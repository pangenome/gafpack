use clap::Parser;
use gafpack::{
    bin_coverages_from_trackers, clip_nodes, compute_coverage, compute_coverage_with_edge_support,
    create_node_bin_trackers, format_coverage_column, parse_gfa_with_edges, cn, ilp,
};
use log::{info, warn, error, debug};

/// Project a GAF alignment file into coverage over GFA graph nodes,
/// optionally estimating copy number using a negative binomial model
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    // ─── Input ───────────────────────────────────────────────────────────────

    /// Input GFA pangenome graph file (supports .gz/.bgz compression)
    #[arg(help_heading = "Input", long)]
    gfa: String,

    /// Input GAF alignment file
    #[arg(help_heading = "Input", short = 'g', long)]
    gaf: String,

    // ─── Output ──────────────────────────────────────────────────────────────

    /// Emit graph coverage vector in a single column
    #[arg(help_heading = "Output", short, long)]
    coverage_column: bool,

    /// Verbosity level: 0=warn, 1=info, 2=debug
    #[arg(help_heading = "Output", short, long, default_value_t = 1)]
    verbose: u8,

    // ─── Coverage Mode ───────────────────────────────────────────────────────

    /// Scale coverage values by node length
    #[arg(help_heading = "Coverage mode", short, long, conflicts_with = "copy_number")]
    len_scale: bool,

    /// Weight coverage by query group occurrences
    #[arg(help_heading = "Coverage mode", short = 'w', long, conflicts_with = "copy_number")]
    weight_queries: bool,

    // ─── Copy Number Estimation ──────────────────────────────────────────────

    /// Enable copy number estimation
    #[arg(help_heading = "Copy number estimation", long)]
    copy_number: bool,

    /// Background CN values to test for parameter estimation
    #[arg(help_heading = "Copy number estimation", long, value_delimiter = ',', default_values_t = vec![1, 2],
          requires = "copy_number")]
    ploidy: Vec<u32>,

    /// Bin size (bp) for negative binomial parameter estimation
    #[arg(help_heading = "Copy number estimation", long, default_value_t = 100,
          requires = "copy_number")]
    bin_size: usize,

    /// CN=0 sensitivity: lower = more aggressive deletion calling
    #[arg(help_heading = "Copy number estimation", long, default_value_t = 0.02,
          requires = "copy_number")]
    epsilon: f64,

    /// Skip ILP flow constraints, use simple ML estimation (faster)
    #[arg(help_heading = "Copy number estimation", long, requires = "copy_number")]
    no_flow: bool,

    /// Disable read deduplication (dedup enabled by default, requires sorted GAF)
    #[arg(help_heading = "Copy number estimation", long, requires = "copy_number")]
    no_dedup: bool,

    // ─── ILP Tuning ──────────────────────────────────────────────────────────

    /// Model complexity: 1=basic, 2=+edge_cov_penalty, 3=+reverse_edge_penalty
    #[arg(help_heading = "ILP tuning", long, default_value_t = 2,
          requires = "copy_number")]
    complexity: u8,

    /// Expensive super-edge penalty (source_prob)
    #[arg(help_heading = "ILP tuning", long, default_value_t = -10000.0,
          requires = "copy_number")]
    source_prob: f64,

    /// Cheap super-edge penalty (cheap_source)
    #[arg(help_heading = "ILP tuning", long, default_value_t = -25.0,
          requires = "copy_number")]
    cheap_source: f64,

    /// Scale factor for log-probabilities (higher = coverage matters more vs flow)
    #[arg(help_heading = "ILP tuning", long, default_value_t = 1.0,
          requires = "copy_number")]
    prob_scale: f64,

    /// CN probability cutoff (auto-computed as 4×|source_prob| if not set)
    #[arg(help_heading = "ILP tuning", long, requires = "copy_number")]
    diff_cutoff: Option<f64>,

    // ─── Solver and Partitioning ─────────────────────────────────────────────

    /// ILP solver: highs (default) or gurobi
    #[cfg(feature = "gurobi")]
    #[arg(help_heading = "Solver", long, default_value = "highs",
          requires = "copy_number")]
    solver: String,

    /// Number of threads for ILP solver (0 = auto)
    #[arg(help_heading = "Solver", short = 't', long, default_value_t = 0,
          requires = "copy_number")]
    threads: u32,

    /// Maximum nodes per ILP partition (0 = no partitioning)
    #[arg(help_heading = "Solver", long, default_value_t = 1_500_000,
          requires = "copy_number")]
    partition_size: usize,

    /// Use METIS graph-aware partitioning
    #[cfg(feature = "metis")]
    #[arg(help_heading = "Solver", long, requires = "copy_number")]
    use_metis: bool,
}

impl Args {
    /// Get solver name (highs is default when gurobi feature not enabled)
    fn solver(&self) -> &str {
        #[cfg(feature = "gurobi")]
        { &self.solver }
        #[cfg(not(feature = "gurobi"))]
        { "highs" }
    }

    /// Check if METIS partitioning is requested
    fn use_metis(&self) -> bool {
        #[cfg(feature = "metis")]
        { self.use_metis }
        #[cfg(not(feature = "metis"))]
        { false }
    }
}

fn main() {
    let args = Args::parse();

    // Initialize logger with verbosity level
    let log_level = match args.verbose {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        _ => log::LevelFilter::Debug,
    };
    env_logger::Builder::new()
        .filter_level(log_level)
        .format_timestamp(None)
        .format_target(false)
        .init();

    if args.copy_number {
        run_copy_number(&args);
    } else {
        run_coverage(&args);
    }
}

/// Original coverage computation mode
fn run_coverage(args: &Args) {
    let (coverage, min_id) =
        compute_coverage(&args.gfa, &args.gaf, args.len_scale, args.weight_queries).unwrap();

    let num_segments = coverage.len();

    if args.coverage_column {
        print!("{}", format_coverage_column(&args.gaf, &coverage));
    } else {
        print!("#sample");
        for n in min_id..min_id + num_segments {
            print!("\tnode.{}", n);
        }
        println!();
        print!("{}", args.gaf);
        for v in coverage {
            print!("\t{}", v);
        }
        println!();
    }
}

/// Copy number estimation mode
fn run_copy_number(args: &Args) {
    info!("Parsing GFA with edges...");

    // Parse GFA with edges (edges will be mutated to track support)
    let (lengths, min_id, mut edges, left_edges, right_edges) =
        parse_gfa_with_edges(&args.gfa).expect("Failed to parse GFA");

    if lengths.is_empty() {
        error!("No segments found in GFA");
        return;
    }

    info!("Found {} segments, {} edges", lengths.len(), edges.len());

    // ─── Node Clipping (Gap 1) ────────────────────────────────────────────────
    // Compute node clipping based on edge overlaps (following floco's graph_processing.py)
    info!("Computing node clipping based on edge overlaps...");
    let clipping = clip_nodes(&lengths, &edges, min_id, &left_edges, &right_edges);

    // Compute clipped lengths for CN calculation
    let clipped_lengths: Vec<usize> = lengths
        .iter()
        .zip(clipping.iter())
        .map(|(&len, clip)| clip.clipped_len(len))
        .collect();

    let total_clipped_bp: usize = clipping.iter()
        .map(|c| c.l_clipping + c.r_clipping)
        .sum();
    let nodes_with_clipping = clipping.iter()
        .filter(|c| c.l_clipping > 0 || c.r_clipping > 0)
        .count();
    debug!("{} nodes have clipping, total {} bp clipped", nodes_with_clipping, total_clipped_bp);

    // Create per-bin trackers using CLIPPED lengths (like floco's bin_nodes())
    debug!("Creating per-bin trackers (bin_size={})...", args.bin_size);
    let mut bin_trackers = create_node_bin_trackers(&clipped_lengths, args.bin_size);
    let nodes_with_bins = bin_trackers.iter().filter(|t| t.is_some()).count();
    debug!("{} nodes have bin trackers", nodes_with_bins);

    // Compute coverage and track edge support + per-bin coverage
    // Deduplication is enabled by default (like floco), use --no-dedup to disable
    let dedup_reads = !args.no_dedup;
    if dedup_reads {
        info!("Computing coverage with read deduplication from GAF...");
    } else {
        info!("Computing coverage from GAF (deduplication disabled)...");
    }
    let (coverage, read_lengths, skipped_alignments) = compute_coverage_with_edge_support(
        &lengths,  // Use original lengths for GAF parsing coordinates
        min_id,
        &args.gaf,
        &mut edges,
        false,  // weight_queries disabled in CN mode (deduplication handles multi-mapping)
        Some(&mut bin_trackers),
        dedup_reads,
    )
    .expect("Failed to compute coverage");

    debug!("Collected {} read lengths for distribution fitting", read_lengths.len());
    if dedup_reads && skipped_alignments > 0 {
        debug!("Skipped {} overlapping alignments during deduplication", skipped_alignments);
    }

    // Log edge support statistics
    let edges_with_support: usize = edges.iter().filter(|e| e.sup_reads > 0).count();
    debug!("{} of {} edges have supporting reads", edges_with_support, edges.len());

    // Extract bins from trackers (actual per-bin coverage, not uniform assumption)
    debug!("Extracting bins from per-bin trackers...");
    let bins = bin_coverages_from_trackers(&bin_trackers);

    if bins.is_empty() {
        warn!("No bins created - nodes may be too short or have insufficient coverage");
        warn!("Using default parameters");
    } else {
        debug!("Created {} bins from nodes with sufficient coverage", bins.len());
    }

    // Estimate NB parameters (logs best ploidy internally)
    info!("Estimating NB parameters (testing ploidies: {:?})...", args.ploidy);
    let (alpha, beta) = cn::estimate_nb_params(&bins, args.bin_size, &args.ploidy);

    // Fit read length distribution (for edge coverage penalty and bin subsampling)
    info!("Fitting read length distribution...");
    let rlen_params = cn::fit_read_length_distribution(&read_lengths);
    let rlen_mean = rlen_params.mean();
    debug!("Read length distribution: mean={:.0}", rlen_mean);

    // ─── Auto-compute diff_cutoff (Gap 5) ─────────────────────────────────────
    // Following floco: diff_cutoff = 4 * |source_prob|
    let source_prob: f64 = args.source_prob;
    let cheap_penalty: f64 = args.cheap_source;
    let diff_cutoff: f64 = args.diff_cutoff.unwrap_or_else(|| {
        let auto_cutoff = 4.0 * source_prob.abs();
        debug!("Auto-computed diff_cutoff = 4 * |{}| = {:.1}", source_prob, auto_cutoff);
        auto_cutoff
    });

    info!("Using complexity level {}, source_prob={:.1}, cheap_source={:.1}, diff_cutoff={:.1}",
          args.complexity, source_prob, cheap_penalty, diff_cutoff);

    // ─── Prepare ILP nodes with per-bin coverage (Gaps 2 & 3) ─────────────────
    debug!("Preparing {} ILP nodes with per-bin coverage...", lengths.len());
    let nodes: Vec<ilp::IlpNode> = (0..lengths.len())
        .map(|i| {
            // Get per-bin coverages for this node (if available)
            let node_bins: Vec<f64> = bin_trackers[i]
                .as_ref()
                .map(|t| t.to_float_bins())
                .unwrap_or_default();

            ilp::IlpNode::new_with_bins(
                min_id + i,
                clipped_lengths[i],  // Use CLIPPED length for CN calculation
                coverage[i],
                &node_bins,
                args.bin_size,
                left_edges[i],
                right_edges[i],
                alpha,
                beta,
                args.epsilon,
                diff_cutoff,
                rlen_mean,  // For bin subsampling
            )
        })
        .collect();

    // Solve

    let cn_calls = if args.no_flow {
        info!("Using simple ML estimation (--no-flow)...");
        ilp::solve_simple(&nodes)
    } else {
        // Decide whether to use partitioning based on graph size
        let use_partitioning = args.partition_size > 0 && nodes.len() > args.partition_size;

        if use_partitioning {
            info!("Solving ILP with flow constraints (partitioned)...");
            match ilp::solve_partitioned(
                &nodes,
                &edges,
                min_id,
                alpha,
                &rlen_params,
                cheap_penalty,
                source_prob,
                args.complexity,
                args.prob_scale,
                args.threads,
                args.partition_size,
                args.use_metis(),
                args.solver(),
            ) {
                Ok(calls) => calls,
                Err(e) => {
                    error!("Partitioned ILP failed: {}", e);
                    warn!("Falling back to simple ML estimation...");
                    ilp::solve_simple(&nodes)
                }
            }
        } else {
            info!("Solving ILP with flow constraints...");
            match ilp::solve(
                &nodes,
                &edges,
                min_id,
                alpha,
                &rlen_params,
                cheap_penalty,
                source_prob,
                args.complexity,
                args.prob_scale,
                args.threads,
                args.solver(),
            ) {
                Ok(calls) => calls,
                Err(e) => {
                    error!("ILP failed: {}", e);
                    warn!("Falling back to simple ML estimation...");
                    ilp::solve_simple(&nodes)
                }
            }
        }
    };

    // Output in same format as coverage mode
    let num_segments = cn_calls.len();

    if args.coverage_column {
        // Column format (PACK format)
        println!("##sample: {}", args.gaf);
        println!("#copy_number");
        for &cn in &cn_calls {
            println!("{}", cn);
        }
    } else {
        // Tab-separated format (same as coverage output)
        print!("#sample");
        for n in min_id..min_id + num_segments {
            print!("\tnode.{}", n);
        }
        println!();
        print!("{}", args.gaf);
        for &cn in &cn_calls {
            print!("\t{}", cn);
        }
        println!();
    }

    info!("Done. Wrote {} CN calls.", cn_calls.len());
}
