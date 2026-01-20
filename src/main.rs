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
    /// Input GFA pangenome graph file (supports .gz/.bgz compression)
    #[arg(long)]
    gfa: String,

    /// Input GAF alignment file
    #[arg(short, long)]
    gaf: String,

    /// Scale coverage values by node length
    #[arg(short, long)]
    len_scale: bool,

    /// Emit graph coverage vector in a single column
    #[arg(short, long)]
    coverage_column: bool,

    /// Weight coverage by query group occurrences
    #[arg(short = 'w', long)]
    weight_queries: bool,

    // ─── Copy Number Options ─────────────────────────────────────────────────

    /// Enable copy number estimation
    #[arg(long)]
    copy_number: bool,

    /// Background CN values to test for parameter estimation [default: 1,2]
    #[arg(long, value_delimiter = ',', default_values_t = vec![1, 2])]
    ploidy: Vec<u32>,

    /// Bin size (bp) for negative binomial parameter estimation
    #[arg(long, default_value_t = 100)]
    bin_size: usize,

    /// CN=0 sensitivity: lower = more aggressive deletion calling
    #[arg(long, default_value_t = 0.02)]
    epsilon: f64,

    /// Skip ILP flow constraints, use simple ML estimation (faster)
    #[arg(long)]
    no_flow: bool,

    /// Model complexity: 1=basic, 2=+edge_cov_penalty, 3=+reverse_edge_penalty
    #[arg(long, default_value_t = 2)]
    complexity: u8,

    /// Expensive super-edge penalty (source_prob)
    #[arg(long, default_value_t = -10000.0)]
    source_prob: f64,

    /// Cheap super-edge penalty (cheap_source)
    #[arg(long, default_value_t = -25.0)]
    cheap_source: f64,

    /// Scale factor for log-probabilities in ILP (higher = coverage matters more vs flow)
    #[arg(long, default_value_t = 1.0)]
    prob_scale: f64,

    /// CN probability cutoff: controls how many CN values are considered per node
    /// If not specified, auto-computed as 4 * |source_prob| (floco-compatible)
    #[arg(long)]
    diff_cutoff: Option<f64>,

    /// Number of threads for ILP solver (0 = auto)
    #[arg(short = 't', long, default_value_t = 0)]
    threads: u32,

    /// Maximum nodes per ILP partition (0 = no partitioning)
    #[arg(long, default_value_t = 1_500_000)]
    partition_size: usize,

    /// Disable read deduplication (dedup is enabled by default for CN mode)
    /// Deduplication filters overlapping multi-mappings per read (requires sorted GAF)
    #[arg(long)]
    no_dedup: bool,

    /// Verbosity level: 0=warn, 1=info, 2=debug [default: 1]
    #[arg(short, long, default_value_t = 1)]
    verbose: u8,
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
    // Warn about ignored options
    if args.len_scale {
        warn!("--len-scale is ignored in copy-number mode (length normalization is handled internally)");
    }
    if args.weight_queries {
        warn!("--weight-queries is ignored in copy-number mode (read deduplication handles multi-mapping)");
    }

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
