//! Flow-constrained copy number optimization via Integer Linear Programming
//!
//! Uses HiGHS solver (MIT licensed) through good_lp crate.
//! Implements network flow constraints to ensure biologically plausible CN assignments.

use crate::cn;
use crate::partition::{create_partitions, create_partitions_metis, extract_partition_nodes, filter_edges_for_partition};
use crate::Edge;
use crate::ReadLengthParams;
use good_lp::{constraint, variable, Expression, ProblemVariables, Solution, SolverModel, Variable};
use good_lp::solvers::highs::{highs, HighsParallelType};
use log::{debug, info, warn};
use std::collections::HashMap;

// ─── Constants ───────────────────────────────────────────────────────────────
// Note: PROB_SCALE is now a parameter passed to solve()

// ─── Data Structures ─────────────────────────────────────────────────────────

/// Simple deterministic RNG (LCG) for reproducible bin subsampling
struct SimpleRng(u64);

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> usize {
        // Linear congruential generator
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
        (self.0 >> 33) as usize
    }

    /// Sample n items from a vector without replacement (Fisher-Yates shuffle variant)
    fn sample<T: Clone>(&mut self, items: &[T], n: usize) -> Vec<T> {
        if n >= items.len() {
            return items.to_vec();
        }
        let mut result: Vec<T> = Vec::with_capacity(n);
        let mut indices: Vec<usize> = (0..items.len()).collect();
        for i in 0..n {
            let j = i + (self.next() % (items.len() - i));
            indices.swap(i, j);
            result.push(items[indices[i]].clone());
        }
        result
    }
}

/// Node data prepared for ILP
pub struct IlpNode {
    pub id: usize,
    pub length: usize,
    pub coverage: f64,
    pub cn_lo: i32,          // Lower bound on CN
    pub cn_hi: i32,          // Upper bound on CN
    pub cn_probs: Vec<f64>,  // Log-probabilities for CN=cn_lo, cn_lo+1, ...
    pub left_edges: usize,   // Number of edges on left side
    pub right_edges: usize,  // Number of edges on right side
}

impl IlpNode {
    /// Create ILP node using aggregate coverage (legacy method)
    pub fn new(
        id: usize,
        length: usize,
        coverage: f64,
        left_edges: usize,
        right_edges: usize,
        alpha: f64,
        beta: f64,
        epsilon: f64,
        diff_cutoff: f64,
    ) -> Self {
        let (cn_lo, cn_probs) = cn::cn_log_probs(coverage, length, alpha, beta, epsilon, diff_cutoff);
        let cn_hi = cn_lo + cn_probs.len() as i32 - 1;

        debug!(
            "Node {}: len={}, cov={:.0}, CN range=[{}, {}]",
            id, length, coverage, cn_lo, cn_hi
        );

        IlpNode {
            id,
            length,
            coverage,
            cn_lo,
            cn_hi,
            cn_probs,
            left_edges,
            right_edges,
        }
    }

    /// Create ILP node using per-bin coverages (floco-compatible method)
    /// Following floco's flow_ilp.py:16-29, 84 for bin subsampling
    ///
    /// # Arguments
    /// * `bin_coverages` - Per-bin coverage values
    /// * `bin_size` - Size of each bin in bp
    /// * `rlen_mean` - Mean read length (for subsampling distance calculation)
    pub fn new_with_bins(
        id: usize,
        length: usize,
        coverage: f64,
        bin_coverages: &[f64],
        bin_size: usize,
        left_edges: usize,
        right_edges: usize,
        alpha: f64,
        beta: f64,
        epsilon: f64,
        diff_cutoff: f64,
        rlen_mean: f64,
    ) -> Self {
        // Compute CN probs based on whether we have bins
        let (cn_lo, cn_probs) = if bin_coverages.is_empty() || bin_coverages.len() == 1 {
            // No bins or single bin - use aggregate coverage with node length
            // Following floco: if nbins <= 1, use single value
            cn::cn_log_probs(coverage, length, alpha, beta, epsilon, diff_cutoff)
        } else {
            // Multiple bins - apply subsampling following floco's flow_ilp.py:16-29
            // subsampling_dist = max(1000, rlen_mean)
            // nbins = floor(n_bins * binsize / subsampling_dist)
            let subsampling_dist = 1000.0_f64.max(rlen_mean);
            let total_len = bin_coverages.len() * bin_size;
            let target_nbins = ((total_len as f64) / subsampling_dist).floor() as usize;

            if target_nbins <= 1 {
                // After subsampling we'd have <= 1 bin, use aggregate
                cn::cn_log_probs(coverage, length, alpha, beta, epsilon, diff_cutoff)
            } else {
                // Subsample bins
                let sampled_bins = if target_nbins >= bin_coverages.len() {
                    bin_coverages.to_vec()
                } else {
                    // Use deterministic RNG with node ID as seed for reproducibility
                    let mut rng = SimpleRng::new(id as u64);
                    rng.sample(bin_coverages, target_nbins)
                };

                debug!(
                    "Node {}: subsampled {} bins to {} (subsampling_dist={:.0})",
                    id, bin_coverages.len(), sampled_bins.len(), subsampling_dist
                );

                cn::cn_log_probs_bins(&sampled_bins, bin_size, alpha, beta, epsilon, diff_cutoff)
            }
        };

        let cn_hi = cn_lo + cn_probs.len() as i32 - 1;

        debug!(
            "Node {}: len={}, cov={:.0}, bins={}, CN range=[{}, {}]",
            id, length, coverage, bin_coverages.len(), cn_lo, cn_hi
        );

        IlpNode {
            id,
            length,
            coverage,
            cn_lo,
            cn_hi,
            cn_probs,
            left_edges,
            right_edges,
        }
    }
}

// ─── ILP Solver ──────────────────────────────────────────────────────────────

/// Solve CN calling with flow constraints
/// Returns vector of CN calls indexed by node position
///
/// # Arguments
/// * `nodes` - ILP node data with coverage and CN probabilities
/// * `edges` - Graph edges with support counts
/// * `min_id` - Minimum node ID (for indexing)
/// * `alpha` - Coverage per bp at CN=1 (for edge penalty calculation)
/// * `rlen_params` - Read length distribution (for edge penalty calculation)
/// * `cheap_penalty` - Penalty for edges with insufficient support
/// * `source_prob` - Expensive super-edge penalty
/// * `complexity` - Model complexity: 1=basic, 2=+edge_cov_pen, 3=+reverse_edge_pen
/// * `prob_scale` - Scale factor for log-probabilities (higher = coverage matters more)
/// * `threads` - Number of threads for parallel solving (0 = auto)
pub fn solve(
    nodes: &[IlpNode],
    edges: &[Edge],
    min_id: usize,
    alpha: f64,
    rlen_params: &ReadLengthParams,
    cheap_penalty: f64,
    source_prob: f64,
    complexity: u8,
    prob_scale: f64,
    threads: u32,
) -> Result<Vec<u32>, String> {
    if nodes.is_empty() {
        warn!("No nodes to solve");
        return Ok(vec![]);
    }

    info!("Setting up ILP with {} nodes, {} edges", nodes.len(), edges.len());

    let mut vars = ProblemVariables::new();
    let n = nodes.len();

    // ─── Decision Variables ──────────────────────────────────────────────────

    // CN variable for each node (integer)
    let cn_vars: Vec<Variable> = nodes
        .iter()
        .map(|node| {
            vars.add(variable().integer().min(node.cn_lo).max(node.cn_hi))
        })
        .collect();

    // Indicator variables z[i][k] = 1 if node i has CN = cn_lo + k
    // This linearizes the piecewise log-probability function
    let mut z_vars: Vec<Vec<Variable>> = Vec::with_capacity(n);
    for node in nodes.iter() {
        let num_cn_vals = node.cn_probs.len();
        let z_node: Vec<Variable> = (0..num_cn_vals)
            .map(|_| vars.add(variable().binary()))
            .collect();
        z_vars.push(z_node);
    }

    debug!("Created {} CN variables with indicator variables", cn_vars.len());

    // Flow variable for each edge (integer >= 0)
    let flow_vars: Vec<Variable> = edges
        .iter()
        .map(|_| vars.add(variable().integer().min(0)))
        .collect();

    debug!("Created {} flow variables", flow_vars.len());

    // Super-edge variables: source and sink for left/right of each node
    let src_left: Vec<Variable> = (0..n).map(|_| vars.add(variable().integer().min(0))).collect();
    let src_right: Vec<Variable> = (0..n).map(|_| vars.add(variable().integer().min(0))).collect();
    let snk_left: Vec<Variable> = (0..n).map(|_| vars.add(variable().integer().min(0))).collect();
    let snk_right: Vec<Variable> = (0..n).map(|_| vars.add(variable().integer().min(0))).collect();

    debug!("Created {} super-edge variables", n * 4);

    // Binary indicator variables for super-edge pair penalties (x1/x2)
    // These activate when src/snk is used, penalizing using both on same side
    let x1_left: Vec<Variable> = (0..n).map(|_| vars.add(variable().binary())).collect();
    let x2_left: Vec<Variable> = (0..n).map(|_| vars.add(variable().binary())).collect();
    let x1_right: Vec<Variable> = (0..n).map(|_| vars.add(variable().binary())).collect();
    let x2_right: Vec<Variable> = (0..n).map(|_| vars.add(variable().binary())).collect();

    debug!("Created {} super-edge indicator variables", n * 4);

    // ─── Build Edge Index Maps ───────────────────────────────────────────────

    // Track which edges connect to each side of each node
    let mut left_in: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut left_out: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut right_in: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut right_out: HashMap<usize, Vec<usize>> = HashMap::new();

    for (i, edge) in edges.iter().enumerate() {
        let from_idx = edge.from - min_id;
        let to_idx = edge.to - min_id;

        // Edge leaves from 'from' node
        if edge.from_rev {
            left_out.entry(from_idx).or_default().push(i);
        } else {
            right_out.entry(from_idx).or_default().push(i);
        }

        // Edge enters 'to' node
        if edge.to_rev {
            right_in.entry(to_idx).or_default().push(i);
        } else {
            left_in.entry(to_idx).or_default().push(i);
        }
    }

    // ─── Build Reverse Edge Pair Index ──────────────────────────────────────

    // Find pairs of edges that are reverse complements of each other
    // Edge (from, to, from_rev, to_rev) has reverse (to, from, !to_rev, !from_rev)
    let mut edge_index: HashMap<(usize, usize, bool, bool), usize> = HashMap::new();
    for (i, edge) in edges.iter().enumerate() {
        edge_index.insert((edge.from, edge.to, edge.from_rev, edge.to_rev), i);
    }

    // For each edge, find its reverse pair (if any)
    let mut edge_pairs: Vec<Option<usize>> = vec![None; edges.len()];
    for (i, edge) in edges.iter().enumerate() {
        let rev_key = (edge.to, edge.from, !edge.to_rev, !edge.from_rev);
        if let Some(&rev_idx) = edge_index.get(&rev_key) {
            if i < rev_idx {  // Only store once per pair
                edge_pairs[i] = Some(rev_idx);
            }
        }
    }

    // Create binary indicator variables for edge pairs
    let num_pairs = edge_pairs.iter().filter(|p| p.is_some()).count();
    let mut pair_x1: Vec<Variable> = Vec::with_capacity(num_pairs);
    let mut pair_x2: Vec<Variable> = Vec::with_capacity(num_pairs);
    let mut pair_indices: Vec<(usize, usize)> = Vec::with_capacity(num_pairs);

    for (i, &maybe_rev) in edge_pairs.iter().enumerate() {
        if let Some(rev_idx) = maybe_rev {
            pair_x1.push(vars.add(variable().binary()));
            pair_x2.push(vars.add(variable().binary()));
            pair_indices.push((i, rev_idx));
        }
    }

    debug!("Found {} reverse edge pairs", num_pairs);

    // ─── Objective Function ──────────────────────────────────────────────────

    let mut objective: Expression = Expression::from(0.0);

    // 1. Add log-probability terms using indicator variables
    //    For each node, sum over k: log_prob[k] * z[k]
    //    Scale by prob_scale to balance with flow penalties
    for (i, node) in nodes.iter().enumerate() {
        for (k, &log_prob) in node.cn_probs.iter().enumerate() {
            objective += (log_prob * prob_scale) * z_vars[i][k];
        }
    }

    // 2. Add super-edge penalties (always enabled)
    // Use CLI-provided penalties instead of constants
    for i in 0..n {
        let has_left = nodes[i].left_edges > 0;
        let has_right = nodes[i].right_edges > 0;

        // Penalty depends on whether real edges exist on each side
        let (left_pen, right_pen) = match (has_left, has_right) {
            (true, true) => (source_prob, source_prob),     // double_sides
            (true, false) => (source_prob, cheap_penalty),  // free_right
            (false, true) => (cheap_penalty, source_prob),  // free_left
            (false, false) => (cheap_penalty, cheap_penalty), // free_both
        };

        objective += left_pen * (src_left[i] + snk_left[i]);
        objective += right_pen * (src_right[i] + snk_right[i]);
    }

    // 3. Add edge coverage penalties (complexity >= 2)
    let mut edges_penalized = 0usize;
    if complexity >= 2 {
        for (i, edge) in edges.iter().enumerate() {
            let penalty = cn::edge_cov_penalty(
                edge.sup_reads,
                alpha,
                edge.overlap,
                rlen_params,
                cheap_penalty,
            );
            if penalty < 0.0 {
                edges_penalized += 1;
                // Penalty is applied proportionally to flow through the edge
                objective += penalty * flow_vars[i];
            }
        }
        debug!("{} edges have insufficient support and will be penalized", edges_penalized);
    }

    debug!("Built objective function with log-probs and penalties");

    // 4. Add super-edge pair penalties (complexity >= 3)
    // Penalty paid when BOTH source and sink are used on the same side
    // pen * (x1 + x2 - 1) = 0 when one is used, = pen when both are used
    if complexity >= 3 {
        let x1x2_penalty = 0.8 * source_prob;  // Following floco: 0.8 * source_prob
        for i in 0..n {
            // Left side: penalize using both src_left and snk_left
            objective += x1x2_penalty * (x1_left[i] + x2_left[i] - 1.0);
            // Right side: penalize using both src_right and snk_right
            objective += x1x2_penalty * (x1_right[i] + x2_right[i] - 1.0);
        }

        // 5. Add reverse edge pair penalties
        // Penalty paid when BOTH directions of an edge pair are used
        for p in 0..pair_indices.len() {
            objective += x1x2_penalty * (pair_x1[p] + pair_x2[p] - 1.0);
        }

        debug!("Added {} super-edge pair penalties and {} edge pair penalties", n * 2, pair_indices.len());
    }

    // ─── Constraints ─────────────────────────────────────────────────────────

    // Build problem first, then configure solver
    let unsolved = vars.maximise(objective);

    // Configure HiGHS solver with parallelism
    let mut problem = if threads > 0 {
        debug!("Using {} threads for ILP solver", threads);
        highs(unsolved)
            .set_parallel(HighsParallelType::On)
            .set_threads(threads)
    } else {
        debug!("Using automatic thread selection for ILP solver");
        highs(unsolved)
            .set_parallel(HighsParallelType::On)
    };

    // 1. Indicator variable constraints: exactly one CN value, and CN = sum(k * z[k])
    for (i, node) in nodes.iter().enumerate() {
        // sum_k z[i][k] = 1
        let mut sum_z: Expression = Expression::from(0.0);
        for k in 0..node.cn_probs.len() {
            sum_z += z_vars[i][k];
        }
        problem = problem.with(constraint!(sum_z == 1.0));

        // CN[i] = sum_k (cn_lo + k) * z[i][k]
        let mut cn_expr: Expression = Expression::from(0.0);
        for k in 0..node.cn_probs.len() {
            let cn_val = node.cn_lo + k as i32;
            cn_expr += cn_val as f64 * z_vars[i][k];
        }
        problem = problem.with(constraint!(cn_expr == cn_vars[i]));
    }

    debug!("Added {} indicator constraints", n * 2);

    // 2. Flow conservation constraints (following floco's model)
    //    Node is a "pipe": flow entering one side exits the other
    for i in 0..n {
        // Collect edge flows for each side
        let l_in: Expression = left_in.get(&i)
            .map(|edges| edges.iter().fold(Expression::from(0.0), |acc, &ei| acc + flow_vars[ei]))
            .unwrap_or_else(|| Expression::from(0.0));

        let l_out: Expression = left_out.get(&i)
            .map(|edges| edges.iter().fold(Expression::from(0.0), |acc, &ei| acc + flow_vars[ei]))
            .unwrap_or_else(|| Expression::from(0.0));

        let r_in: Expression = right_in.get(&i)
            .map(|edges| edges.iter().fold(Expression::from(0.0), |acc, &ei| acc + flow_vars[ei]))
            .unwrap_or_else(|| Expression::from(0.0));

        let r_out: Expression = right_out.get(&i)
            .map(|edges| edges.iter().fold(Expression::from(0.0), |acc, &ei| acc + flow_vars[ei]))
            .unwrap_or_else(|| Expression::from(0.0));

        // Flow balance: left_in → right_out, right_in → left_out
        // source_left + l_edges_in == sink_right + r_edges_out
        problem = problem.with(constraint!(
            src_left[i] + l_in.clone() == snk_right[i] + r_out.clone()
        ));

        // source_right + r_edges_in == sink_left + l_edges_out
        problem = problem.with(constraint!(
            src_right[i] + r_in.clone() == snk_left[i] + l_out.clone()
        ));

        // Total inflow = CN
        problem = problem.with(constraint!(
            src_left[i] + src_right[i] + l_in.clone() + r_in.clone() == cn_vars[i]
        ));

        // Total outflow = CN
        problem = problem.with(constraint!(
            snk_left[i] + snk_right[i] + r_out + l_out == cn_vars[i]
        ));
    }

    debug!("Added {} flow conservation constraints", n * 4);

    // 3. x1/x2 indicator constraints (complexity >= 3)
    // These link binary indicators to flows: C * x >= flow activates x when flow > 0
    if complexity >= 3 {
        const BIG_M: f64 = 10000.0;  // Big-M constant

        // Super-edge indicators: x1 activated by src, x2 activated by snk
        for i in 0..n {
            // Left side
            problem = problem.with(constraint!(BIG_M * x1_left[i] >= src_left[i]));
            problem = problem.with(constraint!(BIG_M * x2_left[i] >= snk_left[i]));
            // CRITICAL: At least one of src or snk must be used (floco's x1 + x2 >= 1)
            problem = problem.with(constraint!(x1_left[i] + x2_left[i] >= 1.0));
            // Right side
            problem = problem.with(constraint!(BIG_M * x1_right[i] >= src_right[i]));
            problem = problem.with(constraint!(BIG_M * x2_right[i] >= snk_right[i]));
            // CRITICAL: At least one of src or snk must be used (floco's x1 + x2 >= 1)
            problem = problem.with(constraint!(x1_right[i] + x2_right[i] >= 1.0));
        }

        // Edge pair indicators: x1 activated by forward flow, x2 activated by reverse flow
        for (p, &(fwd_idx, rev_idx)) in pair_indices.iter().enumerate() {
            problem = problem.with(constraint!(BIG_M * pair_x1[p] >= flow_vars[fwd_idx]));
            problem = problem.with(constraint!(BIG_M * pair_x2[p] >= flow_vars[rev_idx]));
            // CRITICAL: At least one direction must be used (floco's x1 + x2 >= 1)
            problem = problem.with(constraint!(pair_x1[p] + pair_x2[p] >= 1.0));
        }

        debug!("Added {} x1/x2 indicator constraints", n * 6 + pair_indices.len() * 3);
    }

    // ─── Solve ───────────────────────────────────────────────────────────────

    info!("Solving ILP...");
    let solution = problem.solve().map_err(|e| format!("ILP solve failed: {:?}", e))?;

    // Extract CN values
    let cn_calls: Vec<u32> = cn_vars
        .iter()
        .map(|&v| solution.value(v).round() as u32)
        .collect();

    // Log summary statistics
    let cn_sum: u32 = cn_calls.iter().sum();
    let cn_max = cn_calls.iter().max().unwrap_or(&0);
    let cn_min = cn_calls.iter().min().unwrap_or(&0);
    info!(
        "ILP solved: {} nodes, CN range [{}, {}], total CN={}",
        cn_calls.len(), cn_min, cn_max, cn_sum
    );

    // Count CN distribution
    let mut cn_dist: HashMap<u32, usize> = HashMap::new();
    for &cn in &cn_calls {
        *cn_dist.entry(cn).or_insert(0) += 1;
    }
    let mut cn_dist_vec: Vec<_> = cn_dist.into_iter().collect();
    cn_dist_vec.sort_by_key(|&(cn, _)| cn);
    debug!("CN distribution: {:?}", cn_dist_vec);

    // ─── Diagnostics: Super-edge usage and flow-constrained nodes ────────────

    // Count super-edge usage
    let mut src_left_total = 0.0;
    let mut src_right_total = 0.0;
    let mut snk_left_total = 0.0;
    let mut snk_right_total = 0.0;
    let mut nodes_with_super_edges = 0usize;

    for i in 0..n {
        let sl = solution.value(src_left[i]);
        let sr = solution.value(src_right[i]);
        let kl = solution.value(snk_left[i]);
        let kr = solution.value(snk_right[i]);

        src_left_total += sl;
        src_right_total += sr;
        snk_left_total += kl;
        snk_right_total += kr;

        if sl > 0.5 || sr > 0.5 || kl > 0.5 || kr > 0.5 {
            nodes_with_super_edges += 1;
        }
    }

    debug!(
        "Super-edge usage: src_left={:.0}, src_right={:.0}, snk_left={:.0}, snk_right={:.0}",
        src_left_total, src_right_total, snk_left_total, snk_right_total
    );
    debug!("{} nodes use super-edges (flow enters/exits at graph boundaries)", nodes_with_super_edges);

    // Find nodes where ILP CN differs from ML estimate (flow-constrained)
    let mut flow_constrained_nodes: Vec<(usize, u32, u32)> = Vec::new();
    for (i, node) in nodes.iter().enumerate() {
        // ML estimate: argmax of log-probs
        let ml_cn = node.cn_probs
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(k, _)| (node.cn_lo + k as i32).max(0) as u32)
            .unwrap_or(0);

        let ilp_cn = cn_calls[i];

        if ml_cn != ilp_cn {
            flow_constrained_nodes.push((min_id + i, ml_cn, ilp_cn));
        }
    }

    if !flow_constrained_nodes.is_empty() {
        debug!(
            "{} nodes have CN changed by flow constraints:",
            flow_constrained_nodes.len()
        );
        for (node_id, ml_cn, ilp_cn) in &flow_constrained_nodes {
            debug!("  Node {}: ML={} -> ILP={}", node_id, ml_cn, ilp_cn);
        }
    } else {
        debug!("All nodes have CN matching ML estimate (no flow corrections)");
    }

    Ok(cn_calls)
}

/// Solve CN calling with flow constraints using graph partitioning
///
/// Splits large graphs into smaller partitions to avoid solver overflow.
/// Each partition is solved independently, with boundary nodes using the
/// cheap_penalty mechanism for flow entering/exiting the partition.
///
/// # Arguments
/// * `nodes` - ILP node data with coverage and CN probabilities
/// * `edges` - Graph edges with support counts
/// * `min_id` - Minimum node ID (for indexing)
/// * `alpha` - Coverage per bp at CN=1 (for edge penalty calculation)
/// * `rlen_params` - Read length distribution (for edge penalty calculation)
/// * `cheap_penalty` - Penalty for edges with insufficient support
/// * `source_prob` - Expensive super-edge penalty
/// * `complexity` - Model complexity: 1=basic, 2=+edge_cov_pen, 3=+reverse_edge_pen
/// * `prob_scale` - Scale factor for log-probabilities
/// * `threads` - Number of threads for parallel solving (0 = auto)
/// * `partition_size` - Maximum nodes per partition
/// * `use_metis` - Use METIS graph-aware partitioning instead of sequential
pub fn solve_partitioned(
    nodes: &[IlpNode],
    edges: &[Edge],
    min_id: usize,
    alpha: f64,
    rlen_params: &ReadLengthParams,
    cheap_penalty: f64,
    source_prob: f64,
    complexity: u8,
    prob_scale: f64,
    threads: u32,
    partition_size: usize,
    use_metis: bool,
) -> Result<Vec<u32>, String> {
    let num_nodes = nodes.len();

    // Create partitions
    let partitions = if use_metis {
        let num_partitions = (num_nodes + partition_size - 1) / partition_size;
        info!("Using METIS graph-aware partitioning into {} partitions", num_partitions);
        create_partitions_metis(num_nodes, edges, min_id, num_partitions)?
    } else {
        create_partitions(min_id, num_nodes, partition_size)
    };

    if partitions.len() == 1 {
        // No partitioning needed, use regular solve
        info!("Graph fits in single partition, using standard solve");
        return solve(
            nodes, edges, min_id, alpha, rlen_params,
            cheap_penalty, source_prob, complexity, prob_scale, threads,
        );
    }

    info!(
        "Partitioning {} nodes into {} partitions (max {} nodes each)",
        num_nodes, partitions.len(), partition_size
    );

    // Scale cheap_penalty based on partition ratio - smaller partitions need stricter penalty
    // because they have proportionally more boundary nodes (higher surface-to-volume ratio)
    let partition_ratio = num_nodes as f64 / partition_size as f64;
    let scaling_factor = partition_ratio.log2().max(1.0);
    let partition_cheap_penalty = cheap_penalty * scaling_factor;
    info!(
        "Scaling cheap_penalty by {:.2}x for partitioning: {:.1} → {:.1}",
        scaling_factor, cheap_penalty, partition_cheap_penalty
    );

    // Allocate output vector
    let mut cn_calls = vec![0u32; num_nodes];

    // Solve each partition
    for (i, partition) in partitions.iter().enumerate() {
        info!(
            "Solving partition {}/{}: {} nodes (local_min_id={})",
            i + 1, partitions.len(),
            partition.num_nodes(),
            partition.local_min_id()
        );

        // Filter edges to only those internal to this partition
        let (local_edges, local_left_edges, local_right_edges) =
            filter_edges_for_partition(edges, partition, min_id);

        debug!(
            "Partition {}: {} internal edges (dropped {} cross-partition)",
            i + 1, local_edges.len(), edges.len() - local_edges.len()
        );

        // Extract nodes with updated edge counts
        let local_nodes = extract_partition_nodes(
            nodes, partition, &local_left_edges, &local_right_edges,
        );

        // Solve this partition
        let partition_result = solve(
            &local_nodes,
            &local_edges,
            partition.local_min_id(),
            alpha,
            rlen_params,
            partition_cheap_penalty,
            source_prob,
            complexity,
            prob_scale,
            threads,
        )?;

        // Copy results to global output
        // For sparse partitions, we need to use the actual global indices
        let global_indices = partition.global_indices();
        for (local_idx, &cn) in partition_result.iter().enumerate() {
            let global_idx = global_indices[local_idx];
            cn_calls[global_idx] = cn;
        }
    }

    // Log summary statistics
    let cn_sum: u32 = cn_calls.iter().sum();
    let cn_max = cn_calls.iter().max().unwrap_or(&0);
    let cn_min = cn_calls.iter().min().unwrap_or(&0);
    info!(
        "Partitioned ILP complete: {} nodes, CN range [{}, {}], total CN={}",
        cn_calls.len(), cn_min, cn_max, cn_sum
    );

    Ok(cn_calls)
}

/// Simpler solve without flow constraints - just use ML estimate
pub fn solve_simple(nodes: &[IlpNode]) -> Vec<u32> {
    info!("Using simple ML estimation (no flow constraints)");

    let cn_calls: Vec<u32> = nodes
        .iter()
        .map(|node| {
            node.cn_probs
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, _)| (node.cn_lo + i as i32).max(0) as u32)
                .unwrap_or(0)
        })
        .collect();

    // Log summary statistics
    let cn_sum: u32 = cn_calls.iter().sum();
    let cn_max = cn_calls.iter().max().unwrap_or(&0);
    let cn_min = cn_calls.iter().min().unwrap_or(&0);
    info!(
        "ML estimation done: {} nodes, CN range [{}, {}], total CN={}",
        cn_calls.len(), cn_min, cn_max, cn_sum
    );

    cn_calls
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_simple() {
        let nodes = vec![
            IlpNode {
                id: 1,
                length: 100,
                coverage: 200.0,
                cn_lo: 1,
                cn_hi: 4,
                cn_probs: vec![-10.0, -0.1, -5.0, -15.0], // CN=2 is most likely
                left_edges: 0,
                right_edges: 0,
            },
        ];

        let result = solve_simple(&nodes);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 2); // Should pick CN=2 (index 1 in probs, lo=1)
    }
}
