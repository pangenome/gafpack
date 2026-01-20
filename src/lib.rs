use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

pub mod cn;
pub mod ilp;
pub mod partition;

// ─── Edge Structure ──────────────────────────────────────────────────────────

/// GFA edge (L line) representing a link between two segments
#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub from_rev: bool,  // true if from node is reverse strand
    pub to_rev: bool,    // true if to node is reverse strand
    pub overlap: usize,
    pub sup_reads: u32,  // Count of reads supporting this edge (traversing it)
}

// ─── Read Length Distribution ────────────────────────────────────────────────

/// Skew-normal distribution parameters for read lengths
/// Used to compute expected edge support based on overlap and read length distribution
#[derive(Clone, Copy, Debug)]
pub struct ReadLengthParams {
    pub shape: f64,   // Skewness parameter (alpha/a)
    pub loc: f64,     // Location parameter (xi)
    pub scale: f64,   // Scale parameter (omega)
}

impl ReadLengthParams {
    /// Create new read length parameters
    pub fn new(shape: f64, loc: f64, scale: f64) -> Self {
        Self { shape, loc, scale }
    }

    /// Standard normal PDF: phi(x) = (1/sqrt(2*pi)) * exp(-x^2/2)
    #[inline]
    fn std_normal_pdf(x: f64) -> f64 {
        const INV_SQRT_2PI: f64 = 0.3989422804014327; // 1/sqrt(2*pi)
        INV_SQRT_2PI * (-0.5 * x * x).exp()
    }

    /// Standard normal CDF approximation using error function
    /// Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
    #[inline]
    fn std_normal_cdf(x: f64) -> f64 {
        0.5 * (1.0 + erf(x * std::f64::consts::FRAC_1_SQRT_2))
    }

    /// Skew-normal PDF: f(x) = (2/omega) * phi((x-xi)/omega) * Phi(alpha * (x-xi)/omega)
    #[allow(dead_code)]
    pub fn pdf(&self, x: f64) -> f64 {
        let z = (x - self.loc) / self.scale;
        (2.0 / self.scale) * Self::std_normal_pdf(z) * Self::std_normal_cdf(self.shape * z)
    }

    /// Skew-normal CDF (computed numerically via survival function)
    #[allow(dead_code)]
    pub fn cdf(&self, x: f64) -> f64 {
        1.0 - self.sf(x)
    }

    /// Survival function: P(X > x) = 1 - CDF(x)
    /// Uses Owen's T function approximation for skew-normal
    pub fn sf(&self, x: f64) -> f64 {
        let z = (x - self.loc) / self.scale;
        // Approximation: SF ≈ 1 - Phi(z) + 2*T(z, alpha)
        // where T is Owen's T function. For simplicity, use numerical integration
        // or approximation based on standard normal.

        // Simple approximation using standard normal and skewness adjustment
        let base_sf = 1.0 - Self::std_normal_cdf(z);

        // Adjust for skewness (when shape > 0, distribution skews right = heavier right tail)
        // This is an approximation; exact would use Owen's T function
        if self.shape.abs() < 0.01 {
            base_sf
        } else {
            // Use delta parameter: delta = alpha / sqrt(1 + alpha^2)
            let delta = self.shape / (1.0 + self.shape * self.shape).sqrt();
            // Approximate adjustment factor
            let adjustment = 2.0 * Self::std_normal_cdf(-delta * z.abs()) * delta.signum();
            (base_sf + adjustment * 0.5 * base_sf).clamp(0.0, 1.0)
        }
    }

    /// Mean of skew-normal distribution: xi + omega * delta * sqrt(2/pi)
    /// where delta = alpha / sqrt(1 + alpha^2)
    pub fn mean(&self) -> f64 {
        let delta = self.shape / (1.0 + self.shape * self.shape).sqrt();
        self.loc + self.scale * delta * (2.0 / std::f64::consts::PI).sqrt()
    }
}

/// Error function approximation (Abramowitz and Stegun)
fn erf(x: f64) -> f64 {
    let a1 =  0.254829592;
    let a2 = -0.284496736;
    let a3 =  1.421413741;
    let a4 = -1.453152027;
    let a5 =  1.061405429;
    let p  =  0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

/// Node information for CN calling
#[derive(Clone, Debug)]
pub struct NodeInfo {
    pub id: usize,
    pub length: usize,
    pub coverage: f64,
    pub left_edges: usize,   // count of edges on left side
    pub right_edges: usize,  // count of edges on right side
}

// ─── Node Clipping ────────────────────────────────────────────────────────────

/// Per-node clipping information (following floco's graph_processing.py:111-158)
/// Tracks how much of each node's left and right ends should be clipped
/// to account for edge overlaps
#[derive(Clone, Debug, Default)]
pub struct NodeClipping {
    pub l_clipping: usize,  // Amount clipped from left end
    pub r_clipping: usize,  // Amount clipped from right end
}

impl NodeClipping {
    /// Calculate clipped length: original_len - l_clipping - r_clipping
    pub fn clipped_len(&self, original_len: usize) -> usize {
        original_len.saturating_sub(self.l_clipping + self.r_clipping)
    }
}

/// Compute node clipping based on edge overlaps
/// Following floco's graph_processing.py:111-158 exactly:
/// 1. Sort edges by DECREASING overlap size
/// 2. For each edge, decide which node to clip based on:
///    a. If one already has enough clipping, clip the other
///    b. If neither is clipped enough, clip the node with MORE edges on that side
///    c. Tie-breaker: clip the LONGER node
///
/// # Arguments
/// * `lengths` - Node lengths indexed by (node_id - min_id)
/// * `edges` - Graph edges
/// * `min_id` - Minimum node ID
/// * `left_edge_counts` - Number of edges on left side of each node
/// * `right_edge_counts` - Number of edges on right side of each node
///
/// # Returns
/// Vector of NodeClipping, indexed by (node_id - min_id)
pub fn clip_nodes(
    lengths: &[usize],
    edges: &[Edge],
    min_id: usize,
    left_edge_counts: &[usize],
    right_edge_counts: &[usize],
) -> Vec<NodeClipping> {
    let n = lengths.len();
    let mut clipping: Vec<NodeClipping> = vec![NodeClipping::default(); n];

    if edges.is_empty() {
        return clipping;
    }

    // Create sorted edge indices by DECREASING overlap
    let mut edge_indices: Vec<usize> = (0..edges.len()).collect();
    edge_indices.sort_by(|&a, &b| edges[b].overlap.cmp(&edges[a].overlap));

    for &edge_idx in &edge_indices {
        let edge = &edges[edge_idx];
        let overlap = edge.overlap;

        if overlap == 0 {
            continue;
        }

        let from_idx = edge.from - min_id;
        let to_idx = edge.to - min_id;

        // Determine which side of each node this edge affects
        // from_rev: true if edge leaves from left side of 'from' node
        // to_rev: true if edge enters right side of 'to' node
        let from_is_left = edge.from_rev;  // Edge leaves from left side
        let to_is_right = edge.to_rev;     // Edge enters right side

        // Current clipping on the affected sides
        let from_current = if from_is_left {
            clipping[from_idx].l_clipping
        } else {
            clipping[from_idx].r_clipping
        };

        let to_current = if to_is_right {
            clipping[to_idx].r_clipping
        } else {
            clipping[to_idx].l_clipping
        };

        // Check if nodes already have enough clipping
        let from_has_enough = from_current >= overlap;
        let to_has_enough = to_current >= overlap;

        // Decide which node to clip
        let clip_from = if from_has_enough && to_has_enough {
            // Both have enough - no change needed
            continue;
        } else if from_has_enough {
            // 'from' has enough, clip 'to'
            false
        } else if to_has_enough {
            // 'to' has enough, clip 'from'
            true
        } else {
            // Neither has enough - use edge count as tie-breaker
            let from_edge_count = if from_is_left {
                left_edge_counts[from_idx]
            } else {
                right_edge_counts[from_idx]
            };

            let to_edge_count = if to_is_right {
                right_edge_counts[to_idx]
            } else {
                left_edge_counts[to_idx]
            };

            if from_edge_count > to_edge_count {
                // 'from' has more edges on this side, clip 'from'
                true
            } else if to_edge_count > from_edge_count {
                // 'to' has more edges on this side, clip 'to'
                false
            } else {
                // Same edge count - clip the LONGER node
                lengths[from_idx] >= lengths[to_idx]
            }
        };

        // Apply the clipping
        if clip_from {
            if from_is_left {
                clipping[from_idx].l_clipping = clipping[from_idx].l_clipping.max(overlap);
            } else {
                clipping[from_idx].r_clipping = clipping[from_idx].r_clipping.max(overlap);
            }
        } else {
            if to_is_right {
                clipping[to_idx].r_clipping = clipping[to_idx].r_clipping.max(overlap);
            } else {
                clipping[to_idx].l_clipping = clipping[to_idx].l_clipping.max(overlap);
            }
        }
    }

    clipping
}

/// Parse GFA file and extract segment information
/// Returns (segment_lengths, min_id) where segment_lengths[id - min_id] gives the length
pub fn parse_gfa(gfa_path: &str) -> std::io::Result<(Vec<usize>, usize)> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    let mut line = String::new();
    let mut segments_map = HashMap::new();
    let mut min_id = usize::MAX;
    let mut max_id = 0;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();

        // Only process segment lines
        if !line_str.starts_with('S') {
            continue;
        }

        // Parse segment line format: S<tab>id<tab>sequence
        let mut fields = line_str.split('\t');
        let Some((id_str, seq)) = fields.next().and_then(|_type| {
            let id_str = fields.next()?;
            let seq = fields.next()?;
            Some((id_str, seq))
        }) else {
            continue;
        };

        // Parse segment ID
        let id = id_str.parse::<usize>().unwrap();
        min_id = min_id.min(id);
        max_id = max_id.max(id);
        segments_map.insert(id, seq.len());
    }

    // Create a dense vector for O(1) access
    let num_segments = max_id - min_id + 1;
    let mut segment_lengths = vec![0; num_segments];

    for (id, len) in segments_map {
        segment_lengths[id - min_id] = len;
    }

    Ok((segment_lengths, min_id))
}

/// Parse GFA file and extract segments AND edges
/// Returns (segment_lengths, min_id, edges, left_edge_counts, right_edge_counts)
pub fn parse_gfa_with_edges(gfa_path: &str) -> std::io::Result<(Vec<usize>, usize, Vec<Edge>, Vec<usize>, Vec<usize>)> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    let mut line = String::new();

    let mut segments_map: HashMap<usize, usize> = HashMap::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut min_id = usize::MAX;
    let mut max_id = 0;

    // First pass: read all lines
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 { break; }

        let line_str = line.trim();

        if line_str.starts_with('S') {
            // Segment line: S<tab>id<tab>sequence
            let mut fields = line_str.split('\t');
            if let Some((id_str, seq)) = fields.next().and_then(|_| {
                let id_str = fields.next()?;
                let seq = fields.next()?;
                Some((id_str, seq))
            }) {
                if let Ok(id) = id_str.parse::<usize>() {
                    min_id = min_id.min(id);
                    max_id = max_id.max(id);
                    segments_map.insert(id, seq.len());
                }
            }
        } else if line_str.starts_with('L') {
            // Link line: L<tab>from<tab>from_orient<tab>to<tab>to_orient<tab>overlap
            let fields: Vec<&str> = line_str.split('\t').collect();
            if fields.len() >= 6 {
                if let (Ok(from), Ok(to)) = (
                    fields[1].parse::<usize>(),
                    fields[3].parse::<usize>()
                ) {
                    let from_rev = fields[2] == "-";
                    let to_rev = fields[4] == "-";
                    // Parse overlap (e.g., "10M" -> 10)
                    let overlap = fields[5]
                        .trim_end_matches('M')
                        .parse::<usize>()
                        .unwrap_or(0);

                    edges.push(Edge { from, to, from_rev, to_rev, overlap, sup_reads: 0 });
                }
            }
        }
    }

    if min_id == usize::MAX {
        return Ok((vec![], 0, vec![], vec![], vec![]));
    }

    // Create dense vectors
    let num_segments = max_id - min_id + 1;
    let mut segment_lengths = vec![0; num_segments];
    let mut left_edges = vec![0; num_segments];
    let mut right_edges = vec![0; num_segments];

    for (id, len) in segments_map {
        segment_lengths[id - min_id] = len;
    }

    // Count edges per side of each node
    for edge in &edges {
        let from_idx = edge.from - min_id;
        let to_idx = edge.to - min_id;

        // Edge leaves from right side of 'from' if forward, left if reverse
        if edge.from_rev {
            left_edges[from_idx] += 1;
        } else {
            right_edges[from_idx] += 1;
        }

        // Edge enters left side of 'to' if forward, right if reverse
        if edge.to_rev {
            right_edges[to_idx] += 1;
        } else {
            left_edges[to_idx] += 1;
        }
    }

    Ok((segment_lengths, min_id, edges, left_edges, right_edges))
}

/// Create a reader that handles compressed files
pub fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::new(decoder);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(file);
        Ok(Box::new(buf_reader))
    }
}

/// Process each step in a GAF alignment line, calculating coverage for graph nodes
///
/// # Arguments
/// * `line` - A GAF format alignment line
/// * `callback` - Function called for each node with (node_id, coverage_length)
/// * `get_node_len` - Function to get the length of a node by its ID
///
/// # Details
/// Parses GAF alignment lines to extract node coverage information:
/// - Handles both forward (>) and reverse (<) node traversals
/// - Adjusts coverage for partial node alignments at path ends
/// - Accumulates coverage across multi-node paths
fn for_each_step(
    line: &str,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
) {
    //eprintln!("{}", line);
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        //eprintln!("oheunotoeunthoue");
        let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
        let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
        let target_len = target_end - target_start;
        //eprintln!("target_len = {}", target_len);
        let fields = line
            .split('\t')
            .nth(5)
            .unwrap()
            .split(['<', '>'])
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().unwrap())
            .enumerate()
            .collect::<Vec<(usize, usize)>>();
        let mut seen: usize = 0;
        let fields_len = fields.as_slice().len();
        //eprintln!("fields len = {}", fields_len);
        for (i, j) in fields {
            let mut len = get_node_len(j);
            //eprintln!("node {} len = {}", j, len);
            if i == 0 {
                //eprintln!("on first step {} {} {}", len, target_start, seen);
                assert!(len >= target_start);
                len -= target_start;
            }
            if i == fields_len - 1 {
                //eprintln!("on last step {} {} {}", len, target_end, seen);
                assert!(target_len >= seen);
                len = target_len - seen;
            }
            if i == fields_len {
                assert!(false);
            }
            //eprintln!("node {} adj len = {}", j, len);
            seen += len;
            callback(j, len);
        }
        //eprintln!("seen = {}", seen);
    }
    //eprintln!("at end");
}

/// Compute coverage from GAF file using pre-parsed segment data
pub fn compute_coverage_with_segments(
    segment_lengths: &[usize],
    min_id: usize,
    gaf_path: &str,
    len_scale: bool,
    weight_queries: bool,
) -> std::io::Result<Vec<f64>> {
    let num_segments = segment_lengths.len();
    let mut coverage: Vec<f64> = vec![0.0; num_segments];

    // Iterates through each line in a file, applying the provided callback function
    let for_each_line = |callback: &mut dyn FnMut(&str)| -> std::io::Result<()> {
        let file = File::open(gaf_path)?;
        let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
        let buf_reader = BufReader::new(reader);
        for line in buf_reader.lines() {
            callback(&line?);
        }
        Ok(())
    };

    if weight_queries {
        // First pass: count query occurrences
        let mut query_counts: HashMap<String, usize> = HashMap::new();
        for_each_line(&mut |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 4 {
                let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
                *query_counts.entry(query_key).or_insert(0) += 1;
            }
        })?;

        // Second pass: calculate coverage with query count adjustment
        for_each_line(&mut |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
            let count = query_counts.get(&query_key).unwrap_or(&1);

            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64 / *count as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        })?;
    } else {
        // Single pass without weighting
        for_each_line(&mut |l: &str| {
            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        })?;
    }

    if len_scale {
        for (i, cov) in coverage.iter_mut().enumerate() {
            *cov /= segment_lengths[i] as f64;
        }
    }

    Ok(coverage)
}

/// Compute coverage from GAF file on a GFA graph
/// Returns coverage vector indexed by (node_id - min_id)
pub fn compute_coverage(
    gfa_path: &str,
    gaf_path: &str,
    len_scale: bool,
    weight_queries: bool,
) -> std::io::Result<(Vec<f64>, usize)> {
    let (segment_lengths, min_id) = parse_gfa(gfa_path)?;
    let coverage = compute_coverage_with_segments(
        &segment_lengths,
        min_id,
        gaf_path,
        len_scale,
        weight_queries,
    )?;
    Ok((coverage, min_id))
}

/// Format coverage as column output (PACK format)
pub fn format_coverage_column(sample_name: &str, coverage: &[f64]) -> String {
    let mut output = String::new();
    output.push_str(&format!("##sample: {}\n", sample_name));
    output.push_str("#coverage\n");
    for &v in coverage {
        output.push_str(&format!("{}\n", v));
    }
    output
}

// ─── Binning for Parameter Estimation ────────────────────────────────────────

/// Per-node bin data for parameter estimation
pub struct NodeBins {
    pub node_idx: usize,
    pub bins: Vec<f64>,
    pub mean_cov: f64,
}

/// Per-node bin tracker for actual coverage tracking during GAF parsing
/// This tracks coverage per bin as alignments are processed (like floco's update_bins)
#[derive(Clone, Debug)]
pub struct NodeBinTracker {
    pub bins: Vec<u64>,      // Per-bin coverage in bp
    pub bin_size: usize,
}

impl NodeBinTracker {
    /// Create a new bin tracker for a node
    /// Only creates trackers for nodes >= bin_size (matching floco's filtering)
    pub fn new(node_length: usize, bin_size: usize) -> Option<Self> {
        if node_length >= bin_size {
            let n_bins = node_length / bin_size;
            Some(Self {
                bins: vec![0; n_bins],
                bin_size,
            })
        } else {
            // Nodes smaller than bin_size don't get trackers (filtered out like in floco)
            None
        }
    }

    /// Update bins with coverage from an alignment region [start, end)
    /// Following floco's update_bins() logic exactly
    pub fn update(&mut self, start: usize, end: usize) {
        if end <= start || self.bins.is_empty() {
            return;
        }

        let bin_size = self.bin_size;
        let i = start / bin_size;  // First covered bin
        let j = ((end - 1) / bin_size + 1).min(self.bins.len());  // Last covered bin + 1

        for k in i..j {
            let b_start = k * bin_size;
            let b_end = (k + 1) * bin_size;
            // Calculate exact overlap with this bin
            let overlap = end.min(b_end).saturating_sub(start.max(b_start));
            self.bins[k] += overlap as u64;
        }
    }

    /// Get mean coverage across all bins
    pub fn mean_coverage(&self) -> f64 {
        if self.bins.is_empty() {
            return 0.0;
        }
        let sum: u64 = self.bins.iter().sum();
        sum as f64 / self.bins.len() as f64
    }

    /// Check if at least one bin has coverage >= bin_size (floco's quality filter)
    pub fn has_sufficient_coverage(&self) -> bool {
        self.bins.iter().any(|&cov| cov >= self.bin_size as u64)
    }

    /// Convert to f64 bins for parameter estimation
    pub fn to_float_bins(&self) -> Vec<f64> {
        self.bins.iter().map(|&b| b as f64).collect()
    }
}

/// Create bin trackers for all nodes before GAF parsing
pub fn create_node_bin_trackers(
    lengths: &[usize],
    bin_size: usize,
) -> Vec<Option<NodeBinTracker>> {
    lengths.iter()
        .map(|&len| NodeBinTracker::new(len, bin_size))
        .collect()
}

/// Divide large nodes into bins for coverage-based parameter estimation
/// Returns per-node bin data with two-level filtering (following floco):
/// 1. Only include nodes with at least one bin >= bin_size coverage
/// 2. Remove top/bottom 3% of nodes by mean bin coverage
///
/// Note: Currently assumes uniform coverage within nodes. For actual per-bin
/// tracking, use compute_coverage_with_bins() during GAF parsing.
pub fn bin_coverages_by_node(
    coverages: &[f64],
    lengths: &[usize],
    bin_size: usize,
) -> Vec<NodeBins> {
    let bin_size_f = bin_size as f64;
    let mut node_bins: Vec<NodeBins> = Vec::new();

    // First pass: collect bins per node
    for (i, &cov) in coverages.iter().enumerate() {
        let len = lengths[i];
        if len >= bin_size && len > 0 {
            let n_bins = len / bin_size;
            // Coverage per bin (assuming uniform distribution)
            let cov_per_bin = cov * bin_size_f / len as f64;

            // Only include if at least one bin has coverage >= bin_size
            if cov_per_bin >= bin_size_f {
                let bins: Vec<f64> = vec![cov_per_bin; n_bins];
                let mean_cov = cov_per_bin;
                node_bins.push(NodeBins {
                    node_idx: i,
                    bins,
                    mean_cov,
                });
            }
        }
    }

    if node_bins.len() < 10 {
        // Not enough nodes for percentile filtering
        return node_bins;
    }

    // Second pass: remove top/bottom 3% of nodes by mean bin coverage
    let mut mean_covs: Vec<f64> = node_bins.iter().map(|nb| nb.mean_cov).collect();
    mean_covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let lo_idx = (mean_covs.len() as f64 * 0.03) as usize;
    let hi_idx = (mean_covs.len() as f64 * 0.97) as usize;
    let lo_thresh = mean_covs[lo_idx];
    let hi_thresh = mean_covs[hi_idx.min(mean_covs.len() - 1)];

    node_bins.retain(|nb| nb.mean_cov >= lo_thresh && nb.mean_cov <= hi_thresh);

    node_bins
}

/// Extract flat bin vector from per-node bins (for parameter estimation)
/// Also applies bin-level percentile filtering (removes top/bottom 3% of bins)
pub fn flatten_bins_with_filtering(node_bins: &[NodeBins]) -> Vec<f64> {
    // Collect all bins
    let mut all_bins: Vec<f64> = node_bins.iter()
        .flat_map(|nb| nb.bins.iter().copied())
        .collect();

    if all_bins.len() < 10 {
        return all_bins;
    }

    // Sort for percentile calculation
    all_bins.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Remove top/bottom 3%
    let lo_idx = (all_bins.len() as f64 * 0.03) as usize;
    let hi_idx = (all_bins.len() as f64 * 0.97) as usize;

    all_bins[lo_idx..hi_idx].to_vec()
}

/// Legacy function for backwards compatibility
/// Divide large nodes into bins for coverage-based parameter estimation
/// Returns flat vector of bin coverages (coverage per bin_size bp)
pub fn bin_coverages(
    coverages: &[f64],
    lengths: &[usize],
    bin_size: usize,
) -> Vec<f64> {
    let node_bins = bin_coverages_by_node(coverages, lengths, bin_size);
    flatten_bins_with_filtering(&node_bins)
}

// ─── Edge Support Tracking ───────────────────────────────────────────────────

/// Step in a path: (node_id, is_reverse)
#[derive(Clone, Copy, Debug)]
struct PathStep {
    node: usize,
    reverse: bool,
}

/// Parse walk field from GAF line, extracting steps with orientations
fn parse_walk(walk: &str) -> Vec<PathStep> {
    let mut steps = Vec::new();
    let mut current_pos = 0;
    let bytes = walk.as_bytes();

    while current_pos < bytes.len() {
        // Find orientation marker
        let is_reverse = bytes[current_pos] == b'<';
        current_pos += 1; // Skip < or >

        // Find end of node ID (next < or > or end)
        let start = current_pos;
        while current_pos < bytes.len() && bytes[current_pos] != b'<' && bytes[current_pos] != b'>' {
            current_pos += 1;
        }

        if let Ok(node) = walk[start..current_pos].parse::<usize>() {
            steps.push(PathStep { node, reverse: is_reverse });
        }
    }

    steps
}

/// Build edge lookup table for fast edge finding
/// Key: (from, to, from_rev, to_rev) -> edge index
pub fn build_edge_index(edges: &[Edge]) -> HashMap<(usize, usize, bool, bool), usize> {
    edges.iter().enumerate()
        .map(|(i, e)| ((e.from, e.to, e.from_rev, e.to_rev), i))
        .collect()
}

// ─── Read Deduplication ──────────────────────────────────────────────────────

/// Number of bins to divide read length into for deduplication tracking
/// Matches floco's READ_BINS = 64
const READ_BINS: usize = 64;

/// State for read deduplication during GAF parsing
/// Tracks which bins of the current read have been covered by alignments
struct ReadDedup {
    current_name: String,
    bin_size: usize,
    covered: u64,  // 64-bit mask tracking covered bins
}

impl ReadDedup {
    fn new() -> Self {
        Self {
            current_name: String::new(),
            bin_size: 0,
            covered: 0,
        }
    }

    /// Check if an alignment should be included based on deduplication
    /// Returns true if the alignment covers new bins (should be included)
    /// Following floco's filter_gaf() logic exactly
    fn should_include(
        &mut self,
        read_name: &str,
        read_len: usize,
        query_start: usize,  // Start position on query (read)
        query_end: usize,    // End position on query (read)
    ) -> bool {
        if read_name != self.current_name {
            // New read - reset state
            self.current_name = read_name.to_string();
            self.bin_size = (read_len + READ_BINS - 1) / READ_BINS;
            if self.bin_size == 0 {
                self.bin_size = 1;
            }

            let start_bin = query_start / self.bin_size;
            let end_bin = ((query_end - 1) / self.bin_size + 1).min(READ_BINS);
            let nbits = end_bin.saturating_sub(start_bin);

            // Create mask for covered bins
            if nbits > 0 && start_bin < READ_BINS {
                self.covered = if nbits >= 64 {
                    u64::MAX
                } else {
                    ((1u64 << nbits) - 1) << start_bin
                };
            } else {
                self.covered = 0;
            }
            true
        } else {
            // Same read - check for overlap with already-covered bins
            if self.bin_size == 0 {
                return false;
            }

            let start_bin = query_start / self.bin_size;
            let end_bin = ((query_end - 1) / self.bin_size + 1).min(READ_BINS);
            let nbits = end_bin.saturating_sub(start_bin);

            if nbits == 0 || start_bin >= READ_BINS {
                return false;
            }

            let mask = if nbits >= 64 {
                u64::MAX
            } else {
                ((1u64 << nbits) - 1) << start_bin
            };

            // Check if there's NO overlap with already-covered bins
            if (self.covered & mask) == 0 {
                // No overlap - include this alignment and mark bins as covered
                self.covered |= mask;
                true
            } else {
                // Overlap exists - skip this alignment
                false
            }
        }
    }
}

/// Compute coverage and track edge support from GAF file
/// Also collects read lengths for distribution fitting
/// Optionally updates per-bin coverage trackers for accurate bin-level coverage
/// Optionally deduplicates multi-mapped reads (requires sorted GAF by read name)
///
/// Returns (coverage, read_lengths, skipped_alignments)
pub fn compute_coverage_with_edge_support(
    segment_lengths: &[usize],
    min_id: usize,
    gaf_path: &str,
    edges: &mut [Edge],
    weight_queries: bool,
    mut bin_trackers: Option<&mut [Option<NodeBinTracker>]>,
    dedup_reads: bool,
) -> std::io::Result<(Vec<f64>, Vec<u32>, usize)> {
    let num_segments = segment_lengths.len();
    let mut coverage: Vec<f64> = vec![0.0; num_segments];
    let mut read_lengths: Vec<u32> = Vec::new();

    // Build edge lookup for fast matching
    let edge_index = build_edge_index(edges);

    // Read deduplication state (like floco's filter_gaf)
    let mut dedup = ReadDedup::new();
    let mut skipped_alignments: usize = 0;

    let file = File::open(gaf_path)?;
    let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
    let buf_reader = BufReader::new(reader);

    // For weighted queries, we need two passes
    let query_counts: HashMap<String, usize> = if weight_queries {
        let mut counts = HashMap::new();
        let file2 = File::open(gaf_path)?;
        let (reader2, _) = niffler::get_reader(Box::new(file2)).unwrap();
        let buf_reader2 = BufReader::new(reader2);
        for line in buf_reader2.lines() {
            let l = line?;
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 4 {
                let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
                *counts.entry(query_key).or_insert(0) += 1;
            }
        }
        counts
    } else {
        HashMap::new()
    };

    for line in buf_reader.lines() {
        let l = line?;
        let fields: Vec<&str> = l.split('\t').collect();

        if fields.len() < 9 {
            continue;
        }

        let walk = fields[5];
        if walk == "*" {
            continue;
        }

        // Get read info for deduplication
        let read_name = fields[0];
        let read_len: usize = fields[1].parse().unwrap_or(0);
        let query_start: usize = fields[2].parse().unwrap_or(0);
        let query_end: usize = fields[3].parse().unwrap_or(0);

        // Apply read deduplication if enabled (like floco's filter_gaf)
        // Note: requires GAF to be sorted by read name
        if dedup_reads && !dedup.should_include(read_name, read_len, query_start, query_end) {
            skipped_alignments += 1;
            continue;
        }

        // Collect read lengths (only for first alignment per read, or all if no dedup)
        if !dedup_reads || read_name != dedup.current_name || dedup.covered == 0 {
            // For dedup mode, only count once per read
        }
        if read_len > 0 {
            // Count read length for all included alignments
            // (floco counts all, but we'll deduplicate later)
            read_lengths.push(read_len as u32);
        }

        // Parse path coordinates
        let target_start: usize = fields[7].parse().unwrap_or(0);
        let target_end: usize = fields[8].parse().unwrap_or(0);
        let target_len = target_end - target_start;

        // Weight for this alignment
        let weight = if weight_queries {
            let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
            let count = *query_counts.get(&query_key).unwrap_or(&1);
            1.0 / count as f64
        } else {
            1.0
        };

        // Parse walk with orientations
        let steps = parse_walk(walk);

        if steps.is_empty() {
            continue;
        }

        // Process coverage (adapted from for_each_step)
        // Also track per-bin coverage with actual positions
        let mut seen: usize = 0;
        for (i, step) in steps.iter().enumerate() {
            let node_idx = step.node - min_id;
            if node_idx >= num_segments {
                continue;
            }

            let node_len = segment_lengths[node_idx];

            // Calculate alignment region within this node
            // start_in_node: where alignment starts within the node (0-based)
            // end_in_node: where alignment ends within the node (exclusive)
            let (start_in_node, len) = if i == 0 {
                // First node: alignment starts at target_start offset
                if node_len >= target_start {
                    (target_start, node_len - target_start)
                } else {
                    (node_len, 0)
                }
            } else {
                // Middle/last nodes: alignment starts at node beginning
                (0, node_len)
            };

            // For last node, cap the length
            let len = if i == steps.len() - 1 {
                if target_len >= seen {
                    target_len - seen
                } else {
                    0
                }
            } else {
                len
            };

            let end_in_node = start_in_node + len;

            seen += len;
            coverage[node_idx] += len as f64 * weight;

            // Update per-bin coverage if trackers are provided
            // Note: we use integer bp for bin tracking (weight is not applied to bins)
            if let Some(ref mut trackers) = bin_trackers.as_deref_mut() {
                if let Some(ref mut tracker) = trackers.get_mut(node_idx).and_then(|t| t.as_mut()) {
                    tracker.update(start_in_node, end_in_node);
                }
            }
        }

        // Track edge traversals for consecutive steps
        for i in 0..steps.len().saturating_sub(1) {
            let from_step = &steps[i];
            let to_step = &steps[i + 1];

            // Determine edge orientation based on path traversal:
            // - from_rev: true if edge leaves from node's left (reverse) side
            // - to_rev: true if edge enters node's right (reverse) side
            //
            // In walk notation:
            // >A means traverse A forward: enter left, exit right
            // <A means traverse A reverse: enter right, exit left
            //
            // So for >A>B (A fwd, B fwd):
            //   - Exit A from right (from_rev=false)
            //   - Enter B from left (to_rev=false)
            // For >A<B (A fwd, B rev):
            //   - Exit A from right (from_rev=false)
            //   - Enter B from right (to_rev=true)
            // For <A>B (A rev, B fwd):
            //   - Exit A from left (from_rev=true)
            //   - Enter B from left (to_rev=false)
            // For <A<B (A rev, B rev):
            //   - Exit A from left (from_rev=true)
            //   - Enter B from right (to_rev=true)

            let from_rev = from_step.reverse;  // Leaving from left side if reversed
            let to_rev = to_step.reverse;      // Entering right side if reversed

            // Look up edge
            let key = (from_step.node, to_step.node, from_rev, to_rev);
            if let Some(&edge_idx) = edge_index.get(&key) {
                edges[edge_idx].sup_reads += 1;
            }
        }
    }

    Ok((coverage, read_lengths, skipped_alignments))
}

/// Extract bins from trackers for parameter estimation
/// Following floco's filtering: only nodes with at least one bin >= bin_size coverage
/// Then removes top/bottom 3% of nodes by mean coverage, then top/bottom 3% of bins
pub fn bin_coverages_from_trackers(
    trackers: &[Option<NodeBinTracker>],
) -> Vec<f64> {
    // First pass: collect per-node bin data with sufficient coverage
    let mut node_bins: Vec<(usize, Vec<f64>, f64)> = Vec::new();  // (idx, bins, mean)

    for (i, tracker_opt) in trackers.iter().enumerate() {
        if let Some(tracker) = tracker_opt {
            // Only include if at least one bin has sufficient coverage (>= bin_size)
            if tracker.has_sufficient_coverage() {
                let bins = tracker.to_float_bins();
                let mean = tracker.mean_coverage();
                if !bins.is_empty() {
                    node_bins.push((i, bins, mean));
                }
            }
        }
    }

    if node_bins.len() < 10 {
        // Not enough nodes for percentile filtering, return all bins
        return node_bins.into_iter().flat_map(|(_, bins, _)| bins).collect();
    }

    // Second pass: remove top/bottom 3% of nodes by mean bin coverage
    let mut mean_covs: Vec<f64> = node_bins.iter().map(|(_, _, mean)| *mean).collect();
    mean_covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let lo_idx = (mean_covs.len() as f64 * 0.03) as usize;
    let hi_idx = (mean_covs.len() as f64 * 0.97) as usize;
    let lo_thresh = mean_covs[lo_idx];
    let hi_thresh = mean_covs[hi_idx.min(mean_covs.len() - 1)];

    node_bins.retain(|(_, _, mean)| *mean >= lo_thresh && *mean <= hi_thresh);

    // Third pass: collect all remaining bins and apply bin-level filtering
    let mut all_bins: Vec<f64> = node_bins.into_iter()
        .flat_map(|(_, bins, _)| bins)
        .collect();

    if all_bins.len() < 10 {
        return all_bins;
    }

    // Sort for percentile calculation
    all_bins.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Remove top/bottom 3%
    let lo_idx = (all_bins.len() as f64 * 0.03) as usize;
    let hi_idx = (all_bins.len() as f64 * 0.97) as usize;

    all_bins[lo_idx..hi_idx].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_node_clipping_empty() {
        let lengths = vec![100, 200, 300];
        let edges: Vec<Edge> = vec![];
        let left_edges = vec![0, 0, 0];
        let right_edges = vec![0, 0, 0];

        let clipping = clip_nodes(&lengths, &edges, 1, &left_edges, &right_edges);

        assert_eq!(clipping.len(), 3);
        assert_eq!(clipping[0].l_clipping, 0);
        assert_eq!(clipping[0].r_clipping, 0);
    }

    #[test]
    fn test_node_clipping_single_edge() {
        // Node 1 (len 100) -> Node 2 (len 200) with overlap 10
        let lengths = vec![100, 200];
        let edges = vec![Edge {
            from: 1,
            to: 2,
            from_rev: false,  // Leaves from right side of node 1
            to_rev: false,    // Enters left side of node 2
            overlap: 10,
            sup_reads: 0,
        }];
        let left_edges = vec![0, 1];  // Node 2 has 1 left edge
        let right_edges = vec![1, 0]; // Node 1 has 1 right edge

        let clipping = clip_nodes(&lengths, &edges, 1, &left_edges, &right_edges);

        assert_eq!(clipping.len(), 2);
        // The longer node (200) should be clipped
        // Since edge counts are equal (1 each), the longer node gets clipped
        assert_eq!(clipping[1].l_clipping, 10);
        assert_eq!(clipping[0].r_clipping, 0);
    }

    #[test]
    fn test_node_clipping_preserves_length() {
        let clipping = NodeClipping { l_clipping: 10, r_clipping: 15 };
        assert_eq!(clipping.clipped_len(100), 75);
    }

    #[test]
    fn test_node_clipping_saturating() {
        // Ensure clipped_len doesn't go negative
        let clipping = NodeClipping { l_clipping: 60, r_clipping: 60 };
        assert_eq!(clipping.clipped_len(100), 0);
    }

    #[test]
    fn test_read_length_params_mean() {
        let params = ReadLengthParams::new(0.0, 1000.0, 200.0);
        // For shape=0 (symmetric), mean should equal loc
        let mean = params.mean();
        assert!((mean - 1000.0).abs() < 1.0);
    }
}
