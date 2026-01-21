//! Graph partitioning for large-scale ILP solving
//!
//! Splits large graphs into manageable chunks to avoid solver overflow issues.
//! Supports both sequential (contiguous) partitioning and METIS graph-aware partitioning.

use crate::ilp::IlpNode;
use crate::Edge;
use std::collections::HashMap;

/// A partition of the graph - either contiguous range or explicit node list
#[derive(Clone, Debug)]
pub enum Partition {
    /// Contiguous node ID range (sequential partitioning)
    Contiguous {
        /// Start index in global node array (inclusive)
        global_start: usize,
        /// End index in global node array (exclusive)
        global_end: usize,
        /// Number of nodes in this partition
        num_nodes: usize,
        /// The min_id for this partition (global min_id + global_start)
        local_min_id: usize,
    },
    /// Non-contiguous node set (METIS partitioning)
    Sparse {
        /// Global indices of nodes in this partition
        node_indices: Vec<usize>,
        /// Mapping from global index to local index
        global_to_local: HashMap<usize, usize>,
        /// The local min_id for this partition (same as global min_id)
        local_min_id: usize,
    },
}

impl Partition {
    /// Get number of nodes in this partition
    pub fn num_nodes(&self) -> usize {
        match self {
            Partition::Contiguous { num_nodes, .. } => *num_nodes,
            Partition::Sparse { node_indices, .. } => node_indices.len(),
        }
    }

    /// Get the local_min_id for this partition
    pub fn local_min_id(&self) -> usize {
        match self {
            Partition::Contiguous { local_min_id, .. } => *local_min_id,
            Partition::Sparse { local_min_id, .. } => *local_min_id,
        }
    }

    /// Check if a global node index belongs to this partition
    #[inline]
    pub fn contains(&self, global_idx: usize) -> bool {
        match self {
            Partition::Contiguous { global_start, global_end, .. } => {
                global_idx >= *global_start && global_idx < *global_end
            }
            Partition::Sparse { global_to_local, .. } => {
                global_to_local.contains_key(&global_idx)
            }
        }
    }

    /// Convert global node index to local partition index
    #[inline]
    pub fn to_local(&self, global_idx: usize) -> usize {
        match self {
            Partition::Contiguous { global_start, .. } => {
                debug_assert!(self.contains(global_idx));
                global_idx - global_start
            }
            Partition::Sparse { global_to_local, .. } => {
                *global_to_local.get(&global_idx).expect("Node not in partition")
            }
        }
    }

    /// Get the global start index (for contiguous) or minimum index (for sparse)
    pub fn global_start(&self) -> usize {
        match self {
            Partition::Contiguous { global_start, .. } => *global_start,
            Partition::Sparse { node_indices, .. } => {
                *node_indices.iter().min().unwrap_or(&0)
            }
        }
    }

    /// Iterate over global indices in this partition
    pub fn global_indices(&self) -> Vec<usize> {
        match self {
            Partition::Contiguous { global_start, global_end, .. } => {
                (*global_start..*global_end).collect()
            }
            Partition::Sparse { node_indices, .. } => node_indices.clone(),
        }
    }
}

/// Create sequential partitions for a graph with the given number of nodes
///
/// # Arguments
/// * `min_id` - Global minimum node ID
/// * `num_nodes` - Total number of nodes
/// * `partition_size` - Maximum nodes per partition (0 = single partition)
///
/// # Returns
/// Vector of Partition structs covering all nodes
pub fn create_partitions(min_id: usize, num_nodes: usize, partition_size: usize) -> Vec<Partition> {
    if partition_size == 0 || num_nodes <= partition_size {
        // No partitioning needed
        return vec![Partition::Contiguous {
            global_start: 0,
            global_end: num_nodes,
            num_nodes,
            local_min_id: min_id,
        }];
    }

    let mut partitions = Vec::new();
    let mut start = 0;

    while start < num_nodes {
        let end = (start + partition_size).min(num_nodes);
        partitions.push(Partition::Contiguous {
            global_start: start,
            global_end: end,
            num_nodes: end - start,
            local_min_id: min_id + start,
        });
        start = end;
    }

    partitions
}

/// Build CSR (Compressed Sparse Row) adjacency representation for METIS
///
/// METIS requires undirected graph in CSR format:
/// - xadj[i]..xadj[i+1] gives indices into adjncy for node i's neighbors
/// - adjncy contains the actual neighbor node indices
///
/// # Arguments
/// * `num_nodes` - Number of nodes in the graph
/// * `edges` - Graph edges
/// * `min_id` - Minimum node ID (for converting to 0-based indices)
///
/// # Returns
/// (xadj, adjncy) tuple for METIS
#[cfg(feature = "metis")]
fn build_csr(
    num_nodes: usize,
    edges: &[Edge],
    min_id: usize,
) -> (Vec<i32>, Vec<i32>) {
    // Count edges per node (both directions for undirected graph)
    let mut degrees = vec![0usize; num_nodes];
    for e in edges {
        let from_idx = e.from - min_id;
        let to_idx = e.to - min_id;
        // Skip self-loops
        if from_idx != to_idx {
            degrees[from_idx] += 1;
            degrees[to_idx] += 1;
        }
    }

    // Build xadj (prefix sum of degrees)
    let mut xadj = vec![0i32; num_nodes + 1];
    for i in 0..num_nodes {
        xadj[i + 1] = xadj[i] + degrees[i] as i32;
    }

    // Build adjncy
    let mut adjncy = vec![0i32; xadj[num_nodes] as usize];
    let mut pos = xadj[..num_nodes].to_vec();
    for e in edges {
        let u = e.from - min_id;
        let v = e.to - min_id;
        // Skip self-loops
        if u != v {
            adjncy[pos[u] as usize] = v as i32;
            pos[u] += 1;
            adjncy[pos[v] as usize] = u as i32;
            pos[v] += 1;
        }
    }

    (xadj, adjncy)
}

/// Create partitions using METIS graph partitioning
///
/// METIS minimizes edge cuts, which reduces the number of boundary nodes
/// that rely on cheap_penalty super-edges.
///
/// # Arguments
/// * `num_nodes` - Total number of nodes
/// * `edges` - Graph edges
/// * `min_id` - Global minimum node ID
/// * `num_partitions` - Number of partitions to create
///
/// # Returns
/// Vector of Partition structs (Sparse variant)
#[cfg(feature = "metis")]
pub fn create_partitions_metis(
    num_nodes: usize,
    edges: &[Edge],
    min_id: usize,
    num_partitions: usize,
) -> Result<Vec<Partition>, String> {
    use log::debug;

    if num_partitions <= 1 {
        // No partitioning needed
        return Ok(vec![Partition::Contiguous {
            global_start: 0,
            global_end: num_nodes,
            num_nodes,
            local_min_id: min_id,
        }]);
    }

    debug!("Building CSR representation for METIS ({} nodes, {} edges)", num_nodes, edges.len());

    // Build CSR adjacency (undirected - add both directions)
    let (xadj, adjncy) = build_csr(num_nodes, edges, min_id);

    debug!("CSR: xadj len={}, adjncy len={}", xadj.len(), adjncy.len());

    // Call METIS to partition
    let mut part = vec![0i32; num_nodes];

    // METIS requires at least some edges to partition
    if adjncy.is_empty() {
        debug!("Graph has no edges, falling back to sequential partitioning");
        return Ok(create_partitions(min_id, num_nodes, (num_nodes + num_partitions - 1) / num_partitions));
    }

    let result = metis::Graph::new(
        1,                    // numbering starts at 0
        num_partitions as i32,
        &xadj,
        &adjncy,
    );

    let graph = match result {
        Ok(g) => g,
        Err(e) => return Err(format!("METIS graph creation failed: {:?}", e)),
    };

    match graph.part_kway(&mut part) {
        Ok(_) => {}
        Err(e) => return Err(format!("METIS partitioning failed: {:?}", e)),
    }

    // Convert partition vector to Partition structs
    // Group nodes by partition ID
    let mut partition_nodes: Vec<Vec<usize>> = vec![Vec::new(); num_partitions];
    for (node_idx, &part_id) in part.iter().enumerate() {
        partition_nodes[part_id as usize].push(node_idx);
    }

    // Create Partition structs
    let partitions: Vec<Partition> = partition_nodes
        .into_iter()
        .filter(|nodes| !nodes.is_empty())
        .map(|nodes| {
            let global_to_local: HashMap<usize, usize> = nodes
                .iter()
                .enumerate()
                .map(|(local, &global)| (global, local))
                .collect();
            Partition::Sparse {
                node_indices: nodes,
                global_to_local,
                local_min_id: min_id,
            }
        })
        .collect();

    // Count edge cuts
    let mut edge_cuts = 0usize;
    for edge in edges {
        let from_idx = edge.from - min_id;
        let to_idx = edge.to - min_id;
        if part[from_idx] != part[to_idx] {
            edge_cuts += 1;
        }
    }
    debug!("METIS partitioning: {} partitions, {} edge cuts", partitions.len(), edge_cuts);

    Ok(partitions)
}

/// Stub for when METIS feature is not enabled
#[cfg(not(feature = "metis"))]
pub fn create_partitions_metis(
    _num_nodes: usize,
    _edges: &[Edge],
    _min_id: usize,
    _num_partitions: usize,
) -> Result<Vec<Partition>, String> {
    Err("METIS feature not enabled. Compile with --features metis".into())
}

/// Filter edges to only include those internal to a partition
///
/// Cross-partition edges are dropped. Nodes at partition boundaries will have
/// reduced left_edges/right_edges counts, triggering the cheap_penalty mechanism.
///
/// # Arguments
/// * `edges` - All graph edges
/// * `partition` - The partition to filter for
/// * `global_min_id` - Global minimum node ID
///
/// # Returns
/// (filtered_edges, local_left_edges, local_right_edges)
pub fn filter_edges_for_partition(
    edges: &[Edge],
    partition: &Partition,
    global_min_id: usize,
) -> (Vec<Edge>, Vec<usize>, Vec<usize>) {
    let num_nodes = partition.num_nodes();
    let mut local_edges = Vec::new();
    let mut local_left_edges = vec![0usize; num_nodes];
    let mut local_right_edges = vec![0usize; num_nodes];

    for edge in edges {
        let from_global_idx = edge.from - global_min_id;
        let to_global_idx = edge.to - global_min_id;

        // Only include edges where both endpoints are in this partition
        if partition.contains(from_global_idx) && partition.contains(to_global_idx) {
            // For sparse partitions, we need to remap node IDs
            let remapped_edge = match partition {
                Partition::Contiguous { .. } => {
                    // Keep original edge IDs - the solver uses local_min_id to compute indices
                    *edge
                }
                Partition::Sparse { global_to_local, local_min_id, .. } => {
                    // Remap to local indices for sparse partitions
                    let from_local = global_to_local[&from_global_idx];
                    let to_local = global_to_local[&to_global_idx];
                    Edge {
                        from: local_min_id + from_local,
                        to: local_min_id + to_local,
                        from_rev: edge.from_rev,
                        to_rev: edge.to_rev,
                        overlap: edge.overlap,
                        sup_reads: edge.sup_reads,
                    }
                }
            };
            local_edges.push(remapped_edge);

            // Update local edge counts
            let from_local = partition.to_local(from_global_idx);
            let to_local = partition.to_local(to_global_idx);

            if edge.from_rev {
                local_left_edges[from_local] += 1;
            } else {
                local_right_edges[from_local] += 1;
            }

            if edge.to_rev {
                local_right_edges[to_local] += 1;
            } else {
                local_left_edges[to_local] += 1;
            }
        }
    }

    // For boundary nodes, we use the counts from internal edges only.
    // This means boundary nodes will have fewer edges than in the global graph,
    // triggering the cheap_penalty super-edge mechanism for flow in/out.

    (local_edges, local_left_edges, local_right_edges)
}

/// Extract IlpNode data for a partition with updated edge counts
///
/// # Arguments
/// * `global_nodes` - All ILP nodes for the full graph
/// * `partition` - The partition to extract
/// * `local_left_edges` - Local left edge counts (from filter_edges_for_partition)
/// * `local_right_edges` - Local right edge counts (from filter_edges_for_partition)
///
/// # Returns
/// Vector of IlpNode for this partition, with updated edge counts
pub fn extract_partition_nodes(
    global_nodes: &[IlpNode],
    partition: &Partition,
    local_left_edges: &[usize],
    local_right_edges: &[usize],
) -> Vec<IlpNode> {
    match partition {
        Partition::Contiguous { global_start, global_end, .. } => {
            (*global_start..*global_end)
                .enumerate()
                .map(|(local_idx, global_idx)| {
                    let node = &global_nodes[global_idx];
                    IlpNode {
                        id: node.id,  // Keep original global ID for debugging
                        length: node.length,
                        coverage: node.coverage,
                        cn_lo: node.cn_lo,
                        cn_hi: node.cn_hi,
                        cn_probs: node.cn_probs.clone(),
                        left_edges: local_left_edges[local_idx],
                        right_edges: local_right_edges[local_idx],
                    }
                })
                .collect()
        }
        Partition::Sparse { node_indices, .. } => {
            node_indices
                .iter()
                .enumerate()
                .map(|(local_idx, &global_idx)| {
                    let node = &global_nodes[global_idx];
                    IlpNode {
                        id: node.id,  // Keep original global ID for debugging
                        length: node.length,
                        coverage: node.coverage,
                        cn_lo: node.cn_lo,
                        cn_hi: node.cn_hi,
                        cn_probs: node.cn_probs.clone(),
                        left_edges: local_left_edges[local_idx],
                        right_edges: local_right_edges[local_idx],
                    }
                })
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_partitions_no_split() {
        let partitions = create_partitions(100, 500, 1000);
        assert_eq!(partitions.len(), 1);
        match &partitions[0] {
            Partition::Contiguous { global_start, global_end, local_min_id, .. } => {
                assert_eq!(*global_start, 0);
                assert_eq!(*global_end, 500);
                assert_eq!(*local_min_id, 100);
            }
            _ => panic!("Expected Contiguous partition"),
        }
    }

    #[test]
    fn test_create_partitions_split() {
        let partitions = create_partitions(100, 2500, 1000);
        assert_eq!(partitions.len(), 3);

        match &partitions[0] {
            Partition::Contiguous { global_start, global_end, local_min_id, .. } => {
                assert_eq!(*global_start, 0);
                assert_eq!(*global_end, 1000);
                assert_eq!(*local_min_id, 100);  // min_id + 0
            }
            _ => panic!("Expected Contiguous partition"),
        }

        match &partitions[1] {
            Partition::Contiguous { global_start, global_end, local_min_id, .. } => {
                assert_eq!(*global_start, 1000);
                assert_eq!(*global_end, 2000);
                assert_eq!(*local_min_id, 1100);  // min_id + 1000
            }
            _ => panic!("Expected Contiguous partition"),
        }

        match &partitions[2] {
            Partition::Contiguous { global_start, global_end, local_min_id, .. } => {
                assert_eq!(*global_start, 2000);
                assert_eq!(*global_end, 2500);
                assert_eq!(*local_min_id, 2100);  // min_id + 2000
            }
            _ => panic!("Expected Contiguous partition"),
        }
    }

    #[test]
    fn test_partition_contains_contiguous() {
        let partition = Partition::Contiguous {
            global_start: 100,
            global_end: 200,
            num_nodes: 100,
            local_min_id: 100,
        };

        assert!(!partition.contains(99));
        assert!(partition.contains(100));
        assert!(partition.contains(150));
        assert!(partition.contains(199));
        assert!(!partition.contains(200));
    }

    #[test]
    fn test_partition_contains_sparse() {
        let mut global_to_local = HashMap::new();
        global_to_local.insert(5, 0);
        global_to_local.insert(10, 1);
        global_to_local.insert(15, 2);

        let partition = Partition::Sparse {
            node_indices: vec![5, 10, 15],
            global_to_local,
            local_min_id: 0,
        };

        assert!(!partition.contains(0));
        assert!(partition.contains(5));
        assert!(!partition.contains(7));
        assert!(partition.contains(10));
        assert!(partition.contains(15));
        assert!(!partition.contains(20));
    }

    #[test]
    fn test_filter_edges_internal() {
        let edges = vec![
            Edge { from: 10, to: 11, from_rev: false, to_rev: false, overlap: 0, sup_reads: 1 },
            Edge { from: 11, to: 12, from_rev: false, to_rev: false, overlap: 0, sup_reads: 2 },
            Edge { from: 12, to: 15, from_rev: false, to_rev: false, overlap: 0, sup_reads: 3 }, // crosses boundary
        ];

        let partition = Partition::Contiguous {
            global_start: 0,
            global_end: 3, // nodes 10, 11, 12 (indices 0, 1, 2)
            num_nodes: 3,
            local_min_id: 10,
        };

        let (local_edges, local_left, local_right) =
            filter_edges_for_partition(&edges, &partition, 10);

        // Only first 2 edges should be included
        assert_eq!(local_edges.len(), 2);

        // Edges keep their original IDs (not remapped)
        assert_eq!(local_edges[0].from, 10);
        assert_eq!(local_edges[0].to, 11);
        assert_eq!(local_edges[1].from, 11);
        assert_eq!(local_edges[1].to, 12);

        // Edge counts should reflect only internal edges
        assert_eq!(local_left[0], 0);   // node 10: no left edges
        assert_eq!(local_right[0], 1);  // node 10: 1 right edge (to 11)
        assert_eq!(local_left[1], 1);   // node 11: 1 left edge (from 10)
        assert_eq!(local_right[1], 1);  // node 11: 1 right edge (to 12)
        assert_eq!(local_left[2], 1);   // node 12: 1 left edge (from 11)
        assert_eq!(local_right[2], 0);  // node 12: 0 right edges (edge to 15 is external)
    }

    #[test]
    fn test_no_partitioning_with_zero() {
        let partitions = create_partitions(100, 5000, 0);
        assert_eq!(partitions.len(), 1);
        assert_eq!(partitions[0].num_nodes(), 5000);
    }

    #[cfg(feature = "metis")]
    #[test]
    fn test_build_csr() {
        let edges = vec![
            Edge { from: 0, to: 1, from_rev: false, to_rev: false, overlap: 0, sup_reads: 1 },
            Edge { from: 1, to: 2, from_rev: false, to_rev: false, overlap: 0, sup_reads: 1 },
        ];

        let (xadj, adjncy) = build_csr(3, &edges, 0);

        // xadj should have 4 elements (num_nodes + 1)
        assert_eq!(xadj.len(), 4);
        // Node 0 has 1 neighbor (node 1)
        assert_eq!(xadj[1] - xadj[0], 1);
        // Node 1 has 2 neighbors (nodes 0 and 2)
        assert_eq!(xadj[2] - xadj[1], 2);
        // Node 2 has 1 neighbor (node 1)
        assert_eq!(xadj[3] - xadj[2], 1);

        // Total edges in adjncy should be 2 * num_edges (undirected)
        assert_eq!(adjncy.len(), 4);
    }
}
