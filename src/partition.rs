//! Graph partitioning for large-scale ILP solving
//!
//! Splits large graphs into manageable chunks to avoid solver overflow issues.
//! Uses contiguous node ID ranges for simple, dense Vec indexing.

use crate::ilp::IlpNode;
use crate::Edge;

/// A partition of the graph defined by contiguous node ID range
#[derive(Clone, Debug)]
pub struct Partition {
    /// Start index in global node array (inclusive)
    pub global_start: usize,
    /// End index in global node array (exclusive)
    pub global_end: usize,
    /// Number of nodes in this partition
    pub num_nodes: usize,
    /// The min_id for this partition (global min_id + global_start)
    pub local_min_id: usize,
}

impl Partition {
    /// Check if a global node index belongs to this partition
    #[inline]
    pub fn contains(&self, global_idx: usize) -> bool {
        global_idx >= self.global_start && global_idx < self.global_end
    }

    /// Convert global node index to local partition index
    #[inline]
    pub fn to_local(&self, global_idx: usize) -> usize {
        debug_assert!(self.contains(global_idx));
        global_idx - self.global_start
    }
}

/// Create partitions for a graph with the given number of nodes
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
        return vec![Partition {
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
        partitions.push(Partition {
            global_start: start,
            global_end: end,
            num_nodes: end - start,
            local_min_id: min_id + start,
        });
        start = end;
    }

    partitions
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
    let mut local_edges = Vec::new();
    let mut local_left_edges = vec![0usize; partition.num_nodes];
    let mut local_right_edges = vec![0usize; partition.num_nodes];

    for edge in edges {
        let from_global_idx = edge.from - global_min_id;
        let to_global_idx = edge.to - global_min_id;

        // Only include edges where both endpoints are in this partition
        if partition.contains(from_global_idx) && partition.contains(to_global_idx) {
            // Keep original edge IDs - the solver uses local_min_id to compute indices
            local_edges.push(*edge);

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
    (partition.global_start..partition.global_end)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_partitions_no_split() {
        let partitions = create_partitions(100, 500, 1000);
        assert_eq!(partitions.len(), 1);
        assert_eq!(partitions[0].global_start, 0);
        assert_eq!(partitions[0].global_end, 500);
        assert_eq!(partitions[0].local_min_id, 100);
    }

    #[test]
    fn test_create_partitions_split() {
        let partitions = create_partitions(100, 2500, 1000);
        assert_eq!(partitions.len(), 3);

        assert_eq!(partitions[0].global_start, 0);
        assert_eq!(partitions[0].global_end, 1000);
        assert_eq!(partitions[0].local_min_id, 100);  // min_id + 0

        assert_eq!(partitions[1].global_start, 1000);
        assert_eq!(partitions[1].global_end, 2000);
        assert_eq!(partitions[1].local_min_id, 1100);  // min_id + 1000

        assert_eq!(partitions[2].global_start, 2000);
        assert_eq!(partitions[2].global_end, 2500);
        assert_eq!(partitions[2].local_min_id, 2100);  // min_id + 2000
    }

    #[test]
    fn test_partition_contains() {
        let partition = Partition {
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
    fn test_filter_edges_internal() {
        let edges = vec![
            Edge { from: 10, to: 11, from_rev: false, to_rev: false, overlap: 0, sup_reads: 1 },
            Edge { from: 11, to: 12, from_rev: false, to_rev: false, overlap: 0, sup_reads: 2 },
            Edge { from: 12, to: 15, from_rev: false, to_rev: false, overlap: 0, sup_reads: 3 }, // crosses boundary
        ];

        let partition = Partition {
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
        assert_eq!(partitions[0].num_nodes, 5000);
    }
}
