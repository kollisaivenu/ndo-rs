use crate::graph::Graph;

pub(crate) fn initial_bucket_partitioner(graph: &Graph, vertex_weights: &[i64]) -> Vec<usize>{
    // In the partition vector, each node (identified by index) is assigned 3 indicating that it is
    // not yet assigned to any partition (0 or 1).
    let mut partition = vec![3, graph.len()];
    // Stores the weight of each partition
    let mut partition_weight = [0; 2];

    for node in 0..graph.len() {
        // If the current node is already assigned a partition
        if partition[node] != 3 {
            continue;
        }

        // Else find out which partition is lighter and assign the vertex to that partition
        if partition_weight[0] <= partition_weight[1] {
            partition[node] = 0;
            partition_weight[0] += vertex_weights[node];

        } else if partition_weight[1] < partition_weight[0] {
            partition[node] = 1;
            partition_weight[1] += vertex_weights[node];
        }

        // Also ensure that any neighbouring vertex which belongs to the other partition needs
        // to now be a part of the separator vertex (denoted by partition 2)
        for (neighbor,_) in graph.neighbors(node) {
            if partition[neighbor] != 3 || partition[neighbor] == 1 - partition[node] {
                partition[neighbor] = 2;
            }
        }
    }

    partition
}