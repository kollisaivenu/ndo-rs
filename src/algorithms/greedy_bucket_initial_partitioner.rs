use crate::graph::Graph;

pub(crate) fn initial_bucket_partitioner(graph: &Graph, vertex_weights: &[i64]) -> Vec<usize>{
    // In the partition vector, each node (identified by index) is assigned 3 indicating that it is
    // not yet assigned to any partition (0 or 1).
    let mut partition = vec![3; graph.len()];
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

#[cfg(test)]
mod test {
    use crate::algorithms::greedy_bucket_initial_partitioner::*;
    #[test]
    fn test_greedy_bucket_initial_partitioner() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 1, 2.);
        adjacency.insert(1, 2, 1.);
        adjacency.insert(3, 4, 2.);
        adjacency.insert(3, 5, 1.);
        adjacency.insert(1, 3, 3.);
        adjacency.insert(1, 4, 3.);

        adjacency.insert(1, 0, 2.);
        adjacency.insert(2, 1, 1.);
        adjacency.insert(4, 3, 2.);
        adjacency.insert(5, 3, 1.);
        adjacency.insert(3, 1, 3.);
        adjacency.insert(4, 1, 3.);
        let vertex_weights = vec![1, 2, 3, 4, 5, 6];

        let partition = initial_bucket_partitioner(&adjacency, &vertex_weights);
        assert_eq!(partition, vec![2, 2, 0, 2, 0, 1]);
    }
}
