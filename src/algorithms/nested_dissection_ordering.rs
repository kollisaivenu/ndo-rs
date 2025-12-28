use sprs::TriMat;
use crate::algorithms::minimum_degree_order::minimum_degree_ordering;
use crate::algorithms::multilevel_vertex_separator::multilevel_vertex_separator;
use crate::graph::Graph;

fn nested_dissection_ordering(graph: &Graph, nested_dissection_ordering_struct: &NestedDissectionOrdering) -> Vec<usize> {
    // If the graph is small, we run the minimum degree ordering algorithm to find the vertex order.
    if graph.len() <= 50 {
        minimum_degree_ordering(graph.clone());
    }

    // Else we find the vertex separator using multilevel vertex separator algorithm
    // In the partition 0 & 1 indicate the two partitions and 2 represents the vertices in
    // the vertex separator
    let partition = find_separator(graph.clone(),
                                               nested_dissection_ordering_struct.jet_iterations,
                                               nested_dissection_ordering_struct.balance_factor,
                                               nested_dissection_ordering_struct.jet_filter_ratio,
                                               nested_dissection_ordering_struct.jet_tolerance_factor);

    // Create a subgraph out of vertices belonging to partition 0.
    let (subgraph_a, new_vertex_to_old_vertex_mapping_a) = create_subgraph(graph, &partition, 0);
    // Create a subgraph out of vertices belonging to partition 1.
    let (subgraph_b, new_vertex_to_old_vertex_mapping_b) = create_subgraph(graph, &partition, 1);

    // Perform nested dissection ordering for sub graph A (partition 0)
    let order_for_subgraph_a = nested_dissection_ordering(&subgraph_a, nested_dissection_ordering_struct);
    // Perform nested dissection ordering for sub graph B (partition 1)
    let order_for_subgraph_b = nested_dissection_ordering(&subgraph_b, nested_dissection_ordering_struct);

    // Combine the ordering for sub graph A, sub graph B and vertex separator.
    combine_permutation(&order_for_subgraph_a, &order_for_subgraph_b, &new_vertex_to_old_vertex_mapping_a, &new_vertex_to_old_vertex_mapping_b, &partition)
}

fn create_subgraph(graph: &Graph, partition: &Vec<usize>, partition_num: usize) -> (Graph, Vec<usize>) {
    // Find out number of vertices in the partition
    let mut nrows = 0;
    for &part in partition.iter() {
        if part == partition_num {
            nrows += 1;
        }
    }

    // Assign new vertice numbers to the vertices of the sub graph.
    let mut old_to_new_vertex_map = vec![0; graph.len()];
    let mut new_to_old_vertex_map = vec![0; nrows];

    let mut new_vertex = 0usize;

    for vertex in 0..partition.len() {
        if partition[vertex] == partition_num {
            old_to_new_vertex_map[vertex] = new_vertex;
            new_to_old_vertex_map[new_vertex] = vertex;
            new_vertex += 1;
        }
    }

    // Build the new subgraph
    let mut trimat_graph = TriMat::new((nrows, nrows));

    for vertex in 0..graph.len() {
        for (neighbour, edge) in graph.neighbors(vertex) {
            if partition[vertex] == partition[neighbour] && partition[vertex] == partition_num {
                trimat_graph.add_triplet(old_to_new_vertex_map[vertex], old_to_new_vertex_map[neighbour], edge);
            }
        }
    }

    let mut subgraph = Graph::new();
    subgraph.set_matrix(trimat_graph.to_csr());

    // return the new subgraph and the vertex mapping
    (subgraph, new_to_old_vertex_map)
}

fn find_separator(graph: Graph, jet_iterations: u32, balance_factor: f64, jet_filter_ratio: f64, jet_tolerance_factor: f64) -> Vec<usize> {
    // Find the vertex separator
    let mut partition = vec![0; graph.len()];
    let vertex_weights = vec![1; graph.len()];
    multilevel_vertex_separator(&mut partition, &vertex_weights, graph, None, jet_iterations, balance_factor, jet_filter_ratio, jet_tolerance_factor);
    partition
}

fn combine_permutation(ordering_for_subgraph_a: &Vec<usize>, ordering_for_subgraph_b: &Vec<usize>, new_vertex_to_old_vertex_mapping_subgraph_a: &Vec<usize>, new_vertex_to_old_vertex_mapping_subgraph_b: &Vec<usize>, partition: &Vec<usize>) -> Vec<usize> {
    // Combine the ordering of subgraph A, subgraph B and the vertex separator.
    let mut ordering_for_graph = vec![0; partition.len()];
    let mut ordering_index = 0;

    // Find out ordering for sub graph A
    for index in 0..ordering_for_subgraph_a.len() {
        ordering_for_graph[ordering_index] = new_vertex_to_old_vertex_mapping_subgraph_a[ordering_for_subgraph_a[index]];
        ordering_index += 1;
    }

    // Find out ordering for sub graph B
    for index in 0..ordering_for_subgraph_b.len() {
        ordering_for_graph[ordering_index] = new_vertex_to_old_vertex_mapping_subgraph_b[ordering_for_subgraph_b[index]];
        ordering_index += 1;
    }

    // Find the ordering for vertex separator
    for i in 0..partition.len() {
        if partition[i] == 2 {
            ordering_for_graph[ordering_index] = i;
            ordering_index += 1;
        }
    }

    ordering_for_graph
}

fn get_inverse_ordering(ordering: &Vec<usize>) -> Vec<usize> {
    // Determine the inverse ordering
    let mut inverse_ordering = vec![0; ordering.len()];
    for (i, &val) in ordering.iter().enumerate() {
        inverse_ordering[val] = i;
    }

    inverse_ordering
}

pub struct NestedDissectionOrdering {
    pub jet_iterations: u32,
    pub jet_filter_ratio: f64,
    pub balance_factor: f64,
    pub jet_tolerance_factor: f64,
}
impl Default for NestedDissectionOrdering {
    fn default() -> Self {
        Self {
            jet_iterations: 12,
            jet_filter_ratio: 0.75,
            balance_factor: 0.6,
            jet_tolerance_factor: 0.999
        }
    }

}
impl NestedDissectionOrdering {
    pub fn compute_ordering(&self, graph: &Graph) -> (Vec<usize>, Vec<usize>) {
        let ordering = nested_dissection_ordering(graph, &self);
        let inverse_ordering = get_inverse_ordering(&ordering);
        (ordering, inverse_ordering)
    }
}

#[cfg(test)]
mod tests {
    use crate::algorithms::nested_dissection_ordering::{combine_permutation, create_subgraph, get_inverse_ordering};
    use crate::graph::Graph;

    #[test]
    fn test_create_graph() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 3, 2);
        adjacency.insert(1, 3, 1);
        adjacency.insert(2, 3, 4);
        adjacency.insert(2, 5, 2);
        adjacency.insert(5, 4, 1);
        adjacency.insert(5, 6, 1);

        adjacency.insert(3, 0, 2);
        adjacency.insert(3, 1, 1);
        adjacency.insert(3, 2, 4);
        adjacency.insert(5, 2, 2);
        adjacency.insert(4, 5, 1);
        adjacency.insert(6, 5, 1);

        let partition = vec![0, 0, 2, 0, 1, 1, 1];

        let (subgraph_a, order_a) = create_subgraph(&adjacency, &partition, 0);

        assert_eq!(order_a, [0, 1, 3]);
        assert_eq!(subgraph_a.neighbors(0).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [2]);
        assert_eq!(subgraph_a.neighbors(1).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [2]);
        assert_eq!(subgraph_a.neighbors(2).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [0, 1]);

        let (subgraph_b, order_b) = create_subgraph(&adjacency, &partition, 1);

        assert_eq!(order_b, [4, 5, 6]);
        assert_eq!(subgraph_b.neighbors(0).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [1]);
        assert_eq!(subgraph_b.neighbors(1).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [0, 2]);
        assert_eq!(subgraph_b.neighbors(2).into_iter().map(|(neighbor, edge)| neighbor).collect::<Vec<usize>>(), [1]);
    }

    #[test]
    fn test_combine_permutation() {
        let order_subgraph_a = vec![1, 0, 2];
        let order_subgraph_b = vec![2, 1, 0];
        let new_vertex_to_old_vertex_subgraph_a = vec![0usize, 1, 2];
        let new_vertex_to_old_vertex_subgraph_b = vec![4usize, 5, 6];
        let partition = vec![0, 0, 0, 2, 1, 1, 1];
        let final_order = combine_permutation(&order_subgraph_a, &order_subgraph_b, &new_vertex_to_old_vertex_subgraph_a, &new_vertex_to_old_vertex_subgraph_b, &partition);
        assert_eq!(final_order, [1, 0, 2, 6, 5, 4, 3]);
    }

    #[test]
    fn test_get_inverse_ordering() {
        let ordering = vec![3, 2, 4, 1, 0];
        let inverse_ordering = get_inverse_ordering(&ordering);
        assert_eq!(inverse_ordering, vec![4, 3, 1, 0, 2]);
    }
}