use crate::algorithms::minimum_degree_order::minimum_degree_ordering;
use crate::algorithms::multilevel_vertex_separator::multilevel_vertex_separator;
use crate::graph::Graph;

fn nested_dissection_ordering(graph: &Graph) -> Vec<usize> {
    // If the graph is small, we run the minimum degree ordering algorithm to find the vertex order.
    if graph.len() <= 50 {
        return minimum_degree_ordering(graph.clone());
    }

    // Else we find the vertex separator using multilevel vertex separator algorithm
    // In the partition 0 & 1 indicate the two partitions and 2 represents the vertices in
    // the vertex separator
    let partition = find_separator(graph.clone());

    // Create a subgraph out of vertices belonging to partition 0.
    let (subgraph_a, new_vertex_to_old_vertex_mapping_a) = create_subgraph(graph, &partition, 0);
    // Create a subgraph out of vertices belonging to partition 1.
    let (subgraph_b, new_vertex_to_old_vertex_mapping_b) = create_subgraph(graph, &partition, 1);

    // Perform nested dissection ordering for sub graph A (partition 0)
    let order_for_subgraph_a = nested_dissection_ordering(&subgraph_a);
    // Perform nested dissection ordering for sub graph B (partition 1)
    let order_for_subgraph_b = nested_dissection_ordering(&subgraph_b);

    // Combine the ordering for sub graph A, sub graph B and vertex separator.
    combine_permutation(&order_for_subgraph_a, &order_for_subgraph_b, &new_vertex_to_old_vertex_mapping_a, &new_vertex_to_old_vertex_mapping_b, &partition)
}

fn create_subgraph(graph: &Graph, partition: &Vec<usize>, partition_num: usize) -> (Graph, Vec<Option<usize>>) {
    // Find out number of vertices in the partition
    let mut nrows = 0;
    for &part in partition.iter() {
        if part == partition_num {
            nrows += 1;
        }
    }

    // Assign new vertice numbers to the vertices of the sub graph.
    let mut old_to_new_vertex_map = vec![None; graph.len()];
    let mut new_to_old_vertex_map = vec![None; nrows];

    let mut new_vertex = 0usize;

    for vertex in 0..partition.len() {
        if partition[vertex] == partition_num {
            old_to_new_vertex_map[vertex] = Some(new_vertex);
            new_to_old_vertex_map[new_vertex] = Some(vertex);
            new_vertex += 1;
        }
    }

    // Build the new subgraph
    let mut trimat_graph = TriMat::new((nrows, nrows));

    for vertex in 0..graph.len() {
        for (neighbour, edge) in graph.neighbors(vertex) {
            if partition[vertex] == partition[neighbour] && partition[vertex] == partition_num {
                trimat_graph.add_triplet(old_to_new_vertex_map[vertex].unwrap(), old_to_new_vertex_map[neighbour].unwrap(), edge);
            }
        }
    }

    let mut subgraph = Graph::new();
    subgraph.graph_csr = trimat_graph.to_csr();

    // return the new subgraph and the vertex mapping
    (subgraph, new_to_old_vertex_map)
}

fn find_separator(graph: Graph) -> Vec<usize> {
    // Find the vertex separator
    let mut partition = vec![0; graph.len()];
    let vertex_weights = vec![1; graph.len()];
    multilevel_vertex_separator(&mut partition, &vertex_weights, graph, None, 12, 0.5, 0.75, 0.999);
    partition
}

fn combine_permutation(ordering_for_subgraph_a: &Vec<usize>, ordering_for_subgraph_b: &Vec<usize>, new_vertex_to_old_vertex_mapping_subgraph_a: &Vec<Option<usize>>, new_vertex_to_old_vertex_mapping_subgraph_b: &Vec<Option<usize>>, partition: &Vec<usize>) -> Vec<usize> {
    // Combine the ordering of subgraph A, subgraph B and the vertex separator.
    let mut ordering_for_graph = vec![0; partition.len()];
    let mut ordering_index = 0;

    // Find out ordering for sub graph A
    for index in 0..ordering_for_subgraph_a.len() {
        ordering_for_graph[ordering_index] = new_vertex_to_old_vertex_mapping_subgraph_a[ordering_for_subgraph_a[index]].unwrap();
        ordering_index += 1;
    }

    // Find out ordering for sub graph B
    for index in 0..ordering_for_subgraph_b.len() {
        ordering_for_graph[ordering_index] = new_vertex_to_old_vertex_mapping_subgraph_b[ordering_for_subgraph_b[index]].unwrap();
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