use sprs::TriMat;
use crate::graph::Graph;

pub(crate) fn minimum_degree_ordering(mut graph: Graph) -> Vec<usize> {
    let mut ordering = vec![0; graph.len()];
    let mut current_index = 0usize;
    let nrows = graph.len();
    let mut current_nodes = nrows;
    let mut visited = vec![false; nrows];

    while current_nodes > 0 {
        let mut node_with_minimum_vertex = None;
        let mut current_minimum_degree = usize::MAX;

        // Find out which vertex has the lowest degree.
        for vertex in 0..graph.len() {
            if !visited[vertex] && graph.neighbors(vertex).count() < current_minimum_degree {
                current_minimum_degree = graph.neighbors(vertex).count();
                node_with_minimum_vertex = Some(vertex);
            }
        }
        // Set the lowest degree vertex is the first one to be removed.
        ordering[current_index] = node_with_minimum_vertex.unwrap();
        current_index += 1;
        current_nodes -= 1;
        visited[node_with_minimum_vertex.unwrap()] = true;

        // Get the neighbors of the vertex with minimum degree
        let neighbors_of_vertex_minimum_degree:Vec<usize> = graph.neighbors(node_with_minimum_vertex.unwrap()).into_iter().map(|(node, edge)| node).collect();
        // Create a graph without including the vertex with minimum degree,
        let mut trimat_graph = TriMat::new((nrows, nrows));

        // Once the vertex with minimum degree is removed, all its neigbors get connected
        // to each other, i.e new edges are added.
        for &node in neighbors_of_vertex_minimum_degree.iter() {
            for &another_node in neighbors_of_vertex_minimum_degree.iter() {
                if node != another_node {
                    trimat_graph.add_triplet(node, another_node, 1f64);
                }
            }
        }

        // The connections concerned with the other remaining vertices are added back
        for node in 0..nrows {
            if node != node_with_minimum_vertex.unwrap() {
                for (neighbor, _) in graph.neighbors(node) {
                    if neighbor != node_with_minimum_vertex.unwrap() {
                        trimat_graph.add_triplet(node, neighbor, 1f64);
                    }
                }
            }
        }

        graph.graph_csr = trimat_graph.to_csr();
    }

    ordering
}

#[cfg(test)]
mod tests {
    use crate::algorithms::minimum_degree_order::minimum_degree_ordering;
    use crate::graph::Graph;

    #[test]
    fn test_minimum_degree_ordering() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 1, 2.);
        adjacency.insert(1, 2, 1.);
        adjacency.insert(1, 3, 4.);
        adjacency.insert(2, 3, 2.);
        adjacency.insert(3, 4, 1.);

        adjacency.insert(1, 0, 2.);
        adjacency.insert(2, 1, 1.);
        adjacency.insert(3, 1, 4.);
        adjacency.insert(3, 2, 2.);
        adjacency.insert(4, 3, 1.);

        let ordering = minimum_degree_ordering(adjacency.clone());
        assert_eq!(ordering, [0, 4, 1, 2, 3]);
    }
}