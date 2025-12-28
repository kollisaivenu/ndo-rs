use rustc_hash::FxHashSet;
use crate::graph::Graph;

pub(crate) fn minimum_degree_ordering(graph: Graph) -> Vec<usize> {
    let num_of_vertices = graph.len();

    // Create an adjacency list which is a vector of hashsets. This makes it faster to check if
    // an edge exists or not.
    let mut adjacency_list = vec![FxHashSet::default(); num_of_vertices];
    for vertex in 0..num_of_vertices {
        for (neighbor, _) in graph.neighbors(vertex) {
            adjacency_list[neighbor].insert(vertex);
        }
    }

    // Calculate the degree of each vertex.
    let mut degrees = vec![0usize; num_of_vertices];
    for vertex in 0..num_of_vertices {
        degrees[vertex] = adjacency_list[vertex].len();
    }

    let mut ordering = Vec::with_capacity(num_of_vertices);

    // Calculate the ordering
    let mut eliminated = vec![false; num_of_vertices];

    while ordering.len() < num_of_vertices {
        // Find out which vertex has minimum degree currently
        let mut min_degree = usize::MAX;
        let mut min_degree_vertex = None;

        for vertex in 0..graph.len() {
            if degrees[vertex] < min_degree && !eliminated[vertex]{
                min_degree_vertex = Some(vertex);
                min_degree = degrees[vertex];
            }
        }
        // The vertex with the minimum degree will be the one to be eliminated.
        ordering.push(min_degree_vertex.unwrap());
        // This vertex ius now eliminated
        eliminated[min_degree_vertex.unwrap()] = true;

        // Now that the vertex is removed, we need to form a clique to connect it's neighbors.
        let neighbors: Vec<usize> = adjacency_list[min_degree_vertex.unwrap()].iter().cloned().collect();

        for i in 0..neighbors.len() {
            let neighbor1 = neighbors[i];
            // Remove the vertex which has minimum degree.
            adjacency_list[neighbor1].remove(&min_degree_vertex.unwrap());
            degrees[neighbor1] -= 1;

            // Iterate over the neighbors.
            for j in (i + 1)..neighbors.len() {
                let neighbor2 = neighbors[j];
                // If an edge between the two neighbors doesn't exist, we add it.
                if !adjacency_list[neighbor1].contains(&neighbor2) {
                    adjacency_list[neighbor1].insert(neighbor2);
                    adjacency_list[neighbor2].insert(neighbor1);
                    degrees[neighbor1] += 1;
                    degrees[neighbor2] += 1;
                }
            }
        }
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
        adjacency.insert(0, 1, 2);
        adjacency.insert(1, 2, 1);
        adjacency.insert(1, 3, 4);
        adjacency.insert(2, 3, 2);
        adjacency.insert(3, 4, 1);

        adjacency.insert(1, 0, 2);
        adjacency.insert(2, 1, 1);
        adjacency.insert(3, 1, 4);
        adjacency.insert(3, 2, 2);
        adjacency.insert(4, 3, 1);

        let ordering = minimum_degree_ordering(adjacency.clone());
        assert_eq!(ordering, [0, 4, 1, 2, 3]);
    }
}