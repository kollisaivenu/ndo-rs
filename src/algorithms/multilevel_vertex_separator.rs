use rand::rngs::SmallRng;
use rand::SeedableRng;
use crate::algorithms::greedy_bucket_initial_partitioner::initial_bucket_partitioner;
use crate::algorithms::heavy_edge_matching::heavy_edge_matching_coarse;
use crate::algorithms::jet_vertex_separator_refiner::jet_vertex_separator_refiner;
use crate::graph::Graph;

pub(crate) fn multilevel_vertex_separator(
    partition: &mut [usize],
    vertex_weights: &[i64],
    graph: Graph,
    seed: Option<u64>,
    jet_iterations: u32,
    balance_factor: f64,
    jet_filter_ratio: f64,
    jet_tolerance_factor: f64, ) {

    let mut coarse_graphs = Vec::new();
    coarse_graphs.push(graph.clone());
    let mut fine_vertex_to_coarse_vertex_mappings = Vec::new();
    let mut weights_coarse_graphs = Vec::new();
    weights_coarse_graphs.push(vertex_weights.to_vec());

    let mut rng = match seed {
        Some(seed) => SmallRng::seed_from_u64(seed),
        None => SmallRng::from_entropy()
    };

    // Keep coarsening the graph until the graph has less than 100 nodes
    while coarse_graphs.last().unwrap().len() > 50  {

        let (coarse_graph, fine_vertex_to_coarse_vertex_mapping, weights_of_coarse_graph) = heavy_edge_matching_coarse(coarse_graphs.last().unwrap(), &mut rng, weights_coarse_graphs.last().unwrap());
        // Store the coarse graphs at every level
        coarse_graphs.push(coarse_graph);
        // Store the node weights of every coarse graph at each level
        weights_coarse_graphs.push(weights_of_coarse_graph);
        // Store the vertex mapping (coarse node to finer nodes) of the coarse graph at each level.
        fine_vertex_to_coarse_vertex_mappings.push(fine_vertex_to_coarse_vertex_mapping);

    }

    //let mut coarse_graph_partition = get_initial_vertex_separator(coarse_graphs.last().unwrap().len());
    let mut coarse_graph_partition = initial_bucket_partitioner(&coarse_graphs.last().unwrap(), &weights_coarse_graphs.last().unwrap());

    let mut index = coarse_graphs.len() - 1;

    while index >= 0 {
        // Refine using jet vertex separator refiner
        jet_vertex_separator_refiner(
            &mut coarse_graph_partition,
            &weights_coarse_graphs[index],
            coarse_graphs[index].clone(),
            jet_iterations,
            balance_factor,
            jet_filter_ratio,
            jet_tolerance_factor,
        );

        // Uncoarsen the graph till we reach the initial graph.
        if index > 0 {
            coarse_graph_partition = partition_uncoarse(&coarse_graph_partition, &fine_vertex_to_coarse_vertex_mappings[index-1]);
        } else {
            break;
        }

        index -= 1;
    }

    // Copy over the final partition to the partition array which is passed as input.
    partition.copy_from_slice(&coarse_graph_partition);
}

// Refines the partition from a coarse graph back to the original finer graph.
fn partition_uncoarse(partition: &[usize], fine_vertex_to_coarse_vertex_mapping: &Vec<usize>) -> Vec<usize>{
    // Calculate the number of vertices in the uncoarsed graph (1 up level)


    // Create a partition array for the uncoarsed graph (1 level up)
    // If vertex 1 and 2 of the uncoarsed graph were merged into vertex 0 in the coarsed graph
    // and it belonged to partition 0, then vertex 1 and 2 would belong to partition 0 in the uncoarsed graph.
    let mut new_partition: Vec<usize> = vec![0; fine_vertex_to_coarse_vertex_mapping.len()];

    for vertex in 0..fine_vertex_to_coarse_vertex_mapping.len(){
        new_partition[vertex] = partition[fine_vertex_to_coarse_vertex_mapping[vertex]];
    }

    new_partition
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_uncoarse() {
        let fine_vertex_to_coarse_vertex_mapping = vec![0, 2, 1, 0];
        let weights_coarse_graph = [5, 7, 6];
        let coarse_graph_partition = [1, 0, 0];
        let weights_uncoarse_graph = [2, 6, 7, 3];

        let uncoarsed_graph_partition = partition_uncoarse(&coarse_graph_partition, &fine_vertex_to_coarse_vertex_mapping);

        assert_eq!(uncoarsed_graph_partition, vec![1, 0, 0, 1]);
    }
}