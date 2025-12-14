use rustc_hash::FxHashSet;
use crate::graph::Graph;
use crate::imbalance::compute_parts_load;

#[derive(Debug)]
struct Move {
    // Struct to store data about a move that can either lead to better vertex separator or
    // re-balance the weights.

    //The index of the vertex.
    vertex: usize,

    // The partition ID of the partition where the vertex should move to.
    partition_id: usize
}

pub(crate) fn jet_vertex_separator_refiner(
    partition: &mut [usize],
    vertex_weights: &[i64],
    adjacency: Graph,
    iterations: u32,
    balance_factor: f64,
    filter_ratio: f64,
    tolerance_factor: f64) {

    debug_assert!(!partition.is_empty());
    debug_assert_eq!(partition.len(), vertex_weights.len());
    debug_assert_eq!(partition.len(), adjacency.len());

    let mut partition_iter = partition.to_vec();
    let mut current_iteration = 0;
    let mut locked_vertices = vec![false; adjacency.len()];
    let num_of_partitions = 3;
    let mut partition_weights = compute_parts_load(&partition, num_of_partitions, &vertex_weights);
    let mut imbalance_of_current_iter_partition = calculate_imbalance(&partition_weights);
    let mut imbalance_of_best_partition = imbalance_of_current_iter_partition;
    let mut best_vertex_separator_weight = get_weight_of_vertex_separator(&partition_weights);

    while current_iteration < iterations {
        let mut moves = Vec::new();
        if imbalance_of_current_iter_partition < balance_factor {
            // the jetlp subroutine is used to generate a better partition
            moves = jetlp(&adjacency,
                             &partition_iter,
                             &locked_vertices,
                             filter_ratio,
                             vertex_weights);

            // Based on the suggested moves of the jetlp subroutine, the vertices to be moved are locked
            // to ensure that they don't become eligible to move in the next iteration.
            // This prevents oscillation of vertices
            lock_vertices(&moves, &mut locked_vertices);

        } else {
            moves = jetrw(&adjacency,
                             &partition_iter,
                             vertex_weights,
                             balance_factor,
                             &partition_weights);
        }

        apply_moves(&moves, &mut partition_iter, &adjacency, &mut partition_weights, vertex_weights);

        imbalance_of_current_iter_partition = calculate_imbalance(&partition_weights);

        let current_vertex_separator_weight = get_weight_of_vertex_separator(&partition_weights);
        // Check if the current iteration partition is balance
        if imbalance_of_current_iter_partition < balance_factor {
            // Check if the current iteration partition is better than the current best partition
            if current_vertex_separator_weight < best_vertex_separator_weight  {
                // Current iteration partition is chosen as the best partition
                if current_vertex_separator_weight < (tolerance_factor*(best_vertex_separator_weight  as f64)).floor() as i64 {
                    current_iteration = 0;
                }
                partition.copy_from_slice(&partition_iter);
                imbalance_of_best_partition = imbalance_of_current_iter_partition;
                best_vertex_separator_weight = current_vertex_separator_weight;
            } else {
                current_iteration += 1;
            }
        } else if imbalance_of_current_iter_partition < imbalance_of_best_partition {
            // Current iteration is better balanced than the best iteration, hence this is made
            // the best iteration
            partition.copy_from_slice(&partition_iter);
            imbalance_of_best_partition = imbalance_of_current_iter_partition;
            best_vertex_separator_weight = current_vertex_separator_weight;
            current_iteration = 0
        } else {
            current_iteration += 1;
        }
    }
}

fn jetlp(graph: &Graph, partition: &[usize], locked_vertices: &[bool], filter_ratio: f64, vertex_weights: &[i64]) -> Vec<Move> {
    // iterate over all the vertices to find out which vertices provides the best gain (decrease in vertex separator weight)

    let (dest_partition, gain) = calculate_gain(graph, partition, locked_vertices, vertex_weights);

    // First filter is applied to check which of the vertices are eligible for moving from one partition
    // to another. Either the gain should be positive or can be slightly negative (based on the filter ratio).
    // Slightly negative gain vertices are also considered in the hope that they could provide better global solutions
    let first_filter_eligible_moves = gain_conn_ratio_filter(
        locked_vertices,
        partition,
        &gain,
        filter_ratio,
        vertex_weights);

    // A heuristic attempt is made to approximate the true gain that would occur since
    // two positive moves when applied simultaneously can be detrimental.
    calculate_approximate_gain_and_get_positive_moves(graph, first_filter_eligible_moves, partition, vertex_weights, &dest_partition, &gain)
}

fn calculate_gain(graph: &Graph, partition: &[usize], locked_vertices: &[bool], vertex_weights: &[i64]) -> (Vec<usize>, Vec<Option<i64>>){
    let mut dest_partition = vec![0; graph.len()];
    let mut gain = vec![None; graph.len()];
    for vertex in 0..graph.len() {
        // These are values if all the vertex belongs to the same partition as its neighbours (not a boundary vertex).

        if !locked_vertices[vertex]  && partition[vertex] == 2 {
            let mut best_partition = 0;
            let mut gain_of_vertex = i64::MIN;

            for possible_dest in 0..2 {
                let other_dest = 1 - possible_dest;
                let mut strength = 0i64;
                for (neighbor, _) in graph.neighbors(vertex) {
                    if partition[neighbor] == other_dest {
                        strength += vertex_weights[neighbor];
                    }
                }

                if vertex_weights[vertex] - strength > gain_of_vertex {
                    gain_of_vertex = vertex_weights[vertex] - strength;
                    best_partition = possible_dest;
                }

            }

            dest_partition[vertex] = best_partition;
            gain[vertex] = Some(gain_of_vertex);
        }
    }

    (dest_partition, gain)
}

fn calculate_approximate_gain_and_get_positive_moves(graph: &Graph, first_filter_eligible_moves: Vec<usize>, partition: &[usize], vertex_weights: &[i64], dest_partition: &[usize], gain: &[Option<i64>]) -> (Vec<Move>) {
    // A hashset is created as it is faster to check which vertex is eligible to move.
    let first_filter_eligible_vertices: FxHashSet<usize> = first_filter_eligible_moves.clone().into_iter().collect();

    // A heuristic attempt is made to approximate the true gain that would occur since
    // two positive moves when applied simultaneously can be detrimental.
    let mut moves = Vec::with_capacity(first_filter_eligible_moves.len());

    for vertex_index in 0..first_filter_eligible_moves.len() {
        let vertex = first_filter_eligible_moves[vertex_index];

        let mut gain_for_vertex = vertex_weights[vertex];

        for (neighbor_vertex, _) in graph.neighbors(vertex) {
            let dest_partition_of_vertex = dest_partition[vertex];
            let other_dest = 1 - dest_partition_of_vertex;

            if is_higher_placed(neighbor_vertex, vertex, &gain, &first_filter_eligible_vertices) {
                // This is the scenario where a vertex 'x' is moving to partition '0' but one of its
                // neighbor 'y' is moving to partition '1'. This neighbor 'y' can adversely impact
                // the gain of the vertex 'x' as it might have to move vertices into the vertex separator.
                if dest_partition[neighbor_vertex] == other_dest {
                    gain_for_vertex -= vertex_weights[neighbor_vertex];
                }
            } else if partition[neighbor_vertex] == other_dest {
                gain_for_vertex -= vertex_weights[neighbor_vertex];
            }
        }
        if gain_for_vertex > 0 {
            moves.push(Move {vertex: vertex, partition_id: dest_partition[vertex]});
        }
    }

    moves
}

fn lock_vertices(moves: &Vec<Move>, locked_vertices: &mut [bool]) {
    // This function locks the vertices that shouldn't be moved in the subsequent iterations.
    locked_vertices.fill(false);

    for single_move in moves{
        locked_vertices[single_move.vertex] = true;
    }

}

fn gain_conn_ratio_filter(locked_vertices: &[bool], partitions: &[usize], gain: &[Option<i64>], filter_ratio: f64, vertex_weights: &[i64]) -> Vec<usize> {
    // Get a list of vertices that have a positive gain or slightly negative gain value (based on the filter ratio).

    let num_vertices = partitions.len();
    let mut list_of_moveable_vertices  = Vec::new();

    for vertex in 0..num_vertices {
        if (!locked_vertices[vertex])
            &&
            !gain[vertex].is_none()
            &&
            (gain[vertex].unwrap() > 0 || -gain[vertex].unwrap() < (filter_ratio * (vertex_weights[vertex] as f64)).floor() as i64){
            list_of_moveable_vertices.push(vertex);
        }
    }

    list_of_moveable_vertices
}

fn is_higher_placed(vertex1: usize, vertex2: usize, gain: &[Option<i64>], list_of_vertices: &FxHashSet<usize>) -> bool {
    // Checks if vertex1 is better ranked than vertex2 (used in the vertex afterburner). Better ranked
    // indicates that it has a higher gain than vertex2 (and hence should be given more preference)

    if list_of_vertices.contains(&vertex1) && !gain[vertex1].is_none()
        &&
        !gain[vertex2].is_none()
        && ((gain[vertex1].unwrap() > gain[vertex2].unwrap())
        ||
        (gain[vertex1].unwrap() == gain[vertex2].unwrap() && vertex1 < vertex2)){
        return true;
    }

    false
}

fn jetrw(graph: &Graph, partition: &[usize], vertex_weights: &[i64], balance_factor: f64, partition_weights: &[i64]) -> Vec<Move> {
    // Weaker but better rebalancer in terms of the change in edgecut
    let max_slots = 25;

    let heavy_partition= get_heavy_partition(partition_weights);

    let loss = calculate_loss(graph, partition, heavy_partition, vertex_weights);

    // Slot the loss values into different buckets. This is to prevent sorting the loss values
    // which can be expensive.
    let bucket = place_vertices_in_bucket(graph.len(), partition, &loss, max_slots);

    // For each of the heavy partitions, decide the vertices that can be moved from the
    // heavy partitions such that the increase in edge cut is minimized.
    let moves = calculate_moves_to_rebalance(partition_weights, heavy_partition, max_slots, &bucket, &loss, vertex_weights, balance_factor);
    moves
}

fn calculate_loss(graph: &Graph, partition: &[usize], heavy_partition: usize, vertex_weights: &[i64]) -> Vec<Option<i64>> {
    let (loss): (Vec<Option<i64>>) = (0..graph.len()).map(|vertex| {
        let mut calculated_loss_for_vertex = None;

        if partition[vertex] == 2 {
            let mut total_weight = 0i64;
            for (neighbor, _) in graph.neighbors(vertex) {
                if partition[neighbor] == heavy_partition {
                    total_weight += vertex_weights[neighbor];
                }
            }
            calculated_loss_for_vertex = Some(total_weight - vertex_weights[vertex]);
        }
        calculated_loss_for_vertex
    }).collect();

    loss
}

fn place_vertices_in_bucket(num_of_vertices: usize, partition: &[usize], loss: &[Option<i64>], max_slots: usize) -> Vec<Vec<usize>> {
    let mut bucket = init_bucket(1,  max_slots);

    for vertex in 0..num_of_vertices{

        if partition[vertex] == 2 {
            let slot = calculate_slot(loss[vertex].unwrap(), max_slots);
            bucket[get_index_for_bucket(0, slot, max_slots)].push(vertex);
        }
    }

    bucket
}

fn get_heavy_partition(partition_weights: &[i64]) -> usize {
    let heavy_partition:usize;

    if partition_weights[0] > partition_weights[1] {
        heavy_partition = 0;
    } else {
        heavy_partition = 1;
    }

    heavy_partition

}

fn calculate_moves_to_rebalance(partition_weights: &[i64], heavy_partition: usize, max_slots: usize, bucket: &Vec<Vec<usize>>, loss: &[Option<i64>], vertex_weights: &[i64], balance_factor: f64) -> Vec<Move> {
    let mut moves = Vec::new();
    let mut is_still_heavy_partition = true;
    let mut m = 0f64;
    let m_max = partition_weights[heavy_partition] as f64 - (partition_weights[0] as f64 + partition_weights[1] as f64)*balance_factor;

    for slot in 0..max_slots {

        for &vertex in &bucket[get_index_for_bucket(0, slot, max_slots)] {

            if m < m_max {
                m = m + ((loss[vertex].unwrap() + vertex_weights[vertex]) as f64);
                moves.push(Move {vertex, partition_id: 1 - heavy_partition});
            } else {
                is_still_heavy_partition = false;
                break;
            }
        }

        if !is_still_heavy_partition {
            break;
        }
    }

    moves

}

fn calculate_slot(loss: i64, max_slot_size: usize) -> usize {
    // Calculate the slot in which the vertex should be put in based on the loss value.

    if loss < 0 {
        0
    } else if loss == 0 {
        1
    } else {
        ((2 + loss.ilog2()) as usize).min(max_slot_size-1)
    }
}

fn init_bucket(num_heavy_partitions: usize, max_slots: usize) -> Vec<Vec<usize>>{
    // Initialize the bucket where a list of vertices are stored in slots.

    let rows = num_heavy_partitions*max_slots;
    let mut bucket: Vec<Vec<usize>> = Vec::with_capacity(rows);

    for _ in 0..rows {
        bucket.push(Vec::new());
    }

    bucket
}

fn get_index_for_bucket(partition_index: usize, slot: usize, max_slots: usize) -> usize {
    // Gets the index of the bucket based on slot and partition index.

    partition_index * max_slots + slot
}







fn calculate_imbalance(partition_weights: &[i64]) -> f64 {
    // This function calculates the imbalance of the two partitions (0 and 1).
    let total_weight= (partition_weights[0] + partition_weights[1]) as f64;
    let weight_of_heavy_partition = partition_weights[0].max(partition_weights[1]) as f64;

    weight_of_heavy_partition/total_weight
}

fn get_weight_of_vertex_separator(partition_weights: &[i64]) -> i64 {
    // This function retrieves the weight of the vertex separator (denoted by partition 2)
    partition_weights[2]
}

fn apply_moves(moves: &Vec<Move>, partition: &mut [usize], graph: &Graph, partition_weights: &mut Vec<i64>, vertex_weights: &[i64]) {
    // Once the moves are obtained from jetlp/jetrw, the moves are applied.
    for single_move in moves {
        let vertex_id = single_move.vertex;
        let destination_partition = single_move.partition_id;
        // Change the partition of the vertex
        partition[vertex_id] = destination_partition;
        // Decrease the weight of the vertex separator
        partition_weights[2] -= vertex_weights[vertex_id];
        // Increase the weight of the destination partition
        partition_weights[destination_partition] += vertex_weights[vertex_id];

        // Every neighbor of the vertex which is not in the destination partition or the vertex
        // separator move into the vertex separator now
        for (neighbor, _) in graph.neighbors(vertex_id) {
            if partition[neighbor] == 1 - destination_partition {
                // Change the partition of the neighbor
                partition[neighbor] = 2;
                // Decrease the weight of the partition that it belongs to
                partition_weights[1 - destination_partition] -= vertex_weights[neighbor];
                // Increase the weight of the vertex separator partition
                partition_weights[2] += vertex_weights[neighbor];
            }
        }
    }
}

