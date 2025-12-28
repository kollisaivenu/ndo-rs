use rustc_hash::FxHashSet;
use crate::graph::Graph;
use crate::imbalance::{calculate_imbalance, compute_parts_load};

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
            // the jetrs is run balance out the two partitions.
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
    // Iterate over all the vertices to find out which vertices in vertex separator partition provides the best gain (decrease in vertex separator weight)
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
    calculate_approximate_gain_and_get_positive_gain_moves(graph, first_filter_eligible_moves, partition, vertex_weights, &dest_partition, &gain)
}

fn calculate_gain(graph: &Graph, partition: &[usize], locked_vertices: &[bool], vertex_weights: &[i64]) -> (Vec<usize>, Vec<Option<i64>>){
    // This function calculates the gain for the vertices in the vertex separator.
    let mut dest_partition = vec![0; graph.len()];
    let mut gain = vec![None; graph.len()];
    for vertex in 0..graph.len() {
        // Consider only vertices present in the vertex separator partition
        if !locked_vertices[vertex]  && partition[vertex] == 2 {
            let mut best_partition = 0;
            let mut gain_of_vertex = i64::MIN;
            // Check which partition is it best for them to move to (partition 0 or partition 1).
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
            // Assign the best destination partition and the gain associated with the move.
            dest_partition[vertex] = best_partition;
            gain[vertex] = Some(gain_of_vertex);
        }
    }

    (dest_partition, gain)
}

fn calculate_approximate_gain_and_get_positive_gain_moves(graph: &Graph, first_filter_eligible_moves: Vec<usize>, partition: &[usize], vertex_weights: &[i64], dest_partition: &[usize], gain: &[Option<i64>]) -> (Vec<Move>) {
    // A hashset is created as it is faster to check if vertex passed the first filter or not.
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

            if is_more_important(neighbor_vertex, vertex, &gain, &first_filter_eligible_vertices) {
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

fn is_more_important(vertex1: usize, vertex2: usize, gain: &[Option<i64>], list_of_vertices: &FxHashSet<usize>) -> bool {
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
    // This function rebalances the two partitions
    let max_slots = 25;
    // Determine which is the heavier partition
    let heavy_partition= get_heavy_partition(partition_weights);
    // Calculate the loss for vertices in the vertex separator partition
    let loss = calculate_loss(graph, partition, heavy_partition, vertex_weights);

    // Slot the loss values into different buckets. This is to prevent sorting the loss values
    // which can be expensive.
    let bucket = place_vertices_in_bucket(partition, &loss, max_slots);

    // Determine which vertices to move from the vertex separator to the light partition so that the heavy partition becomes lighter.
    determine_moves_to_rebalance(partition_weights, heavy_partition, max_slots, &bucket, &loss, vertex_weights, balance_factor)
}

fn calculate_loss(graph: &Graph, partition: &[usize], heavy_partition: usize, vertex_weights: &[i64]) -> Vec<Option<i64>> {
    // This function calculates the loss of the vertices in the vertex separator. Each vertex is
    // only allowed to move to the lighter partition.
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

fn place_vertices_in_bucket(partition: &[usize], loss: &[Option<i64>], max_slots: usize) -> Vec<Vec<usize>> {
    // This function buckets the vertices based on the loss value.
    let mut bucket = init_bucket(1,  max_slots);
    let num_vertices = partition.len();

    for vertex in 0..num_vertices{

        if partition[vertex] == 2 {
            let slot = calculate_slot(loss[vertex].unwrap(), max_slots);
            bucket[get_index_for_bucket(0, slot, max_slots)].push(vertex);
        }
    }

    bucket
}

fn get_heavy_partition(partition_weights: &[i64]) -> usize {
    // This function determines which of the two partitions is the heavier one.
    let heavy_partition:usize;

    if partition_weights[0] > partition_weights[1] {
        heavy_partition = 0;
    } else {
        heavy_partition = 1;
    }

    heavy_partition

}

fn determine_moves_to_rebalance(partition_weights: &[i64], heavy_partition: usize, max_slots: usize, bucket: &Vec<Vec<usize>>, loss: &[Option<i64>], vertex_weights: &[i64], balance_factor: f64) -> Vec<Move> {
    // This function determines which vertices from the vertex separator should move to reduce the
    // weight of the heavy partition
    let mut moves = Vec::new();
    let mut is_still_heavy_partition = true;
    let mut current_weight = 0f64;
    // Determine how much weight from the heavy partition should be removed.
    let excess_weight = partition_weights[heavy_partition] as f64 - (partition_weights[0] as f64 + partition_weights[1] as f64)*balance_factor;

    for slot in 0..max_slots {

        for &vertex in &bucket[get_index_for_bucket(0, slot, max_slots)] {
            // Keep removing the vertices until the excess weight (m_max
            if current_weight < excess_weight {
                // loss[vertex] + vertex_weights[vertex] is the amount of weight that is lost when
                // vertex moves to the light partition
                current_weight = current_weight + ((loss[vertex].unwrap() + vertex_weights[vertex]) as f64);
                moves.push(Move {vertex, partition_id: 1 - heavy_partition});
            } else {
                // This is for early stopping incase the excess weight is removed.
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

#[cfg(test)]
mod tests {
    use crate::algorithms::jet_vertex_separator_refiner::{apply_moves, calculate_approximate_gain_and_get_positive_gain_moves, calculate_gain, calculate_loss, calculate_slot, determine_moves_to_rebalance, gain_conn_ratio_filter, get_heavy_partition, init_bucket, is_more_important, lock_vertices, place_vertices_in_bucket, Move};
    use crate::graph::Graph;

    #[test]
    fn test_calculate_gain() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 1, 2);
        adjacency.insert(1, 2, 1);
        adjacency.insert(3, 4, 2);
        adjacency.insert(3, 5, 1);
        adjacency.insert(1, 3, 3);
        adjacency.insert(1, 4, 3);

        adjacency.insert(1, 0, 2);
        adjacency.insert(2, 1, 1);
        adjacency.insert(4, 3, 2);
        adjacency.insert(5, 3, 1);
        adjacency.insert(3, 1, 3);
        adjacency.insert(4, 1, 3);

        let partition = [0, 2, 0, 2, 1, 1];
        let vertex_weights = [1, 3, 1, 1, 5, 2];
        let locked_vertices = [false; 6];

        let (dest_partition, gain) = calculate_gain(&adjacency, &partition, &locked_vertices, &vertex_weights);

        assert_eq!(dest_partition, vec![0, 1, 0, 1, 0, 0]);
        assert_eq!(gain, vec![None, Some(1), None, Some(1), None, None]);
    }

    #[test]
    fn test_calculate_approximate_gain_and_get_positive_moves() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 3, 2);
        adjacency.insert(1, 4, 1);
        adjacency.insert(2, 4, 2);
        adjacency.insert(3, 5, 1);
        adjacency.insert(4, 6, 3);
        adjacency.insert(3, 4, 2);

        adjacency.insert(3, 0, 2);
        adjacency.insert(4, 1, 1);
        adjacency.insert(4, 2, 2);
        adjacency.insert(5, 3, 1);
        adjacency.insert(6, 4, 3);
        adjacency.insert(4, 3, 2);

        let partition = [0, 0, 0, 2, 2, 1, 1];
        let vertex_weights = [2, 1, 1, 3, 5, 1, 4];
        let locked_vertices = [false; 7];
        let dest_partition = [0, 0, 0, 0, 1, 0, 0];
        let gain = [None, None, None, Some(2), Some(3), None, None];
        let filtered_vertices = vec![3, 4];

        let moves = calculate_approximate_gain_and_get_positive_gain_moves(&adjacency, filtered_vertices.clone(), &partition, &vertex_weights, &dest_partition, &gain);

        assert_eq!(moves[0].vertex, 4);
        assert_eq!(moves[0].partition_id, 1);
    }

    #[test]
    fn test_gain_conn_ratio_filter() {
        let gain = vec![Some(-8), Some(-1), Some(1)];
        let vertex_weights = vec![1, 5, 3];
        let partition = vec![2, 2, 2];
        let locked_vertices = [false; 3];
        let filtered_vertices = gain_conn_ratio_filter(&locked_vertices, &partition, &gain, 0.75, &vertex_weights);
        assert_eq!(filtered_vertices, vec![1, 2]);
    }

    #[test]
    fn test_apply_moves() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 3, 2);
        adjacency.insert(1, 4, 1);
        adjacency.insert(2, 4, 2);
        adjacency.insert(3, 5, 1);
        adjacency.insert(4, 6, 3);

        adjacency.insert(3, 0, 2);
        adjacency.insert(4, 1, 1);
        adjacency.insert(4, 2, 2);
        adjacency.insert(5, 3, 1);
        adjacency.insert(6, 4, 3);

        let mut partition = [0, 0, 0, 2, 2, 1, 1];
        let vertex_weights = [1, 1, 2, 3, 5, 1, 6];
        let mut partition_weights = vec![4, 7, 8];
        let moves = vec![Move {vertex: 4, partition_id: 1}];

        apply_moves(&moves, &mut partition, &adjacency, &mut partition_weights, &vertex_weights);
        assert_eq!(partition, [0, 2, 2, 2, 1, 1, 1]);
        assert_eq!(partition_weights, [1, 12, 6]);
    }

    #[test]
    fn test_calculate_slot() {
        let slot1 = calculate_slot(-4, 3);
        let slot2 = calculate_slot(0, 3);
        let slot3 = calculate_slot(6, 8);
        let slot4 = calculate_slot(10, 3);

        assert_eq!(slot1, 0);
        assert_eq!(slot2, 1);
        assert_eq!(slot3, 4);
        assert_eq!(slot4, 2);
    }

    #[test]
    fn test_is_more_important(){
        let gain = [Some(4), Some(2), Some(2), Some(1)];
        let list_of_vertices = [0, 1, 2].into_iter().collect();
        let result1 = is_more_important(0, 2, &gain, &list_of_vertices);
        assert_eq!(result1, true);

        let result2 = is_more_important(1, 2, &gain, &list_of_vertices);
        assert_eq!(result2, true);

        let result3 = is_more_important(3, 2, &gain, &list_of_vertices);
        assert_eq!(result3, false);
    }

    #[test]
    fn test_calculate_loss() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 1, 2);
        adjacency.insert(1, 2, 1);
        adjacency.insert(3, 4, 2);
        adjacency.insert(3, 5, 1);
        adjacency.insert(1, 3, 3);
        adjacency.insert(1, 4, 3);

        adjacency.insert(1, 0, 2);
        adjacency.insert(2, 1, 1);
        adjacency.insert(4, 3, 2);
        adjacency.insert(5, 3, 1);
        adjacency.insert(3, 1, 3);
        adjacency.insert(4, 1, 3);

        let partition = [0, 2, 0, 2, 1, 1];
        let vertex_weights = [1, 3, 1, 1, 5, 2];
        let locked_vertices = [false; 6];
        let loss = calculate_loss(&adjacency, &partition, 1, &vertex_weights);

        assert_eq!(loss, vec![None, Some(2), None, Some(6), None, None]);

    }

    #[test]
    fn test_get_heavy_partition() {
        let partition_weight = [20, 15, 5];
        let heavy_partition = get_heavy_partition(&partition_weight);
        assert_eq!(heavy_partition, 0);
    }

    #[test]
    fn test_determine_moves_to_rebalance() {
        let mut adjacency = Graph::new();
        adjacency.insert(0, 1, 2);
        adjacency.insert(0, 2, 1);
        adjacency.insert(0, 3, 2);
        adjacency.insert(1, 4, 1);

        adjacency.insert(1, 0, 2);
        adjacency.insert(2, 0, 1);
        adjacency.insert(3, 0, 2);
        adjacency.insert(4, 1, 1);

        let partition = vec![2, 2, 0, 1, 0];
        let vertex_weights = vec![1, 2, 3, 3, 3];
        let partition_weights = vec![6, 3, 3];
        let loss = calculate_loss(&adjacency, &partition, 0, &vertex_weights);

        assert_eq!(loss, vec![Some(2), Some(1), None, None, None]);

        let bucket = place_vertices_in_bucket(&partition, &loss, 25);

        assert_eq!(bucket[2], [1]);
        assert_eq!(bucket[3], [0]);

        let moves = determine_moves_to_rebalance(&partition_weights, 0, 25, &bucket, &loss, &vertex_weights, 0.5);

        assert_eq!(moves.len(), 1);
        assert_eq!(moves[0].vertex, 1);
        assert_eq!(moves[0].partition_id, 1);
    }

    #[test]
    fn test_place_vertices_in_bucket() {
        let partition = vec![0, 0, 2, 2, 1, 1];
        let loss = vec![None, None, Some(-1), Some(2), None, None];
        let bucket = place_vertices_in_bucket(&partition, &loss, 25);
        assert_eq!(bucket[0], vec![2]);
        assert_eq!(bucket[3], vec![3]);
    }

    #[test]
    fn test_lock_vertices() {
        let moves = vec![Move {vertex: 4, partition_id: 1}, Move {vertex: 1, partition_id: 1}];
        let mut locked_vertices = vec![false; 5];
        lock_vertices(&moves, &mut locked_vertices);
        assert_eq!(locked_vertices, vec![false, true, false, false, true]);
    }
}

