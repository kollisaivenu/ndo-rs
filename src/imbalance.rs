/// Calculates the total weight for each part of a given partition.
pub fn compute_parts_load(partition: &[usize], num_parts: usize, weights: &[i64]) -> Vec<i64> {
    let mut loads = vec![0; num_parts];

    for (&part, w) in partition.iter().zip(weights.into_iter()) {
        if part < num_parts {
            loads[part] += w;
        }
    }

    loads
}

/// Calculates imbalance of the partitions
pub fn calculate_imbalance(partition_weights: &[i64]) -> f64 {
    // This function calculates the imbalance of the two partitions (0 and 1).
    let total_weight= (partition_weights[0] + partition_weights[1]) as f64;
    let weight_of_heavy_partition = partition_weights[0].max(partition_weights[1]) as f64;

    weight_of_heavy_partition/total_weight
}