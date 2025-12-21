/// Calculates the total weight for each part of a given partition.
pub fn compute_parts_load(partition: &[usize], num_parts: usize, weights: &[i64]) -> Vec<i64> {
    // This function computes the weights of partition 0, 1 and 2 (vertex separator).
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

#[cfg(test)]
mod tests {
    use approx::assert_ulps_eq;
    use itertools::assert_equal;
    use crate::imbalance::{compute_parts_load, calculate_imbalance};

    #[test]
    fn test_compute_parts_load() {
        let partition = [0, 0, 1, 1];
        let vtx_weights = vec![4, 7, 5, 2];
        let num_parts = 2;

        let partition_weights = compute_parts_load(&partition, num_parts, &vtx_weights);

        assert_equal(partition_weights, [11, 7]);
    }

    #[test]
    fn test_calculate_imbalance() {
        let partition_weights = vec![3, 3, 1];
        let imb = calculate_imbalance(&partition_weights);
        assert_ulps_eq!(imb, 0.5);
    }
}