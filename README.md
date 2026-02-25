# ndo-rs: Rust Library for performing Nested Dissection Ordering for a Sparse Matrix.

## This project is still a WIP. Use the following code snippet to try it out

```rust
use std::path::Path;
use ndo_rs::io::read_matrix_market_as_graph;
use ndo_rs::algorithms::nested_dissection_ordering::NestedDissectionOrdering;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_path = "./testdata/bcsstk15.mtx";
    let graph = read_matrix_market_as_graph(Path::new(file_path))?;
    let (ordering, inverser_ordering) = NestedDissectionOrdering {..Default::default()}.compute_ordering(&graph);
    println!("{:?}", ordering);
    Ok(())
}
```

### At the moment, the results look as follows

| Graph        | METIS (Number of non-zero values) | ndo-rs (Number of non-zero values) |
|--------------|-----------------------------------|------------------------------------|
| bcsstk15.mtx | 507428                            | 604936                             |
| bcsstk16.mtx | 739059                            | 942803                             |
| bcsstk17.mtx | 1108942                           | 1480873                            |
| cfd1.mtx     | 19008078                          | 29872499                           |
| cfd2.mtx     | 36495160                          | 54408449                           |