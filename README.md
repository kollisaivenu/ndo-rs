# ndo-rs: Rust Library for performing Nested Dissection Ordering for a Sparse Matrix.

## This project is still a WIP. Use the following code snippet to try it out

```rust
use std::path::Path;
use ndo_rs::io::read_matrix_market_as_graph;
use ndo_rs::algorithms::nested_dissection_ordering::NestedDissectionOrdering;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_path = "./testdata/bcsstm12.mtx";
    let graph = read_matrix_market_as_graph(Path::new(file_path))?;
    let ordering = NestedDissectionOrdering {..Default::default()}.compute_ordering(&graph);
    println!("{:?}", ordering);
    
    Ok(())
}
```