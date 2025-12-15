use std::path::Path;
use sprs::io::{read_matrix_market, IoError};
use crate::graph::Graph;

/// Read a matrix market file and output Graph struct.
pub fn read_matrix_market_as_graph(file_path: &Path) -> Result<Graph, IoError> {
    // read the matrix market file as a TriMat with edge lengths.
    let tri_mat = read_matrix_market(file_path);

    match tri_mat {
        Ok(tri_matrix) => {
            // Read was successful, we return it after converting to CSR.
            let csr_matrix = tri_matrix.to_csr();
            Ok(Graph {
                graph_csr: csr_matrix,
            })
        },
        Err(e) => {
            Err(e)
        }
    }
}