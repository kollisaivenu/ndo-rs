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

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Write;
    use std::path::Path;
    use crate::io::read_matrix_market_as_graph;
    use tempfile::tempdir;

    fn create_mock_file(dir: &Path, filename: &str, content: &str) -> String {
        let file_path = dir.join(filename);
        let mut file = File::create(&file_path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file_path.to_str().unwrap().to_string()
    }

    #[test]
    fn test_read_matrix_market() -> Result<(), std::io::Error> {
        let temp_dir = tempdir()?;
        let integer_content = "%%MatrixMarket matrix coordinate real symmetric\n%\n5 5 3\n1 1 1.0\n2 2 2.0\n5 5 5.0\n";
        let integer_matrix_file_path = create_mock_file(temp_dir.path(), "integer_matrix.mtx", integer_content);

        let graph = read_matrix_market_as_graph(&Path::new(&integer_matrix_file_path)).unwrap();

        assert_eq!(graph.graph_csr.rows(), 5);
        assert_eq!(graph.graph_csr.cols(), 5);
        assert_eq!(graph.graph_csr.nnz(), 3);

        Ok(())
    }
}