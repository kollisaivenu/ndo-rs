use std::path::Path;
use faer::perm::Perm;
use faer::{Side};
use faer::sparse::{SymbolicSparseColMat, SymbolicSparseColMatRef};
use faer::sparse::linalg::cholesky::{factorize_symbolic_cholesky, SymmetricOrdering};
use sprs::TriMat;
use ndo_rs::algorithms::nested_dissection_ordering::{ NestedDissectionOrdering};
use ndo_rs::io::read_matrix_market_as_graph;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_paths = vec!["bcsstk15.mtx", "bcsstk16.mtx", "bcsstk17.mtx", "copter2.mtx", "cfd1.mtx", "cfd2.mtx"];

    for path in file_paths {
        println!("Loading test file: {}", path);
        let col_mat_result = create_symbolic_col_matrix(Path::new(&("testdata/".to_owned() + path)));
        match col_mat_result {
            Ok(col_mat)  => {
                println!("Non zero values with no permutation :{}", calculate_num_of_nnz_without_permutation(col_mat.as_ref()));
                println!("Non zero values with AMD :{}", calculate_num_of_nnz_with_amd(col_mat.as_ref()));
                println!("Non zero values with Nested Dissection Ordering :{}", calculate_num_of_nnz_with_ndo(Path::new(&("testdata/".to_owned() + path)), col_mat.as_ref()).unwrap());
            }
            Err(e) => {
                println!("Something went wrong with opening the test file {}", e);
            }
        }
        println!("======================================");
    }

    Ok(())
}
fn create_symbolic_col_matrix(file_path: &Path) -> Result<SymbolicSparseColMat<usize>, Box<dyn std::error::Error>> {
    let triplets: TriMat<f64> = sprs::io::read_matrix_market(Path::new(file_path))?;

    // Change graph to remove self loops
    let mut new_triplets = TriMat::new((triplets.rows(), triplets.cols()));

    for (_, (vertex1, vertex2)) in triplets.triplet_iter() {
        if vertex1 != vertex2 {
            new_triplets.add_triplet(vertex1, vertex2, 1);
        }
    }

    let sprs_csc = new_triplets.to_csc();
    let rows = sprs_csc.rows();
    let cols = sprs_csc.cols();
    let (ind_ptr, indices, _) = sprs_csc.into_raw_storage();

    let faer_symbolic_col_mat = SymbolicSparseColMat::new_checked(
        rows,
        cols,
        ind_ptr,
        None,
        indices,
    );
    Ok(faer_symbolic_col_mat)
}

fn calculate_num_of_nnz_with_amd(symbolic_col_mat_ref: SymbolicSparseColMatRef<usize>) -> usize {
    factorize_symbolic_cholesky(symbolic_col_mat_ref, Side::Lower, SymmetricOrdering::Amd, Default::default()).unwrap().len_val()
}

fn calculate_num_of_nnz_without_permutation(symbolic_col_mat_ref: SymbolicSparseColMatRef<usize>) -> usize {
    factorize_symbolic_cholesky(symbolic_col_mat_ref, Side::Lower, SymmetricOrdering::Identity, Default::default()).unwrap().len_val()
}

fn calculate_num_of_nnz_with_ndo(file_path: &Path, symbolic_col_mat_ref: SymbolicSparseColMatRef<usize>) -> Result<(usize), Box<dyn std::error::Error>> {
    let graph = read_matrix_market_as_graph(Path::new(file_path))?;
    let (ordering, inverse_ordering) = NestedDissectionOrdering {..Default::default()}.compute_ordering(&graph);
    let dim = ordering.len();

    let perm = Perm::new_checked(
        ordering.into_boxed_slice(),
        inverse_ordering.into_boxed_slice(),
        dim,
    );

    let ndo_ordering = SymmetricOrdering::Custom(perm.as_ref());
    Ok(factorize_symbolic_cholesky(symbolic_col_mat_ref, Side::Lower, ndo_ordering, Default::default()).unwrap().len_val())

}