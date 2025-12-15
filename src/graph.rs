use std::iter::{Cloned, Zip};
use std::slice::Iter;
use ::sprs::CsMat;

/// Struct that represents a graph
pub struct Graph{
    /// The CsMat (from sprs) is used to store the graph as a sparse matrix in CSR format
    pub graph_csr: CsMat<f64>
}

impl Graph {

    /// Create a new graph
    pub fn new() -> Self {
        Self {
            graph_csr: CsMat::empty(sprs::CSR, 0)
        }
    }

    /// The number of vertices in the graph.
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.graph_csr.rows(), self.graph_csr.cols());
        self.graph_csr.rows()
    }

    /// An iterator over the neighbors of the given vertex.
    pub fn neighbors(&self, vertex: usize) -> Zip<Cloned<Iter<'_, usize>>, Cloned<Iter<'_, f64>>> {
        let (indices, data) = self.graph_csr.outer_view(vertex).unwrap().into_raw_storage();
        indices.iter().cloned().zip(data.iter().cloned())
    }

    /// Insert an edge with two vertices on either ends.
    pub fn insert(&mut self, vertex1: usize, vertex2: usize, edge_weight: f64) {
        self.graph_csr.insert(vertex1, vertex2, edge_weight);
    }

    /// Get edge weight for a pair of vertices.
    pub fn get_edge_weight(&self, vertex1: usize, vertex2: usize) -> Option<f64> {
        self.graph_csr.get(vertex1, vertex2).cloned()
    }

    /// Clone the graph
    pub fn clone(&self) -> Self {
        Self {
            graph_csr: self.graph_csr.clone()
        }
    }
}