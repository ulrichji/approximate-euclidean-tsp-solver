use kdtree::KdTree;

#[derive(Debug, Copy, Clone)]
pub struct Point<T, const N: usize> {
    pub dims: [T; N]
}

pub struct EuclideanPointDistanceGraph<T, const N: usize> {
    pub point_list: Vec<Point<T, N>>,
    pub kd_tree: KdTree<T, usize, [T; N]>,
    pub edge_counts: Vec<usize>
}

impl<T, const N: usize> EuclideanPointDistanceGraph<T, N> {
    pub fn new(point_list: Vec<Point<T, N>>, kd_tree: KdTree<T, usize, [T; N]>)
            -> EuclideanPointDistanceGraph<T, N> {
        let edge_counts = vec![0; point_list.len()];

        EuclideanPointDistanceGraph {
            point_list,
            kd_tree,
            edge_counts
        }
    }
}
