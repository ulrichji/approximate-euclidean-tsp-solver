use kdtree::KdTree;

#[derive(Debug, Copy, Clone)]
pub struct Point<T, const N: usize> {
    pub dims: [T; N]
}

/*#[derive(PartialEq, PartialOrd)]
pub struct EuclideanEdgeWeight<T> {
    weight: T
}

impl<T> EuclideanEdgeWeight<T> {
    pub fn new(weight: T) -> EuclideanEdgeWeight<T> {
        EuclideanEdgeWeight { weight }
    }
}

impl<T> Eq for EuclideanEdgeWeight<T> where T: PartialOrd { }
impl<T> Ord for EuclideanEdgeWeight<T> where T: PartialOrd {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.weight.partial_cmp(&other.weight).unwrap()
    }
}*/

pub struct EuclideanPointDistanceGraph<T, const N: usize> {
    pub point_list: Vec<Point<T, N>>,
    pub kd_tree: KdTree<T, usize, [T; N]>,
    pub edge_counts: Vec<usize>,
    pub depths: Vec<Option<usize>>
}

impl<T, const N: usize> EuclideanPointDistanceGraph<T, N> {
    pub fn new(point_list: Vec<Point<T, N>>, kd_tree: KdTree<T, usize, [T; N]>)
            -> EuclideanPointDistanceGraph<T, N> {
        let edge_counts = vec![0; point_list.len()];
        let depths = vec![None; point_list.len()];

        EuclideanPointDistanceGraph {
            point_list,
            kd_tree,
            edge_counts,
            depths
        }
    }
}
