use crate::EuclideanPointDistanceGraph;
use crate::Point;
use crate::minimum_spanning_tree;

pub struct TspPointPath<T, const N: usize> {
    points: Vec<Point<T, N>>,
}

impl<T, const N: usize> TspPointPath<T, N> where T: Copy {
    pub fn from_traversed_mst_in_graph(mst: &minimum_spanning_tree::MinimumSpanningTree<usize>,
                                       point_distance_graph: &EuclideanPointDistanceGraph<T, N>)
                                            -> TspPointPath<T, N> {
        let mut points = Vec::new();
        let mut backtrack_vertices = vec![0];

        while backtrack_vertices.len() > 0 {
            let current_index = backtrack_vertices.pop().unwrap();
            let current_vertex = &mst.elements[current_index];

            points.push(point_distance_graph.point_list[current_vertex.vertex_identifier]);

            // The backtrack_vertices is used as a stack, so adding the sibling first will ensure
            // that the graph is traversed preorder (children before sibling)
            match current_vertex.next_sibling_index {
                Some(i) => {
                    backtrack_vertices.push(i);
                },
                None => {}
            };

            // Ref. comment above, the child is added last and will be the next to be evaluated,
            // ensuring preorder traversal.
            match current_vertex.first_child_index {
                Some(i) => {
                    backtrack_vertices.push(i)
                },
                None => {}
            };
        }

        TspPointPath {
            points
        }
    }

    pub fn iter_points(&self) -> PointIterator<T, N> {
        PointIterator {
            points_ref: &self.points,
            current_point_index: 0
        }
    }

    pub fn steps_in_path(&self) -> usize { self.points.len() }
}

pub struct PointIterator<'a, T, const N: usize> {
    points_ref: &'a Vec<Point<T, N>>,
    current_point_index: usize
}

impl<'a, T, const N: usize> Iterator for PointIterator<'a, T, N> {
    type Item = &'a Point<T, N>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_point_index < self.points_ref.len() {
            let vertex = &self.points_ref[self.current_point_index];
            self.current_point_index += 1;

            Some(vertex)
        } else {
            None
        }
    }
}
