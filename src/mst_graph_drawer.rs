
use crate::Point;
use crate::minimum_spanning_tree;
use crate::graph_drawer;

pub fn drawable_vertex_iter<'a, T, const N: usize>(
        mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>,
        point_list: &'a Vec<Point<T, N>>)
            -> impl Iterator<Item=graph_drawer::DrawableVertex2D> + 'a
            where T: Into<f64> + Copy {
    mst.traverse_nodes_preorder()
       .map(move |(n, l)| {
           graph_drawer::DrawableVertex2D::new(point_list[n.vertex_identifier].dims[0].into(),
                                               point_list[n.vertex_identifier].dims[1].into(),
                                               l)
       })
}

pub fn drawable_edge_iter<'a, T, const N: usize>(
        mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>,
        point_list: &'a Vec<Point<T, N>>)
            -> impl Iterator<Item=graph_drawer::DrawableEdge2D> + 'a
            where T: Into<f64> + Copy {
    mst.traverse_edges_preorder()
       .map(move |(f, t, l)| {
           graph_drawer::DrawableEdge2D::new(
               graph_drawer::DrawableVertex2D::new(point_list[f.vertex_identifier].dims[0].into(),
                                                   point_list[f.vertex_identifier].dims[1].into(),
                                                   l),
               graph_drawer::DrawableVertex2D::new(point_list[t.vertex_identifier].dims[0].into(),
                                                   point_list[t.vertex_identifier].dims[1].into(),
                                                   l + 1))
       })
}
