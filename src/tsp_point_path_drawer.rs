use crate::tsp_mst_approximation;
use crate::graph_drawer;

pub fn drawable_vertex_iter<'a, T, const N: usize>(
        tsp_point_path: &'a tsp_mst_approximation::TspPointPath<T, N>)
            -> impl Iterator<Item=graph_drawer::DrawableVertex2D> + 'a
            where T: Into<f64> + Copy {
    tsp_point_path.iter_points()
                  .enumerate()
                  .map(|(i, n)| graph_drawer::DrawableVertex2D::new(n.dims[0].into(),
                                                                    n.dims[1].into(), i))
}

pub fn drawable_edge_iter<'a, T, const N: usize>(
        tsp_point_path: &'a tsp_mst_approximation::TspPointPath<T, N>)
            -> impl Iterator<Item=graph_drawer::DrawableEdge2D> + 'a
            where T: Into<f64> + Copy {
    let from_node_iter = tsp_point_path.iter_points();
    let to_node_iter = tsp_point_path.iter_points()
                                     .skip(1)
                                     .chain(tsp_point_path.iter_points().next());
    from_node_iter
        .zip(to_node_iter)
        .enumerate()
        .map(|(i, (from, to))| graph_drawer::DrawableEdge2D::new(
            graph_drawer::DrawableVertex2D::new(from.dims[0].into(), from.dims[1].into(), i),
            graph_drawer::DrawableVertex2D::new(to.dims[0].into(), to.dims[1].into(), i + 1)))
}
