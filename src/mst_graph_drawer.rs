use crate::Point;
use crate::MstPointEdgeIterator;
use crate::MstPointVertexIterator;

use crate::graph_drawer;

impl <'a, 'b, T, const N: usize> Iterator for MstPointVertexIterator<'a, 'b, T, N>
        where T: Into<f64> + Copy {

    type Item = graph_drawer::DrawableVertex2D;

    fn next(&mut self) -> std::option::Option<<Self as std::iter::Iterator>::Item> {
        if self.current_vertex_index < self.mst.elements.len() {

            let vertex_identifier = self.mst.vertex_identifier_at(self.current_vertex_index);
            let vertex = self.point_distance_graph.get_point_at(*vertex_identifier);
            let vertex_depth = self.point_distance_graph.graph_depth_at(*vertex_identifier)
                                                        .unwrap_or(0);

            self.current_vertex_index += 1;
            Some(graph_drawer::DrawableVertex2D::new(vertex.dims[0].into(),
                                                     vertex.dims[1].into(),
                                                     vertex_depth))
        } else {
            None
        }
    }
}

impl<'a, 'b, T, const N: usize> MstPointEdgeIterator<'a, 'b, T, N> {
    fn current_from_vertex(&self) -> &Point<T, N>{
        let vertex_identifier = self.mst.vertex_identifier_at(self.current_from_vertex_index);
        self.point_distance_graph.get_point_at(*vertex_identifier)
    }

    fn current_to_vertex(&self) -> &Point<T, N> {
        let vertex_identifier = self.mst.vertex_identifier_at(self.current_to_vertex_index);
        self.point_distance_graph.get_point_at(*vertex_identifier)
    }

    fn current_from_depth(&self) -> Option<usize> {
        let vertex_identifier = self.mst.vertex_identifier_at(self.current_from_vertex_index);
        self.point_distance_graph.graph_depth_at(*vertex_identifier)
    }

    fn current_to_depth(&self) -> Option<usize> {
        let vertex_identifier = self.mst.vertex_identifier_at(self.current_to_vertex_index);
        self.point_distance_graph.graph_depth_at(*vertex_identifier)
    }
}

impl <'a, 'b, T, const N: usize> Iterator for MstPointEdgeIterator<'a, 'b, T, N>
        where T: Into<f64> + Copy {

    type Item = graph_drawer::DrawableEdge2D;

    fn next(&mut self) -> std::option::Option<<Self as std::iter::Iterator>::Item> {
        if self.current_from_vertex_index < self.mst.elements.len() {
            let from_vertex = self.current_from_vertex();
            let to_vertex = self.current_to_vertex();

            let from_vertex_depth = self.current_from_depth().unwrap_or(0);
            let to_vertex_depth = self.current_to_depth().unwrap_or(0);

            let result_edge = graph_drawer::DrawableEdge2D::new(
                graph_drawer::DrawableVertex2D::new(from_vertex.dims[0].into(),
                                                    from_vertex.dims[1].into(),
                                                    from_vertex_depth),
                graph_drawer::DrawableVertex2D::new(to_vertex.dims[0].into(),
                                                    to_vertex.dims[1].into(),
                                                    to_vertex_depth));

            self.find_next_edge();

            Some(result_edge)

        } else {
            None
        }
    }
}
