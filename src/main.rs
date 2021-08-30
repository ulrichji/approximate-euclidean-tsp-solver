use num_traits::Float;

use kdtree::KdTree;
use kdtree::distance::squared_euclidean;

mod euclidean_distance_graph;
mod graph_drawer;
mod minimum_spanning_tree;
mod mst_graph_drawer;
mod random_point_generator;
mod tsp_mst_approximation;
mod tsp_point_path_drawer;

use euclidean_distance_graph::{EuclideanPointDistanceGraph, Point};

impl<T, const N: usize> minimum_spanning_tree::Graph<usize, T>
        for EuclideanPointDistanceGraph<T, N>
        where T: Float + std::cmp::PartialOrd {
    fn get_root_vertex_identifier(&self) -> usize { 0usize }
    fn next_neighbour_for_vertex(&mut self, vertex_identifier: &usize)
            -> std::option::Option<minimum_spanning_tree::Edge<usize, T>> {

        let point = &self.point_list[*vertex_identifier];
        let edge_count = &mut self.edge_counts[*vertex_identifier];
        *edge_count += 1;

        let maybe_next_neighbour = self.kd_tree.iter_nearest(&point.dims, &squared_euclidean)
                                               .unwrap()
                                               .skip(*edge_count)
                                               .next();

        match maybe_next_neighbour {
            Some(neighbour) => {
                Some(minimum_spanning_tree::Edge::new(*vertex_identifier,
                                                      *neighbour.1,
                                                      neighbour.0))
            },
            None => None
        }
    }
}


fn create_mst_from_random_points<T, const N: usize>(number_of_points: usize, points_range: (T, T))
        -> (minimum_spanning_tree::MinimumSpanningTree<usize>, EuclideanPointDistanceGraph<T, N>)
        where T: rand::distributions::uniform::SampleUniform +
                 Copy +
                 random_point_generator::Zero +
                 num_traits::Float {
    println!("Generating points");
    let point_list: Vec<Point<T, N>>
        = random_point_generator::RandomPointGenerator::<T, N>::
            new_square_range(points_range.0, points_range.1).take(number_of_points).collect();

    println!("Building tree");
    let mut kdtree = KdTree::new(N);
    for (i, p) in point_list.iter().enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    println!("Building mst");
    let mut euclidean_distance_graph: EuclideanPointDistanceGraph<T, N>
        = EuclideanPointDistanceGraph::new(point_list, kdtree);
    let mst: minimum_spanning_tree::MinimumSpanningTree<usize>
        = minimum_spanning_tree::MinimumSpanningTree::construct::<T>(
            &mut euclidean_distance_graph,
            minimum_spanning_tree::TreeSizeLimit::Finite(number_of_points));

    (mst, euclidean_distance_graph)
}

fn main() {
    draw_mst();
    draw_tsp();
}

fn draw_tsp() {
    type PointNumberType = f64;
    const DIMS: usize = 2;
    let number_of_points = 10000;
    let points_range = (0.0, 1.0);

    let (mst, euclidean_distance_graph) = create_mst_from_random_points::<PointNumberType, DIMS>(
        number_of_points, points_range);

    println!("Computing TSP path from mst");
    let tsp = tsp_mst_approximation::TspPointPath::from_traversed_mst_in_graph(
        &mst, &euclidean_distance_graph);

    println!("Drawing TSP trajectory to image");
    let mut tsp_vertex_iter = tsp_point_path_drawer::drawable_vertex_iter(&tsp);
    let mut tsp_edge_iter = tsp_point_path_drawer::drawable_edge_iter(&tsp);
    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: 512,
        image_height: 512,
        x_min: points_range.0,
        x_max: points_range.1,
        y_min: points_range.0,
        y_max: points_range.1,
        max_depth: tsp.steps_in_path()
    };
    let image = graph_drawer::draw_2d_graph(
        &mut tsp_vertex_iter, &mut tsp_edge_iter, &drawing_properties);

    println!("Saving image to tsp.png");
    image.save("tsp.png").unwrap();
}

fn draw_mst() {
    type PointNumberType = f64;
    const DIMS: usize = 2;
    let number_of_points = 10000;
    let points_range = (0.0, 1.0);

    let (mst, euclidean_distance_graph) = create_mst_from_random_points::<PointNumberType, DIMS>(
        number_of_points, points_range);

    println!("Drawing mst to image");
    let mut mst_vertex_iter = mst_graph_drawer::drawable_vertex_iter(
        &mst, &euclidean_distance_graph.point_list);
    let mut mst_edge_iter = mst_graph_drawer::drawable_edge_iter(
        &mst, &euclidean_distance_graph.point_list);
    let max_depth = mst.traverse_nodes_preorder().map(|(_, l)| l).max().unwrap_or(0);

    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: 512,
        image_height: 512,
        x_min: points_range.0,
        x_max: points_range.1,
        y_min: points_range.0,
        y_max: points_range.1,
        max_depth
    };
    let image = graph_drawer::draw_2d_graph(
        &mut mst_vertex_iter, &mut mst_edge_iter, &drawing_properties);

    println!("Saving image");
    image.save("mst.png").unwrap();
}
