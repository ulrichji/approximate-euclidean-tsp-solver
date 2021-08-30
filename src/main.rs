use num_traits::Float;

use kdtree::KdTree;
use kdtree::distance::squared_euclidean;

use std::path::PathBuf;
use structopt::StructOpt;

mod euclidean_distance_graph;
mod graph_drawer;
mod minimum_spanning_tree;
mod mst_graph_drawer;
mod random_point_generator;
mod tsp_mst_approximation;
mod tsp_point_path_drawer;

mod main_mst_image_maker;
mod main_tsp_image_maker;

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
    let point_list: Vec<Point<T, N>>
        = random_point_generator::RandomPointGenerator::<T, N>::
            new_square_range(points_range.0, points_range.1).take(number_of_points).collect();

    let mut kdtree = KdTree::new(N);
    for (i, p) in point_list.iter().enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    let mut euclidean_distance_graph: EuclideanPointDistanceGraph<T, N>
        = EuclideanPointDistanceGraph::new(point_list, kdtree);
    let mst: minimum_spanning_tree::MinimumSpanningTree<usize>
        = minimum_spanning_tree::MinimumSpanningTree::construct::<T>(
            &mut euclidean_distance_graph,
            minimum_spanning_tree::TreeSizeLimit::Finite(number_of_points));

    (mst, euclidean_distance_graph)
}

#[derive(Debug, StructOpt)]
#[structopt(name = "approximate_tsp_solver",
            about = "Solve the travelling salesman problem with the minimum spanning tree \
                     approximation, and plot the result as pretty looking images. It can also \
                     plot the minimum spanning tree that also can be quite decorative")]
enum ProgramOptions {
    Mst(RandomPointsImageMakerOptions),
    Tsp(RandomPointsImageMakerOptions)
}

#[derive(Debug, StructOpt)]
pub struct RandomPointsImageMakerOptions {
    #[structopt(short, long)]
    number_of_points: usize,

    #[structopt(parse(from_os_str))]
    result_image_path: PathBuf,

    #[structopt(long, default_value = "0.0")]
    points_min_range: f64,
    #[structopt(long, default_value = "1.0")]
    points_max_range: f64,

    #[structopt(long, default_value="512")]
    result_image_width: u32,
    #[structopt(long, default_value="512")]
    result_image_height: u32,

    #[structopt(long, default_value = "0.0")]
    result_image_bounds_min: f64,
    #[structopt(long, default_value = "1.0")]
    result_image_bounds_max: f64,
}

fn main() {
    let program_options = ProgramOptions::from_args();

    match program_options {
        ProgramOptions::Mst(mst_options) => {
            main_mst_image_maker::draw_mst(&mst_options);
        },
        ProgramOptions::Tsp(tsp_options) => {
            main_tsp_image_maker::draw_tsp(&tsp_options);
        }
    }
}
