extern crate num_traits;

use rand::distributions::{Distribution, Uniform, uniform::SampleUniform};

use kdtree::KdTree;
use kdtree::kdtree::NearestIter;
use kdtree::distance::squared_euclidean;

mod minimum_spanning_tree;
mod graph_drawer;

trait Zero {
    fn zero() -> Self;
}

#[derive(Debug, Copy, Clone)]
struct Point<T, const N: usize> {
    dims: [T; N]
}

struct RandomPointGenerator<T, const N: usize> where T: SampleUniform + Copy {
    random_number_generator: rand::rngs::ThreadRng,
    ranges: [(T, T); N]
}

impl<T, const N: usize> RandomPointGenerator<T, N> where T: SampleUniform + Copy {
    pub fn new(ranges: [(T, T); N]) -> RandomPointGenerator<T, N> {
        RandomPointGenerator {
            random_number_generator: rand::thread_rng(),
            ranges
        }
    }

    pub fn new_square_range(min_range: T, max_range: T) -> RandomPointGenerator<T, N> {
        RandomPointGenerator {
            random_number_generator: rand::thread_rng(),
            ranges: [(min_range, max_range); N]
        }
    }
}

impl<T, const N: usize> Iterator for RandomPointGenerator<T, N> where T: SampleUniform + Copy + Zero {
    type Item = Point<T, N>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut dims = [T::zero(); N];
        for (dim, range) in dims.iter_mut().zip(&self.ranges) {
            let distribution = Uniform::new(range.0, range.1);
            *dim = distribution.sample(&mut self.random_number_generator);
        }
        Some(Point{ dims })
    }
}

struct PointWithId<T, const N: usize> {
    point: Point<T, N>,
    id: usize
}

#[derive(PartialEq, PartialOrd)]
struct EuclideanEdgeWeight {
    weight: f64
}

impl EuclideanEdgeWeight {
    fn new(weight: f64) -> EuclideanEdgeWeight {
        EuclideanEdgeWeight { weight }
    }
}

impl Eq for EuclideanEdgeWeight { }
impl Ord for EuclideanEdgeWeight {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if other.weight < other.weight {
            std::cmp::Ordering::Less
        } else if self.weight > other.weight {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

/*
type DistanceCallbackType<'a> = dyn 'a + Fn(&[f64], &[f64]) -> f64;
struct EdgeIterator<'a, 'b> {
    point: Point<f64, 3>,
    nearest_iter: NearestIter<'a, 'b, f64, usize, [f64; 3], DistanceCallbackType<'static>>
}
*/

struct EuclideanPointDistanceGraph<const N: usize> {
    point_list: Vec<Point<f64, N>>,
    kd_tree: KdTree<f64, usize, [f64; N]>,
    edge_counts: Vec<usize>,
    depths: Vec<Option<usize>>
}

impl<const N: usize> EuclideanPointDistanceGraph<N> {
    fn new(point_list: Vec<Point<f64, N>>, kd_tree: KdTree<f64, usize, [f64; N]>) -> EuclideanPointDistanceGraph<N> {
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

impl<const N: usize> minimum_spanning_tree::Graph<usize, EuclideanEdgeWeight> for EuclideanPointDistanceGraph<N> {
    fn get_root_vertex_identifier(&self) -> usize { 0usize }
    fn next_neighbour_for_vertex(&mut self, vertex_identifier: &usize)
            -> std::option::Option<minimum_spanning_tree::Edge<usize, EuclideanEdgeWeight>> {

        let point = &self.point_list[*vertex_identifier];
        let edge_count = &mut self.edge_counts[*vertex_identifier];
        *edge_count += 1;

        let maybe_next_neighbour = self.kd_tree.iter_nearest(&point.dims, &squared_euclidean).unwrap().skip(*edge_count).next();

        match maybe_next_neighbour {
            Some(neighbour) => {
                Some(minimum_spanning_tree::Edge::new(*vertex_identifier,
                                                      *neighbour.1,
                                                      EuclideanEdgeWeight::new(neighbour.0)))
            },
            None => None
        }
    }
}

fn set_graph_depths_from_mst<const N: usize>(graph: &mut EuclideanPointDistanceGraph<N>,
                                             mst: &minimum_spanning_tree::MinimumSpanningTree<usize>,
                                             vertex_index: usize,
                                             depth: usize) {
    let vertex_at_index = &mst.elements[vertex_index];
    let vertex_index_in_point_list = vertex_at_index.vertex_identifier;
    if graph.depths[vertex_index_in_point_list].is_none() {
        graph.depths[vertex_index_in_point_list] = Some(depth);

        if vertex_at_index.next_sibling_index.is_some() {
            set_graph_depths_from_mst(graph, mst, vertex_at_index.next_sibling_index.unwrap(), depth);
        }

        if vertex_at_index.first_child_index.is_some() {
            set_graph_depths_from_mst(graph, mst, vertex_at_index.first_child_index.unwrap(), depth + 1);
        }
    }



}

struct MstPointVertexIterator<'a, 'b, const N: usize> {
    mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>,
    point_distance_graph: &'b EuclideanPointDistanceGraph<N>,
    current_vertex_index: usize
}

impl <'a, 'b, const N: usize> MstPointVertexIterator<'a, 'b, N> {
    fn new(mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>, point_distance_graph: &'b EuclideanPointDistanceGraph<N>) -> MstPointVertexIterator<'a, 'b, N>  {
        MstPointVertexIterator {
            mst,
            point_distance_graph,
            current_vertex_index: 0
        }
    }
}

impl <'a, 'b, const N: usize> Iterator for MstPointVertexIterator<'a, 'b, N> {
    type Item = graph_drawer::DrawableVertex2D;
    fn next(&mut self) -> std::option::Option<<Self as std::iter::Iterator>::Item> {
        if self.current_vertex_index < self.mst.elements.len() {
            let vertex = &self.point_distance_graph.point_list[self.mst.elements[self.current_vertex_index].vertex_identifier];
            let vertex_depth = self.point_distance_graph.depths[self.mst.elements[self.current_vertex_index].vertex_identifier].unwrap_or(0);
            self.current_vertex_index += 1;
            Some(graph_drawer::DrawableVertex2D::new(vertex.dims[0], vertex.dims[1], vertex_depth))
        } else {
            None
        }
    }
}

struct MstPointEdgeIterator<'a, 'b, const N: usize> {
    mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>,
    point_distance_graph: &'b EuclideanPointDistanceGraph<N>,
    current_from_vertex_index: usize,
    current_to_vertex_index: usize
}

impl <'a, 'b, const N: usize> MstPointEdgeIterator<'a, 'b, N> {
    fn new(mst: &'a minimum_spanning_tree::MinimumSpanningTree<usize>, point_distance_graph: &'b EuclideanPointDistanceGraph<N>) -> MstPointEdgeIterator<'a, 'b, N>  {
        MstPointEdgeIterator {
            mst,
            point_distance_graph,
            current_from_vertex_index: 0,
            current_to_vertex_index: mst.elements[0].first_child_index.unwrap_or(0)
        }
    }

    fn find_next_edge(&mut self) {
        let to_vertex = &self.mst.elements[self.current_to_vertex_index];

        self.current_from_vertex_index = match to_vertex.next_sibling_index {
            Some(_) => self.current_from_vertex_index,
            None => {
                let mut next_from_vertex_index = self.current_from_vertex_index + 1;
                while next_from_vertex_index < self.mst.elements.len() && self.mst.elements[next_from_vertex_index].first_child_index.is_none() {
                    next_from_vertex_index += 1;
                };
                next_from_vertex_index
            }
        };

        self.current_to_vertex_index = match to_vertex.next_sibling_index {
            Some(i) => i,
            None => {
                if self.current_from_vertex_index < self.mst.elements.len() {
                    self.mst.elements[self.current_from_vertex_index].first_child_index.unwrap()
                } else {
                    0
                }
            }
        };
    }
}

impl <'a, 'b, const N: usize> Iterator for MstPointEdgeIterator<'a, 'b, N> {
    type Item = graph_drawer::DrawableEdge2D;
    fn next(&mut self) -> std::option::Option<<Self as std::iter::Iterator>::Item> {
        if self.current_from_vertex_index < self.mst.elements.len() {
            let from_vertex = &self.point_distance_graph.point_list[self.mst.elements[self.current_from_vertex_index].vertex_identifier];
            let to_vertex = &self.point_distance_graph.point_list[self.mst.elements[self.current_to_vertex_index].vertex_identifier];
            let from_vertex_depth = self.point_distance_graph.depths[self.mst.elements[self.current_from_vertex_index].vertex_identifier].unwrap_or(0);
            let to_vertex_depth = self.point_distance_graph.depths[self.mst.elements[self.current_to_vertex_index].vertex_identifier].unwrap_or(0);

            self.find_next_edge();

            Some(graph_drawer::DrawableEdge2D::new(graph_drawer::DrawableVertex2D::new(from_vertex.dims[0], from_vertex.dims[1], from_vertex_depth),
                                                   graph_drawer::DrawableVertex2D::new(to_vertex.dims[0], to_vertex.dims[1], to_vertex_depth)))
        } else {
            None
        }
    }
}

struct TspPointPath<T, const N: usize> {
    points: Vec<Point<T, N>>,
}

struct TspPointPathVertexIterator<'a, const N: usize> {
    tsp_point_path: &'a TspPointPath<f64, N>,
    current_index: usize
}

impl <'a, const N: usize> TspPointPathVertexIterator<'a, N> {
    fn new(tsp_point_path: &'a TspPointPath<f64, N>) -> TspPointPathVertexIterator<'a, N> {
        TspPointPathVertexIterator {
            tsp_point_path,
            current_index: 0
        }
    }
}

impl<'a, const N: usize> Iterator for TspPointPathVertexIterator<'a, N> {
    type Item = graph_drawer::DrawableVertex2D;
    fn next(&mut self) -> Option<graph_drawer::DrawableVertex2D> {
        if self.current_index < self.tsp_point_path.points.len() {
            let vertex = self.tsp_point_path.points[self.current_index];
            self.current_index += 1;

            Some(graph_drawer::DrawableVertex2D::new(vertex.dims[0], vertex.dims[1], self.current_index - 1))
        } else {
            None
        }
    }
}

struct TspPointPathEdgeIterator<'a, const N: usize> {
    tsp_point_path: &'a TspPointPath<f64, N>,
    current_index: usize
}

impl <'a, const N: usize> TspPointPathEdgeIterator<'a, N> {
    fn new(tsp_point_path: &'a TspPointPath<f64, N>) -> TspPointPathEdgeIterator<'a, N> {
        TspPointPathEdgeIterator {
            tsp_point_path,
            current_index: 0
        }
    }
}

impl<'a, const N: usize> Iterator for TspPointPathEdgeIterator<'a, N> {
    type Item = graph_drawer::DrawableEdge2D;
    fn next(&mut self) -> Option<graph_drawer::DrawableEdge2D> {
        if self.current_index < self.tsp_point_path.points.len() {
            let from_vertex = self.tsp_point_path.points[self.current_index];
            self.current_index += 1;
            let to_vertex = self.tsp_point_path.points[self.current_index % self.tsp_point_path.points.len()];

            Some(graph_drawer::DrawableEdge2D::new(
                graph_drawer::DrawableVertex2D::new(from_vertex.dims[0], from_vertex.dims[1], self.current_index - 1),
                graph_drawer::DrawableVertex2D::new(to_vertex.dims[0], to_vertex.dims[1], self.current_index)))
        } else {
            None
        }
    }
}

impl<const N: usize> TspPointPath<f64, N> {
    fn from_traversed_mst_in_graph(mst: &minimum_spanning_tree::MinimumSpanningTree<usize>,
                                   point_distance_graph: &EuclideanPointDistanceGraph<N>)
                                        -> TspPointPath<f64, N> {
        let mut points = Vec::new();

        let mut backtrack_vertices: Vec<usize> = Vec::new();
        backtrack_vertices.push(0);

        while backtrack_vertices.len() > 0 {
            let current_index = backtrack_vertices.pop().unwrap();
            let current_vertex = &mst.elements[current_index];

            points.push(point_distance_graph.point_list[current_vertex.vertex_identifier]);

            match current_vertex.next_sibling_index {
                Some(i) => {
                    backtrack_vertices.push(i);
                },
                None => {}
            };

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
}

impl Zero for f64 {
    fn zero() -> f64 { 0f64 }
}

use std::time::{Instant};
fn main() {
    //benchmark_kd_tree();
    //benchmark_mst_construction();
    //draw_mst();
    draw_tsp();
}

fn draw_tsp() {
    const dims: usize = 2;

    let number_of_points = 1000;

    println!("Generating points");
    let point_list: Vec<Point<f64, dims>> = RandomPointGenerator::<f64, dims>::new_square_range(0.0, 2.0).take(number_of_points).collect();

    println!("Building tree");
    let mut kdtree = KdTree::new(dims);
    for (i, p) in point_list.iter().enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    println!("Building mst");
    let mut euclidean_distance_graph: EuclideanPointDistanceGraph<dims> = EuclideanPointDistanceGraph::new(point_list, kdtree);
    let mst: minimum_spanning_tree::MinimumSpanningTree<usize> = minimum_spanning_tree::MinimumSpanningTree::construct::<EuclideanEdgeWeight>(&mut euclidean_distance_graph, minimum_spanning_tree::TreeSize::Finite(number_of_points));

    println!("Computing TSP path from mst");
    let tsp = TspPointPath::from_traversed_mst_in_graph(&mst, &euclidean_distance_graph);

    println!("Drawing TSP trajectory to image");
    let mut tsp_vertex_iter = TspPointPathVertexIterator::new(&tsp);
    let mut tsp_edge_iter = TspPointPathEdgeIterator::new(&tsp);
    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: 512,
        image_height: 512,
        x_min: 0.0,
        x_max: 2.0,
        y_min: 0.0,
        y_max: 2.0,
        max_depth: tsp.points.len()
    };
    let image = graph_drawer::draw_2d_graph(
        &mut tsp_vertex_iter, &mut tsp_edge_iter, &drawing_properties);

    println!("Saving image to tsp.png");
    image.save("tsp.png").unwrap();
}

fn draw_mst() {
    const dims: usize = 2;

    let number_of_points = 6000000;

    println!("Generating points");
    let point_list: Vec<Point<f64, dims>> = RandomPointGenerator::<f64, dims>::new_square_range(0.0, 2.0).take(number_of_points).collect();

    println!("Building tree");
    let mut kdtree = KdTree::new(dims);
    for (i, p) in point_list.iter().enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    println!("Building mst");
    let mut euclidean_distance_graph: EuclideanPointDistanceGraph<dims> = EuclideanPointDistanceGraph::new(point_list, kdtree);
    let mst: minimum_spanning_tree::MinimumSpanningTree<usize> = minimum_spanning_tree::MinimumSpanningTree::construct::<EuclideanEdgeWeight>(&mut euclidean_distance_graph, minimum_spanning_tree::TreeSize::Finite(number_of_points));

    println!("Computing tree depths");
    set_graph_depths_from_mst(&mut euclidean_distance_graph, &mst, 0, 0);

    println!("Drawing mst to image");
    let mut mst_vertex_iter = MstPointVertexIterator::new(&mst, &euclidean_distance_graph);
    let mut mst_edge_iter = MstPointEdgeIterator::new(&mst, &euclidean_distance_graph);
    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: 5900,
        image_height: 5900,
        x_min: 0.0,
        x_max: 2.0,
        y_min: 0.0,
        y_max: 2.0,
        max_depth: euclidean_distance_graph.depths.iter().map(|v| v.unwrap_or(0)).max().unwrap()
    };
    let image = graph_drawer::draw_2d_graph(
        &mut mst_vertex_iter, &mut mst_edge_iter, &drawing_properties);

    println!("Saving image");
    image.save("mst.png").unwrap();
}

fn benchmark_mst_construction() {
    const dims: usize = 3;

    let number_of_points = 100000;

    let pre_build_time_instant = Instant::now();

    println!("Generating points");
    let point_list: Vec<Point<f64, dims>> = RandomPointGenerator::<f64, dims>::new_square_range(0.0, 2.0).take(number_of_points).collect();

    println!("Building tree");
    let mut kdtree = KdTree::new(dims);
    for (i, p) in point_list.iter().enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    let build_duration = pre_build_time_instant.elapsed();
    let pre_mst_build_time_instant = Instant::now();

    println!("Building mst");

    let mut euclidean_distance_graph: EuclideanPointDistanceGraph<dims> = EuclideanPointDistanceGraph::new(point_list, kdtree);
    let mst: minimum_spanning_tree::MinimumSpanningTree<usize> = minimum_spanning_tree::MinimumSpanningTree::construct::<EuclideanEdgeWeight>(&mut euclidean_distance_graph, minimum_spanning_tree::TreeSize::Finite(number_of_points));
    //println!("{:?}", mst);

    let mst_build_duration = pre_mst_build_time_instant.elapsed();
    println!("Done. Building kd-tree with {} points took {} seconds and {} millis, building mst took {} seconds and {} millis",
        number_of_points, build_duration.as_secs(), build_duration.subsec_millis(),
        mst_build_duration.as_secs(), mst_build_duration.subsec_millis());
}

fn benchmark_kd_tree() {
    const dims: usize = 3;

    let number_of_points = 20000000;
    let number_of_searches = 500000;

    let pre_build_time_instant = Instant::now();

    println!("Building tree");
    let mut kdtree = KdTree::new(dims);
    for (i, p) in RandomPointGenerator::<f64, dims>::new_square_range(0.0, 2.0).take(number_of_points).enumerate() {
        kdtree.add(p.dims, i).unwrap();
    }

    let build_duration = pre_build_time_instant.elapsed();
    let pre_search_instant = Instant::now();

    println!("Searching tree");
    let nearest_points: Vec<(f64, &usize)> = RandomPointGenerator::<f64, dims>::new_square_range(0.0, 2.0)
        .take(number_of_searches)
        .map(|p| kdtree.nearest(&p.dims, 1, &squared_euclidean).unwrap()[0])
        .collect();

    let search_duration = pre_search_instant.elapsed();

    for (o, _) in nearest_points.iter().take(1000) {
        println!("Nearest: {:?}", o);
    }

    println!("Done. Building {} points took {} seconds and {} millis, searching {} points took {} seconds and {} millis, which is {} millis per search",
        number_of_points, build_duration.as_secs(), build_duration.subsec_millis(),
        number_of_searches, search_duration.as_secs(), search_duration.subsec_millis(), search_duration.as_millis() as f64 / number_of_searches as f64);
}
