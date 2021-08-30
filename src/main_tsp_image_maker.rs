use crate::create_mst_from_random_points;
use crate::RandomPointsImageMakerOptions;
use crate::graph_drawer;
use crate::tsp_mst_approximation;
use crate::tsp_point_path_drawer;

pub fn draw_tsp(options: &RandomPointsImageMakerOptions) {
    type PointNumberType = f64;
    const DIMS: usize = 2;
    let points_range = (options.points_min_range, options.points_max_range);

    let (mst, euclidean_distance_graph) = create_mst_from_random_points::<PointNumberType, DIMS>(
        options.number_of_points, points_range);

    let tsp = tsp_mst_approximation::TspPointPath::from_traversed_mst_in_graph(
        &mst, &euclidean_distance_graph);

    let mut tsp_vertex_iter = tsp_point_path_drawer::drawable_vertex_iter(&tsp);
    let mut tsp_edge_iter = tsp_point_path_drawer::drawable_edge_iter(&tsp);
    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: options.result_image_width,
        image_height: options.result_image_height,
        x_min: options.result_image_bounds_min,
        x_max: options.result_image_bounds_max,
        y_min: options.result_image_bounds_min,
        y_max: options.result_image_bounds_max,
        max_depth: tsp.steps_in_path()
    };
    let image = graph_drawer::draw_2d_graph(
        &mut tsp_vertex_iter, &mut tsp_edge_iter, &drawing_properties);

    println!("Saving image to {}", &options.result_image_path
                                           .clone()
                                           .into_os_string()
                                           .into_string()
                                           .unwrap_or("".to_string()));
    image.save(&options.result_image_path).unwrap();
}
