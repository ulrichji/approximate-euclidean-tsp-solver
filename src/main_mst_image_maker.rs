use crate::create_mst_from_random_points;
use crate::RandomPointsImageMakerOptions;
use crate::mst_graph_drawer;
use crate::graph_drawer;

pub fn draw_mst(options: &RandomPointsImageMakerOptions) {
    type PointNumberType = f64;
    const DIMS: usize = 2;
    let points_range = (options.points_min_range, options.points_max_range);

    let (mst, euclidean_distance_graph) = create_mst_from_random_points::<PointNumberType, DIMS>(
        options.number_of_points, points_range);

    let mut mst_vertex_iter = mst_graph_drawer::drawable_vertex_iter(
        &mst, &euclidean_distance_graph.point_list);
    let mut mst_edge_iter = mst_graph_drawer::drawable_edge_iter(
        &mst, &euclidean_distance_graph.point_list);
    let max_depth = mst.traverse_nodes_preorder().map(|(_, l)| l).max().unwrap_or(0);

    let drawing_properties = graph_drawer::DrawingProperties {
        image_width: options.result_image_width,
        image_height: options.result_image_height,
        x_min: options.result_image_bounds_min,
        x_max: options.result_image_bounds_max,
        y_min: options.result_image_bounds_min,
        y_max: options.result_image_bounds_max,
        max_depth
    };
    let image = graph_drawer::draw_2d_graph(
        &mut mst_vertex_iter, &mut mst_edge_iter, &drawing_properties);

    println!("Saving image to {}", &options.result_image_path
                                           .clone()
                                           .into_os_string()
                                           .into_string()
                                           .unwrap_or("".to_string()));
    image.save(&options.result_image_path).unwrap();
}
