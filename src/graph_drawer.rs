use imageproc::drawing::draw_line_segment_mut;
use imageproc::drawing::draw_filled_circle_mut;

use image::Rgb;
use image::RgbImage;

use palette::{FromColor, Pixel};

pub struct DrawableVertex2D {
    x: f64,
    y: f64,
    depth: usize
}

impl DrawableVertex2D {
    pub fn new(x: f64, y: f64, depth: usize) -> DrawableVertex2D { DrawableVertex2D { x, y, depth } }
}

pub struct DrawableEdge2D {
    from: DrawableVertex2D,
    to: DrawableVertex2D
}

impl DrawableEdge2D {
    pub fn new(from: DrawableVertex2D, to: DrawableVertex2D) -> DrawableEdge2D { DrawableEdge2D { from, to } }
}

pub struct DrawingProperties {
    pub image_width: u32,
    pub image_height: u32,
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub max_depth: usize
}

pub fn draw_2d_graph<
        'a,
        V: 'a + Iterator<Item=DrawableVertex2D>,
        E: 'a + Iterator<Item=DrawableEdge2D>>(
            vertex_iterator: V,
            edge_iterator: E,
            properties: &DrawingProperties) -> RgbImage {
    let mut img = RgbImage::new(properties.image_width, properties.image_height);

    let width = properties.x_max - properties.x_min;
    let height = properties.y_max - properties.y_min;

    for vertex in vertex_iterator {
        let vertex_center_x = (((vertex.x - properties.x_min) / width) * properties.image_width as f64).round() as i32;
        let vertex_center_y = (((vertex.y - properties.y_min) / height) * properties.image_height as f64).round() as i32;

        let color_hue = palette::Hsv::new(120.0, 0.1, 0.3);
        let color_rgb = palette::Srgb::from_color(color_hue);
        let color_pixel: [u8; 3] = color_rgb.into_format().into_raw();

        draw_filled_circle_mut(&mut img, (vertex_center_x, vertex_center_y), 1, Rgb(color_pixel));
    }

    for edge in edge_iterator {
        let from_vertex_center_x = (((edge.from.x - properties.x_min) / width) * properties.image_width as f64).round() as f32;
        let from_vertex_center_y = (((edge.from.y - properties.y_min) / height) * properties.image_height as f64).round() as f32;

        let to_vertex_center_x = (((edge.to.x - properties.x_min) / width) * properties.image_width as f64).round() as f32;
        let to_vertex_center_y = (((edge.to.y - properties.y_min) / height) * properties.image_height as f64).round() as f32;

        let color_hue_angle = (edge.to.depth as f64 / properties.max_depth as f64) * 360.0;
        let color_hue = palette::Hsv::new(color_hue_angle, 1.0, 1.0);
        let color_rgb = palette::Srgb::from_color(color_hue);
        let color_pixel: [u8; 3] = color_rgb.into_format().into_raw();

        draw_line_segment_mut(&mut img, (from_vertex_center_x, from_vertex_center_y), (to_vertex_center_x, to_vertex_center_y), Rgb(color_pixel));
    }

    img
}
