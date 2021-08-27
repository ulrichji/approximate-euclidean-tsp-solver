use imageproc;
use image;
use palette::{Pixel};

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

struct PixelPosition {
    x: i32,
    y: i32
}

impl PixelPosition {
    fn as_tuple<T>(&self) -> (T, T) where T: std::convert::From<i32> {
        (self.x.into(), self.y.into())
    }
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

impl DrawingProperties {
    fn compute_vertex_pixel_coordinate(&self, vertex: &DrawableVertex2D) -> PixelPosition {
        let vertex_in_image_clamped = DrawableVertex2D {
            x: ((vertex.x - self.x_min) / self.get_drawing_field_width()),
            y: ((vertex.y - self.y_min) / self.get_drawing_field_height()),
            depth: vertex.depth };

        let vertex_pixel_pos = DrawableVertex2D {
            x: (vertex_in_image_clamped.x * self.image_width as f64).round(),
            y: (vertex_in_image_clamped.y * self.image_height as f64).round(),
            depth: vertex_in_image_clamped.depth };

        PixelPosition {
            x: vertex_pixel_pos.x as i32,
            y: vertex_pixel_pos.y as i32
        }
    }

    fn get_drawing_field_width(&self) -> f64 {
        self.x_max - self.x_min
    }

    fn get_drawing_field_height(&self) -> f64 {
        self.y_max - self.y_min
    }
}

pub fn draw_2d_graph<'a, V, E>(
    vertex_iterator: V, edge_iterator: E, properties: &DrawingProperties) -> image::RgbImage
        where V: 'a + Iterator<Item=DrawableVertex2D>,
              E: 'a + Iterator<Item=DrawableEdge2D> {
    let mut img = image::RgbImage::new(properties.image_width, properties.image_height);

    for vertex in vertex_iterator {
        let color_hue_angle = 120.0;
        let color_hsv = palette::Hsv::new(color_hue_angle, 0.1, 0.3);
        draw_vertex(&mut img, &vertex, &properties, color_hsv);
    }

    for edge in edge_iterator {
        let color_hue_angle = (edge.to.depth as f64 / properties.max_depth as f64) * 360.0;
        let color_hvs = palette::Hsv::new(color_hue_angle as f32, 1.0, 1.0);
        draw_edge(&mut img, &edge, &properties, color_hvs);
    }

    img
}

pub fn draw_vertex<C>(img: &mut image::RgbImage,
                      vertex: &DrawableVertex2D,
                      properties: &DrawingProperties,
                      color: C) where C: palette::convert::IntoColor<palette::Srgb> {
    let pixel_position = properties.compute_vertex_pixel_coordinate(vertex).as_tuple();

    let color_rgb: palette::Srgb = color.into_color();
    let color_pixel: [u8; 3] = color_rgb.into_format().into_raw();

    imageproc::drawing::draw_filled_circle_mut(img, pixel_position, 1, image::Rgb(color_pixel));
}

pub fn draw_edge<C>(img: &mut image::RgbImage,
                    edge: &DrawableEdge2D,
                    properties: &DrawingProperties,
                    color: C) where C: palette::convert::IntoColor<palette::Srgb> {
    let from_pixel_position = properties.compute_vertex_pixel_coordinate(&edge.from).as_tuple();
    let to_pixel_position = properties.compute_vertex_pixel_coordinate(&edge.to).as_tuple();

    let color_rgb: palette::Srgb = color.into_color();
    let color_pixel: [u8; 3] = color_rgb.into_format().into_raw();

    imageproc::drawing::draw_antialiased_line_segment_mut(img,
                                                          from_pixel_position,
                                                          to_pixel_position,
                                                          image::Rgb(color_pixel),
                                                          &imageproc::pixelops::interpolate);
}
