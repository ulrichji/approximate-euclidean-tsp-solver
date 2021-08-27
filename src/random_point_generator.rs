use rand::distributions::{Distribution, Uniform, uniform::SampleUniform};

use crate::euclidean_distance_graph::{Point};

pub trait Zero {
    fn zero() -> Self;
}

impl Zero for f64 {
    fn zero() -> f64 { 0f64 }
}

pub struct RandomPointGenerator<T, const N: usize> where T: SampleUniform + Copy {
    random_number_generator: rand::rngs::ThreadRng,
    ranges: [(T, T); N]
}

impl<T, const N: usize> RandomPointGenerator<T, N> where T: SampleUniform + Copy {
    pub fn new_square_range(min_range: T, max_range: T) -> RandomPointGenerator<T, N> {
        RandomPointGenerator {
            random_number_generator: rand::thread_rng(),
            ranges: [(min_range, max_range); N]
        }
    }
}

impl<T, const N: usize> Iterator for RandomPointGenerator<T, N>
        where T: SampleUniform + Copy + Zero {
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
