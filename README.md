# approximate-euclidean-tsp-solver
A minimum spanning tree approximation of the Euclidean travelling salesman problem using a kd-tree in Rust. It also generates cool looking images for the 2D case.

## Requirements
A Rust build environment is required, the following version is was tested:
 - cargo: 1.54.0
 - rustc: 1.54.0

## Building
Standard cargo build commands can be used to build:
```
cargo build --release
```

## Running
The following command will approximate the Euclidean tsp problem using 10000 points, and save it in a file called `tsp.png`
```
cargo run --release -- tsp tsp.png -n 10000
```

Alternatively, you can save the minimum spanning tree in the `mst.png` file using the following command:
```
cargo run --release -- mst mst.png -n 10000
```
