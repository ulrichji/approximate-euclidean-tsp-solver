use std::collections::HashMap;
use std::collections::BinaryHeap;
use std::cmp::Reverse;
pub use std::cmp::Ord;

#[derive(Debug)]
pub struct Vertex<V> {
    pub first_child_index: Option<usize>,
    pub next_sibling_index: Option<usize>,
    pub vertex_identifier: V
}

impl<V> Vertex<V> {
    pub fn new(vertex_identifier: V) -> Vertex<V> {
        Vertex {
            first_child_index: None,
            next_sibling_index: None,
            vertex_identifier
        }
    }

    fn set_first_child_index(&mut self, first_child_index: usize) {
        self.first_child_index = Some(first_child_index);
    }

    fn is_leaf(&self) -> bool {
        self.first_child_index.is_none()
    }

    fn is_stem(&self) -> bool {
        self.first_child_index.is_some()
    }
}

pub struct Edge<V, T> where T: PartialOrd {
    from_vertex_identifier: V,
    to_vertex_identifier: V,
    weight: T
}

impl<V, T> Edge<V, T> where T: PartialOrd {
    pub fn new(from_vertex_identifier: V, to_vertex_identifier: V, weight: T) -> Edge<V, T> {
        Edge::<V, T>{ from_vertex_identifier, to_vertex_identifier, weight }
    }
}

impl<V, T> PartialEq for Edge<V, T> where T: PartialOrd {
    fn eq(&self, other: &Self) -> bool { self.weight == other.weight }
}
impl<V, T> Eq for Edge<V, T> where T: PartialOrd { }

impl<V, T> PartialOrd for Edge<V, T> where T: PartialOrd {
    fn partial_cmp(&self, other: &Self) -> std::option::Option<std::cmp::Ordering> {
        self.weight.partial_cmp(&other.weight)
    }
}
impl<V, T> Ord for Edge<V, T> where T: Ord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.weight.cmp(&other.weight)
    }
}

pub trait Graph<V, T> where T: Ord {
    fn get_root_vertex_identifier(&self) -> V;
    fn next_neighbour_for_vertex(&mut self, vertex_identifier: &V) -> Option<Edge<V, T>>;
}

pub enum TreeSize {
    Finite(usize),
    Infinite
}

#[derive(Debug)]
struct MinimumSpanningTreeElement<T> {
    first_child_index: Option<usize>,
    next_sibling_index: Option<usize>,
    data: T
}

#[derive(Debug)]
pub struct MinimumSpanningTree<T> {
    pub elements: Vec<Vertex<T>>
}

impl<V> MinimumSpanningTree<V> where V: PartialEq + Copy + std::cmp::Eq + std::hash::Hash {
    fn new() -> MinimumSpanningTree<V> {
        MinimumSpanningTree {
            elements: Vec::new()
        }
    }

    pub fn construct<T>(graph: &mut dyn Graph<V, T>, max_tree_size: TreeSize) -> MinimumSpanningTree<V>
            where T: Ord {
        let mut minimum_spanning_tree = Self::new();
        minimum_spanning_tree.elements.push(Vertex::new(graph.get_root_vertex_identifier()));

        let mut vertex_set_in_spanning_tree = HashMap::new();
        vertex_set_in_spanning_tree.insert(graph.get_root_vertex_identifier(), 0);
        let mut edge_priority_queue = BinaryHeap::new();

        let root_vertex_identifier = graph.next_neighbour_for_vertex(&minimum_spanning_tree.elements[0].vertex_identifier);
        match root_vertex_identifier {
            Some(v) => edge_priority_queue.push(Reverse(v)),
            None => return MinimumSpanningTree{ elements: Vec::new() }
        };

        while !edge_priority_queue.is_empty() && match &max_tree_size { TreeSize::Finite(n) => minimum_spanning_tree.elements.len() < *n, TreeSize::Infinite => true } {
            let Reverse(edge) = edge_priority_queue.pop().unwrap();

            if !vertex_set_in_spanning_tree.contains_key(&edge.to_vertex_identifier) {
                minimum_spanning_tree.elements.push(Vertex::new(edge.to_vertex_identifier));
                let child_index = minimum_spanning_tree.elements.len() - 1;
                vertex_set_in_spanning_tree.insert(edge.to_vertex_identifier, child_index);

                let parent_vertex_index = vertex_set_in_spanning_tree.get(&edge.from_vertex_identifier).unwrap();
                minimum_spanning_tree.add_child_index(*parent_vertex_index, child_index);

                let next_edge_from_vertex = graph.next_neighbour_for_vertex(&edge.to_vertex_identifier);
                match next_edge_from_vertex {
                    Some(e) => edge_priority_queue.push(Reverse(e)),
                    None => {}
                }
            }

            // We have taken the edge we know for sure is the shortest edge in/from the mst set.
            // Now, the next might also come from the same vertex, so we add the next edge from
            // that same vertex into the queue.
            let next_edge_from_vertex = graph.next_neighbour_for_vertex(&edge.from_vertex_identifier);
            match next_edge_from_vertex {
                Some(e) => edge_priority_queue.push(Reverse(e)),
                None => {}
            }
        }

        minimum_spanning_tree
    }

    fn add_child_index(&mut self, parent_index: usize, child_index: usize) {
        self.add_child_index_to_parent(parent_index, child_index);
    }

    fn add_child_index_to_parent(&mut self, parent_index: usize, child_index: usize) {
        match self.elements[parent_index].first_child_index {
            Some(i) => { self.add_sibling_index_to_sibling(i, child_index) }
            None => { self.elements[parent_index].set_first_child_index(child_index) }
        }
    }

    fn add_sibling_index_to_sibling(&mut self, vertex_index: usize, sibling_index: usize) {
        match self.elements[vertex_index].next_sibling_index {
            Some(i) => { self.add_sibling_index_to_sibling(i, sibling_index) }
            None => { self.elements[vertex_index].next_sibling_index = Some(sibling_index) }
        }
    }
}
