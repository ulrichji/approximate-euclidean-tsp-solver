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
impl<V, T> Ord for Edge<V, T> where T: PartialOrd {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(&other).unwrap()
    }
}

pub trait Graph<V, T> where T: PartialOrd {
    fn get_root_vertex_identifier(&self) -> V;
    fn next_neighbour_for_vertex(&mut self, vertex_identifier: &V) -> Option<Edge<V, T>>;
}

pub enum TreeSizeLimit {
    Finite(usize),
    Infinite
}

#[derive(Debug)]
pub struct MinimumSpanningTree<T> {
    pub elements: Vec<Vertex<T>>
}

struct MinimumSpanningTreeBuilder<V, T: std::cmp::PartialOrd> {
    tree: MinimumSpanningTree<V>,
    vertices_in_spanning_tree: HashMap<V, usize>,
    edge_priority_queue: BinaryHeap<Reverse<Edge<V, T>>>
}

impl <V, T: std::cmp::PartialOrd> MinimumSpanningTreeBuilder <V, T>
        where V: PartialEq + Copy + std::cmp::Eq + std::hash::Hash,
              T: std::cmp::PartialOrd {

    fn new_with_root(root_vertex_identifier: V) -> MinimumSpanningTreeBuilder <V, T> {
        let mut tree = MinimumSpanningTree::new();
        let mut vertices_in_spanning_tree = HashMap::new();

        let root_vertex = Vertex::new(root_vertex_identifier);
        tree.add_vertex(root_vertex);
        vertices_in_spanning_tree.insert(root_vertex_identifier, 0);

        MinimumSpanningTreeBuilder {
            tree,
            vertices_in_spanning_tree,
            edge_priority_queue: BinaryHeap::new()
        }
    }

    fn construct_from_graph(&mut self, graph: &mut dyn Graph<V, T>, max_tree_size: TreeSizeLimit) {
        self.initialize_edge_priority_queue(graph);

        while !self.is_done(&max_tree_size) {
            let Reverse(edge) = self.edge_priority_queue.pop().unwrap();

            if self.is_vertex_undiscovered(&edge.to_vertex_identifier) {
                let child_vertex = Vertex::new(edge.to_vertex_identifier);
                self.add_vertex_to_mst(child_vertex);
                self.add_child_to_parent(&edge.from_vertex_identifier, &edge.to_vertex_identifier);
                self.push_next_edge_for_vertex_into_queue(graph, &edge.to_vertex_identifier);
            }

            // We have taken the edge we know for sure is the shortest edge in/from the mst set.
            // Now, the next might also come from the same vertex, so we add the next edge from
            // that same vertex into the queue.
            self.push_next_edge_for_vertex_into_queue(graph, &edge.from_vertex_identifier);
        }
    }

    fn initialize_edge_priority_queue(&mut self, graph: &mut dyn Graph<V, T>) {
        let root_element_index = 0;
        let root_vertex_identifier = &self.tree.get_vertex_from_index(root_element_index)
                                               .vertex_identifier;

        let root_neighbour_vertex_identifier =
            graph.next_neighbour_for_vertex(root_vertex_identifier);
        match root_neighbour_vertex_identifier {
            Some(v) => self.edge_priority_queue.push(Reverse(v)),
            None => {  }
        };
    }

    fn is_done(&self, max_tree_size: &TreeSizeLimit) -> bool {
        self.edge_priority_queue.is_empty() ||
        match max_tree_size {
            TreeSizeLimit::Finite(n) => {
                let number_of_found_vertices = self.tree.elements.len();
                number_of_found_vertices >= *n
            },
            TreeSizeLimit::Infinite => false
        }
    }

    fn is_vertex_undiscovered(&self, vertex_identifier: &V) -> bool {
        !self.is_vertex_discovered(vertex_identifier)
    }

    fn is_vertex_discovered(&self, vertex_identifier: &V) -> bool {
        self.vertices_in_spanning_tree.contains_key(vertex_identifier)
    }

    fn add_vertex_to_mst(&mut self, vertex: Vertex<V>) {
        let vertex_index = self.tree.elements.len();
        self.vertices_in_spanning_tree.insert(vertex.vertex_identifier, vertex_index);
        self.tree.add_vertex(vertex);
    }

    fn add_child_to_parent(&mut self, parent_vertex_identifier: &V, child_vertex_identifier: &V) {
        let parent_index = self.vertices_in_spanning_tree.get(parent_vertex_identifier).unwrap();
        let child_index = self.vertices_in_spanning_tree.get(child_vertex_identifier).unwrap();
        self.tree.add_child_index(*parent_index, *child_index);
    }

    fn push_next_edge_for_vertex_into_queue(&mut self,
                                            graph: &mut dyn Graph<V, T>,
                                            vertex_identifier: &V) {
        let next_edge_from_vertex = graph.next_neighbour_for_vertex(vertex_identifier);
        match next_edge_from_vertex {
            Some(e) => self.edge_priority_queue.push(Reverse(e)),
            None => {}
        }
    }

    fn take_spanning_tree(&mut self) -> MinimumSpanningTree<V> {
        let tree = std::mem::replace(&mut self.tree, MinimumSpanningTree::new());
        self.vertices_in_spanning_tree.clear();
        self.edge_priority_queue.clear();

        tree
    }
}

impl<V> MinimumSpanningTree<V> where V: PartialEq + Copy + std::cmp::Eq + std::hash::Hash {
    fn new() -> MinimumSpanningTree<V> {
        MinimumSpanningTree {
            elements: Vec::new()
        }
    }

    pub fn construct<T>(graph: &mut dyn Graph<V, T>, max_tree_size: TreeSizeLimit)
            -> MinimumSpanningTree<V> where T: PartialOrd {
        let mut tree_builder: MinimumSpanningTreeBuilder<V, T> =
            MinimumSpanningTreeBuilder::new_with_root(graph.get_root_vertex_identifier());
        tree_builder.construct_from_graph(graph, max_tree_size);
        tree_builder.take_spanning_tree()
    }

    fn add_vertex(&mut self, vertex: Vertex<V>) {
        self.elements.push(vertex);
    }

    fn add_child_index(&mut self, parent_index: usize, child_index: usize) {
        self.add_child_index_to_parent(parent_index, child_index);
    }

    fn add_child_index_to_parent(&mut self, parent_index: usize, child_index: usize) {
        match self.get_vertex_from_index(parent_index).first_child_index {
            Some(i) => { self.add_sibling_index_to_sibling(i, child_index) }
            None => { self.get_vertex_from_index_mut(parent_index)
                          .set_first_child_index(child_index) }
        }
    }

    fn add_sibling_index_to_sibling(&mut self, vertex_index: usize, sibling_index: usize) {
        match self.get_vertex_from_index(vertex_index).next_sibling_index {
            Some(i) => { self.add_sibling_index_to_sibling(i, sibling_index) }
            None => { self.get_vertex_from_index_mut(vertex_index)
                          .next_sibling_index = Some(sibling_index) }
        }
    }

    fn get_vertex_from_index(&self, index: usize) -> &Vertex<V> {
        &self.elements[index]
    }

    fn get_vertex_from_index_mut(&mut self, index: usize) -> &mut Vertex<V> {
         &mut self.elements[index]
    }
}
