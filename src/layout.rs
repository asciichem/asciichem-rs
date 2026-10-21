//! Deterministic 2D layout for bonded molecules — the Rust port of
//! the TypeScript implementation's `layout.ts` (itself a port of the
//! reference decomposition): a pure walker produces a neutral
//! atom+bond graph (ring-closure pairing mirrors RingBonds), then a
//! chemistry-shaped deterministic placer — rings as regular polygons
//! (fused rings anchored on shared atoms), chains as 120-degree
//! zigzags grown breadth-first in source order.

use crate::model::{self, Node};

pub const BOND_LENGTH: f64 = 56.0;

#[derive(Debug, Clone, PartialEq)]
pub struct RingBond {
    pub digit: String,
    pub from_atom: usize,
    pub to_atom: usize,
}

#[derive(Debug, Clone)]
pub struct LayoutBond {
    pub from: usize,
    pub to: usize,
    pub bond: model::Bond,
}

#[derive(Debug, Clone, Default)]
pub struct LayoutGraph {
    pub atoms: Vec<model::Atom>,
    pub bonds: Vec<LayoutBond>,
    pub ring_bonds: Vec<RingBond>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

// Ring-closure pairing, ported from RingBonds: walk atoms in source
// order; a repeated digit closes the ring between the opener and the
// closer. Multiple digits per atom open in parallel.
pub fn ring_bonds(atoms: &[model::Atom]) -> Vec<RingBond> {
    let mut open: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut rings = Vec::new();
    for (index, atom) in atoms.iter().enumerate() {
        let Some(closures) = &atom.ring_closures else {
            continue;
        };
        for digit in closures.chars() {
            let key: String = digit.to_string();
            if let Some(opener) = open.remove(key.as_str()) {
                rings.push(RingBond {
                    digit: key,
                    from_atom: opener,
                    to_atom: index,
                });
            } else {
                open.insert(key, index);
            }
        }
    }
    rings
}

pub fn unclosed_ring_atoms(atoms: &[model::Atom]) -> Vec<usize> {
    let mut open: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        let Some(closures) = &atom.ring_closures else {
            continue;
        };
        for digit in closures.chars() {
            let key = digit.to_string();
            if open.remove(key.as_str()).is_none() {
                open.insert(key, index);
            }
        }
    }
    let mut unclosed: Vec<usize> = open.into_values().collect();
    unclosed.sort_unstable();
    unclosed
}

// Walks a molecule into a neutral atom+bond graph. Bond tokens
// connect the last atom of the previous unit to the first atom of
// the next; nested groups contribute their walk order.
pub fn walk_molecule(molecule: &model::Molecule) -> LayoutGraph {
    let mut graph = LayoutGraph::default();

    fn walk(node: &Node, atoms: &mut Vec<model::Atom>) -> Option<(usize, usize)> {
        match node {
            Node::Atom(atom) => {
                atoms.push(atom.clone());
                let last = atoms.len() - 1;
                Some((last, last))
            }
            Node::Group(group) => {
                let mut first = None;
                let mut last = None;
                for child in &group.nodes {
                    if let Some((child_first, child_last)) = walk(child, atoms) {
                        first.get_or_insert(child_first);
                        last = Some(child_last);
                    }
                }
                Some((first?, last?))
            }
            _ => None,
        }
    }

    let mut previous = None;
    let mut pending_bond = None;
    for unit in &molecule.nodes {
        match unit {
            Node::Bond(bond) => {
                pending_bond = Some(bond.clone());
            }
            other => {
                if let Some((first, last)) = walk(other, &mut graph.atoms) {
                    if let (Some(prev), Some(bond)) = (previous, pending_bond.take()) {
                        graph.bonds.push(LayoutBond {
                            from: prev,
                            to: first,
                            bond,
                        });
                    }
                    previous = Some(last);
                }
            }
        }
    }

    let rings = ring_bonds(&graph.atoms);
    for ring in &rings {
        graph.bonds.push(LayoutBond {
            from: ring.from_atom,
            to: ring.to_atom,
            bond: model::Bond {
                kind: "single".to_string(),
                direction: None,
            },
        });
    }
    graph.ring_bonds = rings;
    graph
}

fn neighbors_of(graph: &LayoutGraph, node: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for bond in &graph.bonds {
        if bond.from == node {
            out.push(bond.to);
        } else if bond.to == node {
            out.push(bond.from);
        }
    }
    out
}

// Shortest path between the closure endpoints, excluding the closure
// edge itself — the path plus the closure edge is the ring cycle.
fn ring_cycle(graph: &LayoutGraph, closure: &RingBond) -> Option<Vec<usize>> {
    let mut queue: std::collections::VecDeque<Vec<usize>> = std::collections::VecDeque::new();
    queue.push_back(vec![closure.from_atom]);
    let mut visited = std::collections::HashSet::new();
    visited.insert(closure.from_atom);
    while let Some(path) = queue.pop_front() {
        let node = *path.last().expect("seeded");
        if node == closure.to_atom {
            return Some(path);
        }
        for neighbor in neighbors_of(graph, node) {
            let is_closure_edge = (node == closure.from_atom && neighbor == closure.to_atom)
                || (node == closure.to_atom && neighbor == closure.from_atom);
            if is_closure_edge || visited.contains(&neighbor) {
                continue;
            }
            visited.insert(neighbor);
            let mut next = path.clone();
            next.push(neighbor);
            queue.push_back(next);
        }
    }
    None
}

pub fn place_graph(graph: &LayoutGraph) -> Vec<Point> {
    let mut positions: Vec<Option<Point>> = vec![None; graph.atoms.len()];

    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for closure in &graph.ring_bonds {
        if let Some(path) = ring_cycle(graph, closure) {
            if path.len() >= 3 {
                cycles.push(path);
            }
        }
    }

    // Place rings first: each ring is a regular polygon. A ring with
    // an already-placed member anchors on it (fused look); a ring in
    // a new component starts at the next free center to the right.
    let mut free_center_x = 0.0_f64;
    for cycle in &cycles {
        let anchor_index = cycle.iter().position(|member| positions[*member].is_some());
        let radius = polygon_radius(cycle.len());
        match anchor_index {
            None => {
                let center = Point {
                    x: free_center_x,
                    y: 0.0,
                };
                place_polygon(cycle, &center, radius, 0.0, &mut positions);
                free_center_x += radius * 2.0 + BOND_LENGTH;
            }
            Some(anchor_index) => {
                let anchor = cycle[anchor_index];
                let anchor_pos = positions[anchor].expect("anchored on placed");
                // Direction pointing away from any already-placed
                // neighbor of the anchor; the polygon center sits
                // opposite that direction.
                let mut outward = f64::NEG_INFINITY;
                for neighbor in neighbors_of(graph, anchor) {
                    if let Some(npos) = positions[neighbor] {
                        outward = (anchor_pos.y - npos.y).atan2(anchor_pos.x - npos.x);
                        break;
                    }
                }
                let vertex_angle = angle_of_member(cycle.len(), anchor_index);
                let center_angle = if outward == f64::NEG_INFINITY {
                    vertex_angle + std::f64::consts::PI
                } else {
                    outward
                };
                let center = Point {
                    x: anchor_pos.x + radius * center_angle.cos(),
                    y: anchor_pos.y + radius * center_angle.sin(),
                };
                // Rotate the polygon so the anchor lands on its own
                // position.
                let shift = center_angle - vertex_angle;
                place_polygon(cycle, &center, radius, shift, &mut positions);
            }
        }
    }

    // Grow chains breadth-first from every placed atom (zigzag), then
    // from any still-unplaced atom (new component at origin baseline).
    for i in 0..graph.atoms.len() {
        if positions[i].is_none() {
            continue;
        }
        let mut visited = std::collections::HashSet::new();
        visited.insert(i);
        grow(i, &mut positions, graph, &mut visited);
    }
    for i in 0..graph.atoms.len() {
        if positions[i].is_some() {
            continue;
        }
        positions[i] = Some(Point {
            x: free_center_x,
            y: 0.0,
        });
        free_center_x += BOND_LENGTH;
        let mut visited = std::collections::HashSet::new();
        visited.insert(i);
        grow(i, &mut positions, graph, &mut visited);
    }

    positions
        .into_iter()
        .map(|p| p.expect("every atom placed"))
        .collect()
}

fn grow(
    start: usize,
    positions: &mut [Option<Point>],
    graph: &LayoutGraph,
    visited: &mut std::collections::HashSet<usize>,
) {
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(start);
    while let Some(current) = queue.pop_front() {
        let Some(current_pos) = positions[current] else {
            continue;
        };
        // Direction the chain is travelling: away from the origin of
        // the component (0,0) — keeps zigzags deterministic and
        // horizontal-ish.
        let heading = (-current_pos.y).atan2(current_pos.x);
        let children: Vec<usize> = neighbors_of(graph, current)
            .into_iter()
            .filter(|n| !visited.contains(n) && positions[*n].is_none())
            .collect();
        for (i, child) in children.iter().enumerate() {
            if visited.contains(child) {
                continue;
            }
            let spread = match i {
                0 => -std::f64::consts::FRAC_PI_3,
                1 => std::f64::consts::FRAC_PI_3,
                _ => 0.0,
            };
            let angle = heading + spread;
            positions[*child] = Some(Point {
                x: current_pos.x + BOND_LENGTH * angle.cos(),
                y: current_pos.y + BOND_LENGTH * angle.sin(),
            });
            visited.insert(*child);
            queue.push_back(*child);
        }
    }
}

// Circumradius for a regular polygon with the given side length.
fn polygon_radius(sides: usize) -> f64 {
    BOND_LENGTH / (2.0 * (std::f64::consts::PI / sides as f64).sin())
}

fn angle_of_member(sides: usize, member_index: usize) -> f64 {
    std::f64::consts::PI * 2.0 * member_index as f64 / sides as f64 - std::f64::consts::FRAC_PI_2
}

fn place_polygon(
    members: &[usize],
    center: &Point,
    radius: f64,
    shift: f64,
    positions: &mut [Option<Point>],
) {
    for (i, atom) in members.iter().enumerate() {
        let angle = angle_of_member(members.len(), i) + shift;
        positions[*atom] = Some(Point {
            x: center.x + radius * angle.cos(),
            y: center.y + radius * angle.sin(),
        });
    }
}

pub fn layout_molecule(molecule: &model::Molecule) -> (LayoutGraph, Vec<Point>) {
    let graph = walk_molecule(molecule);
    let positions = place_graph(&graph);
    (graph, positions)
}
