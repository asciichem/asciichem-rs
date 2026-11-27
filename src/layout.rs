//! Deterministic 2D layout for bonded molecules — the Rust port of
//! the TypeScript implementation's `layout.ts`: a pure walker produces
//! a neutral atom+bond graph (ring-closure pairing mirrors RingBonds),
//! then a chemistry-shaped deterministic placer — rings as regular
//! polygons (free, fused on placed members, anchored outward from the
//! placed centroid), chains as turn-alternating 120-degree zigzags
//! grown breadth-first with occupied-spot avoidance.

use crate::model::{self, Node};

pub const BOND_LENGTH: f64 = 56.0;
const COLLISION_RADIUS: f64 = BOND_LENGTH / 2.0;
const SPREAD: f64 = std::f64::consts::FRAC_PI_3;

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

// Unclosed ring atoms in first-opened order (the Ruby reference's
// insertion-ordered hash walk, not a sort).
pub fn unclosed_ring_atoms(atoms: &[model::Atom]) -> Vec<usize> {
    let mut open: Vec<(String, usize)> = Vec::new();
    for (index, atom) in atoms.iter().enumerate() {
        let Some(closures) = &atom.ring_closures else {
            continue;
        };
        for digit in closures.chars() {
            let key = digit.to_string();
            if let Some(position) = open.iter().position(|(digit, _)| *digit == key) {
                open.remove(position);
            } else {
                open.push((key, index));
            }
        }
    }
    open.into_iter().map(|(_, index)| index).collect()
}

// Walks a molecule into a neutral atom+bond graph. Bond tokens
// connect the last atom of the previous unit to the first atom of
// the next; nested groups contribute their walk order. A group whose
// first node (unwrapping nested molecules) is a bond is a BRANCH
// group: the leading bond attaches the branch to the atom before the
// group, and the continuation resumes from that root — `C(-O)=O`
// reads C-O and C=O, never O=O.
pub fn walk_molecule(molecule: &model::Molecule) -> LayoutGraph {
    let mut graph = LayoutGraph::default();
    let mut state = WalkState::default();
    walk_nodes(&molecule.nodes, &mut graph, &mut state);

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

#[derive(Default)]
struct WalkState {
    previous: Option<usize>,
    pending_bond: Option<model::Bond>,
}

fn walk_nodes(nodes: &[Node], graph: &mut LayoutGraph, state: &mut WalkState) {
    for unit in nodes {
        match unit {
            Node::Bond(bond) => {
                state.pending_bond = Some(bond.clone());
            }
            Node::Atom(atom) => {
                graph.atoms.push(atom.clone());
                let index = graph.atoms.len() - 1;
                if let (Some(previous), Some(bond)) = (state.previous, state.pending_bond.take()) {
                    graph.bonds.push(LayoutBond {
                        from: previous,
                        to: index,
                        bond,
                    });
                }
                state.previous = Some(index);
                state.pending_bond = None;
            }
            Node::Group(group) => {
                let branch_root = group
                    .nodes
                    .first()
                    .and_then(leading_atom_bond)
                    .and(state.previous);
                walk_nodes(&group.nodes, graph, state);
                if let Some(root) = branch_root {
                    state.previous = Some(root);
                }
            }
            Node::Molecule(molecule) => {
                walk_nodes(&molecule.nodes, graph, state);
            }
            _ => {}
        }
    }
}

// The branch test: the container's first node, unwrapping nested
// molecules, is a bond.
fn leading_atom_bond(first: &Node) -> Option<()> {
    let mut first = first;
    while let Node::Molecule(molecule) = first {
        first = molecule.nodes.first()?;
    }
    match first {
        Node::Bond(_) => Some(()),
        _ => None,
    }
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

struct RingCycle {
    members: Vec<usize>,
}

// Rings first, in source order. Each cycle is placed free (nothing
// placed yet), fused (two or more atoms already placed — fit the
// polygon to them), or anchored (one placed atom — orient away from
// the placed region).
pub fn place_graph(graph: &LayoutGraph) -> Vec<Point> {
    let mut positions: Vec<Option<Point>> = vec![None; graph.atoms.len()];

    let mut cycles: Vec<RingCycle> = Vec::new();
    for closure in &graph.ring_bonds {
        if let Some(path) = ring_cycle(graph, closure) {
            if path.len() >= 3 {
                cycles.push(RingCycle { members: path });
            }
        }
    }

    let mut free_center_x = 0.0_f64;
    for cycle in &cycles {
        let anchored = cycle
            .members
            .iter()
            .filter(|member| positions[**member].is_some())
            .count();
        if anchored == 0 {
            free_center_x = place_free_cycle(&cycle.members, &mut positions, free_center_x);
        } else if anchored >= 2 {
            place_fused_cycle(&cycle.members, &mut positions);
        } else {
            place_anchored_cycle(&cycle.members, &mut positions);
        }
    }

    // Grow chains breadth-first from every placed atom (zigzag), then
    // from any still-unplaced atom (new component at origin baseline).
    for i in 0..graph.atoms.len() {
        if positions[i].is_none() {
            continue;
        }
        grow(i, &mut positions, graph);
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
        grow(i, &mut positions, graph);
    }

    positions
        .into_iter()
        .map(|p| p.expect("every atom placed"))
        .collect()
}

fn place_free_cycle(
    members: &[usize],
    positions: &mut [Option<Point>],
    free_center_x: f64,
) -> f64 {
    let radius = polygon_radius(members.len());
    place_polygon(
        members,
        &Point {
            x: free_center_x,
            y: 0.0,
        },
        radius,
        0.0,
        positions,
    );
    free_center_x + radius * 2.0 + BOND_LENGTH
}

// Fused on shared atoms: derive the rigid transform (rotation +
// translation) that maps the ideal polygon's two placed members onto
// their laid-out positions, then fill only the unplaced members.
// Placed atoms never move — the first ring keeps its geometry.
fn place_fused_cycle(members: &[usize], positions: &mut [Option<Point>]) {
    let a = *members
        .iter()
        .find(|member| positions[**member].is_some())
        .expect("fused ring has placed members");
    let b = *members
        .iter()
        .find(|member| positions[**member].is_some() && **member != a)
        .expect("fused ring has two placed members");
    let ia = members.iter().position(|member| *member == a).expect("member");
    let ib = members.iter().position(|member| *member == b).expect("member");
    let ideal_a = ideal_point(members.len(), ia, 0.0);
    let ideal_b = ideal_point(members.len(), ib, 0.0);
    let a_pos = positions[a].expect("placed");
    let b_pos = positions[b].expect("placed");
    let theta = (b_pos.y - a_pos.y).atan2(b_pos.x - a_pos.x)
        - (ideal_b.y - ideal_a.y).atan2(ideal_b.x - ideal_a.x);
    let cos = theta.cos();
    let sin = theta.sin();
    for (i, member) in members.iter().enumerate() {
        if positions[*member].is_some() {
            continue;
        }
        let ideal = ideal_point(members.len(), i, 0.0);
        let dx = ideal.x - ideal_a.x;
        let dy = ideal.y - ideal_a.y;
        positions[*member] = Some(Point {
            x: a_pos.x + dx * cos - dy * sin,
            y: a_pos.y + dx * sin + dy * cos,
        });
    }
}

// One placed atom: centre the polygon outward, away from the
// centroid of everything already placed (a spiro centre must not
// fold the new ring back onto the old one). The anchor sits opposite
// the centre direction, so the shift must put the anchor's own
// vertex there (+π).
fn place_anchored_cycle(members: &[usize], positions: &mut [Option<Point>]) {
    let anchor = *members
        .iter()
        .find(|member| positions[**member].is_some())
        .expect("anchored ring has a placed member");
    let anchor_index = members.iter().position(|member| *member == anchor).expect("member");
    let anchor_pos = positions[anchor].expect("placed");
    let vertex_angle = angle_of_member(members.len(), anchor_index);
    let center_angle = center_angle_away_from_placed(anchor_pos, positions, vertex_angle);
    let shift = center_angle + std::f64::consts::PI - vertex_angle;
    let radius = polygon_radius(members.len());
    let center = Point {
        x: anchor_pos.x + radius * center_angle.cos(),
        y: anchor_pos.y + radius * center_angle.sin(),
    };
    for (i, member) in members.iter().enumerate() {
        if positions[*member].is_some() {
            continue;
        }
        let ideal = ideal_point(members.len(), i, shift);
        positions[*member] = Some(Point {
            x: center.x + ideal.x,
            y: center.y + ideal.y,
        });
    }
}

fn center_angle_away_from_placed(
    anchor_pos: Point,
    positions: &[Option<Point>],
    fallback: f64,
) -> f64 {
    let placed: Vec<Point> = positions.iter().filter_map(|p| *p).collect();
    if placed.len() < 2 {
        return fallback;
    }
    let count = placed.len() as f64;
    let centroid = Point {
        x: placed.iter().map(|p| p.x).sum::<f64>() / count,
        y: placed.iter().map(|p| p.y).sum::<f64>() / count,
    };
    (anchor_pos.y - centroid.y).atan2(anchor_pos.x - centroid.x)
}

// Zigzag growth: each generation turns the opposite way from its
// parent, so chains alternate instead of curling. Tracked per atom
// as [parent, incomingAngle, turn].
fn grow(
    start: usize,
    positions: &mut [Option<Point>],
    graph: &LayoutGraph,
) {
    let start_pos = positions[start].expect("seed placed");
    let mut visited: std::collections::HashMap<usize, (Option<usize>, f64, f64)> =
        std::collections::HashMap::new();
    visited.insert(
        start,
        (None, start_pos.y.atan2(start_pos.x), -1.0),
    );
    let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    queue.push_back(start);
    while let Some(current) = queue.pop_front() {
        let Some(current_pos) = positions[current] else {
            continue;
        };
        let incoming_angle = visited.get(&current).expect("visited").1;
        let children: Vec<usize> = neighbors_of(graph, current)
            .into_iter()
            .filter(|n| !visited.contains_key(n) && positions[*n].is_none())
            .collect();
        for child in children {
            if visited.contains_key(&child) {
                continue;
            }
            let turn = -visited.get(&current).expect("visited").2;
            let angle = incoming_angle + turn * SPREAD;
            positions[child] = Some(clear_position(current_pos, angle, positions));
            visited.insert(child, (Some(current), angle, turn));
            queue.push_back(child);
        }
    }
}

// Try the preferred spot, then spread out; if the preferred spot is
// occupied (another substituent, a ring atom), fall back to the
// candidate farthest from any placed atom. Deterministic: first best
// wins ties.
fn clear_position(from_pos: Point, angle: f64, positions: &[Option<Point>]) -> Point {
    let spreads = [
        0.0,
        -SPREAD,
        SPREAD,
        -2.0 * SPREAD,
        2.0 * SPREAD,
        std::f64::consts::PI,
    ];
    let candidates: Vec<Point> = spreads
        .iter()
        .map(|spread| {
            let candidate_angle = angle + spread;
            Point {
                x: from_pos.x + BOND_LENGTH * candidate_angle.cos(),
                y: from_pos.y + BOND_LENGTH * candidate_angle.sin(),
            }
        })
        .collect();
    let occupied: Vec<Point> = positions.iter().filter_map(|p| *p).collect();
    for candidate in &candidates {
        if occupied
            .iter()
            .all(|other| distance(*candidate, *other) >= COLLISION_RADIUS)
        {
            return *candidate;
        }
    }
    let mut best = candidates[0];
    let mut best_score = -1.0_f64;
    for candidate in &candidates {
        let score = occupied
            .iter()
            .map(|other| distance(*candidate, *other))
            .sum::<f64>();
        if score > best_score {
            best_score = score;
            best = *candidate;
        }
    }
    best
}

fn distance(a: Point, b: Point) -> f64 {
    ((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)).sqrt()
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

// The ideal polygon vertex position relative to the centre.
fn ideal_point(sides: usize, member_index: usize, shift: f64) -> Point {
    let radius = polygon_radius(sides);
    let angle = angle_of_member(sides, member_index) + shift;
    Point {
        x: radius * angle.cos(),
        y: radius * angle.sin(),
    }
}

pub fn layout_molecule(molecule: &model::Molecule) -> (LayoutGraph, Vec<Point>) {
    let graph = walk_molecule(molecule);
    // Authored coordinates (the CML x2/y2 round-trip, molfile-carried)
    // take precedence over every computed layout — the depiction-hint
    // override: the renderer speaks the molecule's own geometry when it
    // carries one.
    if !graph.atoms.is_empty()
        && graph
            .atoms
            .iter()
            .all(|atom| atom.x2.is_some() && atom.y2.is_some())
    {
        let positions = graph
            .atoms
            .iter()
            .map(|atom| Point {
                x: atom.x2.expect("all placed"),
                y: atom.y2.expect("all placed"),
            })
            .collect();
        return (graph, positions);
    }
    let positions = place_graph(&graph);
    (graph, positions)
}
