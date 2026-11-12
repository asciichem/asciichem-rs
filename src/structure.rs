//! Shared structure plumbing for the interchange formats — the Rust
//! mirror of AsciiChem::Structure (TODO.v2 09).
//!
//! - [`build_graph`] turns a model molecule into a neutral atom array
//!   plus bond list (adjacency) using the same pending-bond walk as
//!   the layout walker, plus ring-closure edges (aromatic between
//!   aromatic atoms — the same default rule the SMILES parser
//!   applies).
//! - [`build_graph_strict`] is the structure-writers' build (the Ruby
//!   `build!`): it declines formula-style adjacency instead of
//!   silently fragmenting the graph.
//! - [`linearize`] turns the adjacency back into linear model nodes:
//!   bond tokens for consecutive edges, ring-closure digits for the
//!   rest, with interval-aware digit reuse (digits pair by occurrence
//!   order, so overlapping intervals never share a digit).

use crate::layout;
use crate::model::{Atom, Bond, Molecule, Node};
use crate::ParseError;

pub const MAX_DIGIT: usize = 9;

#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub kind: String,
    pub direction: Option<String>,
}

impl Edge {
    pub fn ordered(&self) -> (usize, usize) {
        (self.from.min(self.to), self.from.max(self.to))
    }
}

/// Walks a molecule into atoms + edges. Bond tokens connect the last
/// atom of the previous unit to the first atom of the next (carrying
/// kind and E/Z direction); nested groups and molecules contribute
/// their walk order; ring closures pair by digit occurrence order.
pub fn build_graph(molecule: &Molecule) -> (Vec<Atom>, Vec<Edge>) {
    let (atoms, edges, _) = build_graph_with_adjacency(molecule);
    (atoms, edges)
}

/// The structure-writers' strict build — the Ruby `Structure.Graph.build!`.
/// Text adjacency — atoms with no bond token between them (`CCO`, the
/// undashed spacers in `C1=CC=CC=C1`) — is formula notation, not a
/// bond; mixing it with explicit bonds silently fragments the graph
/// (a benzene ring written as a disconnected C4+C2 pair), so the
/// strict build declines when adjacency is present and the final graph
/// is disconnected. SMILES models never trip this (their parser always
/// bonds every continuation — branch reconnection crosses as a
/// ring-closure digit), and an explicit multi-component model (a salt
/// written `A.B`) has no adjacency at all.
pub fn build_graph_strict(molecule: &Molecule) -> Result<(Vec<Atom>, Vec<Edge>), ParseError> {
    let (atoms, edges, adjacency_pairs) = build_graph_with_adjacency(molecule);
    if !adjacency_pairs.is_empty() && !connected(atoms.len(), &edges) {
        return Err(ParseError {
            message: "adjacent atoms have no explicit bond token — a formula \
                      is not a structure. Write explicit bonds (C-C) or ingest \
                      the structure from SMILES/molfile."
                .into(),
        });
    }
    Ok((atoms, edges))
}

/// Breadth-first reachability from atom 0 — the Ruby `connected?`.
fn connected(atom_count: usize, edges: &[Edge]) -> bool {
    if atom_count <= 1 {
        return true;
    }
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); atom_count];
    for edge in edges {
        neighbours[edge.from].push(edge.to);
        neighbours[edge.to].push(edge.from);
    }
    let mut seen = vec![false; atom_count];
    seen[0] = true;
    let mut queue = std::collections::VecDeque::from([0usize]);
    while let Some(current) = queue.pop_front() {
        for &next in &neighbours[current] {
            if seen[next] {
                continue;
            }
            seen[next] = true;
            queue.push_back(next);
        }
    }
    seen.into_iter().all(|visited| visited)
}

fn build_graph_with_adjacency(molecule: &Molecule) -> (Vec<Atom>, Vec<Edge>, Vec<(usize, usize)>) {
    let mut atoms: Vec<Atom> = Vec::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut adjacency_pairs: Vec<(usize, usize)> = Vec::new();
    let mut pending: Option<Bond> = None;
    let mut last: Option<usize> = None;

    fn walk(
        node: &Node,
        atoms: &mut Vec<Atom>,
        edges: &mut Vec<Edge>,
        adjacency_pairs: &mut Vec<(usize, usize)>,
        pending: &mut Option<Bond>,
        last: &mut Option<usize>,
    ) {
        match node {
            Node::Atom(atom) => {
                let index = atoms.len();
                atoms.push(atom.clone());
                match (pending.take(), *last) {
                    (Some(bond), Some(from)) => edges.push(Edge {
                        from,
                        to: index,
                        kind: bond.kind.clone(),
                        direction: bond.direction.clone(),
                    }),
                    (None, Some(from)) => adjacency_pairs.push((from, index)),
                    _ => {}
                }
                *last = Some(index);
            }
            Node::Bond(bond) => *pending = Some(bond.clone()),
            Node::Group(group) => {
                for child in &group.nodes {
                    walk(child, atoms, edges, adjacency_pairs, pending, last);
                }
            }
            Node::Molecule(mol) => {
                for child in &mol.nodes {
                    walk(child, atoms, edges, adjacency_pairs, pending, last);
                }
            }
            _ => {}
        }
    }

    for node in &molecule.nodes {
        walk(
            node,
            &mut atoms,
            &mut edges,
            &mut adjacency_pairs,
            &mut pending,
            &mut last,
        );
    }

    // Ring-closure pairing by digit occurrence order (the same
    // RingBonds walk the layout uses): aromatic between aromatic
    // atoms, single otherwise.
    for ring in layout::ring_bonds(&atoms) {
        let (from, to) = (ring.from_atom, ring.to_atom);
        let kind = if atoms[from].aromatic == Some(true) && atoms[to].aromatic == Some(true) {
            "aromatic"
        } else {
            "single"
        };
        edges.push(Edge {
            from,
            to,
            kind: kind.to_string(),
            direction: None,
        });
    }

    (atoms, edges, adjacency_pairs)
}

/// Adjacency → linear model nodes: a bond token before each atom that
/// continues the chain, ring-closure digits on both endpoints of every
/// non-consecutive edge. Digits are reused interval-aware: a digit is
/// only handed out again once the previous bond carrying it has closed
/// at or before the new bond's first endpoint.
pub fn linearize(atoms: Vec<Atom>, edges: &[Edge]) -> Result<Vec<Node>, ParseError> {
    let mut close_at = [-1i64; MAX_DIGIT];
    let mut digits_for: Vec<String> = vec![String::new(); atoms.len()];
    let mut token_before: Vec<Option<&Edge>> = vec![None; atoms.len()];

    for edge in edges {
        let (first, last) = edge.ordered();
        if last - first == 1 {
            token_before[last] = Some(edge);
            continue;
        }
        let digit_index = close_at
            .iter()
            .position(|&closed| closed <= first as i64)
            .ok_or_else(|| {
                ParseError {
                    message: format!(
                        "more than {MAX_DIGIT} overlapping non-adjacent bonds — beyond the model's ring-closure digit capacity"
                    ),
                }
            })?;
        close_at[digit_index] = last as i64;
        let digit = (digit_index + 1).to_string();
        digits_for[first].push_str(&digit);
        digits_for[last].push_str(&digit);
    }

    let mut owned_atoms = atoms;
    let mut nodes: Vec<Node> = Vec::with_capacity(owned_atoms.len() * 2);
    for (index, mut atom) in owned_atoms.drain(..).enumerate() {
        if index > 0 {
            if let Some(edge) = token_before[index] {
                nodes.push(Node::Bond(Bond {
                    kind: edge.kind.clone(),
                    direction: edge.direction.clone(),
                }));
            }
        }
        let digits = std::mem::take(&mut digits_for[index]);
        if !digits.is_empty() {
            atom.ring_closures = Some(digits);
        }
        nodes.push(Node::Atom(atom));
    }
    Ok(nodes)
}
