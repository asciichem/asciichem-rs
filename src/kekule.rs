//! The text projection of an aromatic ring is its Kekulé form — the
//! Rust mirror of `AsciiChem::Structure::Kekule` (asciichem-ruby#188).
//!
//! Aromatic text bonds are unparseable (the Lewis lone-pair prefix
//! owns `:` in every position an aromatic bond could occupy), so
//! `c1ccccc1` renders as `C1-C=C-C=C-C1`. A Kekulé assignment is a
//! perfect matching on the component's atoms (each atom exactly one
//! double bond), found by deterministic backtracking in Graph edge
//! order (include before skip) so all runtimes derive identical
//! output. Components with no perfect matching (odd rings like furan)
//! stay out of the projection and keep the aromatic rendering (the
//! spec'd text loss). The model itself is never touched — the
//! projection is a render-time view, exactly like perception.

use std::collections::HashMap;

use crate::model::Molecule;
use crate::structure::{build_graph, Edge};

/// The Kekulé projection of `molecule`'s aromatic chain bonds:
/// `((from, to) -> "double" | "single")` keyed by ordered atom-index
/// pair (the same key the positional text walk derives for the bond
/// token it renders). Ring-closure edges carry no bond token (digits
/// are single by grammar) and get no entry.
pub fn bond_kinds(molecule: &Molecule) -> HashMap<(usize, usize), String> {
    let mut kinds = HashMap::new();
    let (_atoms, edges) = build_graph(molecule);
    let aromatic: Vec<&Edge> = edges.iter().filter(|edge| edge.kind == "aromatic").collect();
    if aromatic.is_empty() {
        return kinds;
    }

    for component in components(&aromatic) {
        if let Some(indices) = matching(&component) {
            for (index, edge) in component.iter().enumerate() {
                let (from, to) = edge.ordered();
                if to - from != 1 {
                    continue;
                }
                let kind = if indices.contains(&index) { "double" } else { "single" };
                kinds.insert((from, to), kind.to_string());
            }
        }
    }
    kinds
}

/// Connected components of the aromatic subgraph (edges sharing
/// atoms). Edges leave in shift order — the Ruby port's order.
fn components<'a>(aromatic: &[&'a Edge]) -> Vec<Vec<&'a Edge>> {
    let mut found: Vec<Vec<&Edge>> = Vec::new();
    let mut remaining: Vec<&Edge> = aromatic.to_vec();
    while !remaining.is_empty() {
        let seed = remaining.remove(0);
        let mut component = vec![seed];
        let mut atoms = vec![seed.from, seed.to];
        let mut grew = true;
        while grew {
            grew = false;
            let mut kept = Vec::new();
            for edge in remaining.drain(..) {
                if atoms.contains(&edge.from) || atoms.contains(&edge.to) {
                    component.push(edge);
                    atoms.push(edge.from);
                    atoms.push(edge.to);
                    grew = true;
                } else {
                    kept.push(edge);
                }
            }
            remaining = kept;
        }
        found.push(component);
    }
    found
}

/// Indices of a perfect matching over the component's atoms (each
/// atom incident to exactly one chosen edge), or nil. Closure edges
/// are never chosen: ring-closure digits are single by grammar, so a
/// matching that doubles one would silently lose that bond's order in
/// the text. Include before skip — the same deterministic order as
/// the Ruby port.
fn matching(component: &[&Edge]) -> Option<Vec<usize>> {
    let mut atom_set: Vec<usize> = component.iter().flat_map(|edge| [edge.from, edge.to]).collect();
    atom_set.sort_unstable();
    atom_set.dedup();
    let mut counts: HashMap<usize, i32> = HashMap::new();
    let mut chosen: Vec<usize> = Vec::new();
    let mut found: Option<Vec<usize>> = None;
    search(component, 0, &mut chosen, &mut counts, &atom_set, &mut found);
    found
}

#[allow(clippy::too_many_arguments)]
fn search(
    component: &[&Edge],
    index: usize,
    chosen: &mut Vec<usize>,
    counts: &mut HashMap<usize, i32>,
    atom_set: &[usize],
    found: &mut Option<Vec<usize>>,
) {
    if index == component.len() {
        if atom_set.iter().all(|atom| counts.get(atom).copied().unwrap_or(0) == 1) {
            *found = Some(chosen.clone());
        }
        return;
    }
    let edge = component[index];
    let chain = edge.from.abs_diff(edge.to) == 1;
    let from_count = counts.get(&edge.from).copied().unwrap_or(0);
    let to_count = counts.get(&edge.to).copied().unwrap_or(0);
    if chain && from_count != 1 && to_count != 1 {
        *counts.entry(edge.from).or_insert(0) += 1;
        *counts.entry(edge.to).or_insert(0) += 1;
        chosen.push(index);
        search(component, index + 1, chosen, counts, atom_set, found);
        chosen.pop();
        *counts.entry(edge.from).or_insert(0) -= 1;
        *counts.entry(edge.to).or_insert(0) -= 1;
    }
    if found.is_none() {
        search(component, index + 1, chosen, counts, atom_set, found);
    }
}
