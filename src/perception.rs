//! Perception pass (TODO.v2/13 item 2) — the Rust mirror of
//! AsciiChem::Perception. Chemistry the syntax leaves implicit,
//! derived once and consumed everywhere; semantics are conformance-
//! tested through the shared corpus (fixtures/perception.json), never
//! re-derived. `None` in the corpus means perception declines — this
//! implementation must decline identically.

use std::collections::HashSet;

use crate::model::Atom;
use crate::structure::Edge;

/// Default valences per element (the common ground-state values;
/// multi-valent elements list their allowed valences in order).
fn default_valences(element: &str) -> Option<&'static [f64]> {
    match element {
        "H" => Some(&[1.0]),
        "B" => Some(&[3.0]),
        "C" => Some(&[4.0]),
        "N" => Some(&[3.0]),
        "O" => Some(&[2.0]),
        "F" => Some(&[1.0]),
        "Al" => Some(&[3.0]),
        "Si" => Some(&[4.0]),
        "P" => Some(&[3.0, 5.0]),
        "S" => Some(&[2.0, 4.0, 6.0]),
        "Cl" => Some(&[1.0]),
        "Br" => Some(&[1.0]),
        "I" => Some(&[1.0]),
        _ => None,
    }
}

/// Valence weight of a bond kind; unknown kinds weigh 1.
fn bond_order(kind: &str) -> f64 {
    match kind {
        "double" => 2.0,
        "triple" => 3.0,
        "quadruple" => 4.0,
        _ => 1.0,
    }
}

/// Hydrogen-counting order of an aromatic bond (Hückel): 1.5, so both
/// spellings of benzene sum to 3.0 at every carbon and perceive one
/// implicit hydrogen.
const AROMATIC_H_ORDER: f64 = 1.5;

/// Neighbour indexes of an atom, in edge order.
pub fn neighbours(atom_index: usize, edges: &[Edge]) -> Vec<usize> {
    let mut acc = Vec::new();
    for edge in edges {
        if edge.from == atom_index {
            acc.push(edge.to);
        }
        if edge.to == atom_index {
            acc.push(edge.from);
        }
    }
    acc
}

/// Hydrogen-counting sum: aromatic bonds order as 1.5 (Hückel).
pub fn hydrogen_valence_sum(atom_index: usize, edges: &[Edge]) -> f64 {
    let mut sum = 0.0;
    for edge in edges {
        if edge.from != atom_index && edge.to != atom_index {
            continue;
        }
        sum += if edge.kind == "aromatic" {
            AROMATIC_H_ORDER
        } else {
            bond_order(&edge.kind)
        };
    }
    sum
}

/// Whole-number, non-negative gaps between each default valence and
/// the hydrogen-valence sum; empty when the element is unknown or no
/// valence fits integrally.
pub fn valence_deficits(element: &str, valence_sum: f64) -> Vec<f64> {
    let Some(valences) = default_valences(element) else {
        return Vec::new();
    };
    valences
        .iter()
        .map(|v| v - valence_sum)
        .filter(|d| *d >= 0.0 && d.fract() == 0.0)
        .collect()
}

/// Implicit hydrogens for an atom: the smallest fitting default
/// valence minus the hydrogen-valence sum.
/// - Atoms that STATE hydrogens (`[NH4+]`, `[CH4]`) return `Some(0)` —
///   never add to what the source said.
/// - Charged atoms return `None`: charge algebra is deferred.
/// - Unknown elements, non-integral or over-valent sums return `None`
///   — never invent.
pub fn implicit_hydrogens(atom: &Atom, valence_sum: f64) -> Option<u32> {
    if atom.hydrogens.is_some() {
        return Some(0);
    }
    if atom.charge.is_some() {
        return None;
    }

    let deficits = valence_deficits(&atom.element, valence_sum);
    let smallest = deficits.iter().copied().reduce(f64::min)?;
    Some(smallest as u32)
}

/// pi-electron-donating heteroatoms: a ring atom with no double bond
/// of its own contributes its lone pair (pyrrole N, furan O,
/// thiophene S) — the 2 electrons that make those five-membered rings
/// 4n+2.
const LONE_PAIR_DONORS: [&str; 6] = ["N", "O", "S", "P", "Se", "As"];

/// Aromatic rings by Hückel counting: 5- and 6-membered cycles whose
/// pi contributions total 6 or 10 (4n+2 with n >= 1; n = 0 is a
/// charged three-membered ring, which perception declines). Per atom:
/// an in-system double bond contributes 1; a lone-pair donor with no
/// double contributes 2; a carbon with only single bonds or an
/// exocyclic double contributes 0. Charged atoms make the count
/// undecidable and the ring is declined; aromatic-flagged rings are
/// taken as stated.
///
/// This sees fused systems (naphthalene's shared bond is single, but
/// each atom's double bond lies in the ring system) and keeps quinones
/// out (their carbons' only doubles point at exocyclic oxygens).
pub fn aromatic_rings(atoms: &[Atom], edges: &[Edge]) -> Vec<Vec<usize>> {
    let cycles = ring_cycles(edges, atoms.len());
    let cyclic_atoms: HashSet<usize> = cycles.iter().flatten().copied().collect();
    cycles
        .into_iter()
        .filter(|cycle| aromatic_cycle(cycle, atoms, edges, &cyclic_atoms))
        .collect()
}

fn aromatic_cycle(
    cycle: &[usize],
    atoms: &[Atom],
    edges: &[Edge],
    cyclic_atoms: &HashSet<usize>,
) -> bool {
    if cycle.iter().any(|i| atoms[*i].aromatic == Some(true)) {
        return true;
    }
    if cycle.len() != 5 && cycle.len() != 6 {
        return false;
    }
    let mut total = 0.0;
    for index in cycle {
        let Some(electrons) = ring_pi_electrons(*index, atoms, edges, cyclic_atoms) else {
            return false;
        };
        total += electrons;
    }
    total == 6.0 || total == 10.0
}

fn ring_pi_electrons(
    atom_index: usize,
    atoms: &[Atom],
    edges: &[Edge],
    cyclic_atoms: &HashSet<usize>,
) -> Option<f64> {
    let atom = &atoms[atom_index];
    if atom.charge.is_some() {
        return None;
    }
    if in_system_double(atom_index, edges, cyclic_atoms) {
        return Some(1.0);
    }
    if LONE_PAIR_DONORS.contains(&atom.element.as_str()) {
        return Some(2.0);
    }
    Some(0.0)
}

/// Whether the atom's only doubles reach INTO the ring system: an
/// exocyclic double (quinone carbonyl) does not count, so quinones
/// stay out.
fn in_system_double(atom_index: usize, edges: &[Edge], cyclic_atoms: &HashSet<usize>) -> bool {
    edges.iter().any(|edge| {
        if edge.from != atom_index && edge.to != atom_index {
            return false;
        }
        if edge.kind != "double" {
            return false;
        }
        let partner = if edge.from == atom_index {
            edge.to
        } else {
            edge.from
        };
        cyclic_atoms.contains(&partner)
    })
}

/// Simple cycles up to size 6 via DFS (molecules are small),
/// deduplicated by min/max endpoint (first occurrence kept — the Ruby
/// reference).
fn ring_cycles(edges: &[Edge], atom_count: usize) -> Vec<Vec<usize>> {
    let mut found = Vec::new();
    for start in 0..atom_count {
        dfs_cycles(edges, start, start, vec![start], 6, &mut found);
    }
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for cycle in found {
        let min = cycle.iter().min().copied().expect("cycle is non-empty");
        let max = cycle.iter().max().copied().expect("cycle is non-empty");
        if seen.insert((min, max)) {
            out.push(cycle);
        }
    }
    out
}

fn dfs_cycles(
    edges: &[Edge],
    start: usize,
    current: usize,
    path: Vec<usize>,
    max_size: usize,
    found: &mut Vec<Vec<usize>>,
) {
    // Step straight back to where we came from never closes a ring.
    let came_from = path.len().checked_sub(2).map(|i| path[i]);
    for nb in neighbours(current, edges) {
        if Some(nb) == came_from {
            continue;
        }
        if nb == start {
            if path.len() >= 3 {
                found.push(path.clone());
            }
        } else if !path.contains(&nb) && path.len() < max_size {
            let mut next = path.clone();
            next.push(nb);
            dfs_cycles(edges, start, nb, next, max_size, found);
        }
    }
}
