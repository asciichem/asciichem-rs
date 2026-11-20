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
use crate::model::{Atom, Bond, Group, Molecule, Node};
use crate::ParseError;

pub const MAX_DIGIT: usize = 9;

#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub kind: String,
    pub direction: Option<String>,
    /// The SMILES reader's structural mark — the linearizer
    /// partitions by it instead of index arithmetic (branch atoms
    /// parsed before their continuation are index-consecutive and
    /// would masquerade as chain tokens). None for edges derived from
    /// an already-linearized model.
    pub shape: Option<&'static str>,
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
                        shape: None,
                    }),
                    (None, Some(from)) => adjacency_pairs.push((from, index)),
                    _ => {}
                }
                *last = Some(index);
            }
            Node::Bond(bond) => *pending = Some(bond.clone()),
            Node::Group(group) => {
                // A group whose FIRST node is a bond is a BRANCH group:
                // the leading bond attaches the branch to the atom
                // before the group (the branch root); the continuation
                // resumes from the root — `C(-O)=O` reads C-O and C=O,
                // never O=O.
                let mut first = group.nodes.first();
                while let Some(Node::Molecule(mol)) = first {
                    first = mol.nodes.first();
                }
                let branch_root = if matches!(first, Some(Node::Bond(_))) {
                    *last
                } else {
                    None
                };
                for child in &group.nodes {
                    walk(child, atoms, edges, adjacency_pairs, pending, last);
                }
                if let Some(root) = branch_root {
                    *last = Some(root);
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
            shape: None,
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
    // Partition by the parser's own shape marks — tokens and digits
    // first (so the chain map is complete), then the branch trees. A
    // branch edge inside another branch's tree emits as a NESTED
    // group; only top-level branches render at the flat stream.
    // Producers without shape marks (molfile, CML) ride the flat rule:
    // consecutive indices are chain tokens, every other non-ring edge
    // a digit pair — the shape those producers were recorded against.
    let mut close_at = [-1i64; MAX_DIGIT];
    let mut digits_for: Vec<String> = vec![String::new(); atoms.len()];
    let mut token_before: Vec<Option<&Edge>> = vec![None; atoms.len()];
    let mut children: Vec<Vec<&Edge>> = vec![Vec::new(); atoms.len()];
    let mut in_branch: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut branch_edges: Vec<&Edge> = Vec::new();

    for edge in edges {
        if edge.shape == Some("branch") {
            branch_edges.push(edge);
        } else if edge.shape == Some("chain")
            || (edge.shape.is_none() && edge.ordered().1 - edge.ordered().0 == 1)
        {
            token_before[edge.to] = Some(edge);
            children[edge.from].push(edge);
        } else {
            let (first, last) = edge.ordered();
            let digit_index = close_at
                .iter()
                .position(|&closed| closed <= first as i64)
                .ok_or_else(|| ParseError {
                    message: format!(
                        "more than {MAX_DIGIT} overlapping non-adjacent bonds — beyond the model's ring-closure digit capacity"
                    ),
                })?;
            close_at[digit_index] = last as i64;
            let digit = (digit_index + 1).to_string();
            digits_for[first].push_str(&digit);
            digits_for[last].push_str(&digit);
        }
    }
    for edge in &branch_edges {
        collect_branch(edge.to, &children, &mut in_branch);
        children[edge.from].push(edge);
    }
    let mut branch_at: std::collections::BTreeMap<usize, Vec<&Edge>> =
        std::collections::BTreeMap::new();
    for edge in &branch_edges {
        if !in_branch.contains(&edge.from) {
            branch_at.entry(edge.from).or_default().push(edge);
        }
    }

    fn collect_branch(
        start: usize,
        children: &[Vec<&Edge>],
        in_branch: &mut std::collections::BTreeSet<usize>,
    ) {
        in_branch.insert(start);
        for edge in &children[start] {
            collect_branch(edge.to, children, in_branch);
        }
    }

    fn branch_group(
        edge: &Edge,
        atoms: &[Atom],
        children: &[Vec<&Edge>],
        branch_edges: &[&Edge],
    ) -> Node {
        let mut nodes: Vec<Node> = vec![Node::Bond(Bond {
            kind: edge.kind.clone(),
            direction: edge.direction.clone(),
        })];
        emit_tree(edge.to, &mut nodes, atoms, children, branch_edges);
        Node::Group(Group {
            nodes,
            multiplicity: None,
            bracket: "paren".to_string(),
        })
    }

    fn emit_tree(
        atom_index: usize,
        nodes: &mut Vec<Node>,
        atoms: &[Atom],
        children: &[Vec<&Edge>],
        branch_edges: &[&Edge],
    ) {
        nodes.push(Node::Atom(atoms[atom_index].clone()));
        // Branch groups sit directly after their parent atom (the
        // walkers bind a group's leading bond to the atom before it —
        // emitting a continuation first would steal the binding).
        for edge in &children[atom_index] {
            if branch_edges.contains(edge) {
                nodes.push(branch_group(edge, atoms, children, branch_edges));
            }
        }
        for edge in &children[atom_index] {
            if branch_edges.contains(edge) {
                continue;
            }
            nodes.push(Node::Bond(Bond {
                kind: edge.kind.clone(),
                direction: edge.direction.clone(),
            }));
            emit_tree(edge.to, nodes, atoms, children, branch_edges);
        }
    }

    let mut final_atoms = atoms;
    for index in 0..final_atoms.len() {
        let digits = std::mem::take(&mut digits_for[index]);
        if !digits.is_empty() {
            final_atoms[index].ring_closures = Some(digits);
        }
    }
    let mut nodes: Vec<Node> = Vec::with_capacity(final_atoms.len() * 2);
    for (index, atom) in final_atoms.iter().enumerate() {
        if in_branch.contains(&index) {
            continue;
        }
        if index > 0 {
            if let Some(edge) = token_before[index] {
                nodes.push(Node::Bond(Bond {
                    kind: edge.kind.clone(),
                    direction: edge.direction.clone(),
                }));
            }
        }
        nodes.push(Node::Atom(atom.clone()));
        for edge in branch_at.get(&index).into_iter().flatten() {
            nodes.push(branch_group(edge, &final_atoms, &children, &branch_edges));
        }
    }
    Ok(nodes)
}

#[cfg(test)]
mod branch_projection_tests {
    use crate::model::Molecule;
    use crate::smiles::parse_smiles;

    #[test]
    fn smiles_molecules_text_as_branch_groups() {
        // The branch-projection reference (asciichem-ruby#208): the
        // reader shapes the linearizer, so smiles-ingested molecules
        // text as branch groups instead of carbonyl ring digits —
        // byte-identical with the ruby reference.
        for (source, expected) in [
            ("CC(=O)O", "C-C(=O)-O"),
            (
                "CC(=O)OC1=CC=CC=C1C(=O)O",
                "C-C(=O)-O-C1=C-C=C-C=C1-C(=O)-O",
            ),
            ("N[C@@H](C)C(=O)O", "N-C@TH2(-C)-C(=O)-O"),
            ("C1=CC=CC=C1", "C1=C-C=C-C=C1"),
        ] {
            let nodes = parse_smiles(source).expect("parses");
            let molecule = Molecule {
                nodes,
                ..Molecule::default()
            };
            let rendered = crate::model::Node::Molecule(molecule).to_string();
            assert_eq!(rendered, expected, "{source}");
        }
    }
}
