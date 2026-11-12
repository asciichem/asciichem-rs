//! The identity layer (TODO.v2/17 phase 1, slice 3): standard InChI,
//! AuxInfo and InChIKey derived IN-PROCESS through the
//! [chemomer](https://github.com/chemomer/chemomer) crate — the same
//! native engine the Ruby gem reaches through its C-ABI seam
//! (`NativeEngine`) and the TS/PY runtimes through wasm/PyO3.
//! Same-language means no seam at all: the dependency is version-pinned
//! on crates.io and byte-exact against the shared oracle corpus
//! (chemomer/inchi-tests; this repo's CI gates the identity family at
//! the same pin chemomer's own CI gates).
//!
//! Semantics mirror the Ruby `NativeEngine` contract: a structure goes
//! in as a V2000 molfile, the identity comes out as the standard
//! InChI string, its AuxInfo line, and the InChIKey. The model →
//! molfile emission is the pipe; the writer here is PRIVATE to the
//! identity path and minimal (V2000, no coordinates — InChI never
//! reads them). The full molfile interchange (read/write as a public
//! format) is TODO.v2/09's slice 2, which replaces this writer.

use crate::model::{self, Node};

/// The derived identity of a structure: the standard InChI string,
/// its AuxInfo line, and the InChIKey.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Identity {
    pub inchi: String,
    pub aux: String,
    pub inchikey: String,
}

/// Derives the identity of a V2000 molfile — the engine seam itself.
/// Every runtime's identity seam is this shape: molfile in, the three
/// identity strings out.
pub fn derive_from_molfile(molfile: &str) -> chemomer::Result<Identity> {
    let engine = chemomer::Engine::from_molfile(molfile)?;
    let inchi = engine.inchi_string()?;
    let aux = engine.aux_string()?;
    let inchikey = engine.inchikey()?;
    Ok(Identity {
        inchi,
        aux,
        inchikey,
    })
}

/// Derives the identity of a model molecule: model → V2000 molfile →
/// engine. The molfile emission is the pipe, exactly as the Ruby
/// `NativeEngine#identity` (`Molfile.write(molecule, layout: false)`).
pub fn derive(molecule: &model::Molecule) -> chemomer::Result<Identity> {
    derive_from_molfile(&write_v2000(molecule)?)
}

/// One modelled bond in the built graph, with array indices into the
/// atom list.
struct Edge {
    from: usize,
    to: usize,
    kind: String,
}

/// Model molecule → neutral graph: bond tokens between consecutive
/// atoms become edges (pending-bond semantics, shared across group and
/// nested-molecule recursion — the same walk the Ruby
/// `Structure::Graph` runs); ring-closure digits become edges through
/// the ported `RingBonds` pairing.
///
/// Strict on purpose, exactly as the Ruby `build!`: text adjacency —
/// atoms with no bond token between them (`CCO`) — is formula
/// notation, not a bond, and mixing it with explicit bonds silently
/// fragments the graph, so the build declines when adjacency is
/// present and the final graph is disconnected. The identity of a
/// fragmented molecule is silently WRONG, which is worse than a
/// rejection.
struct GraphBuilder {
    atoms: Vec<model::Atom>,
    edges: Vec<Edge>,
    adjacency_pairs: Vec<(usize, usize)>,
    pending: Option<model::Bond>,
    last: Option<usize>,
}

impl GraphBuilder {
    fn new() -> Self {
        GraphBuilder {
            atoms: Vec::new(),
            edges: Vec::new(),
            adjacency_pairs: Vec::new(),
            pending: None,
            last: None,
        }
    }

    fn walk_nodes(&mut self, nodes: &[Node]) {
        for node in nodes {
            match node {
                Node::Atom(atom) => self.push_atom(atom),
                Node::Bond(bond) => self.pending = Some(bond.clone()),
                Node::Group(group) => self.walk_nodes(&group.nodes),
                Node::Molecule(inner) => self.walk_nodes(&inner.nodes),
                _ => {}
            }
        }
    }

    fn push_atom(&mut self, atom: &model::Atom) {
        let index = self.atoms.len();
        self.atoms.push(atom.clone());
        match (self.last, self.pending.take()) {
            (Some(last), Some(bond)) => self.edges.push(Edge {
                from: last,
                to: index,
                kind: bond.kind.clone(),
            }),
            (Some(last), None) => self.adjacency_pairs.push((last, index)),
            _ => {}
        }
        self.last = Some(index);
    }
}

fn build_graph(molecule: &model::Molecule) -> chemomer::Result<(Vec<model::Atom>, Vec<Edge>)> {
    let mut builder = GraphBuilder::new();
    builder.walk_nodes(&molecule.nodes);
    for ring in crate::layout::ring_bonds(&builder.atoms) {
        builder.edges.push(Edge {
            from: ring.from_atom,
            to: ring.to_atom,
            // Ring-closure digits carry no kind; the default bond rule
            // applies (the Ruby rule speaks aromatic only between
            // aromatic-FLAGGED atoms — a field the SMILES parser sets
            // and the ASCII text never does; the Rust model subset
            // carries no aromatic flag yet, so the closure crosses as
            // the default single — the same verdict the Ruby writer
            // reaches for ASCII input).
            kind: "single".to_string(),
        });
    }
    if !builder.adjacency_pairs.is_empty() && !connected(builder.atoms.len(), &builder.edges) {
        return Err(chemomer::Error::new(
            "adjacent atoms have no explicit bond token — a formula is not a structure. \
             Write explicit bonds (C-C) or ingest the structure from SMILES/molfile.",
        ));
    }
    Ok((builder.atoms, builder.edges))
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

/// The model's number-then-sign charge ("2+", "-") → signed integer —
/// the Ruby `charge_value`.
fn charge_value(charge: &str) -> i32 {
    let digits: String = charge.chars().take_while(|c| c.is_ascii_digit()).collect();
    let count = if digits.is_empty() {
        1
    } else {
        digits.parse().unwrap_or(1)
    };
    if charge.ends_with('-') {
        -count
    } else {
        count
    }
}

/// V2000 bond order code — the Ruby `BOND_TYPES`. Wedge/hash cross as
/// single order carrying the wedge through the stereo field.
fn bond_type(kind: &str) -> chemomer::Result<i32> {
    match kind {
        "single" | "wedge" | "hash" => Ok(1),
        "double" => Ok(2),
        "triple" => Ok(3),
        "aromatic" => Ok(4),
        other => Err(chemomer::Error::new(format!(
            "{other} bonds have no molfile V2000 type"
        ))),
    }
}

/// V2000 bond stereo code — the Ruby `BOND_STEREO`: wedge up (1),
/// hash down (6), everything else unbearing (0).
fn bond_stereo(kind: &str) -> i32 {
    match kind {
        "wedge" => 1,
        "hash" => 6,
        _ => 0,
    }
}

/// The atom-line stereo-parity field: TH1 → 1, TH2 → 2, otherwise 0 —
/// the Ruby `CHIRALITY_PARITY` fetch. The oracle ignores the field for
/// 2D input; enantiomers separate through the WEDGE, which
/// [`wedge_bonds`] lays on the stereocentre's first bond.
fn chirality_parity(atom: &model::Atom) -> i32 {
    match atom.chirality.as_deref() {
        Some("TH1") => 1,
        Some("TH2") => 2,
        _ => 0,
    }
}

/// One wedge per chiral atom, on its first bond, oriented
/// centre → neighbour: TH1 rises (stereo 1), TH2 falls (6). Returns
/// { edge index → (from, to, stereo) } — the Ruby `wedge_bonds`.
fn wedge_bonds(
    atoms: &[model::Atom],
    edges: &[Edge],
) -> std::collections::HashMap<usize, (usize, usize, i32)> {
    let mut wedges = std::collections::HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        let stereo = match atom.chirality.as_deref() {
            Some("TH1") => 1,
            Some("TH2") => 6,
            _ => continue,
        };
        let Some(edge_index) = edges
            .iter()
            .position(|edge| edge.from == index || edge.to == index)
        else {
            continue;
        };
        let edge = &edges[edge_index];
        wedges.insert(
            edge_index,
            if edge.from == index {
                (index, edge.to, stereo)
            } else {
                (index, edge.from, stereo)
            },
        );
    }
    wedges
}

/// One `M  CHG` / `M  ISO` property line; `None` (dropped) when the
/// pair list is empty — the Ruby `property_line`.
fn property_line(prefix: &str, pairs: &[(i32, i32)]) -> Option<String> {
    if pairs.is_empty() {
        return None;
    }
    let mut line = format!("{prefix}{:3}", pairs.len());
    for (index, value) in pairs {
        line.push_str(&format!("{:4}{:4}", index, value));
    }
    Some(line)
}

/// Model molecule → V2000 molfile — the minimal writer the identity
/// pipe needs (TODO.v2/09 slice 2 replaces it with the full
/// interchange). Field formats mirror the Ruby `Molfile::Writer`
/// byte-for-byte; coordinates are all zero (the `layout: false` call
/// the Ruby identity seam makes — InChI never reads coordinates, and
/// tetrahedral stereo crosses through the wedge bonds, not geometry).
fn write_v2000(molecule: &model::Molecule) -> chemomer::Result<String> {
    let (atoms, edges) = build_graph(molecule)?;
    if edges.is_empty() && atoms.len() > 1 {
        return Err(chemomer::Error::new(
            "molecule has no bonds — a formula is not a structure",
        ));
    }

    let mut lines: Vec<String> = vec![
        String::new(),
        "  AsciiChem".to_string(),
        String::new(),
        format!(
            "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000",
            atoms.len(),
            edges.len()
        ),
    ];

    let wedges = wedge_bonds(&atoms, &edges);
    for atom in &atoms {
        // The H-count field stays 0: the model subset carries no
        // explicit hydrogen counts yet (the Ruby writer crosses a
        // stated count as MDL code n + 1).
        lines.push(format!(
            "{:10.4}{:10.4}{:10.4} {:<3} 0  0{:3}{:3}  0  0  0  0  0  0  0  0  0",
            0.0,
            0.0,
            0.0,
            atom.element,
            chirality_parity(atom),
            0,
        ));
    }

    for (index, edge) in edges.iter().enumerate() {
        let stereo = bond_stereo(&edge.kind);
        let kind = bond_type(&edge.kind)?;
        let (from, to, stereo) = match wedges.get(&index) {
            Some(&(from, to, wedge_stereo)) => (from, to, wedge_stereo),
            None => (edge.from, edge.to, stereo),
        };
        lines.push(format!(
            "{:3}{:3}{:3}{:3}  0  0  0  0  0  0  0",
            from + 1,
            to + 1,
            kind,
            stereo
        ));
    }

    let charges: Vec<(i32, i32)> = atoms
        .iter()
        .enumerate()
        .filter_map(|(i, atom)| {
            atom.charge
                .as_ref()
                .map(|charge| (i as i32 + 1, charge_value(charge)))
        })
        .collect();
    if let Some(line) = property_line("M  CHG", &charges) {
        lines.push(line);
    }
    let isotopes: Vec<(i32, i32)> = atoms
        .iter()
        .enumerate()
        .filter_map(|(i, atom)| {
            atom.isotope
                .as_ref()
                .map(|isotope| (i as i32 + 1, isotope.parse().unwrap_or(0)))
        })
        .collect();
    if let Some(line) = property_line("M  ISO", &isotopes) {
        lines.push(line);
    }

    lines.push("M  END".to_string());
    Ok(lines.join("\n") + "\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn molecule(input: &str) -> model::Molecule {
        let nodes = crate::parse_nodes(input).expect("input parses");
        match nodes.first() {
            Some(Node::Molecule(mol)) => mol.clone(),
            other => panic!("not a molecule: {other:?}"),
        }
    }

    fn identity(input: &str) -> Identity {
        derive(&molecule(input)).expect("identity derives")
    }

    #[test]
    fn water_derives_the_recorded_identity() {
        let identity = identity("H-O-H");
        assert_eq!(identity.inchi, "InChI=1S/H2O/h1H2");
        assert_eq!(identity.inchikey, "XLYOFNOQVPJJNP-UHFFFAOYSA-N");
        assert!(identity.aux.starts_with("AuxInfo=1/0/N:"));
    }

    #[test]
    fn acetic_acid_walks_the_branches_and_the_double_bond() {
        // ASCII branch groups carry a leading bond token; the walk's
        // pending-bond semantics chain them linearly — EXACTLY as the
        // Ruby reference's Structure::Graph does (verified against the
        // gem: atoms C,H,H,H,C,O,O with the same edges). Branched
        // structures are SMILES/molfile-ingested in Ruby; ASCII branch
        // input fragments the same way on both sides — a known gap
        // slice 2's interchange owns. The gate pins the PARITY so a
        // silent divergence can never ship.
        let (atoms, edges) = build_graph(&molecule("C(-H)(-H)(-H)-C(=O)-O")).expect("builds");
        assert_eq!(
            atoms.iter().map(|a| a.element.as_str()).collect::<Vec<_>>(),
            ["C", "H", "H", "H", "C", "O", "O"]
        );
        let shape: Vec<(usize, &str, usize)> = edges
            .iter()
            .map(|edge| (edge.from, edge.kind.as_str(), edge.to))
            .collect();
        assert_eq!(
            shape,
            [
                (0, "single", 1),
                (1, "single", 2),
                (2, "single", 3),
                (3, "single", 4),
                (4, "double", 5),
                (5, "single", 6),
            ]
        );
    }

    #[test]
    fn benzene_pairs_the_ring_closure() {
        let identity = identity("C1=C-C=C-C=C1");
        assert_eq!(identity.inchi, "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H");
        assert_eq!(identity.inchikey, "UHOVQNZJYSORNB-UHFFFAOYSA-N");
    }

    #[test]
    fn ethanol_derives_the_linear_chain() {
        let identity = identity("C-C-O");
        assert_eq!(identity.inchi, "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
        assert_eq!(identity.inchikey, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
    }

    #[test]
    fn charges_cross_as_the_chg_block() {
        // The methoxide ion: the charge rides the property block, the
        // engine speaks the /q layer. (Recorded from the engine, which
        // is oracle-gated by the chemomer corpus.)
        let identity = identity("C-O^-");
        assert_eq!(identity.inchi, "InChI=1S/CH3O/c1-2/h1H3/q-1");
        assert_eq!(identity.inchikey, "NBTOZLQBSIZIKS-UHFFFAOYSA-N");
    }

    #[test]
    fn isotopes_cross_as_the_iso_block() {
        // ^13C-labelled methane: the isotope rides M  ISO; the
        // identity speaks the /i layer and the key carries the isotope
        // flag block.
        let identity = identity("^13C-H");
        assert_eq!(identity.inchi, "InChI=1S/CH4/h1H4/i1+1");
        assert_eq!(identity.inchikey, "VNWKTOKETHGBQD-OUBTZVSYSA-N");
    }

    #[test]
    fn wedge_stereo_crosses_through_the_bond_stereo_code() {
        // The TH1/TH2 wedge lands on the stereocentre's first bond
        // (up/down in the bond block's stereo field). Propane has no
        // stereocentre, so the identifiers agree — but the AuxInfo
        // echo records the wedge orientation the engine read (rB:P1
        // up vs rB:N1 down), proving the wedge crossed. Real stereo
        // separation gates through the corpus identity-stereo family,
        // whose recorded wedged molfiles carry the t/m layers.
        let th1 = identity("C@TH1-C-C");
        let th2 = identity("C@TH2-C-C");
        assert_eq!(th1.inchi, th2.inchi);
        assert!(th1.aux.contains("rB:P1"), "TH1 wedge up: {}", th1.aux);
        assert!(th2.aux.contains("rB:N1"), "TH2 wedge down: {}", th2.aux);
    }

    #[test]
    fn a_formula_is_not_a_structure() {
        // H_2O the formula: two atoms, no bond token — the writer
        // declines instead of deriving an identity for a fragmented
        // graph (the Ruby writer's contract).
        let error = derive(&molecule("H_2O")).expect_err("declined");
        assert!(error.message().contains("a formula is not a structure"));
    }

    #[test]
    fn formula_adjacency_with_bonds_is_declined() {
        // `CCO` with one explicit bond: adjacency fragments the graph
        // (C-C bonded pair + lone O), so the strict build declines.
        let error = derive(&molecule("CC-O")).expect_err("declined");
        assert!(error.message().contains("no explicit bond token"));
    }

    #[test]
    fn unknown_bond_kinds_have_no_v2000_type() {
        let error = derive(&molecule("C~>O")).expect_err("declined");
        assert!(error.message().contains("no molfile V2000 type"));
    }

    #[test]
    fn a_lone_atom_is_a_valid_structure() {
        // One atom, zero bonds: the engine's valence perception
        // coalesces a lone neutral oxygen to water (the oracle's
        // lone-atom rule — the same family as a lone H reading H2).
        let identity = identity("O");
        assert_eq!(identity.inchi, "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn the_written_molfile_carries_the_v2000_shape() {
        let molfile = write_v2000(&molecule("C-O^-")).expect("writes");
        let lines: Vec<&str> = molfile.lines().collect();
        assert_eq!(lines[0], "");
        assert_eq!(lines[1], "  AsciiChem");
        assert_eq!(lines[3], "  2  1  0  0  0  0  0  0  0  0999 V2000");
        assert!(
            lines[4].contains(" C   0  0  0  0"),
            "atom line: {}",
            lines[4]
        );
        // The bond line: order 1, unbearing stereo.
        assert!(lines[6].starts_with("  1  2  1  0"));
        // The charge block: one pair, atom 2, -1.
        assert_eq!(lines[7], "M  CHG  1   2  -1");
        assert_eq!(lines[8], "M  END");
        assert!(molfile.ends_with('\n'));
    }

    #[test]
    fn the_wedge_lands_on_the_first_bond_of_the_stereocentre() {
        // TH1 on the first atom, first bond: the bond line's stereo
        // field reads 1 (up), the atom line's parity field 1.
        let molfile = write_v2000(&molecule("C@TH1-C-C")).expect("writes");
        let lines: Vec<&str> = molfile.lines().collect();
        assert!(
            lines[4].contains(" C   0  0  1  0"),
            "atom line: {}",
            lines[4]
        );
        assert!(
            lines[7].starts_with("  1  2  1  1"),
            "bond line: {}",
            lines[7]
        );
        // TH2: parity 2 on the atom line, stereo 6 (down) on the bond.
        let molfile = write_v2000(&molecule("C@TH2-C-C")).expect("writes");
        let lines: Vec<&str> = molfile.lines().collect();
        assert!(
            lines[4].contains(" C   0  0  2  0"),
            "atom line: {}",
            lines[4]
        );
        assert!(
            lines[7].starts_with("  1  2  1  6"),
            "bond line: {}",
            lines[7]
        );
    }
}
