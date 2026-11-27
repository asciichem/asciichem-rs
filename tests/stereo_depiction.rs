//! The enantiomer-pair oracle gate (asciichem-rs#82): every synthesized
//! wedge drawing, re-read as a V2000 molfile with real drawn
//! coordinates through the chemomer engine, must derive the same
//! identity as its source SMILES. The drawing is the contract — if the
//! synthesis polarity were flipped, the engine would derive the OTHER
//! enantiomer and this gate goes red (the same catch that fitted
//! chemomer's phase-2 polarity: mark 1 → dash, mark 2 → up). Runs
//! unconditionally — the engine is a built-in dependency here.

use asciichem::layout;
use asciichem::model::{self, Node};
use asciichem::svg;

fn molecule_from_smiles(smiles: &str) -> model::Molecule {
    let nodes = asciichem::smiles::parse_smiles(smiles).expect("SMILES parses");
    match nodes.first() {
        Some(Node::Molecule(mol)) => mol.clone(),
        other => panic!("not a molecule: {other:?}"),
    }
}

// The drawing, spelled as the V2000 molfile the engine reads: real
// drawn coordinates (molfile y rises, the SVG y falls — mirrored),
// parity field 0 (the wedge codes carry the stereo, exactly as the
// oracle reads 2D input), one wedge/hash code per synthesized bond
// oriented centre → neighbour.
fn drawn_molfile(smiles: &str) -> String {
    let molecule = molecule_from_smiles(smiles);
    let (graph, positions) = layout::layout_molecule(&molecule);
    let wedges = svg::synthesize_stereo_wedges(&graph, &positions);

    let mut out = String::new();
    out.push_str("drawn\n  AsciiChem\n\n");
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        graph.atoms.len(),
        graph.bonds.len()
    ));
    for (index, atom) in graph.atoms.iter().enumerate() {
        let position = positions[index];
        let hcount = atom.hydrogens.map_or(0, |count| count + 1);
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0{:>3}{:>3}  0  0  0  0  0  0  0  0\n",
            position.x,
            -position.y,
            0.0,
            atom.element,
            0,
            hcount
        ));
    }
    for (bond_index, bond) in graph.bonds.iter().enumerate() {
        let (from, to, stereo) = match wedges.get(&bond_index) {
            Some(wedge) => (wedge.center, wedge.neighbor, if wedge.hash { 6 } else { 1 }),
            None => (bond.from, bond.to, 0),
        };
        let bond_type = match bond.bond.kind.as_str() {
            "double" => 2,
            "triple" => 3,
            _ => 1,
        };
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0  0  0  0  0\n",
            from + 1,
            to + 1,
            bond_type,
            stereo
        ));
    }
    out.push_str("M  END\n");
    out
}

const SEEDS: &[&str] = &[
    // The issue's enantiomer pair.
    "C[C@H](O)N",
    "C[C@@H](O)N",
    // The corpus's own l-alanine spelling and its mirror.
    "N[C@@H](C)C(=O)O",
    "N[C@H](C)C(=O)O",
    // Four drawn ligands, no stated hydrogen — the completion-free path.
    "F[C@](Cl)(Br)I",
    "F[C@@](Cl)(Br)I",
];

#[test]
fn synthesized_drawings_read_as_their_source_identities() {
    let mut ran = 0usize;
    for smiles in SEEDS {
        let source = asciichem::identity::derive_from_smiles(smiles).expect("source identity");
        let drawn = asciichem::identity::derive_from_molfile(&drawn_molfile(smiles))
            .unwrap_or_else(|error| panic!("{smiles}: the drawn molfile derives: {error}"));
        assert_eq!(
            source.inchi, drawn.inchi,
            "{smiles}: the drawing reads as a different configuration"
        );
        ran += 1;
    }
    assert!(ran >= 6, "the seed battery shrank: {ran}");

    // The pairs really are enantiomers through the engine (a polarity
    // fit is meaningless if the sources collapse to one identity).
    let l = asciichem::identity::derive_from_smiles("C[C@H](O)N").expect("identity");
    let d = asciichem::identity::derive_from_smiles("C[C@@H](O)N").expect("identity");
    assert_ne!(l.inchi, d.inchi, "the enantiomer pair collapsed");

    // The drawings distinguish the pair byte-wise (the mirror test).
    let render = |smiles: &str| svg::render_formula(&asciichem::parse_nodes(smiles).expect("parses"));
    assert_ne!(
        render("@smiles(\"C[C@H](O)N\")"),
        render("@smiles(\"C[C@@H](O)N\")"),
        "the enantiomers render identically"
    );
}
