//! Unit contract for the SMILES module (TODO.v2 09 + 13 item 2):
//! the deterministic writer's stability fixed point, canonical
//! aromatic perception, and stereo carry-through — the properties the
//! TS and Python implementations spec in their unit suites, mirrored
//! here for the fourth implementation.

use asciichem::model::Node;
use asciichem::smiles::{parse_smiles, write_smiles};
use asciichem::ParseError;

fn round_trip(smiles: &str) -> String {
    let nodes = parse_smiles(smiles).unwrap_or_else(|error| panic!("{smiles:?}: {error}"));
    write_smiles(&nodes).unwrap_or_else(|error| panic!("{smiles:?}: write failed: {error}"))
}

// Canonicalization maps many models to one spelling (Kekulé →
// aromatic), so model equality is not the invariant — composition and
// connectivity must survive and the form must be a fixed point.
#[test]
fn writer_is_stable_and_structure_preserving() {
    for smiles in [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "c1ccc2ccccc2c1",
        "CC(C)(C)C",
        // Ring-closure digits stay out of branches (#103): a
        // closure-bearing child continues the chain.
        "C1=CCC(C)C(C)C1",
        "C1CC2CCC1C2",
        "C1CC(C)CCC1",
    ] {
        let once = round_trip(smiles);
        let twice = round_trip(&once);
        assert_eq!(twice, once, "{smiles:?}: not a fixed point");

        let before = parse_smiles(smiles).expect("parses");
        let after = parse_smiles(&once).expect("re-parses");
        assert_eq!(
            before.len(),
            after.len(),
            "{smiles:?}: component count changed"
        );
        for (before, after) in before.iter().zip(after.iter()) {
            let (Node::Molecule(before), Node::Molecule(after)) = (before, after) else {
                panic!("expected molecule nodes");
            };
            let atom_count = |nodes: &[Node]| {
                nodes
                    .iter()
                    .filter(|node| matches!(node, Node::Atom(_)))
                    .count()
            };
            assert_eq!(
                atom_count(&before.nodes),
                atom_count(&after.nodes),
                "{smiles:?}: atom count changed through {once:?}"
            );
        }
    }
}

// Rings Hückel-perceived aromatic are canonicalized lowercase; rings
// that fail the count (cyclopentadiene's four pi, quinone's exocyclic
// doubles, cyclopropene's three-membered case) stay Kekulé.
#[test]
fn canonicalizes_perceived_aromatic_rings() {
    for (kekule, canonical) in [
        ("C1=CC=CC=C1", "c1ccccc1"),
        ("CC(=O)OC1=CC=CC=C1C(=O)O", "CC(=O)Oc1ccccc1C(=O)O"),
        ("ClC1=CC=CC=C1", "Clc1ccccc1"),
        ("C1=CC=NC=C1", "c1ccncc1"),
        ("C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"),
        // Hückel five-membered rings (corpus 0.8.0).
        ("C1=COC=C1", "c1cocc1"),
        ("C1=CNC=C1", "c1c[nH]cc1"),
        ("C1=CSC=C1", "c1cscc1"),
    ] {
        assert_eq!(round_trip(kekule), canonical, "kekule {kekule:?}");
        assert_eq!(round_trip(canonical), canonical, "canonical {canonical:?}");
    }
}

#[test]
fn keeps_non_aromatic_rings_kekule() {
    for smiles in ["O=C1C=CC=CC1", "O=C1C=CC(=O)C=C1", "C1=CCC=C1", "C1=CC1"] {
        assert_eq!(round_trip(smiles), smiles, "{smiles:?}");
    }
}

// Stereo moved from the v1 rejections to supported features
// (corpus 0.5.0): @/@@ parse into Atom.chirality, / \ into
// Bond.direction — and the writer carries both back out.
#[test]
fn carries_stereo_through_write() {
    // The written @/@@ is corrected against the written neighbour
    // order; both enantiomer spellings surface an @ marker.
    assert!(round_trip("[C@H](C)C").contains('@'));
    assert!(round_trip("[C@@H](C)C").contains('@'));
    assert_eq!(round_trip("C/C=C/C"), "C/C=C/C");
    assert_eq!(round_trip("F/C=C\\F"), "F/C=C\\F");
}

#[test]
fn normalizes_charge_spelling() {
    assert_eq!(round_trip("[Fe++]"), "[Fe+2]");
    assert_eq!(round_trip("[Fe+3]"), "[Fe+3]");
    assert_eq!(round_trip("[O-2]"), "[O-2]");
}

#[test]
fn splits_dot_components() {
    let nodes = parse_smiles("O.CCO").expect("parses");
    assert_eq!(nodes.len(), 2);
    assert_eq!(round_trip("O.CCO"), "O.CCO");
}

#[test]
fn rejects_out_of_subset_forms() {
    for smiles in ["C1CCC", "C=1CCCCC=1", "CC!X", "CC..O", "[C", "C*O"] {
        match parse_smiles(smiles) {
            Err(ParseError { .. }) => {}
            Ok(_) => panic!("{smiles:?}: expected reject"),
        }
    }
}
