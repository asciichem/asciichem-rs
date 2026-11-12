//! Molfile (CTfile V2000) interchange — the corpus gate plus the
//! unit-level properties behind it: stereo codes 1/6 (wedge/hash),
//! chirality crossing as a wedge, charges (property block and the
//! legacy column), isotopes, aromatic bond types, rejects, and
//! byte-stable deterministic re-emission. Semantics mirror the
//! reference implementation's AsciiChem::Molfile.

use asciichem::model::{Atom, Bond, Molecule, Node};
use asciichem::structure::build_graph;

const ETHANOL: &str = "ethanol\n  fixture\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.2500    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4645   -0.0375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1790    0.3750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\nM  END\n";

const AMMONIUM: &str = "ammonium\n  fixture\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9500    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4750    0.8228    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4750   -0.8228    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.9500 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0  0  0  0  0\n  1  4  1  0  0  0  0  0  0  0  0\n  1  5  1  0  0  0  0  0  0  0  0\nM  CHG  1   1   1\nM  END\n";

const WEDGE_STEREO: &str = "2-butanol-wedge\n  fixture\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   -1.5000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -0.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  3  4  1  0  0  0  0  0  0  0  0\n  2  5  1  0  0  0  0  0  0  0  0\nM  END\n";

const BENZENE: &str = "benzene\n  fixture\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0  0  0  0  0\n  2  3  4  0  0  0  0  0  0  0  0\n  3  4  4  0  0  0  0  0  0  0  0\n  4  5  4  0  0  0  0  0  0  0  0\n  5  6  4  0  0  0  0  0  0  0  0\n  6  1  4  0  0  0  0  0  0  0  0\nM  END\n";

fn parse(molfile: &str) -> Molecule {
    asciichem::molfile::parse_molfile(molfile).expect("fixture parses")
}

fn molecule_from_smiles(smiles: &str) -> Molecule {
    match asciichem::smiles::parse_smiles(smiles)
        .expect("smiles parses")
        .first()
    {
        Some(Node::Molecule(molecule)) => molecule.clone(),
        other => panic!("expected a molecule node, got {other:?}"),
    }
}

fn edges_of(molecule: &Molecule) -> Vec<(usize, usize, String)> {
    let (_, edges) = build_graph(molecule);
    let mut shape: Vec<(usize, usize, String)> = edges
        .iter()
        .map(|edge| (edge.from, edge.to, edge.kind.clone()))
        .collect();
    shape.sort();
    shape
}

// -- stereo codes 1/6 (the deliberately-deferred gate) ------------------

#[test]
fn stereo_codes_become_wedge_and_hash() {
    // Bond 1→2 carries stereo code 1 (up/wedge); give bond 2→3 code 6
    // (down/hash) — both ride the single-bond type.
    let molfile = WEDGE_STEREO.replace("  2  3  1  0", "  2  3  1  6");
    let molecule = parse(&molfile);
    assert_eq!(
        edges_of(&molecule),
        vec![
            (0, 1, "wedge".to_string()),
            (1, 2, "hash".to_string()),
            (1, 4, "single".to_string()),
            (2, 3, "single".to_string()),
        ]
    );
}

#[test]
fn writer_emits_stereo_codes_one_and_six() {
    let mut wedge_atom = Atom {
        element: "C".into(),
        ..Atom::default()
    };
    wedge_atom.chirality = Some("TH1".into());
    let hash_atom = Atom {
        element: "C".into(),
        ..Atom::default()
    };
    let molecule = Molecule {
        nodes: vec![
            Node::Atom(wedge_atom),
            Node::Bond(Bond {
                kind: "single".into(),
                direction: None,
            }),
            Node::Atom(hash_atom),
            Node::Bond(Bond {
                kind: "single".into(),
                direction: None,
            }),
            Node::Atom(Atom {
                element: "O".into(),
                ..Atom::default()
            }),
        ],
        ..Molecule::default()
    };

    let written = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    let lines: Vec<&str> = written.lines().collect();
    let counts = lines
        .iter()
        .position(|line| line.ends_with("V2000"))
        .expect("counts line");
    // TH1 → wedge (type 1, stereo 1) on the stereocentre's first bond,
    // oriented centre → neighbour; TH2 would fall (stereo 6).
    assert_eq!(lines[counts + 4], "  1  2  1  1  0  0  0  0  0  0  0");

    // And the hash kind round-trips: a hash edge writes code 6.
    let hash_molecule = Molecule {
        nodes: vec![
            Node::Atom(Atom {
                element: "C".into(),
                ..Atom::default()
            }),
            Node::Bond(Bond {
                kind: "hash".into(),
                direction: None,
            }),
            Node::Atom(Atom {
                element: "C".into(),
                ..Atom::default()
            }),
        ],
        ..Molecule::default()
    };
    let written = asciichem::molfile::write_molfile(&hash_molecule, "").expect("writes");
    assert!(written.contains("  1  2  1  6  0  0  0  0  0  0  0"));
}

#[test]
fn chirality_crosses_as_a_wedge_and_returns_as_parity() {
    // The writer crosses chirality as a wedge on the stereocentre's
    // first bond (the oracle ignores the atom-parity field alone); the
    // reader folds that wedge back to a plain bond once the parity
    // field has restored the chirality.
    let molecule = molecule_from_smiles("[C@H](O)(N)C");
    let written = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    assert!(
        written.contains("  1  2  1  1  0  0  0  0  0  0  0"),
        "{written}"
    );

    let reparsed = parse(&written);
    let atoms = match reparsed.nodes.first() {
        Some(Node::Atom(atom)) => Some(atom),
        _ => None,
    };
    let chirality = atoms.and_then(|atom| atom.chirality.clone());
    assert_eq!(chirality.as_deref(), Some("TH1"));
    assert_eq!(
        edges_of(&reparsed),
        edges_of(&molecule),
        "the wedge folds back to a plain single bond"
    );

    // TH2 falls: stereo code 6, restored as TH2.
    let other = molecule_from_smiles("[C@@H](O)(N)C");
    let written = asciichem::molfile::write_molfile(&other, "").expect("writes");
    assert!(
        written.contains("  1  2  1  6  0  0  0  0  0  0  0"),
        "{written}"
    );
    let reparsed = parse(&written);
    let atoms = match reparsed.nodes.first() {
        Some(Node::Atom(atom)) => Some(atom),
        _ => None,
    };
    assert_eq!(
        atoms.and_then(|atom| atom.chirality.clone()).as_deref(),
        Some("TH2")
    );
}

// -- charges and isotopes ------------------------------------------------

#[test]
fn property_block_charges_apply() {
    let molecule = parse(AMMONIUM);
    match molecule.nodes.first() {
        Some(Node::Atom(atom)) => assert_eq!(atom.charge.as_deref(), Some("+")),
        other => panic!("expected an atom, got {other:?}"),
    }
    let written = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    assert!(written.contains("M  CHG  1   1   1"), "{written}");
}

#[test]
fn legacy_charge_column_applies_and_chg_overrides() {
    // Charge column (0-based 36..39) code 3 → +3; `M  CHG` overrides
    // the column. The atom line: coords, symbol, mass-diff " 0",
    // charge "  3", then zeros.
    let legacy = "charged\n  fixture\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  3  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    let molecule = parse(legacy);
    match molecule.nodes.first() {
        Some(Node::Atom(atom)) => assert_eq!(atom.charge.as_deref(), Some("3+")),
        other => panic!("expected an atom, got {other:?}"),
    }

    let overridden = "charged\n  fixture\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  3  0  0  0  0  0  0  0  0  0  0\nM  CHG  1   1  -2\nM  END\n";
    let molecule = parse(overridden);
    match molecule.nodes.first() {
        Some(Node::Atom(atom)) => assert_eq!(atom.charge.as_deref(), Some("2-")),
        other => panic!("expected an atom, got {other:?}"),
    }
}

#[test]
fn isotope_property_block_applies() {
    let iso = "iso\n  fixture\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  ISO  1   1  14\nM  END\n";
    let molecule = parse(iso);
    match molecule.nodes.first() {
        Some(Node::Atom(atom)) => assert_eq!(atom.isotope.as_deref(), Some("14")),
        other => panic!("expected an atom, got {other:?}"),
    }
    let written = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    assert!(written.contains("M  ISO  1   1  14"), "{written}");
}

// -- aromaticity ---------------------------------------------------------

#[test]
fn aromatic_bond_type_marks_atoms_and_canonicalizes() {
    let molecule = parse(BENZENE);
    let (atoms, edges) = build_graph(&molecule);
    assert_eq!(edges.len(), 6);
    assert!(edges.iter().all(|edge| edge.kind == "aromatic"));
    assert!(atoms.iter().all(|atom| atom.aromatic == Some(true)));

    // The canonical aromatic SMILES falls out of the shared writer.
    let smiles = asciichem::smiles::write_smiles(&[Node::Molecule(molecule)]).expect("writes");
    assert_eq!(smiles, "c1ccccc1");
}

// -- coordinates ---------------------------------------------------------

#[test]
fn authored_coordinates_reemit_byte_identically() {
    // The fixture's atom and bond blocks were emitted by the reference
    // writer's fixed formats; re-emitting the parsed model reproduces
    // them byte-for-byte (only the header — name and program line —
    // is the writer's own).
    let molecule = parse(ETHANOL);
    let written = asciichem::molfile::write_molfile(&molecule, "ethanol").expect("writes");
    let written_blocks: Vec<&str> = written.lines().skip(3).collect();
    let fixture_blocks: Vec<&str> = ETHANOL.lines().skip(3).collect();
    assert_eq!(written_blocks, fixture_blocks);
}

#[test]
fn missing_coordinates_fall_back_to_the_layout() {
    // A model without authored coordinates (text/SMILES origin) gets
    // deterministic layout positions — non-zero for a bonded chain.
    let molecule = molecule_from_smiles("CCC");
    let written = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    let atom_lines: Vec<&str> = written.lines().skip(4).take(3).collect();
    assert!(atom_lines[0] != atom_lines[1], "{written}");
}

// -- determinism ---------------------------------------------------------

#[test]
fn re_emission_is_a_byte_stable_fixed_point() {
    for molfile in [ETHANOL, AMMONIUM, WEDGE_STEREO, BENZENE] {
        let molecule = parse(molfile);
        let once = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
        let twice = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
        assert_eq!(once, twice, "repeated writes diverge for {molfile:?}");

        let round = asciichem::molfile::write_molfile(&parse(&once), "").expect("writes");
        assert_eq!(
            round, once,
            "re-emission is not a fixed point for {molfile:?}"
        );
    }
}

#[test]
fn chirality_write_is_a_fixed_point() {
    let molecule = molecule_from_smiles("[C@H](O)(N)C");
    let once = asciichem::molfile::write_molfile(&molecule, "").expect("writes");
    let round = asciichem::molfile::write_molfile(&parse(&once), "").expect("writes");
    assert_eq!(round, once, "{once}");
}

// -- declines ------------------------------------------------------------

#[test]
fn rejects_are_clean() {
    // Too short.
    assert_eq!(
        asciichem::molfile::parse_molfile("a\nb\nc\n")
            .err()
            .map(|e| e.message),
        Some("molfile too short".into())
    );
    // Malformed counts line.
    assert!(parse_rejects(
        "bad\n  x\n\nzzz  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
    ));
    // Truncated atom block.
    assert!(parse_rejects(
        "trunc\n  fixture\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
    ));
    // Truncated bond block.
    assert!(parse_rejects(
        "trunc\n  fixture\n\n  1  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"
    ));
    // Atom without an element symbol.
    assert!(parse_rejects(
        "noelem\n  fixture\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000     0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n"
    ));
    // Bond index out of range.
    assert!(parse_rejects(
        "oob\n  fixture\n\n  1  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\nM  END\n"
    ));
    // Unsupported bond type.
    assert!(parse_rejects(
        "badtype\n  fixture\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  9  0  0  0  0  0  0  0  0\nM  END\n"
    ));
}

fn parse_rejects(molfile: &str) -> bool {
    matches!(
        asciichem::molfile::parse_molfile(molfile),
        Err(asciichem::ParseError { .. })
    )
}

#[test]
fn formulas_are_not_structures() {
    // A molecule with adjacency but no bonds (the undashed text form)
    // is declined — explicit-bond or SMILES/molfile ingestion only.
    let formula = asciichem::parse_nodes("HHO").expect("parses");
    match formula.first() {
        Some(Node::Molecule(molecule)) => {
            assert!(asciichem::molfile::write_molfile(molecule, "").is_err());
        }
        other => panic!("expected a molecule node, got {other:?}"),
    }
    // So is fragmenting adjacency: one explicit bond, the rest implicit.
    let fragmented = asciichem::parse_nodes("C-C(C)C").expect("parses");
    match fragmented.first() {
        Some(Node::Molecule(molecule)) => {
            assert!(asciichem::molfile::write_molfile(molecule, "").is_err());
        }
        other => panic!("expected a molecule node, got {other:?}"),
    }
}
