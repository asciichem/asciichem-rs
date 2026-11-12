//! Shared-corpus conformance: the parser acceptance/rejection contract
//! and the L1 text round-trip gate.
//!
//! Every fixture case that carries `input` (and is not a linter or
//! convention case, which assert diagnostics rather than parsing) must
//! parse or reject exactly as marked. Round-trip-flagged inputs must
//! render back through the canonical text form byte-identically — or
//! be reported as `Unsupported` skips, never silent mismatches. The
//! fixtures live in the sibling `asciichem-tests` checkout (set
//! `ASCIICHEM_CORPUS` to override). Absent corpus = the test reports a
//! skip rather than failing a fresh clone.

use std::fs;
use std::path::PathBuf;

fn corpus_dir() -> Option<PathBuf> {
    let dir = std::env::var("ASCIICHEM_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../asciichem-tests/corpus/fixtures"));
    dir.is_dir().then_some(dir)
}

struct Case {
    id: String,
    input: String,
    parses: bool,
    round_trip: bool,
}

fn load_cases() -> Option<Vec<Case>> {
    let dir = corpus_dir()?;
    let mut cases = Vec::new();
    for entry in fs::read_dir(&dir).expect("read corpus dir") {
        let path = entry.expect("dir entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
            continue;
        }
        let text = fs::read_to_string(&path).expect("read fixture file");
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&text).expect("valid fixture JSON");
        for case in parsed {
            let Some(input) = case.get("input").and_then(|value| value.as_str()) else {
                continue;
            };
            if case.get("lint").is_some() || case.get("convention").is_some() {
                continue;
            }
            cases.push(Case {
                id: case
                    .get("id")
                    .and_then(|value| value.as_str())
                    .unwrap_or("<unnamed>")
                    .to_string(),
                input: input.to_string(),
                parses: case
                    .get("parses")
                    .and_then(serde_json::Value::as_bool)
                    .unwrap_or(false),
                round_trip: case
                    .get("roundTrip")
                    .and_then(serde_json::Value::as_bool)
                    .unwrap_or(false),
            });
        }
    }
    Some(cases)
}

#[test]
fn parser_acceptance_contract() {
    let Some(cases) = load_cases() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };

    let mut failures = Vec::new();
    let mut checked = 0;
    for case in &cases {
        let actual = asciichem::parse(&case.input).is_ok();
        checked += 1;
        if actual != case.parses {
            failures.push(format!(
                "{}: expected parses={}, got {actual} ({:?})",
                case.id, case.parses, case.input
            ));
        }
    }

    assert!(
        checked > 100,
        "expected a substantial corpus, got {checked}"
    );
    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn text_round_trip_contract() {
    let Some(cases) = load_cases() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };

    let mut pass = 0;
    let mut skipped: Vec<String> = Vec::new();
    let mut mismatches = Vec::new();
    for case in cases.iter().filter(|case| case.parses && case.round_trip) {
        match asciichem::parse_to_text(&case.input) {
            Ok(output) if output == case.input => pass += 1,
            Ok(output) => mismatches.push(format!(
                "{}: round-trip mismatch\n  in: {:?}\n out: {:?}",
                case.id, case.input, output
            )),
            Err(error) => skipped.push(format!("{}: {error}", case.id)),
        }
    }

    assert!(
        pass > 0,
        "expected round-trip passes; got none (skips: {skipped:?})"
    );
    assert!(
        mismatches.is_empty(),
        "{} round-trip mismatch(es):\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!(
        "round-trip: {pass} pass, {} skipped (unsupported subset)",
        skipped.len()
    );
}

#[test]
fn holds_no_state_between_parses() {
    // The grammar singleton is immutable; parsing the same input twice
    // gives the same verdict, and a reject does not poison later calls.
    assert!(asciichem::parse("H_2O").is_ok());
    assert!(asciichem::parse("   ").is_err());
    assert!(asciichem::parse("H_2O").is_ok());
    assert_eq!(asciichem::parse_to_text("H_2O").unwrap(), "H_2O");
}

// -- structure interchange (TODO.v2 09) + perception (item 2) ---------
// The corpus's SMILES and perception families gate the shared
// structure semantics: parse/reject, the deterministic writer (exact
// round trip and canonical aromatic form), and the perception pass
// (implicit hydrogens + Hückel aromatic rings, null = declines).

fn load_keyed_cases(key: &str) -> Option<Vec<serde_json::Value>> {
    let dir = corpus_dir()?;
    let mut cases = Vec::new();
    for entry in fs::read_dir(&dir).expect("read corpus dir") {
        let path = entry.expect("dir entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
            continue;
        }
        let text = fs::read_to_string(&path).expect("read fixture file");
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&text).expect("valid fixture JSON");
        for case in parsed {
            if case.get(key).is_some() {
                cases.push(case);
            }
        }
    }
    Some(cases)
}

fn case_id(case: &serde_json::Value) -> String {
    case.get("id")
        .and_then(|value| value.as_str())
        .unwrap_or("<unnamed>")
        .to_string()
}

#[test]
fn smiles_ingestion_contract() {
    let Some(cases) = load_keyed_cases("smiles") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 40,
        "expected the SMILES family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    for case in &cases {
        let id = case_id(case);
        let smiles = case
            .get("smiles")
            .and_then(|value| value.as_str())
            .expect("smiles string");
        let parses = case
            .get("parses")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false);
        match asciichem::smiles::parse_smiles(smiles) {
            Ok(nodes) => {
                if !parses {
                    failures.push(format!("{id}: expected reject, parsed {smiles:?}"));
                    continue;
                }
                let written = asciichem::smiles::write_smiles(&nodes);
                let round_trip = case
                    .get("smilesRoundTrip")
                    .and_then(serde_json::Value::as_bool)
                    .unwrap_or(false);
                if round_trip {
                    match &written {
                        Ok(output) if output == smiles => {}
                        Ok(output) => failures.push(format!(
                            "{id}: round-trip mismatch\n  in: {smiles:?}\n out: {output:?}"
                        )),
                        Err(error) => failures.push(format!("{id}: write failed: {error}")),
                    }
                }
                if let Some(canonical) = case.get("canonicalSmiles").and_then(|v| v.as_str()) {
                    match &written {
                        Ok(output) if output == canonical => {}
                        Ok(output) => failures.push(format!(
                            "{id}: canonical mismatch\n  expected: {canonical:?}\n  got: {output:?}"
                        )),
                        Err(error) => failures.push(format!("{id}: write failed: {error}")),
                    }
                }
            }
            Err(error) => {
                if parses {
                    failures.push(format!("{id}: expected parse, rejected: {error}"));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} SMILES mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn perception_contract() {
    let Some(cases) = load_keyed_cases("perception") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 14,
        "expected the perception family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    for case in &cases {
        let id = case_id(case);
        let smiles = case
            .get("smiles")
            .and_then(|value| value.as_str())
            .expect("smiles string");
        let expected = case.get("perception").expect("perception object");

        let nodes = asciichem::smiles::parse_smiles(smiles)
            .unwrap_or_else(|error| panic!("{id}: fixture SMILES must parse: {error}"));
        let molecule = match nodes.first() {
            Some(asciichem::model::Node::Molecule(molecule)) => molecule,
            other => panic!("{id}: expected a molecule node, got {other:?}"),
        };
        let (atoms, edges) = asciichem::structure::build_graph(molecule);

        let hydrogens: Vec<Option<u32>> = atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| {
                asciichem::perception::implicit_hydrogens(
                    atom,
                    asciichem::perception::hydrogen_valence_sum(index, &edges),
                )
            })
            .collect();
        let hydrogens_json = serde_json::Value::Array(
            hydrogens
                .iter()
                .map(|h| match h {
                    Some(count) => serde_json::Value::Number((*count).into()),
                    None => serde_json::Value::Null,
                })
                .collect(),
        );
        if &hydrogens_json
            != expected
                .get("implicitHydrogens")
                .expect("implicitHydrogens")
        {
            failures.push(format!(
                "{id}: implicit hydrogens mismatch\n  expected: {}\n  got: {}",
                expected
                    .get("implicitHydrogens")
                    .expect("implicitHydrogens"),
                hydrogens_json
            ));
        }

        let mut rings: Vec<Vec<usize>> = asciichem::perception::aromatic_rings(&atoms, &edges)
            .iter()
            .map(|ring| {
                let mut ring = ring.clone();
                ring.sort_unstable();
                ring
            })
            .collect();
        rings.sort();
        let rings_json = serde_json::Value::Array(
            rings
                .iter()
                .map(|ring| {
                    serde_json::Value::Array(
                        ring.iter()
                            .map(|i| serde_json::Value::Number((*i).into()))
                            .collect(),
                    )
                })
                .collect(),
        );
        if &rings_json != expected.get("aromaticRings").expect("aromaticRings") {
            failures.push(format!(
                "{id}: aromatic rings mismatch\n  expected: {}\n  got: {}",
                expected.get("aromaticRings").expect("aromaticRings"),
                rings_json
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} perception mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// -- molfile interchange (TODO.v2 09, slice 2) --------------------------
// The corpus's molfile family gates V2000 ingest/emission: parse/reject,
// atom and bond counts, and — for `molfileRoundTrip` — a parse/write/
// parse cycle with the same graph (the runner contract's shape
// equality, sorted edges). Stereo codes 1/6 ride the same gate.

// -- CML round-trip (TODO.v2 17 slice 7) --------------------------------
// The corpus's cml family gates the three-way contract: parse ->
// to_cml -> parse_cml -> canonical text, byte-identical to the input.
// Both XML directions ride the linked leptris engine.

#[test]
fn cml_three_way_round_trip_contract() {
    let Some(cases) = load_keyed_cases("cmlRoundTrip") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 20,
        "expected the cml round-trip family, got {} cases",
        cases.len()
    );

    let mut round_trips = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case_id(case);
        let input = case
            .get("input")
            .and_then(|value| value.as_str())
            .expect("cml case input");
        let nodes = match asciichem::parse_nodes(input) {
            Ok(nodes) => nodes,
            Err(error) => {
                failures.push(format!("{id}: model parse failed: {error}"));
                continue;
            }
        };
        let xml = match asciichem::cml::to_cml(&nodes) {
            Ok(xml) => xml,
            Err(error) => {
                failures.push(format!("{id}: CML emit failed: {error}"));
                continue;
            }
        };
        let back = match asciichem::cml::parse_cml(&xml) {
            Ok(back) => back,
            Err(error) => {
                failures.push(format!("{id}: CML read failed: {error}\n  xml: {xml}"));
                continue;
            }
        };
        let text = back
            .iter()
            .map(|node| node.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        if text == input {
            round_trips += 1;
        } else {
            failures.push(format!(
                "{id}: three-way round-trip mismatch\n  in: {input:?}\n out: {text:?}"
            ));
        }
    }
    assert!(round_trips > 0, "expected cml round-trip passes; got none");
    assert!(
        failures.is_empty(),
        "{} cml mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn edge_shape(edges: &[asciichem::structure::Edge]) -> Vec<(usize, usize, String)> {
    let mut shape: Vec<(usize, usize, String)> = edges
        .iter()
        .map(|edge| {
            (
                edge.from.min(edge.to),
                edge.from.max(edge.to),
                edge.kind.clone(),
            )
        })
        .collect();
    shape.sort();
    shape
}

#[test]
fn molfile_ingestion_contract() {
    let Some(cases) = load_keyed_cases("molfile") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 6,
        "expected the molfile family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut round_trips = 0;
    for case in &cases {
        let id = case_id(case);
        let molfile = case
            .get("molfile")
            .and_then(|value| value.as_str())
            .expect("molfile string");
        let parses = case
            .get("parses")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false);
        match asciichem::molfile::parse_molfile(molfile) {
            Ok(molecule) => {
                if !parses {
                    failures.push(format!("{id}: expected reject, parsed"));
                    continue;
                }
                let (atoms, edges) = asciichem::structure::build_graph(&molecule);
                for (actual, expected, label) in [
                    (
                        atoms.len(),
                        case.get("atoms").and_then(|v| v.as_u64()),
                        "atoms",
                    ),
                    (
                        edges.len(),
                        case.get("bonds").and_then(|v| v.as_u64()),
                        "bonds",
                    ),
                ] {
                    match expected {
                        Some(expected) if actual as u64 == expected => {}
                        _ => failures.push(format!(
                            "{id}: {label} count mismatch (expected {expected:?}, got {actual})"
                        )),
                    }
                }
                let round_trip = case
                    .get("molfileRoundTrip")
                    .and_then(serde_json::Value::as_bool)
                    .unwrap_or(false);
                if round_trip {
                    let shape = edge_shape(&edges);
                    match asciichem::molfile::write_molfile(&molecule, "") {
                        Ok(written) => {
                            let (again_atoms, again_edges) = asciichem::structure::build_graph(
                                &asciichem::molfile::parse_molfile(&written).unwrap_or_else(
                                    |error| panic!("{id}: re-parse failed: {error}"),
                                ),
                            );
                            if again_atoms.len() != atoms.len() || edge_shape(&again_edges) != shape
                            {
                                failures
                                    .push(format!("{id}: molfile round-trip changed the graph"));
                            }
                            round_trips += 1;
                        }
                        Err(error) => failures.push(format!("{id}: write failed: {error}")),
                    }
                }
            }
            Err(error) => {
                if parses {
                    failures.push(format!("{id}: expected parse, rejected: {error}"));
                }
            }
        }
    }
    assert!(
        round_trips > 0,
        "expected molfile round-trip passes; got none"
    );
    assert!(
        failures.is_empty(),
        "{} molfile mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
