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
//!
//! The parity families (corpus 0.11.0) gate the remaining shared
//! contracts from the same fixtures: linter-expectations (diagnostics
//! both ways, engine-scaled), resolver (recorded exchanges replayed
//! through the Fetch seam with the clock pinned), citation (byte-equal
//! bibitem goldens), and cli (the offline command contract, invoked
//! in-process). Posture and `pending` cases skip AND report — a
//! missing capability is a visible report, never a silent gap.

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

// -- mechanism family (TODO.v2/13 item 4) --------------------------------
// The `flows` shape gates the electron-flow markers per step label:
// the parsed mechanism's steps must map label -> [{source, target}]
// exactly as the fixture's `flows` object — every step present,
// flow-less steps carrying an empty list, the arrows in order. The
// canonical text round trip rides the generic input contract.

#[test]
fn mechanism_flows_contract() {
    let Some(cases) = load_family("mechanisms.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 3,
        "expected the mechanisms family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut ran = 0;
    for case in &cases {
        let id = case_id(case);
        ran += 1;
        let input = case
            .get("input")
            .and_then(serde_json::Value::as_str)
            .expect("mechanism case input");
        let nodes = match asciichem::parse_nodes(input) {
            Ok(nodes) => nodes,
            Err(error) => {
                failures.push(format!("{id}: mechanism parse failed: {error}"));
                continue;
            }
        };
        let Some(asciichem::model::Node::Mechanism(mechanism)) = nodes.first() else {
            failures.push(format!(
                "{id}: expected a mechanism node, got {:?}",
                nodes.first()
            ));
            continue;
        };
        let actual = serde_json::Value::Object(
            mechanism
                .steps
                .iter()
                .map(|step| {
                    let flows = serde_json::Value::Array(
                        step.flows
                            .iter()
                            .map(|flow| {
                                serde_json::json!({
                                    "source": flow.source,
                                    "target": flow.target,
                                })
                            })
                            .collect(),
                    );
                    (step.label.clone(), flows)
                })
                .collect(),
        );
        let expected = case.get("flows").expect("flows object");
        if &actual != expected {
            failures.push(format!(
                "{id}: flows mismatch\n  expected: {expected}\n  got: {actual}"
            ));
        }
    }
    eprintln!("mechanisms: {ran} run");
    assert!(
        failures.is_empty(),
        "{} mechanism mismatch(es):\n{}",
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

// -- parity families (corpus 0.11.0; TODO.v2/17) ------------------------
// The linter-expectations, resolver, citation, and cli families gate
// the remaining shared contracts from the SHARED fixtures — the same
// match semantics the Ruby/TS/PY runners run. This implementation's
// posture: the chemomer identity engine is IN-PROCESS (an unconditional
// dependency), so the linter's with-engine cases run everywhere CI
// runs and the corpus's engine:"none" cases skip and report. Cases
// `pending` for this implementation skip and report too — a missing
// capability stays visible in the output until it lands.

const IMPLEMENTATION: &str = "asciichem-rs";

/// This run's identity-engine posture: the engine is built in.
const ENGINE_POSTURE: &str = "inchi";

/// The retrieval time the corpus pins so replayed records and the
/// citation goldens that consume them are stable bytes.
const PINNED_RETRIEVED_AT: &str = "2026-01-05T09:30:00Z";

/// The replay key substituted for the fixture's `"{api-key}"`
/// placeholder — never a real credential.
const FAKE_KEY: &str = "corpus-replay-test-key";

/// Loads one named family fixture. The corpus being absent altogether
/// is a skip (the fresh-clone posture); a present corpus missing a
/// family file is a hard failure — the pin must carry suite 0.11.0.
fn load_family(name: &str) -> Option<Vec<serde_json::Value>> {
    let path = corpus_dir()?.join(name);
    assert!(
        path.is_file(),
        "{name} missing from the corpus — the pin does not carry suite 0.11.0"
    );
    let text = fs::read_to_string(&path).expect("read fixture file");
    Some(serde_json::from_str(&text).expect("valid fixture JSON"))
}

/// The implementation ids that cannot run the case yet.
fn pending_for(case: &serde_json::Value) -> Vec<&str> {
    case.get("pending")
        .and_then(serde_json::Value::as_array)
        .map(|list| {
            list.iter()
                .filter_map(serde_json::Value::as_str)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

fn is_pending(case: &serde_json::Value) -> bool {
    pending_for(case).contains(&IMPLEMENTATION)
}

fn engine_marker(case: &serde_json::Value) -> Option<&str> {
    case.get("engine").and_then(serde_json::Value::as_str)
}

fn string_list(case: &serde_json::Value, key: &str) -> Vec<String> {
    case.get(key)
        .and_then(serde_json::Value::as_array)
        .map(|list| {
            list.iter()
                .filter_map(serde_json::Value::as_str)
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

// -- linter-expectations family -----------------------------------------
// The diagnostics contract runs BOTH ways: every expectation matches at
// least one actual diagnostic (severity, message containing
// `messageIncludes`, and — when given — the same node label) AND every
// actual diagnostic matches at least one expectation. `lint: []`
// demands a completely clean run. Molecule-level node labels are
// implementation-local (the schema: the reference prints "Molecule",
// this port prints the hill formula), so the node check only fires
// where the expectation carries one.

fn lint_matches(want: &serde_json::Value, diagnostic: &asciichem::linter::Diagnostic) -> bool {
    let severity = want
        .get("severity")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let fragment = want
        .get("messageIncludes")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let node_ok = want
        .get("node")
        .and_then(serde_json::Value::as_str)
        .map(|node| diagnostic.label.as_deref() == Some(node))
        .unwrap_or(true);
    diagnostic.severity.as_str() == severity && diagnostic.message.contains(fragment) && node_ok
}

fn render_diagnostics(diagnostics: &[asciichem::linter::Diagnostic]) -> String {
    diagnostics
        .iter()
        .map(|diagnostic| {
            format!(
                "({}, {:?}, {:?})",
                diagnostic.severity,
                diagnostic.label,
                diagnostic.message.chars().take(120).collect::<String>()
            )
        })
        .collect::<Vec<_>>()
        .join(", ")
}

#[test]
fn linter_expectations_contract() {
    let Some(cases) = load_family("linter-expectations.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 25,
        "expected the linter-expectations family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut ran = 0;
    let mut skipped_posture: Vec<String> = Vec::new();
    let mut skipped_pending: Vec<String> = Vec::new();
    for case in &cases {
        let id = case_id(case);
        if let Some(marker) = engine_marker(case) {
            if marker != ENGINE_POSTURE {
                eprintln!("SKIP (engine posture {marker} != {ENGINE_POSTURE}): {id}");
                skipped_posture.push(id);
                continue;
            }
        }
        if is_pending(case) {
            eprintln!(
                "SKIP (pending for {IMPLEMENTATION}): {id} — {}",
                pending_for(case).join(", ")
            );
            skipped_pending.push(id);
            continue;
        }
        ran += 1;
        let input = case
            .get("input")
            .and_then(serde_json::Value::as_str)
            .expect("linter case input");
        if !case
            .get("parses")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false)
        {
            if asciichem::parse(input).is_ok() {
                failures.push(format!("{id}: expected reject, parsed {input:?}"));
            }
            continue;
        }
        let nodes = match asciichem::parse_nodes(input) {
            Ok(nodes) => nodes,
            Err(error) => {
                failures.push(format!("{id}: model parse failed: {error}"));
                continue;
            }
        };
        let actual = asciichem::linter::run(&nodes);
        let expected = case
            .get("lint")
            .and_then(serde_json::Value::as_array)
            .expect("lint expectations");
        for want in expected {
            if !actual
                .iter()
                .any(|diagnostic| lint_matches(want, diagnostic))
            {
                failures.push(format!(
                    "{id}: no diagnostic matches {want}; actual: {}",
                    render_diagnostics(&actual)
                ));
            }
        }
        for diagnostic in &actual {
            if !expected.iter().any(|want| lint_matches(want, diagnostic)) {
                failures.push(format!(
                    "{id}: unexpected diagnostic ({}, {:?}, {:?})",
                    diagnostic.severity,
                    diagnostic.label,
                    diagnostic.message.chars().take(120).collect::<String>()
                ));
            }
        }
    }
    eprintln!(
        "linter-expectations: {ran} run, {} skipped (posture: {skipped_posture:?}), \
         {} skipped (pending: {skipped_pending:?})",
        skipped_posture.len(),
        skipped_pending.len()
    );
    assert!(
        failures.is_empty(),
        "{} linter mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn reports_the_posture_skipped_linter_cases() {
    // The skip is part of the report, never silent: the corpus's
    // no-engine cases are listed here by id, and the list changes only
    // when the corpus's posture markers change.
    let Some(cases) = load_family("linter-expectations.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    let skipped: Vec<String> = cases
        .iter()
        .filter(|case| {
            engine_marker(case)
                .map(|marker| marker != ENGINE_POSTURE)
                .unwrap_or(false)
        })
        .map(case_id)
        .collect();
    assert_eq!(
        skipped,
        [
            "lint/inchi-molecule-mismatch",
            "lint/identity-skip-no-engine",
            "lint/inchi-unanalysable-no-engine",
        ],
        "the posture-skip list changed — recheck this runner's posture"
    );
}

// -- resolver family -----------------------------------------------------
// The recorded request/response is replayed through the named source
// adapter — the adapter must build exactly request.url for (source,
// convention, value) and send the adapter-declared request.headers
// ("{api-key}" is the caller's key placeholder). Feeding response.body
// must yield semantically the expectation's substance wire record —
// every identifier and property carrying its provenance with the
// pinned retrievedAt — plus structureParses. The cache markers pin the
// shared caching contract (ttlSeconds, keyFormula, servesOffline). A
// skipped outcome (the key-absent 401) drops the source without
// poisoning the others.

/// One captured request: the URL and the exact header pairs.
type CapturedCall = (String, Vec<(String, String)>);

/// Answers with the recorded response; captures every request for the
/// URL/header assertions (a value fake carrying recorded data).
struct ReplayFetch<'a> {
    response: &'a serde_json::Value,
    calls: std::cell::RefCell<Vec<CapturedCall>>,
}

impl<'a> ReplayFetch<'a> {
    fn new(response: &'a serde_json::Value) -> Self {
        ReplayFetch {
            response,
            calls: std::cell::RefCell::new(Vec::new()),
        }
    }
}

impl asciichem::resolver::Fetch for ReplayFetch<'_> {
    fn get(
        &self,
        url: &str,
        headers: &[(&str, &str)],
    ) -> Result<Option<String>, asciichem::resolver::ResolverError> {
        self.calls.borrow_mut().push((
            url.to_string(),
            headers
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
        ));
        let status = self
            .response
            .get("status")
            .and_then(serde_json::Value::as_i64)
            .unwrap_or_default();
        if status == 404 {
            return Ok(None);
        }
        if status != 200 {
            return Err(asciichem::resolver::ResolverError::Source(format!(
                "HTTP {status}: replay"
            )));
        }
        Ok(Some(match self.response.get("body") {
            Some(serde_json::Value::String(text)) => text.clone(),
            Some(other) => serde_json::to_string(other).expect("recorded body serializes"),
            None => String::new(),
        }))
    }
}

/// Network access from a cache-seeded lookup would be a bug.
struct OfflineFetch;

impl asciichem::resolver::Fetch for OfflineFetch {
    fn get(
        &self,
        url: &str,
        _headers: &[(&str, &str)],
    ) -> Result<Option<String>, asciichem::resolver::ResolverError> {
        Err(asciichem::resolver::ResolverError::Source(format!(
            "network must not be touched: {url}"
        )))
    }
}

/// The headers the adapter must send: the fixture's declared set with
/// the key placeholder substituted (order-insensitive comparison).
fn expected_headers(case: &serde_json::Value) -> std::collections::BTreeMap<String, String> {
    let mut map = std::collections::BTreeMap::new();
    for (name, value) in case
        .get("request")
        .and_then(|request| request.get("headers"))
        .and_then(serde_json::Value::as_object)
        .map(|object| object.iter())
        .into_iter()
        .flatten()
    {
        let value = value
            .as_str()
            .map(|text| if text == "{api-key}" { FAKE_KEY } else { text })
            .unwrap_or_default();
        map.insert(name.clone(), value.to_string());
    }
    map
}

fn captured_headers(call: &CapturedCall) -> std::collections::BTreeMap<String, String> {
    call.1.iter().cloned().collect()
}

/// The named source adapter, with the caller's key when the recorded
/// request carries the `{api-key}` placeholder and no key when it does
/// not (the key-absent case records the request WITHOUT the X-API-KEY
/// header). Building the adapters explicitly — never through the
/// shipped registry's env/file key chain — keeps the family hermetic
/// whatever the machine carries.
fn replay_adapter(case: &serde_json::Value) -> Box<dyn asciichem::resolver::Adapter> {
    let keyed = expected_headers(case)
        .values()
        .any(|value| value == FAKE_KEY);
    match case
        .get("source")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "common_chemistry" => Box::new(if keyed {
            asciichem::resolver::common_chemistry::CommonChemistry::with_key(FAKE_KEY)
        } else {
            asciichem::resolver::common_chemistry::CommonChemistry::keyless()
        }),
        "pubchem" => Box::new(asciichem::resolver::pubchem::PubChem),
        "opsin" => Box::new(asciichem::resolver::opsin::Opsin),
        other => panic!("unknown resolver fixture source: {other}"),
    }
}

/// The substance's canonical wire record, compared semantically (key
/// order cannot fail a case).
fn wire_record(substance: &asciichem::resolver::Substance) -> serde_json::Value {
    serde_json::to_value(substance.to_wire()).expect("wire record serializes")
}

fn temp_corpus_cache(tag: &str) -> asciichem::resolver::Cache {
    use std::sync::atomic::{AtomicU32, Ordering};
    static SEQUENCE: AtomicU32 = AtomicU32::new(0);
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|since| since.as_nanos())
        .unwrap_or_default();
    let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    asciichem::resolver::Cache::at_dir(std::env::temp_dir().join(format!(
        "asciichem-corpus-resolver-{tag}-{unique}-{sequence}"
    )))
}

fn sha256_hex(text: &str) -> String {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(text.as_bytes());
    let mut hex = String::with_capacity(64);
    for byte in digest {
        use std::fmt::Write as _;
        let _ = write!(hex, "{byte:02x}");
    }
    hex
}

#[test]
fn resolver_replay_contract() {
    use asciichem::resolver::adapter::Adapter;
    use asciichem::resolver::cache::DEFAULT_TTL_SECONDS;

    let Some(cases) = load_family("resolver.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 6,
        "expected the resolver family, got {} cases",
        cases.len()
    );

    // The Common Chemistry fan-out never consults the env/file key
    // chain: the shipped slot is replaced with the keyless adapter (the
    // corpus's key-absent case demands the no-key posture).
    let mut hermetic_registry = asciichem::resolver::Registry::shipped();
    hermetic_registry.register("common_chemistry", || {
        Box::new(asciichem::resolver::common_chemistry::CommonChemistry::keyless())
    });

    let mut failures = Vec::new();
    let mut ran = 0;
    let mut skipped: Vec<String> = Vec::new();
    for case in &cases {
        let id = case_id(case);
        if is_pending(case) {
            eprintln!("SKIP (pending for {IMPLEMENTATION}): {id}");
            skipped.push(id);
            continue;
        }
        ran += 1;
        let source = case
            .get("source")
            .and_then(serde_json::Value::as_str)
            .expect("resolver case source");
        let convention = case
            .get("convention")
            .and_then(serde_json::Value::as_str)
            .expect("resolver case convention");
        let value = case
            .get("value")
            .and_then(serde_json::Value::as_str)
            .expect("resolver case value");
        let expectation = case.get("expectation").expect("resolver expectation");
        let outcome = expectation
            .get("outcome")
            .and_then(serde_json::Value::as_str)
            .expect("resolver outcome");
        let cache = temp_corpus_cache("replay");
        let fetch = ReplayFetch::new(case.get("response").expect("recorded response"));

        if outcome == "skipped" {
            // A source-level failure drops the source without
            // poisoning the fan-out: the 401 posture answers empty.
            let results = hermetic_registry
                .resolve_all(
                    value,
                    convention,
                    Some(&[source]),
                    &fetch,
                    &cache,
                    false,
                    Some(PINNED_RETRIEVED_AT),
                )
                .expect("the failure isolates");
            if !results.is_empty() {
                failures.push(format!(
                    "{id}: the skipped source must drop without poisoning the others"
                ));
            }
        } else {
            let substance = replay_adapter(case)
                .resolve(
                    value,
                    convention,
                    &fetch,
                    &cache,
                    false,
                    Some(PINNED_RETRIEVED_AT),
                )
                .expect("replay resolves")
                .unwrap_or_else(|| {
                    panic!("{id}: the source must answer");
                });
            let wire = wire_record(&substance);
            if wire != *expectation.get("substance").expect("expected substance") {
                failures.push(format!(
                    "{id}: substance record mismatch\n  expected: {}\n  got: {}",
                    expectation.get("substance").expect("expected substance"),
                    wire
                ));
            }
            let structure_parses = expectation
                .get("structureParses")
                .and_then(serde_json::Value::as_bool)
                .unwrap_or(false);
            if substance.structure.is_some() != structure_parses {
                failures.push(format!(
                    "{id}: structureParses {} but structure is {}",
                    structure_parses,
                    substance.structure.is_some()
                ));
            }
        }

        // The adapter must build exactly the recorded request.
        let calls = fetch.calls.borrow();
        if calls.len() != 1 {
            failures.push(format!(
                "{id}: expected exactly one request, got {}",
                calls.len()
            ));
            continue;
        }
        let url = &calls[0].0;
        let expected_url = case
            .get("request")
            .and_then(|request| request.get("url"))
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        if url != expected_url {
            failures.push(format!(
                "{id}: adapter URL\n  expected: {expected_url}\n  got: {url}"
            ));
        }
        if captured_headers(&calls[0]) != expected_headers(case) {
            failures.push(format!(
                "{id}: request headers\n  expected: {:?}\n  got: {:?}",
                expected_headers(case),
                calls[0].1
            ));
        }

        // The cache markers pin the shared caching contract.
        let cache_contract = expectation.get("cache").expect("cache contract");
        let ttl = cache_contract
            .get("ttlSeconds")
            .and_then(serde_json::Value::as_u64)
            .expect("ttlSeconds");
        if DEFAULT_TTL_SECONDS != ttl {
            failures.push(format!(
                "{id}: ttlSeconds {} != the cache default {DEFAULT_TTL_SECONDS}",
                ttl
            ));
        }
        let key = cache.key_for(source, convention, value);
        let expected_key = sha256_hex(&format!("{source}:{convention}:{value}"));
        if key != expected_key {
            failures.push(format!(
                "{id}: cache key mismatch (sha256 over source:convention:value)"
            ));
        }

        if !cache_contract
            .get("servesOffline")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false)
        {
            continue;
        }
        // A seeded entry answers without re-fetching — and without a
        // key, for key-gated sources.
        if outcome == "skipped" {
            let seed_cache = temp_corpus_cache("seed");
            let seed_body = cases
                .iter()
                .find(|candidate| case_id(candidate) == "resolver/common-chemistry-cas-aspirin")
                .and_then(|candidate| candidate.get("response"))
                .and_then(|response| response.get("body"))
                .expect("the keyed recorded body for seeding");
            let seed_response = serde_json::json!({ "status": 200, "body": seed_body });
            let seed_fetch = ReplayFetch::new(&seed_response);
            let seeded = asciichem::resolver::common_chemistry::CommonChemistry::with_key(FAKE_KEY)
                .resolve(
                    value,
                    convention,
                    &seed_fetch,
                    &seed_cache,
                    false,
                    Some(PINNED_RETRIEVED_AT),
                )
                .expect("seed resolves")
                .expect("seed answers");
            let served = asciichem::resolver::common_chemistry::CommonChemistry::keyless()
                .resolve(value, convention, &OfflineFetch, &seed_cache, false, None)
                .expect("served offline");
            match served {
                Some(served) if wire_record(&served) == wire_record(&seeded) => {}
                Some(_) => failures.push(format!("{id}: the seeded entry answers differently")),
                None => failures.push(format!("{id}: the seeded entry must serve offline")),
            }
        } else {
            let served = replay_adapter(case)
                .resolve(value, convention, &OfflineFetch, &cache, false, None)
                .expect("served offline");
            match served {
                Some(served)
                    if wire_record(&served)
                        == *expectation.get("substance").expect("expected substance") => {}
                Some(_) => failures.push(format!("{id}: the cached record serves differently")),
                None => failures.push(format!("{id}: the cached record must serve offline")),
            }
        }
    }
    eprintln!(
        "resolver: {ran} run, {} skipped (pending: {skipped:?})",
        skipped.len()
    );
    assert!(
        failures.is_empty(),
        "{} resolver mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// -- citation family ------------------------------------------------------
// For the fixture's substance record (the wire form, retrieval time
// pinned) the citation layer must emit exactly the bibitem bytes —
// byte-equal, no trailing newline, the same captures the Ruby
// reference recorded. A rejected substance (no provenance) is not
// citable: resolve it first.

/// The bibitem's child-element kinds, in first-occurrence order (the
/// corpus pins the documented emission order).
fn child_element_order(xml: &str) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    let mut depth: i32 = 0;
    let mut rest = xml;
    while let Some(start) = rest.find('<') {
        rest = &rest[start + 1..];
        if rest.starts_with('/') {
            depth -= 1;
        } else if rest.starts_with('!') || rest.starts_with('?') {
            // declarations and processing instructions: not elements
        } else {
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_lowercase() || *c == '-')
                .collect();
            if depth == 1 && !names.contains(&name) {
                names.push(name.clone());
            }
            let self_closing = rest
                .find('>')
                .map(|end| rest[..end].ends_with('/'))
                .unwrap_or(false);
            if !self_closing {
                depth += 1;
            }
        }
    }
    names
}

#[test]
fn citation_golden_contract() {
    let Some(cases) = load_family("citation.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 4,
        "expected the citation family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut ran = 0;
    let mut skipped: Vec<String> = Vec::new();
    for case in &cases {
        let id = case_id(case);
        if is_pending(case) {
            eprintln!("SKIP (pending for {IMPLEMENTATION}): {id}");
            skipped.push(id);
            continue;
        }
        ran += 1;
        let substance = asciichem::resolver::Substance::from_model_json(
            &case.get("substance").expect("wire substance").to_string(),
        )
        .unwrap_or_else(|error| panic!("{id}: the wire substance must ingest: {error}"));
        if case
            .get("rejected")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false)
        {
            if asciichem::citation::bibitem(&substance)
                != Err(asciichem::citation::Error::NoProvenance)
            {
                failures.push(format!(
                    "{id}: a substance without provenance is not citable"
                ));
            }
            continue;
        }
        let golden = case
            .get("bibitem")
            .and_then(serde_json::Value::as_str)
            .expect("bibitem golden");
        match asciichem::citation::to_xml(&substance) {
            Ok(xml) if xml == golden => {
                let order = child_element_order(&xml);
                let expected_order = string_list(case, "ordering");
                if order != expected_order {
                    failures.push(format!("{id}: child order {order:?} != {expected_order:?}"));
                }
            }
            Ok(xml) => failures.push(format!(
                "{id}: bibitem bytes differ from the golden (byte-equal, no trailing newline)\n\
                 --- golden\n{golden}\n--- actual\n{xml}"
            )),
            Err(error) => failures.push(format!("{id}: to_xml failed: {error}")),
        }
    }
    eprintln!(
        "citation: {ran} run, {} skipped (pending: {skipped:?})",
        skipped.len()
    );
    assert!(
        failures.is_empty(),
        "{} citation mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// -- cli family -----------------------------------------------------------
// Each case runs `invocation` (the argv after the program name)
// against the in-process CLI dispatch and asserts exit plus the
// stdout/stderr matchers. stdout matchers are exact ("{version}"
// standing for this implementation's own version), substring (every
// entry appears), or semantic JSON (parsed comparison — key order
// cannot fail a case). Error text rides stderr. Invoked in-process:
// no shell, the offline seams injected so a stray network or cache
// access is a bug, not a side effect.

struct CapturedCli {
    code: i32,
    stdout: String,
    stderr: String,
}

/// Runs the CLI dispatch in-process with the offline deps.
fn run_cli(argv: &[&str]) -> CapturedCli {
    struct NoNetwork;

    impl asciichem::resolver::Fetch for NoNetwork {
        fn get(
            &self,
            url: &str,
            _headers: &[(&str, &str)],
        ) -> Result<Option<String>, asciichem::resolver::ResolverError> {
            Err(asciichem::resolver::ResolverError::Source(format!(
                "network must not be touched: {url}"
            )))
        }
    }

    let fetch = NoNetwork;
    let cache = temp_corpus_cache("cli");
    let registry = asciichem::resolver::Registry::shipped();
    let deps = asciichem::cli::Deps {
        fetch: &fetch,
        cache: &cache,
        registry: &registry,
    };
    let mut out = Vec::<u8>::new();
    let mut err = Vec::<u8>::new();
    // The dispatch parses clap-style: argv[0] is the program name (the
    // binary passes env::args_os(), which carries it) — the fixture's
    // invocation is the argv AFTER the program name, so it is restored
    // here, exactly as the Ruby runner's Cli.start does.
    let arguments = std::iter::once("asciichem").chain(argv.iter().copied());
    let code = asciichem::cli::run(arguments, &deps, &mut out, &mut err);
    CapturedCli {
        code,
        stdout: String::from_utf8_lossy(&out).into_owned(),
        stderr: String::from_utf8_lossy(&err).into_owned(),
    }
}

/// Runs one cli case and records any mismatch into `failures` — the
/// shared matcher both the runner contract and the pending-gate test
/// use (exit, then the stdout matchers, then stderr fragments).
fn run_cli_case(case: &serde_json::Value, id: &str, failures: &mut Vec<String>) {
    let invocation: Vec<String> = case
        .get("invocation")
        .and_then(serde_json::Value::as_array)
        .map(|list| {
            list.iter()
                .filter_map(serde_json::Value::as_str)
                .map(str::to_string)
                .collect()
        })
        .expect("cli case invocation");
    let expected_exit = case
        .get("exit")
        .and_then(serde_json::Value::as_i64)
        .expect("cli case exit code");
    let result = run_cli(&invocation.iter().map(String::as_str).collect::<Vec<_>>());
    let detail = format!(
        "exit {} (stdout: {:?}, stderr: {:?})",
        result.code,
        result.stdout.chars().take(120).collect::<String>(),
        result.stderr.chars().take(120).collect::<String>(),
    );
    if result.code as i64 != expected_exit {
        failures.push(format!("{id}: expected exit {expected_exit}; {detail}"));
        return;
    }
    if let Some(stdout) = case.get("stdout").and_then(serde_json::Value::as_str) {
        let expected = stdout.replace("{version}", env!("CARGO_PKG_VERSION"));
        if result.stdout != expected {
            failures.push(format!("{id}: stdout {:?} != {expected:?}", result.stdout));
        }
    } else if let Some(fragments) = case.get("stdoutIncludes") {
        for fragment in fragments.as_array().expect("stdoutIncludes list") {
            let fragment = fragment.as_str().expect("stdoutIncludes fragment");
            if !result.stdout.contains(fragment) {
                failures.push(format!(
                    "{id}: stdout {:?} does not include {fragment:?}",
                    result.stdout
                ));
            }
        }
    } else if let Some(expected) = case.get("stdoutJson") {
        let actual: serde_json::Value = match serde_json::from_str(&result.stdout) {
            Ok(parsed) => parsed,
            Err(error) => {
                failures.push(format!(
                    "{id}: stdout is not JSON ({error}): {:?}",
                    result.stdout
                ));
                return;
            }
        };
        if actual != *expected {
            failures.push(format!(
                "{id}: stdout JSON mismatch\n  expected: {expected}\n  got: {actual}"
            ));
        }
    }
    for fragment in string_list(case, "stderrIncludes") {
        if !result.stderr.contains(&fragment) {
            failures.push(format!(
                "{id}: stderr {:?} does not include {fragment:?}",
                result.stderr
            ));
        }
    }
}

#[test]
fn cli_contract() {
    let Some(cases) = load_family("cli.json") else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    assert!(
        cases.len() >= 20,
        "expected the cli family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut ran = 0;
    let mut skipped: Vec<String> = Vec::new();
    for case in &cases {
        let id = case_id(case);
        if is_pending(case) {
            eprintln!(
                "SKIP (pending for {IMPLEMENTATION}): {id} — {}",
                pending_for(case).join(", ")
            );
            skipped.push(id);
            continue;
        }
        ran += 1;
        run_cli_case(case, &id, &mut failures);
    }
    eprintln!(
        "cli: {ran} run, {} skipped (pending: {skipped:?})",
        skipped.len()
    );
    assert!(
        failures.is_empty(),
        "{} cli mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// -- rinchi family (the reaction-identity oracle) ------------------------
// Each case carries the recorded input edge (RXNfile or RDfile), the
// five oracle lines (identifier, RAuxInfo, the three RInChIKey shapes;
// keys recorded headerless), and the decoded reaction's per-role
// components in wire order — sourced from the inchi-tests RInChI
// tranche's recorded bytes. Derive: the input edge assembles the five
// lines byte-identically (equilibrium cases ride the forced reading —
// the model API can express it, unlike the string-only transports).
// Decode: the identifier plus its RAuxInfo companion reads back into
// the role-resolved components in order, and decode∘derive reproduces
// the oracle bytes — the components' aux/key content rides the
// re-assembled lines. The family rides the corpus pin's next move: an
// absent rinchi.json skips and reports (the pin predates the family),
// the same visible-gap posture as a fresh clone.

/// Strips a key line's header for the corpus's headerless comparison.
fn strip_key_header(line: &str, header: &str) -> String {
    line.strip_prefix(header)
        .map(str::to_owned)
        .unwrap_or_else(|| line.to_owned())
}

#[test]
fn rinchi_family_contract() {
    let Some(dir) = corpus_dir() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    let path = dir.join("rinchi.json");
    let Ok(text) = fs::read_to_string(&path) else {
        eprintln!(
            "rinchi.json not in the corpus pin yet — the family activates at the next pin bump; skipping"
        );
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("valid fixture JSON");
    assert!(
        cases.len() >= 36,
        "expected the rinchi family, got {} cases",
        cases.len()
    );

    let mut failures = Vec::new();
    let mut ran = 0;
    for case in &cases {
        let id = case_id(case);
        let forced = case
            .get("forceEquilibrium")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false);

        // The derive direction: input edge -> reaction -> the five
        // oracle lines.
        let read = if let Some(rxnfile) = case.get("rxnfile").and_then(serde_json::Value::as_str) {
            asciichem::rinchi::reaction_from_rxnfile(rxnfile)
        } else if let Some(rdfile) = case.get("rdfile").and_then(serde_json::Value::as_str) {
            asciichem::rinchi::reaction_from_rdfile(rdfile)
        } else {
            failures.push(format!("{id}: no input edge (rxnfile/rdfile)"));
            continue;
        };
        match read {
            Ok(mut reaction) => {
                if forced {
                    reaction.set_directionality(asciichem::rinchi::Directionality::Equilibrium);
                }
                match asciichem::rinchi::derive(&reaction) {
                    Ok(record) => {
                        ran += 1;
                        let checks = [
                            ("rinchi", record.rinchi, case.get("rinchi")),
                            ("rauxinfo", record.rauxinfo, case.get("rauxinfo")),
                            (
                                "rinchikey",
                                strip_key_header(&record.rinchikey, "Short-RInChIKey="),
                                case.get("rinchikey"),
                            ),
                            (
                                "longKey",
                                strip_key_header(&record.long_key, "Long-RInChIKey="),
                                case.get("longKey"),
                            ),
                            (
                                "webKey",
                                strip_key_header(&record.web_key, "Web-RInChIKey="),
                                case.get("webKey"),
                            ),
                        ];
                        for (what, actual, expected) in checks {
                            let expected = expected
                                .and_then(serde_json::Value::as_str)
                                .unwrap_or_default();
                            if actual != expected {
                                failures.push(format!(
                                    "{id}: derive {what}\n  expected: {expected}\n  actual:   {actual}"
                                ));
                            }
                        }
                    }
                    Err(error) => failures.push(format!("{id}: derive failed: {error}")),
                }
            }
            Err(error) => failures.push(format!("{id}: input edge read failed: {error}")),
        }

        // The decode direction: identifier (+ companion) -> roles in
        // wire order -> the round trip.
        let rinchi = case
            .get("rinchi")
            .and_then(serde_json::Value::as_str)
            .expect("rinchi line");
        let rauxinfo = case
            .get("rauxinfo")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        match asciichem::rinchi::decode_with_auxinfo(rinchi, rauxinfo) {
            Ok(decoded) => {
                ran += 1;
                // The direction flag recovers the real roles for the
                // directional readings — roles compare in wire order.
                // An equilibrium recording is symmetric: either group
                // may output first, so its components compare
                // order-free (the flag/role identity itself is pinned
                // by the round trip below).
                if forced {
                    let mut expected: Vec<String> = ["reactants", "products", "agents"]
                        .into_iter()
                        .flat_map(|role| {
                            case.get("components")
                                .and_then(|components| components.get(role))
                                .and_then(serde_json::Value::as_array)
                                .map(|list| {
                                    list.iter()
                                        .filter_map(serde_json::Value::as_str)
                                        .map(str::to_owned)
                                        .collect::<Vec<_>>()
                                })
                                .unwrap_or_default()
                        })
                        .collect();
                    expected.sort();
                    let mut actual: Vec<String> = decoded
                        .reactants()
                        .iter()
                        .chain(decoded.products().iter())
                        .chain(decoded.agents().iter())
                        .map(|component| component.inchi_string().expect("component inchi"))
                        .collect();
                    actual.sort();
                    if actual != expected {
                        failures.push(format!(
                            "{id}: decoded components\n  expected: {expected:?}\n  actual:   {actual:?}"
                        ));
                    }
                } else {
                    for (role, components) in [
                        ("reactants", decoded.reactants()),
                        ("products", decoded.products()),
                        ("agents", decoded.agents()),
                    ] {
                        let expected: Vec<String> = case
                            .get("components")
                            .and_then(|components| components.get(role))
                            .and_then(serde_json::Value::as_array)
                            .map(|list| {
                                list.iter()
                                    .filter_map(serde_json::Value::as_str)
                                    .map(str::to_owned)
                                    .collect()
                            })
                            .unwrap_or_default();
                        let actual: Vec<String> = components
                            .iter()
                            .map(|component| component.inchi_string().expect("component inchi"))
                            .collect();
                        if actual != expected {
                            failures.push(format!(
                                "{id}: decoded {role}\n  expected: {expected:?}\n  actual:   {actual:?}"
                            ));
                        }
                    }
                }
                let expected_directionality = if forced {
                    asciichem::rinchi::Directionality::Equilibrium
                } else {
                    asciichem::rinchi::Directionality::Directional
                };
                if decoded.directionality() != expected_directionality {
                    failures.push(format!(
                        "{id}: decoded directionality != {:?}",
                        expected_directionality
                    ));
                }
                // decode∘derive is the identity: the re-assembled
                // lines carry every component's aux/key content.
                match asciichem::rinchi::derive(&decoded) {
                    Ok(record) => {
                        for (what, actual, expected) in [
                            ("rinchi", record.rinchi, rinchi.to_owned()),
                            ("rauxinfo", record.rauxinfo, rauxinfo.to_owned()),
                            (
                                "rinchikey",
                                strip_key_header(&record.rinchikey, "Short-RInChIKey="),
                                case.get("rinchikey")
                                    .and_then(serde_json::Value::as_str)
                                    .unwrap_or_default()
                                    .to_owned(),
                            ),
                            (
                                "longKey",
                                strip_key_header(&record.long_key, "Long-RInChIKey="),
                                case.get("longKey")
                                    .and_then(serde_json::Value::as_str)
                                    .unwrap_or_default()
                                    .to_owned(),
                            ),
                            (
                                "webKey",
                                strip_key_header(&record.web_key, "Web-RInChIKey="),
                                case.get("webKey")
                                    .and_then(serde_json::Value::as_str)
                                    .unwrap_or_default()
                                    .to_owned(),
                            ),
                        ] {
                            if actual != expected {
                                failures.push(format!(
                                    "{id}: decode re-encodes {what}\n  expected: {expected}\n  actual:   {actual}"
                                ));
                            }
                        }
                    }
                    Err(error) => failures.push(format!("{id}: re-assembly failed: {error}")),
                }
            }
            Err(error) => failures.push(format!("{id}: decode failed: {error}")),
        }
    }
    eprintln!("rinchi: {ran} comparisons over {} cases", cases.len());
    assert!(
        failures.is_empty(),
        "{} rinchi mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
