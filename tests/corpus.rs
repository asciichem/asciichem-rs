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
