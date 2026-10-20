//! Shared-corpus conformance: the parser acceptance/rejection contract.
//!
//! Every fixture case that carries `input` (and is not a linter or
//! convention case, which assert diagnostics rather than parsing) must
//! parse or reject exactly as marked. The fixtures live in the sibling
//! `asciichem-tests` checkout (set `ASCIICHEM_CORPUS` to override).
//! Absent corpus = the test reports a skip rather than failing a fresh
//! clone.

use std::fs;
use std::path::PathBuf;

fn corpus_dir() -> Option<PathBuf> {
    let dir = std::env::var("ASCIICHEM_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../asciichem-tests/corpus/fixtures"));
    dir.is_dir().then_some(dir)
}

#[test]
fn parser_acceptance_contract() {
    let Some(dir) = corpus_dir() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };

    let mut cases = Vec::new();
    for entry in fs::read_dir(&dir).expect("read corpus dir") {
        let path = entry.expect("dir entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
            continue;
        }
        let text = fs::read_to_string(&path).expect("read fixture file");
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&text).expect("valid fixture JSON");
        cases.extend(parsed);
    }

    let mut failures = Vec::new();
    let mut checked = 0;
    for case in &cases {
        let Some(input) = case.get("input").and_then(|value| value.as_str()) else {
            continue;
        };
        if case.get("lint").is_some() || case.get("convention").is_some() {
            continue;
        }
        let id = case
            .get("id")
            .and_then(|value| value.as_str())
            .unwrap_or("<unnamed>");
        let expected = case
            .get("parses")
            .and_then(serde_json::Value::as_bool)
            .unwrap_or(false);
        let actual = asciichem::parse(input).is_ok();
        checked += 1;
        if actual != expected {
            failures.push(format!(
                "{id}: expected parses={expected}, got {actual} ({:?})",
                input
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
fn holds_no_state_between_parses() {
    // The grammar singleton is immutable; parsing the same input twice
    // gives the same verdict, and a reject does not poison later calls.
    assert!(asciichem::parse("H_2O").is_ok());
    assert!(asciichem::parse("   ").is_err());
    assert!(asciichem::parse("H_2O").is_ok());
}
