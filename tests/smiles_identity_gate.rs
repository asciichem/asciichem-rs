//! The local SMILES identity gate (asciichem-rs#45): every identity
//! fixture's `smiles` field derives the recorded 1.07.5-oracle InChI
//! through the FULL local path — the native smiles parser → the model
//! molecule → the molfile pipe → the chemomer engine. No new oracle,
//! no fixture edits. The recorded exclusions: the 2 tetrahedral cases
//! gate through the pipe since the geometry-aware wedge landed
//! (asciichem-rs#86 — the enantiomer pair rides the equality battery
//! in identity.rs) and the 2 directional-bond cases gate since the
//! E/Z fixtures re-recorded from drawn coordinates (inchi-tests
//! v0.3.43's re-record: the /b layer survives the oracle now) — the
//! remaining skips are the 2 transcription disagreements (the fixture
//! smiles differ from their own recorded molfiles).

use std::fs;
use std::path::PathBuf;

fn suite_dir() -> PathBuf {
    std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"))
}

#[test]
fn the_local_smiles_path_derives_the_recorded_identities() {
    let dir = suite_dir();
    let mut cases: Vec<serde_json::Value> = Vec::new();
    for name in [
        "corpus/fixtures/identity.json",
        "corpus/fixtures/identity-breadth.json",
        "corpus/fixtures/identity-metal.json",
    ] {
        let Ok(text) = fs::read_to_string(dir.join(name)) else {
            return; // the corpus checkout gates in CI
        };
        cases.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).expect("fixtures"));
    }
    if cases.is_empty() {
        return;
    }

    let mut ran = 0usize;
    let mut skipped = Vec::new();
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(smiles) = case["smiles"].as_str() else {
            continue;
        };
        if matches!(id, "identity/sulfuric-acid" | "identity/sulfate-2") {
            skipped.push(id.to_string());
            continue;
        }
        let recorded = case["inchi"].as_str().expect("recorded inchi");
        ran += 1;
        match asciichem::identity::derive_from_smiles(smiles) {
            Ok(identity) if identity.inchi == recorded => {}
            Ok(identity) => {
                failures.push(format!("{id}: {:?} != oracle {recorded:?}", identity.inchi));
            }
            Err(error) => {
                failures.push(format!("{id}: derive errored: {error}"));
            }
        }
    }
    eprintln!(
        "local smiles identity gate: {ran} cases vs the oracle, {} skips, {} failures",
        skipped.len(),
        failures.len()
    );
    for skip in &skipped {
        eprintln!("  skip: {skip}");
    }
    assert!(
        failures.is_empty(),
        "{} local smile identities diverged:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(ran > 1_000, "the fixture corpus shrank: {ran}");
}
