//! Identity-layer conformance over the chemomer oracle corpus
//! (chemomer/inchi-tests): every `identity*.json` case's recorded
//! molfile derives — through the IN-PROCESS chemomer dependency — the
//! recorded InChI, AuxInfo and InChIKey byte-identically. This is the
//! same corpus family and the same pin (v0.3.41) chemomer's own CI
//! gates, run here against the consumed crate rather than the engine
//! checkout. Absent corpus = the test reports a skip rather than
//! failing a fresh clone, the same hermetic shape as the
//! ASCIICHEM_CORPUS gates.

use std::fs;
use std::path::PathBuf;

fn corpus_dir() -> Option<PathBuf> {
    let dir = std::env::var("CHEMOMER_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests/corpus/fixtures"));
    dir.is_dir().then_some(dir)
}

struct Case {
    id: String,
    molfile: String,
    inchi: String,
    aux: Option<String>,
    inchikey: String,
}

fn load_cases() -> Option<Vec<Case>> {
    let dir = corpus_dir()?;
    let mut cases = Vec::new();
    let mut paths: Vec<PathBuf> = fs::read_dir(&dir)
        .expect("read corpus dir")
        .filter_map(|entry| {
            let path = entry.expect("dir entry").path();
            let name = path.file_name()?.to_str()?.to_string();
            name.starts_with("identity")
                .then_some(name.ends_with(".json"))
                .filter(|is_json| *is_json)
                .map(|_| path)
        })
        .collect();
    paths.sort();
    for path in paths {
        let text = fs::read_to_string(&path).expect("read fixture file");
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&text).expect("valid fixture JSON");
        for case in parsed {
            let (Some(molfile), Some(inchi), Some(inchikey)) = (
                case.get("molfile").and_then(serde_json::Value::as_str),
                case.get("inchi").and_then(serde_json::Value::as_str),
                case.get("inchikey").and_then(serde_json::Value::as_str),
            ) else {
                continue;
            };
            cases.push(Case {
                id: case
                    .get("id")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("<unnamed>")
                    .to_string(),
                molfile: molfile.to_string(),
                inchi: inchi.to_string(),
                aux: case
                    .get("aux")
                    .and_then(serde_json::Value::as_str)
                    .map(str::to_string),
                inchikey: inchikey.to_string(),
            });
        }
    }
    Some(cases)
}

#[test]
fn identity_family_is_byte_exact_in_process() {
    let Some(cases) = load_cases() else {
        eprintln!("chemomer corpus not cloned as a sibling (or set CHEMOMER_CORPUS); skipping");
        return;
    };

    let mut mismatches = Vec::new();
    for case in &cases {
        let derived = asciichem::identity::derive_from_molfile(&case.molfile)
            .unwrap_or_else(|error| panic!("{}: derivation failed: {error}", case.id));
        if derived.inchi != case.inchi {
            mismatches.push(format!(
                "{}: inchi\n  recorded: {}\n  derived:  {}",
                case.id, case.inchi, derived.inchi
            ));
        }
        if derived.inchikey != case.inchikey {
            mismatches.push(format!(
                "{}: inchikey\n  recorded: {}\n  derived:  {}",
                case.id, case.inchikey, derived.inchikey
            ));
        }
        if let Some(aux) = &case.aux {
            if &derived.aux != aux {
                mismatches.push(format!(
                    "{}: aux\n  recorded: {}\n  derived:  {}",
                    case.id, aux, derived.aux
                ));
            }
        }
    }

    assert!(
        cases.len() > 5000,
        "expected the full identity family, got {}",
        cases.len()
    );
    assert!(
        mismatches.is_empty(),
        "{} mismatch(es):\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!(
        "identity family: {} cases byte-exact in-process",
        cases.len()
    );
}
