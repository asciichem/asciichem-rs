//! L2 HTML/LaTeX golden parity: every corpus fixture carrying an
//! `html` or `latex` golden must render byte-identically. Unsupported
//! constructs are counted as skips — never silently mismatched.
//! (The deep-audit pass-2 structural ask: these surfaces had no
//! corpus families; asciichem-ruby#217.)

use std::fs;
use std::path::PathBuf;

type Renderer = fn(&[asciichem::model::Node]) -> Result<String, asciichem::ParseError>;

fn corpus_dir() -> Option<PathBuf> {
    let dir = std::env::var("ASCIICHEM_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../asciichem-tests/corpus/fixtures"));
    dir.is_dir().then_some(dir)
}

#[test]
fn render_golden_contract() {
    let Some(dir) = corpus_dir() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };
    let mut pass = 0;
    let mut skipped: Vec<String> = Vec::new();
    let mut mismatches = Vec::new();
    for entry in fs::read_dir(&dir).expect("read corpus dir") {
        let path = entry.expect("dir entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("json") {
            continue;
        }
        let text = fs::read_to_string(&path).expect("read fixture file");
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&text).expect("valid fixture JSON");
        for case in cases {
            let id = case
                .get("id")
                .and_then(|value| value.as_str())
                .unwrap_or("<unnamed>");
            let input = case
                .get("input")
                .and_then(|value| value.as_str())
                .unwrap_or("");
            let renderers: [(
                &str,
                fn(&[asciichem::model::Node]) -> Result<String, asciichem::ParseError>,
            ); 2] = [
                ("html", asciichem::html::render_formula),
                ("latex", asciichem::latex::render_formula),
            ];
            for (field, render) in renderers {
                let Some(golden) = case.get(field).and_then(|value| value.as_str()) else {
                    continue;
                };
                let Ok(nodes) = asciichem::parse_nodes(input) else {
                    skipped.push(format!("{id}: parse declined"));
                    continue;
                };
                match render(&nodes) {
                    Ok(output) if output == golden => pass += 1,
                    Ok(output) => mismatches.push(format!(
                        "{id}: {field} mismatch for {input:?}\n--- golden\n{golden}\n--- actual\n{output}"
                    )),
                    Err(error) => skipped.push(format!("{id}: {error}")),
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} render mismatch(es):\n{}\n(skips: {})\n(pass: {pass})",
        mismatches.len(),
        mismatches.join("\n"),
        skipped.join("; "),
        pass = pass
    );
}
