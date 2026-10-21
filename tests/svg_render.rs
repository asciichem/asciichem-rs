//! SVG contract: byte parity with the TypeScript implementation's
//! structural and linear renderers (the committed reference snapshot
//! was generated from a built asciichem-ts dist), determinism across
//! runs, and graceful degradation (bond-less molecules render
//! typographically, unclosed rings annotate with `?`).

use std::fs;
use std::path::PathBuf;

#[test]
fn svg_matches_the_typescript_reference() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/svg-reference.jsonl");
    let text = fs::read_to_string(&path).expect("reference snapshot");
    let mut pass = 0;
    let mut mismatches = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let case: serde_json::Value = serde_json::from_str(line).expect("valid jsonl");
        let input = case["input"].as_str().expect("input");
        for (field, actual) in [
            ("svg", asciichem::svg::render(input).unwrap()),
            ("linear", asciichem::svg::render_linear(input).unwrap()),
        ] {
            let expected = case[field].as_str().expect(field);
            if actual != expected {
                mismatches.push(format!(
                    "{input:?} ({field}) mismatch\n--- ts\n{expected}\n--- rust\n{actual}"
                ));
            } else {
                pass += 1;
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} SVG mismatch(es):\n{}\n(pass: {pass})",
        mismatches.len(),
        mismatches.join("\n")
    );
    eprintln!("svg reference: {pass} byte-exact renders");
}

#[test]
fn svg_is_deterministic_across_runs() {
    for input in ["C1-C-C-C-C-C1", "C1CCC(CCCCC1)", "C12-C-C1-C2", "A-B -> B"] {
        let first = asciichem::svg::render(input).unwrap();
        let second = asciichem::svg::render(input).unwrap();
        assert_eq!(first, second, "{input:?} is not deterministic");
    }
}

#[test]
fn bondless_molecules_render_typographically() {
    let out = asciichem::svg::render("H_2O").unwrap();
    assert!(!out.contains("<line"), "no bond lines without bonds: {out}");
    assert!(out.contains("<tspan"), "{out}");
}

#[test]
fn unclosed_rings_are_annotated_not_silent() {
    let out = asciichem::svg::render("C1-C-C").unwrap();
    assert!(
        out.contains(">#</text>") || out.contains(">?</text>"),
        "{out}"
    );
}

#[test]
fn structural_molecules_draw_bond_lines_and_labels() {
    let out = asciichem::svg::render("C1-C-C-C-C-C1").unwrap();
    assert!(out.contains("<line"), "{out}");
    assert!(
        out.contains("font-family=\"Helvetica, Arial, sans-serif\""),
        "{out}"
    );
    let bonds = out.matches("<line").count();
    assert_eq!(bonds, 6, "cyclohexane closes all six bonds: {out}");
}
