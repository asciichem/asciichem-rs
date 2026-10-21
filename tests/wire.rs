//! Wire-form contract (L0): the canonical JSON round-trips through
//! emission + ingestion with the model intact, field names match the
//! shared v1 schemas (captured from the reference implementation), and
//! out-of-subset mechanisms report instead of guessing.

use asciichem::wire::{self, WireNode};
use asciichem::ParseError;
use std::fs;
use std::path::PathBuf;

fn corpus_dir() -> Option<PathBuf> {
    let dir = std::env::var("ASCIICHEM_CORPUS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../asciichem-tests/corpus/fixtures"));
    dir.is_dir().then_some(dir)
}

#[test]
fn wire_round_trip_contract() {
    let Some(dir) = corpus_dir() else {
        eprintln!("corpus not cloned as a sibling (or set ASCIICHEM_CORPUS); skipping");
        return;
    };

    let mut pass = 0;
    let skipped: Vec<String> = Vec::new();
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
            if !case
                .get("roundTrip")
                .and_then(serde_json::Value::as_bool)
                .unwrap_or(false)
            {
                continue;
            }
            let Some(input) = case.get("input").and_then(|value| value.as_str()) else {
                continue;
            };
            let id = case
                .get("id")
                .and_then(|value| value.as_str())
                .unwrap_or("<unnamed>");
            let run = (|| -> std::result::Result<(), String> {
                let model = asciichem::parse_nodes(input).map_err(|e: ParseError| e.to_string())?;
                let json = wire::to_model_json(&model).map_err(|e| e.to_string())?;
                let restored = wire::from_model_json(&json).map_err(|e| e.to_string())?;
                if restored != model {
                    return Err(format!(
                        "model drift through wire:\n wire: {json}\n before: {:?}\n  after: {:?}",
                        model, restored
                    ));
                }
                let rendered: Vec<String> = restored.iter().map(|n| n.to_string()).collect();
                if rendered.join(" ") != input {
                    return Err(format!(
                        "text drift through wire: {:?} -> {:?}",
                        input,
                        rendered.join(" ")
                    ));
                }
                Ok(())
            })();
            match run {
                Ok(()) => pass += 1,
                Err(reason) => mismatches.push(format!("{id}: {reason}")),
            }
        }
    }

    assert!(
        mismatches.is_empty(),
        "{} wire round-trip failure(s):\n{}\n(skips: {})\n(pass: {pass})",
        mismatches.len(),
        mismatches.join("\n"),
        skipped.join("; ")
    );
    eprintln!("wire round-trip: {pass} pass");
}

#[test]
fn wire_field_names_match_the_v1_schema() {
    // Captured from the reference implementation's to_model_json —
    // camelCase keys, type tags, empty collections omitted.
    let model = asciichem::parse_nodes(":C.").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    let atom = &json["nodes"][0]["nodes"][0];
    assert_eq!(atom["type"], "atom");
    assert_eq!(atom["element"], "C");
    assert_eq!(atom["lonePairs"], 1);
    assert_eq!(atom["radicalElectrons"], 1);
    assert!(
        json["nodes"][0]["identifiers"].is_null(),
        "empty collections omitted"
    );

    let model = asciichem::parse_nodes("Mn^(VII)").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    assert_eq!(json["nodes"][0]["nodes"][0]["oxidationState"], "VII");
}

#[test]
fn wire_carries_stereo_coefficients_groups_conditions() {
    let model = asciichem::parse_nodes("(a)-C_6H_12O_6").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    assert_eq!(json["nodes"][0]["stereo"], "alpha");

    let model = asciichem::parse_nodes("2H_2O").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    assert_eq!(json["nodes"][0]["coefficient"], "2");

    let model = asciichem::parse_nodes("(OH)_2").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    let group = &json["nodes"][0]["nodes"][0];
    assert_eq!(group["type"], "group");
    assert_eq!(group["multiplicity"], "2");
    assert_eq!(group["bracket"], "paren");

    let model = asciichem::parse_nodes("A ->[cat: Pd] B").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    assert_eq!(json["nodes"][0]["conditions"]["above"], "cat: Pd");
}

#[test]
fn wire_ingests_reference_json() {
    // Reference-emitted wire JSON ingests to the same model the same
    // text parses to.
    let ruby_json = r#"{"type":"formula","nodes":[{"type":"molecule","nodes":[{"type":"atom","element":"H","subscript":"2"},{"type":"atom","element":"O"}]}]}"#;
    let restored = wire::from_model_json(ruby_json).unwrap();
    let parsed = asciichem::parse_nodes("H_2O").unwrap();
    assert_eq!(restored, parsed);
}

#[test]
fn mechanism_steps_carry_typed_reactions_flows_spectators() {
    let model =
        asciichem::parse_nodes("mechanism{\n  s1: A -> B | flow: A > B\n  spectator: Na^+\n}")
            .unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    let mech = &json["nodes"][0];
    assert_eq!(mech["type"], "mechanism");
    assert_eq!(mech["steps"][0]["label"], "s1");
    assert_eq!(mech["steps"][0]["reaction"]["type"], "reaction");
    assert_eq!(mech["steps"][0]["flows"][0]["source"], "A");
    assert_eq!(mech["steps"][0]["flows"][0]["target"], "B");
    assert_eq!(mech["spectators"][0]["nodes"][0]["element"], "Na");
    assert_eq!(mech["spectators"][0]["nodes"][0]["charge"], "+");

    // Steps without flows omit the key entirely.
    let model = asciichem::parse_nodes("mechanism{\n  s1: A -> B\n}").unwrap();
    let json: serde_json::Value =
        serde_json::from_str(&wire::to_model_json(&model).unwrap()).unwrap();
    assert!(json["nodes"][0]["steps"][0].get("flows").is_none());
    assert!(json["nodes"][0].get("spectators").is_none());
}

#[test]
fn cascade_steps_have_no_wire_form() {
    let model = asciichem::parse_nodes("mechanism{\n  s1: A -> B -> C\n}").unwrap();
    let error = wire::to_model_json(&model).unwrap_err();
    assert!(format!("{error}").contains("cascades have no wire form"));
}

#[test]
fn wire_root_must_be_a_formula() {
    let error = wire::from_model_json(r#"{"type":"atom","element":"H"}"#).unwrap_err();
    assert!(format!("{error}").contains("root must be a formula"));
}

#[test]
fn wire_type_names_cover_the_ingest_set() {
    let model = asciichem::parse_nodes("\"heat\"").unwrap();
    let json = wire::to_model_json(&model).unwrap();
    assert!(json.contains(r#""type":"text""#), "{json}");
    let model = asciichem::parse_nodes("`K_c = 1`").unwrap();
    let json = wire::to_model_json(&model).unwrap();
    assert!(json.contains(r#""type":"embedded-math""#), "{json}");
    let _ = WireNode::Text;
}
