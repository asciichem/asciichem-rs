//! Normalisation of parsanol's parslet-shaped `Value` tree into the
//! DSL shapes the model builder consumes — the Rust mirror of the
//! reference engine's TreeAdapter (asciichem-ruby
//! `engine/parg_engine.rb`): single-key wrapper captures unwrap into
//! their payload, `group` collapses with an explicit multiplicity,
//! `embedded_math` re-keys to `math_source`, and sequence captures
//! under the merge keys collapse to one hash (later duplicate keys
//! append, matching parslet's repeated-capture behavior).

use std::collections::HashMap;

use parsanol::portable::transform::Value;

const WRAPPERS: [&str; 9] = [
    "molecule",
    "reaction",
    "annotated_molecule",
    "segment",
    "plus_term",
    "annotation_item",
    "formula_tail",
    "condition",
    "embedded_math",
];

const ARRAY_MERGE_KEYS: [&str; 3] = ["cascade", "molecule", "annotated_molecule"];

pub fn normalise(value: &Value) -> Value {
    match value {
        Value::Array(items) => Value::Array(items.iter().map(normalise).collect()),
        Value::Hash(pairs) => normalise_hash(pairs),
        other => other.clone(),
    }
}

fn normalise_hash(pairs: &HashMap<String, Value>) -> Value {
    let adapted: HashMap<String, Value> = pairs
        .iter()
        .map(|(key, value)| (key.clone(), normalise(value)))
        .collect();

    if adapted.len() == 1 {
        let (key, inner) = adapted.into_iter().next().expect("len checked");
        if key == "embedded_math" {
            let source = match inner {
                Value::Hash(hash) => hash.get("math_source").cloned().unwrap_or(Value::Nil),
                other => other,
            };
            let mut out = HashMap::new();
            out.insert("math_source".to_string(), source);
            return Value::Hash(out);
        }
        if WRAPPERS.contains(&key.as_str()) {
            return inner;
        }
        if key == "group" {
            return merge_multiplicity(merge_items(&inner));
        }
        let mut out = HashMap::new();
        out.insert(key, inner);
        return Value::Hash(out);
    }

    let mut out = adapted;
    for key in ARRAY_MERGE_KEYS {
        if matches!(out.get(key), Some(Value::Array(_))) {
            let merged = merge_items(out.get(key).expect("checked"));
            out.insert(key.to_string(), merged);
        }
    }
    Value::Hash(out)
}

fn merge_multiplicity(mut merged: Value) -> Value {
    if let Value::Hash(hash) = &mut merged {
        hash.entry("multiplicity".to_string()).or_insert(Value::Nil);
    }
    merged
}

// Merges an array of capture hashes into one hash; later duplicate
// keys append to arrays (parslet repeated-capture behavior).
fn merge_items(value: &Value) -> Value {
    let items = match value {
        Value::Array(items) if items.iter().all(|i| matches!(i, Value::Hash(_))) => items,
        _ => return value.clone(),
    };
    let mut merged: HashMap<String, Value> = HashMap::new();
    for item in items {
        if let Value::Hash(pairs) = item {
            for (key, new) in pairs {
                match merged.get_mut(key) {
                    Some(Value::Array(old)) => {
                        let mut combined = old.clone();
                        match new {
                            Value::Array(add) => combined.extend(add.clone()),
                            other => combined.push(other.clone()),
                        }
                        merged.insert(key.clone(), Value::Array(combined));
                    }
                    _ => {
                        merged.insert(key.clone(), new.clone());
                    }
                }
            }
        }
    }
    Value::Hash(merged)
}
