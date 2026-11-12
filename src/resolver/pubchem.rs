//! PubChem source adapter (TODO.impl 39; TODO.v2 07 Layer 2) — the
//! Rust mirror of the Ruby `AsciiChem::Resolver::PubChem`. PubChem is
//! the default source: keyless PUG-REST, permissive terms, and the
//! bulk mirror of Common Chemistry's validated CAS RNs (CAS numbers
//! are indexed as names). Politeness: one request per resolution,
//! 5/10s timeouts, no retries — PubChem asks for max 5 requests per
//! second.
//!
//! Citation profile (TODO.v2 08): publisher NCBI, per-access
//! versioning (retrievedAt provenance), compound link
//! `https://pubchem.ncbi.nlm.nih.gov/compound/{cid}`.

use crate::resolver::adapter::{cgi_escape, parse_structure, Adapter, Fetch};
use crate::resolver::substance::{Identifier, Property, Substance};
use crate::resolver::ResolverError;

const BASE: &str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound";
const PROPERTIES: [&str; 7] = [
    "CanonicalSMILES",
    "IsomericSMILES",
    "MolecularFormula",
    "MolecularWeight",
    "IUPACName",
    "InChI",
    "InChIKey",
];

/// The PubChem source.
#[derive(Debug, Clone, Copy, Default)]
pub struct PubChem;

impl Adapter for PubChem {
    fn source_name(&self) -> &'static str {
        "pubchem"
    }

    fn supported_conventions(&self) -> &'static [&'static str] {
        &[
            "name",
            "cas",
            "pubchem-cid",
            "inchikey",
            "inchi",
            "smiles",
            "canonical-smiles",
        ]
    }

    fn attribution(&self) -> Option<String> {
        Some("PubChem, U.S. National Library of Medicine".to_string())
    }

    fn fetch_record(
        &self,
        value: &str,
        convention: &str,
        fetch: &dyn Fetch,
    ) -> Result<Option<Substance>, ResolverError> {
        let url = url_for(value, convention)?;
        let Some(body) = fetch.get(&url, &[])? else {
            return Ok(None);
        };
        let payload: serde_json::Value = serde_json::from_str(&body).map_err(|error| {
            ResolverError::Source(format!("pubchem returned malformed JSON: {error}"))
        })?;
        let properties = payload
            .get("PropertyTable")
            .and_then(|table| table.get("Properties"))
            .and_then(|list| list.get(0));
        let Some(properties) = properties else {
            return Ok(None);
        };
        Ok(substance_from(self, properties))
    }
}

/// CAS numbers are indexed as PubChem names (Common Chemistry
/// mirror), so both map to the name namespace.
fn url_for(value: &str, convention: &str) -> Result<String, ResolverError> {
    let namespace = match convention {
        "name" | "cas" => "name",
        "pubchem-cid" => "cid",
        "inchikey" => "inchikey",
        "inchi" => "inchi",
        "smiles" | "canonical-smiles" => "smiles",
        other => {
            return Err(ResolverError::Source(format!(
                "pubchem does not answer convention {other:?}"
            )))
        }
    };
    Ok(format!(
        "{BASE}/{namespace}/{}/property/{}/JSON",
        cgi_escape(value),
        PROPERTIES.join(",")
    ))
}

fn json_scalar(value: Option<&serde_json::Value>) -> Option<String> {
    match value? {
        serde_json::Value::String(text) => Some(text.clone()),
        scalar @ (serde_json::Value::Number(_) | serde_json::Value::Bool(_)) => {
            Some(scalar.to_string())
        }
        _ => None,
    }
}

/// The Ruby `String#to_f` posture across the JSON scalar forms the
/// PUG-REST answers with (the fixture's MolecularWeight rides a
/// string).
fn json_scalar_f64(value: Option<&serde_json::Value>) -> Option<f64> {
    match value? {
        serde_json::Value::Number(number) => number.as_f64(),
        serde_json::Value::String(text) => text.trim().parse::<f64>().ok(),
        _ => None,
    }
}

fn substance_from(adapter: &PubChem, properties: &serde_json::Value) -> Option<Substance> {
    let provenance = adapter.default_provenance();
    let mut identifiers = Vec::new();
    if let Some(cid) = json_scalar(properties.get("CID")) {
        identifiers.push(Identifier {
            value: cid,
            convention: "pubchem-cid".to_string(),
            provenance: Some(provenance.clone()),
        });
    }
    for (key, convention) in [
        ("InChIKey", "inchikey"),
        ("InChI", "inchi"),
        ("CanonicalSMILES", "canonical-smiles"),
    ] {
        if let Some(value) = json_scalar(properties.get(key)) {
            identifiers.push(Identifier {
                value,
                convention: convention.to_string(),
                provenance: Some(provenance.clone()),
            });
        }
    }

    let formula = json_scalar(properties.get("MolecularFormula"));
    let molecular_weight_string = json_scalar(properties.get("MolecularWeight"));
    let mut wire_properties = Vec::new();
    if let Some(value) = &formula {
        wire_properties.push(Property {
            name: "molecular-formula".to_string(),
            value: value.clone(),
            units: None,
            provenance: Some(provenance.clone()),
        });
    }
    if let Some(value) = &molecular_weight_string {
        wire_properties.push(Property {
            name: "molecular-weight".to_string(),
            value: value.clone(),
            units: Some("g/mol".to_string()),
            provenance: Some(provenance.clone()),
        });
    }

    let smiles = json_scalar(properties.get("CanonicalSMILES"));
    Some(Substance {
        preferred_name: json_scalar(properties.get("IUPACName")),
        synonyms: Vec::new(),
        formula,
        molecular_weight: json_scalar_f64(properties.get("MolecularWeight")),
        identifiers,
        properties: wire_properties,
        // Best effort: PubChem canonical SMILES may use constructs
        // outside the supported subset (stereo); then the SMILES
        // stays as an identifier and structure is None — no silent loss.
        structure: parse_structure(smiles.as_deref()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Node;
    use crate::resolver::adapter::testing::{fixture, CaptureFetch, RecordedFetch};
    use crate::resolver::cache::temp_cache;
    use crate::resolver::Registry;
    use crate::smiles::write_smiles;

    fn aspirin_fixture() -> String {
        fixture("pubchem/aspirin.json")
    }

    #[test]
    fn registers_itself_as_the_pubchem_source() {
        let registry = Registry::shipped();
        assert_eq!(
            registry.source_names(),
            vec!["pubchem", "opsin", "common_chemistry"]
        );
        let adapter = registry.adapter("pubchem").expect("registered");
        assert_eq!(adapter.source_name(), "pubchem");
    }

    #[test]
    fn supports_the_documented_identifier_conventions() {
        let adapter = PubChem;
        for convention in [
            "name",
            "cas",
            "inchikey",
            "pubchem-cid",
            "inchi",
            "smiles",
            "canonical-smiles",
        ] {
            assert!(adapter.supports(convention), "{convention}");
        }
    }

    #[test]
    fn resolves_aspirin_by_name_from_the_recorded_response() {
        let adapter = PubChem;
        let substance = adapter
            .resolve(
                "aspirin",
                "name",
                &RecordedFetch {
                    body: Some(aspirin_fixture()),
                },
                &temp_cache("pubchem"),
                false,
            )
            .expect("resolves")
            .expect("answers");
        assert_eq!(
            substance.preferred_name.as_deref(),
            Some("2-acetyloxybenzoic acid")
        );
        assert_eq!(substance.formula.as_deref(), Some("C9H8O4"));
        assert_eq!(substance.molecular_weight, Some(180.16));
        assert_eq!(
            substance.identifier_value("inchikey"),
            Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        );
        assert_eq!(substance.identifier_value("pubchem-cid"), Some("2244"));
    }

    #[test]
    fn carries_provenance_with_source_and_attribution_on_every_identifier() {
        let adapter = PubChem;
        let substance = adapter
            .resolve(
                "aspirin",
                "name",
                &RecordedFetch {
                    body: Some(aspirin_fixture()),
                },
                &temp_cache("pubchem-prov"),
                false,
            )
            .expect("resolves")
            .expect("answers");
        let provenance = substance.provenance().expect("carries provenance");
        assert_eq!(provenance.source, "pubchem");
        assert_eq!(
            provenance.attribution.as_deref(),
            Some("PubChem, U.S. National Library of Medicine")
        );
        let retrieved_at = provenance.retrieved_at.as_deref().expect("stamped");
        assert!(
            retrieved_at.len() >= 11 && retrieved_at.chars().nth(10) == Some('T'),
            "{retrieved_at}"
        );
        for identifier in &substance.identifiers {
            assert_eq!(identifier.provenance.as_ref().unwrap().source, "pubchem");
        }
    }

    #[test]
    fn parses_the_structure_when_the_smiles_is_in_the_supported_subset() {
        let adapter = PubChem;
        let substance = adapter
            .resolve(
                "aspirin",
                "name",
                &RecordedFetch {
                    body: Some(aspirin_fixture()),
                },
                &temp_cache("pubchem-struct"),
                false,
            )
            .expect("resolves")
            .expect("answers");
        let structure = substance.structure.as_ref().expect("parses");
        let smiles = write_smiles(&[Node::Molecule(structure.clone())]).expect("writes");
        assert_eq!(smiles, "CC(=O)Oc1ccccc1C(=O)O");
    }

    #[test]
    fn parses_stereo_smiles_into_structures_chirality_left_intact() {
        let stereo = aspirin_fixture().replace("CC(=O)OC1=CC=CC=C1C(=O)O", "C[C@H](N)C(=O)O");
        let adapter = PubChem;
        let substance = adapter
            .resolve(
                "aspirin",
                "name",
                &RecordedFetch { body: Some(stereo) },
                &temp_cache("pubchem-stereo"),
                false,
            )
            .expect("resolves")
            .expect("answers");
        let structure = substance.structure.as_ref().expect("parses");
        let chiral = structure.nodes.iter().find_map(|node| match node {
            Node::Atom(atom) => atom
                .chirality
                .as_deref()
                .filter(|chirality| !chirality.is_empty())
                .map(|c| c.to_string()),
            _ => None,
        });
        assert_eq!(chiral.as_deref(), Some("TH1"));
    }

    #[test]
    fn keeps_the_smiles_as_an_identifier_when_it_is_outside_the_subset() {
        let wildcard =
            aspirin_fixture().replace("CC(=O)OC1=CC=CC=C1C(=O)O", "[*]C(=O)OC1=CC=CC=C1C(=O)O");
        let adapter = PubChem;
        let substance = adapter
            .resolve(
                "aspirin",
                "name",
                &RecordedFetch {
                    body: Some(wildcard),
                },
                &temp_cache("pubchem-wild"),
                false,
            )
            .expect("resolves")
            .expect("answers");
        assert_eq!(substance.structure, None);
        assert!(substance
            .identifier_value("canonical-smiles")
            .expect("kept")
            .contains("[*]"));
    }

    #[test]
    fn returns_none_for_a_not_found_response() {
        let adapter = PubChem;
        let resolved = adapter
            .resolve(
                "nonsense-xyz",
                "name",
                &RecordedFetch { body: None },
                &temp_cache("pubchem-miss"),
                false,
            )
            .expect("not-found is not an error");
        assert_eq!(resolved, None);
    }

    #[test]
    fn builds_pug_rest_urls_per_convention_cas_via_the_name_namespace() {
        let adapter = PubChem;
        let capture = CaptureFetch::new();
        adapter
            .resolve(
                "50-78-2",
                "cas",
                &capture,
                &temp_cache("pubchem-url"),
                false,
            )
            .expect("not-found is not an error");
        let url = capture.calls().pop().expect("requested").url;
        assert!(url.contains("/compound/name/50-78-2/property/"), "{url}");
    }
}
