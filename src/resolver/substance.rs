//! Domain value objects for resolution results — the Rust mirror of
//! the Ruby `Substance`/`Provenance`/`Identifier`/`Property` structs.
//! These are data carriers; (de)serialization to the canonical wire
//! form goes through the shared wire substance-record types
//! ([`crate::wire::SubstanceRecord`], serde derive — the house wire
//! pattern), never hand-rolled JSON.

use crate::wire::{self, IdentifierTag, ProvenanceTag, SubstanceRecordTag, WireError};

/// A retrieval's provenance: which source answered, when, at what
/// version, owing what attribution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Provenance {
    pub source: String,
    pub retrieved_at: Option<String>,
    pub source_version: Option<String>,
    pub attribution: Option<String>,
}

impl Provenance {
    /// The provenance record as one sentence (the citation track's
    /// bibitem note): retrieval source and version plus the source's
    /// required attribution — licence terms (CAS Common Chemistry is
    /// CC BY-NC 4.0) must survive into the rendered document. None
    /// when there is nothing beyond the source name itself.
    pub fn note_text(&self) -> Option<String> {
        let mut details = Vec::new();
        if let Some(version) = &self.source_version {
            details.push(format!("version {version}"));
        }
        if let Some(attribution) = &self.attribution {
            details.push(attribution.clone());
        }
        if details.is_empty() {
            return None;
        }
        Some(format!(
            "Retrieved from {}; {}",
            self.source,
            details.join("; ")
        ))
    }
}

/// One identifier a source answered with, carrying its own provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Identifier {
    pub value: String,
    pub convention: String,
    pub provenance: Option<Provenance>,
}

/// One observed property, carrying its own provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Property {
    pub name: String,
    pub value: String,
    pub units: Option<String>,
    pub provenance: Option<Provenance>,
}

/// A resolved substance: identity + provenance + best-effort
/// structure. `structure` is a model molecule when the source's SMILES
/// parses within the supported subset; otherwise None and the SMILES
/// stays available as an identifier — nothing is dropped. The
/// structure does not ride the wire form (it is derivable from the
/// SMILES identifier), exactly as in the Ruby reference.
#[derive(Debug, Clone, PartialEq)]
pub struct Substance {
    pub preferred_name: Option<String>,
    pub synonyms: Vec<String>,
    pub formula: Option<String>,
    pub molecular_weight: Option<f64>,
    pub identifiers: Vec<Identifier>,
    pub properties: Vec<Property>,
    pub structure: Option<crate::model::Molecule>,
}

impl Substance {
    /// The first identifier carrying `convention`, if any.
    pub fn identifier_value(&self, convention: &str) -> Option<&str> {
        self.identifiers
            .iter()
            .find(|identifier| identifier.convention == convention)
            .map(|identifier| identifier.value.as_str())
    }

    /// The substance-level provenance: the first identifier's.
    pub fn provenance(&self) -> Option<&Provenance> {
        self.identifiers.first().and_then(|i| i.provenance.as_ref())
    }

    /// The canonical wire form (asciichem-model v1 substance-record).
    pub fn to_wire(&self) -> wire::SubstanceRecord {
        wire::SubstanceRecord {
            tag: SubstanceRecordTag::SubstanceRecord,
            preferred_name: self.preferred_name.clone(),
            synonyms: self.synonyms.clone(),
            formula: self.formula.clone(),
            molecular_weight: self.molecular_weight,
            identifiers: self
                .identifiers
                .iter()
                .map(|identifier| wire::ProvenancedIdentifier {
                    identifier: wire::Identifier {
                        tag: IdentifierTag::Identifier,
                        value: identifier.value.clone(),
                        convention: identifier.convention.clone(),
                        dict_ref: None,
                    },
                    provenance: identifier.provenance.as_ref().map(provenance_wire),
                })
                .collect(),
            properties: self
                .properties
                .iter()
                .map(|property| wire::ProvenancedProperty {
                    name: property.name.clone(),
                    value: property.value.clone(),
                    units: property.units.clone(),
                    provenance: property.provenance.as_ref().map(provenance_wire),
                })
                .collect(),
        }
    }

    /// Rebuilds a Substance from the wire form (the cache path). The
    /// structure does not ride the wire — it stays None here, exactly
    /// as in the Ruby reference.
    pub fn from_wire(record: &wire::SubstanceRecord) -> Substance {
        Substance {
            preferred_name: record.preferred_name.clone(),
            synonyms: record.synonyms.clone(),
            formula: record.formula.clone(),
            molecular_weight: record.molecular_weight,
            identifiers: record
                .identifiers
                .iter()
                .map(|provenanced| Identifier {
                    value: provenanced.identifier.value.clone(),
                    convention: provenanced.identifier.convention.clone(),
                    provenance: provenance_domain(provenanced.provenance.as_ref()),
                })
                .collect(),
            properties: record
                .properties
                .iter()
                .map(|provenanced| Property {
                    name: provenanced.name.clone(),
                    value: provenanced.value.clone(),
                    units: provenanced.units.clone(),
                    provenance: provenance_domain(provenanced.provenance.as_ref()),
                })
                .collect(),
            structure: None,
        }
    }

    /// The canonical wire JSON (`to_model_json` in the Ruby
    /// reference).
    pub fn to_model_json(&self) -> Result<String, WireError> {
        serde_json::to_string(&self.to_wire())
            .map_err(|error| WireError::Malformed(error.to_string()))
    }

    /// Parses the canonical wire JSON (`from_model_json` in the Ruby
    /// reference).
    pub fn from_model_json(json: &str) -> Result<Substance, WireError> {
        let record: wire::SubstanceRecord =
            serde_json::from_str(json).map_err(|error| WireError::Malformed(error.to_string()))?;
        Ok(Substance::from_wire(&record))
    }
}

fn provenance_wire(provenance: &Provenance) -> wire::Provenance {
    wire::Provenance {
        tag: ProvenanceTag::Provenance,
        source: provenance.source.clone(),
        retrieved_at: provenance.retrieved_at.clone(),
        source_version: provenance.source_version.clone(),
        attribution: provenance.attribution.clone(),
    }
}

fn provenance_domain(provenance: Option<&wire::Provenance>) -> Option<Provenance> {
    provenance.map(|wire| Provenance {
        source: wire.source.clone(),
        retrieved_at: wire.retrieved_at.clone(),
        source_version: wire.source_version.clone(),
        attribution: wire.attribution.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wire_json() -> &'static str {
        // The shared v1 example's shape: camelCase wire keys, type
        // tags, provenance on every identifier.
        r#"{
            "type": "substance-record",
            "preferredName": "Aspirin",
            "synonyms": ["acetylsalicylic acid"],
            "formula": "C9H8O4",
            "molecularWeight": 180.16,
            "identifiers": [
                {
                    "identifier": {"type": "identifier", "value": "50-78-2", "convention": "cas"},
                    "provenance": {"type": "provenance", "source": "common_chemistry",
                                   "retrievedAt": "2026-09-09T00:00:00Z",
                                   "attribution": "CAS Common Chemistry (CC BY-NC 4.0)"}
                }
            ],
            "properties": [
                {"name": "molecular-formula", "value": "C9H8O4",
                 "provenance": {"type": "provenance", "source": "common_chemistry"}}
            ]
        }"#
    }

    #[test]
    fn ingests_the_canonical_wire_form() {
        let substance = Substance::from_model_json(wire_json()).expect("ingests");
        assert_eq!(substance.preferred_name.as_deref(), Some("Aspirin"));
        assert_eq!(substance.synonyms, vec!["acetylsalicylic acid"]);
        assert_eq!(substance.molecular_weight, Some(180.16));
        assert_eq!(substance.identifier_value("cas"), Some("50-78-2"));
        let provenance = substance.provenance().expect("carries provenance");
        assert_eq!(provenance.source, "common_chemistry");
        assert_eq!(
            provenance.attribution.as_deref(),
            Some("CAS Common Chemistry (CC BY-NC 4.0)")
        );
        assert_eq!(substance.structure, None, "the wire carries no structure");
    }

    #[test]
    fn round_trips_through_the_wire_form() {
        let substance = Substance::from_model_json(wire_json()).expect("ingests");
        let json = substance.to_model_json().expect("emits");
        let restored = Substance::from_model_json(&json).expect("re-ingests");
        assert_eq!(restored, substance);
        // Field names match the v1 schema the other runtimes emit.
        assert!(json.contains("\"type\":\"substance-record\""));
        assert!(json.contains("\"preferredName\":\"Aspirin\""));
        assert!(json.contains("\"retrievedAt\":\"2026-09-09T00:00:00Z\""));
    }

    #[test]
    fn note_text_carries_source_and_attribution() {
        let provenance = Provenance {
            source: "common_chemistry".to_string(),
            retrieved_at: None,
            source_version: None,
            attribution: Some("CAS Common Chemistry (CC BY-NC 4.0)".to_string()),
        };
        assert_eq!(
            provenance.note_text().as_deref(),
            Some("Retrieved from common_chemistry; CAS Common Chemistry (CC BY-NC 4.0)")
        );
        let bare = Provenance {
            source: "pubchem".to_string(),
            retrieved_at: None,
            source_version: None,
            attribution: None,
        };
        assert_eq!(bare.note_text(), None);
    }

    #[test]
    fn malformed_wire_is_reported_not_guessed() {
        let error = Substance::from_model_json("{not json").expect_err("declined");
        assert!(matches!(error, WireError::Malformed(_)));
    }
}
