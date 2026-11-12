//! OPSIN source adapter (TODO.v2/13 item 3) — the Rust mirror of the
//! Ruby `AsciiChem::Resolver::Opsin`: IUPAC name → structure,
//! name-to-structure perception for systematic names, the piece
//! PubChem's name lookup cannot do (PubChem indexes a curated name
//! list; OPSIN *derives* the structure from any systematic IUPAC
//! name). Uses the public OPSIN web service
//! (<https://opsin.ch.cam.ac.uk>), which answers with the derived
//! SMILES as the whole body.
//!
//! Citation profile: publisher University of Cambridge Centre for
//! Molecular Informatics; per-access retrievedAt provenance.

use crate::resolver::adapter::{cgi_escape, parse_structure, Adapter, Fetch};
use crate::resolver::substance::{Identifier, Substance};
use crate::resolver::ResolverError;

const BASE: &str = "https://opsin.ch.cam.ac.uk/opsin";

/// The OPSIN source.
#[derive(Debug, Clone, Copy, Default)]
pub struct Opsin;

impl Adapter for Opsin {
    fn source_name(&self) -> &'static str {
        "opsin"
    }

    /// Derivation only: `name` is the curated-lookup convention
    /// (PubChem/Common Chemistry). OPSIN DERIVES structures from
    /// systematic names — reachable via `iupac-name` so it can never
    /// preempt a lookup source.
    fn supported_conventions(&self) -> &'static [&'static str] {
        &["iupac-name"]
    }

    fn attribution(&self) -> Option<String> {
        Some("OPSIN, University of Cambridge Centre for Molecular Informatics".to_string())
    }

    fn fetch_record(
        &self,
        value: &str,
        _convention: &str,
        fetch: &dyn Fetch,
    ) -> Result<Option<Substance>, ResolverError> {
        let url = format!("{BASE}/{}.smiles", cgi_escape(value));
        let body = fetch
            .get(&url, &[])
            .map_err(|error| ResolverError::Source(format!("opsin request failed: {error}")))?;
        let Some(body) = body else {
            return Ok(None);
        };
        Ok(substance_from(self, body.trim()))
    }
}

fn substance_from(adapter: &Opsin, smiles: &str) -> Option<Substance> {
    if smiles.is_empty() {
        return None;
    }
    let provenance = adapter.default_provenance();
    Some(Substance {
        preferred_name: None,
        synonyms: Vec::new(),
        formula: None,
        molecular_weight: None,
        identifiers: vec![Identifier {
            value: smiles.to_string(),
            convention: "canonical-smiles".to_string(),
            provenance: Some(provenance),
        }],
        properties: Vec::new(),
        // OPSIN answers SMILES that may use features outside the
        // supported subset (wildcards, exotic stereo); then the SMILES
        // stays as an identifier and structure is None — no silent loss.
        structure: parse_structure(Some(smiles)),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::adapter::testing::{fixture, RecordedFetch};
    use crate::resolver::cache::temp_cache;
    use crate::resolver::Registry;

    fn resolve_with(body: &str) -> Option<Substance> {
        Opsin
            .resolve(
                "2,3-dimethylbutane",
                "iupac-name",
                &RecordedFetch {
                    body: Some(body.to_string()),
                },
                &temp_cache("opsin"),
                false,
            )
            .expect("resolves")
    }

    #[test]
    fn registers_itself_as_the_opsin_source() {
        let registry = Registry::shipped();
        assert!(registry.source_names().contains(&"opsin"));
        let adapter = registry.adapter("opsin").expect("registered");
        assert_eq!(adapter.source_name(), "opsin");
    }

    #[test]
    fn supports_the_name_conventions_it_derives_from() {
        assert!(Opsin.supports("iupac-name"));
        assert!(!Opsin.supports("name"));
        assert!(!Opsin.supports("cas"));
    }

    #[test]
    fn derives_the_structure_from_a_systematic_iupac_name() {
        // Recorded OPSIN response body for "2,3-dimethylbutane" — the
        // service answers with the derived SMILES as the whole body.
        let substance = resolve_with(&fixture("opsin/2-3-dimethylbutane.txt")).expect("answers");
        assert!(substance.structure.is_some());
        assert_eq!(
            substance.identifier_value("canonical-smiles"),
            Some("CC(C)C(C)")
        );
    }

    #[test]
    fn attributes_the_cambridge_centre_for_molecular_informatics() {
        let substance = resolve_with("CC(C)C(C)\n").expect("answers");
        let provenance = substance.provenance().expect("carries provenance");
        assert_eq!(provenance.source, "opsin");
        assert!(provenance
            .attribution
            .as_deref()
            .expect("attributed")
            .contains("University of Cambridge"));
    }

    #[test]
    fn keeps_the_smiles_as_an_identifier_when_the_derived_structure_is_outside_the_subset() {
        let substance = resolve_with("C[*]C\n").expect("answers");
        assert_eq!(substance.structure, None);
        assert!(substance
            .identifier_value("canonical-smiles")
            .expect("kept")
            .contains("[*]"));
    }

    #[test]
    fn an_empty_derivation_answers_unknown() {
        assert_eq!(resolve_with("\n"), None);
    }
}
