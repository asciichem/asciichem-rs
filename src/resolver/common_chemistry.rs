//! CAS Common Chemistry source adapter (TODO.impl 40; TODO.v2 07) —
//! the Rust mirror of the Ruby `AsciiChem::Resolver::CommonChemistry`.
//!
//! License posture (maintainer decision 2026-11-11): Common Chemistry
//! data is CC BY-NC 4.0 — acceptable to enable by default; every
//! result carries the required attribution, and the citation track
//! carries it into rendered documents. Documents that must stay
//! NC-free limit sources with the `asciichem-sources` attribute
//! (metanorma-plugin-asciichem) or `Registry::resolve_all(sources)`.
//!
//! API access: CAS key-gated the API (401 without a key; request one
//! at cas.org/services/commonchemistry-api). The access email defines
//! the contract: GET /direct-api/detail with the token as the
//! X-API-KEY header plus x-origin and accept headers. The key comes
//! per instance ([`CommonChemistry::with_key`]), from
//! `COMMON_CHEMISTRY_API_KEY`, or from
//! `~/.config/cas-common-chemistry-api.key` (the token's stored
//! location) — the chain [`CommonChemistry::new`] runs, env first,
//! file second. Tokens are personal and expire after 12 months (CAS
//! emails a renewal notice). Without a key, live lookups answer 401
//! and the resolver skips this source; seeded caches keep working.

use crate::resolver::adapter::{cgi_escape, parse_structure, Adapter, Fetch};
use crate::resolver::substance::{Identifier, Property, Substance};
use crate::resolver::ResolverError;
use std::path::PathBuf;

const BASE: &str = "https://commonchemistry.cas.org/direct-api/detail";
const ORIGIN: &str = "https://commonchemistry.cas.org/api-overview";

/// The CAS Common Chemistry source. The API key rides every request
/// as the X-API-KEY header when present; None sends no key (explicit
/// nil = no key, the Ruby semantics).
#[derive(Debug, Clone)]
pub struct CommonChemistry {
    api_key: Option<String>,
}

impl CommonChemistry {
    /// The default key chain: `COMMON_CHEMISTRY_API_KEY`, else the
    /// stored token file (when readable, trimmed). The lookup runs
    /// now, at construction — exactly when the Ruby
    /// `CommonChemistry.new` default argument evaluates.
    pub fn new() -> Self {
        Self::from_lookup(
            std::env::var("COMMON_CHEMISTRY_API_KEY").ok(),
            Self::default_key_file(),
        )
    }

    /// An explicit key: rides every request as X-API-KEY.
    pub fn with_key(api_key: impl Into<String>) -> Self {
        Self::from_lookup(Some(api_key.into()), Self::default_key_file())
    }

    /// Explicit nil: no key, and the env/file lookup never runs.
    pub fn keyless() -> Self {
        CommonChemistry { api_key: None }
    }

    /// The key chain itself: an explicit env answer first, the stored
    /// token file second (when readable, stripped of surrounding
    /// whitespace). Shared by `new` and the specs — the specs' Rust
    /// shape of the Ruby KEY_FILE const swap.
    pub(crate) fn from_lookup(env_key: Option<String>, key_file: PathBuf) -> Self {
        let api_key = match env_key {
            Some(key) => Some(key),
            None => std::fs::read_to_string(&key_file)
                .ok()
                .map(|text| text.trim().to_string())
                .filter(|key| !key.is_empty()),
        };
        CommonChemistry { api_key }
    }

    /// The token's stored location: `~/.config/cas-common-chemistry-api.key`.
    pub fn default_key_file() -> PathBuf {
        let home = std::env::var("HOME")
            .or_else(|_| std::env::var("USERPROFILE"))
            .unwrap_or_default();
        PathBuf::from(home)
            .join(".config")
            .join("cas-common-chemistry-api.key")
    }

    /// The resolved key — the specs' assertion surface.
    pub fn api_key(&self) -> Option<&str> {
        self.api_key.as_deref()
    }
}

impl Default for CommonChemistry {
    fn default() -> Self {
        Self::new()
    }
}

impl Adapter for CommonChemistry {
    fn source_name(&self) -> &'static str {
        "common_chemistry"
    }

    fn supported_conventions(&self) -> &'static [&'static str] {
        &["cas"]
    }

    /// CC BY-NC 4.0 — REQUIRED on every record this source answers.
    fn attribution(&self) -> Option<String> {
        Some("CAS Common Chemistry (CC BY-NC 4.0)".to_string())
    }

    fn fetch_record(
        &self,
        value: &str,
        _convention: &str,
        fetch: &dyn Fetch,
        now: Option<&str>,
    ) -> Result<Option<Substance>, ResolverError> {
        let mut headers: Vec<(&str, &str)> =
            vec![("Accept", "application/json"), ("x-origin", ORIGIN)];
        if let Some(api_key) = &self.api_key {
            headers.push(("X-API-KEY", api_key));
        }
        let url = format!("{BASE}?cas_rn={}", cgi_escape(value));
        let Some(body) = fetch.get(&url, &headers)? else {
            return Ok(None);
        };
        let data: serde_json::Value = serde_json::from_str(&body).map_err(|error| {
            ResolverError::Source(format!("common_chemistry returned malformed JSON: {error}"))
        })?;
        if !data.is_object() || data.as_object().is_none_or(|object| object.is_empty()) {
            return Ok(None);
        }
        Ok(substance_from(self, &data, now))
    }
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

fn substance_from(
    adapter: &CommonChemistry,
    data: &serde_json::Value,
    now: Option<&str>,
) -> Option<Substance> {
    let provenance = adapter.default_provenance(now);
    let mut identifiers = vec![Identifier {
        value: json_scalar(data.get("rn"))?,
        convention: "cas".to_string(),
        provenance: Some(provenance.clone()),
    }];
    for (key, convention) in [
        ("inchikey", "inchikey"),
        ("smiles", "canonical-smiles"),
        ("inchi", "inchi"),
        ("iupacName", "iupac-name"),
    ] {
        if let Some(value) = json_scalar(data.get(key)) {
            identifiers.push(Identifier {
                value,
                convention: convention.to_string(),
                provenance: Some(provenance.clone()),
            });
        }
    }

    let formula = json_scalar(data.get("formula"));
    let mut properties = Vec::new();
    if let Some(value) = &formula {
        properties.push(Property {
            name: "molecular-formula".to_string(),
            value: value.clone(),
            units: None,
            provenance: Some(provenance.clone()),
        });
    }
    if let Some(value) = json_scalar(data.get("molecularMass")) {
        properties.push(Property {
            name: "molecular-weight".to_string(),
            value,
            units: Some("g/mol".to_string()),
            provenance: Some(provenance.clone()),
        });
    }

    let smiles = json_scalar(data.get("smiles"));
    Some(Substance {
        preferred_name: json_scalar(data.get("name")),
        synonyms: Vec::new(),
        formula,
        molecular_weight: None,
        identifiers,
        properties,
        // Best effort within the supported SMILES subset; the SMILES
        // stays as an identifier when stereo keeps it out.
        structure: parse_structure(smiles.as_deref()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::adapter::testing::{fixture, CaptureFetch, RecordedFetch};
    use crate::resolver::cache::temp_cache;
    use crate::resolver::Registry;
    use std::sync::{Mutex, OnceLock};

    fn cc_fixture() -> String {
        fixture("common_chemistry/aspirin.json")
    }

    fn key_file_fixture() -> PathBuf {
        PathBuf::from(format!(
            "{}/tests/fixtures/resolver/common_chemistry/api.key",
            env!("CARGO_MANIFEST_DIR")
        ))
    }

    /// The one test that mutates the process env holds this lock;
    /// the rest of the chain is exercised hermetically through
    /// `from_lookup`.
    fn env_lock() -> &'static Mutex<()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    #[test]
    fn self_registers_default_on_per_the_maintainer_decision_of_2026_11_11() {
        let registry = Registry::shipped();
        assert!(registry.source_names().contains(&"common_chemistry"));
    }

    #[test]
    fn resolves_by_cas_rn() {
        let registry = Registry::shipped();
        let substance = registry
            .resolve(
                "50-78-2",
                "cas",
                Some(&["common_chemistry"]),
                &RecordedFetch {
                    body: Some(cc_fixture()),
                },
                &temp_cache("cc"),
                false,
                None,
            )
            .expect("resolves")
            .expect("answers");
        assert_eq!(substance.preferred_name.as_deref(), Some("Aspirin"));
        assert_eq!(substance.identifier_value("cas"), Some("50-78-2"));
        assert_eq!(
            substance.identifier_value("inchikey"),
            Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        );
    }

    #[test]
    fn carries_the_cc_by_nc_attribution_on_every_identifier() {
        let substance = CommonChemistry::keyless()
            .resolve(
                "50-78-2",
                "cas",
                &RecordedFetch {
                    body: Some(cc_fixture()),
                },
                &temp_cache("cc-attrib"),
                false,
                None,
            )
            .expect("resolves")
            .expect("answers");
        assert!(!substance.identifiers.is_empty());
        for identifier in &substance.identifiers {
            assert_eq!(
                identifier
                    .provenance
                    .as_ref()
                    .expect("stamped")
                    .attribution
                    .as_deref(),
                Some("CAS Common Chemistry (CC BY-NC 4.0)"),
                "{}",
                identifier.convention
            );
        }
    }

    #[test]
    fn parses_the_structure_from_the_source_smiles() {
        let substance = CommonChemistry::keyless()
            .resolve(
                "50-78-2",
                "cas",
                &RecordedFetch {
                    body: Some(cc_fixture()),
                },
                &temp_cache("cc-struct"),
                false,
                None,
            )
            .expect("resolves")
            .expect("answers");
        let structure = substance.structure.as_ref().expect("parses");
        let smiles =
            crate::smiles::write_smiles(&[crate::model::Node::Molecule(structure.clone())])
                .expect("writes");
        assert_eq!(smiles, "CC(=O)Oc1ccccc1C(=O)O");
    }

    #[test]
    fn queries_the_direct_api_detail_endpoint_by_cas_rn() {
        let capture = CaptureFetch::new();
        CommonChemistry::keyless()
            .resolve(
                "50-78-2",
                "cas",
                &capture,
                &temp_cache("cc-url"),
                false,
                None,
            )
            .expect("not-found is not an error");
        let url = capture.calls().pop().expect("requested").url;
        assert!(url.contains("/direct-api/detail?cas_rn=50-78-2"), "{url}");
        assert!(!url.contains("apikey"), "{url}");
    }

    #[test]
    fn sends_the_x_origin_and_accept_headers_on_every_request() {
        let capture = CaptureFetch::new();
        CommonChemistry::keyless()
            .resolve(
                "50-78-2",
                "cas",
                &capture,
                &temp_cache("cc-headers"),
                false,
                None,
            )
            .expect("not-found is not an error");
        let headers = capture.calls().pop().expect("requested").headers;
        let header = |name: &str| {
            headers
                .iter()
                .find(|(key, _)| key == name)
                .map(|(_, value)| value.as_str())
        };
        assert_eq!(header("x-origin"), Some(ORIGIN));
        assert_eq!(header("Accept"), Some("application/json"));
    }

    #[test]
    fn sends_a_user_provided_key_as_the_x_api_key_header() {
        let capture = CaptureFetch::new();
        CommonChemistry::with_key("cas-key-123")
            .resolve(
                "50-78-2",
                "cas",
                &capture,
                &temp_cache("cc-key"),
                false,
                None,
            )
            .expect("not-found is not an error");
        let headers = capture.calls().pop().expect("requested").headers;
        assert!(
            headers.contains(&("X-API-KEY".to_string(), "cas-key-123".to_string())),
            "{headers:?}"
        );
    }

    #[test]
    fn sends_no_x_api_key_without_a_key() {
        let capture = CaptureFetch::new();
        CommonChemistry::keyless()
            .resolve(
                "50-78-2",
                "cas",
                &capture,
                &temp_cache("cc-nokey"),
                false,
                None,
            )
            .expect("not-found is not an error");
        let headers = capture.calls().pop().expect("requested").headers;
        assert!(
            headers.iter().all(|(name, _)| name != "X-API-KEY"),
            "{headers:?}"
        );
    }

    #[test]
    fn picks_the_key_up_from_common_chemistry_api_key() {
        // env answers before the file is even consulted.
        let adapter =
            CommonChemistry::from_lookup(Some("env-key-456".to_string()), key_file_fixture());
        assert_eq!(adapter.api_key(), Some("env-key-456"));
    }

    #[test]
    fn falls_back_to_the_stored_key_file() {
        let adapter = CommonChemistry::from_lookup(None, key_file_fixture());
        assert_eq!(adapter.api_key(), Some("file-stored-token"));
    }

    #[test]
    fn no_key_anywhere_means_no_key() {
        let adapter = CommonChemistry::from_lookup(None, "/nonexistent/cas-api.key".into());
        assert_eq!(adapter.api_key(), None);
        assert_eq!(CommonChemistry::keyless().api_key(), None);
    }

    #[test]
    fn the_default_chain_reads_the_process_env() {
        let _guard = env_lock().lock().expect("env lock");
        std::env::set_var("COMMON_CHEMISTRY_API_KEY", "live-env-key");
        assert_eq!(CommonChemistry::new().api_key(), Some("live-env-key"));
        std::env::remove_var("COMMON_CHEMISTRY_API_KEY");
    }
}
