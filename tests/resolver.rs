//! The resolver contract (TODO.v2/17 phase 1, slice 5; asciichem-rs#13):
//! the registry's fan-out semantics — first answering source, one
//! record per source, per-source failure isolation, cross-check
//! conflict, the TTL cache serving offline, and the canonical wire
//! substance-record form — mirrored from the Ruby `resolver_spec`,
//! gated offline on the same recorded fixtures (no network access;
//! every fetch rides the injected `Fetch` seam).

use asciichem::resolver::adapter::{Adapter, Fetch};
use asciichem::resolver::cache::Cache;
use asciichem::resolver::substance::{Identifier, Provenance, Substance};
use asciichem::resolver::{Registry, ResolverError};
use std::cell::RefCell;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn fixture(relpath: &str) -> String {
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "fixtures", "resolver"]
        .iter()
        .collect();
    std::fs::read_to_string(path.join(relpath)).expect("recorded fixture present")
}

fn temp_cache(tag: &str) -> Cache {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    Cache::at_dir(std::env::temp_dir().join(format!("asciichem-resolver-it-{tag}-{unique}")))
}

/// Answers every request with the recorded body (the Ruby specs'
/// fetch Struct — a value fake carrying data).
struct RecordedFetch {
    body: Option<String>,
}

impl Fetch for RecordedFetch {
    fn get(&self, _url: &str, _headers: &[(&str, &str)]) -> Result<Option<String>, ResolverError> {
        Ok(self.body.clone())
    }
}

/// Network access from an offline lookup would be a bug.
struct OfflineFetch;

impl Fetch for OfflineFetch {
    fn get(&self, url: &str, _headers: &[(&str, &str)]) -> Result<Option<String>, ResolverError> {
        Err(ResolverError::Source(format!(
            "network must not be touched: {url}"
        )))
    }
}

/// Builds a one-off fixture-backed source registered under `name`
/// (the Ruby spec's `fixture_source`).
#[derive(Clone)]
struct FixtureSource {
    name: &'static str,
    body: String,
}

impl Adapter for FixtureSource {
    fn source_name(&self) -> &'static str {
        self.name
    }

    fn supported_conventions(&self) -> &'static [&'static str] {
        &["name"]
    }

    fn fetch_record(
        &self,
        _value: &str,
        _convention: &str,
        _fetch: &dyn Fetch,
    ) -> Result<Option<Substance>, ResolverError> {
        let payload: serde_json::Value = serde_json::from_str(&self.body)
            .map_err(|error| ResolverError::Source(error.to_string()))?;
        let properties = payload
            .get("PropertyTable")
            .and_then(|table| table.get("Properties"))
            .and_then(|list| list.get(0));
        let Some(properties) = properties else {
            return Ok(None);
        };
        let inchikey = properties
            .get("InChIKey")
            .and_then(|value| value.as_str())
            .map(str::to_string);
        let preferred_name = properties
            .get("IUPACName")
            .and_then(|value| value.as_str())
            .map(str::to_string);
        let Some(inchikey) = inchikey else {
            return Ok(None);
        };
        Ok(Some(Substance {
            preferred_name,
            synonyms: Vec::new(),
            formula: None,
            molecular_weight: None,
            identifiers: vec![Identifier {
                value: inchikey,
                convention: "inchikey".to_string(),
                provenance: Some(Provenance {
                    source: self.name.to_string(),
                    retrieved_at: None,
                    source_version: None,
                    attribution: None,
                }),
            }],
            properties: Vec::new(),
            structure: None,
        }))
    }
}

/// A source-level failure standing in for a key-gated or down
/// endpoint (HTTP 401/500).
struct FailingSource {
    name: &'static str,
}

impl Adapter for FailingSource {
    fn source_name(&self) -> &'static str {
        self.name
    }

    fn supported_conventions(&self) -> &'static [&'static str] {
        &["name"]
    }

    fn fetch_record(
        &self,
        _value: &str,
        _convention: &str,
        _fetch: &dyn Fetch,
    ) -> Result<Option<Substance>, ResolverError> {
        Err(ResolverError::Source("HTTP 401: Unauthorized".to_string()))
    }
}

/// One captured request: the URL and the exact header pairs.
type CapturedCall = (String, Vec<(String, String)>);

/// Captures requested URLs and headers; answers not-found.
struct CaptureFetch {
    captured: RefCell<Vec<CapturedCall>>,
}

impl Fetch for CaptureFetch {
    fn get(&self, url: &str, headers: &[(&str, &str)]) -> Result<Option<String>, ResolverError> {
        self.captured.borrow_mut().push((
            url.to_string(),
            headers
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
        ));
        Ok(None)
    }
}

fn registry_with(extra: impl FnOnce(&mut Registry)) -> Registry {
    let mut registry = Registry::shipped();
    extra(&mut registry);
    registry
}

#[test]
fn resolves_via_the_first_answering_source_by_default() {
    let registry = registry_with(|_| {});
    let substance = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("first"),
            false,
        )
        .expect("resolves")
        .expect("answers");
    assert_eq!(
        substance.provenance().expect("provenance").source,
        "pubchem"
    );
}

#[test]
fn returns_one_substance_per_answering_source() {
    let registry = registry_with(|registry| {
        registry.register("spec_mirror", || {
            Box::new(FixtureSource {
                name: "spec_mirror",
                body: fixture("pubchem/aspirin.json"),
            })
        });
    });
    let results = registry
        .resolve_all(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("fanout"),
            false,
        )
        .expect("resolves");
    // Membership, not a closed world: every registered adapter
    // answering "name" joins the fan-out (opsin did exactly that).
    let mut sources: Vec<&str> = results
        .iter()
        .map(|substance| substance.provenance().expect("provenance").source.as_str())
        .collect();
    assert!(sources.contains(&"pubchem"), "{sources:?}");
    assert!(sources.contains(&"spec_mirror"), "{sources:?}");
    sources.sort_unstable();
    let mut unique = sources.clone();
    unique.dedup();
    assert_eq!(
        unique.len(),
        sources.len(),
        "no duplicate sources: {sources:?}"
    );
}

#[test]
fn skips_a_source_level_failure_instead_of_poisoning_the_fan_out() {
    let registry = registry_with(|registry| {
        registry.register("spec_failing", || {
            Box::new(FailingSource {
                name: "spec_failing",
            })
        });
    });
    let substance = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("isolate"),
            false,
        )
        .expect("the failure isolates")
        .expect("answers");
    assert_eq!(
        substance.provenance().expect("provenance").source,
        "pubchem"
    );
}

#[test]
fn raises_conflict_when_sources_disagree_on_the_inchikey() {
    let disagreement = fixture("pubchem/aspirin.json")
        .replace("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "WRONGKEY-WRONGWRONG-A");
    let registry = registry_with(|registry| {
        registry.register("spec_mirror", move || {
            Box::new(FixtureSource {
                name: "spec_mirror",
                body: disagreement.clone(),
            })
        });
    });
    let error = registry
        .resolve_checked(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("conflict"),
            false,
        )
        .expect_err("the disagreement raises");
    match error {
        ResolverError::Conflict(message) => {
            assert!(message.contains("disagree"), "{message}");
            assert!(message.contains("BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), "{message}");
            assert!(message.contains("WRONGKEY-WRONGWRONG-A"), "{message}");
        }
        other => panic!("expected a conflict: {other:?}"),
    }
}

#[test]
fn agrees_when_sources_cross_check_clean() {
    let registry = registry_with(|registry| {
        registry.register("spec_mirror", || {
            Box::new(FixtureSource {
                name: "spec_mirror",
                body: fixture("pubchem/aspirin.json"),
            })
        });
    });
    let substance = registry
        .resolve_checked(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("crosscheck"),
            false,
        )
        .expect("agreement passes")
        .expect("answers");
    assert_eq!(
        substance.identifier_value("inchikey"),
        Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    );
}

#[test]
fn caches_resolution_results_in_the_user_cache_dir_and_serves_them_offline() {
    let cache = temp_cache("cache");
    let registry = registry_with(|_| {});
    let substance = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &cache,
            false,
        )
        .expect("resolves")
        .expect("answers");
    assert_eq!(
        substance.identifier_value("inchikey"),
        Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    );
    assert!(cache.dir().exists(), "the cache dir came up");

    let cached = registry
        .resolve("aspirin", "name", None, &OfflineFetch, &cache, false)
        .expect("served offline")
        .expect("cache answers");
    assert_eq!(
        cached.identifier_value("inchikey"),
        Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    );

    let fresh = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &cache,
            true,
        )
        .expect("refreshes")
        .expect("answers");
    assert_eq!(fresh.identifier_value("pubchem-cid"), Some("2244"));
}

#[test]
fn serializes_substance_records_to_the_canonical_wire_form() {
    let registry = registry_with(|_| {});
    let substance = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("wire"),
            false,
        )
        .expect("resolves")
        .expect("answers");
    let wire = substance.to_model_json().expect("emits");
    let json: serde_json::Value = serde_json::from_str(&wire).expect("parses");
    assert_eq!(json["type"], "substance-record");
    assert_eq!(json["preferredName"], "2-acetyloxybenzoic acid");
    assert_eq!(
        json["identifiers"][0]["identifier"]["convention"],
        "pubchem-cid"
    );
    assert_eq!(json["identifiers"][0]["provenance"]["source"], "pubchem");
}

#[test]
fn round_trips_a_substance_through_the_wire_form_the_cache_path() {
    let registry = registry_with(|_| {});
    let substance = registry
        .resolve(
            "aspirin",
            "name",
            None,
            &RecordedFetch {
                body: Some(fixture("pubchem/aspirin.json")),
            },
            &temp_cache("roundtrip"),
            false,
        )
        .expect("resolves")
        .expect("answers");
    let restored =
        Substance::from_model_json(&substance.to_model_json().expect("emits")).expect("re-ingests");
    assert_eq!(
        restored.identifier_value("inchikey"),
        substance.identifier_value("inchikey")
    );
    assert_eq!(restored.preferred_name, substance.preferred_name);
}

#[test]
fn raises_an_actionable_error_for_unknown_sources() {
    let registry = registry_with(|_| {});
    let error = match registry.adapter("nope") {
        Err(error) => error,
        Ok(_) => panic!("unknown names are errors"),
    };
    match error {
        ResolverError::Source(message) => {
            assert!(message.contains("unknown resolver source"), "{message}");
            assert!(message.contains("\"nope\""), "{message}");
            assert!(
                message.contains("common_chemistry, opsin, pubchem"),
                "{message}"
            );
        }
        other => panic!("expected a source error: {other:?}"),
    }
}

#[test]
fn sources_limiting_narrows_the_fan_out() {
    // The NC-free posture: limiting to a named subset must keep the
    // CC-licensed source out of the answer set entirely.
    let capture = CaptureFetch {
        captured: RefCell::new(Vec::new()),
    };
    let registry = registry_with(|_| {});
    let answered = registry
        .resolve(
            "50-78-2",
            "cas",
            Some(&["pubchem"]),
            &capture,
            &temp_cache("limiting"),
            false,
        )
        .expect("resolves");
    assert_eq!(
        answered, None,
        "pubchem does not answer the recorded capture"
    );
    for (url, _) in capture.captured.borrow().iter() {
        assert!(url.contains("pubchem.ncbi.nlm.nih.gov"), "{url}");
        assert!(!url.contains("commonchemistry.cas.org"), "{url}");
    }
}
