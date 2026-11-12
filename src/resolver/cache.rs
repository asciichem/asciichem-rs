//! Resolution cache in the USER cache dir (never the crate install
//! dir — a library never writes into its own installed location): the
//! Rust mirror of the Ruby `Resolver::Cache`. Entries carry the
//! retrieval time, the TTL, and the substance in the canonical wire
//! substance-record form; keys are sha256 over `source:convention:value`
//! (the Ruby cache's `Digest::SHA256`, through the `sha2` crate).
//! A cache failure (unreadable entry, malformed JSON) degrades to a
//! miss — it must never kill a live resolution.

use crate::resolver::substance::Substance;
use crate::resolver::timestamp;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

/// Chemistry identifiers are stable — one week, the Ruby default.
pub const DEFAULT_TTL_SECONDS: u64 = 7 * 24 * 60 * 60;

#[derive(Serialize, Deserialize)]
struct Entry {
    cached_at: String,
    ttl_seconds: Option<u64>,
    substance: crate::wire::SubstanceRecord,
}

/// The TTL cache: a directory of wire-form substance records.
#[derive(Debug, Clone)]
pub struct Cache {
    dir: PathBuf,
    ttl: u64,
}

impl Cache {
    /// The default cache: the user cache dir, the default TTL.
    pub fn new() -> Self {
        Cache {
            dir: Self::default_dir(),
            ttl: DEFAULT_TTL_SECONDS,
        }
    }

    /// A cache rooted at `dir` — the specs' isolation seam (and any
    /// caller that wants resolution state out of the user cache).
    pub fn at_dir(dir: impl Into<PathBuf>) -> Self {
        Cache {
            dir: dir.into(),
            ttl: DEFAULT_TTL_SECONDS,
        }
    }

    /// Overrides the TTL (seconds).
    pub fn with_ttl(mut self, ttl: u64) -> Self {
        self.ttl = ttl;
        self
    }

    /// The user cache dir: `XDG_CACHE_HOME`, else `HOME/.cache`, else
    /// `USERPROFILE/.cache`, else the system temp dir (the Ruby
    /// default's `Dir.home` posture, with a portable last resort).
    pub fn default_dir() -> PathBuf {
        let base = std::env::var("XDG_CACHE_HOME")
            .ok()
            .or_else(|| {
                std::env::var("HOME")
                    .ok()
                    .map(|home| format!("{home}/.cache"))
            })
            .or_else(|| {
                std::env::var("USERPROFILE")
                    .ok()
                    .map(|home| format!("{home}/.cache"))
            })
            .unwrap_or_else(|| std::env::temp_dir().to_string_lossy().into_owned());
        Path::new(&base).join("asciichem").join("resolver")
    }

    /// The cache directory this instance reads and writes.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    /// The sha256 hex of `source:convention:value` — file-name-safe
    /// and stable across runs.
    pub fn key_for(&self, source: &str, convention: &str, value: &str) -> String {
        let digest = Sha256::digest(format!("{source}:{convention}:{value}").as_bytes());
        let mut key = String::with_capacity(64);
        for byte in digest {
            use std::fmt::Write as _;
            let _ = write!(key, "{byte:02x}");
        }
        key
    }

    /// The cached substance for `key`, or None on a miss, a stale
    /// entry, or an unreadable/malformed entry (all degrade to a miss).
    pub fn read(&self, key: &str) -> Option<Substance> {
        let path = self.dir.join(format!("{key}.json"));
        let text = std::fs::read_to_string(path).ok()?;
        let entry: Entry = serde_json::from_str(&text).ok()?;
        let cached_at = timestamp::parse_utc(&entry.cached_at)?;
        let now = SystemTime::now().duration_since(UNIX_EPOCH).ok()?.as_secs() as i64;
        let ttl = entry.ttl_seconds.unwrap_or(self.ttl);
        if now - cached_at > ttl as i64 {
            return None;
        }
        Some(Substance::from_wire(&entry.substance))
    }

    /// Writes the substance under `key` as one cache entry. Best
    /// effort: an unwritable cache must not kill a live resolution
    /// (the answer is already in hand).
    pub fn write(&self, key: &str, substance: &Substance) {
        let entry = Entry {
            cached_at: timestamp::now_utc(),
            ttl_seconds: Some(self.ttl),
            substance: substance.to_wire(),
        };
        if let Ok(text) = serde_json::to_string_pretty(&entry) {
            let _ = std::fs::create_dir_all(&self.dir);
            let _ = std::fs::write(self.dir.join(format!("{key}.json")), text);
        }
    }
}

impl Default for Cache {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) fn temp_cache(tag: &str) -> Cache {
    use std::sync::atomic::{AtomicU32, Ordering};
    static SEQUENCE: AtomicU32 = AtomicU32::new(0);
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|since| since.as_nanos())
        .unwrap_or_default();
    let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    Cache::at_dir(
        std::env::temp_dir().join(format!("asciichem-resolver-spec-{tag}-{unique}-{sequence}")),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::substance::{Identifier, Provenance};

    fn substance(inchikey: &str) -> Substance {
        Substance {
            preferred_name: Some("2-acetyloxybenzoic acid".to_string()),
            synonyms: Vec::new(),
            formula: Some("C9H8O4".to_string()),
            molecular_weight: Some(180.16),
            identifiers: vec![Identifier {
                value: inchikey.to_string(),
                convention: "inchikey".to_string(),
                provenance: Some(Provenance {
                    source: "pubchem".to_string(),
                    retrieved_at: Some(timestamp::now_utc()),
                    source_version: None,
                    attribution: Some("PubChem, U.S. National Library of Medicine".to_string()),
                }),
            }],
            properties: Vec::new(),
            structure: None,
        }
    }

    #[test]
    fn writes_reads_and_serves_offline() {
        let cache = temp_cache("cache");
        let key = cache.key_for("pubchem", "name", "aspirin");
        assert_eq!(cache.read(&key), None, "cold cache misses");

        cache.write(&key, &substance("BSYNRYMUTXBXSQ-UHFFFAOYSA-N"));
        assert!(cache.dir().exists());
        let cached = cache.read(&key).expect("warm cache hits");
        assert_eq!(
            cached.identifier_value("inchikey"),
            Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        );
        assert_eq!(
            cached.preferred_name.as_deref(),
            Some("2-acetyloxybenzoic acid")
        );
        assert_eq!(cached.molecular_weight, Some(180.16));
    }

    #[test]
    fn keys_match_the_ruby_sha256() {
        // Digest::SHA256.hexdigest("pubchem:name:aspirin")
        let cache = Cache::at_dir("unused");
        assert_eq!(
            cache.key_for("pubchem", "name", "aspirin"),
            "7ccfaa13c552945e313ba36979762a3c956fb39b13a345818000cdce60c6a1b5"
        );
    }

    #[test]
    fn stale_entries_are_misses() {
        // Ruby semantics: fresh when `now - cached_at <= ttl` (a
        // same-second write+read with ttl 0 is a hit). Backdate the
        // entry instead — an expired record degrades to a miss.
        let cache = temp_cache("stale");
        let key = cache.key_for("pubchem", "name", "aspirin");
        cache.write(&key, &substance("BSYNRYMUTXBXSQ-UHFFFAOYSA-N"));
        assert!(cache.read(&key).is_some(), "fresh entry hits");

        let path = cache.dir().join(format!("{key}.json"));
        let mut entry: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&path).expect("entry file"))
                .expect("entry JSON");
        entry["cached_at"] = "2001-09-09T01:46:40Z".into();
        std::fs::write(&path, serde_json::to_string(&entry).expect("entry JSON")).expect("writes");
        assert_eq!(
            cache.read(&key),
            None,
            "the decades-old entry is stale (the week TTL)"
        );
    }
}
