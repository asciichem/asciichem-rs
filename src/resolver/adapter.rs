//! The source adapter seam — the Rust mirror of the Ruby
//! `AsciiChem::Resolver::Adapter`. A source declares its name, its
//! supported conventions and its attribution, and answers
//! [`Adapter::fetch_record`]; everything around it (the TTL cache,
//! provenance stamping) is the trait's default
//! [`Adapter::resolve`].
//!
//! The transport ([`Fetch`]) is injectable so specs run offline on
//! recorded fixtures — the same seam as the Ruby specs' injected
//! fetch Structs. [`NetFetch`] is the live default: minimal stdlib
//! HTTP over ureq with rustls (the maintainer-vetted dependency),
//! 5s connect / 10s read timeouts, at most three redirects — the
//! Ruby `NetFetch` posture.

use crate::resolver::cache::Cache;
use crate::resolver::substance::{Provenance, Substance};
use crate::resolver::{timestamp, ResolverError};
use std::time::Duration;

/// The transport seam: `Ok(Some(body))` on success, `Ok(None)` for a
/// not-found answer (HTTP 404) — sources legitimately answer "unknown
/// substance" — and [`ResolverError::Source`] for every other outcome.
pub trait Fetch {
    fn get(&self, url: &str, headers: &[(&str, &str)]) -> Result<Option<String>, ResolverError>;
}

/// Minimal HTTP transport over ureq with rustls. No retries: PubChem
/// asks for politeness, and every adapter's failure posture already
/// isolates one dead endpoint.
pub struct NetFetch {
    agent: ureq::Agent,
}

impl Default for NetFetch {
    fn default() -> Self {
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(Duration::from_secs(5))
            .timeout_read(Duration::from_secs(10))
            .redirects(3)
            .build();
        NetFetch { agent }
    }
}

impl Fetch for NetFetch {
    fn get(&self, url: &str, headers: &[(&str, &str)]) -> Result<Option<String>, ResolverError> {
        let mut request = self.agent.get(url).set("User-Agent", "asciichem-resolver");
        for (name, value) in headers {
            request = request.set(name, value);
        }
        match request.call() {
            Ok(response) => response
                .into_string()
                .map(Some)
                .map_err(|error| ResolverError::Source(error.to_string())),
            Err(ureq::Error::Status(404, _)) => Ok(None),
            Err(ureq::Error::Status(code, response)) => Err(ResolverError::Source(format!(
                "HTTP {code}: {}",
                response.status_text()
            ))),
            Err(ureq::Error::Transport(transport)) => {
                let message = if transport.kind() == ureq::ErrorKind::TooManyRedirects {
                    "too many redirects".to_string()
                } else {
                    transport.to_string()
                };
                Err(ResolverError::Source(message))
            }
        }
    }
}

/// Percent-encoding for URL path/query values, byte-for-byte the Ruby
/// `CGI.escape` table: unreserved alphanumerics plus `_`, `.`, `-`
/// pass through, space becomes `+`, everything else (byte-wise, so
/// multi-byte UTF-8 splits) becomes `%XX` uppercase.
pub(crate) fn cgi_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for byte in value.bytes() {
        match byte {
            b'a'..=b'z' | b'A'..=b'Z' | b'0'..=b'9' | b'_' | b'.' | b'-' => {
                escaped.push(byte as char)
            }
            b' ' => escaped.push('+'),
            _ => escaped.push_str(&format!("%{byte:02X}")),
        }
    }
    escaped
}

/// A source adapter. The shipped sources register themselves into
/// [`crate::resolver::Registry::shipped`]; extension sources register
/// the same shape through [`crate::resolver::Registry::register`].
pub trait Adapter {
    /// The source name — the registry key and the provenance `source`.
    fn source_name(&self) -> &'static str;

    /// The identifier conventions this source answers.
    fn supported_conventions(&self) -> &'static [&'static str];

    /// The source's required attribution (licence terms — CAS Common
    /// Chemistry's CC BY-NC 4.0 line must survive into the rendered
    /// document). None when the source name itself is the whole story.
    fn attribution(&self) -> Option<String> {
        None
    }

    /// Whether this source answers `convention`.
    fn supports(&self, convention: &str) -> bool {
        self.supported_conventions().contains(&convention)
    }

    /// The record this source answers for (value, convention), or
    /// `Ok(None)` when it does not know the substance — not-found is
    /// nil, not an error, because sources legitimately answer unknown
    /// during multi-source resolution.
    fn fetch_record(
        &self,
        value: &str,
        convention: &str,
        fetch: &dyn Fetch,
    ) -> Result<Option<Substance>, ResolverError>;

    /// The provenance stamped on identifiers that arrive without one
    /// (the extension seam's backstop; the shipped sources stamp at
    /// construction already).
    fn default_provenance(&self) -> Provenance {
        Provenance {
            source: self.source_name().to_string(),
            retrieved_at: Some(timestamp::now_utc()),
            source_version: None,
            attribution: self.attribution(),
        }
    }

    /// Cache-aware resolution around [`Adapter::fetch_record`]:
    /// serve from the TTL cache unless `refresh`, cache the answer
    /// when there is one, and stamp identifiers that arrived without
    /// provenance.
    fn resolve(
        &self,
        value: &str,
        convention: &str,
        fetch: &dyn Fetch,
        cache: &Cache,
        refresh: bool,
    ) -> Result<Option<Substance>, ResolverError> {
        let key = cache.key_for(self.source_name(), convention, value);
        if !refresh {
            if let Some(cached) = cache.read(&key) {
                return Ok(Some(cached));
            }
        }
        let mut substance = self.fetch_record(value, convention, fetch)?;
        if let Some(substance) = substance.as_mut() {
            let default_provenance = self.default_provenance();
            for identifier in &mut substance.identifiers {
                if identifier.provenance.is_none() {
                    identifier.provenance = Some(default_provenance.clone());
                }
            }
            cache.write(&key, substance);
        }
        Ok(substance)
    }
}

/// Best-effort structure from a source SMILES: parsed when the SMILES
/// sits within the supported subset ([`crate::smiles`]), `None`
/// otherwise — the SMILES stays available as an identifier, nothing is
/// dropped. Shared by all three shipped adapters.
pub(crate) fn parse_structure(smiles: Option<&str>) -> Option<crate::model::Molecule> {
    let smiles = smiles?;
    let nodes = crate::smiles::parse_smiles(smiles).ok()?;
    nodes.into_iter().find_map(|node| match node {
        crate::model::Node::Molecule(molecule) => Some(molecule),
        _ => None,
    })
}

#[cfg(test)]
pub(crate) mod testing {
    //! Recorded-response transports for the offline specs (the Rust
    //! mirror of the Ruby specs' fetch Structs) — value fakes carrying
    //! data, never behaviour mocks.

    use super::{Fetch, ResolverError};
    use std::cell::RefCell;

    /// Answers every request with the recorded body (`None` = 404).
    pub(crate) struct RecordedFetch {
        pub(crate) body: Option<String>,
    }

    impl Fetch for RecordedFetch {
        fn get(
            &self,
            _url: &str,
            _headers: &[(&str, &str)],
        ) -> Result<Option<String>, ResolverError> {
            Ok(self.body.clone())
        }
    }

    /// One captured request: the URL and the exact header pairs.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub(crate) struct CapturedCall {
        pub(crate) url: String,
        pub(crate) headers: Vec<(String, String)>,
    }

    /// Captures requested URLs and headers; answers not-found.
    pub(crate) struct CaptureFetch {
        pub(crate) captured: RefCell<Vec<CapturedCall>>,
    }

    impl CaptureFetch {
        pub(crate) fn new() -> Self {
            CaptureFetch {
                captured: RefCell::new(Vec::new()),
            }
        }

        pub(crate) fn calls(&self) -> Vec<CapturedCall> {
            self.captured.borrow().clone()
        }
    }

    impl Fetch for CaptureFetch {
        fn get(
            &self,
            url: &str,
            headers: &[(&str, &str)],
        ) -> Result<Option<String>, ResolverError> {
            self.captured.borrow_mut().push(CapturedCall {
                url: url.to_string(),
                headers: headers
                    .iter()
                    .map(|(name, value)| (name.to_string(), value.to_string()))
                    .collect(),
            });
            Ok(None)
        }
    }

    /// Reads a recorded fixture from the shared test fixtures.
    pub(crate) fn fixture(relpath: &str) -> String {
        std::fs::read_to_string(format!(
            "{}/tests/fixtures/resolver/{relpath}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .expect("recorded fixture present")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cgi_escape_matches_the_ruby_table() {
        assert_eq!(cgi_escape("aspirin"), "aspirin");
        assert_eq!(cgi_escape("50-78-2"), "50-78-2");
        assert_eq!(cgi_escape("2,3-dimethylbutane"), "2%2C3-dimethylbutane");
        assert_eq!(cgi_escape("acetyl salicylic acid"), "acetyl+salicylic+acid");
        // SMILES punctuation, byte-wise.
        assert_eq!(cgi_escape("CC(=O)O"), "CC%28%3DO%29O");
        assert_eq!(cgi_escape("C[C@H](N)"), "C%5BC%40H%5D%28N%29");
        // Non-ASCII splits byte-wise, like CGI.escape.
        assert_eq!(cgi_escape("é"), "%C3%A9");
    }
}
