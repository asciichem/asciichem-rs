//! Resolution turns identifiers into substance records (TODO.v2 07
//! Layer 2; TODO.impl 38; TODO.v2/17 phase 1 slice 5): PubChem
//! (default), OPSIN name-to-structure, CAS Common Chemistry —
//! multi-source by design. Every answering source produces its own
//! record with provenance, because the citation track (TODO.v2 08)
//! emits one bibitem per (substance, source) the user chooses to
//! cite.
//!
//! Sources are adapters ([`adapter::Adapter`]) registered in a
//! [`Registry`]. The Ruby module-level registry maps to an explicit
//! registry value: [`Registry::shipped()`] gives the three built-in
//! sources in registration order (pubchem, opsin, common_chemistry —
//! Common Chemistry is CC BY-NC 4.0 and default-on per the maintainer
//! decision of 2026-11-11; callers that need NC-free resolution limit
//! sources explicitly), [`Registry::register`] adds opt-in sources,
//! and the `sources` argument narrows a call to a named subset.
//!
//! Failure posture: a source-level failure (HTTP error, rejected
//! request, malformed payload) skips that source — one dead or
//! key-gated endpoint must not poison the others; callers see what
//! answered. [`ResolverError::Conflict`] is the exception: cross-check
//! mode ([`Registry::resolve_checked`], the Ruby `resolve!`) raises
//! when more than one answering source provides an InChIKey and they
//! disagree — the message lists the disagreement, never hidden.
//!
//! The transport is injectable through the [`adapter::Fetch`] trait —
//! the offline seam the recorded-fixture specs run through.
//! [`adapter::NetFetch`] (ureq with rustls — the maintainer-vetted
//! dependency; the engines' zero-dependency policy does not extend to
//! network I/O) is the live default: 5s connect / 10s read timeouts,
//! at most three redirects.

pub mod adapter;
pub mod cache;
pub mod common_chemistry;
pub mod opsin;
pub mod pubchem;
pub mod substance;

pub(crate) mod timestamp;

pub use adapter::{Adapter, Fetch, NetFetch};
pub use cache::Cache;
pub use substance::{Identifier, Property, Provenance, Substance};

use std::sync::Arc;

/// A resolver failure. The Ruby resolver speaks one error type
/// (`AsciiChem::Error`) plus `Resolver::Conflict`; the same two
/// shapes live here as variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResolverError {
    /// A source-level failure: transport error, rejected request,
    /// malformed payload. Isolatable — [`Registry::resolve_all`]
    /// skips the failing source instead of poisoning the fan-out.
    Source(String),
    /// Sources disagree on substance identity (cross-check mode).
    Conflict(String),
}

impl std::fmt::Display for ResolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ResolverError::Source(message) => write!(f, "{message}"),
            ResolverError::Conflict(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for ResolverError {}

/// Builds a source adapter. The factories (not instances) populate the
/// registry because construction can consult the environment — the
/// Common Chemistry key chain reads the env var and the token file at
/// build time, exactly as the Ruby registry's `adapter_class.new` does.
type Factory = Arc<dyn Fn() -> Box<dyn Adapter> + Send + Sync>;

/// The source registry: registration order, `sources:` limiting,
/// extension registration.
#[derive(Clone)]
pub struct Registry {
    adapters: Vec<(String, Factory)>,
}

impl Default for Registry {
    fn default() -> Self {
        Self::shipped()
    }
}

impl Registry {
    /// The three built-in sources in registration order: pubchem,
    /// opsin, common_chemistry. Common Chemistry is CC BY-NC 4.0 and
    /// default-on per the maintainer decision of 2026-11-11.
    pub fn shipped() -> Self {
        Registry {
            adapters: vec![
                (
                    "pubchem".to_string(),
                    Arc::new(|| -> Box<dyn Adapter> { Box::new(pubchem::PubChem) }),
                ),
                (
                    "opsin".to_string(),
                    Arc::new(|| -> Box<dyn Adapter> { Box::new(opsin::Opsin) }),
                ),
                (
                    "common_chemistry".to_string(),
                    Arc::new(|| -> Box<dyn Adapter> {
                        Box::new(common_chemistry::CommonChemistry::new())
                    }),
                ),
            ],
        }
    }

    /// Registers (or replaces, keeping the original registration
    /// position — the Ruby hash-assignment semantics) a source under
    /// `name`.
    pub fn register(
        &mut self,
        name: impl Into<String>,
        factory: impl Fn() -> Box<dyn Adapter> + Send + Sync + 'static,
    ) {
        let name = name.into();
        let factory: Factory = Arc::new(factory);
        match self
            .adapters
            .iter_mut()
            .find(|(existing, _)| *existing == name)
        {
            Some(slot) => slot.1 = factory,
            None => self.adapters.push((name, factory)),
        }
    }

    /// The registered source names, in registration order.
    pub fn source_names(&self) -> Vec<&str> {
        self.adapters
            .iter()
            .map(|(name, _)| name.as_str())
            .collect()
    }

    /// Builds the adapter registered under `name`. Unknown names are
    /// actionable errors listing what IS registered — never silent.
    pub fn adapter(&self, name: &str) -> Result<Box<dyn Adapter>, ResolverError> {
        let (_, factory) = self
            .adapters
            .iter()
            .find(|(existing, _)| existing == name)
            .ok_or_else(|| {
                let mut registered = self.source_names();
                registered.sort_unstable();
                ResolverError::Source(format!(
                    "unknown resolver source {name:?} (registered: {})",
                    registered.join(", ")
                ))
            })?;
        Ok(factory())
    }

    /// Resolves against every registered source (or the selected
    /// subset), returning one Substance per source that answers (None
    /// entries dropped). Sources answer in registration order. A
    /// source-level failure skips that source — one dead or key-gated
    /// endpoint must not poison the others; callers see what answered.
    /// An unknown name in `sources` is an error, not a skip. `now` is
    /// the injectable clock the records' provenance stamps (`None` =
    /// the system clock; replayed fixtures pin the instant).
    #[allow(clippy::too_many_arguments)]
    pub fn resolve_all(
        &self,
        value: &str,
        convention: &str,
        sources: Option<&[&str]>,
        fetch: &dyn Fetch,
        cache: &Cache,
        refresh: bool,
        now: Option<&str>,
    ) -> Result<Vec<Substance>, ResolverError> {
        let selected: Vec<Box<dyn Adapter>> = match sources {
            Some(names) => names
                .iter()
                .map(|name| self.adapter(name))
                .collect::<Result<Vec<_>, _>>()?,
            None => self.adapters.iter().map(|(_, factory)| factory()).collect(),
        };
        let mut answered = Vec::new();
        for adapter in selected {
            if !adapter.supports(convention) {
                continue;
            }
            match adapter.resolve(value, convention, fetch, cache, refresh, now) {
                Ok(Some(substance)) => answered.push(substance),
                Ok(None) => {}
                // Per-source isolation: a dead or key-gated endpoint
                // must not poison the others.
                Err(ResolverError::Source(_)) => {}
                Err(other) => return Err(other),
            }
        }
        Ok(answered)
    }

    /// First answering source wins.
    #[allow(clippy::too_many_arguments)]
    pub fn resolve(
        &self,
        value: &str,
        convention: &str,
        sources: Option<&[&str]>,
        fetch: &dyn Fetch,
        cache: &Cache,
        refresh: bool,
        now: Option<&str>,
    ) -> Result<Option<Substance>, ResolverError> {
        Ok(self
            .resolve_all(value, convention, sources, fetch, cache, refresh, now)?
            .into_iter()
            .next())
    }

    /// Cross-check mode (the Ruby `resolve!`): all sources must agree
    /// on the primary key (InChIKey) when more than one source
    /// provides one.
    #[allow(clippy::too_many_arguments)]
    pub fn resolve_checked(
        &self,
        value: &str,
        convention: &str,
        sources: Option<&[&str]>,
        fetch: &dyn Fetch,
        cache: &Cache,
        refresh: bool,
        now: Option<&str>,
    ) -> Result<Option<Substance>, ResolverError> {
        let substances =
            self.resolve_all(value, convention, sources, fetch, cache, refresh, now)?;
        let mut keys: Vec<&str> = Vec::new();
        for substance in &substances {
            if let Some(key) = substance.identifier_value("inchikey") {
                if !keys.contains(&key) {
                    keys.push(key);
                }
            }
        }
        if keys.len() > 1 {
            return Err(ResolverError::Conflict(format!(
                "sources disagree on substance identity (InChIKeys: {})",
                keys.join(" vs ")
            )));
        }
        Ok(substances.into_iter().next())
    }
}
