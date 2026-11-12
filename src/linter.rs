//! The linter (TODO.v2/17 slice 4) — the Rust mirror of
//! `AsciiChem::Linter`. The parser is total: it accepts any
//! syntactically valid input without judging chemistry correctness.
//! The linter is a separate opt-in pass that walks the model and
//! reports chemistry errors — unbalanced reactions, valence violations,
//! isotope nonsense, broken identifier annotations, structural identity
//! drift.
//!
//! # Registry
//!
//! [`Registry`] holds the checks as `(name, function)` pairs; the
//! built-in set is [`Registry::builtin`] and [`run`] executes it over a
//! parsed formula, returning diagnostics sorted errors-first (then
//! message, for stable output). Extending is [`Registry::new`] +
//! [`Registry::add`] — a caller-composed registry runs any check list.
//!
//! # Checks (the reference's thirteen, and where they stand here)
//!
//! balance, bracket balance, charge balance, element validation,
//! identifier consistency, identifier format, the engine-scaled
//! identity cross-check, isotope sanity, unclosed rings, and valence
//! are ported and gated (the reference's spec semantics mirrored in
//! [`checks`], the shared corpus's `linter-expectations` family in
//! `tests/corpus.rs`).
//!
//! Three reference checks have no Rust counterpart YET, by
//! construction rather than omission: CrystalSanityCheck,
//! SpectrumPeakCheck, and ZMatrixReferenceCheck validate model node
//! kinds (`Crystal`, `Spectrum`, `ZMatrix`) that the Rust model subset
//! does not carry — the model build declines those constructs with
//! `Unsupported` upstream, so no model tree reaching the linter can
//! contain them. They port together with the model subset expansion
//! (the corpus's linter-coverage extension, asciichem-tests#19, tracks
//! the fixtures).
//!
//! The bracket-balance check keeps its registry slot but shrinks to
//! the one property that can drift over Rust models (an unknown
//! bracket kind): the reference's open/close mismatch is
//! unrepresentable because the model derives both characters from the
//! kind.
//!
//! # Totality
//!
//! Every check reports, none panics. Unknown constructs are skipped
//! (conservative `None`), the identity cross-check turns engine
//! rejections into warnings, and diagnostics carry severity + message
//! + the node's diagnostic label ([`Diagnostic::label`]).

pub mod checks;
pub mod identifiers;
pub mod periodic_table;

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::model::Node;

/// How serious a diagnostic is. Ordering is errors first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    Error,
    Warning,
    Info,
}

impl Severity {
    /// The lowercase name the corpus expectations spell.
    pub fn as_str(&self) -> &'static str {
        match self {
            Severity::Error => "error",
            Severity::Warning => "warning",
            Severity::Info => "info",
        }
    }
}

impl std::fmt::Display for Severity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One lint finding. `label` is the offending node's diagnostic label
/// (the reference's `Diagnostic#to_s` context — e.g. `Atom(C)`);
/// `details` carries optional structured data behind the message (e.g.
/// balance's per-element `[reactants, products]` counts) for consumers
/// that want the numbers, not the prose.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub message: String,
    pub label: Option<String>,
    pub details: Option<BTreeMap<String, [i64; 2]>>,
}

impl Diagnostic {
    pub fn error(message: impl Into<String>, label: Option<String>) -> Self {
        Diagnostic::new(Severity::Error, message, label)
    }

    pub fn warning(message: impl Into<String>, label: Option<String>) -> Self {
        Diagnostic::new(Severity::Warning, message, label)
    }

    pub fn info(message: impl Into<String>, label: Option<String>) -> Self {
        Diagnostic::new(Severity::Info, message, label)
    }

    fn new(severity: Severity, message: impl Into<String>, label: Option<String>) -> Self {
        Diagnostic {
            severity,
            message: message.into(),
            label,
            details: None,
        }
    }

    /// Attaches the structured numbers behind the message.
    pub fn with_details(mut self, details: BTreeMap<String, [i64; 2]>) -> Self {
        self.details = Some(details);
        self
    }
}

impl std::fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.label {
            Some(label) => write!(f, "[{}] {label}: {}", self.severity, self.message),
            None => write!(f, "[{}] {}", self.severity, self.message),
        }
    }
}

/// A check: a pure walk over the model tree returning diagnostics.
pub type Check = fn(&[Node]) -> Vec<Diagnostic>;

/// The check registry. The built-in set mirrors the reference's
/// registry order (alphabetical by check name); iteration order is
/// insertion order.
#[derive(Debug, Clone, Default)]
pub struct Registry {
    checks: Vec<(&'static str, Check)>,
}

impl Registry {
    /// An empty registry, ready for caller-composed check lists.
    pub fn new() -> Self {
        Registry::default()
    }

    /// Registers a check under a name (first registration wins — the
    /// reference's dedup-by-name behaviour).
    pub fn add(&mut self, name: &'static str, check: Check) {
        if !self.names().any(|existing| existing == name) {
            self.checks.push((name, check));
        }
    }

    /// The registered check names, in registration order.
    pub fn names(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.checks.iter().map(|(name, _)| *name)
    }

    /// Runs every registered check against the formula, sorted by
    /// severity (errors first) then message for stable output.
    ///
    /// The reference also deduplicates by `(message, node identity)`;
    /// diagnostics here carry labels rather than node identity, and the
    /// built-in checks are disjoint per concern, so no duplicates
    /// arise — dropping same-message-same-label findings could
    /// silently merge distinct nodes (two flagged atoms), so the pass
    /// keeps every diagnostic.
    pub fn run(&self, formula: &[Node]) -> Vec<Diagnostic> {
        let mut diagnostics: Vec<Diagnostic> = self
            .checks
            .iter()
            .flat_map(|(_, check)| check(formula))
            .collect();
        diagnostics.sort_by(|a, b| {
            a.severity
                .cmp(&b.severity)
                .then_with(|| a.message.cmp(&b.message))
        });
        diagnostics
    }

    /// Whether any check produced an error-severity diagnostic.
    pub fn errors(&self, formula: &[Node]) -> bool {
        self.run(formula)
            .iter()
            .any(|diagnostic| diagnostic.severity == Severity::Error)
    }
}

/// The built-in registry: every check this module ports, in the
/// reference's registration order. Names mirror the reference's
/// registry entries (`Registry.names`).
pub fn builtin() -> &'static Registry {
    static CELL: OnceLock<Registry> = OnceLock::new();
    CELL.get_or_init(|| {
        let mut registry = Registry::new();
        registry.add("balance", checks::balance);
        registry.add("bracket_balance", checks::bracket_balance);
        registry.add("charge_balance", checks::charge_balance);
        registry.add("element_validation", checks::element_validation);
        registry.add("identifier_consistency", checks::identifier_consistency);
        registry.add("identifier_format", checks::identifier_format);
        registry.add("identity_cross_check", checks::identity_cross_check);
        registry.add("isotope_sanity", checks::isotope_sanity);
        registry.add("unclosed_ring", checks::unclosed_ring);
        registry.add("valence", checks::valence);
        registry
    })
}

/// Runs the built-in registry over a parsed formula — the linter's
/// front door. Empty output means no findings, not an empty formula.
pub fn run(formula: &[Node]) -> Vec<Diagnostic> {
    builtin().run(formula)
}

/// Convenience: whether any built-in check produced an error.
pub fn errors(formula: &[Node]) -> bool {
    builtin().errors(formula)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_nodes;

    fn lint(input: &str) -> Vec<Diagnostic> {
        run(&parse_nodes(input).expect("input parses"))
    }

    #[test]
    fn the_builtin_registry_carries_every_ported_check() {
        let names: Vec<_> = builtin().names().collect();
        assert_eq!(
            names,
            [
                "balance",
                "bracket_balance",
                "charge_balance",
                "element_validation",
                "identifier_consistency",
                "identifier_format",
                "identity_cross_check",
                "isotope_sanity",
                "unclosed_ring",
                "valence",
            ]
        );
    }

    #[test]
    fn registration_is_by_first_name() {
        let mut registry = Registry::new();
        registry.add("balance", checks::balance);
        registry.add("balance", checks::valence);
        assert_eq!(registry.names().collect::<Vec<_>>(), ["balance"]);
    }

    #[test]
    fn diagnostics_render_severity_then_label_then_message() {
        let diagnostic = Diagnostic::error("bad".to_string(), Some("Atom(C)".to_string()));
        assert_eq!(diagnostic.to_string(), "[error] Atom(C): bad");
        let diagnostic = Diagnostic::info("hello".to_string(), None);
        assert_eq!(diagnostic.to_string(), "[info] hello");
    }

    #[test]
    fn output_is_sorted_errors_first() {
        // `Hx_2O`: an element warning, a valence info, and (no engine
        // annotations) nothing else — warnings precede infos.
        let severities: Vec<_> = lint("Hx_2O")
            .into_iter()
            .map(|diagnostic| diagnostic.severity)
            .collect();
        let mut sorted = severities.clone();
        sorted.sort();
        assert_eq!(severities, sorted);
        assert_eq!(severities.first(), Some(&Severity::Warning));
        assert_eq!(severities.last(), Some(&Severity::Info));
    }

    #[test]
    fn errors_predicate_matches_severity() {
        assert!(errors(&parse_nodes("^5C").expect("parses")));
        assert!(!errors(&parse_nodes("H_2O").expect("parses")));
    }

    #[test]
    fn composed_registries_run_a_subset() {
        let mut registry = Registry::new();
        registry.add("element_validation", checks::element_validation);
        // The unbalanced reaction fires no element diagnostics, so a
        // registry holding only that check stays silent while the
        // builtin one reports.
        let formula = &parse_nodes("H_2 + O_2 -> H_2O").expect("parses");
        assert!(registry.run(formula).is_empty());
        assert!(!run(formula).is_empty());
    }
}
