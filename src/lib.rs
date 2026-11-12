//! AsciiChem — Rust core (TODO.v2/13 item 7).
//!
//! Fourth implementation under the single-contract rule. The grammar
//! is NOT reimplemented: this crate consumes the compiled PARG artifact
//! of the AsciiChem syntax (`artifacts/asciichem.json`, checksummed and
//! verified by `parsanol` at load) and runs it on parsanol's portable
//! engine — the same artifact the reference implementation's `:parg`
//! engine runs.
//!
//! Current scope: parse acceptance/rejection, the L0 wire round-trip,
//! the L1 text round-trip, and L2 MathML goldens over the shared
//! corpus — the parse tree walks through parsanol's parslet-compatible
//! transform into the DSL shapes (`value`), a typed model subset
//! (`model`) built from them (`build`), with the canonical text
//! renderer on the model, the MathML renderer (`mathml`) byte-exact
//! against the corpus L2 goldens, the inline HTML renderer (`html`)
//! and the mhchem LaTeX renderer (`latex`) matching the reference's
//! formatter output, and the canonical wire JSON (`wire`)
//! field-for-field with the shared v1 schemas. The structure
//! interchange layer (`structure`, `smiles`, `molfile`, `perception`)
//! ingests and emits SMILES and molfile V2000 with the shared
//! perception semantics (implicit hydrogens, Hückel aromatic rings,
//! stereo carry-through), gated by the corpus's SMILES, molfile, and
//! perception families. The identity layer (`identity`) derives
//! InChI/AuxInfo/InChIKey in-process through the chemomer dependency.
//! The linter ([`crate::linter`]) walks the model and reports
//! chemistry errors — balance, charge balance, valence, isotope
//! sanity, element validation, unclosed rings, identifier
//! format/consistency, and the identity cross-check through the same
//! in-process engine. The resolver ([`crate::resolver`]) turns
//! identifiers into substance records across PubChem, OPSIN and CAS
//! Common Chemistry — TTL-cached in the user cache dir, provenance-
//! stamped with each source's required attribution, per-source
//! isolated, and serialized through the wire substance-record form.
//! Constructs outside the subset report
//! `Unsupported` rather than guessing.
//!
//! The CML round-trip ([`crate::cml`]) carries the model through the
//! canonical CML wire with the aci: extension channels (atom
//! attributes, top-level elements, molecule groups), parsed and
//! serialized entirely by the linked leptris engine.
//! The citation layer ([`crate::citation`]) turns a resolved substance
//! ([`crate::resolver::Substance`]) into the dataset-type Relaton
//! bibitem XML — one per (substance, source), built and serialized
//! through the leptris engine, byte-exact with the Ruby emission.
//! The CLI ([`crate::cli`]) drives all of it from the `asciichem`
//! binary — the eleven commands the Ruby CLI ships, one dispatch, the
//! contract (stdout shapes + exit codes) testable offline.

use std::sync::OnceLock;

pub mod build;
pub mod citation;
pub mod cli;
pub mod cml;
pub mod html;
pub mod identifiers;
pub mod identity;
pub mod latex;
pub mod layout;
pub mod linter;
pub mod mathml;
pub mod model;
pub mod molfile;
pub mod perception;
pub mod resolver;
pub mod smiles;
pub mod structure;
pub mod svg;
pub mod value;
pub mod wire;

use parsanol::portable::grammar::Grammar;
use parsanol::{AstArena, PargArtifact, PortableParser};

const ARTIFACT: &str = include_str!("../artifacts/asciichem.json");

/// A parse rejection with a human-readable message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// What went wrong (load failure, position, or short capture).
    pub message: String,
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for ParseError {}

pub fn grammar() -> Result<&'static Grammar, ParseError> {
    static CELL: OnceLock<Result<Grammar, String>> = OnceLock::new();
    let loaded = CELL.get_or_init(|| {
        let artifact = PargArtifact::from_json(ARTIFACT)
            .map_err(|error| format!("PARG artifact load failed: {error:?}"))?;
        artifact
            .grammar("formula")
            .map_err(|error| format!("PARG entry extraction failed: {error:?}"))
    });
    loaded.as_ref().map_err(|message| ParseError {
        message: message.clone(),
    })
}

/// Stack reserved for the recursive PEG engine. Deeply nested input
/// (the corpus fuzz cases) exceeds a default thread stack, so parsing
/// runs on a dedicated thread with a generous stack — reserved
/// virtually, committed only as used.
const PARSE_STACK_BYTES: usize = 64 * 1024 * 1024;

/// A persistent parse worker: one thread with the deep stack, spawned
/// once, reused for every parse. Spawning a fresh 64 MiB-reserved
/// thread per call cost ~300 µs — two orders of magnitude above the
/// parse itself.
mod worker {
    use super::ParseError;
    use std::sync::mpsc::{channel, Sender};
    use std::sync::OnceLock;

    type Job = Box<dyn FnOnce() + Send>;

    fn worker() -> &'static Sender<Job> {
        static SENDER: OnceLock<Sender<Job>> = OnceLock::new();
        SENDER.get_or_init(|| {
            let (sender, receiver) = channel::<Job>();
            std::thread::Builder::new()
                .name("asciichem-parse".into())
                .stack_size(super::PARSE_STACK_BYTES)
                .spawn(move || {
                    while let Ok(job) = receiver.recv() {
                        job();
                    }
                })
                .expect("spawn persistent parser thread");
            sender
        })
    }

    /// Runs `f` on the deep-stack worker thread, returning its result.
    pub fn run<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
        let (result_tx, result_rx) = channel();
        worker()
            .send(Box::new(move || {
                result_tx.send(f()).ok();
            }))
            .expect("parser worker alive");
        result_rx.recv().expect("parser worker returned a result")
    }

    /// Convenience for parse results: mirrors the caller's error type.
    pub fn parse<T: Send + 'static>(
        f: impl FnOnce() -> Result<T, ParseError> + Send + 'static,
    ) -> Result<T, ParseError> {
        run(f)
    }
}

/// Parses AsciiChem text into the model subset's formula nodes —
/// the shared front half of every renderer.
/// Inputs at or below this length parse on the CALLER's stack: PEG
/// recursion depth is bounded by input length, so a short input
/// cannot exhaust a default 8 MB stack — and skipping the deep-stack
/// worker hop is worth ~200 µs per parse. Longer inputs go to the
/// worker (the corpus fuzz cases genuinely need the 64 MiB stack).
const INLINE_PARSE_MAX: usize = 16 * 1024;

pub fn parse_nodes(input: &str) -> Result<Vec<model::Node>, ParseError> {
    let source = input.trim();
    if source.len() <= INLINE_PARSE_MAX {
        return parse_nodes_direct(source);
    }
    let owned = source.to_string();
    worker::parse(move || parse_nodes_direct(&owned))
}

/// The parse itself — runs ON the deep-stack worker. Internal: never
/// call from inside another worker job (the worker is single-threaded;
/// nesting self-deadlocks).
pub(crate) fn parse_nodes_direct(source: &str) -> Result<Vec<model::Node>, ParseError> {
    let value = parse_value(source)?;
    build::build_formula(&value).map_err(|error| ParseError {
        message: error.to_string(),
    })
}

pub fn parse_to_text(input: &str) -> Result<String, ParseError> {
    let source = input.trim();
    if source.len() <= INLINE_PARSE_MAX {
        return parse_to_text_inner(source);
    }
    let owned = source.to_string();
    worker::parse(move || parse_to_text_inner(&owned))
}

fn parse_to_text_inner(source: &str) -> Result<String, ParseError> {
    let nodes = parse_nodes_direct(source)?;
    let rendered: Vec<String> = nodes.iter().map(|n| n.to_string()).collect();
    Ok(rendered.join(" "))
}

fn parse_value(source: &str) -> Result<parsanol::portable::transform::Value, ParseError> {
    let grammar = grammar()?;
    let mut arena = AstArena::for_input(source.len());
    let mut parser = PortableParser::new(grammar, source, &mut arena);
    let result = parser.parse_with_end_pos().map_err(|error| ParseError {
        message: error.to_string(),
    })?;
    if result.end_pos < source.len() {
        return Err(ParseError {
            message: format!(
                "input not fully captured: consumed {} of {} bytes",
                result.end_pos,
                source.len()
            ),
        });
    }
    let shaped = parsanol::portable::parslet_transform::to_parslet_compatible(
        &result.value,
        &mut arena,
        source,
    );
    let value = parsanol::portable::transform::ast_to_value(&shaped, &arena, source);
    Ok(value::normalise(&value))
}

/// Parses AsciiChem text through the PARG grammar. `Ok(())` when the
/// whole trimmed input is accepted; a rejection or a short capture is
/// an error (never a silently partial parse).
pub fn parse(input: &str) -> Result<(), ParseError> {
    let source = input.trim().to_string();
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("asciichem-parse".into())
            .stack_size(PARSE_STACK_BYTES)
            .spawn_scoped(scope, || parse_inner(&source))
            .expect("spawn parser thread")
            .join()
            .expect("parser thread panicked")
    })
}

fn parse_inner(source: &str) -> Result<(), ParseError> {
    let result = grammar()?
        .parse_with_pos(source)
        .map_err(|error| ParseError {
            message: error.to_string(),
        })?;
    if result.end_pos < source.len() {
        return Err(ParseError {
            message: format!(
                "input not fully captured: consumed {} of {} bytes",
                result.end_pos,
                source.len()
            ),
        });
    }
    Ok(())
}
