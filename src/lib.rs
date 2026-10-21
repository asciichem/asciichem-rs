//! AsciiChem — Rust core (TODO.v2/13 item 7).
//!
//! Fourth implementation under the single-contract rule. The grammar
//! is NOT reimplemented: this crate consumes the compiled PARG artifact
//! of the AsciiChem syntax (`artifacts/asciichem.json`, checksummed and
//! verified by `parsanol` at load) and runs it on parsanol's portable
//! engine — the same artifact the reference implementation's `:parg`
//! engine runs.
//!
//! Current scope: parse acceptance/rejection plus the L1 text
//! round-trip contract over the shared corpus — the parse tree walks
//! through parsanol's parslet-compatible transform into the DSL shapes
//! (`value`), a typed model subset (`model`) built from them
//! (`build`), and the canonical text renderer. Constructs outside the
//! subset report `Unsupported` rather than guessing; renderers follow
//! the corpus contract in later slices.

use std::sync::OnceLock;

pub mod build;
pub mod model;
pub mod value;

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

/// Parses AsciiChem text and renders the canonical text form:
/// `parse_to_text(s) == s` for canonical `s` (the round-trip
/// contract). Constructs outside the model subset raise
/// `ParseError::Unsupported` — reported as such, never guessed.
pub fn parse_to_text(input: &str) -> Result<String, ParseError> {
    let source = input.trim().to_string();
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("asciichem-parse".into())
            .stack_size(PARSE_STACK_BYTES)
            .spawn_scoped(scope, || parse_to_text_inner(&source))
            .expect("spawn parser thread")
            .join()
            .expect("parser thread panicked")
    })
}

fn parse_to_text_inner(source: &str) -> Result<String, ParseError> {
    let value = parse_value(source)?;
    let nodes = build::build_formula(&value).map_err(|error| ParseError {
        message: error.to_string(),
    })?;
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
