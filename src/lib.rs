//! AsciiChem — Rust core (TODO.v2/13 item 7).
//!
//! Fourth implementation under the single-contract rule. The grammar
//! is NOT reimplemented: this crate consumes the compiled PARG artifact
//! of the AsciiChem syntax (`artifacts/asciichem.json`, checksummed and
//! verified by `parsanol` at load) and runs it on parsanol's portable
//! engine — the same artifact the reference implementation's `:parg`
//! engine runs.
//!
//! Current scope: parse acceptance/rejection over the shared corpus.
//! The model, transform, and renderers follow the corpus contract in
//! later slices.

use std::sync::OnceLock;

use parsanol::portable::grammar::Grammar;
use parsanol::PargArtifact;

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

fn grammar() -> Result<&'static Grammar, ParseError> {
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
