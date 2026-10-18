//! AsciiChem: ASCII syntax for chemistry parsed into a semantic model.
//!
//! Rust port of the [asciichem](https://github.com/asciichem/asciichem-ruby)
//! reference implementation, conforming to the shared
//! [corpus](https://github.com/asciichem/asciichem-tests).

/// Placeholder module marker; the parser and model land with the
/// corpus-conformance implementation.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_semver() {
        assert!(super::VERSION.split('.').count() == 3);
    }
}
