# asciichem-rs

Rust implementation of **AsciiChem** — the ASCII syntax for chemistry
(atoms with isotopes and charges, bonds, reactions, electron
configurations) parsed into a semantic model.

Fourth implementation in the ecosystem, under the single-contract
rule: it parses the same syntax defined by
[asciichem-parg](https://github.com/asciichem/asciichem-parg), builds
the same model shape, and must pass the same shared conformance
corpus ([asciichem-tests](https://github.com/asciichem/asciichem-tests))
that [asciichem-ruby](https://github.com/asciichem/asciichem-ruby)
(the reference),
[asciichem-ts](https://github.com/asciichem/asciichem-ts), and
[asciichem-py](https://github.com/asciichem/asciichem-py) pass.

## Status

**Parse contract at 100%; L0 wire round-trip at 100%; L1 text
round-trip at 100%; L2 MathML goldens at 100%** — the crate consumes
the compiled PARG artifact of
the AsciiChem grammar (`artifacts/asciichem.json`, checksummed and
verified by `parsanol`) and runs it on parsanol's portable engine; the
shared corpus's acceptance/rejection contract passes in CI, every
round-trip fixture survives the canonical wire JSON (field-for-field
with the shared asciichem-model v1 schemas), renders back through the
canonical text form byte-identically (`parse_to_text(s) == s`), and
every corpus MathML golden renders byte-exactly (`mathml::render`).
Deeply nested input runs on a dedicated large-stack thread.

The semantic model subset (atoms, bonds, groups, molecules with
annotations, reactions with conditions, cascades, quoted text,
embedded math, mechanisms) mirrors the reference implementation's
TransformRules; constructs outside the subset report `not yet
supported` instead of guessing.

Next slices follow the corpus: the SVG renderers.

## Development

```sh
cargo test
```
