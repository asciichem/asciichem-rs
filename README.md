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

Scaffold. Implementation follows the
[corpus](https://github.com/asciichem/asciichem-tests) fixture
contract (parse/reject/round-trip), then output formats per the
reference implementation's coverage.

## Development

```sh
cargo test
```
