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
round-trip at 100%; L2 MathML goldens and SVG at 100%; SMILES +
perception + molfile interchange at 100%** — the crate
consumes the compiled PARG artifact of
the AsciiChem grammar (`artifacts/asciichem.json`, checksummed and
verified by `parsanol`) and runs it on parsanol's portable engine; the
shared corpus's acceptance/rejection contract passes in CI, every
round-trip fixture survives the canonical wire JSON (field-for-field
with the shared asciichem-model v1 schemas), renders back through the
canonical text form byte-identically (`parse_to_text(s) == s`), and
every corpus MathML golden renders byte-exactly (`mathml::render`),
and the SVG renderers (`svg::render` structural with linear fallback,
`svg::render_linear`) are byte-exact with the TypeScript
implementation's renderers over a committed reference snapshot.
The SMILES reader/writer (`smiles::parse_smiles`/`smiles::write_smiles`)
passes the corpus's SMILES family (parse/reject, exact round trip,
canonical aromatic form) and the perception family (implicit hydrogens
+ Hückel aromatic rings, declines included) identically to the
TypeScript and Python implementations. The molfile V2000
reader/writer (`molfile::parse_molfile`/`molfile::write_molfile`)
passes the corpus's molfile family — parse/reject, atom/bond counts,
graph-stable round trips, and the stereo-code gate (codes 1/6 →
wedge/hash; chirality crosses as a wedge on the stereocentre's first
bond) — with deterministic re-emission (the written form is a byte
fixed point). The CML round-trip (`cml::to_cml`/`cml::parse_cml`)
passes the corpus's cml family through the three-way contract
(parse → CML → parse → canonical text) with the aci: extension
channels — atom attributes (oxidation state, Lewis markers, ring
closures), top-level elements (electron configuration, embedded math,
quoted text), and molecule groups — with every byte of XML parsing
and serialization performed by the linked
[leptris](https://github.com/leptris/leptris) engine. Deeply nested
input runs on a dedicated large-stack thread.

The semantic model subset (atoms, bonds, groups, molecules with
annotations, reactions with conditions, cascades, quoted text,
embedded math, mechanisms) mirrors the reference implementation's
TransformRules; constructs outside the subset report `not yet
supported` instead of guessing. The 2D layout is the deterministic
chemistry-shaped placer (rings as regular polygons, fused rings
anchored on shared atoms, chains as 120-degree zigzags) shared with
the TypeScript implementation.

Identity is derived IN-PROCESS through the
[chemomer](https://github.com/chemomer/chemomer) crate — the native
Rust InChI engine, byte-exact against the IUPAC oracle, consumed here
as a version-pinned dependency with no FFI seam (the Ruby gem reaches
the same engine through its C ABI, TS/PY through wasm/PyO3).
`identity::derive` pipes the model through a V2000 molfile to the
engine; `identity::derive_from_molfile` is the raw seam (molfile in →
InChI, AuxInfo, InChIKey out). CI gates the oracle corpus's full
recorded identity family (chemomer/inchi-tests, pinned tag) byte-exact
through the dependency.

## Conformance levels

Mapped onto the ecosystem's formal levels — see
[the levels document](https://www.asciichem.org/guides/conformance-levels/):

| Level | Claim |
|---|---|
| L0 notation | **yes** — the parse contract at 100% against the PARG artifact; canonical text round-trip |
| L1 semantics | **yes** — wire round-trip at 100%, field-for-field with the asciichem-model v1 schemas; the model subset mirrors the reference TransformRules — out-of-subset constructs report `not yet supported`, never guessed |
| L2 render | **yes** — MathML goldens at 100%; SVG byte-exact with the TypeScript implementation |
| L3 stereo | **yes** — tetrahedral chirality (`@`/`@@` → TH1/TH2), directional bonds (`/`/`\`), wedge/hash depiction; Hückel aromaticity + implicit-hydrogen perception with the declines; parity-correct canonical SMILES; the molfile stereo-code gate (V2000 codes 1/6 → wedge/hash, chirality crossing as a wedge) |
| L4 identity | **partial** — InChI/AuxInfo/InChIKey derived in-process through the chemomer dependency, byte-exact against the recorded oracle identity family; identifier cross-validation and citation-grade resolution are not yet ported |

## Development

```sh
cargo test
```

The CML round-trip links [libleptris](https://github.com/leptris/leptris)
(the XML engine; the same build leptris-rs's CI runs). Build it once
and point the crate at it:

```sh
cmake -B /path/to/leptris/build-shared -S /path/to/leptris \
  -DCMAKE_BUILD_TYPE=Release \
  -DLEPTRIS_BUILD_SHARED=ON -DLEPTRIS_BUILD_STATIC=OFF \
  -DBUILD_TESTING=OFF -DLEPTRIS_BUILD_CLI=OFF \
  -DLEPTRIS_BUILD_BENCHMARKS=OFF -DLEPTRIS_ENABLE_UTF8PROC=ON
cmake --build /path/to/leptris/build-shared -j
LEPTRIS_LIB_PATH=/path/to/leptris/build-shared/src \
LD_LIBRARY_PATH=/path/to/leptris/build-shared/src \
  cargo test
```

The identity oracle gate needs the chemomer corpus (skipped when
absent): clone [chemomer/inchi-tests](https://github.com/chemomer/inchi-tests)
at the pin as a sibling (`../inchi-tests`) or point `CHEMOMER_CORPUS`
at its `corpus/fixtures`, and run it in release profile — canonical
ranking is debug-bound:

```sh
CHEMOMER_CORPUS=../inchi-tests/corpus/fixtures cargo test --release --test identity_corpus
```
