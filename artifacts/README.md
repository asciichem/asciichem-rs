# Vendored PARG artifact

`asciichem.json` is the compiled PARG grammar of the AsciiChem syntax,
vendored from the contract repository
[asciichem-grammar](https://github.com/asciichem/asciichem-grammar) (the same
artifact the Ruby reference ships at
`lib/asciichem/engine/parg/asciichem.json`). PARG artifacts are
checksummed; `parsanol` refuses a mismatched envelope, so a stale copy
fails loudly rather than parsing a different language.
