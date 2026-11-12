//! Citation goldens (TODO.v2/17 phase 1, slice 6): the bibitem XML
//! for the same substance fixtures the Ruby reference implementation
//! cites, byte-compared against recorded captures of the Ruby
//! emission. The metanorma plugin's assertions are the contract
//! (bibitem type=dataset, docidentifier type=source, the accessed
//! date, keyword nesting, the provenance note); the goldens pin the
//! exact bytes.
//!
//! Offline by construction: substances are built here from the
//! resolver structs — the same shapes the Ruby adapters derive from
//! the recorded resolver fixtures (PubChem CID 2244 / CAS RN 50-78-2
//! aspirin) — with the retrieval time pinned so both sides emit
//! stable bytes. No network.
//!
//! Re-recording a golden (relaton or engine rides that deliberately
//! move the wire shape): resolve the substance through the Ruby gem
//! against its recorded fixture, pin `retrieved_at` to
//! `2026-01-05T09:30:00Z`, and capture `AsciiChem::Citation.to_xml`:
//!
//! ```ruby
//! $LOAD_PATH.unshift "lib"; require "asciichem"; require "asciichem/citation"
//! body = File.read("spec/fixtures/resolver/pubchem/aspirin.json")
//! fetch = Struct.new(:body).new(body)
//! def fetch.get(_url, _headers = {}) = body
//! substance = AsciiChem::Resolver[:pubchem].new.resolve(
//!   value: "50-78-2", convention: "cas", fetch: fetch,
//!   cache: AsciiChem::Resolver::Cache.new(dir: "/tmp/golden"))
//! substance.identifiers.each { |i| i.provenance.retrieved_at = "2026-01-05T09:30:00Z" }
//! File.binwrite("pubchem-aspirin.xml", AsciiChem::Citation.to_xml(substance))
//! ```

use asciichem::citation::{self, Error};
use asciichem::resolver::{Identifier, Provenance, Substance};
use std::fs;
use std::path::PathBuf;

/// The retrieval time pinned in every recorded golden.
const PINNED_RETRIEVED_AT: &str = "2026-01-05T09:30:00Z";

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/citation")
}

fn provenance(source: &str, attribution: Option<&str>) -> Provenance {
    Provenance {
        source: source.to_string(),
        retrieved_at: Some(PINNED_RETRIEVED_AT.to_string()),
        source_version: None,
        attribution: attribution.map(str::to_string),
    }
}

fn identifier(convention: &str, value: &str, provenance: Provenance) -> Identifier {
    Identifier {
        convention: convention.to_string(),
        value: value.to_string(),
        provenance: Some(provenance),
    }
}

fn assert_golden(name: &str, substance: &Substance) {
    let golden = fs::read(fixtures_dir().join(name)).unwrap_or_else(|error| {
        panic!("read golden fixture {name}: {error}");
    });
    let actual = citation::to_xml(substance).unwrap_or_else(|error| panic!("to_xml: {error}"));
    assert_eq!(
        actual.as_bytes(),
        golden,
        "byte mismatch against {name}:\n--- golden\n{}\n--- actual\n{}",
        String::from_utf8_lossy(&golden),
        actual
    );
}

/// The aspirin substance the PubChem adapter derives from its
/// recorded fixture (spec/fixtures/resolver/pubchem/aspirin.json).
fn pubchem_aspirin() -> Substance {
    let provenance = provenance(
        "pubchem",
        Some("PubChem, U.S. National Library of Medicine"),
    );
    Substance {
        preferred_name: Some("2-acetyloxybenzoic acid".to_string()),
        synonyms: vec![],
        formula: Some("C9H8O4".to_string()),
        molecular_weight: Some(180.16),
        identifiers: vec![
            identifier("pubchem-cid", "2244", provenance.clone()),
            identifier(
                "inchikey",
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                provenance.clone(),
            ),
            identifier(
                "inchi",
                "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                provenance.clone(),
            ),
            identifier("canonical-smiles", "CC(=O)OC1=CC=CC=C1C(=O)O", provenance),
        ],
        properties: vec![],
        structure: None,
    }
}

/// The aspirin substance the Common Chemistry adapter derives from its
/// recorded fixture (spec/fixtures/resolver/common_chemistry/aspirin.json).
fn common_chemistry_aspirin() -> Substance {
    let provenance = provenance(
        "common_chemistry",
        Some("CAS Common Chemistry (CC BY-NC 4.0)"),
    );
    Substance {
        preferred_name: Some("Aspirin".to_string()),
        synonyms: vec![],
        formula: Some("C9H8O4".to_string()),
        molecular_weight: None,
        identifiers: vec![
            identifier("cas", "50-78-2", provenance.clone()),
            identifier(
                "inchikey",
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                provenance.clone(),
            ),
            identifier(
                "canonical-smiles",
                "CC(=O)OC1=CC=CC=C1C(=O)O",
                provenance.clone(),
            ),
            identifier(
                "inchi",
                "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                provenance.clone(),
            ),
            identifier("iupac-name", "2-acetyloxybenzoic acid", provenance),
        ],
        properties: vec![],
        structure: None,
    }
}

#[test]
fn pubchem_aspirin_matches_the_ruby_bytes() {
    assert_golden("pubchem-aspirin.xml", &pubchem_aspirin());
}

#[test]
fn common_chemistry_aspirin_matches_the_ruby_bytes() {
    assert_golden("common-chemistry-aspirin.xml", &common_chemistry_aspirin());
}

/// The default-profile path the Ruby specs drive with a bare
/// Substance: provenance with a version and an attribution travels as
/// the provenance note.
#[test]
fn provenance_version_and_attribution_match_the_ruby_bytes() {
    let mut record = provenance("example", Some("Example Data Source"));
    record.source_version = Some("2026.11".to_string());
    let substance = Substance {
        preferred_name: None,
        synonyms: vec![],
        formula: None,
        molecular_weight: None,
        identifiers: vec![identifier("cas", "50-78-2", record)],
        properties: vec![],
        structure: None,
    };
    assert_golden("example-provenance.xml", &substance);
}

/// A source-only provenance record carries no note (it would
/// duplicate the publisher and docidentifier).
#[test]
fn source_only_provenance_matches_the_ruby_bytes() {
    let substance = Substance {
        preferred_name: None,
        synonyms: vec![],
        formula: None,
        molecular_weight: None,
        identifiers: vec![identifier("cas", "50-78-2", provenance("example", None))],
        properties: vec![],
        structure: None,
    };
    assert_golden("example-unattributed.xml", &substance);
}

#[test]
fn a_substance_without_provenance_is_rejected() {
    let bare = Substance {
        preferred_name: None,
        synonyms: vec![],
        formula: None,
        molecular_weight: None,
        identifiers: vec![Identifier {
            convention: "cas".to_string(),
            value: "50-78-2".to_string(),
            provenance: None,
        }],
        properties: vec![],
        structure: None,
    };
    assert_eq!(citation::bibitem(&bare), Err(Error::NoProvenance));
}

/// The plugin-consumer shape, read back through the engine's XPath:
/// the InChIKey keyword, the source docidentifier, the accessed date,
/// and the provenance note are all where the metanorma plugin looks
/// for them.
#[test]
fn the_emission_carries_the_plugin_consumed_shape() {
    let xml = citation::to_xml(&pubchem_aspirin()).expect("to_xml");
    let document = leptris::Document::parse(xml.as_bytes()).expect("parse");
    let keyword = document
        .xpath(
            "//*[local-name()='keyword']/*[local-name()='vocab'][starts-with(text(),'inchikey=')]",
        )
        .expect("xpath");
    assert_eq!(keyword.to_string(), "inchikey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N");

    let docidentifier = document
        .xpath("//*[local-name()='docidentifier'][@type='pubchem']")
        .expect("xpath");
    assert_eq!(docidentifier.to_string(), "PubChem CID 2244");

    let date = document
        .xpath("//*[local-name()='date'][@type='accessed']/*[local-name()='on']")
        .expect("xpath");
    assert_eq!(date.to_string(), "2026-01-05");

    let note = document
        .xpath("//*[local-name()='note'][@type='provenance']")
        .expect("xpath");
    assert_eq!(
        note.to_string(),
        "Retrieved from pubchem; PubChem, U.S. National Library of Medicine"
    );

    let root = document.root().expect("root");
    assert_eq!(root.attribute("type"), Some("dataset"));
    assert_eq!(root.attribute("id"), Some("PubChemCID2244"));
}
