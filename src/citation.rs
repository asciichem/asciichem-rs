//! The citation track (TODO.v2/17 phase 1, slice 6): one dataset-type
//! Relaton bibitem per (substance, source) — two databases agreeing on
//! a substance are still two different documents. Semantics mirror the
//! Ruby `AsciiChem::Citation` contract exactly: source profiles
//! ([`PROFILES`], plus the fallback default), the InChIKey-anchored
//! wire the metanorma plugin consumes, the accessed date, every
//! identifier carried as a keyword for cross-checking, and the
//! provenance record as a bibliographic note (licence terms — CAS
//! Common Chemistry's CC BY-NC 4.0 line — must survive into the
//! rendered document).
//!
//! The bibitem XML is BUILT and SERIALIZED through the leptris engine
//! (the mandated serialization layer): [`relaton`] assembles a typed
//! [`relaton::Bibitem`] payload and emits it via the engine's
//! programmatic construction API ([`engine`]) — no hand-rolled XML
//! strings; escaping, well-formedness, and the serialization bytes all
//! come from libleptris. The goldens in `tests/citation_golden.rs` are
//! byte-captures of the Ruby emission for the same substances.

pub mod engine;
pub mod relaton;

use crate::resolver::Substance;
use relaton::Bibitem;
use std::time::{SystemTime, UNIX_EPOCH};

/// Source-specific citation profile: the publisher string, the source
/// link shape, and the document identifier. Adding a source is one
/// entry in [`PROFILES`] (OCP); unknown sources fall back to
/// [`DEFAULT_PROFILE`], whose fields derive from the substance itself.
struct Profile {
    publisher: Option<&'static str>,
    link_for: fn(&Substance) -> Option<String>,
    identifier_for: fn(&Substance) -> Option<String>,
}

/// The known sources. Keys are the resolver's source names.
const PROFILES: &[(&str, Profile)] = &[
    (
        "pubchem",
        Profile {
            publisher: Some("PubChem, U.S. National Library of Medicine"),
            link_for: |substance| {
                substance
                    .identifier_value("pubchem-cid")
                    .map(|cid| format!("https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"))
            },
            identifier_for: |substance| {
                substance
                    .identifier_value("pubchem-cid")
                    .map(|cid| format!("PubChem CID {cid}"))
            },
        },
    ),
    (
        "common_chemistry",
        Profile {
            publisher: Some("CAS Common Chemistry"),
            link_for: |substance| {
                substance
                    .identifier_value("cas")
                    .map(|cas| format!("https://commonchemistry.cas.org/detail?cas_rn={cas}"))
            },
            identifier_for: |substance| {
                substance
                    .identifier_value("cas")
                    .map(|cas| format!("CAS RN {cas}"))
            },
        },
    ),
];

/// The fallback profile for sources without a dedicated entry: no
/// publisher override, no link, and the document identifier is the
/// substance's first identifier.
const DEFAULT_PROFILE: Profile = Profile {
    publisher: None,
    link_for: |_substance| None,
    identifier_for: |substance| {
        substance
            .identifiers
            .first()
            .map(|first| format!("{}: {}", first.convention, first.value))
    },
};

fn profile_for(source: &str) -> &'static Profile {
    PROFILES
        .iter()
        .find(|(name, _)| *name == source)
        .map(|(_, profile)| profile)
        .unwrap_or(&DEFAULT_PROFILE)
}

/// The citation track speaks one error type, mirroring the Ruby
/// `AsciiChem::Error` shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The substance carries no provenance — it is hand-built, not
    /// resolved. Resolve it first (the resolver, TODO.v2/17 slice 5).
    NoProvenance,
    /// The leptris engine rejected a build or serialization step.
    Engine(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::NoProvenance => {
                write!(
                    f,
                    "substance has no provenance - resolve it first (the resolver)"
                )
            }
            Error::Engine(message) => write!(f, "leptris engine error: {message}"),
        }
    }
}

impl std::error::Error for Error {}

/// Builds a dataset-type bibitem payload for the substance as cited
/// from its recorded source. Raises [`Error::NoProvenance`] when the
/// substance carries no provenance (hand-built, not resolved).
pub fn bibitem(substance: &Substance) -> Result<Bibitem, Error> {
    let provenance = substance.provenance().ok_or(Error::NoProvenance)?;
    let profile = profile_for(&provenance.source);
    let publisher = profile.publisher.unwrap_or(&provenance.source).to_string();
    let docid = (profile.identifier_for)(substance)
        .unwrap_or_else(|| format!("{} substance", provenance.source));

    Ok(Bibitem {
        bibitem_type: "dataset".to_string(),
        title: format!("{} - {} substance record", title_base(substance), publisher),
        docidentifier_type: provenance.source.clone(),
        docidentifier: docid,
        publisher,
        accessed_on: accessed_on(provenance),
        link: (profile.link_for)(substance),
        keywords: substance
            .identifiers
            .iter()
            .map(|identifier| format!("{}={}", identifier.convention, identifier.value))
            .collect(),
        note: provenance.note_text(),
    })
}

/// Convenience: the bibitem XML (what a document pipeline embeds) —
/// the byte-exact Ruby emission shape, serialized through the leptris
/// engine.
pub fn to_xml(substance: &Substance) -> Result<String, Error> {
    bibitem(substance)?.to_xml()
}

/// The title base prefers the resolved name, then the registry keys
/// (CAS, then InChIKey), then the bare "Substance". The lookup
/// results go through `to_string`, so the chain compiles whether the
/// substance hands out owned strings or borrows (the slice-5
/// contract hands out borrows).
fn title_base(substance: &Substance) -> String {
    substance
        .preferred_name
        .clone()
        .or_else(|| {
            substance
                .identifier_value("cas")
                .map(|value| value.to_string())
        })
        .or_else(|| {
            substance
                .identifier_value("inchikey")
                .map(|value| value.to_string())
        })
        .unwrap_or_else(|| "Substance".to_string())
}

/// The accessed date: the provenance's retrieval date (first ten
/// characters — the ISO timestamp's `YYYY-MM-DD`), or today (UTC)
/// when the resolver recorded no retrieval time.
fn accessed_on(provenance: &crate::resolver::Provenance) -> String {
    match &provenance.retrieved_at {
        Some(timestamp) if timestamp.len() >= 10 => timestamp[..10].to_string(),
        _ => today_utc(),
    }
}

/// Today's UTC date as `YYYY-MM-DD`, from the system clock.
fn today_utc() -> String {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before the epoch")
        .as_secs();
    let (year, month, day) = civil_from_days((seconds / 86_400) as i64);
    format!("{year:04}-{month:02}-{day:02}")
}

/// Days since the epoch to a civil (year, month, day) — Howard
/// Hinnant's algorithm, the standard date-arithmetic reference.
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let day_of_era = (z - era * 146_097) as u64; // [0, 146096]
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era as i64 + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let mp = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * mp + 2) / 5 + 1) as u32;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if month <= 2 { year + 1 } else { year }, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::{Identifier, Provenance};

    fn provenance(source: &str) -> Provenance {
        Provenance {
            source: source.to_string(),
            retrieved_at: Some("2026-01-05T09:30:00Z".to_string()),
            source_version: None,
            attribution: None,
        }
    }

    fn substance(source: &str, identifiers: &[(&str, &str)]) -> Substance {
        let provenance = provenance(source);
        Substance {
            preferred_name: None,
            synonyms: Vec::new(),
            formula: None,
            molecular_weight: None,
            identifiers: identifiers
                .iter()
                .map(|(convention, value)| Identifier {
                    value: (*value).to_string(),
                    convention: (*convention).to_string(),
                    provenance: Some(provenance.clone()),
                })
                .collect(),
            properties: Vec::new(),
            structure: None,
        }
    }

    #[test]
    fn builds_the_pubchem_payload() {
        let substance = substance(
            "pubchem",
            &[
                ("pubchem-cid", "2244"),
                ("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"),
            ],
        );
        let item = bibitem(&substance).expect("bibitem");
        assert_eq!(item.bibitem_type, "dataset");
        // The title base falls through the registry keys to the
        // InChIKey when no name or CAS is recorded.
        assert_eq!(
            item.title,
            "BSYNRYMUTXBXSQ-UHFFFAOYSA-N - PubChem, U.S. National Library of Medicine substance record"
        );
        assert_eq!(item.docidentifier_type, "pubchem");
        assert_eq!(item.docidentifier, "PubChem CID 2244");
        assert_eq!(item.accessed_on, "2026-01-05");
        assert_eq!(
            item.link,
            Some("https://pubchem.ncbi.nlm.nih.gov/compound/2244".to_string())
        );
        assert_eq!(
            item.keywords,
            ["pubchem-cid=2244", "inchikey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N"]
        );
    }

    #[test]
    fn prefers_the_name_for_the_title_base() {
        let mut substance = substance("pubchem", &[("pubchem-cid", "2244")]);
        substance.preferred_name = Some("2-acetyloxybenzoic acid".to_string());
        let item = bibitem(&substance).expect("bibitem");
        assert!(item.title.starts_with("2-acetyloxybenzoic acid - "));
    }

    #[test]
    fn falls_back_to_the_cas_then_inchikey_then_substance() {
        let cas = substance("example", &[("cas", "50-78-2")]);
        assert!(bibitem(&cas)
            .expect("bibitem")
            .title
            .starts_with("50-78-2 - "));

        let inchikey = substance("example", &[("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")]);
        assert!(bibitem(&inchikey)
            .expect("bibitem")
            .title
            .starts_with("BSYNRYMUTXBXSQ-UHFFFAOYSA-N - "));

        let nothing = substance("example", &[("note", "x")]);
        assert!(bibitem(&nothing)
            .expect("bibitem")
            .title
            .starts_with("Substance - "));
    }

    #[test]
    fn cites_common_chemistry_with_its_own_profile() {
        let substance = substance("common_chemistry", &[("cas", "50-78-2")]);
        let item = bibitem(&substance).expect("bibitem");
        assert_eq!(item.publisher, "CAS Common Chemistry");
        assert_eq!(item.docidentifier, "CAS RN 50-78-2");
        assert_eq!(
            item.link,
            Some("https://commonchemistry.cas.org/detail?cas_rn=50-78-2".to_string())
        );
    }

    #[test]
    fn unknown_sources_fall_back_to_the_default_profile() {
        let substance = substance("example", &[("cas", "50-78-2")]);
        let item = bibitem(&substance).expect("bibitem");
        // publisher falls back to the source name, docid to the first
        // identifier, no link.
        assert_eq!(item.publisher, "example");
        assert_eq!(item.docidentifier, "cas: 50-78-2");
        assert_eq!(item.docidentifier_type, "example");
        assert_eq!(item.link, None);
    }

    #[test]
    fn a_profile_without_its_identifier_still_cites() {
        // PubChem record with no CID recorded: docid falls back to the
        // "{source} substance" shape, no link.
        let substance = substance("pubchem", &[("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")]);
        let item = bibitem(&substance).expect("bibitem");
        assert_eq!(item.docidentifier, "pubchem substance");
        assert_eq!(item.link, None);
    }

    #[test]
    fn raises_for_substances_without_provenance() {
        let bare = Substance {
            preferred_name: None,
            synonyms: Vec::new(),
            formula: None,
            molecular_weight: None,
            identifiers: vec![Identifier {
                value: "50-78-2".to_string(),
                convention: "cas".to_string(),
                provenance: None,
            }],
            properties: Vec::new(),
            structure: None,
        };
        assert_eq!(bibitem(&bare), Err(Error::NoProvenance));
        assert_eq!(
            Error::NoProvenance.to_string(),
            "substance has no provenance - resolve it first (the resolver)"
        );
    }

    #[test]
    fn the_note_travels_from_the_provenance() {
        let mut provenance = provenance("example");
        provenance.source_version = Some("2026.11".to_string());
        provenance.attribution = Some("Example Data Source".to_string());
        let substance = Substance {
            identifiers: vec![Identifier {
                value: "50-78-2".to_string(),
                convention: "cas".to_string(),
                provenance: Some(provenance),
            }],
            ..substance("example", &[]).clone()
        };
        assert_eq!(
            bibitem(&substance).expect("bibitem").note,
            Some("Retrieved from example; version 2026.11; Example Data Source".to_string())
        );
    }

    #[test]
    fn a_retrieval_time_shorter_than_a_date_falls_back_to_today() {
        let mut provenance = provenance("example");
        provenance.retrieved_at = Some("bad".to_string());
        let substance = Substance {
            identifiers: vec![Identifier {
                value: "50-78-2".to_string(),
                convention: "cas".to_string(),
                provenance: Some(provenance),
            }],
            ..substance("example", &[]).clone()
        };
        let item = bibitem(&substance).expect("bibitem");
        assert_eq!(item.accessed_on, today_utc());
    }

    #[test]
    fn civil_days_land_on_known_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(19_723), (2024, 1, 1)); // 2024-01-01
        assert_eq!(civil_from_days(20_650), (2026, 7, 16)); // 2026-07-16
        assert_eq!(civil_from_days(-1), (1969, 12, 31));
        // The leap day itself.
        assert_eq!(civil_from_days(19_743), (2024, 1, 21));
    }

    #[test]
    fn today_utc_is_the_epoch_day_shape() {
        let today = today_utc();
        assert_eq!(today.len(), 10);
        assert_eq!(today.as_bytes()[4], b'-');
        assert_eq!(today.as_bytes()[7], b'-');
    }
}
