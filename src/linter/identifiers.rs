//! Offline format validation for chemical identifier conventions —
//! the linter's surface over the crate-level validators
//! ([`crate::identifiers`], the Rust mirror of `AsciiChem::Identifiers`).
//!
//! Validation is purely lexical/checksum based. It never touches the
//! network and never decides whether an identifier is *assigned* to a
//! real substance — that is resolution (TODO.v2 07 Layer 2), not
//! linting.

/// The conventions with a registered offline validator. Unknown
/// conventions are not errors — they are simply not checkable offline
/// (e.g. `iupac`, `cid`, `chebi`).
pub fn known(convention: &str) -> bool {
    matches!(convention, "cas" | "inchi" | "inchikey" | "smiles")
}

/// `None` when the value is valid (or the convention is unknown);
/// `Some(reason)` otherwise. The reason is the human-facing message
/// the identifier-format check reports.
pub fn diagnostic(convention: &str, value: &str) -> Option<String> {
    crate::identifiers::diagnostic(convention, value)
}

#[cfg(test)]
mod tests {
    use super::{diagnostic, known};
    use crate::identifiers::{inchi, smiles};

    #[test]
    fn cas_check_digit_matches_the_weighted_sum() {
        assert_eq!(diagnostic("cas", "7732-18-5"), None);
        assert_eq!(diagnostic("cas", "50-78-2"), None);
        assert_eq!(
            diagnostic("cas", "74-82-9"),
            Some("CAS RN \"74-82-9\": check digit 9 does not match weighted sum (8)".into())
        );
    }

    #[test]
    fn cas_format_rejects_wrong_shapes() {
        let expected = "CAS RN \"50782\" is not in DDDDDDD-DD-D form \
                        (2–7 digits, 2 digits, 1 check digit)";
        assert_eq!(diagnostic("cas", "50782"), Some(expected.into()));
        assert!(diagnostic("cas", "773-218-5").is_some());
        assert!(diagnostic("cas", "ab-32-1").is_some());
    }

    #[test]
    fn inchi_checks_the_prefix_and_the_formula_layer() {
        let aspirin = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)";
        assert_eq!(diagnostic("inchi", aspirin), None);
        assert_eq!(
            diagnostic("inchi", "InChI=1S/Cx4"),
            Some("InChI \"InChI=1S/Cx4\": unknown element symbol \"Cx\" in formula layer".into())
        );
        assert_eq!(
            diagnostic("inchi", "garbage"),
            Some("InChI \"garbage\" must start with \"InChI=1/\" or \"InChI=1S/\"".into())
        );
        assert_eq!(
            diagnostic("inchi", "InChI=1S/"),
            Some("InChI \"InChI=1S/\": formula layer is empty".into())
        );
    }

    #[test]
    fn inchi_formula_counts_single_source() {
        let counts = inchi::formula_counts("InChI=1S/H2O/h1H2").expect("analysable");
        assert_eq!(counts.get("H"), Some(&2));
        assert_eq!(counts.get("O"), Some(&1));
        assert_eq!(inchi::formula_counts("garbage"), None);
    }

    #[test]
    fn inchikey_checks_block_shape() {
        assert_eq!(diagnostic("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), None);
        assert_eq!(diagnostic("inchikey", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"), None);
        assert!(diagnostic("inchikey", "bsynrymutxbxsq-uhfffaoysa-n").is_some());
        assert!(diagnostic("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA").is_some());
    }

    #[test]
    fn smiles_flags_corruption() {
        assert_eq!(diagnostic("smiles", "CCO"), None);
        assert_eq!(diagnostic("smiles", "c1ccccc1"), None);
        assert_eq!(diagnostic("smiles", "CC(=O)Oc1ccccc1C(=O)O"), None);
        assert_eq!(
            diagnostic("smiles", ""),
            Some("SMILES must not be empty".into())
        );
        assert_eq!(
            diagnostic("smiles", "C C"),
            Some("SMILES \"C C\" must not contain whitespace".into())
        );
        assert_eq!(
            diagnostic("smiles", "C(C"),
            Some("SMILES \"C(C\": parentheses are not balanced".into())
        );
        assert_eq!(
            diagnostic("smiles", "C[CH"),
            Some("SMILES \"C[CH\": square brackets are not balanced".into())
        );
        assert_eq!(
            diagnostic("smiles", "C1-CC"),
            Some(
                "SMILES \"C1-CC\": ring-bond digit 1 appears an odd number of times \
                 (ring closures must pair)"
                    .into()
            )
        );
        assert_eq!(
            diagnostic("smiles", "Qq"),
            Some("SMILES \"Qq\": contains a malformed or unknown element token".into())
        );
    }

    #[test]
    fn smiles_element_set_normalises_aromatics() {
        let set = smiles::element_set("CC(=O)Oc1ccccc1C(=O)O").expect("analysable");
        assert_eq!(set, ["C", "O"].into_iter().collect::<Vec<_>>());
        let bracketed = smiles::element_set("[NH4+]CCl").expect("analysable");
        assert_eq!(bracketed, ["C", "Cl", "N"]);
        assert_eq!(smiles::element_set("C1-CC"), None);
        assert_eq!(smiles::element_set("C[*]C"), None);
    }

    #[test]
    fn unknown_conventions_are_not_checkable_offline() {
        assert!(!known("iupac"));
        assert_eq!(diagnostic("iupac", "anything"), None);
        assert_eq!(diagnostic("cid", "2244"), None);
    }
}
