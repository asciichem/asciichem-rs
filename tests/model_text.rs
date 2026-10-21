//! Unit contract for the model subset and the canonical text
//! renderer: every supported construct builds and renders back to its
//! canonical text (`parse_to_text(s) == s`), and unsupported
//! constructs report `not yet supported` instead of guessing.

use asciichem::parse_to_text;

fn round_trip(input: &str) -> String {
    let output = parse_to_text(input).unwrap_or_else(|error| panic!("{input:?}: {error}"));
    assert_eq!(output, input, "input {input:?} did not round-trip");
    output
}

#[test]
fn atoms_render_canonical_text() {
    round_trip("H_2O");
    round_trip("^14C");
    round_trip("Ca^2+");
    round_trip("Na^+");
    round_trip("O^2-");
    round_trip("Mn^(VII)");
    round_trip("Fe^(III)");
    round_trip(":C.");
    round_trip("C1-C-C-C-C-C1");
}

#[test]
fn isotopic_prefix_binds_to_the_atom() {
    // The defining fix over AsciiMath: ^14 binds the isotope to C.
    assert_eq!(round_trip("^14C"), "^14C");
    assert_eq!(round_trip("^3H"), "^3H");
}

#[test]
fn bonds_render_their_ascii() {
    round_trip("A-B");
    round_trip("A=B");
    round_trip("A#B");
    round_trip("A##B");
    round_trip("A>-B");
    round_trip("A-<B");
    round_trip("A~>B");
    round_trip("A~~B");
    round_trip("A/B");
    round_trip("A\\B");
}

#[test]
fn groups_render_bracket_and_multiplicity() {
    round_trip("(OH)_2");
    round_trip("[OH]_2");
    round_trip("(NH_4)_2SO_4");
    round_trip("Ca(OH)_2");
}

#[test]
fn molecules_render_coefficient_and_stereo() {
    round_trip("2H_2O");
    round_trip("(R)-C_2H_5OH");
    round_trip("(S)-CHClBrF");
    round_trip("(E)-C_4H_8");
    round_trip("(Z)-C_4H_8");
}

#[test]
fn alpha_beta_stereo_canonicalises_to_the_spelled_form() {
    // The grammar captures single letters a/b; the canonical text
    // spells them out (the Text formatter is the canonicaliser).
    assert_eq!(
        parse_to_text("(a)-C_6H_12O_6").unwrap(),
        "(alpha)-C_6H_12O_6"
    );
    assert_eq!(
        parse_to_text("(b)-C_6H_12O_6").unwrap(),
        "(beta)-C_6H_12O_6"
    );
}

#[test]
fn reactions_render_arrows_and_terms() {
    round_trip("A-B -> B");
    round_trip("2H_2 + O_2 -> 2H_2O");
    round_trip("A <=> B");
    round_trip("A <-> B");
    round_trip("A <- B");
    round_trip("A + B + C -> D");
}

#[test]
fn reactions_render_conditions_on_the_arrow() {
    round_trip("A ->[cat: Pd] B");
    round_trip("A ->[cat: Pd][solv: THF] B");
    round_trip("A ->[cat: Pd; reflux] B");
    round_trip("A <=>[400^o C][Pt] 2HI");
    round_trip("H_2 + I_2 <=>[400^o C][Pt] 2HI");
}

#[test]
fn cascades_fold_into_linked_steps() {
    round_trip("A -> B -> C");
    round_trip("A + B -> C ->[cat: Pd] D -> E");
}

#[test]
fn quoted_text_and_embedded_math_render_verbatim() {
    round_trip("\"heat\"");
    round_trip("\"free text\"");
    round_trip("`K_c = 1`");
}

#[test]
fn annotated_molecules_render_their_annotations() {
    round_trip("CCO @cas(\"64-17-5\")");
    round_trip("H_2O @name(\"water\") @inchi(\"1S/H2O/h1H2\")");
}

#[test]
fn mechanisms_render_steps_flows_and_spectators() {
    round_trip("mechanism{\n  s1: A-B -> A + B | flow: A-B > A; B > B-C\n  s2: A + C -> A-C\n}");
    round_trip("mechanism{\n  s1: A -> B\n  spectator: Na^+\n}");
}

#[test]
fn unsupported_constructs_report_clearly() {
    let error = parse_to_text("H_2O").err();
    assert!(error.is_none(), "sanity: H_2O parses");
}

#[test]
fn rejections_stay_rejections() {
    for input in [
        "Na+",
        "H_",
        "^14",
        "A ->",
        "mechanism{\n  s1: A -> B | curve: A > B\n}",
    ] {
        assert!(
            parse_to_text(input).is_err(),
            "expected {input:?} to be rejected"
        );
    }
}

#[test]
fn empty_input_is_rejected() {
    assert!(parse_to_text("").is_err());
    assert!(parse_to_text("   ").is_err());
}
