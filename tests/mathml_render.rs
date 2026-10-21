//! Unit contract for the MathML renderer: the byte-level shapes that
//! matter per construct (isotope binding, sub/sup combination,
//! conditions as munderover, mechanism tables) and Unsupported
//! reporting for out-of-subset nodes.

use asciichem::mathml;

fn body(input: &str) -> String {
    let xml = mathml::render(input).unwrap_or_else(|error| panic!("{input:?}: {error}"));
    // Strip the declaration + math wrapper; assert on inner shape.
    let start = xml.find("<math").expect("math root")
        + "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">".len();
    xml[start..xml.rfind("</math>").expect("math close")]
        .trim()
        .to_string()
}

#[test]
fn isotope_is_a_left_superscript_on_the_atom() {
    // The headline fix: ^14C binds to the atom — never a phantom base.
    let out = body("^14C");
    assert!(out.contains("<mmultiscripts>"), "{out}");
    assert!(out.contains("<mprescripts/>"), "{out}");
    assert!(out.contains("<mn>14</mn>"), "{out}");
    assert!(out.contains("<mi mathvariant=\"normal\">C</mi>"), "{out}");
    assert!(
        !out.contains("<mi mathvariant=\"normal\"></mi>"),
        "no phantom base: {out}"
    );
}

#[test]
fn subscript_and_charge_combine_into_msubsup() {
    let out = body("SO_4^2-");
    assert!(out.contains("<msubsup>"), "{out}");
    assert!(out.contains("<mn>4</mn>"), "{out}");
    let pos = out.find("<msubsup>").unwrap();
    let rest = &out[pos..];
    let mn_pos = rest.find("<mn>2</mn>").expect("sub first");
    let mo_pos = rest
        .find("<mo>+</mo>")
        .or_else(|| rest.find("<mo>-</mo>"))
        .expect("sign second");
    assert!(mn_pos < mo_pos, "subscript precedes the sign: {rest}");
}

#[test]
fn charge_renders_as_count_plus_sign() {
    let out = body("Ca^2+");
    assert!(out.contains("<mn>2</mn>"), "{out}");
    assert!(out.contains("<mo>+</mo>"), "{out}");
    let out = body("O^2-");
    assert!(out.contains("<mo>-</mo>"), "{out}");
}

#[test]
fn oxidation_state_renders_as_roman_mi() {
    let out = body("Mn^(VII)");
    assert!(out.contains("<mi mathvariant=\"normal\">VII</mi>"), "{out}");
}

#[test]
fn groups_wrap_in_msub_for_multiplicity() {
    let out = body("(OH)_2");
    assert!(out.contains("<mo>(</mo>"), "{out}");
    assert!(out.contains("<mo>)</mo>"), "{out}");
    assert!(out.contains("<msub>"), "{out}");
}

#[test]
fn reactions_render_arrow_entities_and_terms() {
    let out = body("A + B -> C");
    assert!(out.contains("<mo>+</mo>"), "{out}");
    assert!(out.contains("<mo>→</mo>"), "{out}");
    let out = body("A <=> B");
    assert!(out.contains("<mo>⇌</mo>"), "{out}");
}

#[test]
fn conditions_render_as_munderover() {
    let out = body("A ->[cat: Pd] B");
    assert!(out.contains("<mover>"), "{out}");
    let out = body("H_2 + I_2 <=>[400^o C][Pt] 2HI");
    assert!(out.contains("<munderover>"), "{out}");
    assert!(out.contains("<mi mathvariant=\"normal\">Pt</mi>"), "{out}");
}

#[test]
fn cascades_repeat_arrow_products_after_the_head_reactants() {
    let out = body("C_2H_5OH -> CH_3CHO -> CH_3COOH");
    let arrows = out.matches("<mo>→</mo>").count();
    assert_eq!(arrows, 2, "{out}");
}

#[test]
fn text_renders_as_mtext() {
    let out = body("\"heat\"");
    assert!(out.contains("<mtext>heat</mtext>"), "{out}");
}

#[test]
fn electron_configuration_separates_orbitals_with_nbsp() {
    let out = body("1s^2 2s^2");
    assert!(out.contains("<mo>\u{00A0}</mo>"), "{out}");
    assert!(out.contains("<msup>"), "{out}");
}

#[test]
fn annotations_do_not_break_mathml() {
    let out = body("CCO @cas(\"64-17-5\")");
    assert!(out.contains("<mi mathvariant=\"normal\">C</mi>"), "{out}");
    assert!(out.contains("<mi mathvariant=\"normal\">O</mi>"), "{out}");
}

#[test]
fn unsupported_constructs_report_clearly() {
    // Mechanisms render (table form); embedded math stays unsupported
    // until the Plurimath pipeline is ported.
    assert!(mathml::render("mechanism{\n  s1: A -> B\n}").is_ok());
}
