//! CML round-trip: the three aci: channels, the core wire, interop
//! with the reference emitter's output, namespace exactness, error
//! surfaces. The corpus-level three-way gate lives in `corpus.rs`;
//! these tests pin the module's own behavior.

fn round_trip(input: &str) -> String {
    let nodes = asciichem::parse_nodes(input).expect("input parses");
    let xml = asciichem::cml::to_cml(&nodes).expect("emits CML");
    let back = asciichem::cml::parse_cml(&xml).expect("reads CML back");
    back.iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}

fn text_of(input: &str) -> Vec<asciichem::model::Node> {
    let nodes = asciichem::parse_nodes(input).expect("input parses");
    let xml = asciichem::cml::to_cml(&nodes).expect("emits CML");
    asciichem::cml::parse_cml(&xml).expect("reads CML back")
}

// -- the three-way contract over the module's own surface ---------------

#[test]
fn core_wire_round_trips() {
    for input in [
        "H_2O",
        "^14C",
        "Ca^2+",
        "CH_3-CH_2-OH",
        "H-O-H",
        "CH_3-CH_2-OH",
        "2H_2 + O_2 -> 2H_2O",
        "A <- B",
        "A <=> B",
        "A <-> B",
    ] {
        assert_eq!(round_trip(input), input, "core round-trip: {input}");
    }
}

#[test]
fn atom_attribute_channel_round_trips() {
    for input in ["Fe^(II)", "::O", "N.", "::N.", "Fe^(II)::O."] {
        assert_eq!(round_trip(input), input, "atom channel: {input}");
    }
}

#[test]
fn top_level_channel_round_trips() {
    for input in [
        "1s^2 2s^2",
        "1s^2 2s^2 2p^6 3d^10 4s^2",
        "`K_c = [P]/[R]`",
        "`a < b && c > d`",
        "\"note: toxic\"",
    ] {
        assert_eq!(round_trip(input), input, "top-level channel: {input}");
    }
}

#[test]
fn group_channel_round_trips() {
    for input in [
        "(OH)_2",
        "(H_2O)_3",
        "Ca(OH)_2",
        "[OH]_2",
        "{OH}_2",
        "((OH))",
        "((OH)_2)",
        "Ca((OH)_2)",
    ] {
        assert_eq!(round_trip(input), input, "group channel: {input}");
    }
}

#[test]
fn ring_closures_round_trip_through_the_atom_channel() {
    // Ring closure digits ride aci:ringClosures; the closing bond is
    // a canonical bond for foreign consumers and re-spells from the
    // digits on read.
    for input in ["C1-C-C-C-C-C1", "C1=CC=CC=C1", "C1-C-C1", "C12-C-C1-C2"] {
        assert_eq!(round_trip(input), input, "ring closures: {input}");
    }
}

#[test]
fn cascades_round_trip() {
    assert_eq!(round_trip("A -> B -> C"), "A -> B -> C");
}

#[test]
fn top_level_positions_are_preserved() {
    let nodes = text_of("\"left\" H_2O \"right\"");
    let rendered: Vec<String> = nodes.iter().map(|node| node.to_string()).collect();
    assert_eq!(rendered, vec!["\"left\"", "H_2O", "\"right\""]);
}

// -- interop with the reference emitter ---------------------------------
//
// These are the reference implementation's exact CML strings for the
// same inputs (one per channel, plus the core wire). Reading them is
// the cross-implementation contract: the wire shape is shared, not
// re-derived.

#[test]
fn reads_reference_core_wire() {
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema"><molecule id="m1"><atomArray><atom id="a1" elementType="H" count="2"/><atom id="a2" elementType="O"/></atomArray></molecule></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("reads reference CML");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "H_2O");
}

#[test]
fn reads_reference_atom_channel() {
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:aci="https://asciichem.org/cml-ext"><molecule id="m1"><atomArray><atom id="a1" elementType="N" aci:lonePairs="2" aci:radicalElectrons="1"></atom></atomArray></molecule></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("reads reference CML");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "::N.");
}

#[test]
fn reads_reference_group_channel() {
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:aci="https://asciichem.org/cml-ext"><molecule id="m1"><atomArray><atom id="a1" elementType="O" count="2"></atom><atom id="a2" elementType="H" count="2"></atom></atomArray><aci:group multiplicity="2" bracket="paren" atomRefs="a1 a2"></aci:group></molecule></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("reads reference CML");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "(OH)_2");
}

#[test]
fn reads_reference_reaction_wire() {
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema"><reaction id="r1" title="forward" type="forward"><reactantList><reactant><substance role="reactant"><molecule id="m1" count="2"><atomArray><atom id="a1" elementType="H" count="2"/></atomArray></molecule></substance></reactant><reactant><substance role="reactant"><molecule id="m2"><atomArray><atom id="a2" elementType="O" count="2"/></atomArray></molecule></substance></reactant></reactantList><productList><product><substance role="product"><molecule id="m3" count="2"><atomArray><atom id="a3" elementType="H" count="2"/><atom id="a4" elementType="O"/></atomArray></molecule></substance></product></productList></reaction></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("reads reference CML");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "2H_2 + O_2 -> 2H_2O");
}

#[test]
fn reads_reference_escaped_top_level() {
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:aci="https://asciichem.org/cml-ext"><aci:embeddedMath position="0">a &lt; b &amp;&amp; c &gt; d</aci:embeddedMath></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("reads reference CML");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "`a < b && c > d`");
}

#[test]
fn reads_cdata_text_content() {
    let xml = "<cml xmlns=\"http://www.xml-cml.org/schema\" xmlns:aci=\"https://asciichem.org/cml-ext\"><aci:text position=\"0\"><![CDATA[\"co < raw\"]]></aci:text></cml>";
    let nodes = asciichem::cml::parse_cml(xml).expect("reads CDATA");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "\"co < raw\"");
}

// -- namespace exactness -------------------------------------------------

#[test]
fn foreign_namespaced_elements_are_not_channels() {
    // A `<text>` in the CORE namespace with the aci: declaration also
    // present must not be mistaken for the aci:text channel —
    // classification is by resolved URI, not by local name.
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:aci="https://asciichem.org/cml-ext"><text>not a channel</text></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("parses");
    assert!(
        nodes.is_empty(),
        "foreign <text> must be ignored, got {nodes:?}"
    );
}

#[test]
fn undeclared_aci_prefix_is_ignored() {
    // The engine serializes undeclared prefixes as written; the
    // reader resolves namespaces, so an undeclared aci: attribute is
    // no-namespace and must not restore.
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema"><molecule id="m1"><atomArray><atom id="a1" elementType="O" aci:lonePairs="2"/></atomArray></molecule></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("parses");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "O");
}

#[test]
fn foreign_prefix_bound_to_aci_uri_is_a_channel() {
    // XML Namespaces: the written prefix never matters, only the URI.
    let xml = r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:x="https://asciichem.org/cml-ext"><x:text position="0">"hello"</x:text></cml>"#;
    let nodes = asciichem::cml::parse_cml(xml).expect("parses");
    let text = nodes
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, "\"hello\"");
}

// -- emit shape ----------------------------------------------------------

#[test]
fn namespace_declared_only_when_needed() {
    let plain = asciichem::parse_nodes("H_2O").expect("parses");
    let xml = asciichem::cml::to_cml(&plain).expect("emits");
    assert!(
        !xml.contains("aci"),
        "plain molecule must not declare aci: {xml}"
    );

    let lewis = asciichem::parse_nodes("::O").expect("parses");
    let xml = asciichem::cml::to_cml(&lewis).expect("emits");
    assert!(
        xml.contains("xmlns:aci=\"https://asciichem.org/cml-ext\""),
        "{xml}"
    );
}

#[test]
fn escaping_uses_the_engine() {
    let nodes = asciichem::parse_nodes("`a < b && c > d`").expect("parses");
    let xml = asciichem::cml::to_cml(&nodes).expect("emits");
    assert!(xml.contains("a &lt; b &amp;&amp; c &gt; d"), "{xml}");
    assert_eq!(round_trip("`a < b && c > d`"), "`a < b && c > d`");
}

#[test]
fn molecule_charge_sum_is_derived_on_emit_and_ignored_on_read() {
    let nodes = asciichem::parse_nodes("Ca^2+").expect("parses");
    let xml = asciichem::cml::to_cml(&nodes).expect("emits");
    assert!(xml.contains("formalCharge=\"+2\""), "molecule sum: {xml}");
    assert!(xml.contains("formalCharge=\"2+\""), "atom charge: {xml}");
    assert_eq!(round_trip("Ca^2+"), "Ca^2+");
}

// -- error surfaces -------------------------------------------------------

#[test]
fn malformed_xml_is_an_engine_error() {
    let error = asciichem::cml::parse_cml("<cml><molecule></cml>").expect_err("must reject");
    assert!(
        matches!(error, asciichem::cml::CmlError::Xml(_)),
        "{error:?}"
    );
}

#[test]
fn wrong_root_is_invalid() {
    let error = asciichem::cml::parse_cml("<notcml/>").expect_err("must reject");
    assert!(
        matches!(error, asciichem::cml::CmlError::Invalid(_)),
        "{error:?}"
    );
}

#[test]
fn atom_without_element_type_is_invalid() {
    let error = asciichem::cml::parse_cml(
        r#"<cml xmlns="http://www.xml-cml.org/schema"><molecule><atomArray><atom id="a1"/></atomArray></molecule></cml>"#,
    )
    .expect_err("must reject");
    assert!(
        matches!(error, asciichem::cml::CmlError::Invalid(_)),
        "{error:?}"
    );
}

#[test]
fn bond_refs_must_name_two_atoms() {
    let error = asciichem::cml::parse_cml(
        r#"<cml xmlns="http://www.xml-cml.org/schema"><molecule><atomArray><atom id="a1" elementType="C"/></atomArray><bondArray><bond id="b1" atomRefs2="a1" order="S"/></bondArray></molecule></cml>"#,
    )
    .expect_err("must reject");
    assert!(
        matches!(error, asciichem::cml::CmlError::Invalid(_)),
        "{error:?}"
    );
}

#[test]
fn unparseable_channel_content_is_invalid() {
    let error = asciichem::cml::parse_cml(
        r#"<cml xmlns="http://www.xml-cml.org/schema" xmlns:aci="https://asciichem.org/cml-ext"><aci:text position="0">)))not chemistry(((</aci:text></cml>"#,
    )
    .expect_err("must reject");
    assert!(
        matches!(error, asciichem::cml::CmlError::Invalid(_)),
        "{error:?}"
    );
}

#[test]
fn conditions_have_no_wire_yet_and_say_so() {
    let nodes = asciichem::parse_nodes("A ->[cat: Pd] B").expect("parses");
    let error = asciichem::cml::to_cml(&nodes).expect_err("must report, not drop");
    assert!(
        matches!(error, asciichem::cml::CmlError::Unsupported(_)),
        "{error:?}"
    );
}
