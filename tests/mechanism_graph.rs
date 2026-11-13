//! The mechanism step graph: the Rust port of the reference's
//! `MechanismGraph` spec — each step's reaction re-parsed with the same
//! grammar, catalysts (species on both sides of one step), links (an
//! earlier step's product consumed as a later step's reactant),
//! cascade entry/exit sides, coefficient-free matching, and prose
//! entries kept as non-reactions. Also pins the mechanism's renderer
//! bytes (MathML/HTML/LaTeX/SVG element-for-element with the Ruby
//! reference) and the wire round trip through a conditioned step.

fn graph_for(source: &str) -> asciichem::mechanism_graph::Graph {
    let nodes = asciichem::parse_nodes(source).expect("mechanism parses");
    match nodes.first() {
        Some(asciichem::model::Node::Mechanism(mechanism)) => {
            asciichem::mechanism_graph::for_mechanism(mechanism)
        }
        other => panic!("expected a mechanism node, got {other:?}"),
    }
}

fn mechanism_body(lines: &[&str]) -> String {
    format!("mechanism{{\n  {}\n}}", lines.join("\n  "))
}

fn link_shapes(graph: &asciichem::mechanism_graph::Graph) -> Vec<(String, String, Vec<String>)> {
    graph
        .links
        .iter()
        .map(|link| (link.from.clone(), link.to.clone(), link.species.clone()))
        .collect()
}

#[test]
fn parses_each_step_value_into_a_reaction() {
    let graph = graph_for(&mechanism_body(&["step1: A + B -> C"]));
    assert_eq!(graph.steps.len(), 1);
    assert_eq!(graph.steps[0].label, "step1");
    assert!(matches!(
        graph.steps[0].reaction,
        Some(asciichem::model::Node::Reaction(_))
    ));
    assert!(graph.steps[0].catalysts.is_empty());
}

#[test]
fn links_a_later_step_that_consumes_an_earlier_product() {
    let graph = graph_for(&mechanism_body(&["s1: A + B -> C", "s2: C -> D"]));
    assert_eq!(
        link_shapes(&graph),
        vec![("s1".into(), "s2".into(), vec!["C".into()])]
    );
}

#[test]
fn links_chains_transitively_step_to_step_not_across_gaps() {
    let graph = graph_for(&mechanism_body(&["s1: A -> B", "s2: B -> C", "s3: C -> D"]));
    assert_eq!(
        link_shapes(&graph),
        vec![
            ("s1".into(), "s2".into(), vec!["B".into()]),
            ("s2".into(), "s3".into(), vec!["C".into()]),
        ]
    );
}

#[test]
fn does_not_link_steps_with_no_shared_species() {
    let graph = graph_for(&mechanism_body(&["s1: A -> B", "s2: X -> Y"]));
    assert!(graph.links.is_empty());
}

#[test]
fn recognises_a_species_on_both_sides_of_one_step_as_a_catalyst() {
    let graph = graph_for(&mechanism_body(&["s1: A + Pd -> B + Pd"]));
    assert_eq!(graph.steps[0].catalysts, vec!["Pd"]);
}

#[test]
fn links_through_a_shared_catalyst_species_as_well() {
    let graph = graph_for(&mechanism_body(&[
        "s1: A + Pd -> B + Pd",
        "s2: Pd + C -> D",
    ]));
    assert_eq!(
        graph
            .links
            .iter()
            .map(|l| l.species.clone())
            .collect::<Vec<_>>(),
        vec![vec!["Pd".to_string()]]
    );
}

#[test]
fn keeps_prose_entries_as_non_reactions_with_no_links() {
    let graph = graph_for(&mechanism_body(&["s1: proton transfer", "s2: A -> B"]));
    assert_eq!(graph.steps[0].reaction, None);
    assert!(matches!(
        graph.steps[1].reaction,
        Some(asciichem::model::Node::Reaction(_))
    ));
    assert!(graph.links.is_empty());
}

#[test]
fn uses_a_cascade_entry_side_and_exit_side_for_linking() {
    let graph = graph_for(&mechanism_body(&["s1: A -> B -> C", "s2: C -> D"]));
    assert_eq!(
        link_shapes(&graph),
        vec![("s1".into(), "s2".into(), vec!["C".into()])]
    );
}

#[test]
fn ignores_coefficients_when_matching_species() {
    let graph = graph_for(&mechanism_body(&["s1: 2A -> 3B", "s2: 2B -> C"]));
    assert_eq!(
        link_shapes(&graph),
        vec![("s1".into(), "s2".into(), vec!["B".into()])]
    );
}

#[test]
fn carries_flows_into_the_step_graph() {
    let graph = graph_for(&mechanism_body(&["s1: A -> B | flow: A > B"]));
    let flows = &graph.steps[0].flows;
    assert_eq!(flows.len(), 1);
    assert_eq!(flows[0].source, "A");
    assert_eq!(flows[0].target, "B");
}

#[test]
fn leaves_the_canonical_text_and_model_unchanged() {
    let source = mechanism_body(&["s1: A -> B", "s2: B -> C"]);
    assert_eq!(asciichem::parse_to_text(&source).unwrap(), source);
}

// -- renderer bytes (element-for-element with the Ruby reference) -------

#[test]
fn mathml_renders_the_step_table_byte_exact() {
    let out =
        asciichem::mathml::render("mechanism{\n  s1: A -> B | flow: A > B\n  spectator: Na^+\n}")
            .unwrap();
    let expected = r#"<?xml version="1.0" encoding="UTF-8"?>
<math xmlns="http://www.w3.org/1998/Math/MathML">
  <mrow>
    <mrow>
      <mi mathvariant="normal">mechanism</mi>
      <mtable>
        <mtr>
          <mtd>
            <mi mathvariant="normal">s1</mi>
          </mtd>
          <mtd>
            <mtext>A -&gt; B</mtext>
          </mtd>
        </mtr>
        <mtr>
          <mtd>
            <mi mathvariant="normal">spectator</mi>
          </mtd>
          <mtd>
            <mtext>Na^+</mtext>
          </mtd>
        </mtr>
      </mtable>
    </mrow>
  </mrow>
</math>
"#;
    assert_eq!(out, expected);
}

#[test]
fn mathml_renders_a_bare_mechanism_as_the_bare_marker() {
    // The body-less forms build the empty mechanism (the model carries
    // no steps and no spectators), rendering as the bare marker.
    for input in ["mechanism", "mechanism{}"] {
        let out = asciichem::mathml::render(input).unwrap();
        assert!(
            out.contains("<mi mathvariant=\"normal\">mechanism</mi>"),
            "{out}"
        );
        assert!(!out.contains("mtable"), "{out}");
        assert_eq!(asciichem::parse_to_text(input).unwrap(), "mechanism");
    }
}

#[test]
fn html_renders_the_steps_definition_list() {
    let nodes =
        asciichem::parse_nodes("mechanism{\n  s1: A -> B | flow: A > B\n  spectator: Na^+\n}")
            .unwrap();
    let out = asciichem::html::render_formula(&nodes).unwrap();
    let inner = "<span class=\"asciichem-mechanism\"><dl class=\"asciichem-steps\">\
                 <dt>s1</dt><dd>A -&gt; B</dd>\
                 <dt>spectator</dt><dd>Na^+</dd></dl></span>";
    assert!(out.contains(inner), "{out}");
}

#[test]
fn latex_renders_the_step_tabular() {
    let nodes =
        asciichem::parse_nodes("mechanism{\n  s1: A -> B | flow: A > B\n  spectator: Na^+\n}")
            .unwrap();
    let out = asciichem::latex::render_formula(&nodes).unwrap();
    let inner = "\\text{mechanism}\\begin{tabular}{ll}\n\
                s1 & A -> B \\\\\n\
                spectator & Na^+ \\\\\n\
                \\end{tabular}";
    assert!(out.contains(inner), "{out}");
}

#[test]
fn svg_renders_the_full_canonical_text() {
    let input = "mechanism{\n  s1: A -> B | flow: A > B\n  spectator: Na^+\n}";
    for out in [
        asciichem::svg::render(input).unwrap(),
        asciichem::svg::render_linear(input).unwrap(),
    ] {
        // The reference renders constructs without a natural 2D form
        // through the canonical text — the full entry set, spectators
        // included (never just the first line).
        assert!(out.contains("mechanism{"), "{out}");
        assert!(out.contains("spectator: Na^+"), "{out}");
    }
}

// -- wire -----------------------------------------------------------------

#[test]
fn wire_round_trips_flows_spectators_and_conditioned_steps() {
    let source = "mechanism{\n  s1: A-B -> A + B | flow: A-B > A; B > B-C\n  \
                  s2: 2C ->[cat: Pd] D\n  spectator: Na^+\n}";
    let model = asciichem::parse_nodes(source).unwrap();
    let json = asciichem::wire::to_model_json(&model).unwrap();
    let round = asciichem::wire::from_model_json(&json).unwrap();
    let text = round
        .iter()
        .map(|node| node.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    assert_eq!(text, source);

    let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
    let steps = &parsed["nodes"][0]["steps"];
    assert_eq!(steps[0]["flows"][0]["source"], "A-B");
    assert_eq!(steps[0]["flows"][1]["target"], "B-C");
    assert_eq!(steps[1]["reaction"]["conditions"]["above"], "cat: Pd");
    assert_eq!(
        parsed["nodes"][0]["spectators"][0]["nodes"][0]["charge"],
        "+"
    );
}
