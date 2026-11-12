//! Inline HTML renderer — the Rust counterpart of the reference
//! `Formatter::Html`: `<sub>`/`<sup>` scripts inside one
//! `<span class="asciichem">`, no MathJax or MathML dependency.
//! Byte-shaped after the reference visitor (the same four-character
//! XML escaping, the same empty joins); embedded math stays
//! `Unsupported` (the reference routes it through the Plurimath
//! pipeline the Rust core does not embed).

use crate::model::{Crystal, Node};

/// Renders already-built nodes as inline HTML — the
/// grammar-independent seam (`convert -t html`, any `--from` grammar).
pub fn render_formula(nodes: &[Node]) -> Result<String, crate::ParseError> {
    let inner = render_nodes(nodes)?;
    Ok(format!("<span class=\"asciichem\">{inner}</span>"))
}

fn render_nodes(nodes: &[Node]) -> Result<String, crate::ParseError> {
    let mut out = String::new();
    for node in nodes {
        out.push_str(&render_node(node)?);
    }
    Ok(out)
}

fn render_node(node: &Node) -> Result<String, crate::ParseError> {
    match node {
        Node::Atom(atom) => Ok(render_atom(atom)),
        Node::Bond(bond) => Ok(escape(bond.ascii())),
        Node::Group(group) => {
            let (open, close) = group.chars();
            let body = render_nodes(&group.nodes)?;
            let suffix = match &group.multiplicity {
                Some(mult) => format!("<sub>{}</sub>", escape(mult)),
                None => String::new(),
            };
            Ok(format!("{}{body}{}{suffix}", escape(open), escape(close)))
        }
        Node::Molecule(mol) => {
            let prefix = match &mol.coefficient {
                Some(coef) if !coef.is_empty() => escape(coef),
                _ => String::new(),
            };
            let stereo = match mol.stereo_letter() {
                Some(letter) => format!("({})-", escape(letter)),
                None => String::new(),
            };
            let body = render_nodes(&mol.nodes)?;
            Ok(format!("{stereo}{prefix}{body}"))
        }
        Node::Reaction(reaction) => {
            let left = render_terms(&reaction.reactants)?;
            let right = render_terms(&reaction.products)?;
            let arrow = escape(reaction.arrow_entity());
            let Some(conds) = &reaction.conditions else {
                return Ok(format!("{left} {arrow} {right}"));
            };
            let above = match &conds.above {
                Some(text) => format!("<sup>{}</sup>", escape(text)),
                None => String::new(),
            };
            let below = match &conds.below {
                Some(text) => format!("<sub>{}</sub>", escape(text)),
                None => String::new(),
            };
            Ok(format!("{left} {above}{arrow}{below} {right}"))
        }
        Node::ReactionCascade(cascade) => {
            let Some(head) = cascade.steps.first() else {
                return Ok(String::new());
            };
            let mut out = render_terms(&head.reactants)?;
            for step in &cascade.steps {
                let arrow = escape(step.arrow_entity());
                if let Some(conds) = &step.conditions {
                    let above = match &conds.above {
                        Some(text) => format!("<sup>{}</sup>", escape(text)),
                        None => String::new(),
                    };
                    let below = match &conds.below {
                        Some(text) => format!("<sub>{}</sub>", escape(text)),
                        None => String::new(),
                    };
                    out.push_str(&format!(" {above}{arrow}{below}"));
                } else {
                    out.push_str(&format!(" {arrow}"));
                }
                out.push(' ');
                out.push_str(&render_terms(&step.products)?);
            }
            Ok(out)
        }
        Node::ElectronConfiguration(ec) => {
            let orbitals: Vec<String> = ec
                .orbitals
                .iter()
                .map(|(orb, occ)| format!("{}<sup>{}</sup>", escape(orb), escape(occ)))
                .collect();
            Ok(orbitals.join(" "))
        }
        Node::Mechanism(mech) => {
            let mut parts = vec![r#"<span class="asciichem-mechanism">"#.to_string()];
            if mech.steps.is_empty() && mech.spectators.is_empty() {
                parts.push("</span>".to_string());
                return Ok(parts.concat());
            }
            parts.push(r#"<dl class="asciichem-steps">"#.to_string());
            for step in &mech.steps {
                parts.push(format!(
                    "<dt>{}</dt><dd>{}</dd>",
                    escape(&step.label),
                    escape(&step.reaction)
                ));
            }
            for spectator in &mech.spectators {
                parts.push(format!("<dt>spectator</dt><dd>{}</dd>", escape(spectator)));
            }
            parts.push("</dl></span>".to_string());
            Ok(parts.concat())
        }
        Node::Crystal(crystal) => {
            let mut parts = vec![r#"<span class="asciichem-crystal">"#.to_string()];
            if let Some(name) = &crystal.name {
                parts.push(format!(
                    "<span class=\"asciichem-name\">[{}]</span>",
                    escape(name)
                ));
            }
            let pairs: Vec<(String, &str)> = crystal
                .cell_params()
                .into_iter()
                .map(|(label, value)| (Crystal::greek_label(label, "html"), value))
                .collect();
            parts.push(kv_dl(&pairs));
            if !crystal.atoms.is_empty() {
                parts.push(r#"<ol class="asciichem-atoms">"#.to_string());
                for atom in &crystal.atoms {
                    parts.push(format!("<li>{}</li>", render_atom(atom)));
                }
                parts.push("</ol>".to_string());
            }
            parts.push("</span>".to_string());
            Ok(parts.concat())
        }
        Node::EmbeddedMath(_) => Err(crate::ParseError {
            message: "not yet supported in the Rust core: embedded math rendering \
                      (Plurimath pipeline)"
                .into(),
        }),
        Node::Text(content) => Ok(escape(content)),
    }
}

fn render_terms(nodes: &[Node]) -> Result<String, crate::ParseError> {
    let terms: Result<Vec<String>, crate::ParseError> = nodes.iter().map(render_node).collect();
    Ok(terms?.join(" + "))
}

/// One atom: isotope, element, subscript, superscript, charge,
/// oxidation state, ring closures — `<sub>`/`<sup>` scripts (the
/// reference `visit_atom`; lone pairs and radicals stay text-only).
fn render_atom(atom: &crate::model::Atom) -> String {
    let mut parts = Vec::new();
    if let Some(isotope) = &atom.isotope {
        parts.push(format!("<sup>{}</sup>", escape(isotope)));
    }
    parts.push(escape(&atom.element));
    if let Some(sub) = &atom.subscript {
        parts.push(format!("<sub>{}</sub>", escape(sub)));
    }
    if let Some(sup) = &atom.superscript {
        parts.push(format!("<sup>{}</sup>", escape(sup)));
    }
    if let Some(charge) = &atom.charge {
        parts.push(format!("<sup>{}</sup>", escape(charge)));
    }
    if let Some(ox) = &atom.oxidation_state {
        parts.push(format!("<sup>({})</sup>", escape(ox)));
    }
    if let Some(rc) = &atom.ring_closures {
        parts.push(escape(rc));
    }
    parts.concat()
}

/// The reference's minimal escaper: the model content is plain text
/// (digits, element symbols, conditions text), so only the four
/// XML-mandated characters need entities. Absent pieces render "".
fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Generic key-value pairs as `<dl>`: `<dt>key</dt><dd>value</dd>`
/// rows; empty input renders "" (the reference's `kv_dl`).
fn kv_dl(pairs: &[(String, &str)]) -> String {
    if pairs.is_empty() {
        return String::new();
    }
    let mut parts = vec!["<dl>".to_string()];
    for (key, value) in pairs {
        parts.push(format!(
            "<dt>{}</dt><dd>{}</dd>",
            escape(key),
            escape(value)
        ));
    }
    parts.push("</dl>".to_string());
    parts.concat()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn html(input: &str) -> String {
        render_formula(&crate::parse_nodes(input).expect("input parses")).expect("renders")
    }

    #[test]
    fn water_renders_the_subscript_span() {
        assert_eq!(
            html("H_2O"),
            "<span class=\"asciichem\">H<sub>2</sub>O</span>"
        );
    }

    #[test]
    fn isotope_charge_and_oxidation_stack_as_superscripts() {
        assert_eq!(
            html("^14C^2-"),
            "<span class=\"asciichem\"><sup>14</sup>C<sup>2-</sup></span>"
        );
        assert_eq!(
            html("Fe^(III)"),
            "<span class=\"asciichem\">Fe<sup>(III)</sup></span>"
        );
    }

    #[test]
    fn reactions_render_the_entity_arrow_and_condition_scripts() {
        assert_eq!(
            html("H_2 + O_2 -> H_2O"),
            "<span class=\"asciichem\">H<sub>2</sub> + O<sub>2</sub> → H<sub>2</sub>O</span>"
        );
    }

    #[test]
    fn text_is_escaped() {
        assert_eq!(
            html("\"a<b & c\""),
            "<span class=\"asciichem\">a&lt;b &amp; c</span>"
        );
    }

    #[test]
    fn crystal_renders_the_name_cell_table_and_atom_list() {
        assert_eq!(
            html("crystal[NaCl](a=5.64,alpha=90){Na@f(0,0,0)}"),
            "<span class=\"asciichem\">\
             <span class=\"asciichem-crystal\">\
             <span class=\"asciichem-name\">[NaCl]</span>\
             <dl><dt>a</dt><dd>5.64</dd><dt>α</dt><dd>90</dd></dl>\
             <ol class=\"asciichem-atoms\"><li>Na</li></ol>\
             </span></span>"
        );
    }

    #[test]
    fn embedded_math_stays_unsupported() {
        let error =
            render_formula(&crate::parse_nodes("`x^2`").expect("parses")).expect_err("unsupported");
        assert!(error.message.contains("Plurimath"));
    }
}
