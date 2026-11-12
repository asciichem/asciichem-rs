//! mhchem LaTeX renderer — the Rust counterpart of the reference
//! `Formatter::Latex`: the formula wraps in `\ce{...}` (mhchem, version
//! 4), atoms carry bare-digit scripts with braces only where mhchem
//! needs them, and the beyond-formula constructs render as
//! `\text{kind}[label]{ tabular body }`. Embedded math stays
//! `Unsupported` (the reference routes it through the Plurimath
//! pipeline the Rust core does not embed).

use crate::model::{Crystal, Node};

/// Renders already-built nodes as LaTeX — the grammar-independent
/// seam (`convert -t latex`, any `--from` grammar).
pub fn render_formula(nodes: &[Node]) -> Result<String, crate::ParseError> {
    let mut renderer = Latex { inside_ce: false };
    renderer.formula(nodes)
}

/// The renderer state: whether the walk is currently inside the
/// formula-level `\ce{...}` (molecules render bare inside it, wrapped
/// on their own otherwise — the reference's `@inside_ce`).
struct Latex {
    inside_ce: bool,
}

impl Latex {
    fn formula(&mut self, nodes: &[Node]) -> Result<String, crate::ParseError> {
        let previous = self.inside_ce;
        self.inside_ce = true;
        let inner = self.render_nodes(nodes);
        self.inside_ce = previous;
        Ok(format!("\\ce{{{}}}", inner?))
    }

    fn render_nodes(&mut self, nodes: &[Node]) -> Result<String, crate::ParseError> {
        let mut out = String::new();
        for node in nodes {
            out.push_str(&self.render_node(node)?);
        }
        Ok(out)
    }

    fn render_node(&mut self, node: &Node) -> Result<String, crate::ParseError> {
        match node {
            Node::Atom(atom) => Ok(render_atom(atom)),
            Node::Bond(bond) => Ok(match bond.kind.as_str() {
                "single" => "-".to_string(),
                "double" => "=".to_string(),
                "triple" => "\\equiv{}".to_string(),
                _ => bond.ascii().to_string(),
            }),
            Node::Group(group) => {
                let (open, close) = group.chars();
                let body = self.render_nodes(&group.nodes)?;
                let suffix = match &group.multiplicity {
                    Some(mult) => mhchem_subscript(mult),
                    None => String::new(),
                };
                Ok(format!("{open}{body}{close}{suffix}"))
            }
            Node::Molecule(mol) => {
                let prefix = match &mol.coefficient {
                    Some(coef) if !coef.is_empty() => coef.clone(),
                    _ => String::new(),
                };
                let stereo = match mol.stereo_letter() {
                    Some(letter) => format!("({letter})-"),
                    None => String::new(),
                };
                let body = self.render_nodes(&mol.nodes)?;
                let rendered = format!("{stereo}{prefix}{body}");
                if self.inside_ce {
                    Ok(rendered)
                } else {
                    Ok(format!("\\ce{{{rendered}}}"))
                }
            }
            Node::Reaction(reaction) => {
                let left = self.render_terms(&reaction.reactants)?;
                let right = self.render_terms(&reaction.products)?;
                let arrow = arrow_for(&reaction.arrow);
                let Some(conds) = &reaction.conditions else {
                    return Ok(format!("{left} {arrow} {right}"));
                };
                Ok(format!(
                    "{left} {arrow}{}{} {right}",
                    conditions_bracket(conds.above.as_deref()),
                    conditions_bracket(conds.below.as_deref())
                ))
            }
            Node::ReactionCascade(cascade) => {
                let Some(head) = cascade.steps.first() else {
                    return Ok(String::new());
                };
                let mut out = self.render_terms(&head.reactants)?;
                for step in &cascade.steps {
                    let arrow = arrow_for(&step.arrow);
                    if let Some(conds) = &step.conditions {
                        out.push_str(&format!(
                            " {arrow}{}{}",
                            conditions_bracket(conds.above.as_deref()),
                            conditions_bracket(conds.below.as_deref())
                        ));
                    } else {
                        out.push_str(&format!(" {arrow}"));
                    }
                    out.push(' ');
                    out.push_str(&self.render_terms(&step.products)?);
                }
                Ok(out)
            }
            Node::ElectronConfiguration(ec) => {
                let orbitals: Vec<String> = ec
                    .orbitals
                    .iter()
                    .map(|(orb, occ)| format!("{}^{}", wrap(orb), wrap(occ)))
                    .collect();
                Ok(format!("\\ce{{{}}}", orbitals.join(" ")))
            }
            Node::Mechanism(mech) => {
                let mut parts = vec!["\\text{mechanism}".to_string()];
                if mech.steps.is_empty() && mech.spectators.is_empty() {
                    return Ok(parts.concat());
                }
                let mut rows: Vec<Vec<String>> = mech
                    .steps
                    .iter()
                    .map(|step| vec![step.label.clone(), step.reaction.clone()])
                    .collect();
                for spectator in &mech.spectators {
                    rows.push(vec!["spectator".to_string(), spectator.clone()]);
                }
                parts.push(latex_tabular(&rows, None));
                Ok(parts.concat())
            }
            Node::Crystal(crystal) => {
                let mut parts = vec!["\\text{crystal}".to_string()];
                if let Some(name) = &crystal.name {
                    parts.push(format!("[{}]", wrap_latex(name)));
                }
                let pairs: Vec<Vec<String>> = crystal
                    .cell_params()
                    .into_iter()
                    .map(|(label, value)| {
                        vec![Crystal::greek_label(label, "latex"), value.to_string()]
                    })
                    .collect();
                if !pairs.is_empty() {
                    parts.push(latex_tabular(&pairs, None));
                }
                if !crystal.atoms.is_empty() {
                    let rows: Vec<Vec<String>> = crystal
                        .atoms
                        .iter()
                        .map(|atom| vec![render_atom(atom)])
                        .collect();
                    parts.push(latex_tabular(&rows, Some(1)));
                }
                Ok(parts.concat())
            }
            Node::EmbeddedMath(_) => Err(crate::ParseError {
                message: "not yet supported in the Rust core: embedded math rendering \
                     (Plurimath pipeline)"
                    .into(),
            }),
            Node::Text(content) => Ok(content.clone()),
        }
    }

    fn render_terms(&mut self, nodes: &[Node]) -> Result<String, crate::ParseError> {
        let mut terms = Vec::new();
        for node in nodes {
            terms.push(self.render_node(node)?);
        }
        Ok(terms.join(" + "))
    }
}

/// One atom: lone-pair colons, isotope, element, subscript, the
/// superscript-style marker (charge > oxidation state > raw
/// superscript), radical dots, ring closures — mhchem script shapes
/// (the reference `visit_atom`; fractional coordinates stay
/// model-only).
fn render_atom(atom: &crate::model::Atom) -> String {
    let mut parts = Vec::new();
    if let Some(lp) = atom.lone_pairs {
        parts.push(":".repeat(lp as usize));
    }
    if let Some(isotope) = &atom.isotope {
        parts.push(format!("^{}", wrap(isotope)));
    }
    parts.push(wrap(&atom.element));
    if let Some(sub) = &atom.subscript {
        parts.push(mhchem_subscript(sub));
    }
    if let Some(charge) = &atom.charge {
        parts.push(format!("^{}", wrap(charge)));
    } else if let Some(ox) = &atom.oxidation_state {
        parts.push(format!("^{{({})}}", wrap(ox)));
    } else if let Some(sup) = &atom.superscript {
        parts.push(format!("^{}", wrap(sup)));
    }
    if let Some(re) = atom.radical_electrons {
        parts.push(".".repeat(re as usize));
    }
    if let Some(rc) = &atom.ring_closures {
        parts.push(rc.clone());
    }
    parts.concat()
}

/// mhchem accepts bare digits after an element as subscripts: bare
/// form for digit-only values, brace form for anything else.
fn mhchem_subscript(value: &str) -> String {
    if value.bytes().all(|byte| byte.is_ascii_digit()) && !value.is_empty() {
        return value.to_string();
    }
    format!("_{value}")
}

/// mhchem accepts bare alphanumeric values for scripts (`^14`, `_12`);
/// braces only for non-digit-letter values to keep output idiomatic.
fn wrap(value: &str) -> String {
    if !value.is_empty() && value.bytes().all(|byte| byte.is_ascii_alphanumeric()) {
        return value.to_string();
    }
    format!("{{{value}}}")
}

fn arrow_for(arrow: &str) -> &'static str {
    match arrow {
        "reverse" => "<-",
        "equilibrium" => "<=>",
        "resonance" => "<->",
        _ => "->",
    }
}

/// Conditions render as mhchem brackets above/below: `[value]`, ""
/// when absent.
fn conditions_bracket(value: Option<&str>) -> String {
    match value {
        Some(text) if !text.is_empty() => format!("[{text}]"),
        _ => String::new(),
    }
}

/// Wrap content in `\text{}` if it carries spaces or special chars, so
/// tabular cells render cleanly in both pdflatex and xelatex.
fn wrap_latex(content: &str) -> String {
    let allowed = content.bytes().all(|byte| {
        byte.is_ascii_alphanumeric() || matches!(byte, b'+' | b'-' | b'/' | b'.' | b'*')
    }) && !content.is_empty();
    if allowed {
        return content.to_string();
    }
    format!("\\text{{{content}}}")
}

/// A `tabular` environment from rows of cells. Column count is the
/// explicit `columns` or the max row length (rows may be ragged);
/// every row line ends with ` \\` and rows join on newlines.
fn latex_tabular(rows: &[Vec<String>], columns: Option<usize>) -> String {
    if rows.is_empty() {
        return String::new();
    }
    let cols = columns.unwrap_or_else(|| rows.iter().map(Vec::len).max().unwrap_or(0));
    let alignment = "l".repeat(cols);
    let lines: Vec<String> = rows
        .iter()
        .map(|cells| {
            let mut padded: Vec<String> = cells.iter().take(cols).cloned().collect();
            padded.resize(cols, String::new());
            format!("{} \\\\", padded.join(" & "))
        })
        .collect();
    format!(
        "\\begin{{tabular}}{{{}}}\n{}\n\\end{{tabular}}",
        alignment,
        lines.join("\n")
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn latex(input: &str) -> String {
        render_formula(&crate::parse_nodes(input).expect("input parses")).expect("renders")
    }

    #[test]
    fn water_renders_the_mhchem_form() {
        assert_eq!(latex("H_2O"), "\\ce{H2O}");
    }

    #[test]
    fn non_digit_scripts_take_braces() {
        assert_eq!(latex("Fe^2-"), "\\ce{Fe^{2-}}");
        assert_eq!(latex("Fe^(III)"), "\\ce{Fe^{(III)}}");
        // mhchem accepts bare multi-digit scripts (the reference spec).
        assert_eq!(latex("^14C"), "\\ce{^14C}");
    }

    #[test]
    fn bonds_map_to_mhchem_glyphs() {
        assert_eq!(latex("C=C"), "\\ce{C=C}");
        assert_eq!(latex("C#C"), "\\ce{C\\equiv{}C}");
    }

    #[test]
    fn reactions_render_arrow_and_condition_brackets() {
        assert_eq!(
            latex("H_2 + O_2 ->[spark] H_2O"),
            "\\ce{H2 + O2 ->[spark] H2O}"
        );
    }

    #[test]
    fn text_passes_through_unescaped() {
        assert_eq!(latex("\"a & b\""), "\\ce{a & b}");
    }

    #[test]
    fn crystal_renders_the_tabulars() {
        // The formula-level \ce{} wraps the construct form too (the
        // reference's visit_formula is unconditional).
        assert_eq!(
            latex("crystal[NaCl](a=5.64,alpha=90){Na@f(0,0,0) Cl@f(0.5,0.5,0.5)}"),
            "\\ce{\\text{crystal}[NaCl]\\begin{tabular}{ll}\n\
             a & 5.64 \\\\\n\
             \\alpha & 90 \\\\\n\
             \\end{tabular}\
             \\begin{tabular}{l}\n\
             Na \\\\\n\
             Cl \\\\\n\
             \\end{tabular}}"
        );
    }

    #[test]
    fn embedded_math_stays_unsupported() {
        let error =
            render_formula(&crate::parse_nodes("`x^2`").expect("parses")).expect_err("unsupported");
        assert!(error.message.contains("Plurimath"));
    }
}
