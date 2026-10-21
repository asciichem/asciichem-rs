//! SVG renderers — the Rust port of the TypeScript implementation's
//! `formatter/svg.ts` (linear, typographic) and
//! `formatter/structural-svg.ts` (2D skeletal diagram). Bonded
//! molecules render structurally; everything else renders linearly,
//! reactions and beyond-formulas nodes included.

use crate::layout;
use crate::model::{self, Node};

const FONT_SIZE: f64 = 15.0;
const LABEL_RADIUS: f64 = 11.0; // bond lines stop this far from atom centers
const LINE_FONT_SIZE: f64 = 16.0;
const CHAR_WIDTH: f64 = LINE_FONT_SIZE * 0.62; // monospace-ish estimate
const LINE_HEIGHT: f64 = LINE_FONT_SIZE * 1.6;

/// Parses AsciiChem text and renders the structural SVG document
/// (linear typographic fallback inside). Constructs outside the model
/// subset raise `ParseError`.
pub fn render(input: &str) -> Result<String, crate::ParseError> {
    let nodes = crate::parse_nodes(input)?;
    Ok(render_formula(&nodes))
}

/// Parses AsciiChem text and renders the linear (typographic) SVG.
pub fn render_linear(input: &str) -> Result<String, crate::ParseError> {
    let nodes = crate::parse_nodes(input)?;
    let lines: Vec<Vec<Segment>> = nodes.iter().map(segments_for_node).collect();
    Ok(linear_document(&lines))
}

// -- structural -------------------------------------------------------

fn render_formula(nodes: &[Node]) -> String {
    let children: Vec<String> = nodes
        .iter()
        .map(|node| match node {
            Node::Molecule(mol) => render_molecule(mol),
            other => render_linear_node(other),
        })
        .collect();
    let boxes: Vec<(f64, f64, f64, f64)> = children.iter().map(|c| view_box_of(c)).collect();
    let mut bounds = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for (x, y, w, h) in &boxes {
        bounds.0 = bounds.0.min(*x);
        bounds.1 = bounds.1.min(*y);
        bounds.2 = bounds.2.max(x + w);
        bounds.3 = bounds.3.max(y + h);
    }
    if children.is_empty() {
        bounds = (0.0, 0.0, 40.0, 40.0);
    }
    let inner: Vec<String> = children
        .iter()
        .zip(&boxes)
        .map(|(svg, (x, y, w, h))| {
            format!(
                "<svg x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" overflow=\"visible\">{}</svg>",
                round(x - bounds.0),
                round(y - bounds.1),
                round(*w),
                round(*h),
                strip_outer(svg)
            )
        })
        .collect();
    let width = (bounds.2 - bounds.0).max(40.0);
    let height = (bounds.3 - bounds.1).max(40.0);
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\" role=\"img\">\n{}\n</svg>",
        round(width),
        round(height),
        round(width),
        round(height),
        inner
            .iter()
            .map(|s| format!("    {s}"))
            .collect::<Vec<_>>()
            .join("\n")
    )
}

fn render_molecule(mol: &model::Molecule) -> String {
    let has_bond = mol.nodes.iter().any(|n| matches!(n, Node::Bond(_)));
    if !has_bond {
        // No authored connectivity: a formula is not a structure (the
        // isomerism point) — render typographically.
        let owned = Node::Molecule(mol.clone());
        return linear_document(&[segments_for_node(&owned)]);
    }
    let (graph, positions) = layout::layout_molecule(mol);
    let mut elements: Vec<String> = Vec::new();
    let mut drawn = std::collections::HashSet::new();

    for bond in &graph.bonds {
        let key = if bond.from < bond.to {
            (bond.from, bond.to)
        } else {
            (bond.to, bond.from)
        };
        if !drawn.insert(key) {
            continue;
        }
        elements.push(bond_element(
            positions[bond.from],
            positions[bond.to],
            &bond.bond,
        ));
    }

    for (i, atom) in graph.atoms.iter().enumerate() {
        elements.push(atom_element(atom, positions[i]));
    }

    for i in layout::unclosed_ring_atoms(&graph.atoms) {
        let pos = positions[i];
        elements.push(format!(
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#b45309\" font-family=\"monospace\">?</text>",
            round(pos.x + 14.0),
            round(pos.y - 10.0)
        ));
    }

    let mut xs: Vec<f64> = positions.iter().map(|p| p.x).collect();
    xs.push(free_label_extent(&positions));
    let ys: Vec<f64> = positions.iter().map(|p| p.y).collect();
    let min_x = xs.iter().cloned().fold(f64::INFINITY, f64::min) - 24.0;
    let min_y = ys.iter().cloned().fold(f64::INFINITY, f64::min) - 24.0;
    let max_x = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 24.0;
    let max_y = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 24.0;
    wrap(elements, min_x, min_y, max_x, max_y)
}

fn render_linear_node(node: &Node) -> String {
    linear_document(&[segments_for_node(node)])
}

#[allow(clippy::too_many_lines)]
fn bond_element(from: layout::Point, to: layout::Point, bond: &model::Bond) -> String {
    let dx = to.x - from.x;
    let dy = to.y - from.y;
    let len = dx.hypot(dy);
    let len = if len == 0.0 { 1.0 } else { len };
    let ux = dx / len;
    let uy = dy / len;
    let a = layout::Point {
        x: from.x + ux * LABEL_RADIUS,
        y: from.y + uy * LABEL_RADIUS,
    };
    let b = layout::Point {
        x: to.x - ux * LABEL_RADIUS,
        y: to.y - uy * LABEL_RADIUS,
    };
    let nx = -uy;
    let ny = ux;
    let mut lines: Vec<String> = Vec::new();

    match bond.kind.as_str() {
        "single" => lines.push(line_element(&a, &b)),
        "aromatic" => lines.push(line_element(&a, &b).replace(
            "stroke-width=\"1.6\"",
            "stroke-width=\"1.6\" stroke-dasharray=\"4 2.5\"",
        )),
        "double" => {
            lines.push(offset_line(&a, &b, nx, ny, -2.6));
            lines.push(offset_line(&a, &b, nx, ny, 2.6));
        }
        "triple" => {
            lines.push(line_element(&a, &b));
            lines.push(offset_line(&a, &b, nx, ny, -3.4));
            lines.push(offset_line(&a, &b, nx, ny, 3.4));
        }
        "quadruple" => {
            for off in [-4.8, -1.6, 1.6, 4.8] {
                lines.push(offset_line(&a, &b, nx, ny, off));
            }
        }
        "wedge" => lines.push(format!(
            "<polygon points=\"{},{} {},{} {},{}\" fill=\"currentColor\"/>",
            round(from.x),
            round(from.y),
            round(b.x + nx * 4.0),
            round(b.y + ny * 4.0),
            round(b.x - nx * 4.0),
            round(b.y - ny * 4.0)
        )),
        "hash" => {
            let mut t = 0.25;
            while t <= 0.95 {
                let w = 1.5 + 3.5 * t;
                let px = a.x + (b.x - a.x) * t;
                let py = a.y + (b.y - a.y) * t;
                lines.push(format!(
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"currentColor\" stroke-width=\"1.2\"/>",
                    round(px + nx * w),
                    round(py + ny * w),
                    round(px - nx * w),
                    round(py - ny * w)
                ));
                t += 0.18;
            }
        }
        "dative" => {
            lines.push(format!(
                "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"currentColor\" stroke-width=\"1.6\"/>",
                round(a.x),
                round(a.y),
                round(b.x),
                round(b.y)
            ));
            lines.push(arrow_head(b, ux, uy));
        }
        "wavy" => lines.push(wavy_path(&a, &b)),
        _ => lines.push(line_element(&a, &b)),
    }
    lines.join("")
}

fn arrow_head(b: layout::Point, ux: f64, uy: f64) -> String {
    let size = 6.0;
    let left = layout::Point {
        x: b.x - ux * size - uy * size * 0.6,
        y: b.y - uy * size + ux * size * 0.6,
    };
    let right = layout::Point {
        x: b.x - ux * size + uy * size * 0.6,
        y: b.y - uy * size - ux * size * 0.6,
    };
    format!(
        "<polyline points=\"{},{} {},{} {},{}\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.6\"/>",
        round(left.x),
        round(left.y),
        round(b.x),
        round(b.y),
        round(right.x),
        round(right.y)
    )
}

fn wavy_path(a: &layout::Point, b: &layout::Point) -> String {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len = dx.hypot(dy);
    let len = if len == 0.0 { 1.0 } else { len };
    let ux = dx / len;
    let uy = dy / len;
    let nx = -uy;
    let ny = ux;
    let waves = (len / 12.0).round().max(2.0) as usize;
    let mut d = format!("M {} {}", round(a.x), round(a.y));
    for i in 1..=waves {
        let t = i as f64 / waves as f64;
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        let amp = sign * 3.5 * (std::f64::consts::PI * t).sin();
        d.push_str(&format!(
            " L {} {}",
            round(a.x + dx * t + nx * amp),
            round(a.y + dy * t + ny * amp)
        ));
    }
    format!("<path d=\"{d}\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.6\"/>")
}

fn line_element(a: &layout::Point, b: &layout::Point) -> String {
    format!(
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"currentColor\" stroke-width=\"1.6\"/>",
        round(a.x),
        round(a.y),
        round(b.x),
        round(b.y)
    )
}

fn offset_line(a: &layout::Point, b: &layout::Point, nx: f64, ny: f64, off: f64) -> String {
    line_element(
        &layout::Point {
            x: a.x + nx * off,
            y: a.y + ny * off,
        },
        &layout::Point {
            x: b.x + nx * off,
            y: b.y + ny * off,
        },
    )
}

fn atom_element(atom: &model::Atom, pos: layout::Point) -> String {
    let mut label = escape_xml(&atom.element);
    if let Some(isotope) = &atom.isotope {
        label = format!(
            "<tspan font-size=\"10\" dy=\"-4\">{}</tspan><tspan dy=\"4\">{}</tspan>",
            escape_xml(isotope),
            escape_xml(&atom.element)
        );
    }
    if let Some(charge) = &atom.charge {
        label.push_str(&format!(
            "<tspan font-size=\"10\" dy=\"-5\">{}</tspan>",
            escape_xml(charge)
        ));
    }
    format!(
        "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"{}\" font-weight=\"500\" fill=\"currentColor\">{}</text>",
        round(pos.x),
        round(pos.y + 5.0),
        FONT_SIZE,
        label
    )
}

fn free_label_extent(positions: &[layout::Point]) -> f64 {
    positions
        .iter()
        .map(|p| p.x + 12.0)
        .fold(f64::NEG_INFINITY, f64::max)
}

fn view_box_of(svg: &str) -> (f64, f64, f64, f64) {
    let Some(start) = svg.find("viewBox=\"") else {
        return (0.0, 0.0, 40.0, 40.0);
    };
    let rest = &svg[start + "viewBox=\"".len()..];
    let Some(end) = rest.find('"') else {
        return (0.0, 0.0, 40.0, 40.0);
    };
    let parts: Vec<f64> = rest[..end]
        .split_whitespace()
        .filter_map(|p| p.parse().ok())
        .collect();
    match parts.as_slice() {
        [x, y, w, h] => (*x, *y, *w, *h),
        _ => (0.0, 0.0, 40.0, 40.0),
    }
}

fn strip_outer(svg: &str) -> String {
    let Some(open_end) = svg.find('>') else {
        return svg.to_string();
    };
    match svg.rfind("</svg>") {
        Some(close) => svg[open_end + 1..close].trim().to_string(),
        None => svg.to_string(),
    }
}

fn wrap(elements: Vec<String>, min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> String {
    let width = (max_x - min_x).max(40.0);
    let height = (max_y - min_y).max(40.0);
    let body = elements
        .iter()
        .map(|e| format!("    {e}"))
        .collect::<Vec<_>>()
        .join("\n");
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{} {} {} {}\" width=\"{}\" height=\"{}\" role=\"img\">\n  <g fill=\"none\" stroke=\"currentColor\">\n{}\n  </g>\n</svg>",
        round(min_x),
        round(min_y),
        round(width),
        round(height),
        round(width),
        round(height),
        body
    )
}

// -- linear -------------------------------------------------------------

struct Segment {
    text: String,
    dy: f64,
    size: f64,
}

fn seg(text: impl Into<String>, dy: f64) -> Segment {
    Segment {
        text: text.into(),
        dy,
        size: if dy == 0.0 {
            LINE_FONT_SIZE
        } else {
            LINE_FONT_SIZE * 0.72
        },
    }
}

fn segments_for_node(node: &Node) -> Vec<Segment> {
    match node {
        Node::Atom(atom) => atom_segments(atom),
        Node::Group(group) => {
            let (open, close) = group.chars();
            let mut segments = vec![seg(open, 0.0)];
            for child in &group.nodes {
                segments.extend(segments_for_node(child));
            }
            segments.push(seg(close, 0.0));
            if let Some(mult) = &group.multiplicity {
                segments.push(seg(mult, 4.0));
            }
            segments
        }
        Node::Molecule(mol) => {
            let mut segments: Vec<Segment> = Vec::new();
            if let Some(letter) = mol.stereo_letter() {
                segments.push(seg(format!("({letter})-"), 0.0));
            }
            if let Some(coef) = &mol.coefficient {
                segments.push(seg(coef, 0.0));
            }
            for child in &mol.nodes {
                segments.extend(segments_for_node(child));
            }
            segments
        }
        Node::Bond(bond) => vec![seg(bond_glyph(bond.kind.as_str()), 0.0)],
        Node::Reaction(reaction) => {
            let mut segments: Vec<Segment> = Vec::new();
            for term in &reaction.reactants {
                segments.extend(segments_for_node(term));
            }
            segments.push(seg(" ", 0.0));
            segments.push(seg(arrow_glyph(reaction.arrow.as_str()), 0.0));
            push_conditions(&mut segments, reaction.conditions.as_ref());
            segments.push(seg(" ", 0.0));
            for term in &reaction.products {
                segments.extend(segments_for_node(term));
            }
            segments
        }
        Node::ReactionCascade(cascade) => {
            let mut segments: Vec<Segment> = Vec::new();
            let mut steps = cascade.steps.iter();
            let Some(head) = steps.next() else {
                return segments;
            };
            for term in &head.reactants {
                segments.extend(segments_for_node(term));
            }
            for step in std::iter::once(head).chain(steps) {
                segments.push(seg(" ", 0.0));
                segments.push(seg(arrow_glyph(step.arrow.as_str()), 0.0));
                push_conditions(&mut segments, step.conditions.as_ref());
                segments.push(seg(" ", 0.0));
                for term in &step.products {
                    segments.extend(segments_for_node(term));
                }
            }
            segments
        }
        Node::ElectronConfiguration(ec) => {
            let mut segments: Vec<Segment> = Vec::new();
            for (i, (orb, occ)) in ec.orbitals.iter().enumerate() {
                if i > 0 {
                    segments.push(seg(" ", 0.0));
                }
                segments.push(seg(orb, 0.0));
                segments.push(seg(occ, -6.0));
            }
            segments
        }
        Node::EmbeddedMath(source) => vec![seg(format!("`{source}`"), 0.0)],
        Node::Text(content) => vec![seg(format!("\"{content}\""), 0.0)],
        Node::Mechanism(_) => {
            // Linear SVG prints the first line of the canonical text.
            let text = node.to_string();
            vec![seg(text.lines().next().unwrap_or_default(), 0.0)]
        }
    }
}

fn push_conditions(segments: &mut Vec<Segment>, conditions: Option<&model::Conditions>) {
    let Some(conds) = conditions else {
        return;
    };
    if let Some(above) = &conds.above {
        segments.push(seg(format!("[{above}]"), -8.0));
    }
    if let Some(below) = &conds.below {
        segments.push(seg(format!("[{below}]"), 12.0));
    }
}

fn atom_segments(atom: &model::Atom) -> Vec<Segment> {
    let mut segments = Vec::new();
    if let Some(lp) = atom.lone_pairs {
        segments.push(seg(":".repeat(lp as usize), -6.0));
    }
    if let Some(isotope) = &atom.isotope {
        segments.push(seg(isotope, -6.0));
    }
    segments.push(seg(&atom.element, 0.0));
    if let Some(sub) = &atom.subscript {
        segments.push(seg(sub, 5.0));
    }
    if let Some(sup) = &atom.superscript {
        segments.push(seg(sup, -6.0));
    }
    if let Some(charge) = &atom.charge {
        segments.push(seg(charge, -6.0));
    }
    if let Some(ox) = &atom.oxidation_state {
        segments.push(seg(format!("({ox})"), -6.0));
    }
    if let Some(re_) = atom.radical_electrons {
        segments.push(seg("·".repeat(re_ as usize), -6.0));
    }
    if let Some(rc) = &atom.ring_closures {
        segments.push(seg(rc, 0.0));
    }
    segments
}

fn linear_document(lines: &[Vec<Segment>]) -> String {
    let tspans: Vec<String> = lines
        .iter()
        .enumerate()
        .map(|(line_index, segments)| {
            let baseline = 20.0 + line_index as f64 * LINE_HEIGHT;
            let mut x = 4.0;
            let mut out = String::new();
            for s in segments {
                out.push_str(&format!(
                    "<tspan x=\"{}\" y=\"{}\" font-size=\"{}\">{}</tspan>",
                    round(x),
                    round(baseline + s.dy),
                    round(s.size),
                    escape_xml(&s.text)
                ));
                // Width counts characters (the JS reference's
                // .length), never bytes — bond glyphs are multibyte.
                x += s.text.chars().count() as f64 * CHAR_WIDTH * (s.size / LINE_FONT_SIZE);
            }
            format!(
                "    <text x=\"4\" y=\"{}\" font-family=\"monospace\" font-size=\"{}\">{out}</text>",
                round(baseline),
                LINE_FONT_SIZE
            )
        })
        .collect();
    let width = lines
        .iter()
        .map(|segments| {
            segments
                .iter()
                .map(|s| s.text.chars().count() as f64 * CHAR_WIDTH)
                .sum::<f64>()
        })
        .fold(40.0, f64::max);
    let height = lines.len().max(1) as f64 * LINE_HEIGHT + 12.0;
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\" role=\"img\">\n{}\n  </svg>",
        round(width),
        round(height),
        round(width),
        round(height),
        tspans.join("\n")
    )
}

fn bond_glyph(kind: &str) -> &'static str {
    match kind {
        "double" => "═",
        "triple" => "≡",
        "quadruple" => "≣",
        "wedge" => "▲–",
        "hash" => "–▨",
        "dative" => "→",
        "wavy" => "∼",
        "aromatic" => ":",
        _ => "–",
    }
}

fn arrow_glyph(kind: &str) -> &'static str {
    match kind {
        "reverse" => "←",
        "equilibrium" => "⇌",
        "resonance" => "↔",
        _ => "→",
    }
}

fn round(n: f64) -> impl std::fmt::Display {
    let rounded = (n * 10.0).round() / 10.0;
    // Rust prints whole f64 values without a trailing ".0"; the JS
    // reference does the same via number-to-string.
    rounded
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
