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
    Ok(render_linear_formula(&nodes))
}

/// The linear (typographic) SVG for already-built nodes — the
/// grammar-independent seam (see [`render_formula`]).
pub fn render_linear_formula(nodes: &[Node]) -> String {
    let lines: Vec<Vec<Segment>> = nodes.iter().map(segments_for_node).collect();
    linear_document(&lines)
}

/// Renders the structural SVG document for already-built nodes — the
/// grammar-independent seam (`convert --from smiles|molfile`): SMILES
/// and molfile ingestion produce the same model, so the same SVG
/// renderers serve every input grammar.
pub fn render_formula(nodes: &[Node]) -> String {
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
            // The child's content is drawn in its own viewBox coordinates
            // (negative-origin for most layouts). Carrying the child's
            // viewBox maps that content into the placed viewport; without
            // it, raw coordinates paint outside the outer viewBox.
            format!(
                "<svg x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" viewBox=\"{} {} {} {}\" overflow=\"visible\">{}</svg>",
                round(x - bounds.0),
                round(y - bounds.1),
                round(*w),
                round(*h),
                round(*x),
                round(*y),
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
    let wedges = synthesize_stereo_wedges(&graph, &positions);
    let mut elements: Vec<String> = Vec::new();
    let mut drawn = std::collections::HashSet::new();

    for (bond_index, bond) in graph.bonds.iter().enumerate() {
        let key = if bond.from < bond.to {
            (bond.from, bond.to)
        } else {
            (bond.to, bond.from)
        };
        if !drawn.insert(key) {
            continue;
        }
        if let Some(wedge) = wedges.get(&bond_index) {
            // The synthesized glyph points centre → neighbour: narrow
            // at the stereocentre, widening (or hatching) toward the
            // ligand.
            let kind = if wedge.hash { "hash" } else { "wedge" };
            elements.push(bond_element(
                positions[wedge.center],
                positions[wedge.neighbor],
                &model::Bond {
                    kind: kind.to_string(),
                    direction: None,
                },
            ));
        } else {
            elements.push(bond_element(
                positions[bond.from],
                positions[bond.to],
                &bond.bond,
            ));
        }
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

// Tetrahedral wedge synthesis: a chirality-carrying atom (the SMILES
// @/@@ parity, TH1/TH2 in the model) with no authored wedge/hash of
// its own draws one synthesized wedge/hash so enantiomers render as
// distinguishable drawings. The model's parity is defined against the
// reference ligand order — the stated hydrogens first, then the
// neighbours in bond order (the SMILES writer's chirality order) —
// and the synthesized wedge realizes that hand against the laid-out
// neighbour geometry: the wedged ligand rises out of the plane (hash:
// sinks), the stated hydrogen sits at the tetrahedral completion
// opposite the drawn neighbours, and the wedge/hash choice is the one
// whose 3D orientation matches the parity. Synthesis declines (the
// drawing renders exactly as without it) rather than guesses: atoms
// whose ligand set is not tetrahedral, atoms with no claimable plain
// single bond, and degenerate (zero-volume) drawn geometries.
#[derive(Debug, Clone, Copy)]
pub struct SynthesizedWedge {
    pub center: usize,
    pub neighbor: usize,
    pub hash: bool,
}

pub fn synthesize_stereo_wedges(
    graph: &layout::LayoutGraph,
    positions: &[layout::Point],
) -> std::collections::HashMap<usize, SynthesizedWedge> {
    synthesize_stereo_wedges_in(graph, positions, false)
}

/// The mirrored-convention variant: the point set speaks the molfile
/// drawing convention (y rising) instead of the SVG's (y falling).
/// Mirroring the point set flips every orientation sign the core
/// reasons about, so the fitted targets invert with it — the glyphs
/// describe the same hands the unmirrored call would on the vertically
/// mirrored point set, which is what makes write∘parse∘write a byte
/// fixed point for chiral molecules (#86).
pub fn synthesize_stereo_wedges_mirrored(
    graph: &layout::LayoutGraph,
    positions: &[layout::Point],
) -> std::collections::HashMap<usize, SynthesizedWedge> {
    synthesize_stereo_wedges_in(graph, positions, true)
}

fn synthesize_stereo_wedges_in(
    graph: &layout::LayoutGraph,
    positions: &[layout::Point],
    invert: bool,
) -> std::collections::HashMap<usize, SynthesizedWedge> {
    let mut claimed: std::collections::HashMap<usize, SynthesizedWedge> =
        std::collections::HashMap::new();
    let mut used: std::collections::HashSet<usize> = std::collections::HashSet::new();

    for (index, atom) in graph.atoms.iter().enumerate() {
        let Some(parity) = atom.chirality.as_deref() else {
            continue;
        };
        if parity != "TH1" && parity != "TH2" {
            continue;
        }
        let stated_h = matches!(atom.hydrogens, Some(count) if count > 0);

        let mut neighbors: Vec<usize> = Vec::new();
        for bond in &graph.bonds {
            if bond.from == index {
                neighbors.push(bond.to);
            } else if bond.to == index {
                neighbors.push(bond.from);
            }
        }
        let mut distinct = neighbors.clone();
        distinct.sort_unstable();
        distinct.dedup();
        if distinct.len() != neighbors.len() {
            continue;
        }
        if neighbors.len() + usize::from(stated_h) != 4 {
            continue;
        }

        // Claim this centre's wedge on the first incident bond (bond
        // order) that is a plain single and not already carrying
        // another centre's synthesized wedge; authored wedge/hash
        // bonds are rendered verbatim and never claimed.
        let mut bond_index = None;
        let mut neighbor = None;
        for (i, bond) in graph.bonds.iter().enumerate() {
            if used.contains(&i) || bond.bond.kind != "single" {
                continue;
            }
            if bond.from != index && bond.to != index {
                continue;
            }
            bond_index = Some(i);
            neighbor = Some(if bond.from == index {
                bond.to
            } else {
                bond.from
            });
            break;
        }
        let Some(bond_index) = bond_index else {
            continue;
        };
        let neighbor = neighbor.expect("incident bond");
        let Some(wedge_slot) = neighbors.iter().position(|n| *n == neighbor) else {
            continue;
        };

        let center = positions[index];
        // Tetrahedral completion: the stated H sits opposite the drawn
        // neighbours' centroid.
        let hx = -neighbors
            .iter()
            .map(|n| positions[*n].x - center.x)
            .sum::<f64>();
        let hy = -neighbors
            .iter()
            .map(|n| positions[*n].y - center.y)
            .sum::<f64>();

        // Ligand points in reference order under a plain wedge (the
        // wedged ligand at +z toward the viewer, the stated H
        // opposite); a hash flips both out-of-plane signs, which
        // inverts the orientation exactly once.
        let mut points: Vec<(f64, f64, f64)> = Vec::new();
        if stated_h {
            points.push((hx, hy, -1.0));
        }
        for (i, n) in neighbors.iter().enumerate() {
            let z = if i == wedge_slot { 1.0 } else { 0.0 };
            points.push((positions[*n].x - center.x, positions[*n].y - center.y, z));
        }
        let sign = orientation(&points).signum();
        if sign == 0.0 {
            continue;
        }

        // The fitted polarity (tests/stereo_depiction.rs): the
        // stated-H construction (two out-of-plane ligands) and the
        // four-drawn-ligand construction (one apex) are different
        // point sets, so TH1's orientation target differs per
        // construction — each fitted against the engine, TH2 negated.
        // The mirrored convention inverts every target.
        let (th1_stated, th1_drawn) = if invert { (-1.0, 1.0) } else { (1.0, -1.0) };
        let target: f64 = match (parity, stated_h) {
            ("TH1", true) => th1_stated,
            ("TH1", false) => th1_drawn,
            ("TH2", true) => -th1_stated,
            _ => -th1_drawn,
        };
        used.insert(bond_index);
        claimed.insert(
            bond_index,
            SynthesizedWedge {
                center: index,
                neighbor,
                hash: sign != target,
            },
        );
    }

    claimed
}

// The orientation sign of a tetrahedral ligand set: the scalar triple
// product of the first three ligand vectors against the fourth.
fn orientation(points: &[(f64, f64, f64)]) -> f64 {
    let (a, b, c, d) = (points[0], points[1], points[2], points[3]);
    let u = [b.0 - a.0, b.1 - a.1, b.2 - a.2];
    let v = [c.0 - a.0, c.1 - a.1, c.2 - a.2];
    let w = [d.0 - a.0, d.1 - a.1, d.2 - a.2];
    (u[1] * v[2] - u[2] * v[1]) * w[0]
        + (u[2] * v[0] - u[0] * v[2]) * w[1]
        + (u[0] * v[1] - u[1] * v[0]) * w[2]
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
        Node::Adduct(adduct) => {
            let mut segments: Vec<Segment> = Vec::new();
            for (index, component) in adduct.components.iter().enumerate() {
                if index > 0 {
                    segments.push(seg(" · ", 0.0));
                }
                let owned = Node::Molecule(component.clone());
                segments.extend(segments_for_node(&owned));
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
            // The phase spelling, tight to the units.
            if let Some(spelling) = mol.phase_spelling() {
                segments.push(seg(spelling, 0.0));
            }
            // The unicode phase arrow, the glyph MathML typesets.
            if let Some(state) = &mol.state {
                let glyph = if state == "precipitate" { "↓" } else { "↑" };
                segments.push(seg(" ", 0.0));
                segments.push(seg(glyph, 0.0));
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
            // The reference renders constructs without a natural 2D
            // form through the canonical text (Formatter.render(:text,
            // mechanism)) — the full multi-line entry set, as with
            // crystals.
            vec![seg(node.to_string(), 0.0)]
        }
        Node::Crystal(_) => {
            // The reference renders constructs without a natural 2D
            // form through the canonical text (Formatter.render(:text,
            // crystal)) — one line here.
            vec![seg(node.to_string(), 0.0)]
        }
        Node::Spectrum(_) | Node::Calculation(_) | Node::ZMatrix(_) => {
            // The reference renders constructs without a natural 2D
            // form through the canonical text (Formatter.render(:text,
            // spectrum/calculation/z-matrix)) — the full multi-line
            // entry set, as with mechanisms.
            vec![seg(node.to_string(), 0.0)]
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
    if let Some(z) = &atom.atomic_number {
        segments.push(seg(z, -6.0));
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
        "donation" => "⇌⁺",
        "removal" => "⇌⁻",
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
