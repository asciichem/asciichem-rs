//! Presentation MathML renderer — the Rust counterpart of the
//! reference `Formatter::Mathml`, byte-exact against the corpus L2
//! goldens (2-space indent, XML declaration, `<none/>` compaction).
//! Isotope is a LEFT superscript via `mmultiscripts` so `^14C` keeps
//! the binding on the atom.

use crate::model::{Bond, Node, Reaction};

/// Parses AsciiChem text and renders presentation MathML, byte-exact
/// with the corpus L2 goldens. Constructs outside the model subset
/// raise `ParseError` with an `unsupported` message.
pub fn render(input: &str) -> Result<String, crate::ParseError> {
    let source = input.trim().to_string();
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("asciichem-parse".into())
            .stack_size(crate::PARSE_STACK_BYTES)
            .spawn_scoped(scope, || render_inner(&source))
            .expect("spawn parser thread")
            .join()
            .expect("parser thread panicked")
    })
}

fn render_inner(source: &str) -> Result<String, crate::ParseError> {
    let nodes = crate::parse_nodes(source)?;
    let mut math = Element::new("math");
    math.attrs.push((
        "xmlns".to_string(),
        "http://www.w3.org/1998/Math/MathML".to_string(),
    ));
    let mut mrow = Element::new("mrow");
    for node in &nodes {
        mrow.children
            .push(render_node(node).map_err(|error| crate::ParseError {
                message: error.to_string(),
            })?);
    }
    math.children.push(Node_::Element(mrow));
    Ok(compact_empty(&math.serialize()))
}

enum Node_ {
    Element(Element),
    Text(String),
}

struct Element {
    name: String,
    attrs: Vec<(String, String)>,
    children: Vec<Node_>,
}

impl Element {
    fn new(name: &str) -> Self {
        Element {
            name: name.to_string(),
            attrs: vec![],
            children: vec![],
        }
    }

    fn el(name: &str, child: Node_) -> Element {
        let mut element = Element::new(name);
        element.children.push(child);
        element
    }

    fn text(name: &str, content: impl std::fmt::Display) -> Element {
        Element::el(name, Node_::Text(content.to_string()))
    }

    fn mi(content: impl std::fmt::Display) -> Element {
        let mut element = Element::text("mi", content);
        element
            .attrs
            .push(("mathvariant".to_string(), "normal".to_string()));
        element
    }

    fn serialize(&self) -> String {
        let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        write_element(&mut out, self, 0);
        out
    }
}

fn write_element(out: &mut String, element: &Element, depth: usize) {
    out.push_str(&"  ".repeat(depth));
    out.push('<');
    out.push_str(&element.name);
    for (key, value) in &element.attrs {
        out.push_str(&format!(" {key}=\"{value}\""));
    }
    // A lone text child renders inline (`<mi>H</mi>`), matching the
    // reference serializer; element children get indented lines.
    if let [Node_::Text(only)] = &element.children[..] {
        out.push('>');
        out.push_str(&escape_text(only));
        out.push_str("</");
        out.push_str(&element.name);
        out.push_str(">\n");
        return;
    }
    if element.children.is_empty() {
        out.push_str("></");
        out.push_str(&element.name);
        out.push_str(">\n");
        return;
    }
    out.push_str(">\n");
    for child in &element.children {
        match child {
            Node_::Text(text) => {
                out.push_str(&"  ".repeat(depth + 1));
                out.push_str(&escape_text(text));
                out.push('\n');
            }
            Node_::Element(child) => write_element(out, child, depth + 1),
        }
    }
    out.push_str(&"  ".repeat(depth));
    out.push_str("</");
    out.push_str(&element.name);
    out.push_str(">\n");
}

fn escape_text(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

// The reference compacts only these two empty elements back to
// self-closed form; every other empty element stays expanded.
fn compact_empty(xml: &str) -> String {
    xml.replace("<none></none>", "<none/>")
        .replace("<mprescripts></mprescripts>", "<mprescripts/>")
}

type Render<T> = std::result::Result<T, crate::build::BuildError>;

fn render_node(node: &Node) -> Render<Node_> {
    match node {
        Node::Atom(atom) => Ok(Node_::Element(render_atom(atom))),
        Node::Bond(bond) => Ok(Node_::Element(Element::text("mo", bond.entity()))),
        Node::Group(group) => Ok(Node_::Element(render_group(group)?)),
        Node::Molecule(mol) => Ok(Node_::Element(render_molecule(mol)?)),
        Node::Reaction(reaction) => Ok(Node_::Element(render_reaction(reaction)?)),
        Node::ReactionCascade(cascade) => Ok(Node_::Element(render_cascade(cascade)?)),
        Node::ElectronConfiguration(ec) => Ok(Node_::Element(render_electron_configuration(ec))),
        Node::Mechanism(mech) => Ok(Node_::Element(render_mechanism(mech))),
        Node::EmbeddedMath(_) => Err(crate::build::BuildError::Unsupported(
            "embedded math rendering (Plurimath pipeline)".into(),
        )),
        Node::Text(content) => Ok(Node_::Element(Element::text("mtext", content))),
    }
}

fn render_molecule(mol: &crate::model::Molecule) -> Render<Element> {
    let mut mrow = Element::new("mrow");
    if let Some(letter) = mol.stereo_letter() {
        mrow.children.push(Node_::Element(Element::text(
            "mtext",
            format!("({letter})-"),
        )));
    }
    if let Some(coef) = &mol.coefficient {
        mrow.children
            .push(Node_::Element(Element::text("mn", coef)));
    }
    for node in &mol.nodes {
        mrow.children.push(render_node(node)?);
    }
    Ok(mrow)
}

fn render_atom(atom: &crate::model::Atom) -> Element {
    let mut base = Element::mi(&atom.element);
    if let Some(isotope) = &atom.isotope {
        let mut multi = Element::new("mmultiscripts");
        multi.children.push(Node_::Element(base));
        multi.children.push(Node_::Element(Element::new("none")));
        multi.children.push(Node_::Element(Element::new("none")));
        multi
            .children
            .push(Node_::Element(Element::new("mprescripts")));
        multi.children.push(Node_::Element(Element::new("none")));
        multi
            .children
            .push(Node_::Element(Element::text("mn", isotope)));
        base = multi;
    }
    if let Some(lp) = atom.lone_pairs {
        let mut mrow = Element::new("mrow");
        mrow.children.push(Node_::Element(Element::text(
            "mtext",
            ":".repeat(lp as usize),
        )));
        mrow.children.push(Node_::Element(base));
        base = mrow;
    }
    let has_sub = atom.subscript.as_deref().is_some_and(|s| !s.is_empty());
    match (has_sub, super_element(atom)) {
        (true, None) => {
            base = wrap(
                "msub",
                base,
                Element::text("mn", atom.subscript.as_deref().unwrap_or_default()),
            );
        }
        (false, Some(super_node)) => {
            base = wrap("msup", base, super_node);
        }
        (true, Some(super_node)) => {
            let mut msubsup = Element::new("msubsup");
            msubsup.children.push(Node_::Element(base));
            msubsup.children.push(Node_::Element(Element::text(
                "mn",
                atom.subscript.as_deref().unwrap_or_default(),
            )));
            msubsup.children.push(Node_::Element(super_node));
            base = msubsup;
        }
        (false, None) => {}
    }
    if let Some(re_) = atom.radical_electrons {
        let mut mrow = Element::new("mrow");
        mrow.children.push(Node_::Element(base));
        mrow.children.push(Node_::Element(Element::text(
            "mtext",
            ".".repeat(re_ as usize),
        )));
        base = mrow;
    }
    if let Some(rc) = &atom.ring_closures {
        let mut mrow = Element::new("mrow");
        mrow.children.push(Node_::Element(base));
        mrow.children.push(Node_::Element(Element::text("mn", rc)));
        base = mrow;
    }
    base
}

fn wrap(name: &str, base: Element, script: Element) -> Element {
    let mut element = Element::new(name);
    element.children.push(Node_::Element(base));
    element.children.push(Node_::Element(script));
    element
}

// The superscript-style marker: charge > oxidation state > raw
// superscript. None when the atom carries neither.
fn super_element(atom: &crate::model::Atom) -> Option<Element> {
    if let Some(charge) = &atom.charge {
        if let Some(row) = charge_row(charge) {
            return Some(row);
        }
    }
    if let Some(ox) = &atom.oxidation_state {
        let mut mrow = Element::new("mrow");
        mrow.children.push(Node_::Element(Element::mi(ox)));
        return Some(mrow);
    }
    atom.superscript
        .as_ref()
        .map(|sup| Element::text("mn", sup))
}

fn charge_row(charge: &str) -> Option<Element> {
    let bytes = charge.as_bytes();
    let (n, sign) = match (bytes.first(), bytes.last()) {
        (Some(first), Some(last)) if (*last == b'+' || *last == b'-') && first.is_ascii_digit() => {
            (charge.split_at(charge.len() - 1).0, *last)
        }
        (Some(first), Some(last))
            if (*first == b'+' || *first == b'-') && last.is_ascii_digit() =>
        {
            (&charge[1..], *first)
        }
        (Some(only), _) if *only == b'+' || *only == b'-' => ("", *only),
        _ => return None,
    };
    let mut mrow = Element::new("mrow");
    if !n.is_empty() {
        mrow.children.push(Node_::Element(Element::text("mn", n)));
    }
    mrow.children
        .push(Node_::Element(Element::text("mo", sign as char)));
    Some(mrow)
}

fn render_group(group: &crate::model::Group) -> Render<Element> {
    let (open, close) = group.chars();
    let mut mrow = Element::new("mrow");
    mrow.children
        .push(Node_::Element(Element::text("mo", open)));
    for node in &group.nodes {
        mrow.children.push(render_node(node)?);
    }
    mrow.children
        .push(Node_::Element(Element::text("mo", close)));
    Ok(match &group.multiplicity {
        None => mrow,
        Some(mult) => wrap("msub", mrow, Element::text("mn", mult)),
    })
}

fn render_reaction(reaction: &Reaction) -> Render<Element> {
    let mut mrow = Element::new("mrow");
    add_terms(&mut mrow, &reaction.reactants)?;
    mrow.children.push(render_arrow(reaction));
    add_terms(&mut mrow, &reaction.products)?;
    Ok(mrow)
}

fn render_cascade(cascade: &crate::model::ReactionCascade) -> Render<Element> {
    let mut mrow = Element::new("mrow");
    let Some(head) = cascade.steps.first() else {
        return Ok(mrow);
    };
    add_terms(&mut mrow, &head.reactants)?;
    for step in &cascade.steps {
        mrow.children.push(render_arrow(step));
        add_terms(&mut mrow, &step.products)?;
    }
    Ok(mrow)
}

fn render_electron_configuration(ec: &crate::model::ElectronConfiguration) -> Element {
    let mut mrow = Element::new("mrow");
    for (index, (orbital, occupancy)) in ec.orbitals.iter().enumerate() {
        if index > 0 {
            // Orbital separators are non-breaking spaces (IUPAC
            // style), matching the reference formatter.
            mrow.children
                .push(Node_::Element(Element::text("mo", "\u{00A0}")));
        }
        let mut msup = Element::new("msup");
        msup.children.push(Node_::Element(Element::mi(orbital)));
        msup.children
            .push(Node_::Element(Element::text("mn", occupancy)));
        mrow.children.push(Node_::Element(msup));
    }
    mrow
}

fn render_mechanism(mech: &crate::model::Mechanism) -> Element {
    let mut mrow = Element::new("mrow");
    mrow.children.push(Node_::Element(Element::mi("mechanism")));
    if mech.steps.is_empty() && mech.spectators.is_empty() {
        return mrow;
    }
    let mut rows: Vec<Vec<Node_>> = mech
        .steps
        .iter()
        .map(|step| {
            vec![
                Node_::Element(Element::mi(&step.label)),
                Node_::Element(Element::text("mtext", &step.reaction)),
            ]
        })
        .collect();
    for spectator in &mech.spectators {
        rows.push(vec![
            Node_::Element(Element::mi("spectator")),
            Node_::Element(Element::text("mtext", spectator)),
        ]);
    }
    mrow.children.push(Node_::Element(simple_table(rows)));
    mrow
}

fn add_terms(parent: &mut Element, terms: &[Node]) -> Render<()> {
    for (index, term) in terms.iter().enumerate() {
        if index > 0 {
            parent
                .children
                .push(Node_::Element(Element::text("mo", "+")));
        }
        parent.children.push(render_node(term)?);
    }
    Ok(())
}

fn render_arrow(reaction: &Reaction) -> Node_ {
    let op = Element::text("mo", reaction.arrow_entity());
    let Some(conds) = &reaction.conditions else {
        return Node_::Element(op);
    };
    let above = conds.above.as_deref();
    let below = conds.below.as_deref();
    let Some(above_or_below) = above.or(below) else {
        return Node_::Element(op);
    };
    let _ = above_or_below;
    let name = match (above, below) {
        (Some(_), Some(_)) => "munderover",
        (Some(_), None) => "mover",
        _ => "munder",
    };
    let mut element = Element::new(name);
    element.children.push(Node_::Element(op));
    if let Some(above) = above {
        element.children.push(render_condition(above));
    }
    if let Some(below) = below {
        element.children.push(render_condition(below));
    }
    Node_::Element(element)
}

// Conditions are re-parsed as AsciiChem so `_N`/`^N` render as real
// scripts; prose falls back to plain <mtext>.
fn render_condition(text: &str) -> Node_ {
    if text.is_empty() {
        return Node_::Element(Element::text("mtext", ""));
    }
    let rendered = crate::parse_nodes(text).ok().and_then(|nodes| {
        let mut mrow = Element::new("mrow");
        for node in &nodes {
            mrow.children.push(render_node(node).ok()?);
        }
        Some(Node_::Element(mrow))
    });
    rendered.unwrap_or_else(|| Node_::Element(Element::text("mtext", text)))
}

fn simple_table(rows: Vec<Vec<Node_>>) -> Element {
    let mut table = Element::new("mtable");
    for cells in rows {
        let mut tr = Element::new("mtr");
        for cell in cells {
            let mut td = Element::new("mtd");
            td.children.push(cell);
            tr.children.push(Node_::Element(td));
        }
        table.children.push(Node_::Element(tr));
    }
    table
}

impl Bond {
    fn entity(&self) -> &'static str {
        match self.kind.as_str() {
            "double" => "=",
            "triple" => "≡",
            "quadruple" => "≣",
            "wedge" => "↑",
            "hash" => "↓",
            "dative" => "→",
            "wavy" => "∼",
            "aromatic" => ":",
            _ => "-",
        }
    }
}

impl Reaction {
    fn arrow_entity(&self) -> &'static str {
        match self.arrow.as_str() {
            "reverse" => "←",
            "equilibrium" => "⇌",
            "resonance" => "↔",
            _ => "→",
        }
    }
}
