//! Semantic model subset (formulas, molecules, reactions, mechanisms)
//! with the canonical text renderer. The renderer mirrors the
//! reference Text formatter rule for rule; round-trip is the
//! contract.

use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: String,
    pub isotope: Option<String>,
    pub subscript: Option<String>,
    /// Raw superscript — set only when the superscript is neither a
    /// charge nor an oxidation state (those two carry the info).
    pub superscript: Option<String>,
    pub charge: Option<String>,
    pub oxidation_state: Option<String>,
    pub lone_pairs: Option<u32>,
    pub radical_electrons: Option<u32>,
    pub ring_closures: Option<String>,
    pub chirality: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    pub kind: String,
    pub direction: Option<String>,
}

impl Bond {
    fn ascii(&self) -> &'static str {
        match self.kind.as_str() {
            "double" => "=",
            "triple" => "#",
            "quadruple" => "##",
            "wedge" => ">-",
            "hash" => "-<",
            "dative" => "~>",
            "wavy" => "~~",
            "aromatic" => ":",
            _ => "-",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    pub nodes: Vec<Node>,
    pub multiplicity: Option<String>,
    pub bracket: String,
}

impl Group {
    fn chars(&self) -> (&'static str, &'static str) {
        match self.bracket.as_str() {
            "square" => ("[", "]"),
            "brace" => ("{", "}"),
            _ => ("(", ")"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Identifier {
    pub value: String,
    pub convention: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Name {
    pub content: String,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Molecule {
    pub nodes: Vec<Node>,
    pub coefficient: Option<String>,
    pub stereo: Option<String>,
    pub names: Vec<Name>,
    pub identifiers: Vec<Identifier>,
    pub title: Option<String>,
    /// `@formula("...")` concise formulas.
    pub formulas: Vec<String>,
    pub labels: Vec<String>,
    /// `@key("value")` annotations outside the known set.
    pub properties: Vec<(String, String)>,
    /// `@meta("name","content")` pairs.
    pub metadata: Vec<(String, String)>,
}

impl Molecule {
    fn stereo_letter(&self) -> Option<&str> {
        self.stereo.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Conditions {
    pub above: Option<String>,
    pub below: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    pub reactants: Vec<Node>,
    pub products: Vec<Node>,
    pub arrow: String,
    pub conditions: Option<Conditions>,
}

impl Reaction {
    fn arrow_ascii(&self) -> &'static str {
        match self.arrow.as_str() {
            "reverse" => "<-",
            "equilibrium" => "<=>",
            "resonance" => "<->",
            _ => "->",
        }
    }

    fn render_arrow(&self) -> String {
        let mut out = self.arrow_ascii().to_string();
        if let Some(conds) = &self.conditions {
            if let Some(above) = &conds.above {
                out.push_str(&format!("[{above}]"));
            }
            if let Some(below) = &conds.below {
                out.push_str(&format!("[{below}]"));
            }
        }
        out
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionCascade {
    pub steps: Vec<Reaction>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ElectronConfiguration {
    pub orbitals: Vec<(String, String)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Flow {
    pub source: String,
    pub target: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MechanismStep {
    pub label: String,
    pub reaction: String,
    pub flows: Vec<Flow>,
}

impl MechanismStep {
    fn text(&self) -> String {
        if self.flows.is_empty() {
            return self.reaction.clone();
        }
        let arrows: Vec<String> = self
            .flows
            .iter()
            .map(|f| format!("{} > {}", f.source, f.target))
            .collect();
        format!("{} | flow: {}", self.reaction, arrows.join("; "))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mechanism {
    pub steps: Vec<MechanismStep>,
    pub spectators: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Node {
    Atom(Atom),
    Bond(Bond),
    Group(Group),
    Molecule(Molecule),
    Reaction(Reaction),
    ReactionCascade(ReactionCascade),
    ElectronConfiguration(ElectronConfiguration),
    Mechanism(Mechanism),
    EmbeddedMath(String),
    Text(String),
}

fn render_terms(nodes: &[Node]) -> String {
    nodes
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(" + ")
}

fn render_molecule_annotations(mol: &Molecule) -> String {
    let mut parts: Vec<String> = vec![];
    for name in &mol.names {
        parts.push(format!("@name(\"{}\")", name.content));
    }
    for id in &mol.identifiers {
        parts.push(format!("@{}(\"{}\")", id.convention, id.value));
    }
    if let Some(title) = &mol.title {
        parts.push(format!("@title(\"{title}\")"));
    }
    for concise in &mol.formulas {
        parts.push(format!("@formula(\"{concise}\")"));
    }
    for label in &mol.labels {
        parts.push(format!("@label(\"{label}\")"));
    }
    for (title, value) in &mol.properties {
        parts.push(format!("@{title}(\"{value}\")"));
    }
    for (name, content) in &mol.metadata {
        parts.push(format!("@meta(\"{name}\",\"{content}\")"));
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!(" {}", parts.join(" "))
    }
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Node::Atom(atom) => {
                let mut parts: Vec<String> = vec![];
                if let Some(lp) = atom.lone_pairs {
                    parts.push(":".repeat(lp as usize));
                }
                if let Some(iso) = &atom.isotope {
                    parts.push(format!("^{iso}"));
                }
                parts.push(atom.element.clone());
                if let Some(sub) = &atom.subscript {
                    parts.push(format!("_{sub}"));
                }
                if let Some(sup) = &atom.superscript {
                    parts.push(format!("^{sup}"));
                }
                if let Some(chg) = &atom.charge {
                    parts.push(format!("^{chg}"));
                }
                if let Some(ox) = &atom.oxidation_state {
                    parts.push(format!("^({ox})"));
                }
                if let Some(re_) = atom.radical_electrons {
                    parts.push(".".repeat(re_ as usize));
                }
                if let Some(rc) = &atom.ring_closures {
                    parts.push(rc.clone());
                }
                if let Some(ch) = &atom.chirality {
                    parts.push(format!("@{ch}"));
                }
                write!(f, "{}", parts.concat())
            }
            Node::Bond(bond) => {
                if let Some(direction) = &bond.direction {
                    return write!(f, "{}", if direction == "up" { "/" } else { "\\" });
                }
                write!(f, "{}", bond.ascii())
            }
            Node::Group(group) => {
                let (open, close) = group.chars();
                let body: String = group.nodes.iter().map(|n| n.to_string()).collect();
                let suffix = group
                    .multiplicity
                    .as_ref()
                    .map(|m| format!("_{m}"))
                    .unwrap_or_default();
                write!(f, "{open}{body}{close}{suffix}")
            }
            Node::Molecule(mol) => {
                let mut out = String::new();
                if let Some(letter) = mol.stereo_letter() {
                    out.push_str(&format!("({letter})-"));
                }
                if let Some(coef) = &mol.coefficient {
                    out.push_str(coef);
                }
                out.push_str(&mol.nodes.iter().map(|n| n.to_string()).collect::<String>());
                out.push_str(&render_molecule_annotations(mol));
                write!(f, "{out}")
            }
            Node::Reaction(reaction) => {
                let left = render_terms(&reaction.reactants);
                let right = render_terms(&reaction.products);
                write!(f, "{left} {} {right}", reaction.render_arrow())
            }
            Node::ReactionCascade(cascade) => {
                let mut out = String::new();
                for (i, step) in cascade.steps.iter().enumerate() {
                    if i == 0 {
                        out.push_str(&render_terms(&step.reactants));
                        out.push(' ');
                    } else {
                        out.push(' ');
                    }
                    out.push_str(&step.render_arrow());
                    out.push(' ');
                    out.push_str(&render_terms(&step.products));
                }
                write!(f, "{out}")
            }
            Node::ElectronConfiguration(ec) => {
                let mut parts: Vec<String> = ec
                    .orbitals
                    .iter()
                    .map(|(orb, occ)| format!("{orb}^{occ}"))
                    .collect();
                // Term symbols stay unsupported in the Rust core until
                // they enter the corpus round-trip set.
                let _ = &mut parts;
                write!(f, "{}", parts.join(" "))
            }
            Node::Mechanism(mech) => {
                let mut out = String::from("mechanism");
                if !mech.steps.is_empty() || !mech.spectators.is_empty() {
                    let mut lines: Vec<String> = mech
                        .steps
                        .iter()
                        .map(|s| format!("{}: {}", s.label, s.text()))
                        .collect();
                    for sp in &mech.spectators {
                        lines.push(format!("spectator: {sp}"));
                    }
                    out.push_str(&format!("{{\n  {}\n}}", lines.join("\n  ")));
                }
                write!(f, "{out}")
            }
            Node::EmbeddedMath(source) => write!(f, "`{source}`"),
            Node::Text(content) => write!(f, "\"{content}\""),
        }
    }
}
