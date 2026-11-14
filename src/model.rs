//! Semantic model subset (formulas, molecules, reactions, mechanisms)
//! with the canonical text renderer. The renderer mirrors the
//! reference Text formatter rule for rule; round-trip is the
//! contract.

use std::collections::BTreeMap;
use std::fmt;

/// The reference's `String#to_i` for the model's numeric fields:
/// optional sign, leading digits, 0 when none (`"2"` → 2, `"-2"` →
/// -2, `"2a"` → 2, `"a"` → 0). Single parse for every consumer so a
/// count field never gets re-parsed per caller.
pub(crate) fn leading_integer(value: &str) -> i64 {
    let trimmed = value.trim_start();
    let (sign, start) = match trimmed.as_bytes().first() {
        Some(b'+') => (1i64, 1),
        Some(b'-') => (-1i64, 1),
        _ => (1i64, 0),
    };
    let digits: String = trimmed[start..]
        .chars()
        .take_while(|digit| digit.is_ascii_digit())
        .collect();
    if digits.is_empty() {
        return 0;
    }
    sign * digits.parse::<i64>().unwrap_or(0)
}

#[derive(Debug, Clone, Default, PartialEq)]
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
    /// Aromatic atom flag, set by SMILES lowercase input (`c1ccccc1`)
    /// or molfile aromatic contexts. Absent means not aromatic.
    pub aromatic: Option<bool>,
    /// Hydrogen count as stated by the source (SMILES bracket atoms
    /// `[CH4]`/`[nH]`). Absent means "per valence rules" — never
    /// invented.
    pub hydrogens: Option<u32>,
    pub chirality: Option<String>,
    /// Authored 2D/3D coordinates (molfile atom block). Model-only
    /// render data — the canonical wire form carries none; structure
    /// writers use them when present and lay out otherwise.
    pub x2: Option<f64>,
    pub y2: Option<f64>,
    pub z2: Option<f64>,
    /// Fractional (crystallographic) coordinates from the `@f(x,y,z)`
    /// atom annotation. Model-only render data — the canonical wire
    /// atom carries none; the crystal's atoms array references them
    /// through the same atom form.
    pub x_fract: Option<f64>,
    pub y_fract: Option<f64>,
    pub z_fract: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    pub kind: String,
    pub direction: Option<String>,
}

impl Bond {
    /// The plain-text bond glyph (the reference `Bond#ascii`) — the
    /// HTML renderer's bond output and the text canonical form.
    pub(crate) fn ascii(&self) -> &'static str {
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
    pub fn chars(&self) -> (&'static str, &'static str) {
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
    pub fn stereo_letter(&self) -> Option<&str> {
        self.stereo.as_deref()
    }

    /// The integer stoichiometric coefficient: the authored leading
    /// number, or 1 when absent or not positive. Single source for
    /// every consumer (reaction stoichiometry, the balance check,
    /// composition) — coefficients are model data, not per-caller
    /// parsing.
    pub fn stoichiometric_coefficient(&self) -> usize {
        match self.coefficient.as_deref().map(leading_integer) {
            Some(count) if count > 0 => count as usize,
            _ => 1,
        }
    }

    /// Element symbol => total atom count, recursing through groups
    /// and nested molecules with subscripts and multiplicities
    /// applied. Pass `with_coefficient: false` for the composition of
    /// a single molecule instance — identifiers describe the
    /// substance, not the stoichiometric multiplier.
    pub fn element_counts(&self, with_coefficient: bool) -> BTreeMap<String, usize> {
        let multiplier = if with_coefficient {
            self.stoichiometric_coefficient()
        } else {
            1
        };
        let mut tally = BTreeMap::new();
        for node in &self.nodes {
            tally_element(node, multiplier, &mut tally);
        }
        tally
    }

    /// Hill-system canonical formula: C first, then H, then others
    /// alphabetically. Single source of truth for "what is the
    /// canonical formula string of this molecule"
    /// (`H_2SO_4` → `"H2O4S"` — no carbon, alpha sort).
    pub fn hill_formula(&self) -> String {
        let counts = self.element_counts(true);
        let mut ordered: Vec<&String> = Vec::with_capacity(counts.len());
        // Hill order: C first, H second, the rest alphabetical (the
        // BTreeMap iteration).
        if let Some((key, _)) = counts.get_key_value("C") {
            ordered.push(key);
        }
        if let Some((key, _)) = counts.get_key_value("H") {
            ordered.push(key);
        }
        ordered.extend(
            counts
                .keys()
                .filter(|element| element.as_str() != "C" && element.as_str() != "H"),
        );
        ordered
            .iter()
            .map(|element| {
                let count = counts[*element];
                if count == 1 {
                    (*element).clone()
                } else {
                    format!("{element}{count}")
                }
            })
            .collect()
    }
}

fn tally_element(node: &Node, multiplier: usize, tally: &mut BTreeMap<String, usize>) {
    match node {
        Node::Atom(atom) => {
            let sub = atom
                .subscript
                .as_deref()
                .map(leading_integer)
                .filter(|&count| count > 0)
                .unwrap_or(1) as usize;
            *tally.entry(atom.element.clone()).or_insert(0) += sub * multiplier;
        }
        Node::Group(group) => {
            let inner = group
                .multiplicity
                .as_deref()
                .map(leading_integer)
                .filter(|&count| count > 0)
                .unwrap_or(1) as usize;
            for child in &group.nodes {
                tally_element(child, multiplier * inner, tally);
            }
        }
        Node::Molecule(molecule) => {
            let inner = molecule.stoichiometric_coefficient();
            for child in &molecule.nodes {
                tally_element(child, multiplier * inner, tally);
            }
        }
        _ => {}
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

    /// Human-readable label for linter diagnostics (cascade steps are
    /// bare `Reaction`s, so they carry the label too).
    pub fn diagnostic_label(&self) -> String {
        "Reaction".to_string()
    }

    /// The arrow's Unicode glyph (the reference `Reaction#arrow_entity`,
    /// the MathML/HTML entity form).
    pub(crate) fn arrow_entity(&self) -> &'static str {
        match self.arrow.as_str() {
            "reverse" => "\u{2190}",
            "equilibrium" => "\u{21CC}",
            "resonance" => "\u{2194}",
            _ => "\u{2192}",
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

/// A crystal structure: unit-cell parameters, space group, and the
/// asymmetric-unit atoms with fractional coordinates. The cell
/// parameters carry the authored strings (the reference's model keeps
/// the raw captures too); the wire form casts them to numbers.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Crystal {
    pub name: Option<String>,
    pub a: Option<String>,
    pub b: Option<String>,
    pub c: Option<String>,
    pub alpha: Option<String>,
    pub beta: Option<String>,
    pub gamma: Option<String>,
    pub spacegroup: Option<String>,
    pub atoms: Vec<Atom>,
}

impl Crystal {
    /// The six cell parameters in canonical iteration order, paired
    /// with their text labels — the single source of truth for "which
    /// cell fields are set and in what order" (the reference's
    /// `each_cell_param`); formatters consult this instead of inlining
    /// their own field lists.
    pub fn cell_params(&self) -> Vec<(&'static str, &str)> {
        [
            ("a", &self.a),
            ("b", &self.b),
            ("c", &self.c),
            ("alpha", &self.alpha),
            ("beta", &self.beta),
            ("gamma", &self.gamma),
        ]
        .into_iter()
        .filter_map(|(label, value)| value.as_deref().map(|value| (label, value)))
        .collect()
    }

    /// Per-format label overrides for the angle names: MathML, HTML
    /// and LaTeX spell the Greek letters (the reference's
    /// `CELL_LABELS`); text spells them out. `None` = the plain name.
    pub fn greek_label(label: &str, format: &str) -> String {
        match (label, format) {
            ("alpha", "mathml" | "html") => "α".to_string(),
            ("beta", "mathml" | "html") => "β".to_string(),
            ("gamma", "mathml" | "html") => "γ".to_string(),
            ("alpha", "latex") => "\\alpha".to_string(),
            ("beta", "latex") => "\\beta".to_string(),
            ("gamma", "latex") => "\\gamma".to_string(),
            _ => label.to_string(),
        }
    }

    /// Human-readable label for linter diagnostics ("Crystal(NaCl)";
    /// "unnamed" when no name is authored).
    pub fn diagnostic_label(&self) -> String {
        format!("Crystal({})", self.name.as_deref().unwrap_or("unnamed"))
    }
}

/// One spectrum peak: position plus intensity, splitting pattern, and
/// the structural assignment (the reference's `SpectrumPeak`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SpectrumPeak {
    pub position: Option<String>,
    pub intensity: Option<String>,
    pub multiplicity: Option<String>,
    pub assignment: Option<String>,
}

/// A spectroscopy record: technique, acquisition parameters (authored
/// order preserved — the wire's map iterates the same way), and peaks.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Spectrum {
    pub technique: Option<String>,
    pub params: Vec<(String, String)>,
    pub peaks: Vec<SpectrumPeak>,
}

impl Spectrum {
    /// Human-readable label for linter diagnostics
    /// ("Spectrum(1H)"; "unknown" when no technique is authored).
    pub fn diagnostic_label(&self) -> String {
        format!(
            "Spectrum({})",
            self.technique.as_deref().unwrap_or("unknown")
        )
    }
}

/// One computed property line of a calculation record.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CalculationProperty {
    pub title: String,
    pub value: String,
    pub units: Option<String>,
}

/// A computational-chemistry result: method, basis set, properties.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Calculation {
    pub method: Option<String>,
    pub basis: Option<String>,
    pub properties: Vec<CalculationProperty>,
}

impl Calculation {
    /// Human-readable label for linter diagnostics
    /// ("Calculation(B3LYP/6-31G*)").
    pub fn diagnostic_label(&self) -> String {
        format!(
            "Calculation({}/{})",
            self.method.as_deref().unwrap_or_default(),
            self.basis.as_deref().unwrap_or_default()
        )
    }
}

/// One z-matrix row: an atom placed by distance/angle/dihedral to
/// previously defined atoms (the reference's `ZRow`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ZRow {
    pub atom: String,
    pub ref1: Option<String>,
    pub distance: Option<String>,
    pub ref2: Option<String>,
    pub angle: Option<String>,
    pub ref3: Option<String>,
    pub dihedral: Option<String>,
}

/// Internal coordinates: the ordered z-matrix rows.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ZMatrix {
    pub rows: Vec<ZRow>,
}

impl ZMatrix {
    /// Human-readable label for linter diagnostics ("ZMatrix(3 rows)").
    pub fn diagnostic_label(&self) -> String {
        format!("ZMatrix({} rows)", self.rows.len())
    }
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
    Crystal(Crystal),
    Spectrum(Spectrum),
    Calculation(Calculation),
    ZMatrix(ZMatrix),
    EmbeddedMath(String),
    Text(String),
}

impl Molecule {
    /// Human-readable label for linter diagnostics — the Hill formula
    /// when present, "Molecule" otherwise.
    pub fn diagnostic_label(&self) -> String {
        let formula = self.hill_formula();
        if formula.is_empty() {
            "Molecule".to_string()
        } else {
            formula
        }
    }
}

impl Atom {
    /// Human-readable label for linter diagnostics ("Atom(C)").
    pub fn diagnostic_label(&self) -> String {
        format!("Atom({})", self.element)
    }
}

impl Group {
    /// Human-readable label for linter diagnostics ("Group").
    pub fn diagnostic_label(&self) -> String {
        "Group".to_string()
    }
}

impl Node {
    /// Human-readable label for linter diagnostics — the Rust mirror
    /// of `Model::Node#diagnostic_label` (the default snake-case name
    /// with words separated; `Atom` identifies itself by element,
    /// `Mechanism` by step count).
    pub fn diagnostic_label(&self) -> String {
        match self {
            Node::Atom(atom) => format!("Atom({})", atom.element),
            Node::Bond(_) => "Bond".to_string(),
            Node::Group(_) => "Group".to_string(),
            Node::Molecule(_) => "Molecule".to_string(),
            Node::Reaction(reaction) => reaction.diagnostic_label(),
            Node::ReactionCascade(_) => "Reaction Cascade".to_string(),
            Node::ElectronConfiguration(_) => "Electron Configuration".to_string(),
            Node::Mechanism(mechanism) => {
                format!("Mechanism({} steps)", mechanism.steps.len())
            }
            Node::Crystal(crystal) => crystal.diagnostic_label(),
            Node::Spectrum(spectrum) => spectrum.diagnostic_label(),
            Node::Calculation(calculation) => calculation.diagnostic_label(),
            Node::ZMatrix(zmatrix) => zmatrix.diagnostic_label(),
            Node::EmbeddedMath(_) => "Embedded Math".to_string(),
            Node::Text(_) => "Text".to_string(),
        }
    }
}

fn render_terms(nodes: &[Node]) -> String {
    nodes
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(" + ")
}

/// The canonical text form of one atom — the reference Text
/// formatter's `visit_atom` plus its annotation tail. Shared by the
/// Node Display and the crystal body (whose atoms render through the
/// same text rule).
fn atom_text(atom: &Atom) -> String {
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
    if let (Some(x), Some(y), Some(z)) = (atom.x_fract, atom.y_fract, atom.z_fract) {
        parts.push(format!(
            "@f({},{},{})",
            format_coord(x),
            format_coord(y),
            format_coord(z)
        ));
    }
    parts.concat()
}

/// Coordinate text: integral values drop the decimal (the reference's
/// `format_coord` — 0.0 renders "0", 0.5 renders "0.5").
fn format_coord(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{}", value.trunc() as i64)
    } else {
        format!("{value}")
    }
}

/// Float text in the reference's `Float#to_s` shape: integral values
/// keep one decimal (2.0 renders "2.0", 0.5 renders "0.5").
pub(crate) fn format_float(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
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
            Node::Atom(atom) => write!(f, "{}", atom_text(atom)),
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
            Node::Crystal(crystal) => {
                let mut out = String::from("crystal");
                if let Some(name) = &crystal.name {
                    out.push_str(&format!("[{name}]"));
                }
                let mut params: Vec<String> = crystal
                    .cell_params()
                    .into_iter()
                    .map(|(label, value)| format!("{label}={value}"))
                    .collect();
                if let Some(sg) = &crystal.spacegroup {
                    params.push(format!("sg={sg}"));
                }
                if !params.is_empty() {
                    out.push_str(&format!("({})", params.join(",")));
                }
                if !crystal.atoms.is_empty() {
                    let atoms: Vec<String> = crystal.atoms.iter().map(atom_text).collect();
                    out.push_str(&format!("{{{}}}", atoms.join(" ")));
                }
                write!(f, "{out}")
            }
            Node::Spectrum(spectrum) => {
                let mut out = String::from("spectrum");
                if let Some(technique) = &spectrum.technique {
                    out.push_str(&format!("[{technique}]"));
                }
                let params = spectrum
                    .params
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect::<Vec<_>>()
                    .join(",");
                if !params.is_empty() {
                    out.push_str(&format!("({params})"));
                }
                if !spectrum.peaks.is_empty() {
                    let lines = spectrum
                        .peaks
                        .iter()
                        .map(|peak| {
                            let mut line = format!(
                                "{}: {}",
                                peak.position.as_deref().unwrap_or_default(),
                                peak.intensity.as_deref().unwrap_or_default()
                            );
                            if let Some(multiplicity) = &peak.multiplicity {
                                line.push(' ');
                                line.push_str(multiplicity);
                            }
                            if let Some(assignment) = &peak.assignment {
                                line.push_str(&format!(" \"{assignment}\""));
                            }
                            line
                        })
                        .collect::<Vec<_>>()
                        .join("\n  ");
                    out.push_str(&format!("{{\n  {lines}\n}}"));
                }
                write!(f, "{out}")
            }
            Node::Calculation(calc) => {
                let mut out = String::from("calc");
                let mut params: Vec<String> = Vec::new();
                if let Some(method) = &calc.method {
                    params.push(method.clone());
                }
                if let Some(basis) = &calc.basis {
                    params.push(basis.clone());
                }
                if !params.is_empty() {
                    out.push_str(&format!("({})", params.join("/")));
                }
                if !calc.properties.is_empty() {
                    let lines = calc
                        .properties
                        .iter()
                        .map(|property| {
                            let mut line = format!("{}: {}", property.title, property.value);
                            if let Some(units) = &property.units {
                                line.push(' ');
                                line.push_str(units);
                            }
                            line
                        })
                        .collect::<Vec<_>>()
                        .join("\n  ");
                    out.push_str(&format!("{{\n  {lines}\n}}"));
                }
                write!(f, "{out}")
            }
            Node::ZMatrix(zmatrix) => {
                let mut out = String::from("zmatrix");
                if !zmatrix.rows.is_empty() {
                    let lines = zmatrix
                        .rows
                        .iter()
                        .map(|row| {
                            // The reference's token list: presence of a
                            // ref gates its pair, `compact` drops any
                            // nil tail before the two-space join.
                            let mut tokens: Vec<Option<&str>> = vec![Some(&row.atom)];
                            if let Some(ref1) = &row.ref1 {
                                tokens.push(Some(ref1));
                                tokens.push(row.distance.as_deref());
                            }
                            if let Some(ref2) = &row.ref2 {
                                tokens.push(Some(ref2));
                                tokens.push(row.angle.as_deref());
                            }
                            if let Some(ref3) = &row.ref3 {
                                tokens.push(Some(ref3));
                                tokens.push(row.dihedral.as_deref());
                            }
                            tokens.into_iter().flatten().collect::<Vec<_>>().join("  ")
                        })
                        .collect::<Vec<_>>()
                        .join("\n  ");
                    out.push_str(&format!("{{\n  {lines}\n}}"));
                }
                write!(f, "{out}")
            }
            Node::EmbeddedMath(source) => write!(f, "`{source}`"),
            Node::Text(content) => write!(f, "\"{content}\""),
        }
    }
}
