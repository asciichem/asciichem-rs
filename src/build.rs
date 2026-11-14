//! Value tree → semantic model. The Rust counterpart of the
//! reference's TransformRules, scoped to the formula/reaction subset;
//! constructs outside the subset raise `Unsupported` rather than
//! guessing.

use std::collections::HashMap;

use parsanol::portable::transform::Value;

use crate::model::{
    Atom, Bond, Calculation, CalculationProperty, Conditions, Crystal, ElectronConfiguration, Flow,
    Group, Identifier, Mechanism, MechanismStep, Molecule, Name, Node, Reaction, ReactionCascade,
    Spectrum, SpectrumPeak, ZMatrix, ZRow,
};

#[derive(Debug, Clone, PartialEq)]
pub enum BuildError {
    Unsupported(String),
    Malformed(String),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::Unsupported(what) => {
                write!(f, "not yet supported in the Rust core: {what}")
            }
            BuildError::Malformed(what) => write!(f, "malformed capture tree: {what}"),
        }
    }
}

type Result<T> = std::result::Result<T, BuildError>;

/// Atom annotation captures with no model home yet (cartesian
/// coordinates, parity, spin, titles). The fractional triple is NOT
/// here: `@f(x,y,z)` maps into the atom (the crystal node's data).
const UNSUPPORTED_ATOM_KEYS: [&str; 6] = [
    "x2",
    "y2",
    "z2",
    "atom_parity",
    "spin_multiplicity",
    "atom_title",
];

const BOND_KEYS: [(&str, &str); 10] = [
    ("direction_up", "single"),
    ("direction_down", "single"),
    ("quadruple", "quadruple"),
    ("triple", "triple"),
    ("wedge", "wedge"),
    ("hash", "hash"),
    ("dative", "dative"),
    ("wavy", "wavy"),
    ("double", "double"),
    ("single", "single"),
];

const ARROW_KINDS: [(&str, &str); 4] = [
    ("<=>", "equilibrium"),
    ("<->", "resonance"),
    ("->", "forward"),
    ("<-", "reverse"),
];

const STEREO_LETTERS: [(&str, &str); 10] = [
    ("R", "R"),
    ("S", "S"),
    ("E", "E"),
    ("Z", "Z"),
    ("a", "alpha"),
    ("alpha", "alpha"),
    ("α", "alpha"),
    ("b", "beta"),
    ("beta", "beta"),
    ("β", "beta"),
];

const IDENTIFIER_TYPES: [&str; 4] = ["inchi", "smiles", "cas", "chebi"];

pub fn build_formula(value: &Value) -> Result<Vec<Node>> {
    let hash = expect_hash(value)?;
    let subtree = hash
        .get("formula")
        .ok_or_else(|| BuildError::Malformed("tree has no formula capture".into()))?;
    match subtree {
        Value::Array(items) => items.iter().map(build_node).collect(),
        Value::Nil => Ok(vec![]),
        other => Ok(vec![build_node(other)?]),
    }
}

fn expect_hash(value: &Value) -> Result<&HashMap<String, Value>> {
    match value {
        Value::Hash(hash) => Ok(hash),
        _ => Err(BuildError::Malformed(format!(
            "expected hash, got {value:?}"
        ))),
    }
}

fn str_of(value: Option<&Value>) -> Option<String> {
    match value {
        Some(Value::String(s)) => Some(s.clone()),
        _ => None,
    }
}

fn strip_value(value: Option<&Value>) -> Option<String> {
    let s = str_of(value)?;
    let s = s.trim().to_string();
    if s.is_empty() {
        None
    } else {
        Some(s)
    }
}

pub fn build_node(value: &Value) -> Result<Node> {
    let hash = expect_hash(value)?;
    if let Some(text) = hash.get("text_run").or_else(|| hash.get("group_text_run")) {
        return Ok(Node::Text(strip_quotes(
            &str_of(Some(text)).unwrap_or_default(),
        )));
    }
    if let Some(source) = hash.get("math_source") {
        return Ok(Node::EmbeddedMath(str_of(Some(source)).unwrap_or_default()));
    }
    if let Some(atom) = hash.get("atom") {
        return build_atom(atom);
    }
    if let Some(ec) = hash.get("electron_config") {
        return build_electron_configuration(ec);
    }
    if let Some(mech) = hash.get("mechanism_node") {
        return build_mechanism(mech);
    }
    if let Some(crystal) = hash.get("crystal_node") {
        return build_crystal(crystal);
    }
    if let Some(spectrum) = hash.get("spectrum_node") {
        return build_spectrum(spectrum);
    }
    if let Some(calc) = hash.get("calc_node") {
        return build_calculation(calc);
    }
    if let Some(zmatrix) = hash.get("zmatrix_node") {
        return build_zmatrix(zmatrix);
    }
    if hash.contains_key("open_bracket") && hash.contains_key("group_nodes") {
        return build_group(hash);
    }
    if hash.contains_key("mol") && hash.contains_key("annotations") {
        return build_annotated_molecule(hash);
    }
    if hash.contains_key("reactants") && hash.contains_key("arrow") && hash.contains_key("products")
    {
        return build_reaction_from_hash(hash);
    }
    if let Some(cascade) = hash.get("cascade") {
        return build_cascade(cascade);
    }
    if hash.contains_key("units") {
        return build_molecule(hash);
    }
    for (key, kind) in BOND_KEYS {
        if hash.contains_key(key) {
            let direction = match key {
                "direction_up" => Some("up".to_string()),
                "direction_down" => Some("down".to_string()),
                _ => None,
            };
            return Ok(Node::Bond(Bond {
                kind: kind.to_string(),
                direction,
            }));
        }
    }
    Err(BuildError::Unsupported(format!("capture keys {:?}", {
        let mut keys: Vec<String> = hash.keys().cloned().collect();
        keys.sort();
        keys
    })))
}

fn strip_quotes(s: &str) -> String {
    if s.starts_with('"') && s.ends_with('"') && s.len() >= 2 {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

fn build_atom(value: &Value) -> Result<Node> {
    let hash = expect_hash(value)?;
    for key in UNSUPPORTED_ATOM_KEYS {
        if hash.contains_key(key) {
            return Err(BuildError::Unsupported(format!(
                "crystal-coordinate atom fields ({key})"
            )));
        }
    }
    let element = str_of(hash.get("element"))
        .ok_or_else(|| BuildError::Malformed("atom without element".into()))?;
    let superscript = strip_marker(str_of(hash.get("superscript")).as_deref(), Some('^'));
    let charge = detected_charge(superscript.as_deref());
    let oxidation_state = detected_oxidation_state(superscript.as_deref());
    let superscript = if charge.is_none() && oxidation_state.is_none() {
        superscript
    } else {
        None
    };
    Ok(Node::Atom(Atom {
        element,
        isotope: strip_marker(str_of(hash.get("isotope")).as_deref(), None),
        subscript: strip_marker(str_of(hash.get("subscript")).as_deref(), Some('_')),
        superscript,
        charge,
        oxidation_state,
        lone_pairs: lewis_count(hash.get("lone_pairs")),
        radical_electrons: lewis_count(hash.get("radical_electrons")),
        ring_closures: strip_value(hash.get("ring_closures")),
        // Aromatic flags and stated hydrogens are structure-interchange
        // fields (SMILES lowercase atoms, bracket `[nH]`); the AsciiChem
        // text grammar does not capture them.
        aromatic: None,
        hydrogens: None,
        chirality: str_of(hash.get("chirality")),
        x2: None,
        y2: None,
        z2: None,
        x_fract: fract_float(hash.get("x_fract")),
        y_fract: fract_float(hash.get("y_fract")),
        z_fract: fract_float(hash.get("z_fract")),
    }))
}

// Ruby's `String#to_f` for the fractional captures: the leading
// numeric prefix, 0.0 when there is none (the grammar only produces
// clean floats, so the fallback mirrors the reference without
// second-guessing it).
fn fract_float(value: Option<&Value>) -> Option<f64> {
    let s = str_of(value)?;
    let trimmed = s.trim_start();
    let digits: String = trimmed
        .trim_start_matches(['+', '-'])
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '.')
        .collect();
    let parsed: f64 = digits.parse().unwrap_or(0.0);
    Some(match trimmed.starts_with('-') {
        true => -parsed,
        false => parsed,
    })
}

// Strip a leading marker char (`^` or `_`) from the captured value;
// nil/empty → None. Mirrors AtomBuilder.strip_marker.
fn strip_marker(value: Option<&str>, marker: Option<char>) -> Option<String> {
    let mut s = value?.to_string();
    if let Some(marker) = marker {
        if s.starts_with(marker) {
            s.remove(0);
        }
    }
    if s.starts_with('^') || s.starts_with('_') {
        s.remove(0);
    }
    if s.is_empty() {
        None
    } else {
        Some(s)
    }
}

fn detected_charge(superscript: Option<&str>) -> Option<String> {
    let s = strip_marker(superscript, Some('^'))?;
    let s = s.as_str();
    // \A(?<n>\d*)(?<sign>[+-])\z  or  \A(?<sign>[+-])(?<n>\d*)\z
    let bytes = s.as_bytes();
    let (n, sign) = match (bytes.first(), bytes.last()) {
        (Some(first), Some(last)) if (*last == b'+' || *last == b'-') && first.is_ascii_digit() => {
            (s.split_at(s.len() - 1).0, *last)
        }
        (Some(first), Some(last))
            if (*first == b'+' || *first == b'-') && last.is_ascii_digit() =>
        {
            (&s[1..], *first)
        }
        (Some(only), _) if *only == b'+' || *only == b'-' => ("", *only),
        _ => return None,
    };
    let all_digits = n.bytes().all(|b| b.is_ascii_digit());
    if !all_digits {
        return None;
    }
    Some(if n.is_empty() {
        (sign as char).to_string()
    } else {
        format!("{}{}", n, sign as char)
    })
}

fn detected_oxidation_state(superscript: Option<&str>) -> Option<String> {
    let s = strip_marker(superscript, Some('^'))?;
    let inner = s.strip_prefix('(')?.strip_suffix(')')?;
    if !inner.is_empty() && inner.bytes().all(|b| b"IVXLCDM".contains(&b)) {
        Some(inner.to_string())
    } else {
        None
    }
}

fn lewis_count(value: Option<&Value>) -> Option<u32> {
    let s = str_of(value)?;
    let len = s.len() as u32;
    if len > 0 {
        Some(len)
    } else {
        None
    }
}

fn build_molecule(hash: &HashMap<String, Value>) -> Result<Node> {
    let nodes = match hash.get("units") {
        Some(Value::Array(items)) => items.iter().map(build_node).collect::<Result<Vec<_>>>()?,
        Some(Value::Nil) | None => vec![],
        Some(other) => vec![build_node(other)?],
    };
    let coefficient = match hash.get("coefficient") {
        Some(Value::Hash(coef)) => str_of(coef.get("value")),
        Some(Value::String(s)) => Some(s.clone()),
        _ => None,
    };
    let stereo = match str_of(hash.get("stereo")).as_deref() {
        None => None,
        Some(letter) => {
            let symbol = STEREO_LETTERS
                .iter()
                .find(|(l, _)| *l == letter)
                .map(|(_, symbol)| symbol.to_string())
                .ok_or_else(|| BuildError::Unsupported(format!("stereo letter {letter:?}")))?;
            Some(symbol)
        }
    };
    Ok(Node::Molecule(Molecule {
        nodes,
        coefficient,
        stereo,
        names: vec![],
        identifiers: vec![],
        title: None,
        formulas: vec![],
        labels: vec![],
        properties: vec![],
        metadata: vec![],
    }))
}

fn build_annotated_molecule(hash: &HashMap<String, Value>) -> Result<Node> {
    let mut node = build_node(
        hash.get("mol")
            .ok_or_else(|| BuildError::Malformed("annotated molecule without mol".into()))?,
    )?;
    let annotations = match hash.get("annotations") {
        Some(Value::Array(items)) => items.clone(),
        Some(Value::Nil) | None => vec![],
        Some(other) => vec![other.clone()],
    };
    if let Node::Molecule(mol) = &mut node {
        for ann in &annotations {
            apply_annotation(mol, expect_hash(ann)?);
        }
        return Ok(node);
    }
    Err(BuildError::Malformed(
        "annotated molecule body is not a molecule".into(),
    ))
}

fn apply_annotation(mol: &mut Molecule, ann: &HashMap<String, Value>) {
    if let Some(meta_key) = str_of(ann.get("meta_key")) {
        mol.metadata
            .push((meta_key, str_of(ann.get("meta_value")).unwrap_or_default()));
        return;
    }
    let ann_type = str_of(ann.get("ann_type")).unwrap_or_default();
    let value = str_of(ann.get("ann_value")).unwrap_or_default();
    match ann_type.as_str() {
        "name" => mol.names.push(Name { content: value }),
        "title" => mol.title = Some(value),
        "formula" => mol.formulas.push(value),
        "label" => mol.labels.push(value),
        kind if IDENTIFIER_TYPES.contains(&kind) => mol.identifiers.push(Identifier {
            value,
            convention: kind.to_string(),
        }),
        kind => mol.properties.push((kind.to_string(), value)),
    }
}

fn build_group(hash: &HashMap<String, Value>) -> Result<Node> {
    let nodes = match hash.get("group_nodes") {
        Some(Value::Array(items)) => items.iter().map(build_node).collect::<Result<Vec<_>>>()?,
        Some(Value::Nil) | None => vec![],
        Some(other) => vec![build_node(other)?],
    };
    let open = str_of(hash.get("open_bracket")).unwrap_or_default();
    let bracket = match open.as_str() {
        "[" => "square",
        "{" => "brace",
        _ => "paren",
    };
    let multiplicity =
        str_of(hash.get("multiplicity")).map(|m| m.strip_prefix('_').unwrap_or(&m).to_string());
    Ok(Node::Group(Group {
        nodes,
        multiplicity,
        bracket: bracket.to_string(),
    }))
}

fn build_reaction_from_hash(hash: &HashMap<String, Value>) -> Result<Node> {
    let reactants = as_terms(hash.get("reactants"))?;
    let products = as_terms(hash.get("products"))?;
    let arrow = expect_hash(
        hash.get("arrow")
            .ok_or_else(|| BuildError::Malformed("reaction without arrow".into()))?,
    )?;
    Ok(Node::Reaction(build_reaction(reactants, arrow, products)?))
}

fn build_reaction(
    reactants: Vec<Node>,
    arrow: &HashMap<String, Value>,
    products: Vec<Node>,
) -> Result<Reaction> {
    let kind = str_of(arrow.get("kind")).unwrap_or_default();
    let arrow_kind = ARROW_KINDS
        .iter()
        .find(|(literal, _)| *literal == kind)
        .map(|(_, name)| name.to_string())
        .unwrap_or_else(|| "forward".to_string());
    let above = condition_text(arrow.get("above"));
    let below = condition_text(arrow.get("below"));
    let conditions = if above.is_none() && below.is_none() {
        None
    } else {
        Some(Conditions { above, below })
    };
    Ok(Reaction {
        reactants,
        products,
        arrow: arrow_kind,
        conditions,
    })
}

fn condition_text(value: Option<&Value>) -> Option<String> {
    match value? {
        Value::Hash(hash) => str_of(hash.get("text")),
        Value::String(s) => Some(s.clone()),
        _ => None,
    }
}

fn as_terms(value: Option<&Value>) -> Result<Vec<Node>> {
    match value {
        Some(Value::Array(items)) => items.iter().map(build_node).collect(),
        Some(Value::Nil) | None => Ok(vec![]),
        Some(other) => Ok(vec![build_node(other)?]),
    }
}

fn build_cascade(value: &Value) -> Result<Node> {
    let hash = match value {
        Value::Array(items) => items
            .iter()
            .find_map(|item| {
                let hash = expect_hash(item).ok()?;
                if hash.contains_key("first") {
                    Some(hash)
                } else {
                    None
                }
            })
            .ok_or_else(|| BuildError::Malformed("cascade array without first".into()))?,
        Value::Hash(_) => expect_hash(value)?,
        _ => return Err(BuildError::Malformed("cascade is not a hash/array".into())),
    };
    let first = build_node(
        hash.get("first")
            .ok_or_else(|| BuildError::Malformed("cascade without first".into()))?,
    )?;
    let Node::Reaction(reaction) = first else {
        return Err(BuildError::Malformed(
            "cascade first is not a reaction".into(),
        ));
    };
    let mut steps = vec![reaction];
    let segments = match hash.get("segments") {
        Some(Value::Array(items)) => items.clone(),
        Some(Value::Nil) | None => vec![],
        Some(other) => vec![other.clone()],
    };
    for segment in &segments {
        let seg = expect_hash(segment)?;
        let arrow = expect_hash(
            seg.get("arrow")
                .ok_or_else(|| BuildError::Malformed("cascade segment without arrow".into()))?,
        )?;
        let products = as_terms(seg.get("products"))?;
        let prev_products = steps.last().expect("seeded").products.clone();
        steps.push(build_reaction(prev_products, arrow, products)?);
    }
    Ok(Node::ReactionCascade(ReactionCascade { steps }))
}

fn build_electron_configuration(value: &Value) -> Result<Node> {
    let pairs = match value {
        Value::Array(items) => items.clone(),
        Value::Nil => vec![],
        other => vec![other.clone()],
    };
    let mut orbitals = vec![];
    for pair in &pairs {
        let hash = expect_hash(pair)?;
        let orbital = str_of(hash.get("orbital"))
            .ok_or_else(|| BuildError::Malformed("orbital without label".into()))?;
        let occupancy = str_of(hash.get("occupancy"))
            .ok_or_else(|| BuildError::Malformed("orbital without occupancy".into()))?;
        orbitals.push((orbital, occupancy));
    }
    Ok(Node::ElectronConfiguration(ElectronConfiguration {
        orbitals,
    }))
}

fn build_mechanism(value: &Value) -> Result<Node> {
    // The grammar names `mechanism >> body.maybe`: with no body the
    // capture is the bare literal string — the empty mechanism (the
    // reference's MechanismBuilder gets a nil body).
    let hash = match value {
        Value::Hash(hash) => hash,
        Value::String(_) => return Ok(Node::Mechanism(Mechanism::default())),
        other => {
            return Err(BuildError::Malformed(format!(
                "expected hash, got {other:?}"
            )))
        }
    };
    let body = strip_value(hash.get("mechanism_body"));
    let mut steps = vec![];
    let mut spectators = vec![];
    for (label, value) in parse_entries(body.as_deref())? {
        if label == "spectator" {
            for token in value.split_whitespace() {
                spectators.push(token.to_string());
            }
        } else {
            let (reaction, flows) = split_flows(&label, &value)?;
            steps.push(MechanismStep {
                label,
                reaction,
                flows,
            });
        }
    }
    Ok(Node::Mechanism(Mechanism { steps, spectators }))
}

// -- crystals -------------------------------------------------------------

/// Builds a Crystal from the grammar's raw captures (name, the
/// `a=...,sg=...` parameter string, the atoms body) — the reference's
/// CrystalBuilder: the parameters split on `,`/`=`, the body re-parsed
/// through the full formula grammar with its atoms kept.
fn build_crystal(value: &Value) -> Result<Node> {
    let hash = expect_hash(value)?;
    let name = strip_value(hash.get("crystal_name"));
    let params = parse_crystal_params(str_of(hash.get("crystal_params")).as_deref());
    let cell = |key: &str| params.get(key).cloned().flatten();
    Ok(Node::Crystal(Crystal {
        name,
        a: cell("a"),
        b: cell("b"),
        c: cell("c"),
        alpha: cell("alpha"),
        beta: cell("beta"),
        gamma: cell("gamma"),
        spacegroup: cell("sg"),
        atoms: crystal_atoms(str_of(hash.get("crystal_body")).as_deref())?,
    }))
}

/// `a=5.64,sg=Fm-3m` → the parameter map. A pair without `=` carries
/// a nil value (the reference's `val&.strip`).
fn parse_crystal_params(params: Option<&str>) -> HashMap<String, Option<String>> {
    let Some(params) = params else {
        return HashMap::new();
    };
    params
        .split(',')
        .filter_map(|pair| {
            let (key, value) = pair.trim().split_once('=')?;
            let key = key.trim();
            (!key.is_empty()).then(|| (key.to_string(), Some(value.trim().to_string())))
        })
        .collect()
}

/// The crystal body's atoms: the body re-parses through the full
/// formula grammar and its molecules' atoms are kept (the reference's
/// `parse_atoms`). Runs on the CURRENT stack — the deep-stack worker
/// is single-threaded, so routing the inner parse through
/// [`crate::parse_nodes`] from a worker-driven build would deadlock.
fn crystal_atoms(body: Option<&str>) -> Result<Vec<Atom>> {
    let Some(body) = body else {
        return Ok(vec![]);
    };
    let nodes = crate::parse_nodes_direct(body)
        .map_err(|error| BuildError::Malformed(format!("crystal body failed to parse: {error}")))?;
    let mut atoms = Vec::new();
    for node in &nodes {
        let Node::Molecule(molecule) = node else {
            continue;
        };
        for inner in &molecule.nodes {
            if let Node::Atom(atom) = inner {
                atoms.push(atom.clone());
            }
        }
    }
    Ok(atoms)
}

// -- spectra, calculations, z-matrices ------------------------------------

/// Builds a Spectrum from the grammar's raw captures (technique, the
/// `k=v,...` parameter string, the peak-line body) — the reference's
/// SpectrumBuilder: each body line is one peak, a quoted tail is the
/// assignment. A bare `spectrum` literal (no body capture) is the
/// empty record.
fn build_spectrum(value: &Value) -> Result<Node> {
    let (technique, params_raw, body) = match value {
        Value::Hash(hash) => {
            let technique = strip_value(hash.get("spectrum_type"));
            let params_raw = str_of(hash.get("spectrum_params"));
            let body = strip_value(hash.get("spectrum_body"));
            (technique, params_raw, body)
        }
        // The grammar names `spectrum >> body.maybe`: with no body the
        // capture is the bare literal string — the empty record (the
        // reference's builder gets a nil body).
        Value::String(_) => (None, None, None),
        other => {
            return Err(BuildError::Malformed(format!(
                "expected hash, got {other:?}"
            )))
        }
    };
    validate_colon_lines(body.as_deref(), "spectrum peak")?;
    let mut peaks = Vec::new();
    for line in body.iter().flat_map(|body| body.split('\n')) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        // A quoted run anywhere in the line is the assignment (the
        // reference's /"([^"]*)"/ scan); it leaves the token stream.
        let mut assignment = None;
        let mut line = line.to_string();
        if let Some(open) = line.find('"') {
            if let Some(len) = line[open + 1..].find('"') {
                assignment = Some(line[open + 1..][..len].to_string());
                line.replace_range(open..=open + 1 + len, "");
                line = line.trim().to_string();
            }
        }
        let (position, rest) = match line.split_once(':') {
            Some((position, rest)) => (position.trim().to_string(), rest),
            None => (String::new(), line.as_str()),
        };
        let tokens: Vec<&str> = rest.split_whitespace().collect();
        peaks.push(SpectrumPeak {
            position: (!position.is_empty()).then_some(position),
            intensity: tokens.first().map(|token| (*token).to_string()),
            multiplicity: tokens.get(1).map(|token| (*token).to_string()),
            assignment,
        });
    }
    let params = parse_params(params_raw.as_deref());
    Ok(Node::Spectrum(Spectrum {
        technique,
        params,
        peaks,
    }))
}

/// Builds a Calculation from the `method/basis` parameter string and
/// the `title: value units` property lines (the reference's
/// CalculationBuilder). A bare `calc` literal is the empty record.
fn build_calculation(value: &Value) -> Result<Node> {
    let (params_raw, body) = match value {
        Value::Hash(hash) => {
            let params_raw = strip_value(hash.get("calc_params"));
            let body = strip_value(hash.get("calc_body"));
            (params_raw, body)
        }
        Value::String(_) => (None, None),
        other => {
            return Err(BuildError::Malformed(format!(
                "expected hash, got {other:?}"
            )))
        }
    };
    validate_colon_lines(body.as_deref(), "calculation property")?;
    let (method, basis) = match params_raw {
        None => (None, None),
        Some(params) => {
            let mut parts = params.split('/');
            let method = parts.next().unwrap_or_default().trim().to_string();
            let basis = parts.next().unwrap_or_default().trim().to_string();
            (
                (!method.is_empty()).then_some(method),
                (!basis.is_empty()).then_some(basis),
            )
        }
    };
    let mut properties = Vec::new();
    for line in body.iter().flat_map(|body| body.split('\n')) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Some((key, rest)) = line.split_once(':') else {
            continue;
        };
        if key.trim().is_empty() {
            continue;
        }
        let tokens: Vec<&str> = rest.split_whitespace().collect();
        properties.push(CalculationProperty {
            title: key.trim().to_string(),
            value: tokens.first().copied().unwrap_or_default().to_string(),
            units: tokens.get(1).map(|token| (*token).to_string()),
        });
    }
    Ok(Node::Calculation(Calculation {
        method,
        basis,
        properties,
    }))
}

/// Builds a ZMatrix: each body line is one row, whitespace-split into
/// the atom/distance/angle/dihedral token slots (the reference's
/// ZMatrixBuilder). A bare `zmatrix` literal is the empty record.
fn build_zmatrix(value: &Value) -> Result<Node> {
    let body = match value {
        Value::Hash(hash) => strip_value(hash.get("zmatrix_body")),
        Value::String(_) => None,
        other => {
            return Err(BuildError::Malformed(format!(
                "expected hash, got {other:?}"
            )))
        }
    };
    let mut rows = Vec::new();
    for line in body.iter().flat_map(|body| body.split('\n')) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }
        let token = |index: usize| tokens.get(index).map(|token| (*token).to_string());
        rows.push(ZRow {
            atom: tokens[0].to_string(),
            ref1: token(1),
            distance: token(2),
            ref2: token(3),
            angle: token(4),
            ref3: token(5),
            dihedral: token(6),
        });
    }
    Ok(Node::ZMatrix(ZMatrix { rows }))
}

/// `solvent=CDCl3,type=1H` → the ordered parameter map. A pair without
/// `=` keeps its key with the empty value (the reference's
/// `partition('=')` — unlike the crystal's params, spectra keep the
/// authored key).
fn parse_params(params: Option<&str>) -> Vec<(String, String)> {
    let Some(params) = params else {
        return Vec::new();
    };
    params
        .split(',')
        .filter_map(|pair| {
            let pair = pair.trim();
            // The reference's partition('='): no '=' keeps the whole
            // pair as the key with the empty value.
            let (key, value) = match pair.split_once('=') {
                Some((key, value)) => (key, value),
                None => (pair, ""),
            };
            let key = key.trim();
            (!key.is_empty()).then(|| (key.to_string(), value.trim().to_string()))
        })
        .collect()
}

/// The spectrum/calculation bodies are line-oriented `key: rest`
/// records; a nonblank line without the separator is malformed (the
/// reference's `validate_colon_lines`).
fn validate_colon_lines(body: Option<&str>, what: &str) -> Result<()> {
    let Some(body) = body else {
        return Ok(());
    };
    for (index, line) in body.split('\n').enumerate() {
        let line = line.trim();
        if !line.is_empty() && !line.contains(':') {
            return Err(BuildError::Malformed(format!(
                "{what} on line {} is missing ':' separator: {line:?}",
                index + 1
            )));
        }
    }
    Ok(())
}

fn parse_entries(body: Option<&str>) -> Result<Vec<(String, String)>> {
    let body = match body {
        None => return Ok(vec![]),
        Some(body) => body,
    };
    let mut entries = vec![];
    for line in body.split('\n') {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Some((key, value)) = line.split_once(':') else {
            return Err(BuildError::Malformed(format!(
                "mechanism entry is missing ':' separator: {line:?}"
            )));
        };
        entries.push((key.trim().to_string(), value.trim().to_string()));
    }
    Ok(entries)
}

fn split_flows(label: &str, value: &str) -> Result<(String, Vec<Flow>)> {
    let mut clauses = value.split('|');
    let reaction = clauses.next().unwrap_or_default().trim().to_string();
    let mut flows = vec![];
    for clause in clauses {
        let clause = clause.trim();
        let Some((kind, payload)) = clause.split_once(':') else {
            return Err(BuildError::Malformed(format!(
                "mechanism step {label:?} has an unknown clause: {clause:?} (supported: \"flow: SOURCE > TARGET\")"
            )));
        };
        if kind.trim() != "flow" || payload.trim().is_empty() {
            return Err(BuildError::Malformed(format!(
                "mechanism step {label:?} has an unknown clause: {clause:?} (supported: \"flow: SOURCE > TARGET\")"
            )));
        }
        for arrow in payload.split(';') {
            let arrow = arrow.trim();
            if arrow.is_empty() {
                continue;
            }
            let Some((source, target)) = arrow.split_once('>') else {
                return Err(BuildError::Malformed(format!(
                    "mechanism step {label:?} has a malformed flow: {arrow:?} (expected \"SOURCE > TARGET\")"
                )));
            };
            let (source, target) = (source.trim(), target.trim());
            if source.is_empty() || target.is_empty() {
                return Err(BuildError::Malformed(format!(
                    "mechanism step {label:?} has a malformed flow: {arrow:?} (expected \"SOURCE > TARGET\")"
                )));
            }
            flows.push(Flow {
                source: source.to_string(),
                target: target.to_string(),
            });
        }
    }
    Ok((reaction, flows))
}
