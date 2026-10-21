//! The canonical wire form — the Rust mirror of the shared
//! asciichem-model v1 schemas and the reference `WireAdapter`.
//! (De)serialisation goes through serde's derive (camelCase wire
//! keys, `type` tags, empty collections omitted); this module only
//! routes between the model subset and the wire types, field for
//! field, exactly as the reference adapter does.

use serde::{Deserialize, Serialize};

use crate::model::{self, Node};

// Each wire type carries its `type` discriminator as a real field
// (the Rust mirror of the Ruby wire classes' `wire_type`); a
// single-variant marker enum makes the value constant in both
// directions.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FormulaTag {
    #[serde(rename = "formula")]
    Formula,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MoleculeTag {
    #[serde(rename = "molecule")]
    Molecule,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GroupTag {
    #[serde(rename = "group")]
    Group,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AtomTag {
    #[serde(rename = "atom")]
    Atom,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BondTag {
    #[serde(rename = "bond")]
    Bond,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IdentifierTag {
    #[serde(rename = "identifier")]
    Identifier,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NameTag {
    #[serde(rename = "name")]
    Name,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReactionTag {
    #[serde(rename = "reaction")]
    Reaction,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReactionCascadeTag {
    #[serde(rename = "reaction-cascade")]
    ReactionCascade,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElectronConfigurationTag {
    #[serde(rename = "electron-configuration")]
    ElectronConfiguration,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EmbeddedMathTag {
    #[serde(rename = "embedded-math")]
    EmbeddedMath,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TextTag {
    #[serde(rename = "text")]
    Text,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MechanismTag {
    #[serde(rename = "mechanism")]
    Mechanism,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum WireNode {
    Formula(Formula),
    Molecule(Molecule),
    Group(Group),
    Atom(Atom),
    Bond(Bond),
    Identifier(Identifier),
    Name(Name),
    Reaction(Reaction),
    ReactionCascade(ReactionCascade),
    ElectronConfiguration(ElectronConfiguration),
    EmbeddedMath(EmbeddedMath),
    Text(Text),
    Mechanism(Mechanism),
}

impl WireNode {
    pub fn type_name(&self) -> &'static str {
        match self {
            WireNode::Formula(_) => "formula",
            WireNode::Molecule(_) => "molecule",
            WireNode::Group(_) => "group",
            WireNode::Atom(_) => "atom",
            WireNode::Bond(_) => "bond",
            WireNode::Identifier(_) => "identifier",
            WireNode::Name(_) => "name",
            WireNode::Reaction(_) => "reaction",
            WireNode::ReactionCascade(_) => "reaction-cascade",
            WireNode::ElectronConfiguration(_) => "electron-configuration",
            WireNode::EmbeddedMath(_) => "embedded-math",
            WireNode::Text(_) => "text",
            WireNode::Mechanism(_) => "mechanism",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Formula {
    #[serde(rename = "type")]
    pub tag: FormulaTag,
    pub nodes: Vec<WireNode>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    #[serde(rename = "type")]
    pub tag: MoleculeTag,
    pub nodes: Vec<WireNode>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coefficient: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stereo: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub identifiers: Vec<Identifier>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Group {
    #[serde(rename = "type")]
    pub tag: GroupTag,
    pub nodes: Vec<WireNode>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub multiplicity: Option<String>,
    pub bracket: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    #[serde(rename = "type")]
    pub tag: AtomTag,
    pub element: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub isotope: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub charge: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subscript: Option<String>,
    #[serde(rename = "oxidationState", skip_serializing_if = "Option::is_none")]
    pub oxidation_state: Option<String>,
    #[serde(rename = "lonePairs", skip_serializing_if = "Option::is_none")]
    pub lone_pairs: Option<u32>,
    #[serde(rename = "radicalElectrons", skip_serializing_if = "Option::is_none")]
    pub radical_electrons: Option<u32>,
    #[serde(rename = "ringClosures", skip_serializing_if = "Option::is_none")]
    pub ring_closures: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chirality: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    #[serde(rename = "type")]
    pub tag: BondTag,
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub direction: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Identifier {
    #[serde(rename = "type")]
    pub tag: IdentifierTag,
    pub value: String,
    pub convention: String,
    #[serde(rename = "dictRef", skip_serializing_if = "Option::is_none")]
    pub dict_ref: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Name {
    #[serde(rename = "type")]
    pub tag: NameTag,
    pub content: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub convention: Option<String>,
    #[serde(rename = "dictRef", skip_serializing_if = "Option::is_none")]
    pub dict_ref: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reaction {
    #[serde(rename = "type")]
    pub tag: ReactionTag,
    pub reactants: Vec<Molecule>,
    pub products: Vec<Molecule>,
    pub arrow: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conditions: Option<Conditions>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Conditions {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub above: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub below: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReactionCascade {
    #[serde(rename = "type")]
    pub tag: ReactionCascadeTag,
    pub steps: Vec<Reaction>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalOccupancy {
    pub orbital: String,
    pub occupancy: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronConfiguration {
    #[serde(rename = "type")]
    pub tag: ElectronConfigurationTag,
    pub orbitals: Vec<OrbitalOccupancy>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EmbeddedMath {
    #[serde(rename = "type")]
    pub tag: EmbeddedMathTag,
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Text {
    #[serde(rename = "type")]
    pub tag: TextTag,
    pub content: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MechanismFlow {
    pub source: String,
    pub target: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MechanismStep {
    pub label: String,
    pub reaction: Reaction,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub flows: Vec<MechanismFlow>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mechanism {
    #[serde(rename = "type")]
    pub tag: MechanismTag,
    pub steps: Vec<MechanismStep>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub spectators: Vec<Molecule>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WireError {
    Unsupported(String),
    Malformed(String),
}

impl std::fmt::Display for WireError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WireError::Unsupported(what) => write!(f, "no wire form: {what}"),
            WireError::Malformed(what) => write!(f, "malformed wire node: {what}"),
        }
    }
}

impl std::error::Error for WireError {}

type Result<T> = std::result::Result<T, WireError>;

// -- model -> wire --------------------------------------------------

pub fn formula_to_wire(nodes: &[Node]) -> Result<WireNode> {
    Ok(WireNode::Formula(Formula {
        tag: FormulaTag::Formula,
        nodes: nodes.iter().map(node_to_wire).collect::<Result<Vec<_>>>()?,
    }))
}

fn node_to_wire(node: &Node) -> Result<WireNode> {
    Ok(match node {
        Node::Molecule(mol) => WireNode::Molecule(molecule_to_wire(mol)),
        Node::Group(group) => WireNode::Group(Group {
            tag: GroupTag::Group,
            nodes: group
                .nodes
                .iter()
                .map(node_to_wire)
                .collect::<Result<Vec<_>>>()?,
            multiplicity: group.multiplicity.clone(),
            bracket: group.bracket.clone(),
        }),
        Node::Atom(atom) => WireNode::Atom(Atom {
            tag: AtomTag::Atom,
            element: atom.element.clone(),
            isotope: atom.isotope.clone(),
            charge: atom.charge.clone(),
            subscript: atom.subscript.clone(),
            oxidation_state: atom.oxidation_state.clone(),
            lone_pairs: atom.lone_pairs,
            radical_electrons: atom.radical_electrons,
            ring_closures: atom.ring_closures.clone(),
            chirality: atom.chirality.clone(),
        }),
        Node::Bond(bond) => WireNode::Bond(Bond {
            tag: BondTag::Bond,
            kind: bond.kind.clone(),
            direction: bond.direction.clone(),
        }),
        Node::Reaction(reaction) => WireNode::Reaction(reaction_to_wire(reaction)),
        Node::ReactionCascade(cascade) => WireNode::ReactionCascade(ReactionCascade {
            tag: ReactionCascadeTag::ReactionCascade,
            steps: cascade.steps.iter().map(reaction_to_wire).collect(),
        }),
        Node::ElectronConfiguration(ec) => WireNode::ElectronConfiguration(ElectronConfiguration {
            tag: ElectronConfigurationTag::ElectronConfiguration,
            orbitals: ec
                .orbitals
                .iter()
                .map(|(orbital, occupancy)| OrbitalOccupancy {
                    orbital: orbital.clone(),
                    occupancy: occupancy.clone(),
                })
                .collect(),
        }),
        Node::EmbeddedMath(source) => WireNode::EmbeddedMath(EmbeddedMath {
            tag: EmbeddedMathTag::EmbeddedMath,
            source: source.clone(),
        }),
        Node::Text(content) => WireNode::Text(Text {
            tag: TextTag::Text,
            content: content.clone(),
        }),
        Node::Mechanism(mech) => WireNode::Mechanism(mechanism_to_wire(mech)?),
    })
}

fn molecule_to_wire(mol: &model::Molecule) -> Molecule {
    Molecule {
        tag: MoleculeTag::Molecule,
        nodes: mol
            .nodes
            .iter()
            .map(|node| node_to_wire(node).expect("molecule children are wire-expressible"))
            .collect(),
        coefficient: mol.coefficient.clone(),
        stereo: mol.stereo.clone(),
        identifiers: mol
            .identifiers
            .iter()
            .map(|id| Identifier {
                tag: IdentifierTag::Identifier,
                value: id.value.clone(),
                convention: id.convention.clone(),
                dict_ref: None,
            })
            .collect(),
    }
}

fn reaction_to_wire(reaction: &model::Reaction) -> Reaction {
    fn side_terms(nodes: &[Node]) -> Vec<Molecule> {
        nodes
            .iter()
            .map(|node| match node {
                Node::Molecule(mol) => molecule_to_wire(mol),
                other => panic!("reaction term is not a molecule: {other:?}"),
            })
            .collect()
    }
    Reaction {
        tag: ReactionTag::Reaction,
        reactants: side_terms(&reaction.reactants),
        products: side_terms(&reaction.products),
        arrow: reaction.arrow.clone(),
        conditions: reaction.conditions.as_ref().map(|conds| Conditions {
            above: conds.above.clone(),
            below: conds.below.clone(),
        }),
    }
}

// A mechanism step's reaction is TEXT in the model but a typed node
// on the wire, so the same grammar parses it here — exactly like the
// reference adapter. A cascade step has no wire form and says so.
fn mechanism_to_wire(mech: &model::Mechanism) -> Result<Mechanism> {
    let mut steps = Vec::new();
    for step in &mech.steps {
        steps.push(MechanismStep {
            label: step.label.clone(),
            reaction: mechanism_reaction_to_wire(&step.reaction)?,
            flows: step
                .flows
                .iter()
                .map(|flow| MechanismFlow {
                    source: flow.source.clone(),
                    target: flow.target.clone(),
                })
                .collect(),
        });
    }
    let mut spectators = Vec::new();
    for text in &mech.spectators {
        spectators.push(mechanism_spectator_to_wire(text)?);
    }
    Ok(Mechanism {
        tag: MechanismTag::Mechanism,
        steps,
        spectators,
    })
}

// A mechanism step's reaction is TEXT in the model but a typed node
// on the wire, so the same grammar parses it here — exactly like the
// reference adapter. A cascade step has no wire form and says so.
fn mechanism_reaction_to_wire(text: &str) -> Result<Reaction> {
    let nodes = crate::parse_nodes(text).map_err(|error| {
        WireError::Unsupported(format!(
            "mechanism step is not a reaction: {text:?} ({error})"
        ))
    })?;
    match nodes.first() {
        Some(Node::Reaction(reaction)) => Ok(reaction_to_wire(reaction)),
        Some(Node::ReactionCascade(_)) => Err(WireError::Unsupported(format!(
            "mechanism step cascades have no wire form: {text:?}"
        ))),
        _ => Err(WireError::Unsupported(format!(
            "mechanism step is not a reaction: {text:?}"
        ))),
    }
}

fn mechanism_spectator_to_wire(text: &str) -> Result<Molecule> {
    let nodes = crate::parse_nodes(text).map_err(|error| {
        WireError::Unsupported(format!(
            "mechanism spectator is not valid AsciiChem: {text:?} ({error})"
        ))
    })?;
    match nodes.first() {
        Some(Node::Molecule(mol)) => Ok(molecule_to_wire(mol)),
        _ => Err(WireError::Unsupported(format!(
            "mechanism spectator is not a molecule: {text:?}"
        ))),
    }
}

// -- wire -> model ----------------------------------------------------

pub fn node_from_wire(wire: &WireNode) -> Result<Node> {
    Ok(match wire {
        // Formula is the wire root type only — the model has no
        // nested-formula node.
        WireNode::Formula(_) => {
            return Err(WireError::Malformed(
                "formula node outside the wire root".into(),
            ))
        }
        WireNode::Molecule(mol) => Node::Molecule(molecule_from_wire(mol)?),
        WireNode::Group(group) => Node::Group(model::Group {
            nodes: group
                .nodes
                .iter()
                .map(node_from_wire)
                .collect::<Result<Vec<_>>>()?,
            multiplicity: group.multiplicity.clone(),
            bracket: group.bracket.clone(),
        }),
        WireNode::Atom(atom) => Node::Atom(model::Atom {
            element: atom.element.clone(),
            isotope: atom.isotope.clone(),
            subscript: atom.subscript.clone(),
            superscript: None,
            charge: atom.charge.clone(),
            oxidation_state: atom.oxidation_state.clone(),
            lone_pairs: atom.lone_pairs,
            radical_electrons: atom.radical_electrons,
            ring_closures: atom.ring_closures.clone(),
            chirality: atom.chirality.clone(),
        }),
        WireNode::Bond(bond) => Node::Bond(model::Bond {
            kind: bond.kind.clone(),
            direction: bond.direction.clone(),
        }),
        // Identifiers ride on molecules in both the model and the
        // v1 schema; a standalone identifier node has no model form.
        WireNode::Identifier(_) | WireNode::Name(_) => {
            return Err(WireError::Unsupported(
                "standalone identifier/name nodes (they ride on molecules)".into(),
            ))
        }
        WireNode::Reaction(reaction) => Node::Reaction(reaction_from_wire(reaction)?),
        WireNode::ReactionCascade(cascade) => Node::ReactionCascade(model::ReactionCascade {
            steps: cascade
                .steps
                .iter()
                .map(reaction_from_wire)
                .collect::<Result<Vec<_>>>()?,
        }),
        WireNode::ElectronConfiguration(ec) => {
            Node::ElectronConfiguration(model::ElectronConfiguration {
                orbitals: ec
                    .orbitals
                    .iter()
                    .map(|orb| (orb.orbital.clone(), orb.occupancy.clone()))
                    .collect(),
            })
        }
        WireNode::EmbeddedMath(math) => Node::EmbeddedMath(math.source.clone()),
        WireNode::Text(text) => Node::Text(text.content.clone()),
        WireNode::Mechanism(mech) => Node::Mechanism(mechanism_from_wire(mech)?),
    })
}

fn molecule_from_wire(mol: &Molecule) -> Result<model::Molecule> {
    Ok(model::Molecule {
        nodes: mol
            .nodes
            .iter()
            .map(node_from_wire)
            .collect::<Result<Vec<_>>>()?,
        coefficient: mol.coefficient.clone(),
        stereo: mol.stereo.clone(),
        names: vec![],
        identifiers: mol
            .identifiers
            .iter()
            .map(|id| model::Identifier {
                value: id.value.clone(),
                convention: id.convention.clone(),
            })
            .collect(),
        title: None,
        formulas: vec![],
        labels: vec![],
        properties: vec![],
        metadata: vec![],
    })
}

fn reaction_from_wire(reaction: &Reaction) -> Result<model::Reaction> {
    Ok(model::Reaction {
        reactants: reaction
            .reactants
            .iter()
            .map(|mol| Ok(Node::Molecule(molecule_from_wire(mol)?)))
            .collect::<Result<Vec<_>>>()?,
        products: reaction
            .products
            .iter()
            .map(|mol| Ok(Node::Molecule(molecule_from_wire(mol)?)))
            .collect::<Result<Vec<_>>>()?,
        arrow: reaction.arrow.clone(),
        conditions: reaction.conditions.as_ref().map(|conds| model::Conditions {
            above: conds.above.clone(),
            below: conds.below.clone(),
        }),
    })
}

fn mechanism_from_wire(mech: &Mechanism) -> Result<model::Mechanism> {
    Ok(model::Mechanism {
        steps: mech
            .steps
            .iter()
            .map(|step| {
                Ok(model::MechanismStep {
                    label: step.label.clone(),
                    // The model keeps a step's reaction as TEXT; the
                    // wire carries the typed node, so ingestion
                    // renders it back through the canonical form.
                    reaction: Node::Reaction(reaction_from_wire(&step.reaction)?).to_string(),
                    flows: step
                        .flows
                        .iter()
                        .map(|flow| model::Flow {
                            source: flow.source.clone(),
                            target: flow.target.clone(),
                        })
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>>>()?,
        spectators: mech
            .spectators
            .iter()
            .map(|mol| Ok(Node::Molecule(molecule_from_wire(mol)?).to_string()))
            .collect::<Result<Vec<_>>>()?,
    })
}

// -- JSON surface ------------------------------------------------------

pub fn to_model_json(nodes: &[Node]) -> Result<String> {
    let wire = formula_to_wire(nodes)?;
    serde_json::to_string(&wire).map_err(|error| WireError::Malformed(error.to_string()))
}

pub fn from_model_json(json: &str) -> Result<Vec<Node>> {
    let wire: WireNode =
        serde_json::from_str(json).map_err(|error| WireError::Malformed(error.to_string()))?;
    match wire {
        WireNode::Formula(formula) => formula.nodes.iter().map(node_from_wire).collect(),
        _ => Err(WireError::Malformed(
            "wire root must be a formula node".into(),
        )),
    }
}
