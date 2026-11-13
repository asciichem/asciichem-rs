//! The mechanism step graph — the Rust port of the reference
//! implementation's `AsciiChem::MechanismGraph` (TODO.v2/13 item 4).
//!
//! `mechanism{}` entries stay verbatim (label + reaction text) in the
//! model; this service derives the graph the text implies — each
//! step's reaction parsed with the same grammar, species that appear
//! on both sides of a step recognised as catalysts, and steps linked
//! where a product species of an earlier step is a reactant of a
//! later one (a cascade's entry side is its head and its exit side is
//! its last step). Everything is computed; nothing mutates. A step
//! whose value is not a reaction (prose, "proton transfer") keeps a
//! `None` reaction and takes part in no links — never an error.
//!
//! Re-parsing runs through [`crate::parse_nodes`]: fine from ordinary
//! callers, but like every nested parse it must not run INSIDE a
//! deep-stack worker job (the worker is single-threaded; nesting
//! self-deadlocks) — see the crystal build's note.

use crate::model::{self, Node};

/// One derived step: the verbatim label and flows plus the parsed
/// reaction (the model node when the step's value parses as a
/// reaction or cascade) and the catalyst species the reaction's two
/// sides share.
#[derive(Debug, Clone, PartialEq)]
pub struct Step {
    pub label: String,
    pub reaction: Option<Node>,
    pub catalysts: Vec<String>,
    pub flows: Vec<model::Flow>,
}

/// A link between two steps: the earlier step's product species that
/// the later step consumes as reactants (coefficient-free names, in
/// the producing step's product order).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Link {
    pub from: String,
    pub to: String,
    pub species: Vec<String>,
}

/// The derived graph: every step plus the links between them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Graph {
    pub steps: Vec<Step>,
    pub links: Vec<Link>,
}

/// Derives the step graph of a mechanism. Computed on the fly — the
/// canonical text and the model are unchanged.
pub fn for_mechanism(mechanism: &model::Mechanism) -> Graph {
    let steps: Vec<Step> = mechanism.steps.iter().map(graph_step).collect();
    let links = links(&steps);
    Graph { steps, links }
}

fn graph_step(step: &model::MechanismStep) -> Step {
    let reaction = parse_reaction(&step.reaction);
    let catalysts = match sides_of(reaction.as_ref()) {
        Some((reactants, products)) => intersect(&reactants, &products),
        None => Vec::new(),
    };
    Step {
        label: step.label.clone(),
        reaction,
        catalysts,
        flows: step.flows.clone(),
    }
}

fn links(steps: &[Step]) -> Vec<Link> {
    let mut links = Vec::new();
    for (index, from) in steps.iter().enumerate() {
        let Some((_, from_products)) = sides_of(from.reaction.as_ref()) else {
            continue;
        };
        for to in &steps[index + 1..] {
            let Some((to_reactants, _)) = sides_of(to.reaction.as_ref()) else {
                continue;
            };
            let shared = intersect(&from_products, &to_reactants);
            if !shared.is_empty() {
                links.push(Link {
                    from: from.label.clone(),
                    to: to.label.clone(),
                    species: shared,
                });
            }
        }
    }
    links
}

/// The step's reaction text parsed with the same grammar; kept only
/// when it builds a reaction or a cascade (prose stays `None`).
fn parse_reaction(text: &str) -> Option<Node> {
    match crate::parse_nodes(text) {
        Ok(nodes) => match nodes.first() {
            Some(node @ Node::Reaction(_)) | Some(node @ Node::ReactionCascade(_)) => {
                Some(node.clone())
            }
            _ => None,
        },
        Err(_) => None,
    }
}

/// The coefficient-free species names per side: `(reactants,
/// products)`. A cascade's entry side is its head and its exit side is
/// its last step; a reaction-shaped node with no cascade steps (not
/// producible by the build) yields `None` — never a panic.
fn sides_of(reaction: Option<&Node>) -> Option<(Vec<String>, Vec<String>)> {
    match reaction? {
        Node::Reaction(reaction) => Some((
            side_keys(&reaction.reactants),
            side_keys(&reaction.products),
        )),
        Node::ReactionCascade(cascade) => {
            let first = cascade.steps.first()?;
            let last = cascade.steps.last()?;
            Some((side_keys(&first.reactants), side_keys(&last.products)))
        }
        _ => None,
    }
}

/// Species names of one side, coefficient-free (the reference's
/// `side_stoichiometry` keys: the term's canonical text with its
/// leading integer stripped), in first-occurrence order.
fn side_keys(side: &[Node]) -> Vec<String> {
    let mut keys = Vec::new();
    for term in side {
        let name = strip_coefficient(&term.to_string());
        if !keys.contains(&name) {
            keys.push(name);
        }
    }
    keys
}

/// `"2A"` matches `"A"` — the leading integer is the stoichiometric
/// coefficient, never part of the species name.
fn strip_coefficient(text: &str) -> String {
    text.trim_start_matches(|c: char| c.is_ascii_digit())
        .to_string()
}

/// The receiver-order intersection (the reference's `Array#&`): the
/// left side's names that the right side also has, deduplicated.
fn intersect(left: &[String], right: &[String]) -> Vec<String> {
    let mut shared = Vec::new();
    for name in left {
        if right.contains(name) && !shared.contains(name) {
            shared.push(name.clone());
        }
    }
    shared
}
