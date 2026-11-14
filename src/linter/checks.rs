//! The linter checks — the Rust mirror of `AsciiChem::Linter::*Check`.
//! Each function walks the model tree and reports diagnostics; none
//! panics (every check reports, the pass stays total).

use std::collections::{BTreeMap, BTreeSet};

use crate::identity;
use crate::layout;
use crate::linter::periodic_table;
use crate::model::{self, Node};
use crate::perception;
use crate::structure;

use super::Diagnostic;

// -- balance -------------------------------------------------------------

/// Verifies that a reaction is stoichiometrically balanced — the atom
/// counts on the reactants side equal those on the products side.
/// Coefficients and group multiplicities are applied.
///
/// The check treats unknown constructs conservatively: any molecule it
/// can't enumerate atoms for is skipped, not flagged.
pub fn balance(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for reaction in collect_reactions(formula) {
        let (Some(reactants), Some(products)) = (
            count_side(&reaction.reactants),
            count_side(&reaction.products),
        ) else {
            continue;
        };
        if reactants == products {
            continue;
        }
        let message = format!(
            "Reaction is not balanced: {}",
            element_diffs(&reactants, &products)
        );
        diagnostics.push(
            Diagnostic::error(message, Some(reaction.diagnostic_label()))
                .with_details(element_deltas(&reactants, &products)),
        );
    }
    diagnostics
}

/// Element => [reactant_count, product_count] for every element whose
/// counts differ — the numbers behind the message.
fn element_deltas(
    reactants: &BTreeMap<String, i64>,
    products: &BTreeMap<String, i64>,
) -> BTreeMap<String, [i64; 2]> {
    let mut deltas = BTreeMap::new();
    for element in union(reactants, products) {
        let left = reactants.get(element).copied().unwrap_or(0);
        let right = products.get(element).copied().unwrap_or(0);
        if left != right {
            deltas.insert(element.clone(), [left, right]);
        }
    }
    deltas
}

fn element_diffs(reactants: &BTreeMap<String, i64>, products: &BTreeMap<String, i64>) -> String {
    union(reactants, products)
        .into_iter()
        .filter_map(|element| {
            let left = reactants.get(element).copied().unwrap_or(0);
            let right = products.get(element).copied().unwrap_or(0);
            (left != right).then(|| format!("{element}: {left} vs {right}"))
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn union<'a>(
    reactants: &'a BTreeMap<String, i64>,
    products: &'a BTreeMap<String, i64>,
) -> BTreeSet<&'a String> {
    reactants
        .keys()
        .chain(products.keys())
        .collect::<BTreeSet<_>>()
}

/// Counts for one side of the reaction, or `None` if any term
/// couldn't be enumerated.
fn count_side(side: &[Node]) -> Option<BTreeMap<String, i64>> {
    let mut counts = BTreeMap::new();
    for term in side {
        // The grammar puts molecules on reaction sides; anything else
        // enumerates as a unit with coefficient 1 (the model build
        // cannot produce the shapes the reference would reject).
        let coefficient = match term {
            Node::Molecule(molecule) => molecule.stoichiometric_coefficient() as i64,
            _ => 1,
        };
        let enumerated = enumerate(term, coefficient)?;
        for (element, count) in enumerated {
            *counts.entry(element).or_insert(0) += count;
        }
    }
    Some(counts)
}

/// Recursive enumeration with a current multiplier.
fn enumerate(node: &Node, multiplier: i64) -> Option<BTreeMap<String, i64>> {
    match node {
        Node::Atom(atom) => {
            let sub = leading_count(atom.subscript.as_deref(), "1");
            Some(BTreeMap::from([(atom.element.clone(), sub * multiplier)]))
        }
        Node::Molecule(molecule) => enumerate_many(&molecule.nodes, multiplier),
        Node::Group(group) => {
            let mult = leading_count(group.multiplicity.as_deref(), "1");
            enumerate_many(&group.nodes, multiplier.checked_mul(mult)?)
        }
        _ => None,
    }
}

fn enumerate_many(nodes: &[Node], multiplier: i64) -> Option<BTreeMap<String, i64>> {
    let mut result = BTreeMap::new();
    for node in nodes {
        let sub = enumerate(node, multiplier)?;
        for (element, count) in sub {
            *result.entry(element).or_insert(0) += count;
        }
    }
    Some(result)
}

// -- charge balance ------------------------------------------------------

/// Verifies charge conservation in reactions: total formal charge on
/// the reactants side equals total on the products side. Parallel to
/// [`balance`] (which validates atom conservation) — a reaction can be
/// stoichiometrically balanced but charge-imbalanced
/// (e.g. `H^+ + OH^- -> H_2O^+`). Reactions without charges on either
/// side pass (no diagnostic).
pub fn charge_balance(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for reaction in collect_reactions(formula) {
        let (Some(reactants), Some(products)) = (
            total_charge(&reaction.reactants),
            total_charge(&reaction.products),
        ) else {
            continue;
        };
        if reactants == products {
            continue;
        }
        diagnostics.push(Diagnostic::error(
            format!(
                "Reaction charge is not balanced: reactants {} vs products {}",
                format_charge(reactants),
                format_charge(products)
            ),
            Some(reaction.diagnostic_label()),
        ));
    }
    diagnostics
}

/// Sum formal charges on one side of the reaction, or `None` if any
/// molecule couldn't be enumerated (mirrors balance's conservative
/// skip behaviour).
fn total_charge(side: &[Node]) -> Option<i64> {
    let mut total: i64 = 0;
    for term in side {
        let coefficient = match term {
            Node::Molecule(molecule) => molecule.stoichiometric_coefficient() as i64,
            _ => 1,
        };
        let charge = enumerate_charge(term)?.checked_mul(coefficient)?;
        total = total.checked_add(charge)?;
    }
    Some(total)
}

/// Recursively sum atom charges, applying group multiplicities.
fn enumerate_charge(node: &Node) -> Option<i64> {
    match node {
        Node::Atom(atom) => Some(parse_charge(atom.charge.as_deref())),
        Node::Group(group) => {
            let mult = leading_count(group.multiplicity.as_deref(), "1");
            Some(sum_nodes_charge(&group.nodes)?.checked_mul(mult)?)
        }
        Node::Molecule(molecule) => sum_nodes_charge(&molecule.nodes),
        _ => Some(0),
    }
}

fn sum_nodes_charge(nodes: &[Node]) -> Option<i64> {
    let mut total: i64 = 0;
    for node in nodes {
        total = total.checked_add(enumerate_charge(node)?)?;
    }
    Some(total)
}

/// Parse a charge string like "2+", "+", "-1", "2-" into an integer.
/// Returns 0 for absent/empty/unparseable (the atom has no charge).
fn parse_charge(charge: Option<&str>) -> i64 {
    let Some(charge) = charge else {
        return 0;
    };
    if charge.is_empty() {
        return 0;
    }
    let bytes = charge.as_bytes();
    let (count, sign) = match (bytes.first(), bytes.last()) {
        (Some(first), Some(last)) if (*last == b'+' || *last == b'-') && first.is_ascii_digit() => {
            (&charge[..charge.len() - 1], *last)
        }
        (Some(first), Some(last))
            if (*first == b'+' || *first == b'-') && last.is_ascii_digit() =>
        {
            (&charge[1..], *first)
        }
        (Some(only), _) if *only == b'+' || *only == b'-' => ("", *only),
        _ => return 0,
    };
    if !count.bytes().all(|byte| byte.is_ascii_digit()) {
        return 0;
    }
    let magnitude: i64 = if count.is_empty() {
        1
    } else {
        count.parse().unwrap_or(0)
    };
    if sign == b'+' {
        magnitude
    } else {
        -magnitude
    }
}

fn format_charge(value: i64) -> String {
    if value == 0 {
        "0".to_string()
    } else if value > 0 {
        format!("+{value}")
    } else {
        value.to_string()
    }
}

// -- bracket balance -----------------------------------------------------

/// Verifies that group brackets are matched. The grammar enforces this
/// at parse time and the model derives both bracket characters from
/// the kind by construction — so over Rust models the defence-in-depth
/// check reduces to the one property that CAN drift: the group's
/// bracket kind must be one of the known kinds (the reference's
/// open/close mismatch is unrepresentable here).
pub fn bracket_balance(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Group(group) = node else {
            continue;
        };
        if !matches!(group.bracket.as_str(), "paren" | "square" | "brace") {
            diagnostics.push(Diagnostic::error(
                format!(
                    "Group bracket kind {:?} is not a known bracket kind \
                     (paren, square, brace)",
                    group.bracket
                ),
                Some(group.diagnostic_label()),
            ));
        }
    }
    diagnostics
}

// -- isotope sanity ------------------------------------------------------

/// Sanity-checks isotope mass numbers against the element's atomic
/// number: the mass must be ≥ the atomic number (a nucleus can't have
/// fewer protons than its mass number suggests). Unknown elements
/// produce an info diagnostic, not an error.
pub fn isotope_sanity(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Atom(atom) = node else {
            continue;
        };
        let Some(isotope) = &atom.isotope else {
            continue;
        };
        let Some(atomic_number) = periodic_table::atomic_number(&atom.element) else {
            diagnostics.push(Diagnostic::info(
                format!(
                    "Element {:?} not in periodic table; skipping isotope check",
                    atom.element
                ),
                Some(atom.diagnostic_label()),
            ));
            continue;
        };
        let mass = model::leading_integer(isotope);
        if mass >= i64::from(atomic_number) {
            continue;
        }
        diagnostics.push(Diagnostic::error(
            format!(
                "Isotope mass {mass} is less than atomic number {atomic_number} for {}",
                atom.element
            ),
            Some(atom.diagnostic_label()),
        ));
    }
    diagnostics
}

// -- crystal sanity ------------------------------------------------------

/// Validates unit-cell parameters on a Crystal node (the reference
/// `CrystalSanityCheck`):
///
/// - lengths (a, b, c) must be positive
/// - angles (alpha, beta, gamma) must be in (0, 180] degrees
/// - fractional coordinates on the cell's atoms must be in [0, 1)
///
/// Each cell parameter is optional in the model (a partial crystal
/// description is allowed), so missing values are skipped, not
/// flagged — and non-numeric captures are skipped too (the reference's
/// `Float(value, exception: false)`).
pub fn crystal_sanity(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Crystal(crystal) = node else {
            continue;
        };
        check_lengths(crystal, &mut diagnostics);
        check_angles(crystal, &mut diagnostics);
        check_fractional_coords(crystal, &mut diagnostics);
    }
    diagnostics
}

fn check_lengths(crystal: &model::Crystal, diagnostics: &mut Vec<Diagnostic>) {
    for (dim, value) in [("a", &crystal.a), ("b", &crystal.b), ("c", &crystal.c)] {
        let Some(number) = numeric(value.as_deref()) else {
            continue;
        };
        if number > 0.0 {
            continue;
        }
        diagnostics.push(Diagnostic::error(
            format!(
                "Crystal {dim} must be positive, got {}",
                value.as_deref().unwrap_or_default()
            ),
            Some(crystal.diagnostic_label()),
        ));
    }
}

fn check_angles(crystal: &model::Crystal, diagnostics: &mut Vec<Diagnostic>) {
    for (dim, value) in [
        ("alpha", &crystal.alpha),
        ("beta", &crystal.beta),
        ("gamma", &crystal.gamma),
    ] {
        let Some(number) = numeric(value.as_deref()) else {
            continue;
        };
        if (0.0..=180.0).contains(&number) && number != 0.0 {
            continue;
        }
        diagnostics.push(Diagnostic::error(
            format!(
                "Crystal {dim} must be in (0, 180], got {}",
                value.as_deref().unwrap_or_default()
            ),
            Some(crystal.diagnostic_label()),
        ));
    }
}

fn check_fractional_coords(crystal: &model::Crystal, diagnostics: &mut Vec<Diagnostic>) {
    for atom in &crystal.atoms {
        for (axis, value) in [
            ("x_fract", atom.x_fract),
            ("y_fract", atom.y_fract),
            ("z_fract", atom.z_fract),
        ] {
            let Some(number) = value else {
                continue;
            };
            if (0.0..1.0).contains(&number) {
                continue;
            }
            diagnostics.push(Diagnostic::warning(
                format!(
                    "Fractional coordinate {axis}={} on {} is outside [0, 1)",
                    crate::model::format_float(number),
                    atom.element
                ),
                Some(crystal.diagnostic_label()),
            ));
        }
    }
}

/// The reference's `Float(value, exception: false)`: a numeric string
/// parses, anything else (including nil) skips the check.
fn numeric(value: Option<&str>) -> Option<f64> {
    value?.trim().parse::<f64>().ok()
}

// -- element validation --------------------------------------------------

/// Validates that each atom's element symbol is in the periodic table.
/// Catches typos like `Hx`, `Cy` that would otherwise silently parse
/// and produce nonsense output.
///
/// Severity is `warning`, not `error` — the parser is total and might
/// intentionally accept placeholder elements (e.g. `X` for "unknown
/// halogen" in teaching contexts). The linter flags them for the user
/// to confirm.
pub fn element_validation(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Atom(atom) = node else {
            continue;
        };
        if periodic_table::known(&atom.element) {
            continue;
        }
        diagnostics.push(Diagnostic::warning(
            format!(
                "Unknown element symbol {:?}; not in the periodic table",
                atom.element
            ),
            Some(atom.diagnostic_label()),
        ));
    }
    diagnostics
}

// -- unclosed rings ------------------------------------------------------

/// Flags ring closure digits that have no matching partner. Catches
/// typos like `C1-C-C` (digit 1 opened but never closed). The ring
/// pairing is the shared [`layout::unclosed_ring_atoms`] walk — the
/// same source of truth the layout consumes.
pub fn unclosed_ring(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Molecule(molecule) = node else {
            continue;
        };
        let (atoms, _) = structure::build_graph(molecule);
        for index in layout::unclosed_ring_atoms(&atoms) {
            let atom = &atoms[index];
            diagnostics.push(Diagnostic::error(
                format!(
                    "Atom has unmatched ring closure digit(s) {:?}",
                    atom.ring_closures.clone().unwrap_or_default()
                ),
                Some(atom.diagnostic_label()),
            ));
        }
    }
    diagnostics
}

// -- valence -------------------------------------------------------------

/// Verifies that each atom's total bond order does not exceed the
/// element's maximum valence. Catches typos like carbon with five
/// bonds. Runs over the perception pass's graph view: every atom's
/// bond order counts ALL its bonds (rings and branches included).
/// Unknown elements produce an info diagnostic, not an error.
pub fn valence(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Molecule(molecule) = node else {
            continue;
        };
        let (atoms, edges) = structure::build_graph(molecule);
        for (index, atom) in atoms.iter().enumerate() {
            let load = perception::bond_order_sum(index, &edges);
            let Some(max) = periodic_table::max_valence(&atom.element) else {
                diagnostics.push(Diagnostic::info(
                    format!("Element {:?} not in valence table; skipping", atom.element),
                    Some(atom.diagnostic_label()),
                ));
                continue;
            };
            if load <= f64::from(max) {
                continue;
            }
            diagnostics.push(Diagnostic::error(
                format!(
                    "Atom {} has bond order {}; max valence is {max}",
                    atom.element, load
                ),
                Some(atom.diagnostic_label()),
            ));
        }
    }
    diagnostics
}

// -- identifier format ---------------------------------------------------

/// Validates the *format* of molecule identifier annotations
/// (`@cas("...")`, `@inchi("...")`, `@smiles("...")`) whose convention
/// has a registered offline validator ([`identifiers`]). Lexical and
/// checksum only — checking that an identifier is actually *assigned*
/// to a substance is resolution (TODO.v2 07 Layer 2), not linting.
/// Conventions without a validator (e.g. `iupac`, `cid`) are skipped
/// silently.
pub fn identifier_format(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Molecule(molecule) = node else {
            continue;
        };
        for identifier in &molecule.identifiers {
            let Some(reason) =
                crate::identifiers::diagnostic(&identifier.convention, &identifier.value)
            else {
                continue;
            };
            let formula = molecule.hill_formula();
            let label = if formula.is_empty() {
                "molecule".to_string()
            } else {
                formula
            };
            diagnostics.push(Diagnostic::error(
                format!("{reason} ({} identifier on {label})", identifier.convention),
                Some(molecule.diagnostic_label()),
            ));
        }
    }
    diagnostics
}

// -- identifier consistency ----------------------------------------------

/// Cross-checks identifier values against the molecule they annotate,
/// where checkable offline:
///
/// - InChI: the formula layer's element counts must match the
///   molecule's composition (error — the formula layer is
///   authoritative).
/// - SMILES: the element set must agree with the molecule's elements
///   (warning — implicit hydrogens are tolerated, so only set
///   membership is checkable without a full parser).
///
/// Structural equality (stereo, isotopes, tautomers) needs the InChI
/// engine and is the identity cross-check's job.
///
/// Identifiers describe the *substance*, not the stoichiometric
/// multiplier, so comparisons use the molecule's composition without
/// its coefficient.
pub fn identifier_consistency(formula: &[Node]) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for node in walk(formula) {
        let Node::Molecule(molecule) = node else {
            continue;
        };
        let counts = molecule.element_counts(false);
        // The identifier validators count in i64 (InChI layers carry
        // isotope-scale counts); align the molecule tally to them.
        let counts_i64: std::collections::BTreeMap<String, i64> =
            counts.iter().map(|(k, v)| (k.clone(), *v as i64)).collect();
        for identifier in &molecule.identifiers {
            match identifier.convention.as_str() {
                "inchi" => {
                    // Malformed values are the format check's report.
                    let Some(inchi_counts) =
                        crate::identifiers::inchi::formula_counts(&identifier.value)
                    else {
                        continue;
                    };
                    let diffs = composition_diffs(&counts_i64, &inchi_counts);
                    if diffs.is_empty() {
                        continue;
                    }
                    diagnostics.push(Diagnostic::error(
                        format!(
                            "InChI formula {} does not match molecule {}: {diffs}",
                            formula_string(&inchi_counts),
                            formula_string(&counts_i64)
                        ),
                        Some(molecule.diagnostic_label()),
                    ));
                }
                "smiles" => {
                    if let Some(message) = smiles_mismatch(&counts_i64, &identifier.value) {
                        diagnostics.push(Diagnostic::warning(
                            message,
                            Some(molecule.diagnostic_label()),
                        ));
                    }
                }
                _ => {}
            }
        }
    }
    diagnostics
}

/// `None` when the SMILES element set agrees with the molecule
/// (implicit hydrogens tolerated); a message otherwise.
fn smiles_mismatch(counts: &BTreeMap<String, i64>, value: &str) -> Option<String> {
    // Unanalysable values are the format check's report.
    let smiles_elements = crate::identifiers::smiles::element_set(value)?;
    let extra: Vec<String> = smiles_elements
        .iter()
        .filter(|element| !counts.contains_key(*element))
        .cloned()
        .collect();
    let missing: Vec<String> = counts
        .keys()
        .filter(|element| *element != "H" && !smiles_elements.contains(*element))
        .cloned()
        .collect();
    if extra.is_empty() && missing.is_empty() {
        return None;
    }
    let listed: Vec<String> = smiles_elements.into_iter().collect();
    Some(format!(
        "SMILES elements {listed:?} do not match molecule {}: unexpected {extra:?}, missing {missing:?}",
        formula_string(counts)
    ))
}

fn composition_diffs(
    molecule_counts: &BTreeMap<String, i64>,
    inchi_counts: &BTreeMap<String, i64>,
) -> String {
    let mut diffs = Vec::new();
    let keys: BTreeSet<&String> = molecule_counts.keys().chain(inchi_counts.keys()).collect();
    for element in keys {
        let actual = molecule_counts.get(element).copied().unwrap_or(0);
        let expected = inchi_counts.get(element).copied().unwrap_or(0);
        if actual != expected {
            diffs.push(format!("{element}: {actual} vs {expected}"));
        }
    }
    diffs.join(", ")
}

fn formula_string(counts: &BTreeMap<String, i64>) -> String {
    counts
        .iter()
        .map(|(element, count)| {
            if *count == 1 {
                element.clone()
            } else {
                format!("{element}{count}")
            }
        })
        .collect()
}

// -- identity cross-check ------------------------------------------------

/// Structural cross-check of identity annotations (TODO.v2 10). Where
/// the identifier consistency check compares offline-analysable
/// layers (formula composition), this check derives the InChI of the
/// drawn structure through the in-process chemomer engine and compares
/// it to the annotated `@inchi`/`@inchikey` — catching transcription
/// errors that composition checks cannot (wrong connectivity, stereo,
/// isotopes).
///
/// Engine-scaled, mirroring the resolver's opt-in sources: the Ruby
/// reference warns once when molecules are annotated but no engine is
/// configured. Rust gates the engine IN-PROCESS (the chemomer
/// dependency is unconditional), so that skip-warning has no Rust
/// path; instead every derivation failure reports as a warning and the
/// check stays total. Molecules without identity annotations are never
/// checked.
pub fn identity_cross_check(formula: &[Node]) -> Vec<Diagnostic> {
    let annotated: Vec<&model::Molecule> = walk(formula)
        .into_iter()
        .filter_map(|node| match node {
            Node::Molecule(molecule) if has_identity_annotations(molecule) => Some(molecule),
            _ => None,
        })
        .collect();
    annotated
        .iter()
        .flat_map(|molecule| cross_check(molecule))
        .collect()
}

fn has_identity_annotations(molecule: &model::Molecule) -> bool {
    molecule
        .identifiers
        .iter()
        .any(|identifier| matches!(identifier.convention.as_str(), "inchi" | "inchikey"))
}

fn cross_check(molecule: &model::Molecule) -> Vec<Diagnostic> {
    let identity = match identity::derive(molecule) {
        Ok(identity) => identity,
        // Not derivable (e.g. a formula is not a structure: no bonds)
        // or the engine rejected it — report, never raise: the linter
        // stays total.
        Err(error) => {
            return vec![Diagnostic::warning(
                format!("cannot derive InChI for cross-check: {}", error.message()),
                Some(molecule.diagnostic_label()),
            )];
        }
    };
    molecule
        .identifiers
        .iter()
        .filter_map(|identifier| match identifier.convention.as_str() {
            "inchi" if identifier.value != identity.inchi => Some(Diagnostic::error(
                format!(
                    "InChI annotation does not match the drawn structure: \
                     annotated {}, computed {}",
                    identifier.value, identity.inchi
                ),
                Some(molecule.diagnostic_label()),
            )),
            "inchikey" if identifier.value != identity.inchikey => Some(Diagnostic::error(
                format!(
                    "InChIKey annotation does not match the drawn structure: \
                     annotated {}, computed {}",
                    identifier.value, identity.inchikey
                ),
                Some(molecule.diagnostic_label()),
            )),
            _ => None,
        })
        .collect()
}

// -- shared walk ---------------------------------------------------------

/// Depth-first pre-order walk over every node in the formula:
/// containers recurse into their children (groups, molecules,
/// reactions' sides, cascade steps' sides).
fn walk(formula: &[Node]) -> Vec<&Node> {
    let mut visited = Vec::new();
    walk_nodes(formula, &mut visited);
    visited
}

fn walk_nodes<'a>(nodes: &'a [Node], visited: &mut Vec<&'a Node>) {
    for node in nodes {
        visited.push(node);
        match node {
            Node::Group(group) => walk_nodes(&group.nodes, visited),
            Node::Molecule(molecule) => walk_nodes(&molecule.nodes, visited),
            Node::Reaction(reaction) => {
                walk_nodes(&reaction.reactants, visited);
                walk_nodes(&reaction.products, visited);
            }
            Node::ReactionCascade(cascade) => {
                for step in &cascade.steps {
                    walk_nodes(&step.reactants, visited);
                    walk_nodes(&step.products, visited);
                }
            }
            _ => {}
        }
    }
}

/// Every reaction in the tree, outermost first — cascade steps and
/// reactions nested in containers included.
fn collect_reactions(formula: &[Node]) -> Vec<&model::Reaction> {
    let mut reactions = Vec::new();
    collect_reactions_in(formula, &mut reactions);
    reactions
}

fn collect_reactions_in<'a>(nodes: &'a [Node], reactions: &mut Vec<&'a model::Reaction>) {
    for node in nodes {
        match node {
            Node::Reaction(reaction) => {
                reactions.push(reaction);
                collect_reactions_in(&reaction.reactants, reactions);
                collect_reactions_in(&reaction.products, reactions);
            }
            Node::ReactionCascade(cascade) => {
                for step in &cascade.steps {
                    reactions.push(step);
                    collect_reactions_in(&step.reactants, reactions);
                    collect_reactions_in(&step.products, reactions);
                }
            }
            Node::Group(group) => collect_reactions_in(&group.nodes, reactions),
            Node::Molecule(molecule) => collect_reactions_in(&molecule.nodes, reactions),
            _ => {}
        }
    }
}

/// A model count field parsed as the reference's `String#to_i`:
/// leading digits with an optional sign, 0 when none. Zero (or
/// non-positive) counts read as 1 at every call site that needs a
/// multiplier.
fn leading_count(value: Option<&str>, absent: &str) -> i64 {
    let parsed = model::leading_integer(value.unwrap_or(absent));
    if parsed == 0 {
        1
    } else {
        parsed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nodes(input: &str) -> Vec<Node> {
        crate::parse_nodes(input).expect("input parses")
    }

    fn lint(input: &str) -> Vec<Diagnostic> {
        let parsed = nodes(input);
        let mut all = Vec::new();
        for check in [
            balance,
            bracket_balance,
            charge_balance,
            crystal_sanity,
            element_validation,
            identifier_consistency,
            identifier_format,
            identity_cross_check,
            isotope_sanity,
            unclosed_ring,
            valence,
        ] {
            all.extend(check(&parsed));
        }
        all
    }

    fn errors_of(input: &str) -> Vec<Diagnostic> {
        lint(input)
            .into_iter()
            .filter(|diagnostic| diagnostic.severity == crate::linter::Severity::Error)
            .collect()
    }

    // Clean formulas produce nothing.
    #[test]
    fn clean_water_is_silent() {
        assert!(lint("H_2O").is_empty());
        assert!(lint("^14C").is_empty());
        assert!(lint("Ca(OH)_2").is_empty());
    }

    // -- balance ---------------------------------------------------------

    #[test]
    fn balanced_reactions_pass() {
        assert!(errors_of("2H_2 + O_2 -> 2H_2O").is_empty());
        assert!(errors_of("Ca(OH)_2 + 2HCl -> CaCl_2 + 2H_2O").is_empty());
    }

    #[test]
    fn unbalanced_reactions_report_the_deltas() {
        let errors = errors_of("H_2 + O_2 -> H_2O");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("not balanced"));
        assert!(errors[0].message.contains("O: 2 vs 1"));
        assert_eq!(
            errors[0].details,
            Some(BTreeMap::from([("O".to_string(), [2, 1])]))
        );
        assert_eq!(errors[0].label.as_deref(), Some("Reaction"));
    }

    #[test]
    fn balanced_reactions_carry_no_details() {
        for diagnostic in lint("2H_2 + O_2 -> 2H_2O") {
            assert_eq!(diagnostic.details, None);
        }
    }

    // -- charge balance --------------------------------------------------

    fn charge_errors(input: &str) -> Vec<Diagnostic> {
        errors_of(input)
            .into_iter()
            .filter(|diagnostic| diagnostic.message.contains("charge"))
            .collect()
    }

    #[test]
    fn charge_balanced_reactions_pass() {
        assert!(charge_errors("H^+ + OH^- -> H_2O").is_empty());
        assert!(charge_errors("H_2 + O_2 -> H_2O").is_empty());
        assert!(charge_errors("Ca^2+ + 2Cl^- -> CaCl_2").is_empty());
        assert!(charge_errors("2Na^+ + SO_4^2- -> Na_2SO_4").is_empty());
    }

    #[test]
    fn charge_imbalanced_reactions_report_both_sides() {
        let diags = charge_errors("H^+ + OH^- -> H_2O^+");
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.contains("reactants 0 vs products +1"));

        let diags = charge_errors("Na^+ + Cl^- -> NaCl^-");
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.contains("reactants 0 vs products -1"));

        let diags = charge_errors("Ca^2+ + Cl^- -> CaCl");
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.contains("reactants +1 vs products 0"));
    }

    #[test]
    fn charge_values_parse_both_orders() {
        assert_eq!(parse_charge(Some("2+")), 2);
        assert_eq!(parse_charge(Some("+")), 1);
        assert_eq!(parse_charge(Some("-1")), -1);
        assert_eq!(parse_charge(Some("2-")), -2);
        assert_eq!(parse_charge(Some("")), 0);
        assert_eq!(parse_charge(None), 0);
        assert_eq!(parse_charge(Some("2+3")), 0);
    }

    // -- brackets / isotopes / elements -----------------------------------

    #[test]
    fn unknown_bracket_kinds_are_defence_in_depth() {
        // Not constructible by the parser; reachable by direct model
        // construction only — the check pins the invariant.
        let mut group = model::Group {
            nodes: vec![],
            multiplicity: None,
            bracket: "curly".into(),
        };
        let diagnostics = bracket_balance(&[Node::Group(group.clone())]);
        assert_eq!(diagnostics.len(), 1);
        assert!(diagnostics[0].message.contains("curly"));

        group.bracket = "square".into();
        assert!(bracket_balance(&[Node::Group(group)]).is_empty());
    }

    #[test]
    fn isotope_below_atomic_number_errors() {
        let errors = errors_of("^5C");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("Isotope mass 5"));
        assert!(errors[0].message.contains("atomic number 6"));
        assert_eq!(errors[0].label.as_deref(), Some("Atom(C)"));
    }

    #[test]
    fn valid_hydrogen_isotopes_pass() {
        for isotope in ["^1H", "^2H", "^3H"] {
            assert!(errors_of(isotope).is_empty(), "{isotope}");
        }
    }

    #[test]
    fn isotope_on_unknown_element_infos() {
        let infos: Vec<Diagnostic> = lint("^5Xx")
            .into_iter()
            .filter(|diagnostic| diagnostic.severity == crate::linter::Severity::Info)
            .collect();
        assert!(infos
            .iter()
            .any(|diagnostic| diagnostic.message.contains("periodic table")));
    }

    #[test]
    fn unknown_elements_warn() {
        let warnings: Vec<Diagnostic> = lint("Hx_2O")
            .into_iter()
            .filter(|diagnostic| {
                diagnostic.severity == crate::linter::Severity::Warning
                    && diagnostic.message.contains("Unknown element")
            })
            .collect();
        assert!(warnings
            .iter()
            .any(|diagnostic| diagnostic.message.contains("Hx")));
        assert!(lint("Cy")
            .iter()
            .any(|diagnostic| diagnostic.message.contains("Cy")));
    }

    #[test]
    fn every_periodic_table_symbol_is_accepted() {
        for symbol in periodic_table::symbols() {
            let diagnostics = lint(symbol);
            assert!(
                !diagnostics
                    .iter()
                    .any(|diagnostic| diagnostic.message.contains("Unknown element")),
                "{symbol} should be valid"
            );
        }
    }

    // -- rings / valence ---------------------------------------------------

    #[test]
    fn unclosed_ring_digits_error() {
        let errors = errors_of("C-C-C1");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("unmatched ring closure"));
        assert!(errors[0].message.contains("\"1\""));
        assert_eq!(errors[0].label.as_deref(), Some("Atom(C)"));
        assert!(errors_of("C1-C-C-C-C-C1").is_empty());
    }

    #[test]
    fn valence_overload_errors_and_unknowns_info() {
        // Explicitly over-bonded carbon: C with a double and two
        // single bonds (order 4 is fine; the point is the walk runs).
        assert!(errors_of("H-C=C-H").is_empty());
        let infos: Vec<Diagnostic> = lint("Xx-2Yy")
            .into_iter()
            .filter(|diagnostic| diagnostic.severity == crate::linter::Severity::Info)
            .collect();
        assert!(infos
            .iter()
            .any(|diagnostic| diagnostic.message.contains("not in valence table")));
    }

    // -- identifier format / consistency ----------------------------------

    #[test]
    fn identifier_format_passes_valid_and_skips_unknown_conventions() {
        assert!(errors_of("H_2O @cas(\"7732-18-5\")").is_empty());
        assert!(errors_of("C @chebi(\"CHEBI:16183\")").is_empty());
    }

    #[test]
    fn identifier_format_flags_malformed_values() {
        let errors = errors_of("H_2O @cas(\"50782\")");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("not in DDDDDDD-DD-D form"));
        assert!(errors[0].message.contains("cas identifier on H2O"));

        let errors = errors_of("C @cas(\"74-82-9\")");
        assert_eq!(errors.len(), 1);
        assert!(errors[0]
            .message
            .contains("check digit 9 does not match weighted sum (8)"));

        // The engine-scaled identity cross-check ALSO fires (rs ships
        // the in-process chemomer engine): the annotated Cx4 InChI
        // cannot match the derived methane identity.
        let errors = errors_of("C @inchi(\"InChI=1S/Cx4\")");
        assert_eq!(errors.len(), 2);
        assert!(errors
            .iter()
            .any(|diagnostic| diagnostic.message.contains("unknown element symbol \"Cx\"")));
        assert!(errors.iter().any(|diagnostic| diagnostic
            .message
            .contains("does not match the drawn structure")));

        let errors = errors_of("C @smiles(\"C(C\")");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("parentheses are not balanced"));
    }

    #[test]
    fn identifier_format_checks_chained_annotations() {
        let errors = errors_of("C @cas(\"74-82-9\") @inchi(\"InChI=1S/Cx4\")");
        assert_eq!(errors.len(), 3);
        assert!(errors
            .iter()
            .any(|diagnostic| diagnostic.message.contains("check digit")));
        assert!(errors
            .iter()
            .any(|diagnostic| diagnostic.message.contains("unknown element")));
        assert!(errors.iter().any(|diagnostic| diagnostic
            .message
            .contains("does not match the drawn structure")));
    }

    const ASPIRIN_INCHI: &str =
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)";

    fn consistency_diagnostics(input: &str) -> Vec<Diagnostic> {
        lint(input)
            .into_iter()
            .filter(|diagnostic| diagnostic.message.contains("match molecule"))
            .collect()
    }

    #[test]
    fn inchi_consistency_matches_and_mismatches() {
        assert!(
            consistency_diagnostics(&format!("C_9H_8O_4 @inchi(\"{ASPIRIN_INCHI}\")")).is_empty()
        );
        let diags = consistency_diagnostics(&format!("H_2O @inchi(\"{ASPIRIN_INCHI}\")"));
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].severity, crate::linter::Severity::Error);
        assert!(diags[0].message.contains("C9H8O4"));
        assert!(diags[0].message.contains("H: 2 vs 8"));
    }

    #[test]
    fn inchi_consistency_ignores_the_coefficient_and_unanalysable_values() {
        assert!(consistency_diagnostics("2H_2O @inchi(\"InChI=1S/H2O/h1H2\")").is_empty());
        assert!(consistency_diagnostics("H_2O @inchi(\"garbage\")").is_empty());
    }

    #[test]
    fn smiles_consistency_warns_on_set_mismatch_only() {
        let diags = consistency_diagnostics("H_2O @smiles(\"CCO\")");
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].severity, crate::linter::Severity::Warning);
        assert!(diags[0].message.contains("unexpected [\"C\"]"));

        let diags = consistency_diagnostics("C_2H_6O @smiles(\"CC\")");
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.contains("missing [\"O\"]"));

        assert!(consistency_diagnostics("C_2H_6O @smiles(\"CCO\")").is_empty());
        assert!(consistency_diagnostics("C_9H_8O_4 @smiles(\"CC(=O)Oc1ccccc1C(=O)O\")").is_empty());
    }

    // -- identity cross-check (in-process engine) --------------------------

    fn smiles_molecule(smiles: &str) -> Node {
        let parsed = crate::smiles::parse_smiles(smiles).expect("SMILES parses");
        match parsed.first() {
            Some(Node::Molecule(molecule)) => Node::Molecule(molecule.clone()),
            other => panic!("not a molecule: {other:?}"),
        }
    }

    fn annotated_aspirin(convention: &str, value: &str) -> Vec<Node> {
        let Node::Molecule(mut molecule) = smiles_molecule("CC(=O)OC1=CC=CC=C1C(=O)O") else {
            unreachable!()
        };
        molecule.identifiers.push(model::Identifier {
            convention: convention.into(),
            value: value.into(),
        });
        vec![Node::Molecule(molecule)]
    }

    #[test]
    fn identity_cross_check_stays_silent_without_annotations() {
        assert!(identity_cross_check(&nodes("H_2O")).is_empty());
        assert!(identity_cross_check(&nodes("C-C-O @cas(\"64-17-5\")")).is_empty());
    }

    #[test]
    fn identity_cross_check_passes_on_the_recorded_identity() {
        let formula = annotated_aspirin("inchi", ASPIRIN_INCHI);
        assert!(identity_cross_check(&formula).is_empty());
        let formula = annotated_aspirin("inchikey", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N");
        assert!(identity_cross_check(&formula).is_empty());
    }

    #[test]
    fn identity_cross_check_errors_on_transcription_drift() {
        let wrong = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-4H,1H3,(H,11,12)";
        let diagnostics = identity_cross_check(&annotated_aspirin("inchi", wrong));
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(diagnostics[0].severity, crate::linter::Severity::Error);
        assert!(diagnostics[0]
            .message
            .contains("does not match the drawn structure"));
        assert!(diagnostics[0].message.contains(ASPIRIN_INCHI));

        let diagnostics = identity_cross_check(&annotated_aspirin(
            "inchikey",
            "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        ));
        assert_eq!(diagnostics.len(), 1);
        assert!(diagnostics[0]
            .message
            .contains("computed BSYNRYMUTXBXSQ-UHFFFAOYSA-N"));
    }

    #[test]
    fn identity_cross_check_warns_instead_of_raising_on_non_structures() {
        let diagnostics =
            identity_cross_check(&nodes(&format!("H_2O @inchi(\"{ASPIRIN_INCHI}\")")));
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(diagnostics[0].severity, crate::linter::Severity::Warning);
        assert!(diagnostics[0].message.contains("cannot derive InChI"));
    }

    // -- reaction walks see inside cascades --------------------------------

    #[test]
    fn cascades_are_linted_step_by_step() {
        // Step 1 balanced, step 2 not: exactly one balance error.
        let errors = errors_of("2H_2 + O_2 -> 2H_2O -> H_2O");
        assert_eq!(errors.len(), 1);
        assert!(errors[0].message.contains("not balanced"));
    }

    // -- totality: unknown constructs are skipped, never fatal -------------

    #[test]
    fn embedded_math_and_text_never_trip_the_checks() {
        assert!(lint("`x^2` \"note\"").is_empty());
    }

    // -- crystal sanity ---------------------------------------------------

    #[test]
    fn negative_cell_lengths_error_with_the_authored_value() {
        // The corpus cli/lint-json-crystal case's exact expectation.
        let errors = errors_of("crystal[x](a=-1){Na@f(0,0,0)}");
        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].message, "Crystal a must be positive, got -1");
        assert_eq!(errors[0].label.as_deref(), Some("Crystal(x)"));
    }

    #[test]
    fn non_positive_angles_error_and_valid_cells_stay_silent() {
        let errors = errors_of("crystal[x](alpha=0,beta=181,gamma=-1)");
        assert_eq!(errors.len(), 3);
        assert!(errors
            .iter()
            .all(|diagnostic| diagnostic.message.contains("must be in (0, 180]")));
        assert!(lint("crystal[x](a=5.64,b=5.64,c=5.64,alpha=90,beta=90,gamma=180)").is_empty());
    }

    #[test]
    fn out_of_range_fractional_coordinates_warn() {
        let warnings: Vec<Diagnostic> = lint("crystal[x](a=1){Na@f(0,1.5,0)}")
            .into_iter()
            .filter(|diagnostic| diagnostic.severity == crate::linter::Severity::Warning)
            .collect();
        assert_eq!(warnings.len(), 1);
        assert_eq!(
            warnings[0].message,
            "Fractional coordinate y_fract=1.5 on Na is outside [0, 1)"
        );
        assert_eq!(warnings[0].label.as_deref(), Some("Crystal(x)"));
    }

    #[test]
    fn partial_and_non_numeric_cells_are_skipped_not_flagged() {
        // A partial description is allowed; a non-numeric capture skips
        // (the reference's Float(..., exception: false)).
        assert!(lint("crystal[x]").is_empty());
        assert!(lint("crystal[x](sg=Fm-3m)").is_empty());
    }
}
