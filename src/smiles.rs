//! SMILES ingestion and emission — the Rust mirror of
//! AsciiChem::Smiles. SMILES is an input syntax for the one semantic
//! model: parsed molecules are ordinary `Node::Molecule`s with
//! explicit bonds and ring closures, renderable by every formatter.
//!
//! The writer is deterministic (DFS tree from the first atom,
//! single/aromatic continuations, order-independent token tie-breaks,
//! ring digits assigned in encounter order). Aromatic rings perceived
//! from Kekulé input are written in the canonical lowercase form.
//! Stereo carries through: `@`/`@@` tetrahedral parity is corrected
//! against the written neighbour order, and E/Z bond directions ride
//! the directional single bonds.

use std::collections::{BTreeMap, HashMap, HashSet};

use crate::model::{Atom, Molecule, Node};
use crate::perception::{aromatic_rings, hydrogen_valence_sum, implicit_hydrogens};
use crate::structure::{build_graph, linearize, Edge};
use crate::ParseError;

const ORGANIC: [&str; 10] = ["Cl", "Br", "B", "C", "N", "O", "P", "S", "F", "I"];
const AROMATIC_TWO: [&str; 2] = ["se", "as"];
const AROMATIC_ONE: [char; 6] = ['b', 'c', 'n', 'o', 'p', 's'];

type Result<T> = std::result::Result<T, ParseError>;

/// Parses a SMILES string into one `Node::Molecule` per
/// dot-disconnected component.
pub fn parse_smiles(smiles: &str) -> Result<Vec<Node>> {
    SmilesParser::new(smiles).parse()
}

struct SmilesParser<'a> {
    source: &'a str,
    chars: Vec<char>,
    pos: usize,
    atoms: Vec<Atom>,
    adjacency: Vec<BTreeMap<usize, String>>,
    open_rings: HashMap<String, usize>,
    parent: Option<usize>,
    pending_kind: Option<String>,
    pending_direction: Option<String>,
    edge_directions: HashMap<(usize, usize), String>,
}

impl<'a> SmilesParser<'a> {
    fn new(source: &'a str) -> Self {
        SmilesParser {
            source,
            chars: source.chars().collect(),
            pos: 0,
            atoms: Vec::new(),
            adjacency: Vec::new(),
            open_rings: HashMap::new(),
            parent: None,
            pending_kind: None,
            pending_direction: None,
            edge_directions: HashMap::new(),
        }
    }

    fn parse(mut self) -> Result<Vec<Node>> {
        let mut molecules = vec![self.component()?];
        while self.eat('.') {
            molecules.push(self.component()?);
        }
        if self.pos != self.chars.len() {
            return Err(ParseError {
                message: format!(
                    "unexpected character {:?} at position {} in {:?}",
                    self.peek(),
                    self.pos,
                    self.source
                ),
            });
        }
        Ok(molecules)
    }

    fn component(&mut self) -> Result<Node> {
        self.atoms = Vec::new();
        self.adjacency = Vec::new();
        self.open_rings = HashMap::new();
        self.parent = None;
        self.pending_kind = None;
        self.pending_direction = None;

        self.chain()?;

        if !self.open_rings.is_empty() {
            let mut digits: Vec<&String> = self.open_rings.keys().collect();
            digits.sort();
            let listed = digits
                .into_iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                .join(", ");
            return Err(ParseError {
                message: format!("unclosed ring bond digit(s): {listed} in {:?}", self.source),
            });
        }

        let mut edges = Vec::new();
        for (index, neighbours) in self.adjacency.iter().enumerate() {
            for (to, kind) in neighbours {
                if *to > index {
                    let direction = self
                        .edge_directions
                        .get(&(index.min(*to), index.max(*to)))
                        .cloned();
                    edges.push(Edge {
                        from: index,
                        to: *to,
                        kind: kind.clone(),
                        direction,
                    });
                }
            }
        }
        let atoms = std::mem::take(&mut self.atoms);
        Ok(Node::Molecule(Molecule {
            nodes: linearize(atoms, &edges)?,
            ..Molecule::default()
        }))
    }

    fn chain(&mut self) -> Result<()> {
        self.atom_token()?;
        while !self.eoc() {
            if self.peek() == Some('(') {
                self.branch(None, None)?;
            } else if self.bond_start() {
                let (kind, direction) = self.bond_token()?;
                if self.peek() == Some('(') {
                    self.branch(Some(kind), direction)?;
                } else if self.digit_start() {
                    self.ringbond(Some(kind), direction)?;
                } else {
                    self.pending_kind = Some(kind);
                    self.pending_direction = direction;
                    self.atom_token()?;
                }
            } else if self.digit_start() {
                self.ringbond(None, None)?;
            } else {
                self.pending_kind = None;
                self.pending_direction = None;
                self.atom_token()?;
            }
        }
        Ok(())
    }

    fn branch(&mut self, kind: Option<String>, direction: Option<String>) -> Result<()> {
        let mut kind = kind;
        let mut direction = direction;
        self.expect('(')?;
        if kind.is_none() && self.bond_start() {
            let (next_kind, next_direction) = self.bond_token()?;
            kind = Some(next_kind);
            direction = next_direction;
        }
        let saved_parent = self.parent;
        let saved_pending = self.pending_kind.clone();
        let saved_direction = self.pending_direction.clone();
        self.pending_kind = kind;
        self.pending_direction = direction;
        self.chain()?;
        self.expect(')')?;
        self.parent = saved_parent;
        self.pending_kind = saved_pending;
        self.pending_direction = saved_direction;
        Ok(())
    }

    fn ringbond(&mut self, kind: Option<String>, direction: Option<String>) -> Result<()> {
        if direction.is_some() {
            return Err(ParseError {
                message: "directional ring closures are not representable in the model's ring-closure form".into(),
            });
        }
        let digit = self.ring_digit()?;
        let parent = self.parent.expect("ring bond after an atom");
        let opener = self.open_rings.remove(&digit);
        if let Some(opener) = opener {
            if let Some(kind) = &kind {
                if *kind != self.default_kind(opener, parent) {
                    return Err(ParseError {
                        message: format!(
                            "bonded ring closures ({kind}) are not representable in the model's ring-closure form"
                        ),
                    });
                }
            }
            self.add_edge(opener, parent, kind);
        } else {
            self.open_rings.insert(digit, parent);
        }
        Ok(())
    }

    fn atom_token(&mut self) -> Result<()> {
        let atom = self.read_atom()?;
        let index = self.atoms.len();
        self.atoms.push(atom);
        self.adjacency.push(BTreeMap::new());
        if let Some(parent) = self.parent {
            self.add_edge(parent, index, self.pending_kind.clone());
            if let Some(direction) = &self.pending_direction {
                self.edge_directions
                    .insert((parent.min(index), parent.max(index)), direction.clone());
            }
        }
        self.parent = Some(index);
        self.pending_kind = None;
        self.pending_direction = None;
        Ok(())
    }

    fn add_edge(&mut self, from: usize, to: usize, explicit: Option<String>) {
        let kind = explicit.unwrap_or_else(|| self.default_kind(from, to));
        self.adjacency[from].insert(to, kind.clone());
        self.adjacency[to].insert(from, kind);
    }

    fn default_kind(&self, from: usize, to: usize) -> String {
        if self.atoms[from].aromatic == Some(true) && self.atoms[to].aromatic == Some(true) {
            "aromatic".to_string()
        } else {
            "single".to_string()
        }
    }

    // -- tokens ---------------------------------------------------------

    fn read_atom(&mut self) -> Result<Atom> {
        if self.peek() == Some('[') {
            return self.bracket_atom();
        }

        let two: String = self.chars[self.pos..(self.pos + 2).min(self.chars.len())]
            .iter()
            .collect();
        let organic = ORGANIC.iter().find(|s| two.starts_with(*s));
        if let Some(organic) = organic {
            self.pos += organic.chars().count();
            return Ok(Atom {
                element: (*organic).to_string(),
                ..Atom::default()
            });
        }
        let aromatic_two = AROMATIC_TWO.iter().find(|s| two.starts_with(**s));
        let symbol: Option<String> = match aromatic_two {
            Some(s) => Some((*s).to_string()),
            None => match self.peek() {
                Some(ch) if AROMATIC_ONE.contains(&ch) => Some(ch.to_string()),
                _ => None,
            },
        };
        if let Some(symbol) = symbol {
            self.pos += symbol.chars().count();
            return Ok(Atom {
                element: capitalize(&symbol),
                aromatic: Some(true),
                ..Atom::default()
            });
        }
        Err(ParseError {
            message: format!(
                "unexpected character {:?} at position {} in {:?}",
                self.peek(),
                self.pos,
                self.source
            ),
        })
    }

    fn bracket_atom(&mut self) -> Result<Atom> {
        self.expect('[')?;
        let isotope = self.digits();
        let symbol = self.bracket_symbol()?;
        let chirality = self.read_chirality();
        let hydrogens = self.hcount();
        let charge = self.bracket_charge();
        self.skip_class();
        self.expect(']')?;
        Ok(Atom {
            element: capitalize(&symbol),
            isotope,
            charge,
            chirality,
            aromatic: symbol
                .chars()
                .next()
                .is_some_and(|c| c.is_ascii_lowercase())
                .then_some(true),
            hydrogens,
            ..Atom::default()
        })
    }

    fn bracket_symbol(&mut self) -> Result<String> {
        let two: String = self.chars[self.pos..(self.pos + 2).min(self.chars.len())]
            .iter()
            .collect();
        if AROMATIC_TWO.contains(&two.as_str()) {
            self.pos += 2;
            return Ok(two);
        }
        match self.peek() {
            Some(ch) if AROMATIC_ONE.contains(&ch) => {
                self.pos += 1;
                Ok(ch.to_string())
            }
            Some(ch) if ch.is_ascii_uppercase() => {
                let mut symbol = ch.to_string();
                if let Some(next) = self.chars.get(self.pos + 1) {
                    if next.is_ascii_lowercase() {
                        symbol.push(*next);
                    }
                }
                self.pos += symbol.chars().count();
                Ok(symbol)
            }
            other => Err(ParseError {
                message: format!(
                    "expected an element symbol at position {} in {:?} (found {other:?})",
                    self.pos, self.source
                ),
            }),
        }
    }

    // `@` (TH1) / `@@` (TH2): SMILES ordered-neighbour tetrahedral
    // parity, stored verbatim on the atom.
    fn read_chirality(&mut self) -> Option<String> {
        if self.peek() != Some('@') {
            return None;
        }
        let double = self.chars.get(self.pos + 1) == Some(&'@');
        if double {
            self.pos += 2;
            return Some("TH2".to_string());
        }
        self.pos += 1;
        Some("TH1".to_string())
    }

    fn hcount(&mut self) -> Option<u32> {
        if self.peek() != Some('H') {
            return None;
        }
        self.pos += 1;
        let count = self.digits();
        Some(count.and_then(|c| c.parse::<u32>().ok()).unwrap_or(1))
    }

    // "+" | "++" | "+n" | "-" | "--" | "-n" → number-then-sign.
    fn bracket_charge(&mut self) -> Option<String> {
        let sign = self.peek()?;
        if sign != '+' && sign != '-' {
            return None;
        }
        self.pos += 1;
        let second = self.chars.get(self.pos).copied();
        let count = if second == Some(sign) {
            self.pos += 1;
            Some(2)
        } else if second.is_some_and(|c| c.is_ascii_digit()) {
            self.digits().and_then(|d| d.parse::<u32>().ok())
        } else {
            None
        };
        match count {
            Some(count) if count > 0 => Some(format!("{count}{sign}")),
            _ => Some(sign.to_string()),
        }
    }

    fn skip_class(&mut self) {
        if self.peek() != Some(':') {
            return;
        }
        self.pos += 1;
        self.digits();
    }

    // Returns (kind, direction): `/` and `\` are DIRECTIONAL SINGLE
    // bonds (SMILES E/Z stereo).
    fn bond_token(&mut self) -> Result<(String, Option<String>)> {
        let ch = self.peek();
        match ch {
            Some('/') | Some('\\') => {
                self.pos += 1;
                Ok((
                    "single".to_string(),
                    Some(if ch == Some('/') { "up" } else { "down" }.to_string()),
                ))
            }
            Some('-') => {
                self.pos += 1;
                Ok(("single".to_string(), None))
            }
            Some('=') => {
                self.pos += 1;
                Ok(("double".to_string(), None))
            }
            Some('#') => {
                self.pos += 1;
                Ok(("triple".to_string(), None))
            }
            Some('$') => {
                self.pos += 1;
                Ok(("quadruple".to_string(), None))
            }
            Some(':') => {
                self.pos += 1;
                Ok(("aromatic".to_string(), None))
            }
            other => Err(ParseError {
                message: format!(
                    "expected a bond or atom at position {} in {:?} (found {other:?})",
                    self.pos, self.source
                ),
            }),
        }
    }

    fn ring_digit(&mut self) -> Result<String> {
        if self.peek() == Some('%') {
            self.pos += 1;
            let digit: String = self.chars[self.pos..(self.pos + 2).min(self.chars.len())]
                .iter()
                .collect();
            if digit.chars().count() != 2 || !digit.chars().all(|c| c.is_ascii_digit()) {
                return Err(ParseError {
                    message: "malformed %nn ring closure".into(),
                });
            }
            self.pos += 2;
            return Ok(digit);
        }
        match self.peek() {
            Some(ch) if ch.is_ascii_digit() => {
                self.pos += 1;
                Ok(ch.to_string())
            }
            other => Err(ParseError {
                message: format!(
                    "expected ring digit at position {} (found {other:?})",
                    self.pos
                ),
            }),
        }
    }

    // -- character helpers ------------------------------------------------

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn digits(&mut self) -> Option<String> {
        let start = self.pos;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.pos += 1;
        }
        (self.pos != start).then(|| self.chars[start..self.pos].iter().collect())
    }

    fn bond_start(&self) -> bool {
        matches!(
            self.peek(),
            Some('-') | Some('=') | Some('#') | Some('$') | Some(':') | Some('/') | Some('\\')
        )
    }

    fn digit_start(&self) -> bool {
        match self.peek() {
            Some(ch) => ch.is_ascii_digit() || ch == '%',
            None => false,
        }
    }

    fn eoc(&self) -> bool {
        matches!(self.peek(), None | Some('.') | Some(')'))
    }

    fn eat(&mut self, token: char) -> bool {
        if self.peek() == Some(token) {
            self.pos += 1;
            return true;
        }
        false
    }

    fn expect(&mut self, token: char) -> Result<()> {
        if self.peek() != Some(token) {
            return Err(ParseError {
                message: format!(
                    "expected {token:?} at position {} in {:?}",
                    self.pos, self.source
                ),
            });
        }
        self.pos += 1;
        Ok(())
    }
}

fn capitalize(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
        None => String::new(),
    }
}

// -- deterministic writer -----------------------------------------------

const LOWERCASE_AROMATIC: [&str; 8] = ["C", "N", "O", "S", "P", "B", "Se", "As"];
const WRITER_ORGANIC: [&str; 10] = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];
const CONTINUATION_BONDS: [&str; 2] = ["single", "aromatic"];
const HETERO_IMPLICIT_H: [&str; 6] = ["N", "O", "S", "P", "Se", "As"];

fn bond_token_str(kind: &str) -> Option<&'static str> {
    match kind {
        "single" => Some("-"),
        "double" => Some("="),
        "triple" => Some("#"),
        "quadruple" => Some("$"),
        "aromatic" => Some(":"),
        _ => None,
    }
}

/// Emits deterministic SMILES for a formula's molecules, joined with
/// `.` (or for a single molecule, just that molecule).
pub fn write_smiles(nodes: &[Node]) -> Result<String> {
    let molecules: Vec<&Molecule> = nodes
        .iter()
        .filter_map(|node| match node {
            Node::Molecule(mol) => Some(mol),
            _ => None,
        })
        .collect();
    let written = molecules
        .iter()
        .map(|mol| write_molecule(mol))
        .collect::<Result<Vec<_>>>()?;
    Ok(written.join("."))
}

fn write_molecule(molecule: &Molecule) -> Result<String> {
    let (atoms, edges) = build_graph(molecule);
    // A single-atom component is valid SMILES; multiple unbonded atoms
    // are a formula, not a structure.
    if edges.is_empty() && atoms.len() > 1 {
        return Err(ParseError {
            message: "molecule has no bonds — a formula is not a structure".into(),
        });
    }

    let mut adjacency: Vec<BTreeMap<usize, String>> = vec![BTreeMap::new(); atoms.len()];
    let mut directions: HashMap<(usize, usize), String> = HashMap::new();
    let mut model_neighbours: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
    for edge in &edges {
        adjacency[edge.from].insert(edge.to, edge.kind.clone());
        adjacency[edge.to].insert(edge.from, edge.kind.clone());
        if let Some(direction) = &edge.direction {
            directions.insert(edge.ordered(), direction.clone());
        }
        model_neighbours[edge.from].push(edge.to);
        model_neighbours[edge.to].push(edge.from);
    }

    // The canonicalizer consumes the perception pass: atoms in
    // perceived aromatic rings are written lowercase and the ring
    // bonds between them carry no bond symbol, so Kekulé benzene
    // canonicalizes to `c1ccccc1`.
    let mut perceived_aromatic: HashSet<usize> = HashSet::new();
    let mut perceived_ring_edges: HashSet<(usize, usize)> = HashSet::new();
    let mut perceived_hydrogens: HashMap<usize, Option<u32>> = HashMap::new();
    for ring in aromatic_rings(&atoms, &edges) {
        for atom in &ring {
            perceived_aromatic.insert(*atom);
            // Only rings perceived from Kekulé bond orders get an
            // inferred count; source-aromatic atoms already state
            // hydrogens the way SMILES does ([nH]), and aromatic 1.5
            // bond orders do not map onto heteroatom valences.
            if atoms[*atom].aromatic != Some(true) {
                perceived_hydrogens.insert(
                    *atom,
                    implicit_hydrogens(&atoms[*atom], hydrogen_valence_sum(*atom, &edges)),
                );
            }
        }
        for (i, atom) in ring.iter().enumerate() {
            let other = ring[(i + 1) % ring.len()];
            perceived_ring_edges.insert((*atom.min(&other), *atom.max(&other)));
        }
    }

    let atom_count = atoms.len();
    let mut writer = DeterministicWriter {
        atoms,
        adjacency,
        directions,
        model_neighbours,
        perceived_aromatic,
        perceived_ring_edges,
        perceived_hydrogens,
        visited: HashSet::new(),
        digits: vec![Vec::new(); atom_count],
        digit_by_edge: HashMap::new(),
        next_digit: 1,
        tree_children: vec![Vec::new(); atom_count],
        closures: vec![Vec::new(); atom_count],
        tree_size: vec![1; atom_count],
        written_order: vec![Vec::new(); atom_count],
    };
    writer.write()
}

#[derive(Clone)]
enum Unit {
    Literal(&'static str),
    Atom(usize),
}

/// A tetrahedral reference neighbour: the stated implicit hydrogen
/// ("h") or a bonded neighbour by atom index.
#[derive(Clone, PartialEq, Eq)]
enum Item {
    H,
    N(usize),
}

struct DeterministicWriter {
    atoms: Vec<Atom>,
    adjacency: Vec<BTreeMap<usize, String>>,
    directions: HashMap<(usize, usize), String>,
    model_neighbours: Vec<Vec<usize>>,
    perceived_aromatic: HashSet<usize>,
    perceived_ring_edges: HashSet<(usize, usize)>,
    perceived_hydrogens: HashMap<usize, Option<u32>>,

    visited: HashSet<usize>,
    digits: Vec<Vec<String>>,
    digit_by_edge: HashMap<(usize, usize), String>,
    next_digit: u32,
    tree_children: Vec<Vec<usize>>,
    closures: Vec<Vec<usize>>,
    tree_size: Vec<usize>,
    written_order: Vec<Vec<usize>>,
}

impl DeterministicWriter {
    fn write(&mut self) -> Result<String> {
        if self.atoms.is_empty() {
            return Ok(String::new());
        }
        let count = self.atoms.len();
        self.digits = vec![Vec::new(); count];
        self.build_tree(0, None);
        self.compute_tree_size(0);
        self.compute_written_orders(0, None);
        let units = self.emit(0, None)?;
        let mut out = String::new();
        for unit in units {
            match unit {
                Unit::Literal(text) => out.push_str(text),
                Unit::Atom(index) => {
                    out.push_str(&self.atom_token(index));
                    for digit in &self.digits[index] {
                        out.push_str(digit);
                    }
                }
            }
        }
        Ok(out)
    }

    // Deterministic DFS tree: children in index order; edges to
    // already-visited atoms become ring digits. Recording the tree
    // first means a branch can never walk around a ring and swallow
    // the continuation.
    fn build_tree(&mut self, index: usize, parent: Option<usize>) {
        self.visited.insert(index);
        let neighbours: Vec<usize> = self.adjacency[index].keys().copied().collect();
        let mut children = Vec::new();
        let mut closures = Vec::new();
        for nb in neighbours {
            if Some(nb) != parent && self.visited.contains(&nb) {
                closures.push(nb);
            } else if !self.visited.contains(&nb) {
                children.push(nb);
                // Recurse inline: the subtree marks its atoms visited,
                // so a later sibling meeting them closes a ring at this
                // atom rather than growing a second tree.
                self.build_tree(nb, Some(index));
            }
        }
        self.tree_children[index] = children;
        self.closures[index] = closures;
    }

    fn compute_tree_size(&mut self, index: usize) -> usize {
        let mut size = 1;
        for child in self.tree_children[index].clone() {
            size += self.compute_tree_size(child);
        }
        self.tree_size[index] = size;
        size
    }

    // The order neighbours are WRITTEN in: parent first, then ring
    // closures, then branches (rank order), then the continuation.
    // Tetrahedral parity flips with odd permutations of this order
    // against the model's reference order, so the emitted @/@@ is
    // corrected per atom (acyclic-correct; ring closures in the
    // reference order are a later slice).
    fn compute_written_orders(&mut self, index: usize, parent: Option<usize>) {
        let mut order: Vec<usize> = parent.into_iter().collect();
        order.extend(self.closures[index].iter().copied());
        let mut ordered = self.tree_children[index].clone();
        ordered.sort_by(|a, b| self.compare_continuation(index, *a, *b));
        if !ordered.is_empty() {
            let continuation = ordered[0];
            order.extend_from_slice(&ordered[1..]);
            order.push(continuation);
        }
        self.written_order[index] = order;
        for child in self.tree_children[index].clone() {
            self.compute_written_orders(child, Some(index));
        }
    }

    // Correction token for an atom's stated parity: how many
    // transpositions separate the written neighbour order from the
    // model's reference order.
    fn chirality_token(&self, index: usize) -> String {
        let Some(chirality) = &self.atoms[index].chirality else {
            return String::new();
        };
        // The implicit hydrogen participates as a neighbour ("h"),
        // when the atom states one (JS truthiness: a stated 0 does not
        // join the reference order).
        let stated_h = matches!(self.atoms[index].hydrogens, Some(count) if count > 0);
        let mut reference: Vec<Item> = Vec::new();
        if stated_h {
            reference.push(Item::H);
        }
        reference.extend(self.model_neighbours[index].iter().map(|n| Item::N(*n)));
        let mut written: Vec<Item> = Vec::new();
        if stated_h {
            written.push(Item::H);
        }
        written.extend(self.written_order[index].iter().map(|n| Item::N(*n)));
        let mut inversions = 0;
        for (i, item) in written.iter().enumerate() {
            for later in &written[i + 1..] {
                let item_at = reference.iter().position(|candidate| candidate == item);
                let later_at = reference.iter().position(|candidate| candidate == later);
                if let (Some(a), Some(b)) = (item_at, later_at) {
                    if a > b {
                        inversions += 1;
                    }
                }
            }
        }
        let token = if chirality == "TH2" || inversions % 2 == 1 {
            "@@"
        } else {
            "@"
        };
        token.to_string()
    }

    fn emit(&mut self, index: usize, parent: Option<usize>) -> Result<Vec<Unit>> {
        let mut units = Vec::new();
        let incoming_kind = match parent {
            Some(parent) => Some(
                self.adjacency[parent]
                    .get(&index)
                    .cloned()
                    .expect("tree edge exists"),
            ),
            None => None,
        };
        if let Some(kind) = &incoming_kind {
            let parent = parent.expect("kind implies a parent");
            let token = self.bond_symbol(kind, parent, index)?;
            units.push(Unit::Literal(token));
        }
        units.push(Unit::Atom(index));
        // Ring digits attach to per-atom lists (both endpoints); no
        // literal units — they would render twice.
        let closures = self.closures[index].clone();
        for nb in closures {
            self.ring_digit_for(index, nb)?;
        }

        let children = self.tree_children[index].clone();
        if children.is_empty() {
            return Ok(units);
        }

        let mut ordered = children;
        ordered.sort_by(|a, b| self.compare_continuation(index, *a, *b));
        let continuation = ordered[0];
        let branches = &ordered[1..];

        for child in branches {
            units.push(Unit::Literal("("));
            units.extend(self.emit(*child, Some(index))?);
            units.push(Unit::Literal(")"));
        }
        units.extend(self.emit(continuation, Some(index))?);
        Ok(units)
    }

    // A closure-bearing child continues the chain (SMILES writes
    // ring-closure digits at the end of a traversal — and this keeps
    // the written atom order equal to the model's), then
    // single/aromatic bonds, then larger subtrees, then the
    // lexicographically smallest atom token — order-independent, so
    // the canonical form does not depend on how the input numbered its
    // atoms.
    fn compare_continuation(&self, index: usize, a: usize, b: usize) -> std::cmp::Ordering {
        let rank = |child: usize| -> (u8, u8, std::cmp::Reverse<usize>, String) {
            let kind = &self.adjacency[index][&child];
            (
                u8::from(self.closures[child].is_empty()),
                u8::from(!CONTINUATION_BONDS.contains(&kind.as_str())),
                // Larger subtrees rank first (the chain follows the
                // main skeleton; side chains branch off).
                std::cmp::Reverse(self.tree_size[child]),
                self.atom_token(child),
            )
        };
        rank(a).cmp(&rank(b))
    }

    fn ring_digit_for(&mut self, a: usize, b: usize) -> Result<String> {
        let key = (a.min(b), a.max(b));
        if let Some(existing) = self.digit_by_edge.get(&key) {
            return Ok(existing.clone());
        }
        if self.next_digit > 9 {
            return Err(ParseError {
                message: "too many ring closures for SMILES output".into(),
            });
        }
        let digit = self.next_digit.to_string();
        self.next_digit += 1;
        self.digit_by_edge.insert(key, digit.clone());
        self.digits[key.0].push(digit.clone());
        self.digits[key.1].push(digit.clone());
        Ok(digit)
    }

    fn bond_symbol(&self, kind: &str, from: usize, to: usize) -> Result<&'static str> {
        // A bond inside a perceived aromatic ring is implied by the
        // lowercase atoms: writing `=`/`-` there would state Kekulé
        // alternation the canonical form does not carry.
        let key = (from.min(to), from.max(to));
        if self.perceived_ring_edges.contains(&key) && self.aromatic(from) && self.aromatic(to) {
            return Ok("");
        }
        if kind == "single" {
            if let Some(direction) = self.directions.get(&key) {
                return Ok(if direction == "up" { "/" } else { "\\" });
            }
            return Ok(if self.aromatic(from) && self.aromatic(to) {
                "-"
            } else {
                ""
            });
        }
        let token = bond_token_str(kind).ok_or_else(|| ParseError {
            message: format!("{kind} bonds have no SMILES form (v1 subset)"),
        })?;
        Ok(if self.aromatic(from) && self.aromatic(to) {
            ""
        } else {
            token
        })
    }

    fn atom_token(&self, index: usize) -> String {
        let atom = &self.atoms[index];
        let aromatic = self.aromatic(index);
        let lowercase = aromatic && LOWERCASE_AROMATIC.contains(&atom.element.as_str());
        let hyd: Option<Option<u32>> = if atom.hydrogens.is_some() {
            Some(atom.hydrogens)
        } else {
            self.perceived_hydrogens.get(&index).copied()
        };
        // Bare `n`/`o`/`s` imply no hydrogen: a perceived pyrrole N
        // must be written [nH].
        let hetero_needs_h = aromatic
            && lowercase
            && matches!(hyd, Some(Some(count)) if count > 0)
            && HETERO_IMPLICIT_H.contains(&atom.element.as_str());
        let needs_bracket = atom.charge.is_some()
            || atom.isotope.is_some()
            || atom.hydrogens.is_some()
            || atom.chirality.is_some()
            || hetero_needs_h
            || (!WRITER_ORGANIC.contains(&atom.element.as_str()) && !lowercase)
            || (aromatic && !lowercase);
        let symbol = if lowercase {
            atom.element.to_lowercase()
        } else {
            atom.element.clone()
        };
        if !needs_bracket {
            return symbol;
        }

        let mut token = String::from("[");
        if let Some(isotope) = &atom.isotope {
            token.push_str(isotope);
        }
        token.push_str(&symbol);
        token.push_str(&self.chirality_token(index));
        if atom.hydrogens.is_some() || hetero_needs_h {
            token.push('H');
            let count = hyd.flatten();
            if count != Some(1) {
                token.push_str(&count.map(|c| c.to_string()).unwrap_or_default());
            }
        }
        if let Some(charge) = &atom.charge {
            token.push_str(&charge_suffix(charge));
        }
        token.push(']');
        token
    }

    fn aromatic(&self, index: usize) -> bool {
        self.atoms[index].aromatic == Some(true) || self.perceived_aromatic.contains(&index)
    }
}

// The model's number-then-sign charge ("2+") → SMILES ("+2").
fn charge_suffix(charge: &str) -> String {
    let count: u32 = charge
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(1);
    let sign = charge.chars().last().unwrap_or('+');
    if count > 1 {
        format!("{sign}{count}")
    } else {
        sign.to_string()
    }
}
