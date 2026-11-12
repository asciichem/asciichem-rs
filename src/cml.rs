//! CML (Chemical Markup Language) round-trip — TODO.v2/17 slice 7.
//!
//! Bidirectional conversion between the AsciiChem model and the CML
//! XML wire, mirroring the reference implementation's pipeline:
//!
//! ```text
//! model::Node  ──emit──►  canonical wire (molecule/atom/bond/
//!                         reaction + the aci: extension channels)
//!              ──DOM──►   libleptris  ──serialize──►  XML string
//!
//! XML string  ──libleptris SAX──►  namespace-resolved tree
//!              ──read──►  canonical wire  ──restore──►  model::Node
//! ```
//!
//! # The leptris directive
//!
//! Every byte of XML parsing and serialization goes through the
//! leptris engine (the mandated serialization layer for the family —
//! never quick-xml, never hand-rolled strings). The safe `leptris`
//! crate carries the parse half (DOM, SAX, XPath); the engine's
//! public C builder API ([`ffi`]) carries the construction half, the
//! same face the Ruby binding exposes. The engine is linked at build
//! time via `LEPTRIS_LIB_PATH` (see `leptris`'s `build.rs`).
//!
//! # The aci: extension namespace
//!
//! CML's standard wire covers element, isotope, charge, count, and
//! bond order — not oxidation state, Lewis markers, electron
//! configuration, embedded math, free-form text, or group structure.
//! AsciiChem carries these through an `aci:` (AsciiChem extension)
//! namespace, three channels one per scope, each a frozen registry:
//!
//! - **Atom attributes** — `aci:<wire_name>` attributes on `<atom>`
//!   elements (`oxidationState`, `lonePairs`, `radicalElectrons`,
//!   `ringClosures`).
//! - **Top-level elements** — `<aci:<name> position="N">text</…>`
//!   children of `<cml>` for `ElectronConfiguration`,
//!   `EmbeddedMath`, `Text`, `Mechanism`; `position` preserves node
//!   ordering in the formula.
//! - **Molecule groups** — `<aci:group multiplicity="N" bracket="X"
//!   atomRefs="a1 a2"/>` children of `<molecule>`, preserving AsciiChem
//!   group structure through the canonical flattening.
//!
//! The `xmlns:aci` declaration appears only when at least one
//! extension is present; plain atoms and molecules produce clean CML.
//! Unknown elements and attributes are ignored on read; foreign
//! namespaces can never be mistaken for a channel (classification is
//! by resolved namespace URI, not by written prefix).
//!
//! Constructs with no wire representation and no channel
//! (`Crystal`, `Spectrum`, `ZMatrix`, `Calculation`, reaction
//! conditions, molecule annotations) report [`CmlError::Unsupported`]
//! on emit rather than dropping silently — the remaining surface is
//! tracked in asciichem-rs#15.

use std::ffi::{CStr, CString};
use std::fmt;
use std::os::raw::{c_char, c_int};

use crate::model::{leading_integer, Atom, Bond, Group, Molecule, Node, Reaction, ReactionCascade};

/// The CML core namespace (the standard wire).
pub const CML_NAMESPACE: &str = "http://www.xml-cml.org/schema";
/// The AsciiChem extension namespace (the `aci:` channels).
pub const ACI_NAMESPACE: &str = "https://asciichem.org/cml-ext";
/// The written prefix for [`ACI_NAMESPACE`] on emit.
pub const ACI_PREFIX: &str = "aci";

/// An error surfaced by CML emit or read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CmlError {
    /// The XML engine rejected the input, or serialization failed.
    Xml(String),
    /// The XML is well-formed but not the wire this module reads.
    Invalid(String),
    /// The model construct has no wire representation in this slice.
    Unsupported(String),
}

impl fmt::Display for CmlError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CmlError::Xml(message) => write!(f, "CML XML engine error: {message}"),
            CmlError::Invalid(message) => write!(f, "CML wire violation: {message}"),
            CmlError::Unsupported(message) => {
                write!(f, "construct has no CML wire representation: {message}")
            }
        }
    }
}

impl std::error::Error for CmlError {}

// ---------------------------------------------------------------------------
// The engine face: construction + serialization over the linked
// libleptris. The parse half rides the `leptris` crate's safe API
// (`sax::parse`); these declarations mirror the engine's public
// headers (`src/include/leptris.h`) exactly — the same functions the
// Ruby binding calls. Elements are document-owned; freeing the
// document frees every element handle created against it.
// ---------------------------------------------------------------------------

#[allow(non_camel_case_types)]
mod ffi {
    use std::os::raw::{c_char, c_int, c_void};

    #[repr(C)]
    pub struct LeptrisDocument {
        _private: [u8; 0],
    }
    #[repr(C)]
    pub struct LeptrisElement {
        _private: [u8; 0],
    }

    extern "C" {
        pub fn leptris_document_create() -> *mut LeptrisDocument;
        pub fn leptris_document_free(doc: *mut LeptrisDocument);
        pub fn leptris_document_set_root(
            doc: *mut LeptrisDocument,
            elem: *mut LeptrisElement,
        ) -> c_int;
        pub fn leptris_document_serialize(
            doc: *mut LeptrisDocument,
            options: *mut c_void,
        ) -> *mut c_char;
        pub fn leptris_element_create(
            doc: *mut LeptrisDocument,
            name: *const c_char,
        ) -> *mut LeptrisElement;
        pub fn leptris_element_create_child(
            parent: *mut LeptrisElement,
            name: *const c_char,
        ) -> *mut LeptrisElement;
        pub fn leptris_element_append_child(
            parent: *mut LeptrisElement,
            child: *mut LeptrisElement,
        ) -> c_int;
        pub fn leptris_element_set_attribute(
            elem: *mut LeptrisElement,
            name: *const c_char,
            value: *const c_char,
        ) -> c_int;
        pub fn leptris_element_set_text(elem: *mut LeptrisElement, text: *const c_char) -> c_int;
        pub fn leptris_free_string(s: *mut c_char);
    }
}

/// An owned engine document — the emit side's DOM root. Freed on drop.
struct WireDom {
    raw: *mut ffi::LeptrisDocument,
}

impl WireDom {
    fn new() -> Result<WireDom, CmlError> {
        let raw = unsafe { ffi::leptris_document_create() };
        if raw.is_null() {
            return Err(CmlError::Xml("document creation failed".into()));
        }
        Ok(WireDom { raw })
    }

    fn owned_string(ptr: *mut c_char) -> Result<String, CmlError> {
        if ptr.is_null() {
            return Err(CmlError::Xml("serialization failed".into()));
        }
        let s = unsafe {
            let s = CStr::from_ptr(ptr).to_string_lossy().into_owned();
            ffi::leptris_free_string(ptr);
            s
        };
        Ok(s)
    }

    fn create_element(&self, name: &str) -> Result<*mut ffi::LeptrisElement, CmlError> {
        let name = cstr(name)?;
        let raw = unsafe { ffi::leptris_element_create(self.raw, name.as_ptr()) };
        if raw.is_null() {
            return Err(CmlError::Xml(format!("element {name:?} creation failed")));
        }
        Ok(raw)
    }
}

impl Drop for WireDom {
    fn drop(&mut self) {
        unsafe { ffi::leptris_document_free(self.raw) };
    }
}

fn cstr(value: &str) -> Result<CString, CmlError> {
    CString::new(value).map_err(|_| CmlError::Invalid("NUL byte in wire string".into()))
}

fn check(status: c_int, what: &str) -> Result<(), CmlError> {
    if status == 0 {
        Ok(())
    } else {
        Err(CmlError::Xml(format!("{what} failed (status {status})")))
    }
}

fn set_attribute(elem: *mut ffi::LeptrisElement, name: &str, value: &str) -> Result<(), CmlError> {
    let name = cstr(name)?;
    let value = cstr(value)?;
    unsafe {
        check(
            ffi::leptris_element_set_attribute(elem, name.as_ptr(), value.as_ptr()),
            "attribute set",
        )
    }
}

fn set_text(elem: *mut ffi::LeptrisElement, content: &str) -> Result<(), CmlError> {
    let text = cstr(content)?;
    unsafe {
        check(
            ffi::leptris_element_set_text(elem, text.as_ptr()),
            "text set",
        )
    }
}

fn child_element(
    parent: *mut ffi::LeptrisElement,
    name: &str,
) -> Result<*mut ffi::LeptrisElement, CmlError> {
    let name = cstr(name)?;
    let raw = unsafe { ffi::leptris_element_create_child(parent, name.as_ptr()) };
    if raw.is_null() {
        return Err(CmlError::Xml(format!(
            "child element {name:?} creation failed"
        )));
    }
    Ok(raw)
}

// ---------------------------------------------------------------------------
// Frozen channel registries (the aci: extension surface). Adding a
// channel member is one entry here — no other code changes.
// ---------------------------------------------------------------------------

/// Atom-attribute channel: model field name → wire attribute name
/// (without the `aci:` prefix), in wire order.
const ATOM_FIELDS: [(&str, &str); 4] = [
    ("oxidation_state", "oxidationState"),
    ("lone_pairs", "lonePairs"),
    ("radical_electrons", "radicalElectrons"),
    ("ring_closures", "ringClosures"),
];

/// Top-level channel element names.
const TOP_LEVEL_ELEMENTS: [&str; 4] =
    ["electronConfiguration", "embeddedMath", "text", "mechanism"];

/// Group brackets on the wire (the model stores the same names).
const BRACKETS: [&str; 3] = ["paren", "square", "brace"];

/// Bond kind → CML order code (the reference's CML_ORDER_CODES).
fn order_code(kind: &str) -> &'static str {
    match kind {
        "double" => "D",
        "triple" => "T",
        "quadruple" => "Q",
        "wedge" => "W",
        "hash" => "H",
        "dative" => "A",
        "wavy" => "V",
        _ => "S",
    }
}

/// CML order code → bond kind (the reference's KIND_BY_CML_ORDER).
fn kind_from_order(order: &str) -> &'static str {
    match order {
        "D" => "double",
        "T" => "triple",
        "Q" => "quadruple",
        "W" => "wedge",
        "H" => "hash",
        "A" => "dative",
        "V" => "wavy",
        _ => "single",
    }
}

/// Parse an atom charge ("2+", "+2", "-", "2-") into a signed value.
fn charge_value(charge: &str) -> Option<i64> {
    let bytes = charge.as_bytes();
    let (n, sign) = match (bytes.first(), bytes.last()) {
        (Some(first), Some(last)) if (*last == b'+' || *last == b'-') && first.is_ascii_digit() => {
            (&charge[..charge.len() - 1], *last)
        }
        (Some(first), Some(last))
            if (*first == b'+' || *first == b'-') && last.is_ascii_digit() =>
        {
            (&charge[1..], *first)
        }
        (Some(only), _) if *only == b'+' || *only == b'-' => ("", *only),
        _ => return None,
    };
    if !n.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let magnitude: i64 = if n.is_empty() { 1 } else { n.parse().ok()? };
    Some(if sign == b'+' { magnitude } else { -magnitude })
}

// ---------------------------------------------------------------------------
// Canonical IDs — a/b/m/r/g counters shared across one translation
// pass, so re-running on equivalent input yields byte-equal output.
// ---------------------------------------------------------------------------

#[derive(Default)]
struct Ids {
    atoms: usize,
    bonds: usize,
    molecules: usize,
    reactions: usize,
    groups: usize,
}

impl Ids {
    fn next_atom(&mut self) -> String {
        self.atoms += 1;
        format!("a{}", self.atoms)
    }
    fn next_bond(&mut self) -> String {
        self.bonds += 1;
        format!("b{}", self.bonds)
    }
    fn next_molecule(&mut self) -> String {
        self.molecules += 1;
        format!("m{}", self.molecules)
    }
    fn next_reaction(&mut self) -> String {
        self.reactions += 1;
        format!("r{}", self.reactions)
    }
    fn next_group(&mut self) -> String {
        self.groups += 1;
        format!("g{}", self.groups)
    }
}

// ---------------------------------------------------------------------------
// Emit: model → DOM → XML.
// ---------------------------------------------------------------------------

/// Serializes formula nodes to a CML XML string.
///
/// `parse_cml(&to_cml(&nodes)?).unwrap()` rebuilds the model for the
/// gated surface — the three-way round-trip contract.
pub fn to_cml(nodes: &[Node]) -> Result<String, CmlError> {
    let dom = WireDom::new()?;
    let mut ids = Ids::default();

    // Walk every formula node, per-category, like the reference
    // adapter: molecules, reactions, cascades regroup; top-level aci:
    // entries keep their formula position.
    let mut molecules = Vec::new();
    let mut reactions = Vec::new();
    let mut reaction_lists = Vec::new();
    let mut top_level: Vec<(usize, &'static str, String)> = Vec::new();
    let mut needs_aci = false;

    for (position, node) in nodes.iter().enumerate() {
        match node {
            Node::Molecule(molecule) => {
                let (element, extensions) = emit_molecule(&dom, &mut ids, molecule)?;
                needs_aci |= extensions;
                molecules.push(element);
            }
            Node::Reaction(reaction) => {
                let (element, extensions) = emit_reaction(&dom, &mut ids, reaction)?;
                needs_aci |= extensions;
                reactions.push(element);
            }
            Node::ReactionCascade(cascade) => {
                let (element, extensions) = emit_reaction_list(&dom, &mut ids, cascade)?;
                needs_aci |= extensions;
                reaction_lists.push(element);
            }
            Node::Atom(_) | Node::Bond(_) => {
                return Err(CmlError::Unsupported(format!(
                    "bare {} outside a molecule",
                    node.diagnostic_label()
                )))
            }
            other => {
                let (element, content) = top_level_payload(other)?;
                needs_aci = true;
                top_level.push((position, element, content));
            }
        }
    }

    let root = dom.create_element("cml")?;
    set_attribute(root, "xmlns", CML_NAMESPACE)?;
    if needs_aci {
        set_attribute(root, "xmlns:aci", ACI_NAMESPACE)?;
    }
    for (position, element, content) in &top_level {
        let entry = child_element(root, &format!("{ACI_PREFIX}:{element}"))?;
        set_attribute(entry, "position", &position.to_string())?;
        set_text(entry, content)?;
    }
    for molecule in molecules.drain(..) {
        unsafe {
            check(
                ffi::leptris_element_append_child(root, molecule),
                "molecule attach",
            )?
        };
    }
    for reaction in reactions.drain(..) {
        unsafe {
            check(
                ffi::leptris_element_append_child(root, reaction),
                "reaction attach",
            )?
        };
    }
    for list in reaction_lists.drain(..) {
        unsafe {
            check(
                ffi::leptris_element_append_child(root, list),
                "reactionList attach",
            )?
        };
    }
    unsafe { check(ffi::leptris_document_set_root(dom.raw, root), "root set")? };
    let xml = unsafe {
        WireDom::owned_string(ffi::leptris_document_serialize(
            dom.raw,
            std::ptr::null_mut(),
        ))?
    };
    Ok(xml)
}

/// The top-level channel payload: wire element name + text content.
/// Text-shaped constructs serialize as their canonical text render
/// (round-trip-safe by construction); embedded math carries the raw
/// source (the backticks are re-added on read).
fn top_level_payload(node: &Node) -> Result<(&'static str, String), CmlError> {
    match node {
        Node::ElectronConfiguration(_) => Ok(("electronConfiguration", node.to_string())),
        Node::EmbeddedMath(source) => Ok(("embeddedMath", source.clone())),
        Node::Text(content) => Ok(("text", format!("\"{content}\""))),
        Node::Mechanism(_) => Ok(("mechanism", node.to_string())),
        other => Err(CmlError::Unsupported(other.diagnostic_label())),
    }
}

/// Emits one molecule; the bool reports whether any aci: extension
/// was written (atom attributes or group elements), so the root can
/// declare the namespace only when needed.
fn emit_molecule(
    dom: &WireDom,
    ids: &mut Ids,
    molecule: &Molecule,
) -> Result<(*mut ffi::LeptrisElement, bool), CmlError> {
    let mut walker = MoleculeWalker::new(ids);
    walker.walk_nodes(&molecule.nodes, None);
    walker.emit_ring_bonds();
    let MoleculeWalker {
        ids,
        atoms,
        bonds,
        groups,
        ..
    } = walker;

    let mut extensions = false;
    let element = dom.create_element("molecule")?;
    set_attribute(element, "id", &ids.next_molecule())?;
    if let Some(coefficient) = &molecule.coefficient {
        set_attribute(element, "count", coefficient)?;
    }
    if let Some(total) = total_formal_charge(&atoms) {
        set_attribute(element, "formalCharge", &total)?;
    }
    if !atoms.is_empty() {
        let array = child_element(element, "atomArray")?;
        for atom in &atoms {
            let el = child_element(array, "atom")?;
            set_attribute(el, "id", &atom.id)?;
            set_attribute(el, "elementType", &atom.element)?;
            if let Some(isotope) = &atom.isotope {
                set_attribute(el, "isotope", isotope)?;
            }
            if let Some(charge) = &atom.formal_charge {
                set_attribute(el, "formalCharge", charge)?;
            }
            if let Some(count) = &atom.count {
                set_attribute(el, "count", count)?;
            }
            for (field, wire) in ATOM_FIELDS {
                if let Some(value) = atom.extension(field) {
                    set_attribute(el, &format!("{ACI_PREFIX}:{wire}"), &value)?;
                    extensions = true;
                }
            }
        }
    }
    if !bonds.is_empty() {
        let array = child_element(element, "bondArray")?;
        for bond in &bonds {
            let el = child_element(array, "bond")?;
            set_attribute(el, "id", &bond.id)?;
            set_attribute(el, "atomRefs2", &format!("{} {}", bond.from, bond.to))?;
            set_attribute(el, "order", bond.order)?;
        }
    }
    if !groups.is_empty() {
        extensions = true;
        for group in &groups {
            let el = child_element(element, &format!("{ACI_PREFIX}:group"))?;
            if let Some(multiplicity) = &group.multiplicity {
                set_attribute(el, "multiplicity", multiplicity)?;
            }
            set_attribute(el, "bracket", &group.bracket)?;
            set_attribute(el, "atomRefs", &group.atom_ids.join(" "))?;
        }
    }
    Ok((element, extensions))
}

fn emit_reaction(
    dom: &WireDom,
    ids: &mut Ids,
    reaction: &Reaction,
) -> Result<(*mut ffi::LeptrisElement, bool), CmlError> {
    if reaction
        .conditions
        .as_ref()
        .is_some_and(|conditions| conditions.above.is_some() || conditions.below.is_some())
    {
        // The conditions channel (aci:condition, like the reference's
        // ConditionsExtensions) is not in this slice — report rather
        // than drop (asciichem-rs#15 remaining list).
        return Err(CmlError::Unsupported("reaction conditions".into()));
    }
    let element = dom.create_element("reaction")?;
    set_attribute(element, "id", &ids.next_reaction())?;
    set_attribute(element, "title", &reaction.arrow)?;
    set_attribute(element, "type", &reaction.arrow)?;
    let mut extensions = false;
    let reactant_list = child_element(element, "reactantList")?;
    for molecule in &reaction.reactants {
        extensions |= emit_participant(dom, ids, reactant_list, "reactant", molecule)?;
    }
    let product_list = child_element(element, "productList")?;
    for molecule in &reaction.products {
        extensions |= emit_participant(dom, ids, product_list, "product", molecule)?;
    }
    Ok((element, extensions))
}

fn emit_participant(
    dom: &WireDom,
    ids: &mut Ids,
    list: *mut ffi::LeptrisElement,
    role: &str,
    molecule: &Node,
) -> Result<bool, CmlError> {
    let Node::Molecule(molecule) = molecule else {
        return Err(CmlError::Unsupported(format!(
            "reaction participant {}",
            molecule.diagnostic_label()
        )));
    };
    let participant = child_element(list, role)?;
    let substance = child_element(participant, "substance")?;
    set_attribute(substance, "role", role)?;
    let (molecule_el, extensions) = emit_molecule(dom, ids, molecule)?;
    unsafe {
        check(
            ffi::leptris_element_append_child(substance, molecule_el),
            "molecule attach",
        )?
    };
    Ok(extensions)
}

fn emit_reaction_list(
    dom: &WireDom,
    ids: &mut Ids,
    cascade: &ReactionCascade,
) -> Result<(*mut ffi::LeptrisElement, bool), CmlError> {
    let element = dom.create_element("reactionList")?;
    let mut extensions = false;
    for step in &cascade.steps {
        let (reaction_el, step_extensions) = emit_reaction(dom, ids, step)?;
        extensions |= step_extensions;
        unsafe {
            check(
                ffi::leptris_element_append_child(element, reaction_el),
                "reaction attach",
            )?
        };
    }
    Ok((element, extensions))
}

fn total_formal_charge(atoms: &[WireAtom]) -> Option<String> {
    let mut sum = 0i64;
    let mut charged = false;
    for atom in atoms {
        if let Some(value) = atom.formal_charge.as_deref().and_then(charge_value) {
            charged = true;
            sum += value;
        }
    }
    charged.then(|| {
        if sum > 0 {
            format!("+{sum}")
        } else {
            sum.to_string()
        }
    })
}

/// Walks a molecule's nodes, producing the ordered atom list, the
/// bond list with atomRefs, and the group records — the single source
/// of atom ID assignment. Mirrors the reference's MoleculeWalker.
struct MoleculeWalker<'a> {
    ids: &'a mut Ids,
    atoms: Vec<WireAtom>,
    bonds: Vec<WireBond>,
    groups: Vec<WireGroup>,
    group_stack: Vec<usize>,
    pending_bond: Option<String>,
    last_atom: Option<String>,
    open_rings: Vec<(char, String)>,
    ring_pairs: Vec<(String, String)>,
}

struct WireAtom {
    id: String,
    element: String,
    isotope: Option<String>,
    formal_charge: Option<String>,
    count: Option<String>,
    /// The aci: channel values, keyed by model field name.
    extensions: [Option<String>; ATOM_FIELDS.len()],
}

impl WireAtom {
    fn extension(&self, field: &str) -> Option<String> {
        ATOM_FIELDS
            .iter()
            .position(|(name, _)| *name == field)
            .and_then(|index| self.extensions[index].clone())
    }
}

struct WireBond {
    id: String,
    from: String,
    to: String,
    order: &'static str,
}

struct WireGroup {
    multiplicity: Option<String>,
    bracket: String,
    atom_ids: Vec<String>,
}

impl<'a> MoleculeWalker<'a> {
    fn new(ids: &'a mut Ids) -> MoleculeWalker<'a> {
        MoleculeWalker {
            ids,
            atoms: Vec::new(),
            bonds: Vec::new(),
            groups: Vec::new(),
            group_stack: Vec::new(),
            pending_bond: None,
            last_atom: None,
            open_rings: Vec::new(),
            ring_pairs: Vec::new(),
        }
    }

    fn walk_nodes(&mut self, nodes: &[Node], multiplier: Option<i64>) {
        for node in nodes {
            match node {
                Node::Atom(atom) => self.emit_atom(atom, multiplier),
                Node::Bond(bond) => self.pending_bond = Some(bond.kind.clone()),
                Node::Group(group) => self.enter_group(group, multiplier),
                Node::Molecule(nested) => {
                    let inner = combine_multiplier(multiplier, nested.coefficient.as_deref());
                    self.walk_nodes(&nested.nodes, inner);
                }
                _ => {}
            }
        }
    }

    fn enter_group(&mut self, group: &Group, multiplier: Option<i64>) {
        let multiplicity = group.multiplicity.clone();
        let inner = combine_multiplier(multiplier, multiplicity.as_deref());
        let bracket = if BRACKETS.contains(&group.bracket.as_str()) {
            group.bracket.clone()
        } else {
            "paren".to_string()
        };
        // The canonical walk assigns group IDs (`g1`, ...) at the
        // same point the reference's GroupRecord does, keeping the ID
        // stream identical; the wire itself carries groups by
        // atomRefs only.
        let _group_id = self.ids.next_group();
        self.groups.push(WireGroup {
            multiplicity,
            bracket,
            atom_ids: Vec::new(),
        });
        self.group_stack.push(self.groups.len() - 1);
        self.walk_nodes(&group.nodes, inner);
        self.group_stack.pop();
    }

    fn emit_atom(&mut self, atom: &Atom, multiplier: Option<i64>) {
        let id = self.ids.next_atom();
        for index in &self.group_stack {
            self.groups[*index].atom_ids.push(id.clone());
        }
        if let Some(digits) = &atom.ring_closures {
            for digit in digits.chars() {
                match self.open_rings.iter().position(|(open, _)| *open == digit) {
                    Some(index) => {
                        let (_, from) = self.open_rings.remove(index);
                        self.ring_pairs.push((from, id.clone()));
                    }
                    None => self.open_rings.push((digit, id.clone())),
                }
            }
        }
        let count = match multiplier {
            Some(multiplier) => Some(
                combine_multiplier(Some(multiplier), atom.subscript.as_deref())
                    .unwrap_or(multiplier)
                    .to_string(),
            ),
            None => atom.subscript.clone(),
        };
        if self.pending_bond.is_some() {
            if let Some(from) = self.last_atom.clone() {
                let kind = self.pending_bond.take().expect("checked above");
                self.bonds.push(WireBond {
                    id: self.ids.next_bond(),
                    from,
                    to: id.clone(),
                    order: order_code(&kind),
                });
            }
        }
        self.atoms.push(WireAtom {
            id,
            element: atom.element.clone(),
            isotope: atom.isotope.clone(),
            formal_charge: atom.charge.clone(),
            count,
            extensions: [
                atom.oxidation_state.clone(),
                atom.lone_pairs.map(|count| count.to_string()),
                atom.radical_electrons.map(|count| count.to_string()),
                atom.ring_closures.clone(),
            ],
        });
        self.last_atom = self.atoms.last().map(|atom| atom.id.clone());
    }

    /// Emits canonical bonds for each ring-closure pair, skipping
    /// pairs a positional bond already covers (the degenerate `C1-C1`).
    fn emit_ring_bonds(&mut self) {
        for (from, to) in std::mem::take(&mut self.ring_pairs) {
            let duplicate = self.bonds.iter().any(|bond| {
                (bond.from == from && bond.to == to) || (bond.from == to && bond.to == from)
            });
            if !duplicate {
                self.bonds.push(WireBond {
                    id: self.ids.next_bond(),
                    from,
                    to,
                    order: "S",
                });
            }
        }
    }
}

/// The reference walker's multiplier combine: numeric multiply, `None`
/// when neither side carries a usable integer, the usable side
/// otherwise. `None` right with a numeric left yields the multiplier
/// (a group without multiplicity inside a `2`-coefficient molecule
/// still doubles its atom counts).
fn combine_multiplier(left: Option<i64>, right: Option<&str>) -> Option<i64> {
    // Mirror Ruby's Integer(value, exception: false): only a full
    // integer literal combines; "2a" stays a raw string subscript.
    let right = right.and_then(|value| value.parse::<i64>().ok());
    match (left, right) {
        (None, None) => None,
        (left, right) => Some(left.unwrap_or(1) * right.unwrap_or(1)),
    }
}

// ---------------------------------------------------------------------------
// Read: XML → namespace-resolved tree → model. The tree is built from
// the engine's SAX events — the namespace URIs are resolved through
// the declared prefix mappings, so channel classification is exact
// regardless of the written prefix or where declarations appear.
// ---------------------------------------------------------------------------

struct XAttr {
    ns: Option<String>,
    local: String,
    value: String,
}

struct XNode {
    ns: Option<String>,
    local: String,
    attrs: Vec<XAttr>,
    children: Vec<XNode>,
    text: String,
}

impl XNode {
    /// No-namespace attribute by local name.
    fn attr(&self, local: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|attr| attr.ns.is_none() && attr.local == local)
            .map(|attr| attr.value.as_str())
    }

    /// Attribute by namespace URI + local name.
    fn attr_ns(&self, ns: &str, local: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|attr| attr.ns.as_deref() == Some(ns) && attr.local == local)
            .map(|attr| attr.value.as_str())
    }

    fn child(&self, local: &str) -> Option<&XNode> {
        self.children.iter().find(|child| child.local == local)
    }

    fn children_named<'a>(&'a self, local: &'a str) -> impl Iterator<Item = &'a XNode> {
        self.children
            .iter()
            .filter(move |child| child.local == local)
    }

    /// True when this node rides the core CML namespace (or no
    /// namespace — sloppy-but-common undeclared CML).
    fn is_core(&self) -> bool {
        self.ns.as_deref().is_none_or(|ns| ns == CML_NAMESPACE)
    }
}

#[derive(Default)]
struct TreeBuilder {
    scopes: Vec<(String, String)>,
    stack: Vec<XNode>,
    root: Option<XNode>,
    error: Option<String>,
}

/// Innermost binding for `prefix` ("" = the default namespace).
fn resolve_prefix(scopes: &[(String, String)], prefix: &str) -> Option<String> {
    scopes
        .iter()
        .rev()
        .find(|(bound, _)| bound == prefix)
        .map(|(_, uri)| uri.clone())
}

fn qname_split(name: &str) -> (&str, &str) {
    match name.split_once(':') {
        Some((prefix, local)) => (prefix, local),
        None => ("", name),
    }
}

fn build_tree(xml: &str) -> Result<XNode, CmlError> {
    // The SAX callbacks share one builder through a RefCell: the
    // engine fires them one at a time, so no reentrancy is possible.
    let builder = std::rc::Rc::new(std::cell::RefCell::new(TreeBuilder::default()));
    let mut handler = leptris::sax::Handler {
        start_prefix_mapping: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |prefix: &str, uri: &str| {
                builder
                    .borrow_mut()
                    .scopes
                    .push((prefix.to_string(), uri.to_string()));
            }
        })),
        end_prefix_mapping: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |prefix: &str| {
                let mut builder = builder.borrow_mut();
                if let Some(index) = builder
                    .scopes
                    .iter()
                    .rposition(|(bound, _)| bound == prefix)
                {
                    builder.scopes.remove(index);
                }
            }
        })),
        start_element: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |name: &str, attrs: &[(String, String)]| {
                let (prefix, local) = qname_split(name);
                let ns = resolve_prefix(&builder.borrow().scopes, prefix);
                let parsed = attrs
                    .iter()
                    .filter(|(name, _)| {
                        let (prefix, _) = qname_split(name);
                        // xmlns declarations ride the prefix-mapping
                        // events.
                        prefix != "xmlns" && name != "xmlns"
                    })
                    .map(|(name, value)| {
                        let (prefix, local) = qname_split(name);
                        // Unprefixed attributes have no namespace, per
                        // XML Namespaces — the default declaration
                        // never applies.
                        let ns = if prefix.is_empty() {
                            None
                        } else {
                            resolve_prefix(&builder.borrow().scopes, prefix)
                        };
                        XAttr {
                            ns,
                            local: local.to_string(),
                            value: value.clone(),
                        }
                    })
                    .collect();
                let mut builder = builder.borrow_mut();
                builder.stack.push(XNode {
                    ns,
                    local: local.to_string(),
                    attrs: parsed,
                    children: Vec::new(),
                    text: String::new(),
                });
            }
        })),
        characters: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |text: &str| {
                if let Some(node) = builder.borrow_mut().stack.last_mut() {
                    node.text.push_str(text);
                }
            }
        })),
        cdata: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |text: &str| {
                if let Some(node) = builder.borrow_mut().stack.last_mut() {
                    node.text.push_str(text);
                }
            }
        })),
        end_element: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |_: &str| {
                let mut builder = builder.borrow_mut();
                if let Some(done) = builder.stack.pop() {
                    match builder.stack.last_mut() {
                        Some(parent) => parent.children.push(done),
                        None => builder.root = Some(done),
                    }
                }
            }
        })),
        error: Some(Box::new({
            let builder = std::rc::Rc::clone(&builder);
            move |message: &str, _: i32, _: i32| {
                let mut builder = builder.borrow_mut();
                if builder.error.is_none() {
                    builder.error = Some(message.to_string());
                }
            }
        })),
        ..leptris::sax::Handler::default()
    };
    leptris::sax::parse(xml.as_bytes(), &mut handler)
        .map_err(|error| CmlError::Xml(error.to_string()))?;
    drop(handler);
    let mut builder = std::rc::Rc::try_unwrap(builder)
        .map_err(|_| CmlError::Xml("sax callbacks still borrowed".into()))?
        .into_inner();
    if let Some(message) = builder.error {
        return Err(CmlError::Xml(message));
    }
    builder
        .root
        .take()
        .ok_or_else(|| CmlError::Xml("empty document".into()))
}

/// Parses a CML XML string back into formula nodes.
pub fn parse_cml(xml: &str) -> Result<Vec<Node>, CmlError> {
    let root = build_tree(xml)?;
    if root.local != "cml" || !root.is_core() {
        return Err(CmlError::Invalid(format!(
            "expected a <cml> root, got <{}>",
            root.local
        )));
    }

    // Restore order mirrors the reference translator: core wire
    // (atoms carry the atom-attribute channel inline), then group
    // rebuild, then top-level insertion at preserved positions.
    let mut nodes = Vec::new();
    let mut top_level = Vec::new();
    for child in &root.children {
        if child.ns.as_deref() == Some(ACI_NAMESPACE) {
            if TOP_LEVEL_ELEMENTS.contains(&child.local.as_str()) {
                top_level.push(read_top_level(child)?);
            }
            continue;
        }
        if !child.is_core() {
            continue;
        }
        match child.local.as_str() {
            "molecule" => nodes.push(Node::Molecule(read_molecule(child)?)),
            "reaction" => nodes.push(Node::Reaction(read_reaction(child)?)),
            "reactionList" => {
                let mut steps = Vec::new();
                for reaction in child.children_named("reaction") {
                    steps.push(read_reaction(reaction)?);
                }
                nodes.push(Node::ReactionCascade(ReactionCascade { steps }));
            }
            _ => {}
        }
    }

    for entry in top_level {
        let position = entry.position.min(nodes.len());
        nodes.insert(position, entry.node);
    }
    Ok(nodes)
}

struct TopLevelEntry {
    position: usize,
    node: Node,
}

fn read_top_level(element: &XNode) -> Result<TopLevelEntry, CmlError> {
    let position = element
        .attr("position")
        .map(|value| leading_integer(value).max(0) as usize)
        .unwrap_or(0);
    let content = element.text.clone();
    // The wire carries the construct's canonical text render; wrapped
    // constructs get their delimiters back before re-parsing.
    let text = match element.local.as_str() {
        "embeddedMath" => format!("`{content}`"),
        _ => content,
    };
    let node = crate::parse_nodes(&text)
        .map_err(|error| {
            CmlError::Invalid(format!(
                "aci:{} content failed to parse: {error}",
                element.local
            ))
        })?
        .into_iter()
        .next()
        .ok_or_else(|| CmlError::Invalid(format!("aci:{} content is empty", element.local)))?;
    Ok(TopLevelEntry { position, node })
}

fn read_molecule(element: &XNode) -> Result<Molecule, CmlError> {
    let coefficient = element.attr("count").map(str::to_string);
    let mut wire_atoms = Vec::new();
    let mut wire_bonds = Vec::new();
    let mut groups = Vec::new();
    for child in &element.children {
        if child.ns.as_deref() == Some(ACI_NAMESPACE) {
            if child.local == "group" {
                groups.push(read_group(child)?);
            }
            continue;
        }
        match child.local.as_str() {
            "atomArray" => {
                for atom in child.children_named("atom") {
                    wire_atoms.push(read_atom(atom)?);
                }
            }
            "bondArray" => {
                for bond in child.children_named("bond") {
                    wire_bonds.push(read_bond(bond)?);
                }
            }
            _ => {}
        }
    }

    let mut molecule = Molecule {
        coefficient,
        ..Molecule::default()
    };
    rebuild_atoms_and_bonds(&mut molecule, wire_atoms, wire_bonds)?;
    rebuild_groups(&mut molecule, &groups);
    Ok(molecule)
}

/// A wire atom mid-restore: model fields plus the id (for bond
/// refs) and count (for group multiplicity division).
struct RestoredAtom {
    id: Option<String>,
    atom: Atom,
}

fn read_atom(element: &XNode) -> Result<RestoredAtom, CmlError> {
    let element_type = element
        .attr("elementType")
        .ok_or_else(|| CmlError::Invalid("atom without elementType".into()))?;
    let subscript = match element.attr("count") {
        None | Some("1") => None,
        Some(count) => Some(count.to_string()),
    };
    let lone_pairs = element.attr_ns(ACI_NAMESPACE, "lonePairs").map(lewis_count);
    let radical_electrons = element
        .attr_ns(ACI_NAMESPACE, "radicalElectrons")
        .map(lewis_count);
    Ok(RestoredAtom {
        id: element.attr("id").map(str::to_string),
        atom: Atom {
            element: element_type.to_string(),
            isotope: element.attr("isotope").map(str::to_string),
            subscript,
            charge: element.attr("formalCharge").map(str::to_string),
            oxidation_state: element
                .attr_ns(ACI_NAMESPACE, "oxidationState")
                .map(str::to_string),
            lone_pairs,
            radical_electrons,
            ring_closures: element
                .attr_ns(ACI_NAMESPACE, "ringClosures")
                .map(str::to_string),
            ..Atom::default()
        },
    })
}

/// The reference's lenient XML integer cast (Ruby `String#to_i`).
fn lewis_count(value: &str) -> u32 {
    leading_integer(value).max(0) as u32
}

struct RestoredBond {
    refs: Vec<String>,
    order: String,
}

fn read_bond(element: &XNode) -> Result<RestoredBond, CmlError> {
    let refs: Vec<String> = element
        .attr("atomRefs2")
        .unwrap_or_default()
        .split_whitespace()
        .map(str::to_string)
        .collect();
    if refs.len() != 2 {
        return Err(CmlError::Invalid(format!(
            "bond atomRefs2 must reference two atoms, got {:?}",
            refs
        )));
    }
    Ok(RestoredBond {
        refs,
        order: element.attr("order").unwrap_or("S").to_string(),
    })
}

fn read_group(element: &XNode) -> Result<WireGroupRead, CmlError> {
    let bracket = match element.attr("bracket") {
        Some(bracket) if BRACKETS.contains(&bracket) => bracket.to_string(),
        _ => "paren".to_string(),
    };
    Ok(WireGroupRead {
        multiplicity: element.attr("multiplicity").map(str::to_string),
        bracket,
        atom_ids: element
            .attr("atomRefs")
            .unwrap_or_default()
            .split_whitespace()
            .map(str::to_string)
            .collect(),
    })
}

struct WireGroupRead {
    multiplicity: Option<String>,
    bracket: String,
    atom_ids: Vec<String>,
}

/// Rebuilds the molecule's node list: atoms in wire order, bonds
/// inserted just before their later endpoint (the linear-spelling
/// position), ring bonds (spanning atoms) left to the atom-attribute
/// channel's ring-closure digits.
fn rebuild_atoms_and_bonds(
    molecule: &mut Molecule,
    wire_atoms: Vec<RestoredAtom>,
    wire_bonds: Vec<RestoredBond>,
) -> Result<(), CmlError> {
    let ids: Vec<Option<String>> = wire_atoms
        .iter()
        .map(|restored| restored.id.clone())
        .collect();
    let mut nodes: Vec<Node> = wire_atoms
        .into_iter()
        .map(|restored| Node::Atom(restored.atom))
        .collect();
    if wire_bonds.is_empty() {
        molecule.nodes = nodes;
        return Ok(());
    }
    let position = |id: &str| ids.iter().position(|bound| bound.as_deref() == Some(id));
    let mut placed: Vec<(&RestoredBond, usize, usize)> = wire_bonds
        .iter()
        .map(|bond| {
            let mut positions: Vec<usize> =
                bond.refs.iter().filter_map(|id| position(id)).collect();
            positions.sort_unstable();
            let later = positions.last().copied().unwrap_or(0);
            let span = positions
                .last()
                .zip(positions.first())
                .map(|(a, b)| a - b)
                .unwrap_or(0);
            (bond, later, span)
        })
        .collect();
    // Descending later position so earlier insertions don't shift the
    // pending ones.
    placed.sort_by(|(_, a, _), (_, b, _)| b.cmp(a));
    for (bond, later, span) in placed {
        // A ring bond connects non-adjacent atoms; its endpoints span
        // a gap, and the aci:ringClosures digits re-spell it.
        if span > 1 {
            continue;
        }
        nodes.insert(
            later,
            Node::Bond(Bond {
                kind: kind_from_order(&bond.order).to_string(),
                direction: None,
            }),
        );
    }
    molecule.nodes = nodes;
    Ok(())
}

/// Rebuilds Group nodes from the aci:group records. Records are
/// processed in reverse order so inner groups rebuild before their
/// outer parents (the walker records outer groups first).
fn rebuild_groups(molecule: &mut Molecule, groups: &[WireGroupRead]) {
    for record in groups.iter().rev() {
        rebuild_one_group(molecule, record);
    }
}

fn rebuild_one_group(molecule: &mut Molecule, record: &WireGroupRead) {
    let targets: Vec<usize> = record
        .atom_ids
        .iter()
        .filter_map(|id| id.strip_prefix('a')?.parse::<usize>().ok())
        .collect();
    if targets.is_empty() {
        return;
    }
    let Some((first, last)) = find_target_span(&molecule.nodes, &targets) else {
        return;
    };
    apply_group_multiplicity(
        &mut molecule.nodes,
        first,
        last,
        record.multiplicity.as_deref(),
    );
    splice_group(molecule, first, last, record);
}

/// Walks the molecule's nodes in canonical-walk order, tracking which
/// atom id each position covers; the group spans the contiguous run
/// covering any referenced atom (bonds included, preserving the chain
/// inside the group).
fn find_target_span(nodes: &[Node], targets: &[usize]) -> Option<(usize, usize)> {
    let mut first = None;
    let mut last = None;
    let mut canonical = 0usize;
    for (index, node) in nodes.iter().enumerate() {
        let count = count_atoms(node);
        let covers = (1..=count).any(|offset| targets.contains(&(canonical + offset)));
        if covers {
            first = first.or(Some(index));
            last = Some(index);
        }
        canonical += count;
    }
    first.map(|first| (first, last.unwrap_or(first)))
}

/// Counts the atoms transitively contained in a node.
fn count_atoms(node: &Node) -> usize {
    match node {
        Node::Atom(_) => 1,
        Node::Group(group) => group.nodes.iter().map(count_atoms).sum(),
        Node::Molecule(nested) => nested.nodes.iter().map(count_atoms).sum(),
        _ => 0,
    }
}

/// Divides each covered atom's subscript by the group multiplicity
/// (the emit side multiplied it in). Inner groups were already
/// divided when they rebuilt, so only bare atoms in this slice move.
fn apply_group_multiplicity(
    nodes: &mut [Node],
    first: usize,
    last: usize,
    multiplicity: Option<&str>,
) {
    let Some(multiplicity) = multiplicity else {
        return;
    };
    let divisor = leading_integer(multiplicity);
    if divisor <= 0 {
        return;
    }
    for node in &mut nodes[first..=last] {
        let Node::Atom(atom) = node else {
            continue;
        };
        if let Some(subscript) = &atom.subscript {
            let divided = leading_integer(subscript) / divisor;
            atom.subscript = if divided == 1 {
                None
            } else {
                Some(divided.to_string())
            };
        }
    }
}

/// Replaces the target run with a Group node. A bond immediately
/// before the group's first node is its attachment bond (the grammar
/// keeps it inside the group), so it moves in rather than re-spelling
/// as a stray prefix marker.
fn splice_group(molecule: &mut Molecule, first: usize, last: usize, record: &WireGroupRead) {
    let mut inner: Vec<Node> = molecule.nodes.drain(first..=last).collect();
    let mut start = first;
    if start > 0 && matches!(&molecule.nodes[start - 1], Node::Bond(_)) {
        let bond = molecule.nodes.remove(start - 1);
        inner.insert(0, bond);
        start -= 1;
    }
    molecule.nodes.insert(
        start,
        Node::Group(Group {
            nodes: inner,
            multiplicity: record.multiplicity.clone(),
            bracket: record.bracket.clone(),
        }),
    );
}

fn read_reaction(element: &XNode) -> Result<Reaction, CmlError> {
    let arrow = element
        .attr("type")
        .or_else(|| element.attr("title"))
        .map(str::to_string)
        .filter(|wire| {
            matches!(
                wire.as_str(),
                "forward" | "reverse" | "equilibrium" | "resonance"
            )
        })
        .unwrap_or_else(|| "forward".to_string());
    let mut reactants = Vec::new();
    let mut products = Vec::new();
    for child in &element.children {
        if !child.is_core() {
            continue;
        }
        match child.local.as_str() {
            "reactantList" => {
                for participant in child.children_named("reactant") {
                    reactants.push(read_participant(participant)?);
                }
            }
            "productList" => {
                for participant in child.children_named("product") {
                    products.push(read_participant(participant)?);
                }
            }
            _ => {}
        }
    }
    Ok(Reaction {
        reactants,
        products,
        arrow,
        conditions: None,
    })
}

fn read_participant(element: &XNode) -> Result<Node, CmlError> {
    let substance = element
        .child("substance")
        .ok_or_else(|| CmlError::Invalid(format!("{} without substance", element.local)))?;
    let molecule = substance
        .child("molecule")
        .ok_or_else(|| CmlError::Invalid("substance without molecule".into()))?;
    Ok(Node::Molecule(read_molecule(molecule)?))
}
