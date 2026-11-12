//! Molfile (CTfile V2000) ingestion and emission — the Rust mirror of
//! AsciiChem::Molfile (TODO.v2 09). Molfile is the highest-fidelity
//! structure path: atom coordinates are preserved on the model
//! (x2/y2/z2), charges via the `M  CHG` property block (plus the
//! legacy pre-CHG charge column), isotopes via the `M  ISO` block
//! (the mass-difference field is deferred to it — mapping it
//! unambiguously requires isotope tables), bond type 4 is aromatic
//! and marks its atoms, and bond stereo codes 1/6 become wedge/hash.
//!
//! Tetrahedral chirality crosses as a WEDGE on the stereocentre's
//! first bond (TH1 → up, TH2 → down), oriented centre → neighbour:
//! the oracle ignores the V2000 atom-parity field for 2D input, so
//! the parity field alone never carries stereo. Directional bonds
//! (`/`, `\`, E/Z) have no V2000 field: the model's bond direction is
//! documented as a molfile loss rather than silently invented from
//! geometry.

use std::collections::BTreeMap;
use std::collections::HashMap;

use crate::layout;
use crate::model::{Atom, Molecule};
use crate::structure::{build_graph_strict, linearize, Edge};
use crate::ParseError;

type Result<T> = std::result::Result<T, ParseError>;

/// V2000 bond type per model bond kind (`wedge`/`hash` ride the
/// single-bond type with a stereo code).
fn bond_type(kind: &str) -> Option<i64> {
    match kind {
        "single" | "wedge" | "hash" => Some(1),
        "double" => Some(2),
        "triple" => Some(3),
        "aromatic" => Some(4),
        _ => None,
    }
}

/// V2000 bond stereo code per model bond kind (0 = not stereo).
fn bond_stereo(kind: &str) -> i64 {
    match kind {
        "wedge" => 1,
        "hash" => 6,
        _ => 0,
    }
}

/// MDL atom parity → model chirality: 1 odd (TH1), 2 even (TH2),
/// 0/3 not stereo (or either). The parity is relative to the
/// neighbour order in the bond block, which the linearizer preserves —
/// the same ordering the SMILES TH1/TH2 convention uses.
fn chirality_from_parity(parity: i64) -> Option<String> {
    match parity {
        1 => Some("TH1".to_string()),
        2 => Some("TH2".to_string()),
        _ => None,
    }
}

/// MDL H-count code → model hydrogens: 0 default, 1 H0, 2 H1, 3 H2...
fn hydrogens_from_hcount(code: i64) -> Option<u32> {
    u32::try_from(code - 1).ok()
}

// -- fixed-format field readers (Ruby's `String#[start, len]` + lenient
// `to_i`/`to_f` coercion; molfile fields are ASCII fixed columns) ------

fn slice(line: &str, start: usize, len: usize) -> String {
    let bytes = line.as_bytes();
    if start >= bytes.len() {
        return String::new();
    }
    let end = (start + len).min(bytes.len());
    String::from_utf8_lossy(&bytes[start..end]).to_string()
}

fn field(lines: &[String], line_index: usize, start: usize, len: usize) -> String {
    lines
        .get(line_index)
        .map(|line| slice(line, start, len).trim().to_string())
        .unwrap_or_default()
}

/// Ruby's lenient `String#to_i`: leading whitespace, optional sign,
/// leading digit run; otherwise 0.
fn lenient_to_i(text: &str) -> i64 {
    let trimmed = text.trim_start();
    let mut chars = trimmed.chars();
    let mut sign = 1i64;
    if let Some(first) = chars.next() {
        match first {
            '-' => sign = -1,
            '+' => {}
            digit if digit.is_ascii_digit() => {
                let mut value = digit.to_digit(10).expect("digit") as i64;
                for next in chars {
                    match next.to_digit(10) {
                        Some(digit) => value = value * 10 + digit as i64,
                        None => break,
                    }
                }
                return sign * value;
            }
            _ => return 0,
        }
    }
    let mut value = 0i64;
    let mut digits = false;
    for next in chars {
        match next.to_digit(10) {
            Some(digit) => {
                digits = true;
                value = value * 10 + digit as i64;
            }
            None => break,
        }
    }
    if digits {
        sign * value
    } else {
        0
    }
}

/// Ruby's lenient `String#to_f`: leading whitespace, optional sign,
/// digits with an optional fraction and exponent; otherwise 0.
fn lenient_to_f(text: &str) -> f64 {
    let trimmed = text.trim_start();
    let bytes = trimmed.as_bytes();
    let mut end = 0usize;
    if end < bytes.len() && (bytes[end] == b'-' || bytes[end] == b'+') {
        end += 1;
    }
    let digits_start = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    let int_digits = end > digits_start;
    if end < bytes.len() && bytes[end] == b'.' {
        end += 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
    }
    if !int_digits {
        // No digits at all ("abc", "+", ".") → 0, like to_f.
        return 0.0;
    }
    let mut number_end = end;
    if end < bytes.len() && (bytes[end] == b'e' || bytes[end] == b'E') {
        let mut exp_end = end + 1;
        if exp_end < bytes.len() && (bytes[exp_end] == b'-' || bytes[exp_end] == b'+') {
            exp_end += 1;
        }
        let exp_digits = exp_end;
        while exp_end < bytes.len() && bytes[exp_end].is_ascii_digit() {
            exp_end += 1;
        }
        if exp_end > exp_digits {
            number_end = exp_end;
        }
    }
    trimmed[..number_end].parse().unwrap_or(0.0)
}

/// Ruby `text.lines.map(&:chomp)` — \n-terminated segments with the
/// trailing newline (and carriage return) removed; a trailing newline
/// does not add a final empty line.
fn chomp_lines(text: &str) -> Vec<String> {
    text.split_inclusive('\n')
        .map(|line| line.trim_end_matches(['\n', '\r']).to_string())
        .collect()
}

// -- parser --------------------------------------------------------------

/// Parses a V2000 molfile into a `Node::Molecule` (fixed-format
/// blocks: header (3 lines), counts line, atom block, bond block,
/// property block).
pub fn parse_molfile(text: &str) -> Result<Molecule> {
    let lines = chomp_lines(text);
    if lines.len() < 5 {
        return Err(ParseError {
            message: "molfile too short".into(),
        });
    }

    let atom_field = field(&lines, 3, 0, 3);
    let bond_field = field(&lines, 3, 3, 3);
    let parse_count = |text: &str| -> Option<usize> {
        (!text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit()))
            .then(|| text.parse().expect("digits"))
    };
    let (Some(atom_count), Some(bond_count)) = (parse_count(&atom_field), parse_count(&bond_field))
    else {
        return Err(ParseError {
            message: format!("malformed counts line: {:?}", lines[3]),
        });
    };

    let mut atoms = parse_atom_block(&lines, atom_count)?;
    let mut bonds = parse_bond_block(&lines, atom_count, bond_count)?;

    // A wedge/hash bond on an atom whose parity field already carries
    // chirality is that chirality's wedge (our writer emits exactly
    // that pair), not an independent bond kind — fold it back to its
    // plain single bond.
    let chiral: Vec<usize> = atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.chirality.is_some())
        .map(|(index, _)| index)
        .collect();
    for bond in &mut bonds {
        if matches!(bond.kind.as_str(), "wedge" | "hash")
            && (chiral.contains(&bond.from) || chiral.contains(&bond.to))
        {
            bond.kind = "single".to_string();
        }
    }

    apply_legacy_charge(&lines, atom_count, &mut atoms);
    let (charges, isotopes) = parse_properties(&lines);
    for (index, value) in charges {
        if let Some(atom) = atoms.get_mut(index) {
            atom.charge = Some(signed_charge(value));
        }
    }
    for (index, mass) in isotopes {
        if let Some(atom) = atoms.get_mut(index) {
            atom.isotope = Some(mass.to_string());
        }
    }

    // V2000 carries aromaticity on bonds (type 4); the model carries
    // it on atoms and bonds, so atoms touching an aromatic bond are
    // marked aromatic.
    for bond in &bonds {
        if bond.kind == "aromatic" {
            if let Some(atom) = atoms.get_mut(bond.from) {
                atom.aromatic = Some(true);
            }
            if let Some(atom) = atoms.get_mut(bond.to) {
                atom.aromatic = Some(true);
            }
        }
    }

    let mut adjacency: Vec<BTreeMap<usize, String>> = vec![BTreeMap::new(); atoms.len()];
    for bond in &bonds {
        adjacency[bond.from].insert(bond.to, bond.kind.clone());
        adjacency[bond.to].insert(bond.from, bond.kind.clone());
    }
    let mut edges = Vec::new();
    for (index, neighbours) in adjacency.iter().enumerate() {
        for (to, kind) in neighbours {
            if *to > index {
                edges.push(Edge {
                    from: index,
                    to: *to,
                    kind: kind.clone(),
                    direction: None,
                });
            }
        }
    }

    let nodes = linearize(atoms, &edges)?;
    Ok(Molecule {
        nodes,
        ..Molecule::default()
    })
}

struct RawBond {
    from: usize,
    to: usize,
    kind: String,
}

fn parse_atom_block(lines: &[String], count: usize) -> Result<Vec<Atom>> {
    (1..=count)
        .map(|index| {
            let line = lines.get(3 + index).ok_or_else(|| ParseError {
                message: format!("truncated atom block (expected {count} atoms)"),
            })?;
            let element = slice(line, 31, 3).trim().to_string();
            if element.is_empty() {
                return Err(ParseError {
                    message: format!("atom {index} has no element symbol"),
                });
            }
            Ok(Atom {
                element,
                x2: Some(lenient_to_f(&slice(line, 0, 10))),
                y2: Some(lenient_to_f(&slice(line, 10, 10))),
                z2: Some(lenient_to_f(&slice(line, 20, 10))),
                // The mass-difference field encodes isotopes relative
                // to the rounded average mass; mapping it unambiguously
                // requires isotope tables, so v1 defers isotopes to the
                // explicit `M  ISO` block.
                chirality: chirality_from_parity(lenient_to_i(&slice(line, 39, 3))),
                hydrogens: hydrogens_from_hcount(lenient_to_i(&slice(line, 42, 3))),
                ..Atom::default()
            })
        })
        .collect()
}

fn parse_bond_block(lines: &[String], atom_count: usize, count: usize) -> Result<Vec<RawBond>> {
    (1..=count)
        .map(|index| {
            let line = lines
                .get(3 + atom_count + index)
                .ok_or_else(|| ParseError {
                    message: format!("truncated bond block (expected {count} bonds)"),
                })?;
            let from = lenient_to_i(&slice(line, 0, 3)) - 1;
            let to = lenient_to_i(&slice(line, 3, 3)) - 1;
            let bond_type = lenient_to_i(&slice(line, 6, 3));
            let stereo = lenient_to_i(&slice(line, 9, 3));
            if from < 0 || to < 0 || from >= atom_count as i64 || to >= atom_count as i64 {
                return Err(ParseError {
                    message: format!("bond {index} has out-of-range atom indexes"),
                });
            }
            let kind = match stereo {
                1 => "wedge".to_string(),
                6 => "hash".to_string(),
                _ => match bond_type {
                    1 => "single".to_string(),
                    2 => "double".to_string(),
                    3 => "triple".to_string(),
                    4 => "aromatic".to_string(),
                    other => {
                        return Err(ParseError {
                            message: format!("unsupported molfile bond type {other}"),
                        })
                    }
                },
            };
            Ok(RawBond {
                from: from as usize,
                to: to as usize,
                kind,
            })
        })
        .collect()
}

/// Pre-CHG charge column (0-based 36, 3): nonzero values 1..4 mean
/// +1..+4, 5..7 mean -1..-3.
fn apply_legacy_charge(lines: &[String], atom_count: usize, atoms: &mut [Atom]) {
    for index in 1..=atom_count {
        let code = lenient_to_i(&field(lines, 3 + index, 36, 3));
        if code == 0 {
            continue;
        }
        let charge = if code <= 4 { code } else { 4 - code };
        if let Some(atom) = atoms.get_mut(index - 1) {
            atom.charge = Some(signed_charge(charge));
        }
    }
}

/// The signed integer charge → the model's number-then-sign form
/// ("2+", "-").
fn signed_charge(value: i64) -> String {
    let sign = if value < 0 { "-" } else { "+" };
    match value.abs() {
        1 => sign.to_string(),
        magnitude => format!("{magnitude}{sign}"),
    }
}

fn parse_properties(lines: &[String]) -> (HashMap<usize, i64>, HashMap<usize, i64>) {
    let mut charges = HashMap::new();
    let mut isotopes = HashMap::new();
    for line in lines {
        if line.starts_with("M  CHG") || line.starts_with("M  ISO") {
            let is_charge = line.starts_with("M  CHG");
            // "M  CHGnnn8_ppp_vvv..." — first token is the pair count;
            // the pairs follow.
            let tail = slice(line, 6, line.len());
            let parts: Vec<&str> = tail.split_whitespace().collect();
            for pair in parts[1..].chunks(2) {
                if pair.len() != 2 {
                    continue;
                }
                let (Some(index), value) =
                    (lenient_to_i(pair[0]).checked_sub(1), lenient_to_i(pair[1]))
                else {
                    continue;
                };
                if index < 0 {
                    continue;
                }
                if is_charge {
                    charges.insert(index as usize, value);
                } else {
                    isotopes.insert(index as usize, value);
                }
            }
        }
    }
    (charges, isotopes)
}

// -- writer --------------------------------------------------------------

const CHIRALITY_PARITY: [(&str, i64); 2] = [("TH1", 1), ("TH2", 2)];

/// Emits a V2000 molfile. Authored x2/y2 coordinates are used when
/// present; otherwise a deterministic 2D layout is computed (the
/// layout walks atoms in the same order as the structure graph, so
/// positions map by index). Charges and isotopes are emitted as
/// `M  CHG` / `M  ISO` property lines.
pub fn write_molfile(molecule: &Molecule, name: &str) -> Result<String> {
    write_molfile_with(molecule, name, true)
}

/// Emits a V2000 molfile with all coordinates at the origin — the
/// Ruby identity seam's `layout: false` call. Geometry must never
/// invent stereochemistry a downstream perceiver (InChI) would read:
/// the coordinates field is deliberately empty, stereo crossing only
/// through the wedge/hash codes.
pub fn write_molfile_unlaid(molecule: &Molecule, name: &str) -> Result<String> {
    write_molfile_with(molecule, name, false)
}

fn write_molfile_with(molecule: &Molecule, name: &str, layout: bool) -> Result<String> {
    let (atoms, edges) = build_graph_strict(molecule)?;
    if edges.is_empty() && atoms.len() > 1 {
        return Err(ParseError {
            message: "molecule has no bonds — a formula is not a structure".into(),
        });
    }

    // Layout fallback for atoms without authored coordinates. Only
    // trusted when the layout walk saw the same atoms the graph walk
    // did; otherwise unbonded coordinates stay at the origin.
    let layout_positions = layout.then(|| layout::layout_molecule(molecule).1);
    let layout_aligned = layout_positions
        .as_ref()
        .is_some_and(|positions| positions.len() == atoms.len());

    let mut out = String::new();
    out.push_str(name);
    out.push('\n');
    out.push_str("  AsciiChem\n");
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
        atoms.len(),
        edges.len()
    ));
    out.push('\n');

    for (index, atom) in atoms.iter().enumerate() {
        // Laid: authored coordinates win, the deterministic layout
        // covers the rest. Unlaid: the coordinate pair stays at the
        // origin (the Ruby `layout: false` seam behaviour); z passes
        // through in both.
        let placed = layout_positions
            .as_ref()
            .filter(|_| layout_aligned)
            .and_then(|positions| positions.get(index));
        let x = if layout {
            atom.x2.or(placed.map(|p| p.x)).unwrap_or(0.0)
        } else {
            0.0
        };
        let y = if layout {
            atom.y2.or(placed.map(|p| p.y)).unwrap_or(0.0)
        } else {
            0.0
        };
        let z = atom.z2.unwrap_or(0.0);
        let parity = CHIRALITY_PARITY
            .iter()
            .find(|(chirality, _)| Some(&chirality.to_string()) == atom.chirality.as_ref())
            .map_or(0, |(_, parity)| *parity);
        let hcount = atom.hydrogens.map_or(0, |count| count + 1);
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0{:>3}{:>3}  0  0  0  0  0  0  0  0",
            x, y, z, atom.element, parity, hcount
        ));
        out.push('\n');
    }

    // One wedge per chiral atom, on its first bond, oriented
    // centre → neighbour: TH1 rises (stereo 1), TH2 falls (6).
    let mut wedges: HashMap<usize, (usize, usize, i64)> = HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        let Some(stereo) = (match atom.chirality.as_deref() {
            Some("TH1") => Some(1),
            Some("TH2") => Some(6),
            _ => None,
        }) else {
            continue;
        };
        let Some((edge_index, other)) = edges.iter().enumerate().find_map(|(edge_index, edge)| {
            if edge.from == index {
                Some((edge_index, edge.to))
            } else if edge.to == index {
                Some((edge_index, edge.from))
            } else {
                None
            }
        }) else {
            continue;
        };
        wedges.insert(edge_index, (index, other, stereo));
    }

    for (edge_index, edge) in edges.iter().enumerate() {
        let bond_type = bond_type(&edge.kind).ok_or_else(|| ParseError {
            message: format!("{} bonds have no molfile V2000 type", edge.kind),
        })?;
        let stereo = bond_stereo(&edge.kind);
        let (from, to, stereo) = match wedges.get(&edge_index) {
            Some(&(from, to, stereo)) => (from, to, stereo),
            None => (edge.from, edge.to, stereo),
        };
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0  0  0  0  0",
            from + 1,
            to + 1,
            bond_type,
            stereo
        ));
        out.push('\n');
    }

    let charge_pairs: Vec<(usize, i64)> = atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.charge.is_some())
        .map(|(index, atom)| {
            (
                index + 1,
                charge_value(atom.charge.as_deref().unwrap_or("")),
            )
        })
        .collect();
    let isotope_pairs: Vec<(usize, i64)> = atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.isotope.is_some())
        .map(|(index, atom)| {
            (
                index + 1,
                lenient_to_i(atom.isotope.as_deref().unwrap_or("")),
            )
        })
        .collect();
    if let Some(line) = property_line("M  CHG", &charge_pairs) {
        out.push_str(&line);
        out.push('\n');
    }
    if let Some(line) = property_line("M  ISO", &isotope_pairs) {
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  END\n");
    Ok(out)
}

/// The model's number-then-sign charge ("2+") → signed integer.
fn charge_value(charge: &str) -> i64 {
    let count: i64 = charge
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .unwrap_or(1);
    if charge.ends_with('-') {
        -count
    } else {
        count
    }
}

fn property_line(prefix: &str, pairs: &[(usize, i64)]) -> Option<String> {
    if pairs.is_empty() {
        return None;
    }
    let mut line = format!("{prefix}{:>3}", pairs.len());
    for (index, value) in pairs {
        line.push_str(&format!("{index:>4}{value:>4}"));
    }
    Some(line)
}
