//! Identifier validators (TODO.v2/17 slice 4) — the Rust mirror of
//! `AsciiChem::Identifiers`: lexical validation per convention, one
//! [`diagnostic`] dispatcher the linter's format check consumes.
//!
//! These are *lexical* checks: format, check digits, element symbols.
//! Deep verification is the resolution/engine layer's job — a linter
//! stays offline and total.

/// The convention's validator's reason for rejecting `value`, or None
/// when the value is valid or the convention has no validator
/// (skipped silently — resolution is not linting).
pub fn diagnostic(convention: &str, value: &str) -> Option<String> {
    match convention {
        "cas" => cas::diagnostic(value),
        "inchi" => inchi::diagnostic(value),
        "inchikey" => inchikey::diagnostic(value),
        "smiles" => smiles::diagnostic(value),
        _ => None,
    }
}

/// CAS Registry Number: DDDDDDD-DD-D with one check digit equal to the
/// weighted sum of the preceding digits (weights 1..n counted from the
/// right) modulo 10. Reimplementation of the published algorithm
/// (github.com/simonengelke/CAS_Validation, MIT) — not copied.
pub mod cas {
    fn digits(value: &str) -> Vec<u32> {
        value
            .chars()
            .filter(|c| *c != '-')
            .map(|c| c as u32 - '0' as u32)
            .collect()
    }

    pub fn diagnostic(value: &str) -> Option<String> {
        let valid_shape = {
            let parts: Vec<&str> = value.split('-').collect();
            parts.len() == 3
                && parts[0].len() >= 2
                && parts[0].len() <= 7
                && parts[1].len() == 2
                && parts[2].len() == 1
                && value.chars().all(|c| c.is_ascii_digit() || c == '-')
        };
        if !valid_shape {
            return Some(format!(
                "CAS RN {value:?} is not in DDDDDDD-DD-D form (2–7 digits, 2 digits, 1 check digit)"
            ));
        }
        let numbers = digits(value);
        let check = *numbers.last().unwrap_or(&0);
        let weighted: u32 = numbers[..numbers.len() - 1]
            .iter()
            .rev()
            .enumerate()
            .map(|(index, digit)| digit * (index as u32 + 1))
            .sum::<u32>()
            % 10;
        if check == weighted {
            None
        } else {
            Some(format!(
                "CAS RN {value:?}: check digit {check} does not match weighted sum ({weighted})"
            ))
        }
    }
}

/// InChIKey: 27 uppercase letters in XXXXXXXXXXXXXX-XXXXXXXXXX-X
/// (blocks of 14, 10, and 1). Structure verification is the InChI
/// engine's job; format alone catches the common corruptions.
pub mod inchikey {
    pub fn diagnostic(value: &str) -> Option<String> {
        let blocks: Vec<&str> = value.split('-').collect();
        let shaped = blocks.len() == 3
            && blocks[0].len() == 14
            && blocks[1].len() == 10
            && blocks[2].len() == 1
            && value.chars().all(|c| c.is_ascii_uppercase() || c == '-');
        if shaped {
            None
        } else {
            Some(format!(
                "InChIKey {value:?} must be 27 uppercase letters in XXXXXXXXXXXXXX-XXXXXXXXXX-X form (blocks of 14, 10, and 1)"
            ))
        }
    }
}

/// InChI: versioned prefix plus a well-formed formula layer whose
/// element symbols exist in the periodic table. Layers after the
/// formula layer are not validated — full InChI semantics is the
/// InChI software's job (TODO.v2 10), not a linter's.
pub mod inchi {
    use crate::linter::periodic_table;
    use std::collections::BTreeMap;

    /// Element symbol => atom count from the formula layer, or None
    /// when the value is not analysable. Single source of truth for
    /// offline InChI/molecule composition checks.
    pub fn formula_counts(value: &str) -> Option<BTreeMap<String, i64>> {
        let rest = value
            .strip_prefix("InChI=1S/")
            .or(value.strip_prefix("InChI=1/"))?;
        let layer = rest.split('/').next().unwrap_or("");
        parse_formula_layer(layer)?.ok()
    }

    pub fn diagnostic(value: &str) -> Option<String> {
        let Some(layer) = formula_layer(value) else {
            return Some(format!(
                "InChI {value:?} must start with \"InChI=1/\" or \"InChI=1S/\""
            ));
        };
        match parse_formula_layer(layer) {
            None => Some(format!("InChI {value:?}: formula layer is empty")),
            Some(Err(reason)) => Some(format!("InChI {value:?}: {reason}")),
            Some(Ok(_)) => None,
        }
    }

    /// The formula layer: the segment between the prefix and the first
    /// subsequent "/" (e.g. "C9H8O4" in "InChI=1S/C9H8O4/c1-6(10)...").
    fn formula_layer(value: &str) -> Option<&str> {
        let rest = value
            .strip_prefix("InChI=1S/")
            .or(value.strip_prefix("InChI=1/"))?;
        Some(rest.split('/').next().unwrap_or(""))
    }

    /// Ok(counts) on success; Err carries the Ruby-matching reason.
    fn parse_formula_layer(layer: &str) -> Option<Result<BTreeMap<String, i64>, String>> {
        if layer.is_empty() {
            return Some(Err("formula layer is empty".to_string()));
        }
        let mut counts = BTreeMap::new();
        let mut rest = layer;
        while !rest.is_empty() {
            match next_element_token(rest) {
                Err(reason) => return Some(Err(reason)),
                Ok((symbol, count, remainder)) => {
                    *counts.entry(symbol.to_string()).or_insert(0) += count;
                    rest = remainder;
                }
            }
        }
        Some(Ok(counts))
    }

    /// Consumes the element token at the head of `rest`. Err carries
    /// the Ruby-matching reason for a malformed or unknown token.
    fn next_element_token(rest: &str) -> Result<(&str, i64, &str), String> {
        let first = rest.chars().next().unwrap_or(' ');
        if !first.is_ascii_uppercase() {
            return Err(format!("malformed formula layer near {rest:?}"));
        }
        let take = if rest.len() >= 2 && rest.as_bytes()[1].is_ascii_lowercase() {
            2
        } else {
            1
        };
        let symbol = &rest[..take];
        if !periodic_table::known(symbol) {
            return Err(format!(
                "unknown element symbol {symbol:?} in formula layer"
            ));
        }
        let digits: usize = rest[take..]
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .count();
        let (count, consumed) = if digits == 0 {
            (1, take)
        } else {
            match rest[take..take + digits].parse::<i64>() {
                Ok(value) => (value, take + digits),
                Err(_) => return Err(format!("malformed formula layer near {rest:?}")),
            }
        };
        Ok((symbol, count, &rest[consumed..]))
    }
}

/// SMILES sanity validation. Deliberately *not* a SMILES parser — it
/// catches corruption: whitespace, unbalanced brackets and
/// parentheses, unpaired ring-bond digits, and element tokens outside
/// the periodic table. Lowercase aromatic atoms (c1ccccc1) are valid
/// SMILES and pass.
pub mod smiles {
    use crate::linter::periodic_table;

    const ORGANIC: [&str; 10] = ["B", "C", "N", "O", "P", "S", "F", "I", "Cl", "Br"];
    const AROMATIC: [&str; 10] = ["b", "c", "n", "o", "p", "s", "f", "i", "cl", "br"];

    /// Sorted unique element symbols appearing in the SMILES, or None
    /// when the value is not analysable.
    pub fn element_set(value: &str) -> Option<Vec<String>> {
        if structural_reasons(value).is_some() {
            return None;
        }
        let mut symbols: Vec<String> = Vec::new();
        let mut index = 0;
        while index < value.len() {
            let (next, analysable) = advance(value, index, &mut symbols)?;
            if !analysable {
                return None;
            }
            index = next;
        }
        symbols.sort();
        symbols.dedup();
        Some(symbols)
    }

    pub fn diagnostic(value: &str) -> Option<String> {
        if value.is_empty() {
            return Some("SMILES must not be empty".to_string());
        }
        if value.chars().any(|c| c.is_whitespace()) {
            return Some(format!("SMILES {value:?} must not contain whitespace"));
        }
        if let Some(reason) = structural_reasons(value) {
            return Some(reason);
        }
        if element_set(value).is_none() {
            return Some(format!(
                "SMILES {value:?}: contains a malformed or unknown element token"
            ));
        }
        None
    }

    /// None when structurally plausible; a reason otherwise.
    fn structural_reasons(value: &str) -> Option<String> {
        if !balanced(value, '[', ']') {
            return Some(format!(
                "SMILES {value:?}: square brackets are not balanced"
            ));
        }
        if !balanced(value, '(', ')') {
            return Some(format!("SMILES {value:?}: parentheses are not balanced"));
        }
        // Ring digits only exist outside brackets — bracket interiors
        // hold isotope/H-count digits, which must not pair up.
        let mut bare = String::new();
        let mut depth = 0;
        for c in value.chars() {
            match c {
                '[' => depth += 1,
                ']' => depth -= 1,
                _ if depth == 0 => bare.push(c),
                _ => {}
            }
        }
        for digit in b'0'..=b'9' {
            let count = bare.chars().filter(|c| *c as u8 == digit).count();
            if count % 2 == 1 {
                return Some(format!(
                    "SMILES {value:?}: ring-bond digit {} appears an odd number of times (ring closures must pair)",
                    digit as char
                ));
            }
        }
        None
    }

    fn balanced(value: &str, opener: char, closer: char) -> bool {
        let mut depth = 0i32;
        for c in value.chars() {
            if c == opener {
                depth += 1;
            }
            if c == closer {
                depth -= 1;
            }
            if depth < 0 {
                return false;
            }
        }
        depth == 0
    }

    /// Advances past the token at `index`, appending its element
    /// symbol. Returns (new_index, analysable).
    /// None when the token at `index` is unanalysable (an unknown
    /// element) — the whole value is then unanalysable.
    fn advance(value: &str, index: usize, symbols: &mut Vec<String>) -> Option<(usize, bool)> {
        let rest = &value[index..];
        let first = rest.chars().next()?;
        if first == '[' {
            Some(consume_bracket_atom(value, index, symbols))
        } else if first.is_ascii_alphabetic() {
            // An unknown element token makes the whole value
            // unanalysable — None, not a partial set (Ruby semantics;
            // the value-absolute index cannot loop).
            let (consumed, symbol) = organic_token(rest)?;
            symbols.push(symbol);
            Some((index + consumed, true))
        } else {
            Some((index + 1, true))
        }
    }

    fn consume_bracket_atom(value: &str, index: usize, symbols: &mut Vec<String>) -> (usize, bool) {
        let Some(close_rel) = value[index..].find(']') else {
            return (0, false);
        };
        let content = &value[index + 1..index + close_rel];
        match bracket_atom_symbol(content) {
            Some(symbol) => {
                symbols.push(symbol.to_string());
                (index + close_rel + 1, true)
            }
            None => (0, false),
        }
    }

    /// Two-letter element first (Cl/Br, aromatic cl/br), then the
    /// single-letter organic/aromatic subsets. Returns (consumed,
    /// symbol) or None.
    fn organic_token(rest: &str) -> Option<(usize, String)> {
        let two: String = rest.chars().take(2).collect();
        let one: String = rest.chars().take(1).collect();
        if ORGANIC.contains(&two.as_str()) {
            Some((2, two))
        } else if AROMATIC.contains(&two.as_str()) {
            Some((2, capitalize(&two)))
        } else if ORGANIC.contains(&one.as_str()) {
            Some((1, one))
        } else if AROMATIC.contains(&one.as_str()) {
            Some((1, capitalize(&one)))
        } else {
            None
        }
    }

    fn capitalize(symbol: &str) -> String {
        let mut chars = symbol.chars();
        match chars.next() {
            Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
            None => String::new(),
        }
    }

    /// Bracket atom content: optional isotope digits, then the symbol.
    /// Wildcards make composition unknown → None.
    fn bracket_atom_symbol(content: &str) -> Option<&str> {
        let after_digits = content.trim_start_matches(|c: char| c.is_ascii_digit());
        let bytes = after_digits.as_bytes();
        // One uppercase letter, optionally followed by one lowercase
        // (Cl, Br) — "NH4+" reads N, never NH.
        if bytes.is_empty() || !bytes[0].is_ascii_uppercase() {
            return None;
        }
        let take = if bytes.len() >= 2 && bytes[1].is_ascii_lowercase() {
            2
        } else {
            1
        };
        let symbol = &after_digits[..take];
        (symbol != "*" && periodic_table::known(symbol)).then_some(symbol)
    }
}
