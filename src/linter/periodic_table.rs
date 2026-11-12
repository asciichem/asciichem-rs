//! Periodic table as a single source of truth for element data — the
//! Rust mirror of `AsciiChem::PeriodicTable`. Consumed by the linter
//! checks (element validation, isotope sanity, valence) and any future
//! code that needs chemically-grounded defaults.
//!
//! The subset covers the elements chemistry most commonly deals with;
//! adding an element is one table line and every consumer picks it up.

/// One element: symbol, atomic number, and the maximum valence the
/// linter's valence check allows.
pub struct Element {
    pub symbol: &'static str,
    pub atomic_number: u32,
    pub max_valence: u32,
}

/// The element table (the Ruby reference's `ELEMENTS` subset, entry
/// for entry). Sorted by atomic number for readability.
pub static ELEMENTS: [Element; 64] = [
    Element {
        symbol: "H",
        atomic_number: 1,
        max_valence: 1,
    },
    Element {
        symbol: "He",
        atomic_number: 2,
        max_valence: 0,
    },
    Element {
        symbol: "Li",
        atomic_number: 3,
        max_valence: 1,
    },
    Element {
        symbol: "Be",
        atomic_number: 4,
        max_valence: 2,
    },
    Element {
        symbol: "B",
        atomic_number: 5,
        max_valence: 3,
    },
    Element {
        symbol: "C",
        atomic_number: 6,
        max_valence: 4,
    },
    Element {
        symbol: "N",
        atomic_number: 7,
        max_valence: 3,
    },
    Element {
        symbol: "O",
        atomic_number: 8,
        max_valence: 2,
    },
    Element {
        symbol: "F",
        atomic_number: 9,
        max_valence: 1,
    },
    Element {
        symbol: "Ne",
        atomic_number: 10,
        max_valence: 0,
    },
    Element {
        symbol: "Na",
        atomic_number: 11,
        max_valence: 1,
    },
    Element {
        symbol: "Mg",
        atomic_number: 12,
        max_valence: 2,
    },
    Element {
        symbol: "Al",
        atomic_number: 13,
        max_valence: 3,
    },
    Element {
        symbol: "Si",
        atomic_number: 14,
        max_valence: 4,
    },
    Element {
        symbol: "P",
        atomic_number: 15,
        max_valence: 5,
    },
    Element {
        symbol: "S",
        atomic_number: 16,
        max_valence: 6,
    },
    Element {
        symbol: "Cl",
        atomic_number: 17,
        max_valence: 1,
    },
    Element {
        symbol: "Ar",
        atomic_number: 18,
        max_valence: 0,
    },
    Element {
        symbol: "K",
        atomic_number: 19,
        max_valence: 1,
    },
    Element {
        symbol: "Ca",
        atomic_number: 20,
        max_valence: 2,
    },
    Element {
        symbol: "Sc",
        atomic_number: 21,
        max_valence: 3,
    },
    Element {
        symbol: "Ti",
        atomic_number: 22,
        max_valence: 4,
    },
    Element {
        symbol: "V",
        atomic_number: 23,
        max_valence: 5,
    },
    Element {
        symbol: "Cr",
        atomic_number: 24,
        max_valence: 6,
    },
    Element {
        symbol: "Mn",
        atomic_number: 25,
        max_valence: 7,
    },
    Element {
        symbol: "Fe",
        atomic_number: 26,
        max_valence: 6,
    },
    Element {
        symbol: "Co",
        atomic_number: 27,
        max_valence: 4,
    },
    Element {
        symbol: "Ni",
        atomic_number: 28,
        max_valence: 4,
    },
    Element {
        symbol: "Cu",
        atomic_number: 29,
        max_valence: 4,
    },
    Element {
        symbol: "Zn",
        atomic_number: 30,
        max_valence: 2,
    },
    Element {
        symbol: "Ga",
        atomic_number: 31,
        max_valence: 3,
    },
    Element {
        symbol: "Ge",
        atomic_number: 32,
        max_valence: 4,
    },
    Element {
        symbol: "As",
        atomic_number: 33,
        max_valence: 3,
    },
    Element {
        symbol: "Se",
        atomic_number: 34,
        max_valence: 2,
    },
    Element {
        symbol: "Br",
        atomic_number: 35,
        max_valence: 1,
    },
    Element {
        symbol: "Kr",
        atomic_number: 36,
        max_valence: 0,
    },
    Element {
        symbol: "Rb",
        atomic_number: 37,
        max_valence: 1,
    },
    Element {
        symbol: "Sr",
        atomic_number: 38,
        max_valence: 2,
    },
    Element {
        symbol: "Y",
        atomic_number: 39,
        max_valence: 3,
    },
    Element {
        symbol: "Zr",
        atomic_number: 40,
        max_valence: 4,
    },
    Element {
        symbol: "Nb",
        atomic_number: 41,
        max_valence: 5,
    },
    Element {
        symbol: "Mo",
        atomic_number: 42,
        max_valence: 6,
    },
    Element {
        symbol: "Tc",
        atomic_number: 43,
        max_valence: 6,
    },
    Element {
        symbol: "Ru",
        atomic_number: 44,
        max_valence: 6,
    },
    Element {
        symbol: "Rh",
        atomic_number: 45,
        max_valence: 6,
    },
    Element {
        symbol: "Pd",
        atomic_number: 46,
        max_valence: 4,
    },
    Element {
        symbol: "Ag",
        atomic_number: 47,
        max_valence: 4,
    },
    Element {
        symbol: "Cd",
        atomic_number: 48,
        max_valence: 2,
    },
    Element {
        symbol: "In",
        atomic_number: 49,
        max_valence: 3,
    },
    Element {
        symbol: "Sn",
        atomic_number: 50,
        max_valence: 4,
    },
    Element {
        symbol: "Sb",
        atomic_number: 51,
        max_valence: 3,
    },
    Element {
        symbol: "Te",
        atomic_number: 52,
        max_valence: 2,
    },
    Element {
        symbol: "I",
        atomic_number: 53,
        max_valence: 1,
    },
    Element {
        symbol: "Xe",
        atomic_number: 54,
        max_valence: 0,
    },
    Element {
        symbol: "Cs",
        atomic_number: 55,
        max_valence: 1,
    },
    Element {
        symbol: "Ba",
        atomic_number: 56,
        max_valence: 2,
    },
    Element {
        symbol: "W",
        atomic_number: 74,
        max_valence: 6,
    },
    Element {
        symbol: "Pt",
        atomic_number: 78,
        max_valence: 4,
    },
    Element {
        symbol: "Au",
        atomic_number: 79,
        max_valence: 6,
    },
    Element {
        symbol: "Hg",
        atomic_number: 80,
        max_valence: 2,
    },
    Element {
        symbol: "Tl",
        atomic_number: 81,
        max_valence: 3,
    },
    Element {
        symbol: "Pb",
        atomic_number: 82,
        max_valence: 4,
    },
    Element {
        symbol: "Bi",
        atomic_number: 83,
        max_valence: 3,
    },
    Element {
        symbol: "U",
        atomic_number: 92,
        max_valence: 6,
    },
];

fn element(symbol: &str) -> Option<&'static Element> {
    ELEMENTS.iter().find(|element| element.symbol == symbol)
}

/// Whether the symbol is in the table.
pub fn known(symbol: &str) -> bool {
    element(symbol).is_some()
}

/// Atomic number for a symbol, or `None` for unknown elements.
pub fn atomic_number(symbol: &str) -> Option<u32> {
    element(symbol).map(|element| element.atomic_number)
}

/// Maximum valence for a symbol, or `None` for unknown elements.
pub fn max_valence(symbol: &str) -> Option<u32> {
    element(symbol).map(|element| element.max_valence)
}

/// Every known symbol, in table order.
pub fn symbols() -> impl Iterator<Item = &'static str> {
    ELEMENTS.iter().map(|element| element.symbol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn looks_up_the_common_elements() {
        assert_eq!(atomic_number("C"), Some(6));
        assert_eq!(max_valence("C"), Some(4));
        assert_eq!(atomic_number("U"), Some(92));
        assert_eq!(max_valence("He"), Some(0));
    }

    #[test]
    fn unknown_symbols_are_unknown() {
        assert!(!known("Xx"));
        assert_eq!(atomic_number("Hx"), None);
        assert_eq!(max_valence("Cy"), None);
    }

    #[test]
    fn covers_every_entry_of_the_ruby_table() {
        // The mirror contract: same subset, same numbers. Spot the
        // reference's boundary elements per period.
        for (symbol, z, valence) in [
            ("H", 1, 1),
            ("B", 5, 3),
            ("Cl", 17, 1),
            ("Fe", 26, 6),
            ("I", 53, 1),
            ("Ba", 56, 2),
            ("W", 74, 6),
            ("Bi", 83, 3),
            ("U", 92, 6),
        ] {
            assert_eq!(atomic_number(symbol), Some(z), "{symbol}");
            assert_eq!(max_valence(symbol), Some(valence), "{symbol}");
        }
        assert_eq!(symbols().count(), 64);
    }
}
