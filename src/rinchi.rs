//! The reaction-identity layer: RInChI — a reaction (reactants +
//! products + agents, each a molecule) gets its own identifier
//! assembled from the per-component standard InChIs. Derived and
//! decoded IN-PROCESS through the [chemomer](https://github.com/chemomer/chemomer)
//! crate, exactly as the identity layer ([`crate::identity`]) derives
//! per-molecule identifiers — the same oracle-gated native engine, the
//! same byte-exact posture.
//!
//! Two directions, one error type (the engine's own — no second error
//! vocabulary):
//!
//! - **derive** — a reaction's V2000 RXNfile/RDfile text reads into a
//!   [`chemomer::RinchiReaction`] ([`reaction_from_rxnfile`] /
//!   [`reaction_from_rdfile`]; the RDfile edge is the `$DATUM $MFMT`
//!   agent scan), and assembles into the identity record
//!   ([`derive`] — the five lines: identifier, RAuxInfo, the three
//!   RInChIKey shapes). The equilibrium reading rides
//!   [`Directionality`]/`chemomer::RinchiReaction::set_directionality`
//!   (the reference CLI's `/force-equilibrium`).
//! - **decode** — a RInChI string (plus its RAuxInfo companion when
//!   known) reads back into the reaction it identifies
//!   ([`decode`] / [`decode_with_auxinfo`]): the components carry
//!   their identities straight from the string (no structures are
//!   rebuilt), and decode∘derive is the identity — the canonical
//!   string re-assembles byte-for-byte.
//!
//! The record's strings are engine-verbatim (key headers included);
//! the shared corpus records the keys headerless and the runners
//! strip, mirroring the engine's own gate.

pub use chemomer::rinchi::Directionality;

/// The assembled reaction-identity record: the identifier, its
/// RAuxInfo companion line, and the three RInChIKey shapes (the short
/// hashed one under `rinchikey`, per the shared corpus's field
/// naming). Every string is engine-verbatim — key headers included.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RinchiRecord {
    pub rinchi: String,
    pub rauxinfo: String,
    pub rinchikey: String,
    pub long_key: String,
    pub web_key: String,
}

/// Reads a reaction's V2000 RXNfile text into the reaction model —
/// the derive direction's input edge, mirroring
/// [`crate::identity::derive_from_molfile`]'s molfile edge.
pub fn reaction_from_rxnfile(rxnfile: &str) -> chemomer::Result<chemomer::RinchiReaction> {
    chemomer::parse_rxnfile(rxnfile)
}

/// Reads an MDL RDfile (the RXN block verbatim plus the `$DATUM
/// $MFMT` agent scan) into the reaction model — the second input
/// edge.
pub fn reaction_from_rdfile(rdfile: &str) -> chemomer::Result<chemomer::RinchiReaction> {
    chemomer::parse_rdfile(rdfile)
}

/// Assembles the reaction's identity record — the derive direction
/// itself (engine-verbatim lines, key headers included).
pub fn derive(reaction: &chemomer::RinchiReaction) -> chemomer::Result<RinchiRecord> {
    Ok(RinchiRecord {
        rinchi: reaction.rinchi_string()?,
        rauxinfo: reaction.rinchi_auxinfo()?,
        rinchikey: reaction.rinchi_short_key()?,
        long_key: reaction.rinchi_long_key()?,
        web_key: reaction.rinchi_web_key()?,
    })
}

/// Reads + assembles in one call — the one-shot over
/// [`reaction_from_rxnfile`] + [`derive`].
pub fn derive_from_rxnfile(rxnfile: &str) -> chemomer::Result<RinchiRecord> {
    derive(&reaction_from_rxnfile(rxnfile)?)
}

/// Reads + assembles in one call over the RDfile edge.
pub fn derive_from_rdfile(rdfile: &str) -> chemomer::Result<RinchiRecord> {
    derive(&reaction_from_rdfile(rdfile)?)
}

/// Decodes a RInChI string into the reaction it identifies — the
/// components carry their identities straight from the string; the
/// direction flag recovers which group held which role. The
/// components' aux fields stay empty; pair with the RAuxInfo
/// companion through [`decode_with_auxinfo`].
pub fn decode(rinchi: &str) -> chemomer::Result<chemomer::RinchiReaction> {
    chemomer::reaction_from_rinchi(rinchi)
}

/// The two-line decode: the identifier plus its RAuxInfo companion,
/// whose per-component aux bodies distribute over the decoded
/// components in group order.
pub fn decode_with_auxinfo(
    rinchi: &str,
    rauxinfo: &str,
) -> chemomer::Result<chemomer::RinchiReaction> {
    chemomer::rinchi::reaction_from_rinchi_with_auxinfo(rinchi, rauxinfo)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The corpus's water-is-water RXNfile (inchi-tests v0.3.41,
    // rinchi/water-is-water): one water drawn each side.
    const WATER_RXN: &str = "$RXN\nwater-is-water\n  chemomer recording\n\n  1  1\n$MOL\n\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$MOL\n\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const WATER_RAUXINFO_GROUP: &str = "0/N:1/rA:1nO/rB:/rC:12,12,0;";

    #[test]
    fn an_rxnfile_derives_the_recorded_five_lines() {
        let record = derive_from_rxnfile(WATER_RXN).expect("derives");
        assert_eq!(record.rinchi, "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+");
        assert_eq!(
            record.rauxinfo,
            format!("RAuxInfo=1.00.1/{WATER_RAUXINFO_GROUP}<>{WATER_RAUXINFO_GROUP}")
        );
        assert_eq!(
            record.rinchikey,
            "Short-RInChIKey=SA-FUHFF-XLYOFNOQVP-XLYOFNOQVP-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ"
        );
        assert_eq!(
            record.long_key,
            "Long-RInChIKey=SA-FUHFF-XLYOFNOQVPJJNP-UHFFFAOYSA-N--XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        );
        assert_eq!(
            record.web_key,
            "Web-RInChIKey=XLYOFNOQVPJJNPIXP-NUHFFFADPSCTJSA"
        );
    }

    #[test]
    fn the_equilibrium_reading_is_forced_not_parsed() {
        // The RXNfile carries no direction field — parsing yields
        // directional, the canonical water-is-water assembly reads
        // forward, and only the forced reading speaks `/d=` and the E
        // key codes (the reference CLI's /force-equilibrium).
        let mut reaction = reaction_from_rxnfile(WATER_RXN).expect("reads");
        assert_eq!(reaction.directionality(), Directionality::Directional);
        assert_eq!(
            derive(&reaction).expect("derives").rinchi,
            "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+"
        );
        reaction.set_directionality(Directionality::Equilibrium);
        let record = derive(&reaction).expect("derives");
        assert_eq!(record.rinchi, "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d=");
        assert_eq!(
            record.rinchikey,
            "Short-RInChIKey=SA-EUHFF-XLYOFNOQVP-XLYOFNOQVP-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ"
        );
    }

    #[test]
    fn a_rinchi_string_decodes_into_its_roles() {
        // A reversed recording (`/d-`: the products group sorts
        // first): the reader recovers the real roles from the flag.
        let reaction = decode("RInChI=1.00.1S/<>ClH/h1H!H3N/h1H3/d-").expect("decodes");
        assert_eq!(reaction.reactants().len(), 2);
        assert_eq!(reaction.products().len(), 0);
        assert_eq!(
            reaction.reactants()[0].inchi_string().expect("inchi"),
            "InChI=1S/ClH/h1H"
        );
        assert_eq!(reaction.directionality(), Directionality::Directional);

        // The equilibrium flag reads equilibrium.
        let equilibrium = decode("RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d=").expect("decodes");
        assert_eq!(equilibrium.directionality(), Directionality::Equilibrium);
    }

    #[test]
    fn decode_then_derive_is_the_identity() {
        // The RAuxInfo companion distributes over the decoded
        // components; re-assembling the decoded reaction reproduces
        // both lines byte-for-byte — every component's aux/key content
        // rides the assembled lines.
        let rinchi = "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+";
        let rauxinfo = format!("RAuxInfo=1.00.1/{WATER_RAUXINFO_GROUP}<>{WATER_RAUXINFO_GROUP}");
        let decoded = decode_with_auxinfo(rinchi, &rauxinfo).expect("decodes");
        let record = derive(&decoded).expect("re-assembles");
        assert_eq!(record.rinchi, rinchi);
        assert_eq!(record.rauxinfo, rauxinfo);

        // The bare decode (no companion at hand) re-assembles the
        // identifier alone.
        let bare = decode(rinchi).expect("decodes");
        assert_eq!(derive(&bare).expect("re-assembles").rinchi, rinchi);
    }

    #[test]
    fn garbage_is_rejected_through_the_engine_error() {
        let error = decode("RInChI=0.9/ nonsense").expect_err("rejected");
        assert!(error.message().contains("rxn"), "{}", error.message());
        let error = derive_from_rxnfile("not an rxnfile").expect_err("rejected");
        assert!(error.message().contains("$RXN"), "{}", error.message());
    }
}
