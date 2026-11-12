//! The Relaton bibitem emission shape: a typed [`Bibitem`] payload
//! serialized into the XML the metanorma plugin consumes —
//! `<bibitem type="dataset">` with a `<docidentifier type="source">`,
//! the `<date type="accessed">`, the publisher contributor, the
//! provenance note, and every identifier as a `<keyword>` (nested
//! vocab element). The element order, attribute order, and the
//! derived `id` mirror relaton-bib 2.1.9's XML mapping byte for byte;
//! the goldens in `tests/citation_golden.rs` pin the exact bytes.
//!
//! Construction and serialization both run through the leptris engine
//! ([`super::engine`]) — no hand-rolled XML strings.

use super::engine::{Document, Element};
use super::Error;

/// The relaton-models schema version relaton-bib 2.1.9 renders on
/// every bibitem (`Relaton.schema_versions["relaton-models"]`); the
/// goldens pin it, so a relaton ride that bumps it re-records them
/// deliberately.
const SCHEMA_VERSION: &str = "v1.5.6";

/// The dataset bibitem payload: the version-independent fields the
/// citation track assembles, translated to the wire shape here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bibitem {
    /// The bibitem's Relaton type — always `dataset` for substances.
    pub bibitem_type: String,
    /// Main title: `{title base} - {publisher} substance record`.
    pub title: String,
    /// The document identifier scheme — the citing source's name.
    pub docidentifier_type: String,
    /// The document identifier content.
    pub docidentifier: String,
    /// The publisher (the source's organization).
    pub publisher: String,
    /// The accessed date, `YYYY-MM-DD`.
    pub accessed_on: String,
    /// The source record link, when the profile knows the shape.
    pub link: Option<String>,
    /// Every identifier as `convention=value` — cross-checking keys.
    pub keywords: Vec<String>,
    /// The provenance record as a bibliographic note; absent when it
    /// would say nothing beyond the source name itself.
    pub note: Option<String>,
}

impl Bibitem {
    /// Serializes through the leptris engine: the document is built
    /// with the engine's construction API and emitted with its
    /// serializer — the byte shape the goldens pin.
    pub fn to_xml(&self) -> Result<String, Error> {
        let document = Document::new()?;
        let root = Element::create(&document, "bibitem")?;
        if let Some(id) = derived_id(&self.docidentifier) {
            root.set_attribute("id", &id)?;
        }
        root.set_attribute("type", &self.bibitem_type)?;
        root.set_attribute("schema-version", SCHEMA_VERSION)?;

        let title = root.create_child("title")?;
        title.set_attribute("type", "main")?;
        title.set_text(&self.title)?;

        if let Some(link) = &self.link {
            let uri = root.create_child("uri")?;
            uri.set_attribute("type", "src")?;
            uri.set_text(link)?;
        }

        let docidentifier = root.create_child("docidentifier")?;
        docidentifier.set_attribute("type", &self.docidentifier_type)?;
        docidentifier.set_text(&self.docidentifier)?;

        let date = root.create_child("date")?;
        date.set_attribute("type", "accessed")?;
        let on = date.create_child("on")?;
        on.set_text(&self.accessed_on)?;

        let contributor = root.create_child("contributor")?;
        let role = contributor.create_child("role")?;
        role.set_attribute("type", "publisher")?;
        let organization = contributor.create_child("organization")?;
        let name = organization.create_child("name")?;
        name.set_text(&self.publisher)?;

        if let Some(note) = &self.note {
            let note_element = root.create_child("note")?;
            note_element.set_attribute("type", "provenance")?;
            note_element.set_text(note)?;
        }

        for keyword in &self.keywords {
            let keyword_element = root.create_child("keyword")?;
            let vocab = keyword_element.create_child("vocab")?;
            vocab.set_text(keyword)?;
        }

        document.set_root(root)?;
        document.serialize()
    }
}

/// The bibitem `id`: relaton derives it from the document identifier
/// by stripping every non-word character (ASCII letters, digits,
/// underscore). An identifier that reduces to nothing carries no id.
fn derived_id(docidentifier: &str) -> Option<String> {
    let id: String = docidentifier
        .chars()
        .filter(|character| character.is_ascii_alphanumeric() || *character == '_')
        .collect();
    if id.is_empty() {
        None
    } else {
        Some(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bibitem(docidentifier: &str) -> Bibitem {
        Bibitem {
            bibitem_type: "dataset".to_string(),
            title: "Aspirin - example substance record".to_string(),
            docidentifier_type: "example".to_string(),
            docidentifier: docidentifier.to_string(),
            publisher: "example".to_string(),
            accessed_on: "2026-01-05".to_string(),
            link: None,
            keywords: vec!["cas=50-78-2".to_string()],
            note: Some("Retrieved from example".to_string()),
        }
    }

    #[test]
    fn derives_ids_like_relaton() {
        assert_eq!(
            derived_id("PubChem CID 2244"),
            Some("PubChemCID2244".to_string())
        );
        assert_eq!(derived_id("cas: 50-78-2"), Some("cas50782".to_string()));
        assert_eq!(
            derived_id("Test: A-1_b.2/3 x"),
            Some("TestA1_b23x".to_string())
        );
        // Non-ASCII word characters do not survive.
        assert_eq!(derived_id("Ünïcödé"), Some("ncd".to_string()));
        // Nothing survives: no id attribute at all.
        assert_eq!(derived_id("--..--"), None);
    }

    #[test]
    fn serializes_the_relaton_element_order() {
        let xml = bibitem("cas: 50-78-2").to_xml().expect("serialize");
        let id_position = xml.find(r#"id="cas50782""#).expect("id");
        let type_position = xml.find(r#"type="dataset""#).expect("type");
        let schema_position = xml.find(r#"schema-version"#).expect("schema");
        assert!(id_position < type_position && type_position < schema_position);
        let title_position = xml.find("<title").expect("title");
        let docidentifier_position = xml.find("<docidentifier").expect("docidentifier");
        let date_position = xml.find("<date").expect("date");
        let contributor_position = xml.find("<contributor").expect("contributor");
        let note_position = xml.find(r#"<note type="provenance">"#).expect("note");
        let keyword_position = xml.find("<keyword").expect("keyword");
        assert!(title_position < docidentifier_position);
        assert!(docidentifier_position < date_position);
        assert!(date_position < contributor_position);
        assert!(contributor_position < note_position);
        assert!(note_position < keyword_position);
    }

    #[test]
    fn drops_the_note_and_uri_when_absent() {
        let item = Bibitem {
            link: None,
            note: None,
            ..bibitem("cas: 50-78-2")
        };
        let xml = item.to_xml().expect("serialize");
        assert!(!xml.contains("<note"), "{xml}");
        assert!(!xml.contains("<uri"), "{xml}");
    }

    #[test]
    fn engine_parses_back_what_the_builder_wrote() {
        let item = bibitem("PubChem CID 2244");
        let xml = item.to_xml().expect("serialize");
        let parsed = leptris::Document::parse(xml.as_bytes()).expect("parse");
        let root = parsed.root().expect("root");
        assert_eq!(root.name(), "bibitem");
        let docidentifier = root.child("docidentifier").expect("docidentifier");
        assert_eq!(docidentifier.attribute("type"), Some("example"));
        assert_eq!(docidentifier.text(), Some("PubChem CID 2244"));
    }
}
