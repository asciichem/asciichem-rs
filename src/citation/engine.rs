//! The leptris engine's programmatic construction face, as the
//! citation track consumes it. The published `leptris` binding wraps
//! the engine's parse/query/serialize surface; the engine's builder
//! API (`leptris_document_create`, `leptris_element_create_child`,
//! `leptris_element_set_attribute`, `leptris_element_set_text`,
//! `leptris_document_set_root` — all public C API, `src/include/
//! leptris.h`) is declared here so documents can be BUILT through the
//! engine too: the bibitem XML never passes through hand-rolled XML
//! strings — escaping, well-formedness, and the emitted bytes all
//! come from libleptris, the same engine that serializes. The
//! declarations ride the `leptris` dependency's link directives
//! (`links = "leptris"`), so the same `LEPTRIS_LIB_PATH` build is
//! exercised on both sides.
//!
//! Upstream note: this face is a candidate for the binding itself
//! (leptris-rs wraps the same functions in its wasm face); when the
//! native binding grows it, this module shrinks to a re-export.

use super::Error;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

/// Engine status codes (`leptris/types.h`): the zero status is OK.
const LEPTRIS_OK: c_int = 0;

/* Opaque engine handles (leptris/types.h — pointer typedefs). */
#[repr(C)]
struct LeptrisDocument {
    _private: [u8; 0],
}
#[repr(C)]
struct LeptrisElement {
    _private: [u8; 0],
}

/// The frozen serialize options (leptris/types.h
/// `LeptrisSerializeOptions` — the ABI is size-pinned in the header).
#[repr(C)]
struct SerializeOptions {
    indent: c_int,
    xml_declaration: c_int,
    encoding: *const c_char,
}

/// The post-1.9.x serialization extensions (leptris/types.h
/// `LeptrisSerializeExtOptions`); rides the size-aware entry so the
/// header struct can grow independently.
#[repr(C)]
struct SerializeExtOptions {
    indent_text: c_int,
    indent_unit: *const c_char,
    expand_empty: c_int,
    html_method: c_int,
    indent_html: c_int,
}

// The native link rides THIS crate's rlib metadata (the same
// mechanism the `leptris` binding uses via its build script): every
// final binary that links `asciichem` — bins, tests, examples — gets
// `-l leptris` and the propagated `-L` search directory from the
// binding's build script.
#[link(name = "leptris")]
extern "C" {
    fn leptris_document_create() -> *mut LeptrisDocument;
    fn leptris_document_free(doc: *mut LeptrisDocument);
    fn leptris_document_set_root(doc: *mut LeptrisDocument, root: *mut LeptrisElement) -> c_int;
    fn leptris_element_create(
        doc: *mut LeptrisDocument,
        name: *const c_char,
    ) -> *mut LeptrisElement;
    fn leptris_element_create_child(
        parent: *mut LeptrisElement,
        name: *const c_char,
    ) -> *mut LeptrisElement;
    fn leptris_element_set_attribute(
        elem: *mut LeptrisElement,
        name: *const c_char,
        value: *const c_char,
    ) -> c_int;
    fn leptris_element_set_text(elem: *mut LeptrisElement, text: *const c_char) -> c_int;
    fn leptris_document_serialize_ext_sized(
        doc: *mut LeptrisDocument,
        options: *const SerializeOptions,
        ext: *const SerializeExtOptions,
        ext_size: usize,
    ) -> *mut c_char;
    fn leptris_free_string(string: *mut c_char);
    fn leptris_error_message(status: c_int) -> *const c_char;
}

/// A document being built, owned until dropped
/// (`leptris_document_free`).
pub(crate) struct Document {
    raw: *mut LeptrisDocument,
}

impl Document {
    /// An empty document with its own memory pool.
    pub fn new() -> Result<Self, Error> {
        // SAFETY: no arguments; a NULL return is an allocation failure.
        let raw = unsafe { leptris_document_create() };
        if raw.is_null() {
            return Err(Error::Engine(
                "the engine could not allocate a document".to_string(),
            ));
        }
        Ok(Document { raw })
    }

    /// Attach an element as the document's root element. The element
    /// must have been created against this document and must be
    /// unattached.
    pub fn set_root(&self, root: Element) -> Result<(), Error> {
        // SAFETY: both handles are engine-issued and live (the root
        // element is unattached; the document owns the pool).
        let status = unsafe { leptris_document_set_root(self.raw, root.raw) };
        status_to_result(status)
    }

    /// Serialize the document: two-space indentation, no XML
    /// declaration, self-closing empty elements — the Relaton
    /// bibitem emission shape the goldens pin.
    pub fn serialize(&self) -> Result<String, Error> {
        let encoding = b"UTF-8\0";
        let options = SerializeOptions {
            indent: 2,
            xml_declaration: 0,
            encoding: encoding.as_ptr() as *const c_char,
        };
        let ext = SerializeExtOptions {
            indent_text: 0,
            indent_unit: std::ptr::null(),
            expand_empty: 0,
            html_method: 0,
            indent_html: 0,
        };
        // SAFETY: the document handle is live; the option structs
        // borrow only stack memory for the duration of the call; the
        // returned string is engine-allocated and released below.
        let raw = unsafe {
            leptris_document_serialize_ext_sized(
                self.raw,
                &options,
                &ext,
                std::mem::size_of::<SerializeExtOptions>(),
            )
        };
        if raw.is_null() {
            return Err(Error::Engine(
                "the engine could not serialize the document".to_string(),
            ));
        }
        // SAFETY: the engine returned a NUL-terminated string.
        let serialized = unsafe { CStr::from_ptr(raw) }
            .to_string_lossy()
            .into_owned();
        // SAFETY: the string came from the engine's own allocator.
        unsafe { leptris_free_string(raw) };
        Ok(serialized)
    }
}

impl Drop for Document {
    fn drop(&mut self) {
        // SAFETY: the handle is engine-issued and not yet freed.
        unsafe { leptris_document_free(self.raw) };
    }
}

/// An element in a [`Document`] — a borrowed view: the element is
/// pool-owned by its document and lives as long as the document.
#[derive(Clone, Copy)]
pub(crate) struct Element {
    raw: *mut LeptrisElement,
}

impl Element {
    /// Create an element against the document's pool (unattached).
    pub fn create(document: &Document, name: &str) -> Result<Self, Error> {
        let name = c_string(name)?;
        // SAFETY: the document handle is live; the name outlives the
        // call and is pool-copied by the engine.
        let raw = unsafe { leptris_element_create(document.raw, name.as_ptr()) };
        element_from_raw(raw)
    }

    /// Create an element and append it to this element in one call —
    /// the engine's fused builder entry.
    pub fn create_child(&self, name: &str) -> Result<Self, Error> {
        let name = c_string(name)?;
        // SAFETY: the parent handle is live; the name is pool-copied.
        let raw = unsafe { leptris_element_create_child(self.raw, name.as_ptr()) };
        element_from_raw(raw)
    }

    /// Set an attribute (created if absent, replaced if present).
    pub fn set_attribute(&self, name: &str, value: &str) -> Result<(), Error> {
        let name = c_string(name)?;
        let value = c_string(value)?;
        // SAFETY: both handles and strings are live for the call; the
        // engine copies them into the element's pool.
        let status =
            unsafe { leptris_element_set_attribute(self.raw, name.as_ptr(), value.as_ptr()) };
        status_to_result(status)
    }

    /// Set the element's text content (a text child).
    pub fn set_text(&self, text: &str) -> Result<(), Error> {
        let text = c_string(text)?;
        // SAFETY: the element handle is live; the text is pool-copied.
        let status = unsafe { leptris_element_set_text(self.raw, text.as_ptr()) };
        status_to_result(status)
    }
}

fn c_string(value: &str) -> Result<CString, Error> {
    CString::new(value).map_err(|_| Error::Engine("string contains an interior NUL".to_string()))
}

fn element_from_raw(raw: *mut LeptrisElement) -> Result<Element, Error> {
    if raw.is_null() {
        return Err(Error::Engine(
            "the engine could not create the element".to_string(),
        ));
    }
    Ok(Element { raw })
}

fn status_to_result(status: c_int) -> Result<(), Error> {
    if status == LEPTRIS_OK {
        return Ok(());
    }
    // SAFETY: the engine returns a static message string.
    let message = unsafe { CStr::from_ptr(leptris_error_message(status)) }
        .to_string_lossy()
        .into_owned();
    Err(Error::Engine(message))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_and_serializes_a_small_document() {
        let document = Document::new().expect("document");
        let root = Element::create(&document, "bibitem").expect("root");
        root.set_attribute("type", "dataset").expect("attribute");
        let title = root.create_child("title").expect("title");
        title.set_attribute("type", "main").expect("title type");
        title.set_text("Aspirin & friends").expect("title text");
        let role = root.create_child("role").expect("role");
        role.set_attribute("type", "publisher").expect("role type");
        document.set_root(root).expect("set root");

        let xml = document.serialize().expect("serialize");
        // The engine escaped the ampersand, kept the attribute order,
        // and self-closed the empty element.
        assert!(xml.contains(r#"<bibitem type="dataset">"#), "{xml}");
        assert!(
            xml.contains(r#"<title type="main">Aspirin &amp; friends</title>"#),
            "{xml}"
        );
        assert!(xml.contains(r#"<role type="publisher"/>"#), "{xml}");
    }
}
