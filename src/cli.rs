//! The `asciichem` command-line interface (TODO.v2/17 slice 8;
//! asciichem-rs#16): the eleven commands the Ruby CLI ships, on ONE
//! dispatch so the contract — stdout shapes and exit codes — is
//! testable offline. The binary (`src/bin/asciichem.rs`) is a thin
//! main over [`run`]; the same function is the test seam. Resolve and
//! cite take the network through the [`Fetch`] seam the binary wires
//! to [`crate::resolver::NetFetch`] (the real network path); tests
//! inject a fixture fetch and a temporary cache — slice 5's
//! recorded-fixture pattern, now driving the whole dispatch.
//!
//! Contract: mirror the Ruby Thor CLI
//! (`asciichem-ruby/lib/asciichem/cli.rb`) — parse errors exit 1,
//! unknown `convert` output formats 2, resolve errors 3, cite errors
//! 4, identity errors 5, parsing-engine selection 6; `lint` exits 1
//! on error-severity findings, `roundtrip` 1 on a mismatch,
//! `validate` 1 on a known-but-invalid identifier. Documented
//! deviations (version and help text wording, `--engine` accepting
//! only parsanol) live in the command handlers and the PR.

use std::fs;
use std::io::Write;
use std::path::Path;

use clap::{ArgAction, CommandFactory, Parser, Subcommand};

use crate::citation;
use crate::cml;
use crate::identifiers;
use crate::linter;
use crate::model::{Identifier, Node};
use crate::molfile;
use crate::resolver::adapter::{Adapter, Fetch};
use crate::resolver::cache::Cache;
use crate::resolver::{Registry, ResolverError};
use crate::smiles;
use crate::wire;
use crate::{html, identity, latex, mathml, svg, ParseError};

/// The conventions with an offline validator (the reference
/// `Identifiers` registry's keys) — "known" in `validate`'s status
/// column. Conventions without one (`iupac`, `cid`, `chebi`) report
/// `unknown convention` and never affect the exit code.
const KNOWN_CONVENTIONS: [&str; 4] = ["cas", "inchi", "inchikey", "smiles"];

/// The `--engine` value this binary accepts. The Rust core IS the
/// parsanol-powered engine (the compiled PARG artifact on parsanol's
/// portable engine) — the flag is accepted for compatibility with the
/// Ruby CLI's `--engine parsanol`, and every other value is an
/// engine-unavailable error (exit 6), mirroring the reference's
/// missing-gem posture.
const ENGINE: &str = "parsanol";

/// The command environment: the network seam, the resolution cache,
/// and the source registry. The shipped binary takes the real network
/// path (`NetFetch`, the user cache dir, the shipped sources); tests
/// inject the offline equivalents.
pub struct Deps<'a> {
    pub fetch: &'a dyn Fetch,
    pub cache: &'a Cache,
    pub registry: &'a Registry,
}

#[derive(Debug, Parser)]
#[command(
    name = "asciichem",
    version = env!("CARGO_PKG_VERSION"),
    disable_version_flag = true,
    about = "AsciiChem: ASCII syntax for chemistry parsed into a semantic model"
)]
struct Cli {
    /// Print the asciichem version
    #[arg(short = 'v', long = "version", action = ArgAction::Version)]
    _version: (),

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Resolve a substance and emit a dataset-type Relaton bibitem (XML)
    Cite(CiteArgs),
    /// Convert INPUT to FORMAT (mathml|text|html|latex|svg|structural-svg|model-json|cml|smiles|molfile)
    Convert(ConvertArgs),
    /// Derive InChI/InChIKey locally from the structure (offline; the chemomer engine is built in)
    Identity(IdentityArgs),
    /// Run chemistry checks; exit 1 on error, 0 if clean
    Lint(LintArgs),
    /// Parse CML XML and emit AsciiChem text
    #[command(name = "parse-cml", alias = "parse_cml")]
    ParseCml(ParseCmlArgs),
    /// Resolve a substance from a source (network; cached)
    Resolve(ResolveArgs),
    /// Parse and re-emit; exit non-zero if not equal
    Roundtrip(RoundtripArgs),
    /// Print a tree of all available commands
    Tree,
    /// Offline identifier validation
    Validate(ValidateArgs),
    /// Print the asciichem version
    Version,
}

#[derive(Debug, Parser)]
struct CiteArgs {
    /// CAS registry number
    #[arg(long)]
    cas: Option<String>,
    /// Substance name
    #[arg(long)]
    name: Option<String>,
    /// PubChem CID
    #[arg(long)]
    cid: Option<String>,
    /// InChIKey
    #[arg(long)]
    inchikey: Option<String>,
    /// SMILES
    #[arg(long)]
    smiles: Option<String>,
    /// Resolver source
    #[arg(long, default_value = "pubchem")]
    source: String,
    /// Bypass the cache
    #[arg(long)]
    refresh: bool,
}

#[derive(Debug, Parser)]
struct ConvertArgs {
    /// Source text (or '-' for stdin)
    #[arg(short, long)]
    input: Option<String>,
    /// Read source from a file
    #[arg(short = 'f', long)]
    file: Option<String>,
    /// Input grammar: asciichem|smiles|molfile
    #[arg(long, default_value = "asciichem")]
    from: String,
    /// Output format
    #[arg(short = 't', long, default_value = "mathml")]
    format: String,
    /// Parsing engine: parsanol (the only engine — this binary runs the PARG artifact on parsanol's portable engine)
    #[arg(long, default_value = ENGINE)]
    engine: String,
}

#[derive(Debug, Parser)]
struct IdentityArgs {
    /// Source text (or '-' for stdin)
    #[arg(short, long)]
    input: Option<String>,
    /// Read source from a file
    #[arg(short = 'f', long)]
    file: Option<String>,
    /// Input grammar: asciichem|smiles|molfile
    #[arg(long, default_value = "asciichem")]
    from: String,
    /// Path to the inchi-1 binary (not supported: the Rust engine derives in-process)
    #[arg(long)]
    engine_bin: Option<String>,
}

#[derive(Debug, Parser)]
struct LintArgs {
    /// AsciiChem source text
    #[arg(short, long)]
    input: String,
    /// Output format: text or json
    #[arg(short = 'f', long, default_value = "text")]
    format: String,
}

#[derive(Debug, Parser)]
struct ParseCmlArgs {
    /// CML XML source
    #[arg(short, long)]
    input: String,
}

#[derive(Debug, Parser)]
struct ResolveArgs {
    /// CAS registry number
    #[arg(long)]
    cas: Option<String>,
    /// Substance name
    #[arg(long)]
    name: Option<String>,
    /// PubChem CID
    #[arg(long)]
    cid: Option<String>,
    /// InChIKey
    #[arg(long)]
    inchikey: Option<String>,
    /// SMILES
    #[arg(long)]
    smiles: Option<String>,
    /// Resolver source
    #[arg(long, default_value = "pubchem")]
    source: String,
    /// Bypass the cache
    #[arg(long)]
    refresh: bool,
    /// Output: model-json | text | smiles
    #[arg(short = 't', long, default_value = "model-json")]
    format: String,
}

#[derive(Debug, Parser)]
struct RoundtripArgs {
    #[arg(short, long)]
    input: String,
}

#[derive(Debug, Parser)]
struct ValidateArgs {
    #[arg(short, long)]
    input: String,
}

/// A command failure: the stderr line (already prefixed —
/// "Parse error: ", "Resolve error: ", ...) and the exit code the
/// reference maps it to. `line` is empty when the command already
/// wrote its own stderr block (roundtrip's mismatch report).
struct Failure {
    code: i32,
    line: String,
}

fn failure(code: i32, prefix: &str, message: impl std::fmt::Display) -> Failure {
    Failure {
        code,
        line: format!("{prefix}{message}"),
    }
}

/// The mismatch/invalid outcomes that print no Failure line (their
/// stderr block, if any, is already written).
fn silent(code: i32) -> Failure {
    Failure {
        code,
        line: String::new(),
    }
}

/// Ruby `puts` semantics: exactly one trailing newline, never two
/// (the molfile writer's output already ends with one).
fn emit(out: &mut dyn Write, text: &str) {
    let _ = out.write_all(text.as_bytes());
    if !text.ends_with('\n') {
        let _ = out.write_all(b"\n");
    }
}

/// Runs the CLI. The returned i32 is the process exit code. clap's
/// own usage errors (unknown flags, missing required arguments) exit
/// 1 — the Thor posture — and help/version render to stdout with
/// exit 0.
pub fn run<I, T>(arguments: I, deps: &Deps<'_>, out: &mut dyn Write, err: &mut dyn Write) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = match Cli::try_parse_from(arguments) {
        Ok(cli) => cli,
        Err(error) => {
            if error.use_stderr() {
                let _ = write!(err, "{error}");
                return 1;
            }
            let _ = write!(out, "{error}");
            return 0;
        }
    };
    let Some(command) = cli.command else {
        // Bare `asciichem`: the command list, exit 0.
        let _ = write!(out, "{}", Cli::command().render_help());
        return 0;
    };
    let outcome = match command {
        Commands::Cite(args) => cmd_cite(&args, deps, out),
        Commands::Convert(args) => cmd_convert(&args, out),
        Commands::Identity(args) => cmd_identity(&args, out),
        Commands::Lint(args) => cmd_lint(&args, out),
        Commands::ParseCml(args) => cmd_parse_cml(&args, out),
        Commands::Resolve(args) => cmd_resolve(&args, deps, out),
        Commands::Roundtrip(args) => cmd_roundtrip(&args, out, err),
        Commands::Tree => cmd_tree(out),
        Commands::Validate(args) => cmd_validate(&args, out),
        Commands::Version => {
            emit(out, &version_line());
            Ok(())
        }
    };
    match outcome {
        Ok(()) => 0,
        Err(failed) => {
            if !failed.line.is_empty() {
                let _ = writeln!(err, "{}", failed.line);
            }
            failed.code
        }
    }
}

fn version_line() -> String {
    format!("asciichem {}", env!("CARGO_PKG_VERSION"))
}

// -- source ingestion ---------------------------------------------------

enum SourceError {
    /// Neither -i nor -f (convert's guard; identity instead parses
    /// the empty string, mirroring the reference).
    Missing,
    Io(std::io::Error),
}

/// The source text: -f wins over -i, `-` reads stdin (the reference
/// `read_source` order).
fn source_text(input: Option<&str>, file: Option<&str>) -> Result<String, SourceError> {
    if let Some(path) = file {
        return fs::read_to_string(path).map_err(SourceError::Io);
    }
    match input {
        None => Err(SourceError::Missing),
        Some("-") => {
            use std::io::Read;
            let mut buffer = String::new();
            std::io::stdin()
                .read_to_string(&mut buffer)
                .map_err(SourceError::Io)?;
            Ok(buffer)
        }
        Some(text) => Ok(text.to_string()),
    }
}

/// One ingestion point per input grammar (TODO.v2 09): every grammar
/// funnels into the same model, so every output format works
/// regardless of the input language. Molfile input may be a file path.
fn ingest(source: &str, from: &str) -> Result<Vec<Node>, ParseError> {
    match from {
        "asciichem" => crate::parse_nodes(source),
        "smiles" => smiles::parse_smiles(source),
        "molfile" => {
            let text = if Path::new(source).is_file() {
                fs::read_to_string(source).map_err(|error| ParseError {
                    message: error.to_string(),
                })?
            } else {
                source.to_string()
            };
            molfile::parse_molfile(&text).map(|molecule| vec![Node::Molecule(molecule)])
        }
        other => Err(ParseError {
            message: format!("unknown --from grammar: {other}"),
        }),
    }
}

/// Molecules anywhere at the top level: direct formula nodes plus
/// reaction terms. Identity derivation needs exactly one molecule.
fn molecules_in(nodes: &[Node]) -> Vec<&crate::model::Molecule> {
    let mut molecules = Vec::new();
    for node in nodes {
        match node {
            Node::Reaction(reaction) => {
                collect_molecules(reaction.reactants.iter(), &mut molecules);
                collect_molecules(reaction.products.iter(), &mut molecules);
            }
            Node::ReactionCascade(cascade) => {
                for step in &cascade.steps {
                    collect_molecules(step.reactants.iter(), &mut molecules);
                    collect_molecules(step.products.iter(), &mut molecules);
                }
            }
            Node::Molecule(molecule) => molecules.push(molecule),
            _ => {}
        }
    }
    molecules
}

fn collect_molecules<'a>(
    nodes: impl Iterator<Item = &'a Node>,
    molecules: &mut Vec<&'a crate::model::Molecule>,
) {
    for node in nodes {
        if let Node::Molecule(molecule) = node {
            molecules.push(molecule);
        }
    }
}

// -- convert ------------------------------------------------------------

fn cmd_convert(args: &ConvertArgs, out: &mut dyn Write) -> Result<(), Failure> {
    if args.engine != ENGINE {
        return Err(failure(
            6,
            "Engine error: ",
            format!("unknown engine {:?} (available: {ENGINE})", args.engine),
        ));
    }
    let source = match source_text(args.input.as_deref(), args.file.as_deref()) {
        Ok(text) => text,
        Err(SourceError::Missing) => {
            return Err(failure(1, "Parse error: ", "provide -i INPUT or -f FILE"));
        }
        Err(SourceError::Io(error)) => return Err(failure(1, "Error: ", error)),
    };
    let nodes = ingest(&source, &args.from).map_err(|error| failure(1, "Parse error: ", error))?;
    let rendered = render(&nodes, &args.format)?;
    emit(out, &rendered);
    Ok(())
}

/// Emit-side rejections (an unbonded formula's SMILES, a wire
/// construct outside the subset) ride the parse family — the
/// reference's convert maps every non-format failure to exit 1.
fn parse_family(error: impl std::fmt::Display) -> Failure {
    failure(1, "Parse error: ", error)
}

fn render(nodes: &[Node], format: &str) -> Result<String, Failure> {
    match format {
        "mathml" => mathml::render_formula(nodes).map_err(parse_family),
        "text" => Ok(nodes
            .iter()
            .map(Node::to_string)
            .collect::<Vec<_>>()
            .join(" ")),
        "svg" => Ok(svg::render_linear_formula(nodes)),
        "structural-svg" => Ok(svg::render_formula(nodes)),
        "model-json" => wire::to_model_json(nodes).map_err(parse_family),
        "cml" => cml::to_cml(nodes).map_err(parse_family),
        "smiles" => smiles::write_smiles(nodes).map_err(parse_family),
        "molfile" => {
            let molecule = match nodes.first() {
                Some(Node::Molecule(molecule)) => molecule,
                _ => {
                    return Err(failure(
                        1,
                        "Parse error: ",
                        "molfile output needs a molecule",
                    ));
                }
            };
            molfile::write_molfile(molecule, "").map_err(parse_family)
        }
        "html" => html::render_formula(nodes).map_err(parse_family),
        "latex" => latex::render_formula(nodes).map_err(parse_family),
        other => Err(failure(
            2,
            "Format error: ",
            format!("unknown formatter :{other}"),
        )),
    }
}

// -- parse-cml / roundtrip ----------------------------------------------

fn cmd_parse_cml(args: &ParseCmlArgs, out: &mut dyn Write) -> Result<(), Failure> {
    let nodes =
        cml::parse_cml(&args.input).map_err(|error| failure(1, "CML parse error: ", error))?;
    let text = nodes
        .iter()
        .map(Node::to_string)
        .collect::<Vec<_>>()
        .join(" ");
    emit(out, &text);
    Ok(())
}

fn cmd_roundtrip(
    args: &RoundtripArgs,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<(), Failure> {
    let rendered =
        crate::parse_to_text(&args.input).map_err(|error| failure(1, "Parse error: ", error))?;
    if rendered == args.input {
        emit(out, &rendered);
        Ok(())
    } else {
        let _ = writeln!(
            err,
            "round-trip mismatch:\n  input:    {:?}\n  rendered: {:?}",
            args.input, rendered
        );
        Err(silent(1))
    }
}

// -- lint ---------------------------------------------------------------

fn cmd_lint(args: &LintArgs, out: &mut dyn Write) -> Result<(), Failure> {
    let nodes =
        crate::parse_nodes(&args.input).map_err(|error| failure(1, "Parse error: ", error))?;
    let diagnostics = linter::run(&nodes);
    match args.format.as_str() {
        "json" => {
            let payload: Vec<serde_json::Value> = diagnostics
                .iter()
                .map(|diagnostic| {
                    serde_json::json!({
                        "severity": diagnostic.severity.as_str(),
                        "message": diagnostic.message,
                        "node": diagnostic.label,
                    })
                })
                .collect();
            let text = serde_json::to_string_pretty(&payload)
                .map_err(|error| failure(1, "Parse error: ", error))?;
            emit(out, &text);
        }
        _ => {
            for diagnostic in &diagnostics {
                emit(out, &diagnostic.to_string());
            }
        }
    }
    if diagnostics
        .iter()
        .any(|diagnostic| diagnostic.severity == linter::Severity::Error)
    {
        Err(silent(1))
    } else {
        Ok(())
    }
}

// -- validate -----------------------------------------------------------

fn cmd_validate(args: &ValidateArgs, out: &mut dyn Write) -> Result<(), Failure> {
    let nodes =
        crate::parse_nodes(&args.input).map_err(|error| failure(1, "Parse error: ", error))?;
    let annotations: Vec<&Identifier> = nodes
        .iter()
        .filter_map(|node| match node {
            Node::Molecule(molecule) => Some(&molecule.identifiers),
            _ => None,
        })
        .flatten()
        .collect();
    if annotations.is_empty() {
        emit(out, "no identifier annotations found");
        return Ok(());
    }
    let mut invalid = false;
    for identifier in &annotations {
        let known = KNOWN_CONVENTIONS.contains(&identifier.convention.as_str());
        let valid =
            known && identifiers::diagnostic(&identifier.convention, &identifier.value).is_none();
        let status = if !known {
            "unknown convention"
        } else if valid {
            "ok"
        } else {
            "INVALID"
        };
        emit(
            out,
            &format!(
                "{:<12} {:<40} {}",
                identifier.convention, identifier.value, status
            ),
        );
        if known && !valid {
            invalid = true;
        }
    }
    if invalid {
        Err(silent(1))
    } else {
        Ok(())
    }
}

// -- identity -----------------------------------------------------------

fn cmd_identity(args: &IdentityArgs, out: &mut dyn Write) -> Result<(), Failure> {
    if let Some(bin) = &args.engine_bin {
        return Err(failure(
            5,
            "Identity error: ",
            format!(
                "the Rust engine derives identities in-process through chemomer; \
                 --engine-bin ({bin}) is not supported"
            ),
        ));
    }
    // The reference parses whatever read_source gives — no input means
    // the empty string, whose parse failure is an identity error (5).
    let source = match source_text(args.input.as_deref(), args.file.as_deref()) {
        Ok(text) => text,
        Err(SourceError::Missing) => String::new(),
        Err(SourceError::Io(error)) => return Err(failure(1, "Error: ", error)),
    };
    let nodes =
        ingest(&source, &args.from).map_err(|error| failure(5, "Identity error: ", error))?;
    let molecules = molecules_in(&nodes);
    if molecules.len() != 1 {
        return Err(failure(
            5,
            "Identity error: ",
            format!(
                "identity derivation needs exactly one molecule, got {}",
                molecules.len()
            ),
        ));
    }
    let derived = identity::derive(molecules[0])
        .map_err(|error| failure(5, "Identity error: ", error.message()))?;
    emit(out, &derived.inchi);
    emit(out, &derived.inchikey);
    Ok(())
}

// -- resolve / cite -----------------------------------------------------

/// The selected identifier: the first set of
/// --cas/--name/--cid/--inchikey/--smiles wins (the reference's
/// filter_map order); `--cid` answers under the `pubchem-cid`
/// convention.
fn selected_identifier<'a>(
    cas: &'a Option<String>,
    name: &'a Option<String>,
    cid: &'a Option<String>,
    inchikey: &'a Option<String>,
    smiles: &'a Option<String>,
) -> Option<(&'static str, &'a str)> {
    if let Some(value) = cas {
        return Some(("cas", value));
    }
    if let Some(value) = name {
        return Some(("name", value));
    }
    if let Some(value) = cid {
        return Some(("pubchem-cid", value));
    }
    if let Some(value) = inchikey {
        return Some(("inchikey", value));
    }
    if let Some(value) = smiles {
        return Some(("smiles", value));
    }
    None
}

/// The `unknown resolver source` lookup errors immediately (the
/// reference's `Resolver[source]` raises before any request); the
/// per-source isolation of [`Registry::resolve_all`] only ever
/// swallows a source that was successfully looked up.
fn checked_adapter(
    registry: &Registry,
    source: &str,
    prefix: &str,
    code: i32,
) -> Result<(), Failure> {
    let _: Box<dyn Adapter> = registry
        .adapter(source)
        .map_err(|error: ResolverError| failure(code, prefix, error))?;
    Ok(())
}

fn cmd_resolve(args: &ResolveArgs, deps: &Deps<'_>, out: &mut dyn Write) -> Result<(), Failure> {
    let prefix = "Resolve error: ";
    let Some((convention, value)) = selected_identifier(
        &args.cas,
        &args.name,
        &args.cid,
        &args.inchikey,
        &args.smiles,
    ) else {
        return Err(failure(
            3,
            prefix,
            "give one of --cas/--name/--cid/--inchikey/--smiles",
        ));
    };
    checked_adapter(deps.registry, &args.source, prefix, 3)?;
    let sources = [args.source.as_str()];
    let substance = deps
        .registry
        .resolve(
            value,
            convention,
            Some(&sources),
            deps.fetch,
            deps.cache,
            args.refresh,
            // The live CLI path stamps the system clock.
            None,
        )
        .map_err(|error| failure(3, prefix, error))?
        .ok_or_else(|| {
            failure(
                3,
                prefix,
                format!("{} does not know {:?}", args.source, value),
            )
        })?;
    let text = match args.format.as_str() {
        "text" => substance.preferred_name.unwrap_or_default(),
        "smiles" => substance
            .identifier_value("canonical-smiles")
            .unwrap_or_default()
            .to_string(),
        // The reference's else-branch: anything but text/smiles is
        // the model JSON.
        _ => substance
            .to_model_json()
            .map_err(|error| failure(3, prefix, error))?,
    };
    emit(out, &text);
    Ok(())
}

fn cmd_cite(args: &CiteArgs, deps: &Deps<'_>, out: &mut dyn Write) -> Result<(), Failure> {
    let prefix = "Cite error: ";
    let Some((convention, value)) = selected_identifier(
        &args.cas,
        &args.name,
        &args.cid,
        &args.inchikey,
        &args.smiles,
    ) else {
        return Err(failure(
            4,
            prefix,
            "give one of --cas/--name/--cid/--inchikey/--smiles",
        ));
    };
    checked_adapter(deps.registry, &args.source, prefix, 4)?;
    let sources = [args.source.as_str()];
    let substance = deps
        .registry
        .resolve(
            value,
            convention,
            Some(&sources),
            deps.fetch,
            deps.cache,
            args.refresh,
            None,
        )
        .map_err(|error| failure(4, prefix, error))?
        .ok_or_else(|| {
            failure(
                4,
                prefix,
                format!("{} does not know {:?}", args.source, value),
            )
        })?;
    let xml = citation::to_xml(&substance).map_err(|error| failure(4, prefix, error))?;
    emit(out, &xml);
    Ok(())
}

// -- tree ---------------------------------------------------------------

fn cmd_tree(out: &mut dyn Write) -> Result<(), Failure> {
    // Thor 1.5's built-in `tree` task, mirrored: the root label and
    // the connector shape are the reference's; the entries are this
    // binary's own help texts (help and tree excluded, as Thor's tree
    // lists). The entries come from the clap command itself — one
    // source of truth for names and descriptions.
    let mut text = String::from("ascii_chem:cli\n");
    let command = Cli::command();
    let subcommands: Vec<&clap::Command> = command
        .get_subcommands()
        .filter(|sub| sub.get_name() != "help" && sub.get_name() != "tree")
        .collect();
    let last = subcommands.len().saturating_sub(1);
    for (index, sub) in subcommands.iter().enumerate() {
        let connector = if index == last { "└─" } else { "├─" };
        let about = sub
            .get_about()
            .map(|about| about.to_string())
            .unwrap_or_default();
        text.push_str(&format!("  {connector} {} ({about})\n", sub.get_name()));
    }
    emit(out, &text);
    Ok(())
}
