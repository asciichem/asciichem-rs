//! Runtime library path for the leptris engine (the citation track's
//! serialization layer). The `leptris` binding links the engine at
//! build time; its `-rpath` link argument stays on its own artifacts,
//! so THIS crate re-emits the runtime search path for its own
//! binaries (tests included) from the same `LEPTRIS_LIB_PATH` — the
//! test suite loads the exact library it linked against, without a
//! `DYLD_LIBRARY_PATH`/`LD_LIBRARY_PATH` crutch.

use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-env-changed=LEPTRIS_LIB_PATH");
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        // The DLL rides PATH at load time; no rpath equivalent needed.
        return;
    }
    let Ok(path) = std::env::var("LEPTRIS_LIB_PATH") else {
        return;
    };
    let mut candidate = PathBuf::from(&path);
    if !candidate.exists() {
        // Tolerate repo-relative paths: cargo runs build.rs from the
        // crate directory, so a path like `engine/build-shared/src`
        // (natural when set from the repository root) would not
        // resolve.
        if let Some(manifest) = std::env::var_os("CARGO_MANIFEST_DIR") {
            let repo_relative = Path::new(&manifest)
                .parent()
                .map(|root| root.join(&candidate));
            if let Some(found) = repo_relative.filter(|found| found.exists()) {
                candidate = found;
            }
        }
    }
    let directory = if candidate.is_file() {
        candidate
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or(candidate.clone())
    } else {
        candidate.clone()
    };
    // Absolute for -rpath: a relative rpath resolves against the test
    // binary's directory, not the crate.
    let directory = directory.canonicalize().unwrap_or(directory);
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", directory.display());
}
