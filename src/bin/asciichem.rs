//! The `asciichem` binary — a thin main over [asciichem::cli::run].
//! The real network path: the shipped resolver sources, the user
//! cache dir, ureq's NetFetch. All command behaviour (stdout shapes,
//! stderr prefixes, exit codes) lives in [asciichem::cli] so the
//! contract corpus can drive the same dispatch offline.

use asciichem::cli::{self, Deps};
use asciichem::resolver::adapter::NetFetch;
use asciichem::resolver::{Cache, Registry};

fn main() {
    static SHIPPED: std::sync::OnceLock<Registry> = std::sync::OnceLock::new();
    let cache = Cache::new();
    let deps = Deps {
        fetch: &NetFetch::default(),
        cache: &cache,
        registry: SHIPPED.get_or_init(Registry::shipped),
    };
    let code = cli::run(
        std::env::args_os(),
        &deps,
        &mut std::io::stdout().lock(),
        &mut std::io::stderr().lock(),
    );
    std::process::exit(code);
}
