//! Parse + wire round-trip performance: the authoritative runtime's
//! numbers, recorded against the cross-runtime comparison in the Ruby
//! gem's benchmarks/RESULTS.md (TS peggy ~32 µs/input, PY rd ~672 µs,
//! Ruby parslet ~2,960 µs — all per parse input, 10-input batches).

use criterion::{black_box, criterion_group, criterion_main, Criterion};

/// The shared cross-runtime workload (identical to the TS bench and
/// the Ruby benchmarks/README.md contract).
const INPUTS: [&str; 10] = [
    "H_2O",
    "Ca^2+",
    "SO_4^2-",
    "(R)-CH_3CH(OH)COOH",
    "2H_2 + O_2 -> 2H_2O",
    "N_2 + 3H_2 <=>[Fe][400C] 2NH_3",
    "C1-C-C-C-C-C1",
    "CH_3-CH_2-OH",
    "^14C @name(\"carbon-14\") @cas(\"14104-86-4\")",
    "A ->[heat] B ->[cool] C",
];

fn bench_parse(c: &mut Criterion) {
    c.bench_function("parse x10", |b| {
        b.iter(|| {
            for _ in 0..10 {
                for input in INPUTS {
                    let _ = black_box(asciichem::parse_nodes(black_box(input)).expect("parse"));
                }
            }
        })
    });
}

fn bench_parse_to_text(c: &mut Criterion) {
    c.bench_function("parse+text x10", |b| {
        b.iter(|| {
            for _ in 0..10 {
                for input in INPUTS {
                    let _ = black_box(asciichem::parse_to_text(black_box(input)).expect("parse"));
                }
            }
        })
    });
}

criterion_group!(benches, bench_parse, bench_parse_to_text);
criterion_main!(benches);
