use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_seacr(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-seacr");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let exp = manifest.join("tests/golden/exp.bedgraph");
    let igg = manifest.join("tests/golden/igg.bedgraph");
    let out = tempfile::NamedTempFile::new().unwrap();

    c.bench_function("rsomics-seacr golden (control non stringent)", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    exp.to_str().unwrap(),
                    "--control",
                    igg.to_str().unwrap(),
                    "--norm",
                    "non",
                    "--mode",
                    "stringent",
                    "--output",
                    out.path().to_str().unwrap(),
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_seacr);
criterion_main!(benches);
