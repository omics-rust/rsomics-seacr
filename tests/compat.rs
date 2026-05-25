//! Black-box compatibility against the reference SEACR_1.3 binary.
//!
//! Asserts byte-exact output for all six mode combinations:
//!   - `--control igg.bedgraph --norm non --mode stringent`
//!   - `--control igg.bedgraph --norm non --mode relaxed`
//!   - `--control igg.bedgraph --norm norm --mode stringent`
//!   - `--control igg.bedgraph --norm norm --mode relaxed`
//!   - `--fraction 0.3 --mode stringent`
//!   - `--fraction 0.3 --mode relaxed`
//!
//! Requires `SEACR_1.3.sh` on `$SEACR_SH` or at known fallback paths, plus
//! `Rscript` and `bedtools`. Missing any of these causes a loud self-skip.

use std::path::{Path, PathBuf};
use std::process::Command;

fn which(bin: &str) -> bool {
    Command::new(bin)
        .arg("--version")
        .output()
        .map(|o| o.status.success() || !o.stdout.is_empty() || !o.stderr.is_empty())
        .unwrap_or(false)
}

fn seacr_script() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("SEACR_SH") {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let candidates = [
        "SEACR_1.3.sh",
        "/tmp/seacr/SEACR_1.3.sh",
        "/Volumes/KIOXIA/rsomics-scratch/SEACR_1.3.sh",
    ];
    for c in candidates {
        let p = PathBuf::from(c);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

fn bin_path() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-seacr"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run_seacr(script: &Path, exp: &Path, ctrl: &str, norm: &str, mode: &str) -> String {
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("out");
    let status = Command::new("bash")
        .arg(script)
        .arg(exp)
        .arg(ctrl)
        .arg(norm)
        .arg(mode)
        .arg(&prefix)
        .env("LC_ALL", "C")
        .env("LANG", "C")
        .status()
        .expect("failed to launch SEACR");
    assert!(
        status.success(),
        "SEACR exited non-zero (ctrl={ctrl} norm={norm} mode={mode})"
    );
    let out_path = prefix.with_extension(format!("{mode}.bed"));
    std::fs::read_to_string(&out_path)
        .unwrap_or_else(|e| panic!("reading SEACR output {}: {e}", out_path.display()))
}

fn run_ours_control(exp: &Path, igg: &Path, norm: &str, mode: &str) -> String {
    let out = Command::new(bin_path())
        .arg(exp)
        .arg("--control")
        .arg(igg)
        .arg("--norm")
        .arg(norm)
        .arg("--mode")
        .arg(mode)
        .arg("-q")
        .output()
        .expect("failed to launch rsomics-seacr");
    assert!(
        out.status.success(),
        "rsomics-seacr exited non-zero (ctrl norm={norm} mode={mode}): {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn run_ours_fraction(exp: &Path, frac: &str, mode: &str) -> String {
    let out = Command::new(bin_path())
        .arg(exp)
        .arg("--fraction")
        .arg(frac)
        .arg("--mode")
        .arg(mode)
        .arg("-q")
        .output()
        .expect("failed to launch rsomics-seacr");
    assert!(
        out.status.success(),
        "rsomics-seacr exited non-zero (frac={frac} mode={mode}): {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn byte_exact_vs_seacr_all_modes() {
    let Some(script) = seacr_script() else {
        eprintln!("SKIP: SEACR_1.3.sh not found (set $SEACR_SH); cannot run compat oracle");
        return;
    };
    if !which("Rscript") || !which("bedtools") {
        eprintln!("SKIP: Rscript and/or bedtools missing; SEACR oracle cannot run");
        return;
    }

    let exp = golden("exp.bedgraph");
    let igg = golden("igg.bedgraph");

    // Control-based combos: --norm {non,norm} × --mode {stringent,relaxed}
    for norm in ["non", "norm"] {
        for mode in ["stringent", "relaxed"] {
            let oracle = run_seacr(&script, &exp, igg.to_str().unwrap(), norm, mode);
            let ours = run_ours_control(&exp, &igg, norm, mode);
            assert_eq!(
                ours, oracle,
                "byte mismatch: --control igg --norm {norm} --mode {mode}\n\
                 --- ours ---\n{ours}\n--- oracle ---\n{oracle}"
            );
            let n_peaks = oracle.lines().filter(|l| !l.is_empty()).count();
            eprintln!("--control --norm {norm} --mode {mode}: {n_peaks} peaks, byte-exact");
        }
    }

    // Numeric fraction combos: --fraction 0.3 × --mode {stringent,relaxed}
    let frac = "0.3";
    for mode in ["stringent", "relaxed"] {
        let oracle = run_seacr(&script, &exp, frac, "non", mode);
        let ours = run_ours_fraction(&exp, frac, mode);
        assert_eq!(
            ours, oracle,
            "byte mismatch: --fraction {frac} --mode {mode}\n\
             --- ours ---\n{ours}\n--- oracle ---\n{oracle}"
        );
        let n_peaks = oracle.lines().filter(|l| !l.is_empty()).count();
        eprintln!("--fraction {frac} --mode {mode}: {n_peaks} peaks, byte-exact");
    }
}
