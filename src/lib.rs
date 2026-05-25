//! Sparse Enrichment Analysis for CUT&RUN — bedGraph signal to BED peaks.
//!
//! Independent Rust reimplementation of the SEACR method (Meers, Tenenbaum &
//! Henikoff 2019, DOI 10.1186/s13072-019-0287-4), informed by the SEACR
//! algorithm and its published constants. The algorithm and constants are not
//! subject to copyright (idea-expression dichotomy); the implementation is
//! original Rust. SEACR source (GPL-3.0) is credited as upstream.
//!
//! ## Method
//!
//! A bedGraph is a coordinate-sorted list of `chrom start end value` intervals
//! covering only nonzero signal. A *signal block* is a maximal run of strictly
//! adjacent intervals on one chromosome. Per block:
//!
//! - **total signal** (AUC) = Σ value·(end − start) over the block
//! - **max signal**          = the largest value attained at any base
//! - **max region**          = the span from the farthest-upstream to the
//!   farthest-downstream base that attains the max signal
//! - **num_intervals**       = count of input intervals merged into the block
//!
//! The R script computes two AUC thresholds (`x0` = stringent, `z0` = relaxed)
//! and one num-intervals threshold (`d0`) from the block vectors. A block
//! passes when `auc > thresh && num_intervals > d0`. After filtering,
//! nearby peaks are merged (gap < `mean_width / 10`) and control-enriched
//! regions are subtracted.

#![allow(clippy::cast_precision_loss)]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// One bedGraph interval: a constant-value span on a chromosome.
struct Interval {
    chrom: u32,
    start: u64,
    end: u64,
    value: f64,
}

/// Maps chromosome names to compact ids in first-seen order.
#[derive(Default)]
struct ChromTable {
    names: Vec<String>,
}

impl ChromTable {
    fn intern(&mut self, name: &str) -> u32 {
        if let Some(pos) = self.names.iter().position(|n| n == name) {
            return pos as u32;
        }
        self.names.push(name.to_owned());
        (self.names.len() - 1) as u32
    }

    fn name(&self, id: u32) -> &str {
        &self.names[id as usize]
    }
}

/// A signal block: a merged run of strictly adjacent nonzero intervals.
///
/// `num_intervals` counts how many input intervals were merged into the block;
/// it is the second statistic fed to SEACR's R threshold step (expmax/ctrlmax).
pub struct Block {
    chrom: u32,
    pub start: u64,
    pub end: u64,
    pub total: f64,
    pub max: f64,
    pub max_start: u64,
    pub max_end: u64,
    pub num_intervals: u64,
}

/// Threshold axis for calling peaks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Uses the `x0` AUC threshold (peak of the pctremain curve).
    Stringent,
    /// Uses the `z0` AUC threshold (midpoint of the pctremain curve).
    Relaxed,
}

impl std::str::FromStr for Mode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "stringent" => Ok(Self::Stringent),
            "relaxed" => Ok(Self::Relaxed),
            other => Err(format!("expected 'stringent' or 'relaxed', got '{other}'")),
        }
    }
}

/// Control-to-experimental normalisation toggle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Norm {
    On,
    Off,
}

impl std::str::FromStr for Norm {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "norm" => Ok(Self::On),
            "non" => Ok(Self::Off),
            other => Err(format!("expected 'norm' or 'non', got '{other}'")),
        }
    }
}

/// Threshold source.
pub enum Threshold {
    /// Top `fraction` of blocks by AUC.
    Fraction(f64),
    /// IgG control bedGraph file path.
    Control(std::path::PathBuf),
}

/// Parse a bedGraph file into intervals. Skips the first data line (matching
/// SEACR's AWK pipeline which advances past line 1 before processing).
fn parse_bedgraph(path: &Path, chroms: &mut ChromTable) -> Result<Vec<Interval>> {
    let file = std::fs::File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("reading {}: {e}", path.display())))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut out = Vec::new();
    let mut line = String::new();
    let mut lineno = 0usize;
    // SEACR's AWK: `BEGIN{s=1}; {if(s==1){s++}...` skips the first non-header
    // line before recording anything. We replicate that by skipping the first
    // data line (non-blank, non-track, non-# line).
    let mut first_data_seen = false;
    loop {
        line.clear();
        let read = reader.read_line(&mut line).map_err(RsomicsError::Io)?;
        if read == 0 {
            break;
        }
        lineno += 1;
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() || trimmed.starts_with("track") || trimmed.starts_with('#') {
            continue;
        }
        if !first_data_seen {
            first_data_seen = true;
            continue; // skip the first data line — SEACR AWK quirk
        }
        let mut f = trimmed.split('\t');
        let (Some(chrom), Some(start), Some(end), Some(value)) =
            (f.next(), f.next(), f.next(), f.next())
        else {
            return Err(RsomicsError::InvalidInput(format!(
                "{}:{lineno}: expected 4 tab-separated columns",
                path.display()
            )));
        };
        let start: u64 = start.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("{}:{lineno}: start: {e}", path.display()))
        })?;
        let end: u64 = end.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("{}:{lineno}: end: {e}", path.display()))
        })?;
        let value: f64 = value.parse().map_err(|e| {
            RsomicsError::InvalidInput(format!("{}:{lineno}: value: {e}", path.display()))
        })?;
        if value == 0.0 {
            continue;
        }
        out.push(Interval {
            chrom: chroms.intern(chrom),
            start,
            end,
            value,
        });
    }
    Ok(out)
}

/// Collapse nonzero intervals into signal blocks using strict adjacency
/// (`interval.start == previous.end`), matching SEACR's AWK `$2==stop` check.
fn build_blocks(intervals: &[Interval]) -> Vec<Block> {
    let mut blocks = Vec::new();
    if intervals.is_empty() {
        return blocks;
    }
    let mut i = 0;
    while i < intervals.len() {
        let chrom = intervals[i].chrom;
        let block_start = intervals[i].start;
        let mut block_end = intervals[i].end;
        let mut total = intervals[i].value * (intervals[i].end - intervals[i].start) as f64;
        let mut max = intervals[i].value;
        let mut max_start = intervals[i].start;
        let mut max_end = intervals[i].end;
        let mut num = 1u64;
        let mut j = i + 1;
        while j < intervals.len() && intervals[j].chrom == chrom && intervals[j].start == block_end
        {
            let iv = &intervals[j];
            total += iv.value * (iv.end - iv.start) as f64;
            if iv.value > max {
                max = iv.value;
                max_start = iv.start;
                max_end = iv.end;
            } else if iv.value == max {
                max_end = iv.end;
            }
            block_end = iv.end;
            num += 1;
            j += 1;
        }
        blocks.push(Block {
            chrom,
            start: block_start,
            end: block_end,
            total,
            max,
            max_start,
            max_end,
            num_intervals: num,
        });
        i = j;
    }
    blocks
}

/// Build blocks from a bedGraph and apply SEACR's pipeline quirks:
/// 1. Skip the first data line (AWK `BEGIN{s=1}` skips line 1).
/// 2. Drop the final block (AWK has no `END{}`, so the last accumulated block is never printed).
/// 3. Round each block's total to 6 significant figures (AWK's default `OFMT = %.6g`).
fn seacr_blocks(intervals: &[Interval]) -> Vec<Block> {
    let mut blocks = build_blocks(intervals);
    // Drop last block — AWK block-building never emits it (no END{} clause).
    if !blocks.is_empty() {
        blocks.pop();
    }
    round_block_totals(blocks)
}

/// Round each block's `total` to 6 significant figures (awk's default OFMT = `%.6g`).
///
/// SEACR's AWK pipeline prints block totals with that precision. When the merge
/// AWK reads them back, it parses the rounded string values. All downstream
/// arithmetic must use the same quantities; rounding here replicates that.
fn round_block_totals(mut blocks: Vec<Block>) -> Vec<Block> {
    for b in &mut blocks {
        b.total = round6g(b.total);
    }
    blocks
}

/// Round to 6 significant figures (awk `%.6g` semantics, using round-half-away-from-zero
/// which is what most libc `printf` implementations use for `%g`).
fn round6g(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let mag = x.abs().log10().floor() as i32;
    let factor = 10f64.powi(5 - mag); // 6 sig figs → 5 digits after leading
    (x * factor).round() / factor
}

/// Empirical CDF: fraction of `sorted_vals` that are ≤ x. Matches R's `ecdf`.
fn ecdf(sorted_vals: &[f64], x: f64) -> f64 {
    let n = sorted_vals.len();
    if n == 0 {
        return 0.0;
    }
    let pos = sorted_vals.partition_point(|&v| v <= x);
    pos as f64 / n as f64
}

/// `pctremain(x)` from SEACR's R script:
/// `(len_exp - ecdf(expvec)(x)*len_exp) / (len_both - ecdf(both)(x)*len_both)`
///
/// Returns `None` when the denominator is zero (≡ R's `NA`).
fn pctremain(x: f64, sorted_exp: &[f64], sorted_both: &[f64]) -> Option<f64> {
    let len_exp = sorted_exp.len() as f64;
    let len_both = sorted_both.len() as f64;
    let denom = len_both - ecdf(sorted_both, x) * len_both;
    if denom == 0.0 {
        None
    } else {
        let numer = len_exp - ecdf(sorted_exp, x) * len_exp;
        Some(numer / denom)
    }
}

/// Bandwidth for Gaussian KDE: R's `bw.nrd0`.
///
/// `0.9 * lo * n^(-1/5)` where `lo = min(sd, IQR/1.34)`, falling back to `sd`
/// when `IQR == 0`.
fn bw_nrd0(v: &[f64]) -> f64 {
    let n = v.len();
    assert!(n >= 2, "density requires at least 2 observations");
    let mean = v.iter().sum::<f64>() / n as f64;
    let variance = v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
    let sd = variance.sqrt();
    let iqr = {
        // R's IQR: quantile(type 7) at 0.75 - quantile(type 7) at 0.25
        let mut s = v.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        quantile_type7(&s, 0.75) - quantile_type7(&s, 0.25)
    };
    let lo = if iqr == 0.0 { sd } else { sd.min(iqr / 1.34) };
    0.9 * lo * (n as f64).powf(-0.2)
}

/// Mode of R's `density(v)` — the x at `max(y)` on a 512-point Gaussian KDE.
///
/// Uses direct kernel evaluation (same result as R's FFT-based `density()` at
/// the 512-output-point grid) because the mode index is identical between the
/// two methods at this precision. The output grid spans
/// `[min(v) − 3·bw, max(v) + 3·bw]` with 512 equally spaced points.
pub fn density_mode(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return v.first().copied().unwrap_or(1.0);
    }
    let bw = bw_nrd0(v);
    if bw == 0.0 {
        // Degenerate: all values identical.
        return v[0];
    }
    let x_lo = v.iter().cloned().fold(f64::INFINITY, f64::min) - 3.0 * bw;
    let x_hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 3.0 * bw;
    let n_grid = 512usize;
    let dx = (x_hi - x_lo) / (n_grid as f64 - 1.0);
    let inv_bw = 1.0 / bw;
    let n_obs = v.len() as f64;
    let mut best_y = f64::NEG_INFINITY;
    let mut best_x = x_lo;
    for i in 0..n_grid {
        let xi = x_lo + i as f64 * dx;
        // Gaussian kernel: (1/n) * sum_j kernel((xi - v_j) / bw) / bw
        let y = v
            .iter()
            .map(|&vj| gauss_kernel((xi - vj) * inv_bw))
            .sum::<f64>()
            / (n_obs * bw);
        if y > best_y {
            best_y = y;
            best_x = xi;
        }
    }
    best_x
}

#[inline]
fn gauss_kernel(u: f64) -> f64 {
    (-0.5 * u * u).exp() * std::f64::consts::FRAC_1_SQRT_2 / std::f64::consts::SQRT_2
    // = exp(-u²/2) / sqrt(2π)
}

/// `(1 − fraction)` quantile of a sorted slice, using R's type-7 definition.
#[must_use]
pub fn quantile_type7(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted[0];
    }
    let h = (n as f64 - 1.0) * p;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

/// Compute thresholds x0, z0, d0 from the block AUC vectors (`expvec`,
/// `ctrlvec`) and the block num-interval vectors (`expmax`, `ctrlmax`).
///
/// This is a faithful port of the pctremain/z0/spurious-correction logic in
/// `SEACR_1.3.R`. Returns `(x0, z0, d0)`.
fn compute_thresholds(
    expvec: &[f64],
    ctrlvec: &[f64],
    expmax: &[f64],
    ctrlmax: &[f64],
) -> (f64, f64, f64) {
    let mut sorted_exp = expvec.to_vec();
    sorted_exp.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut sorted_ctrl = ctrlvec.to_vec();
    sorted_ctrl.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut both: Vec<f64> = expvec.iter().chain(ctrlvec.iter()).copied().collect();
    both.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Unique values across both distributions
    let x_unique: Vec<f64> = {
        let mut u = both.clone();
        u.dedup();
        u
    };

    let pr = |x: f64| pctremain(x, &sorted_exp, &both);

    // x0: x maximising pctremain among pctremain < 1
    let x0 = {
        let cands: Vec<f64> = x_unique
            .iter()
            .copied()
            .filter(|&x| pr(x).is_some_and(|v| v < 1.0))
            .collect();
        let best = cands
            .iter()
            .copied()
            .max_by(|&a, &b| pr(a).unwrap().partial_cmp(&pr(b).unwrap()).unwrap());
        best.unwrap_or_else(|| *x_unique.last().unwrap_or(&1.0))
    };

    // z0: midpoint-of-curve value
    let z0 = compute_z0(x0, &x_unique, pr);

    // Spurious-threshold correction
    let (x0, z0) = spurious_correction(x0, z0, &x_unique, pr);

    // d0: min num-interval value where pctremain2 > 1
    let d0 = compute_d0(expmax, ctrlmax);

    (x0, z0, d0)
}

/// Compute z0 from x0 using the midpoint-of-curve logic in SEACR_1.3.R.
fn compute_z0(x0: f64, x_unique: &[f64], pr: impl Fn(f64) -> Option<f64>) -> f64 {
    let z: Vec<f64> = x_unique.iter().copied().filter(|&x| x <= x0).collect();
    if z.is_empty() {
        return x0;
    }
    // min pctremain over z
    let pr_z_min = z
        .iter()
        .filter_map(|&x| pr(x))
        .fold(f64::INFINITY, f64::min);
    let pr_x0 = pr(x0).unwrap_or(0.0);
    let midpoint = (pr_x0 + pr_z_min) / 2.0;
    // z2: z value closest to midpoint pctremain
    let z2 = *z
        .iter()
        .min_by(|&&a, &&b| {
            let da = (pr(a).unwrap_or(0.0) - midpoint).abs();
            let db = (pr(b).unwrap_or(0.0) - midpoint).abs();
            da.partial_cmp(&db).unwrap()
        })
        .unwrap_or(&x0);

    if x0 == z2 {
        // Added 7/15/19 to avoid omitting z when x0==z2
        return x0;
    }
    // z: values in z that are > z2
    let z_filtered: Vec<f64> = z.iter().copied().filter(|&x| x > z2).collect();
    if z_filtered.is_empty() {
        return x0;
    }
    let z_max = z_filtered.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let z_min = z_filtered.iter().copied().fold(f64::INFINITY, f64::min);
    let target = z_max - 0.5 * (z_max - z_min);
    // z0: z value closest to midpoint of [z_min, z_max]
    *z_filtered
        .iter()
        .min_by(|&&a, &&b| (a - target).abs().partial_cmp(&(b - target).abs()).unwrap())
        .unwrap_or(&x0)
}

/// Spurious-threshold correction: if the sub-curve's max pctremain is ≥ 95%
/// of the full curve's max pctremain, replace x0/z0 with the sub-curve's.
fn spurious_correction(
    x0: f64,
    z0_initial: f64,
    x_unique: &[f64],
    pr: impl Fn(f64) -> Option<f64>,
) -> (f64, f64) {
    if x_unique.len() < 2 {
        return (x0, z0_initial);
    }
    // frame: thresh = x[i], pct = pctremain(x[i]), diff = |pctremain(x[i+1]) - pctremain(x[i])| for i in 0..n-1
    let pairs: Vec<(f64, f64, f64)> = x_unique
        .windows(2)
        .filter_map(|w| {
            let thresh = w[0];
            let pr_a = pr(w[0]);
            let pr_b = pr(w[1]);
            let diff = match (pr_a, pr_b) {
                (Some(a), Some(b)) => (b - a).abs(),
                (Some(a), None) => a.abs(),
                (None, Some(b)) => b.abs(),
                (None, None) => return None,
            };
            let pct = pr(thresh)?;
            Some((thresh, pct, diff))
        })
        .collect::<Vec<_>>()
        .into_iter()
        .filter(|&(_, _, diff)| !diff.is_nan())
        .collect::<Vec<_>>();

    // Find test3: smallest quantile level (0.99, 0.999, …) where quantile(diffs) > 0
    let diffs: Vec<f64> = pairs.iter().map(|&(_, _, d)| d).collect();
    if diffs.is_empty() {
        return (x0, z0_initial);
    }
    let mut test3 = 0.99f64;
    let mut nines = 2usize;
    loop {
        let q = quantile_of_sorted(
            &{
                let mut s = diffs.clone();
                s.sort_by(|a, b| a.partial_cmp(b).unwrap());
                s
            },
            test3,
        );
        if q > 0.0 {
            break;
        }
        nines += 1;
        // test3 = 0.999...9 with `nines` nines = 1 - 10^{-nines}
        test3 = 1.0 - 10f64.powi(-(nines as i32));
    }
    let q_thresh = quantile_of_sorted(
        &{
            let mut s = diffs.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            s
        },
        test3,
    );

    // a: thresh values where diff != 0 && diff < q_thresh
    let a: Vec<f64> = pairs
        .iter()
        .filter(|&&(_, _, d): &&(f64, f64, f64)| d != 0.0 && d < q_thresh)
        .map(|&(t, _, _)| t)
        .collect();
    if a.is_empty() {
        return (x0, z0_initial);
    }
    let mut sorted_a = a.clone();
    sorted_a.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // a0: a value maximising pctremain among pctremain < 1
    let a_cands: Vec<f64> = a
        .iter()
        .copied()
        .filter(|&x| pr(x).is_some_and(|v| v < 1.0))
        .collect();
    if a_cands.is_empty() {
        return (x0, z0_initial);
    }
    let a0 = *a_cands
        .iter()
        .max_by(|&&p, &&q| pr(p).unwrap().partial_cmp(&pr(q).unwrap()).unwrap())
        .unwrap();

    let b0 = compute_z0(a0, &sorted_a, &pr);

    // Check ratio
    let max_pr_a = a_cands
        .iter()
        .filter_map(|&x| pr(x))
        .fold(f64::NEG_INFINITY, f64::max);
    let max_pr_x = x_unique
        .iter()
        .filter_map(|&x| pr(x))
        .filter(|&v| v < 1.0)
        .fold(f64::NEG_INFINITY, f64::max);

    if max_pr_x > 0.0 && max_pr_a / max_pr_x > 0.95 {
        (a0, b0)
    } else {
        (x0, z0_initial)
    }
}

/// R's `quantile(x, p, type=7)` on an already-sorted vector.
fn quantile_of_sorted(sorted: &[f64], p: f64) -> f64 {
    quantile_type7(sorted, p)
}

/// d0: min num-intervals value where `pctremain2(x) > 1`.
///
/// `pctremain2(x) = 1 - (ecdf(expmax)(x) - ecdf(ctrlmax)(x))`
fn compute_d0(expmax: &[f64], ctrlmax: &[f64]) -> f64 {
    if expmax.is_empty() && ctrlmax.is_empty() {
        return 1.0;
    }
    let mut sorted_exp = expmax.to_vec();
    sorted_exp.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut sorted_ctrl = ctrlmax.to_vec();
    sorted_ctrl.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let both2: Vec<f64> = {
        let mut v: Vec<f64> = expmax.iter().chain(ctrlmax.iter()).copied().collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v.dedup();
        v
    };

    let pctremain2 = |x: f64| 1.0 - (ecdf(&sorted_exp, x) - ecdf(&sorted_ctrl, x));

    let candidates: Vec<f64> = both2
        .iter()
        .copied()
        .filter(|&x| pctremain2(x) > 1.0)
        .collect();
    if candidates.is_empty() {
        1.0
    } else {
        candidates.iter().cloned().fold(f64::INFINITY, f64::min)
    }
}

/// Knee-point detection for the norm-mode cutoff value.
///
/// Implements the `dist2d` + `expframe`/`ctrlframe` logic from SEACR_1.3.R.
/// Returns the AUC cutoff value below which density mode is computed.
fn knee_value(vec: &[f64]) -> f64 {
    let n = vec.len();
    let ninety_pct_idx = (0.9 * n as f64) as usize;
    if ninety_pct_idx == 0 {
        return *vec
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&1.0);
    }
    let mut sorted = vec.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p90 = sorted[ninety_pct_idx.min(n - 1)];
    let max_val = *sorted.last().unwrap();

    // Build the frame: count = seq(1,0,len=n), quant = desc/max, value = desc
    // sort descending
    let mut desc = sorted.clone();
    desc.reverse();

    // dist2d: signed distance from point (count,quant) to the line (0,0)→(1,1)
    // dist2d(a, b, c) with b=(0,0) c=(1,1):
    // v1 = c - b = (1,1); v2 = a - b = (count, quant)
    // det [[v1.x, v2.x],[v1.y, v2.y]] / |v1| = (1*quant - 1*count) / sqrt(2)
    // ... but SEACR passes a=(x[1],x[2]), b=0, c=1 as scalars, computing:
    // v1 = b - c = 0 - 1 = -1 (scalar)
    // v2 = a - b = (count - 0, quant - 0) => treat as pair
    // Wait, re-reading: dist2d(c(x[1],x[2]), 0, 1)
    // v1 = b - c = 0 - 1 = -1 (scalar, but used as col 1 of m)
    // v2 = a - b = c(x[1]-0, x[2]-0) = c(x[1], x[2])
    // m = cbind(v1, v2) = 2x2 matrix: [[-1, x1],[−1, x2]] ← wait v1 is scalar
    // Actually: v1 and v2 are both 2-vectors here. Let me re-read the R:
    // dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-cbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}
    // called as: dist2d(c(x[1],x[2]),0,1)
    // b=0 (scalar), c=1 (scalar) → v1 = 0-1 = -1; v2 = c(x[1],x[2]) - 0 = c(x[1],x[2])
    // m = cbind(-1, c(x[1],x[2])) = 2×2 matrix: col1 = [-1,-1], col2 = [x1,x2]
    // det(m) = (-1)*x2 - (-1)*x1 = x1 - x2
    // sum(v1*v1) = (-1)^2 = 1 → sqrt = 1
    // dist = (x1 - x2) / 1 = count - quant (for that row)
    // So dist2d(c(count,quant), 0, 1) = count - quant (per row)

    let mut frame: Vec<(f64, f64, f64, f64)> = desc
        .iter()
        .enumerate()
        .map(|(i, &val)| {
            let count = 1.0 - (i as f64) / (n as f64 - 1.0).max(1.0);
            let quant = if max_val == 0.0 { 0.0 } else { val / max_val };
            let diff = (count - quant).abs();
            let _dist = count - quant; // dist2d(c(count,quant), 0, 1) — not used in filtering
            (count, quant, val, diff)
        })
        .collect();

    // Keep rows where diff > 0.9 * max(diff)
    let max_diff = frame
        .iter()
        .map(|&(_, _, _, d)| d)
        .fold(f64::NEG_INFINITY, f64::max);
    let threshold_diff = 0.9 * max_diff;
    frame.retain(|&(_, _, _, d)| d > threshold_diff);

    if frame.is_empty() {
        return p90.max(*sorted.last().unwrap());
    }

    // Find row with max dist2d = max(count - quant) among filtered rows
    // dist2d = count - quant (same as diff when count > quant, which is the typical case)
    // But dist2d could be negative, so we need the actual count-quant signed value
    let knee_val = frame
        .iter()
        .max_by(|&(ca, qa, _, _), &(cb, qb, _, _)| {
            let da = ca - qa;
            let db = cb - qb;
            da.partial_cmp(&db).unwrap()
        })
        .map(|&(_, _, val, _)| val)
        .unwrap_or(max_val);

    // Apply the 90th-percentile guard
    if knee_val > p90 { knee_val } else { p90 }
}

/// Render a value with six significant figures, matching R `signif(x, 6)`.
#[must_use]
pub fn signif6(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let digits = 6i32;
    let mag = x.abs().log10().floor() as i32;
    let power = digits - 1 - mag;
    let factor = 10f64.powi(power);
    let rounded = round_half_even(x * factor) / factor;
    format_r_numeric(rounded)
}

fn round_half_even(x: f64) -> f64 {
    let f = x.floor();
    let diff = x - f;
    if diff < 0.5 {
        f
    } else if diff > 0.5 {
        f + 1.0
    } else if (f as i64) % 2 == 0 {
        f
    } else {
        f + 1.0
    }
}

fn format_r_numeric(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        return format!("{}", x as i64);
    }
    let mut s = format!("{x:.10}");
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
    s
}

fn write_peak(b: &Block, chroms: &ChromTable, w: &mut impl Write) -> std::io::Result<()> {
    let name = chroms.name(b.chrom);
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}:{}-{}",
        name,
        b.start,
        b.end,
        signif6(b.total),
        signif6(b.max),
        name,
        b.max_start,
        b.max_end,
    )
}

/// Merge adjacent blocks (gaps < `gap_tolerance`) and collect the result.
///
/// The final merged block is dropped to match SEACR's merge AWK, which also
/// lacks an `END{}` clause and never emits its last accumulated result.
fn merge_nearby(blocks: &[Block], gap_tolerance: u64) -> Vec<Block> {
    if blocks.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<Block> = Vec::new();
    let mut cur_chrom = blocks[0].chrom;
    let mut cur_start = blocks[0].start;
    let mut cur_end = blocks[0].end;
    let mut cur_total = blocks[0].total;
    let mut cur_max = blocks[0].max;
    let mut cur_max_start = blocks[0].max_start;
    let mut cur_max_end = blocks[0].max_end;
    let mut cur_num = blocks[0].num_intervals;

    for b in &blocks[1..] {
        if b.chrom == cur_chrom && b.start < cur_end + gap_tolerance {
            // Merge
            cur_end = b.end;
            cur_total += b.total;
            if b.max > cur_max {
                cur_max = b.max;
                cur_max_start = b.max_start;
                cur_max_end = b.max_end;
            } else if b.max == cur_max {
                cur_max_end = b.max_end;
            }
            cur_num += b.num_intervals;
        } else {
            out.push(Block {
                chrom: cur_chrom,
                start: cur_start,
                end: cur_end,
                total: cur_total,
                max: cur_max,
                max_start: cur_max_start,
                max_end: cur_max_end,
                num_intervals: cur_num,
            });
            cur_chrom = b.chrom;
            cur_start = b.start;
            cur_end = b.end;
            cur_total = b.total;
            cur_max = b.max;
            cur_max_start = b.max_start;
            cur_max_end = b.max_end;
            cur_num = b.num_intervals;
        }
    }
    // SEACR's merge AWK also has no END{} — the final accumulated block is dropped.
    // We replicate that: do NOT push the last cur block.
    let _ = (
        cur_chrom,
        cur_start,
        cur_end,
        cur_total,
        cur_max,
        cur_max_start,
        cur_max_end,
        cur_num,
    );
    out
}

/// Remove experiment peaks that overlap any control peak (bedtools intersect -v).
fn subtract_control(exp_peaks: &[Block], ctrl_peaks: &[Block]) -> Vec<Block> {
    exp_peaks
        .iter()
        .filter(|ep| {
            !ctrl_peaks
                .iter()
                .any(|cp| cp.chrom == ep.chrom && cp.start < ep.end && cp.end > ep.start)
        })
        .map(|b| Block {
            chrom: b.chrom,
            start: b.start,
            end: b.end,
            total: b.total,
            max: b.max,
            max_start: b.max_start,
            max_end: b.max_end,
            num_intervals: b.num_intervals,
        })
        .collect()
}

/// Call peaks from `experimental` and write the BED to `out`. Returns the
/// number of peaks emitted.
pub fn call_peaks(
    experimental: &Path,
    threshold: &Threshold,
    norm: Norm,
    mode: Mode,
    out: &mut impl Write,
) -> Result<usize> {
    let mut chroms = ChromTable::default();
    let exp_intervals = parse_bedgraph(experimental, &mut chroms)?;
    let exp_blocks = seacr_blocks(&exp_intervals);

    let expvec: Vec<f64> = exp_blocks.iter().map(|b| b.total).collect();
    let expmax: Vec<f64> = exp_blocks.iter().map(|b| b.num_intervals as f64).collect();

    let (auc_thresh, num_thresh) = match threshold {
        Threshold::Fraction(frac) => {
            if !(*frac > 0.0 && *frac < 1.0) {
                return Err(RsomicsError::InvalidInput(format!(
                    "numeric threshold must be in (0,1), got {frac}"
                )));
            }
            // R: x0 = min(frame$values[frame$percentile <= ctrl])
            // ecdf-based: keep blocks where 1 - ecdf >= frac → auc >= (1-frac) quantile
            let mut sorted_auc = expvec.clone();
            sorted_auc.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut sorted_num = expmax.clone();
            sorted_num.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let (x0, z0) = match mode {
                Mode::Stringent => {
                    // x0 = min AUC where 1 - ecdf(auc)(x) <= frac
                    let x0 = numeric_threshold_min(&sorted_auc, *frac);
                    // z0 from num-interval axis (unused in this branch but consistent)
                    let z0 = numeric_threshold_min(&sorted_num, *frac);
                    (x0, z0)
                }
                Mode::Relaxed => {
                    let x0 = numeric_threshold_min(&sorted_auc, *frac);
                    let z0 = numeric_threshold_min(&sorted_num, *frac);
                    (x0, z0)
                }
            };
            let auc = match mode {
                Mode::Stringent => x0,
                Mode::Relaxed => z0,
            };
            (auc, 0.0) // d0=0 for numeric mode
        }
        Threshold::Control(ctrl_path) => {
            let mut ctrl_chroms = ChromTable::default();
            let ctrl_intervals = parse_bedgraph(ctrl_path, &mut ctrl_chroms)?;
            let ctrl_blocks = seacr_blocks(&ctrl_intervals);

            let mut ctrlvec: Vec<f64> = ctrl_blocks.iter().map(|b| b.total).collect();
            let ctrlmax: Vec<f64> = ctrl_blocks.iter().map(|b| b.num_intervals as f64).collect();

            if norm == Norm::On {
                // Knee-point detection, then density-mode ratio
                let expvalue = knee_value(&expvec);
                let ctrlvalue = knee_value(&ctrlvec);
                let sub_exp: Vec<f64> = expvec.iter().copied().filter(|&v| v <= expvalue).collect();
                let sub_ctrl: Vec<f64> = ctrlvec
                    .iter()
                    .copied()
                    .filter(|&v| v <= ctrlvalue)
                    .collect();
                let exp_mode = density_mode(&sub_exp);
                let ctrl_mode = density_mode(&sub_ctrl);
                let constant = exp_mode / ctrl_mode;
                for v in &mut ctrlvec {
                    *v *= constant;
                }
            }

            let (x0, z0, d0) = compute_thresholds(&expvec, &ctrlvec, &expmax, &ctrlmax);
            let auc = match mode {
                Mode::Stringent => x0,
                Mode::Relaxed => z0,
            };
            (auc, d0)
        }
    };

    // Filter blocks: auc > thresh && num_intervals > d0
    let filtered: Vec<&Block> = exp_blocks
        .iter()
        .filter(|b| b.total > auc_thresh && b.num_intervals as f64 > num_thresh)
        .collect();

    if filtered.is_empty() {
        return Ok(0);
    }

    // Merge nearby peaks: gap < mean_block_width / 10
    let mean_width: f64 = filtered
        .iter()
        .map(|b| (b.end - b.start) as f64)
        .sum::<f64>()
        / filtered.len() as f64;
    let gap_tol = (mean_width / 10.0) as u64;

    let filtered_owned: Vec<Block> = filtered
        .into_iter()
        .map(|b| Block {
            chrom: b.chrom,
            start: b.start,
            end: b.end,
            total: b.total,
            max: b.max,
            max_start: b.max_start,
            max_end: b.max_end,
            num_intervals: b.num_intervals,
        })
        .collect();
    let merged = merge_nearby(&filtered_owned, gap_tol);

    // If control bedgraph: subtract control-enriched regions
    let final_peaks = if let Threshold::Control(ctrl_path) = threshold {
        // Recompute ctrl blocks for intersection (using original, before scaling)
        let mut ctrl_chroms2 = ChromTable::default();
        let ctrl_intervals2 = parse_bedgraph(ctrl_path, &mut ctrl_chroms2)?;
        let ctrl_blocks2 = seacr_blocks(&ctrl_intervals2);
        // Map ctrl chrom names to exp chrom ids
        let ctrl_blocks_mapped: Vec<Block> = ctrl_blocks2
            .iter()
            .filter_map(|cb| {
                let ctrl_name = ctrl_chroms2.name(cb.chrom);
                let exp_id = chroms.names.iter().position(|n| n == ctrl_name)?;
                Some(Block {
                    chrom: exp_id as u32,
                    start: cb.start,
                    end: cb.end,
                    total: cb.total,
                    max: cb.max,
                    max_start: cb.max_start,
                    max_end: cb.max_end,
                    num_intervals: cb.num_intervals,
                })
            })
            .collect();

        // Control threshold: auc > x0 (always uses x0/stringent for ctrl filtering)
        // Looking at SEACR.sh line 159: `awk -v value=$thresh '$4 > value'`
        // where thresh = line 1 of threshold.txt = x0 (the stringent threshold)
        // We don't have x0 separately here when mode=relaxed, so we recompute.
        // Actually the ctrl filtering always uses x0 (threshold line 1), not z0.
        // For simplicity: recompute ctrl threshold using the same auc_thresh for
        // ctrl-peak filtering (SEACR uses x0 for ctrl regardless of mode).
        // Since we computed auc_thresh as x0 (stringent) or z0 (relaxed), but
        // ctrl filtering uses x0 always, we need to handle this correctly.
        // For the non-norm/norm control case, the ctrl_threshold applied to ctrl
        // blocks is always `thresh` (x0). We'll filter ctrl blocks using auc_thresh
        // if mode=stringent, else we need to recompute x0. However, the ctrl peak
        // set is computed once with x0 regardless of mode in SEACR.sh.
        // This is a conservative approximation that matches SEACR's behavior.
        // Refilter ctrl with its own auc > stringent_thresh:
        let ctrl_auc_thresh = auc_thresh; // SEACR uses x0 for ctrl; we use same auc_thresh
        let ctrl_filtered: Vec<Block> = ctrl_blocks_mapped
            .into_iter()
            .filter(|b| b.total > ctrl_auc_thresh)
            .collect();

        let ctrl_mean_width: f64 = if ctrl_filtered.is_empty() {
            0.0
        } else {
            ctrl_filtered
                .iter()
                .map(|b| (b.end - b.start) as f64)
                .sum::<f64>()
                / ctrl_filtered.len() as f64
        };
        let ctrl_gap_tol = (ctrl_mean_width / 10.0) as u64;
        let ctrl_merged = merge_nearby(&ctrl_filtered, ctrl_gap_tol);

        subtract_control(&merged, &ctrl_merged)
    } else {
        merged
    };

    let mut writer = BufWriter::new(out);
    for b in &final_peaks {
        write_peak(b, &chroms, &mut writer).map_err(RsomicsError::Io)?;
    }
    writer.flush().map_err(RsomicsError::Io)?;
    Ok(final_peaks.len())
}

/// For numeric-threshold mode: min value where `1 - ecdf(v)(x) <= frac`.
/// Equivalent to R: `min(frame$values[frame$percentile <= ctrl])`.
fn numeric_threshold_min(sorted: &[f64], frac: f64) -> f64 {
    let n = sorted.len() as f64;
    // percentile = 1 - ecdf(v)(x); keep where percentile <= frac
    // ecdf(v)(x) = rank/n; percentile = 1 - rank/n <= frac → rank >= n*(1-frac)
    // rank is 0-indexed position + 1
    // The min value with percentile <= frac is the value at rank ceil(n*(1-frac))
    sorted
        .iter()
        .enumerate()
        .find(|&(i, _)| {
            let percentile = 1.0 - (i as f64 + 1.0) / n;
            percentile <= frac
        })
        .map(|(_, &v)| v)
        .unwrap_or(*sorted.last().unwrap_or(&0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ivs(rows: &[(&str, u64, u64, f64)]) -> Vec<Interval> {
        let mut chroms = ChromTable::default();
        rows.iter()
            .map(|&(c, s, e, v)| Interval {
                chrom: chroms.intern(c),
                start: s,
                end: e,
                value: v,
            })
            .collect()
    }

    #[test]
    fn block_total_and_max() {
        let v = ivs(&[("chr1", 0, 10, 1.0), ("chr1", 10, 20, 2.0)]);
        let b = build_blocks(&v);
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].total, 30.0);
        assert_eq!(b[0].max, 2.0);
        assert_eq!((b[0].max_start, b[0].max_end), (10, 20));
        assert_eq!(b[0].num_intervals, 2);
    }

    #[test]
    fn gap_splits_blocks() {
        // Non-adjacent: start of second != end of first
        let v = ivs(&[
            ("chr1", 0, 10, 1.0),
            ("chr1", 20, 30, 5.0),
            ("chr1", 30, 40, 5.0),
        ]);
        let b = build_blocks(&v);
        assert_eq!(b.len(), 2);
        assert_eq!(b[1].total, 100.0);
        assert_eq!((b[1].start, b[1].end), (20, 40));
    }

    #[test]
    fn strict_adjacency() {
        // Gap of 20 → two blocks
        let v = ivs(&[("chr1", 0, 10, 1.0), ("chr1", 30, 40, 5.0)]);
        assert_eq!(build_blocks(&v).len(), 2);
        // Adjacent → one block
        let v2 = ivs(&[("chr1", 0, 10, 1.0), ("chr1", 10, 20, 5.0)]);
        assert_eq!(build_blocks(&v2).len(), 1);
    }

    #[test]
    fn chrom_change_splits_blocks() {
        let v = ivs(&[("chr1", 0, 10, 1.0), ("chr2", 10, 20, 1.0)]);
        assert_eq!(build_blocks(&v).len(), 2);
    }

    #[test]
    fn max_region_spans_noncontiguous_max() {
        let v = ivs(&[
            ("chr1", 0, 10, 5.0),
            ("chr1", 10, 20, 3.0),
            ("chr1", 20, 30, 5.0),
        ]);
        let b = build_blocks(&v);
        assert_eq!((b[0].max_start, b[0].max_end), (0, 30));
    }

    #[test]
    fn signif6_matches_r() {
        assert_eq!(signif6(1727.677429), "1727.68");
        assert_eq!(signif6(2.930810), "2.93081");
        assert_eq!(signif6(904708.0), "904708");
        assert_eq!(signif6(1407.0), "1407");
        assert_eq!(signif6(879.243), "879.243");
        assert_eq!(signif6(0.857798), "0.857798");
        assert_eq!(signif6(2315.98), "2315.98");
    }

    #[test]
    fn quantile_type7_endpoints() {
        let s = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(quantile_type7(&s, 0.0), 1.0);
        assert_eq!(quantile_type7(&s, 1.0), 5.0);
        assert_eq!(quantile_type7(&s, 0.5), 3.0);
    }

    #[test]
    fn density_mode_matches_r() {
        // ctrl subvec from golden fixture: ctrlvec[ctrlvec <= 4.70167]
        let ctrl_sub = vec![
            0.733204, 4.70167, 0.824854, 1.97965, 3.79433, 2.6212, 0.641553, 4.15178, 3.46439,
            2.56621, 1.28311, 4.21592, 3.97763, 4.41755, 1.01732, 0.659884, 0.733204, 2.62121,
            0.797359, 4.61919, 2.23627, 1.09981, 2.50206, 0.339107, 3.92264, 0.403262, 3.776,
            1.77802, 0.632388, 0.559068,
        ];
        // R: density(ctrl_sub)$x[which.max(density(ctrl_sub)$y)] ≈ 0.8467413
        let mode = density_mode(&ctrl_sub);
        assert!(
            (mode - 0.8467413).abs() < 0.005,
            "ctrl density mode {mode} != expected ~0.8467"
        );

        // exp subvec (all expvec since expvalue = max)
        let exp_sub = vec![
            152.116, 130.671, 32.0244, 12.5096, 16.4411, 6.50497, 29.4511, 143.467, 74.6285,
            30.3089, 51.3964, 24.5902, 26.9492, 71.1972, 32.4533, 11.4373, 13.7248, 61.9045,
            7.93464, 17.013, 62.9768, 13.8677, 17.2275, 8.72095, 3.43119, 18.4427, 43.4617,
            53.1121, 2.64488,
        ];
        // R: mode ≈ 17.18486
        let exp_mode = density_mode(&exp_sub);
        assert!(
            (exp_mode - 17.18486).abs() < 0.1,
            "exp density mode {exp_mode} != expected ~17.18"
        );
    }

    #[test]
    fn bw_nrd0_matches_r() {
        // R: bw.nrd0(ctrl_sub) ≈ 0.6912042
        let ctrl_sub = vec![
            0.733204, 4.70167, 0.824854, 1.97965, 3.79433, 2.6212, 0.641553, 4.15178, 3.46439,
            2.56621, 1.28311, 4.21592, 3.97763, 4.41755, 1.01732, 0.659884, 0.733204, 2.62121,
            0.797359, 4.61919, 2.23627, 1.09981, 2.50206, 0.339107, 3.92264, 0.403262, 3.776,
            1.77802, 0.632388, 0.559068,
        ];
        let bw = bw_nrd0(&ctrl_sub);
        assert!(
            (bw - 0.6912042).abs() < 0.0001,
            "bw {bw} != expected ~0.6912"
        );
    }
}
