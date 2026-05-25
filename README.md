# rsomics-seacr

CUT&RUN peak caller: bedGraph signal → BED peaks. Independent Rust reimplementation
of the SEACR algorithm (Sparse Enrichment Analysis for CUT&RUN).

## Usage

```
rsomics-seacr <exp.bedgraph> --fraction 0.01 --mode stringent -o peaks.bed
rsomics-seacr <exp.bedgraph> --control igg.bedgraph --norm norm --mode relaxed -o peaks.bed
```

The experimental bedGraph must be coordinate-sorted and contain only nonzero
signal intervals (as produced by `bedtools genomecov -bg` on fragment BEDs, or
by `rsomics-bam-signal`).

## Method

A *signal block* is a maximal run of strictly adjacent nonzero bedGraph intervals
on one chromosome. Per block:

- **total signal** (AUC) = Σ value·(end − start) — the area under the signal curve
- **max signal** = the largest value attained at any base
- **max region** = the span from the farthest-upstream to the farthest-downstream
  base that attains the max signal
- **num_intervals** = count of input intervals merged into the block

A global threshold separates peaks from background:

- **Numeric mode** (`--fraction n`, `n ∈ (0,1)`) keeps the top-`n` blocks by AUC.
  `stringent` uses the AUC quantile at `1 − n`; `relaxed` uses the num-intervals
  quantile at `1 − n`.
- **Control mode** (`--control igg.bedgraph`) computes empirical AUC thresholds
  (`x0` = stringent, `z0` = relaxed) by maximising the pctremain curve from the
  SEACR algorithm. `--norm norm` first scales the control signal to the experimental
  signal by the ratio of their Gaussian KDE density modes; `--norm non` skips that.

After thresholding, nearby peaks within `mean_block_width / 10` bases are merged,
and peaks overlapping any control-enriched region are removed.

Output columns: `chrom  start  end  total_signal  max_signal  chrom:maxstart-maxend`.
Numeric columns are rendered with six significant figures (`signif(x, 6)`),
matching SEACR's R output formatting exactly.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--fraction` | — | Numeric top-fraction in (0,1) by AUC. Mutually exclusive with `--control`. |
| `--control` | — | IgG control bedGraph for an empirical threshold. |
| `--norm` | norm | `norm` scales control to experimental signal; `non` skips. |
| `--mode` | stringent | `stringent` = x0 (peak of pctremain curve); `relaxed` = z0 (midpoint). |
| `-o` | stdout | Output BED path. |

## Compatibility with SEACR

The output of `rsomics-seacr` is **byte-identical** to `SEACR_1.3.sh` for all six
mode combinations (`--control × {non,norm} × {stringent,relaxed}` and
`--fraction × {stringent,relaxed}`). This is verified by `tests/compat.rs` against
the real SEACR binary on the same golden fixtures.

## Origin

This crate is an independent Rust reimplementation of the SEACR method informed by:

- The published method: Meers MP, Tenenbaum D, Henikoff S. *Peak calling by
  Sparse Enrichment Analysis for CUT&RUN chromatin profiling.* Epigenetics &
  Chromatin 12(1):42, 2019. DOI: 10.1186/s13072-019-0287-4
- The SEACR algorithm and its constants as published in `SEACR_1.3.R` (GPL-3.0),
  read to extract algorithm details (pctremain curve, knee-point detection,
  density-mode normalisation, spurious-threshold correction). The implementation
  is original Rust, not a transcription of the R or shell source.
- Black-box behaviour testing against the upstream `SEACR_1.3.sh` binary.

Test fixtures are slices of the public SEACR test dataset (CTCF / IgG, chr1).

License: MIT OR Apache-2.0.
Upstream credit: SEACR <https://github.com/FredHutch/SEACR> (GPL-3.0).
