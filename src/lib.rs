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

mod blocks;
mod io;
mod model;
mod peaks;
mod types;

pub use io::signif6;
pub use model::{density_mode, quantile_type7};
pub use peaks::call_peaks;
pub use types::{Block, Mode, Norm, Threshold};
