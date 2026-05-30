use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::blocks::seacr_blocks;
use crate::io::{parse_bedgraph, write_peak};
use crate::model::{compute_thresholds, density_mode, knee_value, numeric_threshold_min};
use crate::types::{Block, ChromTable, Mode, Norm, Threshold};

/// Call peaks from `experimental` and write the BED to `out`. Returns peak count.
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
            let mut sorted_auc = expvec.clone();
            sorted_auc.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut sorted_num = expmax.clone();
            sorted_num.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let (x0, z0) = match mode {
                Mode::Stringent => {
                    let x0 = numeric_threshold_min(&sorted_auc, *frac);
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
            (auc, 0.0)
        }
        Threshold::Control(ctrl_path) => {
            let mut ctrl_chroms = ChromTable::default();
            let ctrl_intervals = parse_bedgraph(ctrl_path, &mut ctrl_chroms)?;
            let ctrl_blocks = seacr_blocks(&ctrl_intervals);

            let mut ctrlvec: Vec<f64> = ctrl_blocks.iter().map(|b| b.total).collect();
            let ctrlmax: Vec<f64> = ctrl_blocks.iter().map(|b| b.num_intervals as f64).collect();

            if norm == Norm::On {
                // Knee-point detection, then density-mode ratio normalisation
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

    let filtered: Vec<&Block> = exp_blocks
        .iter()
        .filter(|b| b.total > auc_thresh && b.num_intervals as f64 > num_thresh)
        .collect();

    if filtered.is_empty() {
        return Ok(0);
    }

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

    let final_peaks = if let Threshold::Control(ctrl_path) = threshold {
        // Reparse control using original (pre-normalisation) signal for subtraction.
        let mut ctrl_chroms2 = ChromTable::default();
        let ctrl_intervals2 = parse_bedgraph(ctrl_path, &mut ctrl_chroms2)?;
        let ctrl_blocks2 = seacr_blocks(&ctrl_intervals2);
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

        // SEACR.sh line 159: ctrl peaks filtered with x0 (always stringent threshold).
        let ctrl_auc_thresh = auc_thresh;
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

/// Merge adjacent blocks (gaps < `gap_tolerance`).
///
/// The final merged block is dropped to match SEACR's merge AWK, which lacks
/// an `END{}` clause and never emits its last accumulated result.
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
    // SEACR's merge AWK also lacks END{} — final accumulated block is dropped.
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
