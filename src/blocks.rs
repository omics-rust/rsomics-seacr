use crate::types::{Block, Interval};

/// Collapse nonzero intervals into signal blocks using strict adjacency
/// (`interval.start == previous.end`), matching SEACR's AWK `$2==stop` check.
pub(crate) fn build_blocks(intervals: &[Interval]) -> Vec<Block> {
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

/// Build blocks and apply SEACR's pipeline quirks:
/// 1. Skip the first data line (AWK `BEGIN{s=1}` skips line 1).
/// 2. Drop the final block (AWK lacks `END{}`; the last accumulated block is never printed).
/// 3. Round each block's total to 6 significant figures (AWK's default `OFMT = %.6g`).
pub(crate) fn seacr_blocks(intervals: &[Interval]) -> Vec<Block> {
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
/// AWK reads them back it parses the rounded string values; all downstream
/// arithmetic must use the same quantities.
fn round_block_totals(mut blocks: Vec<Block>) -> Vec<Block> {
    for b in &mut blocks {
        b.total = round6g(b.total);
    }
    blocks
}

/// Round to 6 significant figures (awk `%.6g` semantics, round-half-away-from-zero).
pub(crate) fn round6g(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let mag = x.abs().log10().floor() as i32;
    let factor = 10f64.powi(5 - mag);
    (x * factor).round() / factor
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ChromTable;

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
        let v = ivs(&[("chr1", 0, 10, 1.0), ("chr1", 30, 40, 5.0)]);
        assert_eq!(build_blocks(&v).len(), 2);
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
}
