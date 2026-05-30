use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::types::{Block, ChromTable, Interval};

/// Parse a bedGraph into intervals, skipping the first data line.
///
/// SEACR's AWK pipeline (`BEGIN{s=1}; {if(s==1){s++}...`) skips the first
/// non-header line before recording anything. We replicate that.
pub(crate) fn parse_bedgraph(path: &Path, chroms: &mut ChromTable) -> Result<Vec<Interval>> {
    let file = std::fs::File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("reading {}: {e}", path.display())))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut out = Vec::new();
    let mut line = String::new();
    let mut lineno = 0usize;
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
            continue; // skip first data line — SEACR AWK quirk
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

pub(crate) fn write_peak(
    b: &Block,
    chroms: &ChromTable,
    w: &mut impl Write,
) -> std::io::Result<()> {
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
