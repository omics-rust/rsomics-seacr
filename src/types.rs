use std::path::PathBuf;

/// One bedGraph interval: a constant-value span on a chromosome.
pub(crate) struct Interval {
    pub(crate) chrom: u32,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) value: f64,
}

/// Maps chromosome names to compact ids in first-seen order.
#[derive(Default)]
pub(crate) struct ChromTable {
    pub(crate) names: Vec<String>,
}

impl ChromTable {
    pub(crate) fn intern(&mut self, name: &str) -> u32 {
        if let Some(pos) = self.names.iter().position(|n| n == name) {
            return pos as u32;
        }
        self.names.push(name.to_owned());
        (self.names.len() - 1) as u32
    }

    pub(crate) fn name(&self, id: u32) -> &str {
        &self.names[id as usize]
    }
}

/// A signal block: a merged run of strictly adjacent nonzero intervals.
///
/// `num_intervals` counts input intervals merged into the block; it feeds
/// SEACR's R threshold step (expmax/ctrlmax).
pub struct Block {
    pub(crate) chrom: u32,
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
    /// `x0` AUC threshold (peak of the pctremain curve).
    Stringent,
    /// `z0` AUC threshold (midpoint of the pctremain curve).
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
    Control(PathBuf),
}
