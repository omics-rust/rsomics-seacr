use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_seacr::{Mode, Norm, Threshold, call_peaks};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-seacr",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Experimental (target) bedGraph, coordinate-sorted, nonzero signal only.
    pub experimental: PathBuf,

    /// IgG control bedGraph for an empirical threshold, or a numeric
    /// top-fraction in (0,1) (e.g. 0.01 = top 1% of blocks by total signal).
    #[arg(long = "control")]
    pub control: Option<PathBuf>,

    /// Numeric top-fraction threshold in (0,1). Mutually exclusive with --control.
    #[arg(long = "fraction", conflicts_with = "control")]
    pub fraction: Option<f64>,

    /// Normalise the control track to the experimental track: norm | non.
    #[arg(long = "norm", default_value = "norm")]
    pub norm: Norm,

    /// Threshold axis: stringent (x0 — peak of pctremain curve) | relaxed (z0 — midpoint).
    #[arg(long = "mode", default_value = "stringent")]
    pub mode: Mode,

    /// Output BED file (use `-` for stdout).
    #[arg(short = 'o', long, default_value = "-")]
    pub output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let threshold = match (self.control, self.fraction) {
            (Some(ctrl), None) => Threshold::Control(ctrl),
            (None, Some(frac)) => Threshold::Fraction(frac),
            (None, None) => {
                return Err(RsomicsError::InvalidInput(
                    "supply either --control <igg.bedgraph> or --fraction <n>".into(),
                ));
            }
            (Some(_), Some(_)) => unreachable!("clap conflicts_with rejects this"),
        };

        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };

        let n = call_peaks(
            &self.experimental,
            &threshold,
            self.norm,
            self.mode,
            &mut out,
        )?;

        if !self.common.quiet {
            eprintln!("{n} peaks called");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "CUT&RUN peak caller: bedGraph signal → BED peaks (clean-room SEACR port).",
    origin: Some(Origin {
        upstream: "SEACR (algorithm + constants read; independent Rust reimplementation)",
        upstream_license: "GPL-3.0 (upstream); ours MIT OR Apache-2.0",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1186/s13072-019-0287-4"),
    }),
    usage_lines: &[
        "<exp.bedgraph> --fraction 0.01 --mode stringent -o peaks.bed",
        "<exp.bedgraph> --control igg.bedgraph --norm norm --mode relaxed -o peaks.bed",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "control",
                aliases: &[],
                value: Some("<igg.bedgraph>"),
                type_hint: Some("path"),
                required: false,
                default: None,
                description: "IgG control bedGraph for an empirical threshold.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "fraction",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("f64"),
                required: false,
                default: None,
                description: "Numeric top-fraction in (0,1) by total signal.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "norm",
                aliases: &[],
                value: Some("<norm|non>"),
                type_hint: Some("str"),
                required: false,
                default: Some("norm"),
                description: "Normalise control to experimental signal.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "mode",
                aliases: &[],
                value: Some("<stringent|relaxed>"),
                type_hint: Some("str"),
                required: false,
                default: Some("stringent"),
                description: "stringent uses x0 (peak of pctremain curve); relaxed uses z0 (midpoint).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Top 1% of blocks by total signal, stringent",
            command: "rsomics-seacr exp.bedgraph --fraction 0.01 --mode stringent -o peaks.bed",
        },
        Example {
            description: "Empirical threshold from a normalised IgG control, relaxed",
            command: "rsomics-seacr exp.bedgraph --control igg.bedgraph --norm norm --mode relaxed -o peaks.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }

    #[test]
    fn control_and_fraction_conflict() {
        let r = Cli::try_parse_from([
            "rsomics-seacr",
            "x.bg",
            "--control",
            "c.bg",
            "--fraction",
            "0.01",
        ]);
        assert!(r.is_err());
    }
}
