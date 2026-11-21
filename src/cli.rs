//! Command line interface for the `kronoprime` binary.
//!
//! Supports numbers as positional arguments, streaming from stdin (one
//! number per line when no numbers are given) and the flags `--json`,
//! `--matrix` and `--quadruplets <K>`.

use crate::factor::factor;
use crate::matrix::{self, Placement};
use crate::primality::is_prime;
use std::io::{self, BufRead, IsTerminal, Write};
use std::time::Instant;

/// Upper bound accepted by `--quadruplets` to keep runs bounded.
const MAX_QUADRUPLETS: usize = 100_000;

/// Result of analysing a single input value.
#[derive(Debug, Clone)]
pub struct Analysis {
    /// Value that was analysed.
    pub n: u128,
    /// True when `n` is prime.
    pub prime: bool,
    /// Placement in the 3x3 decade grid (`None` for 0 and multiples of 10).
    pub placement: Option<Placement>,
    /// Prime factors in ascending order (empty for `n < 2`).
    pub factors: Vec<u128>,
    /// Wall-clock analysis time in milliseconds.
    pub elapsed_ms: f64,
}

/// Run the CLI with the given arguments (program arguments, excluding the
/// binary name) and return the process exit code.
pub fn run(args: &[String]) -> u8 {
    let mut json = false;
    let mut show_matrix = false;
    let mut quadruplets: Option<usize> = None;
    let mut inputs: Vec<u128> = Vec::new();
    let mut invalid = 0usize;

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--quadruplets=") {
            match rest.parse::<usize>() {
                Ok(k) => quadruplets = Some(k),
                Err(_) => return reject_quadruplets(),
            }
            continue;
        }
        match arg.as_str() {
            "--json" => json = true,
            "--matrix" => show_matrix = true,
            "--help" | "-h" => {
                print_usage();
                return 0;
            }
            "--version" => {
                println!("kronoprime {}", env!("CARGO_PKG_VERSION"));
                return 0;
            }
            "--quadruplets" => match iter.next().and_then(|v| v.parse::<usize>().ok()) {
                Some(k) => quadruplets = Some(k),
                None => return reject_quadruplets(),
            },
            other => match other.parse::<u128>() {
                Ok(v) => inputs.push(v),
                Err(_) => {
                    eprintln!(
                        "error: invalid input '{other}' (expected an integer in 0..=u128::MAX)"
                    );
                    invalid += 1;
                }
            },
        }
    }

    if let Some(k) = quadruplets {
        return quadruplet_report(k, json);
    }

    if inputs.is_empty() {
        if io::stdin().is_terminal() {
            print_usage();
            return 2;
        }
        for line in io::stdin().lock().lines() {
            match line {
                Ok(line) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    match trimmed.parse::<u128>() {
                        Ok(v) => inputs.push(v),
                        Err(_) => {
                            eprintln!("error: invalid input '{trimmed}' (expected an integer in 0..=u128::MAX)");
                            invalid += 1;
                        }
                    }
                }
                Err(_) => {
                    eprintln!("error: failed to read stdin");
                    return 1;
                }
            }
        }
    }

    if inputs.is_empty() {
        return u8::from(invalid > 0);
    }

    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (i, &n) in inputs.iter().enumerate() {
        let analysis = analyze(n);
        if json {
            let _ = writeln!(out, "{}", json_report(&analysis));
        } else {
            if i > 0 {
                let _ = writeln!(out);
            }
            let _ = write!(out, "{}", text_report(&analysis, show_matrix));
        }
        let _ = out.flush();
    }
    u8::from(invalid > 0)
}

fn reject_quadruplets() -> u8 {
    eprintln!("error: --quadruplets expects a positive integer");
    2
}

/// Analyse one value: primality, 3x3 placement and prime factorization.
pub fn analyze(n: u128) -> Analysis {
    let start = Instant::now();
    let prime = is_prime(n);
    let factors = if n >= 2 { factor(n) } else { Vec::new() };
    let elapsed_ms = start.elapsed().as_secs_f64() * 1_000.0;
    Analysis {
        n,
        prime,
        placement: matrix::try_placement(n),
        factors,
        elapsed_ms,
    }
}

/// Human-readable status label for an analysis.
pub fn status_of(analysis: &Analysis) -> &'static str {
    if analysis.n < 2 {
        "NOT PRIME (0 and 1 have no prime factorization)"
    } else if analysis.prime {
        "PRIME"
    } else if analysis.factors.first().is_some_and(|&f| f > 1_000_000) {
        "COMPOSITE (large prime factors)"
    } else {
        "COMPOSITE"
    }
}

/// Plain-text report for one analysis; appends the 3x3 grid when requested.
pub fn text_report(analysis: &Analysis, show_matrix: bool) -> String {
    let mut out = String::new();
    out.push_str(&format!("Input: {}\n", analysis.n));
    out.push_str(&format!("Status: {}\n", status_of(analysis)));
    match analysis.placement {
        Some(p) => out.push_str(&format!("Matrix N: {} | Slot: +{}\n", p.index, p.slot)),
        None => out.push_str("Matrix N: n/a (value on a decade border) | Slot: n/a\n"),
    }
    let factors = analysis
        .factors
        .iter()
        .map(u128::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    out.push_str(&format!("Factors: [{factors}]\n"));
    out.push_str(&format!("Elapsed: {:.2} ms\n", analysis.elapsed_ms));
    if show_matrix {
        out.push_str(&matrix_block(analysis));
    }
    out
}

/// Structured JSON report for one analysis. Large integers are emitted as
/// strings so no precision is lost by downstream consumers.
pub fn json_report(analysis: &Analysis) -> String {
    let factors = analysis
        .factors
        .iter()
        .map(|f| format!("\"{f}\""))
        .collect::<Vec<_>>()
        .join(",");
    let (matrix_n, slot) = match analysis.placement {
        Some(p) => (format!("\"{}\"", p.index), p.slot.to_string()),
        None => ("null".to_string(), "null".to_string()),
    };
    format!(
        "{{\"input\":\"{}\",\"status\":\"{}\",\"prime\":{},\"matrix_n\":{},\"slot\":{},\"factors\":[{factors}],\"elapsed_ms\":{:.3}}}",
        analysis.n,
        status_of(analysis),
        analysis.prime,
        matrix_n,
        slot,
        analysis.elapsed_ms,
    )
}

/// Render the 3x3 decade grid of the analysed value (skipped on decade
/// borders, which have no grid cell).
fn matrix_block(analysis: &Analysis) -> String {
    let Some(p) = analysis.placement else {
        return String::new();
    };
    let base = 10 * (p.index - 1);
    let width = (base + 9).to_string().len();
    let mut out = String::new();
    out.push_str(&format!(
        "3x3 grid for Matrix N={} ({}..{}):\n",
        p.index,
        base + 1,
        base + 9
    ));
    for row in matrix::decade_cells(p.index) {
        for &value in &row {
            let mark = if is_prime(value) { "*" } else { " " };
            out.push_str(&format!("{value:>width$}{mark} "));
        }
        out.push('\n');
    }
    out.push_str("Active slots:");
    for s in matrix::ACTIVE_SLOTS {
        let value = base + u128::from(s);
        let kind = if is_prime(value) {
            "prime"
        } else {
            "composite"
        };
        let tag = if Some(s) == analysis.placement.map(|q| q.slot) {
            " (input)"
        } else {
            ""
        };
        out.push_str(&format!(" +{s}={value} {kind}{tag}"));
    }
    out.push('\n');
    out
}

/// Report the first `k` quadruplet matrix indices.
fn quadruplet_report(k: usize, json: bool) -> u8 {
    if k == 0 || k > MAX_QUADRUPLETS {
        eprintln!("error: --quadruplets expects an integer in 1..={MAX_QUADRUPLETS}");
        return 2;
    }
    let start = Instant::now();
    let list = matrix::find_quadruplets(k);
    let elapsed_ms = start.elapsed().as_secs_f64() * 1_000.0;
    if json {
        let items = list
            .iter()
            .map(|&n| {
                let base = 10 * (n - 1);
                let values = matrix::ACTIVE_SLOTS
                    .iter()
                    .map(|&s| format!("\"{}\"", base + u128::from(s)))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("{{\"matrix_n\":\"{n}\",\"values\":[{values}]}}")
            })
            .collect::<Vec<_>>()
            .join(",");
        println!("[{items}]");
    } else {
        for &n in &list {
            let base = 10 * (n - 1);
            let values = matrix::ACTIVE_SLOTS
                .iter()
                .map(|&s| (base + u128::from(s)).to_string())
                .collect::<Vec<_>>()
                .join(", ");
            println!("Quadruplet: Matrix N={n} -> [{values}]");
        }
        println!(
            "Searched with the mod-{} wheel in {:.2} ms",
            matrix::WHEEL_MODULUS,
            elapsed_ms
        );
    }
    0
}

fn print_usage() {
    println!(
        "kronoprime {} - analyse and factorize integers up to u128::MAX\n\n\
         Usage: kronoprime [options] [number ...]\n\n\
         With no numbers, reads one integer per line from stdin.\n\n\
         Options:\n\
         \x20 --json               structured JSON output\n\
         \x20 --matrix             also print the 3x3 decade grid of each input\n\
         \x20 --quadruplets <K>    find the first K decade indices whose +1/+3/+7/+9\n\
         \x20                        slots are all prime (mod-3003 wheel)\n\
         \x20 -h, --help           show this help\n\
         \x20 --version            show version",
        env!("CARGO_PKG_VERSION")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn status_labels() {
        assert_eq!(status_of(&analyze(7)), "PRIME");
        assert_eq!(status_of(&analyze(12)), "COMPOSITE");
        assert_eq!(
            status_of(&analyze(100000000520000000627)),
            "COMPOSITE (large prime factors)"
        );
        assert!(status_of(&analyze(1)).starts_with("NOT PRIME"));
    }

    #[test]
    fn text_report_shape() {
        let report = text_report(&analyze(100000000520000000627), true);
        assert!(report.contains("Input: 100000000520000000627"));
        assert!(report.contains("Matrix N: 10000000052000000063 | Slot: +7"));
        assert!(report.contains("Factors: [10000000019, 10000000033]"));
        assert!(report.contains("Elapsed:"));
        assert!(report.contains("+7=100000000520000000627 composite (input)"));
    }

    #[test]
    fn json_report_shape() {
        let json = json_report(&analyze(101));
        assert!(json.contains("\"input\":\"101\""));
        assert!(json.contains("\"prime\":true"));
        assert!(json.contains("\"matrix_n\":\"11\""));
        assert!(json.contains("\"slot\":1"));
        assert!(json.contains("\"factors\":[\"101\"]"));
        // Decade borders have no placement.
        let border = json_report(&analyze(100));
        assert!(border.contains("\"matrix_n\":null"));
        assert!(border.contains("\"slot\":null"));
        assert!(border.contains("\"factors\":[\"2\",\"2\",\"5\",\"5\"]"));
    }

    #[test]
    fn run_reports_invalid_inputs_but_keeps_valid_ones() {
        let args: Vec<String> = ["12", "not-a-number"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        // Exit code 1 flags the invalid input; the valid one is still printed.
        assert_eq!(run(&args), 1);
    }
}
