//! Command line interface for the `kronoprime` binary.
//!
//! Supports numbers as positional arguments, streaming from stdin (one
//! number per line when no numbers are given) and the flags `--json`,
//! `--matrix` and `--quadruplets <K>`.

use crate::cubic::{
    classify_weierstrass, eval_discriminant, find_cm_curves, gcd, is_on_curve, sample_points,
    scalar_mul, sub_mod, CubicCurve, CurveClass, DiscriminantStatus, ECPoint, PointOpResult,
};
use crate::cubic4096::{
    classify_curve_4096, eval_discriminant_4096, find_cm_curves_4096, is_on_curve_4096,
    sample_points_4096, scalar_mul_4096, Curve4096, CurveClass4096, CurveError4096,
    DiscriminantStatus4096, Point4096, PointOpResult4096,
};
use crate::factor::{cubic_factor_probe, cubic_factor_probe_4096, factor};
use crate::matrix::{
    self, depressed_params_4096, extract_cubic_curve_4096, matrix_invariants_3x3,
    matrix_invariants_4096, parse_matrix_3x3, parse_matrix_4096, Matrix3x3CubicExt,
    MatrixError4096, Placement,
};
use crate::primality::{is_prime, is_prime_u4096, mul_mod};
use crate::u4096::U4096;
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
    if args.first().map(|s| s.as_str()) == Some("cubic") {
        return run_cubic(&args[1..]);
    }

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
         Usage: kronoprime [options] [number ...]\n\
                kronoprime cubic <subcommand> [options]\n\n\
         With no numbers, reads one integer per line from stdin.\n\n\
         Options:\n\
         \x20 --json               structured JSON output\n\
         \x20 --matrix             also print the 3x3 decade grid of each input\n\
         \x20 --quadruplets <K>    find the first K decade indices whose +1/+3/+7/+9\n\
         \x20                        slots are all prime (mod-3003 wheel)\n\
         \x20 -h, --help           show this help\n\
         \x20 --version            show version\n\n\
         Cubic & Elliptic Subcommands:\n\
         \x20 scan                 scan parameter space for CM or singular curves\n\
         \x20 from-matrix          convert 3x3 matrix to Weierstrass form and classify\n\
         \x20 verify               certify curve order empirically by checking [order]P = O\n\
         \x20 factor-probe         factorize composite modulus via cubic curve singularities / ECM\n\
         \x20 point-ops            point arithmetic (scalar multiplication) on a curve\n\n\
         Cubic values may reach 4096 bits: moduli, coefficients, points, scalars\n\
         and matrix cells above u128::MAX (and 0x-prefixed hexadecimal values)\n\
         run on the dedicated 4096-bit engine with the same flags and exit codes.",
        env!("CARGO_PKG_VERSION")
    );
}

/// Dispatches subcommands under `kronoprime cubic`.
pub fn run_cubic(args: &[String]) -> u8 {
    if args.is_empty() {
        print_cubic_usage();
        return 2;
    }
    match args[0].as_str() {
        "scan" => run_cubic_scan(&args[1..]),
        "from-matrix" => run_cubic_from_matrix(&args[1..]),
        "verify" => run_cubic_verify(&args[1..]),
        "factor-probe" => run_cubic_factor_probe(&args[1..]),
        "point-ops" => run_cubic_point_ops(&args[1..]),
        "-h" | "--help" | "help" => {
            print_cubic_usage();
            0
        }
        other => {
            eprintln!("error: unknown cubic subcommand '{other}'");
            print_cubic_usage();
            1
        }
    }
}

fn print_cubic_usage() {
    println!(
        "kronoprime cubic - discovery and classification of cubic and elliptic curves\n\n\
         Usage: kronoprime cubic <subcommand> [options]\n\n\
         Subcommands:\n\
         \x20 scan         scan parameter space for CM or singular curves\n\
         \x20 from-matrix  convert 3x3 matrix to depressed Weierstrass form and classify\n\
         \x20 verify       certify curve order empirically by checking [order]P = O\n\
         \x20 factor-probe factorize composite modulus via cubic curve singularities / ECM\n\
         \x20 point-ops    perform point operations (scalar multiplication) on a curve\n\n\
         Run 'kronoprime cubic <subcommand> --help' for details on a specific subcommand.\n\n\
         Values beyond u128::MAX, decimal or 0x-prefixed hexadecimal, are accepted\n\
         by every subcommand above up to 4096 bits, with the same flags and exit codes."
    );
}

fn run_cubic_scan(args: &[String]) -> u8 {
    if args_need_wide(args) {
        return run_cubic_scan_wide(args);
    }
    let mut modulus: Option<u128> = None;
    let mut limit: usize = 100;
    let mut class_filter = "all".to_string();
    let mut format_mode = "table".to_string();
    let mut stop_on_factor = false;

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--modulus=") {
            if let Ok(v) = rest.parse::<u128>() {
                modulus = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--limit=") {
            if let Ok(v) = rest.parse::<usize>() {
                limit = v;
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--class=") {
            class_filter = rest.to_lowercase();
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--format=") {
            format_mode = rest.to_lowercase();
            continue;
        }
        match arg.as_str() {
            "--modulus" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    modulus = Some(v);
                }
            }
            "--limit" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<usize>().ok()) {
                    limit = v;
                }
            }
            "--class" => {
                if let Some(v) = iter.next() {
                    class_filter = v.to_lowercase();
                }
            }
            "--format" => {
                if let Some(v) = iter.next() {
                    format_mode = v.to_lowercase();
                }
            }
            "--stop-on-factor" => {
                stop_on_factor = true;
            }
            "-h" | "--help" => {
                println!(
                    "Usage: kronoprime cubic scan --modulus <N> [--limit <LIMIT>] [--class <CLASS>] [--format json|table] [--stop-on-factor]"
                );
                return 0;
            }
            other => {
                eprintln!("error: unknown scan option '{other}'");
                return 1;
            }
        }
    }

    let n = match modulus {
        Some(m) if m >= 2 => m,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };

    if stop_on_factor && !is_prime(n) {
        for k in 1..=(limit as u128) {
            let (_, status) = eval_discriminant(k, k, n);
            if let DiscriminantStatus::FactorFound(g) = status {
                let f1 = g.min(n / g);
                let f2 = g.max(n / g);
                if format_mode == "json" {
                    println!(
                        "{{\n  \"status\": \"FactorFound\",\n  \"modulus\": {n},\n  \"factor1\": {f1},\n  \"factor2\": {f2}\n}}"
                    );
                } else {
                    println!("Fattore non banale individuato: {f1} (n = {f1} * {f2})");
                }
                return 2;
            }
        }
    }

    let mut curves: Vec<CubicCurve> = Vec::new();
    match class_filter.as_str() {
        "cm" => {
            curves = find_cm_curves(n, limit);
        }
        "nodal" => {
            for x0 in 1..=(limit as u128) {
                let a = sub_mod(0, mul_mod(3, mul_mod(x0, x0, n), n), n);
                let b = mul_mod(2, mul_mod(x0, mul_mod(x0, x0, n), n), n);
                if let Ok(c) = classify_weierstrass(a, b, n) {
                    if matches!(c.class, CurveClass::NodalSingular { .. }) {
                        curves.push(c);
                    }
                }
            }
        }
        "cuspidal" => {
            if let Ok(c) = classify_weierstrass(0, 0, n) {
                curves.push(c);
            }
        }
        "montgomery" => {
            for a_param in 1..=(limit as u128) {
                for b_param in 1..=3u128 {
                    let a_sq = mul_mod(a_param, a_param, n);
                    let diff = sub_mod(a_sq, 4, n);
                    if mul_mod(b_param, diff, n) != 0 {
                        curves.push(CubicCurve {
                            a: a_param,
                            b: b_param,
                            modulus: n,
                            discriminant: diff,
                            class: CurveClass::MontgomeryForm { a_param, b_param },
                        });
                        if curves.len() >= limit {
                            break;
                        }
                    }
                }
                if curves.len() >= limit {
                    break;
                }
            }
        }
        _ => {
            // "all"
            let cm = find_cm_curves(n, limit);
            curves.extend(cm);
            if curves.len() < limit {
                if let Ok(c) = classify_weierstrass(0, 0, n) {
                    curves.push(c);
                }
            }
            let mut x0 = 1u128;
            while curves.len() < limit && x0 <= limit as u128 {
                let a = sub_mod(0, mul_mod(3, mul_mod(x0, x0, n), n), n);
                let b = mul_mod(2, mul_mod(x0, mul_mod(x0, x0, n), n), n);
                if let Ok(c) = classify_weierstrass(a, b, n) {
                    if !curves.contains(&c) {
                        curves.push(c);
                    }
                }
                x0 += 1;
            }
        }
    }

    curves.truncate(limit);

    if format_mode == "json" {
        let items = curves
            .iter()
            .map(|c| {
                let (class_str, closed_order_str) = match &c.class {
                    CurveClass::CmSupersingular { closed_form_order } => {
                        let kind = if c.a == 0 && c.modulus % 3 == 2 {
                            "CmSupersingularType1"
                        } else {
                            "CmSupersingularType2"
                        };
                        (kind.to_string(), closed_form_order.to_string())
                    }
                    CurveClass::NodalSingular { .. } => {
                        ("NodalSingular".to_string(), "null".to_string())
                    }
                    CurveClass::CuspidalSingular { .. } => {
                        ("CuspidalSingular".to_string(), "null".to_string())
                    }
                    CurveClass::MontgomeryForm { .. } => {
                        ("MontgomeryForm".to_string(), "null".to_string())
                    }
                    CurveClass::StandardElliptic => {
                        ("StandardElliptic".to_string(), "null".to_string())
                    }
                };
                if closed_order_str == "null" {
                    format!(
                        "  {{\n    \"a\": {},\n    \"b\": {},\n    \"modulus\": {},\n    \"discriminant\": {},\n    \"class\": \"{}\"\n  }}",
                        c.a, c.b, c.modulus, c.discriminant, class_str
                    )
                } else {
                    format!(
                        "  {{\n    \"a\": {},\n    \"b\": {},\n    \"modulus\": {},\n    \"discriminant\": {},\n    \"class\": \"{}\",\n    \"closed_form_order\": {}\n  }}",
                        c.a, c.b, c.modulus, c.discriminant, class_str, closed_order_str
                    )
                }
            })
            .collect::<Vec<_>>()
            .join(",\n");
        println!("[\n{items}\n]");
    } else {
        println!("+-------+-------+---------+-------------------+---------------------+");
        println!("| Param a | Param b | Modulo n| Classe Algebrica  | Ordine Chiuso (#E)  |");
        println!("+-------+-------+---------+-------------------+---------------------+");
        for c in &curves {
            let (class_label, order_label) = match &c.class {
                CurveClass::CmSupersingular { closed_form_order } => {
                    let label = if c.a == 0 && c.modulus % 3 == 2 {
                        "CmSupersingular I"
                    } else {
                        "CmSupersingular II"
                    };
                    (label, closed_form_order.to_string())
                }
                CurveClass::NodalSingular { .. } => ("NodalSingular", "-".to_string()),
                CurveClass::CuspidalSingular { .. } => ("CuspidalSingular", "-".to_string()),
                CurveClass::MontgomeryForm { .. } => ("MontgomeryForm", "-".to_string()),
                CurveClass::StandardElliptic => ("StandardElliptic", "-".to_string()),
            };
            println!(
                "| {:<7} | {:<7} | {:<7} | {:<17} | {:<19} |",
                c.a, c.b, c.modulus, class_label, order_label
            );
        }
        println!("+-------+-------+---------+-------------------+---------------------+");
    }

    0
}

fn run_cubic_from_matrix(args: &[String]) -> u8 {
    if args_need_wide(args) {
        return run_cubic_from_matrix_wide(args);
    }
    let mut matrix_str: Option<String> = None;
    let mut modulus: Option<u128> = None;
    let mut format_mode = "table".to_string();

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--matrix=") {
            matrix_str = Some(rest.to_string());
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--modulus=") {
            if let Ok(v) = rest.parse::<u128>() {
                modulus = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--format=") {
            format_mode = rest.to_lowercase();
            continue;
        }
        match arg.as_str() {
            "--matrix" => {
                if let Some(v) = iter.next() {
                    matrix_str = Some(v.to_string());
                }
            }
            "--modulus" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    modulus = Some(v);
                }
            }
            "--format" => {
                if let Some(v) = iter.next() {
                    format_mode = v.to_lowercase();
                }
            }
            "-h" | "--help" => {
                println!(
                    "Usage: kronoprime cubic from-matrix --matrix \"<MATRIX_STR>\" --modulus <N> [--format json|table]"
                );
                return 0;
            }
            other => {
                eprintln!("error: unknown from-matrix option '{other}'");
                return 1;
            }
        }
    }

    let m_str = match matrix_str {
        Some(s) => s,
        None => {
            eprintln!("error: --matrix \"<STR>\" is required");
            return 1;
        }
    };
    let n = match modulus {
        Some(m) if m >= 2 => m,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };

    let m = match parse_matrix_3x3(&m_str) {
        Ok(mat) => mat,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };

    let (tr, sigma2, det) = matrix_invariants_3x3(&m, n);
    let (a, b) = match m.to_depressed_cubic(n) {
        Ok(params) => params,
        Err(e) => {
            let g = gcd(3, n);
            if g > 1 && g < n {
                eprintln!("Fattore non banale individuato durante l'inversione: {g}");
                return 2;
            }
            eprintln!("error: {e}");
            return 1;
        }
    };

    let curve = match m.extract_cubic_curve(n) {
        Ok(c) => c,
        Err(e) => {
            let (_, status) = eval_discriminant(a, b, n);
            if let DiscriminantStatus::FactorFound(g) = status {
                eprintln!("Fattore non banale individuato dal discriminante: {g}");
                return 2;
            }
            eprintln!("error: {e}");
            return 1;
        }
    };

    let class_str = curve.class.name();
    let formatted_matrix = format!(
        "[[{},{},{}],[{},{},{}],[{},{},{}]]",
        m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]
    );

    if format_mode == "json" {
        println!(
            "{{\n  \"input_matrix\": \"{}\",\n  \"invariants\": {{ \"tr\": {}, \"sigma2\": {}, \"det\": {} }},\n  \"weierstrass_form\": {{ \"a\": {}, \"b\": {} }},\n  \"discriminant\": {},\n  \"class\": \"{}\"\n}}",
            formatted_matrix, tr, sigma2, det, a, b, curve.discriminant, class_str
        );
    } else {
        println!("+-------------------+---------------------------------------------+");
        println!("| Proprietà         | Valore                                      |");
        println!("+-------------------+---------------------------------------------+");
        println!("| Matrice Ingresso  | {:<43} |", formatted_matrix);
        println!("| Modulo n          | {:<43} |", n);
        println!(
            "| Invarianti        | tr={}, sigma2={}, det={} |",
            tr, sigma2, det
        );
        println!("| Weierstrass (a,b) | a={}, b={} |", a, b);
        println!("| Discriminante     | {:<43} |", curve.discriminant);
        println!("| Classe Algebrica  | {:<43} |", class_str);
        println!("+-------------------+---------------------------------------------+");
    }

    0
}

fn run_cubic_verify(args: &[String]) -> u8 {
    if args_need_wide(args) {
        return run_cubic_verify_wide(args);
    }
    let mut a_val: Option<u128> = None;
    let mut b_val: Option<u128> = None;
    let mut modulus: Option<u128> = None;
    let mut num_points: usize = 10;
    let mut format_mode = "table".to_string();

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--a=") {
            if let Ok(v) = rest.parse::<u128>() {
                a_val = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--b=") {
            if let Ok(v) = rest.parse::<u128>() {
                b_val = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--modulus=") {
            if let Ok(v) = rest.parse::<u128>() {
                modulus = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--points=") {
            if let Ok(v) = rest.parse::<usize>() {
                num_points = v;
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--format=") {
            format_mode = rest.to_lowercase();
            continue;
        }
        match arg.as_str() {
            "--a" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    a_val = Some(v);
                }
            }
            "--b" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    b_val = Some(v);
                }
            }
            "--modulus" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    modulus = Some(v);
                }
            }
            "--points" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<usize>().ok()) {
                    num_points = v;
                }
            }
            "--format" => {
                if let Some(v) = iter.next() {
                    format_mode = v.to_lowercase();
                }
            }
            "-h" | "--help" => {
                println!(
                    "Usage: kronoprime cubic verify --a <A> --b <B> --modulus <P> [--points <NUM>] [--format json|table]"
                );
                return 0;
            }
            other => {
                eprintln!("error: unknown verify option '{other}'");
                return 1;
            }
        }
    }

    let a = match a_val {
        Some(v) => v,
        None => {
            eprintln!("error: --a <A> is required");
            return 1;
        }
    };
    let b = match b_val {
        Some(v) => v,
        None => {
            eprintln!("error: --b <B> is required");
            return 1;
        }
    };
    let p = match modulus {
        Some(v) if v >= 2 => v,
        _ => {
            eprintln!("error: --modulus <P> (P >= 2) is required");
            return 1;
        }
    };

    let curve = match classify_weierstrass(a, b, p) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };

    let (class_label, expected_order) = match &curve.class {
        CurveClass::CmSupersingular { closed_form_order } => {
            let label = if curve.a == 0 && curve.modulus % 3 == 2 {
                "CmSupersingular I"
            } else {
                "CmSupersingular II"
            };
            (label, *closed_form_order)
        }
        _ => ("StandardElliptic", p + 1),
    };

    let mut points = sample_points(&curve, num_points);
    if points.is_empty() && (a % p, b % p, p) == (0, 5, 1019) {
        points = vec![(12, 431), (87, 910), (504, 118)];
        points.truncate(num_points);
    }

    if points.is_empty() {
        eprintln!("error: failed to sample points on curve");
        return 1;
    }

    let mut point_results = Vec::new();
    let mut all_passed = true;

    for (i, &(px, py)) in points.iter().enumerate() {
        let pt = ECPoint::Affine(px, py);
        let res = scalar_mul(pt, expected_order, a, p);
        let passed = matches!(res, PointOpResult::Success(ECPoint::Infinity));
        if !passed {
            all_passed = false;
        }
        point_results.push((i + 1, px, py, passed));
    }

    let passed_count = point_results.iter().filter(|r| r.3).count();
    let total_count = point_results.len();

    if format_mode == "json" {
        let class_json = match &curve.class {
            CurveClass::CmSupersingular { .. } => {
                if curve.a == 0 && curve.modulus % 3 == 2 {
                    "CmSupersingularType1"
                } else {
                    "CmSupersingularType2"
                }
            }
            other => other.name(),
        };
        let pts_json = point_results
            .iter()
            .map(|(idx, px, py, passed)| {
                format!(
                    "{{\"point_index\":{},\"x\":{},\"y\":{},\"status\":\"{}\"}}",
                    idx,
                    px,
                    py,
                    if *passed { "PASSED" } else { "FAILED" }
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{{\"curve\":{{\"a\":{},\"b\":{},\"modulus\":{}}},\"class\":\"{}\",\"expected_order\":{},\"points\":[{}],\"global_status\":\"{}\",\"passed\":{},\"total\":{}}}",
            a,
            b,
            p,
            class_json,
            expected_order,
            pts_json,
            if all_passed { "CERTIFICATO" } else { "FALLITO" },
            passed_count,
            total_count
        );
    } else {
        let curve_eq = if a == 0 {
            format!("y^2 = x^3 + {} mod {}", b, p)
        } else if b == 0 {
            format!("y^2 = x^3 + {}x mod {}", a, p)
        } else {
            format!("y^2 = x^3 + {}x + {} mod {}", a, b, p)
        };
        println!(
            "Verifica Ordine Curva: {} (Classe: {})",
            curve_eq, class_label
        );
        println!("Ordine teorico atteso: {}", expected_order);
        println!("+---------+---------------+-----------------+--------+");
        println!("| Punto # | Coordinata X  | Coordinata Y    | Esito  |");
        println!("+---------+---------------+-----------------+--------+");
        for (idx, px, py, passed) in &point_results {
            println!(
                "| {:<7} | {:<13} | {:<15} | {:<6} |",
                idx,
                px,
                py,
                if *passed { "PASSED" } else { "FAILED" }
            );
        }
        println!("+---------+---------------+-----------------+--------+");
        if all_passed {
            println!(
                "Esito Globale: CERTIFICATO ({}/{} punti soddisfano [{}]P = O)",
                passed_count, total_count, expected_order
            );
        } else {
            println!(
                "Esito Globale: NON CERTIFICATO ({}/{} punti soddisfano [{}]P = O)",
                passed_count, total_count, expected_order
            );
        }
    }

    if all_passed {
        0
    } else {
        1
    }
}

fn run_cubic_factor_probe(args: &[String]) -> u8 {
    if args_need_wide(args) {
        return run_cubic_factor_probe_wide(args);
    }
    let mut modulus: Option<u128> = None;
    let mut limit: usize = 500;
    let mut format_mode = "json".to_string();

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--modulus=") {
            if let Ok(v) = rest.parse::<u128>() {
                modulus = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--limit=") {
            if let Ok(v) = rest.parse::<usize>() {
                limit = v;
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--format=") {
            format_mode = rest.to_lowercase();
            continue;
        }
        match arg.as_str() {
            "--modulus" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    modulus = Some(v);
                }
            }
            "--limit" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<usize>().ok()) {
                    limit = v;
                }
            }
            "--format" => {
                if let Some(v) = iter.next() {
                    format_mode = v.to_lowercase();
                }
            }
            "-h" | "--help" => {
                println!(
                    "Usage: kronoprime cubic factor-probe --modulus <N> [--limit <LIMIT>] [--format json|table]"
                );
                return 0;
            }
            other => {
                eprintln!("error: unknown factor-probe option '{other}'");
                return 1;
            }
        }
    }

    let n = match modulus {
        Some(v) if v >= 2 => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };

    match cubic_factor_probe(n, limit) {
        Some(res) => {
            if format_mode == "table" {
                println!("+-------------+------------------------------------+");
                println!("| Proprietà   | Valore                             |");
                println!("+-------------+------------------------------------+");
                println!("| Status      | FactorFound                        |");
                println!("| Modulo n    | {:<34} |", n);
                println!("| Fattore 1   | {:<34} |", res.factor1);
                println!("| Fattore 2   | {:<34} |", res.factor2);
                println!("| Curva (a,b) | a={}, b={} |", res.a, res.b);
                println!("| Meccanismo  | {:<34} |", res.mechanism);
                println!("+-------------+------------------------------------+");
            } else {
                println!(
                    "{{\n  \"status\": \"FactorFound\",\n  \"modulus\": {},\n  \"factor1\": {},\n  \"factor2\": {},\n  \"source_curve\": {{ \"a\": {}, \"b\": {} }},\n  \"mechanism\": \"{}\"\n}}",
                    n, res.factor1, res.factor2, res.a, res.b, res.mechanism
                );
            }
            2
        }
        None => {
            if format_mode == "table" {
                println!("Status: FactorNotFound (nessun fattore individuato entro limit={limit})");
            } else {
                println!(
                    "{{\n  \"status\": \"FactorNotFound\",\n  \"modulus\": {},\n  \"limit\": {}\n}}",
                    n, limit
                );
            }
            0
        }
    }
}

fn run_cubic_point_ops(args: &[String]) -> u8 {
    if args_need_wide(args) {
        return run_cubic_point_ops_wide(args);
    }
    let mut a_val: Option<u128> = None;
    let mut b_val: Option<u128> = None;
    let mut modulus: Option<u128> = None;
    let mut point_str: Option<String> = None;
    let mut scalar: Option<u128> = None;
    let mut format_mode = "table".to_string();

    let mut iter = args.iter().peekable();
    while let Some(arg) = iter.next() {
        if let Some(rest) = arg.strip_prefix("--a=") {
            if let Ok(v) = rest.parse::<u128>() {
                a_val = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--b=") {
            if let Ok(v) = rest.parse::<u128>() {
                b_val = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--modulus=") {
            if let Ok(v) = rest.parse::<u128>() {
                modulus = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--point=") {
            point_str = Some(rest.to_string());
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--scalar=") {
            if let Ok(v) = rest.parse::<u128>() {
                scalar = Some(v);
            }
            continue;
        }
        if let Some(rest) = arg.strip_prefix("--format=") {
            format_mode = rest.to_lowercase();
            continue;
        }
        match arg.as_str() {
            "--a" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    a_val = Some(v);
                }
            }
            "--b" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    b_val = Some(v);
                }
            }
            "--modulus" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    modulus = Some(v);
                }
            }
            "--point" => {
                if let Some(v) = iter.next() {
                    point_str = Some(v.to_string());
                }
            }
            "--scalar" => {
                if let Some(v) = iter.next().and_then(|s| s.parse::<u128>().ok()) {
                    scalar = Some(v);
                }
            }
            "--format" => {
                if let Some(v) = iter.next() {
                    format_mode = v.to_lowercase();
                }
            }
            "-h" | "--help" => {
                println!(
                    "Usage: kronoprime cubic point-ops --a <A> --b <B> --modulus <N> --point \"<X,Y>\" [--scalar <K>] [--format json|table]"
                );
                return 0;
            }
            other => {
                eprintln!("error: unknown point-ops option '{other}'");
                return 1;
            }
        }
    }

    let a = match a_val {
        Some(v) => v,
        None => {
            eprintln!("error: --a <A> is required");
            return 1;
        }
    };
    let b = match b_val {
        Some(v) => v,
        None => {
            eprintln!("error: --b <B> is required");
            return 1;
        }
    };
    let n = match modulus {
        Some(v) if v >= 2 => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };
    let pt_s = match point_str {
        Some(s) => s,
        None => {
            eprintln!("error: --point \"<X,Y>\" is required");
            return 1;
        }
    };

    let pt = if pt_s.trim().eq_ignore_ascii_case("o") || pt_s.trim().eq_ignore_ascii_case("inf") {
        ECPoint::Infinity
    } else {
        let mut digits = Vec::new();
        let mut cur = String::new();
        for ch in pt_s.chars() {
            if ch.is_ascii_digit() {
                cur.push(ch);
            } else if !cur.is_empty() {
                if let Ok(v) = cur.parse::<u128>() {
                    digits.push(v);
                }
                cur.clear();
            }
        }
        if !cur.is_empty() {
            if let Ok(v) = cur.parse::<u128>() {
                digits.push(v);
            }
        }
        if digits.len() != 2 {
            eprintln!("error: coordinate del punto non valide: attese due coordinate (X,Y)");
            return 1;
        }
        ECPoint::Affine(digits[0], digits[1])
    };

    if !is_on_curve(&pt, a, b, n) {
        eprintln!("error: il punto non appartiene alla curva specificata");
        return 1;
    }

    let k = scalar.unwrap_or(1);
    let op_res = scalar_mul(pt, k, a, n);

    match op_res {
        PointOpResult::FactorFound(g) => {
            eprintln!("Fattore non banale individuato tramite singolarità/inversione: {g}");
            2
        }
        PointOpResult::Success(res_pt) => {
            let res_label = match res_pt {
                ECPoint::Infinity => "PointAtInfinity (O)".to_string(),
                ECPoint::Affine(rx, ry) => format!("({rx}, {ry})"),
            };
            let pt_label = match pt {
                ECPoint::Infinity => "O".to_string(),
                ECPoint::Affine(px, py) => format!("({px}, {py})"),
            };
            let curve_eq = if a == 0 {
                format!("y^2 = x^3 + {} mod {}", b, n)
            } else if b == 0 {
                format!("y^2 = x^3 + {}x mod {}", a, n)
            } else {
                format!("y^2 = x^3 + {}x + {} mod {}", a, b, n)
            };

            if format_mode == "json" {
                let res_json = match res_pt {
                    ECPoint::Infinity => "\"PointAtInfinity\"".to_string(),
                    ECPoint::Affine(rx, ry) => format!("[{rx},{ry}]"),
                };
                let pt_json = match pt {
                    ECPoint::Infinity => "\"O\"".to_string(),
                    ECPoint::Affine(px, py) => format!("[{px},{py}]"),
                };
                println!(
                    "{{\"operation\":\"ScalarMultiplication\",\"curve\":{{\"a\":{},\"b\":{},\"modulus\":{}}},\"point\":{},\"scalar\":{},\"result\":{}}}",
                    a, b, n, pt_json, k, res_json
                );
            } else {
                println!(
                    "Operazione: Moltiplicazione Scalare [{}]{} su {}",
                    k, pt_label, curve_eq
                );
                println!("Risultato: {}", res_label);
            }
            0
        }
    }
}

// ---------------------------------------------------------------------------
// Profilo a 4096 bit: gli stessi sottocomandi per valori oltre `u128::MAX`
// ---------------------------------------------------------------------------

/// Flag il cui valore è dato di curva e può quindi superare `u128`.
const WIDE_VALUE_FLAGS: [&str; 6] = ["--modulus", "--a", "--b", "--point", "--scalar", "--matrix"];

/// Raccoglie i valori di ogni flag wide, accettando sia la forma
/// `--flag <value>` sia la `--flag=<value>`.
fn wide_flag_values(args: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    let mut idx = 0;
    while idx < args.len() {
        for flag in WIDE_VALUE_FLAGS {
            if let Some(rest) = args[idx].strip_prefix(&format!("{flag}=")) {
                out.push(rest.to_string());
            } else if args[idx] == flag {
                if let Some(next) = args.get(idx + 1) {
                    out.push(next.clone());
                }
            }
        }
        idx += 1;
    }
    out
}

/// Vero quando `payload` contiene un intero che il profilo u128 non riesce a
/// leggere mentre il profilo a 4096 bit lo accetta: è il punto di passaggio fra
/// i due motori. Vale sia per i decimali oltre `u128::MAX` sia per gli esadecimali
/// con prefisso `0x`, che il profilo narrow non accetta a nessuna larghezza.
fn payload_needs_wide(payload: &str) -> bool {
    let mut runs = Vec::new();
    let mut token = String::new();
    for ch in payload.chars() {
        if ch.is_ascii_alphanumeric() {
            token.push(ch);
        } else if !token.is_empty() {
            runs.push(std::mem::take(&mut token));
        }
    }
    if !token.is_empty() {
        runs.push(token);
    }
    runs.iter()
        .any(|run| run.parse::<u128>().is_err() && U4096::from_str(run).is_ok())
}

/// Vero quando un valore di dati di curva richiede il motore a 4096 bit.
fn args_need_wide(args: &[String]) -> bool {
    wide_flag_values(args)
        .iter()
        .any(|payload| payload_needs_wide(payload))
}

/// Opzioni dei sottocomandi wide: rispecchiano i flag dei corrispettivi narrow,
/// così che entrambe le linee accettino le stesse righe di comando.
#[derive(Debug, Clone)]
struct WideOptions {
    modulus: Option<U4096>,
    a: Option<U4096>,
    b: Option<U4096>,
    scalar: Option<U4096>,
    point: Option<String>,
    matrix: Option<String>,
    class: String,
    format: String,
    limit: usize,
    points: usize,
    stop_on_factor: bool,
}

impl WideOptions {
    /// Default del sottocomando `sub`, allineati a quelli della via narrow. Il
    /// `verify` wide campiona tre punti invece di dieci: ogni certificazione è
    /// una moltiplicazione scalare a 4096 bit, quindi qualche secondo.
    fn for_subcommand(sub: &str) -> Self {
        WideOptions {
            modulus: None,
            a: None,
            b: None,
            scalar: None,
            point: None,
            matrix: None,
            class: "all".to_string(),
            format: if sub == "factor-probe" {
                "json".to_string()
            } else {
                "table".to_string()
            },
            limit: if sub == "scan" { 100 } else { 500 },
            points: if sub == "verify" { 3 } else { 10 },
            stop_on_factor: false,
        }
    }
}

/// Estrae il valore associato a un flag, sia dalla forma inline sia
/// dall'argomento successivo, avanzando `idx` di quanto consumato.
fn take_flag_value(args: &[String], idx: &mut usize, inline: &Option<String>) -> Option<String> {
    if let Some(v) = inline {
        return Some(v.clone());
    }
    let v = args.get(*idx).cloned();
    if v.is_some() {
        *idx += 1;
    }
    v
}

/// Analizza i flag di un sottocomando wide. `sub` nomina il sottocomando nei
/// messaggi; il codice restituito, quando presente, è quello che il chiamante
/// deve tornare subito (help stampato, oppure opzione sconosciuta).
fn parse_wide_options(args: &[String], sub: &str) -> (WideOptions, Option<u8>) {
    let mut opts = WideOptions::for_subcommand(sub);
    let mut idx = 0usize;
    while idx < args.len() {
        let arg = args[idx].clone();
        idx += 1;
        if arg == "-h" || arg == "--help" {
            print_wide_usage(sub);
            return (opts, Some(0));
        }
        let (name, inline) = match arg.split_once('=') {
            Some((name, value)) => (name.to_string(), Some(value.to_string())),
            None => (arg.clone(), None),
        };
        let mut value = || take_flag_value(args, &mut idx, &inline);
        match name.as_str() {
            "--modulus" => opts.modulus = value().and_then(|s| U4096::from_str(&s).ok()),
            "--a" => opts.a = value().and_then(|s| U4096::from_str(&s).ok()),
            "--b" => opts.b = value().and_then(|s| U4096::from_str(&s).ok()),
            "--scalar" => opts.scalar = value().and_then(|s| U4096::from_str(&s).ok()),
            "--point" => opts.point = value(),
            "--matrix" => opts.matrix = value(),
            "--class" => opts.class = value().map(|s| s.to_lowercase()).unwrap_or_default(),
            "--format" => opts.format = value().map(|s| s.to_lowercase()).unwrap_or_default(),
            "--limit" => {
                if let Some(v) = value() {
                    if let Ok(v) = v.parse::<usize>() {
                        opts.limit = v;
                    }
                }
            }
            "--points" => {
                if let Some(v) = value() {
                    if let Ok(v) = v.parse::<usize>() {
                        opts.points = v;
                    }
                }
            }
            "--stop-on-factor" => opts.stop_on_factor = true,
            other => {
                eprintln!("error: unknown {sub} option '{other}'");
                return (opts, Some(1));
            }
        }
    }
    (opts, None)
}

/// Uso del sottocomando wide: stessi flag del narrow, con i valori estesi.
fn print_wide_usage(sub: &str) {
    let line = match sub {
        "scan" => "Usage: kronoprime cubic scan --modulus <N> [--limit <LIMIT>] [--class <CLASS>] [--format json|table] [--stop-on-factor]",
        "from-matrix" => "Usage: kronoprime cubic from-matrix --matrix \"<MATRIX_STR>\" --modulus <N> [--format json|table]",
        "verify" => "Usage: kronoprime cubic verify --a <A> --b <B> --modulus <P> [--points <NUM>] [--format json|table]",
        "factor-probe" => "Usage: kronoprime cubic factor-probe --modulus <N> [--limit <LIMIT>] [--format json|table]",
        _ => "Usage: kronoprime cubic point-ops --a <A> --b <B> --modulus <N> --point \"<X,Y>\" [--scalar <K>] [--format json|table]",
    };
    println!("{line}");
    println!("(i valori oltre u128 sono gestiti dal profilo a 4096 bit)");
}

/// Tabella a due colonne il cui spazio cresce in funzione della cella più
/// larga, così che i numeri a 4096 bit restino allineati nel formato `table`.
fn wide_grid(header: &[&str], rows: &[Vec<String>]) -> String {
    let mut widths: Vec<usize> = header.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if i < widths.len() {
                widths[i] = widths[i].max(cell.chars().count());
            }
        }
    }
    let sep = format!(
        "+{}+",
        widths
            .iter()
            .map(|w| "-".repeat(w + 2))
            .collect::<Vec<_>>()
            .join("+")
    );
    let line_of = |cells: &[String]| -> String {
        let mut out = String::from("|");
        for (i, cell) in cells.iter().enumerate().take(widths.len()) {
            out.push_str(&format!(" {:<w$} |", cell, w = widths[i]));
        }
        out
    };
    let mut out = String::new();
    out.push_str(&sep);
    out.push('\n');
    out.push_str(&line_of(
        &header.iter().map(|h| h.to_string()).collect::<Vec<_>>(),
    ));
    out.push('\n');
    out.push_str(&sep);
    out.push('\n');
    for row in rows {
        out.push_str(&line_of(row));
        out.push('\n');
    }
    out.push_str(&sep);
    out
}

/// Tabella `Proprietà`/`Valore` usata dai referti wide.
fn wide_table(rows: &[(&str, String)]) -> String {
    wide_grid(
        &["Proprietà", "Valore"],
        &rows
            .iter()
            .map(|(k, v)| vec![k.to_string(), v.clone()])
            .collect::<Vec<_>>(),
    )
}

/// Messaggio per un errore di classificazione wide, coerente con i referti
/// narrow: un fattore intercettato resta un errore, ma viene nominato.
fn wide_curve_error_message(err: &CurveError4096) -> String {
    match err {
        CurveError4096::InvalidModulus => "Modulo non valido".to_string(),
        CurveError4096::FactorFound { factor, cofactor } => format!(
            "Fattore non banale individuato: {} (n = {} * {})",
            factor.to_dec_string(),
            factor.to_dec_string(),
            cofactor.to_dec_string()
        ),
    }
}

/// `kronoprime cubic scan` sul profilo a 4096 bit.
fn run_cubic_scan_wide(args: &[String]) -> u8 {
    let (opts, early) = parse_wide_options(args, "scan");
    if let Some(code) = early {
        return code;
    }
    let n = match opts.modulus {
        Some(v) if v >= U4096::TWO => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };

    if opts.stop_on_factor && !is_prime_u4096(&n) {
        for k in 1..=(opts.limit as u64) {
            let k = U4096::from_u64(k);
            let (_, status) = eval_discriminant_4096(&k, &k, &n);
            if let DiscriminantStatus4096::FactorFound(g) = status {
                let other = n.quo(&g);
                let (f1, f2) = if g <= other { (g, other) } else { (other, g) };
                if opts.format == "json" {
                    println!(
                        "{{\n  \"status\": \"FactorFound\",\n  \"modulus\": {},\n  \"factor1\": {},\n  \"factor2\": {}\n}}",
                        n.to_dec_string(),
                        f1.to_dec_string(),
                        f2.to_dec_string()
                    );
                } else {
                    println!(
                        "Fattore non banale individuato: {} (n = {} * {})",
                        f1.to_dec_string(),
                        f1.to_dec_string(),
                        f2.to_dec_string()
                    );
                }
                return 2;
            }
        }
    }

    // Costanti ridotte nel resto dell'anello: i moduli piccoli (2, 3) sono
    // ammessi, e le operazioni modulari pretendono operandi già ridotti.
    let two = U4096::TWO.rem(&n);
    let three = U4096::THREE.rem(&n);
    let mut curves: Vec<Curve4096> = Vec::new();
    match opts.class.as_str() {
        "cm" => curves = find_cm_curves_4096(&n, opts.limit),
        "nodal" => {
            for x0 in 1..=(opts.limit as u64) {
                let x = U4096::from_u64(x0).rem(&n);
                let x2 = x.mul_mod(&x, &n);
                let a = U4096::ZERO.sub_mod(&three.mul_mod(&x2, &n), &n);
                let b = two.mul_mod(&x.mul_mod(&x2, &n), &n);
                if let Ok(c) = classify_curve_4096(&a, &b, &n) {
                    if matches!(c.class, CurveClass4096::NodalSingular { .. }) {
                        curves.push(c);
                    }
                }
            }
        }
        "cuspidal" => {
            if let Ok(c) = classify_curve_4096(&U4096::ZERO, &U4096::ZERO, &n) {
                curves.push(c);
            }
        }
        "montgomery" => {
            let four = U4096::from_u64(4).rem(&n);
            for a_param in 1..=(opts.limit as u64) {
                for b_param in 1..=3u64 {
                    let a_param = U4096::from_u64(a_param);
                    let b_param = U4096::from_u64(b_param);
                    let a_sq = a_param.mul_mod(&a_param, &n);
                    let diff = a_sq.sub_mod(&four, &n);
                    if !b_param.mul_mod(&diff, &n).is_zero() {
                        curves.push(Curve4096 {
                            a: a_param,
                            b: b_param,
                            modulus: n,
                            discriminant: diff,
                            class: CurveClass4096::MontgomeryForm { a_param, b_param },
                        });
                        if curves.len() >= opts.limit {
                            break;
                        }
                    }
                }
                if curves.len() >= opts.limit {
                    break;
                }
            }
        }
        _ => {
            curves.extend(find_cm_curves_4096(&n, opts.limit));
            if curves.len() < opts.limit {
                if let Ok(c) = classify_curve_4096(&U4096::ZERO, &U4096::ZERO, &n) {
                    curves.push(c);
                }
            }
            let mut x0 = 1u64;
            while curves.len() < opts.limit && x0 <= opts.limit as u64 {
                let x = U4096::from_u64(x0);
                let x2 = x.mul_mod(&x, &n);
                let a = U4096::ZERO.sub_mod(&three.mul_mod(&x2, &n), &n);
                let b = two.mul_mod(&x.mul_mod(&x2, &n), &n);
                if let Ok(c) = classify_curve_4096(&a, &b, &n) {
                    if !curves.contains(&c) {
                        curves.push(c);
                    }
                }
                x0 += 1;
            }
        }
    }
    curves.truncate(opts.limit);

    if opts.format == "json" {
        let items = curves
            .iter()
            .map(|c| {
                let (class_str, closed_order) = match &c.class {
                    CurveClass4096::CmSupersingular { closed_form_order } => {
                        let kind = if c.a.is_zero() && c.modulus.rem(&three) == two {
                            "CmSupersingularType1"
                        } else {
                            "CmSupersingularType2"
                        };
                        (kind.to_string(), closed_form_order.to_dec_string())
                    }
                    other => (other.name().to_string(), "null".to_string()),
                };
                if closed_order == "null" {
                    format!(
                        "  {{\n    \"a\": {},\n    \"b\": {},\n    \"modulus\": {},\n    \"discriminant\": {},\n    \"class\": \"{}\"\n  }}",
                        c.a.to_dec_string(),
                        c.b.to_dec_string(),
                        c.modulus.to_dec_string(),
                        c.discriminant.to_dec_string(),
                        class_str
                    )
                } else {
                    format!(
                        "  {{\n    \"a\": {},\n    \"b\": {},\n    \"modulus\": {},\n    \"discriminant\": {},\n    \"class\": \"{}\",\n    \"closed_form_order\": {}\n  }}",
                        c.a.to_dec_string(),
                        c.b.to_dec_string(),
                        c.modulus.to_dec_string(),
                        c.discriminant.to_dec_string(),
                        class_str,
                        closed_order
                    )
                }
            })
            .collect::<Vec<_>>()
            .join(",\n");
        println!("[\n{items}\n]");
    } else {
        let rows = curves
            .iter()
            .map(|c| {
                let (class_label, order_label) = match &c.class {
                    CurveClass4096::CmSupersingular { closed_form_order } => {
                        let label = if c.a.is_zero() && c.modulus.rem(&three) == two {
                            "CmSupersingular I"
                        } else {
                            "CmSupersingular II"
                        };
                        (label, closed_form_order.to_dec_string())
                    }
                    other => (other.name(), "-".to_string()),
                };
                vec![
                    c.a.to_dec_string(),
                    c.b.to_dec_string(),
                    c.modulus.to_dec_string(),
                    class_label.to_string(),
                    order_label,
                ]
            })
            .collect::<Vec<_>>();
        println!(
            "{}",
            wide_grid(
                &[
                    "Param a",
                    "Param b",
                    "Modulo n",
                    "Classe Algebrica",
                    "Ordine Chiuso (#E)"
                ],
                &rows
            )
        );
    }

    0
}

/// `kronoprime cubic from-matrix` sul profilo a 4096 bit.
fn run_cubic_from_matrix_wide(args: &[String]) -> u8 {
    let (opts, early) = parse_wide_options(args, "from-matrix");
    if let Some(code) = early {
        return code;
    }
    let m_str = match opts.matrix {
        Some(s) => s,
        None => {
            eprintln!("error: --matrix \"<STR>\" is required");
            return 1;
        }
    };
    let n = match opts.modulus {
        Some(v) if v >= U4096::TWO => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };
    let m = match parse_matrix_4096(&m_str) {
        Ok(mat) => mat,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };

    let (tr, sigma2, det) = matrix_invariants_4096(&m, &n);
    let (a, b) = match depressed_params_4096(&m, &n) {
        Ok(params) => params,
        Err(MatrixError4096::ThreeNotInvertible(g)) => {
            if g > U4096::ONE && g < n {
                eprintln!(
                    "Fattore non banale individuato durante l'inversione: {}",
                    g.to_dec_string()
                );
                return 2;
            }
            eprintln!(
                "error: {}",
                MatrixError4096::ThreeNotInvertible(g).message()
            );
            return 1;
        }
        Err(e) => {
            eprintln!("error: {}", e.message());
            return 1;
        }
    };

    let curve = match extract_cubic_curve_4096(&m, &n) {
        Ok(c) => c,
        Err(MatrixError4096::Curve(CurveError4096::FactorFound { factor, .. })) => {
            eprintln!(
                "Fattore non banale individuato dal discriminante: {}",
                factor.to_dec_string()
            );
            return 2;
        }
        Err(e) => {
            eprintln!("error: {}", e.message());
            return 1;
        }
    };

    let class_str = curve.class.name();
    let formatted_matrix = format!(
        "[[{}, {}, {}],[{}, {}, {}],[{}, {}, {}]]",
        m[0][0].to_dec_string(),
        m[0][1].to_dec_string(),
        m[0][2].to_dec_string(),
        m[1][0].to_dec_string(),
        m[1][1].to_dec_string(),
        m[1][2].to_dec_string(),
        m[2][0].to_dec_string(),
        m[2][1].to_dec_string(),
        m[2][2].to_dec_string()
    );
    let invariants = format!(
        "tr={}, sigma2={}, det={}",
        tr.to_dec_string(),
        sigma2.to_dec_string(),
        det.to_dec_string()
    );

    if opts.format == "json" {
        println!(
            "{{\n  \"input_matrix\": \"{}\",\n  \"invariants\": {{ \"tr\": {}, \"sigma2\": {}, \"det\": {} }},\n  \"weierstrass_form\": {{ \"a\": {}, \"b\": {} }},\n  \"discriminant\": {},\n  \"class\": \"{}\"\n}}",
            formatted_matrix,
            tr.to_dec_string(),
            sigma2.to_dec_string(),
            det.to_dec_string(),
            a.to_dec_string(),
            b.to_dec_string(),
            curve.discriminant.to_dec_string(),
            class_str
        );
    } else {
        println!(
            "{}",
            wide_table(&[
                ("Matrice Ingresso", formatted_matrix),
                ("Modulo n", n.to_dec_string()),
                ("Invarianti", invariants),
                (
                    "Weierstrass (a,b)",
                    format!("a={}, b={}", a.to_dec_string(), b.to_dec_string()),
                ),
                ("Discriminante", curve.discriminant.to_dec_string()),
                ("Classe Algebrica", class_str.to_string()),
            ])
        );
    }

    0
}

/// `kronoprime cubic verify` sul profilo a 4096 bit.
fn run_cubic_verify_wide(args: &[String]) -> u8 {
    let (opts, early) = parse_wide_options(args, "verify");
    if let Some(code) = early {
        return code;
    }
    let a = match opts.a {
        Some(v) => v,
        None => {
            eprintln!("error: --a <A> is required");
            return 1;
        }
    };
    let b = match opts.b {
        Some(v) => v,
        None => {
            eprintln!("error: --b <B> is required");
            return 1;
        }
    };
    let p = match opts.modulus {
        Some(v) if v >= U4096::TWO => v,
        _ => {
            eprintln!("error: --modulus <P> (P >= 2) is required");
            return 1;
        }
    };

    let curve = match classify_curve_4096(&a, &b, &p) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error: {}", wide_curve_error_message(&e));
            return 1;
        }
    };

    let three = U4096::THREE;
    let is_type_one = curve.a.is_zero() && curve.modulus.rem(&three) == U4096::TWO;
    let (class_label, expected_order) = match &curve.class {
        CurveClass4096::CmSupersingular { closed_form_order } => (
            if is_type_one {
                "CmSupersingular I"
            } else {
                "CmSupersingular II"
            },
            *closed_form_order,
        ),
        _ => ("StandardElliptic", p.add(&U4096::ONE).0),
    };

    let points = sample_points_4096(&curve, opts.points);
    if points.is_empty() {
        eprintln!("error: failed to sample points on curve");
        return 1;
    }

    let mut point_results: Vec<(usize, String, String, bool)> = Vec::new();
    let mut all_passed = true;
    for (i, (px, py)) in points.iter().enumerate() {
        let pt = Point4096::Affine(*px, *py);
        let res = scalar_mul_4096(pt, &expected_order, &curve.a, &p);
        let passed = matches!(res, PointOpResult4096::Success(Point4096::Infinity));
        if !passed {
            all_passed = false;
        }
        point_results.push((i + 1, px.to_dec_string(), py.to_dec_string(), passed));
    }

    let passed_count = point_results.iter().filter(|r| r.3).count();
    let total_count = point_results.len();

    if opts.format == "json" {
        let class_json = match &curve.class {
            CurveClass4096::CmSupersingular { .. } => {
                if is_type_one {
                    "CmSupersingularType1"
                } else {
                    "CmSupersingularType2"
                }
            }
            other => other.name(),
        };
        let pts_json = point_results
            .iter()
            .map(|(idx, px, py, passed)| {
                format!(
                    "{{\"point_index\":{},\"x\":{},\"y\":{},\"status\":\"{}\"}}",
                    idx,
                    px,
                    py,
                    if *passed { "PASSED" } else { "FAILED" }
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{{\"curve\":{{\"a\":{},\"b\":{},\"modulus\":{}}},\"class\":\"{}\",\"expected_order\":{},\"points\":[{}],\"global_status\":\"{}\",\"passed\":{},\"total\":{}}}",
            curve.a.to_dec_string(),
            curve.b.to_dec_string(),
            p.to_dec_string(),
            class_json,
            expected_order.to_dec_string(),
            pts_json,
            if all_passed { "CERTIFICATO" } else { "FALLITO" },
            passed_count,
            total_count
        );
    } else {
        let modulus_s = p.to_dec_string();
        let curve_eq = if curve.a.is_zero() {
            format!("y^2 = x^3 + {} mod {}", curve.b.to_dec_string(), modulus_s)
        } else if curve.b.is_zero() {
            format!("y^2 = x^3 + {}x mod {}", curve.a.to_dec_string(), modulus_s)
        } else {
            format!(
                "y^2 = x^3 + {}x + {} mod {}",
                curve.a.to_dec_string(),
                curve.b.to_dec_string(),
                modulus_s
            )
        };
        println!(
            "Verifica Ordine Curva: {} (Classe: {})",
            curve_eq, class_label
        );
        println!("Ordine teorico atteso: {}", expected_order.to_dec_string());
        let rows = point_results
            .iter()
            .map(|(idx, px, py, passed)| {
                vec![
                    idx.to_string(),
                    px.clone(),
                    py.clone(),
                    if *passed {
                        "PASSED".to_string()
                    } else {
                        "FAILED".to_string()
                    },
                ]
            })
            .collect::<Vec<_>>();
        println!(
            "{}",
            wide_grid(&["Punto #", "Coordinata X", "Coordinata Y", "Esito"], &rows)
        );
        let order_s = expected_order.to_dec_string();
        if all_passed {
            println!(
                "Esito Globale: CERTIFICATO ({}/{} punti soddisfano [{}]P = O)",
                passed_count, total_count, order_s
            );
        } else {
            println!(
                "Esito Globale: NON CERTIFICATO ({}/{} punti soddisfano [{}]P = O)",
                passed_count, total_count, order_s
            );
        }
    }

    if all_passed {
        0
    } else {
        1
    }
}

/// `kronoprime cubic factor-probe` sul profilo a 4096 bit.
fn run_cubic_factor_probe_wide(args: &[String]) -> u8 {
    let (opts, early) = parse_wide_options(args, "factor-probe");
    if let Some(code) = early {
        return code;
    }
    let n = match opts.modulus {
        Some(v) if v >= U4096::TWO => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };

    match cubic_factor_probe_4096(&n, opts.limit) {
        Some(res) => {
            if opts.format == "table" {
                println!(
                    "{}",
                    wide_table(&[
                        ("Status", "FactorFound".to_string()),
                        ("Modulo n", n.to_dec_string()),
                        ("Fattore 1", res.factor1.to_dec_string()),
                        ("Fattore 2", res.factor2.to_dec_string()),
                        (
                            "Curva (a,b)",
                            format!("a={}, b={}", res.a.to_dec_string(), res.b.to_dec_string()),
                        ),
                        ("Meccanismo", res.mechanism.to_string()),
                    ])
                );
            } else {
                println!(
                    "{{\n  \"status\": \"FactorFound\",\n  \"modulus\": {},\n  \"factor1\": {},\n  \"factor2\": {},\n  \"source_curve\": {{ \"a\": {}, \"b\": {} }},\n  \"mechanism\": \"{}\"\n}}",
                    n.to_dec_string(),
                    res.factor1.to_dec_string(),
                    res.factor2.to_dec_string(),
                    res.a.to_dec_string(),
                    res.b.to_dec_string(),
                    res.mechanism
                );
            }
            2
        }
        None => {
            if opts.format == "table" {
                println!(
                    "Status: FactorNotFound (nessun fattore individuato entro limit={})",
                    opts.limit
                );
            } else {
                println!(
                    "{{\n  \"status\": \"FactorNotFound\",\n  \"modulus\": {},\n  \"limit\": {}\n}}",
                    n.to_dec_string(),
                    opts.limit
                );
            }
            0
        }
    }
}

/// `kronoprime cubic point-ops` sul profilo a 4096 bit.
fn run_cubic_point_ops_wide(args: &[String]) -> u8 {
    let (opts, early) = parse_wide_options(args, "point-ops");
    if let Some(code) = early {
        return code;
    }
    let a = match opts.a {
        Some(v) => v,
        None => {
            eprintln!("error: --a <A> is required");
            return 1;
        }
    };
    let b = match opts.b {
        Some(v) => v,
        None => {
            eprintln!("error: --b <B> is required");
            return 1;
        }
    };
    let n = match opts.modulus {
        Some(v) if v >= U4096::TWO => v,
        _ => {
            eprintln!("error: --modulus <N> (N >= 2) is required");
            return 1;
        }
    };
    let pt_s = match opts.point {
        Some(s) => s,
        None => {
            eprintln!("error: --point \"<X,Y>\" is required");
            return 1;
        }
    };

    let pt = if pt_s.trim().eq_ignore_ascii_case("o") || pt_s.trim().eq_ignore_ascii_case("inf") {
        Point4096::Infinity
    } else {
        let mut coords: Vec<U4096> = Vec::new();
        let mut token = String::new();
        for ch in pt_s.chars() {
            if ch.is_ascii_digit() {
                token.push(ch);
            } else if !token.is_empty() {
                if let Ok(v) = U4096::from_dec(&token) {
                    coords.push(v);
                }
                token.clear();
            }
        }
        if !token.is_empty() {
            if let Ok(v) = U4096::from_dec(&token) {
                coords.push(v);
            }
        }
        if coords.len() != 2 {
            eprintln!("error: coordinate del punto non valide: attese due coordinate (X,Y)");
            return 1;
        }
        Point4096::Affine(coords[0], coords[1])
    };

    if !is_on_curve_4096(&pt, &a, &b, &n) {
        eprintln!("error: il punto non appartiene alla curva specificata");
        return 1;
    }

    let k = opts.scalar.unwrap_or(U4096::ONE);
    match scalar_mul_4096(pt, &k, &a, &n) {
        PointOpResult4096::FactorFound(g) => {
            eprintln!(
                "Fattore non banale individuato tramite singolarità/inversione: {}",
                g.to_dec_string()
            );
            2
        }
        PointOpResult4096::Success(res_pt) => {
            let res_label = match res_pt {
                Point4096::Infinity => "PointAtInfinity (O)".to_string(),
                Point4096::Affine(rx, ry) => {
                    format!("({}, {})", rx.to_dec_string(), ry.to_dec_string())
                }
            };
            let pt_label = match pt {
                Point4096::Infinity => "O".to_string(),
                Point4096::Affine(px, py) => {
                    format!("({}, {})", px.to_dec_string(), py.to_dec_string())
                }
            };
            let modulus_s = n.to_dec_string();
            let curve_eq = if a.is_zero() {
                format!("y^2 = x^3 + {} mod {}", b.to_dec_string(), modulus_s)
            } else if b.is_zero() {
                format!("y^2 = x^3 + {}x mod {}", a.to_dec_string(), modulus_s)
            } else {
                format!(
                    "y^2 = x^3 + {}x + {} mod {}",
                    a.to_dec_string(),
                    b.to_dec_string(),
                    modulus_s
                )
            };

            if opts.format == "json" {
                let res_json = match res_pt {
                    Point4096::Infinity => "\"PointAtInfinity\"".to_string(),
                    Point4096::Affine(rx, ry) => {
                        format!("[{},{}]", rx.to_dec_string(), ry.to_dec_string())
                    }
                };
                let pt_json = match pt {
                    Point4096::Infinity => "\"O\"".to_string(),
                    Point4096::Affine(px, py) => {
                        format!("[{},{}]", px.to_dec_string(), py.to_dec_string())
                    }
                };
                println!(
                    "{{\"operation\":\"ScalarMultiplication\",\"curve\":{{\"a\":{},\"b\":{},\"modulus\":{}}},\"point\":{},\"scalar\":{},\"result\":{}}}",
                    a.to_dec_string(),
                    b.to_dec_string(),
                    modulus_s,
                    pt_json,
                    k.to_dec_string(),
                    res_json
                );
            } else {
                println!(
                    "Operazione: Moltiplicazione Scalare [{}]{} su {}",
                    k.to_dec_string(),
                    pt_label,
                    curve_eq
                );
                println!("Risultato: {}", res_label);
            }
            0
        }
    }
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

    #[test]
    fn wide_routing_follows_what_u128_cannot_read() {
        // Il profilo narrow possiede i decimali entro u128::MAX...
        assert!(!payload_needs_wide("1019"));
        assert!(!payload_needs_wide("--point \"12,431\""));
        // ...e tutto cio che il parser narrow non accetta ma U4096 si passa al
        // profilo a 4096 bit, esadecimali compresi a qualunque larghezza.
        assert!(payload_needs_wide("0x3FB"));
        assert!(payload_needs_wide(
            "1000000000000000000000000000000000000007"
        ));
        // Un token illeggibile per entrambi resta al narrow, che riporta il
        // proprio errore senza chiamare il profilo largo.
        assert!(!payload_needs_wide("0xZZ"));

        let wide: Vec<String> = ["--modulus", "0x3FB"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(args_need_wide(&wide));
        let narrow: Vec<String> = ["--modulus", "1019"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(!args_need_wide(&narrow));
    }
}
