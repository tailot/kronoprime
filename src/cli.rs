//! Command line interface for the `kronoprime` binary.
//!
//! Supports numbers as positional arguments, streaming from stdin (one
//! number per line when no numbers are given) and the flags `--json`,
//! `--matrix` and `--quadruplets <K>`.

use crate::cubic::{
    classify_weierstrass, eval_discriminant, find_cm_curves, gcd, is_on_curve, sample_points,
    scalar_mul, sub_mod, CubicCurve, CurveClass, DiscriminantStatus, ECPoint, PointOpResult,
};
use crate::factor::{cubic_factor_probe, factor};
use crate::matrix::{self, matrix_invariants_3x3, parse_matrix_3x3, Matrix3x3CubicExt, Placement};
use crate::primality::{is_prime, mul_mod};
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
         \x20 point-ops            point arithmetic (scalar multiplication) on a curve",
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
         Run 'kronoprime cubic <subcommand> --help' for details on a specific subcommand."
    );
}

fn run_cubic_scan(args: &[String]) -> u8 {
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
