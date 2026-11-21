//! Integration tests: end-to-end coherence between the matrix grid,
//! primality, factorization and the CLI report layer.

use kronoprime::cli::{self, status_of};
use kronoprime::cubic::{classify_weierstrass, inv_mod, CurveClass};
use kronoprime::factor::factor;
use kronoprime::matrix::{
    self, compute_similar_matrix, extract_depressed_params, slot_of, slot_product,
};
use kronoprime::primality::is_prime;
use std::time::Instant;

/// Every analysed composite must satisfy `product(factors) == n` with all
/// factors prime, and the decade placement must be coherent with the value.
#[test]
fn factorization_and_placement_invariants() {
    let samples: [u128; 10] = [
        2,
        41,
        561,
        1729,
        2u128.pow(40) * 104729,
        u128::from(u64::MAX),
        18446744073709551557,
        10000000019 * 10000000033,
        1000003u128.pow(5),
        170141183460469231731687303715884105727,
    ];
    for &n in &samples {
        let analysis = cli::analyze(n);
        assert_eq!(analysis.n, n);
        assert_eq!(analysis.prime, is_prime(n));
        let product = analysis
            .factors
            .iter()
            .try_fold(1u128, |acc, &f| acc.checked_mul(f));
        assert_eq!(product, Some(n), "factor product for {n}");
        assert!(analysis.factors.iter().all(|&f| is_prime(f)));
        if let Some(p) = analysis.placement {
            assert_eq!(
                matrix::cell_value(p.index, p.row, p.col),
                n,
                "grid cell must hold the value for {n}"
            );
            assert_eq!(p.slot, slot_of(n));
        } else {
            assert!(n == 0 || n % 10 == 0, "only borders lack placement: {n}");
        }
    }
}

/// For values coprime with 10 the slot of the product equals the product of
/// the factor slots modulo 10.
#[test]
fn slot_compatibility_over_factorizations() {
    let samples: [u128; 6] = [
        143,
        221,
        11 * 13 * 17,
        37 * 41,
        10000000019 * 10000000033,
        1000003u128.pow(3),
    ];
    for &n in &samples {
        let combined = factor(n)
            .iter()
            .fold(1u8, |acc, &f| slot_product(acc, slot_of(f)));
        assert_eq!(combined, slot_of(n), "slot product rule for {n}");
    }
}

/// The flagship u128 semiprime must factor well within the time budget.
#[test]
fn fast_large_semiprime() {
    let start = Instant::now();
    let factors = factor(100000000520000000627);
    let elapsed = start.elapsed();
    assert_eq!(factors, vec![10000000019, 10000000033]);
    assert!(elapsed.as_millis() < 5_000, "took {elapsed:?}");
}

/// CLI reports stay coherent for mixed inputs.
#[test]
fn cli_reports_end_to_end() {
    let semiprime = cli::analyze(100000000520000000627);
    assert_eq!(status_of(&semiprime), "COMPOSITE (large prime factors)");
    let text = cli::text_report(&semiprime, true);
    assert!(text.contains("3x3 grid for Matrix N=10000000052000000063"));
    assert!(text.contains("Matrix N: 10000000052000000063 | Slot: +7"));

    let mersenne = cli::analyze(170141183460469231731687303715884105727);
    assert_eq!(status_of(&mersenne), "PRIME");

    let json = cli::json_report(&mersenne);
    assert!(json.contains("\"prime\":true"));
    assert!(json.contains("\"factors\":[\"170141183460469231731687303715884105727\"]"));
}

/// Test 1: Degenerazione Nodale e Cuspidale su Primi Noti (Spec Section 5.1).
#[test]
fn test_singular_curve_classification() {
    let p = 101u128;

    // Caso Cuspidale: a = 0, b = 0 mod p
    let curve_cusp = classify_weierstrass(0, 0, p).unwrap();
    assert_eq!(curve_cusp.discriminant, 0);
    assert!(matches!(
        curve_cusp.class,
        CurveClass::CuspidalSingular { .. }
    ));

    // Caso Nodale: 4*a^3 + 27*b^2 = 0 mod p (es. a = 100 (-1), b = 0 non nodale; a = 100, b calcolato)
    // Scegliamo a = 100 (-1 mod 101), 4*(-1) + 27*b^2 = 0 => 27*b^2 = 4 => b^2 = 4 * 27^-1 mod 101
    // 27^-1 mod 101 = 26. b^2 = 4 * 26 = 104 = 3 mod 101. b = 45 (infatti 45^2 = 2025 = 3 mod 101)
    let curve_nodal = classify_weierstrass(100, 45, p).unwrap();
    assert_eq!(curve_nodal.discriminant, 0);
    assert!(matches!(
        curve_nodal.class,
        CurveClass::NodalSingular { .. }
    ));

    // Verifica analitica nodale supplementare con radici esatte modulo 101: x0 = 1 => a = 98, b = 2
    let curve_nodal_exact = classify_weierstrass(98, 2, p).unwrap();
    assert_eq!(curve_nodal_exact.discriminant, 0);
    assert!(matches!(
        curve_nodal_exact.class,
        CurveClass::NodalSingular { .. }
    ));
}

/// Test 2: Estrazione del Fattore Non Banale su Semiprimi (Spec Section 5.2).
#[test]
#[allow(clippy::erasing_op, clippy::identity_op)]
fn test_factor_extraction_via_discriminant_gcd() {
    let p = 1009u128;
    let q = 2003u128;
    let n = p * q; // Semiprimo

    // Costruiamo parametri tali che 4a^3 + 27b^2 sia un multiplo di p ma non di q
    // Modulo p: a_p = 0, b_p = 0 (singolare)
    // Modulo q: a_q = 1, b_q = 1 (regolare)
    // Ricostruzione CRT per a e b:
    let a = (0 * q * inv_mod(q, p) + 1 * p * inv_mod(p, q)) % n;
    let b = (0 * q * inv_mod(q, p) + 1 * p * inv_mod(p, q)) % n;

    let res = classify_weierstrass(a, b, n);
    assert!(res.is_err()); // Deve fallire restituendo l'intercettazione del fattore
}

/// Test 3: Certificazione Ordine per Curve CM Supersingolari (Spec Section 5.3).
#[test]
fn test_cm_supersingular_closed_form_order() {
    let _p = 1009u128; // 1009 % 4 == 1, usiamo p = 1019 dove 1019 % 3 == 2 e 1019 % 4 == 3
    let p_cm = 1019u128;

    // Tipo I: y^2 = x^3 + 5 (a = 0, p = 2 mod 3)
    let curve_type1 = classify_weierstrass(0, 5, p_cm).unwrap();
    if let CurveClass::CmSupersingular { closed_form_order } = curve_type1.class {
        assert_eq!(closed_form_order, p_cm + 1);
    } else {
        panic!("La curva doveva essere classificata come CM Supersingolare Tipo I");
    }

    // Tipo II: y^2 = x^3 + 7x (b = 0, p = 3 mod 4)
    let curve_type2 = classify_weierstrass(7, 0, p_cm).unwrap();
    if let CurveClass::CmSupersingular { closed_form_order } = curve_type2.class {
        assert_eq!(closed_form_order, p_cm + 1);
    } else {
        panic!("La curva doveva essere classificata come CM Supersingolare Tipo II");
    }
}

/// Test 4: Coerenza della Trasformazione da Matrice a Cubica Ridotta (Spec Section 5.4).
#[test]
fn test_matrix_to_depressed_cubic_invariance() {
    let n = 1009u128;
    // Matrice M generica
    let m = [
        [12u128, 5u128, 7u128],
        [3u128, 19u128, 2u128],
        [8u128, 1u128, 30u128],
    ];

    let (a1, b1) = extract_depressed_params(&m, n);

    // Matrice P * M * P^-1 ha lo stesso polinomio caratteristico
    // Verifichiamo che i parametri della forma depressa rimangano invariati
    let m_similar = compute_similar_matrix(&m, n);
    let (a2, b2) = extract_depressed_params(&m_similar, n);

    assert_eq!(a1, a2);
    assert_eq!(b1, b2);
}

/// Test end-to-end della CLI per i sottocomandi `kronoprime cubic`.
#[test]
fn test_cli_cubic_subcommands_end_to_end() {
    // 1. Scan subcommands
    let args_scan_table: Vec<String> = [
        "cubic",
        "scan",
        "--modulus",
        "1019",
        "--limit",
        "3",
        "--class",
        "cm",
        "--format",
        "table",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_scan_table), 0);

    let args_scan_json: Vec<String> = [
        "cubic",
        "scan",
        "--modulus",
        "1019",
        "--limit",
        "1",
        "--class",
        "cm",
        "--format",
        "json",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_scan_json), 0);

    // 2. From-matrix subcommands
    let args_from_matrix: Vec<String> = [
        "cubic",
        "from-matrix",
        "--matrix",
        "[12,5,7; 3,19,2; 8,1,30]",
        "--modulus",
        "1009",
        "--format",
        "json",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_from_matrix), 0);

    // 3. Verify subcommands
    let args_verify: Vec<String> = [
        "cubic",
        "verify",
        "--a",
        "0",
        "--b",
        "5",
        "--modulus",
        "1019",
        "--points",
        "3",
        "--format",
        "table",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_verify), 0);

    // 4. Factor-probe subcommands (exit code 2 on finding factor)
    let args_factor: Vec<String> = [
        "cubic",
        "factor-probe",
        "--modulus",
        "2021027",
        "--limit",
        "50",
        "--format",
        "json",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_factor), 2);

    // 5. Point-ops subcommands
    let args_point_ops: Vec<String> = [
        "cubic",
        "point-ops",
        "--a",
        "0",
        "--b",
        "5",
        "--modulus",
        "1019",
        "--point",
        "12,431",
        "--scalar",
        "1020",
        "--format",
        "table",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(cli::run(&args_point_ops), 0);
}
