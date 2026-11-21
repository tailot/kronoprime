//! Integration tests: end-to-end coherence between the matrix grid,
//! primality, factorization and the CLI report layer.

use kronoprime::cli::{self, status_of};
use kronoprime::factor::factor;
use kronoprime::matrix::{self, slot_of, slot_product};
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
