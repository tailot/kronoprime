//! Primality testing.
//!
//! Miller-Rabin runs deterministic for the whole standard `u64` profile
//! (`n < 2^64`) thanks to a known compact base set, with two shorter
//! base sets used automatically below `4.76e9` and `2.15e12` (still fully
//! deterministic there, roughly halving the work for the small candidates
//! produced by decade searches). For the extended `u128` profile a fixed
//! set of 40 prime bases is used: the first 13 of them already detect
//! every strong pseudoprime below `3.3e24`, and the full 40-base test has
//! a false-prime probability far below `1e-20` for larger values, which
//! is beyond the practical needs of this engine.
//!
//! The modular helpers here are shared with the factorization stage and
//! are overflow safe for every modulus up to `u128::MAX`.

/// Small primes screened before Miller-Rabin.
const TRIAL_PRIMES: [u128; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Deterministic Miller-Rabin bases for every `n < 2^64`.
const MR_BASES_U64: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];

/// Deterministic bases below [`BASES3_LIMIT`] (classic three-base result).
const MR_BASES_3: [u64; 3] = [2, 7, 61];

/// Deterministic bases below [`BASES5_LIMIT`] (first five prime bases).
const MR_BASES_5: [u64; 5] = [2, 3, 5, 7, 11];

/// Each short bound equals the smallest strong pseudoprime for its base
/// set, so the shorter tests are exactly deterministic strictly below it;
/// at and above the bound the full `u64` set takes over.
const BASES3_LIMIT: u128 = 4_759_123_141;
const BASES5_LIMIT: u128 = 2_152_302_898_747;

/// Fixed bases for the extended `u128` profile (the first 40 primes).
/// Deterministic below `3.3e24`; beyond that the strong probable prime
/// error is bounded by `4^-40` (~`1e-24`).
const MR_BASES_EXTENDED: [u64; 40] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
];

/// Overflow-safe modular addition `(a + b) mod n` for `a, b < n`.
pub(crate) fn add_mod(a: u128, b: u128, n: u128) -> u128 {
    // `n - b` is in `1..=n` because `b < n`; the taken branch never overflows.
    if a >= n - b {
        a - (n - b)
    } else {
        a + b
    }
}

/// Modular multiplication for moduli up to `u128::MAX`.
///
/// When both operands fit in 64 bits their product fits in `u128` and a
/// native remainder is used (the hot path for the whole `u64` profile).
/// Otherwise multiply-and-add over `add_mod` keeps the result overflow free
/// at the cost of at most 128 additions.
pub(crate) fn mul_mod(a: u128, b: u128, n: u128) -> u128 {
    debug_assert!(a < n && b < n);
    if a == 0 || b == 0 {
        return 0;
    }
    if a <= u128::from(u64::MAX) && b <= u128::from(u64::MAX) {
        return (a * b) % n;
    }
    let mut result = 0u128;
    let mut base = a;
    let mut exp = b;
    while exp > 0 {
        if exp & 1 == 1 {
            result = add_mod(result, base, n);
        }
        base = add_mod(base, base, n);
        exp >>= 1;
    }
    result
}

/// Modular exponentiation by square-and-multiply.
fn pow_mod(mut base: u128, mut exp: u128, n: u128) -> u128 {
    debug_assert!(base < n);
    let mut result = 1u128;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, n);
        }
        base = mul_mod(base, base, n);
        exp >>= 1;
    }
    result
}

/// One strong probable prime test for base `a`: returns false only when
/// `a` is a Miller-Rabin witness that `n` is composite.
fn miller_rabin_pass(n: u128, d: u128, s: u32, a: u128) -> bool {
    let a = a % n;
    if a < 2 {
        return true;
    }
    let mut x = pow_mod(a, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 1..s {
        x = mul_mod(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

/// True when `n` is prime.
///
/// Deterministic for every `n < 2^64` (the entire standard `u64` profile,
/// including `u64::MAX`, Mersenne primes and Carmichael numbers). For
/// `2^64 <= n <= u128::MAX` the fixed 40-base test of [`MR_BASES_EXTENDED`]
/// is used, with the error bound documented there.
pub fn is_prime(n: u128) -> bool {
    if n < 2 {
        return false;
    }
    for &p in &TRIAL_PRIMES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let n1 = n - 1;
    let s = n1.trailing_zeros();
    let d = n1 >> s;
    let bases: &[u64] = if n < BASES3_LIMIT {
        // Most decade-search candidates live here: three modular
        // exponentiations instead of seven, still fully deterministic.
        &MR_BASES_3
    } else if n < BASES5_LIMIT {
        &MR_BASES_5
    } else if n < (1u128 << 64) {
        &MR_BASES_U64
    } else {
        &MR_BASES_EXTENDED
    };
    bases
        .iter()
        .all(|&a| miller_rabin_pass(n, d, s, u128::from(a)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_values() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        for p in [2, 3, 5, 7, 11, 13, 97, 101, 104729] {
            assert!(is_prime(p), "{p} should be prime");
        }
        for c in [4, 6, 9, 121, 1000000] {
            assert!(!is_prime(c), "{c} should be composite");
        }
    }

    #[test]
    fn carmichael_numbers_are_detected_as_composite() {
        for c in [561u128, 1105, 1729, 2465, 2821, 8321, 3215031751] {
            assert!(!is_prime(c), "Carmichael/pseudoprime {c} must be composite");
        }
    }

    #[test]
    fn large_primes() {
        assert!(is_prime(18446744073709551557)); // largest prime < u64::MAX region
        assert!(is_prime(2305843009213693951)); // 2^61 - 1 (Mersenne)
        assert!(is_prime(170141183460469231731687303715884105727)); // 2^127 - 1
    }

    #[test]
    fn large_composites() {
        assert!(!is_prime(u128::from(u64::MAX)));
        assert!(!is_prime(18446744073709551614));
        assert!(!is_prime(10000000019u128 * 10000000033));
        assert!(!is_prime(1000003u128 * 1000003));
    }

    #[test]
    fn short_base_boundaries_stay_deterministic() {
        // Smallest strong pseudoprimes for the short base sets: both sit
        // exactly at the bounds, where the full seven-base u64 set runs.
        assert!(!is_prime(4_759_123_141)); // first SPSP for bases 2, 7, 61
        assert!(!is_prime(2_152_302_898_747)); // first SPSP for bases 2..11
                                               // Primes just below the bounds still pass the short tests (the
                                               // largest primes strictly under each bound).
        assert!(is_prime(4_759_123_129));
        assert!(is_prime(2_152_302_898_729));
    }

    #[test]
    fn modular_helpers_are_consistent() {
        // add_mod/mul_mod agreement with checked arithmetic on random-ish values.
        let n = 12345678901234567890u128;
        for k in 1u128..50 {
            let a = (k * 987654321) % n;
            let b = (k * 123456789012) % n;
            assert_eq!(add_mod(a, b, n), (a + b) % n);
            assert_eq!(mul_mod(a, b, n), (a * b) % n);
        }
        // mul_mod doubling path: operands above 2^64.
        let big = u128::MAX - 192;
        assert_eq!(mul_mod(big - 1, big - 2, big), 2);
        assert_eq!(pow_mod(3, 1, big), 3);
    }
}
