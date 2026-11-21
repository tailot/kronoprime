//! Factorization pipeline, from cheapest to strongest stage:
//!
//! 1. Small-prime extraction and trial division up to 1009 (a hard cap, so
//!    numbers with only large factors are not wasted on trial division).
//! 2. Deterministic Miller-Rabin screening of every cofactor (see
//!    [`crate::primality`]).
//! 3. Pollard's rho with Brent's cycle detection and batched `gcd`
//!    (product of 128 differences before each `gcd` call), with several
//!    fallback `(c, y0)` parameter pairs for degenerate cycles.
//! 4. Iterative decomposition through an explicit worklist stack, so
//!    pathological inputs never overflow the call stack.
//!
//! All modular arithmetic is overflow safe for every modulus up to
//! `u128::MAX`. Note that the 3x3 grid index `N(n)` only locates the decade
//! of `n`; it gives no asymptotic shortcut for factorization, which depends
//! on the smallest prime factor `p_min` and is handled by the stages above.

use crate::primality::{add_mod, is_prime, mul_mod};

/// Build the sorted list of primes up to 1009 (exactly 169 of them) at
/// compile time.
const fn build_small_primes() -> [u32; 169] {
    let mut out = [0u32; 169];
    let mut count: usize = 0;
    let mut n: u32 = 2;
    while n <= 1009 {
        let mut prime = true;
        let mut d = 2u32;
        while d * d <= n {
            if n.is_multiple_of(d) {
                prime = false;
                break;
            }
            d += 1;
        }
        if prime {
            out[count] = n;
            count += 1;
        }
        n += 1;
    }
    out
}

/// Primes used by the trial-division stage (up to 1009 inclusive).
const SMALL_PRIMES: [u32; 169] = build_small_primes();

// Compile-time check that the 169-slot table is exactly full.
const _: () = assert!(build_small_primes()[168] == 1009);

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Pollard-Brent for a single `(c, y0)` parameter pair on a composite odd
/// modulus `n` with no prime factors below 1011.
///
/// Returns a non-trivial factor, or `None` for degenerate runs (cycle
/// without a usable difference, or the iteration cap), in which case the
/// caller retries with other parameters.
fn brent_rho(n: u128, c: u128, y0: u128) -> Option<u128> {
    debug_assert!(c < n && y0 < n);
    let f = |x: u128| add_mod(mul_mod(x, x, n), c, n);
    let m: u32 = 128; // batch size for the accumulated product of differences
    let mut y = y0;
    let mut ys = y;
    let mut r: u32 = 1;
    loop {
        let x = y;
        for _ in 0..r {
            y = f(y);
        }
        let mut k: u32 = 0;
        let mut g: u128 = 1;
        while k < r && g == 1 {
            ys = y;
            let upper = m.min(r - k);
            let mut q = 1u128;
            for _ in 0..upper {
                y = f(y);
                q = mul_mod(q, x.abs_diff(y), n);
            }
            g = gcd(q, n);
            k += m;
        }
        if g == 1 {
            if r > (1 << 22) {
                return None; // iteration cap: retry with new parameters
            }
            r *= 2;
            continue;
        }
        if g != n {
            return Some(g);
        }
        // g == n: the batch product is 0 modulo n; replay the saved batch
        // one step at a time to isolate a proper factor.
        for _ in 0..m {
            ys = f(ys);
            let g2 = gcd(x.abs_diff(ys), n);
            if g2 > 1 {
                return if g2 < n { Some(g2) } else { None };
            }
        }
        return None;
    }
}

/// Search a non-trivial factor of a composite odd `n` with no prime factor
/// below 1011, retrying across `(c, y0)` pairs and then with fresh `c`.
fn pollard_rho(n: u128) -> u128 {
    const PARAMS: [(u128, u128); 10] = [
        (1, 2),
        (3, 2),
        (5, 3),
        (7, 2),
        (11, 3),
        (13, 2),
        (17, 4),
        (19, 5),
        (23, 3),
        (29, 7),
    ];
    for &(c, y0) in &PARAMS {
        if c < n && y0 < n {
            if let Some(d) = brent_rho(n, c, y0) {
                return d;
            }
        }
    }
    let mut c = 31u128;
    loop {
        if let Some(d) = brent_rho(n, c, 2) {
            return d;
        }
        c += 2;
    }
}

/// Prime factorization of `n` in ascending order (empty for `n < 2`).
///
/// The decomposition is driven by an explicit stack (worklist) rather than
/// recursion, and the final factor list satisfies `product == n` with every
/// entry prime.
pub fn factor(n: u128) -> Vec<u128> {
    let mut factors = Vec::new();
    if n < 2 {
        return factors;
    }
    let mut stack: Vec<u128> = vec![n];
    while let Some(m) = stack.pop() {
        if m < 2 {
            continue;
        }
        if is_prime(m) {
            factors.push(m);
            continue;
        }
        // Trial-division stage: extract every prime <= 1009.
        let mut rem = m;
        for &p in &SMALL_PRIMES {
            let p = u128::from(p);
            if p * p > rem {
                break;
            }
            while rem.is_multiple_of(p) {
                factors.push(p);
                rem /= p;
            }
        }
        if rem == 1 {
            continue;
        }
        if is_prime(rem) {
            factors.push(rem);
            continue;
        }
        // Pollard-Brent stage for the remaining composite cofactor.
        let d = pollard_rho(rem);
        stack.push(d);
        stack.push(rem / d);
    }
    factors.sort_unstable();
    factors
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_primes_table_is_valid() {
        assert_eq!(SMALL_PRIMES[0], 2);
        assert_eq!(SMALL_PRIMES[168], 1009);
        assert!(SMALL_PRIMES.windows(2).all(|w| w[0] < w[1]));
        assert!(SMALL_PRIMES.iter().all(|&p| is_prime(u128::from(p))));
    }

    #[test]
    fn empty_below_two() {
        assert!(factor(0).is_empty());
        assert!(factor(1).is_empty());
    }

    #[test]
    fn small_values() {
        assert_eq!(factor(2), vec![2]);
        assert_eq!(factor(12), vec![2, 2, 3]);
        assert_eq!(factor(97), vec![97]);
        assert_eq!(factor(97 * 101), vec![97, 101]);
        assert_eq!(factor(2u128.pow(10)), vec![2; 10]);
    }

    #[test]
    fn prime_powers() {
        assert_eq!(factor(59u128.pow(4)), vec![59; 4]);
        assert_eq!(factor(1000003u128.pow(3)), vec![1000003; 3]);
        assert_eq!(factor(10007u128.pow(4)), vec![10007; 4]);
    }

    #[test]
    fn carmichael_factorizations() {
        assert_eq!(factor(561), vec![3, 11, 17]);
        assert_eq!(factor(1105), vec![5, 13, 17]);
        assert_eq!(factor(1729), vec![7, 13, 19]);
        assert_eq!(factor(2465), vec![5, 17, 29]);
    }

    #[test]
    fn u64_edge_cases() {
        assert_eq!(
            factor(u128::from(u64::MAX)),
            vec![3, 5, 17, 257, 641, 65537, 6700417]
        );
        assert_eq!(factor(18446744073709551557), vec![18446744073709551557]);
        assert_eq!(factor(1000003u128 * 1000003), vec![1000003, 1000003]);
    }

    #[test]
    fn u128_semiprime_and_triple() {
        assert_eq!(
            factor(10000000019u128 * 10000000033),
            vec![10000000019, 10000000033]
        );
        let triple = 1000000007u128 * 1000000009 * 1000000033;
        assert_eq!(factor(triple), vec![1000000007, 1000000009, 1000000033]);
    }

    #[test]
    fn product_invariant_and_all_factors_prime() {
        let samples = [
            2u128,
            41,
            561,
            2u128.pow(20) * 3 * 104729,
            u128::from(u64::MAX),
            u128::from(u64::MAX) - 1,
            99999999998 * 37 + 7,
            10000000019u128 * 10000000033,
            170141183460469231731687303715884105727, // 2^127 - 1
        ];
        for n in samples {
            let fs = factor(n);
            let product = fs.iter().try_fold(1u128, |acc, &f| acc.checked_mul(f));
            assert_eq!(product, Some(n), "product of factors must equal {n}");
            assert!(fs.iter().all(|&f| is_prime(f)), "all factors of {n} prime");
        }
    }
}
