//! 3x3 decade grid model: matrix index `N`, slot arithmetic and the mod-3003
//! wheel used to prune candidate decades in quadruplet searches.
//!
//! For a matrix index `N >= 1`, with base `B = 10 * (N - 1)`, the decade
//! spans `[B + 1, B + 9]` laid out in a 3x3 grid: value at row `r` and
//! column `c` is `B + 3r + c + 1`. Active slots (possible primes > 5) are
//! the final digits `+1`, `+3`, `+7`, `+9`.

use crate::cubic::{classify_weierstrass, sub_mod, try_inv_mod, CubicCurve};
use crate::cubic4096::{classify_curve_4096, Curve4096, CurveError4096};
use crate::primality::{add_mod, is_prime, mul_mod};
use crate::u4096::U4096;

/// Final digits of the four active slots of every decade.
pub const ACTIVE_SLOTS: [u8; 4] = [1, 3, 7, 9];

/// Modulus of the quadruplet wheel: `3 * 7 * 11 * 13 = 3003`.
pub const WHEEL_MODULUS: u128 = 3003;

/// Primes captured by the wheel.
pub const WHEEL_PRIMES: [u128; 4] = [3, 7, 11, 13];

/// Position of a value inside the 3x3 decade grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    /// Matrix (decade) index `N = floor((x - 1) / 10) + 1`.
    pub index: u128,
    /// Slot `s(x) = ((x - 1) mod 10) + 1`, in `1..=9`.
    pub slot: u8,
    /// Grid row in `0..=2`.
    pub row: u8,
    /// Grid column in `0..=2`.
    pub col: u8,
}

/// Slot (position inside the decade) of a positive value, per the direct
/// formula `((x - 1) mod 10) + 1`.
///
/// Note: multiples of 10 sit on the border between two decades and are not
/// part of any 3x3 grid; for them the formula yields the boundary value 10.
pub fn slot_of(x: u128) -> u8 {
    ((x - 1) % 10 + 1) as u8
}

/// Matrix index `N(x) = floor((x - 1) / 10) + 1`.
pub fn matrix_index(x: u128) -> u128 {
    (x - 1) / 10 + 1
}

/// Grid row and column for a slot in `1..=9`.
pub fn row_col(slot: u8) -> (u8, u8) {
    let i = slot - 1;
    (i / 3, i % 3)
}

/// Full placement of `x`.
///
/// Returns `None` for `x = 0` and for multiples of 10, which lie between
/// two consecutive decades and have no grid cell.
pub fn try_placement(x: u128) -> Option<Placement> {
    if x < 1 || x.is_multiple_of(10) {
        return None;
    }
    let slot = slot_of(x);
    let (row, col) = row_col(slot);
    Some(Placement {
        index: matrix_index(x),
        slot,
        row,
        col,
    })
}

/// Value stored in cell `(row, col)` of matrix `index`.
pub fn cell_value(index: u128, row: u8, col: u8) -> u128 {
    10 * (index - 1) + 3 * u128::from(row) + u128::from(col) + 1
}

/// The nine values of decade `index` as a 3x3 grid.
pub fn decade_cells(index: u128) -> [[u128; 3]; 3] {
    let base = 10 * (index - 1) + 1;
    let mut cells = [[0u128; 3]; 3];
    for (r, row) in cells.iter_mut().enumerate() {
        for (c, cell) in row.iter_mut().enumerate() {
            *cell = base + 3 * r as u128 + c as u128;
        }
    }
    cells
}

/// Final slot digit produced by multiplying values with slots `a` and `b`.
///
/// For primes `p, q > 5` the product slot is `s_p * s_q (mod 10)`: slot 3
/// comes from the pairs `(1, 3)` and `(7, 9)`, slot 7 from `(1, 7)` and
/// `(3, 9)`; slots 1 and 9 accept every combination of `{1, 3, 7, 9}`.
pub fn slot_product(a: u8, b: u8) -> u8 {
    (a * b) % 10
}

/// True when a value with slot `factor_slot` can divide a number whose slot
/// is `n_slot`, i.e. when the cofactor can carry a compatible active slot.
///
/// This filter is diagnostic and educational: on a single prime factor the
/// modulo-10 constraint alone never excludes a final digit, because every
/// slot in `{1, 3, 7, 9}` is invertible modulo 10. The real pruning power of
/// the grid lives in the modular wheel below, not in slot compatibility.
pub fn matches_slot(n_slot: u8, factor_slot: u8) -> bool {
    ACTIVE_SLOTS
        .iter()
        .any(|&s| slot_product(factor_slot, s) == n_slot)
}

/// True when `M = N - 1` passes the mod-3003 wheel, i.e. no slot value
/// `10M + o` (with `o` in `{1, 3, 7, 9}`) is a proper multiple of 3, 7, 11
/// or 13. A slot value equal to the wheel prime itself stays admissible
/// (that is how decade N=2, holding 11/13, survives the sieve).
// `From` is not a const trait yet, so a widening cast is used inside.
#[allow(clippy::cast_lossless)]
pub const fn wheel_admissible(m: u128) -> bool {
    let mut i = 0;
    while i < WHEEL_PRIMES.len() {
        let p = WHEEL_PRIMES[i];
        let mut j = 0;
        while j < ACTIVE_SLOTS.len() {
            let o = ACTIVE_SLOTS[j] as u128;
            let divisible = (10 * (m % p) + o).is_multiple_of(p);
            // Exception: the divisible value may be the prime p itself.
            // 10m + o == p forces m <= 1, so this stays overflow-free.
            if divisible && !(m <= 1 && 10 * m + o == p) {
                return false;
            }
            j += 1;
        }
        i += 1;
    }
    true
}

/// Number of admissible residue classes modulo 3003 (190 of 3003, ~6.3%).
const ADMISSIBLE_COUNT: usize = 190;

const fn count_admissible() -> usize {
    let mut count = 0usize;
    let mut m = 0u128;
    while m < WHEEL_MODULUS {
        if wheel_admissible(m) {
            count += 1;
        }
        m += 1;
    }
    count
}

const fn build_admissible() -> [u32; ADMISSIBLE_COUNT] {
    let mut out = [0u32; ADMISSIBLE_COUNT];
    let mut count = 0usize;
    let mut m = 0u128;
    while m < WHEEL_MODULUS {
        if wheel_admissible(m) {
            out[count] = m as u32;
            count += 1;
        }
        m += 1;
    }
    out
}

// The hand-set count must match the const sieve exactly.
const _: () = assert!(count_admissible() == ADMISSIBLE_COUNT);

/// Sorted admissible residues of `M = N - 1` modulo 3003, precomputed at
/// compile time (190 classes out of the 3003-class wheel period).
pub const ADMISSIBLE_RESIDUES: [u32; ADMISSIBLE_COUNT] = build_admissible();

/// Admissible residues of `M = N - 1` modulo 3003 (190 classes).
pub fn admissible_residues() -> Vec<u128> {
    ADMISSIBLE_RESIDUES.iter().map(|&r| u128::from(r)).collect()
}

/// First `count` matrix indices `N` (ascending) whose four active slots
/// `10(N-1)+1`, `+3`, `+7`, `+9` are all prime.
///
/// Instead of testing every decade against the wheel, the search jumps
/// directly between admissible residue classes of `N - 1` modulo 3003
/// (190 visits per 3003 decades), so ~93.7% of decades never even reach
/// the filter. Only surviving candidates run primality tests, and the
/// four slot tests short-circuit at the first composite slot.
pub fn find_quadruplets(count: usize) -> Vec<u128> {
    let mut out = Vec::with_capacity(count);
    let mut block = 0u128;
    'blocks: while out.len() < count {
        for &r in &ADMISSIBLE_RESIDUES {
            let m = block + u128::from(r);
            let base = 10 * m;
            if ACTIVE_SLOTS.iter().all(|&o| is_prime(base + u128::from(o))) {
                out.push(m + 1);
                if out.len() == count {
                    break 'blocks;
                }
            }
        }
        block += WHEEL_MODULUS;
    }
    out
}

/// Computes the invariants `(trace, sigma2, determinant)` of a 3x3 matrix modulo `n`.
pub fn matrix_invariants_3x3(m: &[[u128; 3]; 3], n: u128) -> (u128, u128, u128) {
    let tr = add_mod(add_mod(m[0][0] % n, m[1][1] % n, n), m[2][2] % n, n);

    let m00 = sub_mod(
        mul_mod(m[1][1], m[2][2], n),
        mul_mod(m[1][2], m[2][1], n),
        n,
    );
    let m11 = sub_mod(
        mul_mod(m[0][0], m[2][2], n),
        mul_mod(m[0][2], m[2][0], n),
        n,
    );
    let m22 = sub_mod(
        mul_mod(m[0][0], m[1][1], n),
        mul_mod(m[0][1], m[1][0], n),
        n,
    );
    let sigma2 = add_mod(add_mod(m00, m11, n), m22, n);

    let c01 = sub_mod(
        mul_mod(m[1][0], m[2][2], n),
        mul_mod(m[1][2], m[2][0], n),
        n,
    );
    let c02 = sub_mod(
        mul_mod(m[1][0], m[2][1], n),
        mul_mod(m[1][1], m[2][0], n),
        n,
    );

    let term0 = mul_mod(m[0][0], m00, n);
    let term1 = mul_mod(m[0][1], c01, n);
    let term2 = mul_mod(m[0][2], c02, n);
    let det = add_mod(sub_mod(term0, term1, n), term2, n);

    (tr, sigma2, det)
}

/// Multiply two 3x3 matrices modulo `n`.
pub fn mat_mul_3x3(a: &[[u128; 3]; 3], b: &[[u128; 3]; 3], n: u128) -> [[u128; 3]; 3] {
    let mut out = [[0u128; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            let mut sum = 0u128;
            for k in 0..3 {
                sum = add_mod(sum, mul_mod(a[r][k] % n, b[k][c] % n, n), n);
            }
            out[r][c] = sum;
        }
    }
    out
}

/// Generates a similar matrix `M' = P * M * P^-1 mod n` with det(P) = 1.
pub fn compute_similar_matrix(m: &[[u128; 3]; 3], n: u128) -> [[u128; 3]; 3] {
    // Invertible matrix P = [[2, 1, 0], [1, 1, 0], [0, 0, 1]], det(P) = 1
    // P^-1 = [[1, -1, 0], [-1, 2, 0], [0, 0, 1]]
    let p = [[2 % n, 1 % n, 0], [1 % n, 1 % n, 0], [0, 0, 1 % n]];
    let p_inv = [
        [1 % n, (n - 1) % n, 0],
        [(n - 1) % n, 2 % n, 0],
        [0, 0, 1 % n],
    ];
    let pm = mat_mul_3x3(&p, m, n);
    mat_mul_3x3(&pm, &p_inv, n)
}

/// Extracts depressed Weierstrass parameters `(a, b)` directly from a 3x3 matrix.
pub fn extract_depressed_params(m: &[[u128; 3]; 3], modulus: u128) -> (u128, u128) {
    m.to_depressed_cubic(modulus)
        .expect("failed to extract depressed params")
}

/// Parses a 3x3 matrix from a string in formats such as `"[m11,m12,m13; m21,m22,m23; m31,m32,m33]"`
/// or `"[[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]"`.
pub fn parse_matrix_3x3(s: &str) -> Result<[[u128; 3]; 3], &'static str> {
    let mut numbers = Vec::new();
    let mut current_num = String::new();

    for ch in s.chars() {
        if ch.is_ascii_digit() {
            current_num.push(ch);
        } else if !current_num.is_empty() {
            if let Ok(v) = current_num.parse::<u128>() {
                numbers.push(v);
            } else {
                return Err("Invalid matrix number");
            }
            current_num.clear();
        }
    }
    if !current_num.is_empty() {
        if let Ok(v) = current_num.parse::<u128>() {
            numbers.push(v);
        } else {
            return Err("Invalid matrix number");
        }
    }

    if numbers.len() != 9 {
        return Err("Invalid matrix syntax: expected 9 entries for a 3x3 matrix");
    }

    Ok([
        [numbers[0], numbers[1], numbers[2]],
        [numbers[3], numbers[4], numbers[5]],
        [numbers[6], numbers[7], numbers[8]],
    ])
}

/// Extension trait connecting 3x3 matrices to cubic and elliptic curves.
pub trait Matrix3x3CubicExt {
    /// Maps invariants of a 3x3 matrix (trace, minors, det) into Weierstrass parameters `(a, b)`.
    fn to_depressed_cubic(&self, modulus: u128) -> Result<(u128, u128), &'static str>;

    /// Converts a 3x3 matrix directly into its associated `CubicCurve`.
    fn extract_cubic_curve(&self, modulus: u128) -> Result<CubicCurve, &'static str>;
}

impl Matrix3x3CubicExt for [[u128; 3]; 3] {
    fn to_depressed_cubic(&self, modulus: u128) -> Result<(u128, u128), &'static str> {
        let n = modulus;
        if n <= 1 {
            return Err("Invalid modulus");
        }

        let (tr, sigma2, det) = matrix_invariants_3x3(self, n);

        // Coefficients of char poly: P(t) = t^3 + c2*t^2 + c1*t + c0
        let c2 = sub_mod(0, tr, n);
        let c1 = sigma2 % n;
        let c0 = sub_mod(0, det, n);

        // Inverses of 3 and 27
        let inv3 = match try_inv_mod(3, n) {
            Ok(inv) => inv,
            Err(_) => return Err("Modulus not coprime with 3 (cannot invert 3)"),
        };
        let inv27 = mul_mod(mul_mod(inv3, inv3, n), inv3, n);

        let c2_sq = mul_mod(c2, c2, n);
        let c2_cu = mul_mod(c2_sq, c2, n);

        // a = c1 - c2^2 / 3 mod n
        let a = sub_mod(c1, mul_mod(c2_sq, inv3, n), n);

        // b = c0 - c1*c2 / 3 + 2*c2^3 / 27 mod n
        let term_mid = mul_mod(mul_mod(c1, c2, n), inv3, n);
        let term_cubic = mul_mod(mul_mod(2, c2_cu, n), inv27, n);
        let b = add_mod(sub_mod(c0, term_mid, n), term_cubic, n);

        Ok((a, b))
    }

    fn extract_cubic_curve(&self, modulus: u128) -> Result<CubicCurve, &'static str> {
        let (a, b) = self.to_depressed_cubic(modulus)?;
        classify_weierstrass(a, b, modulus)
    }
}

// ---------------------------------------------------------------------------
// The 4096-bit profile
//
// The decade machinery above is tied to decimal digits and therefore to the
// `u128` window it can print. The invariants of a matrix are not: a trace, a sum
// of principal minors and a determinant are just three residue classes, and the
// characteristic polynomial they define can be depressed into Weierstrass
// parameters `(a, b)` at any width. Everything below works on `U4096` cells with
// a `U4096` modulus, keeps every product inside the wide buffer of `U4096`
// itself, and never converts through a machine word.
// ---------------------------------------------------------------------------

/// A 3x3 matrix of `U4096` residues.
pub type Matrix4096 = [[U4096; 3]; 3];

/// What can go wrong while taking invariants or depressed parameters at width.
// The payloads are full-width words and the point of the wide profile is to
// stay stack-allocated, so the variants are kept unboxed on purpose.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MatrixError4096 {
    /// The modulus is `0` or `1`: there is no residue ring to reduce into.
    InvalidModulus,
    /// `3` is not invertible modulo the modulus; the stored value is the gcd,
    /// which for a composite modulus is a factor in disguise.
    ThreeNotInvertible(U4096),
    /// The cubic side rejected the curve; `EllipticCurveMethod`-style failures
    /// carry the factor they found.
    Curve(CurveError4096),
}

impl MatrixError4096 {
    /// The message the command line prints for this failure.
    pub fn message(&self) -> &'static str {
        match self {
            MatrixError4096::InvalidModulus => "Invalid modulus",
            MatrixError4096::ThreeNotInvertible(_) => {
                "Modulus not coprime with 3 (cannot invert 3)"
            }
            MatrixError4096::Curve(CurveError4096::InvalidModulus) => "Invalid modulus",
            MatrixError4096::Curve(CurveError4096::FactorFound { .. }) => {
                "Factor found during curve classification"
            }
        }
    }
}

/// Turns a `u128` matrix into wide cells, without touching their values.
pub fn widen_matrix_3x3(m: &[[u128; 3]; 3]) -> Matrix4096 {
    let mut out = [[U4096::ZERO; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            out[r][c] = U4096::from_u128(m[r][c]);
        }
    }
    out
}

/// Computes `(trace, sigma2, determinant)` of a 3x3 matrix modulo `n`, all wide.
///
/// `sigma2` is the sum of the three principal `2x2` minors, so that
/// `t^3 - trace * t^2 + sigma2 * t - det` is the characteristic polynomial.
pub fn matrix_invariants_4096(m: &Matrix4096, n: &U4096) -> (U4096, U4096, U4096) {
    if n <= &U4096::ONE {
        return (U4096::ZERO, U4096::ZERO, U4096::ZERO);
    }
    let cell = |r: usize, c: usize| m[r][c].rem(n);
    let minor = |a: usize, b: usize, c: usize, d: usize| {
        cell(a, b)
            .mul_mod(&cell(c, d), n)
            .sub_mod(&cell(a, d).mul_mod(&cell(c, b), n), n)
    };

    let tr = cell(0, 0).add_mod(&cell(1, 1), n).add_mod(&cell(2, 2), n);

    let m00 = minor(1, 1, 2, 2);
    let m11 = minor(0, 0, 2, 2);
    let m22 = minor(0, 0, 1, 1);
    let sigma2 = m00.add_mod(&m11, n).add_mod(&m22, n);

    let c01 = minor(1, 0, 2, 2);
    let c02 = minor(1, 0, 2, 1);

    let det = cell(0, 0)
        .mul_mod(&m00, n)
        .sub_mod(&cell(0, 1).mul_mod(&c01, n), n)
        .add_mod(&cell(0, 2).mul_mod(&c02, n), n);

    (tr, sigma2, det)
}

/// Multiply two wide 3x3 matrices modulo `n`.
pub fn mat_mul_4096(a: &Matrix4096, b: &Matrix4096, n: &U4096) -> Matrix4096 {
    let mut out = [[U4096::ZERO; 3]; 3];
    if n <= &U4096::ONE {
        return out;
    }
    for r in 0..3 {
        for c in 0..3 {
            let mut sum = U4096::ZERO;
            for k in 0..3 {
                sum = sum.add_mod(&a[r][k].rem(n).mul_mod(&b[k][c].rem(n), n), n);
            }
            out[r][c] = sum;
        }
    }
    out
}

/// The wide analogue of [`compute_similar_matrix`]: `M' = P M P^-1` with the
/// same unimodular `P`, so the invariants of `M` and `M'` agree by construction.
pub fn compute_similar_matrix_4096(m: &Matrix4096, n: &U4096) -> Matrix4096 {
    if n <= &U4096::ONE {
        return [[U4096::ZERO; 3]; 3];
    }
    let one = U4096::ONE.rem(n);
    let two = U4096::TWO.rem(n);
    let minus_one = U4096::ONE.neg_mod(n);
    let zero = U4096::ZERO;
    let p = [[two, one, zero], [one, one, zero], [zero, zero, one]];
    let p_inv = [
        [one, minus_one, zero],
        [minus_one, two, zero],
        [zero, zero, one],
    ];
    mat_mul_4096(&mat_mul_4096(&p, m, n), &p_inv, n)
}

/// Depresses the characteristic polynomial of a wide matrix into the
/// Weierstrass parameters `(a, b)` of `y^2 = x^3 + a x + b mod n`.
#[allow(clippy::result_large_err)]
pub fn depressed_params_4096(
    m: &Matrix4096,
    modulus: &U4096,
) -> Result<(U4096, U4096), MatrixError4096> {
    let n = modulus;
    if n <= &U4096::ONE {
        return Err(MatrixError4096::InvalidModulus);
    }
    let (tr, sigma2, det) = matrix_invariants_4096(m, n);

    // Characteristic polynomial t^3 + c2 t^2 + c1 t + c0.
    let c2 = tr.neg_mod(n);
    let c1 = sigma2;
    let c0 = det.neg_mod(n);

    let inv3 = U4096::THREE
        .try_inv_mod(n)
        .map_err(MatrixError4096::ThreeNotInvertible)?;
    let inv27 = inv3.mul_mod(&inv3, n).mul_mod(&inv3, n);

    let c2_sq = c2.mul_mod(&c2, n);
    let c2_cu = c2_sq.mul_mod(&c2, n);

    // a = c1 - c2^2 / 3
    let a = c1.sub_mod(&c2_sq.mul_mod(&inv3, n), n);
    // b = c0 - c1 c2 / 3 + 2 c2^3 / 27
    let term_mid = c1.mul_mod(&c2, n).mul_mod(&inv3, n);
    let term_cubic = c2_cu.mul_mod(&U4096::TWO.rem(n), n).mul_mod(&inv27, n);
    let b = c0.sub_mod(&term_mid, n).add_mod(&term_cubic, n);

    Ok((a, b))
}

/// The wide analogue of [`Matrix3x3CubicExt::extract_cubic_curve`]: the matrix
/// becomes a classified curve, and a composite modulus can hand back a factor.
#[allow(clippy::result_large_err)]
pub fn extract_cubic_curve_4096(
    m: &Matrix4096,
    modulus: &U4096,
) -> Result<Curve4096, MatrixError4096> {
    let (a, b) = depressed_params_4096(m, modulus)?;
    classify_curve_4096(&a, &b, modulus).map_err(MatrixError4096::Curve)
}

/// Parses a wide 3x3 matrix in either of the two bracket styles accepted by
/// [`parse_matrix_3x3`], with decimal or `0x` cells up to 4096 bits.
pub fn parse_matrix_4096(s: &str) -> Result<Matrix4096, &'static str> {
    let mut numbers: Vec<U4096> = Vec::new();
    let mut token = String::new();
    let flush = |token: &mut String, numbers: &mut Vec<U4096>| -> Result<(), &'static str> {
        if token.is_empty() {
            return Ok(());
        }
        match U4096::from_str(token) {
            Ok(v) => {
                numbers.push(v);
                token.clear();
                Ok(())
            }
            Err(_) => Err("Invalid matrix number"),
        }
    };
    for ch in s.chars() {
        if ch.is_ascii_hexdigit() || ch == 'x' || ch == 'X' {
            token.push(ch);
        } else {
            flush(&mut token, &mut numbers)?;
        }
    }
    flush(&mut token, &mut numbers)?;

    if numbers.len() != 9 {
        return Err("Invalid matrix syntax: expected 9 entries for a 3x3 matrix");
    }
    let mut out = [[U4096::ZERO; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        row.copy_from_slice(&numbers[i * 3..i * 3 + 3]);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_and_index_roundtrip() {
        for x in 1u128..=1000 {
            let p = try_placement(x);
            if x % 10 == 0 {
                assert_eq!(p, None, "multiple of 10 has no cell: {x}");
                assert_eq!(slot_of(x), 10);
                continue;
            }
            let p = p.expect("placement for x >= 1, x % 10 != 0");
            assert_eq!(p.slot, slot_of(x));
            assert_eq!(p.index, matrix_index(x));
            assert_eq!(cell_value(p.index, p.row, p.col), x);
        }
    }

    #[test]
    fn known_placements() {
        assert_eq!(
            try_placement(100000000520000000627),
            Some(Placement {
                index: 10000000052000000063,
                slot: 7,
                row: 2,
                col: 0,
            })
        );
        assert_eq!(try_placement(11).unwrap().index, 2);
        assert_eq!(try_placement(11).unwrap().slot, 1);
        assert_eq!(
            try_placement(19),
            Some(Placement {
                index: 2,
                slot: 9,
                row: 2,
                col: 2,
            })
        );
        assert_eq!(try_placement(0), None);
        assert_eq!(try_placement(20), None);
    }

    #[test]
    fn decade_cells_match_grid_formula() {
        let cells = decade_cells(2);
        assert_eq!(cells, [[11, 12, 13], [14, 15, 16], [17, 18, 19]]);
        let cells = decade_cells(10000000052000000063);
        assert_eq!(cells[2][0], 100000000520000000627);
    }

    #[test]
    fn wheel_keeps_exactly_190_residues() {
        let residues = admissible_residues();
        assert_eq!(residues.len(), 190);
        assert!(wheel_admissible(1)); // decade 11..19 survives (11 and 13
                                      // are wheel primes themselves)
        assert!(!wheel_admissible(0)); // 3 divides 10*0 + 3
        assert!(!wheel_admissible(2)); // 3 divides 10*2 + 1 and 10*2 + 7
    }

    #[test]
    fn admissible_table_is_sorted_and_complete() {
        assert_eq!(ADMISSIBLE_RESIDUES.len(), 190);
        assert!(ADMISSIBLE_RESIDUES.windows(2).all(|w| w[0] < w[1]));
        for &r in &ADMISSIBLE_RESIDUES {
            assert!(wheel_admissible(u128::from(r)), "residue {r}");
        }
        assert!(ADMISSIBLE_RESIDUES.contains(&1)); // decade N=2 survives
        assert_eq!(admissible_residues()[0], 1);
    }

    #[test]
    fn first_quadruplets() {
        assert_eq!(find_quadruplets(5), vec![2, 11, 20, 83, 149]);
    }

    #[test]
    fn quadruplets_have_four_prime_slots() {
        for &n in &find_quadruplets(30) {
            assert!(wheel_admissible(n - 1));
            let base = 10 * (n - 1);
            for &s in &ACTIVE_SLOTS {
                assert!(is_prime(base + u128::from(s)), "N={n} slot +{s}");
            }
        }
    }

    #[test]
    fn slot_product_rules() {
        assert_eq!(slot_product(3, 7), 1);
        assert_eq!(slot_product(7, 9), 3);
        assert_eq!(slot_product(1, 7), 7);
        assert_eq!(slot_product(3, 9), 7);
        assert!(matches_slot(3, 1));
        assert!(matches_slot(7, 3));
        // 101 * 103 = 10403: slots 1 and 3 -> product slot 3.
        let a = try_placement(101).unwrap().slot;
        let b = try_placement(103).unwrap().slot;
        assert_eq!(slot_product(a, b), try_placement(101 * 103).unwrap().slot);
    }

    #[test]
    fn matrix_parsing_formats() {
        let m1 = parse_matrix_3x3("[12,5,7; 3,19,2; 8,1,30]").unwrap();
        let m2 = parse_matrix_3x3("[[12,5,7],[3,19,2],[8,1,30]]").unwrap();
        assert_eq!(m1, m2);
        assert_eq!(m1[0], [12, 5, 7]);
        assert_eq!(m1[1], [3, 19, 2]);
        assert_eq!(m1[2], [8, 1, 30]);
    }

    #[test]
    fn matrix_depressed_cubic_and_similarity() {
        let n = 1009u128;
        let m = [
            [12u128, 5u128, 7u128],
            [3u128, 19u128, 2u128],
            [8u128, 1u128, 30u128],
        ];
        let (a1, b1) = extract_depressed_params(&m, n);
        let m_similar = compute_similar_matrix(&m, n);
        let (a2, b2) = extract_depressed_params(&m_similar, n);
        assert_eq!(a1, a2);
        assert_eq!(b1, b2);

        let curve = m.extract_cubic_curve(n).unwrap();
        assert_eq!(curve.a, a1);
        assert_eq!(curve.b, b1);
        assert_eq!(curve.modulus, n);
    }

    fn scale_4096(m: &Matrix4096, s: &U4096, n: &U4096) -> Matrix4096 {
        let mut out = [[U4096::ZERO; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                out[r][c] = m[r][c].mul_mod(s, n);
            }
        }
        out
    }

    fn add_4096(a: &Matrix4096, b: &Matrix4096, n: &U4096) -> Matrix4096 {
        let mut out = [[U4096::ZERO; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                out[r][c] = a[r][c].add_mod(&b[r][c], n);
            }
        }
        out
    }

    #[test]
    fn wide_matrix_invariants_match_the_narrow_path() {
        let narrow = [
            [12u128, 5u128, 7u128],
            [3u128, 19u128, 2u128],
            [8u128, 1u128, 30u128],
        ];
        let n128 = 1009u128;
        let n = U4096::from_u128(n128);
        let m = widen_matrix_3x3(&narrow);

        let narrow_inv = matrix_invariants_3x3(&narrow, n128);
        let wide_inv = matrix_invariants_4096(&m, &n);
        assert_eq!(wide_inv.0.to_u128(), Some(narrow_inv.0));
        assert_eq!(wide_inv.1.to_u128(), Some(narrow_inv.1));
        assert_eq!(wide_inv.2.to_u128(), Some(narrow_inv.2));

        let narrow_params = extract_depressed_params(&narrow, n128);
        let (a, b) = depressed_params_4096(&m, &n).expect("1009 is coprime with 3");
        assert_eq!(a.to_u128(), Some(narrow_params.0));
        assert_eq!(b.to_u128(), Some(narrow_params.1));

        let wide_similar = compute_similar_matrix_4096(&m, &n);
        let narrow_similar = widen_matrix_3x3(&compute_similar_matrix(&narrow, n128));
        for r in 0..3 {
            for c in 0..3 {
                assert_eq!(wide_similar[r][c], narrow_similar[r][c].rem(&n));
            }
        }
    }

    #[test]
    fn wide_matrix_cayley_hamilton_and_similarity() {
        // 10^39 + 7: beyond every u128 cell, odd, and coprime with 3, so the
        // whole path runs in wide arithmetic with no narrow escape hatch.
        let n = U4096::from_dec("1000000000000000000000000000000000000007").unwrap();
        assert!(n > U4096::from_u128(u128::MAX));
        assert!(n.is_odd());
        assert_ne!(n.div_rem_u64(3).1, 0, "3 must stay invertible");
        let m: Matrix4096 = [
            [
                U4096::from_dec("999999999999999999999999999999999999999").unwrap(),
                U4096::from_u64(7),
                U4096::from_u128(u128::MAX),
            ],
            [
                U4096::from_u64(12345678901234567890),
                U4096::from_dec("500000000000000000000000000000000000001").unwrap(),
                U4096::from_u64(2),
            ],
            [
                U4096::from_u64(999999999999999999),
                U4096::from_u64(3),
                U4096::from_dec("123456789012345678901234567890123456789").unwrap(),
            ],
        ];

        let (tr, sigma2, det) = matrix_invariants_4096(&m, &n);
        // Cayley-Hamilton: M^3 - trace M^2 + sigma2 M - det I = 0 mod n.
        let m2 = mat_mul_4096(&m, &m, &n);
        let m3 = mat_mul_4096(&m2, &m, &n);
        let identity = [
            [U4096::ONE.rem(&n), U4096::ZERO, U4096::ZERO],
            [U4096::ZERO, U4096::ONE.rem(&n), U4096::ZERO],
            [U4096::ZERO, U4096::ZERO, U4096::ONE.rem(&n)],
        ];
        let mut check = m3;
        check = add_4096(&check, &scale_4096(&m2, &tr.neg_mod(&n), &n), &n);
        check = add_4096(&check, &scale_4096(&m, &sigma2, &n), &n);
        check = add_4096(&check, &scale_4096(&identity, &det.neg_mod(&n), &n), &n);
        for row in check {
            for cell in row {
                assert!(cell.is_zero(), "Cayley-Hamilton broke on {cell}");
            }
        }

        // Similar matrices share their parameters, and the curve agrees.
        let (a, b) = depressed_params_4096(&m, &n).expect("invertible 3");
        let similar = compute_similar_matrix_4096(&m, &n);
        let (a2, b2) = depressed_params_4096(&similar, &n).expect("invertible 3");
        assert_eq!((a, b), (a2, b2));
        let curve = extract_cubic_curve_4096(&m, &n);
        match curve {
            Ok(curve) => {
                assert_eq!(curve.a, a);
                assert_eq!(curve.b, b);
                assert_eq!(curve.modulus, n);
            }
            Err(MatrixError4096::Curve(CurveError4096::FactorFound { factor, cofactor })) => {
                // A composite modulus can surface a factor on the way to the
                // curve, and that factor really divides the wide modulus.
                assert!(factor > U4096::ONE && factor < n);
                assert!(cofactor > U4096::ONE && cofactor < n);
                assert_eq!(n.rem(&factor), U4096::ZERO);
            }
            Err(other) => panic!("unexpected refusal: {other:?}"),
        }
    }

    #[test]
    fn wide_matrix_parsing_and_refusals() {
        let big = "3402823669209384634633746074317682114583";
        let a = parse_matrix_4096(&format!("[1,2,3; 4,5,6; 7,8,{big}]")).unwrap();
        let b = parse_matrix_4096(&format!("[[1,2,3],[4,5,6],[7,8,{big}]]")).unwrap();
        assert_eq!(a, b);
        assert_eq!(a[2][2], U4096::from_dec(big).unwrap());
        assert_eq!(
            parse_matrix_4096("[1,2,3; 4,5,6; 7,8]").err(),
            Some("Invalid matrix syntax: expected 9 entries for a 3x3 matrix")
        );
        // A cell wider than the type is refused rather than silently truncated.
        assert!(U4096::from_hex(&"f".repeat(1024)).is_ok());
        assert!(U4096::from_dec(&"9".repeat(1240)).is_err());
        assert!(parse_matrix_4096(&format!("[{}{}]", "9".repeat(1240), "1")).is_err());

        assert_eq!(
            depressed_params_4096(&a, &U4096::ONE).err(),
            Some(MatrixError4096::InvalidModulus)
        );
        assert_eq!(
            depressed_params_4096(&a, &U4096::ZERO).err(),
            Some(MatrixError4096::InvalidModulus)
        );
        // 3057 = 3 * 1019: 3 is not invertible, and the gcd names the factor.
        let err = depressed_params_4096(&a, &U4096::from_u64(3057))
            .expect_err("3 is not invertible mod 3057");
        assert_eq!(err, MatrixError4096::ThreeNotInvertible(U4096::from_u64(3)));
        assert_eq!(
            err.message(),
            "Modulus not coprime with 3 (cannot invert 3)"
        );
    }
}
