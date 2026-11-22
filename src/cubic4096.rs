//! Cubic and elliptic curve analysis on 4096-bit moduli.
//!
//! This is the wide sibling of [`crate::cubic`]: the same discriminant test,
//! the same classification ladder and the same affine addition law, lifted to
//! [`U4096`] so that a modulus far beyond `u128` can still be classified. Two
//! things change once the modulus leaves the machine word:
//!
//! * Every modular step goes through [`U4096`]'s Knuth/Montgomery cores, so
//!   the cost model is quadratic in the modulus width rather than a single
//!   machine `mulmod`. Point arithmetic therefore keeps the affine shape of the
//!   small module for the *inspectable* operations - a failed inversion is a
//!   factor, which is the whole point of the probe in [`crate::factor`] - and
//!   switches to Jacobian coordinates inside scalar multiplication, where an
//!   inversion per step would dominate the run time.
//! * Point enumeration cannot walk the whole field, so sampling walks `x = 0,
//!   1, 2, ...` and stops as soon as enough residues are found. Roughly half of
//!   the `x` values are residues, so the expected cost tracks the number of
//!   points asked for, not the size of the field.
//!
//! Curves are given in short Weierstrass form `y^2 = x^3 + a x + b mod n`,
//! which needs `2` to be invertible for the group law to behave; callers only
//! build these curves over odd moduli.

use crate::primality::is_prime_u4096_with_bases;
use crate::u4096::{legendre_symbol_4096, mod_sqrt_4096, U4096};

/// How many Miller-Rabin bases a sampling pass spends certifying the modulus.
///
/// Sampling only means anything over a prime field, but a full 40-base
/// certification at 4096 bits costs about a second, which is unreasonable for
/// a helper whose job is to produce example points. Eight bases bound the
/// chance of walking over a pseudoprime by `4^-8`.
const SAMPLING_CERTIFICATES: usize = 8;

/// Result of evaluating the discriminant modulo `n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiscriminantStatus4096 {
    /// Regular curve: `gcd(Delta, n) == 1`.
    Regular(U4096),
    /// The discriminant vanishes modulo all of `n`.
    FullyDegenerate,
    /// The discriminant shares a proper factor with `n`.
    FactorFound(U4096),
}

/// Why a curve could not be built over `Z/nZ`.
// The factor payload is two full-width words, and this crate keeps the number
// path free of heap allocations, so the variants stay unboxed.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CurveError4096 {
    /// `n < 2`: there is no ring to define the curve over.
    InvalidModulus,
    /// `gcd(Delta, n)` exposed a proper factor of `n`.
    FactorFound { factor: U4096, cofactor: U4096 },
}

/// Topological and algebraic classes of cubic and elliptic curves.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CurveClass4096 {
    /// Nodal singular cubic (`E_reg` is a torus over `Z/nZ`).
    NodalSingular { singular_point: (U4096, U4096) },
    /// Cuspidal singular cubic (`E_reg` is the additive group).
    CuspidalSingular { singular_point: (U4096, U4096) },
    /// CM curve with a closed-form point count (`#E = n + 1`).
    CmSupersingular { closed_form_order: U4096 },
    /// Montgomery fast model `y^2 = x^3 + A x^2 + x`.
    MontgomeryForm { a_param: U4096, b_param: U4096 },
    /// Ordinary non-singular elliptic curve.
    StandardElliptic,
}

impl CurveClass4096 {
    pub fn name(&self) -> &'static str {
        match self {
            CurveClass4096::NodalSingular { .. } => "NodalSingular",
            CurveClass4096::CuspidalSingular { .. } => "CuspidalSingular",
            CurveClass4096::CmSupersingular { .. } => "CmSupersingular",
            CurveClass4096::MontgomeryForm { .. } => "MontgomeryForm",
            CurveClass4096::StandardElliptic => "StandardElliptic",
        }
    }
}

/// A cubic `y^2 = x^3 + a x + b mod n` together with its classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Curve4096 {
    pub a: U4096,
    pub b: U4096,
    pub modulus: U4096,
    pub discriminant: U4096,
    pub class: CurveClass4096,
}

/// A point on a short Weierstrass curve, or the identity.
// Copy-by-value of the whole 1024-byte point is deliberate: no heap, no borrow
// bookkeeping around the arithmetic, and the arrays stay on the stack.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Point4096 {
    Infinity,
    Affine(U4096, U4096),
}

/// Outcome of a point operation: a new point, or a factor of the modulus.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointOpResult4096 {
    Success(Point4096),
    FactorFound(U4096),
}

/// A small integer constant taken into the residue ring before it multiplies
/// anything. The arithmetic would be exact either way, but `mul_mod` and
/// `add_mod` keep the standing contract that both operands are already reduced,
/// and that matters as soon as the modulus is smaller than the constant.
fn const_mod(k: u64, n: &U4096) -> U4096 {
    U4096::from_u64(k).rem(n)
}

/// Evaluates `Delta = -16 * (4 a^3 + 27 b^2) mod n` and its gcd with `n`.
pub fn eval_discriminant_4096(a: &U4096, b: &U4096, n: &U4096) -> (U4096, DiscriminantStatus4096) {
    if n <= &U4096::ONE {
        return (U4096::ZERO, DiscriminantStatus4096::FullyDegenerate);
    }
    let a = a.rem(n);
    let b = b.rem(n);
    let four = const_mod(4, n);
    let sixteen = const_mod(16, n);
    let twenty_seven = const_mod(27, n);
    let a_cubed = a.mul_mod(&a, n).mul_mod(&a, n);
    let b_squared = b.mul_mod(&b, n);
    let term1 = four.mul_mod(&a_cubed, n);
    let term2 = twenty_seven.mul_mod(&b_squared, n);
    let inner = term1.add_mod(&term2, n);
    let delta = sixteen.mul_mod(&inner, n).neg_mod(n);

    let g = delta.gcd(n);
    let status = if g == *n {
        DiscriminantStatus4096::FullyDegenerate
    } else if g > U4096::ONE {
        DiscriminantStatus4096::FactorFound(g)
    } else {
        DiscriminantStatus4096::Regular(delta)
    };
    (delta, status)
}

/// Nodal tangent point: the double root `x0 = -3b / 2a` of the cubic.
///
/// Returns `(0, 0)` when `2a` is not invertible, the same fallback the
/// small-width module uses: the singularity is still certified by the
/// vanishing discriminant, only its coordinates are not pinned down.
fn nodal_point_4096(a: &U4096, b: &U4096, n: &U4096) -> (U4096, U4096) {
    if a.is_zero() {
        return (U4096::ZERO, U4096::ZERO);
    }
    let two_a = const_mod(2, n).mul_mod(a, n);
    match two_a.try_inv_mod(n) {
        Ok(inv) => {
            let num = const_mod(3, n).mul_mod(b, n).neg_mod(n);
            (num.mul_mod(&inv, n), U4096::ZERO)
        }
        Err(_) => (U4096::ZERO, U4096::ZERO),
    }
}

/// `n mod 4`, read straight off the two low bits.
fn mod_four(n: &U4096) -> U4096 {
    U4096::from_u64(u64::from(n.bit(0)) + 2 * u64::from(n.bit(1)))
}

/// Classifies `y^2 = x^3 + a x + b` modulo `n`.
// The error holds a full-width factor pair; boxing it would allocate on every
// rejection path, which the stack-only design of this module rules out.
#[allow(clippy::result_large_err)]
pub fn classify_curve_4096(a: &U4096, b: &U4096, n: &U4096) -> Result<Curve4096, CurveError4096> {
    if n <= &U4096::ONE {
        return Err(CurveError4096::InvalidModulus);
    }
    let a = a.rem(n);
    let b = b.rem(n);
    let (delta, status) = eval_discriminant_4096(&a, &b, n);

    match status {
        DiscriminantStatus4096::FactorFound(g) => Err(CurveError4096::FactorFound {
            factor: g.min(*n / g),
            cofactor: g.max(*n / g),
        }),
        DiscriminantStatus4096::FullyDegenerate => {
            let class = if a.is_zero() && b.is_zero() {
                CurveClass4096::CuspidalSingular {
                    singular_point: (U4096::ZERO, U4096::ZERO),
                }
            } else {
                CurveClass4096::NodalSingular {
                    singular_point: nodal_point_4096(&a, &b, n),
                }
            };
            Ok(Curve4096 {
                a,
                b,
                modulus: *n,
                discriminant: U4096::ZERO,
                class,
            })
        }
        DiscriminantStatus4096::Regular(_) => {
            // The closed-form point counts survive at any width: `y^2 = x^3 + b`
            // over a prime `n = 2 mod 3` and `y^2 = x^3 + a x` over a prime
            // `n = 3 mod 4` both have exactly `n + 1` points.
            let cm = (a.is_zero() && n.rem(&U4096::THREE) == U4096::TWO)
                || (b.is_zero() && mod_four(n) == U4096::THREE);
            let class = if cm {
                CurveClass4096::CmSupersingular {
                    closed_form_order: n.add(&U4096::ONE).0,
                }
            } else {
                CurveClass4096::StandardElliptic
            };
            Ok(Curve4096 {
                a,
                b,
                modulus: *n,
                discriminant: delta,
                class,
            })
        }
    }
}

/// Deterministic search for CM curves over a prime `p` (mirrors the small module).
pub fn find_cm_curves_4096(p: &U4096, limit: usize) -> Vec<Curve4096> {
    let mut results = Vec::new();
    if p <= &U4096::ONE {
        return results;
    }
    if p.rem(&U4096::THREE) == U4096::TWO {
        for b in 1..=(limit as u64) {
            if let Ok(curve) = classify_curve_4096(&U4096::ZERO, &U4096::from_u64(b), p) {
                results.push(curve);
            }
        }
    }
    if mod_four(p) == U4096::THREE {
        for a in 1..=(limit as u64) {
            if let Ok(curve) = classify_curve_4096(&U4096::from_u64(a), &U4096::ZERO, p) {
                results.push(curve);
            }
        }
    }
    results
}

/// Checks that `p` satisfies `y^2 = x^3 + a x + b mod n`.
pub fn is_on_curve_4096(p: &Point4096, a: &U4096, b: &U4096, n: &U4096) -> bool {
    match p {
        Point4096::Infinity => true,
        Point4096::Affine(x, y) => {
            if n <= &U4096::ONE {
                return false;
            }
            let x = x.rem(n);
            let y = y.rem(n);
            let a = a.rem(n);
            let b = b.rem(n);
            let lhs = y.mul_mod(&y, n);
            let x3 = x.mul_mod(&x, n).mul_mod(&x, n);
            let ax = a.mul_mod(&x, n);
            let rhs = x3.add_mod(&ax, n).add_mod(&b, n);
            lhs == rhs
        }
    }
}

/// Negates a point: `(x, y) -> (x, -y)`.
pub fn point_neg_4096(p: Point4096, n: &U4096) -> Point4096 {
    match p {
        Point4096::Infinity => Point4096::Infinity,
        Point4096::Affine(x, y) => Point4096::Affine(x.rem(n), y.rem(n).neg_mod(n)),
    }
}

/// Turns a failed inversion into the right report: a proper factor when there
/// is one, otherwise the identity (the denominator vanished modulo all of `n`).
fn inversion_outcome(g: U4096, n: &U4096) -> PointOpResult4096 {
    if g > U4096::ONE && g < *n {
        PointOpResult4096::FactorFound(g)
    } else {
        PointOpResult4096::Success(Point4096::Infinity)
    }
}

/// Affine addition on `y^2 = x^3 + a x + b mod n`.
///
/// A denominator that is not invertible is not an error: `gcd(denominator, n)`
/// is a proper factor whenever the curve was built over a composite modulus,
/// which is exactly the Lenstra mechanism the factor probe relies on.
pub fn point_add_4096(p1: Point4096, p2: Point4096, a: &U4096, n: &U4096) -> PointOpResult4096 {
    if n <= &U4096::ONE {
        return PointOpResult4096::Success(Point4096::Infinity);
    }
    let a = a.rem(n);
    let (x1, y1) = match p1 {
        Point4096::Infinity => return PointOpResult4096::Success(p2),
        Point4096::Affine(x, y) => (x.rem(n), y.rem(n)),
    };
    let (x2, y2) = match p2 {
        Point4096::Infinity => return PointOpResult4096::Success(p1),
        Point4096::Affine(x, y) => (x.rem(n), y.rem(n)),
    };

    let (num, denom) = if x1 == x2 {
        if y1.add_mod(&y2, n).is_zero() || y1 != y2 {
            return PointOpResult4096::Success(Point4096::Infinity);
        }
        // Doubling: lambda = (3 x1^2 + a) / (2 y1).
        let num = const_mod(3, n)
            .mul_mod(&x1, n)
            .mul_mod(&x1, n)
            .add_mod(&a, n);
        (num, const_mod(2, n).mul_mod(&y1, n))
    } else {
        (y2.sub_mod(&y1, n), x2.sub_mod(&x1, n))
    };

    let inv = match denom.try_inv_mod(n) {
        Ok(inv) => inv,
        Err(g) => return inversion_outcome(g, n),
    };
    let lambda = num.mul_mod(&inv, n);
    let x3 = lambda.mul_mod(&lambda, n).sub_mod(&x1, n).sub_mod(&x2, n);
    let y3 = lambda.mul_mod(&x1.sub_mod(&x3, n), n).sub_mod(&y1, n);
    PointOpResult4096::Success(Point4096::Affine(x3, y3))
}

/// Affine doubling on `y^2 = x^3 + a x + b mod n`.
pub fn point_double_4096(p: Point4096, a: &U4096, n: &U4096) -> PointOpResult4096 {
    point_add_4096(p, p, a, n)
}

// ---------------------------------------------------------------------------
// Jacobian coordinates for scalar multiplication.
//
// Affine arithmetic costs one inversion per step; at 4096 bits an inversion
// runs a binary extended gcd over 64 limbs, which is dearer than the whole
// rest of a point operation put together. Jacobian triples `(X : Y : Z)` with
// `x = X / Z^2` and `y = Y / Z^3` need no inversion at all, so a 4096-bit
// scalar multiplication pays a single inversion at the very end - and that one
// inversion is also where a composite modulus confesses a factor.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Jacobian4096 {
    x: U4096,
    y: U4096,
    z: U4096,
}

impl Jacobian4096 {
    const IDENTITY: Self = Self {
        x: U4096::ZERO,
        y: U4096::ZERO,
        z: U4096::ZERO,
    };

    fn is_infinity(&self) -> bool {
        self.z.is_zero()
    }

    fn from_affine(p: Point4096, n: &U4096) -> Self {
        match p {
            Point4096::Infinity => Self::IDENTITY,
            Point4096::Affine(x, y) => Self {
                x: x.rem(n),
                y: y.rem(n),
                z: U4096::ONE,
            },
        }
    }
}

/// `dbl-2009-l`: Jacobian doubling, valid for any short Weierstrass `a`.
fn jac_double_4096(p: Jacobian4096, a: &U4096, n: &U4096) -> Jacobian4096 {
    if p.is_infinity() || p.y.is_zero() {
        return Jacobian4096::IDENTITY;
    }
    let a = a.rem(n);
    let aa = p.x.mul_mod(&p.x, n);
    let bb = p.y.mul_mod(&p.y, n);
    let cc = bb.mul_mod(&bb, n);
    let sum = p.x.add_mod(&bb, n);
    // d = 2 * ((x + b)^2 - aa - cc)
    let dd = sum.mul_mod(&sum, n).sub_mod(&aa, n).sub_mod(&cc, n);
    let dd = dd.add_mod(&dd, n);
    // e = 3 aa + a z^4
    let z2 = p.z.mul_mod(&p.z, n);
    let z4 = z2.mul_mod(&z2, n);
    let e = const_mod(3, n)
        .mul_mod(&aa, n)
        .add_mod(&a.mul_mod(&z4, n), n);
    let f = e.mul_mod(&e, n);
    let x3 = f.sub_mod(&dd, n).sub_mod(&dd, n);
    let mut eight_cc = cc.add_mod(&cc, n);
    eight_cc = eight_cc.add_mod(&eight_cc, n);
    eight_cc = eight_cc.add_mod(&eight_cc, n);
    let y3 = e.mul_mod(&dd.sub_mod(&x3, n), n).sub_mod(&eight_cc, n);
    let z3 = p.y.mul_mod(&p.z, n);
    let z3 = z3.add_mod(&z3, n);
    Jacobian4096 {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// `madd-2007-bl`: generic Jacobian addition, falling back to doubling when the
/// two points coincide and to the identity when they are opposites.
fn jac_add_4096(p1: Jacobian4096, p2: Jacobian4096, a: &U4096, n: &U4096) -> Jacobian4096 {
    if p1.is_infinity() {
        return p2;
    }
    if p2.is_infinity() {
        return p1;
    }
    let z1z1 = p1.z.mul_mod(&p1.z, n);
    let z2z2 = p2.z.mul_mod(&p2.z, n);
    let u1 = p1.x.mul_mod(&z2z2, n);
    let u2 = p2.x.mul_mod(&z1z1, n);
    let s1 = p1.y.mul_mod(&p2.z, n).mul_mod(&z2z2, n);
    let s2 = p2.y.mul_mod(&p1.z, n).mul_mod(&z1z1, n);
    let h = u2.sub_mod(&u1, n);
    let r = s2.sub_mod(&s1, n);
    if h.is_zero() {
        if r.is_zero() {
            return jac_double_4096(p1, a, n);
        }
        return Jacobian4096::IDENTITY;
    }
    let hh = h.mul_mod(&h, n);
    let hhh = hh.mul_mod(&h, n);
    let u1hh = u1.mul_mod(&hh, n);
    let x3 = r
        .mul_mod(&r, n)
        .sub_mod(&hhh, n)
        .sub_mod(&u1hh, n)
        .sub_mod(&u1hh, n);
    let y3 = r
        .mul_mod(&u1hh.sub_mod(&x3, n), n)
        .sub_mod(&s1.mul_mod(&hhh, n), n);
    let z3 = h.mul_mod(&p1.z, n).mul_mod(&p2.z, n);
    Jacobian4096 {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// Converts back to affine, turning a non-invertible `Z` into a factor report.
fn jac_to_affine_4096(p: Jacobian4096, n: &U4096) -> PointOpResult4096 {
    if p.is_infinity() {
        return PointOpResult4096::Success(Point4096::Infinity);
    }
    match p.z.try_inv_mod(n) {
        Ok(inv_z) => {
            let w = inv_z.mul_mod(&inv_z, n);
            let x = p.x.mul_mod(&w, n);
            let y = p.y.mul_mod(&w, n).mul_mod(&inv_z, n);
            PointOpResult4096::Success(Point4096::Affine(x, y))
        }
        Err(g) => {
            if g > U4096::ONE && g < *n {
                PointOpResult4096::FactorFound(g)
            } else {
                // `Z` vanished modulo the whole modulus: the multiple is the
                // identity, the same report the affine code gives when a point
                // with `y = 0` doubles away.
                PointOpResult4096::Success(Point4096::Infinity)
            }
        }
    }
}

/// Scalar multiplication `[k] P` by double-and-add over the bits of `k`.
///
/// The walk is left to right, so the scalar never costs more than
/// `k.bits()` doublings, and the coordinates stay Jacobian until the last step,
/// where a single inversion either recovers an affine point or exposes a factor
/// of a composite modulus.
pub fn scalar_mul_4096(p: Point4096, k: &U4096, a: &U4096, n: &U4096) -> PointOpResult4096 {
    if n <= &U4096::ONE || k.is_zero() {
        return PointOpResult4096::Success(Point4096::Infinity);
    }
    let start = Jacobian4096::from_affine(p, n);
    if start.is_infinity() {
        return PointOpResult4096::Success(Point4096::Infinity);
    }
    let mut acc = start;
    for i in (0..(k.bits() - 1)).rev() {
        acc = jac_double_4096(acc, a, n);
        if k.bit(i) {
            acc = jac_add_4096(acc, start, a, n);
        }
        // No early exit here: an intermediate multiple can be the identity and
        // a later set bit still lifts it off it.
    }
    jac_to_affine_4096(acc, n)
}

/// Sample `count` points on `y^2 = x^3 + a x + b` over a prime modulus.
///
/// `x` walks upward from zero; about half of the candidates are quadratic
/// residues, so the loop terminates after roughly `2 * count` square-root
/// attempts regardless of how large the field is.
pub fn sample_points_4096(curve: &Curve4096, count: usize) -> Vec<(U4096, U4096)> {
    let mut points = Vec::new();
    if count == 0 || curve.modulus <= U4096::ONE {
        return points;
    }
    let p = curve.modulus;
    if !is_prime_u4096_with_bases(&p, SAMPLING_CERTIFICATES) {
        return points;
    }
    // A guard so a field with unusually few points cannot spin: `64 * count +
    // 64` tries is far above the expected cost of collecting `count` residues.
    let guard = 64 * count + 64;
    let mut x = U4096::ZERO;
    let mut tries = 0usize;
    while points.len() < count && tries < guard && x < p {
        tries += 1;
        let x3 = x.mul_mod(&x, &p).mul_mod(&x, &p);
        let ax = curve.a.mul_mod(&x, &p);
        let rhs = x3.add_mod(&ax, &p).add_mod(&curve.b, &p);
        if rhs.is_zero() {
            points.push((x, U4096::ZERO));
        } else if let Some(y) = mod_sqrt_4096(&rhs, &p) {
            points.push((x, y));
            if points.len() < count && !y.is_zero() {
                points.push((x, p.sub(&y).0));
            }
        }
        x = x.add(&U4096::ONE).0;
    }
    points.truncate(count);
    points
}

/// Sample `count` points at seeded abscissae instead of walking up from zero.
///
/// Deterministic but scattered, which is what a curve search wants: a fixed
/// seed replays the same walk, and neighbouring seeds explore different parts
/// of the field.
pub fn sample_points_seeded_4096(
    curve: &Curve4096,
    count: usize,
    seed: u64,
) -> Vec<(U4096, U4096)> {
    let mut points = Vec::new();
    if count == 0 || curve.modulus <= U4096::ONE {
        return points;
    }
    let p = curve.modulus;
    if !is_prime_u4096_with_bases(&p, SAMPLING_CERTIFICATES) {
        return points;
    }
    for i in 0..(8 * count + 8) {
        let x = U4096::sample_below(seed.wrapping_add(i as u64), &p);
        let x3 = x.mul_mod(&x, &p).mul_mod(&x, &p);
        let ax = curve.a.mul_mod(&x, &p);
        let rhs = x3.add_mod(&ax, &p).add_mod(&curve.b, &p);
        let y = if rhs.is_zero() {
            Some(U4096::ZERO)
        } else if legendre_symbol_4096(&rhs, &p) == 1 {
            mod_sqrt_4096(&rhs, &p)
        } else {
            None
        };
        if let Some(y) = y {
            points.push((x, y));
            if points.len() >= count {
                break;
            }
        }
    }
    points.truncate(count);
    points
}

// ---------------------------------------------------------------------------
// Montgomery ladder (x-coordinate only).
//
// The model is `y^2 = x^3 + A x^2 + x mod n`, whose points always come in pairs
// with the 2-torsion abscissa `x = 0`. That is what makes x-only arithmetic
// possible: for any two points,
//
//     x(Q + P) * x(Q - P) = ((x_Q x_P - 1) / (x_Q - x_P))^2
//
// so knowing `x_P` (the fixed difference) and the two consecutive multiples
// `[m] P` and `[m + 1] P` is enough to step the pair forward without ever
// touching a `y` or inverting inside the loop. Coordinates are projective
// pairs `(X : Z)` with `x = X / Z`; doubling is `dbl-2009-l` with
// `a24 = (A + 2) / 4`.
// ---------------------------------------------------------------------------

/// Doubling of a projective Montgomery pair.
fn ladder_double_4096(p: (U4096, U4096), a24: U4096, n: &U4096) -> (U4096, U4096) {
    let (x, z) = p;
    if z.is_zero() {
        return (U4096::ONE, U4096::ZERO);
    }
    let sum = x.add_mod(&z, n);
    let diff = x.sub_mod(&z, n);
    let aa = sum.mul_mod(&sum, n);
    let bb = diff.mul_mod(&diff, n);
    let e = aa.sub_mod(&bb, n);
    (
        aa.mul_mod(&bb, n),
        e.mul_mod(&bb.add_mod(&a24.mul_mod(&e, n), n), n),
    )
}

/// Differential addition `p + q`, where `p - q` is the known point `(xp : 1)`.
///
/// Writing `x = X / Z` in the identity `x(Q + P) * x(Q - P) = ratio^2` puts the
/// known abscissa in the *denominator*, because it is the abscissa of the
/// difference `Q - P`: the projective pair is
///
/// ```text
/// X_out = (X_p X_q - Z_p Z_q)^2 * Z_s
/// Z_out = (X_p Z_q - Z_p X_q)^2 * X_s
/// ```
///
/// for the difference point `(X_s : Z_s)`, and here that is `(xp : 1)`.
fn ladder_diff_add_4096(
    p: (U4096, U4096),
    q: (U4096, U4096),
    xp: U4096,
    n: &U4096,
) -> (U4096, U4096) {
    let (px, pz) = p;
    let (qx, qz) = q;
    // The difference formula divides by `x_p - x_q`, which says nothing about
    // the point at infinity, so the identity is handled by name.
    if pz.is_zero() {
        return (qx, qz);
    }
    if qz.is_zero() {
        return (px, pz);
    }
    let num = px.mul_mod(&qx, n).sub_mod(&pz.mul_mod(&qz, n), n);
    let den = px.mul_mod(&qz, n).sub_mod(&pz.mul_mod(&qx, n), n);
    (num.mul_mod(&num, n), xp.mul_mod(&den.mul_mod(&den, n), n))
}

/// Montgomery x-only ladder: the abscissa of `[k] P` from the abscissa of `P`.
///
/// `None` means the final denominator was not invertible - a factor of a
/// composite `n` in disguise - or that the multiple has no affine abscissa
/// (`k = 0`, or a point killed by the 2-torsion).
pub fn montgomery_ladder_x_4096(
    x_p: &U4096,
    k: &U4096,
    a_param: &U4096,
    n: &U4096,
) -> Option<U4096> {
    if k.is_zero() || n <= &U4096::ONE || n.is_even() {
        return None;
    }
    let xp = x_p.rem(n);
    let a = a_param.rem(n);
    let inv_four = U4096::FOUR.try_inv_mod(n).ok()?;
    let a24 = a.add_mod(&U4096::TWO, n).mul_mod(&inv_four, n);

    // The pair ([m] P, [m + 1] P); the projective identity is (1 : 0).
    let mut r0 = (U4096::ONE, U4096::ZERO);
    let mut r1 = (xp, U4096::ONE);

    for i in (0..k.bits()).rev() {
        if k.bit(i) {
            std::mem::swap(&mut r0, &mut r1);
        }
        r1 = ladder_diff_add_4096(r1, r0, xp, n);
        r0 = ladder_double_4096(r0, a24, n);
        if k.bit(i) {
            std::mem::swap(&mut r0, &mut r1);
        }
    }

    let (x, z) = r0;
    if z.is_zero() {
        return None;
    }
    Some(x.mul_mod(&z.try_inv_mod(n).ok()?, n))
}

/// Affine Montgomery duplication map `x(2P) = (x^2 - 1)^2 / (4x^3 + 4Ax^2 + 4x)`,
/// which follows from `lambda = (3x^2 + 2Ax + 1) / 2y` and the chord sum
/// `x_1 + x_2 + x_3 = -A` of the cubic `x^3 + Ax^2 + x`.
#[cfg(test)]
fn montgomery_double_x_4096(x: &U4096, a: &U4096, n: &U4096) -> Option<U4096> {
    let x = x.rem(n);
    let a = a.rem(n);
    let x2 = x.mul_mod(&x, n);
    let x3 = x2.mul_mod(&x, n);
    let num = x2.sub_mod(&U4096::ONE, n);
    let num = num.mul_mod(&num, n);
    let four = const_mod(4, n);
    let den = four
        .mul_mod(&x3, n)
        .add_mod(&four.mul_mod(&a, n).mul_mod(&x2, n), n)
        .add_mod(&four.mul_mod(&x, n), n);
    Some(num.mul_mod(&den.try_inv_mod(n).ok()?, n))
}

/// The difference identity evaluated in the affine plane, for test ground
/// truth. `x(Q + P) * x(Q - P)` is the squared ratio, so the known abscissa
/// `x(Q - P)` moves to the denominator.
#[cfg(test)]
fn montgomery_add_x_4096(xp: &U4096, xq: &U4096, x_diff: &U4096, n: &U4096) -> Option<U4096> {
    let xp = xp.rem(n);
    let xq = xq.rem(n);
    let num = xp.mul_mod(&xq, n).sub_mod(&U4096::ONE, n);
    let den = xp.sub_mod(&xq, n);
    let ratio = num.mul_mod(&den.try_inv_mod(n).ok()?, n);
    Some(
        ratio
            .mul_mod(&ratio, n)
            .mul_mod(&x_diff.rem(n).try_inv_mod(n).ok()?, n),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small(n: u64) -> U4096 {
        U4096::from_u64(n)
    }

    /// Tiny moduli put the reduced-operand contract under stress: the
    /// discriminant formula multiplies by 4, 16 and 27, every one of which
    /// dwarfs a modulus like 9, and a debug build turns an unreduced operand
    /// into a panic instead of a silent pass.
    #[test]
    fn tiny_moduli_keep_the_reduced_operand_contract() {
        for n in 2u64..=40 {
            let modulus = small(n);
            for a in 0..4u64 {
                for b in 0..4u64 {
                    let (delta, _) = eval_discriminant_4096(&small(a), &small(b), &modulus);
                    assert!(delta < modulus, "Delta is not reduced for n = {n}");
                    if let Ok(curve) = classify_curve_4096(&small(a), &small(b), &modulus) {
                        assert_eq!(curve.modulus, modulus);
                    }
                }
            }
            let point = Point4096::Affine(U4096::ONE, U4096::ONE);
            let _ = point_double_4096(point, &U4096::ONE, &modulus);
            let _ = point_add_4096(point, point, &U4096::ONE, &modulus);
            let _ = scalar_mul_4096(point, &small(6), &U4096::ONE, &modulus);
            let _ = find_cm_curves_4096(&modulus, 3);
        }
        // The ladder rejects even moduli, so 3 is the first width where the
        // constant 4 stops being a reduced residue on its own.
        for n in [3u64, 5, 9, 15] {
            let modulus = small(n);
            let _ = montgomery_ladder_x_4096(&U4096::ONE, &small(5), &U4096::ONE, &modulus);
            let _ = montgomery_double_x_4096(&U4096::ONE, &U4096::ONE, &modulus);
        }
    }

    #[test]
    fn wide_discriminant_matches_the_small_module() {
        // Delta(3, 4) = -16 * (4*27 + 27*16) = -16 * 540 = -8640, and 8640 mod
        // 1019 = 488, so Delta = 1019 - 488 = 531.
        let (delta, status) = eval_discriminant_4096(&small(3), &small(4), &small(1019));
        assert!(matches!(status, DiscriminantStatus4096::Regular(_)));
        assert_eq!(delta.to_u64(), Some(531));
        // Below two there is no ring, and the helper says so instead of
        // dividing by zero.
        let (d0, st0) = eval_discriminant_4096(&small(3), &small(4), &U4096::ONE);
        assert_eq!(d0, U4096::ZERO);
        assert_eq!(st0, DiscriminantStatus4096::FullyDegenerate);
    }

    #[test]
    fn wide_classification_and_cm_detection() {
        let p = small(1019); // 1019 = 2 mod 3 and 3 mod 4
        let curve = classify_curve_4096(&U4096::ZERO, &small(5), &p).unwrap();
        assert_eq!(curve.class.name(), "CmSupersingular");
        // Delta = -16 * 27 * 25 = -10800, and 10800 mod 1019 = 610.
        assert_eq!(curve.discriminant.to_u64(), Some(409));

        // A singular cubic: 4a^3 + 27b^2 = 0 mod 101 for a = -3, b = 2, and the
        // double root sits at x0 = -3b / 2a = -6 / -6 = 1.
        let nodal = classify_curve_4096(&small(98), &small(2), &small(101)).unwrap();
        assert_eq!(nodal.class.name(), "NodalSingular");
        assert_eq!(nodal.discriminant, U4096::ZERO);
        if let CurveClass4096::NodalSingular { singular_point } = nodal.class {
            assert_eq!(singular_point, (small(1), U4096::ZERO));
        } else {
            panic!("expected a nodal class");
        }

        // A cusp needs a = b = 0.
        let cuspidal = classify_curve_4096(&U4096::ZERO, &U4096::ZERO, &small(101)).unwrap();
        assert_eq!(cuspidal.class.name(), "CuspidalSingular");

        assert_eq!(
            classify_curve_4096(&small(1), &small(1), &U4096::ONE),
            Err(CurveError4096::InvalidModulus)
        );
    }

    #[test]
    fn wide_curve_rejects_composite_modulus_with_a_factor() {
        // Delta(0, 1) = -432 = -2^4 * 3^3 shares 9 with n = 9 * 1019.
        let n = small(9 * 1019);
        let err = classify_curve_4096(&U4096::ZERO, &small(1), &n).unwrap_err();
        if let CurveError4096::FactorFound { factor, cofactor } = err {
            assert_eq!(factor, small(9));
            assert_eq!(cofactor, small(1019));
            assert_eq!(factor * cofactor, n);
        } else {
            panic!("expected a factor report, got {err:?}");
        }
    }

    #[test]
    fn wide_point_arithmetic_matches_affine_expectations() {
        let p = small(1019);
        let a = U4096::ZERO;
        let b = small(5);
        // 987^2 = 974169 = 5 mod 1019 and 0^3 + 5 = 5, so (0, 987) is a point
        // on y^2 = x^3 + 5 over F_1019.
        let pt = Point4096::Affine(U4096::ZERO, small(987));
        assert!(is_on_curve_4096(&pt, &a, &b, &p));

        let doubled = point_double_4096(pt, &a, &p);
        let via_add = point_add_4096(pt, pt, &a, &p);
        assert_eq!(doubled, via_add);
        if let PointOpResult4096::Success(q) = doubled {
            assert!(is_on_curve_4096(&q, &a, &b, &p));
            let neg = point_neg_4096(q, &p);
            assert!(matches!(
                point_add_4096(q, neg, &a, &p),
                PointOpResult4096::Success(Point4096::Infinity)
            ));
            assert_eq!(
                point_add_4096(Point4096::Infinity, q, &a, &p),
                PointOpResult4096::Success(q)
            );
        } else {
            panic!("doubling must succeed over a prime field");
        }
        // This point has order three, so doubling lands on its own negative.
        if let PointOpResult4096::Success(q) = point_double_4096(pt, &a, &p) {
            assert_eq!(q, point_neg_4096(pt, &p));
        }
        // A point with y = 0 has order two: doubling reaches the identity.
        // x^3 + 5 = 0 mod 1019 has no solution here, so build the order-two
        // case on a curve that does have one: y^2 = x^3 + x over 5 has (0, 0).
        let two = Point4096::Affine(U4096::ZERO, U4096::ZERO);
        assert!(matches!(
            point_double_4096(two, &U4096::ONE, &small(5)),
            PointOpResult4096::Success(Point4096::Infinity)
        ));
        assert!(is_on_curve_4096(&two, &U4096::ONE, &U4096::ZERO, &small(5)));
    }

    #[test]
    fn jacobian_scalar_mul_agrees_with_affine() {
        let p = small(1019);
        let a = U4096::ZERO;
        let b = small(5);
        let curve = classify_curve_4096(&a, &b, &p).unwrap();
        // Several abscissae, so the walk covers both the small-order point at
        // x = 0 and generic points further along the field.
        for (x, y) in sample_points_4096(&curve, 6) {
            let pt = Point4096::Affine(x, y);
            for k in 1u64..25 {
                let fast = scalar_mul_4096(pt, &small(k), &a, &p);
                let mut slow = PointOpResult4096::Success(Point4096::Infinity);
                for _ in 0..k {
                    slow = match slow {
                        PointOpResult4096::Success(q) => point_add_4096(q, pt, &a, &p),
                        failure => failure,
                    };
                }
                assert_eq!(fast, slow, "[{k}]({x}, {y}) disagrees");
            }
        }
        assert!(matches!(
            scalar_mul_4096(
                Point4096::Affine(U4096::ZERO, small(987)),
                &small(3),
                &a,
                &p
            ),
            PointOpResult4096::Success(Point4096::Infinity)
        ));
        assert!(matches!(
            scalar_mul_4096(Point4096::Infinity, &small(7), &a, &p),
            PointOpResult4096::Success(Point4096::Infinity)
        ));
    }

    #[test]
    fn cm_closed_form_order_annihilates_points_at_width() {
        // y^2 = x^3 + 5 over the prime 1019 has 1020 points, so [1020] P = O
        // for every point; the Jacobian path reports the identity through a Z
        // that vanishes modulo the whole field.
        let p = small(1019);
        let curve = classify_curve_4096(&U4096::ZERO, &small(5), &p).unwrap();
        let order = match curve.class {
            CurveClass4096::CmSupersingular { closed_form_order } => closed_form_order,
            _ => panic!("expected the CM class"),
        };
        assert_eq!(order.to_u64(), Some(1020));
        let sampled = sample_points_4096(&curve, 3);
        assert_eq!(sampled.len(), 3);
        for (x, y) in sampled {
            let res = scalar_mul_4096(Point4096::Affine(x, y), &order, &curve.a, &curve.modulus);
            assert!(
                matches!(res, PointOpResult4096::Success(Point4096::Infinity)),
                "[{order}]({x}, {y}) did not reach O"
            );
        }
    }

    #[test]
    fn affine_inversion_failure_reports_a_factor() {
        // n = 1009 * 1019: a doubling whose denominator carries 1009, and an
        // addition whose denominator carries 1019. Both must come back as
        // proper divisors rather than as arithmetic errors.
        let n = small(1009 * 1019);
        let a = small(1);
        match point_double_4096(Point4096::Affine(small(7), small(2 * 1009)), &a, &n) {
            PointOpResult4096::FactorFound(g) => assert_eq!(g, small(1009)),
            other => panic!("expected the factor 1009, got {other:?}"),
        }
        match point_add_4096(
            Point4096::Affine(small(5), small(11)),
            Point4096::Affine(small(5 + 1019), small(29)),
            &a,
            &n,
        ) {
            PointOpResult4096::FactorFound(g) => assert_eq!(g, small(1019)),
            other => panic!("expected the factor 1019, got {other:?}"),
        }
    }

    #[test]
    fn wide_sampling_and_cm_search() {
        let p = small(1019);
        let curve = classify_curve_4096(&U4096::ZERO, &small(5), &p).unwrap();
        let points = sample_points_4096(&curve, 6);
        assert_eq!(points.len(), 6);
        for (x, y) in &points {
            assert!(is_on_curve_4096(
                &Point4096::Affine(*x, *y),
                &curve.a,
                &curve.b,
                &p
            ));
        }
        let seeded = sample_points_seeded_4096(&curve, 4, 0x5eed);
        assert_eq!(seeded.len(), 4);
        for (x, y) in &seeded {
            assert!(is_on_curve_4096(
                &Point4096::Affine(*x, *y),
                &curve.a,
                &curve.b,
                &p
            ));
        }
        // The same seed replays the same walk.
        assert_eq!(seeded, sample_points_seeded_4096(&curve, 4, 0x5eed));
        // 1019 is 2 mod 3 and 3 mod 4, so both CM families apply.
        assert_eq!(find_cm_curves_4096(&p, 3).len(), 6);
        assert!(sample_points_4096(&curve, 0).is_empty());
        assert!(sample_points_4096(
            &classify_curve_4096(&U4096::ZERO, &small(1), &small(4)).unwrap(),
            3
        )
        .is_empty());
    }

    /// One step of a full affine Montgomery walk on `y^2 = x^3 + A x^2 + x`,
    /// used as ground truth: `lambda = (3x^2 + 2Ax + 1) / 2y` when doubling,
    /// `lambda = (y_m - y_1) / (x_m - x_1)` afterwards, and in both cases
    /// `x_1 + x_2 + x_3 = lambda^2 - A` because the line cuts the cubic in a
    /// quadratic whose `x^2` coefficient is `A - lambda^2`.
    fn mont_step_full_4096(
        xm: (U4096, U4096),
        p1: (U4096, U4096),
        first: bool,
        a: &U4096,
        n: &U4096,
    ) -> Option<(U4096, U4096)> {
        let (xm_x, xm_y) = xm;
        let (x1, y1) = p1;
        let num = if first {
            xm_x.mul_mod(&xm_x, n)
                .mul_mod(&U4096::THREE, n)
                .add_mod(&U4096::TWO.mul_mod(a, n).mul_mod(&xm_x, n), n)
                .add_mod(&U4096::ONE, n)
        } else {
            xm_y.sub_mod(&y1, n)
        };
        let den = if first {
            U4096::TWO.mul_mod(&xm_y, n)
        } else {
            xm_x.sub_mod(&x1, n)
        };
        let lambda = num.mul_mod(&den.try_inv_mod(n).ok()?, n);
        let lam2 = lambda.mul_mod(&lambda, n);
        let out_x = lam2.sub_mod(a, n).sub_mod(&xm_x, n).sub_mod(&x1, n);
        let out_y = lambda
            .mul_mod(&xm_x.sub_mod(&out_x, n), n)
            .sub_mod(&xm_y, n);
        Some((out_x, out_y))
    }

    #[test]
    fn montgomery_ladder_follows_a_real_point_walk() {
        let n = small(1019);
        let a_param = small(4);
        // Take an abscissa that really carries a point over F_1019: the model is
        // `y^2 = x^3 + A x^2 + x`, and x-only arithmetic is only a group walk
        // when the right-hand side is a square.
        let mut start = None;
        for x in 2u64..40 {
            let xs = small(x);
            let x2 = xs.mul_mod(&xs, &n);
            let rhs = x2
                .mul_mod(&xs, &n)
                .add_mod(&a_param.mul_mod(&x2, &n), &n)
                .add_mod(&xs, &n);
            if let Some(y) = mod_sqrt_4096(&rhs, &n) {
                assert_eq!(y.mul_mod(&y, &n), rhs);
                start = Some((xs, y));
                break;
            }
        }
        let p1 = start.expect("a Montgomery abscissa over F_1019");

        // Ground truth: the walk taken with both coordinates, as far as the
        // denominators stay invertible.
        let mut walk = vec![p1];
        while walk.len() < 48 {
            let step = walk.len();
            let Some(next) = mont_step_full_4096(walk[step - 1], p1, step == 1, &a_param, &n)
            else {
                break;
            };
            let x = next.0;
            assert_eq!(
                x.mul_mod(&x, &n)
                    .mul_mod(&x, &n)
                    .add_mod(&a_param.mul_mod(&x, &n).mul_mod(&x, &n), &n)
                    .add_mod(&x, &n),
                next.1.mul_mod(&next.1, &n),
                "step {step} left the curve"
            );
            walk.push(next);
        }
        assert!(
            walk.len() >= 30,
            "the walk died too early to be a useful reference ({} steps)",
            walk.len()
        );

        // The ladder reproduces every abscissa of that walk, and so does the
        // affine recurrence `x_{m+1} = x_{m-1} * ((x_m x_1 - 1)/(x_m - x_1))^2`.
        for (k, point) in walk.iter().enumerate() {
            let scalar = small(k as u64 + 1);
            assert_eq!(
                montgomery_ladder_x_4096(&p1.0, &scalar, &a_param, &n),
                Some(point.0),
                "ladder disagrees with the walk at k = {}",
                k + 1
            );
        }
        for k in 2..walk.len() {
            let recurrence = montgomery_add_x_4096(&walk[k - 1].0, &p1.0, &walk[k - 2].0, &n)
                .expect("invertible denominator");
            assert_eq!(recurrence, walk[k].0, "recurrence breaks at k = {k}");
        }
        assert_eq!(
            montgomery_double_x_4096(&p1.0, &a_param, &n),
            Some(walk[1].0),
            "the duplication map is not the walk's doubling"
        );

        // Composition: [k2]([k1] P) reached by feeding the ladder an abscissa
        // already produced by the walk.
        for (k1, k2) in [(2u64, 10u64), (3, 7), (5, 5), (7, 3)] {
            let product = (k1 * k2) as usize;
            assert!(
                product <= walk.len(),
                "the walk was too short for {k1} * {k2}"
            );
            let inner = walk[(k1 as usize) - 1].0;
            let composed = montgomery_ladder_x_4096(&inner, &small(k2), &a_param, &n)
                .expect("invertible denominator");
            assert_eq!(
                composed,
                walk[product - 1].0,
                "[{k2}]([{k1}] P) differs at ({k1}, {k2})"
            );
        }

        assert_eq!(
            montgomery_ladder_x_4096(&p1.0, &U4096::ZERO, &a_param, &n),
            None
        );
        // An even modulus has no `4` to divide by, and no usable ladder.
        assert_eq!(
            montgomery_ladder_x_4096(&p1.0, &small(3), &a_param, &small(1020)),
            None
        );
    }

    /// Carico reale a piena larghezza: 4096 bit veri, non un campo piccolo.
    /// Misura `mul_mod`, una moltiplicazione scalare e una prova di primalità,
    /// e ancora la correttezza su un punto di ordine tre, il cui esito non
    /// dipende dal fatto che il modulo sia primo. Ignorato perché chiede
    /// secondi: `cargo test --release -- --ignored --nocapture`.
    #[test]
    #[ignore = "4096-bit load: a few seconds of execution"]
    fn full_width_4096_bit_workload_is_measured() {
        use crate::primality::is_prime_u4096_with_bases;
        use std::time::Instant;

        // Modulo dispare a piena larghezza: bit 4095 impostato e LSB dispari.
        let mut n = U4096::from_seed(0x243f_6a88_85a3_08d3);
        n.set_bit(4095, true);
        n.limbs[0] |= 1;
        assert_eq!(n.bits(), 4096);
        assert!(n.is_odd());

        // 512 moltiplicazioni modulari ridotte da un unico operando.
        let x = n.sub(&U4096::ONE).0;
        let t_mul = Instant::now();
        let mut acc = U4096::ONE;
        for _ in 0..512u32 {
            acc = acc.mul_mod(&x, &n);
        }
        let mul_nanos = t_mul.elapsed().as_nanos() as f64 / 512.0;

        // y^2 = x^3 + 4 ammette P = (0, 2) di ordine tre su ogni anello in cui
        // 2 è invertibile: 2P = (0, -2), 3P = O. Vale anche per un modulo
        // composto, quindi è un ancoraggio di correttezza indipendente.
        let a = U4096::ZERO;
        let b = U4096::from_u64(4);
        let p = Point4096::Affine(U4096::ZERO, U4096::from_u64(2));
        assert!(is_on_curve_4096(&p, &a, &b, &n));
        for mult in [3u64, 6, 9, 150] {
            let res = scalar_mul_4096(p, &U4096::from_u64(mult), &a, &n);
            assert!(
                matches!(res, PointOpResult4096::Success(Point4096::Infinity)),
                "[{mult}]P should be the identity on y^2 = x^3 + 4"
            );
        }
        let not_identity = scalar_mul_4096(p, &U4096::from_u64(4), &a, &n);
        assert!(matches!(
            not_identity,
            PointOpResult4096::Success(Point4096::Affine(_, _))
        ));

        // Moltiplicazione scalare con scalare di 192 bit: circa 191 raddoppi
        // più le additions date dai bit impostati.
        // Scalare di esattamente 192 bit, costruito parola per parola perché
        // `from_seed` produce un valore a piena larghezza.
        let mut scalar = U4096::ZERO;
        scalar.limbs[0] = 0x9e37_79b9_7f4a_7c15;
        scalar.limbs[1] = 0xbf58_476d_1ce4_e5b9;
        scalar.limbs[2] = 0x94d0_49bb_1331_11eb | (1u64 << 63);
        assert_eq!(scalar.bits(), 192);
        let triple = U4096::from_u64(3);
        // Prodotto diretto senza riduzione: il risultato resta ben sotto `n`
        // (4096 bit) e conserva la divisibilità per 3, che è ciò che l'asserzione
        // verifica alla fine della moltiplicazione scalare.
        let scalar = scalar * triple;
        let t_scalar = Instant::now();
        let res = scalar_mul_4096(p, &scalar, &a, &n);
        let scalar_ms = t_scalar.elapsed().as_secs_f64() * 1000.0;
        assert!(
            matches!(res, PointOpResult4096::Success(Point4096::Infinity)),
            "[k]P with 3 | k must be the identity"
        );

        // Prova di primalità a 4096 bit: qui il candidato è composto per
        // costruzione, interessa il costo di un giro completo di Miller-Rabin.
        let candidate = n
            .sub(&U4096::ONE)
            .0
            .mul_mod(&U4096::from_u64(65537), &n)
            .add(&U4096::ONE)
            .0;
        let t_mr = Instant::now();
        let probable = is_prime_u4096_with_bases(&candidate, 4);
        let mr_ms = t_mr.elapsed().as_secs_f64() * 1000.0;

        println!("4096-bit mul_mod      : {mul_nanos:.0} ns per operation");
        println!("4096-bit scalar mul   : {scalar_ms:.1} ms for a 192-bit scalar");
        println!(
            "4096-bit Miller-Rabin : {mr_ms:.1} ms with 4 bases (probable prime = {probable})"
        );
    }

    /// Fine a 4096 bit veri: cerca un primo a piena larghezza, lo certifica con
    /// le quaranta basi, ne classifica una curva CM e verifica che l'ordine in
    /// forma chiusa porti davvero all'infinito i punti campionati. La ricerca e
    /// la parte costosa: a questa larghezza compare un primo ogni circa 950
    /// candidati dispari. Ignorato perché chiede qualche decina di secondi.
    #[test]
    #[ignore = "searching for a 4096-bit prime: tens of seconds"]
    fn full_width_4096_bit_curve_order_is_certified() {
        use crate::primality::{is_prime_u4096, is_prime_u4096_with_bases};
        use std::time::Instant;

        let four = U4096::from_u64(4);
        let t_search = Instant::now();
        let mut cand = U4096::from_seed(0x100d_5e37_243f_6a88);
        cand.set_bit(4095, true);
        cand.limbs[0] |= 1;

        // Serve un primo per cui una delle due famiglie CM ha ordine in forma
        // chiusa: j = 0 chiede p = 2 (mod 3), j = 1728 chiede p = 3 (mod 4).
        let mut p: Option<U4096> = None;
        let mut candidates = 0usize;
        for _ in 0..8000usize {
            candidates += 1;
            let usable = cand.rem(&U4096::THREE) == U4096::TWO || cand.rem(&four) == U4096::THREE;
            if usable && is_prime_u4096_with_bases(&cand, 4) {
                p = Some(cand);
                break;
            }
            cand = cand.add(&U4096::TWO).0;
        }
        let p = p.expect("a usable 4096-bit probable prime within 8000 candidates");
        let search_secs = t_search.elapsed().as_secs_f64();

        let hex = p.to_hex_string();
        println!(
            "4096-bit prime search : {candidates} candidates in {search_secs:.2} s (p = {}...{})",
            &hex[..16],
            &hex[hex.len() - 16..]
        );

        let t_mr = Instant::now();
        let certified = is_prime_u4096(&p);
        let mr_ms = t_mr.elapsed().as_secs_f64() * 1000.0;
        println!("4096-bit Miller-Rabin : 40 bases in {mr_ms:.1} ms (prime = {certified})");
        assert!(certified, "the found prime must pass the 40-base run");

        // Famiglia CM in base ai residui del primo trovato.
        let (a, b, family) = if p.rem(&U4096::THREE) == U4096::TWO {
            (U4096::ZERO, U4096::ONE, "j = 0")
        } else {
            (U4096::ONE, U4096::ZERO, "j = 1728")
        };
        let curve = classify_curve_4096(&a, &b, &p).expect("regular curve over a prime");
        let expected = match &curve.class {
            CurveClass4096::CmSupersingular { closed_form_order } => *closed_form_order,
            other => panic!("expected a CM curve, got {}", other.name()),
        };
        assert_eq!(
            expected,
            p.add(&U4096::ONE).0,
            "closed form order must be p + 1"
        );

        let t_walk = Instant::now();
        let points = sample_points_4096(&curve, 2);
        let sample_ms = t_walk.elapsed().as_secs_f64() * 1000.0;
        assert_eq!(points.len(), 2, "sampling must work at 4096 bit ({family})");
        let mut scalars_ms = 0.0;
        for (x, y) in &points {
            let t = Instant::now();
            let res = scalar_mul_4096(Point4096::Affine(*x, *y), &expected, &curve.a, &p);
            scalars_ms += t.elapsed().as_secs_f64() * 1000.0;
            assert!(
                matches!(res, PointOpResult4096::Success(Point4096::Infinity)),
                "[p+1] P must be the identity on the {family} curve over a 4096-bit prime"
            );
        }
        println!(
            "4096-bit curve order  : 2 sampled points in {sample_ms:.0} ms, [p+1]P = O in {scalars_ms:.0} ms total"
        );
    }
}
