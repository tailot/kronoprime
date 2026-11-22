//! 4096-bit unsigned integer arithmetic module (`U4096`).
//!
//! Provides a zero-dependency, stack-allocated 4096-bit integer representation
//! with 64 limbs of 64 bits each (`[u64; 64]`).
//! Supports full arithmetic, bitwise shifts, division, modular arithmetic,
//! Tonelli-Shanks modular square roots, Stein's binary GCD, and binary extended GCD.

use std::cmp::Ordering;
use std::fmt;

/// Number of 64-bit limbs needed for 4096 bits.
pub const LIMBS_4096: usize = 64;

/// Double-width limb count for wide 8192-bit multiplication products.
pub const LIMBS_8192: usize = 128;

/// Unsigned 4096-bit integer, stack-allocated, little-endian limbs.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct U4096 {
    pub limbs: [u64; LIMBS_4096],
}

/// Unsigned 8192-bit intermediate product for overflow-free modular reductions.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct U8192 {
    pub limbs: [u64; LIMBS_8192],
}

// ---------------------------------------------------------------------------
// Limb-level helpers shared by the division core and by Montgomery reduction.
// All of them work on little-endian limb slices and never allocate.
// ---------------------------------------------------------------------------

/// Scratch buffer size for the division and Montgomery cores: a double-width
/// product plus one extra carry limb (`2 * 64 + 1` words).
const WIDE_BUFFER_WORDS: usize = LIMBS_8192 + 1;

/// Lexicographic (most significant limb first) comparison of two equal-length
/// little-endian limb slices.
fn cmp_limbs(a: &[u64], b: &[u64]) -> Ordering {
    debug_assert_eq!(a.len(), b.len());
    for i in (0..a.len()).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => continue,
            other => return other,
        }
    }
    Ordering::Equal
}

/// `a -= b` over equal-length slices, returning the final borrow (0 or 1).
fn sub_limbs(a: &mut [u64], b: &[u64]) -> u64 {
    debug_assert_eq!(a.len(), b.len());
    let mut borrow = 0u64;
    for i in 0..a.len() {
        let (d1, c1) = a[i].overflowing_sub(b[i]);
        let (d2, c2) = d1.overflowing_sub(borrow);
        a[i] = d2;
        borrow = (c1 as u64) + (c2 as u64);
    }
    borrow
}

/// Knuth's Algorithm D, inner loop.
///
/// `u` is an `(m + 1)`-slot buffer holding the normalized dividend (`m`
/// significant words plus the normalization carry in `u[m]`) and `dn` a
/// normalized `n`-word divisor, i.e. `n >= 2` and the top bit of `dn[n - 1]`
/// set. On return the normalized remainder sits in `u[..n]` and the quotient
/// digits have been written to `q` when that buffer is provided. The digit
/// estimate is at most one too large, and that single case is fixed by the
/// standard add-back step.
fn knuth_d(u: &mut [u64], m: usize, dn: &[u64], n: usize, mut q: Option<&mut [u64]>) {
    debug_assert!(n >= 2 && m >= n && dn[n - 1] != 0);
    let d1 = dn[n - 1] as u128;
    let d2 = dn[n - 2] as u128;

    // One digit per position, including the leading position that may come out
    // zero when the normalized dividend is shorter than the divisor shifted up.
    for i in (0..=(m - n)).rev() {
        // Estimate the quotient digit from the two leading dividend words.
        let num = ((u[i + n] as u128) << 64) | (u[i + n - 1] as u128);
        let (mut qhat, mut rhat): (u128, u128) = if (u[i + n] as u128) >= d1 {
            (u128::from(u64::MAX), num - d1 * u128::from(u64::MAX))
        } else {
            (num / d1, num % d1)
        };
        while rhat < d1 && qhat * d2 > ((rhat << 64) | (u[i + n - 2] as u128)) {
            qhat -= 1;
            rhat += d1;
        }

        // u[i..=i+n] -= qhat * dn
        let mut borrow = 0u64;
        for j in 0..n {
            let prod = qhat * (dn[j] as u128) + (borrow as u128);
            let (diff, b1) = u[i + j].overflowing_sub(prod as u64);
            u[i + j] = diff;
            // The owed amount never exceeds qhat, so this addition cannot wrap.
            borrow = (prod >> 64) as u64 + (b1 as u64);
        }
        let (top, b2) = u[i + n].overflowing_sub(borrow);
        u[i + n] = top;

        if b2 {
            // Estimate too large: restore by adding the divisor back.
            let mut carry = 0u64;
            for j in 0..n {
                let (s1, c1) = u[i + j].overflowing_add(dn[j]);
                let (s2, c2) = s1.overflowing_add(carry);
                u[i + j] = s2;
                carry = (c1 as u64) + (c2 as u64);
            }
            u[i + n] = u[i + n].wrapping_add(carry);
            if let Some(qbuf) = q.as_deref_mut() {
                qbuf[i] = (qhat - 1) as u64;
            }
        } else if let Some(qbuf) = q.as_deref_mut() {
            qbuf[i] = qhat as u64;
        }
    }
}

/// Complete normalized division built on [`knuth_d`].
///
/// `dividend` has `m >= n` significant words, `divisor` has `n >= 2` and both
/// buffers may be longer. The remainder is written to `remainder[..n]`.
fn div_knuth(
    dividend: &[u64],
    m: usize,
    divisor: &[u64],
    n: usize,
    q: Option<&mut [u64]>,
    remainder: &mut [u64],
) {
    debug_assert!(n >= 2 && m >= n && divisor[n - 1] != 0 && m < WIDE_BUFFER_WORDS);
    let mut u = [0u64; WIDE_BUFFER_WORDS];
    let mut dn = [0u64; LIMBS_4096];
    let s = divisor[n - 1].leading_zeros();

    if s == 0 {
        u[..m].copy_from_slice(&dividend[..m]);
        dn[..n].copy_from_slice(&divisor[..n]);
    } else {
        let mut carry = 0u64;
        for i in 0..m {
            let w = dividend[i];
            u[i] = (w << s) | carry;
            carry = w >> (64 - s);
        }
        u[m] = carry;
        let mut dcarry = 0u64;
        for i in 0..n {
            let w = divisor[i];
            dn[i] = (w << s) | dcarry;
            dcarry = w >> (64 - s);
        }
        debug_assert_eq!(
            dcarry, 0,
            "normalized divisor must keep the same word count"
        );
    }

    knuth_d(&mut u, m, &dn, n, q);

    if s == 0 {
        remainder[..n].copy_from_slice(&u[..n]);
    } else {
        let c = 64 - s;
        for i in 0..n {
            remainder[i] = (u[i] >> s) | (u[i + 1] << c);
        }
    }
}

/// `(x + m) / 2` for an odd modulus and an odd `x < m`, computed as
/// `m - (m - x) / 2` so that the intermediate sum never leaves the range.
fn half_after_add_mod(x: &U4096, m: &U4096) -> U4096 {
    debug_assert!(m.is_odd() && x.is_odd() && x < m);
    let half = m.sub(x).0.shr1();
    m.sub(&half).0
}

/// Multiplicative inverse of an odd `u64` modulo `2^64`, by Newton iteration
/// (`inv <- inv * (2 - x * inv)`), which doubles the number of correct bits
/// each step: 1 -> 2 -> 4 -> ... -> 64.
fn inverse_mod_2_64(x: u64) -> u64 {
    debug_assert!(x & 1 == 1, "only odd values are invertible modulo 2^64");
    let mut inv = 1u64;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(inv)));
    }
    inv
}

impl U4096 {
    pub const ZERO: Self = Self {
        limbs: [0; LIMBS_4096],
    };

    pub const ONE: Self = {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = 1;
        Self { limbs }
    };

    pub const TWO: Self = {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = 2;
        Self { limbs }
    };

    pub const THREE: Self = {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = 3;
        Self { limbs }
    };

    pub const FOUR: Self = {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = 4;
        Self { limbs }
    };

    pub const TWENTY_SEVEN: Self = {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = 27;
        Self { limbs }
    };

    /// Creates a `U4096` from a `u64`.
    pub const fn from_u64(val: u64) -> Self {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = val;
        Self { limbs }
    }

    /// Creates a `U4096` from a `u128`.
    pub const fn from_u128(val: u128) -> Self {
        let mut limbs = [0; LIMBS_4096];
        limbs[0] = val as u64;
        limbs[1] = (val >> 64) as u64;
        Self { limbs }
    }

    /// Converts to `u64` if it fits.
    pub fn to_u64(&self) -> Option<u64> {
        for &l in &self.limbs[1..] {
            if l != 0 {
                return None;
            }
        }
        Some(self.limbs[0])
    }

    /// Converts to `u128` if it fits.
    pub fn to_u128(&self) -> Option<u128> {
        for &l in &self.limbs[2..] {
            if l != 0 {
                return None;
            }
        }
        Some((self.limbs[0] as u128) | ((self.limbs[1] as u128) << 64))
    }

    /// True if the value is zero.
    pub fn is_zero(&self) -> bool {
        self.limbs.iter().all(|&l| l == 0)
    }

    /// True if the value is one.
    pub fn is_one(&self) -> bool {
        self.limbs[0] == 1 && self.limbs[1..].iter().all(|&l| l == 0)
    }

    /// True if the value is even.
    pub fn is_even(&self) -> bool {
        (self.limbs[0] & 1) == 0
    }

    /// True if the value is odd.
    pub fn is_odd(&self) -> bool {
        (self.limbs[0] & 1) == 1
    }

    /// Bit length of the integer (index of highest set bit + 1).
    pub fn bits(&self) -> usize {
        for i in (0..LIMBS_4096).rev() {
            if self.limbs[i] != 0 {
                return (i * 64) + (64 - self.limbs[i].leading_zeros() as usize);
            }
        }
        0
    }

    /// Number of active 64-bit words.
    pub fn active_words(&self) -> usize {
        for i in (0..LIMBS_4096).rev() {
            if self.limbs[i] != 0 {
                return i + 1;
            }
        }
        0
    }

    /// Value of bit at `index` (0-indexed).
    pub fn bit(&self, index: usize) -> bool {
        if index >= 4096 {
            return false;
        }
        let word = index / 64;
        let shift = index % 64;
        ((self.limbs[word] >> shift) & 1) == 1
    }

    /// Sets or clears bit at `index`.
    pub fn set_bit(&mut self, index: usize, val: bool) {
        if index >= 4096 {
            return;
        }
        let word = index / 64;
        let shift = index % 64;
        if val {
            self.limbs[word] |= 1u64 << shift;
        } else {
            self.limbs[word] &= !(1u64 << shift);
        }
    }

    /// Number of leading zero bits.
    pub fn leading_zeros(&self) -> u32 {
        4096 - self.bits() as u32
    }

    /// Number of trailing zero bits.
    pub fn trailing_zeros(&self) -> u32 {
        let mut count = 0u32;
        for &l in &self.limbs {
            if l == 0 {
                count += 64;
            } else {
                count += l.trailing_zeros();
                break;
            }
        }
        count.min(4096)
    }

    /// Addition with carry `(self + other)`. Returns `(result, carry)`.
    pub fn add(&self, other: &Self) -> (Self, bool) {
        let mut result = Self::ZERO;
        let mut carry = 0u64;
        for i in 0..LIMBS_4096 {
            let (sum1, c1) = self.limbs[i].overflowing_add(other.limbs[i]);
            let (sum2, c2) = sum1.overflowing_add(carry);
            result.limbs[i] = sum2;
            carry = (c1 as u64) + (c2 as u64);
        }
        (result, carry != 0)
    }

    /// Subtraction with borrow `(self - other)`. Returns `(result, borrow)`.
    pub fn sub(&self, other: &Self) -> (Self, bool) {
        let mut result = Self::ZERO;
        let mut borrow = 0u64;
        for i in 0..LIMBS_4096 {
            let (diff1, b1) = self.limbs[i].overflowing_sub(other.limbs[i]);
            let (diff2, b2) = diff1.overflowing_sub(borrow);
            result.limbs[i] = diff2;
            borrow = (b1 as u64) + (b2 as u64);
        }
        (result, borrow != 0)
    }

    /// Shift right by 1 bit.
    pub fn shr1(&self) -> Self {
        let mut res = Self::ZERO;
        let mut carry = 0u64;
        for i in (0..LIMBS_4096).rev() {
            let next_carry = self.limbs[i] & 1;
            res.limbs[i] = (self.limbs[i] >> 1) | (carry << 63);
            carry = next_carry;
        }
        res
    }

    /// Shift left by 1 bit.
    pub fn shl1(&self) -> (Self, bool) {
        let mut res = Self::ZERO;
        let mut carry = 0u64;
        for i in 0..LIMBS_4096 {
            let next_carry = (self.limbs[i] >> 63) & 1;
            res.limbs[i] = (self.limbs[i] << 1) | carry;
            carry = next_carry;
        }
        (res, carry != 0)
    }

    /// General shift right by `shift` bits.
    pub fn shr(&self, shift: usize) -> Self {
        if shift >= 4096 {
            return Self::ZERO;
        }
        let word_shift = shift / 64;
        let bit_shift = shift % 64;
        let mut res = Self::ZERO;

        if bit_shift == 0 {
            for i in word_shift..LIMBS_4096 {
                res.limbs[i - word_shift] = self.limbs[i];
            }
        } else {
            let complement = 64 - bit_shift;
            for i in word_shift..LIMBS_4096 {
                let low = self.limbs[i] >> bit_shift;
                let high = if i + 1 < LIMBS_4096 {
                    self.limbs[i + 1] << complement
                } else {
                    0
                };
                res.limbs[i - word_shift] = low | high;
            }
        }
        res
    }

    /// General shift left by `shift` bits.
    pub fn shl(&self, shift: usize) -> Self {
        if shift >= 4096 {
            return Self::ZERO;
        }
        let word_shift = shift / 64;
        let bit_shift = shift % 64;
        let mut res = Self::ZERO;

        if bit_shift == 0 {
            for i in 0..(LIMBS_4096 - word_shift) {
                res.limbs[i + word_shift] = self.limbs[i];
            }
        } else {
            let complement = 64 - bit_shift;
            for i in 0..(LIMBS_4096 - word_shift) {
                let high = self.limbs[i] << bit_shift;
                let low = if i > 0 {
                    self.limbs[i - 1] >> complement
                } else {
                    0
                };
                res.limbs[i + word_shift] = high | low;
            }
        }
        res
    }

    /// Multiply by single `u64`. Returns `(result, carry)`.
    pub fn mul_u64(&self, other: u64) -> (Self, u64) {
        let mut res = Self::ZERO;
        let mut carry = 0u128;
        let o = other as u128;

        for i in 0..LIMBS_4096 {
            let prod = (self.limbs[i] as u128) * o + carry;
            res.limbs[i] = prod as u64;
            carry = prod >> 64;
        }
        (res, carry as u64)
    }

    /// Division and remainder by single `u64`. Returns `(quotient, remainder)`.
    pub fn div_rem_u64(&self, d: u64) -> (Self, u64) {
        assert!(d != 0, "division by zero in div_rem_u64");
        let mut quot = Self::ZERO;
        let mut rem = 0u128;
        let d128 = d as u128;

        for i in (0..LIMBS_4096).rev() {
            let cur = (rem << 64) | (self.limbs[i] as u128);
            let q = cur / d128;
            rem = cur % d128;
            quot.limbs[i] = q as u64;
        }
        (quot, rem as u64)
    }

    /// Double-width 8192-bit multiplication product.
    pub fn mul_wide(&self, other: &Self) -> U8192 {
        let mut out = [0u64; LIMBS_8192];
        let self_words = self.active_words();
        let other_words = other.active_words();

        for i in 0..self_words {
            let a = self.limbs[i] as u128;
            if a == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..other_words {
                let b = other.limbs[j] as u128;
                let cur = out[i + j] as u128;
                let prod = a * b + cur + carry;
                out[i + j] = prod as u64;
                carry = prod >> 64;
            }
            let mut k = i + other_words;
            while carry > 0 && k < LIMBS_8192 {
                let cur = out[k] as u128;
                let sum = cur + carry;
                out[k] = sum as u64;
                carry = sum >> 64;
                k += 1;
            }
        }
        U8192 { limbs: out }
    }

    /// Division and remainder `(quotient, remainder)`.
    ///
    /// Uses Knuth's Algorithm D over the active words, so the cost is
    /// quadratic in the number of significant limbs (and linear for a
    /// single-limb divisor) instead of one shift-subtract per bit.
    pub fn div_rem(&self, d: &Self) -> Result<(Self, Self), &'static str> {
        if d.is_zero() {
            return Err("Division by zero");
        }
        if self < d {
            return Ok((Self::ZERO, *self));
        }
        if d.active_words() == 1 {
            let (q, r) = self.div_rem_u64(d.limbs[0]);
            return Ok((q, Self::from_u64(r)));
        }

        let m = self.active_words();
        let n = d.active_words();
        let mut quot = [0u64; LIMBS_4096];
        let mut rem = [0u64; LIMBS_4096];
        div_knuth(&self.limbs, m, &d.limbs, n, Some(&mut quot), &mut rem);
        Ok((Self { limbs: quot }, Self { limbs: rem }))
    }

    /// Quotient of the division, panicking on a zero divisor.
    pub fn quo(&self, d: &Self) -> Self {
        self.div_rem(d).expect("division by zero in U4096::quo").0
    }

    /// Remainder of the division, panicking on a zero divisor.
    pub fn rem(&self, d: &Self) -> Self {
        self.div_rem(d).expect("division by zero in U4096::rem").1
    }

    /// Modular addition `(self + other) mod m`.
    ///
    /// Both operands must already be reduced (`self, other < m`), which is the
    /// invariant every modular code path in this crate maintains; a single
    /// conditional subtraction is therefore enough.
    pub fn add_mod(&self, other: &Self, m: &Self) -> Self {
        debug_assert!(self < m && other < m, "add_mod expects reduced operands");
        let (sum, carry) = self.add(other);
        if carry || sum >= *m {
            sum.sub(m).0
        } else {
            sum
        }
    }

    /// Additive inverse `(-self) mod m`.
    pub fn neg_mod(&self, m: &Self) -> Self {
        if self.is_zero() {
            Self::ZERO
        } else {
            m.sub(self).0
        }
    }

    /// Modular subtraction `(self - other) mod m`.
    ///
    /// Correct for arbitrary operands: the difference is reduced when it is not
    /// already below the modulus.
    pub fn sub_mod(&self, other: &Self, m: &Self) -> Self {
        debug_assert!(self < m && other < m, "sub_mod expects reduced operands");
        if self >= other {
            self.sub(other).0
        } else {
            let diff = other.sub(self).0;
            if diff.is_zero() {
                Self::ZERO
            } else {
                m.sub(&diff).0
            }
        }
    }

    /// Modular multiplication `(self * other) mod m` via the double-width
    /// 8192-bit product and one Algorithm D reduction.
    pub fn mul_mod(&self, other: &Self, m: &Self) -> Self {
        debug_assert!(self < m && other < m, "mul_mod expects reduced operands");
        let wide = self.mul_wide(other);
        wide.rem_u4096(m)
    }

    /// Modular exponentiation `self^exp mod m`.
    ///
    /// Odd moduli go through a Montgomery context (one `R = 2^(64n)` setup and
    /// then purely `n^2` multiplications); even moduli fall back to the plain
    /// square-and-multiply path over `mul_mod`.
    pub fn pow_mod(&self, exp: &Self, m: &Self) -> Self {
        if m.is_one() {
            return Self::ZERO;
        }
        if m.is_odd() {
            if let Some(ctx) = MontCtx::new(m) {
                return ctx.pow(self, exp);
            }
        }
        let mut result = Self::ONE;
        let mut base = if *self >= *m { self.rem(m) } else { *self };
        for i in 0..exp.bits() {
            if exp.bit(i) {
                result = result.mul_mod(&base, m);
            }
            base = base.mul_mod(&base, m);
        }
        result
    }

    /// Greatest common divisor via Stein's binary GCD algorithm.
    pub fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() {
            return *other;
        }
        if other.is_zero() {
            return *self;
        }

        let shift_self = self.trailing_zeros();
        let shift_other = other.trailing_zeros();
        let k = shift_self.min(shift_other);

        let mut u = self.shr(shift_self as usize);
        let mut v = other.shr(shift_other as usize);

        loop {
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            v = v.sub(&u).0;
            if v.is_zero() {
                break;
            }
            v = v.shr(v.trailing_zeros() as usize);
        }

        u.shl(k as usize)
    }

    /// Modular inversion.
    ///
    /// Returns `Ok(inv)` with `0 < inv < m` when `gcd(self, m) == 1`, or
    /// `Err(g)` carrying the common divisor `g = gcd(self, m)`, which is the
    /// signal the elliptic-curve factor probe relies on.
    ///
    /// Odd moduli use the binary extended GCD (only shifts, adds and
    /// subtractions). Every coefficient update is kept reduced modulo `m`, and
    /// the "add `m` before halving" step is computed as `m - (m - x) / 2`, so
    /// nothing ever overflows the 4096-bit range even for a modulus with the
    /// top bit set. Even moduli go through the Euclid path in [`Self::inv_mod_euclid`].
    //
    // The error carries the gcd itself, a whole 4096-bit word: boxing it would
    // push a heap allocation into the number path, which this crate avoids.
    #[allow(clippy::result_large_err)]
    pub fn try_inv_mod(&self, m: &Self) -> Result<Self, Self> {
        if *m <= Self::ONE {
            return Err(*m);
        }
        let a = self.rem(m);
        if a.is_zero() {
            return Err(*m);
        }
        if m.is_even() {
            return self.inv_mod_euclid(m);
        }
        let g = a.gcd(m);
        if g != Self::ONE {
            return Err(g);
        }

        let mut u = a;
        let mut v = *m;
        let mut x1 = Self::ONE;
        let mut x2 = Self::ZERO;

        while u != Self::ONE && v != Self::ONE {
            while u.is_even() {
                u = u.shr1();
                // Half x1 modulo m: an odd coefficient first picks up m.
                x1 = if x1.is_odd() {
                    half_after_add_mod(&x1, m)
                } else {
                    x1.shr1()
                };
            }
            while v.is_even() {
                v = v.shr1();
                x2 = if x2.is_odd() {
                    half_after_add_mod(&x2, m)
                } else {
                    x2.shr1()
                };
            }
            if u >= v {
                u = u.sub(&v).0;
                x1 = x1.sub_mod(&x2, m);
            } else {
                v = v.sub(&u).0;
                x2 = x2.sub_mod(&x1, m);
            }
        }

        if u == Self::ONE {
            Ok(x1.rem(m))
        } else {
            Ok(x2.rem(m))
        }
    }

    /// Modular inversion through the Euclidean algorithm, keeping the B\'ezout
    /// coefficient reduced modulo `m` at every step.
    ///
    /// This path costs one full division per Euclid step, but it needs no
    /// parity assumption on the modulus, so it is the fallback for even `m`.
    //
    // Same stack-only reason as [`Self::try_inv_mod`]: the error is a residue.
    #[allow(clippy::result_large_err)]
    fn inv_mod_euclid(&self, m: &Self) -> Result<Self, Self> {
        let a = self.rem(m);
        if a.is_zero() {
            return Err(*m);
        }
        let mut r0 = *m;
        let mut r1 = a;
        let mut t0 = Self::ZERO;
        let mut t1 = Self::ONE;

        while !r1.is_zero() {
            let (q, r2) = r0.div_rem(&r1).expect("non-zero divisor in inv_mod_euclid");
            r0 = r1;
            r1 = r2;
            // q < r0 <= m, so the product stays inside the modular range.
            let step = q.mul_mod(&t1, m);
            t0 = t1;
            t1 = t0.sub_mod(&step, m);
        }

        if r0 != Self::ONE {
            return Err(r0);
        }
        Ok(t0)
    }

    /// Computes modular inverse or returns 0.
    pub fn inv_mod(&self, m: &Self) -> Self {
        self.try_inv_mod(m).unwrap_or(Self::ZERO)
    }

    /// Deterministic pseudo-random value from a 64-bit seed (SplitMix64).
    ///
    /// Curve sampling and the elliptic-curve method both need repeatable
    /// "random" operands without pulling in a dependency, and a hash of the
    /// seed into every limb gives exactly that: the same seed always replays
    /// the same search, which is what makes a failed factor probe reproducible.
    pub fn from_seed(seed: u64) -> Self {
        const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut state = seed;
        let mut limbs = [0u64; LIMBS_4096];
        for w in limbs.iter_mut() {
            state = state.wrapping_add(GOLDEN);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            *w = z ^ (z >> 31);
        }
        Self { limbs }
    }

    /// A pseudo-random residue in `0..m` derived from `seed`.
    pub fn sample_below(seed: u64, m: &Self) -> Self {
        Self::from_seed(seed).rem(m)
    }

    /// Parses from decimal string.
    pub fn from_dec(s: &str) -> Result<Self, &'static str> {
        let s = s.trim();
        if s.is_empty() {
            return Err("Stringa vuota");
        }
        let mut result = Self::ZERO;
        for ch in s.chars() {
            if !ch.is_ascii_digit() {
                return Err("Carattere decimale non valido");
            }
            let digit = (ch as u8 - b'0') as u64;
            let (prod, overflow) = result.mul_u64(10);
            if overflow != 0 {
                return Err("Overflow oltre 4096 bit");
            }
            let (sum, carry) = prod.add(&Self::from_u64(digit));
            if carry {
                return Err("Overflow oltre 4096 bit");
            }
            result = sum;
        }
        Ok(result)
    }

    /// Parses from hex string (with or without optional "0x" prefix).
    pub fn from_hex(s: &str) -> Result<Self, &'static str> {
        let s = s.trim();
        let s = s
            .strip_prefix("0x")
            .or_else(|| s.strip_prefix("0X"))
            .unwrap_or(s);
        if s.is_empty() {
            return Err("Stringa esadecimale vuota");
        }
        if s.len() > 1024 {
            return Err("Stringa esadecimale supera 4096 bit (1024 cifre hex)");
        }
        let mut result = Self::ZERO;
        for (i, ch) in s.chars().rev().enumerate() {
            let val = ch.to_digit(16).ok_or("Carattere esadecimale non valido")? as u64;
            let word = i / 16;
            let shift = (i % 16) * 4;
            if word < LIMBS_4096 {
                result.limbs[word] |= val << shift;
            }
        }
        Ok(result)
    }

    /// Parses from string (auto-detecting hex with "0x" prefix or decimal).
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Result<Self, &'static str> {
        let s = s.trim();
        if s.starts_with("0x") || s.starts_with("0X") {
            Self::from_hex(s)
        } else {
            Self::from_dec(s)
        }
    }

    /// Formats as decimal string.
    pub fn to_dec_string(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let mut chunks = Vec::new();
        let mut val = *self;
        let base = 1_000_000_000_000_000_000u64; // 10^18

        while !val.is_zero() {
            let (q, rem) = val.div_rem_u64(base);
            chunks.push(rem);
            val = q;
        }

        let mut out = String::new();
        for (i, &chunk) in chunks.iter().rev().enumerate() {
            if i == 0 {
                out.push_str(&chunk.to_string());
            } else {
                out.push_str(&format!("{chunk:018}"));
            }
        }
        out
    }

    /// Formats as hex string with "0x" prefix.
    pub fn to_hex_string(&self) -> String {
        if self.is_zero() {
            return "0x0".to_string();
        }
        let mut hex_digits = String::new();
        let mut started = false;
        for i in (0..LIMBS_4096).rev() {
            let word = self.limbs[i];
            if !started {
                if word != 0 {
                    started = true;
                    hex_digits.push_str(&format!("{word:x}"));
                }
            } else {
                hex_digits.push_str(&format!("{word:016x}"));
            }
        }
        format!("0x{hex_digits}")
    }
}

impl PartialOrd for U4096 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for U4096 {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..LIMBS_4096).rev() {
            if self.limbs[i] < other.limbs[i] {
                return Ordering::Less;
            }
            if self.limbs[i] > other.limbs[i] {
                return Ordering::Greater;
            }
        }
        Ordering::Equal
    }
}

impl fmt::Display for U4096 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_dec_string())
    }
}

impl fmt::Debug for U4096 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_dec_string())
    }
}

impl std::str::FromStr for U4096 {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_str(s)
    }
}

impl std::ops::Add for U4096 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        U4096::add(&self, &rhs).0
    }
}

impl std::ops::Sub for U4096 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        U4096::sub(&self, &rhs).0
    }
}

impl std::ops::Mul for U4096 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let wide = self.mul_wide(&rhs);
        let mut out = [0u64; LIMBS_4096];
        out.copy_from_slice(&wide.limbs[..LIMBS_4096]);
        Self { limbs: out }
    }
}

impl std::ops::Div for U4096 {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.div_rem(&rhs)
            .expect("division by zero in U4096::div")
            .0
    }
}

impl std::ops::Rem for U4096 {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self {
        self.div_rem(&rhs)
            .expect("division by zero in U4096::rem")
            .1
    }
}

impl U8192 {
    /// Creates a `U8192` from `U4096`.
    pub fn from_u4096(val: &U4096) -> Self {
        let mut limbs = [0u64; LIMBS_8192];
        limbs[..LIMBS_4096].copy_from_slice(&val.limbs);
        Self { limbs }
    }

    /// Bit length of 8192-bit value.
    pub fn bits(&self) -> usize {
        for i in (0..LIMBS_8192).rev() {
            if self.limbs[i] != 0 {
                return (i * 64) + (64 - self.limbs[i].leading_zeros() as usize);
            }
        }
        0
    }

    /// Number of active 64-bit words.
    pub fn active_words(&self) -> usize {
        for i in (0..LIMBS_8192).rev() {
            if self.limbs[i] != 0 {
                return i + 1;
            }
        }
        0
    }

    /// Shift right by 1 bit.
    pub fn shr1(&self) -> Self {
        let mut res = Self {
            limbs: [0; LIMBS_8192],
        };
        let mut carry = 0u64;
        for i in (0..LIMBS_8192).rev() {
            let next_carry = self.limbs[i] & 1;
            res.limbs[i] = (self.limbs[i] >> 1) | (carry << 63);
            carry = next_carry;
        }
        res
    }

    /// General shift left by `shift` bits.
    pub fn shl(&self, shift: usize) -> Self {
        if shift >= 8192 {
            return Self {
                limbs: [0; LIMBS_8192],
            };
        }
        let word_shift = shift / 64;
        let bit_shift = shift % 64;
        let mut res = Self {
            limbs: [0; LIMBS_8192],
        };

        if bit_shift == 0 {
            for i in 0..(LIMBS_8192 - word_shift) {
                res.limbs[i + word_shift] = self.limbs[i];
            }
        } else {
            let complement = 64 - bit_shift;
            for i in 0..(LIMBS_8192 - word_shift) {
                let high = self.limbs[i] << bit_shift;
                let low = if i > 0 {
                    self.limbs[i - 1] >> complement
                } else {
                    0
                };
                res.limbs[i + word_shift] = high | low;
            }
        }
        res
    }

    /// Subtraction `(self - other)`.
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = Self {
            limbs: [0; LIMBS_8192],
        };
        let mut borrow = 0u64;
        for i in 0..LIMBS_8192 {
            let (diff1, b1) = self.limbs[i].overflowing_sub(other.limbs[i]);
            let (diff2, b2) = diff1.overflowing_sub(borrow);
            result.limbs[i] = diff2;
            borrow = (b1 as u64) + (b2 as u64);
        }
        result
    }

    pub fn cmp_wide(&self, other: &Self) -> Ordering {
        for i in (0..LIMBS_8192).rev() {
            if self.limbs[i] < other.limbs[i] {
                return Ordering::Less;
            }
            if self.limbs[i] > other.limbs[i] {
                return Ordering::Greater;
            }
        }
        Ordering::Equal
    }

    /// Reduces an 8192-bit product modulo the 4096-bit divisor `d`.
    ///
    /// A single Algorithm D pass over the active words replaces the old
    /// one-subtraction-per-bit loop, which is what makes 4096-bit modular
    /// exponentiation practical.
    pub fn rem_u4096(&self, d: &U4096) -> U4096 {
        assert!(!d.is_zero(), "modulo by zero in U8192::rem_u4096");
        let n = d.active_words();
        let m = self.active_words();
        if m == 0 || m < n {
            // Fewer active words than the divisor: the value is already below it.
            let mut out = [0u64; LIMBS_4096];
            let take = m.min(LIMBS_4096);
            out[..take].copy_from_slice(&self.limbs[..take]);
            return U4096 { limbs: out };
        }
        if n == 1 {
            let d0 = d.limbs[0] as u128;
            let mut rem = 0u128;
            for i in (0..LIMBS_8192).rev() {
                let cur = (rem << 64) | (self.limbs[i] as u128);
                rem = cur % d0;
            }
            return U4096::from_u64(rem as u64);
        }
        let mut out = [0u64; LIMBS_4096];
        div_knuth(&self.limbs, m, &d.limbs, n, None, &mut out);
        U4096 { limbs: out }
    }
}

// ---------------------------------------------------------------------------
// Montgomery multiplication
// ---------------------------------------------------------------------------

/// Montgomery context for an odd modulus: precomputes `-m^-1 mod 2^64` and
/// `R^2 mod m` for `R = 2^(64n)`, where `n` is the active word count of `m`.
///
/// Every multiplication then costs one `n x n` schoolbook product plus one
/// `n x n` reduction pass, with no comparison-per-bit and no division. Build
/// the context once per modulus and reuse it across Miller-Rabin bases and
/// across the whole square-and-multiply chain.
#[derive(Clone, Copy, Debug)]
pub struct MontCtx {
    /// Number of significant limbs of the modulus (`R = 2^(64 * n)`).
    n: usize,
    /// Modulus, in the same little-endian layout as [`U4096`].
    m: [u64; LIMBS_4096],
    /// `-m^-1 mod 2^64`, the Montgomery constant.
    m_prime: u64,
    /// `R^2 mod m`, used to move ordinary residues into Montgomery form.
    r2: [u64; LIMBS_4096],
    /// `R mod m`, the Montgomery form of `1`.
    one_mont: [u64; LIMBS_4096],
}

impl MontCtx {
    /// Builds the context, or `None` when the modulus is even, `0` or `1`.
    pub fn new(m: &U4096) -> Option<Self> {
        if m.bits() <= 1 || m.is_even() {
            return None;
        }
        let n = m.active_words();
        let mb = m.limbs;
        let m_prime = inverse_mod_2_64(mb[0]).wrapping_neg();

        // R^2 mod m by doubling: 128n doublings of 1 modulo m.
        let mut r2 = [0u64; LIMBS_4096];
        r2[0] = 1;
        let mut doubled = [0u64; LIMBS_4096];
        for _ in 0..(128 * n) {
            let mut carry = 0u64;
            for i in 0..n {
                let next = r2[i] >> 63;
                doubled[i] = (r2[i] << 1) | carry;
                carry = next;
            }
            if carry != 0 || cmp_limbs(&doubled[..n], &mb[..n]) != Ordering::Less {
                sub_limbs(&mut doubled[..n], &mb[..n]);
            }
            std::mem::swap(&mut r2, &mut doubled);
        }

        let ctx = Self {
            n,
            m: mb,
            m_prime,
            r2,
            one_mont: [0u64; LIMBS_4096],
        };
        let mut ctx = ctx;
        ctx.one_mont = ctx.mul(&ctx.r2, &Self::one_limbs());
        Some(ctx)
    }

    /// The limbs of the integer one.
    fn one_limbs() -> [u64; LIMBS_4096] {
        let mut limbs = [0u64; LIMBS_4096];
        limbs[0] = 1;
        limbs
    }

    /// The modulus this context works modulo.
    pub fn modulus(&self) -> U4096 {
        U4096 { limbs: self.m }
    }

    /// Montgomery form of `1`, i.e. `R mod m`.
    pub fn one(&self) -> [u64; LIMBS_4096] {
        self.one_mont
    }

    /// Compares a Montgomery-form residue against another Montgomery-form
    /// residue without leaving the representation.
    pub fn eq_mont(&self, a: &[u64; LIMBS_4096], b: &[u64; LIMBS_4096]) -> bool {
        cmp_limbs(&a[..self.n], &b[..self.n]) == Ordering::Equal
    }

    /// Montgomery reduction of a `(2n + 1)`-word accumulator.
    fn reduce(&self, t: &mut [u64; WIDE_BUFFER_WORDS]) -> [u64; LIMBS_4096] {
        let n = self.n;
        for i in 0..n {
            let q = ((t[i] as u128) * (self.m_prime as u128)) as u64;
            if q == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..n {
                let cur = (t[i + j] as u128) + (q as u128) * (self.m[j] as u128) + carry;
                t[i + j] = cur as u64;
                carry = cur >> 64;
            }
            let mut k = i + n;
            while carry > 0 {
                let cur = (t[k] as u128) + carry;
                t[k] = cur as u64;
                carry = cur >> 64;
                k += 1;
            }
        }

        let mut out = [0u64; LIMBS_4096];
        out[..n].copy_from_slice(&t[n..2 * n]);
        if t[2 * n] != 0 || cmp_limbs(&out[..n], &self.m[..n]) != Ordering::Less {
            sub_limbs(&mut out[..n], &self.m[..n]);
        }
        out
    }

    /// Montgomery product `a * b * R^-1 mod m` for residues already in
    /// Montgomery form.
    pub fn mul(&self, a: &[u64; LIMBS_4096], b: &[u64; LIMBS_4096]) -> [u64; LIMBS_4096] {
        let n = self.n;
        let mut t = [0u64; WIDE_BUFFER_WORDS];
        for i in 0..n {
            let ai = a[i] as u128;
            if ai == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..n {
                let cur = (t[i + j] as u128) + ai * (b[j] as u128) + carry;
                t[i + j] = cur as u64;
                carry = cur >> 64;
            }
            let mut k = i + n;
            while carry > 0 {
                let cur = (t[k] as u128) + carry;
                t[k] = cur as u64;
                carry = cur >> 64;
                k += 1;
            }
        }
        self.reduce(&mut t)
    }

    /// Moves an ordinary residue into Montgomery form (`x * R mod m`).
    pub fn to_mont(&self, x: &U4096) -> [u64; LIMBS_4096] {
        let m = U4096 { limbs: self.m };
        let reduced = if *x >= m { x.rem(&m) } else { *x };
        self.mul(&reduced.limbs, &self.r2)
    }

    /// Moves a Montgomery-form residue back to the ordinary representation.
    pub fn from_mont(&self, a: &[u64; LIMBS_4096]) -> U4096 {
        let mut t = [0u64; WIDE_BUFFER_WORDS];
        t[..self.n].copy_from_slice(&a[..self.n]);
        U4096 {
            limbs: self.reduce(&mut t),
        }
    }

    /// Square-and-multiply entirely in Montgomery form, returning the residue
    /// still in Montgomery representation (no conversion back).
    pub fn pow_mont_arr(&self, base_mont: &[u64; LIMBS_4096], exp: &U4096) -> [u64; LIMBS_4096] {
        let mut result = self.one_mont;
        let mut base = *base_mont;
        for i in 0..exp.bits() {
            if exp.bit(i) {
                result = self.mul(&result, &base);
            }
            if i + 1 < exp.bits() {
                base = self.mul(&base, &base);
            }
        }
        result
    }

    /// Square-and-multiply in Montgomery form, reduced to an ordinary residue.
    pub fn pow_mont(&self, base_mont: &[u64; LIMBS_4096], exp: &U4096) -> U4096 {
        let result = self.pow_mont_arr(base_mont, exp);
        self.from_mont(&result)
    }

    /// `base^exp mod m` for an ordinary base.
    pub fn pow(&self, base: &U4096, exp: &U4096) -> U4096 {
        let base_mont = self.to_mont(base);
        self.pow_mont(&base_mont, exp)
    }

    /// Montgomery form of a small base, ready to be exponentiated repeatedly
    /// (one `to_mont` and one `pow_mont` per Miller-Rabin base).
    pub fn pow_base_u64(&self, base: u64, exp: &U4096) -> U4096 {
        let base_mont = self.to_mont(&U4096::from_u64(base));
        self.pow_mont(&base_mont, exp)
    }
}

/// Legendre symbol `(a / p)` for an odd prime `p` at 4096 bits.
///
/// Returns `1` (quadratic residue), `-1` (non-residue) or `0` (`a = 0 mod p`).
pub fn legendre_symbol_4096(a: &U4096, p: &U4096) -> i8 {
    let a_mod = a.rem(p);
    if a_mod.is_zero() {
        return 0;
    }
    let p_minus_1 = p.sub(&U4096::ONE).0;
    let exp = p_minus_1.shr1();
    let ls = a_mod.pow_mod(&exp, p);
    if ls == U4096::ONE {
        1
    } else if ls == p_minus_1 {
        -1
    } else {
        0
    }
}

/// Square root modulo an odd prime `p` at 4096 bits (Tonelli-Shanks).
///
/// Returns `Some(r)` with `r^2 = a mod p`, or `None` when `a` is not a
/// quadratic residue. One Montgomery context is built for the modulus and
/// shared by every exponentiation of the algorithm; the `p = 3 mod 4` closed
/// form is tried first, which halves the work for the moduli the curve code
/// usually meets.
pub fn mod_sqrt_4096(a: &U4096, p: &U4096) -> Option<U4096> {
    let a_mod = a.rem(p);
    if a_mod.is_zero() {
        return Some(U4096::ZERO);
    }
    if *p == U4096::TWO {
        return Some(a_mod);
    }
    let ctx = MontCtx::new(p)?;

    // Fast path for p = 3 mod 4: r = a^((p + 1) / 4) mod p, then verify.
    if p.limbs[0] % 4 == 3 {
        let exp = p.add(&U4096::ONE).0.shr(2);
        let r = ctx.pow(&a_mod, &exp);
        return if r.mul_mod(&r, p) == a_mod {
            Some(r)
        } else {
            None
        };
    }

    // General Tonelli-Shanks: p - 1 = q * 2^s with q odd
    let mut q = p.sub(&U4096::ONE).0;
    let mut s = 0u32;
    while q.is_even() {
        q = q.shr1();
        s += 1;
    }

    // Smallest quadratic non-residue, found by trial on the small integers.
    let mut z = U4096::TWO;
    while legendre_symbol_4096(&z, p) != -1 {
        z = z.add(&U4096::ONE).0;
    }

    let mut c = ctx.pow(&z, &q);
    let q_plus_1_half = q.add(&U4096::ONE).0.shr1();
    let mut r = ctx.pow(&a_mod, &q_plus_1_half);
    let mut t = ctx.pow(&a_mod, &q);
    let mut m = s;

    while t != U4096::ONE {
        let mut i = 0u32;
        let mut temp = t;
        while temp != U4096::ONE && i < m {
            temp = temp.mul_mod(&temp, p);
            i += 1;
        }
        if i == m {
            return None;
        }
        let exp = U4096::ONE.shl((m - i - 1) as usize);
        let b = ctx.pow(&c, &exp);
        r = r.mul_mod(&b, p);
        c = b.mul_mod(&b, p);
        t = t.mul_mod(&c, p);
        m = i;
    }
    if r.mul_mod(&r, p) == a_mod {
        Some(r)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u4096_basic_arithmetic() {
        let a = U4096::from_u64(123456789);
        let b = U4096::from_u64(987654321);
        let sum = a + b;
        assert_eq!(sum.to_u64(), Some(123456789 + 987654321));

        let diff = b - a;
        assert_eq!(diff.to_u64(), Some(987654321 - 123456789));

        let prod = a * b;
        assert_eq!(prod.to_u128(), Some(123456789u128 * 987654321));

        let (q, r) = b.div_rem(&a).unwrap();
        assert_eq!(q.to_u64(), Some(987654321 / 123456789));
        assert_eq!(r.to_u64(), Some(987654321 % 123456789));
    }

    #[test]
    fn u4096_hex_and_dec_roundtrip() {
        let hex = "0x123456789abcdef0123456789abcdef0123456789abcdef0";
        let val = U4096::from_hex(hex).unwrap();
        assert_eq!(val.to_hex_string(), hex);

        let dec = "1234567890123456789012345678901234567890";
        let val2 = U4096::from_dec(dec).unwrap();
        assert_eq!(val2.to_dec_string(), dec);
    }

    #[test]
    fn u4096_modular_arithmetic_and_inversion() {
        let m = U4096::from_u64(1019);
        let a = U4096::from_u64(45);
        let inv = a.try_inv_mod(&m).expect("must be invertible");
        assert_eq!(a.mul_mod(&inv, &m), U4096::ONE);

        let g = a.gcd(&m);
        assert_eq!(g, U4096::ONE);
    }

    #[test]
    fn u4096_montgomery_matches_wide_reduction() {
        // A 320-bit odd modulus in normal form (top bit of the leading word set).
        let mut m = U4096::ZERO;
        for (i, limb) in m.limbs.iter_mut().enumerate().take(5) {
            *limb = 0x9e3779b97f4a7c15u64.wrapping_mul((i as u64) + 1) | 1;
        }
        m.limbs[4] |= 1 << 63;
        let ctx = MontCtx::new(&m).expect("odd modulus above one");
        assert_eq!(ctx.n, 5);

        // The Montgomery unit round trips, and so does any residue.
        assert_eq!(ctx.from_mont(&ctx.one()), U4096::ONE);
        let x = U4096::from_u64(0xdeadbeefcafe1234)
            .mul_wide(&U4096::from_u64(0x1234567890abcdef))
            .rem_u4096(&m);
        assert_eq!(ctx.from_mont(&ctx.to_mont(&x)), x);

        // Montgomery exponentiation must agree with the schoolbook path that
        // reduces every step through the full 8192-bit product.
        for seed in 0u64..4 {
            let base = U4096::from_seed(seed).rem(&m);
            let exp = U4096::from_seed(seed.wrapping_add(777));
            let mut acc = U4096::ONE;
            let mut b = base;
            for i in 0..exp.bits() {
                if exp.bit(i) {
                    acc = acc.mul_mod(&b, &m);
                }
                b = b.mul_mod(&b, &m);
            }
            assert_eq!(ctx.pow(&base, &exp), acc, "mismatch for seed {seed}");
        }

        // Montgomery form multiplies cleanly: mul(a~, b~) = (a*b mod m)~.
        let a = U4096::from_seed(31).rem(&m);
        let b = U4096::from_seed(99).rem(&m);
        let prod_mont = ctx.mul(&ctx.to_mont(&a), &ctx.to_mont(&b));
        assert_eq!(ctx.from_mont(&prod_mont), a.mul_mod(&b, &m));

        // Even moduli and the degenerate ones have no Montgomery form.
        assert!(MontCtx::new(&U4096::from_u64(1020)).is_none());
        assert!(MontCtx::new(&U4096::ONE).is_none());
        assert!(MontCtx::new(&U4096::ZERO).is_none());
    }

    #[test]
    fn u4096_division_invariants() {
        // a = q * b + r with r < b, recombined through the 8192-bit product.
        for seed in 0u64..12 {
            let a = U4096::from_seed(seed);
            let b = U4096::from_seed(seed.wrapping_add(0x9e37_79b9_7f4a_7c15));
            if b.is_zero() {
                continue;
            }
            let (q, r) = a.div_rem(&b).unwrap();
            assert!(r < b, "remainder not reduced for seed {seed}");
            let prod = q.mul_wide(&b);
            for w in &prod.limbs[LIMBS_4096..] {
                assert_eq!(*w, 0, "quotient too wide for seed {seed}");
            }
            let mut low = [0u64; LIMBS_4096];
            low.copy_from_slice(&prod.limbs[..LIMBS_4096]);
            let (back, carry) = U4096 { limbs: low }.add(&r);
            assert!(!carry, "recombination overflow for seed {seed}");
            assert_eq!(back, a, "division round trip failed for seed {seed}");
        }

        // Single-limb divisors must agree with the dedicated fast path.
        let d = 0x1234_5678_9abc_def0u64;
        for seed in 0u64..6 {
            let a = U4096::from_seed(seed);
            let (q, r) = a.div_rem(&U4096::from_u64(d)).unwrap();
            let fast = a.div_rem_u64(d);
            assert_eq!(q, fast.0);
            assert_eq!(r.to_u64(), Some(fast.1));
        }

        // A quotient that spans the full width: (2^4096 - 1) / 7.
        // 2^4096 = 2^(3*1365+1) = 2 mod 7, so the remainder is 1.
        let mut ones = [0u64; LIMBS_4096];
        for w in ones.iter_mut() {
            *w = u64::MAX;
        }
        let all_ones = U4096 { limbs: ones };
        let (q, r) = all_ones.div_rem(&U4096::from_u64(7)).unwrap();
        assert_eq!(r.to_u64(), Some(1));
        assert_eq!(q.bits(), 4094);
        let prod = q.mul_wide(&U4096::from_u64(7));
        assert!(prod.limbs[LIMBS_4096..].iter().all(|&w| w == 0));
        let mut low = [0u64; LIMBS_4096];
        low.copy_from_slice(&prod.limbs[..LIMBS_4096]);
        let (back, carry) = U4096 { limbs: low }.add(&r);
        assert!(!carry);
        assert_eq!(back, all_ones);
    }

    #[test]
    fn u4096_wide_reduction_spans_the_full_buffer() {
        // Every 4096-bit curve multiplication squares a 64-limb value, so the
        // dividend really does occupy all 128 slots of the Algorithm D buffer.
        // Pin that shape down: a modulus with the top bit set, squared, reduced
        // by itself, must give zero without walking past the buffer end.
        let mut p = U4096::from_seed(0x5f37_59df_9e37_79b9);
        p.set_bit(4095, true);
        if p.is_even() {
            p = p.add(&U4096::ONE).0;
        }
        assert_eq!(p.bits(), 4096, "modulus must use every limb");
        let prod = p.mul_wide(&p);
        assert_eq!(
            prod.active_words(),
            LIMBS_8192,
            "the product must fill the whole 8192-bit dividend"
        );
        assert_eq!(prod.rem_u4096(&p), U4096::ZERO);

        // Widest dividend against the widest divisor: (2^8192 - 1) is a multiple
        // of (2^4096 - 1).
        let wide = U8192 {
            limbs: [u64::MAX; LIMBS_8192],
        };
        let all_ones = U4096 {
            limbs: [u64::MAX; LIMBS_4096],
        };
        assert_eq!(wide.rem_u4096(&all_ones), U4096::ZERO);

        // A full-width modulus that is not a Mersenne number, checked against
        // the Montgomery path: same residue from both reductions.
        let mut q = U4096::from_seed(0xc0ff_ee00_9e37_79b9);
        q.set_bit(4095, true);
        q.limbs[0] |= 1;
        let ctx = MontCtx::new(&q).expect("odd full-width modulus");
        assert_eq!(ctx.n, LIMBS_4096);
        let x = prod.rem_u4096(&q);
        assert!(x < q);
        let via_mont = ctx.from_mont(&ctx.mul(&ctx.to_mont(&p), &ctx.to_mont(&p)));
        assert_eq!(x, via_mont, "full-width reductions disagree");
    }

    #[test]
    fn u4096_mod_sqrt() {
        let p = U4096::from_u64(1019);
        let a = U4096::from_u64(5);
        let root = mod_sqrt_4096(&a, &p).expect("sqrt should exist");
        assert_eq!(root.mul_mod(&root, &p), a);
    }

    #[test]
    fn u4096_large_scale_2048_bits() {
        // Construct two ~2048-bit numbers and multiply to test 4096-bit capacity
        let mut a = U4096::ZERO;
        let mut b = U4096::ZERO;
        for i in 0..32 {
            a.limbs[i] = 0xdeadbeefcafebabe;
            b.limbs[i] = 0x1234567890abcdef;
        }
        let prod = a.mul_wide(&b);
        assert!(prod.bits() > 4000);
        let m = a;
        let rem = prod.rem_u4096(&m);
        assert_eq!(rem, U4096::ZERO); // since prod = a * b, prod % a must be 0!
    }
}
