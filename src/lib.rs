//! KronoPrime — analysis and factorization engine for large integers.
//!
//! The crate is split into small, independent modules:
//!
//! * [`matrix`] — the 3x3 decade grid: matrix index `N`, slots, slot
//!   compatibility modulo 10 and the modular wheel `mod 3003`.
//! * [`primality`] — deterministic Miller-Rabin (whole `u64` profile) plus a
//!   fixed-base strong probable prime test for the extended `u128` profile.
//! * [`factor`] — factorization pipeline: small-prime trial division,
//!   primality screening and Pollard's rho with Brent's cycle detection.
//! * [`cli`] — command line interface (`--json`, `--matrix`,
//!   `--quadruplets <K>`, arguments and stdin).
//!
//! All arithmetic is overflow safe: `u64` moduli use native `u128`
//! multiplications, `u128` moduli use overflow-free multiply-and-add.

pub mod cli;
pub mod factor;
pub mod matrix;
pub mod primality;
