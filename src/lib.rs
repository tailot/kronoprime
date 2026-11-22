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
//! * [`cubic`] — cubic and elliptic curve discovery, classification,
//!   point arithmetic and modular factorization probing.
//! * [`cubic4096`] — the same curve theory lifted to 4096-bit moduli, with
//!   Jacobian scalar multiplication so a wide walk does not pay an inversion
//!   per step.
//! * [`u4096`] — the stack-allocated 4096-bit unsigned integer: Knuth
//!   division, Montgomery reduction, modular inversion and a wide
//!   Miller-Rabin companion in [`primality`].
//! * [`cli`] — command line interface (`--json`, `--matrix`,
//!   `--quadruplets <K>`, arguments, stdin and `cubic` subcommands).
//!
//! All arithmetic is overflow safe: `u64` moduli use native `u128`
//! multiplications, `u128` moduli use overflow-free multiply-and-add, and the
//! wide profile keeps every intermediate inside its own limb array.

pub mod cli;
pub mod cubic;
pub mod cubic4096;
pub mod factor;
pub mod matrix;
pub mod primality;
pub mod u4096;
