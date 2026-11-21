# KronoPrime

A **Rust** engine for analyzing and factoring large integers (`u64` standard profile, `u128` extended profile), integrating the geometric and modular constraints of the **per-decade 3x3 matrix** grid.

## Project overview

For every input integer, `kronoprime`:

1. decides **primality** (deterministic Miller-Rabin over the whole `u64` profile);
2. computes its **placement in the 3x3 decade matrix** (matrix index `N`, slot, row, column);
3. produces its sorted **prime factorization**, without trial division up to sqrt(n).

The architecture is staged, from cheapest to most powerful:

| Stage | Algorithm | Role |
|---|---|---|
| 1 | Small factors (2, 3, 5) + trial division with a const table of primes up to **1009** | instant extraction of easy factors, with a tight cap so no time is wasted on numbers with only large factors |
| 2 | **Miller-Rabin**, deterministic (short 3/5-base sets below `4.76e9`/`2.15e12`, compact 7-base set for `n < 2^64`; 40 fixed bases for the extended `u128` profile) | certifies remaining cofactors with no false positives on the standard profile, doing the least work possible at each magnitude |
| 3 | **Pollard's rho with Brent's cycle detection** (`f(x) = x^2 + c`, `gcd` accumulated over 128 differences, overflow-free `u128` arithmetic) | large factors up to ~10-12 digits on `u64` inputs in short time; multiple `(c, y0)` pairs with fallback on degenerate cycles |
| 4 | Iterative worklist (explicit stack) | decomposes every non-trivial divisor without stack-overflow risk |

> **Terminology note.** A large *prime* is already indivisible: factorization just returns the number itself after the primality proof. The real computational problem is factoring large **composites** (semiprimes with 12-20+ digits) and certifying primality quickly.

## The 3x3 matrices: visual example

Every decade `N >= 1` (base `B = 10(N - 1)`) is a 3x3 grid covering the range `[B+1, B+9]`. The **active slots** - the final digits `{1, 3, 7, 9}` that can host primes > 5 - are `+1`, `+3`, `+7`, `+9`.

```text
Matrix N=1 (1-9)                Matrix N=2 (11-19)
 +---+---+---+                  +----+----+----+
 | 1 | 2 | 3 |                  | 11 | 12 | 13 |
 +---+---+---+                  +----+----+----+
 | 4 | 5 | 6 |                  | 14 | 15 | 16 |
 +---+---+---+                  +----+----+----+
 | 7 | 8 | 9 |                  | 17 | 18 | 19 |
 +---+---+---+                  +----+----+----+

active slots N=1:  +1=1(-)  +3=3(P)  +7=7(P)  +9=9(-)
active slots N=2:  +1=11(P) +3=13(P) +7=17(P) +9=19(P)   <- first quadruplet
```

Multiples of 10 sit on the border between two decades and occupy no cell.

## Mathematical foundations

**Index and slot.** For `x >= 1`:

- `slot(x) = ((x - 1) mod 10) + 1` in `{1..9}` (it evaluates to 10, on the border, for multiples of 10);
- `N(x) = floor((x - 1)/10) + 1`;
- cell `(r, c)` of matrix `N`: `v = 10(N - 1) + 3r + c + 1`.

**Slot compatibility in products.** For primes `p, q > 5` we have `s_m = s_p * s_q (mod 10)`: slot `3` arises from the pairs `(1,3)` and `(7,9)`, slot `7` from `(1,7)` and `(3,9)`; slots `1` and `9` allow every combination of `{1,3,7,9}`. On a single factor this rule never excludes a final digit (all slots are invertible mod 10): it is a diagnostic filter, not a shortcut.

**The mod-3003 wheel for quadruplets.** To search for decades whose four active slots are all prime, the congruences of `N - 1` modulo `3, 7, 11, 13` (`3003 = 3*7*11*13`) discard in bulk the matrices with at least one trivially composite slot (with the exception that a slot may coincide with the wheel prime itself, like 11 and 13 in decade N=2): only **190 classes out of 3003** survive (~6.3%), before any primality test runs.

**The real role of `N` in factorization.** `N(n)` only identifies the decade of `n` (magnitude ~`n/10`): it reveals nothing about the indices of the factors and does not reduce the complexity, which depends on the smallest prime factor `p_min` and is handled with Miller-Rabin + Pollard rho. The index `N` does give a concrete structural advantage in the quadruplet search (above).

**Beyond Pollard rho (roadmap).** Pollard rho (heuristic `O(p_min^(1/4))`) is excellent up to factors of ~10-12 digits. In the future, beyond KronoPrime's lightweight `u64`/`u128` core: **ECM** (functional in the smallest factor), **SIQS** (balanced semiprimes ~80-100 digits), **NFS** (asymptotically the best, out of scope due to implementation cost).

## Why the 3x3 grid pays off — reasons and measurements

The grid is not decoration: it is a re-indexing of the integers that turns
"is this value composite?" into cheap statements about the **index** `N`.
Every concrete reason, followed by measured numbers.

1. **O(1) placement.** `N = floor((x-1)/10)+1` and `slot = ((x-1) mod 10)+1`
   are one division by a constant (compiler-issued multiply-and-shift).
   Placement of any value, any magnitude, is free and bijective over the
   non-multiples of 10.
2. **A free slot filter.** Only the final digits `{1, 3, 7, 9}` can host
   primes `> 5`, so 5 of the 9 cells die with zero arithmetic — before any
   primality test. This alone is the 3.9x in the table below.
3. **Divisibility becomes a congruence on the index.** A slot value is
   `10(N-1) + o`, so for every prime `p` with `gcd(p, 10) = 1`:
   `p | 10(N-1)+o  <=>  N = 1 - o*10^-1 (mod p)`. One check on `N` kills
   that slot **simultaneously for every decade**, not just the current one.
4. **Periodicity compresses the sieve.** Those congruences for
   `p in {3, 7, 11, 13}` are periodic with period `3003 = 3*7*11*13`: only
   **190 admissible classes out of 3003** (6.3%) can host all-prime slots.
   The whole wheel is therefore a **const 190-entry table** built at
   compile time and resident in L1 — no sieve arrays, so the search works
   even near the `u128` range where no physical sieve could fit in memory.
5. **Amortization.** One congruence on `N` answers for the 4 active slots
   at once; one wheel visit prunes one entire decade (9 values).
6. **Jump instead of test.** Since survivors are known table entries, the
   scanner walks directly from one admissible residue to the next: 190
   visits per 3003 decades, ~93.7% never even reach a filter.
7. **Patterns become indices.** A prime quadruplet `(p, p+2, p+6, p+8)`
   is exactly a decade with all four active slots prime, so a pattern
   search becomes an index iteration with table-driven skips; the first
   ones (N = 2, 11, 20, 83, 149) are compact memory anchors.
8. **Last-digit coherence of products.** `s(p*q) = s(p)*s(q) (mod 10)`
   constrains candidate factor pairs (balanced searches) and validates
   reported factorizations for free.
9. **It stays exact.** The wheel is only a pre-filter: survivors are
   certified by deterministic Miller-Rabin, and the boundary case where a
   slot value equals a wheel prime itself (11 and 13 in N=2) is handled
   exactly, with a `const` assert pinning the class count to 190.

**Measured** (Intel Core Ultra 9 285K, Windows, `cargo build --release`
with LTO; all runs find the same 20,000 quadruplets, the last at
`N = 64,812,215`, sharing the same primality engine):

| Quadruplet-search strategy | Time | Speedup |
|---|---:|---:|
| No grid: scan every odd `x`, test `(x, x+2, x+6, x+8)` | 52,799 ms | 1.0x |
| Decade index, slot filter only (no wheel) | 13,553 ms | 3.9x |
| Decade index + per-decade mod-3003 wheel test | 2,423 ms | 21.8x |
| Decade index + direct residue jumps (**current**) | 2,189 ms | **24.1x** |

Inside the current strategy, the adaptive Miller-Rabin base sets add a
further large cut (same search, same machine):

| Miller-Rabin bases | Time |
|---|---:|
| Fixed 7-base `u64` set at every magnitude | 3,974 ms |
| Adaptive 3/5/7-base sets (**current**) | 2,189 ms (**-45%**) |

**Where it does not pay off — stated honestly.** For one *given* composite
`n`, knowing `N(n)` speeds up nothing: `p | n` is already one remainder,
`p | 10(N-1)+o` is the same computation reworded, and no information about
the smallest factor `p_min` leaks from the decade. Factoring
`100000000520000000627` takes ~80 ms with or without the grid. The grid
earns time only when **generating or discarding candidates in bulk**
(quadruplet searches, tables, range sieves), which is exactly where it is
used in this engine.

## Usage: CLI and Cargo

```sh
cargo build --release
./target/release/kronoprime 561 100000000520000000627
cat numbers.txt | ./target/release/kronoprime --matrix
./target/release/kronoprime --json 104729
./target/release/kronoprime --quadruplets 5
```

Standard output:

```text
Input: 100000000520000000627
Status: COMPOSITE (large prime factors)
Matrix N: 10000000052000000063 | Slot: +7
Factors: [10000000019, 10000000033]
Elapsed: 81.02 ms
```

With `--matrix` the number's 3x3 grid is also shown (prime cells marked with `*`) plus an active-slot summary. With `--json` every number becomes an object with `input`, `status`, `prime`, `matrix_n`, `slot`, `factors`, `elapsed_ms` (big integers as strings, no precision loss).

`--quadruplets 5` example:

```text
Quadruplet: Matrix N=2 -> [11, 13, 17, 19]
Quadruplet: Matrix N=11 -> [101, 103, 107, 109]
Quadruplet: Matrix N=20 -> [191, 193, 197, 199]
Quadruplet: Matrix N=83 -> [821, 823, 827, 829]
Quadruplet: Matrix N=149 -> [1481, 1483, 1487, 1489]
```

Flags: `--json`, `--matrix`, `--quadruplets <K>` (or `--quadruplets=<K>`, max 100000), `-h/--help`, `--version`. Valid inputs: integers in `0..=u128::MAX`; non-numeric or out-of-range inputs are reported on stderr (exit code 1).

## Code structure and quality

```text
src/
  main.rs        entry point
  cli.rs         argument parsing, stdin, text/JSON reports
  matrix.rs      3x3 grid, slots, mod-10 compatibility, mod-3003 wheel
  primality.rs   Miller-Rabin + overflow-safe modular arithmetic (u128)
  factor.rs      trial division, Pollard-Brent, worklist
tests/
  integration.rs end-to-end invariants (factor product, slots, timing)
```

- **No panics on extreme inputs**: overflow-free arithmetic (`add_mod`/`mul_mod` with a native `u128` path below 2^64 and multiply-and-add above it), edge cases covered (`0`, `1`, `u64::MAX`, Mersenne numbers, Carmichael numbers, prime powers).
- **Tests**: `cargo test` and `cargo test --release` (unit + integration).
- **Quality**: `cargo fmt --all -- --check` and `cargo clippy --all-targets -- -D warnings` clean; precomputed tables in `const fn`; zero external dependencies.

## CI/CD (GitHub Actions, Linux)

`.github/workflows/ci.yml` on `ubuntu-latest`:

1. `cargo fmt --all -- --check` + `cargo clippy --all-targets -- -D warnings`;
2. debug and release tests (`cargo test --all-targets`, `cargo test --release --all-targets`);
3. release build for `x86_64-unknown-linux-gnu` with the `kronoprime` artifact (`actions/upload-artifact`);
4. on `v*` tags: the Linux binary is published to GitHub Releases.
