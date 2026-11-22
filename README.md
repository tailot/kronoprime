# KronoPrime

A **Rust** engine for analyzing and factoring large integers (`u64` standard profile, `u128` extended profile), integrating the geometric and modular constraints of the **per-decade 3x3 matrix** grid. The cubic and elliptic-curve layer additionally runs on a dedicated **4096-bit profile** (`U4096`: 64 `u64` limbs, stack-allocated, still zero external dependencies), reached automatically as soon as a value is something `u128` cannot read.

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

**Beyond Pollard rho.** Pollard rho (heuristic `O(p_min^(1/4))`) is excellent up to factors of ~10-12 digits. At 4096 bits `cubic factor-probe` already runs a four-stage attack on a composite modulus - trial division, discriminant gcd, **Lenstra ECM stage 1** (a curve through a chosen point, a smooth multiple of its order via the x-only ladder, the factor read off the non-invertible denominator), then a Pollard-Brent walk sweeping the curve constant - and a 157-bit semiprime falls to it inside the test suite. What is still roadmap: ECM **stage 2** at width, **SIQS** (balanced semiprimes ~80-100 digits), **NFS** (asymptotically the best, out of scope due to implementation cost).

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

## The 4096-bit wide profile

Curve work is not capped at `u128`. `src/u4096.rs` is a fixed-width 4096-bit unsigned integer - 64 `u64` limbs, stack-allocated, no external crates - with Knuth's Algorithm D over a `[u64; 129]` stack buffer, Montgomery multiplication, and full modular arithmetic (`add_mod`, `sub_mod`, `mul_mod`, `div_rem`, `inv_mod`, `gcd`, `pow_mod`, bit accessors, decimal/hex I/O). `src/cubic4096.rs` rebuilds the whole cubic layer on top of it: discriminant `-16(4a^3 + 27b^2)`, classification (nodal, cuspidal, CM supersingular with a closed-form order, Montgomery form, ordinary), affine and Jacobian point arithmetic with Lenstra-style factor extraction, scalar multiplication, a Montgomery x-only ladder, CM curve search and curve-order certification. `3x3` matrix invariants and the factor probe have wide counterparts as well, and Miller-Rabin uses the 40-base set at this width.

Routing is decided per token, not by a flag: any `--modulus`, `--a`, `--b`, `--point`, `--scalar` or `--matrix` cell that the `u128` parser **cannot read** but `U4096` can - a decimal beyond `u128::MAX`, or any `0x`-prefixed hexadecimal - switches that subcommand to the wide engine. Same flags, same table/JSON shapes, same exit codes.

Measured on a desktop release build, by the two load tests (`cargo test --release -- --ignored --nocapture`):

| Operation (4096-bit modulus) | Cost |
|---|---|
| one `mul_mod` (8192-bit product + one Algorithm D reduction) | 4.6 microseconds |
| `[k]P` with a 192-bit scalar | 6.0 ms |
| Miller-Rabin, 4 bases, composite caught early | 0.1 ms |
| Miller-Rabin, 40 bases, on the prime the search finds | 2.14 s |
| finding a random 4096-bit prime (3755 odd candidates) | 24.6 s |
| order certification `[p+1]P = O`, per sampled point | ~135 ms |

The certification test does exactly what it says: it finds a 4096-bit prime `p`, picks the CM family the residue of `p` supports (`j = 0` when `p = 2 mod 3`, `j = 1728` when `p = 3 mod 4`), classifies the curve, and checks that its order really is `p + 1` by annihilating sampled points. Both tests are `#[ignore]`d so the ordinary suite stays fast.

## Usage: CLI and Cargo

```sh
cargo build --release
./target/release/kronoprime 561 100000000520000000627
cat numbers.txt | ./target/release/kronoprime --matrix
./target/release/kronoprime --json 104729
./target/release/kronoprime --quadruplets 5

# Cubic & Elliptic Curves subcommands:
./target/release/kronoprime cubic scan --modulus 1019 --limit 3 --class cm --format table
./target/release/kronoprime cubic from-matrix --matrix "[12,5,7; 3,19,2; 8,1,30]" --modulus 1009 --format json
./target/release/kronoprime cubic verify --a 0 --b 5 --modulus 1019 --points 3 --format table
./target/release/kronoprime cubic factor-probe --modulus 2021027 --limit 50 --format json
./target/release/kronoprime cubic point-ops --a 0 --b 5 --modulus 1019 --point "12,431" --scalar 1020 --format table

# The same subcommands at width: a value u128 cannot read switches engine
./target/release/kronoprime cubic verify --a 0 --b 5 --modulus 0x3FB --points 2 --format table
./target/release/kronoprime cubic factor-probe --modulus 1000000000000000000000000000000000000007 --limit 4 --format json
./target/release/kronoprime cubic point-ops --a 0 --b 4 --modulus 0xFFFFFFFFFFFFFFFF --point "0,2" --scalar 3 --format table
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

Flags: `--json`, `--matrix`, `--quadruplets <K>` (or `--quadruplets=<K>`, max 100000), `-h/--help`, `--version`. Valid inputs: integers in `0..=u128::MAX`; non-numeric or out-of-range inputs are reported on stderr (exit code 1). The `cubic` subcommands use the same flags and accept values up to 4096 bits, decimal or `0x` hexadecimal; their exit codes are width-independent (`2` when `factor-probe` finds a factor, `1` when a `verify` certification fails, `0` otherwise).

## Six worked examples

Run from Linux as `./target/release/kronoprime` (on Windows the same commands carry the `.exe` suffix). Every output below is a real capture, followed by what it actually says.

### 1. The classic Carmichael trap

```sh
./target/release/kronoprime 561
```

```text
Input: 561
Status: COMPOSITE
Matrix N: 57 | Slot: +1
Factors: [3, 11, 17]
Elapsed: 0.01 ms
```

561 is the smallest Carmichael number: it *looks* prime to any Fermat-style probabilistic test, because `b^560 = 1 (mod 561)` for every `b` coprime to it - it passes many probable-prime tests yet is composite. Deterministic Miller-Rabin is exactly what cannot be fooled: `Status: COMPOSITE`, no maybe about it. The placement `Matrix N: 57 | Slot: +1` puts it in the decade 561..569, first active cell, and stage 1 trial division (3 divides it) hands back the full factorization `3 * 11 * 17` in a hundredth of a millisecond.

### 2. The 3x3 grid in high definition

```sh
./target/release/kronoprime --matrix 104729
```

```text
Input: 104729
Status: PRIME
Matrix N: 10473 | Slot: +9
Factors: [104729]
Elapsed: 0.00 ms
3x3 grid for Matrix N=10473 (104721..104729):
104721  104722  104723* 
104724  104725  104726  
104727  104728  104729* 
Active slots: +1=104721 composite +3=104723 prime +7=104727 composite +9=104729 prime (input)
```

104729 is the 10,000th prime, and `--matrix` prints the whole decade it lives in: the cells `104721..104729` laid out as the actual 3x3 grid, with prime cells starred. So the grid answers "where does a prime land?" (row 0, column 2, slot `+9`) and "which cells of this decade are alive?" at a glance: only two of the four active slots are prime here (`+3` = 104723 and `+9` = 104729), while `+1` (divisible by 3) and `+7` are dead. The `(input)` tag marks the queried cell.

### 3. Hunting the first 5 prime quadruplets

```sh
./target/release/kronoprime --quadruplets 5
```

```text
Quadruplet: Matrix N=2 -> [11, 13, 17, 19]
Quadruplet: Matrix N=11 -> [101, 103, 107, 109]
Quadruplet: Matrix N=20 -> [191, 193, 197, 199]
Quadruplet: Matrix N=83 -> [821, 823, 827, 829]
Quadruplet: Matrix N=149 -> [1481, 1483, 1487, 1489]
Searched with the mod-3003 wheel in 0.01 ms
```

A quadruplet `(p, p+2, p+6, p+8)` is exactly a decade whose four active slots `{+1, +3, +7, +9}` are all prime, and the five lines are the first five such decades: N=2 (the only one where the wheel primes 11 and 13 fall inside the window), then N=11, 20, 83, 149. The closing line is the point of the whole grid design: the mod-3003 wheel pre-filter skipped ~93.7% of the decades before any primality test ran, so the hunt cost 0.01 ms.

### 4. The semiprime in ~80 ms: two large primes meet

```sh
./target/release/kronoprime 100000000520000000627
```

```text
Input: 100000000520000000627
Status: COMPOSITE (large prime factors)
Matrix N: 10000000052000000063 | Slot: +7
Factors: [10000000019, 10000000033]
Elapsed: 25.20 ms
```

A 21-digit semiprime: neither factor (both 11-digit primes) is reachable by the trial-division cap, so Miller-Rabin certifies compositeness - `COMPOSITE (large prime factors)` - and Pollard-Brent does the splitting in tens of milliseconds (machine-dependent; the table above measured ~80 ms on a Windows desktop, this capture 25 ms). The slot line even self-validates: the factors sit in slots +9 and +3, and `9 * 3 = 27` ends in 7 - exactly the reported `Slot: +7`, the mod-10 product rule checked for free.

### 5. "Cracking" a composite on the fly with the discriminant GCD

```sh
./target/release/kronoprime cubic factor-probe --modulus 2021027 --limit 50 --format json
```

```json
{
  "status": "FactorFound",
  "modulus": 2021027,
  "factor1": 1009,
  "factor2": 2003,
  "source_curve": { "a": 200300, "b": 200300 },
  "mechanism": "DiscriminantGcd"
}
```

Here the composite is fed to the *curve* layer as a modulus. Before any heavy machinery, the probe computes the discriminant of the constructed cubic `y^2 = x^3 + ax + b` modulo `n` and takes a `gcd` against `n`: one of the seeds hits a non-invertible denominator and the gcd leaks a factor. The JSON says which curve did it - `source_curve: (a, b) = (200300, 200300)`, mechanism `DiscriminantGcd` - and the modulus splits: `2021027 = 1009 * 2003`. Exit code `2` signals "factor found", so scripts can branch on it without parsing the JSON.

### 6. Automatic lane switch to the 4096-bit profile

```sh
./target/release/kronoprime cubic verify --a 0 --b 5 --modulus 0x3FB --points 2 --format table
```

```text
Curve order verification: y^2 = x^3 + 5 mod 1019 (class: CmSupersingular I)
Expected theoretical order: 1020
+---------+--------------+--------------+--------+
| Point # | X coordinate | Y coordinate | Result |
+---------+--------------+--------------+--------+
| 1       | 0            | 987          | PASSED |
| 2       | 0            | 32           | PASSED |
+---------+--------------+--------------+--------+
Overall result: CERTIFIED (2/2 points satisfy [1020]P = O)
```

Nothing here is big - `0x3FB` is just 1019 in hexadecimal - but the `0x` prefix is a token the `u128` parser refuses, so the routing rule kicks in automatically and the *entire* subcommand runs on the 4096-bit engine (`U4096` + `cubic4096`): no flag asked, same table shape, same exit codes. The math: since `1019 = 2 (mod 3)`, `y^2 = x^3 + 5` is a CM supersingular curve (`j = 0` family) with expected order `p + 1 = 1020`, and certification means annihilating sampled points: `[1020]P = O`. The two sampled points share `x = 0` (then `y^2 = 5` has the conjugate roots 32 and `987 = 1019 - 32`), both `PASSED`, so `Overall result: CERTIFIED` and exit code 0.

## Code structure and quality

```text
src/
  main.rs        entry point
  cli.rs         argument parsing, stdin, text/JSON reports, cubic subcommands, wide-profile routing
  matrix.rs      3x3 grid, slots, mod-10 compatibility, mod-3003 wheel, matrix invariants & Tschirnhaus (u128 and 4096)
  primality.rs   Miller-Rabin + overflow-safe modular arithmetic (u128), 40-base Miller-Rabin for U4096
  factor.rs      trial division, Pollard-Brent, worklist, cubic factor probe (u128 and 4096)
  cubic.rs       cubic & elliptic curve discovery, classification, point operations, Tonelli-Shanks
  u4096.rs       4096-bit fixed-width integer: Knuth D on a stack buffer, Montgomery, full modular arithmetic
  cubic4096.rs   the same cubic/curve layer at 4096 bits: discriminant, classes, points, ladder, CM search
tests/
  integration.rs end-to-end invariants (factor product, slots, timing, cubic & matrix validation)
```

- **No panics on extreme inputs**: overflow-free arithmetic (`add_mod`/`mul_mod` with a native `u128` path below 2^64 and multiply-and-add above it), edge cases covered (`0`, `1`, `u64::MAX`, Mersenne numbers, Carmichael numbers, prime powers).
- **Tests**: `cargo test` and `cargo test --release` (unit + integration); `cargo test --release -- --ignored --nocapture` runs the two 4096-bit load and certification tests (~28 s) and prints the timings above.
- **Debug-build contracts**: at width, `add_mod`/`sub_mod`/`mul_mod` assert that both operands are already reduced and `div_knuth` asserts the buffer bound, so the debug suite is a real check rather than a slower copy of the release one. Every constant the formulas multiply by (2, 3, 4, 16, 27) is taken into the residue ring first, which is what keeps moduli as small as 2 honest.
- **Quality**: `cargo fmt --all -- --check` and `cargo clippy --all-targets -- -D warnings` clean; precomputed tables in `const fn`; zero external dependencies - the 4096-bit path stays on the stack, so the handful of lints that would push it onto the heap are allowed with a comment instead of silenced.

## CI/CD (GitHub Actions, Linux)

`.github/workflows/ci.yml` on `ubuntu-latest`:

1. `cargo fmt --all -- --check` + `cargo clippy --all-targets -- -D warnings`;
2. debug and release tests (`cargo test --all-targets`, `cargo test --release --all-targets`);
3. release build for `x86_64-unknown-linux-gnu` with the `kronoprime` artifact (`actions/upload-artifact`);
4. on `v*` tags: the Linux binary is published to GitHub Releases.
