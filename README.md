# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.

---

## What it does

Given an MCMC algorithm and an energy function the checker:

1. **Discovers the state space** by BFS from a seed state: all states reachable by the algorithm from that seed are found automatically.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair using `FullSimplify` with β > 0. Both β and all energy parameters are kept symbolic. The result is an exact PASS or FAIL.

3. **Runs a numerical MCMC simulation** with genuine random numbers and compares the visited-state histogram to the Boltzmann distribution (KL divergence).

Results are printed as a table in the terminal, one row per distinct connected component found.

---

## How to run

```bash
cd ~/Desktop/DetailedBalanceChecker

# Test all bit strings of length 1–4 (both symbolic and numerical)
wolframscript -file check.wls examples2/kawasaki_1d_swap_pass.wl MaxBits=4

# 2D example needs bit strings of length 9 (3×3 grid)
wolframscript -file check.wls examples2/kawasaki_2d_pass.wl MaxBits=9

# Symbolic only (faster, no MCMC)
wolframscript -file check.wls examples2/kawasaki_1d_drift_fail.wl MaxBits=4 Mode=Symbolic

# Options
wolframscript -file check.wls examples2/kawasaki_2d_fail.wl MaxBits=9 NSteps=100000
```

**Options** (key=value, no spaces around `=`):

| Option | Default | Description |
|---|---|---|
| `MaxBits=N` | 8 | Test bit strings of length 1 through N |
| `Mode=Symbolic` | Both | Symbolic DB check only |
| `Mode=Numerical` | Both | Numerical MCMC only |
| `Mode=Both` | Both | Both checks |
| `NSteps=N` | 50000 | MCMC steps per tested system |
| `MaxBitDepth=N` | 20 | BFS tree depth cap per state |
| `Verbose=True` | False | Print per-state BFS progress |

---

## How the batch interface works

The checker generates every bit string of length 1 through `MaxBits`, ordered by length then numerically (e.g. `0`, `1`, `00`, `01`, `10`, `11`, `000`, …).

For each bit string, it calls `BitsToState` (defined in the `.wl` file) to convert the bit string to a seed state. If `BitsToState` returns `None`, the string is skipped. Otherwise:

1. BFS discovers all states reachable from the seed.
2. If the seed is already in a previously-discovered connected component, it is skipped (deduplication — avoids retesting the same system).
3. The checker runs the requested symbolic and/or numerical tests on the discovered component and prints one row.

This allows exhaustive testing of all system configurations up to a given size by simply increasing `MaxBits`.

---

## Writing an algorithm file

A `.wl` file must define exactly four things:

```mathematica
energy[state_]      := ...   (* bare energy, no beta *)
Algorithm[state_]   := ...   (* MCMC move using native random calls *)
BitsToState[bits_]  := ...   (* {0,1,...} list -> state, or None *)
numBeta             = ...    (* numeric inverse temperature *)
```

Optionally, include symbolic energy parameters:

```mathematica
symParams = <|"eps" -> {ε1, ε2, ...}, "couplings" -> {J1, ...}|>
```

### BitsToState

`BitsToState` is the only new requirement compared to a plain seed state. It must:

- Accept a list of 0s and 1s of any length.
- Return a state that `Algorithm` and `energy` can accept, **or `None`** if the bit pattern does not correspond to a valid system configuration.
- Be a **bijection**: different bit strings that return non-`None` should map to different states (or at least to states in distinct connected components).

**Occupation-based example** (3 particles on an L-site ring, sorted-occupancy state):

```mathematica
L = 4
BitsToState[bits_List] :=
  If[Length[bits] =!= L || Total[bits] =!= 3, None,
    Flatten[Position[bits, 1]]]
(* {1,1,1,0} -> {1,2,3}   {0,1,1,1} -> {2,3,4} *)
```

**Labeled-particle example** (L sites, labels assigned left to right):

```mathematica
L = 4
BitsToState[bits_List] :=
  If[Length[bits] =!= L, None,
    bits * Accumulate[bits]]
(* {1,1,1,0} -> {1,2,3,0}   {1,0,1,1} -> {1,0,2,3} *)
```

### Rules for `energy`

- Takes a state, returns a scalar (bare energy, no β).
- Must be deterministic — no random calls.
- May use symbolic (unassigned) parameters listed in `symParams`.
- β must **not** appear inside `energy`; it is injected by the checker.

### Rules for `Algorithm`

- Takes a state, returns a new state (or the same state to stay put).
- Uses `RandomInteger[]`, `RandomReal[]`, and/or `RandomChoice[]` for all randomness.
- Use `MetropolisProb[dE]` for the Metropolis acceptance probability:
  ```mathematica
  MetropolisProb[dE] = Piecewise[{{1, dE ≤ 0}, {Exp[-β dE], dE > 0}}]
  ```
  β remains symbolic during tree-building; the checker assigns it numerically for MCMC.

### Supported random calls

| Call | Behaviour during BFS |
|---|---|
| `RandomReal[]` | deferred token; fires as accept/reject when compared with `< p` |
| `RandomInteger[]` or `RandomInteger[1]` | reads one bit |
| `RandomInteger[{lo, hi}]` | reads k=⌈log₂(hi−lo+1)⌉ bits; out-of-range silently discarded |
| `RandomInteger[n]` | reads k=⌈log₂(n+1)⌉ bits |
| `RandomChoice[list]` | reads k=⌈log₂(n)⌉ bits |

Calls that cannot be intercepted (abort): `RandomVariate`, `AbsoluteTime`, `SessionTime`, `Now`.

---

## examples2/ — self-contained algorithm files

| File | System | Expected |
|---|---|---|
| `kawasaki_1d_swap_pass.wl` | L=4 ring, 3 labeled particles, swap Kawasaki, pairwise coupling | **PASS** |
| `kawasaki_1d_drift_fail.wl` | L=4 ring, 3 particles, rightward drift, flat energy | **FAIL** (symbolic), PASS (numerical) |
| `kawasaki_2d_pass.wl` | 3×3 torus, 3 particles, standard Kawasaki | **PASS** |
| `kawasaki_2d_fail.wl` | 3×3 torus, 3 particles, unbalanced particle selection | **FAIL** |

**`kawasaki_1d_swap_pass.wl`** — True Kawasaki swap: pick a random bond; if both sites are occupied, exchange the two particle labels and apply Metropolis. State = length-4 list of particle labels (0 = empty). With 3 labeled particles and pure swap moves the hole never moves, so each hole position forms a separate connected component of 6 states (3! label permutations). `BitsToState` assigns labels 1, 2, 3, … to occupied sites left to right. With `MaxBits=4`, four connected components are tested (one per hole position), all PASS.

**`kawasaki_1d_drift_fail.wl`** — Deterministic rightward current on the sorted-occupancy representation. Energy is flat (= 0), so the Boltzmann distribution is uniform and KL ≈ 0 for any algorithm (numerical PASS). The symbolic check detects the persistent current (T(i→j) = 1, T(j→i) = 0 for each cyclic pair) and reports FAIL. Demonstrates that the numerical check is necessary but not sufficient.

**`kawasaki_2d_pass.wl`** — Standard Kawasaki on a 3×3 torus with rational site energies. Pick a random particle (rejection sampling over 3), pick a random direction (N/S/E/W), apply Metropolis. With 84 reachable states the symbolic check requires a few minutes; most state pairs have T = 0 trivially. With `MaxBits=9`, the checker finds the first valid 9-bit string with 3 ones, discovers all 84 states, and tests the single connected component. PASS.

**`kawasaki_2d_fail.wl`** — Same 2D system but with a subtly broken particle selection: an unbalanced if-tree consumes 1, 2, or 3 bits per branch. When a hop changes a particle's rank in the sorted state, the forward and backward move draw from different proposal weights, breaking detailed balance. FAIL (symbolic). The violation is small enough that the numerical check may not detect it at moderate NSteps.

---

## How it works internally

The algorithm is called with native Mathematica random functions (`RandomInteger[]`, `RandomReal[]`, etc.), which are intercepted during BFS. Each call is replaced by a deterministic value from a bit tape, making the algorithm's execution path fully reproducible.

- `RandomInteger[]` or `RandomInteger[1]` — reads one bit; contributes ½ to path weight.
- `RandomReal[] < p` — reads one bit; contributes p (accept) or 1−p (reject) to path weight.
- T(i→j) = **sum of path weights** over all bit sequences taking state i to state j.

---

## Practical limits

- **States**: discovered automatically; keep to ≤ 15–20 for fast results. The 2D example with 84 states is feasible (most pairs are non-adjacent and resolve immediately).
- **Bit depth**: `MaxBitDepth=20` covers most algorithms.
- **FullSimplify**: the most expensive step. With symbolic energy parameters, PiecewiseExpand resolves Metropolis conditions branch-by-branch using exact algebra.
- **Flat-energy caveat**: if all energies are equal the Boltzmann distribution is uniform, so KL ≈ 0 for any algorithm — even non-reversible ones. The symbolic DB check catches this; the numerical check is uninformative on flat landscapes.

---

## Files

```
dbc_core.wl                   Core library
check.wls                     Terminal interface (batch bit-string mode)
show_report.py                Python/matplotlib report renderer (legacy)
examples2/
  kawasaki_1d_swap_pass.wl    L=4 ring, 3 labeled particles, swap     [PASS]
  kawasaki_1d_drift_fail.wl   L=4 ring, 3 particles, rightward drift  [FAIL / num-PASS]
  kawasaki_2d_pass.wl         3x3 torus, 3 particles, Kawasaki        [PASS]
  kawasaki_2d_fail.wl         3x3 torus, 3 particles, unbalanced tree [FAIL]
examples/                     Legacy examples (old format, unsupported)
```
