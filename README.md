# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.
Produces a full diagnostic PNG report for each algorithm tested.

---

## What it does

Given an MCMC algorithm and an energy function the checker:

1. **Builds the decision tree** by exhaustive BFS over all random-call sequences the algorithm can make. The state space is discovered automatically from a single seed state.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair using `FullSimplify` with β > 0. Both β and all energy parameters are kept symbolic throughout. The result is an exact PASS or FAIL.

3. **Runs a numerical MCMC simulation** with genuine random numbers and compares the visited-state histogram to the Boltzmann distribution (KL divergence).

4. **Renders a PNG report** showing decision trees, symbolic transition matrix, detailed-balance pair table, and simulated vs Boltzmann scatter plot.

---

## How it works internally

The algorithm is called with native Mathematica random functions (`RandomInteger[]`, `RandomReal[]`, etc.), which are intercepted during BFS. Each call is replaced by a deterministic value from a bit tape, making the algorithm's execution path fully reproducible.

- `RandomInteger[]` or `RandomInteger[1]` — reads one bit; contributes ½ to path weight.
- `RandomReal[] < p` — reads one bit; contributes p (accept) or 1−p (reject) to path weight.
- T(i→j) = **sum of path weights** over all bit sequences taking state i to state j.

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

## Writing an algorithm file

An algorithm file is a plain `.wl` file containing exactly four things:

```mathematica
energy[state_] := ...      (* bare energy, no beta *)

Algorithm[state_] := ...   (* MCMC move using native random calls *)

seedState = ...             (* any valid starting state *)
numBeta   = ...             (* inverse temperature, numeric *)
```

Optionally, include symbolic energy parameters:

```mathematica
symParams = <|"eps" -> {ε1, ε2, ...}, "couplings" -> {J1, ...}|>
```

When `symParams` is present, the checker keeps site energies and coupling strengths as unassigned symbols throughout the symbolic check, then assigns random numerical values for the MCMC run. This allows `FullSimplify` to verify detailed balance across **all** possible energy values simultaneously.

### Rules for `energy`

- Takes a state, returns a scalar (the bare energy, no β).
- Must be deterministic — no random calls.
- May use symbolic (unassigned) parameters listed in `symParams`.
- β must NOT appear inside `energy`; it is injected by the checker via `Block`.

### Rules for `Algorithm`

- Takes a state, returns a new state (or the same state to stay put).
- Uses `RandomInteger[]`, `RandomReal[]`, and/or `RandomChoice[]` for all randomness.
- Use `MetropolisProb[dE]` for the Metropolis acceptance probability:
  ```mathematica
  MetropolisProb[dE] = Piecewise[{{1, dE ≤ 0}, {Exp[-β dE], dE > 0}}]
  ```
  β remains symbolic during tree-building; the checker assigns it numerically for MCMC.
- Custom acceptance criteria work as long as β stays symbolic:
  ```mathematica
  BarkerProb[dE_] := 1 / (1 + Exp[β * dE])
  ```

### Minimal example

```mathematica
(* L=4 ring, 1 particle, Kawasaki hop + Metropolis *)
L = 4
eps = {\[Epsilon]1, \[Epsilon]2, \[Epsilon]3, \[Epsilon]4}
symParams = <|"eps" -> eps|>

energy[s_Integer] := eps[[s]]

Algorithm[s_Integer] := Module[{dir, nbr, dE},
  dir = RandomInteger[];
  nbr = Mod[s + If[dir == 0, -1, 1] - 1, L] + 1;
  dE  = energy[nbr] - energy[s];
  If[RandomReal[] < MetropolisProb[dE], nbr, s]
]

seedState = 1
numBeta   = 1
```

Run with:
```bash
wolframscript -file check.wls my_algorithm.wl
```

---

## Running checks

### Single algorithm

```bash
cd ~/Desktop/DetailedBalanceChecker

# New mode (recommended): file defines Algorithm, energy, seedState, numBeta
wolframscript -file check.wls examples2/kawasaki_1d_swap_pass.wl
wolframscript -file check.wls examples2/kawasaki_1d_drift_fail.wl
wolframscript -file check.wls examples2/kawasaki_2d_pass.wl
wolframscript -file check.wls examples2/kawasaki_2d_fail.wl

# Options (key=value, no spaces around =)
wolframscript -file check.wls my_alg.wl NSteps=200000 Verbose=False
wolframscript -file check.wls my_alg.wl SystemName="My system" OpenWindow=False
```

### Legacy mode (existing examples)

```bash
wolframscript -file check.wls examples/ring_kawasaki.wl 1 KawasakiRing energy$rk 1
wolframscript -file check.wls examples/multi_particle.wl "{1,2}" KawasakiMulti energy$mp 1
```

### Batch scripts (existing examples)

```bash
wolframscript -file run_checks.wls          # core ring Kawasaki examples
wolframscript -file run_variable_bit.wls    # variable-depth examples
wolframscript -file run_new_api.wls         # Barker + asymmetric proposal
wolframscript -file run_extended.wls        # extended feature examples
```

---

## examples2 — self-contained algorithm files

| File | System | Expected |
|---|---|---|
| `kawasaki_1d_swap_pass.wl` | L=4 ring, 3 particles, bond-swap Kawasaki + NN coupling | **PASS** |
| `kawasaki_1d_drift_fail.wl` | L=4 ring, 3 particles, rightward drift, flat energy | **FAIL** (symbolic), PASS (numerical) |
| `kawasaki_2d_pass.wl` | 3×3 torus, 3 particles, standard Kawasaki | **PASS** |
| `kawasaki_2d_fail.wl` | 3×3 torus, 3 particles, unbalanced particle selection | **FAIL** |

### What each example demonstrates

**`kawasaki_1d_swap_pass.wl`** — Bond-swap Kawasaki (pick a random adjacent bond and exchange its contents) on a 4-site ring with 3 particles and nearest-neighbour coupling. Uses fully symbolic site energies and coupling constant. Demonstrates the standard pass case with a different proposal mechanism from the 1-particle hop examples.

**`kawasaki_1d_drift_fail.wl`** — Deterministic rightward current: the particle immediately left of the hole always hops right. Energy is flat, so the Boltzmann distribution is uniform and the MCMC histogram matches it perfectly (KL ≈ 0). Yet the symbolic check detects the persistent current (T(i→j) = 1, T(j→i) = 0 for each cyclic pair). Demonstrates that the numerical check is necessary but not sufficient.

**`kawasaki_2d_pass.wl`** — Standard Kawasaki on a 3×3 torus with rational site energies. Pick a random particle (rejection sampling over 3), pick a random direction (N/S/E/W), apply Metropolis. With 84 states the symbolic check requires several minutes; most state pairs have T = 0 trivially, so only ~300 adjacent pairs need `FullSimplify`. Includes `PlotState` and `PlotTransition` helper functions for visualising the 2D lattice from a Mathematica notebook.

**`kawasaki_2d_fail.wl`** — Same 2D system but with a subtly broken particle selection. Instead of the correct uniform `RandomInteger[{0,2}]`, the algorithm uses a three-level if-tree where rank-1, rank-2, and rank-3 particles consume 1, 2, and 3 bits respectively. This produces unequal proposal probabilities (1/2, 1/4, 1/4 for the three branches — but the rank-3 branch has an extra unused bit that doubles its effective weight relative to rank-2). Whenever a hop changes a particle's rank in the sorted state the forward and backward move draw from different proposal weights, breaking detailed balance. The violation is subtle: all particles can be proposed, direction selection is correct, Metropolis is correct — only the bit-depth imbalance in the selection tree reveals the flaw.

### Visualising the 2D lattice

After running the 2D examples, call the helper functions from a Mathematica notebook:

```mathematica
Get["examples2/kawasaki_2d_pass.wl"]

(* Print current state as ASCII grid *)
PlotState[{1, 5, 9}]

(* Show a transition *)
PlotTransition[{1, 5, 9}, {2, 5, 9}]
```

Grid layout (site numbering):
```
1 2 3
4 5 6
7 8 9
```

---

## Practical limits

- **States**: discovered automatically; keep to ≤ 10–15 for fast results. The 2D example with 84 states is feasible (most pairs are non-adjacent and resolve immediately).
- **Bit depth**: `MaxBitDepth=20` covers most algorithms.
- **FullSimplify**: the most expensive step. With symbolic energy parameters, PiecewiseExpand resolves Metropolis conditions branch-by-branch using exact algebra. For the 2D examples (numeric energies) this is faster.
- **Flat-energy caveat**: if all energies are equal the Boltzmann distribution is uniform, so KL ≈ 0 for any algorithm — even non-reversible ones. The symbolic DB check catches this; the numerical check is uninformative on flat landscapes by design.

---

## Files

```
dbc_core.wl                   Core library
show_report.py                Python/matplotlib report renderer
check.wls                     Terminal interface (new + legacy modes)
run_checks.wls                Core ring Kawasaki examples (3)
run_variable_bit.wls          Variable-depth examples (2)
run_new_api.wls               Barker + asymmetric-proposal examples (4)
run_extended.wls              Extended feature examples (12, A-L)
examples/                     Original examples (legacy format)
examples2/
  kawasaki_1d_swap_pass.wl    L=4 ring, 3 particles, bond-swap Kawasaki   [PASS]
  kawasaki_1d_drift_fail.wl   L=4 ring, 3 particles, rightward drift      [FAIL / num-PASS]
  kawasaki_2d_pass.wl         3x3 torus, 3 particles, Kawasaki            [PASS]
  kawasaki_2d_fail.wl         3x3 torus, 3 particles, unbalanced tree     [FAIL]
```
