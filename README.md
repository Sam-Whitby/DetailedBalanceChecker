# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.
Produces a full diagnostic PNG report for each algorithm tested.

---

## What it does

Given an MCMC algorithm and a bare energy function the checker:

1. **Builds the decision tree** by exhaustive BFS over all bit sequences the algorithm can consume. The state space is **discovered automatically** from a single seed state.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair using `FullSimplify` with β > 0. The result is an exact PASS or FAIL.

3. **Runs a numerical MCMC simulation** with genuine random numbers and compares the visited-state histogram to the Boltzmann distribution (KL divergence). Note: a low KL divergence does **not** guarantee detailed balance — an algorithm can sample the correct stationary distribution without being reversible (e.g. cyclic drift on a flat landscape).

4. **Renders a PNG report** showing the decision trees, symbolic transition matrix, detailed-balance pair table, and simulated vs Boltzmann scatter plot.

### Why bits?

Every random decision flows through `readBit[]` (returns 0 or 1) or `acceptTest[p]` (accept/reject with probability p). Fixing the bit sequence makes the algorithm deterministic and enables exhaustive path enumeration.

- `readBit[]` consumes one fair-coin bit; contributes ½ to the path weight.
- `acceptTest[p]` consumes one bit; contributes p (accept) or 1−p (reject).
- T(i→j) = **sum of path weights** over all bit sequences taking state i to state j.

### Automatic random-call interception

| Call | Converted to | Notes |
|---|---|---|
| `RandomReal[]` | deferred token → `acceptTest[p]` when compared | Works in `If[RandomReal[] < p, ...]` |
| `Random[]` | same as `RandomReal[]` | Deprecated Mathematica form |
| `RandomInteger[1]` or `RandomInteger[]` | `readBit[]` | Returns 0 or 1 |
| `RandomInteger[{lo, hi}]` | reads k=IntegerLength(hi−lo, 2) bits; out-of-range silently discarded | Any integer range |
| `RandomInteger[n]` | reads k=IntegerLength(n, 2) bits; out-of-range silently discarded | Any non-negative upper bound |
| `RandomChoice[list]` | reads k=IntegerLength(n−1, 2) bits; out-of-range silently discarded | Any list length |
| Float literals in `acceptTest` | rationalised: `0.5 → 1/2`, `1.0 → 1`, etc. | See float support below |

**Non-power-of-2 ranges via rejection sampling.** For n values, k = IntegerLength(n−1, 2) bits are read (always an exact integer). Bit strings mapping to values ≥ n are silently discarded. The missing probability fraction is the same from every starting state, so the unnormalised T still satisfies detailed balance exactly.

Calls that **cannot** be intercepted (cause analysis to abort):
- `RandomVariate[dist]` — continuous distributions
- `AbsoluteTime`, `SessionTime`, `Now` — time-dependent
- `RandomSample`, `RandomPermutation`

### Float support

Both float **acceptance probabilities** and float **energy values** are fully supported.

**Float acceptance probabilities** (inside the algorithm): Any `Real` value passed to `acceptTest` is rationalised via `Rationalize[]` at consumption time: `acceptTest[0.5]` → `acceptTest[1/2]`, `acceptTest[Exp[-β*1.0]]` → `acceptTest[Exp[-β]]`. This means `MetropolisProb` with a concrete float `dE` works correctly.

**Float energy values** (in the energy function): Any `Real` returned by `energy[s]` is rationalised in `CheckDetailedBalance` before passing to `FullSimplify`: `energy[s] = 0.5` → `1/2`. This means energy arrays like `{0., 1., 0.5}` work correctly.

**Float energy arrays are still NOT recommended** as primary input — use exact rationals (`{0, 1, 1/2}`) when possible. The rationalisation is a convenience, not a guarantee of correctness for all possible float values (e.g. `0.1` → `1/10`, which may differ from the intended physical value by machine-precision rounding).

### Energy function safety check

Before BFS, the checker scans the DownValues of the energy function for calls that would make the symbolic check meaningless:

- **Unsafe in energy**: `RandomReal`, `RandomInteger`, `RandomVariate`, `AbsoluteTime`, etc.
- **Safe in energy**: exact arithmetic, lookup tables, `Cos`, `Sin`, `Exp` with symbolic arguments, etc.

If unsafe calls are found, `RunFullCheck` aborts immediately with a clear error message. The same check is applied to the algorithm function.

---

## API

```mathematica
RunFullCheck[seedState, alg, energy, numBeta, options...]
```

| Argument | Type | Description |
|---|---|---|
| `seedState` | any | One valid state; others discovered automatically via BFS |
| `alg` | function | `alg[state, readBit, acceptTest]` — see below |
| `energy` | function | `energy[state]` → bare energy (no β). Integers, rationals, or floats. |
| `numBeta` | number | Inverse temperature β for the numerical MCMC run |

### Algorithm signature

```mathematica
myAlg[state_, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = If[b == 0, left[state], right[state]];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
```

The algorithm may also use `RandomInteger[]`, `RandomReal[]`, and `RandomChoice[]` natively. The `readBit_` and `acceptTest_` parameters may then be unused.

**Requirements:**

1. **Return a single next state.** Any Mathematica expression is valid as a state.
2. **Energy must be deterministic.** No random or time-dependent calls inside `energy`.
3. **β stays unassigned** during the symbolic check. Use `MetropolisProb[dE]` or any acceptance function that keeps β symbolic.
4. **States can be anything** — integers, rationals, lists, associations. The BFS discovers all reachable states from the seed automatically.

### `MetropolisProb` and custom acceptance criteria

```
MetropolisProb[dE] = Piecewise[{{1, dE ≤ 0}, {Exp[-β dE], dE > 0}}]
```

Custom criteria work as long as β remains symbolic:
```mathematica
BarkerProb[dE_] := 1 / (1 + Exp[β * dE])
```

### Options

| Option | Default | Description |
|---|---|---|
| `"SystemName"` | `"Unnamed system"` | Label shown in the report |
| `"AlgorithmCode"` | `Automatic` | Source string shown in the report |
| `"MaxBitDepth"` | `20` | Max BFS depth per state |
| `"TimeLimit"` | `60.` | Seconds per state before giving up |
| `"NSteps"` | `100000` | MCMC steps for the numerical check |
| `"WarmupFrac"` | `0.1` | Fraction of steps discarded as warmup |
| `"Verbose"` | `True` | Print BFS progress per state |
| `"OpenWindow"` | `True` | Open PNG report after each check |

---

## Running the examples

```bash
git clone https://github.com/Sam-Whitby/DetailedBalanceChecker.git ~/Desktop/DetailedBalanceChecker
cd ~/Desktop/DetailedBalanceChecker

wolframscript -file run_checks.wls          # core ring Kawasaki examples
wolframscript -file run_variable_bit.wls    # variable-bit-depth examples
wolframscript -file run_new_api.wls         # Barker + asymmetric proposal
wolframscript -file run_extended.wls        # all extended features
```

### Example summary

| Script | Example | Algorithm | Expected |
|---|---|---|---|
| run_checks | 1 | Kawasaki ring L=3, Metropolis | PASS |
| run_checks | 2 | Always-accept ring L=3 | FAIL |
| run_checks | 3 | Wrong acceptance sign | FAIL |
| run_variable_bit | 4 | Two-speed Metropolis L=4 | PASS |
| run_variable_bit | 5 | Biased rightward drift L=4 | FAIL |
| run_new_api | A | Kawasaki ring L=3, Metropolis | PASS |
| run_new_api | B | Always-accept ring L=3 | FAIL |
| run_new_api | C | Barker criterion ring L=4 | PASS |
| run_new_api | D | Asymmetric proposal + Metropolis L=4 | FAIL |
| run_extended | A | Cyclic drift L=4 (flat energy) | DB FAIL, numerical PASS |
| run_extended | B | Two-particle Kawasaki L=4 | PASS |
| run_extended | C | Metropolis on discretised circle | PASS |
| run_extended | D | Kawasaki with native RandomReal[] | PASS |
| run_extended | E | Native RandomReal[], wrong sign | FAIL |
| run_extended | F | RandomInteger[{1,3}] rejection sampling | PASS |
| run_extended | G | Float energy values + float probs | PASS |
| run_extended | H | Unanalyzable: RandomVariate in algorithm | aborts |
| run_extended | I | Unanalyzable: AbsoluteTime in algorithm | aborts |
| run_extended | J | Unanalyzable: AbsoluteTime in energy | aborts |
| run_extended | K | Unanalyzable: RandomVariate in energy | aborts |

---

## Running your own algorithm

### Quick start

```bash
# From the DetailedBalanceChecker directory:
wolframscript -e '
Get["dbc_core.wl"];

L = 4;
eps = {0, 1, 3, 2};    (* exact integers or rationals *)
numBeta = 1;

energy[s_Integer] := eps[[s]];

MyAlg[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = If[b == 0, Mod[state-2,L]+1, Mod[state,L]+1];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
];

RunFullCheck[1, MyAlg, energy, numBeta,
  "SystemName" -> "My algorithm",
  "OpenWindow" -> False]
'
```

### As a script file

Create `my_algorithm.wl`:
```mathematica
L        = 4
eps      = {0, 1, 3, 2}
numBeta  = 1

energy[s_Integer] := eps[[s]]

MyAlg[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = If[b == 0, Mod[s-2,L]+1, Mod[s,L]+1];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
```

Then run:
```bash
wolframscript -e '
Get["dbc_core.wl"];
Get["my_algorithm.wl"];
RunFullCheck[1, MyAlg, energy, numBeta,
  "SystemName" -> "My algorithm",
  "NSteps" -> 80000]
'
```

---

## Practical limits

- **States**: Discovered automatically; keep to ≤ 10–15 for practical run times.
- **Bit depth**: `MaxBitDepth=20` covers most algorithms.
- **FullSimplify**: The most expensive step. Piecewise Metropolis expressions typically simplify in seconds.
- **Float energies**: Supported — rationalised automatically in `CheckDetailedBalance`. Use exact rationals as primary input when possible.
- **Random call support**: `RandomReal[]`, `RandomInteger[]`, and `RandomChoice[]` intercepted for any range. `RandomVariate`, time functions aborted cleanly.
- **Numerical check and flat landscapes**: If all energies are equal the Boltzmann distribution is uniform, so KL ≈ 0 for ANY algorithm (even non-reversible ones). The symbolic DB check catches this — the numerical check is uninformative on flat landscapes by design.

---

## Files

```
dbc_core.wl                   Core library
show_report.py                Python/matplotlib report renderer
run_checks.wls                Core examples (3)
run_variable_bit.wls          Variable-bit-depth examples (2)
run_new_api.wls               Barker + asymmetric-proposal examples (4)
run_extended.wls              Extended feature examples (11, A-K)
examples/
  ring_kawasaki.wl            L=3 ring, Kawasaki Metropolis             [PASS]
  always_accept.wl            L=3 ring, always accept                  [FAIL]
  wrong_sign.wl               L=3 ring, wrong acceptance sign          [FAIL]
  variable_bit_pass.wl        L=4 ring, variable-bit two-speed         [PASS]
  variable_bit_fail.wl        L=4 ring, variable-bit biased drift      [FAIL]
  kawasaki_new.wl             L=3 ring, Metropolis and always-accept   [PASS/FAIL]
  barker_ring.wl              L=4 ring, Barker criterion               [PASS]
  asymmetric_proposal.wl      L=4 ring, biased proposal + Metropolis   [FAIL]
  cyclic_drift.wl             L=4 ring, always-right, flat energy      [DB-FAIL, num-PASS]
  multi_particle.wl           L=4 ring, 2 particles, hard-core         [PASS]
  continuous_metropolis.wl    Discretised circle, rational states       [PASS]
  random_native_pass.wl       Kawasaki using native RandomReal[]       [PASS]
  random_native_fail.wl       Wrong-sign using native RandomReal[]     [FAIL]
  nonpower_of_two.wl          RandomInteger[{1,3}] rejection sampling  [PASS]
  float_energy.wl             Float energy + float probs               [PASS]
  unanalyzable.wl             RandomVariate, AbsoluteTime in alg       [abort]
  unanalyzable_energy.wl      AbsoluteTime, RandomVariate in energy    [abort]
```
