# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.
Produces a full diagnostic PNG report for each algorithm tested.

---

## What it does

Given an MCMC algorithm and a bare energy function the checker:

1. **Builds the decision tree** by exhaustive BFS over all bit sequences the algorithm can consume. The state space is **discovered automatically** from a single seed state — no list of states is required.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair, using `FullSimplify` with β > 0. The result is an exact PASS or FAIL, not a statistical estimate.

3. **Runs a numerical MCMC simulation** with genuine random bits and compares the visited-state histogram to the Boltzmann distribution, reporting the KL divergence as a sanity check. Note: a low KL divergence does **not** guarantee detailed balance — an algorithm can sample the correct stationary distribution without being reversible (e.g. a cyclic drift on a flat energy landscape).

4. **Renders a PNG report** showing the decision trees, symbolic transition matrix, detailed-balance pair table, and simulated vs Boltzmann scatter plot.

### Why bits?

Every random decision in the algorithm flows through `readBit[]` (returns 0 or 1) or `acceptTest[p]` (accept/reject with probability p). By fixing the bit sequence the checker makes the algorithm **deterministic** and enumerates all execution paths exactly.

- `readBit[]` consumes one fair-coin bit; contributes factor ½ to the path weight.
- `acceptTest[p]` consumes one bit; contributes p (accept, bit=1) or 1−p (reject, bit=0).
- T(i→j) = **sum of path weights** over all bit sequences that take state i to state j.

The BFS extends a bit sequence by appending 0 and 1 whenever the algorithm requests more bits. Variable-length paths (different branches consuming different numbers of bits) are handled automatically.

### Automatic random-call interception

Algorithms may use native Mathematica random functions. The checker intercepts these automatically during BFS:

| Call | Converted to | Notes |
|---|---|---|
| `RandomReal[]` | deferred token → `acceptTest[p]` when compared to `p` | Works in `If[RandomReal[] < p, ...]` |
| `Random[]` | same as `RandomReal[]` | Deprecated Mathematica form |
| `RandomInteger[1]` or `RandomInteger[]` | `readBit[]` | Returns 0 or 1 |
| `RandomInteger[{lo, hi}]` | reads k=⌈log₂(hi−lo+1)⌉ bits; out-of-range values silently discarded | Any integer range, power-of-2 or not |
| `RandomInteger[n]` | reads k=⌈log₂(n+1)⌉ bits; out-of-range values silently discarded | Any non-negative upper bound |
| `RandomChoice[list]` | reads k=⌈log₂(n)⌉ bits; out-of-range indices silently discarded | Any list length |
| Float literals in `acceptTest` | rationalised: `0.5 → 1/2`, `0.12 → 3/25`, etc. | Algorithms can use floats without breaking symbolic simplification |

**Non-power-of-2 ranges via rejection sampling.** For a range of n values the checker reads k = ⌈log₂(n)⌉ bits (2ᵏ outcomes). Bit strings that map to values ≥ n are silently discarded from the BFS — they represent the "rejection" branch. The surviving paths each have equal weight 1/2ᵏ and form the correct uniform distribution. The missing probability fraction is the same from every starting state, so the unnormalised transition matrix still satisfies detailed balance exactly.

Calls that **cannot** be intercepted cause the analysis to abort with a clear error message:
- `RandomVariate[dist]` — continuous distributions
- `AbsoluteTime`, `SessionTime`, `Now` — time-dependent values
- `RandomSample`, `RandomPermutation` — unsupported combinatorial sampling

The numerical MCMC run always uses genuine random numbers regardless of interception.

---

## API

```mathematica
RunFullCheck[seedState, alg, energy, numBeta, options...]
```

| Argument | Type | Description |
|---|---|---|
| `seedState` | any | One valid state; others discovered automatically via BFS |
| `alg` | function | `alg[state, readBit, acceptTest]` — see below |
| `energy` | function | `energy[state]` → bare energy (no β). Use exact values (integers or rationals). |
| `numBeta` | number | Inverse temperature β for the numerical MCMC run |

### Algorithm function signature

```mathematica
myAlg[state_, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];                            (* consume 1 fair bit: returns 0 or 1 *)
  nbr = If[b == 0, left[state], right[state]];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
```

The algorithm can also use `RandomInteger[]` and `RandomReal[]` directly (see interception table above). The `readBit_` and `acceptTest_` parameters may be left unused if all randomness comes from native calls.

**Rules:**

1. **Return a single next state.** The algorithm must return one state (any Mathematica expression), not a list.

2. **Use exact energy values** (integers or rationals like `1/2`). Floating-point energies introduce numerical coefficients that break symbolic simplification. Floating-point *probabilities* passed to `acceptTest` are fine — they are rationalised automatically.

3. **β stays global and unassigned** during the symbolic DB check. Assign it numerically via `Block[{β = numBeta}, ...]` if needed inside the algorithm, or use `MetropolisProb[dE]` which handles this automatically.

4. **States can be anything** — integers, rationals, lists, associations. The checker works with any Mathematica expression as a state. This means multi-particle systems, continuous discretisations, and higher-dimensional lattices all work without special configuration.

### `MetropolisProb` and custom acceptance criteria

`MetropolisProb[dE]` is provided by the library:
```
MetropolisProb[dE] = Piecewise[{{1, dE ≤ 0}, {Exp[-β dE], dE > 0}}]
```

Any acceptance probability that keeps β symbolic works, for example the Barker criterion:
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

# Core examples (ring Kawasaki variants)
wolframscript -file run_checks.wls

# Variable-bit-depth examples
wolframscript -file run_variable_bit.wls

# New API examples (Barker, asymmetric proposal)
wolframscript -file run_new_api.wls

# Extended examples (multi-particle, continuous, native random calls,
#                    non-power-of-2, float literals, unanalyzable)
wolframscript -file run_extended.wls
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
| run_extended | G | Float literal probs (0.5, 1.0) | PASS |
| run_extended | H | Unanalyzable: RandomVariate | aborts |
| run_extended | I | Unanalyzable: AbsoluteTime | aborts |

---

## Writing a new algorithm

Create a `.wl` file in `examples/` following this template:

```mathematica
(* ---- system parameters ---- *)
L        = 4
eps      = {0, 1, 3, 2}          (* exact integers or rationals -- no floats *)
numBeta  = 1

rightOf[s_Integer] := Mod[s,     L] + 1
leftOf[s_Integer]  := Mod[s - 2, L] + 1

energy[s_Integer] := eps[[s]]    (* bare energy, no beta *)

(* ---- algorithm ---- *)
MyAlg[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = If[b == 0, leftOf[state], rightOf[state]];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
```

Then in your runner:

```mathematica
Get["dbc_core.wl"]
Get["examples/my_algorithm.wl"]

RunFullCheck[
  1,          (* seed state *)
  MyAlg,
  energy,
  numBeta,
  "SystemName" -> "My algorithm",
  "AlgorithmCode" -> ReadString["examples/my_algorithm.wl"],
  "NSteps" -> 80000
]
```

### Multi-particle and non-integer states

The checker places no restrictions on state types. Some examples:

```mathematica
(* Two-particle state as a sorted list *)
KawasakiMulti[{p1_Integer, p2_Integer}, readBit_, acceptTest_] := ...

(* Continuous position as a rational number *)
ContinuousMetropolis[x_Rational, readBit_, acceptTest_] := ...
```

Starting from a single seed state the BFS discovers all reachable states automatically, whatever type they have.

---

## Practical limits

The symbolic check is exact but grows quickly:

- **States**: Discovered automatically; keep to ≤ 10–15 for practical run times.
- **Bit depth**: Each level doubles paths. `MaxBitDepth=20` covers most algorithms.
- **FullSimplify**: The most expensive step. Piecewise Metropolis expressions with `{β > 0}` assumptions typically simplify in seconds; highly nested expressions may be slow.
- **Float energies**: Must be exact (integers/rationals). Using `0.5` instead of `1/2` for *energy values* makes Piecewise conditions evaluate numerically during the symbolic phase, corrupting the result. Float *acceptance probabilities* passed to `acceptTest` are automatically rationalised.
- **Random call support**: `RandomReal[]`, `RandomInteger[]`, and `RandomChoice[]` are intercepted automatically for any range. `RandomVariate` and time functions cannot be analysed.

---

## Files

```
dbc_core.wl                   Core library — load this first
show_report.py                Python/matplotlib report renderer
run_checks.wls                Core examples (3 ring Kawasaki variants)
run_variable_bit.wls          Variable-bit-depth examples (2)
run_new_api.wls               Barker + asymmetric-proposal examples (4)
run_extended.wls              Extended feature examples (9)
examples/
  kawasaki_new.wl             L=3 ring, Metropolis and always-accept    [PASS/FAIL]
  barker_ring.wl              L=4 ring, Barker (heat-bath) criterion    [PASS]
  asymmetric_proposal.wl      L=4 ring, biased proposal + Metropolis    [FAIL]
  ring_kawasaki.wl            L=3 ring, Kawasaki Metropolis             [PASS]
  always_accept.wl            L=3 ring, always accept                  [FAIL]
  wrong_sign.wl               L=3 ring, wrong acceptance sign          [FAIL]
  variable_bit_pass.wl        L=4 ring, variable-bit two-speed hop     [PASS]
  variable_bit_fail.wl        L=4 ring, variable-bit biased drift      [FAIL]
  cyclic_drift.wl             L=4 ring, always-right, flat energy      [DB-FAIL, num-PASS]
  multi_particle.wl           L=4 ring, 2 particles, hard-core         [PASS]
  continuous_metropolis.wl    Discretised circle, rational states       [PASS]
  random_native_pass.wl       Kawasaki using RandomReal[]/RandomInteger [PASS]
  random_native_fail.wl       Wrong-sign Kawasaki with RandomReal[]    [FAIL]
  nonpower_of_two.wl          RandomInteger[{1,3}] rejection sampling  [PASS]
  float_energy.wl             Float literals (0.5, 1.0) rationalised   [PASS]
  unanalyzable.wl             Two algorithms with unsupported calls    [abort]
```
