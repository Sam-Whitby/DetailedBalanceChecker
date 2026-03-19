# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.
Produces a full diagnostic PNG report for each algorithm tested.

---

## What it does

Given an MCMC algorithm and a bare energy function the checker:

1. **Builds the decision tree** by exhaustive BFS over all bit sequences the algorithm can consume. Each leaf is a (state, probability-weight) pair. The state space is **discovered automatically** from a single seed state — no list of states is required.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair, using `FullSimplify` with β > 0. The result is an exact PASS or FAIL, not a statistical estimate.

3. **Runs a numerical MCMC simulation** with genuine random bits and compares the visited-state histogram to the Boltzmann distribution, reporting the KL divergence as a sanity check.

4. **Renders a PNG report** showing the decision trees, symbolic transition matrix, detailed-balance pair table, and simulated vs Boltzmann scatter plot.

### Why bits?

Every random decision in the algorithm comes from `readBit[]` (returns 0 or 1) or `acceptTest[p]` (accept/reject with probability p). By fixing the bit sequence, the checker makes the algorithm **deterministic** and can enumerate all execution paths exactly.

- `readBit[]` consumes one fair-coin bit; contributes factor ½ to the path weight.
- `acceptTest[p]` consumes one bit; contributes p (accept, bit=1) or 1−p (reject, bit=0) to the path weight.
- The transition probability T(i→j) is the **sum of path weights** over all bit sequences that take state i to state j.

The BFS extends a bit sequence by appending 0 and 1 whenever the algorithm requests more bits than are available. Variable-length paths (where different branches consume different numbers of bits) are handled automatically.

### Halting and safety

- `MaxBitDepth` (default 20) caps the BFS tree depth per state; paths that reach the cap are excluded with a warning.
- `TimeLimit` (seconds per state) provides a second safety net.
- `CheckAlgorithmSafety[alg]` scans the algorithm's definition for forbidden non-deterministic calls (`RandomReal`, `RandomInteger`, etc.) that would break BFS enumeration. It is called automatically by `RunFullCheck`.

---

## Generalised API (recommended)

The primary interface requires only a **single algorithm function** and a **single bare energy function**:

```mathematica
RunFullCheck[seedState, alg, energy, numBeta, options...]
```

| Argument | Type | Description |
|---|---|---|
| `seedState` | any | One valid state; others are discovered automatically |
| `alg` | function | `alg[state, readBit, acceptTest]` — see below |
| `energy` | function | `energy[state]` → bare energy (no β). **Must return exact values (integers or rationals, not floats).** |
| `numBeta` | number | Inverse temperature β for the numerical MCMC run |

### Algorithm function signature

```mathematica
myAlg[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];                            (* consume 1 fair bit: returns 0 or 1 *)
  nbr = If[b == 0, left[state], right[state]];
  dE  = energy[nbr] - energy[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
```

**Rules for algorithm construction:**

1. **Obtain all randomness from `readBit[]` or `acceptTest[p]`.** Never call `RandomReal[]`, `RandomInteger[]`, `RandomChoice[]`, or any other random/time function. The checker will flag these and they will break BFS.

2. **`acceptTest[p]`** takes the acceptance probability p (a number or symbolic expression in β) and returns 1 (accept) or 0 (reject). Use `MetropolisProb[dE]` or any other valid probability expression.

3. **`readBit[]`** returns 0 or 1. It is for discrete choices (which direction to propose, which sub-move to attempt, etc.), **not** for acceptance decisions — use `acceptTest` for those.

4. **Return a single next state.** The algorithm must return one state (integer or any label), not a list of weighted outcomes.

5. **Use exact energy values** (integers or rationals like `1/2`). Floating-point energies (like `0.5`) introduce `1.0 * β` instead of `β` in symbolic expressions, which `FullSimplify` may fail to simplify.

6. **Energy function is called directly inside the algorithm** as a global function. The same function is used for both the symbolic check (where β remains unassigned) and the MCMC run (where `Block[{β = numBeta}, ...]` assigns a numeric value to β).

### `MetropolisProb` and custom acceptance criteria

`MetropolisProb[dE]` is provided by the library:

```
MetropolisProb[dE] = Piecewise[{{1, dE ≤ 0}, {Exp[-β dE], dE > 0}}]
```

Any acceptance probability that keeps β symbolic works. For example the Barker criterion:

```mathematica
BarkerProb[dE_] := 1 / (1 + Exp[β * dE])
```

Both are valid to pass to `acceptTest[...]`.

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
# Clone
git clone https://github.com/Sam-Whitby/DetailedBalanceChecker.git ~/Desktop/DetailedBalanceChecker
cd ~/Desktop/DetailedBalanceChecker

# New generalised API (recommended)
wolframscript -file run_new_api.wls

# Legacy API examples (also still work)
wolframscript -file run_checks.wls
wolframscript -file run_variable_bit.wls
```

### What `run_new_api.wls` tests

| Example | Algorithm | Expected |
|---|---|---|
| A | Kawasaki ring L=3, Metropolis | PASS |
| B | Always-accept ring L=3 (no Metropolis) | FAIL |
| C | Barker criterion ring L=4 | PASS |
| D | Asymmetric proposal (3:1 bias) + Metropolis L=4 | FAIL |

---

## Writing a new algorithm

Create a `.wl` file in `examples/` following this template:

```mathematica
(* ---- system parameters ---- *)
L        = 4
eps      = {0, 1, 3, 2}          (* exact integers or rationals -- no floats *)
numBeta  = 1                      (* exact integer or rational *)

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

---

## Legacy API

The original five-argument form is still supported:

```mathematica
RunFullCheck[allStates, symAlg, numAlg, symEnergy, numEnergy, options...]
```

This requires separate symbolic and numeric versions of the algorithm and energy. See `examples/ring_kawasaki.wl` for a complete example. Use the generalised API for new work.

---

## Files

```
dbc_core.wl                   Core library — load this first
show_report.py                Python/matplotlib report renderer
run_new_api.wls               Generalised API runner (4 examples)
run_checks.wls                Legacy API runner (3 examples)
run_variable_bit.wls          Variable-bit-depth examples
examples/
  kawasaki_new.wl             L=3 ring, Metropolis and always-accept   [new API]
  barker_ring.wl              L=4 ring, Barker (heat-bath) criterion   [new API, PASS]
  asymmetric_proposal.wl      L=4 ring, biased proposal + Metropolis   [new API, FAIL]
  ring_kawasaki.wl            L=3 ring, Metropolis                     [legacy API, PASS]
  always_accept.wl            L=3 ring, always accept                  [legacy API, FAIL]
  wrong_sign.wl               L=3 ring, wrong acceptance sign          [legacy API, FAIL]
  variable_bit_pass.wl        L=4 ring, variable-bit two-speed hop     [legacy API, PASS]
  variable_bit_fail.wl        L=4 ring, variable-bit biased drift      [legacy API, FAIL]
```

---

## Practical limits

The symbolic check is exact but grows quickly:

- **States**: Discovered automatically; keep to ≤ 10 for practical run times.
- **Bit depth**: Each level doubles paths. `MaxBitDepth=20` covers most algorithms; simple algorithms use 2–4 bits.
- **FullSimplify**: Most expensive step. Piecewise Metropolis expressions with `{β > 0}` assumptions typically simplify in seconds. Unusual acceptance functions (nested exponentials, etc.) may be slow.
- **Float warning**: Energy values must be exact (integers/rationals). Using `0.5` instead of `1/2` introduces floating-point coefficients that break symbolic simplification.
