# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.
Produces a full diagnostic PNG report for each algorithm tested.

---

## What it does

Given an MCMC algorithm and an energy function the checker:

1. **Builds the decision tree** by exhaustive BFS over all random-call sequences the algorithm can make. The state space is **discovered automatically** from a single seed state.

2. **Checks detailed balance exactly** — constructs the symbolic transition matrix and verifies T(i→j)·π(i) = T(j→i)·π(j) for every pair using `FullSimplify` with β > 0. Both β and all energy parameters are kept symbolic throughout. The result is an exact PASS or FAIL.

3. **Runs a numerical MCMC simulation** with genuine random numbers and compares the visited-state histogram to the Boltzmann distribution (KL divergence). When symbolic energy parameters are supplied, random numerical values are assigned automatically (with a user-specified seed).

4. **Renders a PNG report** showing the decision trees, symbolic transition matrix, detailed-balance pair table, and simulated vs Boltzmann scatter plot.

### How it works internally

The algorithm is called with native Mathematica random functions (`RandomInteger[]`, `RandomReal[]`, etc.), which are intercepted during BFS. Each call is replaced by a deterministic value from a bit tape, making the algorithm's execution path fully reproducible.

- `RandomInteger[]` or `RandomInteger[1]` — reads one bit (0 or 1); contributes ½ to path weight.
- `RandomReal[] < p` — reads one bit; contributes p (accept) or 1−p (reject) to path weight.
- T(i→j) = **sum of path weights** over all bit sequences taking state i to state j.

### Supported random calls

| Call | Behaviour during BFS | Notes |
|---|---|---|
| `RandomReal[]` | deferred token; fires as accept/reject when compared with `< p` | Works in `If[RandomReal[] < p, ...]` |
| `Random[]` | same as `RandomReal[]` | Deprecated Mathematica form |
| `RandomInteger[]` or `RandomInteger[1]` | reads one bit; returns 0 or 1 | |
| `RandomInteger[{lo, hi}]` | reads k=⌈log₂(hi−lo+1)⌉ bits; out-of-range silently discarded | Any integer range |
| `RandomInteger[n]` | reads k=⌈log₂(n+1)⌉ bits; out-of-range silently discarded | Any non-negative upper bound |
| `RandomChoice[list]` | reads k=⌈log₂(n)⌉ bits; out-of-range silently discarded | Any list length |

**Non-power-of-2 ranges via rejection sampling.** Bit strings mapping to out-of-range values are silently discarded. The missing probability fraction is identical from every starting state, so the unnormalised T still satisfies detailed balance exactly.

Calls that **cannot** be intercepted (cause analysis to abort):
- `RandomVariate[dist]` — continuous distributions
- `AbsoluteTime`, `SessionTime`, `Now` — time-dependent
- `RandomSample`, `RandomPermutation`

---

## Symbolic energy parameters

All energy parameters — site energies **and** pairwise coupling strengths — are kept as unassigned symbols during the symbolic check. The transition matrix entries are therefore symbolic expressions in β, the site energies εᵢ, and the coupling constants Jd. `FullSimplify` with β > 0 verifies detailed balance across all possible energy values simultaneously.

For the numerical MCMC check, random values are assigned automatically to the energy symbols (via `Block`) using a user-specified seed, so results are fully reproducible.

### System parameter helpers

```mathematica
(* Create a symbolic parameter set for an L-site ring *)
params = MakeRingParams[L, "prefix", nCouplingDistances]
(* Returns <|"L" -> L, "eps" -> {ε_prefix_1,...}, "couplings" -> {J_prefix_1,...}|> *)

(* Build the energy function from params *)
energy = BuildRingEnergy[params]
(* energy[s_Integer] for single particle; energy[{p1,p2,...}] for multi-particle *)
```

### Energy with pairwise interactions

The energy of a multi-particle state `{p1, p2, ...}` on an L-site periodic ring:

```
E({p1,...,pn}) = Σᵢ ε_pᵢ + Σᵢ<ⱼ J_{d(pᵢ,pⱼ)}
```

where `d(a,b) = min(|a−b|, L−|a−b|)` is the ring distance and Jd is the coupling at distance d. Hard-sphere exclusion (d=0, same site) is enforced by the algorithm, not the energy function.

---

## API

```mathematica
RunFullCheck[seedState, alg, energy, numBeta, options...]
```

| Argument | Type | Description |
|---|---|---|
| `seedState` | any | One valid state; others discovered automatically via BFS |
| `alg` | function | `alg[state]` — uses native random calls |
| `energy` | function | `energy[state]` → bare energy (no β). May be symbolic in energy parameters. |
| `numBeta` | number | Inverse temperature β for the numerical MCMC run |

### Algorithm signature

```mathematica
MyAlg[state_Integer] := Module[
  {dir, nbr, dE},
  dir = RandomInteger[];
  nbr = If[dir == 0, left[state], right[state]];
  dE  = energy[nbr] - energy[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
]
```

**Requirements:**

1. **Return a single next state.** Any Mathematica expression is valid as a state.
2. **Energy must be deterministic.** No random or time-dependent calls inside `energy`.
3. **β stays unassigned** during the symbolic check. Use `MetropolisProb[dE]` or any acceptance function that keeps β symbolic.
4. **States can be anything** — integers, rationals, lists. The BFS discovers all reachable states from the seed automatically.

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
| `"SysParams"` | `None` | Association from `MakeRingParams` with symbolic energy parameters. When supplied, random numerical values are assigned during the MCMC check. |
| `"NumericSeed"` | `42` | Integer seed passed to `SeedRandom` for reproducible random energy values. |
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
wolframscript -file run_variable_bit.wls    # variable-depth examples
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
| run_extended | B | Two-particle Kawasaki + NN coupling L=4 | PASS |
| run_extended | C | Metropolis on discretised circle | PASS |
| run_extended | D | Kawasaki with native RandomReal[] | PASS |
| run_extended | E | Native RandomReal[], wrong sign | FAIL |
| run_extended | F | RandomInteger[{1,3}] rejection sampling | PASS |
| run_extended | G | Two-particle + NN coupling, wrong sign | FAIL |
| run_extended | H | Three-particle Kawasaki + NN+NNN coupling L=4 | PASS |
| run_extended | I | Unanalyzable: RandomVariate in algorithm | aborts |
| run_extended | J | Unanalyzable: AbsoluteTime in algorithm | aborts |
| run_extended | K | Unanalyzable: AbsoluteTime in energy | aborts |
| run_extended | L | Unanalyzable: RandomVariate in energy | aborts |

---

## Running your own algorithm

### Single-line terminal usage with `check.wls`

```bash
wolframscript -file check.wls <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> [options]
```

| Argument | Description |
|---|---|
| `alg_file` | Path to your `.wl` file (absolute, or relative to the checker directory) |
| `seed_state` | A valid starting state — integer, or `{i,j}` for multi-particle |
| `alg_name` | Name of the algorithm function defined in the file |
| `energy_name` | Name of the energy function defined in the file |
| `num_beta` | Inverse temperature β for the numerical run |

Optional `key=value` flags:

| Option | Default |
|---|---|
| `NSteps=N` | `100000` |
| `MaxBitDepth=N` | `20` |
| `SystemName=label` | algorithm name |
| `Verbose=True/False` | `True` |
| `OpenWindow=True/False` | `True` |

### Writing your algorithm file

Create `my_algorithm.wl` using symbolic energy parameters:

```mathematica
Get["dbc_core.wl"]   (* provides MakeRingParams, BuildRingEnergy, RingDist *)

params  = MakeRingParams[4, "my"]   (* symbolic eps for 4-site ring *)
energy  = BuildRingEnergy[params]
numBeta = 1

MyAlg[state_Integer] := Module[
  {dir, nbr, dE},
  dir = RandomInteger[];
  nbr = If[dir == 0, Mod[state-2, 4]+1, Mod[state, 4]+1];
  dE  = energy[nbr] - energy[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
]
```

Then run:
```bash
wolframscript -file check.wls my_algorithm.wl 1 MyAlg energy 1
```

### Inline usage

```bash
wolframscript -e '
Get["dbc_core.wl"];

params  = MakeRingParams[4, "my", 1];   (* site energies + NN coupling *)
energy  = BuildRingEnergy[params];
numBeta = 1;

MyAlg[state_Integer] := Module[
  {dir, nbr, dE},
  dir = RandomInteger[];
  nbr = If[dir == 0, Mod[state-2, 4]+1, Mod[state, 4]+1];
  dE  = energy[nbr] - energy[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
];

RunFullCheck[1, MyAlg, energy, numBeta,
  "SysParams"   -> params,
  "NumericSeed" -> 42,
  "SystemName"  -> "My algorithm",
  "OpenWindow"  -> False]
'
```

---

## Practical limits

- **States**: Discovered automatically; keep to ≤ 10–15 for practical run times.
- **Bit depth**: `MaxBitDepth=20` covers most algorithms.
- **FullSimplify**: The most expensive step. With symbolic energies, Piecewise conditions are resolved branch-by-branch using exact algebra — typically a few seconds per pair.
- **Random call support**: `RandomReal[]`, `RandomInteger[]`, and `RandomChoice[]` intercepted for any range. `RandomVariate`, time functions aborted cleanly.
- **Numerical check and flat landscapes**: If all energies are equal the Boltzmann distribution is uniform, so KL ≈ 0 for ANY algorithm (even non-reversible ones). The symbolic DB check catches this — the numerical check is uninformative on flat landscapes by design.

---

## Files

```
dbc_core.wl                   Core library
show_report.py                Python/matplotlib report renderer
check.wls                     Single-line terminal interface
run_checks.wls                Core examples (3)
run_variable_bit.wls          Variable-depth examples (2)
run_new_api.wls               Barker + asymmetric-proposal examples (4)
run_extended.wls              Extended feature examples (12, A-L)
examples/
  ring_kawasaki.wl            L=3 ring, Kawasaki Metropolis             [PASS]
  always_accept.wl            L=3 ring, always accept                  [FAIL]
  wrong_sign.wl               L=3 ring, wrong acceptance sign          [FAIL]
  variable_bit_pass.wl        L=4 ring, variable-depth two-speed       [PASS]
  variable_bit_fail.wl        L=4 ring, variable-depth biased drift    [FAIL]
  kawasaki_new.wl             L=3 ring, Metropolis and always-accept   [PASS/FAIL]
  barker_ring.wl              L=4 ring, Barker criterion               [PASS]
  asymmetric_proposal.wl      L=4 ring, biased proposal + Metropolis   [FAIL]
  cyclic_drift.wl             L=4 ring, always-right, flat energy      [DB-FAIL, num-PASS]
  multi_particle.wl           L=4 ring, 2 particles + NN coupling      [PASS]
  coupling_3particle.wl       L=4 ring, 3 particles + NN+NNN coupling  [PASS]
  continuous_metropolis.wl    Discretised circle, rational states       [PASS]
  random_native_pass.wl       Kawasaki, native random calls            [PASS]
  random_native_fail.wl       Wrong-sign, native random calls          [FAIL]
  nonpower_of_two.wl          RandomInteger[{1,3}] rejection sampling  [PASS]
  unanalyzable.wl             RandomVariate, AbsoluteTime in alg       [abort]
  unanalyzable_energy.wl      AbsoluteTime, RandomVariate in energy    [abort]
```
