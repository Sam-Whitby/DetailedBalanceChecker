# DetailedBalanceChecker

Exhaustive symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.

---

## How it works

Given an MCMC algorithm, an energy function, and a way to enumerate states, the checker:

1. **Enumerates states by BFS.** Starting from a seed state, it exhaustively runs the algorithm on every possible sequence of random bits, discovering all reachable states automatically.

2. **Constructs the exact transition matrix.** Each entry T(i→j) is computed as the sum of path weights (products of 1/2 per fair bit and acceptance probability per Metropolis step) over all bit sequences that take state i to state j.

3. **Checks detailed balance symbolically.** For every pair (i, j), verifies T(i→j)·e^{−β·E(i)} = T(j→i)·e^{−β·E(j)} using `FullSimplify` with β > 0 and all energy parameters treated as real symbols. The result is an exact PASS or FAIL — no floating-point approximation.

4. **Validates numerically.** Runs a genuine MCMC simulation and compares the visited-state histogram to the Boltzmann distribution using the KL divergence. A low KL divergence (< 0.02) is a PASS.

The symbolic check is rigorous; the numerical check is a useful sanity check but can give false PASSes (e.g. on flat energy landscapes where any ergodic algorithm gives a uniform histogram).

---

## Running the checker

```bash
cd ~/Desktop/DetailedBalanceChecker

# 1D Kawasaki (all L=1,2,3,4 systems), symbolic + numerical
wolframscript -file check.wls examples3/kawasaki_1d.wl MaxBitString=1111111

# 2D Kawasaki (L=1 and L=2 tori), symbolic only (faster)
wolframscript -file check.wls examples3/kawasaki_2d.wl MaxBitString=1111111 Mode=Symbolic

# Failure examples
wolframscript -file check.wls examples3/kawasaki_1d_fail.wl MaxBitString=1111111
wolframscript -file check.wls examples3/kawasaki_2d_fail.wl MaxBitString=1111111 Mode=Symbolic
```

**Options** (`key=value`, no spaces around `=`):

| Option | Default | Description |
|---|---|---|
| `MaxBitString=XXXX` | `11111111` | Test all bit strings from length 1 up through `len(XXXX)`, stopping at `XXXX` within the final length. Leading zeros are significant. |
| `Mode=Symbolic` | `Both` | Symbolic detailed-balance check only |
| `Mode=Numerical` | `Both` | Numerical MCMC check only |
| `Mode=Both` | `Both` | Both checks |
| `NSteps=N` | `50000` | MCMC steps per tested component |
| `MaxBitDepth=N` | `20` | BFS tree depth cap per state |
| `Verbose=True` | `False` | Print per-state BFS progress |

### How MaxBitString determines which systems are tested

The checker iterates every bit string from length 1 up to and including `MaxBitString`, in order of increasing length then numerically within each length. For each bit string it calls `BitsToState` (defined in your `.wl` file) to get a seed state. If the seed belongs to a previously-discovered connected component it is skipped (deduplication). Otherwise BFS discovers the full component and one row is printed.

For the Kawasaki examples, the bijective encoding maps 7-bit strings to all L=1,2,3,4 labeled-particle arrays, so `MaxBitString=1111111` covers all those systems exhaustively.

---

## Animating an algorithm

```bash
wolframscript -file animate.wls examples3/kawasaki_1d.wl \
  Sites=6 N=3 Steps=300 Beta=1 Delay=0.15

wolframscript -file animate.wls examples3/kawasaki_2d.wl \
  Sites=9 N=4 Steps=200 Delay=0.2

# Supply specific coupling constants
wolframscript -file animate.wls examples3/kawasaki_1d.wl \
  Sites=6 N=3 J12=1.0 J13=0.5 J23=2.0
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `Sites=<n>` | required | Total lattice sites (e.g. 9 for a 3x3 grid) |
| `N=<n>` | required | Number of labeled particles (types 1..N) |
| `Steps=<n>` | `200` | Number of MCMC steps to animate |
| `Beta=<f>` | from `.wl` file | Inverse temperature |
| `Delay=<f>` | `0.1` | Seconds per frame |
| `J<a><b>=<f>` | random | Set coupling constant between types a and b |

---

## Writing your own algorithm file

Copy `template.wl` and fill in the four required definitions. All four live in a single self-contained `.wl` file.

### Required: `energy[state_]`

Bare energy (no beta). Must be deterministic. May use symbolic variables declared in `symParams` or `DynamicSymParams`.

### Required: `Algorithm[state_]`

MCMC move. Use native Mathematica random calls: `RandomInteger[]`, `RandomReal[]`, `RandomChoice[]`. Use `MetropolisProb[dE]` for Metropolis acceptance:

```mathematica
If[RandomReal[] < MetropolisProb[dE], newState, state]
```

The checker intercepts these calls during BFS, replacing them with deterministic values from a bit tape. `MetropolisProb[dE]` expands to `Piecewise[{{1, dE<=0}, {Exp[-beta*dE], dE>0}}]` with beta kept symbolic.

**Supported random calls:**

| Call | Bits consumed |
|---|---|
| `RandomReal[]` | 1 (deferred until compared with `< p`) |
| `RandomInteger[{lo, hi}]` | ceil(log2(hi-lo+1)) bits |
| `RandomInteger[]` or `RandomInteger[1]` | 1 bit |
| `RandomChoice[list]` | ceil(log2(n)) bits |

Out-of-range rejection-sampling paths are silently discarded; the missing probability is uniform across all states, so detailed balance is still checked correctly.

### Required: `BitsToState[bits_List]`

Converts a bit string to a seed state, or returns `None` to skip it.

```mathematica
(* Simple example: 4-site ring, exactly 2 particles *)
L = 4
BitsToState[bits_List] :=
  If[Length[bits] =!= L || Total[bits] =!= 2, None,
    Flatten[Position[bits, 1]]]
```

### Required: `numBeta`

```mathematica
numBeta = 1
```

### Optional: symbolic parameters

**Static** (same for every component):
```mathematica
symParams = <|"eps" -> {eps1, eps2}, "couplings" -> {J12, J13}|>
```

**Dynamic** (auto-generated per component from the discovered states). Use this when parameters depend on which particle types are present — no hardcoded limit needed:
```mathematica
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
        {a, types}, {b, types}]|>]
```

Both may be defined simultaneously; their parameter lists are merged. You can declare any symbolic quantities your energy uses: second-order couplings, particle masses, friction coefficients, time steps, spring constants, etc. The checker keeps them symbolic during the exact check and assigns random numeric values for the MCMC run.

### Optional: `DisplayState[state_]`

If defined, replaces `ToString[state]` in the output table. For 2D states:
```mathematica
DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}], "|"]]
(* {1,2,3,0} -> "{1,2}|{3,0}" *)
```

---

## Examples

| File | System | Result |
|---|---|---|
| `examples3/kawasaki_1d.wl` | 1D ring, labeled particles, dynamic L | **PASS** |
| `examples3/kawasaki_2d.wl` | 2D torus, labeled particles, dynamic L | **PASS** |
| `examples3/kawasaki_1d_fail.wl` | Same as 1D but bond 0 proposed twice as often | **FAIL** |
| `examples3/kawasaki_2d_fail.wl` | Same as 2D but bond 0 proposed twice as often | **FAIL** |

The failure examples differ from the correct ones by exactly one line:

```mathematica
(* Correct *)    dE = energy[newState] - energy[state]
(* Buggy   *)    dE = energy[state] - energy[newState]   (* sign reversed *)
```

With the wrong sign, `MetropolisProb` receives `-ΔE`: uphill moves (to higher energy) are always accepted, while downhill moves are penalised with `exp(-β|ΔE|)`. The algorithm preferentially visits high-energy states, inverting the Boltzmann distribution. For any pair of states (i, j) with E(i) ≠ E(j), T(i→j)·π(i) ≠ T(j→i)·π(j).

**Note on bond-selection bias.** An asymmetric bond proposal — e.g. `Mod[RandomInteger[{0,L}],L]`, which makes bond 0 twice as likely — does *not* break detailed balance for Kawasaki. Both the forward transition i→j and the reverse transition j→i go through the same bond b and therefore have the same proposal probability, which cancels in the DB condition. Only the acceptance step matters.

---

## Checking your own algorithm

The checker is completely agnostic to the algorithm type. To check a custom algorithm:

1. Copy `template.wl`.
2. Implement `energy[state_]` for your system.
3. Implement `Algorithm[state_]` using `RandomInteger`, `RandomReal`, `RandomChoice`, and `MetropolisProb`.
4. Implement `BitsToState[bits_]` to enumerate the configurations you want to test.
5. Declare symbolic parameters in `symParams` or `DynamicSymParams`.

The state can be any Mathematica expression — a list, an association, a matrix. The checker only requires that the algorithm accepts and returns the same type, and that BitsToState returns the same type.

For a **discretised molecular dynamics** integrator (where, say, one lattice dimension represents velocity or momentum): the state is just a list of integers. The Metropolis accept/reject step is handled identically. Symbolic parameters like mass `m`, timestep `dt`, and friction `gamma` go in `symParams`.

---

## Limitations

- **State space size.** BFS discovers the full connected component. Keep components to ~20 states for fast symbolic checks; larger (e.g. 84 states) is feasible but slow.
- **Random call types.** `RandomVariate`, `AbsoluteTime`, `SessionTime`, and `Now` cannot be intercepted and cause an error row in the output.
- **Flat-energy caveat.** If all energies are equal, the Boltzmann distribution is uniform so the numerical check gives PASS for any ergodic algorithm. The symbolic check is not affected.
- **Numerical coupling range.** Random coupling values are drawn from (-1, 1). If your algorithm is only stable for specific parameter ranges, the numerical check may give misleading results; trust the symbolic check.

---

## Files

```
check.wls               Checker entry point
animate.wls             Terminal animation of algorithm execution
template.wl             Template for writing an algorithm file
dbc_core.wl             Core library (BFS, symbolic check, MCMC check)
examples3/
  kawasaki_1d.wl        1D Kawasaki, dynamic L, labeled particles   [PASS]
  kawasaki_2d.wl        2D Kawasaki, dynamic L, labeled particles   [PASS]
  kawasaki_1d_fail.wl   Same with bond-0 bias                       [FAIL]
  kawasaki_2d_fail.wl   Same with bond-0 bias                       [FAIL]
legacy/
  examples/             Old example format (unsupported)
  examples2/            Fixed-L examples (superseded by examples3/)
```
