# DetailedBalanceChecker

Exhaustive symbolic and numerical checker for the detailed balance condition in MCMC algorithms.
Written in Mathematica, run from the terminal via `wolframscript`.

---

## How it works

Given an MCMC algorithm file, the checker:

1. **Enumerates states by BFS.** Starting from a seed state, it exhaustively runs the algorithm on every possible sequence of random bits, discovering all reachable states automatically.

2. **Constructs the exact transition matrix.** Each entry T(i→j) is the sum of path weights (products of acceptance probabilities) over all bit sequences that take state i to state j.

3. **Checks detailed balance symbolically.** For every pair (i, j), verifies T(i→j)·e^{−β·E(i)} = T(j→i)·e^{−β·E(j)} using `FullSimplify` with β and all energy parameters treated as real symbols. The result is an exact PASS or FAIL.

4. **Validates numerically.** Runs a genuine MCMC simulation and compares the visited-state histogram to the Boltzmann distribution using the KL divergence. A KL divergence below 0.02 is a PASS.

---

## Running the checker

```bash
cd ~/Desktop/DetailedBalanceChecker

# 1D Kawasaki (all L=1..4 systems)
wolframscript -file check.wls examples3/kawasaki_1d.wl MaxBitString=1111111

# 2D Kawasaki
wolframscript -file check.wls examples3/kawasaki_2d.wl MaxBitString=1111111 Mode=Symbolic

# 2D VMMC — small systems (L=1 and L=2)
wolframscript -file check.wls examples3/vmmc_2d.wl MaxBitString=1111111 Mode=Symbolic

# 2D VMMC — up to and including L=2 (all 2×2 grid states; ~18-bit ID range)
wolframscript -file check.wls examples3/vmmc_2d.wl MaxBitString=1011000 Mode=Symbolic

# 2D VMMC — report for a specific 3×3 seed state
wolframscript -file report.wls examples3/vmmc_2d.wl BitString=11110101011110011

# Failure examples
wolframscript -file check.wls examples3/kawasaki_1d_fail.wl MaxBitString=1111111
wolframscript -file check.wls examples3/cluster_1d_fail.wl MaxBitString=1111111 Mode=Symbolic
wolframscript -file check.wls examples3/cluster_2d.wl MaxBitString=1111111 Mode=Symbolic
```

**Options** (`key=value`, no spaces around `=`):

| Option | Default | Description |
|---|---|---|
| `MaxBitString=XXXX` | `11111111` | Test all bit strings from length 1 up through `len(XXXX)`, stopping at `XXXX` within the final length |
| `Mode=Symbolic` | `Both` | Symbolic DB check only |
| `Mode=Numerical` | `Both` | Numerical MCMC check only |
| `Mode=Both` | `Both` | Both checks |
| `NSteps=N` | `50000` | MCMC steps per tested component |
| `MaxBitDepth=N` | `20` | BFS tree depth cap per state |
| `Verbose=True` | `False` | Print per-state BFS progress |

### How MaxBitString determines which systems are tested

The checker iterates every candidate bit string up to `MaxBitString`, calling `BitsToState` to get a seed state. If the seed belongs to a previously-discovered connected component it is skipped. Otherwise BFS discovers the full component and one row is printed.

### Fast enumeration with `ValidStateIDs`

If the algorithm file defines `ValidStateIDs[maxId_]`, the checker uses it to enumerate only the integer IDs that decode to valid states, skipping all others without calling `BitsToState`. For `vmmc_2d.wl` this avoids calling `BitsToState` on the ~99% of bit strings that decode to non-square-length arrays, reducing the overhead before reaching L=3 systems from ~500K iterations to ~110K. Previously-seen states are tracked in an `Association` for O(1) lookup, so after the first connected component is found, the remaining IDs in that grid-size range are skipped instantly.

---

## Animating an algorithm

Requires Python 3 with `matplotlib` and `numpy`.

```bash
# 2D VMMC (large system, fast 2-colour mode)
wolframscript -file animate.wls examples3/vmmc_2d.wl \
  Sites=100 N=40 Steps=5000 FPS=30 Simple=1

# 2D Kawasaki torus (3×3 = 9 sites, 4 particles)
wolframscript -file animate.wls examples3/kawasaki_2d.wl \
  Sites=9 N=4 Steps=300 FPS=10 J12=1.5

# 1D Kawasaki ring (6 sites, 3 particles)
wolframscript -file animate.wls examples3/kawasaki_1d.wl \
  Sites=6 N=3 Steps=500 FPS=20

# With attractive coupling and recording every 5th step
wolframscript -file animate.wls examples3/vmmc_2d.wl \
  Sites=100 N=40 Steps=2000 RecordEvery=5 FPS=20 Simple=1 J12=-1.5
```

**Options:**

| Option | Default | Description |
|---|---|---|
| `Sites=<n>` | required | Total lattice sites (e.g. 9 for a 3×3 grid) |
| `N=<n>` | required | Number of labeled particles (types 1..N) |
| `Steps=<n>` | `200` | MCMC steps to run |
| `Beta=<f>` | from `.wl` file | Inverse temperature |
| `FPS=<f>` | `10` | Animation frame rate |
| `Simple=1` | off | Fast 2-colour mode (particles vs holes); enables matplotlib blitting |
| `RecordEvery=<n>` | `1` | Record state every nth step |
| `J<a><b>=<f>` | random | Set coupling constant |

---

## Writing your own algorithm file

Copy `template.wl` and implement the four required functions.

### Required: `energy[state_]`

Bare energy (no beta). Must be deterministic. May use symbolic variables declared in `symParams` or `DynamicSymParams`.

### Required: `Algorithm[state_]`

MCMC move. Use native Mathematica random calls (`RandomInteger[]`, `RandomReal[]`, `RandomChoice[]`). Use `MetropolisProb[dE]` for Metropolis acceptance:

```mathematica
If[RandomReal[] < MetropolisProb[dE], newState, state]
```

The checker intercepts these calls during BFS and replaces them with deterministic values from a bit tape. `MetropolisProb[dE]` expands to `Piecewise[{{1, dE<=0}, {Exp[-β·dE], dE>0}}]` with β kept symbolic.

**Supported random calls:**

| Call | Bits consumed |
|---|---|
| `RandomReal[]` | 1 (deferred until compared with `< p`) |
| `RandomInteger[{lo, hi}]` | ⌈log₂(hi−lo+1)⌉ bits |
| `RandomInteger[]` or `RandomInteger[1]` | 1 bit |
| `RandomChoice[list]` | ⌈log₂(n)⌉ bits |

### Required: `BitsToState[bits_List]`

Converts a bit string to a seed state, or returns `None` to skip it.

### Required: `numBeta`

```mathematica
numBeta = 1
```

### Optional: `ValidStateIDs[maxId_Integer]`

If defined, `check.wls` calls this instead of iterating all bit strings. Return the list of integer IDs (up to `maxId`) that `BitsToState` will accept. For the 2D lattice encoding this is the union of the ID ranges that decode to arrays of length 1, 4, 9, 16, … (perfect squares):

```mathematica
ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2 + 1] - 1, maxId]]];
      L++];
    ids]
```

### Optional: symbolic parameters

**Static** (same for every component):
```mathematica
symParams = <|"eps" -> {eps1, eps2}, "couplings" -> {J12, J13}|>
```

**Dynamic** (auto-generated per component from the discovered particle types):
```mathematica
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
        {a, types}, {b, types}]|>]
```

### Optional: `DisplayState[state_]`

If defined, replaces `ToString[state]` in the output table. For 2D grid states:
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
| `examples3/vmmc_2d.wl` | 2D VMMC, Whitelam–Geissler cluster moves | **PASS** |
| `examples3/kawasaki_1d_fail.wl` | 1D Kawasaki, wrong-sign dE | **FAIL** |
| `examples3/kawasaki_2d_fail.wl` | 2D Kawasaki, wrong-sign dE | **FAIL** |
| `examples3/cluster_1d_fail.wl` | 1D rightward cluster slide (no reverse move) | **FAIL** |
| `examples3/cluster_2d.wl` | 2D rigid cluster move (merging breaks proposal symmetry) | **FAIL** |

---

## 2D Virtual Move Monte Carlo (VMMC)

`examples3/vmmc_2d.wl` implements the Whitelam–Geissler VMMC algorithm (Whitelam & Geissler, *J. Chem. Phys.* 127, 154101, 2007) on a periodic square lattice with nearest-neighbour J couplings.

### State encoding

States are flat arrays of length L² (row-major). `State[[s]] = 0` (empty) or `k ∈ {1,…,N}` (labeled particle of type k). L is inferred as `Sqrt[Length[state]]`; only states whose array length is a perfect square are valid.

The integer ID encoding is the same bijective map used by `kawasaki_2d.wl`: each non-negative integer maps to a unique labeled-occupancy array. `BitsToState` decodes the integer ID from the bit string and checks the length is a perfect square.

### Energy

Energy is the sum of J_ab coupling constants over all unique nearest-neighbour bonds on the L×L torus. Unique bonds are enumerated by `$uniqueBonds2D[L]` (memoised): for each site s, the pairs `{s, right(s)}` and `{s, down(s)}` are sorted and deduplicated. This correctly handles L=2 (where `right(s)=left(s)`) without any special-case divisor.

### Algorithm

Each step:

1. A random occupied site is chosen as the cluster seed; a random direction d ∈ {right, left, down, up} is proposed.

2. A cluster is grown by BFS. For each occupied non-cluster neighbour q of cluster particle p, the neighbour shell tested is the union of neighbours of p, p+d (virtual forward), and p−d (virtual reverse) — this ensures bond-forming moves are correctly link-tested. For each candidate q:

   ```
   eInit = current pair energy between p and q
   eFwd  = pair energy if p is at virtual forward position
   eRev  = pair energy if p is at virtual reverse position
   wFwd  = max(1 − exp(β(eInit − eFwd)), 0)
   wRev  = max(1 − exp(β(eInit − eRev)), 0)
   ```

   Draw r1: if r1 > wFwd, skip q. Otherwise draw r2: if r2 > wRev/wFwd, the link is *frustrated* and the entire move is rejected; otherwise q joins the cluster.

3. If no frustrated link was detected, all cluster particles translate rigidly by one step in direction d.

**Hard-sphere exclusion.** Same-site virtual occupancy (`vI == qSite`) gives `$virtualPairEnergy = Infinity`, making wFwd = 1 (mandatory link attempt). The Piecewise expressions for wFwd, wRev, and the ratio include explicit `=== Infinity` guards to prevent `Exp[β(finite − Infinity)]` with symbolic β, which would generate warnings.

**Why no Metropolis gate.** The Whitelam–Geissler link probabilities enforce superdetailed balance. An additional Metropolis step based on the total energy change would double-count bond contributions and break detailed balance. `MetropolisProb[dE]` is called for checker interface compliance but its return value is not used as an acceptance gate.

**Why Piecewise (not Max).** `Max[1−exp(x), 0]` is opaque to `FullSimplify` for unknown-sign symbolic J. `Piecewise` exposes the case structure and allows `FullSimplify` to verify each sign-case algebraically.

**Why β appears in the exponents.** The checker keeps β as a free symbol. The link weights must contain β so that `T(i→j) ∝ E^{βJ}` and the Boltzmann factor `E^{−βJ}` cancel exactly for all β. Using `Exp[eInit − eFwd]` (implicit β=1) passes numerical checks (since `numBeta=1`) but fails the symbolic check.

---

## Limitations

- **State space size.** Keep components to ~20 states for fast symbolic checks; ~72 states (3×3, 2 particles) is feasible but takes longer.
- **Random call types.** `RandomVariate`, `AbsoluteTime`, `SessionTime`, and `Now` cannot be intercepted.
- **Flat-energy caveat.** If all energies are equal the numerical check passes for any ergodic algorithm; trust the symbolic check.

---

## Files

```
check.wls               Exhaustive checker (iterates bit strings up to MaxBitString)
report.wls              Single-seed report with PNG visualisation
animate.wls             Animation runner
animate_plot.py         Python/matplotlib animation display
show_report.py          Python/matplotlib report renderer
dbc_core.wl             Core library (BFS, symbolic check, MCMC check)
template.wl             Template for writing an algorithm file
examples3/
  kawasaki_1d.wl        1D Kawasaki, dynamic L, labeled particles        [PASS]
  kawasaki_2d.wl        2D Kawasaki, dynamic L, labeled particles        [PASS]
  vmmc_2d.wl            2D VMMC, Whitelam–Geissler cluster moves          [PASS]
  kawasaki_1d_fail.wl   1D Kawasaki, wrong-sign dE                       [FAIL]
  kawasaki_2d_fail.wl   2D Kawasaki, wrong-sign dE                       [FAIL]
  cluster_1d_fail.wl    1D rightward cluster slide (no reverse)          [FAIL]
  cluster_2d.wl         2D rigid cluster move (merging breaks symmetry)  [FAIL]
legacy/
  examples/             Old example format (unsupported)
  examples2/            Fixed-L examples (superseded by examples3/)
```
