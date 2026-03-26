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

# Kawasaki failure examples (wrong-sign dE)
wolframscript -file check.wls examples3/kawasaki_1d_fail.wl MaxBitString=1111111
wolframscript -file check.wls examples3/kawasaki_2d_fail.wl MaxBitString=1111111 Mode=Symbolic

# 1D rightward cluster move (FAILS -- no reverse move)
wolframscript -file check.wls examples3/cluster_1d_fail.wl MaxBitString=1111111 Mode=Symbolic

# 2D rigid cluster move (FAILS -- cluster merging breaks proposal symmetry)
wolframscript -file check.wls examples3/cluster_2d.wl MaxBitString=1111111 Mode=Symbolic
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

Requires Python 3 with `matplotlib` and `numpy`.

```bash
# 1D Kawasaki ring (6 sites, 3 particles)
wolframscript -file animate.wls examples3/kawasaki_1d.wl \
  Sites=6 N=3 Steps=500 FPS=20

# 2D Kawasaki torus (3x3 = 9 sites, 4 particles)
wolframscript -file animate.wls examples3/kawasaki_2d.wl \
  Sites=9 N=4 Steps=300 FPS=10 J12=1.5

# 2D VMMC (large system, fast 2-colour mode)
wolframscript -file animate.wls examples3/vmmc_2d.wl \
  Sites=100 N=40 Steps=5000 FPS=30 Simple=1

# 1D cluster move (also works for fail examples)
wolframscript -file animate.wls examples3/cluster_1d_fail.wl \
  Sites=8 N=4 Steps=400 FPS=15

# 2D cluster move
wolframscript -file animate.wls examples3/cluster_2d.wl \
  Sites=9 N=3 Steps=400 FPS=15

# Specify coupling constants explicitly
wolframscript -file animate.wls examples3/kawasaki_1d.wl \
  Sites=6 N=3 J12=1.0 J13=0.5 J23=2.0

# Record every 5th step for long runs
wolframscript -file animate.wls examples3/kawasaki_2d.wl \
  Sites=9 N=4 Steps=2000 RecordEvery=5 FPS=20
```

The animation opens in a separate matplotlib window with:
- **Left panel**: colour-coded imshow of the lattice (dark = empty, each particle type has a distinct colour). For 2D algorithms the grid is shown as L×L; for 1D as a single row.
- **Right panel**: system energy plotted against step number, growing in real time.

**Fast rendering with `Simple=1`.** By default each particle type gets its own colour, which requires a full figure redraw per frame and limits the achievable frame rate. Passing `Simple=1` switches to a two-colour mode — light grey for empty sites, blue for any occupied site — and enables matplotlib's blitting path (`blit=True`), so only the lattice image is redrawn each frame. Use this when you want high frame rates on large lattices, e.g. `FPS=30 Simple=1` for real-time VMMC visualisation.

**Options:**

| Option | Default | Description |
|---|---|---|
| `Sites=<n>` | required | Total lattice sites (e.g. 9 for a 3×3 grid) |
| `N=<n>` | required | Number of labeled particles (types 1..N) |
| `Steps=<n>` | `200` | MCMC steps to run |
| `Beta=<f>` | from `.wl` file | Inverse temperature |
| `FPS=<f>` | `10` | Animation frame rate in frames per second |
| `Simple=1` | off | Fast 2-colour mode: holes vs particles (enables blitting) |
| `RecordEvery=<n>` | `1` | Record state every nth step |
| `J<a><b>=<f>` | random | Set coupling constant (random otherwise) |

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
| `examples3/vmmc_2d.wl` | 2D VMMC, Whitelam–Geissler cluster moves | **PASS** |
| `examples3/kawasaki_1d_fail.wl` | 1D Kawasaki, wrong-sign dE | **FAIL** |
| `examples3/kawasaki_2d_fail.wl` | 2D Kawasaki, wrong-sign dE | **FAIL** |
| `examples3/cluster_1d_fail.wl` | 1D rightward cluster slide | **FAIL** |
| `examples3/cluster_2d.wl` | 2D rigid cluster move | **FAIL** |

### Kawasaki failure: wrong-sign dE

```mathematica
(* Correct *)    dE = energy[newState] - energy[state]
(* Buggy   *)    dE = energy[state] - energy[newState]   (* sign reversed *)
```

With the wrong sign, `MetropolisProb` receives `-ΔE`: uphill moves are always accepted, downhill moves are penalised. The algorithm preferentially visits high-energy states, inverting the Boltzmann distribution.

**Note on bond-selection bias.** An asymmetric bond proposal — e.g. `Mod[RandomInteger[{0,L}],L]` making bond 0 twice as likely — does *not* break detailed balance for Kawasaki. Both forward and reverse transitions of the same bond have the same proposal probability, which cancels in the DB condition. Only the acceptance step matters.

### 1D rightward cluster slide

Finds all maximal connected runs of particles and slides one randomly chosen run one step clockwise if the site to its right is empty. The reverse (leftward) move is never proposed, so T(j→i) = 0 wherever T(i→j) > 0. The symbolic checker detects this immediately for any component containing a rightward-possible transition.

### 2D rigid cluster move

Picks a connected cluster at random and attempts to translate it one step N/E/S/W if all target sites are empty. This fails because **cluster merging changes the number of available clusters** between the two states:

- State i has n clusters → picks cluster C with probability 1/n
- Cluster C moves east, merging with cluster D → state j has n-1 clusters
- Reverse: state j picks the merged cluster CD with probability 1/(n-1)

So T(i→j) ∝ 1/n but T(j→i) ∝ 1/(n-1), and the ratio 1/n ≠ 1/(n-1) cannot be compensated by the Metropolis acceptance alone. The checker detects this on the 2×2 torus (tested by `MaxBitString=1111111`) where merging transitions exist.

A correct cluster-move algorithm would use a Rosenbluth-weight acceptance criterion that explicitly accounts for the number of cluster choices in the forward and reverse states.

### 2D Virtual Move Monte Carlo (VMMC)

`examples3/vmmc_2d.wl` implements the Whitelam–Geissler VMMC algorithm on a periodic square lattice with nearest-neighbour J couplings (Whitelam & Geissler, *J. Chem. Phys.* 127, 154101, 2007).

**Algorithm summary.** Each step:
1. A random occupied site is chosen as the cluster seed, and a random lattice direction d ∈ {right, left, down, up} is proposed.
2. A cluster is grown by BFS. For each occupied non-cluster neighbour q of cluster particle p, forward and reverse link weights are computed:

   ```
   wFwd = max(1 − exp(eInit − eFwd), 0)
   wRev = max(1 − exp(eInit − eRev), 0)
   ```

   where `eInit` is the current pair energy, `eFwd` is the pair energy if p makes the forward move, and `eRev` is the pair energy if p makes the reverse move. A random number r1 is drawn: if r1 > wFwd the neighbour is skipped; otherwise a second draw r2 determines whether the link is *frustrated* (r2 > wRev/wFwd, reject the whole move) or genuine (q joins the cluster).

3. If no frustrated link was found, all cluster particles translate rigidly by one lattice step in direction d.

**Why no Metropolis gate.** The Whitelam–Geissler link probabilities enforce *superdetailed balance* directly: the cluster-growth mechanism already produces the correct acceptance ratio between forward and reverse moves. Applying an additional Metropolis accept/reject step based on the total energy change would double-count bond contributions and break detailed balance. `MetropolisProb[dE]` is called in the algorithm for interface compliance with the checker, but its return value is not used as an acceptance gate.

**Hard-sphere exclusion.** Same-site occupancy is treated as infinite repulsive energy (wFwd = 1, mandatory link attempt). On small lattices (L = 2), the same site can appear as its own neighbour in both the right and left (or up and down) directions; `DeleteDuplicates` is applied to the neighbour list to prevent the link probability from being applied twice to the same pair.

**Checking VMMC:**

```bash
# Symbolic check on small systems
wolframscript -file check.wls examples3/vmmc_2d.wl MaxBitString=1111111 Mode=Symbolic

# Both symbolic and numerical
wolframscript -file check.wls examples3/vmmc_2d.wl MaxBitString=1111111
```

**Animating VMMC** (use `Simple=1` for fast rendering at high frame rates):

```bash
# 10×10 lattice, 40 particles, 1000 steps at 30 fps
wolframscript -file animate.wls examples3/vmmc_2d.wl \
  Sites=100 N=40 Steps=1000 FPS=30 Simple=1

# With attractive coupling between particle types 1 and 2
wolframscript -file animate.wls examples3/vmmc_2d.wl \
  Sites=100 N=40 Steps=2000 FPS=30 Simple=1 J12=-1.5
```

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
animate.wls             Animation runner (collects data, calls Python)
animate_plot.py         Python/matplotlib animation display
template.wl             Template for writing an algorithm file
dbc_core.wl             Core library (BFS, symbolic check, MCMC check)
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
