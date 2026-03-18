# DetailedBalanceChecker

Symbolic and numerical checker for the detailed balance condition in MCMC algorithms, written in Mathematica and runnable from the terminal via `wolframscript`.

## What it does

Given an MCMC algorithm and an energy function, the checker:

1. **Symbolically** builds the complete decision tree for every starting state by exhaustively enumerating all possible bit sequences the algorithm might consume during execution. It constructs a symbolic transition matrix and verifies that detailed balance holds **exactly**, using Mathematica's `Simplify`.

2. **Numerically** runs the algorithm as a true MCMC chain, then compares the visited-state histogram to the Boltzmann distribution and reports the KL divergence.

### Why bits?

All discrete randomness (site selection, move-type choices, etc.) is drawn from a binary stream. Different execution paths may consume different numbers of bits. The checker uses BFS over bit sequences, extending a path by one bit whenever the algorithm calls `readBit[]` beyond the currently supplied sequence. Each leaf in the resulting tree has probability $(1/2)^k$ where $k$ is the number of bits consumed on that path.

Continuous acceptance steps (Metropolis) are handled by returning a symbolic `Piecewise` probability directly from the algorithm function, rather than consuming bits. This keeps the transition matrix entries as exact symbolic expressions involving $e^{-\beta \Delta E}$, which Mathematica can simplify exactly.

### Limitation: the Halting problem

If the algorithm never returns for some bit sequence, the BFS will eventually hit `MaxBitDepth` and that path will be excluded with a warning. A `TimeLimit` (seconds per state) provides a second safety net.

## Algorithm interface

Provide **two** algorithm functions with the same signature:

```mathematica
symAlg[state_, readBit_]   (* for the symbolic check *)
numAlg[state_, readBit_]   (* for the numerical MCMC check *)
```

- `readBit[]` returns `0` or `1` each time it is called.
- Each function returns **either**:
  - A single new state — if all randomness came from bits.
  - `{{p1, s1}, {p2, s2}, ...}` — explicit probability-weighted outcomes (probabilities sum to 1). Use this for the Metropolis acceptance step.

**`symAlg`**: energy differences should be symbolic (e.g. `eps2 - eps1`). Use `MetropolisProb[deltaE]` which produces a `Piecewise` expression symbolic in `β`. `β` must be unassigned when the symbolic check runs.

**`numAlg`**: all quantities should be fully numeric. Compute `dE` from a numeric energy function that already has `beta` folded in. Use `N@If[dE <= 0, 1, Exp[-dE]]` for Metropolis acceptance.

**`symEnergy[state]`**: bare energy `E(s)`, symbolic in coupling constants, **without** `β`. The checker inserts `β` via `Exp[-β * symEnergy[s]]`.

**`numEnergy[state]`**: fully numeric, **with** `beta` folded in (e.g. `numBeta * E_numeric(s)`). Used to compute Boltzmann weights as `Exp[-numEnergy[s]]`.

## Files

```
dbc_core.wl            Core library (load this)
run_checks.wls         Main runner script
examples/
  ring_kawasaki.wl     Ring system, Kawasaki + Metropolis  (PASS)
  always_accept.wl     Ring system, always accept          (FAIL)
  wrong_sign.wl        Ring system, wrong acceptance sign  (FAIL)
```

## Usage

### Clone and run

```bash
git clone https://github.com/Sam-Whitby/DetailedBalanceChecker.git ~/Desktop/DetailedBalanceChecker
cd ~/Desktop/DetailedBalanceChecker
wolframscript -file run_checks.wls
```

### Use in your own script

```mathematica
Get["path/to/dbc_core.wl"]

allStates  = {1, 2, 3}
numBeta    = 1.5
numEps     = {0., 1., 0.5}
epsSymbols = {eps1, eps2, eps3}   (* unassigned symbols *)

symEnergy[s_] := epsSymbols[[s]]              (* bare E(s), no beta *)
numEnergy[s_] := numBeta * numEps[[s]]        (* beta * E(s), numeric *)

(* Symbolic algorithm -- beta stays unassigned *)
mySymAlg[state_, readBit_] := Module[{dir, nbr, dE, p},
  dir = readBit[];
  nbr = ...;
  dE  = symEnergy[nbr] - symEnergy[state];    (* eps_nbr - eps_state *)
  p   = MetropolisProb[dE];                   (* Piecewise in beta, symbolic *)
  {{p, nbr}, {1 - p, state}}
]

(* Numeric algorithm -- beta already in dE *)
myNumAlg[state_, readBit_] := Module[{dir, nbr, dE, p},
  dir = readBit[];
  nbr = ...;
  dE  = numEnergy[nbr] - numEnergy[state];    (* numeric, beta-scaled *)
  p   = N@If[dE <= 0, 1, Exp[-dE]];
  {{p, nbr}, {1 - p, state}}
]

RunFullCheck[allStates, mySymAlg, myNumAlg, symEnergy, numEnergy,
  "SystemName"  -> "My system",
  "MaxBitDepth" -> 15,
  "NSteps"      -> 100000
]
```

## Practical limits

The symbolic check is exact but expensive. The number of states and the depth of the decision tree grow quickly:

- **States**: $\binom{L}{n}$ for $n$ particles on $L$ sites. Keep $L \leq 5$ for reasonable run times.
- **Bit depth**: Each level doubles the number of paths. `MaxBitDepth=20` covers most simple algorithms.
- **Simplify cost**: `Simplify` on expressions involving multiple `Piecewise` and exponentials can be slow. The check is designed for small illustrative examples, not production-scale systems.
