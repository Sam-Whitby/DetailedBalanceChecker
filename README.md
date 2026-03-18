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

```mathematica
myAlg[state_, readBit_] := ...
```

- `readBit[]` returns `0` or `1` each time it is called (discrete random choices).
- The function returns **either**:
  - A single new state — if all randomness came from bits.
  - A list `{{p1, s1}, {p2, s2}, ...}` — probabilities and new states (for Metropolis acceptance etc.). Probabilities must sum to 1.

Use `MetropolisProb[deltaE]` from the library to produce the standard Metropolis acceptance probability as a symbolic `Piecewise` expression in `β`.

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

(* Define your system *)
allStates = {1, 2, 3}
symEnergy[s_] := eps[[s]]        (* symbolic in beta, J, etc. *)
numEnergy[s_] := 1.5 * neps[[s]] (* fully numeric *)

myAlg[state_, readBit_] := Module[{dir, nbr, dE, p},
  dir = readBit[];               (* consume 1 bit *)
  nbr = ...; dE = symEnergy[nbr] - symEnergy[state];
  p = MetropolisProb[dE];        (* Piecewise Metropolis probability *)
  {{p, nbr}, {1 - p, state}}
]

RunFullCheck[allStates, myAlg, symEnergy, numEnergy,
  "SystemName" -> "My system",
  "MaxBitDepth" -> 15,
  "NSteps" -> 100000
]
```

## Practical limits

The symbolic check is exact but expensive. The number of states and the depth of the decision tree grow quickly:

- **States**: $\binom{L}{n}$ for $n$ particles on $L$ sites. Keep $L \leq 5$ for reasonable run times.
- **Bit depth**: Each level doubles the number of paths. `MaxBitDepth=20` covers most simple algorithms.
- **Simplify cost**: `Simplify` on expressions involving multiple `Piecewise` and exponentials can be slow. The check is designed for small illustrative examples, not production-scale systems.
