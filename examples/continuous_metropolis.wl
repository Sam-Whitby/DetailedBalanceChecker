(* ================================================================
   Example: Metropolis on a discretised continuous ring -- PASSES
   ================================================================

   System: a particle on a circle discretised to N=4 equally-spaced
           rational positions x in {0, 1/4, 1/2, 3/4}.
   State:  a rational number (not an integer).
   Energy: E(x) = 1 - Cos[2*Pi*x]  (cosine well, minimum at x=0).
           Evaluated at the four sites:
             E(0)   = 0
             E(1/4) = 1
             E(1/2) = 2
             E(3/4) = 1
   Move:   propose x +/- 1/4 (wrapping modulo 1), then Metropolis.

   This example demonstrates:
   1. States need not be integers -- any Mathematica expression works.
   2. The energy function can be an analytic formula applied to the
      state rather than a look-up table.
   3. The checker discovers {0, 1/4, 1/2, 3/4} automatically from
      seed 0 and handles rational arithmetic throughout.

   Detailed balance holds because the proposal is symmetric (each
   direction proposed with prob 1/2) and Metropolis is correct.
   ================================================================ *)

step$cm    = 1/4       (* discretisation step size *)
numBeta$cm = 1

(* Bare energy as an analytic function of the rational position *)
energy$cm[x_] := 1 - Cos[2 * Pi * x]

(* Ring arithmetic: wrap position into [0, 1) *)
nextPos$cm[x_] := Mod[x + step$cm, 1]
prevPos$cm[x_] := Mod[x - step$cm, 1]

(* ================================================================
   Algorithm: propose +-step on the discretised ring, Metropolis.
   ================================================================ *)
ContinuousMetropolis[x_, readBit_, acceptTest_] := Module[
  {b, xNew, dE},
  b    = readBit[];
  xNew = If[b == 0, prevPos$cm[x], nextPos$cm[x]];
  dE   = energy$cm[xNew] - energy$cm[x];
  If[acceptTest[MetropolisProb[dE]] == 1, xNew, x]
]
