(* ================================================================
   Example 3: Metropolis with inverted acceptance sign -- FAILS
   ================================================================

   Same L=3 ring system but the acceptance uses +beta*deltaE
   instead of -beta*deltaE.  This preferentially moves to HIGHER
   energy states, sampling exp(+beta*E) rather than exp(-beta*E).
   Detailed balance fails because the acceptance ratio is inverted
   relative to the target Boltzmann distribution.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

(* Buggy acceptance: wrong sign in exponent *)
WrongSignProb[dE_] :=
  Piecewise[{{1, dE >= 0}, {Exp[\[Beta] * dE], dE < 0}}]

WrongSign[state_Integer] := Module[
  {dir, nbr, dE},
  dir = RandomInteger[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = energy$rk[nbr] - energy$rk[state];
  If[RandomReal[] < WrongSignProb[dE], nbr, state]
]
