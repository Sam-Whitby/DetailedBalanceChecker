(* ================================================================
   Example 3: Metropolis with wrong sign in acceptance -- FAILS
   ================================================================

   Same ring system again, but the acceptance probability uses
   +beta*deltaE instead of -beta*deltaE.  This favours moves that
   INCREASE the energy, sampling exp(+beta*E) rather than exp(-beta*E).

   Detailed balance fails because the acceptance ratio is inverted
   relative to the Boltzmann weights.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

allStates$ws  = allStates$rk
symEnergy$ws  = symEnergy$rk
numEnergy$ws  = numEnergy$rk

(* Buggy acceptance: flipped sign *)
WrongSignProb[deltaE_] :=
  Piecewise[{{1, deltaE >= 0}, {Exp[\[Beta] deltaE], deltaE < 0}}]

WrongSign[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = symEnergy$rk[nbr] - symEnergy$rk[state];
  p   = WrongSignProb[dE];
  {{p, nbr}, {1 - p, state}}
]
