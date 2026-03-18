(* ================================================================
   Example 3: Metropolis with inverted acceptance sign -- FAILS
   ================================================================

   Same ring system but the acceptance uses +beta*deltaE instead
   of -beta*deltaE.  This preferentially moves to HIGHER energy
   states, sampling exp(+beta*E) rather than exp(-beta*E).
   Detailed balance fails because the acceptance ratio is inverted
   relative to the target Boltzmann distribution.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

allStates$ws = allStates$rk
symEnergy$ws = symEnergy$rk
numEnergy$ws = numEnergy$rk

(* Buggy symbolic acceptance: wrong sign in exponent *)
WrongSignProb[deltaE_] :=
  Piecewise[{{1, deltaE >= 0}, {Exp[\[Beta] * deltaE], deltaE < 0}}]

WrongSignSym[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = symEnergy$ws[nbr] - symEnergy$ws[state];
  p   = WrongSignProb[dE];
  {{p, nbr}, {1 - p, state}}
]

(* Buggy numerical acceptance: wrong sign *)
WrongSignNum[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = numEnergy$ws[nbr] - numEnergy$ws[state];   (* beta-scaled *)
  p   = N@If[dE >= 0, 1, Exp[dE]];                 (* wrong sign: accept uphill *)
  {{p, nbr}, {1 - p, state}}
]
