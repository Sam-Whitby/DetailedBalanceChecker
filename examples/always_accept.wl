(* ================================================================
   Example 2: Kawasaki with no Metropolis acceptance -- FAILS
   ================================================================

   Same ring system.  The proposed swap is always accepted regardless
   of the energy change.  This samples the UNIFORM distribution, not
   the Boltzmann distribution, so detailed balance fails whenever
   site energies differ.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

allStates$aa = allStates$rk
symEnergy$aa = symEnergy$rk
numEnergy$aa = numEnergy$rk

(* Symbolic: always accept -- returns the new site with prob 1 *)
AlwaysAcceptSym[state_Integer, readBit_] := Module[
  {dir, nbr},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  nbr    (* single state returned; probability 1 *)
]

(* Numerical: identical logic *)
AlwaysAcceptNum[state_Integer, readBit_] := Module[
  {dir, nbr},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  nbr
]
