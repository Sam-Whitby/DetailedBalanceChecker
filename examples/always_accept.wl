(* ================================================================
   Example 2: Kawasaki with no Metropolis acceptance -- FAILS
   ================================================================

   Same ring system as ring_kawasaki.wl but the move always accepts
   the proposed swap, ignoring the energy difference.

   This samples the UNIFORM distribution over states rather than
   the Boltzmann distribution, so detailed balance fails whenever
   the site energies differ.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

(* Reuse same state space and energy functions *)
allStates$aa  = allStates$rk
symEnergy$aa  = symEnergy$rk
numEnergy$aa  = numEnergy$rk  (* same numeric energies as the good example *)

(* ---- Algorithm: always accept (no Metropolis) ---- *)
AlwaysAccept[state_Integer, readBit_] := Module[
  {dir, nbr},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  nbr   (* just return the new state unconditionally *)
]
