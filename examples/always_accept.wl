(* ================================================================
   Example 2: Kawasaki with no Metropolis acceptance -- FAILS
   ================================================================

   Same L=3 ring system.  The proposed hop is always accepted
   regardless of the energy change, sampling the UNIFORM
   distribution rather than the Boltzmann distribution.
   Detailed balance fails whenever site energies differ.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "ring_kawasaki.wl"]

(* Always-accept: same proposal as KawasakiRing but ignores energy *)
AlwaysAccept[state_Integer] := Module[
  {dir, nbr},
  dir = RandomInteger[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  nbr   (* always move, no acceptance test *)
]
