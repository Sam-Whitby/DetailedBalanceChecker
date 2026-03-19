(* ================================================================
   Example: Two-particle Kawasaki on a ring -- PASSES
   ================================================================

   System: L=4 ring, 2 distinguishable particles with hard-core
           exclusion (no two particles on the same site).
   State:  sorted pair {p1, p2} with 1 <= p1 < p2 <= L.
           There are C(4,2) = 6 states:
             {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}.

   Energy: sum of bare site energies of the two occupied sites.
           energy[{p1,p2}] = eps1 + eps2  (no beta).

   Move:   pick a random particle (bit1), pick a random direction
           (bit2), propose a hop to the neighbouring site.
           If the target site is occupied -> stay (hard-core rejection).
           Otherwise apply Metropolis acceptance.

   The proposal is symmetric (each particle proposed with prob 1/2,
   each direction with prob 1/2) and Metropolis is correct, so
   detailed balance holds exactly.

   This example demonstrates that the checker works for multi-particle
   systems; the state is a list, not a single integer.
   ================================================================ *)

L$mp   = 4
eps$mp = {0, 1, 3, 2}   (* bare site energies, no beta *)
numBeta$mp = 1

energy$mp[{p1_Integer, p2_Integer}] :=
  eps$mp[[p1]] + eps$mp[[p2]]

(* Ring neighbours *)
rightOf$mp[s_Integer] := Mod[s,     L$mp] + 1
leftOf$mp[s_Integer]  := Mod[s - 2, L$mp] + 1

(* ================================================================
   Algorithm: single definition for both symbolic and numeric checks.
   ================================================================ *)
KawasakiMulti[{p1_Integer, p2_Integer}, readBit_, acceptTest_] :=
  Module[
    {b1, b2, mover, other, target, dE},
    (* Pick which particle to attempt moving *)
    b1 = readBit[];
    {mover, other} = If[b1 == 0, {p1, p2}, {p2, p1}];
    (* Pick direction *)
    b2     = readBit[];
    target = If[b2 == 0, leftOf$mp[mover], rightOf$mp[mover]];
    (* Hard-core rejection: target occupied *)
    If[target == other,
      {p1, p2},          (* stay put *)
      (* Metropolis acceptance *)
      dE = eps$mp[[target]] - eps$mp[[mover]];
      If[acceptTest[MetropolisProb[dE]] == 1,
        Sort[{target, other}],    (* accept: new sorted state *)
        {p1, p2}]                 (* reject: keep current state *)
    ]
  ]
