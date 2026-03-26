(* ================================================================
   Example: Two-particle Kawasaki with pairwise coupling -- PASSES
   ================================================================

   System: L=4 ring, 2 distinguishable particles with hard-core
           exclusion (no two particles on the same site).
   State:  sorted pair {p1, p2} with 1 <= p1 < p2 <= L.
           There are C(4,2) = 6 states:
             {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}.

   Energy: sum of bare site energies of the two occupied sites
           plus the nearest-neighbour pairwise interaction:

             E({p1,p2}) = \[Epsilon]mp_p1 + \[Epsilon]mp_p2
                        + Jmp1  (if ring-distance between p1 and p2 is 1)

           All energy parameters are symbolic; random numerical values
           are assigned automatically during the numerical MCMC check.

   Move:   pick a random particle (bit1), pick a random direction
           (bit2), propose a hop to the neighbouring site.
           If the target site is occupied -> stay (hard-core rejection).
           Otherwise apply Metropolis acceptance.

   The proposal is symmetric (each particle proposed with prob 1/2,
   each direction with prob 1/2) and Metropolis is correct, so
   detailed balance holds exactly.
   ================================================================ *)

L$mp   = 4
(* Symbolic site energies and nearest-neighbour coupling.
   Kept unassigned during the symbolic check; random numerical values
   are assigned automatically by RunFullCheck via "SysParams" -> params$mp. *)
eps$mp     = {\[Epsilon]mp1, \[Epsilon]mp2, \[Epsilon]mp3, \[Epsilon]mp4}
params$mp  = <|"L" -> L$mp, "eps" -> eps$mp, "couplings" -> {Jmp1}|>
numBeta$mp = 1

(* Ring neighbours *)
rightOf$mp[s_Integer] := Mod[s,     L$mp] + 1
leftOf$mp[s_Integer]  := Mod[s - 2, L$mp] + 1

(* Energy: site energies + nearest-neighbour coupling.
   RingDist is provided by dbc_core.wl. *)
energy$mp[{p1_Integer, p2_Integer}] :=
  eps$mp[[p1]] + eps$mp[[p2]] +
  If[RingDist[p1, p2, L$mp] == 1, Jmp1, 0]

(* ================================================================
   Algorithm: single definition for both symbolic and numeric checks.
   ================================================================ *)
KawasakiMulti[{p1_Integer, p2_Integer}] :=
  Module[
    {b1, b2, mover, other, target, dE},
    (* Pick which particle to attempt moving *)
    b1 = RandomInteger[];
    {mover, other} = If[b1 == 0, {p1, p2}, {p2, p1}];
    (* Pick direction *)
    b2     = RandomInteger[];
    target = If[b2 == 0, leftOf$mp[mover], rightOf$mp[mover]];
    (* Hard-core rejection: target occupied *)
    If[target == other,
      {p1, p2},          (* stay put *)
      (* Metropolis acceptance using full pairwise energy *)
      dE = energy$mp[Sort[{target, other}]] - energy$mp[{p1, p2}];
      If[RandomReal[] < MetropolisProb[dE],
        Sort[{target, other}],    (* accept: new sorted state *)
        {p1, p2}]                 (* reject: keep current state *)
    ]
  ]
