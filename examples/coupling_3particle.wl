(* ================================================================
   Example: Three-particle Kawasaki with pairwise coupling -- PASSES
   ================================================================

   System: L=4 ring, 3 distinguishable particles with hard-core
           exclusion.  With 3 particles on 4 sites there are
           C(4,3) = 4 states:
             {1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}.

   Energy: site energies + pairwise interactions summed over all
           particle pairs.  Two coupling distances are modelled:

             E({p1,p2,p3}) = sum_i \[Epsilon]c3_i
                           + Jc31  for each pair at ring-distance 1 (NN)
                           + Jc32  for each pair at ring-distance 2 (NNN)

           All energy parameters are symbolic; random numerical values
           are assigned automatically during the numerical MCMC check.

   Move:   pick one of the 3 particles using RandomInteger[{0,2}]
           (rejection sampling: reads 2 bits, discards value 3),
           pick a random direction (1 bit), propose a hop.
           Hard-core rejection if target is occupied.
           Metropolis acceptance otherwise.

   The proposal is symmetric (each particle and direction equally
   likely) and Metropolis is correct, so detailed balance holds.
   ================================================================ *)

L$c3   = 4
(* Symbolic site energies and coupling strengths.
   Jc31 = nearest-neighbour coupling (ring-distance 1).
   Jc32 = next-nearest-neighbour coupling (ring-distance 2).
   Kept unassigned during the symbolic check; random values are
   assigned automatically by RunFullCheck via "SysParams" -> params$c3. *)
eps$c3    = {\[Epsilon]c31, \[Epsilon]c32, \[Epsilon]c33, \[Epsilon]c34}
params$c3 = <|"L" -> L$c3, "eps" -> eps$c3, "couplings" -> {Jc31, Jc32}|>
numBeta$c3 = 1

(* Ring neighbours *)
rightOf$c3[s_Integer] := Mod[s,     L$c3] + 1
leftOf$c3[s_Integer]  := Mod[s - 2, L$c3] + 1

(* Energy: site energies + pairwise interactions at distances 1 and 2.
   RingDist is provided by dbc_core.wl. *)
energy$c3[{p1_Integer, p2_Integer, p3_Integer}] :=
  eps$c3[[p1]] + eps$c3[[p2]] + eps$c3[[p3]] +
  (* Nearest-neighbour coupling (d=1) *)
  Jc31 * (If[RingDist[p1, p2, L$c3] == 1, 1, 0] +
          If[RingDist[p1, p3, L$c3] == 1, 1, 0] +
          If[RingDist[p2, p3, L$c3] == 1, 1, 0]) +
  (* Next-nearest-neighbour coupling (d=2) *)
  Jc32 * (If[RingDist[p1, p2, L$c3] == 2, 1, 0] +
          If[RingDist[p1, p3, L$c3] == 2, 1, 0] +
          If[RingDist[p2, p3, L$c3] == 2, 1, 0])

(* ================================================================
   Algorithm: pick a particle, pick a direction, Metropolis.
   RandomInteger[{0,2}] selects one of the 3 particles using the
   checker's built-in rejection sampling (2 bits; value 3 discarded).
   ================================================================ *)
KawasakiCoupled3[{p1_Integer, p2_Integer, p3_Integer}] :=
  Module[
    {b, dir, parts, mover, others, target, newState, dE},
    parts = {p1, p2, p3};
    (* Pick which particle to move: 0 -> p1, 1 -> p2, 2 -> p3 *)
    b     = RandomInteger[{0, 2}];
    mover  = parts[[b + 1]];
    others = Delete[parts, b + 1];
    (* Pick direction *)
    dir    = RandomInteger[];
    target = If[dir == 0, leftOf$c3[mover], rightOf$c3[mover]];
    (* Hard-core rejection: target already occupied *)
    If[MemberQ[others, target],
      {p1, p2, p3},
      (* Metropolis acceptance using full pairwise energy *)
      newState = Sort[Join[{target}, others]];
      dE = energy$c3[newState] - energy$c3[Sort[parts]];
      If[RandomReal[] < MetropolisProb[dE],
        newState,
        Sort[parts]]
    ]
  ]
