(* ================================================================
   1D Kawasaki with rightward drift -- FAILS symbolic check,
   PASSES numerical check
   ================================================================

   System: L=4 ring, 3 particles, 1 hole.
   State:  sorted triple {p1,p2,p3}.
   Energy: FLAT (all energies = 0). Boltzmann distribution is
           uniform over all 4 states.

   Move:   Always hop the particle immediately to the left of
           the hole one step clockwise (into the hole).
           No randomness -- completely deterministic drift.

   Why FAILS detailed balance (symbolic check catches this):
   ---------------------------------------------------------
   The algorithm drives a persistent clockwise current:
     {1,2,3} -> {1,2,4} -> {1,3,4} -> {2,3,4} -> {1,2,3} -> ...
   T(i->next(i)) = 1, but T(next(i)->i) = 0.
   The residual T(i->j)*pi(i) - T(j->i)*pi(j) = 1 != 0 for each
   adjacent pair in the cycle.

   Why PASSES numerically:
   ----------------------
   With flat energy, Boltzmann is uniform (pi = 1/4 for each state).
   The deterministic cycle visits all four states equally often,
   so the simulated histogram matches Boltzmann exactly.
   KL divergence ~ 0 despite the detailed-balance violation.

   This demonstrates that the numerical check is NECESSARY but NOT
   SUFFICIENT: the symbolic check is required for a rigorous guarantee.
   ================================================================ *)

(* ---- System parameters ---- *)
L$df = 4

(* Flat energy: Boltzmann is uniform. No symbolic params needed. *)
energy[{p1_Integer, p2_Integer, p3_Integer}] := 0

(* ================================================================
   Algorithm: deterministic rightward drift.
   Find the hole, hop the particle to its left one step right.
   No random calls -- deterministic algorithm.
   ================================================================ *)
Algorithm[{p1_Integer, p2_Integer, p3_Integer}] :=
  Module[{occ, hole, mover},
    occ   = {p1, p2, p3};
    hole  = Complement[Range[L$df], occ][[1]];
    (* Particle immediately to the left of the hole (periodic) *)
    mover = Mod[hole - 2, L$df] + 1;
    Sort[Append[DeleteCases[occ, mover], hole]]
  ]

(* ---- Checker interface ---- *)
seedState = {1, 2, 3}
numBeta   = 1
