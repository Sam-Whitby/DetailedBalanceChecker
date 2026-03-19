(* ================================================================
   Example: Kawasaki written with native random calls -- PASSES
   ================================================================

   This is functionally identical to kawasaki_new.wl but written
   in "natural" Mathematica style using RandomInteger[] and
   RandomReal[] instead of readBit[]/acceptTest[].

   The checker AUTOMATICALLY intercepts these calls:
     RandomInteger[1]   ->  readBit[]       (0 or 1, weight 1/2 each)
     RandomReal[] < p   ->  acceptTest[p]   (weight p or 1-p)

   so the same BFS and DB verification machinery applies without
   any modification to the algorithm.

   The algorithm correctly implements Metropolis on an L=3 ring and
   therefore passes the detailed balance check.
   ================================================================ *)

L$rn       = 3
eps$rn     = {0, 1, 1/2}   (* exact rationals, no beta *)
numBeta$rn = 3/2

energy$rn[s_Integer] := eps$rn[[s]]

(* ================================================================
   Algorithm written with native Mathematica random calls.
   readBit_ and acceptTest_ appear in the signature but are NOT
   used; the checker intercepts RandomInteger/RandomReal instead.
   ================================================================ *)
KawasakiNative[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE, p},
  dir = RandomInteger[1];    (* 0 or 1: intercepted -> readBit[] *)
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rn] + 1;
  dE  = energy$rn[nbr] - energy$rn[state];
  p   = MetropolisProb[dE];
  (* RandomReal[] < p: intercepted -> acceptTest[p] *)
  If[RandomReal[] < p, nbr, state]
]
