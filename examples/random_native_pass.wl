(* ================================================================
   Example: Kawasaki written with native random calls -- PASSES
   ================================================================

   Functionally identical to kawasaki_new.wl, using RandomInteger[]
   and RandomReal[] natively.  All algorithms now use this style.

   The checker intercepts these calls during BFS:
     RandomInteger[]   ->  0 or 1 from the bit tape (weight 1/2 each)
     RandomReal[] < p  ->  accept with probability p  (weight p or 1-p)

   ================================================================ *)

L$rn       = 3
(* Symbolic site energies -- kept unassigned during the symbolic check.
   Random numerical values are assigned automatically by RunFullCheck
   when "SysParams" -> params$rn is supplied. *)
eps$rn     = {\[Epsilon]rn1, \[Epsilon]rn2, \[Epsilon]rn3}
params$rn  = <|"L" -> L$rn, "eps" -> eps$rn|>
numBeta$rn = 3/2

energy$rn[s_Integer] := eps$rn[[s]]

KawasakiNative[state_Integer] := Module[
  {dir, nbr, dE, p},
  dir = RandomInteger[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rn] + 1;
  dE  = energy$rn[nbr] - energy$rn[state];
  p   = MetropolisProb[dE];
  If[RandomReal[] < p, nbr, state]
]
