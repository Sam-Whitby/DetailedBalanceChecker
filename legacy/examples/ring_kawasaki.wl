(* ================================================================
   Example 1: Kawasaki dynamics on a ring -- PASSES detailed balance
   ================================================================

   System: L=3 sites on a ring, 1 particle.
   State:  integer 1..L (which site the particle occupies).
   Energy: E(i) = eps_i (exact rationals, no beta).
   Move:   RandomInteger[] chooses left (0) or right (1) neighbour,
           then Metropolis acceptance.

   The proposal is symmetric and Metropolis is correct, so
   detailed balance holds exactly.
   ================================================================ *)

(* ---- System parameters ---- *)
L$rk       = 3
(* Symbolic site energies -- kept unassigned during the symbolic check.
   Random numerical values are assigned automatically by RunFullCheck
   when "SysParams" -> params$rk is supplied. *)
eps$rk     = {\[Epsilon]rk1, \[Epsilon]rk2, \[Epsilon]rk3}
params$rk  = <|"L" -> L$rk, "eps" -> eps$rk|>
numBeta$rk = 3/2

energy$rk[s_Integer] := eps$rk[[s]]

(* ================================================================
   Algorithm: uses native Mathematica random calls.
   RandomInteger[] picks direction (0=left, 1=right).
   RandomReal[] < MetropolisProb[dE] is the Metropolis acceptance.
   During the symbolic check \[Beta] is unassigned; during MCMC
   Block[{\[Beta]=numBeta}] makes it numeric.
   ================================================================ *)
KawasakiRing[state_Integer] := Module[
  {dir, nbr, dE},
  dir = RandomInteger[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = energy$rk[nbr] - energy$rk[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
]
