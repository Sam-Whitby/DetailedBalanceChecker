(* ================================================================
   Example 1: Kawasaki dynamics on a ring -- PASSES detailed balance
   ================================================================

   System: L=3 sites on a ring, 1 particle.
   State:  integer 1..L (which site the particle occupies).
   Energy: E(i) = eps_i (exact rationals, no beta).
   Move:   read 1 bit -> choose left (0) or right (1) neighbour,
           then Metropolis acceptance.

   The proposal is symmetric and Metropolis is correct, so
   detailed balance holds exactly.

   Single algorithm used for both the symbolic DB check (beta
   left unassigned) and the numerical MCMC (Block[{beta=numBeta}]).
   ================================================================ *)

(* ---- System parameters ---- *)
L$rk       = 3
eps$rk     = {0, 1, 1/2}   (* exact rationals -- no beta *)
numBeta$rk = 3/2

energy$rk[s_Integer] := eps$rk[[s]]

(* ================================================================
   Algorithm: single definition works symbolically and numerically.
   Uses MetropolisProb[dE] which keeps beta symbolic during the
   DB check, and evaluates numerically via Block[{beta=numBeta}]
   during the MCMC run.
   ================================================================ *)
KawasakiRing[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = energy$rk[nbr] - energy$rk[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]
