(* ================================================================
   Example: Cyclic drift -- FAILS detailed balance, PASSES numerically
   ================================================================

   System: L=4 ring, 1 particle, ALL site energies equal (eps = 0).
   Move:   always hop one step clockwise; no acceptance test.

   Why it FAILS detailed balance
   --------------------------------
   T(i -> right(i)) = 1  for all i,  but
   T(right(i) -> i) = 0,
   so T(i->j)*pi(i) - T(j->i)*pi(j) = pi(i) != 0 for every
   adjacent pair.  The algorithm drives a persistent clockwise
   current -- it is a non-equilibrium steady state, not a
   reversible chain.

   Why it PASSES the numerical test
   -----------------------------------
   With all energies equal the Boltzmann distribution is uniform:
     pi_B(i) = 1/4  for i = 1,2,3,4.
   The cyclic chain visits 1->2->3->4->1->... so every state is
   visited with equal frequency 1/4, matching pi_B exactly.

   This example demonstrates that the NUMERICAL check is necessary
   but NOT sufficient: it can pass even when detailed balance fails.
   The SYMBOLIC check is required for a rigorous guarantee.
   ================================================================ *)

L$cd       = 4
eps$cd     = {0, 0, 0, 0}   (* uniform energy -- Boltzmann is flat *)
numBeta$cd = 1               (* beta irrelevant since all dE = 0 *)

energy$cd[s_Integer] := eps$cd[[s]]

(* ================================================================
   Algorithm: always hop clockwise -- deterministic, no random bits.
   The BFS produces a single leaf per state: {} -> right(state).
   ================================================================ *)
CyclicDrift[state_Integer, readBit_, acceptTest_] :=
  Mod[state, L$cd] + 1   (* always move to next site clockwise *)
