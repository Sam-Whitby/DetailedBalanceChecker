(* ================================================================
   Example: non-power-of-2 RandomInteger via rejection sampling
   ================================================================
   L=3 ring with Metropolis acceptance.  The neighbour is chosen
   by drawing RandomInteger[{1,3}] and skipping the current site,
   giving a uniform proposal over the two neighbours.

   RandomInteger[{1,3}] has 3 outcomes -- not a power of 2.
   The checker reads k=2 bits (4 outcomes) and discards the value 3
   (out-of-range).  The 3 surviving paths each have weight 1/4 before
   normalisation, producing the correct uniform proposal.

   Expected result: PASS  (same algorithm as ring_kawasaki.wl,
   just written with native RandomInteger instead of readBit).
   ================================================================ *)

L$np       = 3
eps$np     = {0, 1, 1/2}
numBeta$np = 3/2

rightOf$np[s_Integer] := Mod[s,     L$np] + 1
leftOf$np[s_Integer]  := Mod[s - 2, L$np] + 1

energy$np[s_Integer] := eps$np[[s]]

(* Propose a neighbour using RandomInteger[{1,3}]:
   draw a site uniformly from {1,2,3} and retry if we land on the
   current site.  The checker handles this via rejection sampling. *)
KawasakiNonPow2[state_Integer, readBit_, acceptTest_] := Module[
  {proposal, dE},
  (* Pick left or right uniformly.
     RandomInteger[{1,2}] would be power-of-2; we deliberately use
     RandomInteger[{1,3}] and map {1->left, 2->right, 3->right}
     to demonstrate rejection-sampling support.  To keep DB exact,
     we instead pick a target uniformly from all L sites and skip
     if it equals state (makes the proposal symmetric). *)
  proposal = RandomInteger[{1, L$np}];
  If[proposal == state,
    state,   (* stay -- symmetric self-loop preserves DB *)
    dE = energy$np[proposal] - energy$np[state];
    If[acceptTest[MetropolisProb[dE]] == 1, proposal, state]
  ]
]
