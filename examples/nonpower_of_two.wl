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
(* Symbolic site energies -- kept unassigned during the symbolic check.
   Random numerical values are assigned automatically by RunFullCheck
   when "SysParams" -> params$np is supplied. *)
eps$np     = {\[Epsilon]np1, \[Epsilon]np2, \[Epsilon]np3}
params$np  = <|"L" -> L$np, "eps" -> eps$np|>
numBeta$np = 3/2

rightOf$np[s_Integer] := Mod[s,     L$np] + 1
leftOf$np[s_Integer]  := Mod[s - 2, L$np] + 1

energy$np[s_Integer] := eps$np[[s]]

(* Propose a neighbour using RandomInteger[{1,3}]:
   draw a site uniformly from {1,2,3} and retry if we land on the
   current site.  The checker handles this via rejection sampling. *)
KawasakiNonPow2[state_Integer] := Module[
  {proposal, dE},
  (* Pick a target site uniformly from all L sites.
     RandomInteger[{1,3}] has 3 outcomes (not a power of 2):
     the checker reads 2 bits and discards the value 3 via
     rejection sampling so the proposal stays symmetric. *)
  proposal = RandomInteger[{1, L$np}];
  If[proposal == state,
    state,   (* stay -- symmetric self-loop preserves DB *)
    dE = energy$np[proposal] - energy$np[state];
    If[RandomReal[] < MetropolisProb[dE], proposal, state]
  ]
]
