(* ================================================================
   1D Kawasaki swap dynamics (labeled particles) -- PASSES DB
   ================================================================

   System: L=4 ring, 3 distinguishable (labeled) particles, 1 hole.
   State:  length-4 list {s1,s2,s3,s4} where s_i in {1,2,3} is the
           label of the particle at site i, or 0 if the site is empty.
           With 3 labeled particles the state space has 4*3*2 = 24 states.

   Energy: pairwise coupling between adjacent occupied pairs.
           Only the IDENTITIES of the two particles matter, not their
           order, so E uses Min/Max to pick the right coupling constant.
           Three coupling constants: J12 (between labels 1 and 2),
                                     J13 (between labels 1 and 3),
                                     J23 (between labels 2 and 3).
           All three are symbolic; random numerical values are assigned
           automatically during the numerical MCMC check.

   Move:   Kawasaki SWAP.
           Pick a random bond (RandomInteger[{0,3}], 2 bits, exact).
           Bond b connects sites (b+1) and Mod[b+1,L]+1.
           - Both occupied:  swap the two particles, Metropolis acceptance.
           - Not both occupied: stay put.

   Proposal is symmetric (bond b chosen with probability 1/4 in both
   directions) and Metropolis is correct, so DB holds exactly.
   ================================================================ *)

(* ---- System parameters ---- *)
L$kswap = 4

(* Symbolic pairwise coupling constants.
   Random values assigned automatically when symParams is provided. *)
symParams = <|"couplings" -> {J12, J13, J23}|>

numBeta = 1

(* ---- Energy ---- *)
(* Coupling between two occupied adjacent sites with labels a and b *)
pairJ[a_Integer, b_Integer] :=
  With[{lo = Min[a, b], hi = Max[a, b]},
    If[lo == 1 && hi == 2, J12,
     If[lo == 1 && hi == 3, J13, J23]]]

(* Total energy: sum over all adjacent bonds of the coupling for
   the pair of labels present (0 if either site is empty) *)
energy[state_List] :=
  Total[Table[
    With[{a = state[[i]], b = state[[Mod[i, L$kswap] + 1]]},
      If[a == 0 || b == 0, 0, pairJ[a, b]]],
    {i, L$kswap}]]

(* ================================================================
   Algorithm: Kawasaki swap.
   Pick bond, swap the two particles if both sites are occupied,
   otherwise stay.  Metropolis acceptance applied to the energy
   difference.
   ================================================================ *)
Algorithm[state_List] :=
  Module[{b, s1, s2, newState, dE},
    b  = RandomInteger[{0, L$kswap - 1}];
    s1 = b + 1;
    s2 = Mod[b + 1, L$kswap] + 1;
    If[state[[s1]] == 0 || state[[s2]] == 0,
      state,   (* one or both sites empty: stay *)
      (* Swap the two particles *)
      newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
      dE = energy[newState] - energy[state];
      If[RandomReal[] < MetropolisProb[dE], newState, state]
    ]
  ]

(* ---- Checker interface ---- *)
seedState = {1, 2, 3, 0}   (* particles 1,2,3 at sites 1,2,3; empty at 4 *)
