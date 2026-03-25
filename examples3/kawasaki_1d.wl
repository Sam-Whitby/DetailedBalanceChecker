(* ================================================================
   1D Kawasaki dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   System: 1D periodic ring of L sites, with labeled (distinguishable)
           particles of types 1, 2, 3, ... and holes (type 0).
           L is NOT hard-coded; it is inferred from the bit string
           supplied to BitsToState.

   State:  A list of L integers.  state[[i]] = 0 means site i is empty;
           state[[i]] = k (k >= 1) means site i holds particle of type k.
           The particle types present in any valid state form a
           consecutive set {1, 2, ..., n} for some n >= 0 (no gaps,
           no duplicates).

   Bit encoding (2 bits per site):
     Each pair of consecutive bits encodes one site value:
       00 -> 0 (empty)    01 -> 1 (type 1)
       10 -> 2 (type 2)   11 -> 3 (type 3)
     A bit string of length N (padded to even length with one leading
     zero if N is odd) encodes an L = N/2 system.
     BitsToState returns None for encodings that violate the
     consecutive-labeling rule (e.g. type 2 present without type 1).

   Move:   Standard local Kawasaki dynamics.
           Pick a bond b uniformly from the L bonds of the periodic
           ring (RandomInteger[{0, L-1}]).  Bond b connects sites
           b+1 and Mod[b+1, L]+1.  Swap the contents of the two sites
           and apply Metropolis acceptance.  This covers all cases:
           particle-particle swap, particle-hole hop, hole-hole (trivial).

   Why DB holds:
     The proposal is symmetric -- bond b is equally likely to be chosen
     in the forward and reverse steps -- and Metropolis acceptance
     satisfies the detailed-balance condition exactly.

   Energy: Nearest-neighbour pairwise coupling between occupied sites.
           Three symbolic coupling constants (J12, J13, J23) describe
           the interaction between each ordered pair of particle types.
           Random numerical values are assigned automatically for the
           numerical MCMC check.

   Run with:
     wolframscript -file check.wls examples3/kawasaki_1d.wl \
                   MaxBitString=111111
   (tests all 1-6 bit strings, covering L=1, L=2, L=3 systems)
   ================================================================ *)


(* ================================================================
   Symbolic coupling constants.
   pairJ[a, b] = coupling between particle types a and b at adjacent
   sites (a, b >= 1).  Symmetric: pairJ[a,b] = pairJ[b,a].
   Types supported up to 3 (the maximum encodable with 2 bits).
   ================================================================ *)
symParams = <|"couplings" -> {J12, J13, J23}|>

pairJ[a_Integer, b_Integer] :=
  With[{lo = Min[a, b], hi = Max[a, b]},
    Which[
      lo == 1 && hi == 2, J12,
      lo == 1 && hi == 3, J13,
      lo == 2 && hi == 3, J23,
      True, 0]]

(* ================================================================
   Energy
   Sum of nearest-neighbour coupling constants over all L bonds of
   the periodic ring.  L is inferred from the state length.
   Empty sites (type 0) contribute zero to the energy.
   ================================================================ *)
energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[
      With[{a = state[[i]], b = state[[Mod[i, L] + 1]]},
        If[a == 0 || b == 0, 0, pairJ[a, b]]],
      {i, L}]]]

(* ================================================================
   Algorithm: standard local Kawasaki on a periodic 1D ring.
   L is computed from the state; no variable is hard-coded.
   ================================================================ *)
Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L = Length[state];
    (* Pick one of the L bonds uniformly *)
    b  = RandomInteger[{0, L - 1}];
    s1 = b + 1;
    s2 = Mod[b + 1, L] + 1;
    (* Swap the two sites (works for any pair: empty/empty, empty/particle,
       particle/particle) *)
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]
  ]

(* ================================================================
   BitsToState: map a bit string to a labeled-particle ring state.

   Encoding: 2 bits per site.
     bit pair -> site value:  00->0, 01->1, 10->2, 11->3
   System size: L = (padded length) / 2
     Odd-length strings are padded with one leading 0 to make L exact.

   Validity rules (returns None if violated):
     (1) The non-zero values in the state must be distinct.
     (2) They must form the set {1, 2, ..., n} for n = count of particles
         (no gaps -- e.g. {1,3} is invalid because type 2 is missing).

   Examples:
     {0}           -> pad to {0,0}  -> L=1, state {0}        (empty ring)
     {1}           -> pad to {0,1}  -> L=1, state {1}        (one particle)
     {0,1,1,0}     -> L=2, state {1,2}   (two-particle sector)
     {1,0,0,1}     -> L=2, state {2,1}   (two-particle sector)
     {0,1,1,0,1,1} -> L=3, state {1,2,3} (fully occupied)
   ================================================================ *)
BitsToState[bits_List] :=
  Module[{paddedBits, L, pairs, state, nonzero, sortedVals},
    (* Pad odd-length strings with one leading zero *)
    paddedBits = If[OddQ[Length[bits]], Prepend[bits, 0], bits];
    L = Length[paddedBits] / 2;
    If[L == 0, Return[None]];
    (* Decode 2-bit pairs into site values *)
    pairs = Partition[paddedBits, 2];
    state = FromDigits[#, 2] & /@ pairs;
    (* Validate: non-zero values must be distinct and form {1,...,n} *)
    nonzero    = Select[state, # != 0 &];
    sortedVals = Sort[nonzero];
    If[Length[nonzero] == 0, Return[state]];          (* all-empty is valid *)
    If[sortedVals =!= Range[Max[sortedVals]],         (* gaps, e.g. {1,3} *)
       Return[None]];
    If[Length[Union[nonzero]] != Length[nonzero],     (* duplicate types   *)
       Return[None]];
    state
  ]

numBeta = 1
