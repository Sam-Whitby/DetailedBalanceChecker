(* ================================================================
   2D Kawasaki dynamics on a periodic square lattice (torus) -- PASSES
   ================================================================

   System: L x L periodic square lattice (torus) with labeled
           (distinguishable) particles of types 1, 2, 3, ... and
           holes (type 0).
           L is NOT hard-coded; it is inferred from the bit string.

   State:  A list of L^2 integers in row-major order.
           state[[s]] = 0 means site s is empty;
           state[[s]] = k means site s holds particle of type k.
           Site numbering (1-indexed, row-major):
             site s -> row r = ceil(s/L), col c = mod(s-1, L)+1
           The particle types present form {1, 2, ..., n} (no gaps,
           no duplicates), as in the 1D file.

   Bit encoding (2 bits per site):
     Each pair of bits encodes one site: 00->0, 01->1, 10->2, 11->3.
     A bit string of length N encodes a system of
       L = ceil(sqrt(ceil(N/2)))
     sites per side.  Short bit strings are padded with leading zeros
     to reach 2*L^2 bits, so every bit string maps to some system.
     BitsToState returns None for encodings that violate the
     consecutive-labeling rule.

     Valid lengths producing a new L:
       N = 2  ->  L=1 (1x1 torus, 1 site)
       N = 8  ->  L=2 (2x2 torus, 4 sites)
       N = 18 ->  L=3 (3x3 torus, 9 sites)
     Shorter strings (e.g. N=3..7) are padded to the next valid length.

   Move:   Standard local Kawasaki dynamics on the torus.
           There are 2*L^2 directed bonds: L^2 horizontal (right) and
           L^2 vertical (down).  Each physical bond appears twice
           (once in each direction), giving a symmetric proposal.
           Pick bond b uniformly from {0, ..., 2*L^2 - 1}
           (RandomInteger[{0, 2*L^2-1}]).
             b < L^2:  horizontal bond -- sites (r,c) and (r, c+1 mod L)
             b >= L^2: vertical bond   -- sites (r,c) and (r+1 mod L, c)
           Swap the two sites and apply Metropolis acceptance.

   Why DB holds:
     The proposal is symmetric: bond b giving sites {s1,s2} and bond
     b' giving sites {s2,s1} both appear in the enumeration with equal
     probability 1/(2*L^2).  Metropolis is correct.

   Energy: Same nearest-neighbour pairwise coupling as kawasaki_1d.wl.
           Coupling constants J12, J13, J23 between particle type pairs.
           Now counted over all horizontal and vertical adjacent bonds.

   Run with:
     wolframscript -file check.wls examples3/kawasaki_2d.wl \
                   MaxBitString=11111111
   (tests all 1-8 bit strings, covering L=1 and L=2 systems;
    for L=3 use MaxBitString=111111111111111111  [18 ones], but note
    that with 84 or more states the symbolic check may take several
    minutes)
   ================================================================ *)


(* ================================================================
   Symbolic coupling constants (shared with 1D convention).
   pairJ[a, b] = coupling between particle types a and b.
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
   2D lattice helpers.
   L is passed as an argument; no global variable is used.
   Site s is at row r = Ceiling[s/L], col c = Mod[s-1, L]+1 (1-indexed).
   ================================================================ *)

(* Right neighbour of site s on an L x L torus *)
rightOf2D[s_Integer, L_Integer] :=
  With[{r = Ceiling[s / L], c = Mod[s - 1, L] + 1},
    (r - 1) * L + Mod[c, L] + 1]

(* Down neighbour of site s on an L x L torus *)
downOf2D[s_Integer, L_Integer] :=
  With[{r = Ceiling[s / L], c = Mod[s - 1, L] + 1},
    Mod[r, L] * L + c]

(* ================================================================
   Energy
   Sum of pairJ over all nearest-neighbour bonds on the L x L torus.
   Each horizontal and each vertical bond is counted once.
   L is inferred from the state length.
   ================================================================ *)
energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[Table[
      With[{a = state[[s]],
            bR = state[[rightOf2D[s, L]]],
            bD = state[[downOf2D[s,  L]]]},
        (* right bond *)
        If[a == 0 || bR == 0, 0, pairJ[a, bR]] +
        (* down bond  *)
        If[a == 0 || bD == 0, 0, pairJ[a, bD]]],
      {s, Length[state]}]]]

(* ================================================================
   getBond2D[L, b]
   Returns {s1, s2}: the two sites connected by bond index b on the
   L x L torus.
     b in {0, ..., L^2-1}          -> horizontal bond
     b in {L^2, ..., 2*L^2-1}      -> vertical bond
   ================================================================ *)
getBond2D[L_Integer, b_Integer] :=
  If[b < L * L,
    (* Horizontal bond b: site at position b+1 in row-major, right neighbour *)
    With[{s = b + 1}, {s, rightOf2D[s, L]}],
    (* Vertical bond b-L^2: site at position b-L^2+1, down neighbour *)
    With[{s = b - L * L + 1}, {s, downOf2D[s, L]}]]

(* ================================================================
   Algorithm: standard local Kawasaki on a periodic 2D torus.
   L is computed from the state; no variable is hard-coded.
   ================================================================ *)
Algorithm[state_List] :=
  Module[{L, nBonds, b, sites, s1, s2, newState, dE},
    L      = Round[Sqrt[Length[state]]];
    nBonds = 2 * L * L;
    (* Pick one bond uniformly from all 2L^2 directed bonds *)
    b      = RandomInteger[{0, nBonds - 1}];
    sites  = getBond2D[L, b];
    s1     = sites[[1]];
    s2     = sites[[2]];
    (* Swap the two sites and apply Metropolis *)
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]
  ]

(* ================================================================
   BitsToState: map a bit string to a labeled-particle 2D state.

   Encoding: 2 bits per site (same as 1D).
   Padding:  L = ceil(sqrt(ceil(N/2))).  The bit string is left-padded
             with zeros to reach 2*L^2 bits, so every input length is
             accepted.  (Short strings describe sparse/small systems.)

   The returned state is a list of L^2 integers.
   Returns None if the decoded site values violate the
   consecutive-labeling rule (same as 1D).

   Examples (L=2, 4 sites, 8-bit encoding):
     {0,0,0,0,0,0,0,0} -> {0,0,0,0}   (empty 2x2)
     {0,0,0,0,0,1,0,0} -> {0,0,1,0}   (one particle at site 3)
     {0,1,0,0,1,0,0,0} -> {1,0,2,0}   (two particles)
     For bit strings shorter than 8: padded with leading zeros.
   ================================================================ *)
BitsToState[bits_List] :=
  Module[{N = Length[bits], nSites, L, paddedBits, pairs, state,
          nonzero, sortedVals},
    (* Find L: smallest integer >= 1 such that L^2 >= ceil(N/2) *)
    nSites     = Ceiling[N / 2];
    L          = Ceiling[Sqrt[Max[nSites, 1]]];
    (* Pad to 2*L^2 bits with leading zeros *)
    paddedBits = PadLeft[bits, 2 * L^2, 0];
    (* Decode 2-bit pairs into site values *)
    pairs = Partition[paddedBits, 2];
    state = FromDigits[#, 2] & /@ pairs;
    (* Validate: non-zero values must be distinct and form {1,...,n} *)
    nonzero    = Select[state, # != 0 &];
    sortedVals = Sort[nonzero];
    If[Length[nonzero] == 0, Return[state]];
    If[sortedVals =!= Range[Max[sortedVals]], Return[None]];
    If[Length[Union[nonzero]] != Length[nonzero], Return[None]];
    state
  ]

numBeta = 1
