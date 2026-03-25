(* ================================================================
   1D Kawasaki dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   System: 1D periodic ring of L sites with labeled (distinguishable)
           particles.  State[[i]] = 0 (empty) or k in {1,...,N}
           (particle of type k).  All N particle types are distinct
           and form a consecutive set {1,...,N}.  L is NOT hard-coded;
           it is inferred from the bit string via BitsToState.

   Bit encoding -- bijective integer map:
     The bit string is interpreted as a binary integer ID.  The ID
     bijectively encodes a (L, N, positions, permutation) tuple:

       ID = countLPrefix[L]           (skip all smaller-L arrays)
          + countLNPrefix[L, N]       (skip smaller-N arrays at this L)
          + rankCombo(positions) * N! (within (L,N): position block)
          + rankPerm(labels)          (within position block: label order)

     where:
       countL[L]       = Sum[C(L,k)*k!, {k,0,L}]  (all labeled arrays of length L)
       countLPrefix[L] = Sum[countL[l], {l,0,L-1}] (cumulative offset)

     This gives a clean bijection.  ID=0 (all-zero bit string) decodes
     to the empty array {} and is mapped to None.  Every non-zero ID
     uniquely identifies one labeled-particle array of some length L,
     which becomes the seed state.  Short bit strings (leading zeros)
     decode to small L systems; longer strings reach larger L.

     First few IDs:
       1 -> {0}      (L=1, empty)
       2 -> {1}      (L=1, one particle)
       3 -> {0,0}    (L=2, empty)
       4 -> {1,0}    (L=2, particle at site 1)
       5 -> {0,1}    (L=2, particle at site 2)
       6 -> {1,2}    (L=2, two particles)
       7 -> {2,1}    (L=2, two particles reversed)
       8 -> {0,0,0}  (L=3 ...)
       ...
       65 -> {1,2,3,4} (L=4, fully occupied)

     With MaxBitString=1111111 (IDs 0-127) all L=1..4 systems
     are covered; with MaxBitString=111111111 (IDs 0-511) L=5 too.

   Move:   Standard local Kawasaki on the periodic ring.
           Pick bond b in {0,...,L-1} uniformly (RandomInteger[{0,L-1}]).
           Swap the two adjacent sites; apply Metropolis acceptance.

   Energy: Nearest-neighbour pairwise coupling.  Coupling between
           particle types a and b is the symbol J<min(a,b)><max(a,b)>
           (e.g. J12, J23, J14 ...).  All coupling symbols listed in
           symParams are kept symbolic during the symbolic check and
           assigned random numerical values for the MCMC check.

   Run with:
     wolframscript -file check.wls examples3/kawasaki_1d.wl \
                   MaxBitString=1111111
   This covers L=1,2,3,4 (all particle-number sectors of each).
   ================================================================ *)


(* ================================================================
   Section 1 -- Bijective integer encoding
   ================================================================ *)

(* Counting: number of valid labeled arrays of given (L, N) *)
$cL[L_Integer]            := $cL[L]   = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_Integer]         := $cLPre[L]  = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_Integer, N_Integer] :=
  $cLNPre[L, N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

(* Rank a sorted 0-indexed combination {p_1,...,p_N} (combinatorial
   number system: rank = Sum[C(p_i, i), {i,1,N}]) *)
$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

(* Inverse: k-th sorted 0-indexed N-combination from {0,...,L-1} *)
$unrankCombo[rank_Integer, L_Integer, N_Integer] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[
      While[Binomial[x, i] > r, x--];
      pos[[i]] = x;
      r -= Binomial[x, i];
      x--,
      {i, N, 1, -1}];
    pos]

(* Rank a permutation of {1,...,N} in Lehmer (factorial number) order *)
$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[
      idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
      rank += idx * Factorial[n - i];
      elems = Delete[elems, idx + 1],
      {i, n}];
    rank]

(* Inverse: k-th permutation of {1,...,N} (0-indexed Lehmer order) *)
$unrankPerm[k_Integer, n_Integer] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[
      idx = Quotient[r, Factorial[i - 1]];
      r   = Mod[r, Factorial[i - 1]];
      AppendTo[perm, elems[[idx + 1]]];
      elems = Delete[elems, idx + 1],
      {i, n, 1, -1}];
    perm]

(* Decode an integer ID to the corresponding labeled-particle array.
   ID=0 -> {} (empty array; BitsToState returns None for this). *)
$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos  = Quotient[r, Factorial[N]];
    rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N];
    perm = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}];
    arr]


(* ================================================================
   Section 2 -- Symbolic coupling constants
   ================================================================ *)

(* Coupling symbol for particle types a and b (a != b, a,b >= 1).
   Uses the naming convention J<lo><hi> where lo=Min[a,b], hi=Max[a,b].
   Additional types can be supported by increasing $maxCouplingN. *)
$maxCouplingN = 4  (* supports particle types 1..4, i.e. L up to 4 fully occupied *)

symParams = <|"couplings" ->
  Flatten @ Table[
    If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
    {a, 1, $maxCouplingN}, {b, 1, $maxCouplingN}]|>

(* pairJ[a, b]: coupling between types a and b at adjacent sites.
   Returns 0 if either site is empty; otherwise the symbolic J symbol. *)
$pairJ[a_Integer, b_Integer] :=
  If[a == 0 || b == 0, 0,
     ToExpression["J" <> ToString[Min[a, b]] <> ToString[Max[a, b]]]]


(* ================================================================
   Section 3 -- Energy function
   ================================================================ *)

(* Sum pairJ over all L nearest-neighbour bonds of the periodic ring.
   L is inferred from the state length. *)
energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[
      $pairJ[state[[i]], state[[Mod[i, L] + 1]]],
      {i, L}]]]


(* ================================================================
   Section 4 -- Algorithm: standard local Kawasaki on a 1D ring
   ================================================================ *)

(* L is computed from the state; no variable is hard-coded.
   Pick bond b uniformly in {0,...,L-1}: sites b+1 and Mod[b+1,L]+1.
   Swap the two sites unconditionally, then apply Metropolis. *)
Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L  = Length[state];
    b  = RandomInteger[{0, L - 1}];
    s1 = b + 1;
    s2 = Mod[b + 1, L] + 1;
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ================================================================
   Section 5 -- BitsToState
   ================================================================ *)

(* Convert a bit string to a seed state for the checker.
   The bit string is read as a binary integer ID; ID=0 (all-zero
   string of any length) returns None.  For any other ID, $decode
   returns the unique labeled-particle array with that index.
   Different lengths of all-zero strings all decode to ID=0, so
   they are all skipped -- the checker's deduplication handles
   remaining duplicates (e.g. "1" and "01" both give ID=1 = {0}). *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
