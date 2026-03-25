(* ================================================================
   2D Kawasaki dynamics on a periodic square lattice (torus) -- PASSES
   ================================================================

   System: L x L periodic square lattice (torus) with labeled
           (distinguishable) particles.  State[[s]] = 0 (empty) or
           k in {1,...,N} (particle of type k).  L is NOT hard-coded;
           it is inferred from the bit string via BitsToState.

   Bit encoding -- same bijective integer map as kawasaki_1d.wl:
     The bit string is a binary integer ID.  $decode[ID] returns a
     labeled-particle array of some length m.  This file accepts only
     those IDs where m is a perfect square, interpreting the state as
     an L x L grid (L = Sqrt[m]).

     Valid 2D systems within short bit strings:
       IDs 1-2:   L=1 (1x1 torus, 1 site)
       IDs 24-88: L=2 (2x2 torus, 4 sites, 65 states)
       With MaxBitString=1111111 (IDs 0-127) both L=1 and L=2 are
       fully covered.

   Site numbering (1-indexed, row-major):
     site s -> row r = Ceiling[s/L], col c = Mod[s-1,L]+1
       1  2  ...  L
      L+1 ...   2L
       ...
     (L-1)L+1 ...  L^2

   Move:   Standard local Kawasaki on the torus.
           Enumerate 2*L^2 directed bonds: first L^2 horizontal
           (each site to its right neighbour), then L^2 vertical
           (each site to its down neighbour).  Each physical bond
           appears twice, giving a symmetric proposal with equal
           weight 1/L^2 per physical bond.
           Pick b uniformly from {0,...,2*L^2-1} via
           RandomInteger[{0, 2*L*L-1}].  Swap the two sites;
           apply Metropolis acceptance.

   Energy: Same nearest-neighbour pairwise coupling as kawasaki_1d.wl
           (J<lo><hi> symbols), now summed over horizontal and vertical
           bonds of the torus.

   Display: States are shown as L rows separated by '|', e.g.
            {1,2}|{0,3} for L=2.

   Run with:
     wolframscript -file check.wls examples3/kawasaki_2d.wl \
                   MaxBitString=1111111
   This covers L=1 and L=2 systems (all particle-number sectors).
   ================================================================ *)


(* ================================================================
   Section 1 -- Bijective integer encoding (identical to 1D file)
   ================================================================ *)

$cL[L_Integer]            := $cL[L]   = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_Integer]         := $cLPre[L]  = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_Integer, N_Integer] :=
  $cLNPre[L, N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_Integer, L_Integer, N_Integer] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[
      While[Binomial[x, i] > r, x--];
      pos[[i]] = x;
      r -= Binomial[x, i];
      x--,
      {i, N, 1, -1}];
    pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[
      idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
      rank += idx * Factorial[n - i];
      elems = Delete[elems, idx + 1],
      {i, n}];
    rank]

$unrankPerm[k_Integer, n_Integer] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[
      idx = Quotient[r, Factorial[i - 1]];
      r   = Mod[r, Factorial[i - 1]];
      AppendTo[perm, elems[[idx + 1]]];
      elems = Delete[elems, idx + 1],
      {i, n, 1, -1}];
    perm]

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
   Section 2 -- Symbolic coupling constants (shared with 1D)
   ================================================================ *)

$maxCouplingN = 4

symParams = <|"couplings" ->
  Flatten @ Table[
    If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
    {a, 1, $maxCouplingN}, {b, 1, $maxCouplingN}]|>

$pairJ[a_Integer, b_Integer] :=
  If[a == 0 || b == 0, 0,
     ToExpression["J" <> ToString[Min[a, b]] <> ToString[Max[a, b]]]]


(* ================================================================
   Section 3 -- 2D lattice helpers (L inferred from state length)
   ================================================================ *)

(* Right neighbour of site s on an L x L torus (row-major, 1-indexed) *)
$right2D[s_Integer, L_Integer] :=
  With[{r = Ceiling[s / L], c = Mod[s - 1, L] + 1},
    (r - 1) * L + Mod[c, L] + 1]

(* Down neighbour of site s on an L x L torus *)
$down2D[s_Integer, L_Integer] :=
  With[{r = Ceiling[s / L], c = Mod[s - 1, L] + 1},
    Mod[r, L] * L + c]

(* The two sites connected by directed bond b on the L x L torus.
   b in {0,...,L^2-1}:   horizontal (right) bonds, site b+1 -> right
   b in {L^2,...,2L^2-1}: vertical (down) bonds, site b-L^2+1 -> down *)
$bond2D[L_Integer, b_Integer] :=
  If[b < L * L,
    With[{s = b + 1}, {s, $right2D[s, L]}],
    With[{s = b - L * L + 1}, {s, $down2D[s, L]}]]


(* ================================================================
   Section 4 -- Energy function
   ================================================================ *)

(* Sum pairJ over all L^2 horizontal and L^2 vertical bonds of the torus.
   L is inferred as Sqrt[Length[state]]. *)
energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    Total[Table[
      $pairJ[state[[s]], state[[$right2D[s, L]]]] +
      $pairJ[state[[s]], state[[$down2D[s, L]]]],
      {s, Length[state]}]]]


(* ================================================================
   Section 5 -- Algorithm: standard local Kawasaki on a 2D torus
   ================================================================ *)

(* L is computed from the state; no variable is hard-coded.
   Pick one of the 2*L^2 directed bonds uniformly; swap the two
   connected sites; apply Metropolis acceptance.
   The proposal is symmetric: each physical bond is enumerated
   twice (as both directed versions), so T_proposal(i->j) = T_proposal(j->i). *)
Algorithm[state_List] :=
  Module[{L, b, sites, s1, s2, newState, dE},
    L = Round[Sqrt[Length[state]]];
    b = RandomInteger[{0, 2 * L * L - 1}];
    sites = $bond2D[L, b];
    s1 = sites[[1]]; s2 = sites[[2]];
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ================================================================
   Section 6 -- BitsToState
   ================================================================ *)

(* Accept only IDs whose decoded array has perfect-square length,
   since only those correspond to valid L x L grids. *)
BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, m, sqrtM},
    If[id == 0, Return[None]];
    state  = $decode[id];
    m      = Length[state];
    sqrtM  = Sqrt[m];
    If[!IntegerQ[sqrtM], Return[None]];
    state]


(* ================================================================
   Section 7 -- DisplayState: show state as an L x L grid
   ================================================================ *)

(* Format the state as rows separated by '|'.
   Example: {1,2,3,0} (L=2) -> "{1,2}|{3,0}" *)
DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table[
        "{" <> StringRiffle[ToString /@ state[[(r - 1) * L + 1 ;; r * L]], ","] <> "}",
        {r, 1, L}],
      "|"]]

numBeta = 1
