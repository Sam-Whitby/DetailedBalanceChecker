(* ================================================================
   1D Kawasaki with rightward bond bias -- FAILS detailed balance
   ================================================================

   Identical to kawasaki_1d.wl except for one line in Algorithm:

     CORRECT:  b = RandomInteger[{0, L - 1}]
     BUGGY:    b = Mod[RandomInteger[{0, L}], L]

   The buggy version draws from {0,...,L} (L+1 values) and wraps the
   extra value L back to bond 0, so bond 0 is proposed with probability
   2/(L+1) while bonds 1,...,L-1 each have probability 1/(L+1).

   This asymmetry breaks detailed balance: transitions involving bond 0
   have unequal forward and reverse proposal rates. The symbolic check
   detects this as a persistent current between sites 1 and 2.
   ================================================================ *)


(* ---- Encoding (identical to kawasaki_1d.wl) ------------------------------ *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]

$pairJ[a_, b_] :=
  If[a == 0 || b == 0, 0,
     ToExpression["J" <> ToString[Min[a,b]] <> ToString[Max[a,b]]]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
        {a, types}, {b, types}]|>]

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[$pairJ[state[[i]], state[[Mod[i, L] + 1]]], {i, L}]]]


(* ---- Buggy algorithm: bond 0 proposed twice as often -------------------- *)

Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L  = Length[state];
    b  = Mod[RandomInteger[{0, L}], L];      (* BUG: bond 0 has double weight *)
    s1 = b + 1; s2 = Mod[b + 1, L] + 1;
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
