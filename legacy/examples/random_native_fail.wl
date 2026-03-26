(* ================================================================
   Example: Inverted acceptance written with native random calls -- FAILS
   ================================================================

   Same L=3 ring system but the acceptance condition is INVERTED:
   the algorithm accepts an UPHILL move with probability 1 and
   a downhill move with probability Exp[+beta*dE] (wrong sign).
   Written using native RandomReal[] instead of acceptTest[].

   The checker intercepts RandomReal[] correctly and detects the
   detailed balance violation exactly as it would for an algorithm
   written with acceptTest[].

   This demonstrates that automatic interception does not mask bugs:
   the DB check still correctly identifies the failure.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "random_native_pass.wl"]

(* Wrong acceptance: prefer uphill moves *)
KawasakiNativeWrongSign[state_Integer] := Module[
  {dir, nbr, dE, pWrong},
  dir    = RandomInteger[];
  nbr    = Mod[state + If[dir == 1, 1, -1] - 1, L$rn] + 1;
  dE     = energy$rn[nbr] - energy$rn[state];
  (* Wrong sign: accept downhill with prob Exp[+beta*dE] < 1 *)
  pWrong = Piecewise[{{1, dE <= 0}, {Exp[\[Beta] * dE], dE > 0}}];
  If[RandomReal[] < pWrong, nbr, state]
]
