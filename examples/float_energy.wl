(* ================================================================
   Example: algorithm written with floating-point literals
   ================================================================
   L=3 ring Kawasaki Metropolis using native RandomReal[]/
   RandomInteger[1] with float literals for probabilities.

   The checker rationalises Real values inside acceptTest at the
   point of consumption: 0.5 → 1/2, Exp[-β*dE] stays symbolic.
   This means algorithms that use floats like `RandomReal[] < 0.5`
   or `acceptTest[0.5]` are handled exactly without the user needing
   to write `1/2` explicitly.

   Expected result: PASS
   ================================================================ *)

L$fe       = 3
eps$fe     = {0, 1, 1/2}    (* energies remain exact rationals *)
numBeta$fe = 3/2

rightOf$fe[s_Integer] := Mod[s,     L$fe] + 1
leftOf$fe[s_Integer]  := Mod[s - 2, L$fe] + 1

energy$fe[s_Integer] := eps$fe[[s]]

(* Direction chosen with float-comparison RandomReal[] < 0.5
   (not readBit[], to demonstrate float literal handling).
   Metropolis acceptance is also called with the float-literal 1.0
   for the downhill branch to further exercise rationalisation. *)
KawasakiFloat[state_Integer, readBit_, acceptTest_] := Module[
  {nbr, dE, p},
  (* RandomReal[] < 0.5 is intercepted: $dbc$rand UpValue fires,
     then acceptTest[0.5] is called; 0.5 is rationalised to 1/2. *)
  nbr = If[RandomReal[] < 0.5, leftOf$fe[state], rightOf$fe[state]];
  dE  = energy$fe[nbr] - energy$fe[state];
  (* Pass a float acceptance probability. 1.0 → 1, Exp[...] stays symbolic. *)
  p   = If[dE <= 0, 1.0, Exp[(-\[Beta]) * dE]];
  If[acceptTest[p] == 1, nbr, state]
]
