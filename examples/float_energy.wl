(* ================================================================
   Example: float energy values and float acceptance probabilities
   ================================================================
   Same L=3 ring Kawasaki physics as ring_kawasaki.wl, but written
   entirely with floating-point literals:

     eps$fe = {0., 1., 0.5}          -- float energy array
     RandomReal[] < 0.5              -- float comparison for direction
     acceptTest[1.0] / acceptTest[p] -- float acceptance probabilities

   The checker handles these as follows:

   1. Energy floats: rationalised in CheckDetailedBalance before
      FullSimplify: 0.5 -> 1/2, 1.0 -> 1, 0. -> 0.
   2. MetropolisProb[dE] with concrete float dE evaluates immediately
      (Piecewise conditions become True/False), giving Exp[-β*dE] or 1.
   3. Reals inside acceptTest are rationalised via
        p /. {r_Real :> Rationalize[r]}
      so Exp[-1.0β] -> Exp[-β], 1.0 -> 1, etc.

   Expected result: PASS
   ================================================================ *)

L$fe       = 3
eps$fe     = {0., 1., 0.5}     (* float energy values *)
numBeta$fe = 1.5               (* float beta also fine for numerical run *)

rightOf$fe[s_Integer] := Mod[s,     L$fe] + 1
leftOf$fe[s_Integer]  := Mod[s - 2, L$fe] + 1

energy$fe[s_Integer] := eps$fe[[s]]

(* Algorithm using float literals throughout.
   RandomReal[] < 0.5 : direction with float comparison
   1.0 / Exp[-β*dE]   : acceptance with float literal *)
KawasakiFloat[state_Integer, readBit_, acceptTest_] := Module[
  {nbr, dE, p},
  nbr = If[RandomReal[] < 0.5, leftOf$fe[state], rightOf$fe[state]];
  dE  = energy$fe[nbr] - energy$fe[state];
  p   = If[dE <= 0, 1.0, Exp[(-\[Beta]) * dE]];
  If[acceptTest[p] == 1, nbr, state]
]
