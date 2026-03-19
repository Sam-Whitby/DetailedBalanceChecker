(* ================================================================
   Examples: Kawasaki dynamics on a ring -- Generalised API
   ================================================================

   System: L=3 ring, 1 particle.  States {1,2,3}.

   Single energy function (bare, no beta):
     energy$kn[s] -> a number (e.g. 0.0, 1.0, 0.5)
   Single algorithm for both symbolic DB check and MCMC.
   numBeta is supplied separately to RunFullCheck.

   KawasakiNew    -- correct Metropolis acceptance    [PASS]
   AlwaysAcceptNew -- always moves, no Metropolis     [FAIL]
   ================================================================ *)

L$kn       = 3
eps$kn     = {0, 1, 1/2}   (* bare site energies -- exact rationals, no beta *)
numBeta$kn = 3/2

rightOf$kn[s_Integer] := Mod[s,     L$kn] + 1   (* clockwise  *)
leftOf$kn[s_Integer]  := Mod[s - 2, L$kn] + 1   (* anti-clock *)

(* Bare energy: returns a plain number, no beta.
   \[Beta] is kept symbolic during the DB check and assigned
   numBeta$kn inside Block[...] during the MCMC run. *)
energy$kn[s_Integer] := eps$kn[[s]]

(* ================================================================
   PASSING algorithm: correct Metropolis
   ================================================================
   Reads ONE bit (direction), then calls acceptTest[MetropolisProb[dE]].
   MetropolisProb keeps \[Beta] symbolic during the symbolic check;
   Block[\[Beta]=numBeta] makes it numeric during MCMC.
   ================================================================ *)
KawasakiNew[state_Integer] := Module[
  {b, nbr, dE},
  b   = RandomInteger[];
  nbr = If[b == 0, leftOf$kn[state], rightOf$kn[state]];
  dE  = energy$kn[nbr] - energy$kn[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
]

(* ================================================================
   FAILING algorithm: always accepts, no Metropolis
   ================================================================
   Draws RandomInteger[] for direction, then always moves there.
   Drives a symmetric random walk regardless of energies ->
   cannot satisfy detailed balance when energies differ.
   ================================================================ *)
AlwaysAcceptNew[state_Integer] := Module[
  {b, nbr},
  b   = RandomInteger[];
  nbr = If[b == 0, leftOf$kn[state], rightOf$kn[state]];
  nbr   (* always move, no energy test *)
]
