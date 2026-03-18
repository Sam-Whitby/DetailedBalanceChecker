(* ================================================================
   Example 1: Kawasaki dynamics on a ring -- PASSES detailed balance
   ================================================================

   System: L sites on a ring, 1 particle.
   State:  integer 1..L giving the particle's site.
   Energy: E(i) = eps_i  (site-dependent, symbolic in eps1,eps2,...,epsL).
   Move:   read 1 bit -> choose left (0) or right (1) neighbour,
           then apply Metropolis acceptance.

   Because the proposal is symmetric (each neighbour chosen with p=1/2)
   and Metropolis is used, detailed balance holds exactly.
   ================================================================ *)

(* ---- Parameters ---- *)
L$rk = 3;   (* ring size *)

(* Symbolic site energies eps1, eps2, eps3 *)
eps$rk = Table[Symbol["eps" <> ToString[i]], {i, L$rk}];

(* Symbolic beta -- must be the symbol \[Beta] used by MetropolisProb *)
(* (already defined as the global symbol in dbc_core.wl) *)

(* All valid states *)
allStates$rk = Range[L$rk];

(* Symbolic energy *)
symEnergy$rk[s_Integer] := eps$rk[[s]]

(* Numerical energy -- use specific values for validation *)
numEps$rk = {0.0, 1.0, 0.5};   (* eps1=0, eps2=1, eps3=0.5 *)
numBeta$rk = 1.5;
numEnergy$rk[s_Integer] := numBeta$rk * numEps$rk[[s]]

(* ---- Algorithm: Kawasaki with Metropolis ---- *)
KawasakiGood[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];                              (* 0 = left, 1 = right *)
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = symEnergy$rk[nbr] - symEnergy$rk[state];
  p   = MetropolisProb[dE];
  {{p, nbr}, {1 - p, state}}                   (* accept / reject *)
]

(* Numerical variant (replaces symbolic energy with numeric) *)
KawasakiGoodNum[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = numBeta$rk * (numEps$rk[[nbr]] - numEps$rk[[state]]);
  p   = N[If[dE <= 0, 1, Exp[-dE]]];
  {{p, nbr}, {1 - p, state}}
]
