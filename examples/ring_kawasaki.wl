(* ================================================================
   Example 1: Kawasaki dynamics on a ring -- PASSES detailed balance
   ================================================================

   System: L sites on a ring, 1 particle.
   State:  integer 1..L (which site the particle occupies).
   Energy: E(i) = eps_i  (symbolic coupling constants eps1, eps2, ...).
   Move:   read 1 bit -> choose left (0) or right (1) neighbour,
           then Metropolis acceptance.

   The proposal is symmetric and Metropolis is correct, so
   detailed balance holds exactly.
   ================================================================ *)

(* ---- System parameters ---- *)
L$rk        = 3
allStates$rk = Range[L$rk]

(* Symbolic site energies: eps1, eps2, eps3 *)
eps$sym$rk = Table[Symbol["eps" <> ToString[i]], {i, L$rk}]

(* Symbolic bare energy (no beta -- beta is inserted by the checker) *)
symEnergy$rk[s_Integer] := eps$sym$rk[[s]]

(* Numeric parameters for validation *)
numBeta$rk  = 1.5
numEps$rk   = {0.0, 1.0, 0.5}

(* numEnergy includes beta so Exp[-numEnergy[s]] gives Boltzmann weight *)
numEnergy$rk[s_Integer] := numBeta$rk * numEps$rk[[s]]

(* ================================================================
   Symbolic algorithm: uses MetropolisProb (keeps beta symbolic).
   Call this with the symbolic checker.
   ================================================================ *)
KawasakiSym$rk[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = symEnergy$rk[nbr] - symEnergy$rk[state];   (* eps_nbr - eps_state *)
  p   = MetropolisProb[dE];                          (* Piecewise, symbolic beta *)
  {{p, nbr}, {1 - p, state}}
]

(* ================================================================
   Numerical algorithm: fully numeric, beta already in dE.
   Call this for the numerical MCMC run.
   ================================================================ *)
KawasakiNum$rk[state_Integer, readBit_] := Module[
  {dir, nbr, dE, p},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$rk] + 1;
  dE  = numEnergy$rk[nbr] - numEnergy$rk[state];   (* already beta-scaled *)
  p   = N@If[dE <= 0, 1, Exp[-dE]];
  {{p, nbr}, {1 - p, state}}
]
