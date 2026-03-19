(* ================================================================
   Example: energy functions that make symbolic checking impossible
   ================================================================
   Two energy functions that CANNOT be used for the symbolic DB check:

   (A) AbsoluteTime in energy: returns a different value each call,
       making path weights time-dependent and the symbolic check
       meaningless.  CheckEnergySafety catches this before BFS starts.

   (B) RandomVariate in energy: stochastic energy makes the DB
       residuals random numbers rather than symbolic expressions.
       CheckEnergySafety catches this too.

   In both cases RunFullCheck returns:
     <|"pass" -> False, "error" -> "energy function contains ...", ...|>
   ================================================================ *)

L$ue       = 3
numBeta$ue = 1

(* Algorithm is a normal Kawasaki ring -- the bug is in the energy. *)
KawasakiUE[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = Mod[state + If[b == 1, 1, -1] - 1, L$ue] + 1;
  dE  = energyTime$ue[nbr] - energyTime$ue[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]

(* ---- (A) AbsoluteTime in energy: non-deterministic ---- *)
energyTime$ue[s_Integer] :=
  Mod[Round[AbsoluteTime[] * 1000] + s, 5]

(* ---- (B) RandomVariate in energy: stochastic ---- *)
energyRand$ue[s_Integer] :=
  s + RandomVariate[NormalDistribution[0, 0.1]]
