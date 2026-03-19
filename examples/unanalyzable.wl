(* ================================================================
   Example: Algorithms that CANNOT be analysed by the checker
   ================================================================

   Three algorithms illustrating different failure modes:

   (A) RandomVariate -- sampling from a continuous distribution.
       Cannot be converted to readBit[]/acceptTest[] because the
       output is a continuous floating-point value, not a discrete
       bit.  The checker throws $dbc$cantHandle immediately.

   (B) Non-power-of-2 RandomInteger range.
       RandomInteger[{1,3}] produces 3 equally-likely values, but
       3 is not a power of 2 so an exact binary representation does
       not exist.  The checker reports the issue and aborts.

   (C) AbsoluteTime -- a time-dependent call.
       Even if the result happens to be used for an acceptance
       decision, the value depends on wall-clock time and cannot
       be reproduced from a bit sequence.  CheckAlgorithmSafety
       flags this before the BFS even starts.

   In each case RunFullCheck returns:
     <|"pass" -> False, "error" -> "reason...", ...|>
   ================================================================ *)

L$ua       = 3
eps$ua     = {0, 1, 1/2}
numBeta$ua = 3/2

energy$ua[s_Integer] := eps$ua[[s]]

(* ---- (A) RandomVariate: unsupported continuous distribution ---- *)
UnanalyzableVariate[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE, u},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$ua] + 1;
  dE  = energy$ua[nbr] - energy$ua[state];
  (* RandomVariate cannot be intercepted: analysis will fail *)
  u   = RandomVariate[UniformDistribution[]];
  If[u < MetropolisProb[dE], nbr, state]
]

(* ---- (B) RandomInteger with non-power-of-2 range ---- *)
UnanalyzableIntRange[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE},
  (* Range {1,3} has 3 values -- not a power of 2 *)
  dir = RandomInteger[{1, 3}];
  nbr = Mod[state + dir - 2, L$ua] + 1;
  dE  = energy$ua[nbr] - energy$ua[state];
  If[acceptTest[MetropolisProb[dE]] == 1, nbr, state]
]

(* ---- (C) AbsoluteTime: time-dependent, flagged by safety check ---- *)
UnanalyzableTime[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$ua] + 1;
  dE  = energy$ua[nbr] - energy$ua[state];
  (* Using wall-clock time as a random source -- completely wrong *)
  If[Mod[Round[AbsoluteTime[] * 1000], 2] == 0,
    nbr, state]
]
