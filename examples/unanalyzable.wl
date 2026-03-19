(* ================================================================
   Example: Algorithms that CANNOT be analysed by the checker
   ================================================================

   Two algorithms illustrating calls that cannot be intercepted:

   (A) RandomVariate -- sampling from a continuous distribution.
       Cannot be converted to readBit[]/acceptTest[] because the
       output is a continuous floating-point value, not a discrete
       bit.  The checker throws $dbc$cantHandle immediately.

   (B) AbsoluteTime -- a time-dependent call.
       Even if the result happens to be used for an acceptance
       decision, the value depends on wall-clock time and cannot
       be reproduced from a bit sequence.  CheckAlgorithmSafety
       flags this before the BFS even starts.

   Note: RandomInteger[{1,3}] (non-power-of-2 range) IS now supported
   via rejection sampling -- see examples/nonpower_of_two.wl.

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

(* ---- (B) AbsoluteTime: time-dependent, flagged by safety check ---- *)
UnanalyzableTime[state_Integer, readBit_, acceptTest_] := Module[
  {dir, nbr, dE},
  dir = readBit[];
  nbr = Mod[state + If[dir == 1, 1, -1] - 1, L$ua] + 1;
  dE  = energy$ua[nbr] - energy$ua[state];
  (* Using wall-clock time as a random source -- completely wrong *)
  If[Mod[Round[AbsoluteTime[] * 1000], 2] == 0,
    nbr, state]
]
