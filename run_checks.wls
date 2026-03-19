#!/usr/bin/env wolframscript
(* ================================================================
   run_checks.wls  --  Detailed Balance Checker  --  main script
   ================================================================
   Run from terminal:
     wolframscript -file run_checks.wls

   Three examples are run in sequence.  Each opens a Mathematica
   window showing:
     - the algorithm source
     - the decision tree for every starting state
     - the symbolic transition matrix
     - the detailed-balance pair table (green = pass, red = fail)
     - the simulated vs Boltzmann bar chart
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 1: Ring Kawasaki + Metropolis  [should PASS]
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE 1: Kawasaki ring + Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/ring_kawasaki.wl"];

result1 = RunFullCheck[
  allStates$rk,
  KawasakiSym$rk,
  KawasakiNum$rk,
  symEnergy$rk,
  numEnergy$rk,
  "SystemName"    -> "Kawasaki ring L=3, 1 particle  [PASS expected]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/ring_kawasaki.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 2: Always accept  [should FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 2: Always accept, no Metropolis  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/always_accept.wl"];

result2 = RunFullCheck[
  allStates$aa,
  AlwaysAcceptSym,
  AlwaysAcceptNum,
  symEnergy$aa,
  numEnergy$aa,
  "SystemName"    -> "Always accept (no Metropolis)  [FAIL expected]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/always_accept.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 3: Wrong acceptance sign  [should FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 3: Wrong acceptance sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/wrong_sign.wl"];

result3 = RunFullCheck[
  allStates$ws,
  WrongSignSym,
  WrongSignNum,
  symEnergy$ws,
  numEnergy$ws,
  "SystemName"    -> "Wrong acceptance sign (exp[+beta*dE])  [FAIL expected]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/wrong_sign.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

Print["\nAll checks complete.  PNG reports have been opened."];
