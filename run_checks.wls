#!/usr/bin/env wolframscript
(* ================================================================
   run_checks.wls  --  Core examples (ring Kawasaki variants)
   ================================================================
   Three examples on an L=3 ring, demonstrating correct and
   incorrect Metropolis acceptance criteria.

   Run from terminal:
     wolframscript -file run_checks.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 1: Ring Kawasaki + Metropolis  [PASS]
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE 1: Kawasaki ring + Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/ring_kawasaki.wl"];

result1 = RunFullCheck[
  1,
  KawasakiRing,
  energy$rk,
  numBeta$rk,
  "SysParams"     -> params$rk,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Kawasaki ring L=3, Metropolis  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/ring_kawasaki.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 2: Always accept (no Metropolis)  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 2: Always accept, no Metropolis  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/always_accept.wl"];

result2 = RunFullCheck[
  1,
  AlwaysAccept,
  energy$rk,
  numBeta$rk,
  "SysParams"     -> params$rk,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Always accept (no Metropolis)  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/always_accept.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 3: Wrong acceptance sign  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 3: Wrong acceptance sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/wrong_sign.wl"];

result3 = RunFullCheck[
  1,
  WrongSign,
  energy$rk,
  numBeta$rk,
  "SysParams"     -> params$rk,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Wrong acceptance sign (Exp[+beta*dE])  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/wrong_sign.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

Print["\nAll checks complete.  PNG reports have been opened."];
