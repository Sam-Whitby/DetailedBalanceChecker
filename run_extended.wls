#!/usr/bin/env wolframscript
(* ================================================================
   run_extended.wls  --  Extended examples
   ================================================================
   Demonstrates the checker on a range of systems and features:

     A. Cyclic drift  -- FAILS DB, PASSES numerically
     B. Multi-particle Kawasaki  -- PASSES
     C. Continuous (rational) states  -- PASSES
     D. Native RandomReal[] usage, correct  -- PASSES
     E. Native RandomReal[] usage, wrong sign  -- FAILS
     F. Unanalyzable (RandomVariate)  -- aborts with clear error
     G. Unanalyzable (non-power-of-2 range)  -- aborts with clear error

   Run from terminal:
     wolframscript -file run_extended.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE A: Cyclic drift -- FAILS DB, PASSES numerically
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE A: Cyclic drift (FAILS DB, PASSES numerically)"];
Print["  Expected: symbolic FAIL, numerical KL ~ 0"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/cyclic_drift.wl"];

resultA = RunFullCheck[
  1,
  CyclicDrift,
  energy$cd,
  numBeta$cd,
  "SystemName"    -> "Cyclic drift L=4  [DB-FAIL, numerically PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/cyclic_drift.wl"],
  "MaxBitDepth"   -> 5,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE B: Multi-particle Kawasaki on L=4 ring  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE B: Two-particle Kawasaki, L=4 ring  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/multi_particle.wl"];

resultB = RunFullCheck[
  {1, 2},             (* seed: particle 1 at site 1, particle 2 at site 2 *)
  KawasakiMulti,
  energy$mp,
  numBeta$mp,
  "SystemName"    -> "Two-particle Kawasaki, L=4  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/multi_particle.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 120000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE C: Metropolis on discretised continuous ring  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE C: Continuous (rational) states  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/continuous_metropolis.wl"];

resultC = RunFullCheck[
  0,                  (* seed: position 0 on circle *)
  ContinuousMetropolis,
  energy$cm,
  numBeta$cm,
  "SystemName"    -> "Metropolis on discretised circle  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/continuous_metropolis.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE D: Native RandomReal[]/RandomInteger[] -- correct  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE D: Native RandomReal[], correct Metropolis  (expected: PASS)"];
Print["  Note: RandomInteger[1] and RandomReal[] < p are intercepted"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/random_native_pass.wl"];

resultD = RunFullCheck[
  1,
  KawasakiNative,
  energy$rn,
  numBeta$rn,
  "SystemName"    -> "Native random calls, Kawasaki ring L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/random_native_pass.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE E: Native RandomReal[] -- wrong acceptance sign  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE E: Native RandomReal[], wrong acceptance sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/random_native_fail.wl"];

resultE = RunFullCheck[
  1,
  KawasakiNativeWrongSign,
  energy$rn,
  numBeta$rn,
  "SystemName"    -> "Native random calls, wrong sign  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/random_native_fail.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE F: Unanalyzable -- RandomVariate (aborts cleanly)
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE F: RandomVariate  (expected: analysis ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable.wl"];

resultF = RunFullCheck[
  1,
  UnanalyzableVariate,
  energy$ua,
  numBeta$ua,
  "SystemName"    -> "Unanalyzable: RandomVariate",
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 10000,
  "Verbose"       -> True,
  "OpenWindow"    -> False
];
Print["  Result[error] = ", resultF["error"]];

(* ================================================================
   EXAMPLE G: Unanalyzable -- non-power-of-2 RandomInteger range
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE G: RandomInteger[{1,3}]  (expected: analysis ABORTED)"];
Print[StringRepeat["#", 62]];

resultG = RunFullCheck[
  1,
  UnanalyzableIntRange,
  energy$ua,
  numBeta$ua,
  "SystemName"    -> "Unanalyzable: RandomInteger non-power-of-2",
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 10000,
  "Verbose"       -> True,
  "OpenWindow"    -> False
];
Print["  Result[error] = ", resultG["error"]];

Print["\n", StringRepeat["=", 62]];
Print["Extended examples complete."];
Print[StringRepeat["=", 62]];
