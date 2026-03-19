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
     F. Non-power-of-2 RandomInteger via rejection sampling  -- PASSES
     G. Float literal acceptance probabilities  -- PASSES
     H. Unanalyzable (RandomVariate)  -- aborts with clear error
     I. Unanalyzable (AbsoluteTime)  -- aborts with clear error

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
   EXAMPLE F: Non-power-of-2 RandomInteger via rejection sampling  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE F: RandomInteger[{1,3}] via rejection sampling  (expected: PASS)"];
Print["  Note: 2 bits read, value=3 paths silently discarded"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/nonpower_of_two.wl"];

resultF = RunFullCheck[
  1,
  KawasakiNonPow2,
  energy$np,
  numBeta$np,
  "SystemName"    -> "Non-power-of-2 RandomInteger, L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/nonpower_of_two.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE G: Float literal acceptance probabilities  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE G: Float literal probs (0.5, 1.0)  (expected: PASS)"];
Print["  Note: Reals rationalised in acceptTest before symbolic check"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/float_energy.wl"];

resultG = RunFullCheck[
  1,
  KawasakiFloat,
  energy$fe,
  numBeta$fe,
  "SystemName"    -> "Float literals rationalised, L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/float_energy.wl"],
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE H: Unanalyzable -- RandomVariate (aborts cleanly)
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE H: RandomVariate  (expected: analysis ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable.wl"];

resultH = RunFullCheck[
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
Print["  Result[error] = ", resultH["error"]];

(* ================================================================
   EXAMPLE I: Unanalyzable -- AbsoluteTime (aborts cleanly)
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE I: AbsoluteTime  (expected: analysis ABORTED)"];
Print[StringRepeat["#", 62]];

resultI = RunFullCheck[
  1,
  UnanalyzableTime,
  energy$ua,
  numBeta$ua,
  "SystemName"    -> "Unanalyzable: AbsoluteTime",
  "MaxBitDepth"   -> 10,
  "NSteps"        -> 10000,
  "Verbose"       -> True,
  "OpenWindow"    -> False
];
Print["  Result[error] = ", resultI["error"]];

Print["\n", StringRepeat["=", 62]];
Print["Extended examples complete."];
Print[StringRepeat["=", 62]];
