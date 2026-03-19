#!/usr/bin/env wolframscript
(* ================================================================
   run_extended.wls  --  Extended examples
   ================================================================
   A. Cyclic drift            -- FAILS DB, PASSES numerically
   B. Multi-particle Kawasaki -- PASSES
   C. Continuous states       -- PASSES
   D. Native RandomReal[], correct  -- PASSES
   E. Native RandomReal[], wrong sign  -- FAILS
   F. RandomInteger[{1,3}] rejection sampling  -- PASSES
   G. Float energy values + float acceptance probs  -- PASSES
   H. Unanalyzable algorithm (RandomVariate)  -- aborts
   I. Unanalyzable algorithm (AbsoluteTime)   -- aborts
   J. Unanalyzable energy (AbsoluteTime)      -- aborts
   K. Unanalyzable energy (RandomVariate)     -- aborts

   Run from terminal:
     wolframscript -file run_extended.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE A: Cyclic drift
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE A: Cyclic drift  (expected: DB FAIL, numerical PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/cyclic_drift.wl"];
resultA = RunFullCheck[
  1, CyclicDrift, energy$cd, numBeta$cd,
  "SystemName"    -> "Cyclic drift L=4  [DB-FAIL, numerically PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/cyclic_drift.wl"],
  "MaxBitDepth"   -> 5, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE B: Two-particle Kawasaki
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE B: Two-particle Kawasaki, L=4  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/multi_particle.wl"];
resultB = RunFullCheck[
  {1, 2}, KawasakiMulti, energy$mp, numBeta$mp,
  "SystemName"    -> "Two-particle Kawasaki, L=4  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/multi_particle.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 120000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE C: Continuous (rational) states
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE C: Continuous (rational) states  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/continuous_metropolis.wl"];
resultC = RunFullCheck[
  0, ContinuousMetropolis, energy$cm, numBeta$cm,
  "SystemName"    -> "Metropolis on discretised circle  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/continuous_metropolis.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE D: Native RandomReal[], correct
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE D: Native RandomReal[], correct Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/random_native_pass.wl"];
resultD = RunFullCheck[
  1, KawasakiNative, energy$rn, numBeta$rn,
  "SystemName"    -> "Native random calls, Kawasaki ring L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/random_native_pass.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE E: Native RandomReal[], wrong acceptance sign
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE E: Native RandomReal[], wrong sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/random_native_fail.wl"];
resultE = RunFullCheck[
  1, KawasakiNativeWrongSign, energy$rn, numBeta$rn,
  "SystemName"    -> "Native random calls, wrong sign  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/random_native_fail.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> False, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE F: Non-power-of-2 RandomInteger via rejection sampling
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE F: RandomInteger[{1,3}] rejection sampling  (expected: PASS)"];
Print["  Note: 2 bits read; value=3 paths silently discarded"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/nonpower_of_two.wl"];
resultF = RunFullCheck[
  1, KawasakiNonPow2, energy$np, numBeta$np,
  "SystemName"    -> "Non-power-of-2 RandomInteger, L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/nonpower_of_two.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE G: Float energy values and float acceptance probabilities
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE G: Float energy + float probs  (expected: PASS)"];
Print["  Note: floats rationalised in acceptTest and CheckDetailedBalance"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/float_energy.wl"];
resultG = RunFullCheck[
  1, KawasakiFloat, energy$fe, numBeta$fe,
  "SystemName"    -> "Float energy + probs, L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/float_energy.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE H: Unanalyzable algorithm -- RandomVariate
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE H: RandomVariate in algorithm  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable.wl"];
resultH = RunFullCheck[
  1, UnanalyzableVariate, energy$ua, numBeta$ua,
  "SystemName" -> "Unanalyzable: RandomVariate",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultH["error"]];

(* ================================================================
   EXAMPLE I: Unanalyzable algorithm -- AbsoluteTime
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE I: AbsoluteTime in algorithm  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

resultI = RunFullCheck[
  1, UnanalyzableTime, energy$ua, numBeta$ua,
  "SystemName" -> "Unanalyzable: AbsoluteTime",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultI["error"]];

(* ================================================================
   EXAMPLE J: Unanalyzable energy -- AbsoluteTime
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE J: AbsoluteTime in energy function  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable_energy.wl"];
resultJ = RunFullCheck[
  1, KawasakiUE, energyTime$ue, numBeta$ue,
  "SystemName" -> "Unanalyzable energy: AbsoluteTime",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultJ["error"]];

(* ================================================================
   EXAMPLE K: Unanalyzable energy -- RandomVariate
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE K: RandomVariate in energy function  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

resultK = RunFullCheck[
  1, KawasakiUE, energyRand$ue, numBeta$ue,
  "SystemName" -> "Unanalyzable energy: RandomVariate",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultK["error"]];

Print["\n", StringRepeat["=", 62]];
Print["Extended examples complete."];
Print[StringRepeat["=", 62]];
