#!/usr/bin/env wolframscript
(* ================================================================
   run_extended.wls  --  Extended examples
   ================================================================
   A. Cyclic drift            -- FAILS DB, PASSES numerically
   B. Two-particle Kawasaki + NN coupling -- PASSES
   C. Continuous states       -- PASSES
   D. Native RandomReal[], correct  -- PASSES
   E. Native RandomReal[], wrong sign  -- FAILS
   F. RandomInteger[{1,3}] rejection sampling  -- PASSES
   G. Two-particle Kawasaki + NN coupling, wrong sign  -- FAILS
   H. Three-particle Kawasaki + NN+NNN coupling  -- PASSES
   I. Unanalyzable algorithm (RandomVariate)  -- aborts
   J. Unanalyzable algorithm (AbsoluteTime)   -- aborts
   K. Unanalyzable energy (AbsoluteTime)      -- aborts
   L. Unanalyzable energy (RandomVariate)     -- aborts

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
   EXAMPLE B: Two-particle Kawasaki + nearest-neighbour coupling
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE B: Two-particle Kawasaki + NN coupling, L=4  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/multi_particle.wl"];
resultB = RunFullCheck[
  {1, 2}, KawasakiMulti, energy$mp, numBeta$mp,
  "SysParams"     -> params$mp,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Two-particle Kawasaki + NN coupling, L=4  [PASS]",
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
  "SysParams"     -> params$rn,
  "NumericSeed"   -> 1,
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
  "SysParams"     -> params$rn,
  "NumericSeed"   -> 1,
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
  "SysParams"     -> params$np,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Non-power-of-2 RandomInteger, L=3  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/nonpower_of_two.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 80000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE G: Two-particle Kawasaki + NN coupling, wrong sign  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE G: Two-particle + NN coupling, wrong sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

(* Reuses params from multi_particle.wl (already loaded above).
   Same system, but acceptance uses +beta*dE (wrong sign). *)
KawasakiMultiWrongSign[{p1_Integer, p2_Integer}] :=
  Module[
    {b1, b2, mover, other, target, dE,
     rightOf = Function[s, Mod[s, L$mp] + 1],
     leftOf  = Function[s, Mod[s - 2, L$mp] + 1]},
    b1 = RandomInteger[];
    {mover, other} = If[b1 == 0, {p1, p2}, {p2, p1}];
    b2     = RandomInteger[];
    target = If[b2 == 0, leftOf[mover], rightOf[mover]];
    If[target == other,
      {p1, p2},
      dE = energy$mp[Sort[{target, other}]] - energy$mp[{p1, p2}];
      (* BUG: wrong sign -- preferentially accepts uphill moves *)
      If[RandomReal[] < Piecewise[{{1, dE >= 0}, {Exp[\[Beta]*dE], dE < 0}}],
        Sort[{target, other}],
        {p1, p2}]
    ]
  ]

resultG = RunFullCheck[
  {1, 2}, KawasakiMultiWrongSign, energy$mp, numBeta$mp,
  "SysParams"     -> params$mp,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Two-particle + NN coupling, wrong sign  [FAIL]",
  "MaxBitDepth"   -> 10, "NSteps" -> 120000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE H: Three-particle Kawasaki + NN+NNN coupling
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE H: Three-particle + NN+NNN coupling, L=4  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/coupling_3particle.wl"];
resultH = RunFullCheck[
  {1, 2, 3}, KawasakiCoupled3, energy$c3, numBeta$c3,
  "SysParams"     -> params$c3,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Three-particle Kawasaki + NN+NNN coupling, L=4  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/coupling_3particle.wl"],
  "MaxBitDepth"   -> 10, "NSteps" -> 120000,
  "Verbose" -> True, "OpenWindow" -> True
];

(* ================================================================
   EXAMPLE I: Unanalyzable algorithm -- RandomVariate
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE I: RandomVariate in algorithm  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable.wl"];
resultI = RunFullCheck[
  1, UnanalyzableVariate, energy$ua, numBeta$ua,
  "SystemName" -> "Unanalyzable: RandomVariate",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultI["error"]];

(* ================================================================
   EXAMPLE J: Unanalyzable algorithm -- AbsoluteTime
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE J: AbsoluteTime in algorithm  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

resultJ = RunFullCheck[
  1, UnanalyzableTime, energy$ua, numBeta$ua,
  "SystemName" -> "Unanalyzable: AbsoluteTime",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultJ["error"]];

(* ================================================================
   EXAMPLE K: Unanalyzable energy -- AbsoluteTime
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE K: AbsoluteTime in energy function  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/unanalyzable_energy.wl"];
resultK = RunFullCheck[
  1, KawasakiUE, energyTime$ue, numBeta$ue,
  "SystemName" -> "Unanalyzable energy: AbsoluteTime",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultK["error"]];

(* ================================================================
   EXAMPLE L: Unanalyzable energy -- RandomVariate
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE L: RandomVariate in energy function  (expected: ABORTED)"];
Print[StringRepeat["#", 62]];

resultL = RunFullCheck[
  1, KawasakiUE, energyRand$ue, numBeta$ue,
  "SystemName" -> "Unanalyzable energy: RandomVariate",
  "MaxBitDepth" -> 10, "NSteps" -> 10000,
  "Verbose" -> True, "OpenWindow" -> False
];
Print["  Result[error] = ", resultL["error"]];

Print["\n", StringRepeat["=", 62]];
Print["Extended examples complete."];
Print[StringRepeat["=", 62]];
