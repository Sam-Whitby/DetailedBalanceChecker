#!/usr/bin/env wolframscript
(* ================================================================
   run_variable_bit.wls  --  Variable-bit-length example runner
   ================================================================
   Demonstrates the checker on two algorithms whose decision trees
   have leaves at DIFFERENT depths (some paths consume 2 bits,
   others 3), something classical analysis cannot handle directly.

   Run from terminal:
     wolframscript -file run_variable_bit.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 4: Two-speed Metropolis  [should PASS]
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE 4: Variable-bit two-speed Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/variable_bit_pass.wl"];

result4 = RunFullCheck[
  allStates$vb,
  TwoSpeedSym,
  TwoSpeedNum,
  symEnergy$vb,
  numEnergy$vb,
  "SystemName"    -> "Variable-bit two-speed Metropolis, L=4 ring  [PASS expected]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/variable_bit_pass.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 120000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 5: Biased rightward drift  [should FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 5: Variable-bit biased drift  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/variable_bit_fail.wl"];

result5 = RunFullCheck[
  allStates$vb,
  BiasedSym,
  BiasedNum,
  symEnergy$vb,
  numEnergy$vb,
  "SystemName"    -> "Variable-bit biased rightward drift, L=4 ring  [FAIL expected]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/variable_bit_fail.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 120000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

Print["\nAll checks complete.  Two PNG reports have been opened."];
