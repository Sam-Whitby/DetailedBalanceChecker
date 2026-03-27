#!/usr/bin/env wolframscript
(* ================================================================
   run_variable_bit.wls  --  Variable-bit-length examples
   ================================================================
   Demonstrates the checker on algorithms whose decision trees have
   leaves at DIFFERENT depths (some paths consume 2 bits, others 3).

   Run from terminal:
     wolframscript -file run_variable_bit.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 4: Two-speed Metropolis  [PASS]
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE 4: Variable-bit two-speed Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/variable_bit_pass.wl"];

result4 = RunFullCheck[
  1,
  TwoSpeed,
  energy$vb,
  numBeta$vb,
  "SysParams"     -> params$vb,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Variable-bit two-speed Metropolis, L=4  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/variable_bit_pass.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 120000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE 5: Biased rightward drift  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 5: Variable-bit biased drift  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/variable_bit_fail.wl"];

result5 = RunFullCheck[
  1,
  Biased,
  energy$vb,
  numBeta$vb,
  "SysParams"     -> params$vb,
  "NumericSeed"   -> 1,
  "SystemName"    -> "Variable-bit biased rightward drift, L=4  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/variable_bit_fail.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 120000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

Print["\nAll checks complete.  Two PNG reports have been opened."];
