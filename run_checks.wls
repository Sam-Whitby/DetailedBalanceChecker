#!/usr/bin/env wolframscript
(* ================================================================
   run_checks.wls  --  Detailed Balance Checker  --  main script
   ================================================================
   Run from terminal:
     wolframscript -file run_checks.wls

   Three examples:
     1. Kawasaki ring with Metropolis          --> PASS
     2. Kawasaki ring, always accept           --> FAIL
     3. Kawasaki ring, wrong acceptance sign   --> FAIL

   For each example, two PNG plots are exported to the same
   directory as this script:
     {SystemName}_frequencies.png       simulated vs Boltzmann bar chart
     {SystemName}_transition_matrix.png numerical T[from,to] heatmap
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 1: Kawasaki ring + Metropolis  [should PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE 1: Kawasaki ring + Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/ring_kawasaki.wl"];

result1 = RunFullCheck[
  allStates$rk,
  KawasakiSym$rk,
  KawasakiNum$rk,
  symEnergy$rk,
  numEnergy$rk,
  "SystemName"    -> "Kawasaki_ring_PASS",
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "ExportPlots"   -> True,
  "PlotDirectory" -> scriptDir
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
  "SystemName"    -> "Always_accept_FAIL",
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "ExportPlots"   -> True,
  "PlotDirectory" -> scriptDir
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
  "SystemName"    -> "Wrong_sign_FAIL",
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "ExportPlots"   -> True,
  "PlotDirectory" -> scriptDir
];

Print["\nAll checks complete."];
Print["Plots saved in: ", scriptDir];
Exit[0]
