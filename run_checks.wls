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
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];

(* Load core library into global context -- no package prefix needed *)
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE 1: Kawasaki ring + Metropolis  [should PASS]
   ================================================================ *)
Print["\n\n", StringRepeat["#", 62]];
Print["EXAMPLE 1: Kawasaki ring + Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/ring_kawasaki.wl"];

result1 = RunFullCheck[
  allStates$rk,
  KawasakiSym$rk,   (* symbolic: uses MetropolisProb, keeps beta symbolic *)
  KawasakiNum$rk,   (* numeric:  uses N[If[dE<=0,1,Exp[-dE]]], beta in dE *)
  symEnergy$rk,
  numEnergy$rk,
  "SystemName"  -> "Ring L=3, 1 particle -- Kawasaki + Metropolis",
  "MaxBitDepth" -> 10,
  "TimeLimit"   -> 60.,
  "NSteps"      -> 80000,
  "Verbose"     -> True
];

(* ================================================================
   EXAMPLE 2: Always accept  [should FAIL]
   ================================================================ *)
Print["\n\n", StringRepeat["#", 62]];
Print["EXAMPLE 2: Always accept, no Metropolis  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/always_accept.wl"];

result2 = RunFullCheck[
  allStates$aa,
  AlwaysAcceptSym,
  AlwaysAcceptNum,
  symEnergy$aa,
  numEnergy$aa,
  "SystemName"  -> "Ring L=3, 1 particle -- always accept",
  "MaxBitDepth" -> 10,
  "TimeLimit"   -> 60.,
  "NSteps"      -> 80000,
  "Verbose"     -> False
];

(* ================================================================
   EXAMPLE 3: Wrong acceptance sign  [should FAIL]
   ================================================================ *)
Print["\n\n", StringRepeat["#", 62]];
Print["EXAMPLE 3: Wrong acceptance sign  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/wrong_sign.wl"];

result3 = RunFullCheck[
  allStates$ws,
  WrongSignSym,
  WrongSignNum,
  symEnergy$ws,
  numEnergy$ws,
  "SystemName"  -> "Ring L=3, 1 particle -- wrong acceptance sign",
  "MaxBitDepth" -> 10,
  "TimeLimit"   -> 60.,
  "NSteps"      -> 80000,
  "Verbose"     -> False
];

Print["\nAll checks complete."];
Exit[0]
