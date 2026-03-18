#!/usr/bin/env wolframscript
(* ================================================================
   run_checks.wls  --  Detailed Balance Checker  --  main script
   ================================================================
   Run from terminal:
     wolframscript -file run_checks.wls

   This loads the core library and runs three examples:
     1. Kawasaki ring with Metropolis          --> PASS
     2. Kawasaki ring, always accept           --> FAIL
     3. Kawasaki ring, wrong acceptance sign   --> FAIL
   ================================================================ *)

(* ---- Locate files relative to this script ---- *)
scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* Make core functions available without package prefix *)
MetropolisProb    = DetailedBalanceChecker`Private`MetropolisProb;
RunWithBits       = DetailedBalanceChecker`Private`RunWithBits;
BuildTransitionMatrix  = DetailedBalanceChecker`Private`BuildTransitionMatrix;
CheckDetailedBalance   = DetailedBalanceChecker`Private`CheckDetailedBalance;
RunNumericalMCMC       = DetailedBalanceChecker`Private`RunNumericalMCMC;
BoltzmannWeights       = DetailedBalanceChecker`Private`BoltzmannWeights;
RunFullCheck           = DetailedBalanceChecker`Private`RunFullCheck;

(* ================================================================
   EXAMPLE 1: Ring Kawasaki with Metropolis  [should PASS]
   ================================================================ *)
Print["\n\n"];
Print[StringRepeat["#", 62]];
Print["EXAMPLE 1: Ring Kawasaki with Metropolis (expected PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/ring_kawasaki.wl"];

(* Run the full check -- symbolic energies keep beta symbolic *)
result1 = RunFullCheck[
  allStates$rk,
  KawasakiGood,          (* symbolic: uses MetropolisProb *)
  symEnergy$rk,
  numEnergy$rk,          (* numeric: specific beta and eps values *)
  "SystemName"   -> "Ring L=3, 1 particle -- Kawasaki + Metropolis",
  "MaxBitDepth"  -> 10,
  "TimeLimit"    -> 60.,
  "NSteps"       -> 80000,
  "Verbose"      -> True
];

(* ================================================================
   EXAMPLE 2: Always accept (no Metropolis)  [should FAIL]
   ================================================================ *)
Print["\n\n"];
Print[StringRepeat["#", 62]];
Print["EXAMPLE 2: Always accept -- no Metropolis (expected FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/always_accept.wl"];

result2 = RunFullCheck[
  allStates$aa,
  AlwaysAccept,
  symEnergy$aa,
  numEnergy$aa,
  "SystemName"   -> "Ring L=3, 1 particle -- always accept",
  "MaxBitDepth"  -> 10,
  "TimeLimit"    -> 60.,
  "NSteps"       -> 80000,
  "Verbose"      -> False
];

(* ================================================================
   EXAMPLE 3: Wrong acceptance sign  [should FAIL]
   ================================================================ *)
Print["\n\n"];
Print[StringRepeat["#", 62]];
Print["EXAMPLE 3: Wrong acceptance sign (expected FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/wrong_sign.wl"];

result3 = RunFullCheck[
  allStates$ws,
  WrongSign,
  symEnergy$ws,
  numEnergy$ws,
  "SystemName"   -> "Ring L=3, 1 particle -- wrong acceptance sign",
  "MaxBitDepth"  -> 10,
  "TimeLimit"    -> 60.,
  "NSteps"       -> 80000,
  "Verbose"      -> False
];

Print["All checks complete."];
Exit[0]
