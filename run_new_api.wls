#!/usr/bin/env wolframscript
(* ================================================================
   run_new_api.wls  --  Generalised API demonstration
   ================================================================
   Four examples using the new three-point interface:

     RunFullCheck[seedState, alg, energy, numBeta, options...]

   States are discovered automatically from the seed; no allStates
   list is needed.  The single algorithm uses acceptTest[p] for
   all accept/reject decisions (no separate sym/num versions).
   A single bare energy function is shared throughout.

   Run from terminal:
     wolframscript -file run_new_api.wls
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ================================================================
   EXAMPLE A: Kawasaki ring, correct Metropolis  [PASS]
   ================================================================ *)
Print[StringRepeat["#", 62]];
Print["EXAMPLE A: Kawasaki ring, Metropolis  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/kawasaki_new.wl"];

resultA = RunFullCheck[
  1,
  KawasakiNew,
  energy$kn,
  numBeta$kn,
  "SystemName"    -> "Kawasaki ring L=3, Metropolis  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/kawasaki_new.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> True,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE B: Always-accept ring (no Metropolis)  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE B: Always-accept ring  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

resultB = RunFullCheck[
  1,
  AlwaysAcceptNew,
  energy$kn,
  numBeta$kn,
  "SystemName"    -> "Always-accept ring L=3  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/kawasaki_new.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE C: Barker criterion on L=4 ring  [PASS]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE C: Barker (heat-bath) criterion, L=4 ring  (expected: PASS)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/barker_ring.wl"];

resultC = RunFullCheck[
  1,
  BarkerRing,
  energy$br,
  numBeta$br,
  "SystemName"    -> "Barker criterion ring L=4  [PASS]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/barker_ring.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

(* ================================================================
   EXAMPLE D: Asymmetric proposal + Metropolis on L=4 ring  [FAIL]
   ================================================================ *)
Print["\n", StringRepeat["#", 62]];
Print["EXAMPLE D: Asymmetric proposal + Metropolis  (expected: FAIL)"];
Print[StringRepeat["#", 62]];

Get[scriptDir <> "examples/asymmetric_proposal.wl"];

resultD = RunFullCheck[
  1,
  AsymmetricProposal,
  energy$br,
  numBeta$br,
  "SystemName"    -> "Asymmetric proposal + Metropolis L=4 ring  [FAIL]",
  "AlgorithmCode" -> ReadString[scriptDir <> "examples/asymmetric_proposal.wl"],
  "MaxBitDepth"   -> 10,
  "TimeLimit"     -> 60.,
  "NSteps"        -> 80000,
  "Verbose"       -> False,
  "OpenWindow"    -> True
];

Print["\nAll checks complete.  PNG reports have been opened."];
