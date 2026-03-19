#!/usr/bin/env wolframscript
(* ================================================================
   check.wls  --  Run a single detailed-balance check from terminal
   ================================================================
   Usage:
     wolframscript -file check.wls \
       <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> \
       [NSteps=100000] [MaxBitDepth=20] [SystemName="My algorithm"] \
       [Verbose=True] [OpenWindow=True]

   Arguments (positional, required):
     alg_file      Path to your .wl file defining the algorithm and energy.
                   Absolute path, or relative to the DetailedBalanceChecker dir.
     seed_state    A valid starting state (integer, list, etc.)
     alg_name      Name of the algorithm function defined in alg_file
     energy_name   Name of the energy function defined in alg_file
     num_beta      Inverse temperature β for the numerical MCMC run

   Named options (optional, key=value, no spaces around =):
     NSteps=N           Number of MCMC steps (default 100000)
     MaxBitDepth=N      Max BFS depth (default 20)
     SystemName=label   Label shown in the PNG report (quote if spaces)
     Verbose=True/False Print BFS progress per state (default True)
     OpenWindow=True/False  Open PNG report when done (default True)

   Examples:
     wolframscript -file check.wls examples/ring_kawasaki.wl 1 KawasakiRing energy 1
     wolframscript -file check.wls my_alg.wl 1 MyAlg energy 2.0 NSteps=80000 Verbose=False
     wolframscript -file check.wls examples/multi_particle.wl "{1,2}" KawasakiMulti energy$mp 1 NSteps=120000
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;

(* Drop the script name itself (first element) *)
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "check.wls"],
  args = Rest[args]
];

If[Length[args] < 5,
  Print["Usage: wolframscript -file check.wls <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> [options]"];
  Print[""];
  Print["Required arguments:"];
  Print["  alg_file     Path to .wl file containing algorithm and energy"];
  Print["  seed_state   Starting state (integer, or {i,j} for multi-particle)"];
  Print["  alg_name     Name of the algorithm function"];
  Print["  energy_name  Name of the energy function"];
  Print["  num_beta     Inverse temperature for numerical MCMC"];
  Print[""];
  Print["Optional key=value options:"];
  Print["  NSteps=N           (default 100000)"];
  Print["  MaxBitDepth=N      (default 20)"];
  Print["  SystemName=label   (default: alg_name)"];
  Print["  Verbose=True/False (default True)"];
  Print["  OpenWindow=True/False (default True)"];
  Exit[1]
];

algFile    = args[[1]];
seedStr    = args[[2]];
algName    = args[[3]];
energyName = args[[4]];
numBetaStr = args[[5]];
extraArgs  = If[Length[args] >= 6, args[[6;;]], {}];

(* Resolve alg_file: try absolute first, then relative to script dir *)
If[!FileExistsQ[algFile],
  algFile = scriptDir <> algFile
];
If[!FileExistsQ[algFile],
  Print["ERROR: Algorithm file not found: ", args[[1]]];
  Exit[1]
];

(* Parse seed state: use ToExpression to handle integers, lists, etc. *)
seedState = ToExpression[seedStr];
If[seedState === $Failed,
  Print["ERROR: Cannot parse seed state: ", seedStr];
  Exit[1]
];

(* Parse num_beta *)
numBeta = ToExpression[numBetaStr];
If[!NumberQ[numBeta],
  Print["ERROR: num_beta must be a number, got: ", numBetaStr];
  Exit[1]
];

(* Parse optional key=value arguments *)
nSteps      = 100000;
maxBitDepth = 20;
systemName  = algName;
verbose     = True;
openWindow  = True;

Do[
  If[StringMatchQ[opt, "NSteps=" ~~ __],
    nSteps = ToExpression[StringDrop[opt, 7]]];
  If[StringMatchQ[opt, "MaxBitDepth=" ~~ __],
    maxBitDepth = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "SystemName=" ~~ __],
    systemName = StringDrop[opt, 11]];
  If[StringMatchQ[opt, "Verbose=" ~~ __],
    verbose = ToExpression[StringDrop[opt, 8]]];
  If[StringMatchQ[opt, "OpenWindow=" ~~ __],
    openWindow = ToExpression[StringDrop[opt, 11]]],
  {opt, extraArgs}
];

(* Load the algorithm file *)
Print["Loading: ", algFile];
Get[algFile];

(* Look up algorithm and energy by name *)
alg    = ToExpression[algName];
energy = ToExpression[energyName];

If[!MatchQ[alg, _Symbol],
  Print["ERROR: Algorithm '", algName, "' not found after loading ", algFile];
  Exit[1]
];
If[!MatchQ[energy, _Symbol],
  Print["ERROR: Energy function '", energyName, "' not found after loading ", algFile];
  Exit[1]
];

(* Run the check *)
Print[""];
Print[StringRepeat["=", 62]];
Print["System:     ", systemName];
Print["Algorithm:  ", algName];
Print["Energy:     ", energyName];
Print["Seed state: ", seedState];
Print["numBeta:    ", numBeta];
Print["NSteps:     ", nSteps];
Print["MaxBitDepth:", maxBitDepth];
Print[StringRepeat["=", 62]];
Print[""];

result = RunFullCheck[
  seedState, alg, energy, numBeta,
  "SystemName"   -> systemName,
  "NSteps"       -> nSteps,
  "MaxBitDepth"  -> maxBitDepth,
  "Verbose"      -> verbose,
  "OpenWindow"   -> openWindow
];

Print[""];
Print[StringRepeat["=", 62]];
If[KeyExistsQ[result, "error"],
  Print["RESULT: ABORTED -- ", result["error"]],
  Print["RESULT: ", If[result["pass"], "PASS", "FAIL"],
        "  (KL divergence: ", result["klDivergence"], ")"]
];
Print[StringRepeat["=", 62]];
