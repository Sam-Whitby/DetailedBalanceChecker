#!/usr/bin/env wolframscript
(* ================================================================
   check.wls  --  Run a single detailed-balance check from terminal
   ================================================================
   TWO MODES:

   NEW (single-file) mode -- recommended:
     wolframscript -file check.wls <alg_file> [options]

     The .wl file must define:
       energy[state_]    -- bare energy function (no beta)
       Algorithm[state_] -- MCMC move function
       seedState         -- a valid starting state
       numBeta           -- inverse temperature (numeric)
     Optionally define:
       symParams = <|"eps" -> {...}, "couplings" -> {...}|>
     for symbolic energy parameters (enables exact symbolic check).

   LEGACY mode -- for existing examples:
     wolframscript -file check.wls \
       <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> \
       [options]

   Named options (key=value, no spaces around =):
     NSteps=N           Number of MCMC steps (default 100000)
     MaxBitDepth=N      Max BFS depth (default 20)
     SystemName=label   Label in the PNG report
     Verbose=True/False Print BFS progress (default True)
     OpenWindow=True/False  Open PNG report when done (default True)
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "check.wls"],
  args = Rest[args]
];

(* Separate positional args from key=value options *)
posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args, StringMatchQ[#,  __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1,
  Print["Usage (new mode):    wolframscript -file check.wls <alg_file> [options]"];
  Print["Usage (legacy mode): wolframscript -file check.wls <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> [options]"];
  Print[""];
  Print["New mode: .wl file must define Algorithm, energy, seedState, numBeta."];
  Print["          Optionally define symParams = <|\"eps\" -> {...}, ...|>."];
  Print[""];
  Print["Options (key=value):"];
  Print["  NSteps=N           (default 100000)"];
  Print["  MaxBitDepth=N      (default 20)"];
  Print["  SystemName=label   (default: algorithm name)"];
  Print["  Verbose=True/False (default True)"];
  Print["  OpenWindow=True/False (default True)"];
  Exit[1]
];

(* ---- Resolve algorithm file path ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile],
  algFile = scriptDir <> algFile
];
If[!FileExistsQ[algFile],
  Print["ERROR: Algorithm file not found: ", posArgs[[1]]];
  Exit[1]
];

(* ---- Parse optional key=value arguments ---- *)
nSteps      = 100000;
maxBitDepth = 20;
systemName  = ""; (* filled in after loading *)
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
  {opt, kvArgs}
];

(* ================================================================
   NEW MODE: single positional argument, standard names in .wl file
   ================================================================ *)
If[Length[posArgs] == 1,

  (* Clear standard names to avoid stale definitions *)
  Remove["Global`Algorithm"];
  Remove["Global`energy"];
  Remove["Global`seedState"];
  Remove["Global`numBeta"];
  Remove["Global`symParams"];

  Print["Loading: ", algFile];
  Get[algFile];

  (* Validate required definitions *)
  If[!ValueQ[Global`seedState],
    Print["ERROR: 'seedState' not defined in ", posArgs[[1]]];
    Print["  Your .wl file must define: seedState = <a valid starting state>"];
    Exit[1]
  ];
  If[!ValueQ[Global`numBeta] || !NumericQ[Global`numBeta],
    Print["ERROR: 'numBeta' (numeric inverse temperature) not defined in ", posArgs[[1]]];
    Exit[1]
  ];
  If[Length[DownValues[Global`Algorithm]] == 0,
    Print["ERROR: 'Algorithm' function not defined in ", posArgs[[1]]];
    Exit[1]
  ];
  If[Length[DownValues[Global`energy]] == 0,
    Print["ERROR: 'energy' function not defined in ", posArgs[[1]]];
    Exit[1]
  ];

  $seedState  = Global`seedState;
  $numBeta    = Global`numBeta;
  $alg        = Global`Algorithm;
  $energy     = Global`energy;
  $sysParams  = If[ValueQ[Global`symParams], Global`symParams, None];
  $algName    = "Algorithm";
  If[systemName == "", systemName = FileBaseName[posArgs[[1]]]];

  ,

(* ================================================================
   LEGACY MODE: five positional arguments
   ================================================================ *)

  If[Length[posArgs] < 5,
    Print["ERROR: Legacy mode requires 5 positional arguments."];
    Print["Usage: wolframscript -file check.wls <alg_file> <seed_state> <alg_name> <energy_name> <num_beta> [options]"];
    Exit[1]
  ];

  seedStr    = posArgs[[2]];
  algName    = posArgs[[3]];
  energyName = posArgs[[4]];
  numBetaStr = posArgs[[5]];

  $seedState = ToExpression[seedStr];
  If[$seedState === $Failed,
    Print["ERROR: Cannot parse seed state: ", seedStr]; Exit[1]
  ];
  $numBeta = ToExpression[numBetaStr];
  If[!NumericQ[$numBeta],
    Print["ERROR: num_beta must be a number, got: ", numBetaStr]; Exit[1]
  ];

  Print["Loading: ", algFile];
  Get[algFile];

  $alg    = ToExpression[algName];
  $energy = ToExpression[energyName];
  $sysParams = None;
  $algName   = algName;
  If[systemName == "", systemName = algName];

  If[!MatchQ[$alg, _Symbol],
    Print["ERROR: Algorithm '", algName, "' not found."]; Exit[1]
  ];
  If[!MatchQ[$energy, _Symbol],
    Print["ERROR: Energy function '", energyName, "' not found."]; Exit[1]
  ];
];

(* ---- Run the check ---- *)
Print[""];
Print[StringRepeat["=", 62]];
Print["System:     ", systemName];
Print["Seed state: ", $seedState];
Print["numBeta:    ", $numBeta];
Print["NSteps:     ", nSteps];
Print["MaxBitDepth:", maxBitDepth];
Print[StringRepeat["=", 62]];
Print[""];

result = RunFullCheck[
  $seedState, $alg, $energy, $numBeta,
  "SysParams"   -> $sysParams,
  "SystemName"  -> systemName,
  "NSteps"      -> nSteps,
  "MaxBitDepth" -> maxBitDepth,
  "Verbose"     -> verbose,
  "OpenWindow"  -> openWindow
];

Print[""];
Print[StringRepeat["=", 62]];
If[KeyExistsQ[result, "error"],
  Print["RESULT: ABORTED -- ", result["error"]],
  Print["RESULT: ", If[result["pass"], "PASS", "FAIL"],
        "  (KL divergence: ", result["kl"], ")"]
];
Print[StringRepeat["=", 62]];
