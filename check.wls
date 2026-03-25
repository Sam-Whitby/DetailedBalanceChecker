#!/usr/bin/env wolframscript
(* ================================================================
   check.wls  --  Exhaustive detailed-balance checker
   ================================================================
   Usage:
     wolframscript -file check.wls <alg_file> [options]

   The .wl file must define:
     energy[state_]      bare energy function (no beta)
     Algorithm[state_]   MCMC move function
     BitsToState[bits_]  {0,1,...} list -> system state, or None if
                         the bit pattern does not correspond to a
                         valid state for this system
     numBeta             numeric inverse temperature
   Optionally:
     symParams = <|"eps" -> {...}, "couplings" -> {...}|>
     (keeps energy parameters symbolic during the symbolic check;
      assigns random numeric values for the MCMC run)

   How it works:
     The checker generates all bit strings of length 1 through
     MaxBits, in order of increasing length and numerically within
     each length.  For each bit string it calls BitsToState to
     obtain a seed state.  If BitsToState returns None the string
     is skipped.  Otherwise BFS discovers the full connected
     component reachable from that seed.  If the component has
     already been seen (deduplication), the string is skipped.
     Otherwise the checker runs the requested symbolic and/or
     numerical tests and prints one result row.

   Options (key=value, no spaces around =):
     MaxBits=N           Largest bit string length to test (default 8)
     Mode=Symbolic       Symbolic detailed-balance check only
     Mode=Numerical      Numerical MCMC check only
     Mode=Both           Both checks (default)
     NSteps=N            MCMC steps per seed state (default 50000)
     MaxBitDepth=N       BFS tree depth cap per state (default 20)
     Verbose=True/False  Print per-state BFS progress (default False)
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "check.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1,
  Print["Usage: wolframscript -file check.wls <alg_file> [options]"];
  Print[""];
  Print["The .wl file must define:"];
  Print["  energy[state_]      bare energy (no beta)"];
  Print["  Algorithm[state_]   MCMC move function"];
  Print["  BitsToState[bits_]  {0,1,...} -> state, or None if invalid"];
  Print["  numBeta             numeric inverse temperature"];
  Print["  symParams           (optional) symbolic parameter association"];
  Print[""];
  Print["Options (key=value, no spaces around =):"];
  Print["  MaxBits=N           Bit string lengths 1..N (default 8)"];
  Print["  Mode=Symbolic|Numerical|Both  (default Both)"];
  Print["  NSteps=N            MCMC steps per system (default 50000)"];
  Print["  MaxBitDepth=N       BFS depth cap (default 20)"];
  Print["  Verbose=True/False  Per-state BFS output (default False)"];
  Exit[1]];

(* ---- Option defaults ---- *)
maxBits     = 8;
mode        = "Both";
nSteps      = 50000;
maxBitDepth = 20;
verbose     = False;

Do[
  If[StringMatchQ[opt, "MaxBits="     ~~ __], maxBits     = ToExpression[StringDrop[opt,  8]]];
  If[StringMatchQ[opt, "Mode="        ~~ __], mode        =              StringDrop[opt,  5]];
  If[StringMatchQ[opt, "NSteps="      ~~ __], nSteps      = ToExpression[StringDrop[opt,  7]]];
  If[StringMatchQ[opt, "MaxBitDepth=" ~~ __], maxBitDepth = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "Verbose="     ~~ __], verbose     = ToExpression[StringDrop[opt,  8]]],
  {opt, kvArgs}];

(* ---- Resolve file path ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

(* ---- Load algorithm file ---- *)
ClearAll[Algorithm, energy, BitsToState, numBeta, symParams];
Print["Loading: ", algFile];
Get[algFile];

(* ---- Validate required symbols ---- *)
If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[energy]] == 0,
  Print["ERROR: 'energy' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[BitsToState]] == 0,
  Print["ERROR: 'BitsToState' not defined in ", posArgs[[1]]]; Exit[1]];
If[!ValueQ[numBeta] || !NumericQ[numBeta],
  Print["ERROR: 'numBeta' (numeric) not defined in ", posArgs[[1]]]; Exit[1]];

$alg       = Algorithm;
$energy    = energy;
$numBeta   = numBeta;
$sysParams = If[ValueQ[symParams], symParams, None];

(* ---- Symbolic energy parameters ---- *)
If[AssociationQ[$sysParams],
  $symVars        = Join[Lookup[$sysParams, "eps",       {}],
                         Lookup[$sysParams, "couplings", {}]];
  $symAssumptions = Map[Element[#, Reals] &, $symVars],
  (* else: energy is purely numeric *)
  $symVars        = {};
  $symAssumptions = {}];

(* ---- Generate all bit strings of length 1..maxBits ---- *)
(* Ordered by length first, then numerically within each length. *)
allBitStrings = Flatten[
  Table[IntegerDigits[i, 2, len],
        {len, 1, maxBits},
        {i, 0, 2^len - 1}],
  1];

(* ---- Column widths ---- *)
colBits  = 14;
colState = 26;
colN     = 9;
colSym   = 10;
colKL    = 12;
colNum   = 10;

totalWidth = colBits + colState + colN +
  If[mode =!= "Numerical", colSym, 0] +
  If[mode =!= "Symbolic",  colKL + colNum, 0];

(* ---- Print header ---- *)
Print[];
Print[StringRepeat["=", 75]];
Print["DETAILED BALANCE CHECKER -- BATCH MODE"];
Print["File:        ", algFile];
Print["MaxBits:     ", maxBits,
      "  (", Length[allBitStrings], " bit strings generated)"];
Print["Mode:        ", mode];
If[mode =!= "Symbolic", Print["NSteps:      ", nSteps]];
Print["MaxBitDepth: ", maxBitDepth];
Print[StringRepeat["=", 75]];
Print[];

Print[StringPadRight["Bits",    colBits],
      StringPadRight["State",   colState],
      StringPadRight["#States", colN],
      If[mode =!= "Numerical", StringPadRight["Symbolic",  colSym], ""],
      If[mode =!= "Symbolic",
         StringPadRight["KL",       colKL] <>
         StringPadRight["Numerical",colNum], ""]];
Print[StringRepeat["-", 75]];

(* ================================================================
   Main loop
   ================================================================ *)
coveredStates = {};
testedCount   = 0;
symPassCount  = 0;
symFailCount  = 0;
numPassCount  = 0;
numFailCount  = 0;

Do[
  (* Convert bit string to system state via user-supplied bijection *)
  state = Check[BitsToState[bitString], $Failed];
  If[state === $Failed || state === None, Continue[]];

  (* Skip if this seed is already in a known connected component *)
  If[MemberQ[coveredStates, state], Continue[]];

  bitsStr  = StringJoin[ToString /@ bitString];
  stateStr = ToString[state];

  (* ---- BFS: discover connected component from seed ---- *)
  treeData = BuildTreeAT[state, $alg,
               "MaxBitDepth" -> maxBitDepth,
               "TimeLimit"   -> 120.,
               "Verbose"     -> verbose];

  (* Algorithm contains an interceptable call -- print error and move on *)
  If[MatchQ[treeData, $dbc$cantHandle[_]],
    Print[StringPadRight[bitsStr,  colBits],
          StringPadRight[stateStr, colState],
          StringPadRight["?",      colN],
          "ERROR: ", treeData[[1]]];
    testedCount++;
    Continue[]];

  allStates = Sort @ Keys[treeData];
  matrix    = TreeATToMatrix[treeData];
  nStates   = Length[allStates];

  (* Record all discovered states so later seeds in this component
     are skipped *)
  coveredStates = Union[coveredStates, allStates];
  testedCount++;

  (* ---- Symbolic detailed-balance check ---- *)
  symResult = "-";
  If[mode =!= "Numerical",
    violations = CheckDetailedBalance[
                   matrix, allStates, $energy, $symAssumptions];
    symResult  = If[violations === {}, "PASS", "FAIL"];
    If[symResult === "PASS", symPassCount++, symFailCount++]
  ];

  (* ---- Numerical MCMC check ---- *)
  (* Skip for trivial single-state systems *)
  klStr     = "-";
  numResult = "-";
  If[mode =!= "Symbolic" && nStates > 1,
    If[$symVars =!= {},
      (* Assign random numeric values to symbolic energy parameters *)
      numSubs = MakeNumericSubs[$sysParams, 42];
      numVals = Last /@ numSubs;
      Block[Evaluate[$symVars],
        MapThread[Set, {$symVars, numVals}];
        mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta,
                       "NSteps" -> nSteps];
        boltzmann  = BoltzmannWeightsAT[allStates, $energy, $numBeta]
      ],
      (* Fully numeric energy: proceed directly *)
      mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta,
                     "NSteps" -> nSteps];
      boltzmann  = BoltzmannWeightsAT[allStates, $energy, $numBeta]
    ];
    simFreq = N[# / Total[mcmcCounts]] & /@ mcmcCounts;
    kl = Total @ Table[
      With[{p = simFreq[s], q = N @ boltzmann[s]},
        If[p > 0 && q > 0, p * Log[p / q], 0.]], {s, allStates}];
    klStr     = ToString[NumberForm[kl, {4, 3}]];
    numResult = If[kl < 0.02, "PASS", "FAIL"];
    If[numResult === "PASS", numPassCount++, numFailCount++]
  ];

  (* ---- Print result row ---- *)
  symCol = If[mode =!= "Numerical", StringPadRight[symResult, colSym], ""];
  numCol = If[mode =!= "Symbolic",
              StringPadRight[klStr,     colKL] <>
              StringPadRight[numResult, colNum], ""];
  Print[StringPadRight[bitsStr,             colBits],
        StringPadRight[stateStr,            colState],
        StringPadRight[ToString[nStates],   colN],
        symCol,
        numCol],

  {bitString, allBitStrings}
];

(* ================================================================
   Summary
   ================================================================ *)
Print[];
Print[StringRepeat["-", 75]];
Print["Systems tested (distinct connected components): ", testedCount];
If[mode =!= "Numerical",
  Print["Symbolic  -- PASS: ", symPassCount, "  FAIL: ", symFailCount]];
If[mode =!= "Symbolic" && (numPassCount + numFailCount) > 0,
  Print["Numerical -- PASS: ", numPassCount, "  FAIL: ", numFailCount]];
If[symFailCount == 0 && numFailCount == 0 && testedCount > 0,
  Print["All tested systems satisfy detailed balance."],
  If[testedCount == 0,
    Print["No valid states found for bit strings of length 1..", maxBits, "."],
    Print["*** DETAILED BALANCE VIOLATION DETECTED ***"]
  ]
];
Print[StringRepeat["=", 75]];
