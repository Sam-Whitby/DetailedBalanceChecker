(* ================================================================
   DetailedBalanceChecker  -  Core Library
   ================================================================
   Load with:  Get["path/to/dbc_core.wl"]

   Entry point:
     RunFullCheck[allStates, symAlg, numAlg, symEnergy, numEnergy, opts]

   See README.md for the full interface description.
   ================================================================ *)


(* ================================================================
   SECTION 1 – PRIMITIVES
   ================================================================ *)

(* ----------------------------------------------------------------
   MetropolisProb
   Standard Metropolis acceptance as a symbolic Piecewise.
   deltaE is the BARE energy difference (no beta); the global
   symbol \[Beta] appears in the result and stays unassigned during
   the symbolic check.
   ---------------------------------------------------------------- *)
MetropolisProb[deltaE_] :=
  Piecewise[{{1, deltaE <= 0}, {Exp[-\[Beta] deltaE], deltaE > 0}}]

(* ----------------------------------------------------------------
   RunWithBits
   Run alg[state, readBit] against a fixed bit list.
   readBit[] returns bits[[1]], bits[[2]], … and throws $OutOfBits
   if the list is exhausted.
   Returns  {outcomes, nBitsConsumed}  or  $OutOfBits.
   outcomes = {{p1,s1},{p2,s2},…}
   ---------------------------------------------------------------- *)
RunWithBits[alg_, state_, bits_List] := Module[
  {pos = 0, readBit, raw},
  readBit[] := (
    pos++;
    If[pos > Length[bits], Throw[$OutOfBits, $dbc$tag], bits[[pos]]]
  );
  raw = Catch[alg[state, readBit], $dbc$tag, ($OutOfBits &)];
  If[raw === $OutOfBits, $OutOfBits,
    {If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
       raw, {{1, raw}}], pos}
  ]
]


(* ================================================================
   SECTION 2 – TREE BUILDING
   ================================================================ *)

(* ----------------------------------------------------------------
   BuildTreeData
   BFS over bit sequences for every starting state.
   Returns  Association[ state -> { {bits, outcomes}, … } ]
   where outcomes = {{p1,s1},…}  (raw leaves of the decision tree).
   ---------------------------------------------------------------- *)
Options[BuildTreeData] = {
  "MaxBitDepth" -> 20,
  "TimeLimit"   -> 60.,
  "Verbose"     -> True
}

BuildTreeData[allStates_List, alg_, OptionsPattern[]] := Module[
  {maxDepth = OptionValue["MaxBitDepth"],
   tlim     = N @ OptionValue["TimeLimit"],
   verbose  = OptionValue["Verbose"],
   result   = <||>,
   queue, bits, res, outcomes, k, t0, timedOut, leaves},

  Do[
    If[verbose, Print["  Tree for state: ", s]];
    queue    = {{}};
    leaves   = {};
    t0       = AbsoluteTime[];
    timedOut = False;

    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[] - t0 > tlim, timedOut = True; Break[]];
      bits  = First[queue]; queue = Rest[queue];
      res   = RunWithBits[alg, s, bits];
      Which[
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],
        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached for state ", s, " at prefix ", bits,
                " -- path excluded (algorithm may not halt on this input)."],
        True,
          {outcomes, k} = res;
          AppendTo[leaves, {bits, outcomes}]
      ]
    ];

    If[timedOut, Print["  WARNING: Time limit reached for state ", s]];
    result[s] = leaves,
    {s, allStates}
  ];
  result
]

(* ----------------------------------------------------------------
   TreeDataToMatrix
   Derive the transition matrix from raw tree leaves.
   Returns  Association[ {from,to} -> symbolicProbability ]
   ---------------------------------------------------------------- *)
TreeDataToMatrix[allStates_List, treeData_Association] := Module[
  {stateSet = Association[# -> True & /@ allStates], matrix = <||>},
  Do[
    Do[
      With[{bits = leaf[[1]], outcomes = leaf[[2]]},
        Do[
          With[{p = out[[1]], ns = out[[2]]},
            If[KeyExistsQ[stateSet, ns],
              matrix[{s, ns}] =
                Lookup[matrix, Key[{s, ns}], 0] + p * (1/2)^Length[bits],
              Print["  WARNING: Invalid state ", ns, " returned from ", s]
            ]
          ],
          {out, outcomes}
        ]
      ],
      {leaf, treeData[s]}
    ],
    {s, allStates}
  ];
  matrix
]

(* BuildTransitionMatrix kept for API compatibility *)
Options[BuildTransitionMatrix] = Options[BuildTreeData]
BuildTransitionMatrix[allStates_List, alg_, opts : OptionsPattern[]] :=
  TreeDataToMatrix[allStates,
    BuildTreeData[allStates, alg,
      "MaxBitDepth" -> OptionValue["MaxBitDepth"],
      "TimeLimit"   -> OptionValue["TimeLimit"],
      "Verbose"     -> OptionValue["Verbose"]]]


(* ================================================================
   SECTION 3 – CHECKERS
   ================================================================ *)

(* ----------------------------------------------------------------
   CheckDetailedBalance
   Verifies T(i->j)*pi(i) = T(j->i)*pi(j) for all i<j pairs,
   where pi(s) = Exp[-beta * symEnergy[s]].
   Uses FullSimplify with beta>0 to resolve Piecewise expressions.
   Returns list of violation records; empty list = PASS.
   ---------------------------------------------------------------- *)
CheckDetailedBalance[matrix_Association, allStates_List, symEnergy_] := Module[
  {n = Length[allStates], violations = {}, si, sj, tij, tji, res},
  Do[
    si = allStates[[i]]; sj = allStates[[j]];
    tij = Lookup[matrix, Key[{si, sj}], 0];
    tji = Lookup[matrix, Key[{sj, si}], 0];
    res = FullSimplify[
      tij * Exp[-\[Beta] * symEnergy[si]] -
      tji * Exp[-\[Beta] * symEnergy[sj]],
      Assumptions -> {\[Beta] > 0}
    ];
    If[res =!= 0,
      AppendTo[violations, <|"pair" -> {si, sj}, "residual" -> res|>]
    ],
    {i, 1, n}, {j, i + 1, n}
  ];
  violations
]

(* ----------------------------------------------------------------
   RunNumericalMCMC
   Run numAlg with true random bits; sample outcomes weighted by
   their returned probabilities.
   Returns  Association[ state -> visitCount ]
   ---------------------------------------------------------------- *)
Options[RunNumericalMCMC] = {
  "NSteps"     -> 100000,
  "WarmupFrac" -> 0.1
}

RunNumericalMCMC[allStates_List, numAlg_, OptionsPattern[]] := Module[
  {nSteps  = OptionValue["NSteps"],
   nWarmup = Round[OptionValue["NSteps"] * OptionValue["WarmupFrac"]],
   state, counts},

  state  = RandomChoice[allStates];
  counts = AssociationThread[allStates -> 0];

  mcmcStep[] := Module[{liveRb, raw, outs, u, cumP, ns},
    liveRb[] := RandomInteger[1];
    raw  = numAlg[state, liveRb];
    outs = If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
               raw, {{1, raw}}];
    u = RandomReal[]; cumP = 0.; ns = outs[[-1, 2]];
    Do[cumP += N[out[[1]]]; If[u < cumP, ns = out[[2]]; Break[]], {out, outs}];
    state = ns
  ];

  Do[mcmcStep[], {nWarmup}];
  Do[mcmcStep[];
     If[KeyExistsQ[counts, state], counts[state]++,
        Print["  WARNING: unexpected state from numAlg: ", state]],
     {nSteps - nWarmup}];
  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeights
   numEnergy[s] must be fully numeric WITH beta already included.
   Returns  Association[ state -> normalisedWeight ]
   ---------------------------------------------------------------- *)
BoltzmannWeights[allStates_List, numEnergy_] := Module[
  {ws, Z},
  ws = N[Exp[-numEnergy[#]] & /@ allStates];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]


(* ================================================================
   SECTION 4 – VISUALISATION PRIMITIVES
   ================================================================ *)

(* ----------------------------------------------------------------
   DrawStateTree
   Render the BFS decision tree for one starting state as a Graph.
   leaves = { {bits_List, outcomes}, ... }
   outcomes = {{p1,s1},{p2,s2},...}
   Green outcome nodes  = algorithm moved to a NEW state.
   Orange outcome nodes = algorithm stayed in the SAME state.
   ---------------------------------------------------------------- *)

(* String vertex IDs -- safer than lists as Graph vertex names *)
$bKey[b_List] := If[b === {}, "root", "b" <> StringJoin[ToString /@ b]]
$oKey[b_List, i_Integer] := $bKey[b] <> "o" <> ToString[i]

(* Compact probability label *)
$probLabel[p_] := Which[
  p === 1 || p === 1.,   Style["p=1", 7, Darker[Green,0.3]],
  p === 0 || p === 0.,   Style["p=0", 7, Red],
  NumericQ[p],           Style["p=" <> ToString[NumberForm[N[p],{3,2}]], 7, GrayLevel[0.3]],
  True,                  Style[TraditionalForm[p], 7, GrayLevel[0.2]]
]

DrawStateTree[startState_, leaves_List] := Module[
  {leafBits, allPfx, intPfx, leafPfx,
   bitEdges, outEdges, allEdges,
   vLabels, eLabels, vStyle, vSize, nLeaves, maxDepth, g},

  If[leaves === {},
    Return @ Framed[Style["(no paths)", 9, Gray],
                   FrameStyle -> LightGray, ImageSize -> 140]
  ];

  leafBits = #[[1]] & /@ leaves;
  nLeaves  = Length[leafBits];
  maxDepth = Max[Length /@ leafBits];

  (* All bit-sequence prefixes appearing as nodes *)
  allPfx = DeleteDuplicates @ Flatten[
    Table[Take[b, k], {b, leafBits}, {k, 0, Length[b]}], 1];
  intPfx = DeleteDuplicates @ Flatten[
    Table[Take[b, k], {b, leafBits}, {k, 0, Length[b]-1}], 1];
  leafPfx = Complement[allPfx, intPfx];

  (* Edges between bit nodes *)
  bitEdges = DeleteDuplicates @ Flatten[
    Table[DirectedEdge[$bKey @ Take[b,k-1], $bKey @ Take[b,k]],
          {b, leafBits}, {k, 1, Length[b]}], 1];

  (* Edges to outcome nodes *)
  outEdges = Flatten @ Table[
    Table[DirectedEdge[$bKey[leaf[[1]]], $oKey[leaf[[1]],i]],
          {i, Length[leaf[[2]]]}],
    {leaf, leaves}];

  allEdges = Join[bitEdges, outEdges];

  (* Vertex labels *)
  vLabels = Flatten @ {
    $bKey[{}] -> Placed[Style["S=" <> ToString[startState], 9, Bold, White], Center],
    Table[$bKey[b] -> Placed[Style["?", 8, White], Center],
          {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> Placed[Style["\[DownArrow]", 9, GrayLevel[0.35]], Center],
          {b, leafPfx}],
    Flatten @ Table[
      With[{bits = leaf[[1]], outs = leaf[[2]]},
        Table[$oKey[bits,i] -> Placed[
          Column[{
            Style["\[RightArrow]" <> ToString[outs[[i,2]]], 8,
                  If[outs[[i,2]] =!= startState, Darker[Green,0.2], Darker[Orange,0.15]]],
            $probLabel[outs[[i,1]]]
          }, Alignment -> Center, Spacings -> 0.1], Center],
        {i, Length[outs]}]
      ],
      {leaf, leaves}]
  };

  (* Edge labels: 0 / 1 on bit edges *)
  eLabels = DeleteDuplicates @ Flatten[
    Table[DirectedEdge[$bKey @ Take[b,k-1], $bKey @ Take[b,k]] ->
            Placed[Style[ToString[b[[k]]], 9, Bold, RGBColor[0.2,0.3,0.7]], Automatic],
          {b, leafBits}, {k, 1, Length[b]}], 1];

  (* Vertex colours *)
  vStyle = Flatten @ {
    $bKey[{}] -> RGBColor[0.22,0.42,0.70],
    Table[$bKey[b] -> GrayLevel[0.52], {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> GrayLevel[0.70], {b, leafPfx}],
    Flatten @ Table[
      With[{bits = leaf[[1]], outs = leaf[[2]]},
        Table[$oKey[bits,i] ->
              If[outs[[i,2]] =!= startState,
                 RGBColor[0.25,0.68,0.38],
                 RGBColor[0.88,0.58,0.20]],
              {i, Length[outs]}]],
      {leaf, leaves}]
  };

  (* Vertex sizes *)
  vSize = Flatten @ {
    $bKey[{}] -> 0.65,
    Table[$bKey[b] -> 0.40, {b, Complement[intPfx, {{}}]}],
    Table[$bKey[b] -> 0.35, {b, leafPfx}],
    Flatten @ Table[Table[$oKey[leaf[[1]],i] -> 0.60, {i, Length[leaf[[2]]]}],
                    {leaf, leaves}]
  };

  Graph[
    allEdges,
    VertexLabels -> vLabels,
    EdgeLabels   -> eLabels,
    VertexStyle  -> vStyle,
    VertexSize   -> vSize,
    GraphLayout  -> {"LayeredDigraphEmbedding",
                     "RootVertex"   -> $bKey[{}],
                     "Orientation"  -> Top},
    ImageSize    -> {Max[220, 140*nLeaves], Max[180, 95*(maxDepth+2)]},
    Background   -> GrayLevel[0.97],
    PlotLabel    -> Style["State " <> ToString[startState], 10, Bold,
                          GrayLevel[0.3]]
  ]
]

(* ----------------------------------------------------------------
   MakeTransitionGrid
   Styled Grid showing the symbolic transition matrix.
   ---------------------------------------------------------------- *)
MakeTransitionGrid[allStates_List, matrix_Association] := Module[
  {n = Length[allStates], hdr, rows},
  hdr = Prepend[
    Style[ToString[#], Bold, 10, GrayLevel[0.2]] & /@ allStates,
    Style["T[i\[Rule]j]", 10, Italic, GrayLevel[0.4]]];
  rows = Table[
    Prepend[
      Table[
        With[{p = Lookup[matrix, Key[{allStates[[i]], allStates[[j]]}], 0]},
          If[p === 0,
             Style["0", 9, GrayLevel[0.7]],
             Style[TraditionalForm @ FullSimplify[p], 9]]],
        {j, n}],
      Style[ToString[allStates[[i]]], Bold, 10, GrayLevel[0.2]]],
    {i, n}];
  Grid[
    Prepend[rows, hdr],
    Frame      -> All,
    FrameStyle -> GrayLevel[0.82],
    Background -> {None, None, Flatten @ {
      Table[{1,j} -> RGBColor[0.84,0.90,1.00], {j,n+1}],
      Table[{i,1} -> RGBColor[0.84,0.90,1.00], {i,n+1}],
      Table[{i+1,i+1} -> RGBColor[1.00,0.98,0.84], {i,n}]}},
    Spacings   -> {2, 1},
    Alignment  -> Center
  ]
]

(* ----------------------------------------------------------------
   MakeDBTable
   Colour-coded table of detailed-balance pair results.
   ---------------------------------------------------------------- *)
MakeDBTable[allStates_List, violations_List] := Module[
  {n = Length[allStates], pairs, violPairs, hdr, rows},
  pairs     = Flatten[Table[{allStates[[i]],allStates[[j]]},
                            {i,n},{j,i+1,n}], 1];
  violPairs = If[violations === {}, {}, #["pair"] & /@ violations];
  hdr = Style[#, Bold, 10] & /@
        {"State i", "State j",
         "T(i\[Rule]j)\[CenterDot]\[Pi](i)  \[Minus]  T(j\[Rule]i)\[CenterDot]\[Pi](j)",
         "Result"};
  rows = Table[
    With[{pass = !MemberQ[violPairs, pair]},
      {Style[ToString[pair[[1]]], 10],
       Style[ToString[pair[[2]]], 10],
       Style[If[pass, "= 0   (FullSimplify, \[Beta] > 0)",
                      ToString[TraditionalForm @
                        First[Select[violations, #["pair"]===pair&],
                              <|"residual"->"?"|>]["residual"]]],
             9, If[pass, Darker[Green,0.2], Darker[Red,0.1]]],
       Style[If[pass, "\[Checkmark] PASS", "\[Times] FAIL"],
             10, Bold, If[pass, Darker[Green], Red]]}],
    {pair, pairs}];
  Grid[
    Prepend[rows, hdr],
    Frame      -> All,
    FrameStyle -> GrayLevel[0.82],
    Background -> {None, None,
      Table[{i+1,4} -> If[!MemberQ[violPairs,pairs[[i]]],
                           RGBColor[0.88,1.00,0.88],
                           RGBColor[1.00,0.88,0.88]],
            {i, Length[pairs]}]},
    Spacings   -> {2, 0.9},
    Alignment  -> {Left, Center}
  ]
]

(* ----------------------------------------------------------------
   MakeFrequencyPanel
   Bar chart of simulated vs Boltzmann frequencies + KL verdict.
   ---------------------------------------------------------------- *)
MakeFrequencyPanel[allStates_List, simFreq_Association,
                   bw_Association, kl_Real] := Module[
  {labels, simD, bwD, pass, chart, verdict},
  labels = ToString /@ allStates;
  simD   = N @ simFreq[#] & /@ allStates;
  bwD    = N @ bw[#]      & /@ allStates;
  pass   = kl < 0.02;
  chart  = BarChart[
    Table[{simD[[i]], bwD[[i]]}, {i, Length[allStates]}],
    ChartLabels  -> {Placed[labels, Below],
                     Placed[{"Simulated", "Boltzmann"}, Above]},
    ChartStyle   -> {Directive[RGBColor[0.20,0.45,0.80], Opacity[0.85]],
                     Directive[RGBColor[0.82,0.22,0.18], Opacity[0.85]]},
    ChartLegends -> Placed[{"Simulated", "Boltzmann"}, {Right, Top}],
    AxesLabel    -> {None, "Probability"},
    PlotLabel    -> Style["Simulated vs Boltzmann state frequencies", 12, Bold],
    ImageSize    -> {520, 300},
    Background   -> White];
  verdict = Style[
    "KL divergence (sim \[DoubleVerticalBar] Boltzmann) = " <>
    ToString[NumberForm[kl,{5,4}]] <> "     " <>
    If[pass, "\[Checkmark] Consistent with Boltzmann",
             "\[Times] Significant deviation from Boltzmann"],
    11, Bold, If[pass, Darker[Green], Red]];
  Column[{chart, Spacer[6], verdict}, Alignment -> Left, Spacings -> 0.3]
]


(* ================================================================
   SECTION 5 – REPORT WINDOW
   ================================================================ *)

(* ----------------------------------------------------------------
   MakeReportWindow
   Assembles all results into a Mathematica notebook and opens it.
   args = Association with keys:
     name, allStates, treeData, matrix, violations,
     simFreq, bw, kl, algCode
   ---------------------------------------------------------------- *)
MakeReportWindow[args_Association] := Module[
  {name, allStates, treeData, matrix, violations,
   simFreq, bw, kl, algCode, pass,
   badge, trees, matGrid, dbTab, freqPanel, cells, nb},

  name       = args["name"];
  allStates  = args["allStates"];
  treeData   = args["treeData"];
  matrix     = args["matrix"];
  violations = args["violations"];
  simFreq    = args["simFreq"];
  bw         = args["bw"];
  kl         = args["kl"];
  algCode    = args["algCode"];
  pass       = violations === {};

  (* ---- Pass/fail banner ---- *)
  badge = Panel[
    Row[{
      Style[name, 18, Bold, GrayLevel[0.1]],
      Spacer[25],
      Style[If[pass,
               "\[FilledCircle]  DETAILED BALANCE:  PASS",
               "\[FilledCircle]  DETAILED BALANCE:  FAIL"],
            17, Bold,
            If[pass, RGBColor[0.04,0.54,0.04], RGBColor[0.74,0.04,0.04]]]
    }],
    Background -> GrayLevel[0.94],
    FrameMargins -> {{16,16},{12,12}}];

  (* ---- Trees: one per state, side by side ---- *)
  trees = Table[
    Framed[
      Column[{DrawStateTree[s, treeData[s]]}, Alignment -> Center],
      FrameStyle -> GrayLevel[0.82], RoundingRadius -> 5,
      Background -> GrayLevel[0.975], FrameMargins -> 8],
    {s, allStates}];

  matGrid  = MakeTransitionGrid[allStates, matrix];
  dbTab    = MakeDBTable[allStates, violations];
  freqPanel = MakeFrequencyPanel[allStates, simFreq, bw, kl];

  (* ---- Notebook cells ---- *)
  cells = {
    Cell[BoxData @ ToBoxes @ badge,
         "Output", CellMargins -> {{8,8},{4,14}}],

    Cell["Algorithm Under Test", "Section"],
    Cell[algCode, "Code"],

    Cell["Decision Trees  (Symbolic Execution Paths)", "Section"],
    Cell[TextData[{
      "Each tree shows all bit sequences the algorithm can consume from a \
given starting state.  ",
      StyleBox["Blue", FontWeight->"Bold",
               FontColor->RGBColor[0.22,0.42,0.70]],
      " = root / internal bit-choice node.  ",
      StyleBox["Green", FontWeight->"Bold", FontColor->Darker[Green,0.2]],
      " outcome = algorithm moved to a new state.  ",
      StyleBox["Orange", FontWeight->"Bold", FontColor->Darker[Orange,0.1]],
      " outcome = algorithm stayed.  Edge labels are bit values (0 / 1).  \
Outcome labels show the destination state and symbolic acceptance probability."
    }], "Text"],
    Cell[BoxData @ ToBoxes @ Row[trees, Spacer[12]],
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Symbolic Transition Matrix  T[ i \[Rule] j ]", "Section"],
    Cell["Entry (i, j) = total probability of transitioning FROM state i \
TO state j, accumulated over all execution paths.  \[Beta] is kept symbolic \
throughout.", "Text"],
    Cell[BoxData @ ToBoxes @ matGrid,
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Detailed Balance Check:  \
T(i\[Rule]j)\[CenterDot]\[Pi](i) = T(j\[Rule]i)\[CenterDot]\[Pi](j)", "Section"],
    Cell["Every pair of distinct states is tested using FullSimplify with \
\[Beta] > 0 so that Piecewise Metropolis expressions are resolved exactly.  \
\[Pi](s) = Exp[\[Minus]\[Beta] E(s)].", "Text"],
    Cell[BoxData @ ToBoxes @ dbTab,
         "Output", CellMargins -> {{8,8},{4,4}}],

    Cell["Numerical MCMC Validation", "Section"],
    Cell["The algorithm is run as a live Markov chain with genuinely random \
bits.  Simulated state frequencies are compared to the analytical Boltzmann \
distribution Exp[\[Minus]\[Beta]E(s)] / Z.", "Text"],
    Cell[BoxData @ ToBoxes @ freqPanel,
         "Output", CellMargins -> {{8,8},{4,18}}]
  };

  (* ---- Open notebook window ---- *)
  nb = Quiet @ Check[
    CreateDocument[
      cells,
      WindowTitle  -> "DetailedBalanceChecker \[LongDash] " <> name,
      WindowSize   -> {1100, 860},
      WindowMargins -> {{Automatic, Automatic}, {Automatic, 0}},
      Editable     -> False,
      Background   -> White,
      StyleDefinitions -> "Default.nb"
    ],
    (Print["  (No Mathematica frontend available; window not opened.)"]; None)
  ];
  nb
]


(* ================================================================
   SECTION 6 – TOP-LEVEL ENTRY POINT
   ================================================================ *)

(* ----------------------------------------------------------------
   RunFullCheck
   Orchestrates BFS, symbolic check, numerical MCMC, and the
   graphical report window.

   Arguments:
     allStates   list of all valid system states
     symAlg      algorithm for symbolic check  (uses MetropolisProb)
     numAlg      algorithm for numerical MCMC  (fully numeric)
     symEnergy   bare energy, symbolic in couplings, no beta
     numEnergy   numeric energy WITH beta folded in

   Options:
     "SystemName"       display name
     "AlgorithmCode"    string shown in the Algorithm section
                        (Automatic = extracted from symAlg DownValues)
     "MaxBitDepth"      BFS depth cap per state (default 20)
     "TimeLimit"        seconds per state BFS (default 60)
     "Verbose"          print BFS progress (default True)
     "NSteps"           MCMC steps (default 100 000)
     "WarmupFrac"       warm-up fraction (default 0.1)
     "OpenWindow"       open graphical report window (default True)
   ---------------------------------------------------------------- *)
Options[RunFullCheck] = Join[
  Options[BuildTreeData],
  Options[RunNumericalMCMC],
  {"SystemName"    -> "Unnamed system",
   "AlgorithmCode" -> Automatic,
   "OpenWindow"    -> True}
]

RunFullCheck[allStates_List, symAlg_, numAlg_,
             symEnergy_, numEnergy_, OptionsPattern[]] := Module[
  {name     = OptionValue["SystemName"],
   n        = Length[allStates],
   algCode, treeData, matrix, violations,
   counts, bw, simFreq, kl, pass},

  algCode = OptionValue["AlgorithmCode"];
  If[algCode === Automatic,
    algCode = StringTrim @
              ToString[InputForm[DownValues[symAlg]], OutputForm]];

  (* ---- Terminal progress ---- *)
  Print[StringRepeat["=", 62]];
  Print["DETAILED BALANCE CHECKER"];
  Print[StringRepeat["=", 62]];
  Print["System : ", name];
  Print["States : ", n, "  --  ", allStates];
  Print[StringRepeat["=", 62]];

  (* ---- 1. Symbolic tree building ---- *)
  Print["\n[1/3] Building decision trees (symAlg) ..."];
  treeData = BuildTreeData[allStates, symAlg,
               "MaxBitDepth" -> OptionValue["MaxBitDepth"],
               "TimeLimit"   -> OptionValue["TimeLimit"],
               "Verbose"     -> OptionValue["Verbose"]];
  matrix   = TreeDataToMatrix[allStates, treeData];
  Print["      Transition matrix: ", Length[matrix], " non-zero entries."];

  (* ---- 2. Detailed balance check ---- *)
  Print["\n[2/3] Checking detailed balance for ",
        Binomial[n,2], " pairs ..."];
  violations = CheckDetailedBalance[matrix, allStates, symEnergy];
  pass       = violations === {};
  If[pass,
    Print["      RESULT: PASS -- all pairs satisfy detailed balance exactly."],
    Print["      RESULT: FAIL -- ", Length[violations], " violation(s) found."]
  ];

  (* ---- 3. Numerical MCMC ---- *)
  Print["\n[3/3] Running ", OptionValue["NSteps"],
        " MCMC steps (numAlg) ..."];
  counts  = RunNumericalMCMC[allStates, numAlg,
              "NSteps"     -> OptionValue["NSteps"],
              "WarmupFrac" -> OptionValue["WarmupFrac"]];
  bw      = BoltzmannWeights[allStates, numEnergy];
  simFreq = N[# / Total[counts]] & /@ counts;
  kl      = Total @ Table[
    With[{p = simFreq[s], q = N@bw[s]},
      If[p > 0 && q > 0, p * Log[p/q], 0.]], {s, allStates}];

  Print["      KL divergence (sim || Boltzmann) = ",
        NumberForm[kl,{5,4}]];
  Print["      Numerical: ",
    If[kl < 0.02, "CONSISTENT with Boltzmann.",
                  "WARNING -- significant deviation from Boltzmann."]];

  Print["\n", StringRepeat["=", 62]];
  Print["OVERALL: ", If[pass, "PASS", "FAIL"]];
  Print[StringRepeat["=", 62]];

  (* ---- Open graphical window ---- *)
  If[TrueQ @ OptionValue["OpenWindow"],
    Print["\nOpening report window ..."];
    MakeReportWindow[<|
      "name"       -> name,
      "allStates"  -> allStates,
      "treeData"   -> treeData,
      "matrix"     -> matrix,
      "violations" -> violations,
      "simFreq"    -> simFreq,
      "bw"         -> bw,
      "kl"         -> kl,
      "algCode"    -> algCode
    |>]
  ];

  <|"matrix"     -> matrix,
    "violations" -> violations,
    "counts"     -> counts,
    "boltzmann"  -> bw,
    "kl"         -> kl,
    "pass"       -> pass|>
]
