(* ================================================================
   DetailedBalanceChecker - Core Library
   ================================================================

   Checks whether an MCMC algorithm satisfies detailed balance
   via exhaustive enumeration of all possible execution paths.

   ALGORITHM INTERFACE
   -------------------
   Provide TWO algorithm functions with the same signature:

     symAlg[state_, readBit_]   -- for the symbolic check.
                                   Energy differences are symbolic.
                                   Use MetropolisProb[dE] for acceptance.

     numAlg[state_, readBit_]   -- for the numerical MCMC check.
                                   All quantities fully numeric.
                                   dE must already include beta.

   Each returns EITHER:
     (a) A single new state     (all randomness from readBit)
     (b) {{p1,s1},{p2,s2},...}  explicit probability-weighted outcomes

   ENERGY INTERFACE
   ----------------
     symEnergy[state]  bare energy E(s), symbolic in couplings, NO beta.
                       The library inserts beta via Exp[-beta*symEnergy[s]].

     numEnergy[state]  fully numeric, WITH beta folded in.
                       Boltzmann weights computed as Exp[-numEnergy[s]].

   HELPER
   ------
     MetropolisProb[deltaE]  standard Metropolis acceptance probability,
                             symbolic Piecewise in beta and deltaE.

   ================================================================ *)

(* ----------------------------------------------------------------
   MetropolisProb
   ---------------------------------------------------------------- *)
MetropolisProb[deltaE_] :=
  Piecewise[{{1, deltaE <= 0}, {Exp[-\[Beta] * deltaE], deltaE > 0}}]

(* ----------------------------------------------------------------
   RunWithBits
   Run alg[state, readBit] with a fixed bit sequence.
   Returns {outcomes, nBitsConsumed} or $OutOfBits.
   ---------------------------------------------------------------- *)
RunWithBits[alg_, state_, bits_List] := Module[
  {pos = 0, readBit, raw},
  readBit[] := (
    pos++;
    If[pos > Length[bits],
      Throw[$OutOfBits, $dbc$outOfBits],
      bits[[pos]]
    ]
  );
  raw = Catch[alg[state, readBit], $dbc$outOfBits, ($OutOfBits &)];
  If[raw === $OutOfBits,
    $OutOfBits,
    {If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
       raw, {{1, raw}}],
     pos}
  ]
]

(* ----------------------------------------------------------------
   BuildTransitionMatrix
   BFS over bit sequences for every starting state.
   Returns Association {fromState, toState} -> symbolicProbability.
   ---------------------------------------------------------------- *)
Options[BuildTransitionMatrix] = {
  "MaxBitDepth" -> 20,
  "TimeLimit"   -> 60.,
  "Verbose"     -> True
}

BuildTransitionMatrix[allStates_List, alg_, OptionsPattern[]] := Module[
  {maxDepth = OptionValue["MaxBitDepth"],
   tlim     = OptionValue["TimeLimit"],
   verbose  = OptionValue["Verbose"],
   stateSet, matrix, queue, bits, res, outcomes, k, t0, timedOut},

  stateSet = Association[# -> True & /@ allStates];
  matrix   = <||>;

  Do[
    If[verbose, Print["  Tree for state: ", s]];
    queue    = {{}};
    t0       = AbsoluteTime[];
    timedOut = False;

    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[] - t0 > tlim, timedOut = True; Break[]];
      bits  = First[queue];
      queue = Rest[queue];
      res   = RunWithBits[alg, s, bits];

      Which[
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],

        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached for state ", s, " on prefix ", bits,
                " -- path excluded (algorithm may not halt on this input)."],

        True,
          {outcomes, k} = res;
          Do[
            With[{p = out[[1]], ns = out[[2]]},
              If[!KeyExistsQ[stateSet, ns],
                Print["  WARNING: Algorithm returned invalid state ", ns,
                      " from ", s],
                matrix[{s, ns}] =
                  Lookup[matrix, Key[{s, ns}], 0] + p * (1/2)^k
              ]
            ],
            {out, outcomes}
          ]
      ]
    ];

    If[timedOut, Print["  WARNING: Time limit reached for state ", s]],
    {s, allStates}
  ];

  matrix
]

(* ----------------------------------------------------------------
   CheckDetailedBalance
   Verifies T(i->j)*pi(i) = T(j->i)*pi(j) for all i<j pairs,
   where pi(s) = Exp[-beta * symEnergy[s]].
   Uses FullSimplify with beta > 0 to resolve Piecewise expressions.
   Returns list of violations; empty = PASS.
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
   Run numAlg with true random bits.
   Returns Association state -> visitCount.
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

  mcmcStep[] := Module[
    {liveRb, raw, outs, u, cumP, ns},
    liveRb[] := RandomInteger[1];
    raw  = numAlg[state, liveRb];
    outs = If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
               raw, {{1, raw}}];
    u    = RandomReal[];
    cumP = 0.;
    ns   = outs[[-1, 2]];
    Do[
      cumP += N[out[[1]]];
      If[u < cumP, ns = out[[2]]; Break[]],
      {out, outs}
    ];
    state = ns
  ];

  Do[mcmcStep[], {nWarmup}];
  Do[
    mcmcStep[];
    If[KeyExistsQ[counts, state],
      counts[state]++,
      Print["  WARNING: numerical MCMC produced unexpected state ", state]
    ],
    {nSteps - nWarmup}
  ];

  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeights
   Normalised Boltzmann probabilities from numEnergy (includes beta).
   ---------------------------------------------------------------- *)
BoltzmannWeights[allStates_List, numEnergy_] := Module[
  {ws, Z},
  ws = N[Exp[-numEnergy[#]] & /@ allStates];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]

(* ----------------------------------------------------------------
   NumericalTransitionMatrix
   Build an n x n numeric matrix from numAlg via BFS.
   Entry [i,j] = probability of going from allStates[[i]] to allStates[[j]].
   ---------------------------------------------------------------- *)
NumericalTransitionMatrix[allStates_List, numAlg_, opts___] := Module[
  {n = Length[allStates], mat, matrix},
  matrix = BuildTransitionMatrix[allStates, numAlg,
             "Verbose" -> False, opts];
  Table[
    N @ Lookup[matrix, Key[{allStates[[i]], allStates[[j]]}], 0],
    {i, n}, {j, n}]
]

(* ----------------------------------------------------------------
   ExportPlots
   Generate and export two PNG files for a given run:
     1. Frequency bar chart   (simulated vs Boltzmann per state)
     2. Transition matrix heatmap (numerical, from numAlg)

   Arguments:
     allStates  list of states
     numAlg     numerical algorithm (for heatmap)
     simFreq    Association state -> simulated frequency
     bw         Association state -> Boltzmann probability
     name       system name string (used in titles and file names)
     outDir     output directory string (include trailing slash)
   ---------------------------------------------------------------- *)
ExportPlots[allStates_List, numAlg_, simFreq_, bw_,
            name_String, outDir_String, algOpts___] := Module[
  {n, labels, safeName, freqData, bwData, freqPlot, numMat, ticks, matPlot, f1, f2},

  n       = Length[allStates];
  labels  = ToString /@ allStates;
  safeName = StringReplace[name, {" " -> "_", "/" -> "-",
                                   "\\" -> "-", "," -> ""}];

  (* ---- 1. Frequency comparison bar chart ---- *)
  freqData = N @ simFreq[#] & /@ allStates;
  bwData   = N @ bw[#]      & /@ allStates;

  freqPlot = BarChart[
    (* one group per state, two bars: simulated then Boltzmann *)
    Table[{freqData[[i]], bwData[[i]]}, {i, n}],
    ChartLabels -> {
      Placed[labels, Below],
      Placed[{"Sim", "Boltz"}, Above]
    },
    ChartStyle  -> {Directive[RGBColor[0.2, 0.4, 0.8], Opacity[0.85]],
                    Directive[RGBColor[0.8, 0.2, 0.2], Opacity[0.85]]},
    ChartLegends -> Placed[{"Simulated", "Boltzmann"}, {Right, Top}],
    AxesLabel    -> {None, "Probability"},
    PlotLabel    -> Style[name <> "\nSimulated vs Boltzmann frequencies", 13, Bold],
    ImageSize    -> {600, 400},
    Background   -> White
  ];

  (* ---- 2. Numerical transition matrix heatmap ---- *)
  numMat = NumericalTransitionMatrix[allStates, numAlg, algOpts];
  ticks  = Table[{i, Rotate[labels[[i]], 0]}, {i, n}];

  matPlot = MatrixPlot[
    numMat,
    ColorFunction  -> "SunsetColors",
    ColorFunctionScaling -> True,
    FrameTicks     -> {{ticks, None}, {ticks, None}},
    FrameLabel     -> {Style["From state", 11], Style["To state", 11]},
    PlotLabel      -> Style[name <> "\nNumerical transition matrix  T[from, to]", 13, Bold],
    ImageSize      -> {500, 500},
    Background     -> White,
    (* Overlay numeric values for small matrices *)
    Epilog -> If[n <= 8,
      Flatten @ Table[
        Text[Style[NumberForm[numMat[[i, j]], {3, 2}], 9, GrayLevel[0.1]],
             {j, n + 1 - i}],   (* MatrixPlot: x=col, y=n+1-row *)
        {i, n}, {j, n}],
      {}]
  ];

  (* ---- Export ---- *)
  f1 = outDir <> safeName <> "_frequencies.png";
  f2 = outDir <> safeName <> "_transition_matrix.png";
  Export[f1, freqPlot,  "PNG"];
  Export[f2, matPlot, "PNG"];

  Print["  Plots exported:"];
  Print["    Frequency chart  -> ", f1];
  Print["    Transition matrix -> ", f2];

  <|"frequencyPlot" -> freqPlot, "matrixPlot" -> matPlot,
    "freqFile" -> f1, "matrixFile" -> f2|>
]

(* ----------------------------------------------------------------
   RunFullCheck
   Top-level entry point.  Runs symbolic + numerical checks and
   prints a self-contained report.  Optionally exports plots.

   Arguments:
     allStates  list of all valid states
     symAlg     algorithm for symbolic check (uses MetropolisProb)
     numAlg     algorithm for numerical check (fully numeric)
     symEnergy  bare energy function, symbolic in couplings
     numEnergy  numeric energy function WITH beta folded in
   ---------------------------------------------------------------- *)
Options[RunFullCheck] = Join[
  Options[BuildTransitionMatrix],
  Options[RunNumericalMCMC],
  {"SystemName"   -> "Unnamed system",
   "ExportPlots"  -> False,
   "PlotDirectory" -> "."}
]

RunFullCheck[allStates_List, symAlg_, numAlg_,
             symEnergy_, numEnergy_, OptionsPattern[]] := Module[
  {name     = OptionValue["SystemName"],
   doPlots  = OptionValue["ExportPlots"],
   plotDir  = OptionValue["PlotDirectory"],
   n        = Length[allStates],
   algOpts, matrix, violations, counts, bw, simFreq, kl, plotDir2},

  (* Options to forward to BuildTransitionMatrix *)
  algOpts = Sequence[
    "MaxBitDepth" -> OptionValue["MaxBitDepth"],
    "TimeLimit"   -> OptionValue["TimeLimit"]
  ];

  (* Ensure plot directory ends with a separator *)
  plotDir2 = If[StringEndsQ[plotDir, $PathnameSeparator],
                plotDir,
                plotDir <> $PathnameSeparator];

  sep[c_] := Print[StringRepeat[c, 62]];

  sep["="];
  Print["DETAILED BALANCE CHECKER"];
  sep["="];
  Print["System : ", name];
  Print["States : ", n, " -- ", allStates];
  sep["="];

  (* ---- 1. Symbolic check ---- *)
  Print["\n[1/2] SYMBOLIC DETAILED BALANCE CHECK"];
  sep["-"];
  Print["Building decision trees for all ", n, " states..."];

  matrix = BuildTransitionMatrix[allStates, symAlg,
    algOpts, "Verbose" -> OptionValue["Verbose"]];

  Print["Transition matrix: ", Length[matrix], " non-zero entries."];
  Print["Checking ", Binomial[n, 2], " pairs for detailed balance..."];

  violations = CheckDetailedBalance[matrix, allStates, symEnergy];

  If[violations === {},
    Print["\n  RESULT: PASS -- all ", Binomial[n, 2],
          " pairs satisfy detailed balance exactly."],
    Print["\n  RESULT: FAIL -- ", Length[violations], " violation(s):"];
    Scan[Function[v,
      Print["    ", v["pair"], "  residual = ", v["residual"]]
    ], violations]
  ];

  (* ---- 2. Numerical check ---- *)
  Print["\n[2/2] NUMERICAL MCMC CHECK"];
  sep["-"];
  Print["Running ", OptionValue["NSteps"], " MCMC steps..."];

  counts  = RunNumericalMCMC[allStates, numAlg,
              "NSteps"     -> OptionValue["NSteps"],
              "WarmupFrac" -> OptionValue["WarmupFrac"]];
  bw      = BoltzmannWeights[allStates, numEnergy];
  simFreq = N[# / Total[counts]] & /@ counts;

  Print["\n  ", PaddedForm["State", 20],
               PaddedForm["Simulated", 12],
               PaddedForm["Boltzmann", 12]];
  Print["  ", StringRepeat["-", 44]];
  Scan[Function[s,
    Print["  ", PaddedForm[ToString[s], 20],
               PaddedForm[NumberForm[simFreq[s], {5, 4}], 12],
               PaddedForm[NumberForm[N @ bw[s],  {5, 4}], 12]]
  ], allStates];

  kl = Total @ Table[
    With[{p = simFreq[s], q = N @ bw[s]},
      If[p > 0 && q > 0, p * Log[p / q], 0.]],
    {s, allStates}];

  Print["\n  KL divergence (sim || Boltzmann) = ", NumberForm[kl, {6, 5}]];
  Print["  Numerical verdict: ",
    If[kl < 0.02,
       "CONSISTENT with Boltzmann.",
       "WARNING -- significant deviation from Boltzmann."]];

  (* ---- 3. Plots (optional) ---- *)
  If[doPlots,
    Print["\n[+] Exporting plots..."];
    sep["-"];
    ExportPlots[allStates, numAlg, simFreq, bw, name, plotDir2, algOpts]
  ];

  sep["="];
  Print["END OF REPORT\n"];

  <|"matrix"     -> matrix,
    "violations" -> violations,
    "counts"     -> counts,
    "boltzmann"  -> bw|>
]
