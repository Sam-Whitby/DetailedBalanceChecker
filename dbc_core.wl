(* ================================================================
   DetailedBalanceChecker - Core Library
   ================================================================

   Checks whether an MCMC algorithm satisfies detailed balance
   via exhaustive enumeration of all possible execution paths.

   ALGORITHM INTERFACE
   -------------------
   Provide TWO algorithm functions with the same signature:

     symAlg[state_, readBit_]   -- for the symbolic check.
                                   Energy differences are symbolic
                                   (e.g. eps2 - eps1).
                                   Use MetropolisProb[dE] for acceptance.

     numAlg[state_, readBit_]   -- for the numerical MCMC check.
                                   All quantities fully numeric.
                                   Use If[dE<=0, 1., Exp[-dE]] etc.,
                                   where dE already includes beta.

   Each function either:
     (a) Returns a single new state  (all randomness from readBit)
     (b) Returns {{p1,s1},{p2,s2},...} -- explicit probability-weighted
         outcomes.  Probabilities must sum to 1.
         Use this for Metropolis acceptance so it stays exact/symbolic.

   ENERGY INTERFACE
   ----------------
     symEnergy[state]  -- bare energy E(s), symbolic in coupling
                          constants.  Beta is NOT included; the library
                          inserts Exp[-beta * symEnergy[s]] for Boltzmann
                          weights and FullSimplify uses beta > 0.

     numEnergy[state]  -- fully numeric, WITH beta folded in, i.e.
                          numBeta * E_numeric(s).  Used to compute
                          Boltzmann weights as Exp[-numEnergy[s]].

   HELPER
   ------
     MetropolisProb[deltaE]  -- Piecewise Metropolis acceptance
                                probability, symbolic in beta and deltaE.

   ================================================================ *)

(* No package wrapper -- load with Get["dbc_core.wl"] *)

(* ----------------------------------------------------------------
   MetropolisProb
   Standard Metropolis probability as an exact symbolic Piecewise.
   deltaE is the bare energy difference (no beta factor).
   The global symbol \[Beta] (beta) is used; it stays symbolic
   during the symbolic check and is substituted by the caller in
   numerical contexts.
   ---------------------------------------------------------------- *)
MetropolisProb[deltaE_] :=
  Piecewise[{{1, deltaE <= 0}, {Exp[-\[Beta] * deltaE], deltaE > 0}}]

(* ----------------------------------------------------------------
   RunWithBits
   Run alg[state, readBit] with a fixed bit sequence.
   readBit[] returns successive bits; throws $dbc$outOfBits when
   the sequence is exhausted.
   Returns:
     {outcomes, nBitsConsumed}   outcomes = {{p,s},...}
     $OutOfBits                  if algorithm needed more bits
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
    (* Normalise to list-of-outcomes form {{p,s},...} *)
    {If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
       raw,
       {{1, raw}}],
     pos}
  ]
]

(* ----------------------------------------------------------------
   BuildTransitionMatrix
   BFS over bit sequences for every starting state.
   A sequence is extended by one bit whenever the algorithm calls
   readBit[] beyond its current length.  Each completed path at
   depth k contributes probability (1/2)^k * p_outcome.

   Returns Association  {fromState, toState} -> symbolicProbability.
   ---------------------------------------------------------------- *)
Options[BuildTransitionMatrix] = {
  "MaxBitDepth" -> 20,   (* hard cap on path length *)
  "TimeLimit"   -> 60.,  (* seconds budget per starting state *)
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
          (* extend tree by one bit level *)
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],

        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached for state ", s, " on bit prefix ", bits,
                " -- path excluded (algorithm may not halt)."],

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
   Verifies T(i->j) * pi(i) = T(j->i) * pi(j) for all pairs i<j,
   where pi(s) = Exp[-beta * symEnergy[s]].

   Uses FullSimplify with Assumptions -> {beta > 0} to handle the
   Piecewise expressions that arise from MetropolisProb.

   Returns a list of violations (empty = PASS).
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
   Run numAlg with true random bits and sample-from-outcomes logic
   for Metropolis branches.  All quantities should be fully numeric.
   Returns Association  state -> visitCount.
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

  (* Single MCMC step: call numAlg with fresh random bits,
     then sample from the returned outcome distribution. *)
  mcmcStep[] := Module[
    {liveRb, raw, outs, u, cumP, ns},
    liveRb[] := RandomInteger[1];
    raw  = numAlg[state, liveRb];
    outs = If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
               raw, {{1, raw}}];
    (* Weighted sample from outcome list *)
    u    = RandomReal[];
    cumP = 0.;
    ns   = outs[[-1, 2]];   (* fallback: last outcome *)
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
      (* state not in allStates -- count as error but don't crash *)
      Print["  WARNING: numerical MCMC produced unexpected state ", state]
    ],
    {nSteps - nWarmup}
  ];

  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeights
   Normalised Boltzmann probabilities.
   numEnergy[s] must be fully numeric and already include beta,
   so the weight is Exp[-numEnergy[s]].
   ---------------------------------------------------------------- *)
BoltzmannWeights[allStates_List, numEnergy_] := Module[
  {ws, Z},
  ws = N[Exp[-numEnergy[#]] & /@ allStates];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]

(* ----------------------------------------------------------------
   RunFullCheck
   Top-level entry point.  Prints a self-contained report.

   Arguments:
     allStates  -- list of all valid states
     symAlg     -- algorithm for symbolic check (uses MetropolisProb)
     numAlg     -- algorithm for numerical check (fully numeric)
     symEnergy  -- bare energy function, symbolic in couplings
     numEnergy  -- numeric energy function WITH beta folded in
   ---------------------------------------------------------------- *)
Options[RunFullCheck] = Join[
  Options[BuildTransitionMatrix],
  Options[RunNumericalMCMC],
  {"SystemName" -> "Unnamed system"}
]

RunFullCheck[allStates_List, symAlg_, numAlg_,
             symEnergy_, numEnergy_, OptionsPattern[]] := Module[
  {name = OptionValue["SystemName"],
   n    = Length[allStates],
   matrix, violations, counts, bw, simFreq, kl},

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
    "MaxBitDepth" -> OptionValue["MaxBitDepth"],
    "TimeLimit"   -> OptionValue["TimeLimit"],
    "Verbose"     -> OptionValue["Verbose"]
  ];

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
               PaddedForm[NumberForm[N@bw[s],    {5, 4}], 12]]
  ], allStates];

  (* KL divergence D(sim || Boltzmann) *)
  kl = Total@Table[
    With[{p = simFreq[s], q = N@bw[s]},
      If[p > 0 && q > 0, p * Log[p / q], 0.]],
    {s, allStates}];

  Print["\n  KL divergence (sim || Boltzmann) = ", NumberForm[kl, {6, 5}]];
  Print["  Numerical verdict: ",
    If[kl < 0.02,
       "CONSISTENT with Boltzmann.",
       "WARNING -- significant deviation from Boltzmann."]];

  sep["="];
  Print["END OF REPORT\n"];

  <|"matrix"     -> matrix,
    "violations" -> violations,
    "counts"     -> counts,
    "boltzmann"  -> bw|>
]
