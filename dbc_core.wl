(* ================================================================
   DetailedBalanceChecker - Core Library
   ================================================================

   Checks whether an MCMC algorithm satisfies detailed balance
   via exhaustive enumeration of all possible execution paths.

   ALGORITHM INTERFACE
   -------------------
   Your algorithm must have signature:

     myAlg[state_, readBit_] := ...

   where readBit[] returns 0 or 1 each time it is called (discrete
   random choices, e.g. site selection).

   The function returns EITHER:
     (a) A single new state  (all randomness consumed via readBit)
     (b) A list {{p1,s1}, {p2,s2}, ...} of symbolic {probability, newState}
         pairs summing to 1.  Use this for continuous acceptance steps
         (Metropolis) so the probability is kept exact and symbolic.

   The helper MetropolisProb[deltaE] produces the standard Metropolis
   acceptance probability as a symbolic Piecewise expression.

   ENERGY FUNCTION INTERFACE
   --------------------------
   Provide two energy functions:
     symEnergy[state]  -- returns expression symbolic in β (and any
                          coupling constants like J).  Used for the
                          Boltzmann weight in the detailed-balance check.
     numEnergy[state]  -- returns a plain number.  Used for the
                          numerical MCMC validation run.

   ================================================================ *)

BeginPackage["DetailedBalanceChecker`"]

(* Public symbols *)
MetropolisProb; RunWithBits; BuildTransitionMatrix;
CheckDetailedBalance; RunNumericalMCMC; BoltzmannWeights; RunFullCheck;

Begin["`Private`"]

(* ----------------------------------------------------------------
   MetropolisProb
   Standard Metropolis acceptance probability as a symbolic Piecewise.
   deltaE should be a symbolic or numeric energy difference.
   beta must be defined in the caller's scope (or passed explicitly).
   ---------------------------------------------------------------- *)
MetropolisProb[deltaE_] :=
  Piecewise[{{1, deltaE <= 0}, {Exp[-\[Beta] deltaE], deltaE > 0}}]

(* ----------------------------------------------------------------
   RunWithBits
   Run alg[state, readBit] with a specific fixed bit sequence.
   readBit[] reads the next bit; throws if the sequence is exhausted.
   Returns:
     {outcomes, nBitsConsumed}   where outcomes is a list {{p,s},...}
     $OutOfBits                  if the algorithm needed more bits
   ---------------------------------------------------------------- *)
RunWithBits[alg_, state_, bits_List] := Module[
  {pos = 0, readBit, raw},
  readBit[] := (
    pos++;
    If[pos > Length[bits],
      Throw[$OutOfBits, $dbc$tag],
      bits[[pos]]
    ]
  );
  raw = Catch[alg[state, readBit], $dbc$tag, ($OutOfBits &)];
  If[raw === $OutOfBits,
    $OutOfBits,
    (* Normalise to list-of-outcomes form *)
    {If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
       raw,                     (* already {{p,s},...} *)
       {{1, raw}}               (* single state -> prob 1 *)
     ], pos}
  ]
]

(* ----------------------------------------------------------------
   BuildTransitionMatrix
   BFS over all possible bit sequences for a single starting state.
   Bits are extended one level at a time whenever the algorithm
   needs more randomness than the current sequence supplies.
   Returns an Association  {fromState, toState} -> symbolicProbability.
   ---------------------------------------------------------------- *)
Options[BuildTransitionMatrix] = {
  "MaxBitDepth" -> 20,   (* Hard limit on path length *)
  "TimeLimit"   -> 60.,  (* Seconds per starting state *)
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
    queue    = {{}};  (* BFS queue of bit sequences to try *)
    t0       = AbsoluteTime[];
    timedOut = False;

    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[] - t0 > tlim, timedOut = True; Break[]];
      bits  = First[queue];
      queue = Rest[queue];
      res   = RunWithBits[alg, s, bits];

      Which[
        (* Algorithm needs more bits - extend by one level *)
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits, 0], Append[bits, 1]}],

        (* Hit depth limit without terminating *)
        res === $OutOfBits,
          Print["  WARNING: MaxBitDepth=", maxDepth,
                " reached for state ", s, " bits ", bits,
                ". Path excluded - algorithm may not halt on this input."],

        (* Algorithm returned outcome(s) *)
        True,
          {outcomes, k} = res;  (* k = bits consumed *)
          Do[
            With[{p = out[[1]], ns = out[[2]]},
              If[!KeyExistsQ[stateSet, ns],
                Print["  WARNING: Invalid state returned: ", ns],
                matrix[{s, ns}] =
                  Lookup[matrix, Key[{s, ns}], 0] + p * (1/2)^k
              ]
            ],
            {out, outcomes}
          ]
      ]
    ];

    If[timedOut, Print["  WARNING: Time limit exceeded for state ", s]],
    {s, allStates}
  ];

  matrix
]

(* ----------------------------------------------------------------
   CheckDetailedBalance
   For every pair {i,j} verify T(i->j)*pi(i) = T(j->i)*pi(j)
   where pi(s) ∝ Exp[-symEnergy[s]].
   Returns list of violations; empty list means PASS.
   ---------------------------------------------------------------- *)
CheckDetailedBalance[matrix_Association, allStates_List, symEnergy_] := Module[
  {n = Length[allStates], violations = {}, si, sj, tij, tji, res},
  Do[
    si = allStates[[i]]; sj = allStates[[j]];
    tij = Lookup[matrix, Key[{si, sj}], 0];
    tji = Lookup[matrix, Key[{sj, si}], 0];
    res = Simplify[
      tij * Exp[-symEnergy[si]] - tji * Exp[-symEnergy[sj]]
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
   Run the algorithm with a true random bit generator.
   Returns an Association  state -> visitCount.
   ---------------------------------------------------------------- *)
Options[RunNumericalMCMC] = {
  "NSteps"     -> 100000,
  "WarmupFrac" -> 0.1
}

RunNumericalMCMC[allStates_List, alg_, numEnergy_, OptionsPattern[]] := Module[
  {nSteps   = OptionValue["NSteps"],
   warmup   = OptionValue["WarmupFrac"],
   counts, state, raw, outcomes, r, nWarmup},

  nWarmup = Round[nSteps * warmup];
  state   = RandomChoice[allStates];
  counts  = AssociationThread[allStates -> 0];

  step[] := Module[{rb, res2, outs, u, cumP, ns},
    rb[] := RandomInteger[1];
    res2 = RunWithBits[alg, state, Table[rb[], {OptionValue["NSteps"]}]];
    (* For numerical run, just call alg with live random bits *)
    With[{liveRb = Function[{}, RandomInteger[1]]},
      raw = alg[state, liveRb];
      (* normalise to outcome list *)
      outs = If[ListQ[raw] && Length[raw] > 0 && ListQ[raw[[1]]],
               raw, {{1, raw}}];
      (* sample from outcome distribution *)
      u = RandomReal[];
      cumP = 0.;
      ns = Last[outs][[2]]; (* fallback *)
      Do[
        cumP += N[out[[1]]];
        If[u < cumP, ns = out[[2]]; Break[]],
        {out, outs}
      ];
      state = ns
    ]
  ];

  Do[step[], {nWarmup}];
  Do[step[]; If[KeyExistsQ[counts, state], counts[state]++], {nSteps - nWarmup}];

  counts
]

(* ----------------------------------------------------------------
   BoltzmannWeights
   Compute normalised Boltzmann probabilities for each state.
   numEnergy[state] must return a plain number.
   ---------------------------------------------------------------- *)
BoltzmannWeights[allStates_List, numEnergy_] := Module[
  {es, ws, Z},
  es = numEnergy /@ allStates;
  ws = N[Exp[-es]];
  Z  = Total[ws];
  AssociationThread[allStates -> ws / Z]
]

(* ----------------------------------------------------------------
   RunFullCheck
   Top-level entry point.  Runs symbolic + numerical checks and
   prints a self-contained report.
   ---------------------------------------------------------------- *)
Options[RunFullCheck] = Join[
  Options[BuildTransitionMatrix],
  Options[RunNumericalMCMC],
  {"SystemName" -> "Unnamed system"}
]

RunFullCheck[allStates_List, alg_, symEnergy_, numEnergy_, OptionsPattern[]] :=
Module[
  {name = OptionValue["SystemName"],
   n    = Length[allStates],
   matrix, violations, counts, bw, simFreq, kl, row},

  printSep[c_] := Print[StringRepeat[c, 62]];

  printSep["="];
  Print["DETAILED BALANCE CHECKER"];
  printSep["="];
  Print["System : ", name];
  Print["States : ", n, " -- ", allStates];
  printSep["="];

  (* ---- 1. Symbolic check ---- *)
  Print["\n[1/2] SYMBOLIC DETAILED BALANCE CHECK"];
  printSep["-"];
  Print["Building decision trees for all ", n, " states..."];

  matrix = BuildTransitionMatrix[allStates, alg,
    "MaxBitDepth" -> OptionValue["MaxBitDepth"],
    "TimeLimit"   -> OptionValue["TimeLimit"],
    "Verbose"     -> OptionValue["Verbose"]
  ];

  Print["Transition matrix built. Non-zero entries: ", Length[matrix]];
  Print["Checking ", Binomial[n, 2], " state pairs for detailed balance..."];

  violations = CheckDetailedBalance[matrix, allStates, symEnergy];

  If[violations === {},
    Print["\n  RESULT: PASS -- all ", Binomial[n,2],
          " pairs satisfy detailed balance exactly."],
    Print["\n  RESULT: FAIL -- ", Length[violations], " violation(s):"];
    Scan[Function[v,
      Print["    ", v["pair"], "  residual = ", v["residual"]]
    ], violations]
  ];

  (* ---- 2. Numerical check ---- *)
  Print["\n[2/2] NUMERICAL MCMC CHECK"];
  printSep["-"];
  Print["Running ", OptionValue["NSteps"], " MCMC steps..."];

  counts = RunNumericalMCMC[allStates, alg, numEnergy,
    "NSteps"     -> OptionValue["NSteps"],
    "WarmupFrac" -> OptionValue["WarmupFrac"]
  ];
  bw      = BoltzmannWeights[allStates, numEnergy];
  simFreq = # / Total[counts] & /@ counts;

  Print["\n  ", PaddedForm["State", 20],
              PaddedForm["Simulated", 12],
              PaddedForm["Boltzmann", 12]];
  Print["  ", StringRepeat["-", 44]];
  Scan[Function[s,
    Print["  ", PaddedForm[ToString[s], 20],
               PaddedForm[NumberForm[N@simFreq[s], {4,4}], 12],
               PaddedForm[NumberForm[N@bw[s],      {4,4}], 12]]
  ], allStates];

  (* KL divergence D(sim || Boltzmann) *)
  kl = Total[Table[
    With[{p = N@simFreq[s], q = N@bw[s]},
      If[p > 0 && q > 0, p Log[p/q], 0.]],
    {s, allStates}]];

  Print["\n  KL divergence (sim || Boltzmann) = ", NumberForm[kl, {5,5}]];
  Print["  Numerical verdict: ",
        If[kl < 0.02, "CONSISTENT with Boltzmann.",
                      "WARNING -- significant deviation from Boltzmann."]];

  printSep["="];
  Print["END OF REPORT\n"];

  <|"matrix" -> matrix, "violations" -> violations,
    "counts" -> counts, "boltzmann" -> bw|>
]

End[]
EndPackage[]
