#!/usr/bin/env wolframscript
(* ================================================================
   animate.wls  --  Terminal animation of an MCMC algorithm
   ================================================================
   Usage:
     wolframscript -file animate.wls <alg_file> \
       Sites=<n> N=<n> [Steps=200] [Beta=1] [Delay=0.1] [J12=0.5 ...]

   Arguments:
     Sites=<n>     Total number of lattice sites (e.g. 6 for a 6-site
                   ring; 9 for a 3x3 grid). Must be >= N.
     N=<n>         Number of labeled particles (types 1..N).
     Steps=<n>     Number of MCMC steps to animate (default 200).
     Beta=<f>      Inverse temperature (overrides numBeta in .wl file).
     Delay=<f>     Pause between frames in seconds (default 0.1).
     J<a><b>=<f>   Set coupling constant J_ab to a specific value.
                   Any J symbol not supplied is assigned a random value.

   Examples:
     wolframscript -file animate.wls examples3/kawasaki_1d.wl \
       Sites=6 N=3 Steps=300 Beta=1 Delay=0.15 J12=1.0 J23=0.5

     wolframscript -file animate.wls examples3/kawasaki_2d.wl \
       Sites=9 N=4 Steps=200 Delay=0.2
   ================================================================ *)

(* MetropolisProb for actual (non-intercepted) execution.
   Uses numBeta which is set below after loading the .wl file. *)
MetropolisProb[dE_] := If[dE <= 0, 1., Exp[-$animBeta * dE]]

scriptDir = DirectoryName[$InputFileName];

(* ---- Parse arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "animate.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1 || !MemberQ[kvArgs, _?(StringMatchQ[#, "Sites=" ~~ __] &)] ||
                          !MemberQ[kvArgs, _?(StringMatchQ[#, "N=" ~~ __] &)],
  Print["Usage: wolframscript -file animate.wls <alg_file> Sites=<n> N=<n> [options]"];
  Print["  Sites=<n>     total lattice sites"];
  Print["  N=<n>         number of labeled particles (types 1..N)"];
  Print["  Steps=<n>     animation steps (default 200)"];
  Print["  Beta=<f>      inverse temperature (default from .wl file)"];
  Print["  Delay=<f>     seconds per frame (default 0.1)"];
  Print["  J<a><b>=<f>   coupling constant (random if omitted)"];
  Exit[1]];

(* ---- Option defaults ---- *)
nSites   = 4;
nPart    = 2;
nSteps   = 200;
animBeta = None;    (* None means use numBeta from the .wl file *)
delay    = 0.1;

Do[
  If[StringMatchQ[opt, "Sites=" ~~ __], nSites   = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "N="     ~~ __], nPart    = ToExpression[StringDrop[opt,  2]]];
  If[StringMatchQ[opt, "Steps=" ~~ __], nSteps   = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "Beta="  ~~ __], animBeta = ToExpression[StringDrop[opt,  5]]];
  If[StringMatchQ[opt, "Delay=" ~~ __], delay    = ToExpression[StringDrop[opt,  6]]],
  {opt, kvArgs}];

(* ---- Load algorithm file ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState];
Get[algFile];

If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined."]; Exit[1]];

(* Use Beta from command line if given, otherwise from the .wl file *)
$animBeta = If[animBeta =!= None, animBeta, numBeta];

(* ---- Build initial state ---- *)
If[nPart > nSites,
  Print["ERROR: N=", nPart, " cannot exceed Sites=", nSites]; Exit[1]];
SeedRandom[1234];
initialState = RandomSample[Join[Range[nPart], ConstantArray[0, nSites - nPart]], nSites];

(* ---- Assign coupling constants ---- *)
(* Parse any J<a><b>=<f> args from the command line *)
jSupplied = Association[
  Map[Function[opt,
    With[{parts = StringSplit[opt, "=", 2]},
      ToExpression[parts[[1]]] -> ToExpression[parts[[2]]]]],
    Select[kvArgs, StringMatchQ[#, "J" ~~ DigitCharacter.. ~~ "="] ||
                   StringMatchQ[#, "J" ~~ DigitCharacter.. ~~ DigitCharacter.. ~~ "="] ||
                   StringMatchQ[#, "J" ~~ DigitCharacter ~~ DigitCharacter ~~ "=" ~~ __] &]]];

(* Assign all J_ab for pairs present in the initial state *)
types = Sort[DeleteCases[initialState, 0]];
SeedRandom[42];
Do[If[a < b,
  sym = ToExpression["J" <> ToString[a] <> ToString[b]];
  Set[Evaluate[sym], Lookup[jSupplied, sym, RandomReal[{0.5, 1.5}]]]],
  {a, types}, {b, types}];

(* ---- Display helper ---- *)
showState[step_, state_] := Module[{stateStr, jStr},
  (* Coupling summary (compact) *)
  jStr = StringRiffle[
    Select[Map[Function[{a},
      Select[Map[Function[{b},
        If[a < b,
          With[{sym = ToExpression["J" <> ToString[a] <> ToString[b]]},
            If[NumericQ[sym], "J" <> ToString[a] <> ToString[b] <> "=" <>
               ToString[NumberForm[N[sym], {3,2}]], Nothing]],
          Nothing]],
        types], StringQ]], types], StringQ],
    "  "];
  stateStr = If[Length[DownValues[DisplayState]] > 0,
    DisplayState[state], StringRiffle[ToString /@ state, " "]];
  WriteString["stdout", "\033[H\033[2J"];   (* clear screen *)
  Print["Algorithm: ", algFile];
  Print["Sites: ", nSites, "  N: ", nPart, "  Beta: ", $animBeta];
  If[jStr =!= "", Print["Couplings: ", jStr]];
  Print[""];
  Print["Step ", step, " / ", nSteps];
  Print[""];
  Print[stateStr];
  Print[""]];

(* ---- Run animation ---- *)
state = initialState;
showState[0, state];
Do[
  state = Algorithm[state];
  showState[i, state];
  Pause[delay],
  {i, 1, nSteps}];

Print["Done."];
