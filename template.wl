(* ================================================================
   template.wl  --  Minimal template for the detailed-balance checker
   ================================================================
   Copy this file, fill in the four required definitions, and run:

     wolframscript -file check.wls your_algorithm.wl MaxBitString=<bits>

   REQUIRED:  energy, Algorithm, BitsToState, numBeta
   OPTIONAL:  symParams or DynamicSymParams (for symbolic parameters)
              DisplayState (for custom state display in output)
   ================================================================ *)


(* ---- Energy function ---------------------------------------------------- *)
(*
   Takes a state, returns the bare energy (a real number, no beta).
   Must be deterministic -- no random calls.
   May use symbolic variables listed in symParams or DynamicSymParams.
*)
energy[state_] :=
  0   (* replace with your energy *)


(* ---- MCMC algorithm ------------------------------------------------------- *)
(*
   Takes a state, returns the next state.
   Use native Mathematica random calls: RandomInteger[], RandomReal[],
   RandomChoice[]. For Metropolis acceptance, use exactly:

     If[RandomReal[] < MetropolisProb[dE], newState, state]

   MetropolisProb[dE] is defined by the checker as:
     Piecewise[{{1, dE<=0}, {Exp[-beta*dE], dE>0}}]
   with beta kept symbolic during the symbolic check.
*)
Algorithm[state_] :=
  state   (* replace with your algorithm *)


(* ---- Symbolic energy parameters (optional) ------------------------------ *)
(*
   Option A -- static parameters (same for every connected component):

     symParams = <|"eps" -> {eps1, eps2}, "couplings" -> {J1, J2}|>

   Option B -- dynamic parameters (auto-generated per component, e.g.
   Kawasaki coupling constants that depend on which particle types
   are present):

     DynamicSymParams[states_List] :=
       <|"couplings" -> ... symbols derived from states ...|>

   Both options may be used together. If neither is defined, the
   checker treats the energy as fully numeric.

   Parameters declared here are kept symbolic during the symbolic
   check (FullSimplify with beta > 0) and assigned random numeric
   values for the MCMC check.

   Example for a 2-particle Kawasaki system:
     DynamicSymParams[states_List] :=
       Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
         <|"couplings" ->
           Flatten @ Table[
             If[a < b, ToExpression["J" <> ToString[a] <> ToString[b]], Nothing],
             {a, types}, {b, types}]|>]
*)


(* ---- BitsToState --------------------------------------------------------- *)
(*
   Converts a bit string (list of 0s and 1s) to a seed state for the
   checker, or returns None to skip this bit string.

   The checker calls BitsToState for every bit string from length 1 up
   to MaxBitString, in order of increasing length. It uses BFS from the
   returned state to discover the full connected component automatically.

   Requirements:
   - Accept a list of any length; return None for invalid patterns.
   - Different bit strings that return non-None should map to states in
     different connected components (to avoid redundant work).
   - The state must be accepted by both energy and Algorithm.

   Simple example: 4-site ring, state = list of site occupancies (0/1)
     BitsToState[bits_List] :=
       If[Length[bits] =!= 4, None, bits]

   Labeled-particle example (kawasaki_1d.wl encoding):
     BitsToState[bits_List] :=
       Module[{id = FromDigits[bits, 2]},
         If[id == 0, None, $decode[id]]]
*)
BitsToState[bits_List] :=
  None   (* replace with your BitsToState *)


(* ---- Inverse temperature ------------------------------------------------ *)
numBeta = 1


(* ---- DisplayState (optional) -------------------------------------------- *)
(*
   If defined, the checker calls DisplayState[state] instead of
   ToString[state] for the State column in the output table.
   Useful for 2D states, e.g.:
     DisplayState[state_List] :=
       With[{L = Round[Sqrt[Length[state]]]},
         StringJoin @ Riffle[
           Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
                 {r, 1, L}], "|"]]
*)
