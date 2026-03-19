(* ================================================================
   Example 4: Variable-bit two-speed Metropolis on a ring -- PASSES
   ================================================================

   System: L=4 sites on a ring, 1 particle.
   State:  integer 1..4 (particle position).
   Energy: E(i) = eps_i (symbolic), or numBeta * numEps[[i]] (numeric).

   Move: the FIRST random bit selects one of two move types, so
   different execution paths consume different numbers of bits:

     bit1 = 0  →  "standard hop"  (2 bits consumed total)
       bit2 = 0 : propose left  neighbour → Metropolis
       bit2 = 1 : propose right neighbour → Metropolis

     bit1 = 1  →  "lazy hop"  (2 OR 3 bits consumed)
       bit2 = 0 : stay put immediately        ← DEPTH-2 LEAF
       bit2 = 1 : read one more bit for direction → Metropolis
         bit3 = 0 : propose left  neighbour   ← DEPTH-3 LEAF
         bit3 = 1 : propose right neighbour   ← DEPTH-3 LEAF

   Probability of proposing each ±1 neighbour:
     via standard hop : (1/2)(1/2)       = 1/4
     via lazy hop     : (1/2)(1/2)(1/2)  = 1/8
     total            : 3/8  (symmetric in both directions)

   Because the proposal is symmetric and Metropolis acceptance is
   used, detailed balance holds exactly.  The self-loop from the
   lazy bit2=0 branch only affects T(i→i) and does not disturb any
   (i,j) balance equation.
   ================================================================ *)

(* ---- System parameters ---- *)
L$vb         = 4
allStates$vb = Range[L$vb]

rightOf$vb[s_Integer] := Mod[s,     L$vb] + 1   (* clockwise  neighbour *)
leftOf$vb[s_Integer]  := Mod[s - 2, L$vb] + 1   (* anti-clock neighbour *)

(* Symbolic bare energies (no beta) *)
eps$sym$vb = Table[Symbol["eps" <> ToString[i]], {i, L$vb}]
symEnergy$vb[s_Integer] := eps$sym$vb[[s]]

(* Numeric energies with beta folded in *)
numBeta$vb = 1.2
numEps$vb  = {0.0, 0.5, 1.5, 1.0}
numEnergy$vb[s_Integer] := numBeta$vb * numEps$vb[[s]]

(* ================================================================
   Symbolic algorithm (used for the exact DB check).
   Uses MetropolisProb so beta stays symbolic.
   ================================================================ *)
TwoSpeedSym[state_Integer, readBit_] := Module[
  {b1, b2, b3, nbr, dE, p},
  b1 = readBit[];
  If[b1 == 0,

    (* ---- Standard hop: 1 more bit picks direction ---- *)
    b2  = readBit[];
    nbr = If[b2 == 0, leftOf$vb[state], rightOf$vb[state]];
    dE  = symEnergy$vb[nbr] - symEnergy$vb[state];
    p   = MetropolisProb[dE];
    {{p, nbr}, {1 - p, state}},

    (* ---- Lazy hop: gatekeeper bit ---- *)
    b2 = readBit[];
    If[b2 == 0,
      (* Stay put -- no further bits needed; this is a depth-2 leaf *)
      {{1, state}},
      (* Willing to move: read direction bit, then Metropolis *)
      b3  = readBit[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = symEnergy$vb[nbr] - symEnergy$vb[state];
      p   = MetropolisProb[dE];
      {{p, nbr}, {1 - p, state}}
    ]
  ]
]

(* ================================================================
   Numerical algorithm (used for the MCMC frequency check).
   beta is already folded into numEnergy.
   ================================================================ *)
TwoSpeedNum[state_Integer, readBit_] := Module[
  {b1, b2, b3, nbr, dE, p},
  b1 = readBit[];
  If[b1 == 0,
    b2  = readBit[];
    nbr = If[b2 == 0, leftOf$vb[state], rightOf$vb[state]];
    dE  = numEnergy$vb[nbr] - numEnergy$vb[state];
    p   = N @ If[dE <= 0, 1, Exp[-dE]];
    {{p, nbr}, {1 - p, state}},
    b2 = readBit[];
    If[b2 == 0,
      {{1, state}},
      b3  = readBit[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = numEnergy$vb[nbr] - numEnergy$vb[state];
      p   = N @ If[dE <= 0, 1, Exp[-dE]];
      {{p, nbr}, {1 - p, state}}
    ]
  ]
]
