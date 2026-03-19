(* ================================================================
   Example 4: Variable-bit two-speed Metropolis on a ring -- PASSES
   ================================================================

   System: L=4 sites on a ring, 1 particle.
   State:  integer 1..4 (particle position).
   Energy: E(i) = eps_i (exact integers, no beta).

   Move: the FIRST random bit selects one of two move types, so
   different execution paths consume different numbers of bits:

     bit1 = 0  ->  "standard hop"  (2 bits consumed total)
       bit2 = 0 : propose left  neighbour -> Metropolis
       bit2 = 1 : propose right neighbour -> Metropolis

     bit1 = 1  ->  "lazy hop"  (2 OR 3 bits consumed)
       bit2 = 0 : stay put immediately        <- DEPTH-2 LEAF
       bit2 = 1 : read one more bit for direction -> Metropolis
         bit3 = 0 : propose left  neighbour   <- DEPTH-3 LEAF
         bit3 = 1 : propose right neighbour   <- DEPTH-3 LEAF

   Probability of proposing each +-1 neighbour:
     via standard hop : (1/2)(1/2)       = 1/4
     via lazy hop     : (1/2)(1/2)(1/2)  = 1/8
     total            : 3/8  (symmetric in both directions)

   Because the proposal is symmetric and Metropolis acceptance is
   used, detailed balance holds exactly.
   ================================================================ *)

(* ---- System parameters ---- *)
L$vb       = 4
eps$vb     = {0, 1, 3, 2}   (* exact integers, no beta *)
numBeta$vb = 1

rightOf$vb[s_Integer] := Mod[s,     L$vb] + 1
leftOf$vb[s_Integer]  := Mod[s - 2, L$vb] + 1

energy$vb[s_Integer] := eps$vb[[s]]

(* ================================================================
   Single algorithm: works symbolically (beta unassigned) and
   numerically (Block[{beta = numBeta$vb}]).
   ================================================================ *)
TwoSpeed[state_Integer] := Module[
  {b1, b2, b3, nbr, dE},
  b1 = RandomInteger[];
  If[b1 == 0,

    (* ---- Standard hop: 1 more random integer picks direction ---- *)
    b2  = RandomInteger[];
    nbr = If[b2 == 0, leftOf$vb[state], rightOf$vb[state]];
    dE  = energy$vb[nbr] - energy$vb[state];
    If[RandomReal[] < MetropolisProb[dE], nbr, state],

    (* ---- Lazy hop: gatekeeper random integer ---- *)
    b2 = RandomInteger[];
    If[b2 == 0,
      (* Stay put -- depth-2 leaf *)
      state,
      (* Willing to move: read direction, then Metropolis *)
      b3  = RandomInteger[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = energy$vb[nbr] - energy$vb[state];
      If[RandomReal[] < MetropolisProb[dE], nbr, state]
    ]
  ]
]
