(* ================================================================
   Example 5: Variable-bit biased rightward drift -- FAILS
   ================================================================

   Same system as variable_bit_pass.wl.

   The ONLY change from the passing version is in the bit1=1, bit2=0
   branch.  Instead of staying put (a harmless self-loop), the
   algorithm now ALWAYS hops one step clockwise with probability 1,
   bypassing Metropolis entirely:

     bit1 = 0  →  standard hop  (unchanged, same as PASS)
     bit1 = 1  →
       bit2 = 0 : hop RIGHT unconditionally, no Metropolis  ← BUG
       bit2 = 1 : read bit3 for direction, apply Metropolis (unchanged)

   Why this breaks detailed balance
   ─────────────────────────────────
   From state i, the probability of transitioning to right(i) is:
     T(i → right(i)) = 1/4 · A(i→r)    [standard hop, going right]
                     + 1/4              [biased hop,   always right, no accept]
                     + 1/8 · A(i→r)    [lazy fair hop, going right]
                     = 1/4 + 3/8 · A(i→r)

   But from right(i) back to i:
     T(right(i) → i) = 3/8 · A(right(i)→i)
   (the biased hop from right(i) goes to right(right(i)), not to i)

   So T(i→r)·π(i)  −  T(r→i)·π(r)
     = [1/4 + 3/8·A(i→r)]·exp(−β E_i)  −  3/8·A(r→i)·exp(−β E_r)

   The extra 1/4 term cannot be cancelled, so every adjacent pair
   violates detailed balance.  Physically the algorithm drives a
   net clockwise current around the ring.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "variable_bit_pass.wl"]

(* ================================================================
   Symbolic version with the biased drift bug
   ================================================================ *)
BiasedSym[state_Integer, readBit_] := Module[
  {b1, b2, b3, nbr, dE, p},
  b1 = readBit[];
  If[b1 == 0,

    (* Standard hop: unchanged from PASS *)
    b2  = readBit[];
    nbr = If[b2 == 0, leftOf$vb[state], rightOf$vb[state]];
    dE  = symEnergy$vb[nbr] - symEnergy$vb[state];
    p   = MetropolisProb[dE];
    {{p, nbr}, {1 - p, state}},

    (* Lazy branch -- bit2=0 now hops right unconditionally *)
    b2 = readBit[];
    If[b2 == 0,
      (* BUG: always move clockwise, no energy check *)
      {{1, rightOf$vb[state]}},
      (* bit2=1: fair hop, unchanged from PASS *)
      b3  = readBit[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = symEnergy$vb[nbr] - symEnergy$vb[state];
      p   = MetropolisProb[dE];
      {{p, nbr}, {1 - p, state}}
    ]
  ]
]

(* ================================================================
   Numerical version with the biased drift bug
   ================================================================ *)
BiasedNum[state_Integer, readBit_] := Module[
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
      {{1, rightOf$vb[state]}},
      b3  = readBit[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = numEnergy$vb[nbr] - numEnergy$vb[state];
      p   = N @ If[dE <= 0, 1, Exp[-dE]];
      {{p, nbr}, {1 - p, state}}
    ]
  ]
]
