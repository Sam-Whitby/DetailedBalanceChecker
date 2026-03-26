(* ================================================================
   Example 5: Variable-bit biased rightward drift -- FAILS
   ================================================================

   Same system as variable_bit_pass.wl.

   The ONLY change from the passing version is in the bit1=1, bit2=0
   branch.  Instead of staying put (a harmless self-loop), the
   algorithm now ALWAYS hops one step clockwise unconditionally,
   bypassing Metropolis entirely:

     bit1 = 0  ->  standard hop  (unchanged, same as PASS)
     bit1 = 1  ->
       bit2 = 0 : hop RIGHT unconditionally, no Metropolis  <- BUG
       bit2 = 1 : read bit3 for direction, apply Metropolis (unchanged)

   Why this breaks detailed balance
   ----------------------------------
   From state i, the probability of transitioning to right(i) is:
     T(i -> right(i)) = 1/4 * A(i->r)    [standard hop, going right]
                      + 1/4              [biased hop,   always right]
                      + 1/8 * A(i->r)    [lazy fair hop, going right]
                      = 1/4 + 3/8 * A(i->r)

   But from right(i) back to i:
     T(right(i) -> i) = 3/8 * A(right(i)->i)

   The extra 1/4 term cannot be cancelled, so every adjacent pair
   violates detailed balance.
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "variable_bit_pass.wl"]

(* ================================================================
   Biased algorithm: bit1=1, bit2=0 hops right unconditionally
   ================================================================ *)
Biased[state_Integer] := Module[
  {b1, b2, b3, nbr, dE},
  b1 = RandomInteger[];
  If[b1 == 0,

    (* Standard hop: unchanged from PASS *)
    b2  = RandomInteger[];
    nbr = If[b2 == 0, leftOf$vb[state], rightOf$vb[state]];
    dE  = energy$vb[nbr] - energy$vb[state];
    If[RandomReal[] < MetropolisProb[dE], nbr, state],

    (* Lazy branch -- second 0 now hops right unconditionally *)
    b2 = RandomInteger[];
    If[b2 == 0,
      (* BUG: always move clockwise, no energy check *)
      rightOf$vb[state],
      (* Second 1: fair hop, unchanged from PASS *)
      b3  = RandomInteger[];
      nbr = If[b3 == 0, leftOf$vb[state], rightOf$vb[state]];
      dE  = energy$vb[nbr] - energy$vb[state];
      If[RandomReal[] < MetropolisProb[dE], nbr, state]
    ]
  ]
]
