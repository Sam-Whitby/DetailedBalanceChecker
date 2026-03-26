(* ================================================================
   Example: Asymmetric proposal with Metropolis -- [FAIL]
   ================================================================

   System: L=4 ring, 1 particle.  States {1,2,3,4}.

   The algorithm proposes moves with a BIASED direction distribution:
     P(propose right) = 3/4,   P(propose left) = 1/4.

   This is encoded in a variable-depth bit tree:
     bit1=0          -> propose right       (prob 1/2)
     bit1=1, bit2=0  -> propose right       (prob 1/4)
     bit1=1, bit2=1  -> propose left        (prob 1/4)

   The Metropolis acceptance is applied afterwards.

   Why this FAILS detailed balance:
     T(i->right(i)) = (3/4) * MetropolisProb[\[CapitalDelta]E_right]
     T(right(i)->i) = (1/4) * MetropolisProb[-\[CapitalDelta]E_right]

   For \[CapitalDelta]E > 0 (uphill move):
     T(i->right(i)) * \[Pi](i) = (3/4) * exp(-\[Beta]*\[CapitalDelta]E) * exp(-\[Beta]*E_i)
     T(right(i)->i) * \[Pi](r) = (1/4) * 1                * exp(-\[Beta]*E_r)
     Ratio = 3/4 * exp(-\[Beta]*\[CapitalDelta]E) * exp(-\[Beta]*E_i)
           / (1/4 * exp(-\[Beta]*(E_i+\[CapitalDelta]E)))
           = (3/4)/(1/4) = 3  != 1.  [FAIL]
   ================================================================ *)

Get[DirectoryName[$InputFileName] <> "barker_ring.wl"]   (* reuse same system *)

(* Asymmetric proposal + Metropolis on L=4 ring *)
AsymmetricProposal[state_Integer] := Module[
  {b1, b2, nbr, dE},
  b1 = RandomInteger[];
  nbr = If[b1 == 0,
    rightOf$br[state],                         (* 0 -> right (50%) *)
    b2 = RandomInteger[];
    If[b2 == 0, rightOf$br[state],             (* 1,0 -> right (25%) *)
                leftOf$br[state]]              (* 1,1 -> left  (25%) *)
  ];
  dE = energy$br[nbr] - energy$br[state];
  If[RandomReal[] < MetropolisProb[dE], nbr, state]
]
