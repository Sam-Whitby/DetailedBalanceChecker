(* ================================================================
   Example: Barker criterion on a ring -- Generalised API  [PASS]
   ================================================================

   System: L=4 ring, 1 particle.  States {1,2,3,4}.

   Barker acceptance (also called the "Glauber" or "heat-bath"
   criterion) uses the acceptance probability:

       p_Barker(i->j) = 1 / (1 + exp(\[Beta] * \[CapitalDelta]E))

   rather than the Metropolis min(1, exp(-\[Beta]*\[CapitalDelta]E)).

   Detailed balance check:
     T(i->j)*\[Pi](i) = p_Barker(i->j) * exp(-\[Beta]*E_i) / 2
                       = exp(-\[Beta]*E_i) / (2*(1+exp(\[Beta]*(E_j-E_i))))
     T(j->i)*\[Pi](j) = p_Barker(j->i) * exp(-\[Beta]*E_j) / 2
                       = exp(-\[Beta]*E_j) / (2*(1+exp(\[Beta]*(E_i-E_j))))
     Ratio = exp(-\[Beta]*(E_i-E_j)) * (1+exp(\[Beta]*(E_i-E_j)))
                                      / (1+exp(-\[Beta]*(E_i-E_j)))
           = exp(-\[Beta]*(E_i-E_j)) * exp(\[Beta]*(E_i-E_j))  =  1  [PASS]

   This example demonstrates that the checker correctly verifies
   acceptance criteria other than Metropolis.
   ================================================================ *)

L$br       = 4
eps$br     = {0, 1, 3, 2}   (* bare site energies -- exact integers *)
numBeta$br = 1

rightOf$br[s_Integer] := Mod[s,     L$br] + 1
leftOf$br[s_Integer]  := Mod[s - 2, L$br] + 1

energy$br[s_Integer] := eps$br[[s]]

(* Barker criterion: p = 1/(1 + exp(\[Beta]*dE))
   Uses the global \[Beta] symbol so it stays symbolic during the
   DB check and evaluates numerically during MCMC via Block. *)
BarkerProb[dE_] := 1 / (1 + Exp[\[Beta] * dE])

BarkerRing[state_Integer, readBit_, acceptTest_] := Module[
  {b, nbr, dE},
  b   = readBit[];
  nbr = If[b == 0, leftOf$br[state], rightOf$br[state]];
  dE  = energy$br[nbr] - energy$br[state];
  If[acceptTest[BarkerProb[dE]] == 1, nbr, state]
]
