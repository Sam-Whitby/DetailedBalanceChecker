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
(* Symbolic site energies -- kept unassigned during the symbolic check.
   Random numerical values are assigned automatically by RunFullCheck
   when "SysParams" -> params$br is supplied. *)
eps$br     = {\[Epsilon]br1, \[Epsilon]br2, \[Epsilon]br3, \[Epsilon]br4}
params$br  = <|"L" -> L$br, "eps" -> eps$br|>
numBeta$br = 1

rightOf$br[s_Integer] := Mod[s,     L$br] + 1
leftOf$br[s_Integer]  := Mod[s - 2, L$br] + 1

energy$br[s_Integer] := eps$br[[s]]

(* Barker criterion: p = 1/(1 + exp(\[Beta]*dE))
   Uses the global \[Beta] symbol so it stays symbolic during the
   DB check and evaluates numerically during MCMC via Block. *)
BarkerProb[dE_] := 1 / (1 + Exp[\[Beta] * dE])

BarkerRing[state_Integer] := Module[
  {b, nbr, dE},
  b   = RandomInteger[];
  nbr = If[b == 0, leftOf$br[state], rightOf$br[state]];
  dE  = energy$br[nbr] - energy$br[state];
  If[RandomReal[] < BarkerProb[dE], nbr, state]
]
