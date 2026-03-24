(* ================================================================
   1D Kawasaki bond-swap dynamics -- PASSES detailed balance
   ================================================================

   System: L=4 ring, 3 particles, 1 hole.
   State:  sorted triple {p1,p2,p3}, 1 <= p1 < p2 < p3 <= 4.
           Four states: {1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}.

   Energy: site energies + nearest-neighbour pairwise coupling
           E({p1,p2,p3}) = eps[[p1]] + eps[[p2]] + eps[[p3]]
                         + J * (number of adjacent occupied pairs)
           All energy parameters are symbolic; random values are
           assigned automatically during the numerical MCMC check.

   Move:   bond-swap Kawasaki.
           Pick a random bond (RandomInteger[{0,3}], 2 bits exact).
           Bond b connects site b+1 to site Mod[b+1,L]+1.
           - particle-hole bond: particle hops, Metropolis acceptance.
           - particle-particle bond: exchange = identity on sorted state.
           - hole-hole bond: impossible with 3 particles on L=4.

   Proposal is symmetric (each bond equally likely in both directions)
   and Metropolis is correct, so detailed balance holds exactly.
   ================================================================ *)

(* ---- System parameters ---- *)
L$sw = 4

(* Symbolic site energies and NN coupling.
   Random numerical values are assigned automatically during the
   numerical MCMC check when symParams is provided. *)
eps$sw    = {\[Epsilon]sw1, \[Epsilon]sw2, \[Epsilon]sw3, \[Epsilon]sw4}
symParams = <|"eps" -> eps$sw, "couplings" -> {Jsw1}|>

numBeta = 1

(* Ring minimum-image distance for site pair (a,b) on L-site ring *)
ringDistSw[a_Integer, b_Integer] := Min[Abs[a - b], L$sw - Abs[a - b]]

(* Energy: site energies + pairwise NN coupling *)
energy[{p1_Integer, p2_Integer, p3_Integer}] :=
  eps$sw[[p1]] + eps$sw[[p2]] + eps$sw[[p3]] +
  Jsw1 * (If[ringDistSw[p1, p2] == 1, 1, 0] +
             If[ringDistSw[p1, p3] == 1, 1, 0] +
             If[ringDistSw[p2, p3] == 1, 1, 0])

(* ================================================================
   Algorithm: bond-swap Kawasaki.
   RandomInteger[{0,3}] picks bond 0..3 (2 bits, no rejection).
   Bond b: sites b+1 and Mod[b+1,L$sw]+1.
   If one site occupied, one empty: propose hop with Metropolis.
   If both occupied: exchange = no state change (identity on sorted state).
   ================================================================ *)
Algorithm[{p1_Integer, p2_Integer, p3_Integer}] :=
  Module[{b, s1, s2, occ, newOcc, dE},
    occ = {p1, p2, p3};
    (* Pick bond *)
    b  = RandomInteger[{0, L$sw - 1}];
    s1 = b + 1;
    s2 = Mod[b + 1, L$sw] + 1;
    Which[
      (* Both occupied: swap = identity on sorted state *)
      MemberQ[occ, s1] && MemberQ[occ, s2],
        Sort[occ],
      (* s1 occupied, s2 empty: particle at s1 hops to s2 *)
      MemberQ[occ, s1],
        newOcc = Sort[Append[DeleteCases[occ, s1], s2]];
        dE = energy[newOcc] - energy[Sort[occ]];
        If[RandomReal[] < MetropolisProb[dE], newOcc, Sort[occ]],
      (* s1 empty, s2 occupied: particle at s2 hops to s1 *)
      MemberQ[occ, s2],
        newOcc = Sort[Append[DeleteCases[occ, s2], s1]];
        dE = energy[newOcc] - energy[Sort[occ]];
        If[RandomReal[] < MetropolisProb[dE], newOcc, Sort[occ]],
      (* Both empty: cannot happen with 3 particles on L=4 *)
      True, Sort[occ]
    ]
  ]

(* ---- Checker interface ---- *)
seedState = {1, 2, 3}
