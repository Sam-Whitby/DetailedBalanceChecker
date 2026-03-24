(* ================================================================
   2D Kawasaki with unbalanced proposal -- FAILS detailed balance
   ================================================================

   System: 3x3 periodic lattice (torus), 3 particles, 6 holes.
   State:  sorted triple {s1,s2,s3}, sites 1..9.
   Energy: same rational energies as kawasaki_2d_pass.wl.

   The bug
   -------
   Particle selection uses an unbalanced if-tree that consumes
   different numbers of random bits per branch:

     b1 = 0 (prob 1/2)  --> select particle 1 (smallest site)
     b1 = 1 (prob 1/2):
       b2 = 0 (prob 1/4) --> select particle 2
       b2 = 1 (prob 1/4):
         b3 = 0 (prob 1/8) --> select particle 3
         b3 = 1 (prob 1/8) --> select particle 3 again (extra branch)

   This gives effective proposal probabilities:
     p1: 1/2,  p2: 1/4,  p3: 3/8  (sums to > 1 -- actually
   the extra b3=1 branch is a "retry" that doubles p3's weight).

   After normalising by the total probability (which is the same
   from every state), the proposal is ASYMMETRIC: particle 1 is
   systematically under-proposed relative to particle 3 whenever
   a hop changes the rank ordering of the particles.

   Why this breaks detailed balance
   ---------------------------------
   Consider a hop that moves a particle from a high-index site to a
   lower-index site, changing that particle's rank from 3 -> 1 (or
   2 -> 1).  The forward move uses the rank-3 proposal weight (3/8)
   while the reverse move uses the rank-1 weight (1/2), so

     T(i->j) * pi(i) != T(j->i) * pi(j).

   The bug is subtle because all particles CAN be proposed (no
   direction is excluded), yet the asymmetric bit consumption
   creates an unequal proposal distribution.
   ================================================================ *)

(* ---- Grid helpers (copied from 2d pass -- self-contained) ---- *)
$Lx2f = 3; $Ly2f = 3;

siteRow2Df[s_Integer] := Ceiling[s / $Lx2f]
siteCol2Df[s_Integer] := Mod[s - 1, $Lx2f] + 1
siteIdx2Df[r_Integer, c_Integer] := (r - 1) * $Lx2f + c

leftOf2Df[s_Integer]  := siteIdx2Df[siteRow2Df[s], Mod[siteCol2Df[s] - 2, $Lx2f] + 1]
rightOf2Df[s_Integer] := siteIdx2Df[siteRow2Df[s], Mod[siteCol2Df[s],     $Lx2f] + 1]
upOf2Df[s_Integer]    := siteIdx2Df[Mod[siteRow2Df[s] - 2, $Ly2f] + 1, siteCol2Df[s]]
downOf2Df[s_Integer]  := siteIdx2Df[Mod[siteRow2Df[s],     $Ly2f] + 1, siteCol2Df[s]]
neighbors2Df[s_Integer] := {leftOf2Df[s], rightOf2Df[s], upOf2Df[s], downOf2Df[s]}
adjacentQ2Df[a_Integer, b_Integer] := MemberQ[neighbors2Df[a], b]

(* ---- Energy (same rational values as the 2d pass example) ---- *)
$siteE2f = {0, 1/2, 0, 1/2, 1, 1/2, 0, 1/2, 0}
$nnJ2f   = 1/2

energy[occ_List] :=
  Total[$siteE2f[[#]] & /@ occ] +
  $nnJ2f * Length[Select[Subsets[occ, {2}],
                         adjacentQ2Df[#[[1]], #[[2]]] &]]

(* ================================================================
   Algorithm: subtly broken particle selection.
   The if-tree below creates unequal proposal probabilities via
   different bit-consumption depths per branch.
   ================================================================ *)
Algorithm[occ_List] :=
  Module[{sorted, b1, b2, b3, mover, dir, target, others, newOcc, dE},
    sorted = Sort[occ];

    (* --- Broken particle selection: unbalanced tree ---
       Correct code would be: RandomInteger[{0,2}]
       This code appears similar but has asymmetric bit depths. *)
    b1 = RandomInteger[];           (* 1 bit *)
    If[b1 == 0,
      mover = sorted[[1]],          (* rank-1: 1 bit consumed (prob 1/2) *)
      b2 = RandomInteger[];         (* 2nd bit *)
      If[b2 == 0,
        mover = sorted[[2]],        (* rank-2: 2 bits consumed (prob 1/4) *)
        (* b1=1, b2=1: use a 3rd bit to pick rank-3 either way *)
        b3 = RandomInteger[];       (* 3rd bit *)
        mover = sorted[[3]]         (* rank-3: 3 bits consumed (prob 1/4) *)
                                    (* N.B. both b3=0 and b3=1 -> rank-3  *)
      ]
    ];

    (* Direction: 0=left, 1=right, 2=up, 3=down (2 bits, exact) *)
    dir    = RandomInteger[{0, 3}];
    target = neighbors2Df[mover][[dir + 1]];

    (* Hard-core rejection: target occupied *)
    If[MemberQ[sorted, target],
      sorted,
      (* Metropolis acceptance *)
      others = DeleteCases[sorted, mover];
      newOcc = Sort[Append[others, target]];
      dE = energy[newOcc] - energy[sorted];
      If[RandomReal[] < MetropolisProb[dE], newOcc, sorted]
    ]
  ]

(* ---- Checker interface ---- *)
seedState = {1, 5, 9}
numBeta   = 1
