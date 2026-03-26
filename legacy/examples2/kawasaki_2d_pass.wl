(* ================================================================
   2D Kawasaki dynamics on a 3x3 torus -- PASSES detailed balance
   ================================================================

   System: 3x3 periodic lattice (torus), 3 particles, 6 holes.
   State:  sorted triple of occupied sites {s1,s2,s3}, sites 1..9.
           84 states total (C(9,3)).

   Site numbering (row-major):
     1 2 3
     4 5 6
     7 8 9

   Energy: sum of site energies (rational) + NN coupling for each
           adjacent pair of occupied sites (horizontal or vertical,
           periodic boundary).  No symbolic params -- energies are
           rational numbers, enabling fast exact simplification.

   Move:   pick a random particle (3 particles, rejection sampling
           via RandomInteger[{0,2}], 2 bits), pick a random direction
           (N/S/E/W via RandomInteger[{0,3}], 2 bits).
           If target site is occupied: hard-core rejection (stay put).
           Otherwise: Metropolis acceptance.

   The proposal is symmetric (each particle equally likely with
   rejection sampling, each direction equally likely) and Metropolis
   is correct, so detailed balance holds exactly.

   NOTE: With 84 states, this check is more computationally intensive
   than 1D examples.  Most state pairs are non-adjacent (T=0 trivially),
   so FullSimplify is only called for the ~300 reachable pairs.
   Expected runtime: a few minutes.
   ================================================================ *)

(* ---- Grid helpers (3x3 torus) ---- *)
$Lx = 3; $Ly = 3;

siteRow2D[s_Integer] := Ceiling[s / $Lx]
siteCol2D[s_Integer] := Mod[s - 1, $Lx] + 1
siteIdx2D[r_Integer, c_Integer] := (r - 1) * $Lx + c

leftOf2D[s_Integer]  := siteIdx2D[siteRow2D[s], Mod[siteCol2D[s] - 2, $Lx] + 1]
rightOf2D[s_Integer] := siteIdx2D[siteRow2D[s], Mod[siteCol2D[s],     $Lx] + 1]
upOf2D[s_Integer]    := siteIdx2D[Mod[siteRow2D[s] - 2, $Ly] + 1, siteCol2D[s]]
downOf2D[s_Integer]  := siteIdx2D[Mod[siteRow2D[s],     $Ly] + 1, siteCol2D[s]]
neighbors2D[s_Integer] := {leftOf2D[s], rightOf2D[s], upOf2D[s], downOf2D[s]}
adjacentQ2D[a_Integer, b_Integer] := MemberQ[neighbors2D[a], b]

(* ---- Energy (rational site energies + NN coupling) ---- *)
(* Site energies: corner=0, edge=1/2, centre=1 *)
$siteE = {0, 1/2, 0, 1/2, 1, 1/2, 0, 1/2, 0}
$nnJ   = 1/2

energy[occ_List] :=
  Total[$siteE[[#]] & /@ occ] +
  $nnJ * Length[Select[Subsets[occ, {2}],
                       adjacentQ2D[#[[1]], #[[2]]] &]]

(* ================================================================
   Algorithm: pick particle, pick direction, Metropolis.
   ================================================================ *)
Algorithm[occ_List] :=
  Module[{sorted, b, mover, dir, target, others, newOcc, dE},
    sorted = Sort[occ];
    (* Pick which particle: RandomInteger[{0,2}] = 0,1,2
       (2 bits; value 3 rejected by the checker's rejection sampling) *)
    b     = RandomInteger[{0, 2}];
    mover = sorted[[b + 1]];
    (* Pick direction: 0=left, 1=right, 2=up, 3=down (2 bits, exact) *)
    dir    = RandomInteger[{0, 3}];
    target = neighbors2D[mover][[dir + 1]];
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

(* BitsToState: bit string -> sorted-occupancy state.
   Accepts exactly ($Lx * $Ly)-bit strings with exactly 3 occupied sites.
   The bit at position i is 1 if site i is occupied (row-major, 1-indexed).
   Example: {1,0,0,0,1,0,0,0,1} -> {1,5,9}
   Returns None for any other length or particle count. *)
BitsToState[bits_List] :=
  If[Length[bits] =!= $Lx * $Ly || Total[bits] =!= 3, None,
    Flatten[Position[bits, 1]]]

numBeta = 1

(* ================================================================
   PlotState: visualise a state as an ASCII grid (for use in a
   Mathematica notebook or Print[] call).

   Example: Print[PlotState[{1, 5, 9}]]
   ================================================================ *)
PlotState[occ_List, label_String : ""] :=
  Module[{rows},
    rows = Table[
      StringJoin[Table[
        If[MemberQ[occ, (r - 1) * $Lx + c], " @ ", " . "],
        {c, 1, $Lx}]],
      {r, 1, $Ly}];
    If[label != "", Print[label]];
    Scan[Print, rows];
    Print["  (@ = particle, . = empty)"]
  ]

(* ================================================================
   PlotTransition: show two states side by side with an arrow.
   Example: PlotTransition[{1,5,9}, {2,5,9}]
   ================================================================ *)
PlotTransition[from_List, to_List] :=
  Module[{makeRow, rowsFrom, rowsTo},
    makeRow[occ_List, r_Integer] :=
      StringJoin[Table[
        If[MemberQ[occ, (r - 1) * $Lx + c], " @ ", " . "],
        {c, 1, $Lx}]];
    rowsFrom = Table[makeRow[from, r], {r, 1, $Ly}];
    rowsTo   = Table[makeRow[to,   r], {r, 1, $Ly}];
    Do[
      Print[rowsFrom[[r]], If[r == 2, "  -->  ", "       "], rowsTo[[r]]],
      {r, 1, $Ly}]
  ]
