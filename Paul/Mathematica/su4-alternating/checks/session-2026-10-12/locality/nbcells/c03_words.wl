(* Word machinery.  A factor is t_a(theta_s + m h), s = 1 (A site of pair 1) or 2 (F site of pair 1), m = -4..4, a = 1,2,3.
   Identically vanishing factors (a zero of the site prefactor g) are dropped: 54 -> 48 factors. *)
lcFactors[ths_] := Select[Flatten[Table[<|"a" -> a, "s" -> s, "m" -> m, "u" -> ths[[s]] + m lcH|>, {s, 2}, {a, 3}, {m, -4, 4}], 2],
  lcGprod[#["a"], #["u"], ths] =!= 0 &];
lcFname[f_] := "t" <> ToString[f["a"]] <> "(th" <> ToString[f["s"]] <> If[f["m"] >= 0, "+", ""] <> ToString[f["m"]] <> "h)";
(* modular reduced row-echelon basis of a set of rows; rows are exact integers mod p *)
lcRR[rows_] := DeleteCases[RowReduce[Round[rows], Modulus -> lcP], {0 ..}];
lcRank[rows_] := MatrixRank[Round[rows], Modulus -> lcP];
lcIn[bas_, x_] := MatrixRank[Append[Round[bas], Round[x]], Modulus -> lcP] == Length[bas];
(* dimension of  span(bas) intersect (C^16 (x) w0):  rank-nullity for  psi -> psi (1 - w0 w0^T/(w0.w0)) on the (pair1) x (rest) reshape *)
lcLocDim[bas_, w0_] := Module[{R = Length[w0], P, proj},
  P = Mod[IdentityMatrix[R] - Outer[Times, w0, w0] PowerMod[Mod[w0 . w0, lcP], -1, lcP], lcP];
  proj = Map[Flatten[Mod[ArrayReshape[#, {16, R}] . P, lcP]] &, Round[bas]];
  Length[bas] - lcRank[proj]];
(* targets, exactly (rational), then modular.  rest = W_a^{L-1} on the other pairs *)
lcTargets[nsites_, a_: lcA0] := Module[{rest = Flatten[Table[lcSiteV[a], {nsites/2 - 1}], 1], v = lcSiteV[a][[1]], u = lcSiteV[a][[2]], T = <||>},
  Do[T[{"T1", a1}] = lcKron[Join[lcSiteV[a1], rest]], {a1, {-2/3, 0, 2/3, 5/4}}];
  Do[T[{"T2", k}] = lcKron[Join[{v (1 - k lcZ), u/(1 - k lcZ)}, rest]], {k, {1/3, -2/5, 5/3}}];
  Do[T[{"T3a", k}] = lcKron[Join[{v (1 - k/lcZ), u}, rest]], {k, {1/3, -2/5, 5/3}}];
  T[{"T3b", "vL^-1 x u(1+L)"}] = lcKron[Join[{v/lcZ, u (1 + lcZ)}, rest]];
  T[{"T3b", "vL^-1 x u(2-L)"}] = lcKron[Join[{v/lcZ, u (2 - lcZ)}, rest]];
  T];
lcTargetsMod[nsites_] := Map[lcModV, lcTargets[nsites]];
lcW0Mod[nsites_] := Round[lcModV[lcKron[Flatten[Table[lcSiteV[lcA0], {nsites/2 - 1}], 1]]]];
(* layered span: V_0 = <psi0>, V_{l+1} = V_l + sum_f V_l . t_f  (f runs over the 48 pair-1 factors, positive powers) *)
lcLayers[ths_, maxl_, chunk_: 40] := Module[{n = Length[ths], F = lcFactors[ths], psi0, bas, prev, T, w0, res = {}, cand, newDim, t0 = AbsoluteTime[]},
  psi0 = lcModV[lcPsi0[n]]; T = lcTargetsMod[n]; w0 = lcW0Mod[n];
  bas = {Round[psi0]};
  AppendTo[res, {0, 1, lcLocDim[bas, w0], Map[lcIn[bas, #] &, T], 0.}];
  Do[
   prev = bas;
   Do[
    cand = Flatten[Table[Round[Transpose[lcPass[f["a"], f["u"], ths, Transpose[blk]]]], {blk, Partition[prev, UpTo[chunk]]}], 1];
    bas = lcRR[Join[bas, cand]],               (* reduce factor by factor to keep the matrices small *)
    {f, F}];
   newDim = Length[bas] - Length[prev];
   AppendTo[res, {l, Length[bas], lcLocDim[bas, w0], Map[lcIn[bas, #] &, T], AbsoluteTime[] - t0}];
   If[newDim == 0, Break[]],
   {l, 1, maxl}];
  res];
lcLayerGrid[res_] := Grid[Prepend[Map[{#[[1]], #[[2]], #[[3]], And @@ KeyValueMap[If[#1[[1]] === "T1", #2, True] &, #[[4]]], And @@ KeyValueMap[If[#1[[1]] === "T2", #2, True] &, #[[4]]], And @@ KeyValueMap[If[#1[[1]] === "T3a", #2, True] &, #[[4]]]} &, res],
   {"layer l (<= l factors)", "rank of span", "dim localised on pair 1", "all T1 in span", "all T2 in span", "all T3a in span"}], Frame -> All, Alignment -> Left];
