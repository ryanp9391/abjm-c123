(* T1 and T2 lie in  span{psi0, psi0 t3(th1-h), psi0 t3(th1-h) t1(th2), psi0 t1(th2)}  (+)  C (S (x) W_a^{L-1}),  S = sum_i e_i (x) e_i on pair 1.
   The four words are the individually localised 'SoV' words; the whole obstruction is the isolated pair singlet. *)
lcFourPlusS[ths_] := Module[{n = Length[ths], psi0, col, w3, w31, w1b, S, T, B4, B5, t12},
  psi0 = Round[lcModV[lcPsi0[n]]]; col = Transpose[{psi0}];
  w3 = Round[Flatten[lcPass[3, ths[[1]] - lcH, ths, col]]];
  w31 = Round[Flatten[lcPass[1, ths[[2]], ths, Transpose[{w3}]]]];
  w1b = Round[Flatten[lcPass[1, ths[[2]], ths, col]]];
  S = Round[lcModV[lcKron[Join[{Flatten[IdentityMatrix[4]]}, Flatten[Table[lcSiteV[lcA0], {n/2 - 1}], 1]]]]];
  T = Select[lcTargetsMod[n], True &]; t12 = KeySelect[T, MemberQ[{"T1", "T2"}, #[[1]]] &];
  B4 = {psi0, w3, w31, w1b}; B5 = Append[B4, S];
  {n, lcRank[B4], lcRank[B5], lcIn[B4, S], And @@ Map[lcIn[B4, #] &, Values[t12]], And @@ Map[lcIn[B5, #] &, Values[t12]]}];
Grid[Prepend[Flatten[Table[Prepend[Rest[lcFourPlusS[Take[ths, n]]], n], {ths, {lcThStd, lcThB, lcThC}}, {n, {4, 6}}], 1],
  {"sites", "rank 4 words", "rank 4 words + S(x)w0", "S(x)w0 in span(4 words)", "all T1,T2 in span(4 words)", "all T1,T2 in span(4 words + S(x)w0)"}], Frame -> All, Alignment -> Left]
