(* One-pass action of t_a(u) on a batch of bras (columns of psi, modular integers stored as machine reals).
   psi . t_a(u): Phi[i,k] = psi delta_ik, then for each site Phi'[i,l] = sum_{k,a} Phi[i,k][..a..] Lax[k,l,a,b]; finally tr with the diagonal twist. *)
lcSite[phi_, s_, n_, d_, B_, mats_] := Module[{nl = 4^(s - 1), nr = 4^(n - s), x},
  x = Transpose[ArrayReshape[phi, {d, d, nl, 4, nr, B}], {1, 5, 2, 6, 3, 4}];   (* -> (i, nl, nr, B, k, a) *)
  x = Mod[ArrayReshape[x, {d nl nr B, 4 d}] . mats, lcP];                        (* exact: entries < 2^22, <= 24 terms *)
  x = Transpose[ArrayReshape[x, {d, nl, nr, B, d, 4}], {1, 3, 5, 6, 2, 4}];     (* -> (i, l, nl, b, nr, B) *)
  ArrayReshape[x, {d, d, nl 4 nr, B}]];
lcPass[which_, u_, ths_, psi_] := Module[{n = Length[ths], d = If[which == 2, 6, 4], B = Dimensions[psi][[2]], phi, types = lcTypes[Length[ths]], tw, out},
  phi = ConstantArray[0., {d, d, 4^n, B}];
  Do[phi[[i, i]] = psi, {i, d}];
  Do[phi = lcSite[phi, s, n, d, B, lcLaxMat[which, types[[s]], u, ths[[s]]]], {s, n}];
  tw = lcMod /@ lcTwistDiag[which];
  out = Mod[Sum[tw[[i]] phi[[i, i]], {i, d}], lcP];
  Mod[lcMod[lcGprod[which, u, ths]] out, lcP]];
(* helpers: Kronecker product of site vectors, exact -> modular *)
lcKron[vs_] := Fold[Flatten[Outer[Times, #1, #2]] &, {1}, vs];
lcModV[v_] := N[Map[lcMod, v]];
lcSiteV[a_] := {{1, 0, 0, a}, {-a, 0, 0, 1}};            (* W_a pair: A site v, F site u *)
lcPsi0[nsites_, a_: lcA0] := lcKron[Flatten[Table[lcSiteV[a], {nsites/2}], 1]];
(* sanity: all t_a(u) commute (modular, exact), 4 and 6 sites *)
lcSanity[ths_] := Module[{n = Length[ths], psi, r, res},
  SeedRandom[11]; psi = Transpose[{N[RandomInteger[{0, lcP - 1}, 4^n]]}];
  res = Table[{a, b, lcPass[a, 7/5, ths, lcPass[b, -3/7, ths, psi]] == lcPass[b, -3/7, ths, lcPass[a, 7/5, ths, psi]]}, {a, 3}, {b, a, 3}];
  Flatten[res, 1]];
Grid[Prepend[Join[Map[Prepend[#, 4] &, lcSanity[Take[lcThStd, 4]]], Map[Prepend[#, 6] &, lcSanity[lcThStd]]], {"sites", "a", "b", "[t_a(7/5), t_b(-3/7)] = 0 mod p"}], Frame -> All]
