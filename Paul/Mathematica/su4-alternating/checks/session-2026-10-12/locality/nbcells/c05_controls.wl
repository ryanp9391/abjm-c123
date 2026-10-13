(* Exact (mod p) identities for the reachable controls and the structure of the F-site step, at 4 and 6 sites.
   rho_s(diag d) multiplies the component with site-s index i by d_i. *)
lcSiteDiag[site_, n_, dvec_] := Flatten[Fold[Outer[Times, #1, #2] &, {1}, Table[If[s == site, dvec, {1, 1, 1, 1}], {s, n}]]];
lcScaleSite[psi_, site_, n_, dvec_] := Mod[psi * (lcMod /@ lcSiteDiag[site, n, dvec]), lcP];     (* psi: flat modular vector *)
lcC1[ths_] := 1/(lcH Product[If[OddQ[s], ths[[1]] - ths[[s]] + lcH, ths[[1]] - ths[[s]]], {s, 2, Length[ths]}]);
lcC3[ths_] := -(3/4) (Times @@ lcZ) Product[With[{d = ths[[1]] - ths[[s]]}, If[OddQ[s], (d - lcH) (d - 2 lcH) (d - 3 lcH), (d - lcH) (d - 2 lcH) (d - 4 lcH)]], {s, 2, Length[ths]}];
lcC12[ths_] := Product[If[OddQ[s], (ths[[2]] - ths[[s]])/(ths[[1]] - ths[[s]] + lcH), (ths[[2]] - ths[[s]] - lcH)/(ths[[1]] - ths[[s]])], {s, 3, Length[ths]}];
lcCtrl[ths_] := Module[{n = Length[ths], psi0, flat, col, p1, r1, r2, r3, x, w0, alpha, sigma, S, vu, fit},
  psi0 = Round[lcModV[lcPsi0[n]]]; col = Transpose[{psi0}]; w0 = lcW0Mod[n];
  p1 = lcScaleSite[psi0, 1, n, 1/lcZ];                                   (* psi0 rho_1(Lambda^-1) *)
  (* (1) calibration  psi0 = c1 (psi0 rho_1(L^-1)) t1(th1) *)
  r1 = Flatten[lcPass[1, ths[[1]], ths, Transpose[{Mod[lcMod[lcC1[ths]] p1, lcP]}]]] == psi0;
  (* (2) (psi0 t1(th1)^-1) t1(th2) = c12 psi0 rho_1(L^-1) rho_2(L) *)
  r2 = Flatten[lcPass[1, ths[[2]], ths, Transpose[{Mod[lcMod[lcC1[ths]] p1, lcP]}]]] == Mod[lcMod[lcC12[ths]] lcScaleSite[p1, 2, n, lcZ], lcP];
  (* (3) psi0 t3(th1-h) = c3 psi0 rho_1(L^-1) *)
  r3 = Flatten[lcPass[3, ths[[1]] - lcH, ths, col]] == Mod[lcMod[lcC3[ths]] p1, lcP];
  (* (4) psi0 t1(th2) = alpha psi0 rho_2(L) + sigma S (x) w0  (S = sum_i e_i (x) e_i on pair 1) *)
  x = Mod[Mod[ArrayReshape[Flatten[lcPass[1, ths[[2]], ths, col]], {16, Length[w0]}] . w0, lcP] PowerMod[Mod[w0 . w0, lcP], -1, lcP], lcP];
  alpha = Mod[x[[4]] PowerMod[lcMod[lcZ[[4]]], -1, lcP], lcP];          (* (A=1,F=4) entry of v (x) u Lam is z4 *)
  sigma = x[[6]];                                                          (* (A=2,F=2) entry: pure singlet *)
  vu = Flatten[KroneckerProduct[lcModV[lcSiteV[lcA0][[1]]], lcModV[lcSiteV[lcA0][[2]] lcZ]]];
  S = Flatten[IdentityMatrix[4]];
  fit = Mod[alpha vu + sigma S, lcP] == x;
  {n, r1, r2, r3, fit}];
Grid[Prepend[Join[Table[lcCtrl[Take[ths, 4]], {ths, {lcThStd, lcThB, lcThC}}], Table[lcCtrl[ths], {ths, {lcThStd, lcThB, lcThC}}]],
  {"sites", "psi0 = c1 psi0 rho1(L^-1) t1(th1)", "t1(th1)^-1 t1(th2) -> c12 rho1(L^-1)rho2(L)", "t3(th1-h) -> c3 rho1(L^-1)", "t1(th2) psi0 = alpha rho2(L) psi0 + sigma S (x) w0"}], Frame -> All, Alignment -> Left]
