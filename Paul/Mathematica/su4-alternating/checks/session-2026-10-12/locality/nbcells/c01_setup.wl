(* Setup: parameters, Lax tensors for t1, t2, t3 (exact rationals), modular reduction.  All symbols prefixed lc. *)
lcP = 4194301;                      (* prime, p < 2^22 so that all products/sums below stay exact in machine doubles *)
lcH = 1/2;
lcThStd = {2/7, -3/5, 8/11, -11/8, 5/13, -7/9};
lcThB = {1/3, -2/5, 3/4, -5/4, 2/9, -7/10};          (* task's second point; NOT generic: theta3 - theta4 = 4h *)
lcThC = {19/23, 34/29, 41/40, -19/12, -7/11, 13/17}; (* generic third point *)
lcZ = {3/2, 5/4, 2/3, 4/7};
lcA0 = 1/3;
lcMod[q_] := Mod[Numerator[q] PowerMod[Denominator[q], -1, lcP], lcP];
lcEm[i_, j_] := Normal[SparseArray[{{i, j} -> 1}, {4, 4}]];
lcRho = <|"A" -> Table[lcEm[k, l], {k, 4}, {l, 4}], "F" -> Table[-lcEm[l, k], {k, 4}, {l, 4}]|>;
lcSubs2 = Subsets[Range[4], {2}];
lcWC[k_, m_, Ir_] := Which[{k, m} == Ir, 1, {m, k} == Ir, -1, True, 0];
lcDw[X_] := Table[Sum[X[[k, Jc[[1]]]] lcWC[k, Jc[[2]], Ir] + X[[k, Jc[[2]]]] lcWC[Jc[[1]], k, Ir], {k, 4}], {Ir, lcSubs2}, {Jc, lcSubs2}];
lcDwE = Table[lcDw[lcEm[a, b]], {a, 4}, {b, 4}];           (* D(e_ab)_{IJ} on wedge^2 C^4, index [[a,b,I,J]] *)
(* Lax tensors Lax[[k,l]] = 4x4 matrix on the site (aux indices k,l) *)
lcLaxEx[1, typ_, u_, th_] := Table[(u - th) KroneckerDelta[k, l] IdentityMatrix[4] + lcH lcRho[typ][[k, l]], {k, 4}, {l, 4}];
lcLaxEx[3, "A", u_, th_] := Table[(u - th + lcH) KroneckerDelta[k, l] IdentityMatrix[4] - lcH lcEm[l, k], {k, 4}, {l, 4}];
lcLaxEx[3, "F", u_, th_] := Table[(u - th - 3 lcH) KroneckerDelta[k, l] IdentityMatrix[4] + lcH lcEm[k, l], {k, 4}, {l, 4}];
lcLaxEx[2, typ_, u_, th_] := Table[(u - th - If[typ === "F", lcH, 0]) KroneckerDelta[ii, jj] IdentityMatrix[4] + lcH Sum[lcDwE[[a, b, ii, jj]] lcRho[typ][[a, b]], {a, 4}, {b, 4}], {ii, 6}, {jj, 6}];
lcTypes[n_] := Table[If[OddQ[s], "A", "F"], {s, n}];
lcGs[1, _, _, _] := 1;
lcGs[2, "A", u_, th_] := u - th - lcH;
lcGs[2, "F", u_, th_] := u - th;
lcGs[3, "A", u_, th_] := (u - th - lcH) (u - th - 2 lcH);
lcGs[3, "F", u_, th_] := (u - th) (u - th - lcH);
lcGprod[which_, u_, ths_] := Times @@ MapThread[lcGs[which, #1, u, #2] &, {lcTypes[Length[ths]], ths}];
(* diagonal twist Lambda = diag(z): t1 uses Lambda, t3 uses det(Lambda) Lambda^{-T}, t2 uses wedge^2 Lambda *)
lcTwistDiag[1] := lcZ;
lcTwistDiag[3] := (Times @@ lcZ)/lcZ;
lcTwistDiag[2] := (Times @@ lcZ[[#]]) & /@ lcSubs2;
(* Lax as a modular (aux,site)x(aux,site) matrix of machine reals: row (k,a), column (l,b) *)
lcLaxMat[which_, typ_, u_, th_] := Module[{ex = lcLaxEx[which, typ, u, th], d = If[which == 2, 6, 4]},
  N[Mod[ArrayReshape[Transpose[Map[lcMod, ex, {4}], {1, 3, 2, 4}], {4 d, 4 d}], lcP]]];
Grid[{{"lcP", lcP}, {"h", lcH}, {"zeros of g (dropped factors) at theta_A, theta_F:", "t2(th_A+h), t2(th_F), t3(th_A+h), t3(th_A+2h), t3(th_F), t3(th_F+h)"}}, Frame -> All, Alignment -> Left]
