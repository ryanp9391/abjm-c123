(* Full cyclic space of psi0 under the commutative algebra of ALL t_a(u) (any shifts, any powers, inverses included, since the
   restriction of an invertible operator to a finite-dimensional invariant subspace has an inverse that is a polynomial in it):
   Krylov space of one generic element M = t1(u1) + c2 t2(u2) + c3 t3(u3) of the algebra, modular exact. *)
lcKrylov[ths_, nvec_, seed_: 1] := Module[{n = Length[ths], psi0, u1, u2, u3, c2, c3, x, vecs},
  SeedRandom[seed];
  {u1, u2, u3} = Table[RandomInteger[{-97, 97}]/RandomInteger[{11, 53}], {3}];
  {c2, c3} = RandomInteger[{2, lcP - 2}, 2];
  psi0 = Round[lcModV[lcPsi0[n]]]; x = Transpose[{psi0}];
  vecs = Table[With[{cur = x}, x = Mod[lcPass[1, u1, ths, cur] + Mod[c2 lcPass[2, u2, ths, cur], lcP] + Mod[c3 lcPass[3, u3, ths, cur], lcP], lcP]; Flatten[cur]], {nvec}];
  lcRR[vecs]];
lcKrylovReport[ths_, nvec_, seed_: 1] := Module[{n = Length[ths], bas = lcKrylov[ths, nvec, seed], T = lcTargetsMod[Length[ths]], w0 = lcW0Mod[Length[ths]]},
  {n, nvec, Length[bas], lcLocDim[bas, w0], And @@ Map[lcIn[bas, #] &, Values[KeySelect[T, #[[1]] === "T1" &]]],
   And @@ Map[lcIn[bas, #] &, Values[KeySelect[T, #[[1]] === "T2" &]]], And @@ Map[lcIn[bas, #] &, Values[KeySelect[T, #[[1]] === "T3a" &]]]}];
lcKry2 = lcKrylovReport[Take[lcThStd, 4], 70];
lcKry3 = lcKrylovReport[lcThStd, 650];
Grid[{{"sites", "# Krylov vectors generated", "cyclic dimension (rank)", "dim localised on pair 1", "all T1", "all T2", "all T3a"}, lcKry2, lcKry3}, Frame -> All, Alignment -> Left]
