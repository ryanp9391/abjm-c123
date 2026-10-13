**Word classes incl. inverses and the 4-factor class (python, exact mod p, in the cyclic space; rank / dim loc / T1,T2 reachable):**

| word set | L = 2 | L = 3 |
|---|---|---|
| 1 + 48 positive singles | 11 / 3 / no | 17 / 3 / no |
| 1 + 48 inverse singles t_f⁻¹ | 49 / 3 / no | 49 / 3 / no |
| 1 + positive singles + inverse singles (97 words) | 54 / 6 / **yes** | 61 / 3 / **no** |
| … + (inverse f)(single g) (2401 words) | 54 / 6 / yes | **588 / 6 / yes** |
| 1 + singles + {t₁(θ₁)⁻¹, t₁(θ₂)⁻¹}×(1, single) (task inverses) | 24 / 4 / no | 42 / 4 / no |
| ≤ 2 positive + {t₁(θ₁)⁻¹, t₁(θ₂)⁻¹}×(≤ 2) | 51 / 6 / yes | 229 / 4 / no |
| ≤ 3 positive + {t₁(θ₁)⁻¹, t₁(θ₂)⁻¹}×(≤ 2) | 53 / 6 / yes | 400 / 4 / no |
| task classes (1225 + A·F·F + A·A·F = 15625 words); any one of A·F·F, A·A·F, A·A·A, F·F·F alone gives the same | 53 / 6 / yes | 400 / 4 / no |
| … + 504 F-site t₁t₂t₃ products | 53 / 6 / yes | 400 / 4 / no |
| … + ×t₁(θ₁), ×t₁(θ₁)⁻¹, both (4-factor ∧³-candidate class) | 54 / 6 / yes | 454, 454, 506 / 4 / **no** |
| ≤ 2 positive + all 48 inverses ×(1, single) (3577 words) | 54 / 6 / yes | **588 / 6 / yes** |
| inverses of t₁ only / t₂ only / t₃ only / A-site only / F-site only / m = 0 only, ×(1, single) | | 296 / 232 / 269 / 385 / 385 / 72, loc 3–4, no |

Caveat: t_f⁻¹ψ₀ is *not* local for generic f (e.g. t₁(θ₂)⁻¹ψ₀ ∝ ψ₀R₂₁⁻¹G₂⁻¹R₂₄⁻¹R₂₃⁻¹ moves the F-site vector onto site 4); only t₁(θ₁)⁻¹ is explicit. Inverse words must be treated exactly: a double-precision run of the L = 2 all-inverse span gave rank 50–52 instead of 54 and spurious non-membership.
