## Exact identities: what is local and universal (controls T3 and the F-site step)

Exact rational check (python `Fraction`, 4 **and 6** sites, θ standard, θ′, θ″; the Wolfram cell checks the same identities mod p). With θ₁ₛ = θ₁−θ_s, ρ_s(d) = multiply the site-s components by d_i:

* ψ₀ t₁(θ₁)⁻¹ = c₁ ψ₀ρ₁(Λ⁻¹), 1/c₁ = h ∏_{F s} θ₁ₛ ∏_{A s≥2}(θ₁ₛ + h);
* ψ₀ t₃(θ₁−h) = c₃ ψ₀ρ₁(Λ⁻¹), c₃ = −(3/4) det Λ ∏_{F s}(θ₁ₛ−h)(θ₁ₛ−2h)(θ₁ₛ−4h) ∏_{A s≥3}(θ₁ₛ−h)(θ₁ₛ−2h)(θ₁ₛ−3h) (h = 1/2; −3/4 = −6h³; closed form fitted from the symbolic 4-site result, verified exactly at 5 (L,θ) points);
* ψ₀ t₁(θ₁)⁻¹ t₁(θ₂) = c₁₂ ψ₀ρ₁(Λ⁻¹)ρ₂(Λ), c₁₂ = ∏_{s≥3}[A: θ₂ₛ/(θ₁ₛ+h); F: (θ₂ₛ−h)/θ₁ₛ];
* ψ₀ t₁(θ₂) = α[ψ₀ρ₂(Λ) + (σ/α) S⊗w₀], S = Σ_i e_i⊗e_i (pair-1 singlet), w₀ = W_a^{⊗(L−1)}; at 4 sites α = hθ₁₂θ₂₃(θ₂₄−h) (times θ₂₅(θ₂₆−h) per extra pair, checked), and **σ/α = h a (z₁−z₄)/θ₁₂ is independent of L**: the F-site step always drags the singlet.

Hence T3a = ψ₀ − (κ/c₃)ψ₀t₃(θ₁−h) (two words), T3b = (α′/c₃)·ψ₀t₃(θ₁−h) + β′(c₁/(c₃c₁₂))·ψ₀t₃(θ₁−h)t₁(θ₂); each added pair multiplies c₁, c₃, c₁₂ by explicit θ-factors — *these coefficients are local/universal*. The 8 individually localised words (≤ 3 factors) are the same at L = 2 and L = 3: 1, t₃(θ₁−h), t₁(θ₂), t₁(θ₁)t₃(θ₁−h), t₃(θ₁−h)t₁(θ₂), t₁(θ₂)t₃(θ₂+3h), t₃(θ₁−h)t₁(θ₂)t₃(θ₂+3h), t₁(θ₁)t₃(θ₁−h)t₁(θ₂); their pair-1 parts span exactly 4 dimensions: ψ₀, ρ₁(Λ⁻¹)ψ₀, ρ₁(Λ⁻¹)ρ₂(Λ)ψ₀ and ρ₂(Λ)ψ₀ + (σ/α)S⊗w₀. The 10 two-word relations among words (parallel bras, e.g. t₁(θ₁−h) ∥ t₂(θ₁)t₃(θ₁−h), t₂(θ₁) ∥ t₁(θ₁−h)t₁(θ₁), t₃(θ₁) ∥ t₁(θ₁)t₂(θ₁−h), t₂(θ₁−h) ∥ t₃(θ₁−h)t₃(θ₁), t₁(θ₁−h)t₃(θ₁) ∥ t₂(θ₁−h)t₂(θ₁) and the five analogues at θ₂+h, θ₂+2h) are fusion relations and hold identically at L = 2 and L = 3.

**Expected output of the cell below** (all True): rows for 4 and 6 sites at the three θ points.
