# G₂/SU(3) = S⁶, its octonion almost complex structure, and the finite-subgroup limit

Three exact integer scripts measured the geometry, and an independent second seat attacked each claim. Its verdicts are in the Status column below; its receipt is [`refute/REFUTE.hbp`](refute/REFUTE.hbp), and the restatements sit beside the claims in [`CORRECTIONS.hbp`](CORRECTIONS.hbp).

| | Claim | Status | Measured or named |
|---|---|---|---|
| C1 | G₂ = Aut(O) acts transitively on the unit imaginary octonions, with stabilizer SU(3), so S⁶ = G₂/SU(3) | supported: measured for the algebra, named for the group | dim der(O) = 14; the stabilizer of e₇ has dimension 8 and is su(3); the orbit map has rank 6 |
| C2 | J_x(v) = x × v is a G₂-invariant almost complex structure, and the tangent space is the 3 of SU(3) | supported, measured | J² = −1, J orthogonal, commuting with su(3); the commutant is span(I, J), so the only invariant almost complex structures are ±J. Which of the two is called 3 is a convention |
| C3 | J is not integrable | supported, measured | the Nijenhuis tensor is nonzero on 24 of 36 basis pairs, measured two independent ways; the Koszul test agrees, with [m¹⁰, m¹⁰] landing purely in m⁰¹ through ε. CP² = SU(3)/U(2) is integrable as a control. Reason restated in [G3](CORRECTIONS.hbp) |
| C4 | (S⁶, round metric, J) is strict nearly Kähler; its cone is flat R⁷ with the G₂ 3-form | supported: measured pointwise, the cone named | (∇_X J)X = 0 with ∇J ≠ 0 at two points, in two octonion tables |
| C5 | W(G₂) = D₆ = S₃ × Z₂, with S₃ = W(SU(3)) from the long roots | supported, measured | \|W\| = 12; W(A₂) has order 6 and does not contain −I; the short roots form two orbits of 3 (the weights of 3 and 3̄); an element of G₂ sends e₇ to −e₇ and J to −J |
| C6 | A Higgs field in the 7 breaks G₂ to SU(3): 8 massless and 6 mass-degenerate vectors | supported, measured (classical, tree level) | 3⟨Xv, Yv⟩ = −tr(XY) on the six-dimensional complement; Casimir −2 on R⁷ |
| C7 | Λ²R⁷ = Λ²₇ ⊕ Λ²₁₄ with Λ²₁₄ = g₂; one chirality of 2-forms on a coassociative R⁴ lies in Λ²₁₄, so R³-invariant pullbacks of 4D instantons of that chirality are G₂-instantons | supported, measured | eigenvalue 2 has dimension 7 and −1 has dimension 14; Λ²₁₄ = der(O); checked on all 7 coassociative planes |
| C8 | Limits of finite subgroups of SU(N) lie at most in N(T) = T ⋊ W | **overstated**, restated in [G4](CORRECTIONS.hbp) | every limit is torus-by-finite, but it need not lie in N(T): Z_m·2I converges to U(1)·2I in SU(3) |
| C9 | The finite repo's S₃ is what maximal-abelian projection leaves of SU(3) | **overstated**, restated in [G5](CORRECTIONS.hbp) | S₃ = Δ(6) ≅ W(SU(3)), and it meets the torus trivially. The finite repo's 8 flat sectors equal the 8 irreps of the quantum double D(S₃), measured. Transfer to continuum SU(N) is OPEN |
| C10 | Routes from finite data to the Clay target | **overstated**, restated in [G6](CORRECTIONS.hbp) | the restatement lists exactly what the Clay problem requires; a gap at fixed physical size is a separate, well-posed target |

## Files

| Folder | Receipt | sha256 | Replay |
|---|---|---|---|
| [`m1/`](m1/) | `M1-OCTONION.hbp`: the octonions, der(O), su(3), J, Nijenhuis, Koszul, CP² control | `0024b12e…` | CI, byte for byte |
| [`m2/`](m2/) | `M2-ROOTS.hbp`: G₂ roots, W(G₂), the Higgs spectrum, the antipodal element | `0503594c…` | CI, byte for byte |
| [`m3/`](m3/) | `M3-TOWERS.hbp`: the Δ(6n²) tower, the Jordan argument, 2I, G₂-instantons | `86b840b1…` | CI, byte for byte |
| [`refute/`](refute/) | `REFUTE.hbp`: 71 checks, 10 controls, 0 failures | `fbd50476…` | evidence, kept verbatim ([`EVIDENCE-ONLY`](refute/EVIDENCE-ONLY)) |

In `m2/m2_roots.py`, the two powers x^dh and (x − p)^dm were rewritten with an integer power helper that refuses a negative exponent, for the gate. Its receipt is byte-identical to the original run.

The refutation folder is pinned by a folder row in [`../../tools/FLOAT-EVIDENCE.txt`](../../tools/FLOAT-EVIDENCE.txt), because its runtime float detector names floats in order to find them.

These are classical results: Cartan, Ehresmann–Libermann, Koszul–Frölicher, Gray, Bryant, Harvey–Lawson. This folder holds exact, reproducible certificates of them, with each claim's scope stated. Whether S⁶ carries any integrable complex structure at all, the Hopf problem, is open.

**Beside the sentence above, November 9, 2026, from a verified literature map.** In the refereed literature the Hopf problem is open. One existence claim is under active, unrefereed scrutiny:

- **The construction.** L. Alpöge, working with Claude, "A compact complex threefold fibred by tori over the projective line, and the six-sphere" (alpo.ge/s6.pdf, August 2026), builds a compact complex threefold fibred by 2-tori over P¹. It is simply connected with the integral homology of S⁶, hence diffeomorphic to S⁶, with algebraic dimension 1.
- **The second account.** P. Engel, "Complex structures on S⁶" ([arXiv:2609.38442](https://arxiv.org/abs/2609.38442), September 2026), recasts the construction as a self-contained account.
- **Expert view.** R. Bryant, quoted in *Scientific American* on September 24, 2026: "There seems to be an emerging consensus that the construction is plausible."
- **Formalization.** A Lean 4 formalization is public (github.com/plby/HopfProblem). It has not been built here (NOT MEASURED).
- **The published theorem it contradicts.** Campana–Demailly–Peternell, *Compositio Math.* 156 (2020), Cor. 2.3, says a complex S⁶ has algebraic dimension 0. The manuscript's §10 names the step of that proof it says fails, at a non-normal fibre.

The results this folder measures hold however the claim resolves. The G₂-invariant J is the unique invariant almost complex structure, and it is not integrable. No complex structure is orthogonal to the round metric (LeBrun 1987).
