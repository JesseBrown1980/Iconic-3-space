# Transfer theorem: finite subgroups of SU(2) reproduce the SU(2) strong-coupling series through their design order

2026-11-09 · seat Opus 5.5 (pid `8467a937cba309f7`) · E=0: local file, no GitHub writes.

**Tags.** **DERIVED**: proved in this document (§3–§5). The proof is written out here and has not yet had a second reader. **MEASURED**: an exact number on disk from the named script (§6), quoted inline. **NAMED**: a cited input, stated and not re-derived here. **CONJECTURE**: unproved. **UNVERIFIED**: a citation or content I could not open.

---

## 0. Result in brief

On every finite lattice whose plaquettes are words in four distinct links (H1), with a mild sharing condition (H2′) that every hypercubic lattice with all side lengths at least 2 satisfies, the following holds. A finite subgroup G ⊂ SU(2) that is a spherical t-design has the same strong-coupling series as SU(2) through β^t, for log Z and for every observable of link-degree d through β^(t−d). The first deviation is an explicit positive integer times a universal constant:

  log Z_G − log Z_SU(2) = N_p · ι_G · β^(t+1) / (2^(t+1) (t+1)!) + O(β^(t+2)),  where ι_G = dim (V_{(t+1)/2})^G = M_{t+1}^G − M_{t+1}^Haar ∈ {1, 2, …}.

Each single plaquette shifts by ⟨f_p⟩_G − ⟨f_p⟩_SU(2) = ι_G β^t / (2^(t+1) t!) + O(β^(t+1)). **DERIVED**, §4. Of 106 predictions on lattices that satisfy the hypotheses, the exact enumeration checks confirmed all 106 (**MEASURED**, §6).

The ladder of subgroups (**MEASURED**, §2.4) runs Dic₃ (t = 3) → 2T (t = 5) → 2O (t = 7) → 2I (t = 11). The value t = 11 is the largest design strength of any finite subgroup of SU(2) (**DERIVED** from the classification, which is **NAMED**).

### Corrections to the statement as posed

Each correction comes with its evidence.

1. **(H1) must be stated.** Each plaquette word must contain each link at most once. "4-link cycle" covers this only if it means four distinct links. A periodic direction with L = 1 breaks it, and then claim (A) itself fails. On the 1×1 periodic torus, Dic₃ already differs from SU(2) at β², with Δlog Z₂ = 7/48 although t = 3. Q8, Dic₄, 2T and 2O first differ at orders 2, 2, 3 and 4 (**MEASURED**). Claims (A) and (C) need **only** H1.
2. **(H2) is sufficient; the condition the proof uses is (H2′):** for every link l and every plaquette p ∋ l, p has a link l′ ≠ l that lies in no other plaquette containing l. H2 implies H2′ (Lemma 1.1). The 2×2 periodic torus breaks H2, because neighbouring plaquettes share two parallel links. It satisfies H2′, and every claim holds on it exactly (**MEASURED**: Q8, Dic₃, Dic₄, 2T, 2O). **So the 2×2 torus is covered by the theorem.**
3. **Counterexamples showing that a sharing condition is needed for (B):**
   - The *pillow* (two plaquettes on the same four links) satisfies H1 and breaks H2′. There Δlog Z_(t+1) = ι/(t+1)!, which is 2^t times the formula: for Dic₃, 1/24 against 1/192 (**MEASURED** and DERIVED, §5.2).
   - The 2×1 periodic strip (one L = 1 direction) breaks H1, H2 and H2′. (A) still holds there, but (B) fails: for Dic₃, 65/2304 against 1/192 (**MEASURED**).
4. **(D) holds as stated, and it also holds for t₁ = t₂.** The formula N_p (M^G1 − M^G2)_(t+1) / (2^(t+1)(t+1)!) is valid then too, and it can be 0. Example: Dic₃ and Dic₄, both with t = 3 and M₄ = 3, agree at β⁴. On the cube surface and the 2×2 torus they first differ at β⁶ (**MEASURED**). Order 6 lies outside the theorem's scope.
5. **Sign and integrality of (B):** M_(t+1)^G − M_(t+1)^Haar = ι_G is a non-negative integer, and it is ≥ 1 exactly when t is the maximal design strength. The finite group's log Z therefore lies *above* SU(2)'s at the first differing order (Corollary 2.5).
6. **Design ⇔ moments:** for a *subgroup* G the two conditions are equivalent, and both directions are proved (Lemma 2.4). For a general point set on S³, only design ⇒ moments holds. The proof of the converse uses dim V_j^G ≥ 0, which needs the group structure.
7. **Which claims use the group structure:** (B), (C′) and the counting in (D) use invariance of the measure (Lemma 4.3). (A) and (C) use only the design property and the product measure, so they hold for any t-design point set X ⊂ SU(2) as well (Remark 8.1).

---

## 1. Setting and hypotheses

**1.1 The complex.** Λ is a finite oriented graph: a finite vertex set, and a finite link set 𝓛 in which each link l has a tail and a head (loops and parallel links are allowed). 𝓟 is a finite set of plaquettes, N_p = |𝓟|. Each plaquette p is a closed walk, given as a word

  w_p = (l₁^ε₁, l₂^ε₂, l₃^ε₃, l₄^ε₄),  εᵢ ∈ {+1, −1},

in which consecutive steps match: the head of a step is the tail of the next, and l^(−1) traverses l from head to tail. For a configuration U ∈ K^𝓛 put U_p = U_{l₁}^ε₁ U_{l₂}^ε₂ U_{l₃}^ε₃ U_{l₄}^ε₄ and f_p = tr(U_p)/2. In SU(2) the eigenvalues are e^(±iθ), so tr U^(−1) = tr U. Hence f_p is unchanged by cyclic rotation of the word and by reversal, f_p ∈ [−1, 1], and the action S = Σ_p f_p is a real function on K^𝓛.

**1.2 Hypotheses.**
- **(H1)** each plaquette word contains four distinct links (each link at most once).
- **(H2)** two distinct plaquettes share at most one link.
- **(H2′)** for every link l and every plaquette p ∋ l there is a link l′ ∈ p, l′ ≠ l, such that l′ ∉ q for every plaquette q ≠ p with l ∈ q. In words: inside the star of l, every plaquette has a link of its own.

**Lemma 1.1** (H1 + H2 ⇒ H2′). *DERIVED.* Let p ∋ l. By H1, p has three links other than l. By H2, every other q ∋ l shares only l with p. So all three are private to p inside the star of l. ∎

**Proposition 1.2** (hypercubic lattices). *DERIVED.* Let Λ be the periodic hypercubic lattice ℤ_{L₁} × … × ℤ_{L_D}, D ≥ 2. It has links μ(x): x → x + μ̂ and plaquettes P_{μν}(x) = μ(x) ν(x+μ̂) μ(x+ν̂)^(−1) ν(x)^(−1) for μ < ν. Alternatively, let Λ be an open box: links and plaquettes are those of ℤ^D lying inside the box.
1. Open box, any sizes: H1 and H2 hold.
2. Periodic with all L_μ ≥ 3: H1 and H2 hold.
3. Periodic with all L_μ ≥ 2: H1 and H2′ hold. If some L_ν = 2, H2 fails: P_{μν}(x) and P_{μν}(x+ν̂) share μ(x) and μ(x+ν̂).
4. Periodic with some L_ν = 1: H1 fails, because P_{μν}(x) = μ(x) ν(x+μ̂) μ(x)^(−1) ν(x)^(−1) contains μ(x) twice.

*Proof.*
- **H1.** The links of P_{μν}(x) are μ(x), μ(x+ν̂), ν(x), ν(x+μ̂). They are distinct iff ν̂ ≢ 0 and μ̂ ≢ 0, that is iff L_μ, L_ν ≥ 2. The links of different directions are always distinct.
- **H2.** Let P = P_{μν}(x) and P′ ≠ P share two links.
  - *Both shared links have the same direction μ.* Then P′ = P_{μν′}(y) and {x, x+ν̂} = {y, y+ν̂′}. If y = x, then ν̂ ≡ ν̂′, which forces ν = ν′ when L ≥ 2, and so P′ = P. Otherwise y = x+ν̂ and ν̂+ν̂′ ≡ 0. If ν′ = ν this needs L_ν | 2. If ν′ ≠ ν it needs L_ν = L_ν′ = 1.
  - *The shared links have directions μ and ν.* Then P′ = P_{μν}(y) with y ∈ {x ± ν̂} ∩ {x ± μ̂}, y ≠ x, which needs L_μ = L_ν = 1.

  So with all L ≥ 3 sharing two links is impossible, and with all L ≥ 2 it happens only through L_ν = 2. In the open box every congruence above is an equality in ℤ^D, and none of them can hold.
- **H2′ for all L ≥ 2.** The star of μ(x) is {P_{μσ}(x), P_{μσ}(x−σ̂) : σ ≠ μ}. P_{μσ}(y) has σ-links only at y and y+μ̂. Take p = P_{μν}(x), with ν-links ν(x) and ν(x+μ̂). Among the other plaquettes of the star, only P_{μν}(x−ν̂) carries ν-links, namely ν(x−ν̂) and ν(x−ν̂+μ̂). A coincidence would need ν̂ ≡ 0 or μ̂ ± ν̂ ≡ 0, that is L_ν = 1 or L_μ = L_ν = 1. So ν(x) is private to p in the star. The same computation applies to p = P_{μν}(x−ν̂), and by symmetry to every link. ∎

**MEASURED.** The script's hypothesis checker reports H1/H2/H2′ = 1/1/1 on the 3×3 torus, 1/0/1 on the 2×2 torus, 0/0/0 on the 2×1 strip and 0/1/0 on the 1×1 torus. These agree with Proposition 1.2 (`TRANSFER-CHECK.hbp`, rows `LAT`).

**1.3 Measures.** K is SU(2) with normalised Haar measure, or a finite subgroup G ⊂ SU(2) with the uniform measure. Write μ_K for either. Both are left- and right-invariant and invariant under g ↦ g^(−1). For G this is immediate, since multiplication and inversion are bijections. For Haar on a compact group it is NAMED (uniqueness and bi-invariance of Haar measure, e.g. Bröcker–tom Dieck, *Representations of Compact Lie Groups*, Ch. I). Configurations carry the product measure μ_K^⊗𝓛, and E_K denotes expectation.

**1.4 The series.** Z_K(β) = E_K[e^(βS)]. Since |S| ≤ N_p and Σ_n |β|^n N_p^n / n! < ∞, dominated convergence gives Z_K(β) = Σ_n z_n^K β^n, with z_n^K = E_K[S^n]/n!. So Z_K is entire, Z_K(0) = 1, and Z_K(β) > 0 for real β. The function log Z_K is then real-analytic on ℝ, and its Taylor series at 0 is the formal logarithm ℓ^K = log(1 + Σ_{n≥1} z_n β^n) (Lemma 4.4). Observables: ⟨O⟩_K(β) = E_K[O e^(βS)] / Z_K(β), with numerator coefficients a_n^K = E_K[O S^n]/n!.

## 2. Designs, moments and the ladder of subgroups

**Definition 2.1** (link-degree). Write a matrix U ∈ SU(2) as U = [[a, b], [−b̄, ā]]. A function F on K^𝓛 has *link-degree ≤ (d_l)_l* if it is a finite linear combination of products ∏_l m_l(U_l), where each m_l is a monomial of degree ≤ d_l in the four entries of U_l. For U ∈ SU(2), U^(−1) = [[U₂₂, −U₁₂], [−U₂₁, U₁₁]], and every complex-conjugated entry is ± an entry. The real coordinates (x₁, x₂, x₃, x₄), with a = x₁ + i x₂ and b = x₃ + i x₄, are related to the entries by an invertible linear map. So "degree ≤ d in the entries of U, of U^(−1), of Ū, or in the four real coordinates" is one and the same class.

**Definition 2.2** (design). A finite G ⊂ SU(2) is a *t-design* if E_G[m] = E_SU(2)[m] for every monomial m of degree ≤ t in the entries of U. Its *strength* t(G) is the largest such t.

**Lemma 2.3** (Definition 2.2 is the spherical t-design of Delsarte–Goethals–Seidel). Under U ↔ x ∈ S³ ⊂ ℝ⁴, Haar measure equals the normalised surface measure σ. Hence G is a t-design in the sense of Definition 2.2 iff the point set G ⊂ S³ is a spherical t-design in the sense of Delsarte, Goethals and Seidel (1977) (NAMED).

*Proof.* Left multiplication by g ∈ SU(2) is ℝ-linear on ℝ⁴. It preserves det U = |x|², so it is orthogonal. σ is O(4)-invariant (NAMED), so σ is a left-invariant probability measure on SU(2), and uniqueness of Haar measure (NAMED) gives σ = Haar. The DGS definition asks for equal averages of all polynomials of degree ≤ t in x, which by Definition 2.1 is the same condition. ∎

Notation for the rest of the document:
- D^j is the irreducible representation of SU(2) of dimension 2j + 1 on V_j ≅ Sym^(2j) ℂ², and χ_j = tr D^j. In particular χ_(1/2)(U) = tr U.
- m(n, j) is the multiplicity of V_j in V_(1/2)^⊗n. It satisfies m(0, 0) = 1, m(n+1, 0) = m(n, ½), and m(n+1, j) = m(n, j−½) + m(n, j+½) for j > 0. This is the Clebsch–Gordan rule V_j ⊗ V_(1/2) ≅ V_(j+1/2) ⊕ V_(j−1/2) (NAMED). Moreover m(n, n/2) = 1.
- M_n^K = E_K[(tr U)^n]. Then M_n^Haar = m(n, 0) = Catalan(n/2) for n even and 0 for n odd (NAMED). **MEASURED:** the recursion gives 2, 5, 14 and 132 at n = 4, 6, 8, 12 (`GROUP` rows, `haar_next`).

**Lemma 2.4** (equivalences for a finite subgroup G). *DERIVED.* For t ≥ 0 the following are equivalent:
- (i) G is a t-design.
- (ii) P_j := |G|^(−1) Σ_g D^j(g) = 0 for all j with 0 < 2j ≤ t.
- (iii) dim V_j^G = 0 for all j with 0 < 2j ≤ t.
- (iv) M_n^G = M_n^Haar for all 0 ≤ n ≤ t.

*Proof.*
- **(ii) ⇔ (iii).** D^j(h) P_j = P_j for every h ∈ G, so P_j² = P_j and the image of P_j consists of G-fixed vectors. Every fixed vector v satisfies P_j v = v. So P_j is a projection onto V_j^G, and P_j = 0 iff V_j^G = 0. Also tr P_j = dim V_j^G = E_G[χ_j].
- **(i) ⇔ (ii).** The entries of D^j(U) are homogeneous polynomials of degree 2j in the entries of U, because D^j is the symmetric power of the defining representation. Conversely, a monomial of degree k in the entries of U is an entry of U^⊗k. By complete reducibility (NAMED), U^⊗k = C (⊕_j D^j(U)^(⊕ m(k,j))) C^(−1) for a fixed invertible C. So that monomial is a linear combination of entries of the D^j(U) with 2j ≤ k and 2j ≡ k mod 2. Hence span{monomials of degree ≤ t} = span{entries of D^j : 2j ≤ t}. Next, E_Haar[D^j] = 0 for j > 0. By invariance it is the projection onto the SU(2)-fixed vectors of V_j, and an irreducible nontrivial V_j has none. Both averages of D⁰ = 1 equal 1. So (i) ⇔ E_G[D^j] = 0 for 0 < 2j ≤ t, which is (ii).
- **(iii) ⇔ (iv).** Expand χ_(1/2)^n = Σ_j m(n, j) χ_j. Averaging gives M_n^G = Σ_j m(n, j) dim V_j^G and M_n^Haar = m(n, 0). Hence

  (∗)  M_n^G − M_n^Haar = Σ_{j>0} m(n, j) dim V_j^G,

  a sum of non-negative terms that contains the term dim V_(n/2)^G with coefficient 1. (iii) ⇒ (iv) follows from (∗) directly. For (iv) ⇒ (iii), take 0 < n ≤ t: the left side of (∗) is 0 and every term is non-negative, so every term vanishes, in particular dim V_(n/2)^G = 0. ∎

**Corollary 2.5** (the first moment defect is an integer). *DERIVED.* If G is a t-design, then M_(t+1)^G − M_(t+1)^Haar = dim V_((t+1)/2)^G =: ι_G ∈ ℤ_(≥0). Moreover ι_G ≥ 1 iff t = t(G). *Proof:* in (∗) with n = t+1, the terms with 0 < 2j ≤ t vanish, and the only remaining j > 0 is 2j = t+1, with multiplicity 1. ∎

**Corollary 2.6.** If −1 ∈ G then t(G) is odd. *Proof:* D^j(−1) = (−1)^(2j), so P_j = −P_j for half-integer j. ∎

**2.4 The ladder.** **MEASURED:** script `check/transfer_check.py`, rows `GROUP`. Each group is generated by exact closure of unit quaternions over ℚ(√d). Three criteria are computed independently:
- (a) every monomial of degree ≤ t in the four real coordinates has the S³ average (the DGS definition, checked directly);
- (b) trace moments;
- (c) dim V_j^G via the Chebyshev recursion χ_(j+1/2) = χ_(1/2) χ_j − χ_(j−1/2).

All three give the same t for every group.

| G | order | t(G) | ι_G = M^G_(t+1) − M^Haar_(t+1) | first log Z deviation per plaquette, ι/(2^(t+1)(t+1)!) | first plaquette deviation, ι/(2^(t+1) t!) |
|---|---|---|---|---|---|
| Z₆ | 6 | 1 | 1 (2 − 1) | 1/8 | 1/4 |
| Q8 = Dic₂ | 8 | 3 | 2 (4 − 2) | 1/192 | 1/48 |
| **Dic₃** | 12 | **3** | **1 (3 − 2)** | **1/384** | **1/96** |
| Dic₄ | 16 | 3 | 1 (3 − 2) | 1/384 | 1/96 |
| 2T | 24 | 5 | 1 (6 − 5) | 1/46080 | 1/7680 |
| 2O | 48 | 7 | 1 (15 − 14) | 1/10321920 | 1/1290240 |
| 2I | 120 | 11 | 1 (133 − 132) | 1/1961990553600 | 1/163499212800 |

The last column was MEASURED for Q8, Dic₃, Dic₄, 2T and 2O on lattices (§6). For Z₆ and 2I it is the value of the DERIVED formula.

*General families* (DERIVED). V_j^G corresponds to the G-invariant binary forms of degree 2j, using V_j ≅ V_j^*.
- **Cyclic** Z_m = ⟨diag(ζ, ζ̄)⟩ with ζ = e^(2πi/m), m ≥ 2: the monomials x and y scale by ζ^(±1) ≠ 1, while xy is invariant. So t = 1, with ι = 1 for m ≥ 3 and ι = 3 for m = 2.
- **Binary dihedral** Dic_n = ⟨a = diag(e^(iπ/n), e^(−iπ/n)), b = [[0, 1], [−1, 0]]⟩, n ≥ 2: −1 = a^n ∈ G. In degree 2, only xy satisfies a − b ≡ 0 mod 2n, and b maps it to −xy, so there is no invariant. In degree 4, x²y² is invariant, and for n = 2 so is x⁴ + y⁴. So t(Dic_n) = 3 for every n ≥ 2, with ι = 1 for n ≥ 3 and ι = 2 for n = 2. This agrees with the MEASURED rows for n = 2, 3, 4.

By the classification of finite subgroups of SU(2), every one is conjugate to a cyclic group, a binary dihedral group, 2T, 2O or 2I (NAMED; Klein 1884). Conjugation by V ∈ SU(2) preserves Haar measure and degree, so it preserves t. **Hence t(G) ≤ 11 for every finite subgroup of SU(2), with equality iff G is conjugate to 2I** (DERIVED from the NAMED classification, the two family computations above, and the MEASURED values for 2T, 2O and 2I).

**The Dic₃ rung and S₃.** The only element of order 2 in SU(2) is −1: U² = 1 with det U = 1 forces U = ±1, and −1 is central. The subgroup of SU(2) covering S₃ ≅ D₃ ⊂ SO(3) is Dic₃ = ⟨a, b | a⁶ = 1, b² = a³, bab^(−1) = a^(−1)⟩. Its quotient Dic₃/{±1} is generated by an element of order 3 and an involution inverting it, so Dic₃/{±1} ≅ S₃ (DERIVED). The SU(2) transfer therefore reaches the S₃ programme through Dic₃:
- t = 3, ι = 1;
- first free-energy deviation N_p/384 · β⁴;
- first plaquette deviation β³/96 (MEASURED in §6).

## 3. The theorem (corrected statement)

Fix a finite complex Λ as in §1.1 satisfying **(H1)**. Let G ⊂ SU(2) be a finite subgroup that is a t-design, t ≥ 0; t need not be maximal. Write ΔX = X^G − X^SU(2) for any quantity X.

- **(A)** *[H1]* z_n^G = z_n^SU(2) and ℓ_n^G = ℓ_n^SU(2) for 1 ≤ n ≤ t. So log Z_G and log Z_SU(2) agree as power series through β^t.
- **(B)** *[H1 + H2′]*

  Δℓ_(t+1) = Δz_(t+1) = N_p (M_(t+1)^G − M_(t+1)^Haar) / (2^(t+1) (t+1)!) = N_p ι_G / (2^(t+1) (t+1)!).

  In E[S^(t+1)] = Σ_{(p₁,…,p_(t+1)) ∈ 𝓟^(t+1)} E[f_p₁ ⋯ f_p_(t+1)], only the N_p constant sequences (p, p, …, p) contribute to the difference.
- **(C)** *[H1]* Let O be bounded with link-degree ≤ d in every link, d ≤ t. Then ⟨O⟩_G and ⟨O⟩_SU(2) have equal Taylor coefficients at β = 0 through order β^(t−d).
- **(C′)** *[H1 + H2′]* For every plaquette p: Δ⟨f_p⟩ = ι_G β^t / (2^(t+1) t!) + O(β^(t+1)). The coefficient is the same for every plaquette and every lattice.
- **(D)** Let G₁ and G₂ be finite subgroups that are t₁- and t₂-designs, and t = min(t₁, t₂).
  - [H1] (A) and (C) hold between G₁ and G₂ with this t.
  - [H1 + H2′] ℓ^G1_(t+1) − ℓ^G2_(t+1) = N_p (M^G1_(t+1) − M^G2_(t+1)) / (2^(t+1)(t+1)!), and ⟨f_p⟩^G1 − ⟨f_p⟩^G2 = (M^G1_(t+1) − M^G2_(t+1)) β^t / (2^(t+1) t!) + O(β^(t+1)).
  - If t₁ < t₂ and t₁ = t(G₁), the coefficient is N_p ι_G1 / (2^(t+1)(t+1)!) > 0.
  - If t₁ = t₂ the formula holds and may vanish.

## 4. Proof

**Lemma 4.1** (degree lemma). If G is a t-design and F has link-degree ≤ t in every link, then E_G[F] = E_SU(2)[F].
*Proof.* By linearity it suffices to take F = ∏_l m_l(U_l) with deg m_l ≤ t. Under the product measure the links are independent, so E_K[F] = ∏_l E_μK[m_l]. Each factor is the same for G and for SU(2) by Definition 2.2. ∎

**Lemma 4.2** (plaquette degree). Assume (H1). For plaquettes p₁, …, p_n (repetition allowed), the product f_p₁ ⋯ f_p_n has link-degree ≤ n_l := #{i : l ∈ p_i} in each link l.
*Proof.* Write f_p = ½ Σ_{a,b,c,e} (U_{l₁}^ε₁)_{ab} (U_{l₂}^ε₂)_{bc} (U_{l₃}^ε₃)_{ce} (U_{l₄}^ε₄)_{ea}. Each factor is ± an entry of U_{lᵢ} (Definition 2.1). The four links are distinct by H1, so every monomial has degree exactly 1 in each link of p and 0 in all other links. Degrees add under products. ∎

**Lemma 4.3** (peeling). Let l′ occur exactly once in the word of p, let F be bounded and independent of U_{l′}, and let a ≥ 0. Then E_K[f_p^a F] = (M_a^K / 2^a) E_K[F].
*Proof.* Cyclically rotate the word, which leaves the trace unchanged, so that U_p = U_{l′}^ε R, where R depends only on links other than l′. By Fubini (finite sums for G; for Haar, bounded measurable functions on a product of compact probability spaces), condition on all links except l′. For fixed R, the map u ↦ u^ε R carries μ_K to μ_K, by inversion and right invariance (§1.3). So E[f_p^a | links ≠ l′] = E_μK[(tr W / 2)^a] = M_a^K / 2^a, a constant. Since F is measurable with respect to the conditioning links, the claim follows. ∎

**Lemma 4.4** (series). Let z, z′ ∈ ℂ[[β]] with z₀ = z′₀ = 1.
- (a) Define ℓ = log z by ℓ₀ = 0 and n ℓ_n = n z_n − Σ_{k=1}^{n−1} k ℓ_k z_(n−k). This is the formal log(1 + X) = Σ_k (−1)^(k+1) X^k / k, and it equals the Taylor series of the analytic logarithm wherever z is a convergent series with z(0) = 1 (composition of convergent power series, NAMED). From the recursion, ℓ_n = z_n + P_n(z₁, …, z_(n−1)) with P_n a fixed polynomial. Hence if z and z′ agree through order m − 1, then ℓ and ℓ′ agree through m − 1 and ℓ_m − ℓ′_m = z_m − z′_m.
- (b) If series a, a′ agree through order m − 1 and z, z′ agree through m − 1, then a/z and a′/z′ agree through m − 1. Moreover (a/z − a′/z′)_m = (a_m − a′_m) − a′₀ (z_m − z′_m). This follows from a/z − a′/z′ = [(a − a′) z′ − a′ (z − z′)] / (z z′), whose denominator has constant term 1. ∎

**Proof of (A).** z_n^K = (1/n!) Σ_{π ∈ 𝓟^n} E_K[∏_{i} f_{p_i}]. For n ≤ t every n_l ≤ n ≤ t, so each term is equal for G and SU(2) by Lemmas 4.2 and 4.1. Lemma 4.4(a) with m = t+1 transfers the equality to ℓ_n, n ≤ t. ∎

**Proof of (B).** Fix π = (p₁, …, p_(t+1)) and put Δ(π) = E_G[∏ f_{p_i}] − E_SU(2)[∏ f_{p_i}].
- *Case 1: n_l ≤ t for every l.* Then Δ(π) = 0 by Lemmas 4.2 and 4.1.
- *Case 2: n_l = t+1 for some l.* Then every p_i contains l. Let q₁, …, q_r be the distinct plaquettes of π, with multiplicities a₁, …, a_r, Σ a_k = t+1. All of them lie in the star of l.
  - *r = 1.* By Lemma 4.3 with F = 1 (any link of q₁ occurs once, by H1), E_K[f_q₁^(t+1)] = M_(t+1)^K / 2^(t+1). So Δ(π) = (M_(t+1)^G − M_(t+1)^Haar) / 2^(t+1).
  - *r ≥ 2.* By **H2′**, q₁ has a link l′ ≠ l that lies in no q_k with k ≥ 2, and by **H1** it occurs once in q₁. Lemma 4.3 with F = ∏_{k≥2} f_{q_k}^{a_k} gives E_K[∏_k f_{q_k}^{a_k}] = (M_a₁^K / 2^a₁) · E_K[F]. Here a₁ ≤ t because r ≥ 2, so M_a₁^G = M_a₁^Haar by Lemma 2.4(iv). F is a product of Σ_{k≥2} a_k ≤ t plaquette factors, so every link has degree ≤ t in it, and E_G[F] = E_SU(2)[F] by Lemmas 4.2 and 4.1. Hence Δ(π) = 0.

Cases 1 and 2 exhaust 𝓟^(t+1). The sequences with r = 1 are exactly the N_p constant ones. Summing, E_G[S^(t+1)] − E_SU(2)[S^(t+1)] = N_p (M_(t+1)^G − M_(t+1)^Haar) / 2^(t+1). Divide by (t+1)! and apply Lemma 4.4(a) with m = t+1, using (A); Corollary 2.5 identifies the bracket as ι_G. ∎

*Where the hypotheses enter.* H1 enters through Lemma 4.2 (degree count) and through the single occurrence of l′ needed in Lemma 4.3. H2′ enters only in Case 2 with r ≥ 2. Under H2, the argument proposed in the task also works: integrate all the other links of every q_k at once, since under H2 they are private to q_k. Each factor then becomes M_a_k / 2^a_k.

**Proof of (C).** a_n^K = (1/n!) Σ_{π ∈ 𝓟^n} E_K[O ∏ f_{p_i}]. The integrand has link-degree ≤ d + n_l ≤ d + n, which is ≤ t for n ≤ t − d, so a_n^G = a_n^SU(2) by Lemmas 4.2 and 4.1. Z agrees through β^t ⊇ β^(t−d) by (A). Lemma 4.4(b) with m = t − d + 1 gives the claim. ∎

**Proof of (C′).** O = f_p has d = 1, so (C) gives agreement through β^(t−1). At order t, a_t^K = (1/t!) Σ_{(p₁,…,p_t)} E_K[f_p f_p₁ ⋯ f_p_t]. This is a sum over sequences of t+1 plaquettes (p, p₁, …, p_t), and the case analysis of (B) applies verbatim: the difference survives only when p₁ = … = p_t = p. So Δa_t = (M_(t+1)^G − M_(t+1)^Haar) / (2^(t+1) t!). By (A), Δz_t = 0 (when t ≥ 1; for t = 0, Δz₀ = 0 trivially). Lemma 4.4(b) gives Δ⟨f_p⟩_t = Δa_t. ∎

**Proof of (D).** Apply (A), (B), (C) and (C′) to each G_i against SU(2) with the common t = min(t₁, t₂). Both G_i are t-designs. Subtracting the two (B) formulas gives N_p (M^G1 − M^G2)_(t+1) / (2^(t+1)(t+1)!), and similarly for (C′). If t₁ < t₂, then G₂ is a (t+1)-design, so M^G2_(t+1) = M^Haar_(t+1) and the coefficient is N_p ι_G1 / (…) ≠ 0 when t₁ = t(G₁). ∎

## 5. Why a sharing condition is needed, and what fails without H1

**5.1 The 1×1 periodic torus: H1 fails, and (A) fails.** The single plaquette is the commutator word U V U^(−1) V^(−1), of degree 2 in each link. Lemma 4.1 then covers only orders n ≤ t/2. Exactly, by Schur orthogonality, E_K[χ_ρ([U, V])] = 1/d_ρ for every irreducible ρ of K (NAMED, Frobenius). So E_K[f^n] depends on how D^j with 2j ≤ n *splits* over G, not only on the trace moments. For Dic₃, D¹|_G = (1-dim) ⊕ (2-dim), which gives E_G[(tr[U,V])²] = 1 + 1 + ½ = 5/2 against 1 + ⅓ = 4/3 for SU(2). So Δz₂ = Δℓ₂ = 7/48 (DERIVED by hand; **MEASURED** `DSERIES_VS_HAAR|lat=torus1x1|group=Dic3|dlogZ_1_to_t2=0,7/48,...`). Measured first differing orders: Q8, Dic₃ and Dic₄ at 2; 2T at 3 (D^(3/2)|_2T = 2′ ⊕ 2″); 2O at 4.

**5.2 The pillow: H1 holds, H2′ fails, and (B) fails by a factor 2^t.** Take four links forming a 4-cycle and two plaquettes P₁ = l₁l₂l₃l₄ and P₂ = P₁ reversed. These are distinct cells sharing all four links; the complex is S² cut into two discs. Since tr U^(−1) = tr U, S = 2 f_P₁ exactly, and U_P₁ is μ_K-distributed (Lemma 4.3). So Z_K = E_μK[e^(β tr W)], and Δz_(t+1) = 2^(t+1) ι / (2^(t+1)(t+1)!) = ι/(t+1)!. The formula of (B) gives 2ι/(2^(t+1)(t+1)!). The ratio is 2^t (DERIVED). **MEASURED**, as measured against formula: Q8 1/12 vs 1/96; Dic₃ 1/24 vs 1/192; 2T 1/720 vs 1/23040; 2O 1/40320 vs 1/5160960. (A) holds there, as it must under H1.

**5.3 The 2×1 periodic strip (one L = 1 direction): H1, H2 and H2′ fail; (A) holds, (B) fails.** Measured against formula: Q8 19/288 vs 1/96; Dic₃ 65/2304 vs 1/192; 2T 41/110592 vs 1/23040; 2O 139/34406400 vs 1/5160960. In each case Δℓ_n = 0 for n ≤ t (**MEASURED**). Mechanism (DERIVED from the character formula of §6.2): Z = Σ_ρ (a_ρ/d_ρ)². The defining representation enters through (a_(1/2)/2)², and the G-irreducible constituents of D¹|_G enter through (a_ρ/d_ρ)². Both contribute at order t+1 = 4 when F = 2.

**5.4 The 2×2 periodic torus: H2 fails, H2′ holds, and every claim holds.** (A), (B) and (C′) hold against SU(2) for Q8, Dic₃, Dic₄, 2T and 2O, and (D) holds for every pair (**MEASURED**, §6.3). This is the case the task proposed as a counterexample. It is covered by the theorem once H2 is replaced by H2′.

## 6. MEASURED checks (independent exact recount)

**6.1 Method.**
- Exact arithmetic only: Python `Fraction` and integers, plus ℚ(√d) pairs, and a Rust 1.81 integer enumerator. No float anywhere.
- Group elements come from exact quaternion closure, which also verifies orders 6, 8, 12, 16, 24, 48 and 120, unit norm, and the multiplication and inverse tables.
- Each lattice is gauge-fixed on a maximal tree. Gauge-fixing lemma (NAMED standard; proof sketch): with g_root = 1 and g_y = g_x U_(x→y) along the tree, U^g has unit tree links. Conditional on the tree links, U ↦ g_tail U g_head^(−1) preserves μ_K^⊗(non-tree), so E[F] = E[F(1_tree, V)] with V i.i.d. μ_K for every gauge-invariant F.
- The enumerator visits **all** |G|^(#free links) configurations and returns exact histograms of T = Σ_p tr U_p, together with the joint histogram with each observable's trace.
- The series z_n, ℓ_n and ⟨O⟩_n are computed exactly from the histograms through β^10.

This recount is a full configuration sum, independent of the peeling mechanism of the proof.

**6.2 The SU(2) side.** On a closed oriented surface with Euler characteristic χ_E and F faces, Z = Σ_ρ d_ρ^(χ_E) (a_ρ/d_ρ)^F, with a_ρ = E[e^(β tr U/2) χ̄_ρ(U)]. This is the character (gluing) formula: Schur orthogonality applied link by link gives d^(V−E) a^F (NAMED; Witten, CMP 141 (1991) 153, bibliographic record checked; Migdal 1975 UNVERIFIED). For SU(2), a_j = Σ_n β^n m(n, j) / (2^n n!). **The formula was validated by direct enumeration for Q8 on all six closed surfaces, through β^10, with exact equality** (rows `GLUE_Q8 … enum_equals_character_sum=1`). For SU(2) on face-transitive surfaces, ⟨f_p⟩ = F^(−1) d log Z/dβ; for the finite groups, ⟨f_p⟩ is computed directly from the joint histograms.

**6.3 Results** (`check/out/TRANSFER-CHECK.hbp`, sha256 `beb54d0a76ad6de855df68ec76d346d0df3ef848b549d21ad92d57a56cc4bab5`; 199 rows; `SUMMARY|hypothesis_lattice_predictions_checked=106|failed=0`)

| lattice | N_p | H1/H2/H2′ | configurations enumerated | measured Δlog Z_(t+1) = formula (B) |
|---|---|---|---|---|
| cube surface (S², χ_E = 2) | 6 | 1/1/1 | up to 254,803,968 (2O) | Q8 1/32 · Dic₃ 1/64 · Dic₄ 1/64 · 2T 1/7680 · 2O 1/1720320 |
| torus 2×2 | 4 | 1/0/1 | up to 254,803,968 (2O) | Q8 1/48 · Dic₃ 1/96 · Dic₄ 1/96 · 2T 1/11520 · 2O 1/2580480 |
| torus 3×3 | 9 | 1/1/1 | 1,073,741,824 (Q8) | Q8 3/64 |
| 3D open stack 1×1×2 (11 plaquettes; middle-layer links in 3 plaquettes) | 11 | 1/1/1 | 134,217,728 (Q8), 5,159,780,352 (Dic₃) | (D) Q8 vs Dic₃: 11/384 = 11·(4−3)/384 |

- In every row above, Δℓ_n = 0 for 1 ≤ n ≤ t (claim A), and Δ⟨f_p⟩ = 0 below β^t with coefficient ι/(2^(t+1) t!) at β^t, for every plaquette of every lattice (claim C′). For example, Dic₃ gives 1/96 on the cube, on the 2×2 torus and, against Q8, on all 11 plaquettes of the stack.
- On the stack, the 1×2 Wilson rectangle (d = 1) agrees between Q8 and Dic₃ through β² as (C) requires; its first difference is at β⁴.
- (D) on the cube surface: Q8–Dic₃ 1/64, Dic₃–2T 1/64, 2T–2O 1/7680, Q8–2O 1/32; the same pairs on the 2×2 torus also match.
- Haar log Z on the cube surface: ℓ = (3/4)β² − (1/64)β⁴ + (5/3072)β⁶ − (17/61440)β⁸ + O(β^10).
- **Observed outside the theorem's scope:** for Dic₃ against Dic₄ (both t = 3, equal M₄), the first difference is at β⁶. On the cube it is 1/7680 and on the 2×2 torus 1/11520, both equal to the single-plaquette value N_p (M₆^Dic3 − M₆^Dic4)/(2⁶ 6!) with M₆ = 11 and 10. On the 2×1 strip it is 19/18432, which differs from that value. This is MEASURED; the theorem's scope ends at order t+1.

**6.4 Build and run.**
- `cargo +1.81 clippy --release -- -D warnings` → exit 0; `cargo +1.81 build --release` → exit 0; rustc 1.81.0 (eeb90cda1 2024-09-04).
- Full run: `python3 -I check/transfer_check.py check/out --big`, wall time 2:41.5, exit 0.
- Source hashes (sha256): `transfer_check.py` 4f301d6c…861fa6, `enum/src/main.rs` a665486a…890df66, `enum/Cargo.toml` e52ba2cb…8e3d6d.
- Per-job histogram hashes are in rows `ENUM … out_sha256=`.

## 7. Scope of the transfer, and open items

**What the transfer carries toward the objective** (finite groups reproducing continuous-group lattice physics exactly up to a computable order, with a computable first deviation). On every finite lattice satisfying H1 and H2′ — any dimension; open boundaries; or periodic with all sides ≥ 2 (Proposition 1.2):
- Each finite subgroup G ⊂ SU(2) reproduces the SU(2) strong-coupling Taylor coefficients of log Z exactly through β^(t(G)).
- Every observable of link-degree d is reproduced through β^(t−d).
- The first free-energy deviation is the closed form N_p ι_G β^(t+1)/(2^(t+1)(t+1)!), and the first plaquette deviation is ι_G β^t/(2^(t+1) t!). Both are positive, and lattice-independent per plaquette.
- Two subgroups compare by the same closed form.
- The ladder is Dic₃ → 2T → 2O → 2I, with t = 3, 5, 7, 11. Dic₃ is the rung that covers S₃.

On a fixed finite lattice this is also a statement about the actual functions. log Z_K is real-analytic on ℝ (§1.4), so log Z_G(β) − log Z_SU(2)(β) = N_p ι_G β^(t+1)/(2^(t+1)(t+1)!) + O(β^(t+2)) as β → 0, with a lattice-dependent constant in the O-term.

**Scope:** the Taylor coefficients at β = 0 on a finite lattice — orders ≤ t+1 for log Z, ≤ t for plaquettes, ≤ t − d for degree-d observables.

**Open:**
1. **Quantitative bounds at a fixed β > 0** (the radius of the series and the constant in O(β^(t+2))): OPEN.
2. **Orders above t+1:** OPEN. The Dic₃/Dic₄ data in §6.3 show that the single-plaquette form can hold at a later order on some lattices (cube, 2×2 torus) and change on others (2×1 strip).
3. **Weak coupling, i.e. the spacing-to-zero target at fixed physical size:** OPEN for this transfer. NAMED: lattice gauge theories on the discrete subgroups 2T, 2O and 2I of SU(2) have a phase structure studied by Monte Carlo. Sources: Petcher and Weingarten, Phys. Rev. D 22, 2465 (1980), APS abstract checked, full text not read. Rebbi, Phys. Rev. D 21, 3350 (1980), finite non-Abelian groups, two-phase structure with a first-order transition, OSTI record checked. Later literature describes the finite groups as freezing at weak coupling; that attribution is secondary (e.g. inspirehep 2016362), UNVERIFIED at the primary-text level. "Bhanot–Rebbi 1980" for SU(2) subgroups: found only as a CERN preprint cited in a Rebbi lecture, UNVERIFIED. The journal Bhanot–Rebbi paper I located is Phys. Rev. D 24, 3319 (1981), on finite subgroups of SU(3) (APS abstract checked). With t(G) ≤ 11 for every finite subgroup (§2.4), the exact transfer from a finite subgroup of SU(2) reaches at most β^11 in log Z.
4. **S₃ with its own action** (not through Dic₃ ⊂ SU(2)): a transfer relative to a compact group containing S₃, such as O(2) or SO(3), is OPEN here. Remark 8.2 states the general form.

## 8. Remarks

**8.1 Non-group designs** (DERIVED by the same proof). Lemmas 4.1 and 4.2 use only the product measure and Definition 2.2. So (A) and (C) hold for any finite X ⊂ SU(2) that is a t-design under its uniform measure, X a subgroup or not. Such X exist for every t (NAMED: Seymour–Zaslavsky 1984, UNVERIFIED here). (B), (C′) and the counting in (D) use the invariance in Lemma 4.3. Gauge invariance requires a group, so for non-group X this extension concerns the formal series alone.

**8.2 Other compact groups** (NAMED: the proof carries over line by line; I did not run it). Take a compact K, a finite G ⊂ K, a unitary representation ρ of dimension N, and the action Σ_p Re tr ρ(U_p)/N. Count degree jointly in the entries of ρ(U) and their conjugates, and call G a t-design for (K, ρ) when all such monomials of total degree ≤ t have equal averages. Then (A) and (C) hold, and, under H1 and H2′, (B) holds in the form Δℓ_(t+1) = N_p (E_G[f^(t+1)] − E_K[f^(t+1)]) / (t+1)!, with f = Re tr ρ / N.

## 9. References and verification status

- P. Delsarte, J. M. Goethals, J. J. Seidel, "Spherical codes and designs", Geom. Dedicata 6 (1977) 363–388, doi:10.1007/BF03187604 — bibliographic record checked.
- Ei. Bannai, Et. Bannai, "A survey on spherical designs and algebraic combinatorics on spheres", European J. Combin. 30 (2009) 1392–1425 — bibliographic record checked.
- D. Petcher, D. H. Weingarten, "Monte Carlo calculations and a model of the phase structure for gauge theories on discrete subgroups of SU(2)", Phys. Rev. D 22 (1980) 2465 — APS abstract checked.
- C. Rebbi, "Phase structure of non-Abelian lattice gauge theories", Phys. Rev. D 21 (1980) 3350 — OSTI record checked.
- G. Bhanot, C. Rebbi, Phys. Rev. D 24 (1981) 3319, "Monte Carlo simulations of lattice models with finite subgroups of SU(3) as gauge groups" — APS abstract checked. An SU(2)-subgroup paper by Bhanot–Rebbi from 1980: UNVERIFIED.
- E. Witten, "On quantum gauge theories in two dimensions", Commun. Math. Phys. 141 (1991) 153–209 — bibliographic record checked. The gluing-formula normalisation used here is the one validated in §6.2.
- K. Osterwalder, E. Seiler, "Gauge field theories on a lattice", Ann. Phys. 110 (1978) 440–471 — bibliographic record checked. Listed for context only.
- Standard inputs, NAMED: Haar measure (uniqueness, bi-invariance), Clebsch–Gordan for SU(2), Schur orthogonality, and the classification of finite subgroups of SU(2) (Klein 1884). These are found in standard texts, e.g. Bröcker–tom Dieck.

## 10. Files

- This theorem: `/mnt/c/tmp/opus55/transfer-20261109/proof/TRANSFER-THEOREM.md`
- Checker: `/mnt/c/tmp/opus55/transfer-20261109/proof/check/transfer_check.py`; enumerator: `/mnt/c/tmp/opus55/transfer-20261109/proof/check/enum/` (Rust 1.81, std only)
- Results: `/mnt/c/tmp/opus55/transfer-20261109/proof/check/out/TRANSFER-CHECK.hbp` with its `.sha256` sidecar; run log `/mnt/c/tmp/opus55/transfer-20261109/proof/check/out/run-big.log`; job inputs and histograms in `/mnt/c/tmp/opus55/transfer-20261109/proof/check/out/jobs/`
