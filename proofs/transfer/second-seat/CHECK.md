# Second-seat check of TRANSFER-THEOREM.md

2026-11-09 · second seat for ACER-CLAUDE-FABLE5 (pid `8467a937cba309f7`), Opus 5.5 · E=0: local files only, no GitHub writes, no Rust build.

- Theorem file: `/mnt/c/tmp/opus55/transfer-20261109/proof/TRANSFER-THEOREM.md`, sha256 `aa93b433168caef4c077d5a2658bb7d8d169ea3453fc29d407daecdac12a1cbe`. I recomputed this hash and it matches (MEASURED).
- Prover's results: `proof/check/out/TRANSFER-CHECK.hbp`, sha256 `beb54d0a76ad6de855df68ec76d346d0df3ef848b549d21ad92d57a56cc4bab5`. Recomputed, matches (MEASURED).
- My hot-path output: `check/out/SECOND-SEAT.hbp`, sha256 `f7bfe8c09ab53dcb8d217b8e7aeb0181b2cb51860d5f5234483eeda18868a785` (sidecar `.sha256`).
- My script: `check/second_seat_check.py`, sha256 `43a5f939aef83230390ea5707ee61216fdb266df625f6097f3142b86ef1e7380`. It was run as `wsl.exe bash check/run_second_seat.sh` (`python3 -I`, wall time 0:40.6, exit 0). It uses Fraction, int and Q(√d) pairs only, with no float. It was written from scratch and neither imports nor reads the prover's checker. Haar moments come from Weyl constant-term extraction, not the Clebsch–Gordan recursion. Design strength is checked against the DGS monomial definition directly.

## Verdict: PROVEN_WITH_CORRECTIONS

I read the proof line by line, §1 through §5, and found no mathematical gap in (A), (B), (C), (C′) or (D). Every exact check I ran agrees with the theorem. The corrections below concern wording and bookkeeping; none of them changes a coefficient.

The scope of what is proven is strong-coupling Taylor coefficients at β = 0 on a fixed finite lattice: log Z through β^(t+1), plaquettes through β^t, and degree-d observables through β^(t−d). Fixed-β bounds and weak coupling remain OPEN, as stated in §7.

## 1. Line-by-line answers

**Design property and function class.** Lemma 4.1 is correct. Under the product measure, a finite sum of products ∏_l m_l(U_l) has expectation ∏_l E[m_l]. Equality therefore needs degree ≤ t *per link*, and other links in the integrand do no harm. Definition 2.1 is right that the entries of U, U⁻¹ and Ū, and the four real coordinates, give the same degree class for SU(2). Lemma 4.2 gives f_p degree exactly 1 in each of its four links under H1. **DERIVED, checked.**

**Order of integration and dependence structure.** Lemma 4.3 conditions on every link except one private link l′ and uses invariance of μ_K under u ↦ u^ε R, with R ∈ K. This is plain Fubini on a finite or compact product space. Independence of F from U_{l′} is exactly what H2′ provides inside the star of l. It is correct. **DERIVED, checked.**

**Order t+1 with a link of degree t+1 shared by two or more distinct plaquettes.** Case 2 with r ≥ 2 is correct. Peeling the private l′ of q₁ gives M_{a₁}/2^{a₁} with a₁ ≤ t, and the remaining F has at most t factors. So Δ(π) = 0.

I tested this directly (MEASURED, rows `MIXED`). E[f_p² f_q²] = 1/16 = (M₂/4)², and E[f_p³ f_q] = E[f_p f_q³] = 0, exactly as peeling gives, for:
- Dic₃ and 2T on the cube surface, where p and q share 1 link;
- Dic₃ and 2T on the 2×2 torus, where p and q share 2 links (H2 fails, H2′ holds).

On the pillow, where H2′ fails, the same moments come out as 3/16 (Dic₃) and 1/8 (2T). That is the mechanism that breaks (B) there. The constant-sequence claim in E[S^(t+1)] also checks out: E_Dic₃[S⁴] − E_Haar[S⁴] = N_p·ι/16, measured as 1/16, 3/8 and 1/4 on the single plaquette, the cube and the 2×2 torus (rows `ES4_DIC3_MINUS_HAAR`).

**Passage from Z to log Z.** Lemma 4.4(a) is correct. From z′ = ℓ′z I re-derived n ℓ_n = n z_n − Σ_{k<n} k ℓ_k z_{n−k}. So ℓ_n − z_n is a polynomial in z₁ … z_{n−1}, and the first differing coefficients are equal. Lemma 4.4(b) is also correct: the numerator vanishes below order m, and the denominator has constant term 1. **DERIVED, checked.**

**Does (C) hold as stated.** Yes, with d the maximum per-link degree, which is what "≤ d in every link" says. I tested it on the cube with Dic₃ against 2T (t = 3), rows `C_DEGREE`, MEASURED:
- O = f_p² (d = 2): (C) predicts agreement through β¹. The series agree through β¹ and first differ at β² by 1/32, which is the f_p⁴ term (M₄ difference 1, divided by 16·2!).
- The Petrie-hexagon Wilson loop W/2 (d = 1, six links): (C) predicts agreement through β². The series agree through β⁴ and first differ at β⁵ by 1/256.

**Also checked and correct** (DERIVED):
- Lemma 1.1.
- Proposition 1.2, including H2′ at L = 2.
- Lemma 2.4, all four equivalences. (iv) ⇒ (iii) correctly uses m(n, j) ≥ 0 and m(n, n/2) = 1.
- Corollaries 2.5 and 2.6.
- The cyclic and binary-dihedral family computations (Dic_n: t = 3, ι = 2 for n = 2 and ι = 1 for n ≥ 3).
- Dic₃/{±1} ≅ S₃ (a³ = −1 and b² = −1 map to the identity, and bab⁻¹ = a⁻¹).
- The Frobenius computation in §5.1 (5/2 against 4/3, giving Δz₂ = 7/48).
- The pillow factor 2^t in §5.2.

## 2. Prediction arithmetic (MEASURED, rows `PAIR_RECOMPUTE`)

I recomputed all 45 frozen `PAIR` rows of `predict/PREDICTION-TRANSFER.hbp` (sha256 `8ada45c8…6d`) from my own moments: 45 rows, **0 mismatches**. Spot entries:
- cube, Dic₃ − 2T: 6·(3−2)/(2⁴·4!) = **1/64** at β⁴, plaquette mean 1/96. This matches the frozen row and my lattice enumeration (`DIC3_VS_2T|lat=cube_surface`, dlogZ = 0,0,0,1/64).
- cube, 2O − SU2: 6/(2⁸·8!) = 1/1720320 at β⁸, plaquette mean 1/1290240. Matches.
- cube, 2I − SU2: 6/(2¹²·12!) = 1/326998425600 at β¹², plaquette mean 1/163499212800. Matches.
- 3×3 torus, Q8 − 2T: 9·2/384 = 3/64 at β⁴. Matches.

**Ladder** (MEASURED, rows `GROUP2`): the moment strength equals the DGS strength for every group.

| Group | Order | t | M_(t+1) | Haar | ι |
|---|---|---|---|---|---|
| Z6 | 6 | 1 | 2 | 1 | 1 |
| Q8 | 8 | 3 | 4 | 2 | 2 |
| Dic₃ | 12 | 3 | 3 | 2 | 1 |
| Dic₄ | 16 | 3 | 3 | 2 | 1 |
| 2T | 24 | 5 | 6 | 5 | 1 |
| 2O | 48 | 7 | 15 | 14 | 1 |
| 2I | 120 | 11 | 133 | 132 | 1 |

Haar moments from the Weyl route are 1, 0, 1, 0, 2, 0, 5, 0, 14, 0, 42, 0, 132 (Catalan).

## 3. Tiny independent test: single plaquette, free boundary (MEASURED)

Four links, word l₀ l₁ l₂⁻¹ l₃⁻¹, with full enumeration of all |G|⁴ configurations (Dic₃ 20,736; 2T 331,776). Gauge-fixed on a 3-link tree, the remaining |G| configurations give a trace histogram equal to the full one divided by |G|³ (`equal=1` for both groups), so the gauge-fixing lemma holds on this lattice.

| Comparison | Δlog Z, β¹…β⁶ | Δ⟨f_p⟩, β⁰…β⁵ | Theorem |
|---|---|---|---|
| Dic₃ − SU2 | 0, 0, 0, **1/384**, 0, −1/5120 | 0, 0, 0, **1/96**, 0, −3/2560 | (A) zero through β³; (B) 1/384; (C′) 1/96. **Holds** |
| 2T − SU2 | 0, 0, 0, 0, 0, **1/46080** | 0, 0, 0, 0, 0, **1/7680** | (A) zero through β⁵; (B) 1/46080; (C′) 1/7680. **Holds** |
| Dic₃ − 2T | 0, 0, 0, **1/384**, 0, −1/4608 | 0, 0, 0, **1/96**, 0, −1/768 | (D) 1/384, 1/96. **Holds** |

Coefficients above order t+1 (for example −1/5120 at β⁶) lie outside the theorem's scope and are recorded only.

**Beyond the task, same script** (MEASURED). For Dic₃ and 2T:

| Lattice | Configurations | Result |
|---|---|---|
| Cube surface (H1/H2/H2′ = 1/1/1) | 248,832 (Dic₃), 7,962,624 (2T) | Δℓ₄ = 1/64 (Dic₃) and Δℓ₆ = 1/7680 (2T), equal to (B). (C′) gives 1/96 and 1/7680. |
| 2×2 torus (1/0/1) | 248,832 (Dic₃), 7,962,624 (2T) | Δℓ₄ = 1/96 and Δℓ₆ = 1/11520, equal to (B) under H2′. Dic₃ − 2T = 1/96, equal to (D). |
| Pillow (1/0/0) | 12 (Dic₃), 24 (2T) | 1/24 against formula 1/192, and 1/720 against 1/23040: (B) fails by 2^t, as §5.2 says. (A) holds. |

The Haar log Z on the cube, from my own Weyl-route character sum, is 3/4, −1/64, 5/3072, −17/61440 at β², β⁴, β⁶, β⁸. This equals §6.3. The character formula also reproduces E[e^(β tr W)] exactly on the pillow.

## 4. Corrections (each with its check)

1. **H2′ wording is weaker than what the checker tests.** §1.2 defines H2′ without "l′ occurs exactly once in p". The prover's `hypotheses()` (transfer_check.py lines ~355–373) requires single occurrence. Under the §1.2 wording, the 2×1 strip and the 1×1 torus *satisfy* H2′; under the single-occurrence form they fail it (MEASURED, rows `HYP_WORDING`: as_written = 1 and single_occurrence = 0 on both). So the statements "the 2×1 strip breaks H1, H2 and H2′" (§0 item 3, §5.3) and the LAT row "1×1 torus 0/1/0" are false as worded. Under H1 the two forms coincide, so (B), (C′) and (D) are unaffected. Fix: add "occurring exactly once in p" to H2′, or say the strip breaks H1 and H2.
2. **Correction 1 overgeneralises.** "A periodic direction with L = 1 breaks [H1], and then claim (A) itself fails" is contradicted by the prover's own rows: on the 2×1 strip, which has an L = 1 direction, `A_holds=1` for Q8, Dic₃, Dic₄, 2T and 2O (MEASURED in TRANSFER-CHECK.hbp). Correct form: without H1, (A) is no longer guaranteed. It fails on the 1×1 torus and holds on the 2×1 strip.
3. **"106 of 106 predictions" mixes frozen predictions with post-freeze checks.** The 106 is a count of checks (transfer_check.py lines 668–785). By my count against the frozen file (sha `8ada45c8…`):
   - About 48 checks correspond to frozen rows, up to 54 if partial rows count. These are the cube rows for Q8, Dic₃, 2T and 2O; the four frozen pairs; the Dic₃–Dic₄ zero at β⁴; the 3×3 Q8 row; and the (A) part on the 2×2 torus.
   - About 52–58 were added after the freeze:
     - every Dic₄ row (Dic₄ appears in the frozen file only in CORRECTION_EXAMPLE);
     - the 1×1×2 stack, which is absent from the frozen file;
     - all (B), (C′) and (D)-coefficient checks on the 2×2 torus, which the frozen file marked `coefficient_prediction=NONE` and `formula_hypothesis=FAILS`.

   All of them pass. Fix: report "106 checks, of which N were frozen predictions" and tag the 2×2-torus coefficient checks as confirmations of the H2′ extension made after the freeze.

   Also NOT_RUN from the frozen file:
   - all 15 rows on the 4×4 torus;
   - 14 of 15 rows on the 3×3 torus (only Q8 − SU2 was measured);
   - every 2I lattice row. The 2I deviations are DERIVED values, not enumerated ones.
4. **"H2 implies H2′"** (task text) needs H1 as well. A word made of a single repeated link has no l′ ≠ l. The theorem's Lemma 1.1 states this correctly as H1 + H2 ⇒ H2′.
5. **(C) clarity.** O should be named as a fixed polynomial on SU(2)^𝓛, restricted to G^𝓛. Every function on a finite G is a polynomial of some degree, so the degree is that of the chosen SU(2) polynomial. The proof already uses it this way.

## 5. What I did not re-measure

- The prover's Rust 1.81 enumerator and its gate `cargo +1.81 clippy --release -- -D warnings` = exit 0: NAMED, from the theorem's own text. I did not re-run it, to keep E=0 and launch no build.
- The 5,159,780,352-configuration stack run, Dic₄, 2O, and the 3×3 torus: NAMED, from TRANSFER-CHECK.hbp. My enumerations cover Dic₃ and 2T on the single plaquette, pillow, cube and 2×2 torus.
- The standard inputs (Haar uniqueness, Clebsch–Gordan, Schur, Klein's classification, the gluing formula in general): NAMED. My checks validate the gluing formula only where they reproduce independent values, namely the pillow (exactly) and orders 2 and 4 on the cube and the 2×2 torus, via (A) against 2T.

## 6. Falsifiers

Any one of the following would refute the theorem:
- A lattice satisfying H1 and the single-occurrence H2′, with a finite subgroup G, where an exact enumeration gives Δℓ_n ≠ 0 for some n ≤ t, or Δℓ_(t+1) ≠ N_p ι_G/(2^(t+1)(t+1)!).
- A product of t+1 plaquette factors from two or more distinct plaquettes in one star where E_G ≠ E_SU2.
- An observable of per-link degree d that differs before β^(t−d+1).

Next falsifier still NOT_RUN: Dic₃ against SU(2) on the 3×3 torus (12¹⁰ ≈ 6.2·10¹⁰ gauge-fixed configurations). The theorem predicts Δℓ₄ = 9/384 = 3/128 there. Open patches are weak tests, because tree gauge-fixing makes their plaquettes independent; use closed surfaces or periodic lattices.
