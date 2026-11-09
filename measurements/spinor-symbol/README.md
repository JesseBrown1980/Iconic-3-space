# The spinor symbol test, measured against its frozen prediction

The prediction and its pre-freeze addendum were frozen in [`predictions/spinor-symbol/`](../../predictions/spinor-symbol/) (merged as #8): `PREDICTION-SPINOR-SYMBOL.hbp` sha256 `c6e86dd4…` and `ADDENDUM-SPINOR-SYMBOL.hbp` sha256 `ad9ae33e…`. Three seats then measured it in exact integers, one per dimension, each writing its measured block before reading any predicted value. A fourth seat refuted all three. Every script is gate-clean (`tools/no_float_gate.py`, 0 findings), and every receipt reproduces byte for byte in a fresh directory against the frozen files on main; CI replays the three measurements on every pull request and push ([`spinor-measurement.yml`](../../.github/workflows/spinor-measurement.yml)).

| Folder | Receipt | sha256 | Compare result |
|---|---|---|---|
| [`4d/`](4d/) | `MEASURE-4D.hbp`: quaternions, Cl(4) = M₂(H), the instanton complex 1, 4, 3, Sp(1) × Sp(1), the obstruction | `c162f553…` | 45 rows: 37 AGREE, 7 AGREE_SCOPED, 1 NOT_MEASURED (the disclosure row), 0 DISAGREE. 39 measured rows pass; 12 of 12 controls fail as required |
| [`7d/`](7d/) | `MEASURE-7D.hbp`: octonions, Cl(7), G₂, the complex 1, 7, 7, 1, the cross product, the 1,344-element table group | `0987062e…` | 96 rows: 91 AGREE, 1 DISAGREE, 2 NOT_MEASURED, 1 REPORTED_REVERSED, 1 SUPPLIED. 10 of 10 controls fail |
| [`8d/`](8d/) | `MEASURE-8D.hbp`: Cl(8) = M₁₆(R), Spin(7), the complex 1, 8, 7, the Cayley form | `08b5936e…` | 57 AGREE, 3 DISAGREE, 2 NOT_MEASURED, 1 NOT_COMPARED, 1 convention followed. All 11 8D prediction rows agree; 11 must-fail controls fail and 1 must-pass control passes |
| [`refute/`](refute/) | `REFUTE-SPINOR.hbp`: reruns, recomputation with its own code, attacks on every disagreement | `d27286af…` | every reported disagreement CONFIRMED; one widened (below). Evidence, kept verbatim ([`EVIDENCE-ONLY`](refute/EVIDENCE-ONLY)) |

## What the measurement found

**Every row of the frozen prediction agrees** in all three dimensions: the Clifford squares γ(ξ)² = −|ξ|²·I as polynomial identities (640, 1,792 and 9,216 coefficients), the algebras M₂(H), M₈(R) ⊕ M₈(R) and M₁₆(R), the symbol sequences 1, 4, 3 and 1, 7, 7, 1 and 1, 8, 7 with ranks 1, 3 and 1, 6, 1 and 1, 7 at every ξ ≠ 0, the projection constants c = 1/2, 1/3, 1/4, the equivariance under Sp(1) × Sp(1), der(O) and spin(7), and the obstruction d_A d_A = [F_A, ·].

**Four rows of the pre-freeze addendum disagree.** They stay as frozen; the measured statements sit beside them here.

| Addendum row | Frozen text says | Measured |
|---|---|---|
| R6, the rolled-up symbols | M (4D), R (7D) and N (8D) "and their doubled self-adjoint forms square to +\|ξ\|²" | The doubled forms square to +\|ξ\|²·I in all three dimensions, and RᵀR, MᵀM, NᵀN equal +\|ξ\|²·I. The literal products differ: R(ξ)R(ξ) = −\|ξ\|²·I₈ (as prediction row 7D.P9 states), and M(ξ)M(ξ) and N(ξ)N(ξ) are non-scalar matrices. The 7D seat reported the disagreement; the refuter found the same in 4D and 8D, where the seats had read "square" as the Gram product |
| R2, the 8D symbol | σ(ξ)a = conj(ξ)·a, with Λ²₇ taken through the Cayley form | In the convention the prediction fixes, the symbol is N(ξ)a = −ξ·conj(a), which row 8D.P9 states and which agrees. conj(ξ)·a fails on 50 of 512 coefficients. It holds only with η′ in S⁻ (Φ₋₁) |
| R2, its note | conj(ξ)·a and −ξ·conj(a) are one symbol under two identifications of S± with O | The space of fixed identifications U with U·L_ξ = L_conj(ξ)·U for all ξ has dimension 0. The note holds if conj(ξ)·a is read as a·conj(ξ), right multiplication, or through the reflection ξ ↦ conj(ξ) of R⁸ (determinant −1, which exchanges Φ₊₁ and Φ₋₁) |
| R3, the Cl(8) block | the generator is [[0, L_conj(ξ)], [L_ξ, 0]] | That block squares to +\|ξ\|²·I₁₆, so it generates Cl(R⁸, +\|ξ\|²), also of dimension 256. The convention fixes Cl(R⁸, −\|ξ\|²), whose generator is [[0, L_ξ], [−L_conj(ξ), 0]] (9,216 coefficients) |

Two further observations, both consistent with the frozen rows:
- **Orientation.** The frozen T₇ split (+2 on 7 dimensions, −1 on 14) holds with vol = −e₁∧…∧e₇ on the frozen octonion table, as convention row 7D.C3 allows; with vol = e₁∧…∧e₇ the eigenvalues are −2 and +1. All four seats found the same reversal, and L_e₁⋯L_e₇ = −I₈, so ε₇ = −1.
- **Format.** Addendum rows R1, R2, R3 and R6 hold the bar character inside values (the norm written with bars), although the prediction's FORMAT row reserves the bar as the field separator. The seats matched quotes against whole rows.

## Independence, as measured

Each seat's measured block was written and hashed before it read any PRED, ADD or DIMS value, and each snapshot is reproduced by its kept pre-read script (`measure_4d.before-compare.py`, `measure_7d.precompare.py`, `measure_8d.before-read.py.txt`, `refute_spinor.before-compare.py`). The refuter adds two limits, carried here: file times cannot show what a seat read through `git show` or its own tool history (CANNOT_SEE), and the task text given to the seats named some figures of the addendum (c = 1/4, the sequence 1, 8, 7, the associator triple), so the blind blocks contain addendum content computed from definitions. The seats worked without GPT-6's independent receipt `c28f8af6…`, by design.

## Scope

These are exact identities of symbols: linear algebra over the integers and Z[i], at the level of the principal symbol and of single gauge-field configurations. They show the Clifford and instanton-complex structure that the owner's "square root of curvature" picture names, in 4, 7 and 8 dimensions. Open here (NOT_MEASURED): the operator-level statements (the complex holding exactly on instantons, the operator part of ADD R5) and the Bott periodicity rows with the plus sign.
