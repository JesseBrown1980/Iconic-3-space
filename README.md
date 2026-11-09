# Iconic 3-space

> I'm trying to prove the illusion of continuous arises from the finite world.
>Iconic and whereby Ionic 3 space, the 24-dimensional world, which is laden with exceptions and these groups, but because the entire subject presumes defining representations. You know, when you have, say, SU3, you imagine complex 3 space as the raison d'etre of that group's existence. The problem here is that these are the tiny number of groups for which this paradigm, that is, deprojectivized linear representations, fails. And so by virtue of the pedagogy, you are led into thinking that somehow you're looking at the platypus and the echidna and saying, well, these are just, you know, sort of exception mammals, when in fact they're the most important mammals of all because they predate the marsupial branching on the phylogenetic tree.
> — Jesse Daniel Brown, November 9, 2026

This repository tests that statement with exact finite mathematics. Its question: how far can a finite set of symmetries stand in for a continuous one before an observer can tell them apart?

The name comes from the owner's phrase of the same day about "the 24-dimensional world", read here as the Leech lattice built as 3-space over the icosians. His full words are in [`OBJECTIVE-QUOTES.md`](OBJECTIVE-QUOTES.md).

**Mathematical foundations:** [`Icosian 3-space, projective lifts and finite-to-continuous tests`](docs/ICOSIAN-3-SPACE-FOUNDATIONS-20261109.md) records the 24-dimensional construction, exact representations/design distinctions, concrete falsification gates and primary literature. The transcript retains the owner's original wording, with editorial interpretations separately labeled. The SU(2) prediction freeze at [`predictions/su2-ladder/`](predictions/su2-ladder/) remains pinned to its original commit and receipt.

## The frame

The owner's frame (OPERATOR-CANON): "we're in a matrix, so the hashes matter too. Like quantum entanglement as a hash." Here that frame is part of the method. A result lives in a receipt pinned by its hash, and a receipt counts only when its script reproduces that hash byte for byte. Every prediction is frozen by its hash before anything is measured, so the record of what was expected can be checked against what was found. The quantum lane will test entanglement as a hash directly. A random enough set of operations, a 2-design, scrambles a system so that each small piece looks completely random on its own, while the information stays whole in the correlations between the pieces: the decoupling theorem. Its prediction will be frozen before it is run.

## The test

A finite set X of points on a sphere, or of elements of a group such as SU(2), **looks exactly continuous up to degree t** when averaging any polynomial of degree at most t over X gives the same number as averaging it over the whole sphere or group. Such a set is called a t-design. Past degree t, some observable sees the difference. That observable, and its exact value, mark the first place the finite world shows through.

The work runs in three steps, each with its own commit and receipt:

1. **Frozen predictions.** For each finite set, the predicted t and the exact deviation at the first break, computed from closed formulas alone, before any point set is built. These are committed first, so their timestamps come before any measurement.
2. **Measurement past the limit.** Exact integer scripts build each finite set and compute every moment through, and well beyond, the predicted break.
3. **Independent refutation.** A second seat reruns every script, checks that each receipt reproduces byte for byte, and attacks each claim.

## What is being tested

- **The SU(2) ladder.** The cyclic groups, the binary dihedral groups (Q8 among them), and the three exceptional groups 2T, 2O and 2I, which the McKay correspondence pairs with E6, E7 and E8.
  - Each group's design degree and first break.
  - The full invariant series out to degree 120.
  - The strong-coupling free energy of 2D lattice gauge theory with that group, compared order by order with SU(2).
- **The exceptional point sets.**
  - The 600-cell, whose 120 vertices are the unit icosians (2I).
  - The 240 roots of E8.
  - The 196,560 shortest vectors of the Leech lattice.
- **The octonion side.**
  - G2 = Aut(O) acting on S⁶ = G2/SU(3), and the almost complex structure the octonions give it.
  - The 16 octonion units, a Moufang loop.
  - The finite frame symmetry inside G2.
- **The quantum side.** 2O is the single-qubit Clifford group. Its predicted break at degree 8 is the statement that it reproduces random single-qubit statistics through the third moment and first differs at the fourth.

## Rules

- **Integer only.** Every tracked file is scanned by [`tools/no_float_gate.py`](tools/no_float_gate.py), the same gate the companion Yang–Mills repositories use, copied verbatim. Rationals are reduced integer pairs, and algebraic numbers are integer vectors over their rings.
- **Hot path first.** Results are `.hbp` rows (`TAG|k=v|...|json=0`) with a `.sha256` sidecar.
- **Every claim carries a tag:**
  - MEASURED: a script here computed it.
  - NAMED: standard literature, with the citation.
  - CONJECTURE: derived but not yet run.
- **Corrections travel beside a claim, never over it.**
- **Predictions are frozen before measurement.**

## Status, November 9, 2026

| Step | Status |
|---|---|
| Frozen prediction: the SU(2) ladder, 61 groups | **frozen**, [`predictions/su2-ladder/`](predictions/su2-ladder/), `PREDICTION.hbp` sha256 `c2601cf0…`, receipt [`FREEZE.hbp`](predictions/su2-ladder/FREEZE.hbp) |
| Frozen prediction: the 600-cell, E8 roots and Leech minimal vectors | **frozen**, [`predictions/exceptional/`](predictions/exceptional/), `PREDICTION-EXCEPTIONAL.hbp` sha256 `e8403cd4…`, receipt [`FREEZE.hbp`](predictions/exceptional/FREEZE.hbp) |
| Frozen prediction: hash-phase quantum states (global, local and quadratic hashes) against continuous random states, 1 to 8 qubits | **frozen**, [`predictions/hash-quantum/`](predictions/hash-quantum/), `PREDICTION-HASH.hbp` sha256 `338df23a…` |
| Measurements past each break | after the freezes |
| Independent refutation | after the measurements |
| Proofs: quaternions give SU(2), Q8, the 16 octonion units, the frame group inside G₂, the ladder's ceiling | **measured and refuted**, [`proofs/quaternion-octonion/`](proofs/quaternion-octonion/): S1, S2, S3 and S6 proven, S4 an analogy made exact, S5 corrected |
| Proofs: G₂/SU(3) = S⁶, its octonion almost complex structure (not integrable), W(G₂), the G₂ Higgs spectrum, G₂-instantons, the finite-subgroup limit | **measured and refuted**, [`proofs/g2-s6/`](proofs/g2-s6/): C1 to C7 supported, C8 to C10 restated |

### The frozen SU(2) ladder prediction

The prediction was computed from the closed invariant (Molien) series alone: Klein 1884, and Springer, *Invariant Theory*, LNM 585, 1977. No group element was built. Each value below is a PREDICTION that the measurement can fail.

| Group | Order | Looks exactly like SU(2) up to degree | First break | Moment there, group : SU(2) | 2D gauge free energy, first difference |
|---|---|---|---|---|---|
| C_n, every n = 3…30 | n | 1 | degree 2 | 2 : 1 | β² term, 1/8 |
| C_2 | 2 | 1 | degree 2 | 4 : 1 | β² term, 3/8 |
| Q8 = Dic_2 | 8 | 3 | degree 4 | 4 : 2 | β⁴ term, 1/192 |
| Dic_n, every n = 3…30 | 4n | 3 | degree 4 | 3 : 2 | β⁴ term, 1/384 |
| 2T | 24 | 5 | degree 6 | 6 : 5 | β⁶ term, 1/46,080 |
| 2O | 48 | 7 | degree 8 | 15 : 14 | β⁸ term, 1/10,321,920 |
| 2I | 120 | **11** | degree 12 | 133 : 132 | β¹² term, 1/1,961,990,553,600 |

The ceiling the prediction states: of the 61 groups, only 2I reaches degree 11. The cyclic and binary dihedral families stay at 1 and 3 at every size. Dic_30 and 2I both have 120 elements: Dic_30 is predicted to stop at degree 3, 2I at degree 11. The file also holds each group's full invariant series to degree 120, its moments to order 60 and its gauge free energy to β³⁰, so the measurement can be checked well past each break.

### The frozen exceptional prediction

Computed from the published distance distributions alone (Conway and Sloane, *Sphere Packings, Lattices and Groups*; Delsarte, Goethals and Seidel 1977), with no vector built. A set is a t-design when its zonal sums S_k vanish for k = 1 to t.

| Point set | Points | Sphere | Looks exactly like the sphere up to degree | First break, S_k | Nonzero S_k through k = 30 |
|---|---|---|---|---|---|
| 600-cell = the 120 unit icosians = 2I | 120 | S³ | 11 | S₁₂ = 120/13 | k = 12, 20, 24, 30 |
| E8 roots | 240 | S⁷ | 7 | S₈ = 720/143 | k = 8, 12, 14, 16, …, 30 (S₁₀ = 0) |
| Leech lattice shortest vectors | 196,560 | S²³ | 11 | S₁₂ = 17955/3596 | k = 12, 16, 18, …, 30 (S₁₄ = 0) |

The two frozen predictions cross-check. The 600-cell is 2I, so its zonal sums must agree with the SU(2) ladder's invariant counts. They do at every degree from 1 to 30: S_k (k+1) = 120 m_k, where m_k is the number of 2I-invariant polynomials of degree k. The first nonzero term is Klein's icosahedral invariant of degree 12.

## Companions

- [Finite-universe Yang–Mills with HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C): finite S₃ lattice gauge theory, spectral certificates and the first-gap proofs.
- [Finite-Universe-SU-N-solved](https://github.com/JesseBrown1980/Finite-Universe-SU-N-solved): the SU(2) finite-subgroup ladder Q8, 2T, 2O and 2I in lattice gauge theory.
- [yang-mills-4-dimensional-continum-proof](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof): the four-dimensional continuum programme.

The design facts tested here are classical. This repository contributes exact, reproducible certificates of them, frozen predictions that the measurements can fail, and the link to the finite-universe programme above.
