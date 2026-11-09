# Icosian 3-space: foundations and finite-to-continuous research gates

**Project:** [Iconic 3-space](https://github.com/JesseBrown1980/Iconic-3-space)  
**Research objective and project attribution:** Jesse Daniel Brown — November 9, 2026  
**Evidence class:** NAMED mathematical background / PROJECT TEST PLAN. The frozen numerical predictions are kept in `predictions/` with their original SHA receipts.

## 1. Reading the source phrase

The submitted transcript says **“Ionic 3 space, the 24-dimensional world”**. The editorial interpretation **icosian 3-space** refers to a rank-three lattice over the icosian quaternionic order, whose associated rational trace lattice has real dimension **24**. The source wording is retained in [OBJECTIVE-QUOTES.md](../OBJECTIVE-QUOTES.md).

The **icosian ring** 𝕀 lives in Hamilton quaternions over `Q(√5)`. Each coordinate contains four quaternionic components with two rational coefficients each, giving rank **8 over Z**. Three coordinates form a rank-24 additive lattice. The Leech lattice arises when the correct integrality, congruence and inner-product conditions are imposed; an arbitrary triple of icosians is a larger ambient module.

For example, John Baez describes an icosian triple construction with congruences `x ≡ y ≡ z (mod h)` and `x+y+z ≡ 0 (mod h*)`. The result, with the indicated lattice normalization, is the **Leech lattice Λ24**.

**Source anchors:** [Baez, *This Week's Finds*, week 59 (1995)](https://math.ucr.edu/home/baez/twf_html/week59.html); [Wilson, *Octonions and the Leech lattice* (2009), introduction](https://webspace.maths.qmul.ac.uk/r.a.wilson/pubs_files/octoLeech1rev.pdf); [Conway–Sloane, *Laminated lattices* (1982)](https://annals.math.princeton.edu/1982/116-3/p07).

## 2. What “exceptional” and “projective” mean here

A defining representation of `SU(3)` acts on `C^3`. A group can also be most naturally presented through a lattice, root system, quaternionic/octonionic geometry, or a projective action. “Exceptional” indicates a distinct mathematical structure rather than the absence of linear representations: exceptional Lie groups have faithful ordinary representations in suitable dimensions.

A **projective representation** `ρ̄: G → PGL(V)` composes only up to a scalar phase when lifted to matrices. If representatives obey `ρ(g)ρ(h) = c(g,h)ρ(gh)`, then `c` is a cocycle. Lifting to a genuine linear action requires checking that its class vanishes, or passing to an appropriate **central extension** of `G`.

The binary polyhedral ladder `2T, 2O, 2I` consists of finite subgroups of `SU(2)` and gives a concrete arena for comparing projective rotations with linear quaternionic double covers. For the icosian Leech lattice, the quaternionic rank-three setting also exposes the Hall–Janko double cover `2.J2`.

**Special recent source:** [Gerald Höhn, *The Hall–Janko Simple Group and the Icosian Leech Lattice* (September 2026), arXiv:2609.22842](https://arxiv.org/abs/2609.22842) studies the rank-three icosian Hermitian lattice, its `2.J2` automorphisms and quaternionic projective configurations. Its abstract identifies three Hermitian components of the Leech minimal shell as **quaternionic projective 5-designs**. Those projective moments are a separate invariant from real-spherical polynomial design moments.

## 3. The finite-to-continuous test

For a finite set `X` on a sphere or in a compact group `K`, compare exact finite averages with invariant continuous averages of the *same* observables. For a polynomial `P`:

```text
A_X(P) = (1/|X|) Σ_{x in X} P(x)
A_K(P) = ∫_K P(g) dμ(g)
Δ_P = A_X(P) − A_K(P)
```

A **t-design** matches every allowed polynomial moment through degree `t`. The smallest degree with a certified `Δ_P ≠ 0` is a **first distinguishing observable**. For compact Lie groups, specify the representation and the precise polynomial/function class; for projective spaces, use projectively meaningful invariant moments.

The first frozen example is `2I` (the 120 vertices of the 600-cell) compared with `SU(2) ≅ S³`. Its [already-frozen prediction](../predictions/su2-ladder/PREDICTION.hbp) is **matching through degree 11**, with the **first distinguishing invariant at degree 12**. This is a prediction awaiting the independent element-by-element measurement.

### Required independent gates

1. **Source gate (frozen):** preserve the existing `2I`, `2O`, `2T`, cyclic and binary-dihedral Molien-series prediction receipts and their commit hashes. Next freeze the separate `E8`, 600-cell and Leech-shell predictions before reconstructing their point sets.
2. **Exact algebra gate:** construct each finite group or minimal shell by an independent method; verify element count, closure where applicable, norm, minimal-shell multiplicity, and conjugation/action laws.
3. **Moment gate:** compute invariant moment sums in reduced integer or algebraic-integer arithmetic, compare with independent continuous integration formulas, and find the first nonzero difference.
4. **Projective gate:** for `2.J2` and the icosian Leech lattice, evaluate the quaternionic projective Gram-angle/line design conditions separately from real-spherical design conditions. Test the projective 5-design theorem with exact witnesses.
5. **Adversarial gate:** inject a broken multiplication, a missing shell vector, a mistaken central sign and a malformed projective lift. Each must have its own deterministic rejected receipt.
6. **Scale gate:** if making a physical continuum-emergence claim, define a family `X_n`, observables, an error metric and a scaling law that proves convergence uniformly over the physically relevant observables. The moment identities above give the finite-level tests for that theorem.

Every execution receives `PREDICTION`, `MEASURED`, `NAMED`, or `CONJECTURE` status and SHA-pinned HBP/HBI receipts. The finite Yang–Mills and 4D companion repositories are separate later applications.

## 4. Literature trail

- J. Tits, *Quaternions over Q(√5), Leech's lattice and the sporadic group of Hall–Janko* (1980), as referenced in Höhn's 2026 paper.
- J. H. Conway and N. J. A. Sloane, [*Laminated lattices*](https://annals.math.princeton.edu/1982/116-3/p07), *Annals of Mathematics* 116 (1982), 593–620.
- John Baez, [*Week 59: Icosians and the Leech lattice*](https://math.ucr.edu/home/baez/twf_html/week59.html), 1995.
- Robert A. Wilson, [*Octonions and the Leech lattice*](https://webspace.maths.qmul.ac.uk/r.a.wilson/pubs_files/octoLeech1rev.pdf), 2009.
- Gerald Höhn, [*The Hall–Janko Simple Group and the Icosian Leech Lattice*](https://arxiv.org/abs/2609.22842), September 2026.

**Provenance:** the research theme and supplied wording are credited to Jesse Daniel Brown. The interpretation “icosian” is a mathematical editorial reading and is tracked apart from the verbatim transcript. The mathematical objects and previously published theorems retain their original scholarly citations.
