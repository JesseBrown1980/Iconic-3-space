# Quaternions, octonions and the SU(2) ladder: exact proofs

The owner asked on November 9, 2026 for proofs of these statements. Each one was measured by an exact integer script and then attacked by an independent second seat, whose verdict appears in the Status column ([`refute/REFUTE2.hbp`](refute/REFUTE2.hbp)). Corrections sit beside the original text in [`CORRECTIONS.hbp`](CORRECTIONS.hbp).

| | Statement | Status | The exact version |
|---|---|---|---|
| S1 | Quaternions give SU(2). | **Proven, measured** | φ(a + bi + cj + dk) = [[a + bi, c + di], [−c + di, a − bi]] is an isomorphism from the unit quaternions S³ onto SU(2). It is multiplicative on all 16 basis pairs, which proves it for every pair by bilinearity, and det φ(q) = N(q) holds as a polynomial identity. Aut(H) = SO(3) = SU(2)/{±1} (Skolem–Noether, NAMED). |
| S2 | The imaginary units of the quaternions generate the group Q8. | **Proven, measured** | i and j generate Q8 = {±1, ±i, ±j, ±k}: order 8, non-abelian, associative on all 512 triples. The set of all unit imaginary quaternions generates the whole of SU(2). |
| S3 | The seven imaginary octonion units don't form a group, because octonions aren't associative. | **Proven, measured** | The 16 units ±1, ±e₁ … ±e₇ are closed under multiplication, with an identity and two-sided inverses. They satisfy all three Moufang identities on every one of the 4,096 triples, so they form a Moufang loop. Associativity holds on 2,752 triples, and the other 1,344 give opposite signs: (e₁e₂)e₃ = −e₆, while e₁(e₂e₃) = +e₆. The units form a group exactly when the algebra is associative. |
| S4 | G₂, the symmetry group of the octonions, takes the role the finite groups play on the quaternion side. | **An analogy, made exact** | See the correspondence table below. G₂ = Aut(O) (Cartan 1914, NAMED) pairs with SO(3) = Aut(H). The finite groups pair with finite subloops of S⁷ and finite subgroups of G₂. |
| S5 | 2I is the largest group on the ladder. | **Corrected** ([C1](CORRECTIONS.hbp)) | 2I, of order 120, is the largest finite subgroup of SU(2) outside every conjugate of the normaliser of a maximal torus. Every finite subgroup of order above 120 is cyclic or binary dihedral (Klein 1884); Dic₃₁, of order 124, is measured in Z[ζ₆₂]. The ladder branches at 2T: 2O sits beside 2I, since 48 does not divide 120, and 2O ∩ 2I = 2T. |
| S6 | No tower of finite subgroups reaches SU(N); it can only approach a torus extended by a finite group. | **Proven, from named theorems and measured inputs** | Every Hausdorff limit of finite subgroups of SU(N), N ≥ 2, has a torus as its identity component. This follows from Jordan's theorem together with Cartan's closed-subgroup theorem ([C4](CORRECTIONS.hbp)), and the result is classical: Turing 1938, "Finite approximations to Lie groups", *Ann. of Math.* 39, 105–111 (NAMED). Measured Jordan indices: Q8 2, 2T 12, 2O 24, 2I 60, Dic₃₁ 2. |

## The exact correspondence behind S4

| Quaternion side | Octonion side |
|---|---|
| Unit quaternions S³ = SU(2), a group | Unit octonions S⁷, a Moufang loop |
| Aut(H) = SO(3) | Aut(O) = G₂ |
| Finite subgroups of S³: Q8, 2T, 2O, 2I | Finite subloops of S⁷: the 16 units; the 240 units of Coxeter's integral octonions |
| Finite subgroups of SO(3): V4, A4, S4, A5, sizes 4, 12, 24, 60, measured as the conjugation images of the four groups | Finite subgroups of G₂: the frame group 2³·GL(3,2) of order 1,344, among others |

The frame group is measured directly. Of the 645,120 signed permutations of e₁ … e₇, exactly 1,344 preserve the multiplication table, and all of them are algebra automorphisms, so they lie in G₂. Its sign-change kernel has order 8. Its image acts on the 7 points of the Fano plane as the full group of 168 collineations, GL(3,2), giving 1,344 = 8 × 168.

The 16 units contain 7 copies of Q8, one for each line of the Fano plane.

## Files

| Folder | Script | Receipt | sha256 |
|---|---|---|---|
| [`m4/`](m4/) | `m4_quaternions.py` | `M4-QUATERNIONS.hbp`: 95 checks, 9 controls, 0 failures | `5add13fe…` |
| [`m5/`](m5/) | `m5_octonion_loop.py` | `M5-OCTONION-LOOP.hbp` | `dbf4e06e…` |
| [`refute/`](refute/) | `refute2_m4m5.py` | `REFUTE2.hbp`: 35 checks, 9 controls, 0 failures; `REFUTE2-EXTERNAL.hbp`, the refuter's run of the repository gate | `6972ff10…` |

Each script is integer only and passes [`../../tools/no_float_gate.py`](../../tools/no_float_gate.py). Each one reproduces its receipt byte for byte from a fresh directory, and the proofs workflow replays all three.

These statements are classical mathematics. What this folder adds is exact, reproducible certificates of them, with the scope of each one stated.
