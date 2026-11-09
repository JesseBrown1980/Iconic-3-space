"""predict_spinor_symbol.py  --  PREDICTION, written before any Clifford matrix is built.

What it freezes.  The measurement that follows builds integer Clifford matrices from the quaternions (4D) and from
the octonions (7D and 8D) and checks, for every covector xi, the symbol identities of three deformation sequences:
the anti-self-dual sequence on R^4, the G2-instanton sequence on R^7 and the Spin(7)-instanton sequence on R^8.
This script fixes, before that run, the conventions the measurement must use and the identities and numbers it
predicts, and it adds the symbol-versus-operator obstruction and the 4D tensor claim 16 = 4 x 4.

Hard rule kept here: no Clifford matrix, no quaternion or octonion multiplication table and no symbol matrix is
built.  The script states and counts.  Its only computations are integer arithmetic on dimensions, binomial
coefficients, Euler counts, ranks solved from dimensions, traces read from eigenvalue multiplicities, and the
Clifford-algebra table of Bott periodicity read as (summands, matrix size, real dimension of the division algebra).
Those self-checks can fail, and each planted wrong dimension must break at least one of them.

Tags.  MEASURED: computed by this script (the self-checks only).  NAMED: stated in the literature, not run here.
DERIVED: derived by hand from NAMED identities in this file's rows, not run here.  CONJECTURE: none is made.
Every identity row carries status PREDICTION and is to be tested by the measurement.

Sources (NAMED): Lawson and Michelsohn 1989, Spin Geometry (Clifford table and Bott periodicity, ch. I sec. 4 and 5);
Harvey 1990, Spinors and Calibrations (octonions, Spin(7), G2, triality); Atiyah, Hitchin and Singer 1978,
Self-duality in four-dimensional Riemannian geometry, Proc. R. Soc. Lond. A 362 (the ASD deformation complex);
Donaldson and Thomas 1998, Gauge theory in higher dimensions, in The Geometric Universe, Oxford (G2 and Spin(7)
instantons); Reyes Carrion 1998, A generalization of the notion of instanton, Diff. Geom. Appl. 8 (the G2 complex);
Baez 2002, The Octonions, Bull. AMS 39 (alternativity, Artin's theorem, Spin(7) and G2); Joyce 2000, Compact
Manifolds with Special Holonomy, ch. 10 (Spin(7) and G2 forms); Bryant 2005, Some remarks on G2-structures.

Integer discipline: integers only, no float, no division operator, no power operator, no fractions, decimal or json;
divmod with a checked remainder where a quotient is needed.  The text writes "one third" and "half" in words.

Run:  python3 -I -B predict_spinor_symbol.py
It writes PREDICTION-SPINOR-SYMBOL.hbp and PREDICTION-SPINOR-SYMBOL.hbp.sha256 beside this file and exits 0 only
when every self-check passes (nonzero otherwise).
"""

import hashlib
import os
import sys

SCHEMA = "predict-spinor-symbol-v1"
DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
SCRIPT_NAME = "predict_spinor_symbol.py"
HBP_NAME = "PREDICTION-SPINOR-SYMBOL.hbp"
SIDECAR_NAME = "PREDICTION-SPINOR-SYMBOL.hbp.sha256"

# ---------------------------------------------------------------- integer helpers


def ipow(b, e):
    if e < 0:
        raise ValueError("negative exponent")
    out = 1
    for _ in range(e):
        out = out * b
    return out


def binom(n, k):
    if k < 0 or k > n:
        return 0
    num = 1
    den = 1
    for i in range(k):
        num = num * (n - i)
        den = den * (i + 1)
    q, r = divmod(num, den)
    if r:
        raise ValueError("binomial did not divide")
    return q


def join(xs, sep=","):
    return sep.join(str(x) for x in xs)


def clean(v):
    s = str(v)
    if "|" in s or "\n" in s or "\r" in s:
        raise ValueError("a field holds a row separator: " + s)
    return s


def row(tag, fields):
    parts = [tag]
    for k, v in fields:
        parts.append(clean(k) + "=" + clean(v))
    parts.append("json=0")
    return "|".join(parts)


# ---------------------------------------------------------------- the predicted dimensions (stated, then checked)
# Suffix 4 = the quaternion case on R^4, 7 = the G2 case on R^7, 8 = the Spin(7) case on R^8.
# small and big are the two summands of Lambda2: 4D Lambda2+ and Lambda2-, 7D Lambda2_7 and Lambda2_14,
# 8D Lambda2_7 and Lambda2_21.  T is the splitting operator: 4D *, 7D *(phi wedge .), 8D *(Phi wedge .).

DIMS = {
    "n_4": 4, "n_7": 7, "n_8": 8,
    "S_plus_4": 4, "S_minus_4": 4, "S_C_half_4": 2,
    "S_7": 8,
    "S_plus_8": 8, "S_minus_8": 8,
    "small_4": 3, "big_4": 3,
    "small_7": 7, "big_7": 14,
    "small_8": 7, "big_8": 21,
    "T_small_4": 1, "T_small_7": 2, "T_small_8": 3, "T_big": -1,
    "k_4": 2, "k_7": 3, "k_8": 4,
    "stab_4": 3, "stab_7": 14, "stab_8": 21,
    "so4": 6, "sp1": 3, "g2": 14, "spin7": 21,
    "seq_4": (1, 4, 3), "seq_7": (1, 7, 7, 1), "seq_8": (1, 8, 7),
    "ranks_4": (1, 3), "ranks_7": (1, 6, 1), "ranks_8": (1, 7),
    "cl_4": 16, "cl_7": 128, "cl_8": 256,
    "alg_4": 16, "alg_7": 64, "alg_8": 256,
    "comm_4": 4, "comm_7": 1, "comm_8": 1,
    "rep_4": 8, "rep_7": 8, "rep_8": 16,
    "gens_4": 4, "gens_7": 7, "gens_8": 8,
    "coef_4": 640, "coef_7": 1792, "coef_8": 9216,
    "spin_C_4": 4, "spin_C_8": 16,
    "g2_meet": 14,
    "su2": 3,
    "fano_points": 7, "fano_lines": 7, "fano_line_size": 3, "fano_lines_per_point": 3,
    "phi_ordered_nonzero": 42, "assoc_pairs": 42, "assoc_rank": 4,
    "L3_split_7": (1, 7, 27),
}

CASES = (("4D", "4"), ("7D", "7"), ("8D", "8"))

# Clifford algebras by Bott periodicity, n = 0..8, as (summands, matrix size m, real dimension d of the division algebra).
# CL_NEG: generators square to -1 (Lawson-Michelsohn Cl_n).  CL_POS: generators square to +1.  NAMED.
CL_NEG = ((1, 1, 1), (1, 1, 2), (1, 1, 4), (2, 1, 4), (1, 2, 4), (1, 4, 2), (1, 8, 1), (2, 8, 1), (1, 16, 1))
CL_NEG_NAMES = ("R", "C", "H", "H+H", "M_2(H)", "M_4(C)", "M_8(R)", "M_8(R)+M_8(R)", "M_16(R)")
CL_POS = ((1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 2, 2), (1, 2, 4), (2, 2, 4), (1, 4, 4), (1, 8, 2), (1, 16, 1))
CL_POS_NAMES = ("R", "R+R", "M_2(R)", "M_2(C)", "M_2(H)", "M_2(H)+M_2(H)", "M_4(H)", "M_8(C)", "M_16(R)")

# Planted wrong dimensions: each must break at least one self-check of the same suite.
PLANTS = (
    ("lambda2_plus_4d_planted_as_4", "small_4", 4),
    ("lambda2_14_g2_planted_as_15", "big_7", 15),
    ("lambda2_7_spin7_planted_as_8", "small_8", 8),
    ("cl8_planted_as_128", "cl_8", 128),
    ("spinor_plus_8d_planted_as_9", "S_plus_8", 9),
    ("ranks_7d_planted_as_1_5_1", "ranks_7", (1, 5, 1)),
)

# What the task text states, compared (never forced) with the table above.
TASK_EXPECT = (
    ("4D_sequence_dims", (1, 4, 3), "seq_4"),
    ("4D_ranks", (1, 3), "ranks_4"),
    ("7D_sequence_dims", (1, 7, 7, 1), "seq_7"),
    ("7D_ranks", (1, 6, 1), "ranks_7"),
    ("8D_sequence_dims", (1, 8, 7), "seq_8"),
    ("8D_ranks", (1, 7), "ranks_8"),
    ("4D_lambda2_plus", 3, "small_4"),
    ("7D_lambda2_7", 7, "small_7"),
    ("7D_lambda2_14", 14, "big_7"),
    ("8D_lambda2_7", 7, "small_8"),
    ("8D_lambda2_21", 21, "big_8"),
    ("8D_spinors_8_plus_8", 16, "rep_8"),
    ("8D_algebra_dim", 256, "alg_8"),
    ("tensor_16", 16, "spin_C_8"),
    ("8D_generator_count_task_text_16", 16, "gens_8"),
)


# ---------------------------------------------------------------- the self-check suite


def solve_ranks(seq):
    # For an exact sequence 0 -> d0 -> d1 -> ... -> 0: r1 = d0, r(m+1) = d(m) - r(m); the last residual must be 0.
    ranks = []
    prev = 0
    for d in seq:
        r = d - prev
        ranks.append(r)
        prev = r
    residual = ranks[-1]
    return tuple(ranks[:-1]), residual


def analyze(D):
    out = []

    def chk(name, ok, detail):
        out.append((name, bool(ok), detail))

    half = {"4": D["S_plus_4"], "7": D["S_7"], "8": D["S_plus_8"]}
    for tag, s in CASES:
        n = D["n_" + s]
        small, big = D["small_" + s], D["big_" + s]
        k, ts, tb = D["k_" + s], D["T_small_" + s], D["T_big"]
        c2 = binom(n, 2)
        chk("lambda2_split_" + tag, small + big == c2, [("small", small), ("big", big), ("binom_n_2", c2)])
        chk("trace_T_zero_" + tag, ts * small + tb * big == 0, [("trace", ts * small + tb * big)])
        chk("k_times_small_is_binom_n_2_" + tag, k * small == c2, [("k", k), ("k_times_small", k * small)])
        chk("k_is_T_small_plus_1_" + tag, k == ts + 1, [("k", k), ("T_small", ts)])
        quad = all(t * t == (k - 2) * t + (k - 1) for t in (ts, tb))
        chk("T_eigenvalues_solve_T2_eq_k_minus_2_T_plus_k_minus_1_" + tag, quad, [("eigenvalues", join((ts, tb)))])
        chk("stabilizer_is_big_summand_" + tag, D["stab_" + s] == big, [("stab", D["stab_" + s]), ("big", big)])
        chk("orbit_is_sphere_of_S_plus_" + tag, c2 - D["stab_" + s] == half[s] - 1,
            [("group_minus_stab", c2 - D["stab_" + s]), ("dim_S_plus_minus_1", half[s] - 1)])
        chk("rank_Psi_is_small_" + tag, small == half[s] - 1, [("small", small), ("dim_S_plus_minus_1", half[s] - 1)])
        seq = D["seq_" + s]
        alt = 0
        for i, d in enumerate(seq):
            q, r = divmod(i, 2)
            alt = alt + (d if r == 0 else -d)
        chk("euler_count_zero_" + tag, alt == 0, [("seq", join(seq)), ("euler", alt)])
        ranks, residual = solve_ranks(seq)
        ok = residual == 0 and all(r >= 0 for r in ranks) and ranks == tuple(D["ranks_" + s])
        chk("exact_ranks_from_dims_" + tag, ok, [("solved", join(ranks)), ("residual", residual),
                                                ("predicted", join(D["ranks_" + s]))])
        even = 0
        odd = 0
        for i, d in enumerate(seq):
            q, r = divmod(i, 2)
            if r == 0:
                even = even + d
            else:
                odd = odd + d
        chk("rolled_up_square_and_spinor_sized_" + tag, even == odd == half[s],
            [("even_sum", even), ("odd_sum", odd), ("dim_S_half", half[s])])
        forms_ok = seq[0] == binom(n, 0) and seq[1] == n and seq[2] == small
        if s == "7":
            forms_ok = forms_ok and len(seq) == 4 and seq[3] == 1
        else:
            forms_ok = forms_ok and len(seq) == 3
        chk("sequence_terms_are_the_forms_" + tag, forms_ok, [("seq", join(seq))])
        cl = D["cl_" + s]
        total = 0
        for j in range(n + 1):
            total = total + binom(n, j)
        summ, m, dd = CL_NEG[n]
        chk("clifford_dim_2_pow_n_" + tag, cl == ipow(2, n) == total == summ * m * m * dd,
            [("cl", cl), ("two_pow_n", ipow(2, n)), ("sum_binom", total), ("bott_entry", CL_NEG_NAMES[n])])
        alg, comm, rep = D["alg_" + s], D["comm_" + s], D["rep_" + s]
        ok = alg == m * m * dd and comm == dd and rep == m * dd and alg * comm == rep * rep and summ * alg == cl
        chk("generated_algebra_and_commutant_" + tag, ok,
            [("alg", alg), ("commutant", comm), ("module", rep), ("alg_times_commutant", alg * comm),
             ("module_squared", rep * rep), ("summands", summ)])
        chk("generator_count_is_n_" + tag, D["gens_" + s] == n, [("gens", D["gens_" + s])])
        pairs = binom(n, 2) + n
        chk("anticommutator_count_" + tag, pairs == binom(n + 1, 2) and D["coef_" + s] == rep * rep * pairs,
            [("pairs", pairs), ("monomials_deg2", binom(n + 1, 2)), ("coefficients", rep * rep * pairs)])

    # Bott tables, both signs, n = 0..8, and the period 8.
    ok = True
    for n in range(9):
        a, m, d = CL_NEG[n]
        b, m2, d2 = CL_POS[n]
        ok = ok and a * m * m * d == ipow(2, n) and b * m2 * m2 * d2 == ipow(2, n)
    ok = ok and CL_NEG[8] == (CL_NEG[0][0], 16 * CL_NEG[0][1], CL_NEG[0][2])
    ok = ok and CL_POS[8] == (CL_POS[0][0], 16 * CL_POS[0][1], CL_POS[0][2])
    chk("bott_tables_dim_2_pow_n_and_period_8", ok, [("n_range", "0..8")])
    chk("cl7_two_simple_factors_image_is_one", D["cl_7"] == 2 * D["alg_7"] and D["alg_7"] == D["rep_7"] * D["rep_7"],
        [("cl7", D["cl_7"]), ("alg7", D["alg_7"])])

    # Spinor decompositions.
    h = D["S_C_half_4"]
    chk("4D_SplusSplus_is_1_plus_lambda2plus_complex", h * h == 1 + D["small_4"], [("two_times_two", h * h)])
    chk("4D_SplusSminus_is_lambda1_complex", h * h == D["n_4"], [("two_times_two", h * h)])
    chk("4D_real_spinors_8", D["S_plus_4"] == 2 * h and D["S_plus_4"] + D["S_minus_4"] == D["rep_4"],
        [("S_plus", D["S_plus_4"]), ("S_minus", D["S_minus_4"])])
    chk("7D_S_is_lambda0_plus_lambda1", 1 + D["n_7"] == D["S_7"] == D["rep_7"], [("one_plus_seven", 1 + D["n_7"])])
    chk("8D_S_plus_is_lambda0_plus_lambda2_7", 1 + D["small_8"] == D["S_plus_8"], [("one_plus_seven", 1 + D["small_8"])])
    chk("8D_S_minus_is_lambda1", D["n_8"] == D["S_minus_8"], [("S_minus", D["S_minus_8"])])
    chk("8D_spinors_16", D["S_plus_8"] + D["S_minus_8"] == D["rep_8"] and D["cl_8"] == D["rep_8"] * D["rep_8"],
        [("sum", D["S_plus_8"] + D["S_minus_8"]), ("cl8", D["cl_8"])])

    # Lie algebras.
    chk("lie_dims", D["g2"] == D["big_7"] and D["spin7"] == D["big_8"] == binom(7, 2)
        and D["so4"] == binom(4, 2) == 2 * D["sp1"] and D["big_4"] == D["sp1"] == D["small_4"],
        [("g2", D["g2"]), ("spin7", D["spin7"]), ("so4", D["so4"]), ("sp1", D["sp1"])])
    chk("spin7_plus_meet_spin7_minus_is_g2", D["big_8"] + D["big_8"] - D["g2_meet"] == binom(8, 2)
        and D["g2_meet"] == D["g2"], [("sum_minus_meet", D["big_8"] + D["big_8"] - D["g2_meet"])])

    # The tensor claim 16 = 4 x 4.
    c4, c8 = D["spin_C_4"], D["spin_C_8"]
    forms4 = 0
    for j in range(5):
        forms4 = forms4 + binom(4, j)
    chk("tensor_16_is_4_times_4", c8 == c4 * c4 == ipow(2, 4) == forms4 and c4 == ipow(2, 2) == 2 * h,
        [("spin_C_8", c8), ("four_times_four", c4 * c4), ("forms_on_R4", forms4)])
    blocks = (1 + D["small_4"]) + (1 + D["big_4"]) + D["n_4"] + D["n_4"]
    chk("tensor_16_by_chirality_blocks", blocks == c8, [("blocks", "(1+3)+(1+3)+4+4"), ("total", blocks)])
    chk("tensor_real_form_is_8D_real_spinor", c8 == D["rep_8"], [("complex_dim", c8), ("real_dim_8D", D["rep_8"])])

    # Octonion incidence (counts only, no table).
    p, ln, sz, lpp = D["fano_points"], D["fano_lines"], D["fano_line_size"], D["fano_lines_per_point"]
    ok = p * lpp == ln * sz and binom(p, 2) == ln * binom(sz, 2) and D["phi_ordered_nonzero"] == ln * 6
    ok = ok and D["assoc_pairs"] == p * (p - 1) and D["assoc_rank"] == D["S_7"] - 4 and p == D["n_7"]
    chk("fano_counts", ok, [("incidences", p * lpp), ("pairs", binom(p, 2)), ("phi_ordered", ln * 6),
                            ("ordered_pairs", p * (p - 1)), ("associator_rank", D["assoc_rank"])])
    chk("7D_lambda3_split_35", sum(D["L3_split_7"]) == binom(7, 3), [("split", join(D["L3_split_7"]))])

    # Obstruction: the su(2)-valued ASD symbol sequence.
    sseq = tuple(D["su2"] * d for d in D["seq_4"])
    sranks, sres = solve_ranks(sseq)
    chk("su2_valued_ASD_symbol_exact", sres == 0 and sranks == tuple(D["su2"] * r for r in D["ranks_4"]),
        [("dims", join(sseq)), ("ranks", join(sranks))])
    return out


# ---------------------------------------------------------------- rows


def conv_rows():
    C = (
        ("C0.1", "all", "clifford_sign", "gamma(u) gamma(v) + gamma(v) gamma(u) = -2 <u, v> I, so gamma(xi)^2 = -nsq(xi) I "
         "(the Lawson-Michelsohn sign: Cl_n = Cl(R^n, -nsq), generators square to -1)", "NAMED"),
        ("C0.2", "all", "symbol", "principal symbols with the factor i dropped: sigma(d)(xi) = xi wedge ., "
         "sigma(d*)(xi) = -iota_xi, sigma(pi d)(xi) = pi(xi wedge .) for a projection pi", "NAMED"),
        ("C0.3", "all", "every_covector", "an identity in xi is checked as an identity of polynomials with integer "
         "coefficients in the coordinates of xi, entry by entry and monomial by monomial; equal coefficients mean "
         "equality for every real covector xi, in finitely many exact integer checks; sampled xi are a control only", "NAMED"),
        ("C0.4", "all", "rank_certificate", "a rank claim for every xi != 0 is certified by a polynomial identity, "
         "such as C(xi) C(xi)^T = nsq(xi) I, never by sampled ranks alone", "NAMED"),
        ("C0.5", "all", "rotation", "R_ij e_i = e_j, R_ij e_j = -e_i, R_ij = 0 on the other basis vectors; for "
         "omega = sum_(i<j) omega_ij e_ij, R_omega = sum omega_ij R_ij and c(omega) = sum_(i<j) omega_ij gamma(e_i) gamma(e_j); "
         "then [c(omega), gamma(v)] = 2 gamma(R_omega v) by the Clifford relations alone", "DERIVED"),
        ("C0.6", "all", "Psi", "Psi(omega) = c(omega) eta, the Clifford action of a 2-form on the fixed spinor eta = 1 "
         "(in S+ for 4D and 8D, in S for 7D); its kernel is the Lie algebra of the stabilizer of eta", "NAMED"),
        ("C0.7", "all", "splitting_operator", "T on Lambda2: 4D T = *, 7D T = *(phi wedge .), 8D T = *(Phi wedge .); "
         "T = k - 1 on the small summand and -1 on the big one, with k = 2, 3, 4; the projection to the small summand is "
         "(T + I) times one over k, and (T + I) is the integer form the measurement uses", "NAMED"),
        ("C0.8", "all", "integer_form", "every identity is stated with integer coefficients on both sides: the factors "
         "sqrt2, sqrt3 and 2 of the isometric normalizations and the projections are cleared as each row states", "DERIVED"),
        ("C0.9", "all", "forms_inner_product", "e_ij = e_i wedge e_j with i < j; the e_I are orthonormal; "
         "<xi wedge a, e_ij> = xi_i a_j - xi_j a_i; nsq(xi wedge a) = nsq(xi) nsq(a) - <xi, a>^2", "NAMED"),
        ("4D.C1", "4D", "basis", "R^4 = H with e1 = 1, e2 = i, e3 = j, e4 = k, Hamilton's i j = k, j k = i, k i = j, "
         "i^2 = j^2 = k^2 = -1; conj the quaternion conjugate; L_x y = x y; orientation e1 wedge e2 wedge e3 wedge e4", "NAMED"),
        ("4D.C2", "4D", "gamma", "S = S+ + S- = H + H, real dimension 4 + 4 = 8; gamma(xi) = [[0, L_xi],[-L_conj(xi), 0]] "
         "for xi in H, S+ the first summand: gamma(xi) maps S- to S+ by L_xi and S+ to S- by -L_conj(xi)", "NAMED"),
        ("4D.C3", "4D", "spin_group", "Spin(4) = Sp(1) x Sp(1): (p, q) acts on R^4 by xi -> p xi conj(q), on S+ by s -> p s, "
         "on S- by s -> q s; gamma(p xi conj(q)) = g gamma(xi) g^-1 with g = diag(L_p, L_q), since conj(p xi conj(q)) = "
         "q conj(xi) conj(p)", "DERIVED"),
        ("4D.C4", "4D", "lambda2_bases", "Lambda2+ basis f1 = e12 + e34, f2 = e13 - e24, f3 = e14 + e23 (self-dual for "
         "e1234: *e12 = e34, *e13 = -e24, *e14 = e23), Gram matrix G = 2 I_3; Lambda2- basis e12 - e34, e13 + e24, "
         "e14 - e23", "NAMED"),
        ("4D.C5", "4D", "c_coordinates", "c_k(xi, a) = <xi wedge a, f_k>: c1 = (xi1 a2 - xi2 a1) + (xi3 a4 - xi4 a3), "
         "c2 = (xi1 a3 - xi3 a1) - (xi2 a4 - xi4 a2), c3 = (xi1 a4 - xi4 a1) + (xi2 a3 - xi3 a2); then "
         "(I + *)(xi wedge a) = 2 (xi wedge a)^+ = c1 f1 + c2 f2 + c3 f3, and c is the coordinate vector of "
         "sqrt2 (xi wedge a)^+ in the orthonormal basis f_k times one over sqrt2", "DERIVED"),
        ("4D.C6", "4D", "M_matrix", "M(xi) is the 4 x 4 integer matrix with row 0 = -xi^T (the -iota_xi part) and row k = "
         "(c_k(xi, e1), c_k(xi, e2), c_k(xi, e3), c_k(xi, e4)) for k = 1, 2, 3; it is the matrix of "
         "sigma(xi) = (-iota_xi, sqrt2 (xi wedge .)^+) : Lambda1 -> Lambda0 + Lambda2+ in orthonormal bases", "DERIVED"),
        ("7D.C1", "7D", "table", "Im O with orthonormal basis e1..e7 from the measurement's octonion table: any alternative "
         "table with e_i^2 = -1 and e_i e_j = -e_j e_i = plus or minus e_k along the 7 lines of a Fano plane; "
         "O = R 1 + Im O; conj fixes 1 and negates Im O; L_x y = x y", "NAMED"),
        ("7D.C2", "7D", "phi_and_cross", "phi(x, y, z) = <x y, z> on Im O, phi_ijk = <e_i e_j, e_k>, totally antisymmetric; "
         "x cross y = Im(x y), (x cross y)_k = sum_ij phi_ijk x_i y_j; (iota_u phi)_ij = sum_k u_k phi_kij; the contraction "
         "of a 2-form omega with phi is the vector w with w_k = sum_(i<j) omega_ij phi_ijk", "NAMED"),
        ("7D.C3", "7D", "orientation", "orient R^7 so that T7 = *(phi wedge .) on Lambda2 has eigenvalue +2 on a "
         "7-dimensional and -1 on a 14-dimensional subspace, the split proofs g2-s6 measured (M3 T_eigen_split_7_14, REFUTE "
         "lambda2_split_7_14); psi = *phi in that orientation; if e1 wedge ... wedge e7 gives -2 and +1, the measurement "
         "reverses it and says so", "NAMED"),
        ("7D.C4", "7D", "spinor", "S = O, real dimension 8 = Lambda0 + Lambda1 by (f, v) <-> f 1 + v; gamma(xi) = L_xi for "
         "xi in Im O; eta = 1", "NAMED"),
        ("7D.C5", "7D", "coordinates", "the Lambda2_7 coordinate of a 2-form omega is w = Psi_7(omega) = the contraction of "
         "omega with phi (integer), and then pi_7 omega = one third of iota_w phi; the Lambda3_1 coordinate of the 3-form "
         "xi wedge pi_7 omega is t = <xi wedge pi_7 omega, phi>", "DERIVED"),
        ("8D.C1", "8D", "basis", "O = R^8 with e0 = 1 and e1..e7 the imaginary units of the same table as 7D.C1; "
         "orientation vol8 = e0 wedge vol7 with vol7 from 7D.C3; e0j = e0 wedge e_j; psi = *7 phi", "NAMED"),
        ("8D.C2", "8D", "gamma", "S+ + S- = O + O, real dimension 8 + 8 = 16; gamma(xi) = [[0, L_xi],[-L_conj(xi), 0]] "
         "for xi in O, S+ the first summand; eta = 1 in S+, eta' = 1 in S-; gamma(e_i) gamma(e_j) = "
         "diag(-L_(e_i) L_conj(e_j), -L_conj(e_i) L_(e_j))", "DERIVED"),
        ("8D.C3", "8D", "cayley_forms", "Phi_s = e0 wedge phi + s psi for s = +1 and s = -1; T_s = *8(Phi_s wedge .) on "
         "Lambda2; *8(e0 wedge alpha) = *7 alpha, so Phi_+1 is self-dual and Phi_-1 anti-self-dual for vol8", "DERIVED"),
        ("8D.C4", "8D", "Psi_maps", "Psi(omega) = c(omega) eta in S+ and Psi'(omega) = c(omega) eta' in S-", "NAMED"),
        ("OB.C1", "obstruction", "pauli", "sigma1 = [[0, 1],[1, 0]], sigma2 = [[0, -i],[i, 0]], sigma3 = [[1, 0],[0, -1]], "
         "entries Gaussian integers written re:im; [sigma_a, sigma_b] = 2 i eps_abc sigma_c; su(2) is spanned by "
         "i sigma1, i sigma2, i sigma3", "NAMED"),
        ("OB.C2", "obstruction", "connection", "d_A = d + [A, .] on su(2)-valued forms, A = sum_i A_i dx^i; "
         "F_A = dA + A wedge A; for constant A, F_A = A wedge A = sum_(i<j) [A_i, A_j] dx^i wedge dx^j and "
         "d_A d_A = [F_A, .]", "NAMED"),
    )
    rows = []
    for cid, case, name, value, claim in C:
        rows.append(row("CONV", [("id", cid), ("case", case), ("name", name), ("value", value),
                                 ("status", "CONVENTION"), ("claim", claim)]))
    return rows


def pred_rows():
    P = (
        # ---------------- 4D
        ("4D.P1", "4D", "gamma(xi)^2 = -nsq(xi) I_8 as a polynomial identity in xi1..xi4",
         "all 640 coefficients of gamma(xi)^2 + nsq(xi) I_8 are 0 (64 entries x 10 monomials of degree 2)",
         "H is associative and x conj(x) = conj(x) x = nsq(x)", "NAMED", "Lawson-Michelsohn I.4"),
        ("4D.P2", "4D", "gamma(e_i) gamma(e_j) + gamma(e_j) gamma(e_i) = -2 delta_ij I_8",
         "10 of 10 relations i <= j hold", "the polarization of 4D.P1", "NAMED", "Lawson-Michelsohn I.1"),
        ("4D.P3", "4D", "the 4 matrices gamma(e1..e4), each 8 x 8 with entries 0, 1, -1, generate a real algebra",
         "dimension 16 = 2^4 = dim Cl_4, Cl_4 = M_2(H), faithful because Cl_4 is simple; commutant in M_8(R) of "
         "dimension 4 (right multiplication by H on both summands); 16 x 4 = 64 = 8^2",
         "Bott table entry n = 4", "NAMED", "Lawson-Michelsohn I.4"),
        ("4D.P4", "4D", "omega4 = gamma(e1) gamma(e2) gamma(e3) gamma(e4)",
         "omega4 = epsilon4 diag(I_4, -I_4) with epsilon4 = -1; omega4^2 = I_8; omega4 anticommutes with every gamma(xi)",
         "block (1,1) = L_(1 conj(i) j conj(k)) = L_(i j k) = L_(-1); block (2,2) = L_(conj(1) i conj(j) k) = L_(1)",
         "DERIVED", "4D.C1, 4D.C2"),
        ("4D.P5", "4D", "ASD sequence Lambda0 -> Lambda1 -> Lambda2+, sigma0(xi) f = f xi, sigma1(xi) a = (xi wedge a)^+",
         "dimensions 1, 4, 3; Euler count 1 - 4 + 3 = 0; sigma1(xi) sigma0(xi) = (xi wedge xi)^+ = 0 as a polynomial "
         "identity; for every xi != 0 rank sigma0 = 1, rank sigma1 = 3 and ker sigma1 = span(xi) = im sigma0 (exact)",
         "the symbol complex of d, d+ (elliptic)", "NAMED", "Atiyah-Hitchin-Singer 1978"),
        ("4D.P6", "4D", "C(xi) = the 3 x 4 matrix of the c-rows of M(xi)",
         "C(xi) xi = 0 and C(xi) C(xi)^T = nsq(xi) I_3 as polynomial identities, so det(C C^T) = nsq(xi)^3 and rank 3 for "
         "every real xi != 0; the signed maximal minors of C(xi) form lambda xi with lambda = plus or minus nsq(xi)",
         "M(xi) M(xi)^T = nsq(xi) I_4 restricted to rows 1..3; Cauchy-Binet", "DERIVED", "4D.P7"),
        ("4D.P7", "4D", "rolled-up symbol sigma(xi) = (-iota_xi, sqrt2 (xi wedge .)^+) : Lambda1 -> Lambda0 + Lambda2+",
         "M(xi)^T M(xi) = M(xi) M(xi)^T = nsq(xi) I_4 as polynomial identities; equivalently, with sigma_u(xi) = "
         "(-iota_xi, c(xi, .)) the f-coordinates of the unscaled projection (I + *)(xi wedge .), Gram G = 2 I_3 on "
         "Lambda2+ and weight 2 on Lambda0 (G' = diag(2, G) = 2 I_4): sigma_u(xi)^T G' sigma_u(xi) = 2 nsq(xi) I_4",
         "nsq((xi wedge a)^+) = half of nsq(xi wedge a) for a decomposable 2-form, because xi wedge a wedge xi wedge a = 0",
         "DERIVED", "4D.C5, 4D.C6, C0.9"),
        ("4D.P8", "4D", "Psi_4(e_ij) = -e_i conj(e_j) for i < j (gamma_i gamma_j on S+ is -L_(e_i) L_conj(e_j))",
         "Psi_4(e12) = i, Psi_4(e13) = j, Psi_4(e14) = k, Psi_4(e23) = k, Psi_4(e24) = -j, Psi_4(e34) = i; so "
         "Psi_4(f1) = 2i, Psi_4(f2) = 2j, Psi_4(f3) = 2k and Psi_4 = 0 on Lambda2-; Psi_4 Psi_4^T = 2 I_3; "
         "Psi_4 (* + I) = 2 Psi_4; c(omega) on S+ equals L_(Psi_4(omega)); ker Psi_4 = Lambda2- = sp(1)_-, dimension 3",
         "hand evaluation with Hamilton's rules", "DERIVED", "4D.C1, 4D.C2, C0.6"),
        ("4D.P9", "4D", "the rolled-up symbol is Clifford multiplication S- -> S+",
         "with Lambda1 -> S-, a -> c(a) eta = -conj(a), and Lambda0 + Lambda2+ -> S+, (f, omega) -> f + Psi_4(omega), the "
         "unnormalized symbol (-iota_xi, (xi wedge .)^+) is L_xi: (M(xi) a)_0 + (M(xi) a)_1 i + (M(xi) a)_2 j + "
         "(M(xi) a)_3 k = -xi conj(a), a bilinear identity (4 components x 16 monomials xi_i a_j = 64 integer "
         "coefficients); tensor form: Lambda1 = S+ (x) S- by (u, v) -> u conj(v), Lambda0 + Lambda2+ = S+ (x) S+ by "
         "(u, w) -> u conj(w), sigma = 1 (x) Clifford multiplication on the second factor, a -> a conj(xi), and with "
         "(s; c) -> -s + c1 i + c2 j + c3 k, M(xi) a goes to a conj(xi)",
         "gamma(xi) gamma(a) = -<xi, a> + c(xi wedge a) in the Clifford algebra, and Psi_4 kills Lambda2-", "DERIVED",
         "4D.P8"),
        ("4D.P10", "4D", "representations",
         "complex: S+ (x) S+ = Lambda2 S+ + Sym2 S+ = 1 + 3 = (Lambda0 + Lambda2+) (x) C and S+ (x) S- = 2 x 2 = 4 = "
         "Lambda1 (x) C; real: Lambda0 + Lambda2+ = H with y -> p y conj(p), Lambda1 = H with a -> p a conj(q); "
         "Lambda2+ = sp(1)_+ acts on S+ and kills S-, Lambda2- = sp(1)_- kills S+",
         "Spin(4) = Sp(1) x Sp(1)", "NAMED", "Atiyah-Hitchin-Singer 1978; Lawson-Michelsohn IV"),
        ("4D.P11", "4D", "equivariance under Sp(1) x Sp(1)",
         "for the 6 rotations R_ij of so(4) = sp(1)_+ + sp(1)_-: sigma(R xi) a + sigma(xi)(R a) = rho(R) sigma(xi) a, rho = 0 "
         "on Lambda0 and the derivation action on 2-forms, which preserves Lambda2+; finite form for (p, q) in Q8 x Q8 "
         "(64 pairs from plus or minus 1, i, j, k, integer signed permutation matrices): M(p xi conj(q)) = "
         "A_p M(xi) B_pq^T with B_pq a = p a conj(q) and A_p = diag(1, Ad_p), Ad_p the 3 x 3 matrix of y -> p y conj(p) on "
         "Im H in the basis i, j, k",
         "(p a conj(q)) conj(p xi conj(q)) = p (a conj(xi)) conj(p)", "DERIVED", "4D.P9"),
        # ---------------- 7D
        ("7D.P1", "7D", "L_xi^2 = -nsq(xi) I_8 as a polynomial identity in xi1..xi7",
         "all 64 x 28 = 1792 coefficients of L_xi^2 + nsq(xi) I_8 are 0",
         "it needs only left alternativity x (x y) = (x x) y with x x = -nsq(x) for x in Im O (Artin: two elements "
         "generate an associative subalgebra); also L_xi^T = -L_xi", "NAMED", "Baez 2002 sec. 2; Harvey 1990"),
        ("7D.P2", "7D", "L_(e_i) L_(e_j) + L_(e_j) L_(e_i) = -2 delta_ij I_8",
         "28 of 28 relations i <= j hold", "the polarization of 7D.P1", "NAMED", "Baez 2002"),
        ("7D.P3", "7D", "the 7 matrices L_(e1..e7), each an 8 x 8 signed permutation matrix, generate a real algebra",
         "dimension 64 = all of M_8(R) = End(O), and not 128 = dim Cl_7: Cl_7 = M_8(R) + M_8(R) has two simple factors and "
         "O carries one of its two irreducible modules; commutant dimension 1; 64 x 1 = 8^2",
         "Bott table entry n = 7", "NAMED", "Lawson-Michelsohn I.4"),
        ("7D.P4", "7D", "L_(e1) L_(e2) ... L_(e7)",
         "epsilon7 I_8 with epsilon7 = +1 or -1; CANNOT_PREDICT the sign here without the table and the ordering; the "
         "measurement reports it and 8D.P4 reuses it",
         "the volume element is central in Cl_7 with square +1 and acts as a scalar on an irreducible module", "NAMED",
         "Lawson-Michelsohn I.3"),
        ("7D.P5", "7D", "nonassociativity witness",
         "for each of the 42 ordered pairs i != j, L_(e_i e_j) - L_(e_i) L_(e_j) is nonzero with exactly 4 nonzero entries, "
         "each plus or minus 2: it is the associator y -> [e_i, e_j, y], zero on the quaternion subalgebra spanned by 1, "
         "e_i, e_j, e_i e_j and twice a signed permutation on the other 4 basis vectors; L is not multiplicative, yet the "
         "Clifford relations of 7D.P2 hold",
         "(e_i e_j) e_k = -e_i (e_j e_k) for basis units on no common line", "NAMED", "Baez 2002 sec. 2"),
        ("7D.P6", "7D", "G2-instanton sequence Lambda0 -> Lambda1 -> Lambda2_7 -> Lambda3_1 in the coordinates of 7D.C5",
         "sigma0(xi) f = f xi, sigma1(xi) v = Psi_7(xi wedge v) = xi cross v, sigma2(xi) w = <xi, w>; dimensions 1, 7, 7, 1; "
         "Euler count 1 - 7 + 7 - 1 = 0; composites xi cross xi = 0 and <xi, xi cross v> = 0 as polynomial identities; "
         "ranks 1, 6, 1 for every xi != 0",
         "<xi wedge iota_w phi, phi> = <iota_w phi, iota_xi phi> = 3 <w, xi> and pi_7 = one third of iota_w phi",
         "NAMED", "Reyes Carrion 1998; Donaldson-Thomas 1998"),
        ("7D.P7", "7D", "X(xi) v = xi cross v (the 7 x 7 integer matrix of sigma1)",
         "X(xi)^T = -X(xi) and X(xi)^2 = -nsq(xi) I_7 + xi xi^T as polynomial identities; so ker X(xi) = span(xi), "
         "im X(xi) = xi^perp, rank 6 for every xi != 0: exactness at Lambda1 and at Lambda2_7",
         "x cross (x cross y) = -nsq(x) y + <x, y> x", "NAMED", "Harvey 1990; Bryant 2005"),
        ("7D.P8", "7D", "the projection to Lambda2_7",
         "(T7 - 2 I)(T7 + I) = 0 on Lambda2 and trace T7 = 2 x 7 - 14 = 0; 3 pi_7 = T7 + I; (T7 + I)(xi wedge v) = "
         "iota_(xi cross v) phi as a bilinear identity; nsq((T7 + I)(xi wedge v)) = 3 nsq(xi wedge v); "
         "Psi_7 (T7 + I) = 3 Psi_7; Psi_7 Psi_7^T = 3 I_7; <iota_u phi, iota_v phi> = 3 <u, v>; ker Psi_7 = Lambda2_14 = "
         "der(O) as rotations, as proofs g2-s6 C7 measured",
         "pi_7(xi wedge v) = one third of iota_(xi cross v) phi", "DERIVED", "Bryant 2005; 7D.C3"),
        ("7D.P9", "7D", "rolled-up symbol on Lambda0 + Lambda1",
         "R(xi)(f, v) = (-<xi, v>, f xi + xi cross v), folding E2 = Lambda2_7 onto Lambda1 by w and E3 = Lambda3_1 onto "
         "Lambda0 by t (the symbol of (f, v) -> (d* v, d f + curl v), curl v = Psi_7(d v), the G2 Dirac operator); "
         "R(xi) = L_xi in the basis 1, e1..e7 for every xi in Im O (linear in xi: 7 checks of 8 x 8); R(xi)^2 = "
         "-nsq(xi) I_8; R(xi)^T R(xi) = nsq(xi) I_8; isometric normalization sqrt3: nsq(Psi_7(omega)) = 3 nsq(omega) on "
         "Lambda2_7",
         "L_xi(f + v) = f xi + xi v and xi v = -<xi, v> + xi cross v for xi, v in Im O", "DERIVED", "7D.P6"),
        ("7D.P10", "7D", "equivariance under G2 = Aut(O)",
         "for each D in der(O) (dimension 14, as proofs g2-s6 m1 measured): D 1 = 0, [D, L_xi] = L_(D xi), so "
         "[D, R(xi)] = R(D xi); the span of L_(e_i) L_(e_j), i < j, has dimension 21 and is closed under commutators "
         "(spin(7) on S); its part killing eta = 1 has dimension 14 and is the image of Lambda2_14: for omega in "
         "Lambda2_14, sum_(i<j) omega_ij L_(e_i) L_(e_j) = 2 R_omega on Im O and 0 on 1, a derivation",
         "D(x y) = (D x) y + x (D y); [half of c(omega), L_v] = L_(R_omega v)", "DERIVED", "C0.5, C0.6"),
        ("7D.P11", "7D", "representations of G2",
         "S = Lambda0 + Lambda1 = 1 + 7 = R + R^7 = O; Lambda2 = Lambda2_7 + Lambda2_14 = 7 + 14 with Lambda2_14 = g2; "
         "Lambda3 = 1 + 7 + 27 = 35", "G2 the stabilizer of eta = 1 in Spin(7)", "NAMED", "Bryant 2005; Harvey 1990"),
        # ---------------- 8D
        ("8D.P1", "8D", "gamma(xi)^2 = -nsq(xi) I_16 as a polynomial identity in xi0..xi7",
         "all 256 x 36 = 9216 coefficients of gamma(xi)^2 + nsq(xi) I_16 are 0",
         "it needs L_xi L_conj(xi) = L_conj(xi) L_xi = nsq(xi) I_8, from alternativity (xi, conj(xi) = 2 Re(xi) - xi and y "
         "lie in the associative subalgebra generated by xi and y) or from the composition law nsq(x y) = nsq(x) nsq(y) "
         "with L_conj(x) = L_x^T", "NAMED", "Baez 2002; Harvey 1990"),
        ("8D.P2", "8D", "gamma(e_i) gamma(e_j) + gamma(e_j) gamma(e_i) = -2 delta_ij I_16",
         "36 of 36 relations 0 <= i <= j <= 7 hold", "the polarization of 8D.P1", "NAMED", "Lawson-Michelsohn I.1"),
        ("8D.P3", "8D", "the 8 matrices gamma(e0), ..., gamma(e7), each 16 x 16, generate a real algebra",
         "dimension 256 = 2^8 = 16^2, all of M_16(R) = Cl_8 (Bott period 8); commutant dimension 1",
         "Bott table entry n = 8", "NAMED", "Lawson-Michelsohn I.4"),
        ("8D.P4", "8D", "gamma(e0) gamma(e1) ... gamma(e7)",
         "epsilon7 diag(I_8, -I_8) with the same epsilon7 as 7D.P4",
         "block (1,1) = L_(e0) L_conj(e1) L_(e2) L_conj(e3) L_(e4) L_conj(e5) L_(e6) L_conj(e7) = L_(e1) ... L_(e7) "
         "(the four signs cancel), block (2,2) = -L_(e1) ... L_(e7)", "DERIVED", "8D.C2"),
        ("8D.P5", "8D", "Psi and the split of Lambda2",
         "Psi(e0j) = e_j and Psi(e_ij) = e_i e_j for 1 <= i < j, so Psi(e0 wedge alpha + beta) = alpha + (beta contracted "
         "with phi) for alpha in R^7, beta in Lambda2 R^7; Psi maps Lambda2 onto Im O (rank 7); Psi Psi^T = 4 I_7 = "
         "I_7 + 3 I_7; ker Psi = (e0 wedge alpha + beta with alpha = -(beta contracted with phi)), dimension 21 = "
         "spin(7)_+; its orthogonal complement Lambda2_7 = (e0 wedge u + iota_u phi), dimension 7; and Psi'(e0 wedge alpha "
         "+ beta) = -alpha + (beta contracted with phi), ker Psi' has alpha = +(beta contracted with phi)",
         "c(e_i wedge e_j) on S+ is -L_(e_i) L_conj(e_j), on S- is -L_conj(e_i) L_(e_j)", "DERIVED", "8D.C2, 8D.C4"),
        ("8D.P6", "8D", "the Cayley form that goes with eta",
         "for both s: (T_s - 3 I)(T_s + I) = 0, trace T_s = 3 x 7 - 21 = 0, eigenspaces of dimension 7 (eigenvalue 3) and "
         "21 (eigenvalue -1); T_s(-3 e0 wedge u + iota_u phi) = (2 - 3 s) iota_u phi + 3 s e0 wedge u, so ker Psi is the "
         "(-1)-eigenspace of T_+1 and is not T_-1-invariant: Psi (T_+1 + I) = 4 Psi holds as a 7 x 28 integer identity "
         "and Psi (T_-1 + I) = 4 Psi fails; ker Psi' is the (-1)-eigenspace of T_-1; eta in S+ goes with "
         "Phi_+1 = e0 wedge phi + psi and eta' in S- with Phi_-1 = e0 wedge phi - psi",
         "NAMED for s = +1 (Joyce 2000 ch. 10); s = -1 is conjugate to s = +1 by e0 -> -e0; the eigenvector computation "
         "uses phi wedge iota_u phi = 2 u wedge psi, psi wedge iota_u phi = 3 *7 u (proofs g2-s6 M3 measured the 3), "
         "*8 of a 6-form gamma on R^7 = e0 wedge *7 gamma", "DERIVED", "8D.C3, 7D.C3"),
        ("8D.P7", "8D", "the two Spin(7) subalgebras",
         "ker Psi intersect ker Psi' = (beta in Lambda2 R^7 with beta contracted with phi = 0) = Lambda2_14 = g2, "
         "dimension 14; ker Psi + ker Psi' = Lambda2, 21 + 21 - 14 = 28",
         "8D.P5", "DERIVED", "8D.P5"),
        ("8D.P8", "8D", "Spin(7)-instanton sequence Lambda0 -> Lambda1 -> Lambda2_7 (of Phi_+1)",
         "sigma0(xi) f = f xi, sigma1(xi) a = pi_7(xi wedge a), read in coordinates by Psi, which kills Lambda2_21: "
         "P(xi) a = Psi(xi wedge a) = -Im(xi conj(a)); dimensions 1, 8, 7; Euler count 1 - 8 + 7 = 0; P(xi) xi = 0 and "
         "P(xi) P(xi)^T = nsq(xi) I_7 as polynomial identities, so ranks 1 and 7 for every xi != 0 and "
         "ker P(xi) = span(xi)",
         "the symbol complex of d, pi_7 d (elliptic)", "NAMED", "Donaldson-Thomas 1998"),
        ("8D.P9", "8D", "the rolled-up symbol is Clifford multiplication S- -> S+",
         "(-iota_xi, pi_7(xi wedge .)) : Lambda1 -> Lambda0 + Lambda2_7; with N(xi) the 8 x 8 integer matrix of "
         "a -> -<xi, a> e0 + Psi(xi wedge a): N(xi) a = -xi conj(a) = L_xi(c(a) eta) with c(a) eta = -L_conj(a) 1 = "
         "-conj(a); so under Lambda1 -> S-, a -> -conj(a), and Lambda0 + Lambda2_7 -> S+, (f, omega) -> f + Psi(omega), "
         "the rolled-up symbol is L_xi : S- -> S+; N(xi)^T N(xi) = N(xi) N(xi)^T = nsq(xi) I_8; isometric normalization 2: "
         "nsq(Psi(omega)) = 4 nsq(omega) on Lambda2_7, so (-iota_xi, 2 pi_7(xi wedge .)) is the norm of xi times an "
         "isometry; integer form nsq((T_+1 + I)(xi wedge a)) = 4 nsq(xi wedge a)",
         "gamma(xi) gamma(a) = -<xi, a> + c(xi wedge a), and Psi kills Lambda2_21", "DERIVED", "8D.P5, 8D.P6"),
        ("8D.P10", "8D", "equivariance under Spin(7)",
         "for omega in ker Psi (21 generators): c(omega) eta = 0 and 2 (N(R_omega xi) + N(xi) R_omega) = c(omega)_S+ N(xi), "
         "with c(omega)_S+ = -sum_(i<j) omega_ij L_(e_i) L_conj(e_j); [c(omega), gamma(xi)] = 2 gamma(R_omega xi) for every "
         "omega in Lambda2; the rotations R_omega, omega in ker Psi, are closed under commutators, a Lie algebra of "
         "dimension 21 = spin(7)",
         "C0.5 with c(omega) eta = 0", "DERIVED", "C0.5, 8D.P9"),
        ("8D.P11", "8D", "representations of Spin(7)",
         "for Spin(7)_+ (the stabilizer of eta): S+ = Lambda0 + Lambda2_7 = 1 + 7, S- = Lambda1 = 8 (the spin "
         "representation of Spin(7)); Lambda2 = Lambda2_7 + Lambda2_21 = 7 + 21 with Lambda2_21 = spin(7)",
         "triality restricted to Spin(7)", "NAMED", "Harvey 1990; Joyce 2000 ch. 10"),
        # ---------------- obstruction
        ("OB.P1", "obstruction", "constant non-abelian A: A1 = i sigma1 dx1, A2 = i sigma2 dx2, A3 = A4 = 0",
         "A1 = [[0:0, 0:1],[0:1, 0:0]], A2 = [[0:0, 1:0],[-1:0, 0:0]]; F_12 = [A1, A2] = -2 i sigma3 = "
         "[[0:-2, 0:0],[0:0, 0:2]], so F_A = -2 i sigma3 dx1 wedge dx2 != 0",
         "[i sigma1, i sigma2] = -[sigma1, sigma2] = -2 i sigma3", "DERIVED", "OB.C1, OB.C2"),
        ("OB.P2", "obstruction", "d_A d_A on Omega0(su(2)) is [F_A, .], zeroth order",
         "for constant f = i sigma1: d_A d_A f = [F_12, i sigma1] dx12 = 4 i sigma2 dx12 = [[0:0, 4:0],[-4:0, 0:0]] dx12 "
         "!= 0, so Lambda0 -> Lambda1 -> Lambda2 with d_A is not a complex as operators; in the ASD sequence "
         "d_A^+ d_A f = [F_A^+, f] with F_A^+ = -i sigma3 (dx12 + dx34) = -i sigma3 f1, and d_A^+ d_A (i sigma1) = "
         "[-i sigma3, i sigma1] f1 = 2 i sigma2 f1 = [[0:0, 2:0],[-2:0, 0:0]] f1 != 0; also F_A^- = -i sigma3 (e12 - e34) "
         "!= 0",
         "[A1, [A2, f]] - [A2, [A1, f]] = [[A1, A2], f] (Jacobi); [A1, f] = 0 and [A2, f] = 2 i sigma3 for f = i sigma1",
         "DERIVED", "OB.C2"),
        ("OB.P3", "obstruction", "the symbol does not see F_A",
         "the principal symbol of d_A is sigma(d) (x) id_su(2), unchanged by the zeroth-order term [A, .]; the su(2)-valued "
         "ASD symbol sequence has dimensions 3, 12, 9, Euler count 0, ranks 3, 9, and stays exact for every xi != 0 while "
         "d_A^+ d_A != 0; the operator sequence is a complex exactly when F_A^+ = 0 (an ASD connection), and here "
         "F_A^+ != 0",
         "the curvature is zeroth order", "NAMED", "Atiyah-Hitchin-Singer 1978"),
        ("OB.P4", "obstruction", "control with commuting constant A: A1 = i sigma3 dx1, A2 = 2 i sigma3 dx2",
         "F_A = 0 and d_A d_A = 0", "[i sigma3, 2 i sigma3] = 0", "DERIVED", "OB.C2"),
        # ---------------- tensor claim
        ("TS.P1", "tensor", "S(X^4 x Y^4) = S(X) (x) S(Y)",
         "complex spinors of R^4 + R^4: C^4 (x) C^4 = C^16, as graded tensor product with Cl(V + W) = Cl(V) graded-(x) "
         "Cl(W); chirality S+(8) = S+ (x) S+ + S- (x) S- (4 + 4 = 8), S-(8) = S+ (x) S- + S- (x) S+ (4 + 4 = 8)",
         "2^4 = 2^2 x 2^2", "NAMED", "Lawson-Michelsohn I.1 and I.5"),
        ("TS.P2", "tensor", "with the spin-invariant form C^4 = (C^4)*, C^4 (x) C^4 = End(C^4) = M_4(C)",
         "complex dimension 16 = 4 x 4; the invariant form on C^4 is skew, a symplectic form on each of S+ and S-; "
         "(g, h) in Spin(4) x Spin(4) acts on End(C^4) by A -> g A h^-1, so (-1, -1) acts trivially and the action "
         "factors through the cover of SO(4) x SO(4) in Spin(8); the identification intertwines Spin(4) x Spin(4)",
         "B(h s, h t) = B(s, t) turns h into the inverse transpose", "NAMED", "Harvey 1990; Lawson-Michelsohn I.5"),
        ("TS.P3", "tensor", "End(C^4) under the diagonal Spin(4)",
         "End(C^4) = Cl_4 (x) C = Lambda C^4 = 1 + 4 + 6 + 4 + 1 = 16; End(S+) = S+ (x) S+ = Lambda0 + Lambda2+ (1 + 3), "
         "End(S-) = Lambda4 + Lambda2- (1 + 3), Hom(S-, S+) and Hom(S+, S-) = Lambda1 and Lambda3 (4 + 4); the "
         "quaternionic structures J (x) J square to +1, giving the real form R^16 = O + O, the real spinors of 8D.C2",
         "Fierz identity", "NAMED", "Lawson-Michelsohn I.5; Harvey 1990"),
    )
    rows = []
    for pid, case, what, predicted, reason, claim, source in P:
        rows.append(row("PRED", [("id", pid), ("case", case), ("what", what), ("predicted", predicted),
                                 ("reason", reason), ("source", source), ("status", "PREDICTION"), ("claim", claim)]))
    return rows


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    rows = []
    rows.append(row("HEAD", [
        ("schema", SCHEMA), ("date", DATE), ("author", AUTHOR), ("kind", "PREDICTION_before_any_measurement"),
        ("measurement_status", "NOT_RUN"), ("matrices_built", 0),
        ("matrices_built_reads", "Clifford matrices, quaternion or octonion multiplication tables and symbol matrices "
                                 "built by this script"),
        ("script", SCRIPT_NAME), ("script_sha256", script_sha), ("run", "python3 -I -B " + SCRIPT_NAME),
        ("float_used", 0),
        ("sources", "NAMED: Lawson and Michelsohn 1989 Spin Geometry; Harvey 1990 Spinors and Calibrations; Atiyah, "
                    "Hitchin and Singer 1978 Proc. R. Soc. Lond. A 362; Donaldson and Thomas 1998 Gauge theory in higher "
                    "dimensions; Reyes Carrion 1998 Diff. Geom. Appl. 8; Baez 2002 The Octonions, Bull. AMS 39; Joyce 2000 "
                    "Compact Manifolds with Special Holonomy; Bryant 2005 Some remarks on G2-structures"),
        ("cross_reference", "proofs g2-s6 of this repository: der(O) dimension 14, the T7 split 7 + 14, the factor 3 in "
                            "psi wedge iota_v phi; measured there, cited here, not rerun"),
    ]))
    rows.append(row("FORMAT", [
        ("nsq", "nsq(x) is the sum of the squares of the coordinates of x, its squared norm; the bar notation is avoided "
                "because the bar separates fields"),
        ("matrix", "[[a, b],[c, d]] is a block or entry matrix, rows first; I_m the identity of size m"),
        ("gaussian", "a Gaussian integer re + im i is written re:im"),
        ("tensor", "(x) is the tensor product"),
        ("ids", "C0.n common conventions, 4D, 7D, 8D, OB (obstruction), TS (tensor claim); P rows are predictions, "
                "C rows conventions"),
        ("claim_tags", "MEASURED computed by this script (self-checks only); NAMED stated in the literature, not run "
                       "here; DERIVED derived by hand from NAMED identities, not run here; CONJECTURE none"),
        ("status", "PREDICTION marks a statement the measurement will test; CONVENTION marks a choice the measurement "
                   "uses or states that it changed"),
    ]))
    rows.extend(conv_rows())
    rows.extend(pred_rows())

    D = DIMS
    for tag, s in CASES:
        rows.append(row("DIMS", [
            ("case", tag), ("n", D["n_" + s]), ("spinor_module", D["rep_" + s]), ("clifford_dim", D["cl_" + s]),
            ("generated_algebra", D["alg_" + s]), ("commutant", D["comm_" + s]), ("generators", D["gens_" + s]),
            ("generator_size", D["rep_" + s]), ("lambda2_small", D["small_" + s]), ("lambda2_big", D["big_" + s]),
            ("T_small", D["T_small_" + s]), ("T_big", D["T_big"]), ("k", D["k_" + s]),
            ("stabilizer_of_eta", D["stab_" + s]), ("sequence_dims", join(D["seq_" + s])),
            ("symbol_ranks", join(D["ranks_" + s])), ("square_coefficients", D["coef_" + s]),
            ("status", "PREDICTION"), ("claim", "NAMED"),
        ]))
    for n in range(9):
        a, m, d = CL_NEG[n]
        b, m2, d2 = CL_POS[n]
        rows.append(row("BOTT", [
            ("n", n), ("cl_square_minus_1", CL_NEG_NAMES[n]), ("dim_minus", a * m * m * d),
            ("irreducible_module_minus", m * d), ("cl_square_plus_1", CL_POS_NAMES[n]), ("dim_plus", b * m2 * m2 * d2),
            ("claim", "NAMED"),
        ]))

    # comparisons with the task text
    cmp_bad = 0
    for item, expect, key in TASK_EXPECT:
        got = D[key]
        agree = 1 if got == expect else 0
        if not agree:
            cmp_bad = cmp_bad + 1
        shown_e = join(expect) if isinstance(expect, tuple) else str(expect)
        shown_g = join(got) if isinstance(got, tuple) else str(got)
        rows.append(row("COMPARE", [("item", item), ("task_text", shown_e), ("this_table", shown_g), ("agree", agree),
                                    ("claim", "MEASURED")]))
    rows.append(row("CORRECTION", [
        ("item", "8D_generator_count_task_text_16"),
        ("task_text_reads", "the 16 matrices gamma(e_i) generate an algebra of real dimension 256"),
        ("restated", "the 8 matrices gamma(e0), ..., gamma(e7), each 16 x 16, generate an algebra of real dimension "
                     "256 = 16^2 = 2^8; 16 is the size of each matrix, 8 the number of generators"),
        ("kept_beside", "the task text is kept as written in COMPARE; this row restates it"),
        ("claim", "DERIVED"),
    ]))

    # self-checks
    checks = analyze(D)
    base_failed = [c[0] for c in checks if not c[1]]
    for name, key, value in PLANTS:
        D2 = dict(D)
        D2[key] = value
        pchecks = analyze(D2)
        broken = [c[0] for c in pchecks if not c[1]]
        shown = join(value, ";") if isinstance(value, tuple) else str(value)
        rows.append(row("PLANTED", [("plant", name), ("key", key), ("planted_value", shown),
                                    ("true_value", join(D[key], ";") if isinstance(D[key], tuple) else D[key]),
                                    ("suite_checks", len(pchecks)), ("checks_broken", len(broken)),
                                    ("broken", join(broken, ";") if broken else "none"), ("claim", "MEASURED")]))
        checks.append(("planted_" + name + "_breaks_a_check", len(broken) >= 1, [("checks_broken", len(broken))]))

    failed = 0
    for name, ok, detail in checks:
        if not ok:
            failed = failed + 1
        rows.append(row("SELF_CHECK", [("name", name), ("status", "PASS" if ok else "FAIL")] + list(detail)
                        + [("claim", "MEASURED")]))
    rows.append(row("END", [
        ("rows", len(rows) + 1), ("checks", len(checks)), ("checks_failed", failed),
        ("base_suite_failed", len(base_failed)), ("compare_mismatches", cmp_bad),
        ("compare_mismatches_reads", "1 expected: the task text's 16 matrices, restated in CORRECTION"),
        ("exit", 1 if failed else 0), ("measurement_status", "NOT_RUN"), ("matrices_built", 0),
    ]))

    data = ("\n".join(rows) + "\n").encode("utf-8")
    with open(os.path.join(here, HBP_NAME), "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(os.path.join(here, SIDECAR_NAME), "wb") as fh:
        fh.write((digest + "  " + HBP_NAME + "\n").encode("utf-8"))
    print(row("DONE", [("hbp_sha256", digest), ("script_sha256", script_sha), ("checks", len(checks)),
                       ("checks_failed", failed), ("compare_mismatches", cmp_bad)]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
