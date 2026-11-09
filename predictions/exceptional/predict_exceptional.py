"""predict_exceptional.py  --  PREDICTION, written before any measurement.

Gegenbauer (zonal) sums S_k, k = 1..30, for three exceptional point sets, computed only from
their published distance distributions (cosine of the angle from one point -> count):

  600-cell vertices = the 120 unit icosians = 2I on S^3 in R^4
  E8 roots, 240 points on S^7 in R^8
  Leech lattice minimal vectors, 196560 points on S^23 in R^24

NAMED inputs: Conway and Sloane, Sphere Packings, Lattices and Groups, ch. 4 and 14 (distance
distributions); Delsarte, Goethals and Seidel 1977, Spherical codes and designs (a distance-invariant
set X is a spherical t-design iff S_k = 0 for k = 1..t).

Hard rule kept here: no lattice vector, root, quaternion or octonion is built, and no sum runs over
a point set.  Every sum below runs over the rows of a distance distribution (9, 5 and 7 rows).

Arithmetic is exact: integers; rationals as reduced (num, den) integer pairs with den > 0; elements
of Q(sqrt5) as pairs (a, b) of such rationals meaning a + b*sqrt5, with
(a + b r)(c + d r) = (ac + 5bd) + (ad + bc) r,  r*r = 5.
Division appears only as divmod with a remainder that is checked or provably zero.

Run:  python3 -I -B predict_exceptional.py
It writes PREDICTION-EXCEPTIONAL.hbp and PREDICTION-EXCEPTIONAL.hbp.sha256 beside this file and
exits 0 only when every self-check passes (nonzero otherwise).
"""

import hashlib
import os
import sys
from math import gcd

KMAX = 30
SCHEMA = "predict-exceptional-v1"
DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
HBP_NAME = "PREDICTION-EXCEPTIONAL.hbp"
SIDECAR_NAME = "PREDICTION-EXCEPTIONAL.hbp.sha256"

# ---------------------------------------------------------------- rationals


def rat(n, d=1):
    if d == 0:
        raise ValueError("rational with zero denominator")
    if d < 0:
        n, d = -n, -d
    g = gcd(n, d)
    qn, rn = divmod(n, g)
    qd, rd = divmod(d, g)
    if rn or rd:
        raise ValueError("gcd did not divide")
    return (qn, qd)


R0 = (0, 1)
R1 = (1, 1)


def r_add(x, y):
    return rat(x[0] * y[1] + y[0] * x[1], x[1] * y[1])


def r_neg(x):
    return (-x[0], x[1])


def r_mul(x, y):
    return rat(x[0] * y[0], x[1] * y[1])


def r_mul_int(x, m):
    return rat(x[0] * m, x[1])


def r_div_int(x, m):
    return rat(x[0], x[1] * m)


def r_sign(x):
    if x[0] > 0:
        return 1
    if x[0] < 0:
        return -1
    return 0


def r_cmp(x, y):
    # sign of x - y, dens positive
    return r_sign(rat(x[0] * y[1] - y[0] * x[1], x[1] * y[1]))


def r_str(x):
    return str(x[0]) + ":" + str(x[1])


# ---------------------------------------------------------------- Q(sqrt5): (a, b) = a + b*sqrt5

Q0 = (R0, R0)
Q1 = (R1, R0)


def q_add(x, y):
    return (r_add(x[0], y[0]), r_add(x[1], y[1]))


def q_neg(x):
    return (r_neg(x[0]), r_neg(x[1]))


def q_sub(x, y):
    return q_add(x, q_neg(y))


def q_mul(x, y):
    a, b = x
    c, d = y
    return (r_add(r_mul(a, c), r_mul_int(r_mul(b, d), 5)), r_add(r_mul(a, d), r_mul(b, c)))


def q_mul_int(x, m):
    return (r_mul_int(x[0], m), r_mul_int(x[1], m))


def q_div_int(x, m):
    return (r_div_int(x[0], m), r_div_int(x[1], m))


def q_is_zero(x):
    return x[0][0] == 0 and x[1][0] == 0


def q_conj(x):
    return (x[0], r_neg(x[1]))


def q_sign(x):
    # exact sign of a + b*sqrt5
    sa = r_sign(x[0])
    sb = r_sign(x[1])
    if sb == 0:
        return sa
    if sa == 0:
        return sb
    if sa == sb:
        return sa
    # opposite signs: compare a*a with 5*b*b
    c = r_cmp(r_mul(x[0], x[0]), r_mul_int(r_mul(x[1], x[1]), 5))
    if c == 0:
        return 0
    return sa if c > 0 else sb


def q_str(x):
    return r_str(x[0]) + "+" + r_str(x[1]) + "r5"


def qr(n, d):
    return (rat(n, d), R0)


def qs(an, ad, bn, bd):
    return (rat(an, ad), rat(bn, bd))


# ---------------------------------------------------------------- the three NAMED distributions
# cosine -> count, from one point; each set is distance-invariant (transitive automorphism group, NAMED).
# 600-cell cosines: phi over 2 = (1 + sqrt5) over 4, (phi - 1) over 2 = (sqrt5 - 1) over 4.

SETS = [
    {
        "name": "600cell_2I_icosians",
        "n": 4,
        "N": 120,
        "expected_t": 11,
        "dist": [
            (qr(1, 1), 1),
            (qs(1, 4, 1, 4), 12),
            (qr(1, 2), 20),
            (qs(-1, 4, 1, 4), 12),
            (qr(0, 1), 30),
            (qs(1, 4, -1, 4), 12),
            (qr(-1, 2), 20),
            (qs(-1, 4, -1, 4), 12),
            (qr(-1, 1), 1),
        ],
    },
    {
        "name": "E8_roots",
        "n": 8,
        "N": 240,
        "expected_t": 7,
        "half_rank": 4,
        "min_norm_half": 1,
        "dist": [
            (qr(1, 1), 1),
            (qr(1, 2), 56),
            (qr(0, 1), 126),
            (qr(-1, 2), 56),
            (qr(-1, 1), 1),
        ],
    },
    {
        "name": "Leech_min_vectors",
        "n": 24,
        "N": 196560,
        "expected_t": 11,
        "half_rank": 12,
        "min_norm_half": 2,
        "dist": [
            (qr(1, 1), 1),
            (qr(1, 2), 4600),
            (qr(1, 4), 47104),
            (qr(0, 1), 93150),
            (qr(-1, 4), 47104),
            (qr(-1, 2), 4600),
            (qr(-1, 1), 1),
        ],
    },
]

# Coxeter degrees (NAMED, Humphreys, Reflection Groups and Coxeter Groups, table 3.1)
COXETER_DEGREES = {"600cell_2I_icosians": [2, 12, 20, 30], "E8_roots": [2, 8, 12, 14, 18, 20, 24, 30]}

# frozen SU(2)-ladder prediction for 2I: first invariant degree 12 (Klein's invariant f12), NAMED input
SU2_LADDER_FIRST_INVARIANT = 12

# ---------------------------------------------------------------- zonal polynomials


def zonal_norm(n, x, kmax):
    # normalised zonal polynomial of S^(n-1), Q_k(1) = 1:
    # Q_0 = 1, Q_1 = x, (k + n - 3) Q_k = (2k + n - 4) x Q_(k-1) - (k - 1) Q_(k-2)
    vals = [Q1, x]
    for k in range(2, kmax + 1):
        num = q_sub(q_mul_int(q_mul(x, vals[k - 1]), 2 * k + n - 4), q_mul_int(vals[k - 2], k - 1))
        vals.append(q_div_int(num, k + n - 3))
    return vals


def zonal_raw(n, x, kmax):
    # unnormalised Gegenbauer C_k^(lambda), 2 lambda = a = n - 2:
    # C_0 = 1, C_1 = a x, k C_k = (2k + a - 2) x C_(k-1) - (k + a - 2) C_(k-2)
    a = n - 2
    vals = [Q1, q_mul_int(x, a)]
    for k in range(2, kmax + 1):
        num = q_sub(q_mul_int(q_mul(x, vals[k - 1]), 2 * k + a - 2), q_mul_int(vals[k - 2], k + a - 2))
        vals.append(q_div_int(num, k))
    return vals


def binom(m, k):
    r = 1
    for i in range(k):
        r, rem = divmod(r * (m - i), i + 1)
        if rem:
            raise ValueError("binomial step not exact")
    return r


def gsums(n, dist, kmax):
    s = [Q0 for _ in range(kmax + 1)]
    for c, m in dist:
        v = zonal_norm(n, c, kmax)
        for k in range(kmax + 1):
            s[k] = q_add(s[k], q_mul_int(v[k], m))
    return s


def series_count(parts, kmax):
    # coefficients of the product over d in parts of 1 over (1 - t^d)
    c = [1] + [0 for _ in range(kmax)]
    for d in parts:
        for k in range(d, kmax + 1):
            c[k] += c[k - d]
    return c


def dim_modular(w):
    # dimension of M_w for SL2(Z)
    if w < 0:
        return 0
    q, r = divmod(w, 2)
    if r:
        return 0
    q, r = divmod(w, 12)
    if r == 2:
        return q
    return q + 1


def count_of(dist, c):
    for x, m in dist:
        if x == c:
            return m
    return None


def is_antipodal(dist):
    return all(count_of(dist, q_neg(c)) == m for c, m in dist)


def is_galois_stable(dist):
    return all(count_of(dist, q_conj(c)) == m for c, m in dist)


def counts_sum(dist):
    return sum(m for _, m in dist)


def dist_str(dist):
    return ",".join(q_str(c) + "#" + str(m) for c, m in dist)


def first_break(s):
    for k in range(1, len(s)):
        if not q_is_zero(s[k]):
            return k
    return None


def row(tag, pairs):
    return tag + "|" + "|".join(k + "=" + str(v) for k, v in pairs) + "|json=0"


# ---------------------------------------------------------------- main


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    rows = []
    checks = []

    def check(name, setname, ok, detail):
        checks.append(row("CHECK", [("name", name), ("set", setname), ("pass", 1 if ok else 0)] + detail))
        return ok

    rows.append(row("HEAD", [
        ("schema", SCHEMA),
        ("date", DATE),
        ("author", AUTHOR),
        ("kind", "PREDICTION_before_any_measurement"),
        ("measurement_status", "NOT_RUN"),
        ("script", "predict_exceptional.py"),
        ("script_sha256", script_sha),
        ("run", "python3 -I -B predict_exceptional.py"),
        ("sources", "NAMED: Conway-Sloane SPLAG ch.4 and ch.14 (distance distributions); "
                    "Delsarte-Goethals-Seidel 1977 Spherical codes and designs (t-design iff S_k zero for k 1..t); "
                    "Klein 1884 Ikosaeder (invariant f12 of 2I, Molien series (1+t^30) over ((1-t^12)(1-t^20))); "
                    "Humphreys 1990 table 3.1 (Coxeter degrees H4, E8); "
                    "Venkov 1984 and 2001 (harmonic theta series of an even unimodular lattice)"),
        ("built_points", 0),
        ("sums_over_point_sets", 0),
        ("float_used", 0),
        ("kmax", KMAX),
    ]))
    rows.append(row("FORMAT", [
        ("rational", "num:den_reduced_den_positive"),
        ("qsqrt5", "a+br5_means_a_plus_b_times_sqrt5_with_a_b_rationals"),
        ("dist_entry", "cosine#count"),
        ("S_k", "sum_over_distribution_rows_of_count_times_Q_k(cosine)"),
        ("pair_sum", "N_times_S_k_equals_sum_over_ordered_pairs_x_y_of_Q_k(x.y)"),
    ]))
    rows.append(row("RECURRENCE", [
        ("start", "Gegenbauer C_k^lambda with lambda equal to (n-2) over 2, so 2 lambda equals n-2"),
        ("raw", "k C_k(x) equals 2(k+lambda-1) x C_(k-1)(x) minus (k+2lambda-2) C_(k-2)(x); C_0 is 1, C_1 is 2lambda x"),
        ("norm", "C_k(1) equals binom(k+2lambda-1,k) equals binom(k+n-3,k); C_(k-1)(1) over C_k(1) is k over (k+n-3);"
                 " C_(k-2)(1) over C_k(1) is k(k-1) over ((k+n-3)(k+n-4))"),
        ("derived", "Q_k equals C_k over C_k(1): (k+n-3) Q_k(x) equals (2k+n-4) x Q_(k-1)(x) minus (k-1) Q_(k-2)(x);"
                    " Q_0 is 1, Q_1 is x"),
        ("n4", "(k+1) Q_k equals 2k x Q_(k-1) minus (k-1) Q_(k-2); Q_k is U_k over (k+1), Chebyshev second kind"),
        ("n8", "(k+5) Q_k equals (2k+4) x Q_(k-1) minus (k-1) Q_(k-2)"),
        ("n24", "(k+21) Q_k equals (2k+20) x Q_(k-1) minus (k-1) Q_(k-2)"),
        ("at_1", "if Q_(k-1)(1) and Q_(k-2)(1) are 1 then Q_k(1) is (2k+n-4-(k-1)) over (k+n-3), which is 1"),
        ("verified_here", "both routes evaluated at every cosine and at 1; see CHECK route_agreement and Q_k_at_1_eq_1"),
    ]))

    results = {}
    failed = 0

    for st in SETS:
        name = st["name"]
        n = st["n"]
        nn = st["N"]
        dist = st["dist"]
        csum = counts_sum(dist)
        anti = is_antipodal(dist)
        gal = is_galois_stable(dist)
        rows.append(row("SET", [
            ("set", name), ("n", n), ("sphere", "S^" + str(n - 1)), ("N", nn), ("rows", len(dist)),
            ("counts_sum", csum), ("antipodal", 1 if anti else 0), ("galois_stable", 1 if gal else 0),
            ("dist", dist_str(dist)), ("basis", "NAMED_distance_distribution"),
        ]))
        s = gsums(n, dist, KMAX)
        results[name] = s
        for k in range(1, KMAX + 1):
            rows.append(row("SK", [
                ("set", name), ("n", n), ("k", k), ("S_a", r_str(s[k][0])), ("S_b", r_str(s[k][1])),
                ("zero", 1 if q_is_zero(s[k]) else 0), ("sign", q_sign(s[k])), ("status", "PREDICTION"),
            ]))
        fb = first_break(s)
        t_pred = (fb - 1) if fb is not None else None
        rows.append(row("TPRED", [
            ("set", name), ("n", n), ("N", nn),
            ("t_pred", t_pred if t_pred is not None else "ge" + str(KMAX)),
            ("rule", "first_k_ge_1_with_S_k_nonzero_minus_1"),
            ("expected_t", st["expected_t"]),
            ("agree_expected", 1 if t_pred == st["expected_t"] else 0),
            ("status", "PREDICTION"),
        ]))
        if fb is not None:
            rows.append(row("FIRST_BREAK", [
                ("set", name), ("k", fb), ("S_a", r_str(s[fb][0])), ("S_b", r_str(s[fb][1])),
                ("pair_sum_a", r_str(r_mul_int(s[fb][0], nn))), ("pair_sum_b", r_str(r_mul_int(s[fb][1], nn))),
                ("sign", q_sign(s[fb])), ("status", "PREDICTION"),
            ]))
        nz = [k for k in range(1, KMAX + 1) if not q_is_zero(s[k])]
        rows.append(row("PROFILE", [
            ("set", name), ("kmax", KMAX), ("nonzero_count", len(nz)),
            ("nonzero_k", ",".join(str(k) for k in nz)),
            ("zero_k", ",".join(str(k) for k in range(1, KMAX + 1) if q_is_zero(s[k]))),
            ("values", ",".join(str(k) + "@" + q_str(s[k]) for k in nz)),
            ("status", "PREDICTION_past_the_limit_profile"),
        ]))

        # ---- self-checks that can fail
        ok = check("counts_sum_eq_N", name, csum == nn, [("counts_sum", csum), ("N", nn)])
        failed += 0 if ok else 1
        ok = check("S_0_eq_N", name, s[0] == qr(nn, 1), [("S_0", q_str(s[0]))])
        failed += 0 if ok else 1
        ones = zonal_norm(n, Q1, KMAX)
        ok = check("Q_k_at_1_eq_1", name, all(v == Q1 for v in ones), [("k_range", "0.." + str(KMAX)), ("n", n)])
        failed += 0 if ok else 1
        agree = True
        for c, _ in [(Q1, 1)] + dist:
            vn = zonal_norm(n, c, KMAX)
            vr = zonal_raw(n, c, KMAX)
            for k in range(KMAX + 1):
                b = binom(k + n - 3, k)
                if vr[k] != q_mul_int(vn[k], b):
                    agree = False
        raw_one = zonal_raw(n, Q1, KMAX)
        norm_ok = all(raw_one[k] == qr(binom(k + n - 3, k), 1) for k in range(KMAX + 1))
        ok = check("route_agreement", name, agree and norm_ok,
                   [("normalised_times_binom_eq_raw_at_every_cosine", 1 if agree else 0),
                    ("raw_at_1_eq_binom", 1 if norm_ok else 0)])
        failed += 0 if ok else 1
        odd_zero = all(q_is_zero(s[k]) for k in range(1, KMAX + 1, 2))
        ok = check("antipodal_odd_k_vanish", name, anti and odd_zero,
                   [("antipodal", 1 if anti else 0), ("odd_k_all_zero", 1 if odd_zero else 0)])
        failed += 0 if ok else 1
        pos = all(q_sign(s[k]) >= 0 for k in range(KMAX + 1))
        ok = check("positivity_S_k_ge_0", name, pos,
                   [("basis", "NAMED_Delsarte-Schoenberg_N_S_k_is_a_sum_of_squares_over_dim_Harm_k")])
        failed += 0 if ok else 1
        if name == "600cell_2I_icosians":
            rational = all(q_is_zero((R0, s[k][1])) for k in range(KMAX + 1))
            ok = check("galois_stable_S_k_rational", name, gal and rational,
                       [("galois_stable", 1 if gal else 0), ("every_S_k_b_zero", 1 if rational else 0),
                        ("basis", "sqrt5_to_minus_sqrt5_permutes_the_distribution_rows")])
            failed += 0 if ok else 1
        if name in COXETER_DEGREES:
            degs = COXETER_DEGREES[name]
            h = series_count([d for d in degs if d != 2], KMAX)
            forced = [k for k in range(1, KMAX + 1) if h[k] == 0]
            ok_imp = all(q_is_zero(s[k]) for k in forced)
            ok = check("coxeter_invariant_zeros", name, ok_imp,
                       [("degrees", ",".join(str(d) for d in degs)),
                        ("no_invariant_harmonic_at_k", ",".join(str(k) for k in forced)),
                        ("rule", "no_invariant_harmonic_of_degree_k_implies_S_k_zero")])
            failed += 0 if ok else 1
            rows.append(row("NAMED_ZEROS", [
                ("set", name), ("source", "Coxeter_degrees_harmonic_invariants"),
                ("forced_zero_k", ",".join(str(k) for k in forced)),
                ("invariant_harmonic_k", ",".join(str(k) for k in range(1, KMAX + 1) if h[k] > 0)),
                ("all_forced_zero_hold", 1 if ok_imp else 0),
            ]))
        if "half_rank" in st:
            forced = [k for k in range(1, KMAX + 1)
                      if dim_modular(st["half_rank"] + k - 12 * st["min_norm_half"]) == 0]
            ok_imp = all(q_is_zero(s[k]) for k in forced)
            ok = check("harmonic_theta_zeros", name, ok_imp,
                       [("rule", "theta_P_is_Delta^" + str(st["min_norm_half"]) + "_times_M_w_with_w_k+"
                         + str(st["half_rank"]) + "-" + str(12 * st["min_norm_half"])
                         + "_so_dim_M_w_zero_implies_S_k_zero"),
                        ("dim_M_w_zero_at_k", ",".join(str(k) for k in forced))])
            failed += 0 if ok else 1
            rows.append(row("NAMED_ZEROS", [
                ("set", name), ("source", "Venkov_harmonic_theta_series_even_unimodular"),
                ("forced_zero_k", ",".join(str(k) for k in forced)),
                ("all_forced_zero_hold", 1 if ok_imp else 0),
            ]))

    # ---- cross-check: 600-cell first break vs frozen SU(2)-ladder prediction for 2I
    s6 = results["600cell_2I_icosians"]
    fb6 = first_break(s6)
    base = series_count([12, 20], KMAX)
    molien = [base[k] + (base[k - 30] if k >= 30 else 0) for k in range(KMAX + 1)]
    molien_first = next((k for k in range(1, KMAX + 1) if molien[k] > 0), None)
    exact = sum(1 for k in range(1, KMAX + 1) if r_mul_int(s6[k][0], k + 1) == rat(120 * molien[k], 1)
                and s6[k][1] == R0)
    agree_cross = (fb6 == SU2_LADDER_FIRST_INVARIANT) and (molien_first == SU2_LADDER_FIRST_INVARIANT)
    rows.append(row("CROSS", [
        ("set", "600cell_2I_icosians"),
        ("pred_first_break_k", fb6),
        ("su2_ladder_first_invariant_degree", SU2_LADDER_FIRST_INVARIANT),
        ("su2_source", "frozen_SU(2)-ladder_prediction_for_2I:Klein_invariant_f12_NAMED"),
        ("molien_first_invariant_degree", molien_first),
        ("reason", "SU(2) is S^3; the spin-j matrix coefficients of SU(2) are exactly the harmonics of degree 2j on S^3"
                   " (Peter-Weyl: Harm_k(S^3) is V_j tensor V_j with j equal to k over 2, dimension (k+1)^2);"
                   " for g with real part cos a the zonal Q_k(cos a) is U_k(cos a) over (k+1), that is chi_j(g) over (k+1);"
                   " so S_k equals 120 m_k over (k+1) with m_k the dimension of the 2I-invariants in V_j, i.e. of the"
                   " 2I-invariant degree-k polynomials on C^2; S_k is nonzero iff an invariant of degree k exists,"
                   " and the first one is Klein's f12"),
        ("per_degree_exact", str(exact) + ":" + str(KMAX)),
        ("per_degree_rule", "S_k times (k+1) equals 120 times m_k, m_k from Molien (1+t^30) over ((1-t^12)(1-t^20))"),
        ("agree", 1 if agree_cross else 0),
        ("status", "PREDICTION_vs_frozen_PREDICTION"),
    ]))
    ok = check("cross_su2_ladder", "600cell_2I_icosians", agree_cross and exact == KMAX,
               [("first_break", fb6), ("su2_ladder", SU2_LADDER_FIRST_INVARIANT),
                ("molien_first", molien_first), ("per_degree_exact", str(exact) + ":" + str(KMAX))])
    failed += 0 if ok else 1

    # ---- control A: E8 with 126 -> 125 must be flagged by the count check
    e8 = SETS[1]
    bad_a = [(c, (125 if m == 126 else m)) for c, m in e8["dist"]]
    flag_a = counts_sum(bad_a) != e8["N"]
    rows.append(row("CONTROL", [
        ("name", "E8_126_to_125"), ("counts_sum", counts_sum(bad_a)), ("N", e8["N"]),
        ("count_check_flags", 1 if flag_a else 0), ("must_flag", 1), ("pass", 1 if flag_a else 0),
    ]))
    ok = check("control_E8_126_to_125_flagged", "E8_roots", flag_a, [("counts_sum", counts_sum(bad_a))])
    failed += 0 if ok else 1

    # ---- control B: Leech 47104 -> 47105 on the +1 over 4 side must break odd-k vanishing and be flagged
    le = SETS[2]
    plus_quarter = qr(1, 4)
    bad_b = [(c, (47105 if c == plus_quarter else m)) for c, m in le["dist"]]
    sb = gsums(le["n"], bad_b, KMAX)
    odd_nz = [k for k in range(1, KMAX + 1, 2) if not q_is_zero(sb[k])]
    anti_b = is_antipodal(bad_b)
    flag_b = (not anti_b) and len(odd_nz) > 0
    rows.append(row("CONTROL", [
        ("name", "Leech_47104_to_47105_at_plus_1:4"), ("counts_sum", counts_sum(bad_b)), ("N", le["N"]),
        ("antipodal", 1 if anti_b else 0), ("odd_k_nonzero_count", len(odd_nz)),
        ("odd_k_nonzero", ",".join(str(k) for k in odd_nz)),
        ("S_1", q_str(sb[1])), ("count_check_flags", 1 if counts_sum(bad_b) != le["N"] else 0),
        ("antipodal_odd_check_flags", 1 if flag_b else 0), ("must_flag", 1), ("pass", 1 if flag_b else 0),
    ]))
    ok = check("control_Leech_47105_flagged", "Leech_min_vectors", flag_b,
               [("odd_k_nonzero", len(odd_nz)), ("antipodal", 1 if anti_b else 0)])
    failed += 0 if ok else 1

    rows.extend(checks)
    rows.append(row("END", [
        ("rows", len(rows) + 1), ("checks", len(checks)), ("checks_failed", failed),
        ("exit", 1 if failed else 0), ("measurement_status", "NOT_RUN"),
    ]))

    data = ("\n".join(rows) + "\n").encode("utf-8")
    hbp_path = os.path.join(here, HBP_NAME)
    with open(hbp_path, "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(os.path.join(here, SIDECAR_NAME), "wb") as fh:
        fh.write((digest + "  " + HBP_NAME + "\n").encode("utf-8"))
    print(row("DONE", [("hbp_sha256", digest), ("script_sha256", script_sha), ("checks", len(checks)),
                       ("checks_failed", failed)]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
