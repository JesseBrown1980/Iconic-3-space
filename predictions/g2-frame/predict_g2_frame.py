"""predict_g2_frame.py  --  PREDICTION, written before any measurement of the frame group.

Three blocks of rows, every value from a closed formula:

  ORBIT  the 14 points +-e_1..+-e_7 on S^6 in R^7 against the uniform measure on S^6: coordinate
         moments E[X_1^a X_2^b] for a + b up to 8, the design degree read over every monomial in all
         seven coordinates up to degree 8, and the first differences.
  G2     the exact moments E_G2[chi_7^n] and E_G2[chi_14^n], n = 0..8, of the continuous compact group
         G2 by the Weyl integration formula, cross-checked by Brauer-Klimyk tensoring with weight
         multisets and dominant-weight peeling; also the G2-fixed dimensions of Sym^k and Wedge^k of 7.
  FRAME  predictions for the group H of order 1344 inside G2 that permutes {+-e_1..+-e_7} and keeps
         the octonion product, to be tested later against H itself.

Hard rule kept here: no octonion multiplication table is built, no element of H and no element of G2
is built.  The only group built is the Weyl group W of G2, a reflection group of 2x2 integer matrices
on the weight lattice, generated from the stated root list.  The only point set built is the 14-point
orbit list, read for its own coordinate moments.

Arithmetic is exact: integers; rationals as reduced num:den integer pairs with den > 0; Laurent
polynomials in two variables as dicts from integer exponent pairs to integer coefficients.  Division
appears only as divmod with a remainder that is checked.

Weight-lattice basis (stated): the simple roots a1 (short) and a2 (long) of G2.  The weight lattice of
G2 equals its root lattice; the script checks this by solving for the fundamental weights, w1 = 2a1 + a2
and w2 = 3a1 + 2a2, and finding the change of basis has determinant 1.  Inner product on the basis
(Gram matrix, short roots of norm 2): (a1,a1) = 2, (a1,a2) = -3, (a2,a2) = 6.

Weyl integration (NAMED): for a class function f, E_G[f] is one over |W| times CT( f Delta conj(Delta) ),
Delta = e^rho times the product over positive roots a of (1 - e^-a), conj(e^m) = e^-m on the maximal
torus, so Delta conj(Delta) is the product over positive a of (1 - e^-a)(1 - e^a).  CT = constant term.

Sphere moments (NAMED): for Y uniform on S^(d-1) in R^d, E[Y^alpha] is 0 when some alpha_i is odd, and
otherwise the product of (alpha_i - 1)!! over d(d+2)...(d+|alpha|-2); here d = 7.

Run:  python3 -I -B predict_g2_frame.py
It writes PREDICTION-G2-FRAME.hbp and PREDICTION-G2-FRAME.hbp.sha256 beside this file and exits 0 only
when every self-check passes (nonzero otherwise).
"""

import hashlib
import itertools
import os
import sys
from math import gcd

SCHEMA = "predict-g2-frame-v1"
DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
SCRIPT_NAME = "predict_g2_frame.py"
HBP_NAME = "PREDICTION-G2-FRAME.hbp"
SIDECAR_NAME = "PREDICTION-G2-FRAME.hbp.sha256"

DIM = 7  # S^6 sits in R^7
ORBIT_DEG = 8  # moments up to total degree 8
NMAX = 8  # tensor powers n = 0..8
ROUTE2_MIN = 6  # the two G2 routes must agree for every n up to this
SYM_MAX = 8  # Sym^k of the 7 for k = 0..8
GROUP_CAP = 2000  # a generated reflection group above this size is reported as not finite here

# Expectations stated in the task: compared, never forced.
TASK_CHI7 = (1, 0, 1, 1, 4, 10, 35)
TASK_DESIGN_DEGREE = 3
TASK_X4_DIFF = (2, 21)
# OEIS A059710 as listed on oeis.org, opened by the prediction agent on 2026-11-09.  Title read there:
# dimension of space of invariants of n-th tensor power of 7-dimensional irreducible representation of G_2.
OEIS_A059710 = (1, 0, 1, 1, 4, 10, 35, 120, 455)

# ---------------------------------------------------------------- G2 root data (stated)
GRAM = ((2, -3), (-3, 6))
SIMPLE = ((1, 0), (0, 1))
SHORT_POS = ((1, 0), (1, 1), (2, 1))
LONG_POS = ((0, 1), (3, 1), (3, 2))
# Planted wrong roots, one coordinate changed in each; each plant must break at least one self-check.
PLANTS = (
    ("short_root_2_1_planted_as_2_2", ((1, 0), (1, 1), (2, 2)), LONG_POS),
    ("long_root_3_2_planted_as_3_3", SHORT_POS, ((0, 1), (3, 1), (3, 3))),
)

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


def r_sub(x, y):
    return rat(x[0] * y[1] - y[0] * x[1], x[1] * y[1])


def r_mul(x, y):
    return rat(x[0] * y[0], x[1] * y[1])


def r_str(x):
    if x is None:
        return "NONE"
    return str(x[0]) + ":" + str(x[1])


def exact_div(n, d):
    q, r = divmod(n, d)
    if r:
        return None
    return q


def ipow(b, e):
    if e < 0:
        raise ValueError("negative exponent")
    out = 1
    for _ in range(e):
        out = out * b
    return out


def dfact(m):
    # double factorial of an odd m >= -1, with (-1)!! = 1
    out = 1
    while m > 1:
        out = out * m
        m = m - 2
    return out


def fact(m):
    out = 1
    for i in range(2, m + 1):
        out = out * i
    return out


def join(xs, sep=","):
    return sep.join(str(x) for x in xs)


def vec_str(v):
    return "(" + str(v[0]) + "," + str(v[1]) + ")"


# ---------------------------------------------------------------- rows


def clean(v):
    s = str(v)
    for bad in ("|", "=", "\n", "\r"):
        if bad in s:
            raise ValueError("value holds a reserved character: " + s)
    return s


def row(tag, pairs):
    parts = [tag]
    for k, v in pairs:
        parts.append(k + "=" + clean(v))
    parts.append("json=0")
    return "|".join(parts)


# ---------------------------------------------------------------- ORBIT block


def orbit_points():
    pts = []
    for i in range(DIM):
        for s in (1, -1):
            v = [0] * DIM
            v[i] = s
            pts.append(tuple(v))
    return pts


def orbit_moment(pts, alpha):
    total = 0
    for p in pts:
        t = 1
        for j in range(DIM):
            t = t * ipow(p[j], alpha[j])
        total = total + t
    return rat(total, len(pts))


def sphere_moment(alpha):
    deg = sum(alpha)
    for a in alpha:
        if a & 1:
            return R0
    num = 1
    for a in alpha:
        num = num * dfact(a - 1)
    den = 1
    k = 0
    while 2 * k < deg:  # factors DIM, DIM + 2, ..., DIM + deg - 2
        den = den * (DIM + 2 * k)
        k = k + 1
    return rat(num, den)


def exponent_vectors(deg):
    out = []
    for combo in itertools.combinations_with_replacement(range(DIM), deg):
        alpha = [0] * DIM
        for i in combo:
            alpha[i] = alpha[i] + 1
        out.append(tuple(alpha))
    return out


def pattern(alpha):
    nz = sorted((a for a in alpha if a), reverse=True)
    return ".".join(str(a) for a in nz) if nz else "1"


def q_power_mean(moment, pts, k):
    # E[(sum_i Y_i^2)^k] by the multinomial expansion; equals 1 on the unit sphere
    total = R0
    for beta in exponent_vectors(k):
        coef = fact(k)
        for b in beta:
            q = exact_div(coef, fact(b))
            if q is None:
                raise ValueError("multinomial not integral")
            coef = q
        alpha = tuple(2 * b for b in beta)
        m = moment(alpha) if pts is None else moment(pts, alpha)
        total = r_add(total, r_mul((coef, 1), m))
    return total


# ---------------------------------------------------------------- lattice and Weyl group


def v_add(u, v):
    return (u[0] + v[0], u[1] + v[1])


def v_neg(u):
    return (-u[0], -u[1])


def ip(u, v):
    return (u[0] * v[0] * GRAM[0][0] + u[0] * v[1] * GRAM[0][1]
            + u[1] * v[0] * GRAM[1][0] + u[1] * v[1] * GRAM[1][1])


def coroot_pair(v, a):
    return exact_div(2 * ip(v, a), ip(a, a))


def reflect(v, a):
    num = 2 * ip(v, a)
    den = ip(a, a)
    out = []
    for j in range(2):
        q = exact_div(num * a[j], den)
        if q is None:
            return None
        out.append(v[j] - q)
    return (out[0], out[1])


def simple_pair(v, i):
    p = coroot_pair(v, SIMPLE[i])
    if p is None:
        raise ValueError("simple coroot pairing not integral")
    return p


IDENT = ((1, 0), (0, 1))


def m_apply(m, v):
    return (m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1])


def m_mul(m, n):
    return ((m[0][0] * n[0][0] + m[0][1] * n[1][0], m[0][0] * n[0][1] + m[0][1] * n[1][1]),
            (m[1][0] * n[0][0] + m[1][1] * n[1][0], m[1][0] * n[0][1] + m[1][1] * n[1][1]))


def m_det(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def reflection_matrix(a):
    c0 = reflect((1, 0), a)
    c1 = reflect((0, 1), a)
    if c0 is None or c1 is None:
        return None
    return ((c0[0], c1[0]), (c0[1], c1[1]))


def preserves_gram(m):
    for u in SIMPLE:
        for v in SIMPLE:
            if ip(m_apply(m, u), m_apply(m, v)) != ip(u, v):
                return False
    return True


def generate_group(gens):
    seen = {IDENT}
    frontier = [IDENT]
    while frontier:
        nxt = []
        for m in frontier:
            for g in gens:
                p = m_mul(g, m)
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
                    if len(seen) > GROUP_CAP:
                        return sorted(seen), False
        frontier = nxt
    return sorted(seen), True


def fundamental_weights():
    # P[i][j] is the pairing of basis vector e_j with simple coroot i; w_i solves P w = e_i
    basis = ((1, 0), (0, 1))
    P = [[simple_pair(e, i) for e in basis] for i in range(2)]
    det = P[0][0] * P[1][1] - P[0][1] * P[1][0]
    ws = []
    for i in range(2):
        r0 = 1 if i == 0 else 0
        r1 = 1 if i == 1 else 0
        x = exact_div(P[1][1] * r0 - P[0][1] * r1, det)
        y = exact_div(-P[1][0] * r0 + P[0][0] * r1, det)
        ws.append((x, y))
    return P, det, ws[0], ws[1]


# ---------------------------------------------------------------- Laurent polynomials


def lp_clean(p):
    return {k: c for k, c in p.items() if c}


def lp_mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            k = (a[0] + b[0], a[1] + b[1])
            out[k] = out.get(k, 0) + ca * cb
    return lp_clean(out)


def lp_add(p, q):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, 0) + c
    return lp_clean(out)


def lp_conj(p):
    return {(-k[0], -k[1]): c for k, c in p.items()}


def lp_from_weights(ws):
    out = {}
    for w in ws:
        out[w] = out.get(w, 0) + 1
    return lp_clean(out)


def lp_ct_pair(f, d):
    # constant term of f times d
    total = 0
    for k, c in f.items():
        total = total + c * d.get((-k[0], -k[1]), 0)
    return total


def lp_alt(group, lam):
    out = {}
    for m in group:
        k = m_apply(m, lam)
        out[k] = out.get(k, 0) + m_det(m)
    return lp_clean(out)


def weyl_dim(lam, rho, pos):
    num = 1
    den = 1
    for a in pos:
        num = num * ip(v_add(lam, rho), a)
        den = den * ip(rho, a)
    if den == 0:
        return None
    return rat(num, den)


def weyl_moments(chi, dd, order, nmax):
    cts = []
    es = []
    p = {(0, 0): 1}
    for n in range(nmax + 1):
        ct = lp_ct_pair(p, dd)
        cts.append(ct)
        es.append(rat(ct, order) if order else None)
        if n < nmax:
            p = lp_mul(p, chi)
    return cts, es


def char_sym(ws, k):
    out = {}
    for combo in itertools.combinations_with_replacement(range(len(ws)), k):
        e = (sum(ws[i][0] for i in combo), sum(ws[i][1] for i in combo))
        out[e] = out.get(e, 0) + 1
    return lp_clean(out)


def char_wedge(ws, k):
    out = {}
    for combo in itertools.combinations(range(len(ws)), k):
        e = (sum(ws[i][0] for i in combo), sum(ws[i][1] for i in combo))
        out[e] = out.get(e, 0) + 1
    return lp_clean(out)


# ---------------------------------------------------------------- Brauer-Klimyk


def dominate(v):
    # reflect v into the closed dominant chamber by simple reflections; the sign is det of the word
    sign = 1
    for _ in range(10000):
        p0 = simple_pair(v, 0)
        p1 = simple_pair(v, 1)
        if p0 < 0:
            v = (v[0] - p0 * SIMPLE[0][0], v[1] - p0 * SIMPLE[0][1])
            sign = -sign
        elif p1 < 0:
            v = (v[0] - p1 * SIMPLE[1][0], v[1] - p1 * SIMPLE[1][1])
            sign = -sign
        else:
            return v, sign, (p0 == 0 or p1 == 0)
    raise ValueError("dominate did not terminate")


def bk_powers(weights, rho, nmax):
    # V(lam) tensor V = sum over weights mu of V of sign(w) V(w(lam + mu + rho) - rho), walls dropped
    cur = {(0, 0): 1}
    trivial = [1]
    negatives = 0
    sizes = [1]
    for _ in range(nmax):
        nxt = {}
        for lam, m in cur.items():
            for mu in weights:
                v = (lam[0] + mu[0] + rho[0], lam[1] + mu[1] + rho[1])
                d, sgn, wall = dominate(v)
                if wall:
                    continue
                nu = (d[0] - rho[0], d[1] - rho[1])
                nxt[nu] = nxt.get(nu, 0) + sgn * m
        cur = lp_clean(nxt)
        negatives = negatives + sum(1 for c in cur.values() if c < 0)
        trivial.append(cur.get((0, 0), 0))
        sizes.append(len(cur))
    return trivial, negatives, sizes


# ---------------------------------------------------------------- the G2 suite (also run on the plants)


def analyze(short_pos, long_pos):
    checks = []

    def note(name, ok, detail):
        checks.append((name, bool(ok), detail))

    pos = list(short_pos) + list(long_pos)
    roots = pos + [v_neg(a) for a in pos]
    root_set = set(roots)
    sn = [ip(a, a) for a in short_pos]
    ln = [ip(a, a) for a in long_pos]
    note("g2_root_count_12_short_norm_2_long_norm_6",
         len(root_set) == 12 and all(x == 2 for x in sn) and all(x == 6 for x in ln),
         [("distinct_roots", len(root_set)), ("short_norms", join(sn)), ("long_norms", join(ln))])

    cartan_bad = 0
    closure_bad = 0
    for a in roots:
        for b in roots:
            if coroot_pair(b, a) is None:
                cartan_bad = cartan_bad + 1
            img = reflect(b, a)
            if img is None or img not in root_set:
                closure_bad = closure_bad + 1
    note("g2_root_axioms_cartan_integers_and_reflection_closure", cartan_bad == 0 and closure_bad == 0,
         [("pairs_read", len(roots) * len(roots)), ("cartan_nonintegral", cartan_bad),
          ("reflection_leaves_root_set", closure_bad)])

    gens = []
    nonint = 0
    for a in pos:
        m = reflection_matrix(a)
        if m is None:
            nonint = nonint + 1
        else:
            gens.append(m)
    group, finite = generate_group(gens)
    order = len(group) if finite else 0
    dets = [m_det(m) for m in group]
    rot = dets.count(1)
    refl = dets.count(-1)
    gram_ok = all(preserves_gram(m) for m in group)
    note("g2_weyl_group_order_12_from_root_reflections",
         finite and nonint == 0 and order == 12 and rot == 6 and refl == 6 and gram_ok,
         [("reflections_used", len(gens)), ("reflections_nonintegral", nonint), ("finite_under_cap", 1 if finite else 0),
          ("order", order), ("det_plus", rot), ("det_minus", refl), ("preserve_gram", 1 if gram_ok else 0)])

    P, pdet, w1, w2 = fundamental_weights()
    lat_ok = None not in w1 and None not in w2
    wdet = (w1[0] * w2[1] - w2[0] * w1[1]) if lat_ok else 0
    note("g2_weight_lattice_equals_root_lattice", lat_ok and (wdet == 1 or wdet == -1),
         [("pairing_matrix", join(P[0]) + ";" + join(P[1])), ("pairing_det", pdet),
          ("w1", vec_str(w1)), ("w2", vec_str(w2)), ("det_w1_w2", wdet)])
    rho_f = v_add(w1, w2)

    s0 = sum(a[0] for a in pos)
    s1 = sum(a[1] for a in pos)
    h0 = exact_div(s0, 2)
    h1 = exact_div(s1, 2)
    rho_h = None if (h0 is None or h1 is None) else (h0, h1)
    note("g2_rho_half_sum_integral_and_equals_w1_plus_w2", rho_h == rho_f,
         [("sum_positive_roots", vec_str((s0, s1))), ("rho_half_sum", "NONE" if rho_h is None else vec_str(rho_h)),
          ("rho_w1_plus_w2", vec_str(rho_f))])
    rho = rho_h if rho_h is not None else rho_f

    delta0 = {(0, 0): 1}
    for a in pos:
        delta0 = lp_mul(delta0, {(0, 0): 1, v_neg(a): -1})
    delta = lp_mul({rho: 1}, delta0)
    alt_rho = lp_alt(group, rho)
    note("g2_weyl_denominator_identity", rho_h is not None and delta == alt_rho,
         [("delta_terms", len(delta)), ("alternating_sum_terms", len(alt_rho))])

    dd = {(0, 0): 1}
    for a in pos:
        dd = lp_mul(dd, lp_mul({(0, 0): 1, v_neg(a): -1}, {(0, 0): 1, a: -1}))
    ddbar = lp_mul(delta, lp_conj(delta))
    ct_dd = dd.get((0, 0), 0)
    note("g2_delta_conj_delta_equals_root_product_and_CT_equals_W", ddbar == dd and order > 0 and ct_dd == order,
         [("root_product_terms", len(dd)), ("CT", ct_dd), ("W", order)])

    w7 = [(0, 0)] + list(short_pos) + [v_neg(a) for a in short_pos]
    w14 = [(0, 0), (0, 0)] + roots
    inv_bad = 0
    for m in group:
        for ws in (w7, w14):
            if sorted(m_apply(m, x) for x in ws) != sorted(ws):
                inv_bad = inv_bad + 1
    note("g2_weight_multisets_7_and_14_W_invariant", inv_bad == 0 and order > 0,
         [("group_elements_read", len(group)), ("failures", inv_bad)])

    chi7 = lp_from_weights(w7)
    chi14 = lp_from_weights(w14)
    c7 = lp_mul(chi7, delta) == lp_alt(group, v_add(w1, rho))
    c14 = lp_mul(chi14, delta) == lp_alt(group, v_add(w2, rho))
    note("g2_weyl_character_formula_7_highest_weight_w1", c7, [("highest_weight", vec_str(w1))])
    note("g2_weyl_character_formula_14_highest_weight_w2", c14, [("highest_weight", vec_str(w2))])

    d7 = weyl_dim(w1, rho_f, pos)
    d14 = weyl_dim(w2, rho_f, pos)
    note("g2_weyl_dimension_formula_7_and_14", d7 == (7, 1) and d14 == (14, 1) and len(w7) == 7 and len(w14) == 14,
         [("dim_w1", r_str(d7)), ("dim_w2", r_str(d14)), ("weights_7", len(w7)), ("weights_14", len(w14))])

    ct7, e7 = weyl_moments(chi7, dd, order, NMAX)
    ct14, e14 = weyl_moments(chi14, dd, order, NMAX)
    note("g2_route1_moments_integral", all(e is not None and e[1] == 1 for e in e7 + e14),
         [("n_read", NMAX + 1)])
    note("g2_route1_normalisation_n0_is_1_n1_is_0",
         e7[0] == R1 and e14[0] == R1 and e7[1] == R0 and e14[1] == R0,
         [("E7_0", r_str(e7[0])), ("E7_1", r_str(e7[1])), ("E14_0", r_str(e14[0])), ("E14_1", r_str(e14[1]))])
    note("g2_E_chi7_squared_is_1_irreducible_real", e7[2] == R1, [("E7_2", r_str(e7[2]))])
    note("g2_E_chi14_squared_is_1_irreducible_real", e14[2] == R1, [("E14_2", r_str(e14[2]))])

    bk7, neg7, sz7 = bk_powers(w7, rho_f, NMAX)
    bk14, neg14, sz14 = bk_powers(w14, rho_f, NMAX)
    ag7 = [1 if e7[n] == (bk7[n], 1) else 0 for n in range(NMAX + 1)]
    ag14 = [1 if e14[n] == (bk14[n], 1) else 0 for n in range(NMAX + 1)]
    note("g2_routes_agree_chi7_n_le_6", all(ag7[:ROUTE2_MIN + 1]),
         [("route1", join(r_str(e) for e in e7[:ROUTE2_MIN + 1])), ("route2", join(bk7[:ROUTE2_MIN + 1]))])
    note("g2_routes_agree_chi14_n_le_6", all(ag14[:ROUTE2_MIN + 1]),
         [("route1", join(r_str(e) for e in e14[:ROUTE2_MIN + 1])), ("route2", join(bk14[:ROUTE2_MIN + 1]))])
    note("g2_brauer_klimyk_multiplicities_nonnegative", neg7 == 0 and neg14 == 0,
         [("negative_7", neg7), ("negative_14", neg14)])

    sq = lp_add(char_sym(w7, 2), char_wedge(w7, 2))
    note("g2_sym2_plus_wedge2_equals_chi7_squared", sq == lp_mul(chi7, chi7), [("terms", len(sq))])

    data = {
        "pos": pos, "roots": roots, "group": group, "order": order, "rho": rho_f, "w1": w1, "w2": w2,
        "w7": w7, "w14": w14, "dd": dd, "delta": delta, "ct7": ct7, "e7": e7, "ct14": ct14, "e14": e14,
        "bk7": bk7, "bk14": bk14, "ag7": ag7, "ag14": ag14, "sz7": sz7, "sz14": sz14,
    }
    return checks, data


# ---------------------------------------------------------------- main


def main():
    if not (sys.flags.isolated and sys.flags.dont_write_bytecode):
        print("HALT|reason=run as python3 -I -B|json=0")
        return 2
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    rows = []
    checks = []

    def check(name, ok, detail):
        checks.append((name, bool(ok), detail))
        return bool(ok)

    rows.append(row("HEAD", [
        ("schema", SCHEMA), ("date", DATE), ("author", AUTHOR), ("kind", "PREDICTION_before_any_measurement"),
        ("measurement_status", "NOT_RUN"), ("group_elements_built", 0),
        ("group_elements_built_reads", "elements of the frame group H and of G2 built by this script"),
        ("octonion_table_built", 0), ("weyl_group_matrices_built", 12), ("points_built", 14),
        ("script", SCRIPT_NAME), ("script_sha256", script_sha), ("run", "python3 -I -B " + SCRIPT_NAME),
        ("float_used", 0),
        ("sources", "NAMED: Weyl integration and character formulas (Bourbaki Lie ch. IX; Brocker and tom Dieck "
                    "ch. VI); Brauer-Klimyk rule (Humphreys 1972 ex. 24.9); sphere moments of the uniform measure "
                    "(Folland 2001, How to integrate a polynomial over a sphere); Cohen and Wales 1983, Finite "
                    "subgroups of G2(C), Comm. Algebra 11; Griess 1995, Basic conjugacy theorems for G2, Invent. "
                    "Math. 121; OEIS A059710"),
    ]))
    rows.append(row("FORMAT", [
        ("rational", "num:den reduced, den positive"),
        ("vector", "(x,y) means x a1 + y a2 in the simple-root basis"),
        ("laurent", "dict from exponent pairs (x,y) to integer coefficients; e^(x,y) is a torus character"),
        ("tensor", "written (x)"), ("claim_tags", "MEASURED computed by this script; NAMED stated and coherent, "
                                                  "not run here; CONJECTURE"),
        ("prediction_tag", "status PREDICTION marks a statement about H to be tested later against H"),
    ]))

    # ---- ORBIT block
    pts = orbit_points()
    pset = set(pts)
    unit = all(sum(x * x for x in p) == 1 for p in pts)
    anti = all(tuple(-x for x in p) in pset for p in pts)
    check("orbit_14_distinct_unit_points_antipodal", len(pts) == 14 and len(pset) == 14 and unit and anti,
          [("points", len(pts)), ("distinct", len(pset)), ("unit", 1 if unit else 0), ("antipodal", 1 if anti else 0)])
    rows.append(row("ORBIT", [
        ("points", "+-e_1..+-e_7 in R^7, each coordinate 0 or +-1, one nonzero per point"), ("count", len(pts)),
        ("measure", "uniform on the 14 points"), ("compare", "uniform measure on S^6"), ("claim", "MEASURED"),
    ]))
    for k in range(5):
        so = q_power_mean(orbit_moment, pts, k)
        ss = q_power_mean(sphere_moment, None, k)
        check("sphere_and_orbit_E_q_power_" + str(k) + "_is_1", so == R1 and ss == R1,
              [("orbit", r_str(so)), ("sphere", r_str(ss)), ("q", "sum of squares")])

    two_var = {}
    for a in range(ORBIT_DEG + 1):
        for b in range(ORBIT_DEG + 1 - a):
            alpha = (a, b, 0, 0, 0, 0, 0)
            o = orbit_moment(pts, alpha)
            s = sphere_moment(alpha)
            d = r_sub(o, s)
            two_var[(a, b)] = (o, s, d)
            rows.append(row("ORBIT_MOMENT", [
                ("a", a), ("b", b), ("degree", a + b), ("E_orbit_X1a_X2b", r_str(o)), ("E_sphere_Y1a_Y2b", r_str(s)),
                ("diff_orbit_minus_sphere", r_str(d)), ("agree", 1 if d == R0 else 0), ("claim", "MEASURED"),
                ("sphere_formula", "NAMED"),
            ]))

    per_degree = []
    for deg in range(ORBIT_DEG + 1):
        total = 0
        bad = 0
        pats = {}
        for alpha in exponent_vectors(deg):
            total = total + 1
            if orbit_moment(pts, alpha) != sphere_moment(alpha):
                bad = bad + 1
                pk = pattern(alpha)
                pats[pk] = pats.get(pk, 0) + 1
        per_degree.append((deg, total, bad))
        rows.append(row("DESIGN", [
            ("degree", deg), ("monomials_in_7_coordinates", total), ("disagree", bad),
            ("failing_patterns", join((pk + "#" + str(pats[pk]) for pk in sorted(pats)), ",") if pats else "none"),
            ("claim", "MEASURED"),
        ]))
    t = -1
    for deg, total, bad in per_degree:
        if bad:
            break
        t = deg
    first_fail = t + 1 if t < ORBIT_DEG else "none_up_to_8"
    read = sum(x[1] for x in per_degree)
    rows.append(row("DESIGN_DEGREE", [
        ("t", t), ("first_failing_degree", first_fail), ("monomials_read", read), ("task_expected", TASK_DESIGN_DEGREE),
        ("match", 1 if t == TASK_DESIGN_DEGREE else 0),
        ("reads", "every moment of degree up to t agrees with S^6 and one of degree t+1 does not"),
        ("claim", "MEASURED"),
    ]))
    check("orbit_design_degree_read_over_all_monomials", read == 6435 and t >= 0,
          [("monomials_read", read), ("expected_monomials", 6435), ("t", t)])

    d40 = two_var[(4, 0)][2]
    d22 = two_var[(2, 2)][2]
    rows.append(row("FIRST_DIFF", [
        ("moment", "E[X_1^4] minus E[Y_1^4]"), ("orbit", r_str(two_var[(4, 0)][0])), ("sphere", r_str(two_var[(4, 0)][1])),
        ("diff", r_str(d40)), ("task_expected", r_str(TASK_X4_DIFF)), ("match", 1 if d40 == TASK_X4_DIFF else 0),
        ("forced", 0), ("claim", "MEASURED"),
    ]))
    rows.append(row("FIRST_DIFF", [
        ("moment", "E[X_1^2 X_2^2] minus E[Y_1^2 Y_2^2]"), ("orbit", r_str(two_var[(2, 2)][0])),
        ("sphere", r_str(two_var[(2, 2)][1])), ("diff", r_str(d22)), ("claim", "MEASURED"),
    ]))
    for (a, b) in ((6, 0), (4, 2), (8, 0), (6, 2), (4, 4)):
        o, s, d = two_var[(a, b)]
        rows.append(row("NEXT_DIFF", [
            ("a", a), ("b", b), ("orbit", r_str(o)), ("sphere", r_str(s)), ("diff", r_str(d)), ("claim", "MEASURED"),
        ]))

    # quartic sum x_i^4 against q^2: ratio at three points of R^7 (one ray each)
    def quartic_ratio(v):
        s4 = sum(x * x * x * x for x in v)
        s2 = sum(x * x for x in v)
        return rat(s4, s2 * s2)

    e1 = (1, 0, 0, 0, 0, 0, 0)
    e12 = (1, 1, 0, 0, 0, 0, 0)
    ones = (1, 1, 1, 1, 1, 1, 1)
    qr1 = quartic_ratio(e1)
    qr2 = quartic_ratio(e12)
    qr3 = quartic_ratio(ones)
    nonconst = qr1 != qr2 or qr2 != qr3
    om4 = R0
    sm4 = R0
    for i in range(DIM):
        alpha = tuple(4 if j == i else 0 for j in range(DIM))
        om4 = r_add(om4, orbit_moment(pts, alpha))
        sm4 = r_add(sm4, sphere_moment(alpha))
    rows.append(row("QUARTIC", [
        ("form", "sum_i x_i^4"), ("ratio_to_q_squared_at_e1", r_str(qr1)), ("at_e1_plus_e2", r_str(qr2)),
        ("at_all_ones", r_str(qr3)), ("ratio_constant_on_S6", 0 if nonconst else 1),
        ("orbit_mean", r_str(om4)), ("sphere_mean", r_str(sm4)),
        ("reading", "a G2-fixed quartic is constant on S^6 since G2 is transitive on S^6 (NAMED); this ratio is not, "
                    "so sum_i x_i^4 is not G2-fixed; H permutes the axes up to sign, so it is H-fixed (NAMED)"),
        ("claim", "MEASURED values; NAMED reading"),
    ]))
    check("quartic_not_constant_on_S6", nonconst, [("ratios", r_str(qr1) + "," + r_str(qr2) + "," + r_str(qr3))])

    # ---- G2 block
    g2checks, g = analyze(SHORT_POS, LONG_POS)
    checks.extend(g2checks)
    rows.append(row("BASIS", [
        ("lattice", "weight lattice of G2, equal to its root lattice (checked)"),
        ("basis", "a1 short simple root, a2 long simple root"), ("gram", "2,-3;-3,6"),
        ("positive_short", join((vec_str(a) for a in SHORT_POS), ";")),
        ("positive_long", join((vec_str(a) for a in LONG_POS), ";")),
        ("rho", vec_str(g["rho"])), ("w1", vec_str(g["w1"])), ("w2", vec_str(g["w2"])),
        ("weights_7", "0 and the 6 short roots"), ("weights_14", "0 twice and the 12 roots"),
        ("W_order", g["order"]), ("W_generated_from", "reflections in the 6 positive roots, closed under products"),
        ("delta_conj_delta_terms", len(g["dd"])), ("delta_conj_delta_CT", g["dd"].get((0, 0), 0)),
        ("claim", "MEASURED"),
    ]))
    for n in range(NMAX + 1):
        rows.append(row("G2_MOMENT", [
            ("rep", 7), ("n", n), ("CT", g["ct7"][n]), ("W", g["order"]), ("E_G2_chi7_pow_n", r_str(g["e7"][n])),
            ("route2_trivial_multiplicity", g["bk7"][n]), ("route2_dominant_weights", g["sz7"][n]),
            ("routes_agree", g["ag7"][n]), ("claim", "MEASURED"),
            ("route1", "Weyl integration (NAMED formula)"), ("route2", "Brauer-Klimyk tensoring (NAMED rule)"),
        ]))
    for n in range(NMAX + 1):
        rows.append(row("G2_MOMENT", [
            ("rep", 14), ("n", n), ("CT", g["ct14"][n]), ("W", g["order"]), ("E_G2_chi14_pow_n", r_str(g["e14"][n])),
            ("route2_trivial_multiplicity", g["bk14"][n]), ("route2_dominant_weights", g["sz14"][n]),
            ("routes_agree", g["ag14"][n]), ("claim", "MEASURED"),
            ("route1", "Weyl integration (NAMED formula)"), ("route2", "Brauer-Klimyk tensoring (NAMED rule)"),
        ]))

    e7 = g["e7"]
    e14 = g["e14"]
    e7_ints = [e[0] for e in e7]
    e14_ints = [e[0] for e in e14]

    sym_dims = []
    for k in range(SYM_MAX + 1):
        cs = char_sym(g["w7"], k)
        ct = lp_ct_pair(cs, g["dd"])
        e = rat(ct, g["order"])
        sym_dims.append(e)
        rows.append(row("G2_FIXED", [
            ("space", "Sym^" + str(k) + " of 7"), ("dim_space", sum(cs.values())), ("CT", ct),
            ("dim_G2_fixed", r_str(e)), ("claim", "MEASURED"),
        ]))
    wedge_dims = []
    for k in range(DIM + 1):
        cw = char_wedge(g["w7"], k)
        ct = lp_ct_pair(cw, g["dd"])
        e = rat(ct, g["order"])
        wedge_dims.append(e)
        rows.append(row("G2_FIXED", [
            ("space", "Wedge^" + str(k) + " of 7"), ("dim_space", sum(cw.values())), ("CT", ct),
            ("dim_G2_fixed", r_str(e)), ("claim", "MEASURED"),
        ]))
    sym_ok = all(sym_dims[k] == (R1 if (k & 1) == 0 else R0) for k in range(SYM_MAX + 1))
    check("g2_fixed_polynomials_on_7_are_powers_of_q_up_to_degree_8", sym_ok,
          [("sym_dims", join(r_str(e) for e in sym_dims)),
           ("reads", "matches the NAMED fact that G2-fixed polynomials on R^7 are polynomials in q")])
    check("g2_fixed_wedge3_is_1_the_3_form", wedge_dims[3] == R1, [("wedge_dims", join(r_str(e) for e in wedge_dims))])

    cmp_bad = 0
    m_task = 1 if e7_ints[:len(TASK_CHI7)] == list(TASK_CHI7) and all(e[1] == 1 for e in e7) else 0
    m_oeis = 1 if e7_ints[:len(OEIS_A059710)] == list(OEIS_A059710) and all(e[1] == 1 for e in e7) else 0
    cmp_bad = cmp_bad + (1 - m_task) + (1 - m_oeis)
    rows.append(row("COMPARE", [
        ("what", "E_G2[chi_7^n], n 0..6, against the values stated in the task"),
        ("computed", join(e7_ints[:len(TASK_CHI7)])), ("stated", join(TASK_CHI7)), ("match", m_task), ("forced", 0),
        ("claim", "MEASURED"),
    ]))
    rows.append(row("COMPARE", [
        ("what", "E_G2[chi_7^n], n 0..8, against OEIS A059710"), ("computed", join(e7_ints)),
        ("listed", join(OEIS_A059710)), ("match", m_oeis), ("forced", 0),
        ("source_status", "OPENED: oeis.org A059710 read by the prediction agent through a web fetch on 2026-11-09; "
                          "title read there: dimension of space of invariants of n-th tensor power of the "
                          "7-dimensional irreducible representation of G_2; terms typed into this script as data"),
        ("claim", "MEASURED comparison; NAMED source"),
    ]))
    rows.append(row("COMPARE", [
        ("what", "E_G2[chi_14^n], n 0..8"), ("computed", join(e14_ints)), ("stated", "none"),
        ("claim", "MEASURED"),
    ]))

    # ---- FRAME block (predictions about H; nothing of H is built)
    n_order = 8
    q_order = 168
    h_order = 1344
    rows.append(row("FRAME_STRUCTURE", [
        ("id", "F0"), ("object", "H, the subgroup of G2 that permutes the 14 points +-e_1..+-e_7 and keeps the "
                                 "octonion product"),
        ("order", h_order), ("order_is_8_times_168", 1 if n_order * q_order == h_order else 0),
        ("N", "H meet the diagonal sign changes, predicted of order 8, isomorphic to 2^3; it acts on axis k by a "
              "character p_k, and p_1..p_7 are the 7 distinct nontrivial characters of N"),
        ("H_mod_N", "predicted GL(3,2) of order 168, acting on the 7 axes as on the points of the Fano plane whose "
                    "lines are the octonion triples"),
        ("orbit_of_e_1", "predicted all 14 points, so the stabiliser of e_1 has order 96"),
        ("stabiliser_check", 1 if 14 * 96 == h_order else 0),
        ("status", "PREDICTION"), ("claim", "NAMED"),
    ]))
    rows.append(row("FRAME_EXTENSION", [
        ("id", "F1"), ("named", "H is the non-split extension 2^3.GL(3,2)"),
        ("sources", "Cohen and Wales 1983, Finite subgroups of G2(C), Comm. Algebra 11; Griess 1995, Basic "
                    "conjugacy theorems for G2, Invent. Math. 121"),
        ("source_status", "UNVERIFIED: not opened by this agent"),
        ("prediction", "no subgroup of H of order 168 meets the sign-change kernel N trivially, so N has no "
                       "complement in H"),
        ("test", "count the subgroups K of H of order 168 with K meet N trivial; predicted count 0"),
        ("status", "PREDICTION"), ("claim", "NAMED"),
    ]))

    reasons7 = {
        0: "normalisation",
        1: "every p_k is nontrivial, so N fixes no nonzero vector of 7 and neither does H",
        2: "7 restricted to N is p_1 + ... + p_7 with the p_k distinct, and H permutes the 7 axes transitively, "
           "so the commutant of H on 7 is the scalars (Clifford); 7 is real, so the mean of chi_7^2 is "
           "the inner product of chi_7 with itself, which is 1",
        3: "the N-fixed tensors in 7 (x) 7 (x) 7 are spanned by e_i (x) e_j (x) e_k with p_i p_j p_k trivial, which "
           "holds exactly when i j k are distinct points of one Fano line: 7 lines times 6 orders gives 42 lines of "
           "tensors; H permutes them transitively, as GL(3,2) is 2-transitive on points and i j fix k, so they form "
           "a transitive monomial representation and the H-fixed space has dimension at most 1; the G2-fixed "
           "3-form lies in it, so it is exactly the G2 value 1",
    }
    for n in range(4):
        rows.append(row("FRAME_CHI7", [
            ("id", "F2." + str(n)), ("n", n), ("E_G2", r_str(e7[n])), ("relation", "EQUAL"),
            ("E_H_predicted", r_str(e7[n])), ("reason", reasons7[n]),
            ("test", "the mean over the 1344 elements of H of chi_7(h)^n"),
            ("status", "PREDICTION"), ("claim", "NAMED"),
        ]))
    # dim U + dim S - dim (U meet S): U the G2-fixed tensors, S the H-fixed symmetric tensors (q^2 and sum x_i^4
    # at least), U meet S the G2-fixed symmetric tensors
    s4 = sym_dims[4][0]
    s_h_lower = s4 + 1
    lb4 = e7_ints[4] + s_h_lower - s4
    rows.append(row("FRAME_CHI7", [
        ("id", "F2.4"), ("n", 4), ("E_G2", r_str(e7[4])), ("relation", "GREATER_OR_EQUAL"),
        ("E_H_lower_bound", lb4), ("E_H_exact", "NOT_PREDICTED"),
        ("reason", "an extra invariant in Sym^4 of the 7 gives an extra invariant in 7 (x)4: Sym^4 of 7 has G2-fixed "
                   "dimension " + r_str(sym_dims[4]) + " (G2_FIXED row, q^2 only); sum_i x_i^4 is H-fixed and not a "
                   "multiple of q^2 (QUARTIC row), so the H-fixed symmetric tensors S have dimension at least "
                   + str(s_h_lower) + "; inside 7 (x)4 the G2-fixed space U meets S in the G2-fixed symmetric "
                   "tensors, so the H-fixed space holds U + S, of dimension at least " + str(e7_ints[4]) + " + "
                   + str(s_h_lower) + " - " + str(s4)),
        ("test", "the mean over H of chi_7(h)^4, predicted at least " + str(lb4)),
        ("status", "PREDICTION"), ("claim", "NAMED"),
    ]))
    for n in range(5, NMAX + 1):
        even = (n & 1) == 0
        lb = e7_ints[n] + (1 if even else 0)
        rows.append(row("FRAME_CHI7", [
            ("id", "F2." + str(n)), ("n", n), ("E_G2", r_str(e7[n])), ("relation", "GREATER_OR_EQUAL"),
            ("E_H_lower_bound", lb), ("E_H_exact", "CANNOT_PREDICT"),
            ("reason", ("G2-fixed tensors are H-fixed; for even n, sum_i x_i^4 times q to the power (n-4) over 2 is "
                        "H-fixed and not a multiple of q to the power n over 2, and Sym^n of 7 has G2-fixed dimension "
                        + r_str(sym_dims[n]) + " (G2_FIXED row), so the count of F2.4 adds 1")
                       if even else "G2-fixed tensors are H-fixed"),
            ("status", "PREDICTION"), ("claim", "NAMED"),
        ]))

    reasons14 = {
        0: "normalisation",
        1: "14 sits in Wedge^2 of 7; e_i wedge e_j carries the N-character p_i p_j, which is p_k for the third point "
           "k of their line and so nontrivial; Wedge^2 of 7 has no N-fixed vector, hence 14 has no H-fixed vector",
        2: "Wedge^2 of 7 is 7 + 14 for G2 (NAMED; REFUTE.hbp lambda2_split_7_14 measured it) and restricts to N as "
           "3 (p_1 + ... + p_7), so 14 restricts as 2 (p_1 + ... + p_7); the p_k part M_k of Wedge^2 is spanned by "
           "the 3 wedges over the 3 lines through k; the inertia group I_k (order 192, over the point stabiliser S_4) "
           "permutes those 3 lines 2-transitively through S_4 onto GL(2,2), so M_k is monomial with 2 double cosets "
           "and its self inner product is at most 2; M_k is (7)_k + (14)_k with (7)_k of dimension 1, which forces "
           "(14)_k irreducible for I_k; by Clifford 14 is induced from it and irreducible for H, and real, so the "
           "mean of chi_14^2 is 1",
    }
    for n in range(3):
        rows.append(row("FRAME_CHI14", [
            ("id", "F3." + str(n)), ("n", n), ("E_G2", r_str(e14[n])), ("relation", "EQUAL"),
            ("E_H_predicted", r_str(e14[n])), ("reason", reasons14[n]),
            ("test", "the mean over H of chi_14(h)^n"),
            ("status", "PREDICTION"), ("claim", "NAMED"),
        ]))
    rows.append(row("FRAME_CHI14_N2_PREDICTABLE", [
        ("answer", "YES"), ("value", r_str(e14[2])),
        ("rests_on", "F0 structure (N of order 8 with characters p_1..p_7, H mod N GL(3,2) on the Fano points) and "
                     "Wedge^2 of 7 equal to 7 + 14"),
        ("status", "PREDICTION"), ("claim", "NAMED"),
    ]))
    for n in range(3, NMAX + 1):
        rows.append(row("FRAME_CHI14", [
            ("id", "F3." + str(n)), ("n", n), ("E_G2", r_str(e14[n])), ("relation", "GREATER_OR_EQUAL"),
            ("E_H_lower_bound", e14_ints[n]), ("E_H_exact", "CANNOT_PREDICT"),
            ("reason", "G2-fixed tensors are H-fixed"), ("status", "PREDICTION"), ("claim", "NAMED"),
        ]))

    # ---- planted wrong roots: each must break at least one check of the same suite
    for name, sp, lp in PLANTS:
        pchecks, _ = analyze(sp, lp)
        broken = [c[0] for c in pchecks if not c[1]]
        rows.append(row("PLANTED", [
            ("plant", name), ("positive_short", join((vec_str(a) for a in sp), ";")),
            ("positive_long", join((vec_str(a) for a in lp), ";")), ("suite_checks", len(pchecks)),
            ("checks_broken", len(broken)), ("broken", join(broken, ";") if broken else "none"), ("claim", "MEASURED"),
        ]))
        check("planted_" + name + "_breaks_a_check", len(broken) >= 1, [("checks_broken", len(broken))])

    failed = 0
    for name, ok, detail in checks:
        if not ok:
            failed = failed + 1
        rows.append(row("SELF_CHECK", [("name", name), ("status", "PASS" if ok else "FAIL")] + list(detail)
                        + [("claim", "MEASURED")]))
    rows.append(row("END", [
        ("rows", len(rows) + 1), ("checks", len(checks)), ("checks_failed", failed),
        ("compare_mismatches", cmp_bad), ("exit", 1 if failed else 0), ("measurement_status", "NOT_RUN"),
        ("group_elements_built", 0),
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
