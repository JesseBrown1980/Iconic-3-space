#!/usr/bin/env python3
"""
measure_8d.py -- the 8D Spin(7) case of the frozen spinor-symbol prediction, measured from definitions.

Seat: opus55 on acer. Integer only: no float, no division operator, no format strings, math.gcd only.
Every identity in xi is an identity of integer polynomials (dict monomial -> int, a monomial being the sorted
tuple of its variable indices), checked coefficient by coefficient over the full monomial basis of its
multidegree; a coefficient outside that basis counts as a failure too. Sampled covectors appear only as controls.

Read before the measured block was written: the CONV and FORMAT rows of PREDICTION-SPINOR-SYMBOL.hbp, nothing else
of the frozen files, no other seat's script or receipt.

Variables: xi_t is variable t (t = 0..7), a_t and x_t are variable 8 + t.
Octonions: e0 = 1, imaginary units e1..e7, e_i e_(i+1) = e_(i+3) on labels 1..7 mod 7, e_i^2 = -1.

Usage: python3 -I -B measure_8d.py OUTDIR                     measured block only
       python3 -I -B measure_8d.py OUTDIR PRED.hbp ADD.hbp    measured block, then COMPARE against the frozen rows
"""
import sys
import os
import hashlib
import itertools
from math import gcd

OUT_NAME = "MEASURE-8D.hbp"
ROWS = []
TALLY = {"pass": 0, "fail": 0, "ctrl_ok": 0, "ctrl_bad": 0}


def clean(s):
    return s.replace("|", ":").replace("\n", " ")


def row(tag, pairs):
    parts = [tag]
    for k, v in pairs:
        parts.append(str(k) + "=" + clean(str(v)))
    parts.append("json=0")
    ROWS.append("|".join(parts))


RESULTS = {}


def meas(rid, name, ok, pairs):
    TALLY["pass" if ok else "fail"] += 1
    RESULTS[rid] = dict(pairs)
    RESULTS[rid]["ok"] = ok
    row("M8", [("id", rid), ("name", name)] + pairs + [("status", "PASS" if ok else "FAIL"), ("claim", "MEASURED")])


def ctrl(rid, name, expect_fail, observed_fail, pairs):
    good = (expect_fail == observed_fail)
    RESULTS[rid] = dict(pairs)
    RESULTS[rid]["observed_fail"] = observed_fail
    RESULTS[rid]["ok"] = good
    TALLY["ctrl_ok" if good else "ctrl_bad"] += 1
    row("CTRL", [("id", rid), ("name", name)] + pairs +
        [("expect", "FAIL" if expect_fail else "PASS"), ("observed", "FAIL" if observed_fail else "PASS"),
         ("control_ok", 1 if good else 0), ("claim", "MEASURED")])


# ---------------------------------------------------------------- integer polynomials
def p_add_into(r, q, k=1):
    for m, c in q.items():
        v = r.get(m, 0) + c * k
        if v:
            r[m] = v
        elif m in r:
            del r[m]
    return r


def p_add(p, q):
    return p_add_into(dict(p), q, 1)


def p_sub(p, q):
    return p_add_into(dict(p), q, -1)


def p_scale(p, k):
    if k == 0:
        return {}
    return {m: c * k for m, c in p.items()}


def p_mul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2))
            v = r.get(m, 0) + c1 * c2
            if v:
                r[m] = v
            elif m in r:
                del r[m]
    return r


def p_var(i):
    return {(i,): 1}


def p_const(c):
    return {(): c} if c else {}


def p_eval(p, point):
    tot = 0
    for m, c in p.items():
        t = c
        for v in m:
            t = t * point[v]
        tot += t
    return tot


def basis(groups):
    parts = [list(itertools.combinations_with_replacement(vs, d)) for vs, d in groups]
    out = []
    for combo in itertools.product(*parts):
        m = ()
        for c in combo:
            m = m + c
        out.append(tuple(sorted(m)))
    return out


def zero_check(entries, monos):
    bset = set(monos)
    checked = 0
    failed = 0
    outside = 0
    for e in entries:
        for m in monos:
            checked += 1
            if e.get(m, 0):
                failed += 1
        for m in e:
            if m not in bset:
                outside += 1
    return checked, failed, outside


def holds(res):
    return res[1] == 0 and res[2] == 0


def cpairs(res):
    return [("coeffs_checked", res[0]), ("coeffs_failed", res[1]), ("coeffs_outside_basis", res[2])]


# ---------------------------------------------------------------- polynomial matrices
def pm_zero(n, m):
    return [[{} for _ in range(m)] for _ in range(n)]


def pm_mul(a, b):
    n = len(a)
    k = len(b)
    m = len(b[0])
    c = pm_zero(n, m)
    for i in range(n):
        for t in range(k):
            x = a[i][t]
            if not x:
                continue
            for j in range(m):
                y = b[t][j]
                if y:
                    p_add_into(c[i][j], p_mul(x, y))
    return c


def pm_lin(a, b, ka, kb):
    return [[p_add_into(p_scale(a[i][j], ka), b[i][j], kb) for j in range(len(a[0]))] for i in range(len(a))]


def pm_T(a):
    return [[dict(a[i][j]) for i in range(len(a))] for j in range(len(a[0]))]


def pm_diag(p, n):
    out = pm_zero(n, n)
    for i in range(n):
        out[i][i] = dict(p)
    return out


def pm_from_int(m):
    return [[p_const(x) for x in r] for r in m]


def flat(a):
    return [e for r in a for e in r]


# ---------------------------------------------------------------- integer matrices
def im_zero(n, m):
    return [[0] * m for _ in range(n)]


def im_ident(n):
    out = im_zero(n, n)
    for i in range(n):
        out[i][i] = 1
    return out


def im_mul(a, b):
    n = len(a)
    k = len(b)
    m = len(b[0])
    c = im_zero(n, m)
    for i in range(n):
        ai = a[i]
        ci = c[i]
        for t in range(k):
            x = ai[t]
            if x:
                bt = b[t]
                for j in range(m):
                    y = bt[j]
                    if y:
                        ci[j] += x * y
    return c


def im_lin(a, b, ka, kb):
    return [[ka * a[i][j] + kb * b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def im_T(a):
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def im_is_zero(a):
    return all(x == 0 for r in a for x in r)


def im_vec(a, v):
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def content_normalize(r):
    g = 0
    for x in r:
        if x:
            g = gcd(g, x)
    if g > 1:
        r = [divmod(x, g)[0] for x in r]
    return r


def rref_int(rows, ncols):
    a = [list(r) for r in rows]
    piv = []
    r = 0
    n = len(a)
    for c in range(ncols):
        if r == n:
            break
        pr = -1
        for i in range(r, n):
            if a[i][c]:
                pr = i
                break
        if pr < 0:
            continue
        a[r], a[pr] = a[pr], a[r]
        prow = a[r]
        pv = prow[c]
        for i in range(n):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = content_normalize([pv * x - f * y for x, y in zip(a[i], prow)])
        piv.append(c)
        r += 1
    return a[:r], piv


def rank_int(rows, ncols):
    if not rows:
        return 0
    return len(rref_int(rows, ncols)[1])


def kernel_int(rows, ncols):
    if not rows:
        return [[1 if t == f else 0 for t in range(ncols)] for f in range(ncols)]
    red, piv = rref_int(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]
    lcm = 1
    for i, c in enumerate(piv):
        p = abs(red[i][c])
        lcm = divmod(lcm * p, gcd(lcm, p))[0]
    out = []
    for f in free:
        v = [0] * ncols
        v[f] = lcm
        for i, c in enumerate(piv):
            v[c] = -red[i][f] * divmod(lcm, red[i][c])[0]
        out.append(content_normalize(v))
    for v in out:
        if any(x for x in im_vec(rows, v)):
            raise SystemExit("kernel_int self-check failed")
    return out


MOD_P = 2305843009213693951


def rank_mod_p(rows):
    pivots = {}
    rank = 0
    for r0 in rows:
        r = {k: v % MOD_P for k, v in r0.items() if v % MOD_P}
        while r:
            c = min(r)
            if c in pivots:
                pr = pivots[c]
                f = r[c]
                pv = pr[c]
                new = {k: (v * pv) % MOD_P for k, v in r.items()}
                for k, v in pr.items():
                    new[k] = (new.get(k, 0) - f * v) % MOD_P
                r = {k: v for k, v in new.items() if v}
            else:
                pivots[c] = r
                rank += 1
                break
    return rank


# ---------------------------------------------------------------- octonions
LINES_FROZEN = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def lines_from_rule():
    return [(i, (i % 7) + 1, ((i + 2) % 7) + 1) for i in range(1, 8)]


def build_table(lines):
    t = [[None] * 8 for _ in range(8)]
    for i in range(8):
        t[0][i] = (1, i)
        t[i][0] = (1, i)
    for i in range(1, 8):
        t[i][i] = (-1, 0)
    for (a, b, c) in lines:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            if t[x][y] is not None or t[y][x] is not None:
                return None
            t[x][y] = (1, z)
            t[y][x] = (-1, z)
    for i in range(8):
        for j in range(8):
            if t[i][j] is None:
                return None
    return t


def ovars(off):
    return [p_var(off + t) for t in range(8)]


def obasis(j):
    return [p_const(1) if t == j else {} for t in range(8)]


def o_mul(tab, x, y):
    z = [{} for _ in range(8)]
    for i in range(8):
        if not x[i]:
            continue
        for j in range(8):
            if not y[j]:
                continue
            s, k = tab[i][j]
            p_add_into(z[k], p_mul(x[i], y[j]), s)
    return z


def o_conj(x):
    return [dict(x[0])] + [p_scale(c, -1) for c in x[1:]]


def o_lin(x, y, kx, ky):
    return [p_add_into(p_scale(x[t], kx), y[t], ky) for t in range(8)]


def o_scalar(p, x):
    return [p_mul(p, c) for c in x]


def o_nsq(x):
    r = {}
    for c in x:
        p_add_into(r, p_mul(c, c))
    return r


def o_dot(x, y):
    r = {}
    for t in range(8):
        p_add_into(r, p_mul(x[t], y[t]))
    return r


def o_text(v):
    terms = []
    for t in range(8):
        if v[t]:
            terms.append(str(v[t]) + "e" + str(t))
    return "+".join(terms) if terms else "0"


def L_matrix(tab, x):
    cols = [o_mul(tab, x, obasis(j)) for j in range(8)]
    return [[cols[j][k] for j in range(8)] for k in range(8)]


def gamma_poly(tab, x, lower=-1):
    lm = L_matrix(tab, x)
    lc = L_matrix(tab, o_conj(x))
    g = pm_zero(16, 16)
    for r in range(8):
        for c in range(8):
            g[r][8 + c] = dict(lm[r][c])
            g[8 + r][c] = p_scale(lc[r][c], lower)
    return g


def coeff_int(g, var):
    return [[e.get((var,), 0) for e in r] for r in g]


# ---------------------------------------------------------------- exterior forms (dict sorted index tuple -> int)
def sort_sign(seq):
    n = len(seq)
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] == seq[j]:
                return 0, None
    s = 1
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                s = -s
    return s, tuple(sorted(seq))


def f_add_term(f, k, v):
    if not v:
        return
    w = f.get(k, 0) + v
    if w:
        f[k] = w
    elif k in f:
        del f[k]


def f_lin(f, g, kf, kg):
    out = {}
    for k, v in f.items():
        f_add_term(out, k, kf * v)
    for k, v in g.items():
        f_add_term(out, k, kg * v)
    return out


def f_wedge(f, g):
    out = {}
    for i, a in f.items():
        for j, b in g.items():
            s, k = sort_sign(i + j)
            if s:
                f_add_term(out, k, s * a * b)
    return out


def f_star(f, universe, sigma):
    out = {}
    for i, a in f.items():
        j = tuple(x for x in universe if x not in i)
        s, _ = sort_sign(i + j)
        f_add_term(out, j, sigma * s * a)
    return out


def f_deriv(rm, f, n):
    out = {}
    for i, a in f.items():
        for m in range(len(i)):
            b = i[m]
            for t in range(n):
                r = rm[t][b]
                if r:
                    s, k = sort_sign(i[:m] + (t,) + i[m + 1:])
                    if s:
                        f_add_term(out, k, s * a * r)
    return out


def f_vec(f, blist):
    idx = {b: n for n, b in enumerate(blist)}
    v = [0] * len(blist)
    for k, c in f.items():
        if k not in idx:
            raise SystemExit("form outside basis")
        v[idx[k]] = c
    return v


def op_matrix(func, in_basis, out_basis):
    cols = [f_vec(func({b: 1}), out_basis) for b in in_basis]
    return [[cols[j][i] for j in range(len(in_basis))] for i in range(len(out_basis))]


PAIRS8 = list(itertools.combinations(range(8), 2))
PAIRS7 = list(itertools.combinations(range(1, 8), 2))
TRIPLES8 = list(itertools.combinations(range(8), 3))
TRIPLES7 = list(itertools.combinations(range(1, 8), 3))
QUADS8 = list(itertools.combinations(range(8), 4))
U7 = tuple(range(1, 8))
U8 = tuple(range(8))


def R_of(w):
    rm = im_zero(8, 8)
    for p, (i, j) in enumerate(PAIRS8):
        if w[p]:
            rm[j][i] += w[p]
            rm[i][j] -= w[p]
    return rm


def w_of(rm):
    return [rm[j][i] for (i, j) in PAIRS8]


def is_skew(rm):
    return all(rm[i][j] == -rm[j][i] for i in range(8) for j in range(8))


def unit(n, k):
    return [1 if t == k else 0 for t in range(n)]


def rot_vec(rm, x):
    return [_rot_entry(rm, x, t) for t in range(8)]


def _rot_entry(rm, x, t):
    acc = {}
    for b in range(8):
        if rm[t][b]:
            p_add_into(acc, x[b], rm[t][b])
    return acc


def poly_text(p):
    if not p:
        return "0"
    out = []
    for m in sorted(p):
        mono = "*".join(("xi" + str(v)) if v < 8 else ("a" + str(v - 8)) for v in m)
        out.append(str(p[m]) + ("*" + mono if mono else ""))
    return "+".join(out)


# ================================================================= measurement
def measure():
    xi = ovars(0)
    av = ovars(8)
    XI = list(range(8))
    AV = list(range(8, 16))
    B_XI1 = basis([(XI, 1)])
    B_XI2 = basis([(XI, 2)])
    B_11 = basis([(XI, 1), (AV, 1)])
    B_21 = basis([(XI, 2), (AV, 1)])
    B_22 = basis([(XI, 2), (AV, 2)])
    nsq = o_nsq(xi)

    # ---------------- A: the octonion table
    rule = lines_from_rule()
    tab = build_table(LINES_FROZEN)
    meas("M8.A1", "fano_lines_from_rule", rule == LINES_FROZEN and tab is not None,
         [("rule", "e_i*e_(i+1)=e_(i+3) mod 7 on 1..7"), ("lines", " ".join("(" + ",".join(str(v) for v in ln) + ")" for ln in rule)),
          ("rule_equals_frozen_lines", 1 if rule == LINES_FROZEN else 0), ("each_imaginary_pair_once", 1 if tab is not None else 0)])
    prod = o_mul(tab, xi, av)
    res = zero_check([p_sub(o_nsq(prod), p_mul(nsq, o_nsq(av)))], B_22)
    meas("M8.A2", "composition_nsq(xi_a)=nsq(xi)nsq(a)", holds(res), cpairs(res) + [("vars", 16), ("bidegree", "2,2")])

    # ---------------- B: gamma(xi) on O + O
    g = gamma_poly(tab, xi, -1)
    lin_ok = all(len(m) == 1 and m[0] < 8 for e in flat(g) for m in e)
    nz = sum(1 for e in flat(g) if e)
    gam = [coeff_int(g, t) for t in range(8)]
    recon_ok = True
    for r in range(16):
        for c in range(16):
            want = {}
            for t in range(8):
                if gam[t][r][c]:
                    want[(t,)] = gam[t][r][c]
            if want != g[r][c]:
                recon_ok = False
    meas("M8.B0", "gamma_structure", lin_ok and recon_ok and nz == 128,
         [("size", "16x16"), ("block", "[[0,L_xi],[-L_conj(xi),0]]"), ("S_plus", "rows 0..7"), ("S_minus", "rows 8..15"),
          ("entries_nonzero", nz), ("entries_single_signed_variable", 1 if lin_ok else 0), ("gamma_equals_sum_xi_t_gamma_t", 1 if recon_ok else 0)])
    g2 = pm_mul(g, g)
    nz2 = sum(len(e) for e in flat(g2))
    e2 = pm_lin(g2, pm_diag(nsq, 16), 1, 1)
    res = zero_check(flat(e2), B_XI2)
    meas("M8.B1", "gamma(xi)^2+nsq(xi)I16=0", holds(res),
         cpairs(res) + [("entries", 256), ("monomials_per_entry", len(B_XI2)), ("nonzero_coeffs_of_gamma_squared", nz2)])
    lm = L_matrix(tab, xi)
    lc = L_matrix(tab, o_conj(xi))
    res1 = zero_check(flat(pm_lin(pm_mul(lm, lc), pm_diag(nsq, 8), 1, -1)), B_XI2)
    res2 = zero_check(flat(pm_lin(pm_mul(lc, lm), pm_diag(nsq, 8), 1, -1)), B_XI2)
    meas("M8.B2", "L_xi*L_conj(xi)=nsq(xi)I8", holds(res1), cpairs(res1))
    meas("M8.B3", "L_conj(xi)*L_xi=nsq(xi)I8", holds(res2), cpairs(res2))
    res = zero_check(flat(pm_lin(lc, pm_T(lm), 1, -1)), B_XI1)
    meas("M8.B4", "L_conj(xi)=L_xi^T", holds(res), cpairs(res) + [("consequence", "gamma(xi)^T=-gamma(xi)")])
    bad = 0
    for i in range(8):
        for j in range(8):
            s = im_lin(im_mul(gam[i], gam[j]), im_mul(gam[j], gam[i]), 1, 1)
            want = im_lin(im_ident(16), im_ident(16), -2 if i == j else 0, 0)
            if s != want:
                bad += 1
    meas("M8.B5", "clifford_relations_gamma_i_gamma_j+gamma_j_gamma_i=-2delta_ij", bad == 0,
         [("ordered_pairs", 64), ("failed", bad), ("sign", "Cl(R^8,-nsq)")])

    # ---------------- C: the algebra generated by gamma(e0..e7)
    def sperm(m):
        rep = []
        for r in range(16):
            nzc = [(c, m[r][c]) for c in range(16) if m[r][c]]
            if len(nzc) != 1 or abs(nzc[0][1]) != 1:
                return None
            rep.append(nzc[0])
        return tuple(rep)

    words = []
    for mask in range(256):
        m = im_ident(16)
        for t in range(8):
            if (mask >> t) & 1:
                m = im_mul(m, gam[t])
        words.append(m)
    reps = [sperm(m) for m in words]
    all_sp = all(r is not None for r in reps)
    gram_bad = 0
    for i in range(256):
        for j in range(i, 256):
            ip = sum(reps[i][r][1] * reps[j][r][1] for r in range(16) if reps[i][r][0] == reps[j][r][0])
            if ip != (16 if i == j else 0):
                gram_bad += 1
    lookup = {}
    for n, r in enumerate(reps):
        lookup[r] = n
    closed = 0
    for n in range(256):
        for t in range(8):
            r = sperm(im_mul(words[n], gam[t]))
            neg = tuple((c, -v) for (c, v) in r) if r is not None else None
            if r is not None and (r in lookup or neg in lookup):
                closed += 1
    dim_alg = 256 if (all_sp and gram_bad == 0 and closed == 2048) else -1
    meas("M8.C1", "generated_algebra_dimension", dim_alg == 256,
         [("ordered_words", 256), ("all_signed_permutations", 1 if all_sp else 0), ("gram_trace_form", "16*I256"),
          ("gram_failures", gram_bad), ("closure_products_word_times_generator", closed), ("closure_expected", 2048),
          ("dimension", dim_alg), ("ambient", "M_16(R) dim 256"), ("algebra_equals_M16R", 1 if dim_alg == 256 else 0)])
    eqs = []
    for t in range(8):
        gm = gam[t]
        for a in range(16):
            for b in range(16):
                e = {}
                for c in range(16):
                    if gm[c][b]:
                        e[a * 16 + c] = e.get(a * 16 + c, 0) + gm[c][b]
                    if gm[a][c]:
                        e[c * 16 + b] = e.get(c * 16 + b, 0) - gm[a][c]
                e = {k: v for k, v in e.items() if v}
                if e:
                    eqs.append(e)
    rk = rank_mod_p(eqs)
    ident_in = all(sum(v * (1 if (k % 17) == 0 else 0) for k, v in e.items()) == 0 for e in eqs)
    dim_comm = 256 - rk
    meas("M8.C2", "commutant_of_gamma_e0..e7", dim_comm == 1 and ident_in,
         [("unknowns", 256), ("equations", len(eqs)), ("rank_mod_p", rk), ("p", "2^61-1"), ("identity_in_commutant", 1 if ident_in else 0),
          ("dimension", dim_comm), ("argument", "rank over Q >= rank mod p and I16 lies in the kernel")])
    eqs = []
    evens = [im_mul(gam[i], gam[j]) for (i, j) in PAIRS8]
    for gm in evens:
        for a in range(16):
            for b in range(16):
                e = {}
                for c in range(16):
                    if gm[c][b]:
                        e[a * 16 + c] = e.get(a * 16 + c, 0) + gm[c][b]
                    if gm[a][c]:
                        e[c * 16 + b] = e.get(c * 16 + b, 0) - gm[a][c]
                e = {k: v for k, v in e.items() if v}
                if e:
                    eqs.append(e)
    rk_even = rank_mod_p(eqs)

    # ---------------- E: the volume element gamma0...gamma7
    om = im_ident(16)
    for t in range(8):
        om = im_mul(om, gam[t])
    eps = om[0][0]
    diag_form = all(om[r][c] == 0 for r in range(16) for c in range(16) if r != c) and \
        all(om[r][r] == eps for r in range(8)) and all(om[r][r] == -eps for r in range(8, 16))
    anti = all(im_is_zero(im_lin(im_mul(om, gam[t]), im_mul(gam[t], om), 1, 1)) for t in range(8))
    meas("M8.E1", "volume_element_gamma0..gamma7", diag_form and anti and im_mul(om, om) == im_ident(16),
         [("value", "diag(" + str(eps) + "*I8," + str(-eps) + "*I8)" if diag_form else "not block scalar"),
          ("square", "I16" if im_mul(om, om) == im_ident(16) else "other"), ("anticommutes_with_each_gamma", 1 if anti else 0),
          ("S_plus_eigenvalue", eps)])
    # the even subalgebra's commutant: R I + R Omega
    om_in = True
    for e in eqs:
        if sum(v * om[divmod(k, 16)[0]][divmod(k, 16)[1]] for k, v in e.items()) != 0:
            om_in = False
    meas("M8.E2", "commutant_of_even_part_gamma_i_gamma_j", 256 - rk_even == 2 and om_in,
         [("equations", len(eqs)), ("rank_mod_p", rk_even), ("dimension", 256 - rk_even), ("volume_element_in_it", 1 if om_in else 0),
          ("reading", "S_plus and S_minus are inequivalent real irreducibles of the even part")])

    # ---------------- F: phi, psi, orientation, the Cayley forms
    phi = {}
    for (i, j, k) in TRIPLES7:
        s, kk = tab[i][j]
        if kk == k:
            phi[(i, j, k)] = s
    anti_ok = True
    for i in range(1, 8):
        for j in range(1, 8):
            for k in range(1, 8):
                s, kk = tab[i][j]
                direct = s if kk == k else 0
                ss, key = sort_sign((i, j, k))
                via = ss * phi.get(key, 0) if ss else 0
                if direct != via:
                    anti_ok = False
    phi_text = "+".join(("" if v > 0 else "-") + "e" + "".join(str(x) for x in k) for k, v in sorted(phi.items()))
    meas("M8.F1", "phi_ijk=<e_i e_j,e_k>", anti_ok and len(phi) == 7,
         [("terms", len(phi)), ("phi", phi_text), ("totally_antisymmetric", 1 if anti_ok else 0)])

    def t7_matrix(sigma):
        def fn(f):
            return f_star(f_wedge(phi, f), U7, sigma)
        return op_matrix(fn, PAIRS7, PAIRS7)

    a7 = t7_matrix(1)
    i21 = im_ident(21)
    q_plus = im_is_zero(im_mul(im_lin(a7, i21, 1, -2), im_lin(a7, i21, 1, 1)))
    q_minus = im_is_zero(im_mul(im_lin(a7, i21, 1, 2), im_lin(a7, i21, 1, -1)))
    r_p1 = rank_int(im_lin(a7, i21, 1, 1), 21)
    r_m2 = rank_int(im_lin(a7, i21, 1, -2), 21)
    r_m1 = rank_int(im_lin(a7, i21, 1, -1), 21)
    r_p2 = rank_int(im_lin(a7, i21, 1, 2), 21)
    if q_plus and r_p1 == 7 and r_m2 == 14:
        sigma7 = 1
        std = "e1..e7 gives +2 on 7 and -1 on 14; kept"
    elif q_minus and r_m1 == 7 and r_p2 == 14:
        sigma7 = -1
        std = "e1..e7 gives -2 on 7 and +1 on 14; reversed, vol7=-e1234567"
    else:
        sigma7 = 0
        std = "neither split"
    meas("M8.F2", "orientation_7D_by_T7_split", sigma7 != 0,
         [("T7", "*7(phi wedge .) on Lambda2(R7)"), ("std_orientation", std), ("sigma7", sigma7),
          ("rank_T7+I_std", r_p1), ("rank_T7-2I_std", r_m2), ("rank_T7-I_std", r_m1), ("rank_T7+2I_std", r_p2)])
    a7s = t7_matrix(sigma7)
    psi = f_star(phi, U7, sigma7)
    psi_text = "+".join(("" if v > 0 else "-") + "e" + "".join(str(x) for x in k) for k, v in sorted(psi.items()))
    meas("M8.F3", "psi=*7phi", len(psi) == 7, [("psi", psi_text), ("vol7_sign", sigma7)])
    sigma8 = sigma7
    star_ok = True
    for tri in TRIPLES7:
        if f_star(f_wedge({(0,): 1}, {tri: 1}), U8, sigma8) != f_star({tri: 1}, U7, sigma7):
            star_ok = False
    meas("M8.F4", "star8(e0 wedge alpha)=star7(alpha)", star_ok, [("alphas", len(TRIPLES7)), ("vol8", "e0 wedge vol7")])
    e0phi = f_wedge({(0,): 1}, phi)
    PHI = {}
    T8 = {}
    i28 = im_ident(28)
    for s in (1, -1):
        PHI[s] = f_lin(e0phi, psi, 1, s)

        def fn(f, s=s):
            return f_star(f_wedge(PHI[s], f), U8, sigma8)
        T8[s] = op_matrix(fn, PAIRS8, PAIRS8)
        sd = f_star(PHI[s], U8, sigma8) == f_lin(PHI[s], {}, s, 0)
        qq = im_is_zero(im_mul(im_lin(T8[s], i28, 1, -3), im_lin(T8[s], i28, 1, 1)))
        rp = rank_int(im_lin(T8[s], i28, 1, 1), 28)
        rm3 = rank_int(im_lin(T8[s], i28, 1, -3), 28)
        ptxt = "+".join(("" if v > 0 else "-") + "e" + "".join(str(x) for x in k) for k, v in sorted(PHI[s].items()))
        meas("M8.F5" + ("p" if s == 1 else "m"), "cayley_form_Phi_s=e0^phi+s*psi_s=" + str(s), sd and qq and rp == 7 and rm3 == 21,
             [("terms", len(PHI[s])), ("Phi", ptxt), ("star8_Phi=s_Phi", 1 if sd else 0),
              ("(T-3I)(T+I)=0", 1 if qq else 0), ("rank_T+I", rp), ("rank_T-3I", rm3), ("eigen", "3 on 7, -1 on 21")])

    # ---------------- G: Psi, Psi', kernels, spin(7), g2
    gg = [im_mul(gam[i], gam[j]) for (i, j) in PAIRS8]
    psi_m = im_zero(7, 28)
    psi_p = im_zero(7, 28)
    stray = 0
    for p in range(28):
        col = [gg[p][r][0] for r in range(16)]
        colp = [gg[p][r][8] for r in range(16)]
        if col[0] or any(col[8:]):
            stray += 1
        if colp[8] or any(colp[:8]):
            stray += 1
        for k in range(1, 8):
            psi_m[k - 1][p] = col[k]
            psi_p[k - 1][p] = colp[8 + k]
    pp = im_mul(psi_m, im_T(psi_m))
    ppp = im_mul(psi_p, im_T(psi_p))
    four7 = im_lin(im_ident(7), im_ident(7), 4, 0)
    meas("M8.G1", "Psi(omega)=c(omega)eta_eta=1_in_S+", stray == 0 and pp == four7,
         [("shape", "7x28 into Im O of S+"), ("real_and_S-_components", stray), ("Psi_PsiT", "4*I7" if pp == four7 else "other"),
          ("nonzero_entries", sum(1 for r in psi_m for x in r if x))])
    meas("M8.G2", "Psi'(omega)=c(omega)eta'_eta'=1_in_S-", ppp == four7,
         [("Psi'_Psi'T", "4*I7" if ppp == four7 else "other")])
    e0j = all(psi_m[j - 1][PAIRS8.index((0, j))] == 1 and sum(abs(psi_m[k][PAIRS8.index((0, j))]) for k in range(7)) == 1 for j in range(1, 8))
    psi7 = im_zero(7, 28)
    for p, (i, j) in enumerate(PAIRS8):
        if i == 0:
            continue
        for k in range(1, 8):
            ss, key = sort_sign((i, j, k))
            psi7[k - 1][p] = ss * phi.get(key, 0) if ss else 0
    so7cols = [p for p, (i, j) in enumerate(PAIRS8) if i > 0]
    restr = all(psi_m[k][p] == psi7[k][p] for k in range(7) for p in so7cols)
    meas("M8.G3", "Psi_on_e0j_and_on_so7", e0j and restr,
         [("Psi(e0j)=e_j", 1 if e0j else 0), ("Psi_on_so7=contraction_with_phi", 1 if restr else 0)])
    ker_psi = kernel_int(psi_m, 28)
    rk_psi = rank_int(psi_m, 28)
    tps = {}
    for s in (1, -1):
        tps[s] = (im_lin(T8[s], i28, 1, 1) == im_mul(im_T(psi_m), psi_m), im_lin(T8[s], i28, 1, 1) == im_mul(im_T(psi_p), psi_p))
    s_match = 1 if tps[1][0] else (-1 if tps[-1][0] else 0)
    sp_match = 1 if tps[1][1] else (-1 if tps[-1][1] else 0)
    meas("M8.G4", "T_s+I=Psi^T_Psi", s_match != 0,
         [("T_+1+I=PsiT_Psi", 1 if tps[1][0] else 0), ("T_-1+I=PsiT_Psi", 1 if tps[-1][0] else 0),
          ("T_+1+I=Psi'T_Psi'", 1 if tps[1][1] else 0), ("T_-1+I=Psi'T_Psi'", 1 if tps[-1][1] else 0),
          ("Psi_matches_s", s_match), ("Psi'_matches_s", sp_match)])
    stab = {}
    for s in (1, -1):
        cols = [f_vec(f_deriv(R_of(unit(28, p)), PHI[s], 8), QUADS8) for p in range(28)]
        stab[s] = [[cols[p][q] for p in range(28)] for q in range(len(QUADS8))]
    dims = {}
    for s in (1, -1):
        dims[s] = (28 - rank_int(stab[s], 28), 28 - rank_int(stab[s] + psi_m, 28), 28 - rank_int(stab[s] + psi_p, 28))
    dm = dims.get(s_match, (0, 0, 0))
    meas("M8.G5", "ker_Psi_dimension_and_stabilizer_of_Phi", 28 - rk_psi == 21 and len(ker_psi) == 21 and
         dm[0] == 21 and dm[1] == 21,
         [("rank_Psi", rk_psi), ("dim_ker_Psi", 28 - rk_psi), ("dim_stab_Phi_+1", dims[1][0]), ("dim_stab_Phi_-1", dims[-1][0]),
          ("dim_kerPsi_cap_stab_Phi_+1", dims[1][1]), ("dim_kerPsi_cap_stab_Phi_-1", dims[-1][1]),
          ("dim_kerPsi'_cap_stab_Phi_+1", dims[1][2]), ("dim_kerPsi'_cap_stab_Phi_-1", dims[-1][2]),
          ("ker_Psi_equals_stab_Phi_s", s_match)])
    closed_bad = 0
    pairs_done = 0
    rms = [R_of(v) for v in ker_psi]
    for a in range(21):
        for b in range(a + 1, 21):
            cm = im_lin(im_mul(rms[a], rms[b]), im_mul(rms[b], rms[a]), 1, -1)
            pairs_done += 1
            if not is_skew(cm) or any(im_vec(psi_m, w_of(cm))):
                closed_bad += 1
    br_bad = 0
    for p in range(28):
        rp_ = R_of(unit(28, p))
        for q in range(28):
            rq = R_of(unit(28, q))
            lhs = f_vec(f_deriv(rp_, {PAIRS8[q]: 1}, 8), PAIRS8)
            rhs = w_of(im_lin(im_mul(rp_, rq), im_mul(rq, rp_), 1, -1))
            if lhs != rhs:
                br_bad += 1
    meas("M8.G6", "ker_Psi_closed_under_commutator", closed_bad == 0 and br_bad == 0,
         [("basis", 21), ("pairs", pairs_done), ("failed", closed_bad), ("bracket", "[R_a,R_b] of C0.5 rotations"),
          ("derivation_on_2forms_equals_commutator", 1 if br_bad == 0 else 0), ("pairs_checked_for_that", 784)])
    rk_both = rank_int(psi_m + psi_p, 28)
    ker_both = kernel_int(psi_m + psi_p, 28)
    so7_rows = [unit(28, PAIRS8.index((0, j))) for j in range(1, 8)]
    cols3 = [f_vec(f_deriv(R_of(unit(28, p)), phi, 8), TRIPLES8) for p in range(28)]
    stab_phi = [[cols3[p][q] for p in range(28)] for q in range(len(TRIPLES8))]
    d_stabphi = 28 - rank_int(stab_phi, 28)
    d_psi7 = 28 - rank_int(psi7 + so7_rows, 28)
    d_all = 28 - rank_int(psi_m + psi_p + stab_phi + psi7 + so7_rows, 28)
    in_so7 = all(v[PAIRS8.index((0, j))] == 0 for v in ker_both for j in range(1, 8))
    rms2 = [R_of(v) for v in ker_both]
    g2_bad = 0
    for a in range(len(rms2)):
        for b in range(a + 1, len(rms2)):
            cm = im_lin(im_mul(rms2[a], rms2[b]), im_mul(rms2[b], rms2[a]), 1, -1)
            w = w_of(cm)
            if any(im_vec(psi_m, w)) or any(im_vec(psi_p, w)):
                g2_bad += 1
    meas("M8.G7", "ker_Psi_cap_ker_Psi'=g2", 28 - rk_both == 14 and d_stabphi == 14 and d_psi7 == 14 and d_all == 14 and in_so7 and g2_bad == 0,
         [("dim_ker_Psi_cap_ker_Psi'", 28 - rk_both), ("dim_stab_so8_phi", d_stabphi), ("dim_ker_Psi7_in_so7", d_psi7),
          ("dim_of_all_three_intersected", d_all), ("lies_in_so7", 1 if in_so7 else 0), ("commutator_pairs", 91), ("closure_failed", g2_bad)])

    # ---------------- H: spin action and equivariance of gamma
    def cliff(w):
        m = im_zero(16, 16)
        for p in range(28):
            if w[p]:
                m = im_lin(m, gg[p], 1, w[p])
        return m

    sum_bad = 0
    for p in range(28):
        rm = R_of(unit(28, p))
        acc = im_zero(16, 16)
        for a in range(8):
            for b in range(8):
                if rm[a][b]:
                    acc = im_lin(acc, im_mul(gam[a], gam[b]), 1, rm[a][b])
        if acc != im_lin(gg[p], gg[p], -2, 0):
            sum_bad += 1
    meas("M8.H1", "matrix_entry_sum_X_ab_gamma_a_gamma_b=-2c(omega)", sum_bad == 0,
         [("forms", 28), ("failed", sum_bad), ("X", "R_omega of C0.5, X_ab its row a column b entry"),
          ("so", "(1/4)sum X_ab gamma_a gamma_b=-(1/2)c(omega); (1/4)sum omega_ab gamma_a gamma_b=+(1/2)c(omega) with omega_ab the form coefficients")])

    def equiv_check(w):
        cm = cliff(w)
        rm = R_of(w)
        lhs = pm_lin(pm_mul(pm_from_int(cm), g), pm_mul(g, pm_from_int(cm)), 1, -1)
        rhs = gamma_poly(tab, rot_vec(rm, xi), -1)
        return zero_check(flat(pm_lin(lhs, rhs, 1, -2)), B_XI1)

    tot = [0, 0, 0]
    for p in range(28):
        r = equiv_check(unit(28, p))
        tot = [tot[0] + r[0], tot[1] + r[1], tot[2] + r[2]]
    meas("M8.H2", "C0.5_[c(omega),gamma(xi)]=2gamma(R_omega xi)_all_28_basis_forms", tot[1] == 0 and tot[2] == 0, cpairs(tot))
    tot = [0, 0, 0]
    for v in ker_psi:
        r = equiv_check(v)
        tot = [tot[0] + r[0], tot[1] + r[1], tot[2] + r[2]]
    meas("M8.H3", "spin7_equivariance_[rho(X),gamma(xi)]=gamma(X xi)_integer_form_[c(omega_X),gamma(xi)]=2gamma(X xi)",
         tot[1] == 0 and tot[2] == 0, cpairs(tot) + [("basis", len(ker_psi)), ("rho", "rho(X)=(1/2)c(omega_X)=(1/4)sum omega_ab gamma_a gamma_b")])
    ann_bad = 0
    eta_p_span = []
    for v in ker_psi:
        cm = cliff(v)
        if any(cm[r][0] for r in range(16)):
            ann_bad += 1
        eta_p_span.append([cm[r][8] for r in range(16)])
    meas("M8.H4", "spin7_annihilates_eta", ann_bad == 0,
         [("basis", len(ker_psi)), ("failed", ann_bad), ("rank_of_spin7_on_eta'", rank_int(eta_p_span, 16))])
    int_bad = 0
    for v in ker_psi:
        cm = cliff(v)
        rv = R_of(v)
        for q in range(28):
            rq = R_of(unit(28, q))
            u = w_of(im_lin(im_mul(rv, rq), im_mul(rq, rv), 1, -1))
            emb = [0] + [psi_m[k][q] for k in range(7)] + [0] * 8
            lhs = im_vec(cm, emb)
            rhs = [0] + [2 * x for x in im_vec(psi_m, u)] + [0] * 8
            if lhs != rhs:
                int_bad += 1
    meas("M8.H5", "Psi_intertwines_c(omega_X)Psi(omega)=2Psi(X.omega)", int_bad == 0,
         [("pairs", 21 * 28), ("failed", int_bad)])

    # ---------------- I: the Spin(7) symbol sequence 1, 8, 7
    def wedge_cols(xv):
        cols = []
        for j in range(8):
            w = [{} for _ in range(28)]
            for i in range(8):
                if i < j:
                    w[PAIRS8.index((i, j))] = p_add(w[PAIRS8.index((i, j))], xv[i])
                elif i > j:
                    w[PAIRS8.index((j, i))] = p_sub(w[PAIRS8.index((j, i))], xv[i])
            cols.append(w)
        return cols

    def P_poly(xv, pm):
        cols = wedge_cols(xv)
        out = pm_zero(7, 8)
        for j in range(8):
            for k in range(7):
                acc = {}
                for p in range(28):
                    if pm[k][p] and cols[j][p]:
                        p_add_into(acc, cols[j][p], pm[k][p])
                out[k][j] = acc
        return out

    P = P_poly(xi, psi_m)
    s1 = [[dict(xi[t])] for t in range(8)]
    r_comp = zero_check(flat(pm_mul(P, s1)), B_XI2)
    r_ppt = zero_check(flat(pm_lin(pm_mul(P, pm_T(P)), pm_diag(nsq, 7), 1, -1)), B_XI2)
    xxt = [[p_mul(xi[i], xi[j]) for j in range(8)] for i in range(8)]
    r_ptp = zero_check(flat(pm_lin(pm_lin(pm_mul(pm_T(P), P), pm_diag(nsq, 8), 1, -1), xxt, 1, 1)), B_XI2)
    r_s1 = zero_check(flat(pm_lin(pm_mul(pm_T(s1), s1), [[dict(nsq)]], 1, -1)), B_XI2)
    seq_ok = holds(r_comp) and holds(r_ppt) and holds(r_ptp) and holds(r_s1)
    meas("M8.I1", "spin7_symbol_sequence", seq_ok,
         [("sequence", "Lambda0 -> Lambda1 -> Lambda2_7"), ("dims", "1,8,7"), ("euler", 1 - 8 + 7), ("ranks", "1,7"),
          ("sigma1", "xi wedge ."), ("sigma2", "P(xi)=Psi(xi wedge .)")])
    meas("M8.I2", "P(xi)xi=0_composition", holds(r_comp), cpairs(r_comp))
    meas("M8.I3", "P(xi)P(xi)^T=nsq(xi)I7_rank7", holds(r_ppt), cpairs(r_ppt))
    meas("M8.I4", "P(xi)^T_P(xi)=nsq(xi)I8-xi_xi^T_exact_at_Lambda1", holds(r_ptp), cpairs(r_ptp))
    meas("M8.I5", "xi^T_xi=nsq(xi)_rank1", holds(r_s1), cpairs(r_s1))
    tpi = im_lin(T8[s_match], i28, 1, 1) if s_match else im_lin(T8[1], i28, 1, 1)
    wc = wedge_cols(xi)
    ent = []
    for j in range(8):
        for q in range(28):
            lhs = {}
            for p in range(28):
                if tpi[q][p] and wc[j][p]:
                    p_add_into(lhs, wc[j][p], tpi[q][p])
            rhs = {}
            for k in range(7):
                if psi_m[k][q] and P[k][j]:
                    p_add_into(rhs, P[k][j], psi_m[k][q])
            ent.append(p_sub(lhs, rhs))
    r_tp = zero_check(ent, B_XI1)
    meas("M8.I6", "(T+I)(xi wedge a)=Psi^T_P(xi)a", holds(r_tp), cpairs(r_tp) + [("s", s_match)])

    # ---------------- J: the rolled-up symbol N(xi)
    def N_poly(xv, pm, sign0):
        pmat = P_poly(xv, pm)
        return [[p_scale(xv[j], sign0) for j in range(8)]] + pmat

    N = N_poly(xi, psi_m, -1)
    r_ntn = zero_check(flat(pm_lin(pm_mul(pm_T(N), N), pm_diag(nsq, 8), 1, -1)), B_XI2)
    r_nnt = zero_check(flat(pm_lin(pm_mul(N, pm_T(N)), pm_diag(nsq, 8), 1, -1)), B_XI2)
    meas("M8.J1", "N(xi)=(-iota_xi,P(xi))_N^T_N=nsq(xi)I8", holds(r_ntn), cpairs(r_ntn) + [("rows", "0: -xi^T, 1..7: P(xi)")])
    meas("M8.J2", "N(xi)N(xi)^T=nsq(xi)I8", holds(r_nnt), cpairs(r_nnt))
    cx = o_conj(xi)
    ca = o_conj(av)
    cands = [("xi*a", o_mul(tab, xi, av)), ("a*xi", o_mul(tab, av, xi)), ("conj(xi)*a", o_mul(tab, cx, av)),
             ("a*conj(xi)", o_mul(tab, av, cx)), ("xi*conj(a)", o_mul(tab, xi, ca)), ("conj(a)*xi", o_mul(tab, ca, xi)),
             ("conj(xi)*conj(a)", o_mul(tab, cx, ca)), ("conj(a)*conj(xi)", o_mul(tab, ca, cx))]
    acol = [[dict(av[t])] for t in range(8)]

    def expr_table(nm):
        na = [r[0] for r in pm_mul(nm, acol)]
        found = []
        for name, ex in cands:
            for sg in (1, -1):
                res = zero_check([p_sub(na[t], p_scale(ex[t], sg)) for t in range(8)], B_11)
                if holds(res):
                    found.append(("" if sg == 1 else "-") + name)
        return na, found

    na, found = expr_table(N)
    meas("M8.J3", "N(xi)a_as_octonion_expression", len(found) == 1,
         [("candidates_tested", 16), ("holding", ",".join(found) if found else "none"), ("coeffs_per_candidate", 8 * len(B_11)),
          ("coordinates", "component 0 = Lambda0, components 1..7 = Psi coordinates of Lambda2_7 on e1..e7")])
    r2 = zero_check([p_sub(na[t], o_mul(tab, cx, av)[t]) for t in range(8)], B_11)
    r2m = zero_check([p_add(na[t], o_mul(tab, cx, av)[t]) for t in range(8)], B_11)
    meas("M8.J4", "addendum_R2_form_N(xi)a=conj(xi)a", holds(r2),
         cpairs(r2) + [("minus_form_N(xi)a=-conj(xi)a_holds", 1 if holds(r2m) else 0),
                       ("N_used", "Psi from eta in S+, row 0 = -iota_xi (C0.2)")])
    variants = []
    for pname, pmat in (("Psi", psi_m), ("Psi'", psi_p)):
        for sign0 in (-1, 1):
            nv = N_poly(xi, pmat, sign0)
            _, fv = expr_table(nv)
            variants.append(pname + "_row0=" + ("-" if sign0 < 0 else "+") + "iota:" + (",".join(fv) if fv else "none"))
    meas("M8.J5", "N_variants_octonion_expression", True, [("variants", " ; ".join(variants))])

    # ---------------- K: R1, the raw projection identity, c = 1/4
    w = [{} for _ in range(28)]
    for p, (i, j) in enumerate(PAIRS8):
        w[p] = p_sub(p_mul(xi[i], av[j]), p_mul(xi[j], av[i]))
    tw = []
    for q in range(28):
        acc = {}
        for p in range(28):
            if tpi[q][p]:
                p_add_into(acc, w[p], tpi[q][p])
        tw.append(acc)
    nsq_w = {}
    for p in range(28):
        p_add_into(nsq_w, p_mul(w[p], w[p]))
    w_tw = {}
    for p in range(28):
        p_add_into(w_tw, p_mul(w[p], tw[p]))
    nsq_tw = {}
    for p in range(28):
        p_add_into(nsq_tw, p_mul(tw[p], tw[p]))
    psw = []
    for k in range(7):
        acc = {}
        for p in range(28):
            if psi_m[k][p]:
                p_add_into(acc, w[p], psi_m[k][p])
        psw.append(acc)
    nsq_psw = {}
    for k in range(7):
        p_add_into(nsq_psw, p_mul(psw[k], psw[k]))
    dot = o_dot(xi, av)
    lagr = p_sub(p_mul(nsq, o_nsq(av)), p_mul(dot, dot))
    r_l = zero_check([p_sub(nsq_w, lagr)], B_22)
    r_a = zero_check([p_sub(w_tw, nsq_w)], B_22)
    r_b = zero_check([p_sub(nsq_tw, p_scale(nsq_w, 4))], B_22)
    r_c = zero_check([p_sub(nsq_psw, nsq_w)], B_22)
    meas("M8.K1", "C0.9_nsq(xi^a)=nsq(xi)nsq(a)-<xi,a>^2", holds(r_l), cpairs(r_l))
    meas("M8.K2", "R1_c=1/4_<xi^a,(T+I)(xi^a)>=nsq(xi^a)", holds(r_a),
         cpairs(r_a) + [("reading", "nsq(pi7(xi^a))=(1/4)nsq(xi^a), pi7=(T+I)(1/4)"), ("c", "1/4")])
    meas("M8.K3", "R1_nsq((T+I)(xi^a))=4nsq(xi^a)", holds(r_b), cpairs(r_b))
    meas("M8.K4", "R1_nsq(Psi(xi^a))=nsq(xi^a)", holds(r_c), cpairs(r_c))
    r_h = zero_check([p_sub(w_tw, p_scale(nsq_w, 2))], B_22)
    r_t = zero_check([p_sub(p_scale(w_tw, 3), p_scale(nsq_w, 4))], B_22)
    ctrl("K8.01", "R1_with_c=1/2_<w,(T+I)w>=2nsq(w)", True, not holds(r_h), cpairs(r_h))
    ctrl("K8.02", "R1_with_c=1/3_3<w,(T+I)w>=4nsq(w)", True, not holds(r_t), cpairs(r_t))
    om_v = [p_var(p) for p in range(28)]
    g_tw = {}
    g_n = {}
    for q in range(28):
        acc = {}
        for p in range(28):
            if tpi[q][p]:
                p_add_into(acc, om_v[p], tpi[q][p])
        p_add_into(g_tw, p_mul(om_v[q], acc))
        p_add_into(g_n, p_mul(om_v[q], om_v[q]))
    r_g = zero_check([p_sub(g_tw, g_n)], basis([(list(range(28)), 2)]))
    ctrl("K8.03", "R1_on_a_general_2form_<omega,(T+I)omega>=nsq(omega)", True, not holds(r_g), cpairs(r_g) + [("vars", 28)])

    # ---------------- L: R3 and alternativity
    def alt(lhs_vec, rhs_vec, mon):
        return zero_check([p_sub(lhs_vec[t], rhs_vec[t]) for t in range(8)], mon)

    nx = o_scalar(nsq, av)
    r3 = alt(o_mul(tab, cx, o_mul(tab, xi, av)), nx, B_21)
    r3b = alt(o_mul(tab, xi, o_mul(tab, cx, av)), nx, B_21)
    r3c = alt(o_mul(tab, o_mul(tab, av, xi), cx), nx, B_21)
    r3d = alt(o_mul(tab, o_mul(tab, av, cx), xi), nx, B_21)
    meas("M8.L1", "R3_conj(xi)(xi x)=nsq(xi)x", holds(r3), cpairs(r3) + [("x", "variables 8..15")])
    meas("M8.L2", "xi(conj(xi)x)=nsq(xi)x", holds(r3b), cpairs(r3b))
    meas("M8.L3", "(x xi)conj(xi)=nsq(xi)x_and_(x conj(xi))xi=nsq(xi)x", holds(r3c) and holds(r3d),
         [("first_failed", r3c[1] + r3c[2]), ("second_failed", r3d[1] + r3d[2]), ("coeffs_checked", r3c[0] + r3d[0])])
    rl = alt(o_mul(tab, xi, o_mul(tab, xi, av)), o_mul(tab, o_mul(tab, xi, xi), av), B_21)
    meas("M8.L4", "left_alternative_xi(xi x)=(xi xi)x", holds(rl), cpairs(rl))
    r_bad = alt(o_mul(tab, xi, o_mul(tab, xi, av)), o_scalar(p_scale(nsq, -1), av), B_21)
    ctrl("K8.04", "R3_with_xi_for_conj(xi)_xi(xi x)=-nsq(xi)x_on_all_of_O", True, not holds(r_bad), cpairs(r_bad))
    im_xi = [{}] + [p_var(t) for t in range(1, 8)]
    nim = o_nsq(im_xi)
    r_im = alt(o_mul(tab, im_xi, o_mul(tab, im_xi, av)), o_scalar(p_scale(nim, -1), av), basis([(list(range(1, 8)), 2), (AV, 1)]))
    meas("M8.L5", "xi(xi x)=-nsq(xi)x_for_imaginary_xi", holds(r_im), cpairs(r_im) + [("xi0", 0)])
    e = [obasis(t) for t in range(8)]
    assoc = o_lin(o_mul(tab, o_mul(tab, e[1], e[2]), e[3]), o_mul(tab, e[1], o_mul(tab, e[2], e[3])), 1, -1)
    assoc_int = [p_eval(c, {}) if c else 0 for c in assoc]
    nz_assoc = 0
    on_line = 0
    for (i, j, k) in TRIPLES7:
        a_ = o_lin(o_mul(tab, o_mul(tab, e[i], e[j]), e[k]), o_mul(tab, e[i], o_mul(tab, e[j], e[k])), 1, -1)
        if any(a_):
            nz_assoc += 1
        else:
            on_line += 1
    ctrl("K8.05", "associativity_on_(e1,e2,e3)", True, any(assoc_int),
         [("(e1e2)e3-e1(e2e3)", o_text(assoc_int)), ("triples_with_nonzero_associator", nz_assoc), ("triples_associative", on_line)])

    # ---------------- M: equivariance of the symbols under spin(7)
    pe_bad = [0, 0, 0]
    ne_bad = [0, 0, 0]
    block_bad = 0
    for v in ker_psi:
        cm = cliff(v)
        rv = R_of(v)
        if any(cm[r][c] for r in range(16) for c in range(16) if (r < 8) != (c < 8)) or any(cm[0][c] for c in range(16)) or any(cm[r][0] for r in range(16)):
            block_bad += 1
        c7 = [[cm[r][c] for c in range(1, 8)] for r in range(1, 8)]
        prx = P_poly(rot_vec(rv, xi), psi_m)
        lhs = pm_lin(prx, pm_mul(P, pm_from_int(rv)), 2, 2)
        res = zero_check(flat(pm_lin(lhs, pm_mul(pm_from_int(c7), P), 1, -1)), B_XI1)
        pe_bad = [pe_bad[0] + res[0], pe_bad[1] + res[1], pe_bad[2] + res[2]]
        c8 = [[0] * 8] + [[0] + r for r in c7]
        nrx = N_poly(rot_vec(rv, xi), psi_m, -1)
        lhs = pm_lin(nrx, pm_mul(N, pm_from_int(rv)), 2, 2)
        res = zero_check(flat(pm_lin(lhs, pm_mul(pm_from_int(c8), N), 1, -1)), B_XI1)
        ne_bad = [ne_bad[0] + res[0], ne_bad[1] + res[1], ne_bad[2] + res[2]]
    meas("M8.M1", "spin7_symbol_equivariance_2P(Xxi)+2P(xi)X=c7(X)P(xi)", holds(pe_bad) and block_bad == 0,
         cpairs(pe_bad) + [("basis", len(ker_psi)), ("c(omega_X)_preserves_ImO_in_S+", 1 if block_bad == 0 else 0)])
    meas("M8.M2", "spin7_rolled_symbol_equivariance_2N(Xxi)+2N(xi)X=diag(0,c7(X))N(xi)", holds(ne_bad), cpairs(ne_bad))

    # ---------------- controls
    gw = gamma_poly(tab, xi, 1)
    gw2 = pm_mul(gw, gw)
    rw = zero_check(flat(pm_lin(gw2, pm_diag(nsq, 16), 1, 1)), B_XI2)
    rwp = zero_check(flat(pm_lin(gw2, pm_diag(nsq, 16), 1, -1)), B_XI2)
    ctrl("K8.06", "wrong_sign_gamma=[[0,L_xi],[+L_conj(xi),0]]_squares_to_-nsq", True, not holds(rw),
         cpairs(rw) + [("squares_to_+nsq_instead", 1 if holds(rwp) else 0), ("symmetric", 1 if holds(zero_check(flat(pm_lin(gw, pm_T(gw), 1, -1)), B_XI1)) else 0)])
    flipped = [(2, 1, 4)] + LINES_FROZEN[1:]
    tf = build_table(flipped)
    lmf = L_matrix(tf, xi)
    lcf = L_matrix(tf, o_conj(xi))
    rf = zero_check(flat(pm_lin(pm_mul(lmf, lcf), pm_diag(nsq, 8), 1, -1)), B_XI2)
    rfc = zero_check([p_sub(o_nsq(o_mul(tf, xi, av)), p_mul(nsq, o_nsq(av)))], B_22)
    ctrl("K8.07", "one_fano_line_reversed_(2,1,4)_L_xi_L_conj(xi)=nsq_I8", True, not holds(rf),
         cpairs(rf) + [("composition_failed_coeffs", rfc[1] + rfc[2])])
    pt = {t: t + 1 for t in range(8)}
    samp_ok = all(p_eval(e2[r][c], pt) == 0 for r in range(16) for c in range(16))
    g2_at = [[p_eval(g2[r][c], pt) for c in range(16)] for r in range(16)]
    ctrl("K8.08", "sampled_covector_xi=(1..8)_gamma^2=-204I_(control_only)", False, not samp_ok,
         [("nsq", sum((t + 1) * (t + 1) for t in range(8))), ("gamma^2_diag", g2_at[0][0]), ("all_entries_match", 1 if samp_ok else 0)])
    fake = pm_lin(g2, pm_diag(p_add(nsq, p_sub(p_mul(xi[0], xi[0]), xi[0])), 16), 1, 1)
    fake_at_sample = all(p_eval(fake[r][c], pt) == 0 for r in range(16) for c in range(16))
    rfk = zero_check(flat(fake), B_XI2)
    ctrl("K8.09", "planted_identity_gamma^2=-(nsq+xi0^2-xi0)I_passes_the_sample_fails_as_polynomial", True, not holds(rfk),
         cpairs(rfk) + [("passes_at_sample", 1 if fake_at_sample else 0)])
    w01 = unit(28, 0)
    c01 = cliff(w01)
    eta_hit = any(c01[r][0] for r in range(16))
    rv = R_of(w01)
    fails = 0
    for q in range(28):
        rq = R_of(unit(28, q))
        u = w_of(im_lin(im_mul(rv, rq), im_mul(rq, rv), 1, -1))
        emb = [0] + [psi_m[k][q] for k in range(7)] + [0] * 8
        lhs = im_vec(c01, emb)
        rhs = [0] + [2 * x for x in im_vec(psi_m, u)] + [0] * 8
        if lhs != rhs:
            fails += 1
    ctrl("K8.10", "X=e01_outside_spin7_annihilates_eta_and_intertwines_Psi", True, eta_hit and fails > 0,
         [("c(e01)eta_nonzero", 1 if eta_hit else 0), ("intertwining_failures_of_28", fails)])
    other = -s_match if s_match else -1
    ctrl("K8.11", "T_(-s)+I=Psi^T_Psi_for_the_other_Cayley_form", True,
         im_lin(T8[other], i28, 1, 1) != im_mul(im_T(psi_m), psi_m), [("s_tested", other)])
    wrong_n = 15
    ctrl("K8.12", "the_other_15_octonion_candidates_for_N(xi)a", True, len(found) == 1, [("failing_candidates", 16 - len(found)), ("expected_failing", wrong_n)])
    return {"tab": tab, "xi": xi, "av": av, "nsq": nsq, "gam": gam, "lm": lm, "lc": lc, "g": g, "om": om, "phi": phi,
            "psi_m": psi_m, "psi_p": psi_p, "psi7": psi7, "T8": T8, "ker_psi": ker_psi, "s_match": s_match, "N": N, "P": P,
            "cliff": cliff, "P_poly": P_poly, "N_poly": N_poly, "wedge_cols": wedge_cols, "tpi": tpi, "PHI": PHI, "sigma7": sigma7}


# ================================================================= after the measured block: supplement and compare
PRE_READ_BLOCK_SHA = "bbd58c4256bc8e88de15a351e3461b2190f209086eeabe5a58e770489735d8f7"
PRE_READ_FILE_SHA = "e737ec713e072c2c12f0393fbc41a124e09554bd61a8a5579ebcef12fa2563cd"
PRE_READ_SCRIPT_SHA = "c50948278e8ba0c5a23b322f965eccb45bcf7f48f50db2167fba09f933e80021"
PRED_SHA = "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e"
ADD_SHA = "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0"
SUPP = {}
VERDICTS = {"AGREE": [], "DISAGREE": [], "NOT_MEASURED": [], "NOT_COMPARED": [], "CONVENTION_FOLLOWED": []}


def supp_row(sid, name, ok, pairs):
    SUPP[sid] = dict(pairs)
    SUPP[sid]["ok"] = ok
    row("SUPP", [("id", sid), ("name", name)] + pairs + [("holds", 1 if ok else 0), ("timing", "POST_READ"), ("claim", "MEASURED")])


def intertwiner_dim(As, Bs, n):
    rows = []
    for A, B in zip(As, Bs):
        for a in range(n):
            for b in range(n):
                r = [0] * (n * n)
                for c in range(n):
                    r[a * n + c] += A[c][b]
                    r[c * n + b] -= B[a][c]
                if any(r):
                    rows.append(r)
    return n * n - rank_int(rows, n * n)


def algebra_of(gens):
    def sperm(m):
        rep = []
        for r in range(16):
            nzc = [(c, m[r][c]) for c in range(16) if m[r][c]]
            if len(nzc) != 1 or abs(nzc[0][1]) != 1:
                return None
            rep.append(nzc[0])
        return tuple(rep)
    words = []
    for mask in range(256):
        m = im_ident(16)
        for t in range(8):
            if (mask >> t) & 1:
                m = im_mul(m, gens[t])
        words.append(m)
    reps = [sperm(m) for m in words]
    if any(r is None for r in reps):
        return -1
    for i in range(256):
        for j in range(i, 256):
            ip = sum(reps[i][r][1] * reps[j][r][1] for r in range(16) if reps[i][r][0] == reps[j][r][0])
            if ip != (16 if i == j else 0):
                return -1
    lookup = set(reps)
    for n in range(256):
        for t in range(8):
            r = sperm(im_mul(words[n], gens[t]))
            if r is None or (r not in lookup and tuple((c, -v) for (c, v) in r) not in lookup):
                return -1
    return 256


def supplement(ctx):
    tab = ctx["tab"]
    xi = ctx["xi"]
    av = ctx["av"]
    nsq = ctx["nsq"]
    gam = ctx["gam"]
    lm = ctx["lm"]
    lc = ctx["lc"]
    psi_m = ctx["psi_m"]
    psi_p = ctx["psi_p"]
    T8 = ctx["T8"]
    phi = ctx["phi"]
    ker_psi = ctx["ker_psi"]
    cliff = ctx["cliff"]
    N = ctx["N"]
    XI = list(range(8))
    AV = list(range(8, 16))
    B_XI1 = basis([(XI, 1)])
    B_XI2 = basis([(XI, 2)])
    i28 = im_ident(28)
    Lint = [[[lm[r][c].get((t,), 0) for c in range(8)] for r in range(8)] for t in range(8)]
    Lcint = [[[lc[r][c].get((t,), 0) for c in range(8)] for r in range(8)] for t in range(8)]
    Rint = []
    for t in range(8):
        rm = im_zero(8, 8)
        for j in range(8):
            s, k = tab[j][t]
            rm[k][j] = s
        Rint.append(rm)
    kint = im_ident(8)
    kint[0][0] = -1

    # S8.01 the 7D volume sign that 8D.P4 reuses, and the blocks of gamma0...gamma7
    p7 = im_ident(8)
    q7 = im_ident(8)
    for t in range(1, 8):
        p7 = im_mul(p7, Lint[t])
        q7 = im_mul(q7, Rint[t])
    eps7 = p7[0][0]
    epsr = q7[0][0]
    scal = p7 == im_lin(im_ident(8), im_ident(8), eps7, 0)
    scalr = q7 == im_lin(im_ident(8), im_ident(8), epsr, 0)
    om = ctx["om"]
    b11 = [[om[r][c] for c in range(8)] for r in range(8)] == p7
    b22 = [[om[8 + r][8 + c] for c in range(8)] for r in range(8)] == im_lin(p7, p7, -1, 0)
    supp_row("S8.01", "L_e1...L_e7=eps7_I8_and_Omega_blocks", scal and b11 and b22,
             [("eps7", eps7 if scal else "not scalar"), ("R_e1...R_e7_eps", epsr if scalr else "not scalar"),
              ("Omega_block11=L_e1...L_e7", 1 if b11 else 0), ("Omega_block22=-L_e1...L_e7", 1 if b22 else 0)])
    # S8.02 is there a fixed U with U L_xi = L_conj(xi) U for all xi, i.e. U anticommuting with L_e1..L_e7
    d_anti = intertwiner_dim([Lint[t] for t in range(1, 8)], [im_lin(Lint[t], Lint[t], -1, 0) for t in range(1, 8)], 8)
    d_lr = intertwiner_dim([Lint[t] for t in range(8)], [Rint[t] for t in range(8)], 8)
    d_lrc = intertwiner_dim([Lint[t] for t in range(8)], [Rint[0]] + [im_lin(Rint[t], Rint[t], -1, 0) for t in range(1, 8)], 8)
    supp_row("S8.02", "only_U=0_solves_U_L_xi=L_conj(xi)U", d_anti == 0,
             [("dim_U_with_U_L_xi=L_conj(xi)_U", d_anti), ("dim_U_with_U_L_xi=R_xi_U", d_lr), ("dim_U_with_U_L_xi=R_conj(xi)_U", d_lrc),
              ("reading", "from xi=e0 any V with U N(xi) V = L_conj(xi) is V = K U^-1, so U L_xi U^-1 = L_conj(xi) is needed; dim 0 means no S+-and-S- preserving identification")])
    nc = ctx["N_poly"](o_conj(xi), psi_m, -1)
    res = zero_check(flat(pm_lin(pm_mul(nc, pm_from_int(kint)), lc, 1, -1)), B_XI1)
    supp_row("S8.03", "N(conj(xi))K=L_conj(xi)_K(a)=-conj(a)", holds(res),
             cpairs(res) + [("reading", "conj(xi)a is the CONV symbol after the reflection xi -> conj(xi) of R^8 (det -1) and a -> -conj(a)")])
    # S8.04 the block form of ADD R3
    gpp = pm_zero(16, 16)
    for r in range(8):
        for c in range(8):
            gpp[r][8 + c] = dict(lc[r][c])
            gpp[8 + r][c] = dict(lm[r][c])
    g2pp = pm_mul(gpp, gpp)
    rp = zero_check(flat(pm_lin(g2pp, pm_diag(nsq, 16), 1, -1)), B_XI2)
    rmn = zero_check(flat(pm_lin(g2pp, pm_diag(nsq, 16), 1, 1)), B_XI2)
    dsg = im_ident(16)
    for r in range(8, 16):
        dsg[r][r] = -1
    rel = zero_check(flat(pm_lin(gpp, pm_mul(pm_from_int(dsg), gamma_poly(tab, o_conj(xi), -1)), 1, -1)), B_XI1)
    gens_pp = [coeff_int(gpp, t) for t in range(8)]
    vol = im_ident(16)
    for t in range(8):
        vol = im_mul(vol, gens_pp[t])
    alg_pp = algebra_of(gens_pp)
    supp_row("S8.04", "ADD_R3_block_[[0,L_conj(xi)],[L_xi,0]]", holds(rp),
             cpairs(rp) + [("squares_to", "+nsq(xi)I16" if holds(rp) else ("-nsq(xi)I16" if holds(rmn) else "other")),
                           ("equals_diag(I8,-I8)gamma(conj(xi))", 1 if holds(rel) else 0),
                           ("volume_element_S_plus_entry", vol[0][0]), ("generated_algebra_dim", alg_pp)])
    # S8.05 the two products of the associator control
    e = [obasis(t) for t in range(8)]
    lhs = [p_eval(c, {}) if c else 0 for c in o_mul(tab, o_mul(tab, e[1], e[2]), e[3])]
    rhs = [p_eval(c, {}) if c else 0 for c in o_mul(tab, e[1], o_mul(tab, e[2], e[3]))]
    supp_row("S8.05", "associator_products_(e1e2)e3_and_e1(e2e3)", True, [("(e1e2)e3", o_text(lhs)), ("e1(e2e3)", o_text(rhs))])
    # S8.06 the doubled self-adjoint rolled-up symbol
    dbl = pm_zero(16, 16)
    for r in range(8):
        for c in range(8):
            dbl[r][8 + c] = dict(N[c][r])
            dbl[8 + r][c] = dict(N[r][c])
    rd = zero_check(flat(pm_lin(pm_mul(dbl, dbl), pm_diag(nsq, 16), 1, -1)), B_XI2)
    supp_row("S8.06", "doubled_[[0,N^T],[N,0]]^2=+nsq(xi)I16", holds(rd), cpairs(rd))
    # S8.07 the block formula of 8D.C2
    bad = 0
    for i in range(8):
        for j in range(8):
            want = im_zero(16, 16)
            a = im_mul(Lint[i], Lcint[j])
            b = im_mul(Lcint[i], Lint[j])
            for r in range(8):
                for c in range(8):
                    want[r][c] = -a[r][c]
                    want[8 + r][8 + c] = -b[r][c]
            if im_mul(gam[i], gam[j]) != want:
                bad += 1
    supp_row("S8.07", "8D.C2_gamma_i_gamma_j=diag(-L_ei_L_conj(ej),-L_conj(ei)_L_ej)", bad == 0, [("ordered_pairs", 64), ("failed", bad)])
    # S8.08 traces, S8.09 the eigenvector formula of 8D.P6, S8.10 and S8.11 Psi against T_s
    tr = {s: sum(T8[s][q][q] for q in range(28)) for s in (1, -1)}
    supp_row("S8.08", "trace_T_s", tr[1] == 0 and tr[-1] == 0, [("trace_T_+1", tr[1]), ("trace_T_-1", tr[-1])])

    def phi_val(i, j, k):
        ss, key = sort_sign((i, j, k))
        return ss * phi.get(key, 0) if ss else 0

    def iota_vec(k):
        return [phi_val(k, i, j) if i > 0 else 0 for (i, j) in PAIRS8]
    e0 = [[1 if (i == 0 and j == k) else 0 for (i, j) in PAIRS8] for k in range(8)]
    fb = {1: 0, -1: 0}
    for s in (1, -1):
        for k in range(1, 8):
            v = [-3 * x + y for x, y in zip(e0[k], iota_vec(k))]
            want = [(2 - 3 * s) * y + 3 * s * x for x, y in zip(e0[k], iota_vec(k))]
            if im_vec(T8[s], v) != want:
                fb[s] += 1
    supp_row("S8.09", "8D.P6_T_s(-3e0^u+iota_u_phi)=(2-3s)iota_u_phi+3s_e0^u", fb[1] == 0 and fb[-1] == 0,
             [("u", "e1..e7"), ("failed_s_+1", fb[1]), ("failed_s_-1", fb[-1])])
    four_psi = im_lin(psi_m, psi_m, 4, 0)
    h1 = im_mul(psi_m, im_lin(T8[1], i28, 1, 1)) == four_psi
    h2 = im_mul(psi_m, im_lin(T8[-1], i28, 1, 1)) == four_psi
    supp_row("S8.10", "Psi(T_+1+I)=4Psi_and_Psi(T_-1+I)=4Psi", h1 and not h2,
             [("Psi(T_+1+I)=4Psi", 1 if h1 else 0), ("Psi(T_-1+I)=4Psi", 1 if h2 else 0)])
    moved = sum(1 for v in ker_psi if any(im_vec(psi_m, im_vec(T8[-1], v))))
    supp_row("S8.11", "ker_Psi_under_T_-1", moved > 0, [("basis_vectors_moved_out_of_ker_Psi", moved), ("of", len(ker_psi))])
    e0j = all(psi_p[j - 1][PAIRS8.index((0, j))] == -1 and sum(abs(psi_p[k][PAIRS8.index((0, j))]) for k in range(7)) == 1 for j in range(1, 8))
    psi7 = ctx["psi7"]
    restr = all(psi_p[k][p] == psi7[k][p] for k in range(7) for p, (i, j) in enumerate(PAIRS8) if i > 0)
    supp_row("S8.12", "Psi'(e0j)=-e_j_and_Psi'_on_so7=contraction_with_phi", e0j and restr,
             [("Psi'(e0j)=-e_j", 1 if e0j else 0), ("Psi'_on_so7", 1 if restr else 0)])
    l27 = all(psi_m[k - 1] == [x + y for x, y in zip(e0[k], iota_vec(k))] for k in range(1, 8))
    supp_row("S8.13", "Lambda2_7=rowspace_Psi=(e0^u+iota_u_phi)", l27, [("rows", 7)])
    ru = zero_check(flat(pm_lin(N, pm_mul(ctx["lm"], pm_from_int(kint)), 1, -1)), B_XI1)
    ga = gamma_poly(tab, av, -1)
    ca = o_conj(av)
    eta_ok = all(not ga[r][0] for r in range(8)) and all(ga[8 + t][0] == p_scale(ca[t], -1) for t in range(8))
    supp_row("S8.14", "N(xi)=L_xi_K_and_gamma(a)eta=(0,-conj(a))", holds(ru) and eta_ok,
             cpairs(ru) + [("gamma(a)eta=(0,-conj(a))", 1 if eta_ok else 0)])
    sm_bad = 0
    c7s = []
    csm = []
    for v in ker_psi:
        cm = cliff(v)
        c7s.append([[cm[r][c] for c in range(1, 8)] for r in range(1, 8)])
        cs = [[cm[8 + r][8 + c] for c in range(8)] for r in range(8)]
        csm.append(cs)
        if im_mul(cs, kint) != im_lin(im_mul(kint, R_of(v)), im_zero(8, 8), 2, 0):
            sm_bad += 1
    d7 = intertwiner_dim(c7s, c7s, 7)
    d8 = intertwiner_dim(csm, csm, 8)
    supp_row("S8.15", "spin7_reps_S-=Lambda1_via_K_and_commutants", sm_bad == 0 and d7 == 1 and d8 == 1,
             [("c(omega)_S-_K=2K_R_omega_failed", sm_bad), ("commutant_on_ImO_of_S+", d7), ("commutant_on_S-", d8),
              ("reading", "commutant dimension 1 for a compact Lie algebra means irreducible of real type")])
    ker_pp = kernel_int(psi_p, 28)
    rk_sum = rank_int(ker_psi + ker_pp, 28)
    supp_row("S8.16", "dim(ker_Psi+ker_Psi')", rk_sum == 28, [("dim", rk_sum)])
    yv = [p_var(k) for k in range(7)]
    om_y = []
    for p in range(28):
        acc = {}
        for k in range(7):
            if psi_m[k][p]:
                p_add_into(acc, yv[k], psi_m[k][p])
        om_y.append(acc)
    ps_y = []
    for k in range(7):
        acc = {}
        for p in range(28):
            if psi_m[k][p]:
                p_add_into(acc, om_y[p], psi_m[k][p])
        ps_y.append(acc)
    n1 = {}
    n2 = {}
    for x in ps_y:
        p_add_into(n1, p_mul(x, x))
    for x in om_y:
        p_add_into(n2, p_mul(x, x))
    r19 = zero_check([p_sub(n1, p_scale(n2, 4))], basis([(list(range(7)), 2)]))
    supp_row("S8.17", "nsq(Psi(omega))=4nsq(omega)_on_Lambda2_7", holds(r19), cpairs(r19) + [("omega", "Psi^T y, y in R^7")])
    cols = ctx["wedge_cols"](xi)
    tpi = ctx["tpi"]
    tw = []
    for j in range(8):
        v = []
        for p in range(28):
            acc = {}
            for q in range(28):
                if tpi[p][q] and cols[j][q]:
                    p_add_into(acc, cols[j][q], tpi[p][q])
            v.append(acc)
        tw.append(v)
    ent = []
    for i in range(8):
        for j in range(8):
            acc = {}
            for p in range(28):
                if cols[i][p] and tw[j][p]:
                    p_add_into(acc, p_mul(cols[i][p], tw[j][p]))
            want = p_sub(nsq if i == j else {}, p_mul(xi[i], xi[j]))
            ent.append(p_sub(acc, want))
    r20 = zero_check(ent, B_XI2)
    supp_row("S8.18", "R1_matrix_form_W^T(T+I)W=nsq(xi)I8-xi_xi^T", holds(r20),
             cpairs(r20) + [("reading", "sigma^T sigma = xi xi^T + (1/4)(nsq I - xi xi^T) with pi = (T+I) times 1/4, so c = 1/4")])


def read_frozen(path, want_sha):
    with open(path, "rb") as fh:
        data = fh.read()
    sha = hashlib.sha256(data).hexdigest()
    rows = []
    for line in data.decode("utf-8").split("\n"):
        if not line:
            continue
        f = line.split("|")
        d = {"_tag": f[0]}
        for x in f[1:]:
            if "=" in x:
                k, v = x.split("=", 1)
                d[k] = v
        rows.append(d)
    return sha, sha == want_sha, rows


def verdict(frozen, says, measured, msays, v, note):
    VERDICTS[v].append(frozen)
    pairs = [("frozen", frozen), ("frozen_says", says), ("measured", measured), ("measured_says", msays), ("verdict", v)]
    if note:
        pairs.append(("note", note))
    row("COMPARE", pairs + [("claim", "MEASURED")])


def compare(pred_path, add_path, ctx):
    psha, pok, prow = read_frozen(pred_path, PRED_SHA)
    asha, aok, arow = read_frozen(add_path, ADD_SHA)
    row("FROZEN", [("file", "PREDICTION-SPINOR-SYMBOL.hbp"), ("sha256", psha), ("matches_freeze", 1 if pok else 0), ("rows", len(prow))])
    row("FROZEN", [("file", "ADDENDUM-SPINOR-SYMBOL.hbp"), ("sha256", asha), ("matches_freeze", 1 if aok else 0), ("rows", len(arow))])
    if not (pok and aok):
        row("COMPARE_SUMMARY", [("status", "FROZEN_SHA_MISMATCH")])
        return
    supplement(ctx)
    R = RESULTS
    S = SUPP
    pred = {d["id"]: d for d in prow if d["_tag"] in ("PRED", "CONV") and "id" in d}
    add = {d["id"]: d for d in arow if d["_tag"] == "ADD"}

    def ok(*ids):
        return all((R[i]["ok"] if i in R else S[i]["ok"]) for i in ids)

    p = pred["8D.P1"]["predicted"]
    v = "AGREE" if (ok("M8.B1") and R["M8.B1"]["coeffs_checked"] == 9216 and "9216" in p) else "DISAGREE"
    verdict("8D.P1", p, "M8.B1", "9216 coefficients checked, 0 nonzero, 0 outside basis", v, "")
    p = pred["8D.P2"]["predicted"]
    v = "AGREE" if (ok("M8.B5") and "36 of 36" in p) else "DISAGREE"
    verdict("8D.P2", p, "M8.B5", "64 ordered pairs (the 36 with i<=j among them) hold exactly", v, "")
    p = pred["8D.P3"]["predicted"]
    v = "AGREE" if (R["M8.C1"]["dimension"] == 256 and R["M8.C2"]["dimension"] == 1 and "256" in p and "commutant dimension 1" in p) else "DISAGREE"
    verdict("8D.P3", p, "M8.C1 M8.C2", "dimension 256 = M_16(R); commutant dimension 1", v, "")
    p = pred["8D.P4"]["predicted"]
    e7 = S["S8.01"]["eps7"]
    v = "AGREE" if (ok("M8.E1", "S8.01") and e7 == R["M8.E1"]["S_plus_eigenvalue"]) else "DISAGREE"
    verdict("8D.P4", p, "M8.E1 S8.01", "gamma0...gamma7 = diag(-I8, I8); L_e1...L_e7 = " + str(e7) + " I8, so epsilon7 = " + str(e7) +
            "; blocks equal L_e1...L_e7 and its negative", v, "epsilon7 measured on the frozen table; 7D.P4 left its sign open")
    p = pred["8D.P5"]["predicted"]
    v = "AGREE" if ok("M8.G1", "M8.G2", "M8.G3", "M8.G5", "S8.12", "S8.13") else "DISAGREE"
    verdict("8D.P5", p, "M8.G1 M8.G2 M8.G3 M8.G5 S8.12 S8.13",
            "Psi(e0j)=e_j, Psi(e_ij)=e_i e_j, rank 7, Psi Psi^T=4 I7, dim ker 21, Lambda2_7 = rows e0^u + iota_u phi, Psi'(e0j)=-e_j", v, "")
    p = pred["8D.P6"]["predicted"]
    g4 = R["M8.G4"]
    v = "AGREE" if (ok("M8.F5p", "M8.F5m", "S8.08", "S8.09", "S8.10", "S8.11") and g4["Psi_matches_s"] == 1 and g4["Psi'_matches_s"] == -1
                    and R["M8.G5"]["dim_kerPsi_cap_stab_Phi_+1"] == 21 and R["M8.G5"]["dim_kerPsi'_cap_stab_Phi_-1"] == 21) else "DISAGREE"
    verdict("8D.P6", p, "M8.F5p M8.F5m M8.G4 M8.G5 S8.08 S8.09 S8.10 S8.11",
            "both s: (T-3I)(T+I)=0, eigen dims 7 and 21, trace 0; eigenvector formula holds; Psi(T_+1+I)=4Psi holds, Psi(T_-1+I)=4Psi fails; "
            "T_+1+I = Psi^T Psi, T_-1+I = Psi'^T Psi'; ker Psi = stab Phi_+1, ker Psi' = stab Phi_-1", v, "")
    p = pred["8D.P7"]["predicted"]
    v = "AGREE" if (ok("M8.G7", "S8.16") and R["M8.G7"]["dim_ker_Psi_cap_ker_Psi'"] == 14) else "DISAGREE"
    verdict("8D.P7", p, "M8.G7 S8.16", "intersection dim 14, inside so(7), equal to beta with beta contracted with phi = 0 and to stab(phi); sum dim 28", v, "")
    p = pred["8D.P8"]["predicted"]
    v = "AGREE" if (ok("M8.I1", "M8.I2", "M8.I3", "M8.I4", "M8.I5") and R["M8.J3"]["holding"] == "-xi*conj(a)") else "DISAGREE"
    verdict("8D.P8", p, "M8.I1..I5 M8.J3", "dims 1,8,7, Euler 0, P xi=0, P P^T=nsq I7, P^T P=nsq I8 - xi xi^T, ranks 1 and 7; components 1..7 of -xi conj(a) = P(xi)a", v, "")
    p = pred["8D.P9"]["predicted"]
    v = "AGREE" if (ok("M8.J1", "M8.J2", "M8.K3", "S8.14", "S8.17") and R["M8.J3"]["holding"] == "-xi*conj(a)" and "N(xi) a = -xi conj(a)" in p) else "DISAGREE"
    verdict("8D.P9", p, "M8.J1 M8.J2 M8.J3 M8.K3 S8.14 S8.17",
            "N(xi)a = -xi conj(a), the only one of 16 candidates; N = L_xi K; gamma(a)eta = (0,-conj(a)); N^T N = N N^T = nsq I8; "
            "nsq(Psi omega) = 4 nsq(omega) on Lambda2_7; nsq((T+I)(xi^a)) = 4 nsq(xi^a)", v, "")
    p = pred["8D.P10"]["predicted"]
    v = "AGREE" if ok("M8.H2", "M8.H4", "M8.G6", "M8.M2", "S8.07") else "DISAGREE"
    verdict("8D.P10", p, "M8.H2 M8.H3 M8.H4 M8.G6 M8.M2 S8.07",
            "c(omega)eta=0 on 21 generators; 2(N(R xi)+N(xi)R) = c(omega)_S+ N(xi); [c(omega),gamma(xi)]=2gamma(R_omega xi) for all 28; closed, dim 21; "
            "c(omega)_S+ = -sum omega_ij L_ei L_conj(ej) from the 8D.C2 blocks", v, "")
    p = pred["8D.P11"]["predicted"]
    v = "AGREE" if ok("M8.H4", "M8.H5", "M8.G5", "S8.15") else "DISAGREE"
    verdict("8D.P11", p, "M8.H4 M8.H5 M8.G5 S8.15",
            "eta fixed (1), Psi intertwines Lambda2_7 with Im O (7), S- = Lambda1 through K (8), commutants 1 and 1, Lambda2_21 = ker Psi = stab(Phi_+1)", v,
            "irreducibility measured by commutant dimension; the triality reason is NAMED, not run")
    dims_row = [d for d in prow if d["_tag"] == "DIMS" and d.get("case") == "8D"][0]
    mv = {"n": 8, "spinor_module": 16, "clifford_dim": 256 if R["M8.C1"]["gram_failures"] == 0 else -1,
          "generated_algebra": R["M8.C1"]["dimension"], "commutant": R["M8.C2"]["dimension"], "generators": 8, "generator_size": 16,
          "lambda2_small": R["M8.F5p"]["rank_T+I"], "lambda2_big": R["M8.F5p"]["rank_T-3I"], "T_small": 3 if ok("M8.F5p") else 0,
          "T_big": -1 if ok("M8.F5p") else 0, "k": 4 if R["M8.G1"]["Psi_PsiT"] == "4*I7" else 0, "stabilizer_of_eta": R["M8.G5"]["dim_ker_Psi"],
          "sequence_dims": R["M8.I1"]["dims"], "symbol_ranks": R["M8.I1"]["ranks"], "square_coefficients": R["M8.B1"]["coeffs_checked"]}
    for k in ("n", "spinor_module", "clifford_dim", "generated_algebra", "commutant", "generators", "generator_size", "lambda2_small",
              "lambda2_big", "T_small", "T_big", "k", "stabilizer_of_eta", "sequence_dims", "symbol_ranks", "square_coefficients"):
        v = "AGREE" if str(mv[k]) == dims_row[k] else "DISAGREE"
        verdict("DIMS.8D." + k, dims_row[k], "measured block", str(mv[k]), v, "")
    bott = [d for d in prow if d["_tag"] == "BOTT" and d.get("n") == "8"][0]
    v = "AGREE" if (bott["cl_square_minus_1"] == "M_16(R)" and bott["dim_minus"] == "256" and bott["irreducible_module_minus"] == "16"
                    and R["M8.C1"]["algebra_equals_M16R"] == 1 and R["M8.C2"]["dimension"] == 1) else "DISAGREE"
    verdict("BOTT.8.minus", bott["cl_square_minus_1"] + " dim " + bott["dim_minus"] + " module " + bott["irreducible_module_minus"],
            "M8.C1 M8.C2", "M_16(R), dim 256, on R^16 with commutant 1", v, "")
    v = "AGREE" if (bott["cl_square_plus_1"] == "M_16(R)" and bott["dim_plus"] == "256" and S["S8.04"]["generated_algebra_dim"] == 256 and S["S8.04"]["ok"]) else "DISAGREE"
    verdict("BOTT.8.plus", bott["cl_square_plus_1"] + " dim " + bott["dim_plus"], "S8.04",
            "generators squaring to +1 generate dim 256 inside M_16(R)", v, "")
    v = "AGREE" if ok("S8.07") else "DISAGREE"
    verdict("8D.C2", "gamma(e_i) gamma(e_j) = diag(-L_(e_i) L_conj(e_j), -L_conj(e_i) L_(e_j))", "S8.07", "64 of 64 ordered pairs", v, "")
    v = "AGREE" if (ok("M8.F4", "M8.F5p", "M8.F5m")) else "DISAGREE"
    verdict("8D.C3", "*8(e0 wedge alpha) = *7 alpha; Phi_+1 self-dual, Phi_-1 anti-self-dual", "M8.F4 M8.F5p M8.F5m",
            "35 of 35 alphas; star8 Phi_s = s Phi_s for both s", v, "")
    verdict("7D.C3", "if e1 wedge ... wedge e7 gives -2 and +1, the measurement reverses it and says so", "M8.F2",
            R["M8.F2"]["std_orientation"], "CONVENTION_FOLLOWED", "vol7 = -e1234567 and vol8 = e0 wedge vol7 are used throughout")
    sc = [d for d in prow if d["_tag"] == "SELF_CHECK" and "8D" in d.get("name", "")]
    rk_seq = [int(x) for x in mv["symbol_ranks"].split(",")]
    smv = {"small": mv["lambda2_small"], "big": mv["lambda2_big"], "binom_n_2": mv["lambda2_small"] + mv["lambda2_big"],
           "trace": S["S8.08"]["trace_T_+1"], "k": mv["k"], "k_times_small": mv["k"] * mv["lambda2_small"],
           "T_small": mv["T_small"], "eigenvalues": str(mv["T_small"]) + "," + str(mv["T_big"]), "stab": mv["stabilizer_of_eta"],
           "group_minus_stab": 28 - mv["stabilizer_of_eta"], "dim_S_plus_minus_1": R["M8.G5"]["rank_Psi"], "seq": mv["sequence_dims"],
           "euler": R["M8.I1"]["euler"], "solved": mv["symbol_ranks"], "residual": 8 - rk_seq[0] - rk_seq[1], "predicted": mv["symbol_ranks"],
           "even_sum": 1 + mv["lambda2_small"], "odd_sum": 8,
           "dim_S_half": 8, "cl": mv["clifford_dim"], "two_pow_n": 256, "sum_binom": 256, "bott_entry": "M_16(R)", "alg": mv["generated_algebra"],
           "commutant": mv["commutant"], "module": 16, "alg_times_commutant": mv["generated_algebra"] * mv["commutant"], "module_squared": 256,
           "summands": 1, "gens": 8, "pairs": 36, "monomials_deg2": 36, "coefficients": mv["square_coefficients"], "one_plus_seven": 1 + mv["lambda2_small"],
           "S_minus": 8, "sum": 16, "cl8": mv["generated_algebra"]}
    for d in sc:
        keys = [k for k in d if k not in ("_tag", "name", "status", "claim", "json")]
        miss = [k for k in keys if k not in smv]
        diff = [k for k in keys if k in smv and str(smv[k]) != d[k]]
        if miss and len(miss) == len(keys):
            v = "NOT_MEASURED"
        elif diff:
            v = "DISAGREE"
        else:
            v = "AGREE"
        verdict("SELF_CHECK." + d["name"], ",".join(k + ":" + d[k] for k in keys), "measured block and SUPP",
                ",".join(k + ":" + str(smv[k]) for k in keys if k in smv), v, ("not measured here: " + ",".join(miss)) if miss else "")
    a = add["R1"]
    v = "AGREE" if (ok("M8.K2", "M8.K3", "S8.18") and R["K8.01"]["observed_fail"] and R["K8.02"]["observed_fail"] and "1:4(Spin(7))" in a["predicted"]) else "DISAGREE"
    verdict("ADD.R1.8D", a["predicted"], "M8.K2 M8.K3 S8.18 K8.01 K8.02",
            "c = 1/4 holds as a polynomial identity (sigma^T sigma = xi xi^T + (1/4)(nsq I - xi xi^T)); c = 1/2 and c = 1/3 fail", v, "")
    a = add["R2"]
    v = "AGREE" if ok("M8.J1", "M8.J2") else "DISAGREE"
    verdict("ADD.R2.8D_isometry", "sigma^T sigma = nsq(xi) I exactly with integer matrices", "M8.J1 M8.J2", "N^T N = N N^T = nsq(xi) I8", v, "")
    v = "AGREE" if R["M8.J4"]["ok"] else "DISAGREE"
    verdict("ADD.R2.8D_form", "8D Lambda2_7 through the Cayley form, so that sigma(xi)a = conj(xi)a", "M8.J4 M8.J5",
            "in the CONV identification (Phi_+1, eta in S+, -iota_xi, Psi coordinates) N(xi)a = -xi conj(a) and conj(xi)a fails on 50 of 512 coefficients; "
            "-conj(xi)a holds with Psi' (eta' in S-, Phi_-1), so conj(xi)a = -N'(xi)a", v,
            "the addendum itself says 8D.P9 governs and asks for this row")
    v = "AGREE" if (S["S8.02"]["dim_U_with_U_L_xi=L_conj(xi)_U"] > 0) else "DISAGREE"
    verdict("ADD.R2.8D_form_note", "conj(xi)a and -xi*conj(a) are the same symbol under two identifications of S+- with O (left versus right multiplication and a -> -conj(a) on one factor)",
            "S8.02 S8.03 M8.J5", "fixed U with U L_xi = L_conj(xi) U: dimension " + str(S["S8.02"]["dim_U_with_U_L_xi=L_conj(xi)_U"]) +
            ", so no fixed U, V give U N(xi) V = L_conj(xi), the map a -> conj(xi)a; U L_xi = R_xi U: dimension " + str(S["S8.02"]["dim_U_with_U_L_xi=R_xi_U"]) +
            "; U L_xi = R_conj(xi) U: dimension " + str(S["S8.02"]["dim_U_with_U_L_xi=R_conj(xi)_U"]) +
            " (conjugation), and flipping the Lambda0 sign of N gives a -> a conj(xi) exactly (M8.J5); conj(xi)a itself is reached only through the reflection "
            "xi -> conj(xi) of R^8 (S8.03, det -1, exchanging Phi_+1 and Phi_-1) or through eta' in S- (M8.J5: -conj(xi)a with Psi')", v,
            "as written the note fails; with conj(xi)a read as a conj(xi) (right multiplication by conj(xi)) it holds under a fixed identification; both forms are isometries with N^T N = nsq I")
    a = add["R3"]
    v = "AGREE" if ok("M8.L1") else "DISAGREE"
    verdict("ADD.R3.relation", "conj(xi)(xi x) = nsq(xi) x for all octonions xi, x", "M8.L1", "2304 coefficients, 0 failed; xi(xi x) = -nsq(xi)x fails (K8.04)", v, "")
    v = "AGREE" if (S["S8.04"]["squares_to"] == "-nsq(xi)I16") else "DISAGREE"
    verdict("ADD.R3.block", "the Cl(8) generator on S+ + S- is [[0, L_conj(xi)],[L_xi, 0]]", "S8.04",
            "it squares to " + S["S8.04"]["squares_to"] + " and equals diag(I8,-I8) gamma(conj(xi)); a generator of Cl(R^8, +nsq), not of the CONV Cl(R^8, -nsq) of C0.1 and 8D.C2",
            v, "a sign convention difference: R6 itself says the positive-square forms are reported in their own rows; under C0.1 this block is not a Clifford generator")
    v = "AGREE" if ok("M8.L4", "M8.L5") else "DISAGREE"
    verdict("ADD.R3.7D", "xi(xi x) = (xi xi)x = -nsq(xi) x for imaginary xi", "M8.L4 M8.L5", "both hold as polynomial identities", v, "")
    a = add["R4"]
    v = "AGREE" if (R["K8.05"]["(e1e2)e3-e1(e2e3)"] == "-2e6" and S["S8.05"]["(e1e2)e3"] == "-1e6" and S["S8.05"]["e1(e2e3)"] == "1e6" and ok("M8.A1")) else "DISAGREE"
    verdict("ADD.R4", a["predicted"], "K8.05 S8.05 M8.A1", "(e1e2)e3 = -e6, e1(e2e3) = +e6, associator -2e6; 28 of 35 triples non-associative, the 7 lines associative", v, "")
    v = "AGREE" if (ok("M8.I1", "M8.I2") and R["M8.I1"]["euler"] == 0) else "DISAGREE"
    verdict("ADD.R5.8D_symbol", "sigma(xi)sigma(xi)=0 whatever F_A is; Euler count 1-8+7=0 (8D)", "M8.I1 M8.I2", "P(xi)xi = 0 as a polynomial identity; Euler 0", v, "")
    verdict("ADD.R5.operator", "the operator sequence is a complex only on instantons, pi_7 F = 0", "none", "outside a symbol measurement", "NOT_MEASURED", "")
    v = "AGREE" if (ok("M8.B1", "M8.J1", "S8.06")) else "DISAGREE"
    verdict("ADD.R6", "gamma(xi)^2 = -nsq(xi) I; the rolled-up N and its doubled self-adjoint form square to +nsq(xi)", "M8.B1 M8.J1 S8.06",
            "gamma^2 = -nsq I16; N^T N = +nsq I8; [[0,N^T],[N,0]]^2 = +nsq I16", v, "")
    verdict("ADD.DISCLOSURE", "GPT-6 reported 18 of 18 symbol checks, receipt c28f8af6", "none", "receipt not read by this seat, by design", "NOT_COMPARED", "")
    row("COMPARE_SUMMARY", [("agree", len(VERDICTS["AGREE"])), ("disagree", len(VERDICTS["DISAGREE"])),
                            ("not_measured", len(VERDICTS["NOT_MEASURED"])), ("not_compared", len(VERDICTS["NOT_COMPARED"])),
                            ("convention_followed", len(VERDICTS["CONVENTION_FOLLOWED"])),
                            ("disagreements", " ".join(VERDICTS["DISAGREE"]) if VERDICTS["DISAGREE"] else "none"),
                            ("not_measured_items", " ".join(VERDICTS["NOT_MEASURED"]) if VERDICTS["NOT_MEASURED"] else "none"),
                            ("measured_block_written_before_reading", 1)])


def main():
    if len(sys.argv) not in (2, 4):
        print("usage: measure_8d.py OUTDIR [PRED.hbp ADD.hbp]")
        return 2
    outdir = sys.argv[1]
    with open(__file__, "rb") as fh:
        self_sha = hashlib.sha256(fh.read()).hexdigest()
    vi = sys.version_info
    row("HEAD", [("seat", "opus55-acer"), ("case", "8D Spin(7)"), ("script", "measure_8d.py"), ("script_sha256", self_sha),
                 ("python", str(vi[0]) + "." + str(vi[1]) + "." + str(vi[2])), ("integer_only", 1), ("float_used", 0),
                 ("read_before_measuring", "CONV and FORMAT rows of PREDICTION-SPINOR-SYMBOL.hbp only"),
                 ("method", "integer polynomial identities over the full monomial basis; sampled xi only in controls")])
    ctx = measure()
    block = "\n".join(ROWS[1:]) + "\n"
    block_sha = hashlib.sha256(block.encode("utf-8")).hexdigest()
    row("MEASURED_BLOCK_END", [("rows", len(ROWS) - 1), ("measured_pass", TALLY["pass"]), ("measured_fail", TALLY["fail"]),
                               ("controls_ok", TALLY["ctrl_ok"]), ("controls_bad", TALLY["ctrl_bad"]),
                               ("block_sha256", block_sha), ("block", "rows between HEAD and this row")])
    if len(sys.argv) == 4:
        same = block_sha == PRE_READ_BLOCK_SHA
        row("PRE_READ", [("block_sha256_before_reading", PRE_READ_BLOCK_SHA), ("reproduced_byte_for_byte", 1 if same else 0),
                         ("file_before_reading", "MEASURED-ONLY-BEFORE-READ-8D.hbp"), ("file_sha256", PRE_READ_FILE_SHA),
                         ("script_sha256_before_reading", PRE_READ_SCRIPT_SHA),
                         ("after_reading", "only the SUPP and COMPARE rows below were added; the measured rows above are unchanged")])
        if not same:
            row("COMPARE_SUMMARY", [("status", "MEASURED_BLOCK_CHANGED_AFTER_READING")])
        else:
            compare(sys.argv[2], sys.argv[3], ctx)
    row("END", [("rows", len(ROWS) + 1)])
    data = ("\n".join(ROWS) + "\n").encode("utf-8")
    path = os.path.join(outdir, OUT_NAME)
    with open(path, "wb") as fh:
        fh.write(data)
    with open(path + ".sha256", "wb") as fh:
        fh.write((hashlib.sha256(data).hexdigest() + "  " + OUT_NAME + "\n").encode("utf-8"))
    print("WROTE|path=" + path + "|rows=" + str(len(ROWS)) + "|sha256=" + hashlib.sha256(data).hexdigest() + "|json=0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
