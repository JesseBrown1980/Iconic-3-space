#!/usr/bin/env python3
"""
measure_4d.py: the 4D case of the spinor-symbol measurement, Iconic-3-space, November 9, 2026.

Written from the CONV rows of predictions/spinor-symbol/PREDICTION-SPINOR-SYMBOL.hbp on branch freeze-spinor-symbol-20261109
(file sha256 c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e). Its HEAD, FORMAT and CONV rows were read first;
no PRED, ADD or DIMS value was read before the measured block below was written and run.

Integers only. Every identity in xi is compared as integer (or Gaussian-integer) polynomials, monomial by monomial; a sampled
covector appears only in a control. Ranks are exact, by fraction-free integer elimination.

Run: python3 -I -B measure_4d.py                  measured block, compare section, summary: MEASURE-4D.hbp and its .sha256
     python3 -I -B measure_4d.py --measure-only   the measured block alone: MEASURE-4D-BEFORE-COMPARE.hbp and its .sha256
"""
import sys
import os
import hashlib
import itertools

if __name__ == "__main__" and not sys.flags.isolated:
    print("MEASURE_4D|needs_python3_-I=1|json=0")
    sys.exit(1)

HERE = os.path.dirname(os.path.abspath(__file__))
MEASURE_ONLY = "--measure-only" in sys.argv[1:]
FROZEN_SHA = "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e"
ADDENDUM_SHA = "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0"

LINES = []


def emit(tag, pairs):
    parts = [tag]
    for k, v in pairs:
        s = str(v)
        if "|" in s or "\n" in s or "|" in k:
            raise SystemExit("bad field " + k)
        parts.append(k + "=" + s)
    parts.append("json=0")
    LINES.append("|".join(parts))


STATS = {"meas": 0, "meas_pass": 0, "ctrl": 0, "ctrl_ok": 0}


def meas(mid, name, ok, pairs):
    STATS["meas"] += 1
    if ok:
        STATS["meas_pass"] += 1
    emit("MEAS", [("id", mid), ("name", name)] + pairs + [("pass", 1 if ok else 0), ("claim", "MEASURED")])


def control(cid, name, failed, pairs):
    STATS["ctrl"] += 1
    if failed:
        STATS["ctrl_ok"] += 1
    emit("CONTROL", [("id", cid), ("name", name)] + pairs
         + [("expected", "FAIL"), ("observed", "FAIL" if failed else "PASS"), ("ok", 1 if failed else 0), ("claim", "MEASURED")])


def lst(v):
    return ",".join(str(x) for x in v)


def mat_str(A):
    return "[" + "; ".join(lst(r) for r in A) + "]"


# ---------------------------------------------------------------- integers

def igcd(a, b):
    a = abs(a)
    b = abs(b)
    while b:
        a, b = b, divmod(a, b)[1]
    return a


def exact_div(a, b):
    q, r = divmod(a, b)
    if r:
        raise SystemExit("inexact division")
    return q


def reduced(n, d):
    g = igcd(n, d)
    n = exact_div(n, g)
    d = exact_div(d, g)
    if d < 0:
        n, d = -n, -d
    return n, d


def rank_int(rows):
    m = [list(r) for r in rows if any(r)]
    if not m:
        return 0
    ncol = len(m[0])
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(m)):
            if m[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        p = m[r][c]
        for i in range(r + 1, len(m)):
            v = m[i][c]
            if v:
                new = [p * m[i][t] - v * m[r][t] for t in range(ncol)]
                g = 0
                for t in new:
                    g = igcd(g, t)
                if g > 1:
                    new = [exact_div(t, g) for t in new]
                m[i] = new
        r += 1
        if r == len(m):
            break
    return r


def izero(n, m):
    return [[0] * m for _ in range(n)]


def ieye(n, k=1):
    return [[k if i == j else 0 for j in range(n)] for i in range(n)]


def imul(A, B):
    k = len(B)
    m = len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(len(A))]


def iadd(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def isub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def iscale(A, k):
    return [[k * x for x in r] for r in A]


def iT(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def iblock(A, B, C, D):
    return [A[i] + B[i] for i in range(len(A))] + [C[i] + D[i] for i in range(len(C))]


def flat(A):
    return [x for r in A for x in r]


def matvec_i(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


# ---------------------------------------------------------------- integer polynomials: dict sorted-var-tuple -> int

def padd(p, q):
    r = dict(p)
    for m, c in q.items():
        v = r.get(m, 0) + c
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def pscale(p, k):
    if k == 0:
        return {}
    return {m: c * k for m, c in p.items()}


def psub(p, q):
    return padd(p, pscale(q, -1))


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2))
            v = r.get(m, 0) + c1 * c2
            if v:
                r[m] = v
            else:
                r.pop(m, None)
    return r


def pvar(k):
    return {(k,): 1}


def pconst(c):
    return {(): c} if c else {}


def peval(p, point):
    s = 0
    for m, c in p.items():
        t = c
        for v in m:
            t = t * point[v]
        s += t
    return s


def psum(ps):
    r = {}
    for p in ps:
        r = padd(r, p)
    return r


def pm_zero(n, m):
    return [[{} for _ in range(m)] for _ in range(n)]


def pm_int(A):
    return [[pconst(x) for x in r] for r in A]


def pm_mul(A, B):
    k = len(B)
    m = len(B[0])
    out = pm_zero(len(A), m)
    for i in range(len(A)):
        for t in range(k):
            a = A[i][t]
            if not a:
                continue
            for j in range(m):
                b = B[t][j]
                if b:
                    out[i][j] = padd(out[i][j], pmul(a, b))
    return out


def pm_add(A, B):
    return [[padd(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def pm_sub(A, B):
    return [[psub(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def pm_scale(A, k):
    return [[pscale(x, k) for x in r] for r in A]


def pm_pscale(A, p):
    return [[pmul(x, p) for x in r] for r in A]


def pm_T(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def pm_diag(n, p):
    return [[dict(p) if i == j else {} for j in range(n)] for i in range(n)]


def pm_outer(u, v):
    return [[pmul(a, b) for b in v] for a in u]


def pm_eval(A, point):
    return [[peval(x, point) for x in r] for r in A]


def pblock(A, B, C, D):
    return [A[i] + B[i] for i in range(len(A))] + [C[i] + D[i] for i in range(len(C))]


def pm_vec(A, v):
    return [psum(pmul(A[i][j], v[j]) for j in range(len(v))) for i in range(len(A))]


def pdet(Mx):
    n = len(Mx)
    if n == 1:
        return Mx[0][0]
    total = {}
    sign = 1
    for j in range(n):
        if Mx[0][j]:
            minor = [row[:j] + row[j + 1:] for row in Mx[1:]]
            total = padd(total, pscale(pmul(Mx[0][j], pdet(minor)), sign))
        sign = -sign
    return total


def mono_space(varlist, deg):
    return set(itertools.combinations_with_replacement(sorted(varlist), deg))


def identity(Lm, Rm, varlist, deg):
    """Compares two polynomial matrices coefficient by coefficient over every monomial of total degree deg in varlist,
    and every monomial outside that space that either side holds. Returns ok, slots compared, mismatching coefficients."""
    if len(Lm) != len(Rm) or len(Lm[0]) != len(Rm[0]):
        raise SystemExit("shape")
    space = mono_space(varlist, deg)
    entries = 0
    extra = 0
    mism = 0
    for i in range(len(Lm)):
        for j in range(len(Lm[0])):
            entries += 1
            l = Lm[i][j]
            r = Rm[i][j]
            for m in set(l) | set(r):
                if m not in space:
                    extra += 1
            mism += len(psub(l, r))
    return mism == 0, entries * len(space) + extra, mism


def idrow(res):
    return [("coeffs_checked", res[1]), ("mismatches", res[2])]


# ---------------------------------------------------------------- Gaussian-integer polynomials: dict -> (re, im)

def gadd(p, q):
    r = dict(p)
    for m, (a, b) in q.items():
        x, y = r.get(m, (0, 0))
        x += a
        y += b
        if x or y:
            r[m] = (x, y)
        else:
            r.pop(m, None)
    return r


def gscale(p, z):
    a, b = z
    out = {}
    for m, (x, y) in p.items():
        u = x * a - y * b
        v = x * b + y * a
        if u or v:
            out[m] = (u, v)
    return out


def gmul(p, q):
    r = {}
    for m1, (a, b) in p.items():
        for m2, (c, d) in q.items():
            m = tuple(sorted(m1 + m2))
            x, y = r.get(m, (0, 0))
            x += a * c - b * d
            y += a * d + b * c
            if x or y:
                r[m] = (x, y)
            else:
                r.pop(m, None)
    return r


def gconjp(p):
    return {m: (a, -b) for m, (a, b) in p.items()}


def gderiv(p, v):
    out = {}
    for m, (a, b) in p.items():
        k = m.count(v)
        if k:
            lm = list(m)
            lm.remove(v)
            out = gadd(out, {tuple(lm): (k * a, k * b)})
    return out


def g_from_real(p):
    return {m: (c, 0) for m, c in p.items()}


def gm_zero(n=2, m=2):
    return [[{} for _ in range(m)] for _ in range(n)]


def gm_const(Z):
    return [[({(): z} if z != (0, 0) else {}) for z in r] for r in Z]


def gm_add(A, B):
    return [[gadd(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def gm_sub(A, B):
    return [[gadd(A[i][j], gscale(B[i][j], (-1, 0))) for j in range(len(A[0]))] for i in range(len(A))]


def gm_mul(A, B):
    k = len(B)
    m = len(B[0])
    out = gm_zero(len(A), m)
    for i in range(len(A)):
        for t in range(k):
            a = A[i][t]
            if not a:
                continue
            for j in range(m):
                b = B[t][j]
                if b:
                    out[i][j] = gadd(out[i][j], gmul(a, b))
    return out


def gm_comm(A, B):
    return gm_sub(gm_mul(A, B), gm_mul(B, A))


def gm_scalar(A, p):
    return [[gmul(x, p) for x in r] for r in A]


def gm_deriv(A, v):
    return [[gderiv(x, v) for x in r] for r in A]


def gm_iszero(A):
    return all(not x for r in A for x in r)


def gm_flat(A):
    return [x for r in A for x in r]


def g_identity(Ls, Rs):
    compared = 0
    mism = 0
    for l, r in zip(Ls, Rs):
        compared += len(set(l) | set(r))
        mism += len(gadd(l, gscale(r, (-1, 0))))
    return mism == 0, compared, mism


# ---------------------------------------------------------------- quaternions (4D.C1): e1 = 1, e2 = i, e3 = j, e4 = k, i j = k

QT = {}
for b0 in range(4):
    QT[(0, b0)] = (1, b0)
    QT[(b0, 0)] = (1, b0)
for a0 in (1, 2, 3):
    QT[(a0, a0)] = (-1, 0)
for (a0, b0, c0) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
    QT[(a0, b0)] = (1, c0)
    QT[(b0, a0)] = (-1, c0)

E_UNIT = [[1 if t == s else 0 for t in range(4)] for s in range(4)]


def qmul_i(x, y):
    out = [0, 0, 0, 0]
    for a in range(4):
        if x[a] == 0:
            continue
        for b in range(4):
            s, c = QT[(a, b)]
            out[c] += s * x[a] * y[b]
    return out


def qconj_i(x):
    return [x[0], -x[1], -x[2], -x[3]]


def qmul_p(x, y):
    out = [{}, {}, {}, {}]
    for a in range(4):
        if not x[a]:
            continue
        for b in range(4):
            if not y[b]:
                continue
            s, c = QT[(a, b)]
            out[c] = padd(out[c], pscale(pmul(x[a], y[b]), s))
    return out


def qconj_p(x):
    return [x[0], pscale(x[1], -1), pscale(x[2], -1), pscale(x[3], -1)]


def Lmat_i(x):
    out = izero(4, 4)
    for a in range(4):
        for b in range(4):
            s, c = QT[(a, b)]
            out[c][b] += s * x[a]
    return out


def Rmat_i(y):
    out = izero(4, 4)
    for a in range(4):
        for b in range(4):
            s, c = QT[(a, b)]
            out[c][a] += s * y[b]
    return out


def Lmat_p(x):
    out = pm_zero(4, 4)
    for a in range(4):
        for b in range(4):
            s, c = QT[(a, b)]
            out[c][b] = padd(out[c][b], pscale(x[a], s))
    return out


def Rmat_p(y):
    out = pm_zero(4, 4)
    for a in range(4):
        for b in range(4):
            s, c = QT[(a, b)]
            out[c][a] = padd(out[c][a], pscale(y[b], s))
    return out


def gamma_i(x):
    Z = izero(4, 4)
    return iblock(Z, Lmat_i(x), iscale(Lmat_i(qconj_i(x)), -1), Z)


def gamma_p(x):
    Z = pm_zero(4, 4)
    return pblock(Z, Lmat_p(x), pm_scale(Lmat_p(qconj_p(x)), -1), Z)


XI = [0, 1, 2, 3]
AV = [4, 5, 6, 7]
X = [pvar(t) for t in XI]
AQ = [pvar(t) for t in AV]


def nsq_of(v):
    return psum(pmul(c, c) for c in v)


NSQ = nsq_of(X)
NSQ_A = nsq_of(AQ)
GENERIC = {0: 2, 1: 3, 2: 5, 3: 7}
E1PT = {0: 1, 1: 0, 2: 0, 3: 0}

# ---------------------------------------------------------------- 2-forms on R^4 (C0.9, 4D.C4)

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PAIR_NAMES = ["e12", "e13", "e14", "e23", "e24", "e34"]


def perm_sign(seq):
    inv = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                inv += 1
    return -1 if inv & 1 else 1


STAR = izero(6, 6)
for col0, (i0, j0) in enumerate(PAIRS):
    rest0 = [t for t in range(4) if t not in (i0, j0)]
    STAR[PAIRS.index(tuple(rest0))][col0] = perm_sign([i0, j0] + rest0)

FB = [[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, -1, 0], [0, 0, 1, 1, 0, 0]]     # f1 = e12 + e34, f2 = e13 - e24, f3 = e14 + e23
GB = [[1, 0, 0, 0, 0, -1], [0, 1, 0, 0, 1, 0], [0, 0, 1, -1, 0, 0]]    # e12 - e34, e13 + e24, e14 - e23


def unit_p(j):
    return [pconst(1) if t == j else {} for t in range(4)]


def wedge_p(u, v):
    return [psub(pmul(u[i], v[j]), pmul(u[j], v[i])) for (i, j) in PAIRS]


def ccoords(w, basis):
    return [psum(pscale(w[I], basis[k][I]) for I in range(6)) for k in range(len(basis))]


def C_of(xq, basis=FB):
    cols = [ccoords(wedge_p(xq, unit_p(j)), basis) for j in range(4)]
    return [[cols[j][k] for j in range(4)] for k in range(len(basis))]


def M_of(xq, basis=FB):
    return [[pscale(xq[j], -1) for j in range(4)]] + C_of(xq, basis)


def colvec(v):
    return [[p] for p in v]


def main_measure():
    emit("HEAD", [("schema", "measure-spinor-symbol-4d-v1"), ("date", "2026-11-09"),
                  ("author", "the local Opus 5.5 seat, pid 8467a937cba309f7 (independent 4D measurement agent)"),
                  ("kind", "MEASUREMENT_independent_of_the_prediction_script"),
                  ("frozen_prediction", "freeze-spinor-symbol-20261109:predictions/spinor-symbol/PREDICTION-SPINOR-SYMBOL.hbp"),
                  ("frozen_prediction_sha256", FROZEN_SHA), ("frozen_addendum_sha256", ADDENDUM_SHA),
                  ("read_before_measuring", "HEAD, FORMAT and CONV rows only; no PRED, ADD or DIMS value"),
                  ("run", "python3 -I -B measure_4d.py"),
                  ("python", lst([sys.version_info[0], sys.version_info[1], sys.version_info[2]])),
                  ("float_used", 0), ("sampling_used", "controls only"),
                  ("index_note", "rows name basis vectors 1-based as in CONV (e1..e4 = 1, i, j, k); the code indexes them 0..3")])

    # ------------------------------------------------ Q: the quaternion table
    fails = 0
    for a in range(4):
        for b in range(4):
            for c in range(4):
                ea, eb, ec = E_UNIT[a], E_UNIT[b], E_UNIT[c]
                if qmul_i(qmul_i(ea, eb), ec) != qmul_i(ea, qmul_i(eb, ec)):
                    fails += 1
    one, qi, qj, qk = E_UNIT
    rules = (qmul_i(qi, qj) == qk and qmul_i(qj, qk) == qi and qmul_i(qk, qi) == qj
             and qmul_i(qi, qi) == [-1, 0, 0, 0] and qmul_i(qj, qj) == [-1, 0, 0, 0] and qmul_i(qk, qk) == [-1, 0, 0, 0])
    meas("M4D.Q1", "quaternion_table", fails == 0 and rules,
         [("table", "Hamilton, i j = k, j k = i, k i = j, i^2 = j^2 = k^2 = -1"), ("assoc_triples", 64), ("assoc_fail", fails),
          ("rules_hold", 1 if rules else 0)])
    XY = qmul_p(X, AQ)
    res = identity([[nsq_of(XY)]], [[pmul(NSQ, NSQ_A)]], XI + AV, 4)
    meas("M4D.Q2", "norm_multiplicative_nsq(xy)=nsq(x)nsq(y)", res[0], idrow(res))
    res1 = identity(pm_mul(pm_T(Lmat_p(X)), Lmat_p(X)), pm_diag(4, NSQ), XI, 2)
    res2 = identity(pm_mul(pm_T(Rmat_p(X)), Rmat_p(X)), pm_diag(4, NSQ), XI, 2)
    res3 = identity(pm_mul(Lmat_p(X), Rmat_p(AQ)), pm_mul(Rmat_p(AQ), Lmat_p(X)), XI + AV, 2)
    meas("M4D.Q3", "L_x^T L_x=nsq I, R_x^T R_x=nsq I, L_x R_y=R_y L_x", res1[0] and res2[0] and res3[0],
         [("coeffs_checked", res1[1] + res2[1] + res3[1]), ("mismatches", res1[2] + res2[2] + res3[2])])

    # ------------------------------------------------ G: gamma(xi) on S+ + S- = H + H (4D.C2)
    GX = gamma_p(X)
    res = identity(pm_mul(GX, GX), pm_diag(8, pscale(NSQ, -1)), XI, 2)
    meas("M4D.G1", "gamma(xi)^2+nsq(xi)I_8=0", res[0],
         idrow(res) + [("entries", 64), ("monomials_per_entry", 10), ("size", "8x8"),
                       ("gamma", "[[0, L_xi],[-L_conj(xi), 0]], S+ first")])
    gens = [gamma_i(E_UNIT[a]) for a in range(4)]
    ok_ac = 0
    for a in range(4):
        for b in range(4):
            lhs = iadd(imul(gens[a], gens[b]), imul(gens[b], gens[a]))
            if lhs == ieye(8, -2 if a == b else 0):
                ok_ac += 1
    meas("M4D.G2", "anticommutators_gamma_a_gamma_b+gamma_b_gamma_a=-2delta_ab", ok_ac == 16,
         [("ordered_pairs", 16), ("hold", ok_ac)])
    basis = [ieye(8)]
    frontier = [ieye(8)]
    cur = 1
    while frontier:
        new = []
        for B in frontier:
            for g in gens:
                P = imul(B, g)
                r = rank_int([flat(x) for x in basis] + [flat(P)])
                if r > cur:
                    basis.append(P)
                    new.append(P)
                    cur = r
        frontier = new
    alg_dim = cur
    subs = []
    for k in range(5):
        for S in itertools.combinations(range(4), k):
            P = ieye(8)
            for a in S:
                P = imul(P, gens[a])
            subs.append(P)
    r16 = rank_int([flat(x) for x in subs])
    meas("M4D.G3", "algebra_generated_by_gamma(e1..e4)", alg_dim == 16 and r16 == 16,
         [("dim", alg_dim), ("ordered_subset_products", len(subs)), ("rank_of_subset_products", r16),
          ("method", "closure of span(I) under right multiplication by the 4 generators, exact integer rank")])

    def commutant_dim(glist):
        rows = []
        for g in glist:
            for r in range(8):
                for c in range(8):
                    coeff = [0] * 64
                    for k in range(8):
                        coeff[8 * r + k] += g[k][c]
                        coeff[8 * k + c] -= g[r][k]
                    rows.append(coeff)
        return 64 - rank_int(rows)

    comm_dim = commutant_dim(gens)
    Z4 = izero(4, 4)
    RH = [iblock(Rmat_i(u), Z4, Z4, Rmat_i(u)) for u in E_UNIT]
    rh_comm = sum(1 for Rh in RH for g in gens if imul(Rh, g) == imul(g, Rh))
    rh_rank = rank_int([flat(x) for x in RH])
    meas("M4D.G4", "commutant_of_the_gamma_algebra_in_End_R(R^8)", comm_dim == 4 and rh_comm == 16 and rh_rank == 4,
         [("dim", comm_dim), ("right_H_multiplications_commuting", rh_comm), ("of", 16), ("right_H_rank", rh_rank),
          ("reading", "commutant = diag(R_u, R_u), u in H, right quaternion multiplication, dimension 4")])
    union_rank = rank_int([flat(x) for x in basis] + [flat(x) for x in RH])
    center = alg_dim + rh_rank - union_rank
    dbl = commutant_dim(RH[1:])
    alg_in = sum(1 for B in basis for Rh in RH if imul(B, Rh) == imul(Rh, B))
    meas("M4D.G5", "center_and_double_commutant", center == 1 and dbl == 16 and alg_in == 4 * len(basis),
         [("center_dim", center), ("commutant_of_right_H_dim", dbl), ("algebra_basis_commuting_with_right_H", alg_in),
          ("of", 4 * len(basis)),
          ("reading", "algebra = End_H(H^2), the commutant of right H, dimension 16, so Cl_4 = M_2(H) (iso type DERIVED from these dims)")])
    Gam = imul(imul(imul(gens[0], gens[1]), gens[2]), gens[3])
    diag = [Gam[i][i] for i in range(8)]
    offdiag = sum(1 for i in range(8) for j in range(8) if i != j and Gam[i][j] != 0)
    sq = imul(Gam, Gam) == ieye(8)
    anti = sum(1 for g in gens if iadd(imul(Gam, g), imul(g, Gam)) == izero(8, 8))
    target = iblock(ieye(4, -1), Z4, Z4, ieye(4))
    meas("M4D.G6", "volume_gamma1_gamma2_gamma3_gamma4", Gam == target and sq and anti == 4,
         [("diag", lst(diag)), ("offdiag_nonzero", offdiag),
          ("value", "diag(-I_4, I_4): -1 on S+, +1 on S-" if Gam == target else "not diag(-I_4, I_4), see diag"),
          ("square_is_I", 1 if sq else 0), ("anticommutes_with_generators", anti)])

    # ------------------------------------------------ S: the ASD symbol sequence (4D.C4, C5, C6)
    st_e12 = [STAR[r][0] for r in range(6)]
    st_e13 = [STAR[r][1] for r in range(6)]
    st_e14 = [STAR[r][2] for r in range(6)]
    star_ok = (st_e12 == [0, 0, 0, 0, 0, 1] and st_e13 == [0, 0, 0, 0, -1, 0] and st_e14 == [0, 0, 0, 1, 0, 0])
    sd = sum(1 for f in FB if matvec_i(STAR, f) == f)
    asd = sum(1 for g in GB if matvec_i(STAR, g) == [-x for x in g])
    gram_f = [[sum(FB[a][t] * FB[b][t] for t in range(6)) for b in range(3)] for a in range(3)]
    gram_g = [[sum(GB[a][t] * GB[b][t] for t in range(6)) for b in range(3)] for a in range(3)]
    cross = [[sum(FB[a][t] * GB[b][t] for t in range(6)) for b in range(3)] for a in range(3)]
    meas("M4D.S1", "hodge_star_and_lambda2_bases", star_ok and sd == 3 and asd == 3 and gram_f == ieye(3, 2)
         and gram_g == ieye(3, 2) and cross == izero(3, 3),
         [("star", "*e12 = e34, *e13 = -e24, *e14 = e23 from the permutation sign, vol = e1234"), ("star_rule_holds", 1 if star_ok else 0),
          ("f_self_dual", sd), ("asd_anti_self_dual", asd), ("gram_f", mat_str(gram_f)), ("gram_asd", mat_str(gram_g)),
          ("f_dot_asd", mat_str(cross))])
    W = wedge_p(X, AQ)
    cdef = ccoords(W, FB)
    x0, x1, x2, x3 = X
    a0, a1, a2, a3 = AQ

    def bil(p, q):
        return pmul(p, q)
    c1 = psum([bil(x0, a1), pscale(bil(x1, a0), -1), bil(x2, a3), pscale(bil(x3, a2), -1)])
    c2 = psum([bil(x0, a2), pscale(bil(x2, a0), -1), pscale(bil(x1, a3), -1), bil(x3, a1)])
    c3 = psum([bil(x0, a3), pscale(bil(x3, a0), -1), bil(x1, a2), pscale(bil(x2, a1), -1)])
    resA = identity(colvec(cdef), colvec([c1, c2, c3]), XI + AV, 2)
    lhs = [padd(W[I], psum(pscale(W[J], STAR[I][J]) for J in range(6))) for I in range(6)]
    rhs = [psum(pscale(cdef[k], FB[k][I]) for k in range(3)) for I in range(6)]
    resB = identity(colvec(lhs), colvec(rhs), XI + AV, 2)
    meas("M4D.S2", "c_k=<xi^a,f_k>_equals_4D.C5_formulas_and_(I+*)(xi^a)=sum_c_k_f_k", resA[0] and resB[0],
         [("coeffs_checked", resA[1] + resB[1]), ("mismatches", resA[2] + resB[2])])
    CX = C_of(X)
    s0 = colvec(X)
    res = identity(pm_mul(CX, s0), pm_zero(3, 1), XI, 2)
    meas("M4D.S3", "composite_C(xi)xi=0_Lambda0->Lambda1->Lambda2+", res[0], idrow(res))
    r_s0 = rank_int(pm_eval(s0, GENERIC))
    r_c = rank_int(pm_eval(CX, GENERIC))
    r_s0e = rank_int(pm_eval(s0, E1PT))
    r_ce = rank_int(pm_eval(CX, E1PT))
    meas("M4D.S4", "ranks_at_integer_xi", r_s0 == 1 and r_c == 3 and r_s0e == 1 and r_ce == 3,
         [("xi_generic", "2,3,5,7"), ("rank_sigma_d0", r_s0), ("rank_C", r_c), ("xi_e1_rank_sigma_d0", r_s0e),
          ("xi_e1_rank_C", r_ce), ("rank_sum", r_s0 + r_c), ("dim_Lambda1", 4)])
    resC = identity(pm_mul(CX, pm_T(CX)), pm_diag(3, NSQ), XI, 2)
    resD = identity(pm_add(pm_mul(pm_T(CX), CX), pm_outer(X, X)), pm_diag(4, NSQ), XI, 2)
    meas("M4D.S5", "C(xi)C(xi)^T=nsq I_3_and_C^T C+xi xi^T=nsq I_4", resC[0] and resD[0],
         [("CCT_coeffs_checked", resC[1]), ("CCT_mismatches", resC[2]), ("CTC_coeffs_checked", resD[1]),
          ("CTC_mismatches", resD[2]),
          ("reading", "rank C(xi) = 3 for every real xi != 0, and ker C(xi) = span(xi), so the sequence is exact for xi != 0")])
    d0 = len(list(itertools.combinations(range(4), 0)))
    d1 = len(list(itertools.combinations(range(4), 1)))
    d2p = rank_int(iadd(STAR, ieye(6)))
    meas("M4D.S6", "dims_and_euler_count_ASD_sequence", d0 - d1 + d2p == 0 and r_s0 + r_c == d1,
         [("dims", lst([d0, d1, d2p])), ("dim_Lambda2+_from", "rank(* + I) on Lambda2"), ("euler", d0 - d1 + d2p)])
    MX = M_of(X)
    resE = identity(pm_mul(pm_T(MX), MX), pm_diag(4, NSQ), XI, 2)
    resF = identity(pm_mul(MX, pm_T(MX)), pm_diag(4, NSQ), XI, 2)
    det = pdet(MX)
    nsq2 = pmul(NSQ, NSQ)
    dsign = 1 if not psub(det, nsq2) else (-1 if not padd(det, nsq2) else 0)
    meas("M4D.S7", "rolled_up_M(xi)^T M(xi)=nsq I_4_and_M M^T=nsq I_4", resE[0] and resF[0] and dsign != 0,
         [("MTM_coeffs_checked", resE[1]), ("MTM_mismatches", resE[2]), ("MMT_coeffs_checked", resF[1]),
          ("MMT_mismatches", resF[2]), ("det_M_over_nsq^2", dsign),
          ("M", "row 0 = -xi^T, rows 1..3 = c_k(xi, e_j)")])
    # R1: raw projection sigma = (-iota_xi, (xi wedge .)^+) in orthonormal bases; its Lambda2+ coordinates are c_k over sqrt2,
    # so 2 sigma^T sigma = 2 xi xi^T + C^T C with every coefficient an integer.
    S2 = pm_add(pm_scale(pm_outer(X, X), 2), pm_mul(pm_T(CX), CX))
    alpha = S2[0][0].get((1, 1), 0)
    beta = S2[0][1].get((0, 1), 0)
    resG = identity(S2, pm_add(pm_diag(4, pscale(NSQ, alpha)), pm_scale(pm_outer(X, X), beta)), XI, 2)
    cn, cd = reduced(alpha, 2)
    dn, dd = reduced(beta, 2)
    meas("M4D.S8", "R1_raw_projection_sigma^T sigma=c nsq I+(1-c)xi xi^T", resG[0] and alpha + beta == 2,
         idrow(resG) + [("cleared_form", "2 sigma^T sigma = " + str(alpha) + " nsq I + " + str(beta) + " xi xi^T"),
                        ("c_num", cn), ("c_den", cd), ("one_minus_c_num", dn), ("one_minus_c_den", dd)])

    # R2: 't Hooft symbols eta^a_(mu nu), a = 1..3, mu, nu = 1..4 with 4 the distinguished index, on the CONV labels e1..e4
    def eps3(a, b, c):
        if len({a, b, c}) < 3:
            return 0
        return perm_sign([a, b, c])
    ETA = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    for a in range(3):
        for m in range(3):
            for n in range(3):
                ETA[a][m][n] = eps3(a, m, n)
        ETA[a][a][3] = 1
        ETA[a][3][a] = -1
    eta_forms = [[ETA[a][i][j] for (i, j) in PAIRS] for a in range(3)]
    eta_sd = sum(1 for e in eta_forms if matvec_i(STAR, e) == e)
    Pm = [[exact_div(sum(eta_forms[a][t] * FB[k][t] for t in range(6)), 2) for k in range(3)] for a in range(3)]
    recon = all([sum(Pm[a][k] * FB[k][t] for k in range(3)) for t in range(6)] == eta_forms[a] for a in range(3))
    ETP = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    for a in range(3):
        for m in range(3):
            for n in range(3):
                ETP[a][1 + m][1 + n] = eps3(a, m, n)
        ETP[a][0][1 + a] = 1
        ETP[a][1 + a][0] = -1
    etp_forms = [[ETP[a][i][j] for (i, j) in PAIRS] for a in range(3)]
    etp_is_f = etp_forms == FB
    TX = [[psum(pscale(X[mu], ETA[a][mu][nu]) for mu in range(4)) for nu in range(4)] for a in range(3)]
    t_bil = [psum(pmul(TX[a][nu], AQ[nu]) for nu in range(4)) for a in range(3)]
    resH = identity(colvec(t_bil), colvec(ccoords(W, eta_forms)), XI + AV, 2)
    resI = identity(TX, pm_mul(pm_int(Pm), CX), XI, 1)
    resJ = identity(pm_mul(TX, pm_T(TX)), pm_diag(3, NSQ), XI, 2)
    lhs = [padd(W[I], psum(pscale(W[J], STAR[I][J]) for J in range(6))) for I in range(6)]
    rhs = [psum(pscale(t_bil[a], eta_forms[a][I]) for a in range(3)) for I in range(6)]
    resK = identity(colvec(lhs), colvec(rhs), XI + AV, 2)
    meas("M4D.S9", "R2_tHooft_eta_symbols", eta_sd == 3 and recon and resH[0] and resI[0] and resJ[0] and resK[0],
         [("eta", "eta^a_bc = eps_abc, eta^a_b4 = delta_ab, eta^a_4b = -delta_ab on e1..e4"), ("eta_self_dual", eta_sd),
          ("eta_in_f_basis", mat_str(Pm)), ("eta_reading", "eta^1 = f3, eta^2 = -f2, eta^3 = f1" if Pm == [[0, 0, 1], [0, -1, 0], [1, 0, 0]] else "see eta_in_f_basis"),
          ("eta_prime_distinguished_e1_equals_f", 1 if etp_is_f else 0),
          ("eta_prime", "eta'^a_1(a+1) = 1, eta'^a_(b+1)(c+1) = eps_abc: the same symbols with e1 = 1 as the distinguished index"),
          ("t_a", "t_a(xi, a) = eta^a_(mu nu) xi_mu a_nu = <xi wedge a, eta^a>"),
          ("TTT_nsq_I3", 1 if resJ[0] else 0), ("(I+*)(xi^a)=sum_t_a_eta^a", 1 if resK[0] else 0),
          ("coeffs_checked", resH[1] + resI[1] + resJ[1] + resK[1]), ("mismatches", resH[2] + resI[2] + resJ[2] + resK[2])])
    sig_t = [[pscale(X[j], -1) for j in range(4)]] + TX
    resL = identity(pm_mul(pm_T(sig_t), sig_t), pm_diag(4, NSQ), XI, 2)
    meas("M4D.S10", "R2.t_coefficients_t_a_as_coordinates", resL[0],
         idrow(resL) + [("identity", "sigma_t^T sigma_t = nsq I_4, sigma_t = (-xi^T; T(xi))"), ("c_num", 1), ("c_den", 1)])
    S4 = pm_add(pm_scale(pm_outer(X, X), 4), pm_mul(pm_T(TX), TX))
    al4 = S4[0][0].get((1, 1), 0)
    be4 = S4[0][1].get((0, 1), 0)
    resM = identity(S4, pm_add(pm_diag(4, pscale(NSQ, al4)), pm_scale(pm_outer(X, X), be4)), XI, 2)
    c4n, c4d = reduced(al4, 4)
    meas("M4D.S11", "R2.s_flat_coefficients_s_a=t_a_over_2_read_as_orthonormal", resM[0] and al4 + be4 == 4,
         idrow(resM) + [("cleared_form", "4 sigma_s^T sigma_s = " + str(al4) + " nsq I + " + str(be4) + " xi xi^T"),
                        ("c_num", c4n), ("c_den", c4d),
                        ("reading", "omega^+ = sum_a s_a eta^a with s_a = t_a over 2; with the Gram 2I of eta kept, c = 1/2 as in R1")])
    S2t = pm_add(pm_scale(pm_outer(X, X), 2), pm_mul(pm_T(TX), TX))
    resN = identity(S2t, pm_add(pm_diag(4, NSQ), pm_outer(X, X)), XI, 2)
    meas("M4D.S12", "R2.s_gram_coefficients_s_a_with_Gram_2I", resN[0],
         idrow(resN) + [("cleared_form", "2 sigma_s^T G sigma_s = 1 nsq I + 1 xi xi^T"), ("c_num", 1), ("c_den", 2)])

    # identification of M(xi) a with a quaternion product
    Ma = pm_vec(MX, AQ)
    matches = []
    tried = 0
    for order in ("xa", "ax"):
        for cx in (0, 1):
            for ca in (0, 1):
                xx = qconj_p(X) if cx else X
                yy = qconj_p(AQ) if ca else AQ
                prod = qmul_p(xx, yy) if order == "xa" else qmul_p(yy, xx)
                nx = "conj(xi)" if cx else "xi"
                na = "conj(a)" if ca else "a"
                for sgn in (1, -1):
                    tried += 1
                    name = ("-" if sgn < 0 else "") + (nx + "*" + na if order == "xa" else na + "*" + nx)
                    if identity(colvec(Ma), colvec([pscale(p, sgn) for p in prod]), XI + AV, 2)[0]:
                        matches.append(name)
    eqv = identity(colvec(qmul_p(X, qconj_p(AQ))), colvec(qconj_p(qmul_p(AQ, qconj_p(X)))), XI + AV, 2)[0]
    meas("M4D.S13", "identification_M(xi)a_as_quaternion_product", len(matches) == 1,
         [("candidates", tried), ("matching", lst(matches) if matches else "none"),
          ("xi*conj(a)=conj(a*conj(xi))", 1 if eqv else 0),
          ("coordinates", "M(xi) a read as the quaternion (row0) 1 + (row1) i + (row2) j + (row3) k")])

    # C0.7: T = * on Lambda2, k = 2
    s2 = imul(STAR, STAR) == ieye(6)
    TpI = iadd(STAR, ieye(6))
    TmI = isub(STAR, ieye(6))
    rk_p = rank_int(TpI)
    rk_m = rank_int(TmI)
    idem = imul(TpI, TpI) == iscale(TpI, 2)
    in_img = rank_int(iT(TpI) + FB) == rk_p
    meas("M4D.S14", "splitting_operator_T=*_k=2", s2 and rk_p == 3 and rk_m == 3 and idem and in_img,
         [("T_squared_I", 1 if s2 else 0), ("rank_T_plus_I", rk_p), ("rank_T_minus_I", rk_m),
          ("eigen_plus1_dim", 6 - rk_m), ("eigen_minus1_dim", 6 - rk_p), ("(T+I)^2=2(T+I)", 1 if idem else 0),
          ("f_basis_spans_image_of_T+I", 1 if in_img else 0)])

    # C0.5: [c(omega), gamma(v)] = 2 gamma(R_omega v)
    OMV = [pvar(10 + t) for t in range(6)]
    VV = [pvar(20 + t) for t in range(4)]
    cOm = pm_zero(8, 8)
    for I, (i, j) in enumerate(PAIRS):
        cOm = pm_add(cOm, pm_pscale(pm_int(imul(gens[i], gens[j])), OMV[I]))
    gv = gamma_p(VV)
    Rv = [{}, {}, {}, {}]
    for I, (i, j) in enumerate(PAIRS):
        Rv[j] = padd(Rv[j], pmul(OMV[I], VV[i]))
        Rv[i] = psub(Rv[i], pmul(OMV[I], VV[j]))
    res = identity(pm_sub(pm_mul(cOm, gv), pm_mul(gv, cOm)), pm_scale(gamma_p(Rv), 2), list(range(10, 16)) + list(range(20, 24)), 2)
    meas("M4D.S15", "C0.5_[c(omega),gamma(v)]=2gamma(R_omega v)", res[0], idrow(res))

    # C0.6: Psi(omega) = c(omega) eta, eta = 1 in S+
    PsiM = [[imul(gens[i], gens[j])[r][0] for (i, j) in PAIRS] for r in range(8)]
    rk_psi = rank_int(PsiM)
    on_asd = sum(1 for g in GB if matvec_i(PsiM, g) == [0] * 8)
    on_sd = sum(1 for f in FB if matvec_i(PsiM, f) != [0] * 8)
    rk_sd = rank_int(iT([matvec_i(PsiM, f) for f in FB]))
    s_minus_zero = all(PsiM[r] == [0] * 6 for r in range(4, 8))
    eta_comp_zero = PsiM[0] == [0] * 6
    meas("M4D.S16", "C0.6_Psi_kernel_and_image_4D", rk_psi == 3 and on_asd == 3 and on_sd == 3 and rk_sd == 3
         and s_minus_zero and eta_comp_zero,
         [("rank", rk_psi), ("kernel_dim", 6 - rk_psi), ("kernel", "Lambda2- (all three anti-self-dual basis forms map to 0)"),
          ("asd_in_kernel", on_asd), ("f_nonzero", on_sd), ("rank_on_Lambda2+", rk_sd),
          ("image", "Im H inside S+ (S- part 0, eta component 0)")])

    # ------------------------------------------------ E: Sp(1) x Sp(1) equivariance on 2T x 2T (4D.C3)
    T2 = []
    for s in (2, -2):
        for t in range(4):
            v = [0, 0, 0, 0]
            v[t] = s
            T2.append(v)
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2_ in (1, -1):
                for s3 in (1, -1):
                    T2.append([s0, s1, s2_, s3])
    T2set = set(tuple(v) for v in T2)
    closed = 0
    for P in T2:
        for Q in T2:
            R = [exact_div(x, 2) for x in qmul_i(P, Q)]
            if tuple(R) in T2set:
                closed += 1
    gen_set = {(0, 2, 0, 0), (1, 1, 1, 1)}
    frontier2 = list(gen_set)
    while frontier2:
        nxt = []
        for P in frontier2:
            for Q in list(gen_set):
                for R in (qmul_i(list(P), list(Q)), qmul_i(list(Q), list(P))):
                    R2 = tuple(exact_div(x, 2) for x in R)
                    if R2 not in gen_set:
                        gen_set.add(R2)
                        nxt.append(R2)
        frontier2 = nxt
    norms = all(sum(x * x for x in P) == 4 for P in T2)
    meas("M4D.E0", "binary_tetrahedral_2T_scaled_by_2", len(T2) == 24 and len(T2set) == 24 and closed == 576 and norms
         and gen_set == T2set,
         [("elements", len(T2set)), ("closed_products", closed), ("of", 576), ("nsq_each", 4),
          ("generated_by_2i_and_1+i+j+k", 1 if gen_set == T2set else 0)])
    ok_g = 0
    slots_g = 0
    mism_g = 0
    ok_m = 0
    slots_m = 0
    mism_m = 0
    Kby_P = {}
    Kset = set()
    blockform = 0
    kkt = 0
    k_is_ad = 0
    swap_fail = 0
    MXe1 = pm_eval(MX, E1PT)
    for P in T2:
        Pp = [pconst(x) for x in P]
        for Q in T2:
            Qb = qconj_i(Q)
            Qbp = [pconst(x) for x in Qb]
            Y = qmul_p(qmul_p(Pp, X), Qbp)
            Gm = iblock(Lmat_i(P), Z4, Z4, Lmat_i(Q))
            res = identity(pm_mul(gamma_p(Y), pm_int(Gm)), pm_scale(pm_mul(pm_int(Gm), GX), 4), XI, 1)
            slots_g += res[1]
            mism_g += res[2]
            if res[0]:
                ok_g += 1
            Gs = iblock(Lmat_i(Q), Z4, Z4, Lmat_i(P))
            if not identity(pm_mul(gamma_p(Y), pm_int(Gs)), pm_scale(pm_mul(pm_int(Gs), GX), 4), XI, 1)[0]:
                swap_fail += 1
            MY = M_of(Y)
            Rho = imul(Lmat_i(P), Rmat_i(Qb))
            K = imul(imul(pm_eval(MY, E1PT), Rho), iT(MXe1))
            res = identity(pm_mul(MY, pm_int(Rho)), pm_mul(pm_int(K), MX), XI, 1)
            slots_m += res[1]
            mism_m += res[2]
            if res[0]:
                ok_m += 1
            Kt = tuple(flat(K))
            Kset.add(Kt)
            Kby_P.setdefault(tuple(P), set()).add(Kt)
            if K[0][0] == 16 and K[0][1:] == [0, 0, 0] and [K[r][0] for r in (1, 2, 3)] == [0, 0, 0]:
                blockform += 1
            if imul(K, iT(K)) == ieye(4, 256):
                kkt += 1
            if K == iscale(imul(Lmat_i(P), Rmat_i(qconj_i(P))), 4):
                k_is_ad += 1
    only_p = all(len(s) == 1 for s in Kby_P.values())
    meas("M4D.E1", "gamma(p xi conj(q))g=g gamma(xi)_on_2Tx2T", ok_g == 576,
         [("pairs", 576), ("hold", ok_g), ("coeffs_checked", slots_g), ("mismatches", mism_g),
          ("integer_form", "gamma(P xi conj(Q)) G = 4 G gamma(xi), G = diag(L_P, L_Q), P = 2p, Q = 2q")])
    meas("M4D.E2", "M_equivariance_on_2Tx2T", ok_m == 576 and blockform == 576 and kkt == 576 and only_p
         and len(Kset) == 12 and k_is_ad == 576,
         [("pairs", 576), ("hold", ok_m), ("coeffs_checked", slots_m), ("mismatches", mism_m),
          ("integer_form", "M(P xi conj(Q)) (L_P R_conj(Q)) = K M(xi), K = M(P conj(Q)) L_P R_conj(Q) M(e1)^T"),
          ("K_block_16_plus_3x3", blockform), ("KKT_256I", kkt), ("K_depends_only_on_P", 1 if only_p else 0),
          ("distinct_K", len(Kset)), ("K_equals_4_Ad(P)", k_is_ad),
          ("reading", "Lambda0 trivial; Lambda2+ carried by SO(3) of the factor p, the one acting on S+")])
    inf_ok = 0
    inf_slots = 0
    for side in (0, 1):
        for u in (qi, qj, qk):
            Lu = Lmat_i(u)
            if side == 0:
                D = iblock(Lu, Z4, Z4, Z4)
                dxi = qmul_p([pconst(x) for x in u], X)
            else:
                D = iblock(Z4, Z4, Z4, Lu)
                dxi = [pscale(c, -1) for c in qmul_p(X, [pconst(x) for x in u])]
            Dp = pm_int(D)
            res = identity(pm_sub(pm_mul(Dp, GX), pm_mul(GX, Dp)), gamma_p(dxi), XI, 1)
            inf_slots += res[1]
            if res[0]:
                inf_ok += 1
    meas("M4D.E3", "infinitesimal_sp(1)+sp(1)_[diag(L_u,L_v),gamma(xi)]=gamma(u xi-xi v)", inf_ok == 6,
         [("generators", 6), ("hold", inf_ok), ("coeffs_checked", inf_slots)])
    control("K4D.E4", "swapped_g=diag(L_q,L_p)_equivariance", swap_fail > 0,
            [("pairs_failing", swap_fail), ("of", 576)])
    kq = _k_by_q(T2, MX, MXe1)
    control("K4D.E5", "claim_K_depends_only_on_Q", not all(len(s) == 1 for s in kq.values()),
            [("distinct_K_for_fixed_Q_max", max(len(s) for s in kq.values())),
             ("distinct_K_for_fixed_P_max", max(len(s) for s in Kby_P.values()))])

    # ------------------------------------------------ OB: the obstruction (OB.C1, OB.C2)
    sig1 = [[(0, 0), (1, 0)], [(1, 0), (0, 0)]]
    sig2 = [[(0, 0), (0, -1)], [(0, 1), (0, 0)]]
    sig3 = [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]
    sigs = [gm_const(s) for s in (sig1, sig2, sig3)]
    Tb = [gm_scalar(s, {(): (0, 1)}) for s in sigs]
    pauli_ok = 0
    for a in range(3):
        for b in range(3):
            lhs = gm_comm(sigs[a], sigs[b])
            rhs = gm_zero()
            for c in range(3):
                e = eps3(a, b, c)
                if e:
                    rhs = gm_add(rhs, gm_scalar(sigs[c], {(): (0, 2 * e)}))
            if gm_iszero(gm_sub(lhs, rhs)):
                pauli_ok += 1
    A_na = [Tb[0], Tb[1], gm_zero(), gm_zero()]
    F12 = gm_comm(A_na[0], A_na[1])
    F12_val = [[F12[i][j].get((), (0, 0)) for j in range(2)] for i in range(2)]

    def su2_coords(Xm):
        y1 = {m: b for m, (a, b) in Xm[0][1].items() if b}
        y2 = {m: a for m, (a, b) in Xm[0][1].items() if a}
        y3 = {m: b for m, (a, b) in Xm[0][0].items() if b}
        R = gm_add(gm_add(gm_scalar(Tb[0], g_from_real(y1)), gm_scalar(Tb[1], g_from_real(y2))), gm_scalar(Tb[2], g_from_real(y3)))
        return [y1, y2, y3], gm_iszero(gm_sub(Xm, R))

    fco, fin = su2_coords(F12)
    fco_val = [p.get((), 0) for p in fco]
    meas("M4D.OB1", "A1=i sigma1_A2=i sigma2_F12=[A1,A2]", pauli_ok == 9 and fin and fco_val == [0, 0, -2],
         [("pauli_commutators_hold", pauli_ok), ("of", 9),
          ("F12", "[[" + lst(F12_val[0][0]) + ", " + lst(F12_val[0][1]) + "],[" + lst(F12_val[1][0]) + ", " + lst(F12_val[1][1]) + "]] as re:im pairs"),
          ("F12_su2_coords_in_i_sigma_basis", lst(fco_val)), ("F12_reading", "F12 = -2 i sigma3, nonzero")])

    def dA0(f, A):
        return [gm_add(gm_deriv(f, i), gm_comm(A[i], f)) for i in range(4)]

    def dA1(w, A):
        out = {}
        for (i, j) in PAIRS:
            out[(i, j)] = gm_sub(gm_add(gm_deriv(w[j], i), gm_comm(A[i], w[j])), gm_add(gm_deriv(w[i], j), gm_comm(A[j], w[i])))
        return out

    def curv(A):
        return {(i, j): gm_add(gm_sub(gm_deriv(A[j], i), gm_deriv(A[i], j)), gm_comm(A[i], A[j])) for (i, j) in PAIRS}

    MONO2 = [()] + [(t,) for t in range(4)] + list(itertools.combinations_with_replacement(range(4), 2))
    fgen = gm_zero()
    for b in range(3):
        gb = {}
        for idx, m in enumerate(MONO2):
            gb[tuple(sorted(m + (100 + 15 * b + idx,)))] = (1, 0)
        fgen = gm_add(fgen, gm_scalar(Tb[b], gb))

    def ddf_check(A):
        d1f = dA0(fgen, A)
        d1f_terms = sum(len(x) for w in d1f for x in gm_flat(w))
        raw_terms = 0
        for (i, j) in PAIRS:
            for part in (gm_deriv(d1f[j], i), gm_comm(A[i], d1f[j]), gm_deriv(d1f[i], j), gm_comm(A[j], d1f[i])):
                raw_terms += sum(len(x) for x in gm_flat(part))
        DF_TERMS.append((d1f_terms, raw_terms))
        lhs = dA1(d1f, A)
        F = curv(A)
        Ls = []
        Rs = []
        nonzero_terms = 0
        nonzero_comp = []
        for I in PAIRS:
            Ls += gm_flat(lhs[I])
            Rs += gm_flat(gm_comm(F[I], fgen))
            t = sum(len(x) for x in gm_flat(lhs[I]))
            nonzero_terms += t
            if t:
                nonzero_comp.append(PAIR_NAMES[PAIRS.index(I)])
        ok, compared, mism = g_identity(Ls, Rs)
        return ok, compared, mism, nonzero_terms, nonzero_comp, lhs, F

    DF_TERMS = []
    ok, compared, mism, nz, nzc, lhs_na, F_na = ddf_check(A_na)
    meas("M4D.OB2", "d_A d_A f=[F_A,f]_for_generic_su2_valued_f", ok and nz > 0 and nzc == ["e12"],
         [("f", "sum_b g_b(x) i sigma_b, g_b generic of degree at most 2 (45 symbolic coefficients)"),
          ("gaussian_coeffs_compared", compared), ("mismatches", mism), ("d_A_d_A_f_nonzero_terms", nz),
          ("d_A_f_terms", DF_TERMS[0][0]), ("terms_before_cancellation_in_d_A_d_A_f", DF_TERMS[0][1]),
          ("nonzero_components", lst(nzc)), ("reading", "d_A d_A = [F_A, .] != 0: the twisted sequence is no complex")])
    control("K4D.OB3", "claim_d_A_d_A_f=0_for_A1=i sigma1_A2=i sigma2", nz > 0, [("nonzero_terms", nz)])
    okab, cab, mab, nzab, _c, _l, _F = ddf_check([Tb[2], Tb[2], gm_zero(), gm_zero()])
    ok0, c0_, m0, nz0, _c0, _l0, _F0 = ddf_check([gm_zero(), gm_zero(), gm_zero(), gm_zero()])
    meas("M4D.OB4", "abelian_A1=A2=i sigma3_and_A=0_give_d_A d_A f=0", okab and nzab == 0 and ok0 and nz0 == 0,
         [("abelian_nonzero_terms", nzab), ("zero_A_nonzero_terms", nz0),
          ("abelian_d_A_f_terms", DF_TERMS[1][0]), ("abelian_terms_before_cancellation", DF_TERMS[1][1]),
          ("zero_A_d_f_terms", DF_TERMS[2][0]), ("zero_A_terms_before_cancellation", DF_TERMS[2][1])])
    cF = [gm_zero() for _ in range(3)]
    cL = [gm_zero() for _ in range(3)]
    for k in range(3):
        for Ii, I in enumerate(PAIRS):
            if FB[k][Ii]:
                cF[k] = gm_add(cF[k], gm_scalar(F_na[I], {(): (FB[k][Ii], 0)}))
                cL[k] = gm_add(cL[k], gm_scalar(lhs_na[I], {(): (FB[k][Ii], 0)}))
    Ls = []
    Rs = []
    for k in range(3):
        Ls += gm_flat(cL[k])
        Rs += gm_flat(gm_comm(cF[k], fgen))
    okp, cp, mp = g_identity(Ls, Rs)
    cF_coords = [lst([p.get((), 0) for p in su2_coords(cF[k])[0]]) for k in range(3)]
    zeroth = []
    for b in range(3):
        col = []
        for k in range(3):
            cc, _ok = su2_coords(gm_comm(cF[k], Tb[b]))
            col += [p.get((), 0) for p in cc]
        zeroth.append(col)
    rk0 = rank_int(zeroth)
    nzp = sum(len(x) for k in range(3) for x in gm_flat(cL[k]))
    meas("M4D.OB5", "self_dual_part_d_A^+ d_A f=[F_A^+,f]", okp and nzp > 0 and rk0 == 2,
         [("F_plus_c_coords_su2", "; ".join(cF_coords)), ("gaussian_coeffs_compared", cp), ("mismatches", mp),
          ("d_A^+_d_A_f_nonzero_terms", nzp), ("rank_of_f->[F^+,f]_on_su2", rk0),
          ("reading", "F_A^+ != 0 for this constant A, so d_A^+ d_A != 0 and the ASD sequence with this A is no complex")])

    ell = {(i, 4 + i): (1, 0) for i in range(4)}
    XIv = [4, 5, 6, 7]

    def symbols(A):
        S0 = [[{} for _ in range(3)] for _ in range(12)]
        xfree = True
        for b in range(3):
            lf = gm_scalar(Tb[b], ell)
            d1 = dA0(lf, A)
            d2 = dA0(Tb[b], A)
            for i in range(4):
                cc, okc = su2_coords(gm_sub(d1[i], gm_scalar(d2[i], ell)))
                xfree = xfree and okc
                for bb in range(3):
                    S0[3 * i + bb][b] = cc[bb]
        S1 = [[{} for _ in range(12)] for _ in range(9)]
        for j in range(4):
            for b in range(3):
                w = [gm_zero() for _ in range(4)]
                w[j] = Tb[b]
                lw = [gm_scalar(x, ell) for x in w]
                o1 = dA1(lw, A)
                o2 = dA1(w, A)
                Om = {I: gm_sub(o1[I], gm_scalar(o2[I], ell)) for I in PAIRS}
                for k in range(3):
                    ck = gm_zero()
                    for Ii, I in enumerate(PAIRS):
                        if FB[k][Ii]:
                            ck = gm_add(ck, gm_scalar(Om[I], {(): (FB[k][Ii], 0)}))
                    cc, okc = su2_coords(ck)
                    xfree = xfree and okc
                    for bb in range(3):
                        S1[3 * k + bb][3 * j + b] = cc[bb]
        for Mx in (S0, S1):
            for r in Mx:
                for p in r:
                    for m in p:
                        if any(v not in XIv for v in m):
                            xfree = False
        return S0, S1, xfree

    S0, S1, xfree = symbols(A_na)
    S0z, S1z, _xf = symbols([gm_zero(), gm_zero(), gm_zero(), gm_zero()])
    a_indep = (S0 == S0z and S1 == S1z)
    NSQX = nsq_of([pvar(t) for t in XIv])
    r1 = identity(pm_mul(S1, S0), pm_zero(9, 3), XIv, 2)
    r2 = identity(pm_mul(pm_T(S0), S0), pm_diag(3, NSQX), XIv, 2)
    r3 = identity(pm_mul(S1, pm_T(S1)), pm_diag(9, NSQX), XIv, 2)
    r4 = identity(pm_add(pm_mul(pm_T(S1), S1), pm_mul(S0, pm_T(S0))), pm_diag(12, NSQX), XIv, 2)
    gpt = {4: 2, 5: 3, 6: 5, 7: 7}
    rk_a = rank_int(pm_eval(S0, gpt))
    rk_b = rank_int(pm_eval(S1, gpt))
    meas("M4D.OB6", "su2_valued_symbol_sequence_3_12_9_exact", xfree and a_indep and r1[0] and r2[0] and r3[0] and r4[0]
         and rk_a == 3 and rk_b == 9,
         [("method", "sigma(D)(xi) u = D(l u) - l D(u), l = sum xi_i x_i, factor i dropped (C0.2), with A = (i sigma1, i sigma2, 0, 0)"),
          ("dims", "3,12,9"), ("euler", 3 - 12 + 9), ("composite_zero", 1 if r1[0] else 0),
          ("S0TS0_nsq_I3", 1 if r2[0] else 0), ("S1S1T_nsq_I9", 1 if r3[0] else 0), ("S1TS1+S0S0T_nsq_I12", 1 if r4[0] else 0),
          ("rank_at_2357", lst([rk_a, rk_b])), ("symbol_independent_of_A", 1 if a_indep else 0), ("x_free", 1 if xfree else 0),
          ("coeffs_checked", r1[1] + r2[1] + r3[1] + r4[1]), ("mismatches", r1[2] + r2[2] + r3[2] + r4[2])])

    # ------------------------------------------------ TS: C^4 (x) C^4 = End(C^4), Spin(4) x Spin(4)
    def zmul(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def zadd(a, b):
        return (a[0] + b[0], a[1] + b[1])

    def zm_mul(A, B):
        n = len(A)
        k = len(B)
        m = len(B[0])
        out = [[(0, 0)] * m for _ in range(n)]
        for i in range(n):
            for j in range(m):
                s = (0, 0)
                for t in range(k):
                    s = zadd(s, zmul(A[i][t], B[t][j]))
                out[i][j] = s
        return out

    def zm_dag(A):
        return [[(A[j][i][0], -A[j][i][1]) for j in range(len(A))] for i in range(len(A[0]))]

    def zm_conj(A):
        return [[(x[0], -x[1]) for x in r] for r in A]

    def zm_scale(A, k):
        return [[(k * x[0], k * x[1]) for x in r] for r in A]

    def zm_eye(n, k=1):
        return [[(k, 0) if i == j else (0, 0) for j in range(n)] for i in range(n)]

    def U_i(P):
        a, b, c, d = P
        return [[(a, b), (c, d)], [(-c, d), (a, -b)]]

    def U_p(xq):
        def z(re, im):
            return gadd(g_from_real(re), gscale(g_from_real(im), (0, 1)))
        return [[z(xq[0], xq[1]), z(xq[2], xq[3])], [z(pscale(xq[2], -1), xq[3]), z(xq[0], pscale(xq[1], -1))]]

    def spin4(P, Q):
        Up = U_i(P)
        Uq = U_i(Q)
        out = zm_eye(4, 0)
        for i in range(2):
            for j in range(2):
                out[i][j] = Up[i][j]
                out[2 + i][2 + j] = Uq[i][j]
        return out

    def realify(A):
        return [x[0] for r in A for x in r] + [x[1] for r in A for x in r]

    def times_i(A):
        return [[(-x[1], x[0]) for x in r] for r in A]

    def Eab(a, b):
        out = zm_eye(4, 0)
        out[a][b] = (1, 0)
        return out

    Ebasis = [Eab(a, b) for a in range(4) for b in range(4)]
    rk_real = rank_int([realify(E) for E in Ebasis] + [realify(times_i(E)) for E in Ebasis])
    rU1 = g_identity(gm_flat(gm_mul(U_p(X), U_p(AQ))), gm_flat(U_p(XY)))
    UXd = [[gconjp(U_p(X)[j][i]) for j in range(2)] for i in range(2)]
    rU2 = g_identity(gm_flat(gm_mul(U_p(X), UXd)), gm_flat([[g_from_real(NSQ), {}], [{}, g_from_real(NSQ)]]))
    meas("M4D.TS1", "U:H->M_2(C)_homomorphism_and_unitarity", rU1[0] and rU2[0],
         [("U", "a + b i + c j + d k -> [[a + b i, c + d i],[-c + d i, a - b i]]"), ("U(x)U(y)=U(xy)", 1 if rU1[0] else 0),
          ("U(x)U(x)^dag=nsq I", 1 if rU2[0] else 0), ("gaussian_coeffs_compared", rU1[1] + rU2[1]),
          ("mismatches", rU1[2] + rU2[2])])

    def gamma_c(xq):
        Z2 = gm_zero()
        Ux = U_p(xq)
        Ub = [[gscale(x, (-1, 0)) for x in r] for r in U_p(qconj_p(xq))]
        top = [Z2[i] + Ux[i] for i in range(2)]
        bot = [Ub[i] + Z2[i] for i in range(2)]
        return top + bot

    GC = gamma_c(X)
    rgc = g_identity(gm_flat(gm_mul(GC, GC)), gm_flat([[g_from_real(pscale(NSQ, -1)) if i == j else {} for j in range(4)] for i in range(4)]))
    gcs = []
    for a in range(4):
        pt = {t: (1 if t == a else 0) for t in range(4)}
        gcs.append([[geval(x, pt) for x in r] for r in GC])
    prods = []
    for k in range(5):
        for S in itertools.combinations(range(4), k):
            P = zm_eye(4)
            for a in S:
                P = zm_mul(P, gcs[a])
            prods.append(P)
    rk_cl = rank_int([realify(P) for P in prods] + [realify(times_i(P)) for P in prods])
    meas("M4D.TS2", "C^4(x)C^4=End(C^4)_dimension_16", rk_real == 32 and rgc[0] and rk_cl == 32,
         [("map", "u (x) v -> u v^T"), ("complex_rank_of_images_of_basis_tensors", exact_div(rk_real, 2)),
          ("complex_gamma_c(xi)^2=-nsq I_4", 1 if rgc[0] else 0), ("complex_rank_of_16_gamma_products", exact_div(rk_cl, 2)),
          ("reading", "dim_C End(C^4) = 16 = 4 x 4, and the complexified Clifford algebra fills it")])
    GEN4 = [((0, 2, 0, 0), (2, 0, 0, 0)), ((1, 1, 1, 1), (2, 0, 0, 0)), ((2, 0, 0, 0), (0, 2, 0, 0)), ((2, 0, 0, 0), (1, 1, 1, 1))]
    GM = [spin4(P, Q) for (P, Q) in GEN4]

    def act(G, H, Am):
        return zm_mul(zm_mul(G, Am), zm_dag(H))

    inter_ok = 0
    inter_n = 0
    for G in GM:
        for H in GM:
            Hc = zm_conj(H)
            for a in range(4):
                for b in range(4):
                    inter_n += 1
                    u = [G[i][a] for i in range(4)]
                    v = [Hc[j][b] for j in range(4)]
                    outer = [[zmul(u[i], v[j]) for j in range(4)] for i in range(4)]
                    if outer == act(G, H, Ebasis[4 * a + b]):
                        inter_ok += 1
    hom_ok = 0
    hom_n = 0
    for G1 in GM:
        for H1 in GM:
            for G2 in GM:
                for H2 in GM:
                    G12 = zm_mul(G1, G2)
                    H12 = zm_mul(H1, H2)
                    for E in Ebasis:
                        hom_n += 1
                        if act(G12, H12, E) == act(G1, H1, act(G2, H2, E)):
                            hom_ok += 1
    meas("M4D.TS3", "Spin(4)xSpin(4)_action_A->gAh^-1_on_generators", inter_ok == inter_n and hom_ok == hom_n,
         [("generators", "(2i,2), (1+i+j+k,2), (2,2i), (2,1+i+j+k) scaled by 2, g = diag(U(p), U(q))"),
          ("intertwining_u(x)v->uv^T", lst([inter_ok, inter_n])), ("homomorphism_checks", lst([hom_ok, hom_n])),
          ("integer_form", "G A H^dag = 4 g A h^-1 for G = 2g, H = 2h")])
    m2 = zm_eye(4, -2)
    p2 = zm_eye(4, 2)
    triv = all(act(m2, m2, E) == zm_scale(E, 4) for E in Ebasis)
    neg1 = all(act(m2, p2, E) == zm_scale(E, -4) for E in Ebasis)
    neg2 = all(act(p2, m2, E) == zm_scale(E, -4) for E in Ebasis)
    meas("M4D.TS4", "(-1,-1)_acts_trivially_(-1,1)_and_(1,-1)_act_as_-1", triv and neg1 and neg2,
         [("(-1,-1)_trivial", 1 if triv else 0), ("(-1,1)_is_-1", 1 if neg1 else 0), ("(1,-1)_is_-1", 1 if neg2 else 0)])
    control("K4D.TS5", "claim_(-1,1)_acts_trivially", not all(act(m2, p2, E) == zm_scale(E, 4) for E in Ebasis), [])
    Ulist = [U_i(P) for P in T2]
    Bt = [[zm_mul(Ulist[i1], zm_dag(Ulist[i2])) == zm_eye(2, 4) for i2 in range(24)] for i1 in range(24)]
    nB = sum(1 for i1 in range(24) for i2 in range(24) if Bt[i1][i2])
    triv_pairs = []
    cand = 0
    for i1 in range(24):
        for i2 in range(24):
            if not Bt[i1][i2]:
                continue
            for j1 in range(24):
                for j2 in range(24):
                    if not Bt[j1][j2]:
                        continue
                    cand += 1
                    G = spin4(T2[i1], T2[j1])
                    H = spin4(T2[i2], T2[j2])
                    if all(act(G, H, E) == zm_scale(E, 4) for E in Ebasis):
                        triv_pairs.append((T2[i1][0], T2[j1][0], T2[i2][0], T2[j2][0]))
    kern_names = ["(" + ("1" if t[0] > 0 else "-1") + "," + ("1" if t[1] > 0 else "-1") + ")" for t in triv_pairs]
    meas("M4D.TS5", "kernel_of_the_action_on_(2Tx2T)x(2Tx2T)", len(triv_pairs) == 2 and cand == 576,
         [("pairs_total", 576 * 576), ("pairs_with_G_H^dag=4I", cand), ("acting_trivially", len(triv_pairs)),
          ("kernel", lst(kern_names) if kern_names else "none"),
          ("reading", "trivial on all E_ab forces G H^dag = 4I by linearity, so the 576 candidates are exhaustive")])

    # ------------------------------------------------ controls that must fail
    def gamma_plus(x):
        Z = pm_zero(4, 4)
        return pblock(Z, Lmat_p(x), Lmat_p(qconj_p(x)), Z)

    def gamma_noconj(x):
        Z = pm_zero(4, 4)
        return pblock(Z, Lmat_p(x), pm_scale(Lmat_p(x), -1), Z)

    gp = gamma_plus(X)
    res = identity(pm_mul(gp, gp), pm_diag(8, pscale(NSQ, -1)), XI, 2)
    control("K4D.G1", "gamma_with_+L_conj(xi)_squares_to_-nsq", not res[0], [("mismatches", res[2])])
    gn = gamma_noconj(X)
    res = identity(pm_mul(gn, gn), pm_diag(8, pscale(NSQ, -1)), XI, 2)
    control("K4D.G2", "gamma_without_conj_squares_to_-nsq", not res[0], [("mismatches", res[2])])
    LiLi = iblock(Lmat_i(qi), Z4, Z4, Lmat_i(qi))
    control("K4D.G3", "diag(L_i,L_i)_in_commutant", not all(imul(LiLi, g) == imul(g, LiLi) for g in gens), [])
    FBw = [FB[0], [0, 1, 0, 0, 1, 0], FB[2]]
    CW = C_of(X, FBw)
    res = identity(pm_mul(CW, pm_T(CW)), pm_diag(3, NSQ), XI, 2)
    sample_ok = imul(pm_eval(CW, E1PT), iT(pm_eval(CW, E1PT))) == ieye(3)
    control("K4D.S1", "wrong_f2=e13+e24_C'C'^T=nsq I", not res[0],
            [("mismatches", res[2]), ("passes_at_sample_xi=e1", 1 if sample_ok else 0),
             ("reading", "a sampled covector would have accepted it; the polynomial comparison does not")])
    res = identity(S2, pm_diag(4, pscale(NSQ, 2)), XI, 2)
    control("K4D.S2", "R1_with_c=1_raw_projection_isometric", not res[0], [("mismatches", res[2])])
    prod = qmul_p(qconj_p(X), AQ)
    res = identity(colvec(Ma), colvec([pscale(p, -1) for p in prod]), XI + AV, 2)
    control("K4D.S3", "M(xi)a=-conj(xi)*a", not res[0], [("mismatches", res[2]), ("other_non_matching_candidates", 15)])
    ETAB = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    for a in range(3):
        for m in range(3):
            for n in range(3):
                ETAB[a][m][n] = eps3(a, m, n)
        ETAB[a][a][3] = -1
        ETAB[a][3][a] = 1
    etab_forms = [[ETAB[a][i][j] for (i, j) in PAIRS] for a in range(3)]
    etab_sd = sum(1 for e in etab_forms if matvec_i(STAR, e) == e)
    TB = [[psum(pscale(X[mu], ETAB[a][mu][nu]) for mu in range(4)) for nu in range(4)] for a in range(3)]
    res = identity(TB, pm_mul(pm_int(Pm), CX), XI, 1)
    control("K4D.S4", "anti_self_dual_eta_bar_equals_P_C", not res[0] and etab_sd == 0,
            [("mismatches", res[2]), ("eta_bar_self_dual", etab_sd)])
    control("K4D.S5", "Psi_kernel=Lambda2+", on_sd == 3, [("f_basis_forms_with_nonzero_Psi", on_sd)])

    emit("MEASURE_SUMMARY", [("meas_rows", STATS["meas"]), ("meas_pass", STATS["meas_pass"]), ("controls", STATS["ctrl"]),
                             ("controls_failed_as_required", STATS["ctrl_ok"]),
                             ("measured_block_sha256", hashlib.sha256(("\n".join(LINES) + "\n").encode("utf-8")).hexdigest())])


def geval(p, pt):
    re_ = 0
    im_ = 0
    for m, (a, b) in p.items():
        t = 1
        for v in m:
            t = t * pt[v]
        re_ += a * t
        im_ += b * t
    return (re_, im_)


def _k_by_q(T2, MX, MXe1):
    out = {}
    for P in T2:
        Pp = [pconst(x) for x in P]
        for Q in T2:
            Qb = qconj_i(Q)
            Y = qmul_p(qmul_p(Pp, X), [pconst(x) for x in Qb])
            Rho = imul(Lmat_i(P), Rmat_i(Qb))
            K = imul(imul(pm_eval(M_of(Y), E1PT), Rho), iT(MXe1))
            out.setdefault(tuple(Q), set()).add(tuple(flat(K)))
    return out


# ================================================================ SUPPLEMENT
# Written after the measured block above was run and frozen (MEASURE-4D-BEFORE-COMPARE.hbp) and after the PRED, ADD and
# DIMS rows were read. It measures predicted details that the measured block did not reach. main_measure above is left
# byte-identical, so the helpers it nests are rebuilt here at module level.

SSTATS = {"rows": 0, "pass": 0}


def strict(pairs):
    for k, _v in pairs:
        if "=" in k or "|" in k:
            raise SystemExit("key holds = or |: " + k)
    return pairs


def supp(sid, name, ok, refs, pairs):
    strict(pairs)
    SSTATS["rows"] += 1
    if ok:
        SSTATS["pass"] += 1
    emit("SUPP", [("id", sid), ("name", name), ("for", refs)] + pairs
         + [("pass", 1 if ok else 0), ("measured_after_reading_predictions", 1), ("claim", "MEASURED")])


def zmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def zm_mul(A, B):
    n = len(A)
    k = len(B)
    m = len(B[0])
    out = [[(0, 0)] * m for _ in range(n)]
    for i in range(n):
        row = [0] * m
        rowi = [0] * m
        for t in range(k):
            a = A[i][t]
            if a == (0, 0):
                continue
            Bt = B[t]
            for j in range(m):
                b = Bt[j]
                if b != (0, 0):
                    row[j] += a[0] * b[0] - a[1] * b[1]
                    rowi[j] += a[0] * b[1] + a[1] * b[0]
        out[i] = [(row[j], rowi[j]) for j in range(m)]
    return out


def zm_add(A, B):
    return [[(A[i][j][0] + B[i][j][0], A[i][j][1] + B[i][j][1]) for j in range(len(A[0]))] for i in range(len(A))]


def zm_sub(A, B):
    return [[(A[i][j][0] - B[i][j][0], A[i][j][1] - B[i][j][1]) for j in range(len(A[0]))] for i in range(len(A))]


def zm_scale(A, k):
    return [[(k * x[0], k * x[1]) for x in r] for r in A]


def zm_eye(n, k=1):
    return [[(k, 0) if i == j else (0, 0) for j in range(n)] for i in range(n)]


def zm_dag(A):
    return [[(A[j][i][0], -A[j][i][1]) for j in range(len(A))] for i in range(len(A[0]))]


def zm_T(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def zm_conj(A):
    return [[(x[0], -x[1]) for x in r] for r in A]


def zm_kron(A, B):
    n = len(A)
    m = len(B)
    out = [[(0, 0)] * (n * m) for _ in range(n * m)]
    for i in range(n):
        for j in range(n):
            a = A[i][j]
            if a == (0, 0):
                continue
            for k in range(m):
                for l in range(m):
                    out[i * m + k][j * m + l] = zmul(a, B[k][l])
    return out


def zm_iszero(A):
    return all(x == (0, 0) for r in A for x in r)


def z_rank(mats):
    rows = []
    for M in mats:
        rows.append([x[0] for r in M for x in r] + [x[1] for r in M for x in r])
        rows.append([-x[1] for r in M for x in r] + [x[0] for r in M for x in r])
    return exact_div(rank_int(rows), 2)


def zmat_rank(M):
    n = len(M)
    m = len(M[0])
    rows = [[M[i][j][0] for j in range(m)] + [-M[i][j][1] for j in range(m)] for i in range(n)]
    rows += [[M[i][j][1] for j in range(m)] + [M[i][j][0] for j in range(m)] for i in range(n)]
    return exact_div(rank_int(rows), 2)


def U_int(P):
    a, b, c, d = P
    return [[(a, b), (c, d)], [(-c, d), (a, -b)]]


def gc_int(xq):
    Ux = U_int(xq)
    Ub = zm_scale(U_int(qconj_i(xq)), -1)
    out = zm_eye(4, 0)
    for i in range(2):
        for j in range(2):
            out[i][2 + j] = Ux[i][j]
            out[2 + i][j] = Ub[i][j]
    return out


def spin4_int(P, Q):
    out = zm_eye(4, 0)
    Up = U_int(P)
    Uq = U_int(Q)
    for i in range(2):
        for j in range(2):
            out[i][j] = Up[i][j]
            out[2 + i][2 + j] = Uq[i][j]
    return out


def s_algebra_dim(gens):
    n = len(gens[0])
    basis = [ieye(n)]
    frontier = [ieye(n)]
    cur = 1
    while frontier:
        new = []
        for B in frontier:
            for g in gens:
                P = imul(B, g)
                r = rank_int([flat(x) for x in basis] + [flat(P)])
                if r > cur:
                    basis.append(P)
                    new.append(P)
                    cur = r
        frontier = new
    return cur, basis


def s_commutant_dim(glist):
    n = len(glist[0])
    rows = []
    for g in glist:
        for r in range(n):
            for c in range(n):
                coeff = [0] * (n * n)
                for k in range(n):
                    coeff[n * r + k] += g[k][c]
                    coeff[n * k + c] -= g[r][k]
                rows.append(coeff)
    return n * n - rank_int(rows)


def s_eps3(a, b, c):
    if len({a, b, c}) < 3:
        return 0
    return perm_sign([a, b, c])


SIG = [gm_const(s) for s in ([[(0, 0), (1, 0)], [(1, 0), (0, 0)]], [[(0, 0), (0, -1)], [(0, 1), (0, 0)]],
                             [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]])]
TBS = [gm_scalar(s, {(): (0, 1)}) for s in SIG]


def ob_dA0(f, A):
    return [gm_add(gm_deriv(f, i), gm_comm(A[i], f)) for i in range(4)]


def ob_dA1(w, A):
    out = {}
    for (i, j) in PAIRS:
        out[(i, j)] = gm_sub(gm_add(gm_deriv(w[j], i), gm_comm(A[i], w[j])), gm_add(gm_deriv(w[i], j), gm_comm(A[j], w[i])))
    return out


def ob_curv(A):
    return {(i, j): gm_add(gm_sub(gm_deriv(A[j], i), gm_deriv(A[i], j)), gm_comm(A[i], A[j])) for (i, j) in PAIRS}


def gm_val(A):
    return [[A[i][j].get((), (0, 0)) for j in range(2)] for i in range(2)]


def gm_halve(V):
    return [[(exact_div(V[i][j][0], 2), exact_div(V[i][j][1], 2)) for j in range(2)] for i in range(2)]


def gv_str(V):
    return "[[" + ", ".join(str(x[0]) + ":" + str(x[1]) for x in V[0]) + "],[" + ", ".join(str(x[0]) + ":" + str(x[1]) for x in V[1]) + "]]"


def qname(v):
    names = ["", "i", "j", "k"]
    parts = []
    for t in range(4):
        c = v[t]
        if c == 0:
            continue
        mag = abs(c)
        s = ("-" if c < 0 else ("+" if parts else "")) + ("" if (mag == 1 and t > 0) else str(mag)) + names[t]
        parts.append(s)
    return "".join(parts) if parts else "0"


def main_supplement():
    gens = [gamma_i(E_UNIT[a]) for a in range(4)]
    Z4 = izero(4, 4)
    CX = C_of(X)
    MX = M_of(X)
    one, qi, qj, qk = E_UNIT

    cnt = 0
    for a in range(4):
        for b in range(a, 4):
            if iadd(imul(gens[a], gens[b]), imul(gens[b], gens[a])) == ieye(8, -2 if a == b else 0):
                cnt += 1
    supp("SUP4D.01", "anticommutator_relations_with_i<=j", cnt == 10, "4D.P2", [("hold", cnt), ("of", 10)])

    rdet = identity([[pdet(pm_mul(CX, pm_T(CX)))]], [[pmul(NSQ, pmul(NSQ, NSQ))]], XI, 6)
    minors = []
    sign = 1
    for j in range(4):
        minors.append(pscale(pdet([[CX[k][t] for t in range(4) if t != j] for k in range(3)]), sign))
        sign = -sign
    lam = 0
    for s in (1, -1):
        if identity(colvec(minors), colvec([pscale(pmul(NSQ, X[j]), s) for j in range(4)]), XI, 3)[0]:
            lam = s
    supp("SUP4D.02", "det(C C^T)=nsq^3_and_signed_maximal_minors_of_C", rdet[0] and lam != 0, "4D.P6",
         [("det_CCT_coeffs_checked", rdet[1]), ("det_CCT_mismatches", rdet[2]),
          ("minor_convention", "m_j = (-1)^(j-1) det(C(xi) without column j), j = 1..4"),
          ("m_equals_lambda_xi_with_lambda", ("+nsq" if lam > 0 else "-nsq") if lam else "neither")])

    rg = identity(pm_mul(pm_mul(pm_T(MX), pm_int(ieye(4, 2))), MX), pm_diag(4, pscale(NSQ, 2)), XI, 2)
    supp("SUP4D.03", "sigma_u^T G' sigma_u=2 nsq I_4_with_G'=2I_4", rg[0], "4D.P7", idrow(rg))

    PsiM = [[imul(gens[i], gens[j])[r][0] for (i, j) in PAIRS] for r in range(8)]
    cols = [[PsiM[r][I] for r in range(4)] for I in range(6)]
    pred_cols = [qi, qj, qk, qk, [-x for x in qj], qi]
    hand = [[-x for x in qmul_i(E_UNIT[i], qconj_i(E_UNIT[j]))] for (i, j) in PAIRS]
    on_f = [[sum(PsiM[r][I] * FB[k][I] for I in range(6)) for r in range(4)] for k in range(3)]
    Psi3 = [PsiM[r] for r in (1, 2, 3)]
    ppt = imul(Psi3, iT(Psi3)) == ieye(3, 2)
    pstar = imul(PsiM, iadd(STAR, ieye(6))) == iscale(PsiM, 2)
    cl_ok = sum(1 for I, (i, j) in enumerate(PAIRS)
                if [row[:4] for row in imul(gens[i], gens[j])[:4]] == Lmat_i(cols[I]))
    asd0 = sum(1 for g in GB if matvec_i(PsiM, g) == [0] * 8)
    supp("SUP4D.04", "Psi_4_values_on_e_ij_and_f_k", cols == pred_cols and cols == hand and on_f == [iscale([qi], 2)[0], iscale([qj], 2)[0], iscale([qk], 2)[0]]
         and ppt and pstar and cl_ok == 6 and asd0 == 3, "4D.P8",
         [("Psi_e12_e13_e14_e23_e24_e34", lst([qname(c) for c in cols])), ("Psi_f1_f2_f3", lst([qname(c) for c in on_f])),
          ("equals_-e_i*conj(e_j)", 1 if cols == hand else 0), ("Psi_Psi^T_is_2I_3", 1 if ppt else 0),
          ("Psi(*+I)_is_2Psi", 1 if pstar else 0), ("c(e_ij)_on_S+_is_L_Psi(e_ij)", lst([cl_ok, 6])), ("Lambda2-_to_0", asd0)])

    Ma = pm_vec(MX, AQ)
    rhs = [pscale(p, -1) for p in qmul_p(X, qconj_p(AQ))]
    bil = set((i, 4 + j) for i in range(4) for j in range(4))
    outside = sum(1 for p in Ma + rhs for m in p if m not in bil)
    mism = sum(len(psub(Ma[t], rhs[t])) for t in range(4))
    w = qmul_p(AQ, qconj_p(X))
    tens = identity(colvec([pscale(Ma[0], -1), Ma[1], Ma[2], Ma[3]]), colvec(w), XI + AV, 2)
    supp("SUP4D.05", "M(xi)a=-xi*conj(a)_bilinear_and_tensor_form", outside == 0 and mism == 0 and tens[0], "4D.P9",
         [("bilinear_coefficients_checked", 4 * len(bil)), ("mismatches", mism), ("monomials_outside_bilinear", outside),
          ("(s;c)->(-s,c)_of_M(xi)a_is_a*conj(xi)", 1 if tens[0] else 0)])

    def cform(vec6):
        out = [[0] * 8 for _ in range(8)]
        for I, (i, j) in enumerate(PAIRS):
            if vec6[I]:
                out = iadd(out, iscale(imul(gens[i], gens[j]), vec6[I]))
        return out

    def blk(Mm, r0, c0):
        return [row[c0:c0 + 4] for row in Mm[r0:r0 + 4]]

    f_kill_minus = sum(1 for f in FB if blk(cform(f), 4, 4) == Z4 and blk(cform(f), 0, 0) != Z4)
    g_kill_plus = sum(1 for g in GB if blk(cform(g), 0, 0) == Z4 and blk(cform(g), 4, 4) != Z4)
    even_diag = sum(1 for v in FB + GB if blk(cform(v), 0, 4) == Z4 and blk(cform(v), 4, 0) == Z4)
    supp("SUP4D.06", "Lambda2+_acts_on_S+_kills_S-_and_Lambda2-_kills_S+", f_kill_minus == 3 and g_kill_plus == 3 and even_diag == 6,
         "4D.P10", [("Lambda2+_forms_killing_S-_acting_on_S+", f_kill_minus), ("Lambda2-_forms_killing_S+_acting_on_S-", g_kill_plus),
                    ("block_diagonal", even_diag)])

    inf_ok = 0
    pres = 0
    inf_slots = 0
    for (i, j) in PAIRS:
        R = izero(4, 4)
        R[j][i] = 1
        R[i][j] = -1
        D = izero(6, 6)
        for J, (p, q) in enumerate(PAIRS):
            ep = [1 if t == p else 0 for t in range(4)]
            eq = [1 if t == q else 0 for t in range(4)]
            Rp = matvec_i(R, ep)
            Rq = matvec_i(R, eq)
            for I, (u, v) in enumerate(PAIRS):
                D[I][J] += (Rp[u] * eq[v] - Rp[v] * eq[u]) + (ep[u] * Rq[v] - ep[v] * Rq[u])
        mcoef = [[exact_div(sum(matvec_i(D, FB[k])[t] * FB[l][t] for t in range(6)), 2) for k in range(3)] for l in range(3)]
        if all(matvec_i(D, FB[k]) == [sum(mcoef[l][k] * FB[l][t] for l in range(3)) for t in range(6)] for k in range(3)):
            pres += 1
        nR = [[-mcoef[l][k] for l in range(3)] for k in range(3)]
        AR = [[0, 0, 0, 0]] + [[0] + nR[k] for k in range(3)]
        Rxi = pm_vec(pm_int(R), X)
        Ra = pm_vec(pm_int(R), AQ)
        lhs = [padd(x, y) for x, y in zip(pm_vec(M_of(Rxi), AQ), pm_vec(MX, Ra))]
        res = identity(colvec(lhs), colvec(pm_vec(pm_int(AR), Ma)), XI + AV, 2)
        inf_slots += res[1]
        if res[0]:
            inf_ok += 1
    Q8 = [[s if t == u else 0 for t in range(4)] for u in range(4) for s in (1, -1)]
    fin_ok = 0
    fin_slots = 0
    for p in Q8:
        Ad = imul(Lmat_i(p), Rmat_i(qconj_i(p)))
        Ap = [[1, 0, 0, 0]] + [[0] + Ad[r][1:] for r in (1, 2, 3)]
        adblock = Ad[0] == [1, 0, 0, 0] and [Ad[r][0] for r in (1, 2, 3)] == [0, 0, 0]
        for q in Q8:
            B = imul(Lmat_i(p), Rmat_i(qconj_i(q)))
            Y = qmul_p(qmul_p([pconst(x) for x in p], X), [pconst(x) for x in qconj_i(q)])
            res = identity(M_of(Y), pm_mul(pm_mul(pm_int(Ap), MX), pm_int(iT(B))), XI, 1)
            fin_slots += res[1]
            if res[0] and adblock:
                fin_ok += 1
    supp("SUP4D.07", "equivariance_6_rotations_R_ij_and_Q8xQ8_finite_form", inf_ok == 6 and pres == 6 and fin_ok == 64, "4D.P11",
         [("R_ij_identities_hold", lst([inf_ok, 6])), ("derivation_preserves_Lambda2+", lst([pres, 6])),
          ("rho_on_Lambda0", 0), ("infinitesimal_coeffs_checked", inf_slots),
          ("Q8xQ8_pairs_M(p xi conj(q))_is_A_p M(xi) B_pq^T", lst([fin_ok, 64])), ("finite_coeffs_checked", fin_slots),
          ("A_p", "diag(1, Ad_p), Ad_p the 3 x 3 block of y -> p y conj(p) on i, j, k")])

    A_na = [TBS[0], TBS[1], gm_zero(), gm_zero()]
    f0 = TBS[0]
    ddf = ob_dA1(ob_dA0(f0, A_na), A_na)
    v12 = gm_val(ddf[(0, 1)])
    others0 = all(gm_iszero(ddf[I]) for I in PAIRS if I != (0, 1))
    F = ob_curv(A_na)
    cdd = [gm_val(gm_zero()) for _ in range(3)]
    cF = [gm_val(gm_zero()) for _ in range(3)]
    gF = [gm_val(gm_zero()) for _ in range(3)]
    for k in range(3):
        acc = gm_zero()
        accF = gm_zero()
        accG = gm_zero()
        for Ii, I in enumerate(PAIRS):
            if FB[k][Ii]:
                acc = gm_add(acc, gm_scalar(ddf[I], {(): (FB[k][Ii], 0)}))
                accF = gm_add(accF, gm_scalar(F[I], {(): (FB[k][Ii], 0)}))
            if GB[k][Ii]:
                accG = gm_add(accG, gm_scalar(F[I], {(): (GB[k][Ii], 0)}))
        cdd[k] = gm_halve(gm_val(acc))
        cF[k] = gm_halve(gm_val(accF))
        gF[k] = gm_halve(gm_val(accG))
    zero2 = [[(0, 0), (0, 0)], [(0, 0), (0, 0)]]
    pred_dd = [[(0, 0), (4, 0)], [(-4, 0), (0, 0)]]
    pred_sd = [[(0, 0), (2, 0)], [(-2, 0), (0, 0)]]
    minus_i_s3 = [[(0, -1), (0, 0)], [(0, 0), (0, 1)]]
    okp2 = (v12 == pred_dd and others0 and cdd == [pred_sd, zero2, zero2] and cF == [minus_i_s3, zero2, zero2]
            and gF == [minus_i_s3, zero2, zero2])
    supp("SUP4D.08", "explicit_d_A_d_A(i sigma1)_and_F^+-", okp2, "OB.P2",
         [("d_A_d_A(i sigma1)_dx12", gv_str(v12)), ("other_components_zero", 1 if others0 else 0),
          ("F^+_f_coefficients", "; ".join(gv_str(V) for V in cF)), ("F^-_asd_coefficients", "; ".join(gv_str(V) for V in gF)),
          ("d_A^+_d_A(i sigma1)_f_coefficients", "; ".join(gv_str(V) for V in cdd)),
          ("normalization", "f-coefficient = <., f_k> over 2 (Gram 2I); M4D.OB5 reports <., f_k> itself")])

    MONO2 = [()] + [(t,) for t in range(4)] + list(itertools.combinations_with_replacement(range(4), 2))
    fgen = gm_zero()
    for b in range(3):
        gb = {}
        for idx, m in enumerate(MONO2):
            gb[tuple(sorted(m + (100 + 15 * b + idx,)))] = (1, 0)
        fgen = gm_add(fgen, gm_scalar(TBS[b], gb))
    A_c = [TBS[2], gm_scalar(TBS[2], {(): (2, 0)}), gm_zero(), gm_zero()]
    Fc = ob_curv(A_c)
    ddc = ob_dA1(ob_dA0(fgen, A_c), A_c)
    d1terms = sum(len(x) for wv in ob_dA0(fgen, A_c) for x in gm_flat(wv))
    okp4 = all(gm_iszero(Fc[I]) for I in PAIRS) and all(gm_iszero(ddc[I]) for I in PAIRS) and d1terms > 0
    supp("SUP4D.09", "commuting_A1=i sigma3_A2=2i sigma3_F=0_d_A_d_A=0", okp4, "OB.P4",
         [("F_zero", 1 if all(gm_iszero(Fc[I]) for I in PAIRS) else 0), ("d_A_d_A_f_zero_for_generic_f", 1 if all(gm_iszero(ddc[I]) for I in PAIRS) else 0),
          ("d_A_f_terms", d1terms)])

    gcs = [gc_int(E_UNIT[a]) for a in range(4)]
    wc = zm_mul(zm_mul(zm_mul(gcs[0], gcs[1]), gcs[2]), gcs[3])
    wdiag = [wc[i][i][0] for i in range(4)]
    G8 = [zm_kron(g, zm_eye(4)) for g in gcs] + [zm_kron(wc, g) for g in gcs]
    ac = 0
    for A in range(8):
        for B in range(8):
            if zm_add(zm_mul(G8[A], G8[B]), zm_mul(G8[B], G8[A])) == zm_eye(16, -2 if A == B else 0):
                ac += 1
    vol = G8[0]
    for A in range(1, 8):
        vol = zm_mul(vol, G8[A])
    vol_is_kron = vol == zm_kron(wc, wc)
    plus_dim = 16 - zmat_rank(zm_sub(vol, zm_eye(16)))
    minus_dim = 16 - zmat_rank(zm_add(vol, zm_eye(16)))
    splus = [i for i in range(4) if wdiag[i] == -1]
    sminus = [i for i in range(4) if wdiag[i] == 1]
    ev = {}
    for name, (X1, X2) in (("S+(x)S+", (splus, splus)), ("S-(x)S-", (sminus, sminus)), ("S+(x)S-", (splus, sminus)), ("S-(x)S+", (sminus, splus))):
        ev[name] = sorted(set(vol[4 * i + j][4 * i + j][0] for i in X1 for j in X2))
    okt1 = (ac == 64 and vol_is_kron and plus_dim == 8 and minus_dim == 8 and ev["S+(x)S+"] == [1] and ev["S-(x)S-"] == [1]
            and ev["S+(x)S-"] == [-1] and ev["S-(x)S+"] == [-1])
    supp("SUP4D.10", "graded_tensor_C^4(x)C^4=C^16_Cl_8_and_chirality", okt1, "TS.P1",
         [("generators", "gamma_a (x) 1 and omega4 (x) gamma_b on C^16, complex 4D gammas built from U"),
          ("omega4_complex_diag", lst(wdiag)), ("anticommutators_hold", lst([ac, 64])), ("vol8_is_omega4(x)omega4", 1 if vol_is_kron else 0),
          ("vol8_plus1_dim", plus_dim), ("vol8_minus1_dim", minus_dim),
          ("vol8_on_blocks", "; ".join(k + " " + lst(v) for k, v in ev.items())),
          ("reading", "vol8 = gamma_1 ... gamma_8 is +1 on S+(x)S+ + S-(x)S- and -1 on the mixed blocks, 8 + 8")])

    eps = [[(0, 0), (1, 0)], [(-1, 0), (0, 0)]]
    J = zm_eye(4, 0)
    for i in range(2):
        for j in range(2):
            J[i][j] = eps[i][j]
            J[2 + i][2 + j] = eps[i][j]
    skew = zm_T(J) == zm_scale(J, -1)
    T2 = [[s if t == u else 0 for t in range(4)] for s in (2, -2) for u in range(4)]
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    T2.append([s0, s1, s2, s3])
    inv = 0
    for P in T2:
        for Q in T2:
            G = spin4_int(P, Q)
            if zm_mul(zm_mul(zm_T(G), J), G) == zm_scale(J, 4):
                inv += 1
    quat = sum(1 for P in T2 if zm_mul(U_int(P), eps) == zm_mul(eps, zm_conj(U_int(P))))
    GEN4 = [((0, 2, 0, 0), (2, 0, 0, 0)), ((1, 1, 1, 1), (2, 0, 0, 0)), ((2, 0, 0, 0), (0, 2, 0, 0)), ((2, 0, 0, 0), (1, 1, 1, 1))]
    GM = [spin4_int(P, Q) for (P, Q) in GEN4]
    JT = zm_T(J)
    jint = 0
    jn = 0
    for G in GM:
        for H in GM:
            for a in range(4):
                for b in range(4):
                    jn += 1
                    Gu = [[G[i][a]] for i in range(4)]
                    Hv = [[H[i][b] for i in range(4)]]
                    lhs = zm_mul(zm_mul(Gu, Hv), JT)
                    E = zm_eye(4, 0)
                    E[a][b] = (1, 0)
                    rhs = zm_mul(zm_mul(G, zm_mul(E, JT)), zm_dag(H))
                    if lhs == rhs:
                        jint += 1
    jj = zm_mul(J, J) == zm_eye(4, -1)
    JJ = zm_kron(J, J)
    real16 = zm_mul(JJ, zm_conj(JJ)) == zm_eye(16)
    supp("SUP4D.11", "spin_invariant_form_J_skew_and_the_J_identification", skew and inv == 576 and jint == jn and quat == 24, "TS.P2",
         [("J", "diag(eps, eps), eps = [[0, 1],[-1, 0]]"), ("J_skew", 1 if skew else 0), ("G^T J G_is_4J_on_2Tx2T", lst([inv, 576])),
          ("det_J_block_S+", eps[0][0][0] * eps[1][1][0] - eps[0][1][0] * eps[1][0][0]),
          ("det_J_block_S-", J[2][2][0] * J[3][3][0] - J[2][3][0] * J[3][2][0]),
          ("u(x)v->u v^T J^T_intertwines_g(x)h_with_A->gAh^-1", lst([jint, jn])),
          ("U(p)eps_is_eps conj(U(p))_on_2T", lst([quat, 24]))])

    prods = {}
    for k in range(5):
        for S in itertools.combinations(range(4), k):
            P = zm_eye(4)
            for a in S:
                P = zm_mul(P, gcs[a])
            prods[S] = P
    deg_counts = [sum(1 for S in prods if len(S) == k) for k in range(5)]
    Pp = zm_sub(zm_eye(4), wc)
    Pm = zm_add(zm_eye(4), wc)

    def sand(L, Xm, R):
        return zm_mul(zm_mul(L, Xm), R)

    def c2(vec6):
        out = zm_eye(4, 0)
        for I, (i, j) in enumerate(PAIRS):
            if vec6[I]:
                out = zm_add(out, zm_scale(prods[(i, j)], vec6[I]))
        return out

    even = [prods[S] for S in prods if len(S) in (0, 2, 4)]
    odd = [prods[S] for S in prods if len(S) in (1, 3)]
    l1 = [prods[S] for S in prods if len(S) == 1]
    l3 = [prods[S] for S in prods if len(S) == 3]
    endp = z_rank([sand(Pp, x, Pp) for x in even])
    endm = z_rank([sand(Pm, x, Pm) for x in even])
    l2p_on_p = z_rank([sand(Pp, c2(f), Pp) for f in FB])
    l2m_on_p = sum(1 for g in GB if zm_iszero(sand(Pp, c2(g), Pp)))
    l2m_on_m = z_rank([sand(Pm, c2(g), Pm) for g in GB])
    l2p_on_m = sum(1 for f in FB if zm_iszero(sand(Pm, c2(f), Pm)))
    endp_from = z_rank([sand(Pp, zm_eye(4), Pp)] + [sand(Pp, c2(f), Pp) for f in FB])
    endm_from = z_rank([sand(Pm, wc, Pm)] + [sand(Pm, c2(g), Pm) for g in GB])
    odd_rank = z_rank(odd)
    odd_offdiag = sum(1 for x in odd if zm_iszero(sand(Pp, x, Pp)) and zm_iszero(sand(Pm, x, Pm)))
    hom_pm = z_rank([sand(Pp, x, Pm) for x in odd])
    hom_mp = z_rank([sand(Pm, x, Pp) for x in odd])
    l1_pm = z_rank([sand(Pp, x, Pm) for x in l1])
    l3_mp = z_rank([sand(Pm, x, Pp) for x in l3])
    l1_mp = z_rank([sand(Pm, x, Pp) for x in l1])
    l1_cap = 4 + hom_pm - z_rank(l1 + [sand(Pp, x, Pm) for x in odd])
    full = z_rank(list(prods.values()))
    okt3 = (deg_counts == [1, 4, 6, 4, 1] and full == 16 and endp == 4 and endm == 4 and l2p_on_p == 3 and l2m_on_p == 3
            and l2m_on_m == 3 and l2p_on_m == 3 and endp_from == 4 and endm_from == 4 and odd_rank == 8 and odd_offdiag == 8
            and hom_pm == 4 and hom_mp == 4 and l1_pm == 4 and l3_mp == 4 and jj and real16)
    supp("SUP4D.12", "End(C^4)=Cl_4(x)C_by_degree_and_chirality_blocks", okt3, "TS.P3",
         [("degree_counts", lst(deg_counts)), ("complex_rank_all_16", full), ("End(S+)_from_even", endp), ("End(S-)_from_even", endm),
          ("End(S+)_from_Lambda0+Lambda2+", endp_from), ("End(S-)_from_Lambda4+Lambda2-", endm_from),
          ("Lambda2+_rank_on_S+", l2p_on_p), ("Lambda2-_zero_on_S+", l2m_on_p), ("Lambda2-_rank_on_S-", l2m_on_m), ("Lambda2+_zero_on_S-", l2p_on_m),
          ("odd_rank", odd_rank), ("odd_off_diagonal", odd_offdiag), ("Hom(S-,S+)_dim", hom_pm), ("Hom(S+,S-)_dim", hom_mp),
          ("Lambda1_onto_Hom(S-,S+)_rank", l1_pm), ("Lambda1_onto_Hom(S+,S-)_rank", l1_mp), ("Lambda3_onto_Hom(S+,S-)_rank", l3_mp),
          ("Lambda1_cap_Hom(S-,S+)_as_subspaces", l1_cap),
          ("J^2_is_-I_4_quaternionic", 1 if jj else 0), ("(J(x)J)conj(J(x)J)_is_I_16_real", 1 if real16 else 0),
          ("P+_scaled", "I - omega4 (omega4 = -1 on S+)")])

    gpl = [iblock(Z4, Lmat_i(E_UNIT[a]), Lmat_i(qconj_i(E_UNIT[a])), Z4) for a in range(4)]
    GP = pblock(pm_zero(4, 4), Lmat_p(X), Lmat_p(qconj_p(X)), pm_zero(4, 4))
    rsq = identity(pm_mul(GP, GP), pm_diag(8, NSQ), XI, 2)
    adim, abasis = s_algebra_dim(gpl)
    cdim = s_commutant_dim(gpl)
    supp("SUP4D.13", "positive_sign_Clifford_Cl(R^4,+nsq)_on_R^8", rsq[0] and adim == 16 and cdim == 4, "BOTT n=4",
         [("gamma'", "[[0, L_xi],[L_conj(xi), 0]]"), ("gamma'(xi)^2_is_+nsq I_8_coeffs_checked", rsq[1]), ("mismatches", rsq[2]),
          ("generated_algebra_dim", adim), ("commutant_dim", cdim)])
    ga_dim, ga_basis = s_algebra_dim(gens)
    gc_dim = s_commutant_dim(gens)
    VV8 = [pvar(40 + t) for t in range(8)]
    acc = pm_zero(8, 8)
    for k in range(5):
        for S in itertools.combinations(range(4), k):
            P = ieye(8)
            for a in S:
                P = imul(P, gens[a])
            Pv = pm_vec(pm_int(P), VV8)
            acc = pm_add(acc, pm_outer(Pv, Pv))
    rirr = identity(acc, pm_diag(8, pscale(nsq_of(VV8), 2)), list(range(40, 48)), 2)
    supp("SUP4D.14", "negative_sign_Clifford_recount_and_irreducibility_certificate", ga_dim == 16 and gc_dim == 4 and rirr[0], "BOTT n=4",
         [("generated_algebra_dim", ga_dim), ("commutant_dim", gc_dim), ("module_dim", 8),
          ("certificate", "sum over the 16 monomials gamma_I v (gamma_I v)^T = 2 nsq(v) I_8 for every v, so the 16 vectors gamma_I v span R^8 for every v != 0 and R^8 has no proper invariant subspace"),
          ("certificate_coeffs_checked", rirr[1]), ("certificate_mismatches", rirr[2])])

    Dm = pblock(pm_zero(4, 4), pm_T(MX), MX, pm_zero(4, 4))
    rD = identity(pm_mul(Dm, Dm), pm_diag(8, NSQ), XI, 2)
    symm = pm_T(Dm) == Dm
    supp("SUP4D.15", "doubled_self_adjoint_rolled_up_symbol_squares_to_+nsq", rD[0] and symm, "ADD R6",
         [("D", "[[0, M^T],[M, 0]]"), ("D_symmetric", 1 if symm else 0)] + idrow(rD))

    ETA = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    ETP = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    for a in range(3):
        for m in range(3):
            for n in range(3):
                ETA[a][m][n] = s_eps3(a, m, n)
                ETP[a][1 + m][1 + n] = s_eps3(a, m, n)
        ETA[a][a][3] = 1
        ETA[a][3][a] = -1
        ETP[a][0][1 + a] = 1
        ETP[a][1 + a][0] = -1
    tb = [psum(pmul(pscale(X[mu], ETA[a][mu][nu]), AQ[nu]) for mu in range(4) for nu in range(4)) for a in range(3)]
    tp = [psum(pmul(pscale(X[mu], ETP[a][mu][nu]), AQ[nu]) for mu in range(4) for nu in range(4)) for a in range(3)]
    imw = [w[1], w[2], w[3]]
    Pm3 = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    r_t = identity(colvec(tb), colvec([psum(pscale(imw[k], Pm3[a][k]) for k in range(3)) for a in range(3)]), XI + AV, 2)
    r_tp = identity(colvec(tp), colvec(imw), XI + AV, 2)
    gram_eta = [[sum(ETA[a][i][j] * ETA[b][i][j] for (i, j) in PAIRS) for b in range(3)] for a in range(3)]
    supp("SUP4D.16", "tHooft_coefficient_vector_and_Im(a*conj(xi))", r_t[0] and r_tp[0] and gram_eta == ieye(3, 2), "ADD R2",
         [("t_a_is_sum_k P_ak Im(a*conj(xi))_k", 1 if r_t[0] else 0), ("P", mat_str(Pm3)),
          ("eta'_(distinguished e1)_t'_is_Im(a*conj(xi))", 1 if r_tp[0] else 0), ("eta_Gram", mat_str(gram_eta)),
          ("coeffs_checked", r_t[1] + r_tp[1]), ("mismatches", r_t[2] + r_tp[2])])

    t_small = [matvec_i(STAR, f) == f for f in FB]
    t_big = [matvec_i(STAR, g) == [-x for x in g] for g in GB]
    supp("SUP4D.17", "eigenvalues_of_T=*_on_Lambda2+-", all(t_small) and all(t_big), "DIMS 4D",
         [("T_on_Lambda2+", 1 if all(t_small) else 0), ("T_on_Lambda2-", -1 if all(t_big) else 0),
          ("k_from_T_small_plus_1", 2 if all(t_small) else 0)])

    eqkeys = []
    for ln in LINES:
        parts = ln.split("|")
        if parts[0] not in ("MEAS", "CONTROL"):
            continue
        mid = parts[1].partition("=")[2]
        for fld in parts[3:]:
            if fld.count("=") > 1 and fld.rpartition("=")[2] in ("0", "1"):
                eqkeys.append(mid + ":" + fld.rpartition("=")[0])
    emit("FORMAT_NOTE", [("about", "the frozen measured block"), ("keys_holding_equals", len(eqkeys)), ("keys", "; ".join(eqkeys)),
                         ("read_as", "for these fields the value is the text after the last =, 0 or 1; every other field splits at the first ="),
                         ("kept", "the measured block stays byte-identical to MEASURE-4D-BEFORE-COMPARE.hbp; SUPP and COMPARE keys hold no =")])
    emit("SUPPLEMENT_SUMMARY", [("rows", SSTATS["rows"]), ("pass", SSTATS["pass"]),
                                ("order", "run after the measured block was frozen and after the PRED, ADD and DIMS rows were read")])


# ================================================================ COMPARE

def find_row(tag, mid):
    for ln in LINES:
        parts = ln.split("|")
        if parts[0] == tag and ("id=" + mid) in parts[1:3]:
            return parts
    return None


def get(mid, key):
    for tag in ("MEAS", "SUPP", "CONTROL"):
        parts = find_row(tag, mid)
        if parts is not None:
            for p in parts[1:]:
                if p.startswith(key + "="):
                    return p[len(key) + 1:]
            return None
    return None


def parse_frozen(path, expect):
    with open(path, "rb") as fh:
        raw = fh.read()
    h = hashlib.sha256(raw).hexdigest()
    if h != expect:
        raise SystemExit("frozen file sha256 differs: " + path)
    rows = []
    for ln in raw.decode("utf-8").split("\n"):
        if not ln:
            continue
        parts = ln.split("|")
        fields = {}
        for p in parts[1:]:
            k, _, v = p.partition("=")
            fields[k] = v
        rows.append((parts[0], fields, hashlib.sha256(ln.encode("utf-8")).hexdigest()[:16]))
    return rows


CSTATS = {"rows": 0}
VERDICTS = {}
DISAGREE = []


def cmp_row(ref, rowsha, predicted, mrows, checks, verdict_ok, note):
    CSTATS["rows"] += 1
    if "=" in ref:
        raise SystemExit("ref holds =: " + ref)
    passes = all(get(m, "pass") == "1" for m in mrows)
    got = []
    ok = passes
    for (mid, key, want) in checks:
        v = get(mid, key)
        got.append(mid + "." + key + "=" + str(v))
        if v != want:
            ok = False
    if verdict_ok == "NOT_MEASURED":
        verdict = "NOT_MEASURED"
    else:
        verdict = verdict_ok if ok else "DISAGREE"
    if verdict.startswith("DISAGREE"):
        DISAGREE.append(ref)
    VERDICTS[verdict] = VERDICTS.get(verdict, 0) + 1
    emit("COMPARE", [("ref", ref), ("frozen_row_sha16", rowsha), ("predicted", predicted), ("measured_rows", lst(mrows) if mrows else "none"),
                     ("measured", "; ".join(got) if got else "none"), ("all_cited_rows_pass", 1 if passes else 0), ("verdict", verdict),
                     ("note", note), ("claim", "MEASURED" if verdict != "NOT_MEASURED" else "NAMED")])


def compare_section():
    prows = parse_frozen(os.path.join(HERE, "frozen", "PREDICTION-SPINOR-SYMBOL.hbp"), FROZEN_SHA)
    arows = parse_frozen(os.path.join(HERE, "frozen", "ADDENDUM-SPINOR-SYMBOL.hbp"), ADDENDUM_SHA)
    emit("COMPARE_HEAD", [("frozen_prediction_sha256", FROZEN_SHA), ("frozen_addendum_sha256", ADDENDUM_SHA),
                          ("scope", "PRED 4D.P1..P11, OB.P1..P4, TS.P1..P3; DIMS case=4D; BOTT n=4; ADD R1, R2, R5, R6, DISCLOSURE (4D parts); the 4D PLANTED row and the frozen COMPARE rows on 4D and the tensor"),
                          ("verdicts", "AGREE: every cited row passes and every quoted value matches; AGREE_SCOPED: the measured part agrees and the named remainder is listed in note; NOT_MEASURED; DISAGREE"),
                          ("first_seen", "PRED, ADD and DIMS values were read after MEASURE-4D-BEFORE-COMPARE.hbp was written; SUPP rows came after")])
    byid = {}
    for tag, f, sh in prows:
        if tag == "PRED":
            byid[f["id"]] = sh
        elif tag == "DIMS":
            byid["DIMS:" + f["case"]] = (sh, f)
        elif tag == "BOTT":
            byid["BOTT:" + f["n"]] = (sh, f)
        elif tag == "PLANTED":
            byid["PLANTED:" + f["plant"]] = (sh, f)
        elif tag == "COMPARE":
            byid["COMPARE:" + f["item"]] = (sh, f)
    for tag, f, sh in arows:
        if tag == "ADD":
            byid["ADD:" + f["id"]] = sh
    need = ["4D.P" + str(n) for n in range(1, 12)] + ["OB.P" + str(n) for n in range(1, 5)] + ["TS.P" + str(n) for n in range(1, 4)]
    for r in need:
        if r not in byid:
            raise SystemExit("frozen row missing " + r)

    cmp_row("4D.P1", byid["4D.P1"], "640 coefficients of gamma(xi)^2 + nsq I_8 all 0 (64 entries x 10 monomials)", ["M4D.G1"],
            [("M4D.G1", "coeffs_checked", "640"), ("M4D.G1", "mismatches", "0"), ("M4D.G1", "entries", "64"), ("M4D.G1", "monomials_per_entry", "10")],
            "AGREE", "same count of coefficients, all zero")
    cmp_row("4D.P2", byid["4D.P2"], "10 of 10 relations i <= j hold", ["M4D.G2", "SUP4D.01"],
            [("SUP4D.01", "hold", "10"), ("M4D.G2", "hold", "16")], "AGREE", "the measured block checked all 16 ordered pairs; the supplement counts the 10 with i <= j")
    cmp_row("4D.P3", byid["4D.P3"], "generated algebra dim 16, Cl_4 = M_2(H), faithful; commutant dim 4 = right H on both summands; 16 x 4 = 64",
            ["M4D.G3", "M4D.G4", "M4D.G5"],
            [("M4D.G3", "dim", "16"), ("M4D.G3", "rank_of_subset_products", "16"), ("M4D.G4", "dim", "4"), ("M4D.G4", "right_H_rank", "4"),
             ("M4D.G5", "center_dim", "1"), ("M4D.G5", "commutant_of_right_H_dim", "16")],
            "AGREE", "faithful: the 16 Clifford monomials stay independent (rank 16); M_2(H) read from center 1, commutant H and double commutant 16, iso type DERIVED; 16 x 4 = 64 is arithmetic")
    cmp_row("4D.P4", byid["4D.P4"], "omega4 = -diag(I_4, -I_4) = diag(-I_4, I_4); omega4^2 = I_8; anticommutes with every gamma(xi)", ["M4D.G6"],
            [("M4D.G6", "diag", "-1,-1,-1,-1,1,1,1,1"), ("M4D.G6", "offdiag_nonzero", "0"), ("M4D.G6", "square_is_I", "1"), ("M4D.G6", "anticommutes_with_generators", "4")],
            "AGREE", "anticommuting with the 4 generators gives every gamma(xi) by linearity")
    cmp_row("4D.P5", byid["4D.P5"], "dims 1,4,3; Euler 0; sigma1 sigma0 = 0 polynomially; rank 1, 3 and ker sigma1 = span(xi) for every xi != 0",
            ["M4D.S3", "M4D.S4", "M4D.S5", "M4D.S6"],
            [("M4D.S6", "dims", "1,4,3"), ("M4D.S6", "euler", "0"), ("M4D.S3", "mismatches", "0"), ("M4D.S4", "rank_sigma_d0", "1"),
             ("M4D.S4", "rank_C", "3"), ("M4D.S5", "CTC_mismatches", "0")],
            "AGREE", "C(xi) = 2 sigma1(xi) in f-coordinates; C^T C + xi xi^T = nsq I_4 gives ker = span(xi) for every xi != 0")
    cmp_row("4D.P6", byid["4D.P6"], "C xi = 0, C C^T = nsq I_3, det(C C^T) = nsq^3, signed maximal minors = lambda xi with lambda = +- nsq",
            ["M4D.S3", "M4D.S5", "SUP4D.02"],
            [("M4D.S5", "CCT_mismatches", "0"), ("SUP4D.02", "det_CCT_mismatches", "0"), ("SUP4D.02", "m_equals_lambda_xi_with_lambda", "+nsq")],
            "AGREE", "lambda = +nsq with m_j = (-1)^(j-1) det(C without column j); the sign depends on that convention, which the prediction leaves open")
    cmp_row("4D.P7", byid["4D.P7"], "M^T M = M M^T = nsq I_4; sigma_u^T G' sigma_u = 2 nsq I_4 with G' = 2 I_4", ["M4D.S7", "SUP4D.03"],
            [("M4D.S7", "MTM_mismatches", "0"), ("M4D.S7", "MMT_mismatches", "0"), ("SUP4D.03", "mismatches", "0")],
            "AGREE", "also measured: det M = -nsq^2 (M4D.S7), not predicted")
    cmp_row("4D.P8", byid["4D.P8"], "Psi(e12..e34) = i, j, k, k, -j, i; Psi(f1..f3) = 2i, 2j, 2k; 0 on Lambda2-; Psi Psi^T = 2 I_3; Psi(*+I) = 2 Psi; c(omega) on S+ = L_Psi(omega); ker = Lambda2-, dim 3",
            ["M4D.S16", "SUP4D.04"],
            [("SUP4D.04", "Psi_e12_e13_e14_e23_e24_e34", "i,j,k,k,-j,i"), ("SUP4D.04", "Psi_f1_f2_f3", "2i,2j,2k"), ("SUP4D.04", "Psi_Psi^T_is_2I_3", "1"),
             ("SUP4D.04", "Psi(*+I)_is_2Psi", "1"), ("SUP4D.04", "c(e_ij)_on_S+_is_L_Psi(e_ij)", "6,6"), ("M4D.S16", "kernel_dim", "3"), ("M4D.S16", "asd_in_kernel", "3")],
            "AGREE", "Psi_4(e_ij) = -e_i conj(e_j) holds for all six (SUP4D.04 equals_-e_i*conj(e_j)=1)")
    cmp_row("4D.P9", byid["4D.P9"], "M(xi) a read as a quaternion = -xi conj(a), 64 bilinear coefficients; (s; c) -> -s + c i: M(xi) a -> a conj(xi)",
            ["M4D.S13", "SUP4D.05"],
            [("M4D.S13", "matching", "-xi*conj(a)"), ("SUP4D.05", "bilinear_coefficients_checked", "64"), ("SUP4D.05", "mismatches", "0"),
             ("SUP4D.05", "(s;c)->(-s,c)_of_M(xi)a_is_a*conj(xi)", "1")],
            "AGREE", "-xi*conj(a) is the only one of 16 sign and conjugation candidates that matches; the measured block compared over all 36 degree-2 monomials in 8 variables (144 slots), the supplement over the 16 bilinear ones (64)")
    cmp_row("4D.P10", byid["4D.P10"], "real: Lambda0 + Lambda2+ = H with y -> p y conj(p), Lambda1 = H with a -> p a conj(q); Lambda2+ acts on S+ and kills S-, Lambda2- kills S+; complex: S+ (x) S+ = 1 + 3, S+ (x) S- = 4",
            ["M4D.E2", "SUP4D.06", "SUP4D.12"],
            [("M4D.E2", "K_equals_4_Ad(P)", "576"), ("M4D.E2", "K_depends_only_on_P", "1"), ("SUP4D.06", "Lambda2+_forms_killing_S-_acting_on_S+", "3"),
             ("SUP4D.06", "Lambda2-_forms_killing_S+_acting_on_S-", "3"), ("SUP4D.12", "End(S+)_from_Lambda0+Lambda2+", "4")],
            "AGREE_SCOPED", "real parts measured on 2T x 2T; complex S+ (x) S+ = Lambda2 S+ + Sym2 S+ = 1 + 3 measured only through End(S+) = Lambda0 + Lambda2+ (1 + 3); S+ (x) S- = 4 = Lambda1 (x) C read from Hom(S-, S+) dim 4 (SUP4D.12)")
    cmp_row("4D.P11", byid["4D.P11"], "6 rotations R_ij: sigma(R xi) a + sigma(xi)(R a) = rho(R) sigma(xi) a, rho = 0 on Lambda0, derivation on 2-forms; Q8 x Q8 (64 pairs): M(p xi conj(q)) = A_p M(xi) B_pq^T, A_p = diag(1, Ad_p)",
            ["SUP4D.07", "M4D.E2", "M4D.E3"],
            [("SUP4D.07", "R_ij_identities_hold", "6,6"), ("SUP4D.07", "derivation_preserves_Lambda2+", "6,6"),
             ("SUP4D.07", "Q8xQ8_pairs_M(p xi conj(q))_is_A_p M(xi) B_pq^T", "64,64"), ("M4D.E2", "hold", "576"), ("M4D.E3", "hold", "6")],
            "AGREE", "the measured block had the finite form on all 576 pairs of 2T x 2T (which contains Q8 x Q8) scaled by 2, and gamma's infinitesimal form; the supplement adds the 6 R_ij on M and the unscaled Q8 x Q8 form")
    cmp_row("OB.P1", byid["OB.P1"], "A1 = [[0:0, 0:1],[0:1, 0:0]], A2 = [[0:0, 1:0],[-1:0, 0:0]], F_12 = -2 i sigma3 = [[0:-2, 0:0],[0:0, 0:2]] != 0", ["M4D.OB1"],
            [("M4D.OB1", "F12", "[[0,-2, 0,0],[0,0, 0,2]] as re:im pairs"), ("M4D.OB1", "F12_su2_coords_in_i_sigma_basis", "0,0,-2")],
            "AGREE", "A1, A2 built from sigma1, sigma2 of OB.C1 times i; the 9 Pauli commutators of OB.C1 hold")
    cmp_row("OB.P2", byid["OB.P2"], "d_A d_A (i sigma1) = 4 i sigma2 dx12 = [[0:0, 4:0],[-4:0, 0:0]]; F^+ = -i sigma3 f1; d_A^+ d_A (i sigma1) = 2 i sigma2 f1 = [[0:0, 2:0],[-2:0, 0:0]] f1; F^- = -i sigma3 (e12 - e34)",
            ["M4D.OB2", "M4D.OB5", "SUP4D.08"],
            [("SUP4D.08", "d_A_d_A(i sigma1)_dx12", "[[0:0, 4:0],[-4:0, 0:0]]"), ("SUP4D.08", "other_components_zero", "1"),
             ("SUP4D.08", "F^+_f_coefficients", "[[0:-1, 0:0],[0:0, 0:1]]; [[0:0, 0:0],[0:0, 0:0]]; [[0:0, 0:0],[0:0, 0:0]]"),
             ("SUP4D.08", "F^-_asd_coefficients", "[[0:-1, 0:0],[0:0, 0:1]]; [[0:0, 0:0],[0:0, 0:0]]; [[0:0, 0:0],[0:0, 0:0]]"),
             ("SUP4D.08", "d_A^+_d_A(i sigma1)_f_coefficients", "[[0:0, 2:0],[-2:0, 0:0]]; [[0:0, 0:0],[0:0, 0:0]]; [[0:0, 0:0],[0:0, 0:0]]"),
             ("M4D.OB5", "F_plus_c_coords_su2", "0,0,-2; 0,0,0; 0,0,0")],
            "AGREE", "M4D.OB5 gives <F, f_k> = (-2 i sigma3, 0, 0); with Gram 2I that is F^+ = -i sigma3 f1, the predicted coefficient; generic f: d_A d_A f = [F_A, f] (M4D.OB2)")
    cmp_row("OB.P3", byid["OB.P3"], "symbol = sigma(d) (x) id, unchanged by [A, .]; dims 3, 12, 9, Euler 0, ranks 3, 9, exact for every xi != 0 while d_A^+ d_A != 0; a complex exactly when F_A^+ = 0", ["M4D.OB6", "M4D.OB5"],
            [("M4D.OB6", "dims", "3,12,9"), ("M4D.OB6", "euler", "0"), ("M4D.OB6", "rank_at_2357", "3,9"), ("M4D.OB6", "symbol_independent_of_A", "1"),
             ("M4D.OB6", "S1TS1+S0S0T_nsq_I12", "1"), ("M4D.OB5", "rank_of_f->[F^+,f]_on_su2", "2")],
            "AGREE_SCOPED", "the iff (a complex exactly when F_A^+ = 0) is NAMED; measured here: d_A^+ d_A = [F_A^+, .] for this A, nonzero since F_A^+ != 0")
    cmp_row("OB.P4", byid["OB.P4"], "A1 = i sigma3, A2 = 2 i sigma3: F_A = 0 and d_A d_A = 0", ["SUP4D.09", "M4D.OB4"],
            [("SUP4D.09", "F_zero", "1"), ("SUP4D.09", "d_A_d_A_f_zero_for_generic_f", "1")],
            "AGREE", "the measured block ran A1 = A2 = i sigma3 and A = 0 (M4D.OB4); the supplement runs the predicted 2 i sigma3")
    cmp_row("TS.P1", byid["TS.P1"], "C^4 (x) C^4 = C^16 as graded tensor product, Cl(V + W) = Cl(V) graded-(x) Cl(W); S+(8) = S+ (x) S+ + S- (x) S- (8), S-(8) = mixed (8)", ["SUP4D.10"],
            [("SUP4D.10", "anticommutators_hold", "64,64"), ("SUP4D.10", "vol8_plus1_dim", "8"), ("SUP4D.10", "vol8_minus1_dim", "8"),
             ("SUP4D.10", "vol8_on_blocks", "S+(x)S+ 1; S-(x)S- 1; S+(x)S- -1; S-(x)S+ -1")],
            "AGREE_SCOPED", "measured: the 8 graded generators satisfy Cl_8 on C^16 and vol8 = omega4 (x) omega4 is +1 on S+(x)S+ + S-(x)S-; calling that eigenspace S+(8) is the orientation convention, which 8D.C1 fixes and this 4D run does not test")
    cmp_row("TS.P2", byid["TS.P2"], "End(C^4) complex dim 16; invariant form skew, symplectic on S+ and S-; (g, h) acts by A -> g A h^-1, (-1, -1) trivial; the identification intertwines Spin(4) x Spin(4)",
            ["M4D.TS2", "M4D.TS3", "M4D.TS4", "M4D.TS5", "SUP4D.11"],
            [("M4D.TS2", "complex_rank_of_images_of_basis_tensors", "16"), ("M4D.TS4", "(-1,-1)_trivial", "1"), ("M4D.TS5", "kernel", "(1,1),(-1,-1)"),
             ("SUP4D.11", "J_skew", "1"), ("SUP4D.11", "G^T J G_is_4J_on_2Tx2T", "576,576"), ("SUP4D.11", "u(x)v->u v^T J^T_intertwines_g(x)h_with_A->gAh^-1", "256,256")],
            "AGREE_SCOPED", "factoring through the cover of SO(4) x SO(4) in Spin(8) is NAMED; measured: the kernel on (2T x 2T) x (2T x 2T) is exactly {(1,1), (-1,-1)}")
    cmp_row("TS.P3", byid["TS.P3"], "End(C^4) = Cl_4 (x) C = 1 + 4 + 6 + 4 + 1; End(S+) = Lambda0 + Lambda2+, End(S-) = Lambda4 + Lambda2- (1 + 3 each); Hom(S-, S+), Hom(S+, S-) = Lambda1, Lambda3 (4 + 4); J (x) J squares to +1",
            ["SUP4D.12", "M4D.TS2"],
            [("SUP4D.12", "degree_counts", "1,4,6,4,1"), ("SUP4D.12", "complex_rank_all_16", "16"), ("SUP4D.12", "End(S+)_from_Lambda0+Lambda2+", "4"),
             ("SUP4D.12", "End(S-)_from_Lambda4+Lambda2-", "4"), ("SUP4D.12", "Hom(S-,S+)_dim", "4"), ("SUP4D.12", "Hom(S+,S-)_dim", "4"),
             ("SUP4D.12", "Lambda1_onto_Hom(S-,S+)_rank", "4"), ("SUP4D.12", "Lambda3_onto_Hom(S+,S-)_rank", "4"), ("SUP4D.12", "(J(x)J)conj(J(x)J)_is_I_16_real", "1")],
            "AGREE_SCOPED", "the identifications hold as projections (P+ X P-), the reading in which End(S+) = Lambda0 + Lambda2+ is meant; as subspaces Lambda1 meets Hom(S-, S+) in 0 (Lambda1_cap_Hom(S-,S+)_as_subspaces) and maps onto both Hom blocks with rank 4; the real form O + O of 8D.C2 is NAMED, not run")

    dsh, dims = byid["DIMS:4D"]
    dmap = [("n", "M4D.G2", None, "sqrt_of_ordered_pairs"), ("spinor_module", "M4D.G1", "size", None), ("clifford_dim", "M4D.G3", "rank_of_subset_products", None),
            ("generated_algebra", "M4D.G3", "dim", None), ("commutant", "M4D.G4", "dim", None), ("generators", "M4D.G2", None, "sqrt_of_ordered_pairs"),
            ("generator_size", "M4D.G1", "size", None), ("lambda2_small", "M4D.S14", "eigen_plus1_dim", None),
            ("lambda2_big", "M4D.S14", "eigen_minus1_dim", None), ("T_small", "SUP4D.17", "T_on_Lambda2+", None),
            ("T_big", "SUP4D.17", "T_on_Lambda2-", None), ("k", "SUP4D.17", "k_from_T_small_plus_1", None), ("stabilizer_of_eta", "M4D.S16", "kernel_dim", None),
            ("sequence_dims", "M4D.S6", "dims", None), ("symbol_ranks", None, None, None), ("square_coefficients", "M4D.G1", "coeffs_checked", None)]
    for key, mid, fld, lit in dmap:
        want = dims.get(key)
        if key == "symbol_ranks":
            have = lst([get("M4D.S4", "rank_sigma_d0"), get("M4D.S4", "rank_C")])
            mids = ["M4D.S4"]
        elif lit is not None:
            op = int(get("M4D.G2", "ordered_pairs"))
            have = str([t for t in range(1, 33) if t * t == op][0])
            mids = [mid]
        else:
            have = get(mid, fld)
            mids = [mid]
        if key in ("spinor_module", "generator_size") and have == "8x8":
            have = "8"
        cmp_row("DIMS:4D." + key, dsh, str(want), mids, [], "AGREE" if have == want else "DISAGREE_VALUE", "measured " + str(have))
    bsh, bott = byid["BOTT:4"]
    cmp_row("BOTT:4", bsh, "Cl(R^4,-nsq) = " + bott["cl_square_minus_1"] + ", dim " + bott["dim_minus"] + ", irreducible module " + bott["irreducible_module_minus"]
            + "; Cl(R^4,+nsq) = " + bott["cl_square_plus_1"] + ", dim " + bott["dim_plus"], ["M4D.G3", "M4D.G5", "SUP4D.13", "SUP4D.14"],
            [("SUP4D.14", "generated_algebra_dim", bott["dim_minus"]), ("SUP4D.14", "module_dim", bott["irreducible_module_minus"]),
             ("SUP4D.13", "generated_algebra_dim", bott["dim_plus"]), ("SUP4D.13", "commutant_dim", "4")],
            "AGREE_SCOPED", "dims measured; the iso types M_2(H) are read from the measured commutant H and center R (DERIVED), not from an explicit isomorphism")
    cmp_row("ADD:R1", byid["ADD:R1"], "raw projection: sigma^T sigma = c nsq I + (1-c) xi xi^T with c = 1:2 in 4D", ["M4D.S8"],
            [("M4D.S8", "c_num", "1"), ("M4D.S8", "c_den", "2"), ("M4D.S8", "mismatches", "0")], "AGREE", "2 sigma^T sigma = nsq I + xi xi^T as an integer polynomial identity")
    cmp_row("ADD:R2", byid["ADD:R2"], "4D: sigma^T sigma = nsq I exactly with the 't Hooft coefficient vector eta^a_mu_nu xi^mu a^nu, equal up to a stated sign convention to Im of the quaternion product; each eta^a has norm^2 2; not in an orthonormal basis of Lambda2+",
            ["M4D.S9", "M4D.S10", "M4D.S8", "SUP4D.16"],
            [("M4D.S10", "mismatches", "0"), ("M4D.S10", "c_num", "1"), ("M4D.S10", "c_den", "1"), ("SUP4D.16", "t_a_is_sum_k P_ak Im(a*conj(xi))_k", "1"),
             ("SUP4D.16", "eta_Gram", "[2,0,0; 0,2,0; 0,0,2]"), ("M4D.S8", "c_den", "2")],
            "AGREE", "the stated sign convention measured: with 4 the distinguished index, t = (c3, -c2, c1) = P Im(a conj(xi)), P = [0,0,1; 0,-1,0; 1,0,0]; with e1 = 1 distinguished, t' = Im(a conj(xi)) exactly; orthonormal coordinates give c = 1/2 (M4D.S8), coefficients s = t over 2 read as orthonormal give c = 1/4 (M4D.S11)")
    cmp_row("ADD:R5", byid["ADD:R5"], "4D Euler count 1 - 4 + 3 = 0; the symbol sequence is a complex whatever F_A is; the operator sequence a complex only on F+ = 0", ["M4D.S6", "M4D.OB6", "M4D.OB5"],
            [("M4D.S6", "euler", "0"), ("M4D.OB6", "symbol_independent_of_A", "1"), ("M4D.OB6", "composite_zero", "1")],
            "AGREE_SCOPED", "the only-if is shown for one A with F+ != 0 (M4D.OB5); the general statement is NAMED")
    cmp_row("ADD:R6", byid["ADD:R6"], "real Clifford gamma(xi)^2 = -nsq I; the rolled-up M and its doubled self-adjoint form square to +nsq", ["M4D.G1", "M4D.S7", "SUP4D.15"],
            [("M4D.G1", "mismatches", "0"), ("M4D.S7", "MTM_mismatches", "0"), ("SUP4D.15", "mismatches", "0"), ("SUP4D.15", "D_symmetric", "1")],
            "AGREE", "both signs reported, each in its own row")
    cmp_row("ADD:DISCLOSURE", byid["ADD:DISCLOSURE"], "a separate implementation that does not read GPT-6's receipt c28f8af6... before measuring", [], [],
            "NOT_MEASURED", "this script was written from the CONV rows alone; neither GPT-6's receipt nor any other agent's script was opened")
    psh, plant = byid["PLANTED:lambda2_plus_4d_planted_as_4"]
    cmp_row("PLANTED:lambda2_plus_4d", psh, "planted " + plant["planted_value"] + ", true " + plant["true_value"], ["M4D.S14"],
            [("M4D.S14", "eigen_plus1_dim", plant["true_value"])], "AGREE", "the measured dim Lambda2+ = 3 is the true value and would have rejected the plant 4")
    for item, mid, fld, conv in (("4D_sequence_dims", "M4D.S6", "dims", None), ("4D_ranks", "M4D.S4", None, None),
                                 ("4D_lambda2_plus", "M4D.S14", "eigen_plus1_dim", None), ("tensor_16", "M4D.TS2", "complex_rank_of_images_of_basis_tensors", None)):
        csh, cf = byid["COMPARE:" + item]
        if item == "4D_ranks":
            have = lst([get("M4D.S4", "rank_sigma_d0"), get("M4D.S4", "rank_C")])
        else:
            have = get(mid, fld)
        cmp_row("FROZEN_COMPARE:" + item, csh, cf["this_table"], [mid], [], "AGREE" if have == cf["this_table"] else "DISAGREE_VALUE", "measured " + str(have))
    emit("COMPARE_SUMMARY", [("compare_rows", CSTATS["rows"]), ("verdict_counts", "; ".join(k + " " + str(v) for k, v in sorted(VERDICTS.items()))),
                             ("disagreements", len(DISAGREE)), ("disagreement_refs", lst(DISAGREE) if DISAGREE else "none"),
                             ("measured_rows_pass", lst([STATS["meas_pass"], STATS["meas"]])), ("controls_failed_as_required", lst([STATS["ctrl_ok"], STATS["ctrl"]])),
                             ("supplement_rows_pass", lst([SSTATS["pass"], SSTATS["rows"]])),
                             ("open", "differential-form projections and gauge consequences beyond these symbol identities are not tested here")])


def write_out(name):
    body = "\n".join(LINES) + "\n"
    path = os.path.join(HERE, name)
    with open(path, "wb") as fh:
        fh.write(body.encode("utf-8"))
    h = hashlib.sha256(body.encode("utf-8")).hexdigest()
    with open(path + ".sha256", "wb") as fh:
        fh.write((h + "  " + name + "\n").encode("utf-8"))
    return h


def script_sha():
    with open(os.path.abspath(__file__), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


if __name__ == "__main__":
    main_measure()
    if MEASURE_ONLY:
        emit("SCRIPT", [("script", "measure_4d.py"), ("script_sha256", script_sha()), ("mode", "measure-only")])
        emit("END", [("rows", len(LINES) + 1)])
        print(write_out("MEASURE-4D-BEFORE-COMPARE.hbp"))
    else:
        main_supplement()
        compare_section()
        emit("SCRIPT", [("script", "measure_4d.py"), ("script_sha256", script_sha()), ("mode", "full"),
                        ("before_compare_script_sha256", "763fef99c27b43c7fa3cf56bc88941d4eae6058367f6ac4ef14efe8f814fc813")])
        emit("END", [("rows", len(LINES) + 1)])
        print(write_out("MEASURE-4D.hbp"))
    print("MEASURE_4D|meas=" + str(STATS["meas"]) + "|meas_pass=" + str(STATS["meas_pass"]) + "|controls=" + str(STATS["ctrl"])
          + "|controls_ok=" + str(STATS["ctrl_ok"]) + "|json=0")
