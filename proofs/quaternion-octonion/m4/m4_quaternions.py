# python3 -I -B m4_quaternions.py
# M4 quaternion step: the algebra H and its norm, H onto SU(2), the group Q8, the binary
# polyhedral groups 2T 2O 2I, their images in SO(3), binary dihedral groups inside the torus
# normaliser N(T), and the finite-subgroup ceiling of SU(2).
# Integer arithmetic only. A ring element is a tuple of integers:
#   Z         (a,)
#   Z[sqrt2]  (a,b) means a + b*sqrt2
#   Z[phi]    (a,b) means a + b*phi, with phi*phi equal to phi + 1
#   Z[zeta_m] coefficient vector of length deg(Phi_m), reduced modulo the cyclotomic polynomial Phi_m
# Half-integer and root coordinates are carried scaled by 2: a stored quaternion X stands for q with X equal to 2q.
# E=0: compute only, no network; writes only M4-QUATERNIONS.hbp and its .sha256 sidecar beside this file.

import sys
import os
import io
import hashlib
import tokenize
import itertools
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "M4-QUATERNIONS.hbp"
OUT = os.path.join(HERE, OUT_NAME)
SHA_OUT = os.path.join(HERE, OUT_NAME + ".sha256")

ROWS = []
STATE = {"checks": 0, "fail": 0, "failed_tags": [], "controls": 0}


# ------------------------------------------------------------------ rows
def clean(v):
    s = str(v)
    s = s.replace("|", ";").replace("=", ":").replace("\n", "_").replace("\r", "_")
    s = s.replace(" ", "_")
    return s


def emit(tag, **kv):
    parts = [tag]
    for k, v in kv.items():
        parts.append("%s=%s" % (k.rstrip("_"), clean(v)))
    parts.append("json=0")
    ROWS.append("|".join(parts))


def check(tag, ok, **kv):
    STATE["checks"] += 1
    if not ok:
        STATE["fail"] += 1
        STATE["failed_tags"].append(tag)
    d = {"result": "PASS" if ok else "FAIL"}
    d.update(kv)
    d["claim"] = "MEASURED"
    emit(tag, **d)
    return ok


def control(tag, observed_ok, expected_ok, **kv):
    # a control passes when the observed outcome equals the expected outcome
    ok = observed_ok == expected_ok
    STATE["checks"] += 1
    STATE["controls"] += 1
    if not ok:
        STATE["fail"] += 1
        STATE["failed_tags"].append(tag)
    d = {
        "result": "PASS" if ok else "FAIL",
        "observed": "PASS" if observed_ok else "FAIL",
        "expected": "PASS" if expected_ok else "FAIL",
    }
    d.update(kv)
    d["claim"] = "MEASURED"
    emit(tag, **d)
    return ok


def proof(tag, claim, statement, argument, checked, relies_on="none"):
    emit(tag, statement=statement, argument=argument, checked=checked, relies_on=relies_on, claim=claim)


def named(tag, statement, ref, **kv):
    d = {"statement": statement, "ref": ref}
    d.update(kv)
    d["citation_opened_in_this_run"] = 0
    d["claim"] = "NAMED"
    emit(tag, **d)


def hist_text(h):
    return ",".join("%d:%d" % (k, h[k]) for k in sorted(h))


# ------------------------------------------------------------------ self lint
ALLOWED_MATH_LINES = {"from math import gcd", "from math import isqrt", "from math import gcd, isqrt"}


def self_lint():
    return lint_source(open(os.path.abspath(__file__), "rb").read())


def lint_source(src):
    bad = []
    sl = chr(47)
    banned_ops = {sl, sl + sl, sl + "=", sl + sl + "="}
    banned_names = {"float", "Fraction", "Decimal", "complex", "fractions", "decimal", "cmath",
                    "exec", "eval", "pickle", "marshal", "__import__", "numpy", "random"}
    for t in tokenize.tokenize(io.BytesIO(src).readline):
        if t.type == tokenize.OP and t.string in banned_ops:
            bad.append(("op", t.start[0]))
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if "." in s or "e" in s or "j" in s or s.startswith("0x"):
                bad.append(("number", t.start[0]))
        elif t.type == tokenize.NAME and t.string in banned_names:
            bad.append(("name", t.start[0]))
        elif t.type == tokenize.NAME and t.string == "math":
            if t.line.strip() not in ALLOWED_MATH_LINES:
                bad.append(("math", t.start[0]))
    slash_chars = src.count(b"\x2f")
    return bad, slash_chars


# ------------------------------------------------------------------ integer polynomials
# a polynomial is a dict: exponent tuple -> nonzero integer coefficient
def p_clean(A):
    return {k: v for k, v in A.items() if v != 0}


def p_add(A, B):
    C = dict(A)
    for k, v in B.items():
        C[k] = C.get(k, 0) + v
    return p_clean(C)


def p_neg(A):
    return {k: -v for k, v in A.items()}


def p_sub(A, B):
    return p_add(A, p_neg(B))


def p_mul(A, B):
    C = {}
    for ka, va in A.items():
        for kb, vb in B.items():
            kk = tuple(x + y for x, y in zip(ka, kb))
            C[kk] = C.get(kk, 0) + va * vb
    return p_clean(C)


def p_scale(c, A):
    return p_clean({k: c * v for k, v in A.items()})


def p_var(t, n):
    e = [0] * n
    e[t] = 1
    return {tuple(e): 1}


def p_const(c, n):
    return p_clean({tuple([0] * n): c})


# Gaussian-integer polynomials: pairs (Re, Im) of integer polynomials
def gp_add(x, y):
    return (p_add(x[0], y[0]), p_add(x[1], y[1]))


def gp_sub(x, y):
    return (p_sub(x[0], y[0]), p_sub(x[1], y[1]))


def gp_neg(x):
    return (p_neg(x[0]), p_neg(x[1]))


def gp_mul(x, y):
    return (p_sub(p_mul(x[0], y[0]), p_mul(x[1], y[1])), p_add(p_mul(x[0], y[1]), p_mul(x[1], y[0])))


def gp_conj(x):
    return (x[0], p_neg(x[1]))


def gpm_mul(A, B):
    return tuple(tuple(gp_add(gp_mul(A[r][0], B[0][c]), gp_mul(A[r][1], B[1][c])) for c in range(2)) for r in range(2))


def gpm_dagger(A):
    return ((gp_conj(A[0][0]), gp_conj(A[1][0])), (gp_conj(A[0][1]), gp_conj(A[1][1])))


# Gaussian integers as pairs
def g_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def g_mul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gm_mul(A, B):
    return tuple(tuple(g_add(g_mul(A[r][0], B[0][c]), g_mul(A[r][1], B[1][c])) for c in range(2)) for r in range(2))


# ------------------------------------------------------------------ quaternion table
NAMES = "1ijk"


def build_qtable(flip=None):
    T = [[None] * 4 for _ in range(4)]
    for a in range(4):
        T[0][a] = (1, a)
        T[a][0] = (1, a)
    for a in range(1, 4):
        T[a][a] = (-1, 0)
    for a, b, c in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        T[a][b] = (1, c)
        T[b][a] = (-1, c)
    if flip is not None:
        a, b = flip
        s, c = T[a][b]
        T[a][b] = (-s, c)
    return T


def table_text(T):
    parts = []
    for a in range(1, 4):
        for b in range(1, 4):
            s, c = T[a][b]
            parts.append("%s%s%s%s" % (NAMES[a], NAMES[b], "+" if s > 0 else "-", NAMES[c]))
    return ",".join(parts)


def unit_int(a):
    return tuple(1 if t == a else 0 for t in range(4))


def qi_mul(T, x, y):
    out = [0, 0, 0, 0]
    for a in range(4):
        if x[a]:
            for b in range(4):
                if y[b]:
                    s, c = T[a][b]
                    out[c] += s * x[a] * y[b]
    return tuple(out)


def qi_conj(x):
    return (x[0], -x[1], -x[2], -x[3])


def assoc_failures(T):
    E = [unit_int(a) for a in range(4)]
    f = 0
    for a in range(4):
        for b in range(4):
            for c in range(4):
                if qi_mul(T, qi_mul(T, E[a], E[b]), E[c]) != qi_mul(T, E[a], qi_mul(T, E[b], E[c])):
                    f += 1
    return f


def qp_mul(T, x, y):
    out = [{}, {}, {}, {}]
    for a in range(4):
        for b in range(4):
            s, c = T[a][b]
            pr = p_mul(x[a], y[b])
            out[c] = p_add(out[c], pr) if s > 0 else p_sub(out[c], pr)
    return out


def norm_mult_poly(T, signs):
    n = 8
    P = [p_var(t, n) for t in range(4)]
    Q = [p_var(4 + t, n) for t in range(4)]
    PQ = qp_mul(T, P, Q)

    def nform(X):
        acc = {}
        for t in range(4):
            sq = p_mul(X[t], X[t])
            acc = p_add(acc, sq) if signs[t] > 0 else p_sub(acc, sq)
        return acc

    NPQ = nform(PQ)
    NPNQ = p_mul(nform(P), nform(Q))
    diff = p_sub(NPQ, NPNQ)
    raw = set()
    for t in range(4):
        raw.update(p_mul(PQ[t], PQ[t]).keys())
    return NPQ, NPNQ, diff, raw


def phi_int(q):
    a, b, c, d = q
    return (((a, b), (c, d)), ((-c, d), (a, -b)))


def phi_int_bad(q):
    # broken control: lower-left entry +c instead of -c
    a, b, c, d = q
    return (((a, b), (c, d)), ((c, d), (a, -b)))


def phi_poly(q):
    return (((q[0], q[1]), (q[2], q[3])), ((p_neg(q[2]), q[3]), (q[0], p_neg(q[1]))))


# ------------------------------------------------------------------ rings
def poly_divmod_monic(num, den):
    num = list(num)
    dd = len(den) - 1
    if den[dd] != 1:
        raise ArithmeticError("divisor not monic")
    if len(num) <= dd:
        return [0], num + [0] * (dd - len(num))
    quo = [0] * (len(num) - dd)
    for k in range(len(num) - 1, dd - 1, -1):
        c = num[k]
        if c:
            quo[k - dd] = c
            for t in range(dd + 1):
                num[k - dd + t] -= c * den[t]
    return quo, num[:dd]


CYC = {}


def cyclotomic(m):
    if m in CYC:
        return CYC[m]
    num = [-1] + [0] * (m - 1) + [1]
    for d in range(1, m):
        if m % d == 0:
            q, r = poly_divmod_monic(num, cyclotomic(d))
            if any(r):
                raise ArithmeticError("cyclotomic division inexact")
            num = q
    while len(num) > 1 and num[-1] == 0:
        num.pop()
    CYC[m] = num
    return num


RZ = {"name": "Z", "kind": "int", "dim": 1}
RS2 = {"name": "Z[sqrt2]", "kind": "quad", "dim": 2, "p": 2, "q": 0}
RPH = {"name": "Z[phi]", "kind": "quad", "dim": 2, "p": 1, "q": 1}


def r_make_cyclo(m):
    P = cyclotomic(m)
    return {"name": "Z[zeta_%d]" % m, "kind": "cyclo", "dim": len(P) - 1, "m": m, "Phi": P}


def r_zero(R):
    return (0,) * R["dim"]


def r_int(R, n):
    return (n,) + (0,) * (R["dim"] - 1)


def r_add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def r_sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def r_neg(u):
    return tuple(-a for a in u)


def r_iszero(u):
    return all(a == 0 for a in u)


def r_half(u):
    out = []
    for a in u:
        h, r = divmod(a, 2)
        if r:
            return None
        out.append(h)
    return tuple(out)


def cyc_reduce(R, poly):
    D = R["dim"]
    _, rem = poly_divmod_monic(poly, R["Phi"])
    rem = list(rem) + [0] * (D - len(rem))
    return tuple(rem[:D])


def r_mul(R, u, v):
    k = R["kind"]
    if k == "int":
        return (u[0] * v[0],)
    if k == "quad":
        t = u[1] * v[1]
        return (u[0] * v[0] + R["p"] * t, u[0] * v[1] + u[1] * v[0] + R["q"] * t)
    D = R["dim"]
    prod = [0] * (2 * D - 1)
    for a in range(D):
        if u[a]:
            for b in range(D):
                if v[b]:
                    prod[a + b] += u[a] * v[b]
    return cyc_reduce(R, prod)


def r_conj(R, u):
    # complex conjugation on Z[zeta_m]: zeta^k goes to zeta^(m-k)
    m = R["m"]
    poly = [0] * m
    for k, c in enumerate(u):
        poly[(m - k) % m] += c
    return cyc_reduce(R, poly)


# quaternions with ring coordinates
def q_from_ints(R, v):
    return tuple(r_int(R, a) for a in v)


def q_mul(R, T, x, y):
    out = [r_zero(R)] * 4
    for a in range(4):
        if r_iszero(x[a]):
            continue
        for b in range(4):
            if r_iszero(y[b]):
                continue
            s, c = T[a][b]
            pr = r_mul(R, x[a], y[b])
            out[c] = r_add(out[c], pr) if s > 0 else r_sub(out[c], pr)
    return tuple(out)


def q_conj(x):
    return (x[0], r_neg(x[1]), r_neg(x[2]), r_neg(x[3]))


def q_sp(R, T, x, y):
    # scaled product: x stands for 2p, y for 2q; returns the stored form 2pq, or None when not in the ring
    z = q_mul(R, T, x, y)
    out = []
    for c in z:
        h = r_half(c)
        if h is None:
            return None
        out.append(h)
    return tuple(out)


def q_norm(R, x):
    s = r_zero(R)
    for c in x:
        s = r_add(s, r_mul(R, c, c))
    return s


def emb(R, X):
    return tuple(r_int(R, c[0]) for c in X)


def perm_parity(p):
    inv = 0
    for a in range(len(p)):
        for b in range(a + 1, len(p)):
            if p[a] > p[b]:
                inv += 1
    return inv % 2


def icosian96(parity):
    base = [(0, 0), (1, 0), (-1, 1), (0, 1)]   # 0, 1, phi^-1 as phi-1, phi (all scaled by 2)
    out = []
    for perm in itertools.permutations(range(4)):
        if perm_parity(perm) != parity:
            continue
        for s1, s2, s3 in itertools.product((1, -1), repeat=3):
            sv = [base[0], (s1 * base[1][0], s1 * base[1][1]), (s2 * base[2][0], s2 * base[2][1]),
                  (s3 * base[3][0], s3 * base[3][1])]
            w = [None] * 4
            for t in range(4):
                w[perm[t]] = sv[t]
            out.append(tuple(w))
    return out


def analyze(R, T, G):
    n = len(G)
    one = q_from_ints(R, [2, 0, 0, 0])
    four = r_int(R, 4)
    info = {"n": n, "distinct": len(set(G)) == n}
    info["norms_ok"] = all(q_norm(R, x) == four for x in G)
    idx = {g: t for t, g in enumerate(G)}
    tab = []
    clo = 0
    for x in G:
        row = []
        for y in G:
            z = q_sp(R, T, x, y)
            t = idx.get(z, -1) if z is not None else -1
            if t < 0:
                clo += 1
            row.append(t)
        tab.append(row)
    info["closure_fail"] = clo
    info["has_one"] = one in idx
    info["tab"] = tab
    info["idx"] = idx
    if not info["has_one"]:
        info["inv_ok"] = False
        info["ident_ok"] = False
        return info
    e = idx[one]
    info["e"] = e
    info["ident_ok"] = all(tab[e][t] == t and tab[t][e] == t for t in range(n))
    inv_ok = True
    for t in range(n):
        c = q_conj(G[t])
        if c not in idx or tab[t][idx[c]] != e or tab[idx[c]][t] != e:
            inv_ok = False
    info["inv_ok"] = inv_ok
    if clo == 0:
        orders = []
        for g in range(n):
            cur = g
            k = 1
            while cur != e and k <= n:
                cur = tab[cur][g]
                k += 1
            orders.append(k)
        hist = {}
        for o in orders:
            hist[o] = hist.get(o, 0) + 1
        info["orders"] = orders
        info["hist"] = hist
        center = [g for g in range(n) if all(tab[g][h] == tab[h][g] for h in range(n))]
        info["center"] = center
        cmax = 0
        for g in range(n):
            if g in center:
                continue
            cs = sum(1 for h in range(n) if tab[g][h] == tab[h][g])
            if cs > cmax:
                cmax = cs
        info["cent_max_noncentral"] = cmax
        info["involutions"] = [g for g in range(n) if orders[g] == 2]
    return info


def rot_matrix(R, T, X, use_conj=True):
    Xc = q_conj(X) if use_conj else X
    M = [[None] * 3 for _ in range(3)]
    imag_ok = True
    for c in range(3):
        e = tuple(r_int(R, 1) if t == c + 1 else r_zero(R) for t in range(4))
        w = q_mul(R, T, q_mul(R, T, X, e), Xc)
        if not r_iszero(w[0]):
            imag_ok = False
        for r in range(3):
            M[r][c] = w[r + 1]
    return tuple(tuple(row) for row in M), imag_ok


def m3_mul(R, A, B):
    out = []
    for r in range(3):
        row = []
        for c in range(3):
            s = r_zero(R)
            for t in range(3):
                s = r_add(s, r_mul(R, A[r][t], B[t][c]))
            row.append(s)
        out.append(tuple(row))
    return tuple(out)


def m3_scale_int(R, k, A):
    return tuple(tuple(r_mul(R, r_int(R, k), A[r][c]) for c in range(3)) for r in range(3))


def m3_vec(R, A, v):
    return tuple(sum_r(R, [r_mul(R, A[r][c], v[c]) for c in range(3)]) for r in range(3))


def sum_r(R, xs):
    s = r_zero(R)
    for x in xs:
        s = r_add(s, x)
    return s


def cross_r(R, u, v):
    return (r_sub(r_mul(R, u[1], v[2]), r_mul(R, u[2], v[1])),
            r_sub(r_mul(R, u[2], v[0]), r_mul(R, u[0], v[2])),
            r_sub(r_mul(R, u[0], v[1]), r_mul(R, u[1], v[0])))


def det3(R, M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    t1 = r_mul(R, a, r_sub(r_mul(R, e, i), r_mul(R, f, h)))
    t2 = r_mul(R, b, r_sub(r_mul(R, d, i), r_mul(R, f, g)))
    t3 = r_mul(R, c, r_sub(r_mul(R, d, h), r_mul(R, e, g)))
    return r_add(r_sub(t1, t2), t3)


def orth_ok(R, M, scale_sq):
    z = r_zero(R)
    target = r_int(R, scale_sq)
    for a in range(3):
        for b in range(3):
            s = sum_r(R, [r_mul(R, M[r][a], M[r][b]) for r in range(3)])
            if s != (target if a == b else z):
                return False
    return True


def rot_analysis(R, T, G, info):
    four = r_int(R, 4)
    z = r_zero(R)
    mats = []
    imag_all = True
    orth_all = True
    det_all = True
    kernel = []
    for X in G:
        M, im = rot_matrix(R, T, X)
        imag_all = imag_all and im
        if not orth_ok(R, M, 16):
            orth_all = False
        if det3(R, M) != r_int(R, 64):
            det_all = False
        if all(M[r][c] == (four if r == c else z) for r in range(3) for c in range(3)):
            kernel.append(X)
        mats.append(M)
    hom_fail = 0
    tab = info["tab"]
    n = len(G)
    for x in range(n):
        for y in range(n):
            if m3_mul(R, mats[x], mats[y]) != m3_scale_int(R, 4, mats[tab[x][y]]):
                hom_fail += 1
    # rotations compatible with the cross product: (M'u) x (M'v) equals 4 M'(u x v)
    cross_fail = 0
    basis3 = [tuple(r_int(R, 1) if t == c else z for t in range(3)) for c in range(3)]
    for M in set(mats):
        for u in basis3:
            for v in basis3:
                lhs = cross_r(R, m3_vec(R, M, u), m3_vec(R, M, v))
                rhs = tuple(r_mul(R, four, x) for x in m3_vec(R, M, cross_r(R, u, v)))
                if lhs != rhs:
                    cross_fail += 1
    return {"mats": mats, "image": set(mats), "imag": imag_all, "orth": orth_all, "det": det_all,
            "kernel": kernel, "hom_fail": hom_fail, "cross_fail": cross_fail}


DIAG = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]


def quarter_int_matrix(R, M):
    Mi = []
    for row in M:
        rr = []
        for e in row:
            h = r_half(e)
            h = r_half(h) if h is not None else None
            if h is None or any(h[1:]):
                return None
            rr.append(h[0])
        Mi.append(tuple(rr))
    return tuple(Mi)


def diag_perm(Mi):
    perm = []
    for d in DIAG:
        v = tuple(sum(Mi[r][c] * d[c] for c in range(3)) for r in range(3))
        hit = -1
        for t, d2 in enumerate(DIAG):
            if v == d2 or v == tuple(-x for x in d2):
                hit = t
        if hit < 0:
            return None
        perm.append(hit)
    return tuple(perm)


def int_m3_mul(A, B):
    return tuple(tuple(sum(A[r][t] * B[t][c] for t in range(3)) for c in range(3)) for r in range(3))


def compose(p, q):
    # apply q first, then p
    return tuple(p[q[t]] for t in range(len(q)))


# Cayley-Dickson form q equals z + w j, z and w in a commutative ring with conjugation
def cd_mul(R, x, y):
    z1, w1 = x
    z2, w2 = y
    return (r_sub(r_mul(R, z1, z2), r_mul(R, w1, r_conj(R, w2))),
            r_add(r_mul(R, z1, w2), r_mul(R, w1, r_conj(R, z2))))


def cd_generate(R, gens, cap=4000):
    one = (r_int(R, 1), r_zero(R))
    S = {one}
    frontier = [one]
    while frontier and len(S) <= cap:
        nxt = []
        for a in frontier:
            for g in gens:
                b = cd_mul(R, a, g)
                if b not in S:
                    S.add(b)
                    nxt.append(b)
        frontier = nxt
    return S


def generate_int(T, gens, cap=4000):
    one = (1, 0, 0, 0)
    S = {one}
    frontier = [one]
    while frontier and len(S) <= cap:
        nxt = []
        for a in frontier:
            for g in gens:
                b = qi_mul(T, a, g)
                if b not in S:
                    S.add(b)
                    nxt.append(b)
        frontier = nxt
    return S


def closure_failures(R, T, G):
    S = set(G)
    f = 0
    for x in G:
        for y in G:
            z = q_sp(R, T, x, y)
            if z is None or z not in S:
                f += 1
    return f


# ------------------------------------------------------------------ main
def main():
    bad, slash = self_lint()
    check("LINT_INTEGER_ONLY", len(bad) == 0 and slash == 0, bad_tokens=len(bad), slash_chars_in_source=slash,
          banned="slash_ops,float_literals,float_fraction_decimal_complex_names,exec_eval_pickle_marshal,math_beyond_gcd_isqrt")
    sample = ("x = 7 " + chr(47) + " 2\ny = 0.5\nz = math.sqrt(2)\nw = eval(s)\n").encode()
    sbad, _ = lint_source(sample)
    control("CONTROL_LINT_CATCHES_BAD_SOURCE", len(sbad) == 0, False, sample="true_division+float_literal+math_sqrt+eval",
            flagged_tokens=len(sbad))
    src = open(os.path.abspath(__file__), "rb").read()
    emit("SCRIPT", name="m4_quaternions.py", sha256=hashlib.sha256(src).hexdigest(), python=sys.version.split()[0],
         E=0, network=0, float_used=0, claim="MEASURED")
    emit("CONVENTIONS", basis="1,i,j,k_as_indices_0..3", table="ij:k,jk:i,ki:j,ji:-k,kj:-i,ik:-j,ii:jj:kk:-1",
         phi_map="a+bi+cj+dk_to_[[a+bi,c+di],[-c+di,a-bi]]", cayley_dickson="q:z+wj_with_jz:conj(z)j",
         scaling="binary_groups_stored_as_2q", rings="Z:(a);Z[sqrt2]:(a,b)_a+b*sqrt2;Z[phi]:(a,b)_a+b*phi_phi^2:phi+1;Z[zeta_m]:coeffs_mod_Phi_m",
         claim="NAMED")

    # ================================================================ 1. H
    T = build_qtable()
    emit("H_TABLE", products=table_text(T), claim="MEASURED")
    E4 = [unit_int(a) for a in range(4)]
    ONE = (1, 0, 0, 0)
    MONE = (-1, 0, 0, 0)
    I_, J_, K_ = E4[1], E4[2], E4[3]
    rel = (qi_mul(T, I_, I_) == MONE and qi_mul(T, J_, J_) == MONE and qi_mul(T, K_, K_) == MONE
           and qi_mul(T, qi_mul(T, I_, J_), K_) == MONE)
    check("H_RELATIONS_I2_J2_K2_IJK_MINUS1", rel, i2=qi_mul(T, I_, I_), j2=qi_mul(T, J_, J_), k2=qi_mul(T, K_, K_),
          ijk=qi_mul(T, qi_mul(T, I_, J_), K_))
    af = assoc_failures(T)
    check("H_ASSOC_BASIS", af == 0, basis_triples=64, failures=af)
    proof("PROOF_H_ASSOCIATIVE", "MEASURED",
          "(xy)z equals x(yz) for all x,y,z in H with coefficients in any commutative ring",
          "the product is bilinear with integer structure constants so the associator (xy)z-x(yz) is trilinear; "
          "a trilinear map that vanishes on every triple of basis vectors vanishes identically",
          "H_ASSOC_BASIS: all 64 basis triples, 0 failures")
    NPQ, NPNQ, diff, raw = norm_mult_poly(T, (1, 1, 1, 1))
    check("H_NORM_MULTIPLICATIVE_POLY", diff == {} and len(NPQ) == 16 and len(NPNQ) == 16,
          variables=8, monomials_N_pq=len(NPQ), monomials_Np_Nq=len(NPNQ), raw_square_monomials=len(raw),
          cancelled_cross_monomials=len(raw) - len(NPQ), difference_terms=len(diff))
    proof("PROOF_H_NORM_MULTIPLICATIVE", "MEASURED",
          "N(pq) equals N(p)N(q) for all p,q, N(q) the sum of the four coordinate squares",
          "N(pq)-N(p)N(q) was expanded as an integer polynomial in 8 indeterminates a0..a3,b0..b3 and every "
          "coefficient is zero; an identity in Z[a,b] specialises to every commutative ring (R, Z[sqrt2], Z[phi], Z[zeta_m])",
          "H_NORM_MULTIPLICATIVE_POLY: full symbolic expansion, no sampling")
    n4 = 4
    Qv = [p_var(t, n4) for t in range(4)]
    Qc = [Qv[0], p_neg(Qv[1]), p_neg(Qv[2]), p_neg(Qv[3])]
    Nq = {}
    for t in range(4):
        Nq = p_add(Nq, p_mul(Qv[t], Qv[t]))
    left = qp_mul(T, Qv, Qc)
    right = qp_mul(T, Qc, Qv)
    check("H_Q_CONJQ_EQUALS_N_POLY", left == [Nq, {}, {}, {}] and right == [Nq, {}, {}, {}],
          identity="q*conj(q):conj(q)*q:N(q)")
    U6 = [p_var(t, 6) for t in range(6)]
    zero6 = {}
    u = [zero6, U6[0], U6[1], U6[2]]
    v = [zero6, U6[3], U6[4], U6[5]]
    uv = qp_mul(T, u, v)
    dotuv = p_add(p_add(p_mul(U6[0], U6[3]), p_mul(U6[1], U6[4])), p_mul(U6[2], U6[5]))
    crs = [p_sub(p_mul(U6[1], U6[5]), p_mul(U6[2], U6[4])), p_sub(p_mul(U6[2], U6[3]), p_mul(U6[0], U6[5])),
           p_sub(p_mul(U6[0], U6[4]), p_mul(U6[1], U6[3]))]
    check("H_IMAG_PRODUCT_POLY", uv == [p_neg(dotuv), crs[0], crs[1], crs[2]],
          identity="uv:-(u.v)+(u_cross_v)_for_imaginary_u,v")

    # ================================================================ 2. phi into 2x2 over Z[i]
    hom_fail = 0
    for a in range(4):
        for b in range(4):
            if phi_int(qi_mul(T, E4[a], E4[b])) != gm_mul(phi_int(E4[a]), phi_int(E4[b])):
                hom_fail += 1
    check("PHI_HOM_BASIS", hom_fail == 0, basis_pairs=16, failures=hom_fail)
    check("PHI_ONE_IS_I", phi_int(ONE) == (((1, 0), (0, 0)), ((0, 0), (1, 0))), phi_1=phi_int(ONE))
    P8 = [p_var(t, 8) for t in range(4)]
    Q8v = [p_var(4 + t, 8) for t in range(4)]
    hom_poly = phi_poly(qp_mul(T, P8, Q8v)) == gpm_mul(phi_poly(P8), phi_poly(Q8v))
    check("PHI_HOM_POLY", hom_poly, variables=8, identity="phi(pq):phi(p)phi(q)_symbolic")
    E = phi_poly(Qv)
    detE = gp_sub(gp_mul(E[0][0], E[1][1]), gp_mul(E[0][1], E[1][0]))
    check("PHI_DET_EQUALS_N_POLY", detE == (Nq, {}), identity="det_phi(q):a^2+b^2+c^2+d^2", imag_part_terms=len(detE[1]))
    check("PHI_DAGGER_EQUALS_PHI_CONJ_POLY", gpm_dagger(phi_poly(Qv)) == phi_poly(Qc),
          identity="phi(q)^dagger:phi(conj_q)")
    rec_ok = all((phi_int(E4[a])[0][0][0], phi_int(E4[a])[0][0][1], phi_int(E4[a])[0][1][0], phi_int(E4[a])[0][1][1]) == E4[a]
                 for a in range(4))
    check("PHI_INJECTIVE", rec_ok, left_inverse="q:(Re_E00,Im_E00,Re_E01,Im_E01)", checked_on="4_basis_vectors_then_linearity")
    X8 = [p_var(t, 8) for t in range(8)]
    Um = (((X8[0], X8[1]), (X8[2], X8[3])), ((X8[4], X8[5]), (X8[6], X8[7])))
    adj = ((Um[1][1], gp_neg(Um[0][1])), (gp_neg(Um[1][0]), Um[0][0]))
    detU = gp_sub(gp_mul(Um[0][0], Um[1][1]), gp_mul(Um[0][1], Um[1][0]))
    zg = ({}, {})
    ok_adj = gpm_mul(Um, adj) == ((detU, zg), (zg, detU)) and gpm_mul(adj, Um) == ((detU, zg), (zg, detU))
    check("ADJUGATE_IDENTITY_POLY", ok_adj, variables=8, identity="U*adj(U):adj(U)*U:det(U)*I_for_general_complex_2x2")
    proof("PROOF_PHI_HOMOMORPHISM", "MEASURED",
          "phi(pq) equals phi(p)phi(q) and phi(1) equals I for all quaternions p,q",
          "both sides are R-bilinear in (p,q); they agree on the 16 basis pairs, hence everywhere; "
          "the 8-variable symbolic expansion confirms it independently",
          "PHI_HOM_BASIS 16 pairs 0 failures; PHI_HOM_POLY; PHI_ONE_IS_I")
    proof("PROOF_S3_ISOMORPHIC_SU2", "MEASURED",
          "phi restricts to a group isomorphism from the unit quaternions S^3 onto SU(2)",
          "(a) N(q):1 gives phi(q)^dagger phi(q):phi(conj q)phi(q):phi(conj(q)q):phi(N(q)):I and det phi(q):N(q):1, so phi(S^3) lies in SU(2). "
          "(b) let U:[[al,be],[ga,de]] be in SU(2); det U:1 and U*adj(U):det(U)I give U^-1:adj(U):[[de,-be],[-ga,al]]; unitarity gives "
          "U^-1:U^dagger:[[conj al,conj ga],[conj be,conj de]]; comparing entries de:conj(al), ga:-conj(be), and det U:|al|^2+|be|^2:1; "
          "with al:a+bi, be:c+di this is U:phi(a+bi+cj+dk) with N:1, so phi(S^3) is all of SU(2). "
          "(c) phi is injective (left inverse reads entries 00 and 01). (d) phi is multiplicative. "
          "Hence a bijective homomorphism; phi is linear so it is a homeomorphism",
          "PHI_HOM_BASIS, PHI_HOM_POLY, PHI_DET_EQUALS_N_POLY, PHI_DAGGER_EQUALS_PHI_CONJ_POLY, H_Q_CONJQ_EQUALS_N_POLY, "
          "ADJUGATE_IDENTITY_POLY, PHI_INJECTIVE")

    # ================================================================ 3. Q8
    Q8 = [tuple(s * x for x in E4[a]) for a in range(4) for s in (1, -1)]
    Q8set = set(Q8)
    clo = sum(1 for x in Q8 for y in Q8 if qi_mul(T, x, y) not in Q8set)
    ident = ONE in Q8set and all(qi_mul(T, ONE, x) == x and qi_mul(T, x, ONE) == x for x in Q8)
    inv = all(qi_mul(T, x, qi_conj(x)) == ONE and qi_mul(T, qi_conj(x), x) == ONE and qi_conj(x) in Q8set for x in Q8)
    assoc = sum(1 for x in Q8 for y in Q8 for z in Q8 if qi_mul(T, qi_mul(T, x, y), z) != qi_mul(T, x, qi_mul(T, y, z)))
    noncomm = sum(1 for x in Q8 for y in Q8 if qi_mul(T, x, y) != qi_mul(T, y, x))

    def order_int(x):
        cur = x
        k = 1
        while cur != ONE and k <= 64:
            cur = qi_mul(T, cur, x)
            k += 1
        return k

    ords = {x: order_int(x) for x in Q8}
    h8 = {}
    for x in Q8:
        h8[ords[x]] = h8.get(ords[x], 0) + 1
    inv2 = [x for x in Q8 if ords[x] == 2]
    gen_ij = generate_int(T, [I_, J_])
    check("Q8_ORDER", len(Q8set) == 8, order=len(Q8set))
    check("Q8_CLOSURE", clo == 0, products=64, outside=clo)
    check("Q8_IDENTITY", ident)
    check("Q8_INVERSES", inv, inverse="conj(x)")
    check("Q8_ASSOC_ALL_TRIPLES", assoc == 0, triples=512, failures=assoc)
    check("Q8_NON_ABELIAN", noncomm > 0, noncommuting_ordered_pairs=noncomm, ij=qi_mul(T, I_, J_), ji=qi_mul(T, J_, I_))
    check("Q8_ORDERS_DIVIDE_4", all(4 % o == 0 for o in ords.values()), order_hist=hist_text(h8))
    check("Q8_UNIQUE_INVOLUTION", inv2 == [MONE], involutions=len(inv2), element=inv2[0] if inv2 else "none")
    check("Q8_GENERATED_BY_I_J", gen_ij == Q8set, generated_size=len(gen_ij))

    # ================================================================ 4. the binary polyhedral groups
    rel_ok = r_mul(RS2, (0, 1), (0, 1)) == (2, 0) and r_mul(RPH, (0, 1), (0, 1)) == (1, 1) and r_mul(RPH, (0, 1), (-1, 1)) == (1, 0)
    check("RING_RELATIONS", rel_ok, sqrt2_sq=r_mul(RS2, (0, 1), (0, 1)), phi_sq=r_mul(RPH, (0, 1), (0, 1)),
          phi_times_phi_minus_1=r_mul(RPH, (0, 1), (-1, 1)))
    Q8s = [q_from_ints(RZ, [2 * s if t == a else 0 for t in range(4)]) for a in range(4) for s in (1, -1)]
    T24 = Q8s + [q_from_ints(RZ, list(sg)) for sg in itertools.product((1, -1), repeat=4)]
    O48 = [emb(RS2, x) for x in T24]
    for x, y in itertools.combinations(range(4), 2):
        for s1 in (1, -1):
            for s2 in (1, -1):
                w = [r_zero(RS2)] * 4
                w[x] = (0, s1)
                w[y] = (0, s2)
                O48.append(tuple(w))
    ico_even = icosian96(0)
    ico_odd = icosian96(1)
    I120 = [emb(RPH, x) for x in T24] + ico_even

    groups = [("Q8", RZ, Q8s, 8), ("2T", RZ, T24, 24), ("2O", RS2, O48, 48), ("2I", RPH, I120, 120)]
    expect_hist = {
        "Q8": {1: 1, 2: 1, 4: 6},
        "2T": {1: 1, 2: 1, 3: 8, 4: 6, 6: 8},
        "2O": {1: 1, 2: 1, 3: 8, 4: 18, 6: 8, 8: 12},
        "2I": {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24},
    }
    expect_cmax = {"Q8": 4, "2T": 6, "2O": 8, "2I": 10}
    INFO = {}
    for name, R, G, order in groups:
        info = analyze(R, T, G)
        INFO[name] = info
        okg = (info["distinct"] and len(G) == order and info["norms_ok"] and info["closure_fail"] == 0
               and info["ident_ok"] and info["inv_ok"])
        check("BIN_%s_GROUP" % name, okg, ring=R["name"], order=len(G), expect=order, distinct=info["distinct"],
              unit_norm=info["norms_ok"], products=len(G) * len(G), closure_outside=info["closure_fail"],
              identity=info["ident_ok"], inverses_conj=info["inv_ok"])
        if info["closure_fail"] == 0:
            ctr = [G[g] for g in info["center"]]
            check("BIN_%s_ORDER_HIST" % name, info["hist"] == expect_hist[name], order_hist=hist_text(info["hist"]))
            check("BIN_%s_CENTER_PM1" % name, len(ctr) == 2 and q_from_ints(R, [2, 0, 0, 0]) in ctr and q_from_ints(R, [-2, 0, 0, 0]) in ctr,
                  center_size=len(ctr))
            check("BIN_%s_UNIQUE_INVOLUTION" % name, len(info["involutions"]) == 1 and G[info["involutions"][0]] == q_from_ints(R, [-2, 0, 0, 0]),
                  involutions=len(info["involutions"]))
            check("BIN_%s_MAX_NONCENTRAL_CENTRALIZER" % name, info["cent_max_noncentral"] == expect_cmax[name],
                  value=info["cent_max_noncentral"], expect=expect_cmax[name])
    proof("PROOF_FINITE_CLOSED_SUBSET_IS_SUBGROUP", "MEASURED",
          "each of Q8, 2T, 2O, 2I is a subgroup of the unit quaternions S^3",
          "unit quaternions form a group (associative by PROOF_H_ASSOCIATIVE, N multiplicative so closed, conj(q) is the inverse); "
          "a finite nonempty subset S of a group closed under the product is a subgroup: left multiplication by s is injective (cancellation), "
          "hence a bijection of the finite S, so s*x:s and s*y:1 have solutions in S; identity and inverses were also checked directly",
          "BIN_*_GROUP: all |G|^2 products in G (64, 576, 2304, 14400), unit norm, identity, conj inverses")
    sub_q8_2t = set(Q8s) <= set(T24)
    sub_2t_2o = set(emb(RS2, x) for x in T24) <= set(O48)
    sub_2t_2i = set(emb(RPH, x) for x in T24) <= set(I120)
    check("CHAIN_Q8_IN_2T", sub_q8_2t)
    check("CHAIN_2T_IN_2O", sub_2t_2o)
    check("CHAIN_2T_IN_2I", sub_2t_2i)
    check("LAGRANGE_2O_NOT_IN_2I", 120 % 48 != 0, remainder_120_mod_48=120 % 48,
          meaning="no_subgroup_of_order_48_in_2I_so_no_copy_or_conjugate_of_2O")

    def bq_s2(x):
        return tuple((c[0], c[1], 0, 0) for c in x)

    def bq_ph(x):
        return tuple((c[0], 0, c[1], 0) for c in x)

    inter = set(bq_s2(x) for x in O48) & set(bq_ph(x) for x in I120)
    t24_bq = set(bq_s2(emb(RS2, x)) for x in T24)
    check("INTERSECTION_2O_2I_IS_2T", inter == t24_bq, common=len(inter), basis="1,sqrt2,phi,sqrt2*phi",
          relies_on="NAMED_Q(sqrt2,sqrt5)_has_degree_4")

    # mirror icosian set: odd permutations, and the j-k swap anti-automorphism
    sig = lambda q: (q[0], q[1], q[3], q[2])
    anti_fail = 0
    for a in range(4):
        for b in range(4):
            if sig(qi_mul(T, E4[a], E4[b])) != qi_mul(T, sig(E4[b]), sig(E4[a])):
                anti_fail += 1
    check("SWAP_JK_IS_ANTIAUTOMORPHISM", anti_fail == 0, basis_pairs=16, failures=anti_fail,
          identity="s(xy):s(y)s(x)_for_s_swapping_j,k")
    I120_mirror = [emb(RPH, x) for x in T24] + ico_odd
    mclo = closure_failures(RPH, T, I120_mirror)
    check("ICOSIAN_MIRROR_ODD_PERMS_ALSO_GROUP", mclo == 0, closure_outside=mclo,
          meaning="2T_plus_odd_perms_is_s(2I);the_even_perm_rule_goes_with_the_orientation_ij:k")
    check("ICOSIAN_EVEN_ODD_DISJOINT", len(set(ico_even) & set(ico_odd)) == 0 and len(set(ico_even)) == 96 and len(set(ico_odd)) == 96,
          even=len(set(ico_even)), odd=len(set(ico_odd)))

    # ================================================================ 5. images in SO(3)
    ROT = {}
    for name, R, G, order in groups:
        info = INFO[name]
        if info["closure_fail"] != 0:
            check("ROT_%s" % name, False, reason="group_closure_failed")
            continue
        ra = rot_analysis(R, T, G, info)
        ROT[name] = ra
        ker = set(ra["kernel"])
        expect_img = {"Q8": 4, "2T": 12, "2O": 24, "2I": 60}[name]
        check("ROT_%s_IN_SO3" % name, ra["imag"] and ra["orth"] and ra["det"], ring=R["name"], matrices=len(G),
              stays_imaginary=ra["imag"], MtM_equals_16I=ra["orth"], det_equals_64=ra["det"], scale="M':4M")
        check("ROT_%s_HOMOMORPHISM" % name, ra["hom_fail"] == 0, pairs=len(G) * len(G), failures=ra["hom_fail"],
              identity="M'(x)M'(y):4M'(xy_scaled)")
        check("ROT_%s_PRESERVES_CROSS" % name, ra["cross_fail"] == 0, tests=len(ra["image"]) * 9, failures=ra["cross_fail"])
        check("ROT_%s_KERNEL_PM1" % name, ker == {q_from_ints(R, [2, 0, 0, 0]), q_from_ints(R, [-2, 0, 0, 0])},
              kernel_size=len(ker))
        check("ROT_%s_IMAGE_SIZE" % name, len(ra["image"]) == expect_img, image=len(ra["image"]), expect=expect_img)

    # identification of the images
    if "Q8" in ROT:
        imgs = list(ROT["Q8"]["image"])
        ab = all(m3_mul(RZ, a, b) == m3_mul(RZ, b, a) for a in imgs for b in imgs)
        sixteenI = tuple(tuple(r_int(RZ, 16) if r == c else r_zero(RZ) for c in range(3)) for r in range(3))
        ex2 = all(m3_mul(RZ, a, a) == sixteenI for a in imgs)
        check("ROT_Q8_IMAGE_IS_KLEIN_FOUR", ab and ex2 and len(imgs) == 4, abelian=ab, exponent_2=ex2, size=len(imgs),
              reason="order_4_exponent_2_group_is_Z2xZ2")
    for name, R in (("2T", RZ), ("2O", RS2)):
        if name not in ROT:
            continue
        imgs = list(ROT[name]["image"])
        ims = [quarter_int_matrix(R, M) for M in imgs]
        perms = [diag_perm(Mi) if Mi is not None else None for Mi in ims]
        okp = all(p is not None for p in perms)
        hom = 0
        if okp:
            pm = {Mi: p for Mi, p in zip(ims, perms)}
            for A in ims:
                for B in ims:
                    AB = int_m3_mul(A, B)
                    if AB not in pm or pm[AB] != compose(pm[A], pm[B]):
                        hom += 1
        nperm = len(set(perms)) if okp else 0
        even = all(perm_parity(p) == 0 for p in perms) if okp else False
        if name == "2T":
            check("ROT_2T_IMAGE_IS_A4", okp and hom == 0 and nperm == 12 and even, integer_entries=okp,
                  action="4_cube_diagonals", hom_failures=hom, distinct_perms=nperm, all_even=even,
                  reason="injective_hom_onto_12_even_perms_of_4_is_A4")
        else:
            check("ROT_2O_IMAGE_IS_S4", okp and hom == 0 and nperm == 24, integer_entries=okp, action="4_cube_diagonals",
                  hom_failures=hom, distinct_perms=nperm, all_even=even, reason="injective_hom_onto_all_24_perms_is_S4")
    if INFO["2I"]["closure_fail"] == 0:
        G = I120
        H = set(emb(RPH, x) for x in T24)
        label = {}
        cosets = []
        for g in G:
            if g in label:
                continue
            C = frozenset(q_sp(RPH, T, g, h) for h in H)
            cid = len(cosets)
            cosets.append(C)
            for c in C:
                label[c] = cid
        reps = [next(iter(C)) for C in cosets]
        perm_of = {}
        for g in G:
            perm_of[g] = tuple(label[q_sp(RPH, T, g, r)] for r in reps)
        idxI = INFO["2I"]["idx"]
        tabI = INFO["2I"]["tab"]
        homc = sum(1 for x in G for y in G if perm_of[G[tabI[idxI[x]][idxI[y]]]] != compose(perm_of[x], perm_of[y]))
        idp = tuple(range(len(cosets)))
        kerc = set(g for g in G if perm_of[g] == idp)
        imgc = set(perm_of.values())
        evenc = all(perm_parity(p) == 0 for p in imgc)
        okc = (len(cosets) == 5 and all(len(C) == 24 for C in cosets) and homc == 0
               and kerc == {q_from_ints(RPH, [2, 0, 0, 0]), q_from_ints(RPH, [-2, 0, 0, 0])} and len(imgc) == 60 and evenc)
        check("ROT_2I_IMAGE_IS_A5", okc, action="left_multiplication_on_5_cosets_of_2T", cosets=len(cosets),
              hom_failures=homc, kernel=len(kerc), image_perms=len(imgc), all_even=evenc,
              reason="2I_mod_pm1_embeds_onto_60_even_perms_of_5_:A5;so_rotation_image_:2I_mod_pm1_:A5")
    proof("PROOF_ROTATION_IMAGES", "MEASURED",
          "q to R_q(v):q v conj(q) maps S^3 into SO(3) with kernel {+1,-1}; images of Q8,2T,2O,2I are V4,A4,S4,A5",
          "Re(q v conj q):Re(v conj(q) q):Re(v)N(q) so imaginary v stays imaginary; N(q v conj q):N(v) so R_q is orthogonal; "
          "R_q:id iff q commutes with i,j,k iff q is real (centre of H, see NT_CHARACTERISATION_POLY), so the kernel in S^3 is {+1,-1}; "
          "det R_q:+1 was computed exactly for every element (and holds on all of S^3 by continuity on the connected S^3). "
          "Image identification: Q8 image is abelian of exponent 2 and order 4; 2T and 2O images act faithfully on the 4 cube diagonals "
          "giving 12 even perms (A4) and all 24 perms (S4); 2I acts on its 5 cosets of 2T with kernel {+1,-1} and image the 60 even perms (A5), "
          "and the rotation image is 2I mod {+1,-1} by the first isomorphism theorem",
          "ROT_*_IN_SO3, ROT_*_HOMOMORPHISM, ROT_*_KERNEL_PM1, ROT_*_IMAGE_SIZE, ROT_Q8_IMAGE_IS_KLEIN_FOUR, ROT_2T_IMAGE_IS_A4, "
          "ROT_2O_IMAGE_IS_S4, ROT_2I_IMAGE_IS_A5")
    named("NAMED_AUT_H_IS_SO3",
          "every R-algebra automorphism of H is conjugation by a unit quaternion and Aut(H) is SO(3)",
          "Skolem-Noether theorem (Th._Skolem_1927;_E._Noether_1933_Math._Z._37) for the central simple R-algebra H",
          reason="(1) an automorphism f fixes R and maps Im(H) to itself (Im(H) is 0 plus the x not in R with x^2 real), and on Im(H) x^2:-N(x), "
                 "so f is orthogonal on Im(H); f(i)f(j):f(k) and uv:-(u.v)+u_cross_v give f(k):f(i)_cross_f(j), so det f:+1. "
                 "(2) conversely every rotation preserves dot and cross product hence, by H_IMAG_PRODUCT_POLY, the product; "
                 "and q to R_q is onto SO(3) (Skolem-Noether, or: S^3 mod {+1,-1} is a compact 3-manifold injected into connected SO(3), open and closed image)",
          measured_support="H_IMAG_PRODUCT_POLY,ROT_*_PRESERVES_CROSS,ROT_*_KERNEL_PM1")

    # ================================================================ 6. binary dihedral groups and N(T)
    # N(T) characterisation: q i - i q : 2d j - 2c k ; q i + i q : 2a i - 2b
    iq = [p_const(0, 4), p_const(1, 4), p_const(0, 4), p_const(0, 4)]
    qi_ = qp_mul(T, Qv, iq)
    iq_ = qp_mul(T, iq, Qv)
    dif = [p_sub(qi_[t], iq_[t]) for t in range(4)]
    sm = [p_add(qi_[t], iq_[t]) for t in range(4)]
    ok_nt = (dif == [{}, {}, p_scale(2, Qv[3]), p_scale(-2, Qv[2])] and sm == [p_scale(-2, Qv[1]), p_scale(2, Qv[0]), {}, {}])
    check("NT_CHARACTERISATION_POLY", ok_nt, identities="qi-iq:2d*j-2c*k;qi+iq:2a*i-2b",
          meaning="qi:iq_iff_q_in_R+Ri;qi:-iq_iff_q_in_Rj+Rk;so_N(T):T_union_Tj")
    jq = [p_const(0, 4), p_const(0, 4), p_const(1, 4), p_const(0, 4)]
    zq = [Qv[0], Qv[1], {}, {}]
    zbar = [Qv[0], p_neg(Qv[1]), {}, {}]
    check("J_Z_EQUALS_ZBAR_J_POLY", qp_mul(T, jq, zq) == qp_mul(T, zbar, jq), identity="j(a+bi):(a-bi)j")
    # Cayley-Dickson product matches the table (Gaussian case m:4)
    R4 = r_make_cyclo(4)

    def to_cd(q):
        return ((q[0], q[1]), (q[2], q[3]))

    cdf = 0
    for a in range(4):
        for b in range(4):
            if cd_mul(R4, to_cd(E4[a]), to_cd(E4[b])) != to_cd(qi_mul(T, E4[a], E4[b])):
                cdf += 1
    check("CD_PRODUCT_MATCHES_TABLE", cdf == 0, basis_pairs=16, failures=cdf, ring=R4["name"],
          formula="(z1+w1j)(z2+w2j):(z1z2-w1conj(w2))+(z1w2+w1conj(z2))j")
    dic_rows = []
    for n in range(2, 9):
        m = 2 * n
        R = r_make_cyclo(m)
        D = R["dim"]
        tot = sum(1 for k in range(1, m + 1) if gcd(k, m) == 1)
        zeta = tuple(1 if t == 1 else 0 for t in range(D))
        one = r_int(R, 1)
        zero = r_zero(R)
        cur = zeta
        zord = 1
        while cur != one and zord <= 4 * m:
            cur = r_mul(R, cur, zeta)
            zord += 1
        zn = zeta
        for _ in range(n - 1):
            zn = r_mul(R, zn, zeta)
        zc = r_conj(R, zeta)
        G = cd_generate(R, [(zeta, zero), (zero, one)])
        Zq = (zeta, zero)
        Zi = (zc, zero)
        unit = all(r_add(r_mul(R, g[0], r_conj(R, g[0])), r_mul(R, g[1], r_conj(R, g[1]))) == one for g in G)
        commute = 0
        invert = 0
        neither = 0
        sq_m1 = 0
        for g in G:
            gz = cd_mul(R, g, Zq)
            if gz == cd_mul(R, Zq, g):
                commute += 1
            elif gz == cd_mul(R, Zi, g):
                invert += 1
            else:
                neither += 1
            if not r_iszero(g[1]) and cd_mul(R, g, g) == (r_neg(one), zero):
                sq_m1 += 1
        ok = (D == tot and zord == m and zn == r_neg(one) and r_mul(R, zeta, zc) == one and zc != zeta
              and len(G) == 4 * n and unit and neither == 0 and commute == 2 * n and invert == 2 * n and sq_m1 == 2 * n)
        check("DIC_%d" % n, ok, ring=R["name"], Phi_degree=D, totient=tot, zeta_order=zord, zeta_n_is_minus1=(zn == r_neg(one)),
              order=len(G), expect=4 * n, unit_norm=unit, commute_with_zeta=commute, invert_zeta=invert, neither=neither,
              non_torus_elements_squaring_to_minus1=sq_m1)
        if n == 2:
            q8cd = set(to_cd(x) for x in Q8)
            check("DIC_2_EQUALS_Q8", G == q8cd, common=len(G & q8cd))
        dic_rows.append("%d:%d" % (n, len(G)))
    proof("PROOF_DIC_N_IN_NT_ALL_N", "MEASURED",
          "for every n>:2 the binary dihedral group Dic_n:<zeta,j>, zeta:exp(pi i n^-1), has order 4n and lies in N(T), T the circle cos t+i sin t",
          "j z:conj(z) j (J_Z_EQUALS_ZBAR_J_POLY) and j^2:-1:zeta^n give Dic_n:{zeta^k, zeta^k j : 0<:k<2n}; zeta^k commutes with zeta and "
          "(zeta^k j)zeta:zeta^-1(zeta^k j); since zeta is not real, commuting with zeta means commuting with i and inverting zeta means qi:-iq, "
          "so by NT_CHARACTERISATION_POLY every element lies in (R+Ri) or (Rj+Rk), i.e. in T or Tj, which is N(T)",
          "DIC_2..DIC_8 built from generators in Z[zeta_2n] (exact, mod Phi_2n), element-by-element commute or invert test, orders 4n; "
          "the all-n statement uses the two polynomial identities")
    # control: the same N(T) test on 2T with the torus through i must fail
    ti = q_from_ints(RZ, [0, 2, 0, 0])
    tii = q_from_ints(RZ, [0, -2, 0, 0])
    neither2t = 0
    for g in T24:
        gt = q_sp(RZ, T, g, ti)
        if gt != q_sp(RZ, T, ti, g) and gt != q_sp(RZ, T, tii, g):
            neither2t += 1
    control("CONTROL_NT_TEST_FAILS_ON_2T", neither2t == 0, False, torus_generator="i", neither=neither2t)
    proof("PROOF_EXCEPTIONAL_NOT_IN_ANY_NT", "MEASURED",
          "none of 2T, 2O, 2I lies in any conjugate of N(T)",
          "a conjugate N(T') is T' union T'j' with T' abelian of index 2, so a finite G inside it has the abelian subgroup G meet T' of index <:2; "
          "an abelian subgroup A that contains a non-central g lies in the centraliser C(g); measured max non-central centralisers are 6, 8, 10, "
          "while index <:2 needs |A|>:12, 24, 60; and A inside the centre {+1,-1} would force |G|<:4",
          "BIN_2T,2O,2I_MAX_NONCENTRAL_CENTRALIZER (6,8,10), BIN_*_CENTER_PM1")
    named("NAMED_KLEIN_CLASSIFICATION",
          "every finite subgroup of SU(2) is conjugate to a cyclic group, a binary dihedral group Dic_n (order 4n), 2T (24), 2O (48) or 2I (120)",
          "F._Klein_1884_Vorlesungen_ueber_das_Ikosaeder_und_die_Aufloesung_der_Gleichungen_vom_fuenften_Grade_(Teubner,Leipzig);"
          "J.H._Conway_and_D.A._Smith_2003_On_Quaternions_and_Octonions_(A_K_Peters)_ch.3;"
          "P._Du_Val_1964_Homographies,_Quaternions_and_Rotations_(Oxford)")
    proof("PROOF_SU2_CEILING", "NAMED",
          "every finite subgroup of SU(2) of order >120 lies in a conjugate of N(T); 2I (order 120) is the largest finite subgroup outside every "
          "conjugate of N(T), and 2I and 2O are maximal finite subgroups",
          "by NAMED_KLEIN_CLASSIFICATION a finite G with |G|>120 is cyclic or Dic_n; a cyclic group is generated by one unit quaternion, "
          "which is conjugate into T (conjugation acts on Im(H) as all of SO(3), so it rotates the imaginary part onto the i axis); "
          "Dic_n is conjugate to the standard form, which lies in N(T) (PROOF_DIC_N_IN_NT_ALL_N). The exceptional groups lie in no conjugate "
          "of N(T) (PROOF_EXCEPTIONAL_NOT_IN_ANY_NT), the largest being 2I. Maximality: a finite G properly containing 2I has order a "
          "multiple of 120 beyond 120, so G is in some N(T'), and then so is 2I, a contradiction; for 2O, 48 does not divide 120 "
          "(LAGRANGE_2O_NOT_IN_2I) and any larger order puts G, hence 2O, in some N(T')",
          "DIC_2..DIC_8, J_Z_EQUALS_ZBAR_J_POLY, NT_CHARACTERISATION_POLY, BIN_*_MAX_NONCENTRAL_CENTRALIZER, LAGRANGE_2O_NOT_IN_2I, ROT_*",
          relies_on="NAMED_KLEIN_CLASSIFICATION")
    proof("PROOF_SU2_TOWER_CANNOT_REACH_SU2", "NAMED",
          "the closure of an increasing union of finite subgroups G_1<G_2<... of SU(2) with unbounded order is not SU(2); "
          "its identity component is abelian (trivial or a circle)",
          "for k large |G_k|>120, so G_k lies in some N(T_k) and x^2 lies in the abelian T_k for every x in G_k (x in T_k j_k gives x^2:-1); "
          "any two x,y of the union lie in one G_k, so x^2 and y^2 commute; by continuity x^2 y^2:y^2 x^2 on the closure; on SU(2) squaring is onto "
          "(cos t+u sin t is the square of cos(t:2)+u sin(t:2)), so if the closure were SU(2) all of SU(2) would commute, but ij is not ji (Q8_NON_ABELIAN); "
          "the same squaring argument on the identity component of the closure, where every element lies in a maximal torus, gives commutativity there",
          "PROOF_SU2_CEILING inputs, Q8_NON_ABELIAN",
          relies_on="NAMED_KLEIN_CLASSIFICATION,NAMED_MAXIMAL_TORUS")
    named("NAMED_JORDAN_THEOREM",
          "for each n there is J(n) such that every finite subgroup of GL(n,C) has an abelian normal subgroup of index at most J(n)",
          "C._Jordan_1878_Memoire_sur_les_equations_differentielles_lineaires_a_integrale_algebrique_J._reine_angew._Math._84_89-215;"
          "optimal_bounds:_M.J._Collins_2007_On_Jordan's_theorem_for_complex_linear_groups_J._Group_Theory_10_411-423")
    named("NAMED_MAXIMAL_TORUS",
          "in a compact connected Lie group every element lies in a maximal torus, and all maximal tori are conjugate",
          "E._Cartan;_e.g._T._Broecker_and_T._tom_Dieck_1985_Representations_of_Compact_Lie_Groups_(GTM_98)_ch.IV")
    proof("PROOF_JORDAN_TOWER_GENERAL_N", "NAMED",
          "if a compact subgroup H of U(n) is the closure of an increasing union of finite subgroups, the identity component of H is a torus, so "
          "H is a torus extended by a finite group and H is not SU(N) for N>:2",
          "let M:J(n)!; in each finite G_k with abelian normal A_k of index <:J(n), x^M lies in A_k for every x, so x^M and y^M commute for x,y in G_k; "
          "any two elements of the union lie in one G_k, so x^M y^M:y^M x^M on the union and by continuity on H; in the identity component H_0 every "
          "element lies in a maximal torus where the M-th power map is onto, so H_0 is abelian, a connected compact abelian Lie group, a torus; "
          "H mod H_0 is finite by compactness; SU(N) is connected and non-abelian",
          "SU(2) case measured: PROOF_SU2_TOWER_CANNOT_REACH_SU2",
          relies_on="NAMED_JORDAN_THEOREM,NAMED_MAXIMAL_TORUS")
    emit("LADDER", orders="Q8:8,2T:24,2O:48,2I:120", chain="Q8<2T<2O_and_2T<2I", branch="2O_not_in_2I",
         dic_orders=",".join(dic_rows), claim="MEASURED")

    # ================================================================ 7. broken controls
    Tb = build_qtable(flip=(1, 2))
    afb = assoc_failures(Tb)
    control("CONTROL_FLIPPED_SIGN_ASSOC_FAILS", afb == 0, False, change="ij:+k_to_ij:-k", failing_basis_triples=afb)
    I119 = I120[:24] + I120[25:]
    c119 = closure_failures(RPH, T, I119)
    control("CONTROL_2I_MINUS_ONE_CLOSURE_FAILS", c119 == 0, False, dropped=I120[24], size=len(I119), products_outside=c119)
    mixed = [emb(RPH, x) for x in T24] + ico_even + ico_odd
    cmix = closure_failures(RPH, T, mixed)
    control("CONTROL_ALL_PERMS_CLOSURE_FAILS", cmix == 0, False, size=len(mixed), products_outside=cmix)
    hb = 0
    for a in range(4):
        for b in range(4):
            if phi_int_bad(qi_mul(T, E4[a], E4[b])) != gm_mul(phi_int_bad(E4[a]), phi_int_bad(E4[b])):
                hb += 1
    control("CONTROL_BROKEN_PHI_HOM_FAILS", hb == 0, False, change="lower_left_-c_to_+c", failing_pairs=hb)
    _, _, diffb, _ = norm_mult_poly(T, (1, 1, 1, -1))
    control("CONTROL_WRONG_NORM_FAILS", diffb == {}, False, form="a^2+b^2+c^2-d^2", difference_terms=len(diffb))
    bad_rot = 0
    for X in T24:
        M, im = rot_matrix(RZ, T, X, use_conj=False)
        if not (im and orth_ok(RZ, M, 16) and det3(RZ, M) == r_int(RZ, 64)):
            bad_rot += 1
    control("CONTROL_NO_CONJUGATE_NOT_ROTATION", bad_rot == 0, False, map="v_to_q_v_q", test="imaginary+orthogonal+det_64",
            failing_elements=bad_rot, of=len(T24))
    O47 = O48[:47]
    c47 = closure_failures(RS2, T, O47)
    control("CONTROL_2O_MINUS_ONE_CLOSURE_FAILS", c47 == 0, False, size=len(O47), products_outside=c47)

    emit("SUMMARY", checks=STATE["checks"], failures=STATE["fail"], controls=STATE["controls"],
         failed=(",".join(STATE["failed_tags"]) if STATE["failed_tags"] else "none"),
         result="PASS" if STATE["fail"] == 0 else "FAIL", claim="MEASURED")


def write_out():
    data = ("\n".join(ROWS) + "\n").encode()
    with open(OUT, "wb") as fh:
        fh.write(data)
    sha = hashlib.sha256(data).hexdigest()
    with open(SHA_OUT, "wb") as fh:
        fh.write(("%s  %s\n" % (sha, OUT_NAME)).encode())
    return sha


if __name__ == "__main__":
    rc = 0
    try:
        main()
    except Exception as ex:
        STATE["fail"] += 1
        emit("ERROR", result="FAIL", kind=type(ex).__name__, detail=str(ex)[:300], claim="MEASURED")
        rc = 3
    sha = write_out()
    for r in ROWS:
        print(r)
    print("HBP_SHA256=%s" % sha)
    if STATE["fail"] != 0 and rc == 0:
        rc = 1
    sys.exit(rc)
