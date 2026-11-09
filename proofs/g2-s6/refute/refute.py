#!/usr/bin/env python3
# REFUTE seat script: independent second-seat attack on M1/M2/M3 and C1..C10.
# Seat: Opus 5.5 second seat for ACER-CLAUDE-FABLE5 pid 8467a937cba309f7. E=0:
# compute only, no network, writes only REFUTE.hbp and REFUTE.hbp.sha256 here.
# Integer only: no true or floor division operators, no float literals, no
# float/fractions/decimal/complex names, math only gcd. Exact division: ediv().
import sys
import os
import io
import hashlib
import tokenize
import itertools
from math import gcd

SELF = os.path.abspath(__file__)
HERE = os.path.dirname(SELF)
BASE = os.path.dirname(HERE)
OUT_NAME = "REFUTE.hbp"
OUT = os.path.join(HERE, OUT_NAME)
SIDE = os.path.join(HERE, OUT_NAME + ".sha256")

ROWS = []
FAILS = []
NPASS = [0]
NCTL = [0]


def fmt(v):
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (list, tuple)):
        return ",".join(fmt(x) for x in v)
    s = str(v)
    return s.replace("|", ";").replace("=", ":").replace("\n", " ")


def row(tag, **kv):
    parts = [tag] + ["%s=%s" % (k, fmt(v)) for k, v in kv.items()] + ["json=0"]
    ROWS.append("|".join(parts))


def check(name, ok, claim="MEASURED", **kv):
    row("CHECK", name=name, status=("PASS" if ok else "FAIL"), claim=claim, **kv)
    if ok:
        NPASS[0] += 1
    else:
        FAILS.append(name)
    return ok


def control(name, observed_ok, **kv):
    # A deliberately broken input. The underlying check must come out FAIL.
    NCTL[0] += 1
    row("CONTROL", name=name, expected="FAIL", observed=("PASS" if observed_ok else "FAIL"),
        status=("FAIL" if observed_ok else "PASS"), claim="MEASURED", **kv)
    if observed_ok:
        FAILS.append("control_" + name)
    else:
        NPASS[0] += 1


# ------------------------------------------------------------------ lint
ALLOWED_IMPORTS = {"sys", "os", "io", "hashlib", "tokenize", "itertools", "math"}
BANNED = {"float", "complex", "Fraction", "fractions", "Decimal", "decimal", "cmath",
          "numpy", "truediv", "floordiv"}


def lint(src):
    bad = []
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    for i, t in enumerate(toks):
        if t.type == tokenize.OP and t.string in ("/", "/=", "//", "//="):
            bad.append("divop_line%d" % t.start[0])
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if not s.startswith(("0x", "0o", "0b")) and ("." in s or "e" in s or "j" in s):
                bad.append("floatlit_line%d" % t.start[0])
        elif t.type == tokenize.NAME and t.string in BANNED:
            bad.append("banned_%s_line%d" % (t.string, t.start[0]))
        elif t.type == tokenize.NAME and t.string == "math" and i + 2 < len(toks) \
                and toks[i + 1].string == ".":
            if toks[i + 2].string not in ("gcd", "isqrt"):
                bad.append("math_%s_line%d" % (toks[i + 2].string, t.start[0]))
        elif t.type == tokenize.NAME and t.string in ("import", "from") and t.start[1] == 0:
            mod = toks[i + 1].string
            if mod not in ALLOWED_IMPORTS:
                bad.append("import_%s_line%d" % (mod, t.start[0]))
            if t.string == "from" and mod == "math":
                j = i + 3
                while j < len(toks) and toks[j].type != tokenize.NEWLINE:
                    if toks[j].type == tokenize.NAME and toks[j].string not in ("gcd", "isqrt"):
                        bad.append("math_import_%s" % toks[j].string)
                    j += 1
    return bad


# ------------------------------------------------------------ integers
def ediv(a, b):
    q, r = divmod(a, b)
    if r:
        raise ArithmeticError("inexact division %d by %d" % (a, b))
    return q


def ilcm(a, b):
    a, b = abs(a), abs(b)
    if a == 0 or b == 0:
        return 0
    return ediv(a * b, gcd(a, b))


def normalize(r):
    g = 0
    for x in r:
        g = gcd(g, x)
    if g == 0:
        return list(r)
    lead = 0
    for x in r:
        if x:
            lead = x
            break
    if lead < 0:
        g = -g
    return [ediv(x, g) for x in r]


def rref(rows, ncols):
    R = [list(r) for r in rows if any(r)]
    piv = []
    r = 0
    for c in range(ncols):
        if r >= len(R):
            break
        p = None
        for i in range(r, len(R)):
            if R[i][c]:
                p = i
                break
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        R[r] = normalize(R[r])
        pr = R[r]
        a = pr[c]
        for i in range(len(R)):
            if i != r and R[i][c]:
                b = R[i][c]
                R[i] = normalize([a * x - b * y for x, y in zip(R[i], pr)])
        piv.append(c)
        r += 1
    return R[:r], piv


def rank(rows, ncols):
    return len(rref(rows, ncols)[1])


def nullspace(rows, ncols):
    R, piv = rref(rows, ncols)
    L = 1
    for k, c in enumerate(piv):
        L = ilcm(L, R[k][c])
    basis = []
    for f in range(ncols):
        if f in piv:
            continue
        v = [0] * ncols
        v[f] = L
        for k, c in enumerate(piv):
            v[c] = ediv(-R[k][f] * L, R[k][c])
        basis.append(normalize(v))
    ok = all(sum(x * y for x, y in zip(rw, v)) == 0 for rw in rows for v in basis)
    return basis, len(piv), ok


def det_bareiss(M):
    n = len(M)
    if n == 0:
        return 1
    A = [list(r) for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if A[i][k]:
                    sw = i
                    break
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = ediv(A[i][j] * A[k][k] - A[i][k] * A[k][j], prev)
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def adjugate(M):
    n = len(M)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor = [[M[r][c] for c in range(n) if c != j] for r in range(n) if r != i]
            adj[j][i] = (-1 if (i + j) & 1 else 1) * det_bareiss(minor)
    return adj


def mmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def madd(A, B, c=1):
    return [[a + c * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, c):
    return [[c * a for a in r] for r in A]


def comm(A, B):
    return madd(mmul(A, B), mmul(B, A), -1)


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def zeros(n):
    return [[0] * n for _ in range(n)]


def flat(A):
    return [x for r in A for x in r]


def perm_sign(seq):
    inv = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                inv += 1
    return -1 if inv & 1 else 1


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vneg(a):
    return [-x for x in a]


def vscale(a, c):
    return [c * x for x in a]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def unit(i, n):
    v = [0] * n
    v[i] = 1
    return v


# ------------------------------------------------- composition algebras
BAEZ = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]
BRYANT = [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 7, 5), (3, 7, 4), (3, 6, 5)]
QUAT = [(1, 2, 3)]


def make_table(triples, dim):
    T = [[None] * dim for _ in range(dim)]
    conflict = [0]

    def put(i, j, v):
        if T[i][j] is not None and T[i][j] != v:
            conflict[0] += 1
        T[i][j] = v
    for i in range(dim):
        put(0, i, (1, i))
        put(i, 0, (1, i))
    for i in range(1, dim):
        put(i, i, (-1, 0))
    for (a, b, c) in triples:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            put(x, y, (1, z))
            put(y, x, (-1, z))
    complete = all(T[i][j] is not None for i in range(dim) for j in range(dim))
    return T, complete, conflict[0]


def amul(T, x, y):
    n = len(x)
    r = [0] * n
    for i in range(n):
        xi = x[i]
        if xi:
            Ti = T[i]
            for j in range(n):
                yj = y[j]
                if yj:
                    s, k = Ti[j]
                    r[k] += s * xi * yj
    return r


def assoc_profile(T):
    dim = len(T)
    E = [unit(i, dim) for i in range(dim)]
    A = {}
    for a in range(dim):
        for b in range(dim):
            ab = amul(T, E[a], E[b])
            for c in range(dim):
                A[(a, b, c)] = vsub(amul(T, ab, E[c]), amul(T, E[a], amul(T, E[b], E[c])))
    fails = 0
    for (a, b, c), v in A.items():
        if A[(b, a, c)] != vneg(v) or A[(a, c, b)] != vneg(v):
            fails += 1
    nonzero = sum(1 for v in A.values() if any(v))
    return fails, nonzero


def lcg(seed):
    s = [seed]

    def nxt(lo, hi):
        s[0] = (s[0] * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return lo + (s[0] >> 33) % (hi - lo + 1)
    return nxt


def norm_mult_fails(T, pairs, seed):
    nx = lcg(seed)
    dim = len(T)
    bad = 0
    for _ in range(pairs):
        x = [nx(-5, 5) for _ in range(dim)]
        y = [nx(-5, 5) for _ in range(dim)]
        if dot(amul(T, x, y), amul(T, x, y)) != dot(x, x) * dot(y, y):
            bad += 1
    return bad


def cross_tab(T):
    m = len(T) - 1
    CR = [[None] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                s, k = T[i + 1][j + 1]
                CR[i][j] = (s, k - 1)
    return CR


def vcross(CR, u, v):
    m = len(CR)
    r = [0] * m
    for i in range(m):
        if u[i]:
            for j in range(m):
                if v[j] and i != j:
                    s, k = CR[i][j]
                    r[k] += s * u[i] * v[j]
    return r


# --------------------------------------------- polynomial vector fields
def p_add(p, q, c=1):
    r = dict(p)
    for m, a in q.items():
        v = r.get(m, 0) + c * a
        if v:
            r[m] = v
        elif m in r:
            del r[m]
    return r


def p_mul(p, q):
    r = {}
    for m1, a in p.items():
        for m2, b in q.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            v = r.get(m, 0) + a * b
            if v:
                r[m] = v
            elif m in r:
                del r[m]
    return r


def p_var(i, nv):
    e = [0] * nv
    e[i] = 1
    return {tuple(e): 1}


def p_const(c, nv):
    return {(0,) * nv: c} if c else {}


def p_diff(p, i):
    r = {}
    for m, a in p.items():
        if m[i]:
            e = list(m)
            e[i] -= 1
            t = tuple(e)
            r[t] = r.get(t, 0) + a * m[i]
    return {k: v for k, v in r.items() if v}


def p_eval(p, pt):
    tot = 0
    for m, a in p.items():
        t = a
        for x, e in zip(pt, m):
            if e:
                t *= x ** e
        tot += t
    return tot


def f_cross(geo, V, W):
    nv = geo["nv"]
    R = [dict() for _ in range(nv)]
    for off, CR in geo["blocks"]:
        m = len(CR)
        for i in range(m):
            Vi = V[off + i]
            if not Vi:
                continue
            for j in range(m):
                if i == j:
                    continue
                Wj = W[off + j]
                if not Wj:
                    continue
                s, k = CR[i][j]
                R[off + k] = p_add(R[off + k], p_mul(Vi, Wj), s)
    return R


def f_J(geo, V, js):
    nv = geo["nv"]
    pos = [p_var(i, nv) for i in range(nv)]
    return [p_add({}, c, js) for c in f_cross(geo, pos, V)]


def f_tangent(geo, p):
    nv = geo["nv"]
    idx = None
    for off, CR in geo["blocks"]:
        if off <= p < off + len(CR):
            idx = list(range(off, off + len(CR)))
    sq = {}
    for i in idx:
        sq = p_add(sq, p_mul(p_var(i, nv), p_var(i, nv)))
    V = [dict() for _ in range(nv)]
    for k in idx:
        t = p_mul(p_var(p, nv), p_var(k, nv))
        V[k] = p_add(sq if k == p else {}, t, -1)
    return V


def jet(V, pt):
    nv = len(V)
    val = [p_eval(c, pt) for c in V]
    jac = [[p_eval(p_diff(c, i), pt) for i in range(nv)] for c in V]
    return val, jac


def br_jet(jv, jw):
    vv, vj = jv
    wv, wj = jw
    nv = len(vv)
    return [sum(vv[i] * wj[k][i] - wv[i] * vj[k][i] for i in range(nv)) for k in range(nv)]


def pt_cross(geo, a, b):
    r = [0] * geo["nv"]
    for off, CR in geo["blocks"]:
        m = len(CR)
        c = vcross(CR, a[off:off + m], b[off:off + m])
        for k in range(m):
            r[off + k] = c[k]
    return r


def pt_J(geo, pt, w, js):
    return vscale(pt_cross(geo, pt, w), js)


def pt_proj(geo, pt, w):
    r = list(w)
    for off, CR in geo["blocks"]:
        m = len(CR)
        pb = pt[off:off + m]
        wb = w[off:off + m]
        d = dot(pb, wb)
        for k in range(m):
            r[off + k] = wb[k] - d * pb[k]
    return r


def nij(geo, pt, jX, jJX, jY, jJY, js, last_sign=-1):
    a = br_jet(jJX, jJY)
    b = pt_J(geo, pt, br_jet(jJX, jY), js)
    c = pt_J(geo, pt, br_jet(jX, jJY), js)
    d = br_jet(jX, jY)
    return [a[k] - b[k] - c[k] + last_sign * d[k] for k in range(len(a))]


def nabla(geo, pt, jX, jY):
    xv = jX[0]
    yj = jY[1]
    n = len(xv)
    D = [sum(yj[k][i] * xv[i] for i in range(n)) for k in range(n)]
    return pt_proj(geo, pt, D)


def sphere_run(geo, pt, tang, js):
    F = {}
    for p in tang:
        X = f_tangent(geo, p)
        JX = f_J(geo, X, js)
        F[p] = (X, JX, jet(X, pt), jet(JX, pt))
    N = {}
    NJ = {}
    for a in tang:
        for b in tang:
            jX, jJX = F[a][2], F[a][3]
            jY, jJY = F[b][2], F[b][3]
            N[(a, b)] = nij(geo, pt, jX, jJX, jY, jJY, js)
            NJ[(a, b)] = vsub(nabla(geo, pt, jX, jJY), pt_J(geo, pt, nabla(geo, pt, jX, jY), js))
    return F, N, NJ


def sphere_block(label, geo, pt, tang, js, expect_integrable):
    nv = geo["nv"]
    F, N, NJ = sphere_run(geo, pt, tang, js)
    vals_ok = all(F[p][2][0] == unit(p, nv) for p in tang)
    tangent_ok = all(dot(v, pt) == 0 for v in N.values())
    nz = sorted(k for k, v in N.items() if any(v))
    zero_pairs = [k for k, v in N.items() if not any(v)]
    lin = {}
    for a, b in N:
        cab = pt_proj(geo, pt, pt_cross(geo, unit(a, nv), unit(b, nv)))
        lin[(a, b)] = cab
    nk_formula = all(NJ[k] == vscale(lin[k], js) for k in N)
    nk_antisym = all(vadd(NJ[(a, b)], NJ[(b, a)]) == [0] * nv for a, b in N)
    nk_strict = any(any(v) for v in NJ.values())
    n_eq_m4JnablaJ = all(N[k] == vscale(pt_J(geo, pt, NJ[k], js), -4) for k in N)
    n_eq_m4Jpi = all(N[k] == vscale(pt_J(geo, pt, lin[k], 1), -4) for k in N)
    row("GEOM_NIJENHUIS", geometry=label, point=pt, js=js, ordered_pairs=len(N),
        nonzero_pairs=len(nz), zero_pairs=len(zero_pairs), fields_restrict_ok=vals_ok,
        tangent=tangent_ok, nablaJ_eq_js_pi_cross=nk_formula, nablaJ_antisym=nk_antisym,
        nablaJ_nonzero=nk_strict, N_eq_minus4_J_nablaJ=n_eq_m4JnablaJ,
        N_eq_minus4_Jplus_pi_cross=n_eq_m4Jpi, method="ambient_polynomial_vector_fields_no_Lie_algebra",
        claim="MEASURED")
    return F, N, NJ, dict(vals_ok=vals_ok, tangent_ok=tangent_ok, nz=len(nz), zero=len(zero_pairs),
                          nk_formula=nk_formula, nk_antisym=nk_antisym, nk_strict=nk_strict,
                          n4=n_eq_m4JnablaJ, n4pi=n_eq_m4Jpi, zero_list=zero_pairs)


def tensoriality(geo, pt, tang, js, q, last_sign):
    nv = geo["nv"]
    f = p_add(p_const(3, nv), p_var(q, nv))
    fv = p_eval(f, pt)
    bad = 0
    for a in tang:
        X = f_tangent(geo, a)
        fX = [p_mul(f, c) for c in X]
        jX, jJX = jet(X, pt), jet(f_J(geo, X, js), pt)
        jfX, jJfX = jet(fX, pt), jet(f_J(geo, fX, js), pt)
        for b in tang:
            Y = f_tangent(geo, b)
            jY, jJY = jet(Y, pt), jet(f_J(geo, Y, js), pt)
            n0 = nij(geo, pt, jX, jJX, jY, jJY, js, last_sign)
            n1 = nij(geo, pt, jfX, jJfX, jY, jJY, js, last_sign)
            n2 = nij(geo, pt, jY, jJY, jfX, jJfX, js, last_sign)
            if n1 != vscale(n0, fv) or n2 != vscale(nij(geo, pt, jY, jJY, jX, jJX, js, last_sign), fv):
                bad += 1
    return fv, bad


# ------------------------------------------------- exterior algebra R^7
def w_wedge(a, b):
    r = {}
    for I, x in a.items():
        sI = set(I)
        for K, y in b.items():
            if sI.intersection(K):
                continue
            seq = I + K
            key = tuple(sorted(seq))
            v = r.get(key, 0) + perm_sign(seq) * x * y
            if v:
                r[key] = v
            elif key in r:
                del r[key]
    return r


def w_hodge(a, n, o):
    full = set(range(1, n + 1))
    r = {}
    for I, x in a.items():
        Ic = tuple(sorted(full - set(I)))
        r[Ic] = r.get(Ic, 0) + o * perm_sign(I + Ic) * x
    return {k: v for k, v in r.items() if v}


def w_interior(i, a):
    r = {}
    for I, x in a.items():
        if i in I:
            pos = I.index(i)
            K = I[:pos] + I[pos + 1:]
            r[K] = r.get(K, 0) + (-1 if pos & 1 else 1) * x
    return {k: v for k, v in r.items() if v}


def phi_form(T):
    r = {}
    for a, b, c in itertools.combinations(range(1, 8), 3):
        s, k = T[b][c]
        if k == a:
            r[(a, b, c)] = s
    return r


PAIRS = list(itertools.combinations(range(1, 8), 2))


def two_vec(form):
    return [form.get(P, 0) for P in PAIRS]


def g2_split(label, T):
    phi = phi_form(T)
    alt_bad = 0
    for a, b, c in itertools.permutations(range(1, 8), 3):
        s, k = T[b][c]
        val = s if k == a else 0
        key = tuple(sorted((a, b, c)))
        # perm_sign((a,b,c)) is the sign relative to sorted order since values are distinct
        if val != perm_sign((a, b, c)) * phi.get(key, 0):
            alt_bad += 1
    vol = (1, 2, 3, 4, 5, 6, 7)
    B = [[w_wedge(w_wedge(w_interior(i, phi), w_interior(j, phi)), phi).get(vol, 0)
          for j in range(1, 8)] for i in range(1, 8)]
    c0 = B[0][0]
    o = 1 if c0 > 0 else -1
    b_diag = (B == mscale(eye(7), c0)) and abs(c0) == 6
    Tm = [[0] * 21 for _ in range(21)]
    for cidx, P in enumerate(PAIRS):
        img = two_vec(w_hodge(w_wedge(phi, {P: 1}), 7, o))
        for ridx in range(21):
            Tm[ridx][cidx] = img[ridx]
    sym = all(Tm[i][j] == Tm[j][i] for i in range(21) for j in range(21))
    T2 = mmul(Tm, Tm)
    minpoly = (madd(madd(T2, Tm, -1), eye(21), -2) == zeros(21))
    r_m2 = rank(madd(Tm, eye(21), -2), 21)
    r_p1 = rank(madd(Tm, eye(21), 1), 21)
    dim7 = 21 - r_m2
    dim14 = 21 - r_p1
    row("G2_FORM_SPLIT", table=label, phi_terms=len(phi), phi=";".join(
        "%s%d" % ("+" if v > 0 else "-", 100 * I[0] + 10 * I[1] + I[2]) for I, v in sorted(phi.items())),
        alternating_failures=alt_bad, B_diag=c0, B_is_6delta_up_to_orientation=b_diag,
        positive_orientation_sign=o, T_symmetric=sym, minpoly_T2_minus_T_minus_2=minpoly,
        dim_eig_plus2=dim7, dim_eig_minus1=dim14, claim="MEASURED")
    return phi, o, Tm, dict(alt=alt_bad, bdiag=b_diag, sym=sym, minpoly=minpoly, d7=dim7, d14=dim14)


def der_basis(T):
    dim = len(T)
    nun = dim * dim
    rows = []
    E = [unit(i, dim) for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            s, t = T[i][j]
            for m in range(dim):
                r = [0] * nun
                # D(e_i e_j)_m = s * D[m][t]
                r[m * dim + t] += s
                # - (D(e_i) e_j)_m = - sum_k D[k][i] (e_k e_j)_m
                for k in range(dim):
                    s2, u = T[k][j]
                    if u == m:
                        r[k * dim + i] -= s2
                # - (e_i D(e_j))_m
                for k in range(dim):
                    s3, u = T[i][k]
                    if u == m:
                        r[k * dim + j] -= s3
                rows.append(r)
    basis, rk, ok = nullspace(rows, nun)
    mats = [[v[r * dim:(r + 1) * dim] for r in range(dim)] for v in basis]
    leib_bad = 0
    for D in mats:
        for i in range(dim):
            for j in range(dim):
                lhs = [sum(D[m][k] * x for k, x in enumerate(amul(T, E[i], E[j]))) for m in range(dim)]
                Di = [D[m][i] for m in range(dim)]
                Dj = [D[m][j] for m in range(dim)]
                rhs = vadd(amul(T, Di, E[j]), amul(T, E[i], Dj))
                if lhs != rhs:
                    leib_bad += 1
    return mats, rk, ok, len(rows), leib_bad


def der_as_2form(D):
    return [D[i][j] for (i, j) in PAIRS]


def chirality(label, T, phi, o, Tm):
    psi = w_hodge(phi, 7, o)
    out = []
    for line in (BAEZ if label == "baez" else BRYANT):
        A = tuple(sorted(line))
        C = tuple(sorted(set(range(1, 8)) - set(A)))
        phiA = phi.get(A, 0)
        eps = psi.get(C, 0)
        orient_agree = (phiA * eps * perm_sign(A + C)) == o
        p0 = C[0]
        sd, asd = [], []
        for x in C[1:]:
            rest = tuple(y for y in C if y not in (p0, x))
            s = eps * perm_sign((p0, x) + rest)
            sd.append({(p0, x): 1, rest: s})
            asd.append({(p0, x): 1, rest: -s})

        def in14(f):
            v = two_vec(f)
            tv = [sum(Tm[i][j] * v[j] for j in range(21)) for i in range(21)]
            return tv == vneg(v)

        def wpsi0(f):
            return not w_wedge(f, psi)
        sd14 = sum(1 for f in sd if in14(f))
        asd14 = sum(1 for f in asd if in14(f))
        asdpsi = all(wpsi0(f) for f in asd)
        sdpsi = all(wpsi0(f) for f in sd)
        out.append((A, phiA, eps, orient_agree, sd14, asd14, asdpsi, sdpsi))
    return out


# --------------------------------------------------------- monomial groups
def mono_mul(A, B, M):
    pa, ea = A
    pb, eb = B
    p = tuple(pb[pa[i]] for i in range(3))
    e = tuple((ea[i] + eb[pa[i]]) % M for i in range(3))
    return (p, e)


def mono_inv(A, M):
    pa, ea = A
    q = [0] * 3
    e = [0] * 3
    for i in range(3):
        q[pa[i]] = i
        e[pa[i]] = (-ea[i]) % M
    return (tuple(q), tuple(e))


def closure(gens, mul, ident):
    G = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = mul(g, h)
                if x not in G:
                    G.add(x)
                    nxt.append(x)
        frontier = nxt
    return G


def classes_of(G, mul, inv):
    left = set(G)
    cls = []
    Gl = list(G)
    while left:
        x = next(iter(left))
        c = set(mul(mul(g, x), inv(g)) for g in Gl)
        cls.append(c)
        left -= c
    return cls


def abelian_normal_subgroups(G, mul, inv, ident):
    cls = classes_of(G, mul, inv)
    idc = [i for i, c in enumerate(cls) if ident in c][0]
    n = len(cls)

    def commute_sets(A, B):
        return all(mul(a, b) == mul(b, a) for a in A for b in B)
    internal = [commute_sets(c, c) for c in cls]
    comp = [[commute_sets(cls[i], cls[j]) for j in range(n)] for i in range(n)]
    others = [i for i in range(n) if i != idc and internal[i]]
    found = []

    def rec(start, chosen):
        S = set().union(*[cls[i] for i in chosen])
        if all(mul(a, b) in S for a in S for b in S):
            found.append(frozenset(S))
        for k in range(start, len(others)):
            i = others[k]
            if all(comp[i][j] for j in chosen):
                rec(k + 1, chosen + [i])
    rec(0, [idc])
    return cls, found


def flat_sector_count(G, mul):
    Gl = list(G)
    idx = {g: i for i, g in enumerate(Gl)}
    n = len(Gl)
    cent = []
    for g in Gl:
        b = 0
        for h in Gl:
            if mul(g, h) == mul(h, g):
                b |= 1 << idx[h]
        cent.append(b)
    tot = 0
    for i in range(n):
        ci = cent[i]
        for j in range(n):
            if (ci >> j) & 1:
                tot += bin(ci & cent[j]).count("1")
    return ediv(tot, n)


# ----------------------------------------------------- Z[phi] quaternions
def zf_mul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def zf_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zf_neg(x):
    return (-x[0], -x[1])


def q_mul(X, Y):
    w1, x1, y1, z1 = X
    w2, x2, y2, z2 = Y
    m = zf_mul
    a = zf_add
    n = zf_neg
    w = a(a(m(w1, w2), n(m(x1, x2))), a(n(m(y1, y2)), n(m(z1, z2))))
    x = a(a(m(w1, x2), m(x1, w2)), a(m(y1, z2), n(m(z1, y2))))
    y = a(a(m(w1, y2), n(m(x1, z2))), a(m(y1, w2), m(z1, x2)))
    z = a(a(m(w1, z2), m(x1, y2)), a(n(m(y1, x2)), m(z1, w2)))
    return (w, x, y, z)


def half(c):
    return (ediv(c[0], 2), ediv(c[1], 2))


def q2mul(X, Y):
    return tuple(half(c) for c in q_mul(X, Y))


def q_norm4(X):
    t = (0, 0)
    for c in X:
        t = zf_add(t, zf_mul(c, c))
    return t


def q_conj(X):
    return (X[0], zf_neg(X[1]), zf_neg(X[2]), zf_neg(X[3]))


def icosians(base):
    Z0 = (0, 0)
    S = set()
    for i in range(4):
        for s in (2, -2):
            v = [Z0] * 4
            v[i] = (s, 0)
            S.add(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        S.add(tuple((s, 0) for s in signs))
    for perm in itertools.permutations(range(4)):
        if perm_sign(perm) != 1:
            continue
        for signs in itertools.product((1, -1), repeat=3):
            vals = [base[0]] + [(s * base[k + 1][0], s * base[k + 1][1]) for k, s in enumerate(signs)]
            v = [None] * 4
            for pos in range(4):
                v[pos] = vals[perm[pos]]
            S.add(tuple(v))
    return S


ONE_Q = ((2, 0), (0, 0), (0, 0), (0, 0))
MONE_Q = ((-2, 0), (0, 0), (0, 0), (0, 0))


def q_order(X):
    P = X
    k = 1
    while P != ONE_Q:
        P = q2mul(P, X)
        k += 1
        if k > 200:
            return -1
    return k


# =================================================================== main
def main_body():
    src = open(SELF, encoding="utf-8").read()
    src_sha = hashlib.sha256(src.encode("utf-8")).hexdigest()
    row("HDR", file=OUT_NAME, script="refute.py", script_sha256=src_sha,
        python=sys.version.split()[0], seat="Opus 5.5 second seat for ACER-CLAUDE-FABLE5 pid 8467a937cba309f7",
        role="REFUTER", E=0, network=0, float_used=0, integer_only=1, date="2026-11-09")
    bad = lint(src)
    check("self_lint_integer_only", not bad, findings=len(bad), detail=(";".join(bad[:8]) or "none"))
    sample = "x = 1 " + chr(47) + " 3\ny = 2" + chr(46) + "5\nfrom fractions import Fraction\nimport math\nz = math.sqrt(2)\n"
    control("lint_catches_bad_source", not lint(sample), flagged=len(lint(sample)))

    # ---------------- reproduction of M1/M2/M3 hbp bytes
    stated = {"m1/M1-OCTONION.hbp": "0024b12eba7a752688ed0847e6a6f05d180c1fe333e95690008cfcd3e6e23fc9",
              "m2/M2-ROOTS.hbp": "0503594c0e3696eb20ff2c84b2bcc7a3f8c5d1d3d9cfe83560e3045694181b98",
              "m3/M3-TOWERS.hbp": "86b840b17f835b54cbc36051e40818c8ca9f9c261bead3d06f0a732d021acd5b"}
    for rel, want in sorted(stated.items()):
        orig = os.path.join(BASE, rel)
        rerun = os.path.join(HERE, "rerun", rel)
        h_o = hashlib.sha256(open(orig, "rb").read()).hexdigest() if os.path.exists(orig) else "missing"
        h_r = hashlib.sha256(open(rerun, "rb").read()).hexdigest() if os.path.exists(rerun) else "missing"
        check("reproduce_" + rel.split("/")[0], h_o == want and h_r == want, stated=want,
              original=h_o, my_rerun=h_r, rerun_cmd="python3 -I -B on a byte-identical copy in refute/rerun")
    aud = os.path.join(HERE, "audit_rows.txt")
    if os.path.exists(aud):
        audsrc = os.path.join(HERE, "audit_float.py")
        ah = hashlib.sha256(open(audsrc, "rb").read()).hexdigest()
        lines = [ln.strip() for ln in open(aud, encoding="utf-8") if ln.strip()]
        allpass = all("status=PASS" in ln for ln in lines)
        for ln in lines:
            ROWS.append(ln.replace("|json=0", "|source=audit_float.py|source_sha256=%s|json=0" % ah))
        check("audit_rows_all_pass", allpass and len(lines) == 8, rows=len(lines))
    else:
        check("audit_rows_present", False, missing=1)

    # ---------------- A. octonion tables, both conventions
    TB, cb, xb = make_table(BAEZ, 8)
    TR, cr, xr = make_table(BRYANT, 8)
    TQ, cq, xq = make_table(QUAT, 4)
    for lab, T, comp, conf in (("baez", TB, cb, xb), ("bryant", TR, cr, xr), ("quat", TQ, cq, xq)):
        af, nz = assoc_profile(T)
        nf = norm_mult_fails(T, 300, 20261109)
        exp_nz = 168 if len(T) == 8 else 0
        check("table_%s_alternative_composition" % lab, comp and conf == 0 and af == 0 and nf == 0 and nz == exp_nz,
              complete=comp, conflicts=conf, assoc_antisym_failures=af, nonzero_basis_associators=nz,
              expect_nonzero=exp_nz, norm_mult_failures_300=nf)
    Tbad, _, _ = make_table([(2, 1, 4)] + BAEZ[1:], 8)
    af, _ = assoc_profile(Tbad)
    control("broken_fano_line_is_alternative", af == 0, assoc_antisym_failures=af)

    # signed-permutation isomorphisms baez -> bryant
    lines_b = [frozenset(l) for l in BAEZ]
    lines_r = set(frozenset(l) for l in BRYANT)
    isos = []
    for perm in itertools.permutations(range(1, 8)):
        pi = dict(zip(range(1, 8), perm))
        if not all(frozenset(pi[x] for x in l) in lines_r for l in lines_b):
            continue
        for signs in itertools.product((1, -1), repeat=7):
            sg = dict(zip(range(1, 8), signs))
            ok = True
            for i in range(1, 8):
                for j in range(1, 8):
                    if i == j:
                        continue
                    s, k = TB[i][j]
                    t, kk = TR[pi[i]][pi[j]]
                    if kk != pi[k] or sg[i] * sg[j] * t != s * sg[k]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                d = perm_sign(perm)
                for x in signs:
                    d *= x
                isos.append((perm, signs, d))
    dets = sorted(set(i[2] for i in isos))
    row("OCT_ISO_BAEZ_TO_BRYANT", count=len(isos), dets=dets,
        example_perm=(isos[0][0] if isos else "none"), example_signs=(isos[0][1] if isos else "none"),
        claim="MEASURED")
    check("baez_and_bryant_tables_isomorphic", len(isos) == 1344, count=len(isos), expect=1344)

    # ---------------- B. geometric Nijenhuis on S^6 (no Lie algebra)
    geo7b = {"nv": 7, "blocks": [(0, cross_tab(TB))]}
    geo7r = {"nv": 7, "blocks": [(0, cross_tab(TR))]}
    e7 = unit(6, 7)
    e1 = unit(0, 7)
    tang7 = [0, 1, 2, 3, 4, 5]
    tang1 = [1, 2, 3, 4, 5, 6]
    res = {}
    for lab, geo, pt, tg in (("S6_baez_at_e7", geo7b, e7, tang7), ("S6_baez_at_e1", geo7b, e1, tang1),
                             ("S6_bryant_at_e7", geo7r, e7, tang7)):
        for js in (1, -1):
            F, N, NJ, st = sphere_block(lab, geo, pt, tg, js, False)
            res[(lab, js)] = (N, NJ, st)
            check("geom_N_%s_js%d" % (lab, js),
                  st["vals_ok"] and st["tangent_ok"] and st["nz"] == 24 and st["zero"] == 12 and st["nk_formula"]
                  and st["nk_antisym"] and st["nk_strict"] and st["n4"],
                  nonzero=st["nz"], zero=st["zero"], expect="24_nonzero_12_zero_nearly_Kahler_strict_N_eq_-4J(nablaJ)")
        same = res[(lab, 1)][0] == res[(lab, -1)][0]
        check("geom_N_minusJ_equals_N_J_%s" % lab, same)
        zl = sorted(res[(lab, 1)][2]["zero_list"])
        Jpairs = set()
        for a in tg:
            Ja = pt_J(geo, pt, unit(a, 7), 1)
            for b in tg:
                if Ja == unit(b, 7) or Ja == vneg(unit(b, 7)) or a == b:
                    Jpairs.add((a, b))
        check("geom_N_zero_set_is_diagonal_plus_J_pairs_%s" % lab, set(zl) == Jpairs, zero_pairs=len(zl))
    N7, NJ7, _ = res[("S6_baez_at_e7", 1)]
    check("geom_N_equals_minus4_J_pi_cross_baez_e7", res[("S6_baez_at_e7", 1)][2]["n4pi"],
          law="N(u,v):-4*(e7 x pi_T(u x v))_with_J:+e7x")
    row("GEOM_N_SAMPLE", table="baez", point="e7", N_e1_e2=N7[(0, 1)][:6], N_e1_e3=N7[(0, 2)][:6],
        N_e2_e4=N7[(1, 3)][:6], N_e1_e4=N7[(0, 3)][:6], units="unit_tangent_vectors", claim="MEASURED")
    # compare with M1's Lie-algebraic N sample (d-scaled, d=2)
    m1p = os.path.join(BASE, "m1", "M1-OCTONION.hbp")
    m1s = None
    if os.path.exists(m1p):
        for ln in open(m1p, encoding="utf-8"):
            if ln.startswith("NIJENHUIS_SAMPLE|"):
                kv = dict(x.split("=", 1) for x in ln.strip().split("|")[1:] if "=" in x)
                m1s = {k: [int(t) for t in kv[k].strip("[]").split(",")] for k in
                       ("N_e1_e2", "N_e1_e3", "N_e2_e4", "N_e1_e4")}
    if m1s:
        mine = {"N_e1_e2": vscale(N7[(0, 1)][:6], 4), "N_e1_e3": vscale(N7[(0, 2)][:6], 4),
                "N_e2_e4": vscale(N7[(1, 3)][:6], 4), "N_e1_e4": vscale(N7[(0, 3)][:6], 4)}
        check("M1_lie_algebraic_N_equals_4x_geometric_N", mine == m1s, m1=";".join("%s:%s" % (k, m1s[k]) for k in sorted(m1s)),
              geometric_times_d2=";".join("%s:%s" % (k, mine[k]) for k in sorted(mine)),
              meaning="M1_formula_sign_and_m_projection_agree_with_vector_field_Nijenhuis")
    else:
        check("M1_sample_readable", False, missing=1)
    # tensoriality (a wrong formula would fail this)
    fv, tb = tensoriality(geo7b, e7, tang7, 1, 0, -1)
    check("geom_N_is_tensorial_f_eq_3_plus_x1", tb == 0, f_at_point=fv, failures=tb, pairs=36)
    fv, tb2 = tensoriality(geo7b, e7, tang7, 1, 0, 1)
    control("broken_Nijenhuis_sign_is_tensorial", tb2 == 0, failures=tb2, change="last_term_+[X,Y]")
    control("S6_N_vanishes_everywhere", res[("S6_baez_at_e7", 1)][2]["nz"] == 0,
            nonzero=res[("S6_baez_at_e7", 1)][2]["nz"])
    # integrable control: S2 x S2 x S2 inside Im H ^3 (real dimension 6)
    CRq = cross_tab(TQ)
    geo9 = {"nv": 9, "blocks": [(0, CRq), (3, CRq), (6, CRq)]}
    pt9 = [0, 0, 1, 0, 0, 1, 0, 0, 1]
    tg9 = [0, 1, 3, 4, 6, 7]
    F9, N9, NJ9, st9 = sphere_block("S2xS2xS2_quat", geo9, pt9, tg9, 1, True)
    check("control_S2cubed_integrable_N_zero_Kahler", st9["nz"] == 0 and not st9["nk_strict"] and st9["tangent_ok"]
          and st9["vals_ok"], nonzero=st9["nz"], nablaJ_nonzero=st9["nk_strict"],
          meaning="same_code_returns_0_on_an_integrable_6dim_product_of_quaternionic_S2")
    fv9, tb9 = tensoriality(geo9, pt9, tg9, 1, 0, -1)
    check("geom_N_tensorial_on_S2cubed", tb9 == 0, failures=tb9)

    # ---------------- C. Lambda^2 split, both tables; chirality
    phis = {}
    for lab, T in (("baez", TB), ("bryant", TR)):
        phi, o, Tm, st = g2_split(lab, T)
        phis[lab] = (phi, o, Tm)
        check("lambda2_split_7_14_%s" % lab, st["alt"] == 0 and st["bdiag"] and st["sym"] and st["minpoly"]
              and st["d7"] == 7 and st["d14"] == 14, orientation=o)
        mats, rk, ok, neq, lb = der_basis(T)
        in14 = all([sum(Tm[i][j] * v[j] for j in range(21)) for i in range(21)] == vneg(v)
                   for v in (der_as_2form(D) for D in mats))
        r14 = rank([der_as_2form(D) for D in mats], 21)
        check("der_%s_dim14_in_eig_minus1" % lab, len(mats) == 14 and rk == 50 and ok and lb == 0 and in14 and r14 == 14,
              equations=neq, rank=rk, dim=len(mats), leibniz_failures=lb, in_eig_minus1=in14, rank_as_2forms=r14)
        Trev = mscale(Tm, -1)
        inrev = all([sum(Trev[i][j] * v[j] for j in range(21)) for i in range(21)] == v
                    for v in (der_as_2form(D) for D in mats))
        row("ORIENTATION_CONVENTION", table=lab, reversed_orientation_eigs="-2(x7),+1(x14)",
            der_in_eig_plus1_when_reversed=inrev,
            meaning="the_7_14_split_is_convention_free;the_signs_2_and_-1_need_phi_positive_orientation",
            claim="MEASURED")
        check("orientation_flip_relabels_eigs_%s" % lab, inrev)
        ch = chirality(lab, T, phi, o, Tm)
        allasd = all(c[5] == 3 and c[4] == 0 and c[6] and c[3] for c in ch)
        row("COASSOC_CHIRALITY", table=lab, planes=len(ch),
            per_plane=";".join("A%d%d%d:phiA%+d,epsC%+d,orient_agree%d,SD_in14:%d,ASD_in14:%d" % (
                c[0][0], c[0][1], c[0][2], c[1], c[2], int(c[3]), c[4], c[5]) for c in ch),
            ASD_wedge_psi_zero=all(c[6] for c in ch), SD_wedge_psi_zero=any(c[7] for c in ch), claim="MEASURED")
        check("coassoc_ASD_half_in_lambda2_14_%s" % lab, allasd, rule="C_oriented_by_psi_on_C")
        if lab == "baez":
            DB = mats
    check("iso_det_equals_product_of_positive_orientations",
          dets == [phis["baez"][1] * phis["bryant"][1]], dets=dets,
          o_baez=phis["baez"][1], o_bryant=phis["bryant"][1])

    # ---------------- D. stabiliser, Higgs Gram, Casimir, commutant (baez)
    D = DB
    rows7 = [[Dm[k][7] for Dm in D] for k in range(8)]
    nd = len(D)
    hc, _, okh = nullspace(rows7, nd)

    def comb(cs):
        M = zeros(8)
        for c, Dm in zip(cs, D):
            if c:
                M = madd(M, Dm, c)
        return M
    H = [comb(c) for c in hc]
    trrows = [[tr(mmul(Dm, Hm)) for Dm in D] for Hm in H]
    mc, _, okm = nullspace(trrows, nd)
    Mm = [comb(c) for c in mc]
    gram_ok = all(3 * dot([X[k][7] for k in range(8)], [Y[k][7] for k in range(8)]) == -tr(mmul(X, Y))
                  for X in Mm for Y in Mm)
    h_kills = all(Hm[k][7] == 0 for Hm in H for k in range(8))
    check("stabiliser_e7_dim8_m_dim6", len(H) == 8 and len(Mm) == 6 and okh and okm and h_kills,
          h=len(H), m=len(Mm))
    check("higgs_gram_3_vv_equals_minus_trace_on_m", gram_ok, identity="3*<Xe7,Ye7>:-tr(XY)",
          meaning="6_broken_masses_equal;8_massless")
    G = [[-tr(mmul(A, B)) for B in D] for A in D]
    dG = det_bareiss(G)
    adjG = adjugate(G)
    Cas = zeros(8)
    for a in range(nd):
        for b in range(nd):
            if adjG[a][b]:
                Cas = madd(Cas, mmul(D[a], D[b]), adjG[a][b])
    target = mscale([[1 if (i == j and i > 0) else 0 for j in range(8)] for i in range(8)], -2 * dG)
    check("casimir_on_ImO_is_minus2_identity", Cas == target and dG > 0, det_gram=dG,
          identity="sum_ab_adj(G)_ab_D_aD_b:-2det(G)I7")
    J6 = [[0] * 6 for _ in range(6)]
    for l in range(1, 7):
        s, k = TB[7][l]
        J6[k - 1][l - 1] = s
    H6 = [[[Hm[i][j] for j in range(1, 7)] for i in range(1, 7)] for Hm in H]
    eq = []
    for Hm in H6:
        for i in range(6):
            for j in range(6):
                r = [0] * 36
                for k in range(6):
                    r[k * 6 + j] += Hm[i][k]
                    r[i * 6 + k] -= Hm[k][j]
                eq.append(r)
    cb_, _, okc = nullspace(eq, 36)
    has_I = rank(cb_ + [flat(eye(6))], 36) == len(cb_)
    has_J = rank(cb_ + [flat(J6)], 36) == len(cb_)
    eq2 = list(eq)
    for i in range(6):
        for j in range(6):
            r = [0] * 36
            for k in range(6):
                r[i * 6 + k] += J6[k][j]
                r[k * 6 + j] += J6[i][k]
            eq2.append(r)
    anti, _, _ = nullspace(eq2, 36)
    check("commutant_su3_on_T_is_span_I_J_complex_type", len(cb_) == 2 and has_I and has_J and len(anti) == 0
          and mmul(J6, J6) == mscale(eye(6), -1) and all(comm(J6, Hm) == zeros(6) for Hm in H6),
          commutant_dim=len(cb_), antilinear_intertwiners=len(anti),
          meaning="T_is_irreducible_complex_type;(T,J)_not_iso_(T,-J);3_vs_3bar_label_is_a_choice_of_iso_h:su(3)")
    # antipodal automorphism and Cartan
    phi_inv = [1, 1, 1, -1, 1, -1, -1, -1]
    aut_ok = True
    for i in range(8):
        for j in range(8):
            s, k = TB[i][j]
            if phi_inv[i] * phi_inv[j] * s != s * phi_inv[k]:
                aut_ok = False
    detphi = 1
    for x in phi_inv[1:]:
        detphi *= x
    P8 = [[phi_inv[i] if i == j else 0 for j in range(8)] for i in range(8)]
    J8 = zeros(8)
    for l in range(1, 7):
        s, k = TB[7][l]
        J8[k][l] = s
    flipsJ = mmul(mmul(P8, J8), P8) == mscale(J8, -1)

    def R(p, q):
        M = zeros(8)
        M[q][p] = 1
        M[p][q] = -1
        return M
    t1 = madd(R(1, 3), R(2, 6), -1)
    t2 = madd(R(2, 6), R(4, 5), -1)
    flatD = [flat(Dm) for Dm in D]
    t_in_der = nd == 14 and rank(flatD + [flat(t1), flat(t2)], 64) == nd
    t_in_h = all(t[k][7] == 0 for t in (t1, t2) for k in range(8))
    t_comm = comm(t1, t2) == zeros(8)
    adneg = mmul(mmul(P8, t1), P8) == mscale(t1, -1) and mmul(mmul(P8, t2), P8) == mscale(t2, -1)
    cen_eq = []
    for t in (t1, t2):
        for i in range(8):
            for j in range(8):
                cen_eq.append([comm(Dm, t)[i][j] for Dm in D])
    cen, _, _ = nullspace(cen_eq, nd)
    check("antipodal_automorphism_reverses_J_and_acts_minus1_on_cartan",
          aut_ok and detphi == 1 and flipsJ and t_in_der and t_in_h and t_comm and adneg and len(cen) == 2,
          phi="diag(1,1,-1,1,-1,-1,-1)_on_e1..e7", det=detphi, e7_to_minus_e7=1, PJP_eq_minusJ=flipsJ,
          cartan_rank=len(cen), Ad_on_t="-1")
    m1on = [1] + [-1] * 7
    aut_m1 = all(m1on[i] * m1on[j] * TB[i][j][0] == TB[i][j][0] * m1on[TB[i][j][1]]
                 for i in range(8) for j in range(8))
    control("minus1_on_ImO_is_automorphism", aut_m1)
    control("R17_is_derivation", rank(flatD + [flat(R(1, 7))], 64) == nd)

    # ---------------- E. W(G2) in simple-root coordinates
    GR = [[2, -3], [-3, 6]]

    def ip(u, v):
        return sum(u[i] * GR[i][j] * v[j] for i in range(2) for j in range(2))

    def refl(beta):
        bb = ip(beta, beta)
        cols = []
        for e in ([1, 0], [0, 1]):
            c = ediv(2 * ip(e, beta), bb)
            cols.append([e[0] - c * beta[0], e[1] - c * beta[1]])
        return ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))

    def m2mul(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)) for i in range(2))

    def m2app(A, v):
        return (A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1])
    I2 = ((1, 0), (0, 1))
    W = closure([refl([1, 0]), refl([0, 1])], m2mul, I2)
    roots = set()
    for w in W:
        roots.add(m2app(w, (1, 0)))
        roots.add(m2app(w, (0, 1)))
    longs = [r for r in roots if ip(r, r) == 6]
    shorts = [r for r in roots if ip(r, r) == 2]

    def order(A):
        P = A
        k = 1
        while P != I2:
            P = m2mul(P, A)
            k += 1
        return k
    hist = {}
    for w in W:
        o_ = order(w)
        hist[o_] = hist.get(o_, 0) + 1
    WA2 = closure([refl(list(r)) for r in longs], m2mul, I2)
    mI = ((-1, 0), (0, -1))
    negWA2 = set(m2mul(mI, w) for w in WA2)
    orbs = []
    left = set(shorts)
    while left:
        x = next(iter(left))
        o_ = set(m2app(w, x) for w in WA2)
        orbs.append(o_)
        left -= o_
    swap = all(set((-r[0], -r[1]) for r in orbs[0]) == orbs[1] for _ in [0]) if len(orbs) == 2 else False
    long_closed = all((a[0] + b[0], a[1] + b[1]) in set(longs) for a in longs for b in longs
                      if (a[0] + b[0], a[1] + b[1]) in roots)
    check("weyl_G2_D6_and_A2_inside", len(W) == 12 and len(roots) == 12 and len(longs) == 6 and len(shorts) == 6
          and hist == {1: 1, 2: 7, 3: 2, 6: 2} and len(WA2) == 6 and mI in W and mI not in WA2
          and (WA2 | negWA2) == W and len(orbs) == 2 and all(len(o_) == 3 for o_ in orbs) and swap and long_closed,
          order=len(W), order_hist=";".join("%d:%d" % kv for kv in sorted(hist.items())), WA2=len(WA2),
          short_orbits_under_WA2=len(orbs), minus_I_swaps_3_3bar=swap)

    # ---------------- F. Delta(6n^2)
    delta = {}
    for n in range(1, 8):
        M = 2 * n
        a = ((0, 1, 2), (2 % M, (-2) % M, 0))
        c = ((1, 2, 0), (0, 0, 0))
        b = ((0, 2, 1), (n % M, n % M, n % M))
        ident = ((0, 1, 2), (0, 0, 0))

        def mul(x, y, M=M):
            return mono_mul(x, y, M)
        Gd = closure([a, b, c], mul, ident)
        det1 = all((sum(g[1]) + (n if perm_sign(g[0]) < 0 else 0)) % M == 0 for g in Gd)
        diag = [g for g in Gd if g[0] == (0, 1, 2)]
        perms = set(g[0] for g in Gd)
        grid = all(all(x % 2 == 0 for x in g[1]) and (sum(g[1]) % M == 0) for g in diag)
        delta[n] = (Gd, mul, ident, M)
        check("delta_6n2_n%d" % n, len(Gd) == 6 * n * n and det1 and len(diag) == n * n and len(perms) == 6 and grid,
              order=len(Gd), expect=6 * n * n, diag_order=len(diag), weyl_image=len(perms), det1=det1)
    bbad = ((0, 2, 1), (0, 0, 0))
    G3 = delta[3]
    detb = (sum(bbad[1]) + (3 if perm_sign(bbad[0]) < 0 else 0)) % 6 == 0
    control("transposition_with_plus_entries_has_det1", detb)
    for n in range(1, 6):
        Gd, mul, ident, M = delta[n]

        def inv(x, M=M):
            return mono_inv(x, M)
        cls, found = abelian_normal_subgroups(Gd, mul, inv, ident)
        best = max(len(S) for S in found)
        nbest = sum(1 for S in found if len(S) == best)
        diagset = frozenset(g for g in Gd if g[0] == (0, 1, 2))
        check("delta_n%d_max_abelian_normal" % n, best == (3 if n == 1 else n * n) and
              (n == 1 or diagset in found), classes=len(cls), abelian_normal=len(found), max_order=best,
              count_max=nbest, index=ediv(len(Gd), best))
    flats = {}
    for n in range(1, 5):
        Gd, mul, ident, M = delta[n]
        flats[n] = flat_sector_count(Gd, mul)
    check("flat_sectors_on_T2_S3_is_8_S4_is_21", flats[1] == 8 and flats[2] == 21,
          counts=";".join("Delta(%d):%d" % (6 * n * n, flats[n]) for n in sorted(flats)),
          meaning="Hom(Z2,G)_mod_G:number_of_flat_connections_on_2-torus:torus_GSD_of_D(G)")
    Gs3 = delta[1][0]
    cls3 = classes_of(Gs3, delta[1][1], lambda x: mono_inv(x, 2))
    k_cent = 0
    for C in cls3:
        g = next(iter(C))
        Cg = set(h for h in Gs3 if delta[1][1](g, h) == delta[1][1](h, g))
        k_cent += len(classes_of(Cg, delta[1][1], lambda x: mono_inv(x, 2)))
    check("S3_flat_count_equals_irreps_of_quantum_double", k_cent == 8 == flats[1], sum_over_classes_k_centraliser=k_cent,
          finite_repo_readme_says="eight-state_flat-vacuum_cluster;D(S3)_torus_ground_space")
    # which S3 character is the plaquette weight? Delta(6) in SU(3) vs the 2-dim irrep
    tr3 = {}
    for g in Gs3:
        fixed = sum(1 for i in range(3) if g[0][i] == i)
        cl = 0 if fixed == 3 else (1 if fixed == 1 else 2)
        tr3[cl] = sum((-1) ** g[1][i] for i in range(3) if g[0][i] == i)
    chi_triv, chi_sign, chi_std = (1, 1, 1), (1, -1, 1), (2, 0, -1)
    t3 = (tr3[0], tr3[1], tr3[2])
    is_sign_plus_std = t3 == tuple(a + b for a, b in zip(chi_sign, chi_std))
    trW = {}
    for w_ in WA2:
        o_ = order(w_)
        trW[o_] = w_[0][0] + w_[1][1]
    w_char = (trW[1], trW[2], trW[3])
    check("plaquette_character_distinction", is_sign_plus_std and t3 != chi_std and w_char == chi_std,
          SU3_fundamental_trace_on_Delta6_classes_e_t_c=t3, decomposes_as="sign+std",
          WA2_reflection_rep_on_cartan_trace=w_char, equals="std_2dim_irrep",
          meaning="README_weight_2-chi(U_p)_needs_chi(e):2_so_it_is_the_2dim_irrep:Weyl_reflection_rep_on_t;Re_tr3_restricted_to_Delta6_is_sign+std")
    control("SU3_trace_on_Delta6_equals_2dim_irrep", t3 == chi_std, t3=t3)

    # ---------------- G. binary icosahedral 2I
    PHI = (0, 1)
    PHIINV = (-1, 1)
    I120 = icosians([(0, 0), (1, 0), PHIINV, PHI])
    norm_ok = all(q_norm4(x) == (4, 0) for x in I120)
    closed = all(q2mul(x, y) in I120 for x in I120 for y in I120)
    hist = {}
    for x in I120:
        k = q_order(x)
        hist[k] = hist.get(k, 0) + 1
    cls2 = classes_of(I120, q2mul, q_conj)
    idc = [i for i, cc in enumerate(cls2) if ONE_Q in cc][0]
    others = [i for i in range(len(cls2)) if i != idc]
    normals = []
    for r in range(len(others) + 1):
        for sub in itertools.combinations(others, r):
            S = set(cls2[idc]).union(*[cls2[i] for i in sub])
            if all(q2mul(x, y) in S for x in S for y in S):
                normals.append(len(S))
    check("binary_icosahedral_120", len(I120) == 120 and norm_ok and closed and
          hist == {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24} and len(cls2) == 9 and sorted(normals) == [1, 2, 120],
          order=len(I120), closed=closed, classes=len(cls2), normal_subgroup_orders=sorted(normals),
          max_abelian_normal=2, jordan_index=60)
    Iwrong = icosians([(0, 0), (1, 0), PHI, PHIINV])
    closed_w = True
    for x in Iwrong:
        for y in Iwrong:
            try:
                if q2mul(x, y) not in Iwrong:
                    closed_w = False
                    break
            except ArithmeticError:
                closed_w = False
                break
        if not closed_w:
            break
    common = len(I120 & Iwrong)
    row("CORRECTION", claim="MEASURED", text="First run used the odd-permutation base (0,1,phi,phi^-1) as a broken control expecting closure to fail; measured: that set is also a closed 120-element group (the Galois conjugate copy, phi to -phi^-1), so the premise was wrong and the control was replaced by a removed-element control. Retained beside the claim.",
        odd_base_closed=closed_w, odd_base_size=len(Iwrong), common_with_2I=common)
    check("second_icosian_copy_is_closed_and_shares_hurwitz_units", closed_w and len(Iwrong) == 120 and common == 24,
          common=common)
    holed = set(I120)
    holed.discard(((1, 0), (1, 0), (1, 0), (1, 0)))
    holed_closed = all(q2mul(x, y) in holed for x in holed for y in holed)
    control("2I_minus_one_element_closed", holed_closed, size=len(holed))
    g2m3 = ((0, 1), (-1, 1), (1, 0), (0, 0))
    g1m3 = ((1, 0), (1, 0), (1, 0), (1, 0))
    gen = closure([g1m3, g2m3], q2mul, ONE_Q)
    check("M3_generators_in_2I_and_generate_it", g2m3 in I120 and g1m3 in I120 and gen == I120, generated=len(gen))
    # C8 counterexample ingredient: 2I has no abelian normal subgroup of index <= 6
    check("2I_not_conjugate_into_N(T)_anywhere", len(I120) == 120 and sorted(normals) == [1, 2, 120], claim="MEASURED+NAMED",
          argument="normal_subgroups_of_2I_are_1,{+-1},2I;so_H:U(1)*2I_in_S(U(2)xU(1))_has_H/H0:2I/{+-1}_with_no_nontrivial_abelian_normal_subgroup;every_abelian_normal_A_of_H_lies_in_H0;[H:A]>=60>6:|W(SU3)|;H_is_a_Hausdorff_limit_of_finite_Z_m*2I_and_is_not_inside_any_conjugate_of_N(T)")

    # ---------------- H. Perron toy: positive kernel => finite gap, always
    S3 = list(itertools.permutations(range(3)))

    def pmul(p, q):
        return tuple(q[p[i]] for i in range(3))

    def pinv(p):
        r = [0] * 3
        for i in range(3):
            r[p[i]] = i
        return tuple(r)

    def pcls(p):
        fixed = sum(1 for i in range(3) if p[i] == i)
        return 0 if fixed == 3 else (1 if fixed == 1 else 2)
    def gapped(w):
        we, wt, wc = w
        Mx = [[w[pcls(pmul(g, pinv(h)))] for h in S3] for g in S3]
        lt = we + 3 * wt + 2 * wc
        ls = we - 3 * wt + 2 * wc
        ld = we - wc
        Pm = mmul(mmul(madd(Mx, eye(6), -lt), madd(Mx, eye(6), -ls)), madd(Mx, eye(6), -ld))
        return Pm == zeros(6) and lt > abs(ls) and lt > abs(ld)
    good = 0
    total = 0
    for we in range(1, 11):
        for wt in range(1, 11):
            for wc in range(1, 11):
                total += 1
                if gapped((we, wt, wc)):
                    good += 1
    check("toy_S3_link_positive_kernel_always_gapped", good == total == 1000, grid=total, gapped=good,
          meaning="Perron-Frobenius:finite_positive_transfer_kernel_gives_finite_gap_for_every_positive_weight")
    control("identity_kernel_has_positive_gap", gapped((1, 0, 0)), weights="1,0,0_not_strictly_positive")

    # ---------------- I. claims verdict rows
    V = [
        ("C1", "SUPPORTED_NAMED", "Measured: der(O) dim 14 (rank 50 of 64, two tables), stab(e7) dim 8 with commutant span(I,J) on T, m dim 6, orbit map rank 6. Group level (transitivity, stabiliser is SU(3) not just su(3)) uses G2 connected and S^6 simply connected: NAMED, Harvey 1990 Spinors and Calibrations ch.6; Baez 2002 Bull AMS 39 sec 4.1; Bryant 1987 Ann Math 126.", "none"),
        ("C2", "SUPPORTED_MEASURED", "J^2:-1, orthogonal, commutes with su(3); commutant span(I,J) and zero antilinear intertwiners, so (T,J) is an irreducible complex 3-dim rep of complex type. Whether it is called 3 or 3bar is a choice of isomorphism h:su(3): with (T,J) as defining rep it is the 3; M2 labelling calls +J the 3bar.", "T_xS^6 with J_x is an irreducible complex 3-dimensional SU(3) representation of complex type; (T,J) and (T,-J) are the 3 and 3bar in some order fixed by convention."),
        ("C3", "SUPPORTED_MEASURED", "Nonintegrability measured two independent ways: M1 Lie-algebraic N plus Koszul (h_C+m10 not closed, [m10,m10] in m01 by eps), and this run's ambient vector-field Nijenhuis on S^6 with no Lie algebra (24 of 36 ordered pairs nonzero, N:-4J(nabla J), equal for +J and -J, both octonion tables, two points), matching M1's sample exactly (4x geometric). Control S2xS2xS2 gives N:0. CP^2 clause: the integrable reason is that m10 is an abelian subalgebra ([m10,m10]:0, g has a Z-grading g_-1+g_0+g_1), while G2/SU(3) has a Z3-grading with [g_1,g_1] in g_2:3bar. Hopf problem open: NAMED, Agricola-Bazzoni-Goertsches-Konstantis-Rollenske 2018 Diff Geom Appl 57; no ROUND-metric-orthogonal integrable J: LeBrun 1987 Proc AMS 101.", "Replace 'because the 2 there is abelian' with 'because m^(1,0):C^2 is an abelian subalgebra, so [m10,m10]:0 lies in h_C+m10'."),
        ("C4", "SUPPORTED_NAMED", "Measured pointwise on S^6 (two tables, two points): (nabla_X J)Y:pi_T(X x Y) antisymmetric so (nabla_X J)X:0, and nabla J nonzero, i.e. strict nearly Kahler, using the Gauss-formula connection (NAMED). Cone statement NAMED: Gray 1976 Math Ann 223; Fukami-Ishihara 1955 Tohoku Math J 7; Bar 1993 Comm Math Phys 154 (nearly Kahler 6-manifolds iff cone holonomy in G2).", "none"),
        ("C5", "SUPPORTED_MEASURED", "Measured here and in M2: |W|:12 with order histogram 1:1,2:7,3:2,6:2 (D6:S3xZ2), W(A2) from long roots order 6 without -I, W:W(A2)x{+-I}, short roots two W(A2)-orbits of 3 swapped by -I; diag(1,1,-1,1,-1,-1,-1) is an automorphism with det +1, sends e7 to -e7, conjugates J to -J and acts as -1 on the Cartan.", "none"),
        ("C6", "SUPPORTED_MEASURED", "Tree-level algebra measured: stabiliser dim 8, 3<Xv,Yv>:-tr(XY) on the 6-dim m so the six broken masses are equal, Casimir -2 on R^7. Vacuum manifold S^6 relies on C1 (NAMED). Mass formula NAMED: Weinberg QFT vol II sec 21.1. Degeneracy is forced by Schur (one invariant symmetric form on an irreducible complex-type m); J supplies the complex structure that makes m the complex 3.", "Scope: classical tree-level Higgs mechanism; no statement about the quantum theory."),
        ("C7", "SUPPORTED_MEASURED", "Measured in both tables: T:*(phi^.) has eigenvalues 2 (dim 7) and -1 (dim 14) when phi is positive for the orientation, der(O) is the -1 eigenspace, and on all 7 coassociative planes the ASD half (C oriented by psi) lies in Lambda2_14 with F^psi:0. Signs 2/-1 and the label ASD flip with orientation conventions; 'one chirality' is the convention-free statement. G2-instanton equivalence NAMED: Donaldson-Thomas 1998 (Gauge theory in higher dimensions, Geometric Universe, OUP); Corrigan-Devchand-Fairlie-Nuyts 1983 Nucl Phys B 214. Pulled-back BPST has infinite action on R^7.", "none"),
        ("C8", "OVERSTATED", "The Jordan-limit argument (M3 LIMIT steps 1-8) is sound: Hausdorff/Chabauty limit H of finite subgroups has an abelian closed subgroup of finite index, so H0 is a torus of dim <= N-1. 'At most T x| W' does not follow as containment: F_m:Z_m*2I in S(U(2)xU(1)) in SU(3) converges to H:U(1)*2I, with H/H0 of order 60 and every abelian normal subgroup inside H0, so H lies in no conjugate of N(T) (whose abelian normal T has index 6). Delta(6n^2):(Z_n^2)x|S3 in N(T) measured n:1..7.", "Every Hausdorff limit of finite subgroups of SU(N) is torus-by-finite with identity component a torus of dimension at most N-1; N(T):T x| W attains that dimension (limit of Delta(6n^2)); limits need not lie in N(T) (example U(1)*2I in SU(3)); so no tower of finite subgroups approaches SU(N)."),
        ("C9", "OVERSTATED", "S3 in the finite repo is Delta(6):{sgn(s)P_s} in N(T), which maps onto W(SU(3)) and meets T trivially (measured diag_order 1). Maximal-abelian projection ('t Hooft 1981 Nucl Phys B 190; Kronfeld-Schierholz-Wiese 1987 Nucl Phys B 293) keeps the torus U(1)^2 and leaves W as residual symmetry; it does not leave S3. Per the repo README (worktree head 6ec2a808, scoped, not the owning gate) the certified gap is E9-E8 above the 8 flat S3 states on 2D spatial tori for r:J_B/J_E>:9, the D(S3) quantum-double phase; this run measures 8:|Hom(Z^2,S3)/S3|:number of D(S3) irreps, and 21 for Delta(24):S4, growing with n (84 for Delta(54), 103 for Delta(96)). The README plaquette weight 2-chi(U_p) uses the 2-dim irrep (the W(A2) reflection character on the Cartan, measured), whereas Re tr_3 of SU(3) restricted to Delta(6) is sign+std (measured), so the model is an S3 gauge theory with its own action, not the SU(3) Wilson action restricted to a subgroup. A finite positive transfer kernel is gapped for every positive weight (toy measured 1000/1000; Perron-Frobenius NAMED). Geometry of G2/SU(3) supplies no gap.", "The S3 lattice theory has gauge group Delta(6):S3 in N(T) in SU(3), isomorphic to W(SU(3)) in W(G2); its certified gaps are statements about finite 2D spatial tori in the D(S3) topological phase above 8 flat sectors; transfer to continuum SU(N) in 3+1 dimensions needs a proof that is not yet available and is OPEN."),
        ("C10", "OVERSTATED", "Routes are real and NAMED: Wilson 1974 PRD 10 2445; Kogut-Susskind 1975 PRD 11 395; q-deformed SU(2)_k: Zache-Gonzalez-Cuadra-Zoller 2023 PRL 131 171902; non-subgroup partitions: Hartung-Jakobs-Jansen-Ostmeyer-Urbach 2022 EPJC 82 237; quantum links: Chandrasekharan-Wiese 1997 Nucl Phys B 492 455. The Clay statement (Jaffe-Witten 2000, CMI) is on R^4: existence for every compact simple G with OS/Wightman axioms and a gap Delta>0 of the infinite-volume Hamiltonian. A gap at fixed physical size is one step; the gap must also stay bounded below as volume goes to infinity, with nontriviality and restored O(4). Closest rigorous fixed-volume continuum results: Balaban 1987-89 CMP (UV stability, 4D lattice YM on a finite torus); Magnen-Rivasseau-Seneor 1993 CMP 155 (YM4 with an infrared cutoff).", "Clay asks for YM on R^4 for each compact simple G, satisfying OS/Wightman axioms incl. Euclidean invariance and nontriviality, with an infinite-volume mass gap. A lattice route needs continuum limit a->0 with asymptotic-freedom renormalisation, OS reconstruction, and a gap bounded below uniformly in both a->0 and L->infinity; the fixed-physical-size finite-universe target is a separate, well-posed statement whose transfer to R^4 is OPEN."),
    ]
    for cid, status, why, restate in V:
        row("CLAIM", id=cid, status=status, why=why, restate=restate,
            claim=("MEASURED" if status == "SUPPORTED_MEASURED" else "NAMED"))
    row("LIMIT_REVIEW", steps="M3_LIMIT_1-8", verdict="sound", notes="Blaschke selection on compact SU(N); finite coset cover passes to limit with bi-invariant metric; limit cosets g_iA cover H with g_i in H; A closed abelian of finite index hence open hence contains H0; H0 torus dim<=N-1; uniform eps_N by contradiction",
        gap_found="none_in_steps_1-8;synthesis_C8_overreads_to_containment_in_T_x|_W", claim="NAMED")


def finalize():
    ok = not FAILS
    row("SUMMARY", checks_passed=NPASS[0], controls=NCTL[0], failed=len(FAILS),
        failed_names=(";".join(FAILS) or "none"), verdict=("ALL_PASS" if ok else "FAIL"))
    body = ("\n".join(ROWS) + "\n").encode("utf-8")
    with open(OUT, "wb") as fh:
        fh.write(body)
    digest = hashlib.sha256(body).hexdigest()
    with open(SIDE, "wb") as fh:
        fh.write(("%s  %s\n" % (digest, OUT_NAME)).encode("utf-8"))
    for r in ROWS:
        if "status=FAIL" in r:
            print(r)
    print("HBP_SHA256|sha256=%s|rows=%d|failed=%d|json=0" % (digest, len(ROWS), len(FAILS)))
    sys.exit(0 if ok else 1)


def main():
    try:
        main_body()
    except Exception as exc:
        FAILS.append("crash")
        row("CRASH", status="FAIL", claim="MEASURED", error=type(exc).__name__, detail=str(exc)[:160])
    finalize()


main()
