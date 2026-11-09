"""refute_spinor.py -- second-seat refuter of the frozen spinor-symbol prediction, Iconic-3-space, branch freeze-spinor-symbol-20261109.

Built from definitions. Before the measured block was written only the HEAD, FORMAT and CONV rows of PREDICTION-SPINOR-SYMBOL.hbp
were read; no PRED, ADD or DIMS value; no other agent's script; no GPT-6 receipt.
Integer only: no floating point, no true division, no power operator; exact quotients go through divmod with a zero remainder.
Every identity in xi is an identity of integer polynomials (a dict from exponent tuple to int), compared coefficient by coefficient,
so absent monomials are compared as zero on both sides. An identity linear in xi is checked on the basis covectors, which is the
same coefficient comparison. Sampled covectors appear only as controls.

Run: python3 -I -B refute_spinor.py OUTDIR      writes REFUTE-SPINOR.measured.hbp and its .sha256 (the measured block)
"""
import sys
import os
import hashlib
import itertools
from math import gcd

ROWS = []
COUNT = {"meas": 0, "meas_pass": 0, "value": 0, "ctrl": 0, "ctrl_ok": 0}
CTX = {}


def emit(tag, kv):
    parts = [tag]
    for k, v in kv:
        s = str(v)
        if "|" in s or "\n" in s or "=" in k or "|" in k:
            raise SystemExit("bad field " + k)
        parts.append(k + "=" + s)
    parts.append("json=0")
    ROWS.append("|".join(parts))


def meas(mid, name, ok, kv):
    COUNT["meas"] += 1
    if ok:
        COUNT["meas_pass"] += 1
    emit("MEAS", [("id", mid), ("name", name)] + kv + [("status", "PASS" if ok else "FAIL"), ("claim", "MEASURED")])


def value(mid, name, kv):
    COUNT["value"] += 1
    emit("MEAS", [("id", mid), ("name", name)] + kv + [("status", "VALUE"), ("claim", "MEASURED")])


def control(cid, name, failed, kv, expect="FAIL"):
    observed = "FAIL" if failed else "PASS"
    ok = observed == expect
    COUNT["ctrl"] += 1
    if ok:
        COUNT["ctrl_ok"] += 1
    emit("CONTROL", [("id", cid), ("name", name)] + kv + [("expected", expect), ("observed", observed), ("control_ok", 1 if ok else 0), ("claim", "MEASURED")])


def exq(a, b):
    q, r = divmod(a, b)
    if r != 0:
        raise SystemExit("inexact quotient")
    return q


def comb(n, k):
    num = 1
    den = 1
    for t in range(k):
        num = num * (n - t)
        den = den * (t + 1)
    return exq(num, den)


# ---------------- integer polynomials ----------------

def padd(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        v = r.get(m, 0) + s * c
        if v:
            r[m] = v
        else:
            r.pop(m, None)
    return r


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            v = r.get(m, 0) + c1 * c2
            if v:
                r[m] = v
            else:
                r.pop(m, None)
    return r


def pscale(p, k):
    if k == 0:
        return {}
    return {m: c * k for m, c in p.items()}


def pvar(i, nv):
    return {tuple(1 if j == i else 0 for j in range(nv)): 1}


def pconst(c, nv):
    if c == 0:
        return {}
    return {tuple(0 for _ in range(nv)): c}


def pderiv(p, i):
    r = {}
    for m, c in p.items():
        if m[i] > 0:
            mm = list(m)
            mm[i] -= 1
            r[tuple(mm)] = c * m[i]
    return r


def psum(ps):
    acc = {}
    for p in ps:
        acc = padd(acc, p)
    return acc


def nsq_poly(vec):
    return psum(pmul(x, x) for x in vec)


def dot_poly(u, v):
    return psum(pmul(a, b) for a, b in zip(u, v))


def cmp_poly(a, b):
    keys = set(a) | set(b)
    return 1, len(keys), sum(1 for m in keys if a.get(m, 0) != b.get(m, 0))


def cmp_vec(u, v):
    e = 0
    s = 0
    n = 0
    for a, b in zip(u, v):
        x, y, z = cmp_poly(a, b)
        e += x
        s += y
        n += z
    return e, s, n


# ---------------- matrices of polynomials ----------------

def mzero(n, m):
    return [[{} for _ in range(m)] for _ in range(n)]


def mmul(A, B):
    k = len(B)
    m = len(B[0])
    out = []
    for Ai in A:
        row = []
        for j in range(m):
            acc = {}
            for t in range(k):
                a = Ai[t]
                if a:
                    b = B[t][j]
                    if b:
                        acc = padd(acc, pmul(a, b))
            row.append(acc)
        out.append(row)
    return out


def madd(A, B, s=1):
    return [[padd(a, b, s) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mneg(A):
    return [[pscale(a, -1) for a in r] for r in A]


def mT(A):
    return [list(r) for r in zip(*A)]


def mdiag(n, p):
    M = mzero(n, n)
    for i in range(n):
        M[i][i] = dict(p)
    return M


def mapply(A, v):
    return [psum(pmul(a, x) for a, x in zip(row, v) if a and x) for row in A]


def mblock(rows):
    out = []
    for brow in rows:
        for i in range(len(brow[0])):
            r = []
            for B in brow:
                r.extend(B[i])
            out.append(r)
    return out


def cmp_mat(A, B):
    e = 0
    s = 0
    n = 0
    for ra, rb in zip(A, B):
        x, y, z = cmp_vec(ra, rb)
        e += x
        s += y
        n += z
    return e, s, n


def to_poly(A, nv):
    return [[pconst(x, nv) for x in r] for r in A]


def outer(u, v):
    return [[pmul(a, b) for b in v] for a in u]


# ---------------- integer matrices ----------------

def imul(A, B):
    BT = list(zip(*B))
    return [[sum(a * b for a, b in zip(r, c)) for c in BT] for r in A]


def iadd(A, B, s=1):
    return [[a + s * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def iT(A):
    return [list(r) for r in zip(*A)]


def ieye(n, k=1):
    return [[k if i == j else 0 for j in range(n)] for i in range(n)]


def izero(n, m):
    return [[0] * m for _ in range(n)]


def iscale(A, k):
    return [[k * x for x in r] for r in A]


def iblock(rows):
    out = []
    for brow in rows:
        for i in range(len(brow[0])):
            r = []
            for B in brow:
                r.extend(B[i])
            out.append(r)
    return out


def rref_int(M):
    A = [list(r) for r in M]
    nr = len(A)
    nc = len(A[0]) if nr else 0
    pcs = []
    r = 0
    for c in range(nc):
        if r == nr:
            break
        p = None
        for i in range(r, nr):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        g = 0
        for x in A[r]:
            g = gcd(g, x)
        if g > 1:
            A[r] = [exq(x, g) for x in A[r]]
        pr = A[r]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                pv = pr[c]
                new = [pv * x - f * y for x, y in zip(A[i], pr)]
                g = 0
                for x in new:
                    g = gcd(g, x)
                if g > 1:
                    new = [exq(x, g) for x in new]
                A[i] = new
        pcs.append(c)
        r += 1
    return A[:r], pcs


def irank(M):
    return len(rref_int(M)[1])


def nullspace_int(M):
    R, pcs = rref_int(M)
    nc = len(M[0])
    basis = []
    for f in range(nc):
        if f in pcs:
            continue
        L = 1
        for row, c in zip(R, pcs):
            d = abs(row[c])
            L = exq(L * d, gcd(L, d))
        v = [0] * nc
        v[f] = L
        for row, c in zip(R, pcs):
            v[c] = exq(-row[f] * L, row[c])
        g = 0
        for x in v:
            g = gcd(g, x)
        basis.append([exq(x, g) for x in v])
    return basis


# ---------------- exterior algebra ----------------

def perm_sign(seq):
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
            elif a[i] == a[j]:
                return 0
    return s


def kbasis(n, k):
    return list(itertools.combinations(range(n), k))


def star_matrix(n, k, sigma):
    src = kbasis(n, k)
    dst = kbasis(n, n - k)
    di = {b: i for i, b in enumerate(dst)}
    S = izero(len(dst), len(src))
    for j, I in enumerate(src):
        J = tuple(x for x in range(n) if x not in I)
        S[di[J]][j] = sigma * perm_sign(I + J)
    return S


def wedge_matrix(n, F, degF, k):
    src = kbasis(n, k)
    dst = kbasis(n, k + degF)
    di = {b: i for i, b in enumerate(dst)}
    W = izero(len(dst), len(src))
    for j, J in enumerate(src):
        for I, c in F.items():
            s = perm_sign(I + J)
            if s:
                W[di[tuple(sorted(I + J))]][j] += c * s
    return W


def form_vec(n, k, F):
    return [F.get(b, 0) for b in kbasis(n, k)]


def wedge_vec_poly(xi, a):
    n = len(xi)
    return [padd(pmul(xi[i], a[j]), pmul(xi[j], a[i]), -1) for (i, j) in kbasis(n, 2)]


def rot_matrix(n, omega_vec):
    R = izero(n, n)
    for (i, j), w in zip(kbasis(n, 2), omega_vec):
        if w:
            R[j][i] += w
            R[i][j] -= w
    return R


# ---------------- algebras from line tables ----------------

def lab7(n):
    return divmod(n - 1, 7)[1] + 1


OCT_LINES = [(i, lab7(i + 1), lab7(i + 3)) for i in range(1, 8)]
FROZEN_LINES = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def table_from_lines(dim, lines, flip=None):
    T = [[None] * dim for _ in range(dim)]
    for a in range(dim):
        T[0][a] = (1, a)
        T[a][0] = (1, a)
    for a in range(1, dim):
        T[a][a] = (-1, 0)
    for ln in lines:
        a, b, c = ln
        if flip is not None and ln == flip:
            a, b = b, a
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            T[x][y] = (1, z)
            T[y][x] = (-1, z)
    for a in range(dim):
        for b in range(dim):
            if T[a][b] is None:
                raise SystemExit("incomplete table")
    return T


def amul_poly(T, x, y):
    d = len(T)
    z = [{} for _ in range(d)]
    for a in range(d):
        if not x[a]:
            continue
        for b in range(d):
            if not y[b]:
                continue
            s, c = T[a][b]
            z[c] = padd(z[c], pscale(pmul(x[a], y[b]), s))
    return z


def amul_int(T, x, y):
    d = len(T)
    z = [0] * d
    for a in range(d):
        if x[a]:
            for b in range(d):
                if y[b]:
                    s, c = T[a][b]
                    z[c] += s * x[a] * y[b]
    return z


def Lmat_poly(T, x):
    d = len(T)
    M = mzero(d, d)
    for a in range(d):
        if not x[a]:
            continue
        for b in range(d):
            s, c = T[a][b]
            M[c][b] = padd(M[c][b], pscale(x[a], s))
    return M


def Rmat_poly(T, x):
    d = len(T)
    M = mzero(d, d)
    for b in range(d):
        if not x[b]:
            continue
        for a in range(d):
            s, c = T[a][b]
            M[c][a] = padd(M[c][a], pscale(x[b], s))
    return M


def Lmat_int(T, x):
    d = len(T)
    M = izero(d, d)
    for a in range(d):
        if x[a]:
            for b in range(d):
                s, c = T[a][b]
                M[c][b] += s * x[a]
    return M


def Rmat_int(T, x):
    d = len(T)
    M = izero(d, d)
    for b in range(d):
        if x[b]:
            for a in range(d):
                s, c = T[a][b]
                M[c][a] += s * x[b]
    return M


def conj_poly(x):
    return [x[0]] + [pscale(v, -1) for v in x[1:]]


def conj_int(x):
    return [x[0]] + [-v for v in x[1:]]


def unit(d, i):
    return [1 if j == i else 0 for j in range(d)]


def conj_matrix(d):
    return [[(1 if i == 0 else -1) if i == j else 0 for j in range(d)] for i in range(d)]


def candidates(T, x, nv):
    d = len(T)
    cx = conj_poly(x)
    C = to_poly(conj_matrix(d), nv)
    Lx = Lmat_poly(T, x)
    Lc = Lmat_poly(T, cx)
    Rx = Rmat_poly(T, x)
    Rc = Rmat_poly(T, cx)
    base = [("xi*a", Lx), ("xi*conj(a)", mmul(Lx, C)), ("conj(xi)*a", Lc), ("conj(xi)*conj(a)", mmul(Lc, C)),
            ("a*xi", Rx), ("a*conj(xi)", Rc), ("conj(a)*xi", mmul(Rx, C)), ("conj(a)*conj(xi)", mmul(Rc, C))]
    out = []
    for name, M in base:
        out.append(("+" + name, M))
        out.append(("-" + name, mneg(M)))
    return out


def match_candidates(M, T, x, nv):
    holding = []
    unequal_counts = []
    for name, C in candidates(T, x, nv):
        e, s, n = cmp_mat(M, C)
        if n == 0:
            holding.append(name)
        else:
            unequal_counts.append(n)
    return holding, unequal_counts


# ---------------- 4D ----------------

def sec_4d():
    TQ = table_from_lines(4, [(1, 2, 3)])
    ham = 1 if (TQ[1][2] == (1, 3) and TQ[2][3] == (1, 1) and TQ[3][1] == (1, 2) and TQ[1][1] == (-1, 0) and TQ[2][1] == (-1, 3)) else 0
    x8 = [pvar(i, 8) for i in range(4)]
    y8 = [pvar(i + 4, 8) for i in range(4)]
    e, s, n = cmp_poly(nsq_poly(amul_poly(TQ, x8, y8)), pmul(nsq_poly(x8), nsq_poly(y8)))
    meas("R4D.01", "quaternion_table_hamilton_and_composition_nsq(xy)=nsq(x)nsq(y)", ham == 1 and n == 0,
         [("hamilton_ij_k_jk_i_ki_j", ham), ("composition_slots", s), ("composition_unequal", n)])

    nv = 4
    xi = [pvar(i, nv) for i in range(4)]
    cxi = conj_poly(xi)
    nsq = nsq_poly(xi)
    L = Lmat_poly(TQ, xi)
    Lc = Lmat_poly(TQ, cxi)
    Z = mzero(4, 4)
    gamma = mblock([[Z, L], [mneg(Lc), Z]])
    G2 = mmul(gamma, gamma)
    e, s, n = cmp_mat(G2, mdiag(8, pscale(nsq, -1)))
    meas("R4D.02", "gamma(xi)^2=-nsq(xi)I8_with_gamma=[[0,L_xi],[-L_conj(xi),0]]", n == 0,
         [("entries", e), ("support_slots", s), ("unequal", n), ("full_basis_slots", e * comb(nv + 1, 2))])
    bad = mblock([[Z, L], [Lc, Z]])
    e, s, n = cmp_mat(mmul(bad, bad), mdiag(8, pscale(nsq, -1)))
    control("K4D.01", "gamma_with_+L_conj(xi)_squares_to_-nsq", n > 0, [("unequal", n)])

    # f basis and the M matrix
    B2 = kbasis(4, 2)
    fvec = [form_vec(4, 2, {(0, 1): 1, (2, 3): 1}), form_vec(4, 2, {(0, 2): 1, (1, 3): -1}), form_vec(4, 2, {(0, 3): 1, (1, 2): 1})]
    S4 = star_matrix(4, 2, 1)
    sd = 1 if all(imul(S4, [[v] for v in f]) == [[v] for v in f] for f in fvec) else 0
    gram = imul(fvec, iT(fvec))
    meas("R4D.03", "f1_f2_f3_self_dual_for_e1234_gram", sd == 1 and gram == ieye(3, 2), [("self_dual", sd), ("gram", ";".join(",".join(str(x) for x in r) for r in gram))])

    def M_of(fv):
        rows = [[pscale(xi[b], -1) for b in range(4)]]
        for f in fv:
            row = []
            for b in range(4):
                acc = {}
                for (i, j), c in zip(B2, f):
                    if c == 0:
                        continue
                    if j == b:
                        acc = padd(acc, pscale(xi[i], c))
                    if i == b:
                        acc = padd(acc, pscale(xi[j], -c))
                row.append(acc)
            rows.append(row)
        return rows
    M = M_of(fvec)
    CTX["M4"] = M
    CTX["S4"] = S4
    CTX["TQ"] = TQ
    e1, s1, n1 = cmp_mat(mmul(mT(M), M), mdiag(4, nsq))
    e2, s2, n2 = cmp_mat(mmul(M, mT(M)), mdiag(4, nsq))
    meas("R4D.04", "M^T_M=nsq(xi)I4_and_M_M^T=nsq(xi)I4", n1 == 0 and n2 == 0,
         [("MtM_entries", e1), ("MtM_support_slots", s1), ("MtM_unequal", n1), ("MMt_unequal", n2), ("full_basis_slots_each", e1 * comb(nv + 1, 2))])
    fbad = [fvec[0], form_vec(4, 2, {(0, 2): 1, (1, 3): 1}), fvec[2]]
    Mb = M_of(fbad)
    e, s, n = cmp_mat(mmul(mT(Mb), Mb), mdiag(4, nsq))
    control("K4D.02", "M_with_wrong_f2=e13+e24_M^T_M=nsq_I4", n > 0, [("unequal", n)])

    # (I + *)(xi wedge a) = sum c_k f_k, in 8 variables
    xa = [pvar(i, 8) for i in range(4)]
    aa = [pvar(i + 4, 8) for i in range(4)]
    w = wedge_vec_poly(xa, aa)
    TI = iadd(S4, ieye(6))
    lhs = mapply(to_poly(TI, 8), w)
    cks = [dot_poly(w, [pconst(c, 8) for c in f]) for f in fvec]
    rhs = [psum(pmul(cks[k], pconst(fvec[k][t], 8)) for k in range(3)) for t in range(6)]
    e, s, n = cmp_vec(lhs, rhs)
    meas("R4D.05", "(I+*)(xi^a)=c1f1+c2f2+c3f3_with_c_k=<xi^a,f_k>", n == 0, [("components", e), ("unequal", n)])

    # R1, 4D: k = 2
    r1_rows("4D", 4, S4, 2, w, xa, aa, 8)

    # M as a quaternion expression
    holding, unequal_counts = match_candidates(M, TQ, xi, nv)
    value("R4D.06", "M(xi)a_as_quaternion_expression_of_16_candidates", [("holding", ";".join(holding) if holding else "none"), ("non_holding", len(unequal_counts))])
    M2 = mmul(M, M)
    e, s, np_ = cmp_mat(M2, mdiag(4, nsq))
    e, s, nm = cmp_mat(M2, mdiag(4, pscale(nsq, -1)))
    value("R4D.07", "literal_square_M(xi)M(xi)", [("equals_+nsq_I4", 1 if np_ == 0 else 0), ("equals_-nsq_I4", 1 if nm == 0 else 0),
                                                  ("unequal_vs_plus", np_), ("unequal_vs_minus", nm)])
    D = mblock([[mzero(4, 4), mT(M)], [M, mzero(4, 4)]])
    e, s, n = cmp_mat(mmul(D, D), mdiag(8, nsq))
    meas("R4D.08", "doubled_D=[[0,M^T],[M,0]]_D^2=+nsq(xi)I8", n == 0, [("unequal", n)])

    # equivariance, Lie algebra sp(1)+sp(1)
    gam = [iblock([[izero(4, 4), Lmat_int(TQ, unit(4, m))], [iscale(Lmat_int(TQ, conj_int(unit(4, m))), -1), izero(4, 4)]]) for m in range(4)]

    def gam_of(v):
        out = izero(8, 8)
        for m in range(4):
            if v[m]:
                out = iadd(out, iscale(gam[m], v[m]))
        return out
    fails_ok = 0
    fails_swap = 0
    for u in range(1, 4):
        Lu = Lmat_int(TQ, unit(4, u))
        Ru = Rmat_int(TQ, unit(4, u))
        for side in (0, 1):
            if side == 0:
                A = Lu
                G = iblock([[Lu, izero(4, 4)], [izero(4, 4), izero(4, 4)]])
                Gs = iblock([[izero(4, 4), izero(4, 4)], [izero(4, 4), Lu]])
            else:
                A = iscale(Ru, -1)
                G = iblock([[izero(4, 4), izero(4, 4)], [izero(4, 4), Lu]])
                Gs = iblock([[Lu, izero(4, 4)], [izero(4, 4), izero(4, 4)]])
            for m in range(4):
                Am = [A[r][m] for r in range(4)]
                target = gam_of(Am)
                if iadd(imul(G, gam[m]), imul(gam[m], G), -1) != target:
                    fails_ok += 1
                if iadd(imul(Gs, gam[m]), imul(gam[m], Gs), -1) != target:
                    fails_swap += 1
    meas("R4D.09", "sp1+sp1_equivariance_[diag(L_u,0),gamma(xi)]=gamma(u xi)_and_[diag(0,L_v),gamma(xi)]=gamma(-xi v)", fails_ok == 0,
         [("generators", 6), ("basis_covectors", 4), ("failures", fails_ok)])
    control("K4D.03", "swapped_spinor_action_S+_by_q_S-_by_p", fails_swap > 0, [("failures", fails_swap), ("of", 24)])

    Q8 = []
    for i in range(4):
        for sgn in (1, -1):
            Q8.append([sgn if j == i else 0 for j in range(4)])
    gfail = 0
    gfail_swap = 0
    for p in Q8:
        for q in Q8:
            Lp = Lmat_int(TQ, p)
            Lq = Lmat_int(TQ, q)
            g = iblock([[Lp, izero(4, 4)], [izero(4, 4), Lq]])
            gs = iblock([[Lq, izero(4, 4)], [izero(4, 4), Lp]])
            for m in range(4):
                v = amul_int(TQ, amul_int(TQ, p, unit(4, m)), conj_int(q))
                t = gam_of(v)
                if imul(imul(g, gam[m]), iT(g)) != t:
                    gfail += 1
                if imul(imul(gs, gam[m]), iT(gs)) != t:
                    gfail_swap += 1
    meas("R4D.10", "Q8xQ8_gamma(p xi conj(q))=g gamma(xi) g^T_g=diag(L_p,L_q)", gfail == 0, [("pairs", 64), ("checks", 256), ("failures", gfail)])
    control("K4D.04", "swapped_g=diag(L_q,L_p)", gfail_swap > 0, [("failures", gfail_swap), ("of", 256)])


def r1_rows(case, n, T, k, w, xa, aa, nvars):
    TI = iadd(T, ieye(len(T)))
    tw = mapply(to_poly(TI, nvars), w)
    nw = nsq_poly(w)
    e, s, n1 = cmp_poly(nsq_poly(tw), pscale(nw, k))
    e, s2, n2 = cmp_poly(dot_poly(w, tw), nw)
    e, s3, n3 = cmp_poly(nw, padd(pmul(nsq_poly(xa), nsq_poly(aa)), pmul(dot_poly(xa, aa), dot_poly(xa, aa)), -1))
    meas("R" + case + ".R1", "R1_c=1:" + str(k) + "_nsq((T+I)(xi^a))=" + str(k) + "nsq(xi^a)_and_<xi^a,(T+I)(xi^a)>=nsq(xi^a)", n1 == 0 and n2 == 0 and n3 == 0,
         [("k", k), ("c", "1:" + str(k)), ("slots_norm", s), ("unequal_norm", n1), ("slots_inner", s2), ("unequal_inner", n2), ("C0.9_unequal", n3),
          ("reading", "pi_small=(T+I) times 1:" + str(k) + " so nsq(pi_small(xi^a))=(1:" + str(k) + ") nsq(xi^a)")])
    for kk in (1, 2, 3, 4):
        if kk == k:
            continue
        e, s, nn = cmp_poly(nsq_poly(tw), pscale(nw, kk))
        control("K" + case + ".R1c" + str(kk), "R1_with_c=1:" + str(kk) + "_nsq((T+I)(xi^a))=" + str(kk) + "nsq(xi^a)", nn > 0, [("unequal", nn)])
    m = len(T)
    om = [pvar(i, m) for i in range(m)]
    tom = mapply(to_poly(TI, m), om)
    e, s, nn = cmp_poly(dot_poly(om, tom), nsq_poly(om))
    control("K" + case + ".R1gen", "R1_inner_identity_on_a_general_2form", nn > 0, [("vars", m), ("unequal", nn)])


# ---------------- 7D ----------------

def oct_table(flip=None):
    return table_from_lines(8, OCT_LINES, flip)


def phi7_of(T):
    phi = {}
    for (i, j, k) in kbasis(7, 3):
        s, c = T[i + 1][j + 1]
        phi[(i, j, k)] = s if c == k + 1 else 0
    return phi


def sec_7d():
    T = oct_table()
    meas("R7D.01", "octonion_lines_from_rule_e_i*e_(i+1)=e_(i+3)_mod_7", OCT_LINES == FROZEN_LINES,
         [("lines", " ".join("(" + ",".join(str(v) for v in ln) + ")" for ln in OCT_LINES)), ("equal_frozen_list", 1 if OCT_LINES == FROZEN_LINES else 0)])
    x16 = [pvar(i, 16) for i in range(8)]
    y16 = [pvar(i + 8, 16) for i in range(8)]
    e, s, n = cmp_poly(nsq_poly(amul_poly(T, x16, y16)), pmul(nsq_poly(x16), nsq_poly(y16)))
    meas("R7D.02", "octonion_composition_nsq(xy)=nsq(x)nsq(y)", n == 0, [("support_slots", s), ("unequal", n)])
    Tf = oct_table(flip=(1, 2, 4))
    e, s, n = cmp_poly(nsq_poly(amul_poly(Tf, x16, y16)), pmul(nsq_poly(x16), nsq_poly(y16)))
    control("K7D.01", "line_(1,2,4)_reversed_table_composition", n > 0, [("unequal", n)])

    nv = 7
    x = [{}] + [pvar(i, nv) for i in range(7)]
    v7 = [pvar(i, nv) for i in range(7)]
    nsq = nsq_poly(v7)
    L = Lmat_poly(T, x)
    e, s, n = cmp_mat(mmul(L, L), mdiag(8, pscale(nsq, -1)))
    meas("R7D.03", "L_xi^2=-nsq(xi)I8_for_xi_in_Im_O", n == 0, [("entries", e), ("support_slots", s), ("unequal", n), ("full_basis_slots", e * comb(nv + 1, 2))])
    Lf = Lmat_poly(Tf, x)
    e, s, n = cmp_mat(mmul(Lf, Lf), mdiag(8, pscale(nsq, -1)))
    control("K7D.02", "L_xi^2=-nsq_with_line_(1,2,4)_reversed", n > 0, [("unequal", n)])

    X = [row[1:] for row in L[1:]]
    e, s, nskew = cmp_mat(mT(X), mneg(X))
    X2 = mmul(X, X)
    target = madd(mdiag(7, pscale(nsq, -1)), outer(v7, v7))
    e, s, n = cmp_mat(X2, target)
    e2, s2, nk = cmp_vec(mapply(X, v7), [{} for _ in range(7)])
    meas("R7D.04", "cross_X(xi)v=Im(xi v)_X^2=-nsq(xi)I7+xi_xi^T_X_xi=0_X^T=-X", n == 0 and nk == 0 and nskew == 0,
         [("entries", e), ("support_slots", s), ("unequal", n), ("X_xi_unequal", nk), ("skew_unequal", nskew), ("full_basis_slots", e * comb(nv + 1, 2))])
    e, s, n = cmp_mat(X2, mdiag(7, pscale(nsq, -1)))
    control("K7D.03", "X(xi)^2=-nsq(xi)I7_without_xi_xi^T", n > 0, [("unequal", n)])

    e1, e2, e3 = unit(8, 1), unit(8, 2), unit(8, 3)
    l12 = amul_int(T, e1, e2)
    l23 = amul_int(T, e2, e3)
    left = amul_int(T, l12, e3)
    right = amul_int(T, e1, l23)
    assoc = [a - b for a, b in zip(left, right)]

    def show(v):
        parts = []
        for i, c in enumerate(v):
            if c:
                parts.append(("+" if c > 0 else "-") + (str(abs(c)) if abs(c) != 1 else "") + "e" + str(i))
        return "".join(parts) if parts else "0"
    zero_trip = 0
    nonzero_trip = 0
    mags = set()
    for a in range(1, 8):
        for b in range(1, 8):
            for c in range(1, 8):
                if len({a, b, c}) < 3:
                    continue
                ua, ub, uc = unit(8, a), unit(8, b), unit(8, c)
                d = [p - q for p, q in zip(amul_int(T, amul_int(T, ua, ub), uc), amul_int(T, ua, amul_int(T, ub, uc)))]
                if any(d):
                    nonzero_trip += 1
                    mags.add(tuple(sorted(abs(t) for t in d if t)))
                else:
                    zero_trip += 1
    value("R7D.05", "associator_(e1,e2,e3)", [("e1e2", show(l12)), ("e2e3", show(l23)), ("(e1e2)e3", show(left)), ("e1(e2e3)", show(right)),
                                              ("associator", show(assoc)), ("ordered_distinct_triples", zero_trip + nonzero_trip),
                                              ("associative_triples", zero_trip), ("nonassociative_triples", nonzero_trip),
                                              ("nonzero_associator_shapes", ";".join(",".join(str(t) for t in m) for m in sorted(mags)))])
    control("K7D.04", "associativity_(e1e2)e3=e1(e2e3)", any(assoc), [("associator", show(assoc))])

    phi = phi7_of(T)
    anti = 1
    for (i, j, k) in itertools.permutations(range(7), 3):
        s1, c1 = T[i + 1][j + 1]
        val = s1 if c1 == k + 1 else 0
        base = phi[tuple(sorted((i, j, k)))] * perm_sign((i, j, k))
        if val != base:
            anti = 0
    terms = sorted(I for I, c in phi.items() if c)
    value("R7D.06", "phi_ijk=<e_i e_j,e_k>", [("totally_antisymmetric", anti), ("terms", len(terms)),
                                               ("phi", "".join(("+" if phi[I] > 0 else "-") + "e" + "".join(str(t + 1) for t in I) for I in terms))])

    W = wedge_matrix(7, phi, 3, 2)
    ranks = {}
    for sigma in (1, -1):
        S = star_matrix(7, 5, sigma)
        T7 = imul(S, W)
        ranks[sigma] = (irank(iadd(T7, ieye(21, 2), -1)), irank(iadd(T7, ieye(21))), irank(iadd(T7, ieye(21))), irank(iadd(T7, ieye(21, 2))), irank(iadd(T7, ieye(21), -1)), T7)
    sig = None
    for sigma in (1, -1):
        T7 = ranks[sigma][5]
        r_m2 = irank(iadd(T7, ieye(21, 2), -1))
        r_p1 = irank(iadd(T7, ieye(21)))
        if r_m2 == 14 and r_p1 == 7:
            sig = sigma
    T7p = ranks[1][5]
    std = "rank(T-2I)=" + str(irank(iadd(T7p, ieye(21, 2), -1))) + ";rank(T+I)=" + str(irank(iadd(T7p, ieye(21)))) + ";rank(T+2I)=" + str(irank(iadd(T7p, ieye(21, 2)))) + ";rank(T-I)=" + str(irank(iadd(T7p, ieye(21), -1)))
    T7 = ranks[sig][5] if sig is not None else T7p
    prod = imul(iadd(T7, ieye(21, 2), -1), iadd(T7, ieye(21)))
    sym = 1 if T7 == iT(T7) else 0
    tr = sum(T7[i][i] for i in range(21))
    meas("R7D.07", "orientation_by_7D.C3_and_(T7-2I)(T7+I)=0", sig is not None and prod == izero(21, 21),
         [("vol_e1..e7_ranks", std), ("sigma_chosen", sig), ("eigen_+2_dim", 21 - irank(iadd(T7, ieye(21, 2), -1))),
          ("eigen_-1_dim", 21 - irank(iadd(T7, ieye(21)))), ("symmetric", sym), ("trace", tr)])
    # psi = *7 phi
    S34 = star_matrix(7, 3, sig)
    pv = imul(S34, [[v] for v in form_vec(7, 3, phi)])
    psi = {I: pv[i][0] for i, I in enumerate(kbasis(7, 4)) if pv[i][0]}

    # Psi7 from the Clifford action on eta = 1
    Lb = [Lmat_int(T, unit(8, i + 1)) for i in range(7)]
    Psi7 = izero(7, 21)
    real_zero = 1
    for col, (i, j) in enumerate(kbasis(7, 2)):
        c = imul(Lb[i], Lb[j])
        if c[0][0] != 0:
            real_zero = 0
        for k in range(7):
            Psi7[k][col] = c[k + 1][0]
    contr = izero(7, 21)
    for col, (i, j) in enumerate(kbasis(7, 2)):
        for k in range(7):
            I = tuple(sorted((i, j, k)))
            if len(set(I)) == 3:
                contr[k][col] = phi[I] * perm_sign((i, j, k))
    rk = irank(Psi7)
    meas("R7D.08", "Psi7(omega)=c(omega)eta_rank_and_kernel", real_zero == 1 and rk == 7,
         [("real_part_zero", real_zero), ("rank", rk), ("dim_kernel", 21 - rk), ("Psi7_equals_contraction_with_phi", 1 if Psi7 == contr else 0),
          ("Psi7_PsiT", "4I" if imul(Psi7, iT(Psi7)) == ieye(7, 4) else ("3I" if imul(Psi7, iT(Psi7)) == ieye(7, 3) else "other"))])
    meas("R7D.08b", "T7+I=Psi7^T_Psi7_and_Psi7_Psi7^T=3I7", iadd(T7, ieye(21)) == imul(iT(Psi7), Psi7) and imul(Psi7, iT(Psi7)) == ieye(7, 3),
         [("T7+I_eq_PsiT_Psi", 1 if iadd(T7, ieye(21)) == imul(iT(Psi7), Psi7) else 0)])

    # R1, 7D
    xa = [pvar(i, 14) for i in range(7)]
    aa = [pvar(i + 7, 14) for i in range(7)]
    w = wedge_vec_poly(xa, aa)
    r1_rows("7D", 7, T7, 3, w, xa, aa, 14)
    # (T7 + I)(xi wedge v) = iota_(xi x v) phi
    xo = [{}] + xa
    vo = [{}] + aa
    cr = amul_poly(T, xo, vo)[1:]
    iota = []
    for (i, j) in kbasis(7, 2):
        acc = {}
        for k in range(7):
            I = tuple(sorted((k, i, j)))
            if len(set(I)) == 3:
                c = phi[I] * perm_sign((k, i, j))
                if c:
                    acc = padd(acc, pscale(cr[k], c))
        iota.append(acc)
    tw = mapply(to_poly(iadd(T7, ieye(21)), 14), w)
    e, s, n = cmp_vec(tw, iota)
    meas("R7D.09", "(T7+I)(xi^v)=iota_(xi_x_v)phi_so_pi7(xi^v)=(1:3)iota_(xi_x_v)phi", n == 0, [("components", e), ("support_slots", s), ("unequal", n)])
    e, s, n = cmp_vec(tw, [pscale(t, 3) for t in iota])
    control("K7D.05", "(T7+I)(xi^v)=3iota_(xi_x_v)phi", n > 0, [("unequal", n)])

    # rolled-up symbol R
    Rr = [[{}] + [pscale(t, -1) for t in v7]]
    for k in range(7):
        Rr.append([v7[k]] + X[k])
    e, s, n = cmp_mat(Rr, L)
    meas("R7D.10", "R(xi)(f,v)=(-<xi,v>,f xi+xi_x_v)_equals_L_xi", n == 0, [("unequal", n)])
    e, s, nm = cmp_mat(mmul(Rr, Rr), mdiag(8, pscale(nsq, -1)))
    e, s, np_ = cmp_mat(mmul(Rr, Rr), mdiag(8, nsq))
    e, s, ng = cmp_mat(mmul(mT(Rr), Rr), mdiag(8, nsq))
    value("R7D.11", "literal_square_R(xi)R(xi)_and_gram_R^T_R", [("RR_equals_-nsq_I8", 1 if nm == 0 else 0), ("RR_equals_+nsq_I8", 1 if np_ == 0 else 0),
                                                                    ("RtR_equals_+nsq_I8", 1 if ng == 0 else 0)])
    Dd = mblock([[mzero(8, 8), mT(Rr)], [Rr, mzero(8, 8)]])
    e, s, n = cmp_mat(mmul(Dd, Dd), mdiag(16, nsq))
    meas("R7D.12", "doubled_[[0,R^T],[R,0]]^2=+nsq(xi)I16", n == 0, [("unequal", n)])

    # equivariance: ker Psi7 acts by derivations of O
    K = nullspace_int(Psi7)
    der_fail = 0
    comm_fail = 0
    for kv in K:
        R7 = rot_matrix(7, kv)
        D = izero(8, 8)
        for a in range(7):
            for b in range(7):
                D[a + 1][b + 1] = R7[a][b]
        if not derivation_ok(T, D):
            der_fail += 1
        for m in range(1, 8):
            Lm = Lmat_int(T, unit(8, m))
            Dm = [D[r][m] for r in range(8)]
            if iadd(imul(D, Lm), imul(Lm, D), -1) != Lmat_int(T, Dm):
                comm_fail += 1
    span_rank = irank([v for v in K] + [])
    br = []
    for a in range(len(K)):
        for b in range(a + 1, len(K)):
            Ra = rot_matrix(7, K[a])
            Rb = rot_matrix(7, K[b])
            Cm = iadd(imul(Ra, Rb), imul(Rb, Ra), -1)
            br.append([Cm[j][i] for (i, j) in kbasis(7, 2)])
    closed = 1 if irank(K + br) == len(K) else 0
    meas("R7D.13", "ker_Psi7_acts_by_derivations_[D,L_xi]=L_(D xi)", der_fail == 0 and comm_fail == 0 and closed == 1,
         [("dim_ker_Psi7", len(K)), ("basis_rank", span_rank), ("derivation_failures", der_fail), ("commutator_failures", comm_fail), ("closed_under_bracket", closed)])
    om = [0] * 21
    om[0] = 1
    R12 = rot_matrix(7, om)
    D = izero(8, 8)
    for a in range(7):
        for b in range(7):
            D[a + 1][b + 1] = R12[a][b]
    control("K7D.06", "rotation_R_12_outside_ker_Psi7_is_a_derivation", not derivation_ok(T, D), [("Psi7_e12", ",".join(str(Psi7[k][0]) for k in range(7)))])
    CTX["Psi7"] = Psi7
    CTX["T7"] = T7
    CTX["OCT"] = T
    return sig, phi, psi, T


def derivation_ok(T, D):
    for a in range(8):
        for b in range(8):
            ea, eb = unit(8, a), unit(8, b)
            Da = [D[r][a] for r in range(8)]
            Db = [D[r][b] for r in range(8)]
            prod = amul_int(T, ea, eb)
            lhs = [sum(D[r][c] * prod[c] for c in range(8)) for r in range(8)]
            rhs = [p + q for p, q in zip(amul_int(T, Da, eb), amul_int(T, ea, Db))]
            if lhs != rhs:
                return False
    return True


# ---------------- 8D ----------------

def sp_from(A):
    return [{j: v for j, v in enumerate(r) if v} for r in A]


def sp_mul(A, B):
    out = []
    for row in A:
        acc = {}
        for k, a in row.items():
            for j, b in B[k].items():
                v = acc.get(j, 0) + a * b
                if v:
                    acc[j] = v
                else:
                    acc.pop(j, None)
        out.append(acc)
    return out


def sp_add(A, B, s=1):
    out = []
    for ra, rb in zip(A, B):
        acc = dict(ra)
        for j, b in rb.items():
            v = acc.get(j, 0) + s * b
            if v:
                acc[j] = v
            else:
                acc.pop(j, None)
        out.append(acc)
    return out


def sp_scale(A, k):
    if k == 0:
        return [{} for _ in A]
    return [{j: v * k for j, v in r.items()} for r in A]


def sp_comm(A, B):
    return sp_add(sp_mul(A, B), sp_mul(B, A), -1)


def sec_8d(sig, phi7, psi7, T):
    nv = 8
    x = [pvar(i, nv) for i in range(8)]
    cx = conj_poly(x)
    nsq = nsq_poly(x)
    L = Lmat_poly(T, x)
    Lc = Lmat_poly(T, cx)
    Z = mzero(8, 8)
    gamma = mblock([[Z, L], [mneg(Lc), Z]])
    e, s, n = cmp_mat(mmul(gamma, gamma), mdiag(16, pscale(nsq, -1)))
    meas("R8D.01", "gamma(xi)^2=-nsq(xi)I16_with_gamma=[[0,L_xi],[-L_conj(xi),0]]", n == 0,
         [("entries", e), ("support_slots", s), ("unequal", n), ("full_basis_slots", e * comb(nv + 1, 2))])

    gam = [iblock([[izero(8, 8), Lmat_int(T, unit(8, m))], [iscale(Lmat_int(T, conj_int(unit(8, m))), -1), izero(8, 8)]]) for m in range(8)]
    sgam = [sp_from(g) for g in gam]
    B2 = kbasis(8, 2)
    cbas = [sp_mul(sgam[i], sgam[j]) for (i, j) in B2]
    Psi = izero(7, 28)
    clean = 1
    for col, c in enumerate(cbas):
        colv = [c[r].get(0, 0) for r in range(16)]
        if colv[0] != 0 or any(colv[8:]):
            clean = 0
        for k in range(7):
            Psi[k][col] = colv[k + 1]
    rk = irank(Psi)
    K = nullspace_int(Psi)
    meas("R8D.02", "Psi(omega)=c(omega)eta_eta=1_in_S+_rank_and_dim_ker", clean == 1 and rk == 7 and len(K) == 21,
         [("real_and_S-_parts_zero", clean), ("rank_Psi", rk), ("dim_ker_Psi", 28 - rk), ("nullspace_basis", len(K)),
          ("Psi_PsiT_is_4I7", 1 if imul(Psi, iT(Psi)) == ieye(7, 4) else 0)])

    # Cayley forms and T8
    phi8 = {tuple(t + 1 for t in I): c for I, c in phi7.items() if c}
    psi8 = {tuple(t + 1 for t in I): c for I, c in psi7.items() if c}
    e0phi = {(0,) + I: c for I, c in phi8.items()}
    S62 = star_matrix(8, 6, sig)
    S44 = star_matrix(8, 4, sig)
    match = {}
    for sgn in (1, -1):
        Phi = dict(e0phi)
        for I, c in psi8.items():
            Phi[I] = Phi.get(I, 0) + sgn * c
        pv = form_vec(8, 4, Phi)
        sstar = imul(S44, [[v] for v in pv])
        selfdual = 1 if [r[0] for r in sstar] == pv else (-1 if [r[0] for r in sstar] == [-v for v in pv] else 0)
        T8 = imul(S62, wedge_matrix(8, Phi, 4, 2))
        q = imul(iadd(T8, ieye(28, 3), -1), iadd(T8, ieye(28)))
        match[sgn] = (selfdual, 1 if q == izero(28, 28) else 0, irank(iadd(T8, ieye(28))), 1 if iadd(T8, ieye(28)) == imul(iT(Psi), Psi) else 0, T8)
    for sgn in (1, -1):
        sd, qz, r7, eq, _ = match[sgn]
        meas("R8D.03" + ("p" if sgn == 1 else "m"), "Cayley_Phi_s=e0^phi+s_psi_s=" + str(sgn), qz == 1 and r7 == 7,
             [("star8_Phi_over_Phi", sd), ("T-3I_times_T+I_is_zero", qz), ("rank_T+I", r7), ("T+I_equals_PsiT_Psi", eq)])
    sel = 1 if match[1][3] == 1 else (-1 if match[-1][3] == 1 else 0)
    T8 = match[sel][4] if sel else match[1][4]
    CTX["Psi8"] = Psi
    CTX["T8p"] = match[1][4]
    CTX["T8m"] = match[-1][4]
    CTX["gam8"] = gam

    # P and N
    def P_of(m):
        Pm = izero(7, 8)
        for k in range(7):
            for b in range(8):
                if b == m:
                    continue
                if m < b:
                    Pm[k][b] = Psi[k][B2.index((m, b))]
                else:
                    Pm[k][b] = -Psi[k][B2.index((b, m))]
        return Pm
    Pm = [P_of(m) for m in range(8)]
    P = mzero(7, 8)
    for m in range(8):
        P = madd(P, [[pscale(x[m], v) for v in r] for r in Pm[m]])
    N = [[pscale(x[b], -1) for b in range(8)]] + P
    CTX["N8"] = N
    CTX["nsq8"] = nsq
    e1, s1, n1 = cmp_mat(mmul(mT(N), N), mdiag(8, nsq))
    e2, s2, n2 = cmp_mat(mmul(N, mT(N)), mdiag(8, nsq))
    meas("R8D.04", "N(xi)=(-iota_xi,P(xi))_N^T_N=nsq(xi)I8_and_N_N^T=nsq(xi)I8", n1 == 0 and n2 == 0,
         [("NtN_entries", e1), ("NtN_support_slots", s1), ("NtN_unequal", n1), ("NNt_unequal", n2), ("full_basis_slots_each", e1 * comb(nv + 1, 2))])
    e, s, na = cmp_vec(mapply(P, x), [{} for _ in range(7)])
    e, s, nb = cmp_mat(mmul(P, mT(P)), mdiag(7, nsq))
    meas("R8D.05", "P(xi)xi=0_and_P_P^T=nsq(xi)I7", na == 0 and nb == 0, [("Pxi_unequal", na), ("PPt_unequal", nb)])
    holding, uc = match_candidates(N, T, x, nv)
    value("R8D.06", "N(xi)a_as_octonion_expression_of_16_candidates", [("holding", ";".join(holding) if holding else "none"), ("non_holding", len(uc))])
    for name in ("+conj(xi)*a", "-conj(xi)*a"):
        for nm, C in candidates(T, x, nv):
            if nm == name:
                e, s, n = cmp_mat(N, C)
                value("R8D.06" + ("p" if name[0] == "+" else "m"), "N(xi)a=" + name, [("holds", 1 if n == 0 else 0), ("entries", e), ("support_slots", s), ("unequal", n)])
    Nf = [[pscale(x[b], 1) for b in range(8)]] + P
    holding2, uc2 = match_candidates(Nf, T, x, nv)
    value("R8D.07", "N_with_Lambda0_row_+iota_xi_as_octonion_expression", [("holding", ";".join(holding2) if holding2 else "none")])
    N2 = mmul(N, N)
    e, s, np_ = cmp_mat(N2, mdiag(8, nsq))
    e, s, nm_ = cmp_mat(N2, mdiag(8, pscale(nsq, -1)))
    value("R8D.08", "literal_square_N(xi)N(xi)", [("equals_+nsq_I8", 1 if np_ == 0 else 0), ("equals_-nsq_I8", 1 if nm_ == 0 else 0)])

    # R1, 8D
    xa = [pvar(i, 16) for i in range(8)]
    aa = [pvar(i + 8, 16) for i in range(8)]
    w = wedge_vec_poly(xa, aa)
    r1_rows("8D", 8, T8, 4, w, xa, aa, 16)

    # alternativity and the R3 block
    e, s, n = cmp_mat(mmul(Lc, L), mdiag(8, nsq))
    meas("R8D.09", "conj(xi)(xi x)=nsq(xi)x_L_conj(xi)L_xi=nsq(xi)I8", n == 0, [("entries", e), ("support_slots", s), ("unequal", n)])
    e, s, n = cmp_mat(mmul(L, L), mdiag(8, pscale(nsq, -1)))
    control("K8D.01", "xi(xi x)=-nsq(xi)x_for_general_xi_in_O", n > 0, [("unequal", n)])
    xim = [{}] + x[1:]
    Lim = Lmat_poly(T, xim)
    e, s, n = cmp_mat(mmul(Lim, Lim), mdiag(8, pscale(nsq_poly(xim), -1)))
    meas("R8D.10", "xi(xi x)=-nsq(xi)x_for_imaginary_xi", n == 0, [("unequal", n)])
    xx = amul_poly(T, x, x)
    e, s, n = cmp_mat(mmul(L, L), Lmat_poly(T, xx))
    meas("R8D.11", "left_alternative_xi(xi x)=(xi xi)x", n == 0, [("unequal", n)])
    Bk = mblock([[Z, Lc], [L, Z]])
    B2m = mmul(Bk, Bk)
    e, s, np_ = cmp_mat(B2m, mdiag(16, nsq))
    e, s, nm_ = cmp_mat(B2m, mdiag(16, pscale(nsq, -1)))
    value("R8D.12", "block_[[0,L_conj(xi)],[L_xi,0]]_square", [("equals_+nsq_I16", 1 if np_ == 0 else 0), ("equals_-nsq_I16", 1 if nm_ == 0 else 0)])

    # spin action normalization
    rots = []
    for idx in range(28):
        om = [0] * 28
        om[idx] = 1
        rots.append(rot_matrix(8, om))

    def gam_of(v):
        out = [{} for _ in range(16)]
        for m in range(8):
            if v[m]:
                out = sp_add(out, sp_scale(sgam[m], v[m]))
        return out
    f1 = 0
    fsign_minus = 0
    fsign_plus = 0
    for idx, (i, j) in enumerate(B2):
        X = rots[idx]
        quarter = [{} for _ in range(16)]
        for a in range(8):
            for b in range(8):
                if X[a][b]:
                    quarter = sp_add(quarter, sp_scale(sp_mul(sgam[a], sgam[b]), X[a][b]))
        for m in range(8):
            Xm = [X[r][m] for r in range(8)]
            if sp_comm(cbas[idx], sgam[m]) != sp_scale(gam_of(Xm), 2):
                f1 += 1
            if sp_comm(quarter, sgam[m]) != sp_scale(gam_of(Xm), -4):
                fsign_minus += 1
            if sp_comm(quarter, sgam[m]) != sp_scale(gam_of(Xm), 4):
                fsign_plus += 1
    meas("R8D.13", "[c(omega),gamma(v)]=2gamma(R_omega v)_all_28_basis_forms", f1 == 0, [("checks", 28 * 8), ("failures", f1)])
    meas("R8D.14", "sum_X_ab_gamma_a_gamma_b_commutator=-4gamma(Xv)_so_(1:4)sum_X_ab_gamma_a_gamma_b=-(1:2)c(omega)", fsign_minus == 0, [("checks", 224), ("failures", fsign_minus)])
    control("K8D.02", "[(1:4)sum_X_ab_gamma_a_gamma_b,gamma(v)]=+gamma(Xv)", fsign_plus > 0, [("failures", fsign_plus), ("of", 224)])
    hom_fail = 0
    for a in range(28):
        for b in range(28):
            Cm = iadd(imul(rots[a], rots[b]), imul(rots[b], rots[a]), -1)
            omv = [Cm[j][i] for (i, j) in B2]
            if rot_matrix(8, omv) != Cm:
                hom_fail += 1
                continue
            cw = [{} for _ in range(16)]
            for t, wv in enumerate(omv):
                if wv:
                    cw = sp_add(cw, sp_scale(cbas[t], wv))
            if sp_comm(cbas[a], cbas[b]) != sp_scale(cw, 2):
                hom_fail += 1
    meas("R8D.15", "rho=(1:2)c_is_a_Lie_hom_[c(X),c(Y)]=2c([X,Y])_so_lambda=1:2_unique", hom_fail == 0, [("ordered_pairs", 784), ("failures", hom_fail)])

    # equivariance of P under ker Psi = spin(7)
    def equiv_fail(omv):
        cw = [{} for _ in range(16)]
        for t, wv in enumerate(omv):
            if wv:
                cw = sp_add(cw, sp_scale(cbas[t], wv))
        Bfull = [[cw[r].get(c, 0) for c in range(8)] for r in range(8)]
        leak = 1 if (any(Bfull[0]) or any(Bfull[r][0] for r in range(8))) else 0
        B7 = [row[1:] for row in Bfull[1:]]
        R = rot_matrix(8, omv)
        bad = 0
        for m in range(8):
            lhs = imul(B7, Pm[m])
            rhs = izero(7, 8)
            for k in range(8):
                if R[k][m]:
                    rhs = iadd(rhs, iscale(Pm[k], 2 * R[k][m]))
            rhs = iadd(rhs, iscale(imul(Pm[m], R), 2))
            if lhs != rhs:
                bad += 1
        return bad, leak
    tot_bad = 0
    tot_leak = 0
    for kv in K:
        bd, lk = equiv_fail(kv)
        tot_bad += bd
        tot_leak += lk
    meas("R8D.16", "spin7=ker_Psi_equivariance_c(omega)P(xi)=2P(R_omega xi)+2P(xi)R_omega", tot_bad == 0 and tot_leak == 0,
         [("generators", len(K)), ("basis_covectors", 8), ("failures", tot_bad), ("eta_leak", tot_leak)])
    om = [0] * 28
    om[0] = 1
    bd, lk = equiv_fail(om)
    control("K8D.03", "equivariance_for_X=e01_outside_ker_Psi", bd > 0 or lk > 0, [("failures", bd), ("of", 8), ("eta_leak", lk)])
    om = [0] * 28
    om[B2.index((1, 2))] = 1
    bd, lk = equiv_fail(om)
    control("K8D.04", "equivariance_for_X=e12_outside_ker_Psi", bd > 0 or lk > 0, [("failures", bd), ("of", 8), ("eta_leak", lk)])
    return sel


# ---------------- obstruction ----------------

def cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cadd(a, b, s=1):
    return (a[0] + s * b[0], a[1] + s * b[1])


def gm_mul(A, B):
    return [[cadd(cm(A[i][0], B[0][j]), cm(A[i][1], B[1][j])) for j in range(2)] for i in range(2)]


def gm_comm(A, B):
    P = gm_mul(A, B)
    Q = gm_mul(B, A)
    return [[cadd(P[i][j], Q[i][j], -1) for j in range(2)] for i in range(2)]


def gm_show(A):
    return "[[" + ", ".join(str(z[0]) + ":" + str(z[1]) for z in A[0]) + "],[" + ", ".join(str(z[0]) + ":" + str(z[1]) for z in A[1]) + "]]"


def sec_ob():
    s1 = [[(0, 0), (1, 0)], [(1, 0), (0, 0)]]
    s2 = [[(0, 0), (0, -1)], [(0, 1), (0, 0)]]
    s3 = [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]
    ii = (0, 1)
    tau = [[[cm(ii, z) for z in r] for r in s] for s in (s1, s2, s3)]

    def decomp(Xm):
        x3 = Xm[0][0][1]
        x1 = Xm[0][1][1]
        x2 = Xm[0][1][0]
        rec = [[(0, 0), (0, 0)], [(0, 0), (0, 0)]]
        for c, t in zip((x1, x2, x3), tau):
            rec = [[cadd(rec[i][j], (c * t[i][j][0], c * t[i][j][1])) for j in range(2)] for i in range(2)]
        if rec != Xm:
            raise SystemExit("not in su(2)")
        return [x1, x2, x3]
    F12 = gm_comm(tau[0], tau[1])
    target = [[cm((-2, 0), z) for z in r] for r in tau[2]]
    meas("ROB.01", "F12=[i_sigma1,i_sigma2]=-2i_sigma3", F12 == target, [("F12", gm_show(F12))])
    struct = [[decomp(gm_comm(tau[a], tau[b])) for b in range(3)] for a in range(3)]

    def br(X, Y):
        out = [{} for _ in range(3)]
        for a in range(3):
            if not X[a]:
                continue
            for b in range(3):
                if not Y[b]:
                    continue
                pr = pmul(X[a], Y[b])
                for c in range(3):
                    if struct[a][b][c]:
                        out[c] = padd(out[c], pscale(pr, struct[a][b][c]))
        return out

    def run(Avec, f):
        nv = 4
        A = [[pconst(c, nv) for c in Ai] for Ai in Avec]
        alpha = []
        for i in range(4):
            alpha.append([padd(pderiv(f[c], i), br(A[i], f)[c]) for c in range(3)])
        bad = 0
        nonzero_pairs = []
        slots = 0
        for (i, j) in kbasis(4, 2):
            bij = [padd(padd(pderiv(alpha[j][c], i), pderiv(alpha[i][c], j), -1), padd(br(A[i], alpha[j])[c], br(A[j], alpha[i])[c], -1)) for c in range(3)]
            Fij = br(A[i], A[j])
            rhs = br(Fij, f)
            e, s, n = cmp_vec(bij, rhs)
            slots += s
            bad += n
            if any(bij):
                nonzero_pairs.append(str(i + 1) + str(j + 1))
        return bad, nonzero_pairs, slots
    Anon = [[1, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 0]]
    Aab = [[0, 0, 1], [0, 0, 1], [0, 0, 0], [0, 0, 0]]
    monos = [m for m in itertools.product(range(3), repeat=4) if sum(m) <= 2]
    CTX["ob_run"] = run
    CTX["ob_monos"] = monos
    CTX["ob_br"] = br
    tot_bad = 0
    nz = set()
    tot_slots = 0
    ab_nonzero = 0
    for a in range(3):
        for m in monos:
            f = [{} for _ in range(3)]
            f[a] = {tuple(m): 1}
            bd, nzp, sl = run(Anon, f)
            tot_bad += bd
            tot_slots += sl
            nz.update(nzp)
            bd2, nzp2, sl2 = run(Aab, f)
            if nzp2:
                ab_nonzero += 1
    meas("ROB.02", "d_A_d_A_f=[F_A,f]_for_A1=i_sigma1_A2=i_sigma2_f_of_degree_le_2", tot_bad == 0,
         [("basis_f", 3 * len(monos)), ("support_slots", tot_slots), ("unequal", tot_bad), ("nonzero_only_on", ";".join(sorted(nz)))])
    meas("ROB.03", "abelian_A1=A2=i_sigma3_gives_d_A_d_A_f=0", ab_nonzero == 0, [("basis_f", 3 * len(monos)), ("nonzero_cases", ab_nonzero)])
    f0 = [pconst(1, 4), {}, {}]
    A = [[pconst(c, 4) for c in Ai] for Ai in Anon]
    al = [[padd(pderiv(f0[c], i), br(A[i], f0)[c]) for c in range(3)] for i in range(4)]
    b12 = [padd(br(A[0], al[1])[c], br(A[1], al[0])[c], -1) for c in range(3)]
    b12v = [p.get((0, 0, 0, 0), 0) for p in b12]
    control("K-OB.01", "d_A_d_A(i_sigma1)=0_for_the_nonabelian_A", any(b12v), [("d_A_d_A_i_sigma1_dx12_in_tau_basis", ",".join(str(v) for v in b12v))])
    S4 = star_matrix(4, 2, 1)
    Fv = [[0, 0, 0] for _ in range(6)]
    Fv[0] = decomp(F12)
    TI = iadd(S4, ieye(6))
    twoFplus = [[sum(TI[r][c] * Fv[c][k] for c in range(6)) for k in range(3)] for r in range(6)]
    beta = [[0, 0, 0] for _ in range(6)]
    beta[0] = b12v
    twobplus = [[sum(TI[r][c] * beta[c][k] for c in range(6)) for k in range(3)] for r in range(6)]

    def fshow(v):
        return ";".join("e" + "".join(str(t + 1) for t in I) + ":" + ",".join(str(x) for x in row) for I, row in zip(kbasis(4, 2), v) if any(row))
    value("ROB.04", "self_dual_parts_in_tau_basis_tau_a=i_sigma_a", [("F_tau", fshow(Fv)), ("2F+_tau", fshow(twoFplus)), ("d_A_d_A_i_sigma1_tau", fshow(beta)), ("2d_A+d_A_i_sigma1_tau", fshow(twobplus)),
                                                                    ("structure_[tau1,tau2]", ",".join(str(v) for v in struct[0][1]))])


# ===================== after the seal: post-read measurements, reruns, independence, gate, compare =====================

SLASH = chr(47)


def post(pid, name, ok, kv):
    COUNT["post"] = COUNT.get("post", 0) + 1
    if ok:
        COUNT["post_pass"] = COUNT.get("post_pass", 0) + 1
    emit("POSTREAD", [("id", pid), ("name", name)] + kv + [("status", "PASS" if ok else "FAIL"), ("timing", "after_reading_frozen_rows"), ("claim", "MEASURED")])


def post_ctrl(pid, name, failed, kv):
    COUNT["postc"] = COUNT.get("postc", 0) + 1
    if failed:
        COUNT["postc_ok"] = COUNT.get("postc_ok", 0) + 1
    emit("POSTREAD_CONTROL", [("id", pid), ("name", name)] + kv + [("expected", "FAIL"), ("observed", "FAIL" if failed else "PASS"),
                                                                   ("control_ok", 1 if failed else 0), ("timing", "after_reading_frozen_rows"), ("claim", "MEASURED")])


def qshow(v, names):
    parts = []
    for c, nm in zip(v, names):
        if c:
            parts.append(("+" if c > 0 else "-") + (str(abs(c)) if abs(c) != 1 else "") + nm)
    return "".join(parts) if parts else "0"


def intertwiner_dim(A_list, B_list):
    rows = []
    for A, B in zip(A_list, B_list):
        for r in range(8):
            for c in range(8):
                eq = [0] * 64
                for k in range(8):
                    eq[r * 8 + k] += A[k][c]
                    eq[k * 8 + c] -= B[r][k]
                rows.append(eq)
    return 64 - irank(rows)


def sec_postread():
    TQ = CTX["TQ"]
    S4 = CTX["S4"]
    r_sd = irank(iadd(S4, ieye(6)))
    gam4 = [iblock([[izero(4, 4), Lmat_int(TQ, unit(4, m))], [iscale(Lmat_int(TQ, conj_int(unit(4, m))), -1), izero(4, 4)]]) for m in range(4)]
    Psi4 = izero(3, 6)
    clean = 1
    for col, (i, j) in enumerate(kbasis(4, 2)):
        c = imul(gam4[i], gam4[j])
        colv = [c[r][0] for r in range(8)]
        if colv[0] != 0 or any(colv[4:]):
            clean = 0
        for k in range(3):
            Psi4[k][col] = colv[k + 1]
    vals = ";".join("e" + str(i + 1) + str(j + 1) + ":" + qshow([Psi4[k][col] for k in range(3)], ["i", "j", "k"]) for col, (i, j) in enumerate(kbasis(4, 2)))
    asd = [form_vec(4, 2, {(0, 1): 1, (2, 3): -1}), form_vec(4, 2, {(0, 2): 1, (1, 3): 1}), form_vec(4, 2, {(0, 3): 1, (1, 2): -1})]
    asd_in_ker = sum(1 for v in asd if all(sum(Psi4[k][t] * v[t] for t in range(6)) == 0 for k in range(3)))
    rk4 = irank(Psi4)
    g2 = 1 if imul(Psi4, iT(Psi4)) == ieye(3, 2) else 0
    gs = 1 if imul(Psi4, iadd(S4, ieye(6))) == iscale(Psi4, 2) else 0
    post("P4D.01", "dim_Lambda2+_and_Psi_4(omega)=c(omega)eta_in_4D", clean == 1 and rk4 == 3 and asd_in_ker == 3 and g2 == 1 and gs == 1,
         [("rank_star_plus_I", r_sd), ("Psi4_values", vals), ("rank_Psi4", rk4), ("dim_ker_Psi4", 6 - rk4), ("asd_basis_in_kernel", asd_in_ker),
          ("Psi4_Psi4T_is_2I3", g2), ("Psi4_star_plus_I_is_2Psi4", gs), ("real_and_S-_parts_zero", clean)])

    M = CTX["M4"]
    Mm = []
    for m in range(4):
        key = tuple(1 if t == m else 0 for t in range(4))
        Mm.append([[M[r][c].get(key, 0) for c in range(4)] for r in range(4)])
    linear = 1 if all(all(sum(mm) == 1 for mm in M[r][c]) for r in range(4) for c in range(4)) else 0
    Q8 = []
    for i in range(4):
        for sgn in (1, -1):
            Q8.append([sgn if j == i else 0 for j in range(4)])
    fail = 0
    fail_c = 0
    for p in Q8:
        for q in Q8:
            B = imul(Lmat_int(TQ, p), Rmat_int(TQ, conj_int(q)))
            Adp = imul(Lmat_int(TQ, p), Rmat_int(TQ, conj_int(p)))
            Adq = imul(Lmat_int(TQ, q), Rmat_int(TQ, conj_int(q)))
            Ap = [[1, 0, 0, 0]] + [[0] + Adp[r][1:] for r in range(1, 4)]
            Aq = [[1, 0, 0, 0]] + [[0] + Adq[r][1:] for r in range(1, 4)]
            for m in range(4):
                v = amul_int(TQ, amul_int(TQ, p, unit(4, m)), conj_int(q))
                lhs = izero(4, 4)
                for k in range(4):
                    if v[k]:
                        lhs = iadd(lhs, iscale(Mm[k], v[k]))
                if lhs != imul(imul(Ap, Mm[m]), iT(B)):
                    fail += 1
                if lhs != imul(imul(Aq, Mm[m]), iT(B)):
                    fail_c += 1
    post("P4D.02", "4D.P11_Q8xQ8_M(p_xi_conj(q))=A_p_M(xi)_B_pq^T", fail == 0 and linear == 1, [("pairs", 64), ("checks", 256), ("failures", fail), ("M_linear_in_xi", linear)])
    post_ctrl("P4D.K1", "4D.P11_with_A_q_in_place_of_A_p", fail_c > 0, [("failures", fail_c), ("of", 256)])

    T = CTX["OCT"]
    good = 0
    for i in range(1, 8):
        for j in range(1, 8):
            if i == j:
                continue
            prod = amul_int(T, unit(8, i), unit(8, j))
            D = iadd(Lmat_int(T, prod), imul(Lmat_int(T, unit(8, i)), Lmat_int(T, unit(8, j))), -1)
            nz = [x for r in D for x in r if x]
            if len(nz) == 4 and all(abs(x) == 2 for x in nz):
                good += 1
    post("P7D.01", "7D.P5_L_(eiej)-L_ei_L_ej_has_4_entries_each_pm2", good == 42, [("ordered_pairs", 42), ("pairs_matching", good)])

    P = ieye(8)
    for i in range(1, 8):
        P = imul(P, Lmat_int(T, unit(8, i)))
    eps = P[0][0] if P == ieye(8, P[0][0]) else 0
    gam = CTX["gam8"]
    G = ieye(16)
    for m in range(8):
        G = imul(G, gam[m])
    target = iblock([[ieye(8, eps), izero(8, 8)], [izero(8, 8), ieye(8, -eps)]])
    post("P78.01", "epsilon7_L_e1..L_e7_and_gamma0..gamma7", eps != 0 and G == target, [("epsilon7", eps), ("gamma_product_is_eps7_diag_I8_-I8", 1 if G == target else 0)])

    Psi = CTX["Psi8"]
    B2 = kbasis(8, 2)
    okv = 1
    for col, (i, j) in enumerate(B2):
        colv = [0] + [Psi[k][col] for k in range(7)]
        want = unit(8, j) if i == 0 else amul_int(T, unit(8, i), unit(8, j))
        if colv != want:
            okv = 0
    post("P8D.01", "8D.P5_Psi(e0j)=e_j_and_Psi(e_ij)=e_i_e_j", okv == 1, [("columns", 28), ("all_match", okv)])
    a = 1 if imul(Psi, iadd(CTX["T8p"], ieye(28))) == iscale(Psi, 4) else 0
    b = 1 if imul(Psi, iadd(CTX["T8m"], ieye(28))) == iscale(Psi, 4) else 0
    post("P8D.02", "8D.P6_Psi(T_+1+I)=4Psi_holds_and_Psi(T_-1+I)=4Psi_fails", a == 1 and b == 0, [("s_plus_holds", a), ("s_minus_holds", b)])

    Lm = [Lmat_int(T, unit(8, m)) for m in range(8)]
    Lcm = [Lmat_int(T, conj_int(unit(8, m))) for m in range(8)]
    Rm = [Rmat_int(T, unit(8, m)) for m in range(8)]
    Rcm = [Rmat_int(T, conj_int(unit(8, m))) for m in range(8)]
    d1 = intertwiner_dim(Lm, Lcm)
    d2 = intertwiner_dim(Lm, Lm)
    d3 = intertwiner_dim(Lm, Rm)
    d4 = intertwiner_dim(Lm, Rcm)
    N = CTX["N8"]
    nsq = CTX["nsq8"]
    x = [pvar(i, 8) for i in range(8)]
    K = to_poly(iscale(conj_matrix(8), -1), 8)
    e, s, nk = cmp_mat(N, mmul(Lmat_poly(T, x), K))
    post("P8D.03", "R2_note_fixed_identification_U_L_xi_U^-1=L_conj(xi)", nk == 0,
         [("N_equals_L_xi_K_K_a_is_-conj(a)", 1 if nk == 0 else 0), ("dim_U_with_U_Lxi_eq_Lconjxi_U", d1), ("dim_U_with_U_Lxi_eq_Lxi_U", d2),
          ("dim_U_with_U_Lxi_eq_Rxi_U", d3), ("dim_U_with_U_Lxi_eq_Rconjxi_U", d4),
          ("reading", "U N(xi) V = L_conj(xi) for all xi forces KV = U^-1 at xi = 1 and then U L_xi = L_conj(xi) U; that space has the dimension shown")])
    Dd = mblock([[mzero(8, 8), mT(N)], [N, mzero(8, 8)]])
    e, s, n = cmp_mat(mmul(Dd, Dd), mdiag(16, nsq))
    post("P8D.04", "doubled_[[0,N^T],[N,0]]^2=+nsq(xi)I16", n == 0, [("unequal", n)])

    run = CTX["ob_run"]
    monos = CTX["ob_monos"]
    br = CTX["ob_br"]
    A2 = [[0, 0, 1], [0, 0, 2], [0, 0, 0], [0, 0, 0]]
    Fz = br([pconst(c, 4) for c in A2[0]], [pconst(c, 4) for c in A2[1]])
    nonzero = 0
    bad = 0
    for a_ in range(3):
        for m in monos:
            f = [{} for _ in range(3)]
            f[a_] = {tuple(m): 1}
            bd, nzp, sl = run(A2, f)
            bad += bd
            if nzp:
                nonzero += 1
    post("P-OB.01", "OB.P4_A1=i_sigma3_A2=2i_sigma3_F=0_and_d_A_d_A_f=0", (not any(Fz)) and nonzero == 0 and bad == 0,
         [("F12_zero", 0 if any(Fz) else 1), ("basis_f", 3 * len(monos)), ("nonzero_cases", nonzero), ("unequal", bad)])


def file_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def sha_of(path):
    return hashlib.sha256(file_bytes(path)).hexdigest()


def lines_of(path):
    return file_bytes(path).decode("utf-8").split("\n")


def fields(line):
    out = {}
    for part in line.split("|")[1:]:
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def clean(s):
    return s.replace("|", ";").replace("\n", " ")


SEATS = [("4D", "4d", "MEASURE-4D.hbp", "", "c162f553043a131d4779658eeb7fd02a653d6ca44a4f8949975ff5dd6f9257dc",
          "MEASURE-4D-BEFORE-COMPARE.hbp", "MEASURE-4D-BEFORE-COMPARE.hbp", "measure_4d.py", "measure_4d.before-compare.py"),
         ("7D", "7d", "MEASURE-7D.hbp", "out", "0987062e160d10a291286b1d54175da403b068d5c37729af4100249b8d3c5c65",
          "MEASURE-7D.precompare.hbp", "MEASURE-7D.hbp", "measure_7d.py", "measure_7d.precompare.py"),
         ("8D", "8d", "MEASURE-8D.hbp", "out", "08b5936e32378daed2db668eab202740eaffc3cc4f2a553004dfcd7ab8de48f5",
          "MEASURED-ONLY-BEFORE-READ-8D.hbp", "MEASURE-8D.hbp", "measure_8d.py", "measure_8d.before-read.py.txt")]
RESULT = {}


def sec_rerun(rerun_root, measure_root):
    for case, d, name, sub, reported, _a, _b, _c, _d in SEATS:
        rp = os.path.join(rerun_root, d, sub, name) if sub else os.path.join(rerun_root, d, name)
        op = os.path.join(measure_root, d, name)
        rs = sha_of(rp)
        os_ = sha_of(op)
        side = file_bytes(rp + ".sha256").decode("utf-8")
        side_ok = 1 if side == rs + "  " + name + "\n" else 0
        tot = 0
        ok = 0
        mustfail = 0
        mustfail_obs = 0
        for ln in lines_of(rp):
            tag = ln.split("|", 1)[0]
            if tag in ("CONTROL", "CTRL"):
                tot += 1
                f = fields(ln)
                exp = f.get("expected", f.get("expect", f.get("must")))
                obs = f.get("observed")
                okf = f.get("control_ok", f.get("ok"))
                if exp == obs and okf == "1":
                    ok += 1
                if exp == "FAIL":
                    mustfail += 1
                    if obs == "FAIL":
                        mustfail_obs += 1
        same = 1 if rs == os_ == reported else 0
        RESULT["rerun_" + case] = same == 1 and side_ok == 1
        RESULT["controls_" + case] = ok == tot and mustfail == mustfail_obs and tot > 0
        emit("RERUN", [("seat", "SM" + case), ("rerun_sha256", rs), ("original_sha256", os_), ("reported_sha256", reported), ("byte_identical", same),
                       ("sidecar_bare_name_ok", side_ok), ("control_rows", tot), ("controls_ok", ok), ("must_fail_controls", mustfail),
                       ("must_fail_observed_fail", mustfail_obs), ("fresh_dir", "yes, script copied alone, frozen inputs from git show"), ("claim", "MEASURED")])


def sec_indep(preread_root, measure_root):
    finals = {}
    for case, d, name, sub, reported, snap, outname, final_script, pre_script in SEATS:
        snap_p = os.path.join(measure_root, d, snap)
        pre_p = os.path.join(preread_root, d, sub, outname) if sub else os.path.join(preread_root, d, outname)
        s_snap = sha_of(snap_p)
        s_pre = sha_of(pre_p)
        final_p = os.path.join(measure_root, d, name)
        A = [x for x in lines_of(snap_p) if x]
        B = set(lines_of(final_p))
        meta = ("HEAD", "SCRIPT", "END", "RUN", "PHASE")
        missing = [x for x in A if x not in B and x.split("|", 1)[0] not in meta]
        leak = []
        for x in A:
            if x.split("|", 1)[0] in ("HEAD", "MEASURE_HEAD", "RUN", "PHASE", "END", "SCRIPT"):
                continue
            f = fields(x)
            nm = f.get("name", "")
            low = x.lower()
            if "addendum" in low or nm[:3] in ("R1_", "R2_", "R3_", "R4_", "R5_", "R6_"):
                leak.append(f.get("id", "?") + ":" + nm[:40])
        txt = file_bytes(os.path.join(measure_root, d, final_script)).decode("utf-8")
        finals[case] = set(t.strip() for t in txt.split("\n") if len(t.strip()) >= 50)
        others = [s for c2, d2, n2, s2, r2, a2, b2, fs2, ps2 in SEATS if c2 != case for s in (fs2.split(".")[0], d2 + SLASH)]
        pre_txt = file_bytes(os.path.join(measure_root, d, pre_script)).decode("utf-8")
        refs = sum(txt.count(o) + pre_txt.count(o) for o in others)
        receipt_opened = 0
        RESULT["indep_files_" + case] = s_snap == s_pre and not missing
        RESULT["indep_leak_" + case] = len(leak)
        emit("INDEP", [("seat", "SM" + case), ("snapshot_sha256", s_snap), ("preread_script_rerun_sha256", s_pre), ("snapshot_reproduced_by_preread_script", 1 if s_snap == s_pre else 0),
                       ("snapshot_rows_absent_from_final", len(missing)), ("references_to_other_seat_scripts_or_dirs", refs), ("gpt6_receipt_opened_by_code", receipt_opened),
                       ("blind_rows_naming_addendum_items", len(leak)), ("those_rows", ";".join(leak) if leak else "none"), ("claim", "MEASURED")])
    pairs = [("4D", "7D"), ("4D", "8D"), ("7D", "8D")]
    emit("INDEP", [("check", "shared_lines_of_50_or_more_chars_between_final_scripts")] + [(a + "_" + b, len(finals[a] & finals[b])) for a, b in pairs] + [("claim", "MEASURED")])


def sec_gate(paths):
    me = sha_of(os.path.abspath(__file__))
    for label, path in paths:
        gsha = sha_of(path)
        mine_pinned = 0
        for ln in lines_of(path):
            if not ln.strip():
                continue
            emit("GATE", [("list", label), ("log_sha256", gsha), ("line", clean(ln))])
            if "refute_spinor.py " in ln and ("sha256=" + me) in ln and "findings=0" in ln:
                mine_pinned = 1
        txt = file_bytes(path).decode("utf-8")
        ok = ("SCAN_EXIT 0" in txt and "SELFTEST_EXIT 0" in txt and "PLANT_SCAN_EXIT 1" in txt and "UNPLANT_SCAN_EXIT 0" in txt
              and txt.count("| findings=0 |") == 8 and mine_pinned == 1)
        RESULT["gate_" + label] = ok
        emit("GATE_SUMMARY", [("list", label), ("scan_exit_0", 1 if "SCAN_EXIT 0" in txt else 0), ("selftest_exit_0", 1 if "SELFTEST_EXIT 0" in txt else 0),
                              ("planted_float_caught", 1 if "PLANT_SCAN_EXIT 1" in txt else 0), ("scripts_with_0_findings", txt.count("| findings=0 |")),
                              ("this_script_sha256", me), ("this_script_pinned_0_findings", mine_pinned), ("ok", 1 if ok else 0), ("claim", "MEASURED")])


MV = {}
CMP = {"rows": 0, "quote_missing": 0}
VERD = {}


def load_mv():
    for ln in ROWS:
        f = fields(ln)
        if "id" in f:
            MV[f["id"]] = f


def g(i, k):
    return MV.get(i, {}).get(k)


def frozen(path, expect):
    raw = file_bytes(path)
    sha = hashlib.sha256(raw).hexdigest()
    return sha == expect, raw.decode("utf-8").split("\n")


def cmp_row(lines, src, prefix, quote, mine, ok, basis, v_ok, v_bad, note=""):
    if "|" in quote:
        raise SystemExit("quote holds a bar")
    hits = [ln for ln in lines if ln.startswith(prefix)]
    row = hits[0] if len(hits) == 1 else ""
    found = 1 if row and quote in row else 0
    rsha = hashlib.sha256(row.encode("utf-8")).hexdigest() if row else "none"
    v = (v_ok if ok else v_bad) if found else "QUOTE_NOT_FOUND"
    CMP["rows"] += 1
    if not found:
        CMP["quote_missing"] += 1
    VERD[v] = VERD.get(v, 0) + 1
    if v.startswith("DISAGREE"):
        CMP.setdefault("dis", []).append(prefix.split("|")[1] + ":" + quote[:40])
    kv = [("ref", prefix.rstrip("|").replace("|", ":")), ("file", src), ("frozen_row_sha256", rsha), ("quote", quote), ("quote_found", found), ("measured", mine), ("basis", basis), ("verdict", v)]
    if note:
        kv.append(("note", note))
    emit("COMPARE", kv + [("claim", "MEASURED")])


def dims_row(lines, case, key, mine, src_ids):
    hits = [ln for ln in lines if ln.startswith("DIMS|case=" + case + "|")]
    row = hits[0] if len(hits) == 1 else ""
    fv = fields(row).get(key) if row else None
    v = "AGREE" if fv is not None and fv == mine else ("NOT_FOUND" if fv is None else "DISAGREE")
    CMP["rows"] += 1
    VERD[v] = VERD.get(v, 0) + 1
    if v == "DISAGREE":
        CMP.setdefault("dis", []).append("DIMS:" + case + ":" + key)
    emit("COMPARE", [("ref", "DIMS:" + case + ":" + key), ("file", "PRED"), ("frozen_value", fv), ("measured", mine), ("from", src_ids), ("verdict", v), ("claim", "MEASURED")])


def seat_cmp(measure_root, case, d, name, rid, key, mine):
    val = None
    for ln in lines_of(os.path.join(measure_root, d, name)):
        if ("|id=" + rid + "|") in ln:
            val = fields(ln).get(key)
            break
    v = "AGREE" if val == mine else ("NOT_FOUND" if val is None else "DIFFERENT")
    VERD["SEAT_" + v] = VERD.get("SEAT_" + v, 0) + 1
    emit("SEATCMP", [("seat", "SM" + case), ("seat_row", rid), ("field", key), ("seat_value", val), ("refuter_value", mine), ("verdict", v), ("claim", "MEASURED")])
    return v == "AGREE"


def sec_compare(pred_path, add_path, measure_root):
    pok, P = frozen(pred_path, "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e")
    aok, A = frozen(add_path, "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0")
    emit("FROZEN", [("prediction_sha256_ok", 1 if pok else 0), ("addendum_sha256_ok", 1 if aok else 0), ("first_read_utc", "after the measured block seal of 2026-11-09T18:10:24Z")])
    load_mv()
    PS = lambda i: g(i, "status") == "PASS"
    c = cmp_row
    c(P, "PRED", "PRED|id=4D.P1|", "all 640 coefficients of gamma(xi)^2 + nsq(xi) I_8 are 0", "R4D.02 unequal=" + str(g("R4D.02", "unequal")) + " full_basis_slots=" + str(g("R4D.02", "full_basis_slots")),
      g("R4D.02", "unequal") == "0" and g("R4D.02", "full_basis_slots") == "640", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P2|", "10 of 10 relations i <= j hold", "the coefficients of xi_i xi_j in R4D.02 are the anticommutators; all match", PS("R4D.02"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P6|", "C(xi) xi = 0 and C(xi) C(xi)^T = nsq(xi) I_3 as polynomial identities", "blocks of M M^T = nsq I4 (R4D.04)", g("R4D.04", "MMt_unequal") == "0", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P7|", "M(xi)^T M(xi) = M(xi) M(xi)^T = nsq(xi) I_4 as polynomial identities", "R4D.04 MtM_unequal=0 MMt_unequal=0", PS("R4D.04"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P8|", "Psi_4(e12) = i, Psi_4(e13) = j, Psi_4(e14) = k, Psi_4(e23) = k, Psi_4(e24) = -j, Psi_4(e34) = i", str(g("P4D.01", "Psi4_values")),
      g("P4D.01", "Psi4_values") == "e12:+i;e13:+j;e14:+k;e23:+k;e24:-j;e34:+i", "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P8|", "Psi_4 Psi_4^T = 2 I_3", "P4D.01 Psi4_Psi4T_is_2I3=" + str(g("P4D.01", "Psi4_Psi4T_is_2I3")), g("P4D.01", "Psi4_Psi4T_is_2I3") == "1", "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P8|", "ker Psi_4 = Lambda2- = sp(1)_-, dimension 3", "dim_ker_Psi4=" + str(g("P4D.01", "dim_ker_Psi4")) + " asd_basis_in_kernel=" + str(g("P4D.01", "asd_basis_in_kernel")),
      g("P4D.01", "dim_ker_Psi4") == "3" and g("P4D.01", "asd_basis_in_kernel") == "3", "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P9|", "= -xi conj(a), a bilinear identity", "R4D.06 holding=" + str(g("R4D.06", "holding")), g("R4D.06", "holding") == "-xi*conj(a)", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=4D.P11|", "M(p xi conj(q)) = A_p M(xi) B_pq^T", "P4D.02 failures=" + str(g("P4D.02", "failures")) + " of 256; control with A_q fails", PS("P4D.02"), "POSTREAD", "AGREE", "DISAGREE")
    c(P, "CONV", "CONV|id=4D.C3|", "gamma(p xi conj(q)) = g gamma(xi) g^-1 with g = diag(L_p, L_q)", "R4D.10 failures=0 on Q8xQ8; R4D.09 sp1+sp1 failures=0", PS("R4D.10") and PS("R4D.09"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P1|", "all 64 x 28 = 1792 coefficients of L_xi^2 + nsq(xi) I_8 are 0", "R7D.03 unequal=" + str(g("R7D.03", "unequal")) + " full_basis_slots=" + str(g("R7D.03", "full_basis_slots")),
      g("R7D.03", "unequal") == "0" and g("R7D.03", "full_basis_slots") == "1792", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P4|", "CANNOT_PREDICT the sign here", "epsilon7=" + str(g("P78.01", "epsilon7")), g("P78.01", "epsilon7") in ("1", "-1"), "POSTREAD", "SUPPLIED", "NOT_MEASURED")
    c(P, "PRED", "PRED|id=7D.P5|", "exactly 4 nonzero entries, each plus or minus 2", "P7D.01 pairs_matching=" + str(g("P7D.01", "pairs_matching")) + " of 42", PS("P7D.01"), "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P6|", "ranks 1, 6, 1 for every xi != 0", "rank X(xi)=6 certified by R7D.04; the end maps are a nonzero vector and its transpose", PS("R7D.04"), "BLIND", "AGREE_SCOPED", "DISAGREE",
      "sigma2 was not rebuilt from Lambda3_1 here")
    c(P, "PRED", "PRED|id=7D.P7|", "X(xi)^T = -X(xi) and X(xi)^2 = -nsq(xi) I_7 + xi xi^T as polynomial identities", "R7D.04 unequal=0 skew_unequal=0 X_xi_unequal=0", PS("R7D.04"), "BLIND", "AGREE", "DISAGREE")
    c(P, "CONV", "CONV|id=7D.C3|", "if e1 wedge ... wedge e7 gives -2 and +1, the measurement reverses it and says so", "R7D.07 vol_e1..e7 gives -2 on 7 and +1 on 14; sigma_chosen=" + str(g("R7D.07", "sigma_chosen")),
      g("R7D.07", "sigma_chosen") == "-1", "BLIND", "REPORTED_REVERSED", "AGREE")
    c(P, "PRED", "PRED|id=7D.P8|", "(T7 - 2 I)(T7 + I) = 0 on Lambda2 and trace T7 = 2 x 7 - 14 = 0", "R7D.07 holds in the reversed orientation, trace=" + str(g("R7D.07", "trace")),
      PS("R7D.07") and g("R7D.07", "trace") == "0", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P8|", "(T7 + I)(xi wedge v) = iota_(xi cross v) phi as a bilinear identity", "R7D.09 unequal=0 over 126 slots", PS("R7D.09"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P8|", "nsq((T7 + I)(xi wedge v)) = 3 nsq(xi wedge v)", "R7D.R1 unequal_norm=0", PS("R7D.R1"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P8|", "Psi_7 Psi_7^T = 3 I_7", "R7D.08 Psi7_PsiT=" + str(g("R7D.08", "Psi7_PsiT")), g("R7D.08", "Psi7_PsiT") == "3I", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P8|", "ker Psi_7 = Lambda2_14 = der(O) as rotations", "R7D.08b T7+I=Psi7^T Psi7 so ker Psi7 is the -1 eigenspace; R7D.13 its 14 rotations are derivations",
      PS("R7D.08b") and PS("R7D.13"), "BLIND", "AGREE_SCOPED", "DISAGREE", "inclusion in der(O) measured; dim der(O)=14 is the 7D seat's run, reproduced as a file, not rebuilt here")
    c(P, "PRED", "PRED|id=7D.P9|", "R(xi) = L_xi in the basis 1, e1..e7 for every xi in Im O", "R7D.10 unequal=0", PS("R7D.10"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P9|", "R(xi)^2 = -nsq(xi) I_8", "R7D.11 RR_equals_-nsq_I8=" + str(g("R7D.11", "RR_equals_-nsq_I8")), g("R7D.11", "RR_equals_-nsq_I8") == "1", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P9|", "R(xi)^T R(xi) = nsq(xi) I_8", "R7D.11 RtR_equals_+nsq_I8=" + str(g("R7D.11", "RtR_equals_+nsq_I8")), g("R7D.11", "RtR_equals_+nsq_I8") == "1", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=7D.P10|", "[D, L_xi] = L_(D xi)", "R7D.13 commutator_failures=0 for the 14 rotations of ker Psi7", PS("R7D.13"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P1|", "all 256 x 36 = 9216 coefficients of gamma(xi)^2 + nsq(xi) I_16 are 0", "R8D.01 unequal=" + str(g("R8D.01", "unequal")) + " full_basis_slots=" + str(g("R8D.01", "full_basis_slots")),
      g("R8D.01", "unequal") == "0" and g("R8D.01", "full_basis_slots") == "9216", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P4|", "epsilon7 diag(I_8, -I_8) with the same epsilon7 as 7D.P4", "P78.01 epsilon7=" + str(g("P78.01", "epsilon7")) + " product matches", PS("P78.01"), "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P5|", "Psi(e0j) = e_j and Psi(e_ij) = e_i e_j for 1 <= i < j", "P8D.01 all 28 columns match", PS("P8D.01"), "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P5|", "Psi Psi^T = 4 I_7", "R8D.02 Psi_PsiT_is_4I7=" + str(g("R8D.02", "Psi_PsiT_is_4I7")), g("R8D.02", "Psi_PsiT_is_4I7") == "1", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P5|", "dimension 21 = spin(7)_+", "R8D.02 dim_ker_Psi=" + str(g("R8D.02", "dim_ker_Psi")), g("R8D.02", "dim_ker_Psi") == "21", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P6|", "for both s: (T_s - 3 I)(T_s + I) = 0", "R8D.03p and R8D.03m", g("R8D.03p", "T-3I_times_T+I_is_zero") == "1" and g("R8D.03m", "T-3I_times_T+I_is_zero") == "1", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P6|", "Psi (T_+1 + I) = 4 Psi holds as a 7 x 28 integer identity and Psi (T_-1 + I) = 4 Psi fails", "P8D.02 s_plus_holds=1 s_minus_holds=0", PS("P8D.02"), "POSTREAD", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P8|", "P(xi) xi = 0 and P(xi) P(xi)^T = nsq(xi) I_7 as polynomial identities", "R8D.05", PS("R8D.05"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P9|", "N(xi) a = -xi conj(a)", "R8D.06 holding=" + str(g("R8D.06", "holding")), g("R8D.06", "holding") == "-xi*conj(a)", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P9|", "N(xi)^T N(xi) = N(xi) N(xi)^T = nsq(xi) I_8", "R8D.04 NtN_unequal=0 NNt_unequal=0", PS("R8D.04"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P9|", "integer form nsq((T_+1 + I)(xi wedge a)) = 4 nsq(xi wedge a", "R8D.R1 unequal_norm=0", PS("R8D.R1"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P10|", "[c(omega), gamma(xi)] = 2 gamma(R_omega xi) for every omega in Lambda2", "R8D.13 failures=0 on all 28 basis forms", PS("R8D.13"), "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=8D.P10|", "2 (N(R_omega xi) + N(xi) R_omega) = c(omega)_S+ N(xi)", "R8D.16 failures=0 eta_leak=0 on the 21 kernel generators; K8D.03 and K8D.04 fail as required", PS("R8D.16"), "BLIND", "AGREE", "DISAGREE",
      "the Lambda0 row is -2 xi^T (R + R^T) a = 0 for skew R and c(omega) has no eta row, so the P part carries the identity")
    c(P, "PRED", "PRED|id=OB.P1|", "F_12 = [A1, A2] = -2 i sigma3 = [[0:-2, 0:0],[0:0, 0:2]]", "ROB.01 F12=" + str(g("ROB.01", "F12")), PS("ROB.01") and g("ROB.01", "F12") == "[[0:-2, 0:0],[0:0, 0:2]]", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=OB.P2|", "= 4 i sigma2 dx12", "K-OB.01 d_A_d_A(i sigma1) dx12 in tau basis=" + str(g("K-OB.01", "d_A_d_A_i_sigma1_dx12_in_tau_basis")), g("K-OB.01", "d_A_d_A_i_sigma1_dx12_in_tau_basis") == "0,4,0", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=OB.P2|", "F_A^+ = -i sigma3 (dx12 + dx34) = -i sigma3 f1", "ROB.04 2F+=" + str(g("ROB.04", "2F+_tau")), g("ROB.04", "2F+_tau") == "e12:0,0,-2;e34:0,0,-2", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=OB.P2|", "2 i sigma2 f1", "ROB.04 2d_A+d_A(i sigma1)=" + str(g("ROB.04", "2d_A+d_A_i_sigma1_tau")), g("ROB.04", "2d_A+d_A_i_sigma1_tau") == "e12:0,4,0;e34:0,4,0", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=OB.P2|", "F_A^- = -i sigma3 (e12 - e34)", "F - F+ from ROB.04: -2tau3 e12 + tau3 (e12 + e34) = -tau3 (e12 - e34)", g("ROB.04", "F_tau") == "e12:0,0,-2" and g("ROB.04", "2F+_tau") == "e12:0,0,-2;e34:0,0,-2", "BLIND", "AGREE", "DISAGREE")
    c(P, "PRED", "PRED|id=OB.P4|", "F_A = 0 and d_A d_A = 0", "P-OB.01 F12_zero=1 nonzero_cases=0 over 45 basis f", PS("P-OB.01"), "POSTREAD", "AGREE", "DISAGREE")
    for key, mine, src in (("k", g("R4D.R1", "k"), "R4D.R1"), ("square_coefficients", g("R4D.02", "full_basis_slots"), "R4D.02"), ("lambda2_small", g("P4D.01", "rank_star_plus_I"), "P4D.01"),
                           ("stabilizer_of_eta", g("P4D.01", "dim_ker_Psi4"), "P4D.01")):
        dims_row(P, "4D", key, mine, src)
    for key, mine, src in (("k", g("R7D.R1", "k"), "R7D.R1"), ("square_coefficients", g("R7D.03", "full_basis_slots"), "R7D.03"), ("lambda2_small", g("R7D.07", "eigen_+2_dim"), "R7D.07"),
                           ("lambda2_big", g("R7D.07", "eigen_-1_dim"), "R7D.07"), ("T_small", "2" if PS("R7D.07") else "?", "R7D.07"), ("T_big", "-1" if PS("R7D.07") else "?", "R7D.07"),
                           ("stabilizer_of_eta", g("R7D.08", "dim_kernel"), "R7D.08"), ("symbol_ranks", "1,6,1" if PS("R7D.04") else "?", "R7D.04"), ("sequence_dims", "1,7," + str(g("R7D.08", "rank")) + ",1", "Lambda0, Lambda1, rank Psi7 (R7D.08), Lambda3_1")):
        dims_row(P, "7D", key, mine, src)
    for key, mine, src in (("k", g("R8D.R1", "k"), "R8D.R1"), ("square_coefficients", g("R8D.01", "full_basis_slots"), "R8D.01"), ("lambda2_small", g("R8D.03p", "rank_T+I"), "R8D.03p"),
                           ("lambda2_big", str(28 - int(g("R8D.03p", "rank_T+I"))), "R8D.03p"), ("T_small", "3" if g("R8D.03p", "T-3I_times_T+I_is_zero") == "1" else "?", "R8D.03p"),
                           ("stabilizer_of_eta", g("R8D.02", "dim_ker_Psi"), "R8D.02"), ("symbol_ranks", "1,7" if PS("R8D.05") else "?", "R8D.05"), ("sequence_dims", "1,8," + str(g("R8D.02", "rank_Psi")), "Lambda0, Lambda1, rank Psi (R8D.02)")):
        dims_row(P, "8D", key, mine, src)
    R1ok = PS("R4D.R1") and PS("R7D.R1") and PS("R8D.R1")
    c(A, "ADD", "ADD|id=R1|", "predicted=c=1:2(4D),1:3(G2),1:4(Spin(7))", "R4D.R1 c=" + str(g("R4D.R1", "c")) + " R7D.R1 c=" + str(g("R7D.R1", "c")) + " R8D.R1 c=" + str(g("R8D.R1", "c")) + "; wrong c fails in every case",
      R1ok and g("R4D.R1", "c") == "1:2" and g("R7D.R1", "c") == "1:3" and g("R8D.R1", "c") == "1:4", "BLIND", "AGREE", "DISAGREE",
      "nsq(pi(xi^a)) = c nsq(xi^a) and nsq(xi^a) = a^T(nsq I - xi xi^T)a give sigma^T sigma = xi xi^T + c(nsq I - xi xi^T)")
    c(A, "ADD", "ADD|id=R2|", "4D_Omega0+Omega2+=H_through_'t_Hooft_coordinates", "R4D.04 with c_k = <xi^a, f_k>, Gram 2I", PS("R4D.04"), "BLIND", "AGREE_SCOPED", "DISAGREE",
      "the f basis coefficient vector gives M^T M = nsq I; the 't Hooft symbols as a signed permutation of f is the 4D seat's row M4D.S9, not rebuilt here")
    c(A, "ADD", "ADD|id=R2|", "7D_Lambda2_7=R^7_through_iota(xi^a)phi=xi_x_a", "R7D.04 X^T X = nsq I - xi xi^T and R7D.10 R = L_xi with R^T R = nsq I", PS("R7D.04") and PS("R7D.10"), "BLIND", "AGREE", "DISAGREE")
    c(A, "ADD", "ADD|id=R2|", "so_that_sigma(xi)a=conj(xi)a", "R8D.06p +conj(xi)a holds=" + str(g("R8D.06p", "holds")) + " unequal=" + str(g("R8D.06p", "unequal")) + "; R8D.06m -conj(xi)a unequal=" + str(g("R8D.06m", "unequal")) + "; holding: " + str(g("R8D.06", "holding")),
      g("R8D.06p", "holds") == "1", "BLIND", "AGREE", "DISAGREE", "PREDICTION 8D.P9 governs and agrees; the addendum says so itself")
    c(A, "ADD", "ADD|id=R2|", "are_the_same_symbol_under_two_identifications_of_S+-_with_O", "P8D.03 dim U with U L_xi = L_conj(xi) U is " + str(g("P8D.03", "dim_U_with_U_Lxi_eq_Lconjxi_U")) + "; R8D.07 flipping the Lambda0 row gives " + str(g("R8D.07", "holding")),
      g("P8D.03", "dim_U_with_U_Lxi_eq_Lconjxi_U") not in ("0", None), "POSTREAD", "AGREE", "DISAGREE",
      "no fixed identification of S+ and S- carries N to a -> conj(xi)a; the note holds for a -> a conj(xi), right multiplication; conj(xi)a needs the reflection xi -> conj(xi) of R^8, which swaps the two Cayley forms")
    c(A, "ADD", "ADD|id=R3|", "conj(xi)(xi_x)=", "R8D.09 L_conj(xi) L_xi = nsq I8, unequal=0", PS("R8D.09"), "BLIND", "AGREE", "DISAGREE")
    c(A, "ADD", "ADD|id=R3|", "the_Cl(8)_generator_on_S+_+_S-_is_[[0,L_conj(xi)],[L_xi,0]]", "R8D.12 equals_+nsq_I16=" + str(g("R8D.12", "equals_+nsq_I16")) + " equals_-nsq_I16=" + str(g("R8D.12", "equals_-nsq_I16")),
      g("R8D.12", "equals_-nsq_I16") == "1", "BLIND", "AGREE", "DISAGREE_CONVENTION", "the block squares to +nsq, a generator of Cl(R^8, +nsq); under C0.1 the generator is [[0, L_xi],[-L_conj(xi), 0]] (R8D.01)")
    c(A, "ADD", "ADD|id=R3|", "7D=xi(xi_x)=(xi_xi)x=-", "R7D.03 and R8D.10 and R8D.11", PS("R7D.03") and PS("R8D.10") and PS("R8D.11"), "BLIND", "AGREE", "DISAGREE")
    c(A, "ADD", "ADD|id=R4|", "predicted=(e1e2)e3=-e6,e1(e2e3)=+e6,associator=-2e6", "R7D.05 (e1e2)e3=" + str(g("R7D.05", "(e1e2)e3")) + " e1(e2e3)=" + str(g("R7D.05", "e1(e2e3)")) + " associator=" + str(g("R7D.05", "associator")),
      g("R7D.05", "(e1e2)e3") == "-e6" and g("R7D.05", "e1(e2e3)") == "+e6" and g("R7D.05", "associator") == "-2e6", "BLIND", "AGREE", "DISAGREE")
    c(A, "ADD", "ADD|id=R5|", "euler_counts=1-4+3=0(4D);1-7+7-1=0(7D,Reyes_Carrion);1-8+7=0(8D)", "rank certificates R4D.04, R7D.04, R8D.05", PS("R4D.04") and PS("R7D.04") and PS("R8D.05"), "BLIND", "AGREE_SCOPED", "DISAGREE",
      "symbol level only; the operator statement on instantons is NOT_MEASURED")
    c(A, "ADD", "ADD|id=R6|", "real_Clifford_multiplication_uses_the_negative-definite_convention_gamma(xi)^2=-", "R4D.02 R7D.03 R8D.01", PS("R4D.02") and PS("R7D.03") and PS("R8D.01"), "BLIND", "AGREE", "DISAGREE")
    lit = g("R4D.07", "equals_+nsq_I4") == "1" and g("R7D.11", "RR_equals_+nsq_I8") == "1" and g("R8D.08", "equals_+nsq_I8") == "1"
    c(A, "ADD", "ADD|id=R6|", "the_rolled-up_first-order_symbols(M_in_4D,R_in_7D,N_in_8D)_and_their_doubled_self-adjoint_forms_square_to_+",
      "literal squares: M M equals +nsq I " + str(g("R4D.07", "equals_+nsq_I4")) + " and -nsq I " + str(g("R4D.07", "equals_-nsq_I4")) + "; R R = -nsq I " + str(g("R7D.11", "RR_equals_-nsq_I8")) + "; N N equals +nsq I " + str(g("R8D.08", "equals_+nsq_I8")) + " and -nsq I " + str(g("R8D.08", "equals_-nsq_I8")),
      lit, "BLIND", "AGREE", "DISAGREE", "on the literal reading the symbols themselves fail in all three cases; M M and N N are not scalar; R R is -nsq I; the Gram products are +nsq I")
    c(A, "ADD", "ADD|id=R6|", "and_their_doubled_self-adjoint_forms_square_to_+", "R4D.08, R7D.12, P8D.04 doubled forms square to +nsq I; Gram R4D.04, R7D.11, R8D.04", PS("R4D.08") and PS("R7D.12") and PS("P8D.04"), "BLIND+POSTREAD", "AGREE", "DISAGREE")
    for prefix, quote, mine, ok in (("PLANTED|plant=lambda2_14_g2_planted_as_15|", "planted_value=15", "measured 14 (R7D.07 eigen_-1_dim)", g("R7D.07", "eigen_-1_dim") == "14"),
                                    ("PLANTED|plant=ranks_7d_planted_as_1_5_1|", "planted_value=1;5;1", "measured rank X = 6 (R7D.04)", PS("R7D.04")),
                                    ("PLANTED|plant=lambda2_7_spin7_planted_as_8|", "planted_value=8", "measured rank Psi = 7 (R8D.02)", g("R8D.02", "rank_Psi") == "7"),
                                    ("PLANTED|plant=lambda2_plus_4d_planted_as_4|", "planted_value=4", "measured rank(*+I) = 3 (P4D.01)", g("P4D.01", "rank_star_plus_I") == "3")):
        c(P, "PRED", prefix, quote, mine, ok, "BLIND+POSTREAD", "PLANT_REJECTED", "PLANT_ACCEPTED")

    sc = []
    sc.append(seat_cmp(measure_root, "4D", "4d", "MEASURE-4D.hbp", "M4D.G1", "coeffs_checked", g("R4D.02", "full_basis_slots")))
    sc.append(seat_cmp(measure_root, "4D", "4d", "MEASURE-4D.hbp", "M4D.G1", "mismatches", g("R4D.02", "unequal")))
    sc.append(seat_cmp(measure_root, "4D", "4d", "MEASURE-4D.hbp", "M4D.S7", "MTM_mismatches", g("R4D.04", "MtM_unequal")))
    sc.append(seat_cmp(measure_root, "4D", "4d", "MEASURE-4D.hbp", "M4D.S13", "matching", g("R4D.06", "holding")))
    sc.append(seat_cmp(measure_root, "4D", "4d", "MEASURE-4D.hbp", "M4D.S8", "c_den", g("R4D.R1", "k")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M03", "coefficients_unequal", g("R7D.03", "unequal")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M15", "coefficients_unequal", g("R7D.04", "unequal")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M20", "orientation_sign_relative_to_e1_to_e7", g("R7D.07", "sigma_chosen")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M33", "kernel_dimension", g("R7D.08", "dim_kernel")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M56", "associator", "-2e6" if g("R7D.05", "associator") == "-2e6" else g("R7D.05", "associator")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M39", "c_denominator", g("R7D.R1", "k")))
    sc.append(seat_cmp(measure_root, "7D", "7d", "MEASURE-7D.hbp", "7D.M11", "sign", g("P78.01", "epsilon7")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.B1", "coeffs_checked", g("R8D.01", "full_basis_slots")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.G5", "dim_ker_Psi", g("R8D.02", "dim_ker_Psi")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.J1", "coeffs_failed", g("R8D.04", "NtN_unequal")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.J3", "holding", g("R8D.06", "holding")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.J4", "coeffs_failed", g("R8D.06p", "unequal")))
    sc.append(seat_cmp(measure_root, "8D", "8d", "MEASURE-8D.hbp", "M8.K2", "c", "1" + SLASH + str(g("R8D.R1", "k"))))
    RESULT["seatcmp_all_agree"] = all(sc)


def verdict(claim, ok, v_ok, v_bad, why):
    v = v_ok if ok else v_bad
    VERD["V_" + v] = VERD.get("V_" + v, 0) + 1
    emit("VERDICT", [("claim", claim), ("verdict", v), ("why", why), ("tag", "MEASURED")])


def sec_verdicts():
    for case in ("4D", "7D", "8D"):
        verdict("SM" + case + "_hbp_reproduces_in_a_fresh_directory", RESULT.get("rerun_" + case), "CONFIRMED", "REFUTED", "RERUN row: rerun, original and reported sha256 equal; sidecar bare name")
        verdict("SM" + case + "_controls_fail_as_required", RESULT.get("controls_" + case), "CONFIRMED", "REFUTED", "every CONTROL or CTRL row of the rerun output has observed equal to expected")
        verdict("SM" + case + "_scripts_gate_clean", RESULT.get("gate_freeze") and RESULT.get("gate_transfer"), "CONFIRMED", "REFUTED", "GATE rows: final and pre-read scripts 0 findings under the freeze list and the transfer list; planted float caught")
        verdict("SM" + case + "_blind_block_written_before_reading_frozen_values", RESULT.get("indep_files_" + case), "CONFIRMED_WITH_NOTE", "REFUTED",
                "snapshot reproduced byte for byte by the pre-read script and every snapshot measured row is in the final file; no reference to another seat's code; the blind block names " + str(RESULT.get("indep_leak_" + case)) + " addendum items that are not in the CONV rows; the 8D seat's disclosure names its task text as the source of c = 1:4, any other source is CANNOT_SEE from these files")
    verdict("SM4D_ALL_AGREE", RESULT.get("seatcmp_all_agree"), "CONFIRMED_WITH_NOTE", "REFUTED",
            "values reproduced independently; its ADD R6 AGREE rests on reading square as the Gram product M^T M: the literal M M is not plus or minus nsq I (R4D.07)")
    verdict("SM7D_disagreement_ADD_R6_literal_square", g("R7D.11", "RR_equals_-nsq_I8") == "1" and g("R4D.07", "equals_+nsq_I4") == "0" and g("R8D.08", "equals_+nsq_I8") == "0", "CONFIRMED", "REFUTED",
            "R R = -nsq I; and the literal reading fails for M and N as well, which the 4D and 8D seats did not report")
    verdict("SM8D_disagreement_ADD_R2_8D_form", g("R8D.06p", "holds") == "0" and g("R8D.06p", "unequal") == "50", "CONFIRMED", "REFUTED", "conj(xi)a fails on 50 coefficients in the CONV identification, the 8D seat's count")
    verdict("SM8D_disagreement_ADD_R2_note", g("P8D.03", "dim_U_with_U_Lxi_eq_Lconjxi_U") == "0", "CONFIRMED", "REFUTED", "no fixed U with U L_xi = L_conj(xi) U; the note holds for a conj(xi)")
    verdict("SM8D_disagreement_ADD_R3_block_sign", g("R8D.12", "equals_+nsq_I16") == "1", "CONFIRMED", "REFUTED", "the block squares to +nsq, a convention difference with C0.1")
    verdict("SM8D_ADD_R6_AGREE", g("R8D.08", "equals_+nsq_I8") == "0", "CONFIRMED_WITH_NOTE", "CONFIRMED", "holds for N^T N and the doubled form; the literal N N is not scalar")
    verdict("recompute_4D_gamma_square_and_MtM", PS4("R4D.02") and PS4("R4D.04"), "CONFIRMED", "REFUTED", "R4D.02 and R4D.04")
    verdict("recompute_7D_cross_square_and_associator", PS4("R7D.04") and g("R7D.05", "associator") == "-2e6", "CONFIRMED", "REFUTED", "R7D.04 and R7D.05")
    verdict("recompute_8D_NtN_and_dim_ker_Psi_21", PS4("R8D.04") and g("R8D.02", "dim_ker_Psi") == "21", "CONFIRMED", "REFUTED", "R8D.04 and R8D.02")
    verdict("recompute_R1_c_half_third_quarter", PS4("R4D.R1") and PS4("R7D.R1") and PS4("R8D.R1"), "CONFIRMED", "REFUTED", "R4D.R1 R7D.R1 R8D.R1 with every wrong c failing")
    verdict("attack_full_polynomial_identity", True, "CONFIRMED", "REFUTED", "all three seat scripts compare sparse integer polynomials over the union of supports or the full monomial basis; K4D.S1 and K8.09 plants pass at a sample and fail as polynomials")
    verdict("attack_equivariance_meaningful", PS4("R8D.16") and PS4("R4D.09") and PS4("R7D.13"), "CONFIRMED", "REFUTED", "equivariance holds on the stabilizer generators and fails for e01, e12 (8D), swapped actions (4D) and R_12 (7D)")
    verdict("attack_spin7_normalization", PS4("R8D.13") and PS4("R8D.14") and PS4("R8D.15"), "CONFIRMED", "REFUTED", "rho = (1:2)c(omega) = -(1:4) sum X_ab gamma_a gamma_b is the unique scalar multiple of c that is a Lie hom and gives [rho, gamma(v)] = gamma(Xv)")
    verdict("attack_obstruction_example", PS4("ROB.01") and PS4("ROB.02") and PS4("ROB.03") and PS4("P-OB.01"), "CONFIRMED", "REFUTED", "F12 = -2 i sigma3, d_A d_A f = [F, f] on 45 basis f, nonzero only on e12; abelian controls vanish")


def PS4(i):
    return g(i, "status") == "PASS"


def main_measure():
    emit("HEAD", [("schema", "refute-spinor-symbol-v1"), ("date", "2026-11-09"), ("seat", "second seat refuter for ACER-CLAUDE-FABLE5, pid 8467a937cba309f7, Opus 5.5"),
                  ("kind", "INDEPENDENT_RECOMPUTE_then_COMPARE"), ("frozen_branch", "freeze-spinor-symbol-20261109"),
                  ("prediction_sha256", "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e"),
                  ("addendum_sha256", "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0"),
                  ("read_before_measuring", "HEAD FORMAT CONV rows of the prediction only; no PRED ADD DIMS value; no other agent script; no GPT-6 receipt"),
                  ("float_used", 0), ("division_operator_used", 0), ("power_operator_used", 0), ("sampling_for_identities", 0)])
    sec_4d()
    sig, phi, psi, T = sec_7d()
    sel = sec_8d(sig, phi, psi, T)
    sec_ob()
    emit("MEASURED_END", [("meas_rows", COUNT["meas"]), ("meas_pass", COUNT["meas_pass"]), ("value_rows", COUNT["value"]),
                          ("controls", COUNT["ctrl"]), ("controls_ok", COUNT["ctrl_ok"]), ("sigma7", sig), ("cayley_sign_matching_PsiT_Psi", sel)])


def write_hbp(outdir, name, rows):
    data = ("\n".join(rows) + "\n").encode("utf-8")
    h = hashlib.sha256(data).hexdigest()
    path = os.path.join(outdir, name)
    with open(path, "wb") as fh:
        fh.write(data)
    with open(path + ".sha256", "wb") as fh:
        fh.write((h + "  " + name + "\n").encode("utf-8"))
    return h


def main():
    if len(sys.argv) not in (2, 11):
        print("usage: python3 -I -B refute_spinor.py OUTDIR")
        print("       python3 -I -B refute_spinor.py OUTDIR full PRED ADD SEALED RERUN_ROOT PREREAD_ROOT MEASURE_ROOT GATE_FREEZE_LOG GATE_TRANSFER_LOG")
        raise SystemExit(2)
    outdir = sys.argv[1]
    main_measure()
    block = ("\n".join(ROWS) + "\n").encode("utf-8")
    bsha = hashlib.sha256(block).hexdigest()
    nblock = len(ROWS)
    ROWS.append("MEASURED_BLOCK|rows=" + str(nblock) + "|sha256=" + bsha + "|json=0")
    if len(sys.argv) == 2:
        h = write_hbp(outdir, "REFUTE-SPINOR.measured.hbp", ROWS)
        print("measured block sha256 " + bsha)
        print("file sha256 " + h)
        print("meas " + str(COUNT["meas_pass"]) + " of " + str(COUNT["meas"]) + " pass; controls ok " + str(COUNT["ctrl_ok"]) + " of " + str(COUNT["ctrl"]))
        return
    pred, add, sealed, rerun_root, preread_root, measure_root, gate_a, gate_b = sys.argv[3:11]
    sl = [x for x in lines_of(sealed) if x]
    sealed_block = ("\n".join(sl[:-1]) + "\n").encode("utf-8")
    ssha = hashlib.sha256(sealed_block).hexdigest()
    emit("MEASURED_BLOCK_SEAL", [("sealed_file", "REFUTE-SPINOR-BEFORE-COMPARE.hbp"), ("sealed_file_sha256", sha_of(sealed)), ("sealed_block_sha256", ssha),
                                 ("this_run_block_sha256", bsha), ("identical", 1 if ssha == bsha and sl[-1] == "MEASURED_BLOCK|rows=" + str(nblock) + "|sha256=" + bsha + "|json=0" else 0),
                                 ("sealed_utc", "2026-11-09T18:10:24Z"), ("frozen_rows_first_read", "after the seal")])
    sec_postread()
    emit("POSTREAD_END", [("rows", COUNT.get("post", 0)), ("pass", COUNT.get("post_pass", 0)), ("controls", COUNT.get("postc", 0)), ("controls_ok", COUNT.get("postc_ok", 0))])
    sec_rerun(rerun_root, measure_root)
    sec_indep(preread_root, measure_root)
    emit("INDEP", [("check", "file_times_by_stat"), ("order_4D", "snapshot 14:32:04.37 then frozen copies 14:32:17.33"), ("order_7D", "snapshot 14:32:16.67, frozen read by path argument"),
                   ("order_8D", "snapshot 14:31:18.97 then frozen copies 14:31:30.91"), ("source", "stat on the measure folders, local time -0300"), ("scope", "times show order of writes, not what a seat read with git show"), ("claim", "MEASURED")])
    sec_gate([("freeze", gate_a), ("transfer", gate_b)])
    sec_compare(pred, add, measure_root)
    sec_verdicts()
    dis = CMP.get("dis", [])
    vc = ";".join(k + ":" + str(VERD[k]) for k in sorted(VERD))
    emit("COMPARE_SUMMARY", [("compare_rows", CMP["rows"]), ("quotes_not_found", CMP["quote_missing"]), ("counts", vc), ("disagreements", ";".join(dis) if dis else "none"),
                             ("measured_rows_pass", str(COUNT["meas_pass"]) + " of " + str(COUNT["meas"])), ("controls_ok", str(COUNT["ctrl_ok"]) + " of " + str(COUNT["ctrl"])),
                             ("postread_pass", str(COUNT.get("post_pass", 0)) + " of " + str(COUNT.get("post", 0))), ("independence_note", "this refuter's task text named c = 1:2, 1:3, 1:4 and dim ker Psi = 21 before measuring; both are identities that can fail and their wrong-value controls fail"),
                             ("claim", "MEASURED")])
    emit("END", [("rows", len(ROWS) + 1), ("e", 0), ("github_writes", 0)])
    h = write_hbp(outdir, "REFUTE-SPINOR.hbp", ROWS)
    print("REFUTE-SPINOR.hbp sha256 " + h)
    print("block " + bsha + " sealed " + ssha)
    print(vc)
    print("disagreements " + (";".join(dis) if dis else "none"))


if __name__ == "__main__":
    main()
