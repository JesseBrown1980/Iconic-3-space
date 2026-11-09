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
    if len(sys.argv) < 2:
        print("usage: python3 -I -B refute_spinor.py OUTDIR")
        raise SystemExit(2)
    outdir = sys.argv[1]
    main_measure()
    block = ("\n".join(ROWS) + "\n").encode("utf-8")
    bsha = hashlib.sha256(block).hexdigest()
    rows = list(ROWS)
    rows.append("MEASURED_BLOCK|rows=" + str(len(ROWS)) + "|sha256=" + bsha + "|json=0")
    h = write_hbp(outdir, "REFUTE-SPINOR.measured.hbp", rows)
    print("measured block sha256 " + bsha)
    print("file sha256 " + h)
    print("meas " + str(COUNT["meas_pass"]) + " of " + str(COUNT["meas"]) + " pass; controls ok " + str(COUNT["ctrl_ok"]) + " of " + str(COUNT["ctrl"]))


if __name__ == "__main__":
    main()
