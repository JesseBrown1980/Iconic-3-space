"""Independent integer measurement of the 7D case of the frozen spinor-symbol prediction.

Octonions come from the frozen table: e_i e_(i+1) = e_(i+3), labels 1..7 taken mod 7, e_i e_i = -1, e0 = 1.
Every identity in a covector is checked as an identity of integer polynomials, coefficient by coefficient.
Integers only, no division operator; exact quotients go through divmod with a zero remainder asserted.
Run: python3 -I -B measure_7d.py OUTDIR [PREDICTION_HBP ADDENDUM_HBP]
"""
import sys
import os
import hashlib
import itertools
import math

ROWS = []
LINES_FROZEN = ((1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3))
LABELS = (1, 2, 3, 4, 5, 6, 7)
PAIRS = tuple((i, j) for i in LABELS for j in LABELS if i < j)
PIDX = {p: n for n, p in enumerate(PAIRS)}
VOL = (1, 2, 3, 4, 5, 6, 7)


def emit(tag, kv):
    parts = [tag]
    for k, v in kv:
        s = str(v)
        if "|" in s or "\n" in s or "=" in k:
            raise ValueError("bad field " + k)
        parts.append(k + "=" + s)
    parts.append("json=0")
    ROWS.append("|".join(parts))


MCOUNT = [0]
KCOUNT = [0]


def two(n):
    return ("0" if n < 10 else "") + str(n)


def mrow(name, kv):
    MCOUNT[0] += 1
    emit("M7", [("id", "7D.M" + two(MCOUNT[0])), ("name", name)] + kv + [("claim", "MEASURED")])


def krow(name, failed, kv):
    KCOUNT[0] += 1
    emit("CONTROL", [("id", "7D.K" + two(KCOUNT[0])), ("name", name), ("must", "FAIL"),
                     ("observed", "FAIL" if failed else "PASS"), ("control_ok", 1 if failed else 0)] + kv + [("claim", "MEASURED")])


def exact_div(a, b):
    q, r = divmod(a, b)
    if r != 0:
        raise ValueError("inexact division")
    return q


def normalize(v):
    g = 0
    for x in v:
        g = math.gcd(g, x)
    if g > 1:
        return [exact_div(x, g) for x in v]
    return list(v)


# ---------- integer polynomials: dict sorted-variable-tuple -> int ----------

def pc(c):
    return {(): c} if c else {}


def pv(i):
    return {(i,): 1}


def padd(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        t = r.get(m, 0) + s * c
        if t:
            r[m] = t
        elif m in r:
            del r[m]
    return r


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2))
            t = r.get(m, 0) + c1 * c2
            if t:
                r[m] = t
            elif m in r:
                del r[m]
    return r


def pscale(p, k):
    if k == 0:
        return {}
    return {m: k * c for m, c in p.items()}


def psum(ps):
    r = {}
    for p in ps:
        r = padd(r, p)
    return r


def peval(p, assign):
    tot = 0
    for m, c in p.items():
        t = c
        for x in m:
            t = t * assign.get(x, 0)
        tot += t
    return tot


# ---------- polynomial matrices ----------

def pm_int(A):
    return [[pc(x) for x in r] for r in A]


def pm_zero(n, m):
    return [[{} for _ in range(m)] for _ in range(n)]


def pm_eye(n, p):
    return [[(dict(p) if i == j else {}) for j in range(n)] for i in range(n)]


def pm_mul(A, B):
    n = len(A)
    k = len(B)
    m = len(B[0])
    out = []
    for i in range(n):
        rv = []
        for j in range(m):
            acc = {}
            for t in range(k):
                a = A[i][t]
                b = B[t][j]
                if a and b:
                    acc = padd(acc, pmul(a, b))
            rv.append(acc)
        out.append(rv)
    return out


def pm_add(A, B, s=1):
    return [[padd(A[i][j], B[i][j], s) for j in range(len(A[0]))] for i in range(len(A))]


def pm_scale(A, k):
    return [[pscale(x, k) for x in r] for r in A]


def pm_T(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def pm_imul(A, B):
    # integer matrix A times polynomial matrix B
    n = len(A)
    k = len(B)
    m = len(B[0])
    out = []
    for i in range(n):
        rv = []
        for j in range(m):
            acc = {}
            for t in range(k):
                if A[i][t] and B[t][j]:
                    acc = padd(acc, B[t][j], A[i][t])
            rv.append(acc)
        out.append(rv)
    return out


def compare(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("shape")
    entries = 0
    coeffs = 0
    bad = 0
    mx = 0
    for i in range(len(A)):
        for j in range(len(A[0])):
            a = A[i][j]
            b = B[i][j]
            entries += 1
            keys = set(a) | set(b)
            coeffs += len(keys)
            for m in keys:
                d = a.get(m, 0) - b.get(m, 0)
                if d:
                    bad += 1
                    if abs(d) > mx:
                        mx = abs(d)
    return entries, coeffs, bad, mx


def ident(name, A, B, variables, statement, extra=()):
    e, c, bad, mx = compare(A, B)
    holds = 1 if bad == 0 else 0
    mrow(name, [("statement", statement), ("holds", holds), ("entries", e), ("coefficients_compared", c),
                ("coefficients_unequal", bad), ("max_abs_residual", mx), ("variables", variables)] + list(extra) +
         [("method", "integer polynomial identity, coefficient by coefficient, no sampling")])
    return holds


def count_terms(A, B):
    # monomial products summed by pm_mul(A, B) before cancellation
    n = 0
    for i in range(len(A)):
        for j in range(len(B[0])):
            for t in range(len(B)):
                n += len(A[i][t]) * len(B[t][j])
    return [("product_terms_summed_before_cancellation", n)]


def fails(A, B):
    e, c, bad, mx = compare(A, B)
    return bad > 0, [("entries", e), ("coefficients_compared", c), ("coefficients_unequal", bad), ("max_abs_residual", mx)]


# ---------- integer matrices ----------

def ieye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def izero(n, m):
    return [[0] * m for _ in range(n)]


def imul(A, B):
    n = len(A)
    k = len(B)
    m = len(B[0])
    out = []
    for i in range(n):
        Ai = A[i]
        rv = [0] * m
        for t in range(k):
            a = Ai[t]
            if a:
                Bt = B[t]
                for j in range(m):
                    if Bt[j]:
                        rv[j] += a * Bt[j]
        out.append(rv)
    return out


def iadd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def iscal(A, k):
    return [[k * x for x in r] for r in A]


def iT(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def iident(name, A, B, statement):
    entries = len(A) * len(A[0])
    bad = 0
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i][j] != B[i][j]:
                bad += 1
    holds = 1 if bad == 0 else 0
    mrow(name, [("statement", statement), ("holds", holds), ("entries", entries), ("entries_unequal", bad),
                ("method", "exact integer matrix identity")])
    return holds


def rank_int(rows):
    M = [normalize(r) for r in rows if any(r)]
    if not M:
        return 0
    ncols = len(M[0])
    pr = 0
    for col in range(ncols):
        piv = -1
        for r in range(pr, len(M)):
            if M[r][col]:
                piv = r
                break
        if piv < 0:
            continue
        M[pr], M[piv] = M[piv], M[pr]
        p = M[pr][col]
        for r in range(pr + 1, len(M)):
            a = M[r][col]
            if a:
                M[r] = normalize([p * x - a * y for x, y in zip(M[r], M[pr])])
        pr += 1
        if pr == len(M):
            break
    return pr


def lcm(a, b):
    a = abs(a)
    b = abs(b)
    return exact_div(a * b, math.gcd(a, b))


def nullspace_int(rows, ncols):
    M = [normalize(r) for r in rows if any(r)]
    pivcols = []
    pr = 0
    for col in range(ncols):
        piv = -1
        for r in range(pr, len(M)):
            if M[r][col]:
                piv = r
                break
        if piv < 0:
            continue
        M[pr], M[piv] = M[piv], M[pr]
        p = M[pr][col]
        for r in range(len(M)):
            if r != pr and M[r][col]:
                a = M[r][col]
                M[r] = normalize([p * x - a * y for x, y in zip(M[r], M[pr])])
        pivcols.append(col)
        pr += 1
    M = M[:pr]
    free = [c for c in range(ncols) if c not in pivcols]
    basis = []
    for f in free:
        L = 1
        for r in range(pr):
            L = lcm(L, M[r][pivcols[r]])
        x = [0] * ncols
        x[f] = L
        for r in range(pr):
            x[pivcols[r]] = -exact_div(M[r][f] * L, M[r][pivcols[r]])
        x = normalize(x)
        for t in x:
            if t:
                if t < 0:
                    x = [-y for y in x]
                break
        basis.append(x)
    return pr, basis


class Echelon:
    def __init__(self):
        self.rows = []

    def reduce(self, v):
        v = list(v)
        for (p, r) in self.rows:
            if v[p]:
                a = v[p]
                b = r[p]
                v = normalize([b * x - a * y for x, y in zip(v, r)])
        return v

    def add(self, v):
        w = self.reduce(v)
        for idx, x in enumerate(w):
            if x:
                self.rows.append((idx, w))
                return 1
        return 0


def flat(A):
    return [x for r in A for x in r]


# ---------- octonions ----------

def rule_lines():
    return tuple((i, (i % 7) + 1, ((i + 2) % 7) + 1) for i in LABELS)


def build_mult(lines, flip=()):
    M = [[None] * 8 for _ in range(8)]
    sets = 0
    for a in range(8):
        M[0][a] = (1, a)
        M[a][0] = (1, a)
    for a in LABELS:
        M[a][a] = (-1, 0)
    for n, (a, b, c) in enumerate(lines):
        s = -1 if n in flip else 1
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            if M[x][y] is not None or M[y][x] is not None:
                raise ValueError("table cell set twice")
            M[x][y] = (s, z)
            M[y][x] = (-s, z)
            sets += 2
    for a in range(8):
        for b in range(8):
            if M[a][b] is None:
                raise ValueError("table hole")
    return M, sets


def left_mats(M):
    Ls = []
    for a in range(8):
        A = izero(8, 8)
        for b in range(8):
            s, c = M[a][b]
            A[c][b] = s
        Ls.append(A)
    return Ls


def omul(M, x, y):
    z = [0] * 8
    for a in range(8):
        if x[a]:
            for b in range(8):
                if y[b]:
                    s, c = M[a][b]
                    z[c] += s * x[a] * y[b]
    return z


def pL(Ls, x):
    out = pm_zero(8, 8)
    for a in range(8):
        if x[a]:
            for r in range(8):
                for c in range(8):
                    if Ls[a][r][c]:
                        out[r][c] = padd(out[r][c], x[a], Ls[a][r][c])
    return out


def unit(a):
    v = [0] * 8
    v[a] = 1
    return v


def vec_str(v):
    parts = []
    for a in range(8):
        if v[a]:
            parts.append(("+" if v[a] > 0 else "") + str(v[a]) + "e" + str(a))
    return "".join(parts) if parts else "0"


# ---------- exterior forms on R^7 with polynomial coefficients ----------

def perm_sign(seq):
    inv = 0
    for a in range(len(seq)):
        for b in range(a + 1, len(seq)):
            if seq[a] > seq[b]:
                inv += 1
    return -1 if inv % 2 else 1


def wedge_basis(I, J):
    if set(I) & set(J):
        return 0, None
    seq = list(I) + list(J)
    return perm_sign(seq), tuple(sorted(seq))


def fclean(C):
    return {K: c for K, c in C.items() if c}


def fwedge(A, B):
    C = {}
    for I, a in A.items():
        for J, b in B.items():
            s, K = wedge_basis(I, J)
            if s:
                C[K] = padd(C.get(K, {}), pmul(a, b), s)
    return fclean(C)


def fhodge(A, orient):
    C = {}
    for I, a in A.items():
        Ic = tuple(x for x in LABELS if x not in I)
        s = perm_sign(list(I) + list(Ic)) * orient
        C[Ic] = padd(C.get(Ic, {}), a, s)
    return fclean(C)


def finner(A, B):
    r = {}
    for I, a in A.items():
        if I in B:
            r = padd(r, pmul(a, B[I]))
    return r


def one_form(vals):
    return fclean({(a,): vals[a - 1] for a in LABELS})


def two_vec(F):
    return [F.get(p, {}) for p in PAIRS]


def form_from_two_vec(vec):
    return fclean({PAIRS[b]: vec[b] for b in range(21)})


def col(vec):
    return [[x] for x in vec]


def measure():
    emit("MEASURE_HEAD", [("schema", "measure-spinor-symbol-7d-v1"), ("date", "2026-11-09"),
                          ("seat", "local Opus 5.5 seat, measurement agent for case 7D, pid 8467a937cba309f7"),
                          ("kind", "MEASUREMENT"), ("case", "7D"),
                          ("independence", "every matrix built here from definitions; the CONV rows of the frozen prediction were read first; no PRED, ADD or DIMS value was read when this block was written; no other agent script or receipt read"),
                          ("frozen_branch", "freeze-spinor-symbol-20261109"),
                          ("prediction_sha256", "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e"),
                          ("addendum_sha256", "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0"),
                          ("conventions_used", "C0.1 C0.2 C0.3 C0.4 C0.5 C0.6 C0.7 C0.8 C0.9 7D.C1 7D.C2 7D.C3 7D.C4 7D.C5"),
                          ("float_used", 0), ("division_operator_used", 0), ("sampling_used_for_identities", 0)])

    # ---- table ----
    lines = rule_lines()
    M, sets = build_mult(lines)
    pair_count = {}
    for (a, b, c) in lines:
        for (x, y) in ((a, b), (b, c), (a, c)):
            q = (min(x, y), max(x, y))
            pair_count[q] = pair_count.get(q, 0) + 1
    fano = sum(1 for q in pair_count if pair_count[q] == 1)
    anti = sum(1 for i in LABELS for j in LABELS if i != j and M[i][j][1] == M[j][i][1] and M[i][j][0] == -M[j][i][0])
    sq = sum(1 for i in LABELS if M[i][i] == (-1, 0))
    mrow("octonion_table", [("rule", "e_i e_(i+1) = e_(i+3) on labels 1..7 mod 7, e_i e_i = -1, e0 = 1"),
                            ("lines", ";".join(",".join(str(x) for x in ln) for ln in lines)),
                            ("lines_equal_frozen_list", 1 if lines == LINES_FROZEN else 0),
                            ("pairs_on_exactly_one_line", str(fano) + " of 21"),
                            ("antisymmetric_pairs", str(anti) + " of 42"), ("squares_minus_one", str(sq) + " of 7"),
                            ("cells_set_by_lines", sets), ("cells_total", 64)])

    Ls = left_mats(M)
    I8 = ieye(8)
    ok = 0
    for i in LABELS:
        for j in LABELS:
            S = iadd(imul(Ls[i], Ls[j]), imul(Ls[j], Ls[i]))
            if S == iscal(I8, -2 if i == j else 0):
                ok += 1
    mrow("clifford_relations_on_basis", [("statement", "L_ei L_ej + L_ej L_ei = -2 delta_ij I_8 for i, j in 1..7 (C0.1, gamma = L by 7D.C4)"),
                                         ("ordered_pairs_checked", 49), ("ordered_pairs_holding", ok), ("holds", 1 if ok == 49 else 0)])

    xi = [{}] + [pv(a) for a in LABELS]
    Lxi = pL(Ls, xi)
    nsq = psum([pmul(pv(a), pv(a)) for a in LABELS])
    ident("Lxi_squared_plus_nsq_is_zero", pm_mul(Lxi, Lxi), pm_eye(8, pscale(nsq, -1)), "xi1..xi7",
          "L_xi L_xi = -nsq(xi) I_8 for imaginary xi, so L_xi L_xi + nsq(xi) I_8 = 0", count_terms(Lxi, Lxi))
    Sd = pm_zero(16, 16)
    LxiT = pm_T(Lxi)
    for r in range(8):
        for c in range(8):
            Sd[r][8 + c] = LxiT[r][c]
            Sd[8 + r][c] = Lxi[r][c]
    ident("doubled_selfadjoint_symbol", pm_mul(Sd, Sd), pm_eye(16, nsq), "xi1..xi7",
          "[[0, L_xi^T],[L_xi, 0]] is symmetric and squares to +nsq(xi) I_16, while L_xi squares to -nsq(xi) I_8 (C0.1 sign)",
          count_terms(Sd, Sd) + [("symmetric", 1 if compare(pm_T(Sd), Sd)[2] == 0 else 0)])
    u = [{}] + [pv(30 + a) for a in LABELS]
    v = [{}] + [pv(10 + a) for a in LABELS]
    Lu = pL(Ls, u)
    Lv = pL(Ls, v)
    uv = psum([pmul(pv(30 + a), pv(10 + a)) for a in LABELS])
    ident("clifford_polarized", pm_add(pm_mul(Lu, Lv), pm_mul(Lv, Lu)), pm_eye(8, pscale(uv, -2)), "u1..u7, v1..v7",
          "L_u L_v + L_v L_u = -2 <u, v> I_8 for imaginary u and v")
    ident("Lxi_skew", pm_T(Lxi), pm_scale(Lxi, -1), "xi1..xi7", "L_xi transpose = -L_xi for imaginary xi")
    x8 = [pv(0)] + [pv(a) for a in LABELS]
    xc = [pv(0)] + [pscale(pv(a), -1) for a in LABELS]
    Lx = pL(Ls, x8)
    Lxc = pL(Ls, xc)
    nsq8 = psum([pmul(pv(a), pv(a)) for a in range(8)])
    ident("composition_LxT_Lx", pm_mul(pm_T(Lx), Lx), pm_eye(8, nsq8), "x0..x7", "L_x transpose L_x = nsq(x) I_8 for every octonion x")
    ident("Lconjx_Lx", pm_mul(Lxc, Lx), pm_eye(8, nsq8), "x0..x7", "L_conj(x) L_x = nsq(x) I_8 for every octonion x")

    # ---- algebra generated, commutant, volume product ----
    E = Echelon()
    E.add(flat(I8))
    frontier = [I8]
    rounds = [1]
    while frontier:
        new = []
        for A in frontier:
            for a in LABELS:
                B = imul(A, Ls[a])
                if E.add(flat(B)):
                    new.append(B)
        frontier = new
        rounds.append(len(E.rows))
    alg_dim = len(E.rows)
    subset_mats = []
    for mask in range(128):
        P = I8
        for a in LABELS:
            if (mask >> (a - 1)) & 1:
                P = imul(P, Ls[a])
        subset_mats.append(flat(P))
    sub_rank = rank_int(subset_mats)
    mrow("algebra_generated_by_L_e1_to_L_e7", [("dimension", alg_dim), ("cumulative_dimension_by_word_length", ",".join(str(x) for x in rounds)),
                                               ("span_of_128_ordered_subset_products", sub_rank), ("full_matrix_algebra_dimension", 64),
                                               ("equals_M8R", 1 if alg_dim == 64 else 0),
                                               ("method", "closure of the span of words in L_e1..L_e7 under right multiplication by the generators; exact rational rank by integer elimination")])
    eqs = []
    for a in LABELS:
        for r in range(8):
            for s in range(8):
                row = [0] * 64
                for c in range(8):
                    row[r * 8 + c] += Ls[a][c][s]
                    row[c * 8 + s] -= Ls[a][r][c]
                eqs.append(row)
    rk, nb = nullspace_int(eqs, 64)
    is_scalar = 1 if (len(nb) == 1 and nb[0] == flat(I8)) else 0
    mrow("commutant_of_L_e1_to_L_e7", [("equations", len(eqs)), ("unknowns", 64), ("rank", rk), ("dimension", 64 - rk),
                                       ("basis_is_identity", is_scalar)])
    P = I8
    for a in LABELS:
        P = imul(P, Ls[a])
    sign = 1 if P == I8 else (-1 if P == iscal(I8, -1) else 0)
    central = sum(1 for a in LABELS if imul(P, Ls[a]) == imul(Ls[a], P))
    mrow("volume_product_L_e1_L_e2_to_L_e7", [("statement", "L_e1 L_e2 L_e3 L_e4 L_e5 L_e6 L_e7 = sign I_8"), ("sign", sign),
                                              ("is_plus_or_minus_identity", 1 if sign else 0), ("square_is_identity", 1 if imul(P, P) == I8 else 0),
                                              ("commutes_with_generators", str(central) + " of 7")])

    # ---- non-associativity ----
    hist = {}
    quat_kernel = 0
    for i in LABELS:
        for j in LABELS:
            if i == j:
                continue
            s, c = M[i][j]
            D = iadd(iscal(Ls[c], s), imul(Ls[i], Ls[j]), -1)
            r = rank_int(D)
            hist[r] = hist.get(r, 0) + 1
            kill = all(all(D[t][b] == 0 for t in range(8)) for b in (0, i, j, c))
            if kill and r == 4:
                quat_kernel += 1
    mrow("rank_L_eiej_minus_L_ei_L_ej", [("ordered_pairs", sum(hist.values())),
                                        ("rank_histogram", ";".join(str(k) + ":" + str(hist[k]) for k in sorted(hist))),
                                        ("kernel_is_quaternion_line_span_1_ei_ej_eiej", str(quat_kernel) + " of 42")])

    # ---- phi and the cross product ----
    PHI = [[[0] * 8 for _ in range(8)] for _ in range(8)]
    for i in LABELS:
        for j in LABELS:
            s, c = M[i][j]
            if c != 0:
                PHI[i][j][c] = s
    nz = 0
    antis = 0
    for i in LABELS:
        for j in LABELS:
            for k in LABELS:
                if PHI[i][j][k]:
                    nz += 1
                if PHI[i][j][k] == -PHI[j][i][k] and PHI[i][j][k] == -PHI[i][k][j]:
                    antis += 1
    phi_form = {}
    terms = []
    for (a, b, c) in itertools.combinations(LABELS, 3):
        if PHI[a][b][c]:
            phi_form[(a, b, c)] = pc(PHI[a][b][c])
            terms.append(("+" if PHI[a][b][c] > 0 else "-") + "e" + str(a) + str(b) + str(c))
    mrow("phi_3form", [("definition", "phi_ijk = <e_i e_j, e_k> on Im O (7D.C2)"), ("nonzero_ordered_components", nz),
                       ("antisymmetric_ordered_components", str(antis) + " of 343"), ("independent_terms", len(phi_form)),
                       ("phi", "".join(terms)), ("phi_e123", PHI[1][2][3])])
    X = [[psum([pscale(pv(i), PHI[i][j][k]) for i in LABELS if PHI[i][j][k]]) for j in LABELS] for k in LABELS]
    Lim = [[Lxi[k][j] for j in LABELS] for k in LABELS]
    ident("cross_matrix_equals_Im_block_of_Lxi", X, Lim, "xi1..xi7", "X(xi) v = xi x v = Im(xi v), (xi x v)_k = sum_ij phi_ijk xi_i v_j")
    xxT = [[pmul(pv(a), pv(b)) for b in LABELS] for a in LABELS]
    ident("cross_matrix_square", pm_mul(X, X), pm_add(pm_eye(7, pscale(nsq, -1)), xxT), "xi1..xi7",
          "X(xi) X(xi) = -nsq(xi) I_7 + xi xi^T")
    ident("cross_matrix_kills_xi", pm_mul(X, col([pv(a) for a in LABELS])), pm_zero(7, 1), "xi1..xi7", "X(xi) xi = 0",
          count_terms(X, col([pv(a) for a in LABELS])))
    ident("cross_matrix_skew", pm_T(X), pm_scale(X, -1), "xi1..xi7", "X(xi) transpose = -X(xi)")
    mrow("cross_matrix_rank_certificate", [("rank_for_every_nonzero_xi", 6),
                                           ("certified_by", "X xi = 0 gives rank at most 6; X X = -nsq I + xi xi^T gives X X w = -nsq w on xi-perp, so X is injective there and rank is at least 6"),
                                           ("identities_used", "cross_matrix_square cross_matrix_kills_xi")])

    # ---- Hodge star, T7, orientation ----
    def T_matrix(orient):
        T = izero(21, 21)
        for b, pb in enumerate(PAIRS):
            F = fhodge(fwedge(phi_form, {pb: pc(1)}), orient)
            for K, cc in F.items():
                T[PIDX[K]][b] = cc.get((), 0)
        return T
    I21 = ieye(21)
    T0 = T_matrix(1)
    qa = imul(iadd(T0, I21, -2), iadd(T0, I21))
    qb = imul(iadd(T0, I21, 2), iadd(T0, I21, -1))
    r0 = {"T+I": rank_int(iadd(T0, I21)), "T-2I": rank_int(iadd(T0, I21, -2)), "T-I": rank_int(iadd(T0, I21, -1)), "T+2I": rank_int(iadd(T0, I21, 2))}
    qa0 = 1 if qa == izero(21, 21) else 0
    qb0 = 1 if qb == izero(21, 21) else 0
    if qa0 and r0["T-2I"] == 14 and r0["T+I"] == 7:
        orient = 1
    elif qb0 and r0["T+2I"] == 14 and r0["T-I"] == 7:
        orient = -1
    else:
        orient = 0
    mrow("T7_with_vol_e1_to_e7", [("T", "*(phi wedge .) on Lambda2 with * from vol = e1 wedge ... wedge e7"),
                                  ("T_minus_2I_times_T_plus_I_is_zero", qa0), ("T_plus_2I_times_T_minus_I_is_zero", qb0),
                                  ("rank_T_plus_I", r0["T+I"]), ("rank_T_minus_2I", r0["T-2I"]), ("rank_T_minus_I", r0["T-I"]), ("rank_T_plus_2I", r0["T+2I"]),
                                  ("eigenvalues_for_e1_to_e7", ("+2 on 7, -1 on 14" if orient == 1 else ("-2 on 7, +1 on 14" if orient == -1 else "neither split"))),
                                  ("trace", sum(T0[i][i] for i in range(21)))])
    mrow("orientation_7D", [("rule", "7D.C3: orient so T7 has +2 on a 7-dim and -1 on a 14-dim subspace"),
                            ("orientation_sign_relative_to_e1_to_e7", orient), ("reversed", 1 if orient == -1 else 0),
                            ("note", "vol = e1 wedge ... wedge e7 kept" if orient == 1 else "vol = -(e1 wedge ... wedge e7) used, the measurement reverses it")])
    T = T_matrix(orient)
    TpI = iadd(T, I21)
    Tm2 = iadd(T, I21, -2)
    iident("T7_minimal_polynomial", imul(Tm2, TpI), izero(21, 21), "(T7 - 2I)(T7 + I) = 0 in the chosen orientation")
    iident("T7_symmetric", iT(T), T, "T7 transpose = T7")
    mrow("T7_eigen_split", [("eigenvalue_plus2_multiplicity", 21 - rank_int(Tm2)), ("eigenvalue_minus1_multiplicity", 21 - rank_int(TpI)),
                            ("trace", sum(T[i][i] for i in range(21))),
                            ("certified_by", "(T-2I)(T+I) = 0 with distinct roots, multiplicities are nullities")])
    psi = fhodge(phi_form, orient)
    pp = fwedge(phi_form, psi)
    mrow("psi_and_phi_wedge_psi", [("psi", "*phi in the chosen orientation"), ("psi_terms", len(psi)),
                                   ("phi_wedge_psi_coefficient_on_e1_to_e7", pp.get(VOL, {}).get((), 0)),
                                   ("phi_wedge_psi_in_chosen_vol", orient * pp.get(VOL, {}).get((), 0)), ("expected_nsq_phi", 7)])
    vvals = [pv(10 + a) for a in LABELS]
    iota_v_phi = fclean({(i, j): psum([pscale(pv(10 + k), PHI[k][i][j]) for k in LABELS if PHI[k][i][j]]) for (i, j) in PAIRS})
    W6 = fwedge(psi, iota_v_phi)
    starv = fhodge(one_form(vvals), orient)
    keys6 = sorted(set(W6) | set(starv))
    found = 0
    for kk in (3, -3, 1, -1):
        A = [[W6.get(K, {})] for K in keys6]
        B = [[pscale(starv.get(K, {}), kk)] for K in keys6]
        if compare(A, B)[2] == 0:
            found = kk
            break
    mrow("psi_wedge_iota_v_phi", [("statement", "psi wedge iota_v phi = k * v as a polynomial identity in v"), ("k", found),
                                  ("variables", "v1..v7")])
    xiv = fwedge(one_form([pv(a) for a in LABELS]), one_form(vvals))
    Y = fhodge(fwedge(xiv, psi), orient)
    Xv = pm_mul(X, col(vvals))
    found1 = 0
    for kk in (1, -1, 3, -3):
        A = [[Y.get((a,), {})] for a in LABELS]
        B = [[pscale(Xv[a - 1][0], kk)] for a in LABELS]
        if compare(A, B)[2] == 0:
            found1 = kk
            break
    mrow("star_xi_v_psi_vs_cross", [("statement", "*(xi wedge v wedge psi) = k (xi x v)"), ("k", found1), ("variables", "xi1..xi7, v1..v7")])

    # ---- Psi_7, Iota, projections ----
    Psi7 = [[PHI[i][j][k] for (i, j) in PAIRS] for k in LABELS]
    Iota = [[PHI[k][i][j] for k in LABELS] for (i, j) in PAIRS]
    cliff_ok = 0
    for b, (i, j) in enumerate(PAIRS):
        colv = [imul(Ls[i], Ls[j])[t][0] for t in range(8)]
        if colv == [0] + [Psi7[k][b] for k in range(7)]:
            cliff_ok += 1
    mrow("Psi7_is_clifford_action_on_eta", [("statement", "Psi(e_ij) = gamma(e_i) gamma(e_j) eta with eta = 1 (C0.6) equals the contraction with phi (7D.C5)"),
                                           ("pairs_holding", str(cliff_ok) + " of 21")])
    iident("TplusI_equals_Iota_Psi7", TpI, imul(Iota, Psi7), "(T7 + I) omega = iota_(Psi7 omega) phi on all of Lambda2")
    iident("Psi7_Iota_is_3I", imul(Psi7, Iota), iscal(ieye(7), 3), "Psi7(iota_u phi) = 3 u")
    iident("Psi7_kills_Lambda2_14", imul(Psi7, Tm2), izero(7, 21), "Psi7 (T7 - 2I) = 0, so the -1 eigenspace lies in ker Psi7")
    iident("Psi7_TplusI_is_3Psi7", imul(Psi7, TpI), iscal(Psi7, 3), "Psi7 (T7 + I) = 3 Psi7")
    iident("TplusI_squared", imul(TpI, TpI), iscal(TpI, 3), "(T7 + I)^2 = 3 (T7 + I), so (T7 + I) times one third is a projection")
    kr, kb = nullspace_int(Psi7, 21)
    mrow("Psi7_rank_and_kernel", [("rank", kr), ("kernel_dimension", len(kb)), ("rank_Iota", rank_int(Iota)),
                                  ("kernel_equals_minus1_eigenspace", 1 if (len(kb) == 21 - rank_int(TpI) and rank_int(Tm2) == len(kb)) else 0)])
    psiwedge = []
    for b, pb in enumerate(PAIRS):
        F = fwedge({pb: pc(1)}, psi)
        psiwedge.append([F.get(tuple(x for x in LABELS if x != m), {}).get((), 0) for m in LABELS])
    wedge_psi_mat = iT(psiwedge)
    iident("wedge_psi_kills_Lambda2_14", imul(wedge_psi_mat, Tm2), izero(7, 21), "(T7 - 2I) omega wedge psi = 0")
    mrow("wedge_psi_rank", [("rank_of_omega_to_omega_wedge_psi", rank_int(wedge_psi_mat))])

    # ---- (T7 + I)(xi wedge v) = iota_(xi x v) phi ----
    xivec = two_vec(xiv)
    lhs = pm_imul(TpI, col(xivec))
    uvec = [Xv[a][0] for a in range(7)]
    rhs = [[psum([pscale(uvec[k - 1], PHI[k][i][j]) for k in LABELS if PHI[k][i][j]])] for (i, j) in PAIRS]
    ident("TplusI_xi_wedge_v_equals_iota_cross_phi", lhs, rhs, "xi1..xi7, v1..v7", "(T7 + I)(xi wedge v) = iota_(xi x v) phi")
    nsq_v = psum([pmul(pv(10 + a), pv(10 + a)) for a in LABELS])
    dotxv = psum([pmul(pv(a), pv(10 + a)) for a in LABELS])
    ident("forms_norm_lagrange", [[psum([pmul(x, x) for x in xivec])]], [[padd(pmul(nsq, nsq_v), pmul(dotxv, dotxv), -1)]], "xi1..xi7, v1..v7",
          "nsq(xi wedge v) = nsq(xi) nsq(v) - <xi, v>^2 (C0.9)")
    ident("R1_cleared", pm_imul(Psi7, lhs), pm_scale(col(uvec), 3), "xi1..xi7, v1..v7",
          "3 Psi7(pi_7(xi wedge v)) = Psi7((T7 + I)(xi wedge v)) = 3 (xi x v); with pi_7 (xi wedge v) = c iota_(xi x v) phi")
    mrow("R1_constant_c", [("c_numerator", 1), ("c_denominator", 3),
                           ("statement", "pi_7(xi wedge v) = c iota_(xi x v) phi with c = one third; cleared: 3 pi_7(xi wedge v) = (T7 + I)(xi wedge v) = iota_(xi x v) phi"),
                           ("from", "TplusI_xi_wedge_v_equals_iota_cross_phi and Psi7_Iota_is_3I and TplusI_squared")])

    # ---- G2 symbol sequence 1, 7, 7, 1 ----
    xi_form = one_form([pv(a) for a in LABELS])
    S0 = col([pv(a) for a in LABELS])
    S1 = [[{} for _ in LABELS] for _ in LABELS]
    for j in LABELS:
        wv = two_vec(fwedge(xi_form, {(j,): pc(1)}))
        for k in range(7):
            S1[k][j - 1] = psum([pscale(wv[b], Psi7[k][b]) for b in range(21) if Psi7[k][b] and wv[b]])
    S2x3 = [[]]
    for k in LABELS:
        iota_ek = fclean({(i, j): pc(PHI[k][i][j]) for (i, j) in PAIRS})
        S2x3[0].append(finner(fwedge(xi_form, iota_ek), phi_form))
    div3 = all(c % 3 == 0 for p in S2x3[0] for c in p.values())
    S2 = [[{m: exact_div(c, 3) for m, c in p.items()} for p in S2x3[0]]]
    mrow("G2_sequence_maps", [("S0", "Lambda0 to Lambda1, f to f xi (symbol of d)"),
                              ("S1", "Lambda1 to Lambda2_7 in w coordinates, v to Psi7(xi wedge v) (symbol of pi_7 d)"),
                              ("S2", "Lambda2_7 in w coordinates to Lambda3_1 in t coordinates, w to t = <xi wedge pi_7 omega, phi> with pi_7 omega = one third iota_w phi"),
                              ("three_S2_coefficients_divisible_by_3", 1 if div3 else 0)])
    ident("G2_S1_equals_cross_matrix", S1, X, "xi1..xi7", "S1(xi) = X(xi)")
    ident("G2_S2_equals_xi_transpose", S2, pm_T(S0), "xi1..xi7", "S2(xi) w = <xi, w> (from 3 S2 = <xi wedge iota_w phi, phi>)")
    tz = [[]]
    wterms = 0
    for b in range(21):
        om = form_from_two_vec([pc(Tm2[a][b]) for a in range(21)])
        w3 = fwedge(xi_form, om)
        wterms += sum(len(p) for p in w3.values())
        tz[0].append(finner(w3, phi_form))
    ident("G2_t_sees_only_Lambda2_7", tz, pm_zero(1, 21), "xi1..xi7", "<xi wedge (T7 - 2I) e_b, phi> = 0 for all 21 b",
          [("three_form_coefficient_terms_paired_with_phi", wterms)])
    h1 = ident("G2_S1_S0_zero", pm_mul(S1, S0), pm_zero(7, 1), "xi1..xi7", "S1 S0 = 0", count_terms(S1, S0))
    h2 = ident("G2_S2_S1_zero", pm_mul(S2, S1), pm_zero(1, 7), "xi1..xi7", "S2 S1 = 0", count_terms(S2, S1))
    h3 = ident("G2_S0T_S0", pm_mul(pm_T(S0), S0), [[nsq]], "xi1..xi7", "S0^T S0 = nsq(xi)")
    h4 = ident("G2_S2_S2T", pm_mul(S2, pm_T(S2)), [[nsq]], "xi1..xi7", "S2 S2^T = nsq(xi)")
    h5 = ident("G2_laplacian_Lambda1", pm_add(pm_mul(S0, pm_T(S0)), pm_mul(pm_T(S1), S1)), pm_eye(7, nsq), "xi1..xi7",
               "S0 S0^T + S1^T S1 = nsq(xi) I_7")
    h6 = ident("G2_laplacian_Lambda2_7", pm_add(pm_mul(S1, pm_T(S1)), pm_mul(pm_T(S2), S2)), pm_eye(7, nsq), "xi1..xi7",
               "S1 S1^T + S2^T S2 = nsq(xi) I_7")
    allh = h1 * h2 * h3 * h4 * h5 * h6
    mrow("G2_symbol_sequence", [("dims", "1,7,7,1"), ("ranks_for_every_nonzero_xi", "1,6,1" if allh else "NOT_CERTIFIED"),
                                ("exact_for_every_nonzero_xi", 1 if allh else 0), ("euler_count", 1 - 7 + 7 - 1),
                                ("certified_by", "S0 injective by S0^T S0 = nsq; S2 onto by S2 S2^T = nsq; S1 S0 = 0 and S0 S0^T + S1^T S1 = nsq I give ker S1 = im S0, rank 6; S2 S1 = 0 and S1 S1^T + S2^T S2 = nsq I give ker S2 = im S1"),
                                ("normalization", "w coordinates on Lambda2_7 and t coordinates on Lambda3_1 as 7D.C5")])
    samp = {1: 1, 2: -1, 3: 2, 4: 0, 5: 3, 6: 1, 7: -2}
    def ev(A):
        return [[peval(x, samp) for x in r] for r in A]
    mrow("sampled_covector_control", [("xi", "1,-1,2,0,3,1,-2"), ("rank_S0", rank_int(iT(ev(S0)))), ("rank_S1", rank_int(ev(S1))),
                                      ("rank_S2", rank_int(ev(S2))), ("role", "control only (C0.3), not a certificate")])

    # ---- rolled-up symbol ----
    R = pm_zero(8, 8)
    for j in LABELS:
        R[0][j] = pscale(pv(j), -1)
        R[j][0] = pv(j)
        for k in LABELS:
            R[k][j] = X[k - 1][j - 1]
    ident("rolled_up_R_equals_Lxi", R, Lxi, "xi1..xi7", "R(xi)(f, v) = (-<xi, v>, f xi + xi x v) equals L_xi under (f, v) to f 1 + v (7D.C4)")
    R2 = pm_zero(8, 8)
    for j in range(7):
        R2[0][j + 1] = pscale(S2[0][j], -1)
        R2[j + 1][0] = S0[j][0]
        for k in range(7):
            R2[k + 1][j + 1] = S1[k][j]
    ident("R2_cross_product_coordinates", R2, Lxi, "xi1..xi7", "[[0, -S2],[S0, S1]] in w coordinates equals L_xi")
    ident("R_transpose_R", pm_mul(pm_T(R), R), pm_eye(8, nsq), "xi1..xi7", "R(xi)^T R(xi) = nsq(xi) I_8")
    ident("R_squared", pm_mul(R, R), pm_eye(8, pscale(nsq, -1)), "xi1..xi7", "R(xi)^2 = -nsq(xi) I_8")

    # ---- R4 associator ----
    e1, e2, e3 = unit(1), unit(2), unit(3)
    p12 = omul(M, e1, e2)
    p23 = omul(M, e2, e3)
    left = omul(M, p12, e3)
    right = omul(M, e1, p23)
    assoc = [left[t] - right[t] for t in range(8)]
    mrow("R4_associator_e1_e2_e3", [("e1e2", vec_str(p12)), ("e2e3", vec_str(p23)), ("(e1e2)e3", vec_str(left)), ("e1(e2e3)", vec_str(right)),
                                    ("associator", vec_str(assoc)), ("phi_e1_e2_e3", PHI[1][2][3])])
    zero_t = 0
    two_t = 0
    alt_t = 0
    for (i, j, k) in itertools.permutations(LABELS, 3):
        a1 = [omul(M, omul(M, unit(i), unit(j)), unit(k))[t] - omul(M, unit(i), omul(M, unit(j), unit(k)))[t] for t in range(8)]
        a2 = [omul(M, omul(M, unit(j), unit(i)), unit(k))[t] - omul(M, unit(j), omul(M, unit(i), unit(k)))[t] for t in range(8)]
        if not any(a1):
            zero_t += 1
        elif sum(1 for t in a1 if t) == 1 and sum(abs(t) for t in a1) == 2:
            two_t += 1
        if a1 == [-t for t in a2]:
            alt_t += 1
    mrow("associator_census", [("ordered_distinct_triples", 210), ("zero", zero_t), ("plus_or_minus_2_e_m", two_t),
                               ("antisymmetric_in_first_two", str(alt_t) + " of 210")])

    # ---- der(O) ----
    deq = []
    for a in range(8):
        La = Ls[a]
        for r in range(8):
            for s in range(8):
                row = [0] * 64
                for c in range(8):
                    row[r * 8 + c] += La[c][s]
                    row[c * 8 + s] -= La[r][c]
                    row[c * 8 + a] -= Ls[c][r][s]
                deq.append(row)
    drank, dbasis = nullspace_int(deq, 64)
    Ds = [[b[r * 8:(r + 1) * 8] for r in range(8)] for b in dbasis]
    kills1 = sum(1 for D in Ds if all(D[r][0] == 0 for r in range(8)) and all(D[0][c] == 0 for c in range(8)))
    skew = sum(1 for D in Ds if iT(D) == iscal(D, -1))
    pair_ok = 0
    prods = {(a, b): omul(M, unit(a), unit(b)) for a in range(8) for b in range(8)}
    for D in Ds:
        good = 1
        for a in range(8):
            for b in range(8):
                ea = unit(a)
                eb = unit(b)
                pab = prods[(a, b)]
                lhs_v = [sum(D[r][t] * pab[t] for t in range(8)) for r in range(8)]
                Dea = [D[r][a] for r in range(8)]
                Deb = [D[r][b] for r in range(8)]
                x1 = omul(M, Dea, eb)
                x2 = omul(M, ea, Deb)
                if lhs_v != [x1[t] + x2[t] for t in range(8)]:
                    good = 0
        pair_ok += good
    mrow("derivations_der_O", [("equations", len(deq)), ("unknowns", 64), ("rank", drank), ("dimension", len(Ds)),
                               ("kill_1_and_preserve_Im_O", str(kills1) + " of " + str(len(Ds))), ("skew", str(skew) + " of " + str(len(Ds))),
                               ("D_of_product_rule_on_all_64_basis_pairs", str(pair_ok) + " of " + str(len(Ds)))])
    for n, D in enumerate(Ds):
        ent = ";".join(str(r) + "," + str(c) + ":" + str(D[r][c]) for r in range(8) for c in range(8) if D[r][c] and r < c)
        mrow("derivation_basis_" + two(n + 1), [("upper_entries_row_col_value", ent), ("lower", "minus the transpose")])
    lie = Echelon()
    for D in Ds:
        lie.add(flat(D))
    base_rank = len(lie.rows)
    grew = 0
    for p in range(len(Ds)):
        for q in range(p + 1, len(Ds)):
            Cm = iadd(imul(Ds[p], Ds[q]), imul(Ds[q], Ds[p]), -1)
            grew += lie.add(flat(Cm))
    mrow("der_O_lie_closure", [("commutators_checked", exact_div(len(Ds) * (len(Ds) - 1), 2)), ("span_rank_before", base_rank),
                               ("span_rank_after", len(lie.rows)), ("closed", 1 if grew == 0 else 0)])
    xi_im = [{}] + [pv(a) for a in LABELS]
    total_c = 0
    total_bad = 0
    hold_im = 0
    hold_full = 0
    for D in Ds:
        Dp = pm_int(D)
        for (vecx, Lv_, tag) in ((xi_im, Lxi, "im"), (x8, Lx, "full")):
            com = pm_add(pm_mul(Dp, Lv_), pm_mul(Lv_, Dp), -1)
            Dx = [psum([pscale(vecx[b], D[a][b]) for b in range(8) if D[a][b] and vecx[b]]) for a in range(8)]
            e, c, bad, mx = compare(com, pL(Ls, Dx))
            total_c += c
            total_bad += bad
            if bad == 0:
                if tag == "im":
                    hold_im += 1
                else:
                    hold_full += 1
    mrow("derivation_equivariance_D_Lxi", [("statement", "[D, L_xi] = L_(D xi) for each basis derivation D"),
                                           ("holding_imaginary_xi", str(hold_im) + " of " + str(len(Ds))),
                                           ("holding_full_octonion_x", str(hold_full) + " of " + str(len(Ds))),
                                           ("coefficients_compared", total_c), ("coefficients_unequal", total_bad),
                                           ("method", "integer polynomial identity, coefficient by coefficient, in xi1..xi7 and in x0..x7")])
    Rk = []
    for kv in kb:
        Rm = izero(8, 8)
        for b, (i, j) in enumerate(PAIRS):
            if kv[b]:
                Rm[j][i] += kv[b]
                Rm[i][j] -= kv[b]
        Rk.append(Rm)
    rk_ok = 0
    for Rm in Rk:
        fv = flat(Rm)
        if all(sum(rw[t] * fv[t] for t in range(64)) == 0 for rw in deq):
            rk_ok += 1
    both = Echelon()
    for D in Ds:
        both.add(flat(D))
    for Rm in Rk:
        both.add(flat(Rm))
    mrow("der_O_equals_R_of_ker_Psi7", [("ker_Psi7_rotations_that_are_derivations", str(rk_ok) + " of " + str(len(Rk))),
                                        ("rank_der", len(Ds)), ("rank_R_ker", rank_int([flat(Rm) for Rm in Rk])),
                                        ("rank_union", len(both.rows)), ("equal", 1 if len(both.rows) == len(Ds) == len(Rk) else 0)])
    om = {}
    cmat = pm_zero(8, 8)
    for (i, j) in PAIRS:
        om[(i, j)] = pv(100 + 10 * i + j)
        cmat = pm_add(cmat, pm_scale(pm_mul(pm_int(imul(Ls[i], Ls[j])), pm_eye(8, om[(i, j)])), 1))
    Rv = [{} for _ in range(8)]
    for (i, j) in PAIRS:
        Rv[j] = padd(Rv[j], pmul(om[(i, j)], pv(10 + i)))
        Rv[i] = padd(Rv[i], pmul(om[(i, j)], pv(10 + j)), -1)
    ident("C0_5_rotation_bracket", pm_add(pm_mul(cmat, Lv), pm_mul(Lv, cmat), -1), pm_scale(pL(Ls, Rv), 2), "omega_ij (21), v1..v7",
          "[c(omega), gamma(v)] = 2 gamma(R_omega v) with c(omega) = sum omega_ij L_ei L_ej")

    # ---- signed permutations preserving the table ----
    auts = []
    nonauts = []
    examined = 0
    for perm in itertools.permutations(LABELS):
        for signs in itertools.product((1, -1), repeat=7):
            examined += 1
            good = True
            for i in LABELS:
                si = signs[i - 1]
                pi_ = perm[i - 1]
                for j in LABELS:
                    if i == j:
                        continue
                    s, k = M[i][j]
                    s2, k2 = M[pi_][perm[j - 1]]
                    if k2 != perm[k - 1] or si * signs[j - 1] * s2 != s * signs[k - 1]:
                        good = False
                        break
                if not good:
                    break
            if good:
                auts.append(tuple(signs[i] * perm[i] for i in range(7)))
            elif len(nonauts) < 1344:
                nonauts.append(tuple(signs[i] * perm[i] for i in range(7)))
    sq_ok = sum(1 for g in auts if all(M[abs(t)][abs(t)] == (-1, 0) for t in g))

    def conj_ok(g):
        p = [0] + [abs(t) for t in g]
        s = [1] + [(1 if t > 0 else -1) for t in g]
        for a in range(8):
            La = Ls[a]
            Lga = Ls[p[a]]
            for r in range(8):
                for c in range(8):
                    if s[r] * s[c] * La[r][c] != s[a] * Lga[p[r]][p[c]]:
                        return False
        return True
    conj_pass = sum(1 for g in auts if conj_ok(g))
    conj_fail_ctrl = sum(1 for g in nonauts if not conj_ok(g))
    aset = set(auts)

    def comp(g, h):
        out = []
        for t in h:
            gt = g[abs(t) - 1]
            out.append(gt if t > 0 else -gt)
        return tuple(out)
    closed = 0
    for g in auts:
        for h in auts:
            if comp(g, h) in aset:
                closed += 1
    ident_t = tuple(LABELS)
    inv_ok = 0
    for g in auts:
        inv = [0] * 7
        for i in range(7):
            t = g[i]
            inv[abs(t) - 1] = (i + 1) if t > 0 else -(i + 1)
        if tuple(inv) in aset and comp(g, tuple(inv)) == ident_t:
            inv_ok += 1
    phi_ok = 0
    for g in auts:
        good = True
        for (a, b, c) in itertools.combinations(LABELS, 3):
            ga, gb, gc = g[a - 1], g[b - 1], g[c - 1]
            val = PHI[abs(ga)][abs(gb)][abs(gc)] * (1 if ga > 0 else -1) * (1 if gb > 0 else -1) * (1 if gc > 0 else -1)
            if val != PHI[a][b][c]:
                good = False
                break
        if good:
            phi_ok += 1
    orbit = len(set(g[0] for g in auts))
    stab = sum(1 for g in auts if g[0] == 1)
    sign_only = sum(1 for g in auts if tuple(abs(t) for t in g) == ident_t)
    perms = len(set(tuple(abs(t) for t in g) for g in auts))
    mrow("table_preserving_signed_permutations", [("candidates_examined", examined), ("preserving_all_42_offdiagonal_products", len(auts)),
                                                  ("squares_preserved", str(sq_ok) + " of " + str(len(auts))),
                                                  ("closed_under_composition", str(closed) + " of " + str(len(auts) * len(auts))),
                                                  ("inverses_in_set", str(inv_ok) + " of " + str(len(auts))),
                                                  ("preserve_phi", str(phi_ok) + " of " + str(len(auts))),
                                                  ("orbit_of_e1", orbit), ("stabilizer_of_e1", stab), ("sign_only", sign_only),
                                                  ("distinct_unsigned_permutations", perms)])
    mrow("signed_permutation_equivariance", [("statement", "g L_xi g^-1 = L_(g xi), g^-1 = g^T, for every table-preserving signed permutation g"),
                                             ("holding", str(conj_pass) + " of " + str(len(auts))),
                                             ("method", "the identity is linear in xi, so its xi_a coefficient is g L_ea g^T = L_(g ea); all 8 coefficients, all 64 entries, exact integers")])

    # ---- dims ----
    mrow("dims_7D", [("dim_S", 8), ("dim_Lambda0", 1), ("dim_Lambda1", 7), ("dim_Lambda2", len(PAIRS)),
                     ("dim_Lambda2_7", 21 - rank_int(Tm2)), ("dim_Lambda2_14", 21 - rank_int(TpI)),
                     ("dim_Lambda3", len(list(itertools.combinations(LABELS, 3)))),
                     ("dim_Lambda3_1", rank_int([[PHI[a][b][c] for (a, b, c) in itertools.combinations(LABELS, 3)]])),
                     ("dim_der_O", len(Ds)), ("dim_ker_Psi7", len(kb)), ("dim_algebra_generated", alg_dim),
                     ("order_signed_automorphisms", len(auts))])
    mrow("dims_7D_note", [("note", "Lambda2_7 = nullity of T7 - 2I and Lambda2_14 = nullity of T7 + I in the chosen orientation; Lambda3_1 = span(phi), its dimension the rank of the coefficient row of phi")])

    # ---- controls that must fail ----
    Mf, _ = build_mult(lines, flip=(0,))
    Lsf = left_mats(Mf)
    Lxif = pL(Lsf, xi)
    f1, kv1 = fails(pm_add(pm_mul(Lxif, Lxif), pm_eye(8, nsq)), pm_zero(8, 8))
    okf = 0
    for i in LABELS:
        for j in LABELS:
            if iadd(imul(Lsf[i], Lsf[j]), imul(Lsf[j], Lsf[i])) == iscal(I8, -2 if i == j else 0):
                okf += 1
    krow("flipped_line_124_table_Lxi_square", f1, [("statement", "L_xi^2 + nsq I = 0 with line (1,2,4) sign flipped"),
                                                   ("clifford_basis_pairs_holding", str(okf) + " of 49")] + kv1)
    f2, kv2 = fails(pm_mul(X, X), pm_eye(7, pscale(nsq, -1)))
    krow("cross_square_without_xi_xiT", f2, [("statement", "X(xi)^2 = -nsq I_7")] + kv2)
    f3, kv3 = fails(lhs, [[pscale(r_[0], 3)] for r_ in rhs])
    krow("c_equal_1_instead_of_one_third", f3, [("statement", "(T7 + I)(xi wedge v) = 3 iota_(xi x v) phi")] + kv3)
    Trev = iscal(T, -1)
    f4, kv4 = fails(pm_imul(iadd(Trev, I21), col(xivec)), rhs)
    krow("reversed_orientation", f4, [("statement", "(-T7 + I)(xi wedge v) = iota_(xi x v) phi")] + kv4)
    Rbad = [list(r_) for r_ in R]
    Rbad[0] = [{}] + [pv(j) for j in LABELS]
    f5, kv5 = fails(Rbad, Lxi)
    krow("rolled_up_sign_flipped", f5, [("statement", "(+<xi, v>, f xi + xi x v) equals L_xi")] + kv5)
    s12, c12 = M[1][2]
    f6 = iscal(Ls[c12], s12) != imul(Ls[1], Ls[2])
    krow("associativity_L_e1e2", f6, [("statement", "L_(e1 e2) = L_e1 L_e2")])
    R12 = izero(8, 8)
    R12[2][1] = 1
    R12[1][2] = -1
    com12 = pm_add(pm_mul(pm_int(R12), Lxi), pm_mul(Lxi, pm_int(R12)), -1)
    Dx12 = [psum([pscale(xi_im[b], R12[a][b]) for b in range(8) if R12[a][b] and xi_im[b]]) for a in range(8)]
    f7, kv7 = fails(com12, pL(Ls, Dx12))
    krow("rotation_R12_not_derivation", f7, [("statement", "[R_12, L_xi] = L_(R_12 xi)"), ("Psi7_e12", ",".join(str(Psi7[k][PIDX[(1, 2)]]) for k in range(7)))] + kv7)
    krow("non_automorphism_signed_permutations", conj_fail_ctrl == len(nonauts) and len(nonauts) == 1344,
         [("statement", "g L_xi g^T = L_(g xi) for the first 1344 non-preserving signed permutations in enumeration order"),
          ("failing", str(conj_fail_ctrl) + " of " + str(len(nonauts)))])
    tzbad = [[]]
    for b in range(21):
        om2 = form_from_two_vec([pc(TpI[a][b]) for a in range(21)])
        tzbad[0].append(finner(fwedge(xi_form, om2), phi_form))
    f9, kv9 = fails(tzbad, pm_zero(1, 21))
    krow("t_on_Lambda2_7_is_not_zero", f9, [("statement", "<xi wedge (T7 + I) e_b, phi> = 0 for all b")] + kv9)
    f10, kv10 = fails(pm_add(pm_mul(S0, pm_T(S0)), pm_mul(pm_T(S1), S1)), pm_eye(7, pscale(nsq, 3)))
    krow("laplacian_wrong_constant", f10, [("statement", "S0 S0^T + S1^T S1 = 3 nsq(xi) I_7")] + kv10)
    emit("MEASURED_END", [("m_rows", MCOUNT[0]), ("control_rows", KCOUNT[0]),
                          ("controls_ok", sum(1 for r_ in ROWS if r_.startswith("CONTROL|") and "|control_ok=1|" in r_))])
    return {"M": M, "Ls": Ls, "PHI": PHI, "TpI": TpI, "Psi7": Psi7, "Iota": Iota, "X": X, "nsq": nsq,
            "phi_form": phi_form, "psi": psi, "kb": kb, "Rk": Rk, "xivec": xivec, "xi_form": xi_form}


def sha_hex(b):
    return hashlib.sha256(b).hexdigest()


def main():
    if len(sys.argv) < 2:
        print("usage: python3 -I -B measure_7d.py OUTDIR [PREDICTION_HBP ADDENDUM_HBP]")
        sys.exit(2)
    outdir = sys.argv[1]
    ctx = measure()
    block = ("\n".join(ROWS) + "\n").encode("ascii")
    mrows_list = list(ROWS)
    measured_rows = len(ROWS)
    bsha = sha_hex(block)
    del ROWS[:]
    emit("MEASURED_BLOCK", [("rows", measured_rows), ("sha256", bsha), ("written_before_reading_PRED_ADD_DIMS", 1)])
    if len(sys.argv) >= 4:
        compare_frozen(sys.argv[2], sys.argv[3], bsha, mrows_list, ctx)
    else:
        emit("PHASE", [("compare", "NOT_RUN"), ("reason", "measure-only run; the frozen PRED, ADD and DIMS rows are unread")])
    with open(os.path.abspath(__file__), "rb") as fh:
        ssha = sha_hex(fh.read())
    emit("RUN", [("script", "measure_7d.py"), ("script_sha256", ssha), ("run", "python3 -I -B measure_7d.py OUTDIR [PREDICTION_HBP ADDENDUM_HBP]"),
                 ("python", str(sys.version_info[0]) + "." + str(sys.version_info[1]) + "." + str(sys.version_info[2])), ("float_used", 0)])
    total = measured_rows + len(ROWS) + 1
    emit("END", [("rows", total)])
    data = block + ("\n".join(ROWS) + "\n").encode("ascii")
    path = os.path.join(outdir, "MEASURE-7D.hbp")
    with open(path, "wb") as fh:
        fh.write(data)
    with open(path + ".sha256", "wb") as fh:
        fh.write((sha_hex(data) + "  MEASURE-7D.hbp\n").encode("ascii"))
    print("MEASURE_7D|rows=" + str(total) + "|measured_block_sha256=" + bsha + "|file_sha256=" + sha_hex(data) + "|json=0")


PRECOMPARE_SEAL = "b716de4b034348b33370b25d2c868c7ee699b340391f5959f512df7ff274aa91"
PRED_SHA = "c6e86dd46662ede6a0bb82a8a8dda133a6fc0c20de54a6d947575e37169b2d3e"
ADD_SHA = "ad9ae33ef3e4ac1befb3528ebfa6faeaa10274cb044a838afdece4859ccc08c0"
QR = {}
QCOUNT = [0]


def qrow(name, kv):
    QCOUNT[0] += 1
    fields = [("id", "7D.Q" + two(QCOUNT[0])), ("name", name)] + kv + [("claim", "MEASURED_POSTREAD")]
    QR[name] = dict((k, str(v)) for k, v in fields)
    emit("POSTREAD", fields)


def q_ident(name, A, B, variables, statement):
    e, c, bad, mx = compare(A, B)
    qrow(name, [("statement", statement), ("holds", 1 if bad == 0 else 0), ("entries", e), ("coefficients_compared", c),
                ("coefficients_unequal", bad), ("max_abs_residual", mx), ("variables", variables),
                ("method", "integer polynomial identity, coefficient by coefficient")])


def q_iident(name, A, B, statement):
    bad = sum(1 for i in range(len(A)) for j in range(len(A[0])) if A[i][j] != B[i][j])
    qrow(name, [("statement", statement), ("holds", 1 if bad == 0 else 0), ("entries", len(A) * len(A[0])), ("entries_unequal", bad),
                ("method", "exact integer matrix identity")])


def iota_form(k, F):
    C = {}
    for I, a in F.items():
        if k in I:
            pos = I.index(k)
            J = tuple(x for x in I if x != k)
            C[J] = padd(C.get(J, {}), a, -1 if pos % 2 else 1)
    return fclean(C)


def postread(ctx):
    M = ctx["M"]
    Ls = ctx["Ls"]
    PHI = ctx["PHI"]
    TpI = ctx["TpI"]
    Psi7 = ctx["Psi7"]
    Iota = ctx["Iota"]
    nsq = ctx["nsq"]
    phi_form = ctx["phi_form"]
    psi = ctx["psi"]
    xivec = ctx["xivec"]
    xi_form = ctx["xi_form"]
    emit("POSTREAD_HEAD", [("reason", "sub-claims of the frozen rows that the blind block did not measure; computed here from the same definitions after the frozen rows were read, so they are not blind"),
                           ("blind_block_unchanged", 1)])
    # 7D.P5 structure of the witness
    n4 = 0
    pm2 = 0
    sp = 0
    asc = 0
    for i in LABELS:
        for j in LABELS:
            if i == j:
                continue
            s, c = M[i][j]
            D = iadd(iscal(Ls[c], s), imul(Ls[i], Ls[j]), -1)
            nzs = [(r, q, D[r][q]) for r in range(8) for q in range(8) if D[r][q]]
            if len(nzs) == 4:
                n4 += 1
            if nzs and all(abs(t[2]) == 2 for t in nzs):
                pm2 += 1
            other = set(range(8)) - {0, i, j, c}
            rs = [t[0] for t in nzs]
            cs = [t[1] for t in nzs]
            if len(nzs) == 4 and set(rs) == other and set(cs) == other and len(set(rs)) == 4 and len(set(cs)) == 4:
                sp += 1
            good = True
            for y in range(8):
                a1 = omul(M, omul(M, unit(i), unit(j)), unit(y))
                a2 = omul(M, unit(i), omul(M, unit(j), unit(y)))
                if [D[r][y] for r in range(8)] != [a1[t] - a2[t] for t in range(8)]:
                    good = False
            if good:
                asc += 1
    qrow("P5_witness_structure", [("ordered_pairs", 42), ("exactly_4_nonzero_entries", str(n4) + " of 42"),
                                  ("all_nonzero_entries_plus_or_minus_2", str(pm2) + " of 42"),
                                  ("twice_a_signed_permutation_of_the_other_4_basis_vectors", str(sp) + " of 42"),
                                  ("column_y_equals_associator_ei_ej_y", str(asc) + " of 42")])
    # 7D.P8 and 7D.P6
    q_iident("Psi7_equals_Iota_transpose", Psi7, iT(Iota), "Psi7 = Iota^T")
    q_iident("Psi7_Psi7T_is_3I", imul(Psi7, iT(Psi7)), iscal(ieye(7), 3), "Psi7 Psi7^T = 3 I_7")
    iu = [psum([pscale(pv(30 + k), PHI[k][i][j]) for k in LABELS if PHI[k][i][j]]) for (i, j) in PAIRS]
    ivv = [psum([pscale(pv(10 + k), PHI[k][i][j]) for k in LABELS if PHI[k][i][j]]) for (i, j) in PAIRS]
    q_ident("iota_u_phi_dot_iota_v_phi", [[psum([pmul(a, b) for a, b in zip(iu, ivv)])]],
            [[pscale(psum([pmul(pv(30 + a), pv(10 + a)) for a in LABELS]), 3)]], "u1..u7, v1..v7", "<iota_u phi, iota_v phi> = 3 <u, v>")
    tv = pm_imul(TpI, col(xivec))
    q_ident("nsq_TplusI_xi_wedge_v", [[psum([pmul(r_[0], r_[0]) for r_ in tv])]], [[pscale(psum([pmul(x, x) for x in xivec]), 3)]],
            "xi1..xi7, v1..v7", "nsq((T7 + I)(xi wedge v)) = 3 nsq(xi wedge v)")
    # 7D.P9
    nu = [pv(200 + b) for b in range(21)]
    om = pm_imul(TpI, col(nu))
    pw = pm_imul(Psi7, om)
    q_ident("nsq_Psi7_on_Lambda2_7", [[psum([pmul(r_[0], r_[0]) for r_ in pw])]], [[pscale(psum([pmul(r_[0], r_[0]) for r_ in om]), 3)]],
            "nu_1..nu_21 with omega = (T7 + I) nu ranging over Lambda2_7", "nsq(Psi7(omega)) = 3 nsq(omega) for omega in Lambda2_7")
    # 7D.P10
    E2 = Echelon()
    mats = []
    for (i, j) in PAIRS:
        Mij = imul(Ls[i], Ls[j])
        mats.append(Mij)
        E2.add(flat(Mij))
    r21 = len(E2.rows)
    grew = 0
    for p in range(21):
        for q in range(p + 1, 21):
            grew += E2.add(flat(iadd(imul(mats[p], mats[q]), imul(mats[q], mats[p]), -1)))
    ok2R = 0
    killeta = 0
    for kv, Rm in zip(ctx["kb"], ctx["Rk"]):
        cm = izero(8, 8)
        for b in range(21):
            if kv[b]:
                cm = iadd(cm, mats[b], kv[b])
        if cm == iscal(Rm, 2):
            ok2R += 1
        if all(cm[r][0] == 0 for r in range(8)):
            killeta += 1
    krank = rank_int([[mats[b][r][0] for b in range(21)] for r in range(8)])
    qrow("P10_spin7_span", [("span_dimension_of_L_ei_L_ej_i_lt_j", r21), ("commutators_checked", 210), ("span_grew", grew),
                            ("closed_under_commutators", 1 if grew == 0 else 0),
                            ("dimension_of_part_killing_eta", 21 - krank),
                            ("c_omega_equals_2_R_omega_on_ker_Psi7_basis", str(ok2R) + " of " + str(len(ctx["kb"]))),
                            ("c_omega_kills_eta_on_ker_Psi7_basis", str(killeta) + " of " + str(len(ctx["kb"])))])
    # 7D.P11
    TRI = list(itertools.combinations(LABELS, 3))
    colsm = []
    for I in TRI:
        F6 = fwedge({I: pc(1)}, phi_form)
        F7 = fwedge({I: pc(1)}, psi)
        colsm.append([F6.get(tuple(x for x in LABELS if x != m), {}).get((), 0) for m in LABELS] + [F7.get(VOL, {}).get((), 0)])
    Amat = iT(colsm)
    rA, k27 = nullspace_int(Amat, 35)
    l37 = []
    for k in LABELS:
        F = iota_form(k, psi)
        l37.append([F.get(I, {}).get((), 0) for I in TRI])
    phi_row = [phi_form.get(I, {}).get((), 0) for I in TRI]
    qrow("P11_Lambda3_split", [("dim_Lambda3", len(TRI)), ("dim_Lambda3_1_span_phi", rank_int([phi_row])),
                               ("dim_Lambda3_7_iota_v_psi", rank_int(l37)),
                               ("rank_of_gamma_to_gamma_wedge_phi_and_gamma_wedge_psi", rA), ("dim_Lambda3_27_kernel", len(k27)),
                               ("rank_of_union", rank_int([phi_row] + l37 + k27))])
    # DIMS clifford_dim and BOTT n = 7: a faithful module
    G = []
    for a in LABELS:
        Gm = izero(16, 16)
        for r in range(8):
            for q in range(8):
                Gm[r][q] = Ls[a][r][q]
                Gm[8 + r][8 + q] = -Ls[a][r][q]
        G.append(Gm)
    I16 = ieye(16)
    okc = 0
    for i in range(7):
        for j in range(7):
            if iadd(imul(G[i], G[j]), imul(G[j], G[i])) == iscal(I16, -2 if i == j else 0):
                okc += 1
    words = []
    for mask in range(128):
        P = I16
        for a in range(7):
            if (mask >> a) & 1:
                P = imul(P, G[a])
        words.append(flat(P))
    P = I16
    for a in range(7):
        P = imul(P, G[a])
    vol_ok = 1 if all(P[r][q] == ((-1 if r < 8 else 1) if r == q else 0) for r in range(16) for q in range(16)) else 0
    qrow("doubled_module_clifford_algebra", [("module", "O + O with gamma(e_i) = diag(L_ei, -L_ei)"), ("clifford_pairs_holding", str(okc) + " of 49"),
                                             ("rank_of_128_words", rank_int(words)), ("volume_product_is_diag_minus_I_plus_I", vol_ok),
                                             ("reading", "a faithful Cl_7 module of dimension 128, the volume element -1 on one block and +1 on the other")])
    # ADD R1 raw projection Gram identity
    Sraw = [[{} for _ in LABELS] for _ in PAIRS]
    for j in LABELS:
        wv = two_vec(fwedge(xi_form, {(j,): pc(1)}))
        cj = pm_imul(TpI, col(wv))
        for b in range(21):
            Sraw[b][j - 1] = cj[b][0]
    G7 = pm_mul(pm_T(Sraw), Sraw)
    base = pm_add(pm_eye(7, nsq), [[pmul(pv(a), pv(b)) for b in LABELS] for a in LABELS], -1)
    kfound = 0
    for kk in range(1, 10):
        if compare(G7, pm_scale(base, kk))[2] == 0:
            kfound = kk
            break
    g = math.gcd(kfound, 9) if kfound else 1
    q_ident("ADD_R1_raw_gram_cleared", G7, pm_scale(base, 3), "xi1..xi7",
            "(T7 + I)(xi wedge .)^T (T7 + I)(xi wedge .) = 3 (nsq(xi) I_7 - xi xi^T); with pi = (T7 + I) times one third, sigma^T sigma = xi xi^T + c (nsq I - xi xi^T)")
    qrow("ADD_R1_c", [("k_in_Gram_of_T_plus_I", kfound), ("c_numerator", exact_div(kfound, g)), ("c_denominator", exact_div(9, g)),
                      ("reading", "c = k over 9 since pi = (T7 + I) times one third")])


def compare_frozen(pred_path, add_path, bsha, mrows_list, ctx):
    emit("MEASURED_BLOCK_SEAL", [("precompare_sha256", PRECOMPARE_SEAL), ("this_run_sha256", bsha), ("identical", 1 if bsha == PRECOMPARE_SEAL else 0),
                                 ("precompare_file", "MEASURE-7D.precompare.hbp"), ("precompare_sealed_utc", "2026-11-09T17:32:13Z")])
    with open(pred_path, "rb") as fh:
        pred = fh.read()
    with open(add_path, "rb") as fh:
        add = fh.read()
    emit("FROZEN_INPUTS", [("prediction_sha256", sha_hex(pred)), ("prediction_expected", 1 if sha_hex(pred) == PRED_SHA else 0),
                           ("addendum_sha256", sha_hex(add)), ("addendum_expected", 1 if sha_hex(add) == ADD_SHA else 0)])
    if sha_hex(pred) != PRED_SHA or sha_hex(add) != ADD_SHA:
        emit("COMPARE_SUMMARY", [("status", "NOT_RUN"), ("reason", "frozen file hash mismatch")])
        return
    prows = {}
    for src, raw, fname in ((pred, "PREDICTION", "PREDICTION-SPINOR-SYMBOL.hbp"), (add, "ADDENDUM", "ADDENDUM-SPINOR-SYMBOL.hbp")):
        for line in src.decode("utf-8").split("\n"):
            if not line:
                continue
            parts = line.split("|")
            f = {}
            for p in parts[1:3]:
                if "=" in p:
                    k, v = p.split("=", 1)
                    f[k] = v
            tag = parts[0]
            if tag == "PRED":
                key = f.get("id")
            elif tag == "CONV":
                key = "CONV:" + f.get("id", "")
            elif tag == "ADD":
                key = "ADD:" + f.get("id", "")
            elif tag == "DIMS":
                key = "DIMS:" + f.get("case", "")
            elif tag == "BOTT":
                key = "BOTT:" + f.get("n", "")
            elif tag == "COMPARE":
                key = "COMPARE:" + f.get("item", "")
            elif tag == "PLANTED":
                key = "PLANTED:" + f.get("plant", "")
            else:
                key = tag
            prows[key] = (line, fname)
    MR = {}
    for r_ in mrows_list:
        parts = r_.split("|")
        f = {}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                f[k] = v
        if "name" in f:
            MR[f["name"]] = f
    postread(ctx)
    tally = {}
    btally = {}
    bars = [0]
    disagreements = []
    quote_missing = [0]

    def ids(*names):
        out = []
        for n in names:
            if n in MR:
                out.append(MR[n]["id"])
            elif n in QR:
                out.append(QR[n]["id"])
            else:
                raise ValueError("unknown measured row " + n)
        return ",".join(out)

    def m(name, key):
        return MR[name][key]

    def q(name, key):
        return QR[name][key]

    def citem(key, sub, quote, measured, mrows, ok, basis="BLIND", verdict=None):
        line, fname = prows[key]
        found = 1 if quote in line else 0
        if not found:
            quote_missing[0] += 1
        qe = quote
        rep = 0
        if "|" in qe:
            qe = qe.replace("|xi|^2", "nsq(xi)").replace("|", "BAR")
            rep = 1
            bars[0] += 1
        if verdict is None:
            verdict = "AGREE" if ok else "DISAGREE"
        tally[verdict] = tally.get(verdict, 0) + 1
        btally[verdict + "@" + basis] = btally.get(verdict + "@" + basis, 0) + 1
        if verdict == "DISAGREE":
            disagreements.append(key + ":" + sub)
        kv = [("frozen", key), ("sub", sub), ("frozen_file", fname), ("frozen_row_sha256", sha_hex(line.encode("utf-8"))),
              ("frozen_quote", qe), ("quote_found_in_frozen_row", found)]
        if rep:
            kv.append(("frozen_quote_bars_replaced", 1))
        kv += [("measured", measured), ("measured_rows", mrows), ("measured_basis", basis), ("verdict", verdict)]
        emit("COMPARE", kv)

    def h(name):
        return m(name, "holds") == "1"

    def hq(name):
        return q(name, "holds") == "1"
    nmono2 = len(list(itertools.combinations_with_replacement(LABELS, 2)))
    # 7D.P1
    slots = int(m("Lxi_squared_plus_nsq_is_zero", "entries")) * nmono2
    citem("7D.P1", "square_coefficients", "all 64 x 28 = 1792 coefficients of L_xi^2 + nsq(xi) I_8 are 0",
          "holds=1 over 64 entries; " + str(nmono2) + " degree-2 monomials in 7 variables give " + str(slots) + " slots, all 0; 56 nonzero coefficients on each side matched",
          ids("Lxi_squared_plus_nsq_is_zero"), h("Lxi_squared_plus_nsq_is_zero") and slots == 1792)
    citem("7D.P1", "skew", "also L_xi^T = -L_xi", "holds=1", ids("Lxi_skew"), h("Lxi_skew"))
    # 7D.P2
    citem("7D.P2", "clifford_basis", "28 of 28 relations i <= j hold", m("clifford_relations_on_basis", "ordered_pairs_holding") + " of 49 ordered pairs, which contain all 28 with i <= j",
          ids("clifford_relations_on_basis", "clifford_polarized"), m("clifford_relations_on_basis", "ordered_pairs_holding") == "49" and h("clifford_polarized"))
    # 7D.P3
    citem("7D.P3", "dimension_64", "dimension 64 = all of M_8(R) = End(O)", "dimension=" + m("algebra_generated_by_L_e1_to_L_e7", "dimension"),
          ids("algebra_generated_by_L_e1_to_L_e7"), m("algebra_generated_by_L_e1_to_L_e7", "dimension") == "64")
    citem("7D.P3", "not_128", "and not 128 = dim Cl_7", "the 128 ordered subset words span " + m("algebra_generated_by_L_e1_to_L_e7", "span_of_128_ordered_subset_products") + " in End(O); in the doubled module they span " + q("doubled_module_clifford_algebra", "rank_of_128_words"),
          ids("algebra_generated_by_L_e1_to_L_e7", "doubled_module_clifford_algebra"),
          m("algebra_generated_by_L_e1_to_L_e7", "span_of_128_ordered_subset_products") == "64" and q("doubled_module_clifford_algebra", "rank_of_128_words") == "128", "BLIND+POSTREAD")
    citem("7D.P3", "two_factors_one_module", "Cl_7 = M_8(R) + M_8(R) has two simple factors and O carries one of its two irreducible modules",
          "O: image M8(R), commutant 1, volume element -I; doubled module diag(L, -L): rank 128, volume diag(-I, +I)",
          ids("commutant_of_L_e1_to_L_e7", "volume_product_L_e1_L_e2_to_L_e7", "doubled_module_clifford_algebra"),
          m("commutant_of_L_e1_to_L_e7", "dimension") == "1" and q("doubled_module_clifford_algebra", "volume_product_is_diag_minus_I_plus_I") == "1", "BLIND+POSTREAD")
    citem("7D.P3", "commutant", "commutant dimension 1", "dimension=" + m("commutant_of_L_e1_to_L_e7", "dimension"), ids("commutant_of_L_e1_to_L_e7"),
          m("commutant_of_L_e1_to_L_e7", "dimension") == "1")
    # 7D.P4
    citem("7D.P4", "scalar", "epsilon7 I_8 with epsilon7 = +1 or -1", "L_e1 L_e2 ... L_e7 = " + m("volume_product_L_e1_L_e2_to_L_e7", "sign") + " I_8",
          ids("volume_product_L_e1_L_e2_to_L_e7"), m("volume_product_L_e1_L_e2_to_L_e7", "is_plus_or_minus_identity") == "1")
    citem("7D.P4", "sign", "CANNOT_PREDICT the sign here without the table and the ordering",
          "epsilon7 = " + m("volume_product_L_e1_L_e2_to_L_e7", "sign") + " for the frozen table, product taken left to right L_e1 L_e2 ... L_e7; for 8D.P4 to reuse",
          ids("volume_product_L_e1_L_e2_to_L_e7"), True, "BLIND", "SUPPLIED")
    # 7D.P5
    citem("7D.P5", "nonzero", "for each of the 42 ordered pairs i != j, L_(e_i e_j) - L_(e_i) L_(e_j) is nonzero",
          "rank histogram " + m("rank_L_eiej_minus_L_ei_L_ej", "rank_histogram"), ids("rank_L_eiej_minus_L_ei_L_ej"),
          m("rank_L_eiej_minus_L_ei_L_ej", "rank_histogram") == "4:42")
    citem("7D.P5", "four_entries_pm2", "with exactly 4 nonzero entries, each plus or minus 2",
          q("P5_witness_structure", "exactly_4_nonzero_entries") + " have 4 nonzero entries; " + q("P5_witness_structure", "all_nonzero_entries_plus_or_minus_2") + " have all entries plus or minus 2",
          ids("P5_witness_structure"), q("P5_witness_structure", "exactly_4_nonzero_entries") == "42 of 42" and q("P5_witness_structure", "all_nonzero_entries_plus_or_minus_2") == "42 of 42", "POSTREAD")
    citem("7D.P5", "associator", "it is the associator y -> [e_i, e_j, y]", q("P5_witness_structure", "column_y_equals_associator_ei_ej_y"),
          ids("P5_witness_structure"), q("P5_witness_structure", "column_y_equals_associator_ei_ej_y") == "42 of 42", "POSTREAD")
    citem("7D.P5", "quaternion_kernel", "zero on the quaternion subalgebra spanned by 1, e_i, e_j, e_i e_j",
          m("rank_L_eiej_minus_L_ei_L_ej", "kernel_is_quaternion_line_span_1_ei_ej_eiej"), ids("rank_L_eiej_minus_L_ei_L_ej"),
          m("rank_L_eiej_minus_L_ei_L_ej", "kernel_is_quaternion_line_span_1_ei_ej_eiej") == "42 of 42")
    citem("7D.P5", "signed_permutation", "twice a signed permutation on the other 4 basis vectors",
          q("P5_witness_structure", "twice_a_signed_permutation_of_the_other_4_basis_vectors"), ids("P5_witness_structure"),
          q("P5_witness_structure", "twice_a_signed_permutation_of_the_other_4_basis_vectors") == "42 of 42", "POSTREAD")
    citem("7D.P5", "clifford_still_holds", "yet the Clifford relations of 7D.P2 hold", "49 of 49", ids("clifford_relations_on_basis"),
          m("clifford_relations_on_basis", "holds") == "1")
    # 7D.P6
    citem("7D.P6", "sigma1", "sigma1(xi) v = Psi_7(xi wedge v) = xi cross v", "S1 = X holds=1", ids("G2_S1_equals_cross_matrix"), h("G2_S1_equals_cross_matrix"))
    citem("7D.P6", "sigma2", "sigma2(xi) w = <xi, w>", "S2 = xi^T holds=1; 3 S2 divisible by 3", ids("G2_S2_equals_xi_transpose", "G2_sequence_maps"),
          h("G2_S2_equals_xi_transpose") and m("G2_sequence_maps", "three_S2_coefficients_divisible_by_3") == "1")
    citem("7D.P6", "dims", "dimensions 1, 7, 7, 1", m("G2_symbol_sequence", "dims"), ids("G2_symbol_sequence"), m("G2_symbol_sequence", "dims") == "1,7,7,1")
    citem("7D.P6", "euler", "Euler count 1 - 7 + 7 - 1 = 0", m("G2_symbol_sequence", "euler_count"), ids("G2_symbol_sequence"), m("G2_symbol_sequence", "euler_count") == "0")
    citem("7D.P6", "composites", "composites xi cross xi = 0 and <xi, xi cross v> = 0 as polynomial identities", "S1 S0 = 0 and S2 S1 = 0 hold",
          ids("G2_S1_S0_zero", "G2_S2_S1_zero"), h("G2_S1_S0_zero") and h("G2_S2_S1_zero"))
    citem("7D.P6", "ranks", "ranks 1, 6, 1 for every xi != 0", m("G2_symbol_sequence", "ranks_for_every_nonzero_xi"),
          ids("G2_symbol_sequence", "G2_laplacian_Lambda1", "G2_laplacian_Lambda2_7", "G2_S0T_S0", "G2_S2_S2T"), m("G2_symbol_sequence", "ranks_for_every_nonzero_xi") == "1,6,1")
    citem("7D.P6", "reason_factor_3", "<xi wedge iota_w phi, phi> = <iota_w phi, iota_xi phi> = 3 <w, xi>",
          "outer equality from 3 S2 = <xi wedge iota_w phi, phi> = 3 <xi, w>; middle from <iota_u phi, iota_v phi> = 3 <u, v>",
          ids("G2_S2_equals_xi_transpose", "iota_u_phi_dot_iota_v_phi"), h("G2_S2_equals_xi_transpose") and hq("iota_u_phi_dot_iota_v_phi"), "BLIND+POSTREAD")
    # 7D.P7
    citem("7D.P7", "skew", "X(xi)^T = -X(xi)", "holds=1", ids("cross_matrix_skew"), h("cross_matrix_skew"))
    citem("7D.P7", "square", "X(xi)^2 = -nsq(xi) I_7 + xi xi^T as polynomial identities", "holds=1", ids("cross_matrix_square"), h("cross_matrix_square"))
    citem("7D.P7", "rank6", "ker X(xi) = span(xi), im X(xi) = xi^perp, rank 6 for every xi != 0", "rank 6 certified; X xi = 0 holds",
          ids("cross_matrix_rank_certificate", "cross_matrix_kills_xi"), m("cross_matrix_rank_certificate", "rank_for_every_nonzero_xi") == "6" and h("cross_matrix_kills_xi"))
    citem("7D.P7", "exactness", "exactness at Lambda1 and at Lambda2_7", "exact=" + m("G2_symbol_sequence", "exact_for_every_nonzero_xi"), ids("G2_symbol_sequence"),
          m("G2_symbol_sequence", "exact_for_every_nonzero_xi") == "1")
    # 7D.P8
    citem("7D.P8", "minimal_polynomial", "(T7 - 2 I)(T7 + I) = 0 on Lambda2", "holds=1 with vol = -(e1 wedge ... wedge e7) as 7D.C3 directs", ids("T7_minimal_polynomial", "orientation_7D"),
          m("T7_minimal_polynomial", "holds") == "1")
    citem("7D.P8", "trace", "trace T7 = 2 x 7 - 14 = 0", "trace=" + m("T7_eigen_split", "trace") + ", multiplicities " + m("T7_eigen_split", "eigenvalue_plus2_multiplicity") + " and " + m("T7_eigen_split", "eigenvalue_minus1_multiplicity"),
          ids("T7_eigen_split"), m("T7_eigen_split", "trace") == "0" and m("T7_eigen_split", "eigenvalue_plus2_multiplicity") == "7")
    citem("7D.P8", "3pi7", "3 pi_7 = T7 + I", "(T7 + I)^2 = 3 (T7 + I) and (T7 + I) = Iota Psi7 hold", ids("TplusI_squared", "TplusI_equals_Iota_Psi7"),
          m("TplusI_squared", "holds") == "1" and m("TplusI_equals_Iota_Psi7", "holds") == "1")
    citem("7D.P8", "bilinear", "(T7 + I)(xi wedge v) = iota_(xi cross v) phi as a bilinear identity", "holds=1", ids("TplusI_xi_wedge_v_equals_iota_cross_phi"),
          h("TplusI_xi_wedge_v_equals_iota_cross_phi"))
    citem("7D.P8", "norm", "nsq((T7 + I)(xi wedge v)) = 3 nsq(xi wedge v)", "holds=" + q("nsq_TplusI_xi_wedge_v", "holds"), ids("nsq_TplusI_xi_wedge_v"),
          hq("nsq_TplusI_xi_wedge_v"), "POSTREAD")
    citem("7D.P8", "Psi_TplusI", "Psi_7 (T7 + I) = 3 Psi_7", "holds=1", ids("Psi7_TplusI_is_3Psi7"), m("Psi7_TplusI_is_3Psi7", "holds") == "1")
    citem("7D.P8", "Psi_PsiT", "Psi_7 Psi_7^T = 3 I_7", "Psi7 Iota = 3 I blind; Psi7 = Iota^T and Psi7 Psi7^T = 3 I postread", ids("Psi7_Iota_is_3I", "Psi7_equals_Iota_transpose", "Psi7_Psi7T_is_3I"),
          m("Psi7_Iota_is_3I", "holds") == "1" and hq("Psi7_Psi7T_is_3I") and hq("Psi7_equals_Iota_transpose"), "BLIND+POSTREAD")
    citem("7D.P8", "iota_inner", "<iota_u phi, iota_v phi> = 3 <u, v>", "holds=" + q("iota_u_phi_dot_iota_v_phi", "holds"), ids("iota_u_phi_dot_iota_v_phi"),
          hq("iota_u_phi_dot_iota_v_phi"), "POSTREAD")
    citem("7D.P8", "kernel", "ker Psi_7 = Lambda2_14 = der(O) as rotations", "ker Psi7 = -1 eigenspace (dim 14) and R(ker Psi7) = der(O)",
          ids("Psi7_rank_and_kernel", "der_O_equals_R_of_ker_Psi7"), m("Psi7_rank_and_kernel", "kernel_equals_minus1_eigenspace") == "1" and m("der_O_equals_R_of_ker_Psi7", "equal") == "1")
    # 7D.P9
    citem("7D.P9", "R_formula", "R(xi)(f, v) = (-<xi, v>, f xi + xi cross v)", "R built from this formula equals L_xi", ids("rolled_up_R_equals_Lxi"), h("rolled_up_R_equals_Lxi"))
    citem("7D.P9", "R_equals_L", "R(xi) = L_xi in the basis 1, e1..e7 for every xi in Im O", "polynomial identity in xi1..xi7 holds", ids("rolled_up_R_equals_Lxi"), h("rolled_up_R_equals_Lxi"))
    citem("7D.P9", "folding", "folding E2 = Lambda2_7 onto Lambda1 by w and E3 = Lambda3_1 onto Lambda0 by t", "[[0, -S2],[S0, S1]] = L_xi holds",
          ids("R2_cross_product_coordinates"), h("R2_cross_product_coordinates"))
    citem("7D.P9", "R_square", "R(xi)^2 = -nsq(xi) I_8", "holds=1", ids("R_squared"), h("R_squared"))
    citem("7D.P9", "R_gram", "R(xi)^T R(xi) = nsq(xi) I_8", "holds=1", ids("R_transpose_R"), h("R_transpose_R"))
    citem("7D.P9", "sqrt3", "nsq(Psi_7(omega)) = 3 nsq(omega) on Lambda2_7", "holds=" + q("nsq_Psi7_on_Lambda2_7", "holds"), ids("nsq_Psi7_on_Lambda2_7"),
          hq("nsq_Psi7_on_Lambda2_7"), "POSTREAD")
    # 7D.P10
    citem("7D.P10", "dim14", "for each D in der(O) (dimension 14", "dimension=" + m("derivations_der_O", "dimension"), ids("derivations_der_O"), m("derivations_der_O", "dimension") == "14")
    citem("7D.P10", "D1_and_bracket", "D 1 = 0, [D, L_xi] = L_(D xi)", m("derivations_der_O", "kill_1_and_preserve_Im_O") + " kill 1; bracket " + m("derivation_equivariance_D_Lxi", "holding_imaginary_xi"),
          ids("derivations_der_O", "derivation_equivariance_D_Lxi"), m("derivations_der_O", "kill_1_and_preserve_Im_O") == "14 of 14" and m("derivation_equivariance_D_Lxi", "holding_imaginary_xi") == "14 of 14")
    citem("7D.P10", "D_R", "so [D, R(xi)] = R(D xi)", "follows from R(xi) = L_xi and [D, L_xi] = L_(D xi), both measured",
          ids("rolled_up_R_equals_Lxi", "derivation_equivariance_D_Lxi"), h("rolled_up_R_equals_Lxi") and m("derivation_equivariance_D_Lxi", "holding_imaginary_xi") == "14 of 14")
    citem("7D.P10", "span21", "the span of L_(e_i) L_(e_j), i < j, has dimension 21 and is closed under commutators (spin(7) on S)",
          "span " + q("P10_spin7_span", "span_dimension_of_L_ei_L_ej_i_lt_j") + ", closed=" + q("P10_spin7_span", "closed_under_commutators"), ids("P10_spin7_span"),
          q("P10_spin7_span", "span_dimension_of_L_ei_L_ej_i_lt_j") == "21" and q("P10_spin7_span", "closed_under_commutators") == "1", "POSTREAD")
    citem("7D.P10", "kill_eta_14", "its part killing eta = 1 has dimension 14 and is the image of Lambda2_14",
          "ker Psi7 = 14 = -1 eigenspace blind; part of the span killing eta = " + q("P10_spin7_span", "dimension_of_part_killing_eta") + " postread",
          ids("Psi7_is_clifford_action_on_eta", "Psi7_rank_and_kernel", "P10_spin7_span"),
          m("Psi7_rank_and_kernel", "kernel_dimension") == "14" and q("P10_spin7_span", "dimension_of_part_killing_eta") == "14", "BLIND+POSTREAD")
    citem("7D.P10", "c_is_2R", "sum_(i<j) omega_ij L_(e_i) L_(e_j) = 2 R_omega on Im O and 0 on 1, a derivation",
          "c(omega) = 2 R_omega on " + q("P10_spin7_span", "c_omega_equals_2_R_omega_on_ker_Psi7_basis") + " kernel basis forms; R_omega derivations " + m("der_O_equals_R_of_ker_Psi7", "ker_Psi7_rotations_that_are_derivations"),
          ids("P10_spin7_span", "der_O_equals_R_of_ker_Psi7"),
          q("P10_spin7_span", "c_omega_equals_2_R_omega_on_ker_Psi7_basis") == "14 of 14" and m("der_O_equals_R_of_ker_Psi7", "ker_Psi7_rotations_that_are_derivations") == "14 of 14", "BLIND+POSTREAD")
    citem("7D.P10", "half_c_bracket", "[half of c(omega), L_v] = L_(R_omega v)", "[c(omega), L_v] = 2 L_(R_omega v) holds", ids("C0_5_rotation_bracket"), h("C0_5_rotation_bracket"))
    # 7D.P11
    citem("7D.P11", "S", "S = Lambda0 + Lambda1 = 1 + 7 = R + R^7 = O", "dim S 8; R(xi) on Lambda0 + Lambda1 equals L_xi on O", ids("dims_7D", "rolled_up_R_equals_Lxi"),
          m("dims_7D", "dim_S") == "8" and h("rolled_up_R_equals_Lxi"))
    citem("7D.P11", "Lambda2", "Lambda2 = Lambda2_7 + Lambda2_14 = 7 + 14 with Lambda2_14 = g2", m("dims_7D", "dim_Lambda2_7") + " + " + m("dims_7D", "dim_Lambda2_14") + ", Lambda2_14 = ker Psi7 = der(O)",
          ids("dims_7D", "der_O_equals_R_of_ker_Psi7"), m("dims_7D", "dim_Lambda2_7") == "7" and m("dims_7D", "dim_Lambda2_14") == "14" and m("der_O_equals_R_of_ker_Psi7", "equal") == "1")
    citem("7D.P11", "Lambda3", "Lambda3 = 1 + 7 + 27 = 35", q("P11_Lambda3_split", "dim_Lambda3_1_span_phi") + " + " + q("P11_Lambda3_split", "dim_Lambda3_7_iota_v_psi") + " + " + q("P11_Lambda3_split", "dim_Lambda3_27_kernel") + ", union rank " + q("P11_Lambda3_split", "rank_of_union"),
          ids("P11_Lambda3_split"), q("P11_Lambda3_split", "dim_Lambda3_27_kernel") == "27" and q("P11_Lambda3_split", "dim_Lambda3_7_iota_v_psi") == "7" and q("P11_Lambda3_split", "rank_of_union") == "35", "POSTREAD")
    # DIMS case=7D, read field by field from the frozen row
    dline = prows["DIMS:7D"][0]
    dfields = []
    for p in dline.split("|")[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            dfields.append((k, v))
    gen_size = math.isqrt(int(m("Lxi_squared_plus_nsq_is_zero", "entries")))
    meas = {"n": (str(len(LABELS)), "Lxi_squared_plus_nsq_is_zero", "BLIND"),
            "spinor_module": (m("dims_7D", "dim_S"), "dims_7D", "BLIND"),
            "clifford_dim": (q("doubled_module_clifford_algebra", "rank_of_128_words"), "doubled_module_clifford_algebra", "POSTREAD"),
            "generated_algebra": (m("algebra_generated_by_L_e1_to_L_e7", "dimension"), "algebra_generated_by_L_e1_to_L_e7", "BLIND"),
            "commutant": (m("commutant_of_L_e1_to_L_e7", "dimension"), "commutant_of_L_e1_to_L_e7", "BLIND"),
            "generators": (str(len(LABELS)), "clifford_relations_on_basis", "BLIND"),
            "generator_size": (str(gen_size), "Lxi_squared_plus_nsq_is_zero", "BLIND"),
            "lambda2_small": (m("dims_7D", "dim_Lambda2_7"), "dims_7D", "BLIND"),
            "lambda2_big": (m("dims_7D", "dim_Lambda2_14"), "dims_7D", "BLIND"),
            "T_small": ("2" if (m("T7_minimal_polynomial", "holds") == "1" and m("T7_eigen_split", "eigenvalue_plus2_multiplicity") == "7") else "NOT_2", "T7_eigen_split", "BLIND"),
            "T_big": ("-1" if (m("T7_minimal_polynomial", "holds") == "1" and m("T7_eigen_split", "eigenvalue_minus1_multiplicity") == "14") else "NOT_-1", "T7_eigen_split", "BLIND"),
            "k": ("3" if m("TplusI_squared", "holds") == "1" else "NOT_3", "TplusI_squared", "BLIND"),
            "stabilizer_of_eta": (m("Psi7_rank_and_kernel", "kernel_dimension"), "Psi7_rank_and_kernel", "BLIND"),
            "sequence_dims": (m("G2_symbol_sequence", "dims"), "G2_symbol_sequence", "BLIND"),
            "symbol_ranks": (m("G2_symbol_sequence", "ranks_for_every_nonzero_xi"), "G2_symbol_sequence", "BLIND"),
            "square_coefficients": (str(slots) if h("Lxi_squared_plus_nsq_is_zero") else "NOT_ZERO", "Lxi_squared_plus_nsq_is_zero", "BLIND")}
    for k, v in dfields:
        if k in ("case", "status", "claim", "json"):
            continue
        if k in meas:
            mvv, mname, basis = meas[k]
            citem("DIMS:7D", k, k + "=" + v, mvv, ids(mname), mvv == v, basis)
        else:
            citem("DIMS:7D", k, k + "=" + v, "no measured counterpart", "none", False, "NONE", "NOT_MEASURED")
    for item, mname, key in (("7D_sequence_dims", "G2_symbol_sequence", "dims"), ("7D_ranks", "G2_symbol_sequence", "ranks_for_every_nonzero_xi"),
                             ("7D_lambda2_7", "dims_7D", "dim_Lambda2_7"), ("7D_lambda2_14", "dims_7D", "dim_Lambda2_14")):
        cline = prows["COMPARE:" + item][0]
        tt = [p.split("=", 1)[1] for p in cline.split("|") if p.startswith("this_table=")][0]
        citem("COMPARE:" + item, "this_table", "this_table=" + tt, m(mname, key), ids(mname), m(mname, key) == tt)
    citem("BOTT:7", "cl_square_minus_1", "cl_square_minus_1=M_8(R)+M_8(R)", "doubled module: rank 128, volume diag(-I, +I); each block generates M8(R) (dimension 64, commutant 1)",
          ids("doubled_module_clifford_algebra", "algebra_generated_by_L_e1_to_L_e7"),
          q("doubled_module_clifford_algebra", "rank_of_128_words") == "128" and m("algebra_generated_by_L_e1_to_L_e7", "dimension") == "64", "BLIND+POSTREAD")
    citem("BOTT:7", "dim_minus", "dim_minus=128", q("doubled_module_clifford_algebra", "rank_of_128_words"), ids("doubled_module_clifford_algebra"),
          q("doubled_module_clifford_algebra", "rank_of_128_words") == "128", "POSTREAD")
    citem("BOTT:7", "irreducible_module_minus", "irreducible_module_minus=8", "O: dimension 8, commutant 1, image all of M8(R)",
          ids("commutant_of_L_e1_to_L_e7", "algebra_generated_by_L_e1_to_L_e7"), m("commutant_of_L_e1_to_L_e7", "dimension") == "1")
    citem("BOTT:7", "cl_square_plus_1", "cl_square_plus_1=M_8(C)", "the plus-sign Clifford algebra is not built in this measurement", "none", False, "NONE", "NOT_MEASURED")
    citem("PLANTED:lambda2_14_g2_planted_as_15", "true_value", "true_value=14", "dim Lambda2_14 = " + m("dims_7D", "dim_Lambda2_14") + "; the planted 15 differs from the measurement",
          ids("dims_7D"), m("dims_7D", "dim_Lambda2_14") == "14")
    citem("PLANTED:ranks_7d_planted_as_1_5_1", "true_value", "true_value=1;6;1", "ranks " + m("G2_symbol_sequence", "ranks_for_every_nonzero_xi") + "; the planted 1;5;1 differs from the measurement",
          ids("G2_symbol_sequence"), m("G2_symbol_sequence", "ranks_for_every_nonzero_xi") == "1,6,1")
    # conventions the measurement was told to test or report
    citem("CONV:7D.C3", "orientation", "if e1 wedge ... wedge e7 gives -2 and +1, the measurement reverses it and says so",
          "e1 wedge ... wedge e7 gave " + m("T7_with_vol_e1_to_e7", "eigenvalues_for_e1_to_e7") + "; reversed=" + m("orientation_7D", "reversed"),
          ids("T7_with_vol_e1_to_e7", "orientation_7D"), m("orientation_7D", "reversed") == "1" and m("T7_with_vol_e1_to_e7", "eigenvalues_for_e1_to_e7") == "-2 on 7, +1 on 14", "BLIND", "REPORTED_REVERSED")
    citem("CONV:C0.5", "derived_bracket", "then [c(omega), gamma(v)] = 2 gamma(R_omega v) by the Clifford relations alone", "holds=1 in 7D", ids("C0_5_rotation_bracket"), h("C0_5_rotation_bracket"))
    citem("CONV:C0.6", "stabilizer", "its kernel is the Lie algebra of the stabilizer of eta", "ker Psi7 = 14-dim, R(ker Psi7) = der(O), closed Lie algebra",
          ids("der_O_equals_R_of_ker_Psi7", "der_O_lie_closure"), m("der_O_equals_R_of_ker_Psi7", "equal") == "1" and m("der_O_lie_closure", "closed") == "1")
    citem("CONV:7D.C5", "pi7", "then pi_7 omega = one third of iota_w phi", "(T7 + I) = Iota Psi7 and (T7 + I)^2 = 3 (T7 + I)", ids("TplusI_equals_Iota_Psi7", "TplusI_squared"),
          m("TplusI_equals_Iota_Psi7", "holds") == "1" and m("TplusI_squared", "holds") == "1")
    # addendum
    citem("ADD:R1", "identity", "sigma^T_sigma=xi_xi^T+c(|xi|^2_I-xi_xi^T)", "(T7 + I)(xi wedge .)^T (T7 + I)(xi wedge .) = 3 (nsq I - xi xi^T) holds=" + q("ADD_R1_raw_gram_cleared", "holds"),
          ids("ADD_R1_raw_gram_cleared"), hq("ADD_R1_raw_gram_cleared"), "POSTREAD")
    citem("ADD:R1", "c_G2", "1:3(G2)", "c = " + q("ADD_R1_c", "c_numerator") + ":" + q("ADD_R1_c", "c_denominator") + " postread; the blind block's pi_7(xi wedge v) = c iota_(xi x v) phi has the same c = " + m("R1_constant_c", "c_numerator") + ":" + m("R1_constant_c", "c_denominator") + " in a different identity",
          ids("ADD_R1_c", "R1_constant_c"), q("ADD_R1_c", "c_numerator") == "1" and q("ADD_R1_c", "c_denominator") == "3", "BLIND+POSTREAD")
    citem("ADD:R2", "integer_isometry", "sigma^T_sigma=|xi|^2_I_exactly,with_integer_matrices", "S0 S0^T + S1^T S1 = nsq I_7 holds (sigma = (-S0^T, S1))",
          ids("G2_laplacian_Lambda1"), h("G2_laplacian_Lambda1"))
    citem("ADD:R2", "cross_coordinates", "7D_Lambda2_7=R^7_through_iota(xi^a)phi=xi_x_a(the_octonion_cross_product)", "Psi7(xi wedge a) = xi x a holds",
          ids("G2_S1_equals_cross_matrix"), h("G2_S1_equals_cross_matrix"))
    citem("ADD:R2", "raw_not_integer", "with_the_raw_projection_the_G2_case_needs_sqrt3_and_cannot_be_integer", "raw Gram has c = 1:3, so an isometry needs a factor sqrt3 on Lambda2_7",
          ids("ADD_R1_c"), q("ADD_R1_c", "c_denominator") == "3", "POSTREAD")
    citem("ADD:R3", "7D_alternativity", "7D=xi(xi_x)=(xi_xi)x=-|xi|^2_x_for_imaginary_xi,alternativity_only", "L_xi L_xi = -nsq I holds",
          ids("Lxi_squared_plus_nsq_is_zero"), h("Lxi_squared_plus_nsq_is_zero"))
    r4 = MR["R4_associator_e1_e2_e3"]
    citem("ADD:R4", "values", "(e1e2)e3=-e6,e1(e2e3)=+e6,associator=-2e6", "(e1e2)e3=" + r4["(e1e2)e3"] + ", e1(e2e3)=" + r4["e1(e2e3)"] + ", associator=" + r4["associator"],
          ids("R4_associator_e1_e2_e3"), r4["(e1e2)e3"] == "-1e6" and r4["e1(e2e3)"] == "+1e6" and r4["associator"] == "-2e6")
    citem("ADD:R4", "off_line", "triple=(xi,xi',x)=(e1,e2,e3),not_on_a_common_line", "phi(e1, e2, e3) = " + r4["phi_e1_e2_e3"] + "; line through 1 and 2 is 1,2,4",
          ids("R4_associator_e1_e2_e3", "octonion_table"), r4["phi_e1_e2_e3"] == "0" and "1,2,4" in m("octonion_table", "lines"))
    citem("ADD:R5", "complex", "the_symbol_sequence_is_a_complex_because_sigma(xi)sigma(xi)=xi^xi^=0", "S1 S0 = 0 and S2 S1 = 0 hold", ids("G2_S1_S0_zero", "G2_S2_S1_zero"),
          h("G2_S1_S0_zero") and h("G2_S2_S1_zero"))
    citem("ADD:R5", "euler_7D", "1-7+7-1=0(7D,Reyes_Carrion)", m("G2_symbol_sequence", "euler_count"), ids("G2_symbol_sequence"), m("G2_symbol_sequence", "euler_count") == "0")
    citem("ADD:R5", "instantons", "the_operator_sequence_is_a_complex_only_on_instantons", "operator level, outside a symbol measurement", "none", False, "NONE", "NOT_MEASURED")
    citem("ADD:R6", "negative_convention", "real_Clifford_multiplication_uses_the_negative-definite_convention_gamma(xi)^2=-|xi|^2_I", "L_xi^2 = -nsq I holds",
          ids("Lxi_squared_plus_nsq_is_zero"), h("Lxi_squared_plus_nsq_is_zero"))
    citem("ADD:R6", "doubled_square", "the_rolled-up_first-order_symbols(M_in_4D,R_in_7D,N_in_8D)_and_their_doubled_self-adjoint_forms_square_to_+|xi|^2",
          "doubled self-adjoint [[0, R^T],[R, 0]] squares to +nsq I_16 (R = L_xi)", ids("doubled_selfadjoint_symbol", "rolled_up_R_equals_Lxi"),
          h("doubled_selfadjoint_symbol") and h("rolled_up_R_equals_Lxi"))
    citem("ADD:R6", "R_7D_square_literal", "the_rolled-up_first-order_symbols(M_in_4D,R_in_7D,N_in_8D)_and_their_doubled_self-adjoint_forms_square_to_+|xi|^2",
          "R(xi) R(xi) = -nsq(xi) I_8 (as 7D.P9 states); +nsq(xi) holds for R^T R and for the doubled form, not for R R",
          ids("R_squared", "R_transpose_R"), not h("R_squared"))
    citem("ADD:DISCLOSURE", "independence", "the_measurement_here_is_a_separate_implementation_that_does_not_read_that_receipt_before_measuring",
          "this seat built every matrix itself and read neither the receipt nor another agent script", "none", True, "PROCESS")
    # rows the task text adds beyond the frozen files
    emit("COMPARE_TASK", [("item", "signed_permutations_1344"), ("task_text", "the 1,344 table-preserving signed permutations g satisfy g L_xi g^-1 = L_(g xi)"),
                          ("measured", m("table_preserving_signed_permutations", "preserving_all_42_offdiagonal_products") + " found among " + m("table_preserving_signed_permutations", "candidates_examined") + "; equivariant " + m("signed_permutation_equivariance", "holding")),
                          ("measured_rows", ids("table_preserving_signed_permutations", "signed_permutation_equivariance")),
                          ("verdict", "AGREE" if (m("table_preserving_signed_permutations", "preserving_all_42_offdiagonal_products") == "1344" and m("signed_permutation_equivariance", "holding") == "1344 of 1344") else "DISAGREE")])
    emit("COMPARE_TASK", [("item", "octonion_table_lines"), ("task_text", "lines (1,2,4)(2,3,5)(3,4,6)(4,5,7)(5,6,1)(6,7,2)(7,1,3)"),
                          ("measured", "rule e_i e_(i+1) = e_(i+3) gives " + m("octonion_table", "lines")), ("measured_rows", ids("octonion_table")),
                          ("verdict", "AGREE" if m("octonion_table", "lines_equal_frozen_list") == "1" else "DISAGREE")])
    emit("FORMAT_FINDING", [("file", "ADDENDUM-SPINOR-SYMBOL.hbp"), ("rows", "ADD R1, R2, R3, R6"),
                            ("finding", "values hold the bar character on both sides of xi for its norm although the prediction FORMAT row says the bar notation is avoided because the bar separates fields; a split on the bar breaks those rows"),
                            ("effect_here", "quotes were matched against the whole row; emitted quotes replace the bar"), ("quotes_with_bars", bars[0]), ("claim", "MEASURED")])
    keys = sorted(tally)
    emit("COMPARE_SUMMARY", [("compare_rows", sum(tally[k] for k in keys)),
                             ("by_verdict", ";".join(k + ":" + str(tally[k]) for k in keys)),
                             ("by_verdict_and_basis", ";".join(k + ":" + str(btally[k]) for k in sorted(btally))),
                             ("disagreements", ";".join(disagreements) if disagreements else "none"),
                             ("quotes_not_found_in_frozen_rows", quote_missing[0]),
                             ("blind_block_identical_to_precompare_seal", 1 if bsha == PRECOMPARE_SEAL else 0),
                             ("postread_rows", QCOUNT[0]),
                             ("reading", "AGREE and DISAGREE compare a measured value with the quoted frozen text; SUPPLIED fills a value the prediction left open; REPORTED_REVERSED is the orientation reversal 7D.C3 provides for; NOT_MEASURED is outside this measurement; PROCESS is a statement about how this run was made")])


if __name__ == "__main__":
    main()
