# python3 -I -B m1_octonion.py
# M1 octonion step: the G2-invariant almost complex structure on G2 mod SU3 = S^6.
# Integer arithmetic only. Gaussian integers are pairs of integers.
# Linear algebra is exact fraction-free elimination over Z (gcd-normalised rows).
# E=0: compute only, no network, writes only the .hbp and its .sha256 sidecar here.
# Run: python3 -I -B m1_octonion.py

import sys
import os
import io
import hashlib
import tokenize
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "M1-OCTONION.hbp"
OUT = os.path.join(HERE, OUT_NAME)
SHA_OUT = os.path.join(HERE, OUT_NAME + ".sha256")

ROWS = []
STATE = {"checks": 0, "fail": 0, "failed_tags": []}


# ------------------------------------------------------------------ rows
def clean(v):
    s = str(v)
    s = s.replace("|", ";").replace("=", ":").replace("\n", "_").replace("\r", "_")
    s = s.replace(" ", "_")
    return s


def emit(tag, **kv):
    parts = [tag]
    for k, v in kv.items():
        parts.append("%s=%s" % (k, clean(v)))
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


def fv(v):
    return "[" + ",".join(str(x) for x in v) + "]"


def fg(z):
    # Gaussian integer pair (re, im) as text
    return "(%d%+di)" % (z[0], z[1])


# ------------------------------------------------------------------ self lint
def self_lint():
    return lint_source(open(os.path.abspath(__file__), "rb").read())


def lint_source(src):
    bad = []
    sl = chr(47)
    banned_ops = {sl, sl + sl, sl + "=", sl + sl + "="}
    banned_names = {"float", "Fraction", "Decimal", "complex", "fractions", "decimal", "cmath"}
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
            line = t.line.strip()
            if line != "from math import gcd":
                bad.append(("math", t.start[0]))
    slash_chars = src.count(b"\x2f")
    return bad, slash_chars


# ------------------------------------------------------------------ integer helpers
def exact_div(a, b):
    q, r = divmod(a, b)
    if r != 0:
        raise ArithmeticError("inexact division %d by %d" % (a, b))
    return q


def lcm(a, b):
    a, b = abs(a), abs(b)
    if a == 0 or b == 0:
        return max(a, b)
    return exact_div(a * b, gcd(a, b))


def vgcd(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    return g


def primitive(v):
    g = vgcd(v)
    if g == 0:
        return list(v)
    w = [exact_div(x, g) for x in v]
    for x in w:
        if x != 0:
            if x < 0:
                w = [-y for y in w]
            break
    return w


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def vscale(k, a):
    return [k * x for x in a]


def is_zero_vec(v):
    return all(x == 0 for x in v)


def rref(rows, ncols):
    A = [list(r) for r in rows if any(r)]
    piv = []
    r = 0
    for c in range(ncols):
        p = -1
        for i in range(r, len(A)):
            if A[i][c] != 0:
                p = i
                break
        if p < 0:
            continue
        A[r], A[p] = A[p], A[r]
        pr = A[r]
        g0 = vgcd(pr)
        if g0 > 1:
            pr = [exact_div(x, g0) for x in pr]
            A[r] = pr
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                a = pr[c]
                b = A[i][c]
                g = gcd(a, b)
                ma = exact_div(a, g)
                mb = exact_div(b, g)
                new = [ma * x - mb * y for x, y in zip(A[i], pr)]
                g2 = vgcd(new)
                if g2 > 1:
                    new = [exact_div(x, g2) for x in new]
                A[i] = new
        piv.append(c)
        r += 1
        A = A[:r] + [row for row in A[r:] if any(row)]
        if r >= len(A):
            break
    return A[:r], piv


def nullspace(rows, ncols):
    R, piv = rref(rows, ncols)
    pset = set(piv)
    free = [c for c in range(ncols) if c not in pset]
    L = 1
    for i, c in enumerate(piv):
        L = lcm(L, R[i][c])
    basis = []
    for f in free:
        x = [0] * ncols
        x[f] = L
        for i, c in enumerate(piv):
            x[c] = -exact_div(R[i][f] * L, R[i][c])
        basis.append(primitive(x))
    return basis, len(piv)


def rank_of(vectors, ncols):
    _, piv = rref(vectors, ncols)
    return len(piv)


def det_bareiss(M):
    A = [list(r) for r in M]
    n = len(A)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = -1
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    sw = i
                    break
            if sw < 0:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = exact_div(A[i][j] * A[k][k] - A[i][k] * A[k][j], prev)
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def adjugate(M):
    n = len(M)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor = [[M[r][c] for c in range(n) if c != i] for r in range(n) if r != j]
            s = 1 if (i + j) % 2 == 0 else -1
            adj[i][j] = s * det_bareiss(minor)
    return adj


# ------------------------------------------------------------------ real matrices
def zmat(n, m=None):
    if m is None:
        m = n
    return [[0] * m for _ in range(n)]


def ident(n):
    I = zmat(n)
    for i in range(n):
        I[i][i] = 1
    return I


def mmul(A, B):
    n = len(A)
    k = len(B)
    m = len(B[0])
    C = zmat(n, m)
    for i in range(n):
        Ai = A[i]
        Ci = C[i]
        for t in range(k):
            a = Ai[t]
            if a == 0:
                continue
            Bt = B[t]
            for j in range(m):
                if Bt[j] != 0:
                    Ci[j] += a * Bt[j]
    return C


def madd(A, B):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def msub(A, B):
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(k, A):
    return [[k * x for x in r] for r in A]


def mneg(A):
    return mscale(-1, A)


def mcomm(A, B):
    return msub(mmul(A, B), mmul(B, A))


def mvec(A, v):
    return [dot(r, v) for r in A]


def mT(A):
    return [list(r) for r in zip(*A)]


def meq(A, B):
    return all(x == y for ra, rb in zip(A, B) for x, y in zip(ra, rb))


def mzero(A):
    return all(x == 0 for r in A for x in r)


def tr_prod(A, B):
    n = len(A)
    return sum(A[i][j] * B[j][i] for i in range(n) for j in range(n))


def flat(A):
    return [x for r in A for x in r]


def unflat(v, n):
    return [list(v[i * n:(i + 1) * n]) for i in range(n)]


def lin_comb(coeffs, mats):
    n = len(mats[0])
    m = len(mats[0][0])
    C = zmat(n, m)
    for c, M in zip(coeffs, mats):
        if c == 0:
            continue
        for i in range(n):
            for j in range(m):
                C[i][j] += c * M[i][j]
    return C


def mprimitive(A):
    n = len(A)
    return unflat(primitive(flat(A)), n)


# ------------------------------------------------------------------ octonions
LINES = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def build_table(lines):
    T = [[None] * 8 for _ in range(8)]
    for i in range(8):
        T[0][i] = (1, i)
        T[i][0] = (1, i)
    for i in range(1, 8):
        T[i][i] = (-1, 0)
    seen = set()
    for (a, b, c) in lines:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            key = (min(x, y), max(x, y))
            if key in seen:
                raise ValueError("pair on two lines")
            seen.add(key)
            T[x][y] = (1, z)
            T[y][x] = (-1, z)
    if len(seen) != 21:
        raise ValueError("not a Fano plane")
    return T


def basis_vec(i, n=8):
    v = [0] * n
    v[i] = 1
    return v


E = [basis_vec(i) for i in range(8)]


def omul(T, x, y):
    r = [0] * 8
    for i in range(8):
        if x[i] == 0:
            continue
        for j in range(8):
            if y[j] == 0:
                continue
            s, k = T[i][j]
            r[k] += s * x[i] * y[j]
    return r


def associator(T, x, y, z):
    return vsub(omul(T, omul(T, x, y), z), omul(T, x, omul(T, y, z)))


def alternativity_scan(T):
    A = {}
    for i in range(8):
        for j in range(8):
            for k in range(8):
                A[(i, j, k)] = associator(T, E[i], E[j], E[k])
    fails = 0
    nonzero = 0
    nonzero_distinct_imag = 0
    for (i, j, k), v in A.items():
        if not is_zero_vec(v):
            nonzero += 1
            if 0 not in (i, j, k) and len({i, j, k}) == 3:
                nonzero_distinct_imag += 1
        if not is_zero_vec(vadd(v, A[(j, i, k)])):
            fails += 1
        if not is_zero_vec(vadd(v, A[(i, k, j)])):
            fails += 1
    return fails, nonzero, nonzero_distinct_imag, A


def norm2(x):
    return sum(t * t for t in x)


def der_equations(T):
    # unknown D[i][j] at index 8*i+j ; D(e_j) = sum_i D[i][j] e_i
    rows = []
    for a in range(8):
        for b in range(8):
            s_ab, k_ab = T[a][b]
            for m in range(8):
                r = [0] * 64
                r[8 * m + k_ab] += s_ab
                for i in range(8):
                    s, k = T[i][b]
                    if k == m:
                        r[8 * i + a] -= s
                    s, k = T[a][i]
                    if k == m:
                        r[8 * i + b] -= s
                rows.append(r)
    return rows


# ------------------------------------------------------------------ Gaussian-integer 3x3 matrices (re, im)
def c3_zero():
    return (zmat(3), zmat(3))


def c3_E(a, b, re=1, im=0):
    R = zmat(3)
    I = zmat(3)
    R[a][b] = re
    I[a][b] = im
    return (R, I)


def c3_add(X, Y):
    return (madd(X[0], Y[0]), madd(X[1], Y[1]))


def c3_sub(X, Y):
    return (msub(X[0], Y[0]), msub(X[1], Y[1]))


def c3_scale(k, X):
    return (mscale(k, X[0]), mscale(k, X[1]))


def c3_mul(X, Y):
    A, B = X
    C, D = Y
    return (msub(mmul(A, C), mmul(B, D)), madd(mmul(A, D), mmul(B, C)))


def c3_comm(X, Y):
    return c3_sub(c3_mul(X, Y), c3_mul(Y, X))


def c3_eq(X, Y):
    return meq(X[0], Y[0]) and meq(X[1], Y[1])


def c3_iszero(X):
    return mzero(X[0]) and mzero(X[1])


def c3_dagger(X):
    return (mT(X[0]), mneg(mT(X[1])))


def c3_in_su3(X):
    antiherm = c3_eq(c3_dagger(X), c3_scale(-1, X))
    trace0 = sum(X[0][i][i] for i in range(3)) == 0 and sum(X[1][i][i] for i in range(3)) == 0
    return antiherm and trace0


def c3_entry_i(X, sgn):
    # entry (a,b) multiplied by sgn(a,b) * i
    R = zmat(3)
    I = zmat(3)
    for a in range(3):
        for b in range(3):
            s = sgn(a, b)
            R[a][b] = -s * X[1][a][b]
            I[a][b] = s * X[0][a][b]
    return (R, I)


def c3_mask(X, keep):
    R = zmat(3)
    I = zmat(3)
    for a in range(3):
        for b in range(3):
            if keep(a, b):
                R[a][b] = X[0][a][b]
                I[a][b] = X[1][a][b]
    return (R, I)


def c3_diag_i(d):
    R = zmat(3)
    I = zmat(3)
    for k in range(3):
        I[k][k] = d[k]
    return (R, I)


def c3_offpair(a, b):
    # real basis of the root plane (a,b): E_ab - E_ba and i(E_ab + E_ba)
    X1 = c3_sub(c3_E(a, b), c3_E(b, a))
    X2 = c3_add(c3_E(a, b, 0, 1), c3_E(b, a, 0, 1))
    return [X1, X2]


def run_homog(name, hb, mb, J, proj, expect_N_zero):
    ok_su3 = all(c3_in_su3(X) for X in hb + mb)
    check(name + "_BASIS_IN_SU3", ok_su3, h_dim=len(hb), m_dim=len(mb))
    neg_one = all(c3_eq(J(J(X)), c3_scale(-1, X)) for X in mb)
    check(name + "_J_SQUARED_MINUS_ONE", neg_one, on="m_basis")
    pres = all(c3_eq(proj(J(X)), J(X)) and c3_in_su3(J(X)) for X in mb)
    check(name + "_J_PRESERVES_M", pres)
    red = all(c3_eq(proj(c3_comm(H, X)), c3_comm(H, X)) for H in hb for X in mb)
    check(name + "_REDUCTIVE_H_M_IN_M", red)
    hclosed = all(c3_iszero(proj(c3_comm(H1, H2))) for H1 in hb for H2 in hb)
    check(name + "_H_SUBALGEBRA", hclosed)
    eqv = all(c3_eq(J(c3_comm(H, X)), c3_comm(H, J(X))) for H in hb for X in mb)
    check(name + "_J_H_EQUIVARIANT", eqv)
    nz = 0
    mm_m_nonzero = 0
    total = 0
    for X in mb:
        for Y in mb:
            total += 1
            JX = J(X)
            JY = J(Y)
            t1 = proj(c3_comm(JX, JY))
            t2 = J(proj(c3_comm(JX, Y)))
            t3 = J(proj(c3_comm(X, JY)))
            t4 = proj(c3_comm(X, Y))
            N = c3_sub(c3_sub(c3_sub(t1, t2), t3), t4)
            if not c3_iszero(N):
                nz += 1
            if not c3_iszero(t4):
                mm_m_nonzero += 1
    if expect_N_zero:
        check(name + "_NIJENHUIS_ZERO", nz == 0, pairs=total, nonzero_pairs=nz,
              mm_bracket_m_part_nonzero_pairs=mm_m_nonzero, meaning="integrable")
    else:
        check(name + "_NIJENHUIS_NONZERO", nz > 0, pairs=total, nonzero_pairs=nz,
              mm_bracket_m_part_nonzero_pairs=mm_m_nonzero, meaning="non_integrable")
    return nz, mm_m_nonzero


def roots_closed(S):
    Sset = set(S)
    for (a, b) in S:
        for (c, d) in S:
            if b == c and a != d and (a, d) not in Sset:
                return False, ((a, b), (c, d), (a, d))
            if d == a and c != b and (c, b) not in Sset:
                return False, ((c, d), (a, b), (c, b))
    return True, None


# ------------------------------------------------------------------ main
def main():
    bad, slash_chars = self_lint()
    check("LINT_INTEGER_ONLY", len(bad) == 0 and slash_chars == 0,
          bad_tokens=len(bad), slash_chars_in_source=slash_chars,
          banned="division_ops,float_literals,float_or_fraction_or_decimal_or_complex_names,math_beyond_gcd")
    bad_src = ("x = 1 " + chr(47) + " 3\ny = 1" + chr(46) + "5\nimport math\nz = math.sqrt(2)\n").encode()
    cbad, cslash = lint_source(bad_src)
    control("CONTROL_LINT_CATCHES_BAD_SOURCE", len(cbad) == 0 and cslash == 0, False,
            sample="true_division+float_literal+math_sqrt", flagged_tokens=len(cbad))
    src_sha = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    emit("SCRIPT", name="m1_octonion.py", sha256=src_sha, python=sys.version.split()[0],
         E=0, network=0, float_used=0, claim="MEASURED")

    # ============================================================== 1. octonions
    T = build_table(LINES)
    emit("OCT_CONVENTION", basis="e0..e7_with_e0_unit",
         lines="(1,2,4)(2,3,5)(3,4,6)(4,5,7)(5,6,1)(6,7,2)(7,1,3)",
         rule="e_i*e_j:+e_k_for_(i,j,k)_cyclic_on_a_line;e_j*e_i:-e_k;e_i*e_i:-e0",
         index_form="(i,i+1,i+3)_mod_7", claim="NAMED",
         ref="Baez_2002_The_Octonions_Bull_AMS_39_145-205_(same_index_rule)")
    fails, nonzero, nz_imag, ASSOC = alternativity_scan(T)
    check("OCT_ASSOCIATOR_TOTALLY_ANTISYMMETRIC", fails == 0, basis_triples=512,
          transposition_tests=1024, failures=fails,
          implies="alternative_by_trilinearity")
    ex = ASSOC[(1, 2, 3)]
    check("OCT_NONASSOCIATIVE", nonzero > 0, nonzero_basis_associators=nonzero,
          nonzero_on_distinct_imaginary_triples=nz_imag, expect_distinct_imag=168,
          example_assoc_e1_e2_e3=fv(ex))
    check("OCT_NONASSOC_COUNT_168", nz_imag == 168, value=nz_imag,
          meaning="ordered_triples_of_distinct_units_off_a_common_line")
    nb_fail = 0
    for i in range(8):
        for j in range(8):
            if norm2(omul(T, E[i], E[j])) != norm2(E[i]) * norm2(E[j]):
                nb_fail += 1
    check("OCT_NORM_MULTIPLICATIVE_BASIS", nb_fail == 0, basis_pairs=64, failures=nb_fail)
    pol_fail = 0
    for i in range(8):
        for j in range(8):
            for k in range(8):
                for l in range(8):
                    lhs = dot(omul(T, E[i], E[j]), omul(T, E[k], E[l])) + \
                        dot(omul(T, E[i], E[l]), omul(T, E[k], E[j]))
                    rhs = 2 * (1 if i == k else 0) * (1 if j == l else 0)
                    if lhs != rhs:
                        pol_fail += 1
    check("OCT_COMPOSITION_POLARIZED", pol_fail == 0, basis_quads=4096, failures=pol_fail,
          identity="<xy,zw>+<xw,zy>:2<x,z><y,w>", implies="N(xy):N(x)N(y)_for_all_x,y")
    seed = 1000081
    rnd_fail = 0
    for t in range(300):
        vals = []
        for q in range(16):
            seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
            vals.append(((seed >> 33) % 11) - 5)
        x = vals[:8]
        y = vals[8:]
        if norm2(omul(T, x, y)) != norm2(x) * norm2(y):
            rnd_fail += 1
    check("OCT_NORM_MULTIPLICATIVE_RANDOM_INT", rnd_fail == 0, pairs=300, coeff_range="-5..5",
          failures=rnd_fail)
    # broken control: reverse orientation of one Fano line
    broken_lines = [(2, 1, 4)] + LINES[1:]
    TB = build_table(broken_lines)
    bf, bnz, _, _ = alternativity_scan(TB)
    control("CONTROL_BROKEN_FANO_ALTERNATIVITY", bf == 0, False,
            change="line_(1,2,4)_reversed_to_(2,1,4)", antisymmetry_failures=bf)

    # ============================================================== 2. der(O)
    DER_ROWS = der_equations(T)
    der_ns, der_rank = nullspace(DER_ROWS, 64)
    check("DER_DIM", len(der_ns) == 14, unknowns=64, equations=len(DER_ROWS), rank=der_rank,
          dim=len(der_ns), expect=14, method="fraction_free_rref_over_Z")
    sol_fail = sum(1 for x in der_ns for r in DER_ROWS if dot(r, x) != 0)
    check("DER_NULLSPACE_SATISFIES_SYSTEM", sol_fail == 0, failures=sol_fail)
    D_list = [unflat(x, 8) for x in der_ns]

    def in_der(M):
        x = flat(M)
        return all(dot(r, x) == 0 for r in DER_ROWS)

    leib_fail = 0
    for D in D_list:
        for a in range(8):
            for b in range(8):
                lhs = mvec(D, omul(T, E[a], E[b]))
                rhs = vadd(omul(T, mvec(D, E[a]), E[b]), omul(T, E[a], mvec(D, E[b])))
                if lhs != rhs:
                    leib_fail += 1
    check("DER_LEIBNIZ_DIRECT", leib_fail == 0, derivations=14, basis_pairs=64, failures=leib_fail)
    kill0 = all(D[i][0] == 0 for D in D_list for i in range(8))
    into_im = all(D[0][j] == 0 for D in D_list for j in range(8))
    check("DER_KILLS_E0", kill0 and into_im, kills_e0=kill0, image_in_ImO=into_im)
    antis = all(D[i][j] == -D[j][i] for D in D_list for i in range(1, 8) for j in range(1, 8))
    check("DER_ANTISYMMETRIC_ON_ImO", antis, meaning="der(O)_inside_so(7)")
    lie_fail = 0
    for a in range(14):
        for b in range(a + 1, 14):
            if not in_der(mcomm(D_list[a], D_list[b])):
                lie_fail += 1
    check("DER_LIE_CLOSED", lie_fail == 0, brackets=91, failures=lie_fail)
    G = [[-tr_prod(D_list[a], D_list[b]) for b in range(14)] for a in range(14)]
    minors = [det_bareiss([row[:k] for row in G[:k]]) for k in range(1, 15)]
    check("DER_TRACE_FORM_NEGATIVE_DEFINITE", all(m > 0 for m in minors),
          leading_minors_of_minus_trace_gram_positive=all(m > 0 for m in minors),
          meaning="compact_real_form_g2")

    # ============================================================== 3. h and m
    h_rows = DER_ROWS + [[1 if c == 8 * i + 7 else 0 for c in range(64)] for i in range(8)]
    h_ns, h_rank = nullspace(h_rows, 64)
    check("H_DIM", len(h_ns) == 8, def_="D_in_der(O)_with_D(e7):0", dim=len(h_ns), expect=8)
    H_list = [unflat(x, 8) for x in h_ns]
    hin = all(in_der(H) and all(H[i][7] == 0 for i in range(8)) for H in H_list)
    check("H_INSIDE_DER_AND_FIXES_E7", hin)
    hl_fail = 0
    for a in range(8):
        for b in range(a + 1, 8):
            C = mcomm(H_list[a], H_list[b])
            if not (in_der(C) and all(C[i][7] == 0 for i in range(8))):
                hl_fail += 1
    check("H_LIE_SUBALGEBRA", hl_fail == 0, brackets=28, failures=hl_fail)
    Gh = [[tr_prod(H_list[a], H_list[b]) for b in range(8)] for a in range(8)]
    dGh = det_bareiss(Gh)
    check("H_TRACE_FORM_NONDEGENERATE", dGh != 0, gram_det=dGh)
    Mtx = [[tr_prod(D_list[k], H_list[l]) for k in range(14)] for l in range(8)]
    m_ns, _ = nullspace(Mtx, 14)
    check("M_DIM", len(m_ns) == 6, def_="trace_orthogonal_complement_of_h_in_der(O)",
          dim=len(m_ns), expect=6)
    X_list = [mprimitive(lin_comb(c, D_list)) for c in m_ns]
    perp = all(tr_prod(X, H) == 0 for X in X_list for H in H_list)
    check("M_ORTHOGONAL_TO_H", perp)
    r14 = rank_of([flat(M) for M in H_list + X_list], 64)
    check("DER_EQUALS_H_PLUS_M", r14 == 14, rank_h_union_m=r14)
    P = [[X_list[c][r + 1][7] for c in range(6)] for r in range(6)]
    side = all(X[0][7] == 0 and X[7][7] == 0 for X in X_list)
    detP = det_bareiss(P)
    check("PHI_M_TO_TANGENT_ISOMORPHISM", detP != 0 and side,
          map="X->X(e7)_in_span(e1..e6)", det=detP, image_in_T=side)
    adj = adjugate(P)
    Y_list = [lin_comb([adj[k][a] for k in range(6)], X_list) for a in range(6)]
    g = 0
    for Y in Y_list:
        g = gcd(g, vgcd(flat(Y)))
    Y_list = [unflat([exact_div(x, g) for x in flat(Y)], 8) for Y in Y_list]
    d = exact_div(detP, g)
    if d < 0:
        Y_list = [mneg(Y) for Y in Y_list]
        d = -d

    def phi(Z):
        return [Z[r + 1][7] for r in range(6)]

    phiok = all(phi(Y_list[a]) == vscale(d, basis_vec(a, 6)) for a in range(6))
    yin = all(in_der(Y) and all(tr_prod(Y, H) == 0 for H in H_list) for Y in Y_list)
    check("M_ADAPTED_BASIS", phiok and yin, rule="Y_a_in_m_with_Y_a(e7):d*e_(a+1)", d=d)
    emit("M_BASIS_SHA256", sha256=hashlib.sha256(
        ",".join(str(x) for Y in Y_list for x in flat(Y)).encode()).hexdigest(),
        H_BASIS_SHA256=hashlib.sha256(",".join(str(x) for H in H_list for x in flat(H)).encode()).hexdigest(),
        claim="MEASURED")

    # ============================================================== 4. J on T
    Jmat = zmat(6)
    im_ok = True
    cross_ok = True
    for c in range(6):
        v = E[c + 1]
        w = omul(T, E[7], v)
        w2 = omul(T, v, E[7])
        if w[0] != 0 or w[7] != 0:
            im_ok = False
        if vscale(2, w) != vsub(w, w2):
            cross_ok = False
        for r in range(6):
            Jmat[r][c] = w[r + 1]
    jmap = []
    for c in range(6):
        col = [Jmat[r][c] for r in range(6)]
        for r in range(6):
            if col[r] != 0:
                jmap.append("e%d->%se%d" % (c + 1, "-" if col[r] < 0 else "", r + 1))
    check("J_DEFINED_ON_T", im_ok and cross_ok, rule="J(v):Im(e7*v):e7_cross_v",
          e7v_in_T=im_ok, equals_cross_product=cross_ok, action=";".join(jmap))
    I6 = ident(6)
    check("J_SQUARED_MINUS_ONE", meq(mmul(Jmat, Jmat), mneg(I6)))
    check("J_ORTHOGONAL", meq(mmul(mT(Jmat), Jmat), I6))
    HT = [[[H[r + 1][c + 1] for c in range(6)] for r in range(6)] for H in H_list]
    pres = all(H[0][c] == 0 and H[7][c] == 0 for H in H_list for c in range(1, 7))
    check("H_PRESERVES_T", pres)
    eq_fail = sum(1 for A in HT if not meq(mmul(A, Jmat), mmul(Jmat, A)))
    check("J_COMMUTES_WITH_H", eq_fail == 0, h_basis=8, failures=eq_fail,
          meaning="T_is_complex_3dim_rep_of_h")
    fidx = []
    for c in range(6):
        if c in fidx or any(Jmat[c][f] != 0 for f in fidx):
            continue
        fidx.append(c)
        if len(fidx) == 3:
            break
    ef = [basis_vec(f, 6) for f in fidx]
    Jef = [mvec(Jmat, v) for v in ef]
    ortho = all(dot(u, v) == (1 if i == j else 0)
                for i, u in enumerate(ef + Jef) for j, v in enumerate(ef + Jef))
    check("COMPLEX_UNITARY_BASIS", ortho,
          f="(" + ",".join("e%d" % (f + 1) for f in fidx) + ")",
          Jf="(" + ",".join("e%d" % (v.index(1) + 1) if 1 in v else "neg" for v in Jef) + ")")

    def ctrace(A):
        re = sum(dot(mvec(A, ef[k]), ef[k]) for k in range(3))
        im = sum(dot(mvec(A, ef[k]), Jef[k]) for k in range(3))
        return (re, im)

    traces = [ctrace(A) for A in HT]
    check("H_COMPLEX_TRACE_ZERO", all(t == (0, 0) for t in traces),
          gaussian_traces=";".join(fg(t) for t in traces),
          meaning="h_preserves_complex_volume_form;h_inside_su(3);dim8:su(3)")
    Jext = zmat(8)
    for r in range(6):
        for c in range(6):
            Jext[r + 1][c + 1] = Jmat[r][c]
    jtr = ctrace(Jmat)
    check("U1_GENERATOR_J_OUTSIDE_DER", (not in_der(Jext)) and jtr != (0, 0),
          J_complex_trace=fg(jtr), meaning="h_is_su(3)_not_u(3)")
    comm_rows = []
    for A in HT:
        for i in range(6):
            for j in range(6):
                r = [0] * 36
                for k in range(6):
                    r[6 * i + k] += A[k][j]
                    r[6 * k + j] -= A[i][k]
                comm_rows.append(r)
    cm_ns, _ = nullspace(comm_rows, 36)
    rI = rank_of(cm_ns + [flat(I6)], 36)
    rJ = rank_of(cm_ns + [flat(Jmat)], 36)
    check("H_COMMUTANT_ON_T", len(cm_ns) == 2 and rI == 2 and rJ == 2,
          dim=len(cm_ns), spanned_by="I,J",
          consequence="h_invariant_almost_complex_structures_are_exactly_+J_and_-J")

    # ============================================================== 5. Nijenhuis on m
    JmY = [lin_comb([Jmat[b][a] for b in range(6)], Y_list) for a in range(6)]
    jm_ok = all(phi(JmY[a]) == vscale(d, mvec(Jmat, basis_vec(a, 6))) and
                all(tr_prod(JmY[a], H) == 0 for H in H_list) for a in range(6))
    check("JM_TRANSPORTED", jm_ok, rule="J_m:phi^-1_J_phi")
    inter_fail = 0
    eqm_fail = 0
    for k, H in enumerate(H_list):
        for a in range(6):
            lhs = mcomm(H, Y_list[a])
            rhs = lin_comb([HT[k][b][a] for b in range(6)], Y_list)
            if not meq(lhs, rhs):
                inter_fail += 1
            JH = mmul(Jmat, HT[k])
            lhs2 = mcomm(H, JmY[a])
            rhs2 = lin_comb([JH[b][a] for b in range(6)], Y_list)
            if not meq(lhs2, rhs2):
                eqm_fail += 1
    check("M_REDUCTIVE_AND_PHI_INTERTWINES", inter_fail == 0, tests=48, failures=inter_fail,
          identity="[H,Y_a]:sum_b_H_T[b][a]Y_b")
    check("JM_H_EQUIVARIANT_ON_M", eqm_fail == 0, tests=48, failures=eqm_fail)
    H_flat = [flat(H) for H in H_list]
    proj_fail = 0
    span_fail = 0
    mm_h_nonzero = 0
    for a in range(6):
        for b in range(6):
            Z = mcomm(Y_list[a], Y_list[b])
            p = phi(Z)
            W = msub(mscale(d, Z), lin_comb(p, Y_list))
            if not (in_der(W) and all(W[i][7] == 0 for i in range(8)) and
                    all(tr_prod(W, Y) == 0 for Y in Y_list)):
                proj_fail += 1
            if rank_of(H_flat + [flat(W)], 64) != 8:
                span_fail += 1
            if not mzero(W):
                mm_h_nonzero += 1
    check("M_PROJECTION_VIA_PHI_EQUALS_TRACE_PROJECTION", proj_fail == 0 and span_fail == 0,
          brackets=36, failures=proj_fail, h_span_failures=span_fail,
          brackets_with_nonzero_h_part=mm_h_nonzero)

    def nij(Jm, JY):
        V = [[None] * 6 for _ in range(6)]
        for a in range(6):
            for b in range(6):
                X = Y_list[a]
                Yb = Y_list[b]
                JX = JY[a]
                JYb = JY[b]
                t1 = phi(mcomm(JX, JYb))
                t2 = mvec(Jm, phi(mcomm(JX, Yb)))
                t3 = mvec(Jm, phi(mcomm(X, JYb)))
                t4 = phi(mcomm(X, Yb))
                V[a][b] = vsub(vsub(vsub(t1, t2), t3), t4)
        return V

    V = nij(Jmat, JmY)
    nzp = sum(1 for a in range(6) for b in range(6) if not is_zero_vec(V[a][b]))
    mm_m_nz = sum(1 for a in range(6) for b in range(6)
                  if not is_zero_vec(phi(mcomm(Y_list[a], Y_list[b]))))
    check("NIJENHUIS_NONZERO", nzp > 0, pairs=36, nonzero_pairs=nzp,
          formula="N(X,Y):[JX,JY]_m-J[JX,Y]_m-J[X,JY]_m-[X,Y]_m",
          units="phi_coords_scale_d", mm_bracket_m_part_nonzero_pairs=mm_m_nz,
          meaning="G2_invariant_J_on_G2_mod_SU3_is_non_integrable")
    emit("NIJENHUIS_SAMPLE", N_e1_e2=fv(V[0][1]), N_e1_e3=fv(V[0][2]), N_e2_e4=fv(V[1][3]),
         N_e1_e4=fv(V[0][3]), basis="Y_a_with_phi(Y_a):d*e_(a+1)", d=d, claim="MEASURED")
    JmYn = [mneg(Y) for Y in JmY]
    Vn = nij(mneg(Jmat), JmYn)
    same = all(Vn[a][b] == V[a][b] for a in range(6) for b in range(6))
    check("NIJENHUIS_MINUS_J_EQUALS_N_J", same,
          consequence="no_G2_invariant_integrable_complex_structure_on_G2_mod_SU3")

    # Koszul: h_C + m10 closed under bracket ?
    def cpair_comm(P1, P2):
        A, B = P1
        C, D2 = P2
        return (msub(mcomm(A, C), mcomm(B, D2)), madd(mcomm(A, D2), mcomm(B, C)))

    Zs = [(Y_list[f], mneg(JmY[f])) for f in fidx]
    hz_fail = 0
    for H in H_list:
        for Zk in Zs:
            A, B = cpair_comm((H, zmat(8)), Zk)
            a = phi(A)
            b = phi(B)
            p01 = (vsub(a, mvec(Jmat, b)), vadd(b, mvec(Jmat, a)))
            hA = msub(mscale(d, A), lin_comb(a, Y_list))
            hB = msub(mscale(d, B), lin_comb(b, Y_list))
            if not (is_zero_vec(p01[0]) and is_zero_vec(p01[1]) and mzero(hA) and mzero(hB)):
                hz_fail += 1
    check("KOSZUL_H_BRACKET_M10_IN_M10", hz_fail == 0, tests=24, failures=hz_fail)
    closed = True
    rel_fail = 0
    for (j, k) in ((0, 1), (1, 2), (2, 0)):
        A, B = cpair_comm(Zs[j], Zs[k])
        a = phi(A)
        b = phi(B)
        two_p01 = (vsub(a, mvec(Jmat, b)), vadd(b, mvec(Jmat, a)))
        two_p10 = (vadd(a, mvec(Jmat, b)), vsub(b, mvec(Jmat, a)))
        hA = msub(mscale(d, A), lin_comb(a, Y_list))
        hB = msub(mscale(d, B), lin_comb(b, Y_list))
        p10_zero = is_zero_vec(two_p10[0]) and is_zero_vec(two_p10[1])
        p01_zero = is_zero_vec(two_p01[0]) and is_zero_vec(two_p01[1])
        h_zero = mzero(hA) and mzero(hB)
        if not p01_zero:
            closed = False
        Nv = V[fidx[j]][fidx[k]]
        if not (two_p01[0] == vscale(-1, Nv) and two_p01[1] == vscale(-1, mvec(Jmat, Nv))):
            rel_fail += 1
        pre = two_p01[0]
        mu = [(dot(pre, ef[l]), -dot(pre, Jef[l])) for l in range(3)]
        emit("KOSZUL_BRACKET", pair="[Z%d,Z%d]" % (j + 1, k + 1),
             m10_part_zero=p10_zero, m01_part_zero=p01_zero, hC_part_zero=h_zero,
             twice_m01_coeffs_on_conjZ1_conjZ2_conjZ3=";".join(fg(z) for z in mu),
             lands_in="m01" if (p10_zero and h_zero and not p01_zero) else "mixed",
             Zk="Y_fk-i*J_m*Y_fk", claim="MEASURED")
    check("KOSZUL_CLOSED_IS_FALSE", not closed, closed=closed,
          meaning="h_C+m10_not_a_subalgebra;Koszul_criterion_agrees_with_N_nonzero")
    check("KOSZUL_M01_PART_EQUALS_MINUS_N", rel_fail == 0, identity="2P01[Zj,Zk]:-(N+iJN)",
          failures=rel_fail)
    pure_fail = 0
    for (j, k) in ((0, 1), (1, 2), (2, 0)):
        A, B = cpair_comm(Zs[j], Zs[k])
        a = phi(A)
        b = phi(B)
        two_p10 = (vadd(a, mvec(Jmat, b)), vsub(b, mvec(Jmat, a)))
        hA = msub(mscale(d, A), lin_comb(a, Y_list))
        hB = msub(mscale(d, B), lin_comb(b, Y_list))
        if not (is_zero_vec(two_p10[0]) and is_zero_vec(two_p10[1]) and mzero(hA) and mzero(hB)):
            pure_fail += 1
    check("M10_BRACKET_LANDS_PURELY_IN_M01", pure_fail == 0, failures=pure_fail,
          meaning="[m10,m10]_in_m01:nearly_Kahler_type_3x3->3bar")

    # ============================================================== 6. structure of N
    anti_fail = sum(1 for a in range(6) for b in range(6) if V[a][b] != vscale(-1, V[b][a]))
    check("NIJENHUIS_ANTISYMMETRIC", anti_fail == 0, failures=anti_fail)

    def Vbil(u, w):
        r = [0] * 6
        for a in range(6):
            if u[a] == 0:
                continue
            for b in range(6):
                if w[b] == 0:
                    continue
                r = vadd(r, vscale(u[a] * w[b], V[a][b]))
        return r

    al1 = 0
    al2 = 0
    for a in range(6):
        for b in range(6):
            ea = basis_vec(a, 6)
            eb = basis_vec(b, 6)
            target = vscale(-1, mvec(Jmat, V[a][b]))
            if Vbil(mvec(Jmat, ea), eb) != target:
                al1 += 1
            if Vbil(ea, mvec(Jmat, eb)) != target:
                al2 += 1
    check("NIJENHUIS_ANTILINEAR_SLOT1", al1 == 0, identity="N(JX,Y):-J_N(X,Y)", failures=al1)
    check("NIJENHUIS_ANTILINEAR_SLOT2", al2 == 0, identity="N(X,JY):-J_N(X,Y)", failures=al2)
    inv_fail = 0
    for A in HT:
        for a in range(6):
            for b in range(6):
                lhs = mvec(A, V[a][b])
                rhs = vadd(Vbil(mvec(A, basis_vec(a, 6)), basis_vec(b, 6)),
                           Vbil(basis_vec(a, 6), mvec(A, basis_vec(b, 6))))
                if lhs != rhs:
                    inv_fail += 1
    check("NIJENHUIS_SU3_INVARIANT", inv_fail == 0, tests=288, failures=inv_fail)

    def zc(t):
        return [(dot(t, ef[k]), dot(t, Jef[k])) for k in range(3)]

    lam = zc(V[fidx[0]][fidx[1]])[2]
    zero = (0, 0)
    e12 = zc(V[fidx[0]][fidx[1]]) == [zero, zero, lam]
    e23 = zc(V[fidx[1]][fidx[2]]) == [lam, zero, zero]
    e31 = zc(V[fidx[2]][fidx[0]]) == [zero, lam, zero]
    check("NIJENHUIS_EPSILON_TYPE", e12 and e23 and e31 and lam != zero,
          law="N(f_a,f_b):lambda*eps_abc*f_c", lam=fg(lam), scale="phi_coords_of_d-scaled_basis",
          N_f1f2=";".join(fg(z) for z in zc(V[fidx[0]][fidx[1]])),
          N_f2f3=";".join(fg(z) for z in zc(V[fidx[1]][fidx[2]])),
          N_f3f1=";".join(fg(z) for z in zc(V[fidx[2]][fidx[0]])),
          meaning="Lambda2(3bar)->3_via_eps_abc;N_antilinear_each_slot")
    # probe: N against the octonion cross product (not a pass or fail check)
    Wc = [[[omul(T, E[a + 1], E[b + 1])[r + 1] for r in range(6)] for b in range(6)] for a in range(6)]
    for label, cand in (("pi_T(u_x_v)", Wc),
                        ("J_pi_T(u_x_v)", [[mvec(Jmat, Wc[a][b]) for b in range(6)] for a in range(6)])):
        ref = None
        for a in range(6):
            for b in range(6):
                for r in range(6):
                    if ref is None and cand[a][b][r] != 0:
                        ref = (cand[a][b][r], V[a][b][r])
        prop = ref is not None and all(
            V[a][b][r] * ref[0] == cand[a][b][r] * ref[1]
            for a in range(6) for b in range(6) for r in range(6))
        if prop:
            gg = gcd(ref[0], ref[1])
            ratio = "%d:%d" % (exact_div(ref[1], gg), exact_div(ref[0], gg))
        else:
            ratio = "none"
        emit("PROBE_N_VS_CROSS_PRODUCT", candidate=label, proportional=prop,
             ratio_N_to_candidate=ratio, claim="MEASURED")

    # ============================================================== 7. controls
    # 7a: CP2 = SU3 mod S(U2xU1), symmetric space
    blk = [0, 0, 1]
    cp_h = [c3_diag_i([1, -1, 0]), c3_diag_i([1, 1, -2])] + c3_offpair(0, 1)
    cp_m = c3_offpair(0, 2) + c3_offpair(1, 2)

    def cp_sgn(a, b):
        if blk[a] == blk[b]:
            return 0
        return 1 if blk[a] < blk[b] else -1

    def cp_J(X):
        return c3_entry_i(X, cp_sgn)

    def cp_proj(X):
        return c3_mask(X, lambda a, b: blk[a] != blk[b])

    Zc = c3_diag_i([1, 1, -2])
    ad_ok = True
    for X in cp_m:
        C = c3_comm(Zc, X)
        try:
            CJ = (unflat([exact_div(x, 3) for x in flat(C[0])], 3),
                  unflat([exact_div(x, 3) for x in flat(C[1])], 3))
        except ArithmeticError:
            ad_ok = False
            continue
        if not c3_eq(CJ, cp_J(X)):
            ad_ok = False
    check("CP2_J_EQUALS_AD_CENTRE", ad_ok, centre="i*diag(1,1,-2)", normalisation="divide_by_3_exact")
    cp_nz, cp_mm = run_homog("CP2", cp_h, cp_m, cp_J, cp_proj, True)
    emit("CP2_NOTE", symmetric_space=(cp_mm == 0), mm_bracket_m_part_nonzero_pairs=cp_mm,
         note="[m,m]_in_h_so_N_vanishes_termwise", claim="MEASURED")
    S_cp = [(0, 1), (1, 0)] + [(a, b) for a in range(3) for b in range(3) if cp_sgn(a, b) == 1]
    cl_cp, _ = roots_closed(S_cp)
    check("CP2_KOSZUL_ROOTS_CLOSED", cl_cp and cp_nz == 0, roots_hC_plus_m10=fv(S_cp),
          closed=cl_cp, agrees_with_N=(cl_cp == (cp_nz == 0)))

    # 7b: flag manifold SU3 mod T2, non-symmetric: integrable J vs nearly-Kahler J
    fl_h = [c3_diag_i([1, -1, 0]), c3_diag_i([0, 1, -1])]
    fl_m = c3_offpair(0, 1) + c3_offpair(0, 2) + c3_offpair(1, 2)

    def fl_proj(X):
        return c3_mask(X, lambda a, b: a != b)

    for label, eps, expect_zero in (("FLAG_INTEGRABLE", {(0, 1): 1, (0, 2): 1, (1, 2): 1}, True),
                                    ("FLAG_NEARLY_KAHLER", {(0, 1): 1, (0, 2): -1, (1, 2): 1}, False)):
        def sg(a, b, eps=eps):
            if a == b:
                return 0
            return eps[(a, b)] if a < b else -eps[(b, a)]

        def Jf(X, sg=sg):
            return c3_entry_i(X, sg)

        fnz, fmm = run_homog(label, fl_h, fl_m, Jf, fl_proj, expect_zero)
        S = [(a, b) for a in range(3) for b in range(3) if a != b and sg(a, b) == 1]
        cl, wit = roots_closed(S)
        check(label + "_KOSZUL_ROOTS", cl == expect_zero and cl == (fnz == 0),
              eps="(%d,%d,%d)" % (eps[(0, 1)], eps[(0, 2)], eps[(1, 2)]),
              roots_m10=fv(S), closed=cl, witness=str(wit), mm_bracket_m_part_nonzero_pairs=fmm,
              agrees_with_N=(cl == (fnz == 0)))

    # 7c: broken J on T (coordinates e2 and e3 swapped), must lose h-equivariance
    S6 = ident(6)
    S6[1][1] = 0
    S6[2][2] = 0
    S6[1][2] = 1
    S6[2][1] = 1
    Jb = mmul(mmul(S6, Jmat), S6)
    check("BROKEN_J_STILL_SQUARES_TO_MINUS_ONE", meq(mmul(Jb, Jb), mneg(I6)))
    check("BROKEN_J_STILL_ORTHOGONAL", meq(mmul(mT(Jb), Jb), I6))
    bfail = sum(1 for A in HT if not meq(mmul(A, Jb), mmul(Jb, A)))
    control("CONTROL_BROKEN_J_EQUIVARIANCE", bfail == 0, False,
            change="J_conjugated_by_swap(e2,e3)", h_basis_elements_failing=bfail)

    # ============================================================== conclusions
    emit("CONCLUSION", claim="MEASURED",
         text="der(O)_dim14_compact;stabiliser_of_e7_dim8_is_su(3);m_dim6_iso_T_e7(S6);"
              "J:e7_cross_is_the_G2_invariant_almost_complex_structure_(unique_up_to_sign);"
              "Nijenhuis_N_nonzero_on_%d_of_36_pairs;N_is_SU3_invariant_antilinear_eps_type_lambda_%s;"
              "Koszul_closure_false;[m10,m10]_lands_purely_in_m01" % (nzp, fg(lam)))
    emit("CONCLUSION", claim="MEASURED",
         text="controls:CP2_N:0_integrable(symmetric);flag_SU3_mod_T2_eps(+,+,+)_N:0;"
              "flag_eps(+,-,+)_N_nonzero;broken_Fano_and_broken_J_fail_as_designed")
    emit("NAMED", claim="NAMED", source="cited_from_memory_not_fetched_E0",
         ref="Kirchhoff_1947_CRAS_225:octonionic_almost_complex_structure_on_S6;"
             "Ehresmann-Libermann_1951_CRAS_232:it_is_not_integrable;"
             "Borel-Serre_1953_AmerJMath_75:S2_and_S6_are_the_only_almost_complex_spheres;"
             "Newlander-Nirenberg_1957_AnnMath_65:integrable_iff_N:0;"
             "Koszul_1955_CanadJMath_7_and_Frolicher_1955_MathAnn_129:invariant_complex_structures_on_G_mod_H;"
             "Kobayashi-Nomizu_vol2_1969_ch10:Nijenhuis_on_reductive_homogeneous_spaces;"
             "Gray_1976_MathAnn_223:nearly_Kahler_S6")
    emit("OPEN", claim="CONJECTURE", status="OPEN",
         text="whether_S6_carries_any_integrable_complex_structure_(not_G2_invariant)_is_the_open_Hopf_problem")
    emit("BRIDGE_NOTE", claim="NAMED",
         text="eps_abc_contraction_3x3->3bar_measured_here_is_the_same_SU3_invariant_as_the_colour_singlet_epsilon;"
              "tie_to_4D_and_finite_Yang-Mills_mass_gap_is_not_computed_in_this_step")



def finalize():
    verdict = "ALL_PASS" if STATE["fail"] == 0 else "SOME_FAIL"
    emit("SUMMARY", checks=STATE["checks"], fail=STATE["fail"], verdict=verdict,
         failed=",".join(STATE["failed_tags"]) if STATE["failed_tags"] else "none", claim="MEASURED")

    body = ("\n".join(ROWS) + "\n").encode()
    with open(OUT, "wb") as fh:
        fh.write(body)
    digest = hashlib.sha256(body).hexdigest()
    with open(SHA_OUT, "wb") as fh:
        fh.write(("%s  %s\n" % (digest, OUT_NAME)).encode())
    for r in ROWS:
        print(r)
    print("HBP_SHA256|sha256=%s|file=%s|json=0" % (digest, OUT_NAME))
    return 0 if STATE["fail"] == 0 else 1


def run():
    try:
        main()
    except Exception as ex:
        STATE["checks"] += 1
        STATE["fail"] += 1
        STATE["failed_tags"].append("CRASH")
        tb = ex.__traceback__
        line = 0
        while tb is not None:
            line = tb.tb_lineno
            tb = tb.tb_next
        emit("CRASH", result="FAIL", error=type(ex).__name__, detail=str(ex)[:120],
             line=line, meaning="a_stage_could_not_complete;downstream_checks_not_run",
             claim="MEASURED")
    return finalize()


if __name__ == "__main__":
    sys.exit(run())
