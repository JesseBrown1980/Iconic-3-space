#!/usr/bin/env python3
# m2_roots.py
# G2 root system in Z^3 on the plane x+y+z=0, W(G2) as 3x3 integer matrices,
# the long-root A2 (the roots of su(3) inside g2), der(O) rebuilt from an explicit
# octonion table, the Higgs breaking G2 -> SU(3) by v = e7 (S^6 is G2/SU(3)),
# the invariant complex structure J = L_e7 on the tangent space, the root system
# read off the actual algebra, the e7 -> -e7 automorphism and its Weyl element,
# and broken controls that must FAIL.
#
# Integer arithmetic only: no floating point, no true division, no rational or
# decimal modules; the math module is used for gcd only. E=0: compute only,
# no network; the only writes are M2-ROOTS.hbp and its .sha256 sidecar here.
# Run: python3 -I -B m2_roots.py
import hashlib
import io
import itertools
import os
import sys
import tokenize
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "M2-ROOTS.hbp"
OUT = os.path.join(HERE, OUT_NAME)
ROWS = []
FAILS = []
COUNTS = {"check": 0, "control": 0}


# ---------------------------------------------------------------- output
def fv(v):
    return ",".join(str(x) for x in v)


def fvl(vs):
    return ";".join(fv(v) for v in vs)


def row(tag, **kv):
    for k in kv:
        s = str(kv[k])
        if "|" in s or "=" in s:
            raise ValueError("separator inside value of " + k)
    line = "|".join([tag] + ["%s=%s" % (k, kv[k]) for k in kv] + ["json=0"])
    ROWS.append(line)
    print(line)


def check(name, cond, **kv):
    COUNTS["check"] += 1
    kv.setdefault("claim", "MEASURED")
    if cond:
        row("CHECK", name=name, status="PASS", **kv)
    else:
        row("FAIL", name=name, status="FAIL", **kv)
        FAILS.append(name)
    return cond


def control(name, broken_passed, **kv):
    # a deliberately broken input: its check must FAIL; if it passes, that is a FAIL
    COUNTS["control"] += 1
    if broken_passed:
        row("FAIL", name=name, status="CONTROL_DID_NOT_FAIL", **kv)
        FAILS.append(name)
    else:
        row("CONTROL", name=name, expected="FAIL", got="FAIL", ok=1, **kv)


# ---------------------------------------------------------------- integer linear algebra
def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def mm(A, B):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(r, c)) for c in Bt] for r in A]


def mv(A, v):
    return [sum(a * b for a, b in zip(r, v)) for r in A]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def tp(A):
    return [list(r) for r in zip(*A)]


def eye(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def neg(A):
    return [[-x for x in r] for r in A]


def madd(A, B, s=1):
    return [[a + s * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def smul(k, A):
    return [[k * x for x in r] for r in A]


def comm(A, B):
    return madd(mm(A, B), mm(B, A), -1)


def flat(A):
    return [x for r in A for x in r]


def key(A):
    return tuple(tuple(r) for r in A)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def exdiv(a, b):
    q, r = divmod(a, b)
    if r != 0:
        raise ArithmeticError("inexact division %d by %d" % (a, b))
    return q


def prim(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    if g > 1:
        v = [exdiv(x, g) for x in v]
    return v


def rref(rows_in):
    # fraction-free reduced row echelon form over Z (each row kept primitive)
    M = [list(r) for r in rows_in if any(r)]
    if not M:
        return [], []
    ncol = len(M[0])
    piv = []
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        M[r] = prim(M[r])
        if M[r][c] < 0:
            M[r] = [-x for x in M[r]]
        a = M[r][c]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                b = M[i][c]
                g = gcd(a, b)
                aa, bb = exdiv(a, g), exdiv(b, g)
                M[i] = prim([aa * x - bb * y for x, y in zip(M[i], M[r])])
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return M[:r], piv


def rank(rows_in):
    return len(rref(rows_in)[1])


def kernel(rows_in, ncol):
    R, piv = rref(rows_in)
    free = [c for c in range(ncol) if c not in piv]
    L = 1
    for i, c in enumerate(piv):
        a = R[i][c]
        L = exdiv(L * a, gcd(L, a))
    basis = []
    for f in free:
        x = [0] * ncol
        x[f] = L
        for i, c in enumerate(piv):
            x[c] = -exdiv(L * R[i][f], R[i][c])
        basis.append(prim(x))
    for x in basis:
        for r in rows_in:
            if dot(r, x) != 0:
                raise ArithmeticError("kernel vector fails a row")
    return basis


def det(A):
    # Bareiss fraction-free determinant
    n = len(A)
    if n == 0:
        return 1
    M = [list(r) for r in A]
    sgn = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sgn = -sgn
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = exdiv(M[i][j] * M[k][k] - M[i][k] * M[k][j], prev)
        prev = M[k][k]
    return sgn * M[n - 1][n - 1]


def minor(A, i, j):
    return [r[:j] + r[j + 1:] for k, r in enumerate(A) if k != i]


def adjugate(A):
    n = len(A)
    return [[(1 if (i + j) % 2 == 0 else -1) * det(minor(A, j, i)) for j in range(n)] for i in range(n)]


def gen_group(gens):
    n = len(gens[0])
    I = eye(n)
    seen = {key(I): I}
    frontier = [I]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                p = mm(g, h)
                k = key(p)
                if k not in seen:
                    seen[k] = p
                    nxt.append(p)
        frontier = nxt
        if len(seen) > 100000:
            raise RuntimeError("group too large")
    return seen


def order(g):
    I = eye(len(g))
    p = g
    k = 1
    while p != I:
        p = mm(p, g)
        k += 1
        if k > 1000:
            return 0
    return k


def commutant(gens, n):
    eqs = []
    for H in gens:
        for r in range(n):
            for c in range(n):
                e = [0] * (n * n)
                for k in range(n):
                    e[r * n + k] += H[k][c]
                    e[k * n + c] -= H[r][k]
                eqs.append(e)
    return kernel(eqs, n * n)


def unflat(v, n):
    return [v[i * n:(i + 1) * n] for i in range(n)]


def lincomb(cs, mats):
    n = len(mats[0])
    out = [[0] * n for _ in range(n)]
    for c, A in zip(cs, mats):
        if c:
            out = madd(out, A, c)
    return out


def reduced_pair(num, den):
    if den < 0:
        num, den = -num, -den
    g = gcd(num, den)
    return exdiv(num, g), exdiv(den, g)


# ---------------------------------------------------------------- self no-float gate
ALLOWED_IMPORTS = {"hashlib", "io", "itertools", "os", "sys", "tokenize"}
BANNED_NAMES = {"float", "complex", "Fraction", "Decimal", "fractions", "decimal"}


def audit(src):
    bad = []
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type == tokenize.OP and t.string in ("/", "/="):
            bad.append("truediv_line%d" % t.start[0])
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if not s.startswith(("0x", "0o", "0b")) and any(ch in s for ch in ".ej"):
                bad.append("floatlit_line%d" % t.start[0])
        elif t.type == tokenize.NAME and t.string in BANNED_NAMES:
            bad.append("name_%s_line%d" % (t.string, t.start[0]))
        elif t.type == tokenize.NAME and t.start[1] == 0 and t.string in ("import", "from"):
            parts = t.line.replace(",", " ").split()
            if parts[0] == "import":
                for mod in parts[1:]:
                    if mod not in ALLOWED_IMPORTS:
                        bad.append("import_" + mod)
            elif parts[0] == "from":
                if parts[1] != "math" or any(nm not in ("gcd", "isqrt") for nm in parts[3:]):
                    bad.append("from_" + parts[1])
    return bad


def part0_gate():
    with open(os.path.abspath(__file__), "r", encoding="utf-8") as fh:
        src = fh.read()
    bad = audit(src)
    check("self_gate_no_float_no_truediv_imports_ok", bad == [], findings=len(bad),
          detail=(fv(bad) if bad else "none"))
    broken = "x " + "/" + " 2\ny " + "=" + " 1" + "." + "5\nfrom fractions import Fraction\n"
    control("ctl_self_gate_catches_truediv_floatlit_fractions", audit(broken) == [],
            caught=len(audit(broken)))
    row("META", script="m2_roots.py", python=fv(sys.version_info[:3]), e=0, float_used=0,
        int_only=1, flags="-I -B", out=OUT_NAME)


# ---------------------------------------------------------------- PART 1: roots
SHORT = [(1, -1, 0), (1, 0, -1), (0, 1, -1)]
LONG = [(2, -1, -1), (-1, 2, -1), (-1, -1, 2)]


def pm(vs):
    return [tuple(v) for v in vs] + [tuple(-x for x in v) for v in vs]


R_SHORT = pm(SHORT)
R_LONG = pm(LONG)
ROOTS = R_SHORT + R_LONG
B1, B2 = (1, -1, 0), (0, 1, -1)


def refl(v, a):
    q, r = divmod(2 * dot(v, a), dot(a, a))
    return r == 0, tuple(x - q * y for x, y in zip(v, a)), q


def rs_checks(roots):
    S = set(roots)
    plane = all(sum(v) == 0 for v in roots)
    integral = True
    closed = True
    cart = set()
    for a in roots:
        for b in roots:
            ok, img, q = refl(b, a)
            if not ok:
                integral = False
                continue
            cart.add(q)
            if img not in S:
                closed = False
    return plane, integral, closed, cart


def part1():
    row("G2_ROOTS", n=len(ROOTS), short=fvl(R_SHORT), long=fvl(R_LONG), claim="MEASURED")
    check("roots_distinct_12", len(set(ROOTS)) == 12, n=len(set(ROOTS)))
    check("short_norm2_x6", sorted(dot(v, v) for v in R_SHORT) == [2] * 6)
    check("long_norm6_x6", sorted(dot(v, v) for v in R_LONG) == [6] * 6)
    nmax = max(dot(v, v) for v in ROOTS)
    nmin = min(dot(v, v) for v in ROOTS)
    check("crystallographic_length_ratio_3", nmax == 3 * nmin, long_sq=nmax, short_sq=nmin)
    plane, integral, closed, cart = rs_checks(ROOTS)
    check("roots_in_plane_x_plus_y_plus_z_0", plane)
    check("reflection_coefficient_exact_integer_all_pairs", integral, pairs=len(ROOTS) * len(ROOTS))
    check("closed_under_all_12_reflections", closed)
    check("cartan_integers_in_range_and_3_occurs", cart <= set(range(-3, 4)) and 3 in cart and -3 in cart,
          values=fv(sorted(cart)))
    box = [v for v in itertools.product(range(-3, 4), repeat=3) if sum(v) == 0]
    check("short_roots_are_all_norm2_plane_lattice_vectors", set(v for v in box if dot(v, v) == 2) == set(R_SHORT))
    check("long_roots_are_all_norm6_plane_lattice_vectors", set(v for v in box if dot(v, v) == 6) == set(R_LONG))
    a1, a2 = (1, -1, 0), (-2, 1, 1)
    cart_m = [[exdiv(2 * dot(x, y), dot(y, y)) for y in (a1, a2)] for x in (a1, a2)]
    check("cartan_matrix_G2_simple_a1short_a2long", cart_m == [[2, -1], [-3, 2]], A=fvl(cart_m))
    coef = {v: (v[2] - v[1], v[2]) for v in ROOTS}
    recon = all(tuple(c[0] * a1[i] + c[1] * a2[i] for i in range(3)) == v for v, c in coef.items())
    check("root_coords_in_simple_basis_exact", recon)
    check("every_root_has_same_sign_simple_coeffs",
          all((c[0] >= 0 and c[1] >= 0) or (c[0] <= 0 and c[1] <= 0) for c in coef.values()))
    pos = sorted(c for c in coef.values() if c[0] >= 0 and c[1] >= 0)
    check("positive_roots_G2", pos == [(0, 1), (1, 0), (1, 1), (2, 1), (3, 1), (3, 2)], pos=fvl(pos))
    hi = tuple(3 * a1[i] + 2 * a2[i] for i in range(3))
    check("highest_root_3a1_2a2_is_long", hi in R_LONG, highest=fv(hi))


# ---------------------------------------------------------------- PART 2: W(G2)
ONE3 = (1, 1, 1)


def lift3(a, eps):
    # 3 * M_eps, where M_eps agrees with s_a on the plane and sends (1,1,1) to eps*(1,1,1)
    k = exdiv(6, dot(a, a))
    return [[3 * (1 if i == j else 0) - k * a[i] * a[j] + (eps - 1) for j in range(3)] for i in range(3)]


def plane2(M):
    def pc(v):
        if sum(v) != 0:
            raise ArithmeticError("vector left the plane")
        return (v[0], -v[2])
    c1 = pc(mv(M, B1))
    c2 = pc(mv(M, B2))
    return [[c1[0], c2[0]], [c1[1], c2[1]]]


def det2(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def underlying_perm(M):
    n = len(M)
    sig = [None] * n
    for i, r in enumerate(M):
        nz = [j for j, x in enumerate(r) if x != 0]
        if len(nz) != 1 or r[nz[0]] not in (1, -1):
            return None
        sig[nz[0]] = i
    if None in sig:
        return None
    return tuple(sig)


def part2():
    LIFT = {}
    EPS = {}
    one_lift = True
    is_refl = True
    plus_integral_long = sum(1 for a in R_LONG if all(x % 3 == 0 for r in lift3(a, 1) for x in r))
    for a in ROOTS:
        found = []
        for eps in (1, -1):
            T3 = lift3(a, eps)
            if all(x % 3 == 0 for r in T3 for x in r):
                found.append((eps, [[exdiv(x, 3) for x in r] for r in T3]))
        if len(found) != 1:
            one_lift = False
            continue
        eps, M = found[0]
        LIFT[a] = M
        EPS[a] = eps
        ok = mm(tp(M), M) == eye(3)
        ok = ok and all(tuple(mv(M, v)) == refl(v, a)[1] for v in (B1, B2))
        ok = ok and mv(M, ONE3) == [eps] * 3
        ok = ok and det2(plane2(M)) == -1
        is_refl = is_refl and ok
    check("exactly_one_orthogonal_integer_lift_per_root", one_lift)
    check("lift_equals_reflection_on_plane_det_minus1", is_refl)
    check("short_lifts_fix_111_long_lifts_negate_111",
          all(EPS[a] == 1 for a in R_SHORT) and all(EPS[a] == -1 for a in R_LONG))
    row("LIFT", convention="M agrees with s_a on x+y+z plane; M(1,1,1) is eps(1,1,1); the R3 reflection of a long root has entries in thirds so its integer lift has eps -1",
        s_short_example=fvl(LIFT[(1, -1, 0)]), s_long_example=fvl(LIFT[(2, -1, -1)]), claim="MEASURED")
    control("ctl_R3_reflection_of_long_root_is_integral", plus_integral_long > 0, integral_count=plus_integral_long)

    W = gen_group([LIFT[a] for a in ROOTS])
    Wl = list(W.values())
    check("W_G2_order_12", len(Wl) == 12, order=len(Wl))
    hist = {}
    for g in Wl:
        o = order(g)
        hist[o] = hist.get(o, 0) + 1
    hs = sorted(hist.items())
    check("W_element_orders_match_D6", hs == [(1, 1), (2, 7), (3, 2), (6, 2)], orders=fvl(hs))
    center = [g for g in Wl if all(mm(g, h) == mm(h, g) for h in Wl)]
    I3 = eye(3)
    mI3 = neg(I3)
    check("W_center_is_I_and_minus_I", sorted(key(c) for c in center) == sorted([key(I3), key(mI3)]), size=len(center))
    check("minus_I_in_W_G2", key(mI3) in W)
    check("W_nonabelian", len(center) < len(Wl))
    r6 = [g for g in Wl if order(g) == 6][0]
    cyc = {}
    p = eye(3)
    for _ in range(6):
        cyc[key(p)] = p
        p = mm(p, r6)
    rinv = tp(r6)
    s2 = [g for g in Wl if order(g) == 2 and key(g) not in cyc][0]
    check("dihedral_relation_s_r_s_is_r_inverse", mm(mm(s2, r6), s2) == rinv and mm(r6, rinv) == I3,
          r=fvl(r6), s=fvl(s2))
    check("minus_I_is_r_cubed", mm(mm(r6, r6), r6) == mI3)
    p2 = {key(plane2(g)) for g in Wl}
    check("W_acts_faithfully_on_plane", len(p2) == 12)
    nrefl = sum(1 for g in Wl if det2(plane2(g)) == -1)
    nrot = sum(1 for g in Wl if det2(plane2(g)) == 1)
    check("W_has_6_plane_reflections_6_rotations", nrefl == 6 and nrot == 6)
    check("all_W_elements_are_signed_permutations", all(underlying_perm(g) is not None for g in Wl))
    check("W_is_all_plus_minus_permutation_matrices",
          {key(g) for g in Wl} == {key(smul(s, [[1 if sig[j] == i else 0 for j in range(3)] for i in range(3)]))
                                    for sig in itertools.permutations(range(3)) for s in (1, -1)})
    row("W_G2", order=12, iso="D6 dihedral of order 12 (also S3 x Z2)", order_hist=fvl(hs),
        minus_I_plane=fvl(plane2(mI3)), claim="MEASURED")
    row("NAMED", statement="W(G2) is dihedral of order 12; G2 roots in x+y+z plane coordinates",
        ref="Humphreys GTM9 1972 sec 12.1; Bourbaki Lie Ch VI Plate IX", claim="NAMED")
    return LIFT, W


# ---------------------------------------------------------------- PART 3: A2 inside
def part3(LIFT, W):
    Wl = list(W.values())
    I3 = eye(3)
    mI3 = neg(I3)
    b1, b2 = (2, -1, -1), (-1, 2, -1)
    cart = [[exdiv(2 * dot(x, y), dot(y, y)) for y in (b1, b2)] for x in (b1, b2)]
    check("long_roots_cartan_matrix_A2", cart == [[2, -1], [-1, 2]], A=fvl(cart))
    plane, integral, closed, cv = rs_checks(R_LONG)
    check("long_roots_closed_root_system_no_3", plane and integral and closed and cv <= set(range(-2, 3)),
          values=fv(sorted(cv)))
    S = set(ROOTS)
    long_sum_pairs = [(a, b) for a in R_LONG for b in R_LONG if tuple(x + y for x, y in zip(a, b)) in S]
    check("long_roots_closed_under_addition_su3_subalgebra",
          all(tuple(x + y for x, y in zip(a, b)) in set(R_LONG) for a, b in long_sum_pairs),
          pairs=len(long_sum_pairs))
    short_to_long = [(a, b) for a in R_SHORT for b in R_SHORT if tuple(x + y for x, y in zip(a, b)) in set(R_LONG)]
    check("short_roots_not_closed_short_plus_short_is_long", len(short_to_long) > 0, pairs=len(short_to_long))
    Lg = gen_group([LIFT[a] for a in R_LONG])
    Ll = list(Lg.values())
    check("W_A2_long_subgroup_order_6", len(Ll) == 6, order=len(Ll))
    perms = [underlying_perm(g) for g in Ll]
    check("W_A2_acts_by_signed_permutations_of_axes", all(p is not None for p in perms))
    hom = all(underlying_perm(mm(g, h)) == tuple(underlying_perm(g)[underlying_perm(h)[j]] for j in range(3))
              for g in Ll for h in Ll)
    check("axis_permutation_map_is_isomorphism_onto_S3", hom and len(set(perms)) == 6)
    def pmat(sig):
        return [[1 if sig[j] == i else 0 for j in range(3)] for i in range(3)]
    check("W_A2_is_sign_of_sigma_times_P_sigma",
          all(g == smul(det(pmat(underlying_perm(g))), pmat(underlying_perm(g))) for g in Ll))
    check("minus_I_not_in_W_A2", key(mI3) not in Lg)
    normal = all(key(mm(mm(g, l), tp(g))) in Lg for g in Wl for l in Ll)
    check("W_A2_normal_in_W_G2", normal)
    check("index_2", len(Wl) == 2 * len(Ll))
    coset = {key(neg(l)) for l in Ll}
    check("W_G2_is_W_A2_times_plus_minus_I", set(Lg) | coset == set(W) and not (set(Lg) & coset))
    Sg = gen_group([LIFT[a] for a in R_SHORT])
    check("short_root_subgroup_is_plain_S3_order_6", len(Sg) == 6 and all(min(flat(g)) >= 0 for g in Sg.values()))
    # short roots under W(A2)
    orbits = []
    seen = set()
    for v in R_SHORT:
        if v in seen:
            continue
        orb = sorted({tuple(mv(g, v)) for g in Ll})
        seen |= set(orb)
        orbits.append(orb)
    check("short_roots_two_W_A2_orbits_of_size_3", sorted(len(o) for o in orbits) == [3, 3])
    box = [v for v in itertools.product(range(-3, 4), repeat=3) if sum(v) == 0]
    w1 = [v for v in box if dot(v, b1) == 3 and dot(v, b2) == 0]
    w2 = [v for v in box if dot(v, b1) == 0 and dot(v, b2) == 3]
    check("fundamental_weights_unique_and_short_roots", len(w1) == 1 and len(w2) == 1 and w1[0] in S and w2[0] in S,
          omega1=fvl(w1), omega2=fvl(w2))
    om1, om2 = w1[0], w2[0]
    O3 = [o for o in orbits if om1 in o][0]
    O3b = [o for o in orbits if om2 in o][0]
    check("orbit_of_omega1_differs_from_orbit_of_omega2", O3 != O3b)
    check("orbit_3bar_is_minus_orbit_3", sorted(tuple(-x for x in v) for v in O3) == sorted(O3b))
    check("weights_of_3_sum_to_zero", [sum(v[i] for v in O3) for i in range(3)] == [0, 0, 0])
    mins = all(divmod(2 * dot(w, b), dot(b, b))[1] == 0 and abs(exdiv(2 * dot(w, b), dot(b, b))) <= 1
               for w in O3 + O3b for b in R_LONG)
    check("short_roots_are_minuscule_su3_weights", mins)
    diffs = all(tuple(x - y for x, y in zip(u, w)) in set(R_LONG) for u in O3 for w in O3 if u != w)
    check("weight_differences_inside_3_are_su3_roots", diffs)
    chain = [om1, tuple(om1[i] - b1[i] for i in range(3)), tuple(om1[i] - b1[i] - b2[i] for i in range(3))]
    check("3_is_omega1_omega1_minus_b1_omega1_minus_b1_minus_b2", sorted(chain) == sorted(O3), chain=fvl(chain))
    check("omega1_plus_omega2_is_highest_su3_root_adjoint", tuple(x + y for x, y in zip(om1, om2)) ==
          tuple(x + y for x, y in zip(b1, b2)))
    swap = all(sorted(tuple(mv(g, v)) for v in O3) == sorted(O3) for g in Ll)
    check("no_W_A2_element_exchanges_3_and_3bar", swap)
    check("minus_I_exchanges_3_and_3bar", sorted(tuple(mv(mI3, v)) for v in O3) == sorted(O3b))
    row("A2_IN_G2", su3_roots=fvl(R_LONG), weights_3=fvl(O3), weights_3bar=fvl(O3b),
        decomposition="g2 is su3 plus 3 plus 3bar: 14 is 8 plus 3 plus 3", claim="MEASURED")
    row("NAMED", statement="long roots of G2 form the A2 subsystem of su(3) inside g2",
        ref="Borel and de Siebenthal 1949 Comment Math Helv 23", claim="NAMED")
    return O3, O3b, Lg, Sg


# ---------------------------------------------------------------- octonions
LINES = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]
MUL = [[None] * 8 for _ in range(8)]
TABLE_CONFLICT = []


def setm(i, j, s, k):
    if MUL[i][j] is not None and MUL[i][j] != (s, k):
        TABLE_CONFLICT.append((i, j))
    MUL[i][j] = (s, k)


for _i in range(8):
    setm(0, _i, 1, _i)
    setm(_i, 0, 1, _i)
for _i in range(1, 8):
    setm(_i, _i, -1, 0)
for (_a, _b, _c) in LINES:
    for (_x, _y, _z) in ((_a, _b, _c), (_b, _c, _a), (_c, _a, _b)):
        setm(_x, _y, 1, _z)
        setm(_y, _x, -1, _z)

E8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]


def omul(x, y):
    out = [0] * 8
    for i in range(8):
        if x[i] == 0:
            continue
        for j in range(8):
            if y[j] == 0:
                continue
            s, k = MUL[i][j]
            out[k] += s * x[i] * y[j]
    return out


def assoc(x, y, z):
    a = omul(omul(x, y), z)
    b = omul(x, omul(y, z))
    return [p - q for p, q in zip(a, b)]


def is_derivation(D):
    for i in range(8):
        for j in range(8):
            x, y = E8[i], E8[j]
            lhs = mv(D, omul(x, y))
            rhs = [a + b for a, b in zip(omul(mv(D, x), y), omul(x, mv(D, y)))]
            if lhs != rhs:
                return False
    return True


def is_automorphism(P):
    for i in range(8):
        for j in range(8):
            if mv(P, omul(E8[i], E8[j])) != omul(mv(P, E8[i]), mv(P, E8[j])):
                return False
    return True


def R8(p, q):
    M = [[0] * 8 for _ in range(8)]
    M[q][p] = 1
    M[p][q] = -1
    return M


def r7(D):
    return [r[1:] for r in D[1:]]


def R7(p, q):
    return r7(R8(p, q))


def octonion_checks():
    check("octonion_table_complete_no_conflict",
          not TABLE_CONFLICT and all(MUL[i][j] is not None for i in range(8) for j in range(8)))
    pairs = {}
    for ln in LINES:
        for u, w in itertools.combinations(ln, 2):
            k = (min(u, w), max(u, w))
            pairs[k] = pairs.get(k, 0) + 1
    check("fano_each_pair_of_units_on_exactly_one_line", len(pairs) == 21 and set(pairs.values()) == {1})
    alt = True
    nonassoc = 0
    for i, j, k in itertools.product(range(8), repeat=3):
        A = assoc(E8[i], E8[j], E8[k])
        if any(A):
            nonassoc += 1
        B = assoc(E8[j], E8[i], E8[k])
        C = assoc(E8[i], E8[k], E8[j])
        if any(a + b for a, b in zip(A, B)) or any(a + c for a, c in zip(A, C)):
            alt = False
    check("associator_alternating_so_algebra_alternative", alt)
    check("algebra_not_associative", nonassoc > 0, nonzero_basis_associators=nonassoc)
    st = 20261109
    normok = True
    for _ in range(300):
        vec = []
        for _k in range(16):
            st = (st * 1103515245 + 12345) % 2147483648
            vec.append(st % 11 - 5)
        x, y = vec[:8], vec[8:]
        if dot(omul(x, y), omul(x, y)) != dot(x, x) * dot(y, y):
            normok = False
    check("norm_multiplicative_300_integer_pairs", normok)
    row("OCTONIONS", table="e_i e_j is e_k for (i,j,k) cyclic in lines 124 235 346 457 561 672 713; e_i e_i is -1; e0 is 1",
        claim="MEASURED")
    row("NAMED", statement="8-dim alternative composition algebra with positive norm is the octonions",
        ref="Hurwitz 1898; Baez Bull AMS 39 (2002) sec 2", claim="NAMED")


# ---------------------------------------------------------------- PART 4: der(O), Higgs
def derivation_kernel():
    eqs = []
    for i in range(8):
        for j in range(8):
            s0, k0 = MUL[i][j]
            for m in range(8):
                e = [0] * 64
                e[8 * m + k0] += s0
                for a in range(8):
                    s2, p = MUL[a][j]
                    if p == m:
                        e[8 * a + i] -= s2
                for b in range(8):
                    s3, p = MUL[i][b]
                    if p == m:
                        e[8 * b + j] -= s3
                eqs.append(e)
    K = kernel(eqs, 64)
    return [unflat(v, 8) for v in K], len(eqs)


def nice_basis():
    out = []
    labels = []
    for m in range(1, 8):
        prs = sorted((a, b) for a in range(1, 8) for b in range(1, 8) if a != b and MUL[a][b] == (1, m))
        if len(prs) != 3:
            raise RuntimeError("bad pair count")
        out.append(madd(R8(*prs[0]), R8(*prs[1]), -1))
        labels.append("R%d%d-R%d%d" % (prs[0] + prs[1]))
        out.append(madd(R8(*prs[1]), R8(*prs[2]), -1))
        labels.append("R%d%d-R%d%d" % (prs[1] + prs[2]))
    return out, labels


def higgs(T7s, v):
    n = len(T7s)
    cols = [mv(Tm, v) for Tm in T7s]
    M = [[dot(cols[a], cols[b]) for b in range(n)] for a in range(n)]
    K = [[tr(mm(T7s[a], T7s[b])) for b in range(n)] for a in range(n)]
    rk = rank(M)
    h = kernel(M, n)
    Kh = [mv(K, c) for c in h]
    m = kernel(Kh, n)
    Mm = [[dot(x, mv(M, y)) for y in m] for x in m]
    Km = [[dot(x, mv(K, y)) for y in m] for x in m]
    ratio = None
    prop = False
    for i in range(len(m)):
        for j in range(len(m)):
            if Km[i][j] != 0 and ratio is None:
                ratio = reduced_pair(Mm[i][j], -Km[i][j])
    if ratio is not None:
        num, den = ratio
        prop = all(Mm[i][j] * den == -Km[i][j] * num for i in range(len(m)) for j in range(len(m)))
    return {"M": M, "K": K, "rank": rk, "h": h, "m": m, "Mm": Mm, "Km": Km, "ratio": ratio, "prop": prop}


def ipow(b, e):
    # integer power for a non-negative exponent; refuses a negative one (no-float gate)
    assert e >= 0
    r = 1
    for _ in range(e):
        r *= b
    return r


def pencil_ok(M, K, p, q, dh, dm):
    dK = det(K)
    for x in range(-7, 8):
        A = [[q * M[i][j] + x * K[i][j] for j in range(len(M))] for i in range(len(M))]
        if det(A) != ipow(x, dh) * ipow(x - p, dm) * dK:
            return False
    return True


def part4():
    octonion_checks()
    KER, neqs = derivation_kernel()
    check("der_O_linear_system_kernel_dim_14", len(KER) == 14, equations=neqs, unknowns=64, dim=len(KER))
    check("der_O_kernel_vectors_kill_1_and_preserve_ImO",
          all(all(D[0][j] == 0 for j in range(8)) and all(D[i][0] == 0 for i in range(8)) for D in KER))
    check("der_O_kernel_vectors_antisymmetric_so_inside_so7", all(tp(D) == neg(D) for D in KER))
    check("der_O_kernel_vectors_pass_direct_derivation_test", all(is_derivation(D) for D in KER))
    NICE, labels = nice_basis()
    check("nice_basis_14_are_derivations", all(is_derivation(D) for D in NICE), basis=";".join(labels))
    check("nice_basis_rank_14_and_spans_der_O",
          rank([flat(D) for D in NICE]) == 14 and rank([flat(D) for D in NICE + KER]) == 14)
    T7 = [r7(D) for D in NICE]
    T7k = [r7(D) for D in KER]
    brk = [flat(comm(T7[a], T7[b])) for a in range(14) for b in range(a + 1, 14)]
    check("g2_closed_under_bracket", rank([flat(X) for X in T7] + brk) == 14)
    check("g2_perfect", rank(brk) == 14)
    c7 = commutant(T7, 7)
    check("g2_irreducible_on_R7_commutant_dim_1", len(c7) == 1 and unflat(c7[0], 7) in (eye(7), neg(eye(7))), dim=len(c7))
    row("NAMED", statement="der(O) is the compact exceptional Lie algebra g2 of dimension 14",
        ref="E Cartan 1914; Schafer Nonassociative Algebras 1966 ch III; Baez 2002 sec 4.1", claim="NAMED")

    v = [0, 0, 0, 0, 0, 0, 1]
    H = higgs(T7, v)
    Hk = higgs(T7k, v)
    row("HIGGS_GRAM", v="e7", M_diag=fv([H["M"][a][a] for a in range(14)]), rank=H["rank"],
        kernel_dim=len(H["h"]), m_dim=len(H["m"]), claim="MEASURED")
    check("gram_rank_6_is_dim_S6", H["rank"] == 6, rank=H["rank"])
    check("gram_kernel_dim_8_is_stabilizer", len(H["h"]) == 8, dim=len(H["h"]))
    check("complement_m_dim_6", len(H["m"]) == 6)
    check("h_plus_m_is_all_g2", rank(H["h"] + H["m"]) == 14)
    hX = [lincomb(c, T7) for c in H["h"]]
    check("kernel_elements_fix_e7", all(mv(X, v) == [0] * 7 for X in hX))
    hb = [flat(comm(hX[i], hX[j])) for i in range(8) for j in range(i + 1, 8)]
    check("stabilizer_closed_and_perfect_dim8_compact_so_su3", rank([flat(X) for X in hX] + hb) == 8 and rank(hb) == 8)
    num, den = H["ratio"]
    check("M_proportional_to_minus_trace_form_on_m", H["prop"], ratio_num=num, ratio_den=den)
    check("ratio_is_1_over_3", (num, den) == (1, 3), ratio="%d over %d" % (num, den))
    check("M_positive_definite_on_m", all(det([r[:k] for r in H["Mm"][:k]]) > 0 for k in range(1, 7)))
    check("same_rank_kernel_ratio_in_solver_basis", Hk["rank"] == 6 and len(Hk["h"]) == 8 and Hk["prop"]
          and Hk["ratio"] == (1, 3), ratio_solver_basis=fv(Hk["ratio"]))
    check("pencil_det_3M_plus_xK_is_x8_x_minus_1_6_detK_15_points",
          pencil_ok(H["M"], H["K"], 1, 3, 8, 6), spectrum="0 x8 and one third x6 relative to minus trace form")
    # Casimir of g2 on R7 with respect to minus trace form
    B = neg(H["K"])
    dB = det(B)
    adjB = adjugate(B)
    S7 = [[0] * 7 for _ in range(7)]
    for a in range(14):
        for b in range(14):
            if adjB[a][b]:
                S7 = madd(S7, mm(T7[a], T7[b]), adjB[a][b])
    check("casimir_on_R7_is_2", neg(S7) == smul(2 * dB, eye(7)), detB=dB)
    check("mass_sum_rule_6_times_ratio_is_casimir_2", 6 * num == 2 * den)
    row("HIGGS_MASS", unbroken="su3 dim 8 massless", broken="6 vector bosons on m", mass_sq_ratio="1 over 3",
        rule="M equals one third of minus tr on m: -tr(XY) is 3 (Xv).(Yv)",
        units="m_W squared is g squared times norm(v) squared over 3 for -tr-orthonormal generators",
        claim="MEASURED")
    row("NAMED", statement="Higgs vector mass matrix is g^2 (T_a v).(T_b v)",
        ref="Weinberg QFT vol II 1996 sec 21.1", claim="NAMED")
    row("NAMED", statement="G2 acts transitively on S6 with stabilizer SU(3)",
        ref="Montgomery and Samelson Ann Math 44 1943; Borel 1949; Baez 2002 sec 4.1; Bryant Ann Math 126 1987", claim="NAMED")

    # su(3) action on the tangent space T (e1..e6) and the invariant complex structure
    H6 = [[r[:6] for r in X[:6]] for X in hX]
    check("stabilizer_preserves_tangent_space", all(all(X[6][j] == 0 and X[j][6] == 0 for j in range(7)) for X in hX))
    J6 = [[0] * 6 for _ in range(6)]
    jok = True
    for j in range(1, 7):
        prod = omul(E8[7], E8[j])
        if prod[0] != 0 or prod[7] != 0:
            jok = False
        for k in range(1, 7):
            J6[k - 1][j - 1] = prod[k]
    check("J_is_left_mult_by_e7_maps_T_to_T", jok)
    check("J_squared_minus_1_and_orthogonal", mm(J6, J6) == neg(eye(6)) and mm(tp(J6), J6) == eye(6))
    Jsum = [[r[:6] for r in X[:6]] for X in [madd(madd(R7(1, 3), R7(2, 6)), R7(4, 5))]][0]
    check("J_is_R13_plus_R26_plus_R45", J6 == Jsum)
    check("J_commutes_with_su3_on_T", all(mm(J6, X) == mm(X, J6) for X in H6))
    c6 = commutant(H6, 6)
    span_ok = rank(c6 + [flat(eye(6)), flat(J6)]) == 2
    check("su3_commutant_on_T_is_span_I_J_so_T_irreducible_complex_type", len(c6) == 2 and span_ok, dim=len(c6))
    row("COMPLEX_STRUCTURE", J="L_e7 on T_e7 S6", commutant="R[J] isomorphic to C",
        consequence="invariant J with J^2 -1 in R[J] is plus or minus J: unique up to sign", claim="MEASURED")
    row("NAMED", statement="the G2-invariant almost complex structure on S6 from octonions is not integrable",
        ref="Kirchhoff 1947 CR Acad Sci 225; Ehresmann and Libermann 1951 CR Acad Sci 232", claim="NAMED")
    Jbad = [[r[:6] for r in X[:6]] for X in [madd(madd(R7(1, 3), R7(2, 6)), R7(4, 5), -1)]][0]
    control("ctl_wrong_complex_structure_commutes_with_su3",
            mm(Jbad, Jbad) == neg(eye(6)) and all(mm(Jbad, X) == mm(X, Jbad) for X in H6))
    return NICE, T7, H, hX, J6


# ---------------------------------------------------------------- roots read off the algebra
def U_of(al):
    a1, a2 = al
    return (2 * a1 + a2, a2 - a1, -a1 - 2 * a2)


def PHI(v):
    return (v[0] - v[1], v[1] - v[2], v[2] - v[0])


def algebra_roots(T7, H, O3, O3b, Lg, Sg):
    t1, t2 = T7[12], T7[13]
    check("cartan_t1_t2_are_R13_minus_R26_and_R26_minus_R45",
          t1 == madd(R7(1, 3), R7(2, 6), -1) and t2 == madd(R7(2, 6), R7(4, 5), -1))
    check("t1_t2_commute_and_fix_e7", comm(t1, t2) == [[0] * 7 for _ in range(7)]
          and mv(t1, [0] * 6 + [1]) == [0] * 7 and mv(t2, [0] * 6 + [1]) == [0] * 7)

    def ad(t, X):
        return comm(t, X)
    A1 = [ad(t1, X) for X in T7]
    A2 = [ad(t2, X) for X in T7]
    A11 = [ad(t1, X) for X in A1]
    A22 = [ad(t2, X) for X in A2]
    A12 = [ad(t1, X) for X in A2]
    fT = [flat(X) for X in T7]

    def rows_of(mats):
        cols = [flat(X) for X in mats]
        return [[cols[b][e] for b in range(14)] for e in range(49)]
    found = {}
    V0 = kernel(rows_of(A1) + rows_of(A2), 14)
    found[(0, 0)] = V0
    for a1 in range(-3, 4):
        for a2 in range(-3, 4):
            if (a1, a2) == (0, 0) or a1 < 0 or (a1 == 0 and a2 < 0):
                continue
            eq = []
            for (mats, c) in ((A11, a1 * a1), (A22, a2 * a2), (A12, a1 * a2)):
                cols = [flat(X) for X in mats]
                eq += [[cols[b][e] + c * fT[b][e] for b in range(14)] for e in range(49)]
            V = kernel(eq, 14)
            if V:
                found[(a1, a2)] = V
    dims = sum(len(V) for V in found.values())
    check("cartan_centralizer_dim_2_rank_2", len(V0) == 2)
    check("root_space_decomposition_fills_g2_14", dims == 14, total=dims)
    nz = [al for al in found if al != (0, 0)]
    check("six_root_pairs_each_real_dim_2", len(nz) == 6 and all(len(found[al]) == 2 for al in nz),
          pairs=fvl(sorted(nz)))
    Us = set()
    for al in nz:
        Us.add(U_of(al))
        Us.add(U_of((-al[0], -al[1])))
    phiR = {PHI(r) for r in ROOTS}
    check("algebra_roots_equal_PHI_of_part1_roots", Us == phiR,
          PHI="(x-y;y-z;z-x) maps short to short and long to long: times sqrt3 and a 30 degree turn")
    in_h = {}
    for al in nz:
        Xs = [lincomb(c, T7) for c in found[al]]
        in_h[al] = all(mv(X, [0] * 6 + [1]) == [0] * 7 for X in Xs)
    long_al = [al for al in nz if dot(U_of(al), U_of(al)) == 18]
    short_al = [al for al in nz if dot(U_of(al), U_of(al)) == 6]
    check("long_root_spaces_lie_in_stabilizer_su3", len(long_al) == 3 and all(in_h[al] for al in long_al))
    check("short_root_spaces_are_broken_m", len(short_al) == 3 and not any(in_h[al] for al in short_al))
    Kmat = H["K"]
    hsp = H["h"]
    orth = all(dot(c, mv(Kmat, x)) == 0 for al in short_al for c in found[al] for x in hsp)
    check("short_root_spaces_trace_orthogonal_to_su3_so_equal_m", orth)
    # weights of T under (t, J)
    planes = [(1, 3), (2, 6), (4, 5)]
    wts = []
    for k, (p, q) in enumerate(planes):
        wts.append((t1[q - 1][p - 1], t2[q - 1][p - 1]))
    recon = all(t == lincomb([w[i] for w in wts], [R7(*pl) for pl in planes]) for i, t in enumerate((t1, t2)))
    check("t_acts_on_each_J_line_by_weight_times_J", recon, weights=fvl(wts))
    Uw = sorted(U_of(w) for w in wts)
    phi3 = sorted(PHI(v) for v in O3)
    phi3b = sorted(PHI(v) for v in O3b)
    which = "3" if Uw == phi3 else ("3bar" if Uw == phi3b else "neither")
    check("tangent_weights_with_J_form_one_W_A2_orbit", which in ("3", "3bar"), matches=which)
    conj = True
    Sl = [key(g) for g in Sg.values()]
    for l in Lg.values():
        hit = [s for s in Sg.values() if all(PHI(mv(l, b)) == tuple(mv(s, PHI(b))) for b in (B1, B2))]
        conj = conj and len(hit) == 1
    check("PHI_conjugates_long_W_A2_to_plain_permutations", conj and len(Sl) == 6)
    row("ALGEBRA_ROOTS", cartan="t1 R13-R26 and t2 R26-R45", roots_on_t=fvl(sorted(nz)),
        long_in_su3=fvl(sorted(long_al)), short_broken=fvl(sorted(short_al)), tangent_weights=fvl(wts),
        tangent_is=which, claim="MEASURED")
    return t1, t2, planes, which


# ---------------------------------------------------------------- PART 5: e7 -> -e7
def sp_matrix(pi, s):
    P = [[0] * 8 for _ in range(8)]
    P[0][0] = 1
    for i in range(1, 8):
        P[pi[i - 1]][i] = s[i - 1]
    return P


def part5(NICE, J6, t1, t2, hX, which):
    lineset = {frozenset(ln) for ln in LINES}
    colls = []
    for pi in itertools.permutations(range(1, 8)):
        if all(frozenset(pi[x - 1] for x in ln) in lineset for ln in LINES):
            colls.append(pi)
    check("fano_collineations_168", len(colls) == 168, n=len(colls))
    autos = []
    per = {}
    for pi in colls:
        cnt = 0
        for s in itertools.product((1, -1), repeat=7):
            ok = True
            for i in range(1, 8):
                for j in range(1, 8):
                    if i == j:
                        continue
                    s0, k = MUL[i][j]
                    s1, k1 = MUL[pi[i - 1]][pi[j - 1]]
                    if k1 != pi[k - 1] or s[i - 1] * s[j - 1] * s1 != s0 * s[k - 1]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                autos.append((pi, s))
                cnt += 1
        per[cnt] = per.get(cnt, 0) + 1
    check("signed_permutation_automorphisms_1344", len(autos) == 1344, n=len(autos), per_collineation=fvl(sorted(per.items())))
    to_minus = [(pi, s) for pi, s in autos if pi[6] == 7 and s[6] == -1]
    to_plus = [(pi, s) for pi, s in autos if pi[6] == 7 and s[6] == 1]
    check("automorphisms_sending_e7_to_minus_e7_exist_and_form_a_coset", len(to_minus) > 0 and len(to_minus) == len(to_plus),
          to_minus=len(to_minus), to_plus=len(to_plus))
    check("full_64_product_check_on_all_signed_perm_automorphisms",
          all(is_automorphism(sp_matrix(pi, s)) for pi, s in autos), checked=len(autos))
    check("all_signed_perm_automorphisms_det_plus1_in_SO7", all(det(r7(sp_matrix(pi, s))) == 1 for pi, s in autos))
    diag = [(pi, s) for pi, s in to_minus if pi == tuple(range(1, 8))]
    check("diagonal_automorphisms_sending_e7_to_minus_e7", len(diag) == 4, n=len(diag))
    pi, s = diag[0]
    P8 = sp_matrix(pi, s)
    check("phi_full_algebra_automorphism_all_64_products", is_automorphism(P8), signs=fv(s))
    P7 = r7(P8)
    check("phi_sends_e7_to_minus_e7", mv(P7, [0] * 6 + [1]) == [0] * 6 + [-1])
    check("phi_is_involution_det_plus1", mm(P7, P7) == eye(7) and det(P7) == 1)
    fixed = [i for i in range(1, 8) if s[i - 1] == 1]
    sub = all(MUL[a][b][1] in [0] + fixed for a in fixed for b in fixed)
    check("phi_fixes_a_quaternion_subalgebra_pointwise", len(fixed) == 3 and sub, fixed=fv(["e%d" % i for i in fixed]))
    mI7 = neg(eye(7))
    mI8 = [[1 if (i == j == 0) else (-1 if i == j else 0) for j in range(8)] for i in range(8)]
    check("antipodal_minus_1_on_ImO_det_minus_1", det(mI7) == -1)
    control("ctl_antipodal_minus_1_on_ImO_is_automorphism", is_automorphism(mI8))
    check("Ad_phi_preserves_der_O", all(is_derivation(mm(mm(P8, D), tp(P8))) for D in NICE))
    check("Ad_phi_preserves_stabilizer_su3", all(mv(mm(mm(P7, X), tp(P7)), [0] * 6 + [1]) == [0] * 7 for X in hX))
    P6 = [r[:6] for r in P7[:6]]
    check("phi_preserves_tangent_space", all(P7[6][j] == 0 and P7[j][6] == 0 for j in range(6)))
    check("phi_J_phi_inverse_is_minus_J_antiholomorphic_on_T", mm(mm(P6, J6), tp(P6)) == neg(J6))
    check("Ad_phi_is_minus_1_on_cartan_t", mm(mm(P7, t1), tp(P7)) == neg(t1) and mm(mm(P7, t2), tp(P7)) == neg(t2))
    # -I on t: in t coordinates (a,b,c) it is -1, which is not a permutation, so not in W(su3)
    check("minus_1_on_t_not_a_permutation_so_not_in_W_su3",
          neg(eye(3)) not in [[[1 if sig[j] == i else 0 for j in range(3)] for i in range(3)]
                              for sig in itertools.permutations(range(3))])
    row("ANTIPODAL", phi="diag on e1..e7 signs " + fv(s), phi_e7="minus e7",
        weyl_element="minus I in W(G2), outside W(A2)",
        meaning="phi normalizes SU(3), acts as minus 1 on its Cartan, sends J to minus J: 3 to 3bar (tangent " + which + " with J)",
        claim="MEASURED")
    return P8


# ---------------------------------------------------------------- PART 6: broken controls
def part6(NICE):
    bad_roots = [(2, -1, 0) if v == (2, -1, -1) else ((-2, 1, 0) if v == (-2, 1, 1) else v) for v in ROOTS]
    pl, it, cl, _ = rs_checks(bad_roots)
    control("ctl_long_root_moved_off_plane", pl and it and cl, plane=int(pl), integral=int(it), closed=int(cl))
    bad_roots2 = [(3, -2, -1) if v == (2, -1, -1) else ((-3, 2, 1) if v == (-2, 1, 1) else v) for v in ROOTS]
    pl, it, cl, _ = rs_checks(bad_roots2)
    control("ctl_long_root_perturbed_in_plane", pl and it and cl, plane=int(pl), integral=int(it), closed=int(cl))
    badT = [madd(NICE[0], R8(1, 7))] + NICE[1:]
    control("ctl_perturbed_T0_is_derivation", is_derivation(badT[0]))
    H = higgs([r7(D) for D in badT], [0, 0, 0, 0, 0, 0, 1])
    good = H["rank"] == 6 and len(H["h"]) == 8 and H["prop"] and H["ratio"] == (1, 3)
    pen = good and pencil_ok(H["M"], H["K"], 1, 3, 8, 6)
    control("ctl_perturbed_T0_keeps_degenerate_third_masses", good and pen,
            rank=H["rank"], kernel=len(H["h"]), prop=int(H["prop"]),
            ratio=(fv(H["ratio"]) if H["ratio"] else "none"))


# ---------------------------------------------------------------- main
def main():
    part0_gate()
    part1()
    LIFT, W = part2()
    O3, O3b, Lg, Sg = part3(LIFT, W)
    NICE, T7, H, hX, J6 = part4()
    t1, t2, planes, which = algebra_roots(T7, H, O3, O3b, Lg, Sg)
    part5(NICE, J6, t1, t2, hX, which)
    part6(NICE)
    verdict = "ALL_PASS" if not FAILS else "SOME_FAIL"
    row("VERDICT", checks=COUNTS["check"], controls=COUNTS["control"], fails=len(FAILS),
        failed=(fv(FAILS) if FAILS else "none"), verdict=verdict, e=0)
    data = ("\n".join(ROWS) + "\n").encode("utf-8")
    with open(OUT, "wb") as fh:
        fh.write(data)
    hx = hashlib.sha256(data).hexdigest()
    with open(OUT + ".sha256", "wb") as fh:
        fh.write(("%s  %s\n" % (hx, OUT_NAME)).encode("utf-8"))
    print("SHA256|file=%s|sha256=%s|json=0" % (OUT_NAME, hx))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
