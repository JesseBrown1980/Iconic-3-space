# m5_octonion_loop.py -- Opus 5.5 seat, lane m5 of workflow wjiq7jhf8, 2026-11-09.
# Exact integer checks on: the octonions O built from the Fano lines (i,i+1,i+3) mod 7, the eight-square
# identity, the octonion loop L of the 16 signed units (closure, inverses, associativity census, Moufang),
# the associator, the signed-permutation frame group inside G2 (order 1344, kernel 2^3, image GL(3,2)),
# and the seven quaternion copies (seven Q8 inside L).
# Integer only: no float, no true division, no fractions, no decimal, no math module, no exec or eval.
# Run as: python3 -I -B m5_octonion_loop.py   E=0: compute only, no network; writes two files beside itself.
import sys
import os
import io
import hashlib
import tokenize
import itertools

SCRIPT_PATH = os.path.abspath(__file__)
HERE = os.path.dirname(SCRIPT_PATH)
HBP_NAME = "M5-OCTONION-LOOP.hbp"
ROWS = []
FAILS = []


# ------------------------------------------------------------------ rows
def clean(v):
    s = str(v)
    for a, b in (("|", "!"), ("\n", "_"), ("=", ":"), (" ", "_")):
        s = s.replace(a, b)
    return s


def row(tag, *pairs):
    parts = [tag]
    for k, v in pairs:
        parts.append(k + "=" + clean(v))
    parts.append("json=0")
    ROWS.append("|".join(parts))


def check(tag, cond, *pairs):
    ok = bool(cond)
    if not ok:
        FAILS.append(tag)
    row(tag, ("result", "PASS" if ok else "FAIL"), *pairs)
    return ok


def control(tag, failed, *pairs):
    # A broken control passes only when the check it is fed FAILS on it.
    return check(tag, failed, ("observed", "FAIL" if failed else "PASS"), ("expected", "FAIL"), *pairs)


def proof(pid, statement, argument, checked, claim):
    row("PROOF", ("id", pid), ("statement", statement), ("argument", argument), ("checked", checked), ("claim", claim))


# ------------------------------------------------------------------ self lint
def lint_source(src):
    sl = chr(47)
    st = chr(42)
    banned_ops = {sl, sl + "=", st + st, st + st + "="}
    banned_names = {"float", "complex", "Fraction", "Decimal", "fractions", "decimal", "cmath", "math",
                    "eval", "exec", "compile", "__import__", "pickle", "marshal", "json", "numpy",
                    "random", "struct", "globals", "vars"}
    bad = []
    for t in tokenize.tokenize(io.BytesIO(src).readline):
        if t.type == tokenize.OP and t.string in banned_ops:
            bad.append("op@" + str(t.start[0]))
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if "." in s or "e" in s or "j" in s or s.startswith("0x"):
                bad.append("number@" + str(t.start[0]))
        elif t.type == tokenize.NAME and t.string in banned_names:
            bad.append("name:" + t.string + "@" + str(t.start[0]))
    return bad, src.count(bytes([47]))


def section0():
    with open(SCRIPT_PATH, "rb") as fh:
        src = fh.read()
    bad, slashes = lint_source(src)
    check("LINT_INTEGER_ONLY", len(bad) == 0 and slashes == 0, ("bad_tokens", len(bad)), ("slash_chars_in_source", slashes),
          ("banned", "true_division,power_operator,float_or_complex_literals,float_complex_fraction_decimal_math_eval_exec_names"),
          ("claim", "MEASURED"))
    bad_src = ("x = 1 " + chr(47) + " 3\ny = 1" + chr(46) + "5\nz = 2 " + chr(42) + chr(42) + " 3\nw = float(x)\nimport math\n").encode()
    bad2, slashes2 = lint_source(bad_src)
    control("CONTROL_LINT_CATCHES_BAD_SOURCE", len(bad2) == 5 and slashes2 == 1, ("flagged_tokens", len(bad2)),
            ("flagged", ",".join(bad2)), ("sample", "true_division+float_literal+power+float_name+math_import"), ("claim", "MEASURED"))
    vi = sys.version_info
    row("SCRIPT", ("name", "m5_octonion_loop.py"), ("sha256", hashlib.sha256(src).hexdigest()),
        ("python", str(vi[0]) + "." + str(vi[1]) + "." + str(vi[2])), ("isolated_flag", sys.flags.isolated),
        ("E", 0), ("network", 0), ("float_used", 0), ("claim", "MEASURED"))


# ------------------------------------------------------------------ octonion core
LINES = ((1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3))
LINESETS = frozenset(frozenset(l) for l in LINES)


def build_table(lines):
    T = [[None] * 8 for _ in range(8)]
    counts = [0, 0]

    def put(i, j, val):
        if T[i][j] is not None:
            counts[1] += 1
        T[i][j] = val
        counts[0] += 1

    for a in range(8):
        put(0, a, (1, a))
    for a in range(1, 8):
        put(a, 0, (1, a))
        put(a, a, (-1, 0))
    for (a, b, c) in lines:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            put(x, y, (1, z))
            put(y, x, (-1, z))
    holes = sum(1 for i in range(8) for j in range(8) if T[i][j] is None)
    return T, counts[0], counts[1], holes


def omul(T, x, y):
    z = [0] * 8
    for i in range(8):
        xi = x[i]
        if xi:
            Ti = T[i]
            for j in range(8):
                yj = y[j]
                if yj:
                    s, k = Ti[j]
                    z[k] += s * xi * yj
    return z


def vadd(x, y):
    return [x[i] + y[i] for i in range(8)]


def vsub(x, y):
    return [x[i] - y[i] for i in range(8)]


def basis(a):
    v = [0] * 8
    v[a] = 1
    return v


def uvec(u):
    v = [0] * 8
    v[u % 8] = 1 if u < 8 else -1
    return v


def as_unit(v):
    nz = [i for i in range(8) if v[i] != 0]
    if len(nz) != 1:
        return None
    i = nz[0]
    if v[i] == 1:
        return i
    if v[i] == -1:
        return i + 8
    return None


def uname(u):
    return ("-" if u >= 8 else "+") + "e" + str(u % 8)


def vname(v):
    return "[" + ",".join(str(c) for c in v) + "]"


def neg(u):
    return (u + 8) % 16


def loop_table(T):
    M = [[-1] * 16 for _ in range(16)]
    nonunit = 0
    for u in range(16):
        for w in range(16):
            r = as_unit(omul(T, uvec(u), uvec(w)))
            if r is None:
                nonunit += 1
            else:
                M[u][w] = r
    return M, nonunit


# ------------------------------------------------------------------ integer polynomials (dict monomial -> int)
def padd(p, q):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + c
    return {m: c for m, c in r.items() if c != 0}


def pneg(p):
    return {m: -c for m, c in p.items()}


def psub(p, q):
    return padd(p, pneg(q))


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2))
            r[m] = r.get(m, 0) + c1 * c2
    return {m: c for m, c in r.items() if c != 0}


def pvar(i):
    return {(i,): 1}


def pomul(T, X, Y):
    Z = [{} for _ in range(8)]
    for i in range(8):
        if X[i]:
            for j in range(8):
                if Y[j]:
                    s, k = T[i][j]
                    t = pmul(X[i], Y[j])
                    Z[k] = padd(Z[k], t if s == 1 else pneg(t))
    return Z


def pnormsq(X):
    r = {}
    for i in range(8):
        r = padd(r, pmul(X[i], X[i]))
    return r


def pconj(X):
    return [X[0]] + [pneg(X[i]) for i in range(1, 8)]


def oct_sub(X, Y):
    return [psub(X[i], Y[i]) for i in range(8)]


def oct_scale(p, X):
    return [pmul(p, X[i]) for i in range(8)]


def oct_monomials(Z):
    return sum(len(z) for z in Z)


def eight_square_diff(T):
    X = [pvar(i) for i in range(8)]
    Y = [pvar(8 + i) for i in range(8)]
    XY = pomul(T, X, Y)
    lhs = pnormsq(XY)
    rhs = pmul(pnormsq(X), pnormsq(Y))
    return psub(lhs, rhs), lhs, rhs, XY


# ------------------------------------------------------------------ section 1: O, Fano rules, eight-square identity
def section1():
    T, assigned, clashes, holes = build_table(LINES)
    check("OCT_TABLE_WELLFORMED", assigned == 64 and clashes == 0 and holes == 0, ("assigned", assigned), ("clashes", clashes),
          ("holes", holes), ("claim", "MEASURED"))
    cover = {}
    for l in LINES:
        for a in l:
            for b in l:
                if a < b:
                    cover[(a, b)] = cover.get((a, b), 0) + 1
    meets = [len(set(LINES[i]) & set(LINES[j])) for i in range(7) for j in range(i + 1, 7)]
    rule_ok = all(LINES[i] == (i + 1, (i + 1) % 7 + 1, (i + 3) % 7 + 1) for i in range(7))
    check("FANO_PLANE", len(cover) == 21 and all(v == 1 for v in cover.values()) and all(m == 1 for m in meets) and rule_ok,
          ("point_pairs_covered_once", sum(1 for v in cover.values() if v == 1)), ("line_pairs_meeting_in_one_point", sum(1 for m in meets if m == 1)),
          ("lines_are_(i,i+1,i+3)_mod_7", rule_ok), ("claim", "MEASURED"))
    anti = sum(1 for a in range(1, 8) for b in range(1, 8) if a != b and T[b][a] == (-T[a][b][0], T[a][b][1]))
    sq = sum(1 for a in range(1, 8) if T[a][a] == (-1, 0))
    cyc = sum(1 for (a, b, c) in LINES if T[a][b] == (1, c) and T[b][c] == (1, a) and T[c][a] == (1, b))
    check("OCT_CONVENTION", anti == 42 and sq == 7 and cyc == 7,
          ("rule", "e_a.e_b:+e_c_along_each_line_in_cyclic_order;e_b.e_a:-e_a.e_b_for_a!:b;e_a.e_a:-e0;e0_unit"),
          ("anticommuting_ordered_pairs", anti), ("squares_to_minus_one", sq), ("lines_cyclic", cyc),
          ("ref", "Baez_2002_The_Octonions_Bull_AMS_39_145-205_sec_2.1_same_index_rule"), ("claim", "MEASURED_table+NAMED_convention"))
    diff, lhs, rhs, XY = eight_square_diff(T)
    raw = sum(len(c) for c in XY)
    check("OCT_EIGHT_SQUARE_IDENTITY", len(diff) == 0, ("variables", "x0..x7,y0..y7"), ("bilinear_terms_in_xy", raw),
          ("N(xy)_monomials", len(lhs)), ("N(x)N(y)_monomials", len(rhs)), ("difference_monomials", len(diff)),
          ("identity", "N(xy)-N(x)N(y)_is_the_zero_polynomial_in_Z[x0..x7,y0..y7]"),
          ("ref", "Degen_1818_eight-square_identity;Hurwitz_1898"), ("claim", "MEASURED"))
    X = [pvar(i) for i in range(8)]
    Y = [pvar(8 + i) for i in range(8)]
    XpY = [padd(X[i], Y[i]) for i in range(8)]
    dadd = psub(pnormsq(XpY), padd(pnormsq(X), pnormsq(Y)))
    control("CONTROL_NORM_NOT_ADDITIVE", len(dadd) > 0, ("identity_tested", "N(x+y)-N(x)-N(y)"), ("difference_monomials", len(dadd)),
            ("meaning", "the polynomial comparator reports a nonzero difference"), ("claim", "MEASURED"))
    NX = pnormsq(X)
    cX = pconj(X)
    NXe0 = [NX] + [{} for _ in range(7)]
    NXY = oct_scale(NX, Y)
    tests = (
        ("x.xbar-N(x)", oct_sub(pomul(T, X, cX), NXe0)),
        ("xbar.x-N(x)", oct_sub(pomul(T, cX, X), NXe0)),
        ("xbar.(x.y)-N(x)y", oct_sub(pomul(T, cX, pomul(T, X, Y)), NXY)),
        ("x.(xbar.y)-N(x)y", oct_sub(pomul(T, X, pomul(T, cX, Y)), NXY)),
        ("(y.x).xbar-N(x)y", oct_sub(pomul(T, pomul(T, Y, X), cX), NXY)),
        ("(y.xbar).x-N(x)y", oct_sub(pomul(T, pomul(T, Y, cX), X), NXY)),
    )
    for name, Z in tests:
        check("OCT_INVERSE_IDENTITY", oct_monomials(Z) == 0, ("identity", name), ("difference_monomials", oct_monomials(Z)), ("claim", "MEASURED"))
    proof("P1_NORM_MULTIPLICATIVE",
          "N(xy) equals N(x)N(y) for all real octonions x,y; the unit octonions S^7 are closed under multiplication",
          "OCT_EIGHT_SQUARE_IDENTITY shows N(xy)-N(x)N(y) is the zero polynomial in Z[x0..x7,y0..y7]; evaluating the 16 variables at any real numbers is a ring map, so the identity holds for every real x,y; N(x)=N(y)=1 gives N(xy)=1",
          "full symbolic expansion of the 8 bilinear product coefficients (64 bilinear terms) and of both sides; 0 difference monomials",
          "MEASURED+PROOF")
    proof("P1B_TWO_SIDED_INVERSE_PROPERTY",
          "for a in S^7 the maps x to ax and x to xa are bijections of S^7 with inverses x to conj(a)x and x to x conj(a)",
          "the six OCT_INVERSE_IDENTITY rows are zero polynomials, so conj(a)(ax)=N(a)x, a(conj(a)x)=N(a)x, (xa)conj(a)=N(a)x, (x conj(a))a=N(a)x hold for all real a,x; with N(a)=1 they are left and right inverse maps",
          "six degree-2 or degree-3 polynomial identities in 16 integer variables, 0 difference monomials each",
          "MEASURED+PROOF")
    return T


# ------------------------------------------------------------------ section 2: the 16 units L
def assoc_census(M, elems):
    ok = 0
    bad = 0
    first = None
    for x in elems:
        for y in elems:
            Mxy = M[x][y]
            for z in elems:
                if M[Mxy][z] == M[x][M[y][z]]:
                    ok += 1
                else:
                    bad += 1
                    if first is None:
                        first = (x, y, z)
    return ok, bad, first


def section2(T):
    M, nonunit = loop_table(T)
    L = list(range(16))
    check("LOOP_CLOSURE", nonunit == 0, ("ordered_pairs", 256), ("products_not_plus_or_minus_a_unit", nonunit), ("claim", "MEASURED"))
    ident = all(M[0][u] == u and M[u][0] == u for u in L)
    check("LOOP_IDENTITY_E0", ident, ("identity", "+e0"), ("claim", "MEASURED"))
    inv_ok = 0
    for u in L:
        invs = [w for w in L if M[u][w] == 0 and M[w][u] == 0]
        expect = u if u % 8 == 0 else neg(u)
        if invs == [expect]:
            inv_ok += 1
    check("LOOP_TWO_SIDED_INVERSES", inv_ok == 16, ("elements_with_unique_two_sided_inverse_as_expected", inv_ok),
          ("rule", "e_a^-1:-e_a_for_a>0;(-e_a)^-1:e_a;e0^-1:e0;(-e0)^-1:-e0"), ("claim", "MEASURED"))
    rows_ok = all(sorted(M[u]) == L for u in L)
    cols_ok = all(sorted(M[u][w] for u in L) == L for w in L)
    check("LOOP_LATIN_SQUARE", rows_ok and cols_ok, ("rows_are_permutations", rows_ok), ("cols_are_permutations", cols_ok),
          ("meaning", "ax:b_and_ya:b_have_unique_solutions_in_L_so_L_is_a_loop"), ("claim", "MEASURED"))
    ok, bad, first = assoc_census(M, L)
    check("LOOP_ASSOCIATIVITY_CENSUS", ok + bad == 4096 and bad > 0, ("ordered_triples", ok + bad), ("associative", ok), ("nonassociative", bad),
          ("first_failure", ",".join(uname(u) for u in first) if first else "none"), ("claim", "MEASURED"))
    bok = sum(1 for i in range(8) for j in range(8) for k in range(8) if M[M[i][j]][k] == M[i][M[j][k]])
    check("BASIS_ASSOCIATIVITY_CENSUS", bok == 344 and ok == 8 * bok and bad == 8 * (512 - bok), ("basis_triples", 512), ("associative", bok),
          ("nonassociative", 512 - bok), ("sign_lift_factor", 8), ("meaning", "signs are central scalars so each basis triple lifts to 8 unit triples with the same verdict"),
          ("claim", "MEASURED"))
    e1, e2, e3 = basis(1), basis(2), basis(3)
    left = omul(T, omul(T, e1, e2), e3)
    right = omul(T, e1, omul(T, e2, e3))
    check("LOOP_EXPLICIT_FAILING_TRIPLE", left != right and as_unit(left) == 14 and as_unit(right) == 6,
          ("e1.e2", uname(as_unit(omul(T, e1, e2)))), ("e2.e3", uname(as_unit(omul(T, e2, e3)))),
          ("(e1.e2).e3", uname(as_unit(left)) + "_" + vname(left)), ("e1.(e2.e3)", uname(as_unit(right)) + "_" + vname(right)),
          ("associator_[e1,e2,e3]", vname(vsub(left, right))), ("claim", "MEASURED"))
    proof("P2_L_NOT_A_GROUP",
          "L = {+-e0,...,+-e7} is not a group under octonion multiplication",
          "a group operation is associative; (e1e2)e3 = -e6 while e1(e2e3) = +e6, so associativity fails on the triple (e1,e2,e3) of L; L is closed with identity and two-sided inverses, so associativity is the only group axiom that fails",
          "all 256 products in L, identity, 16 inverses, all 4096 ordered triples (2752 associative, 1344 not), explicit triple (e1,e2,e3)",
          "MEASURED+PROOF")
    return M


# ------------------------------------------------------------------ section 3: Moufang, alternativity, subloops
def moufang_census(M, elems):
    f1 = f2 = f3 = 0
    for x in elems:
        for y in elems:
            for z in elems:
                if M[z][M[x][M[z][y]]] != M[M[M[z][x]][z]][y]:
                    f1 += 1
                if M[x][M[z][M[y][z]]] != M[M[M[x][z]][y]][z]:
                    f2 += 1
                if M[M[z][x]][M[y][z]] != M[M[z][M[x][y]]][z]:
                    f3 += 1
    return f1, f2, f3


def polarized_moufang_failures(M):
    # F(x,y,z1,z2) is linear in each slot; check F(x,y,ei,ej)+F(x,y,ej,ei) = 0 on all basis x,y,i,j.
    fails = [0, 0, 0]
    for x in range(8):
        for y in range(8):
            for i in range(8):
                for j in range(8):
                    for (k, f) in enumerate((
                            lambda z1, z2: (M[z1][M[x][M[z2][y]]], M[M[M[z1][x]][z2]][y]),
                            lambda z1, z2: (M[x][M[z1][M[y][z2]]], M[M[M[x][z1]][y]][z2]),
                            lambda z1, z2: (M[M[z1][x]][M[y][z2]], M[M[z1][M[x][y]]][z2]))):
                        a1, b1 = f(i, j)
                        a2, b2 = f(j, i)
                        s = vadd(vsub(uvec(a1), uvec(b1)), vsub(uvec(a2), uvec(b2)))
                        if any(s):
                            fails[k] += 1
    return fails


def closure_set(M, gens):
    S = set(gens)
    frontier = list(S)
    while frontier:
        nxt = []
        cur = list(S)
        for x in frontier:
            for y in cur:
                for p in (M[x][y], M[y][x]):
                    if p not in S:
                        S.add(p)
                        nxt.append(p)
        frontier = nxt
    return S


def is_closed(M, S):
    return all(M[x][y] in S for x in S for y in S)


def section3(T, M):
    L = list(range(16))
    f1, f2, f3 = moufang_census(M, L)
    check("LOOP_MOUFANG_IDENTITIES", f1 == 0 and f2 == 0 and f3 == 0, ("triples", 4096),
          ("z(x(zy)):((zx)z)y_failures", f1), ("x(z(yz)):((xz)y)z_failures", f2), ("(zx)(yz):(z(xy))z_failures", f3), ("claim", "MEASURED"))
    la = sum(1 for x in L for y in L if M[M[x][x]][y] != M[x][M[x][y]])
    ra = sum(1 for x in L for y in L if M[M[y][x]][x] != M[y][M[x][x]])
    fl = sum(1 for x in L for y in L if M[M[x][y]][x] != M[x][M[y][x]])
    check("LOOP_ALTERNATIVE_AND_FLEXIBLE", la == 0 and ra == 0 and fl == 0, ("pairs", 256), ("(xx)y:x(xy)_failures", la),
          ("(yx)x:y(xx)_failures", ra), ("(xy)x:x(yx)_failures", fl), ("claim", "MEASURED"))
    comm = [a for a in L if all(M[a][x] == M[x][a] for x in L)]
    nuc = [a for a in L if all(M[M[a][x]][y] == M[a][M[x][y]] and M[M[x][a]][y] == M[x][M[a][y]] and M[M[x][y]][a] == M[x][M[y][a]] for x in L for y in L)]
    squares = sorted(set(M[x][x] for x in L))
    check("LOOP_CENTER_NUCLEUS_SQUARES", comm == [0, 8] and nuc == [0, 8] and squares == [0, 8],
          ("commutant", ",".join(uname(u) for u in comm)), ("nucleus", ",".join(uname(u) for u in nuc)),
          ("squares", ",".join(uname(u) for u in squares)), ("claim", "MEASURED"))
    wd = all(M[a][b] % 8 == M[neg(a)][b] % 8 == M[a][neg(b)] % 8 for a in L for b in L)
    q = [[M[a][b] % 8 for b in range(8)] for a in range(8)]
    q_comm = all(q[a][b] == q[b][a] for a in range(8) for b in range(8))
    q_assoc = all(q[q[a][b]][c] == q[a][q[b][c]] for a in range(8) for b in range(8) for c in range(8))
    q_exp2 = all(q[a][a] == 0 for a in range(8))
    check("LOOP_QUOTIENT_BY_PLUS_MINUS_ONE", wd and q_comm and q_assoc and q_exp2, ("well_defined", wd), ("commutative", q_comm),
          ("associative", q_assoc), ("every_square_trivial", q_exp2), ("order", 8),
          ("meaning", "L mod {+-e0} is an elementary abelian group of order 8; L is a nonassociative central extension of Z2 by it"),
          ("claim", "MEASURED"))
    census = {}
    group_sizes = {}
    nonassoc_closed = []
    for mask in range(1 << 15):
        S = [0] + [u for u in range(1, 16) if (mask >> (u - 1)) & 1]
        Ss = set(S)
        if not is_closed(M, Ss):
            continue
        census[len(S)] = census.get(len(S), 0) + 1
        o, b, _ = assoc_census(M, S)
        if b == 0:
            group_sizes[len(S)] = group_sizes.get(len(S), 0) + 1
        else:
            nonassoc_closed.append(len(S))
    cen = ",".join(str(k) + "x" + str(census[k]) for k in sorted(census))
    gsz = ",".join(str(k) + "x" + str(group_sizes[k]) for k in sorted(group_sizes))
    check("SUBLOOP_CENSUS", census == {1: 1, 2: 1, 4: 7, 8: 7, 16: 1} and group_sizes == {1: 1, 2: 1, 4: 7, 8: 7} and nonassoc_closed == [16],
          ("subsets_containing_e0", 1 << 15), ("closed_subsets_by_size", cen), ("associative_ones_by_size", gsz),
          ("nonassociative_closed_sizes", ",".join(str(s) for s in nonassoc_closed)),
          ("meaning", "every proper subloop of L is a group; L itself is the only nonassociative subloop"), ("claim", "MEASURED"))
    di_bad = 0
    di_max = 0
    for x in L:
        for y in L:
            S = closure_set(M, [x, y])
            di_max = max(di_max, len(S))
            if assoc_census(M, sorted(S))[1] != 0:
                di_bad += 1
    check("LOOP_DIASSOCIATIVE", di_bad == 0 and di_max == 8, ("pairs", 256), ("pairs_generating_a_nonassociative_subloop", di_bad),
          ("largest_two_generated_subloop", di_max), ("ref", "Artin_theorem_Schafer_1966_Intro_Nonassociative_Algebras_Thm_3.1"),
          ("claim", "MEASURED_on_L+NAMED_theorem"))
    g12 = closure_set(M, [1, 2])
    g123 = closure_set(M, [1, 2, 3])
    gall = closure_set(M, list(range(1, 8)))
    check("LOOP_GENERATION", sorted(g12) == [0, 1, 2, 4, 8, 9, 10, 12] and len(g123) == 16 and len(gall) == 16,
          ("generated_by_e1,e2", ",".join(uname(u) for u in sorted(g12))), ("size_generated_by_e1,e2,e3", len(g123)),
          ("size_generated_by_e1..e7", len(gall)),
          ("meaning", "two imaginary units generate a Q8; three off a common line generate all of L, which is not a group"), ("claim", "MEASURED"))
    pf = polarized_moufang_failures(M)
    check("OCT_MOUFANG_ALL_OF_O_POLARIZED", pf == [0, 0, 0], ("basis_quadruples", 4096), ("identities", 3), ("failures", ",".join(str(f) for f in pf)),
          ("form", "F(x,y,ei,ej)+F(x,y,ej,ei):0_for_F_bilinear_in_the_repeated_variable"), ("claim", "MEASURED"))
    row("NAMED_OCTONION_LOOP", ("object", "L_of_order_16_is_the_octonion_loop_a_Moufang_loop"),
        ("ref", "Moufang_1935_Abh_Math_Sem_Hamburg_11_416-430;Conway-Smith_2003_On_Quaternions_and_Octonions_chapter_on_Moufang_loops(ch.7_by_recollection_task_suggested_ch.9_not_opened_here);Griess_1986_Code_loops_J_Algebra_100_224-234"),
        ("claim", "NAMED"))
    proof("P3_L_MOUFANG_LOOP",
          "L is a Moufang loop of order 16 that is not a group",
          "LOOP_LATIN_SQUARE plus LOOP_IDENTITY_E0 make L a loop; LOOP_MOUFANG_IDENTITIES shows all three Moufang identities hold on every triple; P2 shows it is not associative",
          "4096 triples x 3 identities, 0 failures; 256 products Latin; identity e0",
          "MEASURED+PROOF")
    proof("P4_S7_NOT_A_GROUP",
          "the unit octonions S^7 are not a group under octonion multiplication",
          "by P1 S^7 is closed under multiplication; e1,e2,e3 lie in S^7 (norm 1) and (e1e2)e3 = -e6 differs from e1(e2e3) = +e6, so the operation on S^7 is not associative",
          "P1 polynomial identity; explicit triple in LOOP_EXPLICIT_FAILING_TRIPLE",
          "MEASURED+PROOF")
    proof("P5_S7_MOUFANG_LOOP",
          "S^7 is a Moufang loop under octonion multiplication",
          "each Moufang identity M(x,y,z) is linear in x and y and quadratic in z; write F(z1,z2) for the bilinear form with F(z,z) the identity; then 2F(z,z) equals the sum over i,j of z_i z_j S(i,j) with S(i,j)=F(ei,ej)+F(ej,ei); S vanishes on all basis x,y,i,j, so 2F(z,z)=0 and F(z,z)=0 in R^8 for all x,y,z; with P1 (closure) and P1B (left and right multiplication by a unit are bijections of S^7) S^7 is a loop satisfying the Moufang identities",
          "OCT_MOUFANG_ALL_OF_O_POLARIZED on all 4096 basis quadruples for each of the 3 identities, 0 failures",
          "MEASURED+PROOF")
    return census


# ------------------------------------------------------------------ section 4: associator
def section4(M):
    def A(i, j, k):
        return vsub(uvec(M[M[i][j]][k]), uvec(M[i][M[j][k]]))
    t12 = sum(1 for i in range(8) for j in range(8) for k in range(8) if any(vadd(A(i, j, k), A(j, i, k))))
    t23 = sum(1 for i in range(8) for j in range(8) for k in range(8) if any(vadd(A(i, j, k), A(i, k, j))))
    check("ASSOCIATOR_TOTALLY_ANTISYMMETRIC", t12 == 0 and t23 == 0, ("basis_triples", 512), ("swap12_failures", t12), ("swap23_failures", t23),
          ("meaning", "the transpositions (12),(23) generate S3"), ("claim", "MEASURED"))
    nz = [(i, j, k) for i in range(8) for j in range(8) for k in range(8) if any(A(i, j, k))]
    distinct = [(i, j, k) for i in range(1, 8) for j in range(1, 8) for k in range(1, 8) if len({i, j, k}) == 3]
    assoc_unordered = set(frozenset(t) for t in distinct if not any(A(*t)))
    nonassoc_unordered = set(frozenset(t) for t in distinct if any(A(*t)))
    values = set()
    for t in nz:
        v = A(*t)
        nzc = [c for c in v if c != 0]
        values.add((len(nzc), abs(nzc[0]) if nzc else 0))
    check("ASSOCIATOR_BASIS_CENSUS", len(nz) == 168 and all(len(set(t)) == 3 and 0 not in t for t in nz) and assoc_unordered == set(LINESETS)
          and len(nonassoc_unordered) == 28 and values == {(1, 2)},
          ("nonzero_ordered_basis_triples", len(nz)), ("associative_unordered_imaginary_triples", len(assoc_unordered)),
          ("they_equal_the_fano_lines", assoc_unordered == set(LINESETS)), ("nonassociative_unordered_triples", len(nonassoc_unordered)),
          ("nonzero_associator_shape", "plus_or_minus_2_times_one_unit"), ("example_[e1,e2,e3]", vname(A(1, 2, 3))), ("claim", "MEASURED"))
    proof("P6_O_ALTERNATIVE",
          "the associator [x,y,z] is alternating on all of O, so O is alternative: (xx)y = x(xy) and (yx)x = y(xx) for all x,y",
          "[x,y,z] is trilinear because multiplication is bilinear; antisymmetry under (12) and (23) on basis triples (including repeated indices) extends by trilinearity to all x,y,z; then [x,x,y] = -[x,x,y] gives 2[x,x,y] = 0, so [x,x,y] = 0 in R^8",
          "ASSOCIATOR_TOTALLY_ANTISYMMETRIC on all 512 basis triples, 1024 transposition tests, 0 failures",
          "MEASURED+PROOF")


# ------------------------------------------------------------------ section 5: the frame group in G2
def comp(g, h):
    return (0,) + tuple(g[v] if v > 0 else -g[-v] for v in h[1:])


def pcomp(p, q):
    return tuple(p[q[a]] for a in range(8))


def pi(g):
    return tuple(abs(v) for v in g)


def act_unit(g, u):
    # the signed permutation g applied to the unit with index u (index u+8 is the negative of unit u)
    a = u % 8
    if a == 0:
        return u
    t = g[a] if u < 8 else -g[a]
    return t if t > 0 else 8 - t


def inverse_g(g):
    inv = [0] * 8
    for a in range(1, 8):
        t = g[a]
        inv[abs(t)] = a if t > 0 else -a
    return tuple(inv)


def apply_g(g, v):
    w = [0] * 8
    w[0] = v[0]
    for a in range(1, 8):
        if v[a]:
            t = g[a]
            if t > 0:
                w[t] += v[a]
            else:
                w[-t] -= v[a]
    return w


def table_preserving(T, g):
    for a in range(1, 8):
        for b in range(1, 8):
            sig, c = T[a][b]
            ga, gb = g[a], g[b]
            tau, d = T[abs(ga)][abs(gb)]
            sa = 1 if ga > 0 else -1
            sb = 1 if gb > 0 else -1
            gc = g[c] if c > 0 else 0
            if c == 0:
                if not (d == 0 and sig == sa * sb * tau):
                    return False
            else:
                if d != abs(gc) or sig * (1 if gc > 0 else -1) != sa * sb * tau:
                    return False
    return True


def vector_automorphism(T, g):
    for i in range(8):
        for j in range(8):
            if apply_g(g, omul(T, basis(i), basis(j))) != omul(T, apply_g(g, basis(i)), apply_g(g, basis(j))):
                return False
    return True


SIGNS = [tuple(-1 if (sb >> (a - 1)) & 1 else 1 for a in range(1, 8)) for sb in range(128)]


def enumerate_automorphisms(T):
    found = []
    checked = 0
    for p in itertools.permutations(range(1, 8)):
        for S in SIGNS:
            checked += 1
            g = (0,) + tuple(S[a] * p[a] for a in range(7))
            if table_preserving(T, g):
                found.append(g)
    return found, checked


def order_of(x, mul, ident):
    n = 1
    y = x
    while y != ident:
        y = mul(y, x)
        n += 1
    return n


def gen_group(gens, mul, ident):
    seen = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = mul(x, g)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return seen


def hist_str(h):
    return ",".join(str(k) + ":" + str(h[k]) for k in sorted(h))


def section5(T, M):
    G, checked = enumerate_automorphisms(T)
    check("FRAME_ENUMERATION", checked == 645120 and len(G) == 1344, ("candidates_checked", checked), ("table_preserving_signed_permutations", len(G)),
          ("expect", 1344), ("claim", "MEASURED"))
    vec_ok = sum(1 for g in G if vector_automorphism(T, g))
    check("FRAME_ALL_ARE_ALGEBRA_AUTOMORPHISMS", vec_ok == len(G), ("verified_by_vector_multiplication_on_64_basis_pairs", vec_ok), ("claim", "MEASURED"))
    proof("P7_FRAME_IN_G2",
          "each of the 1344 signed permutations is an automorphism of O and so lies in G2 = Aut(O)",
          "g is linear and bijective (a signed permutation fixing e0); g(ei ej) = g(ei) g(ej) on all 64 basis pairs; both sides of g(xy) = g(x)g(y) are bilinear in (x,y), so they agree on all of O; G2 is by definition (Cartan 1914; Baez 2002 sec 4.1) the automorphism group of O",
          "FRAME_ALL_ARE_ALGEBRA_AUTOMORPHISMS: 1344 x 64 basis-pair checks with full vector multiplication",
          "MEASURED+PROOF")
    Gset = set(G)
    ident = tuple(range(8))
    clos_bad = 0
    hom_bad = 0
    PI = {g: pi(g) for g in G}
    for g in G:
        pg = PI[g]
        for h in G:
            c = comp(g, h)
            if c not in Gset:
                clos_bad += 1
            elif PI[c] != pcomp(pg, PI[h]):
                hom_bad += 1
    inv_bad = 0
    for g in G:
        inv = inverse_g(g)
        if inv not in Gset or comp(g, inv) != ident or comp(inv, g) != ident:
            inv_bad += 1
    check("FRAME_IS_A_GROUP", ident in Gset and clos_bad == 0 and inv_bad == 0 and hom_bad == 0, ("ordered_pairs_composed", len(G) * len(G)),
          ("closure_failures", clos_bad), ("inverse_failures", inv_bad), ("identity_present", ident in Gset),
          ("projection_to_point_permutations_homomorphism_failures", hom_bad), ("claim", "MEASURED"))
    K = [g for g in G if PI[g] == ident]
    k_ab = all(comp(a, b) == comp(b, a) for a in K for b in K)
    k_ex2 = all(comp(a, a) == ident for a in K)
    k_cl = all(comp(a, b) in set(K) for a in K for b in K)
    k_char = all(all((1 if k[a] > 0 else -1) * (1 if k[b] > 0 else -1) == (1 if k[c] > 0 else -1) for (a, b, c) in LINES) for k in K)
    inv_of = {g: inverse_g(g) for g in G}
    k_norm = all(comp(comp(g, k), inv_of[g]) in set(K) for g in G for k in K)
    pats = ",".join("".join("+" if k[a] > 0 else "-" for a in range(1, 8)) for k in sorted(K))
    check("FRAME_KERNEL_SIGN_CHANGES", len(K) == 8 and k_ab and k_ex2 and k_cl and k_char and k_norm, ("order", len(K)), ("abelian", k_ab),
          ("exponent_2", k_ex2), ("closed", k_cl), ("normal", k_norm), ("sign_rule", "s_a.s_b:s_c_on_every_line"), ("patterns_e1..e7", pats),
          ("structure", "elementary_abelian_2^3"), ("claim", "MEASURED"))
    Q = sorted(set(PI[g] for g in G))
    coll = [(0,) + p for p in itertools.permutations(range(1, 8))
            if all(frozenset(((0,) + p)[x] for x in l) in LINESETS for l in LINES)]
    lifts = {}
    for g in G:
        lifts[PI[g]] = lifts.get(PI[g], 0) + 1
    check("FRAME_IMAGE_IS_FANO_COLLINEATION_GROUP", len(Q) == 168 and sorted(coll) == Q and set(lifts.values()) == {8} and len(G) == len(K) * len(Q),
          ("image_order", len(Q)), ("line_preserving_permutations_of_7_points_found_independently", len(coll)),
          ("image_equals_full_collineation_group", sorted(coll) == Q), ("lifts_per_collineation", ",".join(str(v) for v in sorted(set(lifts.values())))),
          ("order_check", "1344:8x168"), ("claim", "MEASURED"))
    line_index = {frozenset(l): i for i, l in enumerate(LINES)}
    line_perms = set()
    for p in Q:
        line_perms.add(tuple(line_index[frozenset(p[x] for x in l)] for l in LINES))
    check("FRAME_ACTION_ON_7_LINES", len(line_perms) == 168, ("distinct_induced_line_permutations", len(line_perms)), ("faithful", len(line_perms) == len(Q)),
          ("claim", "MEASURED"))
    labelings = []
    for p in itertools.permutations(range(1, 8)):
        PHI = (0,) + p
        if all(PHI[a] ^ PHI[b] == PHI[c] for (a, b, c) in LINES):
            labelings.append(PHI)
    gl = set()
    for c1 in range(8):
        for c2 in range(8):
            for c3 in range(8):
                imgs = set()
                for m in range(8):
                    v = 0
                    if m & 1:
                        v ^= c1
                    if m & 2:
                        v ^= c2
                    if m & 4:
                        v ^= c3
                    imgs.add(v)
                if len(imgs) == 8:
                    gl.add((c1, c2, c3))
    PHI = labelings[0]
    PHIinv = [0] * 8
    for a in range(8):
        PHIinv[PHI[a]] = a
    fmap = {}
    add_bad = 0
    for p in Q:
        f = tuple(PHI[p[PHIinv[v]]] for v in range(8))
        if any(f[u ^ v] != f[u] ^ f[v] for u in range(8) for v in range(8)):
            add_bad += 1
        fmap[p] = f
    mats = set((f[1], f[2], f[4]) for f in fmap.values())
    hom_bad2 = 0
    for p in Q:
        for q2 in Q:
            fpq = fmap[pcomp(p, q2)]
            if any(fpq[v] != fmap[p][fmap[q2][v]] for v in range(8)):
                hom_bad2 += 1
    check("FRAME_IMAGE_ISOMORPHIC_TO_GL32", len(labelings) == 168 and len(gl) == 168 and add_bad == 0 and mats <= gl and len(mats) == 168 and hom_bad2 == 0,
          ("fano_labelings_by_F2^3_with_lines_u,v,u+v", len(labelings)), ("GL(3,2)_order_by_enumeration", len(gl)),
          ("collineations_conjugated_to_additive_maps", len(Q) - add_bad), ("distinct_matrices", len(mats)), ("homomorphism_failures", hom_bad2),
          ("meaning", "conjugation by the labeling is an injective homomorphism into GL(3,2) and both have order 168, so it is an isomorphism"),
          ("claim", "MEASURED"))
    sl = [(a, b, c, d) for a in range(7) for b in range(7) for c in range(7) for d in range(7) if (a * d - b * c) % 7 == 1]

    def mm7(X, Y):
        return ((X[0] * Y[0] + X[1] * Y[2]) % 7, (X[0] * Y[1] + X[1] * Y[3]) % 7, (X[2] * Y[0] + X[3] * Y[2]) % 7, (X[2] * Y[1] + X[3] * Y[3]) % 7)
    h7 = {}
    for X in sl:
        n = 1
        Y = X
        while Y != (1, 0, 0, 1) and Y != (6, 0, 0, 6):
            Y = mm7(Y, X)
            n += 1
        h7[n] = h7.get(n, 0) + 1
    even = all(v % 2 == 0 for v in h7.values())
    hp = {k: divmod(v, 2)[0] for k, v in h7.items()}
    hq = {}
    for p in Q:
        o = order_of(p, pcomp, ident)
        hq[o] = hq.get(o, 0) + 1
    check("FINGERPRINT_PSL27", len(sl) == 336 and even and hp == hq, ("SL(2,7)_order", len(sl)), ("PSL(2,7)_element_orders", hist_str(hp)),
          ("collineation_group_element_orders", hist_str(hq)), ("match", hp == hq),
          ("scope", "equal_order_statistics_only;the_isomorphism_GL(3,2)_iso_PSL(2,7)_is_NAMED_ATLAS_of_Finite_Groups_1985_L2(7)_iso_L3(2)"), ("claim", "MEASURED_fingerprint+NAMED_isomorphism"))
    hg = {}
    for g in G:
        o = order_of(g, comp, ident)
        hg[o] = hg.get(o, 0) + 1
    row("FRAME_ELEMENT_ORDERS", ("order_histogram", hist_str(hg)), ("total", sum(hg.values())), ("claim", "MEASURED"))
    # complement search: any complement C of K maps isomorphically onto Q, so it holds exactly one lift of each of a and b,
    # and those two lifts generate C. Trying all 8 x 8 lift pairs decides whether a complement exists.
    inv7 = [p for p in Q if order_of(p, pcomp, ident) == 2]
    ord3 = [p for p in Q if order_of(p, pcomp, ident) == 3]
    gen_pair = None
    for a in inv7:
        for b in ord3:
            if len(gen_group([a, b], pcomp, ident)) == 168:
                gen_pair = (a, b)
                break
        if gen_pair:
            break
    by_pi = {}
    for g in G:
        by_pi.setdefault(PI[g], []).append(g)
    sizes = {}
    for la in by_pi[gen_pair[0]]:
        for lb in by_pi[gen_pair[1]]:
            n = len(gen_group([la, lb], comp, ident))
            sizes[n] = sizes.get(n, 0) + 1
    split = 168 in sizes
    check("FRAME_EXTENSION_NONSPLIT", gen_pair is not None and not split, ("generating_pair_orders", "2,3"), ("lift_pairs_tried", sum(sizes.values())),
          ("generated_subgroup_orders", hist_str(sizes)), ("complement_of_order_168_found", split),
          ("meaning", "no subgroup of order 168 meets the kernel trivially, so the extension 2^3.GL(3,2) does not split"), ("claim", "MEASURED"))
    ord7 = [p for p in Q if order_of(p, pcomp, ident) == 7]
    f21 = None
    for s in ord7:
        for t in ord3:
            if len(gen_group([s, t], pcomp, ident)) == 21:
                f21 = (s, t)
                break
        if f21:
            break
    sizes21 = {}
    for la in by_pi[f21[0]]:
        for lb in by_pi[f21[1]]:
            n = len(gen_group([la, lb], comp, ident))
            sizes21[n] = sizes21.get(n, 0) + 1
    check("POSITIVE_CONTROL_COMPLEMENT_FINDER_F21", 21 in sizes21, ("subgroup", "Frobenius_7:3_of_order_21_in_the_image"),
          ("lift_pairs_tried", sum(sizes21.values())), ("generated_subgroup_orders", hist_str(sizes21)),
          ("expected", "a_complement_of_order_21_exists_by_Schur-Zassenhaus_since_gcd(8,21):1"), ("claim", "MEASURED"))
    proof("P8_FRAME_STRUCTURE",
          "the stabilizer of the frame {+-e1,...,+-e7} in G2 is a nonsplit extension 2^3.GL(3,2) of order 1344, acting on the 7 points and the 7 lines of the Fano plane through its full collineation group",
          "every automorphism preserving the frame is a signed permutation fixing e0 (it fixes 1 and maps the unit set {+-e1..+-e7} to itself), so the enumerated 1344 are the whole stabilizer; the projection to point permutations is a homomorphism with kernel the 8 sign changes (first isomorphism theorem: 1344 = 8 x 168); its image is all 168 collineations, isomorphic to GL(3,2) by the F2^3 labeling; the 64-lift search finds no complement",
          "FRAME_ENUMERATION, FRAME_IS_A_GROUP, FRAME_KERNEL_SIGN_CHANGES, FRAME_IMAGE_IS_FANO_COLLINEATION_GROUP, FRAME_IMAGE_ISOMORPHIC_TO_GL32, FRAME_EXTENSION_NONSPLIT",
          "MEASURED+PROOF")
    proof("P8B_AUT_L",
          "the automorphism group of the loop L equals this 1344-element group",
          "a loop automorphism f fixes e0 and maps squares to squares; the squares are {e0,-e0} (LOOP_CENTER_NUCLEUS_SQUARES), so f(-e0) = -e0 and f(-x) = f((-e0)x) = -f(x); f therefore restricts to a signed permutation of e1..e7 preserving the table, which is one of the 1344; each of the 1344 maps L to L and preserves products",
          "LOOP_CENTER_NUCLEUS_SQUARES and FRAME_ENUMERATION",
          "MEASURED+PROOF")
    nonassoc = set((x, y, z) for x in range(16) for y in range(16) for z in range(16) if M[M[x][y]][z] != M[x][M[y][z]])
    orbit = [(act_unit(g, 1), act_unit(g, 2), act_unit(g, 3)) for g in G]
    check("FRAME_SIMPLY_TRANSITIVE_ON_NONASSOCIATIVE_TRIPLES", len(set(orbit)) == len(G) and set(orbit) == nonassoc and len(nonassoc) == 1344,
          ("nonassociative_ordered_unit_triples", len(nonassoc)), ("orbit_of_(e1,e2,e3)", len(set(orbit))), ("stabilizer_order", divmod(len(G), len(set(orbit)))[0]),
          ("orbit_equals_the_nonassociative_triples", set(orbit) == nonassoc),
          ("meaning", "1344 is both the frame group order and the nonassociative triple count because the group acts freely and transitively on those triples"),
          ("ref", "Baez_2002_sec_4.1_G2_acts_simply_transitively_on_basic_triples"), ("claim", "MEASURED_finite_case+NAMED_continuous_case"))
    row("NAMED_FRAME_GROUP_IN_G2", ("statement", "2^3.GL(3,2)_of_order_1344_is_a_maximal_finite_subgroup_of_G2(C)_and_of_compact_G2"),
        ("ref", "Cohen-Wales_1983_Finite_subgroups_of_G2(C)_Comm_Algebra_11_441-459;Griess_1995_Basic_conjugacy_theorems_for_G2_Invent_Math_121_257-277;Wilson_2009_The_Finite_Simple_Groups_sec_4.3"),
        ("checked_here", "order_structure_nonsplit_and_membership_in_Aut(O);maximality_not_checked"), ("claim", "NAMED"))
    return G, Gset


# ------------------------------------------------------------------ section 6: the seven quaternion copies
def gadd(z, w):
    return (padd(z[0], w[0]), padd(z[1], w[1]))


def gsub(z, w):
    return (psub(z[0], w[0]), psub(z[1], w[1]))


def gmul(z, w):
    return (psub(pmul(z[0], w[0]), pmul(z[1], w[1])), padd(pmul(z[0], w[1]), pmul(z[1], w[0])))


def rho(X, a, b, c, s21):
    # q = q0 + qa i + qb j + qc k  ->  [[q0 + qa I, qb + qc I], [-qb + qc I, q0 - qa I]]  (I the Gaussian unit)
    bot = X[b] if s21 == 1 else pneg(X[b])
    return [[(X[0], X[a]), (X[b], X[c])], [(bot, X[c]), (X[0], pneg(X[a]))]]


def m2mul(A, B):
    return [[gadd(gmul(A[i][0], B[0][j]), gmul(A[i][1], B[1][j])) for j in range(2)] for i in range(2)]


def m2zero(A, B):
    return sum(len(gsub(A[i][j], B[i][j])[0]) + len(gsub(A[i][j], B[i][j])[1]) for i in range(2) for j in range(2))


def rho_num(v, a, b, c):
    return [[(v[0], v[a]), (v[b], v[c])], [(-v[b], v[c]), (v[0], -v[a])]]


def gi_mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def m2num_mul(A, B):
    out = []
    for i in range(2):
        r = []
        for j in range(2):
            p = gi_mul(A[i][0], B[0][j])
            q = gi_mul(A[i][1], B[1][j])
            r.append((p[0] + q[0], p[1] + q[1]))
        out.append(r)
    return out


def m2num_dagger(A):
    return [[(A[j][i][0], -A[j][i][1]) for j in range(2)] for i in range(2)]


def section6(T, M, census):
    copies = []
    for (a, b, c) in LINES:
        S = [0, 8, a, a + 8, b, b + 8, c, c + 8]
        Ss = set(S)
        closed = is_closed(M, Ss)
        _, nb, _ = assoc_census(M, S)
        sub_basis = [0, a, b, c]
        span_closed = all(M[x][y] % 8 in sub_basis for x in sub_basis for y in sub_basis)
        span_assoc = all(M[M[x][y]][z] == M[x][M[y][z]] for x in sub_basis for y in sub_basis for z in sub_basis)
        hamilton = (M[a][a] == 8 and M[b][b] == 8 and M[c][c] == 8 and M[a][b] == c and M[b][c] == a and M[c][a] == b
                    and M[b][a] == c + 8 and M[c][b] == a + 8 and M[a][c] == b + 8 and M[M[a][b]][c] == 8)
        orders = {}
        for u in S:
            o = order_of(u, lambda x, y: M[x][y], 0)
            orders[o] = orders.get(o, 0) + 1
        nonab = any(M[x][y] != M[y][x] for x in S for y in S)
        gen2 = closure_set(M, [a, b]) == Ss
        gen3 = closure_set(M, [a, b, c]) == Ss
        P = [{} for _ in range(8)]
        Qp = [{} for _ in range(8)]
        for idx, slot in enumerate((0, a, b, c)):
            P[slot] = pvar(idx)
            Qp[slot] = pvar(4 + idx)
        PQ = pomul(T, P, Qp)
        outside = sum(len(PQ[k]) for k in range(8) if k not in (0, a, b, c))
        rep_bad = m2zero(m2mul(rho(P, a, b, c, -1), rho(Qp, a, b, c, -1)), rho(PQ, a, b, c, -1))
        R = rho(P, a, b, c, -1)
        det = gsub(gmul(R[0][0], R[1][1]), gmul(R[0][1], R[1][0]))
        nP = {}
        for slot in (0, a, b, c):
            nP = padd(nP, pmul(P[slot], P[slot]))
        det_bad = len(psub(det[0], nP)) + len(det[1])
        Rc = rho(pconj(P), a, b, c, -1)
        Rd = [[(R[j][i][0], pneg(R[j][i][1])) for j in range(2)] for i in range(2)]
        dag_bad = m2zero(Rc, Rd)
        su2_bad = 0
        for u in S:
            U = rho_num(uvec(u), a, b, c)
            UU = m2num_mul(U, m2num_dagger(U))
            d = gi_mul(U[0][0], U[1][1])
            e = gi_mul(U[0][1], U[1][0])
            if UU != [[(1, 0), (0, 0)], [(0, 0), (1, 0)]] or (d[0] - e[0], d[1] - e[1]) != (1, 0):
                su2_bad += 1
        hom_num_bad = sum(1 for x in S for y in S if m2num_mul(rho_num(uvec(x), a, b, c), rho_num(uvec(y), a, b, c)) != rho_num(uvec(M[x][y]), a, b, c))
        bad_rep = m2zero(m2mul(rho(P, a, b, c, 1), rho(Qp, a, b, c, 1)), rho(PQ, a, b, c, 1))
        ok = (closed and nb == 0 and span_closed and span_assoc and hamilton and orders == {1: 1, 2: 1, 4: 6} and nonab and gen2 and gen3
              and outside == 0 and rep_bad == 0 and det_bad == 0 and dag_bad == 0 and su2_bad == 0 and hom_num_bad == 0)
        check("QUATERNION_COPY", ok, ("line", "(" + ",".join(str(t) for t in (a, b, c)) + ")"), ("span_1,ea,eb,ec_closed", span_closed),
              ("span_associative_on_64_basis_triples", span_assoc), ("hamilton_table_i:ea,j:eb,k:ec", hamilton),
              ("units_closed", closed), ("unit_triples_nonassociative", nb), ("element_orders", hist_str(orders)), ("nonabelian", nonab),
              ("generated_by_ea,eb", gen2), ("symbolic_product_outside_span_monomials", outside),
              ("rho(p)rho(q)-rho(pq)_monomials", rep_bad), ("det_rho-N_monomials", det_bad), ("rho(conj_q)-rho(q)^dagger_monomials", dag_bad),
              ("Q8_units_not_in_SU(2)_exactly", su2_bad), ("rho_on_Q8_homomorphism_failures", hom_num_bad), ("claim", "MEASURED"))
        control("CONTROL_WRONG_SU2_SIGN", bad_rep > 0, ("line", "(" + ",".join(str(t) for t in (a, b, c)) + ")"),
                ("change", "bottom_left_entry_+qb_instead_of_-qb"), ("difference_monomials", bad_rep), ("claim", "MEASURED"))
        copies.append(frozenset(S))
    order8 = []
    for mask in range(1 << 15):
        if bin(mask).count("1") != 7:
            continue
        S = frozenset([0] + [u for u in range(1, 16) if (mask >> (u - 1)) & 1])
        if is_closed(M, S):
            order8.append(S)
    inter_ok = 0
    for i in range(7):
        for j in range(i + 1, 7):
            common = set(LINES[i]) & set(LINES[j])
            x = min(common)
            if len(common) == 1 and copies[i] & copies[j] == frozenset({0, 8, x, x + 8}):
                inter_ok += 1
    cover = [sum(1 for C in copies if u in C) for u in range(1, 8)]
    union = frozenset().union(*copies)
    check("SEVEN_Q8_IN_L", len(set(copies)) == 7 and set(order8) == set(copies) and inter_ok == 21 and cover == [3] * 7 and union == frozenset(range(16))
          and census.get(8) == 7,
          ("copies", len(set(copies))), ("order_8_subloops_found_by_search", len(order8)), ("they_equal_the_line_copies", set(order8) == set(copies)),
          ("pairs_meeting_in_+-e0,+-e_x_at_the_common_point", inter_ok), ("copies_through_each_imaginary_unit", ",".join(str(c) for c in cover)),
          ("union_is_L", union == frozenset(range(16))), ("claim", "MEASURED"))
    nonline = {0, 8, 1, 9, 2, 10, 3, 11}
    control("CONTROL_NONLINE_SUBSET_NOT_CLOSED", not is_closed(M, nonline), ("subset", "+-e0,+-e1,+-e2,+-e3"), ("witness", "e1.e2:" + uname(M[1][2])),
            ("claim", "MEASURED"))
    proof("P9_QUATERNION_COPIES",
          "for each Fano line (a,b,c), span(1,ea,eb,ec) is a subalgebra isomorphic to H; its unit sphere is a copy of SU(2) inside S^7, and {+-1,+-ea,+-eb,+-ec} is a Q8 subgroup of L",
          "the span is closed (products of its basis stay in it) and its multiplication table is Hamilton's with i,j,k = ea,eb,ec; associativity holds on its 64 basis triples and so, by trilinearity of the associator, on the whole span; rho is linear, injective, and multiplicative as a polynomial identity, with det rho(q) = N(q) and rho(conj q) = rho(q)^dagger, so for N(q)=1 it lands in SU(2) and is a group isomorphism onto its image; that the image is all of SU(2) is NAMED (every SU(2) matrix has the form [[u,v],[-conj v, conj u]] with u conj(u) + v conj(v) equal to 1)",
          "QUATERNION_COPY rows for all 7 lines; CONTROL_WRONG_SU2_SIGN shows the multiplicativity test fails on a sign error",
          "MEASURED+PROOF+NAMED_surjectivity")
    proof("P10_SEVEN_Q8_GLUED",
          "L is the union of exactly seven Q8 subgroups, one per Fano line; any two meet in {+-1,+-e_x} for the common point x; each imaginary unit lies in exactly three",
          "exhaustive search over all 6435 subsets of size 8 containing e0 finds exactly 7 closed ones, equal to the line copies; the 21 pairwise intersections are computed directly; any two Fano lines meet in exactly one point (FANO_PLANE)",
          "SEVEN_Q8_IN_L and SUBLOOP_CENSUS",
          "MEASURED+PROOF")


# ------------------------------------------------------------------ section 6b: the quaternion-side parallel
def section6b(M, G):
    all_ok = True
    for (a, b, c) in LINES:
        pts = (a, b, c)
        auts = []
        for p in itertools.permutations(pts):
            for sb in range(8):
                img = {0: 0, 8: 8}
                for idx, x in enumerate(pts):
                    s = -1 if (sb >> idx) & 1 else 1
                    img[x] = p[idx] if s == 1 else p[idx] + 8
                    img[x + 8] = neg(img[x])
                S = [0, 8, a, a + 8, b, b + 8, c, c + 8]
                if all(img[M[x][y]] == M[img[x]][img[y]] for x in S for y in S):
                    auts.append(tuple(img[x] for x in S))
        S = [0, 8, a, a + 8, b, b + 8, c, c + 8]
        pos = {x: i for i, x in enumerate(S)}

        def cmp_s(f, h):
            return tuple(f[pos[h[i]]] for i in range(8))
        ident8 = tuple(S)
        kern = [f for f in auts if all(f[pos[x]] % 8 == x for x in (a, b, c))]
        imgperm = set(tuple(f[pos[x]] % 8 for x in (a, b, c)) for f in auts)
        hist = {}
        for f in auts:
            o = order_of(f, cmp_s, ident8)
            hist[o] = hist.get(o, 0) + 1
        stab = [g for g in G if {abs(g[a]), abs(g[b]), abs(g[c])} == {a, b, c}]
        restr = set(tuple(act_unit(g, x) for x in S) for g in stab)
        ok = (len(auts) == 24 and len(kern) == 4 and len(imgperm) == 6 and hist == {1: 1, 2: 9, 3: 8, 4: 6}
              and len(stab) == 192 and restr == set(auts))
        all_ok = all_ok and ok
        check("QUATERNION_FRAME_AUTOMORPHISMS", ok, ("line", "(" + ",".join(str(t) for t in pts) + ")"), ("signed_permutation_candidates", 48),
              ("table_preserving", len(auts)), ("sign_kernel", len(kern)), ("axis_permutations", len(imgperm)), ("element_orders", hist_str(hist)),
              ("octonion_frame_elements_stabilizing_the_line", len(stab)), ("their_restrictions_equal_these_24", restr == set(auts)),
              ("claim", "MEASURED"))
    row("NAMED_QUATERNION_SIDE", ("statement", "Aut(Q8)_iso_S4_is_the_frame_stabilizer_in_Aut(H)_iso_SO(3);_unit_quaternions_S^3_iso_SU(2)"),
        ("ref", "Skolem-Noether_(automorphisms_of_H_are_inner);Conway-Smith_2003_chapters_on_quaternions_and_3D_groups"),
        ("measured_here", "order_24_with_S4_order_statistics_and_kernel_4_image_6;_rho_into_SU(2)_on_each_copy"), ("claim", "NAMED"))
    proof("P11_PARALLEL",
          "quaternion side: Q8 sits in the group S^3 iso SU(2), and the table symmetries of its frame form a group of order 24 (S4 order statistics) inside Aut(H) iso SO(3); octonion side: L sits in S^7, which is only a Moufang loop, and the table symmetries of its frame form the group 2^3.GL(3,2) of order 1344 inside Aut(O) = G2; every one of the 24 extends to exactly 8 elements of the 1344",
          "the 24 are found by exhaustive search over 48 signed permutations of each Fano line; the 192 elements of the 1344-group that stabilize that line restrict onto exactly these 24, so restriction is a surjective homomorphism with kernel of order 192 divided by 24, that is 8",
          "QUATERNION_FRAME_AUTOMORPHISMS for all 7 lines; P4, P5, P7, P8",
          "MEASURED+PROOF")
    return all_ok


# ------------------------------------------------------------------ section 7: broken controls
def section7(T, G, Gset):
    broken_lines = ((2, 1, 4),) + LINES[1:]
    TB, assigned, clashes, holes = build_table(broken_lines)
    check("CONTROL_BROKEN_TABLE_WELLFORMED", assigned == 64 and clashes == 0 and holes == 0, ("change", "line_(1,2,4)_reversed_to_(2,1,4)"),
          ("claim", "MEASURED"))
    diff, _, _, _ = eight_square_diff(TB)
    control("CONTROL_BROKEN_LINE_NORM", len(diff) > 0, ("change", "line_(1,2,4)_reversed"), ("difference_monomials", len(diff)), ("claim", "MEASURED"))
    MB, nonunit = loop_table(TB)
    f1, f2, f3 = moufang_census(MB, list(range(16)))
    control("CONTROL_BROKEN_LINE_MOUFANG_ON_L", nonunit == 0 and (f1 + f2 + f3) > 0, ("change", "line_(1,2,4)_reversed"),
            ("failures_per_identity", ",".join(str(f) for f in (f1, f2, f3))), ("claim", "MEASURED"))
    pf = polarized_moufang_failures(MB)
    control("CONTROL_BROKEN_LINE_MOUFANG_ON_O", sum(pf) > 0, ("failures_per_identity", ",".join(str(f) for f in pf)), ("claim", "MEASURED"))
    t12 = 0
    for i in range(8):
        for j in range(8):
            for k in range(8):
                a = vsub(uvec(MB[MB[i][j]][k]), uvec(MB[i][MB[j][k]]))
                b = vsub(uvec(MB[MB[j][i]][k]), uvec(MB[j][MB[i][k]]))
                if any(vadd(a, b)):
                    t12 += 1
    control("CONTROL_BROKEN_LINE_ASSOCIATOR_ANTISYMMETRY", t12 > 0, ("swap12_failures", t12), ("claim", "MEASURED"))
    GB, checked = enumerate_automorphisms(TB)
    control("CONTROL_BROKEN_LINE_FRAME_COUNT", len(GB) != 1344, ("candidates_checked", checked), ("table_preserving_count", len(GB)),
            ("meaning", "the enumerator responds to the table; it does not return 1344 regardless"), ("claim", "MEASURED"))
    g0 = None
    for g in G:
        if g != tuple(range(8)) and g[1] != 1:
            g0 = g
            break
    gs = list(g0)
    gs[1] = -gs[1]
    gs = tuple(gs)
    control("CONTROL_PERTURBED_AUTOMORPHISM_SIGN", not table_preserving(T, gs) and not vector_automorphism(T, gs) and gs not in Gset,
            ("base", ",".join(str(v) for v in g0[1:])), ("perturbed", ",".join(str(v) for v in gs[1:])), ("change", "image_of_e1_negated"), ("claim", "MEASURED"))
    gw = list(g0)
    gw[1], gw[2] = gw[2], gw[1]
    gw = tuple(gw)
    control("CONTROL_PERTURBED_AUTOMORPHISM_SWAP", not table_preserving(T, gw) and not vector_automorphism(T, gw) and gw not in Gset,
            ("perturbed", ",".join(str(v) for v in gw[1:])), ("change", "images_of_e1_and_e2_swapped"), ("claim", "MEASURED"))
    conj = (0, -1, -2, -3, -4, -5, -6, -7)
    control("CONTROL_CONJUGATION_NOT_AUTOMORPHISM", not table_preserving(T, conj) and not vector_automorphism(T, conj) and conj not in Gset,
            ("map", "e_a_to_-e_a_for_all_a>0"), ("meaning", "conjugation reverses products, an anti-automorphism"), ("claim", "MEASURED"))


def main():
    section0()
    T = section1()
    M = section2(T)
    census = section3(T, M)
    section4(M)
    G, Gset = section5(T, M)
    section6(T, M, census)
    section6b(M, G)
    section7(T, G, Gset)
    n_pass = sum(1 for r in ROWS if "|result=PASS|" in r)
    n_fail = sum(1 for r in ROWS if "|result=FAIL|" in r)
    row("SUMMARY", ("checks_pass", n_pass), ("checks_fail", n_fail), ("failed_tags", ",".join(FAILS) if FAILS else "none"),
        ("verdict", "ALL_PASS" if not FAILS else "SOME_FAIL"), ("E", 0), ("claim", "MEASURED"))
    data = ("\n".join(ROWS) + "\n").encode()
    out = os.path.join(HERE, HBP_NAME)
    with open(out, "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(out + ".sha256", "wb") as fh:
        fh.write((digest + "  " + HBP_NAME + "\n").encode())
    for r in ROWS:
        print(r)
    print("HBP_SHA256=" + digest)
    if FAILS:
        for t in FAILS:
            print("FAIL|tag=" + t + "|json=0")
        sys.exit(1)


if __name__ == "__main__":
    main()
