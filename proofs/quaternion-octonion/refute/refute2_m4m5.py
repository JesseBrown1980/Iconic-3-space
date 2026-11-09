# refute2_m4m5.py -- second-seat refuter for lanes m4 (quaternions, SU(2), binary groups) and m5 (octonion loop, frame group)
# of workflow wjiq7jhf8, 2026-11-09. Independent code: nothing is imported from the m4 or m5 scripts.
# Integer only: no float, no true division, no power operator, no fractions, no decimal, no math module, no exec or eval.
# Run as: python3 -I -B refute2_m4m5.py    E=0: compute only, no network; writes REFUTE2.hbp and REFUTE2.hbp.sha256 beside itself.
import sys
import os
import io
import hashlib
import tokenize
import itertools

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
OUT_NAME = "REFUTE2.hbp"
ROWS = []
FAILS = []
COUNT = {"checks": 0, "controls": 0}
SL = chr(47)
ST = chr(42)
PC = chr(37)


# ------------------------------------------------------------------ rows
def clean(v):
    s = str(v)
    for a, b in (("|", "!"), ("\n", "_"), ("\r", "_"), ("=", ":"), (" ", "_")):
        s = s.replace(a, b)
    return s


def row(tag, *pairs):
    parts = [tag]
    for k, v in pairs:
        parts.append(k + "=" + clean(v))
    parts.append("json=0")
    ROWS.append("|".join(parts))


def check(tag, cond, *pairs):
    COUNT["checks"] += 1
    ok = bool(cond)
    if not ok:
        FAILS.append(tag)
    row(tag, ("result", "PASS" if ok else "FAIL"), *(pairs + (("claim", "MEASURED"),)))
    return ok


def control(tag, observed_fail, *pairs):
    # a deliberately broken input; the control passes only when the check fed with it FAILS
    COUNT["controls"] += 1
    return check(tag, observed_fail, ("observed", "FAIL" if observed_fail else "PASS"), ("expected", "FAIL"), *pairs)


def proof(pid, claim, statement, argument, checked):
    row("PROOF", ("id", pid), ("claim", claim), ("statement", statement), ("argument", argument), ("checked", checked))


def hist_str(h):
    return ",".join(str(k) + ":" + str(h[k]) for k in sorted(h))


# ------------------------------------------------------------------ A. static integer-only audit (tokenize, independent of the scripts' own lints)
BANNED_NAMES = {"float", "complex", "Fraction", "Decimal", "fractions", "decimal", "cmath", "eval", "exec", "compile",
                "__import__", "pickle", "marshal", "globals", "vars", "locals", "getattr", "setattr", "delattr", "pow",
                "numpy", "statistics", "random", "time", "struct", "array", "ctypes", "runpy", "importlib", "__builtins__",
                "shelve", "json", "math"}
ALLOWED_MODULES = {"sys", "os", "io", "hashlib", "tokenize", "itertools"}
MATH_OK = {"gcd", "isqrt"}
FMT_FLOAT = set("eEfFgG")


def float_format_hits(text):
    hits = 0
    n = len(text)
    k = 0
    while k < n:
        if text[k] == PC:
            m = k + 1
            while m < n and text[m] in "-+ #0":
                m += 1
            while m < n and text[m].isdigit():
                m += 1
            if m < n and text[m] == chr(46):
                m += 1
                while m < n and text[m].isdigit():
                    m += 1
            if m < n and text[m] in FMT_FLOAT:
                hits += 1
        if text[k] == "{":
            close = text.find("}", k)
            if close > k:
                inner = text[k + 1:close]
                if ":" in inner and len(inner) > 0 and inner[-1] in FMT_FLOAT:
                    hits += 1
        k += 1
    return hits


def token_audit(src):
    r = {"true_div": 0, "floor_div": 0, "pow_op": 0, "kw_unpack": 0, "float_num": 0, "complex_num": 0,
         "banned_names": [], "bad_imports": [], "math_import_lines": 0, "fmt_float": 0,
         "slash_bytes": src.count(bytes([47]))}
    prev = None
    for t in tokenize.tokenize(io.BytesIO(src).readline):
        if t.type in (tokenize.ENCODING, tokenize.NL, tokenize.COMMENT):
            continue
        if t.type in (tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT):
            prev = None
            continue
        if t.type == tokenize.OP:
            s = t.string
            if s in (SL, SL + "="):
                r["true_div"] += 1
            elif s in (SL + SL, SL + SL + "="):
                r["floor_div"] += 1
            elif s == ST + ST + "=":
                r["pow_op"] += 1
            elif s == ST + ST:
                if prev is None or (prev.type == tokenize.OP and prev.string in ("(", ",", "{")) or \
                        (prev.type == tokenize.NAME and prev.string == "lambda"):
                    r["kw_unpack"] += 1
                else:
                    r["pow_op"] += 1
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if s.endswith("j"):
                r["complex_num"] += 1
            elif s.startswith("0x") or s.startswith("0o") or s.startswith("0b"):
                pass
            elif chr(46) in s or "e" in s:
                r["float_num"] += 1
        elif t.type == tokenize.NAME:
            s = t.string
            if (prev is None) and s in ("import", "from"):
                line = t.line.strip()
                words = line.replace(",", " ").split()
                if s == "import":
                    mods = [w.split(chr(46))[0] for w in words[1:] if w != "as"]
                    for m in mods:
                        if m not in ALLOWED_MODULES:
                            r["bad_imports"].append(m + "@" + str(t.start[0]))
                else:
                    mod = words[1].split(chr(46))[0] if len(words) > 1 else ""
                    names = words[3:] if len(words) > 3 else []
                    if mod == "math":
                        r["math_import_lines"] += 1
                        if not set(names) <= MATH_OK:
                            r["bad_imports"].append("math:" + ",".join(names) + "@" + str(t.start[0]))
                    elif mod not in ALLOWED_MODULES:
                        r["bad_imports"].append(mod + "@" + str(t.start[0]))
            elif s in BANNED_NAMES:
                line = t.line.strip()
                if s == "math" and line.startswith("from math import"):
                    pass
                else:
                    r["banned_names"].append(s + "@" + str(t.start[0]))
        elif t.type == tokenize.STRING:
            r["fmt_float"] += float_format_hits(t.string)
        prev = t
    r["violations"] = (r["true_div"] + r["pow_op"] + r["float_num"] + r["complex_num"] + len(r["banned_names"])
                       + len(r["bad_imports"]) + r["fmt_float"])
    return r


def section_audit():
    targets = (("m4_quaternions.py", os.path.join(BASE, "m4", "m4_quaternions.py")),
               ("m5_octonion_loop.py", os.path.join(BASE, "m5", "m5_octonion_loop.py")),
               ("refute2_m4m5.py", os.path.abspath(__file__)))
    for name, path in targets:
        with open(path, "rb") as fh:
            src = fh.read()
        a = token_audit(src)
        row("SOURCE", ("file", name), ("sha256", hashlib.sha256(src).hexdigest()), ("bytes", len(src)), ("claim", "MEASURED"))
        check("AUDIT_INTEGER_ONLY", a["violations"] == 0, ("file", name), ("true_division_ops", a["true_div"]),
              ("floor_division_ops", a["floor_div"]), ("power_ops", a["pow_op"]), ("kwargs_unpack_double_star", a["kw_unpack"]),
              ("float_literals", a["float_num"]), ("complex_literals", a["complex_num"]),
              ("banned_names", ",".join(a["banned_names"]) if a["banned_names"] else "none"),
              ("bad_imports", ",".join(a["bad_imports"]) if a["bad_imports"] else "none"),
              ("math_import_lines_gcd_isqrt_only", a["math_import_lines"]), ("float_format_specs_in_strings", a["fmt_float"]),
              ("slash_bytes_anywhere", a["slash_bytes"]))
    bad = ("x = 7 " + SL + " 2\ny = 1" + chr(46) + "5\nz = 3 " + ST + ST + " 2\nw = float(x)\nv = eval(s)\nimport fractions\n"
           + "from math import sqrt\nu = 2j\nt = '" + PC + chr(46) + "3f' \n").encode()
    b = token_audit(bad)
    control("CONTROL_AUDIT_FLAGS_PLANTED_SOURCE", b["violations"] >= 9 and b["true_div"] == 1 and b["pow_op"] == 1
            and b["float_num"] == 1 and b["complex_num"] == 1 and b["fmt_float"] == 1 and len(b["bad_imports"]) == 2,
            ("planted", "true_div+float_literal+power+float_name+eval+fractions_import+math_sqrt+complex_literal+float_format"),
            ("violations_found", b["violations"]))


# ------------------------------------------------------------------ rings
class QuadRing:
    # elements (a, b) meaning a + b w with w*w equal to p + q w
    def __init__(self, name, p, q):
        self.name = name
        self.p = p
        self.q = q
        self.zero = (0, 0)
        self.one = (1, 0)

    def add(self, u, v):
        return (u[0] + v[0], u[1] + v[1])

    def sub(self, u, v):
        return (u[0] - v[0], u[1] - v[1])

    def neg(self, u):
        return (-u[0], -u[1])

    def mul(self, u, v):
        t = u[1] * v[1]
        return (u[0] * v[0] + self.p * t, u[0] * v[1] + u[1] * v[0] + self.q * t)

    def i(self, n):
        return (n, 0)

    def half(self, u):
        a, ra = divmod(u[0], 2)
        b, rb = divmod(u[1], 2)
        if ra or rb:
            return None
        return (a, b)


class PolyRing:
    # integer polynomials: dict sorted-variable-tuple -> nonzero int
    def __init__(self):
        self.zero = {}

    def add(self, p, q):
        r = dict(p)
        for m, c in q.items():
            r[m] = r.get(m, 0) + c
        return {m: c for m, c in r.items() if c != 0}

    def neg(self, p):
        return {m: -c for m, c in p.items()}

    def sub(self, p, q):
        return self.add(p, self.neg(q))

    def mul(self, p, q):
        r = {}
        for m1, c1 in p.items():
            for m2, c2 in q.items():
                m = tuple(sorted(m1 + m2))
                r[m] = r.get(m, 0) + c1 * c2
        return {m: c for m, c in r.items() if c != 0}


class IntRing:
    def __init__(self):
        self.zero = 0
        self.one = 1

    def add(self, u, v):
        return u + v

    def sub(self, u, v):
        return u - v

    def neg(self, u):
        return -u

    def mul(self, u, v):
        return u * v


ZPHI = QuadRing("Z[phi]", 1, 1)
ZS2 = QuadRing("Z[sqrt2]", 2, 0)
PR = PolyRing()
ZR = IntRing()


def pv(i):
    return {(i,): 1}


# ------------------------------------------------------------------ quaternions by the explicit Hamilton formula (no table)
def hmul(R, x, y, o=1):
    # o=1: ij=k (Hamilton); o=-1: ij=-k (the opposite algebra)
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    m, ad, sb = R.mul, R.add, R.sub

    def cross(u, v):
        return u if o == 1 else R.neg(u)
    r0 = sb(sb(sb(m(a1, a2), m(b1, b2)), m(c1, c2)), m(d1, d2))
    r1 = ad(ad(m(a1, b2), m(b1, a2)), cross(sb(m(c1, d2), m(d1, c2)), None))
    r2 = ad(ad(m(a1, c2), m(c1, a2)), cross(sb(m(d1, b2), m(b1, d2)), None))
    r3 = ad(ad(m(a1, d2), m(d1, a2)), cross(sb(m(b1, c2), m(c1, b2)), None))
    return (r0, r1, r2, r3)


def hconj(R, x):
    return (x[0], R.neg(x[1]), R.neg(x[2]), R.neg(x[3]))


def hnorm(R, x):
    s = R.zero
    for c in x:
        s = R.add(s, R.mul(c, c))
    return s


def sprod(R, X, Y, o=1):
    # X, Y stored as 2x, 2y; returns stored 2xy or None when it is not exact in the ring
    Z = hmul(R, X, Y, o)
    out = []
    for c in Z:
        h = R.half(c)
        if h is None:
            return None
        out.append(h)
    return tuple(out)


def closure_bfs(mul, gens, ident, cap):
    S = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = mul(a, g)
                if b is None:
                    return None
                if b not in S:
                    S.add(b)
                    nxt.append(b)
                    if len(S) > cap:
                        return S
        frontier = nxt
    return S


def is_closed_set(mul, S):
    for x in S:
        for y in S:
            z = mul(x, y)
            if z is None or z not in S:
                return False
    return True


def q_int(R, v):
    return tuple(R.i(a) for a in v)


def section_quaternions():
    R = ZPHI
    one = q_int(R, (1, 0, 0, 0))
    I1 = q_int(R, (0, 1, 0, 0))
    J1 = q_int(R, (0, 0, 1, 0))
    K1 = q_int(R, (0, 0, 0, 1))
    um = lambda a, b: hmul(R, a, b)
    q8_ijk = closure_bfs(um, [I1, J1, K1], one, 100)
    q8_ij = closure_bfs(um, [I1, J1], one, 100)
    expect_q8 = set()
    for v in (I1, J1, K1, one):
        expect_q8.add(v)
        expect_q8.add(tuple(R.neg(c) for c in v))
    check("R2_Q8_GENERATED_BY_I_J_K", q8_ijk == expect_q8 and q8_ij == expect_q8 and len(q8_ijk) == 8,
          ("generated_by_i,j,k", len(q8_ijk)), ("generated_by_i,j", len(q8_ij)), ("equals_+-1,+-i,+-j,+-k", q8_ijk == expect_q8),
          ("method", "explicit_Hamilton_formula_BFS_not_a_table"))
    lit = {I1, J1, K1}
    control("CONTROL_LITERAL_I_J_K_NOT_CLOSED", not is_closed_set(um, lit), ("set", "i,j,k_without_signs_or_1"),
            ("witness", "i.i:-1_not_in_set"))
    # 2T, 2O, 2I generated from two generators each, stored scaled by 2
    s_t = q_int(R, (1, 1, 1, 1))
    i2 = q_int(R, (0, 2, 0, 0))
    one2 = q_int(R, (2, 0, 0, 0))
    t_ico = ((0, 1), (-1, 1), (1, 0), (0, 0))   # (phi + phi^-1 i + j) with phi^-1 = phi - 1, stored as 2q
    mp = lambda a, b: sprod(ZPHI, a, b)
    G2T = closure_bfs(mp, [s_t, i2], one2, 1000)
    G2I = closure_bfs(mp, [s_t, t_ico], one2, 1000)
    s_t2 = tuple(ZS2.i(1) for _ in range(4))
    u_o = ((0, 1), (0, 1), (0, 0), (0, 0))       # sqrt2 (1 + i) is 2 (1+i) over sqrt2
    one2s = (ZS2.i(2), ZS2.zero, ZS2.zero, ZS2.zero)
    ms = lambda a, b: sprod(ZS2, a, b)
    G2O = closure_bfs(ms, [s_t2, u_o], one2s, 1000)
    # independent description sets
    d2t = set()
    for a in range(4):
        for s in (2, -2):
            d2t.add(tuple(R.i(s if t == a else 0) for t in range(4)))
    for sg in itertools.product((1, -1), repeat=4):
        d2t.add(tuple(R.i(c) for c in sg))
    d2o = set(tuple(ZS2.i(c[0]) for c in x) for x in d2t)
    for a, b in itertools.combinations(range(4), 2):
        for s1 in (1, -1):
            for s2 in (1, -1):
                w = [ZS2.zero] * 4
                w[a] = (0, s1)
                w[b] = (0, s2)
                d2o.add(tuple(w))
    vals = [(0, 0), (1, 0), (-1, 1), (0, 1)]

    def ico(parity):
        out = set()
        for perm in itertools.permutations(range(4)):
            inv = sum(1 for x in range(4) for y in range(x + 1, 4) if perm[x] > perm[y])
            if inv % 2 != parity:
                continue
            for sg in itertools.product((1, -1), repeat=3):
                v = [vals[0], (sg[0] * vals[1][0], sg[0] * vals[1][1]), (sg[1] * vals[2][0], sg[1] * vals[2][1]),
                     (sg[2] * vals[3][0], sg[2] * vals[3][1])]
                w = [None] * 4
                for t in range(4):
                    w[perm[t]] = v[t]
                out.add(tuple(w))
        return out
    even = ico(0)
    odd = ico(1)
    d2i_even = d2t | even
    d2i_odd = d2t | odd
    check("R2_2T_FROM_GENERATORS", G2T is not None and G2T == d2t and len(d2t) == 24, ("order", len(G2T) if G2T else "none"),
          ("generators", "(1+i+j+k)_over_2,i"), ("equals_description_Q8_plus_(+-1+-i+-j+-k)_over_2", G2T == d2t))
    check("R2_2O_FROM_GENERATORS", G2O is not None and G2O == d2o and len(d2o) == 48, ("order", len(G2O) if G2O else "none"),
          ("generators", "(1+i+j+k)_over_2,(1+i)_over_sqrt2"), ("equals_description", G2O == d2o))
    check("R2_2I_FROM_GENERATORS", G2I is not None and len(G2I) == 120 and G2I == d2i_even, ("order", len(G2I) if G2I else "none"),
          ("generators", "(1+i+j+k)_over_2,(phi+phi^-1_i+j)_over_2"), ("equals_2T_plus_even_perms", G2I == d2i_even),
          ("equals_2T_plus_odd_perms", G2I == d2i_odd))
    # element order histograms by repeated multiplication (my code) against the m4 rows
    def hist(G, mul, ident):
        h = {}
        for g in G:
            n = 1
            y = g
            while y != ident:
                y = mul(y, g)
                n += 1
            h[n] = h.get(n, 0) + 1
        return h
    h2t = hist(G2T, mp, one2)
    h2o = hist(G2O, ms, one2s)
    h2i = hist(G2I, mp, one2)
    check("R2_ORDER_HISTOGRAMS_MATCH_M4", hist_str(h2t) == "1:1,2:1,3:8,4:6,6:8" and hist_str(h2o) == "1:1,2:1,3:8,4:18,6:8,8:12"
          and hist_str(h2i) == "1:1,2:1,3:20,4:30,5:24,6:20,10:24", ("2T", hist_str(h2t)), ("2O", hist_str(h2o)), ("2I", hist_str(h2i)))
    # icosian parity: both parity classes, under both orientations of the table
    mpo = lambda a, b: sprod(ZPHI, a, b, -1)
    ev_p = is_closed_set(mp, d2i_even)
    od_p = is_closed_set(mp, d2i_odd)
    ev_o = is_closed_set(mpo, d2i_even)
    od_o = is_closed_set(mpo, d2i_odd)
    un_p = is_closed_set(mp, d2t | even | odd)
    un_o = is_closed_set(mpo, d2t | even | odd)
    check("R2_ICOSIAN_PARITY_ORIENTATION_FREE", ev_p and od_p and ev_o and od_o and not un_p and not un_o,
          ("even_closed_ij:k", ev_p), ("odd_closed_ij:k", od_p), ("even_closed_ij:-k", ev_o), ("odd_closed_ij:-k", od_o),
          ("union_closed_either", un_p or un_o),
          ("meaning", "each parity class is a group under BOTH orientations; the parity choice does not depend on the orientation"))
    # conjugacy of the two copies inside 2O: u = w over sqrt2 with w = +-e_a +- e_b, so u x conj(u) = w x conj(w) over 2
    def conj_by(w, X):
        Y = hmul(R, hmul(R, w, X), hconj(R, w))
        out = []
        for c in Y:
            h = R.half(c)
            if h is None:
                return None
            out.append(h)
        return tuple(out)
    swappers = 0
    keepers = 0
    tried = 0
    for a, b in itertools.combinations(range(4), 2):
        for s1 in (1, -1):
            for s2 in (1, -1):
                w = [R.zero] * 4
                w[a] = R.i(s1)
                w[b] = R.i(s2)
                w = tuple(w)
                img = set(conj_by(w, X) for X in d2i_even)
                tried += 1
                if img == d2i_odd:
                    swappers += 1
                elif img == d2i_even:
                    keepers += 1
    check("R2_ICOSIAN_COPIES_CONJUGATE_IN_SU2", swappers > 0, ("elements_of_2O_minus_2T_tried", tried),
          ("conjugate_even_copy_onto_odd_copy", swappers), ("map_even_copy_to_itself", keepers),
          ("meaning", "the two parity copies of 2I are conjugate in SU(2)_by_elements_of_2O"))
    # rotation-free numeric SU(2) test of phi on all of 2I: phi(X) phi(Y) = 2 phi(XY halved), phi(X) phi(X)^dagger = 4 I, det = 4
    def phi(X):
        a, b, c, d = X
        return (((a, b), (c, d)), ((R.neg(c), d), (a, R.neg(b))))

    def gmul(x, y):
        return (R.sub(R.mul(x[0], y[0]), R.mul(x[1], y[1])), R.add(R.mul(x[0], y[1]), R.mul(x[1], y[0])))

    def gadd(x, y):
        return (R.add(x[0], y[0]), R.add(x[1], y[1]))

    def mm(A, B):
        return tuple(tuple(gadd(gmul(A[r][0], B[0][c]), gmul(A[r][1], B[1][c])) for c in range(2)) for r in range(2))

    def dag(A):
        return tuple(tuple((A[c][r][0], R.neg(A[c][r][1])) for c in range(2)) for r in range(2))

    def scal(k, A):
        return tuple(tuple((R.mul(R.i(k), A[r][c][0]), R.mul(R.i(k), A[r][c][1])) for c in range(2)) for r in range(2))
    four_I = (((R.i(4), R.zero), (R.zero, R.zero)), ((R.zero, R.zero), (R.i(4), R.zero)))
    G = sorted(G2I)
    unit_bad = 0
    det_bad = 0
    for X in G:
        P = phi(X)
        if mm(P, dag(P)) != four_I:
            unit_bad += 1
        det = (R.sub(gmul(P[0][0], P[1][1])[0], gmul(P[0][1], P[1][0])[0]), R.sub(gmul(P[0][0], P[1][1])[1], gmul(P[0][1], P[1][0])[1]))
        if det != (R.i(4), R.zero):
            det_bad += 1
    hom_bad = sum(1 for X in G for Y in G if mm(phi(X), phi(Y)) != scal(2, phi(mp(X, Y))))
    check("R2_PHI_ON_2I_IN_SU2_AND_MULTIPLICATIVE", unit_bad == 0 and det_bad == 0 and hom_bad == 0, ("elements", len(G)),
          ("phi_phi_dagger_not_4I", unit_bad), ("det_not_4", det_bad), ("pairs", len(G) * len(G)), ("hom_failures", hom_bad))
    # symbolic re-derivation with my own polynomial engine
    P = tuple(pv(t) for t in range(4))
    Q = tuple(pv(4 + t) for t in range(4))
    PQ = hmul(PR, P, Q)

    def phis(X):
        a, b, c, d = X
        return (((a, b), (c, d)), ((PR.neg(c), d), (a, PR.neg(b))))

    def pgm(x, y):
        return (PR.sub(PR.mul(x[0], y[0]), PR.mul(x[1], y[1])), PR.add(PR.mul(x[0], y[1]), PR.mul(x[1], y[0])))

    def pga(x, y):
        return (PR.add(x[0], y[0]), PR.add(x[1], y[1]))

    def pmm(A, B):
        return tuple(tuple(pga(pgm(A[r][0], B[0][c]), pgm(A[r][1], B[1][c])) for c in range(2)) for r in range(2))
    hom_sym = phis(PQ) == pmm(phis(P), phis(Q))
    NP = hnorm(PR, P)
    E = phis(P)
    detE = (PR.sub(pgm(E[0][0], E[1][1])[0], pgm(E[0][1], E[1][0])[0]), PR.sub(pgm(E[0][0], E[1][1])[1], pgm(E[0][1], E[1][0])[1]))
    norm_mult = PR.sub(hnorm(PR, PQ), PR.mul(hnorm(PR, P), hnorm(PR, Q))) == {}
    check("R2_PHI_SYMBOLIC", hom_sym and detE == (NP, {}) and norm_mult, ("phi(pq):phi(p)phi(q)", hom_sym), ("det_phi:N", detE == (NP, {})),
          ("N(pq):N(p)N(q)", norm_mult), ("variables", 8))
    Pbad = tuple(pv(t) for t in range(4))

    def phibad(X):
        a, b, c, d = X
        return (((a, b), (c, d)), ((c, d), (a, PR.neg(b))))
    control("CONTROL_PHI_WRONG_SIGN_SYMBOLIC", phibad(PQ) != pmm(phibad(Pbad), phibad(Q)), ("change", "lower_left_+c_instead_of_-c"))
    XpY = tuple(PR.add(P[t], Q[t]) for t in range(4))
    control("CONTROL_POLY_ENGINE_NORM_NOT_ADDITIVE", PR.sub(hnorm(PR, XpY), PR.add(hnorm(PR, P), hnorm(PR, Q))) != {},
            ("identity_tested", "N(p+q)-N(p)-N(q)"))
    # opposite orientation is genuinely different on Q8 (so the orientation test above is not vacuous)
    control("CONTROL_OPPOSITE_TABLE_DIFFERS", hmul(R, I1, J1, -1) != hmul(R, I1, J1, 1), ("ij_in_H", "k"), ("ij_in_H_op", "-k"))
    return {"G2T": G2T, "G2O": G2O, "G2I": G2I, "mp": mp, "ms": ms, "one2": one2, "one2s": one2s, "Q8": q8_ijk, "um": um, "one": one}


# ------------------------------------------------------------------ Jordan index on the binary groups; Dic_31 beyond 120
def cayley(G, mul):
    L = sorted(G)
    idx = {g: t for t, g in enumerate(L)}
    tab = [[idx[mul(x, y)] for y in L] for x in L]
    return L, idx, tab


def jordan_index(G, mul, ident):
    L, idx, tab = cayley(G, mul)
    n = len(L)
    e = idx[ident]
    inv = [0] * n
    for x in range(n):
        for y in range(n):
            if tab[x][y] == e:
                inv[x] = y
    classes = []
    seen = set()
    for x in range(n):
        if x in seen:
            continue
        cl = set(tab[tab[g][x]][inv[g]] for g in range(n))
        seen |= cl
        classes.append(frozenset(cl))
    ecls = [c for c in classes if e in c][0]
    others = [c for c in classes if e not in c]
    best = 1
    normals = 0
    for mask in range(1 << len(others)):
        S = set(ecls)
        for t in range(len(others)):
            if (mask >> t) & 1:
                S |= others[t]
        if all(tab[x][y] in S for x in S for y in S):
            normals += 1
            if all(tab[x][y] == tab[y][x] for x in S for y in S):
                if len(S) > best:
                    best = len(S)
    q, r = divmod(n, best)
    return n, len(classes), normals, best, q, r


def cyclotomic(m, cache):
    if m in cache:
        return cache[m]
    num = [-1] + [0] * (m - 1) + [1]
    for d in range(1, m):
        if m % d == 0:
            den = cyclotomic(d, cache)
            num = poly_exact_div(num, den)
    cache[m] = num
    return num


def poly_exact_div(num, den):
    num = list(num)
    dd = len(den) - 1
    quo = [0] * (len(num) - dd)
    for k in range(len(num) - 1, dd - 1, -1):
        c = num[k]
        if c:
            quo[k - dd] = c
            for t in range(dd + 1):
                num[k - dd + t] -= c * den[t]
    if any(num[:dd]):
        raise ArithmeticError("inexact cyclotomic division")
    while len(quo) > 1 and quo[-1] == 0:
        quo.pop()
    return quo


class CycloRing:
    def __init__(self, m):
        self.m = m
        self.Phi = cyclotomic(m, {})
        self.D = len(self.Phi) - 1
        self.zero = tuple([0] * self.D)
        self.one = tuple([1] + [0] * (self.D - 1))

    def reduce(self, poly):
        poly = list(poly)
        D = self.D
        for k in range(len(poly) - 1, D - 1, -1):
            c = poly[k]
            if c:
                for t in range(D + 1):
                    poly[k - D + t] -= c * self.Phi[t]
        poly = poly[:D] + [0] * (D - len(poly[:D]))
        return tuple(poly)

    def add(self, u, v):
        return tuple(a + b for a, b in zip(u, v))

    def neg(self, u):
        return tuple(-a for a in u)

    def sub(self, u, v):
        return tuple(a - b for a, b in zip(u, v))

    def mul(self, u, v):
        prod = [0] * (2 * self.D - 1)
        for a in range(self.D):
            if u[a]:
                for b in range(self.D):
                    if v[b]:
                        prod[a + b] += u[a] * v[b]
        return self.reduce(prod)

    def conj(self, u):
        poly = [0] * self.m
        for k, c in enumerate(u):
            poly[(self.m - k) % self.m] += c
        return self.reduce(poly)


def section_jordan(QD):
    rows = []
    for name, G, mul, ident, expect in (("Q8", QD["Q8"], QD["um"], QD["one"], 2), ("2T", QD["G2T"], QD["mp"], QD["one2"], 12),
                                        ("2O", QD["G2O"], QD["ms"], QD["one2s"], 24), ("2I", QD["G2I"], QD["mp"], QD["one2"], 60)):
        n, ncl, nnorm, best, q, r = jordan_index(G, mul, ident)
        rows.append(name + ":" + str(q))
        check("R2_JORDAN_INDEX_" + name, r == 0 and q == expect, ("order", n), ("conjugacy_classes", ncl), ("normal_subgroups", nnorm),
              ("largest_abelian_normal_subgroup", best), ("index", q), ("expect", expect),
              ("method", "all_unions_of_conjugacy_classes_tested_for_closure_and_commutativity"))
    # Dic_31 in SU(2) over Z[zeta_62]: order 124 above 120
    n = 31
    m = 2 * n
    C = CycloRing(m)
    tot = sum(1 for k in range(1, m + 1) if gcd_int(k, m) == 1)
    zeta = tuple(1 if t == 1 else 0 for t in range(C.D))
    zinv = C.conj(zeta)
    Zm = ((zeta, C.zero), (C.zero, zinv))
    Jm = ((C.zero, C.one), (C.neg(C.one), C.zero))
    Im = ((C.one, C.zero), (C.zero, C.one))

    def mm(A, B):
        return tuple(tuple(C.add(C.mul(A[r][0], B[0][c]), C.mul(A[r][1], B[1][c])) for c in range(2)) for r in range(2))

    def dag(A):
        return tuple(tuple(C.conj(A[c][r]) for c in range(2)) for r in range(2))

    def det(A):
        return C.sub(C.mul(A[0][0], A[1][1]), C.mul(A[0][1], A[1][0]))
    G = closure_bfs(mm, [Zm, Jm], Im, 2000)
    su2 = all(mm(U, dag(U)) == Im and det(U) == C.one for U in G)
    diag = [U for U in G if U[0][1] == C.zero and U[1][0] == C.zero]
    anti = [U for U in G if U[0][0] == C.zero and U[1][1] == C.zero]
    nonab = mm(Zm, Jm) != mm(Jm, Zm)
    dset = set(diag)
    d_closed = all(mm(a, b) in dset for a in diag for b in diag)
    d_ab = all(mm(a, b) == mm(b, a) for a in diag for b in diag)
    d_norm = all(mm(mm(g, a), dag(g)) in dset for g in G for a in diag)
    check("R2_DIC31_ORDER_124_IN_SU2", C.D == tot and len(G) == 124 and su2 and nonab and len(diag) == 62 and len(anti) == 62
          and d_closed and d_ab and d_norm,
          ("ring", "Z[zeta_62]"), ("Phi_degree", C.D), ("totient", tot), ("order", len(G)), ("exceeds_120", len(G) > 120),
          ("all_unitary_det_1", su2), ("nonabelian", nonab), ("diagonal", len(diag)), ("antidiagonal", len(anti)),
          ("diagonal_subgroup_abelian_normal_index_2", d_closed and d_ab and d_norm),
          ("meaning", "a finite subgroup of SU(2)_of_order_124_exists;_it_lies_in_N(T)_so_it_is_toral_type"))
    Jbad = ((C.zero, C.one), (C.one, C.zero))
    Gb = closure_bfs(mm, [Zm, Jbad], Im, 2000)
    control("CONTROL_DIC31_WRONG_J_NOT_SU2", not all(mm(U, dag(U)) == Im and det(U) == C.one for U in Gb),
            ("change", "J:[[0,1],[1,0]]_det_-1"))
    row("R2_JORDAN_INDICES", ("values", ",".join(rows) + ",Dic31:2"), ("max", 60), ("attained_by", "2I"), ("claim", "MEASURED"))


def gcd_int(a, b):
    while b:
        a, b = b, a % b
    return a


# ------------------------------------------------------------------ octonions: Fano table (m5 convention, rebuilt) and Cayley-Dickson doubling
LINES = ((1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3))


def fano_table(lines):
    T = [[None] * 8 for _ in range(8)]
    for a in range(8):
        T[0][a] = (1, a)
        T[a][0] = (1, a)
    for a in range(1, 8):
        T[a][a] = (-1, 0)
    for (a, b, c) in lines:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            T[x][y] = (1, z)
            T[y][x] = (-1, z)
    return T


def tmul(T, x, y):
    z = [0] * 8
    for i in range(8):
        if x[i]:
            for j in range(8):
                if y[j]:
                    s, k = T[i][j]
                    z[k] += s * x[i] * y[j]
    return tuple(z)


def cd_mul(R, x, y, conv):
    a, b = x
    c, d = y
    if conv == "wiki":
        # (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c))
        return (tuple(R.sub(u, v) for u, v in zip(hmul(R, a, c), hmul(R, hconj(R, d), b))),
                tuple(R.add(u, v) for u, v in zip(hmul(R, d, a), hmul(R, b, hconj(R, c)))))
    # baez: (a,b)(c,d) = (ac - d conj(b), conj(a) d + c b)
    return (tuple(R.sub(u, v) for u, v in zip(hmul(R, a, c), hmul(R, d, hconj(R, b)))),
            tuple(R.add(u, v) for u, v in zip(hmul(R, hconj(R, a), d), hmul(R, c, b))))


def cd_table(conv):
    def bvec(t):
        v = [0] * 8
        v[t] = 1
        return (tuple(v[:4]), tuple(v[4:]))
    T = [[None] * 8 for _ in range(8)]
    ok = True
    for a in range(8):
        for b in range(8):
            p = cd_mul(ZR, bvec(a), bvec(b), conv)
            flat = p[0] + p[1]
            nz = [t for t in range(8) if flat[t] != 0]
            if len(nz) != 1 or flat[nz[0]] not in (1, -1):
                ok = False
                T[a][b] = (0, 0)
            else:
                T[a][b] = (flat[nz[0]], nz[0])
    return T, ok


def assoc_vec(T, i, j, k):
    ei = tuple(1 if t == i else 0 for t in range(8))
    ej = tuple(1 if t == j else 0 for t in range(8))
    ek = tuple(1 if t == k else 0 for t in range(8))
    l = tmul(T, tmul(T, ei, ej), ek)
    r = tmul(T, ei, tmul(T, ej, ek))
    return tuple(l[t] - r[t] for t in range(8))


def signed_perm_isos(TA, TB):
    # all maps e_a -> s_a f_p(a) (a = 1..7, e0 -> f0) with g(e_a e_b) = g(e_a) g(e_b) for all a, b in 1..7
    count = 0
    perms_ok = 0
    pairs = [(a, b) for a in range(1, 8) for b in range(1, 8) if a != b]
    for p in itertools.permutations(range(1, 8)):
        P = (0,) + p
        if any(TB[P[a]][P[b]][1] != P[TA[a][b][1]] for (a, b) in pairs):
            continue
        perms_ok += 1
        for sb in range(128):
            S = [1] + [(-1 if (sb >> (a - 1)) & 1 else 1) for a in range(1, 8)]
            good = True
            for (a, b) in pairs:
                sig, c = TA[a][b]
                tau, d = TB[P[a]][P[b]]
                if sig * S[c] != S[a] * S[b] * tau:
                    good = False
                    break
            if good:
                count += 1
    return count, perms_ok


def loop_of(T):
    def uv(u):
        v = [0] * 8
        v[u % 8] = 1 if u < 8 else -1
        return tuple(v)
    M = [[None] * 16 for _ in range(16)]
    for u in range(16):
        for w in range(16):
            v = tmul(T, uv(u), uv(w))
            nz = [t for t in range(8) if v[t] != 0]
            M[u][w] = nz[0] if v[nz[0]] == 1 else nz[0] + 8
    return M


def perm_group_order(gens, npts, cap):
    ident = tuple(range(npts))
    S = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = tuple(g[a[x]] for x in range(npts))
                if b not in S:
                    S.add(b)
                    nxt.append(b)
                    if len(S) > cap:
                        return len(S)
        frontier = nxt
    return len(S)


def section_octonions(QD):
    TF = fano_table(LINES)
    TW, okw = cd_table("wiki")
    TBz, okb = cd_table("baez")
    unit_ok = all(T[0][a] == (1, a) and T[a][0] == (1, a) for T in (TW, TBz) for a in range(8))
    diag_ok = all(T[a][a] == (-1, 0) for T in (TW, TBz) for a in range(1, 8))
    check("R2_CD_TABLES_ARE_SIGNED_UNIT_TABLES", okw and okb and unit_ok and diag_ok, ("wiki_convention", okw), ("baez_convention", okb),
          ("f0_is_unit", unit_ok), ("f_a_squared_is_-1", diag_ok))
    # alternativity (associator antisymmetry on all 512 basis triples) and norm multiplicativity for the CD products
    for name, T, conv in (("wiki", TW, "wiki"), ("baez", TBz, "baez")):
        anti = sum(1 for i in range(8) for j in range(8) for k in range(8)
                   if any(x + y for x, y in zip(assoc_vec(T, i, j, k), assoc_vec(T, j, i, k)))
                   or any(x + y for x, y in zip(assoc_vec(T, i, j, k), assoc_vec(T, i, k, j))))
        X = (tuple(pv(t) for t in range(4)), tuple(pv(4 + t) for t in range(4)))
        Y = (tuple(pv(8 + t) for t in range(4)), tuple(pv(12 + t) for t in range(4)))
        XY = cd_mul(PR, X, Y, conv)

        def n8(Z):
            return PR.add(hnorm(PR, Z[0]), hnorm(PR, Z[1]))
        diff = PR.sub(n8(XY), PR.mul(n8(X), n8(Y)))
        check("R2_CD_OCTONIONS_" + name.upper(), anti == 0 and diff == {}, ("associator_antisymmetry_failures_on_512_triples", anti),
              ("N(xy)-N(x)N(y)_monomials", len(diff)), ("variables", 16))
    # isomorphism census: Fano table (m5) to each CD table, and automorphism recounts
    iso_w, pw = signed_perm_isos(TF, TW)
    iso_b, pb = signed_perm_isos(TF, TBz)
    aut_f, af = signed_perm_isos(TF, TF)
    aut_w, aw = signed_perm_isos(TW, TW)
    check("R2_FANO_TABLE_IS_CAYLEY_DICKSON_OCTONIONS", iso_w == 1344 and iso_b == 1344 and pw == 168 and pb == 168,
          ("candidates_covered_each", 645120), ("isomorphisms_fano_to_cd_wiki", iso_w), ("isomorphisms_fano_to_cd_baez", iso_b),
          ("line_preserving_point_maps", str(pw) + "," + str(pb)),
          ("meaning", "m5_table_is_isomorphic_to_the_doubled_quaternions;_the_iso_set_is_a_torsor_of_Aut"))
    check("R2_FRAME_AUT_RECOUNT_1344", aut_f == 1344 and aut_w == 1344 and af == 168 and aw == 168,
          ("aut_fano_table", aut_f), ("aut_cd_table", aut_w), ("point_permutations", str(af) + "," + str(aw)))
    broken = ((2, 1, 4),) + LINES[1:]
    TX = fano_table(broken)
    iso_x, px = signed_perm_isos(TX, TW)
    control("CONTROL_BROKEN_FANO_NOT_ISOMORPHIC_TO_O", iso_x == 0, ("change", "line_(1,2,4)_reversed"), ("isomorphisms_found", iso_x))
    # the loop of 16 signed units in the CD table
    M = loop_of(TW)
    L = list(range(16))
    clo = all(M[x][y] is not None for x in L for y in L)
    ident = all(M[0][u] == u and M[u][0] == u for u in L)
    inv = all(sum(1 for w in L if M[u][w] == 0 and M[w][u] == 0) == 1 for u in L)
    latin = all(sorted(M[u]) == L for u in L) and all(sorted(M[u][w] for u in L) == L for w in L)
    nas = sum(1 for x in L for y in L for z in L if M[M[x][y]][z] != M[x][M[y][z]])
    mf = [0, 0, 0]
    for x in L:
        for y in L:
            for z in L:
                if M[z][M[x][M[z][y]]] != M[M[M[z][x]][z]][y]:
                    mf[0] += 1
                if M[x][M[z][M[y][z]]] != M[M[M[x][z]][y]][z]:
                    mf[1] += 1
                if M[M[z][x]][M[y][z]] != M[M[z][M[x][y]]][z]:
                    mf[2] += 1
    check("R2_CD_LOOP_NOT_A_GROUP", clo and ident and inv and latin and nas == 1344 and mf == [0, 0, 0],
          ("closed", clo), ("identity", ident), ("unique_two_sided_inverses", inv), ("latin", latin),
          ("nonassociative_ordered_triples", nas), ("associative", 4096 - nas), ("moufang_failures", ",".join(str(f) for f in mf)),
          ("meaning", "associativity_is_the_only_group_axiom_that_fails"))
    seven = set(range(1, 8))
    control("CONTROL_LITERAL_SEVEN_UNITS_NOT_CLOSED", not all(M[x][y] in seven for x in seven for y in seven),
            ("set", "e1..e7_without_signs_or_1"), ("witness", "e1.e1:-e0"))
    q = [[M[a][b] % 8 for b in range(8)] for a in range(8)]
    qgroup = all(q[q[a][b]][c] == q[a][q[b][c]] for a in range(8) for b in range(8) for c in range(8)) and \
        all(q[a][b] == q[b][a] for a in range(8) for b in range(8)) and all(q[a][a] == 0 for a in range(8))
    check("R2_UNITS_MOD_SIGN_FORM_A_GROUP", qgroup, ("structure", "elementary_abelian_order_8"),
          ("meaning", "forgetting_signs_the_8_units_1,e1..e7_DO_form_a_group;_the_obstruction_is_the_sign_cocycle"))
    # multiplication groups: operators x -> u x generate a group even when the units do not
    lg = [tuple(M[u][x] for x in L) for u in L]
    rg = [tuple(M[x][u] for x in L) for u in L]
    lm = perm_group_order(lg, 16, 100000)
    rm = perm_group_order(rg, 16, 100000)
    mlt = perm_group_order(lg + rg, 16, 100000)
    lbad = sum(1 for u in L for w in L if tuple(lg[u][lg[w][x]] for x in L) != lg[M[u][w]])
    q8l = sorted(QD["Q8"])
    qi = {g: t for t, g in enumerate(q8l)}
    qlg = [tuple(qi[QD["um"](u, x)] for x in q8l) for u in q8l]
    qlm = perm_group_order(qlg, 8, 100000)
    qbad = sum(1 for u in range(8) for w in range(8) if tuple(qlg[u][qlg[w][x]] for x in range(8)) != qlg[qi[QD["um"](q8l[u], q8l[w])]])
    check("R2_MULTIPLICATION_GROUPS", qlm == 8 and qbad == 0 and lbad > 0 and lm > 16,
          ("Q8_left_mult_group_order", qlm), ("Q8_pairs_with_L_uL_w_not_L_uw", qbad),
          ("octonion_loop_left_mult_group_order", lm), ("right_mult_group_order", rm), ("full_mult_group_order", mlt),
          ("octonion_loop_pairs_with_L_uL_w_not_L_uw", lbad),
          ("meaning", "for_Q8_the_operators_are_the_group_itself;_for_the_octonion_loop_they_generate_a_larger_group"))
    # the 240 unit integral octonions: test the Kirmse quad system and its seven (0 k) transposes for closure
    quads = [frozenset((0,) + l) for l in LINES] + [frozenset(set(range(1, 8)) - set(l)) for l in LINES]

    def units240(qs):
        S = set()
        for a in range(8):
            for s in (2, -2):
                S.add(tuple(s if t == a else 0 for t in range(8)))
        for qd in qs:
            idxs = sorted(qd)
            for sg in itertools.product((1, -1), repeat=4):
                v = [0] * 8
                for t, ix in enumerate(idxs):
                    v[ix] = sg[t]
                S.add(tuple(v))
        return S

    def smul(X, Y):
        Z = tmul(TF, X, Y)
        out = []
        for c in Z:
            h, r = divmod(c, 2)
            if r:
                return None
            out.append(h)
        return tuple(out)

    def closed_fast(S):
        for x in S:
            for y in S:
                z = smul(x, y)
                if z is None or z not in S:
                    return False
        return True
    results = []
    closed_sets = []
    for k in range(8):
        def sw(i):
            if k == 0:
                return i
            return k if i == 0 else (0 if i == k else i)
        qs = [frozenset(sw(i) for i in qd) for qd in quads]
        S = units240(qs)
        ok = len(S) == 240 and closed_fast(S)
        results.append(str(k) + ":" + ("closed" if ok else "open"))
        if ok:
            closed_sets.append((k, S))
    control("CONTROL_KIRMSE_UNTRANSPOSED_NOT_CLOSED", results[0] == "0:open", ("quad_system", "{0}+lines_and_complements"),
            ("ref", "Coxeter_1946_Integral_Cayley_numbers_Duke_Math_J_13_561-578_(NAMED_not_opened)"))
    if closed_sets:
        k, S = closed_sets[0]
        norms = all(sum(c * c for c in x) == 4 for x in S)
        conj_in = all(tuple([x[0]] + [-c for c in x[1:]]) in S for x in S)
        has_L = all(tuple(s if t == a else 0 for t in range(8)) in S for a in range(8) for s in (2, -2))
        e1, e2, e3 = (0, 2, 0, 0, 0, 0, 0, 0), (0, 0, 2, 0, 0, 0, 0, 0), (0, 0, 0, 2, 0, 0, 0, 0)
        nonassoc = smul(smul(e1, e2), e3) != smul(e1, smul(e2, e3))
        check("R2_OCTAVIAN_240_LOOP", norms and conj_in and has_L and nonassoc, ("systems_tested", ",".join(results)),
              ("closed_systems", len(closed_sets)), ("first_closed_transposition_k", k), ("order", len(S)), ("unit_norm", norms),
              ("conjugate_inverses_present", conj_in), ("contains_L", has_L), ("nonassociative", nonassoc),
              ("meaning", "a_finite_nonassociative_Moufang_subloop_of_S^7_of_order_240_(octonion_side_analogue_of_2I_in_S^3)"))
    else:
        check("R2_OCTAVIAN_240_LOOP", False, ("systems_tested", ",".join(results)), ("closed_systems", 0),
              ("status", "CANNOT_CONSTRUCT_HERE_not_absence"))
    return TF


# ------------------------------------------------------------------ verdict rows
def section_verdicts():
    A = [
        ("M4.PROOF_H_ASSOCIATIVE", "SOUND", "trilinear_associator_vanishing_on_64_basis_triples_vanishes_identically;_valid_over_any_commutative_coefficient_ring"),
        ("M4.PROOF_H_NORM_MULTIPLICATIVE", "SOUND", "full_symbolic_expansion_(Euler_four-square);_rechecked_by_R2_PHI_SYMBOLIC_with_independent_engine"),
        ("M4.PROOF_PHI_HOMOMORPHISM", "SOUND", "phi_is_R-linear_so_both_sides_are_R-bilinear;_16_basis_pairs_suffice;_rechecked_numerically_on_all_14400_pairs_of_2I"),
        ("M4.PROOF_S3_ISOMORPHIC_SU2", "SOUND", "surjectivity_via_U^-1:adj(U):U^dagger_is_complete;_injective;_multiplicative;_homeomorphism_(compact_to_Hausdorff)"),
        ("M4.Q8", "SOUND", "exhaustive_over_8_elements;_R2_rederives_Q8_from_{i,j,k}_and_{i,j}_by_the_explicit_Hamilton_formula"),
        ("M4.PROOF_FINITE_CLOSED_SUBSET_IS_SUBGROUP", "SOUND", "cancellation_in_a_group;_R2_regenerates_2T,2O,2I_from_two_generators_each_and_gets_the_same_sets_and_order_histograms"),
        ("M4.CHAIN_AND_LAGRANGE", "SOUND", "120_mod_48:24_so_no_subgroup_of_order_48_in_2I;_ladder_branches_at_2T"),
        ("M4.INTERSECTION_2O_2I_IS_2T", "SOUND_WITH_NAMED", "needs_1,sqrt2,phi,sqrt2*phi_linearly_independent_over_Q_(NAMED_degree_4)"),
        ("M4.ICOSIAN_PARITY_NOTE", "WRONG_INFERENCE", "measured_facts_correct,_but_'the_even-permutation_rule_depends_on_the_orientation'_does_not_follow:_closure_under_all_ordered_products_is_orientation-free_(x.op_y:yx);_R2_ICOSIAN_PARITY_ORIENTATION_FREE_measures_both_parity_classes_closed_under_both_tables,_and_R2_ICOSIAN_COPIES_CONJUGATE_IN_SU2_shows_they_are_conjugate"),
        ("M4.PROOF_ROTATION_IMAGES", "SOUND_MINOR_GAP", "kernel_{+-1}_on_all_of_S^3_cites_only_the_i-commutation_identity;_the_j_case_follows_by_the_cyclic_automorphism_i->j->k;_inside_each_finite_group_the_kernel_was_measured_directly"),
        ("M4.PROOF_DIC_N_IN_NT_ALL_N", "SOUND", "two_polynomial_identities_give_all_n"),
        ("M4.PROOF_EXCEPTIONAL_NOT_IN_ANY_NT", "SOUND", "index<:2_abelian_subgroup_would_need_12,24,60_elements_vs_max_noncentral_centralisers_6,8,10"),
        ("M4.PROOF_SU2_CEILING", "SOUND_CONCLUSION_GAP_IN_TEXT", "order>120_implies_cyclic_or_Dic_n_correctly_derived_from_Klein;_but_the_2O_maximality_text_says_'any_larger_order_puts_G_in_some_N(T)'_and_a_G_of_order_96_is_not_covered_by_the_>120_lemma;_closes_via_Klein_directly_(order_96_is_cyclic_or_Dic_24)"),
        ("M4.PROOF_SU2_TOWER_CANNOT_REACH_SU2", "SOUND_WITH_NAMED", "squares_commute_on_the_closure;_squaring_onto_SU(2);_relies_on_Klein_and_maximal_torus"),
        ("M4.PROOF_JORDAN_TOWER_GENERAL_N", "SOUND_WITH_NAMED_MISSING_CITATION", "x^M_in_A_k_with_M:J(n)!_correct;_'H_mod_H0_finite'_also_needs_Cartan_closed-subgroup_theorem_(H_is_a_Lie_group);_classical_source_of_the_corollary:_Turing_1938_Ann_Math_39_105-111_(NAMED_not_opened)"),
        ("M5.P1_P1B", "SOUND", "zero_polynomials;_R2_rederives_N(xy):N(x)N(y)_for_two_independent_Cayley-Dickson_products"),
        ("M5.P2_L_NOT_A_GROUP", "SOUND_AIRTIGHT", "explicit_triple_(e1e2)e3:-e6,e1(e2e3):+e6_recomputed;_closure,identity,inverses_hold;_R2_CD_LOOP_NOT_A_GROUP_finds_1344_nonassociative_triples_in_an_independent_table"),
        ("M5.P3_MOUFANG", "SOUND", "the_three_are_the_standard_left,right,middle_Moufang_identities;_0_failures_recomputed_in_the_CD_loop"),
        ("M5.P4_P5_S7", "SOUND", "polarization_correct:_2F(z,z):sum_z_i_z_j_(F(e_i,e_j)+F(e_j,e_i))"),
        ("M5.P6_ALTERNATIVE", "SOUND", "transposition_antisymmetry_is_a_trilinear_identity;_repeated_indices_included"),
        ("M5.P7_P8_P8B", "SOUND", "1344_recounted_independently_on_the_Fano_table_and_on_the_CD_table;_GL(3,2)_via_F2^3_labeling_valid;_non-split_argument_valid_(generating_pair_checked_to_generate_168)"),
        ("M5.P9_P10", "SOUND_UPGRADABLE", "surjectivity_of_rho_onto_SU(2)_tagged_NAMED_in_m5_is_MEASURED_in_m4_(rho_is_m4's_phi_with_i,j,k:e_a,e_b,e_c)"),
        ("M5.P11_PARALLEL", "SOUND", "G2:Aut(O)_parallels_SO(3):Aut(H),_not_SU(2)_and_not_the_finite_groups;_agrees_with_this_refuter"),
        ("M5.NAMED_FRAME_GROUP_MAXIMAL_IN_G2", "NAMED_UNCHECKED", "maximality_of_2^3.GL(3,2)_in_G2_not_measured_by_anyone_in_this_run"),
    ]
    for tid, verdict, why in A:
        row("ATTACK", ("target", tid), ("verdict", verdict), ("why", why), ("claim", "MEASURED+PROOF_REVIEW"))
    proof("R2_PARITY_IS_ORIENTATION_FREE", "MEASURED+PROOF",
          "whether a finite set of unit quaternions is a group does not depend on the orientation of the table",
          "the opposite algebra has x.op_y:yx;_a_set_closed_under_xy_for_all_ordered_pairs_is_closed_under_yx;_so_both_parity_classes_are_groups_in_both_orientations",
          "R2_ICOSIAN_PARITY_ORIENTATION_FREE,R2_ICOSIAN_COPIES_CONJUGATE_IN_SU2,CONTROL_OPPOSITE_TABLE_DIFFERS")
    proof("R2_PURE_UNITS_GENERATE_SU2", "PROOF_ELEMENTARY",
          "the unit imaginary quaternions (all of S^2, not just i,j,k) generate all of S^3:SU(2)",
          "for q:a+w_unit_choose_unit_imaginary_u_perpendicular_to_w;_v:-uq:-au-u_cross_w_is_imaginary_of_norm_1_and_uv:q",
          "elementary_algebra_using_uv:-(u.v)+u_cross_v_(m4_H_IMAG_PRODUCT_POLY)")
    proof("R2_BECAUSE_IS_EXACT", "MEASURED+PROOF",
          "the 16 signed units form a group iff O is associative",
          "L_is_closed_with_identity_and_inverses;_if_L_were_associative_all_512_basis_triples_associate_and_by_trilinearity_O_would_be_associative;_conversely_O_associative_gives_L_associative",
          "R2_CD_LOOP_NOT_A_GROUP,M5_P2,M5_P6")
    proof("R2_FANO_TABLE_IS_O_WITHOUT_HURWITZ", "MEASURED+PROOF",
          "the m5 table defines the octonions",
          "1344_signed-permutation_algebra_isomorphisms_onto_each_Cayley-Dickson_double_of_H_were_found;_isomorphism_on_basis_pairs_extends_bilinearly",
          "R2_FANO_TABLE_IS_CAYLEY_DICKSON_OCTONIONS,R2_CD_OCTONIONS_WIKI,R2_CD_OCTONIONS_BAEZ")
    proof("R2_JORDAN_CONSTANT_SU2_IS_60", "MEASURED+NAMED_Klein",
          "every finite subgroup of SU(2) has an abelian normal subgroup of index at most 60; only conjugates of 2I need 60",
          "by_Klein_(NAMED)_a_finite_subgroup_is_conjugate_to_cyclic_(index_1),Dic_n_(cyclic_index_2),2T_(12),2O_(24),2I_(60);_the_four_binary_indices_were_measured_by_enumerating_all_normal_subgroups",
          "R2_JORDAN_INDEX_Q8,2T,2O,2I;DIC_n_in_m4;R2_DIC31_ORDER_124_IN_SU2")
    proof("R2_HAUSDORFF_STRENGTHENING", "PROOF+NAMED_Jordan,maximal_torus",
          "no sequence of finite subgroups of U(n), nested or not, converges in Hausdorff distance to SU(N) for N>:2",
          "with_M:J(n)!_every_x_in_G_k_has_x^M_in_an_abelian_A_k_so_x^M_and_y^M_commute;_if_G_k->SU(N)_approximate_x,y_by_x_k,y_k;_limit_gives_x^My^M:y^Mx^M_on_SU(N);_the_M-th_power_map_is_onto_(maximal_tori)_so_SU(N)_would_be_abelian",
          "argument_only;_SU(2)_inputs_measured_in_m4")
    S = [
        ("S1", "Quaternions give SU(2)", "PROVEN_MEASURED",
         "the_unit_quaternions_S^3_form_a_group_isomorphic_to_SU(2)_via_phi(a+bi+cj+dk):[[a+bi,c+di],[-c+di,a-bi]]_(m4_PROOF_S3_ISOMORPHIC_SU2_complete;_R2_PHI_SYMBOLIC_and_R2_PHI_ON_2I_recheck);_also_H^x:R_{>0}_x_SU(2)_and_Aut(H):SO(3):SU(2)_mod_{+-1}_(Skolem-Noether_NAMED);_note_the_paired_phrase_'octonions_give_G2'_uses_Aut,_whose_quaternion_counterpart_is_SO(3),_not_SU(2)"),
        ("S2", "The imaginary units of the quaternions generate the group Q8", "PROVEN_MEASURED",
         "precise:_the_basis_units_i,j,k_(any_two_suffice)_generate_Q8:{+-1,+-i,+-j,+-k}_of_order_8;_the_three_units_alone_are_not_a_group_(i.i:-1);_if_'imaginary_units'_means_all_unit_imaginary_quaternions_(S^2)_they_generate_all_of_SU(2),_not_Q8_(R2_PURE_UNITS_GENERATE_SU2)"),
        ("S3", "The seven imaginary octonion units do not form a group because octonions are not associative", "PROVEN_MEASURED",
         "precise:_the_16_signed_units_+-1,+-e1..+-e7_are_closed,_have_identity_and_unique_two-sided_inverses,_form_a_Moufang_loop,_and_fail_only_associativity_((e1e2)e3:-e6,e1(e2e3):+e6;_1344_of_4096_triples);_'because'_is_exact_by_trilinearity_(R2_BECAUSE_IS_EXACT);_the_7_units_without_signs_and_1_fail_closure_already_(same_as_i,j,k);_modulo_sign_the_8_units_DO_form_the_group_Z2^3;_every_two_units_generate_a_Q8"),
        ("S4", "G2, the symmetry group of the octonions, takes the role the finite groups play on the quaternion side", "ANALOGY_ONLY",
         "theorem_part:_G2:Aut(O)_(definition,_Cartan_1914_NAMED);_the_1344_frame_symmetries_are_in_G2_and_equal_Aut(L)_(MEASURED_m5_and_R2_recount);_line_stabiliser_restricts_onto_the_24_quaternion_frame_symmetries_(MEASURED_m5);_interpretation_part:_'takes_the_role_of_the_finite_groups'_is_not_a_theorem;_measured_structure_says_G2:Aut(O)_plays_the_role_of_SO(3):Aut(H);_the_role_of_SU(2):S^3_falls_to_S^7,_only_a_Moufang_loop;_the_role_of_Q8,2T,2O,2I_(finite_subgroups_of_S^3)_falls_to_finite_subloops_of_S^7_(L_of_order_16,_and_the_240-unit_octavian_loop_if_R2_OCTAVIAN_240_LOOP_passed);_the_role_of_V4,A4,S4,A5_in_SO(3)_falls_to_finite_subgroups_of_G2_(e.g._2^3.GL(3,2)_order_1344_MEASURED;_G2(2)_order_12096_as_Aut_of_the_octavians_NAMED_Conway-Smith_2003)"),
        ("S5", "2I is the largest group on the ladder", "WRONG",
         "as_stated_it_fails_on_both_natural_readings:_(a)_inclusion:_2O_is_not_inside_2I_(120_mod_48:24),_so_2I_is_a_maximal_element_of_a_branching_ladder,_not_its_top;_(b)_size_among_finite_subgroups_of_SU(2):_Dic_31_has_order_124_(R2_DIC31_ORDER_124_IN_SU2)_and_cyclic_groups_have_every_order;_only_the_trivial_reading_'largest_order_among_Q8,2T,2O,2I'_holds;_EXACT_TRUE_STATEMENT_(PROVEN_NAMED_Klein_+_MEASURED):_2I_(order_120)_is_the_largest_finite_subgroup_of_SU(2)_not_contained_in_a_conjugate_of_N(T):T_union_Tj;_it_and_2O_are_maximal_finite_subgroups;_every_finite_subgroup_of_order_>120_is_cyclic_or_binary_dihedral,_with_a_cyclic_subgroup_of_index_<:2;_and_2I_is_the_unique_type_attaining_the_SU(2)_Jordan_index_60_(R2_JORDAN_CONSTANT_SU2_IS_60)"),
        ("S6", "Jordan's theorem: no tower of finite subgroups converges to SU(N); it can only approach a torus extended by a finite group", "PROVEN_NAMED",
         "this_is_a_corollary_of_Jordan_(bounded-index_abelian_normal_subgroups),_not_Jordan's_statement_itself;_derivation_in_m4_PROOF_JORDAN_TOWER_GENERAL_N_is_sound_(add_Cartan_closed-subgroup_theorem);_SU(2)_case_measured_inputs_+_Klein;_holds_even_for_non-nested_Hausdorff-convergent_sequences_(R2_HAUSDORFF_STRENGTHENING);_classical_source_Turing_1938_NAMED"),
    ]
    for sid, text, status, why in S:
        row("STATEMENT", ("id", sid), ("text", text), ("status", status), ("why", why))


# ------------------------------------------------------------------ main
def main():
    row("SCRIPT", ("name", "refute2_m4m5.py"), ("python", ".".join(str(v) for v in sys.version_info[:3])), ("isolated_flag", sys.flags.isolated),
        ("E", 0), ("network", 0), ("float_used", 0), ("claim", "MEASURED"))
    section_audit()
    QD = section_quaternions()
    section_jordan(QD)
    section_octonions(QD)
    section_verdicts()
    row("SUMMARY", ("checks", COUNT["checks"]), ("controls", COUNT["controls"]), ("failures", len(FAILS)),
        ("failed_tags", ",".join(FAILS) if FAILS else "none"), ("result", "PASS" if not FAILS else "FAIL"), ("E", 0), ("claim", "MEASURED"))


def write_out():
    data = ("\n".join(ROWS) + "\n").encode()
    out = os.path.join(HERE, OUT_NAME)
    with open(out, "wb") as fh:
        fh.write(data)
    sha = hashlib.sha256(data).hexdigest()
    with open(out + ".sha256", "wb") as fh:
        fh.write((sha + "  " + OUT_NAME + "\n").encode())
    return sha


if __name__ == "__main__":
    rc = 0
    try:
        main()
    except Exception as ex:
        FAILS.append("EXCEPTION")
        row("ERROR", ("result", "FAIL"), ("kind", type(ex).__name__), ("detail", str(ex)[:300]), ("claim", "MEASURED"))
        rc = 3
    sha = write_out()
    for r in ROWS:
        print(r)
    print("HBP_SHA256=" + sha)
    if FAILS:
        for t in FAILS:
            print("FAIL|tag=" + t + "|json=0")
        if rc == 0:
            rc = 1
    sys.exit(rc)
