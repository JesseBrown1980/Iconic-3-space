#!/usr/bin/env python3
"""
measure_ladder.py: MEASUREMENT of the SU(2) ladder. C_n and Dic_n for n = 2..30, then 2T, 2O and 2I, 61 groups.

Every group is built element by element from its definition, by closure under its generators.
  C_n   = < diag(z, z^-1) >, z = x in Z[x]/Phi_n, checked to be a primitive n-th root of unity.
  Dic_n = < diag(w, w^-1), [[0,1],[-1,0]] >, w = x in Z[x]/Phi_2n, checked to be a primitive 2n-th root of unity.
  2T, 2O, 2I = unit quaternions with every coordinate scaled by 2, over Z, Z[sqrt2] = Z[x]/(x^2-2) and Z[phi] = Z[x]/(x^2-x-1).
    Generators: s = (1,1,1,1) for all three, with i = (0,2,0,0) for 2T, (sqrt2,sqrt2,0,0) for 2O, (phi,1,phi-1,0) for 2I.
    Each quaternion is also written as the 2x2 matrix [[A+Bi, C+Di], [-C+Di, A-Bi]] = 2U over Z[zeta_4], Z[zeta_8], Z[zeta_20],
    through sqrt2 = zeta_8 + zeta_8^7 and phi = 1 + zeta_20^4 + zeta_20^16, each image checked against its minimal polynomial.
Checks on every group: order (n, 4n, 24, 48, 120), identity present, closed under the generators and under every pair, det = 1, U U^dagger = I.
Method A, characters: the trace of each element in its ring, spin characters by chi_0 = 1, chi_1/2 = tr, chi_(j+1/2) = tr chi_j - chi_(j-1/2),
  evaluated once per distinct trace value and weighted by how many elements carry it. m(j) = (1/|G|) sum_g chi_j(g) for 2j = 0..120 and
  M_n = (1/|G|) sum_g tr^n for n = 0..60; each sum is required to be an ordinary integer and divisible by |G|.
Method B, invariants: for d = 0..24 the dimension of the degree-d polynomials in (x, y) fixed by every generator, as (d+1) minus the rank of the
  stacked matrices (rho_d(g) - I), rank by fraction-free Bareiss elimination over Z[zeta_m], each division exact through the field norm.
  It is required to equal m(d/2) from method A.
2D gauge: z(beta) = sum_n beta^n M_n / (n! 2^n), f = log z through beta^30 as reduced integer pairs, checked by exp(f) = z, and f_G - f_SU(2)
  with the SU(2) Haar moments Catalan(n/2) for even n and 0 for odd n, Catalan from its own recurrence.
Controls built to fail: 2I with one element dropped, one perturbed 2I trace, method B for 2I with one generator, a non-unitary matrix,
  a perturbed log coefficient, and a perturbed frozen value in the comparison.
Integer only: no float, rationals as reduced (num, den) pairs written num:den. json=0. E=0: no network, nothing written outside OUT_DIR.
Usage: python3 -I -B measure_ladder.py [OUT_DIR [PREDICTION_HBP]]
Exit 0 when every internal check holds and every control fails as built to; exit 1 otherwise.
A disagreement with the frozen prediction is a result: it is written as agree=0 and leaves the exit code alone.
Tag: every value in the output is MEASURED by this script.
"""
import sys
import os
import re
import hashlib
from math import gcd

# The measured rows were written by this script, and their digest recorded, before the prediction file was opened.
# The comparison stage below was added afterwards; it reads PREDICTION.hbp only after every measured row is emitted,
# and requires the measured block to be byte-identical to that first run.
PRECOMPARE_MEASURED_BLOCK = "268d55002ee6c15d7fc1ced73638764f4f9298347436f268d5c0ecd0b785a697"
PRECOMPARE_FILE = "6098cc3459963ba96317f2b85ee4e514173fa3c467170b84f6abbf79d6008c65"
PREDICTION_SHA_PREFIX = "c2601cf0"

TWOJ_MAX = 120
MOMENT_MAX = 60
INV_DMAX = 24
BETA_MAX = 30
NMAX = 30
CLOSURE_CAP = 4096

FAILS = []


class CheckFail(Exception):
    pass


def need(cond, what):
    if not cond:
        raise CheckFail(what)


def check(cond, what):
    if not cond:
        FAILS.append(what)
    return 1 if cond else 0


# ---------------------------------------------------------------- integer polynomials, low degree first

def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pdivmod_monic(a, b):
    a = ptrim(a)
    b = ptrim(b)
    need(b and b[-1] == 1, "divisor not monic")
    db = len(b) - 1
    r = list(a)
    if len(r) - 1 < db:
        return [], ptrim(r)
    q = [0] * (len(r) - db)
    for k in range(len(r) - 1, db - 1, -1):
        c = r[k]
        if c:
            q[k - db] = c
            for i in range(db + 1):
                r[k - db + i] -= c * b[i]
    return ptrim(q), ptrim(r)


CYC = {}


def cyclotomic(m):
    if m in CYC:
        return CYC[m]
    p = [-1] + [0] * (m - 1) + [1]
    for d in range(1, m):
        if m % d == 0:
            q, rem = pdivmod_monic(p, cyclotomic(d))
            need(not rem, "cyclotomic division left a remainder")
            p = q
    CYC[m] = p
    return p


def prime_factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n = n // p
        p += 1
    if n > 1:
        out.append(n)
    return out


# ---------------------------------------------------------------- the ring Z[x]/f, f monic; elements are integer tuples of length deg f

class Ring:
    def __init__(self, name, f, m):
        f = ptrim(f)
        need(f and f[-1] == 1, "ring modulus not monic")
        self.name = name
        self.f = f
        self.r = len(f) - 1
        self.m = m
        r = self.r
        self.zero = tuple([0] * r)
        self.one = tuple([1] + [0] * (r - 1))
        self.red = {}
        v = [-c for c in f[:r]]
        for k in range(r, 2 * r - 1):
            self.red[k] = tuple(v)
            top = v[-1]
            w = [0] + v[:-1]
            if top:
                for i in range(r):
                    w[i] -= top * f[i]
            v = w
        self.units = [k for k in range(1, m + 1) if gcd(k, m) == 1] if m else []

    def from_poly(self, c):
        q, rem = pdivmod_monic(c, self.f)
        rem = list(rem) + [0] * (self.r - len(rem))
        return tuple(rem)

    def from_int(self, n):
        return tuple([n] + [0] * (self.r - 1))

    def add(self, a, b):
        return tuple([x + y for x, y in zip(a, b)])

    def sub(self, a, b):
        return tuple([x - y for x, y in zip(a, b)])

    def neg(self, a):
        return tuple([-x for x in a])

    def smul(self, n, a):
        return tuple([n * x for x in a])

    def is_zero(self, a):
        return not any(a)

    def is_int(self, a):
        return not any(a[1:])

    def mul(self, a, b):
        r = self.r
        if r == 1:
            return (a[0] * b[0],)
        c = [0] * (2 * r - 1)
        for i in range(r):
            x = a[i]
            if x:
                for j in range(r):
                    y = b[j]
                    if y:
                        c[i + j] += x * y
        out = c[:r]
        for k in range(r, 2 * r - 1):
            t = c[k]
            if t:
                v = self.red[k]
                for i in range(r):
                    if v[i]:
                        out[i] += t * v[i]
        return tuple(out)

    def power(self, a, e):
        need(e >= 0, "negative exponent")
        out = self.one
        base = a
        while e:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e = e >> 1
        return out

    def sigma(self, a, k):
        m = self.m
        need(m > 0, "sigma outside a cyclotomic ring")
        v = [0] * m
        for i in range(len(a)):
            c = a[i]
            if c:
                v[(i * k) % m] += c
        return self.from_poly(v)

    def conj(self, a):
        if self.m <= 2:
            return a
        return self.sigma(a, self.m - 1)

    def norm_cofactor(self, b):
        need(not self.is_zero(b), "norm of zero")
        cof = self.one
        for k in self.units:
            if k != 1:
                cof = self.mul(cof, self.sigma(b, k))
        n = self.mul(b, cof)
        need(self.is_int(n) and n[0] != 0, "norm is not a nonzero integer in " + self.name)
        return cof, n[0]

    def div_by(self, a, cofn):
        cof, n = cofn
        t = self.mul(a, cof)
        out = []
        for c in t:
            q, rem = divmod(c, n)
            need(rem == 0, "inexact division in " + self.name)
            out.append(q)
        return tuple(out)

    def halve(self, a):
        out = []
        for c in a:
            q, rem = divmod(c, 2)
            need(rem == 0, "odd coefficient where an even one is required")
            out.append(q)
        return tuple(out)


RINGS = {}


def cyc_ring(m):
    if m not in RINGS:
        RINGS[m] = Ring("Z[zeta_" + str(m) + "]", cyclotomic(m), m)
    return RINGS[m]


def is_primitive_root(R, z, n):
    if R.power(z, n) != R.one:
        return 0
    for p in prime_factors(n):
        if R.power(z, n // p) == R.one:
            return 0
    return 1


# ---------------------------------------------------------------- 2x2 matrices (a, b, c, d) row major

def mat_mul(R, A, B):
    a, b, c, d = A
    e, f, g, h = B
    return (R.add(R.mul(a, e), R.mul(b, g)), R.add(R.mul(a, f), R.mul(b, h)),
            R.add(R.mul(c, e), R.mul(d, g)), R.add(R.mul(c, f), R.mul(d, h)))


def mat_det(R, A):
    a, b, c, d = A
    return R.sub(R.mul(a, d), R.mul(b, c))


def mat_dagger(R, A):
    a, b, c, d = A
    return (R.conj(a), R.conj(c), R.conj(b), R.conj(d))


def mat_scalar(R, s):
    return (s, R.zero, R.zero, s)


def closure(gens, mul, ident, cap):
    seen = {ident}
    order = [ident]
    frontier = [ident]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = mul(x, g)
                if y not in seen:
                    seen.add(y)
                    order.append(y)
                    nxt.append(y)
                    if len(order) > cap:
                        return order
        frontier = nxt
    return order


def group_checks(elems, gens, mul, ident, expected, pairwise):
    S = set(elems)
    distinct = 1 if len(S) == len(elems) else 0
    has_id = 1 if ident in S else 0
    gen_closed = 1
    for x in elems:
        for g in gens:
            if mul(x, g) not in S:
                gen_closed = 0
                break
        if not gen_closed:
            break
    pair_closed = 1
    if pairwise:
        for x in elems:
            for y in elems:
                if mul(x, y) not in S:
                    pair_closed = 0
                    break
            if not pair_closed:
                break
    order_ok = 1 if len(S) == expected else 0
    return {"distinct": distinct, "identity": has_id, "closed_gens": gen_closed, "closed_pairs": pair_closed, "order_ok": order_ok}


# ---------------------------------------------------------------- group builders

def build_cyclic(n):
    R = cyc_ring(n)
    z = R.from_poly([0, 1])
    prim = is_primitive_root(R, z, n)
    zi = R.power(z, n - 1)
    g = (z, R.zero, R.zero, zi)
    ident = mat_scalar(R, R.one)
    mul = lambda A, B: mat_mul(R, A, B)
    elems = closure([g], mul, ident, CLOSURE_CAP)
    return {"name": "C" + str(n), "family": "cyclic", "n": n, "expected": n, "ringA": R, "ringB": R, "elems": elems,
            "gens": [g], "mul": mul, "ident": ident, "kind": "matrix", "prim": prim, "gensB": [g], "scaleB": 1}


def build_dicyclic(n):
    R = cyc_ring(2 * n)
    w = R.from_poly([0, 1])
    prim = is_primitive_root(R, w, 2 * n)
    wi = R.power(w, 2 * n - 1)
    g1 = (w, R.zero, R.zero, wi)
    g2 = (R.zero, R.one, R.neg(R.one), R.zero)
    ident = mat_scalar(R, R.one)
    mul = lambda A, B: mat_mul(R, A, B)
    elems = closure([g1, g2], mul, ident, CLOSURE_CAP)
    return {"name": "Dic" + str(n), "family": "dicyclic", "n": n, "expected": 4 * n, "ringA": R, "ringB": R, "elems": elems,
            "gens": [g1, g2], "mul": mul, "ident": ident, "kind": "matrix", "prim": prim, "gensB": [g1, g2], "scaleB": 1}


def q_mul(R, P, Q):
    a1, b1, c1, d1 = P
    a2, b2, c2, d2 = Q
    m = R.mul
    a = R.sub(R.sub(R.sub(m(a1, a2), m(b1, b2)), m(c1, c2)), m(d1, d2))
    b = R.sub(R.add(R.add(m(a1, b2), m(b1, a2)), m(c1, d2)), m(d1, c2))
    c = R.add(R.add(R.sub(m(a1, c2), m(b1, d2)), m(c1, a2)), m(d1, b2))
    d = R.add(R.sub(R.add(m(a1, d2), m(b1, c2)), m(c1, b2)), m(d1, a2))
    return (R.halve(a), R.halve(b), R.halve(c), R.halve(d))


def embed_fn(R, Rb, theta):
    # Z[x]/f -> Rb, x -> theta; requires f(theta) = 0 in Rb
    acc = Rb.zero
    tp = Rb.one
    for c in R.f:
        acc = Rb.add(acc, Rb.smul(c, tp))
        tp = Rb.mul(tp, theta)
    hom = 1 if Rb.is_zero(acc) else 0
    pw = [Rb.one]
    for k in range(1, R.r):
        pw.append(Rb.mul(pw[-1], theta))

    def emb(a):
        out = Rb.zero
        for k in range(R.r):
            if a[k]:
                out = Rb.add(out, Rb.smul(a[k], pw[k]))
        return out
    return emb, hom


def build_quat(name, R, gens, expected, mb, theta_poly):
    ident = (R.from_int(2), R.zero, R.zero, R.zero)
    mul = lambda P, Q: q_mul(R, P, Q)
    elems = closure(gens, mul, ident, CLOSURE_CAP)
    Rb = cyc_ring(mb)
    theta = Rb.from_poly(theta_poly)
    emb, hom = embed_fn(R, Rb, theta)
    ib = Rb.from_poly([0] * (mb // 4) + [1])
    i_ok = 1 if Rb.mul(ib, ib) == Rb.neg(Rb.one) else 0

    def qmat(q):
        A, B, C, D = [emb(t) for t in q]
        return (Rb.add(A, Rb.mul(B, ib)), Rb.add(C, Rb.mul(D, ib)), Rb.add(Rb.neg(C), Rb.mul(D, ib)), Rb.sub(A, Rb.mul(B, ib)))
    return {"name": name, "family": "exceptional", "n": expected, "expected": expected, "ringA": R, "ringB": Rb, "elems": elems,
            "gens": gens, "mul": mul, "ident": ident, "kind": "quat", "prim": i_ok * hom, "qmat": qmat, "emb": emb,
            "gensB": [qmat(g) for g in gens], "scaleB": 2}


def element_checks(G):
    R = G["ringA"]
    elems = G["elems"]
    det_ok = 0
    unit_ok = 0
    extra_ok = 0
    if G["kind"] == "matrix":
        I = mat_scalar(R, R.one)
        for A in elems:
            if mat_det(R, A) == R.one:
                det_ok += 1
            if mat_mul(R, A, mat_dagger(R, A)) == I:
                unit_ok += 1
        extra_ok = len(elems)
        traces = [R.add(A[0], A[3]) for A in elems]
    else:
        Rb = G["ringB"]
        four = R.from_int(4)
        four_b = mat_scalar(Rb, Rb.from_int(4))
        qmat = G["qmat"]
        emb = G["emb"]
        mats = {}
        for q in elems:
            nq = R.zero
            for t in q:
                nq = R.add(nq, R.mul(t, t))
            M = qmat(q)
            mats[q] = M
            if nq == four and mat_det(Rb, M) == Rb.from_int(4):
                det_ok += 1
            if mat_mul(Rb, M, mat_dagger(Rb, M)) == four_b:
                unit_ok += 1
        # homomorphism and trace agreement: mat(x g) * 2 = mat(x) mat(g), trace(mat(x)) = 2 emb(A)
        for q in elems:
            good = 1
            if Rb.add(mats[q][0], mats[q][3]) != Rb.smul(2, emb(q[0])):
                good = 0
            for g in G["gens"]:
                lhs = tuple([Rb.smul(2, t) for t in qmat(q_mul(R, q, g))])
                if lhs != mat_mul(Rb, mats[q], qmat(g)):
                    good = 0
            extra_ok += good
        traces = [q[0] for q in elems]
    return det_ok, unit_ok, extra_ok, traces


# ---------------------------------------------------------------- method A

def trace_multiset(R, traces):
    counts = {}
    for t in traces:
        counts[t] = counts.get(t, 0) + 1
    return sorted(counts.items())


def char_sums(R, tm, kmax):
    sums = [R.zero] * (kmax + 1)
    for t, cnt in tm:
        prev = R.one
        cur = t
        sums[0] = R.add(sums[0], R.smul(cnt, prev))
        if kmax >= 1:
            sums[1] = R.add(sums[1], R.smul(cnt, cur))
        for k in range(1, kmax):
            nxt = R.sub(R.mul(t, cur), prev)
            prev = cur
            cur = nxt
            sums[k + 1] = R.add(sums[k + 1], R.smul(cnt, cur))
    return sums


def power_sums(R, tm, nmax):
    sums = [R.zero] * (nmax + 1)
    for t, cnt in tm:
        cur = R.one
        for n in range(nmax + 1):
            sums[n] = R.add(sums[n], R.smul(cnt, cur))
            cur = R.mul(cur, t)
    return sums


def averages(R, sums, order):
    # each sum must be an ordinary integer divisible by |G|; returns (values, bad_count)
    vals = []
    bad = 0
    for s in sums:
        if not R.is_int(s):
            bad += 1
            vals.append(None)
            continue
        q, rem = divmod(s[0], order)
        if rem:
            bad += 1
            vals.append(None)
            continue
        vals.append(q)
    return vals, bad


# ---------------------------------------------------------------- method B

def lmul(R, a, b):
    out = [R.zero] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        x = a[i]
        if R.is_zero(x):
            continue
        for j in range(len(b)):
            y = b[j]
            if R.is_zero(y):
                continue
            out[i + j] = R.add(out[i + j], R.mul(x, y))
    return out


def bareiss_rank(R, rows, ncols):
    rows = [dict(r) for r in rows if r]
    prev = None
    rank = 0
    for col in range(ncols):
        piv = -1
        for idx in range(len(rows)):
            if col in rows[idx]:
                piv = idx
                break
        if piv < 0:
            continue
        prow = rows.pop(piv)
        p = prow[col]
        cofn = R.norm_cofactor(prev) if prev is not None else None
        new_rows = []
        for row in rows:
            a = row.get(col)
            keys = set(row.keys())
            if a is not None:
                keys = keys | set(prow.keys())
            keys.discard(col)
            new = {}
            for j in sorted(keys):
                v = R.mul(p, row[j]) if j in row else R.zero
                if a is not None and j in prow:
                    v = R.sub(v, R.mul(a, prow[j]))
                if cofn is not None:
                    v = R.div_by(v, cofn)
                if not R.is_zero(v):
                    new[j] = v
            if new:
                new_rows.append(new)
        rows = new_rows
        prev = p
        rank += 1
    return rank


def invariant_dims(Rb, gens, scale, dmax):
    PQ = []
    for M in gens:
        a, b, c, e = M
        P = [[Rb.one]]
        Q = [[Rb.one]]
        for k in range(dmax):
            P.append(lmul(Rb, P[-1], [b, a]))
            Q.append(lmul(Rb, Q[-1], [e, c]))
        PQ.append((P, Q))
    dims = []
    sd = 1
    for d in range(dmax + 1):
        sde = Rb.from_int(sd)
        rows = []
        for P, Q in PQ:
            cols = [lmul(Rb, P[k], Q[d - k]) for k in range(d + 1)]
            for l in range(d + 1):
                row = {}
                for k in range(d + 1):
                    v = cols[k][l]
                    if l == k:
                        v = Rb.sub(v, sde)
                    if not Rb.is_zero(v):
                        row[k] = v
                if row:
                    rows.append(row)
        rank = bareiss_rank(Rb, rows, d + 1)
        dims.append(d + 1 - rank)
        sd = sd * scale
    return dims


# ---------------------------------------------------------------- rationals as reduced (num, den)

def rnorm(n, d):
    need(d != 0, "zero denominator")
    if d < 0:
        n = -n
        d = -d
    g = gcd(n, d)
    return (n // g, d // g)


def radd(a, b):
    return rnorm(a[0] * b[1] + b[0] * a[1], a[1] * b[1])


def rsub(a, b):
    return rnorm(a[0] * b[1] - b[0] * a[1], a[1] * b[1])


def rmul(a, b):
    return rnorm(a[0] * b[0], a[1] * b[1])


def rstr(a):
    return str(a[0]) + ":" + str(a[1])


def zseries(M, nmax):
    out = []
    fact = 1
    for n in range(nmax + 1):
        if n > 0:
            fact = fact * n
        out.append(rnorm(M[n], fact * (1 << n)))
    return out


def log_series(z, nmax):
    need(z[0] == (1, 1), "z(0) is not 1")
    f = [(0, 1)] * (nmax + 1)
    for n in range(1, nmax + 1):
        acc = (0, 1)
        for k in range(1, n):
            acc = radd(acc, rmul((k, 1), rmul(f[k], z[n - k])))
        f[n] = rsub(z[n], rmul(acc, (1, n)))
    return f


def exp_series(f, nmax):
    e = [(0, 1)] * (nmax + 1)
    e[0] = (1, 1)
    for n in range(1, nmax + 1):
        acc = (0, 1)
        for k in range(1, n + 1):
            acc = radd(acc, rmul((k, 1), rmul(f[k], e[n - k])))
        e[n] = rmul(acc, (1, n))
    return e


def catalan(kmax):
    c = [1]
    for k in range(kmax):
        s = 0
        for i in range(k + 1):
            s += c[i] * c[k - i]
        c.append(s)
    return c


# ---------------------------------------------------------------- output

ROWS = []


def emit(tag, kv):
    parts = [tag]
    for k, v in kv:
        parts.append(k + "=" + str(v))
    parts.append("json=0")
    ROWS.append("|".join(parts))


def ilist(v):
    return ",".join([str(x) for x in v])


def measure_group(G, haar, f_su2):
    name = G["name"]
    R = G["ringA"]
    elems = G["elems"]
    order = G["expected"]
    chk = group_checks(elems, G["gens"], G["mul"], G["ident"], order, True)
    det_ok, unit_ok, extra_ok, traces = element_checks(G)
    built = len(elems)
    ok = 1
    ok &= check(G["prim"] == 1, name + " generator root or embedding check")
    for key in ("distinct", "identity", "closed_gens", "closed_pairs", "order_ok"):
        ok &= check(chk[key] == 1, name + " " + key)
    ok &= check(det_ok == built, name + " det")
    ok &= check(unit_ok == built, name + " unitarity")
    ok &= check(extra_ok == built, name + " matrix form")
    tm = trace_multiset(R, traces)
    emit("GROUP", [("group", name), ("family", G["family"]), ("n", G["n"]), ("order_expected", order), ("elements_built", built),
                   ("ring_traces", R.name), ("ring_invariants", G["ringB"].name), ("generator_ok", G["prim"]),
                   ("distinct", chk["distinct"]), ("identity", chk["identity"]), ("closed_gens", chk["closed_gens"]),
                   ("closed_pairs", chk["closed_pairs"]), ("order_ok", chk["order_ok"]), ("det1", det_ok), ("unitary", unit_ok),
                   ("matrix_form_ok", extra_ok), ("distinct_traces", len(tm)), ("ok", ok), ("tag", "MEASURED")])
    # method A
    csum = char_sums(R, tm, TWOJ_MAX)
    m, bad_m = averages(R, csum, order)
    ok_m = check(bad_m == 0, name + " character sums integer and divisible")
    psum = power_sums(R, tm, MOMENT_MAX)
    M, bad_M = averages(R, psum, order)
    ok_M = check(bad_M == 0, name + " moment sums integer and divisible")
    if not (ok_m and ok_M):
        return None
    emit("MOLIEN", [("group", name), ("twoj", "0.." + str(TWOJ_MAX)), ("m", ilist(m)), ("method", "characters"),
                    ("sums_integer_divisible", len(m) - bad_m), ("tag", "MEASURED")])
    emit("MOMENT", [("group", name), ("n", "0.." + str(MOMENT_MAX)), ("M", ilist(M)), ("sums_integer_divisible", len(M) - bad_M),
                    ("tag", "MEASURED")])
    # method B
    dims = invariant_dims(G["ringB"], G["gensB"], G["scaleB"], INV_DMAX)
    agreeB = 1 if dims == m[:INV_DMAX + 1] else 0
    check(agreeB == 1, name + " method B equals method A")
    emit("INVARIANT", [("group", name), ("d", "0.." + str(INV_DMAX)), ("dim", ilist(dims)), ("method", "kernel_of_rho_minus_I_bareiss"),
                       ("generators", len(G["gensB"])), ("equals_characters", agreeB), ("tag", "MEASURED")])
    # design degree and first break
    fb = 0
    for k in range(1, TWOJ_MAX + 1):
        if m[k] != 0:
            fb = k
            break
    nonzero = [k for k in range(1, TWOJ_MAX + 1) if m[k] != 0]
    if fb:
        emit("DESIGN", [("group", name), ("t", fb - 1), ("first_break_twoj", fb), ("first_break_m", m[fb]),
                        ("nonzero_twoj_count", len(nonzero)), ("tag", "MEASURED")])
    else:
        emit("DESIGN", [("group", name), ("t", "ge" + str(TWOJ_MAX)), ("first_break_twoj", "none"), ("tag", "MEASURED")])
    # moment break against SU(2)
    mb = 0
    for n in range(MOMENT_MAX + 1):
        if M[n] != haar[n]:
            mb = n
            break
    if mb or M[0] != haar[0]:
        consistent = 1 if (mb == fb and M[mb] - haar[mb] == m[fb]) else 0
        check(consistent == 1, name + " moment break matches the first character break")
        emit("MOMENT_BREAK", [("group", name), ("n", mb), ("M", M[mb]), ("su2", haar[mb]), ("diff", M[mb] - haar[mb]),
                              ("equals_first_character_break", consistent), ("tag", "MEASURED")])
    else:
        emit("MOMENT_BREAK", [("group", name), ("n", "none_le" + str(MOMENT_MAX)), ("tag", "MEASURED")])
    # gauge
    z = zseries(M, BETA_MAX)
    f = log_series(z, BETA_MAX)
    check(exp_series(f, BETA_MAX) == z, name + " exp(log z) = z")
    diff = [rsub(f[n], f_su2[n]) for n in range(BETA_MAX + 1)]
    emit("GAUGE_F", [("group", name), ("n", "0.." + str(BETA_MAX)), ("f", ",".join([rstr(x) for x in f])), ("tag", "MEASURED")])
    emit("GAUGE_DIFF", [("group", name), ("n", "0.." + str(BETA_MAX)), ("diff", ",".join([rstr(x) for x in diff])), ("tag", "MEASURED")])
    gb = -1
    for n in range(BETA_MAX + 1):
        if diff[n][0] != 0:
            gb = n
            break
    if gb >= 0:
        expect = rnorm(m[fb], 1) if fb else (0, 1)
        if fb:
            fact = 1
            for k in range(2, fb + 1):
                fact = fact * k
            expect = rnorm(m[fb], fact * (1 << fb))
        consistent = 1 if (gb == fb and diff[gb] == expect) else 0
        check(consistent == 1, name + " gauge first break equals m(first break) / (b! 2^b)")
        emit("GAUGE_FIRST_BREAK", [("group", name), ("n", gb), ("diff", rstr(diff[gb])), ("equals_m_over_factorial", consistent),
                                   ("tag", "MEASURED")])
    else:
        emit("GAUGE_FIRST_BREAK", [("group", name), ("n", "none_le" + str(BETA_MAX)), ("tag", "MEASURED")])
    return {"m": m, "M": M, "dims": dims, "fb": fb, "mb": mb, "gb": gb, "diff": diff, "f": f, "tm": tm, "order": order}


# ---------------------------------------------------------------- controls, each built to fail

def control_rows(G2I, res2I):
    out_ok = 1
    # 1. drop one element of 2I
    elems = list(G2I["elems"])
    dropped = elems[-1]
    kept = elems[:-1]
    chk = group_checks(kept, G2I["gens"], G2I["mul"], G2I["ident"], 120, True)
    detected = 1 if (chk["closed_gens"] == 0 or chk["closed_pairs"] == 0 or chk["order_ok"] == 0) else 0
    out_ok &= detected
    emit("CONTROL", [("name", "drop_one_element_2I"), ("elements", len(kept)), ("closed_gens", chk["closed_gens"]),
                     ("closed_pairs", chk["closed_pairs"]), ("order_ok", chk["order_ok"]), ("must_fail", 1), ("detected", detected),
                     ("tag", "MEASURED")])
    # 2. perturb one trace of 2I by +1
    R = G2I["ringA"]
    tm = list(res2I["tm"])
    t0, c0 = tm[0]
    tm2 = [(t0, c0 - 1), (R.add(t0, R.one), 1)] + tm[1:]
    csum = char_sums(R, tm2, TWOJ_MAX)
    vals, bad = averages(R, csum, 120)
    detected = 1 if bad > 0 else 0
    out_ok &= detected
    emit("CONTROL", [("name", "perturb_one_trace_2I"), ("twoj_failing_divisibility", bad), ("of", TWOJ_MAX + 1), ("must_fail", 1),
                     ("detected", detected), ("tag", "MEASURED")])
    # 3. method B for 2I with only the first generator
    dims1 = invariant_dims(G2I["ringB"], G2I["gensB"][:1], G2I["scaleB"], INV_DMAX)
    mism = sum([1 for d in range(INV_DMAX + 1) if dims1[d] != res2I["m"][d]])
    detected = 1 if mism > 0 else 0
    out_ok &= detected
    emit("CONTROL", [("name", "invariants_one_generator_2I"), ("dim", ilist(dims1)), ("d_disagreeing", mism), ("must_fail", 1),
                     ("detected", detected), ("tag", "MEASURED")])
    # 4. a non-unitary matrix, det 2
    Rc = cyc_ring(4)
    A = (Rc.from_int(2), Rc.zero, Rc.zero, Rc.one)
    u = 1 if mat_mul(Rc, A, mat_dagger(Rc, A)) == mat_scalar(Rc, Rc.one) else 0
    dt = 1 if mat_det(Rc, A) == Rc.one else 0
    detected = 1 if (u == 0 and dt == 0) else 0
    out_ok &= detected
    emit("CONTROL", [("name", "nonunitary_diag_2_1"), ("unitary", u), ("det1", dt), ("must_fail", 1), ("detected", detected), ("tag", "MEASURED")])
    # 5. perturb a log coefficient of 2I, exp must no longer give z
    z = zseries(res2I["M"], BETA_MAX)
    f = list(res2I["f"])
    f[12] = radd(f[12], (1, 1000081))
    same = 1 if exp_series(f, BETA_MAX) == z else 0
    detected = 1 if same == 0 else 0
    out_ok &= detected
    emit("CONTROL", [("name", "perturb_log_coefficient_2I"), ("exp_equals_z", same), ("must_fail", 1), ("detected", detected), ("tag", "MEASURED")])
    return out_ok


# ---------------------------------------------------------------- comparison with the frozen prediction, read only after measurement

CMPSTATE = {"required": [0, 0, 0], "extra": [0, 0, 0]}
DISAGREE = []
MISSING = []


def pred_name(name):
    if name.startswith("Dic"):
        return "Dic_" + name[3:]
    if name.startswith("C"):
        return "C_" + name[1:]
    return name


def norm_tok(tok):
    tok = tok.strip()
    try:
        if ":" in tok:
            a, b = tok.split(":", 1)
            return rstr(rnorm(int(a), int(b)))
        return str(int(tok))
    except (ValueError, CheckFail):
        return tok


def cmp_values(frozen_items, measured_items):
    n = max(len(frozen_items), len(measured_items))
    bad = []
    for i in range(n):
        fv = norm_tok(frozen_items[i]) if i < len(frozen_items) else "absent"
        mv = norm_tok(measured_items[i]) if i < len(measured_items) else "absent"
        if fv != mv:
            bad.append((i, fv, mv))
    return n, bad


def tail(v):
    return v.rsplit("=", 1)[-1] if v is not None else None


def compare(cat, quantity, group, field, frozen, measured, is_list):
    if frozen is None:
        MISSING.append(quantity + "/" + group + "/" + field)
        emit("COMPARE_MISSING", [("set", cat), ("quantity", quantity), ("group", group), ("field", field), ("missing", 1), ("ok", 0),
                                 ("tag", "MEASURED")])
        return 0
    fi = frozen.split(",") if is_list else [frozen]
    mi = [str(x) for x in measured] if is_list else [str(measured)]
    n, bad = cmp_values(fi, mi)
    agree = 1 if not bad else 0
    st = CMPSTATE[cat]
    st[0] += 1
    st[1] += n
    st[2] += n - len(bad)
    emit("COMPARE", [("set", cat), ("quantity", quantity), ("group", group), ("field", field), ("frozen", frozen),
                     ("measured", ",".join(mi)), ("values", n), ("agree_values", n - len(bad)), ("agree", agree), ("tag", "MEASURED")])
    for i, fv, mv in bad:
        DISAGREE.append((cat, quantity, group, field, i, fv, mv))
        emit("COMPARE_DISAGREE", [("set", cat), ("quantity", quantity), ("group", group), ("field", field), ("index", i), ("frozen", fv),
                                  ("measured", mv), ("tag", "MEASURED")])
    return agree


def expand_molien(H, dmax):
    mt = re.fullmatch(r"num\(([^()]*)\)_over_den\(((?:\(1-t\^[0-9]+\))*)\)", H)
    if mt is None:
        return None
    ser = [0] * (dmax + 1)
    for term in mt.group(1).split("+"):
        if term == "1":
            e = 0
        elif term.startswith("t^") and term[2:].isdigit():
            e = int(term[2:])
        else:
            return None
        if e <= dmax:
            ser[e] += 1
    for k in [int(x) for x in re.findall(r"\(1-t\^([0-9]+)\)", mt.group(2))]:
        if k <= 0:
            return None
        for d in range(k, dmax + 1):
            ser[d] += ser[d - k]
    return ser


def parse_hbp(path):
    with open(path, "rb") as fh:
        data = fh.read()
    rows = []
    for line in data.decode("utf-8").split("\n"):
        if not line:
            continue
        parts = line.split("|")
        kv = {}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                kv[k] = v
        rows.append((parts[0], kv))
    return rows, data


def factorial(n):
    out = 1
    for k in range(2, n + 1):
        out = out * k
    return out


def compare_stage(results, haar, f_su2, pred_path):
    present = 1 if os.path.isfile(pred_path) else 0
    side = pred_path + ".sha256"
    if not present:
        check(False, "prediction file missing")
        emit("PREDICTION_FILE", [("path", pred_path), ("missing", 1), ("ok", 0), ("tag", "MEASURED")])
        return 0
    rows, data = parse_hbp(pred_path)
    digest = hashlib.sha256(data).hexdigest()
    side_hex = ""
    if os.path.isfile(side):
        with open(side, "rb") as fh:
            side_hex = fh.read().decode("utf-8").split(" ")[0].strip()
    side_ok = check(side_hex == digest, "prediction sidecar sha256 mismatch")
    frozen_ok = check(digest.startswith(PREDICTION_SHA_PREFIX), "prediction sha256 is not the frozen c2601cf0 file")
    emit("PREDICTION_FILE", [("path", pred_path), ("rows", len(rows)), ("sha256", digest), ("sidecar_equal", side_ok),
                             ("frozen_prefix_c2601cf0", frozen_ok), ("ok", side_ok & frozen_ok), ("tag", "MEASURED")])
    by = {}
    for tag, kv in rows:
        by.setdefault(tag, {})
        if "group" in kv:
            by[tag][kv["group"]] = kv
        else:
            by[tag]["_"] = kv
    names = [k for k in results]
    for name in names:
        r = results[name]
        pn = pred_name(name)
        m = r["m"]
        M = r["M"]
        fb = r["fb"]
        mb = r["mb"]
        gb = r["gb"]
        G = r["G"]
        mol = by.get("MOLIEN", {}).get(pn)
        compare("required", "MOLIEN", pn, "a_d0..120", mol.get("a") if mol else None, m, True)
        des = by.get("DESIGN", {}).get(pn)
        compare("required", "DESIGN", pn, "order", des.get("order") if des else None, len(G["elems"]), False)
        compare("required", "DESIGN", pn, "t_pred", des.get("t_pred") if des else None, fb - 1, False)
        compare("required", "DESIGN", pn, "first_break_degree", des.get("first_break_degree") if des else None, fb, False)
        compare("required", "DESIGN", pn, "m_at_break", des.get("m_at_break") if des else None, m[fb], False)
        mom = by.get("MOMENT", {}).get(pn)
        compare("required", "MOMENT", pn, "M_n0..60", mom.get("M") if mom else None, M, True)
        mbr = by.get("MOMENT_BREAK", {}).get(pn)
        compare("required", "MOMENT_BREAK", pn, "first_n_with_M_ne_Mh", mbr.get("first_n_with_M_ne_Mh") if mbr else None, mb, False)
        compare("required", "MOMENT_BREAK", pn, "M", mbr.get("M") if mbr else None, M[mb], False)
        compare("required", "MOMENT_BREAK", pn, "Mh", mbr.get("Mh") if mbr else None, haar[mb], False)
        compare("required", "MOMENT_BREAK", pn, "diff", mbr.get("diff") if mbr else None, M[mb] - haar[mb], False)
        compare("required", "MOMENT_BREAK", pn, "predicted_diff_value", tail(mbr.get("predicted_diff")) if mbr else None, m[fb], False)
        gd = by.get("GAUGE_DIFF", {}).get(pn)
        compare("required", "GAUGE_DIFF", pn, "diff_beta0..30", gd.get("diff") if gd else None, [rstr(x) for x in r["diff"]], True)
        gfb = by.get("GAUGE_FIRST_BREAK", {}).get(pn)
        compare("required", "GAUGE_FIRST_BREAK", pn, "zero_through_beta_power", gfb.get("zero_through_beta_power") if gfb else None, gb - 1, False)
        compare("required", "GAUGE_FIRST_BREAK", pn, "first_nonzero_power", gfb.get("first_nonzero_power") if gfb else None, gb, False)
        compare("required", "GAUGE_FIRST_BREAK", pn, "value", gfb.get("value") if gfb else None, rstr(r["diff"][gb]), False)
        compare("required", "GAUGE_FIRST_BREAK", pn, "predicted_value", tail(gfb.get("predicted")) if gfb else None,
                rstr(rnorm(m[fb], factorial(fb) * (1 << fb))), False)
        # extra quantities the frozen rows also state
        R = G["ringA"]
        if G["kind"] == "matrix":
            negid = (R.neg(R.one), R.zero, R.zero, R.neg(R.one))
        else:
            negid = (R.from_int(-2), R.zero, R.zero, R.zero)
        has_neg = 1 if negid in set(G["elems"]) else 0
        compare("extra", "MOLIEN", pn, "contains_minus_one", mol.get("contains_minus_one") if mol else None, has_neg, False)
        if mol and "H" in mol:
            ser = expand_molien(mol["H"], TWOJ_MAX)
            if ser is None:
                compare("extra", "MOLIEN", pn, "H_series_vs_measured_m", "unparsed:" + mol["H"], m, True)
            else:
                compare("extra", "MOLIEN", pn, "H_series_vs_measured_m", ",".join([str(x) for x in ser]), m, True)
        gf = by.get("GAUGE_F", {}).get(pn)
        compare("extra", "GAUGE_F", pn, "f_beta0..30", gf.get("f") if gf else None, [rstr(x) for x in r["f"]], True)
    mh = by.get("MOMENT_HAAR", {}).get("SU2_Haar")
    compare("extra", "MOMENT_HAAR", "SU2_Haar", "Mh_n0..60", mh.get("Mh") if mh else None, haar, True)
    gh = by.get("GAUGE_HAAR", {}).get("SU2_Haar")
    compare("extra", "GAUGE_HAAR", "SU2_Haar", "f_beta0..30", gh.get("f") if gh else None, [rstr(x) for x in f_su2], True)
    # ceiling, measured over the 61 groups
    ts = [(results[k]["fb"] - 1, pred_name(k)) for k in names]
    tmax = max([t for t, k in ts])
    att = [k for t, k in ts if t == tmax]
    lower = [t for t, k in ts if t < tmax]
    t2 = max(lower) if lower else -1
    att2 = [k for t, k in ts if t == t2]
    ce = by.get("CEILING", {}).get("_")
    compare("required", "CEILING", "all", "max_t_pred", ce.get("max_t_pred") if ce else None, tmax, False)
    compare("required", "CEILING", "all", "attained_by", ce.get("attained_by") if ce else None, ",".join(att), False)
    compare("required", "CEILING", "all", "count_attaining", ce.get("count_attaining") if ce else None, len(att), False)
    compare("required", "CEILING", "all", "next_highest", ce.get("next_highest") if ce else None, str(t2) + "(" + ",".join(att2) + ")", False)
    compare("required", "CEILING", "all", "groups_listed", ce.get("groups_listed") if ce else None, len(names), False)
    top = [results[k]["fb"] for k in names if results[k]["fb"] - 1 == tmax]
    compare("required", "CEILING", "all", "first_break_at_ceiling", ce.get("first_break_at_ceiling") if ce else None,
            "degree_" + ",".join([str(x) for x in top]), False)
    cf = by.get("CEILING_FAMILY", {}).get("_")
    tc = sorted(set([results[k]["fb"] - 1 for k in names if k.startswith("C")]))
    td = sorted(set([results[k]["fb"] - 1 for k in names if k.startswith("Dic")]))
    compare("extra", "CEILING_FAMILY", "C_2..C_30", "t_values", cf.get("C_n_t_pred_values_n2..30") if cf else None, ",".join([str(x) for x in tc]), False)
    compare("extra", "CEILING_FAMILY", "Dic_2..Dic_30", "t_values", cf.get("Dic_n_t_pred_values_n2..30") if cf else None,
            ",".join([str(x) for x in td]), False)
    compare("extra", "CEILING_FAMILY", "Dic_2", "t_pred", tail(cf.get("Dic_2")) if cf else None, results["Dic2"]["fb"] - 1, False)
    flat = 1 if (len(tc) == 1 and len(td) == 1) else 0
    compare("extra", "CEILING_FAMILY", "C_n_and_Dic_n_n2..30", "growing_order_does_not_raise_t", cf.get("growing_order_does_not_raise_t") if cf else None,
            flat, False)
    emit("SCOPE", [("measured_n", "2..30"), ("ceiling_family_however_large_n", "NAMED_not_measured_here"),
                   ("classification_of_finite_subgroups", "NAMED_Klein_1884_not_measured_here"), ("tag", "NAMED")])
    # predicted groups this run did not measure
    measured_pn = set([pred_name(k) for k in names])
    for tag in ("MOLIEN", "DESIGN", "MOMENT", "MOMENT_BREAK", "GAUGE_DIFF", "GAUGE_FIRST_BREAK"):
        for g in by.get(tag, {}):
            if g not in measured_pn:
                MISSING.append(tag + "/" + g + "/not_measured")
                emit("COMPARE_MISSING", [("set", "required"), ("quantity", tag), ("group", g), ("measured", "absent"), ("missing", 1), ("ok", 0),
                                         ("tag", "MEASURED")])
    return 1


def control_compare(results):
    r = results["2I"]
    t = r["fb"] - 1
    n1, bad1 = cmp_values([str(t + 1)], [str(t)])
    fr = [str(x) for x in r["m"]]
    fr[12] = str(r["m"][12] + 1)
    n2, bad2 = cmp_values(fr, [str(x) for x in r["m"]])
    n3, bad3 = cmp_values(fr[:-1], [str(x) for x in r["m"]])
    fq = [rstr(x) for x in r["diff"]]
    fq[12] = rstr(rnorm(2, 1961990553600))
    n4, bad4 = cmp_values(fq, [rstr(x) for x in r["diff"]])
    detected = 1 if (len(bad1) == 1 and len(bad2) == 1 and bad2[0][0] == 12 and len(bad3) == 2 and len(bad4) == 1) else 0
    emit("CONTROL", [("name", "compare_perturbed_frozen_values_2I"), ("t_plus_1_disagreements", len(bad1)), ("a12_plus_1_disagreements", len(bad2)),
                     ("a12_plus_1_and_truncated_disagreements", len(bad3)), ("gauge12_doubled_disagreements", len(bad4)), ("must_fail", 1),
                     ("detected", detected), ("tag", "MEASURED")])
    return detected


# ---------------------------------------------------------------- main

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = sys.argv[1] if len(sys.argv) > 1 else here
    pred_path = sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/tmp/opus55/iconic-3-space/repo/predictions/su2-ladder/PREDICTION.hbp"
    cat = catalan(MOMENT_MAX)
    haar = [cat[n >> 1] if n % 2 == 0 else 0 for n in range(MOMENT_MAX + 1)]
    f_su2 = log_series(zseries(haar, BETA_MAX), BETA_MAX)
    check(exp_series(f_su2, BETA_MAX) == zseries(haar, BETA_MAX), "SU(2) exp(log z) = z")
    emit("RUN", [("script", "measure_ladder.py"), ("groups", 61), ("twoj_max", TWOJ_MAX), ("moment_max", MOMENT_MAX),
                 ("invariant_dmax", INV_DMAX), ("beta_max", BETA_MAX), ("characters", "once_per_distinct_trace_weighted_by_count"),
                 ("integer_only", 1), ("network", 0), ("tag", "MEASURED")])
    emit("SU2_HAAR", [("n", "0.." + str(MOMENT_MAX)), ("M", ilist(haar)), ("source", "catalan_recurrence"), ("tag", "MEASURED")])
    emit("SU2_GAUGE_F", [("n", "0.." + str(BETA_MAX)), ("f", ",".join([rstr(x) for x in f_su2])), ("tag", "MEASURED")])
    groups = []
    for n in range(2, NMAX + 1):
        groups.append(("C", n))
    for n in range(2, NMAX + 1):
        groups.append(("Dic", n))
    groups.append(("2T", 0))
    groups.append(("2O", 0))
    groups.append(("2I", 0))
    results = {}
    G2I = None
    for fam, n in groups:
        if fam == "C":
            G = build_cyclic(n)
        elif fam == "Dic":
            G = build_dicyclic(n)
        elif fam == "2T":
            Rz = Ring("Z", [0, 1], 0)
            gens = [tuple([Rz.from_int(1)] * 4), (Rz.zero, Rz.from_int(2), Rz.zero, Rz.zero)]
            G = build_quat("2T", Rz, gens, 24, 4, [0])
        elif fam == "2O":
            Rs = Ring("Z[sqrt2]", [-2, 0, 1], 0)
            s2 = Rs.from_poly([0, 1])
            gens = [tuple([Rs.one] * 4), (s2, s2, Rs.zero, Rs.zero)]
            G = build_quat("2O", Rs, gens, 48, 8, [0, 1, 0, 0, 0, 0, 0, 1])
        else:
            Rp = Ring("Z[phi]", [-1, -1, 1], 0)
            ph = Rp.from_poly([0, 1])
            phinv = Rp.sub(ph, Rp.one)
            gens = [tuple([Rp.one] * 4), (ph, Rp.one, phinv, Rp.zero)]
            theta = [1] + [0] * 3 + [1] + [0] * 11 + [1]
            G = build_quat("2I", Rp, gens, 120, 20, theta)
            G2I = G
        res = measure_group(G, haar, f_su2)
        if res is None:
            break
        res["G"] = G
        results[G["name"]] = res
    controls_ok = 0
    if G2I is not None and "2I" in results:
        controls_ok = control_rows(G2I, results["2I"])
    check(controls_ok == 1, "a control did not fail")
    measured_digest = hashlib.sha256(("\n".join(ROWS) + "\n").encode("utf-8")).hexdigest()
    emit("MEASURED_BLOCK", [("rows", len(ROWS)), ("sha256", measured_digest), ("written_before_prediction_read", 1), ("tag", "MEASURED")])
    same_pre = check(measured_digest == PRECOMPARE_MEASURED_BLOCK, "measured block differs from the run made before the prediction was read")
    emit("PRECOMPARE_RUN", [("measured_block_sha256", PRECOMPARE_MEASURED_BLOCK), ("file_sha256", PRECOMPARE_FILE),
                            ("measured_block_equal", same_pre), ("tag", "MEASURED")])
    # comparison: the prediction file is opened here, after every measured row above is fixed
    cmp_ok = 0
    if len(results) == 61:
        cmp_ok = compare_stage(results, haar, f_su2, pred_path)
        controls_ok = controls_ok & control_compare(results)
        check(controls_ok == 1, "the comparison control did not fail")
    else:
        check(False, "not every group was measured, comparison skipped")
    for cat in ("required", "extra"):
        st = CMPSTATE[cat]
        nd = len([d for d in DISAGREE if d[0] == cat])
        emit("COMPARE_SUMMARY", [("set", cat), ("rows_compared", st[0]), ("values_compared", st[1]), ("values_agree", st[2]),
                                 ("values_disagree", st[1] - st[2]), ("disagreement_rows", nd), ("tag", "MEASURED")])
    req_tags = "MOLIEN,DESIGN,MOMENT,MOMENT_BREAK,GAUGE_DIFF,GAUGE_FIRST_BREAK,CEILING"
    emit("COMPARE_SCOPE", [("required_tags", req_tags), ("extra_tags", "MOLIEN.contains_minus_one,MOLIEN.H,GAUGE_F,MOMENT_HAAR,GAUGE_HAAR,CEILING_FAMILY"),
                           ("missing", len(MISSING)), ("compare_ran", cmp_ok), ("tag", "MEASURED")])
    emit("RESULT", [("groups_measured", len(results)), ("internal_failures", len(FAILS)), ("controls_detected_all", controls_ok),
                    ("values_disagree", len(DISAGREE)), ("missing", len(MISSING)),
                    ("status", "PASS" if not FAILS else "FAIL"), ("tag", "MEASURED")])
    for w in FAILS:
        emit("FAIL", [("what", w.replace("|", "/").replace(" ", "_"))])
    text = "\n".join(ROWS) + "\n"
    data = text.encode("utf-8")
    hbp = os.path.join(out_dir, "MEASURE-LADDER.hbp")
    with open(hbp, "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(hbp + ".sha256", "wb") as fh:
        fh.write((digest + "  MEASURE-LADDER.hbp\n").encode("utf-8"))
    for row in ROWS:
        if row.startswith(("RESULT", "CONTROL", "FAIL", "MEASURED_BLOCK", "PRECOMPARE", "PREDICTION_FILE", "COMPARE_SUMMARY", "COMPARE_SCOPE",
                           "COMPARE_DISAGREE", "COMPARE_MISSING")):
            print(row)
    print("SHA256|file=MEASURE-LADDER.hbp|sha256=" + digest + "|json=0")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    try:
        code = main()
    except CheckFail as exc:
        print("CHECK_FAIL|what=" + str(exc).replace(" ", "_") + "|json=0")
        code = 1
    sys.exit(code)
