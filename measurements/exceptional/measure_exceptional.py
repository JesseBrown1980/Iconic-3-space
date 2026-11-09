#!/usr/bin/env python3
"""
measure_exceptional.py, the measuring side for predictions/exceptional/PREDICTION-EXCEPTIONAL.hbp of
JesseBrown1980/Iconic-3-space (main 8f355db). Written without opening predict_*.py. The frozen rows are
read only by the COMPARE stage at the end, after every measured row is computed.

Integers only. Every set is built from its definition:
  1. the 120 unit icosians (600-cell, the group 2I), 4-vectors over Z[phi] scaled by 2, phi^2 = phi + 1,
     phi^-1 = phi - 1; quaternion product (1, i, j, k) with Hamilton's rules;
  2. the 240 E8 roots scaled by 2 in Z^8;
  3. the extended binary Golay code from the cyclic generator of length 23
     g(x) = 1 + x^2 + x^4 + x^5 + x^6 + x^10 + x^11 (rows x^i g(x), i = 0..11, plus an overall parity bit),
     and from it the 196,560 Leech minimal vectors in Conway's construction scaled by sqrt8:
     (+-2^8, 0^16) on the octads with an even number of minus signs, (-+3, +-1^23) with signs on Golay
     codewords, and (+-4, +-4, 0^22).
Measured: counts, norms, closure, every inner product, distance distributions from every point, exact zonal
sums S_k = (1/N) sum_{x,y} Q_k(<x,y>/|x|^2) for k = 0..30 with Q_k built from the Gegenbauer recurrence and
cross-checked by the normalised three-term recurrence, moment identities for 2m = 2..14, t and the nonzero
profile. Controls that must fail are run through the same check functions.

Run: python3 -I -B measure_exceptional.py [--leech-stride K]   (K = 1, every Leech base point, is the record run)
Exit 1 when a check fails, a control fails to fail, the frozen file does not verify against its sidecar and expected sha256
prefix, or the measured block differs from MEASURE-EXCEPTIONAL.precompare.hbp (written before the frozen rows were opened).
Exit 2 when every check passes and some COMPARE row disagrees: a disagreement is a result, reported, never edited away. Exit 0 otherwise.
"""
import sys
import os
import hashlib
import itertools
import operator
import functools
from collections import Counter
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "MEASURE-EXCEPTIONAL.hbp"
KMAX = 30
MAXDEG = 14
ROWS = []
FAILED = []
NCHECK = [0]


def emit(tag, *kv):
    parts = [tag]
    i = 0
    while i < len(kv):
        parts.append(str(kv[i]) + "=" + str(kv[i + 1]))
        i += 2
    parts.append("json=0")
    ROWS.append("|".join(parts))


def check(name, ok, *kv):
    NCHECK[0] += 1
    if not ok:
        FAILED.append(name)
    emit("CHECK", "name", name, "ok", 1 if ok else 0, *kv)
    return ok


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


# ---------------------------------------------------------------- integers, rationals
def ipow(b, e):
    if e < 0:
        raise ValueError("negative exponent")
    r = 1
    while e:
        if e & 1:
            r = r * b
        b = b * b
        e >>= 1
    return r


def binom(a, b):
    if b < 0 or b > a:
        return 0
    r = 1
    for i in range(b):
        q, rem = divmod(r * (a - i), i + 1)
        if rem:
            raise ArithmeticError("binomial step not exact")
        r = q
    return r


def odd_double_factorial(m):
    # 1 * 3 * ... * (2m - 1)
    r = 1
    for i in range(m):
        r = r * (2 * i + 1)
    return r


def rstr(n, d):
    return str(n) + "/" + str(d)


# ---------------------------------------------------------------- Z[phi], phi^2 = phi + 1
ZP0 = (0, 0)
ZP1 = (1, 0)
PHI = (0, 1)
PHIINV = (-1, 1)


def zp_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zp_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def zp_mul(x, y):
    bd = x[1] * y[1]
    return (x[0] * y[0] + bd, x[0] * y[1] + x[1] * y[0] + bd)


def zp_scale(x, k):
    return (x[0] * k, x[1] * k)


def zp_sign(x):
    # sign of a + b*phi = (p + q*sqrt5)/2 with p = 2a + b, q = b
    p = 2 * x[0] + x[1]
    q = x[1]
    if p == 0 and q == 0:
        return 0
    if p >= 0 and q >= 0:
        return 1
    if p <= 0 and q <= 0:
        return -1
    if p > 0:
        return 1 if p * p > 5 * q * q else -1
    return 1 if 5 * q * q > p * p else -1


def zp_str(x):
    if x[1] == 0:
        return str(x[0])
    if x[1] < 0:
        return str(x[0]) + "-" + str(-x[1]) + "*phi"
    return str(x[0]) + "+" + str(x[1]) + "*phi"


def zp_cmp(x, y):
    return zp_sign(zp_sub(x, y))


# ---------------------------------------------------------------- Q(phi) as (A, B, d) = (A + B*phi)/d
def qp_norm(a, b, d):
    if d == 0:
        raise ZeroDivisionError("zero denominator")
    if d < 0:
        a, b, d = -a, -b, -d
    g = gcd(gcd(a, b), d)
    if g > 1:
        a, b, d = a // g, b // g, d // g
    return (a, b, d)


def qp_add(x, y):
    return qp_norm(x[0] * y[2] + y[0] * x[2], x[1] * y[2] + y[1] * x[2], x[2] * y[2])


def qp_mul(x, y):
    a, b = zp_mul((x[0], x[1]), (y[0], y[1]))
    return qp_norm(a, b, x[2] * y[2])


def qp_rat(x, p, q):
    return qp_norm(x[0] * p, x[1] * p, x[2] * q)


def qp_str(x):
    if x[1] == 0:
        return rstr(x[0], x[2])
    return "(" + zp_str((x[0], x[1])) + ")/" + str(x[2])


def qp_is_zero(x):
    return x[0] == 0 and x[1] == 0


# ---------------------------------------------------------------- rational polynomials (coeffs, den)
def pnorm(c, d):
    if d == 0:
        raise ZeroDivisionError("zero denominator")
    if d < 0:
        c = [-a for a in c]
        d = -d
    g = d
    for a in c:
        g = gcd(g, a)
    if g > 1:
        c = [a // g for a in c]
        d = d // g
    while len(c) > 1 and c[-1] == 0:
        c = c[:-1]
    return (tuple(c), d)


def padd(p, q):
    n = max(len(p[0]), len(q[0]))
    c = [0] * n
    for j in range(len(p[0])):
        c[j] += p[0][j] * q[1]
    for j in range(len(q[0])):
        c[j] += q[0][j] * p[1]
    return pnorm(c, p[1] * q[1])


def pscale(p, num, den):
    return pnorm([a * num for a in p[0]], p[1] * den)


def pshift(p):
    return pnorm([0] + list(p[0]), p[1])


def zonal_polys(n, kmax):
    # Gegenbauer C_k^(lam), lam = (n-2)/2: k C_k = 2(k+lam-1) t C_{k-1} - (k+2lam-2) C_{k-2},
    # i.e. k C_k = (2k+n-4) t C_{k-1} - (k+n-4) C_{k-2}, C_0 = 1, C_1 = 2 lam t = (n-2) t. Q_k = C_k / C_k(1).
    C = [pnorm([1], 1), pnorm([0, n - 2], 1)]
    for k in range(2, kmax + 1):
        ac, ad = C[k - 1]
        bc, bd = C[k - 2]
        num = [0] * (k + 1)
        for j in range(len(ac)):
            num[j + 1] += (2 * k + n - 4) * ac[j] * bd
        for j in range(len(bc)):
            num[j] -= (k + n - 4) * bc[j] * ad
        C.append(pnorm(num, ad * bd * k))
    Q = []
    at_one_ok = 0
    for k in range(kmax + 1):
        c, d = C[k]
        val1 = sum(c)
        if val1 == binom(k + n - 3, k) * d:
            at_one_ok += 1
        Q.append(pnorm(list(c), val1))
    # cross-check: (k+n-2) Q_{k+1} = (2k+n-2) t Q_k - k Q_{k-1}
    rec_ok = 0
    for k in range(1, kmax):
        lhs = pscale(Q[k + 1], k + n - 2, 1)
        rhs = padd(pscale(pshift(Q[k]), 2 * k + n - 2, 1), pscale(Q[k - 1], -k, 1))
        if lhs == rhs:
            rec_ok += 1
    q_at_one = sum(1 for q in Q if sum(q[0]) == q[1])
    return Q, at_one_ok, rec_ok, q_at_one


def check_zonal_polys(n):
    Q, at_one_ok, rec_ok, q_at_one = zonal_polys(n, KMAX)
    check("zonal.n" + str(n) + ".C_k(1)_eq_binom(k+n-3,k)", at_one_ok == KMAX + 1, "k", "0.." + str(KMAX), "ok_count", at_one_ok)
    check("zonal.n" + str(n) + ".normalised_recurrence", rec_ok == KMAX - 1, "k", "1.." + str(KMAX - 1), "ok_count", rec_ok)
    check("zonal.n" + str(n) + ".Q_k(1)_eq_1", q_at_one == KMAX + 1, "ok_count", q_at_one)
    check("zonal.n" + str(n) + ".Q_1_eq_t", Q[1] == ((0, 1), 1))
    return Q


def zonal_sums(tally, s, Q, R):
    # S_k = (1/R) sum_v tally[v] Q_k(v/s); v in Z[phi]; Q_k = sum_j c_j t^j / D
    out = []
    for k in range(len(Q)):
        qc, qd = Q[k]
        spow = [ipow(s, k - j) for j in range(len(qc))]
        A = 0
        B = 0
        for v, cnt in tally.items():
            vp = ZP1
            a0 = 0
            a1 = 0
            for j in range(len(qc)):
                c = qc[j] * spow[j]
                a0 += c * vp[0]
                a1 += c * vp[1]
                vp = zp_mul(vp, v)
            A += cnt * a0
            B += cnt * a1
        out.append(qp_norm(A, B, R * qd * ipow(s, k)))
    return out


def zonal_sums_recurrence(tally, s, n, kmax, R):
    # independent path: Q_k evaluated at each cosine by the normalised three-term recurrence
    totals = [(0, 0, 1)] * (kmax + 1)
    for v, cnt in tally.items():
        c = qp_norm(v[0], v[1], s)
        vals = [(1, 0, 1), c]
        for k in range(1, kmax):
            t1 = qp_rat(qp_mul(c, vals[k]), 2 * k + n - 2, 1)
            t2 = qp_rat(vals[k - 1], -k, 1)
            vals.append(qp_rat(qp_add(t1, t2), 1, k + n - 2))
        for k in range(kmax + 1):
            totals[k] = qp_add(totals[k], qp_rat(vals[k], cnt, 1))
    return [qp_rat(t, 1, R) for t in totals]


def design_report(name, n, s, N, tally, R, Q):
    S = zonal_sums(tally, s, Q, R)
    S2 = zonal_sums_recurrence(tally, s, n, KMAX, R)
    check(name + ".zonal_two_paths_agree", S == S2, "k", "0.." + str(KMAX))
    check(name + ".S_0_eq_N", S[0] == (N, 0, 1), "S_0", qp_str(S[0]))
    neg = [k for k in range(KMAX + 1) if zp_sign((S[k][0], S[k][1])) < 0]
    check(name + ".S_k_ge_0", len(neg) == 0, "negative_k", ",".join(str(k) for k in neg) if neg else "none")
    irr = [k for k in range(KMAX + 1) if S[k][1] != 0]
    emit("SK_FIELD", "set", name, "phi_part_nonzero_k", ",".join(str(k) for k in irr) if irr else "none", "tag", "MEASURED")
    for k in range(KMAX + 1):
        emit("SK", "set", name, "n", n, "k", k, "S", qp_str(S[k]), "zero", 1 if qp_is_zero(S[k]) else 0, "tag", "MEASURED")
    t = 0
    while t + 1 <= KMAX and qp_is_zero(S[t + 1]):
        t += 1
    nonzero = [k for k in range(1, KMAX + 1) if not qp_is_zero(S[k])]
    first = nonzero[0] if nonzero else "none"
    emit("T", "set", name, "t", t, "first_nonzero_k", first, "through", KMAX, "tag", "MEASURED")
    emit("PROFILE", "set", name, "nonzero_k", ",".join(str(k) for k in nonzero) if nonzero else "none",
         "zero_even_k", ",".join(str(k) for k in range(2, KMAX + 1, 2) if qp_is_zero(S[k])),
         "odd_all_zero", 1 if all(qp_is_zero(S[k]) for k in range(1, KMAX + 1, 2)) else 0, "through", KMAX, "tag", "MEASURED")
    return S, t, nonzero


def moment_report(name, n, s, N, ipcounters, labels):
    # identity: sum_x <x,v>^(2m) = N s^m |v|^(2m) (2m-1)!! / (n (n+2) ... (n+2m-2)); odd moments vanish
    first_fail = None
    per_v = []
    for vi in range(len(ipcounters)):
        cnt, vv = ipcounters[vi]
        sums = [ZP0] * (MAXDEG + 1)
        for v, c in cnt.items():
            p = ZP1
            for d in range(MAXDEG + 1):
                sums[d] = (sums[d][0] + c * p[0], sums[d][1] + c * p[1])
                p = zp_mul(p, v)
        odd_ok = all(sums[d] == ZP0 for d in range(1, MAXDEG + 1, 2))
        check(name + ".moment.v" + str(vi) + ".odd_vanish", odd_ok, "v", labels[vi])
        vfirst = None
        for m in range(1, MAXDEG // 2 + 1):
            d = 2 * m
            num = N * ipow(s, m) * ipow(vv, m) * odd_double_factorial(m)
            den = 1
            for i in range(m):
                den = den * (n + 2 * i)
            ok = sums[d][1] == 0 and sums[d][0] * den == num
            g = gcd(num, den)
            emit("MOMENT", "set", name, "v", labels[vi], "v_norm", vv, "deg", d, "lhs", zp_str(sums[d]),
                 "identity", rstr(num // g, den // g), "ok", 1 if ok else 0, "tag", "MEASURED")
            if not ok and vfirst is None:
                vfirst = d
        per_v.append(vfirst)
        if vfirst is not None and (first_fail is None or vfirst < first_fail):
            first_fail = vfirst
    emit("MOMENT_BREAK", "set", name, "first_deg_fail_some_v", first_fail if first_fail is not None else "none",
         "per_v", ",".join(str(x) if x is not None else "none" for x in per_v), "degs", "2.." + str(MAXDEG), "tag", "MEASURED")
    return first_fail


def emit_distribution(name, dist, s, phi_ring):
    keys = sorted(dist.keys(), key=functools.cmp_to_key(zp_cmp), reverse=True)
    for v in keys:
        emit("DIST", "set", name, "ip_scaled", zp_str(v), "cos", qp_str(qp_norm(v[0], v[1], s)), "count", dist[v], "tag", "MEASURED")
    emit("DIST_SUM", "set", name, "values", len(keys), "total", sum(dist.values()),
         "profile", ",".join(zp_str(v) + ":" + str(dist[v]) for v in keys), "ring", "Z[phi]" if phi_ring else "Z", "tag", "MEASURED")


# ---------------------------------------------------------------- 600-cell
def perm_parity(p):
    inv = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return inv & 1


def build_icosians(odd_swap_one=False):
    pts = []
    for i in range(4):
        for s in (2, -2):
            v = [ZP0, ZP0, ZP0, ZP0]
            v[i] = (s, 0)
            pts.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        pts.append(tuple((sg, 0) for sg in signs))
    for perm in itertools.permutations(range(4)):
        if perm_parity(perm) != 0:
            continue
        for signs in itertools.product((1, -1), repeat=3):
            vals = [ZP0, zp_scale(ZP1, signs[0]), zp_scale(PHIINV, signs[1]), zp_scale(PHI, signs[2])]
            pts.append(tuple(vals[perm[i]] for i in range(4)))
    if odd_swap_one:
        # control: one 96-type vertex replaced by an odd permutation of its coordinates
        v = pts[24]
        pts[24] = (v[1], v[0], v[2], v[3])
    return pts


def q_mul(x, y):
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    m = zp_mul
    r0 = zp_sub(zp_sub(zp_sub(m(a1, a2), m(b1, b2)), m(c1, c2)), m(d1, d2))
    r1 = zp_sub(zp_add(zp_add(m(a1, b2), m(b1, a2)), m(c1, d2)), m(d1, c2))
    r2 = zp_add(zp_add(zp_sub(m(a1, c2), m(b1, d2)), m(c1, a2)), m(d1, b2))
    r3 = zp_add(zp_sub(zp_add(m(a1, d2), m(b1, c2)), m(c1, b2)), m(d1, a2))
    return (r0, r1, r2, r3)


def halve(x):
    out = []
    for c in x:
        if (c[0] & 1) or (c[1] & 1):
            return None
        out.append((c[0] >> 1, c[1] >> 1))
    return tuple(out)


def q_norm(x):
    acc = ZP0
    for c in x:
        acc = zp_add(acc, zp_mul(c, c))
    return acc


def q_conj(x):
    return (x[0], zp_scale(x[1], -1), zp_scale(x[2], -1), zp_scale(x[3], -1))


def ip4(x, y):
    acc = ZP0
    for i in range(4):
        acc = zp_add(acc, zp_mul(x[i], y[i]))
    return acc


def closure_count(pts, ptset):
    inset = 0
    unhalvable = 0
    total = 0
    for x in pts:
        for y in pts:
            total += 1
            h = halve(q_mul(x, y))
            if h is None:
                unhalvable += 1
            elif h in ptset:
                inset += 1
    return inset, unhalvable, total


def order_histogram(pts):
    ident = ((2, 0), ZP0, ZP0, ZP0)
    hist = Counter()
    for x in pts:
        y = x
        k = 1
        while y != ident and k <= 240:
            y = halve(q_mul(y, x))
            if y is None:
                break
            k += 1
        hist[k if y == ident else 0] += 1
    return dict(hist)


def distance_profiles(pts, ipf):
    dists = []
    total = Counter()
    for x in pts:
        c = Counter(ipf(x, y) for y in pts)
        dists.append(frozenset(c.items()))
        total.update(c)
    distinct = len(set(dists))
    return distinct, dict(dists[0]) if dists else {}, total


# ---------------------------------------------------------------- E8
def build_e8():
    pts = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (2, -2):
                for sj in (2, -2):
                    v = [0] * 8
                    v[i] = si
                    v[j] = sj
                    pts.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if (signs.count(-1) & 1) == 0:
            pts.append(signs)
    return pts


def ipz(x, y):
    return (sum(map(operator.mul, x, y)), 0)


# ---------------------------------------------------------------- Golay and Leech
GOLAY_G = (1 << 0) | (1 << 2) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 10) | (1 << 11)
GOLAY_TEXT = "g(x):1+x^2+x^4+x^5+x^6+x^10+x^11;length23_cyclic;rows_x^i*g(x)_i:0..11;plus_overall_parity_bit_23"
GOLAY_EXPECT = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}


def gf2_rem(a, m):
    dm = m.bit_length() - 1
    while a and a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def golay_rows():
    rows = []
    for i in range(12):
        r = GOLAY_G << i
        r |= (r.bit_count() & 1) << 23
        rows.append(r)
    return rows


def golay_span(rows):
    k = len(rows)
    codes = [0] * (1 << k)
    for m in range(1, 1 << k):
        low = (m & -m).bit_length() - 1
        codes[m] = codes[m ^ (1 << low)] ^ rows[low]
    return codes


def weight_enumerator(codes):
    return dict(Counter(c.bit_count() for c in codes))


def wstr(w):
    return ",".join(str(k) + ":" + str(w[k]) for k in sorted(w))


def build_leech(codes):
    octads = [c for c in codes if c.bit_count() == 8]
    t1 = []
    for o in octads:
        pos = [i for i in range(24) if (o >> i) & 1]
        for signs in itertools.product((1, -1), repeat=8):
            if signs.count(-1) & 1:
                continue
            v = [0] * 24
            for p, sg in zip(pos, signs):
                v[p] = 2 * sg
            t1.append(tuple(v))
    t3 = []
    for c in codes:
        y = [(-1 if (c >> i) & 1 else 1) for i in range(24)]
        for j in range(24):
            v = list(y)
            v[j] = -3 * y[j]
            t3.append(tuple(v))
    t4 = []
    for i in range(24):
        for j in range(i + 1, 24):
            for si in (4, -4):
                for sj in (4, -4):
                    v = [0] * 24
                    v[i] = si
                    v[j] = sj
                    t4.append(tuple(v))
    return octads, t1, t3, t4


def leech_member(v, codeset):
    m = v[0] & 1
    for a in v:
        if (a & 1) != m:
            return False
    if (sum(v) - 4 * m) & 7:
        return False
    target = 2 if m == 0 else 3
    mask = 0
    for i in range(24):
        if (v[i] & 3) == target:
            mask |= 1 << i
    return mask in codeset


LEECH_IPS = (32, 16, 8, 0, -8, -16, -32)


def leech_all_pairs(vecs, stride, shape_of):
    # For base x, T = sum_y (<x,y> + 32) 256^pos(y), packed in one integer; its base-256 digits are the inner products.
    N = len(vecs)
    ones = int.from_bytes(b"\x01" * N, "little")
    P = []
    for i in range(24):
        col = bytes(v[i] + 4 for v in vecs)
        P.append(int.from_bytes(col, "little") - 4 * ones)
    tables = [dict() for _ in range(6)]
    base = 32 * ones
    codes = tuple(ip + 32 for ip in LEECH_IPS)
    zkey = (0, 0, 0, 0)
    ref = None
    tally = [0] * len(LEECH_IPS)
    rows = 0
    rows_bad_set = 0
    rows_bad_self = 0
    rows_diff = 0
    shapes = Counter()
    for a in range(0, N, stride):
        x = vecs[a]
        T = base
        for b in range(6):
            key = x[4 * b:4 * b + 4]
            if key == zkey:
                continue
            t = tables[b].get(key)
            if t is None:
                t = 0
                for o in range(4):
                    if key[o]:
                        t += key[o] * P[4 * b + o]
                tables[b][key] = t
            T += t
        rows += 1
        shapes[shape_of(x)] += 1
        try:
            d = T.to_bytes(N, "little")
        except OverflowError:
            rows_bad_set += 1
            continue
        cnt = tuple(d.count(c) for c in codes)
        if sum(cnt) != N:
            rows_bad_set += 1
        if d[a] != 64:
            rows_bad_self += 1
        if ref is None:
            ref = cnt
        elif cnt != ref:
            rows_diff += 1
        for i in range(len(cnt)):
            tally[i] += cnt[i]
        if (rows & 16383) == 0:
            log("leech rows " + str(rows))
    return rows, rows_bad_set, rows_bad_self, rows_diff, ref, tally, dict(shapes)


def leech_shape(x):
    nz = sum(1 for a in x if a)
    if nz == 8:
        return "2^8"
    if nz == 24:
        return "3.1^23"
    return "4^2"


# ---------------------------------------------------------------- main
def parse_args():
    stride = 1
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--leech-stride" and i + 1 < len(args):
            stride = int(args[i + 1])
            i += 2
        else:
            raise SystemExit("unknown argument " + args[i])
    if stride < 1:
        raise SystemExit("stride must be >= 1")
    return stride


def main():
    stride = parse_args()
    emit("MEASURE", "file", OUT_NAME, "script", "measure_exceptional.py", "repo", "JesseBrown1980/Iconic-3-space", "main", "8f355db",
         "E", 0, "network", 0, "float_used", 0, "kmax", KMAX, "moment_degs", "2.." + str(MAXDEG), "leech_stride", stride,
         "independence", "predict_scripts_not_read", "tag", "MEASURED")

    Q4 = check_zonal_polys(4)
    Q8 = check_zonal_polys(8)
    Q24 = check_zonal_polys(24)
    log("zonal polynomials built")

    # ---------------- 600-cell
    ico = build_icosians()
    S600 = set(ico)
    check("600.count", len(ico) == 120, "count", len(ico))
    check("600.distinct", len(S600) == 120, "distinct", len(S600))
    bad_norm = sum(1 for x in ico if q_norm(x) != (4, 0))
    check("600.norm_scaled_4", bad_norm == 0, "bad", bad_norm)
    inset, unhalv, total = closure_count(ico, S600)
    check("600.closure_quaternion_product", inset == total, "products", total, "in_set", inset, "not_halvable", unhalv)
    check("600.identity_present", ((2, 0), ZP0, ZP0, ZP0) in S600)
    conj_missing = sum(1 for x in ico if q_conj(x) not in S600)
    check("600.inverses_present", conj_missing == 0, "missing", conj_missing)
    hist = order_histogram(ico)
    want_hist = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}
    check("600.order_histogram_of_2I", hist == want_hist, "hist", ",".join(str(k) + ":" + str(hist[k]) for k in sorted(hist)))
    nonab = sum(1 for x in ico for y in ico if q_mul(x, y) != q_mul(y, x))
    emit("GROUP", "set", "600cell", "order", len(S600), "closed", 1 if inset == total else 0, "noncommuting_ordered_pairs", nonab,
         "order_hist", ",".join(str(k) + ":" + str(hist[k]) for k in sorted(hist)),
         "identification", "2I_SL(2,5)_by_closure_and_order_histogram", "tag", "MEASURED;identification:NAMED")
    distinct600, dist600, tally600 = distance_profiles(ico, ip4)
    check("600.distance_invariance_every_point", distinct600 == 1, "base_points", len(ico), "distinct_distributions", distinct600)
    ipset600 = set(tally600.keys())
    emit_distribution("600cell", dist600, 4, True)
    emit("PAIRS", "set", "600cell", "pairs", sum(tally600.values()), "base_points", len(ico), "tag", "MEASURED")
    S600k, t600, nz600 = design_report("600cell", 4, 4, 120, tally600, 120, Q4)
    log("600-cell done")

    # ---------------- E8
    e8 = build_e8()
    SE8 = set(e8)
    check("E8.count", len(e8) == 240, "count", len(e8))
    check("E8.distinct", len(SE8) == 240, "distinct", len(SE8))
    bad8 = sum(1 for x in e8 if sum(a * a for a in x) != 8)
    check("E8.norm_scaled_8", bad8 == 0, "bad", bad8)
    distinct8, dist8, tally8 = distance_profiles(e8, ipz)
    ipset8 = sorted(v[0] for v in tally8)
    check("E8.inner_products_in_{8,4,0,-4,-8}", set(ipset8) <= {8, 4, 0, -4, -8}, "seen", ",".join(str(v) for v in ipset8))
    check("E8.distance_invariance_every_point", distinct8 == 1, "base_points", len(e8), "distinct_distributions", distinct8)
    refl_bad = 0
    for x in e8:
        for r in e8:
            ip = sum(map(operator.mul, x, r))
            if ip & 3:
                refl_bad += 1
                continue
            c = ip >> 2
            y = tuple(x[i] - c * r[i] for i in range(8))
            if y not in SE8:
                refl_bad += 1
    check("E8.closed_under_reflections", refl_bad == 0, "pairs", len(e8) * len(e8), "bad", refl_bad)
    emit_distribution("E8", dist8, 8, False)
    emit("PAIRS", "set", "E8", "pairs", sum(tally8.values()), "base_points", len(e8), "tag", "MEASURED")
    S8k, t8, nz8 = design_report("E8", 8, 8, 240, tally8, 240, Q8)
    log("E8 done")

    # ---------------- Golay
    check("golay.g_divides_x^23+1", gf2_rem((1 << 23) | 1, GOLAY_G) == 0)
    grows = golay_rows()
    codes = golay_span(grows)
    codeset = set(codes)
    wen = weight_enumerator(codes)
    check("golay.codewords_4096", len(codeset) == 4096, "distinct", len(codeset))
    check("golay.weight_enumerator", wen == GOLAY_EXPECT, "weights", wstr(wen))
    selfdual = all((grows[i] & grows[j]).bit_count() % 2 == 0 for i in range(12) for j in range(12))
    check("golay.generator_rows_self_orthogonal", selfdual)
    emit("GOLAY", "generator", GOLAY_TEXT, "codewords", len(codeset), "weights", wstr(wen), "tag", "MEASURED")

    # ---------------- Leech
    octads, t1, t3, t4 = build_leech(codes)
    vecs = t1 + t3 + t4
    N = len(vecs)
    check("leech.octads_759", len(octads) == 759, "octads", len(octads))
    check("leech.shape_2^8_97152", len(t1) == 97152, "count", len(t1))
    check("leech.shape_3_1^23_98304", len(t3) == 98304, "count", len(t3))
    check("leech.shape_4^2_1104", len(t4) == 1104, "count", len(t4))
    SL = set(vecs)
    check("leech.distinct_196560", len(SL) == 196560 and N == 196560, "count", N, "distinct", len(SL))
    badn = sum(1 for x in vecs if sum(a * a for a in x) != 32)
    check("leech.norm_scaled_32", badn == 0, "bad", badn)
    nonmem = sum(1 for x in vecs if not leech_member(x, codeset))
    check("leech.lattice_membership_conway", nonmem == 0, "nonmembers", nonmem)
    antip = sum(1 for x in vecs if tuple(-a for a in x) not in SL)
    check("leech.antipodal", antip == 0, "missing_negatives", antip)
    log("leech built, all-pairs pass")
    rows, bad_set, bad_self, rdiff, ref, ltally, shapes = leech_all_pairs(vecs, stride, leech_shape)
    check("leech.every_inner_product_in_{+-32,+-16,+-8,0}", bad_set == 0, "base_points", rows, "pairs", rows * N, "bad_rows", bad_set)
    check("leech.self_inner_product_32", bad_self == 0, "bad_rows", bad_self)
    check("leech.distance_invariance", rdiff == 0 and rows >= 30, "base_points", rows, "rows_differing", rdiff,
          "shapes", ",".join(k + ":" + str(shapes[k]) for k in sorted(shapes)))
    check("leech.base_points_cover_three_shapes", len(shapes) == 3, "shapes", len(shapes))
    distL = {}
    for i in range(len(LEECH_IPS)):
        distL[(LEECH_IPS[i], 0)] = ref[i]
    emit_distribution("Leech", distL, 32, False)
    emit("PAIRS", "set", "Leech", "pairs", sum(ltally), "base_points", rows, "stride", stride, "full_double_sum", 1 if stride == 1 else 0, "tag", "MEASURED")
    tallyL = {}
    for i in range(len(LEECH_IPS)):
        tallyL[(LEECH_IPS[i], 0)] = ltally[i]
    SLk, tL, nzL = design_report("Leech", 24, 32, N, tallyL, rows, Q24)
    log("leech zonal done")

    # ---------------- moment identities
    tv600 = [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 2, 3, 4), (5, -2, 7, 3)]
    mc = []
    for v in tv600:
        c = Counter()
        for x in ico:
            acc = ZP0
            for i in range(4):
                acc = zp_add(acc, zp_scale(x[i], v[i]))
            c[acc] += 1
        mc.append((c, sum(a * a for a in v)))
    mb600 = moment_report("600cell", 4, 4, 120, mc, [",".join(str(a) for a in v) for v in tv600])
    tv8 = [(1, 0, 0, 0, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0, 0, 0), (2, 2, 0, 0, 0, 0, 0, 0), (1, 2, 3, 4, 5, 6, 7, 8), (3, 1, -4, 1, 5, -9, 2, 6)]
    mc = []
    for v in tv8:
        mc.append((Counter((sum(map(operator.mul, x, v)), 0) for x in e8), sum(a * a for a in v)))
    mb8 = moment_report("E8", 8, 8, 240, mc, [",".join(str(a) for a in v) for v in tv8])
    tvL = [tuple([1] + [0] * 23), tuple([1, 1] + [0] * 22), tuple(range(1, 25)),
           (2, -1, 0, 3, 1, -4, 5, 0, -2, 7, 1, 1, -3, 0, 6, 2, -5, 1, 0, 4, -1, 3, 2, -6),
           tuple([4, 4] + [0] * 22), tuple([-3] + [1] * 23)]
    mc = []
    for v in tvL:
        mc.append((Counter((sum(map(operator.mul, x, v)), 0) for x in vecs), sum(a * a for a in v)))
    mbL = moment_report("Leech", 24, 32, N, mc, [",".join(str(a) for a in v) for v in tvL])
    log("moments done")

    for name, t, nz, mb in (("600cell", t600, nz600, mb600), ("E8", t8, nz8, mb8), ("Leech", tL, nzL, mbL)):
        emit("FIRST_BREAK", "set", name, "zonal_first_nonzero_k", nz[0] if nz else "none", "t", t,
             "moment_first_fail_deg", mb if mb is not None else "none", "tag", "MEASURED")

    # ---------------- controls that must fail
    def control(name, check_passed, *kv):
        ok = not check_passed
        if not ok:
            FAILED.append("control." + name)
        emit("CONTROL", "name", name, "check_passed", 1 if check_passed else 0, "control_failed_as_required", 1 if ok else 0, *kv)

    brows = list(grows)
    brows[0] ^= 1 << 5
    bcodes = golay_span(brows)
    bwen = weight_enumerator(bcodes)
    control("golay_flip_row0_bit5", bwen == GOLAY_EXPECT and len(set(bcodes)) == 4096, "weights", wstr(bwen), "distinct", len(set(bcodes)))

    e8d = e8[1:]
    dd, _, _ = distance_profiles(e8d, ipz)
    control("e8_drop_one_root.count", len(e8d) == 240, "count", len(e8d))
    control("e8_drop_one_root.distance_invariance", dd == 1, "distinct_distributions", dd)

    o = octads[0]
    pos = [i for i in range(24) if (o >> i) & 1]
    vb = [0] * 24
    for p in pos:
        vb[p] = 2
    vb[pos[0]] = -2
    vb = tuple(vb)
    ipb = Counter(sum(map(operator.mul, vb, y)) for y in vecs)
    control("leech_odd_sign_octad.inner_products", set(ipb.keys()) <= set(LEECH_IPS), "seen", ",".join(str(k) for k in sorted(ipb)))
    control("leech_odd_sign_octad.membership", leech_member(vb, codeset))

    Q3 = zonal_polys(3, KMAX)[0]
    Sw = zonal_sums(tally600, 4, Q3, 120)
    tw = 0
    while tw + 1 <= KMAX and qp_is_zero(Sw[tw + 1]):
        tw += 1
    control("zonal_wrong_dimension_n3_on_600cell.t_ge_11", tw >= 11, "t_seen", tw, "S_2", qp_str(Sw[2]))

    icob = build_icosians(odd_swap_one=True)
    ib, _, tb = closure_count(icob, set(icob))
    control("icosian_one_odd_permutation.closure", ib == tb, "in_set", ib, "products", tb)

    cv = Counter((sum(map(operator.mul, x, tvL[2])), 0) for x in vecs)
    s2 = sum(c * v[0] * v[0] for v, c in cv.items())
    vv = sum(a * a for a in tvL[2])
    control("leech_moment_wrong_dimension_n23.deg2", s2 * 23 == N * 32 * vv, "lhs", s2, "identity_n23", rstr(N * 32 * vv, 23))

    emit("SUMMARY", "checks", NCHECK[0], "failed", len(FAILED), "failed_names", ",".join(FAILED) if FAILED else "none", "tag", "MEASURED")
    measured_rows = len(ROWS)
    M = {
        "S": {"600cell": S600k, "E8": S8k, "Leech": SLk},
        "t": {"600cell": t600, "E8": t8, "Leech": tL},
        "nz": {"600cell": nz600, "E8": nz8, "Leech": nzL},
        "N": {"600cell": len(ico), "E8": len(e8), "Leech": N},
        "R": {"600cell": len(ico), "E8": len(e8), "Leech": rows},
        "dist": {"600cell": (dist600, 4), "E8": (dist8, 8), "Leech": (distL, 32)},
        "mb": {"600cell": mb600, "E8": mb8, "Leech": mbL},
        "antipodal": {"600cell": all(tuple(zp_scale(c, -1) for c in x) in S600 for x in ico),
                      "E8": all(tuple(-a for a in x) in SE8 for x in e8), "Leech": antip == 0},
    }
    ok_frozen, n_dis = compare_stage(M, measured_rows)
    with open(os.path.abspath(__file__), "rb") as f:
        script_bytes = f.read()
    emit("SCRIPT", "name", "measure_exceptional.py", "sha256", hashlib.sha256(script_bytes).hexdigest(), "bytes", len(script_bytes), "json_used", 0, "float_used", 0)
    emit("END", "rows", len(ROWS) + 1, "measured_rows", measured_rows, "checks_failed", len(FAILED), "frozen_verified", 1 if ok_frozen else 0,
         "disagreements", n_dis, "exit", 1 if (FAILED or not ok_frozen) else (2 if n_dis else 0))
    write_out()
    if FAILED or not ok_frozen:
        log("FAILED: " + ",".join(FAILED) + (" frozen_not_verified" if not ok_frozen else ""))
        sys.exit(1)
    log("all checks passed; disagreements with frozen rows: " + str(n_dis))
    if n_dis:
        sys.exit(2)


# ---------------------------------------------------------------- COMPARE stage (reads the frozen rows; runs after SUMMARY)
FROZEN_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "iconic-3-space", "repo", "predictions"))
FROZEN_EXC = ("exceptional", "PREDICTION-EXCEPTIONAL.hbp", "e8403cd4")
FROZEN_SU2 = ("su2-ladder", "PREDICTION.hbp", "c2601cf0")
SETMAP = {"600cell_2I_icosians": "600cell", "E8_roots": "E8", "Leech_min_vectors": "Leech"}
CMP_STATS = [0, 0, []]


def rred(n, d):
    if d < 0:
        n, d = -n, -d
    g = gcd(n, d)
    return (n // g, d // g)


def prat(s):
    a, b = s.split(":")
    return rred(int(a), int(b))


def rs(r):
    return str(r[0]) + ":" + str(r[1])


def ps5(s):
    # "a+br5" with a, b as num:den
    if not s.endswith("r5"):
        raise ValueError("not a sqrt5 value " + s)
    a, b = s[:-2].split("+")
    return (prat(a), prat(b))


def s5s(v):
    return rs(v[0]) + "+" + rs(v[1]) + "r5"


def qp_to_s5(x):
    # (A + B*phi)/d with phi = (1 + sqrt5)/2  ->  (2A + B)/(2d) + (B/(2d)) sqrt5
    return (rred(2 * x[0] + x[1], 2 * x[2]), rred(x[1], 2 * x[2]))


def ilist(s):
    if s in ("", "none"):
        return ()
    return tuple(int(a) for a in s.split(","))


def lstr(t):
    return ",".join(str(a) for a in t) if t else "none"


def cmp_emit(row, setname, quantity, fval, mval, fstr, mstr, *extra):
    agree = fval == mval
    CMP_STATS[0] += 1
    if agree:
        CMP_STATS[1] += 1
    else:
        CMP_STATS[2].append(row + ":" + setname + ":" + quantity)
    emit("COMPARE", "row", row, "set", setname, "quantity", quantity, "frozen", fstr, "measured", mstr, "agree", 1 if agree else 0,
         *extra, "tag", "MEASURED_vs_PREDICTION")


def parse_hbp(text):
    out = []
    for line in text.split("\n"):
        if not line:
            continue
        parts = line.split("|")
        d = {}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                d[k] = v
        out.append((parts[0], d))
    return out


def read_frozen(spec):
    path = os.path.join(FROZEN_DIR, spec[0], spec[1])
    try:
        with open(path, "rb") as f:
            data = f.read()
        with open(path + ".sha256", "rb") as f:
            side = f.read().decode("utf-8").split()
    except OSError:
        emit("FROZEN", "file", spec[0] + "/" + spec[1], "missing", 1, "ok", 0)
        return None
    h = hashlib.sha256(data).hexdigest()
    ok = len(side) >= 2 and side[0] == h and side[1] == spec[1] and h.startswith(spec[2])
    emit("FROZEN", "file", spec[0] + "/" + spec[1], "sha256", h, "sidecar_match", 1 if (side and side[0] == h) else 0,
         "prefix_expected", spec[2], "prefix_match", 1 if h.startswith(spec[2]) else 0, "bytes", len(data), "ok", 1 if ok else 0)
    if not ok:
        return None
    return parse_hbp(data.decode("utf-8"))


def molien_2i(k):
    # coefficient of t^k in (1 + t^30) / ((1 - t^12)(1 - t^20)), counted directly
    c = 0
    for kk in (k, k - 30):
        if kk < 0:
            continue
        for a in range(kk // 12 + 1):
            if (kk - 12 * a) % 20 == 0:
                c += 1
    return c


def compare_stage(M, measured_rows):
    block = ("\n".join(ROWS[:measured_rows]) + "\n").encode("utf-8")
    hb = hashlib.sha256(block).hexdigest()
    pre = os.path.join(HERE, "MEASURE-EXCEPTIONAL.precompare.hbp")
    same = None
    try:
        with open(pre, "rb") as f:
            same = f.read() == block
    except OSError:
        same = None
    emit("PRECOMPARE", "measured_rows", measured_rows, "sha256_of_measured_block", hb,
         "precompare_file", "MEASURE-EXCEPTIONAL.precompare.hbp", "precompare_present", 0 if same is None else 1,
         "identical", "unknown" if same is None else (1 if same else 0),
         "meaning", "measured_rows_written_before_the_frozen_rows_were_opened")
    if same is False:
        FAILED.append("precompare_block_changed")
    exc = read_frozen(FROZEN_EXC)
    su2 = read_frozen(FROZEN_SU2)
    if exc is None:
        return False, 0
    S = M["S"]
    for tag, d in exc:
        if tag not in ("SET", "SK", "TPRED", "FIRST_BREAK", "PROFILE", "NAMED_ZEROS", "CROSS"):
            continue
        fset = d.get("set", "")
        ms = SETMAP.get(fset)
        if ms is None:
            emit("COMPARE_SKIP", "row", tag, "set", fset, "reason", "set_name_not_mapped")
            continue
        if tag == "SET":
            cmp_emit("SET", ms, "N", int(d["N"]), M["N"][ms], d["N"], str(M["N"][ms]))
            dist, s = M["dist"][ms]
            md = {}
            for v, c in dist.items():
                md[qp_to_s5(qp_norm(v[0], v[1], s))] = c
            fd = {}
            for ent in d["dist"].split(","):
                val, c = ent.split("#")
                fd[ps5(val)] = int(c)
            order = sorted(dist.keys(), key=functools.cmp_to_key(zp_cmp), reverse=True)
            mtext = ",".join(s5s(qp_to_s5(qp_norm(v[0], v[1], s))) + "#" + str(dist[v]) for v in order)
            cmp_emit("SET", ms, "distance_distribution", fd, md, d["dist"], mtext, "rows_frozen", len(fd), "rows_measured", len(md))
            cmp_emit("SET", ms, "counts_sum", int(d["counts_sum"]), sum(md.values()), d["counts_sum"], str(sum(md.values())))
            cmp_emit("SET", ms, "antipodal", int(d["antipodal"]), 1 if M["antipodal"][ms] else 0, d["antipodal"], "1" if M["antipodal"][ms] else "0")
            gal = all(md.get((k[0], (-k[1][0], k[1][1]))) == md[k] for k in md)
            cmp_emit("SET", ms, "galois_stable", int(d["galois_stable"]), 1 if gal else 0, d["galois_stable"], "1" if gal else "0")
        elif tag == "SK":
            k = int(d["k"])
            fv = (prat(d["S_a"]), prat(d["S_b"]))
            mv = qp_to_s5(S[ms][k])
            fz = int(d["zero"])
            mz = 1 if qp_is_zero(S[ms][k]) else 0
            fsg = int(d["sign"])
            msg = zp_sign((S[ms][k][0], S[ms][k][1]))
            cmp_emit("SK", ms, "S_" + str(k), (fv, fz, fsg), (mv, mz, msg), s5s(fv), s5s(mv), "k", k,
                     "zero_frozen", fz, "zero_measured", mz, "sign_frozen", fsg, "sign_measured", msg)
        elif tag == "TPRED":
            cmp_emit("TPRED", ms, "t", int(d["t_pred"]), M["t"][ms], d["t_pred"], str(M["t"][ms]))
            cmp_emit("TPRED", ms, "N", int(d["N"]), M["N"][ms], d["N"], str(M["N"][ms]))
        elif tag == "FIRST_BREAK":
            k = int(d["k"])
            mk = M["nz"][ms][0] if M["nz"][ms] else -1
            cmp_emit("FIRST_BREAK", ms, "k", k, mk, d["k"], str(mk))
            fv = (prat(d["S_a"]), prat(d["S_b"]))
            mv = qp_to_s5(S[ms][k])
            cmp_emit("FIRST_BREAK", ms, "S_k", fv, mv, s5s(fv), s5s(mv), "k", k)
            R = M["R"][ms]
            ps = qp_rat(S[ms][k], R, 1)
            fp = (prat(d["pair_sum_a"]), prat(d["pair_sum_b"]))
            mp = qp_to_s5(ps)
            cmp_emit("FIRST_BREAK", ms, "pair_sum", fp, mp, s5s(fp), s5s(mp), "k", k, "base_points", R, "full_double_sum", 1 if R == M["N"][ms] else 0)
            cmp_emit("FIRST_BREAK", ms, "sign", int(d["sign"]), zp_sign((S[ms][k][0], S[ms][k][1])), d["sign"], str(zp_sign((S[ms][k][0], S[ms][k][1]))))
            mb = M["mb"][ms]
            cmp_emit("FIRST_BREAK", ms, "k_vs_moment_identity_first_fail_deg", k, mb if mb is not None else -1, d["k"], str(mb),
                     "method", "independent_moment_test_2m_2..14")
        elif tag == "PROFILE":
            kmax = int(d["kmax"])
            mnz = tuple(k for k in range(1, kmax + 1) if not qp_is_zero(S[ms][k]))
            mz = tuple(k for k in range(1, kmax + 1) if qp_is_zero(S[ms][k]))
            cmp_emit("PROFILE", ms, "nonzero_k", ilist(d["nonzero_k"]), mnz, d["nonzero_k"], lstr(mnz))
            cmp_emit("PROFILE", ms, "zero_k", ilist(d["zero_k"]), mz, d["zero_k"], lstr(mz))
            cmp_emit("PROFILE", ms, "nonzero_count", int(d["nonzero_count"]), len(mnz), d["nonzero_count"], str(len(mnz)))
            fvals = {}
            for ent in d["values"].split(","):
                kk, val = ent.split("@")
                fvals[int(kk)] = ps5(val)
            mvals = {}
            for k in mnz:
                mvals[k] = qp_to_s5(S[ms][k])
            cmp_emit("PROFILE", ms, "values", fvals, mvals, d["values"], ",".join(str(k) + "@" + s5s(mvals[k]) for k in mnz))
        elif tag == "NAMED_ZEROS":
            fz = ilist(d["forced_zero_k"])
            held = sum(1 for k in fz if qp_is_zero(S[ms][k]))
            cmp_emit("NAMED_ZEROS", ms, "forced_zero_k_measured_zero[" + d.get("source", "") + "]", len(fz), held,
                     str(len(fz)) + ":" + str(len(fz)), str(held) + ":" + str(len(fz)))
            if "invariant_harmonic_k" in d:
                mnz = tuple(M["nz"][ms])
                cmp_emit("NAMED_ZEROS", ms, "invariant_harmonic_k[" + d.get("source", "") + "]", ilist(d["invariant_harmonic_k"]), mnz,
                         d["invariant_harmonic_k"], lstr(mnz))
        elif tag == "CROSS":
            mk = M["nz"][ms][0] if M["nz"][ms] else -1
            for q in ("pred_first_break_k", "su2_ladder_first_invariant_degree", "molien_first_invariant_degree"):
                if q in d:
                    cmp_emit("CROSS", ms, q, int(d[q]), mk, d[q], str(mk), "measured_as", "first_k_with_S_k_nonzero")
            if ms == "600cell":
                good = 0
                bad_int = 0
                mm = []
                for k in range(1, KMAX + 1):
                    x = qp_rat(S[ms][k], k + 1, 120)
                    if x[1] != 0 or x[2] != 1:
                        bad_int += 1
                        mm.append(None)
                        continue
                    mm.append(x[0])
                    if x[0] == molien_2i(k):
                        good += 1
                pde = d.get("per_degree_exact", "")
                cmp_emit("CROSS", ms, "per_degree_exact[S_k(k+1)_eq_120m_k,m_k_Molien_(1+t^30)over((1-t^12)(1-t^20))]",
                         pde, str(good) + ":" + str(KMAX), pde, str(good) + ":" + str(KMAX), "non_integer_m_k", bad_int,
                         "m_k_measured_k1..30", ",".join("x" if v is None else str(v) for v in mm))
                if su2 is not None:
                    for tg, dd in su2:
                        if dd.get("group") != "2I":
                            continue
                        if tg == "MOLIEN":
                            fa = ilist(dd["a"])
                            mlist = [1] + [v if v is not None else -1 for v in mm]
                            nk = min(len(fa), KMAX + 1)
                            cmp_emit("CROSS_SU2", ms, "MOLIEN_a_d_for_d_0.." + str(nk - 1) + "_vs_measured_S_d(d+1)over120",
                                     tuple(fa[:nk]), tuple(mlist[:nk]), lstr(fa[:nk]), lstr(tuple(mlist[:nk])), "frozen_file", "su2-ladder/PREDICTION.hbp")
                        elif tg == "DESIGN":
                            cmp_emit("CROSS_SU2", ms, "DESIGN_t_pred", int(dd["t_pred"]), M["t"][ms], dd["t_pred"], str(M["t"][ms]), "frozen_file", "su2-ladder/PREDICTION.hbp")
                            cmp_emit("CROSS_SU2", ms, "DESIGN_first_break_degree", int(dd["first_break_degree"]), mk, dd["first_break_degree"], str(mk), "frozen_file", "su2-ladder/PREDICTION.hbp")
                            mab = qp_rat(S[ms][mk], mk + 1, 120)
                            mabv = mab[0] if (mab[1] == 0 and mab[2] == 1) else -1
                            cmp_emit("CROSS_SU2", ms, "DESIGN_m_at_break", int(dd["m_at_break"]), mabv, dd["m_at_break"], str(mabv), "frozen_file", "su2-ladder/PREDICTION.hbp")
    emit("COMPARE_SUMMARY", "compared", CMP_STATS[0], "agree", CMP_STATS[1], "disagree", len(CMP_STATS[2]),
         "disagreements", ",".join(CMP_STATS[2]) if CMP_STATS[2] else "none", "rule", "a_disagreement_is_a_result_not_edited", "tag", "MEASURED_vs_PREDICTION")
    return True, len(CMP_STATS[2])


def write_out():
    path = os.path.join(HERE, OUT_NAME)
    data = ("\n".join(ROWS) + "\n").encode("utf-8")
    with open(path, "wb") as f:
        f.write(data)
    h = hashlib.sha256(data).hexdigest()
    with open(path + ".sha256", "wb") as f:
        f.write((h + "  " + OUT_NAME + "\n").encode("utf-8"))
    log("wrote " + path + " sha256 " + h)


if __name__ == "__main__":
    main()
