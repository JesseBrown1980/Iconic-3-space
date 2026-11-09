#!/usr/bin/env python3
"""
refute_measure.py - second-seat (refuter) independent recompute for JesseBrown1980/Iconic-3-space main 8f355db.

Written without opening predict_*.py or the three measurement scripts' code paths for the values. Every group and point
set is built here from its definition. Integers only: rationals are reduced (num, den) pairs, Z[phi] values are (a, b)
meaning a + b*phi with phi^2 = phi + 1, Z[sqrt2] values are (u, v) meaning u + v*sqrt2, Z[zeta_M] values are coefficient
lists reduced mod the cyclotomic polynomial. No float, no true division, no '**', no json.

Phase 1 (measure) writes every measured row to REFUTE-MEASURE.precompare.hbp before any frozen file is opened.
Phase 2 (compare) reads the frozen PREDICTION rows and appends COMPARE rows. Exit 1 on any failed assertion or any
control that fails to fail; exit 2 if any COMPARE row disagrees; exit 0 otherwise.
Usage: python3 -I -B refute_measure.py [--measure-only]
"""
import sys
import os
import hashlib
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "REFUTE-MEASURE.hbp"
PRE_NAME = "REFUTE-MEASURE.precompare.hbp"
ROWS = []
FAILS = []
KMAX = 30
NMAX = 30


def emit(tag, *kv):
    parts = [tag]
    i = 0
    while i < len(kv):
        v = str(kv[i + 1])
        if "|" in v or "\n" in v:
            FAILS.append("row_format:" + tag + ":" + str(kv[i]))
            v = v.replace("|", "/").replace("\n", " ")
        parts.append(str(kv[i]) + "=" + v)
        i += 2
    parts.append("json=0")
    ROWS.append("|".join(parts))


def check(cond, name):
    if not cond:
        FAILS.append(name)
    return 1 if cond else 0


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


# ---------------------------------------------------------------- rationals
def rq(n, d):
    if d == 0:
        raise ZeroDivisionError("rq")
    if d < 0:
        n, d = -n, -d
    g = gcd(n, d)
    if g > 1:
        n //= g
        d //= g
    return (n, d)


def radd(a, b):
    return rq(a[0] * b[1] + b[0] * a[1], a[1] * b[1])


def rsub(a, b):
    return rq(a[0] * b[1] - b[0] * a[1], a[1] * b[1])


def rmul(a, b):
    return rq(a[0] * b[0], a[1] * b[1])


def rdiv(a, b):
    return rq(a[0] * b[1], a[1] * b[0])


def rstr(a):
    return str(a[0]) + ":" + str(a[1])


def binom(n, k):
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


def ipow(b, e):
    if e < 0:
        raise ValueError("negative exponent")
    r = 1
    for _ in range(e):
        r *= b
    return r


# ---------------------------------------------------------------- rings
# Z
Z_ADD = lambda x, y: x + y
Z_SUB = lambda x, y: x - y
Z_MUL = lambda x, y: x * y


def z_half(x):
    return x >> 1 if x % 2 == 0 else None


# Z[sqrt2]: (u, v) = u + v*sqrt2
def s2_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def s2_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def s2_mul(x, y):
    return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def s2_half(x):
    if x[0] % 2 or x[1] % 2:
        return None
    return (x[0] >> 1, x[1] >> 1)


# Z[phi]: (a, b) = a + b*phi, phi^2 = phi + 1
def ph_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def ph_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def ph_mul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def ph_half(x):
    if x[0] % 2 or x[1] % 2:
        return None
    return (x[0] >> 1, x[1] >> 1)


RING_Z = (0, 1, Z_ADD, Z_SUB, Z_MUL, z_half)
RING_S2 = ((0, 0), (1, 0), s2_add, s2_sub, s2_mul, s2_half)
RING_PH = ((0, 0), (1, 0), ph_add, ph_sub, ph_mul, ph_half)


# ---------------------------------------------------------------- doubled quaternions X = 2q over a ring
def qmul_raw(R, x, y):
    add, sub, mul = R[2], R[3], R[4]
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    a = sub(sub(sub(mul(a1, a2), mul(b1, b2)), mul(c1, c2)), mul(d1, d2))
    b = sub(add(add(mul(a1, b2), mul(b1, a2)), mul(c1, d2)), mul(d1, c2))
    c = add(add(sub(mul(a1, c2), mul(b1, d2)), mul(c1, a2)), mul(d1, b2))
    d = add(sub(add(mul(a1, d2), mul(b1, c2)), mul(c1, b2)), mul(d1, a2))
    return (a, b, c, d)


def dq_mul(R, X, Y):
    # (2p)(2q) = 4pq ; return 2pq or None if not divisible
    P = qmul_raw(R, X, Y)
    out = []
    for c in P:
        h = R[5](c)
        if h is None:
            return None
        out.append(h)
    return tuple(out)


def dq_conj(R, X):
    z = R[0]
    sub = R[3]
    return (X[0], sub(z, X[1]), sub(z, X[2]), sub(z, X[3]))


def dq_norm(R, X):
    add, mul = R[2], R[4]
    s = R[0]
    for c in X:
        s = add(s, mul(c, c))
    return s


def perms4():
    out = []
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if len({a, b, c, d}) == 4:
                        out.append((a, b, c, d))
    return out


def perm_parity(p):
    inv = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return inv % 2


def lift(R, v):
    # integer -> ring element
    if R is RING_Z:
        return v
    return (v, 0)


def hurwitz_units(R):
    out = set()
    for pos in range(4):
        for s in (2, -2):
            v = [lift(R, 0)] * 4
            v[pos] = lift(R, s)
            out.add(tuple(v))
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    out.add(tuple(lift(R, s) for s in (s0, s1, s2, s3)))
    return out


def build_2T():
    return sorted(hurwitz_units(RING_Z))


def build_2O():
    out = set(hurwitz_units(RING_S2))
    # (+-1, +-1, 0, 0)/sqrt2 in all positions; doubled: sqrt2 * (+-1, +-1, 0, 0)
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [(0, 0)] * 4
                    v[i] = (0, si)
                    v[j] = (0, sj)
                    out.add(tuple(v))
    return sorted(out)


def build_2I():
    out = set(hurwitz_units(RING_PH))
    # even permutations of (0, +-1, +-phi^-1, +-phi)/2 ; doubled: (0, +-1, +-(phi-1), +-phi)
    base = [(0, 0), (1, 0), (-1, 1), (0, 1)]
    for p in perms4():
        if perm_parity(p):
            continue
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    vals = [base[0], (s1 * base[1][0], s1 * base[1][1]), (s2 * base[2][0], s2 * base[2][1]),
                            (s3 * base[3][0], s3 * base[3][1])]
                    v = [None] * 4
                    for k in range(4):
                        v[p[k]] = vals[k]
                    out.add(tuple(v))
    return sorted(out)


def check_quat_group(R, G, name, expect_order):
    S = set(G)
    one = (lift(R, 2), lift(R, 0), lift(R, 0), lift(R, 0))
    four = lift(R, 4)
    norms_ok = all(dq_norm(R, X) == four for X in G)
    closed = 0
    for X in G:
        for Y in G:
            P = dq_mul(R, X, Y)
            if P is not None and P in S:
                closed += 1
    inv_ok = all(dq_conj(R, X) in S for X in G)
    res = dict(order=len(S), distinct=1 if len(S) == len(G) else 0, identity=1 if one in S else 0, unit_norm=1 if norms_ok else 0,
               closed_products=closed, pairs=len(G) * len(G), inverses=1 if inv_ok else 0)
    ok = (len(S) == expect_order and res["distinct"] and res["identity"] and norms_ok and closed == len(G) * len(G) and inv_ok)
    return ok, res


# ---------------------------------------------------------------- cyclotomic Z[zeta_M]
PHI_CACHE = {}


def poly_divmod_monic(num, den):
    # integer polynomials as coefficient lists (index = degree); den monic
    num = list(num)
    dd = len(den) - 1
    if len(num) - 1 < dd:
        return [0], num
    q = [0] * (len(num) - dd)
    for i in range(len(num) - 1, dd - 1, -1):
        c = num[i]
        if c:
            q[i - dd] = c
            for j in range(dd + 1):
                num[i - dd + j] -= c * den[j]
    rem = num[:dd] if dd > 0 else [0]
    return q, rem


def cyclotomic(M):
    if M in PHI_CACHE:
        return PHI_CACHE[M]
    p = [-1] + [0] * (M - 1) + [1]  # x^M - 1
    for d in range(1, M):
        if M % d == 0:
            q, r = poly_divmod_monic(p, cyclotomic(d))
            if any(r):
                raise ValueError("cyclotomic division")
            p = q
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    PHI_CACHE[M] = p
    return p


def zeta_reduce(coeffs, M):
    # coeffs: list length M (mod x^M - 1) -> canonical remainder mod Phi_M, trailing zeros stripped
    _, r = poly_divmod_monic(coeffs, cyclotomic(M))
    r = list(r)
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return tuple(r)


def rot_add(vec, a, b, M):
    # vec * (zeta^a + zeta^b) in Z[x]/(x^M - 1)
    out = [0] * M
    for i in range(M):
        c = vec[i]
        if c:
            out[(i + a) % M] += c
            out[(i + b) % M] += c
    return out


# monomial 2x2 matrices over mu_M: ('d', a, b) = diag(z^a, z^b); ('a', a, b) = [[0, z^a], [z^b, 0]]
def mono_mul(x, y, M):
    if x[0] == "d" and y[0] == "d":
        return ("d", (x[1] + y[1]) % M, (x[2] + y[2]) % M)
    if x[0] == "d" and y[0] == "a":
        return ("a", (x[1] + y[1]) % M, (x[2] + y[2]) % M)
    if x[0] == "a" and y[0] == "d":
        return ("a", (x[1] + y[2]) % M, (x[2] + y[1]) % M)
    return ("d", (x[1] + y[2]) % M, (x[2] + y[1]) % M)


def mono_dagger(x, M):
    if x[0] == "d":
        return ("d", (-x[1]) % M, (-x[2]) % M)
    return ("a", (-x[2]) % M, (-x[1]) % M)


def mono_det_is_one(x, M):
    if x[0] == "d":
        return (x[1] + x[2]) % M == 0
    # det [[0, z^a], [z^b, 0]] = -z^(a+b) = 1  <=>  z^(a+b) = -1  <=>  M even and a+b = M/2
    return M % 2 == 0 and (x[1] + x[2]) % M == M >> 1


def bfs_group(gens, mul):
    ident = gens[0][1]
    seen = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for _, h in gens[1:]:
                p = mul(g, h)
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
        frontier = nxt
    return sorted(seen)


def build_mono_group(kind, n):
    if kind == "C":
        M = n
        gens = [("I", ("d", 0, 0)), ("g", ("d", 1 % M, (M - 1) % M))]
    else:
        M = 2 * n
        gens = [("I", ("d", 0, 0)), ("a", ("d", 1, M - 1)), ("b", ("a", 0, n))]
    G = bfs_group(gens, lambda x, y: mono_mul(x, y, M))
    return M, G


def check_mono_group(M, G, expect):
    S = set(G)
    closed = sum(1 for x in G for y in G if mono_mul(x, y, M) in S)
    det_ok = all(mono_det_is_one(x, M) for x in G)
    unitary = all(mono_mul(x, mono_dagger(x, M), M) == ("d", 0, 0) for x in G)
    ok = len(S) == expect and closed == len(G) * len(G) and det_ok and unitary and ("d", 0, 0) in S
    return ok, dict(order=len(S), closed_products=closed, pairs=len(G) * len(G), det1=1 if det_ok else 0, unitary=1 if unitary else 0)


def mono_sums(M, G, nmax):
    # character sums sum_g chi_N(g) and moment sums sum_g tr(g)^N for N = 0..nmax, canonical in Z[zeta_M]
    chi_tot = [[0] * M for _ in range(nmax + 1)]
    mom_tot = [[0] * M for _ in range(nmax + 1)]
    for g in G:
        if g[0] == "a":
            # trace 0: chi = 1, 0, -1, 0, 1, ...; tr^0 = 1, tr^N = 0
            for N in range(nmax + 1):
                v = (1, 0, -1, 0)[N % 4]
                chi_tot[N][0] += v
            mom_tot[0][0] += 1
            continue
        a, b = g[1], g[2]
        prev = [0] * M
        cur = [0] * M
        cur[0] = 1
        pw = [0] * M
        pw[0] = 1
        for N in range(nmax + 1):
            for i in range(M):
                chi_tot[N][i] += cur[i]
                mom_tot[N][i] += pw[i]
            t = rot_add(cur, a, b, M)
            nxt = [t[i] - prev[i] for i in range(M)]
            prev, cur = cur, nxt
            pw = rot_add(pw, a, b, M)
    return [zeta_reduce(c, M) for c in chi_tot], [zeta_reduce(c, M) for c in mom_tot]


def const_div(red, order, name):
    # red must be a constant divisible by order
    ok = len(red) == 1 and red[0] % order == 0
    check(ok, name)
    return red[0] // order if ok else None


def ring_sums(R, G, nmax):
    zero, one, add, sub, mul = R[0], R[1], R[2], R[3], R[4]
    chi_tot = [zero] * (nmax + 1)
    mom_tot = [zero] * (nmax + 1)
    for X in G:
        x = X[0]  # trace of the SU(2) image = 2 Re q = doubled real coordinate
        prev, cur = zero, one
        pw = one
        for N in range(nmax + 1):
            chi_tot[N] = add(chi_tot[N], cur)
            mom_tot[N] = add(mom_tot[N], pw)
            prev, cur = cur, sub(mul(x, cur), prev)
            pw = mul(pw, x)
    return chi_tot, mom_tot


def ring_const_div(R, v, order, name):
    if R is RING_Z:
        ok = v % order == 0
        check(ok, name)
        return v // order if ok else None
    ok = v[1] == 0 and v[0] % order == 0
    check(ok, name)
    return v[0] // order if ok else None


def catalan(n):
    c = [1]
    for k in range(n):
        c.append(sum(c[i] * c[k - i] for i in range(k + 1)))
    return c


# ---------------------------------------------------------------- Gegenbauer zonal
def gegen(lam, t, kmax):
    c = [(1, 1), rmul((2 * lam, 1), t)]
    for n in range(2, kmax + 1):
        a = rmul(rmul((2 * (n + lam - 1), 1), t), c[n - 1])
        b = rmul((n + 2 * lam - 2, 1), c[n - 2])
        c.append(rdiv(rsub(a, b), (n, 1)))
    return c


def sk_from_hist(lam, hist, N, kmax):
    # hist: {t (rational): count over ordered pairs}; S_k = (1/N) sum count * C_k(t)/C_k(1)
    c1 = gegen(lam, (1, 1), kmax)
    for k in range(kmax + 1):
        check(c1[k] == (binom(k + 2 * lam - 1, k), 1), "gegen_C_k(1)_lam" + str(lam) + "_k" + str(k))
    vals = {t: gegen(lam, t, kmax) for t in hist}
    S = []
    for k in range(kmax + 1):
        s = (0, 1)
        for t, cnt in hist.items():
            s = radd(s, rmul((cnt, 1), vals[t][k]))
        S.append(rdiv(s, (N * c1[k][0], 1)))
    return S


def design_profile(S, kmax):
    t = 0
    while t + 1 <= kmax and S[t + 1][0] == 0:
        t += 1
    nz = [k for k in range(1, kmax + 1) if S[k][0] != 0]
    neg = [k for k in range(kmax + 1) if S[k][0] < 0]
    return t, nz, neg


# ---------------------------------------------------------------- PART A: SU(2) finite subgroups (the ladder)
def part_ladder():
    groups = []
    for n in range(2, 31):
        groups.append(("C" + str(n), "mono", ("C", n), n))
    for n in range(2, 31):
        groups.append(("Dic" + str(n), "mono", ("Dic", n), 4 * n))
    groups.append(("2T", "quat", (RING_Z, build_2T), 24))
    groups.append(("2O", "quat", (RING_S2, build_2O), 48))
    groups.append(("2I", "quat", (RING_PH, build_2I), 120))
    results = {}
    for name, kind, spec, order in groups:
        if kind == "mono":
            M, G = build_mono_group(spec[0], spec[1])
            ok, info = check_mono_group(M, G, order)
            chi, mom = mono_sums(M, G, NMAX)
            m = [const_div(chi[N], order, name + "_chi_div_" + str(N)) for N in range(NMAX + 1)]
            mo = [const_div(mom[N], order, name + "_mom_div_" + str(N)) for N in range(NMAX + 1)]
            ring = "Z[zeta_" + str(M) + "]_monomial"
        else:
            R, builder = spec
            G = builder()
            ok, info = check_quat_group(R, G, name, order)
            chi, mom = ring_sums(R, G, NMAX)
            m = [ring_const_div(R, chi[N], order, name + "_chi_div_" + str(N)) for N in range(NMAX + 1)]
            mo = [ring_const_div(R, mom[N], order, name + "_mom_div_" + str(N)) for N in range(NMAX + 1)]
            ring = {id(RING_Z): "Z", id(RING_S2): "Z[sqrt2]", id(RING_PH): "Z[phi]"}[id(R)]
        check(ok, name + "_group_checks")
        check(m[0] == 1, name + "_m0")
        t = 0
        while t + 1 <= NMAX and m[t + 1] == 0:
            t += 1
        brk = t + 1
        emit("GROUP", "name", name, "order", info["order"], "expected_order", order, "ring", ring,
             "closed_products", str(info["closed_products"]) + ":" + str(info["pairs"]), "group_ok", 1 if ok else 0, "tag", "MEASURED")
        emit("MOLIEN_M", "name", name, "twoj", "0.." + str(NMAX), "m", ",".join(str(v) for v in m), "tag", "MEASURED")
        emit("DESIGN_T", "name", name, "t", t, "first_break_twoj", brk, "m_at_break", m[brk] if brk <= NMAX else "beyond",
             "tag", "MEASURED")
        emit("MOMENT_M", "name", name, "n", "0.." + str(NMAX), "M", ",".join(str(v) for v in mo), "tag", "MEASURED")
        results[name] = dict(order=order, m=m, M=mo, t=t, brk=brk)
    cat = catalan(NMAX // 2)
    haar = []
    for n in range(NMAX + 1):
        haar.append(cat[n >> 1] if n % 2 == 0 else 0)
    emit("MOMENT_HAAR_CATALAN", "n", "0.." + str(NMAX), "M", ",".join(str(v) for v in haar), "rule", "C_(n/2)_for_even_n_NAMED_Weyl",
         "tag", "MEASURED")
    for name in results:
        mo = results[name]["M"]
        fb = None
        for n in range(NMAX + 1):
            if mo[n] != haar[n]:
                fb = n
                break
        results[name]["moment_break"] = fb
        emit("MOMENT_BREAK_M", "name", name, "first_n", fb, "M", mo[fb] if fb is not None else "none",
             "Mhaar", haar[fb] if fb is not None else "none", "tag", "MEASURED")
    tmax = max(results[n]["t"] for n in results)
    who = [n for n in results if results[n]["t"] == tmax]
    second = max(results[n]["t"] for n in results if results[n]["t"] < tmax)
    emit("CEILING_M", "groups", len(results), "max_t", tmax, "reached_by", ",".join(who), "next_t", second, "scope",
         "C_n,Dic_n_n=2..30_plus_2T,2O,2I", "tag", "MEASURED")
    # required spot checks
    r2i = results["2I"]
    check(r2i["order"] == 120, "2I_order")
    check(r2i["m"][12] == 1, "2I_m(2j=12)=1")
    check(all(r2i["m"][k] == 0 for k in range(1, 12)), "2I_m_zero_below_12")
    check(results["Dic2"]["M"][4] == 4, "Q8_M4=4")
    check(results["2O"]["M"][8] == 15, "2O_M8=15")
    # Molien series of 2I by Klein's invariants (NAMED formula) as a power series, against the measured m
    ser = [0] * (NMAX + 1)
    for a in range(0, NMAX + 1, 12):
        for b in range(0, NMAX + 1 - a, 20):
            ser[a + b] += 1
            if a + b + 30 <= NMAX:
                ser[a + b + 30] += 1
    check(ser == r2i["m"], "2I_molien_klein_series")
    emit("SPOT_LADDER", "order_2I", r2i["order"], "m_2I_twoj12", r2i["m"][12], "m_2I_twoj1to11", ",".join(str(r2i["m"][k]) for k in range(1, 12)),
         "Q8_M4", results["Dic2"]["M"][4], "O2_M8", results["2O"]["M"][8], "Haar_M8", haar[8], "I2_M12", r2i["M"][12], "Haar_M12", haar[12],
         "klein_series_equal", 1 if ser == r2i["m"] else 0, "tag", "MEASURED")
    return results


# ---------------------------------------------------------------- PART B: point sets
def pair_hist(points, ip):
    total = {}
    per = None
    invariant = 1
    for x in points:
        h = {}
        for y in points:
            v = ip(x, y)
            h[v] = h.get(v, 0) + 1
        hs = tuple(sorted(h.items()))
        if per is None:
            per = hs
        elif hs != per:
            invariant = 0
        for v, c in h.items():
            total[v] = total.get(v, 0) + c
    return total, per, invariant


def ph_ip(x, y):
    s = (0, 0)
    for i in range(4):
        s = ph_add(s, ph_mul(x[i], y[i]))
    return s


def ph_to_s5(v, den):
    # (a + b*phi)/den -> rational parts of r + s*sqrt5 with phi = (1 + sqrt5)/2
    a, b = v
    return rq(2 * a + b, 2 * den), rq(b, 2 * den)


def part_600cell(lad):
    P = build_2I()
    N = len(P)
    total, per, inv = pair_hist(P, ph_ip)
    dist = []
    for v, c in sorted(per, key=lambda vc: (-(vc[0][0] * 2 + vc[0][1]), vc[0])):
        r, s = ph_to_s5(v, 4)
        dist.append(rstr(r) + "+" + rstr(s) + "r5#" + str(c))
    check(inv == 1, "600cell_distance_invariance")
    emit("SET_600CELL", "N", N, "distance_invariant_all_base_points", inv, "base_points_checked", N, "cosine_dist", ",".join(dist), "tag", "MEASURED")
    # S_k with lambda = 1 (Chebyshev U), homogenised: D_k(v) = 4^k C_k(v/4), D_0 = 1, D_1 = 2v, D_n = 2v D_(n-1) - 16 D_(n-2)
    Dv = {}
    for v in total:
        d = [(1, 0), (2 * v[0], 2 * v[1])]
        for n in range(2, KMAX + 1):
            d.append(ph_sub(ph_mul((2 * v[0], 2 * v[1]), d[n - 1]), (16 * d[n - 2][0], 16 * d[n - 2][1])))
        Dv[v] = d
    S = []
    rational_all = 1
    for k in range(KMAX + 1):
        acc = (0, 0)
        for v, c in total.items():
            acc = ph_add(acc, (c * Dv[v][k][0], c * Dv[v][k][1]))
        den = N * ipow(4, k) * (k + 1)
        if acc[1] != 0:
            rational_all = 0
        S.append(rq(acc[0], den))
    check(rational_all == 1, "600cell_S_rational")
    t, nz, neg = design_profile(S, KMAX)
    check(not neg, "600cell_S_nonneg")
    check(S[0] == (N, 1), "600cell_S0")
    emit("SK_600CELL", "lambda", 1, "k", "0.." + str(KMAX), "S", ",".join(rstr(s) for s in S), "tag", "MEASURED")
    emit("TPROFILE", "set", "600cell", "N", N, "t", t, "first_break_k", t + 1, "S_break", rstr(S[t + 1]), "nonzero_k", ",".join(str(k) for k in nz),
         "negative", len(neg), "tag", "MEASURED")
    check(S[12] == (120, 13), "600cell_S12=120/13")
    # 30/30 cross: S_k (k+1)/120 == m_k of 2I
    m = lad["2I"]["m"]
    eq = 0
    for k in range(1, KMAX + 1):
        q = rmul(S[k], (k + 1, 120))
        if q == (m[k], 1):
            eq += 1
    check(eq == KMAX, "600cell_cross_30")
    emit("CROSS_600CELL_2I", "rule", "S_k*(k+1)/120==m_k", "k", "1.." + str(KMAX), "equal", str(eq) + ":" + str(KMAX), "tag", "MEASURED")
    return dict(S=S, t=t, nz=nz, dist=per, N=N)


def build_E8():
    out = set()
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (2, -2):
                for sj in (2, -2):
                    v = [0] * 8
                    v[i] = si
                    v[j] = sj
                    out.add(tuple(v))
    for mask in range(256):
        if bin(mask).count("1") % 2 == 0:
            out.add(tuple(-1 if (mask >> i) & 1 else 1 for i in range(8)))
    return sorted(out)


def z_ip(x, y):
    return sum(a * b for a, b in zip(x, y))


def part_E8(drop_one=False):
    P = build_E8()
    if drop_one:
        P = P[1:]
    N = len(P)
    norms = set(z_ip(x, x) for x in P)
    total, per, inv = pair_hist(P, z_ip)
    if drop_one:
        return dict(N=N, invariant=inv)
    check(N == 240 and norms == {8}, "E8_count_norm")
    check(inv == 1, "E8_distance_invariance")
    dist = ",".join(str(v) + "#" + str(c) for v, c in sorted(per, key=lambda vc: -vc[0]))
    emit("SET_E8", "N", N, "norms", ",".join(str(v) for v in sorted(norms)), "distance_invariant_all_base_points", inv, "base_points_checked", N,
         "ip_doubled_dist", dist, "tag", "MEASURED")
    hist = {}
    for v, c in total.items():
        t = rq(v, 8)
        hist[t] = hist.get(t, 0) + c
    S = sk_from_hist(3, hist, N, KMAX)
    t, nz, neg = design_profile(S, KMAX)
    check(not neg, "E8_S_nonneg")
    emit("SK_E8", "lambda", 3, "k", "0.." + str(KMAX), "S", ",".join(rstr(s) for s in S), "tag", "MEASURED")
    emit("TPROFILE", "set", "E8", "N", N, "t", t, "first_break_k", t + 1, "S_break", rstr(S[t + 1]), "nonzero_k", ",".join(str(k) for k in nz),
         "negative", len(neg), "tag", "MEASURED")
    return dict(S=S, t=t, nz=nz, dist=per, N=N)


# ---------------------------------------------------------------- Golay from the icosahedron, Leech from Conway's test
def icosahedron():
    V = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            base = [(0, 0), (s1, 0), (0, s2)]
            for r in range(3):
                V.append(tuple(base[(k - r) % 3] for k in range(3)))
    return V


def golay_icosahedron(flip=False):
    V = icosahedron()
    n = len(V)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                d = (0, 0)
                for k in range(3):
                    e = ph_sub(V[i][k], V[j][k])
                    d = ph_add(d, ph_mul(e, e))
                if d == (4, 0):
                    adj[i][j] = 1
    degs = [sum(r) for r in adj]
    rows = []
    for i in range(n):
        bits = 1 << i
        for j in range(n):
            if i != j and adj[i][j] == 0:
                bits |= 1 << (12 + j)
        bits |= 1 << (12 + i)
        rows.append(bits)
    # row i = e_i | (J - A)_i : (J - A) has 1 on the diagonal and on the 6 non-neighbours
    if flip:
        rows[0] ^= 1 << 17
    code = [0]
    for r in rows:
        code = code + [c ^ r for c in code]
    return degs, sorted(set(code)), len(V)


def weight_enum(code):
    we = {}
    for c in code:
        w = bin(c).count("1")
        we[w] = we.get(w, 0) + 1
    return we


def conway_member(x, CW):
    m = x[0] & 1
    if any((v & 1) != m for v in x):
        return False
    for nu in range(4):
        s = 0
        for i in range(24):
            if x[i] % 4 == nu:
                s |= 1 << i
        if s not in CW:
            return False
    return sum(x) % 8 == (4 * m) % 8


def combos(n, k):
    # all k-subsets of range(n) as bitmasks
    out = []

    def rec(start, left, mask):
        if left == 0:
            out.append(mask)
            return
        for i in range(start, n - left + 1):
            rec(i + 1, left - 1, mask | (1 << i))
    rec(0, k, 0)
    return out


def bits_of(mask):
    return [i for i in range(24) if (mask >> i) & 1]


def build_leech(CWset, CWlist):
    vecs = []
    shape = {}
    # even, shape (+-4, +-4, 0^22)
    c = 0
    for msk in combos(24, 2):
        i, j = bits_of(msk)
        for si in (4, -4):
            for sj in (4, -4):
                v = [0] * 24
                v[i] = si
                v[j] = sj
                if conway_member(v, CWset):
                    vecs.append(tuple(v))
                    c += 1
    shape["4^2"] = c
    # even, shape (+-4, +-2^4, 0^19): the 2-mod-4 set is the 4-support; test every support
    c = 0
    supports_tested = 0
    for i in range(24):
        for msk in combos(23, 4):
            supports_tested += 1
            sup = [b if b < i else b + 1 for b in bits_of(msk)]
            s = 0
            for b in sup:
                s |= 1 << b
            if s in CWset:
                for sg in range(32):
                    v = [0] * 24
                    v[i] = 4 if sg & 16 else -4
                    for q in range(4):
                        v[sup[q]] = 2 if (sg >> q) & 1 else -2
                    if conway_member(v, CWset):
                        vecs.append(tuple(v))
                        c += 1
    shape["4.2^4"] = c
    shape["4.2^4_supports_tested"] = supports_tested
    # even, shape (+-2^8, 0^16): every 8-support tested
    c = 0
    sup8 = 0
    for msk in combos(24, 8):
        sup8 += 1
        if msk in CWset:
            b = bits_of(msk)
            for sg in range(256):
                v = [0] * 24
                for q in range(8):
                    v[b[q]] = -2 if (sg >> q) & 1 else 2
                if conway_member(v, CWset):
                    vecs.append(tuple(v))
                    c += 1
    shape["2^8"] = c
    shape["2^8_supports_tested"] = sup8
    # odd, shape (+-3, +-1^23): for fixed p, sign pattern <-> A = {i : x_i = 1 mod 4} is a bijection, so the passing
    # candidates are exactly those with A a codeword; the full Conway test is still applied to each
    c = 0
    for p in range(24):
        for A in CWlist:
            v = [0] * 24
            for i in range(24):
                inA = (A >> i) & 1
                if i == p:
                    v[i] = -3 if inA else 3
                else:
                    v[i] = 1 if inA else -1
            if conway_member(v, CWset):
                vecs.append(tuple(v))
                c += 1
    shape["3.1^23"] = c
    return vecs, shape


def leech_ip_hist(vecs):
    # all ordered pairs, by packed big integers: field i of T(x) = <x, y_i> + 32, one byte per field
    N = len(vecs)
    cols = []
    for j in range(24):
        pos = bytearray(N)
        neg = bytearray(N)
        for i, y in enumerate(vecs):
            v = y[j]
            if v > 0:
                pos[i] = v
            elif v < 0:
                neg[i] = -v
        cols.append(int.from_bytes(bytes(pos), "little") - int.from_bytes(bytes(neg), "little"))
    ONES32 = int.from_bytes(bytes([32]) * N, "little")
    LIMIT = 1 << (8 * N)
    # subset-sum tables over 4 chunks of 6 coordinates
    tabs = []
    for ch in range(4):
        t = [0] * 64
        for mask in range(1, 64):
            low = (mask & -mask).bit_length() - 1
            t[mask] = t[mask & (mask - 1)] + cols[6 * ch + low]
        tabs.append(t)
    SUMALL = tabs[0][63] + tabs[1][63] + tabs[2][63] + tabs[3][63]
    keys = [0, 16, 24, 32, 40, 48, 64]
    needles = [bytes([k]) for k in keys]
    total = [0] * 7
    per = None
    invariant = 1
    stray = 0
    done = 0
    for x in vecs:
        P = 0
        Nn = 0
        for i in range(24):
            if x[i] > 0:
                P |= 1 << i
            elif x[i] < 0:
                Nn |= 1 << i
        sp = tabs[0][P & 63] + tabs[1][(P >> 6) & 63] + tabs[2][(P >> 12) & 63] + tabs[3][(P >> 18) & 63]
        if x[0] & 1:
            # odd: x = s + (x_p - s_p) e_p with s = sign vector
            acc = 2 * sp - SUMALL
            for i in range(24):
                if x[i] == 3 or x[i] == -3:
                    acc += (x[i] - (1 if x[i] > 0 else -1)) * cols[i]
        else:
            sn = tabs[0][Nn & 63] + tabs[1][(Nn >> 6) & 63] + tabs[2][(Nn >> 12) & 63] + tabs[3][(Nn >> 18) & 63]
            mags = set(abs(v) for v in x if v)
            if len(mags) == 1:
                acc = (sp - sn) * mags.pop()
            else:
                acc = 0
                for i in range(24):
                    if x[i]:
                        acc += x[i] * cols[i]
        T = acc + ONES32
        if T < 0 or T >= LIMIT:
            stray += 1
            continue
        b = T.to_bytes(N, "little")
        h = tuple(b.count(nd) for nd in needles)
        if sum(h) != N:
            stray += 1
        if per is None:
            per = h
        elif h != per:
            invariant = 0
        for q in range(7):
            total[q] += h[q]
        done += 1
        if done % 20000 == 0:
            log("leech pairs: base points done " + str(done))
    ips = [k - 32 for k in keys]
    return ips, per, total, invariant, stray


def part_leech():
    degs, code, nv = golay_icosahedron()
    we = weight_enum(code)
    check(nv == 12 and all(d == 5 for d in degs), "icosahedron_degree5")
    check(len(code) == 4096 and we == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, "golay_weight_enumerator")
    emit("GOLAY", "construction", "[I12,J-A_icosahedron]_A_from_Z[phi]_vertices", "vertices", nv, "degrees", ",".join(str(d) for d in sorted(set(degs))),
         "codewords", len(code), "weight_enum", ",".join(str(w) + ":" + str(c) for w, c in sorted(we.items())), "tag", "MEASURED")
    CWset = set(code)
    vecs, shape = build_leech(CWset, code)
    N = len(vecs)
    distinct = len(set(vecs))
    norms = set(z_ip(v, v) for v in vecs)
    vs = set(vecs)
    antipodal = all(tuple(-a for a in v) in vs for v in vecs)
    check(N == 196560 and distinct == N, "leech_count")
    check(norms == {32}, "leech_norm")
    check(antipodal, "leech_antipodal")
    emit("LEECH", "N", N, "distinct", distinct, "norms", ",".join(str(v) for v in sorted(norms)), "antipodal", 1 if antipodal else 0,
         "shape_4^2", shape["4^2"], "shape_4.2^4", shape["4.2^4"], "supports_4.2^4_tested", shape["4.2^4_supports_tested"],
         "shape_2^8", shape["2^8"], "supports_2^8_tested", shape["2^8_supports_tested"], "shape_3.1^23", shape["3.1^23"],
         "membership", "Conway_mod4_Cset_test_on_every_vector", "tag", "MEASURED")
    ips, per, total, inv, stray = leech_ip_hist(vecs)
    check(stray == 0, "leech_ip_values_in_set")
    check(inv == 1, "leech_distance_invariance")
    dist = ",".join(str(ips[q]) + "#" + str(per[q]) for q in range(6, -1, -1))
    emit("SET_LEECH", "N", N, "ordered_pairs", N * N, "ip_values_outside_set", stray, "distance_invariant_all_base_points", inv,
         "base_points_checked", N, "ip_dist", dist, "ip_value_set", ",".join(str(ips[q]) for q in range(6, -1, -1) if per[q]), "tag", "MEASURED")
    hist = {}
    for q in range(7):
        if total[q]:
            hist[rq(ips[q], 32)] = total[q]
    S = sk_from_hist(11, hist, N, KMAX)
    t, nz, neg = design_profile(S, KMAX)
    check(not neg, "leech_S_nonneg")
    emit("SK_LEECH", "lambda", 11, "k", "0.." + str(KMAX), "S", ",".join(rstr(s) for s in S), "source", "full_ordered_pair_histogram", "tag", "MEASURED")
    emit("TPROFILE", "set", "Leech", "N", N, "t", t, "first_break_k", t + 1, "S_break", rstr(S[t + 1]), "nonzero_k", ",".join(str(k) for k in nz),
         "negative", len(neg), "tag", "MEASURED")
    return dict(S=S, t=t, nz=nz, N=N, CW=CWset, code=code, ips=ips, per=per, antipodal=1 if antipodal else 0)


# ---------------------------------------------------------------- PART C: bridges
def oct_table(flip=False):
    T = {}
    for i in range(8):
        T[(0, i)] = (1, i)
        T[(i, 0)] = (1, i)
    for i in range(1, 8):
        T[(i, i)] = (-1, 0)
    lines = []
    for i in range(1, 8):
        a = i
        b = i % 7 + 1
        c = (i + 2) % 7 + 1
        lines.append((a, b, c))
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            T[(x, y)] = (1, z)
            T[(y, x)] = (-1, z)
    if flip:
        s, k = T[(1, 2)]
        T[(1, 2)] = (-s, k)
    return T, lines


def omul(T, x, y):
    z = [0] * 8
    for i in range(8):
        if x[i]:
            for j in range(8):
                if y[j]:
                    s, k = T[(i, j)]
                    z[k] += s * x[i] * y[j]
    return z


def onorm(x):
    return sum(a * a for a in x)


def oct_units_from_quads(quads):
    out = set()
    for i in range(8):
        for s in (2, -2):
            v = [0] * 8
            v[i] = s
            out.add(tuple(v))
    for q in quads:
        for sg in range(16):
            v = [0] * 8
            for r, idx in enumerate(q):
                v[idx] = -1 if (sg >> r) & 1 else 1
            out.add(tuple(v))
    return out


def oct_closure(T, U):
    sp = [[(i, x[i]) for i in range(8) if x[i]] for x in U]
    outside = 0
    for xs in sp:
        for ys in sp:
            z = [0] * 8
            for i, xi in xs:
                for j, yj in ys:
                    s, k = T[(i, j)]
                    z[k] += s * xi * yj
            if any(c % 2 for c in z):
                outside += 1
                continue
            if tuple(c >> 1 for c in z) not in U:
                outside += 1
    return outside


def rank_int(rows):
    # fraction-free (Bareiss-like) rank over Q of an integer matrix
    A = [list(r) for r in rows]
    nr = len(A)
    nc = len(A[0]) if A else 0
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                p = A[r][c]
                A[i] = [p * A[i][j] - f * A[r][j] for j in range(nc)]
                g = 0
                for v in A[i]:
                    g = gcd(g, v)
                if g > 1:
                    A[i] = [v // g for v in A[i]]
        r += 1
        if r == nr:
            break
    return r


def root_system_check(V, ip, half_norm):
    # V: list of vectors (tuples), ip integer-valued bilinear form with ip(r, r) = 2*half_norm ... here we pass
    # B with B(r,r) = 2; reflection s_r(x) = x - B(x,r) r must land in V
    S = set(V)
    norms = set(ip(v, v) for v in V)
    outside = 0
    for r in V:
        for x in V:
            b = ip(x, r)
            y = half_norm(x, r, b)
            if y not in S:
                outside += 1
    return norms, outside


def part_octonions():
    T, lines = oct_table()
    # algebra checks over deterministic integer pairs
    seed = 12345
    nm_fail = 0
    alt_fail = 0
    for _ in range(2000):
        x = []
        y = []
        for _k in range(8):
            seed = (1103515245 * seed + 12345) % 2147483648
            x.append((seed >> 16) % 11 - 5)
            seed = (1103515245 * seed + 12345) % 2147483648
            y.append((seed >> 16) % 11 - 5)
        if onorm(omul(T, x, y)) != onorm(x) * onorm(y):
            nm_fail += 1
        if omul(T, omul(T, x, x), y) != omul(T, x, omul(T, x, y)) or omul(T, omul(T, y, x), x) != omul(T, y, omul(T, x, x)):
            alt_fail += 1
    check(nm_fail == 0 and alt_fail == 0, "octonion_algebra")
    # Kirmse quads: {0} + line, and complements
    kirmse = []
    for (a, b, c) in lines:
        q = tuple(sorted((0, a, b, c)))
        kirmse.append(q)
        kirmse.append(tuple(sorted(set(range(8)) - set(q))))
    Uk = oct_units_from_quads(kirmse)
    out_k = oct_closure(T, Uk)
    emit("OCT_TABLE", "lines", ".".join("".join(str(v) for v in L) for L in lines), "norm_mult_fail", str(nm_fail) + ":2000",
         "alternativity_fail", str(alt_fail) + ":2000", "tag", "MEASURED")
    emit("OCT_KIRMSE_NAIVE", "units", len(Uk), "products_outside", out_k, "of", len(Uk) * len(Uk), "closed", 1 if out_k == 0 else 0,
         "tag", "MEASURED")
    # every transposition of coordinates applied to the Kirmse quads
    closed_orders = set()
    closed_count = 0
    for a in range(8):
        for b in range(a + 1, 8):
            def sw(i):
                return b if i == a else (a if i == b else i)
            qs = tuple(sorted(tuple(sorted(sw(i) for i in q)) for q in kirmse))
            U = oct_units_from_quads(qs)
            if oct_closure(T, U) == 0:
                closed_count += 1
                closed_orders.add(qs)
    emit("OCT_TRANSPOSED", "transpositions", 28, "closed", closed_count, "distinct_closed_orders", len(closed_orders), "tag", "MEASURED")
    # the order named by the measuring seat
    ms3 = ["0124", "0137", "0156", "0236", "0257", "0345", "0467", "1235", "1267", "1346", "1457", "2347", "2456", "3567"]
    q3 = tuple(sorted(tuple(int(ch) for ch in s) for s in ms3))
    U = oct_units_from_quads(q3)
    out3 = oct_closure(T, U)
    in_transposed = 1 if q3 in closed_orders else 0
    V = sorted(U)
    norms = set(onorm(v) for v in V)
    total, per, inv = pair_hist(V, z_ip)
    # E8 root system: B(x, y) = <X, Y>/2 for doubled X (unit norm 1 -> B = 2); reflection x - B(x, r) r
    def B(x, y):
        s = z_ip(x, y)
        return s >> 1 if s % 2 == 0 else None
    bvals_int = all(B(x, y) is not None for x in V for y in V)
    _, refl_out = root_system_check(V, lambda x, y: B(x, y), lambda x, r, b: tuple(x[i] - b * r[i] for i in range(8)))
    rk = rank_int(V)
    ok = (len(V) == 240 and norms == {4} and out3 == 0 and inv == 1 and bvals_int and refl_out == 0 and rk == 8)
    check(ok, "octonion_order_E8")
    dist = ",".join(str(v) + "#" + str(c) for v, c in sorted(per, key=lambda vc: -vc[0]))
    emit("OCT_ORDER", "quads", ".".join(ms3), "units", len(V), "norm_doubled", ",".join(str(n) for n in sorted(norms)), "products_outside", out3,
         "is_a_transposed_kirmse_order", in_transposed, "distance_invariant", inv, "ip_doubled_dist", dist, "B_integral", 1 if bvals_int else 0,
         "reflections_outside", refl_out, "rank", rk, "E8_by", "240_roots_rank8_reflection_closed_simply_laced_ADE_NAMED", "tag", "MEASURED")
    # control: one table sign flipped
    Tf, _ = oct_table(flip=True)
    out_f = oct_closure(Tf, U)
    nmf = 0
    seed = 999
    for _ in range(500):
        x = []
        y = []
        for _k in range(8):
            seed = (1103515245 * seed + 12345) % 2147483648
            x.append((seed >> 16) % 7 - 3)
            seed = (1103515245 * seed + 12345) % 2147483648
            y.append((seed >> 16) % 7 - 3)
        if onorm(omul(Tf, x, y)) != onorm(x) * onorm(y):
            nmf += 1
    det_ctl = out_f > 0 and nmf > 0
    check(det_ctl, "control_oct_flip_detected")
    emit("CONTROL", "name", "oct_table_one_sign_flipped", "products_outside", out_f, "norm_mult_fail", str(nmf) + ":500", "must_fail", 1,
         "detected", 1 if det_ctl else 0, "tag", "MEASURED")
    check(out_k > 0, "control_kirmse_naive_not_closed")
    emit("CONTROL", "name", "kirmse_naive_closure", "products_outside", out_k, "must_fail", 1, "detected", 1 if out_k > 0 else 0, "tag", "MEASURED")
    return dict(dist=per, N=len(V), out_k=out_k, closed=closed_count, distinct=len(closed_orders))


def part_icosians():
    G = build_2I()
    # 240 = units and phi^-1 * units ; phi^-1 = phi - 1 = (-1, 1)
    inv_phi = (-1, 1)
    V = sorted(set(G) | set(tuple(ph_mul(inv_phi, c) for c in X) for X in G))

    def Bform(lam_on):
        def B(x, y):
            s = ph_ip(x, y)  # <X, Y> = 4 <x, y> = a + b phi ; 2<x,y> = (a + b phi)/2
            if lam_on:
                tr = s[0] + s[1]  # Tr(lambda (a + b phi)) = a + b, lambda = phi/sqrt5
            else:
                tr = 2 * s[0] + s[1]  # Tr(a + b phi) = 2a + b
            return tr >> 1 if tr % 2 == 0 else None
        return B
    B = Bform(True)
    vals = [[B(x, y) for y in V] for x in V]
    integral = all(v is not None for row in vals for v in row)
    norms = set(vals[i][i] for i in range(len(V)))
    per0 = None
    inv = 1
    for i in range(len(V)):
        h = {}
        for v in vals[i]:
            h[v] = h.get(v, 0) + 1
        hs = tuple(sorted(h.items(), key=lambda kv: -kv[0] if kv[0] is not None else 99))
        if per0 is None:
            per0 = hs
        elif hs != per0:
            inv = 0
    S = set(V)
    refl_out = 0
    for j in range(len(V)):
        r = V[j]
        for i in range(len(V)):
            b = vals[i][j]
            x = V[i]
            y = tuple(ph_sub(x[k], (b * r[k][0], b * r[k][1])) for k in range(4))
            if y not in S:
                refl_out += 1
    flat = [tuple(c for pair in v for c in pair) for v in V]
    rk = rank_int(flat)
    ok = len(V) == 240 and integral and norms == {2} and inv == 1 and refl_out == 0 and rk == 8
    check(ok, "icosian_E8")
    dist = ",".join(str(v) + "#" + str(c) for v, c in per0)
    emit("ICOSIAN_E8", "vectors", len(V), "form", "B(x,y)=Tr(lambda*2<x,y>)_lambda=phi/sqrt5", "B_integral", 1 if integral else 0,
         "norms_B", ",".join(str(n) for n in sorted(x for x in norms if x is not None)), "distance_invariant", inv, "B_dist", dist,
         "reflections_outside", refl_out, "rank_over_Q", rk, "E8_by", "240_roots_rank8_reflection_closed_ADE_NAMED", "tag", "MEASURED")
    # control: plain trace form (lambda = 1)
    B1 = Bform(False)
    norms1 = set(B1(x, x) for x in V)
    det_ctl = norms1 != {2}
    check(det_ctl, "control_icosian_plain_trace_detected")
    emit("CONTROL", "name", "icosian_lambda_1_plain_trace", "norms_B", ",".join(str(n) for n in sorted(x for x in norms1 if x is not None)),
         "uniform_norm_2", 0 if det_ctl else 1, "must_fail", 1, "detected", 1 if det_ctl else 0, "tag", "MEASURED")
    return dict(dist=per0, ok=ok, N=len(V))


# Z[zeta8] as 4-tuples mod x^4 + 1
def z8_add(x, y):
    return tuple(x[i] + y[i] for i in range(4))


def z8_sub(x, y):
    return tuple(x[i] - y[i] for i in range(4))


def z8_mul(x, y):
    r = [0] * 8
    for i in range(4):
        if x[i]:
            for j in range(4):
                r[i + j] += x[i] * y[j]
    return tuple(r[i] - r[i + 4] for i in range(4))


def z8_conj(x):
    # zeta^k -> zeta^-k = -zeta^(4-k) for k = 1..3
    return (x[0], -x[3], -x[2], -x[1])


def z8_zpow(k):
    k %= 8
    v = [0, 0, 0, 0]
    if k < 4:
        v[k] = 1
    else:
        v[k - 4] = -1
    return tuple(v)


Z8_0 = (0, 0, 0, 0)
Z8_1 = (1, 0, 0, 0)
Z8_I = (0, 0, 1, 0)
Z8_S2 = (0, 1, 0, -1)  # sqrt2 = zeta - zeta^3


def m2_mul_half(A, B):
    # A, B are 2U, 2V ; return 2UV = (A B)/2
    out = []
    for i in range(2):
        for j in range(2):
            s = z8_add(z8_mul(A[2 * i], B[j]), z8_mul(A[2 * i + 1], B[2 + j]))
            if any(c % 2 for c in s):
                return None
            out.append(tuple(c >> 1 for c in s))
    return tuple(out)


def m2_det_quarter(A):
    # det(2U) = 4 det U
    d = z8_sub(z8_mul(A[0], A[3]), z8_mul(A[1], A[2]))
    if any(c % 4 for c in d):
        return None
    return tuple(c >> 2 for c in d)


def scal(s, A):
    return tuple(z8_mul(s, a) for a in A)


def part_clifford(S_gen=None, tag_name="clifford"):
    H2 = (Z8_S2, Z8_S2, Z8_S2, tuple(-c for c in Z8_S2))  # 2H = sqrt2 [[1, 1], [1, -1]]
    if S_gen is None:
        S2 = ((2, 0, 0, 0), Z8_0, Z8_0, (0, 0, 2, 0))  # 2S = diag(2, 2i)
    else:
        S2 = S_gen
    I2 = ((2, 0, 0, 0), Z8_0, Z8_0, (2, 0, 0, 0))
    seen = {I2}
    frontier = [I2]
    bad = 0
    while frontier:
        nxt = []
        for g in frontier:
            for h in (H2, S2):
                p = m2_mul_half(g, h)
                if p is None:
                    bad += 1
                    continue
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
        frontier = nxt
        if len(seen) > 5000:
            break
    U2 = sorted(seen)
    # normalise: det U = zeta^(2m) -> +- zeta^(-m) U in SU(2)
    norm_set = set()
    det_exps = set()
    nfail = 0
    for A in U2:
        d = m2_det_quarter(A)
        e = None
        for k in range(8):
            if d == z8_zpow(k):
                e = k
        if e is None or e % 2:
            nfail += 1
            continue
        det_exps.add(e)
        c = z8_zpow(-(e >> 1))
        for s in (c, tuple(-v for v in c)):
            B = scal(s, A)
            norm_set.add(B)
    # moments E |tr U|^(2t) over U2 and over the normalised set
    def moments(Gs, tmax):
        out = []
        for t in range(1, tmax + 1):
            acc = Z8_0
            for A in Gs:
                tr = z8_add(A[0], A[3])  # tr(2U) = 2 tr U
                a2 = z8_mul(tr, z8_conj(tr))
                p = Z8_1
                for _ in range(t):
                    p = z8_mul(p, a2)
                acc = z8_add(acc, p)
            den = ipow(4, t) * len(Gs)
            if acc[1] or acc[2] or acc[3] or acc[0] % den:
                out.append(None)
            else:
                out.append(acc[0] // den)
        return out
    return U2, sorted(norm_set), det_exps, nfail, bad, moments


def quat_to_m2(X):
    # 2q over Z[sqrt2] -> 2U over Z[zeta8]: q = a + b i + c j + d k -> [[a + b i, c + d i], [-c + d i, a - b i]]
    def emb(v):
        return z8_add((v[0], 0, 0, 0), tuple(v[1] * c for c in Z8_S2))
    a, b, c, d = (emb(v) for v in X)
    return (z8_add(a, z8_mul(Z8_I, b)), z8_add(c, z8_mul(Z8_I, d)), z8_add(tuple(-v for v in c), z8_mul(Z8_I, d)), z8_sub(a, z8_mul(Z8_I, b)))


def part_clifford_all(lad):
    U2, NS, det_exps, nfail, bad, moments = part_clifford()
    O = build_2O()
    OM = [quat_to_m2(X) for X in O]
    OMset = set(OM)
    hom_fail = 0
    for i, X in enumerate(O):
        for j, Y in enumerate(O):
            P = dq_mul(RING_S2, X, Y)
            if P is None or m2_mul_half(OM[i], OM[j]) != quat_to_m2(P):
                hom_fail += 1
    eq = 1 if set(NS) == OMset else 0
    mU = moments(U2, 5)
    mN = moments(NS, 5)
    cat = catalan(5)
    haar = cat[1:6]
    fd = None
    for t in range(5):
        if mN[t] != haar[t]:
            fd = t + 1
            break
    ok = len(U2) == 192 and len(NS) == 48 and nfail == 0 and bad == 0 and eq == 1 and hom_fail == 0 and mN == mU and fd == 4
    check(ok, "clifford_equals_2O")
    check(mN[3] == lad["2O"]["M"][8] and mN[4] == lad["2O"]["M"][10], "clifford_moments_vs_my_2O_ladder")
    emit("CLIFFORD", "U2_order", len(U2), "det_exponents_zeta8", ",".join(str(e) for e in sorted(det_exps)), "normalised_SU2", len(NS),
         "normalisation_failures", nfail, "equals_2O_quaternion_set", eq, "quat_to_SU2_hom_failures", str(hom_fail) + ":2304",
         "E_abs_tr^2t_U2", ",".join(str(v) for v in mU), "E_abs_tr^2t_SU2", ",".join(str(v) for v in mN), "Haar_Catalan", ",".join(str(v) for v in haar),
         "first_t_differs", fd, "tag", "MEASURED")
    # control: S replaced by diag(1, -1)
    U2z, NSz, _, nfz, badz, momz = part_clifford(S_gen=((2, 0, 0, 0), Z8_0, Z8_0, (-2, 0, 0, 0)))
    mz = momz(NSz, 5)
    det_ctl = set(NSz) != OMset
    check(det_ctl, "control_clifford_Z_detected")
    emit("CONTROL", "name", "clifford_S_replaced_by_Z", "U2_order", len(U2z), "normalised", len(NSz), "equals_2O", 0 if det_ctl else 1,
         "E_abs_tr^2t", ",".join(str(v) for v in mz), "must_fail", 1, "detected", 1 if det_ctl else 0, "tag", "MEASURED")
    return dict(mN=mN, fd=fd, eq=eq, order=len(NS))


# ---------------------------------------------------------------- controls on the point sets / groups
def controls(lad):
    # 2I minus one element: closure and order must fail
    G = build_2I()[1:]
    ok, info = check_quat_group(RING_PH, G, "2I-1", 120)
    check(not ok, "control_2I_drop_one_detected")
    emit("CONTROL", "name", "2I_drop_one_element", "order", info["order"], "closed_products", str(info["closed_products"]) + ":" + str(info["pairs"]),
         "group_ok", 1 if ok else 0, "must_fail", 1, "detected", 0 if ok else 1, "tag", "MEASURED")
    # Golay with one generator bit flipped
    _, code_f, _ = golay_icosahedron(flip=True)
    we = weight_enum(code_f)
    det = we != {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    check(det, "control_golay_flip_detected")
    emit("CONTROL", "name", "golay_one_bit_flipped", "weight_enum", ",".join(str(w) + ":" + str(c) for w, c in sorted(we.items())), "must_fail", 1,
         "detected", 1 if det else 0, "tag", "MEASURED")
    # E8 minus one root: distance invariance must fail
    r = part_E8(drop_one=True)
    det = r["invariant"] == 0
    check(det, "control_E8_drop_one_detected")
    emit("CONTROL", "name", "E8_drop_one_root", "N", r["N"], "distance_invariant", r["invariant"], "must_fail", 1, "detected", 1 if det else 0,
         "tag", "MEASURED")
    # wrong zonal dimension: lambda = 1/2 (S^2 Legendre) on the 600-cell is not available with integer lam in gegen();
    # use lambda = 2 (S^5 zonals) on the E8 histogram instead: t must change from 7
    P = build_E8()
    total, _, _ = pair_hist(P, z_ip)
    hist = {}
    for v, c in total.items():
        hist[rq(v, 8)] = hist.get(rq(v, 8), 0) + c
    S_wrong = sk_from_hist(2, hist, 240, 12)
    t_wrong, _, neg_w = design_profile(S_wrong, 12)
    det = t_wrong != 7 or len(neg_w) > 0
    check(det, "control_wrong_lambda_detected")
    emit("CONTROL", "name", "E8_zonal_wrong_lambda_2", "t", t_wrong, "S2", rstr(S_wrong[2]), "negative_S_count", len(neg_w), "must_fail", 1,
         "detected", 1 if det else 0, "tag", "MEASURED")
    # Leech membership of an octad vector with one minus sign
    return None


def leech_sign_control(CWset, code):
    octad = None
    for c in code:
        if bin(c).count("1") == 8:
            octad = c
            break
    v = [0] * 24
    first = True
    for i in range(24):
        if (octad >> i) & 1:
            v[i] = -2 if first else 2
            first = False
    mem = conway_member(v, CWset)
    check(not mem, "control_leech_odd_sign_detected")
    emit("CONTROL", "name", "leech_octad_one_minus_sign", "member", 1 if mem else 0, "norm", z_ip(v, v), "must_fail", 1, "detected", 0 if mem else 1,
         "tag", "MEASURED")


# ---------------------------------------------------------------- COMPARE (phase 2; added after phase 1 wrote the measured rows)
FROZEN_ROOT = "/mnt/c/tmp/opus55/iconic-3-space/repo/predictions"
FROZEN_SPECS = (("su2", "su2-ladder/PREDICTION.hbp", "c2601cf08d7c662422e7284cc123250c27cf4943766d4d41baa43b63e3b8c865"),
                ("exc", "exceptional/PREDICTION-EXCEPTIONAL.hbp", "e8403cd4c8058f15e4c78ac6d5b8730901cdeb8bcf30f86a0191c9d74fdac24a"))
CMP = [0, 0, []]


def parse_frozen(text):
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


def read_frozen(rel, want):
    path = os.path.join(FROZEN_ROOT, rel)
    try:
        with open(path, "rb") as f:
            data = f.read()
        with open(path + ".sha256", "rb") as f:
            side = f.read().decode("utf-8").split()
    except OSError:
        emit("FROZEN_FILE", "file", rel, "missing", 1, "ok", 0, "tag", "MEASURED")
        check(False, "frozen_missing:" + rel)
        return None
    h = hashlib.sha256(data).hexdigest()
    ok = h == want and len(side) >= 2 and side[0] == h and side[1] == rel.split("/")[-1]
    emit("FROZEN_FILE", "file", rel, "sha256", h, "sidecar_equal", 1 if side and side[0] == h else 0, "pinned_equal", 1 if h == want else 0,
         "bytes", len(data), "ok", 1 if ok else 0, "tag", "MEASURED")
    check(ok, "frozen_sha:" + rel)
    return parse_frozen(data.decode("utf-8")) if ok else None


def by_key(rows, tag, key):
    out = {}
    for t, d in rows:
        if t == tag and key in d:
            out[d[key]] = d
    return out


def cmpv(src, row, item, quantity, frozen, measured, record=True):
    agree = 1 if str(frozen) == str(measured) else 0
    if record:
        CMP[0] += 1
        if agree:
            CMP[1] += 1
        else:
            CMP[2].append(src + ":" + row + ":" + item + ":" + quantity)
        emit("COMPARE", "src", src, "row", row, "item", item, "quantity", quantity, "frozen", frozen, "measured", measured, "agree", agree,
             "tag", "MEASURED_vs_PREDICTION")
    return agree


def lst(v):
    return ",".join(str(x) for x in v)


def head_list(s, n):
    return ",".join(s.split(",")[:n])


def cos_s5(num, den):
    return rstr(rq(num, den)) + "+0:1r5"


def phase1_equal(pre_rows):
    # the measured block of this run against the phase-1 file written before any frozen file was opened; the only
    # permitted differences are the GOLAY label fix and the assert count it removed
    path = os.path.join(HERE, "REFUTE-MEASURE.phase1.hbp")
    try:
        with open(path, "rb") as f:
            old = f.read().decode("utf-8")
    except OSError:
        emit("PHASE1", "file", "REFUTE-MEASURE.phase1.hbp", "missing", 1, "ok", 0, "tag", "MEASURED")
        check(False, "phase1_missing")
        return
    fixed = old.replace("[I12/J-A_icosahedron]", "[I12,J-A_icosahedron]").replace(
        "assert_failures=1|failed=row_format:GOLAY:construction", "assert_failures=0|failed=none")
    now = "\n".join(pre_rows) + "\n"
    same = 1 if fixed == now else 0
    emit("PHASE1", "file", "REFUTE-MEASURE.phase1.hbp", "phase1_sha256", hashlib.sha256(old.encode("utf-8")).hexdigest(),
         "substitutions", "GOLAY_label_pipe_to_comma,MEASURED_BLOCK_END_assert_count", "measured_block_equal_after_substitution", same,
         "tag", "MEASURED")
    check(same == 1, "phase1_block_changed")


def compare_stage(M, pre_rows):
    phase1_equal(pre_rows)
    su2 = read_frozen(FROZEN_SPECS[0][1], FROZEN_SPECS[0][2])
    exc = read_frozen(FROZEN_SPECS[1][1], FROZEN_SPECS[1][2])
    if su2 is None or exc is None:
        return 0
    lad = M["lad"]

    def pname(n):
        if n.startswith("Dic"):
            return "Dic_" + n[3:]
        if n.startswith("C"):
            return "C_" + n[1:]
        return n
    mol = by_key(su2, "MOLIEN", "group")
    des = by_key(su2, "DESIGN", "group")
    mom = by_key(su2, "MOMENT", "group")
    mbr = by_key(su2, "MOMENT_BREAK", "group")
    upto = NMAX + 1
    for name, r in lad.items():
        g = pname(name)
        if g not in mol or g not in des or g not in mom or g not in mbr:
            cmpv("su2", "ALL", g, "present", 1, 0)
            continue
        cmpv("su2", "MOLIEN", g, "order", mol[g].get("order"), r["order"])
        cmpv("su2", "MOLIEN", g, "a[0.." + str(NMAX) + "]", head_list(mol[g].get("a", ""), upto), lst(r["m"]))
        cmpv("su2", "DESIGN", g, "t_pred", des[g].get("t_pred"), r["t"])
        cmpv("su2", "DESIGN", g, "first_break_degree", des[g].get("first_break_degree"), r["brk"])
        cmpv("su2", "DESIGN", g, "m_at_break", des[g].get("m_at_break"), r["m"][r["brk"]])
        cmpv("su2", "MOMENT", g, "M[0.." + str(NMAX) + "]", head_list(mom[g].get("M", ""), upto), lst(r["M"]))
        fb = r["moment_break"]
        cmpv("su2", "MOMENT_BREAK", g, "first_n", mbr[g].get("first_n_with_M_ne_Mh"), fb)
        cmpv("su2", "MOMENT_BREAK", g, "M", mbr[g].get("M"), r["M"][fb])
    cat = catalan(NMAX // 2)
    haar = [cat[n >> 1] if n % 2 == 0 else 0 for n in range(NMAX + 1)]
    mh = [d for t, d in su2 if t == "MOMENT_HAAR"]
    cmpv("su2", "MOMENT_HAAR", "SU2", "Mh[0.." + str(NMAX) + "]", head_list(mh[0].get("Mh", "") if mh else "", upto), lst(haar))
    ce = [d for t, d in su2 if t == "CEILING"]
    tmax = max(lad[n]["t"] for n in lad)
    who = [n for n in lad if lad[n]["t"] == tmax]
    nxt = max(lad[n]["t"] for n in lad if lad[n]["t"] < tmax)
    nxt_who = [n for n in lad if lad[n]["t"] == nxt]
    c0 = ce[0] if ce else {}
    cmpv("su2", "CEILING", "all", "max_t_pred", c0.get("max_t_pred"), tmax)
    cmpv("su2", "CEILING", "all", "attained_by", c0.get("attained_by"), ",".join(who))
    cmpv("su2", "CEILING", "all", "count_attaining", c0.get("count_attaining"), len(who))
    cmpv("su2", "CEILING", "all", "next_highest", c0.get("next_highest"), str(nxt) + "(" + ",".join(nxt_who) + ")")
    cmpv("su2", "CEILING", "all", "groups_listed", c0.get("groups_listed"), len(lad))
    # exceptional
    setmap = (("600cell", "600cell_2I_icosians"), ("E8", "E8_roots"), ("Leech", "Leech_min_vectors"))
    meas = {"600cell": M["c600"], "E8": M["e8"], "Leech": M["leech"]}
    sk = {}
    for t, d in exc:
        if t == "SK":
            sk[(d.get("set"), d.get("k"))] = d
    tp = by_key(exc, "TPRED", "set")
    fbk = by_key(exc, "FIRST_BREAK", "set")
    prof = by_key(exc, "PROFILE", "set")
    sets = by_key(exc, "SET", "set")
    # measured distance distributions as cosine#count multisets
    dists = {}
    dists["600cell"] = sorted(rstr(ph_to_s5(v, 4)[0]) + "+" + rstr(ph_to_s5(v, 4)[1]) + "r5#" + str(c) for v, c in M["c600"]["dist"])
    dists["E8"] = sorted(cos_s5(v, 8) + "#" + str(c) for v, c in M["e8"]["dist"])
    le = M["leech"]
    dists["Leech"] = sorted(cos_s5(le["ips"][q], 32) + "#" + str(le["per"][q]) for q in range(7) if le["per"][q])
    for ms, fs in setmap:
        S = meas[ms]["S"]
        N = meas[ms]["N"]
        for k in range(1, KMAX + 1):
            d = sk.get((fs, str(k)), {})
            cmpv("exc", "SK", fs, "S_a[k=" + str(k) + "]", d.get("S_a"), rstr(S[k]))
            cmpv("exc", "SK", fs, "S_b[k=" + str(k) + "]", d.get("S_b"), "0:1")
        t = meas[ms]["t"]
        cmpv("exc", "TPRED", fs, "t_pred", tp.get(fs, {}).get("t_pred"), t)
        cmpv("exc", "TPRED", fs, "N", tp.get(fs, {}).get("N"), N)
        f = fbk.get(fs, {})
        cmpv("exc", "FIRST_BREAK", fs, "k", f.get("k"), t + 1)
        cmpv("exc", "FIRST_BREAK", fs, "S_a", f.get("S_a"), rstr(S[t + 1]))
        cmpv("exc", "FIRST_BREAK", fs, "pair_sum_a", f.get("pair_sum_a"), rstr(rmul(S[t + 1], (N, 1))))
        p = prof.get(fs, {})
        nz = meas[ms]["nz"]
        cmpv("exc", "PROFILE", fs, "nonzero_k", p.get("nonzero_k"), lst(nz))
        cmpv("exc", "PROFILE", fs, "zero_k", p.get("zero_k"), lst([k for k in range(1, KMAX + 1) if k not in nz]))
        cmpv("exc", "PROFILE", fs, "nonzero_count", p.get("nonzero_count"), len(nz))
        s = sets.get(fs, {})
        cmpv("exc", "SET", fs, "N", s.get("N"), N)
        cmpv("exc", "SET", fs, "rows", s.get("rows"), len(dists[ms]))
        cmpv("exc", "SET", fs, "counts_sum", s.get("counts_sum"), sum(int(x.split("#")[1]) for x in dists[ms]))
        cmpv("exc", "SET", fs, "dist_multiset", ",".join(sorted(s.get("dist", "").split(","))), ",".join(dists[ms]))
    cmpv("exc", "SET", "Leech_min_vectors", "antipodal", sets.get("Leech_min_vectors", {}).get("antipodal"), le["antipodal"])
    cr = [d for t, d in exc if t == "CROSS"]
    c = cr[0] if cr else {}
    m2i = lad["2I"]["m"]
    first_inv = min(k for k in range(1, NMAX + 1) if m2i[k] > 0)
    eq = 0
    for k in range(1, KMAX + 1):
        if rmul(M["c600"]["S"][k], (k + 1, 120)) == (m2i[k], 1):
            eq += 1
    cmpv("exc", "CROSS", "600cell_2I_icosians", "pred_first_break_k", c.get("pred_first_break_k"), M["c600"]["t"] + 1)
    cmpv("exc", "CROSS", "600cell_2I_icosians", "su2_ladder_first_invariant_degree", c.get("su2_ladder_first_invariant_degree"), lad["2I"]["brk"])
    cmpv("exc", "CROSS", "600cell_2I_icosians", "molien_first_invariant_degree", c.get("molien_first_invariant_degree"), first_inv)
    cmpv("exc", "CROSS", "600cell_2I_icosians", "per_degree_exact", c.get("per_degree_exact"), str(eq) + ":" + str(KMAX))
    # bridges against the frozen rows
    o2 = mom.get("2O", {}).get("M", "").split(",")
    cl = M["cl"]
    for t in range(1, 6):
        cmpv("bridge", "MOMENT", "2O_from_Clifford", "M[" + str(2 * t) + "]", o2[2 * t] if len(o2) > 2 * t else None, cl["mN"][t - 1])
    cmpv("bridge", "MOMENT_BREAK", "2O_from_Clifford", "first_n", mbr.get("2O", {}).get("first_n_with_M_ne_Mh"), 2 * cl["fd"])
    cmpv("bridge", "MOLIEN", "2O_from_Clifford", "order", mol.get("2O", {}).get("order"), cl["order"])
    e8f = ",".join(sorted(sets.get("E8_roots", {}).get("dist", "").split(",")))
    od = ",".join(sorted(cos_s5(v, 4) + "#" + str(c) for v, c in M["octo"]["dist"]))
    idd = ",".join(sorted(cos_s5(v, 2) + "#" + str(c) for v, c in M["ico"]["dist"]))
    cmpv("bridge", "SET", "E8_from_octonion_order", "dist_multiset", e8f, od)
    cmpv("bridge", "SET", "E8_from_icosians", "dist_multiset", e8f, idd)
    # control: perturbed frozen values must be caught
    t2 = des.get("2I", {}).get("t_pred", "0")
    c1 = cmpv("control", "DESIGN", "2I", "t_pred_plus_1", str(int(t2) + 1), lad["2I"]["t"], record=False)
    c2 = cmpv("control", "SK", "600cell", "S12_numerator_plus_1", "121:13", rstr(M["c600"]["S"][12]), record=False)
    half = [x for x in dists["Leech"] if x.startswith("1:2+0:1r5#")]
    c3 = cmpv("control", "SET", "Leech", "one_count_moved", "1:2+0:1r5#4601", half[0] if half else "", record=False)
    det = 1 if (c1 == 0 and c2 == 0 and c3 == 0) else 0
    check(det == 1, "control_compare_perturbed_not_detected")
    emit("CONTROL", "name", "compare_perturbed_frozen_values", "t_plus_1_agree", c1, "S12_plus_agree", c2, "leech_count_moved_agree", c3,
         "must_fail", 1, "detected", det, "tag", "MEASURED")
    emit("COMPARE_SUMMARY", "compared", CMP[0], "agree", CMP[1], "disagree", CMP[0] - CMP[1], "which", ",".join(CMP[2]) if CMP[2] else "none",
         "rule", "a_disagreement_is_a_result_not_edited", "tag", "MEASURED_vs_PREDICTION")
    return CMP[0] - CMP[1]


# ---------------------------------------------------------------- write
def write_file(name, rows):
    path = os.path.join(HERE, name)
    data = ("\n".join(rows) + "\n").encode("utf-8")
    with open(path, "wb") as f:
        f.write(data)
    h = hashlib.sha256(data).hexdigest()
    with open(path + ".sha256", "wb") as f:
        f.write((h + "  " + name + "\n").encode("utf-8"))
    return h


def measure():
    emit("REFUTE_MEASURE", "seat", "second_seat_refuter", "repo", "JesseBrown1980/Iconic-3-space", "main", "8f355db", "E", 0, "network", 0,
         "float_used", 0, "kmax", KMAX, "nmax", NMAX, "independence", "predict_scripts_not_read_measure_scripts_not_imported", "tag", "MEASURED")
    log("ladder")
    lad = part_ladder()
    log("600cell")
    c600 = part_600cell(lad)
    log("E8")
    e8 = part_E8()
    log("octonions")
    octo = part_octonions()
    log("icosians")
    ico = part_icosians()
    log("clifford")
    cl = part_clifford_all(lad)
    log("controls")
    controls(lad)
    log("leech")
    le = part_leech()
    leech_sign_control(le["CW"], le["code"])
    # cross: the octonion and icosian E8 distributions against the coordinate E8
    e8d = tuple((v, c) for v, c in e8["dist"])
    oct_counts = tuple(c for v, c in sorted(octo["dist"], key=lambda vc: -vc[0]))
    e8_counts = tuple(c for v, c in sorted(e8d, key=lambda vc: -vc[0]))
    ico_counts = tuple(c for v, c in ico["dist"])
    same = oct_counts == e8_counts == ico_counts
    check(same, "e8_three_constructions_same_distribution")
    emit("E8_THREE_WAYS", "coordinate", ",".join(str(c) for c in e8_counts), "octonion", ",".join(str(c) for c in oct_counts),
         "icosian", ",".join(str(c) for c in ico_counts), "equal", 1 if same else 0, "tag", "MEASURED")
    return dict(lad=lad, c600=c600, e8=e8, leech=le, octo=octo, cl=cl, ico=ico)


def main():
    measure_only = "--measure-only" in sys.argv
    M = measure()
    emit("MEASURED_BLOCK_END", "rows", len(ROWS) + 1, "assert_failures", len(FAILS), "failed", ",".join(FAILS) if FAILS else "none", "tag", "MEASURED")
    pre_rows = list(ROWS)
    h = write_file(PRE_NAME, pre_rows)
    log("wrote precompare " + h)
    if measure_only:
        sys.exit(1 if FAILS else 0)
    dis = compare_stage(M, pre_rows)
    emit("SUMMARY", "assert_failures", len(FAILS), "failed", ",".join(FAILS) if FAILS else "none", "compare_disagree", dis,
         "status", "PASS" if not FAILS and dis == 0 else ("FAIL" if FAILS else "DISAGREE"), "tag", "MEASURED")
    h2 = write_file(OUT_NAME, ROWS)
    log("wrote " + OUT_NAME + " " + h2)
    if FAILS:
        sys.exit(1)
    if dis:
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
