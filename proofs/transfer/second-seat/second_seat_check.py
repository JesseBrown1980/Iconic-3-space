#!/usr/bin/env python3
"""Second-seat independent exact check of TRANSFER-THEOREM.md.

Exact arithmetic only: fractions.Fraction, integers, and Q(sqrt d) pairs.
No float anywhere.  Written from scratch; it does not import or read the
prover's checker.  Groups are built by quaternion closure, Haar moments by
Weyl integration (constant-term extraction, not the Clebsch-Gordan
recursion), designs by the DGS monomial definition directly.

Usage: python3 -I second_seat_check.py OUTDIR
"""
import hashlib
import itertools
import os
import re
import sys
from fractions import Fraction as Fr
from math import factorial, comb

ROWS = []


def row(tag, **kv):
    parts = [tag] + ["%s=%s" % (k, v) for k, v in kv.items()] + ["json=0"]
    ROWS.append("|".join(parts))


# ---------------------------------------------------------------- Q(sqrt d)
class QF:
    __slots__ = ("a", "b", "d")

    def __init__(self, a, b, d):
        self.a = Fr(a)
        self.b = Fr(b)
        self.d = d

    def __add__(s, o):
        return QF(s.a + o.a, s.b + o.b, s.d)

    def __sub__(s, o):
        return QF(s.a - o.a, s.b - o.b, s.d)

    def __neg__(s):
        return QF(-s.a, -s.b, s.d)

    def __mul__(s, o):
        return QF(s.a * o.a + s.d * s.b * o.b, s.a * o.b + s.b * o.a, s.d)

    def key(s):
        return (s.a, s.b)

    def __eq__(s, o):
        return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((s.a, s.b))


def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qkey(q):
    return tuple(x.key() for x in q)


def closure(gens, d):
    one = (QF(1, 0, d), QF(0, 0, d), QF(0, 0, d), QF(0, 0, d))
    elems = {qkey(one): one}
    frontier = [one]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = qmul(x, g)
                k = qkey(y)
                if k not in elems:
                    elems[k] = y
                    nxt.append(y)
        frontier = nxt
    return list(elems.values())


def build_groups():
    G = {}
    h = Fr(1, 2)
    # Q(sqrt3): Z6, Q8, Dic3, 2T
    d = 3
    z = QF(0, 0, d)
    on = QF(1, 0, d)
    hf = QF(h, 0, d)
    s3h = QF(0, h, d)  # sqrt3/2
    qi = (z, on, z, z)
    qj = (z, z, on, z)
    a6 = (hf, s3h, z, z)  # e^{i pi/3}
    G["Z6"] = (closure([a6], d), d)
    G["Q8"] = (closure([qi, qj], d), d)
    G["Dic3"] = (closure([a6, qj], d), d)
    G["2T"] = (closure([qi, qj, (hf, hf, hf, hf)], d), d)
    # Q(sqrt2): Dic4, 2O
    d = 2
    z = QF(0, 0, d)
    on = QF(1, 0, d)
    hf = QF(h, 0, d)
    r2h = QF(0, h, d)  # sqrt2/2
    qi = (z, on, z, z)
    qj = (z, z, on, z)
    a8 = (r2h, r2h, z, z)
    G["Dic4"] = (closure([a8, qj], d), d)
    G["2O"] = (closure([qi, qj, (hf, hf, hf, hf), a8], d), d)
    # Q(sqrt5): 2I
    d = 5
    z = QF(0, 0, d)
    hf = QF(h, 0, d)
    phi_h = QF(Fr(1, 4), Fr(1, 4), d)       # phi/2 = (1+sqrt5)/4
    iphi_h = QF(Fr(-1, 4), Fr(1, 4), d)     # (1/phi)/2 = (sqrt5-1)/4
    G["2I"] = (closure([(hf, hf, hf, hf), (phi_h, iphi_h, hf, z)], d), d)
    return G


# ---------------------------------------------------------------- Haar via Weyl
def laurent_mul(p, q):
    r = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            r[e1 + e2] = r.get(e1 + e2, 0) + c1 * c2
    return r


def haar_mult(n, twoj):
    """E_Haar[chi_{1/2}^n chi_j] via Weyl: (1/2) CT[(x+1/x)^n chi_j(x) (2-x^2-x^-2)]."""
    p = {0: 1}
    for _ in range(n):
        p = laurent_mul(p, {1: 1, -1: 1})
    chi = {}
    for m2 in range(-twoj, twoj + 1, 2):
        chi[m2] = chi.get(m2, 0) + 1
    p = laurent_mul(p, chi)
    p = laurent_mul(p, {0: 2, 2: -1, -2: -1})
    ct = p.get(0, 0)
    assert ct % 2 == 0
    return ct // 2


def haar_moment(n):
    return haar_mult(n, 0)


# ---------------------------------------------------------------- designs
def sphere_monomial(al):
    if any(a % 2 for a in al):
        return Fr(0)
    num = 1
    for a in al:
        for k in range(a - 1, 0, -2):
            num *= k
    den = 1
    for k in range(sum(al) // 2):
        den *= 4 + 2 * k
    return Fr(num, den)


def group_monomial(elems, d, al):
    tot = QF(0, 0, d)
    for q in elems:
        v = QF(1, 0, d)
        for x, a in zip(q, al):
            for _ in range(a):
                v = v * x
        tot = tot + v
    return QF(tot.a / len(elems), tot.b / len(elems), d)


def dgs_strength(elems, d, maxdeg):
    for deg in range(1, maxdeg + 1):
        for al in itertools.product(range(deg + 1), repeat=4):
            if sum(al) != deg:
                continue
            g = group_monomial(elems, d, al)
            if not (g.b == 0 and g.a == sphere_monomial(al)):
                return deg - 1
    return None


def moments(elems, d, nmax):
    out = []
    for n in range(nmax + 1):
        tot = QF(0, 0, d)
        for q in elems:
            tr = q[0] + q[0]
            v = QF(1, 0, d)
            for _ in range(n):
                v = v * tr
            tot = tot + v
        assert tot.b == 0
        out.append(tot.a / len(elems))
    return out


# ---------------------------------------------------------------- series
def s_mul(a, b, N):
    return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N + 1)]


def s_log(z, N):
    assert z[0] == 1
    l = [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        acc = n * z[n]
        for k in range(1, n):
            acc -= k * l[k] * z[n - k]
        l[n] = acc / n
    return l


def s_div(a, z, N):
    q = [Fr(0)] * (N + 1)
    for n in range(N + 1):
        acc = a[n]
        for k in range(1, n + 1):
            acc -= z[k] * q[n - k]
        q[n] = acc / z[0]
    return q


def s_pow(a, k, N):
    r = [Fr(1)] + [Fr(0)] * N
    for _ in range(k):
        r = s_mul(r, a, N)
    return r


def haar_aj(twoj, N):
    return [Fr(haar_mult(n, twoj), (2 ** n) * factorial(n)) for n in range(N + 1)]


def haar_surface_Z(chiE, F, N):
    Z = [Fr(0)] * (N + 1)
    for twoj in range(0, N + 1):
        dj = twoj + 1
        a = haar_aj(twoj, N)
        r = [x / dj for x in a]
        term = s_pow(r, F, N)
        w = Fr(dj) ** chiE
        Z = [Z[i] + w * term[i] for i in range(N + 1)]
    return Z


# ---------------------------------------------------------------- lattices
def tree_free(nv, links):
    seen = {0}
    tree = set()
    changed = True
    while changed:
        changed = False
        for li, (u, v) in enumerate(links):
            if (u in seen) != (v in seen):
                tree.add(li)
                seen.add(u)
                seen.add(v)
                changed = True
    assert len(seen) == nv
    return tree


def lat_single():
    links = [(0, 1), (1, 2), (3, 2), (0, 3)]
    plaqs = [[(0, 1), (1, 1), (2, -1), (3, -1)]]
    return dict(name="single_plaquette_open", nv=4, links=links, plaqs=plaqs)


def lat_cube():
    verts = list(itertools.product((0, 1), repeat=3))
    vid = {v: i for i, v in enumerate(verts)}
    links = []
    lid = {}
    for v in verts:
        for mu in range(3):
            if v[mu] == 0:
                w = list(v)
                w[mu] = 1
                lid[(v, mu)] = len(links)
                links.append((vid[v], vid[tuple(w)]))
    plaqs = []
    for mu, nu in ((0, 1), (0, 2), (1, 2)):
        for v in verts:
            if v[mu] == 0 and v[nu] == 0:
                vm = list(v); vm[mu] = 1; vm = tuple(vm)
                vn = list(v); vn[nu] = 1; vn = tuple(vn)
                plaqs.append([(lid[(v, mu)], 1), (lid[(vm, nu)], 1),
                              (lid[(vn, mu)], -1), (lid[(v, nu)], -1)])
    return dict(name="cube_surface", nv=8, links=links, plaqs=plaqs, chiE=2, lid=lid)


def lat_torus(L1, L2):
    verts = [(x, y) for x in range(L1) for y in range(L2)]
    vid = {v: i for i, v in enumerate(verts)}
    Ls = (L1, L2)
    links = []
    lid = {}
    for v in verts:
        for mu in range(2):
            w = list(v); w[mu] = (w[mu] + 1) % Ls[mu]
            lid[(v, mu)] = len(links)
            links.append((vid[v], vid[tuple(w)]))
    plaqs = []
    for v in verts:
        vm = ((v[0] + 1) % L1, v[1])
        vn = (v[0], (v[1] + 1) % L2)
        plaqs.append([(lid[(v, 0)], 1), (lid[(vm, 1)], 1), (lid[(vn, 0)], -1), (lid[(v, 1)], -1)])
    return dict(name="torus_%dx%d" % (L1, L2), nv=L1 * L2, links=links, plaqs=plaqs, chiE=0)


def lat_pillow():
    links = [(0, 1), (1, 2), (3, 2), (0, 3)]
    p1 = [(0, 1), (1, 1), (2, -1), (3, -1)]
    p2 = [(l, -e) for (l, e) in reversed(p1)]
    return dict(name="pillow", nv=4, links=links, plaqs=[p1, p2], chiE=2)


def hyp(L):
    ps = [set(l for l, _ in p) for p in L["plaqs"]]
    h1 = all(len(set(l for l, _ in p)) == len(p) == 4 for p in L["plaqs"])
    h2 = all(len(ps[i] & ps[j]) <= 1 for i in range(len(ps)) for j in range(i + 1, len(ps)))
    h2p = True
    for l in range(len(L["links"])):
        star = [i for i in range(len(ps)) if l in ps[i]]
        for i in star:
            others = set().union(*[ps[k] for k in star if k != i]) if len(star) > 1 else set()
            if not any(x != l and x not in others for x in ps[i]):
                h2p = False
    return int(h1), int(h2), int(h2p)


# ---------------------------------------------------------------- group tables
def tables(elems, d):
    idx = {qkey(q): i for i, q in enumerate(elems)}
    n = len(elems)
    mt = [[idx[qkey(qmul(elems[a], elems[b]))] for b in range(n)] for a in range(n)]
    e = idx[qkey((QF(1, 0, d), QF(0, 0, d), QF(0, 0, d), QF(0, 0, d)))]
    inv = [next(b for b in range(n) if mt[a][b] == e) for a in range(n)]
    tr = []
    for q in elems:
        t2 = q[0] + q[0]
        assert t2.b == 0 and t2.a.denominator == 1, "enumerator needs integer traces"
        tr.append(int(t2.a))
    return mt, inv, tr, e


def enumerate_lattice(L, T, track, gauge_fix=True, extra=()):
    """Exact histogram over ALL configurations of the non-tree links (tree links = 1),
    or over all links when gauge_fix=False.  Key: (sum_p tr U_p, tr U_track0, tr U_track1, ...).
    Words in `extra` are tracked observables (indices continue after the plaquettes) and do
    not enter the action.  Innermost free link summed exactly with memoisation of its context."""
    mt, inv, tr, e = T
    n = len(tr)
    tree = tree_free(L["nv"], L["links"]) if gauge_fix else set()
    free = [l for l in range(len(L["links"])) if l not in tree]
    fpos = {l: i for i, l in enumerate(free)}
    words = []
    nplaq = len(L["plaqs"])
    for p in list(L["plaqs"]) + list(extra):
        words.append([(fpos[l], ep) for (l, ep) in p if l not in tree])
    last = len(free) - 1
    inner_p, outer_p = [], []
    for pi, w in enumerate(words):
        cnt = sum(1 for (f, _) in w if f == last)
        if cnt == 0:
            outer_p.append(pi)
        else:
            assert cnt == 1
            k = next(i for i, (f, _) in enumerate(w) if f == last)
            rot = w[k:] + w[:k]
            inner_p.append((pi, rot[0][1], rot[1:]))
    tslot = {p: i for i, p in enumerate(track)}
    memo = {}

    def word_val(wd, cfg):
        g = e
        for (f, ep) in wd:
            x = cfg[f] if ep == 1 else inv[cfg[f]]
            g = mt[g][x]
        return g

    hist = {}
    nt = len(track)
    for outer in itertools.product(range(n), repeat=len(free) - 1):
        cfg = list(outer) + [0]
        base = [0] * (nt + 1)
        for pi in outer_p:
            v = tr[word_val(words[pi], cfg)]
            if pi < nplaq:
                base[0] += v
            if pi in tslot:
                base[1 + tslot[pi]] += v
        ctx = tuple((pi, ep, word_val(rest, cfg)) for (pi, ep, rest) in inner_p)
        ih = memo.get(ctx)
        if ih is None:
            ih = {}
            for x in range(n):
                k = [0] * (nt + 1)
                for (pi, ep, R) in ctx:
                    xx = x if ep == 1 else inv[x]
                    v = tr[mt[xx][R]]
                    if pi < nplaq:
                        k[0] += v
                    if pi in tslot:
                        k[1 + tslot[pi]] += v
                kt = tuple(k)
                ih[kt] = ih.get(kt, 0) + 1
            memo[ctx] = ih
        for kt, c in ih.items():
            key = tuple(b + kk for b, kk in zip(base, kt))
            hist[key] = hist.get(key, 0) + c
    total = n ** len(free)
    assert sum(hist.values()) == total
    return hist, total, len(free)


def series_from_hist(hist, total, N, slot=None):
    """z_n = E[S^n]/n! ; if slot is given, also a_n = E[f_slot S^n]/n!."""
    z = [Fr(0)] * (N + 1)
    a = [Fr(0)] * (N + 1)
    for key, c in hist.items():
        S = Fr(key[0], 2)
        f = Fr(key[1 + slot], 2) if slot is not None else Fr(0)
        pw = Fr(1)
        for k in range(N + 1):
            z[k] += c * pw
            a[k] += c * f * pw
            pw *= S
    z = [z[k] / (total * factorial(k)) for k in range(N + 1)]
    a = [a[k] / (total * factorial(k)) for k in range(N + 1)]
    return z, a


def hist_sha(hist):
    body = "\n".join("%s:%d" % (",".join(map(str, k)), v) for k, v in sorted(hist.items()))
    return hashlib.sha256(body.encode()).hexdigest()


def fstr(x):
    return "%d:%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    N = 8
    G = build_groups()
    expected_order = {"Z6": 6, "Q8": 8, "Dic3": 12, "Dic4": 16, "2T": 24, "2O": 48, "2I": 120}
    haarM = [haar_moment(n) for n in range(14)]
    cat = [comb(n, n // 2) // (n // 2 + 1) if n % 2 == 0 else 0 for n in range(14)]
    row("HAAR_WEYL", M_0_to_13=",".join(map(str, haarM)), equals_catalan=int(haarM == cat))
    info = {}
    for name in ("Z6", "Q8", "Dic3", "Dic4", "2T", "2O", "2I"):
        elems, d = G[name]
        assert len(elems) == expected_order[name], (name, len(elems))
        M = moments(elems, d, 13)
        tm = next(n for n in range(1, 14) if M[n] != haarM[n]) - 1
        td = dgs_strength(elems, d, tm + 1)
        iota = M[tm + 1] - haarM[tm + 1]
        info[name] = dict(t=tm, M=M, iota=iota)
        row("GROUP2", name=name, order=len(elems), field="Q(sqrt%d)" % d, t_moment=tm, t_dgs=td,
            agree=int(tm == td), M_t1=fstr(M[tm + 1]), haar_t1=haarM[tm + 1], iota=fstr(iota),
            M4=fstr(M[4]), M6=fstr(M[6]),
            logZ_per_plaq=fstr(iota / ((2 ** (tm + 1)) * factorial(tm + 1))),
            plaq_coeff=fstr(iota / ((2 ** (tm + 1)) * factorial(tm))))

    # ---- prediction arithmetic: recompute every frozen PAIR row
    predfile = "/mnt/c/tmp/opus55/transfer-20261109/predict/PREDICTION-TRANSFER.hbp"
    alias = {"Dic_3": "Dic3"}
    npair = nbad = 0
    if os.path.exists(predfile):
        for line in open(predfile, encoding="utf-8"):
            if not line.startswith("PAIR|"):
                continue
            kv = dict(x.split("=", 1) for x in line.strip().split("|")[1:] if "=" in x)
            g1, g2 = alias.get(kv["G1"], kv["G1"]), alias.get(kv["G2"], kv["G2"])
            Np = int(kv["N_p"])
            t1 = info[g1]["t"]
            t2 = info[g2]["t"] if g2 != "SU2" else 10 ** 9
            t = min(t1, t2)
            M1 = info[g1]["M"][t + 1]
            M2 = info[g2]["M"][t + 1] if g2 != "SU2" else Fr(haarM[t + 1])
            dl = Np * (M1 - M2) / ((2 ** (t + 1)) * factorial(t + 1))
            dp = (M1 - M2) / ((2 ** (t + 1)) * factorial(t))
            ck = next(k for k in kv if re.match(r"coef_beta\d+_logZ", k))
            pk = next(k for k in kv if re.match(r"plaquette_mean_beta\d+_G1", k))
            fb = int(re.match(r"coef_beta(\d+)_", ck).group(1))
            a, b = kv[ck].split(":")
            pa, pb = kv[pk].split(":")
            ok = (Fr(int(a), int(b)) == dl and Fr(int(pa), int(pb)) == dp and fb == t + 1)
            npair += 1
            nbad += int(not ok)
            if not ok or (kv["lattice"].startswith("cube") and (g1, g2) in
                          (("Dic3", "2T"), ("2O", "SU2"), ("2I", "SU2"))) or \
                    (kv["lattice"].startswith("torus_3x3") and (g1, g2) == ("Q8", "2T")):
                row("PAIR_RECOMPUTE", lattice=kv["lattice"], G1=g1, G2=g2, t=t, beta_power=t + 1,
                    frozen=kv[ck], recomputed=fstr(dl), frozen_pm=kv[pk], recomputed_pm=fstr(dp),
                    match=int(ok))
        row("PAIR_RECOMPUTE_SUMMARY", rows=npair, mismatches=nbad)
    else:
        row("PAIR_RECOMPUTE_SUMMARY", rows=0, missing=1, ok=0)

    # ---- lattice enumerations (Dic3, 2T: integer traces)
    haar_cache = {}
    lats = [lat_single(), lat_pillow(), lat_cube(), lat_torus(2, 2)]
    tabs = {g: tables(*G[g]) for g in ("Dic3", "2T")}
    res = {}
    for L in lats:
        h1, h2, h2p = hyp(L)
        row("LAT2", name=L["name"], links=len(L["links"]), plaqs=len(L["plaqs"]), H1=h1, H2=h2, H2p=h2p)
        # tracked plaquettes: p0 and a neighbour p1 sharing >= 1 link with p0
        p0 = 0
        s0 = set(l for l, _ in L["plaqs"][0])
        p1 = next((i for i in range(1, len(L["plaqs"])) if s0 & set(l for l, _ in L["plaqs"][i])), None)
        track = [p0] + ([p1] if p1 is not None else [])
        shared = len(s0 & set(l for l, _ in L["plaqs"][p1])) if p1 is not None else 0
        for g in ("Dic3", "2T"):
            hist, total, nfree = enumerate_lattice(L, tabs[g], track, gauge_fix=True)
            row("ENUM2", lat=L["name"], group=g, gauge_fixed=1, free_links=nfree, configs=total,
                hist_sha256=hist_sha(hist))
            if L["name"] == "single_plaquette_open":
                hfull, tfull, nf2 = enumerate_lattice(L, tabs[g], track, gauge_fix=False)
                # marginal on the total trace must equal the gauge-fixed one scaled by |G|^3
                m_full, m_fix = {}, {}
                for k, c in hfull.items():
                    m_full[k[0]] = m_full.get(k[0], 0) + c
                for k, c in hist.items():
                    m_fix[k[0]] = m_fix.get(k[0], 0) + c
                scale = tfull // total
                same = all(m_full.get(k, 0) == scale * m_fix.get(k, 0) for k in set(m_full) | set(m_fix))
                row("ENUM2", lat=L["name"], group=g, gauge_fixed=0, free_links=nf2, configs=tfull,
                    hist_sha256=hist_sha(hfull), equals_gauge_fixed_times=scale, equal=int(same))
            z, a = series_from_hist(hist, total, N, slot=0)
            res[(L["name"], g)] = (z, a, hist, total, track, shared)
            if p1 is not None:
                # direct mechanism test: E[f_p0^a f_p1^b] for a+b = 4 (link shared, degree 4 = t+1 for Dic3)
                for (ea, eb) in ((2, 2), (3, 1), (1, 3)):
                    m = Fr(0)
                    for key, c in hist.items():
                        m += c * Fr(key[1], 2) ** ea * Fr(key[2], 2) ** eb
                    m /= total
                    M = info[g]["M"]
                    peel = (M[ea] / 2 ** ea) * (M[eb] / 2 ** eb)
                    row("MIXED", lat=L["name"], group=g, p0=p0, p1=p1, shared_links=shared,
                        moment="E[f_p0^%d f_p1^%d]" % (ea, eb), enum=fstr(m), peeling_value=fstr(peel),
                        equal=int(m == peel))

        # ---- compare with Haar and with each other
        Np = len(L["plaqs"])
        if L["name"] == "single_plaquette_open":
            zH = [Fr(haarM[n], (2 ** n) * factorial(n)) for n in range(N + 1)]
        else:
            zH = haar_surface_Z(L["chiE"], Np, N)
        lH = s_log(zH, N)
        # Haar plaquette mean: single plaquette directly; surfaces via face transitivity (1/F) dlogZ
        if L["name"] == "single_plaquette_open":
            aH = [Fr(haarM[n + 1], (2 ** (n + 1)) * factorial(n)) for n in range(N)] + [Fr(0)]
            fH = s_div(aH, zH, N)
        else:
            fH = [lH[i + 1] * (i + 1) / Np for i in range(N)] + [Fr(0)]
        row("HAAR_SERIES", lat=L["name"], logZ_1_to_8=",".join(fstr(x) for x in lH[1:]))
        h1, h2, h2p = hyp(L)
        for g in ("Dic3", "2T"):
            z, a = res[(L["name"], g)][:2]
            lG = s_log(z, N)
            fG = s_div(a, z, N)
            t = info[g]["t"]
            io = info[g]["iota"]
            dl = [lG[i] - lH[i] for i in range(N + 1)]
            df = [fG[i] - fH[i] for i in range(N + 1)]
            predB = Np * io / ((2 ** (t + 1)) * factorial(t + 1))
            predC = io / ((2 ** (t + 1)) * factorial(t))
            A = all(dl[i] == 0 for i in range(1, t + 1))
            B = dl[t + 1] == predB
            Cp = all(df[i] == 0 for i in range(t)) and df[t] == predC
            row("VS_HAAR", lat=L["name"], group=g, t=t, H1=h1, H2p=h2p,
                dlogZ_1_to_8=",".join(fstr(x) for x in dl[1:]),
                dplaq_0_to_7=",".join(fstr(x) for x in df[:N]),
                predB=fstr(predB), predCprime=fstr(predC), A=int(A), B=int(B), Cprime=int(Cp),
                theorem_applies=int(h1 and h2p))
        z1, a1 = res[(L["name"], "Dic3")][:2]
        z2, a2 = res[(L["name"], "2T")][:2]
        l1, l2 = s_log(z1, N), s_log(z2, N)
        f1, f2 = s_div(a1, z1, N), s_div(a2, z2, N)
        dl = [l1[i] - l2[i] for i in range(N + 1)]
        df = [f1[i] - f2[i] for i in range(N + 1)]
        t = 3
        dM = info["Dic3"]["M"][4] - info["2T"]["M"][4]
        predD = Np * dM / (16 * 24)
        predDp = dM / (16 * 6)
        row("DIC3_VS_2T", lat=L["name"], dlogZ_1_to_8=",".join(fstr(x) for x in dl[1:]),
            dplaq_0_to_7=",".join(fstr(x) for x in df[:N]), predD=fstr(predD), predDp=fstr(predDp),
            D=int(all(dl[i] == 0 for i in range(1, 4)) and dl[4] == predD),
            Dp=int(all(df[i] == 0 for i in range(3)) and df[3] == predDp),
            theorem_applies=int(h1 and h2p))
        # E[S^4] difference, Dic3 vs Haar: theorem says N_p * iota / 2^4
        ES4 = z1[4] * 24 - zH[4] * 24
        row("ES4_DIC3_MINUS_HAAR", lat=L["name"], measured=fstr(ES4),
            constant_sequences_only=fstr(Fr(Np) * info["Dic3"]["iota"] / 16),
            equal=int(ES4 == Fr(Np) * info["Dic3"]["iota"] / 16))

    # ---- claim (C) with d = 2 (O = f_p0^2) and d = 1 (Petrie hexagon Wilson loop), Dic3 vs 2T on the cube
    Lc = lat_cube()
    lid = Lc["lid"]
    hexw = [(lid[((0, 0, 0), 0)], 1), (lid[((1, 0, 0), 1)], 1), (lid[((1, 1, 0), 2)], 1),
            (lid[((0, 1, 1), 0)], -1), (lid[((0, 0, 1), 1)], -1), (lid[((0, 0, 0), 2)], -1)]
    # sanity: the hexagon is a closed walk
    pos = 0
    for (l, ep) in hexw:
        u, v = Lc["links"][l]
        assert (u if ep == 1 else v) == pos
        pos = v if ep == 1 else u
    assert pos == 0
    obs_series = {}
    for g in ("Dic3", "2T"):
        hist, total, _ = enumerate_lattice(Lc, tabs[g], [0, 6], gauge_fix=True, extra=[hexw])
        zc = [Fr(0)] * (N + 1)
        a_sq = [Fr(0)] * (N + 1)
        a_hx = [Fr(0)] * (N + 1)
        for key, c in hist.items():
            S = Fr(key[0], 2)
            fp = Fr(key[1], 2)
            wh = Fr(key[2], 2)
            pw = Fr(1)
            for k in range(N + 1):
                zc[k] += c * pw
                a_sq[k] += c * fp * fp * pw
                a_hx[k] += c * wh * pw
                pw *= S
        zc = [x / (total * factorial(k)) for k, x in enumerate(zc)]
        a_sq = [x / (total * factorial(k)) for k, x in enumerate(a_sq)]
        a_hx = [x / (total * factorial(k)) for k, x in enumerate(a_hx)]
        obs_series[g] = (s_div(a_sq, zc, N), s_div(a_hx, zc, N))
        z_ref = res[("cube_surface", g)][0]
        row("ENUM2_OBS", lat="cube_surface", group=g, configs=total, hist_sha256=hist_sha(hist),
            Z_series_matches_first_run=int(zc == z_ref))
    for oi, (oname, d) in enumerate((("f_p0_squared", 2), ("petrie_hexagon_W_over_2", 1))):
        s1 = obs_series["Dic3"][oi]
        s2 = obs_series["2T"][oi]
        df = [s1[i] - s2[i] for i in range(N + 1)]
        first = next((i for i in range(N + 1) if df[i] != 0), None)
        row("C_DEGREE", lat="cube_surface", g1="Dic3", g2="2T", t=3, observable=oname, d=d,
            predicted_agree_through=3 - d, first_differing_order=first,
            diff_0_to_6=",".join(fstr(x) for x in df[:7]),
            C_holds=int(all(df[i] == 0 for i in range(3 - d + 1))))

    # ---- hypotheses on the H1-breaking tori, under the document's own section 1.2 wording of H2'
    for (a, b) in ((2, 1), (1, 1)):
        Lt = lat_torus(a, b)
        ps = [set(l for l, _ in p) for p in Lt["plaqs"]]
        h1 = int(all(len(ps[i]) == 4 for i in range(len(ps))))
        h2 = int(all(len(ps[i] & ps[j]) <= 1 for i in range(len(ps)) for j in range(i + 1, len(ps))))
        doc = True
        once = True
        for l in range(len(Lt["links"])):
            star = [i for i in range(len(ps)) if l in ps[i]]
            for i in star:
                others = set().union(*[ps[k] for k in star if k != i]) if len(star) > 1 else set()
                if not any(x != l and x not in others for x in ps[i]):
                    doc = False
                cnt = lambda x: sum(1 for (y, _) in Lt["plaqs"][i] if y == x)
                if not any(x != l and x not in others and cnt(x) == 1 for x in ps[i]):
                    once = False
        row("HYP_WORDING", lat=Lt["name"], H1=h1, H2=h2, H2p_as_written_sec1_2=int(doc),
            H2p_with_single_occurrence=int(once))

    # pillow character-formula self-validation: Z = E[e^{beta tr W}] exactly
    zH = haar_surface_Z(2, 2, N)
    direct = [Fr(haarM[n], factorial(n)) for n in range(N + 1)]
    row("CHARFORMULA_VALIDATE", lat="pillow", haar_char_formula_equals_direct=int(zH == direct))

    body = "\n".join(ROWS) + "\n"
    path = os.path.join(outdir, "SECOND-SEAT.hbp")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(body)
    sha = hashlib.sha256(body.encode()).hexdigest()
    with open(path + ".sha256", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("%s  SECOND-SEAT.hbp\n" % sha)
    print(body)
    print("sha256", sha)


if __name__ == "__main__":
    main()
