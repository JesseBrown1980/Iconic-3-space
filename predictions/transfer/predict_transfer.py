# predict_transfer.py -- the PREDICTION for the design-to-lattice transfer theorem, written before any
# measurement is run.
#
# The theorem is stated, corrected and proved in the THEOREM, CORRECTION and PROOF rows this script
# writes. In short: on a finite lattice whose plaquettes each pass through 4 distinct links (H1), with
# link variables uniform on a finite subgroup G of SU(2) or Haar on SU(2) and weight
# exp(beta * sum_p chi(U_p) * (1:2)), a group that is a design of strength t gives log Z_G = log Z_SU2
# through beta^t (A). If also two distinct plaquettes share at most one link (H2), the beta^(t+1)
# coefficients differ by exactly N_p * (M_(t+1)^G - M_(t+1)^Haar) : (2^(t+1) * (t+1)!) (B), and the
# same holds between two finite groups (D).
#
# What it predicts: for every pair among Q8, Dic_3, 2T, 2O, 2I and SU(2), on the 6 faces of one cube
# (free boundary) and on the periodic 3x3 and 4x4 tori, the order through which log Z agrees, the
# first-difference coefficient, and the same for the plaquette mean. On the periodic 2x2 torus, where
# H2 fails, it predicts the agreement order from (A) only, and lists the formula value as a reference.
#
# HARD RULE kept here: no group element, no link variable, no lattice configuration and no sum over
# configurations. Inputs: the closed Molien series of each group and the plaquette count. The only
# lattice object built is the plaquette-link incidence (geometry), to count N_p and to check H1 and H2.
# Sources (NAMED): Molien 1897; Klein 1884; Springer, Invariant Theory, LNM 585, 1977.
#
# Arithmetic: integers only. A rational is a reduced pair (num, den) with den > 0, written num:den,
# reduced by gcd after an exact divisibility assert, with the integer operators // and % only.
# Run: python3 -I -B predict_transfer.py   (writes PREDICTION-TRANSFER.hbp and its .sha256 sidecar
# next to this script; exits 1 if a forced self-check fails)
import hashlib
import os
import sys
from math import gcd

DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
SCHEMA = "TRANSFER_THEOREM_PREDICTION_V1"
OUT_NAME = "PREDICTION-TRANSFER.hbp"
D_MAX = 48  # Molien coefficients a_0..a_48
N_MAX = 24  # moments M_0..M_24
B_LOG = 14  # beta powers of the single-plaquette log route

# ---------------------------------------------------------------- exact rationals


def rat(n, d):
    assert d != 0
    if d < 0:
        n, d = -n, -d
    g = gcd(n, d)
    assert n % g == 0 and d % g == 0
    return (n // g, d // g)


ZERO = (0, 1)
ONE = (1, 1)


def radd(x, y):
    return rat(x[0] * y[1] + y[0] * x[1], x[1] * y[1])


def rsub(x, y):
    return rat(x[0] * y[1] - y[0] * x[1], x[1] * y[1])


def rmul(x, y):
    return rat(x[0] * y[0], x[1] * y[1])


def rint(x, k):
    # x times the integer k
    return rat(x[0] * k, x[1])


def rover(x, k):
    # x times the rational 1:k, k a nonzero integer
    return rat(x[0], x[1] * k)


def rs(x):
    return str(x[0]) + ":" + str(x[1])


def ints(v):
    return ",".join(str(e) for e in v)


def reduced(x):
    return x[1] > 0 and gcd(x[0], x[1]) == 1


# ---------------------------------------------------------------- integer combinatorics

PASCAL = [[1]]
for _n in range(1, 2 * N_MAX + 1):
    _prev = PASCAL[-1]
    PASCAL.append([1] + [_prev[k - 1] + _prev[k] for k in range(1, _n)] + [1])

FACT = [1]
for _n in range(1, N_MAX + 1):
    FACT.append(FACT[-1] * _n)


def binom(n, k):
    if k < 0 or k > n:
        return 0
    return PASCAL[n][k]


def b_closed(n, d):
    # multiplicity of spin j = d:2 in the n-fold tensor power of spin 1:2, B(n,j) = C(n,k) - C(n,k-1), k = (n-d):2
    if d > n or (n - d) % 2 != 0:
        return 0
    k = (n - d) // 2
    return binom(n, k) - binom(n, k - 1)


def b_wrong(n, d):
    # CONTROL: a deliberately wrong closed form, C(n,k) - C(n,k+1)
    if d > n or (n - d) % 2 != 0:
        return 0
    k = (n - d) // 2
    return binom(n, k) - binom(n, k + 1)


def ballot_table(top):
    # second route: tensoring with spin 1:2 moves spin d:2 to (d-1):2 and (d+1):2
    tab = [[1]]
    for n in range(1, top + 1):
        prev = tab[-1]
        new = [0] * (n + 1)
        for d in range(n + 1):
            if d >= 1 and d - 1 <= n - 1:
                new[d] += prev[d - 1]
            if d + 1 <= n - 1:
                new[d] += prev[d + 1]
        tab.append(new)
    return tab


BALLOT = ballot_table(N_MAX)


def b_ballot(n, d):
    if d > n:
        return 0
    return BALLOT[n][d]


# ---------------------------------------------------------------- Molien series (closed formulas only)


def dic(n):
    # binary dihedral Dic_n of order 4n: H(t) = (1 + t^(2n+2)) : ((1 - t^4)(1 - t^(2n)))
    return (2 * n + 2, 4, 2 * n)


def expand(p, q, r, top):
    # (1 + t^p) times the geometric series in t^q times the geometric series in t^r, coefficients of
    # t^0..t^top; multiplying by the geometric series in t^s is the recursion b_d = a_d + b_(d-s)
    a = [0] * (top + 1)
    a[0] += 1
    if p <= top:
        a[p] += 1
    for s in (q, r):
        for d in range(s, top + 1):
            a[d] += a[d - s]
    return a


def enumerate_route(p, q, r, top):
    # second route: count exponent pairs (i, j) with q*i + r*j = d, plus those with q*i + r*j = d - p
    a = [0] * (top + 1)
    for shift in (0, p):
        i = 0
        while shift + q * i <= top:
            j = 0
            while shift + q * i + r * j <= top:
                a[shift + q * i + r * j] += 1
                j += 1
            i += 1
    return a


def first_positive_degree(a):
    for d in range(1, len(a)):
        if a[d] > 0:
            return d
    return -1


# ---------------------------------------------------------------- moments, three routes


def moments_by(a, bfun):
    # M_n = sum_j B(n,j) a_(2j), from chi^n = sum_j B(n,j) chi_j and E[chi_j] = a_(2j)
    return [sum(bfun(n, d) * a[d] for d in range(n + 1)) for n in range(N_MAX + 1)]


def moments_chebyshev(a):
    # third route, free of B: 1 : det(1 - tU) = 1 : (1 - t*chi + t^2) = sum_m (t*chi)^m : (1 + t^2)^(m+1),
    # so H(t) = sum_m M_m t^m (1 + t^2)^-(m+1), a_k = sum_i (-1)^i C(k-i,i) M_(k-2i), solved upward for M_k
    m = []
    for k in range(N_MAX + 1):
        s = a[k]
        i = 1
        while 2 * i <= k:
            sign = -1 if i % 2 == 1 else 1
            s -= sign * binom(k - i, i) * m[k - 2 * i]
            i += 1
        m.append(s)
    return m


def log_route(m):
    # single-plaquette log route: z_n = M_n : (2^n * n!), f = log z by n f_n = n z_n - sum_k k f_k z_(n-k)
    z = [rat(m[n], (1 << n) * FACT[n]) for n in range(B_LOG + 1)]
    assert z[0] == ONE
    f = [ZERO] * (B_LOG + 1)
    for n in range(1, B_LOG + 1):
        s = rint(z[n], n)
        for k in range(1, n):
            s = rsub(s, rmul(rint(f[k], k), z[n - k]))
        f[n] = rover(s, n)
    return f


# ---------------------------------------------------------------- lattice geometry (incidence only)


def cube_lattice():
    # vertices are 3-bit integers; a link is (v, axis) with bit axis of v clear; 6 faces
    plaqs = []
    for fixed in range(3):
        a1, a2 = [ax for ax in range(3) if ax != fixed]
        for side in range(2):
            v = side << fixed
            plaqs.append([(v, a1), (v | (1 << a1), a2), (v | (1 << a2), a1), (v, a2)])
    return plaqs


def cube_ends(link):
    v, ax = link
    return (v, v | (1 << ax))


def torus_lattice(n):
    plaqs = []
    for x in range(n):
        for y in range(n):
            plaqs.append([(x, y, 0), ((x + 1) % n, y, 1), (x, (y + 1) % n, 0), (x, y, 1)])
    return plaqs


def torus_ends(link, n):
    x, y, d = link
    if d == 0:
        return ((x, y), ((x + 1) % n, y))
    return ((x, y), (x, (y + 1) % n))


def geometry(plaqs, ends):
    links = set()
    verts = set()
    distinct4 = True
    cycle4 = True
    for p in plaqs:
        if len(set(p)) != 4:
            distinct4 = False
        deg = {}
        for l in p:
            links.add(l)
            for v in ends(l):
                verts.add(v)
                deg[v] = deg.get(v, 0) + 1
        if len(deg) != 4 or any(c != 2 for c in deg.values()):
            cycle4 = False
    hist = {}
    max_share = 0
    for i in range(len(plaqs)):
        for j in range(i + 1, len(plaqs)):
            s = len(set(plaqs[i]) & set(plaqs[j]))
            hist[s] = hist.get(s, 0) + 1
            if s > max_share:
                max_share = s
    per_link = sorted(set(sum(1 for p in plaqs if l in p) for l in links))
    distinct_plaqs = len(set(frozenset(p) for p in plaqs)) == len(plaqs)
    return {
        "n_p": len(plaqs), "n_links": len(links), "n_verts": len(verts),
        "euler": len(verts) - len(links) + len(plaqs), "distinct4": distinct4, "cycle4": cycle4,
        "distinct_plaqs": distinct_plaqs, "max_share": max_share,
        "hist": ",".join(str(s) + "x" + str(hist[s]) for s in sorted(hist)), "per_link": per_link,
    }


# ---------------------------------------------------------------- the theorem, as frozen text

THEOREM = [
    ("setting", "none",
     "Lambda a finite lattice with a finite list of N_p plaquettes, each a closed path of links; link variables U_l independent and "
     "uniform on K, where K is SU(2) with Haar measure or a finite subgroup G of SU(2) with the uniform measure; U_p the ordered product "
     "of the U_l or their inverses around p, its trace independent of start and orientation since chi(X) = chi(X^-1) is real on SU(2); "
     "S = sum_p chi(U_p)*(1:2); Z_K(beta) = E_K[exp(beta*S)] = sum_n beta^n*E_K[S^n]*(1:n!), a power series with Z_K(0) = 1, so log Z_K "
     "is a power series; its beta^n coefficient is written c_n and the beta^n coefficient of Z_K is written z_n"),
    ("design", "none",
     "a_d = dimension of the G-invariants in spin d:2 = coefficient of t^d in the Molien series of G; M_n^K = E_K[chi^n]; t_G = the largest "
     "t with E_G[D^j] = 0 for every spin j with 0 < 2j <= t, equivalently a_d = 0 for 0 < d <= t, equivalently M_n^G = M_n^Haar for n <= t "
     "(both directions proved, PROOF step 2); t_G + 1 = the smallest positive invariant degree; t_SU2 = inf"),
    ("hypotheses", "none",
     "H1: each plaquette passes through 4 distinct links, each once. H2: two distinct plaquettes share at most one link"),
    ("A", "H1",
     "log Z_G and log Z_SU2 agree as power series through beta^t, t = t_G"),
    ("B", "H1,H2",
     "the beta^(t+1) coefficients of log Z_G and log Z_SU2 differ by exactly N_p*(M_(t+1)^G - M_(t+1)^Haar):(2^(t+1)*(t+1)!) = "
     "N_p*a_(t+1):(2^(t+1)*(t+1)!), which is positive; in z_(t+1) only the N_p single-plaquette terms chi(U_p)^(t+1) differ"),
    ("C", "H1",
     "O a polynomial of degree at most d <= t in the entries of each U_l (entries of U_l^-1 are entries of U_l up to sign on SU(2)): "
     "<O>_G and <O>_SU2 agree through beta^(t-d); under H1 and H2 the plaquette mean S:N_p (d = 1) first differs at beta^t, by "
     "a_(t+1):(2^(t+1)*t!), so the order t-d is attained for d = 1"),
    ("D", "H1,H2",
     "finite G1, G2 with t = min(t1,t2): log Z_G1 and log Z_G2 agree through beta^t and their beta^(t+1) coefficients differ by exactly "
     "N_p*(M_(t+1)^G1 - M_(t+1)^G2):(2^(t+1)*(t+1)!); for t1 != t2 the moment difference is a_(t+1) of the group with the smaller t, at "
     "least 1; for t1 = t2 the same formula holds and is the first difference exactly when the two moments differ"),
]

CORRECTIONS = [
    ("1", "A,C",
     "the task stated A and C under the hypothesis that two distinct plaquettes share at most one link",
     "A and C need H1 only; H2 enters B and the first-difference part of D only (PROOF step 4 uses no H2)",
     "A predicts agreement through beta^t on the periodic 2x2 torus too, where H2 fails; CONTROL_PAIR rows carry that prediction"),
    ("2", "D",
     "the task gave the first difference of D for t1 != t2",
     "the beta^(t+1) formula holds for t1 = t2 as well; it is the first difference exactly when M_(t+1)^G1 != M_(t+1)^G2, which always "
     "holds for t1 != t2 because that difference is a_(t1+1) >= 1",
     "Q8 vs Dic_3 (t = 3 both, a_4 = 2 and 1) differ at beta^4 by the formula; Dic_3 vs Dic_4 (t = 3 both, a_4 = 1 and 1) have beta^4 "
     "difference 0, and the theorem then gives the agreement through beta^4 and no first difference (CORRECTION_EXAMPLE row)"),
    ("3", "setting",
     "the task allowed any boundary condition",
     "any identification of links (free, periodic, any gluing) with H1 and H2 checked on the result; a twisted boundary condition, a "
     "center element z_p in (1,-1) inserted in plaquette p, replaces N_p by sum_p z_p^(t+1), by the same proof with chi(z_p U_p) = z_p chi(U_p)",
     "for the six groups here every group contains -1, t is odd and t+1 even, so twisted and untwisted values coincide"),
    ("4", "C",
     "the task gave C for every d",
     "C holds for d <= t; for d > t it gives nothing",
     "the plaquette mean has d = 1, so it is covered for every group here"),
    ("5", "design equivalence",
     "the task left the design-moment equivalence NAMED",
     "both directions are proved in PROOF step 2: E_G[D^j] = 0 for 0 < 2j <= t iff a_d = 0 for 0 < d <= t iff M_n^G = M_n^Haar for "
     "n <= t; the proof of A to D uses design implies equal polynomial averages (PROOF step 1)",
     "the design strengths here are read from the Molien series directly, t = (smallest positive invariant degree) - 1"),
]

PROOF = [
    ("1", "polynomial averages",
     "On SU(2), U^-1 = U^dagger has the entries of U up to sign, so a polynomial of degree at most k in the entries of U and U^-1 is a "
     "linear combination of entries of the m-fold tensor power of U, m <= k. That tensor power decomposes into spins j with 2j <= m "
     "(Clebsch-Gordan, multiplicity B(m,j)), so the polynomial is a combination of matrix coefficients D^j_ab(U) with 2j <= k. "
     "E_K[D^0] = 1 for every K. If E_G[D^j] = 0 = E_Haar[D^j] for 0 < 2j <= t, then E_G[F] = E_Haar[F] for every polynomial F of degree "
     "at most t. The 4 entries of U and the 4 real coordinates of the point U of S^3 in R^4 are linear combinations of each other, so "
     "this is the spherical t-design property of G on S^3"),
    ("2", "design equals moments",
     "P_j = E_G[D^j(g)] satisfies P_j P_j = P_j and D^j(h) P_j = P_j for h in G (G is a group), so P_j is the projection onto the "
     "G-invariants of spin j and E_G[chi_j] = trace P_j = a_(2j) (Molien); P_j = 0 exactly when a_(2j) = 0. chi^n = sum_j B(n,j) chi_j "
     "gives M_n^K = sum_j B(n,j) a_(2j)^K with a^Haar = (1,0,0,...), so M_n^Haar = B(n,0) = Catalan(n:2) for even n and 0 for odd n. If "
     "a_d = 0 for 0 < d <= t then M_n^G = M_n^Haar for n <= t. Conversely if M_n^G = M_n^Haar for n <= t, then sum over j > 0 of "
     "B(n,j) a_(2j) = 0 with every term >= 0 and B(n,n:2) = 1, so a_n = 0 for 0 < n <= t. With the lower a zero, "
     "M_(t+1)^G - M_(t+1)^Haar = B(t+1,(t+1):2) a_(t+1) = a_(t+1)"),
    ("3", "expansion",
     "z_n = (1:(n!*2^n)) times the sum over ordered sequences (p_1,...,p_n) of plaquettes of E_K[chi(U_(p_1))...chi(U_(p_n))]. Under H1 "
     "chi(U_p) has degree 1 in each of the 4 links of p and degree 0 in the others, so the integrand of a sequence has degree "
     "d_l = number of i with l in p_i, d_l <= n, in U_l"),
    ("4", "hybrid, A and C",
     "Change the measure from G to Haar one link at a time (Fubini on a finite product of probability spaces). With the other links "
     "fixed, the integrand is a polynomial in U_l of degree d_l whose coefficients depend on the other links, and averaging over U_l "
     "leaves a polynomial in the remaining links with no larger degrees. If every d_l <= t, each step keeps the expectation by step 1. "
     "For n <= t every d_l <= n <= t, so z_n^G = z_n^Haar for n <= t. c_n = z_n plus a polynomial in z_1..z_(n-1) (log(1 + (Z - 1))), so "
     "c_n^G = c_n^Haar for n <= t: A. For O of degree at most d per link, E_K[O*S^m] has degree at most d + m per link and agrees for "
     "m <= t - d; <O>_K = E_K[O*exp(beta*S)]:Z_K is a ratio of series agreeing through beta^(t-d) with Z_K(0) = 1: C. H2 is unused"),
    ("5", "order t+1",
     "At n = t+1 a sequence can differ only if some link l has d_l = t+1, that is, l lies in every p_i. If all p_i are one plaquette p "
     "the term is E_K[chi(U_p)^(t+1)]. Otherwise let q_1..q_r, r >= 2, be its distinct plaquettes with multiplicities a_k summing to t+1, "
     "so every a_k <= t; each q_k contains l, and by H2 two of them share only l, so the 3 other links of q_k lie in no other q_m"),
    ("6", "mixed terms agree",
     "Rotate, and invert if needed, each U_(q_k) to the form U_l V_k with V_k the ordered product of the 3 other links of q_k; chi is "
     "unchanged. For fixed U_l in K the V_k depend on disjoint sets of links, so they are independent, and each is uniform on K "
     "(multiplying one of its links by a fixed element of K is a measure-preserving bijection). So the U_l V_k are independent and uniform "
     "on K and the term equals the product over k of M_(a_k)^K whatever U_l is. Since a_k <= t, M_(a_k)^G = M_(a_k)^Haar"),
    ("7", "single plaquette, B",
     "The sequence (p,...,p) gives E_K[chi(U_p)^(t+1)] = M_(t+1)^K because U_p is uniform on K; there are N_p such sequences. So "
     "z_(t+1)^G - z_(t+1)^Haar = N_p*(M_(t+1)^G - M_(t+1)^Haar):((t+1)!*2^(t+1)), and since z_1..z_t agree (step 4) c_(t+1) differs by the "
     "same amount: B. By step 2 the moment difference is a_(t+1) >= 1"),
    ("8", "two groups, D",
     "For t1 < t2, G2 agrees with Haar at beta^(t1+1) (A for G2) and M_(t1+1)^G2 = M_(t1+1)^Haar, so the G1 minus G2 difference at "
     "beta^(t1+1) is the B difference of G1 and equals N_p*(M^G1 - M^G2)_(t+1):(2^(t+1)*(t+1)!). For t1 = t2 = t subtract the two B "
     "statements; the result is the same formula, and it can be 0 (Dic_3 and Dic_4): D"),
    ("9", "plaquette mean",
     "d(log Z):d(beta) = E[S*exp(beta*S)]:Z = N_p times the plaquette mean <S:N_p>, so the beta^t coefficient of the plaquette mean is "
     "(t+1)*c_(t+1):N_p, and under H1 and H2 its G vs Haar difference is a_(t+1):(2^(t+1)*t!) > 0; below beta^t it agrees by A"),
]


# ---------------------------------------------------------------- run


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    checks = []

    def check(name, ok, detail, forced=True):
        checks.append((name, bool(ok), detail, forced))

    groups = [
        {"name": "Q8", "alias": "Dic_2", "family": "binary_dihedral_n=2", "order": 8, "pqr": dic(2)},
        {"name": "Dic_3", "alias": "", "family": "binary_dihedral_n=3", "order": 12, "pqr": dic(3)},
        {"name": "2T", "alias": "", "family": "binary_tetrahedral", "order": 24, "pqr": (12, 6, 8)},
        {"name": "2O", "alias": "", "family": "binary_octahedral", "order": 48, "pqr": (18, 8, 12)},
        {"name": "2I", "alias": "", "family": "binary_icosahedral", "order": 120, "pqr": (30, 12, 20)},
    ]

    a_haar = [1] + [0] * D_MAX
    m_haar = moments_by(a_haar, b_closed)
    T = {"SU2": None}
    MOM = {"SU2": m_haar}
    AV = {"SU2": a_haar}

    two_route_bad, a0_bad, order_bad, feq_bad, parity_bad, tmin_bad = [], [], [], [], [], []
    mroute_bad, mcheb_bad, tmom_bad, jump_bad, wrong_detected = [], [], [], [], []
    rows_g = []

    for g in groups:
        name = g["name"]
        p, q, r = g["pqr"]
        a1 = expand(p, q, r, D_MAX)
        a2 = enumerate_route(p, q, r, D_MAX)
        if a1 != a2:
            two_route_bad.append(name)
        if a1[0] != 1:
            a0_bad.append(name)
        if q * r != 2 * g["order"]:
            order_bad.append(name)
        if p != q + r - 2:
            feq_bad.append(name)
        if any(a1[d] != 0 for d in range(1, D_MAX + 1, 2)):
            parity_bad.append(name)
        d0 = first_positive_degree(a1)
        t = d0 - 1
        if t != min(p, q, r) - 1:
            tmin_bad.append(name)
        m1 = moments_by(a1, b_closed)
        m2 = moments_by(a2, b_ballot)
        m3 = moments_chebyshev(a1)
        if m1 != m2:
            mroute_bad.append(name)
        if m1 != m3:
            mcheb_bad.append(name)
        n_first = next((n for n in range(N_MAX + 1) if m1[n] != m_haar[n]), -1)
        if n_first - 1 != t:
            tmom_bad.append(name)
        if n_first < 0 or m1[n_first] - m_haar[n_first] != a1[d0]:
            jump_bad.append(name)
        if moments_by(a1, b_wrong) != m3:
            wrong_detected.append(name)
        T[name] = t
        MOM[name] = m1
        AV[name] = a1

        label = ["group=" + name] + (["alias=" + g["alias"]] if g["alias"] else [])
        base = label + ["family=" + g["family"], "order=" + str(g["order"]) + "(NAMED)"]
        rows_g.append("|".join(["MOLIEN"] + base + [
            "H=(1+t^" + str(p) + "):((1-t^" + str(q) + ")(1-t^" + str(r) + "))", "d_range=0.." + str(D_MAX),
            "a=" + ints(a1), "json=0"]))
        rows_g.append("|".join(["DESIGN"] + base + [
            "t=" + str(t), "smallest_positive_invariant_degree=" + str(d0), "a_(t+1)=" + str(a1[d0]),
            "rule=t=min(" + str(p) + "," + str(q) + "," + str(r) + ")-1",
            "meaning=E_G[D^j]=0_for_0<2j<=t", "tag=PREDICTION_INPUT", "json=0"]))
        rows_g.append("|".join(["MOMENT"] + base + [
            "quantity=E_G[chi^n]", "n_range=0.." + str(N_MAX), "M=" + ints(m1),
            "first_n_with_M_ne_Haar=" + str(n_first), "M_minus_Haar_there=" + str(m1[n_first] - m_haar[n_first]),
            "routes=binomial_B,ballot_B,chebyshev_inversion", "json=0"]))

    # Haar rows and checks
    mh_ballot = moments_by(a_haar, b_ballot)
    mh_cheb = moments_chebyshev(a_haar)
    if mh_ballot != m_haar:
        mroute_bad.append("SU2")
    if mh_cheb != m_haar:
        mcheb_bad.append("SU2")
    haar_wrong_same = moments_by(a_haar, b_wrong) == m_haar

    cat_closed_ok = True
    for n in range(N_MAX + 1):
        if n % 2 == 1:
            want = 0
        else:
            k = n // 2
            c = binom(2 * k, k)
            assert c % (k + 1) == 0
            want = c // (k + 1)
        if m_haar[n] != want:
            cat_closed_ok = False
    cat = [1]
    for k in range(N_MAX // 2):
        cat.append(sum(cat[i] * cat[k - i] for i in range(k + 1)))
    cat_rec_ok = all(m_haar[n] == (cat[n // 2] if n % 2 == 0 else 0) for n in range(N_MAX + 1))
    rows_h = ["|".join(["MOMENT_HAAR", "group=SU2_Haar", "rule=M_n=B(n,0)=Catalan(n:2)_even_n_0_odd_n",
                        "n_range=0.." + str(N_MAX), "M=" + ints(m_haar), "t=inf", "json=0"])]

    # Dic_4 for correction 2 (closed formula only)
    p4, q4, r4 = dic(4)
    a_d4 = expand(p4, q4, r4, D_MAX)
    m_d4 = moments_by(a_d4, b_closed)
    t_d4 = first_positive_degree(a_d4) - 1
    d4_ok = (a_d4 == enumerate_route(p4, q4, r4, D_MAX) and m_d4 == moments_chebyshev(a_d4) and q4 * r4 == 2 * 16
             and p4 == q4 + r4 - 2)
    row_d4 = "|".join([
        "CORRECTION_EXAMPLE", "for=CORRECTION_2", "groups=Dic_3,Dic_4", "Dic_4_H=(1+t^" + str(p4) + "):((1-t^" + str(q4) + ")(1-t^"
        + str(r4) + "))", "Dic_4_order=16(NAMED)", "t_Dic_3=" + str(T["Dic_3"]), "t_Dic_4=" + str(t_d4),
        "a4_Dic_3=" + str(AV["Dic_3"][4]), "a4_Dic_4=" + str(a_d4[4]), "M4_Dic_3=" + str(MOM["Dic_3"][4]), "M4_Dic_4=" + str(m_d4[4]),
        "beta^4_logZ_difference_any_lattice_with_H1_H2=N_p*(" + str(MOM["Dic_3"][4]) + "-" + str(m_d4[4]) + "):(2^4*4!)="
        + rs(rat(MOM["Dic_3"][4] - m_d4[4], (1 << 4) * FACT[4])) + "_times_N_p",
        "first_difference=NOT_GIVEN_by_the_theorem", "agree_through_beta_power=4", "tag=PREDICTION", "json=0"])

    # lattices
    lats = [
        ("cube_6_faces_free", "free", cube_lattice(), cube_ends, 2, 1, 6, 12, 8),
        ("torus_3x3_periodic", "periodic", torus_lattice(3), lambda l: torus_ends(l, 3), 0, 1, 9, 18, 9),
        ("torus_4x4_periodic", "periodic", torus_lattice(4), lambda l: torus_ends(l, 4), 0, 1, 16, 32, 16),
        ("torus_2x2_periodic", "periodic", torus_lattice(2), lambda l: torus_ends(l, 2), 0, 2, 4, 8, 4),
    ]
    geo_bad = []
    rows_l = []
    GEO = {}
    for lname, bc, plaqs, ends, e_euler, e_share, e_np, e_links, e_verts in lats:
        geo = geometry(plaqs, ends)
        GEO[lname] = geo
        ok = (geo["n_p"] == e_np and geo["n_links"] == e_links and geo["n_verts"] == e_verts and geo["euler"] == e_euler
              and geo["distinct4"] and geo["cycle4"] and geo["distinct_plaqs"] and geo["max_share"] == e_share
              and geo["per_link"] == [2])
        if not ok:
            geo_bad.append(lname)
        h1 = geo["distinct4"] and geo["cycle4"]
        h2 = geo["max_share"] <= 1
        rows_l.append("|".join([
            "LATTICE", "name=" + lname, "boundary=" + bc, "vertices=" + str(geo["n_verts"]), "links=" + str(geo["n_links"]),
            "N_p=" + str(geo["n_p"]), "euler_V-E+F=" + str(geo["euler"]), "plaquettes_per_link=" + ints(geo["per_link"]),
            "plaquette_pairs_by_shared_links=" + geo["hist"], "max_links_shared_by_two_distinct_plaquettes=" + str(geo["max_share"]),
            "H1=" + ("HOLDS" if h1 else "FAILS"), "H2=" + ("HOLDS" if h2 else "FAILS"),
            "built=plaquette_link_incidence_only", "link_variables_assigned=0", "json=0"]))

    # single-plaquette log route per group
    F = {nm: log_route(MOM[nm]) for nm in MOM}

    names = ["Q8", "Dic_3", "2T", "2O", "2I", "SU2"]
    pairs = [(names[i], names[j]) for i in range(len(names)) for j in range(i + 1, len(names))]

    def tstr(x):
        return "inf" if x is None else str(x)

    def pair_core(n1, n2):
        t1, t2 = T[n1], T[n2]
        t = min(x for x in (t1, t2) if x is not None)
        k = t + 1
        mg1 = MOM[n1][k]
        mg2 = MOM[n2][k]
        if n2 == "SU2":
            case = "G_vs_SU2_part_B"
        elif t1 == t2:
            case = "t1_eq_t2_part_D_as_corrected"
        else:
            case = "t1_ne_t2_part_D"
        return t, k, mg1, mg2, mg1 - mg2, case

    rows_p = []
    rows_c = []
    log_bad, pm_bad, red_bad = [], [], []
    PRED = {}
    for lname, bc, plaqs, ends, e_euler, e_share, e_np, e_links, e_verts in lats:
        geo = GEO[lname]
        npl = geo["n_p"]
        h2 = geo["max_share"] <= 1
        for n1, n2 in pairs:
            t, k, mg1, mg2, dm, case = pair_core(n1, n2)
            den = (1 << k) * FACT[k]
            coef = rat(npl * dm, den)
            pm = rat(dm, (1 << k) * FACT[t])
            for x in (coef, pm):
                if not reduced(x):
                    red_bad.append(lname + ":" + n1 + "-" + n2)
            if rat(coef[0] * k, coef[1] * npl) != pm:
                pm_bad.append(lname + ":" + n1 + "-" + n2)
            diff = [rsub(F[n1][i], F[n2][i]) for i in range(B_LOG + 1)]
            zero_through = all(diff[i] == ZERO for i in range(k))
            if not (zero_through and rint(diff[k], npl) == coef):
                log_bad.append(lname + ":" + n1 + "-" + n2)
            PRED[(lname, n1, n2)] = (t, k, coef, pm)
            formula = (str(npl) + "*(" + str(mg1) + "-" + str(mg2) + "):(2^" + str(k) + "*" + str(k) + "!)=" + str(npl * dm) + ":"
                       + str(den))
            head = ["lattice=" + lname, "N_p=" + str(npl), "G1=" + n1, "G2=" + n2, "t1=" + tstr(T[n1]), "t2=" + tstr(T[n2]),
                    "t=min(t1,t2)=" + str(t)]
            if h2:
                first = str(k) if dm != 0 else "NOT_GIVEN_agree_through_" + str(k)
                rows_p.append("|".join(["PAIR"] + head + [
                    "H1=HOLDS", "H2=HOLDS", "case=" + case,
                    "logZ_agree_through_beta_power=" + str(t), "first_diff_beta_power=" + first,
                    "M" + str(k) + "_G1=" + str(mg1), "M" + str(k) + "_G2=" + str(mg2),
                    "coef_beta" + str(k) + "_logZ_G1_minus_G2=" + rs(coef), "formula=" + formula,
                    "plaquette_mean_agree_through_beta_power=" + str(t - 1),
                    "plaquette_mean_beta" + str(t) + "_G1_minus_G2=" + rs(pm),
                    "tag=PREDICTION", "json=0"]))
            else:
                rows_c.append("|".join(["CONTROL_PAIR"] + head + [
                    "H1=HOLDS", "H2=FAILS_max_shared_" + str(geo["max_share"]),
                    "logZ_agree_through_beta_power=" + str(t), "agree_tag=PREDICTION_from_part_A_which_uses_H1_only(CORRECTION_1)",
                    "plaquette_mean_agree_through_beta_power=" + str(t - 1),
                    "coef_beta" + str(k) + "_difference=NOT_PREDICTED_formula_hypothesis_fails",
                    "reference_formula_value=" + rs(coef), "reference_formula=" + formula,
                    "reference_tag=REFERENCE_ONLY_NOT_A_PREDICTION_measurement_decides_hold_or_break", "json=0"]))

    g22 = GEO["torus_2x2_periodic"]
    row_ctl = "|".join([
        "CONTROL", "lattice=torus_2x2_periodic", "N_p=" + str(g22["n_p"]), "links=" + str(g22["n_links"]),
        "plaquette_pairs_by_shared_links=" + g22["hist"], "max_links_shared=" + str(g22["max_share"]),
        "formula_hypothesis=FAILS", "H2=FAILS", "H1=" + ("HOLDS" if g22["distinct4"] and g22["cycle4"] else "FAILS"),
        "coefficient_prediction=NONE", "predicted_only=logZ_agreement_through_beta^t_from_part_A",
        "why_B_proof_breaks=two_adjacent_plaquettes_share_2_links_so_step_5_mixed_terms_carry_two_links_of_degree_t+1",
        "tag=CONTROL", "json=0"])

    # wrong controls
    wd = [nm for nm in ("Q8", "Dic_3", "2T", "2O", "2I") if nm in wrong_detected]
    dim_wrong_fail = [n for n in range(1, N_MAX + 1) if sum(b_wrong(n, d) * (d + 1) for d in range(n + 1)) != (1 << n)]
    t_c, k_c, coef_c, pm_c = PRED[("cube_6_faces_free", "Dic_3", "2T")]
    dm_c = MOM["Dic_3"][k_c] - MOM["2T"][k_c]
    wrong_norm = rat(6 * dm_c, (1 << t_c) * FACT[k_c])
    lr_c = rint(rsub(F["Dic_3"][k_c], F["2T"][k_c]), 6)
    norm_detected = wrong_norm != lr_c
    pw, qw, rw = 12, 4, 8
    aw = expand(pw, qw, rw, D_MAX)
    molien_order_flag = qw * rw != 2 * 24
    molien_feq_flag = pw != qw + rw - 2
    rows_w = [
        "|".join(["WRONG_CONTROL", "name=ballot_closed_form_shifted", "B_wrong(n,j)=C(n,k)-C(n,k+1)",
                  "detected_against_chebyshev_route_in=" + (",".join(wd) or "none"),
                  "haar_row_unchanged=" + ("1" if haar_wrong_same else "0"),
                  "why_haar_unchanged=C(n,n:2-1)=C(n,n:2+1)_so_the_j=0_entry_agrees",
                  "dimension_sum_fails_first_at_n=" + (str(dim_wrong_fail[0]) if dim_wrong_fail else "none"), "json=0"]),
        "|".join(["WRONG_CONTROL", "name=normalization_2^t_for_2^(t+1)", "pair=Dic_3_vs_2T", "lattice=cube_6_faces_free",
                  "wrong_value=" + rs(wrong_norm), "single_plaquette_log_route=" + rs(lr_c),
                  "detected=" + ("1" if norm_detected else "0"), "json=0"]),
        "|".join(["WRONG_CONTROL", "name=2T_molien_with_q=4", "H_wrong=(1+t^12):((1-t^4)(1-t^8))",
                  "q*r=" + str(qw * rw) + "_vs_2*order=48", "p=" + str(pw) + "_vs_q+r-2=" + str(qw + rw - 2),
                  "t_it_would_give=" + str(first_positive_degree(aw) - 1),
                  "order_check_flags=" + ("1" if molien_order_flag else "0"),
                  "functional_equation_flags=" + ("1" if molien_feq_flag else "0"), "json=0"]),
    ]

    # self-checks
    gl = "Q8,Dic_3,2T,2O,2I"
    check("molien_two_routes_agree", not two_route_bad, "geometric_recursion_vs_exponent_pairs_d0.." + str(D_MAX) + "_groups=" + gl
          + "_bad=" + (",".join(two_route_bad) or "none"))
    check("molien_a0_eq_1", not a0_bad, "bad=" + (",".join(a0_bad) or "none"))
    check("order_q*r_eq_2*order", not order_bad, "orders_NAMED=8,12,24,48,120_bad=" + (",".join(order_bad) or "none"))
    check("functional_equation_p_eq_q+r-2", not feq_bad, "H(1:t)=t^2*H(t)_for_G_in_SU2_bad=" + (",".join(feq_bad) or "none"))
    check("contains_minus_one_odd_a_zero", not parity_bad, "bad=" + (",".join(parity_bad) or "none"))
    check("t_eq_min_degree_minus_1", not tmin_bad, "bad=" + (",".join(tmin_bad) or "none"))
    cg_ok = all(BALLOT[n][d] == b_closed(n, d) for n in range(N_MAX + 1) for d in range(n + 1))
    dim_ok = all(sum(b_closed(n, d) * (d + 1) for d in range(n + 1)) == (1 << n) for n in range(N_MAX + 1))
    check("B_binomial_eq_B_ballot_recursion", cg_ok, "n=0.." + str(N_MAX))
    check("B_dimension_sum_eq_2^n", dim_ok, "n=0.." + str(N_MAX))
    check("moments_binomial_route_eq_ballot_route", not mroute_bad, "a_from_recursion_with_binomial_B_vs_a_from_exponent_pairs_with_"
          "ballot_B_groups=" + gl + ",SU2_bad=" + (",".join(mroute_bad) or "none"))
    check("moments_eq_chebyshev_inversion_route", not mcheb_bad, "route_free_of_B_groups=" + gl + ",SU2_bad="
          + (",".join(mcheb_bad) or "none"))
    check("haar_eq_catalan_closed_form", cat_closed_ok, "C(2k,k):(k+1)_exact_n=0.." + str(N_MAX))
    check("haar_eq_catalan_recurrence", cat_rec_ok, "n=0.." + str(N_MAX) + "_M24=" + str(m_haar[N_MAX]))
    check("t_from_molien_eq_t_from_moments", not tmom_bad, "bad=" + (",".join(tmom_bad) or "none"))
    check("moment_jump_at_t+1_eq_a_(t+1)", not jump_bad, "bad=" + (",".join(jump_bad) or "none"))
    hand_t = {"Q8": 3, "Dic_3": 3, "2T": 5, "2O": 7, "2I": 11}
    hand_a = {"Q8": 2, "Dic_3": 1, "2T": 1, "2O": 1, "2I": 1}
    ht_ok = all(T[nm] == hand_t[nm] and AV[nm][hand_t[nm] + 1] == hand_a[nm] for nm in hand_t)
    check("hand_t_and_a_(t+1)", ht_ok, "t=" + ",".join(nm + ":" + str(T[nm]) for nm in hand_t) + "_a_(t+1)="
          + ",".join(nm + ":" + str(AV[nm][T[nm] + 1]) for nm in hand_t))
    hand_m = [("Q8", 4, 4), ("Dic_3", 4, 3), ("SU2", 4, 2), ("2T", 6, 6), ("SU2", 6, 5), ("2O", 8, 15), ("SU2", 8, 14),
              ("2I", 12, 133), ("SU2", 12, 132)]
    check("hand_moments", all(MOM[nm][n] == v for nm, n, v in hand_m),
          "computed=" + ",".join(nm + "_M" + str(n) + ":" + str(MOM[nm][n]) for nm, n, v in hand_m))
    check("lattice_geometry_as_intended", not geo_bad, "cube_V8_E12_F6_euler2,tori_euler0,each_plaquette_a_4-cycle_of_distinct_links,"
          "every_link_in_2_plaquettes_bad=" + (",".join(geo_bad) or "none"))
    h2_ok = all(GEO[nm]["max_share"] == 1 for nm in ("cube_6_faces_free", "torus_3x3_periodic", "torus_4x4_periodic"))
    check("H2_holds_cube_t3_t4", h2_ok, "max_shared=" + ",".join(str(GEO[nm]["max_share"]) for nm in
                                                               ("cube_6_faces_free", "torus_3x3_periodic", "torus_4x4_periodic")))
    check("control_2x2_flagged_H2_fails", GEO["torus_2x2_periodic"]["max_share"] == 2,
          "max_shared=" + str(GEO["torus_2x2_periodic"]["max_share"]) + "_hist=" + GEO["torus_2x2_periodic"]["hist"])
    check("single_plaquette_log_route_eq_formula", not log_bad, "N_p*(f_G1-f_G2)_zero_through_t_and_eq_coef_at_t+1_rows="
          + str(len(lats) * len(pairs)) + "_bad=" + (",".join(log_bad) or "none"))
    check("plaquette_mean_eq_(t+1)*coef:N_p", not pm_bad, "bad=" + (",".join(pm_bad) or "none"))
    check("rationals_reduced_den_positive", not red_bad, "bad=" + (",".join(red_bad) or "none"))
    check("correction_2_Dic_4_routes", d4_ok and t_d4 == 3 and a_d4[4] == 1 and m_d4[4] == MOM["Dic_3"][4],
          "t=" + str(t_d4) + "_a4=" + str(a_d4[4]) + "_M4=" + str(m_d4[4]))
    check("wrong_control_ballot_detected_every_group", len(wd) == 5, "detected_in=" + (",".join(wd) or "none"))
    check("wrong_control_normalization_detected", norm_detected, "wrong=" + rs(wrong_norm) + "_log_route=" + rs(lr_c))
    check("wrong_control_molien_flagged", molien_order_flag and molien_feq_flag, "q*r=32_ne_48_and_12_ne_10")
    ex1 = PRED[("cube_6_faces_free", "Dic_3", "2T")]
    check("example_cube_Dic_3_vs_2T_beta^4_eq_1:64", ex1[1] == 4 and ex1[2] == (1, 64),
          "computed_power=" + str(ex1[1]) + "_value=" + rs(ex1[2]) + "_note=task_example_compared_not_forced", forced=False)
    ex2 = PRED[("cube_6_faces_free", "2T", "2I")]
    check("example_2T_vs_2I_first_at_beta^6", ex2[1] == 6, "computed_power=" + str(ex2[1]) + "_value=" + rs(ex2[2])
          + "_note=task_example_compared_not_forced", forced=False)

    failed = [c for c in checks if c[3] and not c[1]]
    rows_s = ["|".join(["SELF_CHECK", "name=" + nm, "result=" + ("PASS" if ok else "FAIL"), "forced=" + ("1" if forced else "0"),
                        "detail=" + det, "json=0"]) for nm, ok, det, forced in checks]

    head = "|".join([
        "TRANSFER_PREDICT_HEAD", "schema=" + SCHEMA, "date=" + DATE, "author=" + AUTHOR,
        "theorem=design_strength_t_of_G_transfers_to_lattice_log_Z_through_beta^t_with_sharp_beta^(t+1)_difference",
        "theorem_tag=NAMED_proof_written_in_PROOF_rows_machine_check_NOT_RUN",
        "measurement_status=NOT_RUN", "written_before_any_measurement=1", "configurations_enumerated=0",
        "group_elements_built=0", "link_variables_assigned=0", "sums_over_a_group=0",
        "lattice_objects_built=plaquette_link_incidence_only", "integer_only=1", "float_used=0",
        "rational_format=num:den_reduced_den_positive",
        "inputs=Molien_series(Dic_n:(1+t^(2n+2)):((1-t^4)(1-t^(2n))),2T:(1+t^12):((1-t^6)(1-t^8)),2O:(1+t^18):((1-t^8)(1-t^12)),"
        "2I:(1+t^30):((1-t^12)(1-t^20))),N_p",
        "moment_rule=M_n=sum_j_B(n,j)*a_(2j),B(n,j)=C(n,(n-2j):2)-C(n,(n-2j):2-1),Haar=Catalan(n:2)",
        "groups=Q8(Dic_2),Dic_3,2T,2O,2I,SU2_Haar", "pairs=" + str(len(pairs)),
        "lattices=cube_6_faces_free(N_p=6),torus_3x3_periodic(N_p=9),torus_4x4_periodic(N_p=16),control_torus_2x2_periodic(N_p=4)",
        "weight=exp(beta*sum_p_chi(U_p)*(1:2))", "coefficient_formula=N_p*(M_G1-M_G2)_(t+1):(2^(t+1)*(t+1)!)",
        "script=predict_transfer.py", "script_sha256=" + script_sha, "json=0"])
    rows_t = ["|".join(["THEOREM", "part=" + part, "needs=" + needs, "statement=" + text, "tag=NAMED", "json=0"])
              for part, needs, text in THEOREM]
    rows_k = ["|".join(["CORRECTION", "id=" + cid, "on=" + on, "was=" + was, "now=" + now, "consequence=" + cons, "tag=NAMED", "json=0"])
              for cid, on, was, now, cons in CORRECTIONS]
    rows_pf = ["|".join(["PROOF", "step=" + st, "name=" + nm, "text=" + text, "tag=NAMED", "json=0"]) for st, nm, text in PROOF]

    body = [head] + rows_t + rows_k + rows_pf + rows_g + rows_h + rows_l + rows_p + [row_ctl] + rows_c + [row_d4] + rows_w + rows_s
    total = len(body) + 1
    end = "|".join([
        "TRANSFER_PREDICT_END", "rows=" + str(total), "pair_rows=" + str(len(rows_p)), "control_pair_rows=" + str(len(rows_c)),
        "self_checks=" + str(len(checks)), "forced_failed=" + str(len(failed)),
        "status=" + ("ALL_FORCED_SELF_CHECKS_PASS" if not failed else "SELF_CHECK_FAILED"),
        "measurement_status=NOT_RUN", "json=0"])
    text = "\n".join(body + [end]) + "\n"
    for row in body + [end]:
        assert "\n" not in row and row.endswith("|json=0")
    data = text.encode("utf-8")
    out = os.path.join(here, OUT_NAME)
    with open(out, "wb") as fh:
        fh.write(data)
    sha = hashlib.sha256(data).hexdigest()
    with open(out + ".sha256", "wb") as fh:
        fh.write((sha + "  " + OUT_NAME + "\n").encode("ascii"))
    print("TRANSFER_PREDICT|out=" + OUT_NAME + "|sha256=" + sha + "|script_sha256=" + script_sha + "|rows=" + str(total)
          + "|forced_failed=" + str(len(failed)) + "|json=0")
    for c in checks:
        if not c[1]:
            print("TRANSFER_PREDICT_CHECK|name=" + c[0] + "|result=FAIL|forced=" + ("1" if c[3] else "0") + "|detail=" + c[2] + "|json=0")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
