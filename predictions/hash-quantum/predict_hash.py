# predict_hash.py -- the PREDICTION, written before any measurement is run.
#
# The owner's point (2026-11-09): the global hash reactions carry the quantum effects; a discrete
# universe shows the continuous feel of the quantum world only through that hash.
#
# What it predicts, for n = 1..8 qubits, N = 2^n, every split A|B with d_A = 2^k (k = 0..n), and
# t = 1..4, for phase states psi_f(x) = w^f(x) times 1:sqrt(N), w = exp(2 pi i : q):
#   HAAR      continuous Haar-random pure states on C^N
#   GLOBAL_q  f uniform over ALL functions {0,1}^n -> Z_q, q = 2,3,4,5    (the global hash)
#   LOCAL_q   f(x) = h_1(x_1) + ... + h_n(x_n) mod q, each h_i uniform     (a local hash)
#   QUAD_2    f(x) = sum_(i<j) B_ij x_i x_j + sum_i c_i x_i mod 2, every bit uniform
#   1. FRAME   F_t = E |<psi_f|psi_g>|^(2t), the Haar value, the exact difference, the design order
#   2. PURITY  E Tr rho_A^2 for every ensemble and split
#   3. ENTANGLEMENT_GAP at the balanced split k = floor(n:2)
#   4. CONVERGENCE of GLOBAL_3 to Haar, with exact limits from the count polynomials
#
# HARD RULE kept here: no state vector, density matrix or function table is built and nothing is
# summed over functions f. Every value comes from a closed counting formula. The independent checks
# count index tuples of positions x in 0..N-1 (never states, never functions) for small N only.
# Sources (NAMED): Welch 1974; Lubkin 1978; Page 1993; Ji-Liu-Song 2018; Brakerski-Shmueli 2019;
# Dehaene-De Moor 2003; Van den Nest 2010; MacWilliams 1969 (alternating matrices by rank);
# Landsberg 1893 (matrices over GF(q) by rank); MacWilliams-Sloane 1977 ch. 13 (Reed-Muller codes).
#
# Arithmetic: integers only. A rational is a reduced pair (num, den) with den > 0, reduced by gcd
# after an exact divisibility assert, using only the integer operators // and %.
# Run: python3 -I -B predict_hash.py   (writes PREDICTION-HASH.hbp and PREDICTION-HASH.hbp.sha256
# beside this script; exits 1 if any self-check fails)
import hashlib
import os
import sys
from math import gcd

N_MAX = 8
T_MAX = 4
QS = (2, 3, 4, 5)
BRUTE_FRAME_N_MAX = 4      # tuple-pair count by residue class, n <= 4 (N <= 16)
BRUTE_PURITY_N_MAX = 6     # index-quadruple survivor count, n <= 6
FORMAL_LOCAL_N_MAX = 5     # formal h-variable coefficient count, n <= 5
DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
SCHEMA = "ILLUSION_HASH_PREDICTION_V1"
OUT_NAME = "PREDICTION-HASH.hbp"

# ---------------------------------------------------------------- exact integers and rationals


def ipow(b, e):
    assert e >= 0
    r = 1
    for _ in range(e):
        r *= b
    return r


def exact_div(a, b):
    assert b != 0 and a % b == 0
    return a // b


FACT = [1]
for _i in range(1, 41):
    FACT.append(FACT[-1] * _i)


def binom(n, k):
    if k < 0 or n < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = exact_div(r * (n - i), i + 1)
    return r


def falling(n, m):
    r = 1
    for i in range(m):
        r *= n - i
    return r


def rising(n, m):
    r = 1
    for i in range(m):
        r *= n + i
    return r


def dfact_odd(t):
    # (2t-1)!! = 1*3*5*...*(2t-1)
    r = 1
    for i in range(t):
        r *= 2 * i + 1
    return r


def rat(n, d):
    assert d != 0
    if d < 0:
        n, d = -n, -d
    g = gcd(n, d)
    assert n % g == 0 and d % g == 0
    return (n // g, d // g)


ZERO = (0, 1)
ONE = (1, 1)


def ri(n):
    return (n, 1)


def radd(x, y):
    return rat(x[0] * y[1] + y[0] * x[1], x[1] * y[1])


def rsub(x, y):
    return rat(x[0] * y[1] - y[0] * x[1], x[1] * y[1])


def rmul(x, y):
    return rat(x[0] * y[0], x[1] * y[1])


def rdiv(x, y):
    assert y[0] != 0
    return rat(x[0] * y[1], x[1] * y[0])


def rlt(x, y):
    return x[0] * y[1] < y[0] * x[1]


def rle(x, y):
    return x[0] * y[1] <= y[0] * x[1]


def rs(x):
    return str(x[0]) + ":" + str(x[1])


def row(tag, fields):
    return "|".join([tag] + fields + ["json=0"])


# ---------------------------------------------------------------- HAAR (NAMED closed forms)


def haar_frame(N, t):
    return rat(1, binom(N + t - 1, t))


def haar_frame_route2(N, t):
    return rat(FACT[t], rising(N, t))


def haar_purity(dA, dB):
    return rat(dA + dB, dA * dB + 1)


def haar_purity_route2(dA, dB):
    # Tr[(SWAP_A x 1_B)(1 + SWAP)] : (N(N+1)) with Tr SWAP_A x 1_B = d_A d_B^2, Tr SWAP_A SWAP = d_A^2 d_B
    N = dA * dB
    return rat(dA * dB * dB + dA * dA * dB, N * (N + 1))


# ---------------------------------------------------------------- GLOBAL_q frame potential by types

_TYPES = {}


def allowed_pairs(t, q):
    out = []
    for a in range(t + 1):
        for b in range(t + 1):
            if (a, b) != (0, 0) and (a - b) % q == 0:
                out.append((a, b))
    return out


def pair_types(t, q):
    # every multiset of nonzero per-point count pairs (a, b), a = b mod q, sum a = t, sum b = t
    key = (t, q)
    if key in _TYPES:
        return _TYPES[key]
    pairs = allowed_pairs(t, q)
    out = []

    def rec(start, ra, rb, chosen):
        if ra == 0 and rb == 0:
            out.append(tuple(chosen))
            return
        for i in range(start, len(pairs)):
            a, b = pairs[i]
            if a <= ra and b <= rb:
                chosen.append(pairs[i])
                rec(i, ra - a, rb - b, chosen)
                chosen.pop()

    rec(0, t, t, [])
    _TYPES[key] = out
    return out


def type_weight(tp, t):
    # orderings of the two t-tuples inside one placement: (t!)^2 : prod a! b!
    den = 1
    for a, b in tp:
        den *= FACT[a] * FACT[b]
    coef = exact_div(FACT[t] * FACT[t], den)
    mult = {}
    for p in tp:
        mult[p] = mult.get(p, 0) + 1
    sym = 1
    for v in mult.values():
        sym *= FACT[v]
    return coef, sym, len(tp)


def global_count_types(N, t, q):
    # number of pairs of t-tuples in [N]^t x [N]^t whose per-point counts agree mod q
    total = 0
    for tp in pair_types(t, q):
        coef, sym, m = type_weight(tp, t)
        total += coef * exact_div(falling(N, m), sym)
    return total


def falling_poly(m):
    # coefficients of N(N-1)...(N-m+1), index = power of N
    c = [1]
    for j in range(m):
        new = [0] * (len(c) + 1)
        for i, v in enumerate(c):
            new[i + 1] += v
            new[i] -= j * v
        c = new
    return c


def global_count_poly(t, q):
    poly = [ZERO] * (t + 1)
    for tp in pair_types(t, q):
        coef, sym, m = type_weight(tp, t)
        fp = falling_poly(m)
        for i, v in enumerate(fp):
            poly[i] = radd(poly[i], rat(coef * v, sym))
    return poly


def poly_eval(poly, N):
    acc = ZERO
    p = 1
    for c in poly:
        acc = radd(acc, rmul(c, ri(p)))
        p *= N
    return acc


SMALLB = [[binom(a, b) for b in range(T_MAX + 1)] for a in range(T_MAX + 1)]


def global_count_dp(N, t, q):
    # route 2: walk the N points one by one; state (s, r) = x-slots and y-slots already placed.
    # q = 0 means continuous U(1) phases: the counts must agree exactly.
    pairs = []
    for a in range(t + 1):
        for b in range(t + 1):
            if (q == 0 and a == b) or (q > 0 and (a - b) % q == 0):
                pairs.append((a, b))
    dp = {(0, 0): 1}
    for _ in range(N):
        new = {}
        for (s, r), v in dp.items():
            for a, b in pairs:
                if s + a <= t and r + b <= t:
                    k2 = (s + a, r + b)
                    new[k2] = new.get(k2, 0) + v * SMALLB[t - s][a] * SMALLB[t - r][b]
        dp = new
    return dp.get((t, t), 0)


def global_count_brute(N, t, q):
    # route 3: list every t-tuple of positions, class it by its count vector mod q; pairs in the
    # same class are the survivors, so the count is the sum of squared class sizes
    classes = {}
    for idx in range(ipow(N, t)):
        cnt = [0] * N
        v = idx
        for _ in range(t):
            cnt[v % N] += 1
            v = v // N
        key = tuple(c % q for c in cnt)
        classes[key] = classes.get(key, 0) + 1
    return sum(s * s for s in classes.values())


def global_frame(N, t, q):
    return rat(global_count_types(N, t, q), ipow(N, 2 * t))


# ---------------------------------------------------------------- LOCAL_q frame potential


def local_single(t, q):
    return rat(global_count_types(2, t, q), ipow(4, t))


def local_single_route2(t, q):
    # per qubit |overlap|^2 = (2 + z + conj z) : 4, z = w^delta, delta uniform in Z_q
    total = 0
    for i in range(t + 1):
        for j in range(t + 1 - i):
            l = t - i - j
            if (j - l) % q == 0:
                total += exact_div(FACT[t], FACT[i] * FACT[j] * FACT[l]) * ipow(2, i)
    return rat(total, ipow(4, t))


def local_frame(n, t, q):
    one = local_single(t, q)
    r = ONE
    for _ in range(n):
        r = rmul(r, one)
    return r


# ---------------------------------------------------------------- QUAD_2 frame potential


def alt_rank_counts_rec(n):
    # T[r] = number of n x n alternating GF(2) matrices of rank 2r, by bordering with a vector v
    T = [1]
    for m in range(n):
        new = [0] * (len(T) + 1)
        for r, v in enumerate(T):
            new[r] += v * ipow(2, 2 * r)
            new[r + 1] += v * (ipow(2, m) - ipow(2, 2 * r))
        while len(new) > 1 and new[-1] == 0:
            new.pop()
        T = new
    return T


def alt_rank_counts_macwilliams(n):
    out = []
    s = 0
    while 2 * s <= n:
        num = 1
        den = 1
        for i in range(s):
            num *= ipow(2, 2 * i)
            den *= ipow(2, 2 * i + 2) - 1
        for i in range(2 * s):
            num *= ipow(2, n - i) - 1
        out.append(exact_div(num, den))
        s += 1
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def quad_frame(n, t):
    if t == 0:
        return ONE
    T = alt_rank_counts_rec(n)
    total = ZERO
    for r, v in enumerate(T):
        total = radd(total, rat(v, ipow(2, n + binom(n, 2) + 2 * r * (t - 1))))
    return total


def gauss_binom2(n, k):
    num = 1
    den = 1
    for i in range(k):
        num *= ipow(2, n - i) - 1
        den *= ipow(2, i + 1) - 1
    return exact_div(num, den)


def flats3(n):
    # affine 3-flats of GF(2)^n = weight-8 words of RM(n-3, n)
    if n < 3:
        return 0
    return ipow(2, n - 3) * gauss_binom2(n, 3)


def quad_frame_rm(n, t):
    assert t <= 4
    N = 1 << n
    c = global_count_types(N, t, 2)
    if t == 4:
        c += FACT[8] * flats3(n)
    return rat(c, ipow(N, 2 * t))


# ---------------------------------------------------------------- purity


def global_purity(dA, dB):
    surv = dA * dB * dB + dA * dA * dB - dA * dB
    return rat(surv, (dA * dB) * (dA * dB))


def global_purity_brute(n, k, q):
    # count index quadruples (a, a', b, b') whose four-corner net coefficients all vanish mod q
    dA = 1 << k
    dB = 1 << (n - k)
    surv = 0
    for a in range(dA):
        for a2 in range(dA):
            for b in range(dB):
                for b2 in range(dB):
                    net = {}
                    for pa, pb, s in ((a, b, 1), (a2, b, -1), (a2, b2, 1), (a, b2, -1)):
                        net[(pa, pb)] = net.get((pa, pb), 0) + s
                    if all(v % q == 0 for v in net.values()):
                        surv += 1
    N = dA * dB
    return rat(surv, N * N)


def local_purity_formal(n, k, q):
    # formal four-corner coefficient of every local variable h_i(v); a quadruple survives iff every
    # coefficient is 0 mod q. A = the high k bits of x, B = the low n-k bits.
    dA = 1 << k
    dB = 1 << (n - k)
    surv = 0
    for a in range(dA):
        for a2 in range(dA):
            for b in range(dB):
                for b2 in range(dB):
                    coeff = {}
                    for pa, pb, s in ((a, b, 1), (a2, b, -1), (a2, b2, 1), (a, b2, -1)):
                        x = pa * dB + pb
                        for i in range(n):
                            kk = (i, (x >> i) & 1)
                            coeff[kk] = coeff.get(kk, 0) + s
                    if all(v % q == 0 for v in coeff.values()):
                        surv += 1
    N = dA * dB
    return rat(surv, N * N)


def gf2_rank_counts_formula(k, m):
    out = []
    for r in range(min(k, m) + 1):
        num = 1
        den = 1
        for i in range(r):
            num *= (ipow(2, k) - ipow(2, i)) * (ipow(2, m) - ipow(2, i))
            den *= ipow(2, r) - ipow(2, i)
        out.append(exact_div(num, den))
    return out


def gf2_rank_counts_rec(k, m):
    # add rows one at a time: a new row in the current row space keeps rank r (2^r ways)
    R = [1]
    for _ in range(k):
        new = [0] * (len(R) + 1)
        for r, v in enumerate(R):
            new[r] += v * ipow(2, r)
            new[r + 1] += v * (ipow(2, m) - ipow(2, r))
        R = new
    assert all(v == 0 for v in R[min(k, m) + 1:])
    return R[:min(k, m) + 1]


def quad_purity(n, k):
    m = n - k
    R = gf2_rank_counts_formula(k, m)
    total = ZERO
    for r, v in enumerate(R):
        total = radd(total, rat(v, ipow(2, k * m + r)))
    return total


def quad_purity_kernel(n, k):
    dA = 1 << k
    dB = 1 << (n - k)
    return rat(dA + dB - 1, dA * dB)


ENSEMBLES = ["HAAR"] + ["GLOBAL_" + str(q) for q in QS] + ["LOCAL_" + str(q) for q in QS] + ["QUAD_2"]


def ens_q(e):
    return int(e.split("_")[1])


def frame(e, n, t):
    N = 1 << n
    if e == "HAAR":
        return haar_frame(N, t)
    if e.startswith("GLOBAL_"):
        return global_frame(N, t, ens_q(e))
    if e.startswith("LOCAL_"):
        return local_frame(n, t, ens_q(e))
    if e == "QUAD_2":
        return quad_frame(n, t)
    raise ValueError(e)


def purity(e, n, k):
    dA = 1 << k
    dB = 1 << (n - k)
    if e == "HAAR":
        return haar_purity(dA, dB)
    if e.startswith("GLOBAL_"):
        return global_purity(dA, dB)
    if e.startswith("LOCAL_"):
        return ONE
    if e == "QUAD_2":
        return quad_purity(n, k)
    raise ValueError(e)


# ---------------------------------------------------------------- derivation text

DERIVATIONS = [
    ("D01_HAAR_FRAME", "NAMED", "F_t(HAAR)=1:binom(N+t-1,t), the inverse dimension of Sym^t(C^N); every ensemble has F_t>=F_t(HAAR) (Welch bound) with equality iff it is an exact t-design; design order = largest t with equality"),
    ("D02_HAAR_PURITY", "NAMED", "E Tr rho_A^2=(d_A+d_B):(N+1) (Lubkin 1978, Page 1993); route 2: Tr[(SWAP_A x 1_B)(1+SWAP)]:(N(N+1))=(d_A d_B^2+d_A^2 d_B):(N(N+1))"),
    ("D03_GLOBAL_FRAME", "DERIVED", "<psi_f,psi_g>=N^-1 sum_x w^(g(x)-f(x)); h=g-f is uniform over all functions so z_x=w^h(x) are iid uniform q-th roots of unity and F_t=N^(-2t) E abs(sum_x z_x)^(2t); E z^m=[m=0 mod q] so E abs(S)^(2t)=number of pairs of t-tuples X,Y in [N]^t with count_x(X)=count_x(Y) mod q at every x; grouped by the multiset (type) of nonzero per-point pairs (a,b) with a=b mod q: count=sum_types N!:((N-m)! prod mu!) * (t!)^2:prod(a! b!), m=parts, mu=multiplicities; a polynomial in N of degree t; route 2 a point-by-point DP, route 3 a residue-class count of t-tuples for N<=16"),
    ("D04_GLOBAL_Q_LADDER", "DERIVED", "for q>t the condition a=b mod q with 0<=a,b<=t forces a=b, so GLOBAL_q equals uniform continuous U(1) phases (PHASE_U1) for every t<q; for t>=q the types with (q,0) and (0,q) add strictly positive terms when N>=2; leading coefficient: each part has a+b>=2, so t parts need weight 2 each: only (1,1) for q>=3 giving c_t=t!, and also (2,0),(0,2) for q=2 giving c_t=(2t-1)!!"),
    ("D05_LIMITS", "DERIVED", "F_t(HAAR)=t!:(N(N+1)..(N+t-1))=t! N^-t (1-binom(t,2):N+O(N^-2)); with count=c_t N^t+c_(t-1) N^(t-1)+..: F:F_HAAR -> c_t:t!; when c_t=t!: (F-F_HAAR) N^(t+1) -> L_t=c_(t-1)+t! binom(t,2) and (F:F_HAAR-1) N -> L_t:t!; so (F-F_HAAR) N^2 -> 0 for t>=2"),
    ("D06_GLOBAL_PURITY", "DERIVED", "x=(a,b); Tr rho_A^2=N^-2 sum_(a,a2,b,b2) w^(f(a,b)-f(a2,b)+f(a2,b2)-f(a,b2)); over uniform f a term has expectation 1 iff the net coefficient at every point is 0 mod q, else 0; a=a2 or b=b2: corners cancel in pairs; a!=a2 and b!=b2: four distinct points, each net coefficient +1 or -1, and +-1 is nonzero mod q for every q>=2, so q=2 brings no extra cancellation (f(a,b)+f(a2,b)+f(a,b2)+f(a2,b2) mod 2 is a uniform bit, E(-1)^bit=0); survivors d_A d_B^2+d_A^2 d_B-d_A d_B; E Tr rho_A^2=(d_A+d_B-1):N for every q>=2; minus HAAR = -(d_A-1)(d_B-1):(N(N+1)); route 2 counts index quadruples for n<=6"),
    ("D07_LOCAL_FRAME", "DERIVED", "psi_f = tensor_i (w^h_i(0) ket0 + w^h_i(1) ket1):sqrt2; overlaps factor over qubits with independent h_i so F_t(LOCAL_q,n)=F_t(GLOBAL_q,1)^n; route 2 per qubit: abs(overlap)^2=(2+z+conj z):4, z=w^delta, delta uniform in Z_q, factor=4^-t sum_(i+j+l=t) t!:(i!j!l!) 2^i [j=l mod q]"),
    ("D08_LOCAL_PURITY", "DERIVED", "each psi_f is a product state over qubits so rho_A is pure: Tr rho_A^2=1 for every f and every split; PREDICTION min=mean=max=1; route 2: f(a,b)=H_A(a)+H_B(b) so the four-corner exponent vanishes identically: the formal coefficient of every variable h_i(v) is 0 for all N^2 quadruples (counted for n<=5, every q, every split)"),
    ("D09_QUAD_FRAME", "DERIVED", "h=f+g is again a uniform QUAD_2 form; S=sum_x (-1)^h(x); A=B+B^T is a uniform alternating GF(2) matrix; abs(S)^2=N sum_(u in rad A) (-1)^h(u); h is additive on rad A and uniform c makes h restricted to rad A a uniform linear functional, so abs(S:N)^2=2^(-2r) with probability 2^-(n-2r) and 0 otherwise, rank A=2r; F_t=2^-n sum_r P(rank A=2r) 2^(2r(1-t)) for t>=1; alternating counts T(n,2r) by bordering T(n+1,2r)=4^r T(n,2r)+(2^n-4^(r-1)) T(n,2r-2) and by the MacWilliams 1969 product"),
    ("D10_QUAD_VS_GLOBAL2", "DERIVED", "route 2: E(-1)^(sum of h over 2t points)=1 iff the odd-multiplicity support lies in RM(2,n)^perp=RM(n-3,n), minimum weight 8, weight-8 words = affine 3-flats (MacWilliams-Sloane 1977 ch.13); so F_t(QUAD_2)=F_t(GLOBAL_2) for t<=3 and F_4(QUAD_2)-F_4(GLOBAL_2)=8! 2^(n-3) [n,3]_2 : N^8 for n>=3, 0 for n<3"),
    ("D11_QUAD_PURITY", "DERIVED", "four-corner exponent=(a+a2)^T M (b+b2) with M=B_AB the k x (n-k) cross block (Q_A, Q_B and c cancel) so Tr rho_A^2=2^-rank_GF2(M) for every state; mean from rank counts prod_(i<r) (2^k-2^i)(2^(n-k)-2^i):(2^r-2^i) (Landsberg 1893) and a row recurrence; route 3: #(u: u^T M=0)=2^(k-r), its mean is 1+(2^k-1)2^-(n-k), so E 2^-r=2^-k+2^-(n-k)-2^-n=(d_A+d_B-1):N, the GLOBAL_q mean"),
    ("D12_SCOPE", "DERIVED", "splits are the first k qubits against the rest; every ensemble here is invariant under qubit permutations, so any k-subset gives the same mean; GLOBAL_q purity does not depend on q; F_t does not depend on the split; F_0=1 by definition (abs(.)^0=1)"),
]


# ---------------------------------------------------------------- main


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    checks = []

    def check(name, ok, detail):
        checks.append((name, ok, detail))

    rows_d = [row("DERIVATION", ["id=" + i, "tag=" + tg, "text=" + tx]) for i, tg, tx in DERIVATIONS]

    # ------------------------------------------------ frame potentials
    F = {}
    for e in ENSEMBLES:
        for n in range(1, N_MAX + 1):
            for t in range(0, T_MAX + 1):
                F[(e, n, t)] = frame(e, n, t)

    rows_f = []
    welch_bad = []
    for e in ENSEMBLES:
        for n in range(1, N_MAX + 1):
            N = 1 << n
            for t in range(1, T_MAX + 1):
                f = F[(e, n, t)]
                h = F[("HAAR", n, t)]
                d = rsub(f, h)
                if rlt(f, h):
                    welch_bad.append(e + "_n" + str(n) + "_t" + str(t))
                rows_f.append(row("FRAME", [
                    "ensemble=" + e, "n=" + str(n), "N=" + str(N), "t=" + str(t),
                    "F=" + rs(f), "F_haar=" + rs(h), "diff=" + rs(d), "ratio=" + rs(rdiv(f, h)),
                    "equal_haar=" + ("1" if d == ZERO else "0"), "tag=PREDICTION"]))

    # design order
    rows_des = []
    summary = {}
    for e in ENSEMBLES:
        orders = []
        deps = []
        for n in range(1, N_MAX + 1):
            order = 0
            dep = 0
            for t in range(1, T_MAX + 1):
                if F[(e, n, t)] == F[("HAAR", n, t)]:
                    if dep == 0:
                        order = t
                else:
                    if dep == 0:
                        dep = t
            orders.append(order)
            deps.append(dep)
            fd = rs(F[(e, n, dep)]) if dep else "none"
            hd = rs(F[("HAAR", n, dep)]) if dep else "none"
            rows_des.append(row("DESIGN", [
                "ensemble=" + e, "n=" + str(n), "design_order=" + str(order),
                "first_departure_t=" + (str(dep) if dep else "none_through_4"),
                "F_at_departure=" + fd, "F_haar_at_departure=" + hd, "tag=PREDICTION"]))
        summary[e] = (orders, deps)
        rows_des.append(row("DESIGN_SUMMARY", [
            "ensemble=" + e, "n_range=1..8", "design_order_by_n=" + ",".join(str(v) for v in orders),
            "first_departure_by_n=" + ",".join(str(v) for v in deps), "tag=PREDICTION"]))

    # q ladder against continuous U(1) phases, and QUAD_2 against GLOBAL_2
    rows_q = []
    ladder_bad = []
    quad_bad = []
    for n in range(1, N_MAX + 1):
        N = 1 << n
        for t in range(1, T_MAX + 1):
            u1 = rat(global_count_dp(N, t, 0), ipow(N, 2 * t))
            eq = []
            for q in QS:
                same = F[("GLOBAL_" + str(q), n, t)] == u1
                eq.append(str(q) + ":" + ("1" if same else "0"))
                if same != (t < q):
                    ladder_bad.append("n" + str(n) + "t" + str(t) + "q" + str(q))
            rows_q.append(row("Q_LADDER", [
                "n=" + str(n), "t=" + str(t), "F_phase_u1=" + rs(u1),
                "F_by_q=" + ",".join(str(q) + "=" + rs(F[("GLOBAL_" + str(q), n, t)]) for q in QS),
                "equals_u1_by_q=" + ",".join(eq), "rule=equal_iff_t<q", "tag=PREDICTION"]))
            dq = rsub(F[("QUAD_2", n, t)], F[("GLOBAL_2", n, t)])
            pred = ZERO if t <= 3 else rat(FACT[8] * flats3(n), ipow(N, 8))
            if dq != pred:
                quad_bad.append("n" + str(n) + "t" + str(t))
            rows_q.append(row("QUAD_VS_GLOBAL2", [
                "n=" + str(n), "t=" + str(t), "F_quad=" + rs(F[("QUAD_2", n, t)]),
                "F_global2=" + rs(F[("GLOBAL_2", n, t)]), "diff=" + rs(dq),
                "rm_route_diff=" + rs(pred), "flats3=" + str(flats3(n)), "agree=" + ("1" if dq == pred else "0"),
                "tag=PREDICTION"]))

    # ------------------------------------------------ purity
    rows_p = []
    bound_bad = []
    k0_bad = []
    sym_bad = []
    P = {}
    for e in ENSEMBLES:
        for n in range(1, N_MAX + 1):
            for k in range(0, n + 1):
                P[(e, n, k)] = purity(e, n, k)
    for e in ENSEMBLES:
        for n in range(1, N_MAX + 1):
            N = 1 << n
            for k in range(0, n + 1):
                dA = 1 << k
                dB = 1 << (n - k)
                p = P[(e, n, k)]
                h = P[("HAAR", n, k)]
                lo = rat(1, min(dA, dB))
                if not (rle(lo, p) and rle(p, ONE)):
                    bound_bad.append(e + "_n" + str(n) + "_k" + str(k))
                if (k == 0 or k == n) and p != ONE:
                    k0_bad.append(e + "_n" + str(n) + "_k" + str(k))
                if p != P[(e, n, n - k)]:
                    sym_bad.append(e + "_n" + str(n) + "_k" + str(k))
                d = rsub(p, h)
                rows_p.append(row("PURITY", [
                    "ensemble=" + e, "n=" + str(n), "k=" + str(k), "d_A=" + str(dA), "d_B=" + str(dB),
                    "purity=" + rs(p), "purity_haar=" + rs(h), "diff=" + rs(d),
                    "diff_times_N=" + rs(rmul(d, ri(N))), "max_mixed=" + rs(lo),
                    "purity_times_min_d=" + rs(rmul(p, ri(min(dA, dB)))), "tag=PREDICTION"]))

    rows_rank = []
    rank_bad = []
    for n in range(1, N_MAX + 1):
        for k in range(0, n + 1):
            m = n - k
            rf = gf2_rank_counts_formula(k, m)
            rr = gf2_rank_counts_rec(k, m)
            ok = rf == rr and sum(rf) == ipow(2, k * m)
            if not ok:
                rank_bad.append("n" + str(n) + "k" + str(k))
            rows_rank.append(row("QUAD_RANK", [
                "n=" + str(n), "k=" + str(k), "cross_block=" + str(k) + "x" + str(m),
                "count_by_rank=" + ",".join(str(v) for v in rf), "total=" + str(sum(rf)),
                "two_pow_km=" + str(ipow(2, k * m)), "mean_purity=" + rs(quad_purity(n, k)), "tag=PREDICTION"]))
        T = alt_rank_counts_rec(n)
        rows_rank.append(row("QUAD_ALT_RANK", [
            "n=" + str(n), "count_by_rank_2r=" + ",".join(str(v) for v in T), "total=" + str(sum(T)),
            "two_pow_binom_n_2=" + str(ipow(2, binom(n, 2))), "tag=PREDICTION"]))

    # ------------------------------------------------ entanglement gap at the balanced split
    rows_g = []
    for n in range(1, N_MAX + 1):
        k = n // 2
        dA = 1 << k
        mm = rat(1, dA)
        h = P[("HAAR", n, k)]
        for e in ENSEMBLES:
            p = P[(e, n, k)]
            rows_g.append(row("ENTANGLEMENT_GAP", [
                "n=" + str(n), "k=" + str(k), "d_A=" + str(dA), "d_B=" + str(1 << (n - k)), "ensemble=" + e,
                "purity=" + rs(p), "purity_haar=" + rs(h), "minus_haar=" + rs(rsub(p, h)),
                "max_mixed=" + rs(mm), "minus_max_mixed=" + rs(rsub(p, mm)),
                "ratio_to_haar=" + rs(rdiv(p, h)), "ratio_to_max_mixed=" + rs(rdiv(p, mm)), "tag=PREDICTION"]))

    # ------------------------------------------------ polynomials, limits, convergence
    rows_c = []
    limits = {}
    lim_bad = []
    for qq in list(QS) + [0]:
        for t in range(1, T_MAX + 1):
            if qq == 0:
                poly = global_count_poly(t, t + 1)
                name = "PHASE_U1"
            else:
                poly = global_count_poly(t, qq)
                name = "GLOBAL_" + str(qq)
            ct = poly[t]
            ctm = poly[t - 1]
            ratio_lim = rdiv(ct, ri(FACT[t]))
            if ct == ri(FACT[t]):
                L = radd(ctm, ri(FACT[t] * binom(t, 2)))
                lstr = rs(L)
                rel = rs(rdiv(L, ri(FACT[t])))
            else:
                L = None
                lstr = "diverges_ratio_limit_ne_1"
                rel = "diverges"
            limits[(name, t)] = L
            rows_c.append(row("COUNT_POLY", [
                "ensemble=" + name, "t=" + str(t), "count_coeffs_N^0..N^t=" + ",".join(rs(c) for c in poly),
                "lead=" + rs(ct), "t_fact=" + str(FACT[t]), "double_fact_2t-1=" + str(dfact_odd(t)),
                "ratio_limit=" + rs(ratio_lim), "diff_times_N^(t+1)_limit=" + lstr, "rel_times_N_limit=" + rel,
                "tag=PREDICTION"]))
            # the limit must be visible at N = 2^40: abs(diff N^(t+1) - L) * N <= 10^4
            if L is not None:
                NB = 1 << 40
                fb = rat(global_count_types(NB, t, qq if qq else t + 1), ipow(NB, 2 * t))
                x = rsub(rmul(rsub(fb, haar_frame(NB, t)), ri(ipow(NB, t + 1))), L)
                if abs(x[0]) * NB > 10000 * x[1]:
                    lim_bad.append(name + "_t" + str(t))
            for n in range(1, N_MAX + 1):
                N = 1 << n
                if not (poly_eval(poly, N) == ri(global_count_types(N, t, qq if qq else t + 1))):
                    lim_bad.append(name + "_poly_n" + str(n) + "_t" + str(t))

    for t in (2, 3):
        for n in range(1, N_MAX + 1):
            N = 1 << n
            f = F[("GLOBAL_3", n, t)]
            h = F[("HAAR", n, t)]
            d = rsub(f, h)
            rows_c.append(row("CONVERGENCE", [
                "ensemble=GLOBAL_3", "quantity=F_t", "n=" + str(n), "N=" + str(N), "t=" + str(t),
                "diff=" + rs(d), "diff_times_N^2=" + rs(rmul(d, ri(N * N))),
                "diff_times_N^3=" + rs(rmul(d, ri(N * N * N))),
                "diff_times_N^(t+1)=" + rs(rmul(d, ri(ipow(N, t + 1)))),
                "rel_diff_times_N=" + rs(rmul(rsub(rdiv(f, h), ONE), ri(N))),
                "limit_diff_times_N^2=0:1", "limit_diff_times_N^(t+1)=" + rs(limits[("GLOBAL_3", t)]),
                "limit_rel_diff_times_N=" + rs(rdiv(limits[("GLOBAL_3", t)], ri(FACT[t]))), "tag=PREDICTION"]))
    for n in range(1, N_MAX + 1):
        N = 1 << n
        k = n // 2
        dA = 1 << k
        dB = 1 << (n - k)
        d = rsub(P[("GLOBAL_3", n, k)], P[("HAAR", n, k)])
        closed = rat(-(dA - 1) * (dB - 1), N + 1)
        rows_c.append(row("CONVERGENCE", [
            "ensemble=GLOBAL_3", "quantity=purity_balanced", "n=" + str(n), "N=" + str(N), "k=" + str(k),
            "diff=" + rs(d), "diff_times_N=" + rs(rmul(d, ri(N))), "closed_form=-(d_A-1)(d_B-1):(N+1)=" + rs(closed),
            "agree=" + ("1" if rmul(d, ri(N)) == closed else "0"), "limit_diff_times_N=-1:1", "tag=PREDICTION"]))
        if rmul(d, ri(N)) != closed:
            lim_bad.append("purity_closed_n" + str(n))

    # LOCAL growth of F:F_HAAR, per qubit base count_(N=2)(t):2^t
    for qq in QS:
        for t in range(2, T_MAX + 1):
            base = rat(global_count_types(2, t, qq), ipow(2, t))
            rows_c.append(row("LOCAL_GROWTH", [
                "ensemble=LOCAL_" + str(qq), "t=" + str(t), "per_qubit_base=" + rs(base),
                "ratio_n8=" + rs(rdiv(F[("LOCAL_" + str(qq), 8, t)], F[("HAAR", 8, t)])),
                "rule=F:F_HAAR>=base^n:t!_and_base>1_so_ratio_grows_without_bound", "tag=PREDICTION"]))
            if not rlt(ONE, base):
                lim_bad.append("local_base_q" + str(qq) + "_t" + str(t))

    # ------------------------------------------------ self-checks
    f0_bad = [e + "_n" + str(n) for e in ENSEMBLES for n in range(1, N_MAX + 1) if F[(e, n, 0)] != ONE]
    f0_bad += ["QUAD_RM_n" + str(n) for n in range(1, N_MAX + 1) if quad_frame_rm(n, 0) != ONE]
    check("F0_eq_1_all_ensembles", not f0_bad, "bad=" + (",".join(f0_bad) or "none"))
    f1_bad = [e + "_n" + str(n) for e in ENSEMBLES for n in range(1, N_MAX + 1) if F[(e, n, 1)] != rat(1, 1 << n)]
    check("F1_eq_1:N_all_ensembles", not f1_bad, "bad=" + (",".join(f1_bad) or "none"))

    haar_bad = [str(n) + "_" + str(t) for n in range(1, N_MAX + 1) for t in range(T_MAX + 1)
                if haar_frame(1 << n, t) != haar_frame_route2(1 << n, t)]
    haar_bad += ["p" + str(n) + "_" + str(k) for n in range(1, N_MAX + 1) for k in range(n + 1)
                 if haar_purity(1 << k, 1 << (n - k)) != haar_purity_route2(1 << k, 1 << (n - k))]
    check("haar_two_routes", not haar_bad, "frame_binom_vs_rising_and_purity_closed_vs_swap|bad=" + (",".join(haar_bad) or "none"))

    dp_bad = []
    br_bad = []
    for q in QS:
        for n in range(1, N_MAX + 1):
            N = 1 << n
            for t in range(T_MAX + 1):
                c = global_count_types(N, t, q)
                if c != global_count_dp(N, t, q):
                    dp_bad.append("q" + str(q) + "n" + str(n) + "t" + str(t))
                if n <= BRUTE_FRAME_N_MAX and c != global_count_brute(N, t, q):
                    br_bad.append("q" + str(q) + "n" + str(n) + "t" + str(t))
    check("global_frame_types_vs_dp", not dp_bad, "q=2..5|n=1..8|t=0..4|bad=" + (",".join(dp_bad) or "none"))
    check("global_frame_types_vs_tuple_count", not br_bad, "q=2..5|n=1..4|t=0..4|bad=" + (",".join(br_bad) or "none"))

    loc_bad = [str(q) + "_" + str(t) for q in QS for t in range(T_MAX + 1) if local_single(t, q) != local_single_route2(t, q)]
    loc_bad += ["n1_" + str(q) + "_" + str(t) for q in QS for t in range(T_MAX + 1)
                if F[("LOCAL_" + str(q), 1, t)] != F[("GLOBAL_" + str(q), 1, t)]]
    check("local_single_qubit_two_routes", not loc_bad, "types_at_N=2_vs_trinomial|bad=" + (",".join(loc_bad) or "none"))

    alt_bad = [str(n) for n in range(0, N_MAX + 1)
               if alt_rank_counts_rec(n) != alt_rank_counts_macwilliams(n) or sum(alt_rank_counts_rec(n)) != ipow(2, binom(n, 2))]
    check("alternating_rank_counts_recurrence_vs_macwilliams_sum_2^binom(n,2)", not alt_bad, "n=0..8|bad=" + (",".join(alt_bad) or "none"))
    rmq_bad = [str(n) + "_" + str(t) for n in range(1, N_MAX + 1) for t in range(1, T_MAX + 1)
               if quad_frame(n, t) != quad_frame_rm(n, t)]
    check("quad_frame_alt_rank_vs_reed_muller", not rmq_bad, "n=1..8|t=1..4|bad=" + (",".join(rmq_bad) or "none"))
    check("quad_vs_global2_rule", not quad_bad, "bad=" + (",".join(quad_bad) or "none"))
    check("global_q_eq_u1_iff_t<q", not ladder_bad, "bad=" + (",".join(ladder_bad) or "none"))
    check("welch_bound_F_ge_haar", not welch_bad, "bad=" + (",".join(welch_bad) or "none"))

    check("purity_k0_and_kn_eq_1_all_ensembles", not k0_bad, "bad=" + (",".join(k0_bad) or "none"))
    check("purity_between_1:min(d_A,d_B)_and_1", not bound_bad, "bad=" + (",".join(bound_bad) or "none"))
    check("purity_symmetric_k_vs_n-k", not sym_bad, "bad=" + (",".join(sym_bad) or "none"))
    check("gf2_rank_counts_formula_vs_recurrence_sum_2^(k(n-k))", not rank_bad, "n=1..8|bad=" + (",".join(rank_bad) or "none"))
    qk_bad = [str(n) + "_" + str(k) for n in range(1, N_MAX + 1) for k in range(n + 1)
              if quad_purity(n, k) != quad_purity_kernel(n, k)]
    check("quad_purity_rank_mean_vs_kernel_identity", not qk_bad, "bad=" + (",".join(qk_bad) or "none"))
    gp_bad = []
    for q in QS:
        for n in range(1, BRUTE_PURITY_N_MAX + 1):
            for k in range(n + 1):
                if global_purity_brute(n, k, q) != global_purity(1 << k, 1 << (n - k)):
                    gp_bad.append("q" + str(q) + "n" + str(n) + "k" + str(k))
    check("global_purity_four_corner_vs_quadruple_count", not gp_bad,
          "q=2..5|n=1..6|all_k|q2_extra_cancellation_found=" + ("0" if not gp_bad else "1") + "|bad=" + (",".join(gp_bad) or "none"))
    lp_bad = []
    for q in QS:
        for n in range(1, FORMAL_LOCAL_N_MAX + 1):
            for k in range(n + 1):
                if local_purity_formal(n, k, q) != ONE:
                    lp_bad.append("q" + str(q) + "n" + str(n) + "k" + str(k))
    check("local_purity_formal_four_corner_eq_1", not lp_bad, "q=2..5|n=1..5|all_k|bad=" + (",".join(lp_bad) or "none"))
    check("limits_and_poly_consistent", not lim_bad, "N=2^40_bound_and_poly_eval_n=1..8|bad=" + (",".join(lim_bad) or "none"))

    # deliberately wrong control: claim LOCAL_2 purity at n=2, k=1 is 3:4 (below 1)
    claimed = rat(3, 4)
    rederived = local_purity_formal(2, 1, 2)
    caught = claimed != rederived and rlt(claimed, ONE) and rederived == ONE
    row_ctl = row("CONTROL", ["name=wrong_local_purity_below_1", "ensemble=LOCAL_2", "n=2", "k=1",
                              "claimed=" + rs(claimed), "independent_rederivation=" + rs(rederived),
                              "route=formal_four_corner_h_variable_count", "caught=" + ("1" if caught else "0")])
    check("control_wrong_local_purity_caught", caught, "claimed=" + rs(claimed) + "|rederived=" + rs(rederived))

    failed = sum(1 for c in checks if not c[1])
    rows_s = [row("SELF_CHECK", ["name=" + nm, "result=" + ("PASS" if ok else "FAIL"), "detail=" + det])
              for nm, ok, det in checks]

    head = row("ILLUSION_HASH_PREDICT_HEAD", [
        "schema=" + SCHEMA, "date=" + DATE, "author=" + AUTHOR,
        "named_sources=Welch_1974;Lubkin_1978;Page_1993;Ji_Liu_Song_2018;Brakerski_Shmueli_2019;Dehaene_De_Moor_2003;Van_den_Nest_2010;MacWilliams_1969;Landsberg_1893;MacWilliams_Sloane_1977_ch13",
        "formula_tag=NAMED_for_HAAR_DERIVED_for_the_rest", "measurement_status=NOT_RUN", "written_before_any_measurement=1",
        "states_built=0", "density_matrices_built=0", "function_tables_built=0", "sums_over_functions=0",
        "index_tuple_checks=frame_n<=4_purity_n<=6_local_formal_n<=5", "integer_only=1",
        "rational_format=num:den_reduced_den_positive", "ensembles=" + ",".join(ENSEMBLES),
        "n_range=1..8", "k_range=0..n", "t_range=1..4", "split=first_k_qubits_vs_rest",
        "script=predict_hash.py", "script_sha256=" + script_sha])

    body = [head] + rows_d + rows_f + rows_des + rows_q + rows_p + rows_rank + rows_g + rows_c + [row_ctl] + rows_s
    status = "ALL_SELF_CHECKS_PASS" if failed == 0 else "SELF_CHECK_FAILED"
    end = row("ILLUSION_HASH_PREDICT_END", ["rows=" + str(len(body) + 1), "self_checks=" + str(len(checks)),
                                             "failed=" + str(failed), "status=" + status, "measurement_status=NOT_RUN"])
    text = "\n".join(body + [end]) + "\n"
    data = text.encode("ascii")
    out = os.path.join(here, OUT_NAME)
    with open(out, "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(out + ".sha256", "wb") as fh:
        fh.write((digest + "  " + OUT_NAME + "\n").encode("ascii"))

    for nm, ok, det in checks:
        if not ok:
            print("SELF_CHECK_FAIL|name=" + nm + "|" + det + "|json=0")
    print("PREDICT_DONE|rows=" + str(len(body) + 1) + "|self_checks=" + str(len(checks)) + "|failed=" + str(failed)
          + "|script_sha256=" + script_sha + "|hbp_sha256=" + digest + "|json=0")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
