#!/usr/bin/env python3
"""m3_towers.py -- finite gauge towers in SU(3) and SU(2), the Jordan limit
theorem, the icosians over Z[phi], and G2-instantons from 4D ASD 2-forms.

Seat: Opus 5.5 (pid 8467a937cba309f7).  E=0: compute only, no network.
Integer only: no float, no true division, no fractions or decimal modules,
math used only for gcd.  Complex numbers enter as monomials over the ring
Z[x] mod (x^n + 1) (x a primitive 2n-th root of unity) and quaternion
coordinates as integer pairs over Z[phi] or Z[sqrt 2].

Writes M3-TOWERS.hbp and M3-TOWERS.hbp.sha256 next to this file.
Exit code 1 on any failed check and on any broken control that passes.

Items:
 1  Delta(6n^2) in SU(3), n = 1..7, as monomial matrices
 2  Jordan bound: exhaustive search for the largest abelian normal subgroup
 3  limit theorem for finite subgroups (NAMED, argued step by step in rows)
 4  SU(2): 2I (120 icosians) over Z[phi], 2T, 2O over Z[sqrt 2], Dic_m towers
 5  G2: phi(a,b,c) = <a, bc>, Hodge star, T = *(phi ^ .), Lambda2 split,
    der(O), stabiliser of phi, associative and coassociative planes, ASD
 6  broken controls (one sign of phi flipped, bad det, bad icosian, hole)
"""
import hashlib
import io
import itertools
import math
import os
import sys
import tokenize

SELF = os.path.abspath(__file__)
HERE = os.path.dirname(SELF)
OUT_NAME = "M3-TOWERS.hbp"
OUT = os.path.join(HERE, OUT_NAME)
SIDE = os.path.join(HERE, OUT_NAME + ".sha256")

ROWS = []
FAILS = []
COUNT = {"check": 0, "control": 0, "pass": 0}


# ----------------------------------------------------------------- rows
def clean(v):
    s = str(v)
    for ch in ("|", "\n", "\r"):
        s = s.replace(ch, " ")
    return s


def row(rtag, **kv):
    parts = [rtag]
    for k, v in kv.items():
        parts.append(k + "=" + clean(v))
    parts.append("json=0")
    line = "|".join(parts)
    ROWS.append(line)
    print(line)


def check(name, ok, tag="MEASURED", **kv):
    COUNT["check"] += 1
    if ok:
        COUNT["pass"] += 1
    else:
        FAILS.append(name)
    row("CHECK", name=name, status=("PASS" if ok else "FAIL"), tag=tag, **kv)
    return ok


def control(name, observed_ok, **kv):
    """A broken control must FAIL the check it is fed to."""
    COUNT["control"] += 1
    ok = not observed_ok
    if ok:
        COUNT["pass"] += 1
    else:
        FAILS.append("control_" + name)
    row("CONTROL", name=name, expected="FAIL",
        observed=("PASS" if observed_ok else "FAIL"),
        status=("PASS" if ok else "FAIL"), tag="MEASURED", **kv)
    return ok


# ------------------------------------------------------------ self lint
def self_lint():
    with open(SELF, encoding="utf-8") as fh:
        src = fh.read()
    bad = []
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    banned = {"float", "fractions", "decimal", "Fraction", "Decimal",
              "complex", "cmath", "numpy", "truediv"}
    for i, t in enumerate(toks):
        if t.type == tokenize.OP and t.string in ("/", "/="):
            bad.append("true_div_line%d" % t.start[0])
        if t.type == tokenize.NUMBER:
            s = t.string.lower()
            if not s.startswith(("0x", "0o", "0b")) and (
                    "." in s or "e" in s or "j" in s):
                bad.append("float_literal_line%d" % t.start[0])
        if t.type == tokenize.NAME and t.string in banned:
            bad.append("banned_%s_line%d" % (t.string, t.start[0]))
        if (t.type == tokenize.NAME and t.string == "math"
                and i + 2 < len(toks) and toks[i + 1].string == "."):
            if toks[i + 2].string not in ("gcd", "isqrt"):
                bad.append("math_%s_line%d" % (toks[i + 2].string, t.start[0]))
    return src, bad


# ------------------------------------------------------ integer algebra
def perm_parity(seq):
    inv = 0
    L = len(seq)
    for i in range(L):
        si = seq[i]
        for j in range(i + 1, L):
            if si > seq[j]:
                inv += 1
    return inv & 1


def perm_sign(seq):
    return -1 if perm_parity(seq) else 1


def primitive(v):
    g = 0
    for x in v:
        g = math.gcd(g, x)
    if g > 1:
        return [x // g for x in v]
    return list(v)


def rref_int(rows, ncols):
    """Fraction-free Gauss-Jordan over Z. Returns (rows, pivot columns)."""
    M = [list(r) for r in rows if any(r)]
    piv = []
    r = 0
    nrows = len(M)
    for c in range(ncols):
        if r >= nrows:
            break
        p = -1
        for i in range(r, nrows):
            if M[i][c] != 0:
                p = i
                break
        if p < 0:
            continue
        M[r], M[p] = M[p], M[r]
        pr = M[r]
        a = pr[c]
        for i in range(nrows):
            if i != r and M[i][c] != 0:
                b = M[i][c]
                M[i] = primitive([a * x - b * y for x, y in zip(M[i], pr)])
        piv.append(c)
        r += 1
    return M[:r], piv


def rank_int(rows, ncols):
    return len(rref_int(rows, ncols)[1])


def nullspace_int(rows, ncols):
    R, piv = rref_int(rows, ncols)
    pivset = set(piv)
    L = 1
    for k, c in enumerate(piv):
        a = abs(R[k][c])
        L = L * a // math.gcd(L, a)
    basis = []
    for f in range(ncols):
        if f in pivset:
            continue
        v = [0] * ncols
        v[f] = L
        for k, c in enumerate(piv):
            num = -R[k][f] * L
            a = R[k][c]
            if num % a != 0:
                raise ArithmeticError("inexact back substitution")
            v[c] = num // a
        basis.append(primitive(v))
    exact = all(sum(x * y for x, y in zip(rw, v)) == 0 for rw in rows for v in basis)
    return basis, len(piv), exact


def imat_mul(A, B):
    k = len(B)
    m = len(B[0])
    return [[sum(A[i][l] * B[l][j] for l in range(k)) for j in range(m)]
            for i in range(len(A))]


def imat_shift(A, s):
    return [[A[i][j] + (s if i == j else 0) for j in range(len(A[0]))]
            for i in range(len(A))]


def imat_vec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


# ------------------------------------------------- monomial matrix groups
# element = (perm, exps): row i holds zeta^exps[i] in column perm[i];
# zeta a primitive M-th root of unity, M even, -1 = zeta^(M//2).
def mono_mul(x, y, M):
    p, k = x
    q, m = y
    d = len(p)
    return (tuple(q[p[i]] for i in range(d)),
            tuple((k[i] + m[p[i]]) % M for i in range(d)))


def mono_inv(x, M):
    p, k = x
    d = len(p)
    q = [0] * d
    m = [0] * d
    for i in range(d):
        q[p[i]] = i
        m[p[i]] = (-k[i]) % M
    return (tuple(q), tuple(m))


def mono_id(d):
    return (tuple(range(d)), tuple([0] * d))


def mono_det_exp(x, M):
    p, k = x
    return (sum(k) + ((M // 2) if perm_parity(p) else 0)) % M


def mono_closure(gens, M):
    d = len(gens[0][0])
    e = mono_id(d)
    seen = {e}
    frontier = [e]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = mono_mul(x, g, M)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return seen


# ring Z[x] mod (x^n + 1): x maps to zeta = exp(i pi k ... ) of order 2n
def cyc_mono(k, n):
    k %= 2 * n
    v = [0] * n
    if k < n:
        v[k] = 1
    else:
        v[k - n] = -1
    return v


def cyc_mul(a, b, n):
    c = [0] * n
    for i in range(n):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(n):
            bj = b[j]
            if bj == 0:
                continue
            t = i + j
            if t >= n:
                c[t - n] -= ai * bj
            else:
                c[t] += ai * bj
    return c


def cyc_conj(a, n):
    c = [0] * n
    c[0] = a[0]
    for i in range(1, n):
        c[n - i] = -a[i]
    return c


def cmat(x, n):
    p, k = x
    d = len(p)
    A = [[[0] * n for _ in range(d)] for _ in range(d)]
    for i in range(d):
        A[i][p[i]] = cyc_mono(k[i], n)
    return A


def cmat_mul(A, B, n):
    d = len(A)
    C = [[None] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            acc = [0] * n
            for l in range(d):
                t = cyc_mul(A[i][l], B[l][j], n)
                acc = [u + w for u, w in zip(acc, t)]
            C[i][j] = acc
    return C


def cmat_adj(A, n):
    d = len(A)
    return [[cyc_conj(A[j][i], n) for j in range(d)] for i in range(d)]


def cmat_det(A, n):
    d = len(A)
    acc = [0] * n
    for s in itertools.permutations(range(d)):
        t = cyc_mono(0, n)
        for i in range(d):
            t = cyc_mul(t, A[i][s[i]], n)
        sg = perm_sign(s)
        acc = [u + sg * w for u, w in zip(acc, t)]
    return acc


def cmat_identity(d, n):
    return [[cyc_mono(0, n) if i == j else [0] * n for j in range(d)]
            for i in range(d)]


def cmat_is_monomial(A):
    d = len(A)
    r_ok = all(sum(1 for j in range(d) if any(A[i][j])) == 1 for i in range(d))
    c_ok = all(sum(1 for i in range(d) if any(A[i][j])) == 1 for j in range(d))
    return r_ok and c_ok


# ----------------------------------------------- finite group utilities
def group_tables(els, mul):
    els = sorted(els)
    idx = {g: i for i, g in enumerate(els)}
    tab = [[idx[mul(x, y)] for y in els] for x in els]
    return els, idx, tab


def element_orders(tab, e):
    out = {}
    for i in range(len(tab)):
        x = i
        k = 1
        while x != e:
            x = tab[x][i]
            k += 1
        out[k] = out.get(k, 0) + 1
    return dict(sorted(out.items()))


def abelian_normal_search(tab, e):
    """Exhaustive: every abelian normal subgroup is a union of pairwise
    commuting conjugacy classes containing the identity; enumerate every such
    union by backtracking and keep those closed under multiplication."""
    N = len(tab)
    inv = [0] * N
    for i in range(N):
        for j in range(N):
            if tab[i][j] == e:
                inv[i] = j
                break
    seen = [False] * N
    classes = []
    for x in range(N):
        if seen[x]:
            continue
        cl = set()
        for g in range(N):
            cl.add(tab[tab[g][x]][inv[g]])
        for y in cl:
            seen[y] = True
        classes.append(sorted(cl))

    def comm(A, B):
        for x in A:
            tx = tab[x]
            for y in B:
                if tx[y] != tab[y][x]:
                    return False
        return True

    cand = [k for k, cl in enumerate(classes) if e not in cl and comm(cl, cl)]
    compat = {}
    for i in cand:
        for j in cand:
            if i < j:
                compat[(i, j)] = comm(classes[i], classes[j])
    best = {"order": 1, "set": [e]}
    stats = {"enumerated": 0, "subgroups": 0, "max_count": 0}

    def closed(S):
        Ss = set(S)
        for x in S:
            tx = tab[x]
            for y in S:
                if tx[y] not in Ss:
                    return False
        return True

    def rec(pos, chosen):
        stats["enumerated"] += 1
        S = [e]
        for k in chosen:
            S.extend(classes[k])
        if closed(S):
            stats["subgroups"] += 1
            if len(S) > best["order"]:
                best["order"] = len(S)
                best["set"] = sorted(S)
                best["all"] = [sorted(S)]
                stats["max_count"] = 1
            elif len(S) == best["order"]:
                best.setdefault("all", []).append(sorted(S))
                stats["max_count"] += 1
        for q in range(pos, len(cand)):
            k = cand[q]
            if all(compat[(min(k, h), max(k, h))] for h in chosen):
                rec(q + 1, chosen + [k])

    rec(0, [])
    return {"classes": len(classes), "best_order": best["order"],
            "best_set": best["set"], "best_all": best.get("all", [best["set"]]),
            "abelian_normal": stats["subgroups"],
            "unions_enumerated": stats["enumerated"],
            "max_count": stats["max_count"]}


# ============================================================ ITEM 1 + 2
def delta_gens(n, bad_b=False):
    M = 2 * n
    a = ((0, 1, 2), (2 % M, (-2) % M, 0))          # diag(w, w^-1, 1), w = zeta^2
    c = ((1, 2, 0), (0, 0, 0))                      # cyclic permutation matrix
    if bad_b:
        b = ((0, 2, 1), (0, 0, 0))                  # transposition with +1 entries
    else:
        b = ((0, 2, 1), (n % M, n % M, n % M))      # -(transposition), det = +1
    return a, b, c, M


def item1_item2():
    groups = {}
    for n in range(1, 8):
        a, b, c, M = delta_gens(n)
        G = mono_closure([a, b, c], M)
        groups[n] = G
        order = len(G)
        expect = 6 * n * n
        all_det1 = all(mono_det_exp(g, M) == 0 for g in G)
        D = sorted(g for g in G if g[0] == (0, 1, 2))
        Dset = set(D)
        grid = set()
        for j1 in range(n):
            for j2 in range(n):
                j3 = (-j1 - j2) % n
                grid.add(((0, 1, 2), ((2 * j1) % M, (2 * j2) % M, (2 * j3) % M)))
        index = order // len(D) if order % len(D) == 0 else -1
        normal = True
        for g in G:
            gi = mono_inv(g, M)
            for d in D:
                if mono_mul(mono_mul(g, d, M), gi, M) not in Dset:
                    normal = False
                    break
            if not normal:
                break
        H = mono_closure([b, c], M)
        Hperms = set(h[0] for h in H)
        H_nonab = any(mono_mul(x, y, M) != mono_mul(y, x, M) for x in H for y in H)
        s3_ok = len(H) == 6 and len(Hperms) == 6 and H_nonab and H <= G
        weyl_image = len(set(g[0] for g in G))
        # explicit matrices over Z[x] mod (x^n+1)
        one = cyc_mono(0, n)
        I3 = cmat_identity(3, n)
        mono_ok = det_ok = unit_ok = True
        mats = {}
        for g in G:
            A = cmat(g, n)
            mats[g] = A
            if not cmat_is_monomial(A):
                mono_ok = False
            if cmat_det(A, n) != one:
                det_ok = False
            if cmat_mul(A, cmat_adj(A, n), n) != I3:
                unit_ok = False
        mul_ok = True
        pairs = 0
        Gl = sorted(G)
        others = Gl if n <= 3 else [a, b, c]
        for x in Gl:
            for y in others:
                pairs += 1
                if cmat_mul(mats[x], cmat(y, n), n) != cmat(mono_mul(x, y, M), n):
                    mul_ok = False
        is_s3 = (order == 6 and s3_ok and G == H) if n == 1 else "na"
        row("DELTA", n=n, zeta_order=M, order=order, expect_6n2=expect,
            n1_is_S3=is_s3, diag_order=len(D), diag_index=index,
            diag_equals_grid_sum0_mod_n=int(Dset == grid), diag_normal=int(normal),
            S3_subgroup_bc_order=len(H), S3_embeds=int(s3_ok), weyl_image=weyl_image,
            all_det1_symbolic=int(all_det1), explicit_monomial=int(mono_ok),
            explicit_det1=int(det_ok), explicit_unitary=int(unit_ok),
            mult_crosscheck_pairs=pairs, mult_crosscheck=int(mul_ok),
            in_NT="1 (monomial unitary det1; g diag(t) g^-1 is a permuted diag)",
            tag="MEASURED")
        check("delta_order_n%d" % n, order == expect, order=order, expect=expect)
        check("delta_diag_index6_n%d" % n, index == 6 and len(D) == n * n,
              diag_order=len(D), index=index)
        check("delta_diag_grid_n%d" % n, Dset == grid and normal)
        check("delta_S3_embeds_n%d" % n, s3_ok and weyl_image == 6)
        check("delta_SU3_explicit_n%d" % n,
              all_det1 and mono_ok and det_ok and unit_ok and mul_ok,
              pairs=pairs)
        if n == 1:
            check("delta_n1_is_S3", is_s3 is True, order=order)

    # tower inclusions Delta(6n^2) in Delta(6m^2) for n | m
    incl = []
    ok_all = True
    for n in range(1, 8):
        for m in range(n + 1, 8):
            if m % n:
                continue
            f = m // n
            img = set((g[0], tuple((x * f) % (2 * m) for x in g[1])) for g in groups[n])
            ok = img <= groups[m]
            ok_all = ok_all and ok
            incl.append("%d<%d:%d" % (n, m, int(ok)))
    row("DELTA_TOWER", inclusions=" ".join(incl), all_ok=int(ok_all), tag="MEASURED")
    check("delta_tower_inclusions", ok_all, count=len(incl))

    # control: b with +1 entries has det -1
    a, b, c, M = delta_gens(3, bad_b=True)
    Gbad = mono_closure([a, b, c], M)
    bad_det1 = all(mono_det_exp(g, M) == 0 for g in Gbad)
    control("delta_b_plus_entries_det", bad_det1, n=3, order=len(Gbad),
            det_b_exp=mono_det_exp(b, M))

    # ITEM 2: Jordan bound, exhaustive class search, every n = 1..7
    for n in range(1, 8):
        M = 2 * n
        els, idx, tab = group_tables(groups[n], lambda x, y: mono_mul(x, y, M))
        e = idx[mono_id(3)]
        res = abelian_normal_search(tab, e)
        order = len(els)
        best = res["best_order"]
        index = order // best if order % best == 0 else -1
        Didx = sorted(idx[g] for g in groups[n] if g[0] == (0, 1, 2))
        diag_among_max = int(Didx in res["best_all"])
        Dset_i = set(Didx)
        others = []
        for S in res["best_all"]:
            if S == Didx:
                continue
            cyc3 = sum(1 for i in S if els[i][0] in ((1, 2, 0), (2, 0, 1)))
            others.append("capD%d_3cycles%d" % (len(Dset_i & set(S)), cyc3))
        row("JORDAN", group="Delta(6n^2)", n=n, order=order, classes=res["classes"],
            abelian_normal_subgroups=res["abelian_normal"],
            unions_enumerated=res["unions_enumerated"],
            max_abelian_normal_order=best, max_count=res["max_count"],
            min_index=index, bound=6, diag_among_max=diag_among_max,
            other_max=(" ".join(others) if others else "none"), exhaustive=1,
            tag="MEASURED")
        check("jordan_index_le6_n%d" % n, 0 < index <= 6, index=index)
        if n >= 2:
            check("jordan_max_order_n2_diag_among_max_n%d" % n,
                  best == n * n and diag_among_max == 1, max_count=res["max_count"])
    row("CORRECTION", tag="MEASURED",
        text="First run asserted the largest abelian normal subgroup is unique and equals the "
             "diagonal D; measured max_count=4 at n=3 (D plus three order-9 abelian normal "
             "subgroups that meet D in the centre Z3 and hold 3-cycle elements), so the "
             "uniqueness check failed there. The check now asserts max order n^2 (index 6) "
             "with D among the maximal ones; the uniqueness expectation is retracted for n=3.")
    return groups


# ================================================================ ITEM 3
def item3():
    steps = [
        ("1", 0, "Jordan 1878 (J. reine angew. Math. 84): every finite subgroup F of GL(N,C) "
                 "has an abelian normal subgroup A with [F:A] <= J(N); optimal J(2)=60 and "
                 "J(3)=360 (Collins 2007, J. Group Theory 10). Anchors measured in this file: "
                 "2I attains index 60; Delta(6n^2) has index 6."),
        ("2", 0, "Let F_k be finite subgroups of SU(N) with F_k -> H in the Hausdorff metric on "
                 "closed subsets (for compact SU(N) this is the Chabauty topology); H is a closed "
                 "subgroup. Pick A_k abelian normal in F_k, so F_k is a union of m_k <= J(N) "
                 "cosets g_(k,i) A_k; pass to a subsequence with m_k = m constant."),
        ("3", 1, "COMPACTNESS: SU(N) is compact, so its closed subsets form a compact space in "
                 "the Hausdorff metric (Blaschke selection theorem) and SU(N) is sequentially "
                 "compact; pass to a subsequence with A_k -> A and g_(k,i) -> g_i for all i < m."),
        ("4", 0, "A is a closed abelian subgroup: for x, y in A choose x_k, y_k in A_k with "
                 "x_k -> x, y_k -> y; then x_k y_k^-1 lies in A_k and x_k y_k = y_k x_k, and "
                 "continuity of multiplication and inversion gives x y^-1 in A and xy = yx."),
        ("5", 0, "Hausdorff limits commute with finite unions and with translation "
                 "(d(gX, hY) <= d(g,h) + d(X,Y) for a bi-invariant metric), so H = union over "
                 "i < m of g_i A with every g_i in H and A a subgroup of H."),
        ("6", 0, "A has index <= m in H; its complement is a finite union of closed cosets, so A "
                 "is open and closed in H and contains the identity component H0; so H0 is "
                 "abelian. H is a compact Lie group (Cartan closed subgroup theorem), so H0 is a "
                 "torus of dimension <= rank SU(N) = N-1."),
        ("7", 0, "dim H <= N-1 < N^2-1 = dim SU(N) for N >= 2: every limit of finite subgroups "
                 "is torus-by-finite (H0 a torus, H mod H0 finite) of dimension <= N-1, and "
                 "SU(N) lies outside this class."),
        ("8", 1, "COROLLARY, COMPACTNESS again: there is eps_N > 0 with Hausdorff distance "
                 "d(F, SU(N)) >= eps_N for every finite subgroup F of SU(N); a sequence with "
                 "d -> 0 would converge to SU(N), against step 7."),
        ("9", 0, "Attainment of dimension N-1: Delta(6n^2) -> N(T) = T semidirect S3 in SU(3) (dim 2) and "
                 "Dic_m -> N(T) = T u jT in SU(2) (dim 1). Measured ingredients: Delta cap T is "
                 "the full grid diag(w^j1, w^j2, w^j3) with j1+j2+j3 = 0 mod n (n^2 points, angle "
                 "spacing 2 pi over n), the Weyl image is all of S3; Dic_m cap T is all 2m-th "
                 "roots. So the Hausdorff distance is at most the grid spacing and tends to 0."),
    ]
    for sid, comp, text in steps:
        row("LIMIT", step=sid, uses_compactness=comp, tag="NAMED", text=text)
    row("CONSEQUENCE", tag="NAMED",
        text="Finite gauge groups approach at most a torus-by-finite group of dimension N-1, "
             "such as N(T) = T semidirect W; the identity component of any limit is abelian. A "
             "finite-group lattice theory with fixed G in SU(N) tracks the SU(N) theory only "
             "below a freezing coupling (Petcher-Weingarten 1980 Phys Rev D 22 2465 for "
             "finite SU(2) subgroups incl 2I; Bhanot-Rebbi 1981 Phys Rev D 24 3319 for finite "
             "SU(3) subgroups).")
    row("SCOPE", tag="CONJECTURE",
        text="For the finite-universe Yang-Mills program: finitely many sites with the compact "
             "group SU(3) on every link keeps the nonabelian identity component; a tower of "
             "finite subgroups such as Delta(6n^2) carries exactly the maximal-torus plus "
             "Weyl-group sector N(T), the abelian-projection sector, in the limit.")


# ================================================================ ITEM 4
def rmul(x, y, P, Q):
    bd = x[1] * y[1]
    return (x[0] * y[0] + bd * Q, x[0] * y[1] + x[1] * y[0] + bd * P)


def radd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def rsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def rsc(x, s):
    return (x[0] * s, x[1] * s)


def rhalf(x):
    if x[0] % 2 or x[1] % 2:
        raise ArithmeticError("not divisible by 2")
    return (x[0] // 2, x[1] // 2)


def q2mul(x, y, P, Q):
    """Product of quaternions stored scaled by 2: returns 2*(q1 q2)."""
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y

    def m(u, v):
        return rmul(u, v, P, Q)

    r0 = rsub(rsub(rsub(m(a1, a2), m(b1, b2)), m(c1, c2)), m(d1, d2))
    r1 = rsub(radd(radd(m(a1, b2), m(b1, a2)), m(c1, d2)), m(d1, c2))
    r2 = radd(radd(rsub(m(a1, c2), m(b1, d2)), m(c1, a2)), m(d1, b2))
    r3 = radd(rsub(radd(m(a1, d2), m(b1, c2)), m(c1, b2)), m(d1, a2))
    return (rhalf(r0), rhalf(r1), rhalf(r2), rhalf(r3))


def q2norm(x, P, Q):
    acc = (0, 0)
    for c in x:
        acc = radd(acc, rmul(c, c, P, Q))
    return acc


def q2closed(S, P, Q):
    try:
        for x in S:
            for y in S:
                if q2mul(x, y, P, Q) not in S:
                    return False
    except ArithmeticError:
        return False
    return True


def q2closure(gens, P, Q):
    e = ((2, 0), (0, 0), (0, 0), (0, 0))
    seen = {e}
    frontier = [e]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = q2mul(x, g, P, Q)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
        if len(seen) > 100000:
            break
    return seen


ZR = (0, 0)
ONE = (1, 0)
PHI_R = (0, 1)          # phi, with phi^2 = phi + 1   (P=1, Q=1)
PHIINV_R = (-1, 1)      # phi^-1 = phi - 1


def hurwitz_units():
    S = set()
    for pos in range(4):
        for s in (1, -1):
            v = [ZR] * 4
            v[pos] = (2 * s, 0)
            S.add(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        S.add(tuple((s, 0) for s in signs))
    return S


def icosians():
    S = hurwitz_units()
    base = [ZR, ONE, PHIINV_R, PHI_R]
    for perm in itertools.permutations(range(4)):
        if perm_parity(perm):
            continue
        for s1, s2, s3 in itertools.product((1, -1), repeat=3):
            vals = [ZR, rsc(base[1], s1), rsc(base[2], s2), rsc(base[3], s3)]
            v = [None] * 4
            for i in range(4):
                v[perm[i]] = vals[i]
            S.add(tuple(v))
    return S


def octahedral_2O():
    S = hurwitz_units()                 # over Z[sqrt2] (P=0, Q=2) same coordinates
    for i, j in itertools.combinations(range(4), 2):
        for s, t in itertools.product((1, -1), repeat=2):
            v = [ZR] * 4
            v[i] = (0, s)               # 2 * (+-1 over sqrt2) = +-sqrt2
            v[j] = (0, t)
            S.add(tuple(v))
    return S


def item4():
    P, Q = 1, 1
    I120 = icosians()
    norms = all(q2norm(x, P, Q) == (4, 0) for x in I120)
    closed = q2closed(I120, P, Q)
    g1 = ((1, 0), (1, 0), (1, 0), (1, 0))              # (1+i+j+k) over 2
    g2 = (PHI_R, PHIINV_R, ONE, ZR)                     # (phi + phi^-1 i + j) over 2
    gen = q2closure([g1, g2], P, Q)
    els, idx, tab = group_tables(I120, lambda x, y: q2mul(x, y, P, Q))
    e = idx[((2, 0), ZR, ZR, ZR)]
    orders = element_orders(tab, e)
    exp_orders = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}
    res = abelian_normal_search(tab, e)
    center = [els[i] for i in range(len(els)) if all(tab[i][j] == tab[j][i] for j in range(len(els)))]
    row("SU2_2I", ring="Z[phi] phi^2=phi+1", coords="quaternion times 2 as (a,b)=a+b phi",
        order=len(I120), all_norm1=int(norms), closed=int(closed),
        generated_by_g1_g2=int(gen == I120), g1="(1,1,1,1)", g2="(phi,phi-1,1,0)",
        element_orders=orders, center_order=len(center), classes=res["classes"],
        max_abelian_normal_order=res["best_order"],
        jordan_index=len(I120) // res["best_order"], tag="MEASURED")
    check("su2_2I_order120_closed", len(I120) == 120 and norms and closed)
    check("su2_2I_generated_by_two", gen == I120, generated=len(gen))
    check("su2_2I_element_orders", orders == exp_orders, measured=orders)
    check("su2_2I_jordan_index60", res["best_order"] == 2 and len(center) == 2,
          index=len(I120) // res["best_order"])

    T24 = hurwitz_units()
    t_closed = q2closed(T24, P, Q)
    elsT, idxT, tabT = group_tables(T24, lambda x, y: q2mul(x, y, P, Q))
    resT = abelian_normal_search(tabT, idxT[((2, 0), ZR, ZR, ZR)])
    ordT = element_orders(tabT, idxT[((2, 0), ZR, ZR, ZR)])
    row("SU2_2T", order=len(T24), closed=int(t_closed), subset_of_2I=int(T24 <= I120),
        element_orders=ordT, max_abelian_normal_order=resT["best_order"],
        jordan_index=len(T24) // resT["best_order"], tag="MEASURED")
    check("su2_2T_order24", len(T24) == 24 and t_closed and T24 <= I120
          and ordT == {1: 1, 2: 1, 3: 8, 4: 6, 6: 8})

    P2, Q2 = 0, 2
    O48 = octahedral_2O()
    o_norm = all(q2norm(x, P2, Q2) == (4, 0) for x in O48)
    o_closed = q2closed(O48, P2, Q2)
    elsO, idxO, tabO = group_tables(O48, lambda x, y: q2mul(x, y, P2, Q2))
    eO = idxO[((2, 0), ZR, ZR, ZR)]
    ordO = element_orders(tabO, eO)
    resO = abelian_normal_search(tabO, eO)
    row("SU2_2O", ring="Z[sqrt2]", order=len(O48), all_norm1=int(o_norm), closed=int(o_closed),
        element_orders=ordO, max_abelian_normal_order=resO["best_order"],
        jordan_index=len(O48) // resO["best_order"], tag="MEASURED")
    check("su2_2O_order48", len(O48) == 48 and o_norm and o_closed
          and ordO == {1: 1, 2: 1, 3: 8, 4: 18, 6: 8, 8: 12})
    row("ADE", tag="NAMED",
        text="Finite subgroups of SU(2) up to conjugacy: cyclic Z_n (A), binary dihedral "
             "Dic_m of order 4m (D), 2T order 24 (E6), 2O order 48 (E7), 2I order 120 (E8) "
             "(Klein 1884 Vorlesungen ueber das Ikosaeder; McKay 1980). Measured orders "
             "24, 48, 120: 2I is the largest exceptional one.")
    check("su2_ade_exceptional_max_is_2I", max(len(T24), len(O48), len(I120)) == len(I120) == 120)

    # binary dihedral towers Dic_m inside N(T) of SU(2)
    dic = {}
    ok_all = True
    for m in range(2, 9):
        M = 2 * m
        a = ((0, 1), (1, M - 1))         # diag(zeta, zeta^-1), zeta of order 2m
        b = ((1, 0), (0, m))             # [[0, 1], [-1, 0]], b^2 = -1
        G = mono_closure([a, b], M)
        dic[m] = G
        det1 = all(mono_det_exp(g, M) == 0 for g in G)
        Tpart = [g for g in G if g[0] == (0, 1)]
        Jpart = [g for g in G if g[0] == (1, 0)]
        t_res = sorted(set(g[1][0] for g in Tpart))
        j_res = sorted(set(g[1][0] for g in Jpart))
        am = a
        for _ in range(m - 1):
            am = mono_mul(am, a, M)
        a2m = mono_mul(am, am, M)
        rel1 = a2m == mono_id(2)
        rel2 = mono_mul(b, b, M) == am
        rel3 = mono_mul(mono_mul(b, a, M), mono_inv(b, M), M) == mono_inv(a, M)
        mono_ok = all(cmat_is_monomial(cmat(g, m)) for g in G)
        unit_ok = all(cmat_mul(cmat(g, m), cmat_adj(cmat(g, m), m), m) == cmat_identity(2, m) for g in G)
        els, idx, tab = group_tables(G, lambda x, y: mono_mul(x, y, M))
        res = abelian_normal_search(tab, idx[mono_id(2)])
        ok = (len(G) == 4 * m and det1 and len(Tpart) == 2 * m and len(Jpart) == 2 * m
              and t_res == list(range(M)) and j_res == list(range(M))
              and rel1 and rel2 and rel3 and mono_ok and unit_ok)
        ok_all = ok_all and ok
        row("SU2_DIC", m=m, order=len(G), expect_4m=4 * m, det1=int(det1),
            monomial=int(mono_ok), unitary=int(unit_ok), in_NT=int(mono_ok),
            T_part=len(Tpart), T_residues_all=int(t_res == list(range(M))),
            jT_part=len(Jpart), jT_residues_all=int(j_res == list(range(M))),
            relations=int(rel1 and rel2 and rel3),
            max_abelian_normal_order=res["best_order"],
            jordan_index=len(G) // res["best_order"], tag="MEASURED")
        check("su2_dic_m%d" % m, ok, order=len(G))
    incl = []
    inc_ok = True
    for m in range(2, 9):
        for k in range(m + 1, 9):
            if k % m:
                continue
            f = k // m
            img = set((g[0], tuple((x * f) % (2 * k) for x in g[1])) for g in dic[m])
            ok = img <= dic[k]
            inc_ok = inc_ok and ok
            incl.append("%d<%d:%d" % (m, k, int(ok)))
    row("SU2_DIC_TOWER", inclusions=" ".join(incl), all_ok=int(inc_ok),
        dense_in_NT="union over m = 2^k is dense in N(T) = T u jT: each coset holds all "
                    "2m-th roots (measured), spacing 2 pi over 2m -> 0",
        tag="MEASURED+NAMED")
    check("su2_dic_tower", inc_ok)

    # controls
    bad = (ONE, PHIINV_R, ONE, ZR)                      # phi replaced by 1
    control("icosian_norm_bad_generator", q2norm(bad, P, Q) == (4, 0),
            norm=q2norm(bad, P, Q))
    holed = set(I120)
    holed.discard(g2)
    control("icosian_closure_with_hole", q2closed(holed, P, Q), order=len(holed))


# ================================================================ ITEM 5
TRIPLES = ((1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 7, 5), (3, 7, 4), (3, 6, 5))
ALL7 = (1, 2, 3, 4, 5, 6, 7)


def build_omul(triples):
    tab = {}
    dup = False
    for a in range(8):
        tab[(0, a)] = (1, a)
        tab[(a, 0)] = (1, a)
    for a in range(1, 8):
        tab[(a, a)] = (-1, 0)
    used = set()
    for (x, y, z) in triples:
        for (u, v, w) in ((x, y, z), (y, z, x), (z, x, y)):
            if (u, v) in used or (v, u) in used:
                dup = True
            used.add((u, v))
            tab[(u, v)] = (1, w)
            tab[(v, u)] = (-1, w)
    complete = all((a, b) in tab for a in range(8) for b in range(8))
    return tab, complete and not dup


def omul_vec(OM, x, y):
    out = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            s, k = OM[(a, b)]
            out[k] += s * x[a] * y[b]
    return out


def unit(a):
    v = [0] * 8
    v[a] = 1
    return v


def wedge_basis(I, J):
    if set(I) & set(J):
        return None
    seq = tuple(I) + tuple(J)
    return perm_sign(seq), tuple(sorted(seq))


def wedge(F, G):
    out = {}
    for I, a in F.items():
        for J, b in G.items():
            r = wedge_basis(I, J)
            if r is None:
                continue
            s, K = r
            out[K] = out.get(K, 0) + s * a * b
    return {K: v for K, v in out.items() if v}


def star7(F):
    out = {}
    for I, a in F.items():
        J = tuple(x for x in ALL7 if x not in I)
        s = perm_sign(tuple(I) + J)
        out[J] = out.get(J, 0) + s * a
    return {K: v for K, v in out.items() if v}


def interior(v, F):
    out = {}
    for I, a in F.items():
        if v in I:
            p = I.index(v)
            K = tuple(x for x in I if x != v)
            out[K] = out.get(K, 0) + (-a if p & 1 else a)
    return {K: c for K, c in out.items() if c}


L2 = list(itertools.combinations(ALL7, 2))


def vec2(F):
    return [F.get(b, 0) for b in L2]


def form2(v):
    return {L2[i]: v[i] for i in range(len(L2)) if v[i]}


def Tmatrix(phi):
    cols = [vec2(star7(wedge(phi, {b: 1}))) for b in L2]
    return transpose(cols)


def eigen_split(phi):
    T = Tmatrix(phi)
    n = len(T)
    sym = all(T[i][j] == T[j][i] for i in range(n) for j in range(n))
    tr = sum(T[i][i] for i in range(n))
    Pm = imat_mul(imat_shift(T, -2), imat_shift(T, 1))
    minpoly = all(x == 0 for r in Pm for x in r)
    r2 = rank_int(imat_shift(T, -2), n)
    r1 = rank_int(imat_shift(T, 1), n)
    d7 = n - r2
    d14 = n - r1
    ok = sym and minpoly and d7 == 7 and d14 == 14
    return ok, {"symmetric": int(sym), "trace": tr, "minpoly_T2_T1_zero": int(minpoly),
                "dim_eig_plus2": d7, "dim_eig_minus1": d14}, T


def phi_metric(phi):
    G = [[0] * 7 for _ in range(7)]
    for i in ALL7:
        bi = interior(i, phi)
        for j in ALL7:
            f = wedge(wedge(bi, interior(j, phi)), phi)
            G[i - 1][j - 1] = f.get(ALL7, 0)
    return G


def form_str(F):
    parts = []
    for K in sorted(F):
        c = F[K]
        sgn = "+" if c > 0 else "-"
        mag = "" if abs(c) == 1 else str(abs(c))
        parts.append(sgn + mag + "e" + "".join(str(x) for x in K))
    return " ".join(parts)


def item5_item6():
    OM, table_ok = build_omul(TRIPLES)
    E = [unit(a) for a in range(8)]
    # alternative: associator antisymmetric on basis triples
    def assoc(a, b, c):
        l = omul_vec(OM, omul_vec(OM, E[a], E[b]), E[c])
        r = omul_vec(OM, E[a], omul_vec(OM, E[b], E[c]))
        return [x - y for x, y in zip(l, r)]
    alt_ok = True
    for a in range(8):
        for b in range(8):
            for c in range(8):
                A0 = assoc(a, b, c)
                if [-x for x in assoc(b, a, c)] != A0 or [-x for x in assoc(a, c, b)] != A0:
                    alt_ok = False
    nonassoc = any(any(assoc(a, b, c)) for a in range(8) for b in range(8) for c in range(8))
    # composition: polarised norm identity on basis quadruples
    comp_ok = True
    for a in range(8):
        for b in range(8):
            ab = omul_vec(OM, E[a], E[b])
            for c in range(8):
                for d in range(8):
                    cd = omul_vec(OM, E[c], E[d])
                    ad = omul_vec(OM, E[a], E[d])
                    cb = omul_vec(OM, E[c], E[b])
                    lhs = sum(x * y for x, y in zip(ab, cd)) + sum(x * y for x, y in zip(ad, cb))
                    rhs = 2 * (1 if a == c else 0) * (1 if b == d else 0)
                    if lhs != rhs:
                        comp_ok = False
    row("OCT_TABLE", name="Harvey-Lawson-Bryant",
        rule="e_a e_b = e_c for (a,b,c) in 123 145 167 246 275 374 365 and cyclic; "
             "e_b e_a = -e_c; e_a e_a = -1; e_0 = 1",
        complete=int(table_ok), alternative=int(alt_ok), nonassociative=int(nonassoc),
        composition_polarised=int(comp_ok), tag="MEASURED")
    check("oct_table_is_octonions", table_ok and alt_ok and comp_ok and nonassoc)

    # phi(a,b,c) = <e_a, e_b e_c> on Im O
    phit = [[[0] * 8 for _ in range(8)] for _ in range(8)]
    for a in ALL7:
        for b in ALL7:
            for c in ALL7:
                phit[a][b][c] = omul_vec(OM, E[b], E[c])[a]
    alt3 = all(phit[a][b][c] == -phit[b][a][c] == -phit[a][c][b]
               for a in ALL7 for b in ALL7 for c in ALL7)
    PHI = {}
    for t in itertools.combinations(ALL7, 3):
        v = phit[t[0]][t[1]][t[2]]
        if v:
            PHI[t] = v
    PSI = star7(PHI)
    row("PHI", definition="phi(a,b,c)=<a,bc>", form=form_str(PHI), terms=len(PHI),
        alternating=int(alt3), tag="MEASURED")
    row("PSI", definition="psi=*phi, orientation e1234567", form=form_str(PSI),
        terms=len(PSI), tag="MEASURED")
    check("phi_alternating_7_terms", alt3 and len(PHI) == 7)
    g = phi_metric(PHI)
    g_ok = all(g[i][j] == (6 if i == j else 0) for i in range(7) for j in range(7))
    row("PHI_METRIC", formula="(x contract phi)^(y contract phi)^phi = g(x,y) e1234567",
        diag=" ".join(str(g[i][i]) for i in range(7)),
        offdiag_zero=int(all(g[i][j] == 0 for i in range(7) for j in range(7) if i != j)),
        equals_6_delta=int(g_ok), orientation_positive=int(g_ok), tag="MEASURED")
    check("phi_metric_6delta_positive_orientation", g_ok)

    ok, info, T = eigen_split(PHI)
    row("T_EIGEN", operator="T(beta)=*(phi^beta) on Lambda2 R7", dim=21,
        eigenvalues="2 and -1", **info, tag="MEASURED")
    check("T_eigen_split_7_14", ok, **info)

    # Lambda2_7 = { v _| phi }
    L7 = [vec2(interior(v, PHI)) for v in ALL7]
    l7_eig = all(imat_vec(T, w) == [2 * x for x in w] for w in L7)
    l7_rank = rank_int(L7, 21)
    ks = []
    for v in ALL7:
        f6 = wedge(interior(v, PHI), PSI)
        sv = star7({(v,): 1})
        kk = None
        if f6 and set(f6) == set(sv):
            r = [f6[K] * sv[K] for K in sv]   # sv entries are +-1
            if len(set(r)) == 1:
                kk = r[0]
        ks.append(kk)
    row("LAMBDA2_7", basis="v contract phi, v=1..7", eig_2=int(l7_eig), rank=l7_rank,
        beta_wedge_psi_over_star_v=" ".join(str(k) for k in ks),
        note="(v contract phi)^psi = k *e^v with k measured; Bryant 2005 gives 3", tag="MEASURED")
    check("lambda2_7_is_v_contract_phi", l7_eig and l7_rank == 7 and all(k == 3 for k in ks))

    # der(O): 8x8 unknowns D[a][b] = coefficient of e_a in D(e_b)
    def u(a, b):
        return a * 8 + b
    eqs = []
    for i in range(8):
        for j in range(8):
            rows_c = [[0] * 64 for _ in range(8)]
            s, k = OM[(i, j)]
            for c in range(8):
                rows_c[c][u(c, k)] += s
            for a in range(8):
                s2, mm = OM[(a, j)]
                rows_c[mm][u(a, i)] -= s2
                s3, mm = OM[(i, a)]
                rows_c[mm][u(a, j)] -= s3
            eqs.extend(rows_c)
    der, rk, exact = nullspace_int(eqs, 64)
    kills1 = all(all(v[u(a, 0)] == 0 and v[u(0, a)] == 0 for a in range(8)) for v in der)
    antisym = all(all(v[u(a, b)] == -v[u(b, a)] for a in ALL7 for b in ALL7) for v in der)
    der7 = [[v[u(i, j)] for i in ALL7 for j in ALL7] for v in der]
    der_forms = [[v[u(i, j)] for (i, j) in L2] for v in der]
    der_in_14 = all(imat_vec(T, w) == [-x for x in w] for w in der_forms)
    der_rank2 = rank_int(der_forms, 21)
    psi_zero = all(not wedge(form2(w), PSI) for w in der_forms)
    psi_nonzero_on_7 = all(bool(wedge(form2(w), PSI)) for w in L7)
    row("DER_O", unknowns=64, equations=len(eqs), dim=len(der), exact=int(exact),
        kills_1_preserves_ImO=int(kills1), antisymmetric=int(antisym),
        as_2forms_in_eig_minus1=int(der_in_14), rank_as_2forms=der_rank2,
        beta_wedge_psi_zero=int(psi_zero), psi_nonzero_on_Lambda2_7=int(psi_nonzero_on_7),
        tag="MEASURED")
    check("der_O_equals_Lambda2_14", len(der) == 14 and exact and kills1 and antisym
          and der_in_14 and der_rank2 == 14 and psi_zero and psi_nonzero_on_7)

    # stabiliser of phi in gl(7)
    seqs = []
    for (a, b, c) in itertools.combinations(ALL7, 3):
        r = [0] * 49
        for d in ALL7:
            r[(d - 1) * 7 + (a - 1)] += phit[d][b][c]
            r[(d - 1) * 7 + (b - 1)] += phit[a][d][c]
            r[(d - 1) * 7 + (c - 1)] += phit[a][b][d]
        seqs.append(r)
    stab, rks, exs = nullspace_int(seqs, 49)
    st_anti = all(all(v[(i - 1) * 7 + (j - 1)] == -v[(j - 1) * 7 + (i - 1)]
                      for i in ALL7 for j in ALL7) for v in stab)
    joint = rank_int(stab + der7, 49)
    row("STAB_PHI", unknowns=49, equations=len(seqs), dim=len(stab), exact=int(exs),
        antisymmetric=int(st_anti), rank_stab_plus_der=joint, equals_der_O=int(joint == 14),
        tag="MEASURED")
    check("stab_phi_equals_der_O_dim14", len(stab) == 14 and exs and st_anti and joint == 14)

    # associative 3-planes and coassociative 4-planes
    chir = []
    for A in sorted(PHI):
        C = tuple(x for x in ALL7 if x not in A)
        phiA = PHI[A]
        assoc_zero = all(not any(assoc(a, b, c)) for a in A for b in A for c in A)
        phiC_zero = all(phit[a][b][c] == 0 for a in C for b in C for c in C)
        sg, _ = wedge_basis(A, C)
        s1 = phiA * sg                              # vol_A ^ vol_C = vol_7
        s2 = PSI.get(C, 0)                          # psi restricted to C
        # star on Lambda2 C with orientation s2 e^C
        def starC(I):
            J = tuple(x for x in C if x not in I)
            return J, s2 * perm_sign(tuple(I) + J)
        sd, asd = [], []
        star_sq = True
        for J0 in [(C[0], C[1]), (C[0], C[2]), (C[0], C[3])]:
            K, s = starC(J0)
            K2, s_back = starC(K)
            star_sq = star_sq and K2 == J0 and s * s_back == 1
            sd.append(vec2({J0: 1, K: s}))
            asd.append(vec2({J0: 1, K: -s}))

        def dims(B):
            cols = transpose(B)
            in14 = 3 - rank_int(transpose(imat_mul(imat_shift(T, 1), cols)), 21)
            in7 = 3 - rank_int(transpose(imat_mul(imat_shift(T, -2), cols)), 21)
            return in14, in7
        sd14, sd7 = dims(sd)
        asd14, asd7 = dims(asd)
        half = "ASD" if (asd14 == 3 and sd14 == 0) else ("SD" if (sd14 == 3 and asd14 == 0) else "MIXED")
        hv = asd if half == "ASD" else sd
        g2inst = all(not wedge(form2(w), PSI) for w in hv) and half != "MIXED"
        chir.append(half)
        row("COASSOC", A="e" + "".join(map(str, A)), phi_on_A=phiA,
            associator_zero_on_A=int(assoc_zero), C="e" + "".join(map(str, C)),
            phi_on_C_zero=int(phiC_zero), orient_from_A=s1, psi_on_C=s2,
            orientations_agree=int(s1 == s2), starC_sq_1=int(star_sq),
            SD_in_L14=sd14, SD_in_L7=sd7, ASD_in_L14=asd14, ASD_in_L7=asd7,
            half_in_Lambda2_14=half, half_wedge_psi_zero=int(g2inst), tag="MEASURED")
        check("coassoc_%s" % "".join(map(str, C)),
              abs(phiA) == 1 and assoc_zero and phiC_zero and s1 == s2 and abs(s2) == 1
              and star_sq and half != "MIXED" and g2inst and sd7 == 0 and asd7 == 0,
              half=half)
    check("coassoc_same_chirality_all7",
          len(set(chir)) == 1 and len(chir) == 7 and chir[0] in ("SD", "ASD"),
          chirality=chir[0] if chir else "none")
    row("G2_INSTANTON", tag="MEASURED+NAMED", chirality=chir[0] if chir else "none",
        orientation="C oriented by psi on C = vol_A complement orientation (measured equal)",
        text="On R7 = A x C (A associative, C coassociative) a connection pulled back from C, "
             "invariant under translations of A, has curvature F in Lambda2 C tensor g at "
             "every point. If F is anti-self-dual on C (measured: Lambda2_- C lies in "
             "Lambda2_14 and has F^psi = 0), then pi_7(F) = 0, the G2-instanton equation "
             "(Donaldson-Thomas 1998, Gauge theory in higher dimensions). Example: the BPST "
             "instanton (Belavin-Polyakov-Schwarz-Tyupkin 1975) pulled back to R3 x R4 is a "
             "G2-instanton with energy density constant along R3.")

    # ITEM 6: broken controls on phi
    bad = dict(PHI)
    flip = sorted(bad)[0]
    bad[flip] = -bad[flip]
    okb, infob, Tb = eigen_split(bad)
    gb = phi_metric(bad)
    offz = all(gb[i][j] == 0 for i in range(7) for j in range(7) if i != j)
    pos = sum(1 for i in range(7) if gb[i][i] > 0)
    neg = sum(1 for i in range(7) if gb[i][i] < 0)
    control("phi_one_sign_flipped_eigen_split", okb,
            flipped="e" + "".join(map(str, flip)), **infob)
    row("CONTROL_PHI_METRIC", flipped="e" + "".join(map(str, flip)),
        diag=" ".join(str(gb[i][i]) for i in range(7)), offdiag_zero=int(offz),
        positive=pos, negative=neg,
        reading="nondegenerate indefinite B_phi: split-G2 type form (two open GL7 orbits, "
                "Bryant 1987); signature measured from the diagonal",
        tag="MEASURED+NAMED")
    der_bad_ok = False
    try:
        # der(O) check against the flipped T: der forms must leave eig -1
        der_bad_ok = all(imat_vec(Tb, w) == [-x for x in w] for w in der_forms)
    except Exception:
        der_bad_ok = False
    control("der_O_in_eig_minus1_of_flipped_T", der_bad_ok)


# ================================================================== MAIN
def main():
    src, lint = self_lint()
    sha_src = hashlib.sha256(src.encode("utf-8")).hexdigest()
    row("HDR", file=OUT_NAME, script="m3_towers.py", script_sha256=sha_src,
        seat="Opus 5.5 pid 8467a937cba309f7", E=0, network=0, float_used=0,
        integer_only=1, date="2026-11-09")
    check("self_lint_integer_only", not lint, findings=(" ".join(lint) if lint else "none"))
    item1_item2()
    item3()
    item4()
    item5_item6()
    verdict = "ALL_PASS" if not FAILS else "SOME_FAIL"
    row("SUMMARY", checks=COUNT["check"], controls=COUNT["control"], passed=COUNT["pass"],
        failed=len(FAILS), failed_names=(" ".join(FAILS) if FAILS else "none"),
        verdict=verdict)
    body = "\n".join(ROWS) + "\n"
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(body)
    h = hashlib.sha256(body.encode("utf-8")).hexdigest()
    with open(SIDE, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(h + "  " + OUT_NAME + "\n")
    print("HBP_SHA256|sha256=" + h + "|json=0")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
