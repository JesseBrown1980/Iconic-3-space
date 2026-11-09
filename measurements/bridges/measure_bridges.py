#!/usr/bin/env python3
"""
measure_bridges.py: three bridges stated in the Iconic 3-space README, MEASURED in exact integers.

  1. Octonions to E8. A Fano-plane octonion table is stated, every Steiner system S(3,4,8) on the
     eight coordinates is built and tested as an order, the naive Kirmse set is tested, and the
     corrected (Coxeter) order is measured: 240 units, closure, the E8 inner-product distribution
     from every unit, a Z-basis, its Gram matrix, determinant and parity.
     NAMED: Coxeter 1946, Integral Cayley numbers; Conway and Smith ch. 9 (not opened here).
  2. Icosians to E8. The 120 unit icosians are built, closure is tested, the Z[phi]-span is reduced
     to a Z-basis, and the form x -> Tr(lambda N(x)) with lambda = (2 + phi)/5 = phi/sqrt5 gives
     the Gram matrix, determinant and parity. NAMED: Conway and Sloane ch. 8 (not opened here).
  3. The single-qubit Clifford group from H and S, normalised to determinant 1, against 2O built
     independently from unit quaternions, and the moments E|tr U|^(2t), t = 1..5, against Haar SU(2).
     NAMED, citation UNVERIFIED (not opened): Zhu 2015, Webb 2016.

Integer only. A rational is a reduced (num, den) pair. Z[phi] is (a, b) = a + b phi with phi^2 = phi + 1.
Z[sqrt2] is (a, b) = a + b sqrt2. Z[zeta8] is a length-4 integer tuple, a polynomial in zeta8 reduced
mod the cyclotomic polynomial x^4 + 1. Every point set carrying a half is stored doubled, and every
halving is checked for an even remainder.

Order of work: every measured row is computed and written to MEASURE-BRIDGES.hbp first. Only then are
the frozen PREDICTION rows opened, and the COMPARE rows appended.

Run:  python3 -I -B measure_bridges.py [--repo PATH]
Exit 0: every assertion held and every COMPARE row agrees.
Exit 1: an assertion failed; a control that does not fail counts as a failed assertion.
Exit 3: every assertion held and some COMPARE row reads agree=0. That row is a result, kept as written.
"""
import sys
import os
import hashlib
import itertools
import math

if __name__ == "__main__" and not sys.flags.isolated:
    print("MEASURE_BRIDGES|needs_python3_-I=1|json=0")
    sys.exit(1)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_NAME = "MEASURE-BRIDGES.hbp"
REPO_DEFAULT = "/mnt/c/tmp/opus55/iconic-3-space/repo"
FAILS = []
ROWS = []


def check(cond, name):
    if not cond:
        FAILS.append(name)
    return 1 if cond else 0


def row(tag, pairs):
    parts = [tag]
    for k, v in pairs:
        parts.append(str(k) + "=" + str(v))
    parts.append("json=0")
    ROWS.append("|".join(parts))


# ---------------------------------------------------------------- integers and rationals

def rat(n, d):
    if d == 0:
        raise ArithmeticError("zero denominator")
    if d < 0:
        n = -n
        d = -d
    g = math.gcd(n, d)
    return (n // g, d // g)


def radd(a, b):
    return rat(a[0] * b[1] + b[0] * a[1], a[1] * b[1])


def rmul(a, b):
    return rat(a[0] * b[0], a[1] * b[1])


def rdiv(a, b):
    return rat(a[0] * b[1], a[1] * b[0])


def rstr(a):
    if a[1] == 1:
        return str(a[0])
    return str(a[0]) + "/" + str(a[1])


def power_int(b, e):
    if e < 0:
        raise ArithmeticError("negative exponent")
    r = 1
    for _ in range(e):
        r = r * b
    return r


def exact_div(x, d):
    q, r = divmod(x, d)
    if r:
        raise ArithmeticError("inexact division")
    return q


def binom(n, k):
    rowv = [1]
    for _ in range(n):
        nxt = [1]
        for i in range(1, len(rowv)):
            nxt.append(rowv[i - 1] + rowv[i])
        nxt.append(1)
        rowv = nxt
    return rowv[k]


# ---------------------------------------------------------------- integer lattice tools

def hnf(vectors, ncols):
    """Hermite normal form of the Z-span of integer vectors. Returns (basis rows, pivot columns)."""
    rows = [list(v) for v in vectors if any(v)]
    basis = []
    pivots = []
    for c in range(ncols):
        while True:
            nz = [r for r in rows if r[c] != 0]
            if len(nz) <= 1:
                break
            p = min(nz, key=lambda r: abs(r[c]))
            for r in nz:
                if r is p:
                    continue
                q = r[c] // p[c]
                if q:
                    for j in range(c, ncols):
                        r[j] -= q * p[j]
            rows = [r for r in rows if any(r)]
        nz = [r for r in rows if r[c] != 0]
        if nz:
            p = nz[0]
            rows = [r for r in rows if r is not p]
            if p[c] < 0:
                p = [-x for x in p]
            basis.append(p)
            pivots.append(c)
    if rows:
        raise ArithmeticError("hnf left a nonzero row")
    for i in range(len(basis)):
        c = pivots[i]
        for j in range(i):
            q = basis[j][c] // basis[i][c]
            if q:
                basis[j] = [x - q * y for x, y in zip(basis[j], basis[i])]
    return basis, pivots


def in_span(v, basis, pivots):
    w = list(v)
    for b, c in zip(basis, pivots):
        q, r = divmod(w[c], b[c])
        if r:
            return False
        if q:
            w = [x - q * y for x, y in zip(w, b)]
    return not any(w)


def det_bareiss(m):
    a = [list(r) for r in m]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if a[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            a[k], a[sw] = a[sw], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = exact_div(a[i][j] * a[k][k] - a[i][k] * a[k][j], prev)
        prev = a[k][k]
    return sign * a[n - 1][n - 1]


def leading_minors(g):
    return [det_bareiss([r[:k] for r in g[:k]]) for k in range(1, len(g) + 1)]


def mat_str(m):
    return ";".join(",".join(str(x) for x in r) for r in m)


def lattice_report(gram):
    d = det_bareiss(gram)
    mins = leading_minors(gram)
    diag_even = 1 if all(gram[i][i] % 2 == 0 for i in range(len(gram))) else 0
    sym = 1 if all(gram[i][j] == gram[j][i] for i in range(len(gram)) for j in range(len(gram))) else 0
    posdef = 1 if all(x > 0 for x in mins) else 0
    return d, mins, diag_even, sym, posdef


# ---------------------------------------------------------------- Gegenbauer zonal sums

def gegen(lam, kmax, t):
    """C_k^(lam)(t) for k = 0..kmax at rational t, by k C_k = 2(k+lam-1) t C_(k-1) - (k+2lam-2) C_(k-2)."""
    out = [(1, 1), rmul((2 * lam, 1), t)]
    for k in range(2, kmax + 1):
        a = rmul(rmul((2 * (k + lam - 1), 1), t), out[k - 1])
        b = rmul((k + 2 * lam - 2, 1), out[k - 2])
        out.append(rmul(radd(a, (-b[0], b[1])), (1, k)))
    return out[:kmax + 1]


def zonal_sums(dist_counts, npoints, lam, kmax):
    """S_k = (1/N) sum over ordered pairs of P_k(<x,y>), P_k = C_k^(lam)/C_k^(lam)(1). dist_counts: rational -> pair count."""
    at_one = gegen(lam, kmax, (1, 1))
    tot = [(0, 1)] * (kmax + 1)
    for t, cnt in dist_counts.items():
        vals = gegen(lam, kmax, t)
        for k in range(kmax + 1):
            tot[k] = radd(tot[k], rmul((cnt, 1), rdiv(vals[k], at_one[k])))
    return [rdiv(x, (npoints, 1)) for x in tot]


# ================================================================ PART 1: octonions to E8

def fano_lines():
    """Lines (i, i+1, i+3) mod 7 on the labels 1..7: e_i e_(i+1) = e_(i+3)."""
    return [(i, (i % 7) + 1, ((i + 2) % 7) + 1) for i in range(1, 8)]


def octonion_table(lines, flip=None):
    mt = [[None] * 8 for _ in range(8)]
    for j in range(8):
        mt[0][j] = (1, j)
        mt[j][0] = (1, j)
    for i in range(1, 8):
        mt[i][i] = (-1, 0)
    for (a, b, c) in lines:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            if mt[x][y] is not None or mt[y][x] is not None:
                raise ArithmeticError("octonion table pair set twice")
            mt[x][y] = (1, z)
            mt[y][x] = (-1, z)
    for i in range(8):
        for j in range(8):
            if mt[i][j] is None:
                raise ArithmeticError("octonion table pair unset")
    if flip is not None:
        i, j = flip
        s, k = mt[i][j]
        mt[i][j] = (-s, k)
    return mt


def sparse(v):
    return [(i, x) for i, x in enumerate(v) if x]


def omul(mt, u, v):
    w = [0] * 8
    sv = sparse(v)
    for i, a in sparse(u):
        r = mt[i]
        for j, b in sv:
            s, k = r[j]
            w[k] += s * a * b
    return w


def onorm(v):
    return sum(x * x for x in v)


def qmask(quad):
    m = 0
    for i in quad:
        m |= 1 << i
    return m


def is_s348(quads):
    if len(quads) != 14:
        return False
    for tri in itertools.combinations(range(8), 3):
        if sum(1 for q in quads if set(tri) <= q) != 1:
            return False
    return True


def code_masks(quads):
    return {0, 255} | {qmask(q) for q in quads}


def octo_units(quads):
    us = []
    for k in range(8):
        for s in (2, -2):
            v = [0] * 8
            v[k] = s
            us.append(tuple(v))
    for q in sorted(tuple(sorted(x)) for x in quads):
        for signs in itertools.product((1, -1), repeat=4):
            v = [0] * 8
            for i, s in zip(q, signs):
                v[i] = s
            us.append(tuple(v))
    return us


def octo_closure(mt, units, masks, early):
    """Doubled units U = 2u. Products W = U V = 4uv, so 2uv = W/2 must be integer with parity in the code."""
    uset = set(units)
    sp = [sparse(u) for u in units]
    total = 0
    in_order = 0
    in_units = 0
    for su in sp:
        for sv in sp:
            w = [0] * 8
            for i, a in su:
                r = mt[i]
                for j, b in sv:
                    s, k = r[j]
                    w[k] += s * a * b
            total += 1
            ok = True
            x = []
            m = 0
            for idx, c in enumerate(w):
                q, rem = divmod(c, 2)
                if rem:
                    ok = False
                    break
                x.append(q)
                if q & 1:
                    m |= 1 << idx
            if ok and m in masks:
                in_order += 1
                if tuple(x) in uset:
                    in_units += 1
                elif early:
                    return total, in_order, in_units, False
            elif early:
                return total, in_order, in_units, False
    return total, in_order, in_units, (in_units == total)


def lcg_values(seed, count):
    s = seed
    out = []
    for _ in range(count):
        s = (s * 1103515245 + 12345) & 0x7FFFFFFF
        out.append(((s >> 16) % 7) - 3)
    return out


def algebra_checks(mt, trials):
    vals = lcg_values(20261109, 16 * trials)
    nm = 0
    al = 0
    ar = 0
    for t in range(trials):
        x = vals[16 * t:16 * t + 8]
        y = vals[16 * t + 8:16 * t + 16]
        if onorm(omul(mt, x, y)) != onorm(x) * onorm(y):
            nm += 1
        if omul(mt, x, omul(mt, x, y)) != omul(mt, omul(mt, x, x), y):
            al += 1
        if omul(mt, omul(mt, y, x), x) != omul(mt, y, omul(mt, x, x)):
            ar += 1
    nonassoc = 0
    for a in range(1, 8):
        for b in range(1, 8):
            for c in range(1, 8):
                ea = [0] * 8
                eb = [0] * 8
                ec = [0] * 8
                ea[a] = 1
                eb[b] = 1
                ec[c] = 1
                if omul(mt, omul(mt, ea, eb), ec) != omul(mt, ea, omul(mt, eb, ec)):
                    nonassoc += 1
    return nm, al, ar, nonassoc


def part1():
    lines = fano_lines()
    mt = octonion_table(lines)
    row("OCT_TABLE", [("basis", "e0..e7"), ("e0", "1"), ("rule", "e_i*e_(i+1)=e_(i+3)_mod7_on_labels_1..7"),
                      ("lines", ".".join("".join(str(x) for x in ln) for ln in lines)),
                      ("orientation", "e_a*e_b=e_c_e_b*e_c=e_a_e_c*e_a=e_b_reversed_negative"), ("squares", "e_i*e_i=-1"),
                      ("tag", "STATED")])
    trials = 2000
    nm, al, ar, nonassoc = algebra_checks(mt, trials)
    check(nm == 0, "oct_norm_multiplicative")
    check(al == 0 and ar == 0, "oct_alternative")
    check(nonassoc > 0, "oct_nonassociative")
    row("OCT_ALGEBRA", [("lcg_seed", 20261109), ("trials", trials), ("entries", "-3..3"), ("norm_mult_fail", nm),
                        ("left_alt_fail", al), ("right_alt_fail", ar), ("nonassoc_imag_basis_triples", nonassoc),
                        ("of", 343), ("composition_algebra", 1 if nm == 0 else 0), ("tag", "MEASURED")])

    kirmse = [frozenset({0, a, b, c}) for (a, b, c) in lines]
    kirmse = kirmse + [frozenset(set(range(8)) - q) for q in kirmse]
    check(is_s348(kirmse), "kirmse_is_s348")
    systems = set()
    for p in itertools.permutations(range(8)):
        systems.add(frozenset(frozenset(p[i] for i in q) for q in kirmse))
    systems = sorted(systems, key=lambda s: sorted(tuple(sorted(q)) for q in s))
    check(len(systems) == 30, "s348_count_30")
    check(all(is_s348(list(s)) for s in systems), "all_s348")
    closed = []
    for s in systems:
        us = octo_units(s)
        tot, io_, iu, ok = octo_closure(mt, us, code_masks(s), True)
        if ok:
            closed.append(s)
    ktot, kio, kiu, kok = octo_closure(mt, octo_units(kirmse), code_masks(kirmse), False)
    check(not kok, "control_kirmse_not_closed")
    row("OCT_CODES", [("s348_systems_on_8_coords", len(systems)), ("closed_orders", len(closed)),
                      ("kirmse_definition", "quads_{0}+line_and_complements"), ("kirmse_closed", 1 if kok else 0),
                      ("kirmse_products", ktot), ("kirmse_products_in_set", kio), ("kirmse_products_units", kiu),
                      ("tag", "MEASURED")])
    row("CONTROL_OCT_KIRMSE", [("set", "naive_Kirmse"), ("closure", "PASS" if kok else "FAIL"),
                               ("products_outside_set", ktot - kio), ("control_failed_as_required", 1 if not kok else 0),
                               ("tag", "MEASURED")])
    trans = []
    trans_images = set()
    chosen = None
    chosen_name = ""
    for a, b in itertools.combinations(range(8), 2):
        img = frozenset(frozenset((b if i == a else a if i == b else i) for i in q) for q in kirmse)
        if img in closed:
            trans.append(str(a) + str(b))
            trans_images.add(img)
    for k in range(1, 8):
        img = frozenset(frozenset((k if i == 0 else 0 if i == k else i) for i in q) for q in kirmse)
        if img in closed and chosen is None:
            chosen = img
            chosen_name = "Kirmse_with_coordinates_0_and_" + str(k) + "_interchanged"
    if chosen is None and closed:
        chosen = closed[0]
        chosen_name = "first_closed_system_in_sorted_order"
    row("OCT_CORRECTION", [("transpositions_tested", 28), ("transpositions_giving_closed_order", len(trans)),
                           ("which", ".".join(trans) if trans else "none"),
                           ("distinct_closed_orders_reached", len(trans_images)), ("chosen", chosen_name),
                           ("coxeter_corrected", 1), ("tag", "MEASURED")])
    if chosen is None:
        check(False, "no_closed_order")
        return None
    quads = sorted(tuple(sorted(q)) for q in chosen)
    units = octo_units(chosen)
    masks = code_masks(chosen)
    row("OCT_ORDER", [("definition", "O={x/2:x_in_Z^8,x_mod_2_in_C}"), ("C", "extended_Hamming_8_4_4_quads_below"),
                      ("quads", ".".join("".join(str(i) for i in q) for q in quads)),
                      ("coords", "0=real,1..7=e1..e7"), ("tag", "STATED")])
    n1 = sum(1 for u in units if onorm(u) == 4)
    check(len(set(units)) == 240 and n1 == 240, "oct_240_units_norm1")
    tot, io_, iu, ok = octo_closure(mt, units, masks, False)
    check(ok and io_ == tot and tot == 57600, "oct_closure")
    row("OCT_UNITS", [("count", len(set(units))), ("norm_1", n1), ("of_form_+-e_k", 16), ("of_form_half_quad", 224),
                      ("doubled_norm", 4), ("tag", "MEASURED")])
    row("OCT_CLOSURE", [("products", tot), ("in_order", io_), ("in_unit_set", iu), ("closed", 1 if ok else 0),
                        ("tag", "MEASURED")])
    want = {4: 1, 2: 56, 0: 126, -2: 56, -4: 1}
    all_ok = True
    pair_dist = {}
    per_unit = None
    for u in units:
        cnt = {}
        for v in units:
            d = sum(a * b for a, b in zip(u, v))
            cnt[d] = cnt.get(d, 0) + 1
        if cnt != want:
            all_ok = False
        if per_unit is None:
            per_unit = {rat(d, 4): c for d, c in cnt.items()}
        for d, c in cnt.items():
            key = rat(d, 4)
            pair_dist[key] = pair_dist.get(key, 0) + c
    check(all_ok, "oct_e8_distribution_every_unit")
    row("OCT_E8_DIST", [("units_checked", 240), ("cos_1", want[4]), ("cos_1/2", want[2]), ("cos_0", want[0]),
                        ("cos_-1/2", want[-2]), ("cos_-1", want[-4]), ("every_unit_matches", 1 if all_ok else 0),
                        ("measured_pairs", ",".join(rstr(k) + ":" + str(pair_dist[k]) for k in sorted(pair_dist, key=lambda r: (-r[0] * 4 // r[1])))),
                        ("tag", "MEASURED")])
    basis, piv = hnf(units, 8)
    check(len(basis) == 8, "oct_rank8")
    check(all(in_span(u, basis, piv) for u in units), "oct_units_in_basis_span")
    idx = abs(det_bareiss(basis))
    gram = []
    for x in basis:
        gr = []
        for y in basis:
            gr.append(exact_div(sum(a * b for a, b in zip(x, y)), 2))
        gram.append(gr)
    d, mins, de, sym, pd = lattice_report(gram)
    check(d == 1 and de == 1 and sym == 1 and pd == 1, "oct_even_unimodular")
    row("OCT_LATTICE", [("basis_doubled", mat_str(basis)), ("index_of_doubled_order_in_Z^8", idx),
                        ("form", "B(x,y)=2<x,y>=N(x+y)-N(x)-N(y)"), ("gram", mat_str(gram)), ("det", d),
                        ("leading_minors", ",".join(str(x) for x in mins)), ("diag_even", de), ("symmetric", sym),
                        ("positive_definite", pd), ("rank", len(basis)),
                        ("even_unimodular", 1 if (d == 1 and de == 1 and pd == 1) else 0), ("tag", "MEASURED")])
    row("OCT_E8_NAMED", [("claim", "even_unimodular_rank8_positive_definite_lattice_is_E8_unique"),
                         ("cite", "Coxeter_1946_Integral_Cayley_numbers;Conway_Smith_ch9;Conway_Sloane_SPLAG"),
                         ("cite_opened", 0), ("tag", "NAMED")])
    s = zonal_sums(pair_dist, 240, 3, 30)
    tdes = 0
    for k in range(1, 31):
        if s[k][0] == 0:
            tdes = k
        else:
            break
    nonzero = [k for k in range(1, 31) if s[k][0] != 0]
    row("OCT_E8_ZONAL", [("points", 240), ("sphere", "S7"), ("gegenbauer_lambda", 3),
                         ("S_k_def", "(1/N)sum_pairs_C_k(<x,y>)/C_k(1)"), ("t_design", tdes),
                         ("first_break", ("S_" + str(tdes + 1) + "=" + rstr(s[tdes + 1])) if tdes < 30 else "none"),
                         ("nonzero_k_through_30", ",".join(str(k) for k in nonzero)),
                         ("S_1..30", ",".join(rstr(s[k]) for k in range(1, 31))), ("S_0", rstr(s[0])),
                         ("tag", "MEASURED")])

    # Controls on the octonion side.
    flip = (1, 2)
    mtb = octonion_table(lines, flip)
    nmb, alb, arb, _ = algebra_checks(mtb, trials)
    btot, bio, biu, bok = octo_closure(mtb, units, masks, False)
    check(not bok, "control_flip_closure_fails")
    row("CONTROL_OCT_FLIP", [("flip", "e1*e2:+e4->-e4_single_entry"), ("norm_mult_fail", nmb), ("of", trials),
                             ("products", btot), ("in_order", bio), ("in_unit_set", biu),
                             ("closure", "PASS" if bok else "FAIL"), ("det_unchanged", d),
                             ("note", "the_lattice_and_its_norm_form_do_not_read_the_table_so_det_stays"),
                             ("control_failed_as_required", 1 if not bok else 0), ("tag", "MEASURED")])
    graves = [u for u in units if sum(1 for x in u if x) == 1]
    gb, gp = hnf(graves, 8)
    gg = [[exact_div(sum(a * b for a, b in zip(x, y)), 2) for y in gb] for x in gb]
    gd = det_bareiss(gg)
    check(gd != 1, "control_graves_det_not_1")
    row("CONTROL_OCT_GRAVES", [("lattice", "Z-span_of_16_units_+-e_k"), ("rank", len(gb)), ("det", gd),
                               ("unimodular", 1 if gd == 1 else 0), ("control_failed_as_required", 1 if gd != 1 else 0),
                               ("tag", "MEASURED")])
    return {"zonal": s, "t": tdes, "dist": pair_dist, "det": d, "per_unit": per_unit, "every": 1 if all_ok else 0,
            "n": len(set(units)), "antipodal": 1 if (all_ok and per_unit.get((-1, 1)) == 1) else 0}


# ================================================================ PART 2: icosians to E8

def p_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def p_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def phi_mul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def r2_mul(x, y):
    a, b = x
    c, d = y
    return (a * c + 2 * b * d, a * d + b * c)


def qmul(p, q, mul):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    w = p_sub(p_sub(p_sub(mul(w1, w2), mul(x1, x2)), mul(y1, y2)), mul(z1, z2))
    x = p_sub(p_add(p_add(mul(w1, x2), mul(x1, w2)), mul(y1, z2)), mul(z1, y2))
    y = p_add(p_add(p_sub(mul(w1, y2), mul(x1, z2)), mul(y1, w2)), mul(z1, x2))
    z = p_add(p_sub(p_add(mul(w1, z2), mul(x1, y2)), mul(y1, x2)), mul(z1, w2))
    return (w, x, y, z)


def half_pair(x):
    return (exact_div(x[0], 2), exact_div(x[1], 2))


def qmul_doubled(p, q, mul):
    """p = 2a, q = 2b, returns 2ab = pq/2."""
    r = qmul(p, q, mul)
    return tuple(half_pair(c) for c in r)


def qnorm(p, mul):
    t = (0, 0)
    for c in p:
        t = p_add(t, mul(c, c))
    return t


def perm_parity(p):
    inv = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return inv % 2


def icosians():
    """Doubled coordinates over Z[phi]: 2q = (w, x, y, z) for q = w + x i + y j + z k."""
    zero = (0, 0)
    out = []
    for k in range(4):
        for s in (2, -2):
            v = [zero] * 4
            v[k] = (s, 0)
            out.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        out.append(tuple((s, 0) for s in signs))
    tmpl = [(0, 0), (1, 0), (-1, 1), (0, 1)]
    for p in itertools.permutations(range(4)):
        if perm_parity(p):
            continue
        for signs in itertools.product((1, -1), repeat=3):
            v = [zero] * 4
            for t in range(4):
                val = tmpl[t]
                if t >= 1:
                    s = signs[t - 1]
                    val = (s * val[0], s * val[1])
                v[p[t]] = val
            out.append(tuple(v))
    return out


def ico_vec(q):
    return [c for pair in q for c in pair]


def vec_ico(v):
    return tuple((v[2 * i], v[2 * i + 1]) for i in range(4))


def ico_bform(x, y, lam):
    """B(x,y) = Tr(lam * 2<x,y>) for doubled x = 2a, y = 2b; <a,b> = sum a_m b_m in Z[phi]; lam = (p, q, den) = (p + q phi)/den."""
    t = (0, 0)
    for a, b in zip(x, y):
        t = p_add(t, phi_mul(a, b))
    two_inner = half_pair(t)
    u, v = phi_mul((lam[0], lam[1]), two_inner)
    return exact_div(2 * u + v, lam[2])


def part2():
    units = icosians()
    uset = set(units)
    check(len(uset) == 120 and len(units) == 120, "ico_120_distinct")
    n1 = sum(1 for q in units if qnorm(q, phi_mul) == (4, 0))
    check(n1 == 120, "ico_norm1")
    prod_ok = 0
    for a in units:
        for b in units:
            if qmul_doubled(a, b, phi_mul) in uset:
                prod_ok += 1
    check(prod_ok == 14400, "ico_group_closed")
    ident = ((2, 0), (0, 0), (0, 0), (0, 0))
    inv_ok = sum(1 for q in units if (q[0], (-q[1][0], -q[1][1]), (-q[2][0], -q[2][1]), (-q[3][0], -q[3][1])) in uset)
    check(ident in uset and inv_ok == 120, "ico_identity_inverses")
    row("ICO_UNITS", [("definition", "24_Hurwitz_units+96_even_perms_of_(0,+-1,+-phi^-1,+-phi)/2"),
                      ("count", len(uset)), ("norm_1", n1), ("products", 14400), ("products_in_set", prod_ok),
                      ("identity", 1 if ident in uset else 0), ("conjugates_in_set", inv_ok),
                      ("group_2I", 1 if prod_ok == 14400 else 0), ("tag", "MEASURED")])
    phiv = (0, 1)
    gens_units = [ico_vec(q) for q in units]
    gens_phi = [ico_vec(tuple(phi_mul(phiv, c) for c in q)) for q in units]
    basis, piv = hnf(gens_units + gens_phi, 8)
    basis_u, piv_u = hnf(gens_units, 8)
    same = 1 if (basis == basis_u and piv == piv_u) else 0
    check(len(basis) == 8, "ico_rank8")
    check(same == 1, "ico_zspan_equals_zphi_span")
    check(all(in_span(v, basis, piv) for v in gens_units + gens_phi), "ico_gens_in_span")
    ring_closed = 1
    bq = [vec_ico(b) for b in basis]
    for a in bq:
        for b in bq:
            if not in_span(ico_vec(qmul_doubled(a, b, phi_mul)), basis, piv):
                ring_closed = 0
    check(ring_closed == 1, "ico_ring_closed")
    lam = (2, 1, 5)
    gram = [[ico_bform(x, y, lam) for y in bq] for x in bq]
    d, mins, de, sym, pd = lattice_report(gram)
    check(d == 1 and de == 1 and sym == 1 and pd == 1, "ico_even_unimodular")
    row("ICO_RING", [("span", "Z[phi]-span_of_120_icosians"), ("generators", "u_and_phi*u"),
                     ("coords", "2q=(w,x,y,z)_each_a+b*phi_as_a,b"), ("rank_over_Z", len(basis)),
                     ("zspan_of_units_equals_zphi_span", same), ("basis_products_in_span", ring_closed),
                     ("basis_doubled", mat_str(basis)), ("tag", "MEASURED")])
    row("ICO_FORM", [("lambda", "(2+phi)/5=phi/sqrt5"), ("Q", "Q(x)=Tr(lambda*N(x))=a+b_where_N(x)=a+b*phi"),
                     ("B", "B(x,y)=Q(x+y)-Q(x)-Q(y)=Tr(lambda*2<x,y>)"), ("Q_of_unit", 1),
                     ("cite", "Conway_Sloane_SPLAG_ch8;Wilson_Finite_Simple_Groups"), ("cite_opened", 0),
                     ("tag", "STATED")])
    row("ICO_LATTICE", [("gram", mat_str(gram)), ("det", d), ("leading_minors", ",".join(str(x) for x in mins)),
                        ("diag_even", de), ("symmetric", sym), ("positive_definite", pd), ("rank", len(gram)),
                        ("even_unimodular", 1 if (d == 1 and de == 1 and pd == 1) else 0), ("tag", "MEASURED")])
    row("ICO_E8_NAMED", [("claim", "even_unimodular_rank8_hence_E8_by_uniqueness"), ("cite", "Conway_Sloane_SPLAG"),
                         ("cite_opened", 0), ("tag", "NAMED")])
    phinv = (-1, 1)
    roots = units + [tuple(phi_mul(phinv, c) for c in q) for q in units]
    q1 = 0
    for r in roots:
        n = qnorm(r, phi_mul)
        nq = half_pair(half_pair(n))
        u, v = phi_mul((lam[0], lam[1]), nq)
        if exact_div(2 * u + v, lam[2]) == 1:
            q1 += 1
    check(len(set(roots)) == 240 and q1 == 240, "ico_240_Q1")
    want = {2: 1, 1: 56, 0: 126, -1: 56, -2: 1}
    every = 1
    per_root = None
    for x in roots:
        cnt = {}
        for y in roots:
            b = ico_bform(x, y, lam)
            cnt[b] = cnt.get(b, 0) + 1
        if cnt != want:
            every = 0
        if per_root is None:
            per_root = {rat(b, 2): c for b, c in cnt.items()}
    check(every == 1, "ico_roots_e8_distribution")
    row("ICO_ROOTS", [("vectors", "units_and_(phi-1)*units"), ("count", len(set(roots))), ("Q_equals_1", q1),
                      ("B_2", 1), ("B_1", 56), ("B_0", 126), ("B_-1", 56), ("B_-2", 1),
                      ("every_vector_matches", every), ("tag", "MEASURED")])
    lam_c = (1, 0, 1)
    gram_c = [[ico_bform(x, y, lam_c) for y in bq] for x in bq]
    dc = det_bareiss(gram_c)
    check(dc != 1, "control_trace_form_not_unimodular")
    row("CONTROL_ICO_TRACE", [("lambda", "1"), ("form", "Tr(N(x))"), ("det", dc), ("unimodular", 1 if dc == 1 else 0),
                              ("control_failed_as_required", 1 if dc != 1 else 0), ("tag", "MEASURED")])
    return {"det": d, "per_unit": per_root, "every": every, "n": len(set(roots)),
            "antipodal": 1 if (every == 1 and per_root.get((-1, 1)) == 1) else 0}


# ================================================================ PART 3: the Clifford group and 2O

def z_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def z_neg(a):
    return tuple(-x for x in a)


def z_mul(a, b):
    c = [0] * 7
    for i in range(4):
        if a[i]:
            for j in range(4):
                c[i + j] += a[i] * b[j]
    return (c[0] - c[4], c[1] - c[5], c[2] - c[6], c[3])


def z_conj(a):
    """zeta -> zeta^-1 = -zeta^3."""
    return (a[0], -a[3], -a[2], -a[1])


def z_half(a):
    return tuple(exact_div(x, 2) for x in a)


Z0 = (0, 0, 0, 0)
ZETA = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1)]


def m_mul(a, b):
    """2x2 matrices stored doubled (2M), entries in Z[zeta8]: (2A)(2B)/2 = 2AB."""
    a00, a01, a10, a11 = a
    b00, b01, b10, b11 = b
    return (z_half(z_add(z_mul(a00, b00), z_mul(a01, b10))), z_half(z_add(z_mul(a00, b01), z_mul(a01, b11))),
            z_half(z_add(z_mul(a10, b00), z_mul(a11, b10))), z_half(z_add(z_mul(a10, b01), z_mul(a11, b11))))


def m_det(a):
    """det M for doubled a = 2M: det(2M)/4."""
    d = z_add(z_mul(a[0], a[3]), z_neg(z_mul(a[1], a[2])))
    return tuple(exact_div(x, 4) for x in d)


def m_scale(a, z):
    return tuple(z_mul(e, z) for e in a)


IDENT2 = ((2, 0, 0, 0), Z0, Z0, (2, 0, 0, 0))
SQRT2 = (0, 1, 0, -1)


def generate(gens, limit):
    seen = {IDENT2}
    frontier = [IDENT2]
    while frontier:
        nxt = []
        for m in frontier:
            for g in gens:
                p = m_mul(m, g)
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
        frontier = nxt
        if len(seen) > limit:
            raise ArithmeticError("group larger than limit")
    return seen


def normalise_su2(group):
    out = set()
    bad = 0
    for m in group:
        dm = m_det(m)
        if dm not in ZETA:
            bad += 1
            continue
        k = ZETA.index(dm)
        if k % 2:
            bad += 1
            continue
        m0 = (-(k // 2)) % 4
        for mm in (m0, m0 + 4):
            if (2 * mm + k) % 8:
                bad += 1
            out.add(m_scale(m, ZETA[mm]))
    return out, bad


def in_zi_sqrt2(m):
    """2M entries in Z[i, sqrt2] = Z + Z i + Z sqrt2 + Z i sqrt2: zeta and zeta^3 coefficients of equal parity."""
    return all((e[1] - e[3]) % 2 == 0 for e in m)


def binary_octahedral_quats():
    """Doubled coordinates over Z[sqrt2]: 2q, each coordinate (a, b) = a + b sqrt2."""
    zero = (0, 0)
    out = []
    for k in range(4):
        for s in (2, -2):
            v = [zero] * 4
            v[k] = (s, 0)
            out.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        out.append(tuple((s, 0) for s in signs))
    for a, b in itertools.combinations(range(4), 2):
        for sa in (1, -1):
            for sb in (1, -1):
                v = [zero] * 4
                v[a] = (0, sa)
                v[b] = (0, sb)
                out.append(tuple(v))
    return out


def r2_to_z(x):
    """a + b sqrt2 into Z[zeta8], sqrt2 = zeta - zeta^3."""
    return (x[0], x[1], 0, -x[1])


def r2_to_iz(x):
    """i (a + b sqrt2) = a zeta^2 + b (zeta + zeta^3)."""
    return (0, x[1], x[0], x[1])


def quat_to_su2(q):
    """q = w + x i + y j + z k -> [[w + x i, y + z i], [-y + z i, w - x i]] (doubled in, doubled out)."""
    w, x, y, z = q
    return (z_add(r2_to_z(w), r2_to_iz(x)), z_add(r2_to_z(y), r2_to_iz(z)),
            z_add(z_neg(r2_to_z(y)), r2_to_iz(z)), z_add(r2_to_z(w), z_neg(r2_to_iz(x))))


def moments(group, nmax, tmax):
    """E[(tr U)^n], n = 1..nmax, and E|tr U|^(2t), t = 1..tmax, exact rationals over a doubled group."""
    g = len(group)
    tr_pows = [Z0] * (nmax + 1)
    abs_pows = [Z0] * (tmax + 1)
    real = 1
    for m in group:
        t = z_add(m[0], m[3])
        if not (t[2] == 0 and t[1] == -t[3]):
            real = 0
        p = (1, 0, 0, 0)
        for n in range(1, nmax + 1):
            p = z_mul(p, t)
            tr_pows[n] = z_add(tr_pows[n], p)
        a2 = z_mul(t, z_conj(t))
        p = (1, 0, 0, 0)
        for k in range(1, tmax + 1):
            p = z_mul(p, a2)
            abs_pows[k] = z_add(abs_pows[k], p)
    rat_tr = 1 if all(v[1] == 0 and v[2] == 0 and v[3] == 0 for v in tr_pows[1:]) else 0
    rat_abs = 1 if all(v[1] == 0 and v[2] == 0 and v[3] == 0 for v in abs_pows[1:]) else 0
    em = [rat(tr_pows[n][0], g * power_int(2, n)) for n in range(1, nmax + 1)] if rat_tr else []
    ea = [rat(abs_pows[k][0], g * power_int(4, k)) for k in range(1, tmax + 1)] if rat_abs else []
    return em, ea, real, rat_tr, rat_abs


def haar_su2(n):
    """Multiplicity of the trivial representation in V_(1/2)^(tensor n): walks on 2j = 0,1,2,... from 0 to 0."""
    dist = {0: 1}
    for _ in range(n):
        nd = {}
        for pos, c in dist.items():
            for q in (pos + 1, pos - 1):
                if q >= 0:
                    nd[q] = nd.get(q, 0) + c
        dist = nd
    return dist.get(0, 0)


def part3():
    h2 = (SQRT2, SQRT2, SQRT2, z_neg(SQRT2))
    s2 = ((2, 0, 0, 0), Z0, Z0, (0, 0, 2, 0))
    row("CLIFF_GENERATORS", [("H", "(1/sqrt2)[[1,1],[1,-1]]"), ("S", "diag(1,i)"),
                             ("stored", "2M_entries_in_Z[zeta8]_mod_x^4+1"), ("sqrt2", "zeta-zeta^3"),
                             ("tag", "STATED")])
    u2 = generate([h2, s2], 5000)
    dets = set(ZETA.index(m_det(m)) if m_det(m) in ZETA else -1 for m in u2)
    su2, bad = normalise_su2(u2)
    check(bad == 0, "cliff_normalisation_exact")
    det1 = sum(1 for m in su2 if m_det(m) == ZETA[0])
    zi = sum(1 for m in su2 if in_zi_sqrt2(m))
    closed = all(m_mul(a, b) in su2 for a in su2 for b in su2)
    check(len(su2) == 48 and det1 == 48 and zi == 48 and closed, "cliff_su2_48")
    row("CLIFF_GROUP", [("U2_order_generated", len(u2)), ("det_exponents_of_zeta8", ",".join(str(x) for x in sorted(dets))),
                        ("normalised_to_SU2", len(su2)), ("det_1", det1), ("entries_2M_in_Z[i,sqrt2]", zi),
                        ("normalised_set_closed", 1 if closed else 0), ("normalisation_failures", bad),
                        ("tag", "MEASURED")])
    bo = binary_octahedral_quats()
    boset = set(bo)
    bo_closed = sum(1 for a in bo for b in bo if qmul_doubled(a, b, r2_mul) in boset)
    bo_norm = sum(1 for q in bo if qnorm(q, r2_mul) == (4, 0))
    check(len(boset) == 48 and bo_closed == 2304 and bo_norm == 48, "2O_quaternion_group")
    img = {}
    for q in bo:
        img[q] = quat_to_su2(q)
    hom = sum(1 for a in bo for b in bo if quat_to_su2(qmul_doubled(a, b, r2_mul)) == m_mul(img[a], img[b]))
    check(hom == 2304, "2O_map_homomorphism")
    imgset = set(img.values())
    equal = 1 if imgset == su2 else 0
    check(equal == 1, "cliff_equals_2O")
    row("BO_QUATERNIONS", [("definition", "24_Hurwitz_units+24_(+-a+-b)/sqrt2_a!=b_in_{1,i,j,k}"),
                           ("count", len(boset)), ("norm_1", bo_norm), ("products_in_set", bo_closed), ("of", 2304),
                           ("map", "w+xi+yj+zk->[[w+xi,y+zi],[-y+zi,w-xi]]"), ("homomorphism_pairs", hom),
                           ("image_size", len(imgset)), ("tag", "MEASURED")])
    row("CLIFF_EQ_2O", [("clifford_su2", len(su2)), ("2O_image", len(imgset)), ("equal_as_sets", equal),
                        ("tag", "MEASURED")])
    nmax = 10
    tmax = 5
    em, ea, real, rt, ra = moments(su2, nmax, tmax)
    _, ea_u2, _, _, ra_u2 = moments(u2, nmax, tmax)
    check(real == 1 and rt == 1 and ra == 1 and ra_u2 == 1, "cliff_moments_exact_rational")
    haar_n = [haar_su2(n) for n in range(1, nmax + 1)]
    haar_t = [haar_su2(2 * t) for t in range(1, tmax + 1)]
    cat_ok = all(haar_t[t - 1] == exact_div(binom(2 * t, t), t + 1) for t in range(1, tmax + 1))
    check(cat_ok, "haar_catalan")
    first = 0
    for t in range(1, tmax + 1):
        if ea[t - 1] != (haar_t[t - 1], 1):
            first = t
            break
    check(first == 4, "cliff_first_difference_t4")
    check(ea == ea_u2, "cliff_phase_invariance")
    row("CLIFF_MOMENTS", [("quantity", "E_abs_tr_U^(2t)"), ("t", "1,2,3,4,5"), ("clifford_2O", ",".join(rstr(x) for x in ea)),
                          ("clifford_U2_192", ",".join(rstr(x) for x in ea_u2)), ("haar_SU2", ",".join(str(x) for x in haar_t)),
                          ("haar_method", "trivial_multiplicity_in_V_half^(2t)_by_walk_count"), ("haar_equals_catalan", 1 if cat_ok else 0),
                          ("first_t_differs", first), ("unitary_t_design", (first - 1) if first else tmax),
                          ("tag", "MEASURED")])
    row("CLIFF_TRACE_MOMENTS", [("quantity", "M_n=E[(tr_U)^n]_over_2O_in_SU2"), ("n", "1..10"),
                                ("M", ",".join(rstr(x) for x in em)), ("haar_SU2", ",".join(str(x) for x in haar_n)),
                                ("trace_real", real), ("tag", "MEASURED")])
    row("CLIFF_DESIGN_NAMED", [("claim", "single_qubit_Clifford_group_is_unitary_3_design_not_4_design"),
                               ("cite", "Zhu_2015;Webb_2016"), ("cite_opened", 0), ("cite_status", "UNVERIFIED"),
                               ("tag", "NAMED")])
    zc = ((2, 0, 0, 0), Z0, Z0, (-2, 0, 0, 0))
    u2c = generate([h2, zc], 5000)
    suc, badc = normalise_su2(u2c)
    eqc = 1 if suc == su2 else 0
    _, eac, _, _, _ = moments(suc, nmax, tmax)
    firstc = 0
    for t in range(1, tmax + 1):
        if eac[t - 1] != (haar_t[t - 1], 1):
            firstc = t
            break
    check(len(suc) < 48 and eqc == 0, "control_S_to_Z_shrinks")
    row("CONTROL_CLIFF_Z", [("S_replaced_by", "diag(1,-1)"), ("U2_order_generated", len(u2c)),
                            ("normalised_to_SU2", len(suc)), ("normalisation_failures", badc), ("equal_2O", eqc), ("moments_E_abs_tr^2t", ",".join(rstr(x) for x in eac)),
                            ("first_t_differs", firstc), ("control_failed_as_required", 1 if (len(suc) < 48 and eqc == 0) else 0),
                            ("tag", "MEASURED")])
    return {"M": em, "abs": ea, "first": first, "haar_n": haar_n}


# ================================================================ compare (after the measured rows are written)

FROZEN = {
    "su2": ("predictions/su2-ladder/PREDICTION.hbp", "c2601cf0"),
    "exc": ("predictions/exceptional/PREDICTION-EXCEPTIONAL.hbp", "e8403cd4"),
}


def parse_rows(path):
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("|")
            kv = {}
            for p in parts[1:]:
                if "=" in p:
                    k, v = p.split("=", 1)
                    kv[k] = v
            out.append((parts[0], kv))
    return out


def parse_q(t):
    if ":" in t:
        n, d = t.split(":")
        return rat(int(n), int(d))
    return (int(t), 1)


def parse_q5(t):
    """a+br5 with a and b written num:den, returned as the pair of rationals (a, b)."""
    if not t.endswith("r5"):
        raise ArithmeticError("not a Q(sqrt5) value")
    i = t.index("+")
    return parse_q(t[:i]), parse_q(t[i + 1:-2])


def find(rows, tag, key, val):
    hits = [kv for (tg, kv) in rows if tg == tag and kv.get(key) == val]
    if len(hits) != 1:
        raise ArithmeticError("frozen row count for " + tag + " " + val + " is " + str(len(hits)))
    return hits[0]


def compare_section(repo, meas):
    out = []
    dis = []

    def cmp(quantity, frozen_str, measured_str, agree):
        out.append("COMPARE|quantity=" + quantity + "|frozen=" + frozen_str + "|measured=" + measured_str +
                   "|agree=" + str(1 if agree else 0) + "|tag=MEASURED_vs_PREDICTION|json=0")
        if not agree:
            dis.append(quantity)

    files = {}
    for key, (rel, prefix) in FROZEN.items():
        path = os.path.join(repo, rel)
        got = sha256_file(path)
        with open(path + ".sha256", "r", encoding="utf-8") as f:
            side = f.read().split()[0]
        ok = check(got == side and got.startswith(prefix), "frozen_sha_" + key)
        out.append("FROZEN_FILE|file=" + rel + "|sha256=" + got + "|sidecar=" + side + "|expected_prefix=" + prefix +
                   "|match=" + str(ok) + "|opened_after_measured_section_written=1|json=0")
        files[key] = parse_rows(path)
    su2 = files["su2"]
    exc = files["exc"]

    # 1. the frozen 2O moment row, M_2..M_10, against the Clifford group normalised into SU(2)
    mo = find(su2, "MOMENT", "group", "2O")
    check(mo.get("n_range") == "0..60" and mo.get("quantity") == "E_G[chi_half^n]", "frozen_moment_format")
    fm = [parse_q(x) for x in mo["M"].split(",")]
    cl = meas["cliff"]
    for n in range(2, 11):
        mv = cl["M"][n - 1]
        cmp("2O_M_" + str(n) + "_E[(trU)^n]_clifford", rstr(fm[n]), rstr(mv), fm[n] == mv)
    mh = find(su2, "MOMENT_HAAR", "group", "SU2_Haar")
    fh = [parse_q(x) for x in mh["Mh"].split(",")]
    for n in range(2, 11):
        hv = (cl["haar_n"][n - 1], 1)
        cmp("SU2_Haar_M_" + str(n), rstr(fh[n]), rstr(hv), fh[n] == hv)
    mb = find(su2, "MOMENT_BREAK", "group", "2O")
    first_n = 0
    for n in range(1, 11):
        if cl["M"][n - 1] != (cl["haar_n"][n - 1], 1):
            first_n = n
            break
    cmp("2O_first_n_with_M_ne_Mh", mb["first_n_with_M_ne_Mh"], str(first_n), int(mb["first_n_with_M_ne_Mh"]) == first_n)
    if first_n:
        cmp("2O_M_at_first_break", mb["M"], rstr(cl["M"][first_n - 1]), parse_q(mb["M"]) == cl["M"][first_n - 1])
        cmp("SU2_Mh_at_first_break", mb["Mh"], str(cl["haar_n"][first_n - 1]), parse_q(mb["Mh"]) == (cl["haar_n"][first_n - 1], 1))
    de = find(su2, "DESIGN", "group", "2O")
    cmp("2O_first_break_degree_vs_2*first_t_of_E_abs_trU^2t", de["first_break_degree"], str(2 * cl["first"]),
        int(de["first_break_degree"]) == 2 * cl["first"])
    cmp("2O_order_vs_clifford_su2_size", de["order"], "48", int(de["order"]) == 48)

    # 2. the frozen E8 SET row, against the octonion-built and the icosian-built root sets
    st = find(exc, "SET", "set", "E8_roots")
    fdist = {}
    for ent in st["dist"].split(","):
        c, cnt = ent.split("#")
        a, b = parse_q5(c)
        if b != (0, 1):
            raise ArithmeticError("E8 cosine with a sqrt5 part")
        fdist[a] = fdist.get(a, 0) + int(cnt)
    fdist_s = ",".join(rstr(k) + ":" + str(fdist[k]) for k in sorted(fdist, key=lambda r: -r[0] * 4 // r[1]))
    for label, m in (("octonion", meas["oct"]), ("icosian", meas["ico"])):
        md = m["per_unit"]
        md_s = ",".join(rstr(k) + ":" + str(md[k]) for k in sorted(md, key=lambda r: -r[0] * 4 // r[1]))
        cmp("E8_" + label + "_N", st["N"], str(m["n"]), int(st["N"]) == m["n"])
        cmp("E8_" + label + "_distinct_cosines", st["rows"], str(len(md)), int(st["rows"]) == len(md))
        cmp("E8_" + label + "_counts_sum", st["counts_sum"], str(sum(md.values())), int(st["counts_sum"]) == sum(md.values()))
        cmp("E8_" + label + "_distribution_from_every_point", fdist_s, md_s + ("" if m["every"] else "_NOT_EVERY_POINT"),
            fdist == md and m["every"] == 1)
        cmp("E8_" + label + "_antipodal", st["antipodal"], str(m["antipodal"]), int(st["antipodal"]) == m["antipodal"])

    # 3. the frozen E8 zonal sums S_k, k = 1..30, against the octonion-built root set
    oz = meas["oct"]["zonal"]
    for k in range(1, 31):
        sk = [kv for (tg, kv) in exc if tg == "SK" and kv.get("set") == "E8_roots" and kv.get("k") == str(k)]
        if len(sk) != 1:
            raise ArithmeticError("frozen SK row missing at k=" + str(k))
        fa = parse_q(sk[0]["S_a"])
        fb = parse_q(sk[0]["S_b"])
        cmp("E8_octonion_S_" + str(k), sk[0]["S_a"] + "+" + sk[0]["S_b"] + "r5", rstr(oz[k]), fa == oz[k] and fb == (0, 1))
    tp = find(exc, "TPRED", "set", "E8_roots")
    cmp("E8_octonion_t_design", tp["t_pred"], str(meas["oct"]["t"]), int(tp["t_pred"]) == meas["oct"]["t"])
    fbk = find(exc, "FIRST_BREAK", "set", "E8_roots")
    tb = meas["oct"]["t"] + 1
    cmp("E8_octonion_first_break_k", fbk["k"], str(tb), int(fbk["k"]) == tb)
    cmp("E8_octonion_first_break_S", fbk["S_a"], rstr(oz[tb]), parse_q(fbk["S_a"]) == oz[tb] and parse_q(fbk["S_b"]) == (0, 1))
    out.append("COMPARE_END|rows=" + str(len([r for r in out if r.startswith("COMPARE|")])) + "|disagree=" + str(len(dis)) +
               "|which=" + (",".join(dis) if dis else "none") + "|json=0")
    return out, len(dis)


def row_ok(r):
    """TAG|k=v|...|json=0: a tag without '=', then every field holds '=', last field json=0."""
    parts = r.split("|")
    if len(parts) < 2 or "=" in parts[0] or not parts[0] or parts[-1] != "json=0":
        return False
    return all("=" in x and not x.startswith("=") for x in parts[1:])


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    repo = REPO_DEFAULT
    if "--repo" in sys.argv:
        repo = sys.argv[sys.argv.index("--repo") + 1]
    script_sha = sha256_file(os.path.abspath(__file__))
    meas = {}
    try:
        meas["oct"] = part1()
        meas["ico"] = part2()
        meas["cliff"] = part3()
    except ArithmeticError as e:
        check(False, "exception_" + str(e).replace(" ", "_"))
    bad_rows = [r for r in ROWS if not row_ok(r)]
    check(not bad_rows, "row_format_measured")
    measured = list(ROWS)
    sec = hashlib.sha256(("\n".join(measured) + "\n").encode("utf-8")).hexdigest()
    head = ("MEASURE_BRIDGES|script=measure_bridges.py|script_sha256=" + script_sha +
            "|python=" + ".".join(str(x) for x in sys.version_info[:3]) +
            "|integer_only=1|E=0|order=measure_first_then_compare|json=0")
    out = os.path.join(HERE, OUT_NAME)
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write(head + "\n")
        for r in measured:
            f.write(r + "\n")
        f.write("MEASURED_SECTION|rows=" + str(len(measured)) + "|sha256=" + sec + "|assert_failures=" + str(len(FAILS)) +
                "|written_before_predictions_opened=1|json=0\n")
    ROWS.clear()
    cmp_rows = []
    disagree = 0
    if not FAILS and "oct" in meas and "cliff" in meas:
        cmp_rows, disagree = compare_section(repo, meas)
    check(all(row_ok(r) for r in cmp_rows), "row_format_compare")
    n_cmp = sum(1 for r in cmp_rows if r.startswith("COMPARE|"))
    with open(out, "a", encoding="utf-8", newline="\n") as f:
        for r in cmp_rows:
            f.write(r + "\n")
        status = "FAIL" if FAILS else ("DISAGREE" if disagree else "PASS")
        f.write("SUMMARY|assert_failures=" + str(len(FAILS)) + "|failed=" + (",".join(FAILS) if FAILS else "none") +
                "|compare_rows=" + str(n_cmp) + "|compare_disagree=" + str(disagree) + "|status=" + status +
                "|json=0\n")
    digest = sha256_file(out)
    with open(os.path.join(HERE, OUT_NAME + ".sha256"), "w", encoding="utf-8", newline="\n") as f:
        f.write(digest + "  " + OUT_NAME + "\n")
    print("MEASURE_BRIDGES|out=" + OUT_NAME + "|sha256=" + digest + "|measured_section_sha256=" + sec +
          "|assert_failures=" + str(len(FAILS)) + "|compare_disagree=" + str(disagree) + "|json=0")
    if FAILS:
        sys.exit(1)
    if disagree:
        sys.exit(3)
    sys.exit(0)


if __name__ == "__main__":
    main()
