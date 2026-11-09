# predict_ladder.py -- the PREDICTION, written before any measurement is run.
#
# What it predicts, per finite subgroup G of SU(2) acting on C^2:
#   a_d   = number of independent G-invariant polynomials of degree d = multiplicity of the
#           trivial rep in spin j with d = 2j (Molien series, from closed formulas only)
#   t     = largest degree up to which G looks exactly like continuous SU(2)
#   M_n   = E_G[chi_half^n] for n = 0..60, from a_d and Clebsch-Gordan counts
#   f_G   = log of the per-plaquette partition function of 2D Wilson gauge theory with free
#           boundary, exact rational power series through beta^30, and f_G - f_SU(2)
#
# HARD RULE kept here: no group element, matrix or quaternion is built, and nothing is summed
# over a group. Only the closed Molien formulas are used.
# Sources (NAMED): Klein 1884; Springer, Invariant Theory, LNM 585, 1977; the invariant degrees
# of the binary polyhedral groups.
#
# Arithmetic: integers only. A rational is a reduced pair (num, den) with den > 0, reduced by
# gcd after an exact divisibility assert, using only the integer operators // and %.
# Run: python3 -I -B predict_ladder.py   (writes PREDICTION.hbp and PREDICTION.hbp.sha256 next
# to this script; exits 1 if any self-check fails)
import hashlib
import os
import sys
from math import gcd

D_MAX = 120
N_MOM = 60
B_MAX = 30
DATE = "2026-11-09"
AUTHOR = "the local Opus 5.5 seat, pid 8467a937cba309f7 (prediction agent)"
SCHEMA = "ILLUSION_LADDER_PREDICTION_V1"

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


def rover(x, k):
    # x times the rational 1:k, k a nonzero integer
    return rat(x[0], x[1] * k)


def rs(x):
    return str(x[0]) + ":" + str(x[1])


def ints(v):
    return ",".join(str(e) for e in v)


def rats(v):
    return ",".join(rs(e) for e in v)


# ---------------------------------------------------------------- the groups (formulas only)


def group_list():
    out = []
    for n in range(2, 31):
        out.append({
            "name": "C_" + str(n), "alias": "", "family": "cyclic", "n": n, "order": n,
            "p": n, "q": 2, "r": n, "minus_one": n % 2 == 0,
            "num": "1+t^" + str(n), "den": "(1-t^2)(1-t^" + str(n) + ")",
        })
    for n in range(2, 31):
        out.append({
            "name": "Dic_" + str(n), "alias": "Q8" if n == 2 else "", "family": "binary_dihedral", "n": n,
            "order": 4 * n, "p": 2 * n + 2, "q": 4, "r": 2 * n, "minus_one": True,
            "num": "1+t^" + str(2 * n + 2), "den": "(1-t^4)(1-t^" + str(2 * n) + ")",
        })
    for name, order, p, q, r in (("2T", 24, 12, 6, 8), ("2O", 48, 18, 8, 12), ("2I", 120, 30, 12, 20)):
        out.append({
            "name": name, "alias": "", "family": "binary_polyhedral", "n": 0, "order": order,
            "p": p, "q": q, "r": r, "minus_one": True,
            "num": "1+t^" + str(p), "den": "(1-t^" + str(q) + ")(1-t^" + str(r) + ")",
        })
    return out


def expand(p, q, r, top):
    # (1 + t^p) times the geometric series in t^q times the geometric series in t^r,
    # coefficients of t^0..t^top. Multiplying by the geometric series in t^s is the running
    # recursion b_d = a_d + b_(d-s).
    a = [0] * (top + 1)
    a[0] += 1
    if p <= top:
        a[p] += 1
    for s in (q, r):
        for d in range(s, top + 1):
            a[d] += a[d - s]
    return a


def enumerate_route(p, q, r, top):
    # second route: count exponent pairs (i, j) with q*i + r*j = d, plus those with = d - p
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


def cyclic_monomial_route(n, top):
    # third route for C_n only: monomials x^u y^(d-u) with u = d-u mod n
    return [sum(1 for u in range(d + 1) if (u - (d - u)) % n == 0) for d in range(top + 1)]


# ---------------------------------------------------------------- SU(2) combinatorics

BINOM = [[1]]
for _n in range(1, N_MOM + 1):
    _prev = BINOM[-1]
    BINOM.append([1] + [_prev[k - 1] + _prev[k] for k in range(1, _n)] + [1])

FACT = [1]
for _n in range(1, N_MOM + 1):
    FACT.append(FACT[-1] * _n)


def binom(n, k):
    if k < 0 or k > n:
        return 0
    return BINOM[n][k]


def bmult(n, d):
    # multiplicity of spin j = d:2 in the n-fold tensor power of spin 1:2
    if d > n or (n - d) % 2 != 0:
        return 0
    k = (n - d) // 2
    return binom(n, k) - binom(n, k - 1)


def moments(a):
    return [sum(bmult(n, d) * a[d] for d in range(n + 1)) for n in range(N_MOM + 1)]


def z_series(m):
    return [rat(m[n], FACT[n] * 2 ** n) for n in range(B_MAX + 1)]


def smul(x, y):
    out = [ZERO] * (B_MAX + 1)
    for i in range(B_MAX + 1):
        if x[i][0] == 0:
            continue
        for j in range(B_MAX + 1 - i):
            if y[j][0] == 0:
                continue
            out[i + j] = radd(out[i + j], rmul(x[i], y[j]))
    return out


def log_series(z):
    # log(1 + x) with x = z - 1, the series sum over k of (-1)^(k+1) x^k times 1:k
    assert z[0] == ONE
    x = [ZERO] + z[1:]
    f = [ZERO] * (B_MAX + 1)
    power = x[:]
    k = 1
    while k <= B_MAX and any(c[0] != 0 for c in power):
        sign = 1 if k % 2 == 1 else -1
        for i in range(B_MAX + 1):
            if power[i][0] != 0:
                f[i] = radd(f[i], rover((sign * power[i][0], power[i][1]), k))
        power = smul(power, x)
        k += 1
    return f


def log_derivative_route(z):
    # second route: z' = f' z, so n f_n = n z_n - sum_(k=1..n-1) k f_k z_(n-k)
    assert z[0] == ONE
    f = [ZERO] * (B_MAX + 1)
    for n in range(1, B_MAX + 1):
        s = rmul((n, 1), z[n])
        for k in range(1, n):
            s = rsub(s, rmul(rmul((k, 1), f[k]), z[n - k]))
        f[n] = rover(s, n)
    return f


def first_positive_degree(a):
    for d in range(1, len(a)):
        if a[d] > 0:
            return d
    return -1


# ---------------------------------------------------------------- run


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    checks = []

    def check(name, ok, detail, forced=True):
        checks.append((name, bool(ok), detail, forced))

    groups = group_list()

    # Haar SU(2): the same formulas with a_0 = 1 and a_d = 0 for d > 0
    a_haar = [1] + [0] * D_MAX
    m_haar = moments(a_haar)
    z_haar = z_series(m_haar)
    f_haar = log_series(z_haar)
    f_haar_2 = log_derivative_route(z_haar)

    rows_g = []
    results = {}
    all_a0 = True
    two_route_bad = []
    cyc_bad = []
    order_bad = []
    gor_bad = []
    parity_bad = []
    tmin_bad = []
    mono_bad = []
    mbreak_bad = []
    log_bad = []
    gbreak_bad = []

    for g in groups:
        name = g["name"]
        a = expand(g["p"], g["q"], g["r"], D_MAX)
        if a[0] != 1:
            all_a0 = False
        if a != enumerate_route(g["p"], g["q"], g["r"], D_MAX):
            two_route_bad.append(name)
        if g["family"] == "cyclic" and a != cyclic_monomial_route(g["n"], D_MAX):
            cyc_bad.append(name)
        if g["q"] * g["r"] != 2 * g["order"]:
            order_bad.append(name)
        if g["p"] != g["q"] + g["r"] - 2:
            gor_bad.append(name)
        odd_nonzero = [d for d in range(1, D_MAX + 1, 2) if a[d] != 0]
        if g["minus_one"] and odd_nonzero:
            parity_bad.append(name)

        d0 = first_positive_degree(a)
        t_pred = d0 - 1
        m_break = a[d0]
        if t_pred != min(g["p"], g["q"], g["r"]) - 1:
            tmin_bad.append(name)

        m = moments(a)
        if any(m[n] < m_haar[n] for n in range(N_MOM + 1)):
            mono_bad.append(name)
        n_first = next((n for n in range(N_MOM + 1) if m[n] != m_haar[n]), -1)
        m_diff = m[n_first] - m_haar[n_first] if n_first >= 0 else 0
        if not (n_first == d0 and m_diff == m_break):
            mbreak_bad.append(name)

        z = z_series(m)
        f = log_series(z)
        if f != log_derivative_route(z):
            log_bad.append(name)
        diff = [rsub(f[i], f_haar[i]) for i in range(B_MAX + 1)]
        g_first = next((i for i in range(B_MAX + 1) if diff[i][0] != 0), -1)
        predicted_val = rat(m_break, (1 << d0) * FACT[d0])
        g_value = diff[g_first] if g_first >= 0 else ZERO
        zero_through = all(diff[i][0] == 0 for i in range(t_pred + 1))
        g_ok = zero_through and g_first == d0 and g_value == predicted_val
        if not g_ok:
            gbreak_bad.append(name)

        results[name] = {"t": t_pred, "d0": d0, "m": m_break, "gval": g_value, "gfirst": g_first, "pred": predicted_val, "ok": g_ok}

        label = ["group=" + name]
        if g["alias"]:
            label.append("alias=" + g["alias"])
        base = label + ["family=" + g["family"], "order=" + str(g["order"])]
        rows_g.append("|".join(["MOLIEN"] + base + [
            "contains_minus_one=" + ("1" if g["minus_one"] else "0"),
            "H=num(" + g["num"] + ")_over_den(" + g["den"] + ")",
            "d_range=0..120", "a=" + ints(a), "json=0"]))
        rows_g.append("|".join(["DESIGN"] + base + [
            "t_pred=" + str(t_pred), "first_break_degree=" + str(d0), "m_at_break=" + str(m_break),
            "rule=t_pred=min(" + str(g["p"]) + "," + str(g["q"]) + "," + str(g["r"]) + ")-1",
            "meaning=G_matches_SU2_on_every_spin_j_with_2j<=t_pred", "tag=PREDICTION", "json=0"]))
        rows_g.append("|".join(["MOMENT"] + base + [
            "quantity=E_G[chi_half^n]", "n_range=0..60", "M=" + ints(m), "json=0"]))
        rows_g.append("|".join(["MOMENT_BREAK"] + base + [
            "first_n_with_M_ne_Mh=" + str(n_first), "M=" + str(m[n_first]), "Mh=" + str(m_haar[n_first]),
            "diff=" + str(m_diff), "predicted_diff=a_at_first_break_degree=" + str(m_break),
            "agree=" + ("1" if n_first == d0 and m_diff == m_break else "0"), "tag=PREDICTION", "json=0"]))
        rows_g.append("|".join(["GAUGE_F"] + base + [
            "quantity=log_z(beta)", "beta_range=0..30", "rational_format=num:den", "f=" + rats(f), "json=0"]))
        rows_g.append("|".join(["GAUGE_DIFF"] + base + [
            "quantity=f_G-f_SU2", "beta_range=0..30", "rational_format=num:den", "diff=" + rats(diff), "json=0"]))
        rows_g.append("|".join(["GAUGE_FIRST_BREAK"] + base + [
            "zero_through_beta_power=" + str(t_pred) + ("" if zero_through else "_FAILED"),
            "first_nonzero_power=" + str(g_first), "value=" + rs(g_value),
            "predicted=a_" + str(d0) + ":(2^" + str(d0) + "*" + str(d0) + "!)=" + rs(predicted_val),
            "agree=" + ("1" if g_ok else "0"), "tag=PREDICTION", "json=0"]))

    # ---------------- the deliberately wrong control: 2I with numerator 1+t^29
    a_ctl = expand(29, 12, 20, D_MAX)
    ctl_odd = [d for d in range(1, D_MAX + 1, 2) if a_ctl[d] != 0]
    ctl_parity_flag = len(ctl_odd) > 0
    ctl_gor_flag = 29 != 12 + 20 - 2
    row_ctl = "|".join([
        "CONTROL", "name=2I_wrong_numerator", "H=num(1+t^29)_over_den((1-t^12)(1-t^20))", "contains_minus_one=1",
        "odd_degrees_with_nonzero_a=" + ints(ctl_odd), "a_29=" + str(a_ctl[29]),
        "parity_rule=" + ("VIOLATED" if ctl_parity_flag else "not_violated"),
        "gorenstein_degree_rule=" + ("VIOLATED_29_ne_30" if ctl_gor_flag else "not_violated"),
        "flagged=" + ("1" if ctl_parity_flag else "0"), "json=0"])

    # ---------------- Haar rows
    rows_h = [
        "|".join(["MOMENT_HAAR", "group=SU2_Haar", "rule=Mh_n=Catalan(n:2)_for_even_n_0_for_odd_n",
                  "n_range=0..60", "Mh=" + ints(m_haar), "json=0"]),
        "|".join(["GAUGE_HAAR", "group=SU2_Haar", "quantity=log_z(beta)", "beta_range=0..30",
                  "rational_format=num:den", "f=" + rats(f_haar), "json=0"]),
    ]

    # ---------------- ceiling
    t_all = {g["name"]: results[g["name"]]["t"] for g in groups}
    t_max = max(t_all.values())
    attain = [nm for nm in t_all if t_all[nm] == t_max]
    t_next = max(v for v in t_all.values() if v < t_max)
    next_by = [nm for nm in t_all if t_all[nm] == t_next]
    cyc_t = sorted(set(t_all[g["name"]] for g in groups if g["family"] == "cyclic"))
    dic_t = sorted(set(t_all[g["name"]] for g in groups if g["family"] == "binary_dihedral"))
    rows_c = [
        "|".join(["CEILING", "max_t_pred=" + str(t_max), "attained_by=" + ",".join(attain),
                  "count_attaining=" + str(len(attain)), "next_highest=" + str(t_next) + "(" + ",".join(next_by) + ")",
                  "groups_listed=" + str(len(groups)), "first_break_at_ceiling=degree_" + str(t_max + 1),
                  "rule=t_pred=min_positive_invariant_degree-1",
                  "scope=finite_subgroups_of_SU2_are_conjugate_to_C_n,Dic_n,2T,2O,2I(Klein_1884,NAMED)",
                  "tag=PREDICTION", "json=0"]),
        "|".join(["CEILING_FAMILY", "C_n_t_pred_values_n2..30=" + ints(cyc_t),
                  "Dic_n_t_pred_values_n2..30=" + ints(dic_t), "Dic_2=Q8_t_pred=" + str(t_all["Dic_2"]),
                  "statement=every_C_n(n>=2)_has_t_pred=1_and_every_Dic_n(n>=2)_has_t_pred=3_however_large_n",
                  "because=C_n_degrees(2,n,n)_min_2_and_Dic_n_degrees(4,2n,2n+2)_min_4_for_n>=2_and_all_coefficients_nonnegative",
                  "growing_order_does_not_raise_t=1", "tag=PREDICTION", "json=0"]),
    ]

    # ---------------- self-checks (each can fail)
    a2i = expand(30, 12, 20, D_MAX)
    aq8 = expand(6, 4, 4, D_MAX)
    ac4 = expand(4, 2, 4, D_MAX)
    ac3 = expand(3, 2, 3, D_MAX)
    check("a0_equals_1_every_group", all_a0, "groups=" + str(len(groups)))
    check("expansion_two_routes_agree", not two_route_bad, "geometric_recursion_vs_exponent_count_d0..120|bad=" + (",".join(two_route_bad) or "none"))
    check("cyclic_formula_vs_monomial_count", not cyc_bad, "C_2..C_30_d0..120|bad=" + (",".join(cyc_bad) or "none"))
    check("hand_2I_a12_eq_1", a2i[12] == 1, "computed=" + str(a2i[12]))
    check("hand_2I_a20_eq_1", a2i[20] == 1, "computed=" + str(a2i[20]))
    check("hand_2I_a24_eq_1", a2i[24] == 1, "computed=" + str(a2i[24]))
    from_one = sum(1 for i in range(3) for j in range(2) if 12 * i + 20 * j == 30)
    from_t30 = sum(1 for i in range(3) for j in range(2) if 12 * i + 20 * j == 30 - 30)
    check("report_2I_a30", a2i[30] == from_one + from_t30, "computed=" + str(a2i[30]) + "|from_numerator_1=" + str(from_one)
          + "|from_numerator_t^30=" + str(from_t30) + "|note=reported_not_asserted", forced=False)
    check("hand_Q8_a4_eq_2", aq8[4] == 2, "computed=" + str(aq8[4]))
    check("hand_C4_a4_eq_3", ac4[4] == 3, "computed=" + str(ac4[4]) + "|monomials=x^4,y^4,x^2y^2")
    check("hand_C3_a3_eq_2", ac3[3] == 2, "computed=" + str(ac3[3]) + "|monomials=x^3,y^3")
    check("order_from_degrees_q*r_eq_2*order", not order_bad, "bad=" + (",".join(order_bad) or "none"))
    check("gorenstein_degree_p_eq_q+r-2", not gor_bad, "bad=" + (",".join(gor_bad) or "none"))
    check("parity_minus_one_odd_a_zero", not parity_bad, "groups_with_minus_one_checked|bad=" + (",".join(parity_bad) or "none"))
    check("parity_control_2I_numerator_t^29_flagged", ctl_parity_flag, "odd_nonzero_first=" + (str(ctl_odd[0]) if ctl_odd else "none"))
    check("t_pred_eq_min_degree_minus_1", not tmin_bad, "bad=" + (",".join(tmin_bad) or "none"))

    cg_ok = True
    dim_ok = True
    mult = [1]
    for n in range(1, N_MOM + 1):
        new = [0] * (n + 1)
        for d in range(n + 1):
            if d >= 1 and d - 1 <= n - 1:
                new[d] += mult[d - 1]
            if d + 1 <= n - 1:
                new[d] += mult[d + 1]
        mult = new
        if any(mult[d] != bmult(n, d) for d in range(n + 1)):
            cg_ok = False
        if sum(bmult(n, d) * (d + 1) for d in range(n + 1)) != (1 << n):
            dim_ok = False
    check("clebsch_gordan_closed_form_vs_recursion", cg_ok, "n=1..60")
    check("dimension_sum_eq_2^n", dim_ok, "n=0..60")

    cat = [1]
    for k in range(30):
        cat.append(sum(cat[i] * cat[k - i] for i in range(k + 1)))
    haar_ok = all(m_haar[n] == (cat[n // 2] if n % 2 == 0 else 0) for n in range(N_MOM + 1))
    check("haar_moments_eq_catalan_recurrence", haar_ok, "n=0..60|Mh_60=" + str(m_haar[60]))
    check("moments_at_least_haar", not mono_bad, "bad=" + (",".join(mono_bad) or "none"))
    check("moment_break_eq_design_break", not mbreak_bad, "bad=" + (",".join(mbreak_bad) or "none"))
    check("log_two_routes_agree", f_haar == f_haar_2 and not log_bad,
          "log_series_vs_derivative_recursion_beta0..30|bad=" + (",".join(log_bad) or "none"))
    check("gauge_break_eq_design_break", not gbreak_bad, "bad=" + (",".join(gbreak_bad) or "none"))
    check("ceiling_eq_11_by_2I_only", t_max == 11 and attain == ["2I"], "max=" + str(t_max) + "|by=" + ",".join(attain))
    check("family_C_n_t1_Dic_n_t3", cyc_t == [1] and dic_t == [3], "C=" + ints(cyc_t) + "|Dic=" + ints(dic_t))

    failed = sum(1 for c in checks if c[3] and not c[1])
    rows_s = ["|".join(["SELF_CHECK", "name=" + nm, "result=" + ("PASS" if ok else "FAIL"),
                        "forced=" + ("1" if forced else "0"), "detail=" + det, "json=0"])
              for nm, ok, det, forced in checks]

    head = "|".join([
        "ILLUSION_PREDICT_HEAD", "schema=" + SCHEMA, "date=" + DATE, "author=" + AUTHOR,
        "formula_sources=Klein_1884;Springer_Invariant_Theory_LNM_585_1977;binary_polyhedral_invariant_degrees",
        "formula_tag=NAMED", "measurement_status=NOT_RUN", "written_before_any_measurement=1",
        "group_elements_built=0", "sums_over_a_group=0", "integer_only=1", "rational_format=num:den_reduced_den_positive",
        "groups=C_2..C_30,Dic_2..Dic_30(Dic_2=Q8),2T,2O,2I", "group_count=" + str(len(groups)),
        "d_range=0..120", "moment_range=0..60", "beta_range=0..30",
        "gauge_model=2D_Wilson_free_boundary_independent_plaquettes_weight_exp(beta*chi_half(U_p)*(1:2))_U_p_uniform_on_G",
        "z(beta)=sum_n_beta^n*M_n*(1:(n!*2^n))", "f=log_z",
        "script=predict_ladder.py", "script_sha256=" + script_sha, "json=0"])

    body = [head] + rows_g + rows_h + [row_ctl] + rows_c + rows_s
    status = "ALL_SELF_CHECKS_PASS" if failed == 0 else "SELF_CHECK_FAILED"
    end = "|".join(["ILLUSION_PREDICT_END", "rows=" + str(len(body) + 1), "self_checks=" + str(len(checks)),
                    "failed=" + str(failed), "status=" + status, "measurement_status=NOT_RUN", "json=0"])
    text = "\n".join(body + [end]) + "\n"
    data = text.encode("ascii")
    out = os.path.join(here, "PREDICTION.hbp")
    with open(out, "wb") as fh:
        fh.write(data)
    digest = hashlib.sha256(data).hexdigest()
    with open(out + ".sha256", "wb") as fh:
        fh.write((digest + "  PREDICTION.hbp\n").encode("ascii"))

    for nm in ("Dic_2", "Dic_3", "2T", "2O", "2I", "C_3", "C_4"):
        r = results[nm]
        print("FIRST_BREAK|group=" + nm + "|t_pred=" + str(r["t"]) + "|beta_power=" + str(r["gfirst"])
              + "|value=" + rs(r["gval"]) + "|agree=" + ("1" if r["ok"] else "0") + "|json=0")
    print("PREDICT_DONE|rows=" + str(len(body) + 1) + "|self_checks=" + str(len(checks)) + "|failed=" + str(failed)
          + "|script_sha256=" + script_sha + "|hbp_sha256=" + digest + "|json=0")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
