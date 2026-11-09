#!/usr/bin/env python3
# Independent integer-only audit of the M1/M2/M3 scripts (refuter seat, E=0).
# (a) static: tokenize each script and count banned tokens.
# (b) dynamic: run each script copy under sys.settrace and flag any float/complex
#     value returned by a function or held in a function's locals at return time.
# Output rows go to stdout as TAG|k=v|...|json=0.
import sys
import io
import os
import tokenize
import runpy
import hashlib

BASE = "/mnt/c/tmp/opus55/g2-s6-20261109"
AUD = BASE + "/refute/audit"
JOBS = [("m1", "m1_octonion.py"), ("m2", "m2_roots.py"), ("m3", "m3_towers.py")]
BANNED_NAMES = {"float", "complex", "Fraction", "fractions", "Decimal", "decimal",
                "isclose", "fsum", "truediv",
                "numpy", "np", "sympy", "mpmath", "cmath", "statistics", "random"}
ALLOWED_MATH = {"gcd", "isqrt"}


def static_scan(src):
    bad = []
    toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    for i, t in enumerate(toks):
        if t.type == tokenize.OP and t.string in ("/", "//", "/=", "//="):
            # '//' is floor division (integer); flag it separately, it is integer-safe on ints
            bad.append(("op" + t.string, t.start[0]))
        elif t.type == tokenize.NUMBER:
            s = t.string.lower()
            if s.startswith(("0x", "0o", "0b")):
                continue
            if "." in s or "e" in s or s.endswith("j"):
                bad.append(("floatlit_" + t.string, t.start[0]))
        elif t.type == tokenize.NAME and t.string in BANNED_NAMES:
            bad.append(("name_" + t.string, t.start[0]))
        elif t.type == tokenize.NAME and t.string == "math":
            if i + 2 < len(toks) and toks[i + 1].string == ".":
                if toks[i + 2].string not in ALLOWED_MATH:
                    bad.append(("math_" + toks[i + 2].string, t.start[0]))
        elif t.type == tokenize.OP and t.string == "**":
            nxt = toks[i + 1].string if i + 1 < len(toks) else ""
            if nxt == "-":
                bad.append(("negpow", t.start[0]))
        elif t.type == tokenize.NAME and t.string == "pow":
            pass
    return bad


FOUND = []
TARGET = [""]


def _isbad(v):
    tv = type(v)
    return tv is float or tv is complex


def _scan_val(v, depth):
    if _isbad(v):
        return True
    if depth > 0 and type(v) in (list, tuple):
        for w in v:
            if _scan_val(w, depth - 1):
                return True
    if depth > 0 and type(v) is dict:
        for w in v.values():
            if _scan_val(w, depth - 1):
                return True
    return False


def local_tracer(frame, event, arg):
    if event == "return" and frame.f_code.co_filename == TARGET[0]:
        if _scan_val(arg, 2):
            FOUND.append(("ret", frame.f_code.co_name, frame.f_lineno))
        for k, v in frame.f_locals.items():
            if _scan_val(v, 2):
                FOUND.append(("local:" + k, frame.f_code.co_name, frame.f_lineno))
                break
    return local_tracer


def global_tracer(frame, event, arg):
    if event == "call":
        frame.f_trace_lines = False
        return local_tracer
    return None


def _ctl_file():
    return os.path.abspath(__file__)


def main():
    rows = []
    # control: the static scanner must catch a bad sample
    sample = "x = 1 " + chr(47) + " 3\ny = 2" + chr(46) + "5\nimport math\nz = math.sqrt(2)\nw = float(3)\n"
    cb = static_scan(sample)
    rows.append("AUDIT_CONTROL_STATIC|expected=FAIL|observed=%s|flagged=%d|status=%s|json=0" % (
        "FAIL" if cb else "PASS", len(cb), "PASS" if len(cb) >= 4 else "FAIL"))
    ok_all = len(cb) >= 4
    for d, s in JOBS:
        p = BASE + "/" + d + "/" + s
        src = open(p, encoding="utf-8").read()
        sha = hashlib.sha256(src.encode("utf-8")).hexdigest()
        bad = static_scan(src)
        slash_chars = src.count(chr(47))
        floordiv = [b for b in bad if b[0] in ("op//", "op//=")]
        truediv = [b for b in bad if b[0] in ("op/", "op/=")]
        other = [b for b in bad if not b[0].startswith("op")]
        rows.append("AUDIT_STATIC|job=%s|script_sha256=%s|truediv_ops=%d|floordiv_ops=%d|other_findings=%d|slash_chars_total=%d|detail=%s|status=%s|json=0" % (
            d, sha, len(truediv), len(floordiv), len(other), slash_chars,
            ";".join("%s@%d" % b for b in (truediv + other)[:12]) or "none",
            "PASS" if not truediv and not other else "FAIL"))
        if truediv or other:
            ok_all = False
    # dynamic trace
    sys.setrecursionlimit(10000)
    for d, s in JOBS:
        p = AUD + "/" + d + "/" + s
        del FOUND[:]
        TARGET[0] = p
        old_argv = sys.argv
        sys.argv = [p]
        cwd = os.getcwd()
        os.chdir(AUD + "/" + d)
        exit_code = 0
        devnull = open(os.devnull, "w")
        old_out = sys.stdout
        sys.stdout = devnull
        sys.settrace(global_tracer)
        try:
            g = runpy.run_path(p, run_name="__main__")
        except SystemExit as e:
            exit_code = e.code if isinstance(e.code, int) else 1
            g = {}
        finally:
            sys.settrace(None)
            sys.stdout = old_out
            devnull.close()
            os.chdir(cwd)
            sys.argv = old_argv
        gbad = [k for k, v in g.items() if _scan_val(v, 2)]
        n = len(FOUND)
        rows.append("AUDIT_RUNTIME_FLOAT|job=%s|exit=%s|float_or_complex_events=%d|module_globals_bad=%d|first=%s|status=%s|json=0" % (
            d, exit_code, n, len(gbad), (",".join("%s:%s:%d" % f for f in FOUND[:3]) or "none"),
            "PASS" if n == 0 and not gbad and exit_code in (0, None) else "FAIL"))
        if n or gbad or exit_code not in (0, None):
            ok_all = False
    # dynamic control: a function that makes a float must be caught
    del FOUND[:]

    TARGET[0] = _ctl_file()

    def _ctl():
        q = 7
        r = q / 2
        return r
    sys.settrace(global_tracer)
    _ctl()
    sys.settrace(None)
    rows.append("AUDIT_CONTROL_RUNTIME|expected=FAIL|observed=%s|events=%d|status=%s|json=0" % (
        "FAIL" if FOUND else "PASS", len(FOUND), "PASS" if FOUND else "FAIL"))
    if not FOUND:
        ok_all = False
    for r in rows:
        print(r)
    sys.exit(0 if ok_all else 1)


main()
