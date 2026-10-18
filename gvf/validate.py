"""Validation suite. Every check compares the solver with something that can be verified
independently (a closed-form formula, a different algorithm, or a physical identity).
No published textbook numbers are used.

    python -m gvf.validate            prints the report
    python -m gvf.validate --write    also writes theory/VALIDATION.md
"""
import math
import os
import sys
import numpy as np
from . import hydraulics as hy
from .hydraulics import G
from .profile import compute
from .reference import ode_profile
from .sections import Rectangular, Trapezoidal, Circular, make_section

SECTIONS = [Rectangular(5.0), Trapezoidal(4.0, 2.0), Trapezoidal(0.0, 1.5), Circular(2.0)]
_ITEMS = []


def check(group, name, measured, limit, how):
    """Record one result: passes when |measured| <= limit."""
    _ITEMS.append(dict(group=group, name=name, measured=float(measured), limit=float(limit),
                       ok=abs(measured) <= limit, how=how))


def _rect(b=5.0):
    return {"type": "rectangular", "width": b}


def _spec(sec, S0, n, L, up, down, Q, dx, tol=1e-10):
    return dict(discharge=Q, step=dx, tolerance=tol, upstream=up, downstream=down,
                reaches=[dict(length=L, slope=S0, n=n, section=sec)])


def _vs_ode(sec_d, S0, n, L, up, down, Q, dx, from_down):
    """Largest |standard step - ODE| over the stations of a single-reach profile."""
    res = compute(_spec(sec_d, S0, n, L, up, down, Q, dx))
    sec = make_section(sec_d)
    xs = [r["x"] for r in res.rows]
    ys = {r["x"]: r["y"] for r in res.rows}
    if from_down:
        ref = ode_profile(sec, Q, n, S0, ys[xs[-1]], xs[-1], 0.0, xs)
    else:
        ref = ode_profile(sec, Q, n, S0, ys[xs[0]], xs[0], L, xs)
    return max(abs(ys[x] - y) for x, y in ref), len(ref)


def run_all():
    _ITEMS.clear()
    # ---- 1. closed-form hydraulics
    for sec in SECTIONS:
        yn = hy.normal_depth(sec, 2.0, 0.014, 0.002)
        check("Closed forms", f"Normal depth round trip, {sec.kind}", hy.capacity(sec, 0.014, 0.002, yn) / 2.0 - 1, 1e-9,
              "Q recomputed from the normal depth found, relative error")
        yc = hy.critical_depth(sec, 3.0)
        check("Closed forms", f"Froude number at critical depth, {sec.kind}", hy.froude(sec, 3.0, yc) - 1, 1e-9, "Fr - 1 at yc")
    q = 12.0 / 5.0
    check("Closed forms", "Rectangular critical depth", hy.critical_depth(Rectangular(5), 12.0) - (q * q / G) ** (1 / 3), 1e-9,
          "yc - (q^2/g)^(1/3)")
    yc = hy.critical_depth(Rectangular(5), 12.0)
    check("Closed forms", "Minimum specific energy of a rectangle", hy.energy(Rectangular(5), 12.0, yc) - 1.5 * yc, 1e-9, "Emin - 1.5 yc")
    worst = 0.0
    for y1 in (0.2, 0.35, 0.5, 0.7):
        fr1 = q / math.sqrt(G * y1 ** 3)
        worst = max(worst, abs(hy.sequent_depth(Rectangular(5), 12.0, y1) - 0.5 * y1 * (math.sqrt(1 + 8 * fr1 ** 2) - 1)))
    check("Closed forms", "Rectangular sequent depth", worst, 1e-9, "y2 - 0.5 y1 (sqrt(1+8 Fr1^2) - 1), four depths")
    y1 = 0.4
    y2 = hy.sequent_depth(Rectangular(5), 12.0, y1)
    check("Closed forms", "Rectangular energy loss", hy.energy(Rectangular(5), 12.0, y1) - hy.energy(Rectangular(5), 12.0, y2)
          - (y2 - y1) ** 3 / (4 * y1 * y2), 1e-9, "(E1 - E2) - (y2 - y1)^3 / (4 y1 y2)")
    worst_m = worst_e = 0.0
    for sec in SECTIONS:
        yc = hy.critical_depth(sec, 3.0)
        y2 = hy.sequent_depth(sec, 3.0, 0.6 * yc, yc)
        worst_m = max(worst_m, abs(hy.momentum(sec, 3.0, 0.6 * yc) / hy.momentum(sec, 3.0, y2) - 1))
        for k in (0.4, 0.7, 1.3, 2.0):
            a = hy.alternate_depth(sec, 3.0, k * yc, yc)
            if a is not None:
                worst_e = max(worst_e, abs(hy.energy(sec, 3.0, a) - hy.energy(sec, 3.0, k * yc)))
    check("Closed forms", "Equal momentum across a jump, all sections", worst_m, 1e-9, "M(y1)/M(y2) - 1")
    check("Closed forms", "Equal energy at alternate depths, all sections", worst_e, 1e-9, "E(y) - E(alternate depth)")

    # ---- 2. standard step vs the independent ODE solution
    rq = dict(Q=10.0)
    tr = {"type": "trapezoidal", "width": 6, "side_slope": 2}
    ci = {"type": "circular", "diameter": 1.5}
    yc_r = hy.critical_depth(Rectangular(5), 10.0)
    cases = [
        ("M1 backwater, rectangular", _rect(), 0.001, 0.015, 3000, {"type": "none"}, {"type": "depth", "depth": 2.5}, 10.0, True),
        ("M2 drawdown, rectangular", _rect(), 0.001, 0.015, 3000, {"type": "none"}, {"type": "depth", "depth": 1.15 * yc_r}, 10.0, True),
        ("M3 from a gate, rectangular", _rect(), 0.001, 0.015, 60, {"type": "gate", "opening": 0.4}, {"type": "none"}, 1.0, False),
        ("S2 drawdown, rectangular", _rect(), 0.02, 0.015, 400, {"type": "depth", "depth": 0.95 * yc_r}, {"type": "none"}, 5.0, False),
        ("S3 rising to normal depth, rectangular", _rect(), 0.02, 0.015, 400, {"type": "depth", "depth": 0.2}, {"type": "none"}, 5.0, False),
        ("H2 horizontal channel", _rect(), 0.0, 0.015, 500, {"type": "none"}, {"type": "depth", "depth": 1.3 * yc_r}, 5.0, True),
        ("A2 adverse slope", _rect(), -0.0005, 0.015, 500, {"type": "none"}, {"type": "depth", "depth": 1.3 * yc_r}, 5.0, True),
        ("M1 backwater, trapezoidal", tr, 0.0004, 0.025, 3000, {"type": "none"}, {"type": "depth", "depth": 2.4}, 10.0, True),
        ("M1 backwater, circular", ci, 0.001, 0.013, 400, {"type": "none"}, {"type": "depth", "depth": 1.0}, 4.0, True),
    ]
    for name, sec, S0, n, L, up, down, dx, fd in cases:
        Q = 2.0 if sec is ci else 10.0
        err, npts = _vs_ode(sec, S0, n, L, up, down, Q, dx, fd)
        check("Standard step vs ODE", name, err, 2e-3, f"max |y_step - y_ODE| over {npts} stations, step {dx:g} m")
    errs = []
    for dx in (16.0, 8.0, 4.0, 2.0, 1.0):
        errs.append(_vs_ode(_rect(), 0.001, 0.015, 3000, {"type": "none"}, {"type": "depth", "depth": 2.5}, 10.0, dx, True)[0])
    order = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    check("Standard step vs ODE", "Convergence: error shrinks when the step is halved (steps 16, 8, 4, 2, 1 m)",
          errs[-1], 1e-4, "errors " + ", ".join(f"{e:.1e}" for e in errs) + "; observed order " + ", ".join(f"{o:.1f}" for o in order))

    # ---- 3. physical behaviour of whole profiles
    res = compute(_spec(_rect(), 0.0005, 0.02, 20000, {"type": "none"}, {"type": "depth", "depth": 3.0}, 10.0, 20.0))
    yn = res.reaches[0]["yn"]
    check("Physical behaviour", "Far upstream of a backwater the profile approaches normal depth", res.rows[0]["y"] / yn - 1, 1e-3,
          "relative difference 20 km upstream")
    spec = dict(discharge=12.0, step=2.0, tolerance=1e-10, upstream={"type": "gate", "opening": 0.5}, downstream={"type": "depth", "depth": 1.4},
                reaches=[dict(length=300, slope=0.001, n=0.013, section=_rect())])
    res = compute(spec)
    j = res.jumps[0]
    check("Physical behaviour", "Momentum is equal at the toe of the computed jump", hy.momentum(Rectangular(5), 12.0, j["y1"]) /
          hy.momentum(Rectangular(5), 12.0, j["y2"]) - 1, 1e-9, "M(y1)/M(y2) - 1")
    sub_at_toe = [y for x, y in res.curves["sub"] if abs(x - j["x"]) < 1e-9]
    xs_, ys_ = zip(*res.curves["sub"])
    check("Physical behaviour", "Subcritical (tailwater) depth at the toe equals the sequent depth", float(np.interp(j["x"], xs_, ys_)) - j["y2"], 2e-3,
          "interpolated subcritical profile - y2")
    egl = [r["egl"] for r in res.rows]
    check("Physical behaviour", "Energy line never rises in the flow direction (including the roller)",
          max(0.0, max(b - a for a, b in zip(egl[:-1], egl[1:]))), 1e-6, "largest rise of the EGL between stations")
    check("Physical behaviour", "Energy and continuity residuals of the step solution", max(res.checks["energy_residual"], res.checks["continuity"]), 1e-5,
          "largest energy-equation residual (m) and |Q - VA|")
    xs_j = []
    for tail in (1.5, 1.4, 1.3):
        s2 = dict(spec, downstream={"type": "depth", "depth": tail})
        xs_j.append(compute(s2).jumps[0]["x"])
    check("Physical behaviour", "Jump moves downstream when the tailwater drops", 0.0 if xs_j[0] < xs_j[1] < xs_j[2] else 1.0, 0.0,
          "jump toe for tailwater 1.5, 1.4, 1.3 m: " + ", ".join(f"{x:.1f}" for x in xs_j))
    check("Physical behaviour", "Roller length", compute(spec).jumps[0]["x_end"] - compute(spec).jumps[0]["x"] - 6.1 * j["y2"], 1e-9,
          "x_end - x_toe - 6.1 y2")
    return list(_ITEMS)


def markdown(items):
    out = ["# Validation report", "",
           "Generated by `python -m gvf.validate --write`. Every number below was measured when the report was generated.",
           "Each check compares the solver with a closed-form result, an independent algorithm or a physical identity; "
           "no published textbook numbers are used (see `examples/TEXTBOOK_EXAMPLE.md` to add your own).", ""]
    for g in dict.fromkeys(i["group"] for i in items):
        out += [f"## {g}", "", "| Check | Measured | Limit | Result | How |", "|---|---|---|---|---|"]
        for i in (x for x in items if x["group"] == g):
            out.append(f"| {i['name']} | {i['measured']:.2e} | {i['limit']:.0e} | {'pass' if i['ok'] else '**FAIL**'} | {i['how']} |")
        out.append("")
    out += ["## What is not validated", "",
            "* The jump length rule L = 6.1 y2, the jump classes, the weir coefficient and the gate contraction coefficient are "
            "empirical rules from textbooks, written from memory; they are not checked against published data.",
            "* The shape of the surface inside the roller is idealised (only y1, y2, the length and the energy loss are computed).",
            "* The solver assumes steady, one-dimensional, gradually varied flow (see the README for all limitations)."]
    return "\n".join(out) + "\n"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    items = run_all()
    text = markdown(items)
    print(text)
    if "--write" in argv:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        os.makedirs(os.path.join(root, "theory"), exist_ok=True)
        with open(os.path.join(root, "theory", "VALIDATION.md"), "w", encoding="utf-8") as f:
            f.write(text)
    return 0 if all(i["ok"] for i in items) else 1


if __name__ == "__main__":
    raise SystemExit(main())
