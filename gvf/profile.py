"""Water surface profile of a multi-reach channel, with hydraulic jumps.

Algorithm (see theory/THEORY.md for the derivations):
  1. build the reaches and a grid of stations
  2. find the controls: boundary conditions plus critical depth at mild->steep breaks
  3. sweep from every control with the standard step method: subcritical flow is computed
     UPSTREAM from its control, supercritical flow DOWNSTREAM from its control
  4. walk down the channel along the supercritical profile; where its momentum function falls
     to that of the subcritical profile, a hydraulic jump starts (its toe)
  5. draw the jump as a roller of length L = 6.1 y2, then follow the subcritical profile
"""
import math
from types import SimpleNamespace
import pandas as pd
from scipy.optimize import brentq
from . import hydraulics as hy
from .hydraulics import G
from .jump import jump_properties, smooth, LENGTH_FACTOR, LENGTH_NOTE, ROLLER_STATIONS
from .sections import make_section, ChannelError
from .step import Ctx, sweep, solve_step

UP_TYPES = ("normal", "depth", "level", "gate", "critical", "none")
DOWN_TYPES = ("normal", "depth", "level", "weir", "overfall", "none")
MAX_STATIONS = 20000
NAN = float("nan")

COLUMNS = {  # short internal name -> readable name used in tables and CSV
    "x": "Distance (m)", "reach": "Reach", "z": "Bed elevation (m)", "y": "Depth (m)",
    "wse": "Water surface elevation (m)", "v": "Velocity (m/s)", "fr": "Froude number",
    "e": "Specific energy (m)", "sf": "Friction slope", "egl": "Energy grade line elevation (m)",
    "ptype": "Profile type", "flow": "Flow", "note": "Note"}


# --------------------------------------------------------------------------- input
def _num(v, name, positive=False, zero_ok=False):
    try:
        v = float(v)
    except (TypeError, ValueError):
        raise ChannelError(f"{name} must be a number (you gave {v!r}).")
    if math.isnan(v) or math.isinf(v):
        raise ChannelError(f"{name} must be a finite number.")
    if positive and (v < 0 or (v == 0 and not zero_ok)):
        raise ChannelError(f"{name} must be greater than zero (you gave {v:g}).")
    return v


class Reach:
    """One prismatic reach: geometry, roughness, normal/critical depth and slope class."""


def _build_reaches(spec, Q, z_start):
    rl = spec.get("reaches")
    if not isinstance(rl, list) or not rl:
        raise ChannelError("add at least one reach to the channel.")
    out, x0, z0 = [], 0.0, z_start
    for k, rd in enumerate(rl):
        tag = f"Reach {k + 1}"
        r = Reach()
        r.idx, r.name = k, tag
        r.L = _num(rd.get("length"), f"{tag}: length", positive=True)
        r.S0 = _num(rd.get("slope"), f"{tag}: bed slope")
        r.n = _num(rd.get("n", 0.015), f"{tag}: Manning's n", positive=True)
        try:
            r.sec = make_section(rd.get("section"))
            r.yc = hy.critical_depth(r.sec, Q)
            r.yn = hy.normal_depth(r.sec, Q, r.n, r.S0)
        except ChannelError as e:
            raise ChannelError(f"{tag}: {e}")
        r.sc = hy.friction_slope(r.sec, r.n, Q, r.yc)   # critical slope for this Q
        if r.S0 <= 0:
            r.kind, r.cls = "sub", ("H" if r.S0 == 0 else "A")
        else:
            rel = (r.yn - r.yc) / r.yc
            r.kind, r.cls = ("sub", "M") if rel > 0.01 else (("sup", "S") if rel < -0.01 else ("crit", "C"))
        r.x0, r.z0 = x0, z0
        out.append(r)
        x0 += r.L
        z0 -= r.S0 * r.L
    return out


def _bc_depth(bc, side, reaches, Q):
    """Depth implied by a boundary condition: (depth or None, type)."""
    up = side == "upstream"
    r = reaches[0] if up else reaches[-1]
    t = str((bc or {}).get("type", "normal")).lower()
    bc = bc or {}
    if t not in (UP_TYPES if up else DOWN_TYPES):
        raise ChannelError(f"the {side} condition '{t}' is not available. "
                           f"Choose one of: {', '.join(UP_TYPES if up else DOWN_TYPES)}.")
    where = f"{side} condition"
    if t == "none":
        return None, t
    if t == "normal":
        if r.yn is None:
            raise ChannelError(f"the {side} reach has a {'horizontal' if r.S0 == 0 else 'adverse'} slope, so it has no "
                               f"normal depth. Choose another {side} condition (for example a depth).")
        return r.yn, t
    if t in ("critical", "overfall"):
        return r.yc, t
    if t == "depth":
        return _num(bc.get("depth"), f"{where}: depth", positive=True), t
    if t == "level":
        lvl = _num(bc.get("level"), f"{where}: water surface elevation")
        zb = r.z0 if up else r.z0 - r.S0 * r.L
        if lvl - zb <= 0:
            raise ChannelError(f"the {where} water level ({lvl:g} m) is not above the channel bed ({zb:.3f} m).")
        return lvl - zb, t
    if t == "gate":
        a = _num(bc.get("opening"), f"{where}: gate opening", positive=True)
        return _num(bc.get("cc", 0.61), f"{where}: contraction coefficient", positive=True) * a, t
    if t == "weir":
        if r.sec.kind != "rectangular":
            raise ChannelError("the weir option works only for a rectangular last reach. "
                               "For other sections, enter the tailwater depth directly.")
        crest = _num(bc.get("crest"), f"{where}: weir crest height", positive=True, zero_ok=True)
        cw = _num(bc.get("cw", 1.7), f"{where}: weir coefficient", positive=True)
        return crest + (Q / (cw * r.sec.b)) ** (2.0 / 3.0), t


def classify(r, y, flow="subcritical"):
    """Profile type from the depth relative to normal and critical depth.
    A depth exactly at critical depth counts as 'above' for subcritical flow and 'below' for supercritical."""
    yn, yc, c = r.yn, r.yc, r.cls
    if c in "MS" and abs(y - yn) <= 0.01 * yn:
        return "uniform"
    above = y >= yc * (1 - 1e-6) if flow == "subcritical" else y > yc * (1 + 1e-6)
    if c == "M":
        return "M1" if y > yn else ("M2" if above else "M3")
    if c == "S":
        return "S1" if above else ("S2" if y > yn else "S3")
    if c == "C":
        return "C1" if above else "C3"
    return f"{c}{2 if above else 3}"


# --------------------------------------------------------------------------- result
class Result:
    def __init__(self):
        self.warnings, self.jumps, self.reaches, self.rows = [], [], [], []
        self.name, self.Q, self.table, self.checks = "channel", 0.0, None, {}
        self.curves = {"sup": [], "sub": []}   # the two profiles before they are cut at the jump

    def summary(self):
        L = [f"{self.name}   (Q = {self.Q:g} m3/s)", ""]
        for r in self.reaches:
            yn = "none (no normal depth)" if r["yn"] is None else f"{r['yn']:.3f} m"
            L.append(f"{r['name']}: {r['slope_class']} slope | normal depth {yn} | critical depth {r['yc']:.3f} m"
                     f" | profiles: {', '.join(r['types']) or '-'}")
        L.append("")
        if not self.jumps:
            L.append("Hydraulic jump: none")
        for j in self.jumps:
            if j["status"] == "free":
                L.append(f"Hydraulic jump ({j['type']}, Fr1 = {j['fr1']:.2f}) in {j['reach_name']}:")
                L.append(f"  toe at x = {j['x']:.2f} m with y1 = {j['y1']:.3f} m; sequent depth y2 = {j['y2']:.3f} m")
                L.append(f"  roller length ~ {j['length']:.1f} m (approximate), so the jump ends near x = {j['x_end']:.2f} m "
                         f"where the depth is {j['y_end']:.3f} m")
                L.append(f"  energy loss {j['loss']:.3f} m ({j['loss_pct']:.1f}% of E1 = {j['e1']:.3f} m)")
            elif j["status"] == "drowned":
                L.append(f"Jump DROWNED: the tailwater ({j['y_tail']:.3f} m at the upstream end) is higher than the sequent "
                         f"depth ({j['y2']:.3f} m) needed for y1 = {j['y1']:.3f} m, so the jump is pushed up to the upstream end (gate).")
            else:
                L.append(f"Jump SWEPT OUT: flow stays supercritical to the end ({j['y1']:.3f} m). The tailwater "
                         f"({j['y_tail']:.3f} m) is lower than the sequent depth needed ({j['y2']:.3f} m).")
        if any(j["status"] == "free" for j in self.jumps):
            L.append(f"Note: {LENGTH_NOTE}")
        if self.checks:
            L += ["", "Self-checks: " + ", ".join(f"{k} = {v:.2e}" for k, v in self.checks.items())]
        if self.warnings:
            L += ["", "Warnings:"] + [f"  - {w}" for w in self.warnings]
        return "\n".join(L)


# --------------------------------------------------------------------------- main
def compute(spec):
    """Compute the water surface profile for a channel description (a dict, see README)."""
    if not isinstance(spec, dict):
        raise ChannelError("the channel description must be a JSON object.")
    m = SimpleNamespace(res=Result())
    _setup(m, spec)
    _controls(m, spec)
    _sweeps(m)
    _assemble(m)
    _table(m)
    _finish(m)
    return m.res


def _setup(m, spec):
    m.Q = Q = _num(spec.get("discharge"), "Discharge Q", positive=True)
    z_start = _num(spec.get("start_bed_elevation", 100.0), "Starting bed elevation")
    step = _num(spec.get("step", 1.0), "Step size", positive=True)
    m.tol = _num(spec.get("tolerance", 1e-6), "Tolerance", positive=True)
    m.res.name, m.res.Q = str(spec.get("name", "channel")), Q
    m.R = _build_reaches(spec, Q, z_start)
    counts = [max(20, math.ceil(r.L / step - 1e-9)) for r in m.R]       # at least 20 steps per reach
    total = sum(c + 1 for c in counts)
    if total > MAX_STATIONS:
        raise ChannelError(f"that would need {total} stations. Increase the step size (limit {MAX_STATIONS} stations).")
    m.xs, m.rs, m.zs, m.first, m.last = [], [], [], {}, {}
    for r, ns in zip(m.R, counts):   # end of reach k and start of reach k+1 are two stations at the same x
        m.first[r.idx] = len(m.xs)
        for k in range(ns + 1):
            x = r.x0 + r.L * k / ns
            m.xs.append(x), m.rs.append(r.idx), m.zs.append(r.z0 - r.S0 * (x - r.x0))
        m.last[r.idx] = len(m.xs) - 1
    m.N = len(m.xs)
    m.ctx = Ctx(Q, m.R, m.xs, m.rs, m.zs, m.tol)


def _controls(m, spec):
    """Seeds of the sweeps: (station, depth, supercritical?)."""
    R, N, Q, warns = m.R, m.N, m.Q, m.res.warnings
    m.seeds, m.ctrl = [], set()
    yu, m.tu = _bc_depth(spec.get("upstream"), "upstream", R, Q)
    yd, m.td = _bc_depth(spec.get("downstream"), "downstream", R, Q)
    r0, rl = R[0], R[-1]
    if yu is None:
        pass
    elif (yu < r0.yc * (1 - 1e-6)) if m.tu == "normal" else (yu <= r0.yc * (1 + 1e-6)):
        m.seeds.append((0, min(yu, r0.yc), True))               # supercritical start (gate, critical inlet, steep normal)
    elif m.tu not in ("normal", "critical"):
        warns.append(f"The upstream {m.tu} gives a subcritical depth ({yu:.3f} m). Subcritical flow is controlled from "
                     "downstream, so this upstream condition has no effect and was ignored.")
    if yd is None:
        pass
    elif yd >= rl.yc * (1 - 1e-6):
        m.seeds.append((N - 1, max(yd, rl.yc), False))          # subcritical tailwater control
    elif m.td != "normal":
        if rl.kind != "sub":
            warns.append(f"The downstream depth ({yd:.3f} m) is below critical depth ({rl.yc:.3f} m). Supercritical flow "
                         "cannot be affected from downstream, so it was ignored.")
        else:
            yf = rl.yn if rl.yn is not None else rl.yc
            warns.append(f"The downstream depth ({yd:.3f} m) is below critical depth ({rl.yc:.3f} m) in a reach whose flow "
                         f"should be subcritical. It was replaced by {'normal' if rl.yn is not None else 'critical'} depth ({yf:.3f} m).")
            m.seeds.append((N - 1, yf, False))
    rank = {"sub": 0, "crit": 1, "sup": 2}
    for k in range(len(R) - 1):
        if rank[R[k].kind] < rank[R[k + 1].kind]:               # e.g. mild -> steep: critical depth at the break
            i = m.first[k + 1]
            m.ctrl.add(i)
            m.seeds += [(i, R[k + 1].yc, True), (i, R[k + 1].yc, False)]


def _sweeps(m):
    m.sup, m.sub, m.fine = [None] * m.N, [None] * m.N, set()
    for i, y, is_sup in m.seeds:
        ys, _stop, fn = sweep(m.ctx, i, y, 1 if is_sup else -1, is_sup)
        m.fine |= fn
        arr = m.sup if is_sup else m.sub
        for j, v in ys.items():                                  # where sweeps overlap keep the more restrictive depth
            if arr[j] is None or (is_sup and v < arr[j]) or (not is_sup and v > arr[j]):
                arr[j] = v


def _assemble(m):
    """Walk downstream along the supercritical profile; it switches to the subcritical one at a jump."""
    R, Q, N, xs, rs, sup, sub = m.R, m.Q, m.N, m.xs, m.rs, m.sup, m.sub
    res, warns = m.res, m.res.warnings
    M = lambda i, y: hy.momentum(R[rs[i]].sec, Q, y)
    m.y, m.state, m.jb, undefined, st = [NAN] * N, [None] * N, {}, [], None
    for i in range(N):
        su, sb = sup[i], sub[i]
        if i == 0:
            if su is not None and sb is not None and M(0, su) <= M(0, sb) * (1 + 1e-10):     # drowned at the start
                st = "sub"
                y2 = hy.sequent_depth(R[0].sec, Q, su, R[0].yc)
                res.jumps.append(dict(status="drowned", x=0.0, reach=0, y1=su, y2=y2 or sb, y_tail=sb))
            else:
                st = "sup" if su is not None else ("sub" if sb is not None else None)
        elif i in m.ctrl and su is not None:
            st = "sup"                                                                     # critical control at a break
        elif st == "sup":
            if su is None or (sb is not None and M(i, su) <= M(i, sb) * (1 + 1e-10)):
                if sb is None:
                    undefined.append(xs[i])
                else:
                    m.jb[i] = _locate(m, i, forced=su is None)
                    st = "sub"
        elif st == "sub":
            if sb is None and su is not None:
                warns.append(f"Near x = {xs[i]:.1f} m the subcritical profile reaches critical depth going upstream; "
                             "the profile continues as supercritical flow.")
                st = "sup"
        elif st is None:
            st = "sup" if su is not None else ("sub" if sb is not None else None)
        m.state[i] = st
        v = su if st == "sup" else (sb if st == "sub" else None)
        if v is None and st is not None:
            undefined.append(xs[i])
        elif v is not None:
            m.y[i] = v
    if undefined:
        warns.append(f"The profile could not be computed between x = {min(undefined):.1f} m and {max(undefined):.1f} m: "
                     "the flow would have to pass through critical depth in a way this method cannot represent "
                     "(for example choking, or a profile that reaches critical depth upstream of a control).")
    rl = R[-1]
    if m.state[N - 1] == "sup" and sub[N - 1] is not None:      # supercritical to the end although a tailwater exists
        y1e = m.y[N - 1]
        y2 = hy.sequent_depth(rl.sec, Q, min(y1e, rl.yc * (1 - 1e-9)), rl.yc)
        res.jumps.append(dict(status="swept", x=xs[N - 1], reach=len(R) - 1, y1=y1e, y2=y2 or rl.sec.ymax, y_tail=sub[N - 1]))
    res.curves = {"sup": [(xs[i], sup[i]) for i in range(N) if sup[i] is not None],
                  "sub": [(xs[i], sub[i]) for i in range(N) if sub[i] is not None]}


def _locate(m, i, forced):
    """Exact toe of a jump between stations i-1 and i: where M(supercritical) = M(subcritical)."""
    R, Q, xs, rs, sup, sub, ctx = m.R, m.Q, m.xs, m.rs, m.sup, m.sub, m.ctx
    a, b = i - 1, i
    ra, rb = rs[a], rs[b]
    if forced:
        x, ri, y1 = xs[a], ra, sup[a]
        m.res.warnings.append(f"The supercritical profile reaches critical depth near x = {x:.1f} m; "
                              "the jump position there is approximate.")
    elif ra == rb and xs[b] > xs[a] and sub[b] is not None:
        rr = R[ra]
        pa, pb = ctx.pt(a), ctx.pt(b)
        at = lambda t: (t, rr.z0 - rr.S0 * (t - rr.x0), ra)
        y1f = lambda t: solve_step(ctx, pa, sup[a], at(t), True)      # supercritical depth at t (step from a)
        y2f = lambda t: solve_step(ctx, pb, sub[b], at(t), False)     # subcritical depth at t (step from b)

        def g(t):
            u, w = y1f(t), y2f(t)
            return None if (u is None or w is None) else hy.momentum(rr.sec, Q, u) - hy.momentum(rr.sec, Q, w)

        ga, gb = g(xs[a]), g(xs[b])
        x = xs[b]
        if ga is not None and gb is not None and ga > 0 >= gb:
            x = brentq(lambda t: (lambda v: -1.0 if v is None else v)(g(t)), xs[a], xs[b], xtol=1e-9)
        y1, ri = y1f(x), ra
        if y1 is None:
            x, y1 = xs[b], sup[b]
    else:
        x, ri, y1 = xs[b], rb, sup[b]
    r = R[ri]
    y1 = min(y1, r.yc * (1 - 1e-9))
    y2 = hy.sequent_depth(r.sec, Q, y1, r.yc)
    if y2 is None:
        y2 = r.sec.ymax
        m.res.warnings.append("The jump would fill the pipe (pressurised flow); this is not modelled.")
    return dict(status="free", x=x, reach=ri, y1=y1, y2=y2)


def _table(m):
    """Station rows. A free jump becomes: toe row, 12 roller rows (idealised shape), end-of-roller row."""
    R, Q, xs, rs, N = m.R, m.Q, m.xs, m.rs, m.N
    rows, seg, skip_x = [], 0, -math.inf

    def row(x, ri, y, flow, note, sg, e=None, sf=None, ptype=None):
        r = R[ri]
        sec = r.sec
        z = r.z0 - r.S0 * (x - r.x0)
        v = Q / sec.area(y)
        e = y + v * v / (2 * G) if e is None else e
        return dict(x=x, reach=ri + 1, z=z, y=y, wse=z + y, v=v, fr=hy.froude(sec, Q, y), e=e,
                    sf=hy.friction_slope(sec, r.n, Q, y) if sf is None else sf, egl=z + e,
                    ptype=ptype or classify(r, y, flow), flow=flow, note=note, seg=sg, _r=ri, _ctl=False, _fine=False)

    for i in range(N):
        if i in m.jb:
            j = m.jb[i]
            ri, r = j["reach"], R[j["reach"]]
            x_end = min(j["x"] + LENGTH_FACTOR * j["y2"], r.x0 + r.L)    # a roller is not drawn across a junction
            if x_end < j["x"] + LENGTH_FACTOR * j["y2"] - 1e-9:
                m.res.warnings.append(f"The jump (length ~ {LENGTH_FACTOR * j['y2']:.1f} m) would run past the end of "
                                      f"{r.name}; the roller is drawn only to the end of the reach.")
            k = next((q for q in range(i, N) if rs[q] == ri and xs[q] >= x_end - 1e-9 and m.sub[q] is not None), None)
            y_end = j["y2"]
            if k is not None:                                            # subcritical depth where the roller ends
                y_end = m.sub[k] if abs(xs[k] - x_end) < 1e-9 else (
                    solve_step(m.ctx, m.ctx.pt(k), m.sub[k], (x_end, r.z0 - r.S0 * (x_end - r.x0), ri), False) or m.sub[k])
            e1, e_end = hy.energy(r.sec, Q, j["y1"]), hy.energy(r.sec, Q, y_end)
            rows.append(row(j["x"], ri, j["y1"], "supercritical", "jump toe (start of roller)", seg))
            seg += 1
            for q in range(1, ROLLER_STATIONS + 1):
                s = q / ROLLER_STATIONS
                last = q == ROLLER_STATIONS
                rows.append(row(j["x"] + s * (x_end - j["x"]), ri, j["y1"] + (y_end - j["y1"]) * smooth(s),
                                "subcritical" if last else "jump (transition)",
                                "end of roller" if last else "jump roller (idealised shape)", seg,
                                e=e1 + (e_end - e1) * smooth(s), sf=None if last else NAN,
                                ptype=None if last else "jump"))
            seg += 1
            j.update(x_end=x_end, y_end=y_end)
            m.res.jumps.append(j)
            skip_x = x_end
        if math.isnan(m.y[i]) or xs[i] < skip_x - 1e-9:
            continue
        r = rs[i]
        note = ""
        if i == m.first[r] and r > 0:
            note = f"start of reach {r + 1}"
        elif i == m.last[r] and r < len(R) - 1:
            note = f"end of reach {r + 1}"
        if i in m.ctrl:
            note = (note + "; " if note else "") + "critical depth control"
        rw = row(xs[i], r, m.y[i], "supercritical" if m.state[i] == "sup" else "subcritical", note, seg)
        rw["_ctl"] = i in m.ctrl or (i == 0 and m.tu == "critical") or (i == N - 1 and m.td == "overfall")
        rw["_fine"] = i in m.fine
        rows.append(rw)
    m.res.rows = rows


def _finish(m):
    """Jump properties, warnings, reach summaries, DataFrame and self-checks."""
    from .checks import run_checks
    res, R, Q, rows, warns = m.res, m.R, m.Q, m.res.rows, m.res.warnings
    for j in res.jumps:
        r = R[j["reach"]]
        j["reach_name"] = r.name
        j.update({k: v for k, v in jump_properties(r.sec, Q, min(j["y1"], r.yc * (1 - 1e-9)), j["y2"]).items()
                  if k not in ("y1", "y2")})
    if rows:
        real = [r for r in rows if r["flow"] != "jump (transition)"]
        fr = [(r["x"], r["fr"]) for r in real if 0.95 < r["fr"] < 1.05 and not r["_ctl"]]
        if fr:
            warns.append(f"Froude number is within 5% of 1 at {len(fr)} station(s) (x = {fr[0][0]:.1f} to {fr[-1][0]:.1f} m). "
                         "Results near critical depth are very sensitive; an undular jump or waves may occur in reality.")
        worst = (0.0, 0.0)
        for a, b in zip(rows[:-1], rows[1:]):
            if a["seg"] == b["seg"] and a["_r"] == b["_r"] and b["x"] > a["x"] and not (a["_ctl"] or b["_ctl"]) \
                    and "jump" not in a["note"] + b["note"] and a["flow"] != "jump (transition)":
                dy = abs(b["y"] - a["y"]) / max(a["y"], b["y"])
                if dy > worst[0]:
                    worst = (dy, a["x"])
        nfine = sum(1 for r in rows if r["_fine"])
        if nfine:
            warns.append(f"The step size was too coarse at {nfine} station(s) near critical depth; the solver used smaller "
                         "internal sub-steps there. Use a smaller step for a cleaner profile.")
        if worst[0] > 0.10:
            warns.append(f"Step size may be too coarse: depth changes by {100 * worst[0]:.0f}% in one step near x = "
                         f"{worst[1]:.1f} m. Use a smaller step.")
    for r in R:
        if r.sec.kind == "circular":
            if r.yn is not None and hy.capacity(r.sec, r.n, r.S0, r.sec.D) < Q:
                warns.append(f"{r.name}: Q is more than the full-pipe capacity, so the normal depth is above ~82% of the "
                             "diameter and the flow may become pressurised (not modelled).")
            if max([w["y"] for w in rows if w["_r"] == r.idx] or [0]) > 0.95 * r.sec.D:
                warns.append(f"{r.name}: depth exceeds 95% of the diameter; full-pipe flow is not modelled.")
    names = {"M": "mild", "S": "steep", "C": "critical", "H": "horizontal", "A": "adverse"}
    for r in R:
        types = []
        for w in rows:
            if w["_r"] == r.idx and w["ptype"] not in types:
                types.append(w["ptype"])
        res.reaches.append(dict(name=r.name, slope_class=f"{names[r.cls]} ({r.cls})", yn=r.yn, yc=r.yc, sc=r.sc,
                                types=types, x0=r.x0, x1=r.x0 + r.L, S0=r.S0, z0=r.z0))
    res.table = pd.DataFrame([{COLUMNS[k]: v for k, v in rw.items() if k in COLUMNS} for rw in rows])
    res._reach_objs = R
    res.checks = run_checks(res, R)
