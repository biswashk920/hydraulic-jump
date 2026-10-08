"""Builds the channel, finds the controls, runs the standard step sweeps,
puts the subcritical / supercritical profiles together and locates hydraulic jumps."""
import math
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from . import hydraulics as hy
from .hydraulics import G
from .jump import jump_properties, LENGTH_NOTE
from .sections import make_section, ChannelError
from .step import Ctx, sweep, solve_step

UP_TYPES = ("normal", "depth", "level", "gate", "critical", "none")
DOWN_TYPES = ("normal", "depth", "level", "weir", "overfall", "none")
MAX_STATIONS = 20000

COLUMNS = {  # short internal name -> readable name used in tables and CSV
    "x": "Distance (m)", "reach": "Reach", "z": "Bed elevation (m)", "y": "Depth (m)",
    "wse": "Water surface elevation (m)", "v": "Velocity (m/s)", "fr": "Froude number",
    "e": "Specific energy (m)", "sf": "Friction slope", "egl": "Energy grade line elevation (m)",
    "ptype": "Profile type", "flow": "Flow", "note": "Note"}


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
    pass


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
        r.sc = hy.friction_slope(r.sec, r.n, Q, r.yc)  # slope at which uniform flow would be critical
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


def _bc_depth(bc, side, reaches, Q, z_start, warns):
    """Depth implied by a boundary condition and its type. side = 'upstream'/'downstream'."""
    up = side == "upstream"
    r = reaches[0] if up else reaches[-1]
    bc = dict(bc or {"type": "normal"})
    t = str(bc.get("type", "normal")).lower()
    allowed = UP_TYPES if up else DOWN_TYPES
    if t not in allowed:
        raise ChannelError(f"the {side} condition '{t}' is not available. Choose one of: {', '.join(allowed)}.")
    where = f"{side} condition"
    if t == "none":
        return None, t
    if t == "normal":
        if r.yn is None:
            raise ChannelError(f"the {side} reach has a {'horizontal' if r.S0 == 0 else 'adverse'} slope, so it has no "
                               f"normal depth. Choose another {side} condition (for example a depth).")
        return r.yn, t
    if t == "critical" or t == "overfall":
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
        cc = _num(bc.get("cc", 0.61), f"{where}: contraction coefficient", positive=True)
        return cc * a, t
    if t == "weir":
        if r.sec.kind != "rectangular":
            raise ChannelError("the weir option works only for a rectangular last reach. "
                               "For other sections, enter the tailwater depth directly.")
        crest = _num(bc.get("crest"), f"{where}: weir crest height", positive=True, zero_ok=True)
        cw = _num(bc.get("cw", 1.7), f"{where}: weir coefficient", positive=True)
        h = (Q / (cw * r.sec.b)) ** (2.0 / 3.0)
        return crest + h, t


def classify(r, y, flow="subcritical"):
    """Profile type from the depth relative to normal and critical depth.
    A depth exactly at critical depth counts as 'above' for subcritical flow and
    'below' for supercritical flow."""
    yn, yc, c = r.yn, r.yc, r.cls
    if c in "MS" and abs(y - yn) <= 0.01 * yn:
        return "uniform"
    above_c = y >= yc * (1 - 1e-6) if flow == "subcritical" else y > yc * (1 + 1e-6)
    if c == "M":
        return "M1" if y > yn else ("M2" if above_c else "M3")
    if c == "S":
        return "S1" if above_c else ("S2" if y > yn else "S3")
    if c == "C":
        return "C1" if above_c else "C3"
    return "%s%d" % (c, 2 if above_c else 3)


class Result:
    def __init__(self):
        self.warnings, self.jumps, self.reaches, self.rows = [], [], [], []
        self.name, self.Q, self.table, self.checks = "channel", 0.0, None, {}

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
                L.append(f"Hydraulic jump at x = {j['x']:.2f} m ({j['reach_name']}): {j['y1']:.3f} m -> {j['y2']:.3f} m, "
                         f"Fr1 = {j['fr1']:.2f} ({j['type']}), energy loss {j['loss']:.3f} m ({j['loss_pct']:.1f}% of E1), "
                         f"length ~ {j['length']:.1f} m (approximate)")
            elif j["status"] == "drowned":
                L.append(f"Jump DROWNED: the tailwater ({j['y_tail']:.3f} m at the upstream end) is higher than the sequent depth "
                         f"({j['y2']:.3f} m) needed for y1 = {j['y1']:.3f} m, so the jump is pushed up to the upstream end (gate).")
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


def compute(spec):
    """Compute the water surface profile for a channel description (a dict, see README)."""
    from .checks import run_checks
    if not isinstance(spec, dict):
        raise ChannelError("the channel description must be a JSON object.")
    Q = _num(spec.get("discharge"), "Discharge Q", positive=True)
    z_start = _num(spec.get("start_bed_elevation", 100.0), "Starting bed elevation")
    step = _num(spec.get("step", 1.0), "Step size", positive=True)
    tol = _num(spec.get("tolerance", 1e-6), "Tolerance", positive=True)
    res = Result()
    res.name, res.Q = str(spec.get("name", "channel")), Q
    warns = res.warnings
    reaches = _build_reaches(spec, Q, z_start)
    m = len(reaches)

    # ---- station grid (the end of each reach and the start of the next share a distance)
    xs, rs, zs, first_of, last_of = [], [], [], {}, {}
    counts = [max(20, math.ceil(r.L / step - 1e-9)) for r in reaches]
    total = sum(c + 1 for c in counts)
    if total > MAX_STATIONS:
        raise ChannelError(f"that would need {total} stations. Increase the step size (limit {MAX_STATIONS} stations).")
    for r, ns in zip(reaches, counts):
        first_of[r.idx] = len(xs)
        for k in range(ns + 1):
            x = r.x0 + r.L * k / ns
            xs.append(x), rs.append(r.idx), zs.append(r.z0 - r.S0 * (x - r.x0))
        last_of[r.idx] = len(xs) - 1
    N = len(xs)
    if N > MAX_STATIONS:
        raise ChannelError(f"that would need {N} stations. Increase the step size (limit {MAX_STATIONS} stations).")
    ctx = Ctx(Q, reaches, xs, rs, zs, tol)

    # ---- controls
    seeds, ctrl_nodes = [], set()
    yu, tu = _bc_depth(spec.get("upstream"), "upstream", reaches, Q, z_start, warns)
    yd, td = _bc_depth(spec.get("downstream"), "downstream", reaches, Q, z_start, warns)
    r0, rl = reaches[0], reaches[-1]
    if yu is None:
        pass
    elif (yu < r0.yc * (1 - 1e-6)) if tu == "normal" else (yu <= r0.yc * (1 + 1e-6)):
        seeds.append((0, min(yu, r0.yc), True))
    elif tu not in ("normal", "critical"):
        warns.append(f"The upstream {tu} gives a subcritical depth ({yu:.3f} m). Subcritical flow is controlled from "
                     "downstream, so this upstream condition has no effect and was ignored.")
    if yd is None:
        pass
    elif yd >= rl.yc * (1 - 1e-6):
        seeds.append((N - 1, max(yd, rl.yc), False))
    elif td != "normal":
        if rl.kind == "sup" or rl.kind == "crit":
            warns.append(f"The downstream depth ({yd:.3f} m) is below critical depth ({rl.yc:.3f} m). Supercritical flow "
                         "cannot be affected from downstream, so it was ignored.")
        else:
            yf = rl.yn if rl.yn is not None else rl.yc
            warns.append(f"The downstream depth ({yd:.3f} m) is below critical depth ({rl.yc:.3f} m) in a reach whose flow "
                         f"should be subcritical. It was replaced by {'normal' if rl.yn is not None else 'critical'} depth "
                         f"({yf:.3f} m).")
            seeds.append((N - 1, yf, False))
    rank = {"sub": 0, "crit": 1, "sup": 2}
    for k in range(m - 1):
        if rank[reaches[k].kind] < rank[reaches[k + 1].kind]:  # e.g. mild -> steep: critical depth at the break
            i = first_of[k + 1]
            ctrl_nodes.add(i)
            seeds += [(i, reaches[k + 1].yc, True), (i, reaches[k + 1].yc, False)]

    # ---- sweeps
    sup = [None] * N
    sub = [None] * N
    fine = set()
    for i, y, is_sup in seeds:
        ys, _stop, fn = sweep(ctx, i, y, 1 if is_sup else -1, is_sup)
        fine |= fn
        arr = sup if is_sup else sub
        for j, v in ys.items():
            if arr[j] is None or (is_sup and v < arr[j]) or (not is_sup and v > arr[j]):
                arr[j] = v

    M = lambda i, y: hy.momentum(reaches[rs[i]].sec, Q, y)
    # ---- put the two profiles together
    y_final, state, jump_before, st = [math.nan] * N, [None] * N, {}, None
    undefined = []
    for i in range(N):
        su, sb = sup[i], sub[i]
        if i == 0:
            if su is not None and sb is not None and M(0, su) <= M(0, sb) * (1 + 1e-10):
                st = "sub"
                y2n = hy.sequent_depth(reaches[0].sec, Q, su, reaches[0].yc)
                res.jumps.append(dict(status="drowned", x=0.0, reach=0, y1=su, y2=y2n if y2n else sb, y_tail=sb))
            elif su is not None:
                st = "sup"
            elif sb is not None:
                st = "sub"
        elif i in ctrl_nodes and su is not None:
            if st != "sup":
                st = "sup"
        elif st == "sup":
            if su is None or (sb is not None and M(i, su) <= M(i, sb) * (1 + 1e-10)):
                if sb is None:
                    undefined.append(xs[i])
                else:
                    jump_before[i] = _locate(i, ctx, sup, sub, xs, rs, reaches, Q, warns, forced=su is None)
                    st = "sub"
        elif st == "sub":
            if sb is None and su is not None:
                warns.append(f"Near x = {xs[i]:.1f} m the subcritical profile reaches critical depth going upstream; "
                             "the profile continues as supercritical flow.")
                st = "sup"
        elif st is None:
            st = "sup" if su is not None else ("sub" if sb is not None else None)
        state[i] = st
        v = su if st == "sup" else (sb if st == "sub" else None)
        if v is None and st is not None:
            undefined.append(xs[i])
        else:
            y_final[i] = v if v is not None else math.nan
    if undefined:
        warns.append(f"The profile could not be computed between x = {min(undefined):.1f} m and {max(undefined):.1f} m: "
                     "the flow would have to pass through critical depth in a way this method cannot represent "
                     "(for example choking, or a profile that reaches critical depth upstream of a control).")
    if state[N - 1] == "sup" and sub[N - 1] is not None:  # supercritical to the end although a tailwater exists
        sec = rl.sec
        y1e = y_final[N - 1]
        y2n = hy.sequent_depth(sec, Q, min(y1e, rl.yc * (1 - 1e-9)), rl.yc)
        res.jumps.append(dict(status="swept", x=xs[N - 1], reach=m - 1, y1=y1e, y2=y2n if y2n else sec.ymax, y_tail=sub[N - 1]))

    # ---- table rows
    def row(x, ri, y, flow, note, seg):
        r = reaches[ri]
        sec = r.sec
        z = r.z0 - r.S0 * (x - r.x0)
        A = sec.area(y)
        v = Q / A
        e = y + v * v / (2 * G)
        return dict(x=x, reach=ri + 1, z=z, y=y, wse=z + y, v=v, fr=hy.froude(sec, Q, y), e=e,
                    sf=hy.friction_slope(sec, r.n, Q, y), egl=z + e, ptype=classify(r, y, flow), flow=flow,
                    note=note, seg=seg, _r=ri, _ctl=False, _fine=False)

    rows, seg = [], 0
    for i in range(N):
        if i in jump_before:
            j = jump_before[i]
            rows.append(row(j["x"], j["reach"], j["y1"], "supercritical", "jump: before", seg))
            seg += 1
            rows.append(row(j["x"], j["reach"], j["y2"], "subcritical", "jump: after", seg))
            res.jumps.append(j)
        if math.isnan(y_final[i]):
            continue
        r = rs[i]
        note = ""
        if i == first_of[r] and r > 0:
            note = f"start of reach {r + 1}"
        elif i == last_of[r] and r < m - 1:
            note = f"end of reach {r + 1}"
        if i in ctrl_nodes:
            note = (note + "; " if note else "") + "critical depth control"
        rw = row(xs[i], r, y_final[i], "supercritical" if state[i] == "sup" else "subcritical", note, seg)
        rw["_ctl"] = i in ctrl_nodes or (i == 0 and tu == "critical") or (i == N - 1 and td == "overfall")
        rw["_fine"] = i in fine
        rows.append(rw)
    res.rows = rows

    # ---- complete the jump records
    for j in res.jumps:
        r = reaches[j["reach"]]
        j["reach_name"] = r.name
        j.update({k: v for k, v in jump_properties(r.sec, Q, min(j["y1"], r.yc * (1 - 1e-9)), j["y2"]).items()
                  if k not in ("y1", "y2")})
        if j["status"] == "free":
            for q in range(m - 1):
                if abs(j["x"] - (reaches[q].x0 + reaches[q].L)) < j["length"]:
                    warns.append(f"The jump (length ~ {j['length']:.1f} m) is close to the junction between reach {q + 1} "
                                 f"and reach {q + 2}; the jump equations assume one uniform section.")
    # ---- warnings from the table
    if rows:
        fr = [(r["x"], r["fr"]) for r in rows if 0.95 < r["fr"] < 1.05 and not r["_ctl"]]
        if fr:
            warns.append(f"Froude number is within 5% of 1 at {len(fr)} station(s) (x = {fr[0][0]:.1f} to {fr[-1][0]:.1f} m). "
                         "Results near critical depth are very sensitive; an undular jump or waves may occur in reality.")
        worst = (0.0, 0.0)
        for a, b in zip(rows[:-1], rows[1:]):
            if a["seg"] == b["seg"] and a["_r"] == b["_r"] and b["x"] > a["x"] and not (a["_ctl"] or b["_ctl"]):
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
    for r in reaches:
        if r.sec.kind == "circular" and r.yn is not None and hy.capacity(r.sec, r.n, r.S0, r.sec.D) < Q:
            warns.append(f"{r.name}: Q is more than the full-pipe capacity, so the normal depth is above ~82% of the diameter "
                         "and the flow may become pressurised (not modelled).")
        if r.sec.kind == "circular":
            top = max([row_["y"] for row_ in rows if row_["_r"] == r.idx] or [0])
            if top > 0.95 * r.sec.D:
                warns.append(f"{r.name}: depth exceeds 95% of the diameter; full-pipe flow is not modelled.")

    # ---- reach summaries and table
    for r in reaches:
        types = []
        for row_ in rows:
            if row_["_r"] == r.idx and row_["ptype"] not in types:
                types.append(row_["ptype"])
        name = {"M": "mild", "S": "steep", "C": "critical", "H": "horizontal", "A": "adverse"}[r.cls]
        res.reaches.append(dict(name=r.name, slope_class=f"{name} ({r.cls})", yn=r.yn, yc=r.yc, sc=r.sc, types=types,
                                x0=r.x0, x1=r.x0 + r.L, S0=r.S0))
    res.table = pd.DataFrame([{COLUMNS[k]: v for k, v in rw.items() if k in COLUMNS} for rw in rows])
    res.checks = run_checks(res, reaches)
    res._reach_objs = reaches
    return res


def _locate(i, ctx, sup, sub, xs, rs, reaches, Q, warns, forced):
    """Exact position of a jump between stations i-1 and i (momentum balance)."""
    a, b = i - 1, i
    ra, rb = rs[a], rs[b]
    if forced:
        x, ri, y1 = xs[a], ra, sup[a]
        warns.append(f"The supercritical profile reaches critical depth near x = {x:.1f} m; "
                     "the jump position there is approximate.")
    elif ra == rb and xs[b] > xs[a] and sub[b] is not None:
        sec, rr = reaches[ra].sec, reaches[ra]
        pa, pb = ctx.pt(a), ctx.pt(b)
        at = lambda x: (x, rr.z0 - rr.S0 * (x - rr.x0), ra)
        y1f = lambda x: solve_step(ctx, pa, sup[a], at(x), True)    # supercritical depth at x (step from a)
        y2f = lambda x: solve_step(ctx, pb, sub[b], at(x), False)   # subcritical depth at x (step from b)

        def g(x):
            u, w = y1f(x), y2f(x)
            return None if (u is None or w is None) else hy.momentum(sec, Q, u) - hy.momentum(sec, Q, w)

        ga, gb = g(xs[a]), g(xs[b])
        x = xs[b]
        if ga is not None and gb is not None and ga > 0 >= gb:
            x = brentq(lambda t: (lambda v: -1.0 if v is None else v)(g(t)), xs[a], xs[b], xtol=1e-9)
        y1 = y1f(x)
        if y1 is None:
            x, y1 = xs[b], sup[b]
        ri = ra
    else:
        x, ri, y1 = xs[b], rb, sup[b]
    r = reaches[ri]
    y1 = min(y1, r.yc * (1 - 1e-9))
    y2 = hy.sequent_depth(r.sec, Q, y1, r.yc)
    if y2 is None:
        y2 = r.sec.ymax
        warns.append("The jump would fill the pipe (pressurised flow); this is not modelled.")
    return dict(status="free", x=x, reach=ri, y1=y1, y2=y2)
