"""Self-checks on a finished profile: energy equation, continuity, energy line."""
import math
from . import hydraulics as hy
from .hydraulics import G


def run_checks(res, reaches):
    """Returns the largest residuals along the profile (all in metres unless stated).

    energy_residual : |H_up - H_down - mean(Sf)*dx| between neighbouring stations (no jump between)
    continuity      : |Q - V*A| (trivial in 1D steady flow, listed for completeness), in m3/s
    egl_rise        : largest INCREASE of the energy line in the flow direction (should be ~0)
    """
    Q, rows = res.Q, res.rows
    e_res, c_res, rise = 0.0, 0.0, 0.0
    for r in rows:
        sec = reaches[r["_r"]].sec
        c_res = max(c_res, abs(Q - r["v"] * sec.area(r["y"])))
    for a, b in zip(rows[:-1], rows[1:]):
        if a["seg"] != b["seg"] or b["_fine"] or a["_fine"]:
            continue
        # the physical flow direction is always a -> b (downstream), whichever way the step was computed
        hf = 0.5 * (a["sf"] + b["sf"]) * (b["x"] - a["x"])
        H_up, H_dn = a["egl"], b["egl"]
        e_res = max(e_res, abs(H_up - H_dn - hf))
        rise = max(rise, H_dn - H_up)
    return {"energy_residual": e_res, "continuity": c_res, "egl_rise": rise}
