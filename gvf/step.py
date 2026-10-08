"""The standard step method: given the depth at one station, solve the energy
equation for the depth at the next station (upstream or downstream)."""
from scipy.optimize import brentq
from . import hydraulics as hy
from .hydraulics import G


class Ctx:
    """Everything the stepping needs: discharge, reaches and the station grid."""

    def __init__(self, Q, reaches, xs, rs, zs, tol):
        self.Q, self.reaches, self.xs, self.rs, self.zs, self.tol = Q, reaches, xs, rs, zs, tol

    def pt(self, i):
        return (self.xs[i], self.zs[i], self.rs[i])


def solve_step(ctx, pa, ya, pb, sup):
    """Depth at point pb given depth ya at pa. sup=True looks for the supercritical
    root (below critical depth), otherwise the subcritical one. None = no solution."""
    Q = ctx.Q
    ra, rb = ctx.reaches[pa[2]], ctx.reaches[pb[2]]
    L = abs(pb[0] - pa[0])
    d = 1.0 if pb[0] > pa[0] else -1.0
    Aa = ra.sec.area(ya)
    Ha = pa[1] + ya + Q * Q / (2 * G * Aa * Aa)
    Sa = hy.friction_slope(ra.sec, ra.n, Q, ya)
    sec, zb, n = rb.sec, pb[1], rb.n

    def F(y):
        A = sec.area(y)
        sfb = hy.friction_slope(sec, n, Q, y)
        return zb + y + Q * Q / (2 * G * A * A) + d * 0.5 * (Sa + sfb) * L - Ha

    yc = rb.yc
    Fc = F(yc)
    if abs(Fc) <= 1e-9:
        return yc
    if Fc > 0:
        return None  # not enough energy to exist on this branch
    if sup:
        hi, lo = yc, yc * 0.9
        while F(lo) < 0:
            hi, lo = lo, lo * 0.5
            if lo < 1e-9:
                return None
    else:
        lo, hi = yc, yc * 1.1
        cap = sec.ymax
        while F(min(hi, cap)) < 0:
            if hi >= cap or hi > 1e6:
                return None
            lo, hi = hi, hi * 1.6
        hi = min(hi, cap)
    return brentq(F, lo, hi, xtol=ctx.tol)


def try_step(ctx, i, y, j, sup):
    pa, pb = ctx.pt(i), ctx.pt(j)
    r = solve_step(ctx, pa, y, pb, sup)
    if r is not None or pa[0] == pb[0]:
        return r, False
    for k in range(1, 8):  # no root with this step: try smaller sub-steps
        n = 2 ** k
        yy, p, ok = y, pa, True
        for s in range(1, n + 1):
            t = s / n
            q = (pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pb[2])
            yy = solve_step(ctx, p, yy, q, sup)
            if yy is None:
                ok = False
                break
            p = q
        if ok:
            return yy, True
    return None, False


def sweep(ctx, i0, y0, dirn, sup):
    """March from station i0 (known depth y0) in direction dirn (+1 downstream,
    -1 upstream). Returns {station: depth}, the station where it stopped (or None) and the set of
    stations that needed sub-steps."""
    ys, i, y, stop, fine = {i0: y0}, i0, y0, None, set()
    N = len(ctx.xs)
    while 0 <= i + dirn < N:
        j = i + dirn
        r, sub = try_step(ctx, i, y, j, sup)
        if r is None:
            stop = j
            break
        if sub:
            fine.add(j)  # this step needed smaller sub-steps (step size too coarse here)
        ys[j] = y = r
        i = j
    return ys, stop, fine
