"""Basic open-channel hydraulics (SI units): normal depth, critical depth,
Froude number, energy, friction slope, momentum function, sequent depth."""
import math
from scipy.optimize import brentq
from .sections import ChannelError

G = 9.81


def froude(sec, Q, y):
    A = sec.area(y)
    return Q / (A * math.sqrt(G * A / sec.top(y)))


def energy(sec, Q, y):
    A = sec.area(y)
    return y + Q * Q / (2.0 * G * A * A)


def friction_slope(sec, n, Q, y):
    """Manning: Sf = (n Q / (A R^(2/3)))^2"""
    A = sec.area(y)
    R = A / sec.perim(y)
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


def momentum(sec, Q, y):
    """Momentum function M = Q^2/(g A) + A*ybar (per unit weight of water)."""
    return Q * Q / (G * sec.area(y)) + sec.mom(y)


def capacity(sec, n, S0, y):
    """Manning discharge at depth y for bed slope S0."""
    A = sec.area(y)
    return math.sqrt(S0) / n * A * (A / sec.perim(y)) ** (2.0 / 3.0)


def _golden_max(f, a, b, it=80):
    g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(it):
        if f(c) > f(d):
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    return (a + b) / 2


def normal_depth(sec, Q, n, S0):
    """Manning normal depth, or None when none exists (S0 <= 0)."""
    if S0 <= 0:
        return None
    if sec.kind == "circular":
        ypk = _golden_max(lambda y: capacity(sec, n, S0, y), 0.5 * sec.D, sec.D)
        qmax = capacity(sec, n, S0, ypk)
        if Q > qmax * (1 + 1e-9):
            raise ChannelError(
                f"this circular pipe cannot carry {Q:g} m3/s flowing partly full at this slope "
                f"(maximum is about {qmax:.3g} m3/s). Use a larger diameter, a steeper slope or a smaller Q.")
        return brentq(lambda y: capacity(sec, n, S0, y) - Q, 1e-9 * sec.D, ypk, xtol=1e-12)
    hi = 1.0
    while capacity(sec, n, S0, hi) < Q:
        hi *= 2.0
        if hi > 1e6:
            raise ChannelError("no normal depth found (the discharge is far too large for this channel).")
    return brentq(lambda y: capacity(sec, n, S0, y) - Q, 1e-9, hi, xtol=1e-12)


def critical_depth(sec, Q):
    """Depth where Froude number = 1, i.e. Q^2 T / (g A^3) = 1."""
    f = lambda y: Q * Q * sec.top(y) / (G * sec.area(y) ** 3) - 1.0
    if sec.kind == "circular":
        lo, hi = 1e-9 * sec.D, sec.ymax
        if f(hi) > 0:
            raise ChannelError("critical depth would be above 99% of the pipe diameter: "
                               "the pipe is too small for this discharge.")
        return brentq(f, lo, hi, xtol=1e-12)
    hi = 1.0
    while f(hi) > 0:
        hi *= 2.0
        if hi > 1e6:
            raise ChannelError("no critical depth found (the discharge is far too large for this channel).")
    return brentq(f, 1e-9, hi, xtol=1e-12)


def sequent_depth(sec, Q, y1, yc=None):
    """Conjugate (sequent) depth y2 > yc that has the same momentum function as y1 < yc.
    Returns None if the jump would fill a circular pipe."""
    if yc is None:
        yc = critical_depth(sec, Q)
    if y1 > yc * (1 + 1e-9):
        raise ChannelError("a hydraulic jump needs supercritical flow (depth below critical depth).")
    if y1 >= yc * (1 - 1e-7):
        return yc  # Fr = 1: no real jump, the sequent depth is the critical depth
    M1 = momentum(sec, Q, y1)
    f = lambda y: momentum(sec, Q, y) - M1
    lo, hi = yc, max(2.0 * yc, 2.0 * y1)
    cap = sec.ymax
    while f(min(hi, cap)) < 0:
        if hi >= cap:
            return None
        lo, hi = hi, hi * 2.0
    hi = min(hi, cap)
    return brentq(f, lo, hi, xtol=1e-12)
