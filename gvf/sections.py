"""Channel cross-sections.

Every section provides the same small set of functions of the flow depth y:
area, top width, wetted perimeter, and 'mom' = A * (depth of the centroid),
which is the second term of the momentum function used for hydraulic jumps.
"""
import math


class ChannelError(ValueError):
    """Raised with a plain-language message when the input cannot be solved."""


def _num(v, name, positive=True, zero_ok=False):
    try:
        v = float(v)
    except (TypeError, ValueError):
        raise ChannelError(f"{name} must be a number.")
    if math.isnan(v) or math.isinf(v):
        raise ChannelError(f"{name} must be a finite number.")
    if positive and (v < 0 or (v == 0 and not zero_ok)):
        raise ChannelError(f"{name} must be greater than zero (you gave {v:g}).")
    return v


class Rectangular:
    kind = "rectangular"
    ymax = math.inf

    def __init__(self, width):
        self.b = _num(width, "Rectangular channel width")

    def area(self, y):
        return self.b * y

    def top(self, y):
        return self.b

    def perim(self, y):
        return self.b + 2.0 * y

    def mom(self, y):
        return self.b * y * y / 2.0


class Trapezoidal:
    kind = "trapezoidal"
    ymax = math.inf

    def __init__(self, width, side_slope):
        self.b = _num(width, "Trapezoid bottom width", zero_ok=True)
        self.z = _num(side_slope, "Trapezoid side slope (horizontal per 1 vertical)", zero_ok=True)
        if self.b == 0 and self.z == 0:
            raise ChannelError("A trapezoid needs a bottom width or a side slope (both are zero).")
        self.k = math.sqrt(1.0 + self.z * self.z)

    def area(self, y):
        return (self.b + self.z * y) * y

    def top(self, y):
        return self.b + 2.0 * self.z * y

    def perim(self, y):
        return self.b + 2.0 * y * self.k

    def mom(self, y):
        return self.b * y * y / 2.0 + self.z * y ** 3 / 3.0


class Circular:
    """Partly full circular pipe. Depth is limited to 99% of the diameter
    (full-pipe / pressurised flow is not modelled)."""
    kind = "circular"

    def __init__(self, diameter):
        self.D = _num(diameter, "Circular diameter")
        self.r = self.D / 2.0
        self.ymax = 0.99 * self.D

    def _phi(self, y):
        y = min(max(y, 1e-12), self.D)
        c = 1.0 - 2.0 * y / self.D
        return math.acos(max(-1.0, min(1.0, c)))  # half the angle at the centre

    def area(self, y):
        p = self._phi(y)
        return self.r ** 2 * (p - math.sin(p) * math.cos(p))

    def top(self, y):
        return 2.0 * self.r * math.sin(self._phi(y))

    def perim(self, y):
        return 2.0 * self.r * self._phi(y)

    def mom(self, y):
        p = self._phi(y)
        s, c = math.sin(p), math.cos(p)
        return (2.0 / 3.0) * self.r ** 3 * s ** 3 - self.r ** 3 * c * (p - s * c)


def make_section(d):
    if not isinstance(d, dict) or "type" not in d:
        raise ChannelError("each reach needs a 'section' with a 'type' (rectangular, trapezoidal or circular).")
    t = str(d["type"]).lower()
    if t.startswith("rect"):
        return Rectangular(d.get("width"))
    if t.startswith("trap"):
        return Trapezoidal(d.get("width", 0), d.get("side_slope", 0))
    if t.startswith("circ"):
        return Circular(d.get("diameter"))
    raise ChannelError(f"unknown section type '{d['type']}'. Use rectangular, trapezoidal or circular.")
