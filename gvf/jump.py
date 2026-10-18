"""Hydraulic jump: classification, energy loss, length and the idealised roller shape.

Only y1, y2 (momentum equation), the energy loss and the length are computed from
physics/empirical formulas. The surface shape inside the roller is an idealised S-curve."""
from . import hydraulics as hy

# Classes by upstream Froude number Fr1 (commonly quoted USBR / Chow classes)
JUMP_CLASSES = ((1.7, "undular"), (2.5, "weak"), (4.5, "oscillating"), (9.0, "steady"), (float("inf"), "strong"))
LENGTH_FACTOR = 6.1   # L = 6.1 * y2
LENGTH_NOTE = ("Jump length L = 6.1 * y2 is an approximate empirical rule of thumb "
               "(commonly quoted for Fr1 above about 4.5); check it against your textbook.")
ROLLER_STATIONS = 12  # number of intervals used to draw the roller in the station table


def jump_type(fr1):
    for limit, name in JUMP_CLASSES:
        if fr1 < limit:
            return name


def jump_properties(sec, Q, y1, y2):
    """Fr1, Fr2, class, energy before/after, energy loss and length of a jump y1 -> y2."""
    e1, e2 = hy.energy(sec, Q, y1), hy.energy(sec, Q, y2)
    fr1 = hy.froude(sec, Q, y1)
    return dict(y1=y1, y2=y2, fr1=fr1, fr2=hy.froude(sec, Q, y2), type=jump_type(fr1),
                e1=e1, e2=e2, loss=e1 - e2, loss_pct=100.0 * (e1 - e2) / e1,
                length=LENGTH_FACTOR * y2)


def smooth(s):
    """S-curve 0 -> 1 used for the idealised surface and energy line inside the roller."""
    return s * s * (3.0 - 2.0 * s)
