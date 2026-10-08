"""Hydraulic jump properties."""
from . import hydraulics as hy

# Fr1 limits for the jump classes (commonly quoted USBR / Chow classes).
JUMP_CLASSES = ((1.7, "undular"), (2.5, "weak"), (4.5, "oscillating"), (9.0, "steady"), (float("inf"), "strong"))
LENGTH_FACTOR = 6.1
LENGTH_NOTE = ("Jump length L = 6.1 * y2 is an approximate empirical rule of thumb "
               "(commonly quoted for Fr1 above about 4.5); check it against your textbook.")


def jump_type(fr1):
    for limit, name in JUMP_CLASSES:
        if fr1 < limit:
            return name


def jump_properties(sec, Q, y1, y2):
    """Properties of a jump from depth y1 (supercritical) to y2 (sequent depth)."""
    e1 = hy.energy(sec, Q, y1)
    e2 = hy.energy(sec, Q, y2)
    fr1 = hy.froude(sec, Q, y1)
    return dict(y1=y1, y2=y2, fr1=fr1, fr2=hy.froude(sec, Q, y2), type=jump_type(fr1),
                e1=e1, e2=e2, loss=e1 - e2, loss_pct=100.0 * (e1 - e2) / e1,
                length=LENGTH_FACTOR * y2)
