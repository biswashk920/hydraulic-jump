"""Independent reference solutions used for validation (not used by the solver itself).

The standard step method solves the ENERGY EQUATION between two stations. Here the same
profile is obtained a different way: the differential equation of gradually varied flow

        dy/dx = (S0 - Sf) / (1 - Fr^2)

is integrated with an adaptive Runge-Kutta solver (very small error). If both agree, the
step method is implemented correctly."""
from scipy.integrate import solve_ivp
from . import hydraulics as hy


def ode_profile(sec, Q, n, S0, y_start, x_start, x_end, xs):
    """Depths at the stations xs. Integrates from x_start (known depth y_start) towards x_end
    (x_end < x_start means upstream). Stops when the flow gets within 0.1% of critical.
    Returns a list of (x, y) for the stations that were reached."""
    def f(x, y):
        fr = hy.froude(sec, Q, y[0])
        return [(S0 - hy.friction_slope(sec, n, Q, y[0])) / (1.0 - fr * fr)]

    def near_critical(x, y):
        return abs(1.0 - hy.froude(sec, Q, y[0]) ** 2) - 1e-3
    near_critical.terminal = True
    lo, hi = sorted((x_start, x_end))
    pts = sorted([x for x in xs if lo - 1e-9 <= x <= hi + 1e-9], reverse=x_end < x_start)
    sol = solve_ivp(f, (x_start, x_end), [y_start], t_eval=pts, rtol=1e-10, atol=1e-12, method="DOP853",
                    events=near_critical)
    return list(zip(sol.t, sol.y[0]))
