"""Matplotlib figures: longitudinal profile and Froude number."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def box_text(res):
    """Short summary used in the figure box."""
    L = []
    for r in res.reaches:
        yn = "-" if r["yn"] is None else f"{r['yn']:.2f}"
        L.append(f"{r['name']}: {r['slope_class']}; yn = {yn} m, yc = {r['yc']:.2f} m; profiles {', '.join(r['types'])}")
    for j in res.jumps:
        if j["status"] == "free":
            L.append(f"Jump at x = {j['x']:.1f} m: y1 = {j['y1']:.3f} m, y2 = {j['y2']:.3f} m, Fr1 = {j['fr1']:.2f} "
                     f"({j['type']}), energy loss = {j['loss']:.3f} m, length ~ {j['length']:.1f} m (approx.)")
        elif j["status"] == "drowned":
            L.append(f"Jump drowned at the upstream end (tailwater {j['y_tail']:.2f} m > sequent depth {j['y2']:.2f} m)")
        else:
            L.append(f"Jump swept out downstream (tailwater {j['y_tail']:.2f} m < sequent depth {j['y2']:.2f} m)")
    if not res.jumps:
        L.append("No hydraulic jump")
    return "\n".join(L)


def figure_profile(res, egl=True):
    rows, reaches = res.rows, res._reach_objs
    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = [r["x"] for r in rows]
    ax.plot(x, [r["wse"] for r in rows], color="tab:blue", lw=2, label="Water surface")
    if egl:
        ax.plot(x, [r["egl"] for r in rows], color="tab:orange", lw=1.2, label="Energy grade line (EGL)")
    zmin = min(min(r.z0, r.z0 - r.S0 * r.L) for r in reaches)
    for k, r in enumerate(reaches):
        x1, z1 = r.x0 + r.L, r.z0 - r.S0 * r.L
        ax.plot([r.x0, x1], [r.z0, z1], color="black", lw=2, label="Bed" if k == 0 else None)
        ax.fill_between([r.x0, x1], [zmin - 0.05 * (r.z0 - zmin + 1)] * 2, [r.z0, z1], color="0.85")
        if r.yn is not None:
            ax.plot([r.x0, x1], [r.z0 + r.yn, z1 + r.yn], "--", color="tab:green", lw=1,
                    label="Normal depth" if k == 0 else None)
        ax.plot([r.x0, x1], [r.z0 + r.yc, z1 + r.yc], ":", color="tab:red", lw=1.3,
                label="Critical depth" if k == 0 else None)
        if k < len(reaches) - 1:
            ax.axvline(x1, color="0.5", lw=0.8, ls="-.")
            ax.text(x1, ax.get_ylim()[1], " reach boundary", rotation=90, va="top", fontsize=7, color="0.4")
    for j in res.jumps:
        if j["status"] == "free":
            ax.axvline(j["x"], color="tab:purple", lw=1.5, label="Hydraulic jump")
    ax.set_xlabel("Distance along channel (m)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title(f"{res.name}: water surface profile (Q = {res.Q:g} m³/s)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.subplots_adjust(bottom=0.27)
    fig.text(0.01, 0.01, box_text(res), family="monospace", fontsize=7.5, va="bottom",
             bbox=dict(boxstyle="round", fc="white", ec="0.5"))
    return fig


def figure_froude(res):
    rows = res.rows
    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.plot([r["x"] for r in rows], [r["fr"] for r in rows], color="tab:blue", lw=2, label="Froude number")
    ax.axhline(1.0, color="k", ls="--", lw=1, label="Fr = 1 (critical)")
    for r in res.reaches[:-1]:
        ax.axvline(r["x1"], color="0.5", lw=0.8, ls="-.")
    for j in res.jumps:
        if j["status"] == "free":
            ax.axvline(j["x"], color="tab:purple", lw=1.5, label="Hydraulic jump")
    ax.set_xlabel("Distance along channel (m)")
    ax.set_ylabel("Froude number")
    ax.set_title("Froude number along the channel (above 1 = supercritical, below 1 = subcritical)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig
