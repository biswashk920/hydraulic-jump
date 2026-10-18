"""Command line: python -m gvf <channel.json|channel.csv> [options]"""
import argparse
import os
from .io import load_spec, parse_bc
from .profile import compute
from .sections import ChannelError


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m gvf",
                                 description="Water surface profile (standard step method) and hydraulic jump finder.")
    ap.add_argument("channel", help="channel description: .json (everything) or .csv (reaches table)")
    ap.add_argument("--out", default="results", help="output folder (default: results)")
    ap.add_argument("--Q", type=float, help="discharge in m3/s (overrides the file; required for CSV)")
    ap.add_argument("--step", type=float, help="step size in m")
    ap.add_argument("--tol", type=float, help="depth tolerance in m")
    ap.add_argument("--start-z", type=float, help="bed elevation at the upstream end (m)")
    ap.add_argument("--up", help="upstream condition, e.g. normal, critical, none, depth:1.2, level:101.5, gate:0.4[,cc]")
    ap.add_argument("--down", help="downstream condition, e.g. normal, overfall, none, depth:2.5, level:98, weir:1.2[,cw]")
    ap.add_argument("--egl", action="store_true", help="draw the energy grade line")
    ap.add_argument("--no-plots", action="store_true", help="skip the figures")
    a = ap.parse_args(argv)
    try:
        spec = load_spec(a.channel)
        if a.Q is not None:
            spec["discharge"] = a.Q
        if a.step is not None:
            spec["step"] = a.step
        if a.tol is not None:
            spec["tolerance"] = a.tol
        if a.start_z is not None:
            spec["start_bed_elevation"] = a.start_z
        if a.up:
            spec["upstream"] = parse_bc(a.up)
        if a.down:
            spec["downstream"] = parse_bc(a.down)
        res = compute(spec)
    except ChannelError as e:
        print(f"Input problem: {e}")
        return 2
    os.makedirs(a.out, exist_ok=True)
    stem = os.path.splitext(os.path.basename(a.channel))[0]
    csv_path = os.path.join(a.out, f"{stem}_stations.csv")
    res.table.to_csv(csv_path, index=False, float_format="%.8g")
    with open(os.path.join(a.out, f"{stem}_summary.txt"), "w", encoding="utf-8") as f:
        f.write(res.summary() + "\n")
    if not a.no_plots:
        from .plots import figure_profile, figure_froude
        figure_profile(res, egl=a.egl).savefig(os.path.join(a.out, f"{stem}_profile.png"), dpi=130)
        figure_froude(res).savefig(os.path.join(a.out, f"{stem}_froude.png"), dpi=130)
    print(res.summary())
    print(f"\nSaved to {a.out}/: {stem}_stations.csv, {stem}_summary.txt" + ("" if a.no_plots else f", {stem}_profile.png, {stem}_froude.png"))
    return 0
