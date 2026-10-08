"""Tests that can be verified analytically (no published textbook numbers are used).

Run from the project root:   python -m unittest discover -s tests -v    (or: pytest)
"""
import copy
import json
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gvf import compute, ChannelError
from gvf import hydraulics as hy
from gvf.sections import Rectangular, Trapezoidal, Circular

G = hy.G
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def load(name):
    with open(os.path.join(EX, name + ".json"), encoding="utf-8") as f:
        return json.load(f)


def rect(b=5.0):
    return {"type": "rectangular", "width": b}


def chan(reaches, up=None, down=None, Q=12.0, step=2.0, tol=1e-9, **kw):
    d = {"name": "test", "discharge": Q, "step": step, "tolerance": tol, "reaches": reaches}
    if up:
        d["upstream"] = up
    if down:
        d["downstream"] = down
    d.update(kw)
    return d


def reach(L, S0, n=0.013, sec=None):
    return {"length": L, "slope": S0, "n": n, "section": sec or rect()}


def types(res):
    return {t for r in res.reaches for t in r["types"]}


SECTIONS = [Rectangular(5.0), Trapezoidal(4.0, 2.0), Trapezoidal(0.0, 1.5), Circular(2.0)]


class TestNormalAndCritical(unittest.TestCase):
    def test_normal_depth_round_trip(self):
        for sec in SECTIONS:
            for Q in (0.5, 2.0):
                yn = hy.normal_depth(sec, Q, 0.014, 0.002)
                q_back = hy.capacity(sec, 0.014, 0.002, yn)
                self.assertAlmostEqual(q_back, Q, places=9, msg=sec.kind)

    def test_no_normal_depth_for_flat_and_adverse(self):
        for S0 in (0.0, -0.001):
            self.assertIsNone(hy.normal_depth(Rectangular(5), 10, 0.015, S0))

    def test_circular_cannot_carry_q(self):
        with self.assertRaises(ChannelError) as cm:
            hy.normal_depth(Circular(0.5), 5.0, 0.013, 0.001)
        self.assertIn("cannot carry", str(cm.exception))

    def test_rectangular_critical_depth_formula(self):
        for b, Q in ((5.0, 12.0), (2.0, 1.0), (10.0, 80.0)):
            q = Q / b
            self.assertAlmostEqual(hy.critical_depth(Rectangular(b), Q), (q * q / G) ** (1 / 3), places=9)

    def test_froude_is_one_at_critical_depth_all_sections(self):
        for sec in SECTIONS:
            yc = hy.critical_depth(sec, 3.0)
            self.assertAlmostEqual(hy.froude(sec, 3.0, yc), 1.0, places=8, msg=sec.kind)

    def test_slope_classes(self):
        Q = 12.0
        cases = {0.001: "M", 0.02: "S", -0.001: "A", 0.0: "H"}
        for S0, letter in cases.items():
            res = compute(chan([reach(500, S0)], up={"type": "none"}, down={"type": "overfall"}, Q=Q))
            self.assertIn(f"({letter})", res.reaches[0]["slope_class"])
        yc = hy.critical_depth(Rectangular(5), Q)
        Sc = hy.friction_slope(Rectangular(5), 0.013, Q, yc)
        res = compute(chan([reach(500, Sc)], down={"type": "depth", "depth": yc * 1.3}))
        self.assertIn("(C)", res.reaches[0]["slope_class"])


class TestSequentDepth(unittest.TestCase):
    def test_rectangular_closed_form(self):
        sec, Q = Rectangular(5.0), 12.0
        for y1 in (0.2, 0.35, 0.5, 0.7):
            q = Q / 5.0
            fr1 = q / math.sqrt(G * y1 ** 3)
            y2_exact = 0.5 * y1 * (math.sqrt(1 + 8 * fr1 ** 2) - 1)
            self.assertAlmostEqual(hy.sequent_depth(sec, Q, y1), y2_exact, places=9)

    def test_energy_loss_formula_rectangular(self):
        sec, Q, y1 = Rectangular(5.0), 12.0, 0.4
        y2 = hy.sequent_depth(sec, Q, y1)
        loss = hy.energy(sec, Q, y1) - hy.energy(sec, Q, y2)
        self.assertAlmostEqual(loss, (y2 - y1) ** 3 / (4 * y1 * y2), places=9)

    def test_momentum_equal_all_sections(self):
        for sec in SECTIONS:
            Q = 3.0
            yc = hy.critical_depth(sec, Q)
            y1 = 0.6 * yc
            y2 = hy.sequent_depth(sec, Q, y1, yc)
            self.assertGreater(y2, yc)
            self.assertAlmostEqual(hy.momentum(sec, Q, y1) / hy.momentum(sec, Q, y2), 1.0, places=10, msg=sec.kind)

    def test_jump_needs_supercritical_input(self):
        with self.assertRaises(ChannelError):
            hy.sequent_depth(Rectangular(5), 12.0, 2.0)


class TestProfiles(unittest.TestCase):
    def test_approaches_normal_depth_far_from_control(self):
        # long mild channel: high tailwater (M1) and low tailwater/free overfall (M2)
        for down in ({"type": "depth", "depth": 3.0}, {"type": "overfall"}):
            res = compute(chan([reach(20000, 0.0005, 0.02)], down=down, Q=10.0, step=20.0, tol=1e-10))
            yn = res.reaches[0]["yn"]
            far_upstream = res.rows[0]
            self.assertLess(abs(far_upstream["y"] - yn) / yn, 1e-3, down)
        # steep channel entered at critical depth (S2)
        res = compute(chan([reach(3000, 0.02, 0.015)], up={"type": "critical"}, down={"type": "none"}, step=5.0, tol=1e-10))
        self.assertLess(abs(res.rows[-1]["y"] - res.reaches[0]["yn"]) / res.reaches[0]["yn"], 1e-3)

    def test_standard_step_matches_direct_step(self):
        sec, Q, n, S0 = Rectangular(5.0), 10.0, 0.015, 0.001
        yn = hy.normal_depth(sec, Q, n, S0)
        res = compute(chan([reach(4000, S0, n)], down={"type": "depth", "depth": 2.5}, Q=Q, step=10.0, tol=1e-10))
        xs = [r["x"] for r in res.rows]
        ys = [r["y"] for r in res.rows]
        # direct step: choose depths, compute distances (M1 backwater, marching upstream from y=2.5 at x=4000)
        x, y, E = 4000.0, 2.5, hy.energy(sec, Q, 2.5)
        worst = 0.0
        while y > yn * 1.02:
            y1 = y - 0.01
            E1 = hy.energy(sec, Q, y1)
            sf = 0.5 * (hy.friction_slope(sec, n, Q, y) + hy.friction_slope(sec, n, Q, y1))
            dx = (E - E1) / (S0 - sf)  # distance between depth y1 (upstream) and y (downstream)
            x, y, E = x - dx, y1, E1
            y_std = np.interp(x, xs, ys)
            worst = max(worst, abs(y_std - y))
        print(f"[direct step vs standard step: max depth difference {worst:.2e} m]")
        self.assertLess(worst, 2e-3)

    def test_mass_and_energy_consistency(self):
        for f in ("jump_gate_mild", "mild_to_steep", "steep_to_mild_jump", "circular_culvert", "trapezoid_changing"):
            spec = load(f)
            spec["tolerance"] = 1e-10
            res = compute(spec)
            self.assertLess(res.checks["energy_residual"], 1e-6, f)
            self.assertLess(res.checks["continuity"], 1e-10, f)
            self.assertLess(res.checks["egl_rise"], 1e-6, f)

    def test_profile_types(self):
        sup0 = {"type": "none"}
        yc = hy.critical_depth(Rectangular(5), 12.0)
        Sc = hy.friction_slope(Rectangular(5), 0.013, 12.0, yc)
        cases = [
            ("M1", chan([reach(1500, 0.001)], down={"type": "depth", "depth": 3.0}, step=10.0), None),
            ("M2", chan([reach(1500, 0.001)], down={"type": "overfall"}, step=10.0), None),
            ("M3", chan([reach(100, 0.001)], up={"type": "gate", "opening": 0.5}, down={"type": "depth", "depth": 1.15}), None),
            ("S1", chan([reach(300, 0.02)], up={"type": "critical"}, down={"type": "depth", "depth": 2.0}), None),
            ("S2", chan([reach(300, 0.02)], up={"type": "critical"}, down={"type": "none"}), None),
            ("S3", chan([reach(300, 0.02)], up={"type": "depth", "depth": 0.3}, down={"type": "none"}), None),
            ("C1", chan([reach(300, Sc)], down={"type": "depth", "depth": yc * 1.4}), None),
            ("C3", chan([reach(300, Sc)], up={"type": "depth", "depth": yc * 0.6}, down={"type": "none"}), None),
            ("H2", chan([reach(300, 0.0)], up=sup0, down={"type": "overfall"}), None),
            ("H3", chan([reach(300, 0.0)], up={"type": "depth", "depth": 0.4}, down={"type": "none"}), None),
            ("A2", chan([reach(300, -0.001)], up=sup0, down={"type": "overfall"}), None),
            ("A3", chan([reach(300, -0.001)], up={"type": "depth", "depth": 0.4}, down={"type": "none"}), None),
        ]
        for want, spec, _ in cases:
            res = compute(spec)
            self.assertIn(want, types(res), f"{want}: found {types(res)}")

    def test_control_moves_to_critical_depth_at_mild_steep_break(self):
        res = compute(load("mild_to_steep"))
        brk = [r for r in res.rows if "critical depth control" in r["note"]]
        self.assertEqual(len(brk), 1)
        self.assertAlmostEqual(brk[0]["y"], res.reaches[1]["yc"], places=6)
        self.assertEqual(res.jumps, [])

    def test_changing_section_energy_continuous(self):
        spec = chan([reach(500, 0.0008, 0.03, {"type": "trapezoidal", "width": 5, "side_slope": 1.5}),
                     reach(500, 0.0003, 0.022, {"type": "trapezoidal", "width": 8, "side_slope": 1.5})],
                    down={"type": "depth", "depth": 2.4}, Q=15.0, step=25.0)
        res = compute(spec)
        a = [r for r in res.rows if r["note"].startswith("end of reach 1")][0]
        b = [r for r in res.rows if r["note"].startswith("start of reach 2")][0]
        self.assertAlmostEqual(a["egl"], b["egl"], places=6)
        self.assertEqual(a["x"], b["x"])


class TestJump(unittest.TestCase):
    def setUp(self):
        self.res = compute(load("jump_gate_mild"))

    def test_jump_found_with_equal_momentum_and_closed_form(self):
        self.assertEqual(len(self.res.jumps), 1)
        j = self.res.jumps[0]
        self.assertEqual(j["status"], "free")
        sec = Rectangular(5.0)
        self.assertAlmostEqual(hy.momentum(sec, 12.0, j["y1"]) / hy.momentum(sec, 12.0, j["y2"]), 1.0, places=9)
        fr1 = hy.froude(sec, 12.0, j["y1"])
        self.assertAlmostEqual(j["y2"] / j["y1"], 0.5 * (math.sqrt(1 + 8 * fr1 ** 2) - 1), places=9)
        self.assertAlmostEqual(j["loss"], (j["y2"] - j["y1"]) ** 3 / (4 * j["y1"] * j["y2"]), places=9)

    def test_jump_is_where_supercritical_meets_subcritical_profile(self):
        # the table has a 'before' and 'after' row at the jump; downstream of it the water is subcritical
        rows = self.res.rows
        k = [i for i, r in enumerate(rows) if r["note"] == "jump: before"][0]
        self.assertEqual(rows[k]["x"], rows[k + 1]["x"])
        self.assertTrue(all(r["fr"] > 1 for r in rows[:k + 1]))
        self.assertTrue(all(r["fr"] < 1 for r in rows[k + 1:]))

    def test_jump_classes(self):
        from gvf.jump import jump_type
        for fr, name in ((1.2, "undular"), (2.0, "weak"), (3.0, "oscillating"), (6.0, "steady"), (12.0, "strong")):
            self.assertEqual(jump_type(fr), name)

    def test_drowned_and_swept(self):
        d = compute(load("drowned_jump"))
        self.assertEqual([j["status"] for j in d.jumps], ["drowned"])
        self.assertGreater(d.jumps[0]["y_tail"], d.jumps[0]["y2"])
        self.assertTrue(all(r["fr"] < 1 for r in d.rows))  # subcritical all the way
        s = compute(load("swept_jump"))
        self.assertEqual([j["status"] for j in s.jumps], ["swept"])
        self.assertLess(s.jumps[0]["y_tail"], s.jumps[0]["y2"])
        self.assertTrue(all(r["fr"] > 1 for r in s.rows))

    def test_jump_moves_downstream_when_tailwater_drops(self):
        xs = []
        for tail in (1.5, 1.4, 1.3):
            spec = load("jump_gate_mild")
            spec["downstream"] = {"type": "depth", "depth": tail}
            xs.append(compute(spec).jumps[0]["x"])
        self.assertTrue(xs[0] < xs[1] < xs[2], xs)

    def test_circular_jump_momentum(self):
        res = compute(load("circular_culvert"))
        j = res.jumps[0]
        sec = Circular(1.2)
        self.assertAlmostEqual(hy.momentum(sec, 1.6, j["y1"]) / hy.momentum(sec, 1.6, j["y2"]), 1.0, places=9)


class TestValidation(unittest.TestCase):
    def bad(self, spec, text):
        with self.assertRaises(ChannelError) as cm:
            compute(spec)
        self.assertIn(text, str(cm.exception))

    def test_messages(self):
        base = chan([reach(100, 0.001)])
        self.bad({**base, "discharge": -1}, "Discharge")
        self.bad({**base, "discharge": "abc"}, "must be a number")
        self.bad({**base, "reaches": []}, "at least one reach")
        self.bad({**base, "reaches": [reach(-5, 0.001)]}, "length")
        self.bad({**base, "reaches": [reach(100, 0.001, n=0)]}, "Manning")
        self.bad({**base, "reaches": [reach(100, 0.001, sec=rect(-2))]}, "width")
        self.bad({**base, "reaches": [reach(100, 0.001, sec={"type": "oval"})]}, "unknown section")
        self.bad({**base, "downstream": {"type": "depth", "depth": 0}}, "depth")
        self.bad({**base, "downstream": {"type": "teleport"}}, "not available")
        self.bad({**base, "downstream": {"type": "weir", "crest": 1}, "reaches": [reach(100, 0.001, sec={"type": "trapezoidal", "width": 3, "side_slope": 1})]}, "rectangular")
        self.bad(chan([reach(100, 0.0)]), "no normal depth")  # default 'normal' BC on a horizontal reach
        self.bad({**base, "step": 0.0001, "reaches": [reach(100000, 0.001)]}, "Increase the step")
        self.bad(chan([reach(100, 0.001, sec={"type": "circular", "diameter": 0.4})], Q=5.0), "too small")
        self.bad(chan([reach(100, 0.00001, sec={"type": "circular", "diameter": 1.2})], Q=1.6), "cannot carry")

    def test_supercritical_tailwater_ignored_with_warning(self):
        res = compute(chan([reach(300, 0.02)], up={"type": "critical"}, down={"type": "depth", "depth": 0.2}))
        self.assertTrue(any("ignored" in w for w in res.warnings))

    def test_coarse_step_warning(self):
        res = compute(chan([reach(300, 0.001)], up={"type": "gate", "opening": 0.5}, down={"type": "depth", "depth": 1.4}, step=60.0))
        self.assertTrue(any("coarse" in w for w in res.warnings))


if __name__ == "__main__":
    unittest.main()
