"""Behaviour of whole profiles: classification, controls, jumps, warnings and input validation."""
import json
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gvf import compute, ChannelError
from gvf import hydraulics as hy
from gvf.sections import Rectangular, Circular

EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def load(name):
    with open(os.path.join(EX, name + ".json"), encoding="utf-8") as f:
        return json.load(f)


def rect(b=5.0):
    return {"type": "rectangular", "width": b}


def reach(L, S0, n=0.013, sec=None):
    return {"length": L, "slope": S0, "n": n, "section": sec or rect()}


def chan(reaches, up=None, down=None, Q=12.0, step=2.0, tol=1e-9):
    d = {"name": "test", "discharge": Q, "step": step, "tolerance": tol, "reaches": reaches}
    if up:
        d["upstream"] = up
    if down:
        d["downstream"] = down
    return d


def types(res):
    return {t for r in res.reaches for t in r["types"]}


class TestClassification(unittest.TestCase):
    def test_slope_classes(self):
        for S0, letter in {0.001: "M", 0.02: "S", -0.001: "A", 0.0: "H"}.items():
            res = compute(chan([reach(500, S0)], up={"type": "none"}, down={"type": "overfall"}))
            self.assertIn(f"({letter})", res.reaches[0]["slope_class"])
        yc = hy.critical_depth(Rectangular(5), 12.0)
        Sc = hy.friction_slope(Rectangular(5), 0.013, 12.0, yc)
        self.assertIn("(C)", compute(chan([reach(500, Sc)], down={"type": "depth", "depth": yc * 1.3})).reaches[0]["slope_class"])

    def test_every_profile_type_is_produced(self):
        none = {"type": "none"}
        yc = hy.critical_depth(Rectangular(5), 12.0)
        Sc = hy.friction_slope(Rectangular(5), 0.013, 12.0, yc)
        cases = {
            "M1": chan([reach(1500, 0.001)], down={"type": "depth", "depth": 3.0}, step=10.0),
            "M2": chan([reach(1500, 0.001)], down={"type": "overfall"}, step=10.0),
            "M3": chan([reach(100, 0.001)], up={"type": "gate", "opening": 0.5}, down={"type": "depth", "depth": 1.15}),
            "S1": chan([reach(300, 0.02)], up={"type": "critical"}, down={"type": "depth", "depth": 2.0}),
            "S2": chan([reach(300, 0.02)], up={"type": "critical"}, down=none),
            "S3": chan([reach(300, 0.02)], up={"type": "depth", "depth": 0.3}, down=none),
            "C1": chan([reach(300, Sc)], down={"type": "depth", "depth": yc * 1.4}),
            "C3": chan([reach(300, Sc)], up={"type": "depth", "depth": yc * 0.6}, down=none),
            "H2": chan([reach(300, 0.0)], up=none, down={"type": "overfall"}),
            "H3": chan([reach(300, 0.0)], up={"type": "depth", "depth": 0.4}, down=none),
            "A2": chan([reach(300, -0.001)], up=none, down={"type": "overfall"}),
            "A3": chan([reach(300, -0.001)], up={"type": "depth", "depth": 0.4}, down=none),
        }
        for want, spec in cases.items():
            self.assertIn(want, types(compute(spec)), want)


class TestControls(unittest.TestCase):
    def test_critical_depth_control_at_mild_steep_break(self):
        res = compute(load("mild_to_steep"))
        brk = [r for r in res.rows if "critical depth control" in r["note"]]
        self.assertEqual(len(brk), 1)
        self.assertAlmostEqual(brk[0]["y"], res.reaches[1]["yc"], places=6)
        self.assertEqual(res.jumps, [])

    def test_energy_is_continuous_across_a_change_of_section(self):
        spec = chan([reach(500, 0.0008, 0.03, {"type": "trapezoidal", "width": 5, "side_slope": 1.5}),
                     reach(500, 0.0003, 0.022, {"type": "trapezoidal", "width": 8, "side_slope": 1.5})],
                    down={"type": "depth", "depth": 2.4}, Q=15.0, step=25.0)
        res = compute(spec)
        a = [r for r in res.rows if r["note"].startswith("end of reach 1")][0]
        b = [r for r in res.rows if r["note"].startswith("start of reach 2")][0]
        self.assertAlmostEqual(a["egl"], b["egl"], places=6)

    def test_self_checks_on_all_examples(self):
        for f in ("jump_gate_mild", "mild_to_steep", "steep_to_mild_jump", "circular_culvert", "trapezoid_changing"):
            spec = load(f)
            spec["tolerance"] = 1e-10
            c = compute(spec).checks
            self.assertLess(c["energy_residual"], 1e-6, f)
            self.assertLess(c["continuity"], 1e-10, f)
            self.assertLess(c["egl_rise"], 1e-6, f)


class TestJump(unittest.TestCase):
    def setUp(self):
        self.res = compute(load("jump_gate_mild"))
        self.j = self.res.jumps[0]

    def test_jump_has_a_toe_a_roller_and_an_end(self):
        j = self.j
        self.assertEqual(j["status"], "free")
        self.assertAlmostEqual(j["x_end"] - j["x"], 6.1 * j["y2"], places=9)
        notes = [r["note"] for r in self.res.rows]
        self.assertEqual(sum(n.startswith("jump toe") for n in notes), 1)
        self.assertEqual(sum(n == "end of roller" for n in notes), 1)
        roller = [r for r in self.res.rows if r["flow"] == "jump (transition)"]
        self.assertGreaterEqual(len(roller), 10)
        ys = [r["y"] for r in roller]
        self.assertTrue(all(b >= a for a, b in zip(ys, ys[1:])), "depth rises through the roller")
        egl = [r["egl"] for r in self.res.rows if j["x"] <= r["x"] <= j["x_end"]]
        self.assertTrue(all(b <= a + 1e-9 for a, b in zip(egl, egl[1:])), "energy line falls through the roller")

    def test_supercritical_before_subcritical_after(self):
        rows = self.res.rows
        self.assertTrue(all(r["fr"] > 1 for r in rows if r["x"] <= self.j["x"]))
        self.assertTrue(all(r["fr"] < 1 for r in rows if r["x"] >= self.j["x_end"]))
        self.assertTrue(all(r["x"] <= self.j["x"] or r["x"] >= self.j["x_end"] or r["flow"] == "jump (transition)"
                            for r in rows), "no ordinary station inside the roller")

    def test_drowned_and_swept(self):
        d = compute(load("drowned_jump"))
        self.assertEqual([j["status"] for j in d.jumps], ["drowned"])
        self.assertGreater(d.jumps[0]["y_tail"], d.jumps[0]["y2"])
        self.assertTrue(all(r["fr"] < 1 for r in d.rows))
        s = compute(load("swept_jump"))
        self.assertEqual([j["status"] for j in s.jumps], ["swept"])
        self.assertLess(s.jumps[0]["y_tail"], s.jumps[0]["y2"])
        self.assertTrue(all(r["fr"] > 1 for r in s.rows))

    def test_roller_is_cut_at_the_end_of_a_reach(self):
        import numpy as np
        x_cut = self.j["x"] + 3.0                       # the roller is about 8 m long
        xs, ys = zip(*self.res.curves["sub"])
        spec = load("jump_gate_mild")
        spec["reaches"][0]["length"] = x_cut
        spec["downstream"] = {"type": "depth", "depth": float(np.interp(x_cut, xs, ys))}   # keeps the same jump
        res = compute(spec)
        self.assertTrue(any("run past the end" in w for w in res.warnings))
        self.assertAlmostEqual(res.jumps[0]["x_end"], x_cut, places=9)

    def test_circular_jump(self):
        res = compute(load("circular_culvert"))
        j, sec = res.jumps[0], Circular(1.2)
        self.assertAlmostEqual(hy.momentum(sec, 1.6, j["y1"]) / hy.momentum(sec, 1.6, j["y2"]), 1.0, places=9)

    def test_jump_classes(self):
        from gvf.jump import jump_type
        for fr, name in ((1.2, "undular"), (2.0, "weak"), (3.0, "oscillating"), (6.0, "steady"), (12.0, "strong")):
            self.assertEqual(jump_type(fr), name)


class TestInputValidation(unittest.TestCase):
    def bad(self, spec, text):
        with self.assertRaises(ChannelError) as cm:
            compute(spec)
        self.assertIn(text, str(cm.exception))

    def test_friendly_messages(self):
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
        self.bad({**base, "downstream": {"type": "weir", "crest": 1},
                  "reaches": [reach(100, 0.001, sec={"type": "trapezoidal", "width": 3, "side_slope": 1})]}, "rectangular")
        self.bad(chan([reach(100, 0.0)]), "no normal depth")
        self.bad({**base, "step": 0.0001, "reaches": [reach(100000, 0.001)]}, "Increase the step")
        self.bad(chan([reach(100, 0.001, sec={"type": "circular", "diameter": 0.4})], Q=5.0), "too small")
        self.bad(chan([reach(100, 0.00001, sec={"type": "circular", "diameter": 1.2})], Q=1.6), "cannot carry")

    def test_warnings(self):
        res = compute(chan([reach(300, 0.02)], up={"type": "critical"}, down={"type": "depth", "depth": 0.2}))
        self.assertTrue(any("ignored" in w for w in res.warnings))
        res = compute(chan([reach(300, 0.001)], up={"type": "gate", "opening": 0.5}, down={"type": "depth", "depth": 1.4}, step=60.0))
        self.assertTrue(any("coarse" in w for w in res.warnings))


if __name__ == "__main__":
    unittest.main()
