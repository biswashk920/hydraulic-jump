"""The browser version (docs/index.html) must give the same answers as the Python package.
Needs Node.js; the test is skipped if `node` is not installed."""
import json
import os
import shutil
import subprocess
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gvf import compute


@unittest.skipUnless(shutil.which("node"), "Node.js not installed")
class TestJsMatchesPython(unittest.TestCase):
    def test_all_examples(self):
        out = subprocess.run(["node", os.path.join(ROOT, "tests", "js_check.js")], capture_output=True, text=True, check=True)
        js = json.loads(out.stdout)
        self.assertTrue(js["embedded_presets_match"], "presets inside docs/index.html differ from examples/*.json")
        names = [k for k in js if k.endswith(".json")]
        self.assertGreaterEqual(len(names), 10)
        for f in names:
            with open(os.path.join(ROOT, "examples", f), encoding="utf-8") as fh:
                py = compute(json.load(fh))
            j = js[f]
            self.assertEqual(len(py.rows), len(j["rows"]), f)
            for a, b in zip(py.rows, j["rows"]):
                self.assertAlmostEqual(a["x"], b[0], places=6, msg=f)
                self.assertAlmostEqual(a["y"], b[1], delta=1e-4, msg=f)
                self.assertAlmostEqual(a["fr"], b[2], delta=1e-3, msg=f)
            self.assertEqual([x["status"] for x in py.jumps], [x["status"] for x in j["jumps"]], f)
            for a, b in zip(py.jumps, j["jumps"]):
                for k in ("x", "x_end", "y1", "y2", "y_end", "loss"):
                    if k in a:
                        self.assertAlmostEqual(a[k], b[k], delta=1e-3, msg=f"{f} {k}")
            self.assertEqual([r["types"] for r in py.reaches], j["types"], f)

    def test_alternate_and_conjugate_depths(self):
        from gvf import hydraulics as hy
        from gvf.sections import make_section
        out = subprocess.run(["node", os.path.join(ROOT, "tests", "js_check.js")], capture_output=True, text=True, check=True)
        defs = [{"type": "rectangular", "width": 5}, {"type": "trapezoidal", "width": 4, "side_slope": 2},
                {"type": "trapezoidal", "width": 0, "side_slope": 1.5}, {"type": "circular", "diameter": 2}]
        for r in json.loads(out.stdout)["alt"]:
            self.assertLess(r["dE"], 1e-9)   # same specific energy
            self.assertLess(r["dM"], 1e-9)   # same momentum function
            if r["alt"] is not None:
                sec = make_section(defs[r["si"]])
                py = hy.alternate_depth(sec, 3.0, r["y"], r["yc"])
                self.assertAlmostEqual(py, r["alt"], delta=1e-6)


if __name__ == "__main__":
    unittest.main()
