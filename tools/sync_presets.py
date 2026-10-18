"""Copies the example channels into the presets list inside docs/index.html.
Run from the project root after changing an example:   python tools/sync_presets.py
(tests/test_js.py fails if the two are out of step.)"""
import json
import os
import re

ORDER = ["gate_jump_uniform", "jump_gate_mild", "steep_to_mild_jump", "mild_to_steep", "backwater_weir",
         "drowned_jump", "swept_jump", "trapezoid_changing", "circular_culvert", "horizontal_adverse"]
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
presets = []
for name in ORDER:
    with open(os.path.join(root, "examples", name + ".json"), encoding="utf-8") as f:
        presets.append(json.load(f))
path = os.path.join(root, "docs", "index.html")
with open(path, encoding="utf-8") as f:
    html = f.read()
html, n = re.subn(r"window\.PRESETS=.*?;</script>", lambda m: "window.PRESETS=" + json.dumps(presets, separators=(",", ":")) + ";</script>", html, count=1)
assert n == 1, "presets line not found"
with open(path, "w", encoding="utf-8") as f:
    f.write(html)
print(f"{len(presets)} presets written to docs/index.html")
