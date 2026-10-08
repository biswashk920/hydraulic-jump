// Runs the calculation core that is embedded in docs/index.html (the same code the browser runs)
// on every example JSON and prints the results as JSON. Used by tests/test_js.py.
const fs = require("fs"), path = require("path");
const root = path.join(__dirname, "..");
const html = fs.readFileSync(path.join(root, "docs", "index.html"), "utf8");
const core = html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const g = {};
new Function("globalThis", core.replace(/\}\)\(typeof window[^;]*;/, "})(globalThis);"))(g);
const presets = JSON.parse(html.match(/window\.PRESETS=(.*?);<\/script>/)[1]);
const out = {embedded_presets_match: true};
for (const f of fs.readdirSync(path.join(root, "examples")).filter(f => f.endsWith(".json") && !f.startsWith("textbook")).sort()) {
  const spec = JSON.parse(fs.readFileSync(path.join(root, "examples", f), "utf8"));
  if (!spec.reaches) continue;
  const inHtml = presets.find(p => JSON.stringify(p) === JSON.stringify(spec));
  if (!inHtml) out.embedded_presets_match = false;
  const r = g.GVF.compute(spec);
  out[f] = {rows: r.rows.map(w => [w.x, w.y, w.fr]), jumps: r.jumps.map(j => ({status: j.status, x: j.x, y1: j.y1, y2: j.y2, loss: j.loss})),
            types: r.reaches.map(q => q.types), nwarn: r.warnings.length};
}
console.log(JSON.stringify(out));
