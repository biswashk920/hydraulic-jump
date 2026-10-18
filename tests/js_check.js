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
  out[f] = {rows: r.rows.map(w => [w.x, w.y, w.fr]), jumps: r.jumps.map(j => ({status: j.status, x: j.x, x_end: j.x_end, y1: j.y1, y2: j.y2, y_end: j.y_end, loss: j.loss})),
            types: r.reaches.map(q => q.types), nwarn: r.warnings.length};
}
const G = g.GVF, secDefs = [{type:"rectangular",width:5},{type:"trapezoidal",width:4,side_slope:2},{type:"trapezoidal",width:0,side_slope:1.5},{type:"circular",diameter:2}];
out.alt = [];
secDefs.forEach((d, si) => {
  const s = G.mkSec(d), Q = 3.0;
  const yc = (() => { let lo = 1e-6, hi = 50; for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2, A = s.area(m); (Q*Q*s.top(m)/(9.81*A*A*A) > 1 ? lo = m : hi = m); } return (lo + hi) / 2; })();
  [0.4, 0.7, 1.3, 2.0].forEach(k => { const y = k * yc, a = G.alt(s, Q, y, yc), c = G.conj(s, Q, y, yc);
    out.alt.push({si, y, yc, alt: a, conj: c, dE: a === null ? 0 : Math.abs(G.en(s, Q, a) - G.en(s, Q, y)), dM: c === null ? 0 : Math.abs(G.mo(s, Q, c) - G.mo(s, Q, y))}); });
});
console.log(JSON.stringify(out));
