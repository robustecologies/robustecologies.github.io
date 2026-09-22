// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Eigenvalues and equilibria against invariants and closed forms.
import { load, reporter, CORE } from "./harness.mjs";

const DF = load(CORE.concat(["src/core/analysis.js"]));
const R = reporter("analysis");

// Complex determinant of A - lambda I, by elimination with partial pivoting.
function charDet(A, n, lr, li) {
  const re = [], im = [];
  for (let i = 0; i < n; i++) { re.push([]); im.push([]); for (let j = 0; j < n; j++) { re[i].push(A[i * n + j] - (i === j ? lr : 0)); im[i].push(i === j ? -li : 0); } }
  let dr = 1, di = 0;
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let i = c + 1; i < n; i++) if (Math.hypot(re[i][c], im[i][c]) > Math.hypot(re[piv][c], im[piv][c])) piv = i;
    if (piv !== c) { [re[c], re[piv]] = [re[piv], re[c]]; [im[c], im[piv]] = [im[piv], im[c]]; dr = -dr; di = -di; }
    const pr = re[c][c], pi = im[c][c];
    const nr = dr * pr - di * pi, ni = dr * pi + di * pr; dr = nr; di = ni;
    const d2 = pr * pr + pi * pi;
    if (d2 === 0) return 0;
    for (let i = c + 1; i < n; i++) {
      const mr = (re[i][c] * pr + im[i][c] * pi) / d2, mi = (im[i][c] * pr - re[i][c] * pi) / d2;
      for (let j = c; j < n; j++) { re[i][j] -= mr * re[c][j] - mi * im[c][j]; im[i][j] -= mr * im[c][j] + mi * re[c][j]; }
    }
  }
  return Math.hypot(dr, di);
}

// Random matrices: trace, determinant and det(A - lambda I) = 0 for every eigenvalue.
{
  const rng = new DF.RNG(99);
  let worstTrace = 0, worstChar = 0, cases = 0;
  for (const n of [2, 3, 4, 5, 6, 8, 10]) {
    for (let rep = 0; rep < 40; rep++) {
      const A = Float64Array.from({ length: n * n }, () => rng.normal());
      const ev = DF.eigenvalues(A, n);
      let tr = 0; for (let i = 0; i < n; i++) tr += A[i * n + i];
      const sre = ev.reduce((s, l) => s + l.re, 0), sim = ev.reduce((s, l) => s + l.im, 0);
      worstTrace = Math.max(worstTrace, Math.abs(sre - tr) / (1 + Math.abs(tr)), Math.abs(sim));
      // relative characteristic residual |det(A - lambda I)| / prod(1 + |lambda_j - lambda|)
      for (const l of ev) {
        const scale = ev.reduce((p, m) => p * (1 + Math.hypot(m.re - l.re, m.im - l.im)), 1);
        worstChar = Math.max(worstChar, charDet(A, n, l.re, l.im) / scale);
      }
      cases++;
    }
  }
  R.check("eigenvalues of random matrices", worstTrace < 1e-10 && worstChar < 1e-8,
    `${cases} Gaussian matrices of order 2 to 10: worst trace mismatch ${worstTrace.toExponential(1)}, worst |det(A - lambda I)| relative ${worstChar.toExponential(1)}`);
}
{
  const cases = [
    [[0, 1, -2, -3], 2, [[-1, 0], [-2, 0]]],
    [[0, -1, 1, 0], 2, [[0, 1], [0, -1]]],
    [[2, 0, 0, 0, 3, 4, 0, 4, 9], 3, [[1, 0], [2, 0], [11, 0]]],
    [[0, 1, 0, 0, 0, 1, 6, -11, 6], 3, [[1, 0], [2, 0], [3, 0]]]
  ];
  const bad = cases.filter(([A, n, want]) => {
    const got = DF.eigenvalues(Float64Array.from(A), n).map((l) => [l.re, l.im]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const w = want.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    return got.some((g, i) => Math.abs(g[0] - w[i][0]) > 1e-10 || Math.abs(g[1] - w[i][1]) > 1e-10);
  });
  R.check("eigenvalues of known matrices", bad.length === 0, `${cases.length - bad.length} of ${cases.length} (companion of s^2 + 3s + 2, rotation, symmetric block, companion of (s-1)(s-2)(s-3))`);
}

// Equilibria.
{
  const lor = DF.compileSystem("x' = s*(y-x)\ny' = x*(r-z) - y\nz' = x*y - b*z\nparam s = 10\nparam r = 28\nparam b = 2.6666666666666665");
  const eq = DF.findEquilibria(lor, [10, 28, 8 / 3], [[-25, 25], [-30, 30], [0, 50]]);
  const c = Math.sqrt(8 / 3 * 27);
  const ok = eq.length === 3 && eq.every((e) => (Math.abs(e.x[0]) < 1e-8 && e.type === "saddle") || (Math.abs(Math.abs(e.x[0]) - c) < 1e-8 && Math.abs(e.x[2] - 27) < 1e-8 && e.type === "saddle focus"));
  R.check("Lorenz equilibria", ok, `${eq.length} found: ${eq.map((e) => "(" + e.x.map((v) => v.toFixed(4)).join(", ") + ") " + e.type).join("; ")}; C = +/-${c.toFixed(4)}, z = 27`);
  const bru = DF.compileSystem("X' = A - (B+1)*X + X^2*Y\nY' = B*X - X^2*Y\nparam A = 1\nparam B = 3");
  const below = DF.findEquilibria(bru, [1, 1.9], [[0, 4], [0, 6]]), above = DF.findEquilibria(bru, [1, 2.1], [[0, 4], [0, 6]]);
  R.check("Brusselator Hopf at B = 1 + A^2", below.length === 1 && below[0].stable && above.length === 1 && !above[0].stable && Math.abs(below[0].x[1] - 1.9) < 1e-9,
    `B = 1.9: ${below[0].type} at (${below[0].x.map((v) => v.toFixed(4))}); B = 2.1: ${above[0].type}`);
  const logi = DF.compileSystem("x[n+1] = r*x*(1-x)\nparam r = 2.8");
  const f1 = DF.findEquilibria(logi, [2.8], [[0, 1]]), f2 = DF.findEquilibria(logi, [3.2], [[0, 1]]);
  const star = f1.find((e) => e.x[0] > 0.1);
  R.check("logistic fixed points", f1.length === 2 && star && Math.abs(star.x[0] - (1 - 1 / 2.8)) < 1e-10 && star.stable && !f2.find((e) => e.x[0] > 0.1).stable,
    `r = 2.8: x* = ${star.x[0].toFixed(10)} (1 - 1/r = ${(1 - 1 / 2.8).toFixed(10)}) ${star.type}; r = 3.2: ${f2.find((e) => e.x[0] > 0.1).type}`);
  const rm = DF.compileSystem("N' = r*N*(1 - N/K) - a*N*P/(1 + a*h*N)\nP' = e*a*N*P/(1 + a*h*N) - m*P\nparam r = 1\nparam K = 3.7\nparam a = 1\nparam h = 0.4\nparam e = 0.6\nparam m = 0.3");
  const s1 = DF.findEquilibria(rm, [1, 3.7, 1, 0.4, 0.6, 0.3], [[0, 4], [0, 4]]).find((e) => e.x[0] > 0 && e.x[1] > 0);
  const s2 = DF.findEquilibria(rm, [1, 3.8, 1, 0.4, 0.6, 0.3], [[0, 4], [0, 4]]).find((e) => e.x[0] > 0 && e.x[1] > 0);
  R.check("Rosenzweig-MacArthur Hopf at K = 3.75", s1.stable && !s2.stable && Math.abs(s1.x[0] - 0.625) < 1e-9,
    `K = 3.7: ${s1.type} at N* = ${s1.x[0].toFixed(6)} (0.625); K = 3.8: ${s2.type}`);
}

// Branches of equilibria by pseudo-arclength continuation, against closed forms.
{
  const fold = (sys, pars, pi, a, b, box) => DF.continueBranches(DF.compileSystem(sys), Float64Array.from(pars), pi, a, b, box);
  // Cusp x' = r + a x - x^3 at a = 1: folds at r = +/- 2 (1/3)^(3/2), one S-shaped branch.
  const cusp = fold("x' = r + a*x - x^3\nparam r = 0\nparam a = 1", [0, 1], 0, -0.8, 0.8, [[-1.6, 1.6]]);
  const exact = 2 * Math.pow(1 / 3, 1.5), cf = cusp.points.filter((q) => q.kind === "fold").map((q) => q.p).sort((u, v) => u - v);
  const resid = Math.max.apply(null, cusp.branches[0].map((q) => Math.abs(q.p + q.x[0] - Math.pow(q.x[0], 3))));
  R.check("continuation through the folds of the cusp", cusp.branches.length === 1 && cf.length === 2 && Math.abs(cf[0] + exact) < 1e-3 && Math.abs(cf[1] - exact) < 1e-3 && resid < 1e-8,
    `1 branch; folds at r = ${cf.map((v) => v.toFixed(5)).join(", ")} (exact -/+${exact.toFixed(5)}); largest residual ${resid.toExponential(1)}`);
  // Pitchfork x' = a x - x^3: x = 0 stable for a < 0 only, x = +/- sqrt(a) stable, one branch point at a = 0.
  const pf = fold("x' = a*x - x^3\nparam a = 1", [1], 0, -1, 1.5, [[-1.4, 1.4]]);
  const non = pf.branches.flat().filter((q) => Math.abs(q.x[0]) > 1e-3), triv = pf.branches.flat().filter((q) => Math.abs(q.x[0]) < 1e-9 && Math.abs(q.p) > 0.05);
  const bp = pf.points.filter((q) => q.kind === "branch");
  R.check("pitchfork branches and stability", Math.max.apply(null, non.map((q) => Math.abs(q.x[0] * q.x[0] - q.p))) < 1e-9 && non.every((q) => q.stable) && triv.every((q) => q.stable === q.p < 0) && bp.length === 1 && Math.abs(bp[0].p) < 1e-3,
    `${non.length} points with x^2 = a, all stable; trivial branch stable exactly for a < 0; one branch point at a = ${bp.length ? bp[0].p.toExponential(1) : "none"}`);
  // Brusselator: Hopf point B = 1 + A^2 = 2 at A = 1.
  const bru = fold("X' = A - (B + 1)*X + X^2*Y\nY' = B*X - X^2*Y\nparam A = 1\nparam B = 3", [1, 3], 1, 0.3, 6, [[0, 4.5], [0, 5.5]]);
  const hopf = bru.points.filter((q) => q.kind === "hopf");
  R.check("Hopf point of the Brusselator", hopf.length === 1 && Math.abs(hopf[0].p - 2) < 1e-3, `Hopf at B = ${hopf.length ? hopf[0].p.toFixed(5) : "none"} (1 + A^2 = 2)`);
  // Logistic map: period doubling of the fixed point 1 - 1/r at r = 3.
  const lg = fold("x[n+1] = r*x*(1 - x)\nparam r = 3", [3], 0, 2.5, 4, [[0, 1]]);
  const flip = lg.points.filter((q) => q.kind === "flip");
  R.check("period doubling of the logistic map", flip.length === 1 && Math.abs(flip[0].p - 3) < 1e-3, `flip at r = ${flip.length ? flip[0].p.toFixed(5) : "none"} (3)`);
  // May grazing model: folds at the extrema of c(V) = r (1 - V/K) (V^2 + V0^2) / V.
  const gr = fold("V' = r*V*(1 - V/K) - c*V^2/(V^2 + V0^2)\nparam r = 1\nparam K = 10\nparam c = 2\nparam V0 = 1", [1, 10, 2, 1], 2, 1, 3, [[0, 10]]);
  const cV = (V) => (1 - V / 10) * (V * V + 1) / V, ext = [];
  for (let V = 0.05, prev = cV(0.04); V < 10; V += 1e-5) { const c = cV(V), nx = cV(V + 1e-5); if ((c - prev) * (nx - c) < 0) ext.push(c); prev = c; }
  const gf = gr.points.filter((q) => q.kind === "fold").map((q) => q.p).sort((u, v) => u - v), ge = ext.filter((c) => c > 1 && c < 3).sort((u, v) => u - v);
  R.check("folds of the May grazing model", gf.length === 2 && ge.length === 2 && Math.abs(gf[0] - ge[0]) < 2e-3 && Math.abs(gf[1] - ge[1]) < 2e-3, `folds at c = ${gf.map((v) => v.toFixed(5)).join(", ")}; extrema of c(V): ${ge.map((v) => v.toFixed(5)).join(", ")}`);
}

process.exit(R.done().fails ? 1 : 0);
