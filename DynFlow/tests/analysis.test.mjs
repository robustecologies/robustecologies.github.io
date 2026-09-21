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

process.exit(R.done().fails ? 1 : 0);
