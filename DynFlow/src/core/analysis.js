// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Local analysis: Jacobians by central differences, eigenvalues of a real
   matrix, equilibria (flows) and fixed points (maps) by damped Newton
   iteration from many seeds, and their stability. */
(function (DF) {
  "use strict";

  // G(x) = f(x) for flows and F(x) - x for maps, so both reduce to G = 0.
  function residual(sys, t, x, p, out) {
    sys.f(t, x, p, out, function (i) { return x[i]; });
    if (sys.time === "discrete") for (let i = 0; i < x.length; i++) out[i] -= x[i];
    return out;
  }

  // Jacobian of f (flows) or of F (maps), row-major, central differences.
  function jacobian(sys, t, x, p) {
    const n = x.length, J = new Float64Array(n * n), a = new Float64Array(n), b = new Float64Array(n), y = Float64Array.from(x);
    const H = function (i) { return y[i]; };
    for (let j = 0; j < n; j++) {
      const hj = 1e-6 * Math.max(1, Math.abs(x[j]));
      y[j] = x[j] + hj; sys.f(t, y, p, a, H);
      y[j] = x[j] - hj; sys.f(t, y, p, b, H);
      y[j] = x[j];
      for (let i = 0; i < n; i++) J[i * n + j] = (a[i] - b[i]) / (2 * hj);
    }
    return J;
  }

  // Solve A z = r in place (Gaussian elimination, partial pivoting); false if singular.
  function solve(A, r, n) {
    for (let c = 0; c < n; c++) {
      let piv = c;
      for (let i = c + 1; i < n; i++) if (Math.abs(A[i * n + c]) > Math.abs(A[piv * n + c])) piv = i;
      if (Math.abs(A[piv * n + c]) < 1e-300) return false;
      if (piv !== c) {
        for (let j = 0; j < n; j++) { const tmp = A[c * n + j]; A[c * n + j] = A[piv * n + j]; A[piv * n + j] = tmp; }
        const tr = r[c]; r[c] = r[piv]; r[piv] = tr;
      }
      for (let i = c + 1; i < n; i++) {
        const m = A[i * n + c] / A[c * n + c];
        if (m === 0) continue;
        for (let j = c; j < n; j++) A[i * n + j] -= m * A[c * n + j];
        r[i] -= m * r[c];
      }
    }
    for (let i = n - 1; i >= 0; i--) {
      let s = r[i];
      for (let j = i + 1; j < n; j++) s -= A[i * n + j] * r[j];
      r[i] = s / A[i * n + i];
    }
    return true;
  }

  /* Eigenvalues of a real n x n matrix (row-major): reduction to upper
     Hessenberg form by Gaussian similarity transforms, then the shifted QR
     iteration of the EISPACK routine hqr (Wilkinson and Reinsch, Handbook
     for Automatic Computation, vol. 2, 1971). Returns [{re, im}]. */
  function eigenvalues(M, n) {
    const a = [];
    for (let i = 0; i < n; i++) { a.push([]); for (let j = 0; j < n; j++) a[i].push(M[i * n + j]); }
    // elmhes
    for (let m = 1; m < n - 1; m++) {
      let x = 0, i = m;
      for (let j = m; j < n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
      if (i !== m) {
        for (let j = m - 1; j < n; j++) { const t = a[i][j]; a[i][j] = a[m][j]; a[m][j] = t; }
        for (let j = 0; j < n; j++) { const t = a[j][i]; a[j][i] = a[j][m]; a[j][m] = t; }
      }
      if (x !== 0) {
        for (i = m + 1; i < n; i++) {
          let y = a[i][m - 1];
          if (y !== 0) {
            y /= x; a[i][m - 1] = y;
            for (let j = m; j < n; j++) a[i][j] -= y * a[m][j];
            for (let j = 0; j < n; j++) a[j][m] += y * a[j][i];
          }
        }
      }
    }
    for (let i = 2; i < n; i++) for (let j = 0; j < i - 1; j++) a[i][j] = 0;
    // hqr
    const wr = new Array(n).fill(0), wi = new Array(n).fill(0);
    let anorm = 0;
    for (let i = 0; i < n; i++) for (let j = Math.max(i - 1, 0); j < n; j++) anorm += Math.abs(a[i][j]);
    let nn = n - 1, t = 0;
    while (nn >= 0) {
      let its = 0, l;
      do {
        for (l = nn; l >= 1; l--) {
          const s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
          if (Math.abs(a[l][l - 1]) + (s === 0 ? anorm : s) === (s === 0 ? anorm : s)) { a[l][l - 1] = 0; break; }
        }
        const x = a[nn][nn];
        if (l === nn) { wr[nn] = x + t; wi[nn--] = 0; }
        else {
          const y = a[nn - 1][nn - 1], w = a[nn][nn - 1] * a[nn - 1][nn];
          if (l === nn - 1) {
            const p = 0.5 * (y - x), q = p * p + w, z = Math.sqrt(Math.abs(q));
            const xx = x + t;
            if (q >= 0) {
              const zz = p + (p >= 0 ? Math.abs(z) : -Math.abs(z));
              wr[nn - 1] = wr[nn] = xx + zz;
              if (zz) wr[nn] = xx - w / zz;
              wi[nn - 1] = wi[nn] = 0;
            } else {
              wr[nn - 1] = wr[nn] = xx + p;
              wi[nn - 1] = -(wi[nn] = z);
            }
            nn -= 2;
          } else {
            if (its === 60) throw new Error("Eigenvalue iteration did not converge");
            let xs = x, ys = y, ws = w;
            if (its === 10 || its === 20) {
              t += xs;
              for (let i = 0; i <= nn; i++) a[i][i] -= xs;
              const s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
              ys = xs = 0.75 * s; ws = -0.4375 * s * s;
            }
            ++its;
            let m, p, q, r, z;
            for (m = nn - 2; m >= l; m--) {
              z = a[m][m];
              r = xs - z; const s0 = ys - z;
              p = (r * s0 - ws) / a[m + 1][m] + a[m][m + 1];
              q = a[m + 1][m + 1] - z - r - s0;
              r = a[m + 2][m + 1];
              const s = Math.abs(p) + Math.abs(q) + Math.abs(r);
              p /= s; q /= s; r /= s;
              if (m === l) break;
              const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
              const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
              if (u + v === v) break;
            }
            for (let i = m + 2; i <= nn; i++) { a[i][i - 2] = 0; if (i !== m + 2) a[i][i - 3] = 0; }
            for (let k = m; k <= nn - 1; k++) {
              if (k !== m) {
                p = a[k][k - 1]; q = a[k + 1][k - 1]; r = 0;
                if (k !== nn - 1) r = a[k + 2][k - 1];
                xs = Math.abs(p) + Math.abs(q) + Math.abs(r);
                if (xs !== 0) { p /= xs; q /= xs; r /= xs; }
              }
              const s = (p >= 0 ? 1 : -1) * Math.sqrt(p * p + q * q + r * r);
              if (s !== 0) {
                if (k === m) { if (l !== m) a[k][k - 1] = -a[k][k - 1]; } else a[k][k - 1] = -s * xs;
                p += s; xs = p / s; ys = q / s; z = r / s; q /= p; r /= p;
                for (let j = k; j <= nn; j++) {
                  p = a[k][j] + q * a[k + 1][j];
                  if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; }
                  a[k + 1][j] -= p * ys; a[k][j] -= p * xs;
                }
                const mmin = nn < k + 3 ? nn : k + 3;
                for (let i = l; i <= mmin; i++) {
                  p = xs * a[i][k] + ys * a[i][k + 1];
                  if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; }
                  a[i][k + 1] -= p * q; a[i][k] -= p;
                }
              }
            }
          }
        }
      } while (l < nn - 1);
    }
    const out = [];
    for (let i = 0; i < n; i++) out.push({ re: wr[i], im: wi[i] });
    return out;
  }

  /* Equilibria (flows) or fixed points (maps) inside a box. Seeds are a
     Halton sequence over the box plus any given points; each is refined by
     damped Newton iteration; solutions closer than 1e-6 of the box size
     are merged. Stability: all Re(lambda) < 0 for flows, all |lambda| < 1
     for maps; the tolerance on the boundary is 1e-7. */
  function findEquilibria(sys, p, box, opts) {
    opts = opts || {};
    const n = sys.vars.length, t = opts.t || 0, nSeeds = opts.seeds || 60;
    const scale = box.map(function (r) { return Math.max(Math.abs(r[1] - r[0]), 1e-12); });
    const seeds = (opts.extra || []).slice();
    const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    for (let s = 1; s <= nSeeds; s++) {
      seeds.push(box.map(function (r, i) {
        let f = 1, v = 0, k = s;
        const b = primes[i % primes.length];
        while (k > 0) { f /= b; v += f * (k % b); k = Math.floor(k / b); }
        return r[0] + v * (r[1] - r[0]);
      }));
    }
    const found = [];
    const G = new Float64Array(n), G2 = new Float64Array(n), x = new Float64Array(n), y = new Float64Array(n);
    const norm = function (v) { let s = 0; for (let i = 0; i < n; i++) s += (v[i] / scale[i]) * (v[i] / scale[i]); return Math.sqrt(s); };
    seeds.forEach(function (s0) {
      x.set(s0);
      let ok = false;
      for (let it = 0; it < 60; it++) {
        residual(sys, t, x, p, G);
        if (!G.every(isFinite)) return;
        const r0 = norm(G);
        if (r0 < 1e-11) { ok = true; break; }
        const J = jacobian(sys, t, x, p);
        if (sys.time === "discrete") for (let i = 0; i < n; i++) J[i * n + i] -= 1;
        const d = Float64Array.from(G);
        if (!solve(J, d, n)) return;
        let lam = 1;
        for (let ls = 0; ls < 20; ls++) {
          for (let i = 0; i < n; i++) y[i] = x[i] - lam * d[i];
          residual(sys, t, y, p, G2);
          if (G2.every(isFinite) && norm(G2) < (1 - 1e-4 * lam) * r0) break;
          lam *= 0.5;
        }
        x.set(y);
        if (norm(d) * lam < 1e-13) { residual(sys, t, x, p, G); ok = norm(G) < 1e-8; break; }
      }
      if (!ok) return;
      for (let i = 0; i < n; i++) {
        const m = 0.5 * scale[i];
        if (x[i] < box[i][0] - m || x[i] > box[i][1] + m) return;
      }
      if (found.some(function (e) { let d = 0; for (let i = 0; i < n; i++) d = Math.max(d, Math.abs(e.x[i] - x[i]) / scale[i]); return d < 1e-6; })) return;
      found.push({ x: Array.from(x) });
    });
    found.forEach(function (e) {
      const J = jacobian(sys, t, Float64Array.from(e.x), p);
      e.eig = eigenvalues(J, n);
      if (sys.time === "discrete") {
        const rho = Math.max.apply(null, e.eig.map(function (l) { return Math.hypot(l.re, l.im); }));
        e.stable = rho < 1 - 1e-7;
        e.type = rho < 1 - 1e-7 ? "stable" : rho > 1 + 1e-7 ? (e.eig.some(function (l) { return Math.hypot(l.re, l.im) < 1 - 1e-7; }) ? "saddle" : "unstable") : "marginal";
      } else {
        const re = e.eig.map(function (l) { return l.re; });
        const mx = Math.max.apply(null, re), mn = Math.min.apply(null, re);
        const osc = e.eig.some(function (l) { return Math.abs(l.im) > 1e-9; });
        e.stable = mx < -1e-7;
        if (mx < -1e-7) e.type = osc ? "stable focus" : "stable node";
        else if (mn > 1e-7) e.type = osc ? "unstable focus" : "unstable node";
        else if (mn < -1e-7 && mx > 1e-7) e.type = osc ? "saddle focus" : "saddle";
        else e.type = osc ? "centre" : "non-hyperbolic";
      }
    });
    return found;
  }

  DF.residual = residual;
  DF.jacobian = jacobian;
  DF.solveLinear = solve;
  DF.eigenvalues = eigenvalues;
  DF.findEquilibria = findEquilibria;
})(globalThis.DynFlow = globalThis.DynFlow || {});
