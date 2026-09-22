// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Local analysis: Jacobians by central differences, eigenvalues of a real
   matrix, equilibria (flows) and fixed points (maps) by damped Newton
   iteration from many seeds, their stability, and branches of equilibria
   against one parameter by pseudo-arclength continuation. */
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

  /* Branches of equilibria (flows) or fixed points (maps) against parameter
     pi on [from, to], by pseudo-arclength continuation (Allgower and Georg
     2003, ch. 2 and 8). The unknowns are u = (x, p) in coordinates scaled
     by the axis box and the parameter interval; a step predicts along the
     tangent v of G(u) = 0 and corrects by Newton iteration on
       G(x, p) = 0,   v . (u - u_pred) = 0,
     so the branch turns at a fold instead of ending there. Seeds are the
     equilibria found at seven parameter values; a seed that lies on a branch
     already traced is skipped. Returns
       { branches: [[{ x, p, stable, eig }]], points: [{ kind, p, x }] }
     where kind is "fold" (p extremal along the branch, where the number of
     unstable directions changes, located by the vertex of a parabola
     through three points),
     "hopf" (a complex pair crosses the imaginary axis), "torus" (a complex
     pair of a map crosses the unit circle), "flip" (an eigenvalue of a map
     crosses -1, period doubling) or "branch" (a real eigenvalue crosses
     away from a fold, or the branch turns without a change of stability,
     as at a transcritical or pitchfork point). */
  function continueBranches(sys, base, pi, from, to, box, opts) {
    opts = opts || {};
    const n = sys.vars.length, t = opts.t || 0;
    const sx = box.map(function (r) { return Math.max(Math.abs(r[1] - r[0]), 1e-12); }), sp = (to - from) || 1;
    const dsMax = opts.dsMax || 0.02, maxPts = opts.maxPoints || 1500, discrete = sys.time === "discrete";
    const pp = Float64Array.from(base), g = new Float64Array(n), ga = new Float64Array(n), gb = new Float64Array(n), xv = new Float64Array(n);
    const G = function (x, p, out) { pp[pi] = p; xv.set(x); return residual(sys, t, xv, pp, out); };
    // Scaled augmented Jacobian [G_x S_x | G_p s_p] (n rows, n + 1 columns) and the unscaled G_x.
    const jac = function (x, p) {
      pp[pi] = p; xv.set(x);
      const Jx = jacobian(sys, t, xv, pp);
      if (discrete) for (let i = 0; i < n; i++) Jx[i * n + i] -= 1;
      const hp = 1e-6 * Math.max(1, Math.abs(p));
      G(x, p + hp, ga); G(x, p - hp, gb);
      const A = new Float64Array(n * (n + 1));
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) A[i * (n + 1) + j] = Jx[i * n + j] * sx[j];
        A[i * (n + 1) + n] = (ga[i] - gb[i]) / (2 * hp) * sp;
      }
      if (discrete) for (let i = 0; i < n; i++) Jx[i * n + i] += 1;
      return { A: A, Jx: Jx };
    };
    // Solve the bordered system [A; w^T] z = rhs of order n + 1.
    const bordered = function (A, w, rhs) {
      const m = n + 1, M = new Float64Array(m * m), z = Float64Array.from(rhs);
      for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) M[i * m + j] = A[i * m + j];
      for (let j = 0; j < m; j++) M[n * m + j] = w[j];
      return solve(M, z, m) ? z : null;
    };
    const tangent = function (x, p, prev) {
      const A = jac(x, p).A, w = prev || new Float64Array(n + 1).fill(0);
      if (!prev) w[n] = 1;
      const rhs = new Float64Array(n + 1); rhs[n] = 1;
      let v = bordered(A, w, rhs);
      if (!v && !prev) { w[n] = 0; w[0] = 1; v = bordered(A, w, rhs); }
      if (!v) return null;
      let nv = 0; for (let j = 0; j <= n; j++) nv += v[j] * v[j];
      nv = Math.sqrt(nv);
      if (!(nv > 0) || !isFinite(nv)) return null;
      let dot = 0; for (let j = 0; j <= n; j++) { v[j] /= nv; if (prev) dot += v[j] * prev[j]; }
      if (prev && dot < 0) for (let j = 0; j <= n; j++) v[j] = -v[j];
      return v;
    };
    const stability = function (x, p) {
      const eig = eigenvalues(jac(x, p).Jx, n);
      const stable = discrete ? eig.every(function (l) { return Math.hypot(l.re, l.im) < 1 - 1e-7; }) : eig.every(function (l) { return l.re < -1e-7; });
      return { stable: stable, eig: eig };
    };
    const inBox = function (x, p) {
      if (p < from - 0.02 * sp || p > to + 0.02 * sp) return false;
      for (let j = 0; j < n; j++) if (x[j] < box[j][0] - 0.5 * sx[j] || x[j] > box[j][1] + 0.5 * sx[j]) return false;
      return true;
    };
    function run(x0, p0, dir, budget) {
      let x = Float64Array.from(x0), p = p0, v = tangent(x, p, null);
      if (!v) return [];
      if (dir < 0) for (let j = 0; j <= n; j++) v[j] = -v[j];
      const st0 = stability(x, p), out = [{ x: Array.from(x), p: p, stable: st0.stable, eig: st0.eig, tp: v[n] }];
      let ds = dsMax / 2;
      const xn = new Float64Array(n), xp = new Float64Array(n), rhs = new Float64Array(n + 1);
      for (let k = 0; k < budget; k++) {
        let ok = false, pn = p;
        for (let tries = 0; tries < 10 && !ok; tries++) {
          for (let j = 0; j < n; j++) { xp[j] = x[j] + ds * v[j] * sx[j]; xn[j] = xp[j]; }
          const pp0 = p + ds * v[n] * sp; pn = pp0;
          for (let it = 0; it < 12; it++) {
            G(xn, pn, g);
            if (!g.every(isFinite)) break;
            const A = jac(xn, pn).A;
            let arc = (pn - pp0) / sp * v[n];
            for (let j = 0; j < n; j++) arc += (xn[j] - xp[j]) / sx[j] * v[j];
            for (let i = 0; i < n; i++) rhs[i] = g[i];
            rhs[n] = arc;
            const d = bordered(A, v, rhs);
            if (!d) break;
            let nd = 0;
            for (let j = 0; j < n; j++) { xn[j] -= d[j] * sx[j]; nd += d[j] * d[j]; }
            pn -= d[n] * sp; nd = Math.sqrt(nd + d[n] * d[n]);
            if (nd < 1e-10) {
              G(xn, pn, g);
              let gm = 0; for (let i = 0; i < n; i++) gm = Math.max(gm, Math.abs(g[i]));
              ok = gm < 1e-8;
              if (ok) ds = it > 4 ? ds * 0.7 : it < 3 ? Math.min(ds * 1.4, dsMax) : ds;
              break;
            }
          }
          if (!ok) ds *= 0.5;
          if (ds < 1e-7) break;
        }
        if (!ok) break;
        const vn = tangent(xn, pn, v);
        if (!vn) break;
        x = Float64Array.from(xn); p = pn; v = vn;
        const st = stability(x, p);
        out.push({ x: Array.from(x), p: p, stable: st.stable, eig: st.eig, tp: v[n] });
        if (!inBox(x, p)) break;
        // A closed branch (an isola) returns to its first point.
        let back = 0; for (let j = 0; j < n; j++) back += Math.pow((x[j] - out[0].x[j]) / sx[j], 2);
        back = Math.sqrt(back + Math.pow((p - out[0].p) / sp, 2));
        if (out.length > 10 && back < 0.6 * ds) break;
      }
      return out;
    }
    const branches = [];
    let used = 0;
    const onBranch = function (x, p) {
      return branches.some(function (br) {
        return br.some(function (q) {
          let d = Math.pow((p - q.p) / sp, 2);
          for (let j = 0; j < n; j++) d += Math.pow((x[j] - q.x[j]) / sx[j], 2);
          return d < 0.03 * 0.03;
        });
      });
    };
    for (let c = 0; c <= 6 && used < maxPts; c++) {
      const pv = from + sp * c / 6;
      pp[pi] = pv;
      let eqs = [];
      try { eqs = findEquilibria(sys, Float64Array.from(pp), box, { t: t, seeds: opts.seeds || 24 }); } catch (e) { eqs = []; }
      for (let e = 0; e < eqs.length && used < maxPts; e++) {
        if (onBranch(eqs[e].x, pv)) continue;
        const budget = Math.max(20, (maxPts - used) >> 1);
        const back = run(eqs[e].x, pv, -1, budget).reverse();
        back.forEach(function (q) { q.tp = -q.tp; });
        const fwd = run(eqs[e].x, pv, 1, budget);
        const br = back.concat(fwd.slice(1));
        if (br.length > 1) { branches.push(br); used += br.length; }
      }
    }
    // Branches that end where another begins are one branch through a corner
    // of a nonsmooth field (the Stommel model with |T - S|, for instance).
    const gap = function (a, b) {
      let d = Math.pow((a.p - b.p) / sp, 2);
      for (let j = 0; j < n; j++) d += Math.pow((a.x[j] - b.x[j]) / sx[j], 2);
      return Math.sqrt(d);
    };
    for (let merged = true; merged;) {
      merged = false;
      for (let a = 0; a < branches.length && !merged; a++) for (let b = a + 1; b < branches.length && !merged; b++) {
        const A = branches[a], B = branches[b], tol = 1e-4;
        let joined = null;
        if (gap(A[A.length - 1], B[0]) < tol) joined = A.concat(B.slice(1));
        else if (gap(A[A.length - 1], B[B.length - 1]) < tol) joined = A.concat(B.slice(0, -1).reverse());
        else if (gap(A[0], B[B.length - 1]) < tol) joined = B.concat(A.slice(1));
        else if (gap(A[0], B[0]) < tol) joined = B.slice().reverse().concat(A.slice(1));
        if (joined) { branches.splice(b, 1); branches[a] = joined; merged = true; }
      }
    }
    // Special points along each branch. The critical quantity m is the
    // largest real part of an eigenvalue (flows) or the spectral radius minus
    // one (maps); a stability change is located where m, interpolated
    // linearly between two points, vanishes.
    const crit = function (q) {
      return discrete ? Math.max.apply(null, q.eig.map(function (l) { return Math.hypot(l.re, l.im); })) - 1 : Math.max.apply(null, q.eig.map(function (l) { return l.re; }));
    };
    const unstableDim = function (q) {
      return q.eig.filter(function (l) { return discrete ? Math.hypot(l.re, l.im) > 1 + 1e-7 : l.re > 1e-7; }).length;
    };
    const points = [];
    branches.forEach(function (br) {
      const s = [0];
      for (let i = 1; i < br.length; i++) {
        let d = Math.pow((br[i].p - br[i - 1].p) / sp, 2);
        for (let j = 0; j < n; j++) d += Math.pow((br[i].x[j] - br[i - 1].x[j]) / sx[j], 2);
        s.push(s[i - 1] + Math.sqrt(d));
      }
      const turns = [];
      for (let i = 1; i < br.length - 1; i++) {
        const a = br[i - 1].p, b = br[i].p, c = br[i + 1].p;
        if (!((b - a) * (c - b) < 0)) continue;
        const s0 = s[i - 1], s1 = s[i], s2 = s[i + 1];
        const d1 = (b - a) / (s1 - s0), d2 = (c - b) / (s2 - s1), k2 = (d2 - d1) / (s2 - s0);
        const sv = k2 !== 0 ? (s0 + s1) / 2 - d1 / (2 * k2) : s1;
        const w = Math.min(1, Math.max(0, (sv - s0) / (s2 - s0)));
        turns.push(i);
        // A turning point where one eigenvalue crosses (the number of unstable
        // directions changes by one) is a fold; with the same number on both
        // sides the branch meets another one, as at a pitchfork.
        const kind = unstableDim(br[i - 1]) !== unstableDim(br[i + 1]) ? "fold" : "branch";
        points.push({ kind: kind, p: a + d1 * (sv - s0) + k2 * (sv - s0) * (sv - s1), x: br[i - 1].x.map(function (v, j) { return v + w * (br[i + 1].x[j] - v); }) });
      }
      for (let i = 1; i < br.length; i++) {
        if (br[i].stable === br[i - 1].stable) continue;
        if (turns.some(function (f) { return Math.abs(f - i) <= 1; })) continue;
        const m0 = crit(br[i - 1]), m1 = crit(br[i]), w = m1 !== m0 ? Math.min(1, Math.max(0, -m0 / (m1 - m0))) : 0.5;
        // The eigenvalue that crosses: complex (Hopf, or Neimark-Sacker for a
        // map), real and negative for a map (period doubling), otherwise real.
        const q = w < 0.5 ? br[i - 1] : br[i];
        let lead = q.eig[0];
        q.eig.forEach(function (l) { if ((discrete ? Math.hypot(l.re, l.im) : l.re) > (discrete ? Math.hypot(lead.re, lead.im) : lead.re)) lead = l; });
        const complex = Math.abs(lead.im) > 1e-9;
        const kind = complex ? (discrete ? "torus" : "hopf") : discrete && lead.re < 0 ? "flip" : "branch";
        points.push({ kind: kind, p: br[i - 1].p + w * (br[i].p - br[i - 1].p), x: br[i - 1].x.map(function (v, j) { return v + w * (br[i].x[j] - v); }) });
      }
    });
    // One point per bifurcation: a pitchfork is met by two branches.
    const unique = points.filter(function (q, i) {
      return !points.some(function (r, j) { return j < i && r.kind === q.kind && gap(r, q) < 0.01; });
    });
    return { branches: branches, points: unique };
  }

  DF.continueBranches = continueBranches;
  DF.residual = residual;
  DF.jacobian = jacobian;
  DF.solveLinear = solve;
  DF.eigenvalues = eigenvalues;
  DF.findEquilibria = findEquilibria;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
