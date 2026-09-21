// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Coordinates: automatic axis ranges, the 2D and rotating 3D projections
   from state space to the canvas, and axis drawing. */
(function (DF) {
  "use strict";

  /* Axis ranges from a probe run: several members from the initial
     condition, a transient, then the 1st and 99th percentiles of each
     variable with a 6 percent margin. Declared `range` lines win. */
  function autoRanges(sys, params, init, opts) {
    opts = opts || {};
    const dim = sys.vars.length, dt = opts.dt || 0.01, steps = opts.steps || 4000;
    const sim = new DF.Simulator(sys, { dt: dt, params: params, init: init, n: 8, spread: opts.spread || 0.1, seed: 12345, perturbations: opts.perturbations || [] });
    const samples = sys.vars.map(function () { return []; });
    const transient = Math.floor(steps / 4);
    for (let s = 0; s < steps; s++) {
      sim.step();
      if (s < transient || s % 2) continue;
      for (let k = 0; k < sim.n; k++) {
        if (!sim.alive[k]) continue;
        for (let i = 0; i < dim; i++) samples[i].push(sim.X[k * dim + i]);
      }
    }
    const out = {};
    sys.vars.forEach(function (v, i) {
      if (sys.ranges[v]) { out[v] = sys.ranges[v].slice(); return; }
      const a = samples[i].filter(isFinite).sort(function (u, w) { return u - w; });
      let lo, hi;
      if (!a.length) { lo = -1; hi = 1; }
      else { lo = a[Math.floor(0.01 * (a.length - 1))]; hi = a[Math.ceil(0.99 * (a.length - 1))]; }
      if (!(hi > lo)) { const c = isFinite(lo) ? lo : 0, w = Math.max(Math.abs(c) * 0.5, 0.5); lo = c - w; hi = c + w; }
      const m = 0.06 * (hi - lo);
      out[v] = [lo - m, hi + m];
    });
    return out;
  }

  /* Camera. axes holds two or three variable indices; ranges maps each
     axis to [lo, hi]. In 3D the box is centred, rotated by azimuth about the
     vertical axis and tilted by elevation, and projected orthographically. */
  function Camera(axes, ranges, opts) {
    opts = opts || {};
    this.axes = axes; this.ranges = ranges;
    this.azim = opts.azim === undefined ? 0.6 : opts.azim;
    this.elev = opts.elev === undefined ? 0.35 : opts.elev;
    this.zoom = opts.zoom || 1;
    this.pad = opts.pad === undefined ? 0.08 : opts.pad;
    this.upAxis = opts.upAxis === undefined ? 2 : opts.upAxis; // in 3D, which of the three is vertical
    this.simplex = !!opts.simplex && axes.length === 3;           // barycentric triangle for three shares
    this.w = 1; this.h = 1;
    this.inset = { l: 0, r: 0, t: 0, b: 0 };
  }
  Camera.prototype.resize = function (w, h, inset) { this.w = w; this.h = h; if (inset) this.inset = inset; this.tri = null; };
  Camera.prototype.is3D = function () { return this.axes.length === 3 && !this.simplex; };
  // Vertices of the simplex triangle: first axis on top, then lower left, lower right.
  Camera.prototype.triangle = function () {
    const b = this.plotBox(), side = Math.min(b.w * (1 - 2 * this.pad), b.h * (1 - 2 * this.pad) * 2 / Math.sqrt(3)) * this.zoom;
    const cx = b.x + b.w / 2, cy = b.y + b.h / 2 + side * Math.sqrt(3) / 12, hgt = side * Math.sqrt(3) / 2;
    return [[cx, cy - hgt * 2 / 3], [cx - side / 2, cy + hgt / 3], [cx + side / 2, cy + hgt / 3]];
  };
  Camera.prototype.plotBox = function () {
    const i = this.inset;
    return { x: i.l, y: i.t, w: this.w - i.l - i.r, h: this.h - i.t - i.b };
  };
  // Normalised coordinate in [-1, 1] of value v on axis a.
  Camera.prototype.norm = function (a, v) {
    const r = this.ranges[a];
    return 2 * (v - r[0]) / (r[1] - r[0]) - 1;
  };
  Camera.prototype.project = function (x, out) {
    out = out || [0, 0];
    const b = this.plotBox();
    if (this.simplex) {
      const V = this.tri || (this.tri = this.triangle());
      const a = Math.max(0, x[this.axes[0]]), c = Math.max(0, x[this.axes[1]]), d = Math.max(0, x[this.axes[2]]), s = a + c + d || 1;
      out[0] = (a * V[0][0] + c * V[1][0] + d * V[2][0]) / s;
      out[1] = (a * V[0][1] + c * V[1][1] + d * V[2][1]) / s;
      return out;
    }
    if (!this.is3D()) {
      const u = (x[this.axes[0]] - this.ranges[0][0]) / (this.ranges[0][1] - this.ranges[0][0]);
      const v = (x[this.axes[1]] - this.ranges[1][0]) / (this.ranges[1][1] - this.ranges[1][0]);
      const px = this.pad * b.w, py = this.pad * b.h;
      const cx = b.x + b.w / 2, cy = b.y + b.h / 2;
      out[0] = cx + (u - 0.5) * (b.w - 2 * px) * this.zoom;
      out[1] = cy - (v - 0.5) * (b.h - 2 * py) * this.zoom;
      return out;
    }
    // 3D: q = (horizontal a, horizontal b, vertical c) in [-1, 1]^3.
    const up = this.upAxis, hA = up === 0 ? 1 : 0, hB = up === 2 ? 1 : 2;
    const qa = this.norm(hA, x[this.axes[hA]]), qb = this.norm(hB, x[this.axes[hB]]), qc = this.norm(up, x[this.axes[up]]);
    const ca = Math.cos(this.azim), sa = Math.sin(this.azim), ce = Math.cos(this.elev), se = Math.sin(this.elev);
    const X = ca * qa - sa * qb, Y0 = sa * qa + ca * qb;
    const Y = ce * qc - se * Y0;
    const s = Math.min(b.w, b.h) * 0.5 * (1 - this.pad) / 1.5 * this.zoom * 1.25;
    out[0] = b.x + b.w / 2 + s * X;
    out[1] = b.y + b.h / 2 - s * Y;
    out.depth = ce * Y0 + se * qc;
    return out;
  };
  // Inverse of the 2D projection (clicks); null in 3D.
  Camera.prototype.unproject = function (px, py) {
    if (this.is3D()) return null;
    if (this.simplex) {
      const V = this.tri || (this.tri = this.triangle());
      const det = (V[1][1] - V[2][1]) * (V[0][0] - V[2][0]) + (V[2][0] - V[1][0]) * (V[0][1] - V[2][1]);
      const w0 = ((V[1][1] - V[2][1]) * (px - V[2][0]) + (V[2][0] - V[1][0]) * (py - V[2][1])) / det;
      const w1 = ((V[2][1] - V[0][1]) * (px - V[2][0]) + (V[0][0] - V[2][0]) * (py - V[2][1])) / det;
      return [w0, w1, 1 - w0 - w1];
    }
    const b = this.plotBox();
    const padx = this.pad * b.w, pady = this.pad * b.h, cx = b.x + b.w / 2, cy = b.y + b.h / 2;
    const u = (px - cx) / ((b.w - 2 * padx) * this.zoom) + 0.5, v = -(py - cy) / ((b.h - 2 * pady) * this.zoom) + 0.5;
    return [this.ranges[0][0] + u * (this.ranges[0][1] - this.ranges[0][0]), this.ranges[1][0] + v * (this.ranges[1][1] - this.ranges[1][0])];
  };

  // Tick values: 3 to 7 round numbers inside [lo, hi].
  function ticks(lo, hi, target) {
    target = target || 5;
    const span = hi - lo;
    if (!(span > 0)) return [lo];
    const raw = span / target, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map(function (m) { return m * mag; }).find(function (s) { return span / s <= target + 1; }) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);
    return out;
  }
  function fmt(v) {
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(1).replace("e+", "e");
    return String(+v.toPrecision(4));
  }

  /* Axes for a 2D camera: frame, ticks and labels. Labels are variable or
     parameter names; the fonts are those of the page. */
  function drawAxes(ctx, cam, theme, labels, opts) {
    opts = opts || {};
    const th = DF.THEMES[theme] || DF.THEMES["relab-night"];
    const b = cam.plotBox();
    const x0 = cam.project([cam.ranges[0][0], cam.ranges[1][0]].reduce(function (acc, v, j) { acc[cam.axes[j]] = v; return acc; }, []));
    const x1 = cam.project([cam.ranges[0][1], cam.ranges[1][1]].reduce(function (acc, v, j) { acc[cam.axes[j]] = v; return acc; }, []));
    const L = x0[0], R = x1[0], B = x0[1], T = x1[1];
    ctx.save();
    ctx.strokeStyle = th.grid; ctx.lineWidth = 1;
    ctx.font = (opts.fontSize || 11) + "px Jost, system-ui, sans-serif";
    ctx.fillStyle = th.muted;
    const probe = [];
    ticks(cam.ranges[0][0], cam.ranges[0][1]).forEach(function (v) {
      probe[cam.axes[0]] = v; probe[cam.axes[1]] = cam.ranges[1][0];
      const p = cam.project(probe);
      if (opts.grid) { ctx.beginPath(); ctx.moveTo(p[0], B); ctx.lineTo(p[0], T); ctx.stroke(); }
      ctx.beginPath(); ctx.moveTo(p[0], B); ctx.lineTo(p[0], B + 4); ctx.stroke();
      ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(fmt(v), p[0], B + 6);
    });
    ticks(cam.ranges[1][0], cam.ranges[1][1]).forEach(function (v) {
      probe[cam.axes[0]] = cam.ranges[0][0]; probe[cam.axes[1]] = v;
      const p = cam.project(probe);
      if (opts.grid) { ctx.beginPath(); ctx.moveTo(L, p[1]); ctx.lineTo(R, p[1]); ctx.stroke(); }
      ctx.beginPath(); ctx.moveTo(L - 4, p[1]); ctx.lineTo(L, p[1]); ctx.stroke();
      ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText(fmt(v), L - 7, p[1]);
    });
    ctx.strokeStyle = th.muted; ctx.globalAlpha = 0.6;
    ctx.strokeRect(L, T, R - L, B - T);
    ctx.globalAlpha = 1; ctx.fillStyle = th.ink;
    ctx.font = "italic " + ((opts.fontSize || 11) + 2) + "px 'TeX Gyre Pagella', Palatino, serif";
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(labels[0], (L + R) / 2, B + 22);
    ctx.save(); ctx.translate(L - 40, (T + B) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = "bottom"; ctx.fillText(labels[1], 0, 0); ctx.restore();
    ctx.restore();
    return { L: L, R: R, T: T, B: B };
  }

  // Wireframe of the 3D box, faint, for orientation.
  function drawBox3D(ctx, cam, theme) {
    const th = DF.THEMES[theme] || DF.THEMES["relab-night"];
    const c = [0, 1], p = [];
    const corner = function (i, j, k) { const x = []; x[cam.axes[0]] = cam.ranges[0][i]; x[cam.axes[1]] = cam.ranges[1][j]; x[cam.axes[2]] = cam.ranges[2][k]; return cam.project(x, [0, 0]); };
    ctx.save(); ctx.strokeStyle = th.grid; ctx.lineWidth = 1;
    c.forEach(function (i) { c.forEach(function (j) { c.forEach(function (k) { p.push([i, j, k]); }); }); });
    p.forEach(function (a) {
      p.forEach(function (b) {
        const d = Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
        if (d !== 1 || a.join() > b.join()) return;
        const u = corner(a[0], a[1], a[2]), v = corner(b[0], b[1], b[2]);
        ctx.beginPath(); ctx.moveTo(u[0], u[1]); ctx.lineTo(v[0], v[1]); ctx.stroke();
      });
    });
    ctx.restore();
  }

  DF.autoRanges = autoRanges;
  DF.Camera = Camera;
  DF.ticks = ticks;
  DF.fmt = fmt;
  DF.drawAxes = drawAxes;
  DF.drawBox3D = drawBox3D;
  DF.drawSimplex = function (ctx, cam, theme, labels) {
    const th = DF.THEMES[theme] || DF.THEMES["relab-night"], V = cam.triangle();
    ctx.save(); ctx.strokeStyle = th.grid; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(V[0][0], V[0][1]); ctx.lineTo(V[1][0], V[1][1]); ctx.lineTo(V[2][0], V[2][1]); ctx.closePath(); ctx.stroke();
    ctx.fillStyle = th.muted; ctx.font = "italic 13px 'TeX Gyre Pagella', Palatino, serif"; ctx.textAlign = "center";
    ctx.textBaseline = "bottom"; ctx.fillText(labels[0], V[0][0], V[0][1] - 6);
    ctx.textBaseline = "top"; ctx.fillText(labels[1], V[1][0] - 8, V[1][1] + 6); ctx.fillText(labels[2], V[2][0] + 8, V[2][1] + 6);
    ctx.restore();
  };
})(globalThis.DynFlow = globalThis.DynFlow || {});
