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

  /* Playback rate. A scene plays at `rate` units of model time per second
     of wall-clock time (iterations per second for a map), times the speed
     multiplier of the player. When a scene does not fix the rate, it is
     chosen so that the motion reads alike across models: a trajectory
     travels about one plot width per second (1.5 widths of the axis box in
     a rotating 3D view, which draws the box smaller), a particle about 0.6, a fresh
     orbit of the phase plane about 0.5 over its first second; a time
     series scrolls one window in 10 s; a stroboscopic view adds 10 sections
     per second; a sweep crosses its interval in 25 s.

     Speeds are measured in plot widths, each displayed axis scaled by its
     range, and averaged over half-second windows of playback, as the eye
     averages them. A first estimate comes from the median speed of the
     field (over a grid of the axis box) or of a probe ensemble after a
     transient; it is then corrected by replaying, without drawing, what the
     view shows at that rate: the attractor after the view's warm-up for a
     trajectory, particles born and dying as the flow view makes them, fresh
     orbits of the phase plane. Particle and orbit replays are repeated up to
     three times, rate <- rate x target / measured, within fixed bounds of
     the first estimate. */
  const RATE_TARGET = { trajectory: 1.0, flow: 0.6, phase: 0.5, density: 0.5, strobe: 1.0 };
  const PARAM_FORCING = { periodic: 1, quasiperiodic: 1, ramp: 1, step: 1 };
  function median(a) {
    const b = a.filter(isFinite).sort(function (u, v) { return u - v; });
    return b.length ? b[b.length >> 1] : NaN;
  }
  function planeSpeed(d, axes, widths) {
    let s = 0;
    for (let k = 0; k < axes.length; k++) { const u = d[axes[k]] / widths[axes[k]]; s += u * u; }
    return Math.sqrt(s);
  }
  function fieldSpeed(sys, params, init, axes, ranges) {
    const dim = sys.vars.length, x = Float64Array.from(init), d = new Float64Array(dim), p = Float64Array.from(params);
    const widths = ranges.map(function (r) { return r[1] - r[0]; }), G = axes.length >= 3 ? 9 : 14, out = [];
    const H = function (i) { return x[i]; };
    const walk = function (k) {
      if (k === axes.length) {
        sys.f(0, x, p, d, H);
        out.push(planeSpeed(d, axes, widths));
        return;
      }
      const a = axes[k], r = ranges[a];
      for (let g = 0; g < G; g++) { x[a] = r[0] + (g + 0.5) / G * (r[1] - r[0]); walk(k + 1); }
    };
    walk(0);
    return median(out);
  }
  function attractorSpeed(sys, params, init, axes, ranges, opts) {
    const dim = sys.vars.length, widths = ranges.map(function (r) { return r[1] - r[0]; });
    const meanW = widths.reduce(function (s, w) { return s + w; }, 0) / dim;
    const sim = new DF.Simulator(sys, { n: 6, dt: opts.dt, params: params, init: init, spread: 0.01 * meanW, seed: 12345, perturbations: opts.perturbations || [], keepPositive: !!opts.keepPositive });
    const steps = opts.steps || 4000, lagged = sys.kind === "dde", every = 5, out = [];
    const d = new Float64Array(dim), prev = new Float64Array(sim.X.length), x = new Float64Array(dim);
    const H = function (i) { return x[i]; };
    for (let s = 0; s < steps; s++) {
      if (lagged && s % every === every - 1) prev.set(sim.X);
      sim.step();
      if (s < steps / 4 || s % every) continue;
      for (let k = 0; k < sim.n; k++) {
        if (!sim.alive[k]) continue;
        for (let i = 0; i < dim; i++) x[i] = sim.X[k * dim + i];
        if (lagged) for (let i = 0; i < dim; i++) d[i] = (x[i] - prev[k * dim + i]) / sim.h;
        else sys.f(sim.t, x, sim.p, d, H);
        out.push(planeSpeed(d, axes, widths));
      }
    }
    return median(out);
  }
  /* Replay of a view at a given rate, without drawing: the median over
     members of the path length per second in half-second windows (for fresh
     orbits of the phase plane, over their first second, while they move). */
  function replaySpeed(sys, o, rate, kind) {
    const dim = sys.vars.length, h = o.dt, W = o.ranges.map(function (r) { return r[1] - r[0]; }), axes = o.axes, v = o.view || {};
    const perFrame = rate / (60 * h), MAXSTEPS = 8000;
    let acc = 0;
    const stepsNow = function () { acc += perFrame; const n = Math.floor(acc); acc -= n; return n; };
    const disp = function (X, k, prev) { let s = 0; for (let j = 0; j < axes.length; j++) { const i = axes[j], u = (X[k * dim + i] - prev[k * dim + i]) / W[i]; s += u * u; } return Math.sqrt(s); };
    const out = [], WIN = 30;
    const meanW = W.reduce(function (s, w) { return s + w; }, 0) / dim;
    const common = { dt: h, params: o.params, init: o.init, seed: 4321, keepPositive: !!o.keepPositive, box: o.ranges };
    if (kind === "phase") {
      const r = new DF.RNG(99), forcing = (o.perturbations || []).filter(function (q) { return q.kind in PARAM_FORCING; });
      for (let t = 0; t < 8; t++) {
        const x0 = Float64Array.from(o.init);
        if (t) axes.forEach(function (i) { x0[i] = r.range(o.ranges[i][0], o.ranges[i][1]); });
        const sim = new DF.Simulator(sys, Object.assign({}, common, { n: 1, init: x0, deterministic: true, perturbations: forcing }));
        const prev = Float64Array.from(sim.X);
        let path = 0, steps = 0;
        acc = 0;
        const gone = function () { for (let i = 0; i < dim; i++) if (sim.X[i] < o.ranges[i][0] - 2 * W[i] || sim.X[i] > o.ranges[i][1] + 2 * W[i]) return true; return false; };
        for (let f = 0; f < 60 && steps < MAXSTEPS; f++) {
          const n = stepsNow(); steps += n;
          for (let k = 0; k < n; k++) sim.step();
          // The phase view drops an orbit that leaves the box by two widths.
          if (!sim.alive[0] || gone()) break;
          path += disp(sim.X, 0, prev); prev.set(sim.X);
        }
        out.push(path);
      }
      return median(out);
    }
    let sim, n, life = null, age = null, spawn = null;
    if (kind === "flow" && v.life !== "inf" && v.life) {
      n = 24;
      sim = new DF.Simulator(sys, Object.assign({}, common, { n: n, perturbations: o.perturbations || [] }));
      const r = sim.rng, x = new Float64Array(dim);
      spawn = function (k) {
        const mode = v.spawn === "mixed" ? (r.uniform() < (v.spawnMix === undefined ? 0.25 : v.spawnMix) ? "init" : "box") : v.spawn || "box";
        for (let i = 0; i < dim; i++) x[i] = mode === "init" ? o.init[i] + (o.spread || 0.02) * W[i] * r.normal() : r.range(o.ranges[i][0], o.ranges[i][1]);
        sim.setMember(k, x);
      };
      life = new Float64Array(n); age = new Float64Array(n);
      for (let k = 0; k < n; k++) { spawn(k); life[k] = v.life[0] + r.uniform() * (v.life[1] - v.life[0]); age[k] = r.uniform() * life[k]; }
    } else {
      n = 6;
      sim = new DF.Simulator(sys, Object.assign({}, common, { n: n, spread: kind === "flow" ? o.spread : 0.01 * meanW, initMode: kind === "flow" ? o.initMode : "point", perturbations: o.perturbations || [] }));
    }
    const warmFrames = v.warmup === undefined ? (kind === "flow" ? 60 : 0) : v.warmup;
    const warm = Math.min(20000, Math.max(kind === "attractor" ? 1000 : 0, Math.round(warmFrames * perFrame)));
    for (let k = 0; k < warm; k++) sim.step();
    const frames = Math.max(60, Math.min(kind === "flow" ? 150 : 600, Math.floor(MAXSTEPS / Math.max(perFrame, 1e-9))));
    const prev = Float64Array.from(sim.X), wl = new Float64Array(n), wc = new Int32Array(n);
    acc = 0;
    for (let f = 0; f < frames; f++) {
      const m = stepsNow();
      for (let k = 0; k < m; k++) sim.step();
      for (let k = 0; k < n; k++) {
        if (life) {
          age[k] += 1;
          if (!sim.alive[k] || age[k] > life[k]) { spawn(k); age[k] = 0; wl[k] = 0; wc[k] = 0; for (let i = 0; i < dim; i++) prev[k * dim + i] = sim.X[k * dim + i]; continue; }
        } else if (!sim.alive[k]) continue;
        const d = disp(sim.X, k, prev);
        if (!isFinite(d)) continue;
        wl[k] += d;
        if (++wc[k] === WIN) { out.push(wl[k] * 60 / WIN); wl[k] = 0; wc[k] = 0; }
      }
      prev.set(sim.X);
    }
    return median(out);
  }
  function calibrateRate(sys, o) {
    const v = o.view || {}, type = v.type || "flow", discrete = sys.time === "discrete", dt = o.dt;
    const round = function (r) { return +r.toPrecision(2); };
    const windowOf = function (def) { return v.window || def; };
    if (type === "timeseries") return { rate: round(windowOf(discrete ? 100 : 50) / 10), basis: "window" };
    if (type === "density" && (v.mode || (sys.vars.length === 1 ? "carpet" : "map")) === "carpet") return { rate: round(windowOf(60) / 10), basis: "window" };
    if (type === "cobweb") return { rate: 6, basis: "iterations" };
    if (type === "sweep") return { rate: round(1 / (Math.max(1e-9, o.sweepSpeed) * 25)), basis: "sweep" };
    if (type === "strobe" && v.mode !== "section") return { rate: round(10 * o.period), basis: "sections" };
    if (discrete || type === "orbit") return { rate: discrete ? 60 : round(240 * dt), basis: "iterations" };
    // The rotating 3D camera draws the axis box at about 0.77 of the shorter
    // side of the plot, and projection shortens a displacement by about 0.82
    // on average, so a 3D scene aims 1.5 times higher to look as fast on screen.
    const threeD = o.axes.length >= 3 && v.projection !== "simplex" && (type === "trajectory" || type === "flow");
    const target = (RATE_TARGET[type] || 0.6) * (threeD ? 1.5 : 1), axes = o.axes;
    const lo = 60 * dt / 8, hi = 60 * dt * 4000, clamp = function (r) { return Math.min(hi, Math.max(lo, r)); };
    const onAttractor = type === "trajectory" || type === "density" || type === "strobe" || (type === "flow" && (v.life === "inf" || v.spawn === "init"));
    let speed = NaN, basis = "";
    if (onAttractor) { speed = attractorSpeed(sys, o.params, o.init, axes, o.ranges, o); basis = "attractor"; }
    if (!(speed > 1e-12)) { speed = fieldSpeed(sys, o.params, o.init, axes, o.ranges); basis = "field"; }
    if (!(speed > 1e-12)) return { rate: round(120 * dt), basis: "default" };
    const r0 = clamp(target / speed);
    let r = r0;
    try {
      if (type === "trajectory" || type === "density" || type === "strobe") {
        const m = replaySpeed(sys, o, r, "attractor");
        if (m > 1e-9) { r = clamp(r * Math.min(5, Math.max(0.2, target / m))); basis = "replay"; }
      } else if (type === "flow" || type === "phase") {
        const bound = type === "flow" ? 20 : 4;
        for (let it = 0; it < 3; it++) {
          const m = replaySpeed(sys, o, r, type);
          if (!(m > 1e-9)) break;
          const f = target / m;
          r = clamp(Math.min(r0 * bound, Math.max(r0 / bound, r * Math.min(10, Math.max(0.1, f)))));
          basis = "replay";
          if (Math.abs(f - 1) < 0.1) break;
        }
      }
    } catch (e) { r = r0; }
    return { rate: round(r), basis: basis };
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

  // The axes of drawAxes as SVG elements, for vector exports.
  function esc(t) { return String(t).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function axesSVG(cam, theme, labels, opts) {
    opts = opts || {};
    const th = DF.THEMES[theme] || DF.THEMES["relab-night"], fs = opts.fontSize || 11;
    const corner = function (i, j) { const x = []; x[cam.axes[0]] = cam.ranges[0][i]; x[cam.axes[1]] = cam.ranges[1][j]; return cam.project(x, [0, 0]); };
    const c0 = corner(0, 0), c1 = corner(1, 1), L = c0[0], R = c1[0], B = c0[1], T = c1[1];
    const f = function (v) { return v.toFixed(1); };
    let s = '<g font-family="Jost, sans-serif" font-size="' + fs + '" fill="' + th.muted + '" stroke="none">';
    let lines = "";
    const probe = [];
    ticks(cam.ranges[0][0], cam.ranges[0][1]).forEach(function (v) {
      probe[cam.axes[0]] = v; probe[cam.axes[1]] = cam.ranges[1][0];
      const p = cam.project(probe, [0, 0]);
      if (opts.grid) lines += "M" + f(p[0]) + " " + f(B) + "V" + f(T);
      lines += "M" + f(p[0]) + " " + f(B) + "v4";
      s += '<text x="' + f(p[0]) + '" y="' + f(B + 6 + fs * 0.8) + '" text-anchor="middle">' + esc(fmt(v)) + "</text>";
    });
    ticks(cam.ranges[1][0], cam.ranges[1][1]).forEach(function (v) {
      probe[cam.axes[0]] = cam.ranges[0][0]; probe[cam.axes[1]] = v;
      const p = cam.project(probe, [0, 0]);
      if (opts.grid) lines += "M" + f(L) + " " + f(p[1]) + "H" + f(R);
      lines += "M" + f(L - 4) + " " + f(p[1]) + "H" + f(L);
      s += '<text x="' + f(L - 7) + '" y="' + f(p[1] + fs * 0.35) + '" text-anchor="end">' + esc(fmt(v)) + "</text>";
    });
    s += "</g>";
    s = '<path d="' + lines + '" fill="none" stroke="' + th.grid + '" stroke-width="1"/>' + s;
    s += '<rect x="' + f(L) + '" y="' + f(T) + '" width="' + f(R - L) + '" height="' + f(B - T) + '" fill="none" stroke="' + th.muted + '" stroke-opacity="0.6"/>';
    s += '<g font-family="\'TeX Gyre Pagella\', Palatino, serif" font-style="italic" font-size="' + (fs + 2) + '" fill="' + th.ink + '">' +
      '<text x="' + f((L + R) / 2) + '" y="' + f(B + 22 + fs) + '" text-anchor="middle">' + esc(labels[0]) + "</text>" +
      '<text transform="translate(' + f(L - 40) + " " + f((T + B) / 2) + ') rotate(-90)" text-anchor="middle">' + esc(labels[1]) + "</text></g>";
    return s;
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
  DF.calibrateRate = calibrateRate;
  DF.fieldSpeed = fieldSpeed;
  DF.attractorSpeed = attractorSpeed;
  DF.replaySpeed = replaySpeed;
  DF.axesSVG = axesSVG;
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
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
