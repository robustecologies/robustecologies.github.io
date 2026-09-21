// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Views. Each view draws one kind of figure from a running Simulator on the
   three layers of a Player: base (static, redrawn on demand), trail
   (accumulating, faded each frame) and top (cleared each frame).

     flow        particle ensemble with fading trails, 2D or rotating 3D
     trajectory  a few long orbits with a gradient tail, 2D or rotating 3D
     timeseries  scrolling time series of chosen variables
     phase       vector field, nullclines, classified equilibria, orbits
     sweep       slow parameter sweep over the equilibrium branches (hysteresis)
     orbit       bifurcation diagram, built column by column
     density     ensemble density: 2D heat map, or a time carpet in 1D
     strobe      stroboscopic samples every period T, or a Poincare section
     cobweb      cobweb diagram of a one-dimensional map

   Every view implements init(P), frame(P), and optionally drawStatic(P),
   pointer(P, kind, x, y, ev), svg(P) and legend(P). */
(function (DF) {
  "use strict";

  const V = {};

  // ------------------------------------------------------------ helpers
  function colorFor(P, k, x, extra) {
    const s = P.scene.style, pal = P.palette;
    switch (s.colorBy) {
      case "member": return pal[k % pal.length];
      case "dominant": {
        let best = 0;
        for (let i = 1; i < x.length; i++) if (x[i] > x[best]) best = i;
        return pal[best % pal.length];
      }
      case "speed": case "var": case "age": case "time": {
        const u = Math.round(Math.min(1, Math.max(0, extra)) * 23) / 23;
        const c = DF.rampRGB(s.ramp, u);
        return "rgb(" + c[0] + "," + c[1] + "," + c[2] + ")";
      }
      default: return pal[0];
    }
  }
  function Buckets() { this.map = new Map(); }
  Buckets.prototype.seg = function (c, x0, y0, x1, y1) {
    let a = this.map.get(c); if (!a) { a = []; this.map.set(c, a); }
    a.push(x0, y0, x1, y1);
  };
  Buckets.prototype.strokeAll = function (ctx, width, alpha) {
    ctx.lineWidth = width; ctx.globalAlpha = alpha; ctx.lineCap = "round";
    this.map.forEach(function (a, c) {
      ctx.strokeStyle = c; ctx.beginPath();
      for (let i = 0; i < a.length; i += 4) { ctx.moveTo(a[i], a[i + 1]); ctx.lineTo(a[i + 2], a[i + 3]); }
      ctx.stroke();
    });
    ctx.globalAlpha = 1; this.map.clear();
  };
  Buckets.prototype.dotAll = function (ctx, size, alpha) {
    ctx.globalAlpha = alpha;
    this.map.forEach(function (a, c) {
      ctx.fillStyle = c;
      for (let i = 0; i < a.length; i += 4) ctx.fillRect(a[i] - size / 2, a[i + 1] - size / 2, size, size);
    });
    ctx.globalAlpha = 1; this.map.clear();
  };
  function speedOf(P, x) {
    const d = P.tmpDx;
    P.sys.f(P.sim.t, x, P.sim.p, d, function (i) { return x[i]; });
    let s = 0;
    for (let i = 0; i < d.length; i++) { const r = P.fullRanges[i]; const w = r[1] - r[0]; s += (d[i] / w) * (d[i] / w); }
    return Math.sqrt(s);
  }
  function extraFor(P, k, x, age) {
    const s = P.scene.style;
    if (s.colorBy === "speed") { const v = speedOf(P, x); P.speedMax = Math.max(P.speedMax * 0.9995, v); return v / (P.speedMax || 1); }
    if (s.colorBy === "var") { const i = Math.max(0, P.sys.vars.indexOf(s.colorVar)); const r = P.fullRanges[i]; return (x[i] - r[0]) / (r[1] - r[0]); }
    if (s.colorBy === "age") return age;
    return 0;
  }
  function randomInBox(P, out) {
    const r = P.sim.rng;
    for (let i = 0; i < P.sys.vars.length; i++) out[i] = r.range(P.fullRanges[i][0], P.fullRanges[i][1]);
    return out;
  }
  function spawnState(P, out) {
    const v = P.scene.view, r = P.sim.rng, init = P.sim.init;
    const mode = v.spawn === "mixed" ? (r.uniform() < (v.spawnMix === undefined ? 0.25 : v.spawnMix) ? "init" : "box") : v.spawn;
    if (mode === "init") {
      for (let i = 0; i < out.length; i++) { const w = P.fullRanges[i][1] - P.fullRanges[i][0]; out[i] = init[i] + (P.scene.spread || 0.02) * w * r.normal(); }
      return out;
    }
    return randomInBox(P, out);
  }
  function outOfView(P, x) {
    for (let i = 0; i < x.length; i++) {
      const r = P.fullRanges[i], w = r[1] - r[0];
      if (x[i] < r[0] - 2 * w || x[i] > r[1] + 2 * w) return true;
    }
    return false;
  }
  function svgHead(P) {
    return '<svg xmlns="http://www.w3.org/2000/svg" width="' + P.w + '" height="' + P.h + '" viewBox="0 0 ' + P.w + " " + P.h + '">';
  }
  function polyline(pts, color, width, alpha) {
    if (pts.length < 4) return "";
    let d = "M" + pts[0].toFixed(2) + " " + pts[1].toFixed(2);
    for (let i = 2; i < pts.length; i += 2) d += "L" + pts[i].toFixed(2) + " " + pts[i + 1].toFixed(2);
    return '<path d="' + d + '" fill="none" stroke="' + color + '" stroke-width="' + width + '" stroke-opacity="' + alpha + '" stroke-linejoin="round" stroke-linecap="round"/>';
  }
  function rotate3D(P) {
    if (P.cam.is3D() && !P.dragging) P.cam.azim += (P.scene.view.rotate || 0) * 0.01;
  }

  // ------------------------------------------------------------ flow
  V.flow = {
    label: "Flow (particles)",
    init: function (P) {
      const n = P.sim.n, st = new Float64Array(P.sys.vars.length);
      this.age = new Float32Array(n); this.life = new Float32Array(n);
      this.prev = new Float32Array(2 * n); this.has = new Uint8Array(n);
      this.buck = new Buckets(); this.pt = [0, 0];
      for (let k = 0; k < n; k++) {
        spawnState(P, st); P.sim.setMember(k, st);
        this.life[k] = this.newLife(P); this.age[k] = P.sim.rng.uniform() * this.life[k];
      }
      // Warm-up so that the first frame already shows the flow.
      const warm = P.scene.view.warmup === undefined ? 60 : P.scene.view.warmup;
      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();
    },
    newLife: function (P) {
      const L = P.scene.view.life;
      if (!L || L[1] === Infinity || L === "inf") return Infinity;
      return L[0] + P.sim.rng.uniform() * (L[1] - L[0]);
    },
    frame: function (P) {
      const sim = P.sim, n = sim.n, dim = P.sys.vars.length, cam = P.cam, st = P.scene.style;
      rotate3D(P);
      P.fade(st.fade);
      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();
      const x = P.tmpX, discrete = P.sys.time === "discrete", st2 = new Float64Array(dim);
      for (let k = 0; k < n; k++) {
        this.age[k]++;
        if (!sim.alive[k] || this.age[k] > this.life[k]) {
          spawnState(P, st2); sim.setMember(k, st2); this.age[k] = 0; this.life[k] = this.newLife(P); this.has[k] = 0;
          continue;
        }
        sim.member(k, x);
        if (outOfView(P, x)) { this.age[k] = this.life[k] + 1; continue; }
        const q = cam.project(x, this.pt);
        const col = colorFor(P, k, x, extraFor(P, k, x, this.age[k] / (isFinite(this.life[k]) ? this.life[k] : 1)));
        if (discrete) this.buck.seg(col, q[0], q[1], q[0], q[1]);
        else if (this.has[k]) this.buck.seg(col, this.prev[2 * k], this.prev[2 * k + 1], q[0], q[1]);
        this.prev[2 * k] = q[0]; this.prev[2 * k + 1] = q[1]; this.has[k] = 1;
      }
      const ctx = P.ctx.trail;
      ctx.globalCompositeOperation = P.theme.blend;
      if (discrete) this.buck.dotAll(ctx, st.pointSize, st.alpha); else this.buck.strokeAll(ctx, st.lineWidth, st.alpha);
      ctx.globalCompositeOperation = "source-over";
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down" || P.cam.is3D()) return false;
      const u = P.cam.unproject(px, py);
      const s = Float64Array.from(P.sim.init), r = P.sim.rng;
      if (P.cam.simplex) {
        if (u.some(function (c) { return c <= 0.002; })) return false;
        const tot = P.cam.axes.reduce(function (acc, i) { return acc + P.sim.init[i]; }, 0) || 1;
        for (let j = 0; j < 80; j++) {
          const k = Math.floor(r.uniform() * P.sim.n);
          P.cam.axes.forEach(function (i, m) { s[i] = tot * Math.max(1e-6, u[m] + 0.005 * r.normal()); });
          P.sim.setMember(k, s); this.age[k] = 0; this.has[k] = 0;
        }
        return true;
      }
      for (let j = 0; j < 80; j++) {
        const k = Math.floor(r.uniform() * P.sim.n);
        for (let i = 0; i < s.length; i++) s[i] = P.sim.init[i];
        const w0 = P.cam.ranges[0][1] - P.cam.ranges[0][0], w1 = P.cam.ranges[1][1] - P.cam.ranges[1][0];
        s[P.cam.axes[0]] = u[0] + 0.01 * w0 * r.normal(); s[P.cam.axes[1]] = u[1] + 0.01 * w1 * r.normal();
        P.sim.setMember(k, s); this.age[k] = 0; this.has[k] = 0;
      }
      return true;
    }
  };

  // ------------------------------------------------------------ trajectory
  V.trajectory = {
    label: "Trajectory",
    init: function (P) {
      const n = P.sim.n, dim = P.sys.vars.length, L = P.scene.view.tail || 2500;
      this.L = L; this.buf = new Float64Array(n * L * dim); this.len = new Int32Array(n); this.head = new Int32Array(n);
      const warm = P.scene.view.warmup === undefined ? 0 : P.scene.view.warmup;
      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();
    },
    push: function (P) {
      const n = P.sim.n, dim = P.sys.vars.length, L = this.L;
      for (let k = 0; k < n; k++) {
        if (!P.sim.alive[k]) continue;
        const h = this.head[k], o = (k * L + h) * dim;
        for (let i = 0; i < dim; i++) this.buf[o + i] = P.sim.X[k * dim + i];
        this.head[k] = (h + 1) % L; this.len[k] = Math.min(this.len[k] + 1, L);
      }
    },
    points: function (P, k) {
      const dim = P.sys.vars.length, L = this.L, len = this.len[k], out = new Float32Array(2 * len), x = P.tmpX, q = [0, 0];
      const start = (this.head[k] - len + L) % L;
      for (let j = 0; j < len; j++) {
        const o = (k * L + (start + j) % L) * dim;
        for (let i = 0; i < dim; i++) x[i] = this.buf[o + i];
        P.cam.project(x, q); out[2 * j] = q[0]; out[2 * j + 1] = q[1];
      }
      return out;
    },
    frame: function (P) {
      rotate3D(P);
      const every = P.scene.view.sampleEvery || 1;
      for (let s = 0; s < P.scene.stepsPerFrame; s++) { P.sim.step(); if (s % every === 0) this.push(P); }
      const ctx = P.ctx.top, st = P.scene.style, n = P.sim.n, G = 28;
      ctx.clearRect(0, 0, P.w, P.h);
      if (P.cam.is3D() && P.scene.view.box) DF.drawBox3D(ctx, P.cam, P.scene.style.theme);
      ctx.globalCompositeOperation = P.theme.blend; ctx.lineCap = "round"; ctx.lineJoin = "round";
      const discrete = P.sys.time === "discrete";
      for (let k = 0; k < n; k++) {
        const pts = this.points(P, k), m = pts.length / 2;
        if (m < 2) continue;
        for (let g = 0; g < G; g++) {
          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
          if (b <= a) continue;
          const u = (g + 1) / G;
          let col;
          if (st.colorBy === "time" || st.colorBy === "age") { const c = DF.rampRGB(st.ramp, u); col = "rgb(" + c.join(",") + ")"; }
          else col = P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length];
          ctx.strokeStyle = col; ctx.fillStyle = col;
          ctx.globalAlpha = st.alpha * (0.08 + 0.92 * Math.pow(u, 1.4));
          if (discrete) { for (let j = a; j <= b; j++) ctx.fillRect(pts[2 * j] - st.pointSize / 2, pts[2 * j + 1] - st.pointSize / 2, st.pointSize, st.pointSize); continue; }
          ctx.lineWidth = st.lineWidth;
          ctx.beginPath(); ctx.moveTo(pts[2 * a], pts[2 * a + 1]);
          for (let j = a + 1; j <= b; j++) ctx.lineTo(pts[2 * j], pts[2 * j + 1]);
          ctx.stroke();
        }
        ctx.globalAlpha = 1;
        ctx.fillStyle = P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length];
        ctx.beginPath(); ctx.arc(pts[pts.length - 2], pts[pts.length - 1], Math.max(2.5, st.lineWidth * 2), 0, 7); ctx.fill();
      }
      ctx.globalCompositeOperation = "source-over"; ctx.globalAlpha = 1;
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down" || P.cam.is3D()) return false;
      const u = P.cam.unproject(px, py), s = Float64Array.from(P.sim.init);
      s[P.cam.axes[0]] = u[0]; s[P.cam.axes[1]] = u[1];
      const k = (this.nextSeed = ((this.nextSeed || 0) + 1) % P.sim.n);
      P.sim.setMember(k, s); this.len[k] = 0; this.head[k] = 0;
      return true;
    },
    svg: function (P) {
      let out = "";
      const st = P.scene.style;
      for (let k = 0; k < P.sim.n; k++) {
        const pts = Array.from(this.points(P, k));
        out += polyline(pts, P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length], st.lineWidth, st.alpha);
      }
      return out;
    }
  };

  // ------------------------------------------------------------ timeseries
  V.timeseries = {
    label: "Time series",
    axes: true,
    init: function (P) {
      const vars = (P.scene.view.vars && P.scene.view.vars.length ? P.scene.view.vars : P.sys.vars.slice(0, 4)).filter(function (v) { return P.sys.vars.indexOf(v) >= 0; });
      this.vi = vars.map(function (v) { return P.sys.vars.indexOf(v); });
      this.members = Math.min(P.sim.n, P.scene.view.members || 1);
      this.cap = 4000; this.T = new Float64Array(this.cap); this.Y = new Float64Array(this.cap * this.members * this.vi.length);
      this.len = 0; this.head = 0;
      let lo = Infinity, hi = -Infinity;
      this.vi.forEach(function (i) { lo = Math.min(lo, P.fullRanges[i][0]); hi = Math.max(hi, P.fullRanges[i][1]); });
      this.yr = P.scene.view.yRange || [lo, hi];
    },
    frame: function (P) {
      const nv = this.vi.length, M = this.members;
      for (let s = 0; s < P.scene.stepsPerFrame; s++) {
        P.sim.step();
        const h = this.head;
        this.T[h] = P.sim.t;
        for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) this.Y[(h * M + m) * nv + j] = P.sim.X[m * P.sys.vars.length + this.vi[j]];
        this.head = (h + 1) % this.cap; this.len = Math.min(this.len + 1, this.cap);
      }
      const span = P.scene.view.window || (P.sys.time === "discrete" ? 100 : 50);
      const t0 = Math.max(0, P.sim.t - span), t1 = t0 + span;
      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);
      const cam = new DF.Camera([0, 1], [[t0, t1], this.yr], { pad: 0.02 });
      cam.resize(P.w, P.h, P.inset);
      DF.drawAxes(ctx, cam, P.scene.style.theme, [P.sys.time === "discrete" ? "n" : "t", this.vi.map(function (i) { return P.sys.vars[i]; }).join(", ")], { grid: true });
      const st = P.scene.style, q = [0, 0], x = [0, 0];
      ctx.save(); const b = cam.plotBox(); ctx.beginPath(); ctx.rect(b.x, b.y, b.w, b.h); ctx.clip();
      ctx.lineJoin = "round";
      const start = (this.head - this.len + this.cap) % this.cap;
      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) {
        ctx.strokeStyle = P.palette[j % P.palette.length]; ctx.globalAlpha = M > 1 ? Math.max(0.25, st.alpha * 0.6) : st.alpha; ctx.lineWidth = st.lineWidth;
        ctx.beginPath(); let started = false;
        for (let r = 0; r < this.len; r++) {
          const h = (start + r) % this.cap;
          if (this.T[h] < t0) continue;
          x[0] = this.T[h]; x[1] = this.Y[(h * M + m) * nv + j];
          cam.project(x, q);
          if (!started) { ctx.moveTo(q[0], q[1]); started = true; } else ctx.lineTo(q[0], q[1]);
        }
        ctx.stroke();
      }
      ctx.restore(); ctx.globalAlpha = 1;
      this.cam = cam;
    },
    legend: function (P) { return this.vi.map(function (i, j) { return [P.sys.vars[i], P.palette[j % P.palette.length]]; }); },
    svg: function (P) {
      const nv = this.vi.length, M = this.members, cam = this.cam, t0 = cam.ranges[0][0];
      let out = "";
      const start = (this.head - this.len + this.cap) % this.cap, q = [0, 0], x = [0, 0];
      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) {
        const pts = [];
        for (let r = 0; r < this.len; r++) { const h = (start + r) % this.cap; if (this.T[h] < t0) continue; x[0] = this.T[h]; x[1] = this.Y[(h * M + m) * nv + j]; cam.project(x, q); pts.push(q[0], q[1]); }
        out += polyline(pts, P.palette[j % P.palette.length], P.scene.style.lineWidth, P.scene.style.alpha);
      }
      return out;
    }
  };

  // ------------------------------------------------------------ phase
  // Marching squares on a grid of values g (nx by ny), level 0.
  function contour(g, nx, ny, X, Y) {
    const segs = [];
    const lerp = function (a, b) { return a / (a - b); };
    for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
      const a = g[j * nx + i], b = g[j * nx + i + 1], c = g[(j + 1) * nx + i + 1], d = g[(j + 1) * nx + i];
      if (!isFinite(a + b + c + d)) continue;
      const pts = [];
      if ((a > 0) !== (b > 0)) pts.push([X(i + lerp(a, b)), Y(j)]);
      if ((b > 0) !== (c > 0)) pts.push([X(i + 1), Y(j + lerp(b, c))]);
      if ((c > 0) !== (d > 0)) pts.push([X(i + 1 - lerp(c, d)), Y(j + 1)]);
      if ((d > 0) !== (a > 0)) pts.push([X(i), Y(j + 1 - lerp(d, a))]);
      if (pts.length === 2) segs.push(pts[0], pts[1]);
      else if (pts.length === 4) segs.push(pts[0], pts[1], pts[2], pts[3]);
    }
    return segs;
  }

  V.phase = {
    label: "Phase plane",
    axes: true,
    init: function (P) {
      this.trails = []; this.maxTrails = 40;
      this.static = null;
      const n0 = Math.min(P.sim.n, P.scene.view.seeds === undefined ? 6 : P.scene.view.seeds);
      const r = P.sim.rng, s = new Float64Array(P.sys.vars.length);
      for (let k = 0; k < n0; k++) {
        for (let i = 0; i < s.length; i++) s[i] = P.sim.init[i];
        if (k > 0) {
          const a0 = P.cam.axes[0], a1 = P.cam.axes[1];
          s[a0] = r.range(P.cam.ranges[0][0], P.cam.ranges[0][1]); s[a1] = r.range(P.cam.ranges[1][0], P.cam.ranges[1][1]);
        }
        this.addTrail(P, s);
      }
    },
    addTrail: function (P, s) {
      if (this.trails.length >= this.maxTrails) this.trails.shift();
      this.trails.push({ x: Float64Array.from(s), pts: [], life: 0 });
    },
    drawStatic: function (P) {
      const ctx = P.ctx.base, cam = P.cam, sys = P.sys, v = P.scene.view, st = P.scene.style;
      const th = P.theme;
      const frame = DF.drawAxes(ctx, cam, st.theme, [sys.vars[cam.axes[0]], sys.vars[cam.axes[1]]], { grid: false });
      this.frameBox = frame;
      const a0 = cam.axes[0], a1 = cam.axes[1], t = P.sim.t, p = P.sim.p;
      const x = Float64Array.from(P.sim.init), d = new Float64Array(sys.vars.length);
      const H = function (i) { return x[i]; };
      const R0 = cam.ranges[0], R1 = cam.ranges[1];
      ctx.save(); ctx.beginPath(); ctx.rect(frame.L, frame.T, frame.R - frame.L, frame.B - frame.T); ctx.clip();
      // Vector field: arrows on a grid, length by log speed.
      if (v.field !== "none") {
        const nx = v.fieldDensity || 22, ny = Math.round(nx * (frame.B - frame.T) / (frame.R - frame.L));
        const cellW = (frame.R - frame.L) / nx, cellH = (frame.B - frame.T) / ny;
        const q = [0, 0];
        let vmax = 0;
        const vals = [];
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
          x[a0] = R0[0] + (i + 0.5) / nx * (R0[1] - R0[0]); x[a1] = R1[1] - (j + 0.5) / ny * (R1[1] - R1[0]);
          sys.f(t, x, p, d, H);
          const u = d[a0] / (R0[1] - R0[0]) * (frame.R - frame.L), w = -d[a1] / (R1[1] - R1[0]) * (frame.B - frame.T);
          const m = Math.hypot(u, w); vmax = Math.max(vmax, isFinite(m) ? m : 0);
          vals.push([frame.L + (i + 0.5) * cellW, frame.T + (j + 0.5) * cellH, u, w, m]);
        }
        ctx.lineWidth = 1;
        vals.forEach(function (e) {
          if (!(e[4] > 0) || !isFinite(e[4])) return;
          const rel = Math.log1p(9 * e[4] / vmax) / Math.log(10);
          const len = Math.min(cellW, cellH) * 0.42 * (0.35 + 0.65 * rel), ux = e[2] / e[4], uy = e[3] / e[4];
          const c = DF.rampRGB(st.ramp, 0.25 + 0.75 * rel);
          ctx.strokeStyle = th.dark ? "rgba(" + c.join(",") + ",0.55)" : "rgba(" + c.join(",") + ",0.75)";
          const x0 = e[0] - ux * len, y0 = e[1] - uy * len, x1 = e[0] + ux * len, y1 = e[1] + uy * len;
          ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1);
          ctx.moveTo(x1, y1); ctx.lineTo(x1 - 4 * ux + 2.5 * uy, y1 - 4 * uy - 2.5 * ux);
          ctx.moveTo(x1, y1); ctx.lineTo(x1 - 4 * ux - 2.5 * uy, y1 - 4 * uy + 2.5 * ux);
          ctx.stroke();
        });
      }
      // Nullclines f_a0 = 0 and f_a1 = 0 by marching squares.
      if (v.nullclines !== false && P.sys.time === "continuous") {
        const nx = 160, ny = 160, g0 = new Float64Array(nx * ny), g1 = new Float64Array(nx * ny);
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
          x[a0] = R0[0] + i / (nx - 1) * (R0[1] - R0[0]); x[a1] = R1[0] + j / (ny - 1) * (R1[1] - R1[0]);
          sys.f(t, x, p, d, H);
          g0[j * nx + i] = d[a0]; g1[j * nx + i] = d[a1];
        }
        const X = function (i) { return R0[0] + i / (nx - 1) * (R0[1] - R0[0]); }, Y = function (j) { return R1[0] + j / (ny - 1) * (R1[1] - R1[0]); };
        const self = this; self.nullSegs = [];
        [g0, g1].forEach(function (g, idx) {
          const segs = contour(g, nx, ny, X, Y);
          const col = P.palette[(idx + 1) % P.palette.length];
          ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.setLineDash(idx ? [6, 4] : []); ctx.globalAlpha = 0.9;
          ctx.beginPath();
          const q = [0, 0], z = Float64Array.from(P.sim.init);
          const pts = [];
          for (let s = 0; s < segs.length; s += 2) {
            z[a0] = segs[s][0]; z[a1] = segs[s][1]; cam.project(z, q); ctx.moveTo(q[0], q[1]); pts.push(q[0], q[1]);
            z[a0] = segs[s + 1][0]; z[a1] = segs[s + 1][1]; cam.project(z, q); ctx.lineTo(q[0], q[1]); pts.push(q[0], q[1]);
          }
          ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
          self.nullSegs.push({ pts: pts, color: col, dash: idx ? "6 4" : "" });
        });
      }
      ctx.restore();
      // Equilibria in the plane (other variables held at their initial values).
      this.equilibria = [];
      if (v.equilibria !== false && P.sys.time === "continuous") {
        const sub = { vars: [sys.vars[a0], sys.vars[a1]], time: "continuous", f: function (tt, y, pp, out) { x[a0] = y[0]; x[a1] = y[1]; sys.f(tt, x, pp, d, H); out[0] = d[a0]; out[1] = d[a1]; } };
        try { this.equilibria = DF.findEquilibria(sub, p, [R0, R1], { t: t, seeds: 80 }); } catch (e) { this.equilibria = []; }
        const q = [0, 0], z = Float64Array.from(P.sim.init);
        this.equilibria.forEach(function (e) {
          z[a0] = e.x[0]; z[a1] = e.x[1]; cam.project(z, q);
          ctx.lineWidth = 2; ctx.strokeStyle = th.ink; ctx.fillStyle = th.ink;
          ctx.beginPath(); ctx.arc(q[0], q[1], 5.5, 0, 7);
          if (e.stable) ctx.fill();
          else if (/saddle/.test(e.type)) { ctx.fillStyle = P.palette[0]; ctx.fill(); ctx.stroke(); }
          else { ctx.fillStyle = th.dark ? "#0b0620" : "#ffffff"; ctx.fill(); ctx.stroke(); }
        });
      }
    },
    frame: function (P) {
      const sys = P.sys, dim = sys.vars.length, h = P.sim.h, cam = P.cam, st = P.scene.style;
      if (sys.usesTime || P.sim.perturbations.some(function (q) { return q.enabled !== false && DF.PARAM_PERTURBATIONS.indexOf(q.kind) >= 0; })) {
        if ((P.frameCount % 20) === 0) P.redrawStatic();
      }
      for (let s = 0; s < P.scene.stepsPerFrame; s++) P.sim.step();
      const rk = this.rk || (this.rk = DF.makeRK4(dim)), q = [0, 0], tt = P.sim.t;
      const maxPts = P.scene.view.tail || 1500;
      const self = this;
      this.trails.forEach(function (tr) {
        for (let s = 0; s < P.scene.stepsPerFrame; s++) {
          if (sys.time === "discrete") { sys.f(tt, tr.x, P.sim.p, P.tmpDx); tr.x.set(P.tmpDx); }
          else rk(sys.f, tt, tr.x, P.sim.p, h, function (i) { return tr.x[i]; });
          if (!tr.x.every(isFinite) || outOfView(P, tr.x)) { tr.dead = true; break; }
          cam.project(tr.x, q); tr.pts.push(q[0], q[1]);
        }
        if (tr.pts.length > 2 * maxPts) tr.pts.splice(0, tr.pts.length - 2 * maxPts);
      });
      this.trails = this.trails.filter(function (tr) { return !tr.dead || tr.pts.length > 2; });
      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);
      ctx.save();
      if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }
      this.trails.forEach(function (tr, k) {
        const col = P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length];
        const pts = tr.pts, m = pts.length / 2;
        if (m < 2) return;
        ctx.strokeStyle = col; ctx.lineWidth = st.lineWidth; ctx.lineJoin = "round";
        const G = 10;
        for (let g = 0; g < G; g++) {
          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
          if (b <= a) continue;
          ctx.globalAlpha = st.alpha * (0.15 + 0.85 * (g + 1) / G);
          ctx.beginPath(); ctx.moveTo(pts[2 * a], pts[2 * a + 1]);
          for (let j = a + 1; j <= b; j++) ctx.lineTo(pts[2 * j], pts[2 * j + 1]);
          ctx.stroke();
        }
        ctx.globalAlpha = 1; ctx.fillStyle = col;
        if (!tr.dead) { ctx.beginPath(); ctx.arc(pts[pts.length - 2], pts[pts.length - 1], 3.5, 0, 7); ctx.fill(); }
      });
      ctx.restore();
      self.lastTrails = this.trails;
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down") return false;
      const u = P.cam.unproject(px, py), s = Float64Array.from(P.sim.init);
      s[P.cam.axes[0]] = u[0]; s[P.cam.axes[1]] = u[1];
      this.addTrail(P, s);
      return true;
    },
    legend: function (P) {
      const L = [["stable", P.theme.ink, "dot"], ["saddle", P.palette[0], "dot"], ["unstable", P.theme.ink, "ring"]];
      if (P.scene.view.nullclines !== false) L.unshift([P.sys.vars[P.cam.axes[0]] + "-nullcline", P.palette[1 % P.palette.length], "line"], [P.sys.vars[P.cam.axes[1]] + "-nullcline", P.palette[2 % P.palette.length], "dash"]);
      return L;
    },
    svg: function (P) {
      let out = "";
      (this.nullSegs || []).forEach(function (nc) {
        let d = "";
        for (let i = 0; i < nc.pts.length; i += 4) d += "M" + nc.pts[i].toFixed(1) + " " + nc.pts[i + 1].toFixed(1) + "L" + nc.pts[i + 2].toFixed(1) + " " + nc.pts[i + 3].toFixed(1);
        out += '<path d="' + d + '" fill="none" stroke="' + nc.color + '" stroke-width="1.6"' + (nc.dash ? ' stroke-dasharray="' + nc.dash + '"' : "") + "/>";
      });
      (this.trails || []).forEach(function (tr, k) { out += polyline(tr.pts, P.palette[0], P.scene.style.lineWidth, P.scene.style.alpha); });
      const q = [0, 0], z = Float64Array.from(P.sim.init), th = P.theme;
      (this.equilibria || []).forEach(function (e) {
        z[P.cam.axes[0]] = e.x[0]; z[P.cam.axes[1]] = e.x[1]; P.cam.project(z, q);
        out += '<circle cx="' + q[0].toFixed(1) + '" cy="' + q[1].toFixed(1) + '" r="5.5" fill="' + (e.stable ? th.ink : "none") + '" stroke="' + th.ink + '" stroke-width="2"/>';
      });
      return out;
    }
  };

  // ------------------------------------------------------------ sweep
  /* Branches of equilibria (or fixed points) against one parameter, found
     by Newton iteration at each of `cols` parameter values, then a single
     state driven by a slow triangular sweep of that parameter. */
  function branches(P, pi, vi, from, to, cols) {
    const sys = P.sys, p = Float64Array.from(P.sim.base), pts = [];
    const box = P.fullRanges.map(function (r) { return r.slice(); });
    let prev = [];
    for (let c = 0; c <= cols; c++) {
      p[pi] = from + (to - from) * c / cols;
      let eq = [];
      try { eq = DF.findEquilibria(sys, p, box, { seeds: 24, extra: prev }); } catch (e) { eq = []; }
      prev = eq.map(function (e) { return e.x; });
      eq.forEach(function (e) { pts.push({ p: p[pi], v: e.x[vi], stable: e.stable }); });
    }
    return pts;
  }

  V.sweep = {
    label: "Parameter sweep (hysteresis)",
    axes: true,
    init: function (P) {
      const v = P.scene.view, sys = P.sys;
      this.pi = Math.max(0, sys.params.findIndex(function (q) { return q.name === v.param; }));
      this.vi = Math.max(0, sys.vars.indexOf(v.var || sys.vars[0]));
      const q = sys.params[this.pi];
      this.from = v.from === undefined ? q.min : v.from; this.to = v.to === undefined ? q.max : v.to;
      this.pval = this.from; this.dir = 1; this.trail = [];
      this.cam = new DF.Camera([0, 1], [[this.from, this.to], P.fullRanges[this.vi]], { pad: 0.04 });
      this.cam.resize(P.w, P.h, P.inset);
      this.branchPts = null;
    },
    drawStatic: function (P) {
      this.cam.resize(P.w, P.h, P.inset);
      const ctx = P.ctx.base, cam = this.cam, th = P.theme;
      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [P.sys.params[this.pi].name, P.sys.vars[this.vi]], { grid: false });
      this.frameBox = f;
      if (P.scene.view.branches === false) return;
      if (!this.branchPts) this.branchPts = branches(P, this.pi, this.vi, this.from, this.to, Math.min(400, Math.round((f.R - f.L) / 2)));
      const q = [0, 0], x = [0, 0];
      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      this.branchPts.forEach(function (b) {
        x[0] = b.p; x[1] = b.v; cam.project(x, q);
        ctx.fillStyle = b.stable ? th.ink : th.muted;
        const r = b.stable ? 1.7 : 1.1;
        if (b.stable) { ctx.beginPath(); ctx.arc(q[0], q[1], r, 0, 7); ctx.fill(); }
        else { ctx.globalAlpha = 0.8; ctx.fillRect(q[0] - r, q[1] - r, 2 * r, 2 * r); ctx.globalAlpha = 1; }
      });
      ctx.restore();
    },
    frame: function (P) {
      const v = P.scene.view, sim = P.sim;
      if (!P.userParam) {
        const speed = (v.speed || 0.0006) * (this.to - this.from);
        this.pval += this.dir * speed;
        if (this.pval > this.to) { this.pval = this.to; this.dir = -1; }
        if (this.pval < this.from) { this.pval = this.from; this.dir = 1; }
        sim.base[this.pi] = this.pval; sim.updateParams();
        P.emit("param", { name: P.sys.params[this.pi].name, value: this.pval });
      } else this.pval = sim.base[this.pi];
      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();
      if (!sim.alive[0]) sim.setMember(0, sim.init);
      const q = [0, 0];
      this.cam.project([this.pval, sim.X[this.vi]], q);
      this.trail.push(q[0], q[1]);
      if (this.trail.length > 2 * (v.tail || 1400)) this.trail.splice(0, 2);
      const ctx = P.ctx.top, st = P.scene.style;
      ctx.clearRect(0, 0, P.w, P.h);
      ctx.save(); if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }
      const m = this.trail.length / 2, G = 16;
      ctx.strokeStyle = P.palette[1 % P.palette.length]; ctx.lineWidth = st.lineWidth;
      for (let g = 0; g < G; g++) {
        const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
        if (b <= a) continue;
        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * (g + 1) / G);
        ctx.beginPath(); ctx.moveTo(this.trail[2 * a], this.trail[2 * a + 1]);
        for (let j = a + 1; j <= b; j++) ctx.lineTo(this.trail[2 * j], this.trail[2 * j + 1]);
        ctx.stroke();
      }
      ctx.globalAlpha = 1; ctx.fillStyle = P.palette[0];
      ctx.beginPath(); ctx.arc(q[0], q[1], 6, 0, 7); ctx.fill();
      // Parameter marker on the axis.
      if (this.frameBox) { ctx.strokeStyle = P.palette[0]; ctx.globalAlpha = 0.5; ctx.beginPath(); ctx.moveTo(q[0], this.frameBox.T); ctx.lineTo(q[0], this.frameBox.B); ctx.stroke(); ctx.globalAlpha = 1; }
      ctx.restore();
    },
    pointer: function (P, kind, px) {
      if (kind !== "down" && kind !== "drag") return false;
      const u = this.cam.unproject(px, 0);
      P.userParam = true; P.sim.base[this.pi] = Math.min(this.to, Math.max(this.from, u[0])); P.sim.updateParams();
      P.emit("param", { name: P.sys.params[this.pi].name, value: P.sim.base[this.pi] });
      return true;
    },
    legend: function (P) { return [["stable branch", P.theme.ink, "dot"], ["unstable branch", P.theme.muted, "dot"], ["state under the sweep", P.palette[1 % P.palette.length], "line"]]; },
    svg: function (P) {
      let out = "";
      const q = [0, 0], cam = this.cam, th = P.theme;
      (this.branchPts || []).forEach(function (b) { cam.project([b.p, b.v], q); out += '<circle cx="' + q[0].toFixed(1) + '" cy="' + q[1].toFixed(1) + '" r="' + (b.stable ? 1.7 : 1.1) + '" fill="' + (b.stable ? th.ink : th.muted) + '"/>'; });
      out += polyline(this.trail, P.palette[1 % P.palette.length], P.scene.style.lineWidth, P.scene.style.alpha);
      return out;
    }
  };

  // ------------------------------------------------------------ orbit
  V.orbit = {
    label: "Bifurcation diagram",
    axes: true,
    init: function (P) {
      const v = P.scene.view, sys = P.sys;
      this.pi = Math.max(0, sys.params.findIndex(function (q) { return q.name === v.param; }));
      this.vi = Math.max(0, sys.vars.indexOf(v.var || sys.vars[0]));
      const q = sys.params[this.pi];
      this.from = v.from === undefined ? q.min : v.from; this.to = v.to === undefined ? q.max : v.to;
      this.cam = new DF.Camera([0, 1], [[this.from, this.to], P.fullRanges[this.vi]], { pad: 0.02 });
      this.cam.resize(P.w, P.h, P.inset);
      this.col = 0; this.state = Float64Array.from(P.sim.init); this.points = [];
      this.rk = DF.makeRK4(sys.vars.length);
      this.p = Float64Array.from(P.sim.base);
    },
    drawStatic: function (P) {
      this.cam.resize(P.w, P.h, P.inset);
      this.frameBox = DF.drawAxes(P.ctx.base, this.cam, P.scene.style.theme, [P.sys.params[this.pi].name, P.sys.vars[this.vi]], { grid: false });
      this.cols = Math.max(50, Math.round((this.frameBox.R - this.frameBox.L) / (P.scene.view.colWidth || 1.2)));
      // Existing points are redrawn after a resize.
      const ctx = P.ctx.trail, q = [0, 0], st = P.scene.style;
      ctx.clearRect(0, 0, P.w, P.h);
      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;
      for (let i = 0; i < this.points.length; i += 2) { this.cam.project([this.points[i], this.points[i + 1]], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2); }
      ctx.globalAlpha = 1;
    },
    frame: function (P) {
      if (!this.cols || this.col > this.cols) return;
      const v = P.scene.view, sys = P.sys, dim = sys.vars.length, h = P.sim.h, st = P.scene.style;
      const transient = v.transient || (sys.time === "discrete" ? 300 : Math.round(200 / h));
      const samples = v.samples || (sys.time === "discrete" ? 150 : Math.round(400 / h));
      const ctx = P.ctx.trail, q = [0, 0], budget = performance.now() + (v.budget || 12);
      const x = this.state, tmp = P.tmpDx, H = function (i) { return x[i]; };
      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;
      while (this.col <= this.cols && performance.now() < budget) {
        const pv = this.from + (this.to - this.from) * this.col / this.cols;
        this.p[this.pi] = pv;
        if (v.follow === false || !x.every(isFinite)) x.set(P.sim.init);
        let t = 0, prev2 = NaN, prev1 = NaN;
        const adv = function (self) { if (sys.time === "discrete") { sys.f(t, x, self.p, tmp, H); x.set(tmp); t += 1; } else { self.rk(sys.f, t, x, self.p, h, H); t += h; } };
        for (let s = 0; s < transient; s++) adv(this);
        for (let s = 0; s < samples; s++) {
          adv(this);
          const y = x[this.vi];
          if (!isFinite(y)) break;
          if (sys.time === "discrete") { this.points.push(pv, y); this.cam.project([pv, y], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2); }
          else if (prev1 > prev2 && prev1 >= y) {
            // Local maximum: vertex of the parabola through the last three samples.
            const den = prev2 - 2 * prev1 + y, peak = den !== 0 ? prev1 - (prev2 - y) * (prev2 - y) / (8 * den) : prev1;
            this.points.push(pv, peak); this.cam.project([pv, peak], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2);
          }
          prev2 = prev1; prev1 = y;
        }
        this.col++;
      }
      ctx.globalAlpha = 1;
      const top = P.ctx.top; top.clearRect(0, 0, P.w, P.h);
      if (this.col <= this.cols && this.frameBox) {
        const px = this.frameBox.L + (this.frameBox.R - this.frameBox.L) * this.col / this.cols;
        top.strokeStyle = P.palette[1 % P.palette.length]; top.globalAlpha = 0.6; top.beginPath(); top.moveTo(px, this.frameBox.T); top.lineTo(px, this.frameBox.B); top.stroke(); top.globalAlpha = 1;
      }
    },
    legend: function (P) { return [[P.sys.time === "discrete" ? "iterates after a transient" : "local maxima after a transient", P.palette[0], "dot"]]; }
  };

  // ------------------------------------------------------------ density
  V.density = {
    label: "Ensemble density",
    axes: true,
    init: function (P) {
      const v = P.scene.view;
      this.mode = v.mode || (P.sys.vars.length === 1 ? "carpet" : "map");
      this.vi = Math.max(0, P.sys.vars.indexOf(v.var || P.sys.vars[P.cam.axes[1] === undefined ? 0 : P.cam.axes[1]]));
      if (P.sys.vars.length === 1) this.vi = 0;
      this.grid = null; this.img = null;
      const warm = v.warmup || 0;
      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();
    },
    drawStatic: function (P) {
      const st = P.scene.style;
      if (this.mode === "carpet") {
        const span = P.scene.view.window || 60;
        this.cam = new DF.Camera([0, 1], [[-span, 0], P.fullRanges[this.vi]], { pad: 0.0 });
        this.cam.resize(P.w, P.h, P.inset);
        this.frameBox = DF.drawAxes(P.ctx.base, this.cam, st.theme, ["t - t now", P.sys.vars[this.vi]], {});
      } else {
        this.frameBox = DF.drawAxes(P.ctx.base, P.cam, st.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]], {});
      }
      const f = this.frameBox;
      this.nx = Math.max(10, Math.round((f.R - f.L) / (P.scene.view.cell || 3)));
      this.ny = Math.max(10, Math.round((f.B - f.T) / (P.scene.view.cell || 3)));
      this.grid = new Float32Array(this.nx * this.ny);
      if (typeof document !== "undefined") {
        this.off = document.createElement("canvas"); this.off.width = this.nx; this.off.height = this.ny;
        this.octx = this.off.getContext("2d"); this.img = this.octx.createImageData(this.nx, this.ny);
      }
    },
    frame: function (P) {
      const sim = P.sim, dim = P.sys.vars.length, v = P.scene.view, st = P.scene.style;
      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();
      for (let k = 0; k < sim.n; k++) if (!sim.alive[k]) sim.setMember(k, spawnState(P, new Float64Array(dim)));
      if (!this.grid || !this.img) return;
      const nx = this.nx, ny = this.ny, g = this.grid;
      if (this.mode === "carpet") {
        // shift one column left per `shiftEvery` frames, fill the last column
        const r = P.fullRanges[this.vi];
        for (let j = 0; j < ny; j++) { for (let i = 0; i < nx - 1; i++) g[j * nx + i] = g[j * nx + i + 1]; g[j * nx + nx - 1] = 0; }
        for (let k = 0; k < sim.n; k++) {
          const y = sim.X[k * dim + this.vi], j = Math.floor((r[1] - y) / (r[1] - r[0]) * ny);
          if (j >= 0 && j < ny) g[j * nx + nx - 1] += 1;
        }
        let mx = 0; for (let j = 0; j < ny; j++) mx = Math.max(mx, g[j * nx + nx - 1]);
        for (let j = 0; j < ny; j++) g[j * nx + nx - 1] /= (mx || 1);
      } else {
        const decay = v.decay === undefined ? 0.85 : v.decay;
        for (let i = 0; i < g.length; i++) g[i] *= decay;
        const cam = P.cam, a0 = cam.axes[0], a1 = cam.axes[1], R0 = cam.ranges[0], R1 = cam.ranges[1];
        for (let k = 0; k < sim.n; k++) {
          const u = (sim.X[k * dim + a0] - R0[0]) / (R0[1] - R0[0]), w = (R1[1] - sim.X[k * dim + a1]) / (R1[1] - R1[0]);
          const i = Math.floor(u * nx), j = Math.floor(w * ny);
          if (i >= 0 && i < nx && j >= 0 && j < ny) g[j * nx + i] += 1;
        }
      }
      let mx = 0; for (let i = 0; i < g.length; i++) mx = Math.max(mx, g[i]);
      const d = this.img.data, logm = Math.log1p(mx);
      for (let i = 0; i < g.length; i++) {
        const u = mx > 0 ? (v.log === false ? g[i] / mx : Math.log1p(g[i]) / logm) : 0;
        const c = DF.rampRGB(st.ramp, u);
        d[4 * i] = c[0]; d[4 * i + 1] = c[1]; d[4 * i + 2] = c[2]; d[4 * i + 3] = u < 0.004 ? 0 : Math.round(255 * Math.min(1, 0.15 + 1.2 * u));
      }
      this.octx.putImageData(this.img, 0, 0);
      const ctx = P.ctx.top, f = this.frameBox;
      ctx.clearRect(0, 0, P.w, P.h);
      ctx.imageSmoothingEnabled = v.smooth !== false;
      ctx.drawImage(this.off, f.L, f.T, f.R - f.L, f.B - f.T);
    }
  };

  // ------------------------------------------------------------ strobe
  V.strobe = {
    label: "Stroboscopic / Poincare section",
    axes: true,
    init: function (P) {
      const v = P.scene.view;
      const per = P.sim.perturbations.find(function (q) { return q.kind === "periodic" && q.enabled !== false; });
      this.period = v.period || (per ? per.period : 2 * Math.PI);
      this.next = (v.phase || 0) * this.period;
      while (this.next <= P.sim.t) this.next += this.period;
      this.prevSide = null; this.prev = null;
      this.count = 0; this.transient = v.transient === undefined ? 20 : v.transient;
      this.points = [];
    },
    drawStatic: function (P) {
      this.frameBox = DF.drawAxes(P.ctx.base, P.cam, P.scene.style.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]], {});
      const ctx = P.ctx.trail, q = [0, 0], st = P.scene.style, z = Float64Array.from(P.sim.init);
      ctx.clearRect(0, 0, P.w, P.h); ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;
      for (let i = 0; i < this.points.length; i += 2) { z[P.cam.axes[0]] = this.points[i]; z[P.cam.axes[1]] = this.points[i + 1]; P.cam.project(z, q); ctx.fillRect(q[0] - st.pointSize / 2, q[1] - st.pointSize / 2, st.pointSize, st.pointSize); }
      ctx.globalAlpha = 1;
    },
    frame: function (P) {
      const sim = P.sim, dim = P.sys.vars.length, v = P.scene.view, st = P.scene.style, ctx = P.ctx.trail, q = [0, 0];
      const a0 = P.cam.axes[0], a1 = P.cam.axes[1];
      if (st.fade > 0) P.fade(st.fade);
      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha; ctx.globalCompositeOperation = P.theme.blend;
      const plot = (function (self) { return function (x) { self.points.push(x[a0], x[a1]); if (self.points.length > 400000) self.points.splice(0, 2); P.cam.project(x, q); ctx.fillRect(q[0] - st.pointSize / 2, q[1] - st.pointSize / 2, st.pointSize, st.pointSize); }; })(this);
      const x = P.tmpX;
      for (let s = 0; s < P.scene.stepsPerFrame; s++) {
        if (v.mode === "section") {
          const si = Math.max(0, P.sys.vars.indexOf(v.sectionVar || P.sys.vars[dim - 1])), c = v.sectionValue || 0;
          if (!this.prev) this.prev = Float64Array.from(sim.X);
          this.prev.set(sim.X);
          sim.step();
          for (let k = 0; k < sim.n; k++) {
            const y0 = this.prev[k * dim + si] - c, y1 = sim.X[k * dim + si] - c;
            if (y0 < 0 && y1 >= 0) {
              const f = y0 / (y0 - y1);
              for (let i = 0; i < dim; i++) x[i] = this.prev[k * dim + i] + f * (sim.X[k * dim + i] - this.prev[k * dim + i]);
              if (sim.t > this.transient) plot(x);
            }
          }
        } else {
          sim.step();
          if (sim.t + 1e-12 >= this.next) {
            this.count++;
            if (this.count > this.transient) for (let k = 0; k < sim.n; k++) { if (sim.alive[k]) plot(sim.member(k, x)); }
            this.next += this.period;
          }
        }
      }
      ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";
    },
    legend: function (P) { return [[P.scene.view.mode === "section" ? "upward crossings of the section" : "state at t = t0 + kT, T = " + DF.fmt(this.period), P.palette[0], "dot"]]; }
  };

  // ------------------------------------------------------------ cobweb
  V.cobweb = {
    label: "Cobweb (1D maps)",
    axes: true,
    init: function (P) {
      this.x = P.sim.init[0]; this.path = []; this.vi = 0;
      const r = P.fullRanges[0];
      this.cam = new DF.Camera([0, 1], [r, r], { pad: 0.03 });
      this.cam.resize(P.w, P.h, P.inset);
    },
    drawStatic: function (P) {
      this.cam.resize(P.w, P.h, P.inset);
      const ctx = P.ctx.base, cam = this.cam, th = P.theme, v = P.sys.vars[0];
      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [v + "ₙ", v + "ₙ₊₁"], {});
      this.frameBox = f;
      const r = cam.ranges[0], q = [0, 0], out = new Float64Array(P.sys.vars.length), x = Float64Array.from(P.sim.init);
      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      ctx.strokeStyle = th.muted; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
      cam.project([r[0], r[0]], q); ctx.beginPath(); ctx.moveTo(q[0], q[1]); cam.project([r[1], r[1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();
      ctx.setLineDash([]); ctx.strokeStyle = P.palette[2 % P.palette.length]; ctx.lineWidth = 2; ctx.beginPath();
      for (let i = 0; i <= 600; i++) {
        x[0] = r[0] + (r[1] - r[0]) * i / 600; P.sys.f(P.sim.t, x, P.sim.p, out);
        cam.project([x[0], out[0]], q); if (i) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]);
      }
      ctx.stroke(); ctx.restore();
    },
    frame: function (P) {
      const every = Math.max(1, Math.round(8 / Math.max(1, P.scene.stepsPerFrame)));
      if ((P.frameCount % every) !== 0) return;
      const out = P.tmpDx, x = P.tmpX;
      x[0] = this.x; P.sys.f(P.sim.t, x, P.sim.p, out);
      const y = out[0];
      if (!isFinite(y)) { this.x = P.sim.init[0]; this.path = []; return; }
      if (!this.path.length) this.path.push(this.x, this.cam.ranges[1][0] < 0 && this.cam.ranges[1][1] > 0 ? 0 : this.cam.ranges[1][0]);
      this.path.push(this.x, y, y, y);
      if (this.path.length > 2 * (P.scene.view.tail || 120)) this.path.splice(0, 4);
      this.x = y; P.sim.t += 1;
      const ctx = P.ctx.top, st = P.scene.style, q = [0, 0];
      ctx.clearRect(0, 0, P.w, P.h);
      ctx.save(); const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      const m = this.path.length / 2;
      ctx.strokeStyle = P.palette[0]; ctx.lineWidth = st.lineWidth;
      for (let j = 1; j < m; j++) {
        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * j / m);
        ctx.beginPath(); this.cam.project([this.path[2 * j - 2], this.path[2 * j - 1]], q); ctx.moveTo(q[0], q[1]);
        this.cam.project([this.path[2 * j], this.path[2 * j + 1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();
      }
      ctx.globalAlpha = 1; ctx.fillStyle = P.palette[0]; ctx.beginPath(); ctx.arc(q[0], q[1], 4, 0, 7); ctx.fill();
      ctx.restore();
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down") return false;
      this.x = this.cam.unproject(px, py)[0]; this.path = [];
      return true;
    }
  };

  // Which views suit a system.
  // Scalar systems are shown against time; state-space views need two variables.
  function viewsFor(sys) {
    const multi = sys.vars.length >= 2;
    const out = multi ? ["flow", "trajectory", "timeseries", "density"] : ["timeseries", "density"];
    if (sys.time === "continuous" && multi) out.push("phase");
    if (sys.params.length) out.push("sweep", "orbit");
    if (multi) out.push("strobe");
    if (sys.time === "discrete" && !multi) out.push("cobweb");
    return out;
  }

  DF.VIEWS = V;
  DF.viewsFor = viewsFor;
  DF.contour = contour;
})(globalThis.DynFlow = globalThis.DynFlow || {});
