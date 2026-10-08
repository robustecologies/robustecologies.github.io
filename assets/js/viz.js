// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Canvas renderers used by the layouts. Every animation integrates a system
   from flows.js; nothing is a precomputed video. Animations stop when the
   canvas leaves the viewport and draw a single static frame when the reader
   asks for reduced motion. */
(function (root) {
  "use strict";
  const F = root.Flows;
  const reduced = root.matchMedia && root.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function css(name, el) {
    return getComputedStyle(el || document.documentElement).getPropertyValue(name).trim();
  }

  // Resize a canvas to its CSS box at device resolution; returns the 2D context.
  function fit(canvas) {
    const dpr = Math.min(root.devicePixelRatio || 1, 2);
    const r = canvas.getBoundingClientRect();
    canvas.width = Math.max(1, Math.round(r.width * dpr));
    canvas.height = Math.max(1, Math.round(r.height * dpr));
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx: ctx, w: r.width, h: r.height };
  }

  // Run step(dtMs) on animation frames while the element is visible.
  function loop(el, step) {
    let visible = true, last = 0, id = 0;
    function frame(t) {
      const dt = last ? Math.min(t - last, 50) : 16;
      last = t;
      step(dt);
      if (visible) id = requestAnimationFrame(frame);
    }
    if ("IntersectionObserver" in root) {
      new IntersectionObserver(function (entries) {
        const was = visible;
        visible = entries[0].isIntersecting;
        if (visible && !was) { last = 0; id = requestAnimationFrame(frame); }
        if (!visible) cancelAnimationFrame(id);
      }).observe(el);
    }
    id = requestAnimationFrame(frame);
  }

  function hexToRgb(hex) {
    const v = parseInt(hex.replace("#", ""), 16);
    return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
  }

  // Continuous colour by speed: the logarithm of the on-screen displacement
  // of a particle per frame, placed within a running envelope of the whole
  // ensemble (at most five e-folds wide, relaxing slowly inward), is read
  // off a gradient through the palette in 64 levels.
  function SpeedShade(stops) {
    const rgb = stops.map(hexToRgb), cache = new Array(64), env = [0, 0];
    let init = false;
    for (let i = 0; i < 64; i++) {
      const u = i / 63 * (rgb.length - 1), j = Math.min(rgb.length - 2, Math.floor(u)), f = u - j;
      cache[i] = [0, 1, 2].map(function (c) { return Math.round(rgb[j][c] + f * (rgb[j + 1][c] - rgb[j][c])); }).join(",");
    }
    return function (d) {
      const l = Math.log(Math.max(d, 1e-6));
      if (!init) { env[0] = l - 1; env[1] = l; init = true; }
      env[1] = Math.max(l, env[1] - 2e-5); env[0] = Math.max(env[1] - 5, Math.min(l, env[0] + 2e-5));
      const u = Math.min(1, Math.max(0, (l - env[0]) / Math.max(1e-9, env[1] - env[0])));
      return cache[Math.round(u * 63)];
    };
  }
  const GRADIENT = ["#5b2a86", "#A52C60", "#EE6A24", "#FB9E07", "#F6D645"];

  // ---------------------------------------------------------------------------
  // Heteroclinic field: particles follow the May-Leonard flow projected on the
  // simplex. They are born near the interior equilibrium, spiral out, and then
  // visit the three single-species saddles with residence times that grow by
  // the factor (beta - 1)/(1 - alpha) = 1.5 per passage (checks/flows_check.mjs).
  // ---------------------------------------------------------------------------
  function heteroclinicField(canvas, opts) {
    opts = Object.assign({
      n: 1400, speed: 0.22, fade: 0.06, scale: 1.0, dx: 0, dy: 0, rotate: 0,
      gradient: null, background: null, lineWidth: 1.1,
      interactive: true, life: [60, 260], nearEq: 0.25
    }, opts || {});
    const a = F.ML_ALPHA, b = F.ML_BETA;
    const eq = 1 / (1 + a + b);
    let g = fit(canvas);
    const shade = SpeedShade(opts.gradient || GRADIENT);
    const P = [];

    function place(p, x) {
      p.x0 = x[0]; p.x1 = x[1]; p.x2 = x[2];
      p.age = 0;
      p.life = opts.life[0] + Math.random() * (opts.life[1] - opts.life[0]);
      p.px = null;
    }
    // Births: a fraction near the interior equilibrium, which repels slowly
    // (growth rate about 0.017 per unit time, so these orbits spiral), and the
    // rest uniform on the simplex, which reach the heteroclinic cycle quickly.
    function spawn(p) {
      if (Math.random() < opts.nearEq) {
        const e = 0.02 + 0.08 * Math.random();
        place(p, [eq + e * (Math.random() - 0.5), eq + e * (Math.random() - 0.5), eq + e * (Math.random() - 0.5)]);
      } else {
        const u = -Math.log(Math.random()), v = -Math.log(Math.random()), w = -Math.log(Math.random()), t = u + v + w;
        place(p, [0.01 + 0.97 * u / t, 0.01 + 0.97 * v / t, 0.01 + 0.97 * w / t]);
      }
      p.age = -Math.random() * 20;
    }

    function f(x0, x1, x2, out) {
      out[0] = x0 * (1 - x0 - a * x1 - b * x2);
      out[1] = x1 * (1 - x1 - a * x2 - b * x0);
      out[2] = x2 * (1 - x2 - a * x0 - b * x1);
    }
    const k1 = [0, 0, 0], k2 = [0, 0, 0], k3 = [0, 0, 0], k4 = [0, 0, 0];
    function step(p, h) {
      f(p.x0, p.x1, p.x2, k1);
      f(p.x0 + 0.5 * h * k1[0], p.x1 + 0.5 * h * k1[1], p.x2 + 0.5 * h * k1[2], k2);
      f(p.x0 + 0.5 * h * k2[0], p.x1 + 0.5 * h * k2[1], p.x2 + 0.5 * h * k2[2], k3);
      f(p.x0 + h * k3[0], p.x1 + h * k3[1], p.x2 + h * k3[2], k4);
      p.x0 += (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
      p.x1 += (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
      p.x2 += (h / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    }
    for (let i = 0; i < opts.n; i++) { const p = {}; spawn(p); p.age = Math.random() * p.life; P.push(p); }
    // Advance the initial population so the first frame already shows the cycle.
    for (let i = 0; i < P.length; i++) { const p = P[i]; for (let k = 0; k < p.age / 0.25; k++) step(p, 0.25); }
    const cosR = Math.cos(opts.rotate), sinR = Math.sin(opts.rotate);
    function project(p) {
      const q = F.simplexToPlane([p.x0, p.x1, p.x2]);
      const u = q[0] - 0.5, v = q[1] - 0.5;
      const s = Math.max(g.w, g.h * 1.15) * opts.scale;
      return [g.w / 2 + opts.dx * g.w + s * (cosR * u - sinR * v), g.h / 2 + opts.dy * g.h + s * (sinR * u + cosR * v)];
    }
    // Inverse of project for a click: barycentric weights of the point.
    function unproject(X, Y) {
      const s = Math.max(g.w, g.h * 1.15) * opts.scale;
      const u0 = (X - g.w / 2 - opts.dx * g.w) / s, v0 = (Y - g.h / 2 - opts.dy * g.h) / s;
      const u = cosR * u0 + sinR * v0 + 0.5, v = -sinR * u0 + cosR * v0 + 0.5;
      const V = [[0.5, 0.06], [0.04, 0.86], [0.96, 0.86]];
      const det = (V[1][1] - V[2][1]) * (V[0][0] - V[2][0]) + (V[2][0] - V[1][0]) * (V[0][1] - V[2][1]);
      const w0 = ((V[1][1] - V[2][1]) * (u - V[2][0]) + (V[2][0] - V[1][0]) * (v - V[2][1])) / det;
      const w1 = ((V[2][1] - V[0][1]) * (u - V[2][0]) + (V[0][0] - V[2][0]) * (v - V[2][1])) / det;
      const w = [w0, w1, 1 - w0 - w1];
      if (w.some(function (c) { return c <= 0.002; })) return null;
      return w;
    }

    function draw(h) {
      const ctx = g.ctx;
      ctx.globalCompositeOperation = "destination-out";
      ctx.fillStyle = "rgba(0,0,0," + opts.fade + ")";
      ctx.fillRect(0, 0, g.w, g.h);
      ctx.globalCompositeOperation = "lighter";
      ctx.lineWidth = opts.lineWidth;
      for (let i = 0; i < P.length; i++) {
        const p = P[i];
        p.age += h;
        if (p.age < 0) continue;
        if (p.age > p.life || !(p.x0 + p.x1 + p.x2 > 0)) { spawn(p); continue; }
        const before = p.px || project(p);
        step(p, h);
        const now = project(p);
        const c = shade(Math.hypot(now[0] - before[0], now[1] - before[1]));
        const fadeIn = Math.min(1, p.age / 8), fadeOut = Math.min(1, (p.life - p.age) / 10);
        ctx.strokeStyle = "rgba(" + c + "," + (0.55 * fadeIn * fadeOut).toFixed(3) + ")";
        ctx.beginPath(); ctx.moveTo(before[0], before[1]); ctx.lineTo(now[0], now[1]); ctx.stroke();
        p.px = now;
      }
      ctx.globalCompositeOperation = "source-over";
    }

    if (opts.interactive) {
      canvas.addEventListener("pointerdown", function (ev) {
        const r = canvas.getBoundingClientRect();
        const w = unproject(ev.clientX - r.left, ev.clientY - r.top);
        if (!w) return;
        for (let i = 0; i < 60; i++) {
          const p = P[(Math.random() * P.length) | 0];
          const j = 0.01;
          place(p, [w[0] + j * Math.random(), w[1] + j * Math.random(), w[2] + j * Math.random()]);
        }
      });
    }
    root.addEventListener("resize", function () { g = fit(canvas); P.forEach(function (p) { p.px = null; }); });

    if (reduced) { for (let k = 0; k < 60; k++) draw(0.25); return; }
    loop(canvas, function (ms) { draw(opts.speed * ms / 16); });
  }

  // ---------------------------------------------------------------------------
  // Generic particle field for any flow of shared/portraits.js. A system gives
  // its vector field f(x, out), a birth law spawn(), a plane view(x) -> [u, v]
  // with bounds box = [u0, u1, v0, v1], the
  // model time per frame (speed) and the largest RK4 step (hmax). Particles
  // live a random number of frames and are then reborn, so the field keeps
  // showing transients as well as the attractor.
  // ---------------------------------------------------------------------------
  function portraitField(canvas, sys, opts) {
    opts = Object.assign({
      n: 1600, fade: 0.06, gradient: null, lineWidth: 1.1,
      alpha: 0.55, margin: 0.06, interactive: true
    }, opts || {});
    let g = fit(canvas);
    const shade = SpeedShade(opts.gradient || GRADIENT), dim = sys.dim, life = sys.life || [90, 320];
    const k1 = new Array(dim), k2 = new Array(dim), k3 = new Array(dim), k4 = new Array(dim), y = new Array(dim);
    function rk4(x, h) {
      sys.f(x, k1);
      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k1[i];
      sys.f(y, k2);
      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k2[i];
      sys.f(y, k3);
      for (let i = 0; i < dim; i++) y[i] = x[i] + h * k3[i];
      sys.f(y, k4);
      for (let i = 0; i < dim; i++) x[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    function advance(x, T) {
      const m = Math.max(1, Math.ceil(T / sys.hmax)), h = T / m;
      for (let j = 0; j < m; j++) rk4(x, h);
    }
    // Plane to canvas: the box is stretched to the canvas less a margin.
    const B = sys.box;
    function toCanvas(q) {
      const mx = opts.margin * g.w, my = opts.margin * g.h;
      return [mx + ((q[0] - B[0]) / (B[1] - B[0])) * (g.w - 2 * mx), g.h - my - ((q[1] - B[2]) / (B[3] - B[2])) * (g.h - 2 * my)];
    }
    function fromCanvas(X, Y) {
      const mx = opts.margin * g.w, my = opts.margin * g.h;
      return [B[0] + ((X - mx) / (g.w - 2 * mx)) * (B[1] - B[0]), B[2] + ((g.h - my - Y) / (g.h - 2 * my)) * (B[3] - B[2])];
    }
    const P = [];
    function born(p, x) {
      p.x = x; p.age = 0; p.px = null;
      p.life = life[0] + Math.random() * (life[1] - life[0]);
    }
    for (let i = 0; i < opts.n; i++) {
      const p = {};
      born(p, sys.spawn());
      p.age = Math.random() * p.life;
      if (sys.warm) advance(p.x, p.age * sys.speed);
      P.push(p);
    }
    function draw(scale) {
      const ctx = g.ctx;
      ctx.globalCompositeOperation = "destination-out";
      ctx.fillStyle = "rgba(0,0,0," + opts.fade + ")";
      ctx.fillRect(0, 0, g.w, g.h);
      ctx.globalCompositeOperation = "lighter";
      ctx.lineWidth = opts.lineWidth;
      const T = sys.speed * scale;
      for (let i = 0; i < P.length; i++) {
        const p = P[i];
        p.age += scale;
        if (p.age > p.life || !p.x.every(isFinite)) { born(p, sys.spawn()); continue; }
        const before = p.px || toCanvas(sys.view(p.x));
        advance(p.x, T);
        const now = toCanvas(sys.view(p.x));
        const c = shade(Math.hypot(now[0] - before[0], now[1] - before[1]));
        const a = opts.alpha * Math.min(1, p.age / 8) * Math.min(1, (p.life - p.age) / 10);
        ctx.strokeStyle = "rgba(" + c + "," + a.toFixed(3) + ")";
        ctx.beginPath(); ctx.moveTo(before[0], before[1]); ctx.lineTo(now[0], now[1]); ctx.stroke();
        p.px = now;
      }
      ctx.globalCompositeOperation = "source-over";
    }
    if (opts.interactive && sys.seed) {
      canvas.addEventListener("pointerdown", function (ev) {
        const r = canvas.getBoundingClientRect();
        const q = fromCanvas(ev.clientX - r.left, ev.clientY - r.top);
        for (let i = 0; i < 60; i++) born(P[(Math.random() * P.length) | 0], sys.seed(q));
      });
    } else canvas.style.cursor = "default";
    root.addEventListener("resize", function () { g = fit(canvas); P.forEach(function (p) { p.px = null; }); });
    if (reduced) { for (let k = 0; k < 60; k++) draw(1); return; }
    loop(canvas, function (ms) { draw(ms / 16); });
  }

  // ---------------------------------------------------------------------------
  // Static phase portrait of the same flow, drawn as a few long orbits in ink.
  // ---------------------------------------------------------------------------
  function heteroclinicStatic(canvas, opts) {
    opts = Object.assign({ ink: css("--brand"), accent: "#EE6A24", orbits: 7 }, opts || {});
    const g = fit(canvas), ctx = g.ctx;
    const a = F.ML_ALPHA, b = F.ML_BETA, f = F.mayLeonard(a, b), eq = 1 / (1 + a + b);
    const s = Math.min(g.w, g.h * 1.1) * 0.9;
    function P(x) { const q = F.simplexToPlane(x); return [g.w / 2 + s * (q[0] - 0.5), g.h / 2 + s * (q[1] - 0.46)]; }
    const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(P);
    ctx.strokeStyle = opts.ink; ctx.globalAlpha = 0.25; ctx.lineWidth = 0.8;
    ctx.beginPath(); ctx.moveTo(V[0][0], V[0][1]); ctx.lineTo(V[1][0], V[1][1]); ctx.lineTo(V[2][0], V[2][1]); ctx.closePath(); ctx.stroke();
    ctx.globalAlpha = 1;
    for (let o = 0; o < opts.orbits; o++) {
      const e = 0.01 + 0.03 * o;
      let x = [eq + e * Math.cos(o * 0.9), eq + e * Math.sin(o * 0.9), eq];
      ctx.strokeStyle = o === 0 ? opts.accent : opts.ink;
      ctx.globalAlpha = o === 0 ? 0.95 : 0.35;
      ctx.lineWidth = o === 0 ? 1.4 : 0.7;
      ctx.beginPath();
      let q = P(x); ctx.moveTo(q[0], q[1]);
      for (let k = 0; k < 20000; k++) { x = F.rk4(f, x, 0.03); q = P(x); ctx.lineTo(q[0], q[1]); }
      ctx.stroke();
    }
    ctx.globalAlpha = 1; ctx.fillStyle = opts.ink;
    ctx.font = "italic 12px " + css("--font-brand");
    [["e₁", 0, -8], ["e₂", -14, 14], ["e₃", 6, 14]].forEach(function (l, i) { ctx.fillText(l[0], V[i][0] + l[1], V[i][1] + l[2]); });
  }

  // ---------------------------------------------------------------------------
  // Small live traces for the three goals.
  // ---------------------------------------------------------------------------
  function trace(canvas, kind, opts) {
    opts = Object.assign({ colors: ["#EE6A24", "#FB9E07", "#A52C60"], axis: "rgba(127,127,127,0.35)" }, opts || {});
    let g = fit(canvas);
    const N = 240;
    const hist = [];
    let x, f, h, t = 0, kickT = 0;
    if (kind === "transiency") { f = F.mayLeonard(F.ML_ALPHA, F.ML_BETA); x = [0.6, 0.3, 0.1]; h = 0.15; }
    else if (kind === "feedback") { f = F.lotkaVolterra(1.0, 0.5, 0.8, 0.4); x = [3, 1]; h = 0.03; }
    else { f = function (z) { return [0.8 * z[0] * (1 - z[0])]; }; x = [1]; h = 0.05; }

    function advance() {
      if (kind === "transiency") {
        x = F.rk4(f, x, h);
        // A fresh start near coexistence once the cycle has become too slow to watch.
        t += h; if (t > 260) { x = [0.6, 0.3, 0.1]; t = 0; }
        hist.push(x.slice());
      } else if (kind === "feedback") {
        x = F.rk4(f, x, h); hist.push(x.slice());
      } else {
        x = F.rk4(f, x, h); kickT += h;
        if (kickT > 6) { x = [Math.max(0.05, x[0] - (0.35 + 0.5 * Math.random()))]; kickT = 0; }
        hist.push(x.slice());
      }
      if (hist.length > N) hist.shift();
    }
    for (let k = 0; k < N; k++) advance();

    function draw() {
      const ctx = g.ctx, w = g.w, H = g.h, pad = 6;
      ctx.clearRect(0, 0, w, H);
      ctx.strokeStyle = opts.axis; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(pad, H - pad); ctx.lineTo(w - pad, H - pad); ctx.stroke();
      ctx.lineWidth = 1.6; ctx.lineJoin = "round";
      if (kind === "feedback") {
        ctx.beginPath();
        hist.forEach(function (z, i) {
          const X = pad + (z[0] / 6) * (w - 2 * pad), Y = H - pad - (z[1] / 5) * (H - 2 * pad);
          if (i) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
        });
        ctx.strokeStyle = opts.colors[0]; ctx.stroke();
        const z = hist[hist.length - 1];
        ctx.fillStyle = opts.colors[1];
        ctx.beginPath(); ctx.arc(pad + (z[0] / 6) * (w - 2 * pad), H - pad - (z[1] / 5) * (H - 2 * pad), 3, 0, 7); ctx.fill();
        return;
      }
      const series = kind === "transiency" ? 3 : 1;
      for (let s = 0; s < series; s++) {
        ctx.beginPath();
        hist.forEach(function (z, i) {
          const X = pad + (i / (N - 1)) * (w - 2 * pad), Y = H - pad - z[s] * (H - 2 * pad) * 0.92;
          if (i) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
        });
        ctx.strokeStyle = opts.colors[s]; ctx.stroke();
      }
    }
    root.addEventListener("resize", function () { g = fit(canvas); });
    if (reduced) { draw(); return; }
    let acc = 0;
    loop(canvas, function (ms) { acc += ms; while (acc > 33) { advance(); acc -= 33; } draw(); });
  }

  // ---------------------------------------------------------------------------
  // Cusp instrument: bifurcation diagram of dx/dt = r + x - x^3, a state that
  // relaxes under the flow, a slider for r and an automatic slow sweep that
  // traces the hysteresis loop.
  // ---------------------------------------------------------------------------
  function cuspInstrument(host, opts) {
    opts = Object.assign({ a: 1, ink: null, grid: "rgba(127,127,127,0.25)" }, opts || {});
    const a = opts.a, folds = F.cuspFolds(a);
    host.innerHTML =
      '<canvas class="cusp-canvas" aria-label="Bifurcation diagram of the cusp normal form"></canvas>' +
      '<div class="cusp-controls">' +
      '<label>Control parameter r <output>0.00</output></label>' +
      '<input type="range" min="-0.8" max="0.8" step="0.005" value="-0.6">' +
      '<button type="button" class="cusp-sweep" aria-pressed="true">Sweep</button>' +
      '<button type="button" class="cusp-kick">Perturb</button>' +
      "</div>";
    const canvas = host.querySelector("canvas"), slider = host.querySelector("input"), out = host.querySelector("output");
    const sweepBtn = host.querySelector(".cusp-sweep"), kickBtn = host.querySelector(".cusp-kick");
    let g = fit(canvas);
    let r = -0.6, x = F.cuspEquilibria(-0.6, a)[0], sweeping = !reduced, dir = 1;
    const trail = [];
    slider.addEventListener("input", function () { sweeping = false; sweepBtn.setAttribute("aria-pressed", "false"); r = +slider.value; });
    sweepBtn.addEventListener("click", function () { sweeping = !sweeping; sweepBtn.setAttribute("aria-pressed", String(sweeping)); });
    kickBtn.addEventListener("click", function () { x += x > 0 ? -1.1 : 1.1; });

    const R = [-0.8, 0.8], X = [-1.5, 1.5];
    function sx(v) { return 44 + ((v - R[0]) / (R[1] - R[0])) * (g.w - 60); }
    function sy(v) { return 12 + ((X[1] - v) / (X[1] - X[0])) * (g.h - 44); }

    function draw() {
      const ctx = g.ctx, ink = opts.ink || css("--ink", host), muted = css("--muted", host);
      ctx.clearRect(0, 0, g.w, g.h);
      ctx.strokeStyle = opts.grid; ctx.lineWidth = 1;
      folds.forEach(function (fr) { ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(sx(fr), sy(X[1])); ctx.lineTo(sx(fr), sy(X[0])); ctx.stroke(); });
      ctx.setLineDash([]);
      ctx.beginPath(); ctx.moveTo(sx(R[0]), sy(0)); ctx.lineTo(sx(R[1]), sy(0)); ctx.stroke();
      // Equilibrium curve, parametrised by x: r = x^3 - a x. Stable where 3x^2 > a.
      const xc = Math.sqrt(a / 3);
      [[X[0], -xc, true], [-xc, xc, false], [xc, X[1], true]].forEach(function (seg) {
        ctx.beginPath();
        for (let k = 0; k <= 200; k++) {
          const xv = seg[0] + (k / 200) * (seg[1] - seg[0]), rv = xv * xv * xv - a * xv;
          if (k) ctx.lineTo(sx(rv), sy(xv)); else ctx.moveTo(sx(rv), sy(xv));
        }
        ctx.setLineDash(seg[2] ? [] : [5, 5]);
        ctx.strokeStyle = seg[2] ? ink : muted; ctx.lineWidth = seg[2] ? 2 : 1.2; ctx.stroke();
      });
      ctx.setLineDash([]);
      ctx.strokeStyle = "#FB9E07"; ctx.lineWidth = 1.4; ctx.globalAlpha = 0.9;
      ctx.beginPath();
      trail.forEach(function (p, i) { if (i) ctx.lineTo(sx(p[0]), sy(p[1])); else ctx.moveTo(sx(p[0]), sy(p[1])); });
      ctx.stroke(); ctx.globalAlpha = 1;
      ctx.fillStyle = "#EE6A24";
      ctx.beginPath(); ctx.arc(sx(r), sy(x), 6, 0, 7); ctx.fill();
      ctx.fillStyle = muted; ctx.font = "12px " + css("--font-ui");
      ctx.fillText("r", g.w - 16, sy(0) - 6);
      ctx.fillText("x", 30, 22);
      const xf = Math.sqrt(a / 3);
      ctx.fillText("Fold", sx(folds[1]) + 8, sy(-xf) + 4);
      ctx.fillText("Fold", sx(folds[0]) - 34, sy(xf) + 4);
    }
    function advance(dt) {
      if (sweeping) {
        r += dir * 0.0025 * dt;
        if (r > 0.75) dir = -1;
        if (r < -0.75) dir = 1;
        slider.value = r.toFixed(3);
      }
      for (let k = 0; k < 8; k++) x = F.rk4(F.cusp(r, a), [x], 0.05 * dt)[0];
      trail.push([r, x]); if (trail.length > 900) trail.shift();
      out.textContent = r.toFixed(2);
    }
    root.addEventListener("resize", function () { g = fit(canvas); });
    if (reduced) {
      for (let k = 0; k < 1300; k++) { r += dir * 0.0025; if (r > 0.75) dir = -1; if (r < -0.75) dir = 1; for (let j = 0; j < 8; j++) x = F.rk4(F.cusp(r, a), [x], 0.05)[0]; trail.push([r, x]); if (trail.length > 900) trail.shift(); }
      draw(); return;
    }
    loop(canvas, function (ms) { advance(ms / 16); draw(); });
  }

  // ---------------------------------------------------------------------------
  // Lorenz attractor, slowly rotating projection of a long orbit.
  // ---------------------------------------------------------------------------
  function lorenzMini(canvas, opts) {
    opts = Object.assign({ colors: ["#EE6A24", "#FB9E07"], n: 2600 }, opts || {});
    let g = fit(canvas);
    const f = F.lorenz(10, 28, 8 / 3);
    let x = [1, 1, 20];
    for (let k = 0; k < 2000; k++) x = F.rk4(f, x, 0.01);
    const pts = [];
    for (let k = 0; k < opts.n; k++) { x = F.rk4(f, x, 0.01); pts.push(x); }
    let th = 0;
    const c0 = hexToRgb(opts.colors[0]), c1 = hexToRgb(opts.colors[1]);
    function draw() {
      const ctx = g.ctx, s = Math.min(g.w, g.h) / 58;
      ctx.clearRect(0, 0, g.w, g.h);
      const ct = Math.cos(th), st = Math.sin(th);
      ctx.lineWidth = 1;
      for (let i = 1; i < pts.length; i++) {
        const p = pts[i - 1], q = pts[i], u = i / pts.length;
        const col = c0.map(function (v, j) { return Math.round(v + u * (c1[j] - v)); });
        ctx.strokeStyle = "rgba(" + col.join(",") + "," + (0.2 + 0.7 * u).toFixed(2) + ")";
        ctx.beginPath();
        ctx.moveTo(g.w / 2 + s * (ct * p[0] - st * p[1]), g.h / 2 + s * (25 - p[2]) * 0.95);
        ctx.lineTo(g.w / 2 + s * (ct * q[0] - st * q[1]), g.h / 2 + s * (25 - q[2]) * 0.95);
        ctx.stroke();
      }
    }
    function advance() {
      for (let k = 0; k < 3; k++) { x = F.rk4(f, x, 0.01); pts.push(x); pts.shift(); }
      th += 0.003;
    }
    root.addEventListener("resize", function () { g = fit(canvas); });
    if (reduced) { draw(); return; }
    loop(canvas, function () { advance(); draw(); });
  }

  // ---------------------------------------------------------------------------
  // Knowledge graph from the exported graph of RElab_knowledge (vis-network).
  // ---------------------------------------------------------------------------
  function kbGraph(host, opts) {
    opts = Object.assign({ dark: false, highlight: "birkhoff-ergodic-theorem", height: "100%" }, opts || {});
    if (!root.vis || !root.KB_GRAPH) { host.textContent = "Graph data not loaded."; return; }
    const areaColor = { dynamics: "#2a1766", foundations: "#3093CF", ecology: "#118230", stochastics: "#A52C60", philosophy: "#EE6A24", computation: "#FB9E07" };
    if (opts.dark) areaColor.dynamics = "#9a84f0";
    const ink = opts.dark ? "#eceaf4" : "#212529";
    const nodes = root.KB_GRAPH.nodes.map(function (n) {
      const hi = n.id === opts.highlight;
      return {
        id: n.id, label: hi ? n.label : undefined, title: n.label + " (" + n.type + ")",
        shape: n.type === "question" ? "diamond" : n.type === "map" ? "square" : "dot",
        size: hi ? 14 : n.type === "map" ? 9 : 6,
        color: { background: hi ? "#FB9E07" : areaColor[n.area] || "#6c757d", border: hi ? "#EE6A24" : "rgba(0,0,0,0)" },
        font: { color: ink, face: "Jost", size: hi ? 22 : 14, strokeWidth: hi ? 4 : 0, strokeColor: opts.dark ? "#14111f" : "#ffffff" }
      };
    });
    const edges = root.KB_GRAPH.edges.map(function (e) {
      return { from: e.from, to: e.to, color: { color: opts.dark ? "rgba(236,234,244,0.16)" : "rgba(23,12,58,0.14)" }, width: 0.7 };
    });
    const net = new root.vis.Network(host, { nodes: nodes, edges: edges }, {
      physics: { solver: "forceAtlas2Based", forceAtlas2Based: { gravitationalConstant: -34, springLength: 60 }, stabilization: { iterations: 260 } },
      interaction: { hover: true, tooltipDelay: 80, zoomView: false },
      edges: { smooth: false }
    });
    net.once("stabilizationIterationsDone", function () { net.fit({ animation: false }); });
    return net;
  }

  root.Viz = { fit: fit, heteroclinicField: heteroclinicField, portraitField: portraitField, heteroclinicStatic: heteroclinicStatic, trace: trace, cuspInstrument: cuspInstrument, lorenzMini: lorenzMini, kbGraph: kbGraph, reduced: reduced };
})(this);
