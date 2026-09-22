// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Playback rate: every continuous flow, trajectory and phase scene of the
// catalogue moves at a similar speed on screen. The rate chosen by
// calibrateRate is checked against an independent replay of what the view
// shows, with more members, another seed and a longer horizon.
// Run: node tests/playback.test.mjs
import { load, reporter } from "./harness.mjs";

const DF = load(["src/core/rng.js", "src/core/expr.js", "src/core/sim.js", "src/core/analysis.js", "src/render/style.js", "src/render/frame.js", "src/render/mathtype.js", "src/render/views.js", "src/render/player.js", "src/models/catalogue.js"]);
const R = reporter("playback");

// What buildSim does before choosing a rate: ranges, displayed axes, step.
function setup(id) {
  const s = DF.normalizeScene(DF.sceneFor(id)), sys = DF.compileSystem(s.system);
  const params = sys.params.map((q) => (s.params[q.name] !== undefined ? +s.params[q.name] : q.value));
  const init = sys.vars.map((v, i) => (s.init[v] !== undefined ? +s.init[v] : sys.init[i]));
  const dt = s.dt || (sys.time === "discrete" ? 1 : 0.01);
  const need = sys.vars.filter((v) => !(s.view.ranges && s.view.ranges[v]) && !sys.ranges[v]);
  const probe = need.length ? DF.autoRanges(sys, params, init, { dt, steps: Math.min(20000, Math.max(3000, Math.round(60 / dt))), perturbations: s.perturbations }) : {};
  const ranges = sys.vars.map((v) => (s.view.ranges && s.view.ranges[v]) || sys.ranges[v] || probe[v]);
  let axes = (s.view.axes || []).map((v) => sys.vars.indexOf(v)).filter((i) => i >= 0);
  if (axes.length < 2) axes = sys.vars.length >= 3 && (s.view.type === "trajectory" || s.view.type === "flow") && s.view.dim3 !== false ? [0, 1, 2] : [0, 1];
  if (s.view.projection === "simplex") axes = [0, 1, 2];
  axes = Array.from(new Set(axes));
  const o = { view: s.view, dt, params, init, ranges, axes, spread: s.spread, initMode: s.initMode, perturbations: s.perturbations, keepPositive: !!s.keepPositive };
  const threeD = axes.length >= 3 && s.view.projection !== "simplex";
  return { s, sys, o, target: { trajectory: 1, flow: 0.6, phase: 0.5 }[s.view.type] * (threeD ? 1.5 : 1) };
}
// Independent replay: median path length per second over half-second windows
// (phase: over the first second of 16 fresh orbits), 8 trajectory members over
// 20 s, 120 particles over 5 s.
function shown(M, rate) {
  const { s, sys, o } = M, v = s.view, h = o.dt, dim = sys.vars.length, W = o.ranges.map((r) => r[1] - r[0]), perFrame = rate / (60 * h);
  let acc = 0;
  const stepsNow = () => { acc += perFrame; const n = Math.floor(acc); acc -= n; return n; };
  const disp = (X, k, P) => Math.sqrt(o.axes.reduce((a, i) => a + ((X[k * dim + i] - P[k * dim + i]) / W[i]) ** 2, 0));
  const out = [];
  if (v.type === "phase") {
    const rng = new DF.RNG(2024);
    for (let t = 0; t < 16; t++) {
      const x0 = Float64Array.from(o.init);
      if (t) o.axes.forEach((i) => { x0[i] = rng.range(o.ranges[i][0], o.ranges[i][1]); });
      const sim = new DF.Simulator(sys, { n: 1, dt: h, params: o.params, init: x0, deterministic: true, perturbations: o.perturbations.filter((q) => ["periodic", "quasiperiodic", "ramp", "step"].includes(q.kind)) });
      const prev = Float64Array.from(sim.X);
      let path = 0;
      acc = 0;
      for (let f = 0; f < 60; f++) {
        const n = stepsNow(); for (let k = 0; k < n; k++) sim.step();
        if (!sim.alive[0] || sim.X.some((x, i) => x < o.ranges[i][0] - 2 * W[i] || x > o.ranges[i][1] + 2 * W[i])) break;
        path += disp(sim.X, 0, prev); prev.set(sim.X);
      }
      out.push(path);
    }
  } else {
    const flowLives = v.type === "flow" && Array.isArray(v.life);
    const n = flowLives ? 120 : 8, meanW = W.reduce((a, b) => a + b, 0) / dim;
    const sim = new DF.Simulator(sys, { n, dt: h, params: o.params, init: o.init, seed: 2024, box: o.ranges, perturbations: o.perturbations, keepPositive: o.keepPositive, spread: v.type === "flow" ? o.spread : 0.01 * meanW, initMode: v.type === "flow" ? o.initMode : "point" });
    const r = sim.rng, life = new Float64Array(n).fill(Infinity), age = new Float64Array(n);
    const spawn = (k) => {
      const mode = v.spawn === "mixed" ? (r.uniform() < 0.25 ? "init" : "box") : v.spawn;
      const x = new Float64Array(dim);
      for (let i = 0; i < dim; i++) x[i] = mode === "init" ? o.init[i] + (o.spread || 0.02) * W[i] * r.normal() : r.range(o.ranges[i][0], o.ranges[i][1]);
      sim.setMember(k, x);
    };
    if (flowLives) for (let k = 0; k < n; k++) { spawn(k); life[k] = v.life[0] + r.uniform() * (v.life[1] - v.life[0]); age[k] = r.uniform() * life[k]; }
    const warm = Math.min(60000, Math.round((v.warmup === undefined ? (v.type === "flow" ? 60 : 0) : v.warmup) * perFrame));
    for (let k = 0; k < warm; k++) sim.step();
    const prev = Float64Array.from(sim.X), wl = new Float64Array(n), wc = new Int32Array(n), frames = flowLives ? 300 : 1200;
    for (let f = 0; f < frames; f++) {
      const m = stepsNow(); for (let k = 0; k < m; k++) sim.step();
      for (let k = 0; k < n; k++) {
        if (flowLives && ++age[k] > life[k]) { spawn(k); age[k] = 0; wl[k] = 0; wc[k] = 0; for (let i = 0; i < dim; i++) prev[k * dim + i] = sim.X[k * dim + i]; continue; }
        if (!sim.alive[k]) continue;
        wl[k] += disp(sim.X, k, prev);
        if (++wc[k] === 30) { out.push(wl[k] * 2); wl[k] = 0; wc[k] = 0; }
      }
      prev.set(sim.X);
    }
  }
  const a = out.filter(isFinite).sort((x, y) => x - y);
  return a[a.length >> 1];
}

{
  const rows = [];
  for (const m of DF.CATALOGUE) {
    const M = setup(m.id);
    if (!["flow", "trajectory", "phase"].includes(M.s.view.type) || M.sys.time !== "continuous") continue;
    const rate = DF.calibrateRate(M.sys, M.o).rate;
    rows.push({ id: m.id, type: M.s.view.type, rate, rel: shown(M, rate) / M.target });
  }
  for (const t of ["trajectory", "flow", "phase"]) {
    const r = rows.filter((q) => q.type === t), rel = r.map((q) => q.rel), lo = Math.min(...rel), hi = Math.max(...rel);
    const worst = r.slice().sort((a, b) => Math.abs(Math.log(b.rel)) - Math.abs(Math.log(a.rel)))[0];
    // Phase portraits that settle on a node are held within 4 times the first estimate, so they can play slower.
    const floor = t === "phase" ? 0.3 : 0.5;
    R.check("calibrated speed of " + t + " scenes", lo > floor && hi < 2, `${r.length} scenes shown at ${lo.toFixed(2)} to ${hi.toFixed(2)} times the target (spread ${(hi / lo).toFixed(2)}); farthest from it: ${worst.id} at ${worst.rel.toFixed(2)}`);
  }
}
{
  // Scenes saved with a fixed number of steps per frame keep their old pace at 60 frames per second.
  const s = DF.normalizeScene({ system: "x' = r + a*x - x^3\nparam r = 0 [-0.8, 0.8]\nparam a = 1", dt: 0.01, stepsPerFrame: 10, view: { type: "sweep", param: "r", speed: 0.0005 } });
  const m = DF.normalizeScene({ system: "x[n+1] = r*x*(1 - x)\nparam r = 3.7", stepsPerFrame: 2, view: { type: "flow" } });
  R.check("scenes with steps per frame", s.rate === 6 && Math.abs(s.view.speed - 0.005) < 1e-15 && s.stepsPerFrame === undefined && s.version === 2 && m.rate === 120,
    `flow: 10 steps of 0.01 per frame give rate ${s.rate} per second and sweep speed ${s.view.speed} per unit of time; a map at 2 steps per frame gives ${m.rate} iterations per second`);
  const c = DF.calibrateRate(DF.compileSystem("x' = r + a*x - x^3\nparam r = 0\nparam a = 1"), { view: { type: "sweep" }, dt: 0.01, sweepSpeed: 0.005 });
  const w = DF.calibrateRate(DF.compileSystem("x' = -x"), { view: { type: "timeseries", window: 40 }, dt: 0.01 });
  R.check("sweep and time series rates", c.rate === 8 && w.rate === 4, `a sweep at 0.005 of the interval per unit time crosses it in 1 / (0.005 x ${c.rate}) = ${1 / (0.005 * c.rate)} s; a window of 40 scrolls in ${40 / w.rate} s`);
}

process.exit(R.done().fails ? 1 : 0);
