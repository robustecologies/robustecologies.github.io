// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Engine checks against closed forms. Run: node tests/core.test.mjs
import { load, reporter, CORE } from "./harness.mjs";

const DF = load(CORE);
const R = reporter("core");
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const variance = (a) => { const m = mean(a); return a.reduce((s, v) => s + (v - m) * (v - m), 0) / (a.length - 1); };

// ------------------------------------------------------------------ RNG
{
  const r = new DF.RNG(7), N = 200000;
  const u = Array.from({ length: N }, () => r.uniform());
  R.check("uniform moments", Math.abs(mean(u) - 0.5) < 0.003 && Math.abs(variance(u) - 1 / 12) < 0.001 && u.every((v) => v > 0 && v < 1),
    `mean ${mean(u).toFixed(4)} (0.5), variance ${variance(u).toFixed(5)} (${(1 / 12).toFixed(5)}), N = ${N}`);
  const z = Array.from({ length: N }, () => r.normal());
  const k4 = mean(z.map((v) => v ** 4));
  R.check("normal moments", Math.abs(mean(z)) < 0.01 && Math.abs(variance(z) - 1) < 0.01 && Math.abs(k4 - 3) < 0.06,
    `mean ${mean(z).toFixed(4)}, variance ${variance(z).toFixed(4)}, fourth moment ${k4.toFixed(3)} (3)`);
  const s2 = Array.from({ length: N }, () => r.stable(2));
  R.check("stable alpha = 2 is normal with variance 2", Math.abs(variance(s2) - 2) < 0.03, `variance ${variance(s2).toFixed(4)}`);
  const s1 = Array.from({ length: N }, () => Math.abs(r.stable(1))).sort((a, b) => a - b);
  R.check("stable alpha = 1 is Cauchy (median |X| = 1)", Math.abs(s1[N >> 1] - 1) < 0.02, `median |X| = ${s1[N >> 1].toFixed(4)}`);
  const pz = Array.from({ length: 50000 }, () => r.poisson(3.2));
  R.check("Poisson mean and variance", Math.abs(mean(pz) - 3.2) < 0.03 && Math.abs(variance(pz) - 3.2) < 0.08, `mean ${mean(pz).toFixed(3)}, variance ${variance(pz).toFixed(3)} (3.2)`);
  const a = new DF.RNG(42), b = new DF.RNG(42), c = new DF.RNG(43);
  const sa = [a.uniform(), a.normal()], sb = [b.uniform(), b.normal()], sc = [c.uniform(), c.normal()];
  R.check("seeded reproducibility", sa[0] === sb[0] && sa[1] === sb[1] && sa[0] !== sc[0], `seed 42 twice gives ${sa[0].toFixed(6)} both times, seed 43 gives ${sc[0].toFixed(6)}`);
}

// ------------------------------------------------------------ formulas
{
  const ev = (src, x = 3) => {
    const sys = DF.compileSystem(`param a = 2\nq' = ${src}\ninit q = ${x}`);
    const dx = new Float64Array(1);
    sys.f(0.5, Float64Array.from([x]), Float64Array.from([2]), dx);
    return dx[0];
  };
  const cases = [
    ["-q^2", -9], ["2^3^2", 512], ["(1 + 2) * 3 - 4 / 2", 7], ["a*q - 1", 5], ["q > 2", 1], ["q <= 2", 0],
    ["ifelse(q > 1, 10, 20)", 10], ["mod(-1, 3)", 2], ["hill(q, 3, 2)", 0.5], ["min(q, a, 7)", 2], ["max(1, q)", 3],
    ["step(q - 3)", 1], ["clamp(q, 0, 1)", 1], ["sin(pi/2) + exp(0) + t", 2.5], ["2e-1 * 10", 2], ["(q > 1) && (a < 1)", 0], ["!(q > 5)", 1]
  ];
  const bad = cases.filter(([s, v]) => Math.abs(ev(s) - v) > 1e-12);
  R.check("expression semantics", bad.length === 0, `${cases.length - bad.length} of ${cases.length} expressions give the expected value${bad.length ? "; wrong: " + bad.map((c) => c[0] + " = " + ev(c[0])).join(", ") : ""}`);
  const errs = [
    ["x' = y", "Unknown name 'y'"], ["x' = foo(x)", "Unknown function 'foo'"], ["x' = pow(x)", "takes 2 arguments"],
    ["x' = x\ny[n+1] = y", "cannot be mixed"], ["x' = (x", "Expected ')'"], ["x' = x $ 2", "Unexpected character"],
    ["x' = lag(2, 1)", "must be a state variable"], ["param x = 1\nx' = x", "both a variable and a parameter"]
  ];
  const missed = errs.filter(([src, msg]) => { try { DF.compileSystem(src); return true; } catch (e) { return !String(e.message).includes(msg); } });
  R.check("formula errors are reported", missed.length === 0, `${errs.length - missed.length} of ${errs.length} malformed systems rejected with the expected message`);
  const lv = DF.compileSystem("# Lotka-Volterra\ndx/dt = a*x - b*x*y\ny' = d*x*y - c*y\nparam a = 1 [0, 3]\nparam b = 0.5\nparam c = 0.8\nparam d = 0.4\ninit x = 3\ninit y = 1\nrange x = [0, 6]");
  R.check("system metadata", lv.kind === "ode" && lv.vars.join() === "x,y" && lv.params.length === 4 && lv.params[0].max === 3 && lv.init[0] === 3 && lv.ranges.x[1] === 6,
    `kind ${lv.kind}, vars ${lv.vars}, params ${lv.params.map((q) => q.name)}, range of x ${lv.ranges.x}`);
  const kinds = [DF.compileSystem("x[n+1] = r*x*(1-x)\nparam r = 4").kind, DF.compileSystem("x' = -x\nnoise x = 0.1").kind, DF.compileSystem("x' = -lag(x, 1)").kind];
  R.check("kind detection", kinds.join() === "map,sde,dde", `map, sde, dde detected as ${kinds.join(", ")}`);
}

// ------------------------------------------------------------------ ODE
{
  const sys = DF.compileSystem("x' = -x\ninit x = 1");
  const run = (h) => { const s = new DF.Simulator(sys, { dt: h }); for (let k = 0; k < Math.round(1 / h); k++) s.step(); return Math.abs(s.X[0] - Math.exp(-1)); };
  const e1 = run(0.1), e2 = run(0.05);
  R.check("RK4 order", e1 / e2 > 15 && e1 / e2 < 17, `error ${e1.toExponential(3)} at h = 0.1, ${e2.toExponential(3)} at h = 0.05, ratio ${(e1 / e2).toFixed(2)} (16)`);
  const lv = DF.compileSystem("x' = a*x - b*x*y\ny' = d*x*y - c*y\nparam a = 1\nparam b = 0.5\nparam c = 0.8\nparam d = 0.4\ninit x = 3\ninit y = 1");
  const s = new DF.Simulator(lv, { dt: 0.01 });
  const H = (x, y) => 0.4 * x - 0.8 * Math.log(x) + 0.5 * y - 1 * Math.log(y);
  const H0 = H(3, 1);
  let drift = 0;
  for (let k = 0; k < 10000; k++) { s.step(); drift = Math.max(drift, Math.abs(H(s.X[0], s.X[1]) - H0)); }
  R.check("Lotka-Volterra first integral", drift / Math.abs(H0) < 1e-8, `max relative drift ${(drift / Math.abs(H0)).toExponential(2)} over t in [0, 100]`);
}

// ------------------------------------------------------------------ SDE
{
  // Ornstein-Uhlenbeck dx = -theta x dt + sigma dW under Euler-Maruyama has
  // stationary variance sigma^2 / (theta (2 - theta h)) exactly.
  const theta = 1, sigma = 0.5, h = 0.05;
  const sys = DF.compileSystem(`x' = -theta*x\nnoise x = sigma\nparam theta = ${theta}\nparam sigma = ${sigma}\ninit x = 0`);
  const s = new DF.Simulator(sys, { dt: h, n: 20000, seed: 3 });
  for (let k = 0; k < 400; k++) s.step();
  const v = variance(Array.from(s.X)), target = sigma * sigma / (theta * (2 - theta * h));
  R.check("Euler-Maruyama OU stationary variance", Math.abs(v / target - 1) < 0.03, `ensemble variance ${v.toFixed(5)}, exact discrete value ${target.toFixed(5)}, n = 20000`);
  const g = DF.compileSystem("x' = mu*x\nnoise x = s*x\nparam mu = 0.5\nparam s = 0.4\ninit x = 1");
  const sg = new DF.Simulator(g, { dt: 0.01, n: 40000, seed: 5 });
  for (let k = 0; k < 100; k++) sg.step();
  const m = mean(Array.from(sg.X)), exact = Math.pow(1 + 0.5 * 0.01, 100);
  R.check("geometric Brownian motion mean (Ito)", Math.abs(m / exact - 1) < 0.01, `E[x(1)] = ${m.toFixed(4)}, exact under the scheme ${exact.toFixed(4)}, e^0.5 = ${Math.exp(0.5).toFixed(4)}`);
}

// ------------------------------------------------------------------ DDE
{
  // x' = -x(t - 1), x(s) = 1 for s <= 0: x(2) = -1/2 and x(3) = -1/6 by the method of steps.
  const sys = DF.compileSystem("x' = -lag(x, 1)\ninit x = 1");
  const run = (h) => { const s = new DF.Simulator(sys, { dt: h }); const out = {}; for (let k = 1; k <= Math.round(3 / h); k++) { s.step(); if (k === Math.round(2 / h)) out.x2 = s.X[0]; } out.x3 = s.X[0]; return out; };
  const a = run(0.05), b = run(0.025);
  const e1 = Math.abs(a.x3 + 1 / 6), e2 = Math.abs(b.x3 + 1 / 6);
  R.check("DDE method-of-steps solution", Math.abs(b.x2 + 0.5) < 1e-12 && e2 < 1e-12, `x(2) = ${b.x2.toFixed(10)} (-0.5), x(3) = ${b.x3.toFixed(10)} (${(-1 / 6).toFixed(10)})`);
  R.check("DDE exact on piecewise cubics", e1 < 1e-12, `error at t = 3 is ${e1.toExponential(2)} at h = 0.05: the solution is a cubic on [2, 3], which Hermite interpolation and RK4 reproduce exactly`);
  // Mackey-Glass (tau = 17) against a four times finer step: fourth-order convergence.
  const mg = DF.compileSystem("x' = a*lag(x, tau)/(1 + lag(x, tau)^n) - b*x\nparam a = 0.2\nparam b = 0.1\nparam n = 10\nparam tau = 17\ninit x = 0.9");
  const at = (h) => { const s = new DF.Simulator(mg, { dt: h }); for (let k = 0; k < Math.round(50 / h); k++) s.step(); return s.X[0]; };
  const ref = at(0.0125), ea = Math.abs(at(0.1) - ref), eb = Math.abs(at(0.05) - ref);
  R.check("DDE convergence order (Mackey-Glass)", ea / eb > 12, `error at t = 50: ${ea.toExponential(2)} (h = 0.1), ${eb.toExponential(2)} (h = 0.05), ratio ${(ea / eb).toFixed(1)} (16 for order 4)`);
}

// ------------------------------------------------------------ Lyapunov
{
  const logi = DF.compileSystem("x[n+1] = r*x*(1-x)\nparam r = 4\ninit x = 0.3");
  const l1 = DF.largestLyapunov(logi, [4], [0.3], { steps: 200000, every: 1, transient: 100 });
  R.check("logistic r = 4 exponent", Math.abs(l1 - Math.log(2)) < 0.01, `lambda = ${l1.toFixed(4)}, ln 2 = ${Math.log(2).toFixed(4)}`);
  const hen = DF.compileSystem("x[n+1] = 1 - a*x^2 + y\ny[n+1] = b*x\nparam a = 1.4\nparam b = 0.3\ninit x = 0.1\ninit y = 0.1");
  const l2 = DF.largestLyapunov(hen, [1.4, 0.3], [0.1, 0.1], { steps: 200000, every: 1, transient: 1000 });
  R.check("Henon exponent", Math.abs(l2 - 0.419) < 0.01, `lambda = ${l2.toFixed(4)}, reference 0.419`);
  const lor = DF.compileSystem("x' = s*(y-x)\ny' = x*(r-z) - y\nz' = x*y - b*z\nparam s = 10\nparam r = 28\nparam b = 2.6666666666666665\ninit x = 1\ninit y = 1\ninit z = 1");
  const l3 = DF.largestLyapunov(lor, [10, 28, 8 / 3], [1, 1, 1], { dt: 0.005, steps: 200000, every: 20, transient: 10000 });
  R.check("Lorenz exponent", Math.abs(l3 - 0.906) < 0.04, `lambda = ${l3.toFixed(4)}, reference 0.906`);
}

// --------------------------------------------------------- perturbations
{
  const sys = DF.compileSystem("x' = 0\nparam a = 1\ninit x = 0");
  const s = new DF.Simulator(sys, { dt: 0.1, perturbations: [{ kind: "ramp", param: "a", rate: 0.5, t0: 1, span: 2 }] });
  const vals = [];
  for (let k = 0; k < 80; k++) { s.step(); if ([10, 30, 70].includes(s.steps)) vals.push(+s.p[0].toFixed(10)); }
  R.check("ramp", vals.join() === "1,2,3", `a(1) = ${vals[0]}, a(3) = ${vals[1]}, a(7) = ${vals[2]} (1, 2, 3 capped)`);
  const sp = new DF.Simulator(sys, { dt: 0.01, perturbations: [{ kind: "periodic", param: "a", amp: 0.5, period: 4 }] });
  for (let k = 0; k < 100; k++) sp.step();
  R.check("periodic modulation", Math.abs(sp.p[0] - 1.5) < 1e-9, `a(1) = ${sp.p[0].toFixed(10)} (1 + 0.5 sin(pi/2) = 1.5)`);
  const pulse = DF.compileSystem("x' = 0\ninit x = 1");
  const spu = new DF.Simulator(pulse, { dt: 0.1, perturbations: [{ kind: "pulse", var: "x", period: 1, frac: 0.5 }] });
  for (let k = 0; k < 35; k++) spu.step();
  R.check("pulses", Math.abs(spu.X[0] - 0.125) < 1e-12, `x(3.5) = ${spu.X[0]} after pulses at t = 1, 2, 3 removing half (0.125)`);
  const sj = new DF.Simulator(pulse, { dt: 0.01, n: 20000, seed: 9, init: [0], perturbations: [{ kind: "jumps", var: "x", rate: 2, size: 1 }] });
  for (let k = 0; k < 300; k++) sj.step();
  const mj = mean(Array.from(sj.X)), vj = variance(Array.from(sj.X));
  R.check("Poisson jumps", Math.abs(mj - 6) < 0.06 && Math.abs(vj - 6) < 0.25, `mean ${mj.toFixed(3)}, variance ${vj.toFixed(3)} after t = 3 at rate 2 (6 and 6)`);
  const sl = new DF.Simulator(pulse, { dt: 0.01, n: 20000, seed: 11, init: [0], perturbations: [{ kind: "levy", var: "x", sigma: 0.5, alpha: 2 }] });
  for (let k = 0; k < 100; k++) sl.step();
  R.check("alpha-stable noise, alpha = 2", Math.abs(variance(Array.from(sl.X)) / 0.5 - 1) < 0.03, `variance ${variance(Array.from(sl.X)).toFixed(4)} at t = 1 (2 sigma^2 t = 0.5)`);
  const sc = new DF.Simulator(DF.compileSystem("x' = -x\ninit x = 1"), { dt: 0.01, n: 50, seed: 2, spread: 0, perturbations: [{ kind: "additive", var: "x", sigma: 0.3, common: true }] });
  for (let k = 0; k < 200; k++) sc.step();
  const si = new DF.Simulator(DF.compileSystem("x' = -x\ninit x = 1"), { dt: 0.01, n: 50, seed: 2, spread: 0, perturbations: [{ kind: "additive", var: "x", sigma: 0.3 }] });
  for (let k = 0; k < 200; k++) si.step();
  const range = (a) => Math.max.apply(null, a) - Math.min.apply(null, a);
  R.check("common versus independent noise", range(Array.from(sc.X)) === 0 && variance(Array.from(si.X)) > 0.01,
    `range of 50 members after t = 2: ${range(Array.from(sc.X))} with common noise; variance ${variance(Array.from(si.X)).toFixed(4)} with independent noise`);
  const so = new DF.Simulator(sys, { dt: 0.01, seed: 4, perturbations: [{ kind: "ou", param: "a", sigma: 0.4, tau: 2 }] });
  const etas = [];
  for (let k = 0; k < 400000; k++) { so.step(); if (k > 2000 && k % 20 === 0) etas.push(so.p[0] - 1); }
  const target = 0.16 / (1 - 0.01 / 4);
  R.check("OU parameter noise variance", Math.abs(variance(etas) / target - 1) < 0.08, `variance ${variance(etas).toFixed(4)}, stationary value of the scheme ${target.toFixed(4)}`);
}

process.exit(R.done().fails ? 1 : 0);
