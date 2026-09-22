// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Every catalogue model runs, and the claims in its description hold.
import { load, reporter, CORE } from "./harness.mjs";

const DF = load(CORE.concat(["src/core/analysis.js", "src/models/catalogue.js"]));
const R = reporter("models");
const byId = (id) => DF.CATALOGUE.find((m) => m.id === id);
const sys = (id) => DF.compileSystem(byId(id).system);
const pv = (id, over) => { const s = sys(id); return s.params.map((q) => (over && over[q.name] !== undefined ? over[q.name] : q.value)); };
function sim(id, opts) {
  const m = byId(id), s = DF.compileSystem(m.system), sc = m.scene;
  return new DF.Simulator(s, Object.assign({ dt: sc.dt || (s.time === "discrete" ? 1 : 0.01), n: 1, seed: 7, perturbations: sc.perturbations || [], keepPositive: !!sc.keepPositive }, opts || {}));
}
function run(S, steps, each) { for (let k = 0; k < steps; k++) { S.step(); if (each) each(S, k); } return S; }
const eq = (id, over, box) => DF.findEquilibria(sys(id), pv(id, over), box);

// ---------------------------------------------------------------- all run
{
  const failed = [];
  for (const m of DF.CATALOGUE) {
    try {
      const S = sim(m.id, { n: 8, spread: 0.02 });
      run(S, m.id === "seasonal-sir" ? 20000 : 3000);
      // r-tipping escapes and Nicholson-Bailey diverges by design.
      const ok = Array.from(S.alive).some((a) => a) || m.id === "r-tipping" || m.id === "nicholson-bailey";
      if (!ok || !Array.from(S.X).some(isFinite)) failed.push(m.id);
    } catch (e) { failed.push(m.id + " (" + e.message + ")"); }
  }
  R.check("every model runs its default scene", failed.length === 0, `${DF.CATALOGUE.length - failed.length} of ${DF.CATALOGUE.length} models stay finite over 3000 steps${failed.length ? "; failed: " + failed.join(", ") : ""}`);
  const noSource = DF.CATALOGUE.filter((m) => !m.source || !m.about);
  R.check("every model is documented", noSource.length === 0, `${DF.CATALOGUE.length - noSource.length} entries carry a source and a description`);
  const ids = DF.CATALOGUE.map((m) => m.id);
  R.check("model identifiers are unique", new Set(ids).size === ids.length, `${ids.length} identifiers`);
}

// ------------------------------------------------------------- ecology
{
  const S = sim("lotka-volterra", { dt: 0.01 });
  const p = { r: 1, a: 0.05, e: 0.4, m: 0.4 };
  const H = (N, P) => p.e * p.a * N - p.m * Math.log(N) + p.a * P - p.r * Math.log(P);
  const H0 = H(20, 10); let d = 0;
  run(S, 20000, (s) => { d = Math.max(d, Math.abs(H(s.X[0], s.X[1]) - H0)); });
  R.check("Lotka-Volterra: H is conserved", d / Math.abs(H0) < 1e-8, `max relative change ${(d / Math.abs(H0)).toExponential(2)} over t in [0, 200]`);

  const box = [[0.01, 7], [0.01, 5]];
  const i1 = eq("rosenzweig-macarthur", { K: 3.7 }, box).find((e) => e.x[1] > 0.01), i2 = eq("rosenzweig-macarthur", { K: 3.8 }, box).find((e) => e.x[1] > 0.01);
  R.check("Rosenzweig-MacArthur: Hopf at K = 3.75", i1.stable && !i2.stable, `K = 3.7 ${i1.type}, K = 3.8 ${i2.type}`);

  const hp = DF.largestLyapunov(sys("hastings-powell"), pv("hastings-powell"), [0.8, 0.2, 8], { dt: 0.05, steps: 400000, every: 20, transient: 20000 });
  R.check("Hastings-Powell: chaotic", hp > 0.005, `largest Lyapunov exponent ${hp.toFixed(4)} per unit time`);

  // May-Leonard: residence times near the single-species states grow.
  const ml = sim("may-leonard", { dt: 0.02 }); const dw = []; let cur = -1, st = 0;
  run(ml, 150000, (s, k) => { const dom = [0, 1, 2].find((i) => s.X[i] > 0.9); const d0 = dom === undefined ? -1 : dom; if (d0 !== cur) { if (cur >= 0) dw.push((k - st) * 0.02); cur = d0; st = k; } });
  const grows = dw.length >= 5 && dw.every((v, i) => i === 0 || v > dw[i - 1]);
  R.check("May-Leonard: residence times grow", grows, `${dw.length} visits, dwell times ${dw.slice(0, 3).map((v) => v.toFixed(1)).join(", ")} ... ${dw.slice(-2).map((v) => v.toFixed(1)).join(", ")}`);

  const cl = eq("competition-lv", {}, [[0, 110], [0, 110]]);
  const cInt = cl.find((e) => e.x[0] > 1 && e.x[1] > 1);
  R.check("Lotka-Volterra competition: interior saddle", cInt && cInt.type === "saddle" && Math.abs(cInt.x[0] - 40) < 1e-6, `interior equilibrium (${cInt.x.map((v) => v.toFixed(3))}) is a ${cInt.type}; analytic value (40, 40)`);

  const g4 = eq("glv4", {}, [[0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1]]).find((e) => e.x.every((v) => v > 1e-6));
  const g4s = run(sim("glv4", { dt: 0.05 }), 20000);
  const g4d = Math.max(...g4.x.map((v, i) => Math.abs(v - g4s.X[i])));
  R.check("four-species gLV: stable interior equilibrium", g4.stable && g4d < 1e-6, `${g4.type} at (${g4.x.map((v) => v.toFixed(4))}); the orbit from the default start ends within ${g4d.toExponential(1)}`);

  const va = DF.largestLyapunov(sys("vano-lv4"), [], [0.3013, 0.4586, 0.1307, 0.3557], { dt: 0.05, steps: 400000, every: 20, transient: 20000 });
  R.check("Vano four-species LV: chaotic", Math.abs(va - 0.0203) < 0.004, `largest Lyapunov exponent ${va.toFixed(4)}; Vano et al. (2006) report 0.0203`);

  const gz = eq("grazing", {}, [[0, 10]]).map((e) => [e.x[0], e.stable]).sort((a, b) => a[0] - b[0]);
  const gzs = gz.map((e) => e[0].toFixed(3) + (e[1] ? " stable" : " unstable")).join(", ");
  R.check("May grazing: bistable", gz.length === 4 && gz[1][1] && !gz[2][1] && gz[3][1] && Math.abs(gz[3][0] - 7.317) < 1e-3 && Math.abs(gz[1][0] - 0.683) < 1e-3, `equilibria ${gzs}`);

  const rk = sys("ricker");
  const fp = (r) => DF.findEquilibria(rk, [r, 100], [[1, 400]]).find((e) => e.x[0] > 1);
  R.check("Ricker: fixed point K loses stability at r = 2", fp(1.95).stable && !fp(2.05).stable, `r = 1.95 ${fp(1.95).type}, r = 2.05 ${fp(2.05).type}`);

  const co = sim("coleman", { n: 10, initMode: "box", box: [[0.05, 2]], dt: 0.02 });
  run(co, 5000);
  const spread = Math.max(...co.X) - Math.min(...co.X);
  R.check("Coleman: orbits converge to one periodic orbit", spread < 1e-6, `spread of 10 orbits from [0.05, 2] after t = 100: ${spread.toExponential(1)}`);

  const tg = eq("toggle-switch", {}, [[0, 3.5], [0, 3.5]]);
  R.check("toggle switch: two stable states and a saddle", tg.length === 3 && tg.filter((e) => e.stable).length === 2 && tg.some((e) => e.type === "saddle"), tg.map((e) => e.type).join(", "));

  const rp = sim("rock-paper-scissors", { dt: 0.01 }); const P0 = 0.4 * 0.35 * 0.25; let dP = 0;
  run(rp, 20000, (s) => { dP = Math.max(dP, Math.abs(s.X[0] * s.X[1] * s.X[2] - P0)); });
  R.check("rock-paper-scissors: p1 p2 p3 conserved", dP / P0 < 1e-8, `max relative change ${(dP / P0).toExponential(2)}`);

  const hw = run(sim("huisman-weissing", { dt: 0.1 }), 20000);
  R.check("Huisman-Weissing: abundances stay non-negative and bounded", Array.from(hw.X).every((v) => v >= 0 && v < 100), `state after t = 2000: ${Array.from(hw.X).map((v) => v.toFixed(3)).join(" ")}`);
  // With mortality m equal to dilution D, T = R1 + sum_i c_1i N_i obeys T' = D (S - T), so T -> S.
  const tot = hw.X[5] + 0.20 * hw.X[0] + 0.10 * hw.X[1] + 0.10 * hw.X[2] + 0.10 * hw.X[3] + 0.10 * hw.X[4];
  R.check("Huisman-Weissing: resource-1 mass balance", Math.abs(tot - 10) < 1e-6, `R1 + sum c_1i N_i = ${tot.toFixed(6)} (T' = D (S - T) since m = D, so T -> S = 10)`);

  const Hst = 1.5 * Math.log(1.5) / (0.5 * 0.02), Pst = Math.log(1.5) / 0.02;
  const nb = sim("nicholson-bailey", { init: [Hst + 1, Pst] }); const Hs = [];
  run(nb, 60, (s) => Hs.push(s.X[0]));
  const early = Math.max(...Hs.slice(0, 15).map((v) => Math.abs(v - Hst))), late = Math.max(...Hs.slice(40).map((v) => Math.abs(v - Hst)));
  R.check("Nicholson-Bailey: growing oscillations", late > 3 * early, `start 1 unit from H* = ${Hst.toFixed(2)}: max |H - H*| ${early.toFixed(2)} in steps 1-15, ${late.toFixed(2)} in steps 41-60`);
}

// ----------------------------------------------------------- delays
{
  const amp = (r) => { const S = sim("hutchinson", { dt: 0.02, params: [r, 100, 1] }); let lo = Infinity, hi = -Infinity; run(S, 30000, (s, k) => { if (k > 28000) { lo = Math.min(lo, s.X[0]); hi = Math.max(hi, s.X[0]); } }); return hi - lo; };
  const a1 = amp(1.5), a2 = amp(1.65);
  R.check("Hutchinson: Hopf at r tau = pi/2", a1 < 1e-2 && a2 > 10, `oscillation amplitude after t = 560: ${a1.toExponential(1)} at r = 1.5, ${a2.toFixed(1)} at r = 1.65 (threshold pi/2 = 1.571)`);
  const mg = sim("mackey-glass", { dt: 0.1 }); let mn = Infinity, mx = -Infinity;
  run(mg, 20000, (s, k) => { if (k > 5000) { mn = Math.min(mn, s.X[0]); mx = Math.max(mx, s.X[0]); } });
  R.check("Mackey-Glass: bounded irregular oscillation", mn > 0.2 && mx < 1.5 && mx - mn > 0.5, `range [${mn.toFixed(3)}, ${mx.toFixed(3)}] after t = 500`);
}

// ------------------------------------------------------------- chaos
{
  const l = DF.largestLyapunov(sys("lorenz"), pv("lorenz"), [1, 1, 1], { dt: 0.005, steps: 200000, every: 20, transient: 10000 });
  R.check("Lorenz: lambda1 = 0.906", Math.abs(l - 0.906) < 0.04, `lambda1 = ${l.toFixed(4)}`);
  const r = DF.largestLyapunov(sys("rossler"), pv("rossler"), [1, 1, 1], { dt: 0.01, steps: 600000, every: 20, transient: 20000 });
  R.check("Rossler: lambda1 about 0.071", Math.abs(r - 0.0714) < 0.01, `lambda1 = ${r.toFixed(4)}; Sprott (2003) lists 0.0714`);
  const h = DF.largestLyapunov(sys("henon"), pv("henon"), [0.1, 0.1], { steps: 200000, every: 1, transient: 1000 });
  R.check("Henon: lambda1 = 0.419", Math.abs(h - 0.419) < 0.01, `lambda1 = ${h.toFixed(4)}`);
  const lz = DF.largestLyapunov(sys("lozi"), pv("lozi"), [0.1, 0], { steps: 200000, every: 1, transient: 1000 });
  R.check("Lozi: chaotic", Math.abs(lz - 0.47) < 0.02, `lambda1 = ${lz.toFixed(4)}; Sprott (2003) lists 0.4703`);
  const lg = DF.largestLyapunov(sys("logistic"), [4], [0.3], { steps: 200000, every: 1, transient: 100 });
  R.check("logistic: lambda = ln 2 at r = 4", Math.abs(lg - Math.log(2)) < 0.01, `lambda = ${lg.toFixed(4)}`);
  const sm = sys("standard-map"), J = DF.jacobian(sm, 0, Float64Array.from([1.3, 2.1]), Float64Array.from([1.2]));
  const det = J[0] * J[3] - J[1] * J[2];
  R.check("standard map: area preserving", Math.abs(det - 1) < 1e-8, `det J = ${det.toFixed(10)} at (1.3, 2.1)`);
  const ky = sim("kaplan-yorke", { n: 1 }); run(ky, 2000);
  R.check("Kaplan-Yorke: no collapse of the doubling coordinate", ky.X[0] > 1e-6, `x after 2000 steps = ${ky.X[0].toFixed(6)}`);
  const kyNoNoise = DF.compileSystem("x[n+1] = mod(2*x, 1)\ninit x = 0.3678"); const kz = new DF.Simulator(kyNoNoise, {}); run(kz, 60);
  R.check("doubling without re-injected bits collapses (control)", kz.X[0] === 0, `x after 60 steps = ${kz.X[0]}`);
  // Stark forced circle map: fibre exponent <ln|1 + a cos x|> along the orbit.
  const fibre = (a) => { const st = sim("stark-circle", { params: [3.883222077450933, 0.8676521529893011, a, 0.35] }); let sum = 0; run(st, 400000, (s, k) => { if (k >= 2000) sum += Math.log(Math.abs(1 + a * Math.cos(s.X[1]))); }); return sum / 398000; };
  const f9 = fibre(0.9), f25 = fibre(0.25);
  R.check("Stark circle map: fibre exponent", Math.abs(f9 + 0.235) < 0.01 && Math.abs(f25) < 1e-4, `a = 0.9: ${f9.toFixed(4)} (attracting graph); a = 0.25: ${f25.toExponential(1)} (zero, torus)`);
  const za = sim("zaslavsky", { n: 5, spread: 0 }); run(za, 500);
  const zr = Math.max(...Array.from(za.X).filter((v, i) => i % 2 === 1)) - Math.min(...Array.from(za.X).filter((v, i) => i % 2 === 1));
  R.check("Zaslavsky: one phase for the whole ensemble", zr === 0, `range of y over 5 identical members after 500 steps: ${zr}`);
}

// ----------------------------------------------------------- oscillators
{
  const pe = sim("pendulum", { dt: 0.01 }); const E = (x, v) => v * v / 2 - Math.cos(x), E0 = E(1, 0); let dE = 0;
  run(pe, 20000, (s) => { dE = Math.max(dE, Math.abs(E(s.X[0], s.X[1]) - E0)); });
  R.check("pendulum: energy conserved at c = 0", dE < 1e-9, `max |E - E0| = ${dE.toExponential(2)} over t in [0, 200]`);
  const hf = sim("hopf-normal-form", { dt: 0.01 }); run(hf, 20000);
  R.check("Hopf normal form: radius sqrt(mu)", Math.abs(Math.hypot(hf.X[0], hf.X[1]) - Math.sqrt(0.5)) < 1e-6, `radius ${Math.hypot(hf.X[0], hf.X[1]).toFixed(8)}, sqrt(0.5) = ${Math.sqrt(0.5).toFixed(8)}`);
  const kur = (K) => { const S = sim("kuramoto-pair", { dt: 0.01, params: [2 * Math.PI * 1.1, 2 * Math.PI * 0.9, K] }); const psi = []; run(S, 20000, (s, k) => { if (k % 1000 === 0) psi.push(s.X[0] - s.X[1]); }); return psi[psi.length - 1] - psi[psi.length - 5]; };
  const slip = kur(0.5), lock = kur(0.8);
  R.check("coupled phases: locking when 2K > |dw|", Math.abs(lock) < 1e-6 && Math.abs(slip) > 1, `drift of phi1 - phi2 over the last 40 time units: ${slip.toFixed(3)} at K = 0.5 (2K < 1.257), ${lock.toExponential(1)} at K = 0.8`);
  const osc = (id) => { const S = sim(id, { dt: 0.01 }); let lo = Infinity, hi = -Infinity; run(S, 30000, (s, k) => { if (k > 20000) { lo = Math.min(lo, s.X[0]); hi = Math.max(hi, s.X[0]); } }); return hi - lo; };
  const fh = osc("fitzhugh-nagumo"), sk = osc("selkov"), vdp = osc("van-der-pol");
  R.check("FitzHugh-Nagumo, Sel'kov and Van der Pol oscillate", fh > 2 && sk > 0.3 && vdp > 3.5, `late amplitude of the first variable: ${fh.toFixed(2)}, ${sk.toFixed(2)}, ${vdp.toFixed(2)}`);
  const dr = sim("driven-oscillator", { dt: 0.01 }); let hi = 0;
  run(dr, 60000, (s, k) => { if (k > 50000) hi = Math.max(hi, Math.abs(s.X[0])); });
  const Aexact = 0.5 / Math.hypot(1 - 1.69, 2 * 0.1 * 1.3);
  R.check("driven oscillator: response amplitude", Math.abs(hi - Aexact) < 1e-3, `amplitude ${hi.toFixed(5)}, F / sqrt((w0^2 - Om^2)^2 + (2 zeta Om)^2) = ${Aexact.toFixed(5)}`);
}

// ---------------------------------------------------------- stochastic
{
  const ou = sim("ornstein-uhlenbeck", { n: 20000, dt: 0.01, seed: 3 }); run(ou, 1500);
  const m = ou.X.reduce((a, v) => a + v, 0) / ou.n, v = ou.X.reduce((a, x) => a + (x - m) * (x - m), 0) / (ou.n - 1);
  const target = 0.25 / (1 * (2 - 0.01));
  R.check("Ornstein-Uhlenbeck: stationary variance", Math.abs(v / target - 1) < 0.03, `variance ${v.toFixed(4)}, sigma^2 / (2 theta) = 0.125, exact discrete value ${target.toFixed(4)}`);
  const dw = sim("double-well", { n: 2000, dt: 0.02, seed: 5 }); run(dw, 20000);
  const right = Array.from(dw.X).filter((x) => x > 0).length / dw.n;
  R.check("double well: escapes fill both wells", right > 0.4 && right < 0.6, `share in the right well after t = 400: ${right.toFixed(3)} (symmetric potential: 0.5)`);
}

// ------------------------------------------------ tipping and nonautonomous
{
  const rt = (r) => { const S = sim("r-tipping", { dt: 0.005, params: [r, 2, 12] }); let mid = null; run(S, 3000, (s, k) => { if (k === 2600) mid = s.X[0]; }); return { mid: mid, alive: S.alive[0], end: S.X[0] }; };
  const slow = rt(0.8), fast = rt(1.5);
  R.check("rate-induced tipping: threshold r = 1", Math.abs(slow.mid + Math.sqrt(0.2)) < 1e-3 && (!fast.alive || fast.end > 3), `r = 0.8: y = ${slow.mid.toFixed(5)} during the ramp (-sqrt(1 - r) = ${(-Math.sqrt(0.2)).toFixed(5)}); r = 1.5: ${fast.alive ? "y = " + fast.end.toFixed(2) : "escaped to infinity"}`);
  const fo = eq("fold-normal-form", { mu: -0.49 }, [[-1.6, 1.6]]).map((e) => e.x[0]).sort((a, b) => a - b);
  R.check("fold: equilibria +/-sqrt(-mu)", fo.length === 2 && Math.abs(fo[0] + 0.7) < 1e-9 && Math.abs(fo[1] - 0.7) < 1e-9, `mu = -0.49: ${fo.map((x) => x.toFixed(9)).join(", ")}`);
  const cu = (rr) => eq("cusp", { r: rr }, [[-1.6, 1.6]]).length;
  R.check("cusp: three equilibria between the folds", cu(0.38) === 3 && cu(0.39) === 1, `r = 0.38: ${cu(0.38)}, r = 0.39: ${cu(0.39)} (fold at 2/sqrt(27) = 0.3849)`);
  const pf = eq("pitchfork", { a: 0.64 }, [[-1.4, 1.4]]).map((e) => [e.x[0], e.stable]).sort((a, b) => a[0] - b[0]);
  R.check("pitchfork: +/-sqrt(a) stable, 0 unstable", pf.length === 3 && pf[0][1] && !pf[1][1] && pf[2][1] && Math.abs(pf[2][0] - 0.8) < 1e-9, pf.map((e) => e[0].toFixed(4) + (e[1] ? " s" : " u")).join(", "));
  const stom = sys("stommel"); let maxStable = 0, where = [];
  for (let e2 = 0.5; e2 <= 1.5 + 1e-9; e2 += 0.02) { const n = DF.findEquilibria(stom, [3, e2, 0.3], [[0, 4], [0, 4]]).filter((e) => e.stable).length; maxStable = Math.max(maxStable, n); if (n === 2) where.push(e2.toFixed(2)); }
  R.check("Stommel: bistable in eta2", maxStable === 2, `two stable states for eta2 in [${where[0]}, ${where[where.length - 1]}]`);
}

// ------------------------------------------------------------ epidemics
{
  const S = sim("sir", { dt: 0.05 }); run(S, 20000);
  const sInf = S.X[0] / 1000, R0 = 3;
  const resid = Math.log(sInf / 0.999) + R0 * (1 - sInf);
  R.check("SIR: final-size relation", Math.abs(resid) < 1e-3, `S(inf)/N = ${sInf.toFixed(5)}, ln(S_inf/S_0) + R0 (1 - S_inf) = ${resid.toExponential(2)} (0 up to the seed I0 = 1)`);
}

// ------------------------------------------------ scenes corrected by the audit
{
  // Delayed predator and prey: at the equilibrium (5, 3.75) the characteristic
  // equation is lambda^2 + 0.75 lambda exp(-lambda tau) + 0.375 = 0; roots i w
  // need cos(w tau) = 0, so w^2 - 0.75 w - 0.375 = 0 and tau* = (pi / 2) / w.
  const w = (0.75 + Math.sqrt(0.75 * 0.75 + 4 * 0.375)) / 2, tauStar = Math.PI / 2 / w;
  const amp = (tau) => {
    const S = sim("delayed-predator-prey", { dt: 0.01, params: pv("delayed-predator-prey", { tau: tau }), init: [5.2, 3.75] });
    run(S, 5000);
    const win = () => { let lo = Infinity, hi = -Infinity; run(S, 10000, (T) => { lo = Math.min(lo, T.X[0]); hi = Math.max(hi, T.X[0]); }); return [lo, hi]; };
    const a = win(), b = win();
    return { ratio: (b[1] - b[0]) / (a[1] - a[0]), range: b };
  };
  const below = amp(tauStar - 0.03), above = amp(tauStar + 0.03), dflt = amp(1.5);
  R.check("delayed predator and prey: Hopf point in tau", below.ratio < 1 && above.ratio > 1 && Math.abs(tauStar - 1.437) < 5e-4,
    `tau* = ${tauStar.toFixed(5)}; oscillation shrinks by ${below.ratio.toFixed(3)} per 100 time units at tau* - 0.03 and grows by ${above.ratio.toFixed(3)} at tau* + 0.03`);
  const S = sim("delayed-predator-prey", { dt: 0.01 });
  let lo = Infinity, hi = -Infinity, plo = Infinity, phi = -Infinity;
  run(S, 60000, (T, k) => { if (k > 40000) { lo = Math.min(lo, T.X[0]); hi = Math.max(hi, T.X[0]); plo = Math.min(plo, T.X[1]); phi = Math.max(phi, T.X[1]); } });
  R.check("delayed predator and prey: the default cycle fits the plot", lo > 1 && hi < 10.5 && plo > 2 && phi < 6 && hi - lo > 5,
    `at tau = 1.5, N in [${lo.toFixed(2)}, ${hi.toFixed(2)}] and P in [${plo.toFixed(2)}, ${phi.toFixed(2)}], inside N in [0, 14], P in [0, 12]`);
}
{
  // Charney-DeVore: one stable state at F = 2, a second pair of equilibria from
  // a fold near F = 2.5, and two stable states (zonal and blocked) at F = 4.
  const box = [[-10, 10], [-10, 10], [-10, 10]];
  const e2 = eq("charney-devore", { F: 2 }, box), e4 = eq("charney-devore", { F: 4 }, box);
  const st4 = e4.filter((e) => e.stable), sd4 = e4.filter((e) => /saddle/.test(e.type));
  R.check("Charney-DeVore: bistable at F = 4", e2.length === 1 && e2[0].stable && st4.length === 2 && sd4.length === 1,
    `F = 2: ${e2.length} equilibrium (${e2.map((e) => e.type).join(", ")}); F = 4: ${st4.length} stable (x = ${st4.map((e) => e.x[0].toFixed(3)).join(", ")}) and ${sd4.length} saddle`);
  const cd = DF.continueBranches(sys("charney-devore"), Float64Array.from(pv("charney-devore")), 0, 0.2, 5, [[-3, 6], [-4, 4], [-4, 4]]);
  const fd = cd.points.filter((q) => q.kind === "fold"), hp = cd.points.filter((q) => q.kind === "hopf");
  R.check("Charney-DeVore: fold and Hopf points in F", fd.length === 1 && Math.abs(fd[0].p - 2.51) < 0.02 && hp.length === 1 && Math.abs(hp[0].p - 3.41) < 0.02,
    `fold at F = ${fd.map((q) => q.p.toFixed(4)).join(", ")}, Hopf at F = ${hp.map((q) => q.p.toFixed(4)).join(", ")}; the blocked state is stable above the Hopf point`);
}

process.exit(R.done().fails ? 1 : 0);
