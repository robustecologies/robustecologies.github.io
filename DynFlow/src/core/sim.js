// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Simulation engine. A Simulator advances an ensemble of n copies of one
   compiled system, all sharing parameters, time and any common forcing, with
   a fixed step h:

     ode   classical fourth-order Runge-Kutta
     sde   Euler-Maruyama (Ito): x += f h + g dW, dW ~ N(0, h)
     dde   Runge-Kutta 4 with a cubic Hermite interpolant of the stored
           history; constant initial history x(s) = x0 for s <= t0
     map   x[n+1] = F(x[n]), one iteration per step

   Perturbations modify either parameters (periodic, quasiperiodic, ramp,
   step, Ornstein-Uhlenbeck) or states (additive, multiplicative, coloured and
   alpha-stable noise, Poisson jumps, periodic pulses). A perturbation marked
   common uses one realisation for the whole ensemble, as an environmental
   forcing does; otherwise every member receives its own. */
(function (DF) {
  "use strict";

  // ------------------------------------------------------------ RK4 step
  function makeRK4(dim) {
    const k1 = new Float64Array(dim), k2 = new Float64Array(dim), k3 = new Float64Array(dim), k4 = new Float64Array(dim), y = new Float64Array(dim);
    return function (f, t, x, p, h, H, R) {
      f(t, x, p, k1, H, R);
      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k1[i];
      f(t + 0.5 * h, y, p, k2, H, R);
      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k2[i];
      f(t + 0.5 * h, y, p, k3, H, R);
      for (let i = 0; i < dim; i++) y[i] = x[i] + h * k3[i];
      f(t + h, y, p, k4, H, R);
      for (let i = 0; i < dim; i++) x[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    };
  }

  // ----------------------------------------------------- DDE history
  // Ring buffer of (t, x, dx) at the step points; H(i, s) evaluates the
  // cubic Hermite interpolant, exact for cubic solutions.
  function History(dim, cap, x0, t0) {
    this.dim = dim; this.cap = cap; this.x0 = Float64Array.from(x0); this.t0 = t0;
    this.T = new Float64Array(cap); this.X = new Float64Array(cap * dim); this.D = new Float64Array(cap * dim);
    this.len = 0; this.head = 0;
  }
  History.prototype.push = function (t, x, dx) {
    const j = this.head;
    this.T[j] = t;
    for (let i = 0; i < this.dim; i++) { this.X[j * this.dim + i] = x[i]; this.D[j * this.dim + i] = dx[i]; }
    this.head = (j + 1) % this.cap; this.len = Math.min(this.len + 1, this.cap);
  };
  // Enlarge the buffer to newCap points, keeping every stored point in order.
  History.prototype.grow = function (newCap) {
    if (newCap <= this.cap) return;
    const dim = this.dim, T = new Float64Array(newCap), X = new Float64Array(newCap * dim), D = new Float64Array(newCap * dim);
    const oldest = (this.head - this.len + this.cap) % this.cap;
    for (let q = 0; q < this.len; q++) {
      const j = (oldest + q) % this.cap;
      T[q] = this.T[j];
      for (let i = 0; i < dim; i++) { X[q * dim + i] = this.X[j * dim + i]; D[q * dim + i] = this.D[j * dim + i]; }
    }
    this.T = T; this.X = X; this.D = D; this.cap = newCap; this.head = this.len % newCap;
  };
  History.prototype.at = function (i, s) {
    if (s <= this.t0 || this.len === 0) return this.x0[i];
    const cap = this.cap, dim = this.dim;
    const newest = (this.head - 1 + cap) % cap, oldest = (this.head - this.len + cap) % cap;
    if (s >= this.T[newest]) return this.X[newest * dim + i];
    // binary search over the logical order 0..len-1
    let lo = 0, hi = this.len - 1;
    const idx = function (q) { return (oldest + q) % cap; };
    if (s <= this.T[idx(0)]) return this.X[idx(0) * dim + i];
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (this.T[idx(mid)] <= s) lo = mid; else hi = mid; }
    const a = idx(lo), b = idx(hi);
    const ta = this.T[a], hstep = this.T[b] - ta, u = (s - ta) / hstep;
    const ya = this.X[a * dim + i], yb = this.X[b * dim + i], da = this.D[a * dim + i], db = this.D[b * dim + i];
    const u2 = u * u, u3 = u2 * u;
    return (2 * u3 - 3 * u2 + 1) * ya + (u3 - 2 * u2 + u) * hstep * da + (-2 * u3 + 3 * u2) * yb + (u3 - u2) * hstep * db;
  };

  // --------------------------------------------------------- perturbations
  /* Parameter modulators, applied to the base value b of parameter `param`:
       periodic       b + A sin(2 pi t / T + phase)
       quasiperiodic  b + A sin(2 pi t / T) + A2 sin(2 pi t / T2)
       ramp           b + rate (t - t0), clipped to [b, b + span] or [b + span, b]
       step           b before t0, b + A after
       ou             b + eta, d eta = -eta / tau dt + sigma sqrt(2 / tau) dW, common
     State perturbations, on variable `var`:
       additive       dx += sigma dW
       multiplicative dx += sigma x dW
       coloured       dx += eta dt, eta an OU process with time scale tau
       levy           dx += sigma h^(1/alpha) S_alpha, S_alpha symmetric stable
       jumps          at rate lambda, x += size, or x *= (1 - frac)
       pulse          every T time units, x += size, or x *= (1 - frac)  */
  const PARAM_KINDS = { periodic: 1, quasiperiodic: 1, ramp: 1, step: 1, ou: 1 };
  // Memory for the delay histories of one ensemble, in stored points times members times variables.
  const HISTORY_BUDGET = 4e6;

  function Simulator(system, opts) {
    opts = opts || {};
    this.sys = system;
    this.dim = system.vars.length;
    this.n = Math.max(1, opts.n || 1);
    this.h = opts.dt || 0.01;
    this.seed = opts.seed === undefined ? 1 : opts.seed;
    this.base = Float64Array.from(opts.params || system.params.map(function (q) { return q.value; }));
    this.p = Float64Array.from(this.base);
    this.perturbations = (opts.perturbations || []).map(function (q) { return Object.assign({}, q); });
    this.initMode = opts.initMode || "point";
    this.init = Float64Array.from(opts.init || system.init);
    this.spread = opts.spread === undefined ? 0.05 : opts.spread;
    this.box = opts.box || null;
    this.keepPositive = !!opts.keepPositive;
    // Deterministic skeleton: drift only, no state noise and no parameter noise
    // (the orbits of the phase plane and the bifurcation diagram of an SDE).
    this.deterministic = !!opts.deterministic;
    this.t0 = opts.t0 || 0;
    this.maxHistory = opts.maxHistory || Math.max(64, Math.min(200000, Math.floor(HISTORY_BUDGET / (this.n * this.dim))));
    this.historyClamped = false;
    this.reset();
  }

  Simulator.prototype.reset = function (seed) {
    if (seed !== undefined) this.seed = seed;
    const n = this.n, dim = this.dim;
    this.rng = new DF.RNG(this.seed);
    const rng = this.rng;
    this.R = { u: function () { return rng.uniform(); }, n: function () { return rng.normal(); }, U: new Float64Array(8), N: new Float64Array(8) };
    // Time is t0 + steps h, computed from an integer count so that it does not
    // drift through rounding: t equals k T exactly when T is a multiple of h.
    this.t = this.t0; this.steps = 0;
    this.X = new Float64Array(n * dim);
    this.alive = new Uint8Array(n).fill(1);
    for (let k = 0; k < n; k++) this.initMember(k);
    this.rk4 = makeRK4(dim);
    this.dx = new Float64Array(dim); this.gx = new Float64Array(dim); this.xk = new Float64Array(dim);
    const np = this.perturbations.length;
    this.eta = new Float64Array(np);         // common OU states
    this.etaK = new Float64Array(np * n);    // per-member OU states
    this.cDW = new Float64Array(np); this.cJump = new Float64Array(np); this.cStable = new Float64Array(np);
    this.nextPulse = this.perturbations.map(function (q) { return q.kind === "pulse" ? (q.t0 || q.period || 1) : Infinity; });
    this.hist = null;
    this.updateParams();
    if (this.sys.kind === "dde") {
      const cap = this.historyCap(this.lagBound());
      this.hist = [];
      for (let k = 0; k < n; k++) this.hist.push(this.newHistory(k, cap));
    }
  };

  /* Largest delay the scene can ask for: the delays at the current parameters
     and with each parameter at either end of its slider range, so that moving
     a slider does not outrun the stored history. */
  Simulator.prototype.lagBound = function () {
    const sys = this.sys, q = Float64Array.from(this.p);
    let m = sys.maxLag(q);
    sys.params.forEach(function (par, i) {
      const keep = q[i];
      q[i] = par.min; m = Math.max(m, sys.maxLag(q));
      q[i] = par.max; m = Math.max(m, sys.maxLag(q));
      q[i] = keep;
    });
    return m;
  };
  Simulator.prototype.historyCap = function (lag) {
    const need = isFinite(lag) ? Math.ceil(lag / this.h) + 8 : Infinity;
    if (need > this.maxHistory) this.historyClamped = true;
    return Math.min(this.maxHistory, need);
  };
  // Grow every member's history when the current delay needs more points.
  Simulator.prototype.ensureHistory = function () {
    if (!this.hist || !this.hist.length) return;
    const cap = this.historyCap(this.sys.maxLag(this.p));
    if (cap <= this.hist[0].cap) return;
    const c = Math.min(this.maxHistory, Math.max(cap, Math.ceil(this.hist[0].cap * 1.5)));
    this.hist.forEach(function (H) { H.grow(c); });
  };

  // History of member k starting now, holding the initial point with the
  // right derivative f(t0+, x0) so that the first step interpolates correctly.
  Simulator.prototype.newHistory = function (k, cap) {
    const dim = this.dim, x0 = this.X.slice(k * dim, k * dim + dim);
    const hk = new History(dim, cap, x0, this.t);
    const dx = new Float64Array(dim);
    this.sys.f(this.t, x0, this.p, dx, function (i) { return x0[i]; });
    hk.push(this.t, x0, dx);
    return hk;
  };

  Simulator.prototype.initMember = function (k) {
    const dim = this.dim, r = this.rng, x = this.X;
    for (let i = 0; i < dim; i++) {
      let v;
      if (this.initMode === "box" && this.box) v = r.range(this.box[i][0], this.box[i][1]);
      else if (this.initMode === "ball") v = this.init[i] + this.spread * r.normal();
      else v = this.init[i] + (k === 0 ? 0 : this.spread * (r.uniform() - 0.5));
      x[k * dim + i] = v;
    }
    this.alive[k] = 1;
  };

  // Place member k at state s (used by clicks and by respawning views).
  Simulator.prototype.setMember = function (k, s) {
    for (let i = 0; i < this.dim; i++) this.X[k * this.dim + i] = s[i];
    this.alive[k] = 1;
    if (this.hist) this.hist[k] = this.newHistory(k, this.hist[k].cap);
  };

  Simulator.prototype.paramIndex = function (name) {
    return this.sys.params.findIndex(function (q) { return q.name === name; });
  };

  Simulator.prototype.updateParams = function () {
    this.p.set(this.base);
    const t = this.t;
    for (let j = 0; j < this.perturbations.length; j++) {
      const q = this.perturbations[j];
      if (!(q.kind in PARAM_KINDS) || q.enabled === false) continue;
      const i = this.paramIndex(q.param);
      if (i < 0) continue;
      const b = this.base[i];
      switch (q.kind) {
        case "periodic": this.p[i] = b + q.amp * Math.sin(2 * Math.PI * t / q.period + (q.phase || 0)); break;
        case "quasiperiodic": this.p[i] = b + q.amp * Math.sin(2 * Math.PI * t / q.period) + (q.amp2 === undefined ? q.amp : q.amp2) * Math.sin(2 * Math.PI * t / (q.period2 || q.period * (1 + Math.sqrt(5)) / 2)); break;
        case "ramp": {
          const d = q.rate * Math.max(0, t - (q.t0 || 0));
          const span = q.span === undefined ? Infinity : q.span;
          this.p[i] = b + (span >= 0 ? Math.min(d, span) : Math.max(d, span));
          break;
        }
        case "step": this.p[i] = t >= (q.t0 || 0) ? b + q.amp : b; break;
        case "ou": this.p[i] = b + this.eta[j]; break;
      }
    }
    if (this.hist) this.ensureHistory();
  };

  // Advance the whole ensemble by one step of length h.
  Simulator.prototype.step = function () {
    const sys = this.sys, n = this.n, dim = this.dim, h = this.h, t = this.t, r = this.rng;
    const X = this.X, p = this.p, pert = this.perturbations, det = this.deterministic;
    const sqh = Math.sqrt(h), discrete = sys.time === "discrete";
    const tNext = this.t0 + (this.steps + 1) * (discrete ? 1 : h);
    // Common noise increments, drawn once per step for the whole ensemble.
    if (sys.usesRandom && !det) for (let i = 0; i < 8; i++) { this.R.U[i] = r.uniform(); this.R.N[i] = r.normal(); }
    const R = det ? undefined : this.R;
    const commonDW = this.cDW, commonJump = this.cJump, commonStable = this.cStable;
    if (!det) {
      for (let j = 0; j < pert.length; j++) { const q = pert[j]; commonDW[j] = q.common && q.enabled !== false ? r.normal() : 0; }
      for (let j = 0; j < pert.length; j++) { const q = pert[j]; commonJump[j] = q.kind === "jumps" && q.common && q.enabled !== false ? r.poisson(q.rate * h) : 0; }
      for (let j = 0; j < pert.length; j++) { const q = pert[j]; commonStable[j] = q.kind === "levy" && q.common && q.enabled !== false ? r.stable(q.alpha || 1.5) : 0; }
    }

    for (let k = 0; k < n; k++) {
      if (!this.alive[k]) continue;
      const x = X.subarray(k * dim, k * dim + dim);
      const Hk = this.hist ? this.hist[k] : null;
      const H = Hk ? function (i, s) { return Hk.at(i, s); } : undefined;
      if (discrete) {
        sys.f(t, x, p, this.dx, H, R);
        x.set(this.dx);
      } else if (sys.kind === "sde" && !det) {
        sys.f(t, x, p, this.dx, H, R);
        sys.g(t, x, p, this.gx, H, R);
        for (let i = 0; i < dim; i++) x[i] += this.dx[i] * h + (sys.noiseMask[i] ? this.gx[i] * sqh * r.normal() : 0);
      } else {
        this.rk4(sys.f, t, x, p, h, H, R);
      }
      // State perturbations; in discrete time h is 1 and dW has variance 1.
      for (let j = 0; j < pert.length && !det; j++) {
        const q = pert[j];
        if (q.enabled === false || q.kind in PARAM_KINDS || q.kind === "pulse") continue;
        const i = sys.vars.indexOf(q.var);
        if (i < 0) continue;
        const dW = (q.common ? commonDW[j] : r.normal()) * sqh;
        switch (q.kind) {
          case "additive": x[i] += q.sigma * dW; break;
          case "multiplicative": x[i] += q.sigma * x[i] * dW; break;
          case "coloured": {
            const tau = q.tau || 1, ix = q.common ? j : j * n + k, arr = q.common ? this.eta : this.etaK;
            if (!q.common || k === 0) arr[ix] += -arr[ix] / tau * h + q.sigma * Math.sqrt(2 / tau) * (q.common ? commonDW[j] : r.normal()) * sqh;
            x[i] += arr[ix] * h;
            break;
          }
          case "levy": x[i] += q.sigma * Math.pow(h, 1 / (q.alpha || 1.5)) * (q.common ? commonStable[j] : r.stable(q.alpha || 1.5)); break;
          case "jumps": {
            const m = q.common ? commonJump[j] : r.poisson(q.rate * h);
            for (let c = 0; c < m; c++) {
              if (q.frac !== undefined && q.frac !== null && q.frac !== "") x[i] *= (1 - q.frac);
              else x[i] += q.size + (q.sd ? q.sd * r.normal() : 0);
            }
            break;
          }
        }
      }
      if (this.keepPositive) for (let i = 0; i < dim; i++) if (x[i] < 0) x[i] = 0;
      let ok = true;
      for (let i = 0; i < dim; i++) if (!isFinite(x[i]) || Math.abs(x[i]) > 1e12) { ok = false; break; }
      if (!ok) this.alive[k] = 0;
      if (Hk && ok) { sys.f(tNext, x, p, this.dx, H); Hk.push(tNext, x, this.dx); }
    }
    this.steps++;
    this.t = tNext;
    // Pulses at fixed times act on every member at once.
    for (let j = 0; j < pert.length && !det; j++) {
      const q = pert[j];
      if (q.kind !== "pulse" || q.enabled === false) continue;
      const i = sys.vars.indexOf(q.var);
      while (i >= 0 && this.t >= this.nextPulse[j]) {
        for (let k = 0; k < n; k++) {
          const ix = k * dim + i;
          if (q.frac !== undefined && q.frac !== null && q.frac !== "") X[ix] *= (1 - q.frac); else X[ix] += q.size;
        }
        this.nextPulse[j] += q.period;
      }
    }
    // Common Ornstein-Uhlenbeck parameter noise.
    for (let j = 0; j < pert.length && !det; j++) {
      const q = pert[j];
      if (q.kind !== "ou" || q.enabled === false) continue;
      const tau = q.tau || 1;
      this.eta[j] += -this.eta[j] / tau * h + q.sigma * Math.sqrt(2 / tau) * r.normal() * sqh;
    }
    this.updateParams();
  };

  Simulator.prototype.member = function (k, out) {
    out = out || new Float64Array(this.dim);
    for (let i = 0; i < this.dim; i++) out[i] = this.X[k * this.dim + i];
    return out;
  };

  // --------------------------------------------------- Lyapunov estimate
  /* Largest Lyapunov exponent by two nearby orbits renormalised every
     `every` steps (Benettin et al. 1980), after a transient. Continuous time
     gives a rate per unit time; maps give a rate per iteration. */
  function largestLyapunov(system, params, init, opts) {
    opts = opts || {};
    const h = opts.dt || 0.01, steps = opts.steps || 100000, every = opts.every || 10, d0 = opts.d0 || 1e-8;
    const transient = opts.transient || 5000;
    const dim = system.vars.length, rk4 = makeRK4(dim), tmp = new Float64Array(dim);
    const p = Float64Array.from(params);
    const x = Float64Array.from(init), y = new Float64Array(dim);
    const adv = function (z, t) {
      if (system.time === "discrete") { system.f(t, z, p, tmp); z.set(tmp); } else rk4(system.f, t, z, p, h);
    };
    let t = 0;
    for (let s = 0; s < transient; s++) { adv(x, t); t += system.time === "discrete" ? 1 : h; }
    y.set(x); y[0] += d0;
    let sum = 0;
    for (let s = 1; s <= steps; s++) {
      adv(x, t); adv(y, t); t += system.time === "discrete" ? 1 : h;
      if (s % every === 0) {
        let d = 0;
        for (let i = 0; i < dim; i++) d += (y[i] - x[i]) * (y[i] - x[i]);
        d = Math.sqrt(d);
        sum += Math.log(d / d0);
        for (let i = 0; i < dim; i++) y[i] = x[i] + d0 * (y[i] - x[i]) / d;
      }
    }
    return sum / (steps * (system.time === "discrete" ? 1 : h));
  }

  DF.makeRK4 = makeRK4;
  DF.History = History;
  DF.Simulator = Simulator;
  DF.largestLyapunov = largestLyapunov;
  DF.PARAM_PERTURBATIONS = Object.keys(PARAM_KINDS);
  DF.STATE_PERTURBATIONS = ["additive", "multiplicative", "coloured", "levy", "jumps", "pulse"];
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
