// RElabFlow 0.1.0. Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab.
// SPDX-License-Identifier: GPL-3.0-or-later
// Licence text: https://www.gnu.org/licenses/gpl-3.0.html
// Built from 12 source files by tools/build.mjs; edit the sources, not this file.
// ---- src/core/rng.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Seeded random numbers. Every stochastic scene is reproducible from its
   seed: the same seed gives the same noise, the same jumps and the same
   initial ensemble in the studio, in an exported page and in the tests. */
(function (DF) {
  "use strict";

  // splitmix32 expands one 32-bit seed into the four words of xoshiro128**.
  function splitmix32(a) {
    return function () {
      a |= 0; a = (a + 0x9e3779b9) | 0;
      let z = a;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b);
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35);
      return (z ^ (z >>> 16)) >>> 0;
    };
  }

  // xoshiro128** (Blackman and Vigna), period 2^128 - 1.
  function RNG(seed) {
    const sm = splitmix32(seed === undefined ? 1 : seed >>> 0);
    this.s = new Uint32Array([sm(), sm(), sm(), sm()]);
    this._spare = null;
  }
  RNG.prototype.nextU32 = function () {
    const s = this.s;
    const r = Math.imul(rotl(Math.imul(s[1], 5), 7), 9) >>> 0;
    const t = s[1] << 9;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 11);
    return r;
  };
  function rotl(x, k) { return (x << k) | (x >>> (32 - k)); }

  // Uniform on (0, 1), never exactly 0 or 1, with 53 random bits.
  RNG.prototype.uniform = function () {
    const hi = this.nextU32() >>> 5, lo = this.nextU32() >>> 6;
    return (hi * 67108864 + lo + 0.5) / 9007199254740992;
  };
  RNG.prototype.range = function (a, b) { return a + (b - a) * this.uniform(); };

  // Standard normal by the Box-Muller transform, one value cached.
  RNG.prototype.normal = function () {
    if (this._spare !== null) { const v = this._spare; this._spare = null; return v; }
    const u = this.uniform(), v = this.uniform();
    const r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * v;
    this._spare = r * Math.sin(th);
    return r * Math.cos(th);
  };
  RNG.prototype.exponential = function (rate) { return -Math.log(this.uniform()) / rate; };

  // Poisson by inversion for small means and a normal approximation above 60.
  RNG.prototype.poisson = function (mu) {
    if (mu <= 0) return 0;
    if (mu > 60) return Math.max(0, Math.round(mu + Math.sqrt(mu) * this.normal()));
    const L = Math.exp(-mu);
    let k = 0, p = 1;
    do { k++; p *= this.uniform(); } while (p > L);
    return k - 1;
  };

  // Symmetric alpha-stable variate, scale 1, by Chambers, Mallows and Stuck
  // (1976); alpha = 2 gives a normal with variance 2, alpha = 1 a Cauchy.
  RNG.prototype.stable = function (alpha) {
    const V = Math.PI * (this.uniform() - 0.5), W = this.exponential(1);
    if (Math.abs(alpha - 1) < 1e-12) return Math.tan(V);
    return Math.sin(alpha * V) / Math.pow(Math.cos(V), 1 / alpha) *
      Math.pow(Math.cos(V - alpha * V) / W, (1 - alpha) / alpha);
  };

  DF.RNG = RNG;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/core/expr.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Formula language. A system is written as plain text, one statement per
   line, and compiled to JavaScript functions that write into preallocated
   arrays. The parser accepts only numbers, declared names, whitelisted
   functions and operators, so a compiled formula cannot reach anything else.

     # Lotka-Volterra predator and prey
     x' = a*x - b*x*y            differential equation (also dx/dt = ...)
     y' = d*x*y - c*y
     x[n+1] = r*x*(1 - x)        difference equation (map)
     noise x = sigma*x           diffusion coefficient of dW_x (Ito)
     aux h = x/(1 + x)           helper evaluated before the equations
     param a = 1 [0, 3]          parameter, value and slider range
     init x = 0.5                initial condition
     range x = [0, 4]            axis range used by the views
     lag(x, tau)                 delayed state x(t - tau), inside an expression

   Functions: sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10
   log2 sqrt cbrt abs sign floor ceil round min max pow mod step heaviside
   ifelse clamp hill; random draws urand() and nrand() (fresh for every member)
   and ucommon(k), ncommon(k) for k = 0..7 (shared by the ensemble at each
   step, as the common forcing of a random dynamical system); constants pi and e; operators + - * / ^, unary minus,
   comparisons < <= > >= == != (value 1 or 0), && and ||. */
(function (DF) {
  "use strict";

  const FUNCS = {
    sin: 1, cos: 1, tan: 1, asin: 1, acos: 1, atan: 1, sinh: 1, cosh: 1, tanh: 1,
    exp: 1, log: 1, log10: 1, log2: 1, sqrt: 1, cbrt: 1, abs: 1, sign: 1,
    floor: 1, ceil: 1, round: 1, atan2: 2, pow: 2, mod: 2, min: -1, max: -1,
    step: 1, heaviside: 1, ifelse: 3, clamp: 3, hill: 3, lag: 2,
    urand: 0, nrand: 0, ucommon: 1, ncommon: 1
  };
  const CONSTS = { pi: "Math.PI", e: "Math.E" };

  function FormulaError(msg, line) { this.message = msg + (line ? " (line " + line + ")" : ""); this.line = line; }
  FormulaError.prototype = Object.create(Error.prototype);
  FormulaError.prototype.name = "FormulaError";

  // ---------------------------------------------------------------- lexer
  function tokenize(src, line) {
    const toks = [];
    let i = 0;
    while (i < src.length) {
      const c = src[i];
      if (c === " " || c === "\t") { i++; continue; }
      const num = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(src.slice(i));
      if (num) { toks.push({ t: "num", v: num[0] }); i += num[0].length; continue; }
      const id = /^[A-Za-z_][A-Za-z0-9_]*/.exec(src.slice(i));
      if (id) { toks.push({ t: "id", v: id[0] }); i += id[0].length; continue; }
      const op = /^(<=|>=|==|!=|&&|\|\||[-+*/^(),<>!])/.exec(src.slice(i));
      if (op) { toks.push({ t: "op", v: op[0] }); i += op[0].length; continue; }
      throw new FormulaError("Unexpected character '" + c + "'", line);
    }
    toks.push({ t: "end" });
    return toks;
  }

  // ---------------------------------------------------------------- parser
  // Pratt parser; binding powers follow the usual precedence, ^ is right
  // associative and binds tighter than unary minus, so -x^2 = -(x^2).
  const BP = { "||": 1, "&&": 2, "==": 3, "!=": 3, "<": 4, "<=": 4, ">": 4, ">=": 4, "+": 5, "-": 5, "*": 6, "/": 6, "^": 8 };

  function parse(src, line) {
    const toks = tokenize(src, line);
    let k = 0;
    const peek = function () { return toks[k]; };
    const next = function () { return toks[k++]; };
    function expect(v) {
      const tk = next();
      if (tk.t !== "op" || tk.v !== v) throw new FormulaError("Expected '" + v + "'", line);
    }
    function nud(tk) {
      if (tk.t === "num") return { k: "num", v: parseFloat(tk.v) };
      if (tk.t === "id") {
        if (peek().t === "op" && peek().v === "(") {
          next();
          const args = [];
          if (!(peek().t === "op" && peek().v === ")")) {
            for (;;) { args.push(expr(0)); if (peek().t === "op" && peek().v === ",") { next(); continue; } break; }
          }
          expect(")");
          if (!(tk.v in FUNCS)) throw new FormulaError("Unknown function '" + tk.v + "'", line);
          const ar = FUNCS[tk.v];
          if (ar >= 0 && args.length !== ar) throw new FormulaError("Function '" + tk.v + "' takes " + ar + " argument" + (ar > 1 ? "s" : ""), line);
          if (ar < 0 && args.length < 1) throw new FormulaError("Function '" + tk.v + "' needs arguments", line);
          return { k: "call", f: tk.v, args: args };
        }
        return { k: "name", v: tk.v };
      }
      if (tk.t === "op" && tk.v === "(") { const e = expr(0); expect(")"); return e; }
      if (tk.t === "op" && tk.v === "-") return { k: "neg", a: expr(7) };
      if (tk.t === "op" && tk.v === "+") return expr(7);
      if (tk.t === "op" && tk.v === "!") return { k: "not", a: expr(7) };
      throw new FormulaError(tk.t === "end" ? "Incomplete expression" : "Unexpected '" + tk.v + "'", line);
    }
    function expr(rbp) {
      let left = nud(next());
      for (;;) {
        const tk = peek();
        if (tk.t !== "op" || !(tk.v in BP) || BP[tk.v] <= rbp) break;
        next();
        const bp = BP[tk.v];
        const right = expr(tk.v === "^" ? bp - 1 : bp);
        left = { k: "bin", op: tk.v, a: left, b: right };
      }
      return left;
    }
    const ast = expr(0);
    if (peek().t !== "end") throw new FormulaError("Unexpected '" + peek().v + "'", line);
    return ast;
  }

  // ------------------------------------------------------------ code emit
  // scope maps a name to its JavaScript expression; lags collects lag() use.
  function emit(node, scope, line, info) {
    switch (node.k) {
      case "num": return "(" + String(node.v) + ")";
      case "name":
        if (node.v in scope) return scope[node.v];
        if (node.v in CONSTS) return CONSTS[node.v];
        if (node.v === "t") { info.usesTime = true; return "t"; }
        throw new FormulaError("Unknown name '" + node.v + "'", line);
      case "neg": return "(-" + emit(node.a, scope, line, info) + ")";
      case "not": return "(" + emit(node.a, scope, line, info) + " ? 0 : 1)";
      case "bin": {
        const a = emit(node.a, scope, line, info), b = emit(node.b, scope, line, info);
        if (node.op === "^") return "Math.pow(" + a + ", " + b + ")";
        if (node.op === "&&" || node.op === "||") return "((" + a + " " + node.op + " " + b + ") ? 1 : 0)";
        if (["<", "<=", ">", ">=", "==", "!="].indexOf(node.op) >= 0) return "((" + a + " " + (node.op === "==" ? "===" : node.op === "!=" ? "!==" : node.op) + " " + b + ") ? 1 : 0)";
        return "(" + a + " " + node.op + " " + b + ")";
      }
      case "call": {
        const f = node.f;
        if (f === "lag") {
          const v = node.args[0];
          if (v.k !== "name" || !(v.v in info.varIndex)) throw new FormulaError("The first argument of lag() must be a state variable", line);
          info.usesLag = true;
          const tau = emit(node.args[1], scope, line, info);
          info.lagExprs.push(tau);
          return "H(" + info.varIndex[v.v] + ", t - (" + tau + "))";
        }
        const args = node.args.map(function (a) { return emit(a, scope, line, info); });
        switch (f) {
          case "mod": return "__mod(" + args.join(", ") + ")";
          case "step": case "heaviside": return "((" + args[0] + ") >= 0 ? 1 : 0)";
          case "ifelse": return "((" + args[0] + ") ? (" + args[1] + ") : (" + args[2] + "))";
          case "clamp": return "Math.min(Math.max(" + args[0] + ", " + args[1] + "), " + args[2] + ")";
          case "hill": return "__hill(" + args.join(", ") + ")";
          case "urand": info.usesRandom = true; return "R.u()";
          case "nrand": info.usesRandom = true; return "R.n()";
          case "ucommon": info.usesRandom = true; return "R.U[Math.min(7, Math.max(0, (" + args[0] + ") | 0))]";
          case "ncommon": info.usesRandom = true; return "R.N[Math.min(7, Math.max(0, (" + args[0] + ") | 0))]";
          default: return "Math." + f + "(" + args.join(", ") + ")";
        }
      }
    }
    throw new FormulaError("Cannot compile expression", line);
  }

  const PRELUDE =
    "const __mod = function (a, b) { return a - b * Math.floor(a / b); };\n" +
    "const __hill = function (x, K, n) { const u = Math.pow(Math.max(x, 0), n); return u / (Math.pow(K, n) + u); };\n" +
    // Without a random source (analysis, probes) random draws take their mean values.
    "const __R0 = { u: function () { return 0.5; }, n: function () { return 0; }, U: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], N: [0, 0, 0, 0, 0, 0, 0, 0] };\n";

  // ------------------------------------------------------- system parsing
  const RE = {
    ode: /^(?:d\s*([A-Za-z_]\w*)\s*\/\s*dt|([A-Za-z_]\w*)\s*')\s*=\s*(.+)$/,
    map: /^([A-Za-z_]\w*)\s*(?:\[\s*n\s*\+\s*1\s*\]|_\{?\s*n\s*\+\s*1\s*\}?|_next)\s*=\s*(.+)$/,
    noise: /^noise\s+([A-Za-z_]\w*)\s*=\s*(.+)$/,
    aux: /^aux\s+([A-Za-z_]\w*)\s*=\s*(.+)$/,
    param: /^param\s+([A-Za-z_]\w*)\s*=\s*([-+0-9.eE]+)\s*(?:\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\])?\s*$/,
    init: /^init\s+([A-Za-z_]\w*)\s*=\s*([-+0-9.eE]+)\s*$/,
    range: /^range\s+([A-Za-z_]\w*)\s*=\s*\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\]\s*$/
  };
  // Declared names may shadow the constants pi and e (e is a common
  // conversion efficiency); function names and t are reserved.
  const RESERVED = Object.assign({ t: 1, noise: 1, aux: 1, param: 1, init: 1, range: 1 }, FUNCS);

  /* Parse and compile a system. Returns
       { kind: 'ode' | 'map' | 'sde' | 'dde', vars, params, init, ranges,
         f(t, x, p, dx, H), g(t, x, p, gx) or null, noiseMask, maxLag, source } */
  function compileSystem(text) {
    const lines = String(text).split(/\r?\n/);
    const eqs = [], noises = [], auxes = [], params = [], init = {}, ranges = {};
    let kind = null;
    lines.forEach(function (raw, i) {
      const ln = i + 1;
      const s = raw.replace(/#.*$/, "").trim();
      if (!s) return;
      let m;
      if ((m = RE.param.exec(s))) {
        const v = +m[2];
        const lo = m[3] !== undefined ? +m[3] : (v === 0 ? -1 : Math.min(0, 2 * v));
        const hi = m[4] !== undefined ? +m[4] : (v === 0 ? 1 : Math.max(0, 2 * v));
        if (params.some(function (q) { return q.name === m[1]; })) throw new FormulaError("Parameter '" + m[1] + "' is declared twice", ln);
        params.push({ name: m[1], value: v, min: lo, max: hi });
      } else if ((m = RE.init.exec(s))) init[m[1]] = +m[2];
      else if ((m = RE.range.exec(s))) ranges[m[1]] = [+m[2], +m[3]];
      else if ((m = RE.noise.exec(s))) noises.push({ v: m[1], src: m[2], ln: ln });
      else if ((m = RE.aux.exec(s))) auxes.push({ v: m[1], src: m[2], ln: ln });
      else if ((m = RE.ode.exec(s))) {
        if (kind === "map") throw new FormulaError("Differential and difference equations cannot be mixed", ln);
        kind = "ode"; eqs.push({ v: m[1] || m[2], src: m[3], ln: ln });
      } else if ((m = RE.map.exec(s))) {
        if (kind === "ode") throw new FormulaError("Differential and difference equations cannot be mixed", ln);
        kind = "map"; eqs.push({ v: m[1], src: m[2], ln: ln });
      } else throw new FormulaError("Cannot read '" + s + "'", ln);
    });
    if (!eqs.length) throw new FormulaError("No equations: write for example x' = -x");

    const vars = eqs.map(function (e) { return e.v; });
    const varIndex = {};
    vars.forEach(function (v, i) {
      if (v in varIndex) throw new FormulaError("Variable '" + v + "' has two equations", eqs[i].ln);
      if (v in RESERVED) throw new FormulaError("'" + v + "' is a reserved name", eqs[i].ln);
      varIndex[v] = i;
    });
    params.forEach(function (q) {
      if (q.name in varIndex) throw new FormulaError("'" + q.name + "' is both a variable and a parameter");
      if (q.name in RESERVED) throw new FormulaError("'" + q.name + "' is a reserved name");
    });

    // Local names: v_<var>, p_<param>, a_<aux>; the scope maps formula names to them.
    const scope = {};
    vars.forEach(function (v) { scope[v] = "v_" + v; });
    params.forEach(function (q) { scope[q.name] = "p_" + q.name; });
    const info = { varIndex: varIndex, usesLag: false, usesTime: false, usesRandom: false, lagExprs: [] };
    const head = vars.map(function (v, i) { return "const v_" + v + " = x[" + i + "];"; }).join(" ") + "\n" +
      params.map(function (q, i) { return "const p_" + q.name + " = p[" + i + "];"; }).join(" ") + "\n";
    let auxCode = "";
    auxes.forEach(function (a) {
      if (a.v in scope || a.v in RESERVED) throw new FormulaError("'" + a.v + "' is already defined", a.ln);
      auxCode += "const a_" + a.v + " = " + emit(parse(a.src, a.ln), scope, a.ln, info) + ";\n";
      scope[a.v] = "a_" + a.v;
    });
    let body = "";
    eqs.forEach(function (e, i) { body += "dx[" + i + "] = " + emit(parse(e.src, e.ln), scope, e.ln, info) + ";\n"; });
    const usesLag = info.usesLag;

    let gBody = "";
    const noiseMask = vars.map(function () { return 0; });
    noises.forEach(function (nz) {
      if (!(nz.v in varIndex)) throw new FormulaError("noise refers to unknown variable '" + nz.v + "'", nz.ln);
      noiseMask[varIndex[nz.v]] = 1;
      gBody += "gx[" + varIndex[nz.v] + "] = " + emit(parse(nz.src, nz.ln), scope, nz.ln, info) + ";\n";
    });
    if (noises.length && kind === "map") throw new FormulaError("noise lines apply to differential equations; add the noise inside the map instead");

    const f = new Function(PRELUDE + "return function (t, x, p, dx, H, R) {\nR = R || __R0;\n" + head + auxCode + body + "};")();
    const g = noises.length ? new Function(PRELUDE + "return function (t, x, p, gx, H, R) {\nR = R || __R0;\n" + head + auxCode +
      "for (let i = 0; i < gx.length; i++) gx[i] = 0;\n" + gBody + "};")() : null;

    // Constant delays are needed to size the history buffer; a delay that
    // depends on parameters is evaluated at the current parameters.
    const lagFns = info.lagExprs.map(function (src) { return new Function("p", "t", head.split("\n")[1] + "\nreturn " + src + ";"); });

    const k = usesLag ? "dde" : noises.length ? "sde" : kind;
    if (usesLag && kind === "map") throw new FormulaError("lag() is available in differential equations only");
    return {
      kind: k, time: kind === "map" ? "discrete" : "continuous",
      vars: vars, params: params,
      init: vars.map(function (v) { return v in init ? init[v] : 0.1; }),
      ranges: ranges, f: f, g: g, noiseMask: noiseMask, usesTime: info.usesTime, usesRandom: info.usesRandom,
      // A delay that depends on the state or on an aux cannot be bounded in
      // advance; it returns Infinity and the history then uses its full capacity.
      maxLag: function (p) {
        let m = 0;
        lagFns.forEach(function (fn) { let v; try { v = fn(p, 0); } catch (e) { v = Infinity; } m = Math.max(m, isNaN(v) ? Infinity : Math.abs(v)); });
        return m;
      },
      source: String(text)
    };
  }

  // ---------------------------------------------------------------- LaTeX
  const GREEK = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi",
    "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi", "Omega"];
  const GREEK_ALIAS = { eps: "varepsilon", lam: "lambda", gammag: "gamma", sig: "sigma" };

  // x1 -> x_{1}, a12 -> a_{12}, alpha_2 -> \alpha_{2}, K_m -> K_{m}, omega -> \omega
  function texName(name) {
    let base = name, sub = "";
    const us = name.indexOf("_");
    if (us > 0) { base = name.slice(0, us); sub = name.slice(us + 1); }
    else {
      const m = /^([A-Za-z]+?)(\d+)$/.exec(name);
      if (m) { base = m[1]; sub = m[2]; }
    }
    if (base in GREEK_ALIAS) base = GREEK_ALIAS[base];
    let b = GREEK.indexOf(base) >= 0 ? "\\" + base : base.length > 1 ? "\\mathrm{" + base + "}" : base;
    if (sub) b += "_{" + (GREEK.indexOf(sub) >= 0 ? "\\" + sub : sub.length > 1 && !/^\d+$/.test(sub) ? "\\mathrm{" + sub + "}" : sub) + "}";
    return b;
  }

  // \dot over the base symbol only: x1 -> \dot{x}_{1}
  function dotted(name) {
    const t = texName(name), m = /^(\\?[A-Za-z]+|\\mathrm\{[^}]*\})(_\{.*\})?$/.exec(t);
    return m ? "\\dot{" + m[1] + "}" + (m[2] || "") : "\\dot{" + t + "}";
  }

  const TEX_PREC = { "||": 1, "&&": 2, "==": 3, "!=": 3, "<": 4, "<=": 4, ">": 4, ">=": 4, "+": 5, "-": 5, "*": 6, "/": 7, "^": 8 };
  function tex(node, ctx) {
    const wrap = function (child, minPrec) {
      const s = tex(child, ctx);
      const pr = child.k === "bin" ? TEX_PREC[child.op] : child.k === "neg" ? 5.4 : child.k === "call" && child.f === "lag" ? 8.5 : 9;
      return pr < minPrec ? "\\left(" + s + "\\right)" : s;
    };
    switch (node.k) {
      case "num": return String(node.v).replace(/e([+-]?\d+)$/, " \\times 10^{$1}");
      case "name":
        if (node.v === "t") return "t";
        if (node.v === "pi" && !(node.v in ctx.declared)) return "\\pi";
        if (node.v in ctx.aux) return texName(node.v);
        return texName(node.v);
      case "neg": return "-" + wrap(node.a, 6);
      case "not": return "\\neg " + wrap(node.a, 9);
      case "bin": {
        const op = node.op;
        if (op === "/") return "\\frac{" + tex(node.a, ctx) + "}{" + tex(node.b, ctx) + "}";
        if (op === "^") return wrap(node.a, 9) + "^{" + tex(node.b, ctx) + "}";
        if (op === "*") {
          // A negated left factor needs no brackets: -a*x = -(a*x).
          const a = node.a.k === "neg" ? tex(node.a, ctx) : wrap(node.a, 6), b = wrap(node.b, 6);
          const numRight = node.b.k === "num";
          return a + (numRight || /^[\d.]/.test(b) ? " \\cdot " : "\\,") + b;
        }
        const sym = { "+": " + ", "-": " - ", "<": " < ", "<=": " \\le ", ">": " > ", ">=": " \\ge ", "==": " = ", "!=": " \\ne ", "&&": " \\land ", "||": " \\lor " }[op];
        return wrap(node.a, TEX_PREC[op]) + sym + wrap(node.b, op === "-" || op === "+" ? TEX_PREC[op] + 0.5 : TEX_PREC[op]);
      }
      case "call": {
        const a = node.args.map(function (x) { return tex(x, ctx); });
        switch (node.f) {
          case "sqrt": return "\\sqrt{" + a[0] + "}";
          case "cbrt": return "\\sqrt[3]{" + a[0] + "}";
          case "abs": return "\\left|" + a[0] + "\\right|";
          case "exp": return "e^{" + a[0] + "}";
          case "pow": return wrap(node.args[0], 9) + "^{" + a[1] + "}";
          case "lag": return tex(node.args[0], ctx) + "(t - " + a[1] + ")";
          case "hill": return "\\frac{" + wrap(node.args[0], 9) + "^{" + a[2] + "}}{" + wrap(node.args[1], 9) + "^{" + a[2] + "} + " + wrap(node.args[0], 9) + "^{" + a[2] + "}}";
          case "step": case "heaviside": return "\\Theta\\left(" + a[0] + "\\right)";
          case "floor": return "\\lfloor " + a[0] + " \\rfloor";
          case "mod": return a[0] + " \\bmod " + wrap(node.args[1], 9);
          case "ifelse": return "\\begin{cases} " + a[1] + " & " + a[0] + " \\\\ " + a[2] + " & \\text{otherwise} \\end{cases}";
          case "clamp": return "\\mathrm{clamp}\\left(" + a.join(", ") + "\\right)";
          case "sin": case "cos": case "tan": case "sinh": case "cosh": case "tanh": case "log": case "min": case "max":
            return "\\" + node.f + "\\left(" + a.join(", ") + "\\right)";
          case "asin": case "acos": case "atan": return "\\" + node.f.replace("a", "arc") + "\\left(" + a[0] + "\\right)";
          default: return "\\mathrm{" + node.f + "}\\left(" + a.join(", ") + "\\right)";
        }
      }
    }
    return "";
  }

  /* LaTeX for every statement of a system: equations, noise terms and aux
     definitions, in source order. Parameters are listed separately by callers. */
  function systemLatex(text) {
    const out = [], declared = {}, aux = {};
    let kind = "ode";
    const lines = String(text).split(/\r?\n/).map(function (r) { return r.replace(/#.*$/, "").trim(); }).filter(Boolean);
    lines.forEach(function (s) { let m; if ((m = RE.param.exec(s))) declared[m[1]] = 1; if (RE.map.test(s)) kind = "map"; });
    const noiseOf = {};
    lines.forEach(function (s) { const m = RE.noise.exec(s); if (m) noiseOf[m[1]] = m[2]; });
    const ctx = { declared: declared, aux: aux };
    lines.forEach(function (s) {
      let m;
      if ((m = RE.aux.exec(s))) { aux[m[1]] = 1; out.push(texName(m[1]) + " = " + tex(parse(m[2]), ctx)); }
      else if ((m = RE.ode.exec(s))) {
        const v = m[1] || m[2], rhs = tex(parse(m[3]), ctx);
        if (noiseOf[v] !== undefined) out.push("\\mathrm{d}" + texName(v) + " = \\left(" + rhs + "\\right)\\mathrm{d}t + " + tex(parse(noiseOf[v]), ctx) + "\\,\\mathrm{d}W_{" + texName(v) + "}");
        else out.push(dotted(v) + " = " + rhs);
      } else if ((m = RE.map.exec(s))) out.push(texName(m[1]).replace(/_\{(.*)\}$/, "_{$1,\\,n+1}").replace(/^([^_]*)$/, "$1_{n+1}") + " = " + tex(parse(m[2]), ctx));
    });
    return { lines: out, kind: kind };
  }

  DF.FormulaError = FormulaError;
  DF.STATEMENT = RE;
  DF.GREEK = GREEK;
  DF.GREEK_ALIAS = GREEK_ALIAS;
  DF.TEX_PREC = TEX_PREC;
  DF.parseExpression = parse;
  DF.compileSystem = compileSystem;
  DF.systemLatex = systemLatex;
  DF.texName = texName;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/core/sim.js
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

// ---- src/core/analysis.js
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

// ---- src/render/style.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Themes and palettes. Themes set the background and the ink of axes and
   text; palettes set the colours of trajectories. The RElab palettes are
   the brand palettes of the lab (indigo #170C3A and #2a1766, orange #EE6A24,
   amber #FB9E07, pink #A52C60); blackboard follows the ink set used for dark
   slides; the continuous ramps are sampled from the published colour maps. */
(function (DF) {
  "use strict";

  const THEMES = {
    "relab-night": { label: "RElab night", bg: ["radial", "#2a1766", "#170C3A", "#0b0620"], ink: "#eceaf4", muted: "rgba(236,234,244,0.62)", grid: "rgba(236,234,244,0.12)", blend: "lighter", dark: true },
    "blackboard": { label: "Blackboard", bg: ["solid", "#000000"], ink: "#ebebeb", muted: "#bfbfbf", grid: "#404040", blend: "lighter", dark: true },
    "deep-sea": { label: "Deep sea", bg: ["radial", "#0f2a3d", "#07131f", "#02070c"], ink: "#e6f0f5", muted: "rgba(230,240,245,0.6)", grid: "rgba(230,240,245,0.12)", blend: "lighter", dark: true },
    "graphite": { label: "Graphite", bg: ["solid", "#1b1d22"], ink: "#e6e6e6", muted: "#a8a8a8", grid: "#3a3d44", blend: "lighter", dark: true },
    "paper": { label: "Paper", bg: ["solid", "#fffdf8"], ink: "#212529", muted: "#687078", grid: "#e4e1d8", blend: "source-over", dark: false },
    "white": { label: "White", bg: ["solid", "#ffffff"], ink: "#262626", muted: "#666666", grid: "#ebebeb", blend: "source-over", dark: false },
    "transparent": { label: "Transparent", bg: ["none"], ink: "#888888", muted: "#888888", grid: "rgba(128,128,128,0.2)", blend: "source-over", dark: true }
  };

  // Discrete palettes: one colour per variable, member or branch.
  const PALETTES = {
    "relab": { label: "RElab", colors: ["#EE6A24", "#FB9E07", "#A52C60", "#CF4446", "#F6D645", "#764BA2", "#3093CF", "#eceaf4"] },
    "relab-qualitative": { label: "RElab qualitative", colors: ["#FB9E07", "#4777ef", "#009E73", "#764BA2", "#E85D04", "#00B4D8", "#D55E00", "#CC79A7"] },
    "blackboard": { label: "Blackboard ink", colors: ["#5ec5ff", "#ff8c42", "#ff6b6b", "#3ddc97", "#a8d05b", "#ffd166", "#9d8df1", "#ff6ec7", "#d4a373", "#c0c0c0"] },
    "tableau": { label: "Tableau 10", colors: ["#4e79a7", "#f28e2c", "#e15759", "#76b7b2", "#59a14f", "#edc949", "#af7aa1", "#ff9da7", "#9c755f", "#bab0ab"] },
    "okabe-ito": { label: "Okabe-Ito", colors: ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"] },
    "ink": { label: "Ink", colors: ["#170C3A", "#2a1766", "#A52C60", "#EE6A24", "#687078"] },
    "mono-amber": { label: "Amber", colors: ["#FB9E07", "#f5b04a", "#EE6A24", "#ffd699"] },
    "mono-ice": { label: "Ice", colors: ["#9fd8ff", "#5ec5ff", "#d9f1ff", "#3093CF"] }
  };

  // Continuous ramps, as colour stops sampled evenly on [0, 1].
  const RAMPS = {
    "relab-fire": { label: "RElab fire", stops: ["#170C3A", "#2a1766", "#A52C60", "#EE6A24", "#FB9E07", "#F6D645"] },
    "relab-sequential": { label: "RElab sequential", stops: ["#F5F0E6", "#FFE5B4", "#FFD699", "#FBC66A", "#FB9E07", "#E88507", "#D46B07", "#B85507", "#9C4007"] },
    "relab-diverging": { label: "RElab diverging", stops: ["#170C3A", "#2a1766", "#764BA2", "#B8A9C9", "#F5F0E6", "#FFD699", "#FB9E07", "#E85D04", "#C7380B"] },
    "blackboard": { label: "Blackboard ramp", stops: ["#5ec5ff", "#ff6ec7", "#ffd166"] },
    "viridis": { label: "Viridis", stops: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"] },
    "magma": { label: "Magma", stops: ["#000004", "#1c1044", "#4f127b", "#812581", "#b5367a", "#e55064", "#fb8761", "#fec287", "#fcfdbf"] },
    "inferno": { label: "Inferno", stops: ["#000004", "#1f0c48", "#550f6d", "#88226a", "#ba3655", "#e35933", "#f98e09", "#f8c932", "#fcffa4"] },
    "mako": { label: "Mako", stops: ["#0b0405", "#2b1c35", "#3e356b", "#3b5698", "#357ba3", "#38a0ab", "#4bc4ad", "#8ae1b9", "#def5e5"] },
    "sand-red": { label: "Sand red", stops: ["#ffffff", "#fbe3c3", "#f3b27a", "#e0703f", "#b8321f", "#8f0b12"] }
  };

  function hexToRgb(hex) {
    const h = hex.replace("#", "");
    const v = parseInt(h.length === 3 ? h.split("").map(function (c) { return c + c; }).join("") : h, 16);
    return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
  }
  const rgbCache = {};
  function rgb(hex) { return rgbCache[hex] || (rgbCache[hex] = hexToRgb(hex)); }

  // Colour of u in [0, 1] on a ramp, as an [r, g, b] triple.
  function rampRGB(name, u) {
    const st = (RAMPS[name] || RAMPS["relab-fire"]).stops;
    u = Math.min(1, Math.max(0, isFinite(u) ? u : 0)) * (st.length - 1);
    const i = Math.min(st.length - 2, Math.floor(u)), f = u - i;
    const a = rgb(st[i]), b = rgb(st[i + 1]);
    return [Math.round(a[0] + f * (b[0] - a[0])), Math.round(a[1] + f * (b[1] - a[1])), Math.round(a[2] + f * (b[2] - a[2]))];
  }

  // Paint the background of a theme onto a 2D context.
  function paintBackground(ctx, w, h, theme) {
    const th = THEMES[theme] || THEMES["relab-night"], bg = th.bg;
    if (bg[0] === "none") { ctx.clearRect(0, 0, w, h); return; }
    if (bg[0] === "solid") { ctx.fillStyle = bg[1]; ctx.fillRect(0, 0, w, h); return; }
    const g = ctx.createRadialGradient(w * 0.5, h * 0.5, 0, w * 0.5, h * 0.5, Math.hypot(w, h) * 0.6);
    g.addColorStop(0, bg[1]); g.addColorStop(0.55, bg[2]); g.addColorStop(1, bg[3]);
    ctx.fillStyle = g; ctx.fillRect(0, 0, w, h);
  }
  function backgroundCSS(theme) {
    const bg = (THEMES[theme] || THEMES["relab-night"]).bg;
    if (bg[0] === "none") return "transparent";
    if (bg[0] === "solid") return bg[1];
    return "radial-gradient(ellipse at 50% 50%, " + bg[1] + " 0%, " + bg[2] + " 55%, " + bg[3] + " 100%)";
  }

  DF.THEMES = THEMES;
  DF.PALETTES = PALETTES;
  DF.RAMPS = RAMPS;
  DF.rgb = rgb;
  DF.rampRGB = rampRGB;
  DF.paintBackground = paintBackground;
  DF.backgroundCSS = backgroundCSS;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/render/frame.js
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

  /* Sample of the attractor: one orbit from the initial state, without
     noise, run for a transient and then stored at `size` evenly spaced
     points. The run spends at most `budget` steps, 30 percent of them on
     the transient, unless opts fixes `transient` and `span` in model time.
     Returns null when the orbit diverges or settles on an equilibrium, so
     that callers fall back to births over the axis box. */
  function attractorSample(sys, params, init, opts) {
    opts = opts || {};
    const dim = sys.vars.length, dt = opts.dt || 0.01, size = opts.size || 3000, budget = opts.budget || 160000;
    const nT = opts.transient > 0 ? Math.round(opts.transient / dt) : Math.round(0.3 * budget);
    const nS = opts.span > 0 ? Math.round(opts.span / dt) : budget - nT;
    const every = Math.max(1, Math.floor(nS / size));
    let sim;
    try { sim = new DF.Simulator(sys, { dt: dt, params: params, init: init, n: 1, seed: 2468, deterministic: true, perturbations: opts.perturbations || [] }); }
    catch (e) { return null; }
    for (let s = 0; s < nT; s++) { sim.step(); if (!sim.alive[0]) return null; }
    const out = new Float64Array(size * dim), lo = new Float64Array(dim).fill(Infinity), hi = new Float64Array(dim).fill(-Infinity);
    let m = 0;
    for (let s = 0; s < nS && m < size; s++) {
      sim.step();
      if (!sim.alive[0]) return null;
      if (s % every) continue;
      for (let i = 0; i < dim; i++) {
        const v = sim.X[i];
        if (!isFinite(v)) return null;
        out[m * dim + i] = v; if (v < lo[i]) lo[i] = v; if (v > hi[i]) hi[i] = v;
      }
      m++;
    }
    if (m < 50) return null;
    // An equilibrium or a vanishing orbit: every coordinate barely moves
    // against the scale of the variable.
    let moving = false;
    for (let i = 0; i < dim; i++) if (hi[i] - lo[i] > 1e-4 * Math.max(1e-9, Math.abs(hi[i]) + Math.abs(lo[i]))) moving = true;
    if (!moving) return null;
    return { dim: dim, n: m, X: out.subarray(0, m * dim), lo: Array.from(lo), hi: Array.from(hi) };
  }
  // Ranges of the sample, with a margin of `pad` of each width.
  function sampleRanges(sample, pad) {
    pad = pad === undefined ? 0.08 : pad;
    return sample.lo.map(function (l, i) {
      const h = sample.hi[i], w = h - l || Math.max(Math.abs(l), 1) * 0.1;
      return [l - pad * w, h + pad * w];
    });
  }

  /* Eigen-decomposition of a symmetric 3 x 3 matrix by cyclic Jacobi
     rotations; returns eigenvalues in decreasing order and the eigenvectors
     as the rows of V. */
  function symEigen3(A) {
    const a = [A[0].slice(), A[1].slice(), A[2].slice()], V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 50; sweep++) {
      const off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
      if (off < 1e-30) break;
      for (let p = 0; p < 2; p++) for (let q = p + 1; q < 3; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const th = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]), c = Math.cos(th), s = Math.sin(th);
        for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 3; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    const order = [0, 1, 2].sort(function (i, j) { return a[j][j] - a[i][i]; });
    return { values: order.map(function (i) { return a[i][i]; }), vectors: order.map(function (i) { return [V[0][i], V[1][i], V[2][i]]; }) };
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
    } else if (kind === "flow" && v.spawn === "attractor" && o.sample) {
      // Immortal particles born on the attractor: 24 members at evenly
      // spaced points of the sample, so that slow and fast stretches of the
      // attractor are weighted as the ensemble on screen weights them.
      n = 24;
      sim = new DF.Simulator(sys, Object.assign({}, common, { n: n, perturbations: o.perturbations || [] }));
      const S = o.sample;
      for (let k = 0; k < n; k++) { const j = Math.floor((k + 0.5) * S.n / n); sim.setMember(k, S.X.subarray(j * S.dim, (j + 1) * S.dim)); }
    } else {
      n = 6;
      sim = new DF.Simulator(sys, Object.assign({}, common, { n: n, spread: kind === "flow" ? o.spread : 0.01 * meanW, initMode: kind === "flow" ? o.initMode : "point", perturbations: o.perturbations || [] }));
    }
    const warmFrames = v.warmup === undefined ? (kind === "flow" && !(v.spawn === "attractor" && o.sample) ? 60 : 0) : v.warmup;
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
    const onAttractor = type === "trajectory" || type === "density" || type === "strobe" || (type === "flow" && (v.life === "inf" || v.spawn === "init" || v.spawn === "attractor"));
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
      } else if (type === "flow") {
        // The replayed speed m(r) grows with the rate r but need not be
        // proportional to it: when half-second windows reach a fast stretch
        // only above some rate (the spikes of a bursting neuron), m jumps.
        // The target is bracketed by factors of 3 from r0 and then located by
        // bisection on log r, within a factor of 20 of r0.
        const lo0 = r0 / 20, hi0 = r0 * 20, m0 = replaySpeed(sys, o, r0, "flow");
        if (m0 > 1e-9) {
          let a = r0, b = r0, ma = m0, mb = m0;
          if (m0 < target) { while (mb < target && b < hi0) { a = b; ma = mb; b = Math.min(hi0, b * 3); mb = replaySpeed(sys, o, b, "flow"); } }
          else { while (ma > target && a > lo0) { b = a; mb = ma; a = Math.max(lo0, a / 3); ma = replaySpeed(sys, o, a, "flow"); } }
          r = Math.abs(Math.log(ma / target)) < Math.abs(Math.log(mb / target)) ? a : b;
          if ((ma - target) * (mb - target) < 0) {
            for (let it = 0; it < 6; it++) {
              const c = Math.sqrt(a * b), mc = replaySpeed(sys, o, c, "flow");
              if (!(mc > 1e-9)) break;
              if (mc < target) a = c; else b = c;
              r = c;
              if (Math.abs(mc / target - 1) < 0.08) break;
            }
          }
          r = clamp(r); basis = "replay";
        }
      } else if (type === "phase") {
        const bound = 4;
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
    if (this.frame) return this.projectFitted(x, out);
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
  /* Fitted 3D frame. The displayed coordinates, each scaled to [-1, 1] by
     its range, are centred on the mean of an attractor sample and rotated
     onto its principal axes: the widest direction is horizontal, the
     second vertical and the thinnest along the line of sight, so that at
     azimuth and elevation 0 the reader faces the broadest side of the
     attractor. The scale fits the extent of the projected sample over the
     azimuths in [azim - sweep, azim + sweep] to the plot box. */
  Camera.prototype.fitTo = function (sample, opts) {
    opts = opts || {};
    const dim = sample.dim, n = sample.n, ax = this.axes, self = this;
    const c = [0, 0, 0], C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], q = [0, 0, 0];
    const qOf = function (k) { for (let j = 0; j < 3; j++) q[j] = self.norm(j, sample.X[k * dim + ax[j]]); return q; };
    for (let k = 0; k < n; k++) { qOf(k); for (let j = 0; j < 3; j++) c[j] += q[j] / n; }
    for (let k = 0; k < n; k++) {
      qOf(k);
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] += (q[i] - c[i]) * (q[j] - c[j]) / n;
    }
    const E = symEigen3(C), R = E.vectors;
    // Right-handed frame: horizontal x vertical points along the line of sight.
    const cr = [R[0][1] * R[1][2] - R[0][2] * R[1][1], R[0][2] * R[1][0] - R[0][0] * R[1][2], R[0][0] * R[1][1] - R[0][1] * R[1][0]];
    if (cr[0] * R[2][0] + cr[1] * R[2][1] + cr[2] * R[2][2] < 0) R[2] = R[2].map(function (v) { return -v; });
    this.frame = { c: c, R: R, var: E.values, ox: 0, oy: 0, ex: 1, ey: 1 };
    const sweep = opts.sweep === undefined ? 0.6 : opts.sweep, azim = this.azim, elev = this.elev;
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    const steps = sweep >= Math.PI ? 24 : 9, stride = Math.max(1, Math.floor(n / 1500)), u = [0, 0, 0, 0];
    for (let a = 0; a < steps; a++) {
      const az = sweep >= Math.PI ? 2 * Math.PI * a / steps : azim + sweep * (2 * a / (steps - 1) - 1);
      for (let k = 0; k < n; k += stride) {
        this.rawFitted(qOf(k), az, elev, u);
        if (u[0] < x0) x0 = u[0]; if (u[0] > x1) x1 = u[0]; if (u[1] < y0) y0 = u[1]; if (u[1] > y1) y1 = u[1];
      }
    }
    this.frame.ox = (x0 + x1) / 2; this.frame.oy = (y0 + y1) / 2;
    this.frame.ex = Math.max(1e-9, (x1 - x0) / 2); this.frame.ey = Math.max(1e-9, (y1 - y0) / 2);
    return this.frame;
  };
  // Screen-plane coordinates (X right, Y up) and depth of normalised q.
  Camera.prototype.rawFitted = function (q, azim, elev, out) {
    const F = this.frame, d0 = q[0] - F.c[0], d1 = q[1] - F.c[1], d2 = q[2] - F.c[2], R = F.R;
    const qa = R[0][0] * d0 + R[0][1] * d1 + R[0][2] * d2, qc = R[1][0] * d0 + R[1][1] * d1 + R[1][2] * d2, qb = R[2][0] * d0 + R[2][1] * d1 + R[2][2] * d2;
    const ca = Math.cos(azim), sa = Math.sin(azim), ce = Math.cos(elev), se = Math.sin(elev);
    const X = ca * qa - sa * qb, Y0 = sa * qa + ca * qb;
    out[0] = X; out[1] = ce * qc - se * Y0; out[2] = ce * Y0 + se * qc;
    return out;
  };
  Camera.prototype.projectFitted = function (x, out) {
    const F = this.frame, b = this.plotBox(), q = this.fq || (this.fq = [0, 0, 0]), u = this.fu || (this.fu = [0, 0, 0]);
    for (let j = 0; j < 3; j++) q[j] = this.norm(j, x[this.axes[j]]);
    this.rawFitted(q, this.azim, this.elev, u);
    const s = Math.min(b.w / F.ex, b.h / F.ey) * 0.5 * (1 - this.pad) * this.zoom;
    out[0] = b.x + b.w / 2 + s * (u[0] - F.ox);
    out[1] = b.y + b.h / 2 - s * (u[1] - F.oy);
    out.depth = u[2];
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
  DF.attractorSample = attractorSample;
  DF.sampleRanges = sampleRanges;
  DF.symEigen3 = symEigen3;
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

// ---- src/render/mathtype.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Typeset equations without KaTeX. The statements of a system are laid out
   from their parse trees as boxes (width w, ascent a, descent d) that hold
   glyph runs, rules, strokes and dots, and the boxes are painted onto a
   canvas or written as SVG. PNG, GIF and WebM frames, SVG files and pages
   without KaTeX therefore show the same equations as the studio: dotted
   derivatives, fractions, powers, subscripts and Greek letters. The layout
   follows the precedence rules of systemLatex in src/core/expr.js. */
(function (DF) {
  "use strict";

  const GLYPH = {
    alpha: "\u03b1", beta: "\u03b2", gamma: "\u03b3", delta: "\u03b4", epsilon: "\u03f5", varepsilon: "\u03b5", zeta: "\u03b6", eta: "\u03b7",
    theta: "\u03b8", iota: "\u03b9", kappa: "\u03ba", lambda: "\u03bb", mu: "\u03bc", nu: "\u03bd", xi: "\u03be", pi: "\u03c0", rho: "\u03c1",
    sigma: "\u03c3", tau: "\u03c4", upsilon: "\u03c5", phi: "\u03d5", chi: "\u03c7", psi: "\u03c8", omega: "\u03c9",
    Gamma: "\u0393", Delta: "\u0394", Theta: "\u0398", Lambda: "\u039b", Xi: "\u039e", Pi: "\u03a0", Sigma: "\u03a3", Phi: "\u03a6", Psi: "\u03a8", Omega: "\u03a9"
  };
  const MATH_FONT = "KaTeX_Math, 'TeX Gyre Pagella', 'Palatino Linotype', Palatino, 'Times New Roman', serif";
  const MAIN_FONT = "KaTeX_Main, 'TeX Gyre Pagella', 'Palatino Linotype', Palatino, 'Times New Roman', serif";
  const FUNC_NAME = { asin: "arcsin", acos: "arccos", atan: "arctan", log: "ln", log10: "log\u2081\u2080", log2: "log\u2082" };

  let measureCtx = null;
  function fontOf(size, it) { return (it ? "italic " : "") + size.toFixed(2) + "px " + (it ? MATH_FONT : MAIN_FONT); }
  function measure(str, size, it) {
    if (!measureCtx && typeof document !== "undefined") measureCtx = document.createElement("canvas").getContext("2d");
    if (!measureCtx) return 0.55 * size * str.length;
    measureCtx.font = fontOf(size, it);
    return measureCtx.measureText(str).width;
  }

  // ------------------------------------------------------------ boxes
  function text(str, size, it) {
    return { w: measure(str, size, it), a: 0.72 * size, d: 0.22 * size, items: [{ k: "t", s: str, x: 0, y: 0, size: size, it: !!it }] };
  }
  function empty(w) { return { w: w || 0, a: 0, d: 0, items: [] }; }
  // Items carry their offset in x and y; the points of a stroke are relative to it.
  function moved(item, dx, dy) {
    const c = Object.assign({}, item);
    c.x += dx; c.y += dy;
    return c;
  }
  function place(into, box, dx, dy) { box.items.forEach(function (it) { into.items.push(moved(it, dx, dy)); }); }
  function row(boxes) {
    const out = empty(0);
    boxes.forEach(function (b) {
      if (!b) return;
      place(out, b, out.w, 0);
      out.w += b.w; out.a = Math.max(out.a, b.a); out.d = Math.max(out.d, b.d);
    });
    return out;
  }
  function sup(base, e, size) {
    const rise = Math.max(0.42 * size, base.a - 0.55 * e.a), out = empty(0);
    place(out, base, 0, 0); place(out, e, base.w + 0.04 * size, -rise);
    out.w = base.w + 0.04 * size + e.w; out.a = Math.max(base.a, rise + e.a); out.d = Math.max(base.d, e.d - rise);
    return out;
  }
  function sub(base, s, size) {
    const drop = 0.22 * size, out = empty(0);
    place(out, base, 0, 0); place(out, s, base.w + 0.02 * size, drop);
    out.w = base.w + 0.02 * size + s.w; out.a = Math.max(base.a, s.a - drop); out.d = Math.max(base.d, drop + s.d);
    return out;
  }
  function frac(num, den, size) {
    const ax = 0.25 * size, gap = 0.14 * size, th = Math.max(1, 0.05 * size), w = Math.max(num.w, den.w) + 0.24 * size, out = empty(w);
    place(out, num, (w - num.w) / 2, -ax - gap - num.d);
    place(out, den, (w - den.w) / 2, -ax + gap + den.a);
    out.items.push({ k: "r", x: 0.06 * size, y: -ax - th / 2, w: w - 0.12 * size, h: th });
    out.a = ax + gap + num.d + num.a; out.d = Math.max(0, -ax + gap + den.a + den.d);
    return out;
  }
  // Round, bar or floor fences drawn as strokes that span the content.
  function fence(box, open, close, size) {
    const top = -Math.max(box.a, 0.72 * size) - 0.06 * size, bot = Math.max(box.d, 0.22 * size) + 0.06 * size, H = bot - top;
    const lw = Math.max(0.8, 0.06 * size), fw = open === "|" ? 0.22 * size : 0.3 * size + 0.04 * H;
    const glyph = function (kind, left) {
      const g = empty(fw);
      if (kind === "(") {
        const pts = [];
        for (let k = 0; k <= 16; k++) {
          const u = k / 16, ang = Math.PI * (u - 0.5);
          const x = fw * (0.8 - 0.55 * Math.cos(ang));
          pts.push(left ? x : fw - x, top + u * H);
        }
        g.items.push({ k: "p", x: 0, y: 0, pts: pts, lw: lw });
      } else if (kind === "|") g.items.push({ k: "p", x: 0, y: 0, pts: [fw / 2, top, fw / 2, bot], lw: lw });
      else if (kind === "floor") g.items.push({ k: "p", x: 0, y: 0, pts: left ? [fw * 0.35, top, fw * 0.35, bot, fw * 0.8, bot] : [fw * 0.65, top, fw * 0.65, bot, fw * 0.2, bot], lw: lw });
      g.a = -top; g.d = bot;
      return g;
    };
    return row([glyph(open === "|" ? "|" : open === "floor" ? "floor" : "(", true), box, glyph(close === "|" ? "|" : close === "floor" ? "floor" : "(", false)]);
  }
  function radical(box, size, index) {
    const top = -box.a - 0.14 * size, bot = box.d + 0.04 * size, lw = Math.max(0.8, 0.06 * size), lead = 0.55 * size;
    const out = empty(0);
    out.items.push({ k: "p", x: 0, y: 0, pts: [0.02 * size, (top + bot) * 0.5 + 0.1 * size, 0.14 * size, (top + bot) * 0.5, 0.3 * size, bot, 0.5 * size, top, lead + box.w + 0.08 * size, top], lw: lw });
    place(out, box, lead, 0);
    if (index) place(out, index, 0.02 * size, (top + bot) * 0.5 - 0.1 * size);
    out.w = lead + box.w + 0.1 * size; out.a = -top + lw; out.d = bot;
    return out;
  }
  function dot(box, size) {
    const out = empty(box.w);
    place(out, box, 0, 0);
    out.items.push({ k: "c", x: box.w / 2 + 0.07 * size, y: -box.a - 0.14 * size, r: 0.07 * size });
    out.a = box.a + 0.22 * size; out.d = box.d;
    return out;
  }

  // ------------------------------------------------------------ names
  // x1 -> x_1, alpha_2 -> alpha_2, K_m -> K_m, omega -> omega, N0 -> N_0
  function nameParts(name) {
    let base = name, subs = "";
    const us = name.indexOf("_");
    if (us > 0) { base = name.slice(0, us); subs = name.slice(us + 1); }
    else { const m = /^([A-Za-z]+?)(\d+)$/.exec(name); if (m) { base = m[1]; subs = m[2]; } }
    if (base in DF.GREEK_ALIAS) base = DF.GREEK_ALIAS[base];
    return { base: base, sub: subs };
  }
  function symbol(str, size) {
    if (str in GLYPH) return text(GLYPH[str], size, /^[a-z]/.test(str));
    if (str.length === 1) return text(str, size, true);
    if (/^\d+$/.test(str)) return text(str, size, false);
    return text(str, size, false);
  }
  function nameBox(name, size, opts) {
    opts = opts || {};
    const pr = nameParts(name);
    let b = symbol(pr.base, size);
    if (opts.dot) b = dot(b, size);
    const parts = [];
    if (pr.sub) parts.push(symbol(pr.sub, 0.7 * size));
    if (opts.extraSub) parts.push(opts.extraSub(0.7 * size));
    if (parts.length) b = sub(b, parts.length === 1 ? parts[0] : row([parts[0], text(",", 0.7 * size, false), parts[1]]), size);
    return b;
  }

  // ------------------------------------------------------------- trees
  const OPS = { "+": "+", "-": "\u2212", "<": "<", "<=": "\u2264", ">": ">", ">=": "\u2265", "==": "=", "!=": "\u2260", "&&": "\u2227", "||": "\u2228" };
  function prec(node) {
    if (node.k === "bin") return DF.TEX_PREC[node.op];
    if (node.k === "neg") return 5.4;
    if (node.k === "call" && node.f === "lag") return 8.5;
    return 9;
  }
  function number(v, size) {
    const str = String(v), m = /^(-?[\d.]+)e([+-]?\d+)$/.exec(str);
    if (!m) return text(str.replace("-", "\u2212"), size, false);
    const mant = m[1] === "1" ? null : text(m[1] + "\u00d7", size, false);
    return row([mant, sup(text("10", size, false), text(String(+m[2]).replace("-", "\u2212"), 0.7 * size, false), size)]);
  }
  function lay(node, size, ctx) {
    const wrap = function (child, minPrec) { const b = lay(child, size, ctx); return prec(child) < minPrec ? fence(b, "(", ")", size) : b; };
    const op = function (s) { return row([empty(0.22 * size), text(s, size, false), empty(0.22 * size)]); };
    switch (node.k) {
      case "num": return number(node.v, size);
      case "name":
        if (node.v === "pi" && !(node.v in ctx.declared)) return text(GLYPH.pi, size, true);
        return nameBox(node.v, size);
      case "neg": return row([text("\u2212", size, false), wrap(node.a, 6)]);
      case "not": return row([text("\u00ac", size, false), wrap(node.a, 9)]);
      case "bin": {
        const o = node.op;
        if (o === "/") return frac(lay(node.a, 0.8 * size, ctx), lay(node.b, 0.8 * size, ctx), size);
        if (o === "^") return sup(wrap(node.a, 9), lay(node.b, 0.7 * size, ctx), size);
        if (o === "*") {
          // A negated left factor needs no brackets: -a*x = -(a*x).
          const a = node.a.k === "neg" ? lay(node.a, size, ctx) : wrap(node.a, 6), b = wrap(node.b, 6), numeric = node.b.k === "num";
          return row([a, numeric ? op("\u00b7") : empty(0.12 * size), b]);
        }
        return row([wrap(node.a, DF.TEX_PREC[o]), op(OPS[o] || o), wrap(node.b, o === "-" || o === "+" ? DF.TEX_PREC[o] + 0.5 : DF.TEX_PREC[o])]);
      }
      case "call": {
        const f = node.f, A = node.args, arg = function (i) { return lay(A[i], size, ctx); };
        const list = function () { const parts = []; A.forEach(function (x, i) { if (i) parts.push(text(", ", size, false)); parts.push(lay(x, size, ctx)); }); return row(parts); };
        switch (f) {
          case "sqrt": return radical(arg(0), size);
          case "cbrt": return radical(arg(0), size, text("3", 0.5 * size, false));
          case "abs": return fence(arg(0), "|", "|", size);
          case "floor": return fence(arg(0), "floor", "floor", size);
          case "exp": return sup(text("e", size, true), lay(A[0], 0.7 * size, ctx), size);
          case "pow": return sup(wrap(A[0], 9), lay(A[1], 0.7 * size, ctx), size);
          case "lag": return row([lay(A[0], size, ctx), fence(row([text("t", size, true), op("\u2212"), arg(1)]), "(", ")", size)]);
          case "hill": {
            const s8 = 0.8 * size, x = lay(A[0], s8, ctx), K = lay(A[1], s8, ctx), n1 = lay(A[2], 0.56 * size, ctx), n2 = lay(A[2], 0.56 * size, ctx), n3 = lay(A[2], 0.56 * size, ctx);
            return frac(sup(x, n1, s8), row([sup(K, n2, s8), text(" + ", s8, false), sup(lay(A[0], s8, ctx), n3, s8)]), size);
          }
          case "step": case "heaviside": return row([text("\u0398", size, false), fence(arg(0), "(", ")", size)]);
          case "mod": return row([wrap(A[0], 6), op("mod"), wrap(A[1], 9)]);
          case "ifelse": return row([arg(1), text(" if ", size, false), arg(0), text(", else ", size, false), arg(2)]);
          default: return row([text(FUNC_NAME[f] || f, size, false), empty(0.06 * size), fence(list(), "(", ")", size)]);
        }
      }
    }
    return empty(0);
  }

  /* Boxes for every statement of a system, in source order: differential
     equations with a dot, maps with the subscript n + 1, SDEs in differential
     form, and aux definitions. */
  function systemBoxes(src, size) {
    const RE = DF.STATEMENT, out = [], declared = {};
    const lines = String(src).split(/\r?\n/).map(function (r) { return r.replace(/#.*$/, "").trim(); }).filter(Boolean);
    lines.forEach(function (s) { const m = RE.param.exec(s); if (m) declared[m[1]] = 1; });
    const noiseOf = {};
    lines.forEach(function (s) { const m = RE.noise.exec(s); if (m) noiseOf[m[1]] = m[2]; });
    const ctx = { declared: declared };
    const eq = function () { return row([empty(0.28 * size), text("=", size, false), empty(0.28 * size)]); };
    lines.forEach(function (s) {
      let m;
      try {
        if ((m = RE.aux.exec(s))) out.push(row([nameBox(m[1], size), eq(), lay(DF.parseExpression(m[2]), size, ctx)]));
        else if ((m = RE.ode.exec(s))) {
          const v = m[1] || m[2], rhs = lay(DF.parseExpression(m[3]), size, ctx);
          if (noiseOf[v] !== undefined) {
            const d = function () { return text("d", size, false); };
            out.push(row([d(), nameBox(v, size), eq(), fence(rhs, "(", ")", size), empty(0.1 * size), d(), text("t", size, true),
              row([empty(0.22 * size), text("+", size, false), empty(0.22 * size)]), lay(DF.parseExpression(noiseOf[v]), size, ctx), empty(0.1 * size), d(),
              sub(text("W", size, true), nameBox(v, 0.7 * size), size)]));
          } else out.push(row([nameBox(v, size, { dot: true }), eq(), rhs]));
        } else if ((m = RE.map.exec(s))) {
          out.push(row([nameBox(m[1], size, { extraSub: function (z) { return text("n+1", z, false); } }), eq(), lay(DF.parseExpression(m[2]), size, ctx)]));
        }
      } catch (e) { /* a statement that does not parse is left out */ }
    });
    return out;
  }

  // --------------------------------------------------------- painters
  function paint(g, box, x, y, color) {
    g.save();
    g.fillStyle = color; g.strokeStyle = color; g.textBaseline = "alphabetic"; g.textAlign = "left"; g.lineCap = "round"; g.lineJoin = "round";
    box.items.forEach(function (it) {
      if (it.k === "t") { g.font = fontOf(it.size, it.it); g.fillText(it.s, x + it.x, y + it.y); }
      else if (it.k === "r") g.fillRect(x + it.x, y + it.y, it.w, it.h);
      else if (it.k === "c") { g.beginPath(); g.arc(x + it.x, y + it.y, it.r, 0, 2 * Math.PI); g.fill(); }
      else if (it.k === "p") {
        g.lineWidth = it.lw; g.beginPath();
        for (let i = 0; i < it.pts.length; i += 2) { const px = x + it.x + it.pts[i], py = y + it.y + it.pts[i + 1]; if (i) g.lineTo(px, py); else g.moveTo(px, py); }
        g.stroke();
      }
    });
    g.restore();
  }
  function esc(t) { return String(t).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function svg(box, x, y, color) {
    const f = function (v) { return (+v).toFixed(2); };
    let s = "";
    box.items.forEach(function (it) {
      if (it.k === "t") s += '<text x="' + f(x + it.x) + '" y="' + f(y + it.y) + '" font-size="' + f(it.size) + '" font-family="' + esc(it.it ? MATH_FONT : MAIN_FONT) + '"' + (it.it ? ' font-style="italic"' : "") + ' fill="' + color + '">' + esc(it.s) + "</text>";
      else if (it.k === "r") s += '<rect x="' + f(x + it.x) + '" y="' + f(y + it.y) + '" width="' + f(it.w) + '" height="' + f(it.h) + '" fill="' + color + '"/>';
      else if (it.k === "c") s += '<circle cx="' + f(x + it.x) + '" cy="' + f(y + it.y) + '" r="' + f(it.r) + '" fill="' + color + '"/>';
      else if (it.k === "p") {
        let d = "";
        for (let i = 0; i < it.pts.length; i += 2) d += (i ? "L" : "M") + f(x + it.x + it.pts[i]) + " " + f(y + it.y + it.pts[i + 1]);
        s += '<path d="' + d + '" fill="none" stroke="' + color + '" stroke-width="' + f(it.lw) + '" stroke-linecap="round" stroke-linejoin="round"/>';
      }
    });
    return s;
  }

  /* A block of equations, right-aligned at x = right, first line at y = top.
     Returns the height used. Painted on a canvas (g given) or as SVG. */
  function block(src, size, lineGap) {
    const boxes = systemBoxes(src, size), gap = lineGap === undefined ? 0.35 * size : lineGap;
    let h = 0;
    const rows = boxes.map(function (b) { const y = h + b.a; h += b.a + b.d + gap; return { box: b, y: y }; });
    return { rows: rows, height: Math.max(0, h - gap), width: boxes.reduce(function (m, b) { return Math.max(m, b.w); }, 0) };
  }
  function paintBlock(g, src, right, top, size, color) {
    const B = block(src, size);
    B.rows.forEach(function (r) { paint(g, r.box, right - r.box.w, top + r.y, color); });
    return B;
  }
  function svgBlock(src, right, top, size, color) {
    const B = block(src, size);
    let s = "";
    B.rows.forEach(function (r) { s += svg(r.box, right - r.box.w, top + r.y, color); });
    return { svg: s, height: B.height, width: B.width };
  }

  DF.MathType = { systemBoxes: systemBoxes, block: block, paint: paint, svg: svg, paintBlock: paintBlock, svgBlock: svgBlock, nameParts: nameParts };
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/render/views.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Views. Each view draws one kind of figure from a running Simulator on the
   three layers of a Player: base (static, redrawn on demand), trail
   (accumulating, faded each frame) and top (cleared each frame).

     flow        particle ensemble with fading trails, 2D or rotating 3D
     trajectory  a few long orbits with a gradient tail, 2D or rotating 3D
     timeseries  scrolling time series of chosen variables
     phase       vector field, nullclines, classified equilibria, orbits
     sweep       slow parameter sweep over the continued equilibrium branches
     orbit       bifurcation diagram, built column by column
     density     ensemble density: 2D heat map, or a time carpet in 1D
     strobe      stroboscopic samples every period T, or a Poincare section
     cobweb      cobweb diagram of a one-dimensional map

   Every view implements init(P) and frame(P), and optionally drawStatic(P),
   pointer(P, kind, x, y, ev), onParam(P, i), svg(P) and legend(P). A frame
   advances the simulator by P.spf steps, the number the Player's clock
   holds for that frame; P.frameMs is its length in milliseconds. */
(function (DF) {
  "use strict";

  const V = {};
  const FRAME_MS = 1000 / 60;
  const PARAM_FORCING = { periodic: 1, quasiperiodic: 1, ramp: 1, step: 1 };

  // ------------------------------------------------------------ helpers
  // Length of the current frame in frames at 60 per second.
  function frames(P) { return (P.frameMs === undefined ? FRAME_MS : P.frameMs) / FRAME_MS; }
  function warm(P, fr) { const n = P.nominalSteps(fr); for (let s = 0; s < n; s++) P.sim.step(); }
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
        // 64 levels: a continuous gradation to the eye, few stroke batches.
        const u = Math.round(Math.min(1, Math.max(0, extra)) * 63) / 63;
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
    if (s.colorBy === "speed") {
      const v = speedOf(P, x);
      // On a sampled attractor, the logarithm of the speed between its 5th
      // and 95th percentiles spans the ramp from 0.12 to 0.82: slow
      // stretches, where particles crowd, keep a visible colour, and the
      // brightest stops, which additive blending turns white, are left out.
      // Elsewhere the band is a running envelope of the log speed of the
      // ensemble, at most five e-folds wide, that relaxes slowly inward.
      const l = Math.log(Math.max(v, 1e-300));
      if (!P.speedLog) {
        const E = P.speedEnv || (P.speedEnv = [l - 1e-6, l]);
        E[1] = Math.max(l, E[1] - 1e-6); E[0] = Math.max(E[1] - 5, Math.min(l, E[0] + 1e-6));
        return 0.12 + 0.7 * Math.min(1, Math.max(0, (l - E[0]) / Math.max(1e-9, E[1] - E[0])));
      }
      const L = P.speedLog;
      return 0.12 + 0.7 * Math.min(1, Math.max(0, (l - L[0]) / (L[1] - L[0])));
    }
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
    if (mode === "attractor" && P.sample) {
      // A point of the attractor sample, moved by 0.1 percent of each width.
      const S = P.sample, k = Math.floor(r.uniform() * S.n);
      for (let i = 0; i < out.length; i++) { const w = P.fullRanges[i][1] - P.fullRanges[i][0]; out[i] = S.X[k * S.dim + i] + 0.001 * w * r.normal(); }
      return out;
    }
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
  const f1 = function (v) { return (+v).toFixed(1); };
  function pathD(pts) {
    let d = "";
    for (let i = 0; i + 1 < pts.length; i += 2) d += (i ? "L" : "M") + f1(pts[i]) + " " + f1(pts[i + 1]);
    return d;
  }
  function polyline(pts, color, width, alpha, dash) {
    if (pts.length < 4) return "";
    return '<path d="' + pathD(pts) + '" fill="none" stroke="' + color + '" stroke-width="' + width + '" stroke-opacity="' + (+alpha).toFixed(3) + '" stroke-linejoin="round" stroke-linecap="round"' + (dash ? ' stroke-dasharray="' + dash + '"' : "") + "/>";
  }
  function clipOpen(id, f) { return '<clipPath id="' + id + '"><rect x="' + f1(f.L) + '" y="' + f1(f.T) + '" width="' + f1(f.R - f.L) + '" height="' + f1(f.B - f.T) + '"/></clipPath><g clip-path="url(#' + id + ')">'; }
  // A fitted view rocks about its broadest side with a period of `rockPeriod`
  // seconds (24 by default); any other 3D view turns at `rotate`.
  function rotate3D(P) {
    if (!P.cam.is3D()) return;
    const A = P.rockAmplitude ? P.rockAmplitude() : 0;
    if (A > 0) {
      const w = 2 * Math.PI / (P.scene.view.rockPeriod || 24);
      if (P.dragging) { P.rockBase = P.cam.azim - A * Math.sin(P.rockPhase); return; }
      P.rockPhase += w * frames(P) / 60;
      P.cam.azim = P.rockBase + A * Math.sin(P.rockPhase);
      return;
    }
    if (!P.dragging) P.cam.azim += (P.scene.view.rotate || 0) * 0.01 * frames(P);
  }
  // Tail of a trajectory, in stored points, and the stride between stored
  // steps: `tail` counts steps; without it the tail holds `seconds` of playback.
  function tailSpec(P, seconds, perSecond, maxPts) {
    const v = P.scene.view, sps = P.rate / P.sim.h, every = v.sampleEvery || Math.max(1, Math.ceil(sps / perSecond));
    const len = v.tail ? Math.ceil(v.tail / every) : Math.round((v.tailSeconds || seconds) * sps / every);
    return { every: every, len: Math.max(20, Math.min(maxPts, len)) };
  }
  function box3DSVG(P) {
    const cam = P.cam, th = P.theme, c = [0, 1], out = [];
    const corner = function (i, j, k) { const x = []; x[cam.axes[0]] = cam.ranges[0][i]; x[cam.axes[1]] = cam.ranges[1][j]; x[cam.axes[2]] = cam.ranges[2][k]; return cam.project(x, [0, 0]); };
    let d = "";
    c.forEach(function (i) { c.forEach(function (j) { c.forEach(function (k) { out.push([i, j, k]); }); }); });
    out.forEach(function (a) {
      out.forEach(function (b) {
        const dd = Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
        if (dd !== 1 || a.join() > b.join()) return;
        const u = corner(a[0], a[1], a[2]), w = corner(b[0], b[1], b[2]);
        d += "M" + f1(u[0]) + " " + f1(u[1]) + "L" + f1(w[0]) + " " + f1(w[1]);
      });
    });
    return '<path d="' + d + '" fill="none" stroke="' + th.grid + '" stroke-width="1"/>';
  }

  // ------------------------------------------------------------ flow
  V.flow = {
    label: "Flow (particles)",
    init: function (P) {
      const n = P.sim.n, st = new Float64Array(P.sys.vars.length);
      this.age = new Float32Array(n); this.life = new Float32Array(n);
      this.prev = new Float32Array(2 * n); this.has = new Uint8Array(n);
      this.buck = new Buckets(); this.pt = [0, 0]; this.x2 = new Float64Array(P.sys.vars.length);
      for (let k = 0; k < n; k++) {
        spawnState(P, st); P.sim.setMember(k, st);
        this.life[k] = this.newLife(P); this.age[k] = P.sim.rng.uniform() * this.life[k];
      }
      // Warm-up so that the first frame already shows the flow.
      warm(P, P.scene.view.warmup === undefined ? 60 : P.scene.view.warmup);
    },
    // Lifetimes are counted in frames at 60 frames per second.
    newLife: function (P) {
      const L = P.scene.view.life;
      if (!L || L[1] === Infinity || L === "inf") return Infinity;
      return L[0] + P.sim.rng.uniform() * (L[1] - L[0]);
    },
    frame: function (P) {
      const sim = P.sim, n = sim.n, cam = P.cam, st = P.scene.style, df = frames(P);
      rotate3D(P);
      P.fade(st.fade);
      const x = P.tmpX, discrete = P.sys.time === "discrete", st2 = this.x2;
      // A flow draws up to ten segments per particle and frame, one after
      // each share of the frame's steps, so that a fast particle traces its
      // curve and not a chord across it; a map draws its iterates once.
      const passes = discrete || P.spf < 1 ? 1 : Math.min(P.spf, 10);
      let done = 0;
      for (let pass = 0; pass < passes; pass++) {
        const upto = Math.round(P.spf * (pass + 1) / passes);
        for (; done < upto; done++) sim.step();
        for (let k = 0; k < n; k++) {
          if (pass === 0) {
            this.age[k] += df;
            if (!sim.alive[k] || this.age[k] > this.life[k]) {
              spawnState(P, st2); sim.setMember(k, st2); this.age[k] = 0; this.life[k] = this.newLife(P); this.has[k] = 0;
              continue;
            }
          } else if (!sim.alive[k]) continue;
          sim.member(k, x);
          if (outOfView(P, x)) { this.age[k] = this.life[k] + 1; continue; }
          const q = cam.project(x, this.pt);
          const col = colorFor(P, k, x, extraFor(P, k, x, this.age[k] / (isFinite(this.life[k]) ? this.life[k] : 1)));
          if (discrete) { if (P.spf > 0 || !this.has[k]) this.buck.seg(col, q[0], q[1], q[0], q[1]); }
          else if (this.has[k]) this.buck.seg(col, this.prev[2 * k], this.prev[2 * k + 1], q[0], q[1]);
          this.prev[2 * k] = q[0]; this.prev[2 * k + 1] = q[1]; this.has[k] = 1;
        }
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
      const n = P.sim.n, dim = P.sys.vars.length, spec = tailSpec(P, 12, 400, 12000);
      this.L = spec.len; this.every = spec.every; this.k = 0;
      this.buf = new Float64Array(n * this.L * dim); this.len = new Int32Array(n); this.head = new Int32Array(n);
      warm(P, P.scene.view.warmup || 0);
      this.push(P);
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
    colorOf: function (P, k, u) {
      const st = P.scene.style;
      if (st.colorBy === "time" || st.colorBy === "age") { const c = DF.rampRGB(st.ramp, u); return "rgb(" + c.join(",") + ")"; }
      return P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length];
    },
    frame: function (P) {
      rotate3D(P);
      for (let s = 0; s < P.spf; s++) { P.sim.step(); if (++this.k % this.every === 0) this.push(P); }
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
          const u = (g + 1) / G, col = this.colorOf(P, k, u);
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
      const st = P.scene.style, G = 28, discrete = P.sys.time === "discrete";
      if (P.cam.is3D() && (P.scene.view.box || P.scene.view.showAxes)) out += box3DSVG(P);
      else if (!P.cam.is3D() && P.scene.view.showAxes) out += DF.axesSVG(P.cam, st.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]]);
      for (let k = 0; k < P.sim.n; k++) {
        const pts = this.points(P, k), m = pts.length / 2;
        if (m < 2) continue;
        for (let g = 0; g < G; g++) {
          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
          if (b <= a) continue;
          const u = (g + 1) / G, col = this.colorOf(P, k, u), alpha = st.alpha * (0.08 + 0.92 * Math.pow(u, 1.4));
          if (discrete) { for (let j = a; j <= b; j++) out += '<rect x="' + f1(pts[2 * j] - st.pointSize / 2) + '" y="' + f1(pts[2 * j + 1] - st.pointSize / 2) + '" width="' + st.pointSize + '" height="' + st.pointSize + '" fill="' + col + '" fill-opacity="' + alpha.toFixed(3) + '"/>'; continue; }
          out += polyline(Array.from(pts.subarray(2 * a, 2 * b + 2)), col, st.lineWidth, alpha);
        }
        out += '<circle cx="' + f1(pts[pts.length - 2]) + '" cy="' + f1(pts[pts.length - 1]) + '" r="' + Math.max(2.5, st.lineWidth * 2) + '" fill="' + P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length] + '"/>';
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
      this.span = P.scene.view.window || (P.sys.time === "discrete" ? 100 : 50);
      // The buffer holds 4000 samples; one sample every `every` steps keeps
      // 10 percent more than a whole window in it.
      this.cap = 4000; this.every = Math.max(1, Math.ceil(this.span * 1.1 / (P.sim.h * this.cap)));
      this.T = new Float64Array(this.cap); this.Y = new Float64Array(this.cap * this.members * this.vi.length);
      this.len = 0; this.head = 0;
      let lo = Infinity, hi = -Infinity;
      this.vi.forEach(function (i) { lo = Math.min(lo, P.fullRanges[i][0]); hi = Math.max(hi, P.fullRanges[i][1]); });
      this.yr = P.scene.view.yRange || [lo, hi];
      this.record(P);
    },
    record: function (P) {
      const nv = this.vi.length, M = this.members, h = this.head;
      this.T[h] = P.sim.t;
      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) this.Y[(h * M + m) * nv + j] = P.sim.X[m * P.sys.vars.length + this.vi[j]];
      this.head = (h + 1) % this.cap; this.len = Math.min(this.len + 1, this.cap);
    },
    camera: function (P) {
      const t0 = Math.max(0, P.sim.t - this.span), cam = new DF.Camera([0, 1], [[t0, t0 + this.span], this.yr], { pad: 0.02 });
      cam.resize(P.w, P.h, P.inset);
      return cam;
    },
    labels: function (P) { return [P.sys.time === "discrete" ? "n" : "t", this.vi.map(function (i) { return P.sys.vars[i]; }).join(", ")]; },
    frame: function (P) {
      const nv = this.vi.length, M = this.members;
      for (let s = 0; s < P.spf; s++) { P.sim.step(); if (P.sim.steps % this.every === 0) this.record(P); }
      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);
      const cam = this.camera(P), t0 = cam.ranges[0][0];
      this.frameBox = DF.drawAxes(ctx, cam, P.scene.style.theme, this.labels(P), { grid: true });
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
      const nv = this.vi.length, M = this.members, cam = this.cam || this.camera(P), t0 = cam.ranges[0][0], st = P.scene.style;
      let out = DF.axesSVG(cam, st.theme, this.labels(P), { grid: true });
      const b = cam.plotBox();
      out += clipOpen("ts-clip", { L: b.x, T: b.y, R: b.x + b.w, B: b.y + b.h });
      const start = (this.head - this.len + this.cap) % this.cap, q = [0, 0], x = [0, 0];
      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) {
        const pts = [];
        for (let r = 0; r < this.len; r++) { const h = (start + r) % this.cap; if (this.T[h] < t0) continue; x[0] = this.T[h]; x[1] = this.Y[(h * M + m) * nv + j]; cam.project(x, q); pts.push(q[0], q[1]); }
        out += polyline(pts, P.palette[j % P.palette.length], st.lineWidth, M > 1 ? Math.max(0.25, st.alpha * 0.6) : st.alpha);
      }
      return out + "</g>";
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
  /* Join the segments of a contour into polylines, so that a dash pattern
     runs along the whole curve. Two segment ends meet when they agree to
     1e-6 of a grid cell (qx, qy); a cell edge crossed by the curve gives the
     same point to both cells up to rounding. */
  function stitch(segs, qx, qy) {
    const key = function (p) { return Math.round(p[0] / qx) + "," + Math.round(p[1] / qy); };
    const n = segs.length / 2, used = new Uint8Array(n), ends = new Map(), lines = [];
    for (let i = 0; i < 2 * n; i++) { const k = key(segs[i]); let a = ends.get(k); if (!a) { a = []; ends.set(k, a); } a.push(i); }
    const follow = function (tip) {
      const out = [];
      for (;;) {
        const c = (ends.get(key(tip)) || []).find(function (e) { return !used[e >> 1]; });
        if (c === undefined) return out;
        used[c >> 1] = 1;
        tip = segs[c ^ 1];
        out.push(tip);
      }
    };
    for (let s = 0; s < n; s++) {
      if (used[s]) continue;
      used[s] = 1;
      const fwd = follow(segs[2 * s + 1]), back = follow(segs[2 * s]);
      lines.push(back.reverse().concat([segs[2 * s], segs[2 * s + 1]], fwd));
    }
    return lines;
  }
  function equilibriumStyle(P, e) {
    const th = P.theme;
    if (e.stable) return { fill: th.ink, stroke: th.ink };
    if (/saddle/.test(e.type)) return { fill: P.palette[0], stroke: th.ink };
    return { fill: th.dark ? "#0b0620" : "#ffffff", stroke: th.ink };
  }

  V.phase = {
    label: "Phase plane",
    axes: true,
    init: function (P) {
      this.trails = []; this.maxTrails = 40;
      this.spec = tailSpec(P, 10, 300, 6000);
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
    /* Every orbit of the plane is a one-member simulator started now, so it
       keeps its own time (forced systems) and its own delay history (delay
       equations); it follows the deterministic skeleton, without noise, and
       the parameter forcing of the scene. */
    addTrail: function (P, s) {
      if (this.trails.length >= this.maxTrails) this.trails.shift();
      const sim = new DF.Simulator(P.sys, {
        n: 1, dt: P.sim.h, params: P.sim.base, init: s, t0: P.sim.t, deterministic: true, keepPositive: !!P.scene.keepPositive,
        perturbations: P.scene.perturbations.filter(function (q) { return q.kind in PARAM_FORCING; })
      });
      this.trails.push({ sim: sim, pts: [], k: 0 });
    },
    drawStatic: function (P) {
      const ctx = P.ctx.base, cam = P.cam, sys = P.sys, v = P.scene.view, st = P.scene.style;
      const frame = DF.drawAxes(ctx, cam, st.theme, [sys.vars[cam.axes[0]], sys.vars[cam.axes[1]]], { grid: false });
      this.frameBox = frame;
      const a0 = cam.axes[0], a1 = cam.axes[1], t = P.sim.t, p = P.sim.p;
      const x = Float64Array.from(P.sim.init), d = new Float64Array(sys.vars.length);
      const H = function (i) { return x[i]; };
      const R0 = cam.ranges[0], R1 = cam.ranges[1];
      ctx.save(); ctx.beginPath(); ctx.rect(frame.L, frame.T, frame.R - frame.L, frame.B - frame.T); ctx.clip();
      // Vector field: arrows on a grid, length by log speed.
      this.arrows = [];
      if (v.field !== "none") {
        const nx = v.fieldDensity || 22, ny = Math.max(2, Math.round(nx * (frame.B - frame.T) / (frame.R - frame.L)));
        const cellW = (frame.R - frame.L) / nx, cellH = (frame.B - frame.T) / ny;
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
        const self = this;
        vals.forEach(function (e) {
          if (!(e[4] > 0) || !isFinite(e[4])) return;
          const rel = Math.log1p(9 * e[4] / vmax) / Math.log(10);
          const len = Math.min(cellW, cellH) * 0.42 * (0.35 + 0.65 * rel), ux = e[2] / e[4], uy = e[3] / e[4];
          const c = DF.rampRGB(st.ramp, 0.25 + 0.75 * rel);
          const col = "rgba(" + c.join(",") + "," + (P.theme.dark ? 0.55 : 0.75) + ")";
          const x0 = e[0] - ux * len, y0 = e[1] - uy * len, x1 = e[0] + ux * len, y1 = e[1] + uy * len;
          const hd = [x1 - 4 * ux + 2.5 * uy, y1 - 4 * uy - 2.5 * ux, x1 - 4 * ux - 2.5 * uy, y1 - 4 * uy + 2.5 * ux];
          ctx.strokeStyle = col;
          ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1);
          ctx.moveTo(x1, y1); ctx.lineTo(hd[0], hd[1]);
          ctx.moveTo(x1, y1); ctx.lineTo(hd[2], hd[3]);
          ctx.stroke();
          self.arrows.push({ c: col, d: "M" + f1(x0) + " " + f1(y0) + "L" + f1(x1) + " " + f1(y1) + "M" + f1(x1) + " " + f1(y1) + "L" + f1(hd[0]) + " " + f1(hd[1]) + "M" + f1(x1) + " " + f1(y1) + "L" + f1(hd[2]) + " " + f1(hd[3]) });
        });
      }
      // Nullclines f_a0 = 0 and f_a1 = 0 by marching squares, joined into polylines.
      this.nullLines = [];
      if (v.nullclines !== false && P.sys.time === "continuous") {
        const nx = 160, ny = 160, g0 = new Float64Array(nx * ny), g1 = new Float64Array(nx * ny);
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
          x[a0] = R0[0] + i / (nx - 1) * (R0[1] - R0[0]); x[a1] = R1[0] + j / (ny - 1) * (R1[1] - R1[0]);
          sys.f(t, x, p, d, H);
          g0[j * nx + i] = d[a0]; g1[j * nx + i] = d[a1];
        }
        const X = function (i) { return R0[0] + i / (nx - 1) * (R0[1] - R0[0]); }, Y = function (j) { return R1[0] + j / (ny - 1) * (R1[1] - R1[0]); };
        const qx = 1e-6 * (R0[1] - R0[0]) / (nx - 1), qy = 1e-6 * (R1[1] - R1[0]) / (ny - 1);
        const self = this;
        [g0, g1].forEach(function (g, idx) {
          const lines = stitch(contour(g, nx, ny, X, Y), qx, qy);
          const col = P.palette[(idx + 1) % P.palette.length], dash = idx ? [6, 4] : [];
          ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.setLineDash(dash); ctx.globalAlpha = 0.9; ctx.lineJoin = "round";
          const q = [0, 0], z = Float64Array.from(P.sim.init);
          lines.forEach(function (line) {
            const pts = [];
            line.forEach(function (pt) { z[a0] = pt[0]; z[a1] = pt[1]; cam.project(z, q); pts.push(q[0], q[1]); });
            ctx.beginPath(); ctx.moveTo(pts[0], pts[1]);
            for (let i = 2; i < pts.length; i += 2) ctx.lineTo(pts[i], pts[i + 1]);
            ctx.stroke();
            self.nullLines.push({ pts: pts, color: col, dash: idx ? "6 4" : "" });
          });
          ctx.setLineDash([]); ctx.globalAlpha = 1;
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
          const sty = equilibriumStyle(P, e);
          ctx.lineWidth = 2; ctx.strokeStyle = sty.stroke; ctx.fillStyle = sty.fill;
          ctx.beginPath(); ctx.arc(q[0], q[1], 5.5, 0, 7); ctx.fill(); if (!e.stable) ctx.stroke();
        });
      }
    },
    frame: function (P) {
      const sys = P.sys, cam = P.cam, st = P.scene.style, a0 = cam.axes[0], a1 = cam.axes[1];
      if (sys.usesTime || P.sim.perturbations.some(function (q) { return q.enabled !== false && DF.PARAM_PERTURBATIONS.indexOf(q.kind) >= 0; })) {
        if ((P.frameCount % 20) === 0) P.redrawStatic();
      }
      for (let s = 0; s < P.spf; s++) P.sim.step();
      const every = this.spec.every, maxPts = 2 * this.spec.len;
      this.trails.forEach(function (tr) {
        if (tr.dead) return;
        tr.sim.base.set(P.sim.base); tr.sim.updateParams();
        const x = tr.sim.X;
        if (!tr.pts.length) tr.pts.push(x[a0], x[a1]);
        for (let s = 0; s < P.spf; s++) {
          tr.sim.step();
          if (!tr.sim.alive[0] || outOfView(P, x)) { tr.dead = true; break; }
          if (++tr.k % every === 0) tr.pts.push(x[a0], x[a1]);
        }
        if (tr.pts.length > maxPts) tr.pts.splice(0, tr.pts.length - maxPts);
      });
      this.trails = this.trails.filter(function (tr) { return !tr.dead || tr.pts.length > 2; });
      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);
      ctx.save();
      if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }
      const self = this;
      this.trails.forEach(function (tr, k) {
        const col = P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length];
        const pts = self.screen(P, tr), m = pts.length / 2;
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
    },
    // Screen points of a trail, which is stored in state coordinates.
    screen: function (P, tr) {
      const cam = P.cam, z = Float64Array.from(P.sim.init), q = [0, 0], a0 = cam.axes[0], a1 = cam.axes[1], out = new Array(tr.pts.length);
      for (let i = 0; i < tr.pts.length; i += 2) { z[a0] = tr.pts[i]; z[a1] = tr.pts[i + 1]; cam.project(z, q); out[i] = q[0]; out[i + 1] = q[1]; }
      return out;
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down") return false;
      const u = P.cam.unproject(px, py), s = Float64Array.from(P.sim.init);
      s[P.cam.axes[0]] = u[0]; s[P.cam.axes[1]] = u[1];
      this.addTrail(P, s);
      return true;
    },
    onParam: function (P) { P.markStatic(); },
    legend: function (P) {
      const L = [["stable", P.theme.ink, "dot"], ["saddle", P.palette[0], "dot"], ["unstable", P.theme.ink, "ring"]];
      if (P.scene.view.nullclines !== false) L.unshift([P.sys.vars[P.cam.axes[0]] + "-nullcline", P.palette[1 % P.palette.length], "line"], [P.sys.vars[P.cam.axes[1]] + "-nullcline", P.palette[2 % P.palette.length], "dash"]);
      return L;
    },
    svg: function (P) {
      const cam = P.cam, st = P.scene.style, f = this.frameBox, self = this;
      let out = DF.axesSVG(cam, st.theme, [P.sys.vars[cam.axes[0]], P.sys.vars[cam.axes[1]]]);
      if (!f) return out;
      out += clipOpen("phase-clip", f);
      (this.arrows || []).forEach(function (a) { out += '<path d="' + a.d + '" fill="none" stroke="' + a.c + '" stroke-width="1"/>'; });
      (this.nullLines || []).forEach(function (nc) { out += polyline(nc.pts, nc.color, 1.6, 0.9, nc.dash); });
      this.trails.forEach(function (tr, k) {
        const col = P.palette[(st.colorBy === "member" ? k : 0) % P.palette.length], pts = self.screen(P, tr), m = pts.length / 2, G = 10;
        for (let g = 0; g < G; g++) {
          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
          if (b > a) out += polyline(pts.slice(2 * a, 2 * b + 2), col, st.lineWidth, st.alpha * (0.15 + 0.85 * (g + 1) / G));
        }
        if (!tr.dead && m >= 1) out += '<circle cx="' + f1(pts[pts.length - 2]) + '" cy="' + f1(pts[pts.length - 1]) + '" r="3.5" fill="' + col + '"/>';
      });
      out += "</g>";
      const q = [0, 0], z = Float64Array.from(P.sim.init);
      (this.equilibria || []).forEach(function (e) {
        z[cam.axes[0]] = e.x[0]; z[cam.axes[1]] = e.x[1]; cam.project(z, q);
        const sty = equilibriumStyle(P, e);
        out += '<circle cx="' + f1(q[0]) + '" cy="' + f1(q[1]) + '" r="5.5" fill="' + sty.fill + '"' + (e.stable ? "" : ' stroke="' + sty.stroke + '" stroke-width="2"') + "/>";
      });
      return out;
    }
  };

  // ------------------------------------------------------------ sweep
  /* Fallback when continuation finds no branch: equilibria found by Newton
     iteration at each of `cols` parameter values, drawn as points. */
  function scanBranches(P, pi, vi, from, to, cols) {
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
  const POINT_LABEL = { fold: "fold", hopf: "Hopf", branch: "branch point", flip: "period doubling", torus: "torus" };

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
      this.data = null;
      P.sim.base[this.pi] = this.pval; P.sim.updateParams();
    },
    // Branches as runs of equal stability, in data coordinates [p, x, p, x, ...].
    compute: function (P) {
      const sys = P.sys, pi = this.pi, vi = this.vi;
      let res = null;
      try { res = DF.continueBranches(sys, Float64Array.from(P.sim.base), pi, this.from, this.to, P.fullRanges.map(function (r) { return r.slice(); })); } catch (e) { res = null; }
      const runs = [];
      if (res) res.branches.forEach(function (br) {
        let cur = null;
        br.forEach(function (q) {
          if (!cur || cur.stable !== q.stable) {
            const last = cur ? cur.pts.slice(-2) : null;
            cur = { stable: q.stable, pts: last ? last : [] };
            runs.push(cur);
          }
          cur.pts.push(q.p, q.x[vi]);
        });
      });
      if (!runs.length) return { runs: [], points: [], dots: scanBranches(P, pi, vi, this.from, this.to, 200) };
      return { runs: runs, points: res.points.map(function (q) { return { kind: q.kind, p: q.p, v: q.x[vi] }; }), dots: null };
    },
    drawStatic: function (P) {
      this.cam.resize(P.w, P.h, P.inset);
      const ctx = P.ctx.base, cam = this.cam, th = P.theme, pname = P.sys.params[this.pi].name;
      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [pname, P.sys.vars[this.vi]], { grid: false });
      this.frameBox = f;
      if (P.scene.view.branches === false) return;
      if (!this.data) this.data = this.compute(P);
      const q = [0, 0], x = [0, 0];
      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      ctx.lineJoin = "round"; ctx.lineCap = "round";
      this.data.runs.forEach(function (run) {
        ctx.strokeStyle = run.stable ? th.ink : th.muted; ctx.lineWidth = run.stable ? 2.2 : 1.6; ctx.setLineDash(run.stable ? [] : [6, 5]);
        ctx.beginPath();
        for (let i = 0; i < run.pts.length; i += 2) { x[0] = run.pts[i]; x[1] = run.pts[i + 1]; cam.project(x, q); if (i) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]); }
        ctx.stroke();
      });
      ctx.setLineDash([]);
      (this.data.dots || []).forEach(function (b) {
        x[0] = b.p; x[1] = b.v; cam.project(x, q);
        ctx.fillStyle = b.stable ? th.ink : th.muted; ctx.beginPath(); ctx.arc(q[0], q[1], b.stable ? 1.7 : 1.1, 0, 7); ctx.fill();
      });
      ctx.restore();
      ctx.font = "10.5px Jost, system-ui, sans-serif"; ctx.textBaseline = "middle";
      this.data.points.forEach(function (sp) {
        x[0] = sp.p; x[1] = sp.v; cam.project(x, q);
        if (q[0] < f.L - 1 || q[0] > f.R + 1 || q[1] < f.T - 1 || q[1] > f.B + 1) return;
        ctx.lineWidth = 1.6; ctx.strokeStyle = th.ink; ctx.fillStyle = th.dark ? "#0b0620" : "#ffffff";
        ctx.beginPath(); ctx.arc(q[0], q[1], 4.5, 0, 7); ctx.fill(); ctx.stroke();
        const right = q[0] < (f.L + f.R) / 2;
        ctx.fillStyle = th.muted; ctx.textAlign = right ? "left" : "right";
        ctx.fillText(POINT_LABEL[sp.kind] + ", " + pname + " = " + DF.fmt(sp.p), q[0] + (right ? 9 : -9), q[1]);
      });
    },
    frame: function (P) {
      const sim = P.sim, span = this.to - this.from, h = sim.h;
      if (!P.userParam) {
        // The parameter moves by speed x span per unit of model time, a little at every step.
        const dp = P.sweepSpeed * span * (P.sys.time === "discrete" ? 1 : h);
        for (let s = 0; s < P.spf; s++) {
          this.pval += this.dir * dp;
          if (this.pval > this.to) { this.pval = this.to; this.dir = -1; }
          if (this.pval < this.from) { this.pval = this.from; this.dir = 1; }
          sim.base[this.pi] = this.pval; sim.updateParams();
          sim.step();
        }
        if (P.spf) P.emit("param", { name: P.sys.params[this.pi].name, value: this.pval });
      } else {
        this.pval = sim.base[this.pi];
        for (let s = 0; s < P.spf; s++) sim.step();
      }
      if (!sim.alive[0]) sim.setMember(0, sim.init);
      if (P.spf || !this.trail.length) {
        this.trail.push(this.pval, sim.X[this.vi]);
        if (this.trail.length > 2 * (P.scene.view.tail || 1400)) this.trail.splice(0, 2);
      }
      const ctx = P.ctx.top, st = P.scene.style, q = [0, 0], x = [0, 0], cam = this.cam;
      ctx.clearRect(0, 0, P.w, P.h);
      ctx.save(); if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }
      const pts = [];
      for (let i = 0; i < this.trail.length; i += 2) { x[0] = this.trail[i]; x[1] = this.trail[i + 1]; cam.project(x, q); pts.push(q[0], q[1]); }
      const m = pts.length / 2, G = 16;
      ctx.strokeStyle = P.palette[1 % P.palette.length]; ctx.lineWidth = st.lineWidth; ctx.lineJoin = "round";
      for (let g = 0; g < G; g++) {
        const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);
        if (b <= a) continue;
        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * (g + 1) / G);
        ctx.beginPath(); ctx.moveTo(pts[2 * a], pts[2 * a + 1]);
        for (let j = a + 1; j <= b; j++) ctx.lineTo(pts[2 * j], pts[2 * j + 1]);
        ctx.stroke();
      }
      x[0] = this.pval; x[1] = sim.X[this.vi]; cam.project(x, q);
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
    // Moving the swept parameter by hand stops the sweep; any other parameter changes the branches.
    onParam: function (P, i) {
      if (i === this.pi) { P.userParam = true; return; }
      this.data = null; P.markStatic();
    },
    legend: function (P) {
      const L = [["stable branch", P.theme.ink, "line"], ["unstable branch", P.theme.muted, "dash"]];
      if (this.data && this.data.points.length) L.push([Array.from(new Set(this.data.points.map(function (q) { return POINT_LABEL[q.kind]; }))).join(", "), P.theme.ink, "ring"]);
      L.push(["state under the sweep", P.palette[1 % P.palette.length], "line"]);
      return L;
    },
    svg: function (P) {
      const cam = this.cam, th = P.theme, st = P.scene.style, f = this.frameBox, pname = P.sys.params[this.pi].name, q = [0, 0], x = [0, 0];
      let out = DF.axesSVG(cam, st.theme, [pname, P.sys.vars[this.vi]]);
      if (!f) return out;
      const toScreen = function (arr) { const pts = []; for (let i = 0; i < arr.length; i += 2) { x[0] = arr[i]; x[1] = arr[i + 1]; cam.project(x, q); pts.push(q[0], q[1]); } return pts; };
      out += clipOpen("sweep-clip", f);
      if (this.data) {
        this.data.runs.forEach(function (run) { out += polyline(toScreen(run.pts), run.stable ? th.ink : th.muted, run.stable ? 2.2 : 1.6, 1, run.stable ? "" : "6 5"); });
        (this.data.dots || []).forEach(function (b) { const s = toScreen([b.p, b.v]); out += '<circle cx="' + f1(s[0]) + '" cy="' + f1(s[1]) + '" r="' + (b.stable ? 1.7 : 1.1) + '" fill="' + (b.stable ? th.ink : th.muted) + '"/>'; });
      }
      out += polyline(toScreen(this.trail), P.palette[1 % P.palette.length], st.lineWidth, st.alpha);
      out += "</g>";
      if (this.data) this.data.points.forEach(function (sp) {
        const s = toScreen([sp.p, sp.v]);
        if (s[0] < f.L - 1 || s[0] > f.R + 1 || s[1] < f.T - 1 || s[1] > f.B + 1) return;
        const right = s[0] < (f.L + f.R) / 2;
        out += '<circle cx="' + f1(s[0]) + '" cy="' + f1(s[1]) + '" r="4.5" fill="' + (th.dark ? "#0b0620" : "#ffffff") + '" stroke="' + th.ink + '" stroke-width="1.6"/>';
        out += '<text x="' + f1(s[0] + (right ? 9 : -9)) + '" y="' + f1(s[1] + 3.5) + '" text-anchor="' + (right ? "start" : "end") + '" font-family="Jost, sans-serif" font-size="10.5" fill="' + th.muted + '">' + POINT_LABEL[sp.kind] + ", " + pname + " = " + DF.fmt(sp.p) + "</text>";
      });
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
      // Delay equations need their history and random maps their draws, so
      // they run in a simulator; flows and maps without randomness use a
      // bare integrator, and stochastic equations their drift alone.
      this.useSim = sys.kind === "dde" || sys.usesRandom;
      this.runner = null;
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
      if (!this.cols || this.col > this.cols) { this.drawCursor(P); return; }
      const v = P.scene.view, sys = P.sys, h = P.sim.h, st = P.scene.style, discrete = sys.time === "discrete";
      const transient = v.transient || (discrete ? 300 : Math.round(200 / h));
      const samples = v.samples || (discrete ? 150 : Math.round(400 / h));
      const ctx = P.ctx.trail, q = [0, 0], budget = (typeof performance !== "undefined" ? performance.now() : Date.now()) + (v.budget || 12);
      const x = this.state, tmp = P.tmpDx, H = function (i) { return x[i]; };
      const clock = function () { return typeof performance !== "undefined" ? performance.now() : Date.now(); };
      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;
      while (this.col <= this.cols && clock() < budget) {
        const pv = this.from + (this.to - this.from) * this.col / this.cols;
        this.p[this.pi] = pv;
        let adv, read;
        // Following the attractor, a state that sits on an equilibrium to the
        // last bit stays there after the equilibrium turns unstable (with a
        // constant delay history nothing seeds the oscillation), so each new
        // column starts from the previous state nudged by one part in 1e6.
        const kick = function (z) { for (let i = 0; i < z.length; i++) z[i] += 1e-6 * (Math.abs(z[i]) + 1e-6) * (i % 2 ? -1 : 1); };
        if (this.useSim) {
          if (!this.runner || v.follow === false || !this.runner.alive[0]) this.runner = new DF.Simulator(sys, { n: 1, dt: h, params: this.p, init: P.sim.init, seed: P.scene.seed + this.col });
          else { this.runner.base[this.pi] = pv; this.runner.updateParams(); kick(this.runner.X); }
          const R = this.runner;
          adv = function () { R.step(); };
          read = function () { return R.alive[0] ? R.X[this.vi] : NaN; }.bind(this);
        } else {
          if (v.follow === false || !x.every(isFinite)) x.set(P.sim.init); else if (this.col) kick(x);
          let t = 0;
          const self = this;
          adv = function () { if (discrete) { sys.f(t, x, self.p, tmp, H); x.set(tmp); t += 1; } else { self.rk(sys.f, t, x, self.p, h, H); t += h; } };
          read = function () { return x[self.vi]; };
        }
        let prev2 = NaN, prev1 = NaN;
        for (let s = 0; s < transient; s++) adv();
        for (let s = 0; s < samples; s++) {
          adv();
          const y = read();
          if (!isFinite(y)) break;
          if (discrete) { this.points.push(pv, y); this.cam.project([pv, y], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2); }
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
      this.drawCursor(P);
    },
    drawCursor: function (P) {
      const top = P.ctx.top; top.clearRect(0, 0, P.w, P.h);
      if (this.cols && this.col <= this.cols && this.frameBox) {
        const px = this.frameBox.L + (this.frameBox.R - this.frameBox.L) * this.col / this.cols;
        top.strokeStyle = P.palette[1 % P.palette.length]; top.globalAlpha = 0.6; top.beginPath(); top.moveTo(px, this.frameBox.T); top.lineTo(px, this.frameBox.B); top.stroke(); top.globalAlpha = 1;
      }
    },
    onParam: function (P) { P.restart(); },
    legend: function (P) {
      if (P.sys.time === "discrete") return [["iterates after a transient", P.palette[0], "dot"]];
      return [[P.sys.kind === "sde" ? "local maxima of the drift (noise off) after a transient" : "local maxima after a transient", P.palette[0], "dot"]];
    }
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
      warm(P, v.warmup || 0);
    },
    drawStatic: function (P) {
      const st = P.scene.style;
      if (this.mode === "carpet") {
        this.span = P.scene.view.window || 60;
        this.cam = new DF.Camera([0, 1], [[-this.span, 0], P.fullRanges[this.vi]], { pad: 0.0 });
        this.cam.resize(P.w, P.h, P.inset);
        this.frameBox = DF.drawAxes(P.ctx.base, this.cam, st.theme, ["t - t now", P.sys.vars[this.vi]], {});
      } else {
        this.frameBox = DF.drawAxes(P.ctx.base, P.cam, st.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]], {});
      }
      const f = this.frameBox;
      this.nx = Math.max(10, Math.round((f.R - f.L) / (P.scene.view.cell || 3)));
      this.ny = Math.max(10, Math.round((f.B - f.T) / (P.scene.view.cell || 3)));
      this.grid = new Float32Array(this.nx * this.ny);
      // A carpet column covers window / nx units of time, so the axis spans the window.
      this.colDt = (this.span || 1) / this.nx; this.nextCol = P.sim.t + this.colDt;
      if (typeof document !== "undefined") {
        this.off = document.createElement("canvas"); this.off.width = this.nx; this.off.height = this.ny;
        this.octx = this.off.getContext("2d"); this.img = this.octx.createImageData(this.nx, this.ny);
      }
    },
    column: function (P) {
      const sim = P.sim, dim = P.sys.vars.length, nx = this.nx, ny = this.ny, g = this.grid, r = P.fullRanges[this.vi];
      for (let j = 0; j < ny; j++) { g.copyWithin(j * nx, j * nx + 1, j * nx + nx); g[j * nx + nx - 1] = 0; }
      for (let k = 0; k < sim.n; k++) {
        const y = sim.X[k * dim + this.vi], j = Math.floor((r[1] - y) / (r[1] - r[0]) * ny);
        if (j >= 0 && j < ny) g[j * nx + nx - 1] += 1;
      }
      let mx = 0; for (let j = 0; j < ny; j++) mx = Math.max(mx, g[j * nx + nx - 1]);
      for (let j = 0; j < ny; j++) g[j * nx + nx - 1] /= (mx || 1);
    },
    frame: function (P) {
      const sim = P.sim, dim = P.sys.vars.length, v = P.scene.view, st = P.scene.style;
      const respawn = function () { for (let k = 0; k < sim.n; k++) if (!sim.alive[k]) sim.setMember(k, spawnState(P, new Float64Array(dim))); };
      if (this.mode === "carpet" && this.grid) {
        for (let s = 0; s < P.spf; s++) {
          sim.step();
          for (let c = 0; sim.t >= this.nextCol && c < this.nx; c++) { this.column(P); this.nextCol += this.colDt; }
          if (sim.t >= this.nextCol) this.nextCol = sim.t + this.colDt;
        }
        respawn();
      } else {
        for (let s = 0; s < P.spf; s++) sim.step();
        respawn();
      }
      if (!this.grid || !this.img) return;
      const nx = this.nx, ny = this.ny, g = this.grid;
      if (this.mode !== "carpet") {
        const decay = Math.pow(v.decay === undefined ? 0.85 : v.decay, frames(P));
        for (let i = 0; i < g.length; i++) g[i] *= decay;
        const cam = P.cam, a0 = cam.axes[0], a1 = cam.axes[1], R0 = cam.ranges[0], R1 = cam.ranges[1];
        for (let k = 0; k < sim.n; k++) {
          const u = (sim.X[k * dim + a0] - R0[0]) / (R0[1] - R0[0]), w = (R1[1] - sim.X[k * dim + a1]) / (R1[1] - R1[0]);
          const i = Math.floor(u * nx), j = Math.floor(w * ny);
          if (i >= 0 && i < nx && j >= 0 && j < ny) g[j * nx + i] += frames(P);
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
      this.period = P.strobePeriod();
      // The Player chose h = T / N, so the state is sampled every N steps, at t = k T exactly.
      this.N = Math.max(1, Math.round(this.period / (P.sys.time === "discrete" ? 1 : P.sim.h)));
      this.offset = Math.round((((v.phase || 0) % 1) + 1) % 1 * this.N);
      this.prev = null;
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
      const self = this;
      const plot = function (x) { self.points.push(x[a0], x[a1]); if (self.points.length > 400000) self.points.splice(0, 2); P.cam.project(x, q); ctx.fillRect(q[0] - st.pointSize / 2, q[1] - st.pointSize / 2, st.pointSize, st.pointSize); };
      const x = P.tmpX;
      for (let s = 0; s < P.spf; s++) {
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
          if ((sim.steps - this.offset) % this.N === 0) {
            this.count++;
            if (this.count > this.transient) for (let k = 0; k < sim.n; k++) { if (sim.alive[k]) plot(sim.member(k, x)); }
          }
        }
      }
      ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";
    },
    // Points sampled under another parameter belong to another attractor.
    onParam: function (P) { this.points = []; this.count = 0; P.clearTrail(); },
    legend: function (P) { return [[P.scene.view.mode === "section" ? "upward crossings of the section" : "state at t = t0 + kT, T = " + DF.fmt(this.period), P.palette[0], "dot"]]; }
  };

  // ------------------------------------------------------------ cobweb
  V.cobweb = {
    label: "Cobweb (1D maps)",
    axes: true,
    init: function (P) {
      this.x = P.sim.init[0]; this.path = []; this.vi = 0; this.n = 0;
      const r = P.fullRanges[0];
      this.cam = new DF.Camera([0, 1], [r, r], { pad: 0.03 });
      this.cam.resize(P.w, P.h, P.inset);
    },
    curve: function (P) {
      const r = this.cam.ranges[0], out = new Float64Array(P.sys.vars.length), x = Float64Array.from(P.sim.init), pts = [];
      for (let i = 0; i <= 600; i++) { x[0] = r[0] + (r[1] - r[0]) * i / 600; P.sys.f(this.n, x, P.sim.p, out); pts.push(x[0], out[0]); }
      return pts;
    },
    drawStatic: function (P) {
      this.cam.resize(P.w, P.h, P.inset);
      const ctx = P.ctx.base, cam = this.cam, th = P.theme, v = P.sys.vars[0];
      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [v + "\u2099", v + "\u2099\u208a\u2081"], {});
      this.frameBox = f;
      const r = cam.ranges[0], q = [0, 0], c = this.curve(P);
      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      ctx.strokeStyle = th.muted; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
      cam.project([r[0], r[0]], q); ctx.beginPath(); ctx.moveTo(q[0], q[1]); cam.project([r[1], r[1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();
      ctx.setLineDash([]); ctx.strokeStyle = P.palette[2 % P.palette.length]; ctx.lineWidth = 2; ctx.beginPath();
      for (let i = 0; i < c.length; i += 2) { cam.project([c[i], c[i + 1]], q); if (i) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]); }
      ctx.stroke(); ctx.restore();
    },
    frame: function (P) {
      const out = P.tmpDx, x = P.tmpX;
      for (let s = 0; s < P.spf; s++) {
        x[0] = this.x; P.sys.f(this.n, x, P.sim.p, out);
        const y = out[0];
        if (!isFinite(y)) { this.x = P.sim.init[0]; this.path = []; continue; }
        if (!this.path.length) this.path.push(this.x, this.cam.ranges[1][0] < 0 && this.cam.ranges[1][1] > 0 ? 0 : this.cam.ranges[1][0]);
        this.path.push(this.x, y, y, y);
        if (this.path.length > 2 * (P.scene.view.tail || 120)) this.path.splice(0, 4);
        this.x = y; this.n += 1;
      }
      const ctx = P.ctx.top, st = P.scene.style, q = [0, 0];
      ctx.clearRect(0, 0, P.w, P.h);
      if (!this.frameBox) return;
      ctx.save(); const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();
      const m = this.path.length / 2;
      ctx.strokeStyle = P.palette[0]; ctx.lineWidth = st.lineWidth;
      for (let j = 1; j < m; j++) {
        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * j / m);
        ctx.beginPath(); this.cam.project([this.path[2 * j - 2], this.path[2 * j - 1]], q); ctx.moveTo(q[0], q[1]);
        this.cam.project([this.path[2 * j], this.path[2 * j + 1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();
      }
      if (m) { ctx.globalAlpha = 1; ctx.fillStyle = P.palette[0]; ctx.beginPath(); ctx.arc(q[0], q[1], 4, 0, 7); ctx.fill(); }
      ctx.restore();
    },
    pointer: function (P, kind, px, py) {
      if (kind !== "down") return false;
      this.x = this.cam.unproject(px, py)[0]; this.path = [];
      return true;
    },
    onParam: function (P) { P.markStatic(); },
    svg: function (P) {
      const cam = this.cam, th = P.theme, st = P.scene.style, f = this.frameBox, v = P.sys.vars[0], r = cam.ranges[0], q = [0, 0];
      let out = DF.axesSVG(cam, st.theme, [v + "\u2099", v + "\u2099\u208a\u2081"]);
      if (!f) return out;
      const toScreen = function (arr) { const pts = []; for (let i = 0; i < arr.length; i += 2) { cam.project([arr[i], arr[i + 1]], q); pts.push(q[0], q[1]); } return pts; };
      out += clipOpen("cobweb-clip", f);
      out += polyline(toScreen([r[0], r[0], r[1], r[1]]), th.muted, 1, 1, "4 4");
      out += polyline(toScreen(this.curve(P)), P.palette[2 % P.palette.length], 2, 1);
      const pts = toScreen(this.path), m = pts.length / 2;
      for (let j = 1; j < m; j++) out += polyline(pts.slice(2 * j - 2, 2 * j + 2), P.palette[0], st.lineWidth, st.alpha * (0.1 + 0.9 * j / m));
      return out + "</g>";
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
  DF.stitchContour = stitch;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/render/player.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Player: renders one scene into a host element. The same Player runs in
   the studio, in the <relab-flow> element and in exported standalone pages.

   A scene is plain JSON and carries its own system text, so it is complete
   without the catalogue:

     { version: 2, name, system, params: {a: 1}, init: {x: 0.1},
       dt, rate, speed, seed, n, initMode, spread,
       perturbations: [{kind, ...}],
       view:  { type, axes: ['x','y'], ranges: {x: [lo, hi]}, ... },
       style: { theme, palette, ramp, colorBy, fade, lineWidth, alpha, pointSize },
       overlay: { title, subtitle, caption, equations, readout, legend } }

   Time runs on the wall clock: a scene plays `rate` units of model time per
   second (iterations per second for a map) times the multiplier `speed`.
   Without a rate the Player chooses one from the model (DF.calibrateRate),
   so that the motion reads alike across models and on any screen refresh
   rate. Each frame advances the whole number of steps the clock has
   accumulated; when those steps would take longer than the frame budget,
   the Player advances fewer and reports the slowdown. */
(function (DF) {
  "use strict";

  const PARAM_FORCING = { periodic: 1, quasiperiodic: 1, ramp: 1, step: 1 };
  const VIEW_DEFAULTS = {
    flow: { n: 1500, life: [80, 320], spawn: "box", fade: 0.06, lineWidth: 1.1, alpha: 0.55, colorBy: "speed" },
    trajectory: { n: 1, fade: 0, lineWidth: 1.3, alpha: 0.95, colorBy: "time", rotate: 0.25 },
    timeseries: { n: 1, fade: 0, lineWidth: 1.6, alpha: 0.95, colorBy: "solid" },
    phase: { n: 6, fade: 0, lineWidth: 1.5, alpha: 0.9, colorBy: "solid" },
    sweep: { n: 1, fade: 0, lineWidth: 1.6, alpha: 0.9, colorBy: "solid" },
    orbit: { n: 1, fade: 0, lineWidth: 1, alpha: 0.35, colorBy: "solid" },
    density: { n: 3000, fade: 0, lineWidth: 1, alpha: 1, colorBy: "solid" },
    strobe: { n: 200, fade: 0, lineWidth: 1, alpha: 0.6, colorBy: "solid", pointSize: 1.4 },
    cobweb: { n: 1, fade: 0, lineWidth: 1.3, alpha: 0.9, colorBy: "solid" }
  };
  const FRAME_MS = 1000 / 60, BUDGET_MS = 12;
  const now = function () { return typeof performance !== "undefined" ? performance.now() : Date.now(); };

  function clone(o) { return JSON.parse(JSON.stringify(o === undefined ? null : o)); }
  function isMapText(src) {
    return String(src || "").split(/\r?\n/).some(function (l) { return DF.STATEMENT.map.test(l.replace(/#.*$/, "").trim()); });
  }

  // Fill defaults; never throws on a partial scene.
  function normalizeScene(input) {
    const s = clone(input || {}) || {};
    s.name = s.name || "Untitled scene";
    s.system = s.system || "x' = -y\ny' = x\ninit x = 1";
    s.params = s.params || {};
    s.init = s.init || {};
    s.view = s.view || {};
    s.view.type = s.view.type && DF.VIEWS[s.view.type] ? s.view.type : "flow";
    const d = VIEW_DEFAULTS[s.view.type];
    // Scenes written before the wall clock advanced a fixed number of steps
    // per frame; at 60 frames per second that is the rate they played at,
    // and a sweep speed given per frame becomes a speed per unit of time.
    if (s.version !== 2 && s.stepsPerFrame > 0) {
      const step = (isMapText(s.system) ? 1 : s.dt || 0.01) * s.stepsPerFrame;
      if (s.rate === undefined && s.view.type !== "orbit") s.rate = +(60 * step).toPrecision(6);
      if (s.view.type === "sweep") s.view.speed = (s.view.speed === undefined ? 0.0006 : s.view.speed) / step;
    }
    delete s.stepsPerFrame;
    s.version = 2;
    if (!(s.rate > 0)) delete s.rate;
    s.speed = s.speed > 0 ? s.speed : 1;
    s.n = s.n || d.n;
    s.seed = s.seed === undefined ? 1 : s.seed;
    s.initMode = s.initMode || "point";
    s.spread = s.spread === undefined ? 0.05 : s.spread;
    s.perturbations = s.perturbations || [];
    // Births and lives of the flow view depend on the dimension of the
    // system, and are settled in buildSim once the system is compiled.
    if (s.view.rotate === undefined && d.rotate !== undefined) s.view.rotate = d.rotate;
    s.style = Object.assign({ theme: "relab-night", palette: "relab", ramp: "relab-fire", colorBy: d.colorBy, fade: d.fade, lineWidth: d.lineWidth, alpha: d.alpha, pointSize: d.pointSize || 1.6, renderScale: 1 }, s.style || {});
    s.overlay = Object.assign({ title: "", subtitle: "", caption: "", equations: false, readout: false, legend: true, position: "top-left" }, s.overlay || {});
    return s;
  }

  /* Ranges, displayed axes and attractor sample of a scene, with the births
     and lives of the flow view settled once the system is known. Shared by
     the player and the playback tests, so that a rate is checked on the
     view as drawn. */
  function resolveView(sys, s, params, init, dt) {
    const discrete = sys.time === "discrete", out = { ranges: null, axes: null, sample: null, fitted: false };
    // Full ranges for every variable: scene, then declared, then a probe run.
    const need = sys.vars.filter(function (v) { return !(s.view.ranges && s.view.ranges[v]) && !sys.ranges[v]; });
    const probe = need.length ? DF.autoRanges(sys, params, init, { dt: dt, steps: discrete ? 3000 : Math.min(20000, Math.max(3000, Math.round(60 / dt))), perturbations: s.perturbations }) : {};
    out.ranges = sys.vars.map(function (v) { return (s.view.ranges && s.view.ranges[v]) || sys.ranges[v] || probe[v]; });
    let axes = (s.view.axes || []).map(function (v) { return sys.vars.indexOf(v); }).filter(function (i) { return i >= 0; });
    if (axes.length < 2) axes = sys.vars.length >= 3 && (s.view.type === "trajectory" || s.view.type === "flow") && s.view.dim3 !== false ? [0, 1, 2] : sys.vars.length >= 2 ? [0, 1] : [0, 0];
    if (s.view.projection === "simplex" && sys.vars.length >= 3 && axes.length < 3) axes = [0, 1, 2];
    // A rotating 3D view of a flow is drawn from a sample of its attractor:
    // flow particles are born on it and live on, the camera faces its
    // broadest side and rocks about it, and the ranges enclose it.
    const v = s.view, solid = axes.length === 3 && axes[0] !== axes[1] && v.projection !== "simplex" && !discrete && (v.type === "flow" || v.type === "trajectory");
    out.sample = null;
    if (solid && (v.fit !== false || v.spawn === "attractor")) {
      out.sample = DF.attractorSample(sys, params, init, { dt: dt, transient: v.sampleTransient, span: v.sampleSpan, perturbations: s.perturbations.filter(function (q) { return q.kind in PARAM_FORCING; }) });
    }
    if (v.type === "flow") {
      if (v.spawn === undefined) v.spawn = out.sample ? "attractor" : VIEW_DEFAULTS.flow.spawn;
      if (v.spawn === "attractor" && !out.sample) v.spawn = "box";
      if (v.life === undefined) v.life = v.spawn === "attractor" ? "inf" : VIEW_DEFAULTS.flow.life.slice();
    }
    // Without chosen axes, a system of more than three variables shows the
    // three that move most on the attractor, each against its own range.
    if (out.sample && !(v.axes && v.axes.length) && sys.vars.length > 3) {
      const S = out.sample, fr = out.ranges, spread = sys.vars.map(function (q, i) {
        let m = 0, m2 = 0;
        for (let k = 0; k < S.n; k++) { const z = S.X[k * S.dim + i]; m += z; m2 += z * z; }
        m /= S.n;
        return Math.sqrt(Math.max(0, m2 / S.n - m * m)) / ((fr[i][1] - fr[i][0]) || 1);
      });
      axes = sys.vars.map(function (q, i) { return i; }).sort(function (i, j) { return spread[j] - spread[i]; }).slice(0, 3).sort(function (i, j) { return i - j; });
    }
    out.fitted = !!(out.sample && v.fit !== false);
    if (out.fitted) {
      // A variable that barely moves on the attractor (a species gone
      // extinct, say) keeps its declared range, so it does not blow up.
      const sr = DF.sampleRanges(out.sample);
      out.ranges = out.ranges.map(function (r, i) {
        if (v.ranges && v.ranges[sys.vars[i]]) return r;
        return sr[i][1] - sr[i][0] > 1e-3 * (r[1] - r[0]) ? sr[i] : r;
      });
    }
    out.axes = axes;
    return out;
  }

  function Player(host, scene, opts) {
    this.opts = opts || {};
    this.host = host;
    this.listeners = {};
    this.running = false; this.frameCount = 0; this.fps = 0;
    this.speedMax = 1e-9; this.speedEnv = null;
    this.frameMs = FRAME_MS; this.spf = 0; this.slow = 1;
    this.buildDom();
    this.watchFonts();
    this.load(scene);
  }

  Player.prototype.on = function (ev, cb) { (this.listeners[ev] = this.listeners[ev] || []).push(cb); return this; };
  Player.prototype.emit = function (ev, data) { (this.listeners[ev] || []).forEach(function (cb) { try { cb(data); } catch (e) { if (typeof console !== "undefined") console.error(e); } }); };

  Player.prototype.buildDom = function () {
    const h = this.host, doc = h.ownerDocument;
    if (getComputedStyle(h).position === "static") h.style.position = "relative";
    h.style.overflow = "hidden";
    const mk = function (tag, cls, css) { const e = doc.createElement(tag); e.className = cls; e.style.cssText = css; h.appendChild(e); return e; };
    const fill = "position:absolute;inset:0;width:100%;height:100%;";
    this.bgEl = mk("div", "df-bg", fill);
    this.cv = { base: mk("canvas", "df-base", fill), trail: mk("canvas", "df-trail", fill), top: mk("canvas", "df-top", fill + "touch-action:none;") };
    this.ctx = { base: this.cv.base.getContext("2d"), trail: this.cv.trail.getContext("2d"), top: this.cv.top.getContext("2d") };
    this.overlayEl = mk("div", "df-overlay", "position:absolute;inset:0;pointer-events:none;font-family:Jost,system-ui,sans-serif;");
    const self = this;
    if (typeof ResizeObserver !== "undefined") { this.ro = new ResizeObserver(function () { self.resize(); }); this.ro.observe(h); }
    if (typeof IntersectionObserver !== "undefined") {
      this.io = new IntersectionObserver(function (e) { self.visible = e[0].isIntersecting; if (self.visible && self.running) self.kick(); });
      this.io.observe(h);
    }
    this.visible = true;
    this.bindPointer();
  };

  // Text drawn on a canvas before its web font arrives stays in the fallback
  // font; the static layer and the overlay are drawn again once it arrives.
  Player.prototype.watchFonts = function () {
    const fonts = this.host.ownerDocument.fonts, self = this;
    if (!fonts || !fonts.load || !fonts.check) return;
    const need = ["11px Jost", "italic 13px 'TeX Gyre Pagella'"];
    let missing;
    try { missing = need.filter(function (f) { return !fonts.check(f); }); } catch (e) { return; }
    if (!missing.length) return;
    Promise.all(missing.map(function (f) { return fonts.load(f).catch(function () { return []; }); })).then(function () {
      if (self.disposed || !self.view) return;
      self.renderOverlay(); self.markStatic();
    });
  };

  Player.prototype.bindPointer = function () {
    const self = this, el = this.cv.top;
    let down = null;
    const pos = function (ev) { const r = el.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    el.addEventListener("pointerdown", function (ev) {
      if (!self.view) return;
      const p = pos(ev); down = { x: p[0], y: p[1], azim: self.cam ? self.cam.azim : 0, elev: self.cam ? self.cam.elev : 0, moved: false };
      el.setPointerCapture(ev.pointerId);
      if (!(self.cam && self.cam.is3D()) && self.view.pointer && self.view.pointer(self, "down", p[0], p[1], ev)) self.emit("interact", {});
    });
    el.addEventListener("pointermove", function (ev) {
      if (!down || !self.view) return;
      const p = pos(ev), dx = p[0] - down.x, dy = p[1] - down.y;
      if (Math.abs(dx) + Math.abs(dy) > 3) down.moved = true;
      if (self.cam && self.cam.is3D()) { self.dragging = true; self.cam.azim = down.azim + dx * 0.01; self.cam.elev = Math.max(-1.5, Math.min(1.5, down.elev + dy * 0.01)); self.clearTrail(); }
      else if (self.view.pointer) self.view.pointer(self, "drag", p[0], p[1], ev);
    });
    const up = function () { down = null; self.dragging = false; };
    el.addEventListener("pointerup", up); el.addEventListener("pointercancel", up);
    el.addEventListener("wheel", function (ev) {
      if (!self.cam || !self.scene || !self.scene.view.zoomable) return;
      ev.preventDefault(); self.cam.zoom = Math.max(0.2, Math.min(20, self.cam.zoom * Math.exp(-ev.deltaY * 0.001))); self.redrawStatic(); self.clearTrail();
    }, { passive: false });
  };

  /* Load (or reload) a scene: compile, build the simulator, the camera and
     the view. A scene whose system does not compile is refused as a whole:
     the Player keeps the scene it was showing, sets `error` and emits it,
     and load returns false. */
  Player.prototype.load = function (scene) {
    const s = normalizeScene(scene);
    let sys;
    try { sys = DF.compileSystem(s.system); }
    catch (e) { this.error = e; this.emit("error", e); return false; }
    const prev = { scene: this.scene, sys: this.sys };
    this.scene = s; this.sys = sys;
    try { this.buildSim(); }
    catch (e) {
      this.scene = prev.scene; this.sys = prev.sys;
      if (prev.scene) { try { this.buildSim(); } catch (e2) { this.view = null; } } else this.view = null;
      this.error = e; this.emit("error", e); return false;
    }
    this.error = null;
    this.emit("scene", s);
    return true;
  };

  // Period of the stroboscopic view: the scene's, else that of a periodic forcing, else 2 pi.
  Player.prototype.strobePeriod = function () {
    const v = this.scene.view, per = (this.scene.perturbations || []).find(function (q) { return q.kind === "periodic" && q.enabled !== false; });
    return v.period || (per ? per.period : 2 * Math.PI);
  };

  /* Rocking of a fitted 3D view: the azimuth swings by this amplitude, in
     radians, about the broadest side. Zero turns rocking off, and a scene
     with a rotation speed and no rocking turns full circles instead. */
  Player.prototype.rockAmplitude = function () {
    const v = this.scene.view;
    if (!this.fitted) return 0;
    if (v.rock !== undefined) return Math.max(0, +v.rock || 0);
    return 0.6;
  };

  /* Default sweep speed, as a fraction of the parameter interval per unit
     of model time: slow against the slowest relaxation at the start, so
     that the state follows its branch until the branch ends at a fold. */
  Player.prototype.defaultSweepSpeed = function (params) {
    const sys = this.sys, discrete = sys.time === "discrete";
    let rate = Infinity;
    try {
      DF.findEquilibria(sys, Float64Array.from(params), this.fullRanges, { seeds: 30 }).forEach(function (e) {
        if (!e.stable) return;
        e.eig.forEach(function (l) { const r = discrete ? -Math.log(Math.max(1e-12, Math.hypot(l.re, l.im))) : -l.re; if (r > 0) rate = Math.min(rate, r); });
      });
    } catch (e) { rate = Infinity; }
    return isFinite(rate) ? Math.min(0.02, Math.max(1e-5, rate / 200)) : 0.005;
  };

  Player.prototype.buildSim = function () {
    const s = this.scene, sys = this.sys, discrete = sys.time === "discrete";
    const params = sys.params.map(function (q) { return s.params[q.name] !== undefined ? +s.params[q.name] : q.value; });
    const init = sys.vars.map(function (v, i) { return s.init[v] !== undefined ? +s.init[v] : sys.init[i]; });
    const dt = s.dt || (discrete ? 1 : 0.01);
    s.dt = dt;
    if (sys.vars.length === 1 && ["cobweb", "density", "timeseries", "sweep", "orbit"].indexOf(s.view.type) < 0) {
      // A scalar system is shown as x against time through the timeseries view.
      s.view.type = "timeseries";
    }
    const R = resolveView(sys, s, params, init, dt);
    this.fullRanges = R.ranges; this.sample = R.sample; this.fitted = R.fitted;
    const axes = R.axes;
    // A stroboscopic view samples every N steps; its step divides the period
    // exactly, h = T / round(T / dt), so that samples fall on t = k T.
    let h = dt;
    const strobe = s.view.type === "strobe" && s.view.mode !== "section" && !discrete;
    if (strobe) {
      const T = this.strobePeriod();
      if (T > 0) h = T / Math.max(1, Math.round(T / dt));
    }
    // Playback rate: the scene's, or one measured from the model.
    this.sweepSpeed = s.view.type === "sweep" ? (s.view.speed > 0 ? s.view.speed : this.defaultSweepSpeed(params)) : 0;
    this.autoRate = null;
    if (!(s.rate > 0)) {
      const vi = Math.max(0, sys.vars.indexOf(s.view.var || sys.vars[0]));
      const shown = s.view.type === "timeseries" || (s.view.type === "density" && sys.vars.length === 1) ? [vi] : axes.filter(function (a, i, arr) { return arr.indexOf(a) === i; });
      this.autoRate = DF.calibrateRate(sys, {
        view: s.view, dt: h, params: params, init: init, ranges: this.fullRanges, axes: shown, sample: this.sample, spread: s.spread, initMode: s.initMode,
        perturbations: s.perturbations, keepPositive: !!s.keepPositive, period: this.strobePeriod(), sweepSpeed: this.sweepSpeed
      });
    }
    this.rate = s.rate > 0 ? s.rate : this.autoRate.rate;
    // At least one step per frame at speed 1: a slower clock would move the
    // figure every second or third frame only. The step is shortened instead.
    if (!discrete && !strobe && this.rate / (60 * h) < 1) h = this.rate / 60;
    this.sim = new DF.Simulator(sys, {
      n: s.n, dt: h, seed: s.seed, params: params, init: init, initMode: s.initMode, spread: s.spread,
      perturbations: s.perturbations, keepPositive: !!s.keepPositive, box: this.fullRanges
    });
    this.tmpX = new Float64Array(sys.vars.length); this.tmpDx = new Float64Array(sys.vars.length);
    this.speedLog = null;
    if (this.sample) {
      // Percentiles of the log speed over the attractor sample, in plot
      // widths per unit time, for the speed colouring.
      const S = this.sample, F = sys.f, d = this.tmpDx, x = new Float64Array(S.dim), fr = this.fullRanges, L = [];
      const H = function (i) { return x[i]; };
      for (let k = 0; k < S.n; k += 2) {
        for (let i = 0; i < S.dim; i++) x[i] = S.X[k * S.dim + i];
        F(0, x, this.sim.p, d, H);
        let q = 0;
        for (let i = 0; i < S.dim; i++) { const w = fr[i][1] - fr[i][0]; q += (d[i] / w) * (d[i] / w); }
        if (q > 0 && isFinite(q)) L.push(0.5 * Math.log(q));
      }
      L.sort(function (a, b) { return a - b; });
      if (L.length > 20) { const lo = L[Math.floor(0.05 * L.length)], hi = L[Math.floor(0.95 * L.length)]; if (hi > lo) this.speedLog = [lo, hi]; }
    }
    const self = this;
    this.cam = new DF.Camera(axes, axes.map(function (i) { return self.fullRanges[i]; }), {
      azim: s.view.azim === undefined && this.fitted ? 0 : s.view.azim, elev: s.view.elev === undefined && this.fitted ? 0.35 : s.view.elev,
      zoom: s.view.zoom, pad: s.view.pad, upAxis: s.view.upAxis, simplex: s.view.projection === "simplex"
    });
    if (this.fitted) this.cam.fitTo(this.sample, { sweep: this.rockAmplitude() > 0 ? this.rockAmplitude() : (s.view.rotate > 0 ? Math.PI : 0) });
    this.rockPhase = 0; this.rockBase = this.cam.azim;
    this.theme = DF.THEMES[s.style.theme] || DF.THEMES["relab-night"];
    this.palette = (DF.PALETTES[s.style.palette] || DF.PALETTES.relab).colors;
    this.acc = 0; this.nCap = Infinity; this.slow = 1;
    this.view = Object.create(DF.VIEWS[s.view.type]);
    this.userParam = false;
    this.frameCount = 0;
    this.speedMax = 1e-9; this.speedEnv = null;
    this.staticDirty = false;
    this.resize(true);
    // A new scene starts on empty layers; otherwise trails of the previous
    // scene remain under views that never draw on the trail layer.
    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);
    this.view.init(this);
    this.redrawStatic();
    this.renderOverlay();
    if (!this.running) this.drawOnce();
  };

  // ------------------------------------------------------------ clock
  // Steps that `ms` milliseconds of playback hold; the fraction carries over.
  Player.prototype.stepsFor = function (ms) {
    this.acc += this.rate * (this.scene.speed || 1) / this.sim.h * ms / 1000;
    const n = Math.floor(this.acc);
    this.acc -= n;
    return n;
  };
  // Steps in `frames` frames at 60 frames per second and speed 1 (warm-ups).
  Player.prototype.nominalSteps = function (frames) { return Math.round(frames * this.rate / (60 * this.sim.h)); };
  Player.prototype.stepsPerSecond = function () { return this.rate * (this.scene.speed || 1) / this.sim.h; };

  Player.prototype.resize = function (skipStatic) {
    const r = this.host.getBoundingClientRect();
    const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
    const dpr = Math.min(3, (typeof devicePixelRatio !== "undefined" ? devicePixelRatio : 1) * (this.scene ? this.scene.style.renderScale || 1 : 1));
    const changed = w !== this.w || h !== this.h || dpr !== this.dpr;
    this.w = w; this.h = h; this.dpr = dpr;
    if (changed) {
      for (const k in this.cv) {
        this.cv[k].width = Math.round(w * dpr); this.cv[k].height = Math.round(h * dpr);
        this.ctx[k].setTransform(dpr, 0, 0, dpr, 0, 0);
      }
    }
    if (!this.scene) return;
    this.bgEl.style.background = DF.backgroundCSS(this.scene.style.theme);
    const axesView = this.view && this.view.axes || (this.scene.view.showAxes && this.cam && !this.cam.is3D());
    const topPad = (this.scene.overlay.title || this.scene.overlay.subtitle) && this.scene.overlay.position !== "none" ? 70 : this.scene.overlay.legend ? 34 : 18;
    this.inset = axesView ? { l: 66, r: 22, t: topPad, b: 52 } : { l: 0, r: 0, t: 0, b: 0 };
    if (this.cam) this.cam.resize(w, h, this.inset);
    if (changed && !skipStatic && this.view) { this.clearTrail(); this.redrawStatic(); this.renderOverlay(); if (!this.running) this.drawOnce(); }
  };

  Player.prototype.clearTrail = function () { this.ctx.trail.clearRect(0, 0, this.w, this.h); };
  // Fade the trail layer by `a` per frame at 60 frames per second, whatever the refresh rate.
  Player.prototype.fade = function (a) {
    if (!(a > 0)) return;
    const eff = 1 - Math.pow(1 - Math.min(1, a), this.frameMs / FRAME_MS);
    const c = this.ctx.trail;
    c.globalCompositeOperation = "destination-out"; c.fillStyle = "rgba(0,0,0," + eff + ")";
    c.fillRect(0, 0, this.w, this.h); c.globalCompositeOperation = "source-over";
  };
  Player.prototype.redrawStatic = function () {
    if (!this.view) return;
    this.ctx.base.clearRect(0, 0, this.w, this.h);
    if (this.scene.view.showAxes && this.cam && !this.cam.is3D() && !this.view.axes) DF.drawAxes(this.ctx.base, this.cam, this.scene.style.theme, [this.sys.vars[this.cam.axes[0]], this.sys.vars[this.cam.axes[1]]], {});
    if (this.scene.view.showAxes && this.cam && this.cam.is3D()) DF.drawBox3D(this.ctx.base, this.cam, this.scene.style.theme);
    if (this.cam && this.cam.simplex && this.scene.view.showAxes !== false) DF.drawSimplex(this.ctx.base, this.cam, this.scene.style.theme, this.cam.axes.map(function (i) { return this.sys.vars[i]; }, this));
    if (this.view.drawStatic) this.view.drawStatic(this);
  };
  // Ask for the static layer to be drawn again: at the next frame while
  // running (so that a dragged slider redraws once per frame), at once when paused.
  Player.prototype.markStatic = function () {
    if (!this.view) return;
    if (this.running) { this.staticDirty = true; return; }
    this.staticDirty = false; this.redrawStatic(); this.drawOnce();
  };

  // ------------------------------------------------------------ overlay
  Player.prototype.textSizes = function () {
    return {
      title: Math.round(Math.max(15, Math.min(26, this.w * 0.024))),
      subtitle: Math.round(Math.max(12, Math.min(16, this.w * 0.015))),
      eq: Math.round(Math.max(11, Math.min(15, this.w * 0.014)))
    };
  };
  Player.prototype.eqTop = function () { return this.view && this.view.axes ? 70 : 16; };
  Player.prototype.bottomPad = function () { return this.view && this.view.axes ? 6 : 12; };
  Player.prototype.legendItems = function () {
    const o = this.scene.overlay;
    if (!o.legend || !this.view || !this.view.legend || o.position === "none") return [];
    return this.view.legend(this) || [];
  };
  // Where the legend goes, in CSS pixels: inside the plot frame when the view
  // has one, at its lower right corner when the equations take the upper right.
  Player.prototype.legendAnchor = function () {
    const o = this.scene.overlay, fb = this.view && this.view.frameBox, ins = this.inset || { r: 0, t: 0, b: 0 };
    const atBottom = !!o.equations;
    if (fb) return atBottom ? { right: this.w - fb.R + 8, bottom: this.h - fb.B + 8 } : { right: this.w - fb.R + 8, top: fb.T + 8 };
    if (atBottom) return { right: ins.r + 10, bottom: this.view && this.view.axes ? ins.b + 8 : (o.readout ? 30 : 12) };
    return { right: ins.r + 10, top: this.view && this.view.axes ? ins.t + 8 : (o.title ? 64 : 12) };
  };
  Player.prototype.renderOverlay = function () {
    const el = this.overlayEl;
    el.innerHTML = "";
    this.readoutEl = null;
    if (!this.scene || !this.view) return;
    const o = this.scene.overlay, th = this.theme, doc = el.ownerDocument, sz = this.textSizes();
    if (o.position === "none") return;
    const box = doc.createElement("div");
    box.style.cssText = "position:absolute;left:18px;top:14px;right:18px;color:" + th.ink + ";";
    if (o.title) { const t = doc.createElement("div"); t.textContent = o.title; t.style.cssText = "font-weight:600;font-size:" + sz.title + "px;line-height:1.15;"; box.appendChild(t); }
    if (o.subtitle) { const t = doc.createElement("div"); t.textContent = o.subtitle; t.style.cssText = "font-weight:300;font-size:" + sz.subtitle + "px;color:" + th.muted + ";margin-top:2px;"; box.appendChild(t); }
    el.appendChild(box);
    if (o.equations) {
      const eq = doc.createElement("div");
      eq.style.cssText = "position:absolute;right:18px;top:" + this.eqTop() + "px;color:" + th.ink + ";font-size:" + sz.eq + "px;text-align:right;opacity:0.92;";
      const k = typeof katex !== "undefined" ? katex : (typeof window !== "undefined" ? window.katex : undefined);
      if (k) {
        DF.systemLatex(this.scene.system).lines.forEach(function (L) {
          const d = doc.createElement("div"); d.style.margin = "2px 0";
          try { k.render(L, d, { throwOnError: false, displayMode: false }); } catch (e) { d.textContent = L; }
          eq.appendChild(d);
        });
      } else {
        // Without KaTeX (a standalone page) the built-in typesetter draws them as SVG.
        const size = sz.eq * 1.21, B = DF.MathType.svgBlock(this.scene.system, 0, 0, size, th.ink);
        const W = Math.ceil(B.width + 4), H = Math.ceil(B.height + size * 0.5);
        eq.innerHTML = '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="' + (-W + 2) + " 0 " + W + " " + H + '">' + B.svg + "</svg>";
      }
      el.appendChild(eq);
    }
    if (o.caption) {
      const c = doc.createElement("div"); c.textContent = o.caption;
      c.style.cssText = "position:absolute;left:18px;bottom:" + this.bottomPad() + "px;max-width:min(560px,70%);font-size:11px;line-height:16px;color:" + th.muted + ";";
      el.appendChild(c);
    }
    if (o.readout) {
      const r = doc.createElement("div");
      r.style.cssText = "position:absolute;right:18px;bottom:" + this.bottomPad() + "px;font:11px ui-monospace,monospace;color:" + th.muted + ";text-align:right;";
      el.appendChild(r); this.readoutEl = r;
    }
    const items = this.legendItems();
    if (items.length) {
      const lg = doc.createElement("div"), A = this.legendAnchor();
      const panelBg = th.dark ? "rgba(11,6,32,0.55)" : "rgba(255,255,255,0.8)";
      lg.style.cssText = "position:absolute;right:" + A.right + "px;" + (A.bottom !== undefined ? "bottom:" + A.bottom + "px;" : "top:" + A.top + "px;") +
        "max-width:48%;display:flex;flex-direction:column;align-items:flex-start;gap:1px;padding:5px 9px;border-radius:6px;background:" + panelBg + ";font-size:11px;line-height:16px;color:" + th.muted + ";";
      items.forEach(function (it) {
        const d = doc.createElement("span"), shape = it[2] || "line", c = it[1];
        const glyph = shape === "dot" ? '<span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:' + c + ';margin-right:5px;vertical-align:-1px"></span>'
          : shape === "ring" ? '<span style="display:inline-block;width:8px;height:8px;border-radius:50%;border:1.5px solid ' + c + ';margin-right:5px;vertical-align:-1px;box-sizing:border-box"></span>'
          : shape === "dash" ? '<span style="display:inline-block;width:16px;border-top:2px dashed ' + c + ';margin-right:5px;vertical-align:middle"></span>'
          : '<span style="display:inline-block;width:16px;height:2px;background:' + c + ';margin-right:5px;vertical-align:middle"></span>';
        d.innerHTML = glyph; d.appendChild(doc.createTextNode(it[0])); lg.appendChild(d);
      });
      el.appendChild(lg);
    }
    if (this.readoutEl) this.readoutEl.textContent = this.readoutText();
  };
  Player.prototype.readoutText = function () {
    const sim = this.sim, sys = this.sys;
    let s = (sys.time === "discrete" ? "n = " + sim.t : "t = " + sim.t.toFixed(2));
    const mod = sim.perturbations.filter(function (q) { return q.enabled !== false && q.param; });
    mod.forEach(function (q) { const i = sim.paramIndex(q.param); if (i >= 0) s += "   " + q.param + " = " + DF.fmt(sim.p[i]); });
    if (this.view.pval !== undefined) s += "   " + sys.params[this.view.pi].name + " = " + DF.fmt(this.view.pval);
    return s;
  };

  // ------------------------------------------------------------ loop
  Player.prototype.drawOnce = function () {
    if (!this.view) return;
    const ms = this.frameMs;
    this.spf = 0; this.frameMs = 0;
    try { this.view.frame(this); } catch (e) { /* first frame of some views needs steps */ }
    this.frameMs = ms;
  };
  /* One frame covering `ms` milliseconds of playback (1000/60 by default).
     While playing live, the steps of a frame are capped so that the frame
     stays within the budget; `slow` is then the fraction of real time kept. */
  Player.prototype.tick = function (ms, live) {
    if (!this.view) return;
    ms = ms === undefined ? FRAME_MS : Math.max(0, Math.min(100, ms));
    this.frameMs = ms;
    let n = this.stepsFor(ms);
    if (live) {
      let slow = 1;
      if (n > this.nCap) { slow = this.nCap / n; n = this.nCap; this.acc = 0; }
      this.slow = 0.9 * this.slow + 0.1 * slow;
    }
    this.spf = n;
    if (this.staticDirty) { this.staticDirty = false; this.redrawStatic(); }
    const t0 = now();
    try { this.view.frame(this); }
    catch (e) { this.pause(); this.error = e; this.emit("error", e); return; }
    if (live) {
      const cost = now() - t0;
      if (cost > BUDGET_MS && n > 1) this.nCap = Math.max(1, Math.floor(n * BUDGET_MS / cost));
      else if (cost < 0.6 * BUDGET_MS && this.nCap < Infinity) this.nCap = this.nCap > 1e7 ? Infinity : Math.ceil(this.nCap * 1.25 + 1);
    }
    this.frameCount++;
    if (this.readoutEl && (this.frameCount % 4) === 0) this.readoutEl.textContent = this.readoutText();
    this.emit("frame", this.frameCount);
  };
  Player.prototype.kick = function () {
    if (this.raf) return;
    const self = this;
    let last = 0, acc = 0, frames = 0;
    const loop = function (ts) {
      self.raf = 0;
      if (!self.running || !self.visible) return;
      if (last) { acc += ts - last; frames++; if (acc > 500) { self.fps = frames * 1000 / acc; acc = 0; frames = 0; } }
      self.tick(last ? ts - last : FRAME_MS, true);
      last = ts;
      self.raf = requestAnimationFrame(loop);
    };
    this.raf = requestAnimationFrame(loop);
  };
  Player.prototype.play = function () { if (this.running || !this.view) return; this.running = true; this.emit("state", "play"); this.kick(); };
  Player.prototype.pause = function () { this.running = false; if (this.raf) cancelAnimationFrame(this.raf); this.raf = 0; this.emit("state", "pause"); };
  Player.prototype.toggle = function () { if (this.running) this.pause(); else this.play(); };
  Player.prototype.step = function () { this.tick(FRAME_MS); };
  Player.prototype.restart = function (seed) {
    if (!this.view) return;
    if (seed !== undefined) this.scene.seed = seed;
    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);
    this.buildSim();
    this.emit("restart", this.scene.seed);
  };
  // Advance n frames of `ms` milliseconds each without waiting for the screen (static renders, exports, tests).
  Player.prototype.advance = function (n, ms) { for (let i = 0; i < n; i++) this.tick(ms); };

  // Change one parameter live, without restarting.
  Player.prototype.setParam = function (name, value) {
    const i = this.sys.params.findIndex(function (q) { return q.name === name; });
    if (i < 0) return;
    this.scene.params[name] = value;
    this.sim.base[i] = value; this.sim.updateParams();
    if (this.view.onParam) this.view.onParam(this, i);
  };
  // Playback: `rate` fixes model time per second (a non-positive rate returns to the automatic one), `speed` multiplies it.
  Player.prototype.setRate = function (rate) {
    if (rate > 0) { this.scene.rate = rate; this.rate = rate; }
    else { delete this.scene.rate; this.rate = this.autoRate ? this.autoRate.rate : this.rate; }
  };
  Player.prototype.setSpeed = function (m) { this.scene.speed = m > 0 ? m : 1; };
  Player.prototype.setStyle = function (patch) {
    Object.assign(this.scene.style, patch);
    this.theme = DF.THEMES[this.scene.style.theme] || DF.THEMES["relab-night"];
    this.palette = (DF.PALETTES[this.scene.style.palette] || DF.PALETTES.relab).colors;
    this.bgEl.style.background = DF.backgroundCSS(this.scene.style.theme);
    if (patch.renderScale !== undefined) { this.w = 0; this.resize(); }
    this.redrawStatic(); this.renderOverlay();
  };
  Player.prototype.setOverlay = function (patch) { Object.assign(this.scene.overlay, patch); this.resize(true); this.redrawStatic(); this.renderOverlay(); };
  Player.prototype.getScene = function () {
    const s = clone(this.scene), sim = this.sim, sys = this.sys;
    sys.params.forEach(function (q, i) { s.params[q.name] = sim.base[i]; });
    if (this.cam && this.cam.is3D()) { s.view.azim = +this.cam.azim.toFixed(4); s.view.elev = +this.cam.elev.toFixed(4); }
    return s;
  };

  // ------------------------------------------------------------ export
  /* Composite of background, layers and overlay text on one canvas, at the
     current resolution, for PNG export and video frames. */
  Player.prototype.composite = function (target) {
    const W = this.cv.top.width, H = this.cv.top.height, c = target || this.host.ownerDocument.createElement("canvas");
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    const g = c.getContext("2d");
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.clearRect(0, 0, W, H);
    g.save(); g.scale(this.dpr, this.dpr); DF.paintBackground(g, this.w, this.h, this.scene.style.theme); g.restore();
    g.drawImage(this.cv.base, 0, 0); g.drawImage(this.cv.trail, 0, 0); g.drawImage(this.cv.top, 0, 0);
    g.save(); g.scale(this.dpr, this.dpr); this.paintText(g); g.restore();
    return c;
  };
  // Lines of the caption wrapped at the width of the overlay caption.
  Player.prototype.captionLines = function (measure) {
    const words = String(this.scene.overlay.caption).split(/\s+/).filter(Boolean), maxW = Math.min(560, this.w * 0.7), lines = [];
    let line = "";
    words.forEach(function (w) { const tt = line ? line + " " + w : w; if (measure(tt) > maxW && line) { lines.push(line); line = w; } else line = tt; });
    if (line) lines.push(line);
    return lines;
  };
  // The legend box, laid out once for the canvas and the SVG painters.
  Player.prototype.legendLayout = function (measure) {
    const items = this.legendItems();
    if (!items.length) return null;
    const A = this.legendAnchor(), rowH = 17, padX = 9, padY = 5, glyphW = 21;
    const w = padX * 2 + glyphW + Math.max.apply(null, items.map(function (it) { return measure(it[0]); }));
    const h = padY * 2 + items.length * rowH;
    const x = this.w - A.right - w, y = A.bottom !== undefined ? this.h - A.bottom - h : A.top;
    return { x: x, y: y, w: w, h: h, rowH: rowH, padX: padX, padY: padY, items: items, bg: this.theme.dark ? "rgba(11,6,32,0.55)" : "rgba(255,255,255,0.8)" };
  };
  Player.prototype.paintText = function (g) {
    const o = this.scene.overlay, th = this.theme, sz = this.textSizes(), self = this;
    if (o.position === "none") return;
    let y = 14;
    g.textBaseline = "top"; g.textAlign = "left";
    if (o.title) { g.fillStyle = th.ink; g.font = "600 " + sz.title + "px Jost, system-ui, sans-serif"; g.fillText(o.title, 18, y); y += Math.round(sz.title * 1.15) + 2; }
    if (o.subtitle) { g.fillStyle = th.muted; g.font = "300 " + sz.subtitle + "px Jost, system-ui, sans-serif"; g.fillText(o.subtitle, 18, y); }
    const bottom = this.bottomPad();
    if (o.caption) {
      g.fillStyle = th.muted; g.font = "11px Jost, system-ui, sans-serif"; g.textBaseline = "bottom";
      const lines = this.captionLines(function (t) { return g.measureText(t).width; });
      lines.forEach(function (L, i) { g.fillText(L, 18, self.h - bottom - (lines.length - 1 - i) * 16); });
    }
    if (o.equations) DF.MathType.paintBlock(g, this.scene.system, this.w - 18, this.eqTop() + 2, sz.eq * 1.21, th.ink);
    if (o.readout) { g.fillStyle = th.muted; g.font = "11px ui-monospace, monospace"; g.textAlign = "right"; g.textBaseline = "bottom"; g.fillText(this.readoutText(), this.w - 18, this.h - bottom); }
    g.font = "11px Jost, system-ui, sans-serif";
    const L = this.legendLayout(function (t) { return g.measureText(t).width; });
    if (L) {
      g.save();
      g.fillStyle = L.bg; g.beginPath();
      if (g.roundRect) g.roundRect(L.x, L.y, L.w, L.h, 6); else g.rect(L.x, L.y, L.w, L.h);
      g.fill();
      g.textAlign = "left"; g.textBaseline = "middle";
      L.items.forEach(function (it, i) {
        const cy = L.y + L.padY + (i + 0.5) * L.rowH, gx = L.x + L.padX, c = it[1], shape = it[2] || "line";
        g.fillStyle = c; g.strokeStyle = c;
        if (shape === "dot") { g.beginPath(); g.arc(gx + 4, cy, 4, 0, 2 * Math.PI); g.fill(); }
        else if (shape === "ring") { g.lineWidth = 1.5; g.beginPath(); g.arc(gx + 4, cy, 3.25, 0, 2 * Math.PI); g.stroke(); }
        else { g.lineWidth = 2; g.setLineDash(shape === "dash" ? [4, 3] : []); g.beginPath(); g.moveTo(gx, cy); g.lineTo(gx + 16, cy); g.stroke(); g.setLineDash([]); }
        g.fillStyle = th.muted; g.fillText(it[0], gx + 21, cy);
      });
      g.restore();
    }
  };

  /* SVG of the current frame. Views that keep vector geometry (trajectory,
     time series, phase plane, sweep, cobweb) return a complete drawing,
     axes included; the others are embedded as one raster image. The text of
     the overlay (title, equations, legend, caption, readout) is SVG text. */
  Player.prototype.toSVG = function () {
    const th = this.theme, o = this.scene.overlay, w = this.w, h = this.h, sz = this.textSizes(), self = this;
    const esc = function (s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); };
    let s = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="' + w + '" height="' + h + '" viewBox="0 0 ' + w + " " + h + '">';
    const bg = th.bg;
    if (bg[0] === "solid") s += '<rect width="100%" height="100%" fill="' + bg[1] + '"/>';
    else if (bg[0] === "radial") s += '<defs><radialGradient id="bg" cx="50%" cy="50%" r="72%"><stop offset="0" stop-color="' + bg[1] + '"/><stop offset="0.55" stop-color="' + bg[2] + '"/><stop offset="1" stop-color="' + bg[3] + '"/></radialGradient></defs><rect width="100%" height="100%" fill="url(#bg)"/>';
    const vec = this.view.svg ? this.view.svg(this) : null;
    if (vec !== null && vec !== undefined) s += '<g id="figure">' + vec + "</g>";
    else {
      const tmp = this.host.ownerDocument.createElement("canvas"); tmp.width = this.cv.top.width; tmp.height = this.cv.top.height;
      const g = tmp.getContext("2d"); g.drawImage(this.cv.base, 0, 0); g.drawImage(this.cv.trail, 0, 0); g.drawImage(this.cv.top, 0, 0);
      s += '<image width="' + w + '" height="' + h + '" xlink:href="' + tmp.toDataURL("image/png") + '"/>';
    }
    if (o.position !== "none") {
      const mc = this.host.ownerDocument.createElement("canvas").getContext("2d");
      const measure = function (font) { return function (t) { mc.font = font; return mc.measureText(t).width; }; };
      let y = 14;
      if (o.title) { s += '<text x="18" y="' + (y + sz.title * 0.82).toFixed(1) + '" font-family="Jost, sans-serif" font-weight="600" font-size="' + sz.title + '" fill="' + th.ink + '">' + esc(o.title) + "</text>"; y += Math.round(sz.title * 1.15) + 2; }
      if (o.subtitle) s += '<text x="18" y="' + (y + sz.subtitle * 0.82).toFixed(1) + '" font-family="Jost, sans-serif" font-weight="300" font-size="' + sz.subtitle + '" fill="' + th.muted + '">' + esc(o.subtitle) + "</text>";
      if (o.equations) s += '<g id="equations">' + DF.MathType.svgBlock(this.scene.system, w - 18, this.eqTop() + 2, sz.eq * 1.21, th.ink).svg + "</g>";
      const bottom = this.bottomPad();
      if (o.caption) {
        const lines = this.captionLines(measure("11px Jost, sans-serif"));
        lines.forEach(function (L, i) { s += '<text x="18" y="' + (h - bottom - 3 - (lines.length - 1 - i) * 16) + '" font-family="Jost, sans-serif" font-size="11" fill="' + th.muted + '">' + esc(L) + "</text>"; });
      }
      if (o.readout) s += '<text x="' + (w - 18) + '" y="' + (h - bottom - 3) + '" text-anchor="end" font-family="ui-monospace, monospace" font-size="11" fill="' + th.muted + '">' + esc(this.readoutText()) + "</text>";
      const L = this.legendLayout(measure("11px Jost, sans-serif"));
      if (L) {
        s += '<g id="legend"><rect x="' + L.x.toFixed(1) + '" y="' + L.y.toFixed(1) + '" width="' + L.w.toFixed(1) + '" height="' + L.h.toFixed(1) + '" rx="6" fill="' + L.bg + '"/>';
        L.items.forEach(function (it, i) {
          const cy = L.y + L.padY + (i + 0.5) * L.rowH, gx = L.x + L.padX, c = it[1], shape = it[2] || "line";
          if (shape === "dot") s += '<circle cx="' + (gx + 4).toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="4" fill="' + c + '"/>';
          else if (shape === "ring") s += '<circle cx="' + (gx + 4).toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="3.25" fill="none" stroke="' + c + '" stroke-width="1.5"/>';
          else s += '<path d="M' + gx.toFixed(1) + " " + cy.toFixed(1) + "h16" + '" stroke="' + c + '" stroke-width="2"' + (shape === "dash" ? ' stroke-dasharray="4 3"' : "") + "/>";
          s += '<text x="' + (gx + 21).toFixed(1) + '" y="' + (cy + 4).toFixed(1) + '" font-family="Jost, sans-serif" font-size="11" fill="' + th.muted + '">' + esc(it[0]) + "</text>";
        });
        s += "</g>";
      }
    }
    void self;
    return s + "</svg>";
  };

  Player.prototype.dispose = function () {
    this.pause();
    this.disposed = true;
    if (this.ro) this.ro.disconnect();
    if (this.io) this.io.disconnect();
    for (const k in this.cv) this.cv[k].remove();
    this.bgEl.remove(); this.overlayEl.remove();
  };

  DF.VIEW_DEFAULTS = VIEW_DEFAULTS;
  DF.resolveView = resolveView;
  DF.normalizeScene = normalizeScene;
  DF.Player = Player;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/render/export.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Exports: still images (PNG, SVG), moving images (WebM, GIF), and the
   scene itself (JSON, share link, embed snippet, standalone page). */
(function (DF) {
  "use strict";

  const E = {};

  E.download = function (blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 4000);
  };
  E.slug = function (s) { return String(s || "scene").toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "scene"; };

  E.png = function (player) {
    return new Promise(function (resolve) { player.composite().toBlob(resolve, "image/png"); });
  };
  E.svg = function (player) { return new Blob([player.toSVG()], { type: "image/svg+xml" }); };

  // ------------------------------------------------------------- WebM
  /* Records the composite of the player for `seconds` at `fps` with the
     browser's MediaRecorder; resolves to a WebM blob. Every video frame
     advances the scene by 1/fps s of playback, so the film runs at the speed
     of the live figure whatever the frame rate. */
  E.webm = function (player, seconds, fps, onProgress) {
    fps = fps || 30;
    const c = player.composite();
    const stream = c.captureStream(fps);
    const type = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm"].find(function (t) { return typeof MediaRecorder !== "undefined" && MediaRecorder.isTypeSupported(t); });
    if (!type) return Promise.reject(new Error("This browser cannot record WebM"));
    const rec = new MediaRecorder(stream, { mimeType: type, videoBitsPerSecond: 12e6 });
    const chunks = [];
    rec.ondataavailable = function (e) { if (e.data && e.data.size) chunks.push(e.data); };
    const wasRunning = player.running;
    player.pause();
    return new Promise(function (resolve, reject) {
      rec.onstop = function () { resolve(new Blob(chunks, { type: "video/webm" })); if (wasRunning) player.play(); };
      rec.onerror = function (e) { reject(e.error || e); };
      rec.start(250);
      const total = Math.round(seconds * fps);
      let k = 0;
      const frame = function () {
        if (k >= total) { rec.stop(); return; }
        player.tick(1000 / fps); player.composite(c); k++;
        if (onProgress) onProgress(k / total);
        setTimeout(frame, 1000 / fps);
      };
      frame();
    });
  };

  // -------------------------------------------------------------- GIF
  /* GIF89a encoder. The palette holds the 256 most frequent colours of a
     5-5-5 bit histogram over all frames; pixels map to the nearest palette
     entry; frames are LZW-compressed with the variable-width code of the
     GIF specification (W3C, GIF89a, 1990, appendix F). */
  function buildPalette(frames) {
    const hist = new Uint32Array(32768);
    frames.forEach(function (f) {
      const d = f.data;
      for (let i = 0; i < d.length; i += 4 * 3) hist[((d[i] >> 3) << 10) | ((d[i + 1] >> 3) << 5) | (d[i + 2] >> 3)]++;
    });
    const idx = [];
    for (let i = 0; i < 32768; i++) if (hist[i]) idx.push(i);
    idx.sort(function (a, b) { return hist[b] - hist[a]; });
    const pal = idx.slice(0, 256).map(function (i) { return [((i >> 10) & 31) * 8 + 4, ((i >> 5) & 31) * 8 + 4, (i & 31) * 8 + 4]; });
    while (pal.length < 256) pal.push([0, 0, 0]);
    return pal;
  }
  function mapper(pal) {
    const cache = new Int16Array(32768).fill(-1);
    return function (r, g, b) {
      const key = ((r >> 3) << 10) | ((g >> 3) << 5) | (b >> 3);
      let c = cache[key];
      if (c >= 0) return c;
      let best = 0, bd = Infinity;
      for (let i = 0; i < 256; i++) { const p = pal[i], d = (p[0] - r) * (p[0] - r) + (p[1] - g) * (p[1] - g) + (p[2] - b) * (p[2] - b); if (d < bd) { bd = d; best = i; } }
      cache[key] = best; return best;
    };
  }
  function lzw(indices, minCode) {
    const out = [];
    let cur = 0, curBits = 0;
    const clear = 1 << minCode, eoi = clear + 1;
    let codeSize = minCode + 1, next = eoi + 1;
    let dict = new Map();
    const emit = function (code) {
      cur |= code << curBits; curBits += codeSize;
      while (curBits >= 8) { out.push(cur & 255); cur >>>= 8; curBits -= 8; }
    };
    emit(clear);
    let prefix = indices[0];
    for (let i = 1; i < indices.length; i++) {
      const k = indices[i], key = prefix * 256 + k, hit = dict.get(key);
      if (hit !== undefined) { prefix = hit; continue; }
      emit(prefix);
      if (next < 4096) {
        dict.set(key, next++);
        if (next > (1 << codeSize) && codeSize < 12) codeSize++;
      } else { emit(clear); dict = new Map(); codeSize = minCode + 1; next = eoi + 1; }
      prefix = k;
    }
    emit(prefix); emit(eoi);
    if (curBits > 0) out.push(cur & 255);
    return out;
  }
  function encodeGIF(frames, w, h, delayCs) {
    const pal = buildPalette(frames), map = mapper(pal);
    const bytes = [];
    const str = function (s) { for (let i = 0; i < s.length; i++) bytes.push(s.charCodeAt(i)); };
    const u16 = function (v) { bytes.push(v & 255, (v >> 8) & 255); };
    str("GIF89a"); u16(w); u16(h); bytes.push(0xf7, 0, 0);
    pal.forEach(function (c) { bytes.push(c[0], c[1], c[2]); });
    bytes.push(0x21, 0xff, 11); str("NETSCAPE2.0"); bytes.push(3, 1, 0, 0, 0);
    frames.forEach(function (f) {
      bytes.push(0x21, 0xf9, 4, 0x04); u16(delayCs); bytes.push(0, 0);
      bytes.push(0x2c); u16(0); u16(0); u16(w); u16(h); bytes.push(0);
      const d = f.data, idx = new Uint8Array(w * h);
      for (let i = 0, j = 0; j < idx.length; i += 4, j++) idx[j] = map(d[i], d[i + 1], d[i + 2]);
      bytes.push(8);
      const data = lzw(idx, 8);
      for (let i = 0; i < data.length; i += 255) { const n = Math.min(255, data.length - i); bytes.push(n); for (let j = 0; j < n; j++) bytes.push(data[i + j]); }
      bytes.push(0);
    });
    bytes.push(0x3b);
    return new Uint8Array(bytes);
  }
  E.encodeGIF = encodeGIF;

  E.gif = function (player, seconds, fps, maxWidth, onProgress) {
    fps = fps || 20; maxWidth = maxWidth || 640;
    const src = player.composite();
    const scale = Math.min(1, maxWidth / src.width);
    const w = Math.max(2, Math.round(src.width * scale)), h = Math.max(2, Math.round(src.height * scale));
    const small = document.createElement("canvas"); small.width = w; small.height = h;
    const g = small.getContext("2d", { willReadFrequently: true });
    const frames = [], total = Math.round(seconds * fps);
    const wasRunning = player.running;
    player.pause();
    return new Promise(function (resolve) {
      let k = 0;
      const step = function () {
        if (k >= total) {
          if (onProgress) onProgress(1, "encoding");
          setTimeout(function () {
            const bytes = encodeGIF(frames, w, h, Math.round(100 / fps));
            if (wasRunning) player.play();
            resolve(new Blob([bytes], { type: "image/gif" }));
          }, 20);
          return;
        }
        player.tick(1000 / fps);
        player.composite(src);
        g.drawImage(src, 0, 0, w, h);
        frames.push(g.getImageData(0, 0, w, h));
        k++;
        if (onProgress) onProgress(k / total, "capturing");
        setTimeout(step, 0);
      };
      step();
    });
  };

  // ---------------------------------------------------------- scenes
  E.sceneJSON = function (scene) { return JSON.stringify(scene, null, 2); };

  function b64url(bytes) {
    let s = "";
    for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
    return btoa(s).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
  }
  function unb64url(s) {
    s = s.replace(/-/g, "+").replace(/_/g, "/"); while (s.length % 4) s += "=";
    const bin = atob(s), out = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }
  // Scene to a compact URL fragment: "z" + deflate-raw + base64url, or "j" + base64url of JSON.
  E.encodeScene = function (scene) {
    const bytes = new TextEncoder().encode(JSON.stringify(scene));
    if (typeof CompressionStream === "undefined") return Promise.resolve("j" + b64url(bytes));
    const cs = new Blob([bytes]).stream().pipeThrough(new CompressionStream("deflate-raw"));
    return new Response(cs).arrayBuffer().then(function (buf) { return "z" + b64url(new Uint8Array(buf)); });
  };
  E.decodeScene = function (code) {
    const kind = code[0], bytes = unb64url(code.slice(1));
    if (kind === "j") return Promise.resolve(JSON.parse(new TextDecoder().decode(bytes)));
    const ds = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
    return new Response(ds).text().then(JSON.parse);
  };

  function escAttr(s) { return String(s).replace(/&/g, "&amp;").replace(/'/g, "&#39;").replace(/</g, "&lt;"); }
  E.embedSnippet = function (scene, src) {
    return '<script src="' + (src || "relabflow.js") + '"></script>\n' +
      "<relab-flow style=\"display:block;width:100%;height:420px\" controls scene='" + escAttr(JSON.stringify(scene)) + "'></relab-flow>";
  };
  // A page that needs nothing else: the engine source, the scene and a full-window player.
  E.standaloneHTML = function (scene, opts) {
    opts = opts || {};
    if (!DF.SOURCE) throw new Error("The standalone export needs the built relabflow.js (run node tools/build.mjs)");
    const th = DF.THEMES[scene.style && scene.style.theme] || DF.THEMES["relab-night"];
    const bg = th.bg[0] === "solid" ? th.bg[1] : th.bg[0] === "radial" ? th.bg[3] : "#000";
    const title = (scene.overlay && scene.overlay.title) || scene.name || "RElabFlow scene";
    return "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<title>" +
      title.replace(/</g, "&lt;") + "</title>\n<link rel=\"icon\" href=\"data:,\">\n<style>html,body{margin:0;height:100%;background:" + bg + "}relab-flow{display:block;width:100vw;height:100vh}</style>\n</head>\n<body>\n" +
      "<relab-flow" + (opts.controls === false ? "" : " controls") + " scene='" + escAttr(JSON.stringify(scene)) + "'></relab-flow>\n" +
      "<script>\n" + DF.SOURCE.replace(/<\/script/gi, "<\\/script") + "\n</script>\n</body>\n</html>\n";
  };

  DF.Export = E;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/models/catalogue.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Model catalogue. Every entry is a system written in the formula language
   of src/core/expr.js plus a default scene; `source` names where the model
   and its parameter values were taken from (an RElab package, file and
   line, or the original literature), so that each entry can be traced.
   The numerical claims in `about` are checked by tests/models.test.mjs. */
(function (DF) {
  "use strict";

  const M = [];
  function add(id, name, group, source, about, system, scene) {
    M.push({ id: id, name: name, group: group, source: source, about: about, system: system.trim().replace(/^ +/gm, ""), scene: scene || {} });
  }
  const PI2 = 2 * Math.PI;

  // =================================================================== ecology
  add("lotka-volterra", "Lotka-Volterra predator and prey", "Ecology",
    "janos R/shiny_app.R:404 (lotka_volterra)",
    "Neutral cycles around the coexistence point: H = eaN - m ln N + aP - r ln P is conserved, so every orbit is closed.",
    `N' = r*N - a*N*P
     P' = e*a*N*P - m*P
     param r = 1 [0.2, 2]
     param a = 0.05 [0.01, 0.2]
     param e = 0.4 [0.1, 1]
     param m = 0.4 [0.05, 1.5]
     init N = 20
     init P = 10
     range N = [0, 80]
     range P = [0, 60]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02, overlay: { equations: true } });

  add("rosenzweig-macarthur", "Rosenzweig-MacArthur", "Ecology",
    "janos R/shiny_app.R:417 (rosenzweig)",
    "Paradox of enrichment: the coexistence equilibrium loses stability in a Hopf bifurcation at K = 3.75, and a limit cycle grows with K.",
    `N' = r*N*(1 - N/K) - a*N*P/(1 + a*h*N)
     P' = e*a*N*P/(1 + a*h*N) - m*P
     param r = 1 [0.2, 2]
     param K = 6 [1, 12]
     param a = 1 [0.2, 2]
     param h = 0.4 [0.1, 1]
     param e = 0.6 [0.1, 1]
     param m = 0.3 [0.05, 0.6]
     init N = 5
     init P = 1
     range N = [0, 7]
     range P = [0, 4.5]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.02, overlay: { equations: true } });

  add("hastings-powell", "Hastings-Powell food chain", "Ecology",
    "kaRma R/demo_system.R:87 and janos R/shiny_app.R:532 (Hastings and Powell 1991)",
    "Three-level food chain with type II responses; chaotic 'teacup' attractor at the classic parameters.",
    `X' = X*(1 - X) - a1*X*Y/(1 + b1*X)
     Y' = a1*X*Y/(1 + b1*X) - a2*Y*Z/(1 + b2*Y) - d1*Y
     Z' = a2*Y*Z/(1 + b2*Y) - d2*Z
     param a1 = 5 [3, 6]
     param b1 = 3 [2, 6.2]
     param a2 = 0.1 [0.05, 0.2]
     param b2 = 2 [1, 3]
     param d1 = 0.4 [0.2, 0.6]
     param d2 = 0.01 [0.001, 0.06]
     init X = 0.8
     init Y = 0.2
     init Z = 8
     range X = [0, 1]
     range Y = [0, 0.5]
     range Z = [7, 10.5]`,
    { view: { type: "flow" }, n: 2000, dt: 0.05, overlay: { equations: true } });

  add("may-leonard", "May-Leonard cyclic competition", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1241 and HiRsch R/systems.R:469 (May and Leonard 1975)",
    "Rock-paper-scissors competition among three species. With alpha > 1 > beta and alpha + beta > 2 orbits approach a heteroclinic cycle through the single-species states and linger ever longer near each.",
    `N1' = N1*(1 - N1 - alpha*N2 - beta*N3)
     N2' = N2*(1 - beta*N1 - N2 - alpha*N3)
     N3' = N3*(1 - alpha*N1 - beta*N2 - N3)
     param alpha = 1.5 [0, 2.5]
     param beta = 0.7 [0, 2.5]
     init N1 = 0.5
     init N2 = 0.3
     init N3 = 0.2
     range N1 = [0, 1]
     range N2 = [0, 1]
     range N3 = [0, 1]`,
    { view: { type: "flow", projection: "simplex", spawn: "mixed", life: [120, 420] }, n: 1800, dt: 0.05, style: { colorBy: "speed", ramp: "relab-fire" } });

  add("competition-lv", "Lotka-Volterra competition", "Ecology",
    "janos R/analysis_phase_portrait.R:207",
    "Two competitors with strong interspecific competition (a12 = a21 = 1.5): the coexistence point is a saddle whose stable manifold separates the two exclusion outcomes.",
    `N1' = r1*N1*(1 - N1/K1 - a12*N2/K1)
     N2' = r2*N2*(1 - N2/K2 - a21*N1/K2)
     param r1 = 1 [0.1, 2]
     param r2 = 1 [0.1, 2]
     param K1 = 100 [50, 150]
     param K2 = 100 [50, 150]
     param a12 = 1.5 [0, 2]
     param a21 = 1.5 [0, 2]
     init N1 = 10
     init N2 = 80
     range N1 = [0, 110]
     range N2 = [0, 110]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.02 });

  add("glv4", "Four-species generalised Lotka-Volterra", "Ecology",
    "janos R/shiny_app.R:491 (glv4)",
    "Competitive community converging to a stable interior equilibrium.",
    `x1' = x1*(1.0 - 1.0*x1 - 0.6*x2 - 0.3*x3 - 0.1*x4)
     x2' = x2*(0.9 - 0.2*x1 - 1.0*x2 - 0.5*x3 - 0.2*x4)
     x3' = x3*(0.8 - 0.1*x1 - 0.3*x2 - 1.0*x3 - 0.4*x4)
     x4' = x4*(0.7 - 0.2*x1 - 0.1*x2 - 0.2*x3 - 1.0*x4)
     init x1 = 0.4
     init x2 = 0.3
     init x3 = 0.2
     init x4 = 0.1
     range x1 = [0, 1]
     range x2 = [0, 1]
     range x3 = [0, 1]
     range x4 = [0, 1]`,
    { view: { type: "timeseries", window: 60, members: 6 }, n: 6, spread: 0.5, dt: 0.05 });

  add("vano-lv4", "Chaotic four-species Lotka-Volterra", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1118 (Vano et al. 2006)",
    "Four competitors with a chaotic attractor, the smallest competitive Lotka-Volterra community known to be chaotic.",
    `x1' = x1*(1 - x1 - 1.09*x2 - 1.52*x3)
     x2' = 0.72*x2*(1 - x2 - 0.44*x3 - 1.36*x4)
     x3' = 1.53*x3*(1 - 2.33*x1 - x3 - 0.47*x4)
     x4' = 1.27*x4*(1 - 1.21*x1 - 0.51*x2 - 0.35*x3 - x4)
     init x1 = 0.3013
     init x2 = 0.4586
     init x3 = 0.1307
     init x4 = 0.3557
     range x1 = [0, 1]
     range x2 = [0, 1]
     range x3 = [0, 0.5]
     range x4 = [0, 1]`,
    { view: { type: "flow", axes: ["x1", "x2", "x3"] }, n: 2000, dt: 0.05 });

  add("act-lv", "Arneodo-Coullet-Tresser Lotka-Volterra", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1088 (Arneodo, Coullet and Tresser 1980)",
    "Three-species Lotka-Volterra system with a chaotic attractor near mu = 1.5.",
    `N1' = N1*(0.5*(1 - N1) + 0.5*(1 - N2) + 0.1*(1 - N3))
     N2' = N2*(-0.5*(1 - N1) - 0.1*(1 - N2) + 0.1*(1 - N3))
     N3' = N3*(mu*(1 - N1) + 0.1*(1 - N2) + 0.1*(1 - N3))
     param mu = 1.52 [1.3, 1.6]
     init N1 = 0.5
     init N2 = 0.3
     init N3 = 0.2`,
    { view: { type: "flow" }, n: 2000, dt: 0.05 });

  add("huisman-weissing", "Huisman-Weissing resource competition", "Ecology",
    "wadaR R/multispecies_competition.R:232 (Huisman and Weissing 2001, Am. Nat. 157: 488)",
    "Five phytoplankton species compete for three resources with Liebig growth; the surviving community depends sensitively on the initial abundances.",
    `aux q1 = max(R1, 0)
     aux q2 = max(R2, 0)
     aux q3 = max(R3, 0)
     aux mu1 = r*min(q1/(0.20 + q1), q2/(0.25 + q2), q3/(0.15 + q3))
     aux mu2 = r*min(q1/(0.05 + q1), q2/(0.10 + q2), q3/(0.95 + q3))
     aux mu3 = r*min(q1/(1.00 + q1), q2/(0.05 + q2), q3/(0.35 + q3))
     aux mu4 = r*min(q1/(0.05 + q1), q2/(1.00 + q2), q3/(0.10 + q3))
     aux mu5 = r*min(q1/(1.20 + q1), q2/(0.40 + q2), q3/(0.05 + q3))
     N1' = N1*(mu1 - m)
     N2' = N2*(mu2 - m)
     N3' = N3*(mu3 - m)
     N4' = N4*(mu4 - m)
     N5' = N5*(mu5 - m)
     R1' = D*(S - R1) - (0.20*mu1*N1 + 0.10*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.10*mu5*N5)
     R2' = D*(S - R2) - (0.10*mu1*N1 + 0.20*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.20*mu5*N5)
     R3' = D*(S - R3) - (0.10*mu1*N1 + 0.10*mu2*N2 + 0.20*mu3*N3 + 0.20*mu4*N4 + 0.10*mu5*N5)
     param r = 1 [0.5, 1.5]
     param m = 0.25 [0.1, 0.4]
     param D = 0.25 [0.1, 0.5]
     param S = 10 [5, 15]
     init N1 = 0.1
     init N2 = 0.1
     init N3 = 0.1
     init N4 = 0.1
     init N5 = 0.1
     init R1 = 10
     init R2 = 10
     init R3 = 10`,
    { view: { type: "timeseries", vars: ["N1", "N2", "N3", "N4", "N5"], window: 800 }, dt: 0.1, keepPositive: true, style: { palette: "relab-qualitative" } });

  add("grazing", "May grazing model", "Ecology",
    "nonautonomeR R/systems.R:1018 (May 1977)",
    "Vegetation under grazing with a type III response: two stable states separated by an unstable one, and hysteresis when the grazing rate c is swept.",
    `V' = r*V*(1 - V/K) - c*V^2/(V^2 + V0^2)
     param r = 1 [0.5, 2]
     param K = 10 [5, 15]
     param c = 2 [1, 3]
     param V0 = 1 [0.5, 2]
     init V = 7
     range V = [0, 10]`,
    { view: { type: "sweep", param: "c", var: "V", from: 1, to: 3, speed: 0.0015 }, dt: 0.05 });

  add("allee", "Strong Allee effect", "Ecology",
    "janos R/analysis_fokker_planck.R:854",
    "Populations below the threshold A decline to extinction; above it they grow to K.",
    `x' = r*x*(x/A - 1)*(1 - x/K)
     param r = 1 [0.2, 2]
     param A = 0.3 [0.05, 0.6]
     param K = 1 [0.7, 1.5]
     init x = 0.5
     range x = [0, 1.2]`,
    { view: { type: "timeseries", window: 20, members: 12 }, n: 12, initMode: "box", dt: 0.02 });

  add("nicholson-bailey", "Nicholson-Bailey host and parasitoid", "Ecology",
    "janos R/shiny_app.R:719 (Nicholson and Bailey 1935)",
    "Host-parasitoid map whose equilibrium is always unstable: oscillations of growing amplitude.",
    `H[n+1] = lambda*H*exp(-a*P)
     P[n+1] = c*H*(1 - exp(-a*P))
     param lambda = 1.5 [1.05, 2]
     param a = 0.02 [0.005, 0.05]
     param c = 1 [0.5, 2]
     init H = 25
     init P = 10
     range H = [0, 200]
     range P = [0, 200]`,
    { view: { type: "trajectory", tail: 60, dim3: false }, style: { pointSize: 4 } });

  add("ricker", "Ricker map", "Ecology",
    "janos R/shiny_app.R:697 (Ricker 1954)",
    "Density-dependent growth of a single population: period doubling from r = 2 and chaos above r of about 2.69.",
    `N[n+1] = N*exp(r*(1 - N/K))
     param r = 2.7 [1.5, 3.5]
     param K = 100 [50, 150]
     init N = 10
     range N = [0, 400]`,
    { view: { type: "orbit", param: "r", var: "N", from: 1.5, to: 3.5 }, style: { alpha: 0.3 } });

  add("seasonal-rm", "Seasonally forced Rosenzweig-MacArthur", "Ecology",
    "nonautonomeR vignettes/pullback-attraction.Rmd:655",
    "Prey growth modulated with period T: the limit cycle entrains at weak forcing and becomes chaotic at eps = 0.8; the stroboscopic section shows the attractor.",
    `x' = x*(1 + eps*sin(2*pi*t/T))*(1 - x) - a*x*y/(b + x)
     y' = y*(a*x/(b + x) - d)
     param a = 1 [0.5, 1.5]
     param b = 0.3 [0.1, 0.5]
     param d = 0.35 [0.2, 0.5]
     param eps = 0.8 [0, 1]
     param T = 10 [5, 20]
     init x = 0.5
     init y = 0.3
     range x = [0, 1]
     range y = [0, 0.8]`,
    { view: { type: "strobe", period: 10, transient: 10 }, n: 600, spread: 0.3, dt: 0.02, style: { pointSize: 1.3, alpha: 0.7 } });

  add("coleman", "Coleman logistic with seasonal carrying capacity", "Ecology",
    "nonautonomeR R/systems.R:1119 and vignettes/coleman-model.Rmd:58",
    "Logistic growth with periodic K(t): all positive orbits converge to one periodic orbit, the pullback attractor.",
    `x' = r*x*(1 - x/(1 + A*sin(w*t)))
     param r = 1 [0.2, 3]
     param A = 0.3 [0, 0.8]
     param w = 0.2 [0.05, 1]
     init x = 0.5
     range x = [0, 2]`,
    { view: { type: "timeseries", window: 60, members: 10 }, n: 10, initMode: "box", dt: 0.02 });

  add("toggle-switch", "Genetic toggle switch", "Ecology",
    "janos vignettes/qualitative-analysis.Rmd:641 (Gardner, Cantor and Collins 2000)",
    "Two mutually repressing genes: two stable states separated by a saddle; noise drives switches between them.",
    `u' = alpha/(1 + v^beta) - u
     v' = alpha/(1 + u^gamma) - v
     param alpha = 3 [1, 6]
     param beta = 2.5 [1, 4]
     param gamma = 2.5 [1, 4]
     init u = 2.5
     init v = 0.5
     range u = [0, 3.5]
     range v = [0, 3.5]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.02 });

  add("rock-paper-scissors", "Rock-paper-scissors replicator", "Ecology",
    "janos R/shiny_app.R:1133 (rps)",
    "Replicator dynamics of the zero-sum rock-paper-scissors game: neutral cycles on the simplex around the mixed equilibrium.",
    `aux f1 = -p2 + p3
     aux f2 = p1 - p3
     aux f3 = -p1 + p2
     aux fbar = p1*f1 + p2*f2 + p3*f3
     p1' = p1*(f1 - fbar)
     p2' = p2*(f2 - fbar)
     p3' = p3*(f3 - fbar)
     init p1 = 0.4
     init p2 = 0.35
     init p3 = 0.25
     range p1 = [0, 1]
     range p2 = [0, 1]
     range p3 = [0, 1]`,
    { view: { type: "flow", projection: "simplex", life: [200, 500] }, n: 1200, dt: 0.02 });

  add("gilpin", "Gilpin spiral chaos", "Ecology",
    "Gilpin (1979), Am. Nat. 113: 306",
    "Two competing prey and one predator: the predator repeatedly falls to very low densities and recovers, and orbits wind around a chaotic spiral (largest Lyapunov exponent about 0.006 per unit time).",
    `x1' = x1*(1 - 0.001*x1 - 0.001*x2 - c*x3)
     x2' = x2*(1 - 0.0015*x1 - 0.001*x2 - 0.001*x3)
     x3' = x3*(-1 + 0.005*x1 + 0.0005*x2)
     param c = 0.01 [0.005, 0.02]
     init x1 = 73.02
     init x2 = 894.2
     init x3 = 1
     range x1 = [0, 1000]
     range x2 = [0, 1000]
     range x3 = [0, 300]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("mccann-yodzis", "McCann-Yodzis food chain", "Ecology",
    "McCann and Yodzis (1994), Ecology 75: 561",
    "Resource, consumer and predator with bioenergetic rates from body size and metabolism; at K = 0.99 the chain is chaotic (largest Lyapunov exponent about 0.009 per unit time).",
    `R' = R*(1 - R/K) - xc*yc*C*R/(R + R0)
     C' = xc*C*(-1 + yc*R/(R + R0)) - xp*yp*P*C/(C + C0)
     P' = xp*P*(-1 + yp*C/(C + C0))
     param K = 0.99 [0.8, 1.1]
     param xc = 0.4 [0.2, 0.6]
     param yc = 2.009 [1.5, 2.5]
     param xp = 0.08 [0.04, 0.12]
     param yp = 2.876 [2, 3.5]
     param R0 = 0.16129 [0.1, 0.3]
     param C0 = 0.5 [0.3, 0.7]
     init R = 0.664
     init C = 0.3832
     init P = 0.7149
     range R = [0.1, 0.95]
     range C = [0.1, 0.6]
     range P = [0.55, 1.1]`,
    { view: { type: "flow" }, n: 2000, dt: 0.02, style: { ramp: "magma" } });

  add("blasius-upca", "Blasius vegetation-herbivore-predator", "Ecology",
    "Blasius, Huppert and Stone (1999), Nature 399: 354",
    "Plants, herbivores and predators with uniform phase and chaotic amplitude: the cycles keep a regular rhythm while their peaks vary chaotically (largest Lyapunov exponent about 0.06 per unit time).",
    `u' = a*u - al1*u*v/(1 + k1*u)
     v' = -b*v + al1*u*v/(1 + k1*u) - al2*v*w
     w' = -c*(w - ws) + al2*v*w
     param a = 1 [0.5, 1.5]
     param b = 1 [0.5, 1.5]
     param c = 10 [5, 15]
     param al1 = 0.2 [0.1, 0.3]
     param al2 = 1 [0.5, 1.5]
     param k1 = 0.05 [0, 0.1]
     param ws = 0.006 [0.002, 0.01]
     init u = 10.11
     init v = 5.216
     init w = 0.01188
     range u = [3, 16]
     range v = [3, 15]
     range w = [0, 1.3]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, style: { ramp: "viridis" } });

  // ============================================================== delays
  add("mackey-glass", "Mackey-Glass", "Delay equations",
    "janos R/shiny_app.R:733 (Mackey and Glass 1977)",
    "Blood-cell production with a delay of 17 time units: chaotic oscillations.",
    `x' = a*lag(x, tau)/(1 + lag(x, tau)^n) - b*x
     param a = 0.2 [0.1, 0.3]
     param b = 0.1 [0.05, 0.2]
     param n = 10 [4, 12]
     param tau = 17 [2, 30]
     init x = 0.9
     range x = [0.2, 1.5]`,
    { view: { type: "timeseries", window: 600 }, dt: 0.1 });

  add("hutchinson", "Hutchinson delayed logistic", "Delay equations",
    "janos R/shiny_app.R:743 (Hutchinson 1948)",
    "Logistic growth with delayed feedback: the equilibrium K loses stability when r tau exceeds pi/2.",
    `N' = r*N*(1 - lag(N, tau)/K)
     param r = 1.6 [0.5, 2.5]
     param K = 100 [50, 150]
     param tau = 1 [0.2, 2]
     init N = 20
     range N = [0, 320]`,
    { view: { type: "timeseries", window: 40 }, dt: 0.02 });

  add("nicholson-blowflies", "Nicholson blowflies", "Delay equations",
    "janos vignettes/introduction.Rmd:265 and symplectoR R/data.R:123 (Gurney, Blythe and Nisbet 1980)",
    "Delayed recruitment with a hump-shaped birth function: large irregular population cycles.",
    `N' = P*lag(N, tau)*exp(-lag(N, tau)/N0) - delta*N
     param P = 8 [2, 12]
     param N0 = 1 [0.5, 2]
     param delta = 0.175 [0.1, 0.4]
     param tau = 15 [5, 20]
     init N = 3
     range N = [0, 20]`,
    { view: { type: "timeseries", window: 300 }, dt: 0.1 });

  add("delayed-predator-prey", "Delayed predator and prey", "Delay equations",
    "janos vignettes/qualitative-analysis.Rmd:768",
    "Prey self-regulation acts with a delay tau, which destabilises coexistence into cycles in a Hopf bifurcation at tau = 1.437. At tau = 1.5 the prey cycles between about 1.7 and 10; at the value 3 of the source the cycles reach N of about 330, far outside the plot, so tau = 1.5 is the default here.",
    `N' = r*N*(1 - lag(N, tau)/K) - a*N*P
     P' = b*N*P - d*P
     param r = 1.5 [0.5, 2.5]
     param K = 10 [5, 15]
     param a = 0.2 [0.1, 0.4]
     param b = 0.1 [0.05, 0.2]
     param d = 0.5 [0.2, 1]
     param tau = 1.5 [0.5, 5]
     init N = 5
     init P = 2
     range N = [0, 14]
     range P = [0, 12]`,
    { view: { type: "trajectory" }, dt: 0.02 });

  // ================================================================ chaos
  add("lorenz", "Lorenz", "Chaotic flows",
    "janos R/shiny_app.R:506 and tuRbulence R/dynamical_systems.R:61 (Lorenz 1963)",
    "Convection model with the butterfly attractor; largest Lyapunov exponent 0.906 at the classical parameters.",
    `x' = sigma*(y - x)
     y' = x*(rho - z) - y
     z' = x*y - beta*z
     param sigma = 10 [1, 20]
     param rho = 28 [0.5, 50]
     param beta = 2.6666666666666665 [0.5, 4]
     init x = 1
     init y = 1
     init z = 1
     range x = [-22, 22]
     range y = [-28, 28]
     range z = [0, 52]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, overlay: { equations: true } });

  add("rossler", "Rossler", "Chaotic flows",
    "janos R/shiny_app.R:520 and tuRbulence R/dynamical_systems.R:254 (Rossler 1976)",
    "Band chaos from a single folded band; largest Lyapunov exponent about 0.07 at a = b = 0.2, c = 5.7.",
    `x' = -y - z
     y' = x + a*y
     z' = b + z*(x - c)
     param a = 0.2 [0, 0.4]
     param b = 0.2 [0, 1]
     param c = 5.7 [2, 12]
     init x = 1
     init y = 1
     init z = 1
     range x = [-12, 14]
     range y = [-14, 11]
     range z = [0, 24]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01 });

  add("chua", "Chua double scroll", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:167 (Chua, Komuro and Matsumoto 1986)",
    "Electronic circuit with a piecewise-linear diode: the double-scroll attractor.",
    `x' = alpha*(y - x - (m1*x + 0.5*(m0 - m1)*(abs(x + 1) - abs(x - 1))))
     y' = x - y + z
     z' = -beta*y
     param alpha = 15.6 [8, 20]
     param beta = 28 [20, 35]
     param m0 = -1.143 [-1.5, -0.8]
     param m1 = -0.714 [-1, -0.4]
     init x = 0.1
     init y = 0
     init z = 0
     range x = [-2.6, 2.6]
     range y = [-0.5, 0.5]
     range z = [-4, 4]`,
    { view: { type: "flow", axes: ["x", "z", "y"] }, n: 2000, dt: 0.005 });

  add("chen", "Chen", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:204 (Chen and Ueta 1999)",
    "A Lorenz-like system with a double-wing attractor of different topology.",
    `x' = a*(y - x)
     y' = (c - a)*x - x*z + c*y
     z' = x*y - b*z
     param a = 35 [30, 40]
     param b = 3 [1, 5]
     param c = 28 [20, 30]
     init x = -10
     init y = 0
     init z = 37
     range x = [-30, 30]
     range y = [-32, 32]
     range z = [0, 60]`,
    { view: { type: "flow" }, n: 2000, dt: 0.001 });

  add("lu", "Lu", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:231 (Lu and Chen 2002)",
    "Intermediate between the Lorenz and Chen attractors.",
    `x' = a*(y - x)
     y' = -x*z + c*y
     z' = x*y - b*z
     param a = 36 [30, 40]
     param b = 3 [1, 5]
     param c = 20 [12, 28]
     init x = 0.1
     init y = 0.2
     init z = 0.3
     range x = [-25, 25]
     range y = [-28, 28]
     range z = [0, 45]`,
    { view: { type: "flow" }, n: 2000, dt: 0.002 });

  add("shimizu-morioka", "Shimizu-Morioka", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:269 (Shimizu and Morioka 1980)",
    "Lorenz-type attractor of a laser model.",
    `x' = y
     y' = x - lambda*y - x*z
     z' = -alpha*z + x^2
     param alpha = 0.45 [0.3, 0.6]
     param lambda = 0.75 [0.5, 1]
     init x = 0.1
     init y = 0.1
     init z = 0.1
     range x = [-2, 2]
     range y = [-1.6, 1.6]
     range z = [0, 2.4]`,
    { view: { type: "flow" }, n: 2000, dt: 0.02 });

  add("nose-hoover", "Nose-Hoover (Sprott A)", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:301 (Hoover 1985)",
    "Thermostatted oscillator: a conservative flow in which a chaotic sea coexists with invariant tori.",
    `x' = y
     y' = -x - z*y
     z' = y^2 - a
     param a = 1 [0.5, 2]
     init x = 0
     init y = 5
     init z = 0
     range x = [-4, 4]
     range y = [-5, 5]
     range z = [-4, 4]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01 });

  add("sprott-jerk", "Sprott minimal jerk", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:343 (Sprott 1997)",
    "Third-order equation x''' = -a x'' + x'^2 - x, one of the simplest chaotic flows.",
    `x' = y
     y' = z
     z' = -a*z + y^2 - x
     param a = 2.017 [1.9, 2.1]
     init x = 0
     init y = 0
     init z = 1
     range x = [-7, 5]
     range y = [-3, 3]
     range z = [-4, 3]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01 });

  add("thomas", "Thomas cyclically symmetric", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:375 (Thomas 1999)",
    "Symmetric flow driven by sines; small damping b gives a labyrinth of chaotic motion.",
    `x' = sin(y) - b*x
     y' = sin(z) - b*y
     z' = sin(x) - b*z
     param b = 0.18 [0.05, 0.3]
     init x = 2.4
     init y = 2.5
     init z = 2.6
     range x = [-5, 5]
     range y = [-5, 5]
     range z = [-5, 5]`,
    { view: { type: "flow", life: [200, 600], rotate: 0.15 }, n: 1600, dt: 0.05, style: { colorBy: "speed", ramp: "relab-fire" } });

  add("halvorsen", "Halvorsen", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:424",
    "Cyclically symmetric attractor with three lobes.",
    `x' = -a*x - 4*y - 4*z - y^2
     y' = -a*y - 4*z - 4*x - z^2
     z' = -a*z - 4*x - 4*y - x^2
     param a = 1.4 [1.2, 1.6]
     init x = 1
     init y = 0
     init z = 0
     range x = [-12, 8]
     range y = [-12, 8]
     range z = [-12, 8]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005 });

  add("aizawa", "Aizawa", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:457",
    "A sphere-like attractor with a tube along its axis.",
    `x' = (z - b)*x - d*y
     y' = d*x + (z - b)*y
     z' = c + a*z - z^3/3 - (x^2 + y^2)*(1 + e*z) + f*z*x^3
     param a = 0.95 [0.7, 1]
     param b = 0.7 [0.5, 0.9]
     param c = 0.6 [0.4, 0.8]
     param d = 3.5 [2, 5]
     param e = 0.25 [0, 0.5]
     param f = 0.1 [0, 0.3]
     init x = 0.1
     init y = 0
     init z = 0
     range x = [-1.6, 1.6]
     range y = [-1.6, 1.6]
     range z = [-0.6, 2]`,
    { view: { type: "flow", life: [150, 500], rotate: 0.2 }, n: 1500, dt: 0.01, style: { colorBy: "speed", ramp: "mako" } });

  add("rabinovich-fabrikant", "Rabinovich-Fabrikant", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:490 (Rabinovich and Fabrikant 1979)",
    "Modulation instability in a non-equilibrium medium; strongly stretched attractor.",
    `x' = y*(z - 1 + x^2) + gamma*x
     y' = x*(3*z + 1 - x^2) + gamma*y
     z' = -2*z*(alpha + x*y)
     param gamma = 0.87 [0.1, 1]
     param alpha = 1.1 [0.9, 1.3]
     init x = -1
     init y = 0
     init z = 0.5
     range x = [-2.5, 2.5]
     range y = [-3, 3]
     range z = [0, 2]`,
    { view: { type: "flow" }, n: 2000, dt: 0.002 });

  add("lorenz-84", "Lorenz-84 atmosphere", "Chaotic flows",
    "tuRbulence R/dynamical_systems.R:447 and nonautonomeR R/systems.R:651 (Lorenz 1984)",
    "Low-order model of the westerlies (X) and a travelling wave (Y, Z), chaotic at F = 8, G = 1.",
    `X' = -Y^2 - Z^2 - a*X + a*F
     Y' = X*Y - b*X*Z - Y + G
     Z' = b*X*Y + X*Z - Z
     param a = 0.25 [0.1, 0.5]
     param b = 4 [2, 6]
     param F = 8 [4, 10]
     param G = 1 [0, 2]
     init X = 1
     init Y = 1
     init Z = 1
     range X = [-1.5, 2.8]
     range Y = [-2.5, 2.8]
     range Z = [-2.6, 2.6]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01 });

  add("lorenz-96", "Lorenz-96, five sites", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:653 (Lorenz 1996)",
    "Ring of five sites with advection, damping and forcing F: spatiotemporal chaos.",
    `x1' = (x2 - x4)*x5 - x1 + F
     x2' = (x3 - x5)*x1 - x2 + F
     x3' = (x4 - x1)*x2 - x3 + F
     x4' = (x5 - x2)*x3 - x4 + F
     x5' = (x1 - x3)*x4 - x5 + F
     param F = 8 [2, 12]
     init x1 = 8.01
     init x2 = 8
     init x3 = 8
     init x4 = 8
     init x5 = 8
     range x1 = [-8, 13]
     range x2 = [-8, 13]
     range x3 = [-8, 13]
     range x4 = [-8, 13]
     range x5 = [-8, 13]`,
    { view: { type: "flow", axes: ["x1", "x2", "x3"] }, n: 2000, dt: 0.005 });

  add("hyperchaotic-rossler", "Hyperchaotic Rossler", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:974 (Rossler 1979)",
    "Four-dimensional flow with two positive Lyapunov exponents.",
    `x' = -y - z
     y' = x + a*y + w
     z' = b + x*z
     w' = -c*z + d*w
     param a = 0.25 [0.2, 0.3]
     param b = 3 [2, 4]
     param c = 0.5 [0.3, 0.7]
     param d = 0.05 [0.02, 0.08]
     init x = -10
     init y = -6
     init z = 0
     init w = 10
     range x = [-60, 40]
     range y = [-40, 40]
     range z = [0, 100]
     range w = [0, 120]`,
    { view: { type: "flow", axes: ["x", "y", "w"] }, n: 2000, dt: 0.005 });

  add("newton-leipnik", "Newton-Leipnik", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:1283 (Leipnik and Newton 1981)",
    "Rigid-body motion with feedback: two coexisting strange attractors, reached from z0 = -0.16 and z0 = -0.18.",
    `x' = -a*x + y + 10*y*z
     y' = -x - a*y + 5*x*z
     z' = alpha*z - 5*x*y
     param a = 0.4 [0.3, 0.5]
     param alpha = 0.175 [0.1, 0.25]
     init x = 0.349
     init y = 0
     init z = -0.16
     range x = [-0.6, 0.6]
     range y = [-0.6, 0.6]
     range z = [-0.6, 0.1]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01 });

  add("charney-devore", "Charney-DeVore three modes", "Chaotic flows",
    "tuRbulence R/charney_devore.R:84 (Charney and DeVore 1979)",
    "Truncated barotropic flow over topography. At F = 4 a zonal state (weak waves) and a blocked state (strong waves) are both stable, separated by a saddle; particles released over the box settle into one regime or the other. Below F of about 3.4 the blocked state is unstable, and below about 2.5 it does not exist.",
    `x' = k*(F - x) - alpha*y*z + beta*y
     y' = -k*y + alpha*x*z - beta*x - delta*z
     z' = -k*z + delta*y
     param F = 4 [0.5, 5]
     param k = 0.1 [0.05, 0.3]
     param alpha = 1 [0.5, 1.5]
     param beta = 0.5 [0.2, 1]
     param delta = 1 [0.5, 1.5]
     init x = 1
     init y = 0.1
     init z = 0.1
     range x = [0, 4.5]
     range y = [-1, 1]
     range z = [-2.5, 2.8]`,
    { view: { type: "flow", life: [240, 720], rotate: 0.15 }, n: 1200, dt: 0.02, style: { colorBy: "speed", ramp: "relab-fire" } });

  add("rucklidge", "Rucklidge double convection", "Chaotic flows",
    "Rucklidge (1992), J. Fluid Mech. 237: 209",
    "Truncated model of two-dimensional convection in a magnetic field or a rotating layer; a Lorenz-like attractor with largest Lyapunov exponent about 0.19 at k = 2, l = 6.7.",
    `x' = -k*x + l*y - y*z
     y' = x
     z' = -z + y^2
     param k = 2 [1, 3]
     param l = 6.7 [5, 8]
     init x = 6.708
     init y = -3.317
     init z = 13.07
     range x = [-11, 11]
     range y = [-6, 6]
     range z = [0, 16]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("burke-shaw", "Burke-Shaw", "Chaotic flows",
    "Shaw (1981), Z. Naturforsch. 36a: 80",
    "A Lorenz-type flow with the symmetry (x, y, z) -> (-x, -y, z); largest Lyapunov exponent about 0.7 at s = 10, v = 4.272.",
    `x' = -s*(x + y)
     y' = -y - s*x*z
     z' = s*x*y + v
     param s = 10 [8, 12]
     param v = 4.272 [3.5, 5]
     init x = -0.329
     init y = 0.5766
     init z = 1.617
     range x = [-1.8, 1.8]
     range y = [-2.4, 2.4]
     range z = [-2, 2]`,
    { view: { type: "flow" }, n: 2000, dt: 0.002, style: { ramp: "inferno" } });

  add("dadras", "Dadras-Momeni", "Chaotic flows",
    "Dadras and Momeni (2009), Phys. Lett. A 373: 3637",
    "Three-dimensional quadratic flow whose attractor has two, three or four scrolls as the parameters change; chaotic at the default values (largest Lyapunov exponent about 0.4).",
    `x' = y - a*x + b*y*z
     y' = c*y - x*z + z
     z' = d*x*y - e*z
     param a = 3 [2, 4]
     param b = 2.7 [2, 3.5]
     param c = 1.7 [1, 2.5]
     param d = 2 [1, 3]
     param e = 9 [6, 12]
     init x = 0.4958
     init y = 2.025
     init z = 0.1494
     range x = [-20, 22]
     range y = [-13, 9]
     range z = [-16, 13]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, style: { ramp: "magma" } });

  add("chen-lee", "Chen-Lee rigid body", "Chaotic flows",
    "Chen and Lee (2004), Chaos Solitons Fractals 21: 957",
    "Rigid body under feedback torques: two symmetric wings around z > 0; largest Lyapunov exponent about 0.14.",
    `x' = a*x - y*z
     y' = b*y + x*z
     z' = c*z + x*y/3
     param a = 5 [3, 6]
     param b = -10 [-12, -8]
     param c = -0.38 [-0.6, -0.2]
     init x = -0.5446
     init y = -0.373
     init z = 9.408
     range x = [-16, 16]
     range y = [-13, 13]
     range z = [3, 15]`,
    { view: { type: "flow" }, n: 2000, dt: 0.002, style: { ramp: "relab-fire" } });

  add("sprott-b", "Sprott B", "Chaotic flows",
    "Sprott (1994), Phys. Rev. E 50: R647",
    "One of nineteen chaotic flows with five terms and two quadratic nonlinearities found by computer search; largest Lyapunov exponent about 0.21.",
    `x' = y*z
     y' = x - y
     z' = 1 - x*y
     init x = 1.338
     init y = 1.307
     init z = -0.4326
     range x = [-8, 9]
     range y = [-4, 4.5]
     range z = [-8.5, 8]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "blackboard" } });

  add("sprott-c", "Sprott C", "Chaotic flows",
    "Sprott (1994), Phys. Rev. E 50: R647",
    "Sprott's case C, x' = yz, y' = x - y, z' = 1 - x^2; largest Lyapunov exponent about 0.16.",
    `x' = y*z
     y' = x - y
     z' = 1 - x^2
     init x = -0.6523
     init y = -0.5747
     init z = 0.8349
     range x = [-4.2, 4.2]
     range y = [-2.5, 2.5]
     range z = [-8.5, 6]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "mako" } });

  add("sprott-d", "Sprott D", "Chaotic flows",
    "Sprott (1994), Phys. Rev. E 50: R647",
    "Sprott's case D, x' = -y, y' = x + z, z' = xz + 3y^2; largest Lyapunov exponent about 0.10.",
    `x' = -y
     y' = x + z
     z' = x*z + 3*y^2
     init x = -1.975
     init y = -1.611
     init z = 3.625
     range x = [-5, 0]
     range y = [-2.6, 2.1]
     range z = [0, 5]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, style: { ramp: "inferno" } });

  add("sprott-e", "Sprott E", "Chaotic flows",
    "Sprott (1994), Phys. Rev. E 50: R647",
    "Sprott's case E, x' = yz, y' = x^2 - y, z' = 1 - 4x; largest Lyapunov exponent about 0.08.",
    `x' = y*z
     y' = x^2 - y
     z' = 1 - 4*x
     init x = 0.02179
     init y = 0.06937
     init z = 2.994
     range x = [-1.6, 2.4]
     range y = [0, 2.5]
     range z = [-2.8, 6]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "viridis" } });

  add("genesio-tesi", "Genesio-Tesi jerk", "Chaotic flows",
    "Genesio and Tesi (1992), Automatica 28: 531",
    "Third-order equation x''' = -a x'' - b x' - c x + x^2 with a single quadratic term; a folded band attractor, largest Lyapunov exponent about 0.11.",
    `x' = y
     y' = z
     z' = -c*x - b*y - a*z + x^2
     param a = 0.44 [0.3, 0.6]
     param b = 1.1 [0.9, 1.3]
     param c = 1 [0.8, 1.2]
     init x = 0.169
     init y = -0.03774
     init z = -0.1851
     range x = [-0.55, 1.25]
     range y = [-0.7, 1]
     range z = [-0.85, 1]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("arneodo", "Arneodo-Coullet-Tresser spiral", "Chaotic flows",
    "Arneodo, Coullet and Tresser (1981), Commun. Math. Phys. 79: 573",
    "Jerk flow x''' = 5.5 x - 3.5 x' - x'' - x^3 with a Shilnikov spiral around two symmetric foci; largest Lyapunov exponent about 0.24.",
    `x' = y
     y' = z
     z' = a*x - b*y - z - x^3
     param a = 5.5 [4, 7]
     param b = 3.5 [3, 4]
     init x = 0.8461
     init y = 0.8614
     init z = 0.9103
     range x = [-3.8, 3.8]
     range y = [-6.5, 6.5]
     range z = [-12, 12]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "magma" } });

  add("liu-chen", "Liu-Chen", "Chaotic flows",
    "Liu, Liu, Liu and Zhou (2004), Chaos Solitons Fractals 22: 1031",
    "Lorenz-type system with a quadratic term x^2 in the third equation; a two-wing attractor with largest Lyapunov exponent about 1.6.",
    `x' = a*(y - x)
     y' = b*x - k*x*z
     z' = -c*z + h*x^2
     param a = 10 [5, 15]
     param b = 40 [30, 50]
     param c = 2.5 [1.5, 3.5]
     param k = 1 [0.5, 1.5]
     param h = 4 [2, 6]
     init x = 1.672
     init y = -9.284
     init z = 77
     range x = [-18, 18]
     range y = [-34, 34]
     range z = [5, 110]`,
    { view: { type: "flow" }, n: 2000, dt: 0.002, style: { ramp: "inferno" } });

  add("moore-spiegel", "Moore-Spiegel oscillator", "Chaotic flows",
    "Moore and Spiegel (1966), Astrophys. J. 143: 871",
    "Thermally excited fluid element in a stratified, rotating star; chaotic at T = 6, R = 20 (largest Lyapunov exponent about 0.10).",
    `x' = y
     y' = z
     z' = -z - (T - R + R*x^2)*y - T*x
     param T = 6 [4, 8]
     param R = 20 [15, 25]
     init x = 0.2981
     init y = 0.6634
     init z = 1.866
     range x = [-2.6, 2.5]
     range y = [-7, 6.2]
     range z = [-48, 61]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, style: { ramp: "mako" } });

  add("finance", "Ma-Chen finance system", "Chaotic flows",
    "Ma and Chen (2001), Appl. Math. Mech. 22: 1375",
    "Interest rate x, investment demand y and price index z of a macroeconomic model; chaotic at a = 0.00001, b = 0.1, c = 1 (largest Lyapunov exponent about 0.13).",
    `x' = z + (y - a)*x
     y' = 1 - b*y - x^2
     z' = -x - c*z
     param a = 0.00001 [0, 1]
     param b = 0.1 [0.05, 0.3]
     param c = 1 [0.5, 1.5]
     init x = 0.9893
     init y = -0.6592
     init z = -1.385
     range x = [-3.8, 3.6]
     range y = [-2.7, 3.7]
     range z = [-1.8, 1.8]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("rikitake", "Rikitake dynamo", "Chaotic flows",
    "Rikitake (1958), Math. Proc. Camb. Philos. Soc. 54: 89",
    "Two coupled disk dynamos; the currents reverse irregularly, as the geomagnetic field does (largest Lyapunov exponent about 0.13 at m = 2, a = 5).",
    `x' = -m*x + y*z
     y' = -m*y + x*(z - a)
     z' = 1 - x*y
     param m = 2 [1, 3]
     param a = 5 [3, 7]
     init x = 0.1385
     init y = 0.06071
     init z = 6.326
     range x = [-6, 7]
     range y = [-3, 3.7]
     range z = [2, 10]`,
    { view: { type: "flow" }, n: 2000, dt: 0.005, style: { ramp: "blackboard" } });

  add("windmi", "WINDMI magnetosphere", "Chaotic flows",
    "Horton and Doxas (1998), J. Geophys. Res. 103: 4561",
    "Jerk reduction of the solar-wind-driven magnetosphere and ionosphere model, x''' = -a x'' - x' + b - exp(x); largest Lyapunov exponent about 0.08.",
    `x' = y
     y' = z
     z' = -a*z - y + b - exp(x)
     param a = 0.7 [0.5, 0.9]
     param b = 2.5 [2, 3]
     init x = 1.516
     init y = 2.936
     init z = -3.169
     range x = [-27, 4]
     range y = [-16, 9]
     range z = [-20, 12]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "viridis" } });

  add("willamowski-rossler", "Willamowski-Rossler chemistry", "Chaotic flows",
    "Willamowski and Rossler (1980), Z. Naturforsch. 35a: 317",
    "Mass-action network of three species with only bimolecular reactions, the first realistic chemical scheme shown to be chaotic (largest Lyapunov exponent about 0.53).",
    `x' = k1*x - km1*x^2 - k2*x*y + km2*y^2 - k4*x*z + km4
     y' = k2*x*y - km2*y^2 - k3*y + km3
     z' = -k4*x*z + km4 + k5*z - km5*z^2
     param k1 = 30 [25, 35]
     param km1 = 0.25 [0.2, 0.3]
     param k2 = 1
     param km2 = 0.0001
     param k3 = 10 [8, 12]
     param km3 = 0.001
     param k4 = 1
     param km4 = 0.001
     param k5 = 16.5 [14, 19]
     param km5 = 0.5 [0.4, 0.6]
     init x = 0.4932
     init y = 57.25
     init z = 0.05962
     range x = [0, 95]
     range y = [0, 100]
     range z = [0, 30]`,
    { view: { type: "flow" }, n: 2000, dt: 0.0005, style: { ramp: "inferno" } });

  // ============================================================== forced
  add("duffing", "Forced Duffing oscillator", "Forced oscillators",
    "janos R/shiny_app.R:547 (duffing)",
    "Double-well oscillator driven with period 2 pi / 1.2: the stroboscopic map reveals a strange attractor.",
    `x' = y
     y' = -delta*y - alpha*x - beta*x^3 + gamma*cos(omega*t)
     param delta = 0.3 [0.1, 0.5]
     param alpha = -1 [-1.5, 1]
     param beta = 1 [0.5, 1.5]
     param gamma = 0.5 [0.2, 0.6]
     param omega = 1.2 [0.8, 1.6]
     init x = 0.1
     init y = 0
     range x = [-2, 2]
     range y = [-1.6, 1.6]`,
    { view: { type: "strobe", period: PI2 / 1.2, transient: 5 }, n: 500, spread: 1, dt: 0.02, style: { pointSize: 1.2, alpha: 0.6 } });

  add("forced-van-der-pol", "Forced Van der Pol", "Forced oscillators",
    "janos vignettes/chaotic-systems.Rmd:532",
    "Relaxation oscillator driven at a frequency near its own: chaos between locking regimes.",
    `x' = y
     y' = mu*(1 - x^2)*y - x + A*sin(omega*t)
     param mu = 3 [0.5, 5]
     param A = 5 [0, 8]
     param omega = 1.788 [1, 3]
     init x = 0.1
     init y = 0
     range x = [-3, 3]
     range y = [-8, 8]`,
    { view: { type: "trajectory" }, dt: 0.005 });

  add("forced-pendulum", "Forced damped pendulum", "Forced oscillators",
    "wadaR R/basins.R:73 and R/wada_detection.R:584",
    "Periodically driven pendulum x'' + gamma x' + sin x = F cos t; at gamma = 0.2, F = 1.66 its basins have the Wada property.",
    `x' = v
     v' = -gamma*v - sin(x) + F*cos(t)
     param gamma = 0.2 [0.05, 0.5]
     param F = 1.66 [0.5, 2.5]
     init x = 0
     init v = 0
     range x = [-3.1416, 3.1416]
     range v = [-4, 4]`,
    { view: { type: "strobe", period: PI2, transient: 10 }, n: 600, spread: 3, dt: 0.02, style: { pointSize: 1.3 } });

  add("driven-oscillator", "Driven damped linear oscillator", "Forced oscillators",
    "janos R/analysis_stroboscopic.R:75",
    "Every orbit converges to the unique periodic response; the stroboscopic points collapse onto one fixed point.",
    `x' = v
     v' = -w0^2*x - 2*zeta*v + F*cos(Om*t)
     param w0 = 1 [0.5, 2]
     param zeta = 0.1 [0.02, 0.5]
     param F = 0.5 [0, 1]
     param Om = 1.3 [0.5, 2]
     init x = 1
     init v = 0
     range x = [-2, 2]
     range v = [-2, 2]`,
    { view: { type: "strobe", period: PI2 / 1.3, transient: 0 }, n: 300, spread: 2, dt: 0.01, style: { fade: 0.02, pointSize: 2 } });

  // ========================================================== oscillators
  add("van-der-pol", "Van der Pol", "Oscillators and excitable media",
    "janos R/shiny_app.R:466 (van der Pol 1926)",
    "Self-sustained oscillation with nonlinear damping; relaxation oscillations for large mu.",
    `x' = y
     y' = mu*(1 - x^2)*y - x
     param mu = 2 [0.1, 8]
     init x = 2
     init y = 0
     range x = [-3, 3]
     range y = [-6, 6]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.01 });

  add("fitzhugh-nagumo", "FitzHugh-Nagumo", "Oscillators and excitable media",
    "janos R/shiny_app.R:453 (FitzHugh 1961; Nagumo et al. 1962)",
    "Slow-fast model of a neuron: excitable at low input I, oscillating above a Hopf threshold.",
    `v' = v - v^3/3 - w + I
     w' = eps*(v + a - b*w)
     param I = 0.5 [0, 2]
     param eps = 0.08 [0.01, 0.3]
     param a = 0.7 [0.3, 1]
     param b = 0.8 [0.3, 1]
     init v = -1
     init w = -0.5
     range v = [-2.5, 2.5]
     range w = [-1, 2]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02 });

  add("brusselator", "Brusselator", "Oscillators and excitable media",
    "janos R/shiny_app.R:437 (Prigogine and Lefever 1968)",
    "Autocatalytic reaction scheme: the equilibrium (A, B/A) undergoes a Hopf bifurcation at B = 1 + A^2.",
    `X' = A - (B + 1)*X + X^2*Y
     Y' = B*X - X^2*Y
     param A = 1 [0.5, 2]
     param B = 3 [0.3, 6]
     init X = 1
     init Y = 1
     range X = [0, 4.5]
     range Y = [0, 5.5]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.01 });

  add("selkov", "Sel'kov glycolysis", "Oscillators and excitable media",
    "janos vignettes/introduction.Rmd:209 (Sel'kov 1968)",
    "Glycolytic oscillations: a Hopf bifurcation in the flux b creates a limit cycle.",
    `x' = -x + a*y + x^2*y
     y' = b - a*y - x^2*y
     param a = 0.1 [0, 0.2]
     param b = 0.5 [0.1, 1.2]
     init x = 0.5
     init y = 0.5
     range x = [0, 3]
     range y = [0, 3]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.02 });

  add("pendulum", "Pendulum", "Oscillators and excitable media",
    "classical mechanics",
    "Frictionless pendulum x'' = -sin x: librations inside the separatrix through the saddles at x = +/-pi, rotations outside.",
    `x' = v
     v' = -sin(x) - c*v
     param c = 0 [0, 0.5]
     init x = 1
     init v = 0
     range x = [-7, 7]
     range v = [-3.2, 3.2]`,
    { view: { type: "phase", seeds: 14 }, dt: 0.01 });

  add("kuramoto-pair", "Two coupled phase oscillators", "Oscillators and excitable media",
    "nonautonomeR R/systems.R:1453 (kuramoto_pair)",
    "Phase difference psi = phi1 - phi2 obeys psi' = dw - 2K sin psi: phase locking when 2K > |dw|, phase slips otherwise.",
    `phi1' = w1 + K*sin(phi2 - phi1)
     phi2' = w2 + K*sin(phi1 - phi2)
     param w1 = 6.911503837897544 [5, 8]
     param w2 = 5.654866776461628 [5, 8]
     param K = 0.5 [0, 1.5]
     init phi1 = 0
     init phi2 = 1
     range phi1 = [0, 60]
     range phi2 = [0, 60]`,
    { view: { type: "timeseries", window: 30 }, dt: 0.01 });

  add("hindmarsh-rose", "Hindmarsh-Rose bursting neuron", "Oscillators and excitable media",
    "Hindmarsh and Rose (1984), Proc. R. Soc. Lond. B 221: 87",
    "Fast spiking (x, y) slaved to a slow adaptation current z; at I = 3.25 the number of spikes per burst varies chaotically (largest Lyapunov exponent about 0.013).",
    `x' = y - a*x^3 + b*x^2 - z + I
     y' = c - d*x^2 - y
     z' = r*(s*(x - xr) - z)
     param a = 1
     param b = 3
     param c = 1
     param d = 5
     param r = 0.006 [0.001, 0.02]
     param s = 4
     param xr = -1.6
     param I = 3.25 [2.5, 3.5]
     init x = 0.6241
     init y = 0.4467
     init z = 3.321
     range x = [-1.3, 1.9]
     range y = [-6.5, 0.8]
     range z = [2.9, 3.45]`,
    { view: { type: "flow" }, n: 2000, dt: 0.01, style: { ramp: "relab-fire" } });

  // ================================================================= maps
  add("logistic", "Logistic map", "Maps",
    "janos R/shiny_app.R:680 and kaRma R/demo_system.R:164 (May 1976)",
    "Period-doubling cascade to chaos; at r = 4 the Lyapunov exponent is ln 2.",
    `x[n+1] = r*x*(1 - x)
     param r = 3.9 [2.5, 4]
     init x = 0.2
     range x = [0, 1]`,
    { view: { type: "orbit", param: "r", from: 2.5, to: 4 }, style: { alpha: 0.25 } });

  add("logistic-cobweb", "Logistic map, cobweb", "Maps",
    "janos vignettes/qualitative-analysis.Rmd:346",
    "Graphical iteration of x[n+1] = r x (1 - x); click to restart from another x.",
    `x[n+1] = r*x*(1 - x)
     param r = 3.7 [2.5, 4]
     init x = 0.2
     range x = [0, 1]`,
    { view: { type: "cobweb", tail: 80 } });

  add("henon", "Henon map", "Maps",
    "janos R/shiny_app.R:688 and nonautonomeR R/systems.R:542 (Henon 1976)",
    "Stretch-and-fold map with a strange attractor; largest Lyapunov exponent 0.419 at a = 1.4, b = 0.3.",
    `x[n+1] = 1 - a*x^2 + y
     y[n+1] = b*x
     param a = 1.4 [0.8, 1.42]
     param b = 0.3 [0, 0.4]
     init x = 0.1
     init y = 0.1
     range x = [-1.5, 1.5]
     range y = [-0.45, 0.45]`,
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "relab-fire" } });

  add("lozi", "Lozi map", "Maps",
    "janos vignettes/chaotic-systems.Rmd:810 (Lozi 1978)",
    "Piecewise-linear analogue of the Henon map with a strange attractor.",
    `x[n+1] = 1 - a*abs(x) + y
     y[n+1] = b*x
     param a = 1.7 [1.2, 1.8]
     param b = 0.5 [0.2, 0.6]
     init x = 0.1
     init y = 0
     range x = [-1.4, 1.4]
     range y = [-0.7, 0.7]`,
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age" } });

  add("ikeda", "Ikeda map", "Maps",
    "janos vignettes/chaotic-systems.Rmd:842 (Ikeda 1979)",
    "Light in a nonlinear optical ring cavity; spiralling strange attractor.",
    `aux th = kappa - alpha/(1 + x^2 + y^2)
     x[n+1] = a + b*(x*cos(th) - y*sin(th))
     y[n+1] = b*(x*sin(th) + y*cos(th))
     param a = 1 [0.5, 1.5]
     param b = 0.9 [0.6, 0.95]
     param kappa = 0.4 [0, 1]
     param alpha = 6 [4, 8]
     init x = 0.1
     init y = 0.1
     range x = [-0.6, 2.2]
     range y = [-2.4, 1]`,
    { view: { type: "flow", life: [40, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "mako" } });

  add("standard-map", "Chirikov standard map", "Maps",
    "janos R/shiny_app.R:705 (Chirikov 1979)",
    "Area-preserving kicked rotor on the torus: KAM curves break up near K = 0.9716 and a chaotic sea spreads.",
    `p[n+1] = mod(p + K*sin(th), 2*pi)
     th[n+1] = mod(th + p + K*sin(th), 2*pi)
     param K = 1.2 [0, 3]
     init p = 0.5
     init th = 0.5
     range p = [0, 6.2832]
     range th = [0, 6.2832]`,
    { view: { type: "flow", life: [400, 2000], axes: ["th", "p"], dim3: false }, n: 1500, style: { fade: 0, pointSize: 1, alpha: 0.5, colorBy: "member", palette: "relab" } });

  add("zaslavsky", "Zaslavsky random map", "Maps",
    "nonautonomeR R/systems.R:741 (Namenson, Ott and Antonsen 1996)",
    "Dissipative kicked rotor with a random phase c_n drawn afresh at each step and shared by the ensemble: a snapshot attractor that changes shape every step.",
    `aux xn = mod(x + y*(1 - exp(-alpha))/alpha, 2*pi)
     x[n+1] = xn
     y[n+1] = kappa*sin(xn + 2*pi*ucommon(0)) + exp(-alpha)*y
     param alpha = 0.09 [0.02, 0.3]
     param kappa = 0.5 [0.1, 1]
     init x = 1
     init y = 0.1
     range x = [0, 6.2832]
     range y = [-4, 4]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, spread: 0, initMode: "box", style: { fade: 0.6, pointSize: 1.4, alpha: 0.8, colorBy: "solid", palette: "mono-amber" } });

  add("random-baker", "Random baker's map", "Maps",
    "nonautonomeR R/systems.R:801",
    "Baker's map with a random cut c_n shared by the ensemble: a fractal snapshot attractor with known generalised dimensions.",
    `aux c = ucommon(0)
     x[n+1] = ifelse(y <= c, lambda*x, 0.5 + lambda*x)
     y[n+1] = ifelse(y <= c, y/c, (y - c)/(1 - c))
     param lambda = 0.4 [0.1, 0.5]
     init x = 0.3
     init y = 0.6
     range x = [0, 1]
     range y = [0, 1]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, initMode: "box", style: { fade: 0.7, pointSize: 1.3, alpha: 0.8, colorBy: "solid", palette: "mono-ice" } });

  add("stark-circle", "Quasiperiodically forced circle map", "Maps",
    "nonautonomeR R/demo_stark_skew.R:47 (Stark 1999)",
    "Circle diffeomorphism driven by an irrational rotation theta. At a = 0.9 the fibre Lyapunov exponent is negative (about -0.24) and the attractor is the graph of a function of theta; at the package default a = 0.25 it is zero and orbits fill the torus.",
    `th[n+1] = mod(th + Om, 2*pi)
     x[n+1] = mod(x + nu + a*sin(x) + b*sin(th), 2*pi)
     param Om = 3.883222077450933 [0, 6.2832]
     param nu = 0.8676521529893011 [0, 6.2832]
     param a = 0.9 [0, 0.95]
     param b = 0.35 [0, 1]
     init th = 0.7
     init x = 0.2
     range th = [0, 6.2832]
     range x = [0, 6.2832]`,
    { view: { type: "flow", life: [300, 1200], dim3: false }, n: 2000, initMode: "box", style: { fade: 0.02, pointSize: 1.3, alpha: 0.6, colorBy: "var", colorVar: "x", ramp: "blackboard" } });

  add("kaplan-yorke", "Kaplan-Yorke map", "Maps",
    "nonautonomeR R/systems.R:915 (Kaplan and Yorke 1979)",
    "Doubling map driving a contracting coordinate: a strange attractor of Kaplan-Yorke dimension 1 + ln 2 / |ln alpha|. Floating-point doubling discards one bit per step and would reach x = 0 within 53 steps, so a noise of size 1e-9 re-injects the low-order bits.",
    `x[n+1] = mod(2*x + 0.000000001*nrand(), 1)
     y[n+1] = alpha*y + cos(4*pi*x)
     param alpha = 0.2 [0.05, 0.6]
     init x = 0.3678
     init y = 0.6677
     range x = [0, 1]
     range y = [-1.3, 1.3]`,
    { view: { type: "flow", life: [20, 60], dim3: false }, n: 3000, initMode: "box", style: { fade: 0.05, pointSize: 1.2, alpha: 0.5 } });

  add("de-jong", "Peter de Jong attractor", "Maps",
    "Pickover (1990), Computers, Pattern, Chaos and Beauty",
    "Trigonometric map whose attractors are used in generative art; many parameter sets are chaotic.",
    `x[n+1] = sin(a*y) - cos(b*x)
     y[n+1] = sin(c*x) - cos(d*y)
     param a = 1.4 [-3, 3]
     param b = -2.3 [-3, 3]
     param c = 2.4 [-3, 3]
     param d = -2.1 [-3, 3]
     init x = 0
     init y = 0
     range x = [-2.2, 2.2]
     range y = [-2.2, 2.2]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "age", ramp: "relab-fire" } });

  add("clifford", "Clifford attractor", "Maps",
    "Pickover (1990), Computers, Pattern, Chaos and Beauty",
    "Trigonometric map related to the de Jong map, with flowing filamentary attractors.",
    `x[n+1] = sin(a*y) + c*cos(a*x)
     y[n+1] = sin(b*x) + d*cos(b*y)
     param a = -1.4 [-3, 3]
     param b = 1.6 [-3, 3]
     param c = 1 [-3, 3]
     param d = 0.7 [-3, 3]
     init x = 0.1
     init y = 0.1
     range x = [-2.2, 2.2]
     range y = [-2.2, 2.2]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "var", colorVar: "y", ramp: "mako" } });

  // =========================================================== stochastic
  add("ornstein-uhlenbeck", "Ornstein-Uhlenbeck process", "Stochastic",
    "janos R/shiny_app.R:765 (ou)",
    "Mean-reverting diffusion with stationary variance sigma^2 / (2 theta) = 0.125.",
    `x' = -theta*(x - mu)
     noise x = sigma
     param theta = 1 [0.1, 3]
     param mu = 0 [-1, 1]
     param sigma = 0.5 [0, 1.5]
     init x = 2
     range x = [-1.5, 2.2]`,
    { view: { type: "density", window: 20 }, n: 3000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("double-well", "Noisy double well", "Stochastic",
    "janos R/shiny_app.R:793 and nonautonomeR vignettes/melancholia-states.Rmd:298",
    "Bistable potential V = x^4/4 - x^2/2 with additive noise: Kramers escapes between the wells at x = -1 and x = 1.",
    `x' = x - x^3
     noise x = sigma
     param sigma = 0.45 [0, 1]
     init x = -1
     range x = [-2, 2]`,
    { view: { type: "density", window: 200 }, n: 3000, dt: 0.02, style: { ramp: "magma" } });

  add("stochastic-resonance", "Stochastic resonance", "Stochastic",
    "janos vignettes/noise-in-dynamical-systems.Rmd:372 (Benzi, Sutera and Vulpiani 1981)",
    "A weak periodic tilt, too small to push the state over the barrier alone, becomes visible in the switching when the noise is tuned.",
    `x' = x - x^3 + A*cos(w*t)
     noise x = sigma
     param A = 0.12 [0, 0.4]
     param w = 0.1 [0.02, 0.5]
     param sigma = 0.35 [0, 1]
     init x = -1
     range x = [-2, 2]`,
    { view: { type: "timeseries", window: 400, members: 1 }, dt: 0.02 });

  add("verhulst-noise", "Noise-induced transition (stochastic Verhulst)", "Stochastic",
    "janos R/shiny_app.R:784 (Horsthemke and Lefever 1984)",
    "Logistic growth with multiplicative noise: the stationary density changes shape (a P-bifurcation) as sigma crosses sqrt(a).",
    `x' = a*x - x^2
     noise x = sigma*x
     param a = 1 [0.2, 2]
     param sigma = 0.8 [0, 2]
     init x = 1
     range x = [0, 3]`,
    { view: { type: "density", window: 60 }, n: 3000, dt: 0.005, keepPositive: true, style: { ramp: "mako" } });

  add("stochastic-lv", "Stochastic Lotka-Volterra", "Stochastic",
    "janos R/shiny_app.R:802 (sde_lv)",
    "Environmental noise pushes orbits off the conserved cycles of the Lotka-Volterra model.",
    `N' = r*N - a*N*P
     P' = e*a*N*P - m*P
     noise N = s*N
     noise P = s*P
     param r = 1 [0.2, 2]
     param a = 0.05 [0.01, 0.2]
     param e = 0.4 [0.1, 1]
     param m = 0.4 [0.05, 1.5]
     param s = 0.05 [0, 0.3]
     init N = 20
     init P = 10
     range N = [0, 80]
     range P = [0, 60]`,
    { view: { type: "density", decay: 0.97 }, n: 3000, dt: 0.01, style: { ramp: "inferno" } });

  add("coherence-resonance", "Coherence resonance", "Stochastic",
    "janos vignettes/noise-in-dynamical-systems.Rmd:390 (Pikovsky and Kurths 1997)",
    "Excitable FitzHugh-Nagumo unit below threshold: noise alone triggers spikes, most regular at an intermediate noise level.",
    `v' = (v - v^3/3 - w)/eps
     w' = v + a
     noise v = sigma
     param eps = 0.05 [0.01, 0.2]
     param a = 1.05 [0.9, 1.3]
     param sigma = 0.5 [0, 2]
     init v = -1
     init w = -0.6
     range v = [-2.5, 2.5]
     range w = [-1.2, 1.2]`,
    { view: { type: "timeseries", vars: ["v"], window: 40 }, dt: 0.002 });

  add("maier-stein", "Maier-Stein", "Stochastic",
    "nonautonomeR R/systems.R:992 (Maier and Stein 1993)",
    "Two attractors at (+/-1, 0) and a saddle at the origin; for alpha different from mu the drift is not a gradient and escape paths bend.",
    `x' = x - x^3 - alpha*x*y^2
     y' = -mu*(1 + x^2)*y
     noise x = sigma
     noise y = sigma
     param alpha = 3 [0, 6]
     param mu = 1 [0.2, 3]
     param sigma = 0.25 [0, 0.6]
     init x = -1
     init y = 0
     range x = [-1.8, 1.8]
     range y = [-1, 1]`,
    { view: { type: "density", decay: 0.96 }, n: 3000, dt: 0.01, style: { ramp: "magma" } });

  add("geometric-brownian", "Geometric Brownian motion", "Stochastic",
    "janos R/shiny_app.R:775 (gbm)",
    "Multiplicative noise with drift mu: the mean grows as exp(mu t) while the median grows as exp((mu - sigma^2/2) t).",
    `S' = mu*S
     noise S = sigma*S
     param mu = 0.08 [-0.2, 0.3]
     param sigma = 0.3 [0, 0.8]
     init S = 100
     range S = [0, 400]`,
    { view: { type: "timeseries", window: 10, members: 20 }, n: 20, dt: 0.002 });

  add("noisy-van-der-pol", "Noisy Van der Pol", "Stochastic",
    "janos vignettes/qualitative-analysis.Rmd:689",
    "Limit cycle blurred by velocity noise: the ensemble diffuses along the cycle and loses its phase.",
    `x' = y
     y' = mu*(1 - x^2)*y - x
     noise y = sigma
     param mu = 1.5 [0.2, 4]
     param sigma = 0.5 [0, 1.5]
     init x = 2
     init y = 0
     range x = [-3, 3]
     range y = [-5, 5]`,
    { view: { type: "density", decay: 0.9 }, n: 3000, spread: 0.01, dt: 0.01, style: { ramp: "relab-fire" } });

  // ================================================= tipping and nonautonomous
  add("r-tipping", "Rate-induced tipping", "Tipping and nonautonomous",
    "normal form of Ashwin, Wieczorek, Vitolo and Cox (2012), Phil. Trans. R. Soc. A 370: 1166",
    "In the frame y = x + lambda of x' = (x + lambda)^2 - 1, a ramp of lambda at rate r adds r to the drift. For r < 1 the state shifts to y = -sqrt(1 - r) and recovers when the ramp ends; for r > 1 no equilibrium exists and the state escapes (R-tipping) if the ramp outlasts the escape time, although every frozen lambda is safe.",
    `y' = y^2 - 1 + r*step(t - t0)*step(t0 + L - t)
     param r = 1.2 [0, 3]
     param t0 = 2 [0, 10]
     param L = 12 [1, 30]
     init y = -1
     range y = [-2, 3]`,
    { view: { type: "timeseries", window: 25, members: 1 }, dt: 0.005, overlay: { equations: true } });

  add("fold-normal-form", "Saddle-node (fold) normal form", "Tipping and nonautonomous",
    "janos vignettes/advanced-dynamics.Rmd:213",
    "x' = mu + x^2: two equilibria for mu < 0 that collide and vanish at mu = 0.",
    `x' = mu + x^2
     param mu = -1 [-1, 0.5]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "mu", from: -1, to: 0.3, speed: 0.00625 }, dt: 0.01 });

  add("cusp", "Cusp catastrophe", "Tipping and nonautonomous",
    "normal form of the cusp catastrophe (Thom 1972; Zeeman 1977)",
    "x' = r + a x - x^3: for a > 0 two stable branches coexist between the folds at r = +/-2 (a/3)^(3/2), so a slow sweep of r jumps and shows hysteresis.",
    `x' = r + a*x - x^3
     param r = 0 [-0.8, 0.8]
     param a = 1 [-0.5, 2]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "r", from: -0.8, to: 0.8, speed: 0.005 }, dt: 0.01, overlay: { equations: true } });

  add("pitchfork", "Pitchfork normal form", "Tipping and nonautonomous",
    "janos R/analysis_bifurcation_sweep.R:80",
    "x' = a x - x^3: the origin splits into two symmetric stable states at a = 0.",
    `x' = a*x - x^3
     param a = 1 [-1, 1.5]
     init x = 0.05
     range x = [-1.4, 1.4]`,
    { view: { type: "sweep", param: "a", from: -1, to: 1.5, speed: 0.00625 }, dt: 0.01 });

  add("hopf-normal-form", "Hopf normal form", "Tipping and nonautonomous",
    "janos vignettes/advanced-dynamics.Rmd:241",
    "Supercritical Hopf bifurcation: for mu > 0 a stable limit cycle of radius sqrt(mu) surrounds the origin.",
    `x' = mu*x - y - x*(x^2 + y^2)
     y' = x + mu*y - y*(x^2 + y^2)
     param mu = 0.5 [-0.5, 1]
     init x = 0.1
     init y = 0.1
     range x = [-1.3, 1.3]
     range y = [-1.3, 1.3]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02 });

  add("bogdanov-takens", "Bogdanov-Takens normal form", "Tipping and nonautonomous",
    "normal form (Kuznetsov 2004, Elements of Applied Bifurcation Theory, ch. 8)",
    "Codimension-two unfolding with fold, Hopf and homoclinic bifurcation curves meeting at the origin of (b1, b2).",
    `x' = y
     y' = b1 + b2*x + x^2 + x*y
     param b1 = -0.1 [-0.5, 0.2]
     param b2 = 0.2 [-0.5, 0.5]
     init x = 0
     init y = 0.1
     range x = [-1, 1]
     range y = [-1, 1]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.01 });

  add("stommel", "Stommel two-box ocean", "Tipping and nonautonomous",
    "tuRbulence R/stommel.R:73 (Stommel 1961)",
    "Thermohaline circulation with temperature T and salinity S: bistability between strong and weak overturning, with hysteresis in the freshwater forcing eta2.",
    `T' = eta1 - T*(1 + abs(T - S))
     S' = eta2 - S*(eta3 + abs(T - S))
     param eta1 = 3 [2, 4]
     param eta2 = 1 [0.5, 1.5]
     param eta3 = 0.3 [0.1, 0.6]
     init T = 2
     init S = 1
     range T = [0, 3.5]
     range S = [0, 3.5]`,
    { view: { type: "sweep", param: "eta2", var: "S", from: 0.5, to: 1.5, speed: 0.002 }, dt: 0.02 });

  add("lorenz84-forced", "Lorenz-84 under climate change", "Tipping and nonautonomous",
    "nonautonomeR R/systems.R:651 and vignettes/ergodicity.Rmd:174 (J\u00e1nosi, T\u00e9l and co-authors)",
    "Annual forcing F(t) = F0 + 2 sin(2 pi t / 73) with F0 ramping down: a cloud of 1500 members, all driven alike, traces a snapshot attractor that deforms as the climate changes.",
    `X' = -Y^2 - Z^2 - a*X + a*(F0 + AF*sin(2*pi*t/73))
     Y' = X*Y - b*X*Z - Y + G
     Z' = b*X*Y + X*Z - Z
     param a = 0.25 [0.1, 0.5]
     param b = 4 [2, 6]
     param G = 1 [0, 2]
     param F0 = 9.5 [6, 10]
     param AF = 2 [0, 3]
     init X = 1
     init Y = 0
     init Z = 0
     range X = [-1.5, 3.2]
     range Y = [-3, 3]
     range Z = [-3, 3]`,
    { view: { type: "flow", life: "inf", rotate: 0.1 }, n: 1500, spread: 1.2, initMode: "ball", dt: 0.02, perturbations: [{ kind: "ramp", param: "F0", rate: -0.000274, t0: 0, span: -2 }], style: { fade: 0.25, colorBy: "speed", ramp: "relab-fire" }, overlay: { readout: true } });

  add("lorenz-drift", "Lorenz with drifting rho", "Tipping and nonautonomous",
    "nonautonomeR vignettes/pullback-scenes.Rmd:145",
    "Rayleigh number rho(t) = 28 + 0.02 t: the attractor inflates slowly while the ensemble spreads over it.",
    `x' = sigma*(y - x)
     y' = x*(rho - z) - y
     z' = x*y - beta*z
     param sigma = 10 [5, 15]
     param rho = 28 [20, 40]
     param beta = 2.6666666666666665 [2, 3]
     init x = 1
     init y = 1
     init z = 20
     range x = [-26, 26]
     range y = [-34, 34]
     range z = [0, 62]`,
    { view: { type: "flow", life: "inf", rotate: 0.12 }, n: 1500, initMode: "ball", spread: 6, dt: 0.005, perturbations: [{ kind: "ramp", param: "rho", rate: 0.02, t0: 0, span: 12 }], style: { fade: 0.2 }, overlay: { readout: true } });

  add("duffing-drift", "Duffing with drifting forcing", "Tipping and nonautonomous",
    "nonautonomeR R/systems.R:963 (Janosi and Tel 2024)",
    "Forcing amplitude eps(t) = 0.4 + 0.00045 t: the stroboscopic cloud of an ensemble is a snapshot attractor whose shape follows the drift.",
    `x' = v
     v' = x - x^3 - 2*beta*v + eps*cos(omega*t)
     param beta = 0.2 [0.1, 0.4]
     param eps = 0.4 [0.2, 0.6]
     param omega = 1 [0.8, 1.2]
     init x = 0.5
     init v = 0
     range x = [-2, 2]
     range v = [-1.6, 1.6]`,
    { view: { type: "strobe", period: PI2, transient: 3 }, n: 1500, spread: 1.5, initMode: "ball", dt: 2 * Math.PI / 128, perturbations: [{ kind: "ramp", param: "eps", rate: 0.00045, t0: 0, span: 0.3 }], style: { fade: 0.35, pointSize: 1.6, alpha: 0.8 }, overlay: { readout: true } });

  add("tilted-well", "Double well with a slow tilt", "Tipping and nonautonomous",
    "nonautonomeR vignettes/melancholia-states.Rmd:524",
    "The tilt lambda is ramped slowly: the occupied well loses stability at a fold, and the ensemble switches, earlier when noise is present.",
    `x' = x - x^3 + lambda
     noise x = sigma
     param lambda = -0.6 [-0.6, 0.6]
     param sigma = 0.15 [0, 0.5]
     init x = -1.2
     range x = [-1.8, 1.8]`,
    { view: { type: "density", window: 120 }, n: 3000, dt: 0.02, perturbations: [{ kind: "ramp", param: "lambda", rate: 0.01, t0: 0, span: 1.2 }], style: { ramp: "magma" }, overlay: { readout: true } });

  // ============================================================ epidemics
  add("sir", "SIR epidemic", "Epidemics",
    "janos R/shiny_app.R:477 (Kermack and McKendrick 1927)",
    "Closed epidemic with basic reproduction number R0 = beta / gamma = 3.",
    `S' = -beta*S*I/N0
     I' = beta*S*I/N0 - gamma*I
     R' = gamma*I
     param beta = 0.3 [0.05, 1]
     param gamma = 0.1 [0.02, 0.5]
     param N0 = 1000 [100, 2000]
     init S = 999
     init I = 1
     init R = 0
     range S = [0, 1000]
     range I = [0, 1000]
     range R = [0, 1000]`,
    { view: { type: "timeseries", window: 160 }, dt: 0.05 });

  add("seasonal-sir", "Seasonally forced SIR", "Epidemics",
    "janos vignettes/chaotic-systems.Rmd:1204 (Olsen and Schaffer 1990)",
    "Measles-like epidemic with seasonal transmission: irregular outbreaks from a smooth forcing.",
    `S' = mu - beta0*(1 + beta1*cos(2*pi*t))*S*I - mu*S
     I' = beta0*(1 + beta1*cos(2*pi*t))*S*I - (gam + mu)*I + eps
     param mu = 0.02 [0.01, 0.04]
     param beta0 = 1800 [1000, 2500]
     param beta1 = 0.08 [0, 0.3]
     param gam = 100 [50, 150]
     param eps = 0.000001 [0, 0.00001]
     init S = 0.065
     init I = 0.0002
     range S = [0.045, 0.085]
     range I = [0, 0.002]`,
    { view: { type: "timeseries", window: 20, vars: ["I"] }, dt: 0.0005 });

  DF.CATALOGUE = M;
  DF.catalogueGroups = function () {
    const seen = [];
    M.forEach(function (m) { if (seen.indexOf(m.group) < 0) seen.push(m.group); });
    return seen;
  };
  // Scene for a catalogue entry: the entry's defaults with its system text.
  DF.sceneFor = function (id, overrides) {
    const m = M.find(function (e) { return e.id === id; });
    if (!m) throw new Error("No model '" + id + "' in the catalogue");
    const s = JSON.parse(JSON.stringify(m.scene));
    s.name = m.name; s.model = m.id; s.system = m.system;
    s.overlay = Object.assign({ title: m.name }, s.overlay || {});
    return Object.assign(s, overrides || {});
  };
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

// ---- src/component/relab-flow.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* <relab-flow>: a scene as an HTML element. One script tag and one element
   place a live figure in a web page, a Quarto or R Markdown document, a
   reveal.js slide or a pkgdown article.

     <script src="relabflow.js"></script>
     <relab-flow scene='{"system": "...", "view": {"type": "trajectory"}}'></relab-flow>
     <relab-flow model="lorenz" theme="blackboard" controls></relab-flow>
     <relab-flow src="figures/fold.json" paused></relab-flow>

   Attributes: scene (JSON), src (URL of a scene file), model (catalogue
   identifier), theme, view, title, controls (show play, restart and
   fullscreen buttons), paused (do not start automatically), static (draw a
   still after a number of frames given by the attribute, default 240).
   Changing scene, src, model, theme, view or title loads the scene again. The
   element exposes .player, .play(), .pause(), .restart(), .setParam(name, v)
   and .scene. A figure starts only when it is visible and shows a still frame
   when the reader asks for reduced motion. */
(function (DF) {
  "use strict";
  if (typeof HTMLElement === "undefined" || typeof customElements === "undefined") return;

  const CSS = ":host{display:block;position:relative;min-height:200px;height:360px;contain:content}" +
    ".stage{position:absolute;inset:0}" +
    ".ctl{position:absolute;right:10px;bottom:10px;display:flex;gap:6px;opacity:0;transition:opacity .2s;z-index:5}" +
    ":host(:hover) .ctl,.ctl:focus-within{opacity:1}" +
    "button{all:unset;cursor:pointer;width:30px;height:30px;border-radius:50%;display:grid;place-items:center;background:rgba(23,12,58,.72);color:#fff;font:14px system-ui}" +
    "button:hover{background:#EE6A24}button:focus-visible{outline:2px solid #FB9E07}" +
    ".err{position:absolute;inset:0;display:grid;place-items:center;padding:16px;font:13px system-ui;color:#CF4446;background:#fff8f5;text-align:center}";

  class RElabFlowElement extends HTMLElement {
    static get observedAttributes() { return ["scene", "src", "model", "theme", "view", "title"]; }
    /* The shadow tree and the player outlive a move of the element: moving it
       to another parent (a slide framework does this) pauses the figure and
       resumes it where it was, and does not build a second shadow root. */
    connectedCallback() {
      if (this._root) { if (this.player && this._wasRunning) this.player.play(); return; }
      this._root = this.shadowRoot || this.attachShadow({ mode: "open" });
      this._root.innerHTML = "";
      const style = document.createElement("style"); style.textContent = CSS;
      this._stage = document.createElement("div"); this._stage.className = "stage";
      this._root.append(style, this._stage);
      if (this.hasAttribute("controls")) this.buildControls();
      this.loadScene();
    }
    disconnectedCallback() { if (this.player) { this._wasRunning = this.player.running; this.player.pause(); } }
    attributeChangedCallback(name, oldV, newV) { if (this._root && oldV !== newV && (this.player || this._failed)) this.loadScene(); }

    buildControls() {
      const c = document.createElement("div"); c.className = "ctl";
      const mk = (label, title, fn) => { const b = document.createElement("button"); b.textContent = label; b.title = title; b.setAttribute("aria-label", title); b.addEventListener("click", fn); c.appendChild(b); return b; };
      this._playBtn = mk("\u275a\u275a", "Pause", () => this.player && this.player.toggle());
      mk("\u21bb", "Restart", () => this.restart());
      mk("\u26f6", "Full screen", () => { if (document.fullscreenElement) document.exitFullscreen(); else this.requestFullscreen && this.requestFullscreen(); });
      this._root.appendChild(c);
    }

    async loadScene() {
      let scene;
      try {
        if (this.hasAttribute("scene")) scene = JSON.parse(this.getAttribute("scene"));
        else if (this.hasAttribute("src")) {
          const src = this.getAttribute("src");
          let res;
          try { res = await fetch(src); } catch (e) { throw new Error("Could not load " + src + ". A page opened from disk cannot fetch files: serve the folder, or give the scene attribute instead"); }
          if (!res.ok) throw new Error("Could not load " + src + " (HTTP " + res.status + ")");
          try { scene = await res.json(); } catch (e) { throw new Error(src + " is not a scene file (JSON)"); }
        }
        else if (this.hasAttribute("model")) {
          if (!DF.sceneFor) throw new Error("The model catalogue is not loaded");
          scene = DF.sceneFor(this.getAttribute("model"));
        } else throw new Error("Give a scene, src or model attribute");
      } catch (e) { this._failed = true; this.showError(e); return; }
      if (this.hasAttribute("theme")) scene.style = Object.assign({}, scene.style, { theme: this.getAttribute("theme") });
      if (this.hasAttribute("view")) scene.view = Object.assign({}, scene.view, { type: this.getAttribute("view") });
      if (this.hasAttribute("title")) scene.overlay = Object.assign({}, scene.overlay, { title: this.getAttribute("title") });
      if (this.player) { this.player.dispose(); this.player = null; }
      this._stage.innerHTML = "";
      try {
        this.player = new DF.Player(this._stage, scene);
        if (this.player.error) throw this.player.error;
      } catch (e) { if (this.player) { this.player.dispose(); this.player = null; } this._failed = true; this.showError(e); return; }
      this._failed = false;
      this.player.on("state", (s) => { if (this._playBtn) { this._playBtn.textContent = s === "play" ? "\u275a\u275a" : "\u25b6"; this._playBtn.title = s === "play" ? "Pause" : "Play"; } });
      this.player.on("error", (e) => this.showError(e));
      const reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      if (this.hasAttribute("static") || reduced) this.player.advance(parseInt(this.getAttribute("static"), 10) || 240);
      else if (!this.hasAttribute("paused")) this.player.play();
      this.dispatchEvent(new CustomEvent("ready", { detail: this.player }));
    }
    showError(e) {
      const d = document.createElement("div"); d.className = "err";
      d.textContent = "RElabFlow: " + (e && e.message ? e.message : String(e));
      this._stage.innerHTML = ""; this._stage.appendChild(d);
    }
    get scene() { return this.player ? this.player.getScene() : null; }
    play() { if (this.player) this.player.play(); }
    pause() { if (this.player) this.player.pause(); }
    restart(seed) { if (this.player) { this.player.restart(seed); if (!this.hasAttribute("paused")) this.player.play(); } }
    setParam(name, value) { if (this.player) this.player.setParam(name, value); }
  }

  if (!customElements.get("relab-flow")) customElements.define("relab-flow", RElabFlowElement);
  DF.RElabFlowElement = RElabFlowElement;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});

globalThis.RElabFlow.VERSION = "0.1.0";
