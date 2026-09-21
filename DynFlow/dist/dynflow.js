// DynFlow 0.1.0. Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab.
// SPDX-License-Identifier: GPL-3.0-or-later. https://www.gnu.org/licenses/gpl-3.0.html
// Built from 11 source files by tools/build.mjs; edit the sources, not this file.
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
})(globalThis.DynFlow = globalThis.DynFlow || {});

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
      maxLag: function (p) { let m = 0; lagFns.forEach(function (fn) { m = Math.max(m, fn(p, 0)); }); return m; },
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
          const a = wrap(node.a, 6), b = wrap(node.b, 6);
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
  DF.parseExpression = parse;
  DF.compileSystem = compileSystem;
  DF.systemLatex = systemLatex;
  DF.texName = texName;
})(globalThis.DynFlow = globalThis.DynFlow || {});

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
    this.maxHistory = opts.maxHistory || 20000;
    this.reset();
  }

  Simulator.prototype.reset = function (seed) {
    if (seed !== undefined) this.seed = seed;
    const n = this.n, dim = this.dim;
    this.rng = new DF.RNG(this.seed);
    const rng = this.rng;
    this.R = { u: function () { return rng.uniform(); }, n: function () { return rng.normal(); }, U: new Float64Array(8), N: new Float64Array(8) };
    this.t = 0; this.steps = 0;
    this.X = new Float64Array(n * dim);
    this.alive = new Uint8Array(n).fill(1);
    for (let k = 0; k < n; k++) this.initMember(k);
    this.rk4 = makeRK4(dim);
    this.dx = new Float64Array(dim); this.gx = new Float64Array(dim); this.xk = new Float64Array(dim);
    this.eta = new Float64Array(this.perturbations.length);         // common OU states
    this.etaK = new Float64Array(this.perturbations.length * n);    // per-member OU states
    this.nextPulse = this.perturbations.map(function (q) { return q.kind === "pulse" ? (q.t0 || q.period || 1) : Infinity; });
    this.hist = null;
    this.updateParams();
    if (this.sys.kind === "dde") {
      const lag = this.sys.maxLag(this.p);
      const cap = Math.min(this.maxHistory, Math.ceil(lag / this.h) + 8);
      this.hist = [];
      for (let k = 0; k < n; k++) this.hist.push(this.newHistory(k, cap));
    }
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
  };

  // Advance the whole ensemble by one step of length h.
  Simulator.prototype.step = function () {
    const sys = this.sys, n = this.n, dim = this.dim, h = this.h, t = this.t, r = this.rng;
    const X = this.X, p = this.p, pert = this.perturbations;
    const sqh = Math.sqrt(h);
    // Common noise increments, drawn once per step for the whole ensemble.
    if (sys.usesRandom) for (let i = 0; i < 8; i++) { this.R.U[i] = r.uniform(); this.R.N[i] = r.normal(); }
    const R = this.R;
    const commonDW = pert.map(function (q) { return q.common && q.enabled !== false ? r.normal() : 0; });
    const commonJump = pert.map(function (q) { return q.kind === "jumps" && q.common && q.enabled !== false ? r.poisson(q.rate * h) : 0; });
    const commonStable = pert.map(function (q) { return q.kind === "levy" && q.common && q.enabled !== false ? r.stable(q.alpha || 1.5) : 0; });

    for (let k = 0; k < n; k++) {
      if (!this.alive[k]) continue;
      const x = X.subarray(k * dim, k * dim + dim);
      const Hk = this.hist ? this.hist[k] : null;
      const H = Hk ? function (i, s) { return Hk.at(i, s); } : undefined;
      if (sys.time === "discrete") {
        sys.f(t, x, p, this.dx, H, R);
        x.set(this.dx);
      } else if (sys.kind === "sde") {
        sys.f(t, x, p, this.dx, H, R);
        sys.g(t, x, p, this.gx, H, R);
        for (let i = 0; i < dim; i++) x[i] += this.dx[i] * h + (sys.noiseMask[i] ? this.gx[i] * sqh * r.normal() : 0);
      } else {
        this.rk4(sys.f, t, x, p, h, H, R);
      }
      // State perturbations; in discrete time h is 1 and dW has variance 1.
      for (let j = 0; j < pert.length; j++) {
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
      if (Hk && ok) { sys.f(t + h, x, p, this.dx, H); Hk.push(t + h, x, this.dx); }
    }
    this.t = t + (sys.time === "discrete" ? 1 : h);
    this.steps++;
    // Pulses at fixed times act on every member at once.
    for (let j = 0; j < pert.length; j++) {
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
    for (let j = 0; j < pert.length; j++) {
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
})(globalThis.DynFlow = globalThis.DynFlow || {});

// ---- src/core/analysis.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Local analysis: Jacobians by central differences, eigenvalues of a real
   matrix, equilibria (flows) and fixed points (maps) by damped Newton
   iteration from many seeds, and their stability. */
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

  DF.residual = residual;
  DF.jacobian = jacobian;
  DF.solveLinear = solve;
  DF.eigenvalues = eigenvalues;
  DF.findEquilibria = findEquilibria;
})(globalThis.DynFlow = globalThis.DynFlow || {});

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
})(globalThis.DynFlow = globalThis.DynFlow || {});

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

// ---- src/render/player.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Player: renders one scene into a host element. The same Player runs in
   the studio, in the <dyn-flow> element and in exported standalone pages.

   A scene is plain JSON and carries its own system text, so it is complete
   without the catalogue:

     { version, name, system, params: {a: 1}, init: {x: 0.1},
       dt, stepsPerFrame, seed, n, initMode, spread,
       perturbations: [{kind, ...}],
       view:  { type, axes: ['x','y'], ranges: {x: [lo, hi]}, ... },
       style: { theme, palette, ramp, colorBy, fade, lineWidth, alpha, pointSize },
       overlay: { title, subtitle, caption, equations, readout, legend } } */
(function (DF) {
  "use strict";

  const VIEW_DEFAULTS = {
    flow: { n: 1500, stepsPerFrame: 2, life: [80, 320], spawn: "box", fade: 0.06, lineWidth: 1.1, alpha: 0.55, colorBy: "dominant" },
    trajectory: { n: 1, stepsPerFrame: 4, tail: 2500, fade: 0, lineWidth: 1.3, alpha: 0.95, colorBy: "time", rotate: 0.25 },
    timeseries: { n: 1, stepsPerFrame: 2, fade: 0, lineWidth: 1.6, alpha: 0.95, colorBy: "solid" },
    phase: { n: 6, stepsPerFrame: 2, fade: 0, lineWidth: 1.5, alpha: 0.9, colorBy: "solid" },
    sweep: { n: 1, stepsPerFrame: 8, fade: 0, lineWidth: 1.6, alpha: 0.9, colorBy: "solid" },
    orbit: { n: 1, stepsPerFrame: 1, fade: 0, lineWidth: 1, alpha: 0.35, colorBy: "solid" },
    density: { n: 3000, stepsPerFrame: 2, fade: 0, lineWidth: 1, alpha: 1, colorBy: "solid" },
    strobe: { n: 200, stepsPerFrame: 20, fade: 0, lineWidth: 1, alpha: 0.6, colorBy: "solid", pointSize: 1.4 },
    cobweb: { n: 1, stepsPerFrame: 1, fade: 0, lineWidth: 1.3, alpha: 0.9, colorBy: "solid" }
  };

  function clone(o) { return JSON.parse(JSON.stringify(o === undefined ? null : o)); }

  // Fill defaults; never throws on a partial scene.
  function normalizeScene(input) {
    const s = clone(input || {}) || {};
    s.version = 1;
    s.name = s.name || "Untitled scene";
    s.system = s.system || "x' = -y\ny' = x\ninit x = 1";
    s.params = s.params || {};
    s.init = s.init || {};
    s.view = s.view || {};
    s.view.type = s.view.type && DF.VIEWS[s.view.type] ? s.view.type : "flow";
    const d = VIEW_DEFAULTS[s.view.type];
    s.n = s.n || d.n;
    s.stepsPerFrame = s.stepsPerFrame || d.stepsPerFrame;
    s.seed = s.seed === undefined ? 1 : s.seed;
    s.initMode = s.initMode || "point";
    s.spread = s.spread === undefined ? 0.05 : s.spread;
    s.perturbations = s.perturbations || [];
    if (s.view.life === undefined && s.view.type === "flow") s.view.life = d.life;
    if (s.view.spawn === undefined && s.view.type === "flow") s.view.spawn = d.spawn;
    if (s.view.rotate === undefined && d.rotate !== undefined) s.view.rotate = d.rotate;
    if (s.view.tail === undefined && d.tail !== undefined) s.view.tail = d.tail;
    s.style = Object.assign({ theme: "relab-night", palette: "relab", ramp: "relab-fire", colorBy: d.colorBy, fade: d.fade, lineWidth: d.lineWidth, alpha: d.alpha, pointSize: d.pointSize || 1.6, renderScale: 1 }, s.style || {});
    s.overlay = Object.assign({ title: "", subtitle: "", caption: "", equations: false, readout: false, legend: true, position: "top-left" }, s.overlay || {});
    return s;
  }

  function Player(host, scene, opts) {
    this.opts = opts || {};
    this.host = host;
    this.listeners = {};
    this.running = false; this.frameCount = 0; this.fps = 0;
    this.speedMax = 1e-9;
    this.buildDom();
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

  Player.prototype.bindPointer = function () {
    const self = this, el = this.cv.top;
    let down = null;
    const pos = function (ev) { const r = el.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    el.addEventListener("pointerdown", function (ev) {
      const p = pos(ev); down = { x: p[0], y: p[1], azim: self.cam ? self.cam.azim : 0, elev: self.cam ? self.cam.elev : 0, moved: false };
      el.setPointerCapture(ev.pointerId);
      if (!(self.cam && self.cam.is3D()) && self.view.pointer && self.view.pointer(self, "down", p[0], p[1], ev)) self.emit("interact", {});
    });
    el.addEventListener("pointermove", function (ev) {
      if (!down) return;
      const p = pos(ev), dx = p[0] - down.x, dy = p[1] - down.y;
      if (Math.abs(dx) + Math.abs(dy) > 3) down.moved = true;
      if (self.cam && self.cam.is3D()) { self.dragging = true; self.cam.azim = down.azim + dx * 0.01; self.cam.elev = Math.max(-1.5, Math.min(1.5, down.elev + dy * 0.01)); self.clearTrail(); }
      else if (self.view.pointer) self.view.pointer(self, "drag", p[0], p[1], ev);
    });
    const up = function () { down = null; self.dragging = false; };
    el.addEventListener("pointerup", up); el.addEventListener("pointercancel", up);
    el.addEventListener("wheel", function (ev) {
      if (!self.cam || !self.scene.view.zoomable) return;
      ev.preventDefault(); self.cam.zoom = Math.max(0.2, Math.min(20, self.cam.zoom * Math.exp(-ev.deltaY * 0.001))); self.redrawStatic(); self.clearTrail();
    }, { passive: false });
  };

  // Load (or reload) a scene: compile, build the simulator, the camera and the view.
  Player.prototype.load = function (scene) {
    const s = normalizeScene(scene);
    let sys;
    try { sys = DF.compileSystem(s.system); }
    catch (e) { this.error = e; this.emit("error", e); if (!this.sys) return; sys = this.sys; }
    this.error = null;
    this.scene = s; this.sys = sys;
    this.buildSim();
    this.emit("scene", s);
  };

  Player.prototype.buildSim = function () {
    const s = this.scene, sys = this.sys;
    const params = sys.params.map(function (q) { return s.params[q.name] !== undefined ? +s.params[q.name] : q.value; });
    const init = sys.vars.map(function (v, i) { return s.init[v] !== undefined ? +s.init[v] : sys.init[i]; });
    const dt = s.dt || (sys.time === "discrete" ? 1 : 0.01);
    s.dt = dt;
    // Full ranges for every variable: scene, then declared, then a probe run.
    const need = sys.vars.filter(function (v) { return !(s.view.ranges && s.view.ranges[v]) && !sys.ranges[v]; });
    const probe = need.length ? DF.autoRanges(sys, params, init, { dt: dt, steps: sys.time === "discrete" ? 3000 : Math.min(20000, Math.max(3000, Math.round(60 / dt))), perturbations: s.perturbations }) : {};
    this.fullRanges = sys.vars.map(function (v) { return (s.view.ranges && s.view.ranges[v]) || sys.ranges[v] || probe[v]; });
    this.sim = new DF.Simulator(sys, {
      n: s.n, dt: dt, seed: s.seed, params: params, init: init, initMode: s.initMode, spread: s.spread,
      perturbations: s.perturbations, keepPositive: !!s.keepPositive, box: this.fullRanges
    });
    this.tmpX = new Float64Array(sys.vars.length); this.tmpDx = new Float64Array(sys.vars.length);
    let axes = (s.view.axes || []).map(function (v) { return sys.vars.indexOf(v); }).filter(function (i) { return i >= 0; });
    if (axes.length < 2) axes = sys.vars.length >= 3 && (s.view.type === "trajectory" || s.view.type === "flow") && s.view.dim3 !== false ? [0, 1, 2] : sys.vars.length >= 2 ? [0, 1] : [0, 0];
    if (sys.vars.length === 1 && s.view.type !== "cobweb" && s.view.type !== "density" && s.view.type !== "timeseries" && s.view.type !== "sweep" && s.view.type !== "orbit") {
      // A scalar system is shown as x against time through the timeseries view.
      s.view.type = "timeseries";
    }
    const self = this;
    if (s.view.projection === "simplex" && sys.vars.length >= 3 && axes.length < 3) axes = [0, 1, 2];
    this.cam = new DF.Camera(axes, axes.map(function (i) { return self.fullRanges[i]; }), { azim: s.view.azim, elev: s.view.elev, zoom: s.view.zoom, pad: s.view.pad, upAxis: s.view.upAxis, simplex: s.view.projection === "simplex" });
    this.theme = DF.THEMES[s.style.theme] || DF.THEMES["relab-night"];
    this.palette = (DF.PALETTES[s.style.palette] || DF.PALETTES.relab).colors;
    this.view = Object.create(DF.VIEWS[s.view.type]);
    this.userParam = false;
    this.frameCount = 0;
    this.speedMax = 1e-9;
    this.resize(true);
    // A new scene starts on empty layers; otherwise trails of the previous
    // scene remain under views that never draw on the trail layer.
    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);
    this.view.init(this);
    this.redrawStatic();
    this.renderOverlay();
    if (!this.running) this.drawOnce();
  };

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
    if (changed && !skipStatic && this.view) { this.clearTrail(); this.redrawStatic(); if (!this.running) this.drawOnce(); }
  };

  Player.prototype.clearTrail = function () { this.ctx.trail.clearRect(0, 0, this.w, this.h); };
  Player.prototype.fade = function (a) {
    if (!(a > 0)) return;
    const c = this.ctx.trail;
    c.globalCompositeOperation = "destination-out"; c.fillStyle = "rgba(0,0,0," + a + ")";
    c.fillRect(0, 0, this.w, this.h); c.globalCompositeOperation = "source-over";
  };
  Player.prototype.redrawStatic = function () {
    this.ctx.base.clearRect(0, 0, this.w, this.h);
    if (this.scene.view.showAxes && this.cam && !this.cam.is3D() && !this.view.axes) DF.drawAxes(this.ctx.base, this.cam, this.scene.style.theme, [this.sys.vars[this.cam.axes[0]], this.sys.vars[this.cam.axes[1]]], {});
    if (this.scene.view.showAxes && this.cam && this.cam.is3D()) DF.drawBox3D(this.ctx.base, this.cam, this.scene.style.theme);
    if (this.cam && this.cam.simplex && this.scene.view.showAxes !== false) DF.drawSimplex(this.ctx.base, this.cam, this.scene.style.theme, this.cam.axes.map(function (i) { return this.sys.vars[i]; }, this));
    if (this.view.drawStatic) this.view.drawStatic(this);
  };

  // ------------------------------------------------------------ overlay
  Player.prototype.renderOverlay = function () {
    const o = this.scene.overlay, el = this.overlayEl, th = this.theme, doc = el.ownerDocument;
    el.innerHTML = "";
    if (o.position === "none") return;
    const box = doc.createElement("div");
    box.style.cssText = "position:absolute;left:18px;top:14px;right:18px;color:" + th.ink + ";";
    if (o.title) { const t = doc.createElement("div"); t.textContent = o.title; t.style.cssText = "font-weight:600;font-size:clamp(15px,2.4vw,26px);line-height:1.15;"; box.appendChild(t); }
    if (o.subtitle) { const t = doc.createElement("div"); t.textContent = o.subtitle; t.style.cssText = "font-weight:300;font-size:clamp(12px,1.5vw,16px);color:" + th.muted + ";margin-top:2px;"; box.appendChild(t); }
    el.appendChild(box);
    if (o.equations) {
      const eq = doc.createElement("div");
      eq.style.cssText = "position:absolute;right:18px;top:" + (this.view.axes ? 70 : 16) + "px;color:" + th.ink + ";font-size:clamp(11px,1.4vw,15px);text-align:right;opacity:0.92;";
      const lines = DF.systemLatex(this.scene.system).lines;
      const k = typeof katex !== "undefined" ? katex : (typeof window !== "undefined" ? window.katex : undefined);
      if (k) {
        lines.forEach(function (L) {
          const d = doc.createElement("div"); d.style.margin = "2px 0";
          try { k.render(L, d, { throwOnError: false, displayMode: false }); } catch (e) { d.textContent = L; }
          eq.appendChild(d);
        });
      } else {
        // Without KaTeX (a standalone page), the equations are shown as written.
        String(this.scene.system).split(/\r?\n/).map(function (l) { return l.replace(/#.*$/, "").trim(); })
          .filter(function (l) { return /'\s*=|\[n\+1\]|dt\s*=|^noise/.test(l); })
          .forEach(function (L) { const d = doc.createElement("div"); d.textContent = L; d.style.cssText = "margin:2px 0;font-family:'TeX Gyre Pagella',Palatino,serif;font-style:italic;"; eq.appendChild(d); });
      }
      el.appendChild(eq);
    }
    if (o.caption) {
      const c = doc.createElement("div"); c.textContent = o.caption;
      c.style.cssText = "position:absolute;left:18px;bottom:" + (this.view.axes ? 6 : 12) + "px;max-width:min(560px,70%);font-size:11px;line-height:1.45;color:" + th.muted + ";";
      el.appendChild(c);
    }
    this.readoutEl = null;
    if (o.readout) {
      const r = doc.createElement("div");
      r.style.cssText = "position:absolute;right:18px;bottom:" + (this.view.axes ? 6 : 12) + "px;font:11px ui-monospace,monospace;color:" + th.muted + ";text-align:right;";
      el.appendChild(r); this.readoutEl = r;
    }
    if (o.legend && this.view.legend) {
      const items = this.view.legend(this);
      if (items && items.length) {
        const lg = doc.createElement("div");
        const ins = this.inset || { r: 0, t: 0 };
        const panelBg = th.dark ? "rgba(11,6,32,0.55)" : "rgba(255,255,255,0.8)";
        // With equations in the top-right corner the legend goes to the bottom of the plot.
        const atBottom = !!o.equations;
        const bottomPx = this.view.axes ? ins.b + 8 : (o.readout ? 30 : 12);
        lg.style.cssText = "position:absolute;right:" + (ins.r + 10) + "px;" + (atBottom ? "bottom:" + bottomPx + "px;" : "top:" + (this.view.axes ? ins.t + 8 : (o.title ? 64 : 12)) + "px;") +
          "max-width:48%;display:flex;flex-direction:column;align-items:flex-start;gap:1px;padding:5px 9px;border-radius:6px;background:" + panelBg + ";font-size:11px;color:" + th.muted + ";";
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
    }
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
    const spf = this.scene.stepsPerFrame;
    this.scene.stepsPerFrame = 0;
    try { this.view.frame(this); } catch (e) { /* first frame of some views needs steps */ }
    this.scene.stepsPerFrame = spf;
  };
  Player.prototype.tick = function () {
    try { this.view.frame(this); }
    catch (e) { this.pause(); this.error = e; this.emit("error", e); return; }
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
      last = ts;
      self.tick();
      self.raf = requestAnimationFrame(loop);
    };
    this.raf = requestAnimationFrame(loop);
  };
  Player.prototype.play = function () { if (this.running) return; this.running = true; this.emit("state", "play"); this.kick(); };
  Player.prototype.pause = function () { this.running = false; if (this.raf) cancelAnimationFrame(this.raf); this.raf = 0; this.emit("state", "pause"); };
  Player.prototype.toggle = function () { if (this.running) this.pause(); else this.play(); };
  Player.prototype.step = function () { this.tick(); };
  Player.prototype.restart = function (seed) {
    if (seed !== undefined) this.scene.seed = seed;
    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);
    this.buildSim();
    this.emit("restart", this.scene.seed);
  };
  // Advance n frames without waiting for the screen (static renders, tests).
  Player.prototype.advance = function (n) { for (let i = 0; i < n; i++) this.tick(); };

  // Change one parameter live, without restarting.
  Player.prototype.setParam = function (name, value) {
    const i = this.sys.params.findIndex(function (q) { return q.name === name; });
    if (i < 0) return;
    this.scene.params[name] = value;
    this.sim.base[i] = value; this.sim.updateParams();
    if (this.view.pi === i) this.userParam = true;
    if (this.view.branchPts) this.view.branchPts = null;
    if (this.scene.view.type === "phase" || this.scene.view.type === "cobweb") this.redrawStatic();
    if (this.scene.view.type === "orbit") { this.restart(); }
  };
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
  Player.prototype.paintText = function (g) {
    const o = this.scene.overlay, th = this.theme;
    if (o.position === "none") return;
    let y = 14;
    g.textBaseline = "top"; g.textAlign = "left";
    if (o.title) { g.fillStyle = th.ink; g.font = "600 " + Math.round(Math.max(15, Math.min(26, this.w * 0.024))) + "px Jost, system-ui, sans-serif"; g.fillText(o.title, 18, y); y += Math.max(18, Math.min(30, this.w * 0.028)); }
    if (o.subtitle) { g.fillStyle = th.muted; g.font = "300 " + Math.round(Math.max(12, Math.min(16, this.w * 0.015))) + "px Jost, system-ui, sans-serif"; g.fillText(o.subtitle, 18, y); }
    const bottom = this.view.axes ? 8 : 14;
    if (o.caption) {
      g.fillStyle = th.muted; g.font = "11px Jost, system-ui, sans-serif"; g.textBaseline = "bottom";
      const words = o.caption.split(" "), maxW = Math.min(560, this.w * 0.7), lines = [];
      let line = "";
      words.forEach(function (w) { const tt = line ? line + " " + w : w; if (g.measureText(tt).width > maxW && line) { lines.push(line); line = w; } else line = tt; });
      if (line) lines.push(line);
      lines.forEach(function (L, i) { g.fillText(L, 18, this.h - bottom - (lines.length - 1 - i) * 15); }, this);
    }
    if (o.equations) {
      g.fillStyle = th.ink; g.font = "italic 13px 'TeX Gyre Pagella', Palatino, serif"; g.textAlign = "right"; g.textBaseline = "top";
      const eqs = String(this.scene.system).split(/\r?\n/).map(function (l) { return l.replace(/#.*$/, "").trim(); }).filter(function (l) { return /'|\[n\+1\]|dt\s*=|^noise/.test(l); });
      eqs.forEach(function (L, i) { g.fillText(L, this.w - 18, (this.view.axes ? 70 : 16) + 18 * i); }, this);
    }
    if (o.readout) { g.fillStyle = th.muted; g.font = "11px ui-monospace, monospace"; g.textAlign = "right"; g.textBaseline = "bottom"; g.fillText(this.readoutText(), this.w - 18, this.h - bottom); }
  };

  /* SVG of the current frame: vector paths for the views that keep them
     (trajectory, time series, phase plane, sweep), the rendered raster
     embedded for the others; text as SVG text either way. */
  Player.prototype.toSVG = function () {
    const th = this.theme, o = this.scene.overlay, w = this.w, h = this.h;
    const esc = function (s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); };
    let s = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="' + w + '" height="' + h + '" viewBox="0 0 ' + w + " " + h + '">';
    const bg = th.bg;
    if (bg[0] === "solid") s += '<rect width="100%" height="100%" fill="' + bg[1] + '"/>';
    else if (bg[0] === "radial") s += '<defs><radialGradient id="bg" cx="50%" cy="50%" r="72%"><stop offset="0" stop-color="' + bg[1] + '"/><stop offset="0.55" stop-color="' + bg[2] + '"/><stop offset="1" stop-color="' + bg[3] + '"/></radialGradient></defs><rect width="100%" height="100%" fill="url(#bg)"/>';
    const base = this.cv.base.toDataURL("image/png");
    s += '<image width="' + w + '" height="' + h + '" xlink:href="' + base + '"/>';
    const vec = this.view.svg ? this.view.svg(this) : null;
    if (vec !== null) s += '<g id="vector">' + vec + "</g>";
    else {
      const tmp = this.host.ownerDocument.createElement("canvas"); tmp.width = this.cv.top.width; tmp.height = this.cv.top.height;
      const g = tmp.getContext("2d"); g.drawImage(this.cv.trail, 0, 0); g.drawImage(this.cv.top, 0, 0);
      s += '<image width="' + w + '" height="' + h + '" xlink:href="' + tmp.toDataURL("image/png") + '"/>';
    }
    if (o.position !== "none") {
      if (o.title) s += '<text x="18" y="36" font-family="Jost, sans-serif" font-weight="600" font-size="22" fill="' + th.ink + '">' + esc(o.title) + "</text>";
      if (o.subtitle) s += '<text x="18" y="58" font-family="Jost, sans-serif" font-weight="300" font-size="14" fill="' + th.muted + '">' + esc(o.subtitle) + "</text>";
      if (o.caption) s += '<text x="18" y="' + (h - 12) + '" font-family="Jost, sans-serif" font-size="11" fill="' + th.muted + '">' + esc(o.caption) + "</text>";
    }
    return s + "</svg>";
  };

  Player.prototype.dispose = function () {
    this.pause();
    if (this.ro) this.ro.disconnect();
    if (this.io) this.io.disconnect();
    for (const k in this.cv) this.cv[k].remove();
    this.bgEl.remove(); this.overlayEl.remove();
  };

  DF.VIEW_DEFAULTS = VIEW_DEFAULTS;
  DF.normalizeScene = normalizeScene;
  DF.Player = Player;
})(globalThis.DynFlow = globalThis.DynFlow || {});

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
     browser's MediaRecorder; resolves to a WebM blob. */
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
        player.tick(); player.composite(c); k++;
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
    const per = Math.max(1, Math.round(60 / fps));
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
        for (let i = 0; i < per; i++) player.tick();
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
    return '<script src="' + (src || "dynflow.js") + '"></script>\n' +
      "<dyn-flow style=\"display:block;width:100%;height:420px\" controls scene='" + escAttr(JSON.stringify(scene)) + "'></dyn-flow>";
  };
  // A page that needs nothing else: the engine source, the scene and a full-window player.
  E.standaloneHTML = function (scene, opts) {
    opts = opts || {};
    if (!DF.SOURCE) throw new Error("The standalone export needs the built dynflow.js (run node tools/build.mjs)");
    const th = DF.THEMES[scene.style && scene.style.theme] || DF.THEMES["relab-night"];
    const bg = th.bg[0] === "solid" ? th.bg[1] : th.bg[0] === "radial" ? th.bg[3] : "#000";
    const title = (scene.overlay && scene.overlay.title) || scene.name || "DynFlow scene";
    return "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<title>" +
      title.replace(/</g, "&lt;") + "</title>\n<link rel=\"icon\" href=\"data:,\">\n<style>html,body{margin:0;height:100%;background:" + bg + "}dyn-flow{display:block;width:100vw;height:100vh}</style>\n</head>\n<body>\n" +
      "<dyn-flow" + (opts.controls === false ? "" : " controls") + " scene='" + escAttr(JSON.stringify(scene)) + "'></dyn-flow>\n" +
      "<script>\n" + DF.SOURCE.replace(/<\/script/gi, "<\\/script") + "\n</script>\n</body>\n</html>\n";
  };

  DF.Export = E;
})(globalThis.DynFlow = globalThis.DynFlow || {});

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
    { view: { type: "phase", seeds: 8 }, dt: 0.02, stepsPerFrame: 3, overlay: { equations: true } });

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
    { view: { type: "phase", seeds: 6 }, dt: 0.02, stepsPerFrame: 4, overlay: { equations: true } });

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
    { view: { type: "trajectory", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 12, overlay: { equations: true } });

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
    { view: { type: "flow", projection: "simplex", spawn: "mixed", life: [120, 420] }, n: 1800, dt: 0.05, stepsPerFrame: 3, style: { colorBy: "dominant" } });

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
    { view: { type: "phase", seeds: 10 }, dt: 0.02, stepsPerFrame: 4 });

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
    { view: { type: "timeseries", window: 60, members: 6 }, n: 6, spread: 0.5, dt: 0.05, stepsPerFrame: 2 });

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
    { view: { type: "trajectory", axes: ["x1", "x2", "x3"], warmup: 200, tail: 5000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 10 });

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
    { view: { type: "trajectory", warmup: 300, tail: 5000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 8 });

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
    { view: { type: "timeseries", vars: ["N1", "N2", "N3", "N4", "N5"], window: 800 }, dt: 0.1, stepsPerFrame: 20, keepPositive: true, style: { palette: "relab-qualitative" } });

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
    { view: { type: "sweep", param: "c", var: "V", from: 1, to: 3 }, dt: 0.05, stepsPerFrame: 8 });

  add("allee", "Strong Allee effect", "Ecology",
    "janos R/analysis_fokker_planck.R:854",
    "Populations below the threshold A decline to extinction; above it they grow to K.",
    `x' = r*x*(x/A - 1)*(1 - x/K)
     param r = 1 [0.2, 2]
     param A = 0.3 [0.05, 0.6]
     param K = 1 [0.7, 1.5]
     init x = 0.5
     range x = [0, 1.2]`,
    { view: { type: "timeseries", window: 20, members: 12 }, n: 12, initMode: "box", dt: 0.02, stepsPerFrame: 2 });

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
    { view: { type: "trajectory", tail: 60, dim3: false }, stepsPerFrame: 1, style: { pointSize: 4 } });

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
    { view: { type: "strobe", period: 10, transient: 10 }, n: 600, spread: 0.3, dt: 0.02, stepsPerFrame: 50, style: { pointSize: 1.3, alpha: 0.7 } });

  add("coleman", "Coleman logistic with seasonal carrying capacity", "Ecology",
    "nonautonomeR R/systems.R:1119 and vignettes/coleman-model.Rmd:58",
    "Logistic growth with periodic K(t): all positive orbits converge to one periodic orbit, the pullback attractor.",
    `x' = r*x*(1 - x/(1 + A*sin(w*t)))
     param r = 1 [0.2, 3]
     param A = 0.3 [0, 0.8]
     param w = 0.2 [0.05, 1]
     init x = 0.5
     range x = [0, 2]`,
    { view: { type: "timeseries", window: 60, members: 10 }, n: 10, initMode: "box", dt: 0.02, stepsPerFrame: 3 });

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
    { view: { type: "phase", seeds: 10 }, dt: 0.02, stepsPerFrame: 3 });

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
    { view: { type: "flow", projection: "simplex", life: [200, 500] }, n: 1200, dt: 0.02, stepsPerFrame: 4 });

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
    { view: { type: "timeseries", window: 600 }, dt: 0.1, stepsPerFrame: 10 });

  add("hutchinson", "Hutchinson delayed logistic", "Delay equations",
    "janos R/shiny_app.R:743 (Hutchinson 1948)",
    "Logistic growth with delayed feedback: the equilibrium K loses stability when r tau exceeds pi/2.",
    `N' = r*N*(1 - lag(N, tau)/K)
     param r = 1.6 [0.5, 2.5]
     param K = 100 [50, 150]
     param tau = 1 [0.2, 2]
     init N = 20
     range N = [0, 320]`,
    { view: { type: "timeseries", window: 40 }, dt: 0.02, stepsPerFrame: 4 });

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
    { view: { type: "timeseries", window: 300 }, dt: 0.1, stepsPerFrame: 6 });

  add("delayed-predator-prey", "Delayed predator and prey", "Delay equations",
    "janos vignettes/qualitative-analysis.Rmd:768",
    "Prey self-regulation acts with a delay tau, which destabilises coexistence into cycles.",
    `N' = r*N*(1 - lag(N, tau)/K) - a*N*P
     P' = b*N*P - d*P
     param r = 1.5 [0.5, 2.5]
     param K = 10 [5, 15]
     param a = 0.2 [0.1, 0.4]
     param b = 0.1 [0.05, 0.2]
     param d = 0.5 [0.2, 1]
     param tau = 3 [0.5, 5]
     init N = 5
     init P = 2
     range N = [0, 14]
     range P = [0, 12]`,
    { view: { type: "trajectory", tail: 3000 }, dt: 0.02, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", warmup: 500, tail: 3500, rotate: 0.25 }, dt: 0.005, stepsPerFrame: 5, overlay: { equations: true } });

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
    { view: { type: "trajectory", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", warmup: 400, tail: 5000, rotate: 0.2, axes: ["x", "z", "y"] }, dt: 0.005, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", warmup: 300, tail: 4000, rotate: 0.25 }, dt: 0.001, stepsPerFrame: 12 });

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
    { view: { type: "trajectory", warmup: 300, tail: 4000, rotate: 0.25 }, dt: 0.002, stepsPerFrame: 8 });

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
    { view: { type: "trajectory", warmup: 400, tail: 5000, rotate: 0.2 }, dt: 0.02, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", tail: 6000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", warmup: 400, tail: 5000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });

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
    { view: { type: "flow", life: [200, 600], rotate: 0.15 }, n: 1600, dt: 0.05, stepsPerFrame: 2, style: { colorBy: "speed", ramp: "relab-fire" } });

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
    { view: { type: "trajectory", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.005, stepsPerFrame: 8 });

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
    { view: { type: "flow", life: [150, 500], rotate: 0.2 }, n: 1500, dt: 0.01, stepsPerFrame: 3, style: { colorBy: "speed", ramp: "mako" } });

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
    { view: { type: "trajectory", warmup: 400, tail: 6000, rotate: 0.2 }, dt: 0.002, stepsPerFrame: 10 });

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
    { view: { type: "trajectory", warmup: 300, tail: 5000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", axes: ["x1", "x2", "x3"], warmup: 400, tail: 4000, rotate: 0.2 }, dt: 0.005, stepsPerFrame: 6 });

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
    { view: { type: "trajectory", axes: ["x", "y", "w"], warmup: 200, tail: 4000, rotate: 0.15 }, dt: 0.005, stepsPerFrame: 10 });

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
    { view: { type: "trajectory", warmup: 200, tail: 6000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 8 });

  add("charney-devore", "Charney-DeVore three modes", "Chaotic flows",
    "tuRbulence R/charney_devore.R:84 (Charney and DeVore 1979)",
    "Truncated barotropic flow over topography: zonal and blocked regimes of the midlatitude atmosphere.",
    `x' = k*(F - x) - alpha*y*z + beta*y
     y' = -k*y + alpha*x*z - beta*x - delta*z
     z' = -k*z + delta*y
     param F = 1.5 [0.5, 3]
     param k = 0.1 [0.05, 0.3]
     param alpha = 1 [0.5, 1.5]
     param beta = 0.5 [0.2, 1]
     param delta = 1 [0.5, 1.5]
     init x = 1
     init y = 0.1
     init z = 0.1`,
    { view: { type: "trajectory", warmup: 200, tail: 4000, rotate: 0.2 }, dt: 0.02, stepsPerFrame: 6 });

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
    { view: { type: "strobe", period: PI2 / 1.2, transient: 5 }, n: 500, spread: 1, dt: 0.02, stepsPerFrame: 60, style: { pointSize: 1.2, alpha: 0.6 } });

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
    { view: { type: "trajectory", tail: 2500 }, dt: 0.005, stepsPerFrame: 8 });

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
    { view: { type: "strobe", period: PI2, transient: 10 }, n: 600, spread: 3, dt: 0.02, stepsPerFrame: 60, style: { pointSize: 1.3 } });

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
    { view: { type: "strobe", period: PI2 / 1.3, transient: 0 }, n: 300, spread: 2, dt: 0.01, stepsPerFrame: 30, style: { fade: 0.02, pointSize: 2 } });

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
    { view: { type: "phase", seeds: 8 }, dt: 0.01, stepsPerFrame: 5 });

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
    { view: { type: "phase", seeds: 8 }, dt: 0.02, stepsPerFrame: 5 });

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
    { view: { type: "phase", seeds: 6 }, dt: 0.01, stepsPerFrame: 5 });

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
    { view: { type: "phase", seeds: 6 }, dt: 0.02, stepsPerFrame: 5 });

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
    { view: { type: "phase", seeds: 14 }, dt: 0.01, stepsPerFrame: 5 });

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
    { view: { type: "timeseries", window: 30 }, dt: 0.01, stepsPerFrame: 3 });

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
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "relab-fire" } });

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
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age" } });

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
    { view: { type: "flow", life: [40, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "mako" } });

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
    { view: { type: "flow", life: [400, 2000], axes: ["th", "p"], dim3: false }, n: 1500, stepsPerFrame: 1, style: { fade: 0, pointSize: 1, alpha: 0.5, colorBy: "member", palette: "relab" } });

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
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, spread: 0, initMode: "box", stepsPerFrame: 1, style: { fade: 0.6, pointSize: 1.4, alpha: 0.8, colorBy: "solid", palette: "mono-amber" } });

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
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, initMode: "box", stepsPerFrame: 1, style: { fade: 0.7, pointSize: 1.3, alpha: 0.8, colorBy: "solid", palette: "mono-ice" } });

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
    { view: { type: "flow", life: [300, 1200], dim3: false }, n: 2000, initMode: "box", stepsPerFrame: 1, style: { fade: 0.02, pointSize: 1.3, alpha: 0.6, colorBy: "var", colorVar: "x", ramp: "blackboard" } });

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
    { view: { type: "flow", life: [20, 60], dim3: false }, n: 3000, initMode: "box", stepsPerFrame: 1, style: { fade: 0.05, pointSize: 1.2, alpha: 0.5 } });

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
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", stepsPerFrame: 1, style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "age", ramp: "relab-fire" } });

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
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", stepsPerFrame: 1, style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "var", colorVar: "y", ramp: "mako" } });

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
    { view: { type: "density", window: 20 }, n: 3000, dt: 0.01, stepsPerFrame: 2, style: { ramp: "relab-fire" } });

  add("double-well", "Noisy double well", "Stochastic",
    "janos R/shiny_app.R:793 and nonautonomeR vignettes/melancholia-states.Rmd:298",
    "Bistable potential V = x^4/4 - x^2/2 with additive noise: Kramers escapes between the wells at x = -1 and x = 1.",
    `x' = x - x^3
     noise x = sigma
     param sigma = 0.45 [0, 1]
     init x = -1
     range x = [-2, 2]`,
    { view: { type: "density", window: 200 }, n: 3000, dt: 0.02, stepsPerFrame: 3, style: { ramp: "magma" } });

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
    { view: { type: "timeseries", window: 400, members: 1 }, dt: 0.02, stepsPerFrame: 15 });

  add("verhulst-noise", "Noise-induced transition (stochastic Verhulst)", "Stochastic",
    "janos R/shiny_app.R:784 (Horsthemke and Lefever 1984)",
    "Logistic growth with multiplicative noise: the stationary density changes shape (a P-bifurcation) as sigma crosses sqrt(a).",
    `x' = a*x - x^2
     noise x = sigma*x
     param a = 1 [0.2, 2]
     param sigma = 0.8 [0, 2]
     init x = 1
     range x = [0, 3]`,
    { view: { type: "density", window: 60 }, n: 3000, dt: 0.005, stepsPerFrame: 6, keepPositive: true, style: { ramp: "mako" } });

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
    { view: { type: "density", decay: 0.97 }, n: 3000, dt: 0.01, stepsPerFrame: 4, style: { ramp: "inferno" } });

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
    { view: { type: "timeseries", vars: ["v"], window: 40 }, dt: 0.002, stepsPerFrame: 20 });

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
    { view: { type: "density", decay: 0.96 }, n: 3000, dt: 0.01, stepsPerFrame: 4, style: { ramp: "magma" } });

  add("geometric-brownian", "Geometric Brownian motion", "Stochastic",
    "janos R/shiny_app.R:775 (gbm)",
    "Multiplicative noise with drift mu: the mean grows as exp(mu t) while the median grows as exp((mu - sigma^2/2) t).",
    `S' = mu*S
     noise S = sigma*S
     param mu = 0.08 [-0.2, 0.3]
     param sigma = 0.3 [0, 0.8]
     init S = 100
     range S = [0, 400]`,
    { view: { type: "timeseries", window: 10, members: 20 }, n: 20, dt: 0.002, stepsPerFrame: 10 });

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
    { view: { type: "density", decay: 0.9 }, n: 3000, spread: 0.01, dt: 0.01, stepsPerFrame: 3, style: { ramp: "relab-fire" } });

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
    { view: { type: "timeseries", window: 25, members: 1 }, dt: 0.005, stepsPerFrame: 4, overlay: { equations: true } });

  add("fold-normal-form", "Saddle-node (fold) normal form", "Tipping and nonautonomous",
    "janos vignettes/advanced-dynamics.Rmd:213",
    "x' = mu + x^2: two equilibria for mu < 0 that collide and vanish at mu = 0.",
    `x' = mu + x^2
     param mu = -1 [-1, 0.5]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "mu", from: -1, to: 0.3, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 8 });

  add("cusp", "Cusp catastrophe", "Tipping and nonautonomous",
    "normal form of the cusp catastrophe (Thom 1972; Zeeman 1977)",
    "x' = r + a x - x^3: for a > 0 two stable branches coexist between the folds at r = +/-2 (a/3)^(3/2), so a slow sweep of r jumps and shows hysteresis.",
    `x' = r + a*x - x^3
     param r = 0 [-0.8, 0.8]
     param a = 1 [-0.5, 2]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "r", from: -0.8, to: 0.8, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 10, overlay: { equations: true } });

  add("pitchfork", "Pitchfork normal form", "Tipping and nonautonomous",
    "janos R/analysis_bifurcation_sweep.R:80",
    "x' = a x - x^3: the origin splits into two symmetric stable states at a = 0.",
    `x' = a*x - x^3
     param a = 1 [-1, 1.5]
     init x = 0.05
     range x = [-1.4, 1.4]`,
    { view: { type: "sweep", param: "a", from: -1, to: 1.5, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 8 });

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
    { view: { type: "phase", seeds: 8 }, dt: 0.02, stepsPerFrame: 3 });

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
    { view: { type: "phase", seeds: 10 }, dt: 0.01, stepsPerFrame: 4 });

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
    { view: { type: "sweep", param: "eta2", var: "S", from: 0.5, to: 1.5, speed: 0.0004 }, dt: 0.02, stepsPerFrame: 10 });

  add("lorenz84-forced", "Lorenz-84 under climate change", "Tipping and nonautonomous",
    "nonautonomeR R/systems.R:651 and vignettes/ergodicity.Rmd:174 (Jánosi, Tél and co-authors)",
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
    { view: { type: "flow", life: "inf", rotate: 0.1 }, n: 1500, spread: 1.2, initMode: "ball", dt: 0.02, stepsPerFrame: 3, perturbations: [{ kind: "ramp", param: "F0", rate: -0.000274, t0: 0, span: -2 }], style: { fade: 0.25, colorBy: "speed", ramp: "relab-fire" }, overlay: { readout: true } });

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
    { view: { type: "flow", life: "inf", rotate: 0.12 }, n: 1500, initMode: "ball", spread: 6, dt: 0.005, stepsPerFrame: 4, perturbations: [{ kind: "ramp", param: "rho", rate: 0.02, t0: 0, span: 12 }], style: { fade: 0.2 }, overlay: { readout: true } });

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
    { view: { type: "strobe", period: PI2, transient: 3 }, n: 1500, spread: 1.5, initMode: "ball", dt: 2 * Math.PI / 128, stepsPerFrame: 128, perturbations: [{ kind: "ramp", param: "eps", rate: 0.00045, t0: 0, span: 0.3 }], style: { fade: 0.35, pointSize: 1.6, alpha: 0.8 }, overlay: { readout: true } });

  add("tilted-well", "Double well with a slow tilt", "Tipping and nonautonomous",
    "nonautonomeR vignettes/melancholia-states.Rmd:524",
    "The tilt lambda is ramped slowly: the occupied well loses stability at a fold, and the ensemble switches, earlier when noise is present.",
    `x' = x - x^3 + lambda
     noise x = sigma
     param lambda = -0.6 [-0.6, 0.6]
     param sigma = 0.15 [0, 0.5]
     init x = -1.2
     range x = [-1.8, 1.8]`,
    { view: { type: "density", window: 120 }, n: 3000, dt: 0.02, stepsPerFrame: 3, perturbations: [{ kind: "ramp", param: "lambda", rate: 0.01, t0: 0, span: 1.2 }], style: { ramp: "magma" }, overlay: { readout: true } });

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
    { view: { type: "timeseries", window: 160 }, dt: 0.05, stepsPerFrame: 4 });

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
    { view: { type: "timeseries", window: 20, vars: ["I"] }, dt: 0.0005, stepsPerFrame: 40 });

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
})(globalThis.DynFlow = globalThis.DynFlow || {});

// ---- src/component/dyn-flow.js
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* <dyn-flow>: a scene as an HTML element. One script tag and one element
   place a live figure in a web page, a Quarto or R Markdown document, a
   reveal.js slide or a pkgdown article.

     <script src="dynflow.js"></script>
     <dyn-flow scene='{"system": "...", "view": {"type": "trajectory"}}'></dyn-flow>
     <dyn-flow model="lorenz" theme="blackboard" controls></dyn-flow>
     <dyn-flow src="figures/fold.json" paused></dyn-flow>

   Attributes: scene (JSON), src (URL of a scene file), model (catalogue
   identifier), theme, view, title, controls (show play, restart and
   fullscreen buttons), paused (do not start automatically), static (draw a
   still after a number of frames given by the attribute, default 240). The
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

  class DynFlowElement extends HTMLElement {
    static get observedAttributes() { return ["scene", "model", "theme", "view"]; }
    connectedCallback() {
      if (this._root) return;
      this._root = this.attachShadow({ mode: "open" });
      const style = document.createElement("style"); style.textContent = CSS;
      this._stage = document.createElement("div"); this._stage.className = "stage";
      this._root.append(style, this._stage);
      if (this.hasAttribute("controls")) this.buildControls();
      this.loadScene();
    }
    disconnectedCallback() { if (this.player) { this.player.dispose(); this.player = null; } this._root = null; }
    attributeChangedCallback(name, oldV, newV) { if (this._root && oldV !== newV && this.player) this.loadScene(); }

    buildControls() {
      const c = document.createElement("div"); c.className = "ctl";
      const mk = (label, title, fn) => { const b = document.createElement("button"); b.textContent = label; b.title = title; b.setAttribute("aria-label", title); b.addEventListener("click", fn); c.appendChild(b); return b; };
      this._playBtn = mk("❚❚", "Pause", () => this.player && this.player.toggle());
      mk("↻", "Restart", () => this.restart());
      mk("⛶", "Full screen", () => { if (document.fullscreenElement) document.exitFullscreen(); else this.requestFullscreen && this.requestFullscreen(); });
      this._root.appendChild(c);
    }

    async loadScene() {
      let scene;
      try {
        if (this.hasAttribute("scene")) scene = JSON.parse(this.getAttribute("scene"));
        else if (this.hasAttribute("src")) scene = await (await fetch(this.getAttribute("src"))).json();
        else if (this.hasAttribute("model")) {
          if (!DF.sceneFor) throw new Error("The model catalogue is not loaded");
          scene = DF.sceneFor(this.getAttribute("model"));
        } else throw new Error("Give a scene, src or model attribute");
      } catch (e) { this.showError(e); return; }
      if (this.hasAttribute("theme")) scene.style = Object.assign({}, scene.style, { theme: this.getAttribute("theme") });
      if (this.hasAttribute("view")) scene.view = Object.assign({}, scene.view, { type: this.getAttribute("view") });
      if (this.hasAttribute("title")) scene.overlay = Object.assign({}, scene.overlay, { title: this.getAttribute("title") });
      if (this.player) { this.player.dispose(); this.player = null; }
      try {
        this.player = new DF.Player(this._stage, scene);
        if (this.player.error) throw this.player.error;
      } catch (e) { this.showError(e); return; }
      this.player.on("state", (s) => { if (this._playBtn) { this._playBtn.textContent = s === "play" ? "❚❚" : "▶"; this._playBtn.title = s === "play" ? "Pause" : "Play"; } });
      this.player.on("error", (e) => this.showError(e));
      const reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      if (this.hasAttribute("static") || reduced) this.player.advance(parseInt(this.getAttribute("static"), 10) || 240);
      else if (!this.hasAttribute("paused")) this.player.play();
      this.dispatchEvent(new CustomEvent("ready", { detail: this.player }));
    }
    showError(e) {
      const d = document.createElement("div"); d.className = "err";
      d.textContent = "DynFlow: " + (e && e.message ? e.message : String(e));
      this._stage.innerHTML = ""; this._stage.appendChild(d);
    }
    get scene() { return this.player ? this.player.getScene() : null; }
    play() { if (this.player) this.player.play(); }
    pause() { if (this.player) this.player.pause(); }
    restart(seed) { if (this.player) { this.player.restart(seed); if (!this.hasAttribute("paused")) this.player.play(); } }
    setParam(name, value) { if (this.player) this.player.setParam(name, value); }
  }

  if (!customElements.get("dyn-flow")) customElements.define("dyn-flow", DynFlowElement);
  DF.DynFlowElement = DynFlowElement;
})(globalThis.DynFlow = globalThis.DynFlow || {});

globalThis.DynFlow.VERSION = "0.1.0";
globalThis.DynFlow.SOURCE = "// DynFlow 0.1.0. Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab.\n// SPDX-License-Identifier: GPL-3.0-or-later. https://www.gnu.org/licenses/gpl-3.0.html\n// Built from 11 source files by tools/build.mjs; edit the sources, not this file.\n// ---- src/core/rng.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Seeded random numbers. Every stochastic scene is reproducible from its\n   seed: the same seed gives the same noise, the same jumps and the same\n   initial ensemble in the studio, in an exported page and in the tests. */\n(function (DF) {\n  \"use strict\";\n\n  // splitmix32 expands one 32-bit seed into the four words of xoshiro128**.\n  function splitmix32(a) {\n    return function () {\n      a |= 0; a = (a + 0x9e3779b9) | 0;\n      let z = a;\n      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b);\n      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35);\n      return (z ^ (z >>> 16)) >>> 0;\n    };\n  }\n\n  // xoshiro128** (Blackman and Vigna), period 2^128 - 1.\n  function RNG(seed) {\n    const sm = splitmix32(seed === undefined ? 1 : seed >>> 0);\n    this.s = new Uint32Array([sm(), sm(), sm(), sm()]);\n    this._spare = null;\n  }\n  RNG.prototype.nextU32 = function () {\n    const s = this.s;\n    const r = Math.imul(rotl(Math.imul(s[1], 5), 7), 9) >>> 0;\n    const t = s[1] << 9;\n    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];\n    s[2] ^= t; s[3] = rotl(s[3], 11);\n    return r;\n  };\n  function rotl(x, k) { return (x << k) | (x >>> (32 - k)); }\n\n  // Uniform on (0, 1), never exactly 0 or 1, with 53 random bits.\n  RNG.prototype.uniform = function () {\n    const hi = this.nextU32() >>> 5, lo = this.nextU32() >>> 6;\n    return (hi * 67108864 + lo + 0.5) / 9007199254740992;\n  };\n  RNG.prototype.range = function (a, b) { return a + (b - a) * this.uniform(); };\n\n  // Standard normal by the Box-Muller transform, one value cached.\n  RNG.prototype.normal = function () {\n    if (this._spare !== null) { const v = this._spare; this._spare = null; return v; }\n    const u = this.uniform(), v = this.uniform();\n    const r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * v;\n    this._spare = r * Math.sin(th);\n    return r * Math.cos(th);\n  };\n  RNG.prototype.exponential = function (rate) { return -Math.log(this.uniform()) / rate; };\n\n  // Poisson by inversion for small means and a normal approximation above 60.\n  RNG.prototype.poisson = function (mu) {\n    if (mu <= 0) return 0;\n    if (mu > 60) return Math.max(0, Math.round(mu + Math.sqrt(mu) * this.normal()));\n    const L = Math.exp(-mu);\n    let k = 0, p = 1;\n    do { k++; p *= this.uniform(); } while (p > L);\n    return k - 1;\n  };\n\n  // Symmetric alpha-stable variate, scale 1, by Chambers, Mallows and Stuck\n  // (1976); alpha = 2 gives a normal with variance 2, alpha = 1 a Cauchy.\n  RNG.prototype.stable = function (alpha) {\n    const V = Math.PI * (this.uniform() - 0.5), W = this.exponential(1);\n    if (Math.abs(alpha - 1) < 1e-12) return Math.tan(V);\n    return Math.sin(alpha * V) / Math.pow(Math.cos(V), 1 / alpha) *\n      Math.pow(Math.cos(V - alpha * V) / W, (1 - alpha) / alpha);\n  };\n\n  DF.RNG = RNG;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/core/expr.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Formula language. A system is written as plain text, one statement per\n   line, and compiled to JavaScript functions that write into preallocated\n   arrays. The parser accepts only numbers, declared names, whitelisted\n   functions and operators, so a compiled formula cannot reach anything else.\n\n     # Lotka-Volterra predator and prey\n     x' = a*x - b*x*y            differential equation (also dx/dt = ...)\n     y' = d*x*y - c*y\n     x[n+1] = r*x*(1 - x)        difference equation (map)\n     noise x = sigma*x           diffusion coefficient of dW_x (Ito)\n     aux h = x/(1 + x)           helper evaluated before the equations\n     param a = 1 [0, 3]          parameter, value and slider range\n     init x = 0.5                initial condition\n     range x = [0, 4]            axis range used by the views\n     lag(x, tau)                 delayed state x(t - tau), inside an expression\n\n   Functions: sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10\n   log2 sqrt cbrt abs sign floor ceil round min max pow mod step heaviside\n   ifelse clamp hill; random draws urand() and nrand() (fresh for every member)\n   and ucommon(k), ncommon(k) for k = 0..7 (shared by the ensemble at each\n   step, as the common forcing of a random dynamical system); constants pi and e; operators + - * / ^, unary minus,\n   comparisons < <= > >= == != (value 1 or 0), && and ||. */\n(function (DF) {\n  \"use strict\";\n\n  const FUNCS = {\n    sin: 1, cos: 1, tan: 1, asin: 1, acos: 1, atan: 1, sinh: 1, cosh: 1, tanh: 1,\n    exp: 1, log: 1, log10: 1, log2: 1, sqrt: 1, cbrt: 1, abs: 1, sign: 1,\n    floor: 1, ceil: 1, round: 1, atan2: 2, pow: 2, mod: 2, min: -1, max: -1,\n    step: 1, heaviside: 1, ifelse: 3, clamp: 3, hill: 3, lag: 2,\n    urand: 0, nrand: 0, ucommon: 1, ncommon: 1\n  };\n  const CONSTS = { pi: \"Math.PI\", e: \"Math.E\" };\n\n  function FormulaError(msg, line) { this.message = msg + (line ? \" (line \" + line + \")\" : \"\"); this.line = line; }\n  FormulaError.prototype = Object.create(Error.prototype);\n  FormulaError.prototype.name = \"FormulaError\";\n\n  // ---------------------------------------------------------------- lexer\n  function tokenize(src, line) {\n    const toks = [];\n    let i = 0;\n    while (i < src.length) {\n      const c = src[i];\n      if (c === \" \" || c === \"\\t\") { i++; continue; }\n      const num = /^(\\d+\\.?\\d*|\\.\\d+)([eE][+-]?\\d+)?/.exec(src.slice(i));\n      if (num) { toks.push({ t: \"num\", v: num[0] }); i += num[0].length; continue; }\n      const id = /^[A-Za-z_][A-Za-z0-9_]*/.exec(src.slice(i));\n      if (id) { toks.push({ t: \"id\", v: id[0] }); i += id[0].length; continue; }\n      const op = /^(<=|>=|==|!=|&&|\\|\\||[-+*/^(),<>!])/.exec(src.slice(i));\n      if (op) { toks.push({ t: \"op\", v: op[0] }); i += op[0].length; continue; }\n      throw new FormulaError(\"Unexpected character '\" + c + \"'\", line);\n    }\n    toks.push({ t: \"end\" });\n    return toks;\n  }\n\n  // ---------------------------------------------------------------- parser\n  // Pratt parser; binding powers follow the usual precedence, ^ is right\n  // associative and binds tighter than unary minus, so -x^2 = -(x^2).\n  const BP = { \"||\": 1, \"&&\": 2, \"==\": 3, \"!=\": 3, \"<\": 4, \"<=\": 4, \">\": 4, \">=\": 4, \"+\": 5, \"-\": 5, \"*\": 6, \"/\": 6, \"^\": 8 };\n\n  function parse(src, line) {\n    const toks = tokenize(src, line);\n    let k = 0;\n    const peek = function () { return toks[k]; };\n    const next = function () { return toks[k++]; };\n    function expect(v) {\n      const tk = next();\n      if (tk.t !== \"op\" || tk.v !== v) throw new FormulaError(\"Expected '\" + v + \"'\", line);\n    }\n    function nud(tk) {\n      if (tk.t === \"num\") return { k: \"num\", v: parseFloat(tk.v) };\n      if (tk.t === \"id\") {\n        if (peek().t === \"op\" && peek().v === \"(\") {\n          next();\n          const args = [];\n          if (!(peek().t === \"op\" && peek().v === \")\")) {\n            for (;;) { args.push(expr(0)); if (peek().t === \"op\" && peek().v === \",\") { next(); continue; } break; }\n          }\n          expect(\")\");\n          if (!(tk.v in FUNCS)) throw new FormulaError(\"Unknown function '\" + tk.v + \"'\", line);\n          const ar = FUNCS[tk.v];\n          if (ar >= 0 && args.length !== ar) throw new FormulaError(\"Function '\" + tk.v + \"' takes \" + ar + \" argument\" + (ar > 1 ? \"s\" : \"\"), line);\n          if (ar < 0 && args.length < 1) throw new FormulaError(\"Function '\" + tk.v + \"' needs arguments\", line);\n          return { k: \"call\", f: tk.v, args: args };\n        }\n        return { k: \"name\", v: tk.v };\n      }\n      if (tk.t === \"op\" && tk.v === \"(\") { const e = expr(0); expect(\")\"); return e; }\n      if (tk.t === \"op\" && tk.v === \"-\") return { k: \"neg\", a: expr(7) };\n      if (tk.t === \"op\" && tk.v === \"+\") return expr(7);\n      if (tk.t === \"op\" && tk.v === \"!\") return { k: \"not\", a: expr(7) };\n      throw new FormulaError(tk.t === \"end\" ? \"Incomplete expression\" : \"Unexpected '\" + tk.v + \"'\", line);\n    }\n    function expr(rbp) {\n      let left = nud(next());\n      for (;;) {\n        const tk = peek();\n        if (tk.t !== \"op\" || !(tk.v in BP) || BP[tk.v] <= rbp) break;\n        next();\n        const bp = BP[tk.v];\n        const right = expr(tk.v === \"^\" ? bp - 1 : bp);\n        left = { k: \"bin\", op: tk.v, a: left, b: right };\n      }\n      return left;\n    }\n    const ast = expr(0);\n    if (peek().t !== \"end\") throw new FormulaError(\"Unexpected '\" + peek().v + \"'\", line);\n    return ast;\n  }\n\n  // ------------------------------------------------------------ code emit\n  // scope maps a name to its JavaScript expression; lags collects lag() use.\n  function emit(node, scope, line, info) {\n    switch (node.k) {\n      case \"num\": return \"(\" + String(node.v) + \")\";\n      case \"name\":\n        if (node.v in scope) return scope[node.v];\n        if (node.v in CONSTS) return CONSTS[node.v];\n        if (node.v === \"t\") { info.usesTime = true; return \"t\"; }\n        throw new FormulaError(\"Unknown name '\" + node.v + \"'\", line);\n      case \"neg\": return \"(-\" + emit(node.a, scope, line, info) + \")\";\n      case \"not\": return \"(\" + emit(node.a, scope, line, info) + \" ? 0 : 1)\";\n      case \"bin\": {\n        const a = emit(node.a, scope, line, info), b = emit(node.b, scope, line, info);\n        if (node.op === \"^\") return \"Math.pow(\" + a + \", \" + b + \")\";\n        if (node.op === \"&&\" || node.op === \"||\") return \"((\" + a + \" \" + node.op + \" \" + b + \") ? 1 : 0)\";\n        if ([\"<\", \"<=\", \">\", \">=\", \"==\", \"!=\"].indexOf(node.op) >= 0) return \"((\" + a + \" \" + (node.op === \"==\" ? \"===\" : node.op === \"!=\" ? \"!==\" : node.op) + \" \" + b + \") ? 1 : 0)\";\n        return \"(\" + a + \" \" + node.op + \" \" + b + \")\";\n      }\n      case \"call\": {\n        const f = node.f;\n        if (f === \"lag\") {\n          const v = node.args[0];\n          if (v.k !== \"name\" || !(v.v in info.varIndex)) throw new FormulaError(\"The first argument of lag() must be a state variable\", line);\n          info.usesLag = true;\n          const tau = emit(node.args[1], scope, line, info);\n          info.lagExprs.push(tau);\n          return \"H(\" + info.varIndex[v.v] + \", t - (\" + tau + \"))\";\n        }\n        const args = node.args.map(function (a) { return emit(a, scope, line, info); });\n        switch (f) {\n          case \"mod\": return \"__mod(\" + args.join(\", \") + \")\";\n          case \"step\": case \"heaviside\": return \"((\" + args[0] + \") >= 0 ? 1 : 0)\";\n          case \"ifelse\": return \"((\" + args[0] + \") ? (\" + args[1] + \") : (\" + args[2] + \"))\";\n          case \"clamp\": return \"Math.min(Math.max(\" + args[0] + \", \" + args[1] + \"), \" + args[2] + \")\";\n          case \"hill\": return \"__hill(\" + args.join(\", \") + \")\";\n          case \"urand\": info.usesRandom = true; return \"R.u()\";\n          case \"nrand\": info.usesRandom = true; return \"R.n()\";\n          case \"ucommon\": info.usesRandom = true; return \"R.U[Math.min(7, Math.max(0, (\" + args[0] + \") | 0))]\";\n          case \"ncommon\": info.usesRandom = true; return \"R.N[Math.min(7, Math.max(0, (\" + args[0] + \") | 0))]\";\n          default: return \"Math.\" + f + \"(\" + args.join(\", \") + \")\";\n        }\n      }\n    }\n    throw new FormulaError(\"Cannot compile expression\", line);\n  }\n\n  const PRELUDE =\n    \"const __mod = function (a, b) { return a - b * Math.floor(a / b); };\\n\" +\n    \"const __hill = function (x, K, n) { const u = Math.pow(Math.max(x, 0), n); return u / (Math.pow(K, n) + u); };\\n\" +\n    // Without a random source (analysis, probes) random draws take their mean values.\n    \"const __R0 = { u: function () { return 0.5; }, n: function () { return 0; }, U: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], N: [0, 0, 0, 0, 0, 0, 0, 0] };\\n\";\n\n  // ------------------------------------------------------- system parsing\n  const RE = {\n    ode: /^(?:d\\s*([A-Za-z_]\\w*)\\s*\\/\\s*dt|([A-Za-z_]\\w*)\\s*')\\s*=\\s*(.+)$/,\n    map: /^([A-Za-z_]\\w*)\\s*(?:\\[\\s*n\\s*\\+\\s*1\\s*\\]|_\\{?\\s*n\\s*\\+\\s*1\\s*\\}?|_next)\\s*=\\s*(.+)$/,\n    noise: /^noise\\s+([A-Za-z_]\\w*)\\s*=\\s*(.+)$/,\n    aux: /^aux\\s+([A-Za-z_]\\w*)\\s*=\\s*(.+)$/,\n    param: /^param\\s+([A-Za-z_]\\w*)\\s*=\\s*([-+0-9.eE]+)\\s*(?:\\[\\s*([-+0-9.eE]+)\\s*,\\s*([-+0-9.eE]+)\\s*\\])?\\s*$/,\n    init: /^init\\s+([A-Za-z_]\\w*)\\s*=\\s*([-+0-9.eE]+)\\s*$/,\n    range: /^range\\s+([A-Za-z_]\\w*)\\s*=\\s*\\[\\s*([-+0-9.eE]+)\\s*,\\s*([-+0-9.eE]+)\\s*\\]\\s*$/\n  };\n  // Declared names may shadow the constants pi and e (e is a common\n  // conversion efficiency); function names and t are reserved.\n  const RESERVED = Object.assign({ t: 1, noise: 1, aux: 1, param: 1, init: 1, range: 1 }, FUNCS);\n\n  /* Parse and compile a system. Returns\n       { kind: 'ode' | 'map' | 'sde' | 'dde', vars, params, init, ranges,\n         f(t, x, p, dx, H), g(t, x, p, gx) or null, noiseMask, maxLag, source } */\n  function compileSystem(text) {\n    const lines = String(text).split(/\\r?\\n/);\n    const eqs = [], noises = [], auxes = [], params = [], init = {}, ranges = {};\n    let kind = null;\n    lines.forEach(function (raw, i) {\n      const ln = i + 1;\n      const s = raw.replace(/#.*$/, \"\").trim();\n      if (!s) return;\n      let m;\n      if ((m = RE.param.exec(s))) {\n        const v = +m[2];\n        const lo = m[3] !== undefined ? +m[3] : (v === 0 ? -1 : Math.min(0, 2 * v));\n        const hi = m[4] !== undefined ? +m[4] : (v === 0 ? 1 : Math.max(0, 2 * v));\n        if (params.some(function (q) { return q.name === m[1]; })) throw new FormulaError(\"Parameter '\" + m[1] + \"' is declared twice\", ln);\n        params.push({ name: m[1], value: v, min: lo, max: hi });\n      } else if ((m = RE.init.exec(s))) init[m[1]] = +m[2];\n      else if ((m = RE.range.exec(s))) ranges[m[1]] = [+m[2], +m[3]];\n      else if ((m = RE.noise.exec(s))) noises.push({ v: m[1], src: m[2], ln: ln });\n      else if ((m = RE.aux.exec(s))) auxes.push({ v: m[1], src: m[2], ln: ln });\n      else if ((m = RE.ode.exec(s))) {\n        if (kind === \"map\") throw new FormulaError(\"Differential and difference equations cannot be mixed\", ln);\n        kind = \"ode\"; eqs.push({ v: m[1] || m[2], src: m[3], ln: ln });\n      } else if ((m = RE.map.exec(s))) {\n        if (kind === \"ode\") throw new FormulaError(\"Differential and difference equations cannot be mixed\", ln);\n        kind = \"map\"; eqs.push({ v: m[1], src: m[2], ln: ln });\n      } else throw new FormulaError(\"Cannot read '\" + s + \"'\", ln);\n    });\n    if (!eqs.length) throw new FormulaError(\"No equations: write for example x' = -x\");\n\n    const vars = eqs.map(function (e) { return e.v; });\n    const varIndex = {};\n    vars.forEach(function (v, i) {\n      if (v in varIndex) throw new FormulaError(\"Variable '\" + v + \"' has two equations\", eqs[i].ln);\n      if (v in RESERVED) throw new FormulaError(\"'\" + v + \"' is a reserved name\", eqs[i].ln);\n      varIndex[v] = i;\n    });\n    params.forEach(function (q) {\n      if (q.name in varIndex) throw new FormulaError(\"'\" + q.name + \"' is both a variable and a parameter\");\n      if (q.name in RESERVED) throw new FormulaError(\"'\" + q.name + \"' is a reserved name\");\n    });\n\n    // Local names: v_<var>, p_<param>, a_<aux>; the scope maps formula names to them.\n    const scope = {};\n    vars.forEach(function (v) { scope[v] = \"v_\" + v; });\n    params.forEach(function (q) { scope[q.name] = \"p_\" + q.name; });\n    const info = { varIndex: varIndex, usesLag: false, usesTime: false, usesRandom: false, lagExprs: [] };\n    const head = vars.map(function (v, i) { return \"const v_\" + v + \" = x[\" + i + \"];\"; }).join(\" \") + \"\\n\" +\n      params.map(function (q, i) { return \"const p_\" + q.name + \" = p[\" + i + \"];\"; }).join(\" \") + \"\\n\";\n    let auxCode = \"\";\n    auxes.forEach(function (a) {\n      if (a.v in scope || a.v in RESERVED) throw new FormulaError(\"'\" + a.v + \"' is already defined\", a.ln);\n      auxCode += \"const a_\" + a.v + \" = \" + emit(parse(a.src, a.ln), scope, a.ln, info) + \";\\n\";\n      scope[a.v] = \"a_\" + a.v;\n    });\n    let body = \"\";\n    eqs.forEach(function (e, i) { body += \"dx[\" + i + \"] = \" + emit(parse(e.src, e.ln), scope, e.ln, info) + \";\\n\"; });\n    const usesLag = info.usesLag;\n\n    let gBody = \"\";\n    const noiseMask = vars.map(function () { return 0; });\n    noises.forEach(function (nz) {\n      if (!(nz.v in varIndex)) throw new FormulaError(\"noise refers to unknown variable '\" + nz.v + \"'\", nz.ln);\n      noiseMask[varIndex[nz.v]] = 1;\n      gBody += \"gx[\" + varIndex[nz.v] + \"] = \" + emit(parse(nz.src, nz.ln), scope, nz.ln, info) + \";\\n\";\n    });\n    if (noises.length && kind === \"map\") throw new FormulaError(\"noise lines apply to differential equations; add the noise inside the map instead\");\n\n    const f = new Function(PRELUDE + \"return function (t, x, p, dx, H, R) {\\nR = R || __R0;\\n\" + head + auxCode + body + \"};\")();\n    const g = noises.length ? new Function(PRELUDE + \"return function (t, x, p, gx, H, R) {\\nR = R || __R0;\\n\" + head + auxCode +\n      \"for (let i = 0; i < gx.length; i++) gx[i] = 0;\\n\" + gBody + \"};\")() : null;\n\n    // Constant delays are needed to size the history buffer; a delay that\n    // depends on parameters is evaluated at the current parameters.\n    const lagFns = info.lagExprs.map(function (src) { return new Function(\"p\", \"t\", head.split(\"\\n\")[1] + \"\\nreturn \" + src + \";\"); });\n\n    const k = usesLag ? \"dde\" : noises.length ? \"sde\" : kind;\n    if (usesLag && kind === \"map\") throw new FormulaError(\"lag() is available in differential equations only\");\n    return {\n      kind: k, time: kind === \"map\" ? \"discrete\" : \"continuous\",\n      vars: vars, params: params,\n      init: vars.map(function (v) { return v in init ? init[v] : 0.1; }),\n      ranges: ranges, f: f, g: g, noiseMask: noiseMask, usesTime: info.usesTime, usesRandom: info.usesRandom,\n      maxLag: function (p) { let m = 0; lagFns.forEach(function (fn) { m = Math.max(m, fn(p, 0)); }); return m; },\n      source: String(text)\n    };\n  }\n\n  // ---------------------------------------------------------------- LaTeX\n  const GREEK = [\"alpha\", \"beta\", \"gamma\", \"delta\", \"epsilon\", \"zeta\", \"eta\", \"theta\", \"iota\", \"kappa\", \"lambda\", \"mu\", \"nu\", \"xi\",\n    \"pi\", \"rho\", \"sigma\", \"tau\", \"upsilon\", \"phi\", \"chi\", \"psi\", \"omega\", \"Gamma\", \"Delta\", \"Theta\", \"Lambda\", \"Xi\", \"Pi\", \"Sigma\", \"Phi\", \"Psi\", \"Omega\"];\n  const GREEK_ALIAS = { eps: \"varepsilon\", lam: \"lambda\", gammag: \"gamma\", sig: \"sigma\" };\n\n  // x1 -> x_{1}, a12 -> a_{12}, alpha_2 -> \\alpha_{2}, K_m -> K_{m}, omega -> \\omega\n  function texName(name) {\n    let base = name, sub = \"\";\n    const us = name.indexOf(\"_\");\n    if (us > 0) { base = name.slice(0, us); sub = name.slice(us + 1); }\n    else {\n      const m = /^([A-Za-z]+?)(\\d+)$/.exec(name);\n      if (m) { base = m[1]; sub = m[2]; }\n    }\n    if (base in GREEK_ALIAS) base = GREEK_ALIAS[base];\n    let b = GREEK.indexOf(base) >= 0 ? \"\\\\\" + base : base.length > 1 ? \"\\\\mathrm{\" + base + \"}\" : base;\n    if (sub) b += \"_{\" + (GREEK.indexOf(sub) >= 0 ? \"\\\\\" + sub : sub.length > 1 && !/^\\d+$/.test(sub) ? \"\\\\mathrm{\" + sub + \"}\" : sub) + \"}\";\n    return b;\n  }\n\n  // \\dot over the base symbol only: x1 -> \\dot{x}_{1}\n  function dotted(name) {\n    const t = texName(name), m = /^(\\\\?[A-Za-z]+|\\\\mathrm\\{[^}]*\\})(_\\{.*\\})?$/.exec(t);\n    return m ? \"\\\\dot{\" + m[1] + \"}\" + (m[2] || \"\") : \"\\\\dot{\" + t + \"}\";\n  }\n\n  const TEX_PREC = { \"||\": 1, \"&&\": 2, \"==\": 3, \"!=\": 3, \"<\": 4, \"<=\": 4, \">\": 4, \">=\": 4, \"+\": 5, \"-\": 5, \"*\": 6, \"/\": 7, \"^\": 8 };\n  function tex(node, ctx) {\n    const wrap = function (child, minPrec) {\n      const s = tex(child, ctx);\n      const pr = child.k === \"bin\" ? TEX_PREC[child.op] : child.k === \"neg\" ? 5.4 : child.k === \"call\" && child.f === \"lag\" ? 8.5 : 9;\n      return pr < minPrec ? \"\\\\left(\" + s + \"\\\\right)\" : s;\n    };\n    switch (node.k) {\n      case \"num\": return String(node.v).replace(/e([+-]?\\d+)$/, \" \\\\times 10^{$1}\");\n      case \"name\":\n        if (node.v === \"t\") return \"t\";\n        if (node.v === \"pi\" && !(node.v in ctx.declared)) return \"\\\\pi\";\n        if (node.v in ctx.aux) return texName(node.v);\n        return texName(node.v);\n      case \"neg\": return \"-\" + wrap(node.a, 6);\n      case \"not\": return \"\\\\neg \" + wrap(node.a, 9);\n      case \"bin\": {\n        const op = node.op;\n        if (op === \"/\") return \"\\\\frac{\" + tex(node.a, ctx) + \"}{\" + tex(node.b, ctx) + \"}\";\n        if (op === \"^\") return wrap(node.a, 9) + \"^{\" + tex(node.b, ctx) + \"}\";\n        if (op === \"*\") {\n          const a = wrap(node.a, 6), b = wrap(node.b, 6);\n          const numRight = node.b.k === \"num\";\n          return a + (numRight || /^[\\d.]/.test(b) ? \" \\\\cdot \" : \"\\\\,\") + b;\n        }\n        const sym = { \"+\": \" + \", \"-\": \" - \", \"<\": \" < \", \"<=\": \" \\\\le \", \">\": \" > \", \">=\": \" \\\\ge \", \"==\": \" = \", \"!=\": \" \\\\ne \", \"&&\": \" \\\\land \", \"||\": \" \\\\lor \" }[op];\n        return wrap(node.a, TEX_PREC[op]) + sym + wrap(node.b, op === \"-\" || op === \"+\" ? TEX_PREC[op] + 0.5 : TEX_PREC[op]);\n      }\n      case \"call\": {\n        const a = node.args.map(function (x) { return tex(x, ctx); });\n        switch (node.f) {\n          case \"sqrt\": return \"\\\\sqrt{\" + a[0] + \"}\";\n          case \"cbrt\": return \"\\\\sqrt[3]{\" + a[0] + \"}\";\n          case \"abs\": return \"\\\\left|\" + a[0] + \"\\\\right|\";\n          case \"exp\": return \"e^{\" + a[0] + \"}\";\n          case \"pow\": return wrap(node.args[0], 9) + \"^{\" + a[1] + \"}\";\n          case \"lag\": return tex(node.args[0], ctx) + \"(t - \" + a[1] + \")\";\n          case \"hill\": return \"\\\\frac{\" + wrap(node.args[0], 9) + \"^{\" + a[2] + \"}}{\" + wrap(node.args[1], 9) + \"^{\" + a[2] + \"} + \" + wrap(node.args[0], 9) + \"^{\" + a[2] + \"}}\";\n          case \"step\": case \"heaviside\": return \"\\\\Theta\\\\left(\" + a[0] + \"\\\\right)\";\n          case \"floor\": return \"\\\\lfloor \" + a[0] + \" \\\\rfloor\";\n          case \"mod\": return a[0] + \" \\\\bmod \" + wrap(node.args[1], 9);\n          case \"ifelse\": return \"\\\\begin{cases} \" + a[1] + \" & \" + a[0] + \" \\\\\\\\ \" + a[2] + \" & \\\\text{otherwise} \\\\end{cases}\";\n          case \"clamp\": return \"\\\\mathrm{clamp}\\\\left(\" + a.join(\", \") + \"\\\\right)\";\n          case \"sin\": case \"cos\": case \"tan\": case \"sinh\": case \"cosh\": case \"tanh\": case \"log\": case \"min\": case \"max\":\n            return \"\\\\\" + node.f + \"\\\\left(\" + a.join(\", \") + \"\\\\right)\";\n          case \"asin\": case \"acos\": case \"atan\": return \"\\\\\" + node.f.replace(\"a\", \"arc\") + \"\\\\left(\" + a[0] + \"\\\\right)\";\n          default: return \"\\\\mathrm{\" + node.f + \"}\\\\left(\" + a.join(\", \") + \"\\\\right)\";\n        }\n      }\n    }\n    return \"\";\n  }\n\n  /* LaTeX for every statement of a system: equations, noise terms and aux\n     definitions, in source order. Parameters are listed separately by callers. */\n  function systemLatex(text) {\n    const out = [], declared = {}, aux = {};\n    let kind = \"ode\";\n    const lines = String(text).split(/\\r?\\n/).map(function (r) { return r.replace(/#.*$/, \"\").trim(); }).filter(Boolean);\n    lines.forEach(function (s) { let m; if ((m = RE.param.exec(s))) declared[m[1]] = 1; if (RE.map.test(s)) kind = \"map\"; });\n    const noiseOf = {};\n    lines.forEach(function (s) { const m = RE.noise.exec(s); if (m) noiseOf[m[1]] = m[2]; });\n    const ctx = { declared: declared, aux: aux };\n    lines.forEach(function (s) {\n      let m;\n      if ((m = RE.aux.exec(s))) { aux[m[1]] = 1; out.push(texName(m[1]) + \" = \" + tex(parse(m[2]), ctx)); }\n      else if ((m = RE.ode.exec(s))) {\n        const v = m[1] || m[2], rhs = tex(parse(m[3]), ctx);\n        if (noiseOf[v] !== undefined) out.push(\"\\\\mathrm{d}\" + texName(v) + \" = \\\\left(\" + rhs + \"\\\\right)\\\\mathrm{d}t + \" + tex(parse(noiseOf[v]), ctx) + \"\\\\,\\\\mathrm{d}W_{\" + texName(v) + \"}\");\n        else out.push(dotted(v) + \" = \" + rhs);\n      } else if ((m = RE.map.exec(s))) out.push(texName(m[1]).replace(/_\\{(.*)\\}$/, \"_{$1,\\\\,n+1}\").replace(/^([^_]*)$/, \"$1_{n+1}\") + \" = \" + tex(parse(m[2]), ctx));\n    });\n    return { lines: out, kind: kind };\n  }\n\n  DF.FormulaError = FormulaError;\n  DF.parseExpression = parse;\n  DF.compileSystem = compileSystem;\n  DF.systemLatex = systemLatex;\n  DF.texName = texName;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/core/sim.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Simulation engine. A Simulator advances an ensemble of n copies of one\n   compiled system, all sharing parameters, time and any common forcing, with\n   a fixed step h:\n\n     ode   classical fourth-order Runge-Kutta\n     sde   Euler-Maruyama (Ito): x += f h + g dW, dW ~ N(0, h)\n     dde   Runge-Kutta 4 with a cubic Hermite interpolant of the stored\n           history; constant initial history x(s) = x0 for s <= t0\n     map   x[n+1] = F(x[n]), one iteration per step\n\n   Perturbations modify either parameters (periodic, quasiperiodic, ramp,\n   step, Ornstein-Uhlenbeck) or states (additive, multiplicative, coloured and\n   alpha-stable noise, Poisson jumps, periodic pulses). A perturbation marked\n   common uses one realisation for the whole ensemble, as an environmental\n   forcing does; otherwise every member receives its own. */\n(function (DF) {\n  \"use strict\";\n\n  // ------------------------------------------------------------ RK4 step\n  function makeRK4(dim) {\n    const k1 = new Float64Array(dim), k2 = new Float64Array(dim), k3 = new Float64Array(dim), k4 = new Float64Array(dim), y = new Float64Array(dim);\n    return function (f, t, x, p, h, H, R) {\n      f(t, x, p, k1, H, R);\n      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k1[i];\n      f(t + 0.5 * h, y, p, k2, H, R);\n      for (let i = 0; i < dim; i++) y[i] = x[i] + 0.5 * h * k2[i];\n      f(t + 0.5 * h, y, p, k3, H, R);\n      for (let i = 0; i < dim; i++) y[i] = x[i] + h * k3[i];\n      f(t + h, y, p, k4, H, R);\n      for (let i = 0; i < dim; i++) x[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);\n    };\n  }\n\n  // ----------------------------------------------------- DDE history\n  // Ring buffer of (t, x, dx) at the step points; H(i, s) evaluates the\n  // cubic Hermite interpolant, exact for cubic solutions.\n  function History(dim, cap, x0, t0) {\n    this.dim = dim; this.cap = cap; this.x0 = Float64Array.from(x0); this.t0 = t0;\n    this.T = new Float64Array(cap); this.X = new Float64Array(cap * dim); this.D = new Float64Array(cap * dim);\n    this.len = 0; this.head = 0;\n  }\n  History.prototype.push = function (t, x, dx) {\n    const j = this.head;\n    this.T[j] = t;\n    for (let i = 0; i < this.dim; i++) { this.X[j * this.dim + i] = x[i]; this.D[j * this.dim + i] = dx[i]; }\n    this.head = (j + 1) % this.cap; this.len = Math.min(this.len + 1, this.cap);\n  };\n  History.prototype.at = function (i, s) {\n    if (s <= this.t0 || this.len === 0) return this.x0[i];\n    const cap = this.cap, dim = this.dim;\n    const newest = (this.head - 1 + cap) % cap, oldest = (this.head - this.len + cap) % cap;\n    if (s >= this.T[newest]) return this.X[newest * dim + i];\n    // binary search over the logical order 0..len-1\n    let lo = 0, hi = this.len - 1;\n    const idx = function (q) { return (oldest + q) % cap; };\n    if (s <= this.T[idx(0)]) return this.X[idx(0) * dim + i];\n    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (this.T[idx(mid)] <= s) lo = mid; else hi = mid; }\n    const a = idx(lo), b = idx(hi);\n    const ta = this.T[a], hstep = this.T[b] - ta, u = (s - ta) / hstep;\n    const ya = this.X[a * dim + i], yb = this.X[b * dim + i], da = this.D[a * dim + i], db = this.D[b * dim + i];\n    const u2 = u * u, u3 = u2 * u;\n    return (2 * u3 - 3 * u2 + 1) * ya + (u3 - 2 * u2 + u) * hstep * da + (-2 * u3 + 3 * u2) * yb + (u3 - u2) * hstep * db;\n  };\n\n  // --------------------------------------------------------- perturbations\n  /* Parameter modulators, applied to the base value b of parameter `param`:\n       periodic       b + A sin(2 pi t / T + phase)\n       quasiperiodic  b + A sin(2 pi t / T) + A2 sin(2 pi t / T2)\n       ramp           b + rate (t - t0), clipped to [b, b + span] or [b + span, b]\n       step           b before t0, b + A after\n       ou             b + eta, d eta = -eta / tau dt + sigma sqrt(2 / tau) dW, common\n     State perturbations, on variable `var`:\n       additive       dx += sigma dW\n       multiplicative dx += sigma x dW\n       coloured       dx += eta dt, eta an OU process with time scale tau\n       levy           dx += sigma h^(1/alpha) S_alpha, S_alpha symmetric stable\n       jumps          at rate lambda, x += size, or x *= (1 - frac)\n       pulse          every T time units, x += size, or x *= (1 - frac)  */\n  const PARAM_KINDS = { periodic: 1, quasiperiodic: 1, ramp: 1, step: 1, ou: 1 };\n\n  function Simulator(system, opts) {\n    opts = opts || {};\n    this.sys = system;\n    this.dim = system.vars.length;\n    this.n = Math.max(1, opts.n || 1);\n    this.h = opts.dt || 0.01;\n    this.seed = opts.seed === undefined ? 1 : opts.seed;\n    this.base = Float64Array.from(opts.params || system.params.map(function (q) { return q.value; }));\n    this.p = Float64Array.from(this.base);\n    this.perturbations = (opts.perturbations || []).map(function (q) { return Object.assign({}, q); });\n    this.initMode = opts.initMode || \"point\";\n    this.init = Float64Array.from(opts.init || system.init);\n    this.spread = opts.spread === undefined ? 0.05 : opts.spread;\n    this.box = opts.box || null;\n    this.keepPositive = !!opts.keepPositive;\n    this.maxHistory = opts.maxHistory || 20000;\n    this.reset();\n  }\n\n  Simulator.prototype.reset = function (seed) {\n    if (seed !== undefined) this.seed = seed;\n    const n = this.n, dim = this.dim;\n    this.rng = new DF.RNG(this.seed);\n    const rng = this.rng;\n    this.R = { u: function () { return rng.uniform(); }, n: function () { return rng.normal(); }, U: new Float64Array(8), N: new Float64Array(8) };\n    this.t = 0; this.steps = 0;\n    this.X = new Float64Array(n * dim);\n    this.alive = new Uint8Array(n).fill(1);\n    for (let k = 0; k < n; k++) this.initMember(k);\n    this.rk4 = makeRK4(dim);\n    this.dx = new Float64Array(dim); this.gx = new Float64Array(dim); this.xk = new Float64Array(dim);\n    this.eta = new Float64Array(this.perturbations.length);         // common OU states\n    this.etaK = new Float64Array(this.perturbations.length * n);    // per-member OU states\n    this.nextPulse = this.perturbations.map(function (q) { return q.kind === \"pulse\" ? (q.t0 || q.period || 1) : Infinity; });\n    this.hist = null;\n    this.updateParams();\n    if (this.sys.kind === \"dde\") {\n      const lag = this.sys.maxLag(this.p);\n      const cap = Math.min(this.maxHistory, Math.ceil(lag / this.h) + 8);\n      this.hist = [];\n      for (let k = 0; k < n; k++) this.hist.push(this.newHistory(k, cap));\n    }\n  };\n\n  // History of member k starting now, holding the initial point with the\n  // right derivative f(t0+, x0) so that the first step interpolates correctly.\n  Simulator.prototype.newHistory = function (k, cap) {\n    const dim = this.dim, x0 = this.X.slice(k * dim, k * dim + dim);\n    const hk = new History(dim, cap, x0, this.t);\n    const dx = new Float64Array(dim);\n    this.sys.f(this.t, x0, this.p, dx, function (i) { return x0[i]; });\n    hk.push(this.t, x0, dx);\n    return hk;\n  };\n\n  Simulator.prototype.initMember = function (k) {\n    const dim = this.dim, r = this.rng, x = this.X;\n    for (let i = 0; i < dim; i++) {\n      let v;\n      if (this.initMode === \"box\" && this.box) v = r.range(this.box[i][0], this.box[i][1]);\n      else if (this.initMode === \"ball\") v = this.init[i] + this.spread * r.normal();\n      else v = this.init[i] + (k === 0 ? 0 : this.spread * (r.uniform() - 0.5));\n      x[k * dim + i] = v;\n    }\n    this.alive[k] = 1;\n  };\n\n  // Place member k at state s (used by clicks and by respawning views).\n  Simulator.prototype.setMember = function (k, s) {\n    for (let i = 0; i < this.dim; i++) this.X[k * this.dim + i] = s[i];\n    this.alive[k] = 1;\n    if (this.hist) this.hist[k] = this.newHistory(k, this.hist[k].cap);\n  };\n\n  Simulator.prototype.paramIndex = function (name) {\n    return this.sys.params.findIndex(function (q) { return q.name === name; });\n  };\n\n  Simulator.prototype.updateParams = function () {\n    this.p.set(this.base);\n    const t = this.t;\n    for (let j = 0; j < this.perturbations.length; j++) {\n      const q = this.perturbations[j];\n      if (!(q.kind in PARAM_KINDS) || q.enabled === false) continue;\n      const i = this.paramIndex(q.param);\n      if (i < 0) continue;\n      const b = this.base[i];\n      switch (q.kind) {\n        case \"periodic\": this.p[i] = b + q.amp * Math.sin(2 * Math.PI * t / q.period + (q.phase || 0)); break;\n        case \"quasiperiodic\": this.p[i] = b + q.amp * Math.sin(2 * Math.PI * t / q.period) + (q.amp2 === undefined ? q.amp : q.amp2) * Math.sin(2 * Math.PI * t / (q.period2 || q.period * (1 + Math.sqrt(5)) / 2)); break;\n        case \"ramp\": {\n          const d = q.rate * Math.max(0, t - (q.t0 || 0));\n          const span = q.span === undefined ? Infinity : q.span;\n          this.p[i] = b + (span >= 0 ? Math.min(d, span) : Math.max(d, span));\n          break;\n        }\n        case \"step\": this.p[i] = t >= (q.t0 || 0) ? b + q.amp : b; break;\n        case \"ou\": this.p[i] = b + this.eta[j]; break;\n      }\n    }\n  };\n\n  // Advance the whole ensemble by one step of length h.\n  Simulator.prototype.step = function () {\n    const sys = this.sys, n = this.n, dim = this.dim, h = this.h, t = this.t, r = this.rng;\n    const X = this.X, p = this.p, pert = this.perturbations;\n    const sqh = Math.sqrt(h);\n    // Common noise increments, drawn once per step for the whole ensemble.\n    if (sys.usesRandom) for (let i = 0; i < 8; i++) { this.R.U[i] = r.uniform(); this.R.N[i] = r.normal(); }\n    const R = this.R;\n    const commonDW = pert.map(function (q) { return q.common && q.enabled !== false ? r.normal() : 0; });\n    const commonJump = pert.map(function (q) { return q.kind === \"jumps\" && q.common && q.enabled !== false ? r.poisson(q.rate * h) : 0; });\n    const commonStable = pert.map(function (q) { return q.kind === \"levy\" && q.common && q.enabled !== false ? r.stable(q.alpha || 1.5) : 0; });\n\n    for (let k = 0; k < n; k++) {\n      if (!this.alive[k]) continue;\n      const x = X.subarray(k * dim, k * dim + dim);\n      const Hk = this.hist ? this.hist[k] : null;\n      const H = Hk ? function (i, s) { return Hk.at(i, s); } : undefined;\n      if (sys.time === \"discrete\") {\n        sys.f(t, x, p, this.dx, H, R);\n        x.set(this.dx);\n      } else if (sys.kind === \"sde\") {\n        sys.f(t, x, p, this.dx, H, R);\n        sys.g(t, x, p, this.gx, H, R);\n        for (let i = 0; i < dim; i++) x[i] += this.dx[i] * h + (sys.noiseMask[i] ? this.gx[i] * sqh * r.normal() : 0);\n      } else {\n        this.rk4(sys.f, t, x, p, h, H, R);\n      }\n      // State perturbations; in discrete time h is 1 and dW has variance 1.\n      for (let j = 0; j < pert.length; j++) {\n        const q = pert[j];\n        if (q.enabled === false || q.kind in PARAM_KINDS || q.kind === \"pulse\") continue;\n        const i = sys.vars.indexOf(q.var);\n        if (i < 0) continue;\n        const dW = (q.common ? commonDW[j] : r.normal()) * sqh;\n        switch (q.kind) {\n          case \"additive\": x[i] += q.sigma * dW; break;\n          case \"multiplicative\": x[i] += q.sigma * x[i] * dW; break;\n          case \"coloured\": {\n            const tau = q.tau || 1, ix = q.common ? j : j * n + k, arr = q.common ? this.eta : this.etaK;\n            if (!q.common || k === 0) arr[ix] += -arr[ix] / tau * h + q.sigma * Math.sqrt(2 / tau) * (q.common ? commonDW[j] : r.normal()) * sqh;\n            x[i] += arr[ix] * h;\n            break;\n          }\n          case \"levy\": x[i] += q.sigma * Math.pow(h, 1 / (q.alpha || 1.5)) * (q.common ? commonStable[j] : r.stable(q.alpha || 1.5)); break;\n          case \"jumps\": {\n            const m = q.common ? commonJump[j] : r.poisson(q.rate * h);\n            for (let c = 0; c < m; c++) {\n              if (q.frac !== undefined && q.frac !== null && q.frac !== \"\") x[i] *= (1 - q.frac);\n              else x[i] += q.size + (q.sd ? q.sd * r.normal() : 0);\n            }\n            break;\n          }\n        }\n      }\n      if (this.keepPositive) for (let i = 0; i < dim; i++) if (x[i] < 0) x[i] = 0;\n      let ok = true;\n      for (let i = 0; i < dim; i++) if (!isFinite(x[i]) || Math.abs(x[i]) > 1e12) { ok = false; break; }\n      if (!ok) this.alive[k] = 0;\n      if (Hk && ok) { sys.f(t + h, x, p, this.dx, H); Hk.push(t + h, x, this.dx); }\n    }\n    this.t = t + (sys.time === \"discrete\" ? 1 : h);\n    this.steps++;\n    // Pulses at fixed times act on every member at once.\n    for (let j = 0; j < pert.length; j++) {\n      const q = pert[j];\n      if (q.kind !== \"pulse\" || q.enabled === false) continue;\n      const i = sys.vars.indexOf(q.var);\n      while (i >= 0 && this.t >= this.nextPulse[j]) {\n        for (let k = 0; k < n; k++) {\n          const ix = k * dim + i;\n          if (q.frac !== undefined && q.frac !== null && q.frac !== \"\") X[ix] *= (1 - q.frac); else X[ix] += q.size;\n        }\n        this.nextPulse[j] += q.period;\n      }\n    }\n    // Common Ornstein-Uhlenbeck parameter noise.\n    for (let j = 0; j < pert.length; j++) {\n      const q = pert[j];\n      if (q.kind !== \"ou\" || q.enabled === false) continue;\n      const tau = q.tau || 1;\n      this.eta[j] += -this.eta[j] / tau * h + q.sigma * Math.sqrt(2 / tau) * r.normal() * sqh;\n    }\n    this.updateParams();\n  };\n\n  Simulator.prototype.member = function (k, out) {\n    out = out || new Float64Array(this.dim);\n    for (let i = 0; i < this.dim; i++) out[i] = this.X[k * this.dim + i];\n    return out;\n  };\n\n  // --------------------------------------------------- Lyapunov estimate\n  /* Largest Lyapunov exponent by two nearby orbits renormalised every\n     `every` steps (Benettin et al. 1980), after a transient. Continuous time\n     gives a rate per unit time; maps give a rate per iteration. */\n  function largestLyapunov(system, params, init, opts) {\n    opts = opts || {};\n    const h = opts.dt || 0.01, steps = opts.steps || 100000, every = opts.every || 10, d0 = opts.d0 || 1e-8;\n    const transient = opts.transient || 5000;\n    const dim = system.vars.length, rk4 = makeRK4(dim), tmp = new Float64Array(dim);\n    const p = Float64Array.from(params);\n    const x = Float64Array.from(init), y = new Float64Array(dim);\n    const adv = function (z, t) {\n      if (system.time === \"discrete\") { system.f(t, z, p, tmp); z.set(tmp); } else rk4(system.f, t, z, p, h);\n    };\n    let t = 0;\n    for (let s = 0; s < transient; s++) { adv(x, t); t += system.time === \"discrete\" ? 1 : h; }\n    y.set(x); y[0] += d0;\n    let sum = 0;\n    for (let s = 1; s <= steps; s++) {\n      adv(x, t); adv(y, t); t += system.time === \"discrete\" ? 1 : h;\n      if (s % every === 0) {\n        let d = 0;\n        for (let i = 0; i < dim; i++) d += (y[i] - x[i]) * (y[i] - x[i]);\n        d = Math.sqrt(d);\n        sum += Math.log(d / d0);\n        for (let i = 0; i < dim; i++) y[i] = x[i] + d0 * (y[i] - x[i]) / d;\n      }\n    }\n    return sum / (steps * (system.time === \"discrete\" ? 1 : h));\n  }\n\n  DF.makeRK4 = makeRK4;\n  DF.History = History;\n  DF.Simulator = Simulator;\n  DF.largestLyapunov = largestLyapunov;\n  DF.PARAM_PERTURBATIONS = Object.keys(PARAM_KINDS);\n  DF.STATE_PERTURBATIONS = [\"additive\", \"multiplicative\", \"coloured\", \"levy\", \"jumps\", \"pulse\"];\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/core/analysis.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Local analysis: Jacobians by central differences, eigenvalues of a real\n   matrix, equilibria (flows) and fixed points (maps) by damped Newton\n   iteration from many seeds, and their stability. */\n(function (DF) {\n  \"use strict\";\n\n  // G(x) = f(x) for flows and F(x) - x for maps, so both reduce to G = 0.\n  function residual(sys, t, x, p, out) {\n    sys.f(t, x, p, out, function (i) { return x[i]; });\n    if (sys.time === \"discrete\") for (let i = 0; i < x.length; i++) out[i] -= x[i];\n    return out;\n  }\n\n  // Jacobian of f (flows) or of F (maps), row-major, central differences.\n  function jacobian(sys, t, x, p) {\n    const n = x.length, J = new Float64Array(n * n), a = new Float64Array(n), b = new Float64Array(n), y = Float64Array.from(x);\n    const H = function (i) { return y[i]; };\n    for (let j = 0; j < n; j++) {\n      const hj = 1e-6 * Math.max(1, Math.abs(x[j]));\n      y[j] = x[j] + hj; sys.f(t, y, p, a, H);\n      y[j] = x[j] - hj; sys.f(t, y, p, b, H);\n      y[j] = x[j];\n      for (let i = 0; i < n; i++) J[i * n + j] = (a[i] - b[i]) / (2 * hj);\n    }\n    return J;\n  }\n\n  // Solve A z = r in place (Gaussian elimination, partial pivoting); false if singular.\n  function solve(A, r, n) {\n    for (let c = 0; c < n; c++) {\n      let piv = c;\n      for (let i = c + 1; i < n; i++) if (Math.abs(A[i * n + c]) > Math.abs(A[piv * n + c])) piv = i;\n      if (Math.abs(A[piv * n + c]) < 1e-300) return false;\n      if (piv !== c) {\n        for (let j = 0; j < n; j++) { const tmp = A[c * n + j]; A[c * n + j] = A[piv * n + j]; A[piv * n + j] = tmp; }\n        const tr = r[c]; r[c] = r[piv]; r[piv] = tr;\n      }\n      for (let i = c + 1; i < n; i++) {\n        const m = A[i * n + c] / A[c * n + c];\n        if (m === 0) continue;\n        for (let j = c; j < n; j++) A[i * n + j] -= m * A[c * n + j];\n        r[i] -= m * r[c];\n      }\n    }\n    for (let i = n - 1; i >= 0; i--) {\n      let s = r[i];\n      for (let j = i + 1; j < n; j++) s -= A[i * n + j] * r[j];\n      r[i] = s / A[i * n + i];\n    }\n    return true;\n  }\n\n  /* Eigenvalues of a real n x n matrix (row-major): reduction to upper\n     Hessenberg form by Gaussian similarity transforms, then the shifted QR\n     iteration of the EISPACK routine hqr (Wilkinson and Reinsch, Handbook\n     for Automatic Computation, vol. 2, 1971). Returns [{re, im}]. */\n  function eigenvalues(M, n) {\n    const a = [];\n    for (let i = 0; i < n; i++) { a.push([]); for (let j = 0; j < n; j++) a[i].push(M[i * n + j]); }\n    // elmhes\n    for (let m = 1; m < n - 1; m++) {\n      let x = 0, i = m;\n      for (let j = m; j < n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }\n      if (i !== m) {\n        for (let j = m - 1; j < n; j++) { const t = a[i][j]; a[i][j] = a[m][j]; a[m][j] = t; }\n        for (let j = 0; j < n; j++) { const t = a[j][i]; a[j][i] = a[j][m]; a[j][m] = t; }\n      }\n      if (x !== 0) {\n        for (i = m + 1; i < n; i++) {\n          let y = a[i][m - 1];\n          if (y !== 0) {\n            y /= x; a[i][m - 1] = y;\n            for (let j = m; j < n; j++) a[i][j] -= y * a[m][j];\n            for (let j = 0; j < n; j++) a[j][m] += y * a[j][i];\n          }\n        }\n      }\n    }\n    for (let i = 2; i < n; i++) for (let j = 0; j < i - 1; j++) a[i][j] = 0;\n    // hqr\n    const wr = new Array(n).fill(0), wi = new Array(n).fill(0);\n    let anorm = 0;\n    for (let i = 0; i < n; i++) for (let j = Math.max(i - 1, 0); j < n; j++) anorm += Math.abs(a[i][j]);\n    let nn = n - 1, t = 0;\n    while (nn >= 0) {\n      let its = 0, l;\n      do {\n        for (l = nn; l >= 1; l--) {\n          const s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);\n          if (Math.abs(a[l][l - 1]) + (s === 0 ? anorm : s) === (s === 0 ? anorm : s)) { a[l][l - 1] = 0; break; }\n        }\n        const x = a[nn][nn];\n        if (l === nn) { wr[nn] = x + t; wi[nn--] = 0; }\n        else {\n          const y = a[nn - 1][nn - 1], w = a[nn][nn - 1] * a[nn - 1][nn];\n          if (l === nn - 1) {\n            const p = 0.5 * (y - x), q = p * p + w, z = Math.sqrt(Math.abs(q));\n            const xx = x + t;\n            if (q >= 0) {\n              const zz = p + (p >= 0 ? Math.abs(z) : -Math.abs(z));\n              wr[nn - 1] = wr[nn] = xx + zz;\n              if (zz) wr[nn] = xx - w / zz;\n              wi[nn - 1] = wi[nn] = 0;\n            } else {\n              wr[nn - 1] = wr[nn] = xx + p;\n              wi[nn - 1] = -(wi[nn] = z);\n            }\n            nn -= 2;\n          } else {\n            if (its === 60) throw new Error(\"Eigenvalue iteration did not converge\");\n            let xs = x, ys = y, ws = w;\n            if (its === 10 || its === 20) {\n              t += xs;\n              for (let i = 0; i <= nn; i++) a[i][i] -= xs;\n              const s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);\n              ys = xs = 0.75 * s; ws = -0.4375 * s * s;\n            }\n            ++its;\n            let m, p, q, r, z;\n            for (m = nn - 2; m >= l; m--) {\n              z = a[m][m];\n              r = xs - z; const s0 = ys - z;\n              p = (r * s0 - ws) / a[m + 1][m] + a[m][m + 1];\n              q = a[m + 1][m + 1] - z - r - s0;\n              r = a[m + 2][m + 1];\n              const s = Math.abs(p) + Math.abs(q) + Math.abs(r);\n              p /= s; q /= s; r /= s;\n              if (m === l) break;\n              const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));\n              const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));\n              if (u + v === v) break;\n            }\n            for (let i = m + 2; i <= nn; i++) { a[i][i - 2] = 0; if (i !== m + 2) a[i][i - 3] = 0; }\n            for (let k = m; k <= nn - 1; k++) {\n              if (k !== m) {\n                p = a[k][k - 1]; q = a[k + 1][k - 1]; r = 0;\n                if (k !== nn - 1) r = a[k + 2][k - 1];\n                xs = Math.abs(p) + Math.abs(q) + Math.abs(r);\n                if (xs !== 0) { p /= xs; q /= xs; r /= xs; }\n              }\n              const s = (p >= 0 ? 1 : -1) * Math.sqrt(p * p + q * q + r * r);\n              if (s !== 0) {\n                if (k === m) { if (l !== m) a[k][k - 1] = -a[k][k - 1]; } else a[k][k - 1] = -s * xs;\n                p += s; xs = p / s; ys = q / s; z = r / s; q /= p; r /= p;\n                for (let j = k; j <= nn; j++) {\n                  p = a[k][j] + q * a[k + 1][j];\n                  if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; }\n                  a[k + 1][j] -= p * ys; a[k][j] -= p * xs;\n                }\n                const mmin = nn < k + 3 ? nn : k + 3;\n                for (let i = l; i <= mmin; i++) {\n                  p = xs * a[i][k] + ys * a[i][k + 1];\n                  if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; }\n                  a[i][k + 1] -= p * q; a[i][k] -= p;\n                }\n              }\n            }\n          }\n        }\n      } while (l < nn - 1);\n    }\n    const out = [];\n    for (let i = 0; i < n; i++) out.push({ re: wr[i], im: wi[i] });\n    return out;\n  }\n\n  /* Equilibria (flows) or fixed points (maps) inside a box. Seeds are a\n     Halton sequence over the box plus any given points; each is refined by\n     damped Newton iteration; solutions closer than 1e-6 of the box size\n     are merged. Stability: all Re(lambda) < 0 for flows, all |lambda| < 1\n     for maps; the tolerance on the boundary is 1e-7. */\n  function findEquilibria(sys, p, box, opts) {\n    opts = opts || {};\n    const n = sys.vars.length, t = opts.t || 0, nSeeds = opts.seeds || 60;\n    const scale = box.map(function (r) { return Math.max(Math.abs(r[1] - r[0]), 1e-12); });\n    const seeds = (opts.extra || []).slice();\n    const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];\n    for (let s = 1; s <= nSeeds; s++) {\n      seeds.push(box.map(function (r, i) {\n        let f = 1, v = 0, k = s;\n        const b = primes[i % primes.length];\n        while (k > 0) { f /= b; v += f * (k % b); k = Math.floor(k / b); }\n        return r[0] + v * (r[1] - r[0]);\n      }));\n    }\n    const found = [];\n    const G = new Float64Array(n), G2 = new Float64Array(n), x = new Float64Array(n), y = new Float64Array(n);\n    const norm = function (v) { let s = 0; for (let i = 0; i < n; i++) s += (v[i] / scale[i]) * (v[i] / scale[i]); return Math.sqrt(s); };\n    seeds.forEach(function (s0) {\n      x.set(s0);\n      let ok = false;\n      for (let it = 0; it < 60; it++) {\n        residual(sys, t, x, p, G);\n        if (!G.every(isFinite)) return;\n        const r0 = norm(G);\n        if (r0 < 1e-11) { ok = true; break; }\n        const J = jacobian(sys, t, x, p);\n        if (sys.time === \"discrete\") for (let i = 0; i < n; i++) J[i * n + i] -= 1;\n        const d = Float64Array.from(G);\n        if (!solve(J, d, n)) return;\n        let lam = 1;\n        for (let ls = 0; ls < 20; ls++) {\n          for (let i = 0; i < n; i++) y[i] = x[i] - lam * d[i];\n          residual(sys, t, y, p, G2);\n          if (G2.every(isFinite) && norm(G2) < (1 - 1e-4 * lam) * r0) break;\n          lam *= 0.5;\n        }\n        x.set(y);\n        if (norm(d) * lam < 1e-13) { residual(sys, t, x, p, G); ok = norm(G) < 1e-8; break; }\n      }\n      if (!ok) return;\n      for (let i = 0; i < n; i++) {\n        const m = 0.5 * scale[i];\n        if (x[i] < box[i][0] - m || x[i] > box[i][1] + m) return;\n      }\n      if (found.some(function (e) { let d = 0; for (let i = 0; i < n; i++) d = Math.max(d, Math.abs(e.x[i] - x[i]) / scale[i]); return d < 1e-6; })) return;\n      found.push({ x: Array.from(x) });\n    });\n    found.forEach(function (e) {\n      const J = jacobian(sys, t, Float64Array.from(e.x), p);\n      e.eig = eigenvalues(J, n);\n      if (sys.time === \"discrete\") {\n        const rho = Math.max.apply(null, e.eig.map(function (l) { return Math.hypot(l.re, l.im); }));\n        e.stable = rho < 1 - 1e-7;\n        e.type = rho < 1 - 1e-7 ? \"stable\" : rho > 1 + 1e-7 ? (e.eig.some(function (l) { return Math.hypot(l.re, l.im) < 1 - 1e-7; }) ? \"saddle\" : \"unstable\") : \"marginal\";\n      } else {\n        const re = e.eig.map(function (l) { return l.re; });\n        const mx = Math.max.apply(null, re), mn = Math.min.apply(null, re);\n        const osc = e.eig.some(function (l) { return Math.abs(l.im) > 1e-9; });\n        e.stable = mx < -1e-7;\n        if (mx < -1e-7) e.type = osc ? \"stable focus\" : \"stable node\";\n        else if (mn > 1e-7) e.type = osc ? \"unstable focus\" : \"unstable node\";\n        else if (mn < -1e-7 && mx > 1e-7) e.type = osc ? \"saddle focus\" : \"saddle\";\n        else e.type = osc ? \"centre\" : \"non-hyperbolic\";\n      }\n    });\n    return found;\n  }\n\n  DF.residual = residual;\n  DF.jacobian = jacobian;\n  DF.solveLinear = solve;\n  DF.eigenvalues = eigenvalues;\n  DF.findEquilibria = findEquilibria;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/render/style.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Themes and palettes. Themes set the background and the ink of axes and\n   text; palettes set the colours of trajectories. The RElab palettes are\n   the brand palettes of the lab (indigo #170C3A and #2a1766, orange #EE6A24,\n   amber #FB9E07, pink #A52C60); blackboard follows the ink set used for dark\n   slides; the continuous ramps are sampled from the published colour maps. */\n(function (DF) {\n  \"use strict\";\n\n  const THEMES = {\n    \"relab-night\": { label: \"RElab night\", bg: [\"radial\", \"#2a1766\", \"#170C3A\", \"#0b0620\"], ink: \"#eceaf4\", muted: \"rgba(236,234,244,0.62)\", grid: \"rgba(236,234,244,0.12)\", blend: \"lighter\", dark: true },\n    \"blackboard\": { label: \"Blackboard\", bg: [\"solid\", \"#000000\"], ink: \"#ebebeb\", muted: \"#bfbfbf\", grid: \"#404040\", blend: \"lighter\", dark: true },\n    \"deep-sea\": { label: \"Deep sea\", bg: [\"radial\", \"#0f2a3d\", \"#07131f\", \"#02070c\"], ink: \"#e6f0f5\", muted: \"rgba(230,240,245,0.6)\", grid: \"rgba(230,240,245,0.12)\", blend: \"lighter\", dark: true },\n    \"graphite\": { label: \"Graphite\", bg: [\"solid\", \"#1b1d22\"], ink: \"#e6e6e6\", muted: \"#a8a8a8\", grid: \"#3a3d44\", blend: \"lighter\", dark: true },\n    \"paper\": { label: \"Paper\", bg: [\"solid\", \"#fffdf8\"], ink: \"#212529\", muted: \"#687078\", grid: \"#e4e1d8\", blend: \"source-over\", dark: false },\n    \"white\": { label: \"White\", bg: [\"solid\", \"#ffffff\"], ink: \"#262626\", muted: \"#666666\", grid: \"#ebebeb\", blend: \"source-over\", dark: false },\n    \"transparent\": { label: \"Transparent\", bg: [\"none\"], ink: \"#888888\", muted: \"#888888\", grid: \"rgba(128,128,128,0.2)\", blend: \"source-over\", dark: true }\n  };\n\n  // Discrete palettes: one colour per variable, member or branch.\n  const PALETTES = {\n    \"relab\": { label: \"RElab\", colors: [\"#EE6A24\", \"#FB9E07\", \"#A52C60\", \"#CF4446\", \"#F6D645\", \"#764BA2\", \"#3093CF\", \"#eceaf4\"] },\n    \"relab-qualitative\": { label: \"RElab qualitative\", colors: [\"#FB9E07\", \"#4777ef\", \"#009E73\", \"#764BA2\", \"#E85D04\", \"#00B4D8\", \"#D55E00\", \"#CC79A7\"] },\n    \"blackboard\": { label: \"Blackboard ink\", colors: [\"#5ec5ff\", \"#ff8c42\", \"#ff6b6b\", \"#3ddc97\", \"#a8d05b\", \"#ffd166\", \"#9d8df1\", \"#ff6ec7\", \"#d4a373\", \"#c0c0c0\"] },\n    \"tableau\": { label: \"Tableau 10\", colors: [\"#4e79a7\", \"#f28e2c\", \"#e15759\", \"#76b7b2\", \"#59a14f\", \"#edc949\", \"#af7aa1\", \"#ff9da7\", \"#9c755f\", \"#bab0ab\"] },\n    \"okabe-ito\": { label: \"Okabe-Ito\", colors: [\"#E69F00\", \"#56B4E9\", \"#009E73\", \"#F0E442\", \"#0072B2\", \"#D55E00\", \"#CC79A7\", \"#000000\"] },\n    \"ink\": { label: \"Ink\", colors: [\"#170C3A\", \"#2a1766\", \"#A52C60\", \"#EE6A24\", \"#687078\"] },\n    \"mono-amber\": { label: \"Amber\", colors: [\"#FB9E07\", \"#f5b04a\", \"#EE6A24\", \"#ffd699\"] },\n    \"mono-ice\": { label: \"Ice\", colors: [\"#9fd8ff\", \"#5ec5ff\", \"#d9f1ff\", \"#3093CF\"] }\n  };\n\n  // Continuous ramps, as colour stops sampled evenly on [0, 1].\n  const RAMPS = {\n    \"relab-fire\": { label: \"RElab fire\", stops: [\"#170C3A\", \"#2a1766\", \"#A52C60\", \"#EE6A24\", \"#FB9E07\", \"#F6D645\"] },\n    \"relab-sequential\": { label: \"RElab sequential\", stops: [\"#F5F0E6\", \"#FFE5B4\", \"#FFD699\", \"#FBC66A\", \"#FB9E07\", \"#E88507\", \"#D46B07\", \"#B85507\", \"#9C4007\"] },\n    \"relab-diverging\": { label: \"RElab diverging\", stops: [\"#170C3A\", \"#2a1766\", \"#764BA2\", \"#B8A9C9\", \"#F5F0E6\", \"#FFD699\", \"#FB9E07\", \"#E85D04\", \"#C7380B\"] },\n    \"blackboard\": { label: \"Blackboard ramp\", stops: [\"#5ec5ff\", \"#ff6ec7\", \"#ffd166\"] },\n    \"viridis\": { label: \"Viridis\", stops: [\"#440154\", \"#482878\", \"#3e4989\", \"#31688e\", \"#26828e\", \"#1f9e89\", \"#35b779\", \"#6ece58\", \"#b5de2b\", \"#fde725\"] },\n    \"magma\": { label: \"Magma\", stops: [\"#000004\", \"#1c1044\", \"#4f127b\", \"#812581\", \"#b5367a\", \"#e55064\", \"#fb8761\", \"#fec287\", \"#fcfdbf\"] },\n    \"inferno\": { label: \"Inferno\", stops: [\"#000004\", \"#1f0c48\", \"#550f6d\", \"#88226a\", \"#ba3655\", \"#e35933\", \"#f98e09\", \"#f8c932\", \"#fcffa4\"] },\n    \"mako\": { label: \"Mako\", stops: [\"#0b0405\", \"#2b1c35\", \"#3e356b\", \"#3b5698\", \"#357ba3\", \"#38a0ab\", \"#4bc4ad\", \"#8ae1b9\", \"#def5e5\"] },\n    \"sand-red\": { label: \"Sand red\", stops: [\"#ffffff\", \"#fbe3c3\", \"#f3b27a\", \"#e0703f\", \"#b8321f\", \"#8f0b12\"] }\n  };\n\n  function hexToRgb(hex) {\n    const h = hex.replace(\"#\", \"\");\n    const v = parseInt(h.length === 3 ? h.split(\"\").map(function (c) { return c + c; }).join(\"\") : h, 16);\n    return [(v >> 16) & 255, (v >> 8) & 255, v & 255];\n  }\n  const rgbCache = {};\n  function rgb(hex) { return rgbCache[hex] || (rgbCache[hex] = hexToRgb(hex)); }\n\n  // Colour of u in [0, 1] on a ramp, as an [r, g, b] triple.\n  function rampRGB(name, u) {\n    const st = (RAMPS[name] || RAMPS[\"relab-fire\"]).stops;\n    u = Math.min(1, Math.max(0, isFinite(u) ? u : 0)) * (st.length - 1);\n    const i = Math.min(st.length - 2, Math.floor(u)), f = u - i;\n    const a = rgb(st[i]), b = rgb(st[i + 1]);\n    return [Math.round(a[0] + f * (b[0] - a[0])), Math.round(a[1] + f * (b[1] - a[1])), Math.round(a[2] + f * (b[2] - a[2]))];\n  }\n\n  // Paint the background of a theme onto a 2D context.\n  function paintBackground(ctx, w, h, theme) {\n    const th = THEMES[theme] || THEMES[\"relab-night\"], bg = th.bg;\n    if (bg[0] === \"none\") { ctx.clearRect(0, 0, w, h); return; }\n    if (bg[0] === \"solid\") { ctx.fillStyle = bg[1]; ctx.fillRect(0, 0, w, h); return; }\n    const g = ctx.createRadialGradient(w * 0.5, h * 0.5, 0, w * 0.5, h * 0.5, Math.hypot(w, h) * 0.6);\n    g.addColorStop(0, bg[1]); g.addColorStop(0.55, bg[2]); g.addColorStop(1, bg[3]);\n    ctx.fillStyle = g; ctx.fillRect(0, 0, w, h);\n  }\n  function backgroundCSS(theme) {\n    const bg = (THEMES[theme] || THEMES[\"relab-night\"]).bg;\n    if (bg[0] === \"none\") return \"transparent\";\n    if (bg[0] === \"solid\") return bg[1];\n    return \"radial-gradient(ellipse at 50% 50%, \" + bg[1] + \" 0%, \" + bg[2] + \" 55%, \" + bg[3] + \" 100%)\";\n  }\n\n  DF.THEMES = THEMES;\n  DF.PALETTES = PALETTES;\n  DF.RAMPS = RAMPS;\n  DF.rgb = rgb;\n  DF.rampRGB = rampRGB;\n  DF.paintBackground = paintBackground;\n  DF.backgroundCSS = backgroundCSS;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/render/frame.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Coordinates: automatic axis ranges, the 2D and rotating 3D projections\n   from state space to the canvas, and axis drawing. */\n(function (DF) {\n  \"use strict\";\n\n  /* Axis ranges from a probe run: several members from the initial\n     condition, a transient, then the 1st and 99th percentiles of each\n     variable with a 6 percent margin. Declared `range` lines win. */\n  function autoRanges(sys, params, init, opts) {\n    opts = opts || {};\n    const dim = sys.vars.length, dt = opts.dt || 0.01, steps = opts.steps || 4000;\n    const sim = new DF.Simulator(sys, { dt: dt, params: params, init: init, n: 8, spread: opts.spread || 0.1, seed: 12345, perturbations: opts.perturbations || [] });\n    const samples = sys.vars.map(function () { return []; });\n    const transient = Math.floor(steps / 4);\n    for (let s = 0; s < steps; s++) {\n      sim.step();\n      if (s < transient || s % 2) continue;\n      for (let k = 0; k < sim.n; k++) {\n        if (!sim.alive[k]) continue;\n        for (let i = 0; i < dim; i++) samples[i].push(sim.X[k * dim + i]);\n      }\n    }\n    const out = {};\n    sys.vars.forEach(function (v, i) {\n      if (sys.ranges[v]) { out[v] = sys.ranges[v].slice(); return; }\n      const a = samples[i].filter(isFinite).sort(function (u, w) { return u - w; });\n      let lo, hi;\n      if (!a.length) { lo = -1; hi = 1; }\n      else { lo = a[Math.floor(0.01 * (a.length - 1))]; hi = a[Math.ceil(0.99 * (a.length - 1))]; }\n      if (!(hi > lo)) { const c = isFinite(lo) ? lo : 0, w = Math.max(Math.abs(c) * 0.5, 0.5); lo = c - w; hi = c + w; }\n      const m = 0.06 * (hi - lo);\n      out[v] = [lo - m, hi + m];\n    });\n    return out;\n  }\n\n  /* Camera. axes holds two or three variable indices; ranges maps each\n     axis to [lo, hi]. In 3D the box is centred, rotated by azimuth about the\n     vertical axis and tilted by elevation, and projected orthographically. */\n  function Camera(axes, ranges, opts) {\n    opts = opts || {};\n    this.axes = axes; this.ranges = ranges;\n    this.azim = opts.azim === undefined ? 0.6 : opts.azim;\n    this.elev = opts.elev === undefined ? 0.35 : opts.elev;\n    this.zoom = opts.zoom || 1;\n    this.pad = opts.pad === undefined ? 0.08 : opts.pad;\n    this.upAxis = opts.upAxis === undefined ? 2 : opts.upAxis; // in 3D, which of the three is vertical\n    this.simplex = !!opts.simplex && axes.length === 3;           // barycentric triangle for three shares\n    this.w = 1; this.h = 1;\n    this.inset = { l: 0, r: 0, t: 0, b: 0 };\n  }\n  Camera.prototype.resize = function (w, h, inset) { this.w = w; this.h = h; if (inset) this.inset = inset; this.tri = null; };\n  Camera.prototype.is3D = function () { return this.axes.length === 3 && !this.simplex; };\n  // Vertices of the simplex triangle: first axis on top, then lower left, lower right.\n  Camera.prototype.triangle = function () {\n    const b = this.plotBox(), side = Math.min(b.w * (1 - 2 * this.pad), b.h * (1 - 2 * this.pad) * 2 / Math.sqrt(3)) * this.zoom;\n    const cx = b.x + b.w / 2, cy = b.y + b.h / 2 + side * Math.sqrt(3) / 12, hgt = side * Math.sqrt(3) / 2;\n    return [[cx, cy - hgt * 2 / 3], [cx - side / 2, cy + hgt / 3], [cx + side / 2, cy + hgt / 3]];\n  };\n  Camera.prototype.plotBox = function () {\n    const i = this.inset;\n    return { x: i.l, y: i.t, w: this.w - i.l - i.r, h: this.h - i.t - i.b };\n  };\n  // Normalised coordinate in [-1, 1] of value v on axis a.\n  Camera.prototype.norm = function (a, v) {\n    const r = this.ranges[a];\n    return 2 * (v - r[0]) / (r[1] - r[0]) - 1;\n  };\n  Camera.prototype.project = function (x, out) {\n    out = out || [0, 0];\n    const b = this.plotBox();\n    if (this.simplex) {\n      const V = this.tri || (this.tri = this.triangle());\n      const a = Math.max(0, x[this.axes[0]]), c = Math.max(0, x[this.axes[1]]), d = Math.max(0, x[this.axes[2]]), s = a + c + d || 1;\n      out[0] = (a * V[0][0] + c * V[1][0] + d * V[2][0]) / s;\n      out[1] = (a * V[0][1] + c * V[1][1] + d * V[2][1]) / s;\n      return out;\n    }\n    if (!this.is3D()) {\n      const u = (x[this.axes[0]] - this.ranges[0][0]) / (this.ranges[0][1] - this.ranges[0][0]);\n      const v = (x[this.axes[1]] - this.ranges[1][0]) / (this.ranges[1][1] - this.ranges[1][0]);\n      const px = this.pad * b.w, py = this.pad * b.h;\n      const cx = b.x + b.w / 2, cy = b.y + b.h / 2;\n      out[0] = cx + (u - 0.5) * (b.w - 2 * px) * this.zoom;\n      out[1] = cy - (v - 0.5) * (b.h - 2 * py) * this.zoom;\n      return out;\n    }\n    // 3D: q = (horizontal a, horizontal b, vertical c) in [-1, 1]^3.\n    const up = this.upAxis, hA = up === 0 ? 1 : 0, hB = up === 2 ? 1 : 2;\n    const qa = this.norm(hA, x[this.axes[hA]]), qb = this.norm(hB, x[this.axes[hB]]), qc = this.norm(up, x[this.axes[up]]);\n    const ca = Math.cos(this.azim), sa = Math.sin(this.azim), ce = Math.cos(this.elev), se = Math.sin(this.elev);\n    const X = ca * qa - sa * qb, Y0 = sa * qa + ca * qb;\n    const Y = ce * qc - se * Y0;\n    const s = Math.min(b.w, b.h) * 0.5 * (1 - this.pad) / 1.5 * this.zoom * 1.25;\n    out[0] = b.x + b.w / 2 + s * X;\n    out[1] = b.y + b.h / 2 - s * Y;\n    out.depth = ce * Y0 + se * qc;\n    return out;\n  };\n  // Inverse of the 2D projection (clicks); null in 3D.\n  Camera.prototype.unproject = function (px, py) {\n    if (this.is3D()) return null;\n    if (this.simplex) {\n      const V = this.tri || (this.tri = this.triangle());\n      const det = (V[1][1] - V[2][1]) * (V[0][0] - V[2][0]) + (V[2][0] - V[1][0]) * (V[0][1] - V[2][1]);\n      const w0 = ((V[1][1] - V[2][1]) * (px - V[2][0]) + (V[2][0] - V[1][0]) * (py - V[2][1])) / det;\n      const w1 = ((V[2][1] - V[0][1]) * (px - V[2][0]) + (V[0][0] - V[2][0]) * (py - V[2][1])) / det;\n      return [w0, w1, 1 - w0 - w1];\n    }\n    const b = this.plotBox();\n    const padx = this.pad * b.w, pady = this.pad * b.h, cx = b.x + b.w / 2, cy = b.y + b.h / 2;\n    const u = (px - cx) / ((b.w - 2 * padx) * this.zoom) + 0.5, v = -(py - cy) / ((b.h - 2 * pady) * this.zoom) + 0.5;\n    return [this.ranges[0][0] + u * (this.ranges[0][1] - this.ranges[0][0]), this.ranges[1][0] + v * (this.ranges[1][1] - this.ranges[1][0])];\n  };\n\n  // Tick values: 3 to 7 round numbers inside [lo, hi].\n  function ticks(lo, hi, target) {\n    target = target || 5;\n    const span = hi - lo;\n    if (!(span > 0)) return [lo];\n    const raw = span / target, mag = Math.pow(10, Math.floor(Math.log10(raw)));\n    const step = [1, 2, 2.5, 5, 10].map(function (m) { return m * mag; }).find(function (s) { return span / s <= target + 1; }) || 10 * mag;\n    const out = [];\n    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);\n    return out;\n  }\n  function fmt(v) {\n    const a = Math.abs(v);\n    if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(1).replace(\"e+\", \"e\");\n    return String(+v.toPrecision(4));\n  }\n\n  /* Axes for a 2D camera: frame, ticks and labels. Labels are variable or\n     parameter names; the fonts are those of the page. */\n  function drawAxes(ctx, cam, theme, labels, opts) {\n    opts = opts || {};\n    const th = DF.THEMES[theme] || DF.THEMES[\"relab-night\"];\n    const b = cam.plotBox();\n    const x0 = cam.project([cam.ranges[0][0], cam.ranges[1][0]].reduce(function (acc, v, j) { acc[cam.axes[j]] = v; return acc; }, []));\n    const x1 = cam.project([cam.ranges[0][1], cam.ranges[1][1]].reduce(function (acc, v, j) { acc[cam.axes[j]] = v; return acc; }, []));\n    const L = x0[0], R = x1[0], B = x0[1], T = x1[1];\n    ctx.save();\n    ctx.strokeStyle = th.grid; ctx.lineWidth = 1;\n    ctx.font = (opts.fontSize || 11) + \"px Jost, system-ui, sans-serif\";\n    ctx.fillStyle = th.muted;\n    const probe = [];\n    ticks(cam.ranges[0][0], cam.ranges[0][1]).forEach(function (v) {\n      probe[cam.axes[0]] = v; probe[cam.axes[1]] = cam.ranges[1][0];\n      const p = cam.project(probe);\n      if (opts.grid) { ctx.beginPath(); ctx.moveTo(p[0], B); ctx.lineTo(p[0], T); ctx.stroke(); }\n      ctx.beginPath(); ctx.moveTo(p[0], B); ctx.lineTo(p[0], B + 4); ctx.stroke();\n      ctx.textAlign = \"center\"; ctx.textBaseline = \"top\"; ctx.fillText(fmt(v), p[0], B + 6);\n    });\n    ticks(cam.ranges[1][0], cam.ranges[1][1]).forEach(function (v) {\n      probe[cam.axes[0]] = cam.ranges[0][0]; probe[cam.axes[1]] = v;\n      const p = cam.project(probe);\n      if (opts.grid) { ctx.beginPath(); ctx.moveTo(L, p[1]); ctx.lineTo(R, p[1]); ctx.stroke(); }\n      ctx.beginPath(); ctx.moveTo(L - 4, p[1]); ctx.lineTo(L, p[1]); ctx.stroke();\n      ctx.textAlign = \"right\"; ctx.textBaseline = \"middle\"; ctx.fillText(fmt(v), L - 7, p[1]);\n    });\n    ctx.strokeStyle = th.muted; ctx.globalAlpha = 0.6;\n    ctx.strokeRect(L, T, R - L, B - T);\n    ctx.globalAlpha = 1; ctx.fillStyle = th.ink;\n    ctx.font = \"italic \" + ((opts.fontSize || 11) + 2) + \"px 'TeX Gyre Pagella', Palatino, serif\";\n    ctx.textAlign = \"center\"; ctx.textBaseline = \"top\";\n    ctx.fillText(labels[0], (L + R) / 2, B + 22);\n    ctx.save(); ctx.translate(L - 40, (T + B) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = \"bottom\"; ctx.fillText(labels[1], 0, 0); ctx.restore();\n    ctx.restore();\n    return { L: L, R: R, T: T, B: B };\n  }\n\n  // Wireframe of the 3D box, faint, for orientation.\n  function drawBox3D(ctx, cam, theme) {\n    const th = DF.THEMES[theme] || DF.THEMES[\"relab-night\"];\n    const c = [0, 1], p = [];\n    const corner = function (i, j, k) { const x = []; x[cam.axes[0]] = cam.ranges[0][i]; x[cam.axes[1]] = cam.ranges[1][j]; x[cam.axes[2]] = cam.ranges[2][k]; return cam.project(x, [0, 0]); };\n    ctx.save(); ctx.strokeStyle = th.grid; ctx.lineWidth = 1;\n    c.forEach(function (i) { c.forEach(function (j) { c.forEach(function (k) { p.push([i, j, k]); }); }); });\n    p.forEach(function (a) {\n      p.forEach(function (b) {\n        const d = Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);\n        if (d !== 1 || a.join() > b.join()) return;\n        const u = corner(a[0], a[1], a[2]), v = corner(b[0], b[1], b[2]);\n        ctx.beginPath(); ctx.moveTo(u[0], u[1]); ctx.lineTo(v[0], v[1]); ctx.stroke();\n      });\n    });\n    ctx.restore();\n  }\n\n  DF.autoRanges = autoRanges;\n  DF.Camera = Camera;\n  DF.ticks = ticks;\n  DF.fmt = fmt;\n  DF.drawAxes = drawAxes;\n  DF.drawBox3D = drawBox3D;\n  DF.drawSimplex = function (ctx, cam, theme, labels) {\n    const th = DF.THEMES[theme] || DF.THEMES[\"relab-night\"], V = cam.triangle();\n    ctx.save(); ctx.strokeStyle = th.grid; ctx.lineWidth = 1;\n    ctx.beginPath(); ctx.moveTo(V[0][0], V[0][1]); ctx.lineTo(V[1][0], V[1][1]); ctx.lineTo(V[2][0], V[2][1]); ctx.closePath(); ctx.stroke();\n    ctx.fillStyle = th.muted; ctx.font = \"italic 13px 'TeX Gyre Pagella', Palatino, serif\"; ctx.textAlign = \"center\";\n    ctx.textBaseline = \"bottom\"; ctx.fillText(labels[0], V[0][0], V[0][1] - 6);\n    ctx.textBaseline = \"top\"; ctx.fillText(labels[1], V[1][0] - 8, V[1][1] + 6); ctx.fillText(labels[2], V[2][0] + 8, V[2][1] + 6);\n    ctx.restore();\n  };\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/render/views.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Views. Each view draws one kind of figure from a running Simulator on the\n   three layers of a Player: base (static, redrawn on demand), trail\n   (accumulating, faded each frame) and top (cleared each frame).\n\n     flow        particle ensemble with fading trails, 2D or rotating 3D\n     trajectory  a few long orbits with a gradient tail, 2D or rotating 3D\n     timeseries  scrolling time series of chosen variables\n     phase       vector field, nullclines, classified equilibria, orbits\n     sweep       slow parameter sweep over the equilibrium branches (hysteresis)\n     orbit       bifurcation diagram, built column by column\n     density     ensemble density: 2D heat map, or a time carpet in 1D\n     strobe      stroboscopic samples every period T, or a Poincare section\n     cobweb      cobweb diagram of a one-dimensional map\n\n   Every view implements init(P), frame(P), and optionally drawStatic(P),\n   pointer(P, kind, x, y, ev), svg(P) and legend(P). */\n(function (DF) {\n  \"use strict\";\n\n  const V = {};\n\n  // ------------------------------------------------------------ helpers\n  function colorFor(P, k, x, extra) {\n    const s = P.scene.style, pal = P.palette;\n    switch (s.colorBy) {\n      case \"member\": return pal[k % pal.length];\n      case \"dominant\": {\n        let best = 0;\n        for (let i = 1; i < x.length; i++) if (x[i] > x[best]) best = i;\n        return pal[best % pal.length];\n      }\n      case \"speed\": case \"var\": case \"age\": case \"time\": {\n        const u = Math.round(Math.min(1, Math.max(0, extra)) * 23) / 23;\n        const c = DF.rampRGB(s.ramp, u);\n        return \"rgb(\" + c[0] + \",\" + c[1] + \",\" + c[2] + \")\";\n      }\n      default: return pal[0];\n    }\n  }\n  function Buckets() { this.map = new Map(); }\n  Buckets.prototype.seg = function (c, x0, y0, x1, y1) {\n    let a = this.map.get(c); if (!a) { a = []; this.map.set(c, a); }\n    a.push(x0, y0, x1, y1);\n  };\n  Buckets.prototype.strokeAll = function (ctx, width, alpha) {\n    ctx.lineWidth = width; ctx.globalAlpha = alpha; ctx.lineCap = \"round\";\n    this.map.forEach(function (a, c) {\n      ctx.strokeStyle = c; ctx.beginPath();\n      for (let i = 0; i < a.length; i += 4) { ctx.moveTo(a[i], a[i + 1]); ctx.lineTo(a[i + 2], a[i + 3]); }\n      ctx.stroke();\n    });\n    ctx.globalAlpha = 1; this.map.clear();\n  };\n  Buckets.prototype.dotAll = function (ctx, size, alpha) {\n    ctx.globalAlpha = alpha;\n    this.map.forEach(function (a, c) {\n      ctx.fillStyle = c;\n      for (let i = 0; i < a.length; i += 4) ctx.fillRect(a[i] - size / 2, a[i + 1] - size / 2, size, size);\n    });\n    ctx.globalAlpha = 1; this.map.clear();\n  };\n  function speedOf(P, x) {\n    const d = P.tmpDx;\n    P.sys.f(P.sim.t, x, P.sim.p, d, function (i) { return x[i]; });\n    let s = 0;\n    for (let i = 0; i < d.length; i++) { const r = P.fullRanges[i]; const w = r[1] - r[0]; s += (d[i] / w) * (d[i] / w); }\n    return Math.sqrt(s);\n  }\n  function extraFor(P, k, x, age) {\n    const s = P.scene.style;\n    if (s.colorBy === \"speed\") { const v = speedOf(P, x); P.speedMax = Math.max(P.speedMax * 0.9995, v); return v / (P.speedMax || 1); }\n    if (s.colorBy === \"var\") { const i = Math.max(0, P.sys.vars.indexOf(s.colorVar)); const r = P.fullRanges[i]; return (x[i] - r[0]) / (r[1] - r[0]); }\n    if (s.colorBy === \"age\") return age;\n    return 0;\n  }\n  function randomInBox(P, out) {\n    const r = P.sim.rng;\n    for (let i = 0; i < P.sys.vars.length; i++) out[i] = r.range(P.fullRanges[i][0], P.fullRanges[i][1]);\n    return out;\n  }\n  function spawnState(P, out) {\n    const v = P.scene.view, r = P.sim.rng, init = P.sim.init;\n    const mode = v.spawn === \"mixed\" ? (r.uniform() < (v.spawnMix === undefined ? 0.25 : v.spawnMix) ? \"init\" : \"box\") : v.spawn;\n    if (mode === \"init\") {\n      for (let i = 0; i < out.length; i++) { const w = P.fullRanges[i][1] - P.fullRanges[i][0]; out[i] = init[i] + (P.scene.spread || 0.02) * w * r.normal(); }\n      return out;\n    }\n    return randomInBox(P, out);\n  }\n  function outOfView(P, x) {\n    for (let i = 0; i < x.length; i++) {\n      const r = P.fullRanges[i], w = r[1] - r[0];\n      if (x[i] < r[0] - 2 * w || x[i] > r[1] + 2 * w) return true;\n    }\n    return false;\n  }\n  function svgHead(P) {\n    return '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' + P.w + '\" height=\"' + P.h + '\" viewBox=\"0 0 ' + P.w + \" \" + P.h + '\">';\n  }\n  function polyline(pts, color, width, alpha) {\n    if (pts.length < 4) return \"\";\n    let d = \"M\" + pts[0].toFixed(2) + \" \" + pts[1].toFixed(2);\n    for (let i = 2; i < pts.length; i += 2) d += \"L\" + pts[i].toFixed(2) + \" \" + pts[i + 1].toFixed(2);\n    return '<path d=\"' + d + '\" fill=\"none\" stroke=\"' + color + '\" stroke-width=\"' + width + '\" stroke-opacity=\"' + alpha + '\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>';\n  }\n  function rotate3D(P) {\n    if (P.cam.is3D() && !P.dragging) P.cam.azim += (P.scene.view.rotate || 0) * 0.01;\n  }\n\n  // ------------------------------------------------------------ flow\n  V.flow = {\n    label: \"Flow (particles)\",\n    init: function (P) {\n      const n = P.sim.n, st = new Float64Array(P.sys.vars.length);\n      this.age = new Float32Array(n); this.life = new Float32Array(n);\n      this.prev = new Float32Array(2 * n); this.has = new Uint8Array(n);\n      this.buck = new Buckets(); this.pt = [0, 0];\n      for (let k = 0; k < n; k++) {\n        spawnState(P, st); P.sim.setMember(k, st);\n        this.life[k] = this.newLife(P); this.age[k] = P.sim.rng.uniform() * this.life[k];\n      }\n      // Warm-up so that the first frame already shows the flow.\n      const warm = P.scene.view.warmup === undefined ? 60 : P.scene.view.warmup;\n      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();\n    },\n    newLife: function (P) {\n      const L = P.scene.view.life;\n      if (!L || L[1] === Infinity || L === \"inf\") return Infinity;\n      return L[0] + P.sim.rng.uniform() * (L[1] - L[0]);\n    },\n    frame: function (P) {\n      const sim = P.sim, n = sim.n, dim = P.sys.vars.length, cam = P.cam, st = P.scene.style;\n      rotate3D(P);\n      P.fade(st.fade);\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();\n      const x = P.tmpX, discrete = P.sys.time === \"discrete\", st2 = new Float64Array(dim);\n      for (let k = 0; k < n; k++) {\n        this.age[k]++;\n        if (!sim.alive[k] || this.age[k] > this.life[k]) {\n          spawnState(P, st2); sim.setMember(k, st2); this.age[k] = 0; this.life[k] = this.newLife(P); this.has[k] = 0;\n          continue;\n        }\n        sim.member(k, x);\n        if (outOfView(P, x)) { this.age[k] = this.life[k] + 1; continue; }\n        const q = cam.project(x, this.pt);\n        const col = colorFor(P, k, x, extraFor(P, k, x, this.age[k] / (isFinite(this.life[k]) ? this.life[k] : 1)));\n        if (discrete) this.buck.seg(col, q[0], q[1], q[0], q[1]);\n        else if (this.has[k]) this.buck.seg(col, this.prev[2 * k], this.prev[2 * k + 1], q[0], q[1]);\n        this.prev[2 * k] = q[0]; this.prev[2 * k + 1] = q[1]; this.has[k] = 1;\n      }\n      const ctx = P.ctx.trail;\n      ctx.globalCompositeOperation = P.theme.blend;\n      if (discrete) this.buck.dotAll(ctx, st.pointSize, st.alpha); else this.buck.strokeAll(ctx, st.lineWidth, st.alpha);\n      ctx.globalCompositeOperation = \"source-over\";\n    },\n    pointer: function (P, kind, px, py) {\n      if (kind !== \"down\" || P.cam.is3D()) return false;\n      const u = P.cam.unproject(px, py);\n      const s = Float64Array.from(P.sim.init), r = P.sim.rng;\n      if (P.cam.simplex) {\n        if (u.some(function (c) { return c <= 0.002; })) return false;\n        const tot = P.cam.axes.reduce(function (acc, i) { return acc + P.sim.init[i]; }, 0) || 1;\n        for (let j = 0; j < 80; j++) {\n          const k = Math.floor(r.uniform() * P.sim.n);\n          P.cam.axes.forEach(function (i, m) { s[i] = tot * Math.max(1e-6, u[m] + 0.005 * r.normal()); });\n          P.sim.setMember(k, s); this.age[k] = 0; this.has[k] = 0;\n        }\n        return true;\n      }\n      for (let j = 0; j < 80; j++) {\n        const k = Math.floor(r.uniform() * P.sim.n);\n        for (let i = 0; i < s.length; i++) s[i] = P.sim.init[i];\n        const w0 = P.cam.ranges[0][1] - P.cam.ranges[0][0], w1 = P.cam.ranges[1][1] - P.cam.ranges[1][0];\n        s[P.cam.axes[0]] = u[0] + 0.01 * w0 * r.normal(); s[P.cam.axes[1]] = u[1] + 0.01 * w1 * r.normal();\n        P.sim.setMember(k, s); this.age[k] = 0; this.has[k] = 0;\n      }\n      return true;\n    }\n  };\n\n  // ------------------------------------------------------------ trajectory\n  V.trajectory = {\n    label: \"Trajectory\",\n    init: function (P) {\n      const n = P.sim.n, dim = P.sys.vars.length, L = P.scene.view.tail || 2500;\n      this.L = L; this.buf = new Float64Array(n * L * dim); this.len = new Int32Array(n); this.head = new Int32Array(n);\n      const warm = P.scene.view.warmup === undefined ? 0 : P.scene.view.warmup;\n      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();\n    },\n    push: function (P) {\n      const n = P.sim.n, dim = P.sys.vars.length, L = this.L;\n      for (let k = 0; k < n; k++) {\n        if (!P.sim.alive[k]) continue;\n        const h = this.head[k], o = (k * L + h) * dim;\n        for (let i = 0; i < dim; i++) this.buf[o + i] = P.sim.X[k * dim + i];\n        this.head[k] = (h + 1) % L; this.len[k] = Math.min(this.len[k] + 1, L);\n      }\n    },\n    points: function (P, k) {\n      const dim = P.sys.vars.length, L = this.L, len = this.len[k], out = new Float32Array(2 * len), x = P.tmpX, q = [0, 0];\n      const start = (this.head[k] - len + L) % L;\n      for (let j = 0; j < len; j++) {\n        const o = (k * L + (start + j) % L) * dim;\n        for (let i = 0; i < dim; i++) x[i] = this.buf[o + i];\n        P.cam.project(x, q); out[2 * j] = q[0]; out[2 * j + 1] = q[1];\n      }\n      return out;\n    },\n    frame: function (P) {\n      rotate3D(P);\n      const every = P.scene.view.sampleEvery || 1;\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) { P.sim.step(); if (s % every === 0) this.push(P); }\n      const ctx = P.ctx.top, st = P.scene.style, n = P.sim.n, G = 28;\n      ctx.clearRect(0, 0, P.w, P.h);\n      if (P.cam.is3D() && P.scene.view.box) DF.drawBox3D(ctx, P.cam, P.scene.style.theme);\n      ctx.globalCompositeOperation = P.theme.blend; ctx.lineCap = \"round\"; ctx.lineJoin = \"round\";\n      const discrete = P.sys.time === \"discrete\";\n      for (let k = 0; k < n; k++) {\n        const pts = this.points(P, k), m = pts.length / 2;\n        if (m < 2) continue;\n        for (let g = 0; g < G; g++) {\n          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);\n          if (b <= a) continue;\n          const u = (g + 1) / G;\n          let col;\n          if (st.colorBy === \"time\" || st.colorBy === \"age\") { const c = DF.rampRGB(st.ramp, u); col = \"rgb(\" + c.join(\",\") + \")\"; }\n          else col = P.palette[(st.colorBy === \"member\" ? k : 0) % P.palette.length];\n          ctx.strokeStyle = col; ctx.fillStyle = col;\n          ctx.globalAlpha = st.alpha * (0.08 + 0.92 * Math.pow(u, 1.4));\n          if (discrete) { for (let j = a; j <= b; j++) ctx.fillRect(pts[2 * j] - st.pointSize / 2, pts[2 * j + 1] - st.pointSize / 2, st.pointSize, st.pointSize); continue; }\n          ctx.lineWidth = st.lineWidth;\n          ctx.beginPath(); ctx.moveTo(pts[2 * a], pts[2 * a + 1]);\n          for (let j = a + 1; j <= b; j++) ctx.lineTo(pts[2 * j], pts[2 * j + 1]);\n          ctx.stroke();\n        }\n        ctx.globalAlpha = 1;\n        ctx.fillStyle = P.palette[(st.colorBy === \"member\" ? k : 0) % P.palette.length];\n        ctx.beginPath(); ctx.arc(pts[pts.length - 2], pts[pts.length - 1], Math.max(2.5, st.lineWidth * 2), 0, 7); ctx.fill();\n      }\n      ctx.globalCompositeOperation = \"source-over\"; ctx.globalAlpha = 1;\n    },\n    pointer: function (P, kind, px, py) {\n      if (kind !== \"down\" || P.cam.is3D()) return false;\n      const u = P.cam.unproject(px, py), s = Float64Array.from(P.sim.init);\n      s[P.cam.axes[0]] = u[0]; s[P.cam.axes[1]] = u[1];\n      const k = (this.nextSeed = ((this.nextSeed || 0) + 1) % P.sim.n);\n      P.sim.setMember(k, s); this.len[k] = 0; this.head[k] = 0;\n      return true;\n    },\n    svg: function (P) {\n      let out = \"\";\n      const st = P.scene.style;\n      for (let k = 0; k < P.sim.n; k++) {\n        const pts = Array.from(this.points(P, k));\n        out += polyline(pts, P.palette[(st.colorBy === \"member\" ? k : 0) % P.palette.length], st.lineWidth, st.alpha);\n      }\n      return out;\n    }\n  };\n\n  // ------------------------------------------------------------ timeseries\n  V.timeseries = {\n    label: \"Time series\",\n    axes: true,\n    init: function (P) {\n      const vars = (P.scene.view.vars && P.scene.view.vars.length ? P.scene.view.vars : P.sys.vars.slice(0, 4)).filter(function (v) { return P.sys.vars.indexOf(v) >= 0; });\n      this.vi = vars.map(function (v) { return P.sys.vars.indexOf(v); });\n      this.members = Math.min(P.sim.n, P.scene.view.members || 1);\n      this.cap = 4000; this.T = new Float64Array(this.cap); this.Y = new Float64Array(this.cap * this.members * this.vi.length);\n      this.len = 0; this.head = 0;\n      let lo = Infinity, hi = -Infinity;\n      this.vi.forEach(function (i) { lo = Math.min(lo, P.fullRanges[i][0]); hi = Math.max(hi, P.fullRanges[i][1]); });\n      this.yr = P.scene.view.yRange || [lo, hi];\n    },\n    frame: function (P) {\n      const nv = this.vi.length, M = this.members;\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) {\n        P.sim.step();\n        const h = this.head;\n        this.T[h] = P.sim.t;\n        for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) this.Y[(h * M + m) * nv + j] = P.sim.X[m * P.sys.vars.length + this.vi[j]];\n        this.head = (h + 1) % this.cap; this.len = Math.min(this.len + 1, this.cap);\n      }\n      const span = P.scene.view.window || (P.sys.time === \"discrete\" ? 100 : 50);\n      const t0 = Math.max(0, P.sim.t - span), t1 = t0 + span;\n      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);\n      const cam = new DF.Camera([0, 1], [[t0, t1], this.yr], { pad: 0.02 });\n      cam.resize(P.w, P.h, P.inset);\n      DF.drawAxes(ctx, cam, P.scene.style.theme, [P.sys.time === \"discrete\" ? \"n\" : \"t\", this.vi.map(function (i) { return P.sys.vars[i]; }).join(\", \")], { grid: true });\n      const st = P.scene.style, q = [0, 0], x = [0, 0];\n      ctx.save(); const b = cam.plotBox(); ctx.beginPath(); ctx.rect(b.x, b.y, b.w, b.h); ctx.clip();\n      ctx.lineJoin = \"round\";\n      const start = (this.head - this.len + this.cap) % this.cap;\n      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) {\n        ctx.strokeStyle = P.palette[j % P.palette.length]; ctx.globalAlpha = M > 1 ? Math.max(0.25, st.alpha * 0.6) : st.alpha; ctx.lineWidth = st.lineWidth;\n        ctx.beginPath(); let started = false;\n        for (let r = 0; r < this.len; r++) {\n          const h = (start + r) % this.cap;\n          if (this.T[h] < t0) continue;\n          x[0] = this.T[h]; x[1] = this.Y[(h * M + m) * nv + j];\n          cam.project(x, q);\n          if (!started) { ctx.moveTo(q[0], q[1]); started = true; } else ctx.lineTo(q[0], q[1]);\n        }\n        ctx.stroke();\n      }\n      ctx.restore(); ctx.globalAlpha = 1;\n      this.cam = cam;\n    },\n    legend: function (P) { return this.vi.map(function (i, j) { return [P.sys.vars[i], P.palette[j % P.palette.length]]; }); },\n    svg: function (P) {\n      const nv = this.vi.length, M = this.members, cam = this.cam, t0 = cam.ranges[0][0];\n      let out = \"\";\n      const start = (this.head - this.len + this.cap) % this.cap, q = [0, 0], x = [0, 0];\n      for (let m = 0; m < M; m++) for (let j = 0; j < nv; j++) {\n        const pts = [];\n        for (let r = 0; r < this.len; r++) { const h = (start + r) % this.cap; if (this.T[h] < t0) continue; x[0] = this.T[h]; x[1] = this.Y[(h * M + m) * nv + j]; cam.project(x, q); pts.push(q[0], q[1]); }\n        out += polyline(pts, P.palette[j % P.palette.length], P.scene.style.lineWidth, P.scene.style.alpha);\n      }\n      return out;\n    }\n  };\n\n  // ------------------------------------------------------------ phase\n  // Marching squares on a grid of values g (nx by ny), level 0.\n  function contour(g, nx, ny, X, Y) {\n    const segs = [];\n    const lerp = function (a, b) { return a / (a - b); };\n    for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {\n      const a = g[j * nx + i], b = g[j * nx + i + 1], c = g[(j + 1) * nx + i + 1], d = g[(j + 1) * nx + i];\n      if (!isFinite(a + b + c + d)) continue;\n      const pts = [];\n      if ((a > 0) !== (b > 0)) pts.push([X(i + lerp(a, b)), Y(j)]);\n      if ((b > 0) !== (c > 0)) pts.push([X(i + 1), Y(j + lerp(b, c))]);\n      if ((c > 0) !== (d > 0)) pts.push([X(i + 1 - lerp(c, d)), Y(j + 1)]);\n      if ((d > 0) !== (a > 0)) pts.push([X(i), Y(j + 1 - lerp(d, a))]);\n      if (pts.length === 2) segs.push(pts[0], pts[1]);\n      else if (pts.length === 4) segs.push(pts[0], pts[1], pts[2], pts[3]);\n    }\n    return segs;\n  }\n\n  V.phase = {\n    label: \"Phase plane\",\n    axes: true,\n    init: function (P) {\n      this.trails = []; this.maxTrails = 40;\n      this.static = null;\n      const n0 = Math.min(P.sim.n, P.scene.view.seeds === undefined ? 6 : P.scene.view.seeds);\n      const r = P.sim.rng, s = new Float64Array(P.sys.vars.length);\n      for (let k = 0; k < n0; k++) {\n        for (let i = 0; i < s.length; i++) s[i] = P.sim.init[i];\n        if (k > 0) {\n          const a0 = P.cam.axes[0], a1 = P.cam.axes[1];\n          s[a0] = r.range(P.cam.ranges[0][0], P.cam.ranges[0][1]); s[a1] = r.range(P.cam.ranges[1][0], P.cam.ranges[1][1]);\n        }\n        this.addTrail(P, s);\n      }\n    },\n    addTrail: function (P, s) {\n      if (this.trails.length >= this.maxTrails) this.trails.shift();\n      this.trails.push({ x: Float64Array.from(s), pts: [], life: 0 });\n    },\n    drawStatic: function (P) {\n      const ctx = P.ctx.base, cam = P.cam, sys = P.sys, v = P.scene.view, st = P.scene.style;\n      const th = P.theme;\n      const frame = DF.drawAxes(ctx, cam, st.theme, [sys.vars[cam.axes[0]], sys.vars[cam.axes[1]]], { grid: false });\n      this.frameBox = frame;\n      const a0 = cam.axes[0], a1 = cam.axes[1], t = P.sim.t, p = P.sim.p;\n      const x = Float64Array.from(P.sim.init), d = new Float64Array(sys.vars.length);\n      const H = function (i) { return x[i]; };\n      const R0 = cam.ranges[0], R1 = cam.ranges[1];\n      ctx.save(); ctx.beginPath(); ctx.rect(frame.L, frame.T, frame.R - frame.L, frame.B - frame.T); ctx.clip();\n      // Vector field: arrows on a grid, length by log speed.\n      if (v.field !== \"none\") {\n        const nx = v.fieldDensity || 22, ny = Math.round(nx * (frame.B - frame.T) / (frame.R - frame.L));\n        const cellW = (frame.R - frame.L) / nx, cellH = (frame.B - frame.T) / ny;\n        const q = [0, 0];\n        let vmax = 0;\n        const vals = [];\n        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {\n          x[a0] = R0[0] + (i + 0.5) / nx * (R0[1] - R0[0]); x[a1] = R1[1] - (j + 0.5) / ny * (R1[1] - R1[0]);\n          sys.f(t, x, p, d, H);\n          const u = d[a0] / (R0[1] - R0[0]) * (frame.R - frame.L), w = -d[a1] / (R1[1] - R1[0]) * (frame.B - frame.T);\n          const m = Math.hypot(u, w); vmax = Math.max(vmax, isFinite(m) ? m : 0);\n          vals.push([frame.L + (i + 0.5) * cellW, frame.T + (j + 0.5) * cellH, u, w, m]);\n        }\n        ctx.lineWidth = 1;\n        vals.forEach(function (e) {\n          if (!(e[4] > 0) || !isFinite(e[4])) return;\n          const rel = Math.log1p(9 * e[4] / vmax) / Math.log(10);\n          const len = Math.min(cellW, cellH) * 0.42 * (0.35 + 0.65 * rel), ux = e[2] / e[4], uy = e[3] / e[4];\n          const c = DF.rampRGB(st.ramp, 0.25 + 0.75 * rel);\n          ctx.strokeStyle = th.dark ? \"rgba(\" + c.join(\",\") + \",0.55)\" : \"rgba(\" + c.join(\",\") + \",0.75)\";\n          const x0 = e[0] - ux * len, y0 = e[1] - uy * len, x1 = e[0] + ux * len, y1 = e[1] + uy * len;\n          ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1);\n          ctx.moveTo(x1, y1); ctx.lineTo(x1 - 4 * ux + 2.5 * uy, y1 - 4 * uy - 2.5 * ux);\n          ctx.moveTo(x1, y1); ctx.lineTo(x1 - 4 * ux - 2.5 * uy, y1 - 4 * uy + 2.5 * ux);\n          ctx.stroke();\n        });\n      }\n      // Nullclines f_a0 = 0 and f_a1 = 0 by marching squares.\n      if (v.nullclines !== false && P.sys.time === \"continuous\") {\n        const nx = 160, ny = 160, g0 = new Float64Array(nx * ny), g1 = new Float64Array(nx * ny);\n        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {\n          x[a0] = R0[0] + i / (nx - 1) * (R0[1] - R0[0]); x[a1] = R1[0] + j / (ny - 1) * (R1[1] - R1[0]);\n          sys.f(t, x, p, d, H);\n          g0[j * nx + i] = d[a0]; g1[j * nx + i] = d[a1];\n        }\n        const X = function (i) { return R0[0] + i / (nx - 1) * (R0[1] - R0[0]); }, Y = function (j) { return R1[0] + j / (ny - 1) * (R1[1] - R1[0]); };\n        const self = this; self.nullSegs = [];\n        [g0, g1].forEach(function (g, idx) {\n          const segs = contour(g, nx, ny, X, Y);\n          const col = P.palette[(idx + 1) % P.palette.length];\n          ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.setLineDash(idx ? [6, 4] : []); ctx.globalAlpha = 0.9;\n          ctx.beginPath();\n          const q = [0, 0], z = Float64Array.from(P.sim.init);\n          const pts = [];\n          for (let s = 0; s < segs.length; s += 2) {\n            z[a0] = segs[s][0]; z[a1] = segs[s][1]; cam.project(z, q); ctx.moveTo(q[0], q[1]); pts.push(q[0], q[1]);\n            z[a0] = segs[s + 1][0]; z[a1] = segs[s + 1][1]; cam.project(z, q); ctx.lineTo(q[0], q[1]); pts.push(q[0], q[1]);\n          }\n          ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;\n          self.nullSegs.push({ pts: pts, color: col, dash: idx ? \"6 4\" : \"\" });\n        });\n      }\n      ctx.restore();\n      // Equilibria in the plane (other variables held at their initial values).\n      this.equilibria = [];\n      if (v.equilibria !== false && P.sys.time === \"continuous\") {\n        const sub = { vars: [sys.vars[a0], sys.vars[a1]], time: \"continuous\", f: function (tt, y, pp, out) { x[a0] = y[0]; x[a1] = y[1]; sys.f(tt, x, pp, d, H); out[0] = d[a0]; out[1] = d[a1]; } };\n        try { this.equilibria = DF.findEquilibria(sub, p, [R0, R1], { t: t, seeds: 80 }); } catch (e) { this.equilibria = []; }\n        const q = [0, 0], z = Float64Array.from(P.sim.init);\n        this.equilibria.forEach(function (e) {\n          z[a0] = e.x[0]; z[a1] = e.x[1]; cam.project(z, q);\n          ctx.lineWidth = 2; ctx.strokeStyle = th.ink; ctx.fillStyle = th.ink;\n          ctx.beginPath(); ctx.arc(q[0], q[1], 5.5, 0, 7);\n          if (e.stable) ctx.fill();\n          else if (/saddle/.test(e.type)) { ctx.fillStyle = P.palette[0]; ctx.fill(); ctx.stroke(); }\n          else { ctx.fillStyle = th.dark ? \"#0b0620\" : \"#ffffff\"; ctx.fill(); ctx.stroke(); }\n        });\n      }\n    },\n    frame: function (P) {\n      const sys = P.sys, dim = sys.vars.length, h = P.sim.h, cam = P.cam, st = P.scene.style;\n      if (sys.usesTime || P.sim.perturbations.some(function (q) { return q.enabled !== false && DF.PARAM_PERTURBATIONS.indexOf(q.kind) >= 0; })) {\n        if ((P.frameCount % 20) === 0) P.redrawStatic();\n      }\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) P.sim.step();\n      const rk = this.rk || (this.rk = DF.makeRK4(dim)), q = [0, 0], tt = P.sim.t;\n      const maxPts = P.scene.view.tail || 1500;\n      const self = this;\n      this.trails.forEach(function (tr) {\n        for (let s = 0; s < P.scene.stepsPerFrame; s++) {\n          if (sys.time === \"discrete\") { sys.f(tt, tr.x, P.sim.p, P.tmpDx); tr.x.set(P.tmpDx); }\n          else rk(sys.f, tt, tr.x, P.sim.p, h, function (i) { return tr.x[i]; });\n          if (!tr.x.every(isFinite) || outOfView(P, tr.x)) { tr.dead = true; break; }\n          cam.project(tr.x, q); tr.pts.push(q[0], q[1]);\n        }\n        if (tr.pts.length > 2 * maxPts) tr.pts.splice(0, tr.pts.length - 2 * maxPts);\n      });\n      this.trails = this.trails.filter(function (tr) { return !tr.dead || tr.pts.length > 2; });\n      const ctx = P.ctx.top; ctx.clearRect(0, 0, P.w, P.h);\n      ctx.save();\n      if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }\n      this.trails.forEach(function (tr, k) {\n        const col = P.palette[(st.colorBy === \"member\" ? k : 0) % P.palette.length];\n        const pts = tr.pts, m = pts.length / 2;\n        if (m < 2) return;\n        ctx.strokeStyle = col; ctx.lineWidth = st.lineWidth; ctx.lineJoin = \"round\";\n        const G = 10;\n        for (let g = 0; g < G; g++) {\n          const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);\n          if (b <= a) continue;\n          ctx.globalAlpha = st.alpha * (0.15 + 0.85 * (g + 1) / G);\n          ctx.beginPath(); ctx.moveTo(pts[2 * a], pts[2 * a + 1]);\n          for (let j = a + 1; j <= b; j++) ctx.lineTo(pts[2 * j], pts[2 * j + 1]);\n          ctx.stroke();\n        }\n        ctx.globalAlpha = 1; ctx.fillStyle = col;\n        if (!tr.dead) { ctx.beginPath(); ctx.arc(pts[pts.length - 2], pts[pts.length - 1], 3.5, 0, 7); ctx.fill(); }\n      });\n      ctx.restore();\n      self.lastTrails = this.trails;\n    },\n    pointer: function (P, kind, px, py) {\n      if (kind !== \"down\") return false;\n      const u = P.cam.unproject(px, py), s = Float64Array.from(P.sim.init);\n      s[P.cam.axes[0]] = u[0]; s[P.cam.axes[1]] = u[1];\n      this.addTrail(P, s);\n      return true;\n    },\n    legend: function (P) {\n      const L = [[\"stable\", P.theme.ink, \"dot\"], [\"saddle\", P.palette[0], \"dot\"], [\"unstable\", P.theme.ink, \"ring\"]];\n      if (P.scene.view.nullclines !== false) L.unshift([P.sys.vars[P.cam.axes[0]] + \"-nullcline\", P.palette[1 % P.palette.length], \"line\"], [P.sys.vars[P.cam.axes[1]] + \"-nullcline\", P.palette[2 % P.palette.length], \"dash\"]);\n      return L;\n    },\n    svg: function (P) {\n      let out = \"\";\n      (this.nullSegs || []).forEach(function (nc) {\n        let d = \"\";\n        for (let i = 0; i < nc.pts.length; i += 4) d += \"M\" + nc.pts[i].toFixed(1) + \" \" + nc.pts[i + 1].toFixed(1) + \"L\" + nc.pts[i + 2].toFixed(1) + \" \" + nc.pts[i + 3].toFixed(1);\n        out += '<path d=\"' + d + '\" fill=\"none\" stroke=\"' + nc.color + '\" stroke-width=\"1.6\"' + (nc.dash ? ' stroke-dasharray=\"' + nc.dash + '\"' : \"\") + \"/>\";\n      });\n      (this.trails || []).forEach(function (tr, k) { out += polyline(tr.pts, P.palette[0], P.scene.style.lineWidth, P.scene.style.alpha); });\n      const q = [0, 0], z = Float64Array.from(P.sim.init), th = P.theme;\n      (this.equilibria || []).forEach(function (e) {\n        z[P.cam.axes[0]] = e.x[0]; z[P.cam.axes[1]] = e.x[1]; P.cam.project(z, q);\n        out += '<circle cx=\"' + q[0].toFixed(1) + '\" cy=\"' + q[1].toFixed(1) + '\" r=\"5.5\" fill=\"' + (e.stable ? th.ink : \"none\") + '\" stroke=\"' + th.ink + '\" stroke-width=\"2\"/>';\n      });\n      return out;\n    }\n  };\n\n  // ------------------------------------------------------------ sweep\n  /* Branches of equilibria (or fixed points) against one parameter, found\n     by Newton iteration at each of `cols` parameter values, then a single\n     state driven by a slow triangular sweep of that parameter. */\n  function branches(P, pi, vi, from, to, cols) {\n    const sys = P.sys, p = Float64Array.from(P.sim.base), pts = [];\n    const box = P.fullRanges.map(function (r) { return r.slice(); });\n    let prev = [];\n    for (let c = 0; c <= cols; c++) {\n      p[pi] = from + (to - from) * c / cols;\n      let eq = [];\n      try { eq = DF.findEquilibria(sys, p, box, { seeds: 24, extra: prev }); } catch (e) { eq = []; }\n      prev = eq.map(function (e) { return e.x; });\n      eq.forEach(function (e) { pts.push({ p: p[pi], v: e.x[vi], stable: e.stable }); });\n    }\n    return pts;\n  }\n\n  V.sweep = {\n    label: \"Parameter sweep (hysteresis)\",\n    axes: true,\n    init: function (P) {\n      const v = P.scene.view, sys = P.sys;\n      this.pi = Math.max(0, sys.params.findIndex(function (q) { return q.name === v.param; }));\n      this.vi = Math.max(0, sys.vars.indexOf(v.var || sys.vars[0]));\n      const q = sys.params[this.pi];\n      this.from = v.from === undefined ? q.min : v.from; this.to = v.to === undefined ? q.max : v.to;\n      this.pval = this.from; this.dir = 1; this.trail = [];\n      this.cam = new DF.Camera([0, 1], [[this.from, this.to], P.fullRanges[this.vi]], { pad: 0.04 });\n      this.cam.resize(P.w, P.h, P.inset);\n      this.branchPts = null;\n    },\n    drawStatic: function (P) {\n      this.cam.resize(P.w, P.h, P.inset);\n      const ctx = P.ctx.base, cam = this.cam, th = P.theme;\n      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [P.sys.params[this.pi].name, P.sys.vars[this.vi]], { grid: false });\n      this.frameBox = f;\n      if (P.scene.view.branches === false) return;\n      if (!this.branchPts) this.branchPts = branches(P, this.pi, this.vi, this.from, this.to, Math.min(400, Math.round((f.R - f.L) / 2)));\n      const q = [0, 0], x = [0, 0];\n      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();\n      this.branchPts.forEach(function (b) {\n        x[0] = b.p; x[1] = b.v; cam.project(x, q);\n        ctx.fillStyle = b.stable ? th.ink : th.muted;\n        const r = b.stable ? 1.7 : 1.1;\n        if (b.stable) { ctx.beginPath(); ctx.arc(q[0], q[1], r, 0, 7); ctx.fill(); }\n        else { ctx.globalAlpha = 0.8; ctx.fillRect(q[0] - r, q[1] - r, 2 * r, 2 * r); ctx.globalAlpha = 1; }\n      });\n      ctx.restore();\n    },\n    frame: function (P) {\n      const v = P.scene.view, sim = P.sim;\n      if (!P.userParam) {\n        const speed = (v.speed || 0.0006) * (this.to - this.from);\n        this.pval += this.dir * speed;\n        if (this.pval > this.to) { this.pval = this.to; this.dir = -1; }\n        if (this.pval < this.from) { this.pval = this.from; this.dir = 1; }\n        sim.base[this.pi] = this.pval; sim.updateParams();\n        P.emit(\"param\", { name: P.sys.params[this.pi].name, value: this.pval });\n      } else this.pval = sim.base[this.pi];\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();\n      if (!sim.alive[0]) sim.setMember(0, sim.init);\n      const q = [0, 0];\n      this.cam.project([this.pval, sim.X[this.vi]], q);\n      this.trail.push(q[0], q[1]);\n      if (this.trail.length > 2 * (v.tail || 1400)) this.trail.splice(0, 2);\n      const ctx = P.ctx.top, st = P.scene.style;\n      ctx.clearRect(0, 0, P.w, P.h);\n      ctx.save(); if (this.frameBox) { const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip(); }\n      const m = this.trail.length / 2, G = 16;\n      ctx.strokeStyle = P.palette[1 % P.palette.length]; ctx.lineWidth = st.lineWidth;\n      for (let g = 0; g < G; g++) {\n        const a = Math.floor(g * (m - 1) / G), b = Math.floor((g + 1) * (m - 1) / G);\n        if (b <= a) continue;\n        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * (g + 1) / G);\n        ctx.beginPath(); ctx.moveTo(this.trail[2 * a], this.trail[2 * a + 1]);\n        for (let j = a + 1; j <= b; j++) ctx.lineTo(this.trail[2 * j], this.trail[2 * j + 1]);\n        ctx.stroke();\n      }\n      ctx.globalAlpha = 1; ctx.fillStyle = P.palette[0];\n      ctx.beginPath(); ctx.arc(q[0], q[1], 6, 0, 7); ctx.fill();\n      // Parameter marker on the axis.\n      if (this.frameBox) { ctx.strokeStyle = P.palette[0]; ctx.globalAlpha = 0.5; ctx.beginPath(); ctx.moveTo(q[0], this.frameBox.T); ctx.lineTo(q[0], this.frameBox.B); ctx.stroke(); ctx.globalAlpha = 1; }\n      ctx.restore();\n    },\n    pointer: function (P, kind, px) {\n      if (kind !== \"down\" && kind !== \"drag\") return false;\n      const u = this.cam.unproject(px, 0);\n      P.userParam = true; P.sim.base[this.pi] = Math.min(this.to, Math.max(this.from, u[0])); P.sim.updateParams();\n      P.emit(\"param\", { name: P.sys.params[this.pi].name, value: P.sim.base[this.pi] });\n      return true;\n    },\n    legend: function (P) { return [[\"stable branch\", P.theme.ink, \"dot\"], [\"unstable branch\", P.theme.muted, \"dot\"], [\"state under the sweep\", P.palette[1 % P.palette.length], \"line\"]]; },\n    svg: function (P) {\n      let out = \"\";\n      const q = [0, 0], cam = this.cam, th = P.theme;\n      (this.branchPts || []).forEach(function (b) { cam.project([b.p, b.v], q); out += '<circle cx=\"' + q[0].toFixed(1) + '\" cy=\"' + q[1].toFixed(1) + '\" r=\"' + (b.stable ? 1.7 : 1.1) + '\" fill=\"' + (b.stable ? th.ink : th.muted) + '\"/>'; });\n      out += polyline(this.trail, P.palette[1 % P.palette.length], P.scene.style.lineWidth, P.scene.style.alpha);\n      return out;\n    }\n  };\n\n  // ------------------------------------------------------------ orbit\n  V.orbit = {\n    label: \"Bifurcation diagram\",\n    axes: true,\n    init: function (P) {\n      const v = P.scene.view, sys = P.sys;\n      this.pi = Math.max(0, sys.params.findIndex(function (q) { return q.name === v.param; }));\n      this.vi = Math.max(0, sys.vars.indexOf(v.var || sys.vars[0]));\n      const q = sys.params[this.pi];\n      this.from = v.from === undefined ? q.min : v.from; this.to = v.to === undefined ? q.max : v.to;\n      this.cam = new DF.Camera([0, 1], [[this.from, this.to], P.fullRanges[this.vi]], { pad: 0.02 });\n      this.cam.resize(P.w, P.h, P.inset);\n      this.col = 0; this.state = Float64Array.from(P.sim.init); this.points = [];\n      this.rk = DF.makeRK4(sys.vars.length);\n      this.p = Float64Array.from(P.sim.base);\n    },\n    drawStatic: function (P) {\n      this.cam.resize(P.w, P.h, P.inset);\n      this.frameBox = DF.drawAxes(P.ctx.base, this.cam, P.scene.style.theme, [P.sys.params[this.pi].name, P.sys.vars[this.vi]], { grid: false });\n      this.cols = Math.max(50, Math.round((this.frameBox.R - this.frameBox.L) / (P.scene.view.colWidth || 1.2)));\n      // Existing points are redrawn after a resize.\n      const ctx = P.ctx.trail, q = [0, 0], st = P.scene.style;\n      ctx.clearRect(0, 0, P.w, P.h);\n      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;\n      for (let i = 0; i < this.points.length; i += 2) { this.cam.project([this.points[i], this.points[i + 1]], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2); }\n      ctx.globalAlpha = 1;\n    },\n    frame: function (P) {\n      if (!this.cols || this.col > this.cols) return;\n      const v = P.scene.view, sys = P.sys, dim = sys.vars.length, h = P.sim.h, st = P.scene.style;\n      const transient = v.transient || (sys.time === \"discrete\" ? 300 : Math.round(200 / h));\n      const samples = v.samples || (sys.time === \"discrete\" ? 150 : Math.round(400 / h));\n      const ctx = P.ctx.trail, q = [0, 0], budget = performance.now() + (v.budget || 12);\n      const x = this.state, tmp = P.tmpDx, H = function (i) { return x[i]; };\n      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;\n      while (this.col <= this.cols && performance.now() < budget) {\n        const pv = this.from + (this.to - this.from) * this.col / this.cols;\n        this.p[this.pi] = pv;\n        if (v.follow === false || !x.every(isFinite)) x.set(P.sim.init);\n        let t = 0, prev2 = NaN, prev1 = NaN;\n        const adv = function (self) { if (sys.time === \"discrete\") { sys.f(t, x, self.p, tmp, H); x.set(tmp); t += 1; } else { self.rk(sys.f, t, x, self.p, h, H); t += h; } };\n        for (let s = 0; s < transient; s++) adv(this);\n        for (let s = 0; s < samples; s++) {\n          adv(this);\n          const y = x[this.vi];\n          if (!isFinite(y)) break;\n          if (sys.time === \"discrete\") { this.points.push(pv, y); this.cam.project([pv, y], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2); }\n          else if (prev1 > prev2 && prev1 >= y) {\n            // Local maximum: vertex of the parabola through the last three samples.\n            const den = prev2 - 2 * prev1 + y, peak = den !== 0 ? prev1 - (prev2 - y) * (prev2 - y) / (8 * den) : prev1;\n            this.points.push(pv, peak); this.cam.project([pv, peak], q); ctx.fillRect(q[0] - 0.6, q[1] - 0.6, 1.2, 1.2);\n          }\n          prev2 = prev1; prev1 = y;\n        }\n        this.col++;\n      }\n      ctx.globalAlpha = 1;\n      const top = P.ctx.top; top.clearRect(0, 0, P.w, P.h);\n      if (this.col <= this.cols && this.frameBox) {\n        const px = this.frameBox.L + (this.frameBox.R - this.frameBox.L) * this.col / this.cols;\n        top.strokeStyle = P.palette[1 % P.palette.length]; top.globalAlpha = 0.6; top.beginPath(); top.moveTo(px, this.frameBox.T); top.lineTo(px, this.frameBox.B); top.stroke(); top.globalAlpha = 1;\n      }\n    },\n    legend: function (P) { return [[P.sys.time === \"discrete\" ? \"iterates after a transient\" : \"local maxima after a transient\", P.palette[0], \"dot\"]]; }\n  };\n\n  // ------------------------------------------------------------ density\n  V.density = {\n    label: \"Ensemble density\",\n    axes: true,\n    init: function (P) {\n      const v = P.scene.view;\n      this.mode = v.mode || (P.sys.vars.length === 1 ? \"carpet\" : \"map\");\n      this.vi = Math.max(0, P.sys.vars.indexOf(v.var || P.sys.vars[P.cam.axes[1] === undefined ? 0 : P.cam.axes[1]]));\n      if (P.sys.vars.length === 1) this.vi = 0;\n      this.grid = null; this.img = null;\n      const warm = v.warmup || 0;\n      for (let s = 0; s < warm * P.scene.stepsPerFrame; s++) P.sim.step();\n    },\n    drawStatic: function (P) {\n      const st = P.scene.style;\n      if (this.mode === \"carpet\") {\n        const span = P.scene.view.window || 60;\n        this.cam = new DF.Camera([0, 1], [[-span, 0], P.fullRanges[this.vi]], { pad: 0.0 });\n        this.cam.resize(P.w, P.h, P.inset);\n        this.frameBox = DF.drawAxes(P.ctx.base, this.cam, st.theme, [\"t - t now\", P.sys.vars[this.vi]], {});\n      } else {\n        this.frameBox = DF.drawAxes(P.ctx.base, P.cam, st.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]], {});\n      }\n      const f = this.frameBox;\n      this.nx = Math.max(10, Math.round((f.R - f.L) / (P.scene.view.cell || 3)));\n      this.ny = Math.max(10, Math.round((f.B - f.T) / (P.scene.view.cell || 3)));\n      this.grid = new Float32Array(this.nx * this.ny);\n      if (typeof document !== \"undefined\") {\n        this.off = document.createElement(\"canvas\"); this.off.width = this.nx; this.off.height = this.ny;\n        this.octx = this.off.getContext(\"2d\"); this.img = this.octx.createImageData(this.nx, this.ny);\n      }\n    },\n    frame: function (P) {\n      const sim = P.sim, dim = P.sys.vars.length, v = P.scene.view, st = P.scene.style;\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) sim.step();\n      for (let k = 0; k < sim.n; k++) if (!sim.alive[k]) sim.setMember(k, spawnState(P, new Float64Array(dim)));\n      if (!this.grid || !this.img) return;\n      const nx = this.nx, ny = this.ny, g = this.grid;\n      if (this.mode === \"carpet\") {\n        // shift one column left per `shiftEvery` frames, fill the last column\n        const r = P.fullRanges[this.vi];\n        for (let j = 0; j < ny; j++) { for (let i = 0; i < nx - 1; i++) g[j * nx + i] = g[j * nx + i + 1]; g[j * nx + nx - 1] = 0; }\n        for (let k = 0; k < sim.n; k++) {\n          const y = sim.X[k * dim + this.vi], j = Math.floor((r[1] - y) / (r[1] - r[0]) * ny);\n          if (j >= 0 && j < ny) g[j * nx + nx - 1] += 1;\n        }\n        let mx = 0; for (let j = 0; j < ny; j++) mx = Math.max(mx, g[j * nx + nx - 1]);\n        for (let j = 0; j < ny; j++) g[j * nx + nx - 1] /= (mx || 1);\n      } else {\n        const decay = v.decay === undefined ? 0.85 : v.decay;\n        for (let i = 0; i < g.length; i++) g[i] *= decay;\n        const cam = P.cam, a0 = cam.axes[0], a1 = cam.axes[1], R0 = cam.ranges[0], R1 = cam.ranges[1];\n        for (let k = 0; k < sim.n; k++) {\n          const u = (sim.X[k * dim + a0] - R0[0]) / (R0[1] - R0[0]), w = (R1[1] - sim.X[k * dim + a1]) / (R1[1] - R1[0]);\n          const i = Math.floor(u * nx), j = Math.floor(w * ny);\n          if (i >= 0 && i < nx && j >= 0 && j < ny) g[j * nx + i] += 1;\n        }\n      }\n      let mx = 0; for (let i = 0; i < g.length; i++) mx = Math.max(mx, g[i]);\n      const d = this.img.data, logm = Math.log1p(mx);\n      for (let i = 0; i < g.length; i++) {\n        const u = mx > 0 ? (v.log === false ? g[i] / mx : Math.log1p(g[i]) / logm) : 0;\n        const c = DF.rampRGB(st.ramp, u);\n        d[4 * i] = c[0]; d[4 * i + 1] = c[1]; d[4 * i + 2] = c[2]; d[4 * i + 3] = u < 0.004 ? 0 : Math.round(255 * Math.min(1, 0.15 + 1.2 * u));\n      }\n      this.octx.putImageData(this.img, 0, 0);\n      const ctx = P.ctx.top, f = this.frameBox;\n      ctx.clearRect(0, 0, P.w, P.h);\n      ctx.imageSmoothingEnabled = v.smooth !== false;\n      ctx.drawImage(this.off, f.L, f.T, f.R - f.L, f.B - f.T);\n    }\n  };\n\n  // ------------------------------------------------------------ strobe\n  V.strobe = {\n    label: \"Stroboscopic / Poincare section\",\n    axes: true,\n    init: function (P) {\n      const v = P.scene.view;\n      const per = P.sim.perturbations.find(function (q) { return q.kind === \"periodic\" && q.enabled !== false; });\n      this.period = v.period || (per ? per.period : 2 * Math.PI);\n      this.next = (v.phase || 0) * this.period;\n      while (this.next <= P.sim.t) this.next += this.period;\n      this.prevSide = null; this.prev = null;\n      this.count = 0; this.transient = v.transient === undefined ? 20 : v.transient;\n      this.points = [];\n    },\n    drawStatic: function (P) {\n      this.frameBox = DF.drawAxes(P.ctx.base, P.cam, P.scene.style.theme, [P.sys.vars[P.cam.axes[0]], P.sys.vars[P.cam.axes[1]]], {});\n      const ctx = P.ctx.trail, q = [0, 0], st = P.scene.style, z = Float64Array.from(P.sim.init);\n      ctx.clearRect(0, 0, P.w, P.h); ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha;\n      for (let i = 0; i < this.points.length; i += 2) { z[P.cam.axes[0]] = this.points[i]; z[P.cam.axes[1]] = this.points[i + 1]; P.cam.project(z, q); ctx.fillRect(q[0] - st.pointSize / 2, q[1] - st.pointSize / 2, st.pointSize, st.pointSize); }\n      ctx.globalAlpha = 1;\n    },\n    frame: function (P) {\n      const sim = P.sim, dim = P.sys.vars.length, v = P.scene.view, st = P.scene.style, ctx = P.ctx.trail, q = [0, 0];\n      const a0 = P.cam.axes[0], a1 = P.cam.axes[1];\n      if (st.fade > 0) P.fade(st.fade);\n      ctx.fillStyle = P.palette[0]; ctx.globalAlpha = st.alpha; ctx.globalCompositeOperation = P.theme.blend;\n      const plot = (function (self) { return function (x) { self.points.push(x[a0], x[a1]); if (self.points.length > 400000) self.points.splice(0, 2); P.cam.project(x, q); ctx.fillRect(q[0] - st.pointSize / 2, q[1] - st.pointSize / 2, st.pointSize, st.pointSize); }; })(this);\n      const x = P.tmpX;\n      for (let s = 0; s < P.scene.stepsPerFrame; s++) {\n        if (v.mode === \"section\") {\n          const si = Math.max(0, P.sys.vars.indexOf(v.sectionVar || P.sys.vars[dim - 1])), c = v.sectionValue || 0;\n          if (!this.prev) this.prev = Float64Array.from(sim.X);\n          this.prev.set(sim.X);\n          sim.step();\n          for (let k = 0; k < sim.n; k++) {\n            const y0 = this.prev[k * dim + si] - c, y1 = sim.X[k * dim + si] - c;\n            if (y0 < 0 && y1 >= 0) {\n              const f = y0 / (y0 - y1);\n              for (let i = 0; i < dim; i++) x[i] = this.prev[k * dim + i] + f * (sim.X[k * dim + i] - this.prev[k * dim + i]);\n              if (sim.t > this.transient) plot(x);\n            }\n          }\n        } else {\n          sim.step();\n          if (sim.t + 1e-12 >= this.next) {\n            this.count++;\n            if (this.count > this.transient) for (let k = 0; k < sim.n; k++) { if (sim.alive[k]) plot(sim.member(k, x)); }\n            this.next += this.period;\n          }\n        }\n      }\n      ctx.globalAlpha = 1; ctx.globalCompositeOperation = \"source-over\";\n    },\n    legend: function (P) { return [[P.scene.view.mode === \"section\" ? \"upward crossings of the section\" : \"state at t = t0 + kT, T = \" + DF.fmt(this.period), P.palette[0], \"dot\"]]; }\n  };\n\n  // ------------------------------------------------------------ cobweb\n  V.cobweb = {\n    label: \"Cobweb (1D maps)\",\n    axes: true,\n    init: function (P) {\n      this.x = P.sim.init[0]; this.path = []; this.vi = 0;\n      const r = P.fullRanges[0];\n      this.cam = new DF.Camera([0, 1], [r, r], { pad: 0.03 });\n      this.cam.resize(P.w, P.h, P.inset);\n    },\n    drawStatic: function (P) {\n      this.cam.resize(P.w, P.h, P.inset);\n      const ctx = P.ctx.base, cam = this.cam, th = P.theme, v = P.sys.vars[0];\n      const f = DF.drawAxes(ctx, cam, P.scene.style.theme, [v + \"ₙ\", v + \"ₙ₊₁\"], {});\n      this.frameBox = f;\n      const r = cam.ranges[0], q = [0, 0], out = new Float64Array(P.sys.vars.length), x = Float64Array.from(P.sim.init);\n      ctx.save(); ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();\n      ctx.strokeStyle = th.muted; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);\n      cam.project([r[0], r[0]], q); ctx.beginPath(); ctx.moveTo(q[0], q[1]); cam.project([r[1], r[1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();\n      ctx.setLineDash([]); ctx.strokeStyle = P.palette[2 % P.palette.length]; ctx.lineWidth = 2; ctx.beginPath();\n      for (let i = 0; i <= 600; i++) {\n        x[0] = r[0] + (r[1] - r[0]) * i / 600; P.sys.f(P.sim.t, x, P.sim.p, out);\n        cam.project([x[0], out[0]], q); if (i) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]);\n      }\n      ctx.stroke(); ctx.restore();\n    },\n    frame: function (P) {\n      const every = Math.max(1, Math.round(8 / Math.max(1, P.scene.stepsPerFrame)));\n      if ((P.frameCount % every) !== 0) return;\n      const out = P.tmpDx, x = P.tmpX;\n      x[0] = this.x; P.sys.f(P.sim.t, x, P.sim.p, out);\n      const y = out[0];\n      if (!isFinite(y)) { this.x = P.sim.init[0]; this.path = []; return; }\n      if (!this.path.length) this.path.push(this.x, this.cam.ranges[1][0] < 0 && this.cam.ranges[1][1] > 0 ? 0 : this.cam.ranges[1][0]);\n      this.path.push(this.x, y, y, y);\n      if (this.path.length > 2 * (P.scene.view.tail || 120)) this.path.splice(0, 4);\n      this.x = y; P.sim.t += 1;\n      const ctx = P.ctx.top, st = P.scene.style, q = [0, 0];\n      ctx.clearRect(0, 0, P.w, P.h);\n      ctx.save(); const f = this.frameBox; ctx.beginPath(); ctx.rect(f.L, f.T, f.R - f.L, f.B - f.T); ctx.clip();\n      const m = this.path.length / 2;\n      ctx.strokeStyle = P.palette[0]; ctx.lineWidth = st.lineWidth;\n      for (let j = 1; j < m; j++) {\n        ctx.globalAlpha = st.alpha * (0.1 + 0.9 * j / m);\n        ctx.beginPath(); this.cam.project([this.path[2 * j - 2], this.path[2 * j - 1]], q); ctx.moveTo(q[0], q[1]);\n        this.cam.project([this.path[2 * j], this.path[2 * j + 1]], q); ctx.lineTo(q[0], q[1]); ctx.stroke();\n      }\n      ctx.globalAlpha = 1; ctx.fillStyle = P.palette[0]; ctx.beginPath(); ctx.arc(q[0], q[1], 4, 0, 7); ctx.fill();\n      ctx.restore();\n    },\n    pointer: function (P, kind, px, py) {\n      if (kind !== \"down\") return false;\n      this.x = this.cam.unproject(px, py)[0]; this.path = [];\n      return true;\n    }\n  };\n\n  // Which views suit a system.\n  // Scalar systems are shown against time; state-space views need two variables.\n  function viewsFor(sys) {\n    const multi = sys.vars.length >= 2;\n    const out = multi ? [\"flow\", \"trajectory\", \"timeseries\", \"density\"] : [\"timeseries\", \"density\"];\n    if (sys.time === \"continuous\" && multi) out.push(\"phase\");\n    if (sys.params.length) out.push(\"sweep\", \"orbit\");\n    if (multi) out.push(\"strobe\");\n    if (sys.time === \"discrete\" && !multi) out.push(\"cobweb\");\n    return out;\n  }\n\n  DF.VIEWS = V;\n  DF.viewsFor = viewsFor;\n  DF.contour = contour;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/render/player.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Player: renders one scene into a host element. The same Player runs in\n   the studio, in the <dyn-flow> element and in exported standalone pages.\n\n   A scene is plain JSON and carries its own system text, so it is complete\n   without the catalogue:\n\n     { version, name, system, params: {a: 1}, init: {x: 0.1},\n       dt, stepsPerFrame, seed, n, initMode, spread,\n       perturbations: [{kind, ...}],\n       view:  { type, axes: ['x','y'], ranges: {x: [lo, hi]}, ... },\n       style: { theme, palette, ramp, colorBy, fade, lineWidth, alpha, pointSize },\n       overlay: { title, subtitle, caption, equations, readout, legend } } */\n(function (DF) {\n  \"use strict\";\n\n  const VIEW_DEFAULTS = {\n    flow: { n: 1500, stepsPerFrame: 2, life: [80, 320], spawn: \"box\", fade: 0.06, lineWidth: 1.1, alpha: 0.55, colorBy: \"dominant\" },\n    trajectory: { n: 1, stepsPerFrame: 4, tail: 2500, fade: 0, lineWidth: 1.3, alpha: 0.95, colorBy: \"time\", rotate: 0.25 },\n    timeseries: { n: 1, stepsPerFrame: 2, fade: 0, lineWidth: 1.6, alpha: 0.95, colorBy: \"solid\" },\n    phase: { n: 6, stepsPerFrame: 2, fade: 0, lineWidth: 1.5, alpha: 0.9, colorBy: \"solid\" },\n    sweep: { n: 1, stepsPerFrame: 8, fade: 0, lineWidth: 1.6, alpha: 0.9, colorBy: \"solid\" },\n    orbit: { n: 1, stepsPerFrame: 1, fade: 0, lineWidth: 1, alpha: 0.35, colorBy: \"solid\" },\n    density: { n: 3000, stepsPerFrame: 2, fade: 0, lineWidth: 1, alpha: 1, colorBy: \"solid\" },\n    strobe: { n: 200, stepsPerFrame: 20, fade: 0, lineWidth: 1, alpha: 0.6, colorBy: \"solid\", pointSize: 1.4 },\n    cobweb: { n: 1, stepsPerFrame: 1, fade: 0, lineWidth: 1.3, alpha: 0.9, colorBy: \"solid\" }\n  };\n\n  function clone(o) { return JSON.parse(JSON.stringify(o === undefined ? null : o)); }\n\n  // Fill defaults; never throws on a partial scene.\n  function normalizeScene(input) {\n    const s = clone(input || {}) || {};\n    s.version = 1;\n    s.name = s.name || \"Untitled scene\";\n    s.system = s.system || \"x' = -y\\ny' = x\\ninit x = 1\";\n    s.params = s.params || {};\n    s.init = s.init || {};\n    s.view = s.view || {};\n    s.view.type = s.view.type && DF.VIEWS[s.view.type] ? s.view.type : \"flow\";\n    const d = VIEW_DEFAULTS[s.view.type];\n    s.n = s.n || d.n;\n    s.stepsPerFrame = s.stepsPerFrame || d.stepsPerFrame;\n    s.seed = s.seed === undefined ? 1 : s.seed;\n    s.initMode = s.initMode || \"point\";\n    s.spread = s.spread === undefined ? 0.05 : s.spread;\n    s.perturbations = s.perturbations || [];\n    if (s.view.life === undefined && s.view.type === \"flow\") s.view.life = d.life;\n    if (s.view.spawn === undefined && s.view.type === \"flow\") s.view.spawn = d.spawn;\n    if (s.view.rotate === undefined && d.rotate !== undefined) s.view.rotate = d.rotate;\n    if (s.view.tail === undefined && d.tail !== undefined) s.view.tail = d.tail;\n    s.style = Object.assign({ theme: \"relab-night\", palette: \"relab\", ramp: \"relab-fire\", colorBy: d.colorBy, fade: d.fade, lineWidth: d.lineWidth, alpha: d.alpha, pointSize: d.pointSize || 1.6, renderScale: 1 }, s.style || {});\n    s.overlay = Object.assign({ title: \"\", subtitle: \"\", caption: \"\", equations: false, readout: false, legend: true, position: \"top-left\" }, s.overlay || {});\n    return s;\n  }\n\n  function Player(host, scene, opts) {\n    this.opts = opts || {};\n    this.host = host;\n    this.listeners = {};\n    this.running = false; this.frameCount = 0; this.fps = 0;\n    this.speedMax = 1e-9;\n    this.buildDom();\n    this.load(scene);\n  }\n\n  Player.prototype.on = function (ev, cb) { (this.listeners[ev] = this.listeners[ev] || []).push(cb); return this; };\n  Player.prototype.emit = function (ev, data) { (this.listeners[ev] || []).forEach(function (cb) { try { cb(data); } catch (e) { if (typeof console !== \"undefined\") console.error(e); } }); };\n\n  Player.prototype.buildDom = function () {\n    const h = this.host, doc = h.ownerDocument;\n    if (getComputedStyle(h).position === \"static\") h.style.position = \"relative\";\n    h.style.overflow = \"hidden\";\n    const mk = function (tag, cls, css) { const e = doc.createElement(tag); e.className = cls; e.style.cssText = css; h.appendChild(e); return e; };\n    const fill = \"position:absolute;inset:0;width:100%;height:100%;\";\n    this.bgEl = mk(\"div\", \"df-bg\", fill);\n    this.cv = { base: mk(\"canvas\", \"df-base\", fill), trail: mk(\"canvas\", \"df-trail\", fill), top: mk(\"canvas\", \"df-top\", fill + \"touch-action:none;\") };\n    this.ctx = { base: this.cv.base.getContext(\"2d\"), trail: this.cv.trail.getContext(\"2d\"), top: this.cv.top.getContext(\"2d\") };\n    this.overlayEl = mk(\"div\", \"df-overlay\", \"position:absolute;inset:0;pointer-events:none;font-family:Jost,system-ui,sans-serif;\");\n    const self = this;\n    if (typeof ResizeObserver !== \"undefined\") { this.ro = new ResizeObserver(function () { self.resize(); }); this.ro.observe(h); }\n    if (typeof IntersectionObserver !== \"undefined\") {\n      this.io = new IntersectionObserver(function (e) { self.visible = e[0].isIntersecting; if (self.visible && self.running) self.kick(); });\n      this.io.observe(h);\n    }\n    this.visible = true;\n    this.bindPointer();\n  };\n\n  Player.prototype.bindPointer = function () {\n    const self = this, el = this.cv.top;\n    let down = null;\n    const pos = function (ev) { const r = el.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };\n    el.addEventListener(\"pointerdown\", function (ev) {\n      const p = pos(ev); down = { x: p[0], y: p[1], azim: self.cam ? self.cam.azim : 0, elev: self.cam ? self.cam.elev : 0, moved: false };\n      el.setPointerCapture(ev.pointerId);\n      if (!(self.cam && self.cam.is3D()) && self.view.pointer && self.view.pointer(self, \"down\", p[0], p[1], ev)) self.emit(\"interact\", {});\n    });\n    el.addEventListener(\"pointermove\", function (ev) {\n      if (!down) return;\n      const p = pos(ev), dx = p[0] - down.x, dy = p[1] - down.y;\n      if (Math.abs(dx) + Math.abs(dy) > 3) down.moved = true;\n      if (self.cam && self.cam.is3D()) { self.dragging = true; self.cam.azim = down.azim + dx * 0.01; self.cam.elev = Math.max(-1.5, Math.min(1.5, down.elev + dy * 0.01)); self.clearTrail(); }\n      else if (self.view.pointer) self.view.pointer(self, \"drag\", p[0], p[1], ev);\n    });\n    const up = function () { down = null; self.dragging = false; };\n    el.addEventListener(\"pointerup\", up); el.addEventListener(\"pointercancel\", up);\n    el.addEventListener(\"wheel\", function (ev) {\n      if (!self.cam || !self.scene.view.zoomable) return;\n      ev.preventDefault(); self.cam.zoom = Math.max(0.2, Math.min(20, self.cam.zoom * Math.exp(-ev.deltaY * 0.001))); self.redrawStatic(); self.clearTrail();\n    }, { passive: false });\n  };\n\n  // Load (or reload) a scene: compile, build the simulator, the camera and the view.\n  Player.prototype.load = function (scene) {\n    const s = normalizeScene(scene);\n    let sys;\n    try { sys = DF.compileSystem(s.system); }\n    catch (e) { this.error = e; this.emit(\"error\", e); if (!this.sys) return; sys = this.sys; }\n    this.error = null;\n    this.scene = s; this.sys = sys;\n    this.buildSim();\n    this.emit(\"scene\", s);\n  };\n\n  Player.prototype.buildSim = function () {\n    const s = this.scene, sys = this.sys;\n    const params = sys.params.map(function (q) { return s.params[q.name] !== undefined ? +s.params[q.name] : q.value; });\n    const init = sys.vars.map(function (v, i) { return s.init[v] !== undefined ? +s.init[v] : sys.init[i]; });\n    const dt = s.dt || (sys.time === \"discrete\" ? 1 : 0.01);\n    s.dt = dt;\n    // Full ranges for every variable: scene, then declared, then a probe run.\n    const need = sys.vars.filter(function (v) { return !(s.view.ranges && s.view.ranges[v]) && !sys.ranges[v]; });\n    const probe = need.length ? DF.autoRanges(sys, params, init, { dt: dt, steps: sys.time === \"discrete\" ? 3000 : Math.min(20000, Math.max(3000, Math.round(60 / dt))), perturbations: s.perturbations }) : {};\n    this.fullRanges = sys.vars.map(function (v) { return (s.view.ranges && s.view.ranges[v]) || sys.ranges[v] || probe[v]; });\n    this.sim = new DF.Simulator(sys, {\n      n: s.n, dt: dt, seed: s.seed, params: params, init: init, initMode: s.initMode, spread: s.spread,\n      perturbations: s.perturbations, keepPositive: !!s.keepPositive, box: this.fullRanges\n    });\n    this.tmpX = new Float64Array(sys.vars.length); this.tmpDx = new Float64Array(sys.vars.length);\n    let axes = (s.view.axes || []).map(function (v) { return sys.vars.indexOf(v); }).filter(function (i) { return i >= 0; });\n    if (axes.length < 2) axes = sys.vars.length >= 3 && (s.view.type === \"trajectory\" || s.view.type === \"flow\") && s.view.dim3 !== false ? [0, 1, 2] : sys.vars.length >= 2 ? [0, 1] : [0, 0];\n    if (sys.vars.length === 1 && s.view.type !== \"cobweb\" && s.view.type !== \"density\" && s.view.type !== \"timeseries\" && s.view.type !== \"sweep\" && s.view.type !== \"orbit\") {\n      // A scalar system is shown as x against time through the timeseries view.\n      s.view.type = \"timeseries\";\n    }\n    const self = this;\n    if (s.view.projection === \"simplex\" && sys.vars.length >= 3 && axes.length < 3) axes = [0, 1, 2];\n    this.cam = new DF.Camera(axes, axes.map(function (i) { return self.fullRanges[i]; }), { azim: s.view.azim, elev: s.view.elev, zoom: s.view.zoom, pad: s.view.pad, upAxis: s.view.upAxis, simplex: s.view.projection === \"simplex\" });\n    this.theme = DF.THEMES[s.style.theme] || DF.THEMES[\"relab-night\"];\n    this.palette = (DF.PALETTES[s.style.palette] || DF.PALETTES.relab).colors;\n    this.view = Object.create(DF.VIEWS[s.view.type]);\n    this.userParam = false;\n    this.frameCount = 0;\n    this.speedMax = 1e-9;\n    this.resize(true);\n    // A new scene starts on empty layers; otherwise trails of the previous\n    // scene remain under views that never draw on the trail layer.\n    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);\n    this.view.init(this);\n    this.redrawStatic();\n    this.renderOverlay();\n    if (!this.running) this.drawOnce();\n  };\n\n  Player.prototype.resize = function (skipStatic) {\n    const r = this.host.getBoundingClientRect();\n    const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));\n    const dpr = Math.min(3, (typeof devicePixelRatio !== \"undefined\" ? devicePixelRatio : 1) * (this.scene ? this.scene.style.renderScale || 1 : 1));\n    const changed = w !== this.w || h !== this.h || dpr !== this.dpr;\n    this.w = w; this.h = h; this.dpr = dpr;\n    if (changed) {\n      for (const k in this.cv) {\n        this.cv[k].width = Math.round(w * dpr); this.cv[k].height = Math.round(h * dpr);\n        this.ctx[k].setTransform(dpr, 0, 0, dpr, 0, 0);\n      }\n    }\n    if (!this.scene) return;\n    this.bgEl.style.background = DF.backgroundCSS(this.scene.style.theme);\n    const axesView = this.view && this.view.axes || (this.scene.view.showAxes && this.cam && !this.cam.is3D());\n    const topPad = (this.scene.overlay.title || this.scene.overlay.subtitle) && this.scene.overlay.position !== \"none\" ? 70 : this.scene.overlay.legend ? 34 : 18;\n    this.inset = axesView ? { l: 66, r: 22, t: topPad, b: 52 } : { l: 0, r: 0, t: 0, b: 0 };\n    if (this.cam) this.cam.resize(w, h, this.inset);\n    if (changed && !skipStatic && this.view) { this.clearTrail(); this.redrawStatic(); if (!this.running) this.drawOnce(); }\n  };\n\n  Player.prototype.clearTrail = function () { this.ctx.trail.clearRect(0, 0, this.w, this.h); };\n  Player.prototype.fade = function (a) {\n    if (!(a > 0)) return;\n    const c = this.ctx.trail;\n    c.globalCompositeOperation = \"destination-out\"; c.fillStyle = \"rgba(0,0,0,\" + a + \")\";\n    c.fillRect(0, 0, this.w, this.h); c.globalCompositeOperation = \"source-over\";\n  };\n  Player.prototype.redrawStatic = function () {\n    this.ctx.base.clearRect(0, 0, this.w, this.h);\n    if (this.scene.view.showAxes && this.cam && !this.cam.is3D() && !this.view.axes) DF.drawAxes(this.ctx.base, this.cam, this.scene.style.theme, [this.sys.vars[this.cam.axes[0]], this.sys.vars[this.cam.axes[1]]], {});\n    if (this.scene.view.showAxes && this.cam && this.cam.is3D()) DF.drawBox3D(this.ctx.base, this.cam, this.scene.style.theme);\n    if (this.cam && this.cam.simplex && this.scene.view.showAxes !== false) DF.drawSimplex(this.ctx.base, this.cam, this.scene.style.theme, this.cam.axes.map(function (i) { return this.sys.vars[i]; }, this));\n    if (this.view.drawStatic) this.view.drawStatic(this);\n  };\n\n  // ------------------------------------------------------------ overlay\n  Player.prototype.renderOverlay = function () {\n    const o = this.scene.overlay, el = this.overlayEl, th = this.theme, doc = el.ownerDocument;\n    el.innerHTML = \"\";\n    if (o.position === \"none\") return;\n    const box = doc.createElement(\"div\");\n    box.style.cssText = \"position:absolute;left:18px;top:14px;right:18px;color:\" + th.ink + \";\";\n    if (o.title) { const t = doc.createElement(\"div\"); t.textContent = o.title; t.style.cssText = \"font-weight:600;font-size:clamp(15px,2.4vw,26px);line-height:1.15;\"; box.appendChild(t); }\n    if (o.subtitle) { const t = doc.createElement(\"div\"); t.textContent = o.subtitle; t.style.cssText = \"font-weight:300;font-size:clamp(12px,1.5vw,16px);color:\" + th.muted + \";margin-top:2px;\"; box.appendChild(t); }\n    el.appendChild(box);\n    if (o.equations) {\n      const eq = doc.createElement(\"div\");\n      eq.style.cssText = \"position:absolute;right:18px;top:\" + (this.view.axes ? 70 : 16) + \"px;color:\" + th.ink + \";font-size:clamp(11px,1.4vw,15px);text-align:right;opacity:0.92;\";\n      const lines = DF.systemLatex(this.scene.system).lines;\n      const k = typeof katex !== \"undefined\" ? katex : (typeof window !== \"undefined\" ? window.katex : undefined);\n      if (k) {\n        lines.forEach(function (L) {\n          const d = doc.createElement(\"div\"); d.style.margin = \"2px 0\";\n          try { k.render(L, d, { throwOnError: false, displayMode: false }); } catch (e) { d.textContent = L; }\n          eq.appendChild(d);\n        });\n      } else {\n        // Without KaTeX (a standalone page), the equations are shown as written.\n        String(this.scene.system).split(/\\r?\\n/).map(function (l) { return l.replace(/#.*$/, \"\").trim(); })\n          .filter(function (l) { return /'\\s*=|\\[n\\+1\\]|dt\\s*=|^noise/.test(l); })\n          .forEach(function (L) { const d = doc.createElement(\"div\"); d.textContent = L; d.style.cssText = \"margin:2px 0;font-family:'TeX Gyre Pagella',Palatino,serif;font-style:italic;\"; eq.appendChild(d); });\n      }\n      el.appendChild(eq);\n    }\n    if (o.caption) {\n      const c = doc.createElement(\"div\"); c.textContent = o.caption;\n      c.style.cssText = \"position:absolute;left:18px;bottom:\" + (this.view.axes ? 6 : 12) + \"px;max-width:min(560px,70%);font-size:11px;line-height:1.45;color:\" + th.muted + \";\";\n      el.appendChild(c);\n    }\n    this.readoutEl = null;\n    if (o.readout) {\n      const r = doc.createElement(\"div\");\n      r.style.cssText = \"position:absolute;right:18px;bottom:\" + (this.view.axes ? 6 : 12) + \"px;font:11px ui-monospace,monospace;color:\" + th.muted + \";text-align:right;\";\n      el.appendChild(r); this.readoutEl = r;\n    }\n    if (o.legend && this.view.legend) {\n      const items = this.view.legend(this);\n      if (items && items.length) {\n        const lg = doc.createElement(\"div\");\n        const ins = this.inset || { r: 0, t: 0 };\n        const panelBg = th.dark ? \"rgba(11,6,32,0.55)\" : \"rgba(255,255,255,0.8)\";\n        // With equations in the top-right corner the legend goes to the bottom of the plot.\n        const atBottom = !!o.equations;\n        const bottomPx = this.view.axes ? ins.b + 8 : (o.readout ? 30 : 12);\n        lg.style.cssText = \"position:absolute;right:\" + (ins.r + 10) + \"px;\" + (atBottom ? \"bottom:\" + bottomPx + \"px;\" : \"top:\" + (this.view.axes ? ins.t + 8 : (o.title ? 64 : 12)) + \"px;\") +\n          \"max-width:48%;display:flex;flex-direction:column;align-items:flex-start;gap:1px;padding:5px 9px;border-radius:6px;background:\" + panelBg + \";font-size:11px;color:\" + th.muted + \";\";\n        items.forEach(function (it) {\n          const d = doc.createElement(\"span\"), shape = it[2] || \"line\", c = it[1];\n          const glyph = shape === \"dot\" ? '<span style=\"display:inline-block;width:8px;height:8px;border-radius:50%;background:' + c + ';margin-right:5px;vertical-align:-1px\"></span>'\n            : shape === \"ring\" ? '<span style=\"display:inline-block;width:8px;height:8px;border-radius:50%;border:1.5px solid ' + c + ';margin-right:5px;vertical-align:-1px;box-sizing:border-box\"></span>'\n            : shape === \"dash\" ? '<span style=\"display:inline-block;width:16px;border-top:2px dashed ' + c + ';margin-right:5px;vertical-align:middle\"></span>'\n            : '<span style=\"display:inline-block;width:16px;height:2px;background:' + c + ';margin-right:5px;vertical-align:middle\"></span>';\n          d.innerHTML = glyph; d.appendChild(doc.createTextNode(it[0])); lg.appendChild(d);\n        });\n        el.appendChild(lg);\n      }\n    }\n  };\n  Player.prototype.readoutText = function () {\n    const sim = this.sim, sys = this.sys;\n    let s = (sys.time === \"discrete\" ? \"n = \" + sim.t : \"t = \" + sim.t.toFixed(2));\n    const mod = sim.perturbations.filter(function (q) { return q.enabled !== false && q.param; });\n    mod.forEach(function (q) { const i = sim.paramIndex(q.param); if (i >= 0) s += \"   \" + q.param + \" = \" + DF.fmt(sim.p[i]); });\n    if (this.view.pval !== undefined) s += \"   \" + sys.params[this.view.pi].name + \" = \" + DF.fmt(this.view.pval);\n    return s;\n  };\n\n  // ------------------------------------------------------------ loop\n  Player.prototype.drawOnce = function () {\n    const spf = this.scene.stepsPerFrame;\n    this.scene.stepsPerFrame = 0;\n    try { this.view.frame(this); } catch (e) { /* first frame of some views needs steps */ }\n    this.scene.stepsPerFrame = spf;\n  };\n  Player.prototype.tick = function () {\n    try { this.view.frame(this); }\n    catch (e) { this.pause(); this.error = e; this.emit(\"error\", e); return; }\n    this.frameCount++;\n    if (this.readoutEl && (this.frameCount % 4) === 0) this.readoutEl.textContent = this.readoutText();\n    this.emit(\"frame\", this.frameCount);\n  };\n  Player.prototype.kick = function () {\n    if (this.raf) return;\n    const self = this;\n    let last = 0, acc = 0, frames = 0;\n    const loop = function (ts) {\n      self.raf = 0;\n      if (!self.running || !self.visible) return;\n      if (last) { acc += ts - last; frames++; if (acc > 500) { self.fps = frames * 1000 / acc; acc = 0; frames = 0; } }\n      last = ts;\n      self.tick();\n      self.raf = requestAnimationFrame(loop);\n    };\n    this.raf = requestAnimationFrame(loop);\n  };\n  Player.prototype.play = function () { if (this.running) return; this.running = true; this.emit(\"state\", \"play\"); this.kick(); };\n  Player.prototype.pause = function () { this.running = false; if (this.raf) cancelAnimationFrame(this.raf); this.raf = 0; this.emit(\"state\", \"pause\"); };\n  Player.prototype.toggle = function () { if (this.running) this.pause(); else this.play(); };\n  Player.prototype.step = function () { this.tick(); };\n  Player.prototype.restart = function (seed) {\n    if (seed !== undefined) this.scene.seed = seed;\n    this.clearTrail(); this.ctx.top.clearRect(0, 0, this.w, this.h);\n    this.buildSim();\n    this.emit(\"restart\", this.scene.seed);\n  };\n  // Advance n frames without waiting for the screen (static renders, tests).\n  Player.prototype.advance = function (n) { for (let i = 0; i < n; i++) this.tick(); };\n\n  // Change one parameter live, without restarting.\n  Player.prototype.setParam = function (name, value) {\n    const i = this.sys.params.findIndex(function (q) { return q.name === name; });\n    if (i < 0) return;\n    this.scene.params[name] = value;\n    this.sim.base[i] = value; this.sim.updateParams();\n    if (this.view.pi === i) this.userParam = true;\n    if (this.view.branchPts) this.view.branchPts = null;\n    if (this.scene.view.type === \"phase\" || this.scene.view.type === \"cobweb\") this.redrawStatic();\n    if (this.scene.view.type === \"orbit\") { this.restart(); }\n  };\n  Player.prototype.setStyle = function (patch) {\n    Object.assign(this.scene.style, patch);\n    this.theme = DF.THEMES[this.scene.style.theme] || DF.THEMES[\"relab-night\"];\n    this.palette = (DF.PALETTES[this.scene.style.palette] || DF.PALETTES.relab).colors;\n    this.bgEl.style.background = DF.backgroundCSS(this.scene.style.theme);\n    if (patch.renderScale !== undefined) { this.w = 0; this.resize(); }\n    this.redrawStatic(); this.renderOverlay();\n  };\n  Player.prototype.setOverlay = function (patch) { Object.assign(this.scene.overlay, patch); this.resize(true); this.redrawStatic(); this.renderOverlay(); };\n  Player.prototype.getScene = function () {\n    const s = clone(this.scene), sim = this.sim, sys = this.sys;\n    sys.params.forEach(function (q, i) { s.params[q.name] = sim.base[i]; });\n    if (this.cam && this.cam.is3D()) { s.view.azim = +this.cam.azim.toFixed(4); s.view.elev = +this.cam.elev.toFixed(4); }\n    return s;\n  };\n\n  /* Composite of background, layers and overlay text on one canvas, at the\n     current resolution, for PNG export and video frames. */\n  Player.prototype.composite = function (target) {\n    const W = this.cv.top.width, H = this.cv.top.height, c = target || this.host.ownerDocument.createElement(\"canvas\");\n    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }\n    const g = c.getContext(\"2d\");\n    g.setTransform(1, 0, 0, 1, 0, 0);\n    g.clearRect(0, 0, W, H);\n    g.save(); g.scale(this.dpr, this.dpr); DF.paintBackground(g, this.w, this.h, this.scene.style.theme); g.restore();\n    g.drawImage(this.cv.base, 0, 0); g.drawImage(this.cv.trail, 0, 0); g.drawImage(this.cv.top, 0, 0);\n    g.save(); g.scale(this.dpr, this.dpr); this.paintText(g); g.restore();\n    return c;\n  };\n  Player.prototype.paintText = function (g) {\n    const o = this.scene.overlay, th = this.theme;\n    if (o.position === \"none\") return;\n    let y = 14;\n    g.textBaseline = \"top\"; g.textAlign = \"left\";\n    if (o.title) { g.fillStyle = th.ink; g.font = \"600 \" + Math.round(Math.max(15, Math.min(26, this.w * 0.024))) + \"px Jost, system-ui, sans-serif\"; g.fillText(o.title, 18, y); y += Math.max(18, Math.min(30, this.w * 0.028)); }\n    if (o.subtitle) { g.fillStyle = th.muted; g.font = \"300 \" + Math.round(Math.max(12, Math.min(16, this.w * 0.015))) + \"px Jost, system-ui, sans-serif\"; g.fillText(o.subtitle, 18, y); }\n    const bottom = this.view.axes ? 8 : 14;\n    if (o.caption) {\n      g.fillStyle = th.muted; g.font = \"11px Jost, system-ui, sans-serif\"; g.textBaseline = \"bottom\";\n      const words = o.caption.split(\" \"), maxW = Math.min(560, this.w * 0.7), lines = [];\n      let line = \"\";\n      words.forEach(function (w) { const tt = line ? line + \" \" + w : w; if (g.measureText(tt).width > maxW && line) { lines.push(line); line = w; } else line = tt; });\n      if (line) lines.push(line);\n      lines.forEach(function (L, i) { g.fillText(L, 18, this.h - bottom - (lines.length - 1 - i) * 15); }, this);\n    }\n    if (o.equations) {\n      g.fillStyle = th.ink; g.font = \"italic 13px 'TeX Gyre Pagella', Palatino, serif\"; g.textAlign = \"right\"; g.textBaseline = \"top\";\n      const eqs = String(this.scene.system).split(/\\r?\\n/).map(function (l) { return l.replace(/#.*$/, \"\").trim(); }).filter(function (l) { return /'|\\[n\\+1\\]|dt\\s*=|^noise/.test(l); });\n      eqs.forEach(function (L, i) { g.fillText(L, this.w - 18, (this.view.axes ? 70 : 16) + 18 * i); }, this);\n    }\n    if (o.readout) { g.fillStyle = th.muted; g.font = \"11px ui-monospace, monospace\"; g.textAlign = \"right\"; g.textBaseline = \"bottom\"; g.fillText(this.readoutText(), this.w - 18, this.h - bottom); }\n  };\n\n  /* SVG of the current frame: vector paths for the views that keep them\n     (trajectory, time series, phase plane, sweep), the rendered raster\n     embedded for the others; text as SVG text either way. */\n  Player.prototype.toSVG = function () {\n    const th = this.theme, o = this.scene.overlay, w = this.w, h = this.h;\n    const esc = function (s) { return String(s).replace(/[&<>\"]/g, function (c) { return { \"&\": \"&amp;\", \"<\": \"&lt;\", \">\": \"&gt;\", '\"': \"&quot;\" }[c]; }); };\n    let s = '<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"' + w + '\" height=\"' + h + '\" viewBox=\"0 0 ' + w + \" \" + h + '\">';\n    const bg = th.bg;\n    if (bg[0] === \"solid\") s += '<rect width=\"100%\" height=\"100%\" fill=\"' + bg[1] + '\"/>';\n    else if (bg[0] === \"radial\") s += '<defs><radialGradient id=\"bg\" cx=\"50%\" cy=\"50%\" r=\"72%\"><stop offset=\"0\" stop-color=\"' + bg[1] + '\"/><stop offset=\"0.55\" stop-color=\"' + bg[2] + '\"/><stop offset=\"1\" stop-color=\"' + bg[3] + '\"/></radialGradient></defs><rect width=\"100%\" height=\"100%\" fill=\"url(#bg)\"/>';\n    const base = this.cv.base.toDataURL(\"image/png\");\n    s += '<image width=\"' + w + '\" height=\"' + h + '\" xlink:href=\"' + base + '\"/>';\n    const vec = this.view.svg ? this.view.svg(this) : null;\n    if (vec !== null) s += '<g id=\"vector\">' + vec + \"</g>\";\n    else {\n      const tmp = this.host.ownerDocument.createElement(\"canvas\"); tmp.width = this.cv.top.width; tmp.height = this.cv.top.height;\n      const g = tmp.getContext(\"2d\"); g.drawImage(this.cv.trail, 0, 0); g.drawImage(this.cv.top, 0, 0);\n      s += '<image width=\"' + w + '\" height=\"' + h + '\" xlink:href=\"' + tmp.toDataURL(\"image/png\") + '\"/>';\n    }\n    if (o.position !== \"none\") {\n      if (o.title) s += '<text x=\"18\" y=\"36\" font-family=\"Jost, sans-serif\" font-weight=\"600\" font-size=\"22\" fill=\"' + th.ink + '\">' + esc(o.title) + \"</text>\";\n      if (o.subtitle) s += '<text x=\"18\" y=\"58\" font-family=\"Jost, sans-serif\" font-weight=\"300\" font-size=\"14\" fill=\"' + th.muted + '\">' + esc(o.subtitle) + \"</text>\";\n      if (o.caption) s += '<text x=\"18\" y=\"' + (h - 12) + '\" font-family=\"Jost, sans-serif\" font-size=\"11\" fill=\"' + th.muted + '\">' + esc(o.caption) + \"</text>\";\n    }\n    return s + \"</svg>\";\n  };\n\n  Player.prototype.dispose = function () {\n    this.pause();\n    if (this.ro) this.ro.disconnect();\n    if (this.io) this.io.disconnect();\n    for (const k in this.cv) this.cv[k].remove();\n    this.bgEl.remove(); this.overlayEl.remove();\n  };\n\n  DF.VIEW_DEFAULTS = VIEW_DEFAULTS;\n  DF.normalizeScene = normalizeScene;\n  DF.Player = Player;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/render/export.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Exports: still images (PNG, SVG), moving images (WebM, GIF), and the\n   scene itself (JSON, share link, embed snippet, standalone page). */\n(function (DF) {\n  \"use strict\";\n\n  const E = {};\n\n  E.download = function (blob, name) {\n    const a = document.createElement(\"a\");\n    a.href = URL.createObjectURL(blob); a.download = name;\n    document.body.appendChild(a); a.click(); a.remove();\n    setTimeout(function () { URL.revokeObjectURL(a.href); }, 4000);\n  };\n  E.slug = function (s) { return String(s || \"scene\").toLowerCase().replace(/[^a-z0-9]+/g, \"-\").replace(/^-|-$/g, \"\") || \"scene\"; };\n\n  E.png = function (player) {\n    return new Promise(function (resolve) { player.composite().toBlob(resolve, \"image/png\"); });\n  };\n  E.svg = function (player) { return new Blob([player.toSVG()], { type: \"image/svg+xml\" }); };\n\n  // ------------------------------------------------------------- WebM\n  /* Records the composite of the player for `seconds` at `fps` with the\n     browser's MediaRecorder; resolves to a WebM blob. */\n  E.webm = function (player, seconds, fps, onProgress) {\n    fps = fps || 30;\n    const c = player.composite();\n    const stream = c.captureStream(fps);\n    const type = [\"video/webm;codecs=vp9\", \"video/webm;codecs=vp8\", \"video/webm\"].find(function (t) { return typeof MediaRecorder !== \"undefined\" && MediaRecorder.isTypeSupported(t); });\n    if (!type) return Promise.reject(new Error(\"This browser cannot record WebM\"));\n    const rec = new MediaRecorder(stream, { mimeType: type, videoBitsPerSecond: 12e6 });\n    const chunks = [];\n    rec.ondataavailable = function (e) { if (e.data && e.data.size) chunks.push(e.data); };\n    const wasRunning = player.running;\n    player.pause();\n    return new Promise(function (resolve, reject) {\n      rec.onstop = function () { resolve(new Blob(chunks, { type: \"video/webm\" })); if (wasRunning) player.play(); };\n      rec.onerror = function (e) { reject(e.error || e); };\n      rec.start(250);\n      const total = Math.round(seconds * fps);\n      let k = 0;\n      const frame = function () {\n        if (k >= total) { rec.stop(); return; }\n        player.tick(); player.composite(c); k++;\n        if (onProgress) onProgress(k / total);\n        setTimeout(frame, 1000 / fps);\n      };\n      frame();\n    });\n  };\n\n  // -------------------------------------------------------------- GIF\n  /* GIF89a encoder. The palette holds the 256 most frequent colours of a\n     5-5-5 bit histogram over all frames; pixels map to the nearest palette\n     entry; frames are LZW-compressed with the variable-width code of the\n     GIF specification (W3C, GIF89a, 1990, appendix F). */\n  function buildPalette(frames) {\n    const hist = new Uint32Array(32768);\n    frames.forEach(function (f) {\n      const d = f.data;\n      for (let i = 0; i < d.length; i += 4 * 3) hist[((d[i] >> 3) << 10) | ((d[i + 1] >> 3) << 5) | (d[i + 2] >> 3)]++;\n    });\n    const idx = [];\n    for (let i = 0; i < 32768; i++) if (hist[i]) idx.push(i);\n    idx.sort(function (a, b) { return hist[b] - hist[a]; });\n    const pal = idx.slice(0, 256).map(function (i) { return [((i >> 10) & 31) * 8 + 4, ((i >> 5) & 31) * 8 + 4, (i & 31) * 8 + 4]; });\n    while (pal.length < 256) pal.push([0, 0, 0]);\n    return pal;\n  }\n  function mapper(pal) {\n    const cache = new Int16Array(32768).fill(-1);\n    return function (r, g, b) {\n      const key = ((r >> 3) << 10) | ((g >> 3) << 5) | (b >> 3);\n      let c = cache[key];\n      if (c >= 0) return c;\n      let best = 0, bd = Infinity;\n      for (let i = 0; i < 256; i++) { const p = pal[i], d = (p[0] - r) * (p[0] - r) + (p[1] - g) * (p[1] - g) + (p[2] - b) * (p[2] - b); if (d < bd) { bd = d; best = i; } }\n      cache[key] = best; return best;\n    };\n  }\n  function lzw(indices, minCode) {\n    const out = [];\n    let cur = 0, curBits = 0;\n    const clear = 1 << minCode, eoi = clear + 1;\n    let codeSize = minCode + 1, next = eoi + 1;\n    let dict = new Map();\n    const emit = function (code) {\n      cur |= code << curBits; curBits += codeSize;\n      while (curBits >= 8) { out.push(cur & 255); cur >>>= 8; curBits -= 8; }\n    };\n    emit(clear);\n    let prefix = indices[0];\n    for (let i = 1; i < indices.length; i++) {\n      const k = indices[i], key = prefix * 256 + k, hit = dict.get(key);\n      if (hit !== undefined) { prefix = hit; continue; }\n      emit(prefix);\n      if (next < 4096) {\n        dict.set(key, next++);\n        if (next > (1 << codeSize) && codeSize < 12) codeSize++;\n      } else { emit(clear); dict = new Map(); codeSize = minCode + 1; next = eoi + 1; }\n      prefix = k;\n    }\n    emit(prefix); emit(eoi);\n    if (curBits > 0) out.push(cur & 255);\n    return out;\n  }\n  function encodeGIF(frames, w, h, delayCs) {\n    const pal = buildPalette(frames), map = mapper(pal);\n    const bytes = [];\n    const str = function (s) { for (let i = 0; i < s.length; i++) bytes.push(s.charCodeAt(i)); };\n    const u16 = function (v) { bytes.push(v & 255, (v >> 8) & 255); };\n    str(\"GIF89a\"); u16(w); u16(h); bytes.push(0xf7, 0, 0);\n    pal.forEach(function (c) { bytes.push(c[0], c[1], c[2]); });\n    bytes.push(0x21, 0xff, 11); str(\"NETSCAPE2.0\"); bytes.push(3, 1, 0, 0, 0);\n    frames.forEach(function (f) {\n      bytes.push(0x21, 0xf9, 4, 0x04); u16(delayCs); bytes.push(0, 0);\n      bytes.push(0x2c); u16(0); u16(0); u16(w); u16(h); bytes.push(0);\n      const d = f.data, idx = new Uint8Array(w * h);\n      for (let i = 0, j = 0; j < idx.length; i += 4, j++) idx[j] = map(d[i], d[i + 1], d[i + 2]);\n      bytes.push(8);\n      const data = lzw(idx, 8);\n      for (let i = 0; i < data.length; i += 255) { const n = Math.min(255, data.length - i); bytes.push(n); for (let j = 0; j < n; j++) bytes.push(data[i + j]); }\n      bytes.push(0);\n    });\n    bytes.push(0x3b);\n    return new Uint8Array(bytes);\n  }\n  E.encodeGIF = encodeGIF;\n\n  E.gif = function (player, seconds, fps, maxWidth, onProgress) {\n    fps = fps || 20; maxWidth = maxWidth || 640;\n    const src = player.composite();\n    const scale = Math.min(1, maxWidth / src.width);\n    const w = Math.max(2, Math.round(src.width * scale)), h = Math.max(2, Math.round(src.height * scale));\n    const small = document.createElement(\"canvas\"); small.width = w; small.height = h;\n    const g = small.getContext(\"2d\", { willReadFrequently: true });\n    const frames = [], total = Math.round(seconds * fps);\n    const wasRunning = player.running;\n    player.pause();\n    const per = Math.max(1, Math.round(60 / fps));\n    return new Promise(function (resolve) {\n      let k = 0;\n      const step = function () {\n        if (k >= total) {\n          if (onProgress) onProgress(1, \"encoding\");\n          setTimeout(function () {\n            const bytes = encodeGIF(frames, w, h, Math.round(100 / fps));\n            if (wasRunning) player.play();\n            resolve(new Blob([bytes], { type: \"image/gif\" }));\n          }, 20);\n          return;\n        }\n        for (let i = 0; i < per; i++) player.tick();\n        player.composite(src);\n        g.drawImage(src, 0, 0, w, h);\n        frames.push(g.getImageData(0, 0, w, h));\n        k++;\n        if (onProgress) onProgress(k / total, \"capturing\");\n        setTimeout(step, 0);\n      };\n      step();\n    });\n  };\n\n  // ---------------------------------------------------------- scenes\n  E.sceneJSON = function (scene) { return JSON.stringify(scene, null, 2); };\n\n  function b64url(bytes) {\n    let s = \"\";\n    for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);\n    return btoa(s).replace(/\\+/g, \"-\").replace(/\\//g, \"_\").replace(/=+$/, \"\");\n  }\n  function unb64url(s) {\n    s = s.replace(/-/g, \"+\").replace(/_/g, \"/\"); while (s.length % 4) s += \"=\";\n    const bin = atob(s), out = new Uint8Array(bin.length);\n    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);\n    return out;\n  }\n  // Scene to a compact URL fragment: \"z\" + deflate-raw + base64url, or \"j\" + base64url of JSON.\n  E.encodeScene = function (scene) {\n    const bytes = new TextEncoder().encode(JSON.stringify(scene));\n    if (typeof CompressionStream === \"undefined\") return Promise.resolve(\"j\" + b64url(bytes));\n    const cs = new Blob([bytes]).stream().pipeThrough(new CompressionStream(\"deflate-raw\"));\n    return new Response(cs).arrayBuffer().then(function (buf) { return \"z\" + b64url(new Uint8Array(buf)); });\n  };\n  E.decodeScene = function (code) {\n    const kind = code[0], bytes = unb64url(code.slice(1));\n    if (kind === \"j\") return Promise.resolve(JSON.parse(new TextDecoder().decode(bytes)));\n    const ds = new Blob([bytes]).stream().pipeThrough(new DecompressionStream(\"deflate-raw\"));\n    return new Response(ds).text().then(JSON.parse);\n  };\n\n  function escAttr(s) { return String(s).replace(/&/g, \"&amp;\").replace(/'/g, \"&#39;\").replace(/</g, \"&lt;\"); }\n  E.embedSnippet = function (scene, src) {\n    return '<script src=\"' + (src || \"dynflow.js\") + '\"></script>\\n' +\n      \"<dyn-flow style=\\\"display:block;width:100%;height:420px\\\" controls scene='\" + escAttr(JSON.stringify(scene)) + \"'></dyn-flow>\";\n  };\n  // A page that needs nothing else: the engine source, the scene and a full-window player.\n  E.standaloneHTML = function (scene, opts) {\n    opts = opts || {};\n    if (!DF.SOURCE) throw new Error(\"The standalone export needs the built dynflow.js (run node tools/build.mjs)\");\n    const th = DF.THEMES[scene.style && scene.style.theme] || DF.THEMES[\"relab-night\"];\n    const bg = th.bg[0] === \"solid\" ? th.bg[1] : th.bg[0] === \"radial\" ? th.bg[3] : \"#000\";\n    const title = (scene.overlay && scene.overlay.title) || scene.name || \"DynFlow scene\";\n    return \"<!doctype html>\\n<html lang=\\\"en\\\">\\n<head>\\n<meta charset=\\\"utf-8\\\">\\n<meta name=\\\"viewport\\\" content=\\\"width=device-width, initial-scale=1\\\">\\n<title>\" +\n      title.replace(/</g, \"&lt;\") + \"</title>\\n<link rel=\\\"icon\\\" href=\\\"data:,\\\">\\n<style>html,body{margin:0;height:100%;background:\" + bg + \"}dyn-flow{display:block;width:100vw;height:100vh}</style>\\n</head>\\n<body>\\n\" +\n      \"<dyn-flow\" + (opts.controls === false ? \"\" : \" controls\") + \" scene='\" + escAttr(JSON.stringify(scene)) + \"'></dyn-flow>\\n\" +\n      \"<script>\\n\" + DF.SOURCE.replace(/<\\/script/gi, \"<\\\\/script\") + \"\\n</script>\\n</body>\\n</html>\\n\";\n  };\n\n  DF.Export = E;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/models/catalogue.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* Model catalogue. Every entry is a system written in the formula language\n   of src/core/expr.js plus a default scene; `source` names where the model\n   and its parameter values were taken from (an RElab package, file and\n   line, or the original literature), so that each entry can be traced.\n   The numerical claims in `about` are checked by tests/models.test.mjs. */\n(function (DF) {\n  \"use strict\";\n\n  const M = [];\n  function add(id, name, group, source, about, system, scene) {\n    M.push({ id: id, name: name, group: group, source: source, about: about, system: system.trim().replace(/^ +/gm, \"\"), scene: scene || {} });\n  }\n  const PI2 = 2 * Math.PI;\n\n  // =================================================================== ecology\n  add(\"lotka-volterra\", \"Lotka-Volterra predator and prey\", \"Ecology\",\n    \"janos R/shiny_app.R:404 (lotka_volterra)\",\n    \"Neutral cycles around the coexistence point: H = eaN - m ln N + aP - r ln P is conserved, so every orbit is closed.\",\n    `N' = r*N - a*N*P\n     P' = e*a*N*P - m*P\n     param r = 1 [0.2, 2]\n     param a = 0.05 [0.01, 0.2]\n     param e = 0.4 [0.1, 1]\n     param m = 0.4 [0.05, 1.5]\n     init N = 20\n     init P = 10\n     range N = [0, 80]\n     range P = [0, 60]`,\n    { view: { type: \"phase\", seeds: 8 }, dt: 0.02, stepsPerFrame: 3, overlay: { equations: true } });\n\n  add(\"rosenzweig-macarthur\", \"Rosenzweig-MacArthur\", \"Ecology\",\n    \"janos R/shiny_app.R:417 (rosenzweig)\",\n    \"Paradox of enrichment: the coexistence equilibrium loses stability in a Hopf bifurcation at K = 3.75, and a limit cycle grows with K.\",\n    `N' = r*N*(1 - N/K) - a*N*P/(1 + a*h*N)\n     P' = e*a*N*P/(1 + a*h*N) - m*P\n     param r = 1 [0.2, 2]\n     param K = 6 [1, 12]\n     param a = 1 [0.2, 2]\n     param h = 0.4 [0.1, 1]\n     param e = 0.6 [0.1, 1]\n     param m = 0.3 [0.05, 0.6]\n     init N = 5\n     init P = 1\n     range N = [0, 7]\n     range P = [0, 4.5]`,\n    { view: { type: \"phase\", seeds: 6 }, dt: 0.02, stepsPerFrame: 4, overlay: { equations: true } });\n\n  add(\"hastings-powell\", \"Hastings-Powell food chain\", \"Ecology\",\n    \"kaRma R/demo_system.R:87 and janos R/shiny_app.R:532 (Hastings and Powell 1991)\",\n    \"Three-level food chain with type II responses; chaotic 'teacup' attractor at the classic parameters.\",\n    `X' = X*(1 - X) - a1*X*Y/(1 + b1*X)\n     Y' = a1*X*Y/(1 + b1*X) - a2*Y*Z/(1 + b2*Y) - d1*Y\n     Z' = a2*Y*Z/(1 + b2*Y) - d2*Z\n     param a1 = 5 [3, 6]\n     param b1 = 3 [2, 6.2]\n     param a2 = 0.1 [0.05, 0.2]\n     param b2 = 2 [1, 3]\n     param d1 = 0.4 [0.2, 0.6]\n     param d2 = 0.01 [0.001, 0.06]\n     init X = 0.8\n     init Y = 0.2\n     init Z = 8\n     range X = [0, 1]\n     range Y = [0, 0.5]\n     range Z = [7, 10.5]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 12, overlay: { equations: true } });\n\n  add(\"may-leonard\", \"May-Leonard cyclic competition\", \"Ecology\",\n    \"janos vignettes/chaotic-systems.Rmd:1241 and HiRsch R/systems.R:469 (May and Leonard 1975)\",\n    \"Rock-paper-scissors competition among three species. With alpha > 1 > beta and alpha + beta > 2 orbits approach a heteroclinic cycle through the single-species states and linger ever longer near each.\",\n    `N1' = N1*(1 - N1 - alpha*N2 - beta*N3)\n     N2' = N2*(1 - beta*N1 - N2 - alpha*N3)\n     N3' = N3*(1 - alpha*N1 - beta*N2 - N3)\n     param alpha = 1.5 [0, 2.5]\n     param beta = 0.7 [0, 2.5]\n     init N1 = 0.5\n     init N2 = 0.3\n     init N3 = 0.2\n     range N1 = [0, 1]\n     range N2 = [0, 1]\n     range N3 = [0, 1]`,\n    { view: { type: \"flow\", projection: \"simplex\", spawn: \"mixed\", life: [120, 420] }, n: 1800, dt: 0.05, stepsPerFrame: 3, style: { colorBy: \"dominant\" } });\n\n  add(\"competition-lv\", \"Lotka-Volterra competition\", \"Ecology\",\n    \"janos R/analysis_phase_portrait.R:207\",\n    \"Two competitors with strong interspecific competition (a12 = a21 = 1.5): the coexistence point is a saddle whose stable manifold separates the two exclusion outcomes.\",\n    `N1' = r1*N1*(1 - N1/K1 - a12*N2/K1)\n     N2' = r2*N2*(1 - N2/K2 - a21*N1/K2)\n     param r1 = 1 [0.1, 2]\n     param r2 = 1 [0.1, 2]\n     param K1 = 100 [50, 150]\n     param K2 = 100 [50, 150]\n     param a12 = 1.5 [0, 2]\n     param a21 = 1.5 [0, 2]\n     init N1 = 10\n     init N2 = 80\n     range N1 = [0, 110]\n     range N2 = [0, 110]`,\n    { view: { type: \"phase\", seeds: 10 }, dt: 0.02, stepsPerFrame: 4 });\n\n  add(\"glv4\", \"Four-species generalised Lotka-Volterra\", \"Ecology\",\n    \"janos R/shiny_app.R:491 (glv4)\",\n    \"Competitive community converging to a stable interior equilibrium.\",\n    `x1' = x1*(1.0 - 1.0*x1 - 0.6*x2 - 0.3*x3 - 0.1*x4)\n     x2' = x2*(0.9 - 0.2*x1 - 1.0*x2 - 0.5*x3 - 0.2*x4)\n     x3' = x3*(0.8 - 0.1*x1 - 0.3*x2 - 1.0*x3 - 0.4*x4)\n     x4' = x4*(0.7 - 0.2*x1 - 0.1*x2 - 0.2*x3 - 1.0*x4)\n     init x1 = 0.4\n     init x2 = 0.3\n     init x3 = 0.2\n     init x4 = 0.1\n     range x1 = [0, 1]\n     range x2 = [0, 1]\n     range x3 = [0, 1]\n     range x4 = [0, 1]`,\n    { view: { type: \"timeseries\", window: 60, members: 6 }, n: 6, spread: 0.5, dt: 0.05, stepsPerFrame: 2 });\n\n  add(\"vano-lv4\", \"Chaotic four-species Lotka-Volterra\", \"Ecology\",\n    \"janos vignettes/chaotic-systems.Rmd:1118 (Vano et al. 2006)\",\n    \"Four competitors with a chaotic attractor, the smallest competitive Lotka-Volterra community known to be chaotic.\",\n    `x1' = x1*(1 - x1 - 1.09*x2 - 1.52*x3)\n     x2' = 0.72*x2*(1 - x2 - 0.44*x3 - 1.36*x4)\n     x3' = 1.53*x3*(1 - 2.33*x1 - x3 - 0.47*x4)\n     x4' = 1.27*x4*(1 - 1.21*x1 - 0.51*x2 - 0.35*x3 - x4)\n     init x1 = 0.3013\n     init x2 = 0.4586\n     init x3 = 0.1307\n     init x4 = 0.3557\n     range x1 = [0, 1]\n     range x2 = [0, 1]\n     range x3 = [0, 0.5]\n     range x4 = [0, 1]`,\n    { view: { type: \"trajectory\", axes: [\"x1\", \"x2\", \"x3\"], warmup: 200, tail: 5000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 10 });\n\n  add(\"act-lv\", \"Arneodo-Coullet-Tresser Lotka-Volterra\", \"Ecology\",\n    \"janos vignettes/chaotic-systems.Rmd:1088 (Arneodo, Coullet and Tresser 1980)\",\n    \"Three-species Lotka-Volterra system with a chaotic attractor near mu = 1.5.\",\n    `N1' = N1*(0.5*(1 - N1) + 0.5*(1 - N2) + 0.1*(1 - N3))\n     N2' = N2*(-0.5*(1 - N1) - 0.1*(1 - N2) + 0.1*(1 - N3))\n     N3' = N3*(mu*(1 - N1) + 0.1*(1 - N2) + 0.1*(1 - N3))\n     param mu = 1.52 [1.3, 1.6]\n     init N1 = 0.5\n     init N2 = 0.3\n     init N3 = 0.2`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 5000, rotate: 0.2 }, dt: 0.05, stepsPerFrame: 8 });\n\n  add(\"huisman-weissing\", \"Huisman-Weissing resource competition\", \"Ecology\",\n    \"wadaR R/multispecies_competition.R:232 (Huisman and Weissing 2001, Am. Nat. 157: 488)\",\n    \"Five phytoplankton species compete for three resources with Liebig growth; the surviving community depends sensitively on the initial abundances.\",\n    `aux q1 = max(R1, 0)\n     aux q2 = max(R2, 0)\n     aux q3 = max(R3, 0)\n     aux mu1 = r*min(q1/(0.20 + q1), q2/(0.25 + q2), q3/(0.15 + q3))\n     aux mu2 = r*min(q1/(0.05 + q1), q2/(0.10 + q2), q3/(0.95 + q3))\n     aux mu3 = r*min(q1/(1.00 + q1), q2/(0.05 + q2), q3/(0.35 + q3))\n     aux mu4 = r*min(q1/(0.05 + q1), q2/(1.00 + q2), q3/(0.10 + q3))\n     aux mu5 = r*min(q1/(1.20 + q1), q2/(0.40 + q2), q3/(0.05 + q3))\n     N1' = N1*(mu1 - m)\n     N2' = N2*(mu2 - m)\n     N3' = N3*(mu3 - m)\n     N4' = N4*(mu4 - m)\n     N5' = N5*(mu5 - m)\n     R1' = D*(S - R1) - (0.20*mu1*N1 + 0.10*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.10*mu5*N5)\n     R2' = D*(S - R2) - (0.10*mu1*N1 + 0.20*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.20*mu5*N5)\n     R3' = D*(S - R3) - (0.10*mu1*N1 + 0.10*mu2*N2 + 0.20*mu3*N3 + 0.20*mu4*N4 + 0.10*mu5*N5)\n     param r = 1 [0.5, 1.5]\n     param m = 0.25 [0.1, 0.4]\n     param D = 0.25 [0.1, 0.5]\n     param S = 10 [5, 15]\n     init N1 = 0.1\n     init N2 = 0.1\n     init N3 = 0.1\n     init N4 = 0.1\n     init N5 = 0.1\n     init R1 = 10\n     init R2 = 10\n     init R3 = 10`,\n    { view: { type: \"timeseries\", vars: [\"N1\", \"N2\", \"N3\", \"N4\", \"N5\"], window: 800 }, dt: 0.1, stepsPerFrame: 20, keepPositive: true, style: { palette: \"relab-qualitative\" } });\n\n  add(\"grazing\", \"May grazing model\", \"Ecology\",\n    \"nonautonomeR R/systems.R:1018 (May 1977)\",\n    \"Vegetation under grazing with a type III response: two stable states separated by an unstable one, and hysteresis when the grazing rate c is swept.\",\n    `V' = r*V*(1 - V/K) - c*V^2/(V^2 + V0^2)\n     param r = 1 [0.5, 2]\n     param K = 10 [5, 15]\n     param c = 2 [1, 3]\n     param V0 = 1 [0.5, 2]\n     init V = 7\n     range V = [0, 10]`,\n    { view: { type: \"sweep\", param: \"c\", var: \"V\", from: 1, to: 3 }, dt: 0.05, stepsPerFrame: 8 });\n\n  add(\"allee\", \"Strong Allee effect\", \"Ecology\",\n    \"janos R/analysis_fokker_planck.R:854\",\n    \"Populations below the threshold A decline to extinction; above it they grow to K.\",\n    `x' = r*x*(x/A - 1)*(1 - x/K)\n     param r = 1 [0.2, 2]\n     param A = 0.3 [0.05, 0.6]\n     param K = 1 [0.7, 1.5]\n     init x = 0.5\n     range x = [0, 1.2]`,\n    { view: { type: \"timeseries\", window: 20, members: 12 }, n: 12, initMode: \"box\", dt: 0.02, stepsPerFrame: 2 });\n\n  add(\"nicholson-bailey\", \"Nicholson-Bailey host and parasitoid\", \"Ecology\",\n    \"janos R/shiny_app.R:719 (Nicholson and Bailey 1935)\",\n    \"Host-parasitoid map whose equilibrium is always unstable: oscillations of growing amplitude.\",\n    `H[n+1] = lambda*H*exp(-a*P)\n     P[n+1] = c*H*(1 - exp(-a*P))\n     param lambda = 1.5 [1.05, 2]\n     param a = 0.02 [0.005, 0.05]\n     param c = 1 [0.5, 2]\n     init H = 25\n     init P = 10\n     range H = [0, 200]\n     range P = [0, 200]`,\n    { view: { type: \"trajectory\", tail: 60, dim3: false }, stepsPerFrame: 1, style: { pointSize: 4 } });\n\n  add(\"ricker\", \"Ricker map\", \"Ecology\",\n    \"janos R/shiny_app.R:697 (Ricker 1954)\",\n    \"Density-dependent growth of a single population: period doubling from r = 2 and chaos above r of about 2.69.\",\n    `N[n+1] = N*exp(r*(1 - N/K))\n     param r = 2.7 [1.5, 3.5]\n     param K = 100 [50, 150]\n     init N = 10\n     range N = [0, 400]`,\n    { view: { type: \"orbit\", param: \"r\", var: \"N\", from: 1.5, to: 3.5 }, style: { alpha: 0.3 } });\n\n  add(\"seasonal-rm\", \"Seasonally forced Rosenzweig-MacArthur\", \"Ecology\",\n    \"nonautonomeR vignettes/pullback-attraction.Rmd:655\",\n    \"Prey growth modulated with period T: the limit cycle entrains at weak forcing and becomes chaotic at eps = 0.8; the stroboscopic section shows the attractor.\",\n    `x' = x*(1 + eps*sin(2*pi*t/T))*(1 - x) - a*x*y/(b + x)\n     y' = y*(a*x/(b + x) - d)\n     param a = 1 [0.5, 1.5]\n     param b = 0.3 [0.1, 0.5]\n     param d = 0.35 [0.2, 0.5]\n     param eps = 0.8 [0, 1]\n     param T = 10 [5, 20]\n     init x = 0.5\n     init y = 0.3\n     range x = [0, 1]\n     range y = [0, 0.8]`,\n    { view: { type: \"strobe\", period: 10, transient: 10 }, n: 600, spread: 0.3, dt: 0.02, stepsPerFrame: 50, style: { pointSize: 1.3, alpha: 0.7 } });\n\n  add(\"coleman\", \"Coleman logistic with seasonal carrying capacity\", \"Ecology\",\n    \"nonautonomeR R/systems.R:1119 and vignettes/coleman-model.Rmd:58\",\n    \"Logistic growth with periodic K(t): all positive orbits converge to one periodic orbit, the pullback attractor.\",\n    `x' = r*x*(1 - x/(1 + A*sin(w*t)))\n     param r = 1 [0.2, 3]\n     param A = 0.3 [0, 0.8]\n     param w = 0.2 [0.05, 1]\n     init x = 0.5\n     range x = [0, 2]`,\n    { view: { type: \"timeseries\", window: 60, members: 10 }, n: 10, initMode: \"box\", dt: 0.02, stepsPerFrame: 3 });\n\n  add(\"toggle-switch\", \"Genetic toggle switch\", \"Ecology\",\n    \"janos vignettes/qualitative-analysis.Rmd:641 (Gardner, Cantor and Collins 2000)\",\n    \"Two mutually repressing genes: two stable states separated by a saddle; noise drives switches between them.\",\n    `u' = alpha/(1 + v^beta) - u\n     v' = alpha/(1 + u^gamma) - v\n     param alpha = 3 [1, 6]\n     param beta = 2.5 [1, 4]\n     param gamma = 2.5 [1, 4]\n     init u = 2.5\n     init v = 0.5\n     range u = [0, 3.5]\n     range v = [0, 3.5]`,\n    { view: { type: \"phase\", seeds: 10 }, dt: 0.02, stepsPerFrame: 3 });\n\n  add(\"rock-paper-scissors\", \"Rock-paper-scissors replicator\", \"Ecology\",\n    \"janos R/shiny_app.R:1133 (rps)\",\n    \"Replicator dynamics of the zero-sum rock-paper-scissors game: neutral cycles on the simplex around the mixed equilibrium.\",\n    `aux f1 = -p2 + p3\n     aux f2 = p1 - p3\n     aux f3 = -p1 + p2\n     aux fbar = p1*f1 + p2*f2 + p3*f3\n     p1' = p1*(f1 - fbar)\n     p2' = p2*(f2 - fbar)\n     p3' = p3*(f3 - fbar)\n     init p1 = 0.4\n     init p2 = 0.35\n     init p3 = 0.25\n     range p1 = [0, 1]\n     range p2 = [0, 1]\n     range p3 = [0, 1]`,\n    { view: { type: \"flow\", projection: \"simplex\", life: [200, 500] }, n: 1200, dt: 0.02, stepsPerFrame: 4 });\n\n  // ============================================================== delays\n  add(\"mackey-glass\", \"Mackey-Glass\", \"Delay equations\",\n    \"janos R/shiny_app.R:733 (Mackey and Glass 1977)\",\n    \"Blood-cell production with a delay of 17 time units: chaotic oscillations.\",\n    `x' = a*lag(x, tau)/(1 + lag(x, tau)^n) - b*x\n     param a = 0.2 [0.1, 0.3]\n     param b = 0.1 [0.05, 0.2]\n     param n = 10 [4, 12]\n     param tau = 17 [2, 30]\n     init x = 0.9\n     range x = [0.2, 1.5]`,\n    { view: { type: \"timeseries\", window: 600 }, dt: 0.1, stepsPerFrame: 10 });\n\n  add(\"hutchinson\", \"Hutchinson delayed logistic\", \"Delay equations\",\n    \"janos R/shiny_app.R:743 (Hutchinson 1948)\",\n    \"Logistic growth with delayed feedback: the equilibrium K loses stability when r tau exceeds pi/2.\",\n    `N' = r*N*(1 - lag(N, tau)/K)\n     param r = 1.6 [0.5, 2.5]\n     param K = 100 [50, 150]\n     param tau = 1 [0.2, 2]\n     init N = 20\n     range N = [0, 320]`,\n    { view: { type: \"timeseries\", window: 40 }, dt: 0.02, stepsPerFrame: 4 });\n\n  add(\"nicholson-blowflies\", \"Nicholson blowflies\", \"Delay equations\",\n    \"janos vignettes/introduction.Rmd:265 and symplectoR R/data.R:123 (Gurney, Blythe and Nisbet 1980)\",\n    \"Delayed recruitment with a hump-shaped birth function: large irregular population cycles.\",\n    `N' = P*lag(N, tau)*exp(-lag(N, tau)/N0) - delta*N\n     param P = 8 [2, 12]\n     param N0 = 1 [0.5, 2]\n     param delta = 0.175 [0.1, 0.4]\n     param tau = 15 [5, 20]\n     init N = 3\n     range N = [0, 20]`,\n    { view: { type: \"timeseries\", window: 300 }, dt: 0.1, stepsPerFrame: 6 });\n\n  add(\"delayed-predator-prey\", \"Delayed predator and prey\", \"Delay equations\",\n    \"janos vignettes/qualitative-analysis.Rmd:768\",\n    \"Prey self-regulation acts with a delay tau, which destabilises coexistence into cycles.\",\n    `N' = r*N*(1 - lag(N, tau)/K) - a*N*P\n     P' = b*N*P - d*P\n     param r = 1.5 [0.5, 2.5]\n     param K = 10 [5, 15]\n     param a = 0.2 [0.1, 0.4]\n     param b = 0.1 [0.05, 0.2]\n     param d = 0.5 [0.2, 1]\n     param tau = 3 [0.5, 5]\n     init N = 5\n     init P = 2\n     range N = [0, 14]\n     range P = [0, 12]`,\n    { view: { type: \"trajectory\", tail: 3000 }, dt: 0.02, stepsPerFrame: 6 });\n\n  // ================================================================ chaos\n  add(\"lorenz\", \"Lorenz\", \"Chaotic flows\",\n    \"janos R/shiny_app.R:506 and tuRbulence R/dynamical_systems.R:61 (Lorenz 1963)\",\n    \"Convection model with the butterfly attractor; largest Lyapunov exponent 0.906 at the classical parameters.\",\n    `x' = sigma*(y - x)\n     y' = x*(rho - z) - y\n     z' = x*y - beta*z\n     param sigma = 10 [1, 20]\n     param rho = 28 [0.5, 50]\n     param beta = 2.6666666666666665 [0.5, 4]\n     init x = 1\n     init y = 1\n     init z = 1\n     range x = [-22, 22]\n     range y = [-28, 28]\n     range z = [0, 52]`,\n    { view: { type: \"trajectory\", warmup: 500, tail: 3500, rotate: 0.25 }, dt: 0.005, stepsPerFrame: 5, overlay: { equations: true } });\n\n  add(\"rossler\", \"Rossler\", \"Chaotic flows\",\n    \"janos R/shiny_app.R:520 and tuRbulence R/dynamical_systems.R:254 (Rossler 1976)\",\n    \"Band chaos from a single folded band; largest Lyapunov exponent about 0.07 at a = b = 0.2, c = 5.7.\",\n    `x' = -y - z\n     y' = x + a*y\n     z' = b + z*(x - c)\n     param a = 0.2 [0, 0.4]\n     param b = 0.2 [0, 1]\n     param c = 5.7 [2, 12]\n     init x = 1\n     init y = 1\n     init z = 1\n     range x = [-12, 14]\n     range y = [-14, 11]\n     range z = [0, 24]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });\n\n  add(\"chua\", \"Chua double scroll\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:167 (Chua, Komuro and Matsumoto 1986)\",\n    \"Electronic circuit with a piecewise-linear diode: the double-scroll attractor.\",\n    `x' = alpha*(y - x - (m1*x + 0.5*(m0 - m1)*(abs(x + 1) - abs(x - 1))))\n     y' = x - y + z\n     z' = -beta*y\n     param alpha = 15.6 [8, 20]\n     param beta = 28 [20, 35]\n     param m0 = -1.143 [-1.5, -0.8]\n     param m1 = -0.714 [-1, -0.4]\n     init x = 0.1\n     init y = 0\n     init z = 0\n     range x = [-2.6, 2.6]\n     range y = [-0.5, 0.5]\n     range z = [-4, 4]`,\n    { view: { type: \"trajectory\", warmup: 400, tail: 5000, rotate: 0.2, axes: [\"x\", \"z\", \"y\"] }, dt: 0.005, stepsPerFrame: 6 });\n\n  add(\"chen\", \"Chen\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:204 (Chen and Ueta 1999)\",\n    \"A Lorenz-like system with a double-wing attractor of different topology.\",\n    `x' = a*(y - x)\n     y' = (c - a)*x - x*z + c*y\n     z' = x*y - b*z\n     param a = 35 [30, 40]\n     param b = 3 [1, 5]\n     param c = 28 [20, 30]\n     init x = -10\n     init y = 0\n     init z = 37\n     range x = [-30, 30]\n     range y = [-32, 32]\n     range z = [0, 60]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 4000, rotate: 0.25 }, dt: 0.001, stepsPerFrame: 12 });\n\n  add(\"lu\", \"Lu\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:231 (Lu and Chen 2002)\",\n    \"Intermediate between the Lorenz and Chen attractors.\",\n    `x' = a*(y - x)\n     y' = -x*z + c*y\n     z' = x*y - b*z\n     param a = 36 [30, 40]\n     param b = 3 [1, 5]\n     param c = 20 [12, 28]\n     init x = 0.1\n     init y = 0.2\n     init z = 0.3\n     range x = [-25, 25]\n     range y = [-28, 28]\n     range z = [0, 45]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 4000, rotate: 0.25 }, dt: 0.002, stepsPerFrame: 8 });\n\n  add(\"shimizu-morioka\", \"Shimizu-Morioka\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:269 (Shimizu and Morioka 1980)\",\n    \"Lorenz-type attractor of a laser model.\",\n    `x' = y\n     y' = x - lambda*y - x*z\n     z' = -alpha*z + x^2\n     param alpha = 0.45 [0.3, 0.6]\n     param lambda = 0.75 [0.5, 1]\n     init x = 0.1\n     init y = 0.1\n     init z = 0.1\n     range x = [-2, 2]\n     range y = [-1.6, 1.6]\n     range z = [0, 2.4]`,\n    { view: { type: \"trajectory\", warmup: 400, tail: 5000, rotate: 0.2 }, dt: 0.02, stepsPerFrame: 6 });\n\n  add(\"nose-hoover\", \"Nose-Hoover (Sprott A)\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:301 (Hoover 1985)\",\n    \"Thermostatted oscillator: a conservative flow in which a chaotic sea coexists with invariant tori.\",\n    `x' = y\n     y' = -x - z*y\n     z' = y^2 - a\n     param a = 1 [0.5, 2]\n     init x = 0\n     init y = 5\n     init z = 0\n     range x = [-4, 4]\n     range y = [-5, 5]\n     range z = [-4, 4]`,\n    { view: { type: \"trajectory\", tail: 6000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });\n\n  add(\"sprott-jerk\", \"Sprott minimal jerk\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:343 (Sprott 1997)\",\n    \"Third-order equation x''' = -a x'' + x'^2 - x, one of the simplest chaotic flows.\",\n    `x' = y\n     y' = z\n     z' = -a*z + y^2 - x\n     param a = 2.017 [1.9, 2.1]\n     init x = 0\n     init y = 0\n     init z = 1\n     range x = [-7, 5]\n     range y = [-3, 3]\n     range z = [-4, 3]`,\n    { view: { type: \"trajectory\", warmup: 400, tail: 5000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });\n\n  add(\"thomas\", \"Thomas cyclically symmetric\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:375 (Thomas 1999)\",\n    \"Symmetric flow driven by sines; small damping b gives a labyrinth of chaotic motion.\",\n    `x' = sin(y) - b*x\n     y' = sin(z) - b*y\n     z' = sin(x) - b*z\n     param b = 0.18 [0.05, 0.3]\n     init x = 2.4\n     init y = 2.5\n     init z = 2.6\n     range x = [-5, 5]\n     range y = [-5, 5]\n     range z = [-5, 5]`,\n    { view: { type: \"flow\", life: [200, 600], rotate: 0.15 }, n: 1600, dt: 0.05, stepsPerFrame: 2, style: { colorBy: \"speed\", ramp: \"relab-fire\" } });\n\n  add(\"halvorsen\", \"Halvorsen\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:424\",\n    \"Cyclically symmetric attractor with three lobes.\",\n    `x' = -a*x - 4*y - 4*z - y^2\n     y' = -a*y - 4*z - 4*x - z^2\n     z' = -a*z - 4*x - 4*y - x^2\n     param a = 1.4 [1.2, 1.6]\n     init x = 1\n     init y = 0\n     init z = 0\n     range x = [-12, 8]\n     range y = [-12, 8]\n     range z = [-12, 8]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 4000, rotate: 0.2 }, dt: 0.005, stepsPerFrame: 8 });\n\n  add(\"aizawa\", \"Aizawa\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:457\",\n    \"A sphere-like attractor with a tube along its axis.\",\n    `x' = (z - b)*x - d*y\n     y' = d*x + (z - b)*y\n     z' = c + a*z - z^3/3 - (x^2 + y^2)*(1 + e*z) + f*z*x^3\n     param a = 0.95 [0.7, 1]\n     param b = 0.7 [0.5, 0.9]\n     param c = 0.6 [0.4, 0.8]\n     param d = 3.5 [2, 5]\n     param e = 0.25 [0, 0.5]\n     param f = 0.1 [0, 0.3]\n     init x = 0.1\n     init y = 0\n     init z = 0\n     range x = [-1.6, 1.6]\n     range y = [-1.6, 1.6]\n     range z = [-0.6, 2]`,\n    { view: { type: \"flow\", life: [150, 500], rotate: 0.2 }, n: 1500, dt: 0.01, stepsPerFrame: 3, style: { colorBy: \"speed\", ramp: \"mako\" } });\n\n  add(\"rabinovich-fabrikant\", \"Rabinovich-Fabrikant\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:490 (Rabinovich and Fabrikant 1979)\",\n    \"Modulation instability in a non-equilibrium medium; strongly stretched attractor.\",\n    `x' = y*(z - 1 + x^2) + gamma*x\n     y' = x*(3*z + 1 - x^2) + gamma*y\n     z' = -2*z*(alpha + x*y)\n     param gamma = 0.87 [0.1, 1]\n     param alpha = 1.1 [0.9, 1.3]\n     init x = -1\n     init y = 0\n     init z = 0.5\n     range x = [-2.5, 2.5]\n     range y = [-3, 3]\n     range z = [0, 2]`,\n    { view: { type: \"trajectory\", warmup: 400, tail: 6000, rotate: 0.2 }, dt: 0.002, stepsPerFrame: 10 });\n\n  add(\"lorenz-84\", \"Lorenz-84 atmosphere\", \"Chaotic flows\",\n    \"tuRbulence R/dynamical_systems.R:447 and nonautonomeR R/systems.R:651 (Lorenz 1984)\",\n    \"Low-order model of the westerlies (X) and a travelling wave (Y, Z), chaotic at F = 8, G = 1.\",\n    `X' = -Y^2 - Z^2 - a*X + a*F\n     Y' = X*Y - b*X*Z - Y + G\n     Z' = b*X*Y + X*Z - Z\n     param a = 0.25 [0.1, 0.5]\n     param b = 4 [2, 6]\n     param F = 8 [4, 10]\n     param G = 1 [0, 2]\n     init X = 1\n     init Y = 1\n     init Z = 1\n     range X = [-1.5, 2.8]\n     range Y = [-2.5, 2.8]\n     range Z = [-2.6, 2.6]`,\n    { view: { type: \"trajectory\", warmup: 300, tail: 5000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 6 });\n\n  add(\"lorenz-96\", \"Lorenz-96, five sites\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:653 (Lorenz 1996)\",\n    \"Ring of five sites with advection, damping and forcing F: spatiotemporal chaos.\",\n    `x1' = (x2 - x4)*x5 - x1 + F\n     x2' = (x3 - x5)*x1 - x2 + F\n     x3' = (x4 - x1)*x2 - x3 + F\n     x4' = (x5 - x2)*x3 - x4 + F\n     x5' = (x1 - x3)*x4 - x5 + F\n     param F = 8 [2, 12]\n     init x1 = 8.01\n     init x2 = 8\n     init x3 = 8\n     init x4 = 8\n     init x5 = 8\n     range x1 = [-8, 13]\n     range x2 = [-8, 13]\n     range x3 = [-8, 13]\n     range x4 = [-8, 13]\n     range x5 = [-8, 13]`,\n    { view: { type: \"trajectory\", axes: [\"x1\", \"x2\", \"x3\"], warmup: 400, tail: 4000, rotate: 0.2 }, dt: 0.005, stepsPerFrame: 6 });\n\n  add(\"hyperchaotic-rossler\", \"Hyperchaotic Rossler\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:974 (Rossler 1979)\",\n    \"Four-dimensional flow with two positive Lyapunov exponents.\",\n    `x' = -y - z\n     y' = x + a*y + w\n     z' = b + x*z\n     w' = -c*z + d*w\n     param a = 0.25 [0.2, 0.3]\n     param b = 3 [2, 4]\n     param c = 0.5 [0.3, 0.7]\n     param d = 0.05 [0.02, 0.08]\n     init x = -10\n     init y = -6\n     init z = 0\n     init w = 10\n     range x = [-60, 40]\n     range y = [-40, 40]\n     range z = [0, 100]\n     range w = [0, 120]`,\n    { view: { type: \"trajectory\", axes: [\"x\", \"y\", \"w\"], warmup: 200, tail: 4000, rotate: 0.15 }, dt: 0.005, stepsPerFrame: 10 });\n\n  add(\"newton-leipnik\", \"Newton-Leipnik\", \"Chaotic flows\",\n    \"janos vignettes/chaotic-systems.Rmd:1283 (Leipnik and Newton 1981)\",\n    \"Rigid-body motion with feedback: two coexisting strange attractors, reached from z0 = -0.16 and z0 = -0.18.\",\n    `x' = -a*x + y + 10*y*z\n     y' = -x - a*y + 5*x*z\n     z' = alpha*z - 5*x*y\n     param a = 0.4 [0.3, 0.5]\n     param alpha = 0.175 [0.1, 0.25]\n     init x = 0.349\n     init y = 0\n     init z = -0.16\n     range x = [-0.6, 0.6]\n     range y = [-0.6, 0.6]\n     range z = [-0.6, 0.1]`,\n    { view: { type: \"trajectory\", warmup: 200, tail: 6000, rotate: 0.2 }, dt: 0.01, stepsPerFrame: 8 });\n\n  add(\"charney-devore\", \"Charney-DeVore three modes\", \"Chaotic flows\",\n    \"tuRbulence R/charney_devore.R:84 (Charney and DeVore 1979)\",\n    \"Truncated barotropic flow over topography: zonal and blocked regimes of the midlatitude atmosphere.\",\n    `x' = k*(F - x) - alpha*y*z + beta*y\n     y' = -k*y + alpha*x*z - beta*x - delta*z\n     z' = -k*z + delta*y\n     param F = 1.5 [0.5, 3]\n     param k = 0.1 [0.05, 0.3]\n     param alpha = 1 [0.5, 1.5]\n     param beta = 0.5 [0.2, 1]\n     param delta = 1 [0.5, 1.5]\n     init x = 1\n     init y = 0.1\n     init z = 0.1`,\n    { view: { type: \"trajectory\", warmup: 200, tail: 4000, rotate: 0.2 }, dt: 0.02, stepsPerFrame: 6 });\n\n  // ============================================================== forced\n  add(\"duffing\", \"Forced Duffing oscillator\", \"Forced oscillators\",\n    \"janos R/shiny_app.R:547 (duffing)\",\n    \"Double-well oscillator driven with period 2 pi / 1.2: the stroboscopic map reveals a strange attractor.\",\n    `x' = y\n     y' = -delta*y - alpha*x - beta*x^3 + gamma*cos(omega*t)\n     param delta = 0.3 [0.1, 0.5]\n     param alpha = -1 [-1.5, 1]\n     param beta = 1 [0.5, 1.5]\n     param gamma = 0.5 [0.2, 0.6]\n     param omega = 1.2 [0.8, 1.6]\n     init x = 0.1\n     init y = 0\n     range x = [-2, 2]\n     range y = [-1.6, 1.6]`,\n    { view: { type: \"strobe\", period: PI2 / 1.2, transient: 5 }, n: 500, spread: 1, dt: 0.02, stepsPerFrame: 60, style: { pointSize: 1.2, alpha: 0.6 } });\n\n  add(\"forced-van-der-pol\", \"Forced Van der Pol\", \"Forced oscillators\",\n    \"janos vignettes/chaotic-systems.Rmd:532\",\n    \"Relaxation oscillator driven at a frequency near its own: chaos between locking regimes.\",\n    `x' = y\n     y' = mu*(1 - x^2)*y - x + A*sin(omega*t)\n     param mu = 3 [0.5, 5]\n     param A = 5 [0, 8]\n     param omega = 1.788 [1, 3]\n     init x = 0.1\n     init y = 0\n     range x = [-3, 3]\n     range y = [-8, 8]`,\n    { view: { type: \"trajectory\", tail: 2500 }, dt: 0.005, stepsPerFrame: 8 });\n\n  add(\"forced-pendulum\", \"Forced damped pendulum\", \"Forced oscillators\",\n    \"wadaR R/basins.R:73 and R/wada_detection.R:584\",\n    \"Periodically driven pendulum x'' + gamma x' + sin x = F cos t; at gamma = 0.2, F = 1.66 its basins have the Wada property.\",\n    `x' = v\n     v' = -gamma*v - sin(x) + F*cos(t)\n     param gamma = 0.2 [0.05, 0.5]\n     param F = 1.66 [0.5, 2.5]\n     init x = 0\n     init v = 0\n     range x = [-3.1416, 3.1416]\n     range v = [-4, 4]`,\n    { view: { type: \"strobe\", period: PI2, transient: 10 }, n: 600, spread: 3, dt: 0.02, stepsPerFrame: 60, style: { pointSize: 1.3 } });\n\n  add(\"driven-oscillator\", \"Driven damped linear oscillator\", \"Forced oscillators\",\n    \"janos R/analysis_stroboscopic.R:75\",\n    \"Every orbit converges to the unique periodic response; the stroboscopic points collapse onto one fixed point.\",\n    `x' = v\n     v' = -w0^2*x - 2*zeta*v + F*cos(Om*t)\n     param w0 = 1 [0.5, 2]\n     param zeta = 0.1 [0.02, 0.5]\n     param F = 0.5 [0, 1]\n     param Om = 1.3 [0.5, 2]\n     init x = 1\n     init v = 0\n     range x = [-2, 2]\n     range v = [-2, 2]`,\n    { view: { type: \"strobe\", period: PI2 / 1.3, transient: 0 }, n: 300, spread: 2, dt: 0.01, stepsPerFrame: 30, style: { fade: 0.02, pointSize: 2 } });\n\n  // ========================================================== oscillators\n  add(\"van-der-pol\", \"Van der Pol\", \"Oscillators and excitable media\",\n    \"janos R/shiny_app.R:466 (van der Pol 1926)\",\n    \"Self-sustained oscillation with nonlinear damping; relaxation oscillations for large mu.\",\n    `x' = y\n     y' = mu*(1 - x^2)*y - x\n     param mu = 2 [0.1, 8]\n     init x = 2\n     init y = 0\n     range x = [-3, 3]\n     range y = [-6, 6]`,\n    { view: { type: \"phase\", seeds: 8 }, dt: 0.01, stepsPerFrame: 5 });\n\n  add(\"fitzhugh-nagumo\", \"FitzHugh-Nagumo\", \"Oscillators and excitable media\",\n    \"janos R/shiny_app.R:453 (FitzHugh 1961; Nagumo et al. 1962)\",\n    \"Slow-fast model of a neuron: excitable at low input I, oscillating above a Hopf threshold.\",\n    `v' = v - v^3/3 - w + I\n     w' = eps*(v + a - b*w)\n     param I = 0.5 [0, 2]\n     param eps = 0.08 [0.01, 0.3]\n     param a = 0.7 [0.3, 1]\n     param b = 0.8 [0.3, 1]\n     init v = -1\n     init w = -0.5\n     range v = [-2.5, 2.5]\n     range w = [-1, 2]`,\n    { view: { type: \"phase\", seeds: 8 }, dt: 0.02, stepsPerFrame: 5 });\n\n  add(\"brusselator\", \"Brusselator\", \"Oscillators and excitable media\",\n    \"janos R/shiny_app.R:437 (Prigogine and Lefever 1968)\",\n    \"Autocatalytic reaction scheme: the equilibrium (A, B/A) undergoes a Hopf bifurcation at B = 1 + A^2.\",\n    `X' = A - (B + 1)*X + X^2*Y\n     Y' = B*X - X^2*Y\n     param A = 1 [0.5, 2]\n     param B = 3 [0.3, 6]\n     init X = 1\n     init Y = 1\n     range X = [0, 4.5]\n     range Y = [0, 5.5]`,\n    { view: { type: \"phase\", seeds: 6 }, dt: 0.01, stepsPerFrame: 5 });\n\n  add(\"selkov\", \"Sel'kov glycolysis\", \"Oscillators and excitable media\",\n    \"janos vignettes/introduction.Rmd:209 (Sel'kov 1968)\",\n    \"Glycolytic oscillations: a Hopf bifurcation in the flux b creates a limit cycle.\",\n    `x' = -x + a*y + x^2*y\n     y' = b - a*y - x^2*y\n     param a = 0.1 [0, 0.2]\n     param b = 0.5 [0.1, 1.2]\n     init x = 0.5\n     init y = 0.5\n     range x = [0, 3]\n     range y = [0, 3]`,\n    { view: { type: \"phase\", seeds: 6 }, dt: 0.02, stepsPerFrame: 5 });\n\n  add(\"pendulum\", \"Pendulum\", \"Oscillators and excitable media\",\n    \"classical mechanics\",\n    \"Frictionless pendulum x'' = -sin x: librations inside the separatrix through the saddles at x = +/-pi, rotations outside.\",\n    `x' = v\n     v' = -sin(x) - c*v\n     param c = 0 [0, 0.5]\n     init x = 1\n     init v = 0\n     range x = [-7, 7]\n     range v = [-3.2, 3.2]`,\n    { view: { type: \"phase\", seeds: 14 }, dt: 0.01, stepsPerFrame: 5 });\n\n  add(\"kuramoto-pair\", \"Two coupled phase oscillators\", \"Oscillators and excitable media\",\n    \"nonautonomeR R/systems.R:1453 (kuramoto_pair)\",\n    \"Phase difference psi = phi1 - phi2 obeys psi' = dw - 2K sin psi: phase locking when 2K > |dw|, phase slips otherwise.\",\n    `phi1' = w1 + K*sin(phi2 - phi1)\n     phi2' = w2 + K*sin(phi1 - phi2)\n     param w1 = 6.911503837897544 [5, 8]\n     param w2 = 5.654866776461628 [5, 8]\n     param K = 0.5 [0, 1.5]\n     init phi1 = 0\n     init phi2 = 1\n     range phi1 = [0, 60]\n     range phi2 = [0, 60]`,\n    { view: { type: \"timeseries\", window: 30 }, dt: 0.01, stepsPerFrame: 3 });\n\n  // ================================================================= maps\n  add(\"logistic\", \"Logistic map\", \"Maps\",\n    \"janos R/shiny_app.R:680 and kaRma R/demo_system.R:164 (May 1976)\",\n    \"Period-doubling cascade to chaos; at r = 4 the Lyapunov exponent is ln 2.\",\n    `x[n+1] = r*x*(1 - x)\n     param r = 3.9 [2.5, 4]\n     init x = 0.2\n     range x = [0, 1]`,\n    { view: { type: \"orbit\", param: \"r\", from: 2.5, to: 4 }, style: { alpha: 0.25 } });\n\n  add(\"logistic-cobweb\", \"Logistic map, cobweb\", \"Maps\",\n    \"janos vignettes/qualitative-analysis.Rmd:346\",\n    \"Graphical iteration of x[n+1] = r x (1 - x); click to restart from another x.\",\n    `x[n+1] = r*x*(1 - x)\n     param r = 3.7 [2.5, 4]\n     init x = 0.2\n     range x = [0, 1]`,\n    { view: { type: \"cobweb\", tail: 80 } });\n\n  add(\"henon\", \"Henon map\", \"Maps\",\n    \"janos R/shiny_app.R:688 and nonautonomeR R/systems.R:542 (Henon 1976)\",\n    \"Stretch-and-fold map with a strange attractor; largest Lyapunov exponent 0.419 at a = 1.4, b = 0.3.\",\n    `x[n+1] = 1 - a*x^2 + y\n     y[n+1] = b*x\n     param a = 1.4 [0.8, 1.42]\n     param b = 0.3 [0, 0.4]\n     init x = 0.1\n     init y = 0.1\n     range x = [-1.5, 1.5]\n     range y = [-0.45, 0.45]`,\n    { view: { type: \"flow\", life: [30, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: \"age\", ramp: \"relab-fire\" } });\n\n  add(\"lozi\", \"Lozi map\", \"Maps\",\n    \"janos vignettes/chaotic-systems.Rmd:810 (Lozi 1978)\",\n    \"Piecewise-linear analogue of the Henon map with a strange attractor.\",\n    `x[n+1] = 1 - a*abs(x) + y\n     y[n+1] = b*x\n     param a = 1.7 [1.2, 1.8]\n     param b = 0.5 [0.2, 0.6]\n     init x = 0.1\n     init y = 0\n     range x = [-1.4, 1.4]\n     range y = [-0.7, 0.7]`,\n    { view: { type: \"flow\", life: [30, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: \"age\" } });\n\n  add(\"ikeda\", \"Ikeda map\", \"Maps\",\n    \"janos vignettes/chaotic-systems.Rmd:842 (Ikeda 1979)\",\n    \"Light in a nonlinear optical ring cavity; spiralling strange attractor.\",\n    `aux th = kappa - alpha/(1 + x^2 + y^2)\n     x[n+1] = a + b*(x*cos(th) - y*sin(th))\n     y[n+1] = b*(x*sin(th) + y*cos(th))\n     param a = 1 [0.5, 1.5]\n     param b = 0.9 [0.6, 0.95]\n     param kappa = 0.4 [0, 1]\n     param alpha = 6 [4, 8]\n     init x = 0.1\n     init y = 0.1\n     range x = [-0.6, 2.2]\n     range y = [-2.4, 1]`,\n    { view: { type: \"flow\", life: [40, 200], dim3: false }, n: 3000, stepsPerFrame: 1, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: \"age\", ramp: \"mako\" } });\n\n  add(\"standard-map\", \"Chirikov standard map\", \"Maps\",\n    \"janos R/shiny_app.R:705 (Chirikov 1979)\",\n    \"Area-preserving kicked rotor on the torus: KAM curves break up near K = 0.9716 and a chaotic sea spreads.\",\n    `p[n+1] = mod(p + K*sin(th), 2*pi)\n     th[n+1] = mod(th + p + K*sin(th), 2*pi)\n     param K = 1.2 [0, 3]\n     init p = 0.5\n     init th = 0.5\n     range p = [0, 6.2832]\n     range th = [0, 6.2832]`,\n    { view: { type: \"flow\", life: [400, 2000], axes: [\"th\", \"p\"], dim3: false }, n: 1500, stepsPerFrame: 1, style: { fade: 0, pointSize: 1, alpha: 0.5, colorBy: \"member\", palette: \"relab\" } });\n\n  add(\"zaslavsky\", \"Zaslavsky random map\", \"Maps\",\n    \"nonautonomeR R/systems.R:741 (Namenson, Ott and Antonsen 1996)\",\n    \"Dissipative kicked rotor with a random phase c_n drawn afresh at each step and shared by the ensemble: a snapshot attractor that changes shape every step.\",\n    `aux xn = mod(x + y*(1 - exp(-alpha))/alpha, 2*pi)\n     x[n+1] = xn\n     y[n+1] = kappa*sin(xn + 2*pi*ucommon(0)) + exp(-alpha)*y\n     param alpha = 0.09 [0.02, 0.3]\n     param kappa = 0.5 [0.1, 1]\n     init x = 1\n     init y = 0.1\n     range x = [0, 6.2832]\n     range y = [-4, 4]`,\n    { view: { type: \"flow\", life: \"inf\", dim3: false }, n: 6000, spread: 0, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.6, pointSize: 1.4, alpha: 0.8, colorBy: \"solid\", palette: \"mono-amber\" } });\n\n  add(\"random-baker\", \"Random baker's map\", \"Maps\",\n    \"nonautonomeR R/systems.R:801\",\n    \"Baker's map with a random cut c_n shared by the ensemble: a fractal snapshot attractor with known generalised dimensions.\",\n    `aux c = ucommon(0)\n     x[n+1] = ifelse(y <= c, lambda*x, 0.5 + lambda*x)\n     y[n+1] = ifelse(y <= c, y/c, (y - c)/(1 - c))\n     param lambda = 0.4 [0.1, 0.5]\n     init x = 0.3\n     init y = 0.6\n     range x = [0, 1]\n     range y = [0, 1]`,\n    { view: { type: \"flow\", life: \"inf\", dim3: false }, n: 6000, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.7, pointSize: 1.3, alpha: 0.8, colorBy: \"solid\", palette: \"mono-ice\" } });\n\n  add(\"stark-circle\", \"Quasiperiodically forced circle map\", \"Maps\",\n    \"nonautonomeR R/demo_stark_skew.R:47 (Stark 1999)\",\n    \"Circle diffeomorphism driven by an irrational rotation theta. At a = 0.9 the fibre Lyapunov exponent is negative (about -0.24) and the attractor is the graph of a function of theta; at the package default a = 0.25 it is zero and orbits fill the torus.\",\n    `th[n+1] = mod(th + Om, 2*pi)\n     x[n+1] = mod(x + nu + a*sin(x) + b*sin(th), 2*pi)\n     param Om = 3.883222077450933 [0, 6.2832]\n     param nu = 0.8676521529893011 [0, 6.2832]\n     param a = 0.9 [0, 0.95]\n     param b = 0.35 [0, 1]\n     init th = 0.7\n     init x = 0.2\n     range th = [0, 6.2832]\n     range x = [0, 6.2832]`,\n    { view: { type: \"flow\", life: [300, 1200], dim3: false }, n: 2000, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.02, pointSize: 1.3, alpha: 0.6, colorBy: \"var\", colorVar: \"x\", ramp: \"blackboard\" } });\n\n  add(\"kaplan-yorke\", \"Kaplan-Yorke map\", \"Maps\",\n    \"nonautonomeR R/systems.R:915 (Kaplan and Yorke 1979)\",\n    \"Doubling map driving a contracting coordinate: a strange attractor of Kaplan-Yorke dimension 1 + ln 2 / |ln alpha|. Floating-point doubling discards one bit per step and would reach x = 0 within 53 steps, so a noise of size 1e-9 re-injects the low-order bits.\",\n    `x[n+1] = mod(2*x + 0.000000001*nrand(), 1)\n     y[n+1] = alpha*y + cos(4*pi*x)\n     param alpha = 0.2 [0.05, 0.6]\n     init x = 0.3678\n     init y = 0.6677\n     range x = [0, 1]\n     range y = [-1.3, 1.3]`,\n    { view: { type: \"flow\", life: [20, 60], dim3: false }, n: 3000, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.05, pointSize: 1.2, alpha: 0.5 } });\n\n  add(\"de-jong\", \"Peter de Jong attractor\", \"Maps\",\n    \"Pickover (1990), Computers, Pattern, Chaos and Beauty\",\n    \"Trigonometric map whose attractors are used in generative art; many parameter sets are chaotic.\",\n    `x[n+1] = sin(a*y) - cos(b*x)\n     y[n+1] = sin(c*x) - cos(d*y)\n     param a = 1.4 [-3, 3]\n     param b = -2.3 [-3, 3]\n     param c = 2.4 [-3, 3]\n     param d = -2.1 [-3, 3]\n     init x = 0\n     init y = 0\n     range x = [-2.2, 2.2]\n     range y = [-2.2, 2.2]`,\n    { view: { type: \"flow\", life: \"inf\", dim3: false }, n: 4000, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: \"age\", ramp: \"relab-fire\" } });\n\n  add(\"clifford\", \"Clifford attractor\", \"Maps\",\n    \"Pickover (1990), Computers, Pattern, Chaos and Beauty\",\n    \"Trigonometric map related to the de Jong map, with flowing filamentary attractors.\",\n    `x[n+1] = sin(a*y) + c*cos(a*x)\n     y[n+1] = sin(b*x) + d*cos(b*y)\n     param a = -1.4 [-3, 3]\n     param b = 1.6 [-3, 3]\n     param c = 1 [-3, 3]\n     param d = 0.7 [-3, 3]\n     init x = 0.1\n     init y = 0.1\n     range x = [-2.2, 2.2]\n     range y = [-2.2, 2.2]`,\n    { view: { type: \"flow\", life: \"inf\", dim3: false }, n: 4000, initMode: \"box\", stepsPerFrame: 1, style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: \"var\", colorVar: \"y\", ramp: \"mako\" } });\n\n  // =========================================================== stochastic\n  add(\"ornstein-uhlenbeck\", \"Ornstein-Uhlenbeck process\", \"Stochastic\",\n    \"janos R/shiny_app.R:765 (ou)\",\n    \"Mean-reverting diffusion with stationary variance sigma^2 / (2 theta) = 0.125.\",\n    `x' = -theta*(x - mu)\n     noise x = sigma\n     param theta = 1 [0.1, 3]\n     param mu = 0 [-1, 1]\n     param sigma = 0.5 [0, 1.5]\n     init x = 2\n     range x = [-1.5, 2.2]`,\n    { view: { type: \"density\", window: 20 }, n: 3000, dt: 0.01, stepsPerFrame: 2, style: { ramp: \"relab-fire\" } });\n\n  add(\"double-well\", \"Noisy double well\", \"Stochastic\",\n    \"janos R/shiny_app.R:793 and nonautonomeR vignettes/melancholia-states.Rmd:298\",\n    \"Bistable potential V = x^4/4 - x^2/2 with additive noise: Kramers escapes between the wells at x = -1 and x = 1.\",\n    `x' = x - x^3\n     noise x = sigma\n     param sigma = 0.45 [0, 1]\n     init x = -1\n     range x = [-2, 2]`,\n    { view: { type: \"density\", window: 200 }, n: 3000, dt: 0.02, stepsPerFrame: 3, style: { ramp: \"magma\" } });\n\n  add(\"stochastic-resonance\", \"Stochastic resonance\", \"Stochastic\",\n    \"janos vignettes/noise-in-dynamical-systems.Rmd:372 (Benzi, Sutera and Vulpiani 1981)\",\n    \"A weak periodic tilt, too small to push the state over the barrier alone, becomes visible in the switching when the noise is tuned.\",\n    `x' = x - x^3 + A*cos(w*t)\n     noise x = sigma\n     param A = 0.12 [0, 0.4]\n     param w = 0.1 [0.02, 0.5]\n     param sigma = 0.35 [0, 1]\n     init x = -1\n     range x = [-2, 2]`,\n    { view: { type: \"timeseries\", window: 400, members: 1 }, dt: 0.02, stepsPerFrame: 15 });\n\n  add(\"verhulst-noise\", \"Noise-induced transition (stochastic Verhulst)\", \"Stochastic\",\n    \"janos R/shiny_app.R:784 (Horsthemke and Lefever 1984)\",\n    \"Logistic growth with multiplicative noise: the stationary density changes shape (a P-bifurcation) as sigma crosses sqrt(a).\",\n    `x' = a*x - x^2\n     noise x = sigma*x\n     param a = 1 [0.2, 2]\n     param sigma = 0.8 [0, 2]\n     init x = 1\n     range x = [0, 3]`,\n    { view: { type: \"density\", window: 60 }, n: 3000, dt: 0.005, stepsPerFrame: 6, keepPositive: true, style: { ramp: \"mako\" } });\n\n  add(\"stochastic-lv\", \"Stochastic Lotka-Volterra\", \"Stochastic\",\n    \"janos R/shiny_app.R:802 (sde_lv)\",\n    \"Environmental noise pushes orbits off the conserved cycles of the Lotka-Volterra model.\",\n    `N' = r*N - a*N*P\n     P' = e*a*N*P - m*P\n     noise N = s*N\n     noise P = s*P\n     param r = 1 [0.2, 2]\n     param a = 0.05 [0.01, 0.2]\n     param e = 0.4 [0.1, 1]\n     param m = 0.4 [0.05, 1.5]\n     param s = 0.05 [0, 0.3]\n     init N = 20\n     init P = 10\n     range N = [0, 80]\n     range P = [0, 60]`,\n    { view: { type: \"density\", decay: 0.97 }, n: 3000, dt: 0.01, stepsPerFrame: 4, style: { ramp: \"inferno\" } });\n\n  add(\"coherence-resonance\", \"Coherence resonance\", \"Stochastic\",\n    \"janos vignettes/noise-in-dynamical-systems.Rmd:390 (Pikovsky and Kurths 1997)\",\n    \"Excitable FitzHugh-Nagumo unit below threshold: noise alone triggers spikes, most regular at an intermediate noise level.\",\n    `v' = (v - v^3/3 - w)/eps\n     w' = v + a\n     noise v = sigma\n     param eps = 0.05 [0.01, 0.2]\n     param a = 1.05 [0.9, 1.3]\n     param sigma = 0.5 [0, 2]\n     init v = -1\n     init w = -0.6\n     range v = [-2.5, 2.5]\n     range w = [-1.2, 1.2]`,\n    { view: { type: \"timeseries\", vars: [\"v\"], window: 40 }, dt: 0.002, stepsPerFrame: 20 });\n\n  add(\"maier-stein\", \"Maier-Stein\", \"Stochastic\",\n    \"nonautonomeR R/systems.R:992 (Maier and Stein 1993)\",\n    \"Two attractors at (+/-1, 0) and a saddle at the origin; for alpha different from mu the drift is not a gradient and escape paths bend.\",\n    `x' = x - x^3 - alpha*x*y^2\n     y' = -mu*(1 + x^2)*y\n     noise x = sigma\n     noise y = sigma\n     param alpha = 3 [0, 6]\n     param mu = 1 [0.2, 3]\n     param sigma = 0.25 [0, 0.6]\n     init x = -1\n     init y = 0\n     range x = [-1.8, 1.8]\n     range y = [-1, 1]`,\n    { view: { type: \"density\", decay: 0.96 }, n: 3000, dt: 0.01, stepsPerFrame: 4, style: { ramp: \"magma\" } });\n\n  add(\"geometric-brownian\", \"Geometric Brownian motion\", \"Stochastic\",\n    \"janos R/shiny_app.R:775 (gbm)\",\n    \"Multiplicative noise with drift mu: the mean grows as exp(mu t) while the median grows as exp((mu - sigma^2/2) t).\",\n    `S' = mu*S\n     noise S = sigma*S\n     param mu = 0.08 [-0.2, 0.3]\n     param sigma = 0.3 [0, 0.8]\n     init S = 100\n     range S = [0, 400]`,\n    { view: { type: \"timeseries\", window: 10, members: 20 }, n: 20, dt: 0.002, stepsPerFrame: 10 });\n\n  add(\"noisy-van-der-pol\", \"Noisy Van der Pol\", \"Stochastic\",\n    \"janos vignettes/qualitative-analysis.Rmd:689\",\n    \"Limit cycle blurred by velocity noise: the ensemble diffuses along the cycle and loses its phase.\",\n    `x' = y\n     y' = mu*(1 - x^2)*y - x\n     noise y = sigma\n     param mu = 1.5 [0.2, 4]\n     param sigma = 0.5 [0, 1.5]\n     init x = 2\n     init y = 0\n     range x = [-3, 3]\n     range y = [-5, 5]`,\n    { view: { type: \"density\", decay: 0.9 }, n: 3000, spread: 0.01, dt: 0.01, stepsPerFrame: 3, style: { ramp: \"relab-fire\" } });\n\n  // ================================================= tipping and nonautonomous\n  add(\"r-tipping\", \"Rate-induced tipping\", \"Tipping and nonautonomous\",\n    \"normal form of Ashwin, Wieczorek, Vitolo and Cox (2012), Phil. Trans. R. Soc. A 370: 1166\",\n    \"In the frame y = x + lambda of x' = (x + lambda)^2 - 1, a ramp of lambda at rate r adds r to the drift. For r < 1 the state shifts to y = -sqrt(1 - r) and recovers when the ramp ends; for r > 1 no equilibrium exists and the state escapes (R-tipping) if the ramp outlasts the escape time, although every frozen lambda is safe.\",\n    `y' = y^2 - 1 + r*step(t - t0)*step(t0 + L - t)\n     param r = 1.2 [0, 3]\n     param t0 = 2 [0, 10]\n     param L = 12 [1, 30]\n     init y = -1\n     range y = [-2, 3]`,\n    { view: { type: \"timeseries\", window: 25, members: 1 }, dt: 0.005, stepsPerFrame: 4, overlay: { equations: true } });\n\n  add(\"fold-normal-form\", \"Saddle-node (fold) normal form\", \"Tipping and nonautonomous\",\n    \"janos vignettes/advanced-dynamics.Rmd:213\",\n    \"x' = mu + x^2: two equilibria for mu < 0 that collide and vanish at mu = 0.\",\n    `x' = mu + x^2\n     param mu = -1 [-1, 0.5]\n     init x = -1\n     range x = [-1.6, 1.6]`,\n    { view: { type: \"sweep\", param: \"mu\", from: -1, to: 0.3, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 8 });\n\n  add(\"cusp\", \"Cusp catastrophe\", \"Tipping and nonautonomous\",\n    \"normal form of the cusp catastrophe (Thom 1972; Zeeman 1977)\",\n    \"x' = r + a x - x^3: for a > 0 two stable branches coexist between the folds at r = +/-2 (a/3)^(3/2), so a slow sweep of r jumps and shows hysteresis.\",\n    `x' = r + a*x - x^3\n     param r = 0 [-0.8, 0.8]\n     param a = 1 [-0.5, 2]\n     init x = -1\n     range x = [-1.6, 1.6]`,\n    { view: { type: \"sweep\", param: \"r\", from: -0.8, to: 0.8, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 10, overlay: { equations: true } });\n\n  add(\"pitchfork\", \"Pitchfork normal form\", \"Tipping and nonautonomous\",\n    \"janos R/analysis_bifurcation_sweep.R:80\",\n    \"x' = a x - x^3: the origin splits into two symmetric stable states at a = 0.\",\n    `x' = a*x - x^3\n     param a = 1 [-1, 1.5]\n     init x = 0.05\n     range x = [-1.4, 1.4]`,\n    { view: { type: \"sweep\", param: \"a\", from: -1, to: 1.5, speed: 0.0005 }, dt: 0.01, stepsPerFrame: 8 });\n\n  add(\"hopf-normal-form\", \"Hopf normal form\", \"Tipping and nonautonomous\",\n    \"janos vignettes/advanced-dynamics.Rmd:241\",\n    \"Supercritical Hopf bifurcation: for mu > 0 a stable limit cycle of radius sqrt(mu) surrounds the origin.\",\n    `x' = mu*x - y - x*(x^2 + y^2)\n     y' = x + mu*y - y*(x^2 + y^2)\n     param mu = 0.5 [-0.5, 1]\n     init x = 0.1\n     init y = 0.1\n     range x = [-1.3, 1.3]\n     range y = [-1.3, 1.3]`,\n    { view: { type: \"phase\", seeds: 8 }, dt: 0.02, stepsPerFrame: 3 });\n\n  add(\"bogdanov-takens\", \"Bogdanov-Takens normal form\", \"Tipping and nonautonomous\",\n    \"normal form (Kuznetsov 2004, Elements of Applied Bifurcation Theory, ch. 8)\",\n    \"Codimension-two unfolding with fold, Hopf and homoclinic bifurcation curves meeting at the origin of (b1, b2).\",\n    `x' = y\n     y' = b1 + b2*x + x^2 + x*y\n     param b1 = -0.1 [-0.5, 0.2]\n     param b2 = 0.2 [-0.5, 0.5]\n     init x = 0\n     init y = 0.1\n     range x = [-1, 1]\n     range y = [-1, 1]`,\n    { view: { type: \"phase\", seeds: 10 }, dt: 0.01, stepsPerFrame: 4 });\n\n  add(\"stommel\", \"Stommel two-box ocean\", \"Tipping and nonautonomous\",\n    \"tuRbulence R/stommel.R:73 (Stommel 1961)\",\n    \"Thermohaline circulation with temperature T and salinity S: bistability between strong and weak overturning, with hysteresis in the freshwater forcing eta2.\",\n    `T' = eta1 - T*(1 + abs(T - S))\n     S' = eta2 - S*(eta3 + abs(T - S))\n     param eta1 = 3 [2, 4]\n     param eta2 = 1 [0.5, 1.5]\n     param eta3 = 0.3 [0.1, 0.6]\n     init T = 2\n     init S = 1\n     range T = [0, 3.5]\n     range S = [0, 3.5]`,\n    { view: { type: \"sweep\", param: \"eta2\", var: \"S\", from: 0.5, to: 1.5, speed: 0.0004 }, dt: 0.02, stepsPerFrame: 10 });\n\n  add(\"lorenz84-forced\", \"Lorenz-84 under climate change\", \"Tipping and nonautonomous\",\n    \"nonautonomeR R/systems.R:651 and vignettes/ergodicity.Rmd:174 (Jánosi, Tél and co-authors)\",\n    \"Annual forcing F(t) = F0 + 2 sin(2 pi t / 73) with F0 ramping down: a cloud of 1500 members, all driven alike, traces a snapshot attractor that deforms as the climate changes.\",\n    `X' = -Y^2 - Z^2 - a*X + a*(F0 + AF*sin(2*pi*t/73))\n     Y' = X*Y - b*X*Z - Y + G\n     Z' = b*X*Y + X*Z - Z\n     param a = 0.25 [0.1, 0.5]\n     param b = 4 [2, 6]\n     param G = 1 [0, 2]\n     param F0 = 9.5 [6, 10]\n     param AF = 2 [0, 3]\n     init X = 1\n     init Y = 0\n     init Z = 0\n     range X = [-1.5, 3.2]\n     range Y = [-3, 3]\n     range Z = [-3, 3]`,\n    { view: { type: \"flow\", life: \"inf\", rotate: 0.1 }, n: 1500, spread: 1.2, initMode: \"ball\", dt: 0.02, stepsPerFrame: 3, perturbations: [{ kind: \"ramp\", param: \"F0\", rate: -0.000274, t0: 0, span: -2 }], style: { fade: 0.25, colorBy: \"speed\", ramp: \"relab-fire\" }, overlay: { readout: true } });\n\n  add(\"lorenz-drift\", \"Lorenz with drifting rho\", \"Tipping and nonautonomous\",\n    \"nonautonomeR vignettes/pullback-scenes.Rmd:145\",\n    \"Rayleigh number rho(t) = 28 + 0.02 t: the attractor inflates slowly while the ensemble spreads over it.\",\n    `x' = sigma*(y - x)\n     y' = x*(rho - z) - y\n     z' = x*y - beta*z\n     param sigma = 10 [5, 15]\n     param rho = 28 [20, 40]\n     param beta = 2.6666666666666665 [2, 3]\n     init x = 1\n     init y = 1\n     init z = 20\n     range x = [-26, 26]\n     range y = [-34, 34]\n     range z = [0, 62]`,\n    { view: { type: \"flow\", life: \"inf\", rotate: 0.12 }, n: 1500, initMode: \"ball\", spread: 6, dt: 0.005, stepsPerFrame: 4, perturbations: [{ kind: \"ramp\", param: \"rho\", rate: 0.02, t0: 0, span: 12 }], style: { fade: 0.2 }, overlay: { readout: true } });\n\n  add(\"duffing-drift\", \"Duffing with drifting forcing\", \"Tipping and nonautonomous\",\n    \"nonautonomeR R/systems.R:963 (Janosi and Tel 2024)\",\n    \"Forcing amplitude eps(t) = 0.4 + 0.00045 t: the stroboscopic cloud of an ensemble is a snapshot attractor whose shape follows the drift.\",\n    `x' = v\n     v' = x - x^3 - 2*beta*v + eps*cos(omega*t)\n     param beta = 0.2 [0.1, 0.4]\n     param eps = 0.4 [0.2, 0.6]\n     param omega = 1 [0.8, 1.2]\n     init x = 0.5\n     init v = 0\n     range x = [-2, 2]\n     range v = [-1.6, 1.6]`,\n    { view: { type: \"strobe\", period: PI2, transient: 3 }, n: 1500, spread: 1.5, initMode: \"ball\", dt: 2 * Math.PI / 128, stepsPerFrame: 128, perturbations: [{ kind: \"ramp\", param: \"eps\", rate: 0.00045, t0: 0, span: 0.3 }], style: { fade: 0.35, pointSize: 1.6, alpha: 0.8 }, overlay: { readout: true } });\n\n  add(\"tilted-well\", \"Double well with a slow tilt\", \"Tipping and nonautonomous\",\n    \"nonautonomeR vignettes/melancholia-states.Rmd:524\",\n    \"The tilt lambda is ramped slowly: the occupied well loses stability at a fold, and the ensemble switches, earlier when noise is present.\",\n    `x' = x - x^3 + lambda\n     noise x = sigma\n     param lambda = -0.6 [-0.6, 0.6]\n     param sigma = 0.15 [0, 0.5]\n     init x = -1.2\n     range x = [-1.8, 1.8]`,\n    { view: { type: \"density\", window: 120 }, n: 3000, dt: 0.02, stepsPerFrame: 3, perturbations: [{ kind: \"ramp\", param: \"lambda\", rate: 0.01, t0: 0, span: 1.2 }], style: { ramp: \"magma\" }, overlay: { readout: true } });\n\n  // ============================================================ epidemics\n  add(\"sir\", \"SIR epidemic\", \"Epidemics\",\n    \"janos R/shiny_app.R:477 (Kermack and McKendrick 1927)\",\n    \"Closed epidemic with basic reproduction number R0 = beta / gamma = 3.\",\n    `S' = -beta*S*I/N0\n     I' = beta*S*I/N0 - gamma*I\n     R' = gamma*I\n     param beta = 0.3 [0.05, 1]\n     param gamma = 0.1 [0.02, 0.5]\n     param N0 = 1000 [100, 2000]\n     init S = 999\n     init I = 1\n     init R = 0\n     range S = [0, 1000]\n     range I = [0, 1000]\n     range R = [0, 1000]`,\n    { view: { type: \"timeseries\", window: 160 }, dt: 0.05, stepsPerFrame: 4 });\n\n  add(\"seasonal-sir\", \"Seasonally forced SIR\", \"Epidemics\",\n    \"janos vignettes/chaotic-systems.Rmd:1204 (Olsen and Schaffer 1990)\",\n    \"Measles-like epidemic with seasonal transmission: irregular outbreaks from a smooth forcing.\",\n    `S' = mu - beta0*(1 + beta1*cos(2*pi*t))*S*I - mu*S\n     I' = beta0*(1 + beta1*cos(2*pi*t))*S*I - (gam + mu)*I + eps\n     param mu = 0.02 [0.01, 0.04]\n     param beta0 = 1800 [1000, 2500]\n     param beta1 = 0.08 [0, 0.3]\n     param gam = 100 [50, 150]\n     param eps = 0.000001 [0, 0.00001]\n     init S = 0.065\n     init I = 0.0002\n     range S = [0.045, 0.085]\n     range I = [0, 0.002]`,\n    { view: { type: \"timeseries\", window: 20, vars: [\"I\"] }, dt: 0.0005, stepsPerFrame: 40 });\n\n  DF.CATALOGUE = M;\n  DF.catalogueGroups = function () {\n    const seen = [];\n    M.forEach(function (m) { if (seen.indexOf(m.group) < 0) seen.push(m.group); });\n    return seen;\n  };\n  // Scene for a catalogue entry: the entry's defaults with its system text.\n  DF.sceneFor = function (id, overrides) {\n    const m = M.find(function (e) { return e.id === id; });\n    if (!m) throw new Error(\"No model '\" + id + \"' in the catalogue\");\n    const s = JSON.parse(JSON.stringify(m.scene));\n    s.name = m.name; s.model = m.id; s.system = m.system;\n    s.overlay = Object.assign({ title: m.name }, s.overlay || {});\n    return Object.assign(s, overrides || {});\n  };\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\n// ---- src/component/dyn-flow.js\n// SPDX-License-Identifier: GPL-3.0-or-later\n// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab\n/* <dyn-flow>: a scene as an HTML element. One script tag and one element\n   place a live figure in a web page, a Quarto or R Markdown document, a\n   reveal.js slide or a pkgdown article.\n\n     <script src=\"dynflow.js\"></script>\n     <dyn-flow scene='{\"system\": \"...\", \"view\": {\"type\": \"trajectory\"}}'></dyn-flow>\n     <dyn-flow model=\"lorenz\" theme=\"blackboard\" controls></dyn-flow>\n     <dyn-flow src=\"figures/fold.json\" paused></dyn-flow>\n\n   Attributes: scene (JSON), src (URL of a scene file), model (catalogue\n   identifier), theme, view, title, controls (show play, restart and\n   fullscreen buttons), paused (do not start automatically), static (draw a\n   still after a number of frames given by the attribute, default 240). The\n   element exposes .player, .play(), .pause(), .restart(), .setParam(name, v)\n   and .scene. A figure starts only when it is visible and shows a still frame\n   when the reader asks for reduced motion. */\n(function (DF) {\n  \"use strict\";\n  if (typeof HTMLElement === \"undefined\" || typeof customElements === \"undefined\") return;\n\n  const CSS = \":host{display:block;position:relative;min-height:200px;height:360px;contain:content}\" +\n    \".stage{position:absolute;inset:0}\" +\n    \".ctl{position:absolute;right:10px;bottom:10px;display:flex;gap:6px;opacity:0;transition:opacity .2s;z-index:5}\" +\n    \":host(:hover) .ctl,.ctl:focus-within{opacity:1}\" +\n    \"button{all:unset;cursor:pointer;width:30px;height:30px;border-radius:50%;display:grid;place-items:center;background:rgba(23,12,58,.72);color:#fff;font:14px system-ui}\" +\n    \"button:hover{background:#EE6A24}button:focus-visible{outline:2px solid #FB9E07}\" +\n    \".err{position:absolute;inset:0;display:grid;place-items:center;padding:16px;font:13px system-ui;color:#CF4446;background:#fff8f5;text-align:center}\";\n\n  class DynFlowElement extends HTMLElement {\n    static get observedAttributes() { return [\"scene\", \"model\", \"theme\", \"view\"]; }\n    connectedCallback() {\n      if (this._root) return;\n      this._root = this.attachShadow({ mode: \"open\" });\n      const style = document.createElement(\"style\"); style.textContent = CSS;\n      this._stage = document.createElement(\"div\"); this._stage.className = \"stage\";\n      this._root.append(style, this._stage);\n      if (this.hasAttribute(\"controls\")) this.buildControls();\n      this.loadScene();\n    }\n    disconnectedCallback() { if (this.player) { this.player.dispose(); this.player = null; } this._root = null; }\n    attributeChangedCallback(name, oldV, newV) { if (this._root && oldV !== newV && this.player) this.loadScene(); }\n\n    buildControls() {\n      const c = document.createElement(\"div\"); c.className = \"ctl\";\n      const mk = (label, title, fn) => { const b = document.createElement(\"button\"); b.textContent = label; b.title = title; b.setAttribute(\"aria-label\", title); b.addEventListener(\"click\", fn); c.appendChild(b); return b; };\n      this._playBtn = mk(\"❚❚\", \"Pause\", () => this.player && this.player.toggle());\n      mk(\"↻\", \"Restart\", () => this.restart());\n      mk(\"⛶\", \"Full screen\", () => { if (document.fullscreenElement) document.exitFullscreen(); else this.requestFullscreen && this.requestFullscreen(); });\n      this._root.appendChild(c);\n    }\n\n    async loadScene() {\n      let scene;\n      try {\n        if (this.hasAttribute(\"scene\")) scene = JSON.parse(this.getAttribute(\"scene\"));\n        else if (this.hasAttribute(\"src\")) scene = await (await fetch(this.getAttribute(\"src\"))).json();\n        else if (this.hasAttribute(\"model\")) {\n          if (!DF.sceneFor) throw new Error(\"The model catalogue is not loaded\");\n          scene = DF.sceneFor(this.getAttribute(\"model\"));\n        } else throw new Error(\"Give a scene, src or model attribute\");\n      } catch (e) { this.showError(e); return; }\n      if (this.hasAttribute(\"theme\")) scene.style = Object.assign({}, scene.style, { theme: this.getAttribute(\"theme\") });\n      if (this.hasAttribute(\"view\")) scene.view = Object.assign({}, scene.view, { type: this.getAttribute(\"view\") });\n      if (this.hasAttribute(\"title\")) scene.overlay = Object.assign({}, scene.overlay, { title: this.getAttribute(\"title\") });\n      if (this.player) { this.player.dispose(); this.player = null; }\n      try {\n        this.player = new DF.Player(this._stage, scene);\n        if (this.player.error) throw this.player.error;\n      } catch (e) { this.showError(e); return; }\n      this.player.on(\"state\", (s) => { if (this._playBtn) { this._playBtn.textContent = s === \"play\" ? \"❚❚\" : \"▶\"; this._playBtn.title = s === \"play\" ? \"Pause\" : \"Play\"; } });\n      this.player.on(\"error\", (e) => this.showError(e));\n      const reduced = window.matchMedia && window.matchMedia(\"(prefers-reduced-motion: reduce)\").matches;\n      if (this.hasAttribute(\"static\") || reduced) this.player.advance(parseInt(this.getAttribute(\"static\"), 10) || 240);\n      else if (!this.hasAttribute(\"paused\")) this.player.play();\n      this.dispatchEvent(new CustomEvent(\"ready\", { detail: this.player }));\n    }\n    showError(e) {\n      const d = document.createElement(\"div\"); d.className = \"err\";\n      d.textContent = \"DynFlow: \" + (e && e.message ? e.message : String(e));\n      this._stage.innerHTML = \"\"; this._stage.appendChild(d);\n    }\n    get scene() { return this.player ? this.player.getScene() : null; }\n    play() { if (this.player) this.player.play(); }\n    pause() { if (this.player) this.player.pause(); }\n    restart(seed) { if (this.player) { this.player.restart(seed); if (!this.hasAttribute(\"paused\")) this.player.play(); } }\n    setParam(name, value) { if (this.player) this.player.setParam(name, value); }\n  }\n\n  if (!customElements.get(\"dyn-flow\")) customElements.define(\"dyn-flow\", DynFlowElement);\n  DF.DynFlowElement = DynFlowElement;\n})(globalThis.DynFlow = globalThis.DynFlow || {});\n\nglobalThis.DynFlow.VERSION = \"0.1.0\";\n";
