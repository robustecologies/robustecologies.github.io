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
