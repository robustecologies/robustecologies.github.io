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
