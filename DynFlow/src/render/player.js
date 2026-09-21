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
