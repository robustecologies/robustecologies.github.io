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

  const VIEW_DEFAULTS = {
    flow: { n: 1500, life: [80, 320], spawn: "box", fade: 0.06, lineWidth: 1.1, alpha: 0.55, colorBy: "dominant" },
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
    if (s.view.life === undefined && s.view.type === "flow") s.view.life = d.life;
    if (s.view.spawn === undefined && s.view.type === "flow") s.view.spawn = d.spawn;
    if (s.view.rotate === undefined && d.rotate !== undefined) s.view.rotate = d.rotate;
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
    // Full ranges for every variable: scene, then declared, then a probe run.
    const need = sys.vars.filter(function (v) { return !(s.view.ranges && s.view.ranges[v]) && !sys.ranges[v]; });
    const probe = need.length ? DF.autoRanges(sys, params, init, { dt: dt, steps: discrete ? 3000 : Math.min(20000, Math.max(3000, Math.round(60 / dt))), perturbations: s.perturbations }) : {};
    this.fullRanges = sys.vars.map(function (v) { return (s.view.ranges && s.view.ranges[v]) || sys.ranges[v] || probe[v]; });
    let axes = (s.view.axes || []).map(function (v) { return sys.vars.indexOf(v); }).filter(function (i) { return i >= 0; });
    if (axes.length < 2) axes = sys.vars.length >= 3 && (s.view.type === "trajectory" || s.view.type === "flow") && s.view.dim3 !== false ? [0, 1, 2] : sys.vars.length >= 2 ? [0, 1] : [0, 0];
    if (s.view.projection === "simplex" && sys.vars.length >= 3 && axes.length < 3) axes = [0, 1, 2];
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
        view: s.view, dt: h, params: params, init: init, ranges: this.fullRanges, axes: shown, spread: s.spread, initMode: s.initMode,
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
    const self = this;
    this.cam = new DF.Camera(axes, axes.map(function (i) { return self.fullRanges[i]; }), { azim: s.view.azim, elev: s.view.elev, zoom: s.view.zoom, pad: s.view.pad, upAxis: s.view.upAxis, simplex: s.view.projection === "simplex" });
    this.theme = DF.THEMES[s.style.theme] || DF.THEMES["relab-night"];
    this.palette = (DF.PALETTES[s.style.palette] || DF.PALETTES.relab).colors;
    this.acc = 0; this.nCap = Infinity; this.slow = 1;
    this.view = Object.create(DF.VIEWS[s.view.type]);
    this.userParam = false;
    this.frameCount = 0;
    this.speedMax = 1e-9;
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
  DF.normalizeScene = normalizeScene;
  DF.Player = Player;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
