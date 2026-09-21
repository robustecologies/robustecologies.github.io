// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* DynFlow studio: choose or write a model, perturb it, choose a view and a
   style, and export the scene as an image, a film, a page or an element. */
(function () {
  "use strict";
  const DF = window.DynFlow;
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
  const h = (tag, attrs, ...kids) => {
    const e = document.createElement(tag);
    for (const k in attrs || {}) {
      const v = attrs[k];
      if (k === "class") e.className = v;
      else if (k === "html") e.innerHTML = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (v === true) e.setAttribute(k, "");
      else if (v !== false && v !== null && v !== undefined) e.setAttribute(k, v);
    }
    kids.flat().forEach((c) => { if (c !== null && c !== undefined && c !== false) e.append(c.nodeType ? c : document.createTextNode(String(c))); });
    return e;
  };
  const tex = (el, src, display) => { try { if (window.katex) window.katex.render(src, el, { throwOnError: false, displayMode: !!display }); else el.textContent = src; } catch (e) { el.textContent = src; } return el; };
  const clone = (o) => JSON.parse(JSON.stringify(o));
  const fmt = (v) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? (+v).toExponential(3) : String(+(+v).toPrecision(6)));

  const ICON = {
    play: '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M7 4.5v15l13-7.5z"/></svg>',
    pause: '<svg viewBox="0 0 24 24" fill="currentColor"><rect x="6" y="4.5" width="4" height="15" rx="1"/><rect x="14" y="4.5" width="4" height="15" rx="1"/></svg>',
    restart: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 12a8 8 0 1 0 2.4-5.7"/><path d="M4 4v5h5"/></svg>',
    step: '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M6 5v14l9-7z"/><rect x="16" y="5" width="2.5" height="14" rx="1"/></svg>',
    dice: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8"><rect x="4" y="4" width="16" height="16" rx="3"/><circle cx="9" cy="9" r="1.2" fill="currentColor"/><circle cx="15" cy="15" r="1.2" fill="currentColor"/><circle cx="15" cy="9" r="1.2" fill="currentColor"/><circle cx="9" cy="15" r="1.2" fill="currentColor"/></svg>',
    full: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/></svg>',
    search: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="6.5"/><path d="m16 16 4.5 4.5" stroke-linecap="round"/></svg>',
    export: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M12 15V3M7 8l5-5 5 5M5 14v5h14v-5"/></svg>',
    link: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M10 14a4 4 0 0 0 5.7 0l3-3a4 4 0 0 0-5.7-5.7l-1 1"/><path d="M14 10a4 4 0 0 0-5.7 0l-3 3a4 4 0 0 0 5.7 5.7l1-1"/></svg>',
    help: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="9"/><path d="M9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.7.3-1 .9-1 1.7M12 17h.01"/></svg>',
    lib: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 6h16M4 12h16M4 18h10"/></svg>',
    sliders: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 7h10M18 7h2M4 17h4M12 17h8"/><circle cx="16" cy="7" r="2"/><circle cx="10" cy="17" r="2"/></svg>'
  };
  const VIEW_ICON = {
    flow: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M3 18c5-12 12-12 26-10"/><path d="M3 12c7-8 14-7 26-3"/><path d="M5 22c6-8 12-9 24-8"/></svg>',
    trajectory: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.4"><path d="M16 12c-6-9-13-4-9 2s10 3 9-2c1 5 8 8 11 2s-3-11-11-2"/></svg>',
    timeseries: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.4"><path d="M2 20h28M3 14c3-10 5 8 8-2s5 9 8 0 5-8 8 1"/></svg>',
    phase: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M5 5l3 2M13 4l3 1M22 5l3 0M5 13l2 3M24 14l-1 3M9 20l3-1M19 20l3-2"/><circle cx="16" cy="12" r="2.2" fill="currentColor"/></svg>',
    sweep: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.4"><path d="M3 19c8-1 12-2 14-6-2-5-8-5-6 0 3 5 8 3 18 1"/></svg>',
    orbit: '<svg viewBox="0 0 32 24" fill="currentColor"><circle cx="4" cy="12" r="1"/><circle cx="9" cy="12" r="1"/><circle cx="13" cy="8" r="1"/><circle cx="13" cy="16" r="1"/><circle cx="17" cy="6" r="1"/><circle cx="17" cy="10" r="1"/><circle cx="17" cy="15" r="1"/><circle cx="17" cy="18" r="1"/><circle cx="21" cy="5" r=".8"/><circle cx="22" cy="9" r=".8"/><circle cx="22" cy="14" r=".8"/><circle cx="23" cy="19" r=".8"/><circle cx="26" cy="7" r=".8"/><circle cx="27" cy="12" r=".8"/><circle cx="27" cy="17" r=".8"/></svg>',
    density: '<svg viewBox="0 0 32 24"><rect x="2" y="3" width="28" height="18" rx="2" fill="none" stroke="currentColor" stroke-width="1.2"/><ellipse cx="11" cy="12" rx="5" ry="4" fill="currentColor" opacity=".6"/><ellipse cx="22" cy="12" rx="4" ry="3" fill="currentColor" opacity=".3"/></svg>',
    strobe: '<svg viewBox="0 0 32 24" fill="currentColor"><circle cx="6" cy="14" r="1"/><circle cx="9" cy="10" r="1"/><circle cx="13" cy="8" r="1"/><circle cx="17" cy="9" r="1"/><circle cx="20" cy="12" r="1"/><circle cx="23" cy="15" r="1"/><circle cx="26" cy="13" r="1"/><circle cx="11" cy="15" r=".8"/><circle cx="16" cy="16" r=".8"/></svg>',
    cobweb: '<svg viewBox="0 0 32 24" fill="none" stroke="currentColor" stroke-width="1.3"><path d="M3 21 29 3" stroke-dasharray="2 2"/><path d="M3 21c6-22 20-22 26 0"/><path d="M8 21V10h6v4h-6"/></svg>'
  };
  const VIEW_LABEL = { flow: "Flow", trajectory: "Trajectory", timeseries: "Time series", phase: "Phase plane", sweep: "Sweep", orbit: "Bifurcation", density: "Density", strobe: "Strobe", cobweb: "Cobweb" };

  // ------------------------------------------------------------ showcase
  const SHOWCASE = [
    { label: "Transiency begets persistence", id: "may-leonard", over: { view: { type: "flow", projection: "simplex", spawn: "mixed", life: [140, 460] }, n: 2200, overlay: { title: "Robust Ecologies Lab", subtitle: "Exploring how transiency begets persistence", caption: "Three competing species (May-Leonard): orbits leave coexistence and linger ever longer near each single-species state. Click to seed new initial conditions." } } },
    { label: "Fold catastrophe", id: "cusp", over: { style: { theme: "relab-night" }, overlay: { title: "Fold catastrophe", subtitle: "Hysteresis under a slow sweep of r", equations: true, readout: true } } },
    { label: "Lorenz attractor", id: "lorenz", over: { overlay: { title: "Lorenz system", subtitle: "sigma = 10, rho = 28, beta = 8/3", equations: true } } },
    { label: "Snapshot attractor", id: "zaslavsky", over: { overlay: { title: "A snapshot attractor", subtitle: "6000 states under one random forcing" } } },
    { label: "Rate-induced tipping", id: "r-tipping", over: { overlay: { title: "Rate-induced tipping", subtitle: "Ramp faster than r = 1 and the state escapes" } } },
    { label: "Kramers escape", id: "double-well", over: { overlay: { title: "Noise-driven escape", subtitle: "Ensemble density in a double well" } } },
    { label: "Duffing section", id: "duffing", over: { style: { theme: "blackboard", palette: "blackboard" }, overlay: { title: "Forced Duffing oscillator", subtitle: "Stroboscopic section" } } },
    { label: "Route to chaos", id: "logistic", over: { overlay: { title: "Logistic map", subtitle: "Period doubling to chaos" } } },
    { label: "Paradox of enrichment", id: "rosenzweig-macarthur", over: { overlay: { title: "Rosenzweig-MacArthur", subtitle: "Nullclines, equilibria and orbits", equations: true } } },
    { label: "Labyrinth chaos", id: "thomas", over: { overlay: { title: "Thomas attractor" } } },
    { label: "Climate change attractor", id: "lorenz84-forced", over: { overlay: { title: "Lorenz-84 under a changing climate", subtitle: "One forcing, 1500 members" } } },
    { label: "De Jong attractor", id: "de-jong", over: { overlay: { title: "" } } }
  ];

  const S = { scene: null, model: null, player: null, tab: "model", aspect: "fill", dirty: false };

  // ----------------------------------------------------------- layout
  function buildTop() {
    const top = $("#top");
    top.append(
      h("button", { class: "tbtn iconbtn only-narrow", title: "Models", "aria-label": "Show the model library", onclick: () => document.body.classList.toggle("show-lib"), html: ICON.lib }),
      h("a", { class: "brand", href: "#", onclick: (e) => { e.preventDefault(); loadShowcase(0); } },
        h("img", { src: "vendor/relab-logo.png", alt: "" }),
        h("span", { class: "wordmark", html: "<b>Dyn</b>Flow<small>studio</small>" })),
      h("input", { class: "scene-name", id: "sceneName", "aria-label": "Scene name", oninput: (e) => { S.scene.name = e.target.value; } }),
      h("div", { class: "spacer" }),
      h("button", { class: "tbtn hide-narrow", onclick: () => newModel(), title: "Write a model from formulas" }, "New model"),
      h("button", { class: "tbtn hide-narrow", onclick: openFile }, "Open"),
      h("button", { class: "tbtn hide-narrow", onclick: shareLink, html: ICON.link + "<span>Share</span>" }),
      h("button", { class: "tbtn primary", id: "exportBtn", onclick: (e) => { e.stopPropagation(); $("#exportMenu").classList.toggle("on"); }, html: ICON.export + "<span>Export</span>" }),
      h("button", { class: "tbtn iconbtn", title: "Help", "aria-label": "Help", onclick: () => $("#help").classList.add("on"), html: ICON.help }),
      h("button", { class: "tbtn iconbtn only-narrow", title: "Inspector", "aria-label": "Show the inspector", onclick: () => document.body.classList.toggle("show-insp"), html: ICON.sliders })
    );
    const m = $("#exportMenu");
    const item = (label, key, fn) => h("button", { onclick: () => { m.classList.remove("on"); fn(); } }, label, key ? h("small", {}, key) : null);
    m.append(
      item("PNG image", "P", () => exportPNG()),
      item("SVG image", null, () => exportSVG()),
      h("hr"),
      item("WebM video (8 s)", null, () => exportVideo("webm", 8)),
      item("GIF animation (5 s)", null, () => exportVideo("gif", 5)),
      h("hr"),
      item("Standalone HTML page", null, () => exportHTML()),
      item("Copy embed snippet", null, () => copyEmbed()),
      item("Scene file (JSON)", "S", () => exportJSON()),
      item("Copy share link", null, () => shareLink())
    );
    document.addEventListener("click", (e) => { if (!m.contains(e.target)) m.classList.remove("on"); });
  }

  function buildLibrary() {
    const lib = $("#library");
    const search = h("input", { type: "search", placeholder: "Search 80 models", "aria-label": "Search models", oninput: (e) => filterLibrary(e.target.value) });
    lib.append(h("div", { class: "lib-head" }, h("div", { class: "search", html: ICON.search }, search),
      h("button", { class: "btn block", onclick: () => newModel() }, "+ Write a model")));
    const list = h("div", { class: "lib-list", id: "libList" });
    lib.append(list);
    const sc = h("details", { class: "grp", open: true }, h("summary", {}, "Showcase", h("span", { class: "count" }, SHOWCASE.length)));
    const grid = h("div", { class: "showcase" });
    SHOWCASE.forEach((s, i) => {
      const cv = h("canvas", { width: 240, height: 128, "data-i": i });
      grid.append(h("button", { class: "card", title: s.label, onclick: () => loadShowcase(i) }, cv, h("span", {}, s.label)));
    });
    sc.append(grid); list.append(sc);
    DF.catalogueGroups().forEach((g) => {
      const items = DF.CATALOGUE.filter((m) => m.group === g);
      const d = h("details", { class: "grp", open: true, "data-group": g }, h("summary", {}, g, h("span", { class: "count" }, items.length)));
      items.forEach((m) => {
        const sys = DF.compileSystem(m.system);
        const tag = { ode: "ODE", sde: "SDE", map: "map", dde: "DDE" }[sys.kind] + " · " + sys.vars.length + "D" + (sys.usesTime ? " · forced" : "") + (sys.usesRandom ? " · random" : "");
        d.append(h("button", { class: "item", "data-id": m.id, "data-text": (m.name + " " + m.about + " " + m.group + " " + tag).toLowerCase(), onclick: () => loadModel(m.id) },
          h("span", { class: "nm" }, m.name), h("span", { class: "tg" }, tag)));
      });
      list.append(d);
    });
    renderThumbs();
  }
  function filterLibrary(q) {
    q = q.trim().toLowerCase();
    $$(".item").forEach((it) => { it.hidden = q && !it.dataset.text.includes(q); });
    $$(".grp[data-group]").forEach((g) => { g.hidden = !$$(".item", g).some((it) => !it.hidden); if (q) g.open = true; });
  }
  // Small previews of the showcase scenes, rendered once each, lazily.
  function renderThumbs() {
    const io = new IntersectionObserver((entries) => {
      entries.forEach((en) => {
        if (!en.isIntersecting) return;
        io.unobserve(en.target);
        const i = +en.target.dataset.i;
        // Previews render one at a time when the page is idle, after the live figure has started.
        queueThumb(() => thumb(en.target, SHOWCASE[i]));
      });
    });
    $$(".showcase canvas").forEach((c) => io.observe(c));
  }
  const thumbQueue = [];
  let thumbBusy = false;
  function queueThumb(fn) { thumbQueue.push(fn); pumpThumbs(); }
  function pumpThumbs() {
    if (thumbBusy || !thumbQueue.length) return;
    thumbBusy = true;
    const idle = window.requestIdleCallback || ((f) => setTimeout(f, 120));
    setTimeout(() => idle(() => { const fn = thumbQueue.shift(); try { fn(); } finally { thumbBusy = false; pumpThumbs(); } }, { timeout: 2000 }), 350);
  }
  function thumb(canvas, sc) {
    const host = h("div", { style: "position:fixed;left:-10000px;top:0;width:240px;height:128px" });
    document.body.append(host);
    try {
      const s = sceneOf(sc);
      s.n = Math.min(s.n || 400, 400); s.overlay = { position: "none" };
      const p = new DF.Player(host, s);
      p.advance(s.view && s.view.type === "orbit" ? 30 : 110);
      const c = p.composite();
      canvas.getContext("2d").drawImage(c, 0, 0, canvas.width, canvas.height);
      p.dispose();
    } catch (e) { /* a preview is optional */ }
    host.remove();
  }
  function sceneOf(sc) {
    const base = DF.sceneFor(sc.id);
    const over = clone(sc.over || {});
    ["view", "style", "overlay"].forEach((k) => { if (over[k]) { base[k] = Object.assign({}, base[k] || {}, over[k]); delete over[k]; } });
    return Object.assign(base, over);
  }

  // ------------------------------------------------------------- stage
  function buildStage() {
    const wrap = $("#stageWrap");
    const area = h("div", { class: "stage-area", id: "stageArea" });
    const stage = h("div", { class: "stage", id: "stage" });
    area.append(stage, h("div", { class: "err-bar", id: "errBar" }), h("div", { class: "stage-hint", id: "hint" }, "Click the figure to seed trajectories; drag to rotate 3D views"));
    const play = h("button", { class: "play", id: "playBtn", title: "Play or pause (space)", "aria-label": "Play or pause", onclick: () => S.player && S.player.toggle(), html: ICON.pause });
    const speed = h("input", { type: "range", class: "speed", id: "speed", min: 0, max: 100, "aria-label": "Speed", oninput: (e) => { const v = Math.max(1, Math.round(Math.pow(10, e.target.value / 40) - 0.5)); S.scene.stepsPerFrame = v; if (S.player) S.player.scene.stepsPerFrame = v; $("#speedVal").textContent = v + "×"; } });
    const seed = h("input", { class: "seed", id: "seed", type: "number", "aria-label": "Random seed", title: "Random seed", onchange: (e) => { S.scene.seed = +e.target.value || 1; restart(); } });
    const aspect = h("select", { "aria-label": "Aspect ratio", onchange: (e) => { S.aspect = e.target.value; fitStage(); } },
      ...["fill", "16:9", "4:3", "3:2", "1:1", "4:5", "9:16"].map((a) => h("option", { value: a }, a === "fill" ? "Fill" : a)));
    const tr = h("div", { class: "transport" },
      play,
      h("button", { class: "tool", title: "Restart (R)", "aria-label": "Restart", onclick: () => restart(), html: ICON.restart }),
      h("button", { class: "tool", title: "Step one frame (S)", "aria-label": "Step", onclick: () => { S.player.pause(); S.player.step(); }, html: ICON.step }),
      h("button", { class: "tool", title: "New random seed (N)", "aria-label": "New seed", onclick: () => newSeed(), html: ICON.dice }),
      seed,
      h("span", { class: "tlabel hide-narrow" }, "Speed"), speed, h("span", { class: "tlabel", id: "speedVal" }, ""),
      h("div", { class: "spacer" }),
      h("span", { class: "readout", id: "readout" }, ""),
      aspect,
      h("button", { class: "tool", title: "Full screen (F)", "aria-label": "Full screen", onclick: () => { const st = $("#stage"); if (document.fullscreenElement) document.exitFullscreen(); else st.requestFullscreen && st.requestFullscreen(); }, html: ICON.full })
    );
    wrap.append(area, tr);
    new ResizeObserver(fitStage).observe(area);
  }
  function fitStage() {
    const area = $("#stageArea"), st = $("#stage");
    if (S.aspect === "fill") { st.style.width = "100%"; st.style.height = "100%"; return; }
    const [a, b] = S.aspect.split(":").map(Number);
    const W = area.clientWidth - 36, H = area.clientHeight - 36;
    const w = Math.min(W, H * a / b);
    st.style.width = Math.floor(w) + "px"; st.style.height = Math.floor(w * b / a) + "px";
  }

  // --------------------------------------------------------- scene flow
  function loadShowcase(i) { S.model = SHOWCASE[i].id; setScene(sceneOf(SHOWCASE[i])); }
  function loadModel(id) {
    const s = DF.sceneFor(id);
    if (S.scene && S.scene.style && $("#keepTheme") && $("#keepTheme").checked) s.style = Object.assign({}, s.style || {}, { theme: S.scene.style.theme });
    S.model = id; setScene(s);
    document.body.classList.remove("show-lib");
  }
  function setScene(scene, keepTab) {
    const s = DF.normalizeScene(scene);
    if (!S.player) {
      S.player = new DF.Player($("#stage"), s);
      S.player.on("state", (st) => { $("#playBtn").innerHTML = st === "play" ? ICON.pause : ICON.play; });
      S.player.on("error", (e) => showError(e));
      S.player.on("frame", onFrame);
      S.player.on("param", (d) => syncParam(d.name, d.value));
      S.player.on("interact", () => { $("#hint").style.opacity = 0; });
    } else S.player.load(s);
    if (S.player.error) { showError(S.player.error); return; }
    hideError();
    S.scene = S.player.scene;
    S.model = s.model || S.model;
    $("#sceneName").value = S.scene.name || "";
    $("#seed").value = S.scene.seed;
    const spf = S.scene.stepsPerFrame;
    $("#speed").value = Math.round(40 * Math.log10(spf + 0.5)); $("#speedVal").textContent = spf + "×";
    $$(".item").forEach((it) => it.classList.toggle("active", it.dataset.id === S.model));
    renderInspector();
    S.player.play();
    save();
  }
  // Rebuild the player from the current scene (after structural changes).
  function restart() {
    if (!S.player) return;
    const wasRunning = S.player.running;
    S.player.load(clone(S.scene));
    if (S.player.error) { showError(S.player.error); return; }
    hideError(); S.scene = S.player.scene;
    if (wasRunning || !S.player.running) S.player.play();
    save();
  }
  function newSeed() { S.scene.seed = 1 + Math.floor(Math.random() * 99999); $("#seed").value = S.scene.seed; restart(); }
  let lastRead = 0;
  function onFrame() {
    const now = performance.now();
    if (now - lastRead < 200) return;
    lastRead = now;
    const p = S.player, sim = p.sim;
    $("#readout").textContent = (p.sys.time === "discrete" ? "n " + sim.t : "t " + sim.t.toFixed(2)) + "  " + Math.round(p.fps || 0) + " fps";
    if (S.tab === "perturb") $$("[data-live-param]").forEach((el) => { const i = sim.paramIndex(el.dataset.liveParam); if (i >= 0) el.textContent = fmt(sim.p[i]); });
  }
  function showError(e) { const b = $("#errBar"); b.textContent = e && e.message ? e.message : String(e); b.classList.add("on"); }
  function hideError() { $("#errBar").classList.remove("on"); }
  let saveTimer = 0;
  function save() {
    clearTimeout(saveTimer);
    saveTimer = setTimeout(() => { try { localStorage.setItem("dynflow:scene", JSON.stringify(S.player.getScene())); } catch (e) { /* storage unavailable */ } }, 400);
  }

  // ---------------------------------------------------------- inspector
  const TABS = [["model", "Model"], ["perturb", "Perturb"], ["view", "View"], ["style", "Style"], ["export", "Export"]];
  function buildInspector() {
    const ins = $("#inspector");
    const tabs = h("div", { class: "tabs", role: "tablist" });
    TABS.forEach(([k, label]) => tabs.append(h("button", { class: "tab", role: "tab", "aria-selected": k === S.tab ? "true" : "false", "data-tab": k, onclick: () => { S.tab = k; renderInspector(); } }, label)));
    ins.append(tabs);
    TABS.forEach(([k]) => ins.append(h("div", { class: "pane", id: "pane-" + k, role: "tabpanel", hidden: k !== S.tab })));
  }
  function renderInspector() {
    $$(".tab").forEach((t) => t.setAttribute("aria-selected", t.dataset.tab === S.tab ? "true" : "false"));
    TABS.forEach(([k]) => { $("#pane-" + k).hidden = k !== S.tab; });
    ({ model: renderModel, perturb: renderPerturb, view: renderView, style: renderStyle, export: renderExport })[S.tab]();
  }
  const sec = (title, act, ...kids) => h("div", { class: "sec" }, h("h3", {}, title, act ? h("span", { class: "act" }, act) : null), ...kids);
  function numInput(value, onchange, opts) {
    opts = opts || {};
    return h("input", { type: "number", value: fmt(value), step: opts.step || "any", min: opts.min, max: opts.max, "aria-label": opts.label || "", onchange: (e) => { const v = parseFloat(e.target.value); if (isFinite(v)) onchange(v); } });
  }

  // Model
  function renderModel() {
    const pane = $("#pane-model"); pane.innerHTML = "";
    const p = S.player, sys = p.sys, m = DF.CATALOGUE.find((e) => e.id === S.model && e.system === S.scene.system);
    if (m) pane.append(sec(m.name, null, h("p", { class: "about" }, m.about), h("p", { class: "source" }, "Source: " + m.source)));
    const eq = h("div", { class: "eqs" });
    DF.systemLatex(S.scene.system).lines.forEach((L) => eq.append(tex(h("div"), L, true)));
    pane.append(sec("Equations", h("button", { class: "btn small", onclick: () => { const d = $("#formulaBox"); d.open = !d.open; } }, "Edit"), eq));
    const ta = h("textarea", { class: "code", id: "formula", spellcheck: "false", rows: Math.min(22, S.scene.system.split("\n").length + 2) });
    ta.value = S.scene.system;
    const err = h("div", { class: "formula-err", id: "formulaErr" });
    const apply = () => {
      try { DF.compileSystem(ta.value); } catch (e) { err.textContent = e.message; err.className = "formula-err"; return; }
      err.textContent = "Compiled"; err.className = "formula-err formula-ok";
      S.scene.system = ta.value; S.scene.params = {}; S.scene.init = {};
      if (S.scene.view) delete S.scene.view.ranges;
      const sys2 = DF.compileSystem(ta.value);
      if (!DF.viewsFor(sys2).includes(S.scene.view.type)) S.scene.view.type = sys2.time === "discrete" ? "flow" : "trajectory";
      if (S.scene.view.axes && !S.scene.view.axes.every((v) => sys2.vars.includes(v))) delete S.scene.view.axes;
      S.model = null; restart(); renderInspector();
    };
    ta.addEventListener("keydown", (e) => { if ((e.ctrlKey || e.metaKey) && e.key === "Enter") { e.preventDefault(); apply(); } });
    const box = h("details", { id: "formulaBox", class: "sec" }, h("summary", { class: "hint", style: "cursor:pointer" }, "Formula source"), ta, err,
      h("div", { class: "btns" }, h("button", { class: "btn primary", onclick: apply }, "Apply"), h("span", { class: "hint" }, "Ctrl + Enter; see Help for the language")));
    pane.append(box);

    if (sys.params.length) {
      const list = h("div");
      sys.params.forEach((q, i) => {
        const val = p.sim.base[i];
        const lo = Math.min(q.min, val), hi = Math.max(q.max, val);
        const range = h("input", { type: "range", min: lo, max: hi, step: "any", value: val, "data-param": q.name });
        const num = numInput(val, (v) => { range.value = v; setParam(q.name, v); }, { label: q.name });
        num.dataset.paramNum = q.name;
        range.addEventListener("input", () => { num.value = fmt(+range.value); setParam(q.name, +range.value); });
        list.append(h("div", { class: "row" }, tex(h("label", { title: q.name }), DF.texName(q.name)), range, num));
      });
      pane.append(sec("Parameters", h("button", { class: "btn small", title: "Back to the model's values", onclick: () => { S.scene.params = {}; restart(); renderInspector(); } }, "Reset"), list));
    }
    const ic = h("div");
    sys.vars.forEach((v, i) => ic.append(h("div", { class: "row two" }, tex(h("label", { title: v }), DF.texName(v) + "(0)"), numInput(p.sim.init[i], (x) => { S.scene.init[v] = x; restart(); }, { label: v }))));
    pane.append(sec("Initial state", null, ic));
    const kindName = { ode: "Runge-Kutta 4", sde: "Euler-Maruyama (Ito)", dde: "Runge-Kutta 4 with Hermite history", map: "iteration" }[sys.kind];
    pane.append(sec("Integration", null,
      h("p", { class: "hint", style: "margin:0 0 6px" }, "Scheme: " + kindName + (sys.time === "discrete" ? "" : ", fixed step")),
      sys.time === "continuous" ? h("div", { class: "row two" }, h("label", {}, "Step dt"), numInput(S.scene.dt, (v) => { S.scene.dt = Math.max(1e-6, v); restart(); }, { label: "dt" })) : null,
      h("div", { class: "row two" }, h("label", {}, "Steps per frame"), numInput(S.scene.stepsPerFrame, (v) => { S.scene.stepsPerFrame = Math.max(1, Math.round(v)); S.player.scene.stepsPerFrame = S.scene.stepsPerFrame; }, { label: "steps per frame" })),
      h("div", { class: "row two" }, h("label", {}, "Ensemble size"), numInput(S.scene.n, (v) => { S.scene.n = Math.max(1, Math.min(20000, Math.round(v))); restart(); }, { label: "ensemble size" })),
      h("div", { class: "row two" }, h("label", {}, "Initial spread"), h("select", { onchange: (e) => { S.scene.initMode = e.target.value; restart(); } },
        ...[["point", "Around the initial state"], ["ball", "Gaussian ball"], ["box", "Uniform over the axes"]].map(([v, t]) => h("option", { value: v, selected: S.scene.initMode === v }, t)))),
      h("div", { class: "row two" }, h("label", {}, "Spread"), numInput(S.scene.spread, (v) => { S.scene.spread = Math.max(0, v); restart(); }, { label: "spread" })),
      h("label", { class: "check" }, h("input", { type: "checkbox", checked: !!S.scene.keepPositive, onchange: (e) => { S.scene.keepPositive = e.target.checked; restart(); } }), "Keep states non-negative")
    ));
  }
  function setParam(name, v) { S.scene.params[name] = v; S.player.setParam(name, v); save(); }
  function syncParam(name, v) {
    if (S.tab !== "model") return;
    const r = $('input[data-param="' + name + '"]'), n = $('input[data-param-num="' + name + '"]');
    if (r) r.value = v; if (n && document.activeElement !== n) n.value = fmt(v);
  }

  // Perturb
  const PKIND = {
    periodic: { label: "Periodic forcing", target: "param", fields: [["amp", "Amplitude", 0.1], ["period", "Period", 10], ["phase", "Phase", 0]] },
    quasiperiodic: { label: "Quasiperiodic forcing", target: "param", fields: [["amp", "Amplitude 1", 0.1], ["period", "Period 1", 10], ["amp2", "Amplitude 2", 0.1], ["period2", "Period 2", 16.18]] },
    ramp: { label: "Ramp (rate-induced)", target: "param", fields: [["rate", "Rate", 0.01], ["t0", "Start", 0], ["span", "Total change", 1]] },
    step: { label: "Step change", target: "param", fields: [["amp", "Jump", 0.5], ["t0", "At time", 10]] },
    ou: { label: "Coloured noise on a parameter", target: "param", fields: [["sigma", "Std deviation", 0.1], ["tau", "Correlation time", 1]] },
    additive: { label: "Additive noise", target: "var", fields: [["sigma", "Intensity", 0.1]], common: true },
    multiplicative: { label: "Multiplicative noise", target: "var", fields: [["sigma", "Intensity", 0.1]], common: true },
    coloured: { label: "Coloured (OU) noise", target: "var", fields: [["sigma", "Std deviation", 0.2], ["tau", "Correlation time", 1]], common: true },
    levy: { label: "Levy (alpha-stable) noise", target: "var", fields: [["sigma", "Scale", 0.05], ["alpha", "Alpha (0, 2]", 1.5]], common: true },
    jumps: { label: "Random jumps", target: "var", fields: [["rate", "Rate", 0.2], ["size", "Size", 0.5], ["sd", "Size std", 0], ["frac", "Or fraction lost", ""]], common: true },
    pulse: { label: "Periodic pulses (harvest, kick)", target: "var", fields: [["period", "Period", 10], ["t0", "First at", 10], ["size", "Size", 0], ["frac", "Or fraction lost", 0.3]] }
  };
  function renderPerturb() {
    const pane = $("#pane-perturb"); pane.innerHTML = "";
    const sys = S.player.sys, list = S.scene.perturbations;
    const add = h("select", { "aria-label": "Add a perturbation", onchange: (e) => {
      const k = e.target.value; if (!k) return;
      const def = PKIND[k], q = { kind: k };
      if (def.target === "param") { if (!sys.params.length) { toast("This model has no parameters"); e.target.value = ""; return; } q.param = sys.params[0].name; }
      else q.var = sys.vars[0];
      def.fields.forEach(([f, , d]) => { if (d !== "") q[f] = d; });
      if (k === "ramp") { const pr = sys.params[0]; q.rate = +((pr.max - pr.min) / 100).toPrecision(2); q.span = +((pr.max - pr.min) / 2).toPrecision(2); }
      list.push(q); restart(); renderPerturb();
    } },
      h("option", { value: "" }, "Add a perturbation..."),
      h("optgroup", { label: "Forcing of a parameter" }, ...["periodic", "quasiperiodic", "ramp", "step", "ou"].map((k) => h("option", { value: k }, PKIND[k].label))),
      h("optgroup", { label: "Noise and shocks on a state" }, ...["additive", "multiplicative", "coloured", "levy", "jumps", "pulse"].map((k) => h("option", { value: k }, PKIND[k].label))));
    pane.append(sec("Perturbations", null, add, h("p", { class: "hint" }, "Parameter forcing makes the system nonautonomous; state noise makes it stochastic. Common noise drives the whole ensemble with one realisation, as an environment does, which reveals snapshot and pullback attractors.")));
    list.forEach((q, j) => {
      const def = PKIND[q.kind];
      const card = h("div", { class: "pcard" + (q.enabled === false ? " off" : "") });
      const live = (f, v) => { q[f] = v; const sp = S.player.sim.perturbations[j]; if (sp) { sp[f] = v; S.player.sim.updateParams(); } save(); };
      card.append(h("header", {}, h("input", { type: "checkbox", checked: q.enabled !== false, title: "Enabled", onchange: (e) => { q.enabled = e.target.checked; card.classList.toggle("off", !q.enabled); live("enabled", q.enabled); } }),
        h("strong", {}, def.label), def.target === "param" ? h("span", { class: "kind-tag", "data-live-param": q.param }, "") : null,
        h("button", { class: "x", title: "Remove", "aria-label": "Remove", onclick: () => { list.splice(j, 1); restart(); renderPerturb(); } }, "×")));
      const targets = def.target === "param" ? sys.params.map((p) => p.name) : sys.vars;
      card.append(h("div", { class: "row" }, h("label", {}, def.target === "param" ? "Parameter" : "Variable"),
        h("select", { onchange: (e) => { q[def.target] = e.target.value; restart(); renderPerturb(); } }, ...targets.map((t) => h("option", { value: t, selected: q[def.target] === t }, t)))));
      def.fields.forEach(([f, label]) => {
        const inp = h("input", { type: "number", step: "any", value: q[f] === undefined || q[f] === null ? "" : q[f], onchange: (e) => { const v = e.target.value === "" ? undefined : parseFloat(e.target.value); if (v === undefined) delete q[f]; if (v === undefined || isFinite(v)) { if (f === "period" || f === "t0" || f === "tau") { q[f] = v; restart(); } else live(f, v); } } });
        card.append(h("div", { class: "row" }, h("label", {}, label), inp));
      });
      if (def.common) card.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: !!q.common, onchange: (e) => live("common", e.target.checked) }), "Common to the ensemble"));
      pane.append(card);
    });
  }

  // View
  function renderView() {
    const pane = $("#pane-view"); pane.innerHTML = "";
    const sys = S.player.sys, v = S.scene.view, ok = DF.viewsFor(sys);
    const grid = h("div", { class: "views" });
    Object.keys(VIEW_LABEL).forEach((k) => grid.append(h("button", { class: "vbtn", "aria-pressed": v.type === k ? "true" : "false", disabled: !ok.includes(k), title: DF.VIEWS[k].label, onclick: () => switchView(k), html: VIEW_ICON[k] + "<span>" + VIEW_LABEL[k] + "</span>" })));
    pane.append(sec("View", null, grid));
    const set = (k, val, noRestart) => { v[k] = val; if (!noRestart) restart(); save(); };
    const selectOf = (cur, opts, fn) => h("select", { onchange: (e) => fn(e.target.value) }, ...opts.map((o) => Array.isArray(o) ? h("option", { value: o[0], selected: String(cur) === String(o[0]) }, o[1]) : h("option", { value: o, selected: cur === o }, o)));
    const row = (label, el) => h("div", { class: "row two" }, h("label", {}, label), el);
    const opts = h("div");
    const axesViews = ["flow", "trajectory", "phase", "density", "strobe"];
    if (axesViews.includes(v.type) && sys.vars.length >= 2) {
      const cam = S.player.cam, cur = cam.axes.map((i) => sys.vars[i]);
      const three = cur.length === 3 && !cam.simplex;
      const setAxis = (i, val) => { const a = cur.slice(); a[i] = val; v.axes = a; restart(); renderView(); };
      opts.append(row("Horizontal", selectOf(cur[0], sys.vars, (x) => setAxis(0, x))), row("Vertical", selectOf(cur[1], sys.vars, (x) => setAxis(1, x))));
      if ((v.type === "flow" || v.type === "trajectory") && sys.vars.length >= 3) {
        opts.append(row("Projection", selectOf(cam.simplex ? "simplex" : three ? "3d" : "2d", [["2d", "Plane"], ["3d", "Rotating 3D"], ["simplex", "Simplex (shares)"]], (x) => {
          if (x === "2d") { v.axes = cur.slice(0, 2); v.dim3 = false; delete v.projection; }
          else if (x === "3d") { v.axes = cur.length === 3 ? cur : cur.concat(sys.vars.find((q) => !cur.includes(q))); v.dim3 = true; delete v.projection; }
          else { v.projection = "simplex"; v.axes = sys.vars.slice(0, 3); }
          restart(); renderView();
        })));
        if (three) {
          opts.append(row("Depth axis", selectOf(cur[2], sys.vars, (x) => setAxis(2, x))));
          const rot = h("input", { type: "range", min: 0, max: 1.5, step: 0.01, value: v.rotate || 0, oninput: (e) => { v.rotate = +e.target.value; } });
          opts.append(row("Rotation", rot));
        }
      }
    }
    const pnames = sys.params.map((q) => q.name);
    switch (v.type) {
      case "flow":
        opts.append(row("Particle life", h("div", { style: "display:flex;gap:6px" },
          numInput(Array.isArray(v.life) ? v.life[0] : 0, (x) => { v.life = [x, Array.isArray(v.life) ? v.life[1] : x * 3]; restart(); }, { label: "minimum life in frames" }),
          numInput(Array.isArray(v.life) ? v.life[1] : 0, (x) => { v.life = [Array.isArray(v.life) ? v.life[0] : x / 3, x]; restart(); }, { label: "maximum life in frames" }))));
        opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.life === "inf", onchange: (e) => { v.life = e.target.checked ? "inf" : [80, 320]; restart(); renderView(); } }), "Immortal particles (snapshot of the ensemble)"));
        opts.append(row("Births", selectOf(v.spawn || "box", [["box", "Uniform over the axes"], ["init", "Near the initial state"], ["mixed", "Mixed"]], (x) => set("spawn", x))));
        break;
      case "trajectory":
        opts.append(row("Tail length", numInput(v.tail || 2500, (x) => set("tail", Math.max(10, Math.round(x))), { label: "tail" })));
        opts.append(row("Warm-up frames", numInput(v.warmup || 0, (x) => set("warmup", Math.max(0, Math.round(x))), { label: "warm-up" })));
        break;
      case "timeseries":
        opts.append(row("Window", numInput(v.window || (sys.time === "discrete" ? 100 : 50), (x) => set("window", x), { label: "window" })));
        opts.append(row("Members shown", numInput(v.members || 1, (x) => { v.members = Math.max(1, Math.round(x)); S.scene.n = Math.max(S.scene.n, v.members); restart(); }, { label: "members" })));
        const vs = h("div", { style: "display:flex;flex-wrap:wrap;gap:4px 12px" });
        const cur = v.vars && v.vars.length ? v.vars : sys.vars.slice(0, 4);
        sys.vars.forEach((q) => vs.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: cur.includes(q), onchange: (e) => { const s2 = new Set(cur); if (e.target.checked) s2.add(q); else s2.delete(q); v.vars = sys.vars.filter((z) => s2.has(z)); restart(); } }), q)));
        opts.append(row("Variables", vs));
        break;
      case "phase":
        opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.field !== "none", onchange: (e) => set("field", e.target.checked ? "arrows" : "none") }), "Vector field"));
        opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.nullclines !== false, onchange: (e) => set("nullclines", e.target.checked) }), "Nullclines"));
        opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.equilibria !== false, onchange: (e) => set("equilibria", e.target.checked) }), "Equilibria and their stability"));
        opts.append(row("Arrow grid", numInput(v.fieldDensity || 22, (x) => set("fieldDensity", Math.max(6, Math.min(60, Math.round(x)))), { label: "arrow grid" })));
        opts.append(h("p", { class: "hint" }, "Click the plane to launch an orbit. With more than two variables the others are held at their initial values."));
        break;
      case "sweep": case "orbit": {
        const cur = v.param || pnames[0], q = sys.params.find((z) => z.name === cur) || sys.params[0];
        opts.append(row("Parameter", selectOf(cur, pnames, (x) => { v.param = x; const qq = sys.params.find((z) => z.name === x); v.from = qq.min; v.to = qq.max; restart(); renderView(); })));
        opts.append(row("Variable", selectOf(v.var || sys.vars[0], sys.vars, (x) => set("var", x))));
        opts.append(row("From", numInput(v.from === undefined ? q.min : v.from, (x) => set("from", x), { label: "from" })), row("To", numInput(v.to === undefined ? q.max : v.to, (x) => set("to", x), { label: "to" })));
        if (v.type === "sweep") {
          opts.append(row("Sweep speed", h("input", { type: "range", min: 0.00005, max: 0.004, step: 0.00005, value: v.speed || 0.0006, oninput: (e) => { v.speed = +e.target.value; } })));
          opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.branches !== false, onchange: (e) => set("branches", e.target.checked) }), "Equilibrium branches"));
          opts.append(h("p", { class: "hint" }, "Click the plot to set the parameter by hand; Restart resumes the sweep. A faster sweep delays the jump past the fold."));
        } else {
          opts.append(row("Transient", numInput(v.transient || (sys.time === "discrete" ? 300 : Math.round(200 / S.scene.dt)), (x) => set("transient", Math.round(x)), { label: "transient" })));
          opts.append(row("Samples", numInput(v.samples || (sys.time === "discrete" ? 150 : Math.round(400 / S.scene.dt)), (x) => set("samples", Math.round(x)), { label: "samples" })));
          opts.append(h("label", { class: "check" }, h("input", { type: "checkbox", checked: v.follow !== false, onchange: (e) => set("follow", e.target.checked) }), "Follow the attractor across parameter values"));
        }
        break;
      }
      case "density":
        opts.append(row("Mode", selectOf(v.mode || (sys.vars.length === 1 ? "carpet" : "map"), [["map", "Heat map of the ensemble"], ["carpet", "Density of one variable over time"]], (x) => set("mode", x))));
        opts.append(row("Carpet variable", selectOf(v.var || sys.vars[0], sys.vars, (x) => set("var", x))));
        opts.append(row("Memory", h("input", { type: "range", min: 0, max: 0.995, step: 0.005, value: v.decay === undefined ? 0.85 : v.decay, oninput: (e) => { v.decay = +e.target.value; } })));
        opts.append(row("Cell size (px)", numInput(v.cell || 3, (x) => set("cell", Math.max(1, x)), { label: "cell" })));
        break;
      case "strobe":
        opts.append(row("Section", selectOf(v.mode || "stroboscopic", [["stroboscopic", "Every period T"], ["section", "Poincare plane"]], (x) => { set("mode", x); renderView(); })));
        if (v.mode === "section") {
          opts.append(row("Plane variable", selectOf(v.sectionVar || sys.vars[sys.vars.length - 1], sys.vars, (x) => set("sectionVar", x))));
          opts.append(row("Plane value", numInput(v.sectionValue || 0, (x) => set("sectionValue", x), { label: "plane value" })));
        } else opts.append(row("Period T", numInput(v.period || 2 * Math.PI, (x) => set("period", x), { label: "period" })));
        opts.append(row("Transient", numInput(v.transient === undefined ? 20 : v.transient, (x) => set("transient", x), { label: "transient" })));
        break;
      case "cobweb":
        opts.append(row("Tail", numInput(v.tail || 120, (x) => set("tail", Math.round(x)), { label: "tail" })));
        break;
    }
    pane.append(sec("Options", null, opts));
    // Ranges
    const rg = h("div");
    sys.vars.forEach((q, i) => {
      const r = S.player.fullRanges[i];
      rg.append(h("div", { class: "row" }, tex(h("label"), DF.texName(q)),
        numInput(r[0], (x) => { v.ranges = Object.assign({}, v.ranges || {}); v.ranges[q] = [x, S.player.fullRanges[i][1]]; restart(); }, { label: q + " minimum" }),
        numInput(r[1], (x) => { v.ranges = Object.assign({}, v.ranges || {}); v.ranges[q] = [S.player.fullRanges[i][0], x]; restart(); }, { label: q + " maximum" })));
    });
    pane.append(sec("Axis ranges", h("button", { class: "btn small", title: "Measure the ranges from a probe run at the current parameters", onclick: () => fitRanges() }, "Fit"), rg,
      h("label", { class: "check" }, h("input", { type: "checkbox", checked: !!v.showAxes, onchange: (e) => set("showAxes", e.target.checked) }), "Show axes on flow and trajectory views")));
  }
  function switchView(k) {
    const v = S.scene.view, d = DF.VIEW_DEFAULTS[k];
    const keep = { axes: v.axes, ranges: v.ranges, projection: k === "flow" || k === "trajectory" ? v.projection : undefined };
    S.scene.view = Object.assign({ type: k }, keep);
    S.scene.n = d.n; S.scene.stepsPerFrame = Math.max(d.stepsPerFrame, k === "orbit" ? 1 : S.scene.stepsPerFrame);
    Object.assign(S.scene.style, { fade: d.fade, lineWidth: d.lineWidth, alpha: d.alpha, colorBy: d.colorBy });
    if (k === "sweep" || k === "orbit") { S.scene.view.param = S.player.sys.params[0].name; }
    restart(); renderView();
  }
  function fitRanges() {
    const p = S.player, r = DF.autoRanges(p.sys, Array.from(p.sim.base), Array.from(p.sim.init), { dt: S.scene.dt, steps: p.sys.time === "discrete" ? 3000 : Math.min(40000, Math.max(3000, Math.round(80 / S.scene.dt))), perturbations: S.scene.perturbations });
    S.scene.view.ranges = {};
    p.sys.vars.forEach((v) => { if (!p.sys.ranges[v]) S.scene.view.ranges[v] = r[v]; else S.scene.view.ranges[v] = r[v]; });
    restart(); renderView(); toast("Ranges fitted to a probe run");
  }

  // Style
  function renderStyle() {
    const pane = $("#pane-style"); pane.innerHTML = "";
    const st = S.scene.style, ov = S.scene.overlay;
    const th = h("div", { class: "themes" });
    Object.keys(DF.THEMES).forEach((k) => {
      const t = DF.THEMES[k];
      th.append(h("button", { class: "tswatch", "aria-pressed": st.theme === k ? "true" : "false", title: t.label, style: "background:" + (k === "transparent" ? "repeating-conic-gradient(#555 0% 25%, #333 0% 50%) 50% / 12px 12px" : DF.backgroundCSS(k)) + ";color:" + t.ink, onclick: () => { S.player.setStyle({ theme: k }); S.player.clearTrail(); renderStyle(); save(); } }, h("span", {}, t.label)));
    });
    pane.append(sec("Theme", null, th, h("label", { class: "check" }, h("input", { type: "checkbox", id: "keepTheme", checked: true }), "Keep the theme when changing model")));
    const palSel = h("select", { onchange: (e) => { S.player.setStyle({ palette: e.target.value }); renderStyle(); save(); } }, ...Object.keys(DF.PALETTES).map((k) => h("option", { value: k, selected: st.palette === k }, DF.PALETTES[k].label)));
    const pr = h("div", { class: "palrow" }, ...DF.PALETTES[st.palette].colors.map((c) => h("i", { style: "background:" + c })));
    const rampSel = h("select", { onchange: (e) => { S.player.setStyle({ ramp: e.target.value }); renderStyle(); save(); } }, ...Object.keys(DF.RAMPS).map((k) => h("option", { value: k, selected: st.ramp === k }, DF.RAMPS[k].label)));
    const gr = h("div", { class: "grad", style: "background:linear-gradient(90deg," + DF.RAMPS[st.ramp].stops.join(",") + ")" });
    const colorBy = h("select", { onchange: (e) => { S.player.setStyle({ colorBy: e.target.value }); renderStyle(); save(); } },
      ...[["solid", "One colour"], ["dominant", "Largest variable"], ["member", "Ensemble member"], ["speed", "Speed"], ["var", "Value of a variable"], ["time", "Position along the tail"], ["age", "Age of the particle"]].map(([k, t]) => h("option", { value: k, selected: st.colorBy === k }, t)));
    pane.append(sec("Colour", null,
      h("div", { class: "row two" }, h("label", {}, "Palette"), h("div", {}, palSel, pr)),
      h("div", { class: "row two" }, h("label", {}, "Ramp"), h("div", {}, rampSel, gr)),
      h("div", { class: "row two" }, h("label", {}, "Colour by"), colorBy),
      st.colorBy === "var" ? h("div", { class: "row two" }, h("label", {}, "Variable"), h("select", { onchange: (e) => { S.player.setStyle({ colorVar: e.target.value }); save(); } }, ...S.player.sys.vars.map((q) => h("option", { value: q, selected: st.colorVar === q }, q)))) : null));
    const slider = (label, key, min, max, step) => h("div", { class: "row" }, h("label", {}, label),
      h("input", { type: "range", min: min, max: max, step: step, value: st[key], oninput: (e) => { st[key] = +e.target.value; e.target.nextSibling && (e.target.parentNode.lastChild.textContent = fmt(+e.target.value)); save(); } }),
      h("span", { class: "lab" }, fmt(st[key])));
    pane.append(sec("Strokes", null,
      slider("Trail fade", "fade", 0, 0.6, 0.005), slider("Line width", "lineWidth", 0.2, 5, 0.1), slider("Opacity", "alpha", 0.02, 1, 0.01), slider("Point size", "pointSize", 0.4, 6, 0.1),
      h("div", { class: "row two" }, h("label", {}, "Resolution"), h("select", { onchange: (e) => { S.player.setStyle({ renderScale: +e.target.value }); restart(); } },
        ...[[1, "Screen"], [2, "2× (print)"], [3, "3× (poster)"]].map(([v, t]) => h("option", { value: v, selected: +st.renderScale === v }, t))))));
    const txt = (label, key, area) => h("div", { class: "row two" }, h("label", {}, label),
      h(area ? "textarea" : "input", { type: "text", rows: 2, value: ov[key] || "", oninput: (e) => { S.player.setOverlay({ [key]: e.target.value }); save(); } }));
    const tgl = (label, key) => h("label", { class: "check" }, h("input", { type: "checkbox", checked: !!ov[key], onchange: (e) => { S.player.setOverlay({ [key]: e.target.checked }); save(); } }), label);
    const capBox = txt("Caption", "caption", true); $("textarea", capBox) && ($("textarea", capBox).value = ov.caption || "");
    pane.append(sec("Text on the figure", null, txt("Title", "title"), txt("Subtitle", "subtitle"), capBox,
      tgl("Equations", "equations"), tgl("Time and forcing readout", "readout"), tgl("Legend", "legend"),
      h("label", { class: "check" }, h("input", { type: "checkbox", checked: ov.position === "none", onchange: (e) => { S.player.setOverlay({ position: e.target.checked ? "none" : "top-left" }); save(); } }), "Hide all text")));
  }

  // Export
  function renderExport() {
    const pane = $("#pane-export"); pane.innerHTML = "";
    const prog = h("div", { class: "progress", id: "prog" }, h("i"));
    const secs = h("input", { type: "number", value: 8, min: 1, max: 60, "aria-label": "seconds" }), fps = h("input", { type: "number", value: 30, min: 5, max: 60, "aria-label": "frames per second" });
    const gifW = h("input", { type: "number", value: 640, min: 120, max: 1600, "aria-label": "GIF width" });
    pane.append(sec("Image", null, h("div", { class: "btns" }, h("button", { class: "btn primary", onclick: exportPNG }, "PNG"), h("button", { class: "btn", onclick: exportSVG }, "SVG")),
      h("p", { class: "hint" }, "PNG at the current resolution (Style, Resolution for print). SVG keeps orbits, nullclines and branches as vector paths and embeds particle clouds as an image.")));
    pane.append(sec("Moving image", null,
      h("div", { class: "row two" }, h("label", {}, "Seconds"), secs), h("div", { class: "row two" }, h("label", {}, "Frames per second"), fps), h("div", { class: "row two" }, h("label", {}, "GIF width (px)"), gifW),
      h("div", { class: "btns" }, h("button", { class: "btn primary", onclick: () => exportVideo("webm", +secs.value, +fps.value) }, "WebM video"), h("button", { class: "btn", onclick: () => exportVideo("gif", +secs.value, Math.min(+fps.value, 25), +gifW.value) }, "GIF")), prog,
      h("p", { class: "hint" }, "WebM plays in browsers, PowerPoint 365 and Keynote; convert to MP4 with ffmpeg -i scene.webm scene.mp4.")));
    const snip = h("textarea", { class: "snippet", readonly: true });
    snip.value = DF.Export.embedSnippet(S.player.getScene());
    pane.append(sec("Live figure", null,
      h("div", { class: "btns" }, h("button", { class: "btn primary", onclick: exportHTML }, "Standalone page"), h("button", { class: "btn", onclick: copyEmbed }, "Copy embed code")), snip,
      h("p", { class: "hint" }, "The standalone page runs offline and contains everything; use it in reveal.js or as an iframe. The embed code needs dynflow.js from the dist folder next to your page.")));
    const rcode = h("textarea", { class: "snippet", readonly: true });
    rcode.value = "# In an R Markdown or Quarto document, or a pkgdown article\nhtmltools::tagList(\n  htmltools::tags$script(src = \"dynflow.js\"),\n  htmltools::HTML('<dyn-flow style=\"display:block;height:420px\" src=\"" + DF.Export.slug(S.scene.name) + ".json\" controls></dyn-flow>')\n)";
    pane.append(sec("R Markdown, Quarto, pkgdown", null, rcode, h("p", { class: "hint" }, "Save the scene file (JSON) next to the document and copy dist/dynflow.js into the site.")));
    pane.append(sec("Scene", null, h("div", { class: "btns" }, h("button", { class: "btn", onclick: exportJSON }, "Save scene (JSON)"), h("button", { class: "btn", onclick: openFile }, "Open scene"), h("button", { class: "btn", onclick: shareLink }, "Copy share link"))));
  }

  // ------------------------------------------------------------ exports
  const fileBase = () => DF.Export.slug(S.scene.name || S.model || "scene");
  async function exportPNG() { DF.Export.download(await DF.Export.png(S.player), fileBase() + ".png"); toast("PNG saved"); }
  function exportSVG() { DF.Export.download(DF.Export.svg(S.player), fileBase() + ".svg"); toast("SVG saved"); }
  function exportJSON() { DF.Export.download(new Blob([DF.Export.sceneJSON(S.player.getScene())], { type: "application/json" }), fileBase() + ".json"); toast("Scene saved"); }
  function exportHTML() {
    try { DF.Export.download(new Blob([DF.Export.standaloneHTML(S.player.getScene())], { type: "text/html" }), fileBase() + ".html"); toast("Standalone page saved"); }
    catch (e) { toast(e.message); }
  }
  async function copyEmbed() { await copy(DF.Export.embedSnippet(S.player.getScene())); toast("Embed code copied"); }
  async function shareLink() {
    const code = await DF.Export.encodeScene(S.player.getScene());
    const url = location.href.split("#")[0] + "#s=" + code;
    history.replaceState(null, "", "#s=" + code);
    await copy(url); toast("Share link copied (" + url.length + " characters)");
  }
  async function exportVideo(kind, seconds, fps, width) {
    const prog = $("#prog"), bar = prog ? $("i", prog) : null;
    if (prog) prog.classList.add("on");
    const onP = (f) => { if (bar) bar.style.width = Math.round(100 * f) + "%"; };
    toast(kind === "gif" ? "Capturing the GIF..." : "Recording the video...");
    try {
      const blob = kind === "gif" ? await DF.Export.gif(S.player, seconds, fps || 20, width || 640, onP) : await DF.Export.webm(S.player, seconds, fps || 30, onP);
      DF.Export.download(blob, fileBase() + "." + kind);
      toast((kind === "gif" ? "GIF" : "Video") + " saved, " + (blob.size / 1e6).toFixed(1) + " MB");
    } catch (e) { toast(e.message || String(e)); }
    if (prog) setTimeout(() => prog.classList.remove("on"), 600);
  }
  async function copy(text) { try { await navigator.clipboard.writeText(text); } catch (e) { const t = h("textarea"); t.value = text; document.body.append(t); t.select(); document.execCommand("copy"); t.remove(); } }
  function openFile() {
    const inp = h("input", { type: "file", accept: ".json,application/json", onchange: async (e) => {
      const f = e.target.files[0]; if (!f) return;
      try { const sc = JSON.parse(await f.text()); S.model = sc.model || null; setScene(sc); toast("Scene loaded"); } catch (err) { toast("Not a scene file: " + err.message); }
    } });
    inp.click();
  }
  function newModel() {
    S.model = null;
    setScene({ name: "My model", system: "# One equation per line; see Help for the language.\n# Holling type II predator and prey:\nN' = r*N*(1 - N/K) - a*N*P/(1 + a*h*N)\nP' = e*a*N*P/(1 + a*h*N) - m*P\nparam r = 1 [0, 2]\nparam K = 5 [1, 10]\nparam a = 1 [0, 2]\nparam h = 0.4 [0, 1]\nparam e = 0.6 [0, 1]\nparam m = 0.3 [0, 1]\ninit N = 3\ninit P = 1\nrange N = [0, 6]\nrange P = [0, 4]", view: { type: "phase" }, overlay: { title: "My model", equations: true } });
    S.tab = "model"; renderInspector();
    const b = $("#formulaBox"); if (b) { b.open = true; $("#formula").focus(); }
    document.body.classList.add("show-insp");
  }
  let toastTimer = 0;
  function toast(msg) { const t = $("#toast"); t.textContent = msg; t.classList.add("on"); clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.remove("on"), 2600); }

  // ------------------------------------------------------------- keys
  function keys() {
    document.addEventListener("keydown", (e) => {
      if (e.target.closest("input, textarea, select") || e.ctrlKey || e.metaKey || e.altKey) return;
      const k = e.key.toLowerCase();
      if (k === " ") { e.preventDefault(); S.player.toggle(); }
      else if (k === "r") restart();
      else if (k === "s") { S.player.pause(); S.player.step(); }
      else if (k === "n") newSeed();
      else if (k === "f") { const st = $("#stage"); if (document.fullscreenElement) document.exitFullscreen(); else st.requestFullscreen && st.requestFullscreen(); }
      else if (k === "p") exportPNG();
      else if (k === "escape") { $("#help").classList.remove("on"); document.body.classList.remove("show-lib", "show-insp"); }
    });
  }

  // ------------------------------------------------------------- start
  async function start() {
    buildTop(); buildLibrary(); buildStage(); buildInspector(); keys();
    $("#help").addEventListener("click", (e) => { if (e.target.id === "help" || e.target.closest("[data-close]")) $("#help").classList.remove("on"); });
    let scene = null;
    const hash = location.hash;
    if (hash.startsWith("#s=")) { try { scene = await DF.Export.decodeScene(hash.slice(3)); } catch (e) { toast("The link does not hold a valid scene"); } }
    else if (hash.startsWith("#model=")) { try { scene = DF.sceneFor(hash.slice(7)); } catch (e) { toast(e.message); } }
    if (!scene) { try { const saved = localStorage.getItem("dynflow:scene"); if (saved) scene = JSON.parse(saved); } catch (e) { scene = null; } }
    if (scene) { S.model = scene.model || null; setScene(scene); }
    else loadShowcase(0);
    setTimeout(() => { const hint = $("#hint"); if (hint) hint.style.opacity = 0; }, 7000);
  }
  document.readyState === "loading" ? document.addEventListener("DOMContentLoaded", start) : start();
  window.DynFlowStudio = S;
})();
