// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Headless-browser checks of the studio, the exports and the <relab-flow>
// element, through the Chrome DevTools protocol. Needs google-chrome (or
// CHROME=/path/to/chrome) and python3 for a static server.
// Run: node tests/browser.test.mjs
import { spawn } from "child_process";
import { writeFileSync, mkdirSync, rmSync } from "fs";
import { join } from "path";
import { tmpdir } from "os";
import { ROOT, reporter } from "./harness.mjs";

const R = reporter("browser");
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const PORT = 8790 + Math.floor(Math.random() * 100), DBG = 9700 + Math.floor(Math.random() * 200);
const OUT = join(tmpdir(), "relabflow-browser-test");
rmSync(OUT, { recursive: true, force: true }); mkdirSync(OUT, { recursive: true });

const server = spawn("python3", ["-m", "http.server", String(PORT), "--bind", "127.0.0.1"], { cwd: ROOT, stdio: "ignore" });
const chrome = spawn(process.env.CHROME || "google-chrome", ["--headless=new", "--disable-gpu", `--remote-debugging-port=${DBG}`, "--window-size=1500,950", `--user-data-dir=${join(OUT, "profile")}`, "--autoplay-policy=no-user-gesture-required", "about:blank"], { stdio: "ignore" });
const cleanup = () => { try { chrome.kill(); } catch (e) {} try { server.kill(); } catch (e) {} };
process.on("exit", cleanup);

async function page(url, init) {
  let target;
  for (let i = 0; i < 60 && !target; i++) { await sleep(200); try { target = await (await fetch(`http://127.0.0.1:${DBG}/json/new?about:blank`, { method: "PUT" })).json(); } catch (e) {} }
  const ws = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((r) => (ws.onopen = r));
  let id = 0; const pending = {}; const errors = [];
  ws.onmessage = (m) => {
    const d = JSON.parse(m.data);
    if (d.method === "Runtime.exceptionThrown") errors.push(d.params.exceptionDetails.exception ? d.params.exceptionDetails.exception.description : d.params.exceptionDetails.text);
    if (d.method === "Runtime.consoleAPICalled" && d.params.type === "error") errors.push(d.params.args.map((a) => a.value || a.description).join(" "));
    if (pending[d.id]) { pending[d.id](d); delete pending[d.id]; }
  };
  const send = (method, params = {}) => new Promise((r) => { pending[++id] = r; ws.send(JSON.stringify({ id, method, params })); });
  await send("Runtime.enable"); await send("Page.enable");
  await send("Emulation.setDeviceMetricsOverride", { width: 1500, height: 950, deviceScaleFactor: 1, mobile: false });
  if (init) await send("Page.addScriptToEvaluateOnNewDocument", { source: init });
  await send("Page.navigate", { url });
  await sleep(2500);
  const ev = async (expr) => {
    const r = await send("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.result.exceptionDetails) throw new Error(r.result.exceptionDetails.exception ? r.result.exceptionDetails.exception.description : r.result.exceptionDetails.text);
    return r.result.result.value;
  };
  return { ev, errors, close: () => ws.close(), send };
}

const base = `http://127.0.0.1:${PORT}`;
try {
  // --------------------------------------------------------------- studio
  const P = await page(base + "/index.html?test");
  await P.ev("localStorage.clear(); true");
  const models = await P.ev("RElabFlow.CATALOGUE.map(m => m.id)");
  const perModel = await P.ev(`(async () => {
    const S = window.RElabFlowStudio, out = [];
    for (const id of ${JSON.stringify(models)}) {
      document.querySelector('.item[data-id="' + id + '"]').click();
      const f0 = S.player.frameCount;
      await new Promise(r => setTimeout(r, 120));
      const err = document.querySelector('#errBar').classList.contains('on') ? document.querySelector('#errBar').textContent : '';
      const c = S.player.composite(), g = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
      let lit = 0; const bg = [g[0], g[1], g[2]];
      for (let i = 0; i < g.length; i += 4 * 37) if (Math.abs(g[i] - bg[0]) + Math.abs(g[i + 1] - bg[1]) + Math.abs(g[i + 2] - bg[2]) > 30) lit++;
      out.push({ id, err, frames: S.player.frameCount - f0, lit, view: S.player.scene.view.type });
    }
    return out;
  })()`);
  const bad = perModel.filter((m) => m.err || m.frames < 1 || m.lit < 3);
  R.check("studio: every catalogue model loads, animates and draws", bad.length === 0, `${perModel.length - bad.length} of ${perModel.length}${bad.length ? "; problems: " + bad.map((m) => m.id + "(" + (m.err || "frames " + m.frames + ", lit " + m.lit) + ")").join(", ") : ""}`);

  const viewRuns = await P.ev(`(async () => {
    const S = window.RElabFlowStudio, out = [];
    for (const id of ["lorenz", "rosenzweig-macarthur", "logistic", "double-well", "duffing", "cusp", "mackey-glass", "zaslavsky"]) {
      document.querySelector('.item[data-id="' + id + '"]').click();
      await new Promise(r => setTimeout(r, 60));
      const views = RElabFlow.viewsFor(S.player.sys);
      for (const v of views) {
        S.tab = 'view'; document.querySelector('[data-tab="view"]').click();
        const b = [...document.querySelectorAll('.vbtn')].find(x => x.textContent.trim() === { flow: 'Flow', trajectory: 'Trajectory', timeseries: 'Time series', phase: 'Phase plane', sweep: 'Sweep', orbit: 'Bifurcation', density: 'Density', strobe: 'Strobe', cobweb: 'Cobweb' }[v]);
        b.click();
        const f0 = S.player.frameCount;
        await new Promise(r => setTimeout(r, 150));
        const err = document.querySelector('#errBar').classList.contains('on') ? document.querySelector('#errBar').textContent : '';
        out.push({ id, v, ok: !err && S.player.frameCount > f0 && S.player.scene.view.type === v, err });
      }
    }
    return out;
  })()`);
  const vbad = viewRuns.filter((r) => !r.ok);
  R.check("studio: every applicable view runs", vbad.length === 0, `${viewRuns.length - vbad.length} of ${viewRuns.length} model-view pairs${vbad.length ? "; failed: " + vbad.map((r) => r.id + "/" + r.v + " " + r.err).join(", ") : ""}`);

  // Regression: switching model must not leave the previous figure on the layers.
  const leftover = await P.ev(`(async () => {
    const S = window.RElabFlowStudio, out = [];
    const lit = (c) => { const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data; let n = 0; for (let i = 3; i < d.length; i += 4 * 7) if (d[i] > 0) n++; return n; };
    for (const [a, b] of [["may-leonard", "lorenz"], ["thomas", "rosenzweig-macarthur"], ["henon", "cusp"], ["duffing", "sir"]]) {
      document.querySelector('.item[data-id="' + a + '"]').click(); await new Promise(r => setTimeout(r, 700));
      document.querySelector('.item[data-id="' + b + '"]').click();
      // the new view has drawn one frame at most; trail pixels left from the old scene show up here
      out.push(a + ' to ' + b + ': ' + lit(S.player.cv.trail));
    }
    return out;
  })()`);
  R.check("studio: switching model clears the previous figure", leftover.every((x) => / 0$/.test(x)), leftover.join("; "));

  // Tabs render without errors for a model with perturbations.
  const tabs = await P.ev(`(async () => {
    document.querySelector('.item[data-id="lorenz84-forced"]').click();
    await new Promise(r => setTimeout(r, 100));
    const out = [];
    for (const t of ['model', 'perturb', 'view', 'style', 'export']) { document.querySelector('[data-tab="' + t + '"]').click(); await new Promise(r => setTimeout(r, 40)); out.push(t + ':' + document.querySelector('#pane-' + t).children.length); }
    return out.join(' ');
  })()`);
  R.check("studio: inspector tabs render", !/:0/.test(tabs), tabs);

  // Formula editor: a new model compiles and runs; a broken one reports the line.
  const fe = await P.ev(`(async () => {
    const S = window.RElabFlowStudio;
    document.querySelector('[data-tab="model"]').click();
    const ta = document.querySelector('#formula');
    ta.value = "x' = y\\ny' = -x - c*y\\nparam c = 0.2 [0, 1]\\ninit x = 1";
    [...document.querySelectorAll('#pane-model .btn.primary')].find(b => b.textContent === 'Apply').click();
    await new Promise(r => setTimeout(r, 150));
    const good = S.player.sys.vars.join() + ' ' + S.player.frameCount;
    const ta2 = document.querySelector('#formula');
    ta2.value = "x' = y\\ny' = -x - c*\\nparam c = 0.2";
    [...document.querySelectorAll('#pane-model .btn.primary')].find(b => b.textContent === 'Apply').click();
    return { good, err: document.querySelector('#formulaErr').textContent };
  })()`);
  R.check("studio: formula editor compiles and reports errors", /^x,y \d+/.test(fe.good) && /line 2/.test(fe.err), `valid model runs (${fe.good}); broken model: "${fe.err}"`);

  // Parameter slider changes the running simulation live.
  const live = await P.ev(`(async () => {
    const S = window.RElabFlowStudio;
    document.querySelector('.item[data-id="lorenz"]').click(); await new Promise(r => setTimeout(r, 80));
    document.querySelector('[data-tab="model"]').click();
    const r = document.querySelector('input[data-param="rho"]'); r.value = 15; r.dispatchEvent(new Event('input'));
    return S.player.sim.p[S.player.sim.paramIndex('rho')];
  })()`);
  R.check("studio: parameter slider acts live", Math.abs(live - 15) < 1e-9, `rho in the simulator after moving the slider: ${live}`);

  // --------------------------------------------------------------- exports
  const ex = await P.ev(`(async () => {
    const S = window.RElabFlowStudio, E = RElabFlow.Export;
    document.querySelector('.item[data-id="lorenz"]').click(); await new Promise(r => setTimeout(r, 300));
    const png = await E.png(S.player);
    const pngHead = new Uint8Array(await png.slice(0, 8).arrayBuffer());
    const svg = S.player.toSVG();
    const gif = await E.gif(S.player, 0.5, 10, 320);
    const gifBytes = Array.from(new Uint8Array(await gif.arrayBuffer()));
    let webm = null; try { const w = await E.webm(S.player, 1, 20); webm = Array.from(new Uint8Array(await w.arrayBuffer())); } catch (e) { webm = String(e.message || e); }
    const scene = S.player.getScene();
    const code = await E.encodeScene(scene), back = await E.decodeScene(code);
    const html = E.standaloneHTML(scene);
    return { pngHead: Array.from(pngHead), pngSize: png.size, svgPaths: (svg.match(/<path/g) || []).length, svgOk: svg.startsWith('<svg'), gifBytes, webm, round: JSON.stringify(back) === JSON.stringify(scene), codeLen: code.length, jsonLen: JSON.stringify(scene).length, html };
  })()`);
  R.check("export: PNG", ex.pngHead.join(",") === "137,80,78,71,13,10,26,10" && ex.pngSize > 10000, `PNG signature present, ${(ex.pngSize / 1024).toFixed(0)} KiB`);
  R.check("export: SVG keeps the orbit as vector paths", ex.svgOk && ex.svgPaths >= 1, `${ex.svgPaths} path elements`);
  writeFileSync(join(OUT, "t.gif"), Buffer.from(ex.gifBytes));
  if (Array.isArray(ex.webm)) writeFileSync(join(OUT, "t.webm"), Buffer.from(ex.webm));
  const round = ex.round;
  R.check("export: share link round trip", round, `scene of ${ex.jsonLen} characters encoded in ${ex.codeLen}`);
  writeFileSync(join(OUT, "standalone.html"), ex.html);

  // GIF and WebM decoded by independent tools.
  const run = (cmd, args) => new Promise((res) => { const p = spawn(cmd, args); let o = ""; p.stdout.on("data", (d) => (o += d)); p.stderr.on("data", (d) => (o += d)); p.on("close", (c) => res({ code: c, out: o })); p.on("error", () => res({ code: -1, out: "" })); });
  const gifInfo = await run("python3", ["-c", `from PIL import Image, ImageSequence; im = Image.open(${JSON.stringify(join(OUT, "t.gif"))}); print(sum(1 for _ in ImageSequence.Iterator(im)), im.size[0], im.size[1])`]);
  R.check("export: GIF decodes", gifInfo.code === 0 && +gifInfo.out.split(" ")[0] === 5, `Pillow reads ${gifInfo.out.trim()} (frames, width, height)`);
  if (Array.isArray(ex.webm)) {
    const probe = await run("ffprobe", ["-v", "error", "-show_entries", "stream=codec_name,width,height", "-of", "csv=p=0", join(OUT, "t.webm")]);
    R.check("export: WebM decodes", probe.code === 0 && /vp[89]/.test(probe.out), `ffprobe: ${probe.out.trim() || "(no output)"}`);
  } else R.check("export: WebM", false, "recording failed: " + ex.webm);

  // ------------------------------------------------ defects found by the audit
  const click = (id) => `document.querySelector('.item[data-id="${id}"]').click(); await new Promise(r => setTimeout(r, 80));`;
  const fx = await P.ev(`(async () => {
    const S = window.RElabFlowStudio, out = {};
    // Time series: the buffer covers the window (seasonal SIR covered 10 percent of it).
    ${click("seasonal-sir")} let p = S.player; p.pause(); while (p.sim.t < 25) p.tick();
    let v = p.view; out.cover = (p.sim.t - v.T[(v.head - v.len + v.cap) % v.cap]) / v.span;
    // Density carpet: nx columns of window / nx time units (the axis said 200 over 14.7 units).
    ${click("double-well")} p = S.player; v = p.view; out.carpet = [v.span, v.colDt * v.nx];
    // Strobe: T = N h exactly, samples at t = k T (34 of 60 sections were one step late).
    ${click("seasonal-rm")} p = S.player; p.pause(); p.advance(30); v = p.view; out.strobe = [Math.abs(v.N * p.sim.h - v.period), Math.abs(p.sim.t - p.sim.steps * p.sim.h)];
    // Refresh rate: 1 s of playback in 144 or 60 frames covers the same model time.
    ${click("lorenz")} p = S.player; p.pause(); let t0 = p.sim.t; p.advance(144, 1000 / 144); const a = p.sim.t - t0; t0 = p.sim.t; p.advance(60, 1000 / 60); out.refresh = [a, p.sim.t - t0, p.rate, p.sim.h];
    // WebM: one video second is one second of playback (it was half a second).
    t0 = p.sim.t; await RElabFlow.Export.webm(p, 1, 30); out.webm = [p.sim.t - t0, p.rate]; p.pause();
    // Mackey-Glass in tau: the bifurcation diagram keeps the delay (it was empty).
    ${click("mackey-glass")} S.scene.view = { type: 'orbit', param: 'tau', from: 2, to: 30 }; S.player.load(JSON.parse(JSON.stringify(S.scene))); p = S.player; p.pause();
    for (let i = 0; i < 4000 && p.view.col <= p.view.cols; i++) p.tick();
    const pts = p.view.points, at17 = [], early = [];
    for (let i = 0; i < pts.length; i += 2) { if (Math.abs(pts[i] - 17) < 0.05) at17.push(pts[i + 1]); if (pts[i] < 4.6 && Math.abs(pts[i + 1] - 1) > 1e-3) early.push(pts[i + 1]); }
    out.mackey = [pts.length / 2, Math.min(...at17), Math.max(...at17), early.length];
    // Phase plane of a delay equation: the orbit is the simulator's solution (it ignored the delay).
    ${click("delayed-predator-prey")} S.scene.view = { type: 'phase', seeds: 1 }; S.scene.n = 1; S.player.load(JSON.parse(JSON.stringify(S.scene))); p = S.player; p.pause();
    let diff = 0; for (let f = 0; f < 400; f++) { p.tick(); const x = p.view.trails[0].sim.X; diff = Math.max(diff, Math.abs(x[0] - p.sim.X[0]), Math.abs(x[1] - p.sim.X[1])); } out.phaseDDE = diff;
    // Phase plane of a forced oscillator: the orbit keeps its own time (it was frozen within a frame).
    ${click("driven-oscillator")} S.scene.view = { type: 'phase', seeds: 1 }; S.scene.n = 1; S.player.load(JSON.parse(JSON.stringify(S.scene))); p = S.player; p.pause();
    while (p.sim.t < 150) p.tick(); let amp = 0; const t1 = p.sim.t; while (p.sim.t < t1 + 30) { p.tick(); amp = Math.max(amp, Math.abs(p.view.trails[0].sim.X[0])); } out.forced = amp;
    // Sweep: continued branches with their folds, recomputed when another parameter moves (they went stale).
    ${click("cusp")} p = S.player; v = p.view; const f1 = v.data.points.filter(q => q.kind === 'fold').map(q => q.p).sort((x, y) => x - y);
    p.pause(); p.setParam('a', 0.5); const f2 = v.data.points.filter(q => q.kind === 'fold').map(q => q.p).sort((x, y) => x - y);
    out.sweep = [v.data.runs.length, f1, f2];
    // Nullclines: joined polylines, so the dash pattern shows (each tiny segment restarted it).
    ${click("rosenzweig-macarthur")} p = S.player; out.nullLines = p.view.nullLines.length;
    // PNG: legend and typeset equations are painted (the legend was missing, the equations were source text).
    p.pause(); S.scene.overlay.equations = true; p.setOverlay({ equations: true });
    const px = (legend) => { p.setOverlay({ legend: legend }); const c = p.composite(), L = p.legendLayout((t) => 6 * t.length), g = c.getContext('2d'), d = L ? g.getImageData(Math.round(L.x * p.dpr), Math.round(L.y * p.dpr), Math.round(L.w * p.dpr), Math.round(L.h * p.dpr)).data : null; return { L: !!L, sum: d ? d.reduce((s, q) => s + q, 0) : 0 }; };
    const withL = px(true), without = px(false); p.setOverlay({ legend: true });
    out.png = [withL.L, withL.sum !== without.sum, RElabFlow.MathType.systemBoxes(p.scene.system, 15).map(b => b.items.filter(i => i.k === 'r').length)];
    // SVG: vector figure with axes, no embedded raster (a time series lost its axes).
    ${click("sir")} await new Promise(r => setTimeout(r, 300)); const svg = S.player.toSVG(); out.svg = [(svg.match(/<image/g) || []).length, (svg.match(/<text/g) || []).length];
    return out;
  })()`);
  R.check("audit: time series fill their window", fx.cover > 0.99, `seasonal SIR: the buffer covers ${fx.cover.toFixed(2)} of the window`);
  R.check("audit: density carpet axis", Math.abs(fx.carpet[1] - fx.carpet[0]) < 1e-9, `axis ${fx.carpet[0]} time units; columns cover ${fx.carpet[1].toFixed(6)}`);
  R.check("audit: stroboscopic samples on t = kT", fx.strobe[0] < 1e-12 && fx.strobe[1] < 1e-9, `|N h - T| = ${fx.strobe[0].toExponential(1)}, |t - steps h| = ${fx.strobe[1].toExponential(1)}`);
  R.check("audit: speed independent of the refresh rate", Math.abs(fx.refresh[0] - fx.refresh[1]) <= fx.refresh[3] + 1e-9 && Math.abs(fx.refresh[1] - fx.refresh[2]) <= fx.refresh[3] + 1e-9, `1 s of playback: ${fx.refresh[0].toFixed(4)} time units in 144 frames, ${fx.refresh[1].toFixed(4)} in 60 (rate ${fx.refresh[2]}, step ${fx.refresh[3]})`);
  R.check("audit: WebM plays at the live speed", Math.abs(fx.webm[0] - fx.webm[1]) / fx.webm[1] < 0.02, `one video second holds ${fx.webm[0].toFixed(4)} time units; live rate ${fx.webm[1]} per second`);
  R.check("audit: bifurcation diagram of a delay equation", fx.mackey[0] > 1000 && fx.mackey[1] > 0.95 && fx.mackey[1] < 1 && fx.mackey[2] > 1.3 && fx.mackey[2] < 1.34 && fx.mackey[3] === 0, `Mackey-Glass in tau: ${fx.mackey[0]} points; peaks at tau = 17 in [${fx.mackey[1].toFixed(3)}, ${fx.mackey[2].toFixed(3)}] (a direct run gives [0.971, 1.317]); none below the Hopf point tau = 4.71`);
  R.check("audit: phase-plane orbits of delay and forced equations", fx.phaseDDE === 0 && Math.abs(fx.forced - 0.5 / Math.hypot(1 - 1.69, 0.26)) < 5e-4, `delayed predator and prey: the orbit equals the simulator's solution (max difference ${fx.phaseDDE}); driven oscillator: amplitude ${fx.forced.toFixed(4)} (exact ${(0.5 / Math.hypot(1 - 1.69, 0.26)).toFixed(4)})`);
  const e1 = 2 * Math.pow(1 / 3, 1.5), e2 = 2 * Math.pow(0.5 / 3, 1.5);
  R.check("audit: sweep branches as curves with their folds", fx.sweep[0] >= 3 && Math.abs(fx.sweep[1][0] + e1) < 1e-3 && Math.abs(fx.sweep[1][1] - e1) < 1e-3 && Math.abs(fx.sweep[2][1] - e2) < 1e-3, `${fx.sweep[0]} runs of equal stability; folds at r = ${fx.sweep[1].map((q) => q.toFixed(4)).join(", ")} (exact -/+${e1.toFixed(4)}); after a = 0.5, ${fx.sweep[2].map((q) => q.toFixed(4)).join(", ")} (exact -/+${e2.toFixed(4)})`);
  R.check("audit: nullclines are joined", fx.nullLines <= 6, `${fx.nullLines} polylines for the two nullclines of Rosenzweig-MacArthur`);
  R.check("audit: PNG carries the legend and typeset equations", fx.png[0] && fx.png[1] && fx.png[2].join() === "2,1", `legend painted: ${fx.png[1]}; fraction bars per equation: ${fx.png[2].join(", ")} (N' has N/K and aNP/(1 + ahN), P' has eaNP/(1 + ahN))`);
  R.check("audit: SVG of a time series is vector with axes", fx.svg[0] === 0 && fx.svg[1] >= 8, `${fx.svg[0]} embedded images, ${fx.svg[1]} text elements`);

  const leaked = P.errors.filter((e) => !/favicon/.test(e));
  R.check("studio: no uncaught errors", leaked.length === 0, leaked.length ? leaked.slice(0, 3).join(" | ") : "console clean");
  P.close();

  // ------------------------------------------- opening the studio (audit)
  const LT = "window.__long = []; new PerformanceObserver((l) => l.getEntries().forEach((e) => window.__long.push(e.duration))).observe({ type: 'longtask', buffered: true });";
  const O = await page(base + "/index.html", LT);
  await sleep(8000);
  const lt = await O.ev("window.__long");
  R.check("audit: no freeze while the studio opens", Math.max(0, ...lt) < 400, `longest task in the first 10 s: ${Math.round(Math.max(0, ...lt))} ms (a preview froze the page for 2 s)`);
  const ui = await O.ev(`({ logo: document.querySelector('.brand').href, target: document.querySelector('.brand').target, search: document.querySelector('.lib-head input').placeholder, menu: [...document.querySelectorAll('#exportMenu button small')].map(x => x.textContent).join(''), n: RElabFlow.CATALOGUE.length })`);
  R.check("audit: logo link and interface text", ui.logo === "https://robustecologies.github.io/" && ui.target === "_blank" && ui.search === "Search " + ui.n + " models" && ui.menu === "P", `logo opens ${ui.logo} in a new tab; search box "${ui.search}"; keys in the export menu: ${ui.menu}`);
  O.close();
  const brokenLink = "j" + Buffer.from(JSON.stringify({ name: "Broken", system: "x' = y\ny' = -x - c*\nparam c = 0.2" })).toString("base64").replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
  const L1 = await page(base + "/index.html#s=" + brokenLink);
  await sleep(800);
  const bl = await L1.ev(`({ scene: RElabFlowStudio.scene && RElabFlowStudio.scene.name, toast: document.querySelector('#toast').textContent, saved: (JSON.parse(localStorage.getItem('relabflow:scene') || '{}')).name })`);
  L1.close();
  const L2 = await page(base + "/index.html");
  const re = await L2.ev(`({ scene: RElabFlowStudio.scene && RElabFlowStudio.scene.name, inspector: document.querySelector('#pane-model').children.length })`);
  await L2.ev(`(localStorage.setItem('relabflow:scene', JSON.stringify({ name: 'Bad', system: "x' = (" })), 1)`);
  L2.close();
  const L3 = await page(base + "/index.html");
  const rb = await L3.ev(`({ scene: RElabFlowStudio.scene && RElabFlowStudio.scene.name, toast: document.querySelector('#toast').textContent })`);
  await L3.ev("localStorage.clear(), 1");
  L3.close();
  R.check("audit: a scene that does not compile never blocks the studio", bl.scene && /does not compile/.test(bl.toast) && bl.saved !== "Broken" && re.scene && re.inspector > 0 && rb.scene && /could not be opened/.test(rb.toast),
    `broken link: "${bl.toast}", showing ${bl.scene}, saved ${bl.saved}; next start: ${re.scene}; a broken saved scene: "${rb.toast}", showing ${rb.scene}`);

  // ------------------------------------------------------ standalone page
  writeFileSync(join(ROOT, "examples", ".standalone-test.html"), ex.html);
  const Q = await page(base + "/examples/.standalone-test.html");
  await sleep(1500);
  const st = await Q.ev(`(() => { const el = document.querySelector('relab-flow'); return { frames: el.player ? el.player.frameCount : -1, type: el.player ? el.player.scene.view.type : null, net: performance.getEntriesByType('resource').map(e => e.name) }; })()`);
  R.check("standalone page runs with no other file", st.frames > 20 && st.net.length === 0, `${st.frames} frames of a ${st.type} view; network requests: ${st.net.length ? st.net.join(", ") : "none"}`);
  Q.close();
  rmSync(join(ROOT, "examples", ".standalone-test.html"), { force: true });

  // ------------------------------------------------------ <relab-flow> element
  const C = await page(base + "/examples/embed.html");
  await sleep(1500);
  const els = await C.ev(`(() => [...document.querySelectorAll('relab-flow')].map(e => ({ how: e.getAttribute('data-how'), frames: e.player ? e.player.frameCount : -1, err: e.shadowRoot && e.shadowRoot.querySelector('.err') ? e.shadowRoot.querySelector('.err').textContent : '' })))()`);
  const ebad = els.filter((e) => e.frames < 1 && e.how !== "error" && e.how !== "paused");
  const errEl = els.find((e) => e.how === "error"), pausedEl = els.find((e) => e.how === "paused");
  R.check("element: scene, src and model attributes", ebad.length === 0 && els.length >= 4, els.map((e) => e.how + " " + e.frames).join(", "));
  R.check("element: paused attribute and error message", pausedEl && pausedEl.frames === 0 && errEl && /Unknown name/.test(errEl.err), `paused element at ${pausedEl && pausedEl.frames} frames; broken scene shows "${errEl && errEl.err}"`);
  const mv = await C.ev(`(async () => { const el = document.querySelector('relab-flow'), p = el.player; const div = document.createElement('div'); document.body.prepend(div); div.append(el); el.scrollIntoView(); const f0 = el.player.frameCount; await new Promise(r => setTimeout(r, 500)); return { same: el.player === p, frames: el.player.frameCount - f0 }; })()`);
  R.check("element: moving it to another parent keeps the figure", mv.same && mv.frames > 5 && C.errors.length === 0, `same player after the move: ${mv.same}; ${mv.frames} frames in 0.5 s; errors: ${C.errors.length ? C.errors[0] : "none"}`);
  C.close();
} catch (e) {
  R.check("browser run", false, e.stack || String(e));
}
cleanup();
process.exit(R.done().fails ? 1 : 0);
