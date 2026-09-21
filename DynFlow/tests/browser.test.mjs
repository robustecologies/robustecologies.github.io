// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Headless-browser checks of the studio, the exports and the <dyn-flow>
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
const OUT = join(tmpdir(), "dynflow-browser-test");
rmSync(OUT, { recursive: true, force: true }); mkdirSync(OUT, { recursive: true });

const server = spawn("python3", ["-m", "http.server", String(PORT), "--bind", "127.0.0.1"], { cwd: ROOT, stdio: "ignore" });
const chrome = spawn(process.env.CHROME || "google-chrome", ["--headless=new", "--disable-gpu", `--remote-debugging-port=${DBG}`, "--window-size=1500,950", `--user-data-dir=${join(OUT, "profile")}`, "--autoplay-policy=no-user-gesture-required", "about:blank"], { stdio: "ignore" });
const cleanup = () => { try { chrome.kill(); } catch (e) {} try { server.kill(); } catch (e) {} };
process.on("exit", cleanup);

async function page(url) {
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
  const models = await P.ev("DynFlow.CATALOGUE.map(m => m.id)");
  const perModel = await P.ev(`(async () => {
    const S = window.DynFlowStudio, out = [];
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
    const S = window.DynFlowStudio, out = [];
    for (const id of ["lorenz", "rosenzweig-macarthur", "logistic", "double-well", "duffing", "cusp", "mackey-glass", "zaslavsky"]) {
      document.querySelector('.item[data-id="' + id + '"]').click();
      await new Promise(r => setTimeout(r, 60));
      const views = DynFlow.viewsFor(S.player.sys);
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
    const S = window.DynFlowStudio, out = [];
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
    const S = window.DynFlowStudio;
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
    const S = window.DynFlowStudio;
    document.querySelector('.item[data-id="lorenz"]').click(); await new Promise(r => setTimeout(r, 80));
    document.querySelector('[data-tab="model"]').click();
    const r = document.querySelector('input[data-param="rho"]'); r.value = 15; r.dispatchEvent(new Event('input'));
    return S.player.sim.p[S.player.sim.paramIndex('rho')];
  })()`);
  R.check("studio: parameter slider acts live", Math.abs(live - 15) < 1e-9, `rho in the simulator after moving the slider: ${live}`);

  // --------------------------------------------------------------- exports
  const ex = await P.ev(`(async () => {
    const S = window.DynFlowStudio, E = DynFlow.Export;
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

  const leaked = P.errors.filter((e) => !/favicon/.test(e));
  R.check("studio: no uncaught errors", leaked.length === 0, leaked.length ? leaked.slice(0, 3).join(" | ") : "console clean");
  P.close();

  // ------------------------------------------------------ standalone page
  writeFileSync(join(ROOT, "examples", ".standalone-test.html"), ex.html);
  const Q = await page(base + "/examples/.standalone-test.html");
  await sleep(1500);
  const st = await Q.ev(`(() => { const el = document.querySelector('dyn-flow'); return { frames: el.player ? el.player.frameCount : -1, type: el.player ? el.player.scene.view.type : null, net: performance.getEntriesByType('resource').map(e => e.name) }; })()`);
  R.check("standalone page runs with no other file", st.frames > 20 && st.net.length === 0, `${st.frames} frames of a ${st.type} view; network requests: ${st.net.length ? st.net.join(", ") : "none"}`);
  Q.close();
  rmSync(join(ROOT, "examples", ".standalone-test.html"), { force: true });

  // ------------------------------------------------------ <dyn-flow> element
  const C = await page(base + "/examples/embed.html");
  await sleep(1500);
  const els = await C.ev(`(() => [...document.querySelectorAll('dyn-flow')].map(e => ({ how: e.getAttribute('data-how'), frames: e.player ? e.player.frameCount : -1, err: e.shadowRoot && e.shadowRoot.querySelector('.err') ? e.shadowRoot.querySelector('.err').textContent : '' })))()`);
  const ebad = els.filter((e) => e.frames < 1 && e.how !== "error" && e.how !== "paused");
  const errEl = els.find((e) => e.how === "error"), pausedEl = els.find((e) => e.how === "paused");
  R.check("element: scene, src and model attributes", ebad.length === 0 && els.length >= 4, els.map((e) => e.how + " " + e.frames).join(", "));
  R.check("element: paused attribute and error message", pausedEl && pausedEl.frames === 0 && errEl && /Unknown name/.test(errEl.err), `paused element at ${pausedEl && pausedEl.frames} frames; broken scene shows "${errEl && errEl.err}"`);
  C.close();
} catch (e) {
  R.check("browser run", false, e.stack || String(e));
}
cleanup();
process.exit(R.done().fails ? 1 : 0);
