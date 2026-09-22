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
