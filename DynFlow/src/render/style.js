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
