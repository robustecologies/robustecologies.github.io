/* Interface of the knowledge base site: routing, views, graph, search and recall.
   It reads window.KB_DATA, which tools/kb.py build writes into data/kb-data.js. */
(function () {
  "use strict";

  const DATA = window.KB_DATA;
  const main = document.getElementById("main");
  if (!DATA) {
    main.innerHTML = '<div class="empty"><h1>No data</h1><p>The file data/kb-data.js is missing. Run <code>python3 tools/kb.py build</code> and open the built site.</p></div>';
    return;
  }

  // ============================ Utilities
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));
  const esc = v => String(v == null ? "" : v).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

  // Highlighting for the check scripts, which are R. The text is escaped first and the spans
  // are added to the escaped text, so no markup from a script can reach the page.
  const CODE_KEYWORDS = /\b(function|for|while|repeat|if|else|in|next|break|return|TRUE|FALSE|NULL|NA|NaN|Inf|library|require|def|import|from|lambda|None|True|False)\b/g;
  function highlightCode(src) {
    return esc(src).split("\n").map(line => {
      const hash = line.indexOf("#");
      let code = line, tail = "";
      if (hash >= 0) {
        const before = line.slice(0, hash);
        if ((before.match(/&quot;/g) || []).length % 2 === 0) {
          code = before;
          tail = `<span class="c-com">${line.slice(hash)}</span>`;
        }
      }
      code = code
        .replace(/(&quot;[^&]*?&quot;|&#39;[^&]*?&#39;)/g, '<span class="c-str">$1</span>')
        .replace(/\b(\d+\.?\d*([eE][-+]?\d+)?)\b/g, '<span class="c-num">$1</span>')
        .replace(CODE_KEYWORDS, '<span class="c-kw">$1</span>')
        .replace(/([A-Za-z._][A-Za-z0-9._]*)(\()/g, '<span class="c-fn">$1</span>$2');
      return code + tail;
    }).join("\n");
  }

  const enc = encodeURIComponent;
  const debounce = (fn, ms) => { let t; return (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); }; };
  const stripTags = html => { const d = document.createElement("div"); d.innerHTML = html || ""; return d.textContent || ""; };
  const plainTitle = t => String(t || "").replace(/\$([^$]*)\$/g, (m, tex) => tex
    .replace(/\\(?:mathrm|mathbf|mathcal|mathbb|operatorname|text)\{([^}]*)\}/g, "$1")
    .replace(/\\([A-Za-z]+)/g, "$1").replace(/[{}]/g, ""));
  const store = {
    get(key, fallback) { try { const v = localStorage.getItem("kb:" + key); return v == null ? fallback : JSON.parse(v); } catch (e) { return fallback; } },
    set(key, value) { try { localStorage.setItem("kb:" + key, JSON.stringify(value)); } catch (e) { /* storage unavailable */ } },
  };
  const cssVar = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  function isDark() {
    const t = document.documentElement.getAttribute("data-theme");
    return t ? t === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  const isTyping = el => el && (el.tagName === "INPUT" || el.tagName === "TEXTAREA" || el.tagName === "SELECT" || el.isContentEditable);

  // ============================ Data
  const META = DATA.meta || {};
  const S = DATA.schema || {};
  const A = DATA.analytics || {};
  const NOTES = DATA.notes || [];
  const EDGES = DATA.edges || [];
  const REL = S.relations || {};
  const LEVELS = S.levels || {};
  const MACROS = S.macros || {};
  const PHILO = S.philosophy || { dimensions: {}, themes: {} };
  const SYMBOLS = DATA.symbols || [];
  const symbolByKey = new Map(SYMBOLS.map(s => [s.key, s]));
  const byId = new Map(NOTES.map(n => [n.id, n]));
  const outE = new Map();
  const inE = new Map();
  for (const e of EDGES) {
    if (!outE.has(e.s)) outE.set(e.s, []);
    outE.get(e.s).push(e);
    if (!inE.has(e.t)) inE.set(e.t, []);
    inE.get(e.t).push(e);
  }
  const CONCEPT_TYPES = ["concept", "theorem", "method", "model", "example"];
  const isTyped = e => Object.prototype.hasOwnProperty.call(REL, e.k);
  const LOCAL = META.mode !== "public";

  const typeLabel = t => ((S.types || {})[t] || {}).label || t;
  const typePlural = t => ((S.types || {})[t] || {}).plural || t;
  const level = (kind, id) => ((LEVELS[kind] || {})[id]) || {};
  const levelLabel = (kind, id) => level(kind, id).label || id;
  // Level ids in increasing rank, whatever the key order of the data file.
  const levelIds = kind => Object.keys(LEVELS[kind] || {}).sort((a, b) => (level(kind, a).rank || 0) - (level(kind, b).rank || 0));
  const levelIdsDesc = kind => levelIds(kind).reverse();
  const areaIds = () => Object.keys(S.areas || {}).sort((a, b) => (S.areas[a].rank || 0) - (S.areas[b].rank || 0));
  const areaColor = area => { const a = (S.areas || {})[area]; return a ? (isDark() ? a.dark : a.light) : cssVar("--neutral-node"); };
  const domainLabel = d => ((S.domains || {})[d] || {}).label || d;
  const noteTitle = id => { const n = byId.get(id); return n ? n.title_html : esc(id); };

  // ============================ Math
  function typeset(root) {
    if (!window.katex || !root) return;
    root.querySelectorAll(".math:not([data-typeset])").forEach(el => {
      const tex = el.textContent;
      try {
        window.katex.render(tex, el, { displayMode: el.classList.contains("math-display"), throwOnError: false, strict: "ignore", macros: Object.assign({}, MACROS) });
      } catch (err) {
        el.classList.add("math-error");
        el.title = String(err && err.message || err);
      }
      el.setAttribute("data-typeset", "1");
    });
  }

  // ============================ Small components
  const SHAPES = {
    dot: '<circle cx="8" cy="8" r="5.5"/>',
    diamond: '<path d="M8 1.5 14.5 8 8 14.5 1.5 8Z"/>',
    square: '<rect x="2.5" y="2.5" width="11" height="11" rx="1"/>',
    triangle: '<path d="M8 1.8 14.6 13.6H1.4Z"/>',
    star: '<path d="m8 1.4 2 4.2 4.6.6-3.4 3.2.9 4.6L8 11.8 3.9 14l.9-4.6L1.4 6.2 6 5.6Z"/>',
    hexagon: '<path d="M4.4 2h7.2L15 8l-3.4 6H4.4L1 8Z"/>',
    triangleDown: '<path d="M1.4 2.4h13.2L8 14.2Z"/>',
    box: '<rect x="1.5" y="3.5" width="13" height="9" rx="2"/>',
  };
  const shapeOf = type => ((S.types || {})[type] || {}).shape || "dot";
  const shapeIcon = (type, color) => `<svg class="shape" viewBox="0 0 16 16" aria-hidden="true" style="fill:${color || "currentColor"}">${SHAPES[shapeOf(type)] || SHAPES.dot}</svg>`;
  const pips = (k, n) => `<span class="pips" aria-hidden="true">${Array.from({ length: n }, (_, i) => `<i class="pip${i < k ? " on" : ""}"></i>`).join("")}</span>`;

  function typeChip(n) { return `<span class="chip">${shapeIcon(n.type)}${esc(typeLabel(n.type))}</span>`; }
  function maturityChip(n) {
    if (!n.maturity) return "";
    const lv = level("maturity", n.maturity);
    return `<span class="chip" data-tip="${esc(lv.definition || "")}">${pips(lv.rank || 0, 3)}${esc(lv.label || n.maturity)}</span>`;
  }
  function confidenceChip(n) {
    if (!n.confidence) return "";
    const lv = level("confidence", n.confidence);
    return `<span class="chip chip-conf conf-${esc(n.confidence)}" data-tip="${esc(lv.definition || "")}">${esc(lv.label || n.confidence)} confidence</span>`;
  }
  function statusChip(n) {
    if (n.type !== "question" || !n.status) return "";
    const icon = { open: "○", "in-progress": "◐", resolved: "●" }[n.status] || "";
    const cls = { open: "chip-warn", "in-progress": "", resolved: "chip-good" }[n.status] || "";
    return `<span class="chip ${cls}" data-tip="${esc(level("status", n.status).definition || "")}">${icon} ${esc(levelLabel("status", n.status))}</span>`;
  }
  function verificationChip(n) {
    const v = n.verification || {};
    if (!v.status) return "";
    const icon = { verified: "✓", manual: "✓", unverified: "!", mismatch: "✕" }[v.status] || "";
    const cls = { verified: "chip-good", manual: "chip-good", unverified: "chip-warn", mismatch: "chip-bad" }[v.status] || "";
    return `<span class="chip ${cls}">${icon} ${esc(levelLabel("verification", v.status))}</span>`;
  }
  function readingChip(n) {
    if (n.type !== "source" || !n.reading) return "";
    return `<span class="chip" data-tip="${esc(level("reading", n.reading).definition || "")}">${pips(level("reading", n.reading).rank || 0, 4)}${esc(levelLabel("reading", n.reading))}</span>`;
  }
  function understandingChip(n) {
    if (!n.understanding) return "";
    return `<span class="chip" data-tip="The owner's own rating of understanding, from 1 to 5">${pips(n.understanding, 5)}Understanding</span>`;
  }
  function levelChips(n) {
    return [maturityChip(n), confidenceChip(n), statusChip(n), verificationChip(n), readingChip(n), understandingChip(n),
      n.visibility === "private" ? '<span class="chip chip-warn">Private</span>' : ""].join("");
  }
  // The marks follow the logo: the existential quantifier for what there is, the
  // turnstile for what can be derived and known.
  const DIM_MARK = { ontology: "∃", epistemology: "⊢" };
  const themeOf = id => (PHILO.themes || {})[id] || null;
  function themeChip(id) {
    const t = themeOf(id);
    if (!t) return "";
    return `<a class="chip chip-theme dim-${esc(t.dimension)}" href="#/philosophy/${enc(id)}" data-tip="${esc(t.question || "")}"><span class="dim-mark">${DIM_MARK[t.dimension] || ""}</span>${esc(t.label)}</a>`;
  }
  const kpiHtml = (v, k) => `<div class="kpi"><div class="v">${v}</div><div class="k">${k}</div></div>`;
  const noteLink = (id, html) => byId.has(id)
    ? `<a class="wikilink" href="#/n/${enc(id)}" data-id="${esc(id)}">${html || noteTitle(id)}</a>`
    : `<span class="wikilink hidden">${esc(id)}</span>`;
  const datalist = () => `<datalist id="note-titles">${NOTES.filter(n => n.type !== "source").map(n => `<option value="${esc(plainTitle(n.title))}"></option>`).join("")}</datalist>`;

  function findNote(text) {
    const q = String(text || "").trim();
    if (!q) return null;
    if (byId.has(q)) return q;
    const low = q.toLowerCase();
    for (const n of NOTES) {
      if (plainTitle(n.title).toLowerCase() === low || (n.aliases || []).some(a => a.toLowerCase() === low)) return n.id;
    }
    const hits = searchNotes(q);
    return hits.length ? hits[0].id : null;
  }

  function relationGroups(id) {
    const groups = new Map();
    const add = (label, other, note) => { if (!groups.has(label)) groups.set(label, []); groups.get(label).push({ id: other, note }); };
    for (const e of outE.get(id) || []) if (isTyped(e)) add(REL[e.k].label, e.t, e.n);
    for (const e of inE.get(id) || []) if (isTyped(e)) add(REL[e.k].symmetric ? REL[e.k].label : REL[e.k].inverse, e.s, e.n);
    return groups;
  }

  // ============================ Popover and tooltip
  const popover = $("#popover");
  const tip = $("#tip");
  let popTimer = null;
  function previewHtml(id) {
    const n = byId.get(id);
    if (!n) return "";
    if (n.type === "source") {
      return `<div class="chips">${typeChip(n)}${verificationChip(n)}</div><h4>${esc(n.label || plainTitle(n.title))}</h4><p class="reference">${n.reference_html || ""}</p>`;
    }
    return `<div class="chips">${typeChip(n)}${maturityChip(n)}</div><h4>${n.title_html}</h4><p>${n.summary_html || ""}</p>`;
  }
  function showPopover(anchor, id) {
    const html = previewHtml(id);
    if (!html) return;
    popover.innerHTML = html;
    popover.hidden = false;
    typeset(popover);
    const r = anchor.getBoundingClientRect();
    const w = popover.offsetWidth;
    const h = popover.offsetHeight;
    let top = r.bottom + 8;
    if (top + h > window.innerHeight - 8) top = Math.max(8, r.top - h - 8);
    const left = Math.min(Math.max(8, r.left), window.innerWidth - w - 8);
    popover.style.top = top + "px";
    popover.style.left = left + "px";
  }
  function hidePopover() { clearTimeout(popTimer); popover.hidden = true; }
  function hideTip() { tip.hidden = true; }
  const canHover = window.matchMedia("(hover: hover)").matches;
  document.addEventListener("mouseover", ev => {
    const t = ev.target;
    if (!(t instanceof Element)) return;
    const link = t.closest("a.wikilink[data-id], a.cite-link[data-id]");
    if (link && canHover && !link.closest(".popover")) {
      clearTimeout(popTimer);
      popTimer = setTimeout(() => showPopover(link, link.getAttribute("data-id")), 300);
    }
    const tipEl = t.closest("[data-tip]");
    if (tipEl && tipEl.getAttribute("data-tip")) {
      tip.textContent = tipEl.getAttribute("data-tip");
      tip.hidden = false;
    }
  });
  document.addEventListener("mousemove", ev => {
    if (!tip.hidden) {
      const w = tip.offsetWidth;
      tip.style.left = Math.min(ev.clientX + 12, window.innerWidth - w - 8) + "px";
      tip.style.top = (ev.clientY + 16) + "px";
    }
  });
  document.addEventListener("mouseout", ev => {
    const t = ev.target;
    if (!(t instanceof Element)) return;
    if (t.closest("a.wikilink[data-id], a.cite-link[data-id]")) { clearTimeout(popTimer); popTimer = setTimeout(() => { popover.hidden = true; }, 120); }
    if (t.closest("[data-tip]")) hideTip();
  });

  // ============================ Search
  let miniSearch = null;
  function searchIndex() {
    if (miniSearch || !window.MiniSearch) return miniSearch;
    miniSearch = new window.MiniSearch({
      fields: ["title", "aliases", "summary", "notation", "text"],
      storeFields: ["id"],
      processTerm: term => term.toLowerCase().normalize("NFD").replace(/[̀-ͯ]/g, ""),
      searchOptions: { boost: { title: 5, aliases: 4, summary: 2, notation: 1.5, text: 1 }, prefix: true, fuzzy: 0.2 },
    });
    miniSearch.addAll(NOTES.map(n => ({
      id: n.id, title: plainTitle(n.title), aliases: (n.aliases || []).join(" "), summary: n.summary || "",
      notation: (n.notation || []).map(x => stripTags(x.meaning_html)).join(" "), text: n.text || "",
    })));
    return miniSearch;
  }
  function searchNotes(q) {
    const index = searchIndex();
    if (!index || !q.trim()) return [];
    let hits = index.search(q, { combineWith: "AND" });
    if (!hits.length) hits = index.search(q, { combineWith: "OR" });
    return hits.map(h => byId.get(h.id)).filter(Boolean);
  }
  const modal = $("#search-modal");
  const input = $("#search-input");
  const resultsBox = $("#search-results");
  let results = [];
  let selected = 0;
  function renderResults() {
    const q = input.value;
    let groupTitle = "";
    if (q.trim()) {
      results = searchNotes(q).slice(0, 30);
    } else {
      results = NOTES.filter(n => n.type !== "source").slice().sort((a, b) => String(b.updated).localeCompare(String(a.updated))).slice(0, 8);
      groupTitle = results.length ? '<div class="group">Recently updated</div>' : "";
    }
    selected = Math.min(selected, Math.max(0, results.length - 1));
    resultsBox.innerHTML = groupTitle + (results.length ? results.map((n, i) => `
      <a class="result${i === selected ? " sel" : ""}" role="option" aria-selected="${i === selected}" href="#/n/${enc(n.id)}" data-i="${i}">
        ${shapeIcon(n.type, n.area ? areaColor(n.area) : null)}
        <span class="r-title">${n.title_html}<span class="r-meta">${esc(typeLabel(n.type))}${n.domains && n.domains.length ? " · " + esc(domainLabel(n.domains[0])) : ""}</span></span>
        <span class="r-sum">${n.type === "source" ? esc(n.label || "") : (n.summary_html || "")}</span>
      </a>`).join("") : `<div class="group">No note matches "${esc(q)}".</div>`);
    typeset(resultsBox);
  }
  function openSearch() { modal.hidden = false; input.value = ""; selected = 0; renderResults(); input.focus(); }
  function closeSearch() { modal.hidden = true; }
  $("#search-open").addEventListener("click", openSearch);
  input.addEventListener("input", () => { selected = 0; renderResults(); });
  input.addEventListener("keydown", ev => {
    if (ev.key === "ArrowDown") { ev.preventDefault(); selected = Math.min(selected + 1, results.length - 1); renderResults(); scrollSel(); }
    else if (ev.key === "ArrowUp") { ev.preventDefault(); selected = Math.max(selected - 1, 0); renderResults(); scrollSel(); }
    else if (ev.key === "Enter") { ev.preventDefault(); if (results[selected]) { closeSearch(); location.hash = "#/n/" + enc(results[selected].id); } }
    else if (ev.key === "Escape") { closeSearch(); }
  });
  function scrollSel() { const el = resultsBox.querySelector(".result.sel"); if (el) el.scrollIntoView({ block: "nearest" }); }
  resultsBox.addEventListener("click", ev => { if (ev.target.closest(".result")) closeSearch(); });
  modal.addEventListener("click", ev => { if (ev.target === modal) closeSearch(); });
  document.addEventListener("keydown", ev => {
    if ((ev.key === "/" && !isTyping(ev.target)) || ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === "k")) {
      ev.preventDefault();
      openSearch();
    } else if (ev.key === "?" && !isTyping(ev.target)) {
      location.hash = "#/guide";
    } else if (ev.key === "Escape") {
      if (!modal.hidden) closeSearch();
      const card = $("#node-card");
      if (card) card.hidden = true;
      hidePopover();
    }
  });

  // ============================ Graph drawing
  const activeNets = [];
  function destroyNets() { while (activeNets.length) { try { activeNets.pop().destroy(); } catch (e) { /* already destroyed */ } } }
  function graphTheme() {
    return {
      surface: cssVar("--graph-surface"), ink: cssVar("--ink"), ink2: cssVar("--ink-2"), muted: cssVar("--muted"),
      line: cssVar("--line-strong"), accent: cssVar("--accent"), analytic: cssVar("--edge-analytic"), concern: cssVar("--edge-concern"),
      neutral: cssVar("--neutral-node"), matSeed: cssVar("--mat-seed"), matDev: cssVar("--mat-developing"), matMature: cssVar("--mat-mature"),
    };
  }
  function shade(hex, amount) {
    const m = /^#?([0-9a-f]{6})$/i.exec(hex || "");
    if (!m) return hex;
    const v = parseInt(m[1], 16);
    const f = c => Math.max(0, Math.min(255, Math.round(c + (amount < 0 ? c * amount : (255 - c) * amount))));
    return "#" + [f(v >> 16 & 255), f(v >> 8 & 255), f(v & 255)].map(x => x.toString(16).padStart(2, "0")).join("");
  }
  function wrapLabel(text, width) {
    const words = plainTitle(text).split(/\s+/);
    const lines = [];
    let line = "";
    for (const w of words) {
      if ((line + " " + w).trim().length > (width || 18) && line) { lines.push(line); line = w; } else { line = (line + " " + w).trim(); }
    }
    if (line) lines.push(line);
    return lines.join("\n");
  }
  function nodeColor(n, mode, T) {
    if (mode === "maturity") return ({ seed: T.matSeed, developing: T.matDev, mature: T.matMature })[n.maturity] || T.neutral;
    if (mode === "confidence") { const c = level("confidence", n.confidence); return c.background || T.neutral; }
    if (mode === "understanding") { const k = n.understanding || 0; return k ? shade(T.analytic, 0.7 - 0.15 * k) : T.neutral; }
    return n.area ? areaColor(n.area) : T.neutral;
  }
  function visNode(n, mode, T, degree, focus) {
    const shape = shapeOf(n.type);
    const color = nodeColor(n, mode, T);
    const node = {
      id: n.id, shape, label: wrapLabel(n.title),
      color: { background: color, border: shade(color, -0.3), highlight: { background: color, border: T.ink }, hover: { background: color, border: T.ink } },
      font: { color: T.ink, size: 13, face: "Libre Franklin, system-ui, sans-serif", strokeWidth: 4, strokeColor: T.surface },
      borderWidth: focus ? 4 : (n.maturity === "mature" ? 2.5 : 1.5),
    };
    if (focus) node.color.border = T.accent;
    if (shape === "box") {
      node.font = Object.assign({}, node.font, { color: "#ffffff", strokeWidth: 0, size: 12 });
      node.margin = 7;
      node.color.background = T.neutral;
      node.color.highlight.background = T.neutral;
      node.color.hover.background = T.neutral;
    } else {
      node.size = (focus ? 15 : 9) + 3 * Math.sqrt(degree || 0);
    }
    if (n.maturity === "seed") node.shapeProperties = { borderDashes: [3, 3] };
    return node;
  }
  function edgeStyle(e, T) {
    if (isTyped(e)) {
      const r = REL[e.k];
      const fam = {
        structure: { color: T.ink2, dashes: false, width: 1.4 },
        analytic: { color: T.analytic, dashes: false, width: 1.8 },
        application: { color: T.ink2, dashes: [7, 4], width: 1.3 },
        contrast: { color: T.ink2, dashes: [2, 5], width: 1.8 },
        history: { color: T.muted, dashes: [11, 6], width: 1.2 },
      }[r.family] || { color: T.ink2, dashes: false, width: 1.3 };
      return Object.assign(fam, { arrows: !r.symmetric });
    }
    if (e.k === "mentions") return { color: T.line, dashes: false, width: 0.9, arrows: false };
    if (e.k === "cites") return { color: T.muted, dashes: [1, 4], width: 1, arrows: true };
    return { color: T.concern, dashes: [4, 4], width: 1.3, arrows: true };
  }
  const edgeText = e => isTyped(e) ? REL[e.k].label : ({ mentions: "mentions", cites: "cites", concerns: "concerns", "resolved-by": "resolved by" }[e.k] || e.k);
  function visEdge(e, i, T) {
    const st = edgeStyle(e, T);
    return {
      id: "e" + i, from: e.s, to: e.t, _e: e, width: st.width, dashes: st.dashes,
      color: { color: st.color, highlight: T.accent, hover: T.accent, inherit: false, opacity: 0.9 },
      arrows: { to: { enabled: st.arrows, scaleFactor: 0.55 } },
      smooth: { enabled: true, type: "continuous", roundness: 0.3 },
      font: { size: 11, color: T.ink2, strokeWidth: 4, strokeColor: T.surface, align: "horizontal", face: "Jost, system-ui, sans-serif" },
    };
  }
  function networkOptions(hierarchical, count, compact) {
    const base = {
      autoResize: true,
      interaction: { hover: true, tooltipDelay: 400, dragView: true, zoomView: true, multiselect: false, navigationButtons: false, keyboard: false },
      edges: { selectionWidth: 1.2, hoverWidth: 0.6 },
    };
    if (hierarchical) {
      return Object.assign(base, {
        layout: { hierarchical: { enabled: true, direction: "DU", sortMethod: "directed", shakeTowards: "roots", nodeSpacing: compact ? 130 : 180, levelSeparation: compact ? 85 : 115 } },
        physics: { enabled: true, solver: "hierarchicalRepulsion", hierarchicalRepulsion: { nodeDistance: compact ? 120 : 160, avoidOverlap: 0.8 }, stabilization: { enabled: true, iterations: 300 } },
      });
    }
    return Object.assign(base, {
      layout: { improvedLayout: count <= 120, randomSeed: 11 },
      physics: {
        enabled: true, solver: "forceAtlas2Based", minVelocity: 0.75,
        forceAtlas2Based: { gravitationalConstant: compact ? -45 : -75, centralGravity: 0.015, springLength: compact ? 85 : 120, springConstant: 0.07, damping: 0.55, avoidOverlap: 0.7 },
        stabilization: { enabled: true, iterations: 400, fit: true },
      },
    });
  }
  function makeNetwork(container, nodeList, edgeList, opts) {
    opts = opts || {};
    if (!container) return null;
    if (!window.vis) { container.innerHTML = '<p class="empty">The graph library did not load.</p>'; return null; }
    const T = graphTheme();
    const degree = {};
    edgeList.forEach(e => { degree[e.s] = (degree[e.s] || 0) + 1; degree[e.t] = (degree[e.t] || 0) + 1; });
    const nodes = new window.vis.DataSet(nodeList.map(n => visNode(n, opts.color || "area", T, degree[n.id], opts.focus === n.id)));
    const edges = new window.vis.DataSet(edgeList.map((e, i) => visEdge(e, i, T)));
    const net = new window.vis.Network(container, { nodes, edges }, networkOptions(!!opts.hierarchical, nodeList.length, !!opts.compact));
    net.kbData = { nodes, edges };
    net.once("stabilizationIterationsDone", () => {
      net.setOptions({ physics: { enabled: false } });
      net.fit({ animation: false });
      // fit() frames node centres and shapes but not every label, so small
      // graphs are zoomed out slightly to keep labels inside the frame.
      if (opts.compact) net.moveTo({ scale: net.getScale() * 0.82 });
    });
    if (opts.navigate) {
      net.on("click", p => { if (p.nodes.length) location.hash = "#/n/" + enc(p.nodes[0]); });
      net.on("hoverNode", () => { container.style.cursor = "pointer"; });
      net.on("blurNode", () => { container.style.cursor = ""; });
    }
    activeNets.push(net);
    return net;
  }
  function neighbourhood(start, edges, depth) {
    const adj = new Map();
    for (const e of edges) {
      if (!adj.has(e.s)) adj.set(e.s, []);
      if (!adj.has(e.t)) adj.set(e.t, []);
      adj.get(e.s).push(e.t);
      adj.get(e.t).push(e.s);
    }
    const dist = new Map([[start, 0]]);
    const queue = [start];
    while (queue.length) {
      const cur = queue.shift();
      if (dist.get(cur) >= depth) continue;
      for (const nb of adj.get(cur) || []) if (!dist.has(nb)) { dist.set(nb, dist.get(cur) + 1); queue.push(nb); }
    }
    return dist;
  }
  function shortestPaths(a, b, k, includeMentions) {
    const adj = new Map();
    const add = (x, y, e, fwd) => { if (!adj.has(x)) adj.set(x, []); adj.get(x).push([y, e, fwd]); };
    for (const e of EDGES) if (isTyped(e) || (includeMentions && e.k === "mentions")) { add(e.s, e.t, e, true); add(e.t, e.s, e, false); }
    if (a === b) return [[]];
    const dist = new Map([[a, 0]]);
    const preds = new Map();
    const queue = [a];
    while (queue.length) {
      const cur = queue.shift();
      for (const [nb, e, fwd] of adj.get(cur) || []) {
        if (!dist.has(nb)) { dist.set(nb, dist.get(cur) + 1); queue.push(nb); }
        if (dist.get(nb) === dist.get(cur) + 1) { if (!preds.has(nb)) preds.set(nb, []); preds.get(nb).push([cur, e, fwd]); }
      }
    }
    if (!dist.has(b)) return [];
    const found = [];
    const walk = (node, acc) => {
      if (found.length >= k) return;
      if (node === a) { found.push(acc.slice().reverse()); return; }
      for (const [prev, e, fwd] of preds.get(node) || []) walk(prev, acc.concat([[prev, node, e, fwd]]));
    };
    walk(b, []);
    return found;
  }

  // ============================ Chrome
  const TABS = [["explore", "Explore"], ["glossary", "Glossary"], ["symbols", "Symbols"], ["maps", "Maps"], ["paths", "Paths"],
    ["philosophy", "Philosophy"], ["notes", "Notes"], ["references", "References"], ["state", "State"], ["timeline", "Timeline"],
    ["questions", "Questions"], ["recall", "Recall"], ["guide", "Guide"]];
  function renderTabs(active) {
    $("#tabs").innerHTML = TABS.map(([k, label]) => `<a class="tab${k === active ? " active" : ""}" role="tab" aria-selected="${k === active}" href="#/${k}">${label}</a>`).join("");
  }
  function setupChrome() {
    const siteTitle = $("#site-title");
    if (!siteTitle.querySelector(".wordmark")) siteTitle.textContent = META.title || "Knowledge base";
    $("#site-subtitle").textContent = META.subtitle || "";
    // The logo and the home button open the lab website in the same tab.
    const home = $("#home-link");
    const brandHome = $("#brand-home");
    if (META.home) {
      home.href = META.home;
      brandHome.href = META.home;
    } else {
      home.hidden = true;
      brandHome.removeAttribute("href");
    }
    $("#theme-toggle").addEventListener("click", () => {
      const next = isDark() ? "light" : "dark";
      store.set("theme", next);
      document.documentElement.setAttribute("data-theme", next);
      route();
    });
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { if (!store.get("theme", null)) route(); });
    const measure = () => document.documentElement.style.setProperty("--header-h", $(".topbar").offsetHeight + "px");
    measure();
    window.addEventListener("resize", debounce(measure, 150));
    window.addEventListener("load", measure);
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(measure);
    const concepts = NOTES.filter(n => CONCEPT_TYPES.includes(n.type)).length;
    const year = String(META.built || "").slice(0, 4) || new Date().getFullYear();
    $("#footer").innerHTML = [
      `<span>© ${esc(year)} ${esc(META.organisation || META.author || "")}</span>`,
      META.license ? `<span>${META.license_url ? `<a href="${esc(META.license_url)}" target="_blank" rel="noopener">${esc(META.license)}</a>` : esc(META.license)}</span>` : "",
      `<span>${concepts} notes, ${NOTES.filter(n => n.type === "source").length} references</span>`,
      `<span>Built ${esc(String(META.built || "").slice(0, 10))}${LOCAL ? " (local build)" : ""}</span>`,
    ].join("");
    if (LOCAL) {
      const v = A.validation || {};
      const banner = $("#banner");
      banner.innerHTML = `<div class="banner"><strong>Local build</strong><span>Private notes and owner ratings are shown.</span><span>Validation: ${v.errors || 0} errors, ${v.warnings || 0} warnings.</span><a href="#/state">Open the state page</a></div>`;
      banner.hidden = false;
    }
  }

  let teardown = [];
  function clearView() {
    teardown.forEach(fn => { try { fn(); } catch (e) { /* ignore */ } });
    teardown = [];
    destroyNets();
    hidePopover();
    hideTip();
  }
  function setView(tab, html, title, opts) {
    opts = opts || {};
    renderTabs(tab);
    main.classList.toggle("full", !!opts.full);
    $("#banner").style.display = opts.full ? "none" : "";
    main.innerHTML = html;
    document.title = (title ? title + " · " : "") + (META.title || "");
    typeset(main);
    if (!opts.keepScroll) window.scrollTo(0, 0);
  }

  // ============================ Explore
  const exploreDefaults = () => ({ color: "area", hiddenAreas: [], hiddenTypes: ["source"], hiddenRelations: [], hiddenMaturity: [], mentions: true, layout: "network", depth: 0, focus: null });
  const exploreState = () => Object.assign(exploreDefaults(), store.get("explore", {}));
  const saveExplore = st => store.set("explore", st);
  const nodeVisible = (n, st) => !st.hiddenTypes.includes(n.type) && !(n.area && st.hiddenAreas.includes(n.area)) && !(n.maturity && st.hiddenMaturity.includes(n.maturity));
  const edgeVisible = (e, st) => isTyped(e) ? !st.hiddenRelations.includes(e.k) : (e.k === "mentions" ? st.mentions : true);
  let explore = null;

  function showExplore(focusId) {
    const st = exploreState();
    if (focusId && byId.has(focusId)) st.focus = focusId;
    setView("explore", `
      <section class="explore" id="explore">
        <aside class="filters" id="filters" aria-label="Graph filters"></aside>
        <div class="canvas-wrap">
          <div class="graph" id="graph" role="img" aria-label="Graph of notes and their relations"></div>
          <div class="canvas-toolbar">
            <button class="btn filters-toggle" type="button" id="filters-toggle">Filters</button>
            <input type="search" id="find-node" list="note-titles" placeholder="Find a note in the graph" aria-label="Find a note in the graph">
            <button class="btn" type="button" id="g-fit">Fit</button>
            <button class="btn" type="button" id="g-relayout">Rearrange</button>
            <button class="btn" type="button" id="g-png">Save image</button>
            <span class="count" id="g-count"></span>
          </div>
          <div class="legend" id="legend"></div>
          <div class="hint">Click a note for its summary and relations. Double-click opens it. Drag to move, scroll to zoom.</div>
          <div class="node-card" id="node-card" hidden></div>
          <div class="graph-empty" id="graph-empty" hidden>No notes match the filters.</div>
        </div>
      </section>${datalist()}`, "Explore", { full: true });
    renderFilters(st);
    drawExplore(st);
    $("#filters-toggle").addEventListener("click", () => $("#explore").classList.toggle("filters-open"));
    $("#g-fit").addEventListener("click", () => explore && explore.net && explore.net.fit({ animation: { duration: 400 } }));
    $("#g-relayout").addEventListener("click", () => {
      if (!explore || !explore.net) return;
      explore.net.setOptions({ physics: { enabled: true } });
      explore.net.stabilize(300);
      explore.net.once("stabilizationIterationsDone", () => explore.net.setOptions({ physics: { enabled: false } }));
    });
    $("#g-png").addEventListener("click", () => {
      const canvas = $("#graph canvas");
      if (!canvas) return;
      const a = document.createElement("a");
      a.href = canvas.toDataURL("image/png");
      a.download = "knowledge-graph.png";
      document.body.appendChild(a);
      a.click();
      a.remove();
    });
    const find = $("#find-node");
    find.addEventListener("change", () => { const id = findNote(find.value); if (id) focusExplore(id); });
    find.addEventListener("keydown", ev => { if (ev.key === "Enter") { const id = findNote(find.value); if (id) focusExplore(id); } });
  }

  function renderFilters(st) {
    const count = (list, key) => list.reduce((m, n) => { const k = key(n); if (k) m[k] = (m[k] || 0) + 1; return m; }, {});
    const areaCounts = count(NOTES, n => n.area);
    const typeCounts = count(NOTES, n => n.type);
    const matCounts = count(NOTES, n => n.maturity);
    const relCounts = EDGES.reduce((m, e) => { m[e.k] = (m[e.k] || 0) + 1; return m; }, {});
    const row = (kind, id, label, lead, n, on) => `<label class="check-row"><input type="checkbox" data-kind="${kind}" data-id="${esc(id)}"${on ? " checked" : ""}>${lead}<span>${esc(label)}</span>${kind === "area" || kind === "type" ? `<button type="button" class="linkish only" data-only="${kind}:${esc(id)}">Only</button>` : ""}<span class="count">${n || 0}</span></label>`;
    const seg = (name, options, value) => `<div class="seg" role="group">${options.map(([v, l]) => `<button type="button" data-${name}="${v}" class="${v === value ? "on" : ""}">${l}</button>`).join("")}</div>`;
    const T = graphTheme();
    const sample = e => { const s = edgeStyle(e, T); return `<svg class="edge-sample" viewBox="0 0 28 10" aria-hidden="true"><line x1="1" y1="5" x2="27" y2="5" stroke="${s.color}" stroke-width="${Math.max(1.2, s.width)}" ${s.dashes ? `stroke-dasharray="${s.dashes.join(" ")}"` : ""}/></svg>`; };
    const colourModes = [["area", "Area"], ["maturity", "Maturity"], ["confidence", "Confidence"]];
    if (NOTES.some(n => n.understanding)) colourModes.push(["understanding", "Understanding"]);
    const families = ["structure", "analytic", "application", "contrast", "history"];
    const relRows = families.map(f => Object.keys(REL).filter(k => REL[k].family === f && relCounts[k]).map(k => row("relation", k, REL[k].label, sample({ k }), relCounts[k], !st.hiddenRelations.includes(k))).join("")).join("");
    $("#filters").innerHTML = `
      <div class="filter-group"><h3>Colour by</h3>${seg("color", colourModes, st.color)}</div>
      <div class="filter-group"><h3>Layout</h3>${seg("layout", [["network", "Network"], ["prereq", "Prerequisites"]], st.layout)}
        <p class="small muted">Prerequisites draws the requires relation as layers, with foundations at the top.</p></div>
      <div class="filter-group"><h3>Around the selected note</h3>${seg("depth", [["0", "Whole graph"], ["1", "1 step"], ["2", "2 steps"], ["3", "3 steps"]], String(st.depth))}</div>
      <div class="filter-group"><h3>Areas <button type="button" class="linkish" data-all="area">Show all</button></h3>${areaIds().filter(a => areaCounts[a]).map(a => row("area", a, S.areas[a].label, `<span class="swatch" style="background:${areaColor(a)}"></span>`, areaCounts[a], !st.hiddenAreas.includes(a))).join("")}</div>
      <div class="filter-group"><h3>Note types <button type="button" class="linkish" data-all="type">Show all</button></h3>${Object.keys(S.types || {}).filter(t => typeCounts[t]).map(t => row("type", t, typePlural(t), shapeIcon(t), typeCounts[t], !st.hiddenTypes.includes(t))).join("")}</div>
      <div class="filter-group"><h3>Maturity</h3>${levelIds("maturity").filter(m => matCounts[m]).map(m => row("maturity", m, levelLabel("maturity", m), pips(level("maturity", m).rank || 0, 3), matCounts[m], !st.hiddenMaturity.includes(m))).join("")}</div>
      <div class="filter-group"><h3>Relations</h3>${relRows || '<p class="small muted">No typed relations yet.</p>'}
        <label class="check-row"><input type="checkbox" data-kind="mentions"${st.mentions ? " checked" : ""}>${sample({ k: "mentions" })}<span>Mentions in the text</span><span class="count">${relCounts.mentions || 0}</span></label></div>
      <div class="filter-group"><button type="button" class="btn" id="filters-reset">Reset filters</button></div>`;
    const filters = $("#filters");
    const redraw = () => { saveExplore(st); renderFilters(st); drawExplore(st); };
    filters.querySelectorAll("input[type=checkbox]").forEach(cb => cb.addEventListener("change", () => {
      const { kind, id } = cb.dataset;
      const toggle = (list, value, visible) => { const i = list.indexOf(value); if (visible && i >= 0) list.splice(i, 1); if (!visible && i < 0) list.push(value); };
      if (kind === "area") toggle(st.hiddenAreas, id, cb.checked);
      if (kind === "type") toggle(st.hiddenTypes, id, cb.checked);
      if (kind === "maturity") toggle(st.hiddenMaturity, id, cb.checked);
      if (kind === "relation") toggle(st.hiddenRelations, id, cb.checked);
      if (kind === "mentions") st.mentions = cb.checked;
      redraw();
    }));
    filters.querySelectorAll("[data-only]").forEach(b => b.addEventListener("click", ev => {
      ev.preventDefault();
      const [kind, id] = b.dataset.only.split(":");
      if (kind === "area") st.hiddenAreas = areaIds().filter(a => a !== id);
      if (kind === "type") st.hiddenTypes = Object.keys(S.types || {}).filter(t => t !== id);
      redraw();
    }));
    filters.querySelectorAll("[data-all]").forEach(b => b.addEventListener("click", () => {
      if (b.dataset.all === "area") st.hiddenAreas = [];
      if (b.dataset.all === "type") st.hiddenTypes = [];
      redraw();
    }));
    filters.querySelectorAll("[data-color]").forEach(b => b.addEventListener("click", () => { st.color = b.dataset.color; redraw(); }));
    filters.querySelectorAll("[data-layout]").forEach(b => b.addEventListener("click", () => { st.layout = b.dataset.layout; redraw(); }));
    filters.querySelectorAll("[data-depth]").forEach(b => b.addEventListener("click", () => { st.depth = Number(b.dataset.depth); redraw(); }));
    $("#filters-reset").addEventListener("click", () => { Object.assign(st, exploreDefaults()); redraw(); });
  }

  function drawExplore(st) {
    let nodes = NOTES.filter(n => nodeVisible(n, st));
    const ids = new Set(nodes.map(n => n.id));
    let edges = EDGES.filter(e => ids.has(e.s) && ids.has(e.t) && edgeVisible(e, st));
    if (st.layout === "prereq") {
      edges = edges.filter(e => e.k === "requires");
      const involved = new Set();
      edges.forEach(e => { involved.add(e.s); involved.add(e.t); });
      nodes = nodes.filter(n => involved.has(n.id));
    }
    if (st.depth > 0 && st.focus && nodes.some(n => n.id === st.focus)) {
      const keep = neighbourhood(st.focus, edges, st.depth);
      nodes = nodes.filter(n => keep.has(n.id));
      edges = edges.filter(e => keep.has(e.s) && keep.has(e.t));
    }
    $("#g-count").textContent = `${nodes.length} notes · ${edges.length} links`;
    $("#graph-empty").hidden = nodes.length > 0;
    $("#node-card").hidden = true;
    destroyNets();
    const net = makeNetwork($("#graph"), nodes, edges, { color: st.color, hierarchical: st.layout === "prereq", focus: st.focus });
    renderLegend(st, nodes, edges);
    if (!net) return;
    explore = { net, st, labelled: [] };
    net.on("click", p => { if (p.nodes.length) selectExplore(p.nodes[0], p.pointer.DOM); else closeNodeCard(); });
    net.on("doubleClick", p => { if (p.nodes.length) location.hash = "#/n/" + enc(p.nodes[0]); });
    net.on("hoverNode", () => { $("#graph").style.cursor = "pointer"; });
    net.on("blurNode", () => { $("#graph").style.cursor = ""; });
    net.once("stabilizationIterationsDone", () => {
      if (st.focus && nodes.some(n => n.id === st.focus)) {
        net.selectNodes([st.focus]);
        labelEdges(st.focus);
      }
    });
  }

  function labelEdges(id) {
    if (!explore) return;
    const data = explore.net.kbData;
    if (explore.labelled.length) data.edges.update(explore.labelled.map(eid => ({ id: eid, label: "" })));
    const ids = explore.net.getConnectedEdges(id);
    data.edges.update(ids.map(eid => ({ id: eid, label: edgeText(data.edges.get(eid)._e) })));
    explore.labelled = ids;
  }
  function closeNodeCard() {
    const card = $("#node-card");
    if (card) card.hidden = true;
    if (explore && explore.labelled.length) {
      explore.net.kbData.edges.update(explore.labelled.map(eid => ({ id: eid, label: "" })));
      explore.labelled = [];
    }
  }
  function nodeCardHtml(n) {
    const groups = [...relationGroups(n.id).entries()];
    const rels = groups.slice(0, 10).map(([label, items]) => `<div class="nc-rel"><span class="nc-label">${esc(label)}:</span> ${items.map(it => `<button type="button" class="linkish" data-focus="${esc(it.id)}">${esc(plainTitle((byId.get(it.id) || {}).title || it.id))}</button>`).join(", ")}</div>`).join("");
    const prereq = (A.prerequisites || {})[n.id];
    return `<div class="nc-head"><div class="chips">${typeChip(n)}${n.area ? `<span class="chip"><span class="swatch" style="background:${areaColor(n.area)}"></span>${esc(S.areas[n.area].label)}</span>` : ""}</div><button class="icon-btn nc-close" type="button" aria-label="Close">×</button></div>
      <h3>${n.title_html}</h3>
      <div class="nc-summary">${n.type === "source" ? (n.reference_html || "") : (n.summary_html || "")}</div>
      <div class="chips">${levelChips(n)}</div>
      ${rels ? `<div class="nc-rels">${rels}</div>` : ""}
      <div class="nc-actions"><a class="btn btn-primary" href="#/n/${enc(n.id)}">Open note</a>${prereq ? `<a class="btn" href="#/paths/prereq/${enc(n.id)}">Prerequisites (${prereq.length})</a>` : ""}</div>`;
  }
  function selectExplore(id, pos) {
    const n = byId.get(id);
    if (!n || !explore) return;
    explore.st.focus = id;
    saveExplore(explore.st);
    labelEdges(id);
    const card = $("#node-card");
    card.innerHTML = nodeCardHtml(n);
    card.hidden = false;
    typeset(card);
    const wrap = card.parentElement.getBoundingClientRect();
    const w = card.offsetWidth;
    const h = card.offsetHeight;
    let left = pos ? pos.x + 18 : wrap.width - w - 12;
    let top = pos ? pos.y - 20 : 60;
    if (left + w > wrap.width - 12) left = (pos ? pos.x : wrap.width) - w - 18;
    left = Math.max(12, left);
    top = Math.max(56, Math.min(top, wrap.height - h - 12));
    card.style.left = left + "px";
    card.style.top = top + "px";
    card.querySelector(".nc-close").addEventListener("click", closeNodeCard);
    card.querySelectorAll("[data-focus]").forEach(b => b.addEventListener("click", () => focusExplore(b.dataset.focus)));
  }
  function focusExplore(id) {
    if (!explore) return;
    if (!explore.net.kbData.nodes.get(id)) { location.hash = "#/n/" + enc(id); return; }
    explore.net.selectNodes([id]);
    explore.net.focus(id, { scale: 1.15, animation: { duration: 350 } });
    selectExplore(id, null);
  }
  function renderLegend(st, nodes, edges) {
    const legend = $("#legend");
    if (!legend) return;
    const T = graphTheme();
    let colour = "";
    if (st.color === "area") {
      const present = new Set(nodes.map(n => n.area).filter(Boolean));
      colour = areaIds().filter(a => present.has(a)).map(a => `<button type="button" class="item" data-legend-area="${esc(a)}" title="Hide or show this area"><span class="swatch" style="background:${areaColor(a)}"></span>${esc(S.areas[a].label)}</button>`).join("");
    } else if (st.color === "maturity") {
      colour = [["seed", T.matSeed], ["developing", T.matDev], ["mature", T.matMature]].map(([k, c]) => `<span class="item"><span class="swatch" style="background:${c}"></span>${esc(levelLabel("maturity", k))}</span>`).join("");
    } else if (st.color === "confidence") {
      colour = Object.keys(LEVELS.confidence || {}).map(k => `<span class="item"><span class="swatch" style="background:${level("confidence", k).background}"></span>${esc(levelLabel("confidence", k))}</span>`).join("");
    } else {
      colour = '<span class="item">Lighter means lower understanding</span>';
    }
    const types = [...new Set(nodes.map(n => n.type))].map(t => `<span class="item">${shapeIcon(t)}${esc(typeLabel(t))}</span>`).join("");
    legend.innerHTML = colour + types + (nodes.some(n => n.maturity === "seed") ? '<span class="item muted">Dashed outline: seed</span>' : "");
    legend.querySelectorAll("[data-legend-area]").forEach(b => b.addEventListener("click", () => {
      const a = b.dataset.legendArea;
      const i = st.hiddenAreas.indexOf(a);
      if (i >= 0) st.hiddenAreas.splice(i, 1); else st.hiddenAreas.push(a);
      saveExplore(st);
      renderFilters(st);
      drawExplore(st);
    }));
  }

  // ============================ Note view
  function showNote(id, anchor) {
    const n = byId.get(id);
    if (!n) return showWanted(id);
    const tab = { map: "maps", source: "references", question: "questions" }[n.type] || "notes";
    const crumbs = [
      n.area ? `<span class="crumb-area"><span class="swatch" style="background:${areaColor(n.area)}"></span>${esc(S.areas[n.area].label)}</span>` : "",
      ...(n.domains || []).map(d => `<span>${esc(domainLabel(d))}</span>`),
    ].filter(Boolean).join('<span class="sep">/</span>');
    const aliases = (n.aliases || []).length ? `<p class="aliases">Also called ${n.aliases.map(esc).join(", ")}</p>` : "";
    const notation = (n.notation || []).length ? `<div class="box"><h2>Notation <a class="small" href="#/symbols">All symbols</a></h2><div class="table-wrap"><table class="data-table notation-table"><tbody>${n.notation.map(x => `<tr><td><a class="sym-link" href="#/symbols/${enc(x.key || "")}" data-tip="Open in the symbol atlas"><span class="math math-inline">${esc(x.symbol)}</span></a></td><td>${x.meaning_html}</td></tr>`).join("")}</tbody></table></div></div>` : "";
    let special = "";
    if (n.type === "source") {
      const v = n.verification || {};
      const contexts = (n.contexts || []).length ? `<h2 class="box-sub">Citing sentences</h2><ul class="ctx-list">${n.contexts.map(c => `<li><div class="small muted">${noteLink(c.note)}${c.locator ? `, ${esc(c.locator)}` : ""}</div><div class="ctx">${c.html}</div></li>`).join("")}</ul>` : "";
      special = `<div class="box"><h2>Reference</h2><p class="reference">${n.reference_html || ""}</p>${v.note ? `<p class="small muted">${esc(v.note)}</p>` : ""}<div class="ref-actions"><a class="btn" href="#/references/${enc(n.id)}">Show in the references</a><button type="button" class="btn" data-copy="${esc(n.id)}">Copy BibTeX</button></div>${contexts}</div>`;
    }
    if (n.type === "question") {
      special = `<div class="box"><h2>Concerns</h2><div class="chips">${(n.concerns || []).map(c => noteLink(c)).join(" · ") || "None listed"}</div>${(n.resolved_by || []).length ? `<h2 style="margin-top:10px">Answered in</h2><div class="chips">${n.resolved_by.map(c => noteLink(c)).join(" · ")}</div>` : ""}</div>`;
    }
    const checks = checksSection(n);
    const cites = citationsHtml(n);
    const meta = [
      n.created ? `<span>Created ${esc(n.created)}</span>` : "",
      n.updated ? `<span>Updated ${esc(n.updated)}</span>` : "",
      n.reviewed ? `<span>Reviewed by the owner ${esc(n.reviewed)}</span>` : "",
      `<span>File <code>${esc(n.path)}</code></span>`,
    ].join("");
    const html = `
      <div class="note-layout">
        <article class="note">
          <nav class="crumbs" aria-label="Area and domains">${crumbs}</nav>
          <h1 class="note-title">${n.title_html}</h1>
          ${aliases}
          <div class="chips">${typeChip(n)}${levelChips(n)}</div>
          ${n.summary_html ? `<p class="lead">${n.summary_html}</p>` : ""}
          ${special}${notation}
          <div class="prose">${n.html || ""}</div>
          ${checks}${cites}
          <footer class="note-meta">${meta}</footer>
        </article>
        <aside class="note-aside">${asideHtml(n)}</aside>
      </div>`;
    setView(tab, html, plainTitle(n.title));
    decoratePhilosophy(n);
    bindCopy(main);
    const graphBox = $("#local-graph");
    if (graphBox) drawLocalGraph(n, 1);
    $$("[data-local-depth]").forEach(b => b.addEventListener("click", () => {
      $$("[data-local-depth]").forEach(x => x.classList.toggle("on", x === b));
      destroyNets();
      drawLocalGraph(n, Number(b.dataset.localDepth));
    }));
    if (anchor) {
      const target = document.getElementById(anchor);
      if (target) setTimeout(() => target.scrollIntoView({ behavior: "smooth", block: "start" }), 60);
    }
  }
  function decoratePhilosophy(n) {
    const heading = main.querySelector(".prose #ontology-and-epistemology");
    if (!heading) return;
    const wrap = document.createElement("section");
    wrap.className = "philo";
    heading.parentNode.insertBefore(wrap, heading);
    let el = heading;
    while (el) {
      const next = el.nextElementSibling;
      wrap.appendChild(el);
      if (!next || next.tagName === "H2") break;
      el = next;
    }
    if ((n.themes || []).length) heading.insertAdjacentHTML("afterend", `<div class="chips philo-themes">${n.themes.map(themeChip).join("")}</div>`);
  }
  function copyText(text, button) {
    if (!button.dataset.label) button.dataset.label = button.textContent;
    const done = label => { button.textContent = label; setTimeout(() => { button.textContent = button.dataset.label; }, 1600); };
    const fallback = () => {
      const area = document.createElement("textarea");
      area.value = text;
      document.body.appendChild(area);
      area.select();
      let ok = false;
      try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
      area.remove();
      done(ok ? "Copied" : "Copy failed");
    };
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(text).then(() => done("Copied"), fallback);
    else fallback();
  }
  function bindCopy(root) {
    root.querySelectorAll("[data-copy]").forEach(b => b.addEventListener("click", () => {
      const s = byId.get(b.dataset.copy);
      copyText((s && s.bibtex) || "", b);
    }));
  }
  function bindCopyTex(root) {
    root.querySelectorAll("[data-copytex]").forEach(b => {
      if (b.dataset.bound) return;
      b.dataset.bound = "1";
      b.addEventListener("click", ev => {
        ev.preventDefault();
        const s = symbolByKey.get(b.dataset.copytex);
        copyText((s && s.tex) || "", b);
      });
    });
  }
  // One Validation section at the end of a note: the prose of each check, written in the body
  // as a callout and moved here by the build, followed by the values that the run recorded, the
  // figure it drew and the script itself.
  function checksSection(n) {
    const cards = n.checks || [];
    const notes = (n.check_notes || []).slice();
    if (!cards.length && !notes.length) return "";
    const parts = [];
    for (const c of cards) {
      let prose = "";
      const k = notes.findIndex(x => x.script && c.script && x.script === c.script);
      if (k >= 0) prose = notes.splice(k, 1)[0].html;
      parts.push(`<div class="check-block">${prose}${checkHtml(c)}</div>`);
    }
    for (const rest of notes) parts.push(rest.html);
    return `<section class="box" id="validation"><h2>Validation</h2>${parts.join("")}</section>`;
  }
  function checkHtml(c) {
    const cls = c.stale ? "stale" : c.status;
    const label = c.stale ? "Stale: the script changed after the last run" : ({ pass: "✓ Pass", fail: "✕ Fail", error: "✕ Error", "not-run": "! Not run" }[c.status] || c.status);
    const metrics = Object.entries(c.metrics || {}).map(([k, v]) => `<tr><td>${esc(k.replace(/_/g, " "))}</td><td>${esc(v)}</td></tr>`).join("");
    return `<div class="check ${esc(cls)}"><div class="check-head"><span class="check-status">${esc(label)}</span><code>${esc(c.id)}</code>${c.date ? `<span class="muted small">${esc(String(c.date).slice(0, 10))}</span>` : ""}${c.runtime_s != null ? `<span class="muted small">${esc(c.runtime_s)} s</span>` : ""}</div>
      ${c.summary ? `<div class="small">${esc(c.summary)}</div>` : ""}
      ${metrics ? `<details><summary class="small">Measured values</summary><table class="metrics"><tbody>${metrics}</tbody></table></details>` : ""}
      ${c.figure ? `<figure class="check-figure"><a class="figure-download" href="${esc(c.figure)}" download="${esc(c.id)}.svg" title="Click the figure to download it"><img src="${esc(c.figure)}" alt="Figure drawn by the check ${esc(c.id)}" loading="lazy"></a><figcaption class="small muted">Drawn by the run recorded above. Click the figure to download it as <code>${esc(c.id)}.svg</code>.</figcaption></figure>` : ""}
      ${c.script || c.code ? `<details class="check-code"><summary class="code-summary">Show the script${c.code ? ` <span class="muted small">${esc(String(c.code).split("\n").length)} lines</span>` : ""}</summary><pre class="code-block"><code>${c.code ? highlightCode(c.code) : ""}</code></pre></details>` : ""}
      ${c.script ? `<div class="code-actions"><a class="btn btn-small" href="${esc(c.script)}" download="${esc(c.script.split("/").pop())}">Download the script</a><code class="small muted">${esc(c.script)}</code></div>` : ""}</div>`;
  }
  function citationsHtml(n) {
    const grouped = new Map();
    for (const c of n.citations || []) {
      if (!grouped.has(c.key)) grouped.set(c.key, []);
      if (c.locator && !grouped.get(c.key).includes(c.locator)) grouped.get(c.key).push(c.locator);
    }
    if (!grouped.size) return "";
    return `<section class="box"><h2>References</h2><ul class="cited-list">${[...grouped.entries()].map(([k, locs]) => {
      const s = byId.get(k);
      return `<li>${s ? `<a class="cite-link" href="#/n/${enc(k)}" data-id="${esc(k)}">${esc(s.label || k)}</a> ${verificationChip(s)}${locs.length ? ` <span class="muted small">${locs.map(esc).join("; ")}</span>` : ""}<div class="small">${s.reference_html || ""}</div>` : esc(k)}</li>`;
    }).join("")}</ul></section>`;
  }
  function asideHtml(n) {
    const parts = [];
    const toc = (n.toc || []).filter(t => t.level <= 3 && t.id !== "validation");
    // The Validation section is rendered by the site at the end of the note, so its entry is
    // added here rather than taken from the body, where the heading may no longer stand.
    if ((n.checks || []).length || (n.check_notes || []).length) toc.push({ level: 2, id: "validation", html: "Validation" });
    if (toc.length > 1) parts.push(`<section class="panel toc"><h2>On this page</h2>${toc.map(t => `<a class="l${t.level}" href="#/n/${enc(n.id)}/${enc(t.id)}">${t.html}</a>`).join("")}</section>`);
    const groups = [...relationGroups(n.id).entries()];
    if (groups.length) {
      parts.push(`<section class="panel"><h2>Relations</h2><dl class="rel-groups">${groups.map(([label, items]) => `<dt>${esc(label)}</dt>${items.map(it => `<dd>${noteLink(it.id)}${it.note ? `<div class="rel-note">${it.note}</div>` : ""}</dd>`).join("")}`).join("")}</dl></section>`);
    }
    const hasNeighbours = (outE.get(n.id) || []).length + (inE.get(n.id) || []).length > 0;
    if (hasNeighbours) parts.push(`<section class="panel"><h2>Local graph <span class="seg"><button type="button" class="on" data-local-depth="1">1 step</button><button type="button" data-local-depth="2">2 steps</button></span></h2><div class="mini-graph" id="local-graph"></div></section>`);
    const prereq = (A.prerequisites || {})[n.id];
    if (prereq && prereq.length) {
      parts.push(`<section class="panel"><h2>Prerequisites <a class="small" href="#/paths/prereq/${enc(n.id)}">Open path</a></h2><ol>${prereq.map(([p, lv]) => `<li>${noteLink(p)} <span class="muted small">level ${lv}</span></li>`).join("")}</ol></section>`);
    }
    const mentions = (inE.get(n.id) || []).filter(e => e.k === "mentions").map(e => e.s);
    if (mentions.length) parts.push(`<section class="panel"><h2>Mentioned in</h2><ul>${mentions.map(m => `<li>${noteLink(m)}</li>`).join("")}</ul></section>`);
    if (n.type === "map") {
      const members = (outE.get(n.id) || []).filter(e => e.k === "mentions").map(e => e.t);
      if (members.length) parts.push(`<section class="panel"><h2>Notes in this map</h2><ul>${members.map(m => `<li>${noteLink(m)}</li>`).join("")}</ul></section>`);
    }
    const questions = (inE.get(n.id) || []).filter(e => e.k === "concerns").map(e => byId.get(e.s)).filter(Boolean);
    if (questions.length) parts.push(`<section class="panel"><h2>Questions</h2><ul>${questions.map(q => `<li>${noteLink(q.id)} ${statusChip(q)}</li>`).join("")}</ul></section>`);
    if (n.type === "source") {
      const citing = (inE.get(n.id) || []).filter(e => e.k === "cites").map(e => byId.get(e.s)).filter(Boolean);
      parts.push(`<section class="panel"><h2>Cited by</h2>${citing.length ? `<ul>${citing.map(c => {
        const locs = (c.citations || []).filter(x => x.key === n.id && x.locator).map(x => x.locator);
        return `<li>${noteLink(c.id)}${locs.length ? ` <span class="muted small">${locs.map(esc).join("; ")}</span>` : ""}</li>`;
      }).join("")}</ul>` : '<p class="small muted">No note cites this source yet.</p>'}</section>`);
    }
    return parts.join("");
  }
  function drawLocalGraph(n, depth) {
    const box = $("#local-graph");
    if (!box) return;
    const usable = EDGES.filter(e => isTyped(e) || e.k === "mentions" || e.k === "concerns" || (n.type === "source" && e.k === "cites"));
    const dist = neighbourhood(n.id, usable, depth);
    const nodes = NOTES.filter(x => dist.has(x.id));
    const edges = usable.filter(e => dist.has(e.s) && dist.has(e.t));
    makeNetwork(box, nodes, edges, { focus: n.id, compact: true, navigate: true });
  }

  function showWanted(id) {
    const from = ((A.wanted || []).find(w => w.id === id) || {}).from || [];
    setView("notes", `<div class="empty"><h1>Not written yet</h1><p>No note has the id <code>${esc(id)}</code>.</p>
      ${from.length ? `<p>It is linked from ${from.map(f => noteLink(f)).join(", ")}.</p>` : ""}
      ${LOCAL ? `<p class="small">Create it with <code>python3 tools/kb.py new concept "Title" --id ${esc(id)}</code>.</p>` : ""}
      <p><a class="btn" href="#/notes">All notes</a></p></div>`, "Not written yet");
  }

  // ============================ Notes index
  const notesState = () => Object.assign({ q: "", types: [], areas: [], maturity: [], sort: "title" }, store.get("notes", {}));
  function showNotes() {
    const st = notesState();
    const concepts = NOTES.filter(n => CONCEPT_TYPES.includes(n.type));
    setView("notes", `<div class="page-head"><div><h1>Notes</h1><p>${concepts.length} notes on concepts, theorems, methods, models and worked examples. The <a href="#/glossary">glossary</a> lists the same notes by name and the <a href="#/symbols">symbols</a> list their notation.</p></div></div><div id="notes-body"></div>`, "Notes");
    const body = $("#notes-body");
    const chip = (kind, id, label, lead) => `<button type="button" class="filter-chip${st[kind].includes(id) ? "" : " off"}" data-kind="${kind}" data-id="${esc(id)}">${lead || ""}${esc(label)}</button>`;
    const typesPresent = CONCEPT_TYPES.filter(t => concepts.some(n => n.type === t));
    if (!st.types.length) st.types = typesPresent.slice();
    const areasPresent = areaIds().filter(a => concepts.some(n => n.area === a));
    if (!st.areas.length) st.areas = areasPresent.slice();
    if (!st.maturity.length) st.maturity = Object.keys(LEVELS.maturity || {});
    body.innerHTML = `<div class="toolbar">
        <input type="search" class="grow" id="notes-q" placeholder="Filter by title, alias or summary" value="${esc(st.q)}" aria-label="Filter notes">
        <select id="notes-sort" aria-label="Sort notes">
          ${[["title", "Title"], ["updated", "Recently updated"], ["degree", "Most connected"], ["maturity", "Maturity"]].map(([v, l]) => `<option value="${v}"${st.sort === v ? " selected" : ""}>${l}</option>`).join("")}
        </select>
      </div>
      <div class="toolbar">${typesPresent.map(t => chip("types", t, typePlural(t), shapeIcon(t))).join("")}</div>
      <div class="toolbar">${areasPresent.map(a => chip("areas", a, S.areas[a].label, `<span class="swatch" style="background:${areaColor(a)}"></span>`)).join("")}
        ${Object.keys(LEVELS.maturity || {}).map(m => chip("maturity", m, levelLabel("maturity", m), pips(level("maturity", m).rank, 3))).join("")}</div>
      <div class="grid" id="cards"></div>`;
    const draw = () => {
      store.set("notes", st);
      const q = st.q.trim().toLowerCase();
      let list = concepts.filter(n => st.types.includes(n.type) && (!n.area || st.areas.includes(n.area)) && (!n.maturity || st.maturity.includes(n.maturity)));
      if (q) {
        const hits = new Set(searchNotes(q).map(n => n.id));
        list = list.filter(n => hits.has(n.id) || plainTitle(n.title).toLowerCase().includes(q));
      }
      const degree = A.degree || {};
      const sorters = {
        title: (a, b) => plainTitle(a.title).localeCompare(plainTitle(b.title)),
        updated: (a, b) => String(b.updated).localeCompare(String(a.updated)),
        degree: (a, b) => (degree[b.id] || 0) - (degree[a.id] || 0),
        maturity: (a, b) => (level("maturity", b.maturity).rank || 0) - (level("maturity", a.maturity).rank || 0),
      };
      list.sort(sorters[st.sort] || sorters.title);
      $("#cards").innerHTML = list.length ? list.map(cardHtml).join("") : '<p class="muted">No note matches the filters.</p>';
      typeset($("#cards"));
    };
    $("#notes-q").addEventListener("input", debounce(ev => { st.q = ev.target.value; draw(); }, 120));
    $("#notes-sort").addEventListener("change", ev => { st.sort = ev.target.value; draw(); });
    body.querySelectorAll(".filter-chip").forEach(b => b.addEventListener("click", () => {
      const list = st[b.dataset.kind];
      const i = list.indexOf(b.dataset.id);
      if (i >= 0) list.splice(i, 1); else list.push(b.dataset.id);
      b.classList.toggle("off", i >= 0);
      draw();
    }));
    draw();
  }
  function cardHtml(n) {
    const degree = (A.degree || {})[n.id] || 0;
    return `<a class="card" href="#/n/${enc(n.id)}" style="--card-color:${n.area ? areaColor(n.area) : "var(--line-strong)"}">
      <div class="chips">${typeChip(n)}${maturityChip(n)}</div>
      <h3>${n.title_html}</h3>
      <div class="summary">${n.summary_html || ""}</div>
      <div class="foot">${confidenceChip(n)}<span>${degree} ${degree === 1 ? "link" : "links"}</span>${(n.domains || []).slice(0, 2).map(d => `<span>${esc(domainLabel(d))}</span>`).join("")}</div>
    </a>`;
  }
  function showGlossary() {
    const concepts = NOTES.filter(n => CONCEPT_TYPES.includes(n.type));
    setView("glossary", `<div class="page-head"><div><h1>Glossary</h1><p>Every concept, theorem, method, model and worked example by name, with the other names it goes by and the sentence that defines it. The <a href="#/notes">notes</a> page filters the same set and the <a href="#/symbols">symbols</a> page lists their notation.</p></div></div>
      <div class="toolbar"><input type="search" class="grow" id="gloss-q" placeholder="Filter by name, alias or summary" aria-label="Filter the glossary"></div>
      <div id="gloss-body"></div>`, "Glossary");
    const body = $("#gloss-body");
    const draw = q => {
      const needle = q.trim().toLowerCase();
      const shown = !needle ? concepts : concepts.filter(n =>
        [plainTitle(n.title), (n.aliases || []).join(" "), stripTags(n.summary_html || "")].join(" ").toLowerCase().includes(needle));
      renderGlossary(body, shown);
    };
    draw("");
    const q = $("#gloss-q");
    if (q) q.addEventListener("input", () => draw(q.value));
  }

  function renderGlossary(body, concepts) {
    const sorted = concepts.slice().sort((a, b) => plainTitle(a.title).localeCompare(plainTitle(b.title)));
    const groups = new Map();
    for (const n of sorted) {
      const letter = (plainTitle(n.title)[0] || "#").toUpperCase().normalize("NFD").replace(/[̀-ͯ]/g, "");
      if (!groups.has(letter)) groups.set(letter, []);
      groups.get(letter).push(n);
    }
    body.innerHTML = `<div class="letters">${[...groups.keys()].map(l => `<a href="#/glossary" data-letter="${esc(l)}">${esc(l)}</a>`).join("")}</div>
      <div class="gloss">${[...groups.entries()].map(([l, list]) => `<h2 id="letter-${esc(l)}">${esc(l)}</h2><dl>${list.map(n => `<dt>${noteLink(n.id)}${(n.aliases || []).length ? ` <span class="muted small">also ${n.aliases.map(esc).join(", ")}</span>` : ""}</dt><dd>${n.summary_html || ""}</dd>`).join("")}</dl>`).join("") || '<p class="muted">No notes yet.</p>'}</div>`;
    body.querySelectorAll("[data-letter]").forEach(a => a.addEventListener("click", ev => { ev.preventDefault(); const el = document.getElementById("letter-" + a.dataset.letter); if (el) el.scrollIntoView({ behavior: "smooth" }); }));
    typeset(body);
  }
  // ============================ Symbols
  const FAMILY_ORDER = ["Latin", "Greek", "Calligraphic", "Blackboard and bold", "Structures", "Relations", "Operators and accents"];
  function showSymbols(focusKey) {
    const st = Object.assign({ q: "", view: "atlas", group: "family", hiddenFamilies: [], selected: null }, store.get("symbols", {}));
    if (focusKey && symbolByKey.has(focusKey)) Object.assign(st, { view: "atlas", q: "", hiddenFamilies: [], selected: focusKey });
    let selected = st.selected && symbolByKey.has(st.selected) ? st.selected : ((SYMBOLS[0] || {}).key || null);
    const areaOf = s => (s.note && byId.has(s.note) ? byId.get(s.note).area : null);
    const colorOf = s => (areaOf(s) ? areaColor(areaOf(s)) : cssVar("--accent"));
    const noteTitleOf = s => (s.note && byId.has(s.note) ? plainTitle(byId.get(s.note).title) : "");
    const conventionLabel = "Conventions without a defining note";
    const families = FAMILY_ORDER.filter(f => SYMBOLS.some(s => s.family === f));
    const uses = SYMBOLS.reduce((a, s) => a + (s.uses || []).length, 0);
    const macros = Object.entries(MACROS);
    setView("symbols", `<div class="page-head"><div><h1>Symbols</h1><p>Every symbol of the knowledge base: what it denotes, the convention that governs it, the note that defines it and the formulas of other notes in which it appears.</p></div></div>
      <div class="kpis">${kpiHtml(SYMBOLS.length, "Symbols")}${kpiHtml(SYMBOLS.filter(s => s.note).length, "Defined in notes")}${kpiHtml(SYMBOLS.filter(s => s.section).length, "Governed by a convention")}${kpiHtml(uses, "Notes that reuse them")}</div>
      <div class="toolbar"><input type="search" class="grow" id="sym-q" placeholder="Filter by meaning, TeX or note" value="${esc(st.q)}" aria-label="Filter symbols">
        <select id="sym-group" aria-label="Group symbols">${[["family", "Group by typography"], ["area", "Group by area"], ["note", "Group by defining note"], ["none", "No groups"]].map(([v, l]) => `<option value="${v}"${st.group === v ? " selected" : ""}>${l}</option>`).join("")}</select>
        <div class="seg" role="group">${[["atlas", "Atlas"], ["table", "Table"]].map(([v, l]) => `<button type="button" data-symview="${v}" class="${st.view === v ? "on" : ""}">${l}</button>`).join("")}</div></div>
      <div class="toolbar">${families.map(f => `<button type="button" class="filter-chip${st.hiddenFamilies.includes(f) ? " off" : ""}" data-family="${esc(f)}">${esc(f)} <span class="muted">${SYMBOLS.filter(s => s.family === f).length}</span></button>`).join("")}</div>
      <div id="sym-body"></div>
      ${macros.length ? `<section class="panel macro-panel"><h2>TeX macros</h2><p class="small muted">Notes may write these commands, and the site and the TeX check expand them in the same way.</p><div class="macro-grid">${macros.map(([name, expansion]) => `<div class="macro-card"><span class="macro-glyph"><span class="math math-inline">${esc(expansion.includes("#") ? `${name}{x}` : name)}</span></span><code>${esc(expansion.includes("#") ? `${name}{x}` : name)}</code><code class="muted">${esc(expansion)}</code></div>`).join("")}</div></section>` : ""}`, "Symbols");

    const matches = s => {
      if (st.hiddenFamilies.includes(s.family)) return false;
      const q = st.q.trim().toLowerCase();
      if (!q) return true;
      return [s.meaning, s.tex, stripTags(s.object_html || ""), s.section || "", noteTitleOf(s)].join(" ").toLowerCase().includes(q);
    };
    const groupOf = s => st.group === "family" ? s.family
      : st.group === "area" ? (areaOf(s) ? S.areas[areaOf(s)].label : conventionLabel)
      : st.group === "note" ? (noteTitleOf(s) || conventionLabel) : "";
    const orderGroups = keys => st.group === "family" ? FAMILY_ORDER.filter(f => keys.includes(f))
      : keys.sort((a, b) => (a === conventionLabel) - (b === conventionLabel) || a.localeCompare(b));
    // The size of a glyph is chosen from an estimate of how wide it will render, not from the
    // length of its TeX: a command such as \mathcal draws one letter and a plain string such as
    // "L(T), B(T)" draws ten, so counting the source overflows the tile for the second kind.
    // Each command counts as one character and the braces and the spacing commands as none.
    const glyphWeight = tex => String(tex).replace(/\\[A-Za-z]+/g, "x").replace(/[{}\\\s]/g, "").length;
    const glyphSize = s => { const w = glyphWeight(s.tex); return w > 14 ? " tiny" : w > 7 ? " long" : ""; };
    const isLong = s => glyphSize(s) !== "";
    const specimen = s => {
      if (!s) return '<p class="muted">No symbol matches.</p>';
      const note = s.note ? byId.get(s.note) : null;
      const list = s.uses || [];
      return `<div class="specimen" style="--sym-color:${colorOf(s)}">
        <div class="specimen-stage"><span class="specimen-family">${esc(s.family)}</span><span class="specimen-glyph${glyphSize(s)}"><span class="math math-inline">${esc(s.tex)}</span></span></div>
        <div class="specimen-body">
          <p class="specimen-meaning">${s.meaning_html || ""}</p>
          ${s.note && s.object_html ? `<p class="small muted">Listed in the notation conventions (${esc(s.section || "")}) as ${s.object_html}.</p>` : ""}
          ${s.rule_html ? `<div class="specimen-rule"><span class="label">Rule</span>${s.rule_html}</div>` : ""}
          <div class="specimen-tex"><code>${esc(s.tex)}</code><button type="button" class="btn" data-copytex="${esc(s.key)}">Copy TeX</button>${(s.macros || []).map(m => `<span class="chip" data-tip="A macro that gives the same symbol"><code>${esc(m)}</code></span>`).join("")}</div>
          <h3>Defined in</h3>
          ${note ? `<div class="specimen-note"><div class="chips">${typeChip(note)}${noteLink(note.id)}</div><div class="small">${note.summary_html || ""}</div></div>` : `<p class="small">No note defines it; the notation conventions (${esc(s.section || "")}) fix its meaning.</p>`}
          <h3>Reused in ${list.length ? `<span class="muted small">${list.length} ${list.length === 1 ? "note" : "notes"}</span>` : ""}</h3>
          ${list.length ? `<ul class="use-list">${list.map(u => `<li>${noteLink(u.note)}${u.tex ? `<div class="use-formula"><span class="math math-inline">${esc(u.tex)}</span></div>` : ""}</li>`).join("")}</ul>` : `<p class="small muted">${s.throughout ? "A short symbol that the conventions fix throughout the knowledge base." : "No other note uses it in a formula yet."}</p>`}
          ${s.note && s.tex.length ? '<p class="small muted specimen-foot">A symbol of fewer than four TeX tokens is searched only in the notes linked to its defining note.</p>' : ""}
        </div></div>`;
    };
    const tile = s => `<button type="button" class="glyph-tile${s.key === selected ? " on" : ""}" data-key="${esc(s.key)}" style="--sym-color:${colorOf(s)}" aria-pressed="${s.key === selected}">
        <span class="glyph${glyphSize(s)}"><span class="math math-inline">${esc(s.tex)}</span></span>
        <span class="glyph-cap">${s.object_html || s.meaning_html || ""}</span>
        ${(s.uses || []).length ? `<span class="glyph-uses" data-tip="Notes that reuse the symbol">${s.uses.length}</span>` : ""}</button>`;
    const row = s => {
      const list = s.uses || [];
      return `<tr style="--sym-color:${colorOf(s)}"><td class="sym-cell"><a class="sym-link" href="#/symbols/${enc(s.key)}"><span class="math math-inline">${esc(s.tex)}</span></a></td>
        <td>${s.meaning_html || ""}</td>
        <td class="small">${s.rule_html || (s.section ? esc(s.section) : '<span class="muted">None</span>')}</td>
        <td>${s.note ? noteLink(s.note) : '<span class="small muted">Convention</span>'}</td>
        <td class="small">${list.length ? list.slice(0, 5).map(u => noteLink(u.note)).join(", ") + (list.length > 5 ? ` and ${list.length - 5} more` : "") : `<span class="muted">${s.throughout ? "Throughout" : "None"}</span>`}</td>
        <td class="nowrap"><code class="small">${esc(s.tex)}</code> <button type="button" class="linkish small" data-copytex="${esc(s.key)}">Copy</button></td></tr>`;
    };
    const select = key => {
      selected = key;
      st.selected = key;
      store.set("symbols", st);
      $$(".glyph-tile").forEach(x => { const on = x.dataset.key === key; x.classList.toggle("on", on); x.setAttribute("aria-pressed", String(on)); });
      const panel = $("#sym-specimen");
      if (!panel) return;
      panel.innerHTML = specimen(symbolByKey.get(key));
      typeset(panel);
      bindCopyTex(panel);
      history.replaceState(null, "", "#/symbols/" + enc(key));
    };
    const draw = () => {
      store.set("symbols", st);
      $$("[data-symview]").forEach(b => b.classList.toggle("on", b.dataset.symview === st.view));
      const body = $("#sym-body");
      const list = SYMBOLS.filter(matches);
      if (!list.length) { body.innerHTML = '<p class="muted">No symbol matches.</p>'; return; }
      if (!list.some(s => s.key === selected)) selected = list[0].key;
      const groups = new Map();
      list.forEach(s => { const g = groupOf(s); if (!groups.has(g)) groups.set(g, []); groups.get(g).push(s); });
      const keys = orderGroups([...groups.keys()]);
      const heading = g => (g ? `<h2>${esc(g)} <span class="muted small">${groups.get(g).length}</span></h2>` : "");
      if (st.view === "table") {
        body.innerHTML = keys.map(g => `<section class="sym-group">${heading(g)}<div class="panel table-wrap"><table class="data-table sym-table"><thead><tr><th>Symbol</th><th>Meaning</th><th>Convention</th><th>Defined in</th><th>Reused in</th><th>TeX</th></tr></thead><tbody>${groups.get(g).map(row).join("")}</tbody></table></div></section>`).join("");
      } else {
        body.innerHTML = `<div class="atlas"><div class="atlas-grid">${keys.map(g => `<section class="sym-group">${heading(g)}<div class="glyph-grid">${groups.get(g).map(tile).join("")}</div></section>`).join("")}</div><aside class="atlas-specimen" id="sym-specimen" aria-live="polite">${specimen(symbolByKey.get(selected))}</aside></div>`;
        body.querySelectorAll(".glyph-tile").forEach(b => b.addEventListener("click", () => {
          select(b.dataset.key);
          if (window.matchMedia("(max-width: 1080px)").matches) $("#sym-specimen").scrollIntoView({ behavior: "smooth", block: "start" });
        }));
      }
      typeset(body);
      bindCopyTex(body);
    };
    $("#sym-q").addEventListener("input", debounce(ev => { st.q = ev.target.value; draw(); }, 120));
    $("#sym-group").addEventListener("change", ev => { st.group = ev.target.value; draw(); });
    $$("[data-symview]").forEach(b => b.addEventListener("click", () => { st.view = b.dataset.symview; draw(); }));
    $$(".filter-chip[data-family]").forEach(b => b.addEventListener("click", () => {
      const i = st.hiddenFamilies.indexOf(b.dataset.family);
      if (i >= 0) st.hiddenFamilies.splice(i, 1); else st.hiddenFamilies.push(b.dataset.family);
      b.classList.toggle("off", i < 0);
      draw();
    }));
    draw();
    typeset($(".macro-panel"));
    if (focusKey && symbolByKey.has(focusKey)) {
      const tileEl = main.querySelector(`.glyph-tile[data-key="${CSS.escape(focusKey)}"]`);
      if (tileEl) setTimeout(() => tileEl.scrollIntoView({ block: "center" }), 60);
    }
  }

  // ============================ Maps
  function showMaps() {
    const maps = NOTES.filter(n => n.type === "map");
    const cards = maps.map(m => {
      const members = (outE.get(m.id) || []).filter(e => e.k === "mentions").map(e => byId.get(e.t)).filter(Boolean);
      const byArea = members.reduce((acc, x) => { if (x.area) acc[x.area] = (acc[x.area] || 0) + 1; return acc; }, {});
      const total = Object.values(byArea).reduce((a, b) => a + b, 0) || 1;
      const stack = areaIds().filter(a => byArea[a]).map(a => `<span style="width:${100 * byArea[a] / total}%;background:${areaColor(a)}" data-tip="${esc(S.areas[a].label)}: ${byArea[a]}"></span>`).join("");
      return `<a class="card" href="#/n/${enc(m.id)}" style="--card-color:${m.area ? areaColor(m.area) : "var(--accent)"}"><div class="chips">${typeChip(m)}${maturityChip(m)}</div><h3>${m.title_html}</h3><div class="summary">${m.summary_html || ""}</div><div class="stack">${stack}</div><div class="foot"><span>${members.length} linked notes</span></div></a>`;
    }).join("");
    setView("maps", `<div class="page-head"><div><h1>Maps</h1><p>Curated overviews that order the notes of a domain and explain how they connect.</p></div></div>
      ${maps.length ? `<div class="grid">${cards}</div>` : '<div class="empty"><p>No maps yet. A map is a note of type map that links the notes of a domain in a reading order.</p></div>'}`, "Maps");
  }

  // ============================ Paths
  function showPaths(mode, a, b) {
    const sub = [["prereq", "Prerequisites"], ["connect", "Connection"], ["near", "Neighbourhood"]].map(([k, l]) => `<a class="btn${k === mode ? " btn-primary" : ""}" href="#/paths/${k}">${l}</a>`).join("");
    const nameOf = id => id && byId.has(id) ? plainTitle(byId.get(id).title) : "";
    const intro = {
      prereq: "Everything a note requires, directly or through other notes, in an order that never meets a note before its prerequisites.",
      connect: "The shortest chains of relations between two notes, with the direction in which each relation is written.",
      near: "The notes within one, two or three steps of a note.",
    }[mode];
    let picker = "";
    if (mode === "prereq") picker = `<input type="search" id="p-a" list="note-titles" placeholder="Choose a note" value="${esc(nameOf(a))}" aria-label="Note"><button class="btn btn-primary" type="button" id="p-go">Show prerequisites</button>`;
    if (mode === "connect") picker = `<input type="search" id="p-a" list="note-titles" placeholder="From" value="${esc(nameOf(a))}" aria-label="From"><input type="search" id="p-b" list="note-titles" placeholder="To" value="${esc(nameOf(b))}" aria-label="To"><label class="check-row"><input type="checkbox" id="p-mentions"${store.get("paths-mentions", false) ? " checked" : ""}>Follow mentions</label><button class="btn btn-primary" type="button" id="p-go">Connect</button>`;
    if (mode === "near") picker = `<input type="search" id="p-a" list="note-titles" placeholder="Choose a note" value="${esc(nameOf(a))}" aria-label="Note"><select id="p-depth" aria-label="Steps">${[1, 2, 3].map(d => `<option value="${d}"${String(b || 2) === String(d) ? " selected" : ""}>${d} ${d === 1 ? "step" : "steps"}</option>`).join("")}</select><button class="btn btn-primary" type="button" id="p-go">Show</button>`;
    setView("paths", `<div class="page-head"><div><h1>Paths</h1><p>${intro}</p></div><div class="subnav">${sub}</div></div>
      <div class="picker">${picker}</div><div class="paths-grid"><div id="p-list"></div><div class="graph-box" id="p-graph"></div></div>${datalist()}`, "Paths");
    const go = () => {
      const ida = findNote(($("#p-a") || {}).value);
      if (mode === "prereq" && ida) location.hash = `#/paths/prereq/${enc(ida)}`;
      if (mode === "near" && ida) location.hash = `#/paths/near/${enc(ida)}/${$("#p-depth").value}`;
      if (mode === "connect") {
        const idb = findNote($("#p-b").value);
        store.set("paths-mentions", $("#p-mentions").checked);
        if (ida && idb) location.hash = `#/paths/connect/${enc(ida)}/${enc(idb)}`;
      }
    };
    $("#p-go").addEventListener("click", go);
    $$(".picker input[type=search]").forEach(i => i.addEventListener("keydown", ev => { if (ev.key === "Enter") go(); }));
    const list = $("#p-list");
    const gbox = $("#p-graph");
    if (!a || !byId.has(a)) { list.innerHTML = '<p class="muted">Choose a note to begin.</p>'; gbox.hidden = true; return; }
    if (mode === "prereq") {
      const chain = (A.prerequisites || {})[a] || [];
      if (!chain.length) { list.innerHTML = `<p>No prerequisites are recorded for ${noteLink(a)}. Prerequisites come from the relation <em>requires</em>.</p>`; gbox.hidden = true; typeset(list); return; }
      let html = "";
      let current = -1;
      chain.forEach(([id, lv], i) => {
        if (lv !== current) { if (current >= 0) html += "</ol>"; html += `<div class="level-head">Level ${lv}</div><ol class="steps">`; current = lv; }
        const n = byId.get(id);
        html += `<li><span class="n">${i + 1}</span><span class="t">${noteLink(id)} ${maturityChip(n)}${understandingChip(n)}</span><span class="s">${n.summary_html || ""}</span></li>`;
      });
      html += `</ol><div class="level-head">Then</div><ol class="steps"><li><span class="n">${chain.length + 1}</span><span class="t">${noteLink(a)}</span><span class="s">${byId.get(a).summary_html || ""}</span></li></ol>`;
      list.innerHTML = html;
      typeset(list);
      const ids = new Set(chain.map(c => c[0]).concat([a]));
      makeNetwork(gbox, NOTES.filter(n => ids.has(n.id)), EDGES.filter(e => e.k === "requires" && ids.has(e.s) && ids.has(e.t)), { hierarchical: true, focus: a, compact: true, navigate: true });
    }
    if (mode === "connect") {
      if (!b || !byId.has(b)) { list.innerHTML = '<p class="muted">Choose the second note.</p>'; gbox.hidden = true; return; }
      const paths = shortestPaths(a, b, 3, store.get("paths-mentions", false));
      if (!paths.length) { list.innerHTML = `<p>No chain of ${store.get("paths-mentions", false) ? "relations or mentions" : "typed relations"} connects ${noteLink(a)} and ${noteLink(b)}.</p>`; gbox.hidden = true; typeset(list); return; }
      list.innerHTML = paths.map(p => `<div class="path-chain"><span class="path-node">${noteLink(a)}</span>${p.map(([, to, e, fwd]) => `<span class="path-edge">${fwd ? `${esc(edgeText(e))} →` : `← ${esc(edgeText(e))}`}</span><span class="path-node">${noteLink(to)}</span>`).join("")}</div>`).join("")
        + '<p class="small muted">An arrow to the right means the relation is written from the left note to the right note.</p>';
      typeset(list);
      const ids = new Set([a]);
      const used = new Set();
      paths.forEach(p => p.forEach(([from, to, e]) => { ids.add(from); ids.add(to); used.add(e); }));
      makeNetwork(gbox, NOTES.filter(n => ids.has(n.id)), [...used], { focus: a, compact: true, navigate: true });
    }
    if (mode === "near") {
      const depth = Math.max(1, Math.min(3, Number(b) || 2));
      const usable = EDGES.filter(e => isTyped(e) || e.k === "mentions");
      const dist = neighbourhood(a, usable, depth);
      let html = "";
      for (let d = 1; d <= depth; d++) {
        const ids = [...dist.entries()].filter(([, v]) => v === d).map(([k]) => k).filter(k => byId.has(k));
        html += `<div class="level-head">${d} ${d === 1 ? "step" : "steps"}: ${ids.length} notes</div><div class="chips">${ids.map(id => noteLink(id)).join(" · ") || '<span class="muted small">None</span>'}</div>`;
      }
      list.innerHTML = html;
      typeset(list);
      makeNetwork(gbox, NOTES.filter(n => dist.has(n.id)), usable.filter(e => dist.has(e.s) && dist.has(e.t)), { focus: a, compact: true, navigate: true });
    }
  }

  // ============================ References
  function showReferences(focusId) {
    const sources = NOTES.filter(n => n.type === "source");
    const st = Object.assign({ q: "", sort: "oldest", group: "decade", hiddenKinds: [], view: "list" }, store.get("references", {}));
    const citingIds = id => [...new Set((inE.get(id) || []).filter(e => e.k === "cites").map(e => e.s))];
    const yearOf = s => Number((s.bib || {}).year) || null;
    const years = sources.map(yearOf).filter(Boolean);
    const minYear = years.length ? Math.min(...years) : null;
    const maxYear = years.length ? Math.max(...years) : null;
    const citing = new Set(sources.flatMap(s => citingIds(s.id)));
    const sentences = sources.reduce((a, s) => a + (s.contexts || []).length, 0);
    const kinds = [...new Set(sources.map(s => s.kind || "article"))].sort();
    const kindLabel = k => String(k || "article").replace(/^./, c => c.toUpperCase());
    const familyOf = s => String((((s.bib || {}).authors || [])[0]) || s.label || s.id).split(",")[0];
    const areasPresent = areaIds().filter(a => sources.some(s => s.area === a));
    if (focusId && byId.has(focusId)) Object.assign(st, { view: "list", q: "", hiddenKinds: [] });
    setView("references", `<div class="page-head"><div><h1>References</h1><p>${sources.length} works cited across ${citing.size} notes${years.length ? `, published from ${minYear} to ${maxYear}` : ""}. Each work shows the notes that cite it and the sentences in which they cite it.</p></div>
        <div class="subnav"><a class="btn" href="references.bib" download="references.bib">BibTeX file</a></div></div>
      <div class="kpis">${kpiHtml(sources.length, "Works")}${kpiHtml(citing.size, "Citing notes")}${kpiHtml(sentences, "Citing sentences")}${kpiHtml(years.length ? maxYear - minYear : 0, "Years spanned")}</div>
      ${years.length ? `<section class="panel timeline-panel"><h2>Timeline</h2><div class="timeline" id="ref-timeline"></div><div class="legend-inline">${areasPresent.map(a => `<span><span class="swatch" style="background:${areaColor(a)}"></span>${esc(S.areas[a].label)}</span>`).join("")}<span>The area of a circle grows with the number of citing notes; click a circle to find the work.</span></div></section>` : ""}
      <div class="toolbar"><input type="search" class="grow" id="ref-q" placeholder="Filter by author, title, venue or citing note" value="${esc(st.q)}" aria-label="Filter references">
        <select id="ref-sort" aria-label="Sort references">${[["oldest", "Oldest first"], ["newest", "Newest first"], ["author", "Author"], ["cited", "Most cited"]].map(([v, l]) => `<option value="${v}"${st.sort === v ? " selected" : ""}>${l}</option>`).join("")}</select>
        <select id="ref-group" aria-label="Group references">${[["decade", "Group by decade"], ["kind", "Group by kind"], ["area", "Group by area"], ["none", "No groups"]].map(([v, l]) => `<option value="${v}"${st.group === v ? " selected" : ""}>${l}</option>`).join("")}</select>
        <div class="seg" role="group">${[["list", "List"], ["map", "Citation map"]].map(([v, l]) => `<button type="button" data-view="${v}" class="${st.view === v ? "on" : ""}">${l}</button>`).join("")}</div></div>
      <div class="toolbar">${kinds.map(k => `<button type="button" class="filter-chip${st.hiddenKinds.includes(k) ? " off" : ""}" data-kind="${esc(k)}">${esc(kindLabel(k))}</button>`).join("")}</div>
      <div id="ref-body"></div>`, "References");
    const focusEntry = id => {
      const el = document.getElementById("ref-" + id);
      if (!el) return;
      el.scrollIntoView({ behavior: "smooth", block: "center" });
      el.classList.add("flash");
      setTimeout(() => el.classList.remove("flash"), 1800);
    };
    const entryHtml = s => {
      const notes = citingIds(s.id);
      const locators = id => [...new Set(((byId.get(id) || {}).citations || []).filter(c => c.key === s.id && c.locator).map(c => c.locator))];
      const ctx = s.contexts || [];
      const doi = (s.bib || {}).doi;
      return `<article class="ref-card" id="ref-${esc(s.id)}" style="--card-color:${s.area ? areaColor(s.area) : "var(--line-strong)"}">
        <div class="ref-side"><div class="ref-year">${esc(yearOf(s) || "n.d.")}</div><div class="ref-author">${esc(familyOf(s))}</div><div class="chips"><span class="chip">${esc(kindLabel(s.kind))}</span>${verificationChip(s)}</div></div>
        <div class="ref-main"><div class="reference">${s.reference_html || ""}</div>
          <div class="ref-cited">${notes.length ? `<span class="muted small">Cited in</span> ${notes.map(id => `<span class="ref-note">${noteLink(id)}${locators(id).length ? ` <span class="muted small">${locators(id).map(esc).join("; ")}</span>` : ""}</span>`).join("")}` : '<span class="muted small">No note cites this work yet.</span>'}</div>
          ${ctx.length ? `<details class="ref-contexts"><summary>${ctx.length} citing ${ctx.length === 1 ? "sentence" : "sentences"}</summary><ul class="ctx-list">${ctx.map(c => `<li><div class="small muted">${noteLink(c.note)}${c.locator ? `, ${esc(c.locator)}` : ""}</div><div class="ctx">${c.html}</div></li>`).join("")}</ul></details>` : ""}
          <div class="ref-actions"><a class="btn" href="#/n/${enc(s.id)}">Open note</a><button type="button" class="btn" data-copy="${esc(s.id)}">Copy BibTeX</button>${doi ? `<a class="btn" href="https://doi.org/${esc(doi)}" target="_blank" rel="noopener">DOI</a>` : ""}</div></div></article>`;
    };
    const draw = () => {
      store.set("references", st);
      $$("[data-view]").forEach(b => b.classList.toggle("on", b.dataset.view === st.view));
      const body = $("#ref-body");
      destroyNets();
      const q = st.q.trim().toLowerCase();
      let list = sources.filter(s => !st.hiddenKinds.includes(s.kind || "article"));
      if (q) list = list.filter(s => (stripTags(s.reference_html) + " " + citingIds(s.id).map(id => plainTitle((byId.get(id) || {}).title || id)).join(" ")).toLowerCase().includes(q));
      if (st.view === "map") {
        body.innerHTML = list.length ? '<div class="graph-box ref-map" id="ref-map"></div><p class="small muted">Boxes are works and the other shapes are notes; an arrow runs from a note to each work that it cites. Click a node to open it.</p>' : '<p class="muted">No reference matches.</p>';
        if (list.length) {
          const ids = new Set(list.map(s => s.id));
          list.forEach(s => citingIds(s.id).forEach(id => ids.add(id)));
          const works = new Set(list.map(s => s.id));
          makeNetwork($("#ref-map"), NOTES.filter(n => ids.has(n.id)), EDGES.filter(e => e.k === "cites" && ids.has(e.s) && works.has(e.t)), { navigate: true });
        }
        return;
      }
      const sorters = {
        oldest: (a, b) => (yearOf(a) || 9999) - (yearOf(b) || 9999) || familyOf(a).localeCompare(familyOf(b)),
        newest: (a, b) => (yearOf(b) || 0) - (yearOf(a) || 0) || familyOf(a).localeCompare(familyOf(b)),
        author: (a, b) => familyOf(a).localeCompare(familyOf(b)) || (yearOf(a) || 0) - (yearOf(b) || 0),
        cited: (a, b) => citingIds(b.id).length - citingIds(a.id).length || familyOf(a).localeCompare(familyOf(b)),
      };
      list.sort(sorters[st.sort] || sorters.oldest);
      const groupOf = s => st.group === "decade" ? (yearOf(s) ? `${Math.floor(yearOf(s) / 10) * 10}s` : "Undated")
        : st.group === "kind" ? kindLabel(s.kind) : st.group === "area" ? (s.area ? S.areas[s.area].label : "No area") : "";
      const groups = new Map();
      list.forEach(s => { const g = groupOf(s); if (!groups.has(g)) groups.set(g, []); groups.get(g).push(s); });
      body.innerHTML = list.length ? [...groups.entries()].map(([g, items]) => `<section class="ref-group">${g ? `<h2>${esc(g)} <span class="muted small">${items.length}</span></h2>` : ""}${items.map(entryHtml).join("")}</section>`).join("") : '<p class="muted">No reference matches.</p>';
      typeset(body);
      bindCopy(body);
    };
    const drawTimeline = () => {
      const box = $("#ref-timeline");
      if (!box) return;
      const lo = Math.floor(minYear / 10) * 10;
      const hi = Math.max(lo + 10, Math.ceil((maxYear + 1) / 10) * 10);
      const W = 1000;
      const base = 150;
      const pad = 40;
      const x = y => pad + (W - 2 * pad) * (y - lo) / (hi - lo);
      const placed = [];
      const dots = sources.filter(yearOf).sort((a, b) => yearOf(a) - yearOf(b) || familyOf(a).localeCompare(familyOf(b))).map(s => {
        const y = yearOf(s);
        const px = x(y);
        let level = 0;
        while (placed.some(p => p.level === level && Math.abs(p.x - px) < 82)) level += 1;
        placed.push({ x: px, level });
        const n = citingIds(s.id).length;
        const r = 6 + 2.4 * Math.sqrt(n);
        const cy = base - 26 - level * 34;
        const color = s.area ? areaColor(s.area) : cssVar("--neutral-node");
        return `<g class="tl-dot" data-ref="${esc(s.id)}" tabindex="0" role="link" aria-label="${esc(s.label)}"><line x1="${px}" y1="${base}" x2="${px}" y2="${cy + r}" class="tl-stem"/><circle cx="${px}" cy="${cy}" r="${r}" fill="${color}" data-tip="${esc(s.label)}: cited in ${n} ${n === 1 ? "note" : "notes"}"/><text x="${px}" y="${cy - r - 5}" text-anchor="middle" class="tl-label">${esc(familyOf(s))} ${y}</text></g>`;
      }).join("");
      const ticks = [];
      for (let y = lo; y <= hi; y += 10) ticks.push(`<g class="tl-tick"><line x1="${x(y)}" y1="${base}" x2="${x(y)}" y2="${base + 6}"/><text x="${x(y)}" y="${base + 22}" text-anchor="middle">${y}</text></g>`);
      const levels = placed.reduce((m, p) => Math.max(m, p.level), 0);
      const top = Math.min(0, base - 26 - levels * 34 - 36);
      box.innerHTML = `<svg viewBox="0 ${top} ${W} ${base + 32 - top}" role="img" aria-label="Publication years of the works cited"><line x1="${pad}" y1="${base}" x2="${W - pad}" y2="${base}" class="tl-axis"/>${ticks.join("")}${dots}</svg>`;
      box.querySelectorAll("[data-ref]").forEach(g => {
        const go = () => { if (st.view !== "list") { st.view = "list"; draw(); } focusEntry(g.dataset.ref); };
        g.addEventListener("click", go);
        g.addEventListener("keydown", ev => { if (ev.key === "Enter") go(); });
      });
    };
    $("#ref-q").addEventListener("input", debounce(ev => { st.q = ev.target.value; draw(); }, 120));
    $("#ref-sort").addEventListener("change", ev => { st.sort = ev.target.value; draw(); });
    $("#ref-group").addEventListener("change", ev => { st.group = ev.target.value; draw(); });
    $$("[data-view]").forEach(b => b.addEventListener("click", () => { st.view = b.dataset.view; draw(); }));
    $$(".filter-chip[data-kind]").forEach(b => b.addEventListener("click", () => {
      const i = st.hiddenKinds.indexOf(b.dataset.kind);
      if (i >= 0) st.hiddenKinds.splice(i, 1); else st.hiddenKinds.push(b.dataset.kind);
      b.classList.toggle("off", i < 0);
      draw();
    }));
    drawTimeline();
    draw();
    if (focusId && byId.has(focusId)) setTimeout(() => focusEntry(focusId), 80);
  }

  // ============================ Philosophy
  function showPhilosophy(focusTheme) {
    const dims = PHILO.dimensions || {};
    const themes = PHILO.themes || {};
    const concepts = NOTES.filter(n => CONCEPT_TYPES.includes(n.type));
    const analysed = concepts.filter(n => (n.toc || []).some(t => t.id === "ontology-and-epistemology"));
    const byTitle = (a, b) => plainTitle(a.title).localeCompare(plainTitle(b.title));
    const notesFor = id => concepts.filter(n => (n.themes || []).includes(id)).sort(byTitle);
    const dimIds = Object.keys(dims);
    const themesOf = dim => Object.keys(themes).filter(k => themes[k].dimension === dim);
    const allThemes = dimIds.flatMap(themesOf);
    const focus = focusTheme && themes[focusTheme] ? focusTheme : null;
    const sectionHtml = n => {
      const tmp = document.createElement("div");
      tmp.innerHTML = n.html || "";
      const h = tmp.querySelector("#ontology-and-epistemology");
      if (!h) return "";
      const parts = [];
      for (let el = h.nextElementSibling; el && el.tagName !== "H2"; el = el.nextElementSibling) parts.push(el.outerHTML);
      return parts.join("");
    };
    const card = id => {
      const t = themes[id];
      const list = notesFor(id);
      return `<article class="theme-card dim-${esc(t.dimension)}${focus === id ? " focus" : ""}" id="theme-${esc(id)}">
        <header><span class="dim-mark">${DIM_MARK[t.dimension] || ""}</span><h3><a href="#/philosophy/${enc(id)}">${esc(t.label)}</a></h3><span class="count" data-tip="Notes that raise this theme">${list.length}</span></header>
        <p class="theme-q">${esc(t.question)}</p><p class="small">${esc(t.definition)}</p>
        <div class="theme-notes">${list.map(n => noteLink(n.id)).join('<span class="sep"> · </span>') || '<span class="small muted">No note raises this theme yet.</span>'}</div></article>`;
    };
    const focusHtml = focus ? `<section class="panel focus-panel dim-${esc(themes[focus].dimension)}" id="theme-focus"><h2><span><span class="dim-mark">${DIM_MARK[themes[focus].dimension] || ""}</span> ${esc(themes[focus].label)}</span><a class="small" href="#/philosophy">All themes</a></h2><p class="theme-q">${esc(themes[focus].question)}</p>
        ${notesFor(focus).map(n => `<article class="excerpt"><h3>${noteLink(n.id)}<span class="chips inline">${(n.themes || []).filter(x => x !== focus).map(themeChip).join("")}</span></h3><div class="prose">${sectionHtml(n)}</div></article>`).join("") || '<p class="small muted">No note raises this theme yet.</p>'}</section>` : "";
    const rows = concepts.filter(n => (n.themes || []).length).sort(byTitle);
    const matrix = rows.length ? `<div class="table-wrap"><table class="data-table theme-matrix"><thead><tr><th>Note</th>${allThemes.map(id => `<th class="rot dim-${esc(themes[id].dimension)}"><a href="#/philosophy/${enc(id)}" data-tip="${esc(themes[id].question)}">${esc(themes[id].label)}</a></th>`).join("")}</tr></thead><tbody>${rows.map(n => `<tr><td>${noteLink(n.id)}</td>${allThemes.map(id => `<td class="cell dim-${esc(themes[id].dimension)}">${(n.themes || []).includes(id) ? `<span class="dot" data-tip="${esc(themes[id].label)}"></span>` : ""}</td>`).join("")}</tr>`).join("")}</tbody></table></div>` : '<p class="small muted">No note lists a theme yet.</p>';
    setView("philosophy", `<div class="page-head"><div><h1>Philosophy</h1><p>What the concepts of this knowledge base take to exist, and how claims about them can be known. Each developing note lists its philosophical themes and argues them in a section on ontology and epistemology.</p></div></div>
      <div class="kpis">${kpiHtml(`${analysed.length}/${concepts.length}`, "Notes with an analysis")}${dimIds.map(d => kpiHtml(themesOf(d).length, `${esc(dims[d].label)} themes`)).join("")}${kpiHtml(concepts.reduce((a, n) => a + (n.themes || []).length, 0), "Themes raised in notes")}</div>
      ${focusHtml}
      <div class="dims">${dimIds.map(d => `<section class="dim dim-${esc(d)}"><div class="dim-head"><span class="dim-mark big">${DIM_MARK[d] || ""}</span><div><h2>${esc(dims[d].label)}</h2><p class="theme-q">${esc(dims[d].question)}</p><p class="small muted">${esc(dims[d].description)}</p></div></div><div class="theme-grid">${themesOf(d).map(card).join("")}</div></section>`).join("")}</div>
      <section class="panel"><h2>Themes by note</h2>${matrix}</section>`, focus ? `${themes[focus].label} · Philosophy` : "Philosophy");
  }


  // ============================ Timeline
  const timelineDefaults = () => ({ order: "new", area: null, sessions: true });
  const timelineState = () => Object.assign(timelineDefaults(), store.get("timeline", {}));

  function timelineDays(area) {
    const byDate = new Map();
    const add = date => { if (!byDate.has(date)) byDate.set(date, { date, notes: [], sessions: [] }); return byDate.get(date); };
    NOTES.forEach(n => {
      if (!n.created) return;
      if (area && n.area !== area) return;
      add(n.created).notes.push(n);
    });
    ((DATA.timeline || {}).sessions || []).forEach(s => add(s.date).sessions.push(s));
    return [...byDate.values()].sort((a, b) => (a.date < b.date ? -1 : a.date > b.date ? 1 : 0));
  }

  function timelineChart(days) {
    if (!days.length) return "";
    const W = 900, H = 200, padL = 34, padR = 12, padT = 10, padB = 26;
    const cols = days.length;
    const step = (W - padL - padR) / Math.max(cols, 1);
    const cell = Math.max(3, Math.min(11, step * 0.5));
    const tallest = Math.max(1, ...days.map(d => d.notes.length));
    const perColumn = Math.max(1, Math.floor((H - padT - padB) / (cell + 2)));
    const scale = tallest > perColumn ? perColumn / tallest : 1;   // squeeze a tall day into the panel
    let total = 0;
    const cumulative = days.map(d => { total += d.notes.length; return total; });
    const maxTotal = Math.max(1, total);
    const squares = days.map((d, i) => {
      const x = padL + i * step + (step - cell) / 2;
      return d.notes.map((n, k) => {
        const y = H - padB - (k + 1) * (cell + 2) * scale;
        const colour = n.area ? areaColor(n.area) : cssVar("--neutral-node");
        return `<a href="#/n/${enc(n.id)}"><rect class="tl-unit" x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${cell.toFixed(1)}" height="${(cell * scale).toFixed(1)}" rx="1.5" fill="${colour}"><title>${esc(plainTitle(n.title))} (${esc(levelLabel("type", n.type) || n.type)}, ${esc(d.date)})</title></rect></a>`;
      }).join("");
    }).join("");
    const line = days.map((d, i) => {
      const x = padL + i * step + step / 2;
      const y = H - padB - (cumulative[i] / maxTotal) * (H - padT - padB);
      return `${i ? "L" : "M"}${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(" ");
    const labels = days.map((d, i) => {
      if (cols > 12 && i % Math.ceil(cols / 12) !== 0 && i !== cols - 1) return "";
      const x = padL + i * step + step / 2;
      return `<text class="tl-axis" x="${x.toFixed(1)}" y="${H - 8}" text-anchor="middle">${esc(d.date.slice(5))}</text>`;
    }).join("");
    return `<figure class="tl-chart"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Notes added on each day, one square per note, with the cumulative total as a line">
      <line class="tl-base" x1="${padL - 6}" y1="${H - padB}" x2="${W - padR}" y2="${H - padB}"/>
      ${squares}
      <path class="tl-cumulative" d="${line}"/>
      <text class="tl-axis" x="${padL - 10}" y="${padT + 8}" text-anchor="end">${total}</text>
      <text class="tl-axis" x="${padL - 10}" y="${H - padB}" text-anchor="end">0</text>
      ${labels}
    </svg><figcaption>One square for each note, coloured by area and stacked over the day it was created; the line is the running total, which reaches ${total}.</figcaption></figure>`;
  }

  function timelineDayCard(day, showSessions) {
    const groups = new Map();
    day.notes.forEach(n => { if (!groups.has(n.type)) groups.set(n.type, []); groups.get(n.type).push(n); });
    const order = ["concept", "theorem", "method", "model", "example", "map", "question", "source"];
    const sorted = [...groups.entries()].sort((a, b) => order.indexOf(a[0]) - order.indexOf(b[0]));
    const chips = sorted.map(([type, list]) => `<div class="tl-group"><span class="tl-kind">${esc(levelLabel("type", type) || type)}</span>${list.map(n => `<a class="tl-note" href="#/n/${enc(n.id)}">${shapeIcon(n.type, n.area ? areaColor(n.area) : null)}<span>${n.title_html}</span></a>`).join("")}</div>`).join("");
    const areas = [...new Set(day.notes.map(n => n.area).filter(Boolean))];
    const sessions = showSessions ? day.sessions.map(s => `<details class="tl-session"><summary>${esc(s.title)}</summary><div class="tl-session-body">${s.html}</div></details>`).join("") : "";
    const heads = day.sessions.map(s => esc(s.title)).join(" · ");
    return `<li class="tl-day"><div class="tl-dot" aria-hidden="true"></div>
      <div class="tl-card">
        <div class="tl-head"><h2>${esc(day.date)}</h2><span class="chip">${day.notes.length} note${day.notes.length === 1 ? "" : "s"}</span>${day.sessions.length ? `<span class="chip">${day.sessions.length} session${day.sessions.length === 1 ? "" : "s"}</span>` : ""}${areas.map(a => `<span class="chip" style="border-color:${areaColor(a)}">${esc(((S.areas || {})[a] || {}).label || a)}</span>`).join("")}</div>
        ${heads && !showSessions ? `<p class="small muted">${heads}</p>` : ""}
        ${sessions}
        ${chips || '<p class="small muted">No new notes on this day.</p>'}
      </div></li>`;
  }

  function showTimeline() {
    const st = timelineState();
    const tl = DATA.timeline || {};
    const days = timelineDays(st.area);
    const shown = st.order === "new" ? [...days].reverse() : days;
    const areas = [...new Set(NOTES.map(n => n.area).filter(Boolean))].sort();
    const notes = NOTES.filter(n => n.type !== "source").length;
    const sources = NOTES.filter(n => n.type === "source").length;
    const controls = `<div class="tl-controls">
      <div class="seg" role="tablist" aria-label="Order">
        <button class="${st.order === "new" ? "on" : ""}" data-order="new">Newest first</button>
        <button class="${st.order === "old" ? "on" : ""}" data-order="old">As it happened</button>
      </div>
      <div class="chips"><button class="chip${st.area ? "" : " on"}" data-area="">All areas</button>${areas.map(a => `<button class="chip${st.area === a ? " on" : ""}" data-area="${esc(a)}" style="border-color:${areaColor(a)}">${esc(((S.areas || {})[a] || {}).label || a)}</button>`).join("")}</div>
      ${(tl.sessions || []).length ? `<label class="tl-toggle"><input type="checkbox" id="tl-sessions"${st.sessions ? " checked" : ""}> Session notes</label>` : ""}
    </div>`;
    const span = tl.first && tl.last ? `from ${esc(tl.first)} to ${esc(tl.last)}` : "";
    const privacy = (tl.sessions || []).length
      ? `<p class="small muted">The session notes are the entries of the changelog, and they appear on this page only in the local build.</p>`
      : `<p class="small muted">The public build shows what was added and when, and leaves out the working notes of each session.</p>`;
    setView("timeline", `<div class="page-head"><div><h1>Timeline</h1>
        <p>${notes} notes and ${sources} sources, added over ${(tl.days || days.length)} working days ${span}. Every square is one note, placed on the day it was written; follow the spine downwards to read the order in which the ideas arrived.</p>
        ${privacy}</div></div>
      ${controls}
      ${timelineChart(days)}
      <ol class="tl-spine">${shown.map(d => timelineDayCard(d, st.sessions)).join("")}</ol>`, "Timeline");
    main.querySelectorAll("[data-order]").forEach(b => b.addEventListener("click", () => { store.set("timeline", Object.assign(st, { order: b.dataset.order })); showTimeline(); }));
    main.querySelectorAll("[data-area]").forEach(b => b.addEventListener("click", () => { store.set("timeline", Object.assign(st, { area: b.dataset.area || null })); showTimeline(); }));
    const toggle = $("#tl-sessions");
    if (toggle) toggle.addEventListener("change", () => { store.set("timeline", Object.assign(st, { sessions: toggle.checked })); showTimeline(); });
  }

  // ============================ Questions
  function showQuestions() {
    const qs = NOTES.filter(n => n.type === "question");
    const column = (status, title) => {
      const list = qs.filter(q => q.status === status);
      return `<div class="column"><h2>${esc(title)} <span class="chip">${list.length}</span></h2>${list.map(q => `<a class="qcard" href="#/n/${enc(q.id)}"><h3>${q.title_html}</h3><p>${q.summary_html || ""}</p><div class="chips">${(q.concerns || []).slice(0, 4).map(c => `<span class="chip">${esc(plainTitle((byId.get(c) || {}).title || c))}</span>`).join("")}</div></a>`).join("") || '<p class="small muted" style="margin:6px">None.</p>'}</div>`;
    };
    setView("questions", `<div class="page-head"><div><h1>Questions</h1><p>Open conceptual questions mark where the knowledge base is still incomplete. A question is resolved when the notes that answer it exist.</p></div></div>
      <div class="board">${column("open", levelLabel("status", "open"))}${column("in-progress", levelLabel("status", "in-progress"))}${column("resolved", levelLabel("status", "resolved"))}</div>`, "Questions");
  }

  // ============================ State
  function showState() {
    const concepts = NOTES.filter(n => CONCEPT_TYPES.includes(n.type));
    const sources = NOTES.filter(n => n.type === "source");
    const typedCount = EDGES.filter(isTyped).length;
    const checks = concepts.flatMap(n => n.checks || []);
    const uniqueChecks = [...new Map(checks.map(c => [c.id, c])).values()];
    const passing = uniqueChecks.filter(c => c.status === "pass" && !c.stale).length;
    const checkedStatus = new Set(["verified", "manual"]);
    const verified = sources.filter(s => checkedStatus.has((s.verification || {}).status)).length;
    const kpi = (v, k) => `<div class="kpi"><div class="v">${v}</div><div class="k">${k}</div></div>`;
    const maxType = Math.max(1, ...CONCEPT_TYPES.map(t => concepts.filter(n => n.type === t).length));
    const matColors = { seed: "var(--mat-seed)", developing: "var(--mat-developing)", mature: "var(--mat-mature)" };
    const maturityRows = CONCEPT_TYPES.filter(t => concepts.some(n => n.type === t)).map(t => {
      const list = concepts.filter(n => n.type === t);
      const segs = levelIds("maturity").map(m => { const c = list.filter(n => n.maturity === m).length; return c ? `<span style="width:${100 * c / maxType}%;background:${matColors[m]}" data-tip="${esc(typePlural(t))}, ${esc(levelLabel("maturity", m).toLowerCase())}: ${c}"></span>` : ""; }).join("");
      return `<span>${esc(typePlural(t))}</span><div class="stack">${segs}</div><span class="muted">${list.length}</span>`;
    }).join("");
    const confCounts = levelIdsDesc("confidence").map(k => ({ k, c: concepts.filter(n => n.confidence === k).length }));
    const maxConf = Math.max(1, ...confCounts.map(x => x.c));
    const confRows = confCounts.map(({ k, c }) => `<span>${esc(levelLabel("confidence", k))}</span><div class="bar-track"><div class="bar" style="width:${100 * c / maxConf}%;background:${level("confidence", k).background};outline:1px solid ${level("confidence", k).foreground}" data-tip="${esc(levelLabel("confidence", k))} confidence: ${c}"></div></div><span class="muted">${c}</span>`).join("");
    const domCounts = {};
    concepts.forEach(n => (n.domains || []).forEach(d => { domCounts[d] = (domCounts[d] || 0) + 1; }));
    const maxDom = Math.max(1, ...Object.values(domCounts));
    const domRows = Object.entries(domCounts).sort((x, y) => y[1] - x[1]).map(([d, c]) => { const area = ((S.domains || {})[d] || {}).area; return `<span>${esc(domainLabel(d))}</span><div class="bar-track"><div class="bar" style="width:${100 * c / maxDom}%;background:${areaColor(area)}" data-tip="${esc(domainLabel(d))}: ${c} notes"></div></div><span class="muted">${c}</span>`; }).join("");
    const ranked = (items, fmt) => items.length ? `<ol class="rank-list">${items.map(fmt).join("")}</ol>` : '<p class="small muted">None.</p>';
    const develop = (A.develop_next || []).filter(d => byId.has(d.id)).slice(0, 10);
    const hubs = Object.entries(A.degree || {}).filter(([id]) => byId.has(id) && CONCEPT_TYPES.includes(byId.get(id).type)).sort((x, y) => y[1] - x[1]).slice(0, 8);
    const bridges = Object.entries(A.betweenness || {}).filter(([id, v]) => byId.has(id) && v > 0).sort((x, y) => y[1] - x[1]).slice(0, 8);
    const bottlenecks = Object.entries(A.dependents || {}).filter(([, v]) => v > 0).sort((x, y) => y[1] - x[1]).slice(0, 8);
    const months = {};
    NOTES.forEach(n => { const m = String(n.updated || "").slice(0, 7); if (m) months[m] = (months[m] || 0) + 1; });
    const monthKeys = Object.keys(months).sort().slice(-18);
    const maxMonth = Math.max(1, ...monthKeys.map(m => months[m]));
    const activity = monthKeys.length ? `<div class="columns">${monthKeys.map(m => `<span style="height:${Math.max(4, 100 * months[m] / maxMonth)}%" data-tip="${m}: ${months[m]} notes updated"></span>`).join("")}</div><div class="legend-inline"><span>${esc(monthKeys[0])}</span><span style="margin-left:auto">${esc(monthKeys[monthKeys.length - 1])}</span></div>` : '<p class="small muted">No dates yet.</p>';
    const v = A.validation || {};
    const findings = (A.findings || []);
    const byCode = {};
    findings.forEach(f => { (byCode[f.code] = byCode[f.code] || []).push(f); });
    const findingsHtml = LOCAL ? (findings.length ? Object.entries(byCode).map(([code, list]) => `<details><summary>${esc(code)} ${esc((S.codes || {})[code] || "")} (${list.length})</summary><ul>${list.slice(0, 60).map(f => `<li><code>${esc(f.path)}${f.line ? ":" + f.line : ""}</code> ${esc(f.message)}</li>`).join("")}</ul></details>`).join("") : '<p class="small">No findings.</p>') : "";
    const checkRows = uniqueChecks.map(c => `<tr><td><code>${esc(c.id)}</code></td><td class="${c.status === "pass" && !c.stale ? "chip-good" : "chip-bad"}">${c.stale ? "! Stale" : ({ pass: "✓ Pass", fail: "✕ Fail", error: "✕ Error", "not-run": "! Not run" }[c.status] || esc(c.status))}</td><td class="small">${esc(c.summary || "")}</td><td class="small muted">${esc(String(c.date || "").slice(0, 10))}</td></tr>`).join("");
    const unverified = sources.filter(s => !checkedStatus.has((s.verification || {}).status));
    const analysed = concepts.filter(n => (n.toc || []).some(t => t.id === "ontology-and-epistemology")).length;
    const themeCounts = Object.keys(PHILO.themes || {}).map(id => ({ id, t: PHILO.themes[id], c: concepts.filter(n => (n.themes || []).includes(id)).length })).filter(x => x.c).sort((a, b) => b.c - a.c);
    const maxTheme = Math.max(1, ...themeCounts.map(x => x.c));
    const themeRows = themeCounts.map(({ id, t, c }) => `<a href="#/philosophy/${enc(id)}"><span class="dim-mark dim-${esc(t.dimension)}">${DIM_MARK[t.dimension] || ""}</span> ${esc(t.label)}</a><div class="bar-track"><div class="bar dim-bar dim-${esc(t.dimension)}" style="width:${100 * c / maxTheme}%" data-tip="${esc(t.label)}: ${c} notes"></div></div><span class="muted">${c}</span>`).join("");
    const html = `<div class="page-head"><div><h1>Knowledge state</h1><p>How much is written, how mature and well supported it is, and where to work next.</p></div></div>
      <div class="kpis">
        ${kpi(concepts.length, "Notes")}${kpi(NOTES.filter(n => n.type === "map").length, "Maps")}
        ${kpi(NOTES.filter(n => n.type === "question" && n.status !== "resolved").length, "Open questions")}
        ${kpi(`${verified}/${sources.length}`, "References checked")}${kpi(typedCount, "Typed relations")}
        ${kpi(`${passing}/${uniqueChecks.length}`, "Checks passing")}${kpi(`${analysed}/${concepts.length}`, "Philosophical analyses")}
      </div>
      <div class="dash">
        <section class="panel"><h2>Maturity by type</h2>${maturityRows ? `<div class="bars">${maturityRows}</div><div class="legend-inline">${levelIds("maturity").map(m => `<span><span class="swatch" style="background:${matColors[m]}"></span>${esc(levelLabel("maturity", m))}</span>`).join("")}</div>` : '<p class="small muted">No notes yet.</p>'}</section>
        <section class="panel"><h2>Confidence</h2><div class="bars">${confRows}</div></section>
        <section class="panel"><h2>Philosophical themes</h2>${themeRows ? `<div class="bars">${themeRows}</div><div class="legend-inline"><span><span class="dim-mark dim-ontology">${DIM_MARK.ontology}</span> Ontology</span><span><span class="dim-mark dim-epistemology">${DIM_MARK.epistemology}</span> Epistemology</span></div>` : '<p class="small muted">No note lists a theme yet.</p>'}</section>
        <section class="panel"><h2>Notes by domain</h2>${domRows ? `<div class="bars">${domRows}</div>` : '<p class="small muted">No notes yet.</p>'}</section>
        <section class="panel"><h2>Develop next</h2><p class="small muted">Seed and developing notes ranked by how many notes relate to them, mention them or require them.</p>${ranked(develop, d => `<li>${noteLink(d.id)} ${maturityChip(byId.get(d.id))}</li>`)}</section>
        <section class="panel"><h2>Most connected</h2>${ranked(hubs, ([id, c]) => `<li>${noteLink(id)} <span class="muted small">${c} links</span></li>`)}</section>
        <section class="panel"><h2>Bridges</h2><p class="small muted">Notes that lie on many shortest paths between other notes (normalised betweenness).</p>${ranked(bridges, ([id, b]) => `<li>${noteLink(id)} <span class="muted small">${Number(b).toFixed(3)}</span></li>`)}</section>
        <section class="panel"><h2>Prerequisite bottlenecks</h2><p class="small muted">Notes that the largest number of other notes require, directly or indirectly.</p>${ranked(bottlenecks, ([id, c]) => `<li>${noteLink(id)} <span class="muted small">${c} dependent notes</span></li>`)}</section>
        <section class="panel"><h2>Wanted notes</h2>${ranked((A.wanted || []).slice(0, 12), w => `<li><a class="wikilink wanted" href="#/wanted/${enc(w.id)}">${esc(w.id)}</a> <span class="muted small">linked from ${w.from.length}</span></li>`)}</section>
        <section class="panel"><h2>Isolated notes</h2>${ranked((A.isolated || []).filter(id => byId.has(id)), id => `<li>${noteLink(id)}</li>`)}</section>
        <section class="panel"><h2>References to check</h2>${ranked(unverified, s => `<li>${noteLink(s.id, esc(s.label || s.id))} ${verificationChip(s)}</li>`)}</section>
        <section class="panel"><h2>Activity</h2><p class="small muted">Notes by the month of their last update.</p>${activity}</section>
        <section class="panel"><h2>Numerical checks</h2>${uniqueChecks.length ? `<div class="table-wrap"><table class="data-table"><tbody>${checkRows}</tbody></table></div>` : '<p class="small muted">No checks yet.</p>'}</section>
        ${LOCAL ? `<section class="panel findings"><h2>Validation</h2><p class="small">${v.errors || 0} errors and ${v.warnings || 0} warnings at build time.</p>${findingsHtml}</section>` : ""}
      </div>`;
    setView("state", html, "Knowledge state");
  }

  // ============================ Recall
  function showRecall() {
    const st = Object.assign({ area: "", type: "", missedOnly: false }, store.get("recall-filter", {}));
    const record = store.get("recall", {});
    const pool = () => NOTES.filter(n => CONCEPT_TYPES.includes(n.type) && (!st.area || n.area === st.area) && (!st.type || n.type === st.type) && (!st.missedOnly || ((record[n.id] || {}).m || 0) > ((record[n.id] || {}).k || 0)));
    setView("recall", `<div class="recall"><div class="page-head"><div><h1>Recall</h1><p>Test yourself on the notes: read the title, recall the definition, then reveal it. Your answers stay in this browser.</p></div></div>
      <div class="toolbar"><select id="rc-area" aria-label="Area"><option value="">All areas</option>${areaIds().map(a => `<option value="${esc(a)}"${st.area === a ? " selected" : ""}>${esc(S.areas[a].label)}</option>`).join("")}</select>
      <select id="rc-type" aria-label="Type"><option value="">All types</option>${CONCEPT_TYPES.map(t => `<option value="${t}"${st.type === t ? " selected" : ""}>${esc(typePlural(t))}</option>`).join("")}</select>
      <label class="check-row"><input type="checkbox" id="rc-missed"${st.missedOnly ? " checked" : ""}>Only notes I missed more often than I knew</label></div>
      <div id="rc-card"></div><div class="progress" id="rc-progress"></div></div>`, "Recall");
    let current = null;
    const pick = () => {
      const list = pool();
      if (!list.length) return null;
      const score = n => { const r = record[n.id]; if (!r) return -1 + Math.random() * 0.1; return (r.k - r.m) + (Date.now() - (r.t || 0)) / -8.64e7 + Math.random() * 0.3; };
      const candidates = list.filter(n => !current || n.id !== current.id);
      return (candidates.length ? candidates : list).sort((a, b) => score(a) - score(b))[0];
    };
    const excerpt = n => {
      const tmp = document.createElement("div");
      tmp.innerHTML = n.html || "";
      const names = ["definition", "statement", "procedure", "setup", "problem"];
      const h = $$("h2", tmp).find(x => names.includes(x.textContent.trim().toLowerCase()));
      if (!h) return "";
      const parts = [];
      let el = h.nextElementSibling;
      while (el && el.tagName !== "H2" && parts.length < 4) { parts.push(el.outerHTML); el = el.nextElementSibling; }
      return parts.join("");
    };
    const progress = () => {
      const list = pool();
      const seen = list.filter(n => record[n.id]).length;
      const knew = list.reduce((s, n) => s + ((record[n.id] || {}).k || 0), 0);
      const missed = list.reduce((s, n) => s + ((record[n.id] || {}).m || 0), 0);
      $("#rc-progress").innerHTML = `<span>${seen} of ${list.length} notes seen</span><span>${knew} answers known</span><span>${missed} missed</span>${seen ? '<button type="button" class="linkish" id="rc-reset">Forget my answers</button>' : ""}`;
      const reset = $("#rc-reset");
      if (reset) reset.addEventListener("click", () => { Object.keys(record).forEach(k => delete record[k]); store.set("recall", record); next(); });
    };
    const show = (revealed) => {
      const box = $("#rc-card");
      if (!current) { box.innerHTML = '<div class="empty"><p>No notes match these settings.</p></div>'; progress(); return; }
      const n = current;
      const symbols = (n.notation || []).map(x => `<span class="math math-inline">${esc(x.symbol)}</span>`).join(", ");
      box.innerHTML = `<div class="flash"><div class="chips">${typeChip(n)}${n.area ? `<span class="chip"><span class="swatch" style="background:${areaColor(n.area)}"></span>${esc(S.areas[n.area].label)}</span>` : ""}</div>
        <h2>${n.title_html}</h2>${symbols ? `<p class="prompt">Notation: ${symbols}</p>` : ""}
        ${revealed ? `<div class="answer"><p class="lead">${n.summary_html || ""}</p><div class="prose">${excerpt(n)}</div><p><a href="#/n/${enc(n.id)}">Open the note</a></p></div>` : '<p class="prompt">State the definition or the result in your own words, then reveal the note.</p>'}</div>
        <div class="recall-actions">${revealed ? '<button type="button" class="btn btn-primary" id="rc-knew">I knew it <kbd>1</kbd></button><button type="button" class="btn" id="rc-missed-btn">Not yet <kbd>2</kbd></button>' : '<button type="button" class="btn btn-primary" id="rc-reveal">Reveal <kbd>Space</kbd></button><button type="button" class="btn" id="rc-skip">Skip</button>'}</div>`;
      typeset(box);
      if (revealed) {
        $("#rc-knew").addEventListener("click", () => answer(true));
        $("#rc-missed-btn").addEventListener("click", () => answer(false));
      } else {
        $("#rc-reveal").addEventListener("click", () => show(true));
        $("#rc-skip").addEventListener("click", next);
      }
      box.dataset.revealed = revealed ? "1" : "";
      progress();
    };
    const answer = knew => {
      const r = record[current.id] || { k: 0, m: 0, t: 0 };
      if (knew) r.k += 1; else r.m += 1;
      r.t = Date.now();
      record[current.id] = r;
      store.set("recall", record);
      next();
    };
    const next = () => { current = pick(); show(false); };
    const onKey = ev => {
      if (isTyping(ev.target) || !current) return;
      const revealed = $("#rc-card").dataset.revealed === "1";
      if (ev.key === " " && !revealed) { ev.preventDefault(); show(true); }
      else if (ev.key === "1" && revealed) answer(true);
      else if (ev.key === "2" && revealed) answer(false);
    };
    document.addEventListener("keydown", onKey);
    teardown.push(() => document.removeEventListener("keydown", onKey));
    const refilter = () => { st.area = $("#rc-area").value; st.type = $("#rc-type").value; st.missedOnly = $("#rc-missed").checked; store.set("recall-filter", st); next(); };
    ["#rc-area", "#rc-type", "#rc-missed"].forEach(sel => $(sel).addEventListener("change", refilter));
    next();
  }

  // ============================ Guide
  function showGuide() {
    const types = Object.entries(S.types || {}).map(([k, t]) => `<tr><td>${shapeIcon(k)}</td><td><strong>${esc(t.label)}</strong></td><td>${esc(t.description || "")}</td></tr>`).join("");
    const rels = Object.entries(REL).map(([k, r]) => `<tr><td><strong>${esc(r.label)}</strong><div class="small muted">${esc(r.symmetric ? "symmetric" : "inverse: " + r.inverse)}</div></td><td>${esc(r.definition)}</td><td class="small"><code>${esc(r.example || "")}</code></td></tr>`).join("");
    const levelTable = (kind, chip) => `<div class="table-wrap"><table class="data-table"><tbody>${Object.entries(LEVELS[kind] || {}).map(([k, lv]) => `<tr><td style="white-space:nowrap">${chip(k, lv)}</td><td>${esc(lv.definition)}</td></tr>`).join("")}</tbody></table></div>`;
    const areas = areaIds().map(a => `<tr><td><span class="swatch" style="background:${areaColor(a)}"></span></td><td><strong>${esc(S.areas[a].label)}</strong><div class="small muted">${esc(S.areas[a].description || "")}</div></td><td class="small">${Object.entries(S.domains || {}).filter(([, d]) => d.area === a).map(([, d]) => esc(d.label)).join(", ")}</td></tr>`).join("");
    const flow = `<svg viewBox="0 0 900 150" role="img" aria-label="From a source to this site in five steps">
      <defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 10 5 0 10Z" fill="currentColor"/></marker></defs>
      ${[["Source", "PDF, DOI, arXiv or ISBN"], ["Ingest", "identify the work"], ["Write", "concepts and relations"], ["Check", "rules and numbers"], ["Build", "this site"]].map(([t, s], i) => `
        <g transform="translate(${10 + i * 180},30)"><rect width="150" height="80" rx="12" fill="none" stroke="currentColor" stroke-width="1.5"/>
        <text x="75" y="36" text-anchor="middle" font-family="Jost, sans-serif" font-size="18" font-weight="600" fill="currentColor">${t}</text>
        <text x="75" y="58" text-anchor="middle" font-family="Libre Franklin, sans-serif" font-size="12.5" fill="currentColor" opacity=".8">${s}</text></g>
        ${i < 4 ? `<line x1="${162 + i * 180}" y1="70" x2="${188 + i * 180}" y2="70" stroke="currentColor" stroke-width="1.5" marker-end="url(#arr)"/>` : ""}`).join("")}
    </svg>`;
    setView("guide", `<div class="page-head"><div><h1>Guide</h1><p>How to read the notes, the graph and the levels that qualify each note.</p></div></div>
      <div class="guide">
        <section><h2>What this site is</h2>
          <p>${esc(META.title || "This knowledge base")} is a personal body of scientific knowledge kept as plain text files. Each note holds one concept, theorem, method, model or worked example. Notes are joined by typed relations, cite the sources that support their claims, and state their maturity and their confidence, so that a reader can judge how far to rely on each one.</p>
          <p>Numerical claims are reproduced by scripts whose last run is shown on the note, and every work cited is listed on the References page with the sentences that cite it.</p>
          <div class="flow">${flow}</div></section>
        <section><h2>Note types</h2><div class="table-wrap"><table class="data-table"><tbody>${types}</tbody></table></div></section>
        <section><h2>Relations</h2><p>A relation is read from the note in which it is written: <em>pushforward measure requires measurable map</em>. The other note shows the inverse reading. In the graph, structural relations are solid grey, analytic relations solid violet, applications dashed, contrasts dotted and historical links long-dashed. Mentions are thin lines for links inside the text.</p>
          <div class="table-wrap"><table class="data-table"><thead><tr><th>Relation</th><th>Meaning</th><th>Example</th></tr></thead><tbody>${rels}</tbody></table></div></section>
        <section><h2>Maturity</h2>${levelTable("maturity", (k, lv) => `<span class="chip">${pips(lv.rank, 3)}${esc(lv.label)}</span>`)}</section>
        <section><h2>Confidence</h2>${levelTable("confidence", k => confidenceChip({ confidence: k }))}</section>
        <section><h2>Questions</h2>${levelTable("status", (k, lv) => esc(lv.label))}</section>
        <section><h2>Ontology and epistemology</h2><p>Each developing note has a section on ontology and epistemology. The section asks what the concept takes to exist and how claims about it can be known, and the note lists the themes that it raises: the mark ${DIM_MARK.ontology} labels an ontological theme and the mark ${DIM_MARK.epistemology} an epistemological one. The Philosophy page gathers the themes with the notes that raise them.</p>
          <div class="table-wrap"><table class="data-table"><tbody>${Object.entries(PHILO.themes || {}).map(([id, t]) => `<tr><td style="white-space:nowrap">${themeChip(id)}</td><td>${esc(t.definition)}</td></tr>`).join("")}</tbody></table></div></section>
        <section><h2>Areas and domains</h2><div class="table-wrap"><table class="data-table"><tbody>${areas}</tbody></table></div></section>
        <section><h2>Using the site</h2><div class="table-wrap"><table class="data-table"><tbody>
          <tr><td><kbd>/</kbd> or <kbd>Ctrl</kbd> <kbd>K</kbd></td><td>Search titles, aliases, summaries, notation and text.</td></tr>
          <tr><td><kbd>?</kbd></td><td>Open this guide.</td></tr>
          <tr><td><kbd>Esc</kbd></td><td>Close the search, the preview or the note card.</td></tr>
          <tr><td>Hover a link</td><td>Preview the note or the reference without leaving the page.</td></tr>
          <tr><td>Explore</td><td>Click a node for its card, double-click to open it, and use the filters to show one area, one type or a neighbourhood.</td></tr>
          <tr><td>Paths</td><td>List the prerequisites of a note, connect two notes, or show the neighbourhood of a note.</td></tr>
          <tr><td>Symbols</td><td>Browse every symbol as a tile, read its meaning, its convention, its defining note and the formulas that reuse it, copy its TeX, or switch to the table.</td></tr>
          <tr><td>Philosophy</td><td>Browse the ontological and epistemological themes, and read the analyses of the notes that raise a theme.</td></tr>
          <tr><td>References</td><td>Follow the works on a timeline, filter and group them, read the sentences that cite each work, copy its BibTeX entry, or draw the citation map.</td></tr>
          <tr><td>Recall</td><td>Test yourself on definitions and results; the record stays in this browser.</td></tr>
        </tbody></table></div></section>
      </div>`, "Guide");
  }

  // ============================ Router
  function route() {
    clearView();
    const hash = location.hash || "#/explore";
    let m;
    if ((m = hash.match(/^#\/n\/([^/]+)(?:\/(.+))?$/))) return showNote(decodeURIComponent(m[1]), m[2] ? decodeURIComponent(m[2]) : null);
    if ((m = hash.match(/^#\/wanted\/(.+)$/))) return showWanted(decodeURIComponent(m[1]));
    // The glossary and the symbols were pages inside the notes index and are now tabs of their
    // own; the old addresses still work and go to the new pages.
    if ((m = hash.match(/^#\/notes(?:\/(cards|glossary|notation))?$/))) {
      if (m[1] === "glossary") return showGlossary();
      if (m[1] === "notation") return showSymbols(null);
      return showNotes();
    }
    if ((m = hash.match(/^#\/symbols(?:\/(.+))?$/))) return showSymbols(m[1] ? decodeURIComponent(m[1]) : null);
    if ((m = hash.match(/^#\/paths(?:\/(prereq|connect|near))?(?:\/([^/]*))?(?:\/([^/]*))?$/))) return showPaths(m[1] || "prereq", m[2] ? decodeURIComponent(m[2]) : "", m[3] ? decodeURIComponent(m[3]) : "");
    if ((m = hash.match(/^#\/explore\/(.+)$/))) return showExplore(decodeURIComponent(m[1]));
    if ((m = hash.match(/^#\/references(?:\/(.+))?$/))) return showReferences(m[1] ? decodeURIComponent(m[1]) : null);
    if ((m = hash.match(/^#\/philosophy(?:\/(.+))?$/))) return showPhilosophy(m[1] ? decodeURIComponent(m[1]) : null);
    const views = { "#/explore": showExplore, "#/glossary": showGlossary, "#/notes": showNotes, "#/timeline": showTimeline, "#/maps": showMaps, "#/sources": () => showReferences(null), "#/questions": showQuestions, "#/state": showState, "#/recall": showRecall, "#/guide": showGuide };
    return (views[hash] || showExplore)();
  }

  setupChrome();
  window.addEventListener("hashchange", route);
  route();
})();
