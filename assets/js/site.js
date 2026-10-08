// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Behaviour of the lab site: phone menu, live title bands, the RElabFlow
   hero and goal figures, category filters, copy buttons for code and the
   search of the site. Every page reads without it. */
(function () {
  "use strict";
  const V = window.Viz, DF = window.RElabFlow;

  // ------------------------------------------------------------ phone menu
  const top = document.querySelector(".top"), toggle = document.querySelector(".nav-toggle");
  if (top && toggle) {
    toggle.addEventListener("click", function () {
      const open = top.classList.toggle("open");
      toggle.setAttribute("aria-expanded", String(open));
    });
  }

  // ------------------------------------------------------------ menus
  // People and Projects open their list; they do not lead to a page.
  const menus = Array.from(document.querySelectorAll(".has-menu"));
  menus.forEach(function (li) {
    const b = li.querySelector(".menu-btn");
    b.addEventListener("click", function (e) {
      e.stopPropagation();
      const open = !li.classList.contains("open");
      menus.forEach(function (o) { o.classList.remove("open"); o.querySelector(".menu-btn").setAttribute("aria-expanded", "false"); });
      li.classList.toggle("open", open); b.setAttribute("aria-expanded", String(open));
    });
  });
  document.addEventListener("click", function () { menus.forEach(function (o) { o.classList.remove("open"); o.querySelector(".menu-btn").setAttribute("aria-expanded", "false"); }); });

  // ------------------------------------------------------------ title bands
  document.querySelectorAll(".band-canvas").forEach(function (c) {
    if (V) V.heteroclinicField(c, { n: 700, scale: 1.3, dx: 0.28, dy: 0.1, fade: 0.05, interactive: false, lineWidth: 1, nearEq: 0 });
  });

  // ------------------------------------------------------------ home hero
  const stage = document.querySelector(".hero-canvas");
  if (stage && window.Portraits) {
    const m = window.Portraits.choose();
    if (m) {
      stage.setAttribute("aria-label", "Live phase portrait: " + m.name);
      const cap = document.querySelector(".hero-caption");
      if (cap) cap.textContent = window.Portraits.caption(m);
      window.Portraits.render(stage, m);
    }
  }
  const over = document.querySelector(".top.over-hero");
  if (over) {
    const onScroll = function () { over.classList.toggle("solid", window.scrollY > window.innerHeight * 0.6); };
    window.addEventListener("scroll", onScroll, { passive: true }); onScroll();
  }

  // ------------------------------------------------------------ live figures
  // Exact scenes; every number in the legends beside them was checked
  // numerically (the torus by its stroboscopic section, the focus by its
  // eigenvalues, the snapshot attractor by its exponent).
  // Every live figure carries its RElabFlow scene as JSON inside the element.
  function play(host, scene) {
    if (!DF || !DF.Player || !scene) return;
    const sc = JSON.parse(JSON.stringify(scene));
    sc.overlay = { title: "", legend: false, equations: false, readout: false };
    sc.style = Object.assign({ theme: "relab-night" }, sc.style);
    const el = document.createElement("relab-flow");
    el.setAttribute("scene", JSON.stringify(sc));
    el.style.height = "100%";
    host.appendChild(el);
  }
  document.querySelectorAll(".live-fig").forEach(function (host) {
    const s = host.querySelector('script[type="application/json"]');
    if (s) play(host, JSON.parse(s.textContent));
  });
  // The cusp instrument: bifurcation diagram, slider for r and a slow sweep.
  document.querySelectorAll('[data-instrument="cusp"]').forEach(function (host) {
    if (V && V.cuspInstrument) V.cuspInstrument(host, { ink: "#eceaf4" });
  });

  // ------------------------------------------------------------ licensing
  // The panel opens on hover (CSS) and on click or Enter; Escape or a click
  // elsewhere closes it.
  const lic = document.querySelector(".lic"), licBtn = lic && lic.querySelector(".lic-btn");
  if (licBtn) {
    const setLic = function (open) { lic.classList.toggle("open", open); licBtn.setAttribute("aria-expanded", String(open)); };
    licBtn.addEventListener("click", function (e) { e.stopPropagation(); setLic(!lic.classList.contains("open")); });
    lic.querySelector(".lic-panel").addEventListener("click", function (e) { e.stopPropagation(); });
    document.addEventListener("click", function () { setLic(false); });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape" && lic.classList.contains("open")) { setLic(false); licBtn.focus(); } });
  }

  // ------------------------------------------------------------ theme
  // The switch moves between the light and dark themes and remembers the choice.
  const themeBtn = document.querySelector(".theme-toggle");
  if (themeBtn) {
    themeBtn.addEventListener("click", function () {
      const root = document.documentElement;
      const dark = root.dataset.theme ? root.dataset.theme === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
      root.dataset.theme = dark ? "light" : "dark";
      try { localStorage.setItem("relab:theme", root.dataset.theme); } catch (e) { /* storage unavailable */ }
    });
  }

  // ------------------------------------------------------------ cite
  document.querySelectorAll(".copy-cite").forEach(function (b) {
    b.addEventListener("click", function () {
      const block = b.closest(".cite-block"), t = block.querySelector("code, .cite-text");
      if (!navigator.clipboard || !t) return;
      navigator.clipboard.writeText(t.innerText).then(function () { b.textContent = "Copied"; setTimeout(function () { b.textContent = "Copy"; }, 1600); });
    });
  });

  // ------------------------------------------------------------ categories
  // Buttons of the category column and the chips of each item filter the
  // list; #category=Name in the address selects one on arrival.
  const buttons = Array.from(document.querySelectorAll(".cats [data-filter]"));
  function filter(name) {
    buttons.forEach(function (b) { b.setAttribute("aria-pressed", String(b.dataset.filter === name)); });
    document.querySelectorAll("[data-cats]").forEach(function (it) {
      it.hidden = !!name && it.dataset.cats.split("|").indexOf(name) < 0;
    });
    document.querySelectorAll(".listing section.level2").forEach(function (sec) {
      const items = sec.querySelectorAll("[data-cats]");
      if (items.length) sec.hidden = Array.from(items).every(function (i) { return i.hidden; });
    });
  }
  // On a phone the category column starts folded, above the list.
  const catBox = document.querySelector(".cats details");
  if (catBox && window.matchMedia("(max-width: 860px)").matches) catBox.open = false;
  if (catBox) catBox.querySelector("summary").addEventListener("click", function (e) { if (!window.matchMedia("(max-width: 860px)").matches) e.preventDefault(); });
  if (buttons.length) {
    buttons.forEach(function (b) { b.addEventListener("click", function () { filter(b.dataset.filter); }); });
    const fromHash = function () {
      const m = decodeURIComponent(location.hash).match(/^#category=(.+)$/);
      if (m) filter(m[1].trim());
    };
    window.addEventListener("hashchange", fromHash); fromHash();
  }

  // ------------------------------------------------------------ code
  document.querySelectorAll(".doc-body pre, main > pre").forEach(function (pre) {
    if (!navigator.clipboard) return;
    const b = document.createElement("button");
    b.type = "button"; b.className = "copy-btn"; b.textContent = "Copy";
    b.addEventListener("click", function () {
      const code = pre.querySelector("code") || pre;
      navigator.clipboard.writeText(code.innerText).then(function () {
        b.textContent = "Copied"; setTimeout(function () { b.textContent = "Copy"; }, 1600);
      });
    });
    pre.appendChild(b);
  });

  // ------------------------------------------------------------ search
  // One index (assets/js/search-index.js) serves the search panel of the
  // header and the search of the 404 page. Addresses in the index are
  // relative to the site root; data-root on <html> leads from this page to it.
  const pages = window.SEARCH_INDEX || [];
  const root = document.documentElement.dataset.root || "";
  const norm = function (s) { return s.normalize("NFD").replace(/[\u0300-\u036f]/g, "").toLowerCase(); };
  const escH = function (s) { return s.replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); };
  pages.forEach(function (p) { p.nt = norm(p.t); p.nc = norm(p.c.join(" ")); p.nx = norm(p.d + " " + p.x); });
  function snippet(p, terms) {
    const text = p.d + " " + p.x, low = norm(text);
    let at = -1;
    terms.forEach(function (t) { const i = low.indexOf(t); if (i >= 0 && (at < 0 || i < at)) at = i; });
    if (at < 0) return escH(p.d || text.slice(0, 160));
    const a = Math.max(0, at - 70), b = Math.min(text.length, at + 110);
    let s = escH((a ? "\u2026" : "") + text.slice(a, b) + (b < text.length ? "\u2026" : ""));
    terms.forEach(function (t) {
      const re = new RegExp("(" + t.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + ")", "gi");
      s = s.replace(re, "<mark>$1</mark>");
    });
    return s;
  }
  function search(input, out) {
    const terms = norm(input.value).split(/\s+/).filter(function (t) { return t.length > 1; });
    out.innerHTML = "";
    if (!terms.length) return;
    const hits = pages.map(function (p) {
      let score = 0;
      const all = terms.every(function (t) {
        const inT = p.nt.indexOf(t) >= 0, inC = p.nc.indexOf(t) >= 0, n = p.nx.split(t).length - 1;
        score += (inT ? 10 : 0) + (inC ? 4 : 0) + Math.min(n, 5);
        return inT || inC || n > 0;
      });
      return { p: p, s: all ? score : 0 };
    }).filter(function (h) { return h.s > 0; }).sort(function (a, b) { return b.s - a.s; }).slice(0, 8);
    if (!hits.length) { out.innerHTML = '<li class="none">No page lies in that basin; try another word.</li>'; return; }
    out.innerHTML = hits.map(function (h) {
      const where = "/" + h.p.u.replace(/index\.html$/, "");
      return '<li><a href="' + root + h.p.u + '">' + escH(h.p.t) + '</a><span class="where">' + escH(where) + '</span><span class="snip">' + snippet(h.p, terms) + "</span></li>";
    }).join("");
  }
  document.querySelectorAll(".site-search-input").forEach(function (input) {
    const out = input.closest("form").parentNode.querySelector(".search-results");
    input.addEventListener("input", function () { search(input, out); });
  });
  // Header panel: the magnifying glass, or the key "/", opens it; Escape closes it.
  const panel = document.getElementById("search-panel"), opener = document.querySelector(".search-toggle");
  if (panel && opener) {
    const field = panel.querySelector("input");
    const show = function (on) {
      panel.hidden = !on; opener.setAttribute("aria-expanded", String(on));
      if (on) field.focus();
    };
    opener.addEventListener("click", function () { show(panel.hidden); });
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && !panel.hidden) { show(false); opener.focus(); }
      if (e.key === "/" && panel.hidden && !/^(INPUT|TEXTAREA)$/.test(document.activeElement.tagName)) { e.preventDefault(); show(true); }
    });
  }
  // 404: the last part of the missing address is the first query.
  const lost = document.getElementById("site-search");
  if (lost) {
    const last = location.pathname.split("/").filter(Boolean).pop();
    if (last && !/^404/.test(last)) {
      lost.value = decodeURIComponent(last).replace(/\.html?$/, "").replace(/[-_]+/g, " ");
      lost.dispatchEvent(new Event("input"));
    }
  }
})();
