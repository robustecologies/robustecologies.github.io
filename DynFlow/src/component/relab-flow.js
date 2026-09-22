// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* <relab-flow>: a scene as an HTML element. One script tag and one element
   place a live figure in a web page, a Quarto or R Markdown document, a
   reveal.js slide or a pkgdown article.

     <script src="relabflow.js"></script>
     <relab-flow scene='{"system": "...", "view": {"type": "trajectory"}}'></relab-flow>
     <relab-flow model="lorenz" theme="blackboard" controls></relab-flow>
     <relab-flow src="figures/fold.json" paused></relab-flow>

   Attributes: scene (JSON), src (URL of a scene file), model (catalogue
   identifier), theme, view, title, controls (show play, restart and
   fullscreen buttons), paused (do not start automatically), static (draw a
   still after a number of frames given by the attribute, default 240).
   Changing scene, src, model, theme, view or title loads the scene again. The
   element exposes .player, .play(), .pause(), .restart(), .setParam(name, v)
   and .scene. A figure starts only when it is visible and shows a still frame
   when the reader asks for reduced motion. */
(function (DF) {
  "use strict";
  if (typeof HTMLElement === "undefined" || typeof customElements === "undefined") return;

  const CSS = ":host{display:block;position:relative;min-height:200px;height:360px;contain:content}" +
    ".stage{position:absolute;inset:0}" +
    ".ctl{position:absolute;right:10px;bottom:10px;display:flex;gap:6px;opacity:0;transition:opacity .2s;z-index:5}" +
    ":host(:hover) .ctl,.ctl:focus-within{opacity:1}" +
    "button{all:unset;cursor:pointer;width:30px;height:30px;border-radius:50%;display:grid;place-items:center;background:rgba(23,12,58,.72);color:#fff;font:14px system-ui}" +
    "button:hover{background:#EE6A24}button:focus-visible{outline:2px solid #FB9E07}" +
    ".err{position:absolute;inset:0;display:grid;place-items:center;padding:16px;font:13px system-ui;color:#CF4446;background:#fff8f5;text-align:center}";

  class RElabFlowElement extends HTMLElement {
    static get observedAttributes() { return ["scene", "src", "model", "theme", "view", "title"]; }
    /* The shadow tree and the player outlive a move of the element: moving it
       to another parent (a slide framework does this) pauses the figure and
       resumes it where it was, and does not build a second shadow root. */
    connectedCallback() {
      if (this._root) { if (this.player && this._wasRunning) this.player.play(); return; }
      this._root = this.shadowRoot || this.attachShadow({ mode: "open" });
      this._root.innerHTML = "";
      const style = document.createElement("style"); style.textContent = CSS;
      this._stage = document.createElement("div"); this._stage.className = "stage";
      this._root.append(style, this._stage);
      if (this.hasAttribute("controls")) this.buildControls();
      this.loadScene();
    }
    disconnectedCallback() { if (this.player) { this._wasRunning = this.player.running; this.player.pause(); } }
    attributeChangedCallback(name, oldV, newV) { if (this._root && oldV !== newV && (this.player || this._failed)) this.loadScene(); }

    buildControls() {
      const c = document.createElement("div"); c.className = "ctl";
      const mk = (label, title, fn) => { const b = document.createElement("button"); b.textContent = label; b.title = title; b.setAttribute("aria-label", title); b.addEventListener("click", fn); c.appendChild(b); return b; };
      this._playBtn = mk("\u275a\u275a", "Pause", () => this.player && this.player.toggle());
      mk("\u21bb", "Restart", () => this.restart());
      mk("\u26f6", "Full screen", () => { if (document.fullscreenElement) document.exitFullscreen(); else this.requestFullscreen && this.requestFullscreen(); });
      this._root.appendChild(c);
    }

    async loadScene() {
      let scene;
      try {
        if (this.hasAttribute("scene")) scene = JSON.parse(this.getAttribute("scene"));
        else if (this.hasAttribute("src")) {
          const src = this.getAttribute("src");
          let res;
          try { res = await fetch(src); } catch (e) { throw new Error("Could not load " + src + ". A page opened from disk cannot fetch files: serve the folder, or give the scene attribute instead"); }
          if (!res.ok) throw new Error("Could not load " + src + " (HTTP " + res.status + ")");
          try { scene = await res.json(); } catch (e) { throw new Error(src + " is not a scene file (JSON)"); }
        }
        else if (this.hasAttribute("model")) {
          if (!DF.sceneFor) throw new Error("The model catalogue is not loaded");
          scene = DF.sceneFor(this.getAttribute("model"));
        } else throw new Error("Give a scene, src or model attribute");
      } catch (e) { this._failed = true; this.showError(e); return; }
      if (this.hasAttribute("theme")) scene.style = Object.assign({}, scene.style, { theme: this.getAttribute("theme") });
      if (this.hasAttribute("view")) scene.view = Object.assign({}, scene.view, { type: this.getAttribute("view") });
      if (this.hasAttribute("title")) scene.overlay = Object.assign({}, scene.overlay, { title: this.getAttribute("title") });
      if (this.player) { this.player.dispose(); this.player = null; }
      this._stage.innerHTML = "";
      try {
        this.player = new DF.Player(this._stage, scene);
        if (this.player.error) throw this.player.error;
      } catch (e) { if (this.player) { this.player.dispose(); this.player = null; } this._failed = true; this.showError(e); return; }
      this._failed = false;
      this.player.on("state", (s) => { if (this._playBtn) { this._playBtn.textContent = s === "play" ? "\u275a\u275a" : "\u25b6"; this._playBtn.title = s === "play" ? "Pause" : "Play"; } });
      this.player.on("error", (e) => this.showError(e));
      const reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      if (this.hasAttribute("static") || reduced) this.player.advance(parseInt(this.getAttribute("static"), 10) || 240);
      else if (!this.hasAttribute("paused")) this.player.play();
      this.dispatchEvent(new CustomEvent("ready", { detail: this.player }));
    }
    showError(e) {
      const d = document.createElement("div"); d.className = "err";
      d.textContent = "RElabFlow: " + (e && e.message ? e.message : String(e));
      this._stage.innerHTML = ""; this._stage.appendChild(d);
    }
    get scene() { return this.player ? this.player.getScene() : null; }
    play() { if (this.player) this.player.play(); }
    pause() { if (this.player) this.player.pause(); }
    restart(seed) { if (this.player) { this.player.restart(seed); if (!this.hasAttribute("paused")) this.player.play(); } }
    setParam(name, value) { if (this.player) this.player.setParam(name, value); }
  }

  if (!customElements.get("relab-flow")) customElements.define("relab-flow", RElabFlowElement);
  DF.RElabFlowElement = RElabFlowElement;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
