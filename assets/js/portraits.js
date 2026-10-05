/* Phase portraits for the home hero, played by RElabFlow
   (shared/vendor/relabflow/relabflow.js, a copy of the RElabFlow bundle
   without its embedded source text). One model is drawn at random on each
   visit, never the same one twice in a row, from every chaotic flow of the
   catalogue and its chaotic ecological models; ?portrait=<id> forces a
   choice. The caption is the catalogue's own description, whose numbers
   RElabFlow's model tests check. */
(function (root) {
  "use strict";
  const DF = root.RElabFlow;
  const ECOLOGY = ["may-leonard", "hastings-powell", "vano-lv4", "act-lv", "gilpin", "mccann-yodzis", "blasius-upca", "hindmarsh-rose"];
  // Flows whose default scene is not an attractor seen in 3D.
  const SKIP = ["charney-devore", "nose-hoover"];

  function list() {
    if (!DF || !DF.CATALOGUE) return [];
    // A pool chosen in RElabSite (window.HERO_MODELS) replaces the default one.
    const chosen = root.HERO_MODELS;
    if (chosen && chosen.length) return DF.CATALOGUE.filter(function (m) { return chosen.indexOf(m.id) >= 0; });
    return DF.CATALOGUE.filter(function (m) {
      return (m.group === "Chaotic flows" || ECOLOGY.indexOf(m.id) >= 0) && SKIP.indexOf(m.id) < 0;
    });
  }

  // Short source for a caption: the authors and year of the original work.
  function credit(m) {
    const hit = m.source.match(/\(([^()]*\d{4})\)\s*$/) || m.source.match(/^(.*?\(\d{4}\))/);
    if (!hit) return "";
    return hit[1].indexOf("(") >= 0 ? hit[1].replace(/\s*\((\d{4})\)/, " $1") : hit[1];
  }

  function choose() {
    const L = list();
    if (!L.length) return null;
    let forced = null, last = null;
    try { forced = new URLSearchParams(root.location.search).get("portrait"); } catch (e) { forced = null; }
    const hit = L.filter(function (m) { return m.id === forced; })[0];
    if (hit) return hit;
    try { last = root.localStorage.getItem("relab:portrait"); } catch (e) { last = null; }
    const pool = L.filter(function (m) { return m.id !== last; });
    const pick = pool[Math.floor(Math.random() * pool.length)];
    try { root.localStorage.setItem("relab:portrait", pick.id); } catch (e) { /* storage unavailable */ }
    return pick;
  }

  function caption(m) {
    const c = credit(m);
    return "Live: " + m.name + (c ? " (" + c + ")" : "") + ". " + m.about + (m.scene.view && m.scene.view.projection === "simplex" ? " Click to seed new initial conditions." : " Drag to turn the attractor.");
  }

  // Play model m in host, with the lab's night theme and no overlay text.
  function render(host, m) {
    const scene = DF.sceneFor(m.id);
    scene.overlay = { title: "", legend: false, equations: false, readout: false };
    scene.style = Object.assign({}, scene.style, { theme: "relab-night", alpha: 0.5 });
    scene.n = Math.min(scene.n || 1400, 1400);
    const player = new DF.Player(host, scene);
    player.play();
    // The hero pauses once less than a quarter of it shows, so that it does
    // not compete with the figures further down the page.
    if ("IntersectionObserver" in root) {
      new IntersectionObserver(function (e) {
        if (e[0].intersectionRatio >= 0.25) player.play(); else player.pause();
      }, { threshold: [0, 0.25, 0.5] }).observe(host);
    }
    return player;
  }

  root.Portraits = { list: list, choose: choose, caption: caption, render: render };
})(this);
