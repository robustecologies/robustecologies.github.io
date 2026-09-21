// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Concatenates the sources into dist/dynflow.js, one classic script that
// works from file:// and from any server. The bundle carries its own source
// text as DynFlow.SOURCE, which the standalone-page export embeds.
// Run: node tools/build.mjs
import { readFileSync, writeFileSync, mkdirSync } from "fs";
import { dirname, join } from "path";
import { fileURLToPath } from "url";

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");
const files = JSON.parse(readFileSync(join(ROOT, "tools", "files.json"), "utf8"))
  .concat(["src/render/export.js", "src/models/catalogue.js", "src/component/dyn-flow.js"]);
const version = JSON.parse(readFileSync(join(ROOT, "package.json"), "utf8")).version;
const head = `// DynFlow ${version}. Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab.\n` +
  "// SPDX-License-Identifier: GPL-3.0-or-later. https://www.gnu.org/licenses/gpl-3.0.html\n" +
  "// Built from " + files.length + " source files by tools/build.mjs; edit the sources, not this file.\n";
const body = files.map((f) => "// ---- " + f + "\n" + readFileSync(join(ROOT, f), "utf8")).join("\n");
const bundle = head + body + "\nglobalThis.DynFlow.VERSION = " + JSON.stringify(version) + ";\n";
const out = bundle + "globalThis.DynFlow.SOURCE = " + JSON.stringify(bundle) + ";\n";
mkdirSync(join(ROOT, "dist"), { recursive: true });
writeFileSync(join(ROOT, "dist", "dynflow.js"), out);
console.log(`dist/dynflow.js: ${files.length} files, ${(out.length / 1024).toFixed(0)} KiB (engine ${(bundle.length / 1024).toFixed(0)} KiB)`);
