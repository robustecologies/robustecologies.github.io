// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
// Loads the engine sources into a fresh context and reports checks.
import { readFileSync } from "fs";
import { dirname, join } from "path";
import { fileURLToPath } from "url";
import vm from "vm";

export const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");
export const CORE = ["src/core/rng.js", "src/core/expr.js", "src/core/sim.js"];

export function load(files) {
  const ctx = { console, Math, Float64Array, Uint8Array, Uint32Array, Array, Object, JSON };
  ctx.globalThis = ctx;
  vm.createContext(ctx);
  for (const f of files) vm.runInContext(readFileSync(join(ROOT, f), "utf8"), ctx, { filename: f });
  return ctx.RElabFlow;
}

export function reporter(title) {
  const rows = [];
  let fails = 0;
  return {
    check(name, pass, detail) {
      if (!pass) fails++;
      const row = `${pass ? "PASS" : "FAIL"}  ${name}: ${detail}`;
      rows.push(row);
      console.log(row);
    },
    done() {
      const s = `${title}: ${rows.length - fails} of ${rows.length} checks pass`;
      console.log(s);
      return { rows: rows.concat([s]), fails };
    }
  };
}
