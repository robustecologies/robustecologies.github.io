// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Seeded random numbers. Every stochastic scene is reproducible from its
   seed: the same seed gives the same noise, the same jumps and the same
   initial ensemble in the studio, in an exported page and in the tests. */
(function (DF) {
  "use strict";

  // splitmix32 expands one 32-bit seed into the four words of xoshiro128**.
  function splitmix32(a) {
    return function () {
      a |= 0; a = (a + 0x9e3779b9) | 0;
      let z = a;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b);
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35);
      return (z ^ (z >>> 16)) >>> 0;
    };
  }

  // xoshiro128** (Blackman and Vigna), period 2^128 - 1.
  function RNG(seed) {
    const sm = splitmix32(seed === undefined ? 1 : seed >>> 0);
    this.s = new Uint32Array([sm(), sm(), sm(), sm()]);
    this._spare = null;
  }
  RNG.prototype.nextU32 = function () {
    const s = this.s;
    const r = Math.imul(rotl(Math.imul(s[1], 5), 7), 9) >>> 0;
    const t = s[1] << 9;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 11);
    return r;
  };
  function rotl(x, k) { return (x << k) | (x >>> (32 - k)); }

  // Uniform on (0, 1), never exactly 0 or 1, with 53 random bits.
  RNG.prototype.uniform = function () {
    const hi = this.nextU32() >>> 5, lo = this.nextU32() >>> 6;
    return (hi * 67108864 + lo + 0.5) / 9007199254740992;
  };
  RNG.prototype.range = function (a, b) { return a + (b - a) * this.uniform(); };

  // Standard normal by the Box-Muller transform, one value cached.
  RNG.prototype.normal = function () {
    if (this._spare !== null) { const v = this._spare; this._spare = null; return v; }
    const u = this.uniform(), v = this.uniform();
    const r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * v;
    this._spare = r * Math.sin(th);
    return r * Math.cos(th);
  };
  RNG.prototype.exponential = function (rate) { return -Math.log(this.uniform()) / rate; };

  // Poisson by inversion for small means and a normal approximation above 60.
  RNG.prototype.poisson = function (mu) {
    if (mu <= 0) return 0;
    if (mu > 60) return Math.max(0, Math.round(mu + Math.sqrt(mu) * this.normal()));
    const L = Math.exp(-mu);
    let k = 0, p = 1;
    do { k++; p *= this.uniform(); } while (p > L);
    return k - 1;
  };

  // Symmetric alpha-stable variate, scale 1, by Chambers, Mallows and Stuck
  // (1976); alpha = 2 gives a normal with variance 2, alpha = 1 a Cauchy.
  RNG.prototype.stable = function (alpha) {
    const V = Math.PI * (this.uniform() - 0.5), W = this.exponential(1);
    if (Math.abs(alpha - 1) < 1e-12) return Math.tan(V);
    return Math.sin(alpha * V) / Math.pow(Math.cos(V), 1 / alpha) *
      Math.pow(Math.cos(V - alpha * V) / W, (1 - alpha) / alpha);
  };

  DF.RNG = RNG;
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
