/* Dynamical systems behind every animation of the layout study.
   Loaded as a classic script in the browser (window.Flows) and through
   require() in node for checks/flows_check.mjs. Every system is integrated
   with the classical fourth-order Runge-Kutta scheme at a fixed step. */
(function (root) {
  "use strict";

  // One RK4 step of dx/dt = f(x) for a state array x of length n.
  function rk4(f, x, h) {
    const n = x.length;
    const k1 = f(x);
    const y = new Array(n);
    for (let i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k1[i];
    const k2 = f(y);
    for (let i = 0; i < n; i++) y[i] = x[i] + 0.5 * h * k2[i];
    const k3 = f(y);
    for (let i = 0; i < n; i++) y[i] = x[i] + h * k3[i];
    const k4 = f(y);
    const out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = x[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return out;
  }

  // May-Leonard competition among three species (May and Leonard 1975):
  //   dx_i/dt = x_i (1 - x_i - alpha x_{i+1} - beta x_{i+2}),  indices mod 3.
  // For alpha < 1 < beta and alpha + beta > 2 the three single-species
  // equilibria e_i are saddles joined by an attracting heteroclinic cycle, and
  // orbits spend geometrically longer times near each saddle.
  function mayLeonard(alpha, beta) {
    return function (x) {
      return [
        x[0] * (1 - x[0] - alpha * x[1] - beta * x[2]),
        x[1] * (1 - x[1] - alpha * x[2] - beta * x[0]),
        x[2] * (1 - x[2] - alpha * x[0] - beta * x[1])
      ];
    };
  }

  // Lotka-Volterra predator and prey:
  //   dx/dt = a x - b x y,  dy/dt = d x y - c y.
  // H = d x - c ln x + b y - a ln y is constant along orbits.
  function lotkaVolterra(a, b, c, d) {
    return function (z) {
      return [a * z[0] - b * z[0] * z[1], d * z[0] * z[1] - c * z[1]];
    };
  }
  function lvInvariant(a, b, c, d, z) {
    return d * z[0] - c * Math.log(z[0]) + b * z[1] - a * Math.log(z[1]);
  }

  // Cusp normal form dx/dt = r + a x - x^3. For a > 0 the equilibria fold at
  // r = +/- 2 (a/3)^(3/2); between the folds two stable branches coexist.
  function cusp(r, a) {
    return function (x) { return [r + a * x[0] - x[0] * x[0] * x[0]]; };
  }
  function cuspFolds(a) {
    const r = 2 * Math.pow(a / 3, 1.5);
    return [-r, r];
  }
  // Real equilibria of r + a x - x^3 = 0, sorted, found from the cubic.
  function cuspEquilibria(r, a) {
    // x^3 - a x - r = 0, depressed cubic t^3 + p t + q with p = -a, q = -r.
    const p = -a, q = -r;
    const disc = -(4 * p * p * p + 27 * q * q);
    if (disc > 0) {
      const m = 2 * Math.sqrt(-p / 3);
      const th = Math.acos((3 * q) / (p * m)) / 3;
      return [0, 1, 2].map(function (k) { return m * Math.cos(th - (2 * Math.PI * k) / 3); })
        .sort(function (u, v) { return u - v; });
    }
    const s = Math.sqrt((q * q) / 4 + (p * p * p) / 27);
    return [Math.cbrt(-q / 2 + s) + Math.cbrt(-q / 2 - s)];
  }

  // Rosenzweig-MacArthur predator and prey with a Holling type II response,
  // in units where the half-saturation constant and the prey growth rate are 1:
  //   dx/dt = x (1 - x/K) - x y/(1 + x),  dy/dt = y (b x/(1 + x) - d).
  // The equilibrium x* = d/(b - d) loses stability in a Hopf bifurcation at
  // K = 1 + 2 x*; above it every interior orbit but x* tends to a limit cycle.
  function rosenzweigMacArthur(K, b, d) {
    return function (z) {
      const g = z[0] / (1 + z[0]);
      return [z[0] * (1 - z[0] / K) - g * z[1], z[1] * (b * g - d)];
    };
  }

  // Lotka-Volterra competition between two species, symmetric coefficient a:
  //   dx/dt = x (1 - x - a y),  dy/dt = y (1 - y - a x).
  // For a > 1 the interior point (1, 1)/(1 + a) is a saddle whose stable
  // manifold is the diagonal x = y; (1, 0) and (0, 1) are both stable.
  function competition(a) {
    return function (z) {
      return [z[0] * (1 - z[0] - a * z[1]), z[1] * (1 - z[1] - a * z[0])];
    };
  }

  // Hastings-Powell food chain (Hastings and Powell 1991): resource x,
  // consumer y, top predator z, with f_i(u) = a_i u/(1 + b_i u):
  //   dx/dt = x (1 - x) - f1(x) y
  //   dy/dt = f1(x) y - f2(y) z - d1 y
  //   dz/dt = f2(y) z - d2 z.
  // At a1 = 5, b1 = 3, a2 = 0.1, b2 = 2, d1 = 0.4, d2 = 0.01 orbits settle on
  // the chaotic teacup attractor.
  function hastingsPowell(a1, b1, a2, b2, d1, d2) {
    return function (u) {
      const f1 = (a1 * u[0]) / (1 + b1 * u[0]), f2 = (a2 * u[1]) / (1 + b2 * u[1]);
      return [u[0] * (1 - u[0]) - f1 * u[1], f1 * u[1] - f2 * u[2] - d1 * u[1], f2 * u[2] - d2 * u[2]];
    };
  }
  const HP = [5, 3, 0.1, 2, 0.4, 0.01];

  // Lorenz system at the classical parameters sigma = 10, rho = 28, beta = 8/3.
  function lorenz(sigma, rho, beta) {
    return function (x) {
      return [sigma * (x[1] - x[0]), x[0] * (rho - x[2]) - x[1], x[0] * x[1] - beta * x[2]];
    };
  }

  // Barycentric projection of a point of the positive orthant onto an
  // equilateral triangle with vertices e_1 (top), e_2 (lower left), e_3
  // (lower right); returns coordinates in the unit square.
  function simplexToPlane(x) {
    const s = x[0] + x[1] + x[2];
    const w0 = x[0] / s, w1 = x[1] / s, w2 = x[2] / s;
    const V = [[0.5, 0.06], [0.04, 0.86], [0.96, 0.86]];
    return [w0 * V[0][0] + w1 * V[1][0] + w2 * V[2][0], w0 * V[0][1] + w1 * V[1][1] + w2 * V[2][1]];
  }

  const api = {
    rk4: rk4,
    mayLeonard: mayLeonard,
    lotkaVolterra: lotkaVolterra,
    lvInvariant: lvInvariant,
    cusp: cusp,
    cuspFolds: cuspFolds,
    cuspEquilibria: cuspEquilibria,
    rosenzweigMacArthur: rosenzweigMacArthur,
    competition: competition,
    hastingsPowell: hastingsPowell,
    HP: HP,
    lorenz: lorenz,
    simplexToPlane: simplexToPlane,
    ML_ALPHA: 0.8,
    ML_BETA: 1.3
  };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Flows = api;
})(this);
