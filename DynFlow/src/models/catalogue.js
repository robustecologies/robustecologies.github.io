// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab
/* Model catalogue. Every entry is a system written in the formula language
   of src/core/expr.js plus a default scene; `source` names where the model
   and its parameter values were taken from (an RElab package, file and
   line, or the original literature), so that each entry can be traced.
   The numerical claims in `about` are checked by tests/models.test.mjs. */
(function (DF) {
  "use strict";

  const M = [];
  function add(id, name, group, source, about, system, scene) {
    M.push({ id: id, name: name, group: group, source: source, about: about, system: system.trim().replace(/^ +/gm, ""), scene: scene || {} });
  }
  const PI2 = 2 * Math.PI;

  // =================================================================== ecology
  add("lotka-volterra", "Lotka-Volterra predator and prey", "Ecology",
    "janos R/shiny_app.R:404 (lotka_volterra)",
    "Neutral cycles around the coexistence point: H = eaN - m ln N + aP - r ln P is conserved, so every orbit is closed.",
    `N' = r*N - a*N*P
     P' = e*a*N*P - m*P
     param r = 1 [0.2, 2]
     param a = 0.05 [0.01, 0.2]
     param e = 0.4 [0.1, 1]
     param m = 0.4 [0.05, 1.5]
     init N = 20
     init P = 10
     range N = [0, 80]
     range P = [0, 60]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02, overlay: { equations: true } });

  add("rosenzweig-macarthur", "Rosenzweig-MacArthur", "Ecology",
    "janos R/shiny_app.R:417 (rosenzweig)",
    "Paradox of enrichment: the coexistence equilibrium loses stability in a Hopf bifurcation at K = 3.75, and a limit cycle grows with K.",
    `N' = r*N*(1 - N/K) - a*N*P/(1 + a*h*N)
     P' = e*a*N*P/(1 + a*h*N) - m*P
     param r = 1 [0.2, 2]
     param K = 6 [1, 12]
     param a = 1 [0.2, 2]
     param h = 0.4 [0.1, 1]
     param e = 0.6 [0.1, 1]
     param m = 0.3 [0.05, 0.6]
     init N = 5
     init P = 1
     range N = [0, 7]
     range P = [0, 4.5]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.02, overlay: { equations: true } });

  add("hastings-powell", "Hastings-Powell food chain", "Ecology",
    "kaRma R/demo_system.R:87 and janos R/shiny_app.R:532 (Hastings and Powell 1991)",
    "Three-level food chain with type II responses; chaotic 'teacup' attractor at the classic parameters.",
    `X' = X*(1 - X) - a1*X*Y/(1 + b1*X)
     Y' = a1*X*Y/(1 + b1*X) - a2*Y*Z/(1 + b2*Y) - d1*Y
     Z' = a2*Y*Z/(1 + b2*Y) - d2*Z
     param a1 = 5 [3, 6]
     param b1 = 3 [2, 6.2]
     param a2 = 0.1 [0.05, 0.2]
     param b2 = 2 [1, 3]
     param d1 = 0.4 [0.2, 0.6]
     param d2 = 0.01 [0.001, 0.06]
     init X = 0.8
     init Y = 0.2
     init Z = 8
     range X = [0, 1]
     range Y = [0, 0.5]
     range Z = [7, 10.5]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.2 }, dt: 0.05, overlay: { equations: true } });

  add("may-leonard", "May-Leonard cyclic competition", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1241 and HiRsch R/systems.R:469 (May and Leonard 1975)",
    "Rock-paper-scissors competition among three species. With alpha > 1 > beta and alpha + beta > 2 orbits approach a heteroclinic cycle through the single-species states and linger ever longer near each.",
    `N1' = N1*(1 - N1 - alpha*N2 - beta*N3)
     N2' = N2*(1 - beta*N1 - N2 - alpha*N3)
     N3' = N3*(1 - alpha*N1 - beta*N2 - N3)
     param alpha = 1.5 [0, 2.5]
     param beta = 0.7 [0, 2.5]
     init N1 = 0.5
     init N2 = 0.3
     init N3 = 0.2
     range N1 = [0, 1]
     range N2 = [0, 1]
     range N3 = [0, 1]`,
    { view: { type: "flow", projection: "simplex", spawn: "mixed", life: [120, 420] }, n: 1800, dt: 0.05, style: { colorBy: "dominant" } });

  add("competition-lv", "Lotka-Volterra competition", "Ecology",
    "janos R/analysis_phase_portrait.R:207",
    "Two competitors with strong interspecific competition (a12 = a21 = 1.5): the coexistence point is a saddle whose stable manifold separates the two exclusion outcomes.",
    `N1' = r1*N1*(1 - N1/K1 - a12*N2/K1)
     N2' = r2*N2*(1 - N2/K2 - a21*N1/K2)
     param r1 = 1 [0.1, 2]
     param r2 = 1 [0.1, 2]
     param K1 = 100 [50, 150]
     param K2 = 100 [50, 150]
     param a12 = 1.5 [0, 2]
     param a21 = 1.5 [0, 2]
     init N1 = 10
     init N2 = 80
     range N1 = [0, 110]
     range N2 = [0, 110]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.02 });

  add("glv4", "Four-species generalised Lotka-Volterra", "Ecology",
    "janos R/shiny_app.R:491 (glv4)",
    "Competitive community converging to a stable interior equilibrium.",
    `x1' = x1*(1.0 - 1.0*x1 - 0.6*x2 - 0.3*x3 - 0.1*x4)
     x2' = x2*(0.9 - 0.2*x1 - 1.0*x2 - 0.5*x3 - 0.2*x4)
     x3' = x3*(0.8 - 0.1*x1 - 0.3*x2 - 1.0*x3 - 0.4*x4)
     x4' = x4*(0.7 - 0.2*x1 - 0.1*x2 - 0.2*x3 - 1.0*x4)
     init x1 = 0.4
     init x2 = 0.3
     init x3 = 0.2
     init x4 = 0.1
     range x1 = [0, 1]
     range x2 = [0, 1]
     range x3 = [0, 1]
     range x4 = [0, 1]`,
    { view: { type: "timeseries", window: 60, members: 6 }, n: 6, spread: 0.5, dt: 0.05 });

  add("vano-lv4", "Chaotic four-species Lotka-Volterra", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1118 (Vano et al. 2006)",
    "Four competitors with a chaotic attractor, the smallest competitive Lotka-Volterra community known to be chaotic.",
    `x1' = x1*(1 - x1 - 1.09*x2 - 1.52*x3)
     x2' = 0.72*x2*(1 - x2 - 0.44*x3 - 1.36*x4)
     x3' = 1.53*x3*(1 - 2.33*x1 - x3 - 0.47*x4)
     x4' = 1.27*x4*(1 - 1.21*x1 - 0.51*x2 - 0.35*x3 - x4)
     init x1 = 0.3013
     init x2 = 0.4586
     init x3 = 0.1307
     init x4 = 0.3557
     range x1 = [0, 1]
     range x2 = [0, 1]
     range x3 = [0, 0.5]
     range x4 = [0, 1]`,
    { view: { type: "trajectory", axes: ["x1", "x2", "x3"], warmup: 200, rotate: 0.2 }, dt: 0.05 });

  add("act-lv", "Arneodo-Coullet-Tresser Lotka-Volterra", "Ecology",
    "janos vignettes/chaotic-systems.Rmd:1088 (Arneodo, Coullet and Tresser 1980)",
    "Three-species Lotka-Volterra system with a chaotic attractor near mu = 1.5.",
    `N1' = N1*(0.5*(1 - N1) + 0.5*(1 - N2) + 0.1*(1 - N3))
     N2' = N2*(-0.5*(1 - N1) - 0.1*(1 - N2) + 0.1*(1 - N3))
     N3' = N3*(mu*(1 - N1) + 0.1*(1 - N2) + 0.1*(1 - N3))
     param mu = 1.52 [1.3, 1.6]
     init N1 = 0.5
     init N2 = 0.3
     init N3 = 0.2`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.2 }, dt: 0.05 });

  add("huisman-weissing", "Huisman-Weissing resource competition", "Ecology",
    "wadaR R/multispecies_competition.R:232 (Huisman and Weissing 2001, Am. Nat. 157: 488)",
    "Five phytoplankton species compete for three resources with Liebig growth; the surviving community depends sensitively on the initial abundances.",
    `aux q1 = max(R1, 0)
     aux q2 = max(R2, 0)
     aux q3 = max(R3, 0)
     aux mu1 = r*min(q1/(0.20 + q1), q2/(0.25 + q2), q3/(0.15 + q3))
     aux mu2 = r*min(q1/(0.05 + q1), q2/(0.10 + q2), q3/(0.95 + q3))
     aux mu3 = r*min(q1/(1.00 + q1), q2/(0.05 + q2), q3/(0.35 + q3))
     aux mu4 = r*min(q1/(0.05 + q1), q2/(1.00 + q2), q3/(0.10 + q3))
     aux mu5 = r*min(q1/(1.20 + q1), q2/(0.40 + q2), q3/(0.05 + q3))
     N1' = N1*(mu1 - m)
     N2' = N2*(mu2 - m)
     N3' = N3*(mu3 - m)
     N4' = N4*(mu4 - m)
     N5' = N5*(mu5 - m)
     R1' = D*(S - R1) - (0.20*mu1*N1 + 0.10*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.10*mu5*N5)
     R2' = D*(S - R2) - (0.10*mu1*N1 + 0.20*mu2*N2 + 0.10*mu3*N3 + 0.10*mu4*N4 + 0.20*mu5*N5)
     R3' = D*(S - R3) - (0.10*mu1*N1 + 0.10*mu2*N2 + 0.20*mu3*N3 + 0.20*mu4*N4 + 0.10*mu5*N5)
     param r = 1 [0.5, 1.5]
     param m = 0.25 [0.1, 0.4]
     param D = 0.25 [0.1, 0.5]
     param S = 10 [5, 15]
     init N1 = 0.1
     init N2 = 0.1
     init N3 = 0.1
     init N4 = 0.1
     init N5 = 0.1
     init R1 = 10
     init R2 = 10
     init R3 = 10`,
    { view: { type: "timeseries", vars: ["N1", "N2", "N3", "N4", "N5"], window: 800 }, dt: 0.1, keepPositive: true, style: { palette: "relab-qualitative" } });

  add("grazing", "May grazing model", "Ecology",
    "nonautonomeR R/systems.R:1018 (May 1977)",
    "Vegetation under grazing with a type III response: two stable states separated by an unstable one, and hysteresis when the grazing rate c is swept.",
    `V' = r*V*(1 - V/K) - c*V^2/(V^2 + V0^2)
     param r = 1 [0.5, 2]
     param K = 10 [5, 15]
     param c = 2 [1, 3]
     param V0 = 1 [0.5, 2]
     init V = 7
     range V = [0, 10]`,
    { view: { type: "sweep", param: "c", var: "V", from: 1, to: 3, speed: 0.0015 }, dt: 0.05 });

  add("allee", "Strong Allee effect", "Ecology",
    "janos R/analysis_fokker_planck.R:854",
    "Populations below the threshold A decline to extinction; above it they grow to K.",
    `x' = r*x*(x/A - 1)*(1 - x/K)
     param r = 1 [0.2, 2]
     param A = 0.3 [0.05, 0.6]
     param K = 1 [0.7, 1.5]
     init x = 0.5
     range x = [0, 1.2]`,
    { view: { type: "timeseries", window: 20, members: 12 }, n: 12, initMode: "box", dt: 0.02 });

  add("nicholson-bailey", "Nicholson-Bailey host and parasitoid", "Ecology",
    "janos R/shiny_app.R:719 (Nicholson and Bailey 1935)",
    "Host-parasitoid map whose equilibrium is always unstable: oscillations of growing amplitude.",
    `H[n+1] = lambda*H*exp(-a*P)
     P[n+1] = c*H*(1 - exp(-a*P))
     param lambda = 1.5 [1.05, 2]
     param a = 0.02 [0.005, 0.05]
     param c = 1 [0.5, 2]
     init H = 25
     init P = 10
     range H = [0, 200]
     range P = [0, 200]`,
    { view: { type: "trajectory", tail: 60, dim3: false }, style: { pointSize: 4 } });

  add("ricker", "Ricker map", "Ecology",
    "janos R/shiny_app.R:697 (Ricker 1954)",
    "Density-dependent growth of a single population: period doubling from r = 2 and chaos above r of about 2.69.",
    `N[n+1] = N*exp(r*(1 - N/K))
     param r = 2.7 [1.5, 3.5]
     param K = 100 [50, 150]
     init N = 10
     range N = [0, 400]`,
    { view: { type: "orbit", param: "r", var: "N", from: 1.5, to: 3.5 }, style: { alpha: 0.3 } });

  add("seasonal-rm", "Seasonally forced Rosenzweig-MacArthur", "Ecology",
    "nonautonomeR vignettes/pullback-attraction.Rmd:655",
    "Prey growth modulated with period T: the limit cycle entrains at weak forcing and becomes chaotic at eps = 0.8; the stroboscopic section shows the attractor.",
    `x' = x*(1 + eps*sin(2*pi*t/T))*(1 - x) - a*x*y/(b + x)
     y' = y*(a*x/(b + x) - d)
     param a = 1 [0.5, 1.5]
     param b = 0.3 [0.1, 0.5]
     param d = 0.35 [0.2, 0.5]
     param eps = 0.8 [0, 1]
     param T = 10 [5, 20]
     init x = 0.5
     init y = 0.3
     range x = [0, 1]
     range y = [0, 0.8]`,
    { view: { type: "strobe", period: 10, transient: 10 }, n: 600, spread: 0.3, dt: 0.02, style: { pointSize: 1.3, alpha: 0.7 } });

  add("coleman", "Coleman logistic with seasonal carrying capacity", "Ecology",
    "nonautonomeR R/systems.R:1119 and vignettes/coleman-model.Rmd:58",
    "Logistic growth with periodic K(t): all positive orbits converge to one periodic orbit, the pullback attractor.",
    `x' = r*x*(1 - x/(1 + A*sin(w*t)))
     param r = 1 [0.2, 3]
     param A = 0.3 [0, 0.8]
     param w = 0.2 [0.05, 1]
     init x = 0.5
     range x = [0, 2]`,
    { view: { type: "timeseries", window: 60, members: 10 }, n: 10, initMode: "box", dt: 0.02 });

  add("toggle-switch", "Genetic toggle switch", "Ecology",
    "janos vignettes/qualitative-analysis.Rmd:641 (Gardner, Cantor and Collins 2000)",
    "Two mutually repressing genes: two stable states separated by a saddle; noise drives switches between them.",
    `u' = alpha/(1 + v^beta) - u
     v' = alpha/(1 + u^gamma) - v
     param alpha = 3 [1, 6]
     param beta = 2.5 [1, 4]
     param gamma = 2.5 [1, 4]
     init u = 2.5
     init v = 0.5
     range u = [0, 3.5]
     range v = [0, 3.5]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.02 });

  add("rock-paper-scissors", "Rock-paper-scissors replicator", "Ecology",
    "janos R/shiny_app.R:1133 (rps)",
    "Replicator dynamics of the zero-sum rock-paper-scissors game: neutral cycles on the simplex around the mixed equilibrium.",
    `aux f1 = -p2 + p3
     aux f2 = p1 - p3
     aux f3 = -p1 + p2
     aux fbar = p1*f1 + p2*f2 + p3*f3
     p1' = p1*(f1 - fbar)
     p2' = p2*(f2 - fbar)
     p3' = p3*(f3 - fbar)
     init p1 = 0.4
     init p2 = 0.35
     init p3 = 0.25
     range p1 = [0, 1]
     range p2 = [0, 1]
     range p3 = [0, 1]`,
    { view: { type: "flow", projection: "simplex", life: [200, 500] }, n: 1200, dt: 0.02 });

  // ============================================================== delays
  add("mackey-glass", "Mackey-Glass", "Delay equations",
    "janos R/shiny_app.R:733 (Mackey and Glass 1977)",
    "Blood-cell production with a delay of 17 time units: chaotic oscillations.",
    `x' = a*lag(x, tau)/(1 + lag(x, tau)^n) - b*x
     param a = 0.2 [0.1, 0.3]
     param b = 0.1 [0.05, 0.2]
     param n = 10 [4, 12]
     param tau = 17 [2, 30]
     init x = 0.9
     range x = [0.2, 1.5]`,
    { view: { type: "timeseries", window: 600 }, dt: 0.1 });

  add("hutchinson", "Hutchinson delayed logistic", "Delay equations",
    "janos R/shiny_app.R:743 (Hutchinson 1948)",
    "Logistic growth with delayed feedback: the equilibrium K loses stability when r tau exceeds pi/2.",
    `N' = r*N*(1 - lag(N, tau)/K)
     param r = 1.6 [0.5, 2.5]
     param K = 100 [50, 150]
     param tau = 1 [0.2, 2]
     init N = 20
     range N = [0, 320]`,
    { view: { type: "timeseries", window: 40 }, dt: 0.02 });

  add("nicholson-blowflies", "Nicholson blowflies", "Delay equations",
    "janos vignettes/introduction.Rmd:265 and symplectoR R/data.R:123 (Gurney, Blythe and Nisbet 1980)",
    "Delayed recruitment with a hump-shaped birth function: large irregular population cycles.",
    `N' = P*lag(N, tau)*exp(-lag(N, tau)/N0) - delta*N
     param P = 8 [2, 12]
     param N0 = 1 [0.5, 2]
     param delta = 0.175 [0.1, 0.4]
     param tau = 15 [5, 20]
     init N = 3
     range N = [0, 20]`,
    { view: { type: "timeseries", window: 300 }, dt: 0.1 });

  add("delayed-predator-prey", "Delayed predator and prey", "Delay equations",
    "janos vignettes/qualitative-analysis.Rmd:768",
    "Prey self-regulation acts with a delay tau, which destabilises coexistence into cycles in a Hopf bifurcation at tau = 1.437. At tau = 1.5 the prey cycles between about 1.7 and 10; at the value 3 of the source the cycles reach N of about 330, far outside the plot, so tau = 1.5 is the default here.",
    `N' = r*N*(1 - lag(N, tau)/K) - a*N*P
     P' = b*N*P - d*P
     param r = 1.5 [0.5, 2.5]
     param K = 10 [5, 15]
     param a = 0.2 [0.1, 0.4]
     param b = 0.1 [0.05, 0.2]
     param d = 0.5 [0.2, 1]
     param tau = 1.5 [0.5, 5]
     init N = 5
     init P = 2
     range N = [0, 14]
     range P = [0, 12]`,
    { view: { type: "trajectory" }, dt: 0.02 });

  // ================================================================ chaos
  add("lorenz", "Lorenz", "Chaotic flows",
    "janos R/shiny_app.R:506 and tuRbulence R/dynamical_systems.R:61 (Lorenz 1963)",
    "Convection model with the butterfly attractor; largest Lyapunov exponent 0.906 at the classical parameters.",
    `x' = sigma*(y - x)
     y' = x*(rho - z) - y
     z' = x*y - beta*z
     param sigma = 10 [1, 20]
     param rho = 28 [0.5, 50]
     param beta = 2.6666666666666665 [0.5, 4]
     init x = 1
     init y = 1
     init z = 1
     range x = [-22, 22]
     range y = [-28, 28]
     range z = [0, 52]`,
    { view: { type: "trajectory", warmup: 500, rotate: 0.25 }, dt: 0.005, overlay: { equations: true } });

  add("rossler", "Rossler", "Chaotic flows",
    "janos R/shiny_app.R:520 and tuRbulence R/dynamical_systems.R:254 (Rossler 1976)",
    "Band chaos from a single folded band; largest Lyapunov exponent about 0.07 at a = b = 0.2, c = 5.7.",
    `x' = -y - z
     y' = x + a*y
     z' = b + z*(x - c)
     param a = 0.2 [0, 0.4]
     param b = 0.2 [0, 1]
     param c = 5.7 [2, 12]
     init x = 1
     init y = 1
     init z = 1
     range x = [-12, 14]
     range y = [-14, 11]
     range z = [0, 24]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.2 }, dt: 0.01 });

  add("chua", "Chua double scroll", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:167 (Chua, Komuro and Matsumoto 1986)",
    "Electronic circuit with a piecewise-linear diode: the double-scroll attractor.",
    `x' = alpha*(y - x - (m1*x + 0.5*(m0 - m1)*(abs(x + 1) - abs(x - 1))))
     y' = x - y + z
     z' = -beta*y
     param alpha = 15.6 [8, 20]
     param beta = 28 [20, 35]
     param m0 = -1.143 [-1.5, -0.8]
     param m1 = -0.714 [-1, -0.4]
     init x = 0.1
     init y = 0
     init z = 0
     range x = [-2.6, 2.6]
     range y = [-0.5, 0.5]
     range z = [-4, 4]`,
    { view: { type: "trajectory", warmup: 400, rotate: 0.2, axes: ["x", "z", "y"] }, dt: 0.005 });

  add("chen", "Chen", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:204 (Chen and Ueta 1999)",
    "A Lorenz-like system with a double-wing attractor of different topology.",
    `x' = a*(y - x)
     y' = (c - a)*x - x*z + c*y
     z' = x*y - b*z
     param a = 35 [30, 40]
     param b = 3 [1, 5]
     param c = 28 [20, 30]
     init x = -10
     init y = 0
     init z = 37
     range x = [-30, 30]
     range y = [-32, 32]
     range z = [0, 60]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.25 }, dt: 0.001 });

  add("lu", "Lu", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:231 (Lu and Chen 2002)",
    "Intermediate between the Lorenz and Chen attractors.",
    `x' = a*(y - x)
     y' = -x*z + c*y
     z' = x*y - b*z
     param a = 36 [30, 40]
     param b = 3 [1, 5]
     param c = 20 [12, 28]
     init x = 0.1
     init y = 0.2
     init z = 0.3
     range x = [-25, 25]
     range y = [-28, 28]
     range z = [0, 45]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.25 }, dt: 0.002 });

  add("shimizu-morioka", "Shimizu-Morioka", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:269 (Shimizu and Morioka 1980)",
    "Lorenz-type attractor of a laser model.",
    `x' = y
     y' = x - lambda*y - x*z
     z' = -alpha*z + x^2
     param alpha = 0.45 [0.3, 0.6]
     param lambda = 0.75 [0.5, 1]
     init x = 0.1
     init y = 0.1
     init z = 0.1
     range x = [-2, 2]
     range y = [-1.6, 1.6]
     range z = [0, 2.4]`,
    { view: { type: "trajectory", warmup: 400, rotate: 0.2 }, dt: 0.02 });

  add("nose-hoover", "Nose-Hoover (Sprott A)", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:301 (Hoover 1985)",
    "Thermostatted oscillator: a conservative flow in which a chaotic sea coexists with invariant tori.",
    `x' = y
     y' = -x - z*y
     z' = y^2 - a
     param a = 1 [0.5, 2]
     init x = 0
     init y = 5
     init z = 0
     range x = [-4, 4]
     range y = [-5, 5]
     range z = [-4, 4]`,
    { view: { type: "trajectory", rotate: 0.2 }, dt: 0.01 });

  add("sprott-jerk", "Sprott minimal jerk", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:343 (Sprott 1997)",
    "Third-order equation x''' = -a x'' + x'^2 - x, one of the simplest chaotic flows.",
    `x' = y
     y' = z
     z' = -a*z + y^2 - x
     param a = 2.017 [1.9, 2.1]
     init x = 0
     init y = 0
     init z = 1
     range x = [-7, 5]
     range y = [-3, 3]
     range z = [-4, 3]`,
    { view: { type: "trajectory", warmup: 400, rotate: 0.2 }, dt: 0.01 });

  add("thomas", "Thomas cyclically symmetric", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:375 (Thomas 1999)",
    "Symmetric flow driven by sines; small damping b gives a labyrinth of chaotic motion.",
    `x' = sin(y) - b*x
     y' = sin(z) - b*y
     z' = sin(x) - b*z
     param b = 0.18 [0.05, 0.3]
     init x = 2.4
     init y = 2.5
     init z = 2.6
     range x = [-5, 5]
     range y = [-5, 5]
     range z = [-5, 5]`,
    { view: { type: "flow", life: [200, 600], rotate: 0.15 }, n: 1600, dt: 0.05, style: { colorBy: "speed", ramp: "relab-fire" } });

  add("halvorsen", "Halvorsen", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:424",
    "Cyclically symmetric attractor with three lobes.",
    `x' = -a*x - 4*y - 4*z - y^2
     y' = -a*y - 4*z - 4*x - z^2
     z' = -a*z - 4*x - 4*y - x^2
     param a = 1.4 [1.2, 1.6]
     init x = 1
     init y = 0
     init z = 0
     range x = [-12, 8]
     range y = [-12, 8]
     range z = [-12, 8]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.2 }, dt: 0.005 });

  add("aizawa", "Aizawa", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:457",
    "A sphere-like attractor with a tube along its axis.",
    `x' = (z - b)*x - d*y
     y' = d*x + (z - b)*y
     z' = c + a*z - z^3/3 - (x^2 + y^2)*(1 + e*z) + f*z*x^3
     param a = 0.95 [0.7, 1]
     param b = 0.7 [0.5, 0.9]
     param c = 0.6 [0.4, 0.8]
     param d = 3.5 [2, 5]
     param e = 0.25 [0, 0.5]
     param f = 0.1 [0, 0.3]
     init x = 0.1
     init y = 0
     init z = 0
     range x = [-1.6, 1.6]
     range y = [-1.6, 1.6]
     range z = [-0.6, 2]`,
    { view: { type: "flow", life: [150, 500], rotate: 0.2 }, n: 1500, dt: 0.01, style: { colorBy: "speed", ramp: "mako" } });

  add("rabinovich-fabrikant", "Rabinovich-Fabrikant", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:490 (Rabinovich and Fabrikant 1979)",
    "Modulation instability in a non-equilibrium medium; strongly stretched attractor.",
    `x' = y*(z - 1 + x^2) + gamma*x
     y' = x*(3*z + 1 - x^2) + gamma*y
     z' = -2*z*(alpha + x*y)
     param gamma = 0.87 [0.1, 1]
     param alpha = 1.1 [0.9, 1.3]
     init x = -1
     init y = 0
     init z = 0.5
     range x = [-2.5, 2.5]
     range y = [-3, 3]
     range z = [0, 2]`,
    { view: { type: "trajectory", warmup: 400, rotate: 0.2 }, dt: 0.002 });

  add("lorenz-84", "Lorenz-84 atmosphere", "Chaotic flows",
    "tuRbulence R/dynamical_systems.R:447 and nonautonomeR R/systems.R:651 (Lorenz 1984)",
    "Low-order model of the westerlies (X) and a travelling wave (Y, Z), chaotic at F = 8, G = 1.",
    `X' = -Y^2 - Z^2 - a*X + a*F
     Y' = X*Y - b*X*Z - Y + G
     Z' = b*X*Y + X*Z - Z
     param a = 0.25 [0.1, 0.5]
     param b = 4 [2, 6]
     param F = 8 [4, 10]
     param G = 1 [0, 2]
     init X = 1
     init Y = 1
     init Z = 1
     range X = [-1.5, 2.8]
     range Y = [-2.5, 2.8]
     range Z = [-2.6, 2.6]`,
    { view: { type: "trajectory", warmup: 300, rotate: 0.2 }, dt: 0.01 });

  add("lorenz-96", "Lorenz-96, five sites", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:653 (Lorenz 1996)",
    "Ring of five sites with advection, damping and forcing F: spatiotemporal chaos.",
    `x1' = (x2 - x4)*x5 - x1 + F
     x2' = (x3 - x5)*x1 - x2 + F
     x3' = (x4 - x1)*x2 - x3 + F
     x4' = (x5 - x2)*x3 - x4 + F
     x5' = (x1 - x3)*x4 - x5 + F
     param F = 8 [2, 12]
     init x1 = 8.01
     init x2 = 8
     init x3 = 8
     init x4 = 8
     init x5 = 8
     range x1 = [-8, 13]
     range x2 = [-8, 13]
     range x3 = [-8, 13]
     range x4 = [-8, 13]
     range x5 = [-8, 13]`,
    { view: { type: "trajectory", axes: ["x1", "x2", "x3"], warmup: 400, rotate: 0.2 }, dt: 0.005 });

  add("hyperchaotic-rossler", "Hyperchaotic Rossler", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:974 (Rossler 1979)",
    "Four-dimensional flow with two positive Lyapunov exponents.",
    `x' = -y - z
     y' = x + a*y + w
     z' = b + x*z
     w' = -c*z + d*w
     param a = 0.25 [0.2, 0.3]
     param b = 3 [2, 4]
     param c = 0.5 [0.3, 0.7]
     param d = 0.05 [0.02, 0.08]
     init x = -10
     init y = -6
     init z = 0
     init w = 10
     range x = [-60, 40]
     range y = [-40, 40]
     range z = [0, 100]
     range w = [0, 120]`,
    { view: { type: "trajectory", axes: ["x", "y", "w"], warmup: 200, rotate: 0.15 }, dt: 0.005 });

  add("newton-leipnik", "Newton-Leipnik", "Chaotic flows",
    "janos vignettes/chaotic-systems.Rmd:1283 (Leipnik and Newton 1981)",
    "Rigid-body motion with feedback: two coexisting strange attractors, reached from z0 = -0.16 and z0 = -0.18.",
    `x' = -a*x + y + 10*y*z
     y' = -x - a*y + 5*x*z
     z' = alpha*z - 5*x*y
     param a = 0.4 [0.3, 0.5]
     param alpha = 0.175 [0.1, 0.25]
     init x = 0.349
     init y = 0
     init z = -0.16
     range x = [-0.6, 0.6]
     range y = [-0.6, 0.6]
     range z = [-0.6, 0.1]`,
    { view: { type: "trajectory", warmup: 200, rotate: 0.2 }, dt: 0.01 });

  add("charney-devore", "Charney-DeVore three modes", "Chaotic flows",
    "tuRbulence R/charney_devore.R:84 (Charney and DeVore 1979)",
    "Truncated barotropic flow over topography. At F = 4 a zonal state (weak waves) and a blocked state (strong waves) are both stable, separated by a saddle; particles released over the box settle into one regime or the other. Below F of about 3.4 the blocked state is unstable, and below about 2.5 it does not exist.",
    `x' = k*(F - x) - alpha*y*z + beta*y
     y' = -k*y + alpha*x*z - beta*x - delta*z
     z' = -k*z + delta*y
     param F = 4 [0.5, 5]
     param k = 0.1 [0.05, 0.3]
     param alpha = 1 [0.5, 1.5]
     param beta = 0.5 [0.2, 1]
     param delta = 1 [0.5, 1.5]
     init x = 1
     init y = 0.1
     init z = 0.1
     range x = [0, 4.5]
     range y = [-1, 1]
     range z = [-2.5, 2.8]`,
    { view: { type: "flow", life: [240, 720], rotate: 0.15 }, n: 1200, dt: 0.02, style: { colorBy: "speed", ramp: "relab-fire" } });

  // ============================================================== forced
  add("duffing", "Forced Duffing oscillator", "Forced oscillators",
    "janos R/shiny_app.R:547 (duffing)",
    "Double-well oscillator driven with period 2 pi / 1.2: the stroboscopic map reveals a strange attractor.",
    `x' = y
     y' = -delta*y - alpha*x - beta*x^3 + gamma*cos(omega*t)
     param delta = 0.3 [0.1, 0.5]
     param alpha = -1 [-1.5, 1]
     param beta = 1 [0.5, 1.5]
     param gamma = 0.5 [0.2, 0.6]
     param omega = 1.2 [0.8, 1.6]
     init x = 0.1
     init y = 0
     range x = [-2, 2]
     range y = [-1.6, 1.6]`,
    { view: { type: "strobe", period: PI2 / 1.2, transient: 5 }, n: 500, spread: 1, dt: 0.02, style: { pointSize: 1.2, alpha: 0.6 } });

  add("forced-van-der-pol", "Forced Van der Pol", "Forced oscillators",
    "janos vignettes/chaotic-systems.Rmd:532",
    "Relaxation oscillator driven at a frequency near its own: chaos between locking regimes.",
    `x' = y
     y' = mu*(1 - x^2)*y - x + A*sin(omega*t)
     param mu = 3 [0.5, 5]
     param A = 5 [0, 8]
     param omega = 1.788 [1, 3]
     init x = 0.1
     init y = 0
     range x = [-3, 3]
     range y = [-8, 8]`,
    { view: { type: "trajectory" }, dt: 0.005 });

  add("forced-pendulum", "Forced damped pendulum", "Forced oscillators",
    "wadaR R/basins.R:73 and R/wada_detection.R:584",
    "Periodically driven pendulum x'' + gamma x' + sin x = F cos t; at gamma = 0.2, F = 1.66 its basins have the Wada property.",
    `x' = v
     v' = -gamma*v - sin(x) + F*cos(t)
     param gamma = 0.2 [0.05, 0.5]
     param F = 1.66 [0.5, 2.5]
     init x = 0
     init v = 0
     range x = [-3.1416, 3.1416]
     range v = [-4, 4]`,
    { view: { type: "strobe", period: PI2, transient: 10 }, n: 600, spread: 3, dt: 0.02, style: { pointSize: 1.3 } });

  add("driven-oscillator", "Driven damped linear oscillator", "Forced oscillators",
    "janos R/analysis_stroboscopic.R:75",
    "Every orbit converges to the unique periodic response; the stroboscopic points collapse onto one fixed point.",
    `x' = v
     v' = -w0^2*x - 2*zeta*v + F*cos(Om*t)
     param w0 = 1 [0.5, 2]
     param zeta = 0.1 [0.02, 0.5]
     param F = 0.5 [0, 1]
     param Om = 1.3 [0.5, 2]
     init x = 1
     init v = 0
     range x = [-2, 2]
     range v = [-2, 2]`,
    { view: { type: "strobe", period: PI2 / 1.3, transient: 0 }, n: 300, spread: 2, dt: 0.01, style: { fade: 0.02, pointSize: 2 } });

  // ========================================================== oscillators
  add("van-der-pol", "Van der Pol", "Oscillators and excitable media",
    "janos R/shiny_app.R:466 (van der Pol 1926)",
    "Self-sustained oscillation with nonlinear damping; relaxation oscillations for large mu.",
    `x' = y
     y' = mu*(1 - x^2)*y - x
     param mu = 2 [0.1, 8]
     init x = 2
     init y = 0
     range x = [-3, 3]
     range y = [-6, 6]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.01 });

  add("fitzhugh-nagumo", "FitzHugh-Nagumo", "Oscillators and excitable media",
    "janos R/shiny_app.R:453 (FitzHugh 1961; Nagumo et al. 1962)",
    "Slow-fast model of a neuron: excitable at low input I, oscillating above a Hopf threshold.",
    `v' = v - v^3/3 - w + I
     w' = eps*(v + a - b*w)
     param I = 0.5 [0, 2]
     param eps = 0.08 [0.01, 0.3]
     param a = 0.7 [0.3, 1]
     param b = 0.8 [0.3, 1]
     init v = -1
     init w = -0.5
     range v = [-2.5, 2.5]
     range w = [-1, 2]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02 });

  add("brusselator", "Brusselator", "Oscillators and excitable media",
    "janos R/shiny_app.R:437 (Prigogine and Lefever 1968)",
    "Autocatalytic reaction scheme: the equilibrium (A, B/A) undergoes a Hopf bifurcation at B = 1 + A^2.",
    `X' = A - (B + 1)*X + X^2*Y
     Y' = B*X - X^2*Y
     param A = 1 [0.5, 2]
     param B = 3 [0.3, 6]
     init X = 1
     init Y = 1
     range X = [0, 4.5]
     range Y = [0, 5.5]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.01 });

  add("selkov", "Sel'kov glycolysis", "Oscillators and excitable media",
    "janos vignettes/introduction.Rmd:209 (Sel'kov 1968)",
    "Glycolytic oscillations: a Hopf bifurcation in the flux b creates a limit cycle.",
    `x' = -x + a*y + x^2*y
     y' = b - a*y - x^2*y
     param a = 0.1 [0, 0.2]
     param b = 0.5 [0.1, 1.2]
     init x = 0.5
     init y = 0.5
     range x = [0, 3]
     range y = [0, 3]`,
    { view: { type: "phase", seeds: 6 }, dt: 0.02 });

  add("pendulum", "Pendulum", "Oscillators and excitable media",
    "classical mechanics",
    "Frictionless pendulum x'' = -sin x: librations inside the separatrix through the saddles at x = +/-pi, rotations outside.",
    `x' = v
     v' = -sin(x) - c*v
     param c = 0 [0, 0.5]
     init x = 1
     init v = 0
     range x = [-7, 7]
     range v = [-3.2, 3.2]`,
    { view: { type: "phase", seeds: 14 }, dt: 0.01 });

  add("kuramoto-pair", "Two coupled phase oscillators", "Oscillators and excitable media",
    "nonautonomeR R/systems.R:1453 (kuramoto_pair)",
    "Phase difference psi = phi1 - phi2 obeys psi' = dw - 2K sin psi: phase locking when 2K > |dw|, phase slips otherwise.",
    `phi1' = w1 + K*sin(phi2 - phi1)
     phi2' = w2 + K*sin(phi1 - phi2)
     param w1 = 6.911503837897544 [5, 8]
     param w2 = 5.654866776461628 [5, 8]
     param K = 0.5 [0, 1.5]
     init phi1 = 0
     init phi2 = 1
     range phi1 = [0, 60]
     range phi2 = [0, 60]`,
    { view: { type: "timeseries", window: 30 }, dt: 0.01 });

  // ================================================================= maps
  add("logistic", "Logistic map", "Maps",
    "janos R/shiny_app.R:680 and kaRma R/demo_system.R:164 (May 1976)",
    "Period-doubling cascade to chaos; at r = 4 the Lyapunov exponent is ln 2.",
    `x[n+1] = r*x*(1 - x)
     param r = 3.9 [2.5, 4]
     init x = 0.2
     range x = [0, 1]`,
    { view: { type: "orbit", param: "r", from: 2.5, to: 4 }, style: { alpha: 0.25 } });

  add("logistic-cobweb", "Logistic map, cobweb", "Maps",
    "janos vignettes/qualitative-analysis.Rmd:346",
    "Graphical iteration of x[n+1] = r x (1 - x); click to restart from another x.",
    `x[n+1] = r*x*(1 - x)
     param r = 3.7 [2.5, 4]
     init x = 0.2
     range x = [0, 1]`,
    { view: { type: "cobweb", tail: 80 } });

  add("henon", "Henon map", "Maps",
    "janos R/shiny_app.R:688 and nonautonomeR R/systems.R:542 (Henon 1976)",
    "Stretch-and-fold map with a strange attractor; largest Lyapunov exponent 0.419 at a = 1.4, b = 0.3.",
    `x[n+1] = 1 - a*x^2 + y
     y[n+1] = b*x
     param a = 1.4 [0.8, 1.42]
     param b = 0.3 [0, 0.4]
     init x = 0.1
     init y = 0.1
     range x = [-1.5, 1.5]
     range y = [-0.45, 0.45]`,
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "relab-fire" } });

  add("lozi", "Lozi map", "Maps",
    "janos vignettes/chaotic-systems.Rmd:810 (Lozi 1978)",
    "Piecewise-linear analogue of the Henon map with a strange attractor.",
    `x[n+1] = 1 - a*abs(x) + y
     y[n+1] = b*x
     param a = 1.7 [1.2, 1.8]
     param b = 0.5 [0.2, 0.6]
     init x = 0.1
     init y = 0
     range x = [-1.4, 1.4]
     range y = [-0.7, 0.7]`,
    { view: { type: "flow", life: [30, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age" } });

  add("ikeda", "Ikeda map", "Maps",
    "janos vignettes/chaotic-systems.Rmd:842 (Ikeda 1979)",
    "Light in a nonlinear optical ring cavity; spiralling strange attractor.",
    `aux th = kappa - alpha/(1 + x^2 + y^2)
     x[n+1] = a + b*(x*cos(th) - y*sin(th))
     y[n+1] = b*(x*sin(th) + y*cos(th))
     param a = 1 [0.5, 1.5]
     param b = 0.9 [0.6, 0.95]
     param kappa = 0.4 [0, 1]
     param alpha = 6 [4, 8]
     init x = 0.1
     init y = 0.1
     range x = [-0.6, 2.2]
     range y = [-2.4, 1]`,
    { view: { type: "flow", life: [40, 200], dim3: false }, n: 3000, style: { fade: 0.03, pointSize: 1.2, alpha: 0.5, colorBy: "age", ramp: "mako" } });

  add("standard-map", "Chirikov standard map", "Maps",
    "janos R/shiny_app.R:705 (Chirikov 1979)",
    "Area-preserving kicked rotor on the torus: KAM curves break up near K = 0.9716 and a chaotic sea spreads.",
    `p[n+1] = mod(p + K*sin(th), 2*pi)
     th[n+1] = mod(th + p + K*sin(th), 2*pi)
     param K = 1.2 [0, 3]
     init p = 0.5
     init th = 0.5
     range p = [0, 6.2832]
     range th = [0, 6.2832]`,
    { view: { type: "flow", life: [400, 2000], axes: ["th", "p"], dim3: false }, n: 1500, style: { fade: 0, pointSize: 1, alpha: 0.5, colorBy: "member", palette: "relab" } });

  add("zaslavsky", "Zaslavsky random map", "Maps",
    "nonautonomeR R/systems.R:741 (Namenson, Ott and Antonsen 1996)",
    "Dissipative kicked rotor with a random phase c_n drawn afresh at each step and shared by the ensemble: a snapshot attractor that changes shape every step.",
    `aux xn = mod(x + y*(1 - exp(-alpha))/alpha, 2*pi)
     x[n+1] = xn
     y[n+1] = kappa*sin(xn + 2*pi*ucommon(0)) + exp(-alpha)*y
     param alpha = 0.09 [0.02, 0.3]
     param kappa = 0.5 [0.1, 1]
     init x = 1
     init y = 0.1
     range x = [0, 6.2832]
     range y = [-4, 4]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, spread: 0, initMode: "box", style: { fade: 0.6, pointSize: 1.4, alpha: 0.8, colorBy: "solid", palette: "mono-amber" } });

  add("random-baker", "Random baker's map", "Maps",
    "nonautonomeR R/systems.R:801",
    "Baker's map with a random cut c_n shared by the ensemble: a fractal snapshot attractor with known generalised dimensions.",
    `aux c = ucommon(0)
     x[n+1] = ifelse(y <= c, lambda*x, 0.5 + lambda*x)
     y[n+1] = ifelse(y <= c, y/c, (y - c)/(1 - c))
     param lambda = 0.4 [0.1, 0.5]
     init x = 0.3
     init y = 0.6
     range x = [0, 1]
     range y = [0, 1]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 6000, initMode: "box", style: { fade: 0.7, pointSize: 1.3, alpha: 0.8, colorBy: "solid", palette: "mono-ice" } });

  add("stark-circle", "Quasiperiodically forced circle map", "Maps",
    "nonautonomeR R/demo_stark_skew.R:47 (Stark 1999)",
    "Circle diffeomorphism driven by an irrational rotation theta. At a = 0.9 the fibre Lyapunov exponent is negative (about -0.24) and the attractor is the graph of a function of theta; at the package default a = 0.25 it is zero and orbits fill the torus.",
    `th[n+1] = mod(th + Om, 2*pi)
     x[n+1] = mod(x + nu + a*sin(x) + b*sin(th), 2*pi)
     param Om = 3.883222077450933 [0, 6.2832]
     param nu = 0.8676521529893011 [0, 6.2832]
     param a = 0.9 [0, 0.95]
     param b = 0.35 [0, 1]
     init th = 0.7
     init x = 0.2
     range th = [0, 6.2832]
     range x = [0, 6.2832]`,
    { view: { type: "flow", life: [300, 1200], dim3: false }, n: 2000, initMode: "box", style: { fade: 0.02, pointSize: 1.3, alpha: 0.6, colorBy: "var", colorVar: "x", ramp: "blackboard" } });

  add("kaplan-yorke", "Kaplan-Yorke map", "Maps",
    "nonautonomeR R/systems.R:915 (Kaplan and Yorke 1979)",
    "Doubling map driving a contracting coordinate: a strange attractor of Kaplan-Yorke dimension 1 + ln 2 / |ln alpha|. Floating-point doubling discards one bit per step and would reach x = 0 within 53 steps, so a noise of size 1e-9 re-injects the low-order bits.",
    `x[n+1] = mod(2*x + 0.000000001*nrand(), 1)
     y[n+1] = alpha*y + cos(4*pi*x)
     param alpha = 0.2 [0.05, 0.6]
     init x = 0.3678
     init y = 0.6677
     range x = [0, 1]
     range y = [-1.3, 1.3]`,
    { view: { type: "flow", life: [20, 60], dim3: false }, n: 3000, initMode: "box", style: { fade: 0.05, pointSize: 1.2, alpha: 0.5 } });

  add("de-jong", "Peter de Jong attractor", "Maps",
    "Pickover (1990), Computers, Pattern, Chaos and Beauty",
    "Trigonometric map whose attractors are used in generative art; many parameter sets are chaotic.",
    `x[n+1] = sin(a*y) - cos(b*x)
     y[n+1] = sin(c*x) - cos(d*y)
     param a = 1.4 [-3, 3]
     param b = -2.3 [-3, 3]
     param c = 2.4 [-3, 3]
     param d = -2.1 [-3, 3]
     init x = 0
     init y = 0
     range x = [-2.2, 2.2]
     range y = [-2.2, 2.2]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "age", ramp: "relab-fire" } });

  add("clifford", "Clifford attractor", "Maps",
    "Pickover (1990), Computers, Pattern, Chaos and Beauty",
    "Trigonometric map related to the de Jong map, with flowing filamentary attractors.",
    `x[n+1] = sin(a*y) + c*cos(a*x)
     y[n+1] = sin(b*x) + d*cos(b*y)
     param a = -1.4 [-3, 3]
     param b = 1.6 [-3, 3]
     param c = 1 [-3, 3]
     param d = 0.7 [-3, 3]
     init x = 0.1
     init y = 0.1
     range x = [-2.2, 2.2]
     range y = [-2.2, 2.2]`,
    { view: { type: "flow", life: "inf", dim3: false }, n: 4000, initMode: "box", style: { fade: 0.004, pointSize: 0.8, alpha: 0.25, colorBy: "var", colorVar: "y", ramp: "mako" } });

  // =========================================================== stochastic
  add("ornstein-uhlenbeck", "Ornstein-Uhlenbeck process", "Stochastic",
    "janos R/shiny_app.R:765 (ou)",
    "Mean-reverting diffusion with stationary variance sigma^2 / (2 theta) = 0.125.",
    `x' = -theta*(x - mu)
     noise x = sigma
     param theta = 1 [0.1, 3]
     param mu = 0 [-1, 1]
     param sigma = 0.5 [0, 1.5]
     init x = 2
     range x = [-1.5, 2.2]`,
    { view: { type: "density", window: 20 }, n: 3000, dt: 0.01, style: { ramp: "relab-fire" } });

  add("double-well", "Noisy double well", "Stochastic",
    "janos R/shiny_app.R:793 and nonautonomeR vignettes/melancholia-states.Rmd:298",
    "Bistable potential V = x^4/4 - x^2/2 with additive noise: Kramers escapes between the wells at x = -1 and x = 1.",
    `x' = x - x^3
     noise x = sigma
     param sigma = 0.45 [0, 1]
     init x = -1
     range x = [-2, 2]`,
    { view: { type: "density", window: 200 }, n: 3000, dt: 0.02, style: { ramp: "magma" } });

  add("stochastic-resonance", "Stochastic resonance", "Stochastic",
    "janos vignettes/noise-in-dynamical-systems.Rmd:372 (Benzi, Sutera and Vulpiani 1981)",
    "A weak periodic tilt, too small to push the state over the barrier alone, becomes visible in the switching when the noise is tuned.",
    `x' = x - x^3 + A*cos(w*t)
     noise x = sigma
     param A = 0.12 [0, 0.4]
     param w = 0.1 [0.02, 0.5]
     param sigma = 0.35 [0, 1]
     init x = -1
     range x = [-2, 2]`,
    { view: { type: "timeseries", window: 400, members: 1 }, dt: 0.02 });

  add("verhulst-noise", "Noise-induced transition (stochastic Verhulst)", "Stochastic",
    "janos R/shiny_app.R:784 (Horsthemke and Lefever 1984)",
    "Logistic growth with multiplicative noise: the stationary density changes shape (a P-bifurcation) as sigma crosses sqrt(a).",
    `x' = a*x - x^2
     noise x = sigma*x
     param a = 1 [0.2, 2]
     param sigma = 0.8 [0, 2]
     init x = 1
     range x = [0, 3]`,
    { view: { type: "density", window: 60 }, n: 3000, dt: 0.005, keepPositive: true, style: { ramp: "mako" } });

  add("stochastic-lv", "Stochastic Lotka-Volterra", "Stochastic",
    "janos R/shiny_app.R:802 (sde_lv)",
    "Environmental noise pushes orbits off the conserved cycles of the Lotka-Volterra model.",
    `N' = r*N - a*N*P
     P' = e*a*N*P - m*P
     noise N = s*N
     noise P = s*P
     param r = 1 [0.2, 2]
     param a = 0.05 [0.01, 0.2]
     param e = 0.4 [0.1, 1]
     param m = 0.4 [0.05, 1.5]
     param s = 0.05 [0, 0.3]
     init N = 20
     init P = 10
     range N = [0, 80]
     range P = [0, 60]`,
    { view: { type: "density", decay: 0.97 }, n: 3000, dt: 0.01, style: { ramp: "inferno" } });

  add("coherence-resonance", "Coherence resonance", "Stochastic",
    "janos vignettes/noise-in-dynamical-systems.Rmd:390 (Pikovsky and Kurths 1997)",
    "Excitable FitzHugh-Nagumo unit below threshold: noise alone triggers spikes, most regular at an intermediate noise level.",
    `v' = (v - v^3/3 - w)/eps
     w' = v + a
     noise v = sigma
     param eps = 0.05 [0.01, 0.2]
     param a = 1.05 [0.9, 1.3]
     param sigma = 0.5 [0, 2]
     init v = -1
     init w = -0.6
     range v = [-2.5, 2.5]
     range w = [-1.2, 1.2]`,
    { view: { type: "timeseries", vars: ["v"], window: 40 }, dt: 0.002 });

  add("maier-stein", "Maier-Stein", "Stochastic",
    "nonautonomeR R/systems.R:992 (Maier and Stein 1993)",
    "Two attractors at (+/-1, 0) and a saddle at the origin; for alpha different from mu the drift is not a gradient and escape paths bend.",
    `x' = x - x^3 - alpha*x*y^2
     y' = -mu*(1 + x^2)*y
     noise x = sigma
     noise y = sigma
     param alpha = 3 [0, 6]
     param mu = 1 [0.2, 3]
     param sigma = 0.25 [0, 0.6]
     init x = -1
     init y = 0
     range x = [-1.8, 1.8]
     range y = [-1, 1]`,
    { view: { type: "density", decay: 0.96 }, n: 3000, dt: 0.01, style: { ramp: "magma" } });

  add("geometric-brownian", "Geometric Brownian motion", "Stochastic",
    "janos R/shiny_app.R:775 (gbm)",
    "Multiplicative noise with drift mu: the mean grows as exp(mu t) while the median grows as exp((mu - sigma^2/2) t).",
    `S' = mu*S
     noise S = sigma*S
     param mu = 0.08 [-0.2, 0.3]
     param sigma = 0.3 [0, 0.8]
     init S = 100
     range S = [0, 400]`,
    { view: { type: "timeseries", window: 10, members: 20 }, n: 20, dt: 0.002 });

  add("noisy-van-der-pol", "Noisy Van der Pol", "Stochastic",
    "janos vignettes/qualitative-analysis.Rmd:689",
    "Limit cycle blurred by velocity noise: the ensemble diffuses along the cycle and loses its phase.",
    `x' = y
     y' = mu*(1 - x^2)*y - x
     noise y = sigma
     param mu = 1.5 [0.2, 4]
     param sigma = 0.5 [0, 1.5]
     init x = 2
     init y = 0
     range x = [-3, 3]
     range y = [-5, 5]`,
    { view: { type: "density", decay: 0.9 }, n: 3000, spread: 0.01, dt: 0.01, style: { ramp: "relab-fire" } });

  // ================================================= tipping and nonautonomous
  add("r-tipping", "Rate-induced tipping", "Tipping and nonautonomous",
    "normal form of Ashwin, Wieczorek, Vitolo and Cox (2012), Phil. Trans. R. Soc. A 370: 1166",
    "In the frame y = x + lambda of x' = (x + lambda)^2 - 1, a ramp of lambda at rate r adds r to the drift. For r < 1 the state shifts to y = -sqrt(1 - r) and recovers when the ramp ends; for r > 1 no equilibrium exists and the state escapes (R-tipping) if the ramp outlasts the escape time, although every frozen lambda is safe.",
    `y' = y^2 - 1 + r*step(t - t0)*step(t0 + L - t)
     param r = 1.2 [0, 3]
     param t0 = 2 [0, 10]
     param L = 12 [1, 30]
     init y = -1
     range y = [-2, 3]`,
    { view: { type: "timeseries", window: 25, members: 1 }, dt: 0.005, overlay: { equations: true } });

  add("fold-normal-form", "Saddle-node (fold) normal form", "Tipping and nonautonomous",
    "janos vignettes/advanced-dynamics.Rmd:213",
    "x' = mu + x^2: two equilibria for mu < 0 that collide and vanish at mu = 0.",
    `x' = mu + x^2
     param mu = -1 [-1, 0.5]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "mu", from: -1, to: 0.3, speed: 0.00625 }, dt: 0.01 });

  add("cusp", "Cusp catastrophe", "Tipping and nonautonomous",
    "normal form of the cusp catastrophe (Thom 1972; Zeeman 1977)",
    "x' = r + a x - x^3: for a > 0 two stable branches coexist between the folds at r = +/-2 (a/3)^(3/2), so a slow sweep of r jumps and shows hysteresis.",
    `x' = r + a*x - x^3
     param r = 0 [-0.8, 0.8]
     param a = 1 [-0.5, 2]
     init x = -1
     range x = [-1.6, 1.6]`,
    { view: { type: "sweep", param: "r", from: -0.8, to: 0.8, speed: 0.005 }, dt: 0.01, overlay: { equations: true } });

  add("pitchfork", "Pitchfork normal form", "Tipping and nonautonomous",
    "janos R/analysis_bifurcation_sweep.R:80",
    "x' = a x - x^3: the origin splits into two symmetric stable states at a = 0.",
    `x' = a*x - x^3
     param a = 1 [-1, 1.5]
     init x = 0.05
     range x = [-1.4, 1.4]`,
    { view: { type: "sweep", param: "a", from: -1, to: 1.5, speed: 0.00625 }, dt: 0.01 });

  add("hopf-normal-form", "Hopf normal form", "Tipping and nonautonomous",
    "janos vignettes/advanced-dynamics.Rmd:241",
    "Supercritical Hopf bifurcation: for mu > 0 a stable limit cycle of radius sqrt(mu) surrounds the origin.",
    `x' = mu*x - y - x*(x^2 + y^2)
     y' = x + mu*y - y*(x^2 + y^2)
     param mu = 0.5 [-0.5, 1]
     init x = 0.1
     init y = 0.1
     range x = [-1.3, 1.3]
     range y = [-1.3, 1.3]`,
    { view: { type: "phase", seeds: 8 }, dt: 0.02 });

  add("bogdanov-takens", "Bogdanov-Takens normal form", "Tipping and nonautonomous",
    "normal form (Kuznetsov 2004, Elements of Applied Bifurcation Theory, ch. 8)",
    "Codimension-two unfolding with fold, Hopf and homoclinic bifurcation curves meeting at the origin of (b1, b2).",
    `x' = y
     y' = b1 + b2*x + x^2 + x*y
     param b1 = -0.1 [-0.5, 0.2]
     param b2 = 0.2 [-0.5, 0.5]
     init x = 0
     init y = 0.1
     range x = [-1, 1]
     range y = [-1, 1]`,
    { view: { type: "phase", seeds: 10 }, dt: 0.01 });

  add("stommel", "Stommel two-box ocean", "Tipping and nonautonomous",
    "tuRbulence R/stommel.R:73 (Stommel 1961)",
    "Thermohaline circulation with temperature T and salinity S: bistability between strong and weak overturning, with hysteresis in the freshwater forcing eta2.",
    `T' = eta1 - T*(1 + abs(T - S))
     S' = eta2 - S*(eta3 + abs(T - S))
     param eta1 = 3 [2, 4]
     param eta2 = 1 [0.5, 1.5]
     param eta3 = 0.3 [0.1, 0.6]
     init T = 2
     init S = 1
     range T = [0, 3.5]
     range S = [0, 3.5]`,
    { view: { type: "sweep", param: "eta2", var: "S", from: 0.5, to: 1.5, speed: 0.002 }, dt: 0.02 });

  add("lorenz84-forced", "Lorenz-84 under climate change", "Tipping and nonautonomous",
    "nonautonomeR R/systems.R:651 and vignettes/ergodicity.Rmd:174 (J\u00e1nosi, T\u00e9l and co-authors)",
    "Annual forcing F(t) = F0 + 2 sin(2 pi t / 73) with F0 ramping down: a cloud of 1500 members, all driven alike, traces a snapshot attractor that deforms as the climate changes.",
    `X' = -Y^2 - Z^2 - a*X + a*(F0 + AF*sin(2*pi*t/73))
     Y' = X*Y - b*X*Z - Y + G
     Z' = b*X*Y + X*Z - Z
     param a = 0.25 [0.1, 0.5]
     param b = 4 [2, 6]
     param G = 1 [0, 2]
     param F0 = 9.5 [6, 10]
     param AF = 2 [0, 3]
     init X = 1
     init Y = 0
     init Z = 0
     range X = [-1.5, 3.2]
     range Y = [-3, 3]
     range Z = [-3, 3]`,
    { view: { type: "flow", life: "inf", rotate: 0.1 }, n: 1500, spread: 1.2, initMode: "ball", dt: 0.02, perturbations: [{ kind: "ramp", param: "F0", rate: -0.000274, t0: 0, span: -2 }], style: { fade: 0.25, colorBy: "speed", ramp: "relab-fire" }, overlay: { readout: true } });

  add("lorenz-drift", "Lorenz with drifting rho", "Tipping and nonautonomous",
    "nonautonomeR vignettes/pullback-scenes.Rmd:145",
    "Rayleigh number rho(t) = 28 + 0.02 t: the attractor inflates slowly while the ensemble spreads over it.",
    `x' = sigma*(y - x)
     y' = x*(rho - z) - y
     z' = x*y - beta*z
     param sigma = 10 [5, 15]
     param rho = 28 [20, 40]
     param beta = 2.6666666666666665 [2, 3]
     init x = 1
     init y = 1
     init z = 20
     range x = [-26, 26]
     range y = [-34, 34]
     range z = [0, 62]`,
    { view: { type: "flow", life: "inf", rotate: 0.12 }, n: 1500, initMode: "ball", spread: 6, dt: 0.005, perturbations: [{ kind: "ramp", param: "rho", rate: 0.02, t0: 0, span: 12 }], style: { fade: 0.2 }, overlay: { readout: true } });

  add("duffing-drift", "Duffing with drifting forcing", "Tipping and nonautonomous",
    "nonautonomeR R/systems.R:963 (Janosi and Tel 2024)",
    "Forcing amplitude eps(t) = 0.4 + 0.00045 t: the stroboscopic cloud of an ensemble is a snapshot attractor whose shape follows the drift.",
    `x' = v
     v' = x - x^3 - 2*beta*v + eps*cos(omega*t)
     param beta = 0.2 [0.1, 0.4]
     param eps = 0.4 [0.2, 0.6]
     param omega = 1 [0.8, 1.2]
     init x = 0.5
     init v = 0
     range x = [-2, 2]
     range v = [-1.6, 1.6]`,
    { view: { type: "strobe", period: PI2, transient: 3 }, n: 1500, spread: 1.5, initMode: "ball", dt: 2 * Math.PI / 128, perturbations: [{ kind: "ramp", param: "eps", rate: 0.00045, t0: 0, span: 0.3 }], style: { fade: 0.35, pointSize: 1.6, alpha: 0.8 }, overlay: { readout: true } });

  add("tilted-well", "Double well with a slow tilt", "Tipping and nonautonomous",
    "nonautonomeR vignettes/melancholia-states.Rmd:524",
    "The tilt lambda is ramped slowly: the occupied well loses stability at a fold, and the ensemble switches, earlier when noise is present.",
    `x' = x - x^3 + lambda
     noise x = sigma
     param lambda = -0.6 [-0.6, 0.6]
     param sigma = 0.15 [0, 0.5]
     init x = -1.2
     range x = [-1.8, 1.8]`,
    { view: { type: "density", window: 120 }, n: 3000, dt: 0.02, perturbations: [{ kind: "ramp", param: "lambda", rate: 0.01, t0: 0, span: 1.2 }], style: { ramp: "magma" }, overlay: { readout: true } });

  // ============================================================ epidemics
  add("sir", "SIR epidemic", "Epidemics",
    "janos R/shiny_app.R:477 (Kermack and McKendrick 1927)",
    "Closed epidemic with basic reproduction number R0 = beta / gamma = 3.",
    `S' = -beta*S*I/N0
     I' = beta*S*I/N0 - gamma*I
     R' = gamma*I
     param beta = 0.3 [0.05, 1]
     param gamma = 0.1 [0.02, 0.5]
     param N0 = 1000 [100, 2000]
     init S = 999
     init I = 1
     init R = 0
     range S = [0, 1000]
     range I = [0, 1000]
     range R = [0, 1000]`,
    { view: { type: "timeseries", window: 160 }, dt: 0.05 });

  add("seasonal-sir", "Seasonally forced SIR", "Epidemics",
    "janos vignettes/chaotic-systems.Rmd:1204 (Olsen and Schaffer 1990)",
    "Measles-like epidemic with seasonal transmission: irregular outbreaks from a smooth forcing.",
    `S' = mu - beta0*(1 + beta1*cos(2*pi*t))*S*I - mu*S
     I' = beta0*(1 + beta1*cos(2*pi*t))*S*I - (gam + mu)*I + eps
     param mu = 0.02 [0.01, 0.04]
     param beta0 = 1800 [1000, 2500]
     param beta1 = 0.08 [0, 0.3]
     param gam = 100 [50, 150]
     param eps = 0.000001 [0, 0.00001]
     init S = 0.065
     init I = 0.0002
     range S = [0.045, 0.085]
     range I = [0, 0.002]`,
    { view: { type: "timeseries", window: 20, vars: ["I"] }, dt: 0.0005 });

  DF.CATALOGUE = M;
  DF.catalogueGroups = function () {
    const seen = [];
    M.forEach(function (m) { if (seen.indexOf(m.group) < 0) seen.push(m.group); });
    return seen;
  };
  // Scene for a catalogue entry: the entry's defaults with its system text.
  DF.sceneFor = function (id, overrides) {
    const m = M.find(function (e) { return e.id === id; });
    if (!m) throw new Error("No model '" + id + "' in the catalogue");
    const s = JSON.parse(JSON.stringify(m.scene));
    s.name = m.name; s.model = m.id; s.system = m.system;
    s.overlay = Object.assign({ title: m.name }, s.overlay || {});
    return Object.assign(s, overrides || {});
  };
})(globalThis.RElabFlow = globalThis.RElabFlow || {});
