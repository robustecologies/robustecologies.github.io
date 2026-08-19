# Changelog

## symplectoR 0.1.0

### Added

- Objective specification layer:
  [`sym_objective()`](https://robustecologies.github.io/symplectoR/reference/sym_objective.md)
  (R closures or compiled function pointers, box constraints with
  logarithmic barrier or reflection, automatic central finite-difference
  gradients) and
  [`sym_compile()`](https://robustecologies.github.io/symplectoR/reference/sym_compile.md)
  (thread-safe compiled objectives enabling every OpenMP path).
- Family A solvers in
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md):
  `"rgd"` (relativistic gradient descent with the `1/sqrt(delta)` step
  bound), its exact parameter presets `"nag"`, `"heavy_ball"`, `"cm"`,
  and `"leapfrog"` (dissipative presymplectic leapfrog with constant,
  Nesterov and mixed damping schedules in the overflow-safe
  finite-difference form).
- Family B solvers: `"slc_poly"` and `"slc_expo"` (fully explicit
  Bregman symplectic leapfrog compositions with gradient momentum
  restarting and temporal looping, both on by default) and `"wibisono"`
  (the p = 2 rate-matching three-sequence scheme with its
  estimate-sequence guarantee).
- Family C solver: `"qhd"` (classical split-step Fourier simulation of
  quantum Hamiltonian descent on up to three dimensions; zeroth order,
  requires only function values; seeded measurement, exact grid
  incumbent, density marginals).
- Multi-start ensembles (`n_starts` in
  [`sym_optim()`](https://robustecologies.github.io/symplectoR/reference/sym_optim.md))
  and parameter sweeps
  ([`sym_sweep()`](https://robustecologies.github.io/symplectoR/reference/sym_sweep.md)),
  OpenMP-parallel for compiled objectives with bit-identical results
  across thread counts.
- Inverse modeling:
  [`sym_inverse()`](https://robustecologies.github.io/symplectoR/reference/sym_inverse.md)
  builds parameter estimation objectives from prediction functions or
  janos `system_spec` trajectory matching (least squares, Gaussian
  likelihood, custom losses; named parameter boxes).
- Ecosystem bridge:
  [`as_incumbent_solver()`](https://robustecologies.github.io/symplectoR/reference/as_incumbent_solver.md)
  adapts any trajectory method to the `optim`-style
  `(par, fn, lower, upper)` local-solver contract.
- Benchmarks with exact ground truth:
  [`sym_benchmark()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark.md)
  (conditioned quadratic factory, Rosenbrock, damped oscillator and
  Nesterov-damping problems with closed-form continuous trajectories)
  and the 12-function non-smooth global suite of Leng et al. (2025) via
  [`sym_benchmark_suite()`](https://robustecologies.github.io/symplectoR/reference/sym_benchmark_suite.md).
- S3 classes `sym_fit`, `sym_ensemble`, `sym_sweep`, `sym_objective`,
  `sym_benchmark`, each with print, summary and plot methods. Single
  views cover traces, iterate paths, phase portraits, Lyapunov energy,
  gradient norms, per-step increments, quantum densities and measured
  samples, start maps, value histograms, stability heatmaps,
  cost-quality frontiers, marginal success profiles, axis profiles and
  analytic trajectories; `type = "dashboard"` composes the views
  available for the object at hand into a single annotated figure.
- Empirical data for the inverse-modeling front end:
  `paramecium_didinium` (Veilleux microcosm predator-prey record) and
  `nicholson_blowfly` (four blowfly cages), each documented from the
  primary literature.
