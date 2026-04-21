# Package index

## Model specification

Unified formula interface to every dynamical system family.

- [`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md)
  : Define a dynamical system model
- [`dynamical_system()`](https://robustecologies.github.io/janos/reference/dynamical_system.md)
  : Specify a dynamical system for Lyapunov analysis
- [`observe()`](https://robustecologies.github.io/janos/reference/observe.md)
  : Add observation noise to a simulation
- [`equilibrium()`](https://robustecologies.github.io/janos/reference/equilibrium.md)
  : Find an equilibrium of an ODE system

## Simulation dispatch

Top-level entry points that hand the model to a solver or an ensemble.

- [`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md)
  : Simulate a dynamical system
- [`ensemble_sim()`](https://robustecologies.github.io/janos/reference/ensemble_sim.md)
  : Simulate an ensemble of independent replicates

## Noise processes

Drop-in noise specifications for SDE and jump-diffusion solvers.

- [`correlated_noise()`](https://robustecologies.github.io/janos/reference/correlated_noise.md)
  : Correlated Gaussian noise specification
- [`levy_noise()`](https://robustecologies.github.io/janos/reference/levy_noise.md)
  : Levy alpha-stable noise specification
- [`fbm_noise()`](https://robustecologies.github.io/janos/reference/fbm_noise.md)
  : Fractional Brownian motion noise specification
- [`colored_noise()`](https://robustecologies.github.io/janos/reference/colored_noise.md)
  : Colored (1/f^beta) noise specification

## Solvers - deterministic

ODE, DDE and map integrators; stiff and non-stiff.

- [`solver_rk4()`](https://robustecologies.github.io/janos/reference/solver_rk4.md)
  : Fixed-step fourth-order Runge-Kutta solver
- [`solver_rk45()`](https://robustecologies.github.io/janos/reference/solver_rk45.md)
  : Adaptive Dormand-Prince RK4(5) solver
- [`solver_rosenbrock()`](https://robustecologies.github.io/janos/reference/solver_rosenbrock.md)
  : Rosenbrock (Rodas3) implicit solver for stiff systems
- [`solver_map()`](https://robustecologies.github.io/janos/reference/solver_map.md)
  : Discrete map solver
- [`solver_dde()`](https://robustecologies.github.io/janos/reference/solver_dde.md)
  : DDE solver with fixed-step RK4 and history buffer

## Solvers - SDE and jump-diffusion

Weak and strong Ito schemes plus Poisson-driven jump channels.

- [`solver_euler_maruyama()`](https://robustecologies.github.io/janos/reference/solver_euler_maruyama.md)
  : Euler-Maruyama SDE solver
- [`solver_milstein()`](https://robustecologies.github.io/janos/reference/solver_milstein.md)
  : Milstein SDE solver
- [`solver_jump_diffusion()`](https://robustecologies.github.io/janos/reference/solver_jump_diffusion.md)
  : Jump-diffusion SDE solver

## Solvers - CTMC

Exact and approximate stochastic simulation of reaction networks.

- [`solver_ssa_direct()`](https://robustecologies.github.io/janos/reference/solver_ssa_direct.md)
  : Gillespie direct method SSA solver
- [`solver_ssa_nrm()`](https://robustecologies.github.io/janos/reference/solver_ssa_nrm.md)
  : Next-reaction method SSA solver
- [`solver_ssa_mnrm()`](https://robustecologies.github.io/janos/reference/solver_ssa_mnrm.md)
  : Modified next-reaction method SSA solver
- [`solver_tau_leap()`](https://robustecologies.github.io/janos/reference/solver_tau_leap.md)
  : Adaptive tau-leap solver
- [`solver_tau_leap_midpoint()`](https://robustecologies.github.io/janos/reference/solver_tau_leap_midpoint.md)
  : Midpoint tau-leap solver
- [`solver_hybrid()`](https://robustecologies.github.io/janos/reference/solver_hybrid.md)
  : Hybrid SSA/CLE solver

## Solvers - spatial and switched

PDE, RDME and piecewise deterministic Markov processes.

- [`solver_mol()`](https://robustecologies.github.io/janos/reference/solver_mol.md)
  : Method of lines (MOL) solver for PDE systems
- [`solver_mol2d()`](https://robustecologies.github.io/janos/reference/solver_mol2d.md)
  : Specify a 2D method-of-lines solver
- [`solver_rdme()`](https://robustecologies.github.io/janos/reference/solver_rdme.md)
  : RDME solver (reaction-diffusion master equation)
- [`solver_pdmp()`](https://robustecologies.github.io/janos/reference/solver_pdmp.md)
  : PDMP solver (piecewise deterministic Markov process)

## Graph constructors

Adjacency builders for graph-structured RDME and metapopulations.

- [`lattice_graph()`](https://robustecologies.github.io/janos/reference/lattice_graph.md)
  : Create a 2D rectangular lattice adjacency matrix
- [`ring_graph()`](https://robustecologies.github.io/janos/reference/ring_graph.md)
  : Create a ring (cycle) graph adjacency matrix
- [`star_graph()`](https://robustecologies.github.io/janos/reference/star_graph.md)
  : Create a star graph adjacency matrix
- [`complete_graph()`](https://robustecologies.github.io/janos/reference/complete_graph.md)
  : Create a complete graph adjacency matrix
- [`random_graph()`](https://robustecologies.github.io/janos/reference/random_graph.md)
  : Create an Erdos-Renyi random graph adjacency matrix

## Qualitative portraits

Vector fields, nullclines, and family-specific phase portraits.

- [`phase_portrait()`](https://robustecologies.github.io/janos/reference/phase_portrait.md)
  : Compute phase portrait of an ODE system
- [`map_portrait()`](https://robustecologies.github.io/janos/reference/map_portrait.md)
  : Compute portrait of a discrete-time map
- [`sde_portrait()`](https://robustecologies.github.io/janos/reference/sde_portrait.md)
  : Compute stochastic phase portrait of an SDE system
- [`dde_portrait()`](https://robustecologies.github.io/janos/reference/dde_portrait.md)
  : Compute phase portrait of a delay differential equation system
- [`pdmp_portrait()`](https://robustecologies.github.io/janos/reference/pdmp_portrait.md)
  : Compute portrait of a piecewise deterministic Markov process

## Lyapunov stack

Unified advisor plus the eight construction algorithms.

- [`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
  : Unified Lyapunov analysis for a model_spec
- [`lyapunov_advisor()`](https://robustecologies.github.io/janos/reference/lyapunov_advisor.md)
  : Advise on feasible Lyapunov techniques for a model
- [`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
  : Construct a Lyapunov function for a dynamical system
- [`lyapunov_discrete()`](https://robustecologies.github.io/janos/reference/lyapunov_discrete.md)
  : Discrete Lyapunov function for a map
- [`lyapunov_stochastic()`](https://robustecologies.github.io/janos/reference/lyapunov_stochastic.md)
  : Stochastic Lyapunov function for an SDE
- [`lyapunov_krasovskii()`](https://robustecologies.github.io/janos/reference/lyapunov_krasovskii.md)
  : Lyapunov-Krasovskii functional for a DDE
- [`lyapunov_pdmp()`](https://robustecologies.github.io/janos/reference/lyapunov_pdmp.md)
  : Common-quadratic Lyapunov function for a PDMP
- [`lyapunov_foster()`](https://robustecologies.github.io/janos/reference/lyapunov_foster.md)
  : Foster-Lyapunov drift function for a CTMC
- [`lyapunov_functional()`](https://robustecologies.github.io/janos/reference/lyapunov_functional.md)
  : Energy Lyapunov functional for a reaction-diffusion PDE

## Fokker-Planck and rare events

Stationary densities, quasi-potentials, Kramers rates and QSD.

- [`fp_stationary()`](https://robustecologies.github.io/janos/reference/fp_stationary.md)
  : Compute the stationary density of the Fokker-Planck equation
- [`fp_stationary_2d()`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md)
  : Compute the 2D stationary density of the Fokker-Planck equation
- [`fp_potential()`](https://robustecologies.github.io/janos/reference/fp_potential.md)
  : Reconstruct the drift potential and quasi-potential
- [`fp_kramers_rate()`](https://robustecologies.github.io/janos/reference/fp_kramers_rate.md)
  : Compute Kramers escape rates for noise-induced transitions
- [`estimate_qsd()`](https://robustecologies.github.io/janos/reference/estimate_qsd.md)
  : Estimate the quasi-stationary distribution via Fleming-Viot
  particles
- [`estimate_extinction()`](https://robustecologies.github.io/janos/reference/estimate_extinction.md)
  : Estimate rare event probability via importance sampling
- [`density_landscape_2d()`](https://robustecologies.github.io/janos/reference/density_landscape_2d.md)
  : Compute an empirical 2D probability density landscape from ensemble
  data

## Sensitivity, continuation, spectral

Parameter gradients, bifurcation tracking and power spectra.

- [`continuation()`](https://robustecologies.github.io/janos/reference/continuation.md)
  : Numerical bifurcation continuation
- [`bifurcation_sweep()`](https://robustecologies.github.io/janos/reference/bifurcation_sweep.md)
  : Compute a bifurcation diagram by parameter sweep
- [`adjoint_sensitivity()`](https://robustecologies.github.io/janos/reference/adjoint_sensitivity.md)
  : Adjoint sensitivity analysis for ODE models
- [`spectral_analysis()`](https://robustecologies.github.io/janos/reference/spectral_analysis.md)
  : Power spectral density estimation
- [`mlmc_estimate()`](https://robustecologies.github.io/janos/reference/mlmc_estimate.md)
  : Estimate expectations of stochastic models via multi-level Monte
  Carlo

## Chaos diagnostics

Numerical Lyapunov spectrum, fractal dimension, Poincare sections,
Gottwald-Melbourne chaos test and attractor-level bifurcation diagrams.

- [`lyapunov_spectrum()`](https://robustecologies.github.io/janos/reference/lyapunov_spectrum.md)
  : Lyapunov spectrum via QR renormalisation
- [`correlation_dimension()`](https://robustecologies.github.io/janos/reference/correlation_dimension.md)
  : Correlation dimension
- [`poincare_section()`](https://robustecologies.github.io/janos/reference/poincare_section.md)
  : Poincare section
- [`zero_one_test()`](https://robustecologies.github.io/janos/reference/zero_one_test.md)
  : 0-1 test for chaos
- [`bifurcation_diagram()`](https://robustecologies.github.io/janos/reference/bifurcation_diagram.md)
  : Bifurcation diagram

## Miscellaneous

Demos and easter eggs.

- [`VonNeumanns_elephant()`](https://robustecologies.github.io/janos/reference/VonNeumanns_elephant.md)
  : With four parameters I can fit an elephant, and with five I can make
  him wiggle his trunk

## S3 methods

Print, summary and plot methods attached to the 24 exported S3 classes.
Users normally do not call these directly; they are listed here for
completeness and cross-reference.

- [`print(`*`<continuation_result>`*`)`](https://robustecologies.github.io/janos/reference/print.continuation_result.md)
  : Print method for continuation_result objects
- [`print(`*`<dde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/print.dde_portrait.md)
  : Print method for dde_portrait objects
- [`print(`*`<dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/print.dyn_sim.md)
  : Concise textual report for a simulated dynamical system
- [`print(`*`<dynamical_system>`*`)`](https://robustecologies.github.io/janos/reference/print.dynamical_system.md)
  : Print method for a dynamical_system object
- [`print(`*`<ensemble_sim>`*`)`](https://robustecologies.github.io/janos/reference/print.ensemble_sim.md)
  : Print method for ensemble_sim objects
- [`print(`*`<equilibrium_point>`*`)`](https://robustecologies.github.io/janos/reference/print.equilibrium_point.md)
  : Print method for equilibrium_point objects
- [`print(`*`<fp_kramers>`*`)`](https://robustecologies.github.io/janos/reference/print.fp_kramers.md)
  : Print method for fp_kramers objects
- [`print(`*`<fp_potential>`*`)`](https://robustecologies.github.io/janos/reference/print.fp_potential.md)
  : Print method for fp_potential objects
- [`print(`*`<fp_stationary>`*`)`](https://robustecologies.github.io/janos/reference/print.fp_stationary.md)
  : Print method for fp_stationary objects
- [`print(`*`<fp_stationary_2d>`*`)`](https://robustecologies.github.io/janos/reference/print.fp_stationary_2d.md)
  : Print method for fp_stationary_2d objects
- [`print(`*`<lyapunov_advisor>`*`)`](https://robustecologies.github.io/janos/reference/print.lyapunov_advisor.md)
  : Print method for a lyapunov_advisor object
- [`print(`*`<lyapunov_function>`*`)`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md)
  : Print method for a lyapunov_function object
- [`print(`*`<lyapunov_report>`*`)`](https://robustecologies.github.io/janos/reference/print.lyapunov_report.md)
  : Print method for a lyapunov_report object
- [`print(`*`<map_portrait>`*`)`](https://robustecologies.github.io/janos/reference/print.map_portrait.md)
  : Print method for map_portrait objects
- [`print(`*`<mlmc_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.mlmc_estimate.md)
  : Print method for mlmc_estimate objects
- [`print(`*`<model_spec>`*`)`](https://robustecologies.github.io/janos/reference/print.model_spec.md)
  : Print method for model_spec objects
- [`print(`*`<noise_spec>`*`)`](https://robustecologies.github.io/janos/reference/print.noise_spec.md)
  : Print method for noise_spec objects
- [`print(`*`<observed_dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/print.observed_dyn_sim.md)
  : Print method for observed_dyn_sim objects
- [`print(`*`<pdmp_portrait>`*`)`](https://robustecologies.github.io/janos/reference/print.pdmp_portrait.md)
  : Print method for pdmp_portrait objects
- [`print(`*`<phase_portrait>`*`)`](https://robustecologies.github.io/janos/reference/print.phase_portrait.md)
  : Print method for phase_portrait objects
- [`print(`*`<qsd_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.qsd_estimate.md)
  : Print method for qsd_estimate objects
- [`print(`*`<rare_event_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.rare_event_estimate.md)
  : Print method for rare_event_estimate objects
- [`print(`*`<sde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/print.sde_portrait.md)
  : Print method for sde_portrait objects
- [`print(`*`<sensitivity_result>`*`)`](https://robustecologies.github.io/janos/reference/print.sensitivity_result.md)
  : Print method for sensitivity_result objects
- [`print(`*`<solver_spec>`*`)`](https://robustecologies.github.io/janos/reference/print.solver_spec.md)
  : Print method for solver_spec objects
- [`print(`*`<spectral_analysis>`*`)`](https://robustecologies.github.io/janos/reference/print.spectral_analysis.md)
  : Print method for spectral_analysis objects
- [`print(`*`<summary.dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.dyn_sim.md)
  : Print method for summary.dyn_sim objects
- [`print(`*`<summary.mlmc_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.mlmc_estimate.md)
  : Print method for summary.mlmc_estimate objects
- [`print(`*`<summary.observed_dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.observed_dyn_sim.md)
  : Print method for summary.observed_dyn_sim objects
- [`print(`*`<summary.qsd_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.qsd_estimate.md)
  : Print method for summary.qsd_estimate objects
- [`print(`*`<summary.rare_event_estimate>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.rare_event_estimate.md)
  : Print method for summary.rare_event_estimate objects
- [`print(`*`<summary.spectral_analysis>`*`)`](https://robustecologies.github.io/janos/reference/print.summary.spectral_analysis.md)
  : Print method for summary.spectral_analysis objects
- [`print(`*`<zero_one_test>`*`)`](https://robustecologies.github.io/janos/reference/print.zero_one_test.md)
  : Print method for zero_one_test
- [`summary(`*`<continuation_result>`*`)`](https://robustecologies.github.io/janos/reference/summary.continuation_result.md)
  : Summary method for continuation_result objects
- [`summary(`*`<dde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/summary.dde_portrait.md)
  : Summary method for dde_portrait objects
- [`summary(`*`<dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/summary.dyn_sim.md)
  : Numerical summary of a simulated dynamical system
- [`summary(`*`<ensemble_sim>`*`)`](https://robustecologies.github.io/janos/reference/summary.ensemble_sim.md)
  : Summary method for ensemble_sim objects
- [`summary(`*`<equilibrium_point>`*`)`](https://robustecologies.github.io/janos/reference/summary.equilibrium_point.md)
  : Summary method for equilibrium_point objects
- [`summary(`*`<fp_kramers>`*`)`](https://robustecologies.github.io/janos/reference/summary.fp_kramers.md)
  : Summary method for fp_kramers objects
- [`summary(`*`<fp_potential>`*`)`](https://robustecologies.github.io/janos/reference/summary.fp_potential.md)
  : Summary method for fp_potential objects
- [`summary(`*`<fp_stationary>`*`)`](https://robustecologies.github.io/janos/reference/summary.fp_stationary.md)
  : Summary method for fp_stationary objects
- [`summary(`*`<fp_stationary_2d>`*`)`](https://robustecologies.github.io/janos/reference/summary.fp_stationary_2d.md)
  : Summary method for fp_stationary_2d objects
- [`summary(`*`<lyapunov_advisor>`*`)`](https://robustecologies.github.io/janos/reference/summary.lyapunov_advisor.md)
  : Summary method for a lyapunov_advisor object
- [`summary(`*`<lyapunov_function>`*`)`](https://robustecologies.github.io/janos/reference/summary.lyapunov_function.md)
  : Summary method for a lyapunov_function object
- [`summary(`*`<lyapunov_report>`*`)`](https://robustecologies.github.io/janos/reference/summary.lyapunov_report.md)
  : Summary method for a lyapunov_report object
- [`summary(`*`<map_portrait>`*`)`](https://robustecologies.github.io/janos/reference/summary.map_portrait.md)
  : Summary method for map_portrait objects
- [`summary(`*`<mlmc_estimate>`*`)`](https://robustecologies.github.io/janos/reference/summary.mlmc_estimate.md)
  : Summary method for mlmc_estimate objects
- [`summary(`*`<observed_dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/summary.observed_dyn_sim.md)
  : Summary method for observed_dyn_sim objects
- [`summary(`*`<pdmp_portrait>`*`)`](https://robustecologies.github.io/janos/reference/summary.pdmp_portrait.md)
  : Summary method for pdmp_portrait objects
- [`summary(`*`<phase_portrait>`*`)`](https://robustecologies.github.io/janos/reference/summary.phase_portrait.md)
  : Summary method for phase_portrait objects
- [`summary(`*`<qsd_estimate>`*`)`](https://robustecologies.github.io/janos/reference/summary.qsd_estimate.md)
  : Summary method for qsd_estimate objects
- [`summary(`*`<rare_event_estimate>`*`)`](https://robustecologies.github.io/janos/reference/summary.rare_event_estimate.md)
  : Summary method for rare_event_estimate objects
- [`summary(`*`<sde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/summary.sde_portrait.md)
  : Summary method for sde_portrait objects
- [`summary(`*`<sensitivity_result>`*`)`](https://robustecologies.github.io/janos/reference/summary.sensitivity_result.md)
  : Summary method for sensitivity_result objects
- [`summary(`*`<spectral_analysis>`*`)`](https://robustecologies.github.io/janos/reference/summary.spectral_analysis.md)
  : Summary method for spectral_analysis objects
- [`plot(`*`<bifurcation_sweep>`*`)`](https://robustecologies.github.io/janos/reference/plot.bifurcation_sweep.md)
  : Plot method for bifurcation_sweep objects
- [`plot(`*`<continuation_result>`*`)`](https://robustecologies.github.io/janos/reference/plot.continuation_result.md)
  : Plot method for continuation_result objects
- [`plot(`*`<dde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/plot.dde_portrait.md)
  : Plot method for dde_portrait objects
- [`plot(`*`<density_landscape_2d>`*`)`](https://robustecologies.github.io/janos/reference/plot.density_landscape_2d.md)
  : Plot method for density_landscape_2d objects
- [`plot(`*`<dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/plot.dyn_sim.md)
  : Visualise a simulated dynamical system
- [`plot(`*`<ensemble_sim>`*`)`](https://robustecologies.github.io/janos/reference/plot.ensemble_sim.md)
  : Plot method for ensemble_sim objects
- [`plot(`*`<fp_kramers>`*`)`](https://robustecologies.github.io/janos/reference/plot.fp_kramers.md)
  : Plot method for fp_kramers objects
- [`plot(`*`<fp_potential>`*`)`](https://robustecologies.github.io/janos/reference/plot.fp_potential.md)
  : Plot method for fp_potential objects
- [`plot(`*`<fp_stationary>`*`)`](https://robustecologies.github.io/janos/reference/plot.fp_stationary.md)
  : Plot method for fp_stationary objects
- [`plot(`*`<fp_stationary_2d>`*`)`](https://robustecologies.github.io/janos/reference/plot.fp_stationary_2d.md)
  : Plot method for fp_stationary_2d objects
- [`plot(`*`<lyapunov_advisor>`*`)`](https://robustecologies.github.io/janos/reference/plot.lyapunov_advisor.md)
  : Plot method for a lyapunov_advisor object
- [`plot(`*`<lyapunov_function>`*`)`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
  : Plot method for a lyapunov_function object
- [`plot(`*`<lyapunov_report>`*`)`](https://robustecologies.github.io/janos/reference/plot.lyapunov_report.md)
  : Plot method for a lyapunov_report object
- [`plot(`*`<map_portrait>`*`)`](https://robustecologies.github.io/janos/reference/plot.map_portrait.md)
  : Plot method for map_portrait objects
- [`plot(`*`<mlmc_estimate>`*`)`](https://robustecologies.github.io/janos/reference/plot.mlmc_estimate.md)
  : Plot method for mlmc_estimate objects
- [`plot(`*`<observed_dyn_sim>`*`)`](https://robustecologies.github.io/janos/reference/plot.observed_dyn_sim.md)
  : Plot method for observed_dyn_sim objects
- [`plot(`*`<pdmp_portrait>`*`)`](https://robustecologies.github.io/janos/reference/plot.pdmp_portrait.md)
  : Plot method for pdmp_portrait objects
- [`plot(`*`<phase_portrait>`*`)`](https://robustecologies.github.io/janos/reference/plot.phase_portrait.md)
  : Plot method for phase_portrait objects
- [`plot(`*`<qsd_estimate>`*`)`](https://robustecologies.github.io/janos/reference/plot.qsd_estimate.md)
  : Plot method for qsd_estimate objects
- [`plot(`*`<rare_event_estimate>`*`)`](https://robustecologies.github.io/janos/reference/plot.rare_event_estimate.md)
  : Plot method for rare_event_estimate objects
- [`plot(`*`<sde_portrait>`*`)`](https://robustecologies.github.io/janos/reference/plot.sde_portrait.md)
  : Plot method for sde_portrait objects
- [`plot(`*`<sensitivity_result>`*`)`](https://robustecologies.github.io/janos/reference/plot.sensitivity_result.md)
  : Plot method for sensitivity_result objects
- [`plot(`*`<spectral_analysis>`*`)`](https://robustecologies.github.io/janos/reference/plot.spectral_analysis.md)
  : Plot method for spectral_analysis objects
- [`plot(`*`<zero_one_test>`*`)`](https://robustecologies.github.io/janos/reference/plot.zero_one_test.md)
  : Plot method for zero_one_test
