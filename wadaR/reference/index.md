# Package index

## Package overview

Top-level documentation entry.

- [`wadaR`](https://robustecologies.github.io/wadaR/reference/wadaR-package.md)
  [`wadaR-package`](https://robustecologies.github.io/wadaR/reference/wadaR-package.md)
  : wadaR: Wada basins of attraction and bifurcations in chaotic
  dynamical systems

## Attractor specification

Constructors describing what the integrator should recognise as an
asymptotic outcome.

- [`attractor_point()`](https://robustecologies.github.io/wadaR/reference/attractor_point.md)
  : Create a point attractor specification
- [`attractor_cycle()`](https://robustecologies.github.io/wadaR/reference/attractor_cycle.md)
  : Create a limit cycle attractor specification
- [`attractor_exit()`](https://robustecologies.github.io/wadaR/reference/attractor_exit.md)
  : Create an escape channel specification
- [`attractor_outcome()`](https://robustecologies.github.io/wadaR/reference/attractor_outcome.md)
  : Create a discrete outcome specification

## System catalogue

Canonical Wada systems and the user-extensible compiled-system
interface.

- [`forced_damped_pendulum()`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
  : Forced damped pendulum system
- [`henon_heiles_system()`](https://robustecologies.github.io/wadaR/reference/henon_heiles_system.md)
  : Henon-Heiles Hamiltonian system
- [`newton_fractal_system()`](https://robustecologies.github.io/wadaR/reference/newton_fractal_system.md)
  : Newton fractal system
- [`multispecies_competition()`](https://robustecologies.github.io/wadaR/reference/multispecies_competition.md)
  : Huisman-Weissing multispecies competition system
- [`compiled_system()`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)
  : Create a compiled dynamical system with parallel basin computation
- [`is.compiled_system()`](https://robustecologies.github.io/wadaR/reference/is.compiled_system.md)
  : Check if an object is a compiled_system
- [`compile_basin_function()`](https://robustecologies.github.io/wadaR/reference/compile_basin_function.md)
  : Compile a parallel basin computation function for user-defined
  dynamics
- [`initial_state_fixed_energy()`](https://robustecologies.github.io/wadaR/reference/initial_state_fixed_energy.md)
  : Create initial state function for fixed energy Hamiltonian systems

## Basin computation

Trajectory integration over a grid of initial conditions; basin
labelling.

- [`compute_basins()`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
  : Compute basins of attraction for a 2D dynamical system
- [`compute_basins_3d()`](https://robustecologies.github.io/wadaR/reference/compute_basins_3d.md)
  : Compute 3D basins of attraction
- [`compute_newton_basins()`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md)
  : Compute Newton fractal basins
- [`compute_competition_basins()`](https://robustecologies.github.io/wadaR/reference/compute_competition_basins.md)
  : Compute basins of attraction for multispecies competition
- [`simulate_competition()`](https://robustecologies.github.io/wadaR/reference/simulate_competition.md)
  : Simulate Huisman-Weissing competition dynamics
- [`detect_attractors()`](https://robustecologies.github.io/wadaR/reference/detect_attractors.md)
  : Detect attractors automatically from trajectory analysis
- [`extract_basin()`](https://robustecologies.github.io/wadaR/reference/extract_basin.md)
  : Extract basin at specific parameter value
- [`get_boundary()`](https://robustecologies.github.io/wadaR/reference/get_boundary.md)
  : Get basin boundary points
- [`merge_basins()`](https://robustecologies.github.io/wadaR/reference/merge_basins.md)
  : Merge basins for Wada analysis

## Wada detection

Three published methods plus the unified dispatcher. All require \>= 3
attractors.

- [`detect_wada()`](https://robustecologies.github.io/wadaR/reference/detect_wada.md)
  : Detect Wada basins using multiple methods
- [`wada_grid_method()`](https://robustecologies.github.io/wadaR/reference/wada_grid_method.md)
  : Detect Wada basins using the grid method
- [`wada_merging_method()`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
  : Detect Wada basins using the merging method
- [`wada_straddle_method()`](https://robustecologies.github.io/wadaR/reference/wada_straddle_method.md)
  : Detect Wada basins using the saddle-straddle method
- [`hausdorff_distance()`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md)
  : Compute Hausdorff distance between two point sets

## Basin entropy

Information-theoretic boundary diagnostic (Daza 2016).

- [`basin_entropy()`](https://robustecologies.github.io/wadaR/reference/basin_entropy.md)
  : Compute basin entropy

## Bifurcation analysis

Parameter sweeps with regime classification and animation.

- [`bifurcation_basins()`](https://robustecologies.github.io/wadaR/reference/bifurcation_basins.md)
  : Bifurcation analysis of basins of attraction
- [`animate()`](https://robustecologies.github.io/wadaR/reference/animate.md)
  : Animate an object
- [`animate(`*`<bifurcation_result>`*`)`](https://robustecologies.github.io/wadaR/reference/animate.bifurcation_result.md)
  : Animate bifurcation results

## 3D visualization

plotly-based volumetric, slice and rotating-camera visualizations.

- [`plot_3d_basins()`](https://robustecologies.github.io/wadaR/reference/plot_3d_basins.md)
  : Plot 3D basins using plotly
- [`slice_3d_basins()`](https://robustecologies.github.io/wadaR/reference/slice_3d_basins.md)
  : Extract 2D slice from 3D basins
- [`animate_3d_rotation()`](https://robustecologies.github.io/wadaR/reference/animate_3d_rotation.md)
  : Create animated rotation of 3D basins

## Visualization helpers

Discrete colour palettes used across the plot methods.

- [`palettes()`](https://robustecologies.github.io/wadaR/reference/palettes.md)
  : Color palettes (delegated to koloRo)
- [`palette_ramp()`](https://robustecologies.github.io/wadaR/reference/palette_ramp.md)
  : Interpolate a palette to N colors (delegated to koloRo)

## Interactive dashboard

- [`shinywadaR()`](https://robustecologies.github.io/wadaR/reference/shinywadaR.md)
  : Launch interactive wadaR Shiny application

## S3 methods

Print, summary and plot methods for the 12 exported S3 classes.

- [`print(`*`<attractor_detection>`*`)`](https://robustecologies.github.io/wadaR/reference/print.attractor_detection.md)
  : Print method for attractor_detection objects
- [`print(`*`<basin_entropy_result>`*`)`](https://robustecologies.github.io/wadaR/reference/print.basin_entropy_result.md)
  : Print method for basin_entropy_result
- [`print(`*`<basin_result_3d>`*`)`](https://robustecologies.github.io/wadaR/reference/print.basin_result_3d.md)
  : Print method for basin_result_3d
- [`print(`*`<bifurcation_result>`*`)`](https://robustecologies.github.io/wadaR/reference/print.bifurcation_result.md)
  : Print method for bifurcation_result
- [`print(`*`<compiled_system>`*`)`](https://robustecologies.github.io/wadaR/reference/print.compiled_system.md)
  : Print method for compiled_system objects
- [`print(`*`<merged_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/print.merged_basins.md)
  : Print method for merged_basins
- [`print(`*`<multispecies_competition>`*`)`](https://robustecologies.github.io/wadaR/reference/print.multispecies_competition.md)
  : Print method for multispecies competition system
- [`print(`*`<wada_analysis>`*`)`](https://robustecologies.github.io/wadaR/reference/print.wada_analysis.md)
  : Print method for wada_analysis
- [`print(`*`<wada_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/print.wada_basins.md)
  : Print method for wada_basins
- [`print(`*`<wada_grid_result>`*`)`](https://robustecologies.github.io/wadaR/reference/print.wada_grid_result.md)
  : Print method for wada_grid_result
- [`print(`*`<wada_merging_result>`*`)`](https://robustecologies.github.io/wadaR/reference/print.wada_merging_result.md)
  : Print method for wada_merging_result
- [`print(`*`<wada_straddle_result>`*`)`](https://robustecologies.github.io/wadaR/reference/print.wada_straddle_result.md)
  : Print method for wada_straddle_result
- [`summary(`*`<attractor_detection>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.attractor_detection.md)
  : Summary method for attractor_detection objects
- [`summary(`*`<basin_entropy_result>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.basin_entropy_result.md)
  : Summary method for basin_entropy_result objects
- [`summary(`*`<basin_result_3d>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.basin_result_3d.md)
  : Summary method for basin_result_3d objects
- [`summary(`*`<bifurcation_result>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.bifurcation_result.md)
  : Summary method for bifurcation_result objects
- [`summary(`*`<compiled_system>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.compiled_system.md)
  : Summary method for compiled_system objects
- [`summary(`*`<merged_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.merged_basins.md)
  : Summary method for merged_basins objects
- [`summary(`*`<multispecies_competition>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.multispecies_competition.md)
  : Summary method for multispecies_competition objects
- [`summary(`*`<wada_analysis>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.wada_analysis.md)
  : Summary method for wada_analysis objects
- [`summary(`*`<wada_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.wada_basins.md)
  : Summary method for wada_basins objects
- [`summary(`*`<wada_grid_result>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.wada_grid_result.md)
  : Summary method for wada_grid_result objects
- [`summary(`*`<wada_merging_result>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.wada_merging_result.md)
  : Summary method for wada_merging_result objects
- [`summary(`*`<wada_straddle_result>`*`)`](https://robustecologies.github.io/wadaR/reference/summary.wada_straddle_result.md)
  : Summary method for wada_straddle_result objects
- [`plot(`*`<attractor_detection>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.attractor_detection.md)
  : Plot method for attractor_detection objects
- [`plot(`*`<basin_entropy_result>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.basin_entropy_result.md)
  : Plot method for basin_entropy_result
- [`plot(`*`<basin_result_3d>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.basin_result_3d.md)
  : Plot method for 3D basin results
- [`plot(`*`<bifurcation_result>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.bifurcation_result.md)
  : Plot method for bifurcation_result
- [`plot(`*`<compiled_system>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.compiled_system.md)
  : Plot method for compiled system objects
- [`plot(`*`<merged_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.merged_basins.md)
  : Plot method for merged_basins
- [`plot(`*`<multispecies_competition>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.multispecies_competition.md)
  : Plot method for multispecies competition objects
- [`plot(`*`<wada_analysis>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.wada_analysis.md)
  : Plot method for Wada analysis results
- [`plot(`*`<wada_basins>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.wada_basins.md)
  : Plot method for wada_basins (S3)
- [`plot(`*`<wada_grid_result>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.wada_grid_result.md)
  : Plot grid method results
- [`plot(`*`<wada_merging_result>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.wada_merging_result.md)
  : Plot merging method results
- [`plot(`*`<wada_straddle_result>`*`)`](https://robustecologies.github.io/wadaR/reference/plot.wada_straddle_result.md)
  : Plot saddle-straddle method results
- [`plot_3d_basins()`](https://robustecologies.github.io/wadaR/reference/plot_3d_basins.md)
  : Plot 3D basins using plotly
