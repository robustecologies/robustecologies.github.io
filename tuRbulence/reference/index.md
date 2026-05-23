# Package index

## Package overview

Top-level documentation entry.

- [`tuRbulence`](https://robustecologies.github.io/tuRbulence/reference/tuRbulence-package.md)
  [`tuRbulence-package`](https://robustecologies.github.io/tuRbulence/reference/tuRbulence-package.md)
  : tuRbulence: Stochastic chaos and turbulent attractors in dynamical
  systems

## Unified simulation interface

Single dispatcher and observable extractor.

- [`turbulence_sim()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_sim.md)
  : Simulate a turbulent dynamical system
- [`simulate_system()`](https://robustecologies.github.io/tuRbulence/reference/simulate_system.md)
  : Simulate a dynamical system (legacy interface)
- [`get_primary_series()`](https://robustecologies.github.io/tuRbulence/reference/get_primary_series.md)
  : Extract primary time series from simulation

## Canonical chaotic systems

Lorenz-63, Lorenz-84 and Rossler integrators with optional
Ornstein-Uhlenbeck forcing.

- [`lorenz_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz_sim.md)
  : Simulate the Lorenz system
- [`lorenz84_sim()`](https://robustecologies.github.io/tuRbulence/reference/lorenz84_sim.md)
  : Simulate the Lorenz-84 atmospheric model
- [`rossler_sim()`](https://robustecologies.github.io/tuRbulence/reference/rossler_sim.md)
  : Simulate the Rössler system

## Geophysical and climate systems

Charney-DeVore (3- and 6-mode), Stommel two-box thermohaline, Hasselmann
stochastic climate, von-Karman swirling-flow Duffing.

- [`charney_devore_sim()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_sim.md)
  : Simulate Charney-DeVore atmospheric blocking model
- [`charney_devore_6mode()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_6mode.md)
  : Simulate the full 6-mode Charney-DeVore model with bistability
- [`charney_devore_batch()`](https://robustecologies.github.io/tuRbulence/reference/charney_devore_batch.md)
  : Batch simulation over forcing values
- [`stommel_sim()`](https://robustecologies.github.io/tuRbulence/reference/stommel_sim.md)
  : Simulate Stommel box model for thermohaline circulation
- [`stommel_batch()`](https://robustecologies.github.io/tuRbulence/reference/stommel_batch.md)
  : Batch simulation over freshwater flux values
- [`hasselmann_sim()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_sim.md)
  : Simulate Hasselmann stochastic climate model
- [`hasselmann_batch()`](https://robustecologies.github.io/tuRbulence/reference/hasselmann_batch.md)
  : Batch simulation over forcing values
- [`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md)
  : Simulate stochastic Duffing oscillator for von Kármán flow
- [`vonkarman_batch()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_batch.md)
  : Batch simulation over multiple control parameter values
- [`vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_peaks.md)
  : Extract peaks from time series for attractor embedding
- [`vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_attractor.md)
  : Create embedded attractor from peaks

## Phase-space reconstruction

Takens delay embedding and Cao false-nearest-neighbours dimension.

- [`turbulence_embed()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_embed.md)
  : Create delay embedding from time series
- [`turbulence_cao()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_cao.md)
  : Estimate embedding dimension using Cao's method

## Lyapunov spectra

Wolf and Rosenstein algorithms for the largest exponent.

- [`turbulence_lyapunov()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_lyapunov.md)
  : Estimate largest Lyapunov exponent

## Bifurcation analysis

Parameter sweeps with regime classification.

- [`turbulence_bifurcation()`](https://robustecologies.github.io/tuRbulence/reference/turbulence_bifurcation.md)
  : Compute Lyapunov exponents across parameter range
- [`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)
  : Create bifurcation panel showing attractor evolution

## 3D visualization

plotly-based attractors, trajectories, animations.

- [`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md)
  : Create static 3D embedded attractor
- [`create_trajectory_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_trajectory_3d.md)
  : Create static 3D phase space trajectory
- [`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md)
  : Create animated 3D attractor with playback controls
- [`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md)
  : Create animated attractor showing progressive formation
- [`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md)
  : Save plotly attractor as interactive HTML
- [`export_animation_frames()`](https://robustecologies.github.io/tuRbulence/reference/export_animation_frames.md)
  : Export animation frames as images for video creation

## Interactive dashboard

- [`shiny_tuRbulence()`](https://robustecologies.github.io/tuRbulence/reference/shiny_tuRbulence.md)
  : Interactive Shiny application for tuRbulence

## S3 methods

Print, summary and plot methods for the 18 exported S3 classes.

- [`print(`*`<charney_devore_6mode>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_6mode.md)
  :

  Print method for `charney_devore_6mode` objects

- [`print(`*`<charney_devore_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_batch.md)
  : Print method for Charney-DeVore batch results

- [`print(`*`<charney_devore_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.charney_devore_sim.md)
  :

  Print method for `charney_devore_sim` objects

- [`print(`*`<hasselmann_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.hasselmann_batch.md)
  : Print method for Hasselmann batch results

- [`print(`*`<hasselmann_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.hasselmann_sim.md)
  :

  Print method for `hasselmann_sim` objects

- [`print(`*`<lorenz84_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz84_sim.md)
  :

  Print method for `lorenz84_sim` objects

- [`print(`*`<lorenz_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.lorenz_sim.md)
  :

  Print method for `lorenz_sim` objects

- [`print(`*`<rossler_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.rossler_sim.md)
  :

  Print method for `rossler_sim` objects

- [`print(`*`<stommel_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.stommel_batch.md)
  : Print method for Stommel batch results

- [`print(`*`<stommel_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.stommel_sim.md)
  :

  Print method for `stommel_sim` objects

- [`print(`*`<turbulence_bifurcation>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_bifurcation.md)
  :

  Print method for `turbulence_bifurcation` objects

- [`print(`*`<turbulence_cao>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_cao.md)
  :

  Print method for `turbulence_cao` objects

- [`print(`*`<turbulence_embedding>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_embedding.md)
  :

  Print method for `turbulence_embedding` objects

- [`print(`*`<turbulence_lyapunov>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.turbulence_lyapunov.md)
  : Print method for Lyapunov-exponent estimation results

- [`print(`*`<vonkarman_attractor>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_attractor.md)
  :

  Print method for `vonkarman_attractor` objects

- [`print(`*`<vonkarman_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_batch.md)
  : Print method for von Karman batch results

- [`print(`*`<vonkarman_peaks>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_peaks.md)
  : Print method for von Karman peak extraction

- [`print(`*`<vonkarman_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_sim.md)
  : Print method for von Karman stochastic Duffing simulations

- [`summary(`*`<charney_devore_6mode>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_6mode.md)
  :

  Summary method for `charney_devore_6mode` objects

- [`summary(`*`<charney_devore_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_batch.md)
  : Summary method for Charney-DeVore batch results

- [`summary(`*`<charney_devore_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.charney_devore_sim.md)
  :

  Summary method for `charney_devore_sim` objects

- [`summary(`*`<hasselmann_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.hasselmann_batch.md)
  : Summary method for Hasselmann batch results

- [`summary(`*`<hasselmann_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.hasselmann_sim.md)
  :

  Summary method for `hasselmann_sim` objects

- [`summary(`*`<lorenz84_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz84_sim.md)
  :

  Summary method for `lorenz84_sim` objects

- [`summary(`*`<lorenz_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.lorenz_sim.md)
  :

  Summary method for `lorenz_sim` objects

- [`summary(`*`<rossler_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.rossler_sim.md)
  :

  Summary method for `rossler_sim` objects

- [`summary(`*`<stommel_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.stommel_batch.md)
  : Summary method for Stommel batch results

- [`summary(`*`<stommel_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.stommel_sim.md)
  :

  Summary method for `stommel_sim` objects

- [`summary(`*`<turbulence_bifurcation>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_bifurcation.md)
  :

  Summary method for `turbulence_bifurcation` objects

- [`summary(`*`<turbulence_cao>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_cao.md)
  : Summary method for Cao's method results

- [`summary(`*`<turbulence_embedding>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_embedding.md)
  : Summary method for turbulence embedding objects

- [`summary(`*`<turbulence_lyapunov>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.turbulence_lyapunov.md)
  : Summary method for Lyapunov-exponent estimation results

- [`summary(`*`<vonkarman_attractor>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_attractor.md)
  : Summary method for von Karman attractor objects

- [`summary(`*`<vonkarman_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_batch.md)
  : Summary method for von Karman batch results

- [`summary(`*`<vonkarman_peaks>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_peaks.md)
  : Summary method for von Karman peak extraction

- [`summary(`*`<vonkarman_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_sim.md)
  :

  Summary method for `vonkarman_sim` objects

- [`plot(`*`<charney_devore_6mode>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_6mode.md)
  :

  Plot method for `charney_devore_6mode` objects

- [`plot(`*`<charney_devore_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_batch.md)
  : Plot Charney-DeVore batch simulation results

- [`plot(`*`<charney_devore_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.charney_devore_sim.md)
  : Plot a Charney-DeVore atmospheric blocking simulation

- [`plot(`*`<hasselmann_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.hasselmann_batch.md)
  : Plot Hasselmann batch simulation results

- [`plot(`*`<hasselmann_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.hasselmann_sim.md)
  : Plot a Hasselmann stochastic climate simulation

- [`plot(`*`<lorenz84_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.lorenz84_sim.md)
  : Plot a Lorenz-84 atmospheric model simulation

- [`plot(`*`<lorenz_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.lorenz_sim.md)
  : Plot a Lorenz system simulation

- [`plot(`*`<rossler_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.rossler_sim.md)
  : Plot a Rössler system simulation

- [`plot(`*`<stommel_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.stommel_batch.md)
  : Plot Stommel batch simulation results

- [`plot(`*`<stommel_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.stommel_sim.md)
  : Plot a Stommel box model simulation

- [`plot(`*`<turbulence_bifurcation>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_bifurcation.md)
  : Plot bifurcation diagram

- [`plot(`*`<turbulence_cao>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_cao.md)
  : Plot Cao's method results for embedding dimension

- [`plot(`*`<turbulence_embedding>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_embedding.md)
  : Plot method for turbulence embedding objects

- [`plot(`*`<turbulence_lyapunov>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.turbulence_lyapunov.md)
  : Plot Lyapunov exponent estimation results

- [`plot(`*`<vonkarman_attractor>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_attractor.md)
  : Plot a von Kármán embedded attractor

- [`plot(`*`<vonkarman_batch>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_batch.md)
  : Plot von Kármán batch simulation results

- [`plot(`*`<vonkarman_peaks>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_peaks.md)
  : Plot method for von Karman peak extraction

- [`plot(`*`<vonkarman_sim>`*`)`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_sim.md)
  : Plot a von Kármán turbulent flow simulation

- [`plot_bifurcation_panel()`](https://robustecologies.github.io/tuRbulence/reference/plot_bifurcation_panel.md)
  : Create bifurcation panel showing attractor evolution
