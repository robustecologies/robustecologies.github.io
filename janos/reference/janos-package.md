# janos: General-Purpose Simulation Framework for Spatio-Temporal Stochastic Dynamical Systems

A unified framework for specifying, simulating and analysing dynamical
systems from a single formula-based interface. Covers ordinary, delay
and stochastic differential equations; 1D and 2D partial differential
equations via the method of lines; discrete maps; jump-diffusion
processes; piecewise deterministic Markov processes; continuous-time
Markov chains; and the reaction-diffusion master equation on regular
grids and arbitrary graphs. Stochastic models support correlated, Levy
alpha-stable, fractional Brownian and colored (1/f^beta) noise. Model
formulas are compiled to C++ on demand and dispatched to a suite of
solvers (Dormand-Prince RK4(5), Rosenbrock Rodas3 with symbolic Jacobian
generation, method of lines, Gillespie Direct and Next-Reaction,
adaptive tau-leaping, hybrid SSA/CLE, Euler-Maruyama, Milstein, and
Lewis-Shedler thinning). The analysis layer provides phase and
qualitative portraits for ODE, SDE, DDE and PDMP systems; numerical
bifurcation continuation with fold and Hopf detection; Lyapunov function
construction via eight algebraic and numerical methods with a
family-aware advisor dispatching over any model specification;
Fokker-Planck stationary densities, quasi-potential landscapes and
Kramers escape rates; continuous adjoint sensitivity; quasi-stationary
distribution and rare-event probability estimation; multi-level Monte
Carlo variance reduction; and ensemble simulation with two-tier C++
OpenMP and R-level parallelism.

## See also

Useful links:

- <https://github.com/robustecologies/janos>

- <https://robustecologies.github.io/janos>

- Report bugs at <https://github.com/robustecologies/janos/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))
