# symplectoR: Symplectic and Accelerated-Gradient Optimization via Dissipative Hamiltonian Dynamics

Accelerated gradient methods implemented as structure-preserving
discretizations of continuous dissipative Hamiltonian systems. Provides
the relativistic gradient descent master kernel of Franca et al (2020)
<https://arxiv.org/abs/1903.04100>, which subsumes classical momentum
and Nesterov acceleration as exact parameter settings; dissipative
presymplectic leapfrog integrators with general damping schedules
following Franca et al. (2021)
[doi:10.1088/1742-5468/abf5d4](https://doi.org/10.1088/1742-5468/abf5d4)
; the fully explicit symplectic leapfrog compositions of the polynomial
and exponential Bregman subfamilies with gradient momentum restarting
and temporal looping following Duruisseaux and Leok (2023)
<https://arxiv.org/abs/2207.11460>; the rate-matching three-sequence
discretization of Wibisono et al (2016)
<https://doi.org/10.1073/pnas.1614734113>; and a classical split-step
Fourier simulation of quantum Hamiltonian descent for non-smooth global
optimization following Leng et al (2025)
<https://arxiv.org/abs/2503.15878>. An inverse-modeling front end
estimates model parameters from data by trajectory matching or
likelihood minimization. Performance-critical kernels are implemented in
C++ via 'Rcpp' and 'RcppArmadillo' with OpenMP parallelization over
multi-start ensembles, parameter sweeps and pseudospectral grids.

## See also

Useful links:

- <https://github.com/robustecologies/symplectoR>

- <https://robustecologies.github.io/symplectoR>

- Report bugs at <https://github.com/robustecologies/symplectoR/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))

Authors:

- Pablo Almaraz <pablo.almaraz@csic.es>
  ([ORCID](https://orcid.org/0000-0003-1416-2695))
