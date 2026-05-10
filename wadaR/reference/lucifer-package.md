# lucifer: complete environment for Bayesian inference

Provides a complete environment for Bayesian inference using a variety
of different samplers. The package implements 43 MCMC algorithms,
approximately 50 probability distributions, numerous matrix utilities,
and comprehensive diagnostic tools. Computationally intensive operations
are implemented in C++ via Rcpp and RcppArmadillo with OpenMP
parallelization for multicore acceleration.

## Details

The main function is
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
which performs MCMC sampling with any of 43 algorithms. Multi-chain MCMC
is supported directly via the `Chains`, `CPUs`, and `Type` arguments.
The recommended setting for production inference is `Chains = 4`.

Additional inference methods include
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md)
(GP-based probabilistic integration),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)
(population Monte Carlo), and
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

The package provides extensive MCMC diagnostics
([`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
[`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
among others), model comparison tools
([`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`WAIC`](https://robustecologies.github.io/lucifer/reference/WAIC.md),
[`LOO`](https://robustecologies.github.io/lucifer/reference/LOO.md)),
posterior predictive checks
([`PosteriorChecks`](https://robustecologies.github.io/lucifer/reference/PosteriorChecks.md)),
and visualization methods for all fitted objects.

## References

The vignettes bundled with this package provide extensive background:
the lucifer Tutorial, the Examples compendium, and the Bayesian
Inference introduction.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`BayesianQuadrature`](https://robustecologies.github.io/lucifer/reference/BayesianQuadrature.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md),
[`PMC`](https://robustecologies.github.io/lucifer/reference/PMC.md)

## Author

Pablo Almaraz (creator, maintainer), Byron Hall (creator), Martina Hall
(creator)
