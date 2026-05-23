# Neural ODE inference via Bayesian gradient matching

Fits neural ordinary differential equations (NODEs) to multivariate time
series using the two-stage Bayesian Neural Gradient Matching (BNGM)
approach of Bonnaffe (2023). Stage 1 interpolates each variable with a
sinusoidal neural network and extracts analytical temporal derivatives
via the chain rule. Stage 2 fits a dual-pathway single-layer perceptron
(linear + exponential activations) to explain per-capita growth rates as
functions of all state variables, recovering the interaction structure
(Jacobian) and causal attribution (Geber contributions) of the dynamical
system without specifying its functional form.

## Usage

``` r
NODE(
  data,
  times = NULL,
  var.names = NULL,
  obs.neurons = "auto",
  obs.ensemble = 100L,
  obs.retain = 0.1,
  obs.sigma2 = NULL,
  proc.neurons = 5L,
  proc.ensemble = 100L,
  proc.retain = 0.1,
  proc.sd1 = NULL,
  proc.sd2 = NULL,
  proc.train.split = 2/3,
  cv.grid = NULL,
  cv.folds = 2L,
  bfgs.maxiter = 500L,
  threshold = 0.1,
  forecast = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  A numeric matrix or data.frame with variables in columns and time
  points in rows. Observations should be strictly positive because the
  method applies a log transform internally; values \<= 0 are clamped to
  0.005. Accepts `matrix`, `data.frame`, or `ts` objects.

- times:

  Optional numeric vector of observation times (length equal to
  `nrow(data)`). If `NULL`, sequential integers `1:nrow(data)` are used.
  The spacing `diff(times[1:2])` is used to de-standardize the temporal
  derivatives.

- var.names:

  Optional character vector of variable names. Defaults to
  `colnames(data)` or `V1, V2, ...`.

- obs.neurons:

  Integer, `"auto"`, or integer vector. Number of sinusoidal neurons per
  variable in the observation model. Each neuron contributes three
  learnable parameters (weight, frequency, phase), so the total
  parameter count per variable is `3 * obs.neurons`. If `"auto"`
  (default), a grid of candidate values is evaluated via a fast pilot
  process model fit, and the value yielding the highest minimum process
  R-squared across variables is selected. If an integer vector of length
  \> 1, it is used as the candidate grid for auto-selection. If a single
  integer, it is used directly without auto-selection. A reasonable
  manual range is `T/10` to `T/3` where `T` is the number of time
  points.

- obs.ensemble:

  Integer. Number of random initializations for the observation model's
  anchored ensemble (Pearce et al. 2018). Each initialization is
  optimized independently via BFGS. Larger values improve coverage of
  the loss surface at the cost of computation time (linear scaling).
  Default 100.

- obs.retain:

  Numeric in (0, 1\]. Fraction of best ensemble members to retain after
  ranking by log marginal likelihood. The retained members provide
  uncertainty quantification through their spread. Default 0.1.

- obs.sigma2:

  Numeric, `NULL`, or `"cv"`. Regularization coefficient \\c\\ for the
  observation model's Foresee-Hagan marginal posterior approximation. If
  `NULL` (default), uses the original formula \\c = 1/H\\ from Bonnaffe
  (2023), which is data-independent and well-calibrated across systems.
  If `"cv"`, selects via cross-validation (slower, can be noisy with
  small ensembles). If numeric, uses the specified value directly.
  Smaller values produce smoother interpolations; larger values allow
  the network to track higher-frequency fluctuations.

- proc.neurons:

  Integer. Number of neurons per pathway (linear and exponential) in the
  dual-pathway process model. The total parameter count per target
  variable is `2 * proc.neurons * (2 + n_vars)`. Default 5.

- proc.ensemble:

  Integer. Number of random initializations for the process model.
  Default 100.

- proc.retain:

  Numeric in (0, 1\]. Fraction of best process model ensemble members to
  retain, ranked by test-set prediction error. Default 0.1.

- proc.sd1:

  Numeric or `NULL`. Standard deviation \\\sigma_1\\ of the Gaussian
  likelihood in the process model posterior. Controls the tolerance for
  mismatch between observed and predicted per-capita growth rates. If
  `NULL` (default), set adaptively as
  `max(0.1, median(sd(growth_rates)) * 0.5)` after the observation model
  is fitted, ensuring the regularization scales with the data. For clean
  simulated data, 0.1 works well; for noisy data or systems with large
  growth rate variation, increase to 0.15-0.3.

- proc.sd2:

  Numeric or `NULL`. Standard deviation \\\sigma_2\\ of the spherical
  Gaussian prior on process model parameters. Controls weight
  regularization. If `NULL` (default), set equal to `proc.sd1`. Smaller
  values produce more constrained (sparser) interaction estimates.

- proc.train.split:

  Numeric in (0, 1). Fraction of time points used for training in the
  process model; the remainder is held out as a test set for ranking
  ensemble members. Default 2/3.

- cv.grid:

  Numeric vector or `NULL`. Grid of regularization coefficients for
  cross-validation (used only when `obs.sigma2 = "cv"`). Default: 12
  log-spaced values from 0.005 to 0.5.

- cv.folds:

  Integer. Number of cross-validation folds. Default 2.

- bfgs.maxiter:

  Integer. Maximum BFGS iterations per optimization run. The optimizer
  uses a self-contained C++ BFGS with backtracking Wolfe line search,
  avoiding R callback overhead. Default 500.

- threshold:

  Numeric. Geber contribution threshold for pruning weak interactions.
  Interactions whose normalized contribution \\\bar{c}\_{ij}\\ falls
  below this value are zeroed in the thresholded effects matrix and
  adjacency matrix. Default 0.1.

- forecast:

  Optional list for forward simulation via RK4 integration. Must contain
  `horizon` (numeric, time units to forecast beyond the last
  observation) and optionally `n.steps` (integer, number of output time
  steps; default 100). If `NULL`, no forecast is computed.

- verbose:

  Logical. Print progress messages with timing information.

## Value

An object of class `node_fit`, a list containing:

- fitted:

  Matrix (n_time x n_vars) of interpolated values on the original
  (untransformed) scale.

- derivatives:

  Matrix (n_time x n_vars) of temporal derivatives dY/dt.

- growth_rates:

  Matrix (n_time x n_vars) of per-capita growth rates \\r_i =
  (1/Y_i)(dY_i/dt)\\.

- proc.fitted:

  Matrix (n_time x n_vars) of process model predictions for the growth
  rates.

- obs.r2:

  Named numeric vector of observation model R-squared per variable.

- proc.r2:

  Named numeric vector of process model R-squared per variable.

- effects:

  Mean Jacobian matrix (n_vars x n_vars) after threshold pruning. Entry
  (i,j) is the mean effect of variable j on variable i.

- effects_raw:

  Mean Jacobian before thresholding.

- effects_sd:

  Ensemble standard deviation of the Jacobian (n_vars x n_vars).

- effects_ts:

  List of n_vars matrices, each (n_time x n_vars), giving the
  time-varying Jacobian rows.

- contributions:

  Geber contribution matrix (n_vars x n_vars), normalized so rows sum to
  1.

- adjacency:

  Binary adjacency matrix (n_vars x n_vars) after thresholding
  contributions.

- threshold:

  Numeric threshold used for pruning.

- obs_networks:

  Serialized observation model ensemble (for prediction).

- proc_networks:

  Serialized process model ensemble (for prediction and forecasting).

- cv:

  Cross-validation results (NULL if `obs.sigma2 != "cv"`).

- forecast:

  Forecast results (NULL if not requested). Contains `times`, `states`
  (matrix), and `x0`.

- normalization:

  Internal standardization coefficients for prediction/forecasting.

- call:

  The matched call.

- Minutes:

  Runtime in minutes.

## Details

The BNGM method decouples the inverse problem for dynamical systems into
two independent regression tasks, avoiding the costly numerical ODE
integration required by standard NODE methods (Chen et al. 2018). This
makes it 100-600 times faster than standard NODE fitting while achieving
higher accuracy in interaction recovery (Bonnaffe 2023, Table 1).

**Stage 1 (observation model):** Each variable is interpolated with a
sinusoidal neural network \\\hat{y}\_i(t) = \sum\_{j=1}^{H} w_j
\sin(\pi(t \cdot a_j + \phi_j))\\, where \\H\\ is the number of neurons
and the parameters \\(w_j, a_j, \phi_j)\\ are weights, learnable
frequencies, and phases. Unlike a truncated Fourier series, the
frequencies are optimized to concentrate representational capacity on
the frequencies actually present in the signal. The temporal derivative
\\d\hat{y}\_i/dt = \sum_j w_j \pi a_j \cos(\pi(t a_j + \phi_j))\\ is
available analytically, avoiding finite-difference noise amplification.
Regularization uses the Foresee and Hagan (1997) marginal posterior
approximation.

**Stage 2 (process model):** The per-capita growth rate \\r_i =
(1/Y_i)(dY_i/dt)\\ for each variable is modeled as a function of all
state variables using a dual-pathway single-layer perceptron:
\\f\_{p,i}(\mathbf{x}) = \mathbf{v}\_{\mathrm{lin}}^\top
(\mathbf{b}\_{\mathrm{lin}} + \mathbf{W}\_{\mathrm{lin}} \mathbf{x}) +
\mathbf{v}\_{\exp}^\top \exp(\mathbf{b}\_{\exp} + \mathbf{W}\_{\exp}
\mathbf{x})\\. The linear pathway captures density-independent dynamics;
the exponential pathway captures nonlinear density-dependent
interactions (predation saturation, competitive exclusion, Allee
effects).

**Uncertainty quantification:** Both stages use anchored Bayesian
ensembles (Pearce et al. 2018): \\K\\ random parameter initializations
are optimized independently via BFGS (parallelized with OpenMP), and the
best \\\rho K\\ fraction is retained. Ensemble spread provides
approximate posterior uncertainty without MCMC sampling.

**Interaction inference:** The Jacobian matrix \\e\_{ij} = \partial
f\_{p,i} / \partial x_j\\ gives the effect of variable \\j\\ on the
growth rate of variable \\i\\ (positive = facilitation, negative =
suppression). The Geber method (Hairston et al. 2005) weights these
effects by the temporal variation they explain, yielding normalized
contributions \\\bar{c}\_{ij}\\ that sum to 1 per target variable.

**Limitations:** The method assumes all dynamically important variables
are observed. Systems that converge monotonically to equilibrium without
oscillations provide insufficient temporal variation for reliable
interaction recovery. Numerical differentiation amplifies observation
noise, so noisy data degrades the process model more than the
observation model.

## References

Bonnaffe, W. (2023). Fast fitting of neural ordinary differential
equations by Bayesian neural gradient matching. *Methods in Ecology and
Evolution*, 14(6), 1456–1468.
[doi:10.1111/2041-210X.14121](https://doi.org/10.1111/2041-210X.14121)

Chen, R.T.Q., Rubanova, Y., Bettencourt, J. & Duvenaud, D.K. (2018).
Neural ordinary differential equations. *Advances in Neural Information
Processing Systems*, 31. <https://arxiv.org/abs/1806.07366>

Pearce, T., Leibfried, F. & Brintrup, A. (2020). Uncertainty in neural
networks: approximately Bayesian ensembling. *Proceedings of the
International Conference on Artificial Intelligence and Statistics
(AISTATS)*, PMLR 108, 234–244.

Foresee, F.D. & Hagan, M.T. (1997). Gauss-Newton approximation to
Bayesian learning. *Proceedings of the International Joint Conference on
Neural Networks (IJCNN)*, 1930–1935.
[doi:10.1109/ICNN.1997.614194](https://doi.org/10.1109/ICNN.1997.614194)

Hairston, N.G., Ellner, S.P., Geber, M.A., Yoshida, T. & Fox, J.A.
(2005). Rapid evolution and the convergence of ecological and
evolutionary time. *Ecology Letters*, 8(10), 1114–1127.
[doi:10.1111/j.1461-0248.2005.00812.x](https://doi.org/10.1111/j.1461-0248.2005.00812.x)

## See also

[`NODE_predict`](https://robustecologies.github.io/lucifer/reference/NODE_predict.md)
for forecasting,
[`SDE`](https://robustecologies.github.io/lucifer/reference/SDE.md) for
parametric stochastic differential equation models,
[`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md) for
simulation-based inference.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate a Lotka-Volterra predator-prey system
set.seed(42)
n <- 100; dt <- 0.1
prey <- predator <- numeric(n)
prey[1] <- 1.0; predator[1] <- 0.5
for (i in 2:n) {
  prey[i] <- prey[i-1] + dt * prey[i-1] * (1 - 0.5 * predator[i-1])
  predator[i] <- predator[i-1] + dt * predator[i-1] * (-0.5 + 0.3 * prey[i-1])
}
data <- cbind(prey = prey, predator = predator)
times <- seq(dt, n * dt, by = dt)

# Fit with automatic obs.neurons selection (default)
fit <- NODE(data, times = times, obs.ensemble = 50, proc.ensemble = 50)

# Inspect recovered interaction structure
print(fit)
summary(fit)

# Auto-selection metadata (which H was chosen and why)
fit$obs.neurons.selection

# Visualize all plot types
plot(fit)                          # interpolation (default)
plot(fit, type = "dynamics")       # per-capita growth rates
plot(fit, type = "effects")        # Jacobian heatmap
plot(fit, type = "contributions")  # Geber causal attribution
plot(fit, type = "network")        # circular interaction network
plot(fit, type = "phase")          # phase space trajectories

# Compare against ground truth Jacobian
# True: prey->pred = +0.3, pred->prey = -0.5
print(round(fit$effects_raw, 3))

# Forecast via RK4 and visualize
pred <- NODE_predict(fit, horizon = 5, n.steps = 60)
print(pred)
plot(pred)                         # trajectory with forecast zone
plot(pred, type = "phase")         # phase space forecast
} # }
```
