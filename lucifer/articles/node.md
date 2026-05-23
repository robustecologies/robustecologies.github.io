# Neural ODE inference via Bayesian gradient matching

## Introduction

A central problem in ecology and dynamical systems science is inferring
interaction structure from observed time series: given measurements of
multiple interacting variables over time, which variables drive which,
and how strongly? Parametric approaches, whether Lotka-Volterra models,
state-space formulations, or mechanistic ODEs, require specifying the
functional form of interactions in advance. When the true functional
form is unknown or misspecified, parameter estimates absorb the model
error and can lead to qualitatively wrong conclusions about interaction
structure.

Neural ordinary differential equations (NODEs; Chen et al. 2018
[\[3\]](#ref3)) offer a model-free alternative by using neural networks
to represent the unknown right-hand side of the ODE system. However,
standard NODE fitting requires numerically integrating the entire system
at every parameter update, making the computational cost proportional to
the number of time steps multiplied by the number of solver sub-steps.
For ecological datasets with complex, potentially stiff dynamics, this
cost is prohibitive. The adjoint method reduces memory but not
computation, and backpropagation through the solver introduces
additional numerical instabilities.

Bayesian Neural Gradient Matching (BNGM; Bonnaffé 2023 [\[1\]](#ref1))
solves this problem by decoupling ODE fitting into two independent
regression stages, avoiding numerical integration entirely. Stage 1
interpolates each time series with a sinusoidal neural network whose
temporal derivatives are available analytically via the chain rule.
Stage 2 fits a dual-pathway neural network to predict per-capita growth
rates from the interpolated state variables. Because both stages reduce
to standard regression problems with known gradients, the method is
100-600 times faster than standard NODE fitting while recovering
interaction structure with significantly higher accuracy (Bonnaffé 2023
reports ~5-10% error vs ~50% for standard NODEs on benchmark systems).

The
[`NODE()`](https://robustecologies.github.io/lucifer/reference/NODE.md)
function in lucifer implements BNGM with a C++/Armadillo backend
parallelized via OpenMP, providing further speedups of 10-50x over the
reference R implementation. This vignette develops the mathematical
foundations, demonstrates the method through four examples of increasing
complexity, and provides practical guidelines for parameter selection.

## Mathematical foundations

### The inverse problem for dynamical systems

Consider a system of \\N\\ interacting variables \\\mathbf{y}(t) =
(y_1(t), \ldots, y_N(t))\\ governed by an unknown ODE system:

\\\frac{dy_i}{dt} = g_i(\mathbf{y}, t), \qquad i = 1, \ldots, N\\

where the functions \\g_i\\ encode the ecological interactions. The
inverse problem is to recover the \\g_i\\ (or at least their qualitative
structure, i.e., which partial derivatives \\\partial g_i / \partial
y_j\\ are non-zero and what sign they carry) from discrete, noisy
observations \\\\(t_n, Y\_{n,1}, \ldots, Y\_{n,N})\\\_{n=1}^{T}\\.

Standard approaches parameterize \\g_i\\ as a known functional form
(e.g., \\g_i = r_i y_i (1 - \sum_j \alpha\_{ij} y_j / K_i)\\ for
competitive Lotka-Volterra) and estimate the parameters by maximum
likelihood or Bayesian inference. BNGM instead approximates each \\g_i\\
with a neural network \\f\_{\theta,i}\\, requiring no prior knowledge of
the functional form.

### Stage 1: sinusoidal observation model

The first stage interpolates the observed time series to obtain smooth
estimates \\\hat{y}\_i(t)\\ and their temporal derivatives
\\\dot{\hat{y}}\_i(t)\\. Bonnaffé (2023) uses a sum-of-sinusoids
architecture rather than the more common sigmoid or ReLU activations.
The sinusoidal basis has two advantages: it naturally captures the
oscillatory dynamics common in ecological systems, and its derivatives
are analytically available without numerical differentiation.

For each variable \\i\\, the interpolation network has \\H\\ neurons:

\\\hat{y}\_i(t, \omega) = \sum\_{j=1}^{H} w_j \sin\bigl(\pi(t \cdot
a_j + \phi_j)\bigr)\\

where \\\omega = (w_1, \ldots, w_H, a_1, \ldots, a_H, \phi_1, \ldots,
\phi_H)\\ contains \\3H\\ parameters: output weights \\w_j\\,
frequencies \\a_j\\, and phases \\\phi_j\\. The temporal derivative
follows directly:

\\\frac{\partial \hat{y}\_i}{\partial t}(t, \omega) = \sum\_{j=1}^{H}
w_j \cdot \pi a_j \cdot \cos\bigl(\pi(t \cdot a_j + \phi_j)\bigr)\\

This is a Fourier-like basis expansion, but unlike a truncated Fourier
series, the frequencies \\a_j\\ are learnable parameters rather than
fixed harmonics. This gives the network the flexibility to adapt its
frequency content to the data, concentrating representational capacity
on the frequencies actually present in the signal.

#### Bayesian regularization via the marginal posterior

Fitting the observation model uses the Bayesian framework of Foresee &
Hagan (1997) [\[4\]](#ref4), which replaces the standard log-posterior
with an approximation to the log marginal posterior that integrates out
the noise variance:

\\\log P(\omega \mid Y) \approx -\frac{T}{2}
\log\Bigl(\frac{1}{2}\sum\_{n=1}^{T} \varepsilon_n^2 + 1\Bigr) - c \cdot
\frac{M}{2} \log\Bigl(\frac{1}{2}\sum\_{j=1}^{M} \omega_j^2 + 1\Bigr)\\

where \\\varepsilon_n = Y_n - \hat{y}(t_n, \omega)\\ are residuals, \\M
= 3H\\ is the number of parameters, and \\c\\ controls the
regularization strength. The first term is a robust alternative to the
sum-of-squares likelihood (it grows logarithmically rather than linearly
with the residual magnitude, making it less sensitive to outliers). The
second term is a spherical Gaussian prior on the parameters, preventing
overfitting.

The default regularization \\c = 1/H\\ (inversely proportional to the
number of neurons) is data-independent and works well across a wide
range of systems (Bonnaffé 2023). Stronger regularization (smaller
\\c\\) produces smoother interpolations at the cost of fitting fidelity;
weaker regularization (larger \\c\\) allows the network to track
high-frequency fluctuations that may be noise.

### Stage 2: dual-pathway process model

After interpolation, the per-capita growth rate for each variable is
computed:

\\r_i(t) = \frac{1}{\hat{y}\_i(t)} \cdot \frac{d\hat{y}\_i}{dt}(t) =
\frac{d}{dt}\bigl\[\log \hat{y}\_i(t)\bigr\]\\

This transformation removes the scaling effect of population size and
isolates the intrinsic dynamics. The process model then predicts \\r_i\\
as a function of all state variables using a dual-pathway single-layer
perceptron (SLP):

\\f\_{p,i}(\mathbf{x}) = \underbrace{\mathbf{v}\_{\text{lin}}^\top
\bigl(\mathbf{b}\_{\text{lin}} + \mathbf{W}\_{\text{lin}}
\mathbf{x}\bigr)}\_{\text{linear pathway}} +
\underbrace{\mathbf{v}\_{\exp}^\top \exp\bigl(\mathbf{b}\_{\exp} +
\mathbf{W}\_{\exp} \mathbf{x}\bigr)}\_{\text{exponential pathway}}\\

Each pathway has \\H_p\\ neurons with \\N\\ inputs (one per state
variable), giving \\H_p(2 + N)\\ parameters per pathway and \\2 H_p(2 +
N)\\ total parameters per target variable. The linear pathway captures
density-independent effects and weak linear interactions, while the
exponential pathway captures nonlinear density-dependent interactions
such as predation saturation, competitive exclusion, and Allee effects.
This separation improves both interpretability (one can examine each
pathway’s contribution independently) and fitting stability (the linear
pathway provides a baseline that the exponential pathway refines).

The process model uses a standard Gaussian log-posterior:

\\\log P(\theta \mid \mathbf{X}, \mathbf{r}) = -\sum\_{t}
\frac{(r_i(t) - f\_{p,i}(\mathbf{x}(t)))^2}{\sigma_1^2} - \sum_j
\frac{\theta_j^2}{\sigma_2^2}\\

where \\\sigma_1\\ controls the likelihood width (data fit tolerance)
and \\\sigma_2\\ controls the prior width (parameter regularization).
Unlike the observation model, the process model uses a train/test split:
parameters are optimized on the training fraction, and test-set
prediction error is used to rank ensemble members.

### Anchored Bayesian ensembles

Both stages use anchored ensemble sampling (Pearce et al. 2018
[\[2\]](#ref2)) for uncertainty quantification without MCMC. The
procedure is:

1.  Draw \\K\\ random parameter initializations from \\\mathcal{N}(0,
    0.001^2 I)\\
2.  Optimize each initialization to a (local) posterior maximum via BFGS
3.  Retain the best \\\rho K\\ members ranked by marginal likelihood
    (stage 1) or test loss (stage 2)

The retained ensemble provides approximate posterior uncertainty through
the spread of the members. This is computationally cheaper than MCMC but
provides only a point-estimate approximation to the posterior. For the
observation model, an acceptance criterion ensures quality: only
initializations that improve the log-posterior by at least 10 units are
retained, discarding poor local optima.

### Interaction inference: Jacobian and Geber method

The Jacobian matrix \\\mathbf{J}\\ with entries

\\e\_{ij}(t) = \frac{\partial f\_{p,i}}{\partial x_j}(\mathbf{x}(t))\\

quantifies the instantaneous effect of variable \\j\\ on the per-capita
growth rate of variable \\i\\ at time \\t\\. Positive values indicate
facilitation (mutualism, predation benefit to the consumer); negative
values indicate suppression (competition, predation cost to the prey).
The time-averaged Jacobian \\\bar{e}\_{ij} = T^{-1} \sum_t e\_{ij}(t)\\
gives the mean effect.

However, a strong effect does not necessarily translate into a strong
contribution to the observed dynamics. A large \\e\_{ij}\\ has little
impact if variable \\j\\ barely changes over the observation period. The
Geber method (Hairston et al. 2005 [\[5\]](#ref5)) addresses this by
weighting effects by the temporal variation they explain:

\\c\_{ij} = \frac{\sum_t \bigl(\dot{x}\_j(t) \cdot
e\_{ij}(t)\bigr)^2}{\sum_k \sum_t \bigl(\dot{x}\_k(t) \cdot
e\_{ik}(t)\bigr)^2}\\

These contributions are normalized so that \\\sum_j c\_{ij} = 1\\ for
each target variable \\i\\, giving a proportional decomposition of the
drivers of each variable’s dynamics. Contributions below a
user-specified threshold (default 0.1) are set to zero, yielding the
thresholded adjacency matrix that defines the inferred interaction
network.

### Forecasting via RK4 integration

Once the process model ensemble is fitted, forward simulation is
possible via classical fourth-order Runge-Kutta (RK4) integration. The
integration operates in log space (\\Z_i = \log Y_i\\) where the
dynamics are:

\\\frac{dZ_i}{dt} = r_i = \mu\_{r,i} + \sigma\_{r,i} \cdot
f\_{p,i}\Bigl(\frac{Z - \mu_x}{\sigma_x}\Bigr)\\

using the standardization coefficients from the fitting stage. Forecasts
inherit uncertainty from the process model ensemble: each ensemble
member produces a different trajectory, and the spread indicates
forecast uncertainty. A clamping mechanism prevents divergence when the
fitted dynamics have unrealistic eigenvalues at states far from the
training data.

## Example 1: Lotka-Volterra predator-prey

The classic two-species predator-prey model provides a ground-truth
benchmark where the interaction structure is known analytically:

\\\frac{dN}{dt} = rN(1 - \alpha P), \qquad \frac{dP}{dt} = P(-\delta +
\beta N)\\

The Jacobian of per-capita growth rates is:

\\\mathbf{J} = \begin{pmatrix} \partial r_N/\partial N & \partial
r_N/\partial P \\ \partial r_P/\partial N & \partial r_P/\partial P
\end{pmatrix} = \begin{pmatrix} 0 & -r\alpha \\ \beta & 0
\end{pmatrix}\\

We expect: zero diagonal (no self-regulation in this simple form),
negative predator-to-prey effect, and positive prey-to-predator effect.

``` r

library(lucifer)
set.seed(666)

n <- 100; dt <- 0.1
prey <- predator <- numeric(n)
prey[1] <- 1.0; predator[1] <- 0.5
for (i in 2:n) {
    prey[i] <- prey[i-1] + dt * prey[i-1] * (1 - 0.5 * predator[i-1])
    predator[i] <- predator[i-1] + dt * predator[i-1] * (-0.5 + 0.3 * prey[i-1])
}
data_lv <- cbind(prey = prey, predator = predator)
times_lv <- seq(dt, n * dt, by = dt)
```

Fitting with the default `obs.neurons = "auto"` lets the method
automatically select the number of sinusoidal neurons by evaluating a
grid of candidates via a fast pilot process model fit. The `times`
argument must match the actual observation spacing so that temporal
derivatives and per-capita growth rates are computed in the correct
units:

``` r

fit_lv <- NODE(
    data_lv,
    times = times_lv,
    obs.neurons = 10L,
    # obs.neurons = "auto",
    obs.ensemble = 50L,
    obs.retain = 0.2,
    proc.neurons = 5L,
    proc.ensemble = 50L,
    proc.retain = 0.2,
    bfgs.maxiter = 300L
)
print(fit_lv)
```

The auto-selection metadata shows which candidates were evaluated and
their pilot process R-squared:

``` r

sel <- fit_lv$obs.neurons.selection
if (!is.null(sel)) {
    knitr::kable(data.frame(
        H = sel$grid,
        `min proc R2` = round(sel$min_scores, 3),
        selected = ifelse(sel$grid == sel$selected, "\u2714", "")
    ), align = "rcl", caption = "Obs.neurons auto-selection scores")
}
```

The interpolation plot shows the sinusoidal network’s fit to the raw
data (left column) and the extracted temporal derivatives (right
column):

``` r

plot(fit_lv, type = "interpolation")
```

The dynamics plot compares observed per-capita growth rates (grey
points, computed from the interpolated derivatives) with the process
model’s fitted values (colored lines). High R-squared values indicate
the dual-pathway SLP successfully captures the functional relationship
between states and growth rates:

``` r

plot(fit_lv, type = "dynamics")
```

The effects heatmap reveals the recovered Jacobian. The expected pattern
(prey suppressed by predator, predator facilitated by prey) should be
clearly visible with the correct signs:

``` r

plot(fit_lv, type = "effects")
```

We can compare the recovered Jacobian against the analytical ground
truth. With \\r = 1\\, \\\alpha = 0.5\\, \\\beta = 0.3\\, the true
Jacobian is \\\partial r_N / \partial P = -0.5\\ and \\\partial r_P /
\partial N = 0.3\\:

``` r

true_J <- matrix(c(0, -0.5, 0.3, 0), 2, 2,
                 dimnames = list(c("prey", "predator"), c("prey", "predator")))
recovered_J <- round(fit_lv$effects_raw, 4)

# Percent error on non-zero entries
mask <- true_J != 0
pct_err <- round(abs(recovered_J[mask] - true_J[mask]) / abs(true_J[mask]) * 100, 1)

comparison <- data.frame(
    Interaction = c("prey self", "pred \u2192 prey", "prey \u2192 pred", "pred self"),
    True = as.vector(true_J),
    Recovered = as.vector(t(recovered_J)),
    stringsAsFactors = FALSE
)
comparison$Error <- ifelse(comparison$True == 0, "\u2014",
    paste0(round(abs(comparison$Recovered - comparison$True) /
                  pmax(abs(comparison$True), 1e-10) * 100, 1), "%"))

knitr::kable(comparison, align = "lrrr",
             caption = "Ground truth vs recovered Jacobian entries")
```

The Geber contributions decompose the drivers of each variable’s
dynamics. For prey, the predator contribution should dominate (the
prey’s growth rate varies mainly because predator density changes). For
predator, prey contribution should dominate:

``` r

plot(fit_lv, type = "contributions")
```

The network diagram provides a visual summary. Green edges indicate
positive effects (facilitation), red edges indicate negative effects
(suppression), and edge width is proportional to the Geber contribution:

``` r

plot(fit_lv, type = "network")
```

Phase portraits overlay observed and fitted trajectories in state space:

``` r

plot(fit_lv, type = "phase")
```

Forecasting via RK4 integration produces trajectories that continue the
oscillatory dynamics beyond the observed window. `NODE_predict` returns
a `node_forecast` object with its own `print` and `plot` methods:

``` r

pred_lv <- NODE_predict(fit_lv, horizon = 3, n.steps = 60)
print(pred_lv)
```

The trajectory plot overlays observed data, the fitted interpolation,
and the forecast extension. A shaded region marks the forecast zone, and
a diamond marker indicates the forecast origin:

``` r

plot(pred_lv)
```

The phase space view shows how the forecast continues the fitted orbit.
The diamond marks the forecast start; the triangle marks its endpoint:

``` r

plot(pred_lv, type = "phase")
```

## Example 2: three-species food web

Ecological communities typically involve more than two species arranged
in trophic levels. The tri-trophic food chain, with a resource (R),
consumer (C), and predator (P), is a canonical system whose interaction
structure is well understood:

\\\dot{R} = R(r_R - a\_{RR}R - a\_{RC}C), \quad \dot{C} = C(-d_C +
a\_{CR}R - a\_{CP}P), \quad \dot{P} = P(-d_P + a\_{PC}C)\\

The expected Jacobian has R self-regulation (negative diagonal),
top-down control of C on R (negative), bottom-up flow from R to C
(positive), top-down control of P on C (negative), and bottom-up flow
from C to P (positive). The predator P has no direct effect on R (only
indirect, mediated through C).

``` r

set.seed(666)
n <- 200; dt <- 0.1
R <- C <- P <- numeric(n)
R[1] <- 2.0; C[1] <- 1.0; P[1] <- 0.5
for (i in 2:n) {
    R[i] <- max(0.01, R[i-1] + dt * R[i-1] * (1.5 - 0.3*R[i-1] - 0.8*C[i-1]))
    C[i] <- max(0.01, C[i-1] + dt * C[i-1] * (-0.5 + 0.5*R[i-1] - 0.6*P[i-1]))
    P[i] <- max(0.01, P[i-1] + dt * P[i-1] * (-0.3 + 0.3*C[i-1]))
}
data_fw <- cbind(R = R, C = C, P = P)
times_fw <- seq(dt, n * dt, by = dt)

fit_fw <- NODE(
    data_fw,
    times = times_fw,
    obs.ensemble = 80L,
    obs.retain = 0.15,
    proc.neurons = 6L,
    proc.ensemble = 80L,
    proc.retain = 0.15,
    threshold = 0.05,
    bfgs.maxiter = 400L
)
print(fit_fw)
```

Ground truth comparison for the food web. The true per-capita growth
rate Jacobian has entries \\\partial r_R/\partial R = -0.3\\, \\\partial
r_R/\partial C = -0.8\\, \\\partial r_C/\partial R = 0.5\\, \\\partial
r_C/\partial P = -0.6\\, \\\partial r_P/\partial C = 0.3\\:

``` r

true_fw <- matrix(c(-0.3, -0.8, 0, 0.5, 0, -0.6, 0, 0.3, 0), 3, 3,
                   dimnames = list(c("R","C","P"), c("R","C","P")))
rec_fw <- round(fit_fw$effects_raw, 4)

entries <- data.frame(
    Interaction = c("R self", "C \u2192 R", "R \u2192 C",
                     "P \u2192 C", "C \u2192 P"),
    True = c(-0.3, -0.8, 0.5, -0.6, 0.3),
    Recovered = c(rec_fw["R","R"], rec_fw["R","C"], rec_fw["C","R"],
                   rec_fw["C","P"], rec_fw["P","C"]),
    stringsAsFactors = FALSE
)
entries$Error <- paste0(round(abs(entries$Recovered - entries$True) /
                               pmax(abs(entries$True), 1e-10) * 100, 1), "%")
knitr::kable(entries, align = "lrrr",
             caption = "Food web: ground truth vs recovered Jacobian")
```

The network visualization reveals the trophic cascade: R feeds C, C
feeds P, and top-down control propagates back through the chain. The
absence of a direct R-P link tests the method’s ability to distinguish
direct from indirect interactions.

``` r

plot(fit_fw, type = "network")
```

``` r

plot(fit_fw, type = "effects")
```

The dynamics plot for three species shows the quality of the process
model fit across all trophic levels:

``` r

plot(fit_fw, type = "dynamics")
```

Forecasting the tri-trophic system forward in time tests whether the
fitted NODE captures the qualitative long-term behavior (damped
oscillations toward the coexistence equilibrium):

``` r

pred_fw <- NODE_predict(fit_fw, horizon = 5, n.steps = 80)
print(pred_fw)
plot(pred_fw)
```

## Example 3: four-species community with trophic structure

To test the method’s scalability, we construct a four-species community
with two prey species (grass G and herb H) and two predator species (a
generalist predator F1 that feeds on both prey, and a specialist
predator F2 that feeds exclusively on H). This system generates the
following expected interaction matrix:

- G is limited by self-regulation, competition with H, and predation by
  F1
- H is limited by competition with G, predation by F1, and heavy
  predation by F2
- F1 benefits from both G and H (generalist)
- F2 benefits only from H (specialist)

``` r

set.seed(666)
n <- 400; dt <- 0.05
G <- H <- F1 <- F2 <- numeric(n)
G[1] <- 3.0; H[1] <- 2.5; F1[1] <- 1.0; F2[1] <- 0.8

for (i in 2:n) {
    dG  <- G[i-1] * (1.2 - 0.2*G[i-1] - 0.15*H[i-1] - 0.4*F1[i-1])
    dH  <- H[i-1] * (1.0 - 0.1*G[i-1] - 0.25*H[i-1] - 0.2*F1[i-1] - 0.5*F2[i-1])
    dF1 <- F1[i-1] * (-0.4 + 0.2*G[i-1] + 0.15*H[i-1])
    dF2 <- F2[i-1] * (-0.5 + 0.4*H[i-1])
    G[i]  <- max(0.01, G[i-1] + dt * dG)
    H[i]  <- max(0.01, H[i-1] + dt * dH)
    F1[i] <- max(0.01, F1[i-1] + dt * dF1)
    F2[i] <- max(0.01, F2[i-1] + dt * dF2)
}
data_4sp <- cbind(G = G, H = H, F1 = F1, F2 = F2)
times_4sp <- seq(dt, n * dt, by = dt)

fit_4sp <- NODE(
    data_4sp,
    times = times_4sp,
    obs.ensemble = 80L,
    obs.retain = 0.1,
    proc.neurons = 6L,
    proc.ensemble = 80L,
    proc.retain = 0.1,
    threshold = 0.05,
    bfgs.maxiter = 400L
)
print(fit_4sp)
```

Ground truth comparison for the 4-species system. The true Jacobian has
10 non-zero entries encoding self-regulation, competition, and
predation:

``` r

true_4sp <- matrix(c(
    -0.2, -0.15, -0.4,  0,
    -0.1, -0.25, -0.2, -0.5,
     0.2,  0.15,  0,    0,
     0,    0.4,   0,    0
), 4, 4, byrow = TRUE,
dimnames = list(c("G","H","F1","F2"), c("G","H","F1","F2")))
rec_4sp <- fit_4sp$effects_raw

entries_4sp <- data.frame(
    Interaction = c("G self", "H \u2192 G", "F1 \u2192 G",
                     "G \u2192 H", "H self", "F1 \u2192 H", "F2 \u2192 H",
                     "G \u2192 F1", "H \u2192 F1", "H \u2192 F2"),
    True = c(-0.2, -0.15, -0.4, -0.1, -0.25, -0.2, -0.5, 0.2, 0.15, 0.4),
    Recovered = c(rec_4sp["G","G"], rec_4sp["G","H"], rec_4sp["G","F1"],
                   rec_4sp["H","G"], rec_4sp["H","H"], rec_4sp["H","F1"],
                   rec_4sp["H","F2"], rec_4sp["F1","G"], rec_4sp["F1","H"],
                   rec_4sp["F2","H"]),
    stringsAsFactors = FALSE
)
entries_4sp$Recovered <- round(entries_4sp$Recovered, 4)
entries_4sp$`Sign OK` <- ifelse(sign(entries_4sp$True) == sign(entries_4sp$Recovered),
                                 "\u2714", "\u2718")
knitr::kable(entries_4sp, align = "lrrrc",
             caption = "4-species community: ground truth vs recovered Jacobian")
```

The four-species network diagram displays the trophic structure, with
bottom-up positive edges (green) from prey to predators and top-down
negative edges (red) from predators to prey:

``` r

plot(fit_4sp, type = "network")
```

The effects heatmap reveals the full interaction matrix, including weak
indirect effects that may be below the contribution threshold:

``` r

plot(fit_4sp, type = "effects")
```

``` r

plot(fit_4sp, type = "contributions")
```

The interpolation quality for all four species confirms that the
sinusoidal network captures the oscillatory dynamics with high fidelity:

``` r

plot(fit_4sp, type = "interpolation")
```

The four-species phase space forecast shows the multi-dimensional
trajectory continuing from the fitted orbit. With four variables, phase
plots are generated for all unique pairs:

``` r

pred_4sp <- NODE_predict(fit_4sp, horizon = 3, n.steps = 60)
print(pred_4sp)
plot(pred_4sp)
```

``` r

plot(pred_4sp, type = "phase")
```

## Example 4: noisy observations

Real ecological time series are contaminated by observation error,
measurement imprecision, and stochastic environmental fluctuations.
Adding multiplicative log-normal noise to the Lotka-Volterra system
tests the method’s robustness. Because the derivative computation
amplifies noise (numerical differentiation is an ill-conditioned
operation), even moderate observation error can substantially degrade
the process model’s ability to recover interaction structure. This
example uses 5% coefficient of variation, which is at the lower end of
ecological monitoring data:

``` r

set.seed(666)
noise_level <- 0.05
data_noisy <- data_lv * exp(matrix(rnorm(prod(dim(data_lv)), 0, noise_level),
                                    nrow = nrow(data_lv)))

fit_noisy <- NODE(
    data_noisy,
    times = times_lv,
    obs.ensemble = 60L,
    obs.retain = 0.15,
    proc.neurons = 5L,
    proc.ensemble = 60L,
    proc.retain = 0.15,
    bfgs.maxiter = 300L
)
print(fit_noisy)
```

The ground truth comparison reveals how observation noise degrades
interaction recovery. The signs may still be correct, but magnitudes
diverge more than in the clean case:

``` r

rec_noisy <- round(fit_noisy$effects_raw, 4)
comparison_noisy <- data.frame(
    Interaction = c("prey self", "pred \u2192 prey", "prey \u2192 pred", "pred self"),
    True = c(0, -0.5, 0.3, 0),
    Clean = as.vector(t(round(fit_lv$effects_raw, 4))),
    Noisy = as.vector(t(rec_noisy)),
    stringsAsFactors = FALSE
)
comparison_noisy$`Sign OK` <- ifelse(
    comparison_noisy$True == 0, "\u2014",
    ifelse(sign(comparison_noisy$True) == sign(comparison_noisy$Noisy), "\u2714", "\u2718"))
knitr::kable(comparison_noisy, align = "lrrrr",
             caption = "Effect of 5% observation noise on Jacobian recovery")
```

The interpolation shows smooth fits that capture the underlying dynamics
while filtering observation noise. The sinusoidal basis acts as a
natural low-pass filter, with the regularization parameter controlling
the cutoff frequency:

``` r

plot(fit_noisy, type = "interpolation")
```

Despite the noise, the recovered interaction structure should preserve
the correct sign pattern, though with weaker magnitudes and larger
ensemble uncertainty compared to the noise-free case:

``` r

plot(fit_noisy, type = "effects")
```

The summary provides ensemble spread information to assess confidence in
the recovered interactions:

``` r

summary(fit_noisy)
```

## Practical guidelines

### Time scale and the `times` argument

The `times` argument controls the units in which derivatives and
per-capita growth rates are computed. If omitted,
[`NODE()`](https://robustecologies.github.io/lucifer/reference/NODE.md)
defaults to `times = 1:nrow(data)`, treating each row as one time unit
apart. This means the recovered Jacobian entries will be in units of
“per row” rather than “per second” or “per year”. For correct physical
interpretation of interaction effects, always pass the actual
observation times; the spacing `diff(times[1:2])` determines the
derivative scaling factor. For instance, if the simulation uses
`dt = 0.1`, passing `times = seq(0.1, 10, by = 0.1)` yields effects in
the correct continuous-time units, whereas omitting `times` would shrink
them by a factor of 10.

### Choosing observation model neurons

The sinusoidal network’s capacity must match the temporal complexity of
the signal. Too few neurons produce overly smooth interpolations that
miss real dynamics; too many allow the network to fit noise. A
reasonable starting point is `obs.neurons` \\\approx T/5\\ to \\T/3\\,
but the built-in regularization (\\c = 1/H\\) provides automatic
complexity control regardless of the neuron count. For oscillatory
systems, ensure enough neurons to represent the dominant frequencies: a
signal with \\k\\ distinct periodicities needs at least \\2k\\ neurons.

### Choosing process model neurons

The dual-pathway architecture already separates linear from nonlinear
effects, so moderate neuron counts (5-10 per pathway) suffice for most
ecological systems. The total parameter count per target variable is
\\2H_p(2 + N)\\, which grows linearly with the number of state variables
\\N\\. For a 4-species system with \\H_p = 6\\, this gives \\2 \times 6
\times 6 = 72\\ parameters, well within the capacity of BFGS to
optimize.

### Ensemble sizing

Larger ensembles improve coverage of the loss surface but increase
computation linearly. The key tradeoff is between ensemble size \\K\\
and retention fraction \\\rho\\. A conservative choice is \\K = 100\\,
\\\rho = 0.1\\ (keeping 10 members). For quick exploration, \\K = 30\\,
\\\rho = 0.3\\ (keeping 9) works well. The retained members should have
comparable loss values; if the best and worst retained members differ
substantially, increase \\K\\ to improve coverage.

### Process model regularization

The `proc.sd1` parameter controls the likelihood width: larger values
tolerate more process error, smaller values demand tighter fits. For
clean simulated data, `proc.sd1 = 0.1` works well. For noisy real data,
increase to 0.15-0.3 to prevent overfitting the observation model’s
interpolation artifacts. The prior width `proc.sd2` controls parameter
shrinkage; setting it equal to `proc.sd1` (the default) provides
balanced regularization.

### Interpreting results

The effects matrix (Jacobian) gives the direction and magnitude of each
pairwise interaction. However, magnitude comparison across variables
requires caution: effects are in the units of per-capita growth rate per
unit change in the (log-transformed, standardized) state variable. The
Geber contributions provide a more robust basis for comparison because
they account for how much temporal variation each interaction actually
explains.

The threshold parameter controls network sparsity. Lower thresholds
(0.01-0.05) retain more edges, including weak indirect effects; higher
thresholds (0.1-0.2) produce sparser networks focusing on dominant
interactions. For exploratory analysis, start with a low threshold and
increase gradually.

### When to use NODE vs SDE vs SBI

[`NODE()`](https://robustecologies.github.io/lucifer/reference/NODE.md)
is the right choice when the functional form of interactions is unknown
and the primary goal is recovering the qualitative interaction structure
(who affects whom, with what sign) from time series data. For systems
where a mechanistic model is available and the goal is parameter
estimation with proper uncertainty quantification,
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md) or
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
with MCMC are more appropriate.
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md) is
best when the likelihood is intractable but a forward simulator is
available.

### Limitations

The BNGM method assumes all dynamically important state variables are
observed. Unobserved latent variables (e.g., nutrient concentrations,
temperature) can bias the recovered interaction structure by forcing the
network to attribute their effects to observed variables. Systems that
converge monotonically to equilibrium (e.g., purely competitive
communities approaching a stable coexistence point) provide little
temporal variation for the gradient matching to distinguish effects,
leading to lower process model R-squared and less reliable interaction
recovery. The method works best on systems with sustained oscillations
or transient dynamics that generate sufficient variation in all state
variables.

## Computational performance

The C++ backend with OpenMP parallelization provides substantial
speedups. On a modern multi-core system, typical runtimes are:

- 2-species, 100 time points, \\K = 50\\: \< 1 second
- 3-species, 200 time points, \\K = 80\\: 2-5 seconds
- 4-species, 400 time points, \\K = 80\\: 30-90 seconds

The dominant cost is the process model fitting, which scales as \\O(K
\cdot N \cdot I \cdot H_p^2 N^2)\\ where \\I\\ is the number of BFGS
iterations. For systems with more than 5-6 variables, consider reducing
the ensemble size or neuron count.

## References

**\[1\]** Bonnaffé, W. (2023). Fast fitting of neural ordinary
differential equations by Bayesian neural gradient matching. *Methods in
Ecology and Evolution*, 14(6), 1456-1468. [DOI:
10.1111/2041-210X.14121](https://doi.org/10.1111/2041-210X.14121)

**\[2\]** Pearce, T., Leibfried, F., & Brintrup, A. (2018). Uncertainty
in neural networks: approximately Bayesian ensembling. *Proceedings of
the International Conference on Artificial Intelligence and Statistics
(AISTATS)*.

**\[3\]** Chen, R.T.Q., Rubanova, Y., Bettencourt, J., & Duvenaud, D.K.
(2018). Neural ordinary differential equations. *Advances in Neural
Information Processing Systems*, 31.

**\[4\]** Foresee, F.D. & Hagan, M.T. (1997). Gauss-Newton approximation
to Bayesian learning. *Proceedings of the International Joint Conference
on Neural Networks (IJCNN)*.

**\[5\]** Hairston, N.G., Ellner, S.P., Geber, M.A., Yoshida, T., & Fox,
J.A. (2005). Rapid evolution and the convergence of ecological and
evolutionary time. *Ecology Letters*, 8(10), 1114-1127. [DOI:
10.1111/j.1461-0248.2005.00812.x](https://doi.org/10.1111/j.1461-0248.2005.00812.x)
