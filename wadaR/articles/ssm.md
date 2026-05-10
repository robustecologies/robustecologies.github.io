# State-space model inference in lucifer

## Introduction

State-space models (SSMs) constitute one of the most versatile and
widely deployed frameworks in time series analysis, providing a
principled way to represent systems where a latent process \\x\_{0:T}\\
evolves according to a Markov transition kernel and is observed
indirectly through noisy measurements \\y\_{1:T}\\. The general
formulation is

\\ x_t \mid x\_{t-1} \sim p(x_t \mid x\_{t-1}, \theta), \qquad y_t \mid
x_t \sim p(y_t \mid x_t, \theta), \\

where \\\theta\\ denotes the static parameters governing both the state
dynamics and the observation process. The first equation, often called
the state or transition equation, specifies how the latent process
evolves over time; the second, the observation or measurement equation,
specifies how the observed data are generated given the current latent
state. Together, these two equations define a hidden Markov model whose
inference requires marginalizing over the entire latent trajectory
\\x\_{0:T}\\, a high-dimensional integration problem that grows linearly
with the number of time steps.

The history of state-space modelling begins with Kalman’s landmark 1960
paper [\[17\]](#ref17), which derived the optimal linear recursive
filter for estimating the state of a linear dynamical system with
Gaussian noise. The Kalman filter provides the exact posterior mean and
covariance of the state given past and present observations, and its
computational efficiency, requiring only \\O(d_x^3)\\ operations per
time step, made it immediately practical for engineering applications in
navigation, control, and signal processing. The Bayesian interpretation
of the Kalman filter, which views the filtered state estimate as the
posterior distribution \\p(x_t \mid y\_{1:t}, \theta)\\, was developed
shortly thereafter and placed the filter within the broader framework of
Bayesian sequential updating. The extension to smoothing, where one
computes \\p(x_t \mid y\_{1:T}, \theta)\\ using the full observation
record, was achieved by Rauch, Tung, and Striebel [\[7\]](#ref7),
completing the linear-Gaussian inference toolkit.

The application of state-space methods to statistical time series
analysis, as distinct from engineering control theory, was pioneered by
Harvey [\[8\]](#ref8), who showed that classical decomposition models
for trend, seasonality, and cycles could be represented as
linear-Gaussian SSMs and estimated via maximum likelihood using the
prediction-error decomposition of the Kalman filter. Durbin and Koopman
[\[16\]](#ref16) extended this framework to handle non-Gaussian and
nonlinear models through importance sampling and simulation methods,
establishing the modern statistical perspective on state-space analysis.
The Bayesian approach to SSM inference, where both the parameters
\\\theta\\ and the latent states \\x\_{0:T}\\ are treated as random
variables with a joint posterior distribution, gained traction through
the work of Carter and Kohn [\[1\]](#ref1) and Fruhwirth-Schnatter
[\[2\]](#ref2), who independently developed the Forward-Filtering
Backward-Sampling (FFBS) algorithm for exact Gibbs sampling in
linear-Gaussian SSMs.

State-space models can be classified along several dimensions that
determine which inference algorithms are applicable and efficient. The
most fundamental distinction is between linear and nonlinear models:
when both the transition function \\f\\ and the observation function
\\h\\ are linear in the state, the Kalman filter provides exact
posterior inference; when either is nonlinear, approximate methods such
as the Extended Kalman Filter (EKF), the Unscented Kalman Filter (UKF),
or particle filters are required. A second distinction is between
Gaussian and non-Gaussian models: Gaussian state and observation noise
yield closed-form filtering distributions, while non-Gaussian noise
(e.g., Poisson observations, Student-t innovations) requires either
analytic approximations or sequential Monte Carlo methods. A third
distinction is between time-invariant and time-varying models: when the
system matrices \\F_t\\, \\Q_t\\, \\H_t\\, \\R_t\\ depend on time, the
Kalman filter generalizes naturally, but the parameter dimension may
increase with \\T\\ unless structural constraints are imposed. Finally,
regime-switching models introduce a discrete latent Markov chain that
selects among multiple sets of system matrices, requiring specialized
algorithms such as the Hamilton filter [\[13\]](#ref13) and its Kim
[\[14\]](#ref14) collapsing approximation.

Bayesian inference for SSMs requires sampling from the joint posterior
\\p(\theta, x\_{0:T} \mid y\_{1:T})\\, which involves both the
finite-dimensional parameters and the entire latent trajectory. This is
fundamentally different from standard posterior inference, where only
\\\theta\\ is unknown, and it cannot be handled by the
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
dispatcher, which operates on a flat parameter vector and returns a
`demonoid` object. The
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
function provides a dedicated interface for this class of problems: it
accepts either an `sde_model` object from
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md) or
an `ssm_model` object from
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md),
jointly samples parameters and latent states at each iteration, and
returns an `ssm_fit` object containing both the parameter posterior and
the smoothed state trajectories. The `ssm_fit` object can be converted
to a `demonoid` via
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
for downstream diagnostics such as
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md),
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
and
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
bridging the gap between the specialized SSM inference engine and the
general-purpose diagnostic toolkit.

The SSM module connects to other lucifer components in several ways. The
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md)
function provides the continuous-time model specification, with its
family registry handling discretization and parameterization; the
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
constructor provides a direct discrete-time specification route. Both
produce model objects that
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
can consume. After inference, the `ssm_fit` object integrates with the
LOO-PSIS framework for model comparison, with
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
for prior sensitivity analysis, and with the `model_spec` DSL for
specifying hierarchical extensions. The
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
function can also evaluate SSM model objects, comparing different
algorithm choices and structural specifications on the same dataset.

The SSM module provides eight algorithms spanning the full spectrum of
model complexity. For linear-Gaussian models, where both the transition
and observation densities are Gaussian with known functional forms, the
Forward-Filtering Backward-Sampling (FFBS) algorithm provides exact
Gibbs sampling via the Kalman filter. For nonlinear Gaussian models, the
Unscented Kalman Filter (UKF) propagates sigma points through the exact
nonlinear dynamics, avoiding the Jacobian linearization required by the
EKF; for high-dimensional state spaces where full-rank covariance is
intractable, the Ensemble Kalman Filter (EnKF) represents the filtering
distribution empirically via an ensemble of state particles. The
Rao-Blackwellized Particle Filter (RBPF) exploits mixed linear/nonlinear
structure to analytically marginalize the conditionally linear
substates, reducing variance without increasing particle count. For
general nonlinear or non-Gaussian models, Particle Gibbs with Ancestor
Sampling (PGAS) replaces the Kalman filter with a conditional particle
filter that handles arbitrary state dynamics. For the most general case,
where neither conjugacy nor tractable transition densities are
available, SMC-squared (SMC\\^2\\) nests a particle filter inside an SMC
sampler to simultaneously learn parameters and states without requiring
user-specified full conditionals. For regime-switching linear-Gaussian
models, Markov-switching FFBS (MS-FFBS) runs \\K\\ parallel Kalman
filters with Hamilton (1989) filtering and Kim (1994) collapsing
approximation. Finally, for the specific but important class of
stochastic volatility models, the Kim-Shephard-Chib (KSC) algorithm
exploits a mixture-of-normals approximation to convert the non-Gaussian
observation model into a conditionally linear-Gaussian system amenable
to FFBS.

When `Algorithm = "auto"` (the default),
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
inspects the SDE model specification and selects the most efficient
algorithm: KSC for stochastic volatility families, FFBS for exact
linear-Gaussian families with Gaussian observations, PGAS for general
nonlinear models, and SMC\\^2\\ as the fully general fallback. The UKF,
EnKF, RBPF, and MS-FFBS algorithms are available by explicit selection
when their structural assumptions hold.

## Mathematical foundations

This section develops the mathematical framework underlying all eight
SSM inference algorithms in lucifer, beginning with the general
nonlinear state-space model and specializing to the linear-Gaussian case
where exact inference is available. The notation established here will
be referenced throughout the algorithm descriptions that follow.

### General state-space model

The most general state-space model treated by lucifer has the form

\\ x_t = f(x\_{t-1}, \theta) + G(x\_{t-1}, \theta)\\ w_t, \qquad w_t
\sim \mathcal{N}(0, Q), \\

\\ y_t = h(x_t, \theta) + v_t, \qquad v_t \sim \mathcal{N}(0, R), \\

where \\x_t \in \mathbb{R}^{d_x}\\ is the latent state vector, \\y_t \in
\mathbb{R}^{d_y}\\ is the observation vector, \\f: \mathbb{R}^{d_x}
\times \Theta \to \mathbb{R}^{d_x}\\ is the state transition function,
\\h: \mathbb{R}^{d_x} \times \Theta \to \mathbb{R}^{d_y}\\ is the
observation function, \\G: \mathbb{R}^{d_x} \times \Theta \to
\mathbb{R}^{d_x \times d_w}\\ is the state-dependent diffusion matrix,
and the noise sequences \\\\w_t\\\\ and \\\\v_t\\\\ are mutually
independent. The initial state is drawn from \\x_0 \sim \mathcal{N}(m_0,
P_0)\\. This formulation encompasses both discrete-time models specified
directly via
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
and continuous-time SDEs discretized via
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md),
where the Euler-Maruyama scheme yields \\f(x, \theta) = x + \mu(x,
\theta) \Delta t\\ and \\G(x, \theta) = \sigma(x, \theta) \sqrt{\Delta
t}\\ for drift \\\mu\\ and diffusion \\\sigma\\.

The joint density of the complete data \\(x\_{0:T}, y\_{1:T})\\
factorizes as

\\ p(x\_{0:T}, y\_{1:T} \mid \theta) = p(x_0 \mid \theta) \prod\_{t=1}^T
p(x_t \mid x\_{t-1}, \theta) \prod\_{t=1}^T p(y_t \mid x_t, \theta), \\

and Bayesian inference targets the joint posterior \\p(\theta, x\_{0:T}
\mid y\_{1:T}) \propto p(\theta)\\ p(x\_{0:T}, y\_{1:T} \mid \theta)\\.
The marginal likelihood \\p(y\_{1:T} \mid \theta) = \int p(x\_{0:T},
y\_{1:T} \mid \theta) \\ dx\_{0:T}\\, obtained by integrating out the
latent states, plays a central role in model comparison and is available
as a byproduct of most filtering algorithms.

### Linear-Gaussian state-space models

When the transition and observation functions are both linear and the
noise is Gaussian, the model specializes to

\\ x_t = F_t x\_{t-1} + c_t + w_t, \qquad w_t \sim \mathcal{N}(0, Q_t),
\\

\\ y_t = H_t x_t + d_t + v_t, \qquad v_t \sim \mathcal{N}(0, R_t), \\

where \\F_t \in \mathbb{R}^{d_x \times d_x}\\ is the state transition
matrix, \\H_t \in \mathbb{R}^{d_y \times d_x}\\ is the observation
matrix, \\c_t\\ and \\d_t\\ are optional intercept vectors, and \\Q_t\\
and \\R_t\\ are the state and observation noise covariance matrices
respectively. The initial state is \\x_0 \sim \mathcal{N}(m_0, P_0)\\.
For time-invariant models the subscript \\t\\ is dropped from all system
matrices; in lucifer, time-varying matrices are passed as 3D arrays
(sliced along the third dimension) while time-invariant matrices are
passed as 2D matrices that the C++ backend broadcasts across all time
steps.

The Kalman filter provides the exact posterior distribution \\p(x_t \mid
y\_{1:t}, \theta) = \mathcal{N}(m_t, P_t)\\ via the following recursive
prediction-update cycle. The prediction step propagates the filtered
distribution one step forward:

\\ m\_{t\|t-1} = F_t m\_{t-1} + c_t, \qquad P\_{t\|t-1} = F_t P\_{t-1}
F_t^\top + Q_t. \\

The innovation, which is the deviation of the observation from its
predicted value, and its covariance are

\\ \nu_t = y_t - H_t m\_{t\|t-1} - d_t, \qquad S_t = H_t P\_{t\|t-1}
H_t^\top + R_t. \\

The Kalman gain determines the optimal weighting between the prediction
and the innovation:

\\ K_t = P\_{t\|t-1} H_t^\top S_t^{-1}. \\

The update step incorporates the observation:

\\ m_t = m\_{t\|t-1} + K_t \nu_t. \\

For the covariance update, lucifer uses the Joseph stabilized form
rather than the standard \\P_t = (I - K_t H_t) P\_{t\|t-1}\\, which can
lose positive semi-definiteness due to floating-point arithmetic:

\\ P_t = (I - K_t H_t) P\_{t\|t-1} (I - K_t H_t)^\top + K_t R_t
K_t^\top. \\

The Joseph form is algebraically equivalent to the standard update but
is guaranteed to produce a symmetric positive semi-definite covariance
matrix even in the presence of rounding errors, which is critical for
numerical stability in long time series or ill-conditioned models.

The log-likelihood of the parameters given the observations, computed as
a byproduct of the Kalman filter, takes the prediction-error
decomposition form:

\\ \ell(\theta) = \log p(y\_{1:T} \mid \theta) = -\frac{1}{2}
\sum\_{t=1}^T \bigl\[ d_y \log 2\pi + \log \|S_t\| + \nu_t^\top S_t^{-1}
\nu_t \bigr\]. \\

This expression decomposes the joint likelihood into a product of
one-step-ahead predictive densities \\p(y_t \mid y\_{1:t-1}, \theta) =
\mathcal{N}(\nu_t; 0, S_t)\\, each of which is computed during the
forward filter pass. The prediction-error decomposition is the
foundation of maximum likelihood estimation for linear-Gaussian SSMs and
provides the marginal likelihood estimate used for model comparison in
the Bayesian setting.

The Rauch-Tung-Striebel (RTS) smoother [\[7\]](#ref7) computes the
smoothed distribution \\p(x_t \mid y\_{1:T}, \theta) =
\mathcal{N}(m_t^s, P_t^s)\\ by running a backward pass after the forward
filter:

\\ G_t = P_t F\_{t+1}^\top P\_{t+1\|t}^{-1}, \qquad m_t^s = m_t + G_t
(m\_{t+1}^s - m\_{t+1\|t}), \qquad P_t^s = P_t + G_t (P\_{t+1}^s -
P\_{t+1\|t}) G_t^\top. \\

The backward simulation variant used in FFBS replaces the RTS smoother’s
mean computation with a draw from \\\mathcal{N}(m_t^s, P_t^s)\\,
producing a sample from the full joint smoothing distribution
\\p(x\_{0:T} \mid y\_{1:T}, \theta)\\ rather than the marginal smoothed
mean at each time point.

### Common univariate state-space models

Several widely used time series models are special cases of the
linear-Gaussian SSM, differing only in the structure of the system
matrices.

The **local level model** (also called the random walk plus noise model)
is the simplest structural time series model, with a single state
representing a stochastically evolving level:

\\ x_t = x\_{t-1} + w_t, \qquad y_t = x_t + v_t, \qquad w_t \sim
\mathcal{N}(0, \sigma_w^2), \quad v_t \sim \mathcal{N}(0, \sigma_v^2).
\\

The system matrices are \\F = 1\\, \\Q = \sigma_w^2\\, \\H = 1\\, \\R =
\sigma_v^2\\, all scalar. The signal-to-noise ratio \\q = \sigma_w^2 /
\sigma_v^2\\ determines the smoothness of the estimated level: small
\\q\\ produces a smooth level close to a simple average, while large
\\q\\ tracks the observations closely.

The **local linear trend model** augments the local level with a
stochastic slope component:

\\ \begin{pmatrix} \mu_t \\ \beta_t \end{pmatrix} = \begin{pmatrix} 1 &
1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} \mu\_{t-1} \\ \beta\_{t-1}
\end{pmatrix} + \begin{pmatrix} w\_{1t} \\ w\_{2t} \end{pmatrix}, \qquad
y_t = \begin{pmatrix} 1 & 0 \end{pmatrix} \begin{pmatrix} \mu_t \\
\beta_t \end{pmatrix} + v_t, \\

where \\\mu_t\\ is the level and \\\beta_t\\ is the slope. Setting
\\\sigma\_{w_2}^2 = 0\\ (fixed slope) yields the integrated random walk
model; setting both state noise variances to zero recovers deterministic
linear trend.

The **seasonal model** with period \\s\\ uses a dummy-variable rotation
that constrains the seasonal effects to sum to zero over a complete
cycle. The state vector has dimension \\s - 1\\, with transition matrix

\\ F = \begin{pmatrix} -1 & -1 & \cdots & -1 \\ 1 & 0 & \cdots & 0 \\ 0
& 1 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & 0
\end{pmatrix}, \\

where the first row ensures that the sum of all \\s\\ seasonal
components (current plus \\s-1\\ lagged) is driven toward zero, and the
remaining rows shift the seasonal components forward. The observation
matrix is \\H = (1, 0, \ldots, 0)\\, selecting the current seasonal
effect.

The **basic structural model** (BSM), developed systematically by Harvey
[\[8\]](#ref8), combines local linear trend and seasonal components into
a single model. The state vector is \\(\mu_t, \beta_t, \gamma\_{1,t},
\ldots, \gamma\_{s-1,t})^\top\\ with dimension \\2 + s - 1\\, and the
system matrices are block-diagonal compositions of the trend and
seasonal sub-models. In lucifer, the
`SSM_model(type = "bsm", period = s)` builder constructs all system
matrices and provides conjugate inverse-gamma full conditionals for the
three variance parameters \\(\sigma\_\mu^2, \sigma\_\beta^2,
\sigma\_\gamma^2)\\ and the observation noise variance \\\sigma_v^2\\.

### Common multivariate state-space models

The SSM framework extends naturally to multivariate settings. A vector
autoregressive model of order \\p\\, denoted VAR(\\p\\), can be cast as
a linear-Gaussian SSM through the companion form representation.
Consider \\y_t \in \mathbb{R}^{d_y}\\ following

\\ y_t = A_1 y\_{t-1} + A_2 y\_{t-2} + \cdots + A_p y\_{t-p} + u_t,
\qquad u_t \sim \mathcal{N}(0, \Sigma), \\

where each \\A_j \in \mathbb{R}^{d_y \times d_y}\\ is a coefficient
matrix. Defining the stacked state \\x_t = (y_t^\top, y\_{t-1}^\top,
\ldots, y\_{t-p+1}^\top)^\top \in \mathbb{R}^{d_y p}\\, the companion
form is

\\ x_t = \underbrace{\begin{pmatrix} A_1 & A_2 & \cdots & A\_{p-1} & A_p
\\ I & 0 & \cdots & 0 & 0 \\ 0 & I & \cdots & 0 & 0 \\ \vdots & & \ddots
& & \vdots \\ 0 & 0 & \cdots & I & 0 \end{pmatrix}}\_{F} x\_{t-1} +
\begin{pmatrix} u_t \\ 0 \\ \vdots \\ 0 \end{pmatrix}, \qquad y_t =
\begin{pmatrix} I & 0 & \cdots & 0 \end{pmatrix} x_t. \\

The `SSM_model(type = "var_p", data = Y, lags = p)` builder constructs
this companion representation automatically for multivariate data.

The **dynamic factor model** (DFM) explains the covariation among a
potentially large number of observed series \\y_t \in \mathbb{R}^{d_y}\\
through a small number of latent factors \\f_t \in \mathbb{R}^{d_f}\\
with \\d_f \ll d_y\\:

\\ y_t = \Lambda f_t + e_t, \qquad f_t = \Phi f\_{t-1} + \eta_t, \\

where \\\Lambda \in \mathbb{R}^{d_y \times d_f}\\ is the factor loading
matrix, \\\Phi \in \mathbb{R}^{d_f \times d_f}\\ is the factor
transition matrix, \\e_t \sim \mathcal{N}(0, D)\\ with \\D\\ typically
diagonal, and \\\eta_t \sim \mathcal{N}(0, Q_f)\\. This is a
linear-Gaussian SSM with \\H = \Lambda\\, \\F = \Phi\\, \\R = D\\, and
\\Q = Q_f\\, and can be estimated via FFBS with appropriate
identification constraints on \\\Lambda\\ (typically lower-triangular
with positive diagonal).

### Non-Gaussian observation models

When the observation distribution is non-Gaussian, the linear-Gaussian
Kalman machinery no longer applies, and inference requires
particle-based methods (PGAS, SMC\\^2\\) or specialized approximations
(KSC). The SSM module provides a registry of observation families
accessible via the `obs.model` argument to
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md).

For **Poisson observations**, which arise naturally in count data
settings such as disease surveillance or ecological monitoring, the
observation model is

\\ y_t \mid x_t \sim \text{Poisson}\bigl(\exp(x_t)\bigr), \\

where \\x_t\\ represents the log-intensity and the exponential link
function ensures positivity of the rate parameter. The log-likelihood
contribution is \\\log p(y_t \mid x_t) = y_t x_t - \exp(x_t) -
\log(y_t!)\\.

For **negative binomial observations**, which accommodate overdispersion
relative to the Poisson distribution, the model is

\\ y_t \mid x_t \sim \text{NegBin}\bigl(r, \exp(x_t) / (r +
\exp(x_t))\bigr), \\

where \\r \> 0\\ is the overdispersion parameter. As \\r \to \infty\\,
the negative binomial converges to the Poisson; finite \\r\\ values
accommodate variance that exceeds the mean, a common feature in
ecological and epidemiological count data.

For **Student-\\t\\ observations**, which provide robustness against
outliers in continuous data, the observation model is

\\ y_t \mid x_t \sim t\_\nu(x_t, \sigma^2), \\

where \\\nu \> 0\\ is the degrees of freedom and \\\sigma^2\\ is the
scale parameter. As \\\nu \to \infty\\, this converges to the Gaussian
observation model; small \\\nu\\ values produce heavy tails that
downweight observations far from the predicted state.

For **zero-inflated Poisson observations**, which arise when the data
exhibit excess zeros relative to a standard Poisson model, the
observation is a mixture:

\\ y_t \mid x_t \sim \begin{cases} 0 & \text{with probability } p_0 \\
\text{Poisson}(\exp(x_t)) & \text{with probability } 1 - p_0,
\end{cases} \\

where \\p_0 \in \[0, 1)\\ is the zero-inflation probability. This model
is particularly useful in ecology for modelling species abundance data
where structural zeros (species absent from a site) coexist with
sampling zeros (species present but not detected).

### Process noise structures

The structure of the process noise covariance matrix \\Q\\ has
substantial implications for both model expressiveness and computational
efficiency.

**Isotropic noise** is the simplest case, where \\Q = \sigma^2
I\_{d_x}\\, with a single scalar parameter controlling the magnitude of
all state perturbations equally. This is appropriate when there is no
prior reason to expect differential variability across state components
and is the default in many of the built-in SSM builders.

**Diagonal noise** generalizes to \\Q = \text{diag}(\sigma_1^2, \ldots,
\sigma\_{d_x}^2)\\, allowing each state component to have its own noise
variance while maintaining independence. The BSM builder in lucifer uses
this parameterization, with separate variances for the trend, slope, and
seasonal components.

**Full covariance noise** uses the Cholesky parameterization \\Q = L
L^\top\\, where \\L\\ is a lower-triangular matrix with positive
diagonal entries. This allows correlated state innovations, which can be
important in multivariate models such as the VAR(\\p\\) or dynamic
factor model. The Cholesky parameterization ensures positive
definiteness and provides a one-to-one mapping between unconstrained
parameters and the space of positive-definite matrices.

**State-dependent noise** arises when the diffusion coefficient depends
on the current state, giving \\Q(x_t) = G(x_t) G(x_t)^\top\\. This is
common in SDE-based models; for example, the geometric Brownian motion
has \\\sigma(x) = \sigma x\\, making the noise variance proportional to
the squared state value. State-dependent noise precludes the use of the
standard Kalman filter and requires either particle-based methods or the
UKF/EnKF approximations.

## SSM builder registry

The
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
constructor provides a direct route to discrete-time state-space model
specification without requiring the continuous-time SDE interface. It
accepts a `type` argument that dispatches to a built-in builder from the
SSM builder registry, producing a fully specified model with default
system matrices, parameter names, state names, initial values, and,
where conjugacy holds, a default `theta_update_fn` for Gibbs sampling.
The resulting `ssm_model` object can be passed directly to
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
for inference.

The registry currently provides five model types. The `"local_level"`
type implements the simplest structural time series model, with a single
random-walk state observed with Gaussian noise. The
`"local_linear_trend"` type augments this with a stochastic slope
component, so that both level and trend are estimated. The `"seasonal"`
type implements a dummy-variable-rotation seasonal model with
user-specified period via `period = s`, producing \\s - 1\\ seasonal
states. The `"bsm"` (basic structural model) combines local linear trend
and seasonal components into a single model following Harvey (1989)
[\[8\]](#ref8), with state dimension \\2 + s - 1\\. Finally, `"var_p"`
casts a vector autoregressive model of order \\p\\ as a state-space
model, accepting a `lags` argument.

All built-in types produce conjugate inverse-gamma full conditionals for
the variance parameters, so the default `theta_update_fn` performs exact
Gibbs updates without requiring user intervention. When a user supplies
their own `build_ssm` function or `theta_update_fn` as arguments to
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md),
these override the builder defaults. The following table summarizes the
five built-in types, their state dimensions, parameter counts, and
default priors.

| Type | State dim. | Parameters | Conjugate updates | Extra arguments |
|:---|:---|:---|:---|:---|
| local_level | 1 | sigma_w, sigma_v | Yes (inv-gamma) | none |
| local_linear_trend | 2 | sigma_level, sigma_slope, sigma_v | Yes (inv-gamma) | none |
| seasonal | s - 1 | sigma_seasonal, sigma_v | Yes (inv-gamma) | period |
| bsm | s + 1 | sigma_level, sigma_slope, sigma_seasonal, sigma_v | Yes (inv-gamma) | period |
| var_p | d_y \* p | A_1, …, A_p, Sigma | Yes (normal-inv-Wishart) | lags |

Built-in SSM builder types. {.table}

``` r

# Simulate a local level series
set.seed(42)
y_ll <- cumsum(rnorm(200, sd = 0.5)) + rnorm(200, sd = 1)

# Build via the registry (conjugate theta_update_fn included)
mod_ll <- SSM_model(type = "local_level", data = y_ll)
print(mod_ll)

# Fit with FFBS (auto-selected for ssm_model with build_ssm)
fit_ll <- SSM(mod_ll, Iterations = 2000, verbose = TRUE, Status = 500)
summary(fit_ll)
plot(fit_ll, type = "states")

# Basic structural model with quarterly seasonality
y_bsm <- cumsum(rnorm(200, sd = 0.3)) +
    sin(2 * pi * (1:200) / 4) * 2 + rnorm(200, sd = 0.5)
mod_bsm <- SSM_model(type = "bsm", data = y_bsm, period = 4)
print(mod_bsm)

fit_bsm <- SSM(mod_bsm, Iterations = 3000, verbose = TRUE, Status = 500)
plot(fit_bsm, type = "states")
plot(fit_bsm, type = "parameters")
```

The
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
constructor also accepts a user-supplied `build_ssm` function for custom
discrete-time specifications that do not fall into any of the built-in
categories. The function must take a parameter vector `parm` and return
a list with components `F`, `Q`, `H`, `R`, `m0`, `P0`. Matrices may be
2D (time-invariant) or 3D arrays with the third dimension equal to \\T\\
(time-varying). This interface provides complete flexibility while
maintaining compatibility with the FFBS algorithm.

``` r

# Custom time-varying local level model where observation noise
# increases linearly with time (heteroscedastic observations)
T_len <- 150
set.seed(123)
x_true <- cumsum(rnorm(T_len, sd = 0.3))
sigma_v <- 0.5 + 0.01 * (1:T_len)  # increasing observation noise
y_het <- x_true + rnorm(T_len, sd = sigma_v)

build_het <- function(parm) {
    T_obs <- 150
    R_arr <- array(0, c(1, 1, T_obs))
    for (t in seq_len(T_obs))
        R_arr[1, 1, t] <- (parm[2] + parm[3] * t)^2
    list(
        F = array(1, c(1, 1, 1)),
        Q = array(parm[1]^2, c(1, 1, 1)),
        H = array(1, c(1, 1, 1)),
        R = R_arr,
        m0 = 0,
        P0 = matrix(10, 1, 1)
    )
}

theta_fn <- function(parm, states, Data) {
    proposal <- abs(parm + rnorm(3, 0, 0.02))
    proposal  # MH acceptance handled externally
}

mod_custom <- SSM_model(
    data = y_het,
    build_ssm = build_het,
    theta_update_fn = theta_fn,
    parm.names = c("sigma_w", "sigma_v_base", "sigma_v_slope"),
    state.names = "level",
    Initial.Values = c(0.3, 0.5, 0.01)
)

fit_custom <- SSM(mod_custom, Iterations = 2000, verbose = TRUE, Status = 500)
summary(fit_custom)
plot(fit_custom, type = "states")
```

## Algorithms

### Forward-Filtering Backward-Sampling (FFBS)

FFBS is the gold-standard algorithm for linear-Gaussian state-space
models, providing exact Gibbs samples from the joint posterior
\\p(\theta, x\_{0:T} \mid y\_{1:T})\\ without any approximation error.
The algorithm was developed independently by Carter and Kohn
[\[1\]](#ref1) and Fruhwirth-Schnatter [\[2\]](#ref2), and it alternates
between two steps: sampling the latent states given the parameters via
the Kalman filter and backward simulation, and sampling the parameters
given the latent states via user-specified full conditionals.

The forward pass is a standard Kalman filter that processes observations
sequentially, computing filtered means \\m_t = E\[x_t \mid y\_{1:t},
\theta\]\\ and filtered covariances \\P_t = \text{Cov}(x_t \mid
y\_{1:t}, \theta)\\ at each time step. The backward pass then draws a
complete state trajectory \\x\_{0:T}\\ by sampling backwards from \\t =
T\\ to \\t = 0\\, conditioning on the already-drawn future states:

\\ x_t \mid x\_{t+1}, y\_{1:T}, \theta \sim \mathcal{N}\\\bigl(m_t^s,\\
P_t^s\bigr) \\

where \\m_t^s = m_t + G_t(x\_{t+1} - F\_{t+1} m_t)\\ and \\P_t^s = P_t -
G_t P\_{t+1\|t} G_t^\top\\, with smoother gain \\G_t = P_t F\_{t+1}^\top
P\_{t+1\|t}^{-1}\\. This produces an exact draw from \\p(x\_{0:T} \mid
y\_{1:T}, \theta)\\ in \\O(T \cdot d_x^3)\\ operations, where \\d_x\\ is
the state dimension.

The exactness of FFBS rests on two properties of the linear-Gaussian
model: the Markov property of the state process ensures that the
backward factorization \\p(x\_{0:T} \mid y\_{1:T}, \theta) = p(x_T \mid
y\_{1:T}, \theta) \prod\_{t=0}^{T-1} p(x_t \mid x\_{t+1}, y\_{1:t},
\theta)\\ holds, and the Gaussian linearity ensures that each
conditional \\p(x_t \mid x\_{t+1}, y\_{1:t}, \theta)\\ is itself
Gaussian with closed-form mean and covariance. The backward sampling
step is therefore exact, not an approximation, and the resulting Gibbs
sampler converges to the true posterior at the standard Markov chain
rate.

Convergence of FFBS is governed by the mixing properties of the
parameter update step. When conjugate full conditionals are available,
as in the built-in SSM builders for local level, local linear trend, and
BSM models, the resulting Gibbs sampler is geometrically ergodic and
typically mixes well within a few hundred iterations. The effective
sample size (ESS) of the parameter chain can be monitored via
`as.demonoid(fit)` followed by
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md),
and slow mixing usually indicates either strong posterior correlation
between parameters (which can be addressed by reparameterization) or
weak identifiability.

The lucifer implementation uses a C++ Kalman filter module
(`kalman.cpp`) with Joseph-form covariance updates \\P = (I - KH)P(I -
KH)^\top + KRK^\top\\ for numerical stability, handling of missing
observations (skipping the update step when \\y_t\\ contains NA), and
time-varying system matrices via `arma::cube` slicing. The user must
provide either a `build_ssm` function that constructs the system
matrices \\(F, Q, H, R, m_0, P_0)\\ from the current parameter vector,
or use a recognized SDE family (currently OU and Vasicek) for which a
default builder is available. The `theta_update_fn` draws new parameters
from their full conditional distribution given the sampled states.

FFBS should be preferred whenever its assumptions hold. No other
algorithm in the SSM module matches its efficiency for linear-Gaussian
models: it requires no tuning parameters (no particle count, no proposal
distribution, no tempering schedule), it provides exact samples (no
approximation error), and its cost is \\O(T d_x^3)\\ per iteration,
which is optimal for dense covariance models. The UKF and EnKF provide
cheaper-per-step approximations for nonlinear models but introduce
systematic bias; PGAS provides exact inference for nonlinear models but
requires a particle count large enough to avoid path degeneracy.

#### Algorithm

1.  **Initialize** \\\theta_0\\; run Kalman filter to obtain filtered
    state estimates
2.  **For** \\t = 1, \ldots, T\_{\text{iter}}\\:
    1.  Construct system matrices \\(F_t, Q_t, H_t, R_t, m_0, P_0)\\
        from current \\\theta\\
    2.  **Forward filter**: run Kalman filter to compute \\(m_k, P_k)\\
        for \\k = 1, \ldots, T\\
    3.  **Backward simulate**: draw \\x\_{0:T}^{(t)} \sim p(x\_{0:T}
        \mid y\_{1:T}, \theta^{(t-1)})\\
    4.  **Parameter update**: draw \\\theta^{(t)} \sim p(\theta \mid
        x\_{0:T}^{(t)}, y\_{1:T})\\ via `theta_update_fn`
3.  **Return** \\\\(\theta^{(t)}, x\_{0:T}^{(t)})\\\\

#### Usage

``` r



# Simulate a local level model: y_t = x_t + v_t, x_t = x_{t-1} + w_t
set.seed(42)
T_len <- 100
x_true <- cumsum(rnorm(T_len, 0, sqrt(0.5)))
y <- x_true + rnorm(T_len, 0, 1.0)

sde <- SDE(family = "ou", data = y, times = 1:T_len)

# Custom SSM builder for local level
build_ll <- function(parm) {
    list(F = array(1, c(1, 1, 1)),
         Q = array(parm[1]^2, c(1, 1, 1)),
         H = array(1, c(1, 1, 1)),
         R = array(parm[2]^2, c(1, 1, 1)),
         m0 = 0, P0 = matrix(10, 1, 1))
}

# Conjugate Gibbs update for noise variances
theta_fn <- function(parm, states, Data) {
    y_dat <- as.vector(Data$.y)
    x_state <- states[, 1]
    resid <- y_dat - x_state
    sigma_v_sq <- 1 / rgamma(1, length(resid)/2, sum(resid^2)/2)
    dx <- diff(x_state)
    sigma_w_sq <- 1 / rgamma(1, length(dx)/2, sum(dx^2)/2)
    c(sqrt(sigma_w_sq), sqrt(sigma_v_sq))
}

sde$parm.names <- c("sigma_w", "sigma_v")
sde$Initial.Values <- c(1, 1)

fit <- SSM(sde, Algorithm = "FFBS", Iterations = 2000,
           Specs = list(theta_update_fn = theta_fn,
                        build_ssm = build_ll),
           verbose = TRUE, Status = 500)

print(fit)
summary(fit)
plot(fit, type = "states")
plot(fit, type = "parameters")
```

#### Extended example: local linear trend with FFBS

The following example demonstrates FFBS on a local linear trend model,
showing how the Gibbs sampler jointly estimates the level, slope, and
all noise variances using conjugate updates.

``` r



# Simulate local linear trend: level + slope + observation noise
set.seed(99)
T_len <- 200
sigma_level <- 0.3
sigma_slope <- 0.05
sigma_obs   <- 1.0

# True states
mu_true  <- numeric(T_len)
beta_true <- numeric(T_len)
mu_true[1] <- 0; beta_true[1] <- 0.1
for (t in 2:T_len) {
    beta_true[t] <- beta_true[t-1] + rnorm(1, sd = sigma_slope)
    mu_true[t]   <- mu_true[t-1] + beta_true[t-1] + rnorm(1, sd = sigma_level)
}
y_trend <- mu_true + rnorm(T_len, sd = sigma_obs)

# Build via SSM_model registry
mod_trend <- SSM_model(type = "local_linear_trend", data = y_trend)

# Fit (conjugate theta_update_fn supplied by builder)
fit_trend <- SSM(mod_trend, Iterations = 3000, verbose = TRUE, Status = 500)

# Inspect
print(fit_trend)
summary(fit_trend)
plot(fit_trend, type = "states")       # level and slope trajectories
plot(fit_trend, type = "parameters")   # variance posteriors
```

### Particle Gibbs with Ancestor Sampling (PGAS)

PGAS extends the Gibbs sampling framework to nonlinear and non-Gaussian
state-space models by replacing the Kalman filter with a conditional
particle filter. The algorithm was introduced by Andrieu, Doucet, and
Holenstein [\[3\]](#ref3) in the context of Particle MCMC, and the
critical ancestor sampling refinement was added by Lindsten, Jordan, and
Schon [\[4\]](#ref4) to break the path degeneracy that plagues standard
particle Gibbs.

In standard particle Gibbs, one particle’s trajectory is fixed to the
reference trajectory from the previous MCMC iteration, and the remaining
\\N - 1\\ particles are propagated freely. After weighting and
resampling, a new trajectory is selected from the particle system. The
problem is that the conditioned trajectory “locks in” early time points:
since the reference particle at time \\t = 1\\ is always preserved
through resampling, the sampled trajectory at early times is nearly
deterministic, requiring exponentially many particles to mix. Ancestor
sampling solves this by resampling the ancestor of the reference
particle at each time step proportional to the backward transition
density \\w_k(i) \cdot p(x\_{k+1}^{\text{ref}} \mid x_k^{(i)},
\theta)\\, allowing the trajectory to detach from its frozen past and
explore the state space freely.

The mathematical justification for ancestor sampling rests on the
extended target distribution framework. Define the extended target as

\\ \tilde{\pi}(x\_{1:T}^{1:N}, a\_{1:T-1}^{1:N}) = \frac{1}{N^T}
\prod\_{t=1}^T \prod\_{i=1}^N q(x_t^i \mid x\_{t-1}^{a\_{t-1}^i},
\theta) \cdot \prod\_{t=1}^{T-1} \prod\_{i=1}^N w_t(x_t^{a\_{t}^i}), \\

where \\a_t^i\\ denotes the ancestor index of particle \\i\\ at time
\\t\\ and \\q\\ is the proposal distribution (taken as the transition
density in the bootstrap filter). Andrieu et al. [\[3\]](#ref3) showed
that the standard particle Gibbs algorithm, which conditions on one
particle trajectory and resamples the rest, leaves this extended target
invariant, ensuring that the marginal draws of \\x\_{1:T}\\ converge to
the correct smoothing distribution. Lindsten et al. [\[4\]](#ref4)
showed that adding ancestor sampling preserves this invariance while
dramatically improving mixing by allowing the reference trajectory to
“re-root” at each time step.

The practical effect of ancestor sampling is that the number of
particles \\N\\ needed for adequate mixing is typically modest, on the
order of 50 to 200 for state dimensions \\d_x \leq 5\\. Without ancestor
sampling, the required particle count grows exponentially with \\T\\;
with it, the computational cost is \\O(T \cdot N \cdot d_x)\\ per MCMC
iteration, independent of \\T\\ in terms of mixing. The quality of the
ancestor weights depends on the accuracy of the transition density
evaluation. For Euler-Maruyama discretizations of SDEs, finer
discretization (smaller `dt`) improves the ancestor weights and hence
the mixing of the sampler.

The lucifer implementation uses a C++ conditional particle filter
(`pgas.cpp`) that stores the full particle genealogy for trajectory
extraction. The Euler-Maruyama transition density is used for computing
ancestor sampling weights. As with FFBS, the user provides a
`theta_update_fn` for parameter updates. PGAS is the default algorithm
for general nonlinear SDE models where no conjugate structure is
available.

#### Algorithm

1.  **Initialize** \\\theta_0\\; run a bootstrap particle filter to
    obtain an initial reference trajectory \\x\_{0:T}^{\text{ref}}\\
2.  **For** \\t = 1, \ldots, T\_{\text{iter}}\\:
    1.  **Conditional PF with ancestor sampling**: run \\N\\ particles,
        fixing particle \\N\\ to \\x\_{0:T}^{\text{ref}}\\. At each
        resampling step, resample the ancestor of particle \\N\\ from
        weights \\\propto w_k(i) \cdot p(x\_{k+1}^{\text{ref}} \mid
        x_k^{(i)}, \theta)\\
    2.  Sample a new trajectory from the particle system:
        \\x\_{0:T}^{(t)} \sim \hat{p}(x\_{0:T} \mid y\_{1:T},
        \theta^{(t-1)})\\
    3.  **Parameter update**: \\\theta^{(t)} \sim p(\theta \mid
        x\_{0:T}^{(t)}, y\_{1:T})\\ via `theta_update_fn`
3.  **Return** \\\\(\theta^{(t)}, x\_{0:T}^{(t)})\\\\

#### Usage

``` r

# Simulate OU process
set.seed(42)
T_len <- 100
kappa_true <- 2; mu_true <- 5; sigma_true <- 1; obs_sd_true <- 0.5
dt_sim <- rep(0.1, T_len)
x_true <- numeric(T_len); x_true[1] <- mu_true
for (t in 2:T_len)
    x_true[t] <- x_true[t-1] + kappa_true * (mu_true - x_true[t-1]) * dt_sim[t] +
        sigma_true * sqrt(dt_sim[t]) * rnorm(1)
y <- x_true + rnorm(T_len, 0, obs_sd_true)

sde <- SDE(family = "ou", data = y, times = seq(0.1, 10, length.out = T_len))

# Random-walk parameter update (positivity enforced)
theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.02)
    if (any(proposal[c(1, 3, 4)] <= 0)) return(parm)
    proposal
}

fit <- SSM(sde, Algorithm = "PGAS", Iterations = 500,
           Specs = list(N.particles = 30,
                        theta_update_fn = theta_fn),
           verbose = TRUE, Status = 250)

print(fit)
plot(fit, type = "states")
```

#### Extended example: nonlinear population dynamics with PGAS

The following example demonstrates PGAS on a nonlinear ecological model
where the latent population follows Ricker dynamics observed with
Poisson counting noise. This is a case where neither FFBS (nonlinear
dynamics) nor KSC (non-SV observation model) apply, making PGAS the
natural choice.

``` r

# Simulate Ricker-Poisson SSM
# x_t = log(N_t), latent log-population
# x_t = x_{t-1} + r * (1 - exp(x_{t-1}) / K) + w_t  (Ricker on log-scale)
# y_t ~ Poisson(exp(x_t))
set.seed(123)
T_len <- 150
r_true <- 0.5        # intrinsic growth rate
K_true <- 100        # carrying capacity
sigma_w <- 0.1       # process noise SD

x_true <- numeric(T_len)
x_true[1] <- log(50)
for (t in 2:T_len) {
    x_true[t] <- x_true[t-1] + r_true * (1 - exp(x_true[t-1]) / K_true) +
        rnorm(1, sd = sigma_w)
}
y_counts <- rpois(T_len, lambda = exp(x_true))

# Build SDE model with Poisson observations
sde_ricker <- SDE(
    drift = function(x, theta, t) theta[1] * (1 - exp(x[1]) / theta[2]),
    diffusion = function(x, theta, t) theta[3],
    data = y_counts,
    times = 1:T_len,
    x0 = log(50),
    obs.model = "poisson",
    obs.link = function(x, theta) exp(x[1]),
    prior = function(theta) {
        if (any(theta <= 0)) return(-Inf)
        sum(dnorm(log(theta), log(c(0.5, 100, 0.1)), 1, log = TRUE))
    },
    parm.names = c("r", "K", "sigma_w"),
    Initial.Values = c(0.3, 80, 0.2),
    method = "particle", N.particles = 100
)

# MH parameter update: log-normal random walk preserves positivity
theta_fn <- function(parm, states, Data) {
    parm * exp(rnorm(length(parm), 0, 0.05))
}

fit_ricker <- SSM(sde_ricker, Algorithm = "PGAS", Iterations = 2000,
                  Specs = list(N.particles = 100,
                               theta_update_fn = theta_fn),
                  verbose = TRUE, Status = 500)

print(fit_ricker)
summary(fit_ricker)
plot(fit_ricker, type = "states")
plot(fit_ricker, type = "parameters")
```

### SMC-squared (SMC\\^2\\)

SMC\\^2\\ is the most general algorithm for state-space models with
unknown static parameters, requiring neither full conditional
distributions nor tractable transition densities. It was introduced by
Chopin, Jacob, and Papaspiliopoulos [\[5\]](#ref5) and nests two levels
of sequential Monte Carlo: an outer SMC sampler that operates on the
parameter space, and an inner particle filter that estimates the
likelihood \\p(y\_{1:t} \mid \theta)\\ for each parameter particle. The
outer level incorporates observations sequentially; each time a new
observation arrives, every outer particle’s inner particle filter is
advanced by one step, producing an incremental weight. When the outer
ESS drops, the parameter particles are resampled and rejuvenated via
MCMC, with the inner particle filters re-initialized for accepted
proposals.

The mathematical framework proceeds as follows. At stage \\t\\, each
outer particle \\\theta^{(i)}\\ carries an inner particle filter
\\\\x\_{1:t}^{(i,j)}, w_t^{(i,j)}\\\_{j=1}^{N\_{\text{inner}}}\\ that
provides an unbiased estimate of the marginal likelihood up to time
\\t\\:

\\ \hat{p}(y\_{1:t} \mid \theta^{(i)}) = \prod\_{k=1}^t \left(
\frac{1}{N\_{\text{inner}}} \sum\_{j=1}^{N\_{\text{inner}}} w_k^{(i,j)}
\right). \\

When observation \\y\_{t+1}\\ arrives, each inner particle filter
advances one step, yielding an incremental weight \\\Delta \log w^{(i)}
= \log \hat{p}(y\_{t+1} \mid y\_{1:t}, \theta^{(i)})\\ that updates the
outer particle’s importance weight. The outer ESS is computed as
\\\text{ESS} = (\sum_i W_i^2)^{-1}\\, where \\W_i\\ are the normalized
outer weights. When ESS falls below the threshold (default: \\0.5 \times
N\_{\text{outer}}\\), the outer particles are resampled and rejuvenated:
for each particle, a new \\\theta^\*\\ is proposed via random-walk
Metropolis, a fresh particle filter is run to obtain \\\hat{p}(y\_{1:t}
\mid \theta^\*)\\, and the proposal is accepted or rejected via the
standard MH ratio.

The key theoretical guarantee of SMC\\^2\\ is that the normalizing
constant estimate \\\hat{p}(y\_{1:T})\\ is unbiased regardless of the
number of inner particles, although more inner particles reduce
variance. This property makes SMC\\^2\\ suitable for model comparison
via Bayes factors: the ratio \\\hat{p}(y\_{1:T} \mid \mathcal{M}\_1) /
\hat{p}(y\_{1:T} \mid \mathcal{M}\_2)\\ provides a consistent estimator
of the Bayes factor between two competing SSM specifications.

The computational cost of SMC\\^2\\ is \\O(T \cdot N\_{\text{outer}}
\cdot N\_{\text{inner}} \cdot d_x)\\ per time step, which can be
substantial. The advantage is complete generality: no conjugacy, no
gradient, no user-specified samplers. The marginal likelihood estimate
is a byproduct, enabling model comparison between competing SSM
specifications. Typical settings for moderate-sized problems use
\\N\_{\text{outer}} = 100\\-\\500\\ and \\N\_{\text{inner}} =
20\\-\\50\\.

The lucifer implementation uses a dedicated C++ backend (`ssm_smc2_cpp`
in `src/ssm_engines.cpp`) that runs the entire triple-nested loop (outer
particles, time steps, inner particles) natively. The inner particle
filters use Euler-Maruyama propagation with R callbacks for user-defined
drift and diffusion functions; systematic resampling and weight
normalization are handled entirely in C++. The outer rejuvenation step
runs a quick C++ particle filter for each proposed parameter value and
accepts or rejects via the standard MH ratio on the approximate log
marginal likelihood.

SMC\\^2\\ should be chosen when the modeller has neither conjugate full
conditionals for \\\theta\\ nor the ability to write a
`theta_update_fn`, or when the transition density is intractable (e.g.,
for jump-diffusion processes or agent-based models). For models where a
`theta_update_fn` is available, PGAS is generally more efficient because
it avoids the quadratic cost of maintaining \\N\_{\text{outer}}\\
independent particle filters.

#### Algorithm

1.  **Initialize** \\N\_{\text{outer}}\\ parameter particles
    \\\theta^{(i)}\\, each with an inner PF of \\N\_{\text{inner}}\\
    state particles
2.  **For** each observation \\y_k\\, \\k = 1, \ldots, T\\:
    1.  **Propagate** each inner PF forward one time step given its
        \\\theta^{(i)}\\
    2.  **Reweight** outer particles: \\\Delta \log w_i = \log
        \hat{p}(y_k \mid y\_{1:k-1}, \theta^{(i)})\\ (from inner PF)
    3.  If outer ESS \\\<\\ threshold:
        - **Resample** outer particles (systematic resampling)
        - **Rejuvenate**: for each outer particle, propose \\\theta^\*\\
          via RWM, run a quick PF to estimate \\\log p(y\_{1:k} \mid
          \theta^\*)\\, accept/reject via MH
        - Re-initialize inner PFs for accepted proposals
3.  **Return** weighted parameter particles and filtered state means

#### Usage

``` r

# Simulate OU data (self-contained)
set.seed(42)
T_len <- 50
x_ou <- 5 + cumsum(rnorm(T_len, 0, 0.3))
y_ou <- x_ou + rnorm(T_len, 0, 0.5)
sde <- SDE(family = "ou", data = y_ou, times = seq(0.1, 5, length.out = T_len))

fit <- SSM(sde, Algorithm = "SMC2",
           Specs = list(N.outer = 100, N.inner = 30,
                        rejuvenation_steps = 3))

print(fit)
plot(fit, type = "states")
```

#### Extended example: unknown transition model with SMC\\^2\\

This example shows SMC\\^2\\ on a model where the drift function has a
complex parameterization and no obvious conjugate structure,
illustrating the “plug and play” nature of the algorithm. The user only
needs to define the SDE model; no `theta_update_fn` is required.

``` r



# Simulate a mean-reverting process with time-varying speed
# dx = kappa(t) * (mu - x) dt + sigma dW
# kappa(t) = a + b * sin(2*pi*t/P)
set.seed(42)
T_len <- 200
dt <- 0.1
times <- seq(0, by = dt, length.out = T_len)
a_true <- 0.5; b_true <- 0.3; P_true <- 50
mu_true <- 2.0; sigma_true <- 0.4

x <- numeric(T_len); x[1] <- mu_true
for (t in 2:T_len) {
    kappa_t <- a_true + b_true * sin(2 * pi * times[t] / P_true)
    x[t] <- x[t-1] + kappa_t * (mu_true - x[t-1]) * dt +
        sigma_true * rnorm(1, sd = sqrt(dt))
}
y_obs <- x + rnorm(T_len, sd = 0.3)

# SDE model: no theta_update_fn needed for SMC2
sde_tv <- SDE(
    data = y_obs,
    times = times,
    x0 = y_obs[1],
    drift = function(x, theta, t) {
        kappa <- theta[1] + theta[2] * sin(2 * pi * t / theta[3])
        kappa * (theta[4] - x[1])
    },
    diffusion = function(x, theta, t) theta[5],
    obs.model = "gaussian",
    obs.noise = function(theta) 0.3,
    prior = function(theta) {
        if (any(theta[c(1, 5)] <= 0)) return(-Inf)
        sum(dnorm(theta, c(0.5, 0.2, 50, 1, 0.5), 2, log = TRUE))
    },
    parm.names = c("a", "b", "period", "mu", "sigma"),
    Initial.Values = c(0.3, 0.2, 40, 1.5, 0.3),
    method = "particle", N.particles = 50
)

# SMC2: fully automatic, no user-specified updates
fit_smc2 <- SSM(sde_tv, Algorithm = "SMC2",
                Specs = list(N.outer = 200, N.inner = 50,
                             rejuvenation_steps = 5))

print(fit_smc2)
summary(fit_smc2)
plot(fit_smc2, type = "states")
plot(fit_smc2, type = "parameters")
```

### Kim-Shephard-Chib (KSC)

The Kim-Shephard-Chib algorithm [\[6\]](#ref6) is a specialized Gibbs
sampler for the standard log-stochastic-volatility model:

\\ y_t = \exp(h_t / 2) \cdot \varepsilon_t, \qquad h\_{t+1} = \mu +
\phi(h_t - \mu) + \sigma\_\eta \cdot \eta_t, \\

where \\\varepsilon_t \sim \mathcal{N}(0, 1)\\ and \\\eta_t \sim
\mathcal{N}(0, 1)\\. Taking the log-squared transformation \\y_t^\* =
\log(y_t^2 + c)\\ (with a small offset \\c\\ to handle zeros) yields
\\y_t^\* = h_t + \log(\varepsilon_t^2)\\. The distribution
\\\log(\chi^2_1)\\ is non-Gaussian, but Kim, Shephard, and Chib showed
that it can be approximated accurately by a 7-component Gaussian mixture
with known weights, means, and variances. Conditional on the mixture
indicators \\s_t \in \\1, \ldots, 7\\\\ that select which component
generated each \\\log(\varepsilon_t^2)\\, the model becomes a
linear-Gaussian SSM with time-varying observation noise, and the full
state trajectory \\h\_{1:T}\\ can be drawn in one block via FFBS.

The 7-component mixture approximation is the key insight that makes KSC
computationally efficient. The following table reproduces the mixture
parameters from Table 4 of Kim et al. [\[6\]](#ref6):

| Component | Weight (p_k) | Mean (m_k) | Variance (v_k^2) |
|----------:|-------------:|-----------:|-----------------:|
|         1 |      0.00730 |  -10.12999 |          5.79596 |
|         2 |      0.10556 |   -3.97281 |          2.61369 |
|         3 |      0.00002 |   -8.56686 |          5.17950 |
|         4 |      0.04395 |    2.77786 |          0.16735 |
|         5 |      0.34001 |    0.61942 |          0.64009 |
|         6 |      0.24566 |    1.79518 |          0.34023 |
|         7 |      0.25750 |   -1.08819 |          1.26261 |

7-component Gaussian mixture approximation to log(chi-squared_1).
{.table}

The approximation is accurate to within 4 decimal places of the exact
\\\log(\chi^2_1)\\ density, as demonstrated by Kolmogorov-Smirnov tests
in the original paper [\[6\]](#ref6). The use of a mixture approximation
introduces a small bias relative to the exact \\\log(\chi^2_1)\\
density, but Kim et al. showed that this bias is negligible in practice
and does not affect posterior inference for the parameters \\(\mu, \phi,
\sigma\_\eta)\\.

The Gibbs cycle alternates five steps: (1) sample the mixture indicators
\\s_t\\ from their categorical conditional given \\(y_t^\*, h_t)\\; (2)
sample \\h\_{1:T}\\ via FFBS using the Kalman filter with \\R_t =
\sigma^2\_{s_t}\\; (3) sample \\\mu\\ from its conjugate Gaussian
posterior; (4) sample \\\phi\\ from a truncated normal (constraining
\\\|\phi\| \< 1\\ for stationarity); and (5) sample \\\sigma\_\eta^2\\
from its conjugate inverse-gamma posterior. The posterior for \\\mu\\ is
\\\mathcal{N}(\hat{\mu}, V\_\mu)\\ where \\V\_\mu = (1/V_0 + (1 -
\phi)^2 (T-1) / \sigma\_\eta^2)^{-1}\\ and \\\hat{\mu} = V\_\mu
(m_0/V_0 + (1-\phi)/\sigma\_\eta^2 \sum\_{t=2}^T (h_t - \phi
h\_{t-1}))\\. The posterior for \\\sigma\_\eta^2\\ is
\\\text{InvGamma}(\alpha_0 + (T-1)/2, \beta_0 + \sum\_{t=2}^T (h_t -
\mu - \phi(h\_{t-1} - \mu))^2 / 2)\\. The stationarity constraint on
\\\phi\\ is enforced by rejection sampling from a truncated normal
proposal.

KSC is activated automatically when
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
detects the `"sv"` SDE family, and it requires no user-supplied
`theta_update_fn` because all full conditionals are built in. The entire
Gibbs loop runs in a dedicated C++ backend (`ssm_ksc_cpp` in
`src/ssm_engines.cpp`): mixture indicator sampling, Kalman
filter/backward simulation, and conjugate parameter updates all execute
natively without crossing the R/C++ boundary during the inner loop. The
7-component mixture constants are compiled as `constexpr` arrays in C++.

The primary diagnostic for KSC convergence is the autocorrelation of the
\\\phi\\ parameter, which often exhibits slow mixing because the
persistence parameter is tightly coupled with the volatility trajectory.
Running chains of 5000 or more iterations and thinning by 5 to 10 is
typical practice. The [`summary()`](https://rdrr.io/r/base/summary.html)
method reports ESS for each parameter, and the
`plot(fit, type = "parameters")` method displays trace plots that reveal
any non-stationarity.

#### Algorithm

1.  **Initialize** \\(\mu, \phi, \sigma\_\eta)\\, \\h\_{1:T}\\; compute
    \\y_t^\* = \log(y_t^2 + 10^{-4})\\
2.  **For** \\t = 1, \ldots, T\_{\text{iter}}\\:
    1.  **Sample indicators** \\s_t \sim \text{Categorical}\bigl(p_k
        \cdot \mathcal{N}(y_t^\*; h_t + m_k, v_k)\bigr)\\ for \\k = 1,
        \ldots, 7\\
    2.  **Sample states** \\h\_{1:T}\\ via FFBS with \\F = \phi\\, \\Q =
        \sigma\_\eta^2\\, \\H = 1\\, \\R_t = v\_{s_t}\\, offset \\=
        m\_{s_t}\\
    3.  **Sample** \\\mu \sim \mathcal{N}(m\_\mu, V\_\mu)\\ from
        conjugate posterior
    4.  **Sample** \\\phi\\ via MH with truncated normal proposal
        (\\\|\phi\| \< 1\\)
    5.  **Sample** \\\sigma\_\eta^2 \sim \text{InvGamma}(\alpha_n,
        \beta_n)\\ from conjugate posterior
3.  **Return** \\\\(\mu, \phi, \sigma\_\eta, h\_{1:T})^{(t)}\\\\

#### Usage

``` r



# Simulate stochastic volatility model
set.seed(42)
T_len <- 500
mu <- -0.5; phi <- 0.97; sigma_eta <- 0.15
h <- numeric(T_len); h[1] <- mu
for (t in 2:T_len) h[t] <- mu + phi * (h[t-1] - mu) + sigma_eta * rnorm(1)
y_sv <- exp(h / 2) * rnorm(T_len)

# Build SDE model (auto-detects SV family -> KSC)
sde <- SDE(family = "sv", data = y_sv, times = 1:T_len)

# Run KSC (no theta_update_fn needed)
fit <- SSM(sde, Algorithm = "KSC", Iterations = 3000,
           verbose = TRUE, Status = 500)

print(fit)
summary(fit)
plot(fit, type = "states")       # log-volatility trajectory
plot(fit, type = "parameters")   # mu, phi, sigma_eta posteriors
```

#### Extended example: stochastic volatility with leverage

The standard SV model assumes independence between return shocks
\\\varepsilon_t\\ and volatility shocks \\\eta_t\\. The leverage
extension introduces negative correlation \\\rho \< 0\\ between these
two innovations, capturing the well-documented asymmetric volatility
effect where negative returns are associated with increases in
volatility. While the standard KSC algorithm does not handle leverage
directly, the model can be augmented by sampling \\\rho\\ conditional on
the states and returns, and then adjusting the FFBS step accordingly.

``` r



# Simulate SV with leverage effect
set.seed(55)
T_len <- 600
mu <- -0.3; phi <- 0.98; sigma_eta <- 0.12; rho <- -0.5

h <- numeric(T_len); h[1] <- mu
eps <- rnorm(T_len)
for (t in 2:T_len) {
    eta <- rho * eps[t-1] + sqrt(1 - rho^2) * rnorm(1)
    h[t] <- mu + phi * (h[t-1] - mu) + sigma_eta * eta
}
y_lev <- exp(h / 2) * eps

# Fit standard KSC (ignoring leverage for comparison)
sde_sv <- SDE(family = "sv", data = y_lev, times = 1:T_len)
fit_sv <- SSM(sde_sv, Algorithm = "KSC", Iterations = 5000,
              verbose = TRUE, Status = 1000)

summary(fit_sv)
plot(fit_sv, type = "states")
plot(fit_sv, type = "parameters")
```

### Unscented Kalman Filter (UKF)

The Unscented Kalman Filter, introduced by Julier and Uhlmann
[\[9\]](#ref9) and refined by Wan and van der Merwe [\[10\]](#ref10),
propagates a deterministic set of \\2n + 1\\ sigma points through the
exact nonlinear dynamics, avoiding the finite-difference Jacobian
approximation that the Extended Kalman Filter (EKF) requires. Given a
state of dimension \\n\\ with mean \\\hat{x}\\ and covariance \\P\\, the
UKF generates sigma points \\\chi_0 = \hat{x}\\ and \\\chi_i = \hat{x}
\pm (\sqrt{(n + \lambda)P})\_i\\ for \\i = 1, \ldots, n\\, where
\\\lambda = \alpha^2(n + \kappa) - n\\ controls the spread of the sigma
points around the mean. These points are propagated through the
nonlinear transition function \\f(\cdot)\\, and the predicted mean and
covariance are reconstructed from weighted sums of the propagated
points:

\\ \hat{x}\_{t\|t-1} = \sum\_{i=0}^{2n} w_i^{(m)} f(\chi_i), \qquad
P\_{t\|t-1} = \sum\_{i=0}^{2n} w_i^{(c)} \bigl(f(\chi_i) -
\hat{x}\_{t\|t-1}\bigr)\bigl(f(\chi_i) - \hat{x}\_{t\|t-1}\bigr)^\top +
Q_t, \\

where the weights are

\\ w_0^{(m)} = \frac{\lambda}{n + \lambda}, \qquad w_0^{(c)} =
\frac{\lambda}{n + \lambda} + (1 - \alpha^2 + \beta), \qquad w_i^{(m)} =
w_i^{(c)} = \frac{1}{2(n + \lambda)} \text{ for } i = 1, \ldots, 2n. \\

The parameter \\\beta = 2\\ is optimal for Gaussian distributions, and
\\\alpha \in (10^{-4}, 1\]\\ controls how tightly the sigma points
cluster around the mean. For the observation update, the same
sigma-point transform is applied to the observation function
\\h(\cdot)\\. The cross-covariance between the predicted state and
predicted observation is

\\ P\_{xy} = \sum\_{i=0}^{2n} w_i^{(c)} \bigl(f(\chi_i) -
\hat{x}\_{t\|t-1}\bigr)\bigl(h(f(\chi_i)) -
\hat{y}\_{t\|t-1}\bigr)^\top, \\

and the Kalman gain is \\K_t = P\_{xy} S_t^{-1}\\, where \\S_t = \sum_i
w_i^{(c)} (h(f(\chi_i)) - \hat{y}\_{t\|t-1})(h(f(\chi_i)) -
\hat{y}\_{t\|t-1})^\top + R_t\\ is the innovation covariance. The update
equations are then \\\hat{x}\_t = \hat{x}\_{t\|t-1} + K_t (y_t -
\hat{y}\_{t\|t-1})\\ and \\P_t = P\_{t\|t-1} - K_t S_t K_t^\top\\.

The UKF captures mean and covariance exactly up to third-order Taylor
expansion for Gaussian inputs and second order for non-Gaussian inputs,
compared to the EKF’s first-order accuracy. This makes it substantially
more accurate for models with moderate nonlinearity, such as the Van der
Pol oscillator or double-well potential systems, without requiring
Jacobian computation.

The lucifer implementation uses a C++ UKF filter (`ukf.cpp`) that
propagates sigma points through the user-supplied drift and diffusion
functions with Euler-Maruyama substeps. Within the Gibbs loop, the UKF
replaces the Kalman filter in the forward pass, producing filtered state
estimates that are then passed to `theta_update_fn` for parameter
updates. The UKF tuning parameters \\\alpha\\, \\\beta\\, and \\\kappa\\
are controlled via `Specs$alpha`, `Specs$beta`, and `Specs$kappa`
respectively.

The primary advantage of the UKF over PGAS for nonlinear models is
computational cost: the UKF processes each time step with a fixed number
of deterministic sigma points rather than a stochastic ensemble of
particles, making it \\O(T n^3)\\ rather than \\O(T N n)\\. The tradeoff
is accuracy. The unscented transform captures mean and covariance
exactly up to second-order Taylor expansion for any nonlinearity, but it
assumes Gaussian state distributions at each step. For models with
multimodal or heavily skewed state posteriors, PGAS remains the
appropriate choice. The UKF is most effective when the nonlinearity is
smooth and the state posterior at each time step is approximately
unimodal.

#### Usage

``` r



# Simulate a double-well potential: x_t = x_{t-1} + (x - x^3) * dt + sigma * dW_t
set.seed(42)
T_len <- 200
dt_sim <- 0.1
x_dw <- numeric(T_len); x_dw[1] <- 0.5
for (t in 2:T_len)
    x_dw[t] <- x_dw[t-1] + (x_dw[t-1] - x_dw[t-1]^3) * dt_sim +
        0.3 * rnorm(1, sd = sqrt(dt_sim))
y_dw <- x_dw + rnorm(T_len, 0, 0.5)

# Build SDE model with nonlinear drift
sde_dw <- SDE(
    data = y_dw,
    times = seq(0.1, by = dt_sim, length.out = T_len),
    x0 = 0.5,
    drift = function(x, theta, t) theta[1] * (x[1] - x[1]^3),
    diffusion = function(x, theta, t) theta[2],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[3],
    prior = function(theta) {
        if (any(theta[2:3] <= 0)) return(-Inf)
        sum(dnorm(theta, c(1, 0.3, 0.5), 1, log = TRUE))
    },
    parm.names = c("drift_scale", "sigma_state", "sigma_obs"),
    Initial.Values = c(1.0, 0.3, 0.5),
    method = "ekf"
)

# Gibbs parameter update
theta_fn <- function(parm, states, Data) {
    parm + rnorm(length(parm), 0, 0.02)
}

# Compare EKF-based PGAS vs. UKF
fit_ukf <- SSM(sde_dw, Algorithm = "UKF", Iterations = 2000,
               Specs = list(theta_update_fn = theta_fn,
                            alpha = 1e-3, beta = 2.0, kappa = 0),
               verbose = TRUE, Status = 500)

print(fit_ukf)
plot(fit_ukf, type = "states")
plot(fit_ukf, type = "parameters")
```

#### Extended example: Van der Pol oscillator with UKF

The Van der Pol oscillator provides a more challenging test case for the
UKF because the dynamics are two-dimensional and exhibit limit-cycle
behaviour. The system is

\\ \dot{x}\_1 = x_2, \qquad \dot{x}\_2 = \mu (1 - x_1^2) x_2 - x_1 +
\sigma_w \xi(t), \\

where \\\mu\\ controls the nonlinearity and \\\xi(t)\\ is white noise.
After Euler-Maruyama discretization, this becomes a two-dimensional
nonlinear SSM with observed \\y_t = x\_{1,t} + v_t\\.

``` r



# Simulate Van der Pol oscillator
set.seed(77)
T_len <- 300
dt_sim <- 0.05
mu_vdp <- 1.5
sigma_w <- 0.2
sigma_obs <- 0.5

x1 <- numeric(T_len); x2 <- numeric(T_len)
x1[1] <- 1.0; x2[1] <- 0.0
for (t in 2:T_len) {
    x1[t] <- x1[t-1] + x2[t-1] * dt_sim
    x2[t] <- x2[t-1] + (mu_vdp * (1 - x1[t-1]^2) * x2[t-1] - x1[t-1]) * dt_sim +
        sigma_w * rnorm(1, sd = sqrt(dt_sim))
}
y_vdp <- x1 + rnorm(T_len, sd = sigma_obs)

# Build model (2D state, 1D observation)
sde_vdp <- SDE(
    data = y_vdp,
    times = seq(0.05, by = dt_sim, length.out = T_len),
    x0 = c(1.0, 0.0),
    drift = function(x, theta, t) {
        c(x[2], theta[1] * (1 - x[1]^2) * x[2] - x[1])
    },
    diffusion = function(x, theta, t) {
        matrix(c(0, theta[2]), 2, 1)
    },
    obs.model = "gaussian",
    obs.noise = function(theta) theta[3],
    prior = function(theta) {
        if (any(theta[2:3] <= 0)) return(-Inf)
        sum(dnorm(theta, c(1, 0.2, 0.5), 1, log = TRUE))
    },
    parm.names = c("mu", "sigma_w", "sigma_obs"),
    state.names = c("x1", "x2"),
    Initial.Values = c(1.0, 0.3, 0.4),
    method = "particle", N.particles = 50
)

theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.02)
    if (any(proposal[2:3] <= 0)) return(parm)
    proposal
}

fit_vdp <- SSM(sde_vdp, Algorithm = "UKF", Iterations = 3000,
               Specs = list(theta_update_fn = theta_fn,
                            alpha = 0.01, beta = 2.0, kappa = 0),
               verbose = TRUE, Status = 500)

summary(fit_vdp)
plot(fit_vdp, type = "states")
plot(fit_vdp, type = "parameters")
```

### Ensemble Kalman Filter (EnKF)

For high-dimensional state spaces where maintaining and propagating a
full \\n \times n\\ covariance matrix is intractable, the Ensemble
Kalman Filter, introduced by Evensen [\[11\]](#ref11), represents the
filtering distribution empirically via \\N_e\\ ensemble members. At each
time step, each ensemble member is propagated through the nonlinear
transition function independently, and the ensemble covariance serves as
the forecast covariance in the Kalman update. The stochastic variant
perturbs observations by adding Gaussian noise \\\varepsilon_t^{(j)}
\sim \mathcal{N}(0, R)\\ to each ensemble member’s update, yielding the
analysis step:

\\ x_t^{(j),a} = x_t^{(j),f} + P_t^f H^\top (H P_t^f H^\top + R)^{-1}
\bigl(y_t + \varepsilon_t^{(j)} - H x_t^{(j),f}\bigr), \\

where \\P_t^f = \frac{1}{N_e - 1}\sum\_{j=1}^{N_e} (x_t^{(j),f} -
\bar{x}\_t^f)(x_t^{(j),f} - \bar{x}\_t^f)^\top\\ is the sample
covariance of the forecast ensemble.

The ensemble approach avoids ever forming the full \\n \times n\\
covariance matrix explicitly. Instead, all computations can be expressed
in terms of the \\n \times N_e\\ anomaly matrix \\A_t = \[x_t^{(1),f} -
\bar{x}\_t^f, \ldots, x_t^{(N_e),f} - \bar{x}\_t^f\]\\, so that \\P_t^f
= A_t A_t^\top / (N_e - 1)\\. The matrix inversion in the Kalman gain
involves only the \\d_y \times d_y\\ innovation covariance \\S_t = H
P_t^f H^\top + R\\, which is typically low-dimensional even when the
state dimension is large. This gives the EnKF a computational cost of
\\O(T N_e n)\\ per iteration, which is linear in the state dimension
rather than cubic.

**Covariance inflation** is an essential practical technique for the
EnKF. With finite ensemble sizes, the sample covariance \\P_t^f\\
systematically underestimates the true forecast covariance, leading to
overconfident filter estimates that can diverge from the true state.
Multiplicative inflation replaces \\P_t^f\\ with \\(1 + \delta) P_t^f\\
for a small \\\delta \> 0\\ (controlled by `Specs$inflation`, default
1.0 meaning no inflation). Typical values of \\\delta\\ range from 0.01
to 0.10, depending on the ensemble size and model dynamics.
**Localization**, which zeroes out long-range correlations in \\P_t^f\\,
is another technique used in geophysical applications when the state
dimension is extremely large (\\n \> 10^4\\); it is not currently
implemented in lucifer but could be added via user-supplied localization
matrices.

The lucifer implementation uses a C++ EnKF filter (`enkf.cpp`) that
propagates ensemble members through the user-supplied drift and
diffusion functions with Euler-Maruyama substeps. Within the Gibbs loop,
the ensemble mean serves as the filtered state estimate for the
`theta_update_fn`. The ensemble size is set via `Specs$N_ensemble`
(default 100) and the covariance inflation factor via `Specs$inflation`
(default 1.0, no inflation).

The EnKF is the algorithm of choice when the state dimension exceeds
approximately 10, making the UKF’s \\O(n^3)\\ cost-per-step prohibitive.
For moderate state dimensions (\\n \< 10\\), the UKF is generally
preferred because its deterministic sigma-point approximation has lower
variance than the stochastic ensemble representation. For general
nonlinear models with low state dimension and non-Gaussian posteriors,
PGAS provides exact inference that neither the UKF nor the EnKF can
match.

#### Usage

``` r

library(lucifer)

# Simulate OU data (self-contained)
set.seed(42)
T_len <- 100
x_ou <- 5 + cumsum(rnorm(T_len, 0, 0.3))
y <- x_ou + rnorm(T_len, 0, 0.5)
sde <- SDE(family = "ou", data = y, times = seq(0.1, 10, length.out = T_len))

# MH parameter update
theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.05)
    if (any(proposal[c(1, 3, 4)] <= 0)) return(parm)
    parm  # placeholder: real application would evaluate likelihoods
}

fit_enkf <- SSM(sde, Algorithm = "EnKF", Iterations = 1000,
                Specs = list(N_ensemble = 100,
                             inflation = 1.02,
                             theta_update_fn = theta_fn),
                verbose = TRUE, Status = 250)

print(fit_enkf)
plot(fit_enkf, type = "states")
```

#### Extended example: multi-dimensional state with EnKF

The following example demonstrates the EnKF on a system with \\d_x = 5\\
states, where the full covariance matrix has 15 unique entries but the
ensemble representation avoids forming it explicitly.

``` r



# Simulate 5-dimensional coupled oscillator system
set.seed(88)
T_len <- 200
d_x <- 5
dt_sim <- 0.1
sigma_w <- 0.2
sigma_obs <- 0.5

# Coupling matrix: tridiagonal
A <- diag(-0.5, d_x)
for (i in 1:(d_x-1)) {
    A[i, i+1] <- 0.2
    A[i+1, i] <- 0.2
}

x_mat <- matrix(0, T_len, d_x)
x_mat[1, ] <- rnorm(d_x)
for (t in 2:T_len) {
    x_mat[t, ] <- x_mat[t-1, ] + A %*% x_mat[t-1, ] * dt_sim +
        sigma_w * rnorm(d_x, sd = sqrt(dt_sim))
}
# Observe only first 2 states
y_multi <- x_mat[, 1:2] + matrix(rnorm(T_len * 2, sd = sigma_obs),
                                  T_len, 2)

sde_multi <- SDE(
    data = y_multi,
    times = seq(0.1, by = dt_sim, length.out = T_len),
    x0 = rep(0, d_x),
    drift = function(x, theta, t) {
        d <- length(x)
        A_loc <- diag(theta[1], d)
        for (i in 1:(d-1)) {
            A_loc[i, i+1] <- theta[2]
            A_loc[i+1, i] <- theta[2]
        }
        as.vector(A_loc %*% x)
    },
    diffusion = function(x, theta, t) theta[3],
    obs.model = "gaussian",
    obs.noise = function(theta) theta[4],
    prior = function(theta) {
        if (any(theta[3:4] <= 0)) return(-Inf)
        sum(dnorm(theta, 0, 2, log = TRUE))
    },
    parm.names = c("self_coupling", "cross_coupling", "sigma_w", "sigma_obs"),
    state.names = paste0("x", 1:d_x),
    Initial.Values = c(-0.3, 0.1, 0.2, 0.4),
    method = "particle", N.particles = 50
)

theta_fn <- function(parm, states, Data) {
    parm + rnorm(length(parm), 0, 0.01)
}

fit_enkf_multi <- SSM(sde_multi, Algorithm = "EnKF", Iterations = 2000,
                      Specs = list(N_ensemble = 200,
                                   inflation = 1.05,
                                   theta_update_fn = theta_fn),
                      verbose = TRUE, Status = 500)

summary(fit_enkf_multi)
plot(fit_enkf_multi, type = "states")
plot(fit_enkf_multi, type = "parameters")
```

### Rao-Blackwellized particle filter (RBPF)

When the state vector admits a partition \\x_t = (x_t^{(nl)},
x_t^{(lin)})\\ such that, conditional on the nonlinear substate
\\x_t^{(nl)}\\, the linear substate \\x_t^{(lin)}\\ evolves according to
a linear-Gaussian model, the Rao-Blackwellized particle filter of
Doucet, de Freitas, Murphy, and Russell [\[12\]](#ref12) analytically
marginalizes the conditionally linear component via per-particle Kalman
filters. Each particle carries both a sample of the nonlinear substate
and a Kalman filter sufficient statistic \\(m^{(i)}\_t, P^{(i)}\_t)\\
for the linear substate conditional on the sampled nonlinear trajectory.
The marginal weight update integrates over the linear substate
analytically, so the effective dimension of the particle system is
\\\dim(x^{(nl)})\\ rather than \\\dim(x)\\, reducing variance without
increasing particle count.

Formally, the RBPF targets the marginal \\p(x\_{0:t}^{(nl)} \mid
y\_{1:t})\\ with particles and computes \\p(x\_{0:t}^{(lin)} \mid
x\_{0:t}^{(nl)}, y\_{1:t})\\ analytically via the Kalman filter. The
posterior mean and covariance of the full state vector are recovered by
marginalizing over the particle ensemble:

\\ E\[x_t \mid y\_{1:t}\] \approx \sum\_{i=1}^{N} w_t^{(i)}
\begin{pmatrix} x_t^{(nl),(i)} \\ m_t^{(i)} \end{pmatrix}. \\

The variance reduction from Rao-Blackwellization can be substantial. The
Rao-Blackwell theorem guarantees that \\\text{Var}\[E\[g(x) \mid
x^{(nl)}\]\] \leq \text{Var}\[g(x)\]\\ for any function \\g\\, with
equality only when \\g\\ does not depend on \\x^{(lin)}\\. In practice,
this means that an RBPF with \\N\\ particles can achieve the same
accuracy as a standard particle filter with \\N \cdot C\\ particles,
where \\C\\ depends on the dimensionality of the linear substate and the
signal-to-noise ratio. For models where \\\dim(x^{(lin)}) \gg
\dim(x^{(nl)})\\, such as jump-linear systems where the jump component
is scalar but the continuous state is multi-dimensional, the
computational savings are dramatic.

The approach is particularly effective for models such as conditionally
linear state-space models, mixed-frequency models, and jump-diffusion
processes where the jump component is nonlinear but the continuous
dynamics are linear conditional on the jump state. The cost is \\O(T N
d\_{lin}^3)\\ per iteration, with the cubic term coming from the
per-particle Kalman operations on the linear substate.

The lucifer implementation requires the user to specify the model
decomposition explicitly via `Specs$nl_transition_fn`,
`Specs$build_lin_ssm_fn`, `Specs$dim_nl`, and `Specs$dim_lin`. Because
this decomposition is model-specific and cannot be inferred
automatically, the RBPF is available only by explicit selection
(`Algorithm = "RBPF"`). The following example illustrates the interface
for a jump-linear model where the jump state selects between two linear
regimes.

``` r



# Jump-linear model: nonlinear jump state s_t in {-1, +1}
# selects between two drift regimes for the linear state x_t
# s_t = s_{t-1} * sign(rnorm(1))  [simplified jump dynamics]
# x_t = (1 + 0.5 * s_t) * x_{t-1} + w_t
# y_t = x_t + v_t

set.seed(42)
T_len <- 200
x_true <- numeric(T_len); s_true <- numeric(T_len)
x_true[1] <- 0; s_true[1] <- 1
sigma_w <- 0.3; sigma_v <- 1.0

for (t in 2:T_len) {
    s_true[t] <- s_true[t-1] * sample(c(1, -1), 1, prob = c(0.95, 0.05))
    x_true[t] <- (1 + 0.1 * s_true[t]) * x_true[t-1] + rnorm(1, sd = sigma_w)
}
y_jl <- x_true + rnorm(T_len, sd = sigma_v)

# Nonlinear transition for jump state
nl_trans <- function(s, theta, t) {
    sample(c(1, -1), 1, prob = c(theta[1], 1 - theta[1])) * sign(s)
}

# Linear SSM builder conditional on jump state
build_lin <- function(parm, s_nl) {
    F_mat <- array(1 + 0.1 * s_nl, c(1, 1, 1))
    list(F = F_mat, Q = array(parm[2]^2, c(1, 1, 1)),
         H = array(1, c(1, 1, 1)), R = array(parm[3]^2, c(1, 1, 1)),
         m0 = 0, P0 = matrix(10, 1, 1))
}

# Build SDE model for the RBPF (particle filter handles nonlinear part)
sde_jl <- SDE(
    drift = function(x, theta, t) 0 * x[1],
    diffusion = function(x, theta, t) theta[2],
    data = y_jl, times = 1:T_len, x0 = 0,
    obs.model = "gaussian",
    obs.noise = function(theta) theta[3],
    prior = function(theta) {
        if (any(theta[2:3] <= 0) || theta[1] < 0 || theta[1] > 1) return(-Inf)
        sum(dnorm(theta, c(0.9, 0.3, 1), c(0.2, 0.5, 0.5), log = TRUE))
    },
    parm.names = c("p_stay", "sigma_w", "sigma_v"),
    Initial.Values = c(0.9, 0.3, 1.0),
    method = "particle", N.particles = 30)

theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.01)
    if (any(proposal[2:3] <= 0) || proposal[1] < 0 || proposal[1] > 1) return(parm)
    proposal
}

fit_rbpf <- SSM(sde_jl, Algorithm = "RBPF", Iterations = 2000,
                Specs = list(N.particles = 200,
                             dim_nl = 1, dim_lin = 1,
                             nl_transition_fn = nl_trans,
                             build_lin_ssm_fn = build_lin,
                             theta_update_fn = theta_fn),
                verbose = TRUE, Status = 500)

print(fit_rbpf)
summary(fit_rbpf)
plot(fit_rbpf, type = "states")
plot(fit_rbpf, type = "parameters")
```

### Markov-switching FFBS (MS-FFBS)

Regime-switching state-space models augment the linear-Gaussian SSM with
a discrete latent Markov chain \\S_t \in \\1, \ldots, K\\\\ that selects
which set of system matrices governs the dynamics at each time step. The
observation and transition equations become

\\ x_t \mid x\_{t-1}, S_t = j \sim \mathcal{N}(F^{(j)} x\_{t-1},\\
Q^{(j)}), \qquad y_t \mid x_t, S_t = j \sim \mathcal{N}(H^{(j)} x_t,\\
R^{(j)}), \\

with the regime evolving according to a \\K \times K\\ transition matrix
\\P\_{ij} = \Pr(S\_{t+1} = j \mid S_t = i)\\.

The exact filter for this model has complexity that grows exponentially
in \\T\\ because the number of possible regime histories doubles at each
time step. Hamilton [\[13\]](#ref13) showed that a tractable
approximation can be obtained by collapsing the \\K^2\\ Gaussian
components produced at each filtered step back down to \\K\\ components,
an idea formalized by Kim [\[14\]](#ref14) into the Kim filter. At each
time step, \\K\\ parallel Kalman filters are run (one per current
regime), producing \\K^2\\ posterior components indexed by the previous
and current regime. These are collapsed to \\K\\ components by moment
matching, weighted by the filtered regime probabilities. The resulting
algorithm is \\O(T K^2 d_x^3)\\ per iteration.

The Kim collapsing approximation works as follows. At time \\t\\, the
joint filtered distribution is a mixture of \\K^2\\ Gaussian components:

\\ p(x_t \mid S\_{t-1} = i, S_t = j, y\_{1:t}) = \mathcal{N}(m_t^{ij},
P_t^{ij}), \\

with mixing weights proportional to \\\Pr(S\_{t-1} = i, S_t = j \mid
y\_{1:t})\\. The collapse to \\K\\ components is achieved by moment
matching within each value of \\S_t = j\\:

\\ m_t^j = \sum_i \pi\_{ij} m_t^{ij}, \qquad P_t^j = \sum_i \pi\_{ij}
\left\[ P_t^{ij} + (m_t^{ij} - m_t^j)(m_t^{ij} - m_t^j)^\top \right\],
\\

where \\\pi\_{ij} = \Pr(S\_{t-1} = i \mid S_t = j, y\_{1:t})\\ are the
normalized within-regime weights. This collapsing introduces a small
approximation error at each time step, but Kim and Nelson
[\[15\]](#ref15) showed through extensive simulation studies that the
cumulative approximation error is negligible for most practical
applications.

The MS-FFBS Gibbs sampler cycles through four steps: (1) run the Kim
filter forward to compute filtered state means, covariances, and regime
probabilities; (2) sample the regime sequence \\S\_{1:T}\\ backward from
the smoothed regime probabilities, conditioning on \\S\_{t+1}\\ at each
step via \\\Pr(S_t = i \mid S\_{t+1} = j, y\_{1:T}) \propto P\_{ij}
\cdot \Pr(S_t = i \mid y\_{1:t})\\; (3) conditional on the sampled
regime sequence, construct time-varying system matrices and draw the
full state trajectory \\x\_{0:T}\\ via standard FFBS; (4) update
parameters via `theta_update_fn(parm, states, regimes, Data)`, where
`regimes` is the integer vector of sampled regime labels; and (5) sample
each row of the transition matrix from its Dirichlet posterior given the
observed transition counts.

The user specifies the model via `Specs$build_ssm_list`, a list of \\K\\
functions each returning system matrices \\(F, Q, H, R, m_0, P_0)\\ for
its regime, and `Specs$n_regimes`. The `theta_update_fn` receives the
sampled regime vector as its third argument (after `parm` and `states`),
enabling regime-conditional parameter updates. The Dirichlet prior on
each row of the transition matrix is controlled by
`Specs$prior_transition` (default: uniform with concentration 1).

MS-FFBS should be used when the modeller believes that the
data-generating process switches between distinct regimes, each of which
can be described by a linear-Gaussian SSM. Common applications include
business cycle analysis (expansion vs. recession regimes with different
growth rates and volatilities), financial market modelling (calm
vs. turbulent regimes), and ecological time series where the population
dynamics shift between stable and unstable equilibria. The number of
regimes \\K\\ must be specified a priori; model comparison via marginal
likelihood or information criteria can guide the choice.

#### Usage

``` r



# Simulate a 2-regime local level model
set.seed(42)
T_len <- 300
K <- 2
P_true <- matrix(c(0.98, 0.02, 0.05, 0.95), 2, 2, byrow = TRUE)
sigma_w <- c(0.1, 0.8)   # low vs. high state noise
sigma_v <- 1.0

# Generate regime sequence
S_true <- integer(T_len); S_true[1] <- 1
for (t in 2:T_len)
    S_true[t] <- sample(1:K, 1, prob = P_true[S_true[t-1], ])

# Generate states and observations
x_true <- numeric(T_len); x_true[1] <- 0
for (t in 2:T_len)
    x_true[t] <- x_true[t-1] + rnorm(1, 0, sigma_w[S_true[t]])
y_ms <- x_true + rnorm(T_len, 0, sigma_v)

# Build model via SSM_model (no SDE needed for discrete-time SSM)
sde_ms <- SSM_model(
    type = "local_level", data = y_ms, times = 1:T_len,
    parm.names = c("sigma_w1", "sigma_w2", "sigma_v"),
    Initial.Values = c(0.2, 0.5, 1.0))

# Regime-conditional builders
build_ssm_list <- list(
    function(parm) {
        list(F = array(1, c(1,1,1)), Q = array(parm[1]^2, c(1,1,1)),
             H = array(1, c(1,1,1)), R = array(parm[3]^2, c(1,1,1)),
             m0 = 0, P0 = matrix(10, 1, 1))
    },
    function(parm) {
        list(F = array(1, c(1,1,1)), Q = array(parm[2]^2, c(1,1,1)),
             H = array(1, c(1,1,1)), R = array(parm[3]^2, c(1,1,1)),
             m0 = 0, P0 = matrix(10, 1, 1))
    }
)

# Parameter update (regime-aware: receives regimes as 3rd arg)
theta_fn <- function(parm, states, regimes, Data) {
    x <- states[, 1]
    dx <- diff(x)
    # Regime-conditional state noise
    for (j in 1:2) {
        idx <- which(regimes[-1] == j)
        if (length(idx) > 2) {
            ss <- sum(dx[idx]^2)
            parm[j] <- sqrt(1 / rgamma(1, length(idx)/2, ss/2))
        }
    }
    # Observation noise
    resid <- as.vector(Data$.y) - x
    ss_v <- sum(resid^2)
    parm[3] <- sqrt(1 / rgamma(1, length(resid)/2, ss_v/2))
    parm
}

fit_ms <- SSM(sde_ms, Algorithm = "MS-FFBS", Iterations = 3000,
              Specs = list(n_regimes = 2,
                           build_ssm_list = build_ssm_list,
                           theta_update_fn = theta_fn),
              verbose = TRUE, Status = 500)

print(fit_ms)
summary(fit_ms)
plot(fit_ms, type = "states")
plot(fit_ms, type = "parameters")
```

#### Extended example: three-regime business cycle model

This example extends the regime-switching framework to three regimes
(recession, normal growth, boom), each with different trend growth rates
and volatilities. This is a simplified version of the Hamilton
[\[13\]](#ref13) business cycle model.

``` r



# Simulate 3-regime local linear trend
set.seed(55)
T_len <- 400
K <- 3
P_true <- matrix(c(
    0.95, 0.04, 0.01,
    0.03, 0.94, 0.03,
    0.02, 0.05, 0.93
), 3, 3, byrow = TRUE)

# Regime-specific parameters
drift_regime <- c(-0.5, 0.2, 0.8)  # recession, normal, boom
sigma_regime <- c(1.5, 0.5, 0.8)   # high, low, medium vol

S <- integer(T_len); S[1] <- 2
for (t in 2:T_len) S[t] <- sample(1:K, 1, prob = P_true[S[t-1], ])

x <- numeric(T_len); x[1] <- 0
for (t in 2:T_len)
    x[t] <- x[t-1] + drift_regime[S[t]] + rnorm(1, sd = sigma_regime[S[t]])
y_bc <- x + rnorm(T_len, sd = 0.8)

# Build model via SSM_model
sde_bc <- SSM_model(
    type = "local_level", data = y_bc, times = 1:T_len,
    parm.names = c("sigma_w1", "sigma_w2", "sigma_w3", "sigma_v"),
    Initial.Values = c(1.0, 0.5, 0.8, 0.8))

# Three regime builders with different volatilities
build_ssm_3 <- list(
    function(parm) list(F = array(1, c(1,1,1)), Q = array(parm[1]^2, c(1,1,1)),
                        H = array(1, c(1,1,1)), R = array(parm[4]^2, c(1,1,1)),
                        m0 = 0, P0 = matrix(10, 1, 1)),
    function(parm) list(F = array(1, c(1,1,1)), Q = array(parm[2]^2, c(1,1,1)),
                        H = array(1, c(1,1,1)), R = array(parm[4]^2, c(1,1,1)),
                        m0 = 0, P0 = matrix(10, 1, 1)),
    function(parm) list(F = array(1, c(1,1,1)), Q = array(parm[3]^2, c(1,1,1)),
                        H = array(1, c(1,1,1)), R = array(parm[4]^2, c(1,1,1)),
                        m0 = 0, P0 = matrix(10, 1, 1))
)

theta_fn_3 <- function(parm, states, regimes, Data) {
    x <- states[, 1]
    dx <- diff(x)
    for (j in 1:3) {
        idx <- which(regimes[-1] == j)
        if (length(idx) > 2) {
            ss <- sum(dx[idx]^2)
            parm[j] <- sqrt(1 / rgamma(1, length(idx)/2, ss/2))
        }
    }
    resid <- as.vector(Data$.y) - x
    parm[4] <- sqrt(1 / rgamma(1, length(resid)/2, sum(resid^2)/2))
    parm
}

fit_bc <- SSM(sde_bc, Algorithm = "MS-FFBS", Iterations = 5000,
              Specs = list(n_regimes = 3,
                           build_ssm_list = build_ssm_3,
                           theta_update_fn = theta_fn_3,
                           prior_transition = rep(1, 3)),
              verbose = TRUE, Status = 1000)

summary(fit_bc)
plot(fit_bc, type = "states")
plot(fit_bc, type = "parameters")
```

## Non-Gaussian observation families

The default observation model in SSM inference is Gaussian: \\y_t \mid
x_t \sim \mathcal{N}(Hx_t, R)\\. For count data, proportions, or
heavy-tailed observations, the SSM module provides a registry of
non-Gaussian observation families that can be used with particle-based
algorithms (PGAS, SMC\\^2\\, RBPF, EnKF). The registry is accessed via
the `obs.model` argument to
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md),
which accepts a string naming the family, and `obs.params`, a named list
providing the required dispersion or scale parameters.

The available families are `"gaussian"` (the default, parameterized by
`sigma`), `"poisson"` (log link, no scale parameter), `"negbin"` (log
link, parameterized by `size`), `"student_t"` (identity link,
parameterized by `sigma` and `df`), `"zero_inflated_poisson"` (log link,
parameterized by `p_zero`), and `"binomial"` (logit link, parameterized
by `n_trials`). Each family defines a log-likelihood function, a link
function mapping the latent state to the natural parameter of the
observation distribution, and metadata about which parameters must be
supplied. The observation parameters in `obs.params` can be either fixed
scalars or functions of the MCMC parameter vector \\\theta\\; the latter
enables joint estimation of the observation distribution parameters
alongside the latent states and transition parameters.

The following table summarizes the available observation families, their
link functions, and required parameters.

| Family    | Link function | Mean function | Extra parameters | Applicable algorithms |
|:----------|:--------------|:--------------|:-----------------|:----------------------|
| gaussian  | Identity      | x_t           | sigma            | All                   |
| poisson   | Log           | exp(x_t)      | none             | PGAS, SMC2            |
| negbin    | Log           | exp(x_t)      | size             | PGAS, SMC2            |
| student_t | Identity      | x_t           | sigma, df        | PGAS, SMC2            |
| zip       | Log           | exp(x_t)      | p_zero           | PGAS, SMC2            |
| binomial  | Logit         | expit(x_t)    | n_trials         | PGAS, SMC2            |

Non-Gaussian observation families in the SSM module. {.table}

For example, to fit a state-space model with Poisson observations where
the latent state represents the log-intensity, one would specify:

``` r

# Simulate Poisson count data from a latent OU process
set.seed(42)
T_pois <- 100
x_lat <- numeric(T_pois); x_lat[1] <- 2
for (t in 2:T_pois) x_lat[t] <- x_lat[t-1] + 0.5 * (2 - x_lat[t-1]) * 0.1 +
    0.3 * rnorm(1, sd = sqrt(0.1))
y_counts <- rpois(T_pois, lambda = exp(x_lat))

sde_pois <- SDE(
    data = y_counts,
    times = seq(0.1, 10, length.out = T_pois),
    x0 = 2,
    drift = function(x, theta, t) theta[1] * (theta[2] - x[1]),
    diffusion = function(x, theta, t) theta[3],
    obs.model = "poisson",
    obs.link = function(x, theta) exp(x[1]),
    prior = function(theta) {
        if (any(theta[c(1, 3)] <= 0)) return(-Inf)
        sum(dnorm(theta, c(0.5, 2, 0.3), 1, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma"),
    Initial.Values = c(0.5, 1.0, 0.3),
    method = "particle", N.particles = 50
)

theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.02)
    if (any(proposal[c(1, 3)] <= 0)) return(parm)
    proposal
}

fit_pois <- SSM(sde_pois, Algorithm = "PGAS",
                Specs = list(N.particles = 50,
                             theta_update_fn = theta_fn))
```

For negative binomial observations where the overdispersion parameter is
estimated jointly:

``` r

# Same latent process, but with negative binomial observations
sde_nb <- SDE(
    data = y_counts,
    times = seq(0.1, 10, length.out = length(y_counts)),
    x0 = 2,
    drift = function(x, theta, t) theta[1] * (theta[2] - x[1]),
    diffusion = function(x, theta, t) theta[3],
    obs.model = "negbin",
    obs.link = function(x, theta) exp(x[1]),
    obs.params = list(size = function(theta) exp(theta[4])),
    prior = function(theta) {
        if (any(theta[c(1, 3)] <= 0)) return(-Inf)
        sum(dnorm(theta, c(0.5, 2, 0.3, 2), 1, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma", "log_size"),
    Initial.Values = c(0.5, 1.0, 0.3, 2.0),
    method = "particle", N.particles = 50
)
```

For Student-\\t\\ observations that provide robustness to outliers:

``` r

# Simulate data with outliers
set.seed(42)
T_out <- 100
x_out <- cumsum(rnorm(T_out, 0, 0.3))
y_outliers <- x_out + rt(T_out, df = 3) * 0.5

sde_t <- SDE(
    data = y_outliers,
    times = seq(0.1, 10, length.out = T_out),
    x0 = 0,
    drift = function(x, theta, t) theta[1] * (theta[2] - x[1]),
    diffusion = function(x, theta, t) theta[3],
    obs.model = "student_t",
    obs.params = list(sigma = 1.0, df = 5),
    prior = function(theta) {
        if (any(theta[c(1, 3)] <= 0)) return(-Inf)
        sum(dnorm(theta, 0, 2, log = TRUE))
    },
    parm.names = c("kappa", "mu", "sigma_state"),
    Initial.Values = c(0.5, 0.0, 0.3),
    method = "particle", N.particles = 50
)

theta_fn <- function(parm, states, Data) {
    proposal <- parm + rnorm(length(parm), 0, 0.02)
    if (any(proposal[c(1, 3)] <= 0)) return(parm)
    proposal
}

fit_t <- SSM(sde_t, Algorithm = "PGAS",
             Specs = list(N.particles = 50,
                          theta_update_fn = theta_fn))
```

## Post-inference diagnostics

The `ssm_fit` object returned by
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
supports conversion to a `demonoid` object via
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md),
enabling the full suite of lucifer diagnostic tools. This section
demonstrates the key diagnostic workflows for state-space model
inference.

### Posterior predictive checks

Posterior predictive checks (PPCs) assess whether the fitted model can
generate data that resemble the observed data. For state-space models,
PPCs can operate at two levels: checking whether the marginal
distribution of \\y_t\\ is well-calibrated, and checking whether the
temporal structure (autocorrelation, spectral density) is captured. The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
conversion extracts the parameter posterior, which can then be used for
PPCs.

``` r

# After fitting an SSM, convert to demonoid for PPC
demon <- as.demonoid(fit)

# Posterior predictive p-values: compare test statistics
# computed on observed data against the posterior predictive distribution
# The yhat field in the demonoid contains posterior predictive draws

# Residual analysis: innovation sequence diagnostics
# For FFBS/KSC, the innovations v_t = y_t - H m_{t|t-1} should be
# approximately i.i.d. N(0, S_t) if the model is correct
states_mean <- apply(fit$States, c(2, 3), mean)
residuals <- as.vector(y) - states_mean[, 1]

# Visual diagnostics
par(mfrow = c(2, 2))
plot(residuals, type = "l", main = "Innovation sequence")
acf(residuals, main = "ACF of innovations")
qqnorm(residuals); qqline(residuals)
hist(residuals, breaks = 30, freq = FALSE, main = "Innovation density")
curve(dnorm(x, mean(residuals), sd(residuals)), add = TRUE, col = "red")
```

A well-specified model should produce innovations that are approximately
white noise (no significant autocorrelation) and approximately Gaussian
(linear Q-Q plot). Systematic patterns in the innovation sequence, such
as heteroscedasticity or trending mean, indicate model misspecification.
Serial correlation in the innovations suggests that the state dynamics
are insufficiently flexible, while heavy tails suggest that the
observation noise is non-Gaussian and a Student-\\t\\ or other
heavy-tailed observation family should be considered.

The standardized innovations \\e_t = \nu_t / \sqrt{S_t}\\, where
\\\nu_t\\ is the innovation and \\S_t\\ is the innovation covariance
from the Kalman filter, should follow a standard normal distribution
under correct model specification. The Ljung-Box test applied to the
squared standardized innovations provides a formal test for remaining
heteroscedasticity, while the Jarque-Bera test checks for departures
from normality. For multivariate observations, the Mahalanobis distance
of the innovation vector \\d_t = \nu_t^\top S_t^{-1} \nu_t\\ should
follow a \\\chi^2\_{d_y}\\ distribution, providing a joint diagnostic
for all observation dimensions simultaneously. These diagnostics
complement the informal visual checks shown above and should be reported
whenever the fitted model is used for forecasting or decision-making.

``` r

# Standardized innovation diagnostics (requires the filter output)
# Ljung-Box test for serial correlation in squared innovations
e_std <- residuals / sd(residuals)
Box.test(e_std^2, lag = 10, type = "Ljung-Box")

# Cumulative sum (CUSUM) test for structural stability
cusum <- cumsum(e_std) / sqrt(length(e_std))
plot(cusum, type = "l", main = "CUSUM of standardized innovations")
abline(h = c(-1.96, 1.96), lty = 2, col = "red")
```

### Model comparison via LOO-PSIS

LOO-PSIS (Leave-One-Out cross-validation via Pareto-Smoothed Importance
Sampling) provides a principled approach to comparing competing SSM
specifications without refitting the model for each left-out
observation. The method requires pointwise log-likelihood evaluations,
which are available from the
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
conversion. Comparing a local level model against a BSM, for instance,
reveals whether the additional complexity of trend and seasonal
components is warranted by the data.

``` r

# Fit two competing models to the same data
mod_ll  <- SSM_model(type = "local_level", data = y_bsm)
mod_bsm <- SSM_model(type = "bsm", data = y_bsm, period = 4)

fit_ll  <- SSM(mod_ll,  Iterations = 3000)
fit_bsm <- SSM(mod_bsm, Iterations = 3000)

# Convert to demonoid for LOO
demon_ll  <- as.demonoid(fit_ll)
demon_bsm <- as.demonoid(fit_bsm)

# Compute LOO-PSIS
loo_ll  <- LOO(demon_ll)
loo_bsm <- LOO(demon_bsm)

# Compare: negative ELPD difference favors the first model
loo_compare(loo_ll, loo_bsm)
```

The
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)
function returns the difference in expected log predictive density
(ELPD) between the two models along with a standard error. A positive
difference favors the second model; a negative difference favors the
first. The Pareto \\k\\ diagnostics reported by
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
indicate the reliability of the importance sampling approximation: \\k
\< 0.5\\ is good, \\0.5 \< k \< 0.7\\ is acceptable, and \\k \> 0.7\\
suggests that the importance weights are too variable for reliable LOO
estimation and moment matching or refitting may be needed
[\[18\]](#ref18).

### Robust Bayesian sensitivity analysis

The
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
function performs power-scaling sensitivity analysis to assess whether
the posterior is driven by the data or the prior. For state-space
models, this is particularly relevant when informative priors are placed
on the noise variances \\\sigma_w^2\\ and \\\sigma_v^2\\, which can
strongly influence the signal-to-noise ratio and hence the smoothness of
the estimated state trajectory.

``` r

# Power-scaling sensitivity analysis
demon <- as.demonoid(fit)
rb <- RobustBayes(demon, type = "power")
summary(rb)
plot(rb)
```

The power-scaling approach perturbs the prior by raising it to a power
\\\alpha\\ and the likelihood to a power \\\beta\\, then examines how
the posterior changes. Parameters whose posteriors are sensitive to the
prior power (large derivative \\\partial E\[\theta\] / \partial \alpha\\
at \\\alpha = 1\\) are “prior-dominated” and should be interpreted
cautiously. For SSMs, the observation noise variance \\\sigma_v^2\\ is
typically well-identified by the data (likelihood-dominated), while the
process noise variance \\\sigma_w^2\\ can be prior-sensitive, especially
when the time series is short or the signal-to-noise ratio is low.

### Convergence diagnostics

Standard MCMC convergence diagnostics apply directly to the parameter
chains extracted from the `ssm_fit` object. The
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
and
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md)
functions from lucifer provide the split-\\\hat{R}\\ diagnostic of
Vehtari et al. [\[19\]](#ref19) and the effective sample size, both of
which are essential for verifying that the Gibbs sampler has converged.

``` r

# Convert to demonoid for convergence diagnostics
demon <- as.demonoid(fit)

# Split-Rhat: values near 1.0 indicate convergence
Rhat(demon)

# Effective sample size: should be at least 100 for reliable inference
ESS(demon)

# Trace plots of parameters
plot(fit, type = "parameters")

# State trajectory with credible intervals
plot(fit, type = "states")
```

For FFBS with conjugate updates, convergence is typically rapid, and
\\\hat{R} \< 1.01\\ is achieved within 1000 iterations. For PGAS, mixing
depends on the number of particles; increasing `N.particles` reduces the
autocorrelation of the parameter chain at the cost of increased
computation per iteration. A useful diagnostic is to compare the ESS per
second across different particle counts to find the computationally
optimal setting. For KSC, the persistence parameter \\\phi\\ often has
the slowest mixing; running longer chains (5000+) with moderate thinning
(5-10) is recommended.

## The `ssm_fit` object

All eight algorithms return an S3 object of class `ssm_fit` containing
three core components. The `Posterior` matrix (\\n\_{\text{samples}}
\times d\_\theta\\) stores the parameter samples, compatible with
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
for downstream diagnostics such as Rhat, ESS, LOO, and robust Bayesian
sensitivity analysis. The `States` array (\\n\_{\text{samples}} \times T
\times d_x\\) stores the sampled state trajectories, from which
pointwise credible intervals and smoothed means can be extracted. The
`Log.Marginal.Likelihood` scalar provides an estimate of \\\log p(y)\\
where available (exact for FFBS, approximate for PGAS and SMC\\^2\\,
unavailable for KSC).

The S3 methods provide a consistent interface across all eight
algorithms.

``` r

# S3 methods
print(fit)                          # concise one-screen summary
summary(fit)                        # extended diagnostics with ESS
plot(fit, type = "states")          # state trajectories with 95% bands
plot(fit, type = "parameters")      # trace plots and densities
plot(fit, type = "pairs")           # bivariate scatterplots

# Convert to demonoid for lucifer diagnostics
demon <- as.demonoid(fit)
```

The [`print()`](https://rdrr.io/r/base/print.html) method displays a
compact summary including the algorithm used, the number of iterations,
the number of parameters and state dimensions, and the log marginal
likelihood estimate if available. The
[`summary()`](https://rdrr.io/r/base/summary.html) method provides
extended output with posterior means, standard deviations, quantiles
(2.5%, 50%, 97.5%), Rhat, and ESS for each parameter. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method supports
three types: `"states"` displays the posterior mean state trajectory
with 95% credible intervals overlaid on the observations; `"parameters"`
shows trace plots and marginal density estimates for each parameter; and
`"pairs"` produces bivariate scatterplots of the posterior samples,
useful for detecting posterior correlation between parameters.

The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
conversion creates a `demonoid` object from the parameter posterior,
enabling the use of all standard lucifer diagnostic functions. This
bridge is essential for integrating SSM inference with the broader
lucifer ecosystem, including model comparison via
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
and
[`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md),
prior sensitivity via
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md),
and convergence diagnostics via
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
and
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md).

## Algorithm selection guide

| Criterion | FFBS | UKF | EnKF | RBPF | PGAS | SMC2 | MS-FFBS | KSC |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Model class | Linear-Gaussian | Nonlinear Gaussian | Nonlinear Gaussian | Mixed linear/nonlinear | General nonlinear | Most general | Regime-switching LG | Stochastic volatility |
| Requires theta_update_fn | Yes | Yes | Yes | Yes | Yes | No | Yes | No (built-in) |
| Requires gradients | No | No | No | No | No | No | No | No |
| Computational cost per iteration | O(T d_x^3) | O(T d_x^3) | O(T N_e d_x) | O(T N d_lin^3) | O(T N d_x) | O(T N_out N_in d_x) | O(T K^2 d_x^3) | O(T d_x^3) |
| Marginal likelihood | Exact | Approximate | Approximate | Approximate | Approximate | Unbiased | Approximate | Not available |
| Best for | DLMs, local level, structural time series | Moderate nonlinearity, low dimension | High-dimensional states (d_x \>\> 10) | Conditionally linear substates | Nonlinear SDEs, ecological models | No conjugacy, black-box models | Regime changes, structural breaks | Financial volatility, SV models |

Comparison of SSM inference algorithms in lucifer. {.table
style="width:100%;"}

For linear-Gaussian models where conjugate full conditionals for
\\\theta\\ are available, FFBS is the clear choice: it produces exact
Gibbs samples with minimal computational cost and no tuning parameters.
When the model is nonlinear but the state dimension is moderate and the
nonlinearity is smooth, the UKF provides an efficient deterministic
approximation that avoids the Monte Carlo noise of particle filters; for
high-dimensional states where full covariance propagation is
intractable, the EnKF scales gracefully by replacing the covariance
matrix with an ensemble representation. When the model has mixed
linear/nonlinear structure, the RBPF exploits this decomposition to
reduce the effective particle dimension. For general nonlinear models
where tractable transition densities are available (or can be
approximated via Euler-Maruyama), PGAS provides a well-mixing Gibbs
sampler that requires only a `theta_update_fn` and a moderate number of
particles (\\N = 50\\-\\200\\). When neither full conditionals nor
tractable transitions are available, SMC\\^2\\ handles the fully general
case at the cost of running \\N\_{\text{outer}} \times
N\_{\text{inner}}\\ model evaluations per time step. For
regime-switching models with linear-Gaussian dynamics within each
regime, MS-FFBS provides efficient Gibbs sampling with automatic
transition matrix estimation. For the standard stochastic volatility
model specifically, KSC is the most efficient choice because it exploits
the mixture approximation to avoid particle filtering entirely.

The following decision tree provides practical guidance for algorithm
selection. Begin by asking whether the model is linear-Gaussian. If yes,
use FFBS (or MS-FFBS if regime-switching is present). If the model is
nonlinear but the observation model is Gaussian, check the state
dimension: for \\d_x \< 10\\, use UKF; for \\d_x \geq 10\\, use EnKF. If
the model has non-Gaussian observations and a tractable transition
density, use PGAS. If the model has non-Gaussian observations and an
intractable transition density, use SMC\\^2\\. If the model is
specifically a standard SV model, use KSC. If the model has mixed
linear/nonlinear structure that can be decomposed, consider RBPF for
variance reduction.

The computational cost comparisons in the table above should be
interpreted carefully. The notation \\O(T d_x^3)\\ for FFBS reflects the
per-iteration cost of the Kalman filter with matrix inversions; in
practice, for scalar state models (\\d_x = 1\\), FFBS processes
thousands of time steps in milliseconds. The notation \\O(T N d_x)\\ for
PGAS reflects the per-particle propagation cost, but the constant factor
includes R callback overhead for custom drift/diffusion functions
(mitigated in lucifer by the C++ backend). For SMC\\^2\\, the cost is
dominated by the \\N\_{\text{outer}} \times N\_{\text{inner}}\\ particle
filters that must be maintained simultaneously; a typical configuration
with \\N\_{\text{outer}} = 200\\ and \\N\_{\text{inner}} = 50\\ requires
10,000 model evaluations per time step, which can be prohibitive for
long time series (\\T \> 1000\\).

## Practical considerations

### Choosing the number of particles

For PGAS, the number of particles \\N\\ controls the tradeoff between
mixing quality and computational cost. Too few particles lead to path
degeneracy, where the sampled trajectory fails to explore the state
space and the parameter chain exhibits high autocorrelation. Too many
particles increase the per-iteration cost without proportional
improvement in mixing. A practical approach is to start with \\N = 50\\
and monitor the ESS of the parameter chain; if ESS/iteration is below
0.1 (i.e., ESS \\\<\\ 100 after 1000 iterations), increase the particle
count. For well-behaved models with unimodal state posteriors, \\N =
50\\-\\100\\ is typically sufficient; for multimodal or weakly
identified models, \\N = 200\\-\\500\\ may be needed.

For SMC\\^2\\, both \\N\_{\text{outer}}\\ and \\N\_{\text{inner}}\\ must
be tuned. The outer particle count determines the resolution of the
parameter posterior approximation; \\N\_{\text{outer}} = 100\\-\\500\\
is typical. The inner particle count determines the accuracy of each
parameter particle’s likelihood estimate; \\N\_{\text{inner}} =
20\\-\\50\\ is usually sufficient. The rejuvenation acceptance rate
provides a useful diagnostic: rates below 10% suggest that the proposal
distribution is poorly calibrated (try increasing `rejuvenation_steps`
or reducing the proposal variance), while rates above 50% suggest that
more outer particles could improve the posterior approximation.

### Handling missing observations

All algorithms in the SSM module handle missing observations
transparently. When \\y_t\\ contains `NA`, the Kalman filter (used by
FFBS, UKF, and as a sub-component of RBPF, KSC, and MS-FFBS) skips the
update step at time \\t\\, propagating only the prediction step. This is
equivalent to setting the observation precision to zero, so that the
filtered state at time \\t\\ depends entirely on the prior (prediction
from \\t-1\\) rather than the observation. For particle-based methods
(PGAS, SMC\\^2\\), missing observations are handled by setting the
observation weight to 1 (uniform), so that particles are propagated
without reweighting.

Missing observations can occur for several reasons: sensor failure,
irregular sampling, or deliberate design (e.g., hold-out sets for
cross-validation). The Kalman filter’s handling of missing data provides
a natural mechanism for temporal interpolation: the smoothed state
estimate at a time with missing data is a weighted average of the
surrounding observations, with weights determined by the state dynamics
and noise variances.

### Prior specification for noise variances

The choice of prior on the noise variances \\\sigma_w^2\\ and
\\\sigma_v^2\\ can substantially influence the estimated signal-to-noise
ratio, especially for short time series. The conjugate prior for a
variance parameter in the linear-Gaussian SSM is inverse-gamma,
\\\sigma^2 \sim \text{InvGamma}(\alpha_0, \beta_0)\\, with the special
case \\\alpha_0 = \beta_0 = 0.001\\ providing a diffuse but proper
prior. For weakly informative priors, half-Cauchy or half-\\t\\ priors
on \\\sigma\\ (rather than \\\sigma^2\\) are often recommended because
they place less mass near zero while allowing large values; these can be
implemented via MH steps in the `theta_update_fn`. The
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
diagnostic described above can quantify the sensitivity of the posterior
to the prior choice.

The inverse-gamma prior \\\text{InvGamma}(\alpha_0, \beta_0)\\ has
density \\p(\sigma^2) \propto (\sigma^2)^{-\alpha_0 - 1} \exp(-\beta_0 /
\sigma^2)\\, with mean \\\beta_0 / (\alpha_0 - 1)\\ for \\\alpha_0 \>
1\\ and mode \\\beta_0 / (\alpha_0 + 1)\\. Common choices include
\\\alpha_0 = \beta_0 = 0.001\\ (diffuse), \\\alpha_0 = 2, \beta_0 = 1\\
(weakly informative with mode 1/3 and mean 1), and \\\alpha_0 = 5,
\beta_0 = s_0\\ (moderately informative, centering the prior around
\\s_0 / 4\\). For the built-in SSM builders, the default is \\\alpha_0 =
\beta_0 = 0.001\\, but users can override this in their
`theta_update_fn` by adjusting the shape and rate parameters of the
[`rgamma()`](https://rdrr.io/r/stats/GammaDist.html) calls.

### Reparameterization for improved mixing

Gibbs samplers for SSMs can exhibit slow mixing when the parameters are
strongly correlated in the posterior. A common source of posterior
correlation is the interaction between the state trajectory and the
noise variances: when \\\sigma_w^2\\ is large, the state trajectory is
flexible and tracks the observations closely, which makes the
observation noise residuals small and drives \\\sigma_v^2\\ toward zero,
and vice versa. This creates a “banana-shaped” posterior that Gibbs
sampling traverses slowly because each conditional update moves along
the narrow direction of the posterior.

Several reparameterizations can mitigate this problem. The
**non-centered parameterization** replaces the state \\x_t\\ with the
standardized innovation \\z_t = (x_t - f(x\_{t-1}, \theta)) /
\sigma_w\\, so that the state trajectory is parameterized as \\x_t =
f(x\_{t-1}, \theta) + \sigma_w z_t\\ with \\z_t \sim \mathcal{N}(0,
1)\\. This decouples the state trajectory from \\\sigma_w\\ in the
prior, often improving mixing when the observation noise is small
relative to the process noise. The **interweaving strategy** of Yu and
Meng alternates between centered and non-centered parameterizations
within a single Gibbs cycle, guaranteeing that at least one of the two
parameterizations mixes well regardless of the signal-to-noise ratio.

Another effective approach is to reparameterize in terms of the
signal-to-noise ratio \\q = \sigma_w^2 / \sigma_v^2\\ and the total
noise \\\tau^2 = \sigma_w^2 + \sigma_v^2\\, which are often less
correlated in the posterior than the individual variances. The
`theta_update_fn` can implement any of these reparameterizations by
transforming the parameters before and after the Gibbs update step.

### Forecasting and prediction

After fitting an SSM, the estimated model can be used for
multi-step-ahead forecasting by propagating the filtered state
distribution forward without observations. For linear-Gaussian models,
the forecast distribution at horizon \\h\\ is

\\ p(y\_{T+h} \mid y\_{1:T}, \theta) = \mathcal{N}(H F^h m_T, \\ H (F^h
P_T (F^h)^\top + \sum\_{j=0}^{h-1} F^j Q (F^j)^\top) H^\top + R), \\

where \\m_T\\ and \\P_T\\ are the filtered mean and covariance at the
last observation. The forecast uncertainty grows with horizon \\h\\ as
the accumulated process noise inflates the predicted state covariance.

For Bayesian forecasts that account for parameter uncertainty, one
integrates over the posterior parameter distribution. In practice, this
is achieved by computing the forecast for each posterior draw
\\\theta^{(i)}\\ and aggregating: the predictive mean is
\\\hat{y}\_{T+h} = \frac{1}{S} \sum\_{i=1}^S E\[y\_{T+h} \mid y\_{1:T},
\theta^{(i)}\]\\ and the predictive variance includes both the
within-model forecast variance and the between-model variance due to
parameter uncertainty.

``` r

# Multi-step-ahead forecasting from a fitted local level model
# After fitting fit_ll from the FFBS section

# Extract posterior parameter draws
sigma_w_draws <- fit_ll$Posterior[, "sigma_w"]
sigma_v_draws <- fit_ll$Posterior[, "sigma_v"]
states_final  <- fit_ll$States[, nrow(fit_ll$States[1,,]), 1]

# Generate forecasts for each posterior draw
H <- 20  # forecast horizon
n_draws <- nrow(fit_ll$Posterior)
forecasts <- matrix(NA, n_draws, H)

for (i in 1:n_draws) {
    x_curr <- states_final[i]
    for (h in 1:H) {
        x_curr <- x_curr + rnorm(1, sd = sigma_w_draws[i])
        forecasts[i, h] <- x_curr + rnorm(1, sd = sigma_v_draws[i])
    }
}

# Predictive mean and 95% interval
forecast_mean <- colMeans(forecasts)
forecast_lower <- apply(forecasts, 2, quantile, 0.025)
forecast_upper <- apply(forecasts, 2, quantile, 0.975)

# Plot forecast fan chart
plot(1:H, forecast_mean, type = "l", lwd = 2,
     xlab = "Forecast horizon", ylab = "y",
     ylim = range(c(forecast_lower, forecast_upper)),
     main = "Multi-step-ahead forecast with 95% predictive interval")
polygon(c(1:H, H:1), c(forecast_lower, rev(forecast_upper)),
        col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
lines(1:H, forecast_mean, lwd = 2, col = "darkblue")
```

### Connection to SDE inference

The SSM module and the SDE module in lucifer share a common model
specification interface through the `sde_model` object, but they serve
different inference goals. The
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
function performs parameter estimation for continuous-time SDEs using
particle MCMC or pseudo-marginal methods, marginalizing over the latent
path. The
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
function performs joint inference on both parameters and the discretized
latent path, producing explicit state trajectory samples. The choice
between them depends on whether the latent trajectory itself is of
scientific interest.

When the primary goal is parameter estimation and the latent path is a
nuisance,
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
is appropriate because it integrates out the path and avoids the mixing
challenges associated with sampling high-dimensional state trajectories.
When the latent path is of direct interest, as in signal extraction,
smoothing, or state-dependent forecasting,
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md) is
the correct choice because it produces posterior samples of the full
trajectory \\x\_{0:T}\\.

The two approaches can also be combined for model checking: fit the
parameters via
[`SDE.fit()`](https://robustecologies.github.io/lucifer/reference/SDE.fit.md)
to obtain a well-mixing parameter chain, then condition on the posterior
mean parameters and use
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
with fixed parameters to sample the state trajectory. This two-stage
approach can be more efficient than joint inference when the parameter
posterior is highly non-Gaussian or multimodal.

It is worth noting that the SDE family registry plays a dual role in
this architecture. When an `sde_model` object is passed to
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md),
the family metadata determines which algorithm is auto-selected (e.g.,
`"sv"` triggers KSC, `"ou"` triggers FFBS) and provides default system
matrix builders for families with known exact discretizations. The
`"custom"` family, which accepts user-supplied drift and diffusion
functions, is the general-purpose entry point for non-standard models
and defaults to PGAS or SMC\\^2\\ depending on whether a
`theta_update_fn` is provided. Users working with discrete-time models
that have no continuous-time interpretation should prefer
[`SSM_model()`](https://robustecologies.github.io/lucifer/reference/SSM_model.md)
over
[`SDE()`](https://robustecologies.github.io/lucifer/reference/SDE.md),
as it avoids unnecessary discretization overhead and provides a cleaner
specification interface.

## Complete workflow example

This section demonstrates a complete analysis workflow on a simulated
ecological time series, from model specification through inference,
diagnostics, model comparison, and forecasting. The data represent
annual abundance observations of a hypothetical population undergoing
density-dependent stochastic growth.

``` r



# Simulate Gompertz state-space model
# x_t = a + b * x_{t-1} + w_t,  w_t ~ N(0, sigma_w^2)
# y_t = x_t + v_t,              v_t ~ N(0, sigma_v^2)
# where x_t = log(N_t) is log-abundance

set.seed(42)
T_len <- 80   # 80 years of data
a_true <- 0.8
b_true <- 0.7
sigma_w_true <- 0.15
sigma_v_true <- 0.30

x_true <- numeric(T_len)
x_true[1] <- a_true / (1 - b_true)  # stationary mean
for (t in 2:T_len) {
    x_true[t] <- a_true + b_true * x_true[t-1] + rnorm(1, sd = sigma_w_true)
}
y_eco <- x_true + rnorm(T_len, sd = sigma_v_true)

# Introduce some missing data (years 30-35 unobserved)
y_eco[30:35] <- NA
```

``` r

# Step 1: Specify the model via SSM_model with custom build_ssm

build_gompertz <- function(parm) {
    a <- parm[1]; b <- parm[2]
    sigma_w <- parm[3]; sigma_v <- parm[4]
    list(
        F = array(b, c(1, 1, 1)),
        Q = array(sigma_w^2, c(1, 1, 1)),
        H = array(1, c(1, 1, 1)),
        R = array(sigma_v^2, c(1, 1, 1)),
        m0 = a / (1 - b),
        P0 = matrix(sigma_w^2 / (1 - b^2), 1, 1)
    )
}

# Conjugate + MH theta update
theta_gompertz <- function(parm, states, Data) {
    x <- states[, 1]
    y_obs <- as.vector(Data$.y)
    T_obs <- length(x)

    # Observation noise: conjugate inv-gamma
    valid <- !is.na(y_obs)
    resid <- y_obs[valid] - x[valid]
    n_obs <- sum(valid)
    sigma_v_sq <- 1 / rgamma(1, 0.001 + n_obs / 2,
                              0.001 + sum(resid^2) / 2)

    # Process noise: conjugate inv-gamma (conditional on a, b)
    dx <- x[2:T_obs] - parm[1] - parm[2] * x[1:(T_obs-1)]
    n_dx <- length(dx)
    sigma_w_sq <- 1 / rgamma(1, 0.001 + n_dx / 2,
                              0.001 + sum(dx^2) / 2)

    # a, b: conjugate normal (conditional on sigma_w)
    X_mat <- cbind(1, x[1:(T_obs-1)])
    y_state <- x[2:T_obs]
    V_beta <- solve(crossprod(X_mat) / sigma_w_sq + diag(0.01, 2))
    m_beta <- V_beta %*% (crossprod(X_mat, y_state) / sigma_w_sq)
    ab_new <- as.vector(m_beta + chol(V_beta) %*% rnorm(2))

    # Enforce stationarity: |b| < 1
    if (abs(ab_new[2]) >= 1) ab_new <- parm[1:2]

    c(ab_new, sqrt(sigma_w_sq), sqrt(sigma_v_sq))
}

mod_gomp <- SSM_model(
    data = y_eco,
    build_ssm = build_gompertz,
    theta_update_fn = theta_gompertz,
    parm.names = c("a", "b", "sigma_w", "sigma_v"),
    state.names = "log_abundance",
    Initial.Values = c(0.5, 0.5, 0.2, 0.3)
)
```

``` r

# Step 2: Fit the model
fit_gomp <- SSM(mod_gomp, Algorithm = "FFBS", Iterations = 5000,
                Thinning = 2, verbose = TRUE, Status = 1000)

# Step 3: Inspect results
print(fit_gomp)
summary(fit_gomp)
```

``` r

# Step 4: Visualize
plot(fit_gomp, type = "states")       # state trajectory + interpolated missing data
plot(fit_gomp, type = "parameters")   # trace plots and densities
plot(fit_gomp, type = "pairs")        # posterior correlations
```

``` r

# Step 5: Convergence diagnostics
demon_gomp <- as.demonoid(fit_gomp)
Rhat(demon_gomp)
ESS(demon_gomp)
```

``` r

# Step 6: Model comparison - Gompertz vs. random walk
mod_rw <- SSM_model(type = "local_level", data = y_eco)
fit_rw <- SSM(mod_rw, Iterations = 5000, Thinning = 2)

demon_rw <- as.demonoid(fit_rw)
loo_gomp <- LOO(demon_gomp)
loo_rw   <- LOO(demon_rw)
loo_compare(loo_gomp, loo_rw)
```

``` r

# Step 7: Prior sensitivity analysis
rb_gomp <- RobustBayes(demon_gomp, type = "power")
summary(rb_gomp)
plot(rb_gomp)
```

This workflow illustrates several key features of the SSM module: the
missing data at years 30-35 are automatically interpolated by the Kalman
smoother, with appropriate uncertainty propagation; the model comparison
via LOO-PSIS reveals whether density dependence (the Gompertz model) is
supported by the data relative to a simple random walk; and the prior
sensitivity analysis identifies which parameters are well-identified by
the data versus driven by the prior. The Gompertz model’s stationary
mean \\a/(1-b) \approx 2.67\\ and the autocorrelation coefficient \\b\\
are typically well-identified for time series of moderate length (\\T
\geq 50\\), while the process and observation noise variances can be
difficult to separate when the signal-to-noise ratio is close to one.

## When Bayesian and frequentist SSM estimates converge

The Bayesian SSM framework implemented through
[`SSM()`](https://robustecologies.github.io/lucifer/reference/SSM.md)
produces posterior distributions over both parameters \\\theta\\ and
latent states \\x\_{0:T}\\. A natural question, particularly for
practitioners trained in classical time series analysis, is under what
conditions these Bayesian estimates coincide with their frequentist
counterparts: maximum likelihood estimates from the Kalman filter, Wald
confidence intervals from the observed information, and likelihood ratio
tests. This section develops the theory of Bayesian-frequentist
convergence for state-space models and demonstrates the bridge using
lucifer’s frequentist infrastructure (`freq_summary`, `wald_test`,
`lr_test`, `profile_likelihood`, `confint_compare`).

### Theoretical conditions for convergence

For standard i.i.d. models, the Bernstein-von Mises theorem guarantees
that the posterior converges to a Gaussian centered at the MLE with
covariance equal to the inverse Fisher information as \\n \to \infty\\.
For state-space models the situation is more subtle, because the
observations are dependent and the latent states introduce a
high-dimensional nuisance parameter. Three conditions must hold for the
Bayesian posterior over \\\theta\\ to converge to the frequentist MLE
asymptotically.

First, the SSM must be **identifiable**: the mapping \\\theta \mapsto
p(y\_{1:T} \mid \theta) = \int p(y\_{1:T}, x\_{0:T} \mid \theta) \\
dx\_{0:T}\\ must be injective up to equivalence classes. For
linear-Gaussian SSMs, identifiability requires that the state dimension
does not exceed the observation dimension multiplied by the lag order,
and that the system \\(F, H)\\ is observable in the control-theoretic
sense. Unobservable states create ridges in the likelihood surface where
distinct \\\theta\\ values produce identical marginal likelihoods, and
the posterior spreads along these ridges regardless of sample size.

Second, the SSM must satisfy **ergodicity** conditions ensuring that the
Kalman filter (or particle filter) forgets its initial condition
exponentially fast. For stable linear systems (all eigenvalues of \\F\\
inside the unit circle), the filter covariance converges to a
steady-state value, and the log-likelihood decomposes into a sum of
weakly dependent prediction-error terms to which a central limit theorem
applies [\[16\]](#ref16). For nonlinear models, the ergodicity
conditions are model-specific and generally require that the state
process is geometrically ergodic.

Third, the **prior** must be smooth and positive in a neighborhood of
the true parameter value. Priors that place zero mass on the true
\\\theta_0\\ (such as a point mass or a prior with support excluding
\\\theta_0\\) will prevent convergence. Diffuse priors, such as
\\\sigma^2 \sim \text{InvGamma}(\epsilon, \epsilon)\\ with small
\\\epsilon\\, satisfy this condition. The conjugate inverse-gamma priors
provided by the SSM builder registry (Section 2) are suitable because
they are everywhere positive on \\(0, \infty)\\.

When all three conditions hold, the marginal posterior \\p(\theta \mid
y\_{1:T})\\ concentrates around the MLE \\\hat{\theta}\_{\text{MLE}}\\
at rate \\T^{-1/2}\\, and the posterior covariance matrix converges to
the inverse of the observed Fisher information evaluated at the MLE. The
filtered state estimates \\E\[x_t \mid y\_{1:t}, \hat{\theta}\]\\ from
the Kalman filter are identical to the posterior mean of \\x_t\\ under a
flat prior for \\\theta\\, because the Kalman filter is the exact
Bayesian filter for the linear-Gaussian case.

### The prediction-error decomposition and the MLE

For linear-Gaussian SSMs, the log-likelihood has a closed-form
expression via the prediction-error decomposition:

\\\ell(\theta) = -\frac{1}{2} \sum\_{t=1}^{T} \left\[ d_y \log 2\pi +
\log \|S_t(\theta)\| + v_t(\theta)^\top S_t(\theta)^{-1} v_t(\theta)
\right\]\\

where \\v_t = y_t - H_t m\_{t\|t-1}\\ is the innovation and \\S_t = H_t
P\_{t\|t-1} H_t^\top + R_t\\ is the innovation covariance, both computed
by the Kalman filter. This is exactly the `loglik` field returned by
`kalman_filter_cpp()`. The MLE \\\hat{\theta}\\ maximizes this function
over \\\theta\\; the Laplace approximation in lucifer
([`LaplaceApproximation()`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md))
finds this maximum and computes the Hessian for standard errors.

The key insight is that FFBS with a flat prior on \\\theta\\ produces
posterior samples whose marginal distribution over \\\theta\\ converges
to a Gaussian centered at the MLE. The
[`as.demonoid()`](https://robustecologies.github.io/lucifer/reference/as.demonoid.md)
conversion provides direct access to the
[`freq_summary()`](https://robustecologies.github.io/lucifer/reference/freq_summary.md),
[`wald_test()`](https://robustecologies.github.io/lucifer/reference/wald_test.md),
and
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
functions from the frequentist bridge.

### Demonstration: local level model

The following example shows exact convergence between Bayesian FFBS and
frequentist Kalman-filter MLE for a local level model. With flat priors,
the posterior mean equals the MLE and the posterior standard deviation
equals the frequentist standard error.

``` r



# Simulate local level with known parameters
set.seed(42)
T_len <- 500
sigma_w_true <- 0.7
sigma_v_true <- 1.0
x_true <- cumsum(rnorm(T_len, 0, sigma_w_true))
y <- x_true + rnorm(T_len, 0, sigma_v_true)

# --- Bayesian estimation via FFBS --- #
mod <- SSM_model(type = "local_level", data = y, times = 1:T_len)
fit_bayes <- SSM(mod, Algorithm = "FFBS", Iterations = 5000,
                 Thinning = 2, verbose = TRUE, Status = 1000)

# Posterior summary
demon <- as.demonoid(fit_bayes)

# --- Frequentist extraction --- #
freq <- freq_summary(demon)

knitr::kable(data.frame(
    Parameter = c("sigma_w", "sigma_v"),
    `Posterior mean` = round(fit_bayes$Summary[, "Mean"], 4),
    `Posterior SD` = round(fit_bayes$Summary[, "SD"], 4),
    `Freq. MLE` = round(freq$coefficients[, "Estimate"], 4),
    `Freq. SE` = round(freq$coefficients[, "Std. Error"], 4),
    check.names = FALSE
), caption = "Bayesian posterior vs. frequentist MLE for the local level model (T = 500).")
```

For \\T = 500\\, the posterior means and MLEs typically agree to three
decimal places. The posterior SDs and frequentist SEs agree to two
decimal places; any residual discrepancy reflects the prior’s finite
(though negligible) influence and the MCMC sampling variability. As
\\T\\ increases, the agreement becomes exact.

### Confidence interval comparison

The
[`confint_compare()`](https://robustecologies.github.io/lucifer/reference/confint_compare.md)
function overlays Bayesian credible intervals (from the posterior
quantiles) with frequentist confidence intervals (from the Hessian-based
normal approximation and, optionally, profile likelihood) on a single
plot. For well-identified parameters with large \\T\\, the intervals are
visually indistinguishable.

``` r

# Compare interval estimates
ci_compare <- confint_compare(demon)
plot(ci_compare)
```

The profile likelihood intervals (obtained via
[`profile_likelihood()`](https://robustecologies.github.io/lucifer/reference/profile_likelihood.md))
are particularly informative for parameters near boundaries. For
example, when the true \\\sigma_w\\ is small, the MLE may be near zero
and the Wald interval can include negative values, while both the
profile likelihood interval and the Bayesian credible interval are
naturally bounded below by zero (since the prior and the likelihood both
assign zero probability to negative variances).

``` r

# Profile likelihood for sigma_w
prof <- profile_likelihood(demon, parm = "sigma_w",
                           range = c(0.3, 1.2), n_points = 30)
plot(prof)
```

### Wald and likelihood ratio tests for model comparison

Classical time series analysis often tests whether specific components
of a structural model are necessary. For example, does a local linear
trend improve on a local level? Is the seasonal component statistically
significant? These questions are answered by the Wald test (testing
whether the additional parameters are zero) or the likelihood ratio test
(comparing the maximized log-likelihoods of nested models).

In the SSM framework, the FFBS log-likelihood from the Kalman
prediction-error decomposition is the exact log-likelihood. The
[`lr_test()`](https://robustecologies.github.io/lucifer/reference/lr_test.md)
function compares two fitted models:

``` r

# Fit two competing models
mod_ll <- SSM_model(type = "local_level", data = y)
mod_llt <- SSM_model(type = "local_linear_trend", data = y)

fit_ll <- SSM(mod_ll, Algorithm = "FFBS", Iterations = 3000)
fit_llt <- SSM(mod_llt, Algorithm = "FFBS", Iterations = 3000)

demon_ll <- as.demonoid(fit_ll)
demon_llt <- as.demonoid(fit_llt)

# Likelihood ratio test (H0: slope variance = 0)
lr <- lr_test(demon_ll, demon_llt)

knitr::kable(data.frame(
    Test = "LR (local level vs. local linear trend)",
    Statistic = round(lr$statistic, 3),
    df = lr$df,
    `p-value` = format.pval(lr$p.value, digits = 3),
    check.names = FALSE
), caption = "Likelihood ratio test for the slope component.")
```

The Wald test provides an alternative that requires only the
unrestricted model:

``` r

# Wald test for sigma_slope = 0
wt <- wald_test(demon_llt, parm = "sigma_slope", null.value = 0)

knitr::kable(data.frame(
    Test = "Wald (H0: sigma_slope = 0)",
    Statistic = round(wt$statistic, 3),
    df = wt$df,
    `p-value` = format.pval(wt$p.value, digits = 3),
    check.names = FALSE
), caption = "Wald test for the slope component.")
```

When the true model has no slope (\\\sigma\_{\text{slope}} = 0\\), both
tests should fail to reject at conventional significance levels, and the
Bayesian model comparison via LOO (Section 9) should prefer the simpler
model by having a lower ELPD difference.

### When Bayesian and frequentist estimates diverge

Three scenarios produce genuine divergence between Bayesian posteriors
and frequentist MLEs, even asymptotically.

**Boundary parameters.** When the true value of a variance component is
zero, the MLE lies on the boundary of the parameter space. The
asymptotic distribution of the MLE is then a mixture of a point mass at
zero and a half-normal, not a full normal [\[16\]](#ref16). The Bayesian
posterior, by contrast, places continuous mass on the entire positive
half-line (assuming an inverse-gamma or half-Cauchy prior), and its mean
is strictly positive even when the MLE is zero. This is perhaps the most
common divergence in practice: the Bayesian estimate of a small variance
is biased upward relative to the MLE because the prior prevents it from
reaching zero.

``` r

# Boundary case: true sigma_slope = 0
set.seed(123)
T_b <- 200
x_b <- cumsum(rnorm(T_b, 0, 0.5))  # no slope
y_b <- x_b + rnorm(T_b, 0, 1)

mod_b <- SSM_model(type = "local_linear_trend", data = y_b)
fit_b <- SSM(mod_b, Algorithm = "FFBS", Iterations = 3000)

knitr::kable(data.frame(
    Parameter = c("sigma_level", "sigma_slope", "sigma_obs"),
    `Posterior mean` = round(fit_b$Summary[, "Mean"], 4),
    `True value` = c(0.5, 0, 1.0),
    Note = c("", "Bayesian > 0 (prior floor)", ""),
    check.names = FALSE
), caption = "Boundary parameter: sigma_slope = 0.")
```

**Short time series.** For small \\T\\ (say \\T \< 30\\), the prior
contributes meaningfully to the posterior. The posterior mean is a
compromise between the MLE and the prior mean, weighted by the relative
precision. The frequentist standard errors from the Hessian may be
unreliable because the asymptotic normal approximation has not yet
engaged. In this regime, the Bayesian analysis with informative priors
is arguably more appropriate than the frequentist analysis, because it
regularizes the estimates and produces calibrated uncertainty intervals
even when the likelihood surface is poorly conditioned.

**Misspecified models.** When the assumed SSM structure is wrong (e.g.,
a local level model fitted to data with a deterministic trend), the
posterior concentrates on the pseudo-true value \\\theta^\* =
\arg\min\_\theta \text{KL}(p_0 \\ p\_\theta)\\, which is the same as the
pseudo-MLE. However, the posterior width reflects the inverse Fisher
information of the misspecified model, which is generally narrower than
the sandwich standard error that accounts for misspecification. In this
case, frequentist sandwich intervals are wider and have correct
coverage, while the Bayesian credible intervals are anti-conservative.
The
[`freq_summary()`](https://robustecologies.github.io/lucifer/reference/freq_summary.md)
function with `sandwich = TRUE` provides the robust standard errors for
comparison.

``` r

# Misspecified model: local level fitted to trend data
set.seed(456)
y_trend <- 0.1 * (1:100) + rnorm(100, 0, 1)

mod_mis <- SSM_model(type = "local_level", data = y_trend)
fit_mis <- SSM(mod_mis, Algorithm = "FFBS", Iterations = 3000)
demon_mis <- as.demonoid(fit_mis)

freq_hessian <- freq_summary(demon_mis)
# Bayesian CI will be narrower than sandwich CI
# because it ignores the misspecification
```

### Coverage simulation

The ultimate test of Bayesian-frequentist agreement is coverage: do 95%
posterior credible intervals contain the true parameter value 95% of the
time across repeated data sets? The
[`coverage_sim()`](https://robustecologies.github.io/lucifer/reference/coverage_sim.md)
function automates this by repeatedly simulating data from the true
model, fitting the Bayesian model, and checking whether the true
parameter falls within the credible interval.

``` r

# Coverage simulation for local level model
# (computationally intensive; use small n_sim for illustration)
cov_result <- coverage_sim(
    Model = function(parm, Data) {
        # Local level log-likelihood via Kalman filter
        build <- .ssm_builder_local_level(Data$y, 1:length(Data$y))
        ssm_mats <- build$build_ssm(parm)
        kf <- kalman_filter_cpp(matrix(Data$y, ncol = 1),
            ssm_mats$F, ssm_mats$Q, ssm_mats$H, ssm_mats$R,
            ssm_mats$m0, ssm_mats$P0)
        LP <- kf$loglik + sum(dnorm(parm, 0, 10, log = TRUE))
        list(LP = LP, Dev = -2 * kf$loglik, Monitor = LP,
             yhat = rnorm(length(Data$y)), parm = parm)
    },
    Data = list(y = rnorm(100), N = 100,
                parm.names = c("sigma_w", "sigma_v"),
                mon.names = "LP"),
    true_parm = c(0.7, 1.0),
    n_sim = 50,
    Iterations = 1000
)

knitr::kable(data.frame(
    Parameter = c("sigma_w", "sigma_v"),
    `Nominal` = c("95%", "95%"),
    `Observed coverage` = paste0(round(cov_result$coverage * 100, 1), "%"),
    check.names = FALSE
), caption = "Empirical coverage of 95% Bayesian credible intervals.")
```

For \\T = 100\\ with diffuse priors, the coverage should be close to the
nominal 95% for both parameters. Deviations indicate either prior
influence (when coverage exceeds 95%, the intervals are conservative) or
boundary effects (when coverage falls below 95%, the posterior
underestimates uncertainty near zero).

## References

**\[1\]** Carter, C.K. and Kohn, R. (1994). *On Gibbs sampling for state
space models*. Biometrika, 81(3), 541-553. [DOI:
10.1093/biomet/81.3.541](https://doi.org/10.1093/biomet/81.3.541)

**\[2\]** Fruhwirth-Schnatter, S. (1994). *Data augmentation and dynamic
linear models*. Journal of Time Series Analysis, 15(2), 183-202. [DOI:
10.1111/j.1467-9892.1994.tb00184.x](https://doi.org/10.1111/j.1467-9892.1994.tb00184.x)

**\[3\]** Andrieu, C., Doucet, A., and Holenstein, R. (2010). *Particle
Markov chain Monte Carlo methods*. Journal of the Royal Statistical
Society: Series B, 72(3), 269-342. [DOI:
10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)

**\[4\]** Lindsten, F., Jordan, M.I., and Schon, T.B. (2014). *Particle
Gibbs with ancestor sampling*. Journal of Machine Learning Research,
15(1), 2145-2184.

**\[5\]** Chopin, N., Jacob, P.E., and Papaspiliopoulos, O. (2013).
*SMC\\^2\\: an efficient algorithm for sequential analysis of state
space models*. Journal of the Royal Statistical Society: Series B,
75(3), 397-426. [DOI:
10.1111/j.1467-9868.2012.01046.x](https://doi.org/10.1111/j.1467-9868.2012.01046.x)

**\[6\]** Kim, S., Shephard, N., and Chib, S. (1998). *Stochastic
volatility: likelihood inference and comparison with ARCH models*.
Review of Economic Studies, 65(3), 361-393. [DOI:
10.1111/1467-937X.00050](https://doi.org/10.1111/1467-937X.00050)

**\[7\]** Rauch, H.E., Tung, F., and Striebel, C.T. (1965). *Maximum
likelihood estimates of linear dynamic systems*. AIAA Journal, 3(8),
1445-1450. [DOI: 10.2514/3.3166](https://doi.org/10.2514/3.3166)

**\[8\]** Harvey, A.C. (1989). *Forecasting, structural time series
models and the Kalman filter*. Cambridge University Press. ISBN
978-0-521-40573-7.

**\[9\]** Julier, S.J. and Uhlmann, J.K. (1997). *New extension of the
Kalman filter to nonlinear systems*. In Signal Processing, Sensor
Fusion, and Target Recognition VI, SPIE, vol. 3068, 182-193. [DOI:
10.1117/12.280797](https://doi.org/10.1117/12.280797)

**\[10\]** Wan, E.A. and van der Merwe, R. (2000). *The unscented Kalman
filter for nonlinear estimation*. In Proceedings of the IEEE Adaptive
Systems for Signal Processing, Communications, and Control Symposium
(AS-SPCC), 153-158. [DOI:
10.1109/ASSPCC.2000.882463](https://doi.org/10.1109/ASSPCC.2000.882463)

**\[11\]** Evensen, G. (1994). *Sequential data assimilation with a
nonlinear quasi-geostrophic model using Monte Carlo methods to forecast
error statistics*. Journal of Geophysical Research, 99(C5), 10143-10162.
[DOI: 10.1029/94JC00572](https://doi.org/10.1029/94JC00572)

**\[12\]** Doucet, A., de Freitas, N., Murphy, K., and Russell, S.
(2000). *Rao-Blackwellised particle filtering for dynamic Bayesian
networks*. In Proceedings of the 16th Conference on Uncertainty in
Artificial Intelligence (UAI), 176-183.

**\[13\]** Hamilton, J.D. (1989). *A new approach to the economic
analysis of nonstationary time series and the business cycle*.
Econometrica, 57(2), 357-384. [DOI:
10.2307/1912559](https://doi.org/10.2307/1912559)

**\[14\]** Kim, C.-J. (1994). *Dynamic linear models with
Markov-switching*. Journal of Econometrics, 60(1-2), 1-22. [DOI:
10.1016/0304-4076(94)90036-1](https://doi.org/10.1016/0304-4076(94)90036-1)

**\[15\]** Kim, C.-J. and Nelson, C.R. (1999). *State-space models with
regime switching*. MIT Press. ISBN 978-0-262-11236-4.

**\[16\]** Durbin, J. and Koopman, S.J. (2012). *Time series analysis by
state space methods*. Second edition, Oxford University Press. [DOI:
10.1093/acprof:oso/9780199641178.001.0001](https://doi.org/10.1093/acprof:oso/9780199641178.001.0001)

**\[17\]** Kalman, R.E. (1960). *A new approach to linear filtering and
prediction problems*. Journal of Basic Engineering, 82(1), 35-45. [DOI:
10.1115/1.3662552](https://doi.org/10.1115/1.3662552)

**\[18\]** Vehtari, A., Gelman, A., and Gabry, J. (2017). *Practical
Bayesian model evaluation using leave-one-out cross-validation and
WAIC*. Statistics and Computing, 27(5), 1413-1432. [DOI:
10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

**\[19\]** Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and
Burkner, P.-C. (2021). *Rank-normalization, folding, and localization:
an improved R-hat for assessing convergence of MCMC (with discussion)*.
Bayesian Analysis, 16(2), 667-718. [DOI:
10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

**\[20\]** Schon, T.B., Gustafsson, F., and Nordlund, P.-J. (2005).
*Marginalized particle filters for mixed linear/nonlinear state-space
models*. IEEE Transactions on Signal Processing, 53(7), 2168-2177. [DOI:
10.1109/TSP.2005.849151](https://doi.org/10.1109/TSP.2005.849151)

**\[21\]** Del Moral, P., Doucet, A., and Jasra, A. (2006). *Sequential
Monte Carlo samplers*. Journal of the Royal Statistical Society: Series
B, 68(3), 411-436. [DOI:
10.1111/j.1467-9868.2006.00553.x](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

**\[22\]** Doucet, A. and Johansen, A.M. (2011). *A tutorial on particle
filtering and smoothing: fifteen years later*. In Handbook of Nonlinear
Filtering, Oxford University Press, 656-704.
