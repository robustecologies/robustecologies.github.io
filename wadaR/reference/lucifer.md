# lucifer

The `lucifer` function is the main function of the lucifer package.
Given data, a model specification, and initial values, `lucifer`
maximizes the logarithm of the unnormalized joint posterior density with
MCMC and provides samples of the marginal posterior distributions,
deviance, and other monitored variables. `lucifer` supports multi-chain
MCMC directly via the `Chains`, `CPUs`, and `Type` arguments. When
`Chains > 1`, multiple chains are run (sequentially or in parallel
depending on `CPUs`), combined via
[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md),
and the Gelman-Rubin diagnostic is automatically computed. The
recommended setting for production inference is `Chains=4`.

After sampling, burn-in removal splits the output into two posterior
matrices: `Posterior1` retains all thinned samples for diagnostics and
trace plots, while `Posterior2` contains only the post-burn-in
(stationary) samples used for inference. By default, the burn-in
boundary is estimated automatically via the BMK stationarity diagnostic;
the `BurnIn` argument allows the user to override this with a fixed
number of thinned samples to discard. The `Summary2`, `DIC2`, and `LML`
fields, as well as
[`predict.demonoid`](https://robustecologies.github.io/lucifer/reference/predict.demonoid.md),
all operate on `Posterior2`.

## Usage

``` r
lucifer(
  Model,
  Data,
  Initial.Values,
  Covar = NULL,
  Iterations = 10000,
  Status = 100,
  Thinning = 10,
  Algorithm = "MWG",
  Specs = list(B = NULL),
  Debug = list(DB.chol = FALSE, DB.eigen = FALSE, DB.MCSE = FALSE, DB.Model = TRUE),
  LogFile = "",
  Chains = 3L,
  CPUs = 3L,
  Type = NULL,
  Grad = NULL,
  BurnIn = NULL,
  ...
)
```

## Arguments

- Model:

  This required argument receives the model from a user-defined function
  that must be named Model. The user-defined function is where the model
  is specified. `lucifer` passes two arguments to the model function,
  `parms` and `Data`, and receives five arguments from the model
  function: `LP` (the logarithm of the unnormalized joint posterior),
  `Dev` (the deviance), `Monitor` (the monitored variables), `yhat` (the
  variables for posterior predictive checks), and `parm`, the vector of
  parameters, which may be constrained in the model function. More
  information on the Model specification function may be found in the
  "lucifer tutorial" vignette and the
  [`is.model`](https://robustecologies.github.io/lucifer/reference/is.model.md)
  function. Many examples of model specification functions may be found
  in the "Examples" vignette.

- Data:

  This required argument accepts a list of data. The list of data must
  contain `mon.names` which contains monitored variable names, and must
  contain `parm.names` which contains parameter names. The
  [`as.parm.names`](https://robustecologies.github.io/lucifer/reference/as.parm.names.md)
  function may be helpful for preparing the data, and the
  [`is.data`](https://robustecologies.github.io/lucifer/reference/is.data.md)
  function may be helpful for checking data.

- Initial.Values:

  This argument requires a vector of initial values equal in length to
  the number of parameters. Each initial value will be the starting
  point for an adaptive chain or a non-adaptive Markov chain of a
  parameter. Parameters are assumed to be continuous, unless specified
  to be discrete (see `dparm`), which is not accepted by all algorithms
  (see
  [`dcrmrf`](https://robustecologies.github.io/lucifer/reference/dist.ContinuousRelaxation.md)
  for an alternative). If all initial values are set to zero, then
  Lucifer will attempt to optimize the initial values with the
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
  function. After Lucifer finishes updating, it may be desired to
  continue updating from where it left off. To continue, this argument
  should receive the last iteration of the previous update, for example
  `Initial.Values=as.initial.values(Fit)`. Initial values may be
  generated randomly with the
  [`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md)
  function. When `Chains > 1` and `Initial.Values` is a vector, chain 1
  uses the original values and chains 2..n receive Gaussian jitter. When
  `Initial.Values` is a matrix with `nrow == Chains`, each row is used
  directly.

- Covar:

  This argument defaults to `NULL`, but may otherwise accept a \\K
  \times K\\ proposal covariance matrix (where \\K\\ is the number of
  dimensions or parameters), a variance vector, or a list of covariance
  matrices (for blockwise sampling in some algorithms). When the model
  is updated for the first time and prior variance or covariance is
  unknown, then `Covar=NULL` should be used. Some algorithms require
  covariance, some only require variance, and some require neither.
  Lucifer automatically converts the user input to the required form.
  Once Lucifer has finished updating, it may be desired to continue
  updating where it left off, in which case the proposal covariance
  matrix from the last run can be input into the next run. The
  covariance matrix may also be input from the
  [`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md)
  function, if used.

- Iterations:

  This required argument accepts integers larger than 10, and determines
  the number of iterations that Lucifer will update the parameters while
  searching for target distributions. The required amount of computer
  memory will increase with `Iterations`. If computer memory is
  exceeded, then all will be lost. The
  [`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md)
  function can be used later to combine multiple updates.

- Status:

  This argument accepts an integer between 1 and the number of
  iterations, and indicates how often, in iterations, the user would
  like the status printed to the screen or log file.

- Thinning:

  This argument accepts integers between 1 and the number of iterations,
  and indicates that every nth iteration will be retained, while the
  other iterations are discarded. Thinning is performed to reduce
  autocorrelation and the number of marginal posterior samples.

- Algorithm:

  This argument accepts the abbreviated name of the MCMC algorithm,
  which must appear in quotes. A list of MCMC algorithms appears in the
  details section, and the abbreviated name is in parenthesis.

- Specs:

  This argument defaults to `NULL`, and accepts a list of specifications
  for the MCMC algorithm declared in the `Algorithm` argument. The
  specifications associated with each algorithm are described in the
  details section. For the Gibbs sampler, when `Specs` is `NULL`,
  automatic full conditionals are constructed internally using a mixed
  strategy: continuous parameters are sampled via univariate slice
  sampling (Neal 2003), and discrete parameters are sampled via
  exhaustive enumeration over their finite support. Discrete parameters
  are declared through `Data$dparm` (integer vector of parameter
  indices) and optionally `Data$dsupport` (list of support vectors,
  parallel to `dparm`; defaults to `c(0, 1)` for Bernoulli indicators).
  In auto-FC mode the acceptance rate is 100% and no MWG updates are
  performed. To use manual full conditionals, supply
  `Specs = list(FC = myFC)` as before. See
  [`vignette("gibbs", package = "lucifer")`](https://robustecologies.github.io/lucifer/articles/gibbs.md)
  for details.

- Debug:

  This argument accepts a list of logical scalars that control whether
  or not errors or warnings are reported due to a `try` function or
  non-finite values. List components include `DB.chol` regarding `chol`,
  `DB.eigen` regarding `eigen`, `DB.MCSE` regarding
  [`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md),
  and `DB.Model` regarding the Model specification function.

- LogFile:

  This argument is used to specify a log file name in quotes in the
  working directory as a destination, rather than the console, for the
  output messages.

- Chains:

  Number of MCMC chains to run. Defaults to 3. When `Chains > 1`,
  multiple chains are executed (in parallel if `CPUs > 1`), combined via
  [`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md),
  and the Gelman-Rubin diagnostic is automatically computed. When
  `Initial.Values` is a vector and `Chains > 1`, chain 1 uses the
  original values and chains 2..n receive Gaussian jitter. When
  `Initial.Values` is a matrix with `nrow == Chains`, each row is used
  directly. Set `Chains=1` for single-chain runs. Recommended:
  `Chains=4` for production inference.

- CPUs:

  Number of CPUs for parallel execution. Defaults to 3. Automatically
  capped at the number of detected cores.

- Type:

  Parallelization type. Defaults to `NULL`, which auto-detects: `"FORK"`
  on Unix/macOS (faster, shared memory) and `"PSOCK"` on Windows. Can be
  set explicitly to `"PSOCK"`, `"FORK"`, or `"MPI"`.

- Grad:

  Controls gradient computation for gradient-based MCMC algorithms
  (AHMC, HMC, HMCDA, MALA, NUTS, Refractive, SGLD, THMC). Three modes
  are supported: `NULL` (default) uses the built-in forward
  finite-difference approximation via `partial_cpp`, which costs LIV+1
  model evaluations per gradient; a user-supplied function with
  signature `function(parm, Data)` returning a numeric vector of length
  LIV, which is validated at startup against the initial values; or the
  string `"numDeriv"`, which wraps
  [`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html) with
  Richardson extrapolation for higher accuracy (requires the numDeriv
  package). Ignored for non-gradient algorithms.

- BurnIn:

  Controls burn-in (warmup) removal from the posterior samples. Three
  modes are available. When `NULL` (the default), burn-in is estimated
  automatically using the BMK stationarity diagnostic (Hellinger
  distance \> 0.5 on 10 batches); if no stationary region is detected,
  the first 50% of samples are discarded as a conservative fallback.
  When a non-negative integer, that many thinned samples are discarded
  from the beginning of each chain, bypassing the automatic diagnostic
  entirely. For example, `BurnIn = 200` discards the first 200 thinned
  samples. The burn-in affects `Posterior2`, `Summary2`, and `DIC2` in
  the returned object; `Posterior1` always contains all thinned samples.
  When `Chains > 1`, the burn-in is applied after chains are
  concatenated by
  [`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md).
  Post-hoc burn-in adjustment is available via
  [`deburn`](https://robustecologies.github.io/lucifer/reference/deburn.md).

- ...:

  Additional arguments are unused.

## Value

`lucifer` returns an object of class `demonoid`. When `Chains > 1`, the
returned object contains combined posteriors from all chains (via
[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md))
plus per-chain results and the Gelman-Rubin diagnostic. Each object of
class `demonoid` is a list with the following components:

- Acceptance.Rate:

  The acceptance rate of the MCMC algorithm, indicating the percentage
  of iterations in which the proposals were accepted.

- Algorithm:

  The specific algorithm used.

- Call:

  The matched call of `lucifer`.

- Chains:

  When `Chains > 1`, a list of per-chain `demonoid` objects. `NULL` for
  single-chain runs.

- Chains.n:

  Integer number of chains (1 for single-chain runs).

- CPUs:

  Integer number of CPUs used.

- Covar:

  The \\K \times K\\ proposal covariance matrix, variance vector, or
  list of covariance matrices.

- CovarDHis:

  An \\N \times K\\ matrix storing the diagonal of the proposal
  covariance matrix of each adaptation.

- Deviance:

  A vector of the deviance of the model, with a length equal to the
  number of thinned samples that were retained.

- DIC1:

  A vector of three values: Dbar, pD, and DIC, calculated over all
  retained samples.

- DIC2:

  Identical to `DIC1` but calculated over only the samples considered
  stationary by `BMK.Diagnostic`.

- Gelman:

  When `Chains >= 2`, contains the output of
  [`Gelman.Diagnostic`](https://robustecologies.github.io/lucifer/reference/Gelman.Diagnostic.md),
  including PSRF per parameter and MPSRF. `NULL` for single-chain runs.

- Initial.Values:

  The vector of initial values, which may have been optimized.

- Iterations:

  The number of iterations for updating.

- LML:

  An approximation of the logarithm of the marginal likelihood of the
  data (see
  [`LML`](https://robustecologies.github.io/lucifer/reference/LML.md)).

- Minutes:

  The number of minutes that `lucifer` was running.

- Model:

  The model specification function.

- Monitor:

  A vector or matrix of monitored variables.

- Parameters:

  The number of parameters.

- Posterior1:

  A matrix of marginal posterior distributions composed of all thinned
  samples.

- Posterior2:

  A matrix equal to `Posterior1`, except that only stationary samples
  are included.

- Rec.BurnIn.Thinned:

  The recommended burn-in for the thinned samples.

- Rec.BurnIn.UnThinned:

  The recommended burn-in for all samples.

- Rec.Thinning:

  The recommended value for the `Thinning` argument.

- Specs:

  An optional list of algorithm specifications.

- Status:

  The value in the `Status` argument.

- Summary1:

  A matrix summarizing the marginal posterior distributions over all
  samples in `Posterior1`, including mean, SD, MCSE, ESS, and quantiles
  (2.5%, 50%, 97.5%).

- Summary2:

  Identical to `Summary1` but calculated only on the stationary samples
  in `Posterior2`.

- Thinned.Samples:

  The number of thinned samples retained.

- Thinning:

  The value of the `Thinning` argument.

## Details

lucifer

`lucifer` offers numerous MCMC algorithms for numerical approximation in
Bayesian inference, organized below by methodological family.

### Available algorithms

`lucifer` dispatches more than eighty MCMC samplers, organized below by
methodological family. Each algorithm is identified by an abbreviation
accepted by the `Algorithm` argument; the full description, default
specifications, and references for every algorithm are documented in
[`vignette("mcmc-algorithms", package = "lucifer")`](https://robustecologies.github.io/lucifer/articles/mcmc-algorithms.md).

- **Random-walk and adaptive Metropolis:** Random-Walk Metropolis (RWM),
  Adaptive Metropolis (AM), Adaptive-Mixture Metropolis (AMM), Robust
  Adaptive Metropolis (RAM), Delayed Rejection Metropolis (DRM), Delayed
  Rejection Adaptive Metropolis (DRAM), Delayed Acceptance MCMC (DA),
  Independence Metropolis (IM), Multiple-Try Metropolis (MTM), Random
  Dive Metropolis-Hastings (RDMH), Interchain Adaptation (INCA),
  Pseudo-Marginal MCMC (PMCMC).

- **Componentwise and Gibbs samplers:** Gibbs Sampler (Gibbs),
  Metropolis-within-Gibbs (MWG), Adaptive Metropolis-within-Gibbs
  (AMWG), Adaptive Directional Metropolis-within-Gibbs (ADMG),
  Griddy-Gibbs (GG), Adaptive Griddy-Gibbs (AGG), Sequential
  Metropolis-within-Gibbs (SMWG), Sequential Adaptive
  Metropolis-within-Gibbs (SAMWG), Updating Sequential
  Metropolis-within-Gibbs (USMWG), Updating Sequential Adaptive
  Metropolis-within-Gibbs (USAMWG), Zanella locally balanced proposal
  (Zanella), Polya-Gamma data augmentation (PG).

- **Hamiltonian Monte Carlo:** Hamiltonian Monte Carlo (HMC), Adaptive
  Hamiltonian Monte Carlo (AHMC), No-U-Turn Sampler (NUTS), Hamiltonian
  Monte Carlo with Dual-Averaging (HMCDA), Tempered Hamiltonian Monte
  Carlo (THMC), Microcanonical Hamiltonian Monte Carlo (MCHMC), MEADS,
  autoMALA, Metropolis Adjusted Langevin Trajectories (MALT),
  Apogee-to-Apogee Path Sampler (AAPS), Gibbs Self-Tuning NUTS (GIST).

- **Langevin and stochastic-gradient samplers:** Metropolis-Adjusted
  Langevin Algorithm (MALA), Stochastic Gradient Langevin Dynamics
  (SGLD), Stochastic Gradient Hamiltonian Monte Carlo (SGHMC), Barker
  proposal (Barker), Non-Reversible Overdamped Langevin (NROLangevin).

- **Slice samplers:** Slice Sampler (Slice), Automated Factor Slice
  Sampler (AFSS), Elliptical Slice Sampler (ESS), Oblique Hyperrectangle
  Slice Sampler (OHSS), Reflective Slice Sampler (RSS), Univariate
  Eigenvector Slice Sampler (UESS), Quantile Slice Sampler (QSS), Latent
  Slice Sampler (LSS), Gibbsian Polar Slice Sampler (GPSS).

- **Ensemble and population samplers:** Affine-Invariant Ensemble
  Sampler (AIES), Differential Evolution Markov Chain (DEMC),
  Differential Evolution Adaptive Metropolis (DREAM), Ensemble Slice
  Sampler (Zeus), t-walk (twalk).

- **Quadratic-Monte Carlo ensemble:** Quadratic Monte Carlo (QMC),
  Third-Order Monte Carlo (QMC3), Nth-Order Monte Carlo (QMCN), Directed
  Quadratic Monte Carlo (DQMC), Quadratic Simplex Monte Carlo (SQMC),
  Modified Affine Monte Carlo (MAMC), Affine Simplex Monte Carlo (SAMC),
  Walk Monte Carlo (WMC).

- **Hit-and-run and directional samplers:** Hit-And-Run Metropolis
  (HARM), Componentwise Hit-And-Run Metropolis (CHARM), Refractive
  Sampler (Refractive).

- **Multimodal and tempering:** Metropolis-Coupled Markov Chain Monte
  Carlo (MCMCMC), Reversible-Jump (RJ), Non-Reversible Parallel
  Tempering (NRPT), Simulated Tempering (SimTemp), Non-Reversible
  Simulated Tempering (NRST), Wang-Landau (WL).

- **Piecewise-deterministic Markov processes:** Zig-Zag Sampler
  (ZigZag), Bouncy Particle Sampler (BPS), Boomerang Sampler
  (Boomerang), Randomized HMC (RHMC).

- **Riemannian and geometric MCMC:** Riemannian Manifold HMC (RMHMC),
  Lagrangian Monte Carlo (LMC), Magnetic HMC (MHMC), Relativistic Monte
  Carlo (Relativistic).

- **Constraint-handling samplers:** Projected Langevin (ProjLang),
  proximal MCMC (ProxMCMC).

- **Flow-enhanced samplers:** Neural Transport MCMC (NeuTra), adaptive
  normalizing-flow MCMC (flowMC). Both require the torch package.

- **Infinite-dimensional targets:** Preconditioned Crank-Nicolson (pCN).

### Algorithm specifications

The `Specs` argument is a named list whose contents depend on the chosen
algorithm. When `Specs = NULL` (the default), an algorithm-specific
defaults function is invoked, except for the eleven algorithms in which
at least one specification has no sensible default value (`AGG`, `GG`,
`IM`, `RJ`, `RSS`, `SAMWG`, `SMWG`, `USAMWG`, `USMWG`, `SGLD`, and
`SGHMC`). The recognized specification elements are:

- `A`:

  Length of the warm-up window in AAPS, AFSS, Barker, GIST, HMCDA,
  MCHMC, MEADS, NUTS, OHSS, and UESS. In MALA it is the maximum
  acceptable Euclidean norm of the adaptive drift.

- `Adaptive`:

  Iteration at which adaptation begins, used in AM, AMM, DRAM, INCA, and
  Refractive.

- `a`:

  Quadratic step amplitude in QMC, QMC3, QMCN, DQMC, SQMC, MAMC, SAMC,
  WMC.

- `alpha`:

  Friction parameter in MCHMC and SGHMC.

- `alpha_gist`:

  Sharpness of the softmax used to draw a step size from the GIST
  candidate grid.

- `alpha.star`:

  Target acceptance rate in MALA and RAM, optional in CHARM and HARM.
  Recommended values: 0.234 for multivariate proposals, 0.44 for
  componentwise proposals, 0.574 for MALA.

- `at`, `aw`:

  Traverse and walk moves in twalk (`at = 6`, `aw = 1.5` recommended).

- `B`:

  List of blocked parameters, optional in AFSS, AMM, AMWG, ESS, HARM,
  MWG, PG, PMCMC, RAM, RWM, Slice, UESS, Zanella; see
  [`Blocks`](https://robustecologies.github.io/lucifer/reference/Blocks.md).

- `b`:

  Drift bias in DQMC.

- `base_algorithm`:

  Sampler used by Delayed Acceptance MCMC after the surrogate-stage
  check.

- `balancing_fn`:

  Locally balanced function in Zanella (`"sqrt"`, `"min"`, or
  `"barker"`).

- `Begin`:

  Time-period at which to begin updating in USAMWG and USMWG.

- `beta`:

  Scale parameter for AIES (default 2) or autoregressive parameter for
  pCN.

- `bin.n`, `bin.p`:

  Size and probability parameters for the binomial prior on model size
  in RJ.

- `Bounds`:

  Vector of length two with the lower and upper boundary of the slice in
  the Slice algorithm.

- `c`:

  Speed-of-light cap in Relativistic Monte Carlo.

- `constraint`:

  Constraint type in Projected Langevin (`"box"` or `"simplex"`).

- `CR`:

  Crossover probability in DREAM.

- `b_star`:

  Snooker noise variance in DREAM.

- `delta`:

  Target acceptance rate in AAPS, Barker, GIST, HMCDA, MEADS, NUTS
  (recommended values 0.65, 0.40, 0.6, 0.65, 0.6 and 0.6 respectively).
  In MALA it is a constant in the bounded drift function.

- `Dist`:

  Proposal distribution in RAM (`"t"` or `"N"`).

- `dparm`:

  Vector of integers indicating discrete parameters, for use with AGG or
  GG.

- `Dyn`:

  \\T \times K\\ matrix of dynamic parameters used by SAMWG, SMWG,
  USAMWG, and USMWG.

- `epsilon`:

  Step-size used in AAPS, AHMC, autoMALA (`epsilon_max`), Barker, GIST,
  HMC, HMCDA, MALA, MALT, MCHMC, MEADS, NRST, NROLangevin, NUTS,
  ProjLang, ProxMCMC, Relativistic, RHMC, RMHMC, SGHMC, SGLD, SimTemp,
  THMC, and Zig-Zag (`T_max`).

- `epsilon_max`:

  Upper bound on the proposed step size searched by autoMALA.

- `excess_rate`:

  Refresh-rate boost for piecewise deterministic samplers (BPS,
  Boomerang, ZigZag).

- `f_init`, `f_min`, `flatness`:

  Adaptation schedule for the Wang-Landau histogram.

- `FC`:

  Used in Gibbs; accepts either a function receiving
  `(parameters, data)` and returning updated parameters, or the literal
  string `"auto"` to construct full conditionals automatically (slice
  sampling for continuous parameters and exhaustive enumeration for
  discrete parameters declared via `Data$dparm`).

- `field_strength`:

  Magnetic-field magnitude in MHMC.

- `file`:

  Quoted name of a .csv data file for SGLD and SGHMC.

- `Fit`:

  Object of class `demonoid` in USAMWG and USMWG.

- `fixedpoint_maxiter`, `fixedpoint_tol`:

  Convergence controls for the fixed-point iteration of the generalized
  leapfrog used by RMHMC.

- `Gaussian`:

  Logical indicating whether the QMC family should sample with Gaussian
  proposals (`TRUE`) or default quadratic proposals (`FALSE`).

- `gamma`:

  Step size in DEMC, friction in MALT, or decay of adaptation in MALA
  and RAM.

- `Grid`:

  Vector or list of vectors of evenly-spaced grid points for AGG or GG.

- `grid_range`:

  Half-width (log2) of the GIST candidate grid of step sizes.

- `K`:

  Number of proposals in MTM, or number of step-size candidates in GIST.

- `L`:

  Number of leapfrog steps in AHMC, HMC, LMC, MALT, MCHMC, MEADS, MHMC,
  Relativistic, RMHMC, SGHMC, and THMC.

- `lambda`:

  Scalar trajectory length in HMCDA, ladder spacing in MCMCMC and NRPT,
  geometric tempering exponent in SimTemp and NRST.

- `lambda_reg`:

  Regularization strength in ProxMCMC.

- `Lmax`:

  Maximum for `L` in AAPS, GIST, HMCDA, and NUTS.

- `lower`, `upper`:

  Box-constraint vectors in Projected Langevin.

- `m`:

  Context-dependent: mass matrix, maximum steps, or steps per iteration
  in AFSS, AHMC, HMC, Refractive, RSS, Slice, THMC, and UESS.

- `max_search`:

  Maximum number of doubling/halving steps in autoMALA's step-size
  search.

- `max_steps`:

  Maximum number of stepping-out iterations in GPSS.

- `mean_L`:

  Mean of the geometric distribution of leapfrog lengths in Randomized
  HMC.

- `method`:

  Variant of the Zig-Zag sampler (`"zigzag"` or `"local"`).

- `metric`:

  Metric tensor used by RMHMC, LMC, and MHMC (`"identity"`, `"hessian"`,
  or `"softabs"`).

- `mu`:

  Vector equal in length to initial values, used as the mean of the
  proposal distribution in IM and the mode in Zeus.

- `MWG`:

  Vector of parameter positions receiving Metropolis-within-Gibbs
  updates in Gibbs (ignored when `FC = "auto"`).

- `n`:

  Number of previous iterations in ADMG, AFSS, AMM, AMWG, OHSS, RAM, and
  UESS.

- `n1`:

  Subset size for each set of points in twalk (`n1 = 4` recommended).

- `n_apogees`:

  Number of apogees explored per iteration in AAPS.

- `n_bins`:

  Number of energy bins in Wang-Landau.

- `n_flow_epochs`:

  Number of training epochs for the normalizing flow used by NeuTra and
  flowMC.

- `n_local`:

  Number of local Metropolis steps per round in flowMC.

- `n_pairs`:

  Number of difference pairs sampled per DREAM update.

- `n_rounds`:

  Number of train-then-sample rounds in flowMC.

- `n_steps`:

  Number of within-temperature steps between NRPT swap moves.

- `n_temps`:

  Number of temperatures in NRPT, SimTemp, and NRST.

- `n_warmup`:

  Length of the warm-up phase used to train the NeuTra normalizing flow.

- `nGuidePoints`:

  Number of guide points in SQMC, SAMC, and WMC.

- `nOrder`:

  Polynomial order in QMCN.

- `Nc`:

  Number of walkers in AIES, DEMC, DREAM, Zeus, QMC, QMC3, QMCN, DQMC,
  SQMC, MAMC, SAMC, and WMC, or number of columns of big data in SGLD
  and SGHMC.

- `Nr`:

  Number of rows of big data in SGLD and SGHMC.

- `parm.p`:

  Vector of probabilities for parameter selection in RJ.

- `Periodicity`:

  Adaptation frequency (in iterations) for AHMC, AM, AMM, AMWG, DRAM,
  INCA, NRPT, SAMWG, and USAMWG.

- `prior.mean`, `prior.var`:

  Gaussian prior moments for Polya-Gamma data augmentation.

- `prox`:

  Proximity operator in ProxMCMC (`"l1"`, `"tv"`, or `"box"`).

- `r`:

  Ratio between r1 and r2 in the Refractive algorithm.

- `refresh_rate`:

  Velocity refreshment rate in BPS, Boomerang, and ZigZag.

- `s_init`:

  Initial latent slice width in LSS.

- `selectable`:

  Vector of indicators for variable selection eligibility in RJ.

- `selected`:

  Vector of indicators for initial variable selection state in RJ.

- `SIV`:

  Secondary initial values for twalk.

- `size`:

  Number of rows of big data read per iteration in SGLD and SGHMC.

- `smax`:

  Maximum allowable tuning parameter sigma in AGG.

- `softabs_alpha`:

  Sharpness of the SoftAbs metric used by RMHMC, LMC, and MHMC.

- `subset_size`:

  Fraction of the dataset evaluated by the surrogate stage of Delayed
  Acceptance MCMC.

- `surrogate`:

  Surrogate-likelihood strategy in Delayed Acceptance MCMC (`"subset"`
  or a user function).

- `target_accept`:

  Target acceptance rate in autoMALA.

- `Temperature`:

  Heats up momentum in the first half of leapfrog steps in THMC.

- `T_max`:

  Maximum simulated time per Zig-Zag iteration.

- `Type`:

  Within Specs, used by Slice: `"Continuous"`, `"Nominal"`, or
  `"Ordinal"`.

- `vorticity`:

  Strength of the antisymmetric drift in Non-Reversible Overdamped
  Langevin.

- `w`:

  Mixture weight in AMM, step-size interval in Slice, AFSS, RSS,
  Refractive, GPSS, QSS, or stretch parameter in DEMC.

- `weightRange`:

  Range of guide-point weights in SQMC and SAMC.

- `Z`:

  \\T \times J\\ matrix or \\T \times J \times Nc\\ array of thinned
  samples for DEMC and the QMC family.

## References

Atchade, Y.F. (2006). "An Adaptive Version for the Metropolis Adjusted
Langevin Algorithm with a Truncated Drift". *Methodology and Computing
in Applied Probability*, 8, p. 235–254.
[doi:10.1007/s11009-006-8550-0](https://doi.org/10.1007/s11009-006-8550-0)

Bai, Y. (2009). "An Adaptive Directional Metropolis-within-Gibbs
Algorithm". Technical Report in Department of Statistics at the
University of Toronto.

Beskos, A., Roberts, G.O., Stuart, A.M., and Voss, J. (2008). "MCMC
Methods for Diffusion Bridges". *Stochastic Dynamics*, 8, p. 319–350.
[doi:10.1142/S0219493708002378](https://doi.org/10.1142/S0219493708002378)

Boyles, L.B. and Welling, M. (2012). "Refractive Sampling".

Craiu, R.V., Rosenthal, J., and Yang, C. (2009). "Learn From Thy
Neighbor: Parallel-Chain and Regional Adaptive MCMC". *Journal of the
American Statistical Association*, 104(488), p. 1454–1466.
[doi:10.1198/jasa.2009.tm08393](https://doi.org/10.1198/jasa.2009.tm08393)

Christen, J.A. and Fox, C. (2010). "A General Purpose Sampling Algorithm
for Continuous Distributions (the t-walk)". *Bayesian Analysis*, 5(2),
p. 263–282. [doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

Dutta, S. (2012). "Multiplicative Random Walk Metropolis-Hastings on the
Real Line". *Sankhya B*, 74(2), p. 315–342.
[doi:10.1007/s13571-012-0040-5](https://doi.org/10.1007/s13571-012-0040-5)

Duane, S., Kennedy, A.D., Pendleton, B.J., and Roweth, D. (1987).
"Hybrid Monte Carlo". *Physics Letters B*, 195, p. 216–222.
[doi:10.1016/0370-2693(87)91197-X](https://doi.org/10.1016/0370-2693%2887%2991197-X)

Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian Data
Analysis, Texts in Statistical Science, 2nd ed.". Chapman and Hall,
London.

Geman, S. and Geman, D. (1984). "Stochastic Relaxation, Gibbs
Distributions, and the Bayesian Restoration of Images". *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 6(6), p.
721–741.
[doi:10.1109/TPAMI.1984.4767596](https://doi.org/10.1109/TPAMI.1984.4767596)

Geyer, C.J. (1991). "Markov Chain Monte Carlo Maximum Likelihood". In
Keramidas, E.M. Computing Science and Statistics: Proceedings of the
23rd Symposium of the Interface. Fairfax Station VA: Interface
Foundation. p. 156–163.

Goodman, J. and Weare, J. (2010). "Ensemble Samplers with Affine
Invariance". *Communications in Applied Mathematics and Computational
Science*, 5(1), p. 65–80.
[doi:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)

Green, P.J. (1995). "Reversible Jump Markov Chain Monte Carlo
Computation and Bayesian Model Determination". *Biometrika*, 82, p.
711–732.
[doi:10.1093/biomet/82.4.711](https://doi.org/10.1093/biomet/82.4.711)

Haario, H., Laine, M., Mira, A., and Saksman, E. (2006). "DRAM:
Efficient Adaptive MCMC". *Statistics and Computing*, 16, p. 339–354.
[doi:10.1007/s11222-006-9438-0](https://doi.org/10.1007/s11222-006-9438-0)

Haario, H., Saksman, E., and Tamminen, J. (2001). "An Adaptive
Metropolis Algorithm". *Bernoulli*, 7, p. 223–242.
[doi:10.2307/3318737](https://doi.org/10.2307/3318737)

Hoffman, M.D. and Gelman, A. (2014). "The No-U-Turn Sampler: Adaptively
Setting Path Lengths in Hamiltonian Monte Carlo". *Journal of Machine
Learning Research*, 15, p. 1593–1623.

Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". *Journal of the
American Statistical Association*, 90(430), p. 773–795.
[doi:10.1080/01621459.1995.10476572](https://doi.org/10.1080/01621459.1995.10476572)

Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
Posterior Simulation with the Laplace-Metropolis Estimator". *Journal of
the American Statistical Association*, 92, p. 648–655.
[doi:10.1080/01621459.1997.10474016](https://doi.org/10.1080/01621459.1997.10474016)

Liu, J., Liang, F., and Wong, W. (2000). "The Multiple-Try Method and
Local Optimization in Metropolis Sampling". *Journal of the American
Statistical Association*, 95, p. 121–134.
[doi:10.1080/01621459.2000.10473908](https://doi.org/10.1080/01621459.2000.10473908)

Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller, E.
(1953). "Equation of State Calculations by Fast Computing Machines".
*Journal of Chemical Physics*, 21, p. 1087–1092.
[doi:10.1063/1.1699114](https://doi.org/10.1063/1.1699114)

Mira, A. (2001). "On Metropolis-Hastings Algorithms with Delayed
Rejection". *Metron*, Vol. LIX, n. 3-4, p. 231–241.

Murray, I., Adams, R.P., and MacKay, D.J. (2010). "Elliptical Slice
Sampling". *Journal of Machine Learning Research*, 9, p. 541–548.

Neal, R.M. (2003). "Slice Sampling" (with discussion). *Annals of
Statistics*, 31(3), p. 705–767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

Ritter, C. and Tanner, M. (1992). "Facilitating the Gibbs Sampler: the
Gibbs Stopper and the Griddy-Gibbs Sampler". *Journal of the American
Statistical Association*, 87, p. 861–868.
[doi:10.1080/01621459.1992.10476243](https://doi.org/10.1080/01621459.1992.10476243)

Roberts, G.O. and Rosenthal, J.S. (2009). "Examples of Adaptive MCMC".
*Computational Statistics and Data Analysis*, 18, p. 349–367.
[doi:10.1016/j.csda.2008.09.028](https://doi.org/10.1016/j.csda.2008.09.028)

Roberts, G.O. and Tweedie, R.L. (1996). "Exponential Convergence of
Langevin Distributions and Their Discrete Approximations". *Bernoulli*,
2(4), p. 341–363. [doi:10.2307/3318418](https://doi.org/10.2307/3318418)

Rosenthal, J.S. (2007). "AMCMC: An R interface for adaptive MCMC".
*Computational Statistics and Data Analysis*, 51, p. 5467–5470.
[doi:10.1016/j.csda.2006.11.012](https://doi.org/10.1016/j.csda.2006.11.012)

Smith, R.L. (1984). "Efficient Monte Carlo Procedures for Generating
Points Uniformly Distributed Over Bounded Region". *Operations
Research*, 32, p. 1296–1308.
[doi:10.1287/opre.32.6.1296](https://doi.org/10.1287/opre.32.6.1296)

Ter Braak, C.J.F. and Vrugt, J.A. (2008). "Differential Evolution Markov
Chain with Snooker Updater and Fewer Chains". *Statistics and
Computing*, 18(4), p. 435–446.
[doi:10.1007/s11222-008-9104-9](https://doi.org/10.1007/s11222-008-9104-9)

Tibbits, M., Groendyke, C., Haran, M., and Liechty, J. (2014).
"Automated Factor Slice Sampling". *Journal of Computational and
Graphical Statistics*, 23(2), p. 543–563.
[doi:10.1080/10618600.2013.791193](https://doi.org/10.1080/10618600.2013.791193)

Thompson, M.D. (2011). "Slice Sampling with Multivariate Steps".
<http://hdl.handle.net/1807/31955>

Vihola, M. (2012). "Robust Adaptive Metropolis Algorithm with Coerced
Acceptance Rate". *Statistics and Computing*, 22(5), p. 997–1008.
[doi:10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5)

Welling, M. and Teh, Y.W. (2011). "Bayesian Learning via Stochastic
Gradient Langevin Dynamics". *Proceedings of the 28th International
Conference on Machine Learning (ICML)*, p. 681–688.

## See also

[`AcceptanceRate`](https://robustecologies.github.io/lucifer/reference/AcceptanceRate.md),
[`as.initial.values`](https://robustecologies.github.io/lucifer/reference/as.initial.values.md),
[`as.parm.names`](https://robustecologies.github.io/lucifer/reference/as.parm.names.md),
[`BayesFactor`](https://robustecologies.github.io/lucifer/reference/BayesFactor.md),
[`Blocks`](https://robustecologies.github.io/lucifer/reference/Blocks.md),
[`BMK.Diagnostic`](https://robustecologies.github.io/lucifer/reference/BMK.Diagnostic.md),
[`Combine`](https://robustecologies.github.io/lucifer/reference/Combine.md),
[`Consort`](https://robustecologies.github.io/lucifer/reference/Consort.md),
[`dcrmrf`](https://robustecologies.github.io/lucifer/reference/dist.ContinuousRelaxation.md),
[`ESS`](https://robustecologies.github.io/lucifer/reference/ESS.md),
[`GIV`](https://robustecologies.github.io/lucifer/reference/GIV.md),
[`is.data`](https://robustecologies.github.io/lucifer/reference/is.data.md),
[`is.model`](https://robustecologies.github.io/lucifer/reference/is.model.md),
[`IterativeQuadrature`](https://robustecologies.github.io/lucifer/reference/IterativeQuadrature.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`lucifer.RAM`](https://robustecologies.github.io/lucifer/reference/lucifer.RAM.md),
[`LML`](https://robustecologies.github.io/lucifer/reference/LML.md),
[`MCSE`](https://robustecologies.github.io/lucifer/reference/MCSE.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(lucifer)
set.seed(666)

## ----- Synthetic Gaussian linear regression ---------------------------
## The simplest Bayesian model that every MCMC sampler in lucifer can
## fit: an unconstrained intercept, an unconstrained slope, and a
## log-transformed scale so that all three parameters live in R and
## gradient-based samplers (HMC, NUTS, MALA, RMHMC, ...) never collide
## with a positivity boundary.
N <- 100
ground_truth <- c(beta0 = 1, beta1 = 2, log.sigma = log(0.5))
X <- cbind(1, rnorm(N))
y <- as.vector(X %*% ground_truth[c("beta0", "beta1")]) +
     rnorm(N, 0, exp(ground_truth["log.sigma"]))

## Data list with PGF for random initial values via GIV()
PGF <- function(Data) c(rnorm(2, 0, 1), log(runif(1, 0.1, 1)))
Data <- list(N = N, J = 2, X = X, y = y,
             mon.names = "LP",
             parm.names = c("beta0", "beta1", "log.sigma"),
             PGF = PGF)

## Model specification (Gaussian likelihood, vague Gaussian priors)
Model <- function(parm, Data) {
     beta  <- parm[1:2]
     sigma <- exp(parm[3])
     mu    <- as.vector(Data$X %*% beta)
     LL    <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
     LP    <- LL +
              sum(dnorm(beta, 0, 10, log = TRUE)) +
              dnorm(parm[3], 0, 1, log = TRUE)
     list(LP = LP, Dev = -2 * LL, Monitor = LP,
          yhat = rnorm(length(mu), mu, sigma), parm = parm)
}

IV <- GIV(Model, Data, PGF = TRUE)

## ----- Canonical workflow (illustrated with NUTS) ---------------------
fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "NUTS", Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)

## ----- Sweep over every MCMC algorithm with default specifications ----
## The vector below enumerates every sampler in lucifer that runs with
## Specs = NULL on a continuous, unconstrained Gaussian target.
## Excluded:
##   * algorithms whose Specs are intrinsically mandatory (AGG, GG, IM,
##     RJ, RSS, SAMWG, SMWG, USAMWG, USMWG, SGLD, SGHMC) and MCMCMC
##     (which needs CPUs > 1 inside Specs);
##   * Polya-Gamma data augmentation (PG), which is restricted to
##     logistic / binomial regression and cannot fit a Gaussian
##     likelihood;
##   * NeuTra and flowMC, which require the torch package.
## All of these are documented in
## vignette("mcmc-algorithms", package = "lucifer").
algos <- c(
     ## Random-walk and adaptive Metropolis
     "RWM", "AM", "AMM", "RAM", "DRM", "DRAM", "DA", "MTM", "RDMH",
     "INCA", "PMCMC",
     ## Componentwise and Gibbs samplers
     "MWG", "AMWG", "ADMG", "Gibbs", "Zanella",
     ## Hamiltonian Monte Carlo
     "HMC", "AHMC", "NUTS", "HMCDA", "THMC", "MCHMC", "MEADS",
     "autoMALA", "MALT", "AAPS", "GIST",
     ## Langevin and stochastic gradient
     "MALA", "Barker", "NROLangevin",
     ## Slice samplers
     "Slice", "AFSS", "ESS", "OHSS", "UESS", "QSS", "LSS", "GPSS",
     ## Ensemble and population
     "AIES", "DEMC", "DREAM", "Zeus", "twalk",
     ## Quadratic Monte Carlo
     "QMC", "QMC3", "QMCN", "DQMC", "SQMC", "MAMC", "SAMC", "WMC",
     ## Hit-and-run and directional
     "HARM", "CHARM", "Refractive",
     ## Multimodal and tempering
     "NRPT", "SimTemp", "NRST", "WL",
     ## Piecewise-deterministic Markov processes
     "ZigZag", "BPS", "Boomerang", "RHMC",
     ## Riemannian and geometric
     "RMHMC", "LMC", "MHMC", "Relativistic",
     ## Constraint-handling and infinite-dimensional
     "ProjLang", "ProxMCMC", "pCN")

for (alg in algos) {
     fit <- lucifer(Model, Data, IV,
                    Iterations = 2000, Status = 1000, Thinning = 2,
                    Algorithm = alg, Specs = NULL, Chains = 3)
     plot(fit)
     ppc <- predict(fit, Model = Model, Data = Data)
     plot(ppc, Style = "Fitted", Data = Data)
     ppc_dens_overlay(ppc)
     plot(fit, ground_truth = ground_truth)
}
} # }
```
