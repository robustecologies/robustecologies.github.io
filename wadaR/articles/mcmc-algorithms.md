# A bestiary of the MCMC algorithms provided by lucifer

### Algorithms

The MCMC algorithms in the lucifer package are organized by algorithmic
family:

[Random-walk Metropolis and adaptive proposals](#rwmap)

- [Random-Walk Metropolis](#rwm)
- [Adaptive Metropolis](#am)
- [Adaptive-Mixture Metropolis](#amm)
- [Robust Adaptive Metropolis](#ram)
- [Delayed Rejection Adaptive Metropolis](#dram)
- [Delayed Rejection Metropolis](#drm)
- [Delayed Acceptance MCMC](#da) – *New*
- [Independence Metropolis](#im)
- [Multiple-Try Metropolis](#mtm)
- [Interchain Adaptation](#inca)

[Componentwise and Gibbs samplers](#cwgs)

- [Gibbs Sampler](#gibbs)
- [Metropolis-within-Gibbs](#mwg)
- [Adaptive Metropolis-within-Gibbs](#amwg)
- [Adaptive Directional Metropolis-within-Gibbs](#admg)
- [Griddy-Gibbs](#gg)
- [Adaptive Griddy-Gibbs](#agg)
- [Sequential Metropolis-within-Gibbs](#smwg)
- [Sequential Adaptive Metropolis-within-Gibbs](#samwg)
- [Updating Sequential Metropolis-within-Gibbs](#usmwg)
- [Updating Sequential Adaptive Metropolis-within-Gibbs](#usamwg)
- [Random Dive Metropolis-Hastings](#rdmh)
- [Zanella Locally Balanced](#zanella) – *New*
- [Polya-Gamma Data Augmentation](#pg) – *New*

[Hamiltonian Monte Carlo](#hmcalg)

- [Hamiltonian Monte Carlo](#hmc)
- [Adaptive Hamiltonian Monte Carlo](#ahmc)
- [No-U-Turn Sampler](#nuts)
- [Hamiltonian Monte Carlo with Dual-Averaging](#hmcda)
- [Tempered Hamiltonian Monte Carlo](#thmc)
- [Microcanonical Hamiltonian Monte Carlo](#mchmc) – *New*
- [MEADS](#meads) – *New*
- [autoMALA](#automala) – *New*
- [MALT](#malt) – *New*
- [AAPS](#aaps) – *New*
- [GIST](#gist) – *New*

[Langevin and stochastic gradient methods](#lsgm)

- [Metropolis-Adjusted Langevin Algorithm](#mala)
- [Stochastic Gradient Langevin Dynamics](#sgld)
- [Stochastic Gradient Hamiltonian Monte Carlo](#sghmc) – *New*
- [Barker Proposal](#barker) – *New*
- [Non-Reversible Overdamped Langevin](#nrolangevin) – *New*

[Slice samplers](#ss)

- [Slice Sampler](#slice)
- [Automated Factor Slice Sampler](#afss)
- [Elliptical Slice Sampler](#ess)
- [Oblique Hyperrectangle Slice Sampler](#ohss)
- [Reflective Slice Sampler](#rss)
- [Univariate Eigenvector Slice Sampler](#uess)
- [Quantile Slice Sampler](#qss) – *New*
- [Latent Slice Sampler](#lss) – *New*
- [Gibbsian Polar Slice Sampler](#gpss) – *New*

[Ensemble and population methods](#epm)

- [Affine-Invariant Ensemble Sampler](#aies)
- [Differential Evolution Markov Chain](#demc)
- [t-walk](#twalk)
- [Quadratic Monte Carlo](#qmc) – *New*
- [Third-Order Monte Carlo](#qmc3) – *New*
- [Nth-Order Monte Carlo](#qmcn) – *New*
- [Directed Quadratic Monte Carlo](#dqmc) – *New*
- [Quadratic Simplex Monte Carlo](#sqmc) – *New*
- [Modified Affine Monte Carlo](#mamc) – *New*
- [Affine Simplex Monte Carlo](#samc) – *New*
- [Walk Monte Carlo](#wmc) – *New*
- [Ensemble Slice Sampler](#zeus) – *New*
- [Differential Evolution Adaptive Metropolis](#dream) – *New*

[Hit-and-run and directional samplers](#hrds)

- [Hit-And-Run Metropolis](#harm)
- [Componentwise Hit-And-Run Metropolis](#charm)
- [Refractive Sampler](#refractive)

[Multimodal and transdimensional methods](#mmtm)

- [Metropolis-Coupled Markov Chain Monte Carlo](#mcmcmc)
- [Reversible-Jump](#rj)
- [Non-Reversible Parallel Tempering](#nrpt) – *New*
- [Simulated Tempering](#simtemp) – *New*
- [Non-Reversible Simulated Tempering](#nrst) – *New*
- [Wang-Landau](#wanglandau) – *New*

[Piecewise-deterministic samplers](#pdmps)

- [Zig-Zag Sampler](#zigzag) – *New*
- [Bouncy Particle Sampler](#bps) – *New*
- [Boomerang Sampler](#boomerang) – *New*
- [Randomized HMC](#rhmc) – *New*

[Riemannian and geometric MCMC](#riemannian) – *New*

- [Riemannian Manifold HMC](#rmhmc) – *New*
- [Lagrangian Monte Carlo](#lmc) – *New*
- [Magnetic HMC](#mhmc) – *New*
- [Relativistic Monte Carlo](#relativistic) – *New*

[Constraint-handling samplers](#constraints) – *New*

- [Projected Langevin](#projlang) – *New*
- [ProxMCMC](#proxmcmc) – *New*

[Intractable likelihoods and infinite-dimensional targets](#ilidt) –
*New*

- [Pseudo-Marginal MCMC](#pmcmc) – *New*
- [Preconditioned Crank-Nicolson](#pcn)

## Markov Chain Monte Carlo

Markov chain Monte Carlo (MCMC) algorithms, also called samplers, are
numerical approximation algorithms. There are a large number of MCMC
algorithms, too many to review here. Popular families (which are often
non-distinct) include Gibbs sampling, Metropolis-Hastings, slice
sampling, Hamiltonian Monte Carlo, and many others. Though the name is
misleading, Metropolis-within-Gibbs (MWG) was developed first
(Metropolis, Rosenbluth, and Teller 1953), and Metropolis-Hastings was a
generalization of MWG (Hastings 1970). All MCMC algorithms are known as
special cases of the Metropolis-Hastings algorithm. Regardless of the
algorithm, the goal in Bayesian inference is to maximize the
unnormalized joint posterior distribution and collect samples of the
target distributions, which are marginal posterior distributions, later
to be used for inference.

The most generalizable MCMC algorithm is the Metropolis-Hastings (MH)
generalization of the MWG algorithm. The MH algorithm extended MWG to
include asymmetric proposal distributions. For years, the main
disadvantage of the MWG and MH algorithms was that the proposal variance
(see below) had to be tuned manually, and therefore other MCMC
algorithms have become popular because they do not need to be tuned.

Gibbs sampling became popular for Bayesian inference, though it requires
conditional sampling of conjugate distributions, so it is precluded from
non-conjugate sampling in its purest form. Gibbs sampling also suffers
under high correlations (Gilks and Roberts 1996). Due to these
limitations, Gibbs sampling is less generalizable than random-walk
Metropolis (RWM), though RWM and other algorithms are not immune to
problems with correlation. The original slice sampling algorithm of Neal
(1997) is a special case of Gibbs sampling that samples a distribution
by sampling uniformly from the region under the plot of its density
function, and is more appropriate with bounded distributions that cannot
approach infinity, though the improved Slice sampler of Neal (2003) is
available here.

### Blockwise Sampling

Usually, there is more than one target distribution, in which case it
must be determined whether it is best to sample from target
distributions individually, in groups, or all at once. Block updating
refers to splitting a multivariate vector into groups called blocks, and
each block is sampled separately. A block may contain one or more
parameters.

Parameters are usually grouped into blocks such that parameters within a
block are as correlated as possible, and parameters between blocks are
as independent as possible. This strategy retains as much of the
parameter correlation as possible for blockwise sampling, as opposed to
componentwise sampling where parameter correlation is often ignored. The
PosteriorChecks function can be used on the output of previous runs to
find highly correlated parameters, and the Blocks function may be used
to create blocks based on posterior correlation.

Advantages of blockwise sampling are that a different MCMC algorithm may
be used for each block (or parameter, for that matter), creating a more
specialized approach (though different algorithms by block are not
supported here), the acceptance of a newly proposed state is likely to
be higher than sampling from all target distributions at once in high
dimensions, and large proposal covariance matrices can be reduced in
size, which is most helpful again in high dimensions.

Disadvantages of blockwise sampling are that correlations probably exist
between parameters between blocks, and each block is updated while
holding the other blocks constant, ignoring these correlations of
parameters between blocks. Without simultaneously taking everything into
account, the algorithm may converge slowly or never arrive at the proper
solution. However, there are instances when it may be best when
everything is not taken into account at once, such as in state-space
models. Also, as the number of blocks increases, more computation is
required, which slows the algorithm. In general, blockwise sampling
allows a more specialized approach at the expense of accuracy,
generalization, and speed.

Blockwise sampling is offered in the following algorithms: Adaptive
Metropolis-within-Gibbs (AMWG), Adaptive-Mixture Metropolis (AMM),
Automated Factor Slice Sampler (AFSS), Elliptical Slice Sampler (ESS),
Hit-And-Run Metropolis (HARM), Metropolis-within-Gibbs (MWG),
Random-Walk Metropolis (RWM), Robust Adaptive Metropolis (RAM), and the
Univariate Eigenvector Slice Sampler (UESS).

### Markov chain Properties

Only the basics of Markov chain properties are introduced here. A Markov
chain is Markovian when the current iteration depends only on the
previous iteration. Many (but not all) adaptive algorithms are merely
chains but not Markov chains when the adaptation is based on the history
of the chains, not just the previous iteration. A Markov chain is said
to be aperiodic when it is not repeating a cycle. A Markov chain is
considered irreducible when it is possible to go from any state to any
other state, though not necessarily in one iteration. A Markov chain is
said to be recurrent if it will eventually return to a given state with
probability 1, and it is positive recurrent if the expected return time
is finite, and null recurrent otherwise. The ergodic theorem states that
a Markov chain is ergodic when it is aperiodic, irreducible, and
positive recurrent.

The non-Markovian chains of an adaptive algorithm that adapt based on
the history of the chains should have two conditions: containment and
diminishing adaptation. Containment is difficult to implement. The
condition of diminishing adaptation is fulfilled when the amount of
adaptation diminishes with the length of the chain. Diminishing
adaptation can be achieved when the proposal variances become smaller or
by decreasing the probability of performing adaptations with more
iterations (Roberts and Rosenthal 2007).

### Algorithm Summary Table

Each algorithm, presented above, has an algorithm summary table that is
meant to provide a quick overview and facilitate comparisons between
algorithms. Below is a description of the contents of these tables.

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is a consideration of the number of accepted proposals and the number of rejected proposals. It is usually a ratio of the number of proposals to the number of iterations. When known, the optimal acceptance rate and recommended interval is reported. |
| Applications | Since no algorithm is suitable to every problem, this is meant to provide guidance to the applicability of the algorithm. |
| Difficulty | The difficulty of using the algorithm is considered for a beginner. Considerations include the number of algorithm specifications, tuning the algorithm, whether or not it may be used as a final algorithm, and general experiences with use. |
| Final Algorithm? | This description is usually ‘Yes’, ‘No’, or ‘User Discretion’, followed by an explanation. When a ‘final algorithm’ seems to have converged, the samples may be used for inference. When it is not a final algorithm, then one is suggested. |
| Proposal | This describes how the proposal is generated, and is usually ‘Multivariate’ or ‘Componentwise’, though other descriptions exist as well. In general, a multivariate proposal means that each iteration, a proposal is generated for all parameters that takes correlation into account, usually from a multivariate normal distribution and a proposal covariance matrix. Componentwise proposals usually indicate that a proposal is made for each parameter, without considering correlation. A componentwise algorithm must evaluate the model a number of times equal to the number of parameters, per iteration. Componentwise algorithms are consierably slower per iteration, and speed per iteration decreases as the number of parameters increases. |

### Stopping Rules

After a MCMC algorithm is selected, a model is specified and the priors
are updated given the likelihood, and numerous iterations. A question
with MCMC is: when to stop iterating/updating? Many stopping rules have
been proposed. The most popular single-chain stopping rule involves
Monte Carlo Standard Error (MCSE). An alternative when this is too
difficult may be when effective sample size (ESS) reaches a target. The
most popular parallel-chain stopping rule uses the Gelman-Rubin
Diagnostic. MCMC Diagnostics

After the user decides to stop iterating MCMC and updating the model,
MCMC diagnostics are performed to assess non-convergence. Markov (and
non-Markovian) Chains cannot be confirmed to have converged, so
diagnostics are used to look for signs of non-convergence. Numerous MCMC
diagnostics have been proposed. It is customary to observe trace plots
and autocorrelation plots of the chains, assess stationarity with the
BMK Diagnostic, and search for signs of non-convergence with other
diagnostics. The lucifer software package offers a variety of MCMC
diagnostics. References for the introductory material above and for
every algorithm-specific section that follows are collected at the end
of this vignette under [References](#references).

## Random-walk Metropolis and adaptive proposals learn from chain history

These methods generate proposals by perturbing the current state with
symmetric or nearly symmetric jumps, typically drawn from a multivariate
normal centered on the current position. Adaptive variants estimate the
proposal covariance from chain history to improve mixing without
requiring gradient information. The family ranges from the foundational
random-walk Metropolis to sophisticated delayed-rejection and
multiple-try extensions that reduce the cost of rejection.

### Random-Walk Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is moderately difficult for a beginner because the proposal covariance matrix must be tuned. |
| Final Algorithm? | Yes. |
| Proposal | Blockwise or Multivariate. |

Random-Walk Metropolis (RWM) is a multivariate extension of
Metropolis-within-Gibbs (MWG). Multiple parameters usually exist, and
therefore correlations may occur between the parameters. Many MCMC
algorithms attempt to estimate multivariate proposals, usually taking
correlations into account through a covariance matrix. Each iteration, a
multivariate proposal is made by taking a draw from a multivariate
normal distribution, given a proposal covariance matrix.

RWM does not have any required algorithm specifications, though
blockwise sampling may be specified with `B`, which accepts a list of
proposal covariance matrices equal in length to the number of blocks. By
default, blockwise sampling is not performed, so all parameters are
updated with one multivariate proposal.

Given the number of dimensions (`K`) or parameters, the optimal scale of
the proposal covariance, also called the jumping kernel, has been
reported as \\2.4^2/K\\ based on the asymptotic limit of
infinite-dimensional Gaussian target distributions that are independent
and identically-distributed [\[1\]](#ref1). In applied settings, each
problem is different, so the amount of correlation varies between
variables, target distributions may be non-Gaussian, the target
distributions may be non-IID, and the scale should be optimized. There
are algorithms in statistical literature that attempt to optimize this
scale, such as the Robust Adaptive Metropolis (RAM) algorithm.

There have been numerous methods introduced for tuning the proposal
covariance matrix. Many adaptive MCMC algorithms such as Adaptive
Metropolis (AM), Adaptive-Mixture Metropolis (AMM), and RAM will tune
the proposal covariance matrix. Alternatively, initially specify an
identity matrix with small-scale diagonal values such as `1e-3` as the
proposal covariance matrix, update with RWM, and then update again but
this time with the observed covariance of the historical samples. Done
repeatedly, this may arrive at an acceptable proposal covariance matrix,
suitable for longer and more serious updates.

Since RWM is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal covariance
    \\\Sigma\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Propose \\\theta^\* \sim \mathcal{N}(\theta\_{t-1}, \Sigma)\\
    2.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    3.  Draw \\u \sim \text{Uniform}(0,1)\\; if \\u \< \alpha\\, accept
        \\\theta_t = \theta^\*\\; otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

The optimal scaling of \\\Sigma\\ is \\2.4^2/d \cdot
\hat{\Sigma}\_{\text{target}}\\, where \\d\\ is the parameter dimension
and \\\hat{\Sigma}\_{\text{target}}\\ approximates the posterior
covariance.

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| B | list() | No blocking; the default uses a single multivariate Gaussian proposal. Supply a list of disjoint parameter index vectors to enable blockwise updates. |

Default `Specs` for RWM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RWM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AM](#am)
- [AMM](#amm)
- [MWG](#mwg)
- [RAM](#ram)

### Adaptive Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications, and these are easy to specify. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The RWM algorithm is recommended as the final algorithm, though AM may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Multivariate. |

The Adaptive Metropolis (AM) algorithm of [\[2\]](#ref2) is an extension
of Random-Walk Metropolis (RWM) that adapts based on the observed
covariance matrix from the history of the chains. AM is historically
significant as the first adaptive MCMC algorithm.

AM has two algorithm specifications: `Adaptive` is the iteration at
which adaptation begins, and `Periodicity` is the frequency in
iterations of adaptation. The user should not allow AM to adapt
immediately, since AM adapts based on the observed covariance matrix of
historical and accepted samples. Since enough samples are needed to
obtain a valid covariance matrix before adaptation, a small covariance
matrix is often used initially to encourage a high acceptance rate.

As recommended by [\[2\]](#ref2), there are two tricks that may be used
to assist the AM algorithm in the beginning. Although the suggested
“greedy start” method is not used here, the second suggested trick is
used, which is to shrink the proposal as long as the acceptance rate is
less than 5%, and there have been at least five acceptances.
[\[2\]](#ref2) suggest loosely that if “it has not moved enough during
some number of iterations, the proposal could be shrunk by a constant
factor”. For each iteration that the acceptance rate is less than 5% and
that the AM algorithm is used but the current iteration is prior to
adaptation, Laplace’s Demon multiplies the proposal covariance or
variance by (1 - 1/Iterations). Over pre-adaptive time, this encourages
a smaller proposal covariance or variance to increase the acceptance
rate so that when adaptation begins, the observed covariance or variance
of the chains will not be constant, and then shrinkage will cease and
adaptation will take it from there.

AM is best suited for a model with small or medium dimensions (number of
parameters). The Adaptive-Mixture Metropolis (AMM) of [\[3\]](#ref3) and
Robust Adaptive Metropolis (RAM) of Vihola (2011) are extensions of the
AM algorithm. If AM is used for adaptation, then the final, non-adaptive
algorithm should be RWM.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, initial proposal covariance
    \\\Sigma_0\\, scaling \\s_d = 2.4^2 / d\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  If \\t \<\\ `Adaptive`: use \\\Sigma_0\\; if acceptance rate
        \\\< 5\\\\, shrink \\\Sigma_0 \leftarrow (1 - 1/t) \cdot
        \Sigma_0\\
    2.  If \\t \geq\\ `Adaptive` and \\t \bmod\\ `Periodicity` \\= 0\\:
        set \\\Sigma_t \leftarrow s_d \cdot \text{Cov}(\theta_0, \ldots,
        \theta_t) + s_d \cdot \delta \cdot I_d\\
    3.  Propose \\\theta^\* \sim \mathcal{N}(\theta\_{t-1}, \Sigma_t)\\
    4.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Adaptive | floor(Iterations/2) (= 1000) | Lets the empirical covariance accumulate over the first half of the run before adaptation begins, avoiding pathological proposals from a noisy early sample covariance. |
| Periodicity | 1 | Re-estimate the proposal covariance every iteration; cheap when the dimension is small and accelerates burn-in. |

Default `Specs` for AM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMM](#amm)
- [INCA](#inca)
- [RAM](#ram)
- [RWM](#rwm)

### Adaptive-Mixture Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters, or a larger number of blocked parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. The algorithm specifications are easy to specify. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The RWM algorithm is recommended as the final algorithm, though AMM may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Blockwise or Multivariate. |

The Adaptive-Mixture Metropolis (AMM) algorithm is an extension by
[\[3\]](#ref3) of the Adaptive Metropolis (AM) algorithm of
[\[2\]](#ref2). AMM differs from the AM algorithm in two respects.
First, AMM updates a scatter matrix based on the cumulative current
parameters and the cumulative associated outer-products, and these are
used to generate a multivariate normal proposal. This is more efficient
with large numbers of parameters adapting over many iterations,
especially with frequent adaptations, and results in a much faster
algorithm. The second (and main) difference, is that the proposal is a
mixture. The two mixture components are adaptive multivariate and
static/symmetric univariate proposals. The mixture is determined at each
iteration with a mixture weight. The mixture weight must be in the
interval (0,1\], and it defaults to 0.05, as in [\[3\]](#ref3). A higher
value of the mixture weight is associated with more static/symmetric
univariate proposals, and a lower weight is associated with more
adaptive multivariate proposals. The algorithm will be unable to include
the multivariate mixture component until it has accumulated some
history, and models with more parameters will take longer to be able to
use adaptive multivariate proposals.

AMM has five algorithm specifications: `Adaptive` is the iteration in
which adaptation begins, `B` optionally accepts a list of blocked
parameters and defaults to NULL, `n` is the number of previous
iterations, `Periodicity` is the frequency in iterations of adaptation,
and `w` is the weight of the small-variance mixture component. `B`
allows the user to organize parameters into blocks. `B` accepts a list,
in which each component is a block and accepts a vector that consists of
numbers that point to the associated parameters in `parm.names`. `B`
defaults to NULL, in which case blocking does not occur. When blocking
does occur, the proposal covariance matrix may be either NULL or a list
in which each component is the covariance matrix for a block. As more
blocks are added, the algorithm becomes closer to Adaptive
Metropolis-within-Gibbs (AMWG).

The advantages of AMM over AMWG (when `B`=NULL) are that it takes
correlation into account as it adapts, and is much faster to update each
iteration. The disadvantages are that AMWG does not require a burn-in
period before it can begin to adapt, and more information must be
learned in the covariance matrix to adapt properly. Disadvantages of AMM
compared to Robust Adaptive Metropolis (RAM) are that RAM does not
require a burn-in period before it can begin to adapt, RAM is more
likely to better handle multimodal or heavy-tailed targets, and RAM also
adapts to the shape of the target distributions and coerces the
acceptance rate. If AMM is used for adaptation, then the final,
non-adaptive algorithm should be Random-Walk Metropolis (RWM).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, scatter matrix \\S = 0\\,
    mixture weight \\w\\, iteration count \\n\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  If \\t \geq\\ `Adaptive` and \\t \bmod\\ `Periodicity` \\= 0\\:
        update proposal covariance \\\Sigma_t\\ from cumulative scatter
        matrix
    2.  Draw \\u \sim \text{Uniform}(0,1)\\
    3.  If \\u \> w\\ and \\\Sigma_t\\ is available: propose \\\theta^\*
        \sim \mathcal{N}(\theta\_{t-1}, \Sigma_t)\\ (adaptive
        multivariate)
    4.  Otherwise: propose \\\theta^\* \sim \mathcal{N}(\theta\_{t-1},
        0.01/d \cdot I_d)\\ (static symmetric)
    5.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    6.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Adaptive | floor(Iterations/2) (= 1000) | Same rationale as AM: defer the adaptive multivariate component until enough history has accumulated. |
| B | NULL | Single-block update by default; supply a list of index vectors to enable blockwise sampling. |
| n | 0 | Initial counter for the cumulative scatter matrix. |
| Periodicity | 1 | Update the cumulative scatter matrix every iteration. |
| w | 0.05 | Mixture weight on the static symmetric component, exactly the value advocated by Roberts and Rosenthal (2009). |

Default `Specs` for AMM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AMM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AM](#am)
- [AMWG](#amwg)
- [RAM](#ram)
- [RWM](#rwm)

### Robust Adaptive Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The user specifies the target acceptance rate, α\* (23.4% is recommended). The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters, or a larger number of blocked parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications, and these are easy to specify. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The RWM algorithm is recommended as the final algorithm, though RAM may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Blockwise or Multivariate. |

The Adaptive Metropolis (AM) and Adaptive-Mixture Metropolis (AMM)
algorithms adapt the scale of the proposal distribution to attain a
theoretical acceptance rate. However, these algorithms are unable to
adapt to the shape of the target distribution. The Robust Adaptive
Metropolis (RAM) algorithm estimates the shape of the target
distribution and simultaneously coerces the acceptance rate
[\[4\]](#ref4). If the acceptance probability, α, is less (or greater)
than an acceptance rate target, α∗, then the proposal distribution is
shrunk (or expanded). Matrix S is computed as a rank one Cholesky
update. Therefore, the algorithm is computationally efficient up to a
relatively high dimension. The AM and AMM algorithms require a burn-in
period prior to adaptation, so that these algorithms can adapt to the
sample covariance. RAM does not require a burn-in period prior to
adaptation. RAM allows the user the option of using the traditional
normally-distributed proposals, or t-distributed proposals for
heavier-tailed target densities. Unlike AM and AMM, RAM can cope with
targets having arbitrarily heavy tails, and handles multimodal targets
better than AM. The user is still assumed to know and specify the target
acceptance rate.

RAM has five algorithm specifications: `alpha.star` is the target
acceptance rate, `B` optionally accepts a list of blocked parameters and
defaults to NULL, `Dist` is the target distribution as either “N” for
normal or “t” for the Student t with 5 degrees of freedom, `gamma`
accepts a scalar in the interval (0.5,1\] and controls the decay of
adaptation (0.66 is recommended), and `n` is the number of previous
iterations. RAM adapts only when the variance-covariance matrix is
positive-definite.

The advantages of RAM over AMM are that RAM does not require a burn-in
period before it can begin to adapt. RAM is more likely to better handle
multimodal or heavy-tailed targets, adapts to the shape of the target
distributions, and attempts to coerce the acceptance rate. The
advantages of RAM over Adaptive Metropolis-within-Gibbs (AMWG) are that
RAM takes correlations into account, and is much faster to update each
iteration. The disadvantage of RAM compared to AMWG is that more
information must be learned in the covariance matrix to adapt properly,
and frequent adaptation may be desirable, but slow. If RAM is used for
adaptation, then the final, non-adaptive algorithm should be Random-Walk
Metropolis (RWM).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, lower-triangular Cholesky
    factor \\S\\ (from initial covariance), target acceptance
    \\\alpha^\*\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw \\u \sim \mathcal{N}(0, I)\\ (or \\t\\-distributed if
        `Dist = "t"`)
    2.  Propose \\\theta^\* = \theta\_{t-1} + S \\ u\\
    3.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    4.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
    5.  **Adapt** \\S\\ via rank-one Cholesky update: \\S\_{t+1}
        S\_{t+1}^\top = S_t \bigl(I + \eta_t (\alpha - \alpha^\*)
        \frac{u u^\top}{\\u\\^2}\bigr) S_t^\top\\, where \\\eta_t =
        t^{-\gamma}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| alpha.star | 0.234 | Roberts and Rosenthal optimal acceptance rate for high-dimensional Gaussian targets; raise to 0.44 in one dimension. |
| B | NULL | Multivariate proposal by default; supply blocks for very high dimensions. |
| Dist | “N” | Gaussian proposal; switch to “t” for heavy-tailed targets. |
| gamma | 0.66 | Decay exponent in (0.5, 1\] suggested by Vihola (2012); guarantees diminishing adaptation. |
| n | 0 | Adaptation starts immediately because RAM does not need a burn-in to estimate covariance. |

Default `Specs` for RAM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RAM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AM](#am)
- [AMM](#amm)
- [AMWG](#amwg)
- [RWM](#rwm)

### Delayed Rejection Adaptive Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications, and these are easy to specify. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The RWM algorithm is recommended as the final algorithm, though DRAM may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Multivariate. Whenever a proposal is rejected, an alternate proposal is attempted. |

The Delayed Rejection Adaptive Metropolis (DRAM) algorithm is merely the
combination of both Delayed Rejection Metropolis (DRM) and Adaptive
Metropolis (AM) [\[5\]](#ref5). DRAM has been demonstrated to be robust
in extreme situations where DRM or AM fail separately. [\[5\]](#ref5)
present an example involving ordinary differential equations in which
least squares could not find a stable solution, and DRAM did well.

DRAM has two algorithm specifications: `Adaptive` is the iteration in
which DRAM becomes adaptive, and `Periodicity` is the frequency in
iterations of adaptation.

The DRAM algorithm is useful to assist AM when the acceptance rate is
low. As an alternative, the Adaptive-Mixture Metropolis (AMM) is an
extension of the AM algorithm that includes a mixture of proposals, and
one mixture component has a small proposal standard deviation to assist
in overcoming initially low acceptance rates. If DRAM is used for
adaptation, then the final algorithm should be Random-Walk Metropolis
(RWM).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal covariance
    \\\Sigma\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Stage 1**: propose \\\theta^{(1)} \sim
        \mathcal{N}(\theta\_{t-1}, \Sigma)\\
    2.  Compute \\\alpha_1 = \min\\\bigl(1,\\ p(\theta^{(1)} \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    3.  With probability \\\alpha_1\\, accept \\\theta_t =
        \theta^{(1)}\\ and go to step 2
    4.  **Stage 2 (delayed rejection)**: propose \\\theta^{(2)} \sim
        \mathcal{N}(\theta\_{t-1}, \Sigma/4)\\ (halved scale)
    5.  Compute \\\alpha_2 = \min\\\bigl(1,\\ \frac{p(\theta^{(2)}
        \mid y) \cdot q(\theta^{(1)} \mid \theta^{(2)}) \cdot (1 -
        \alpha(\theta^{(2)}, \theta^{(1)}))}{p(\theta\_{t-1} \mid y)
        \cdot q(\theta^{(1)} \mid \theta\_{t-1}) \cdot (1 -
        \alpha_1)}\bigr)\\
    6.  With probability \\\alpha_2\\, accept \\\theta_t =
        \theta^{(2)}\\; otherwise \\\theta_t = \theta\_{t-1}\\
    7.  If \\t \geq\\ `Adaptive` and \\t \bmod\\ `Periodicity` \\= 0\\:
        update \\\Sigma\\ from sample history (as in AM)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Adaptive | floor(Iterations/2) (= 1000) | Same as AM: defer covariance adaptation until enough history is available. |
| Periodicity | 1 | Refresh the proposal covariance every iteration. |

Default `Specs` for DRAM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DRAM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AM](#am)
- [AMM](#amm)
- [DRM](#drm)
- [RWM](#rwm)

### Delayed Rejection Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. The proposal covariance matrix should have been solved with an adaptive algorithm such as DRAM. |
| Difficulty | This algorithm is easy for a beginner when the proposal covariance has been tuned with another algorithm. Otherwise, it may be tedious for the user to tune the proposal covariance matrix. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. Whenever a proposal is rejected, an alternate proposal is attempted. |

The Delayed Rejection Metropolis (DRM or DR) algorithm is a Random-Walk
Metropolis (RWM) with one, small twist. Whenever a proposal is rejected,
the DRM algorithm will try one or more alternate proposals, and correct
for the probability of this conditional acceptance. By delaying
rejection, autocorrelation in the chains may be decreased, and the
algorithm is encouraged to move. The additional calculations may slow
each iteration of the algorithm in which the first set of proposals is
rejected, but it may also converge faster. For more information on DRM,
see [\[6\]](#ref6).

DRM does not have any algorithm specifications.

DRM may be considered to be an adaptive MCMC algorithm, because it
adapts the proposal based on a rejection. However, DRM does not violate
the Markov property, because the proposal is based on the current state.
DRM is not considered to be an adaptive MCMC algorithm here, because it
is not adapting to the target distribution by considering previous
states in the Markov chain, but merely makes more attempts from the
current state.

lucifer also temporarily shrinks the proposal covariance arbitrarily by
50% for delayed rejection. A smaller proposal covariance is more likely
to be accepted, and the goal of delayed rejection is to increase
acceptance. In the long-term, a proposal covariance that is too small is
undesirable, and so it is only used in this case to assist acceptance.

Since DRM is non-adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal covariance
    \\\Sigma\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Stage 1**: propose \\\theta^{(1)} \sim
        \mathcal{N}(\theta\_{t-1}, \Sigma)\\
    2.  Compute \\\alpha_1 = \min\\\bigl(1,\\ p(\theta^{(1)} \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    3.  With probability \\\alpha_1\\, accept \\\theta_t =
        \theta^{(1)}\\ and go to step 2
    4.  **Delayed rejection**: propose \\\theta^{(2)} \sim
        \mathcal{N}(\theta\_{t-1}, \Sigma/4)\\ (reduced scale)
    5.  Compute corrected acceptance \\\alpha_2\\ accounting for the
        conditional probability of reaching stage 2
    6.  With probability \\\alpha_2\\, accept \\\theta_t =
        \theta^{(2)}\\; otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| (none) | (none) | DRM has no algorithm-specific specifications. The proposal scale is taken from `Covar` and the secondary delayed proposal is generated automatically by halving the primary scale. |

Default `Specs` for DRM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DRM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [RWM](#rwm)

### Delayed Acceptance MCMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The overall acceptance rate is typically lower than RWM because proposals must pass two stages. The observed acceptance rate may be suitable in the interval \[10%,40%\]. |
| Applications | This algorithm is designed for models where the full posterior is expensive to evaluate. A cheap surrogate (such as a subsampled likelihood or a Gaussian process emulator) filters out obviously bad proposals before the full model is evaluated. |
| Difficulty | This algorithm is moderately difficult for a beginner. The user must either provide a surrogate Model function or set surrogate=“subset” for automatic subset-data approximation. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Delayed Acceptance MCMC (DA) was introduced by [\[7\]](#ref7) as a
two-stage Metropolis-Hastings algorithm that reduces the number of
expensive full-model evaluations. The first stage evaluates a cheap
surrogate of the log-posterior and applies the standard MH accept/reject
criterion. Only proposals that survive stage 1 proceed to stage 2, where
the full (expensive) model is evaluated and a correction factor is
applied to ensure that the chain targets the exact posterior despite
using the surrogate in stage 1. The correction ratio is
\\\exp\\\bigl\[(LP(\theta^\*) - LP_s(\theta^\*)) - (LP(\theta) -
LP_s(\theta))\bigr\]\\, which accounts for the discrepancy between the
surrogate and the true posterior at both the proposed and current
positions.

DA has three algorithm specifications: `surrogate` is either a function
with the same signature as Model that returns a cheaper approximation,
or the string `"subset"` which causes the algorithm to construct an
automatic subset-data surrogate by evaluating the Model on a random
fraction of the data (controlled by `subset_size`, default 0.1) and
scaling the log-likelihood accordingly. The parameter `base_algorithm`
(default `"RWM"`) selects the proposal mechanism. The proposal
covariance is handled identically to RWM.

The advantage of DA is that it can dramatically reduce computational
cost when the surrogate is much cheaper than the full model and
reasonably accurate: most bad proposals are rejected at stage 1 without
ever touching the expensive likelihood. The disadvantage is that the
overall acceptance rate is the product of the two stage acceptance
rates, so if the surrogate is a poor approximation, the chain mixes
slowly. Since DA is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal covariance
    \\\Sigma\\, surrogate function \\\tilde{p}\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Propose \\\theta^\* \sim \mathcal{N}(\theta\_{t-1}, \Sigma)\\
    2.  **Stage 1** (cheap): compute \\\alpha_1 = \min\\\bigl(1,\\
        \tilde{p}(\theta^\* \mid y) / \tilde{p}(\theta\_{t-1} \mid
        y)\bigr)\\
    3.  With probability \\\alpha_1\\, proceed to stage 2; otherwise
        retain \\\theta\_{t-1}\\
    4.  **Stage 2** (expensive): compute \\\alpha_2 = \min\\\bigl(1,\\
        \exp\bigl\[(LP^\* - \tilde{LP}^\*) - (LP\_{t-1} -
        \tilde{LP}\_{t-1})\bigr\]\bigr)\\
    5.  With probability \\\alpha_2\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| surrogate | “subset” | Cheap surrogate built from a random subset of the data; the most generic option requiring no user-supplied function. |
| subset_size | 0.1 | 10% of the rows of `Data$y` are evaluated by the surrogate stage. Small enough to provide a meaningful speed-up, large enough to keep correlation with the full likelihood. |
| base_algorithm | “RWM” | Random-walk Metropolis is the simplest base sampler that delayed acceptance can wrap; switch to gradient-based samplers when gradients are available. |

Default `Specs` for DA. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DA", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [RWM](#rwm)
- [DRAM](#dram)

### Independence Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small to medium number of parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications, and these are easy to specify. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Proposed by [\[9\]](#ref9) and popularized by [\[10\]](#ref10), the
Independence Metropolis (IM) algorithm (also called the independence
sampler) is an algorithm in which the proposal distribution does not
depend on the previous state or iteration. The proposal distribution
must be a good approximation of the target distribution for IM to
perform well, and the proposal distribution should have slightly heavier
tails than the target distribution. IM is used most often to obtain
additional posterior samples given an algorithm that has already
converged.

IM has one algorithm specification: `mu`. The usual approach to IM is to
update the model with Laplace Approximation, and then supply the
posterior mode and covariance to IM. The posterior mode vector of
Laplace Approximation becomes the `mu` argument in the algorithm
specifications for IM. The covariance matrix from Laplace Approximation
is expanded by multiplying it by 1.1 so that it has heavier tails. Each
iteration, IM draws from a multivariate normal distribution as the
proposal distribution. Alternatively, posterior means and covariances
may be used from other algorithms, such as other MCMC algorithms.

Since IM is non-adaptive and uses a proposal distribution that remains
fixed for all iterations, it may be used as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal mean \\\mu\\ (e.g.,
    from Laplace Approximation), proposal covariance \\\Sigma\\
    (slightly inflated, e.g., \\1.1 \times\\ posterior covariance)
2.  **For** \\t = 1, \ldots, T\\:
    1.  Propose \\\theta^\* \sim \mathcal{N}(\mu, \Sigma)\\ (independent
        of \\\theta\_{t-1}\\)
    2.  Compute importance weight ratio \\w = \frac{p(\theta^\* \mid y)
        / q(\theta^\*)}{p(\theta\_{t-1} \mid y) / q(\theta\_{t-1})}\\
    3.  Compute \\\alpha = \min(1, w)\\
    4.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| mu | IV | Independence Metropolis has no defaults; the proposal mean must be supplied. We use the initial values, which is a sensible mode-centered choice when no Laplace fit is available. |

Default `Specs` for IM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "IM", Specs = list(mu = IV),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- none

### Multiple-Try Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is above 44%, and increases with the number of tries. The observed acceptance rate may be suitable in the interval \[10%,90%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is moderate in difficulty because the number of proposals and user-specified proposal variance may need tuning. Additional optional algorithm specifications control parallelization. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Multiple-Try Metropolis (MTM) algorithm was introduced in
[\[11\]](#ref11) as a componentwise Metropolis algorithm with an
improved acceptance rate due to an increased proposal variance and
number of proposals. The user specifies that K proposals will be made.
For each parameter at each iteration, K normally-distributed proposals
are made around the current parameter, scaled according to the
user-specified proposal variance. Each proposal is weighted according to
its unnormalized joint posterior density. One proposal is selected
randomly, with probability proportional to the weights. A reference set
of length K is created, in which the first K-1 elements are draws from a
normal distribution centered around the selected proposal, again
according to proposal variance. The last element, K, is the selected
proposal itself. A Metropolis step is performed for the weighted
reference set. If the weighted reference set is accepted, then the
selected proposal becomes the new value for the parameter.

MTM has four algorithm specifications: `K` is the number of proposals,
`CPUs` is the number of central processing units, `Packages` accepts any
package required for the model function, and `Dyn.libs` accepts dynamic
libraries for parallelization, if required.

Liu (2000) demonstrate the combination of MTM with the Snooker algorithm
from Adaptive Directional Sampling (ADS), Conjugate Gradient Monte Carlo
(CGMC), Hit-And-Run modified as Random-Ray Monte Carlo (RRMC), and
Griddy-Gibbs (GG). MTM has since been extended to multivariate
proposals, proposals with different scales, and more.

Advantages of MTM over Metropolis-within-Gibbs (MWG) is that the
acceptance rate is higher, and multiple evaluations of the model
specification function are parallelized each iteration. The advantage of
MTM over Griddy-Gibbs (GG) is exact rather than approximate estimation
and that an equilibrium distribution cannot be guaranteed from an
approximation of the conditional such as in GG. A disadvantage of MTM
compared to MWG is that MTM must evaluate the model specification
function multiple times per parameter per iteration, resulting in an
algorithm that is slower per iteration. Since MTM is not adaptive, it is
suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, number of tries \\K\\,
    proposal variance \\\sigma^2\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\:
        - Generate \\K\\ proposals: \\\theta^{(k)}\_j \sim
          \mathcal{N}(\theta\_{t,j}, \sigma^2)\\ for \\k = 1, \ldots,
          K\\
        - Compute weights \\w_k = p(\theta^{(k)} \mid y)\\ for each
          proposal
        - Select one proposal \\\theta^\*\_j\\ with probability
          \\\propto w_k\\
        - Generate reference set: \\\tilde{\theta}^{(k)}\_j \sim
          \mathcal{N}(\theta^\*\_j, \sigma^2)\\ for \\k = 1, \ldots,
          K-1\\; set \\\tilde{\theta}^{(K)}\_j = \theta\_{t,j}\\
        - Compute \\\alpha = \min\\\bigl(1,\\ \sum w_k / \sum
          \tilde{w}\_k\bigr)\\
        - With probability \\\alpha\\, accept \\\theta\_{t+1,j} =
          \theta^\*\_j\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| K | 4 | Four trial proposals per iteration: enough to materially raise the acceptance rate without quadrupling the cost per iteration. |
| CPUs | 1 | Sequential evaluation of the K proposals. |
| Packages | NULL | No additional packages needed for the base example. |
| Dyn.libs | NULL | No dynamic libraries to attach. |

Default `Specs` for MTM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MTM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [GG](#gg)
- [MWG](#mwg)

### Interchain Adaptation

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm for parallel chains only. It is best suited to models with a small to medium number of parameters. The proposal covariance matrix must be solved, and this matrix grows with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications, and these are easy to specify. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The RWM algorithm is recommended as the final algorithm, though INCA may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Multivariate. |

The Interchain Adaptation (INCA) algorithm of [\[12\]](#ref12) is an
extension of the Adaptive Metropolis (AM) algorithm of [\[2\]](#ref2).
[\[12\]](#ref12) refer to INCA as inter-chain adaptation and inter-chain
adaptive MCMC. INCA uses parallel chains that are independent, except
that they share the adaptive component, and this sharing speeds
convergence. Since parallel chains are a defining feature of INCA, this
algorithm requires `Chains > 1` and `CPUs > 1` in
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

INCA has two algorithm specifications: `Adaptive` is the iteration in
which it becomes adaptive, and `Periodicity` is the frequency in
iterations of adaptation.

As with AM, the proposal covariance matrix is set equal to the
historical (or sample) covariance matrix. Ample pre-adaptive iterations
are recommended [\[12\]](#ref12), and initial values should be dispersed
to aid in discovering multimodal marginal posterior distributions. After
adaptation begins, INCA combines the historical covariance matrices from
all parallel chains during each adaptation. Each chain learns from
experience as in AM, and in INCA, each chain also learns from the other
parallel chains.

[\[13\]](#ref13) found a dramatic reduction in the number of iterations
to convergence when INCA used 10 parallel chains, compared against a
single-chain AM algorithm. Similar improvements have been noted in the
lucifer package, though more time is required per iteration.

The Gelman.Diagnostic is recommended by [\[12\]](#ref12) to determine
when the parallel chains have stopped sharing different information
about the target distributions. The exchange of information between
chains decreases as the number of iterations increases.

This implementation of INCA uses a server function that is built into
lucifer.hpc. If the connection to this server fails, then the user must
kill the process and then close all open connections with the
closeAllConnections function.

Since INCA is an adaptive algorithm, the final algorithm should be
Random-Walk Metropolis (RWM).

#### Algorithm

1.  **Initialize** \\M\\ parallel chains with parameters
    \\\\\theta^{(1)}\_0, \ldots, \theta^{(M)}\_0\\\\, individual
    covariances \\\Sigma^{(m)}\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each chain \\m = 1, \ldots, M\\ (independently):
        - Propose \\\theta^{(m)\*} \sim \mathcal{N}(\theta^{(m)}\_{t-1},
          \Sigma^{(m)})\\
        - Accept or reject via standard Metropolis step
    2.  If \\t \geq\\ `Adaptive` and \\t \bmod\\ `Periodicity` \\= 0\\:
        - Each chain updates its own \\\Sigma^{(m)}\\ from its sample
          history (as in AM)
        - **Inter-chain sharing**: combine covariance estimates across
          all chains: \\\Sigma\_{\text{shared}} = \frac{1}{M}
          \sum\_{m=1}^M \Sigma^{(m)}\\
        - Set \\\Sigma^{(m)} \leftarrow \Sigma\_{\text{shared}}\\ for
          all chains
3.  **Return** samples from all chains

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Adaptive | floor(Iterations/2) (= 1000) | Defer adaptation until enough samples are available across chains. |
| Periodicity | 1 | Re-estimate the joint proposal every iteration; INCA pools information across chains so the cost is amortized. |

Default `Specs` for INCA. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "INCA", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AM](#am)
- [RWM](#rwm)

## Componentwise and Gibbs samplers update parameters one at a time

Componentwise methods update one parameter (or one block) per
sub-iteration while holding the remaining parameters fixed, cycling
through all parameters to complete a single full iteration. Pure Gibbs
sampling draws each parameter from its exact full conditional
distribution, which requires conjugacy; Metropolis-within-Gibbs relaxes
this requirement by using a Metropolis step for each component. Adaptive
variants tune the componentwise proposal scales or directions
automatically, while sequential and updating-sequential variants control
the schedule and scope of adaptation.

### Gibbs Sampler

> *New automatic full conditionals generation in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm that generally requires conjugacy. |
| Difficulty | This algorithm is difficult for a beginner when the conditional distribution must be specified. Otherwise, it is fully automatic. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

Gibbs sampling was introduced by [\[16\]](#ref16), and later by brothers
[\[14\]](#ref14). The Geman brothers named the algorithm after the
physicist J. W. Gibbs, some eight decades after his death, in reference
to an analogy between the sampling algorithm and statistical physics.
Geman and Geman introduced Gibbs sampling in the context of image
restoration.

Gibbs has two algorithm specifications: `FC` accepts a user-specified
function to calculate full-conditionals, and `MWG` defaults to `NULL`
and otherwise accepts a numeric vector that indicates which parameters
are updated with Metropolis-within-Gibbs. FC accepts two arguments,
`parm` and `Data`, just like the Model specification function, and
returns the full vector of parameters `parm`. Each iteration, the
full-conditional distributions are completed first, then MWG updates, if
any.

In its basic version, Gibbs sampling is a special case of the
Metropolis-Hastings algorithm. However, in its extended versions, Gibbs
sampling can be considered a general framework for sampling from a large
set of variables by sampling each variable (or in some cases, each group
of variables) in turn, and can incorporate the Metropolis-Hastings
algorithm (such as Metropolis-within-Gibbs or similar methods such as
Slice sampling) to implement one or more of the sampling steps.
Currently, Metropolis-within-Gibbs is included.

Gibbs sampling is applicable when the joint distribution is not known
explicitly or is difficult to sample from directly, but the conditional
distribution of each variable is known and easy to sample from. A Gibbs
sampler generates a draw from the distribution of each parameter or
variable in turn, conditional on the current values of the other
parameters or variables. Therefore, a Gibbs sampler is a componentwise
algorithm. As a simplified example, given a joint probability
distribution \\p(θ_1, θ_2 \| y)\\, a Gibbs sampler would draw \\p(θ_1 \|
θ_2, y)\\, then \\p(θ_2 \| θ_1, y)\\.

For a user to determine each conditional distribution, the joint
distribution must be known first, then all parameters are held constant
except the current parameter of interest, and the equation is simplified
to the conditional distribution.

There are numerous variations of Gibbs sampling, such as blocked Gibbs
sampling, collapsed Gibbs sampling, and ordered overrelaxation.

The advantage of Gibbs sampling to other MCMC algorithms is that it is
often more efficient when it is appropriate, due to a 100% acceptance
rate. The disadvantages are that a Gibbs sampler is appropriate only
with conjugate distributions and low correlations between parameters,
and therefore Gibbs sampling is less generally applicable than other
MCMC algorithms. Since Gibbs is not adaptive, it is suitable as a final
algorithm.

#### Automatic full conditionals

When `Specs` is `NULL` (or omitted), the Gibbs sampler automatically
constructs full-conditional draws using a mixed strategy that handles
both continuous and discrete parameters. Continuous parameters are
sampled via stepping-out slice sampling [\[15\]](#ref15), which requires
only log-posterior evaluations and produces exact draws from
\\p(\theta_j \mid \theta\_{-j}, y)\\. Discrete parameters, declared via
`Data$dparm`, are sampled by exhaustive enumeration of their finite
support, computing the log-posterior for each value and drawing from the
resulting categorical distribution. Both methods are exact, and the
acceptance rate is 100%.

The cost is approximately 3–10 Model evaluations per continuous
parameter and \\\|\text{support}\|\\ evaluations per discrete parameter
per iteration. For models where conjugate distributions are available
and speed is critical, providing a manual `FC` function remains
preferable. For all other cases, including models with spike-and-slab
priors (SSVS), auto-FC provides a correct, zero-configuration
alternative.

Usage:

``` r

# Continuous parameters only
fit <- lucifer(Model, Data, IV, Algorithm = "Gibbs")

# With discrete parameters (e.g., Bernoulli indicators for SSVS)
Data$dparm <- pos.Gamma         # integer vector of discrete parameter indices
Data$dsupport <- NULL            # defaults to {0, 1}; or list of support vectors
fit <- lucifer(Model, Data, IV, Algorithm = "Gibbs")
```

Missing data can be handled by treating unobserved values as additional
latent parameters in the model, exactly as in JAGS or BUGS. Auto-FC
samples these imputed values via slice sampling alongside the model
parameters. See
[`vignette("gibbs", package = "lucifer")`](https://robustecologies.github.io/lucifer/articles/gibbs.md)
for a comprehensive treatment including theory, examples, and manual FC
derivations.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Full-conditional updates** (via user-supplied `FC` function):
        - **For** each parameter \\j\\ with a known conditional
          distribution:
          - Draw \\\theta\_{t+1,j} \sim p(\theta_j \mid \theta\_{-j},
            y)\\ (the full conditional)
    2.  **MWG updates** (for parameters without conjugate conditionals):
        - **For** each parameter \\j\\ in the `MWG` set:
          - Propose \\\theta^\*\_j \sim q(\theta^\*\_j \mid
            \theta\_{t,j})\\
          - Accept or reject via Metropolis step
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| FC | “auto” | When `Specs = NULL` lucifer constructs full conditionals automatically: continuous parameters are updated via univariate slice sampling (Neal 2003) and discrete parameters via exhaustive enumeration over `Data$dparm`. Replace with a user function `function(parm, Data)` for analytic Gibbs steps. |
| MWG | NULL | No Metropolis-within-Gibbs fallback positions are needed when FC is automatic. |

Default `Specs` for Gibbs. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Gibbs", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [GG](#gg)
- [MWG](#mwg)
- [Slice](#slice)

### Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner when the proposal variance has been tuned with another algorithm. Otherwise, it may be tedious for the user to tune the proposal variance. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

Metropolis-within-Gibbs (MWG) is the original MCMC algorithm, introduced
in [\[18\]](#ref18). Since it was the original MCMC algorithm, it
pre-dated Gibbs sampling (Gibbs), and was not known as
Metropolis-within-Gibbs. MWG was later proposed as a hybrid algorithm
that combines Metropolis-Hastings and Gibbs sampling, and was suggested
in [\[10\]](#ref10). The idea was to substitute a Metropolis step when
Gibbs sampling fails. MWG is available in this respect in the Gibbs
sampling algorithm. However, Gibbs sampling is not included in this
particular version in lucifer (or most versions), making it an algorithm
with a misleading name. Without Gibbs sampling, the more honest name
would be componentwise random-walk Metropolis, but the more common name
is MWG. MWG is also referred to as Metropolis within Gibbs,
Metropolis-in-Gibbs, Random-Walk Metropolis-within-Gibbs, single-site
Metropolis, the Metropolis algorithm, or Variable-at-a-Time Metropolis.

##### Componentwise

MWG is a componentwise algorithm, meaning that each parameter is updated
individually each iteration. This implies that the model specification
function is evaluated a number of times equal to the number of
parameters, per iteration. Componentwise sampling usually ignores
parameter correlation. The order of the parameters for updating is
randomized each iteration (random-scan MWG), as opposed to sequential
updating (deterministic-scan MWG). MWG often uses blocks, and blockwise
sampling may be specified with the algorithm specification `B`, in which
a list is supplied, and each list component consists of a vector of
parameters. Blocks are updated sequentially, but the order of parameters
is randomized within each block.

##### Metropolis Step

The update of each parameter occurs in a Metropolis step, otherwise
called an accept/reject step or a Metropolis accept/reject step. A
componentwise proposal is generated randomly and the model is evaluated
with the proposed parameter. If the proposal is an improvement in the
logarithm of the unnormalized joint posterior density, then the proposal
is accepted. If the proposal is not an improvement, then the proposal is
accepted or rejected according to a probability.

##### Acceptance Rate

The acceptance rate is the total number of proposals accepted in the
Metropolis step divided by the total number of iterations. If the
acceptance rate is too high, then the proposal variance is too small. In
this case, the algorithm will take longer than necessary to find the
target distribution and the samples will be highly autocorrelated. If
the acceptance rate is too low, then the proposal variance is too large,
and the algorithm is ineffective at exploration. In the worst case
scenario, no proposals are accepted and the algorithm fails to move.
Under theoretical conditions, the optimal acceptance rate for a sole,
independent and identically distributed (IID), Gaussian, marginal
posterior distribution is 0.44 or 44%. The optimal acceptance rate for
an infinite number of distributions that are IID and Gaussian is 0.234
or 23.4%. Since MWG is a componentwise algorithm, it is most efficient
when the acceptance rate of each parameter is 0.44.

##### Random-Walk Behavior

MWG is a componentwise Random-Walk Metropolis (RWM) algorithm.
Random-walk behavior is desirable because it allows the algorithm to
explore, and hopefully avoid getting trapped in undesirable regions. On
the other hand, random-walk behavior is undesirable because it takes
longer to converge to the target distribution while the algorithm
explores. The algorithm generally progresses in the right direction, but
may periodically wander away. Such exploration may uncover multimodal
target distributions, which other algorithms may fail to recognize, and
then converge incorrectly. With enough iterations, MWG is guaranteed
theoretically to converge to the correct target distribution, regardless
of the starting point of each parameter, provided the proposal variance
for each proposal of a target distribution is sensible.

Historically, MWG first ran on an early computer that was built
specifically for it, called MANIAC I. Metropolis discarded the first 16
iterations as burn-in, and updated an additional 48-64 iterations, which
required 4-5 hours on MANIAC I.

[\[19\]](#ref19) introduced an adaptive version of
Metropolis-within-Gibbs, called Adaptive Metropolis-within-Gibbs (AMWG).
[\[17\]](#ref17) extended this to the Adaptive Directional
Metropolis-within-Gibbs (ADMG) algorithm.

The advantage of MWG over the multivariate version, RWM, is that it is
more efficient with information per iteration, so convergence is faster
in iterations. The disadvantages of MWG are that covariance is not
included in proposals, and it is more time-consuming due to the
evaluation of the model specification function for each parameter per
iteration. As the number of parameters increases, and especially as
model complexity increases, the run-time per iteration increases. Since
fewer iterations are completed in a given time-interval, the possible
amount of thinning is also at a disadvantage. MWG is non-adaptive, and
is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal standard deviations
    \\\sigma_j\\ for \\j = 1, \ldots, J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Propose \\\theta^\*\_j \sim \mathcal{N}(\theta\_{t,j},
          \sigma_j^2)\\, holding all other parameters fixed
        - Evaluate model: compute \\\log p(\theta^\* \mid y)\\
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta_t \mid y)\bigr)\\
        - Draw \\u \sim \text{Uniform}(0,1)\\; if \\u \< \alpha\\, set
          \\\theta\_{t+1,j} = \theta^\*\_j\\; otherwise
          \\\theta\_{t+1,j} = \theta\_{t,j}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| B | NULL | Componentwise update over each parameter; supply a list of index vectors to enable group updates. |

Default `Specs` for MWG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MWG", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [ADMG](#admg)
- [AMWG](#amwg)
- [Gibbs](#gibbs)
- [RWM](#rwm)

### Adaptive Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is relatively easy for a beginner. It has three algorithm specifications, which are easy to use. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The MWG algorithm is recommended as the final algorithm, though AMWG may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Componentwise. |

The Adaptive Metropolis-within-Gibbs (AMWG) algorithm is presented in
(Roberts and Rosenthal 2009; Rosenthal 2007). It is an adaptive version
of Metropolis-within-Gibbs (MWG).

AMWG has three algorithm specifications: `B` is optional and allows
blockwise sampling. It defaults to NULL, but may accept a list in which
each list component is a block of parameters. If blockwise sampling is
used, then blocks are updated sequentially, and the order of parameters
within blocks is randomized. The `n` specification defaults to 0 and is
used to keep track of how many previous adaptive iterations were run.
The `Periodicity` specification indicates the frequency in iterations at
which the algorithm adapts. If Periodicity is set too low, such as
`Periodicity=1`, the algorithm may adapt too quickly to poor information
in the beginning, and become unstable. `Periodicity=50` is recommended.

In AMWG, the standard deviation of the proposal of each parameter is
manipulated to optimize the associated acceptance rate toward 0.44. This
is much simpler than other adaptive methods that adapt based on sample
covariance in large dimensions. Large covariance matrices require a
large number of elements to adapt, which takes exponentially longer to
adapt as the dimension increases. Regardless of dimension, the AMWG
optimizes each parameter to a univariate acceptance rate, and a sample
covariance matrix does not need to be estimated for adaptation, which
consumes time and memory. The order of the parameters for updating is
randomized each iteration (random-scan AMWG), as opposed to sequential
updating (deterministic-scan AMWG).

Compared to other adaptive algorithms with multivariate proposals, a
disadvantage is the time to complete each iteration increases as a
function of parameters and model complexity, as noted in MWG. For
example, in a 100-parameter model, AMWG completes its first iteration as
the Adaptive-Mixture Metropolis (AMM) algorithm completes its 100^(th).
However, to adapt accurately, the AMM algorithm must correctly estimate
5,050 elements of a sample covariance matrix, while AMWG must correctly
estimate only 100 proposal standard deviations. [\[3\]](#ref3) have
shown an example model with 500 parameters that had a burn-in of around
25,000 iterations.

The advantages of AMWG over AMM are that AMWG does not require a burn-in
period before it can begin to adapt, and that AMWG does not need to
estimate a covariance matrix to adapt properly. The disadvantages of
AMWG compared to AMM are that correlation can be problematic since it is
not taken into account with a proposal covariance matrix, and AMWG
solves the model function once per parameter per iteration, which can be
unacceptably slow with large or complicated models. The advantage of
AMWG over Robust Adaptive Metropolis (RAM) is that AMWG does not need to
estimate a covariance matrix to adapt properly. The disadvantages of
AMWG compared to RAM are AMWG is less likely to handle multimodal or
heavy-tailed targets, and AMWG solves the model function once per
parameter per iteration, which can be unacceptably slow with large or
complicated models. If AMWG is used for adaptation, then the final,
non-adaptive algorithm should be MWG.

[\[17\]](#ref17) extended AMWG to the Adaptive Directional
Metropolis-within-Gibbs (ADMG) algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, log-scale proposal standard
    deviations \\\log \sigma_j\\ for \\j = 1, \ldots, J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Propose \\\theta^\*\_j \sim \mathcal{N}(\theta\_{t,j}, \exp(2
          \log \sigma_j))\\
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta_t \mid y)\bigr)\\
        - With probability \\\alpha\\, accept \\\theta\_{t+1,j} =
          \theta^\*\_j\\
    2.  If \\t \bmod\\ `Periodicity` \\= 0\\ and \\t \> n\\, **adapt**
        each \\\sigma_j\\:
        - If acceptance rate for parameter \\j \> 0.44\\: \\\log
          \sigma_j \leftarrow \log \sigma_j + \delta(t)\\
        - If acceptance rate for parameter \\j \< 0.44\\: \\\log
          \sigma_j \leftarrow \log \sigma_j - \delta(t)\\
        - where \\\delta(t) \to 0\\ ensures diminishing adaptation
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| B | NULL | Componentwise updates by default. |
| n | 0 | Adaptation history starts empty. |
| Periodicity | 50 | Variance adaptation is performed every 50 iterations: frequent enough to converge during warm-up, rare enough that the variances stabilize and diminishing adaptation is preserved. |

Default `Specs` for AMWG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AMWG", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [ADMG](#admg)
- [MWG](#mwg)

### Adaptive Directional Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is relatively easy for a beginner. It has only two algorithm specifications. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The MWG algorithm is recommended as the final algorithm, though ADMG may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Componentwise. |

The Adaptive Directional Metropolis-within-Gibbs (ADMG) algorithm was
proposed by [\[17\]](#ref17) as an extension of the Adaptive
Metropolis-within-Gibbs (AMWG) algorithm in which the direction and
jumping distance of the componentwise proposal is affected by an
estimate of the historical covariance matrix.

ADMG has two algorithm specifications: `n` is the number of previous
iterations and `Periodicity` is the frequency in iterations of
adaptation. The minimum Periodicity is equal to the number of
parameters, and should probably be greater.

Although ADMG is a componentwise algorithm, it does not ignore
correlation of parameters, as is typical with componentwise algorithms.
ADMG makes each componentwise proposal based on an estimate of the
historical covariance matrix. Adaptation is based on the singular value
decomposition (SVD) of the estimate of the historical covariance matrix.
In large dimensions, SVD is computationally expensive, and larger
integers are recommended for the algorithm specifications. If SVD
results in a singularity, then the algorithm defaults to a simple
Metropolis-within-Gibbs MWG proposal. Otherwise, when SVD is successful,
then MWG is performed under the rotated coordinates.

The author notes that the MWG part of the algorithm may proceed by
deterministic-scan (which he refers to as system-scan) and random-scan.
The algorithm implemented here uses random-scan, and is denoted by the
author as ADRSMG. Random-scan updates means that the order of the update
of the parameters is randomized each iteration.

Compared to other adaptive algorithms with multivariate proposals, a
disadvantage is the time to complete each iteration increases as a
function of parameters and model complexity, as noted in MWG. For
example, in a 100-parameter model, ADMG completes its first iteration as
the AMM algorithm completes its 100^(th). However, ADMG is more
efficient with information per iteration.

The advantage of ADMG to other componentwise algorithms (such as AMWG,
MWG, RDMH, and Slice) is that ADMG makes each proposal while accounting
for parameter correlation, while other componentwise algorithms
traditionally ignore parameter correlation. The author asserts that
SVD-based adaptation is more efficient than adapting directly to the
historical covariance matrix or an estimate thereof. The disadvantage of
ADMG to other componentwise algorithms is that SVD becomes
computationally expensive in large dimensions. Since ADMG is adaptive,
MWG should be used as the final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, proposal standard deviations
    \\\sigma_j\\ for \\j = 1, \ldots, J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  If \\t \bmod\\ `Periodicity` \\= 0\\ and \\t \> n\\, compute SVD
        of the estimated historical covariance \\\hat{\Sigma}\\
    2.  **For** each parameter \\j\\ in random order:
        - If SVD succeeded, propose \\\theta^\*\_j\\ along the rotated
          coordinate axis from the SVD basis
        - Otherwise, fall back to standard MWG: \\\theta^\*\_j \sim
          \mathcal{N}(\theta\_{t,j}, \sigma_j^2)\\
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta_t \mid y)\bigr)\\
        - With probability \\\alpha\\, set \\\theta\_{t+1,j} =
          \theta^\*\_j\\; otherwise retain \\\theta\_{t,j}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| n | 0 | Adaptation history starts empty. |
| Periodicity | 1 | Re-estimate the directional proposal every iteration. |

Default `Specs` for ADMG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "ADMG", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMWG](#amwg)
- [MWG](#mwg)

### Griddy-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is moderate in difficulty because the user-specified grid may need tuning. With an appropriate grid, the method is fully automatic. There is only one required algorithm specification for defining the grid. Additional optional algorithm specifications control parallelization or allow discrete parameters. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

Introduced by [\[20\]](#ref20), the Griddy-Gibbs (GG) sampler is a
componentwise algorithm that approximates the conditional distribution
by evaluating the model at a discrete set of points, the user-specified
grid for each parameter. The proposal is a random draw from the
approximated distribution. In this implementation, splinetic
interpolation is used to approximate the distribution for continuous
parameters with 1,000 points, given the evaluated points. The acceptance
rate is 100%.

GG has five algorithm specifications: `Grid` is a vector or list of
evenly-spaced points in the grid, `dparm` is a vector that indicates
discrete parameters and defaults to NULL, `CPUs` is the number of
available central processing units (CPUs), `Packages` is a vector of
package names, and `Dyn.lib` is a vector of shared libraries. The
default for `Grid` is `seq(from=-0.1, to=0.1, len=5)`, which creates a
grid with the values -0.1, -0.05, 0, 0.05, and 0.1. For each continuous
parameter in each iteration, the grid values are added to the latest
value of the parameter, and the model is evaluated with the parameter at
each point on the grid. For each discrete parameter, the model is
evaluated at each grid point and a sample is taken.

At least five grid points are recommended for a continuous parameter,
and a grid with more points will better approximate the distribution,
but requires more evaluations. However, the approximation in GG may be
parallelized (within lucifer, not lucifer.hpc), so a large computer
environment may approach an excellent approximation with little
inefficiency. It is natural to desire a grid with a larger range, but in
practice this becomes problematic, so it is recommended to keep the
range of the grid relatively small, say within \[-0.1,0.1\] or
\[-0.2,0.2\], and may require experimentation. After observing a Markov
chain, the user may adjust the range of the grid to decrease
autocorrelation in a future update.

When an odd number of grid points is used for a continuous parameter,
the current position is evaluated. If there are too few grid points,
then the current point may be the only point with appreciable density,
and the parameter may never move. This can be checked afterward with the
AcceptanceRate function. If the acceptance rate is less than one for any
parameter, then there are too few grid points. This failure may be
harder to find when there are numerous parameters, an even number of
grid points, and too few grid points, because the acceptance rate may be
100% for a parameter, yet it may be oscillating between two values.

GG is one of the slower algorithms per iteration, since the model must
be evaluated multiple times per parameter per iteration. This may mostly
be alleviated in a parallel environment. GG seems appropriate when the
problem must be solved with a componentwise algorithm, and excessive
thinning is required. GG may help reduce thinning by making proposals
from the approximated conditional distribution, though parameter
correlation may increase autocorrelation.

Advantages of parallelized GG over most componentwise algorithms is that
it yields a 100% acceptance rate, and draws from the approximated
distribution should be less autocorrelated. A disadvantage is that more
model evaluations are required, and even if a parallel environment had
zero overhead, GG would be twice as slow per iteration as other
componentwise algorithms. However, if the user tunes the grid, it may be
more efficient than other componentwise algorithms. The Adaptive
Griddy-Gibbs (AGG) sampler is available so the user can avoid tuning the
grid. Since GG is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, grid \\G\\ of evenly-spaced
    points
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Construct evaluation grid: \\G_j = \theta\_{t,j} + G\\
        - Evaluate \\p(\theta \mid \text{data})\\ at each grid point \\g
          \in G_j\\, holding other parameters fixed
        - For continuous parameters, interpolate the conditional
          distribution via cubic splines over 1000 points
        - For discrete parameters, normalize the evaluations directly as
          probabilities
        - Draw \\\theta\_{t+1,j}\\ from the (interpolated) conditional
          distribution
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Grid | seq(-3, 3, length = 11) | Mandatory: an evenly spaced grid centered at zero with eleven nodes is wide enough to cover the marginal posterior of `beta0`, `beta1` and `log.sigma` in the example. |
| dparm | NULL | All parameters treated as continuous in this example. |
| CPUs | 1 | Sequential evaluation of the conditional density on the grid. |
| Packages | NULL | No extra packages. |
| Dyn.libs | NULL | No dynamic libraries. |

Default `Specs` for GG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "GG", Specs = list(Grid = seq(-3, 3, length = 11), dparm = NULL,
             CPUs = 1L, Packages = NULL, Dyn.libs = NULL),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AGG](#agg)
- [Gibbs](#gibbs)
- [MTM](#mtm)
- [Slice](#slice)

### Adaptive Griddy-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy to use because the scale of the grid is adaptively tuned. There is only one required algorithm specification for defining the grid. Additional optional algorithm specifications control parallelization, allow discrete parameters, or prevent extreme values. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Adaptive Griddy-Gibbs (AGG) sampler is an extension of the
Griddy-Gibbs (GG) sampler in which the scale of the centered grid is
tuned each iteration for continuous parameters. The scale is calculated
as the standard deviation of the conditional distribution. If discrete
parameters are used, then the grid of discrete parameters is not tuned,
since all discrete values are evaluated.

AGG has six algorithm specifications: `Grid` is a vector or list of
points in the grid, `dparm` is a vector that indicates discrete
parameters and defaults to NULL, `smax` is the maximum allowable
conditional standard deviation, `CPUs` is the number of available
central processing units (CPUs), `Packages` is a vector of package
names, and `Dyn.lib` is a vector of shared libraries. The default for
Grid is GaussHermiteQuadRule(3)\$nodes. For each continuous parameter in
each iteration, the scaled grid values are added to the latest value of
the parameter, and the model is evaluated with the parameter at each
point on the grid and a sample is taken. For each discrete parameter,
the model is evaluated at each grid point and a sample is taken.

The points in the grid for continuous parameters are selected by the
user as the nodes of a Gauss-Hermite quadrature rule (with the
GaussHermiteQuadRule function), and these points are not evenly-spaced.
Many problems require as few as 3 nodes, while others require perhaps as
many as 99 nodes. When the number of nodes is larger than necessary,
time is wasted in computation, fewer points occur within most of the
density, and occasional extreme values of LP are observed, but the
chains reach the target distributions in fewer iterations due to the
wider grid. When the number of nodes is too small, it may not converge
on the correct distribution. It is therefore suggested to begin with a
small, odd number of grid points, such as 3, to reduce computation time.
An odd number of grid points is preferred. If occasional extreme values
of LP are observed, set `smax` to something reasonable. If `smax` is too
small, then higher autocorrelation will occur and more thinning is
necessary.

Advantages of parallelized AGG over most componentwise algorithms is
that it yields a 100% acceptance rate, and draws from the approximated
distribution should be less autocorrelated. A disadvantage is that more
model evaluations are required, and even if a parallel environment had
zero overhead, AGG would be twice as slow per iteration as other
componentwise algorithms. AGG is adaptive in the sense of
self-adjusting, but not in the sense of being non-Markovian, so it is
suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, grid \\G\\ of Gauss-Hermite
    quadrature nodes
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Estimate conditional standard deviation \\\hat{\sigma}\_j\\
          from the current state
        - Scale grid: \\G_j^{(t)} = \theta\_{t,j} + \hat{\sigma}\_j
          \times G\\
        - Evaluate \\p(\theta \mid \text{data})\\ at each grid point \\g
          \in G_j^{(t)}\\, holding other parameters fixed
        - Interpolate the conditional CDF via cubic splines over 1000
          points
        - Draw \\\theta\_{t+1,j}\\ from the interpolated conditional
          distribution via inverse CDF sampling
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Grid | GaussHermiteQuadRule(3)\$nodes | Mandatory: a 3-node Gauss–Hermite rule produces an adaptive grid that follows the local curvature. |
| dparm | NULL | All parameters continuous. |
| smax | 5 | Upper bound on the adaptive scale parameter sigma; large enough to absorb prior dispersion without diverging. |
| CPUs | 1 | Sequential evaluation. |
| Packages | NULL | No extra packages. |
| Dyn.libs | NULL | No dynamic libraries. |

Default `Specs` for AGG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AGG", Specs = list(Grid = GaussHermiteQuadRule(3)$nodes, dparm = NULL,
             smax = 5, CPUs = 1L, Packages = NULL, Dyn.libs = NULL),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [GG](#gg)
- [Slice](#slice)

### Sequential Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This algorithm is applicable with state-space models (SSMs), including dynamic linear models (DLMs). |
| Difficulty | This algorithm is relatively easy for a beginner when the proposal variance has been tuned with the SAMWG algorithm. Otherwise, it may be tedious for the user to tune the proposal variance. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Sequential Metropolis-within-Gibbs (SMWG) algorithm is the
non-adaptive version of the Sequential Adaptive Metropolis-within-Gibbs
(SAMWG) algorithm, and is used for final sampling of state-space models
(SSMs).

#### Algorithm

1.  **Initialize** parameters \\\theta_0 = (\theta\_{\text{static}},
    \theta\_{\text{dynamic},1:T_s})\\, tuned proposal SDs \\\sigma_j\\
    (from SAMWG)
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Update static parameters** (random order): for each static
        \\j\\, propose \\\theta^\*\_j \sim \mathcal{N}(\theta\_{t,j},
        \sigma_j^2)\\; accept/reject via Metropolis step
    2.  **Update dynamic parameters sequentially** through time: for
        each \\s = 1, \ldots, T_s\\ and dynamic parameter \\j\\, propose
        and accept/reject (no adaptation)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Dyn | matrix(parm.names, nrow = 1) | Mandatory: a single static block over all parameters is the safest default when no temporal blocking is required. |

Default `Specs` for SMWG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "SMWG", Specs = list(Dyn = matrix(Data$parm.names, nrow = 1)),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [SAMWG](#samwg)

### Sequential Adaptive Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This algorithm is applicable with state-space models (SSMs), including dynamic linear models (DLMs). |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The SMWG algorithm is recommended as the final algorithm, though SAMWG may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Componentwise. |

The Sequential Adaptive Metropolis-within-Gibbs (SAMWG) algorithm is for
state-space models (SSMs), including dynamic linear models (DLMs). It is
identical to the Adaptive Metropolis-within-Gibbs (AMWG) algorithm,
except with regard to the order of updating parameters (and here,
sequential does not refer to deterministic-scan). Parameters are grouped
into two blocks: static and dynamic. At each iteration, static
parameters are updated first, followed by dynamic parameters, which are
updated sequentially through the time-periods of the model. The order of
the static parameters is randomly selected at each iteration, and if
there are multiple dynamic parameters for each time-period, then the
order of the dynamic parameters is also randomly selected. The SAMWG
algorithm is adapted from [\[21\]](#ref21) for lucifer. The SAMWG is a
single-site update algorithm that is more efficient in terms of
iterations, though convergence can be slow with high intercorrelations
in the state vector [\[22\]](#ref22). If SAMWG is used for adaptation,
then the final, non-adaptive algorithm should be Sequential
Metropolis-within-Gibbs (SMWG).

#### Algorithm

1.  **Initialize** parameters \\\theta_0 = (\theta\_{\text{static}},
    \theta\_{\text{dynamic},1:T_s})\\, proposal SDs \\\sigma_j\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Update static parameters** (random order): for each static
        parameter \\j\\, propose \\\theta^\*\_j \sim
        \mathcal{N}(\theta\_{t,j}, \sigma_j^2)\\; accept/reject via
        Metropolis step
    2.  **Update dynamic parameters sequentially** through time periods
        \\s = 1, \ldots, T_s\\: for each dynamic parameter \\j\\ at time
        \\s\\ (random order within time), propose and accept/reject via
        Metropolis step
    3.  If \\t \bmod\\ `Periodicity` \\= 0\\, adapt each \\\sigma_j\\
        toward 44% acceptance (as in AMWG)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Dyn | matrix(parm.names, nrow = 1) | Mandatory: a single block. |
| Periodicity | 50 | Variance adaptation every 50 iterations. |

Default `Specs` for SAMWG. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "SAMWG", Specs = list(Dyn = matrix(Data$parm.names, nrow = 1), Periodicity = 50L),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMWG](#amwg)
- [SMWG](#smwg)

### Updating Sequential Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This algorithm is applicable with state-space models (SSMs), including dynamic linear models (DLMs). |
| Difficulty | This algorithm is relatively easy for a beginner when the proposal variance has been tuned with the USAMWG algorithm. Otherwise, it may be tedious for the user to tune the proposal variance. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Updating Sequential Metropolis-within-Gibbs (USMWG) algorithm is the
non-adaptive version of the USAMWG algorithm, and is used for final
sampling when updating state-space models (SSMs).

For example, suppose a time-series of daily values was fit with the
SAMWG algorithm, and SMWG may have been used for final samples. Let’s
also suppose the last day in the model sample was a Monday, and that
one-step ahead forecasts are produced, so the model predicted Tuesday.
After the actual value for Tuesday is known, it may be included in the
model. Using USAMWG, the latest Tuesday is filtered and Wednesday is
forecast, while the days of the original model sample are not estimated
again. USMWG may then be used for final samples.

#### Algorithm

1.  **Given** a previously fitted SSM with tuned proposal SDs
    \\\sigma_j\\ (from USAMWG)
2.  **Receive** new observation \\y\_{T_s+1}\\
3.  **For** \\t = 1, \ldots, T\\:
    1.  **Filter**: update dynamic state \\\theta\_{\text{dyn},T_s+1}\\
        using non-adaptive MWG proposals with fixed \\\sigma_j\\
    2.  **Predict**: forecast \\\theta\_{\text{dyn},T_s+2}\\ from the
        state transition model
4.  **Return** updated parameter samples and one-step-ahead forecasts

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Dyn | (supplied by user) | Mandatory. |
| Fit | (supplied by user) | Mandatory: USMWG continues sampling from a previous `demonoid` fit. |
| Begin | (supplied by user) | Mandatory: time-step at which to resume. |

Default `Specs` for USMWG. {.table}

#### Example with default specifications

``` r

## USMWG requires user-supplied specifications that have no sensible
## defaults. See the table above and the help page ?lucifer for details.
```

#### See Also

- [USAMWG](#usamwg)

### Updating Sequential Adaptive Metropolis-within-Gibbs

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This algorithm is applicable with state-space models (SSMs), including dynamic linear models (DLMs). |
| Difficulty | This algorithm is relatively easy for a beginner. It has few algorithm specifications. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The USMWG algorithm is recommended as the final algorithm, though USAMWG may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Componentwise. |

The Updating Sequential Adaptive Metropolis-within-Gibbs (USAMWG) is for
state-space models (SSMs), including dynamic linear models (DLMs). After
a model is fit with Sequential Adaptive Metropolis-within-Gibbs (SAMWG)
and Sequential Metropolis-within-Gibbs (SMWG), and information is later
obtained regarding the first future state predicted by the model, the
USAMWG algorithm may be applied to update the model given the new
information. In SSM terminology, updating is filtering and predicting.
This is more efficient than re-estimating the entire model each time new
information is obtained.

For example, suppose a time-series of daily values was fit with the
SAMWG algorithm. Let’s also suppose the last day in the model sample was
a Monday, and that one-step ahead forecasts are produced, so the model
predicted Tuesday. After the actual value for Tuesday is known, it may
be included in the model. Using USAMWG, the latest Tuesday is filtered
and Wednesday is forecast, while the days of the original model sample
are not estimated again.

#### Algorithm

1.  **Given** a previously fitted SSM with parameters \\\theta\\ and
    dynamic states \\\theta\_{\text{dyn},1:T_s}\\
2.  **Receive** new observation \\y\_{T_s+1}\\
3.  **For** \\t = 1, \ldots, T\\:
    1.  **Filter**: update dynamic state \\\theta\_{\text{dyn},T_s+1}\\
        given \\y\_{T_s+1}\\ and current static parameters, using
        adaptive MWG proposals
    2.  **Predict**: forecast \\\theta\_{\text{dyn},T_s+2}\\ from the
        state transition model
    3.  Adapt proposal standard deviations \\\sigma_j\\ toward 44%
        acceptance (as in AMWG)
4.  **Return** updated parameter samples and one-step-ahead forecasts

#### Default specifications

| Specification | Default            | Justification |
|:--------------|:-------------------|:--------------|
| Dyn           | (supplied by user) | Mandatory.    |
| Periodicity   | (supplied by user) | Mandatory.    |
| Fit           | (supplied by user) | Mandatory.    |
| Begin         | (supplied by user) | Mandatory.    |

Default `Specs` for USAMWG. {.table}

#### Example with default specifications

``` r

## USAMWG requires user-supplied specifications that have no sensible
## defaults. See the table above and the help page ?lucifer for details.
```

#### See Also

- [SAMWG](#samwg)
- [SMWG](#smwg)

### Random Dive Metropolis-Hastings

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is especiallly useful with multimodal or fat-tailed target distributions. |
| Difficulty | This algorithm is easy for a beginner because it is fully automatic; tuning is not required. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

Random Dive Metropolis-Hastings (RDMH) was introduced by
[\[23\]](#ref23). RDMH is a componentwise Metropolis-Hastings algorithm
in which the proposal is the product or ratio of a current parameter and
a random quantity that does not require any tuning parameters. The
random quantity is in the interval (-1,1). Although RDMH is a very
simple algorithm, it has excellent convergence and mixing properties.
However, if it is reasonable that the origin (0) has positive
probability, then the model should be reparameterized, because RDMH
fails in the obscure case when the origin has positive probability
(which can arise if the state-space is the set of integers).

RDMH does not have any algorithm specifications.

RDMH allows the proposed state to be far away from the current state,
and yet has a good acceptance rate. Therefore, RDMH can explore the
state space faster than many other MCMC algorithms.

RDMH is geometrically ergodic for a class of densities that is much
larger than most other MCMC algorithms. RDMH excels at representing
target densities that are multimodal or fat-tailed, and has been
compared with Gibbs sampling (Gibbs), Metropolis-Adjusted Langevin
Algorithm (MALA), and Metropolis-within-Gibbs (MWG). Each of these other
algorithms became stuck at local modes in multimodal target densities
with large distances between the modes. RDMH explored the multimodal
densities, or the fat tails, with ease.

As a componentwise algorithm, the model specification function is
evaluated a number of times equal to the number of parameters, per
iteration. The order of the parameters for updating is randomized each
iteration (random-scan), as opposed to sequential updating
(deterministic-scan).

The advantages of RDMH over MWG are that RDMH does not require tuning,
better explores multimodal and fat-tailed target distributions, is
better able to find acceptable proposals that are distant, and that it
may therefore explore the state-space faster. Compared to multivariate
MCMC algorithms, RDMH shares common disadvantages that it is slower per
iteration, and correlated parameters are not taken into account. Since
RDMH is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Draw \\v \sim \text{Uniform}(-1, 1)\\
        - If \\v \> 0\\: propose \\\theta^\*\_j = \theta\_{t,j} \cdot
          (1 + v)\\ (multiplicative step)
        - If \\v \< 0\\: propose \\\theta^\*\_j = \theta\_{t,j} / (1 -
          v)\\ (division step)
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta_t \mid y) \cdot \|1/(1-v)\|\bigr)\\ (with Jacobian)
        - With probability \\\alpha\\, accept \\\theta\_{t+1,j} =
          \theta^\*\_j\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| (none) | (none) | RDMH has no algorithm specifications. The dive proposal is generated from `Covar` directly. |

Default `Specs` for RDMH. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RDMH", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Gibbs](#gibbs)
- [MALA](#mala)
- [MWG](#mwg)

### Zanella Locally Balanced

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | Not directly comparable to continuous samplers. Acceptance rates vary widely depending on the discrete support size and balancing function. |
| Applications | This algorithm targets models with discrete parameters (variable selection indicators, change-point locations, hidden Markov model states, network edges). It dramatically improves mixing over exhaustive enumeration when the discrete support is large. |
| Difficulty | This algorithm is easy for a beginner. It has two algorithm specifications. The default balancing function (sqrt) works well in most settings. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

[\[24\]](#ref24) introduced locally balanced proposals for
discrete-space MCMC that are Peskun-optimal: they maximize the
probability of moving to a new state among all reversible kernels with
the same neighborhood structure. In standard discrete Gibbs sampling (as
in the auto-FC mechanism), the sampler enumerates all possible values of
a discrete parameter and samples proportionally to the posterior. This
works well when the support is small but becomes expensive when the
support is large. Zanella’s method instead proposes from a weighted
distribution over the neighborhood \\\mathcal{N}(x)\\ of the current
state, where each neighbor \\y\\ receives weight proportional to a
balancing function \\g\\ applied to the posterior ratio: \\q(y \mid x)
\propto g\\\bigl(\pi(y)/\pi(x)\bigr)\\. The two recommended balancing
functions are \\g(r) = \sqrt{r}\\ (which minimizes asymptotic variance)
and \\g(r) = \min(1, r)\\ (the Barker choice, which is more robust).

Zanella has two algorithm specifications: `balancing_fn` is either
`"sqrt"` (default) or `"min"`, and `B` accepts an optional list of
blocks. For discrete parameters (identified by `Data$dparm` indices,
with support defined in `Data$dsupport`), the locally balanced proposal
is used. For continuous parameters in the same model, the sampler falls
back to standard Metropolis-within-Gibbs with Gaussian proposals. The MH
acceptance probability for the locally balanced proposal corrects for
the asymmetry: \\\alpha = \min\\\bigl(1,\\ Z(x) / Z(y)\bigr)\\, where
\\Z(x) = \sum\_{z \in \mathcal{N}(x)} g\\\bigl(\pi(z)/\pi(x)\bigr)\\ is
the normalizing constant of the proposal at \\x\\.

The advantage of Zanella over exhaustive enumeration is that it requires
fewer model evaluations per update (only the neighborhood, not the full
support) while achieving better mixing. The advantage over uniform
random proposals (standard MH for discrete spaces) is that the locally
balanced weighting concentrates proposals on high-probability neighbors,
dramatically improving acceptance rates in high-dimensional discrete
spaces.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ (random scan):
        - If \\j\\ is discrete with neighborhood
          \\\mathcal{N}(\theta_j)\\:
          1.  For each \\y \in \mathcal{N}(\theta_j)\\: compute \\r_y =
              \pi(\theta\_{-j}, y) / \pi(\theta)\\, weight \\w_y =
              g(r_y)\\
          2.  Sample \\\theta_j^\* \sim \text{Categorical}(\\w_y / \sum
              w\\)\\
          3.  Compute reverse sum \\Z^\* = \sum\_{z \in
              \mathcal{N}(\theta_j^\*)} g\\\bigl(\pi(\theta\_{-j}, z) /
              \pi(\theta\_{-j}, \theta_j^\*)\bigr)\\
          4.  Accept with probability \\\min(1,\\ Z(\theta_j) / Z^\*)\\
        - If \\j\\ is continuous: standard MWG Gaussian proposal
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| balancing_fn | “sqrt” | Square-root locally balanced function: the choice that gave the best mixing in Zanella (2020). |
| B | NULL | Single-block update. |

Default `Specs` for Zanella. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Zanella", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MWG](#mwg)
- [Gibbs](#gibbs)
- [AMWG](#amwg)

### Polya-Gamma Data Augmentation

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 100% because this is a Gibbs sampler with exact full conditionals. |
| Applications | This algorithm is specifically designed for logistic regression, binomial regression, and negative-binomial regression models with a logit link. It requires the design matrix Data\\X and binary/count response Data\\y. |
| Difficulty | This algorithm is easy for a beginner. The algorithm specifications (prior mean and variance) have sensible defaults. The model must have the logistic/binomial structure. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise (block Gibbs). |

Polya-Gamma data augmentation was introduced by [\[26\]](#ref26) as a
method for exact Bayesian inference in logistic and binomial regression
models without tuning parameters. The key insight is that the logistic
likelihood \\\prod_i \[e^{\psi_i}/(1 + e^{\psi_i})\]^{y_i} \[1/(1 +
e^{\psi_i})\]^{1-y_i}\\, where \\\psi_i = x_i^\top \beta\\, can be
represented as a Gaussian scale mixture by introducing auxiliary
Polya-Gamma random variables \\\omega_i \sim \text{PG}(b_i, \psi_i)\\.
Conditional on \\\omega\\, the posterior for \\\beta\\ is multivariate
Gaussian: \\\beta \mid \omega, y \sim \mathcal{N}(m\_\omega,
V\_\omega)\\, where \\V\_\omega = (X^\top \Omega X +
\Sigma_0^{-1})^{-1}\\ and \\m\_\omega = V\_\omega(X^\top \kappa +
\Sigma_0^{-1}\mu_0)\\, with \\\Omega = \text{diag}(\omega)\\ and
\\\kappa_i = y_i - b_i/2\\.

PG has three algorithm specifications: `prior.mean` (default 0),
`prior.var` (default 100), and `B` for optional blocks. The model’s
`Data` list must contain `X` (the \\N \times p\\ design matrix), `y`
(binary 0/1 response or counts), and optionally `n` (number of trials
for binomial, defaulting to 1 for Bernoulli). The algorithm alternates
between drawing \\\omega_i \sim \text{PG}(n_i, x_i^\top\beta)\\ for all
observations (using the C++ Polya-Gamma random variate generator based
on the Devroye 2009 method) and drawing \\\beta\\ from the resulting
Gaussian conditional.

The advantage of PG is that it produces an exact Gibbs sampler for
logistic models with no tuning parameters and 100% acceptance rate,
leading to fast mixing and simple diagnostics. The disadvantage is that
it only applies to models with the specific logistic/binomial structure
and requires the design matrix to be available as `Data$X`. For
non-logistic models, use a general-purpose algorithm such as NUTS or
RWM. When using the Gibbs Sampler with `FC = "auto"` and the data
contains `pg_augment = TRUE`, the auto-FC mechanism automatically
detects the logistic structure and uses PG augmentation.

The C++ backend avoids the \\N \times N\\ diagonal matrix allocation by
computing \\X^\top \Omega X\\ as `X.t() * (X.each_col() % omega)` in
Armadillo, and calls the Polya-Gamma RNG directly via the internal C
function `rpg_vec()` rather than crossing the R/C++ boundary for each
draw.

#### Algorithm

1.  **Initialize** coefficients \\\beta_0\\; set prior \\\beta \sim
    \mathcal{N}(\mu_0, \Sigma_0)\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute linear predictor \\\psi_i = x_i^\top \beta\_{t-1}\\ for
        \\i = 1, \ldots, N\\
    2.  Draw \\\omega_i \sim \text{PG}(n_i, \psi_i)\\ independently for
        each observation
    3.  Form \\\Omega = \text{diag}(\omega_1, \ldots, \omega_N)\\ and
        \\\kappa_i = y_i - n_i/2\\
    4.  Compute \\V\_\omega = (X^\top \Omega X + \Sigma_0^{-1})^{-1}\\
        and \\m\_\omega = V\_\omega(X^\top\kappa + \Sigma_0^{-1}\mu_0)\\
    5.  Draw \\\beta_t \sim \mathcal{N}(m\_\omega, V\_\omega)\\
3.  **Return** \\\\\beta_1, \ldots, \beta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| prior.mean | rep(0, LIV) | Vague Gaussian prior centered at zero, matching the model used in this vignette. |
| prior.var | rep(100, LIV) | Variance of 100 corresponds to a vague prior on logistic regression coefficients (sd = 10). |
| B | NULL | Single block. |

Default `Specs` for PG. {.table}

#### Example with default specifications

``` r

## PG is restricted to logistic / binomial / negative-binomial regression
## with a logit link, and so cannot fit the Gaussian linear regression
## model used by every other example in this vignette. The snippet below
## defines a self-contained logistic regression and runs PG on it.
set.seed(42)
N_pg  <- 200L
X_pg  <- cbind(1, rnorm(N_pg), rnorm(N_pg))
beta_true <- c(-0.5, 1.2, -0.8)
y_pg  <- rbinom(N_pg, 1, 1 / (1 + exp(-as.vector(X_pg %*% beta_true))))

Data_pg <- list(N = N_pg, J = ncol(X_pg), X = X_pg, y = y_pg,
                mon.names = "LP", parm.names = paste0("beta", 1:3),
                PGF = function(Data) rnorm(Data$J))

Model_pg <- function(parm, Data) {
    eta <- as.vector(Data$X %*% parm)
    LL  <- sum(Data$y * eta - log1p(exp(eta)))
    LP  <- LL + sum(dnorm(parm, 0, 10, log = TRUE))
    list(LP = LP, Dev = -2 * LL, Monitor = LP,
         yhat = plogis(eta), parm = parm)
}

IV_pg <- GIV(Model_pg, Data_pg, PGF = TRUE)

fit <- lucifer(Model_pg, Data_pg, IV_pg,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "PG", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model_pg, Data = Data_pg)
plot(ppc, Style = "Fitted", Data = Data_pg)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = setNames(beta_true, Data_pg$parm.names))
```

#### See Also

- [Gibbs](#gibbs)
- [MWG](#mwg)

## Hamiltonian Monte Carlo exploits target geometry through gradient flow

Hamiltonian Monte Carlo augments the parameter space with auxiliary
momentum variables and simulates Hamiltonian dynamics using the gradient
of the log-posterior as a force field. The resulting proposals follow
the curvature of the target distribution, producing long-range moves
with high acceptance probability. Variants in this family automate
trajectory length selection (NUTS, HMCDA), add tempering for
multimodality (THMC), eliminate the momentum resampling step entirely
(MCHMC), or adapt the mass matrix without manual tuning (MEADS).

### Hamiltonian Monte Carlo

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 65% when L \> 1, or 57.4% when L = 1. The observed acceptance rate may be suitable in the interval \[60%,70%\] when L \> 1, or \[40%,80%\] when L = 1. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small number of parameters. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is difficult for a beginner. It has three algorithm specifications, all require tuning, and tuning is difficult. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. Each iteration involves numerous proposals, due to partial derivatives and L. |

Introduced under the name of hybrid Monte Carlo [\[28\]](#ref28), the
name Hamiltonian Monte Carlo (HMC) surpasses it in popularity in
statistics literature. HMC introduces auxiliary momentum variables with
independent, Gaussian proposals. Momentum variables receive alternate
updates, from simple updates to Metropolis updates. Metropolis updates
result in the proposal of a new state by computing a trajectory
according to Hamiltonian dynamics, from physics. Hamiltonian dynamics is
discretized with the leapfrog method. In this way, distant jumps can be
proposed and random-walk behavior avoided.

HMC has three algorithm specifications: a vector of the step size of the
leapfrog steps, `epsilon` (ε), that is equal in length to the number of
parameters, the number of leapfrog steps, `L`, and an optional mass
matrix `M`. When `L = 1`, HMC reduces to Langevin Monte Carlo (LMC),
also called the Metropolis-Adjusted Langevin Algorithm (MALA),
introduced by [\[30\]](#ref30). These tuning parameters must be adjusted
until the acceptance rate is appropriate. The optimal acceptance rate of
HMC is 65%, or in the case of LMC or MALA 57.4% is optimal. Tuning ε and
L, however, is very difficult. The trajectory length, εL must also be
considered.

Suggestions for tuning ε and L are found in [\[29\]](#ref29). When ε is
too large, the algorithm becomes unstable and suffers from a low
acceptance rate. When ε is too small, the algorithm takes too many small
steps and is inefficient. When L is too large, trajectory lengths (εL)
result in double-back behavior and become computationally
self-defeating. When L is too small, more random-walk behavior occurs
and mixing becomes slower. Even if ε and L are tuned optimally, a
well-tuned mass matrix may be necessary for efficient sampling.

If a user is new to tuning HMC algorithms, then good advice may be to
leave `L = 1` and begin with small values for ε, say 0.1 or smaller. It
is easy to experience problems when inexperienced, but HMC is a
rewarding algorithm once proficiency is acquired. As can be expected,
the adaptive extensions (AHMC, HMCDA, and NUTS), will also be easier,
since ε is adapted and does not require tuning (and in the case of NUTS,
\\L\\ does not require tuning).

Although HMC generates multivariate proposals, parameter correlation is
not taken into account unless a mass matrix is used. Some sources refer
to a mass matrix in HMC as a covariance matrix, and some as a precision
matrix. Here, the mass matrix is a covariance matrix. It is difficult to
tune for several reasons. The optimal mass matrix is different with
different configurations of the parameters, so what was optimal at one
point is sub-optimal at another. It is not recommended to substitute the
historical covariance matrix of samples while pursuing equilibrium.

Partial derivatives are required, and hence the parameters must be
differentiable everywhere the algorithm explores. Partial derivatives
are computationally intensive, and computational expense increases with
the number of parameters. For \\K\\ parameters and \\L\\ leapfrog steps,
there are \\L + KL\\ evaluations of the model specification function per
iteration. In practice, starting any of the algorithms in the HMC family
(AHMC, HMC, HMCDA, or THMC) in a region that is distant from density may
result in failure due to differentiation, unless manipulated with
priors.

Since HMC requires the approximation of partial derivatives, it is
slower per iteration than most algorithms, and much slower in higher
dimensions. Tuned well, HMC is an excellent algorithm, but tuning can be
very difficult. The AHMC algorithm and HMCDA are adaptive versions of
HMC in which ε is adapted based on recent history of acceptance and
rejection. The NUTS algorithm is a fully adaptive version that does not
require tuning of ε or \\L\\. Since HMC is not adaptive, it is suitable
as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    leapfrog steps \\L\\, mass matrix \\M\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Sample momentum \\r \sim \mathcal{N}(0, M)\\
    2.  Set \\(\theta^\*, r^\*) \leftarrow (\theta_t, r)\\; compute
        \\H_0 = U(\theta_t) + \tfrac{1}{2} r^\top M^{-1} r\\
    3.  **Leapfrog** for \\\ell = 1, \ldots, L\\:
        - \\r^\* \leftarrow r^\* + \tfrac{\epsilon}{2} \nabla \log
          p(\theta^\* \mid y)\\
        - \\\theta^\* \leftarrow \theta^\* + \epsilon \\ M^{-1} r^\*\\
        - \\r^\* \leftarrow r^\* + \tfrac{\epsilon}{2} \nabla \log
          p(\theta^\* \mid y)\\
    4.  Negate momentum: \\r^\* \leftarrow -r^\*\\
    5.  Compute \\\alpha = \min\\\bigl(1,\\ \exp(H_0 - H(\theta^\*,
        r^\*))\bigr)\\
    6.  With probability \\\alpha\\, accept \\\theta\_{t+1} =
        \theta^\*\\; otherwise \\\theta\_{t+1} = \theta_t\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | rep(1/LIV, LIV) (= rep(1/3, 3)) | Step-size scaled inversely with dimension; works well for unit-variance Gaussian targets. |
| L | 2 | Two leapfrog steps per iteration: a deliberately conservative default to prevent divergences before tuning. |
| m | diag(LIV) | Identity mass matrix as the no-information default. |

Default `Specs` for HMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "HMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AHMC](#ahmc)
- [HMCDA](#hmcda)
- [MALA](#mala)
- [NUTS](#nuts)
- [THMC](#thmc)

### Adaptive Hamiltonian Monte Carlo

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 65% when `L` \> 1, or 57.4% when `L` = 1. The observed acceptance rate may be suitable in the interval \[60%,70%\] when `L` \> 1, or \[40%,80%\] when `L` = 1. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small number of parameters. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is relatively difficult for a beginner. It has four algorithm specifications. Even though it is adaptive, it is difficult to tune. Since it is an adaptive algorithm, the user must regard diminishing adaptation. |
| Final Algorithm? | User Discretion. The HMC algorithm is recommended as the final algorithm, though AHMC may be used if diminishing adaptation occurs and adaptation ceases effectively. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. However, proposals are not generated with a multivariate distribution and a proposal covariance matrix. Each iteration involves numerous proposals, due to partial derivatives and `L`. |

This is an adaptive form of Hamiltonian Monte Carlo (HMC) called
Adaptive Hamiltonian Monte Carlo (AHMC). In AHMC, the ε parameter is
adapted, though `L` and `M` are not. When adapting, and considering K
parameters, AHMC multiplies εk by 0.8 when a proposal for parameter sk
has not been accepted in the last 10 iterations, or multiplies it by 1.2
when a proposal has been accepted at least 8 of the last 10 iterations,
as suggested by [\[29\]](#ref29).

AHMC has four algorithm specifications: `epsilon` or ε is the step-size,
`L` is the number of leapfrog steps, `M` is an optional mass matrix, and
`Periodicity` is the frequency in iterations for adaptation of ε.

Although AHMC generates multivariate proposals, parameter correlation is
not taken into account unless a mass matrix is used. Some sources refer
to a mass matrix in HMC as a covariance matrix, and some as a precision
matrix. Here, the mass matrix is a covariance matrix. It is difficult to
tune for several reasons. The optimal mass matrix is different with
different configurations of the parameters, so what was optimal at one
point is sub-optimal at another. It is not recommended to substitute the
historical covariance matrix of samples while pursuing equilibrium.

As with HMC, the AHMC algorithm is slower per iteration than many other
algorithms, but often produces chains with good mixing when well-tuned.
An alternative to AHMC that should perform better is HMCDA. AHMC is more
consistent with respect to time per iteration, because `L` remains
constant, than HMCDA and NUTS, which may have some iterations that are
much slower than others. If AHMC is used for adaptation, then the final,
non-adaptive algorithm should be HMC.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step sizes \\\epsilon_j\\
    for \\j = 1, \ldots, J\\, leapfrog count \\L\\, mass matrix \\M\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Sample momentum \\r \sim \mathcal{N}(0, M)\\
    2.  Set \\(\theta^\*, r^\*) \leftarrow (\theta_t, r)\\
    3.  **Leapfrog** for \\\ell = 1, \ldots, L\\: half-step \\r^\*\\,
        full-step \\\theta^\*\\ with step sizes \\\epsilon\\, half-step
        \\r^\*\\
    4.  Compute \\\alpha = \min\\\bigl(1,\\ \exp(-H(\theta^\*, r^\*) +
        H(\theta_t, r))\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta\_{t+1} =
        \theta^\*\\; otherwise \\\theta\_{t+1} = \theta_t\\
    6.  If \\t \bmod\\ `Periodicity` \\= 0\\, **adapt** \\\epsilon_j\\:
        multiply by \\0.8\\ if parameter \\j\\ not accepted in last 10
        iterations, or by \\1.2\\ if accepted \\\geq 8\\ of last 10
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | rep(1/LIV, LIV) | Initial step-size; will be adapted dynamically. |
| L | 2 | Initial number of leapfrog steps. |
| m | diag(LIV) | Identity mass matrix. |
| Periodicity | 1 | Adapt at every iteration. |

Default `Specs` for AHMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AHMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [HMCDA](#hmcda)
- [MALA](#mala)
- [NUTS](#nuts)
- [THMC](#thmc)

### No-U-Turn Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The user specifies the target acceptance rate as delta (60% is recommended). A recommended, suitable, observed acceptance rate may be in the interval \[delta-5%,delta+5%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small number of parameters. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is relatively easy for a beginner. It has only a few algorithm specifications. However, since it is adaptive, the user must regard diminishing adaptation. |
| Final Algorithm? | Yes, if not adaptive, otherwise: User Discretion. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. However, proposals are not generated with a multivariate distribution and a proposal covariance matrix. Each iteration involves numerous proposals, due to partial derivatives and L. |

The No-U-Turn Sampler (NUTS) is an extension of Hamiltonian Monte Carlo
(HMC) that adapts both the scalar step-size ε and scalar number of
leapfrog steps L. This is algorithm \#6 in Hoffman and Gelman (2012).

NUTS has four algorithm specifications: the number of adaptive
iterations `A`, the target acceptance rate `delta` or δ (and 0.6 is
recommended), a scalar step-size `epsilon` or ε, and the maximum number
of leapfrog steps `Lmax`.

Each iteration \\i ≤ A\\, NUTS adapts both ε and \\L\\, and coerces the
target acceptance rate δ. \\L\\ continues to change after adaptation
ends, but is not an adaptive parameter in the sense of destroying
ergodicity. The adaptive samples are discarded and only the thinned
non-adaptive samples are returned.

The main advantage of NUTS over other HMC algorithms is that NUTS is the
algorithm most likely to produce approximately independent samples, in
the sense of low autocorrelation. Due to computational complexity, NUTS
is slower per iteration than HMC, and the HMC family is among the
slowest. Despite this, NUTS often produces chains with excellent mixing,
and should outperform other adaptive versions of HMC, such as AHMC and
HMCDA. Per iteration, NUTS should generally perform better than other
HMC algorithms. Per minute, however, is another story.

NUTS has been extended elsewhere to allow for a non-diagonal mass matrix
(proposal covariance matrix for momentum). This extension are not yet
included here.

In complex and high-dimensional models, NUTS may produce approximately
independent samples much more slowly in minutes than other MCMC
algorithms, such as Adaptive Metropolis-within-Gibbs (AMWG). This is
because the combination of calculating partial derivatives and the
search each iteration for L is computationally intensive.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    target acceptance \\\delta\\, maximum tree depth \\L\_{\max}\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Sample momentum \\r \sim \mathcal{N}(0, I)\\; set \\\theta^- =
        \theta^+ = \theta_t\\, \\r^- = r^+ = r\\, tree depth \\j = 0\\,
        candidate set size \\n = 1\\
    2.  **Build tree** (recursive doubling):
        - Draw direction \\v \sim \\-1, +1\\\\ uniformly
        - Double the tree in direction \\v\\ by taking \\2^j\\ leapfrog
          steps with step size \\\epsilon\\
        - **U-turn check**: stop if \\(\theta^+ - \theta^-) \cdot r^- \<
          0\\ or \\(\theta^+ - \theta^-) \cdot r^+ \< 0\\
        - Select candidate from the new subtree with probability
          \\\propto \exp(-H(\theta^\*, r^\*))\\
        - Increment \\j\\; stop if \\j = L\_{\max}\\ or U-turn detected
    3.  If \\t \leq A\\, **adapt** \\\epsilon\\ via dual-averaging
        toward \\\delta\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | min(500, floor(Iterations/2)) (= 500) | Length of the dual-averaging warm-up window. 500 iterations are typically enough to converge the step size for low-dimensional targets. |
| delta | 0.6 | Target Metropolis acceptance probability; lower than Stan’s 0.8 default to favour longer trajectories on simple targets. |
| epsilon | NULL | Step size automatically initialized by find-reasonable-epsilon (Hoffman and Gelman 2014). |
| Lmax | 10 | Maximum tree depth of 2^10 = 1024 leapfrog steps; safety cap against runaway trajectories. |

Default `Specs` for NUTS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "NUTS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AHMC](#ahmc)
- [HMC](#hmc)
- [HMCDA](#hmcda)
- [MALA](#mala)

### Hamiltonian Monte Carlo with Dual-Averaging

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The user specifies the target acceptance rate as delta (65% is recommended). A recommended, suitable, observed acceptance rate may be in the interval \[delta-5%,delta+5%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small number of parameters. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is relatively difficult for a beginner. It has several algorithm specifications. Since it is an adaptive algorithm, the user must regard diminishing adaptation. |
| Final Algorithm? | Yes, when non-adaptive, otherwise: User Discretion. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. However, proposals are not generated with a multivariate distribution and a proposal covariance matrix. Each iteration involves numerous proposals, due to partial derivatives and L. |

The Hamiltonian Monte Carlo with Dual-Averaging (HMCDA) algorithm is an
extension of Hamiltonian Monte Carlo (HMC) that adapts the scalar
step-size parameter, ε, according to dual-averaging. This is algorithm
\#5 in Hoffman and Gelman (2012).

HMCDA has five algorithm specifications: the number of adaptive
iterations `A`, the target acceptance rate `delta` or δ (and 0.65 is
recommended), a scalar step-size `epsilon` or ε, and the trajectory
length `lambda` or λ. The trajectory length is scalar λ = εL, where
\\L\\ is the unspecified number of leapfrog steps that is determined
from ε and λ. Each iteration i ≤ A, HMCDA adapts and coerces the target
acceptance rate δ.

As with HMC, the HMCDA algorithm is slower per iteration than many other
algorithms, but often produces chains with good mixing. HMCDA should
outperform Adaptive Hamiltonian Monte Carlo (AHMC), and iterates faster
as well, unless L becomes large. When mixing is inadequate, consider
switching to the MALA or NUTS algorithm. When parameters are highly
correlated, another algorithm should be preferred in which correlation
is taken into account, such as AMM, or in which the algorithm is
generally invariant to correlation, such as twalk. When adaptive, it is
not suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    trajectory length \\\lambda = \epsilon L\\, target acceptance
    \\\delta\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute \\L = \max(1, \lfloor \lambda / \epsilon \rfloor)\\
        leapfrog steps
    2.  Sample momentum \\r \sim \mathcal{N}(0, I)\\; set \\(\theta^\*,
        r^\*) \leftarrow (\theta_t, r)\\
    3.  **Leapfrog** for \\\ell = 1, \ldots, L\\: half-step \\r^\*\\,
        full-step \\\theta^\*\\, half-step \\r^\*\\
    4.  Compute \\\alpha = \min\\\bigl(1,\\ \exp(H_0 - H(\theta^\*,
        r^\*))\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta\_{t+1} =
        \theta^\*\\
    6.  If \\t \leq A\\, **adapt** \\\epsilon\\ via dual-averaging:
        update \\\bar{H}\_t\\ toward \\\delta - \alpha\\, set \\\log
        \epsilon_t = \mu - \sqrt{t}/\gamma \cdot \bar{H}\_t\\
3.  After warmup, fix \\\epsilon \leftarrow \bar{\epsilon}\\ (averaged
    step size)
4.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | min(500, floor(Iterations/2)) (= 500) | Dual-averaging warm-up window. |
| delta | 0.65 | Target Metropolis acceptance rate recommended by Hoffman and Gelman (2014) for HMCDA. |
| epsilon | NULL | Step size auto-initialized. |
| Lmax | 10 | Upper bound on the trajectory length budget. |
| lambda | 1 | Total trajectory length; combined with epsilon to derive L. |

Default `Specs` for HMCDA. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "HMCDA", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AHMC](#ahmc)
- [HMC](#hmc)
- [MALA](#mala)
- [NUTS](#nuts)

### Tempered Hamiltonian Monte Carlo

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 65% when L \> 1, or 57.4% when L = 1. The observed acceptance rate may be suitable in the interval \[60%,70%\] when L \> 1, or \[40%,80%\] when L = 1. |
| Applications | This is a widely applicable, general-purpose algorithm that is best suited to models with a small number of parameters. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is difficult for a beginner. It has a several algorithm specifications, and these are difficult to tune. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. Each iteration involves numerous proposals, due to partial derivatives and L. |

The Tempered Hamiltonian Monte Carlo (THMC) algorithm is an extension of
Hamiltonian Monte Carlo (HMC) to include another algorithm
specification: Temperature, which must be positive. When greater than 1,
the algorithm should explore more diffuse distributions, and may be
helpful with multimodal distributions.

THMC has four algorithm specifications: step-size `epsilon`, the number
`L` of leapfrog steps, an optional mass matrix `M`, and `Temperature`.
Algorithm specifications are the same as for HMC, with the exception of
temperature, which is described below.

There are a variety of ways to include tempering in HMC, and this
algorithm, named here as THMC, uses “tempered trajectory”, as described
by [\[29\]](#ref29). When L \> 1 and during the first half of the
leapfrog steps, the momentum is increased (heated) by multiplying it by
√T, where T is temperature, each leapfrog step. In the last half of the
leapfrog steps, the momentum decreases (is cooled down) by dividing it
by √T. The momentum is largest in the middle of the leapfrog steps,
where mode-switching behavior becomes most likely to occur. This
preserves the trajectory, εL.

As with HMC, THMC is a difficult algorithm to tune. Since THMC is
non-adaptive, it is sufficient as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    leapfrog steps \\L\\, mass matrix \\M\\, temperature \\T\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Sample momentum \\r \sim \mathcal{N}(0, M)\\; set \\(\theta^\*,
        r^\*) \leftarrow (\theta_t, r)\\
    2.  **Tempered leapfrog** for \\\ell = 1, \ldots, L\\:
        - If \\\ell \leq L/2\\: **heat** momentum \\r^\* \leftarrow r^\*
          \cdot \sqrt{T}\\ (increasing kinetic energy)
        - Standard leapfrog step: half-step \\r^\*\\, full-step
          \\\theta^\*\\, half-step \\r^\*\\
        - If \\\ell \> L/2\\: **cool** momentum \\r^\* \leftarrow r^\* /
          \sqrt{T}\\ (decreasing kinetic energy)
    3.  Negate momentum: \\r^\* \leftarrow -r^\*\\
    4.  Compute \\\alpha = \min\\\bigl(1,\\ \exp(H_0 - H(\theta^\*,
        r^\*))\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta\_{t+1} =
        \theta^\*\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | rep(1/LIV, LIV) | Step size scaled inversely with dimension. |
| L | 2 | Two leapfrog steps. |
| m | diag(LIV) | Identity mass matrix. |
| Temperature | 2 | Doubles momentum variance during the first half of each trajectory to encourage barrier crossing. |

Default `Specs` for THMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "THMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [MALA](#mala)

### Microcanonical Hamiltonian Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is effectively 100%, since MCHMC does not use a Metropolis accept/reject step. Invalid trajectories (where the potential energy exceeds total energy) are rejected by momentum reflection. |
| Applications | This is a widely applicable, general-purpose algorithm. It is particularly effective in high-dimensional problems where standard HMC suffers from step size sensitivity. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has four algorithm specifications. The step size and trajectory length require tuning, though step size adaptation is performed automatically during warmup. |
| Final Algorithm? | User Discretion. MCHMC may be used as a final algorithm after adaptation ceases, but the user should verify that the adapted step size has stabilized. |
| Proposal | Multivariate. |

Microcanonical Hamiltonian Monte Carlo (MCHMC) was proposed by
[\[32\]](#ref32) as an alternative to canonical HMC that operates on the
microcanonical manifold where the total Hamiltonian \\H = U(\theta) +
K(r)\\ is held constant, rather than resampling momentum from a fixed
distribution at each iteration. In canonical HMC, the kinetic energy is
refreshed by drawing new momenta from a Gaussian distribution,
discarding all momentum information from the previous trajectory. MCHMC
instead performs a partial momentum refresh via a rotation:
\\r\_{\text{new}} = \cos(\alpha) \\ r\_{\text{old}} + \sin(\alpha) \\
\eta\\, where \\\eta \sim \mathcal{N}(0, I)\\ and \\\alpha\\ controls
the refresh rate. The momentum is then rescaled to satisfy the energy
constraint \\K(r) = H - U(\theta)\\, ensuring that the sampler remains
on the microcanonical manifold.

MCHMC has four algorithm specifications: `epsilon` is the leapfrog step
size, `L` is the number of leapfrog steps per iteration, `alpha` is the
partial momentum refresh angle in radians (smaller values preserve more
momentum and improve autocorrelation, while larger values inject more
randomness), and `A` is the number of warmup iterations during which the
step size is adapted. During warmup, the step size is adjusted to
maintain a target acceptance rate of 90%; if the acceptance rate falls
below 72%, the step size is reduced by 5%, and if it exceeds 90%, the
step size is increased by 5%.

The leapfrog integration follows the standard Stormer-Verlet scheme.
After \\L\\ leapfrog steps, the algorithm checks energy conservation: if
\\\|H\_{\text{new}} - H\_{\text{old}}\| \> 10\\ nats, the trajectory is
flagged as divergent (reported in diagnostics). When the potential
energy at the proposed position exceeds the total energy \\E\\, the
kinetic energy would need to be negative, which is unphysical; in this
case the momentum is reflected and the proposal is rejected, returning
to the previous state.

The advantage of MCHMC over canonical HMC is that the persistent
momentum produces longer effective trajectories with less wasted kinetic
energy, leading to better mixing per gradient evaluation in many
problems. The microcanonical constraint also provides a natural energy
diagnostic: if the sampler cannot maintain the energy constraint, this
signals numerical integration issues. The disadvantage is that MCHMC is
more sensitive to the initial energy level and requires the target
distribution to be sufficiently smooth for the microcanonical manifold
to be well-defined.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, momentum \\r_0 \sim
    \mathcal{N}(0, I)\\; compute total energy \\E = U(\theta_0) +
    K(r_0)\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Partial momentum refresh**: draw \\\eta \sim \mathcal{N}(0,
        I)\\; set \\r \leftarrow \cos(\alpha)\\ r\_{t-1} +
        \sin(\alpha)\\ \eta\\
    2.  **Rescale** momentum to satisfy energy constraint: \\r
        \leftarrow r \cdot \sqrt{(E - U(\theta\_{t-1})) / K(r)}\\ (if
        \\E - U(\theta\_{t-1}) \> 0\\)
    3.  **Leapfrog** for \\\ell = 1, \ldots, L\\ with step size
        \\\epsilon\\:
        - \\r \leftarrow r + \tfrac{\epsilon}{2} \nabla \log p(\theta
          \mid y)\\; \\\theta \leftarrow \theta + \epsilon\\ r\\; \\r
          \leftarrow r + \tfrac{\epsilon}{2} \nabla \log p(\theta \mid
          y)\\
    4.  If \\U(\theta^\*) \> E\\ (kinetic energy would be negative),
        **reflect** momentum and retain \\\theta\_{t-1}\\
    5.  If \\t \leq A\\, **adapt** \\\epsilon\\: reduce by 5% if
        acceptance \\\< 72\\\\, increase by 5% if \\\> 90\\\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 | Conservative step size for the constant-energy trajectory. |
| L | 10 | Ten partial-momentum-refreshment steps per iteration. |
| alpha | 0.1 | Partial momentum refreshment rate; 0.1 keeps trajectories long while breaking periodicity. |
| A | min(500, floor(Iterations/2)) (= 500) | Step-size adaptation window. |

Default `Specs` for MCHMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MCHMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [HMCDA](#hmcda)
- [NUTS](#nuts)
- [MEADS](#meads)
- [SGHMC](#sghmc)

### MEADS

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The user specifies the target acceptance rate as delta (65% is recommended). The observed acceptance rate may be suitable in the interval \[50%,80%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. MEADS is specifically designed to be tuning-free: the step size and mass matrix are adapted automatically from dual chains. |
| Difficulty | This algorithm is easy for a beginner despite having four algorithm specifications, because the algorithm is designed to be fully automatic. The user primarily needs to set the number of warmup iterations. |
| Final Algorithm? | User Discretion. MEADS may be used as a final algorithm after the warmup period, when dual-averaging has converged and the mass matrix has stabilized. |
| Proposal | Multivariate. |

MEADS (Metropolis-adjusted Evolving-Adaptation Dual-averaging Sampler)
was proposed by [\[33\]](#ref33) as a tuning-free generalization of
Hamiltonian Monte Carlo. The central idea is to run two parallel chains
simultaneously and use the inter-chain statistics to estimate a diagonal
mass matrix, while the step size is adapted via the dual-averaging
scheme of [\[34\]](#ref34). This eliminates the need for a separate
warmup phase where the mass matrix is estimated from a single chain’s
history, which can be unreliable in the early stages of sampling.

MEADS has four algorithm specifications: `delta` is the target
acceptance rate (0.65 is recommended), `epsilon` is the initial step
size (this is adapted during warmup via dual-averaging), `L` is the
number of leapfrog steps per iteration, and `A` is the number of warmup
iterations for adaptation. The dual-averaging mechanism adjusts the step
size to achieve the target acceptance rate \\\delta\\, using the average
acceptance probability across both chains. The adapted step size
converges to a stable value \\\bar{\epsilon}\\ that is used for all
post-warmup iterations.

The mass matrix estimation proceeds as follows. Each iteration, the
difference \\d_t = \theta^{(1)}\_t - \theta^{(2)}\_t\\ between the two
chains’ positions is computed, and Welford’s online algorithm maintains
running estimates of the mean and variance of these differences. After
an initial burn-in of 10 iterations, the diagonal of the mass matrix is
set to the estimated variance of inter-chain differences, clipped at a
minimum of \\10^{-8}\\. This approach is more robust than estimating the
mass matrix from a single chain’s sample covariance, because the
inter-chain differences are approximately Gaussian even when the
marginal distributions are not, and the estimate stabilizes faster.

Both chains undergo independent HMC transitions using the same step size
and mass matrix, each with a standard Metropolis accept/reject step. The
acceptance probabilities from both chains are averaged for the
dual-averaging update. Divergent transitions (where \\\Delta H \> 10\\
nats) are tracked and reported. Only the first chain’s samples are
retained; the second chain exists solely to provide information for the
mass matrix estimate.

The advantage of MEADS is that it requires essentially no manual tuning:
the step size adapts automatically via dual-averaging, and the mass
matrix adapts from inter-chain statistics. This makes it an attractive
default for users who do not wish to experiment with tuning parameters.
The disadvantage is the computational overhead of running two chains
(roughly doubling the cost per iteration compared to single-chain HMC),
and the restriction to a diagonal mass matrix, which does not capture
posterior correlations. For highly correlated posteriors, NUTS with a
dense mass matrix may be preferable.

#### Algorithm

1.  **Initialize** two chains \\\theta^{(1)}\_0, \theta^{(2)}\_0\\; step
    size \\\epsilon\\; diagonal mass matrix \\M = I\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each chain \\c \in \\1, 2\\\\:
        - Sample \\r^{(c)} \sim \mathcal{N}(0, M)\\; leapfrog \\L\\
          steps with step size \\\epsilon\\ and mass \\M\\
        - Accept or reject via Metropolis step; record acceptance
          probability \\\alpha^{(c)}\\
    2.  **Dual-averaging** (if \\t \leq A\\): update \\\epsilon\\ toward
        target \\\delta\\ using average \\\bar{\alpha} = (\alpha^{(1)} +
        \alpha^{(2)})/2\\
    3.  **Mass matrix update** (if \\t \> 10\\): compute inter-chain
        difference \\d_t = \theta^{(1)}\_t - \theta^{(2)}\_t\\; update
        running variance of \\\\d_1, \ldots, d_t\\\\ via Welford’s
        algorithm; set \\M\_{\text{diag}} \leftarrow
        \max\\\bigl(\text{Var}(d), 10^{-8}\bigr)\\
3.  After warmup, fix \\\epsilon \leftarrow \bar{\epsilon}\\ and \\M\\
4.  **Return** samples from the first chain

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| delta | 0.6 | Target acceptance rate for the dual-averaging stage. |
| epsilon | NULL | Auto-initialized. |
| L | 10 | Ten ensemble HMC integration steps. |
| A | min(500, floor(Iterations/2)) (= 500) | Adaptation window. |

Default `Specs` for MEADS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MEADS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [HMCDA](#hmcda)
- [NUTS](#nuts)
- [MCHMC](#mchmc)

### autoMALA

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is approximately 57.4%, which is the optimal rate for MALA in the high-dimensional limit. Step size adaptation targets this rate automatically. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is designed as a fully automated alternative to MALA that requires no user tuning, adapting the step size locally at each iteration. |
| Difficulty | This algorithm is easy for a beginner. It has three algorithm specifications, and the key feature is that no warmup adaptation is needed because the step size is tuned on the fly. |
| Final Algorithm? | Yes. autoMALA may be used directly as a final algorithm without warmup. |
| Proposal | Multivariate. |

autoMALA was introduced by [\[35\]](#ref35) as a gradient-based MCMC
algorithm that automatically tunes its step size at every iteration
without requiring a separate warmup phase. Standard MALA requires
careful step size tuning (typically via dual-averaging during warmup);
if the step size is wrong, acceptance rates collapse. autoMALA instead
performs a local search at each iteration: starting from an initial
candidate step size \\\epsilon\_{\max}\\, it repeatedly halves the step
size until the Metropolis acceptance probability exceeds a threshold
derived from the target acceptance rate. This local adaptation means the
step size varies across the parameter space, automatically using larger
steps in flat regions and smaller steps in regions of high curvature,
without any global adaptation scheme.

autoMALA has three algorithm specifications: `epsilon_max` is the
initial maximum step size (default auto-scaled from the gradient
magnitude), `target_accept` is the target acceptance probability
(default 0.574), and `max_search` is the maximum number of step size
halvings per iteration (default 10).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, maximum step size
    \\\epsilon\_{\max}\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Set \\\epsilon = \epsilon\_{\max}\\
    2.  **Repeat** up to `max_search` times:
        - Propose \\\theta^\* = \theta\_{t-1} +
          \frac{\epsilon^2}{2}\nabla \log \pi(\theta\_{t-1}) + \epsilon
          z\\, where \\z \sim \mathcal{N}(0, I_d)\\
        - Compute Metropolis-Hastings acceptance probability \\\alpha\\
        - If \\\alpha \geq\\ target: accept this step size; **break**
        - Otherwise: \\\epsilon \leftarrow \epsilon / 2\\
    3.  Accept \\\theta^\*\\ with probability \\\min(1, \alpha)\\;
        otherwise retain \\\theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon_max | 1.0 | Upper bound on the proposal scale; the doubling/halving search starts here and contracts as needed. |
| target_accept | 0.574 | Roberts/Rosenthal optimum for MALA-style proposals. |
| max_search | 10 | At most ten doubling/halving steps per proposal: bounds worst-case cost while still allowing 1024-fold scale changes. |

Default `Specs` for autoMALA. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "autoMALA", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [MALT](#malt)
- [Barker Proposal](#barker)

### MALT

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on the step size and friction parameter. Typical target rates are 60-80%, interpolating between MALA-like and HMC-like behavior. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It unifies MALA and HMC in a single framework by introducing a friction parameter that controls the dissipation of kinetic energy. |
| Difficulty | This algorithm is moderately easy for a beginner. It has three algorithm specifications. The friction parameter gamma interpolates smoothly between MALA (gamma=1) and HMC (gamma=0). |
| Final Algorithm? | Yes. MALT may be used as a final algorithm. |
| Proposal | Multivariate. |

Metropolis Adjusted Langevin Trajectories (MALT) was introduced by
[\[36\]](#ref36) as a framework that unifies MALA and HMC by introducing
a friction (damping) parameter \\\gamma \in \[0, 1\]\\ into the leapfrog
integrator. When \\\gamma = 0\\, the dynamics are purely Hamiltonian and
MALT reduces to standard HMC; when \\\gamma = 1\\ and \\L = 1\\, the
dynamics are fully dissipative and MALT reduces to MALA. Intermediate
values of \\\gamma\\ produce a partially damped trajectory that combines
the local gradient information of MALA with the long-range ballistic
exploration of HMC. At each leapfrog step, the momentum is partially
refreshed: \\p \leftarrow \sqrt{1 - \gamma}\\p + \sqrt{\gamma}\\\xi\\
where \\\xi \sim \mathcal{N}(0, I_d)\\. This partial refreshment
controls how much kinetic energy survives between steps, trading off
between the ergodic properties of Langevin dynamics and the
trajectory-following properties of Hamiltonian dynamics.

MALT has three algorithm specifications: `epsilon` is the leapfrog step
size, `L` is the number of leapfrog steps, and `gamma` is the friction
parameter (default 0.5).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    trajectory length \\L\\, friction \\\gamma\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum \\p \sim \mathcal{N}(0, I_d)\\
    2.  **For** \\\ell = 1, \ldots, L\\:
        - Partial momentum refreshment: \\p \leftarrow \sqrt{1 -
          \gamma}\\p + \sqrt{\gamma}\\\xi\\, \\\xi \sim \mathcal{N}(0,
          I_d)\\
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\nabla \log \pi(\theta)\\
        - Full-step position: \\\theta \leftarrow \theta + \epsilon p\\
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\nabla \log \pi(\theta)\\
    3.  Compute acceptance probability \\\alpha\\ from the augmented
        Hamiltonian
    4.  Accept with probability \\\min(1, \alpha)\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step-size scaling for diffusive samplers. |
| L | 10 | Ten leapfrog steps per trajectory. |
| gamma | 1.0 | Friction coefficient; gamma = 1 corresponds to underdamped Langevin trajectories. |

Default `Specs` for MALT. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MALT", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [HMC](#hmc)
- [autoMALA](#automala)

### AAPS

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is adapted during warmup via the delta parameter. Typical targets are 60-80%. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It automatically selects trajectory lengths by detecting apogee points (turning points in the Hamiltonian trajectory), providing a principled alternative to NUTS for trajectory termination. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has five algorithm specifications. The apogee-based termination criterion is conceptually simpler than the U-turn criterion of NUTS but requires more leapfrog steps per iteration. |
| Final Algorithm? | User Discretion. AAPS may be used as a final algorithm after warmup adaptation stabilizes. |
| Proposal | Multivariate (Hamiltonian dynamics). |

The Apogee-to-Apogee Path Sampler (AAPS) was introduced by
[\[37\]](#ref37) as an alternative to NUTS that determines trajectory
length by detecting apogee points, defined as local maxima of the
potential energy \\U(\theta) = -\log \pi(\theta)\\ along the Hamiltonian
trajectory. At an apogee, the particle has converted most of its kinetic
energy into potential energy and is about to turn back; this is a
natural point to stop. AAPS runs the leapfrog integrator forward until
\\n\_{\text{apogees}}\\ apogee points are encountered, then selects the
proposal from the states between the first and last apogee using a
multinomial weighting scheme that preserves detailed balance. The apogee
criterion \\U(\theta\_\ell) \> U(\theta\_{\ell-1})\\ followed by
\\U(\theta\_\ell) \> U(\theta\_{\ell+1})\\ is simpler to evaluate than
the U-turn criterion of NUTS and provides a different (sometimes
superior) notion of “the trajectory has gone far enough.”

AAPS has five algorithm specifications: `A` is the number of warmup
iterations for step size adaptation (default 500), `delta` is the target
acceptance rate (default 0.70), `epsilon` is the initial step size
(auto-scaled if NULL), `n_apogees` is the number of apogee points to
detect before stopping (default 2), and `Lmax` is the maximum number of
leapfrog steps as a safety limit (default 1000).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum \\p \sim \mathcal{N}(0, I_d)\\
    2.  Run leapfrog forward from \\(\theta\_{t-1}, p)\\, recording
        trajectory states
    3.  After each leapfrog step, check if an apogee occurred:
        \\U(\theta\_\ell) \> U(\theta\_{\ell-1})\\ and
        \\U(\theta\_{\ell+1}) \< U(\theta\_\ell)\\
    4.  Stop when \\n\_{\text{apogees}}\\ apogees are detected or
        \\L\_{\max}\\ steps reached
    5.  Select proposal from trajectory states between first and last
        apogee with multinomial weights \\\propto \exp(-H(\theta\_\ell,
        p\_\ell))\\
    6.  Accept with Metropolis correction
    7.  If \\t \leq A\\: adapt \\\epsilon\\ via dual-averaging toward
        target \\\delta\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | min(500, floor(Iterations/2)) (= 500) | Warm-up window for step-size adaptation. |
| delta | 0.65 | Target acceptance rate. |
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Initial step size for the apogee-to-apogee trajectory. |
| n_apogees | 2 | Two apogees per iteration: enough to span the modal region without excessive computation. |
| Lmax | 1000 | Upper bound on the number of leapfrog steps per trajectory. |

Default `Specs` for AAPS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AAPS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [NUTS](#nuts)
- [HMC](#hmc)
- [GIST](#gist)

### GIST

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on the Gibbs-selected step size. Typical acceptance rates are 60-80%. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It augments NUTS with Gibbs-based step size selection, eliminating the need for warmup adaptation entirely. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has seven algorithm specifications. The Gibbs step size selection adds complexity but removes the warmup requirement. |
| Final Algorithm? | Yes. GIST is designed to be used without warmup, making it immediately usable as a final algorithm. |
| Proposal | Multivariate (Hamiltonian dynamics). |

GIST (Gibbs Self-Tuning NUTS) was introduced by [\[38\]](#ref38) as an
extension of the No-U-Turn Sampler that eliminates the need for
warmup-based step size adaptation by treating the step size as an
auxiliary variable sampled via Gibbs. At each iteration, instead of
using a single fixed step size, GIST evaluates the NUTS trajectory at
\\K\\ candidate step sizes drawn from a discrete grid centered (in log
space) on the current step size, and selects among them with
probabilities proportional to \\\exp(-\alpha\_{\text{gist}} \cdot
H\_{\text{error}})\\, where \\H\_{\text{error}}\\ is the Hamiltonian
error at each candidate. This Gibbs selection automatically concentrates
on step sizes that produce low integration error, adapting locally
without any global warmup scheme. The sharpness parameter
\\\alpha\_{\text{gist}}\\ controls how strongly the selection favors
low-error candidates.

GIST has seven algorithm specifications: `A` is the number of warmup
iterations (default 0, since GIST does not require warmup), `delta` is
the target acceptance rate for optional dual-averaging warmup (default
0.80), `epsilon` is the initial step size (auto-scaled if NULL), `Lmax`
is the maximum tree depth (default 10), `K` is the number of candidate
step sizes (default 5), `alpha_gist` is the Gibbs selection sharpness
(default 1.0), and `grid_range` is the half-width of the log2 step size
grid (default 2.0).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum \\p \sim \mathcal{N}(0, I_d)\\
    2.  Construct candidate grid: \\\epsilon_k = \epsilon_0 \cdot
        2^{u_k}\\ for \\u_k\\ uniformly spaced in \\\[-r, r\]\\, \\k =
        1, \ldots, K\\
    3.  **For** each candidate \\\epsilon_k\\: run NUTS trajectory
        (build tree until U-turn) and record Hamiltonian error \\\Delta
        H_k\\
    4.  Select \\\epsilon^\*\\ with probability \\\propto
        \exp(-\alpha\_{\text{gist}} \cdot \|\Delta H_k\|)\\
    5.  Use the NUTS proposal from the selected step size trajectory
    6.  Accept with the standard NUTS multinomial selection criterion
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | min(500, floor(Iterations/2)) (= 500) | Step-size warm-up window. |
| delta | 0.6 | Target acceptance rate (matches NUTS). |
| epsilon | NULL | Reference step size auto-initialized. |
| Lmax | 10 | Maximum tree depth of 2^10 leapfrog steps. |
| K | 5 | Five candidate step sizes on the log-uniform grid: a sweet spot between exploration and per-iteration cost. |
| alpha_gist | 1.0 | Sharpness of the softmax used to sample step sizes; 1.0 gives a moderate selection bias. |
| grid_range | 2.0 | Half-width (log2) of the candidate grid; spans factor-4 around the reference epsilon. |

Default `Specs` for GIST. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "GIST", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [NUTS](#nuts)
- [AAPS](#aaps)
- [autoMALA](#automala)

## Langevin and stochastic gradient methods scale to large datasets

Langevin dynamics adds gradient-driven drift to random-walk proposals,
yielding faster convergence than pure random walk at the cost of
requiring first-order derivatives. Stochastic gradient variants (SGLD,
SGHMC) replace the full-data gradient with a minibatch estimate, making
each iteration computationally cheap for large datasets. The
Metropolis-adjusted variant (MALA) corrects for discretization error
through an accept-reject step, while the stochastic gradient variants
typically forgo correction and rely on decreasing step sizes for
asymptotic exactness.

### Metropolis-Adjusted Langevin Algorithm

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 57.4%. The observed acceptance rate may be suitable in the interval \[40%,80%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. The number of model evaluations per iteration increases with the number of parameters. |
| Difficulty | This algorithm is easy for a beginner to use. Although it has five algorithm specifications, it is fully automatic at the default values. |
| Final Algorithm? | Yes, when not adaptive. |
| Proposal | Multivariate. |

Also called Langevin Monte Carlo (LMC), the Metropolis-Adjusted Langevin
Algorithm (MALA) was proposed in [\[40\]](#ref40), an adaptive version
in [\[39\]](#ref39), and an alternative adaptive version in
[\[41\]](#ref41). MALA was inspired by stochastic models of molecular
dynamics. MALA is an extension of the multivariate random-walk
Metropolis (RWM) algorithm that includes partial derivatives to improve
mixing. [\[40\]](#ref40) presented ULA, MALA, and MALTA, and recommended
MALTA. MALTA is a refinement of MALA that uses a truncated drift, where
the drift parameter is a step-size parameter for the partial
derivatives. The non-adaptive version of MALA is nearly equivalent to
the Hamiltonian Monte Carlo (HMC) algorithm with L=1 leapfrog steps,
except MALA also includes a proposal covariance matrix.

The original, non-adaptive MALA family is difficult to tune, and optimal
tuning differs between transience and stationarity. For this reason, the
MALA presented here is the adaptive version presented in
[\[39\]](#ref39). This adaptive MALA uses stochastic approximation to
update an expectation, scale, and covariance every iteration. The scale
of the proposal is adapted to obtain a target acceptance rate.

This version of MALA has five algorithm specifications: `A` defaults to
`1e7` as the maximum acceptable value of the Euclidean norm of the
adaptive parameters `mu` and `sigma`, and the Frobenius norm of the
covariance matrix, `alpha.star` is the target acceptance rate and
defaults to 57.4%, `gamma` defaults to 1 and accepts a constant in
\[0,T\] where \\T\\ is iterations and controls decay in adaptation,
`delta` defaults to 1 and is a constant in the bounded drift function,
and `epsilon` is a vector of length two that defaults to `c(1e-6,1e-7)`,
in which the first element is the acceptable minimum of adaptive scale
`sigma`, and the second element is added to the diagonal of the
covariance matrix for regularization. The expectation, scale, and
covariance adapt each iteration, and the amount of adaptation is a
decreasing function `gamma/t` for \\T\\ iterations. When `gamma=0`, the
algorithm does not adapt. Otherwise, `gamma` is the iteration through
which full adaptation occurs and after which adaptation decays. The
optimal acceptance rate is 57.4%, and is acceptable in the interval
\[40%, 80%\].

MALA approximates partial derivatives for multivariate proposals.
Approximating partial derivatives is computationally expensive, and
requires \\J + 1\\ model evaluations per \\J\\ parameters per iteration.
This results in a multivariate algorithm that is slightly slower per
iteration than a traditional componentwise algorithm, though the higher
acceptance rate and additional information from partial derivatives may
allow it to mix better and approach convergence faster. Unlike most
componentwise algorithms, this version of MALA accounts for parameter
correlation.

When non-adaptive, MALA is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, scale \\\sigma\\, drift
    parameter \\\delta\\, covariance \\\Sigma = I\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute gradient \\g = \nabla \log p(\theta\_{t-1} \mid y)\\;
        truncate if \\\\g\\ \> \delta\\: \\g \leftarrow \delta \cdot g /
        \\g\\\\
    2.  Propose \\\theta^\* \sim \mathcal{N}\\\bigl(\theta\_{t-1} +
        \tfrac{\sigma^2}{2} \Sigma \\ g,\\ \sigma^2 \Sigma\bigr)\\
    3.  Compute \\\alpha = \min\\\bigl(1,\\ \frac{p(\theta^\* \mid y)
        \cdot q(\theta\_{t-1} \mid \theta^\*)}{p(\theta\_{t-1} \mid y)
        \cdot q(\theta^\* \mid \theta\_{t-1})}\bigr)\\
        (Metropolis-Hastings ratio with asymmetric proposal)
    4.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\
    5.  If adaptive (\\\gamma \> 0\\), update \\\sigma\\, expectation
        \\\mu\\, and \\\Sigma\\ via stochastic approximation with decay
        \\\gamma / t\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | 1e7 | Maximum norm allowed for the bounded drift; effectively disables clipping for well-behaved targets. |
| alpha.star | 0.574 | Roberts/Rosenthal optimum acceptance rate for MALA. |
| gamma | 0.66 | Decay exponent for the Robbins-Monro adaptation of the step size. |
| delta | 1 | Constant in the bounded drift function. |
| epsilon | c(1e-6, 1e-7) | Lower and upper bounds on the adaptive step size. |

Default `Specs` for MALA. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MALA", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AHMC](#ahmc)
- [HMC](#hmc)
- [HMCDA](#hmcda)
- [NUTS](#nuts)
- [THMC](#thmc)

### Stochastic Gradient Langevin Dynamics

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 100%. |
| Applications | This is a widely applicable, general-purpose algorithm that is specifically designed for big data. |
| Difficulty | This algorithm is easy for a beginner, though the step size must be specified and tuned. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

The Stochastic Gradient Langevin Dynamics (SGLD) algorithm of
[\[42\]](#ref42) is the first MCMC algorithm designed specifically for
big data. Traditional MCMC algorithms require the entire data set to be
included in the model evaluation each iteration. In contrast, SGLD reads
and processes only a small, randomly selected batch of records each
iteration. In addition to saving computation time, the entire data set
does not need to be loaded into memory at once.

SGLD expands the stochastic gradient descent optimization algorithm to
include Gaussian noise with Langevin dynamics. The multivariate proposal
is merely the vector of current values plus a step size times the
gradient, plus Gaussian noise. The scale of the Gaussian noise is
determined by the step size, ε.

SGLD has five algorithm specifications: `epsilon` or ε is the step size,
`file` accepts the quoted name of a .csv file that is the big data set,
`Nr` is the number of rows in the big data set, `Nc` is the number of
columns in the big data set, and `size` is the number of rows to be read
and processed each iteration.

Since SGLD is designed for big data, the entire data set is not included
in the Data list, but one small batch must be included and named \\X\\.
All data must be included. For example, both the dependent variable
\\y\\ and design matrix \\X\\ in linear regression are included. The
requirement for the small batch to be in `Data` is so that numerous
checks may be passed after
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
is called and before the SGLD algorithm begins. Each iteration, SGLD
uses the scan function, without headers, to read a random block of rows
from, say, `X.csv`, stores it in `Data$X`, and passes it to the Model
specification function. The Model function must differ from the other
examples found in the lucifer package in that multiple objects, such as
\\X\\ and \\y\\ must be read from `Data$X`, where usually there is both
`Data$X` and `Data$y`.

The user tunes SGLD with step size ε via the `epsilon` argument, which
accepts either a scalar for a constant step size or a vector equal in
length to the number of iterations for an annealing schedule. The step
size must remain in the interval (0,1). The user may define an annealing
schedule with any function desired. Examples are given in
[\[42\]](#ref42) of decreasing schedules, as well as for using a
constant. When `ε = 0`, both the stochastic gradient and Langevin
dynamics components of the equation are reduced to zero and the
algorithm will not move. When ε is too large, degenerate results occur.
A good recommendation seems to be to begin with ε set to `1/Nr`. From
testing on numerous examples, it seems preferable to perform several
short runs and experiment with a constant, rather than using a
decreasing schedule, but your mileage may vary.

The SGLD algorithm presented here is the simplest one, which is equation
4 in [\[42\]](#ref42). Other components were also proposed such as a
preconditioning matrix and covariance matrix, though these are not
currently included here.

SGLD does not include a Metropolis step for acceptance and rejection, so
the acceptance rate is 100%. SGLD is slower than a componentwise
algorithm for two reasons: first it must read data from an external file
each iteration, and second, gradients with numerical differencing are
computationally expensive, requiring many model evaluations per
iteration. At least `Nr / size` iterations are suggested. SGLD has great
promise for the application of MCMC to massive data sets.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step-size schedule
    \\\epsilon_t\\, mini-batch size \\n\\, total data size \\N\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Read a random mini-batch of \\n\\ rows from the external data
        file into `Data$X`
    2.  Compute stochastic gradient: \\\hat{g}\_t = \nabla \log
        p(\theta\_{t-1}) + \frac{N}{n} \sum\_{i \in \text{batch}} \nabla
        \log p(x_i \mid \theta\_{t-1})\\
    3.  Draw noise \\\eta_t \sim \mathcal{N}(0, \epsilon_t \cdot I)\\
    4.  Update: \\\theta_t = \theta\_{t-1} + \frac{\epsilon_t}{2}
        \hat{g}\_t + \eta_t\\ (no accept/reject step)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 1 / Nr | Step-size schedule for stochastic gradient Langevin dynamics. |
| file | (supplied by user) | Mandatory: a CSV file with the rows of the big-data design matrix. |
| Nr | (supplied by user) | Mandatory: number of rows of the big-data file. |
| Nc | (supplied by user) | Mandatory: number of columns. |
| size | (supplied by user) | Mandatory: rows read per iteration. |

Default `Specs` for SGLD. {.table}

#### Example with default specifications

``` r

## SGLD requires user-supplied specifications that have no sensible
## defaults. See the table above and the help page ?lucifer for details.
```

#### See Also

- [HMC](#hmc)

### Stochastic Gradient Hamiltonian Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 100%, since SGHMC does not include a Metropolis accept/reject step. |
| Applications | This algorithm is designed for big data problems where the full dataset does not fit in memory. Data are streamed from an external file in mini-batches. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has seven algorithm specifications, and requires the user to save data to a file and specify file dimensions. Tuning the step size and friction coefficient requires care. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Stochastic Gradient Hamiltonian Monte Carlo (SGHMC) was proposed by
[\[43\]](#ref43) as an extension of Stochastic Gradient Langevin
Dynamics (SGLD) that incorporates second-order momentum dynamics. Where
SGLD uses only gradient information, discarding curvature, SGHMC
introduces a momentum variable and a friction term that corrects for the
noise introduced by stochastic gradient estimates. The resulting
dynamics closely resemble those of Hamiltonian Monte Carlo (HMC), but
applied to mini-batches rather than the full dataset. This makes SGHMC
suitable for large-scale Bayesian inference where evaluating the full
likelihood at every iteration is prohibitively expensive.

SGHMC has seven algorithm specifications: `epsilon` is the step size
(learning rate), `alpha` is the friction coefficient in (0,1) that
controls the balance between momentum persistence and noise injection,
`file` is the path to the external data file in CSV format, `Nr` is the
number of rows in the file, `Nc` is the number of columns in the file,
`size` is the mini-batch size (number of rows to read per iteration),
and `L` is the number of leapfrog steps per iteration. The data file
must be in CSV format. Each iteration, a random contiguous block of
`size` rows is read from the file and placed in `Data[["X"]]` for model
evaluation. At least `Nr / size` iterations are suggested. Larger values
of `L` produce longer trajectories but increase computational cost
proportionally.

The momentum update follows the Langevin diffusion with friction:
\\r\_{t+1} = (1 - \alpha) r_t + \epsilon \nabla \log p(\theta \| x) +
\mathcal{N}(0, 2\alpha\epsilon)\\, where \\r\\ is the momentum,
\\\alpha\\ is the friction coefficient, and \\\epsilon\\ is the step
size. The noise term compensates for the stochastic gradient noise, so
that the stationary distribution remains the correct posterior. Unlike
standard HMC, SGHMC does not require a Metropolis correction step,
yielding a 100% acceptance rate.

The advantage of SGHMC over SGLD is improved mixing through momentum,
which allows the sampler to traverse the parameter space more
efficiently. Compared to full-data HMC, the computational cost per
iteration is dramatically lower when the dataset is large. The
disadvantage is that the friction coefficient and step size require
careful tuning, and the stochastic gradient noise is only approximately
corrected; in practice this means SGHMC is asymptotically exact but may
exhibit bias in finite samples, particularly with large step sizes.
Since SGHMC is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, momentum \\r_0 = 0\\, step
    size \\\epsilon\\, friction \\\alpha\\, mini-batch size \\n\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Read a random mini-batch of \\n\\ rows from the external data
        file
    2.  Compute stochastic gradient: \\\hat{g}\_t = \nabla \log
        p(\theta\_{t-1}) + \frac{N}{n} \sum\_{i \in \text{batch}} \nabla
        \log p(x_i \mid \theta\_{t-1})\\
    3.  **For** \\\ell = 1, \ldots, L\\ leapfrog steps:
        - \\r \leftarrow (1 - \alpha)\\ r + \epsilon\\ \hat{g}\_t +
          \mathcal{N}(0, 2\alpha\epsilon\\ I)\\
        - \\\theta \leftarrow \theta + \epsilon\\ r\\
    4.  Set \\\theta_t = \theta\\ (no accept/reject step)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | (supplied by user) | Mandatory: step size of the stochastic gradient HMC trajectory. |
| alpha | (supplied by user) | Mandatory: friction coefficient. |
| file | (supplied by user) | Mandatory: big-data CSV file. |
| Nr | (supplied by user) | Mandatory. |
| Nc | (supplied by user) | Mandatory. |
| size | (supplied by user) | Mandatory. |
| L | (supplied by user) | Mandatory: number of leapfrog steps. |

Default `Specs` for SGHMC. {.table}

#### Example with default specifications

``` r

## SGHMC requires user-supplied specifications that have no sensible
## defaults. See the table above and the help page ?lucifer for details.
```

#### See Also

- [HMC](#hmc)
- [SGLD](#sgld)
- [MCHMC](#mchmc)

### Barker Proposal

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is 40% (Livingstone and Zanella 2022). The observed acceptance rate may be suitable in the interval \[25%,55%\]. |
| Applications | This is a widely applicable, gradient-based algorithm that serves as a robust alternative to MALA. It is particularly useful when the optimal step size is unknown, because its acceptance rate degrades gracefully with step size misspecification, unlike MALA which collapses catastrophically. |
| Difficulty | This algorithm is easy for a beginner. It has three algorithm specifications, and the step size is adapted automatically during warmup. |
| Final Algorithm? | User Discretion. The Barker proposal may be used as a final algorithm after the warmup period. |
| Proposal | Multivariate. |

The Barker proposal was introduced by [\[44\]](#ref44) as a
gradient-informed MCMC algorithm that achieves the same \\O(d^{-1/3})\\
optimal scaling as the Metropolis-Adjusted Langevin Algorithm (MALA)
while being as robust to step size misspecification as random-walk
Metropolis (RWM). The proposal draws each component independently:
\\\theta_j^\* = \theta_j + \|\varepsilon_j\| \cdot s_j\\, where
\\\|\varepsilon_j\| \sim \text{Exp}(1/\sigma_j)\\ is the magnitude and
\\s_j \in \\-1, +1\\\\ is a random sign drawn from
\\\text{Bernoulli}(\text{sigmoid}(\varepsilon_j \cdot \partial \log \pi
/ \partial \theta_j))\\. This construction creates a skewed proposal
that tends to move in the direction of the gradient but with exponential
(heavy) tails, making it far more robust than MALA’s Gaussian proposal.
When the step size is too large, MALA’s acceptance rate drops to zero
exponentially fast; the Barker proposal’s acceptance rate degrades only
polynomially.

Barker has three algorithm specifications: `epsilon` is the step size
(scalar or vector of length \\d\\; if NULL, it is initialized
automatically), `A` is the number of warmup iterations for step size
adaptation (default \\\min(500, T/2)\\), and `delta` is the target
acceptance rate (default 0.40). During warmup, the step size is adapted
using Nesterov dual-averaging (the same scheme used by NUTS and HMCDA)
to achieve the target acceptance rate.

The Metropolis-Hastings acceptance ratio accounts for the asymmetry of
the Barker proposal. The forward proposal probability includes the
sigmoid-weighted direction, so the acceptance ratio contains correction
terms involving the sigmoid function evaluated at both the current and
proposed gradients.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\sigma\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute gradient \\g = \nabla \log \pi(\theta\_{t-1})\\
    2.  **For** each dimension \\j\\: draw \\\varepsilon_j \sim
        \text{Exp}(1/\sigma_j)\\; draw \\u_j \sim \text{Uniform}(0,1)\\;
        if \\u_j \< \text{sigmoid}(\varepsilon_j \cdot g_j)\\, set \\d_j
        = +1\\, else \\d_j = -1\\; propose \\\theta_j^\* = \theta_j +
        d_j \cdot \varepsilon_j\\
    3.  Compute gradient at proposal \\g^\* = \nabla \log
        \pi(\theta^\*)\\
    4.  Compute \\\log \alpha = \bigl\[\log \pi(\theta^\*) - \log
        \pi(\theta\_{t-1})\bigr\] + \sum_j \bigl\[\log
        \text{sigmoid}(-d_j g^\*\_j \varepsilon_j) - \log
        \text{sigmoid}(d_j g_j \varepsilon_j)\bigr\]\\
    5.  With probability \\\min(1, \exp(\log \alpha))\\, accept;
        otherwise retain
    6.  If \\t \leq A\\: adapt \\\sigma\\ via dual-averaging toward
        target \\\delta\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | NULL | Step size auto-initialized via find-reasonable-epsilon. |
| A | min(500, floor(Iterations/2)) (= 500) | Adaptation window. |
| delta | 0.40 | Target acceptance rate; 0.40 is the Barker optimum reported by Livingstone and Zanella (2022). |

Default `Specs` for Barker. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Barker", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [NUTS](#nuts)
- [HMCDA](#hmcda)

### Non-Reversible Overdamped Langevin

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | There is no acceptance step; all proposals are accepted. The asymptotic bias from discretization vanishes as the step size decreases. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is particularly useful when non-reversibility can accelerate convergence, such as for targets with elongated level sets where the vorticity drives circulation around the contours. |
| Difficulty | This algorithm is easy for a beginner. It has two algorithm specifications. There is no Metropolis correction step, so the output has asymptotic bias controlled by the step size. |
| Final Algorithm? | User Discretion. The algorithm introduces discretization bias, so it should be used with small step sizes or as a fast exploratory tool. |
| Proposal | Multivariate. |

Non-Reversible Overdamped Langevin was introduced by [\[45\]](#ref45) as
a modification of the overdamped Langevin diffusion that adds a
divergence-free (antisymmetric) drift term to break reversibility and
accelerate convergence. The standard overdamped Langevin diffusion
\\d\theta = \frac{1}{2}\nabla \log \pi(\theta)\\dt + dW_t\\ satisfies
detailed balance; adding an antisymmetric component \\J\nabla \log
\pi(\theta)\\ where \\J = -J^\top\\ produces \\d\theta = \frac{1}{2}(I +
J)\nabla \log \pi(\theta)\\dt + dW_t\\. Because \\\nabla \cdot (J \nabla
\log \pi \cdot \pi) = 0\\ for antisymmetric \\J\\, the stationary
distribution \\\pi(\theta)\\ is unchanged, but the dynamics are no
longer reversible. The vorticity induced by \\J\\ creates circulation
around the contours of \\\pi\\, which can dramatically accelerate mixing
when the target has elongated or multimodal structure. The discretized
version uses an Euler-Maruyama scheme without a Metropolis correction.

Non-Reversible Overdamped Langevin has two algorithm specifications:
`epsilon` is the step size, and `vorticity` controls the magnitude of
the antisymmetric matrix \\J\\ (default 1.0). The matrix \\J\\ is
constructed as \\J = \omega \cdot (A - A^\top)\\ where \\A\\ is a fixed
random matrix and \\\omega\\ is the vorticity parameter.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    antisymmetric matrix \\J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute gradient \\g = \nabla \log \pi(\theta\_{t-1})\\
    2.  Compute drift: \\\mu = \theta\_{t-1} + \frac{\epsilon}{2}(I +
        J)g\\
    3.  Draw \\\theta_t = \mu + \sqrt{\epsilon}\\\xi\\, where \\\xi \sim
        \mathcal{N}(0, I_d)\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.01 / sqrt(LIV) (~ 0.0058) | Conservative step size for the overdamped Langevin discretisation. |
| vorticity | 1.0 | Magnitude of the antisymmetric drift; 1.0 generates moderate non-reversibility. |

Default `Specs` for NROLangevin. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "NROLangevin", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [SGLD](#sgld)
- [Barker Proposal](#barker)

## Slice samplers avoid tuning by sampling uniformly under the density

Slice sampling introduces an auxiliary variable that defines a
horizontal “slice” through the target density, then samples uniformly
from the region where the density exceeds this level. The approach
guarantees acceptance (the acceptance rate is always 1) and requires no
proposal variance tuning, since the slice boundaries adapt automatically
to the local scale of the distribution. Extensions in this family handle
multivariate targets through factor rotations (AFSS), elliptical
geometry for Gaussian priors (ESS), oblique hyperrectangles for
correlated parameters (OHSS), reflections for bounded domains (RSS), or
eigenvector-aligned coordinates (UESS).

### Slice Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner. There are several algorithm specifications, and the defaults are recommended for automatic use with continuous parameters. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The origin of slice sampling dates back to [\[46\]](#ref46), and the
current algorithm was introduced in [\[47\]](#ref47) and improved in
[\[15\]](#ref15). In slice sampling, a distribution is sampled by
sampling uniformly from the region under the plot of its density
function. Here, slice sampling uses two phases as follows. First, an
interval is created for the slice, and second, rejection sampling is
performed within this interval.

Slice has five algorithm specifications: `B` defaults to NULL in which
case blockwise updating is not performed, or `B` accepts a list in which
each component contains a vector that indicates the position number of
parameters to be sampled per block. Within each block, the order of
evaluations is randomized each iteration. The `Bounds` specification
accepts either a vector of length two for the lower and upper bound, or
a list of vectors, one for each block. The `m` specification is the
maximum number of steps (an integer in \[1, ∞\], and defaults to Inf),
and may be a scalar or a list with a scalar per block. The `Type`
specification accepts the following strings: “`Continuous`”,
“`Nominal`”, or “`Ordinal`”. This string may be a scalar, or a list of
scalars, one for each block. This specification refers to the type of
parameter, which may be continuous or discrete. Discrete slice sampling
is performed either for nominal or ordinal parameters. Finally, the `w`
specification accepts a scalar step size (in (0, ∞), and defaults to 1),
or accepts a list of scalars, one for each block. Ideally, `w` is 3
standard deviations of the target distribution for continuous
parameters, and is usually 1 for discrete parameters. All parameters in
a block receive the same specifications, but specifications may differ
by block.

The Slice algorithm implemented here is componentwise, and the algorithm
for continuous parameters is based on figures 3 and 5 in
[\[15\]](#ref15), in which the slice is replaced with an interval \\I\\
that contains most of the slice. For each parameter, an interval \\I\\
is created and expanded by the “stepping out” procedure with step size
`w` until the interval is larger than the slice, and the number of steps
`m` may be limited by the user. The original slice sampler is
inappropriate for the unbounded interval (-∞, ∞), and this improved
version allows this unbounded interval by replacing the slice with the
created interval \\I\\arkov property.. This algorithm adaptively changes
the scale, though it is not an adaptive algorithm in the sense that it
retains the Markov property.

Blockwise sampling of a componentwise algorithm allows different
specifications per block, and allows the user to control the order of
evaluations, because blocks are processed in order, though parameters
per block are selected in a randomized order.

The lower and upper bounds of interval \\I\\ default to (-∞, ∞), and may
be constrained for constrained parameters. Once the interval is
constructed, a proposal is drawn randomly from within the interval until
the proposal is accepted, and the interval is otherwise shrunk. The
acceptance rate of this Slice algorithm is 1, though multiple model
evaluations occur per iteration.

When this Slice sampler uses the default specifications and begins far
from the target distributions, the time per iteration should decrease as
the algorithm approaches the target distributions. Considerable time may
be spent in the first iterations. One strategy may be to limit `m`.

This componentwise Slice algorithm has been noted to be more efficient
than Metropolis updates, may be easier to implement than Gibbs sampling
(Gibbs), and is attractive for routine and automated use
[\[15\]](#ref15). An adaptive version is AFSS. As a componentwise
algorithm, Slice is slower per iteration than algorithms that use
multivariate proposals. Multivariate slice samplers have been proposed,
such as ESS, OHSS, and UESS. Since Slice is not an adaptive algorithm,
it is acceptable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\w\\, maximum
    steps \\m\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order within each block:
        - Draw slice level \\y \sim \text{Uniform}\bigl(0,\\ p(\theta_t
          \mid \text{data})\bigr)\\
        - **Step out**: set \\L = \theta\_{t,j} - w \cdot u\\, \\R = L +
          w\\ (\\u \sim \text{Uniform}(0,1)\\); expand \\L\\ leftward
          and \\R\\ rightward by \\w\\ while \\p(\theta \mid
          \text{data}) \> y\\ at boundaries, up to \\m\\ steps
        - **Shrink and sample**: draw \\\theta^\*\_j \sim
          \text{Uniform}(L, R)\\; if \\p(\theta^\* \mid \text{data}) \>
          y\\, accept \\\theta\_{t+1,j} = \theta^\*\_j\\; otherwise
          shrink the interval toward \\\theta\_{t,j}\\ and repeat
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| B | list(1:LIV) | Single-block update over all parameters. |
| Bounds | c(-Inf, Inf) | Unbounded slices, appropriate for unconstrained parameters. |
| m | Inf | No upper limit on the stepping-out interval; slice widening continues until the level set is bracketed. |
| Type | “Continuous” | Continuous parameters. |
| w | 1 | Initial slice width. |

Default `Specs` for Slice. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Slice", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AFSS](#afss)
- [ESS](#ess)
- [Gibbs](#gibss)
- [OHSS](#ohss)
- [RSS](#rss)
- [UESS](#uess)

### Automated Factor Slice Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner. There are several algorithm specifications, and the defaults are recommended for automatic use. |
| Final Algorithm? | Yes, when `A`=0. See below. |
| Proposal | Multivariate proposals with componentwise evaluation. |

The Automated Factor Slice Sampler of [\[48\]](#ref48) is an extension
of the componentwise slice sampler in which linear correlation between
parameters is accounted for and the interval widths are tuned. This
results in a componentwise slice sampler that is faster than a
traditional componentwise slice sampler because less evaluations are
necessary to estimate the slice boundaries, and it samples efficiently
in high-dimensional models with continuous parameters.

AFSS has five algorithm specifications: `A` accepts a scalar that
indicates the iteration in which adaptation will cease, and defaults to
Inf. `B` is for blockwise sampling, defaults to NULL, and otherwise
accepts a list of parameter positions, in which each list component is a
block. Blocks are processed sequentially, but the order of parameters is
selected randomly within each block, and each block must contain at
least two parameters. The `m` specification accepts either a scalar
integer or a vector equal in length to the number of parameters that
indicates the maximum number of steps. It defaults to 100 and must be in
\[1, ∞\]. The `n` specification accepts a scalar that indicates the
total previous number of adaptive iterations, and defaults to 0.
Finally, w accepts either a scalar or a vector equal in length to the
number of parameters, and is initial step size, which defaults to 1 and
must be in (0, ∞).

AFSS performs two kinds of adaptation. While adaptive, each iteration
AFSS collects an approximation of the covariance of the parameters via a
scatter matrix. Periodically, eigenvectors are estimated from the
approximated covariance matrix. More parameters results in less frequent
adaptation. Each iteration, componentwise evaluation is performed as per
a componentwise sampler. However, after the first eigenvalue estimation,
all parameters are updated for each componentwise evaluation, which
results in a sampler that is not truly multivariate or componentwise.
Each parameter is updated in turn, but all other parameters are also
updated due to linear correlations with the current parameter. Also,
after each estimation of eigenvalues, the interval widths are tuned
toward an equal probability of interval expansion during the step-out
phase and interval shrinkage during rejection sampling. As the number of
parameters increases, so does the computational cost of estimating
eigenvalues, and blockwise sampling becomes more practical in
high-dimensional models.

While AFSS is efficient with linear relationships, it may be inefficient
due to nonlinear relationships. The AFSS algorithm is slower per
iteration than multivariate samplers, but typically improves more per
iteration than most algorithms. To date, testing indicates AFSS is
remarkably efficient. AFSS is acceptable as a final algorithm when
non-adaptive samples are drawn.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step sizes \\w_j\\, scatter
    matrix \\S = 0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Update scatter matrix: \\S \leftarrow S + \theta_t
        \theta_t^\top\\
    2.  Periodically estimate eigenvectors \\V\\ from \\S\\ (less
        frequently as \\J\\ increases)
    3.  **For** each parameter \\j\\ in random order within each block:
        - Draw slice level \\y \sim \text{Uniform}\bigl(0,\\ p(\theta_t
          \mid \text{data})\bigr)\\
        - **Step out**: expand interval \\I\\ by \\w_j\\ until both ends
          fall below \\y\\, limited by \\m\\ steps
        - **Shrink and sample**: draw \\\theta^\*\_j\\ uniformly from
          \\I\\; if \\p(\theta^\* \mid \text{data}) \< y\\, shrink \\I\\
          and repeat
        - If eigenvectors available, update all parameters via linear
          transformation: \\\theta \leftarrow \theta + V \Delta\\
    4.  If \\t \leq A\\, tune \\w_j\\ toward equal step-out and
        shrinkage probabilities
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | Inf | Disable forced step-size truncation. |
| B | NULL | Single-block update. |
| m | Inf (per-coordinate) | No upper limit on the slice-widening loop. |
| n | 0 | Adaptation history starts empty. |
| w | 1 (per-coordinate) | Unit initial slice width per factor. |

Default `Specs` for AFSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AFSS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice](#slice)

### Elliptical Slice Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This algorithm is applicable only to models in which the prior mean of all parameters is zero. |
| Difficulty | This algorithm is easy for a beginner. There is only one optional algorithm specification for block updating, and tuning is unnecessary, though optimal performance is gained with a user-specified prior covariance matrix for the ellipse. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Elliptical Slice Sampling (ESS) was introduced by [\[49\]](#ref49) as a
multivariate extension of componentwise Slice sampling (Slice) that
utilizes ellipse ν and angle θ. Each iteration, ellipse ν is drawn from
N(0, Σ), where Σ is a user-specified prior covariance for the ellipse.
The authors recommend using a form of the covariance of the data. Even
though other parameters are not discussed, an identity matrix, or a
combination thereof, performs well in practice. The prior covariance
matrix is accepted as `VarCov` in the lucifer function, rather than the
usual proposal covariance matrix.

ESS has one algorithm specification `B` for blocks of parameters, and
defaults to `B=NULL`. For large-dimensional problems, block updating may
be used. To specify block updating, the user gives the `B` specification
a list in which each component is a block, and each component is a
vector of integers pointing to the position of the parameters for that
block. Block updating also requires the `Covar` argument to receive a
list of prior covariance matrices of the appropriate dimension.

Each proposal is an additive, trigonometric function of the current
state, ellipse ν, and angle θ. Angle θ is originally in the interval (0,
2π\] each iteration, and is shrunk until acceptance occurs. This results
in a 100% acceptance rate, and usually multiple model evaluations per
iteration.

This algorithm is applicable only to models in which the prior mean of
all parameters is zero. It is often possible to apply ESS when the
variables are centered and scaled.

An advantage of ESS over the (componentwise) Slice algorithm is the
computational efficiency of a multivariate proposal. A disadvantage is
that the user must specify the prior covariance Σ for optimal ESS
performance, and the algorithm can be very sensitive to this prior.
Since ESS is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, prior covariance \\\Sigma\\
    for the ellipse
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw ellipse \\\nu \sim \mathcal{N}(0, \Sigma)\\
    2.  Set log-likelihood threshold: \\\log y = \log p(\theta_t \mid
        \text{data}) + \log u\\, where \\u \sim \text{Uniform}(0,1)\\
    3.  Draw angle \\\vartheta \sim \text{Uniform}(0, 2\pi)\\; set
        bracket \\\[\vartheta\_{\min}, \vartheta\_{\max}\] =
        \[\vartheta - 2\pi, \vartheta\]\\
    4.  **Loop**:
        - Propose \\\theta^\* = \theta_t \cos \vartheta + \nu \sin
          \vartheta\\
        - If \\\log p(\theta^\* \mid \text{data}) \> \log y\\, accept
          \\\theta\_{t+1} = \theta^\*\\ and **break**
        - Otherwise, shrink the bracket around \\\vartheta\\ and draw a
          new \\\vartheta\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification                                |
|:--------------|:--------|:---------------------------------------------|
| B             | list()  | Single elliptical slice over all parameters. |

Default `Specs` for ESS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "ESS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice](#slice)

### Oblique Hyperrectangle Slice Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy to use and does not require tuning. |
| Final Algorithm? | Yes, when it is non-adaptive. |
| Proposal | Multivariate. |

Oblique Hyperrectangle Slice Sampler (OHSS) is an adaptive algorithm
that approximates the sample covariance matrix. Introduced by
[\[50\]](#ref50) where it is referred to as “Oblique Hyperrect”, OHSS is
an extension of the hyperrectangle method in [\[15\]](#ref15), which in
turn is a multivariate extension of the slice sampler (Slice). The
initial hyperrectangle is oriented to have edges that are parallel to
the eigenvectors of the sample covariance matrix. The initial
hyperrectangle is also scaled to have lengths proportional to the square
roots of the eigenvalues of the sample covariance matrix. Rejection
sampling is performed, and when a sample is rejected, the slice
approximation is shrunk. Aside from specifying the initial
hyperrectangle, this is the hyperrectangle method in [\[15\]](#ref15).

OHSS has two algorithm specifications: `A` is the number of adaptive
iterations, `n` is the number of previous iterations, and defaults to
zero. The number of previous iterations, if any, is used to weight the
sample covariance matrix.

Each iteration, there is a 5% probability that a non-adaptive step is
taken, rather than an adaptive step. The first ten iterations or so are
non-adaptive. Eigenvectors are estimated no more than every tenth
iteration. When adaptive, the sample covariance matrix is updated along
with the eigenvectors.

Once the slice interval is estimated, a sample is drawn uniformly with
rejection sampling from within this interval. If rejected, then the
interval is shrunk until an acceptable sample is drawn. OHSS has a 100%
acceptance rate.

The time per iteration varies, since rejection sampling often requires
more than one evaluation. When OHSS is adaptive, it is not suitable as a
final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, sample covariance \\\Sigma =
    I\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  If adaptive (with 95% probability) and \\t \> 10\\: update
        \\\Sigma\\ from sample history; estimate eigenvectors \\V\\ and
        eigenvalues \\\Lambda\\ periodically
    2.  Draw slice level \\y \sim \text{Uniform}\bigl(0,\\ p(\theta_t
        \mid \text{data})\bigr)\\
    3.  Construct initial hyperrectangle \\H\\ with edges along
        eigenvectors of \\\Sigma\\, lengths \\\propto \sqrt{\Lambda}\\
    4.  **Rejection sampling**: draw \\\theta^\*\\ uniformly from \\H\\
        - If \\p(\theta^\* \mid \text{data}) \> y\\, accept
          \\\theta\_{t+1} = \theta^\*\\
        - Otherwise, shrink \\H\\ and repeat
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | Iterations + 1 | Disable adaptation by setting the warm-up window beyond the run length. |
| n | 0 | Adaptation history starts empty. |

Default `Specs` for OHSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "OHSS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice](#slice)
- [UESS](#uess)

### Reflective Slice Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is 1. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is difficult to tune. There are two tuning parameters: the number of steps to take, and the step size. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Introduced in [\[47\]](#ref47), the Reflective Slice Sampler (RSS) is a
multivariate slice sampler that uses partial derivatives and reflects
into a slice interval. As described in [\[15\]](#ref15), each iteration,
direction and momentum within the slice interval are determined
randomly, partial derivatives are taken once, and a number of steps is
taken according to a step size. The distance of each step is the product
of the step size and momentum. Direction is changed by reflecting off
the boundaries of the slice interval. The reflected momentum direction
is a function of the incident momentum and the partial derivatives.
Random-walk behavior is suppressed, because a large number of steps may
be taken in the same direction. The acceptance rate is 100%.

RSS has two algorithm specifications: `m` and `w`. Each iteration, `m`
steps are taken with step size `w`. While `m` must be a scalar, `w` may
be a scalar or vector equal in length to the number of parameters. A
step size that is too small is inefficient, but too small is better than
too large, which can entirely avoid the target distribution.

RSS is difficult to tune, but is consistent in time per iteration.
[\[15\]](#ref15) states that Hamiltonian Monte Carlo (HMC) performs
better than RSS when both HMC and RSS are optimally tuned. Since RSS is
not adaptive, it may be used as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, number of steps \\m\\, step
    size \\w\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw slice level \\y \sim \text{Uniform}\bigl(0,\\ p(\theta_t
        \mid \text{data})\bigr)\\
    2.  Compute gradient \\g = \nabla \log p(\theta_t \mid y)\\; draw
        random momentum \\r\\
    3.  **For** \\s = 1, \ldots, m\\:
        - Move \\\theta^\* \leftarrow \theta^\* + w \cdot r\\ (step in
          current direction)
        - If \\p(\theta^\* \mid \text{data}) \< y\\ (outside the slice),
          **reflect**: update \\r\\ using gradient \\g\\ at boundary
    4.  Accept \\\theta\_{t+1} = \theta^\*\\ (100% acceptance)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| m | 100 | Mandatory: maximum number of reflections per iteration; 100 is enough for low-dimensional targets. |
| w | rep(1, LIV) | Mandatory: per-coordinate initial widths of the reflection box. |

Default `Specs` for RSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RSS", Specs = list(m = 100L, w = rep(1, length(IV))),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [Slice](#slice)

### Univariate Eigenvector Slice Sampler

| Aspect           | Description                                             |
|:-----------------|:--------------------------------------------------------|
| Acceptance Rate  | The acceptance rate is 1.                               |
| Applications     | This is a widely applicable, general-purpose algorithm. |
| Difficulty       | This algorithm is easy to use.                          |
| Final Algorithm? | Yes, when non-adaptive.                                 |
| Proposal         | Multivariate.                                           |

The Univariate Eigenvector Slice Sampler (UESS) is an adaptive algorithm
that approximates the sample covariance matrix in the same manner as
Adaptive-Mixture Metropolis (AMM). Described as “Univar Eigen” in
[\[50\]](#ref50), UESS is a multivariate variation of slice sampling
(Slice) that makes univariate updates along the eigenvectors of the
sample covariance matrix.

UESS has four algorithm specifications: `A` is the number of adaptive
iterations, `B` accepts a list of blocks and defaults to NULL, `m` is
the number of steps, and `n` is the number of previous iterations, and
defaults to zero. Each iteration, a slice interval is estimated with up
to `m` steps. The default is `m=100` steps. The step size is affected by
the eigenvectors of the sample covariance matrix. The number of previous
iterations, if any, is used to weight the sample covariance matrix.

The contours of the target distribution are approximately ellipsoidal
and the eigenvectors of its covariance coincide roughly with the axes of
these ellipsoids, provided that the shape of the target distribution is
approximately convex. Steps taken along these eigenvectors include steps
along the long axes of a slice.

Each iteration, there is a 5% probability that a non-adaptive step is
taken, rather than an adaptive step. The first ten iterations or so are
non-adaptive. Eigenvectors are estimated no more than every tenth
iteration. When adaptive, the sample covariance matrix is updated each
iteration.

Once the slice interval is estimated, a sample is drawn uniformly with
rejection sampling from within this interval. If rejected, then the
interval is shrunk until an acceptable sample is drawn. UESS has a 100%
acceptance rate.

The time per iteration varies, since building the slice area requires up
to `m` steps, and rejection sampling often requires more than one
evaluation. Although UESS is a combination of the AMM and Slice
algorithms, it usually performs as well or better than either. When UESS
is adaptive, it is not suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, sample covariance \\\Sigma =
    I\\, step count \\m\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  If adaptive (with 95% probability) and \\t \> 10\\: update
        \\\Sigma\\ from sample history; estimate eigenvectors \\V\\ and
        eigenvalues \\\Lambda\\
    2.  **For** each eigenvector direction \\v_j\\:
        - Draw slice level \\y \sim \text{Uniform}\bigl(0,\\ p(\theta_t
          \mid \text{data})\bigr)\\
        - **Step out** along \\v_j\\: expand interval by steps scaled to
          \\\sqrt{\Lambda_j}\\, up to \\m\\ steps
        - **Shrink and sample**: draw \\\theta^\*\\ uniformly along
          \\v_j\\ within the interval; if \\p(\theta^\* \mid
          \text{data}) \> y\\, accept; otherwise shrink and repeat
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| A | Inf | Disable forced truncation of the slice-widening loop. |
| B | NULL | Single block. |
| m | 100 | Maximum number of stepping-out iterations per coordinate. |
| n | 0 | Adaptation history starts empty. |

Default `Specs` for UESS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "UESS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMM](#amm)
- [Slice](#slice)

### Quantile Slice Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is effectively 100% because slice sampling always accepts. |
| Applications | This algorithm is applicable to univariate full conditionals within a Gibbs scheme. It is particularly useful when the target density has heavy tails or unbounded support, because the quantile transformation maps the slice interval to a bounded domain. |
| Difficulty | This algorithm is easy for a beginner. It has one algorithm specification. |
| Final Algorithm? | Yes. The Quantile Slice Sampler may be used as a final algorithm. |
| Proposal | Univariate (componentwise). |

The Quantile Slice Sampler (QSS) was introduced by [\[51\]](#ref51) as a
slice sampling method that uses a quantile transformation to convert the
stepping-out and shrinking procedures of Neal (2003) into a bounded
interval. Standard slice sampling on unbounded distributions requires
stepping-out procedures that can be slow when the tails are heavy; QSS
instead maps the real line to \\(0, 1)\\ via the CDF of a reference
distribution (typically Cauchy), performs slice sampling on the
transformed scale, and maps back. The Cauchy reference is heavy-tailed
enough to handle most targets without excessive stepping out. On the
transformed scale, the slice interval is always within \\\[0, 1\]\\, and
the stepping-out procedure is replaced by a simple shrinkage scheme on
this bounded interval.

QSS has one algorithm specification: `w` is the scale parameter of the
Cauchy reference distribution (default 1.0). Larger values spread the
quantile mapping, which is appropriate when the full conditional has
wider tails.

#### Algorithm

1.  **For** each component \\j\\ (within Gibbs):
    1.  Draw vertical level \\y \sim \text{Uniform}(0, \pi(\theta_j \mid
        \theta\_{-j}))\\
    2.  Map current value to quantile scale: \\u =
        F\_{\text{Cauchy}}(\theta_j; w)\\
    3.  Draw interval \\\[L, R\]\\ uniformly containing \\u\\ within
        \\\[0, 1\]\\
    4.  **Repeat** (shrinkage):
        - Draw \\u^\* \sim \text{Uniform}(L, R)\\
        - Map back: \\\theta_j^\* = F\_{\text{Cauchy}}^{-1}(u^\*; w)\\
        - If \\\pi(\theta_j^\* \mid \theta\_{-j}) \> y\\: accept;
          **break**
        - Shrink \\\[L, R\]\\ around \\u\\
2.  **Return** updated \\\theta\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| w | 1.0 | Initial width of the quantile slice; the algorithm is insensitive to this choice for unit-scale targets. |

Default `Specs` for QSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "QSS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice Sampler](#slice)
- [Latent Slice Sampler](#lss)
- [Gibbsian Polar Slice Sampler](#gpss)

### Latent Slice Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is effectively 100% because slice sampling always accepts. |
| Applications | This algorithm is applicable to general targets. It eliminates the need for stepping-out or doubling procedures by introducing a latent shrinkage variable that controls the slice width adaptively. |
| Difficulty | This algorithm is easy for a beginner. It has one algorithm specification, and the latent variable s automatically adapts the slice width. |
| Final Algorithm? | Yes. The Latent Slice Sampler may be used as a final algorithm. |
| Proposal | Multivariate. |

The Latent Slice Sampler (LSS) was introduced by [\[52\]](#ref52) as a
slice sampling method that avoids the stepping-out and shrinking
procedures of Neal (2003) by augmenting the state space with a latent
“shrinkage” variable \\s \> 0\\. At each iteration, the joint target is
\\\pi(\theta, s) \propto \pi(\theta) \cdot \mathbf{1}\[s \<
g(\theta)\]\\ for some monotonically related function \\g\\, and the
update alternates between sampling \\s\\ given \\\theta\\ (a uniform
draw) and \\\theta\\ given \\s\\ (a uniform draw from the slice
\\\\\theta : g(\theta) \> s\\\\). The key insight is that \\s\\ acts as
an adaptive slice width: when the current \\\theta\\ is in a
high-density region, \\s\\ can be large, narrowing the slice; when
\\\theta\\ is in a low-density region, \\s\\ is small, widening it. This
eliminates all tuning parameters associated with stepping-out procedures
while maintaining the exact sampling guarantee of slice sampling.

LSS has one algorithm specification: `s_init` is the initial value of
the latent shrinkage variable (default 1.0).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, shrinkage \\s_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw vertical level: \\s_t \sim \text{Uniform}(0,
        \pi(\theta\_{t-1}))\\
    2.  Draw direction \\d \sim \text{Uniform}(\mathcal{S}^{d-1})\\
        (unit sphere)
    3.  Draw step size \\\ell \sim \text{Uniform}(-s\_{t-1}, s\_{t-1})\\
    4.  Propose \\\theta^\* = \theta\_{t-1} + \ell \cdot d\\
    5.  If \\\pi(\theta^\*) \> s_t\\: accept \\\theta_t = \theta^\*\\,
        update \\s_t = \|\ell\|\\
    6.  Otherwise: shrink \\s_t\\ and repeat from (c)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| s_init | 1.0 | Initial value of the latent slice variable; updated automatically. |

Default `Specs` for LSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "LSS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice Sampler](#slice)
- [Quantile Slice Sampler](#qss)
- [Gibbsian Polar Slice Sampler](#gpss)

### Gibbsian Polar Slice Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is effectively 100% because slice sampling always accepts. |
| Applications | This algorithm is applicable to multivariate targets. It is particularly effective for targets with strong correlations or banana-shaped geometries, because the polar decomposition separates radial and angular components. |
| Difficulty | This algorithm is moderately easy for a beginner. It has two algorithm specifications. The polar coordinate transformation adds conceptual complexity but the implementation is straightforward. |
| Final Algorithm? | Yes. GPSS may be used as a final algorithm. |
| Proposal | Multivariate (polar coordinates). |

The Gibbsian Polar Slice Sampler (GPSS) was introduced by
[\[53\]](#ref53) as a slice sampling method that operates in polar
coordinates. The parameter vector \\\theta \in \mathbb{R}^d\\ is
decomposed into a radius \\r = \\\theta\\\\ and a direction \\\omega =
\theta / \\\theta\\ \in \mathcal{S}^{d-1}\\. The Gibbs update alternates
between updating the radius (given the direction) using standard
univariate slice sampling on \\\pi(r\omega) \cdot r^{d-1}\\, and
updating the direction (given the radius) using a slice sampling step on
the sphere. The directional update uses a stepping-out procedure along a
great circle, which naturally respects the spherical geometry. The polar
decomposition is advantageous because correlations in Cartesian
coordinates often manifest as smooth, easy-to-sample radial profiles
combined with concentrated directional distributions, rather than the
entangled structure that makes Cartesian slice sampling inefficient.

GPSS has two algorithm specifications: `w` is the initial slice width
for the radial stepping-out procedure (default 1.0), and `max_steps` is
the maximum number of stepping-out steps (default 32).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, decompose into \\r_0 =
    \\\theta_0\\\\, \\\omega_0 = \theta_0 / r_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Update radius** (given direction \\\omega\_{t-1}\\):
        - Draw vertical level \\y \sim \text{Uniform}(0, \pi(r\_{t-1}
          \omega\_{t-1}) \cdot r\_{t-1}^{d-1})\\
        - Find slice interval \\\[L, R\]\\ via stepping-out with width
          \\w\\
        - Draw \\r_t\\ uniformly from \\\\r \in \[L, R\] : \pi(r
          \omega\_{t-1}) \cdot r^{d-1} \> y\\\\ via shrinkage
    2.  **Update direction** (given radius \\r_t\\):
        - Draw random great circle through \\\omega\_{t-1}\\
        - Perform univariate slice sampling along the great circle with
          target \\\pi(r_t \omega)\\
        - Set \\\omega_t\\ to the accepted point
    3.  Reconstruct \\\theta_t = r_t \omega_t\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| w | 1.0 | Initial slice width along the polar direction. |
| max_steps | 32 | Cap on the number of stepping-out iterations: prevents runaway expansion in pathological cases. |

Default `Specs` for GPSS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "GPSS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Slice Sampler](#slice)
- [Quantile Slice Sampler](#qss)
- [Latent Slice Sampler](#lss)

## Ensemble and population methods use walker collectives to explore the target

Ensemble methods maintain a population of walkers (also called particles
or chains) that interact to generate proposals, inheriting global
information about the target geometry from the collective walker
positions without requiring gradient evaluations or explicit covariance
estimation. The affine-invariant ensemble sampler (AIES) uses stretch
moves between walker pairs; differential evolution (DEMC) uses vector
differences; the t-walk uses four complementary move types. The
Quadratic Monte Carlo (QMC) family, introduced by Militzer (2023, 2025),
generalizes stretch moves via Lagrange polynomial interpolation through
multiple guide walkers, enabling proposals that follow curved manifolds
rather than straight lines. See
[`vignette("qmc", package = "lucifer")`](https://robustecologies.github.io/lucifer/articles/qmc.md)
for a detailed treatment of the QMC methods.

### Affine-Invariant Ensemble Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is relatively easy for a beginner. The user must select the number of walkers, though the recommended default is suitable, must supply a vector of initial values for each walker, and specify a scale parameter, which again has a suggested default. If the user experiences difficulty, the recommendation is to increase the number of walkers. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

The Affine-Invariant Ensemble Sampler (AIES) of [\[55\]](#ref55) uses a
complementary ensemble of at least 2J walkers for J parameters. Each
walker receives J initial values, and initial values must differ for
each walker. At each iteration, AIES makes a multivariate proposal for
each walker, given a scaled difference in position by parameter between
the current walker and another randomly-selected walker.

AIES has six algorithm specifications: `Nc` is the number of walkers,
`Z` is a \\N_c\\ x \\J\\ matrix of initial values, `β` is a scale
parameter, `CPUs` is the number of central processing units (CPUs),
`Packages` is a vector of package names, and `Dyn.lib` is a vector of
shared libraries. The recommended number of walkers is at least \\2J\\.
If separate sets of initial values are not supplied in `Z`, since `Z`
defaults to NULL, then the GIV function is used to generate initial
values. The original article referred to the scale parameter as α,
though it has been renamed here to β to avoid a conflict with the
acceptance probability α in the Metropolis step. The β parameter may be
manipulated to affect the desired acceptance rate, though in practice,
the acceptance rate may potentially be better affected by increasing the
number of walkers. It is recommended to specify `CPUs=1` and leave the
remaining arguments to `NULL`, unless needed.

This version returns the samples from one walker, and the other walkers
assist the main walker. A disadvantage of this approach is that all
samples from all walkers are not returned. An advantage of this approach
is that if a particular walker is an outlier, then it does not affect
the main walker, unless of course it is the main walker! Multiple sets
of samples are best returned in a list, such as with parallel chains in
the lucifer.hpc function, though it is not applicable in lucifer.

AIES has been parallelized by Foreman-Mackey et al. (2012), and this
style of parallelization is available here as well. The user is
cautioned to prefer CPUs=1, because parallelizing may result in a slower
algorithm due to communication between the master and slave CPUs each
iteration. This communication is costly, and is best overcome with a
large number of CPUs available, and when the `Model` function is slow to
evaluate in comparison to network communication time.

Both the parallelized (CPUs \> 1) and un-parallelized (CPUs=1) versions
should be called from lucifer, not lucifer.hpc. When parallelized, the
number of walkers must be an even number (odd numbers are not
permitted), and the walkers are split into two equal-sized batches. Each
iteration, each walker moves in relation to a randomly-selected walker
in the other batch. This retains detailed balance.

AIES is attractive for offering affine-invariance, and therefore being
generally robust to badly scaled posteriors, such as with highly
correlated parameters. It is also attractive for making a multivariate
proposal without a proposal covariance matrix. However, since at least
2J walkers are recommended, the number of model evaluations per
iteration exceeds most componentwise algorithms by at least twice,
making AIES a slow algorithm per iteration, and computation scales
poorly with model dimension. Large-scale computing environments may
overcome this limitation with parallelization, but parallelization is
probably not very helpful (and may be detrimental) in small-scale
computing environments when evaluating the model function is not slow in
comparison with network communication time. AIES is not an adaptive
algorithm, and is therefore suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c\\ walkers with positions \\\\z_1, \ldots,
    z\_{N_c}\\\\, each of dimension \\J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\k = 1, \ldots, N_c\\:
        - Select a complementary walker \\j \neq k\\ uniformly at random
        - Draw scale factor \\Z \sim g(z) \propto 1/\sqrt{z}\\ for \\z
          \in \[1/\beta,\\ \beta\]\\
        - Propose \\\theta^\* = z_j + Z \cdot (z_k - z_j)\\
        - Compute \\\alpha = \min\\\bigl(1,\\ Z^{J-1} \cdot p(\theta^\*
          \mid y) / p(z_k \mid y)\bigr)\\
        - With probability \\\alpha\\, set \\z_k \leftarrow \theta^\*\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers in the ensemble: at least 2\*LIV+1 to avoid degeneracy, padded to 10 for variety. |
| Z | NULL | Walker positions auto-initialized via PGF. |
| beta | 2 | Stretch parameter; 2 is the Goodman/Weare default. |
| CPUs | 1 | Sequential evaluation of walkers. |
| Packages | NULL | No extra packages. |
| Dyn.libs | NULL | No dynamic libraries. |

Default `Specs` for AIES. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "AIES", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [twalk](#twalk)
- [QMC](#qmc)
- [DEMC](#demc)

### Differential Evolution Markov Chain

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is moderate in difficulty. It has few algorithm specifications that are easy to specify, but proficiency comes with practice. This algorithm requires an additional matrix of initial values, and the number of chains depends on the number of parameters. |
| Final Algorithm? | No. |
| Proposal | Multivariate. |

The original Differential Evolution Markov Chain (DEMC)
[\[56\]](#ref56), referred to in the literature as DE-MC, included a
Metropolis step on a genetic algorithm called Differential Evolution
with multiple chains (per parameter), in which the chains learn from
each other. It could be considered as parallel adaptive direction
sampling (ADS) with the Gibbs sampling (Gibbs) step replaced by a
Metropolis step, or as a non-parametric form of Random-Walk Metropolis
(RWM). However, for a model with J parameters, the original DEMC
required at least 2J chains, and often required as many as 20J chains.
Hence, from 2J to 20J model evaluations were required per iteration,
whereas a typical multivariate sampler such as Adaptive-Mixture
Metropolis (AMM) requires one evaluation, or a componentwise sampler
such as Adaptive Metropolis-within-Gibbs (AMWG) requires J evaluations
per iteration.

The version included here was presented in [\[57\]](#ref57), and the
required number of chains is drastically reduced by adapting based on
historical, thinned samples (the parallel direction move), and
periodically using a snooker move instead. In the article, three chains
were used to update as many as 50-100 parameters. Larger models may
require blocks (as suggested in the article, and blocking is not
included here), or more chains (see below). In testing, here, a few
200-dimensional models have been solved with 5-10 chains.

DEMC has four algorithm specifications: `Nc` is the number of chains (at
least 3), `Z` is a \\T\\ x \\J\\ matrix or \\T\\ x \\J\\ x \\N_c\\ array
of \\T\\ thinned samples for \\J\\ parameters and \\N_c\\ chains,
`gamma` is a scale parameter, and `w` is the probability of a snooker
move for each iteration. When `gamma=NULL`, the scale parameter defaults
to the recommended 2.381204 / sqrt(2J), though for snooker moves, it is
2.381204 / sqrt(2) regardless of the algorithm specification. The
default, recommended probability of a snooker move is `w = 0.1`.

The parallel direction move consists of a multivariate proposal for each
chain in which two randomly selected previous iterations from two
randomly selected other chains are differenced and scaled, and a small
jitter, `U ∼ (-0.001, 0.001)^J`, is added. The snooker move differences
these with another randomly selected (current) chain in the current
iteration, and with a fixed scale. Another variation is to use snooker
with past chain time-periods. The snooker move facilitates mode-jumping
behavior for multimodal posteriors.

For the first update, `Z` is usually NULL. Internally, lucifer populates
`Z` with GIV, and it is strongly recommended that PGF is specified by
the user. As the sampler iterates, `Z` is used for adaptation, and
elements of `Z` are replaced with thinned samples. A short, first run is
recommended to obtain thinned samples, such as in `Fit$Posterior1`. For
consecutive updates, this posterior matrix is used as `Z`.

In this implementation, samples are returned of the main chain, for
which `Initial.Values` are specified. Samples for other chains
(associated with PCIV) are not returned, but are used to assist the main
chain. The authors note that too many chains can be problematic when an
outlier chain exists. Here, samples of other chains are not returned,
outlier or not. If an outlier chain exists, it simply does not help the
main chain much and wastes computational resources, but does not
negatively affect the main chain.

An attractive property is that DEMC does not require a proposal
covariance matrix, but adapts instead based on past (thinned) samples.
In the output of lucifer, the thinned samples are also stored in
CovarDHis, though they are thinned samples, not the history of the
diagonal of the covariance matrix. This is done so the absolute
differences can be observed for diminishing adaptation. Another
attractive property is that the chains may be parallelized, such as
across CPUs, in the future, though the current version is not
parallelized.

DEMC has been considered to perform like a version of
Metropolis-within-Gibbs (MWG) that updates by chain, rather than by
component. DEMC may be one form of a compromise between a one-evaluation
multivariate proposal and J componentwise proposals. Since it is
adaptive, DEMC is not recommended as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ chains with positions
    \\\\\theta^{(1)}, \ldots, \theta^{(N_c)}\\\\, past samples matrix
    \\Z\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each chain \\k = 1, \ldots, N_c\\:
        - Draw \\u \sim \text{Uniform}(0,1)\\
        - If \\u \< w\\ (**snooker move**): select three random elements
          from \\Z\\ and current chains; propose along the line
          connecting two of them with fixed scale \\\gamma =
          2.38/\sqrt{2}\\
        - Otherwise (**parallel direction move**): select two random
          past samples \\z_a, z_b\\ from different chains in \\Z\\;
          propose \\\theta^\* = \theta^{(k)}\_t + \gamma \cdot (z_a -
          z_b) + \text{jitter}\\, where \\\gamma = 2.38/\sqrt{2J}\\
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta^{(k)}\_t \mid y)\bigr)\\ (with Jacobian correction
          for snooker)
        - With probability \\\alpha\\, accept \\\theta^{(k)}\_{t+1} =
          \theta^\*\\
    2.  Periodically store thinned samples in \\Z\\
3.  **Return** samples from the main chain

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers. |
| Z | NULL | Auto-initialized. |
| gamma | 2.38 / sqrt(2 \* LIV) | ter Braak (2006) optimal differential-evolution step size. |
| w | 0.05 | Mixture weight on the static Gaussian noise component. |

Default `Specs` for DEMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DEMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMM](#amm)
- [AMWG](#amwg)
- [Gibbs](#gibbs)
- [MWG](#mwg)
- [RWM](#rwm)

### t-walk

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner. It has few algorithm specifications, and these are recommended to remain at the default values. This algorithm requires an additional vector of initial values, but is otherwise fully automatic. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. However, proposals are not generated with a multivariate distribution and a proposal covariance matrix. Proposals are generated for only a subset of parameters at each iteration, and for only one of two sets of initial values. See below. |

The t-walk (twalk) algorithm of [\[58\]](#ref58) is a general-purpose
MCMC algorithm that requires no tuning, is scale-invariant, is
technically non-adaptive (but self-adjusting), and can sample from
target distributions with arbitrary scale and correlation structures. A
random subset of one of two vectors is moved around the state-space to
influence one of two chains, per iteration.

twalk has four algorithm specifications: `SIV` is a vector secondary
initial values and the default is `NULL`, `n1` affects the size of the
subset of each set of points to adjust and the default is 4, `at`
affects the traverse move and the default is 6, and `aw` affects the
walk move and the default is 1.5. The vector of secondary initial values
may be left to its default, `NULL`, in which case it is generated with
the GIV function and it is recommended that `PGF` is specified. `SIV`
should be similar to, but unique from, `Initial.Values`. The secondary
initial values are used for a second chain, which is merely used here to
help the first chain, and its results are not reported. The `n1`
specification relates to the number of parameters. For example, if
`n1=4` and a model has \\J=100\\ parameters, then there is a \\p(0.04) =
4 / 100\\ probability that a point is moved that affects each parameter,
though this affects only one of two chains per iteration. Put another
way, there is a 4% chance that each parameter changes each iteration,
and a 50% chance each iteration that the observed chain is selected. The
traverse specification argument, `at`, affects the traverse move, which
helps when some parameters are highly correlated, and the correlation
structure may change through the state-space. The traverse move is
associated with an acceptance rate that decreases as the number of
parameters increases, and is the reason that n1 is used to select a
subset of parameters each iteration. Finally, the walk specification
argument, `aw`, affects the walk move. The authors recommend keeping
these specification arguments in `n1` in \[2,20\], `at`’ in \[2,10\],
and `aw` in \[0.3,2\]. The ‘hop’ and ‘blow’ moves do not have
specifications, but help with multimodality, ensure irreducibility, and
prevent the two chains from collapsing together. The hop move is
centered on the primary chain, and the blow move is centered on the
secondary chain.

The authors have provided the t-walk algorithm in R code as well as
other languages. It is called the “t-walk” for “traverse” or
“thoughtful” walk, as opposed to Random-Walk Metropolis (RWM). Where
adaptive algorithms are designed to adapt to the scale and correlation
structure of target distributions, the t-walk is invariant to this
structure. The step-size and direction continuously “adjust” to the
local structure. The t-walk uses one of four proposal distributions or
‘moves’ per iteration, with the following probabilities: traverse
(p=0.4918), walk (p=0.4918), hop (p=0.0082), and blow (p=0.0082).

Testing in lucifer with the default specifications suggests the t-walk
is very promising, but due to the subset of proposals, it is important
to note that the reported acceptance rate indicates the proportion of
iterations in which moves were accepted, but that only a subset of
parameters changed, and only one chain is selected each iteration.
Therefore, a user who updates a high-dimensional model should find that
parameter values change much less frequently, and this requires more
iterations.

The main advantage of t-walk, like the Hit-And-Run Metropolis (HARM) and
Metropolis-within-Gibbs (MWG) families, over multivariate adaptive
algorithms such as Adaptive-Mixture Metropolis (AMM) and Robust Adaptive
Metropolis (RAM) is that t-walk does not adapt to a proposal covariance
matrix, which can be limiting in random-access memory (RAM) and other
respects in large dimensions, making t-walk suitable for
high-dimensional exploration. Other advantages are that t-walk is
invariant to all but the most extreme correlation structures, does not
need to burn-in before adapting since it technically is non-adaptive
(though it ‘adjusts’ continuously), and continuous adjustment is an
advantage, so `Periodicity` does not need to be adjusted. The advantage
of t-walk over componentwise algorithms such as the MWG family, is that
the model specification does not have to be evaluated a number of times
equal to the number of parameters in each iteration, allowing the t-walk
algorithm to iterate significantly faster in high dimension. The
disadvantage of t-walk, compared to these algorithms, is that more
iterations are required because only a subset of parameters can change
at each iteration (though it still updates twice the number of
parameters per iteration, on average, than the MWG family).

Since twalk is technically non-adaptive, it is suitable as a final
algorithm.

#### Algorithm

1.  **Initialize** two sets of parameters \\\theta_0, \theta'\_0\\
    (primary and secondary chains)
2.  **For** \\t = 1, \ldots, T\\:
    1.  Select one of the two chains at random; select a move type:
        - **Traverse** (\\p = 0.4918\\): for a random subset of
          parameters (size \\\propto n_1/J\\), propose \\\theta^\*\_j =
          \theta'\_j + \beta (\theta'\_j - \theta_j)\\ where \\\beta\\
          depends on `at`
        - **Walk** (\\p = 0.4918\\): propose \\\theta^\*\_j = \theta_j +
          u_j (\theta_j - \theta'\_j)\\ where \\u_j \sim
          \text{Uniform}(-\text{aw}, \text{aw})\\
        - **Hop** (\\p = 0.0082\\): propose from \\\mathcal{N}(\theta_j,
          \|\theta_j - \theta'\_j\|/3)\\ (centered on primary)
        - **Blow** (\\p = 0.0082\\): propose from
          \\\mathcal{N}(\theta'\_j, \|\theta_j - \theta'\_j\|)\\
          (centered on secondary)
    2.  Compute \\\alpha\\ (with appropriate Jacobian for traverse move)
    3.  With probability \\\alpha\\, accept the proposal for the
        selected chain
3.  **Return** samples from the primary chain

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| SIV | (auto-generated) | Secondary initial values; lucifer regenerates them via PGF when SIV is NULL. |
| n1 | 4 | Subset size for the walk move; n1 = 4 recommended by Christen and Fox (2010). |
| at | 6 | Traverse-move parameter; 6 recommended. |
| aw | 1.5 | Walk-move parameter; 1.5 recommended. |

Default `Specs` for twalk. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "twalk", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AMM](#amm)
- [HARM](#harm)
- [MWG](#mwg)
- [RAM](#ram)
- [RWM](#rwm)

### Quadratic Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and the knot range parameter a. Typical rates are 20-50%. |
| Applications | General-purpose ensemble sampler suitable for low to moderate dimensions. Particularly effective on posteriors with curved ridges where linear proposals fail. |
| Difficulty | Easy. Requires setting ensemble size Nc and knot range a. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Quadratic Monte Carlo (QMC) was introduced by [\[59\]](#ref59) as a
generalization of the affine-invariant ensemble sampler (AIES) that
replaces linear stretch moves with second-order Lagrange polynomial
interpolation through three walkers. The target walker \\i\\ and two
randomly selected guide walkers \\j\\ and \\k\\ from the complementary
ensemble define a quadratic arc in parameter space. By drawing a
proposal knot \\t'\\ and evaluating the Lagrange weights, the proposal
follows the curvature of the fitness landscape rather than interpolating
linearly between two points.

QMC has four algorithm specifications: `Nc` is the number of walkers
(minimum 3, recommended \\2J\\ to \\3J\\ for \\J\\ parameters), `a` is
the knot range parameter controlling proposal aggressiveness
(recommended 0.5-1.0 for low dimensions, smaller for high dimensions),
`Z` is an optional \\N_c \times J\\ matrix of initial walker positions,
and `Gaussian` is a logical flag selecting Gaussian (\\t \sim
\mathcal{N}(0, a^2)\\) rather than uniform (\\t \sim \text{Unif}\[-a,
a\]\\) knot sampling. The Jacobian correction \\\|w_i\|^N\\ ensures
detailed balance. Because the Jacobian penalty grows exponentially with
dimension, the knot range `a` must be reduced as \\N\\ increases.

QMC is the foundation of the QMC family and serves as a general-purpose
replacement for AIES on curved posteriors. Compared to AIES, QMC
produces higher acceptance rates on banana-shaped distributions because
the quadratic arc tracks the posterior ridge. The cost is marginally
higher per proposal (three walkers instead of two). Since QMC is not
adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers with positions \\\\z_1,
    \ldots, z\_{N_c}\\\\, each of dimension \\J\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i = 1, \ldots, N_c\\:
        - Select two guide walkers \\j, k \neq i\\ uniformly at random
        - Draw target knot \\t_i \sim P_S(t)\\ and proposal knot \\t'
          \sim P_S(t)\\; fix guide knots \\t_j = -1\\, \\t_k = +1\\
        - Compute Lagrange weights \\w_i, w_j, w_k\\ from the three
          knots
        - Propose \\z^\* = w_i z_i + w_j z_j + w_k z_k\\
        - Compute \\\alpha = \min\\\bigl(1,\\ \|w_i\|^J \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
        - With probability \\\alpha\\, set \\z_i \leftarrow z^\*\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers in the quadratic ensemble. |
| Z | NULL | Auto-initialized. |
| a | 1.0 | Quadratic step amplitude. |
| Gaussian | FALSE | Use the quadratic proposal rather than the Gaussian fallback. |

Default `Specs` for QMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "QMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AIES](#aies)
- [DEMC](#demc)
- [QMC3](#qmc3)
- [DQMC](#dqmc)

### Third-Order Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and the knot range parameter a. Typical rates are 20-50%. |
| Applications | General-purpose ensemble sampler suitable for low to moderate dimensions. The cubic interpolant can match S-shaped or asymmetrically curved manifolds that a quadratic cannot. |
| Difficulty | Easy. Requires setting ensemble size Nc and knot range a. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Third-Order Monte Carlo (QMC3) extends QMC from quadratic to cubic
Lagrange polynomial interpolation by using three guide walkers with
fixed knots at \\t_j = -1\\, \\t_k = 0\\, and \\t_l = +1\\. The
four-point cubic interpolant can follow S-shaped or asymmetrically
curved manifolds that a second-order parabola cannot match. The cost is
an additional guide walker per proposal and slightly more extreme
Lagrange weights for the same knot range, so somewhat smaller `a` may be
needed compared to QMC.

QMC3 has four algorithm specifications: `Nc` is the number of walkers
(minimum 4), `a` is the knot range parameter, `Z` is an optional matrix
of initial walker positions, and `Gaussian` selects Gaussian knot
sampling. The Jacobian correction and acceptance probability are
identical in form to QMC, using \\\|w_i\|^J\\ where \\w_i\\ is the
Lagrange weight of the target walker.

Since QMC3 is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 4\\ walkers with positions \\\\z_1,
    \ldots, z\_{N_c}\\\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select three guide walkers \\j, k, l \neq i\\ at random; fix
          knots \\t_j = -1\\, \\t_k = 0\\, \\t_l = +1\\
        - Draw \\t_i \sim P_S(t)\\ and proposal knot \\t' \sim P_S(t)\\
        - Compute degree-3 Lagrange weights \\(w_i, w_j, w_k, w_l)\\
        - Propose \\z^\* = w_i z_i + w_j z_j + w_k z_k + w_l z_l\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^J \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default                  | Justification           |
|:--------------|:-------------------------|:------------------------|
| Nc            | max(2\*LIV+1, 10) (= 10) | Ten walkers.            |
| Z             | NULL                     | Auto-initialized.       |
| a             | 1.0                      | Cubic step amplitude.   |
| Gaussian      | FALSE                    | Use the cubic proposal. |

Default `Specs` for QMC3. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "QMC3", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [QMC](#qmc)
- [QMCN](#qmcn)
- [AIES](#aies)

### Nth-Order Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension, knot range a, and polynomial order. Higher orders produce more extreme weights, typically requiring smaller a. |
| Applications | General-purpose ensemble sampler that can match complex posterior geometry. Most useful when QMC or QMC3 underperform due to insufficient polynomial flexibility. |
| Difficulty | Moderate. Requires setting ensemble size Nc, knot range a, and polynomial order nOrder. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Nth-Order Monte Carlo (QMCN) generalizes QMC to arbitrary polynomial
order controlled by the `nOrder` parameter. The method selects `nOrder`
guide walkers placed at equispaced knots on \\\[-1, +1\]\\ and
constructs the degree-`nOrder` Lagrange interpolant. This flexibility
allows matching complex posterior geometry, but diminishing returns set
in beyond order 4 or 5 because the Lagrange weights oscillate
increasingly at high polynomial orders, analogous to Runge’s phenomenon
in polynomial interpolation.

QMCN has five algorithm specifications: `Nc` is the number of walkers
(minimum `nOrder + 1`), `a` is the knot range parameter, `nOrder` is the
polynomial order (default 4), `Z` is an optional matrix of initial
walker positions, and `Gaussian` selects Gaussian knot sampling.
Higher-order interpolation requires more guide walkers per proposal and
generally more conservative knot ranges to maintain acceptable
acceptance rates.

Since QMCN is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq n\_{\text{order}} + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n\_{\text{order}}\\ guide walkers at random; place at
          equispaced knots on \\\[-1, +1\]\\
        - Draw \\t_i \sim P_S(t)\\ and proposal knot \\t' \sim P_S(t)\\
        - Compute degree-\\n\_{\text{order}}\\ Lagrange weights
        - Propose \\z^\* = \sum_m w_m z_m\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^J \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers. |
| Z | NULL | Auto-initialized. |
| a | 1.0 | Polynomial step amplitude. |
| Gaussian | FALSE | Use the polynomial proposal. |
| nOrder | 4 | Quartic polynomial fit; default in Militzer (2025). |

Default `Specs` for QMCN. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "QMCN", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [QMC](#qmc)
- [QMC3](#qmc3)
- [AIES](#aies)

### Directed Quadratic Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension, knot range a, and the directed width b. The energy-directed step can improve acceptance on unimodal posteriors. |
| Applications | General-purpose ensemble sampler with energy-guided proposals. Best suited to unimodal posteriors with strong curvature. |
| Difficulty | Moderate. Requires setting ensemble size Nc, knot range a, and directed width b. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Directed Quadratic Monte Carlo (DQMC) augments the standard QMC proposal
with energy guidance [\[60\]](#ref60). After computing the three QMC
knots and evaluating the negative log-posterior \\E = -\log \pi\\ at the
corresponding positions, DQMC fits a parabola \\E(t) = At^2 + Bt + C\\
through the three energy values. When \\A \> 0\\ (the parabola opens
upward), the minimum is at \\t\_{\min} = -B/(2A)\\, and the proposal
knot is drawn from \\t' \sim \mathcal{N}(t\_{\min}, b^2)\\ rather than
from the undirected distribution. The parameter `b` controls how
aggressively the proposal targets the energy minimum: small `b`
concentrates proposals near \\t\_{\min}\\, while large `b` diffuses them
and recovers undirected QMC behavior. When \\A \leq 0\\, the method
falls back to standard undirected QMC.

Because the forward and reverse proposal densities differ, DQMC requires
an asymmetric Metropolis-Hastings correction. The acceptance probability
includes a ratio of the undirected and directed proposal densities
evaluated at both the forward and reverse knots.

DQMC has five algorithm specifications: `Nc` is the number of walkers
(minimum 3), `a` is the knot range parameter, `b` is the directed
proposal width (default 0.5), `Z` is an optional matrix of initial
walker positions, and `Gaussian` selects Gaussian knot sampling. Since
DQMC is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two guide walkers \\j, k\\; draw target knot \\t_i\\;
          fix \\t_j = -1\\, \\t_k = +1\\
        - Evaluate energies \\E_i, E_j, E_k\\ at the three walker
          positions
        - Fit parabola \\E(t) = At^2 + Bt + C\\; if \\A \> 0\\, set
          \\t\_{\min} = -B/(2A)\\
        - If directed: draw \\t' \sim \mathcal{N}(t\_{\min}, b^2)\\;
          else: draw \\t' \sim P_S(t)\\
        - Compute Lagrange weights and propose \\z^\* = w_i z_i + w_j
          z_j + w_k z_k\\
        - Accept with asymmetric Metropolis-Hastings ratio including
          Jacobian \\\|w_i\|^J\\ and proposal density corrections
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default                  | Justification                     |
|:--------------|:-------------------------|:----------------------------------|
| Nc            | max(2\*LIV+1, 10) (= 10) | Ten walkers.                      |
| Z             | NULL                     | Auto-initialized.                 |
| a             | 1.0                      | Quadratic amplitude.              |
| Gaussian      | FALSE                    | Use directed quadratic proposals. |
| b             | 1.0                      | Drift bias parameter.             |

Default `Specs` for DQMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DQMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [QMC](#qmc)
- [AIES](#aies)

### Quadratic Simplex Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and knot range a. Simplex averaging can improve proposal quality when the ensemble is large. |
| Applications | General-purpose ensemble sampler. The averaging suppresses individual walker noise and can improve proposals for noisy or multimodal ensembles. |
| Difficulty | Moderate. Requires setting ensemble size Nc, knot range a, and the number of guide points nGuidePoints. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Quadratic Simplex Monte Carlo (SQMC) replaces the individual guide
walkers in QMC with weighted averages of groups of walkers
[\[60\]](#ref60). Two groups of `nGuidePoints` walkers are selected from
the complementary ensemble, and each group is collapsed to a randomly
weighted centroid using Dirichlet-like weights controlled by
`weightRange`. The standard QMC interpolation then uses these averaged
positions instead of individual walkers. This averaging suppresses
individual walker noise and can improve proposal quality when the
ensemble is large, at the cost of requiring more walkers: \\N_c \geq 2
\cdot n_G + 1\\ where \\n_G\\ is `nGuidePoints`.

SQMC has six algorithm specifications: `Nc` is the number of walkers,
`a` is the knot range parameter, `nGuidePoints` is the number of walkers
per guide group (default 3), `weightRange` controls the Dirichlet weight
variability (default 1.0), `Z` is an optional matrix of initial walker
positions, and `Gaussian` selects Gaussian knot sampling. The Jacobian
correction and acceptance probability follow the same form as QMC.

Since SQMC is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 2 n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two groups of \\n_G\\ guide walkers from the
          complementary ensemble
        - For each group, draw Dirichlet-like weights and compute
          weighted centroid \\\bar{z}\_j\\, \\\bar{z}\_k\\
        - Perform standard QMC interpolation using \\z_i\\,
          \\\bar{z}\_j\\, \\\bar{z}\_k\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ \|w_i\|^J \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*3 + LIV, 10) (= 10) | Ten guide-supported walkers. |
| Z | NULL | Auto-initialized. |
| a | 1.0 | Quadratic amplitude. |
| Gaussian | FALSE | Use the quadratic simplex proposal. |
| nGuidePoints | 3 | Three simplex guide points provide adequate coverage in low dimensions. |
| weightRange | 1.0 | Equal-weight guide points. |

Default `Specs` for SQMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "SQMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [QMC](#qmc)
- [SAMC](#samc)
- [AIES](#aies)

### Modified Affine Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and the scale parameter a (must be \> 1). Typical rates are 20-50%, similar to AIES. |
| Applications | General-purpose ensemble sampler. The two-guide stretch can capture correlations between guide positions that single-guide AIES misses. |
| Difficulty | Easy. Requires setting ensemble size Nc and scale parameter a (default 2.0). |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Modified Affine Monte Carlo (MAMC) uses two guide walkers to define both
the stretch direction and the reference point, unlike standard AIES
which stretches from a single companion [\[60\]](#ref60). The proposal
is \\z' = z_j + Z(z_i - z_k)\\ where \\Z\\ is drawn from the stretch
distribution \\g(z) \propto 1/\sqrt{z}\\ on \\\[1/a, a\]\\ and walkers
\\j\\ and \\k\\ are selected at random from the complementary ensemble.
By using two guides to define the direction, MAMC can capture
correlations between guide positions that AIES misses. The acceptance
correction is \\Z^{J-1}\\ as in standard AIES.

MAMC has four algorithm specifications: `Nc` is the number of walkers
(minimum 3), `a` is the scale parameter (must be \> 1, default 2.0), `Z`
is an optional matrix of initial walker positions, and `Gaussian` is
ignored (MAMC always uses the stretch distribution). Since MAMC is not
adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq 3\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select two guide walkers \\j, k \neq i\\ at random
        - Draw stretch factor \\Z \sim g(z)\\ on \\\[1/a,\\ a\]\\
        - Propose \\z^\* = z_j + Z (z_i - z_k)\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ Z^{J-1} \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers. |
| Z | NULL | Auto-initialized. |
| a | 2.0 | Stretch parameter \> 1 (mandatory by validator). |

Default `Specs` for MAMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MAMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AIES](#aies)
- [SAMC](#samc)
- [QMC](#qmc)

### Affine Simplex Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and the scale parameter a (must be \> 1). The centroid reference provides more stable proposals than single-walker AIES. |
| Applications | General-purpose ensemble sampler. The centroid-based stretch reduces proposal variance compared to AIES. |
| Difficulty | Easy. Requires setting ensemble size Nc, scale parameter a, and nGuidePoints. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Affine Simplex Monte Carlo (SAMC) combines the AIES stretch move with
simplex averaging [\[60\]](#ref60): instead of stretching from a single
companion, it stretches from the weighted centroid of `nGuidePoints`
walkers. The proposal is \\z' = \bar{z} + Z(z_i - \bar{z})\\ where
\\\bar{z}\\ is the guide centroid and \\Z\\ follows the stretch
distribution \\g(z) \propto 1/\sqrt{z}\\ on \\\[1/a, a\]\\. The centroid
provides a more stable reference point than a single walker, which can
reduce proposal variance. The acceptance correction is \\Z^{J-1}\\.

SAMC has five algorithm specifications: `Nc` is the number of walkers
(minimum \\n_G + 1\\), `a` is the scale parameter (must be \> 1, default
2.0), `nGuidePoints` is the number of walkers used to compute the
centroid (default 3), `Z` is an optional matrix of initial walker
positions, and `Gaussian` is ignored (SAMC always uses the stretch
distribution). Since SAMC is not adaptive, it is suitable as a final
algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n_G\\ guide walkers from the complementary ensemble;
          compute centroid \\\bar{z}\\
        - Draw stretch factor \\Z \sim g(z)\\ on \\\[1/a,\\ a\]\\
        - Propose \\z^\* = \bar{z} + Z(z_i - \bar{z})\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ Z^{J-1} \cdot p(z^\*
          \mid y) / p(z_i \mid y)\bigr)\\
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default                 | Justification         |
|:--------------|:------------------------|:----------------------|
| Nc            | max(3 + LIV, 10) (= 10) | Ten walkers.          |
| Z             | NULL                    | Auto-initialized.     |
| a             | 2.0                     | Stretch \> 1.         |
| nGuidePoints  | 3                       | Three simplex guides. |
| weightRange   | 1.0                     | Equal weights.        |

Default `Specs` for SAMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "SAMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AIES](#aies)
- [MAMC](#mamc)
- [SQMC](#sqmc)

### Walk Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on dimension and scale parameter a. No Jacobian correction is needed because the proposal is symmetric. |
| Applications | General-purpose ensemble sampler. Serves as a robust baseline when the other QMC methods encounter Jacobian penalties that are too severe. |
| Difficulty | Easy. Requires setting ensemble size Nc, scale parameter a, and nGuidePoints. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate by Walker. |

Walk Monte Carlo (WMC) constructs a symmetric random walk perturbation
from the ensemble covariance structure [\[60\]](#ref60). The proposal is
\\z' = z_i + \sum\_{g=1}^{n_G} a \\ Z_g \\ (z_g - \bar{z})\\ where \\Z_g
\sim \mathcal{N}(0,1)\\, \\z_g\\ are guide walker positions, and
\\\bar{z}\\ is their mean. Because the Gaussian perturbation is
symmetric around \\z_i\\, no Jacobian correction is needed; the
acceptance reduces to the standard Metropolis ratio \\\min\[1,
\pi(z')/\pi(z_i)\]\\. WMC is the simplest method in the QMC family and
serves as a robust baseline when the other methods’ Jacobian penalties
become problematic in high dimensions.

WMC has five algorithm specifications: `Nc` is the number of walkers
(minimum \\n_G + 1\\), `a` is the scale parameter controlling
perturbation magnitude (default 1.0), `nGuidePoints` is the number of
guide walkers (default 3), `Z` is an optional matrix of initial walker
positions, and `Gaussian` is ignored. Since WMC is not adaptive, it is
suitable as a final algorithm.

#### Algorithm

1.  **Initialize** \\N_c \geq n_G + 1\\ walkers
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\i\\:
        - Select \\n_G\\ guide walkers; compute their mean \\\bar{z}\\
        - Draw \\Z_g \sim \mathcal{N}(0,1)\\ for \\g = 1, \ldots, n_G\\
        - Propose \\z^\* = z_i + a \sum\_{g=1}^{n_G} Z_g (z_g -
          \bar{z})\\
        - Accept with \\\alpha = \min\\\bigl(1,\\ p(z^\* \mid y) / p(z_i
          \mid y)\bigr)\\ (no Jacobian correction)
3.  **Return** samples from the main walker

#### Default specifications

| Specification | Default                 | Justification     |
|:--------------|:------------------------|:------------------|
| Nc            | max(3 + LIV, 10) (= 10) | Ten walkers.      |
| Z             | NULL                    | Auto-initialized. |
| a             | 1.0                     | Walk amplitude.   |
| nGuidePoints  | 3                       | Three guides.     |

Default `Specs` for WMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "WMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AIES](#aies)
- [QMC](#qmc)
- [MAMC](#mamc)

### Ensemble Slice Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate is effectively 100% because slice sampling always accepts. |
| Applications | This is a widely applicable, gradient-free algorithm that inherits the affine invariance of ensemble methods. It is particularly useful for correlated posteriors where gradient-based methods are unavailable or expensive. |
| Difficulty | This algorithm is easy for a beginner. It has three algorithm specifications, and the main tuning parameter (slice width mu) adapts automatically. |
| Final Algorithm? | Yes. |
| Proposal | Ensemble. |

The Ensemble Slice Sampler (ESS, marketed as “zeus”) was introduced by
[\[61\]](#ref61) as a gradient-free, affine-invariant ensemble MCMC
method that uses slice sampling instead of the stretch move (used by
AIES) or differential evolution (used by DEMC). An ensemble of \\N_c\\
walkers is maintained, split into two complementary subsets \\S_a\\ and
\\S_b\\. For each walker in \\S_a\\, two walkers from \\S_b\\ define a
direction vector \\d = x_k - x_j\\; univariate slice sampling is then
performed along the line \\x_i + t \cdot d\\ for \\t \in \mathbb{R}\\.
The slice sampling uses the doubling-out procedure of [\[15\]](#ref15)
to find the slice interval, followed by shrinkage to locate the accepted
position. After all walkers in \\S_a\\ are updated, the roles swap and
walkers in \\S_b\\ are updated using \\S_a\\.

ESS has three algorithm specifications: `Nc` is the number of walkers
(default \\\max(2d+1, 10)\\), `Z` is an optional \\N_c \times d\\ matrix
of initial walker positions (generated via GIV if NULL), and `mu` is the
initial slice width (default 1.0). Only the first walker’s samples are
thinned into the output.

The advantage of ESS over AIES is that slice sampling eliminates the
need for step size tuning entirely: the doubling procedure automatically
finds the appropriate scale. The advantage over DEMC is that slice
sampling always accepts, so there is no wasted computation from rejected
proposals. The disadvantage, shared with all ensemble methods, is the
cost of maintaining \\N_c\\ walkers (each requiring a model evaluation),
and the restriction that \\N_c \geq 2d + 1\\ for the direction vectors
to span the parameter space.

#### Algorithm

1.  **Initialize** \\N_c\\ walkers \\\\x_1, \ldots, x\_{N_c}\\\\; split
    into \\S_a\\ (even indices) and \\S_b\\ (odd indices)
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each walker \\x_i \in S_a\\:
        - Pick two walkers \\x_j, x_k\\ from \\S_b\\; set direction \\d
          = x_k - x_j\\
        - Draw slice level \\\log y = \log \pi(x_i) - \text{Exp}(1)\\
        - **Doubling**: initialize interval \\\[L, R\] = \[-\mu,
          \mu\]\\; while \\\log \pi(x_i + L \cdot d) \geq \log y\\ or
          \\\log \pi(x_i + R \cdot d) \geq \log y\\, double the interval
        - **Shrinkage**: sample \\t^\* \sim \text{Uniform}(L, R)\\; if
          \\\log \pi(x_i + t^\* d) \geq \log y\\, accept \\x_i
          \leftarrow x_i + t^\* d\\; otherwise shrink interval toward 0
    2.  **Swap roles**: repeat for each walker in \\S_b\\ using \\S_a\\
3.  **Return** thinned samples from the first walker

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers in the ensemble slice. |
| Z | NULL | Auto-initialized. |
| mu | 1.0 | Initial slice mode used by the doubling procedure. |

Default `Specs` for Zeus. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Zeus", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [AIES](#aies)
- [DEMC](#demc)
- [Slice](#slice)

### Differential Evolution Adaptive Metropolis

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate depends on the number of chain pairs and the crossover probability. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, gradient-free ensemble algorithm. It is particularly effective for multimodal distributions and high-dimensional problems where standard DE methods struggle. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has five algorithm specifications, though the defaults work well for most problems. |
| Final Algorithm? | User Discretion. DREAM may be used as a final algorithm, though it includes periodic snooker updates that introduce non-stationarity. |
| Proposal | Ensemble (multivariate). |

Differential Evolution Adaptive Metropolis (DREAM) was introduced by
[\[62\]](#ref62) as an enhancement of the Differential Evolution Markov
Chain (DEMC) algorithm. DREAM maintains \\N_c\\ chains and generates
proposals via differential evolution, but improves on DEMC in three
ways: (1) multiple chain pairs are used for the differential evolution
proposal instead of a single pair, reducing the variance of the proposal
direction; (2) a crossover operation randomly selects a subset of
dimensions to update (with probability CR per dimension), which is
beneficial in high dimensions where most dimensions may be nearly
independent; and (3) periodic snooker updates (every 5th iteration,
using \\\gamma = 1\\ instead of the standard scaling) improve
exploration of the tails.

DREAM has five algorithm specifications: `Nc` is the number of chains
(default \\\max(2d+1, 10)\\), `Z` is an optional initial position
matrix, `n_pairs` is the number of chain pairs for the DE proposal
(default 3), `CR` is the crossover probability (default 0.9), and
`b_star` is the scale of additive noise for ergodicity (default
\\10^{-6}\\). The proposal for chain \\i\\ is \\\theta^\* = \theta_i +
\gamma \sum\_{p=1}^{\delta} (\theta\_{r\_{1p}} - \theta\_{r\_{2p}}) +
\varepsilon\\, where \\\gamma = 2.38 / \sqrt{2 \delta d\_{\text{CR}}}\\,
\\\delta\\ is the number of pairs, \\d\_{\text{CR}}\\ is the number of
dimensions being updated, and \\\varepsilon \sim \mathcal{N}(0, b^\*
I)\\.

The advantage of DREAM over DEMC is faster convergence in high
dimensions due to the multiple pairs and crossover mechanism. The
advantage over AIES and ESS is that DREAM’s DE proposals can be more
effective for multimodal distributions because multiple chain pairs
provide better directional information. The disadvantage is the
computational overhead of maintaining \\N_c\\ chains and the additional
complexity of tuning the crossover probability.

#### Algorithm

1.  **Initialize** \\N_c\\ chains \\\\x_1, \ldots, x\_{N_c}\\\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each chain \\i\\:
        - Sample \\\delta\\ pairs of distinct chains \\(r\_{1p},
          r\_{2p})\\ with \\r \neq i\\
        - Generate crossover mask: include dimension \\j\\ with
          probability CR (at least one dimension)
        - Compute \\\gamma = 2.38/\sqrt{2\delta \cdot d\_{\text{CR}}}\\;
          every 5th iteration, set \\\gamma = 1\\ (snooker)
        - Propose \\\theta_i^\* = \theta_i + \gamma \sum_p
          (\theta\_{r\_{1p}} - \theta\_{r\_{2p}}) + \mathcal{N}(0, b^\*
          I)\\, applied only to selected dimensions
        - Accept/reject via standard MH ratio
3.  **Return** thinned samples from the first chain

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Nc | max(2\*LIV+1, 10) (= 10) | Ten walkers. |
| Z | NULL | Auto-initialized. |
| n_pairs | 3 | Three difference pairs per update; the value recommended by Vrugt (2016). |
| CR | 0.9 | Crossover probability. |
| b_star | 1e-6 | Snooker noise variance. |

Default `Specs` for DREAM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "DREAM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [DEMC](#demc)
- [AIES](#aies)
- [Zeus](#zeus)

## Hit-and-run and directional samplers propose along random directions

Hit-and-run methods choose a random direction in parameter space and
then propose a new point along that direction, effectively reducing the
multivariate sampling problem to a series of univariate problems along
random lines. The componentwise variant (CHARM) restricts directions to
coordinate axes, while the full hit-and-run (HARM) samples directions
uniformly on the unit sphere. The refractive sampler extends this idea
by using gradient information to bend the proposal direction at density
boundaries, analogously to Snell’s law in optics, achieving better
penetration of narrow regions.

### Hit-And-Run Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner. There are no required algorithm specifications (though two are optional) and tuning is unnecessary. |
| Final Algorithm? | Yes, if not adaptive. |
| Proposal | Blockwise or Multivariate. Proposals are multivariate only in the sense that proposals for multiple parameters are generated at once. However, proposals are not generated with a multivariate distribution and a proposal covariance matrix. |

The Hit-And-Run algorithm is a variation of Random-Walk Metropolis (RWM)
that has been around as long as Gibbs sampling (Gibbs). Hit-And-Run
randomly samples a direction on the unit sphere as in [\[65\]](#ref65),
and a proposal is made for each parameter in its randomly-selected
direction for a uniformly-distributed distance. This version of
Hit-And-Run, called Hit-And-Run Metropolis (HARM), includes multivariate
proposals and a Metropolis step for rejection. Introduced by
[\[16\]](#ref16), along with the original Gibbs sampler, and popularized
by [\[66\]](#ref66), Hit-And-Run was given its name later due to its
ability to run across the state-space and arrive at a distant
“hit-point”. It is related to other algorithms with interesting names,
such as Hide-And-Seek and Shake-And-Bake.

HARM has two optional algorithm specifications: `alpha.star` is the
target acceptance rate, and `B` is a list of blocked parameters. When
`alpha.star=NULL`, non-adaptive HARM is used. Otherwise, the
Robbins-Monro stochastic approximation of [\[64\]](#ref64) is applied to
the proposal distance. `alpha.star=0.234` is recommended. For blockwise
sampling, each component of the list is a block and consists of a vector
indicating the position of the parameters per block.

HARM is the fastest MCMC algorithm per iteration in the lucifer package,
because it is very simple. For example, HARM does not use a proposal
covariance matrix, and there are no tuning parameters, with one optional
exception discussed below. The proposal is multivariate in the sense
that all parameters are proposed at once, though from univariate draws.
HARM often mixes better than the Gibbs sampler [\[65\]](#ref65),
especially with correlated parameters (Chen and Schmeiser 1992).

Adaptive HAR (without the Metropolis step) with a multivariate proposal
is available in the LaplaceApproximation function.

The HARM algorithm is able to traverse complex spaces with bounded sets
in one iteration. As such, HARM may work well with multimodal posteriors
due to potentially good mode-switching behavior. However, HARM may have
difficulty sampling regions of high probability that are spike-shaped or
near constraints, and this difficulty is likely to be more problematic
in high dimensions. When HARM is non-adaptive, it may be used as a final
algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Sample a random direction \\d\\ uniformly on the unit sphere in
        \\\mathbb{R}^J\\
    2.  Draw distance \\\lambda \sim \text{Uniform}(-1, 1)\\ (or
        Robbins-Monro-scaled if adaptive)
    3.  Propose \\\theta^\* = \theta\_{t-1} + \lambda \cdot d\\
        (multivariate proposal along random direction)
    4.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
    6.  If adaptive, update step-size via Robbins-Monro toward
        `alpha.star`
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| alpha.star | NA | When NA, the chain runs without acceptance-rate coercion. Set to 0.234 (or 0.44) to coerce. |
| B | NULL | Single block. |

Default `Specs` for HARM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "HARM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [CHARM](#charm)
- [Gibbs](#gibbs)

### Componentwise Hit-And-Run Metropolis

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is 44%, and is based on the univariate normality of each marginal posterior distribution. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is easy for a beginner. There are no algorithm specifications (although one is optional), and tuning is unnecessary. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Hit-And-Run algorithm is a variation of RWM that has been around as
long as Gibbs sampling. Hit-And-Run randomly samples a direction on the
unit sphere, and a proposal is made for each parameter in its
randomly-selected direction for a uniformly-distributed distance
[\[65\]](#ref65). This version of Hit-And-Run, called Componentwise
Hit-And-Run Metropolis (CHARM), includes componentwise proposals and a
Metropolis step for rejection. Introduced by [\[16\]](#ref16) along with
Gibbs sampling, and popularized by [\[66\]](#ref66), Hit-And-Run was
given its name later due to its ability to run across the state-space
and arrive at a distant “hit-point”. It is related to other algorithms
with interesting names, such as Hide-And-Seek and Shake-And-Bake.

CHARM has one optional algorithm specification: `alpha.star`. When
`Specs=NULL`, CHARM is non-adaptive. Otherwise, `alpha.star` is the
target acceptance rate, and is recommended to be 0.44. When a target
acceptance rate is declared, CHARM applies the Robbins-Monro stochastic
approximation of [\[64\]](#ref64) to the proposal distance to attain the
target acceptance rate.

As a componentwise algorithm, the model is evaluated after a proposal is
made for each parameter, which results in an algorithm that takes more
time per iteration. As with the Metropolis-within-Gibbs (MWG) family,
the time to complete each iteration grows with the number of parameters.
Compared to other algorithms with multivariate proposals, a disadvantage
is the time to complete each iteration increases as a function of
parameters and model complexity. For example, in a 100-parameter model,
CHARM completes its first iteration as HARM completes its 100^(th).

CHARM enjoys many of the advantages of HARM, such as having no tuning
parameters (unless the adaptive form is used), traversing complex spaces
with bounded sets in one iteration, not being adaptive (unless specified
as adaptive), handling high correlations well, and having the potential
to work well with multimodal distributions. When non-adaptive, CHARM may
be used as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **For** each parameter \\j\\ in random order:
        - Sample a random direction \\d_j\\ on the unit sphere
        - Draw distance \\\lambda \sim \text{Uniform}(-1, 1)\\ (or
          Robbins-Monro-scaled if adaptive)
        - Propose \\\theta^\*\_j = \theta\_{t,j} + \lambda \cdot d_j\\
        - Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
          p(\theta_t \mid y)\bigr)\\
        - With probability \\\alpha\\, accept \\\theta\_{t+1,j} =
          \theta^\*\_j\\
    2.  If adaptive, update step-size via Robbins-Monro toward target
        acceptance rate `alpha.star`
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification                               |
|:--------------|:--------|:--------------------------------------------|
| alpha.star    | NA      | Componentwise hit-and-run without coercion. |

Default `Specs` for CHARM. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "CHARM", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Gibbs](#gibbs)
- [HARM](#harm)
- [MWG](#mwg)
- [RWM](#rwm)

### Refractive Sampler

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is between 0.6-0.7. |
| Applications | This is a widely applicable, general-purpose algorithm. |
| Difficulty | This algorithm is moderately difficult because it requires tuning of the number of steps, as well as the step size. However, the step size may be tuned automatically with the adaptive version, leaving only the number of steps to be tuned by the user. |
| Final Algorithm? | Yes, if non-adaptive, or no if adaptive. |
| Proposal | Multivariate. |

Refractive sampling was introduced by [\[67\]](#ref67) as an alternative
to Hamiltonian Monte Carlo (HMC) and the No-U-Turn Sampler (NUTS). While
the refractive sampler, HMC, and NUTS all use partial derivatives, the
refractive sampler uses partial derivatives only for direction, not
magnitude. The partial derivatives always point toward the side with the
higher index of refraction. Like HMC and unlike NUTS, the refractive
sampler requires tuning.

The Refractive algorithm has four algorithm specifications: `Adaptive`,
the number `m` of steps to take per iteration, step size `w`, and `r` is
the ratio between indices of refraction. The Adaptive argument does not
appear in [\[67\]](#ref67). When Adaptive is less than the number of
iterations, an optional Robbins-Monro stochastic approximation of
[\[64\]](#ref64) is applied to step size `w`. All specifications must be
scalar, and it is recommended that `r=1.3`. The number of steps `m` is
similar to the number of leapfrog steps `L` in HMC, and step size `w` is
similar to `epsilon`. Typically, as the number of parameters increases,
it is often better for m to be larger and `w` smaller.

Compared to other MCMC algorithms, a higher posterior density is often
found. An advantage over NUTS is that iterative speed is consistent. An
advantage over HMC is that Refractive is easier to tune. When `Adaptive`
is greater than the number of iterations, the Refractive sampler is not
adaptive and is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\w\\, number of
    steps \\m\\, refraction ratio \\r\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute gradient \\g = \nabla \log p(\theta\_{t-1} \mid y)\\
        (direction only, not magnitude)
    2.  Initialize velocity \\v\\ from random direction; set \\\theta^\*
        = \theta\_{t-1}\\
    3.  **For** \\s = 1, \ldots, m\\:
        - Move: \\\theta^\* \leftarrow \theta^\* + w \cdot v\\
        - Compute new gradient \\g^\*\\; if the density boundary is
          crossed, **refract** \\v\\ using Snell’s law with ratio \\r\\
    4.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    5.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\
    6.  If adaptive, update \\w\\ via Robbins-Monro toward target
        acceptance
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| Adaptive | 1 | Adaptation begins at iteration 1; refractive samplers are robust to early adaptation. |
| m | 2 | Two refractive bounces per iteration. |
| w | 0.1 | Step-size for the bounce. |
| r | 1.3 | Refractive ratio between dense and sparse media. |

Default `Specs` for Refractive. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Refractive", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [NUTS](#nuts)
- [RSS](#rss)

## Methods for multimodal posteriors and transdimensional targets

Standard MCMC methods can become trapped in a single mode of a
multimodal posterior, producing biased inference. Metropolis-coupled
MCMC (MC3, also called parallel tempering) runs multiple chains at
different temperatures, with hot chains exploring freely across modes
and cold chains sampling accurately, periodically exchanging states to
transfer modal information. Reversible-jump MCMC (RJ) extends the
framework to variable-dimension problems where the number of parameters
itself is unknown, enabling model comparison and model averaging within
a single MCMC run.

### Metropolis-Coupled Markov Chain Monte Carlo

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that specializes in representing multimodal distributions. |
| Difficulty | This algorithm is moderately difficult for a beginner, because the proposal covariance must be tuned, the swap acceptance rate should be tuned, and the number of CPUs selected. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

The Metropolis-Coupled Markov Chain Monte Carlo (MCMCMC) algorithm –
also referred to as MC3, (MC)3, Metropolis-Coupled MCMC, Parallel
Tempering, or Replica Exchange – was introduced by [\[71\]](#ref71) for
multimodal distributions. The name “parallel tempering” was introduced
in [\[69\]](#ref69).

MCMCMC consists of parallel chains that use different temperatures and
combine to form a mixture distribution, where each chain is a mixture
component. The chains or mixture components are updated in parallel,
each with a Metropolis step that is adjusted by temperature. After this
first set of Metropolis steps, two parallel chains are selected at
random, and another adjusted Metropolis step is used to accept or reject
a swap between chains. The act of swapping is referred to as coupling.
Each swap may be a jump between modes. MCMCMC has both within-chain and
between-chain proposals each iteration, but only the coolest chain is
retained.

MCMCMC has four algorithm specifications: `lambda` is either a
positive-only scalar that controls the temperature spacing or a vector
of temperature spacings, `CPUs` is the number of central processing
units in the computer, `Packages` is a vector of any package names that
are required, and `Dyn.libs` are any dynamic libraries that are
required. A larger scalar value of λ results in greater differences in
temperature between chains, and often a lower swap acceptance rate (see
below). Given a scalar λ and m = 1,…,M CPUs, the current or proposal
distribution is exponentiated to this power: 1/\[1 + λ(m - 1)\]. The
number of CPUs must be at least two; as programmed, MCMCMC will not
function on a single-core computer.

MCMCMC is called with
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
using `Chains > 1`. Each iteration,
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
communicates with all CPUs, collects the latest for all the chains, and
then the master process calculates the Metropolis steps for each chain,
as well as the swap.

As with Random-Walk Metropolis (RWM), the proposal covariance matrix
must be tuned. When nothing is known about this matrix, as is often the
case, it is suggested to begin with an identity matrix with small-scale
diagonal values, such as 1e-3. After a short run that hopefully has some
acceptances, another short run can begin with the observed historical
covariance matrix. Eventually, the proposal covariance matrix is usually
satisfactory.

[\[68\]](#ref68) demonstrate that MCMCMC performance improves as the
acceptance rate of proposed swaps approaches 0.234. In lucifer, the swap
acceptance rate is printed to the console at the end of each update. The
swap acceptance rate is affected by the temperature spacing between
parallel chains. If the default temperature spacing that results from a
scalar λ is unacceptable, then a different scalar value may be
attempted, or the user may enter their own temperature spacing directly
as a vector for λ.

Coupling induces dependence among the chains, and the chains are no
longer Markovian. The whole stochastic process of m chains together does
form a Markov chain.

Along with the J-walking algorithm of [\[70\]](#ref70), MCMCMC was one
of the first extensions of Metropolis-Hastings for multimodal
distributions. Many additional MCMC algorithms, such as serial tempering
or simulated tempering, were influenced by or variations of MCMCMC.

The advantage of MCMCMC over RWM is that MCMCMC is better able to
approximate multimodal distributions, and that successful coupling
(swaps) improves mixing. A disadvantage is that most of the information
in the warmer chains is lost, speed per iteration is slower than RWM due
to communication with CPUs, and distributions with many modes require at
least as many CPUs. Since MCMCMC is not adaptive, it is suitable as a
final algorithm.

#### Algorithm

1.  **Initialize** \\M\\ chains at temperatures \\T_1 = 1 \< T_2 \<
    \cdots \< T_M\\, where \\T_m = 1 + \lambda(m-1)\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Within-chain updates**: for each chain \\m\\, propose
        \\\theta^{(m)\*} \sim \mathcal{N}(\theta^{(m)}\_{t-1}, \Sigma)\\
        and accept with \\\alpha_m = \min\\\bigl(1,\\ \[p(\theta^{(m)\*}
        \mid y) / p(\theta^{(m)}\_{t-1} \mid y)\]^{1/T_m}\bigr)\\
    2.  **Swap proposal**: select two chains \\i, j\\ at random; propose
        swapping their states
    3.  Compute swap acceptance \\\alpha\_{\text{swap}} =
        \min\\\bigl(1,\\ \[p(\theta^{(j)} \mid y) / p(\theta^{(i)} \mid
        y)\]^{1/T_i - 1/T_j}\bigr)\\
    4.  With probability \\\alpha\_{\text{swap}}\\, exchange
        \\\theta^{(i)} \leftrightarrow \theta^{(j)}\\
3.  **Return** samples from the coldest chain (\\T_1 = 1\\)

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| lambda | 1 | Geometric temperature spacing exponent; 1 corresponds to evenly spaced inverse temperatures. |
| CPUs | 2 | MCMCMC is intrinsically parallel: at least one cold and one hot chain are needed. |
| Packages | NULL | No additional packages. |
| Dyn.libs | NULL | No dynamic libraries. |

Default `Specs` for MCMCMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MCMCMC", Specs = list(lambda = 1, CPUs = 2L, Packages = NULL, Dyn.libs = NULL),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [RWM](#rwm)

### Reversible-Jump

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is unknown (not for reversible-jump in general, but for the CHARM algorithm within this version of reversible-jump). As such, the observed acceptance rate may be suitable in the interval \[10%,90%\]. |
| Applications | This version of reversible-jump is intended for variable selection and Bayesian Model Averaging (BMA). |
| Difficulty | This algorithm is difficult for a beginner. There are several algorithm specifications. Convergence will be difficult to assess. |
| Final Algorithm? | Yes. |
| Proposal | Componentwise. |

Reversible-jump Markov chain Monte Carlo (RJMCMC) was introduced by
[\[72\]](#ref72) as an extension to MCMC in which the dimension of the
model is uncertain and to be estimated. Reversible-jump belongs to a
family of trans-dimensional algorithms, and it has many applications,
including variable selection, model selection, mixture component
selection, and more. The RJ algorithm, here, is one of the simplest
possible implementations and is intended for variable selection and
Bayesian Model Averaging (BMA).

The Componentwise Hit-And-Run Metropolis (CHARM) algorithm was selected
in lucifer to be extended with reversible-jump. CHARM was selected
because it does not have tuning parameters, it is not adaptive (which
simplifies things with RJ), and it performs well. Even though it is a
componentwise algorithm, it does not evaluate all potential predictors
each iteration, but only those that are included. This novel combination
is Reversible-Jump (RJ) with Componentwise Hit-And-Run Metropolis
(CHARM). The optimal acceptance rate, and a good suggested acceptance
rate range, are unknowns.

RJ proceeds by making two proposals during each iteration. First, a
within-model move is proposed. This means that the model dimension does
not change, and the algorithm proceeds like a traditional CHARM
algorithm. Next, a between-models move is proposed, where a selectable
parameter is sampled, and its status in the model is changed. If this
selected parameter is in the current model, then RJ proposes a model
that excludes it. If this selected parameter is not in the current
model, then RJ proposes a model that includes it. RJ also includes a
user-specified binomial prior distribution for the expected size of the
model (the number of included, selectable parameters), as well as
user-specified prior probabilities for the inclusion of each of the
selectable parameters.

Behind the scenes, RJ keeps track of the most recent non-zero value for
each selectable parameter. If a parameter is currently excluded, then
its value is currently set to zero. When this parameter is proposed to
be included, the most recent non-zero value is used as the basis of the
proposal, rather than zero. In this way, an excluded parameter does not
have to travel back toward its previously included density, which may be
far from zero. However, if RJ begins updating after a previous run had
ended, then it will not begin again with this most recent non-zero
value. Please keep this in mind with this implementation of RJ.

RJ has five algorithm specifications. `bin.n` is the scalar size
parameter of the binomial prior distribution for model size, and is the
maximum size that RJ will explore. `bin.p` is the scalar probability
parameter of the binomial prior distribution for model size, and the
mean or median expected model size is `bin.n` x `bin.p.` `parm.p` is a
vector of parameter-inclusion probabilities that must be equal in length
to the number of initial values. `selectable` is a vector of indicators
(0 or 1) that indicate whether or not a parameter is selectable by
reversible-jump. When an element is zero, it is always in the model.
Finally, the `selected` vector contains indicators of whether or not
each parameter is in the model when RJ begins to update. All of these
specifications must receive an argument with exactly that name (such as
bin.n=bin.n, for the Consort function to recognize it, with the
exception of the `selected` specification.

RJ provides a sampler-based alternative to variable selection, compared
with the Bayesian Adaptive Lasso (BAL) or Stochastic Search Variable
Selection (SSVS), as two of many other approaches. Examples of BAL and
SSVS are in the Examples vignette. Advantages of RJ are that it is
easier for the user to apply to a given model than writing the
variable-selection code into the model, and RJ requires fewer
parameters, because variable-inclusion is handled by the specialized
sampler, rather than the model specification function. A disadvantage of
RJ is that other methods allow the user to freely change to other MCMC
algorithms, if the current algorithm is unsatisfactory.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, model indicators `selected`,
    prior inclusion probabilities `parm.p`
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Within-model move** (CHARM): for each included parameter \\j\\
        with `selected[j] = 1`:
        - Propose \\\theta^\*\_j\\ via Hit-And-Run in random direction;
          accept/reject via Metropolis step
    2.  **Between-models move**: select a parameter \\j\\ from the
        `selectable` set
        - If \`selected\[j\] = 1\$: propose exclusion (set \\\theta_j =
          0\\)
        - If \`selected\[j\] = 0\$: propose inclusion (restore last
          non-zero value)
        - Compute \\\alpha\\ including binomial model-size prior and
          parameter-inclusion prior `parm.p[j]`
        - With probability \\\alpha\\, update `selected[j]`
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\ and model indicators

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| bin.n | LIV | Maximum model size equals the number of available parameters. |
| bin.p | 0.5 | Equal prior probability of inclusion for each parameter. |
| parm.p | rep(0.5, LIV) | Equal selection probability for every parameter. |
| selectable | rep(1, LIV) | All parameters are eligible for inclusion/exclusion. |
| selected | rep(1, LIV) | Start from the full model. |

Default `Specs` for RJ. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RJ", Specs = list(bin.n = length(IV), bin.p = 0.5,
             parm.p = rep(0.5, length(IV)),
             selectable = rep(1L, length(IV)),
             selected = rep(1L, length(IV))),
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [CHARM](#charm)

### Non-Reversible Parallel Tempering

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The swap acceptance rate depends on the temperature spacing. Adaptive spacing targets equalized swap rates across all adjacent pairs. Within-chain acceptance follows RWM guidelines. |
| Applications | This algorithm targets multimodal posterior distributions where single-chain methods get trapped in local modes. It is the modern replacement for standard parallel tempering (MCMCMC), with provably faster mixing across the temperature ladder. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has four algorithm specifications. The temperature ladder is adapted automatically. |
| Final Algorithm? | User Discretion. NRPT may be used as a final algorithm after the temperature ladder has stabilized. |
| Proposal | Multivariate (parallel chains). |

Non-Reversible Parallel Tempering (NRPT) was introduced by
[\[73\]](#ref73) as a modification of standard parallel tempering that
achieves \\O(1)\\ mixing across the temperature ladder, compared to
\\O(K^2)\\ for reversible random-swap parallel tempering (where \\K\\ is
the number of temperatures). The key innovation is a deterministic
even-odd sweep schedule combined with direction variables. Standard
parallel tempering (MCMCMC in lucifer) selects random pairs of adjacent
chains and proposes swaps; this random selection creates a diffusive
process where information about the target distribution must perform a
random walk across the ladder. NRPT instead sweeps deterministically
through odd pairs \\(0,1), (2,3), \ldots\\ then even pairs \\(1,2),
(3,4), \ldots\\, maintaining a direction variable \\d_k \in \\+1, -1\\\\
for each adjacent pair. When a swap is rejected, the direction is
reversed. This non-reversibility creates a persistent drift that pushes
information through the ladder in \\O(K)\\ steps.

NRPT has four algorithm specifications: `n_temps` is the number of
temperatures (default \\\max(4, \text{CPUs})\\), `lambda` controls the
initial geometric spacing \\\beta_k = 1/(1 + \lambda k)\\ (default 1.0),
`Periodicity` is the adaptation frequency (default 100), and `n_steps`
is the number of within-temperature RWM steps between swap rounds
(default 1). The temperature ladder is adapted every `Periodicity`
iterations to equalize swap acceptance rates across all adjacent pairs,
following the strategy of [\[74\]](#ref74).

Each temperature chain targets the tempered distribution \\\pi_k(\theta)
\propto p(\theta) \cdot L(\theta \mid y)^{\beta_k}\\, where \\\beta_0 =
1\\ (the posterior) and \\\beta\_{K-1}\\ is close to 0 (the prior).
Within-chain moves use RWM with the proposal covariance scaled by
\\1/\sqrt{\beta_k}\\ to account for the flatter landscape at higher
temperatures. Only the cold chain (\\\beta = 1\\) is recorded. The swap
acceptance ratio for adjacent chains \\k\\ and \\k+1\\ is \\\alpha =
\exp\\\bigl\[(\beta_k - \beta\_{k+1})(LL\_{k+1} - LL_k)\bigr\]\\, where
\\LL\\ is the log-likelihood.

#### Algorithm

1.  **Initialize** \\K\\ chains at temperatures \\\beta_0 = 1 \> \beta_1
    \> \ldots \> \beta\_{K-1}\\; direction variables \\d_k = +1\\ for
    all pairs
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Within-chain RWM**: for each temperature \\k\\, run
        \\n\_{\text{steps}}\\ of RWM targeting \\\pi_k\\
    2.  **Odd-pair sweep**: for \\k = 0, 2, 4, \ldots\\, propose swap of
        chains \\k\\ and \\k+1\\; accept with probability \\\alpha\\; if
        rejected, flip \\d_k\\
    3.  **Even-pair sweep**: for \\k = 1, 3, 5, \ldots\\, same procedure
    4.  If \\t \bmod \text{Periodicity} = 0\\: adapt temperature ladder
        to equalize swap acceptance rates
3.  **Return** samples from the cold chain (\\\beta_0 = 1\\)

#### Default specifications

| Specification | Default | Justification                                      |
|:--------------|:--------|:---------------------------------------------------|
| n_temps       | 4       | Four temperatures: cold + three tempered replicas. |
| lambda        | 1.0     | Geometric tempering exponent.                      |
| Periodicity   | 100     | Adapt the temperature ladder every 100 iterations. |
| n_steps       | 1       | Single within-temperature step between swap moves. |

Default `Specs` for NRPT. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "NRPT", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MCMCMC](#mcmcmc)
- [RWM](#rwm)

### Simulated Tempering

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The within-temperature acceptance rate follows RWM guidelines. The temperature swap acceptance rate depends on the spacing and normalizing constant estimates. |
| Applications | This algorithm targets multimodal posterior distributions. Unlike parallel tempering which runs multiple chains simultaneously, simulated tempering runs a single chain that moves through a temperature ladder, requiring less memory and computation. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has two algorithm specifications. The normalizing constants at each temperature must be estimated adaptively, which adds complexity. |
| Final Algorithm? | User Discretion. Simulated Tempering may be used as a final algorithm once the temperature ladder and normalizing constant estimates have stabilized. |
| Proposal | Multivariate (single chain, multiple temperatures). |

Simulated Tempering was introduced by [\[75\]](#ref75) as a method for
sampling from multimodal distributions by augmenting the state space
with a temperature index \\k \in \\0, 1, \ldots, K-1\\\\. The chain
alternates between within-temperature moves (RWM targeting the tempered
distribution \\\pi_k(\theta) \propto \pi(\theta)^{\beta_k}\\) and
temperature moves (proposing \\k \to k \pm 1\\). The key challenge is
that the temperature move acceptance probability involves the ratio of
normalizing constants \\Z_k / Z\_{k+1}\\, which are unknown. These must
be estimated on the fly, typically by tracking the fraction of time
spent at each temperature and adjusting the weights \\\\w_k\\\\ to
achieve uniform occupancy. Unlike parallel tempering (which runs \\K\\
chains simultaneously), simulated tempering runs a single chain, making
it more memory-efficient but requiring the normalizing constant
estimates.

Simulated Tempering has two algorithm specifications: `n_temps` is the
number of temperatures in the ladder (default 5), and `lambda` is the
temperature spacing exponent controlling the geometric schedule
\\\beta_k = (1 - k/(K-1))^{\lambda}\\ (default 1.5).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, temperature index \\k = 0\\,
    weights \\w_k = 0\\ for all \\k\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Within-temperature move**: propose \\\theta^\* \sim q(\cdot
        \mid \theta\_{t-1})\\; accept with probability \\\min(1,
        \pi(\theta^\*)^{\beta_k} / \pi(\theta\_{t-1})^{\beta_k})\\
    2.  **Temperature move**: propose \\k^\* = k \pm 1\\ (uniformly);
        accept with probability \\\min\bigl(1,
        \exp\bigl\[(\beta\_{k^\*} - \beta_k) \log \pi(\theta) +
        w\_{k^\*} - w_k\bigr\]\bigr)\\
    3.  Adapt weights: \\w_k \leftarrow w_k - \eta_t \cdot
        (\text{indicator}\[k_t = k\] - 1/K)\\
    4.  Record \\\theta_t\\ only when \\k = 0\\ (the target temperature)
3.  **Return** samples collected at \\\beta_0 = 1\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| n_temps | 5 | Five temperature levels. |
| lambda | 1.5 | Geometric spacing exponent; 1.5 produces a ladder denser at the cold end. |

Default `Specs` for SimTemp. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "SimTemp", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MCMCMC](#mcmcmc)
- [Non-Reversible Parallel Tempering](#nrpt)
- [Non-Reversible Simulated Tempering](#nrst)

### Non-Reversible Simulated Tempering

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The within-temperature acceptance rate follows RWM guidelines. Temperature traversal is faster than reversible simulated tempering due to the persistent direction variable. |
| Applications | This algorithm targets multimodal posterior distributions. It combines the single-chain memory efficiency of simulated tempering with the non-reversible acceleration of NRPT, achieving O(1) round-trip times through the temperature ladder. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has two algorithm specifications. The non-reversible direction variable and normalizing constant adaptation add complexity over standard simulated tempering. |
| Final Algorithm? | User Discretion. NRST may be used as a final algorithm once the temperature ladder stabilizes. |
| Proposal | Multivariate (single chain, multiple temperatures). |

Non-Reversible Simulated Tempering (NRST) was introduced by
[\[35\]](#ref35) as a modification of simulated tempering that adds a
direction variable \\d \in \\+1, -1\\\\ to create persistent motion
through the temperature ladder. In standard simulated tempering, the
temperature index performs a random walk, requiring \\O(K^2)\\ steps for
a round trip; NRST’s direction variable creates a persistent drift,
reducing the round-trip time to \\O(K)\\. At each step, the chain
proposes to move in the current direction: \\k^\* = k + d\\. If the
temperature move is rejected, the direction is reversed (\\d \leftarrow
-d\\), just as in Non-Reversible Parallel Tempering (NRPT). The
difference from NRPT is that NRST uses a single chain rather than \\K\\
parallel chains, trading parallelism for memory efficiency. The
normalizing constant estimates \\\\w_k\\\\ are adapted identically to
standard simulated tempering.

NRST has two algorithm specifications: `n_temps` is the number of
temperatures (default 5), and `lambda` is the temperature spacing
exponent (default 1.5).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, temperature index \\k = 0\\,
    direction \\d = +1\\, weights \\w_k = 0\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  **Within-temperature move**: RWM step targeting \\\pi_k(\theta)
        \propto \pi(\theta)^{\beta_k}\\
    2.  **Temperature move**: propose \\k^\* = k + d\\
        - If \\k^\* \in \\0, \ldots, K-1\\\\: accept with probability
          \\\min\bigl(1, \exp\bigl\[(\beta\_{k^\*} - \beta_k)\log
          \pi(\theta) + w\_{k^\*} - w_k\bigr\]\bigr)\\
        - If accepted: \\k \leftarrow k^\*\\
        - If rejected or out of bounds: \\d \leftarrow -d\\
    3.  Adapt weights \\w_k\\
    4.  Record \\\theta_t\\ only when \\k = 0\\
3.  **Return** samples collected at \\\beta_0 = 1\\

#### Default specifications

| Specification | Default | Justification            |
|:--------------|:--------|:-------------------------|
| n_temps       | 5       | Five temperature levels. |
| lambda        | 1.5     | Geometric ladder.        |

Default `Specs` for NRST. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "NRST", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Simulated Tempering](#simtemp)
- [Non-Reversible Parallel Tempering](#nrpt)
- [Wang-Landau](#wanglandau)

### Wang-Landau

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate varies; the algorithm targets a flat histogram across energy bins, so the effective acceptance depends on the current weight estimates. |
| Applications | This algorithm targets multimodal distributions and is particularly useful for estimating the density of states (normalizing constant as a function of energy). It originated in statistical physics for computing partition functions. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has four algorithm specifications. The flatness criterion and modification factor schedule require understanding. |
| Final Algorithm? | User Discretion. Wang-Landau should be used in two stages: a learning phase where the weights are adapted, followed by a production phase with fixed weights. |
| Proposal | Multivariate. |

The Wang-Landau algorithm was introduced by [\[76\]](#ref76) as an
adaptive method for estimating the density of states in statistical
physics. In the MCMC context, it provides a flat-histogram approach to
sampling multimodal distributions. The energy space is partitioned into
bins, and a weight function \\g(E)\\ (an estimate of the log density of
states) is adapted so that the modified target \\\pi(\theta) /
g(E(\theta))\\ has approximately uniform energy distribution. At each
iteration, a Metropolis move targeting the modified distribution is
performed, and the weight of the current bin is updated: \\g(E)
\leftarrow g(E) \cdot f\\, where \\f \> 1\\ is a modification factor. A
histogram \\H(E)\\ tracks visits to each bin; when the histogram is
“flat” (all entries within a tolerance of the mean), the modification
factor is reduced (\\f \leftarrow \sqrt{f}\\) and the histogram is
reset. This process is repeated until \\f\\ is close to 1, at which
point the weight estimates have converged and the sampler explores all
energy levels uniformly.

Wang-Landau has four algorithm specifications: `n_bins` is the number of
energy bins (default 50), `f_init` is the initial modification factor
(default \\e \approx 2.718\\), `f_min` is the final modification factor
threshold for convergence (default 1.001), and `flatness` is the
flatness criterion for the histogram (default 0.8, meaning all bins must
have at least 80% of the mean count).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, weight estimates \\g(E) =
    0\\ for all bins, histogram \\H(E) = 0\\, modification factor \\f =
    f\_{\text{init}}\\
2.  **While** \\f \> f\_{\min}\\:
    1.  Propose \\\theta^\* \sim q(\cdot \mid \theta)\\
    2.  Accept with probability \\\min\bigl(1, \exp\bigl\[\log
        \pi(\theta^\*) - \log \pi(\theta) - g(E^\*) +
        g(E)\bigr\]\bigr)\\
    3.  Update weight: \\g(E\_{\text{current}}) \leftarrow
        g(E\_{\text{current}}) + \log f\\
    4.  Update histogram: \\H(E\_{\text{current}}) \leftarrow
        H(E\_{\text{current}}) + 1\\
    5.  If histogram is flat (\$H / (H) \> \$ flatness): \\f \leftarrow
        \sqrt{f}\\, reset \\H\\
3.  **Production phase**: run MCMC with fixed \\g(E)\\ weights
4.  **Return** production samples reweighted to recover \\\pi(\theta)\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| n_bins | 50 | Energy histogram resolution; 50 bins balance variance and visit time. |
| f_init | 2.718 | Initial multiplicative weight = e. |
| f_min | 1.001 | Termination threshold for the multiplicative weight. |
| flatness | 0.8 | Histogram flatness criterion of Wang and Landau (2001). |

Default `Specs` for WL. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "WL", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Simulated Tempering](#simtemp)
- [Non-Reversible Simulated Tempering](#nrst)
- [MCMCMC](#mcmcmc)

## Piecewise-deterministic samplers explore without random walks

Piecewise-deterministic Markov processes (PDMPs) represent a
fundamentally different approach to MCMC. Instead of proposing random
jumps and accepting or rejecting them, PDMPs move deterministically
along straight lines (or other curves) and stochastically change
direction at random times. The trajectory between direction changes is
continuous and deterministic, producing a piecewise-linear (or
piecewise-smooth) sample path that can be subsampled at any resolution.
PDMPs require gradient information to determine the rate of direction
changes.

### Zig-Zag Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | Not applicable in the traditional sense. The sampler accepts all moves; the stochastic element is the timing and selection of velocity flips (bounces). |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is particularly efficient when componentwise gradient evaluations are cheap, because each bounce event requires only one partial derivative. |
| Difficulty | This algorithm is difficult for a beginner. It has four algorithm specifications and operates on a fundamentally different principle than standard MCMC. The output requires subsampling from a continuous trajectory. |
| Final Algorithm? | User Discretion. The zig-zag sampler produces asymptotically exact samples, but the discretized output requires care in choosing the subsampling rate. |
| Proposal | Continuous trajectory (event-driven). |

The Zig-Zag sampler was introduced by [\[77\]](#ref77) as a
piecewise-deterministic Markov process (PDMP) for Bayesian computation.
A particle moves along straight lines with velocity \\v \in \\-1,
+1\\^d\\ (each component is \\\pm 1\\). At random times, one velocity
component flips sign (a “bounce”). The rate of bouncing in component
\\j\\ is \\\lambda_j(x, v) = \max(0, -v_j \cdot \partial \log \pi(x) /
\partial x_j)\\: the particle bounces when it is moving downhill in that
component. Between bounces, the particle traces a piecewise-linear
trajectory that is a valid sample from the target distribution.
Additionally, velocity refreshment events (at rate `refresh_rate`)
randomize individual velocity components to ensure ergodicity.

Zig-Zag has four algorithm specifications: `T_max` is the total
trajectory simulation time (auto-scaled if NULL), `refresh_rate` is the
rate of velocity refreshment events (default 0.1), `excess_rate` is the
upper bound multiplier for Poisson thinning (default 1.5), and `method`
is either `"zigzag"` (coordinate-wise flips, default) or `"bps"` (bouncy
particle sampler with full velocity reflection).

The zig-zag sampler uses Poisson thinning to simulate the inhomogeneous
Poisson process governing bounce times. At each step, a candidate event
time is drawn from a homogeneous Poisson process with rate equal to the
upper bound of the total bounce rate times `excess_rate`. The event is
then accepted or rejected based on the ratio of the actual rate to the
upper bound at the candidate time. This is exact (no approximation error
from the thinning), but efficiency depends on how tight the upper bound
is.

The advantage of the zig-zag sampler over standard MCMC is that between
bounces the sampler moves deterministically, producing
autocorrelation-free samples when subsampled at the correct rate. Each
bounce requires only one componentwise gradient evaluation, making it
potentially very efficient for problems where full gradient computation
is expensive but individual partial derivatives are cheap. The
disadvantage is the fundamentally different output format (a skeleton of
event times and positions rather than a sequence of states), which
requires interpolation for standard MCMC diagnostics.

The C++ backend uses `compute_gradient()` from `sampler_common.h`, which
dispatches to either a user-supplied gradient function or the built-in
finite-difference engine, eliminating the per-component R loop overhead.
Bounce component selection uses a cumulative probability scan in C++
rather than R’s [`sample.int()`](https://rdrr.io/r/base/sample.html).

#### Algorithm

1.  **Initialize** position \\x_0\\, velocity \\v = (\pm 1, \ldots, \pm
    1)\\ uniformly at random
2.  **Repeat** until trajectory length exceeds \\T\_{\max}\\:
    1.  Compute componentwise bounce rates \\\lambda_j = \max(0, -v_j
        \cdot \partial_j \log \pi(x))\\
    2.  Set total rate \\\Lambda = \sum_j \lambda_j + d \cdot
        r\_{\text{refresh}}\\; draw \\\tau \sim \text{Exp}(c \cdot
        \Lambda)\\ where \\c\\ = excess_rate
    3.  Move deterministically: \\x \leftarrow x + v \cdot \tau\\
    4.  Compute actual rates at new position; accept event with
        probability \\\Lambda\_{\text{actual}} / (c \cdot \Lambda)\\
    5.  If accepted: select component \\j\\ proportional to rates; if
        bounce, flip \\v_j\\; if refresh, randomize \\v_j\\
3.  **Subsample** the piecewise-linear trajectory at equally-spaced time
    points
4.  **Return** subsampled positions

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| T_max | NULL | Auto-derived per iteration from the local rate. |
| refresh_rate | 0.1 | Velocity refreshment rate; 0.1 keeps the process ergodic without too many refreshments. |
| excess_rate | 1.5 | Excess Poisson rate for thinning. |
| method | “zigzag” | Standard Zig-Zag dynamics. |

Default `Specs` for ZigZag. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "ZigZag", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [NUTS](#nuts)
- [Barker](#barker)

### Bouncy Particle Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | Not applicable. All events are accepted; the stochastic element is the timing and geometry of velocity reflections. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is especially efficient in high dimensions where the gradient is available, because the bounce rate depends on the alignment between velocity and gradient. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has three algorithm specifications and operates on PDMP principles requiring trajectory subsampling. |
| Final Algorithm? | User Discretion. BPS produces asymptotically exact samples; output quality depends on subsampling rate and trajectory length. |
| Proposal | Continuous trajectory (event-driven). |

The Bouncy Particle Sampler (BPS) was introduced by [\[78\]](#ref78) as
a piecewise-deterministic Markov process for Bayesian computation.
Unlike the Zig-Zag sampler, which flips individual velocity components,
BPS operates with a continuous velocity vector \\v \in \mathbb{R}^d\\.
The particle moves along the ray \\x(t) = x_0 + tv\\ until a bounce
event occurs, at which point the velocity is reflected off the gradient
of the negative log-target: \\v \leftarrow v - 2\frac{\langle \nabla
U(x), v \rangle}{\\\nabla U(x)\\^2} \nabla U(x)\\, where \\U(x) = -\log
\pi(x)\\. The bounce rate is \\\lambda(x, v) = \max(0, \langle \nabla
U(x), v \rangle)\\, so the particle bounces only when it moves into
regions of decreasing probability. Velocity refreshment events at rate
`refresh_rate` replace \\v\\ with a draw from \\\mathcal{N}(0, I_d)\\,
ensuring ergodicity and preventing the particle from cycling in
subspaces.

BPS has three algorithm specifications: `refresh_rate` is the rate of
velocity refreshment events (default 0.1), `excess_rate` is the upper
bound multiplier for Poisson thinning (default 1.5), and `T_max` is the
total trajectory simulation time (auto-scaled if NULL). The bounce
events are simulated via Poisson thinning with the same machinery as the
Zig-Zag sampler.

#### Algorithm

1.  **Initialize** position \\x_0\\, velocity \\v \sim \mathcal{N}(0,
    I_d)\\
2.  **Repeat** until trajectory length exceeds \\T\_{\max}\\:
    1.  Compute bounce rate \\\lambda(x, v) = \max(0, \langle \nabla
        U(x), v \rangle)\\
    2.  Set total rate \\\Lambda = \lambda + r\_{\text{refresh}}\\; draw
        \\\tau \sim \text{Exp}(c \cdot \Lambda)\\ where \\c\\ =
        excess_rate
    3.  Move deterministically: \\x \leftarrow x + v \cdot \tau\\
    4.  Accept event with probability \\\Lambda\_{\text{actual}} / (c
        \cdot \Lambda)\\
    5.  If bounce event: reflect \\v \leftarrow v - 2 \frac{\langle
        \nabla U(x), v \rangle}{\\\nabla U(x)\\^2} \nabla U(x)\\
    6.  If refresh event: draw \\v \sim \mathcal{N}(0, I_d)\\
3.  **Subsample** the piecewise-linear trajectory at equally-spaced time
    points
4.  **Return** subsampled positions

#### Default specifications

| Specification | Default | Justification          |
|:--------------|:--------|:-----------------------|
| refresh_rate  | 0.1     | Velocity refresh rate. |
| excess_rate   | 1.5     | Excess thinning rate.  |

Default `Specs` for BPS. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "BPS", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Zig-Zag Sampler](#zigzag)
- [Boomerang Sampler](#boomerang)
- [Randomized HMC](#rhmc)

### Boomerang Sampler

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | Not applicable. All events are accepted; the stochastic element is the timing of velocity refreshments and bounces along elliptical arcs. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is especially useful when the target is approximately Gaussian, because the deterministic dynamics follow elliptical orbits that naturally match such geometry. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has three algorithm specifications and requires understanding of Hamiltonian-like elliptical dynamics within the PDMP framework. |
| Final Algorithm? | User Discretion. The Boomerang sampler produces asymptotically exact samples and can outperform BPS on near-Gaussian targets. |
| Proposal | Continuous trajectory (elliptical arcs). |

The Boomerang Sampler was introduced by [\[79\]](#ref79) as a
piecewise-deterministic Markov process that generalizes the Bouncy
Particle Sampler by replacing straight-line deterministic dynamics with
elliptical (Hamiltonian) dynamics. Between bounce events, the particle
follows the orbit of a quadratic Hamiltonian \\H(x, v) =
\frac{1}{2}\\v\\^2 + \frac{1}{2}x^\top \Sigma^{-1} x\\, which traces
ellipses in \\(x, v)\\ space. This means the deterministic trajectory
curves back toward the origin (hence “boomerang”), which is advantageous
when the target is approximately Gaussian. The bounce rate is
\\\lambda(x, v) = \max(0, \langle \nabla U(x) - \Sigma^{-1} x, v
\rangle)\\, where \\\Sigma^{-1} x\\ is the gradient of the reference
quadratic that is already accounted for by the elliptical dynamics. When
the target is exactly \\\mathcal{N}(0, \Sigma)\\, the bounce rate is
identically zero and the sampler moves without any bounces, achieving
perfect efficiency.

Boomerang has three algorithm specifications: `refresh_rate` is the rate
of velocity refreshment events (default 0.1), `excess_rate` is the upper
bound multiplier for Poisson thinning (default 1.5), and `T_max` is the
total trajectory simulation time (auto-scaled if NULL).

#### Algorithm

1.  **Initialize** position \\x_0\\, velocity \\v \sim \mathcal{N}(0,
    I_d)\\, reference precision \\\Sigma^{-1}\\
2.  **Repeat** until trajectory length exceeds \\T\_{\max}\\:
    1.  Compute bounce rate \\\lambda = \max(0, \langle \nabla U(x) -
        \Sigma^{-1} x, v \rangle)\\
    2.  Set total rate \\\Lambda = \lambda + r\_{\text{refresh}}\\; draw
        candidate time \\\tau\\ via Poisson thinning
    3.  Evolve deterministically along ellipse: \\(x, v) \leftarrow
        (\cos(\tau)x + \sin(\tau)\Sigma v,\\ -\Sigma^{-1}\sin(\tau)x +
        \cos(\tau)v)\\
    4.  Accept event with probability \\\Lambda\_{\text{actual}} / (c
        \cdot \Lambda)\\
    5.  If bounce: reflect \\v \leftarrow v - 2 \frac{\langle \nabla
        U(x) - \Sigma^{-1}x, v \rangle}{\\\nabla U(x) -
        \Sigma^{-1}x\\^2} (\nabla U(x) - \Sigma^{-1}x)\\
    6.  If refresh: draw \\v \sim \mathcal{N}(0, I_d)\\
3.  **Subsample** the trajectory at equally-spaced time points
4.  **Return** subsampled positions

#### Default specifications

| Specification | Default | Justification          |
|:--------------|:--------|:-----------------------|
| refresh_rate  | 0.1     | Velocity refresh rate. |
| excess_rate   | 1.5     | Excess thinning rate.  |

Default `Specs` for Boomerang. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Boomerang", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Bouncy Particle Sampler](#bps)
- [Zig-Zag Sampler](#zigzag)
- [Randomized HMC](#rhmc)

### Randomized HMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on the integration accuracy. With a symplectic integrator and moderate step sizes, acceptance rates above 65% are typical. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It bridges the gap between standard HMC (fixed trajectory length) and PDMP samplers (continuous-time), and is particularly useful when trajectory length tuning is difficult. |
| Difficulty | This algorithm is moderately easy for a beginner. It has two algorithm specifications and the randomized trajectory length reduces sensitivity to periodic orbits. |
| Final Algorithm? | Yes. Randomized HMC may be used as a final algorithm. |
| Proposal | Multivariate (Hamiltonian dynamics). |

Randomized HMC was introduced by [\[80\]](#ref80) as a modification of
standard Hamiltonian Monte Carlo that draws the trajectory length from
an exponential distribution rather than fixing it. In standard HMC, the
number of leapfrog steps \\L\\ is a fixed hyperparameter; poor choices
of \\L\\ can cause near-periodic orbits that degrade mixing. Randomized
HMC eliminates this problem by drawing \\L \sim
\text{Poisson}(\bar{L})\\ (or equivalently, running the trajectory for a
random time \\T \sim \text{Exp}(1/\bar{L})\\). This randomization breaks
periodicity and produces a sampler that, in the continuous-time limit,
converges to a PDMP that is ergodic without velocity refreshment. The
randomization also provides a theoretical advantage: Bou-Rabee and
Sanz-Serna show that the randomized variant achieves faster convergence
to equilibrium than the deterministic variant for certain target
classes.

Randomized HMC has two algorithm specifications: `epsilon` is the
leapfrog step size, and `mean_L` is the mean of the Poisson distribution
from which the number of leapfrog steps is drawn (default 10). The
Metropolis acceptance step uses the standard HMC accept/reject criterion
based on the Hamiltonian error.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\, mean
    trajectory length \\\bar{L}\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum \\p \sim \mathcal{N}(0, I_d)\\
    2.  Draw number of leapfrog steps \\L \sim \max(1,
        \text{Poisson}(\bar{L}))\\
    3.  Integrate \\L\\ leapfrog steps: \\(\theta^\*, p^\*) =
        \text{Leapfrog}(\theta\_{t-1}, p, \epsilon, L)\\
    4.  Compute \\\log \alpha = \log \pi(\theta^\*) -
        \frac{1}{2}\\p^\*\\^2 - \log \pi(\theta\_{t-1}) +
        \frac{1}{2}\\p\\^2\\
    5.  With probability \\\min(1, \exp(\log \alpha))\\, accept
        \\\theta^\*\\; otherwise retain \\\theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step size. |
| mean_L | 10 | Mean of the geometric distribution of leapfrog lengths. |

Default `Specs` for RHMC. {.table style="width:100%;"}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RHMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [NUTS](#nuts)
- [Bouncy Particle Sampler](#bps)

## Riemannian and geometric MCMC explore the curvature of the target distribution

Riemannian and geometric MCMC methods exploit the curvature of the
target distribution’s parameter space to achieve more efficient
sampling. Instead of treating the parameter space as flat Euclidean
space (as standard HMC does with a constant mass matrix), these methods
equip the space with a position-dependent metric tensor \\G(\theta)\\
derived from the Fisher information matrix or the Hessian of the
negative log-posterior. The resulting dynamics follow geodesics on the
statistical manifold, automatically adapting step sizes and directions
to the local geometry. This is particularly advantageous for targets
with strong correlations, funnel geometries, or varying curvature across
the parameter space.

### Riemannian Manifold HMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is typically 60-80%. The generalized leapfrog integrator introduces additional numerical error relative to standard leapfrog, requiring smaller step sizes. |
| Applications | This algorithm is applicable to targets with strong parameter correlations, funnel geometries, or varying curvature. It excels on hierarchical models where the Neal funnel phenomenon causes standard HMC to fail. |
| Difficulty | This algorithm is difficult for a user. It has six algorithm specifications and requires the metric tensor (Fisher information or Hessian) to be computable at every point. The implicit generalized leapfrog integrator involves fixed-point iterations. |
| Final Algorithm? | User Discretion. RMHMC may be used as a final algorithm once the metric mode and step size are tuned. |
| Proposal | Multivariate (Riemannian Hamiltonian dynamics). |

Riemannian Manifold Hamiltonian Monte Carlo (RMHMC) was introduced by
[\[81\]](#ref81) as an extension of HMC that replaces the constant mass
matrix \\M\\ with a position-dependent metric tensor \\G(\theta)\\,
typically chosen as the Fisher information matrix or the negative
Hessian of the log-posterior. The Hamiltonian becomes \\H(\theta, p) =
-\log \pi(\theta) + \frac{1}{2}\log\|G(\theta)\| + \frac{1}{2}p^\top
G(\theta)^{-1} p\\, and the resulting dynamics follow geodesics on the
statistical manifold. Because \\G(\theta)\\ varies with position, the
standard leapfrog integrator is no longer sufficient; RMHMC uses a
generalized leapfrog (implicit midpoint) integrator that requires
fixed-point iterations at each step. The SoftAbs metric [\[82\]](#ref82)
provides a smooth, positive-definite regularization of the Hessian that
avoids the singularities that can occur with the raw Hessian, using
\\G(\theta) = H \cdot \text{coth}(\alpha H)\\ where \\H\\ is the Hessian
and \\\alpha\\ controls the softness of the absolute value.

RMHMC has six algorithm specifications: `epsilon` is the leapfrog step
size, `L` is the number of leapfrog steps, `metric` selects the metric
tensor type (`"softabs"`, `"hessian"`, or `"identity"`, default
`"softabs"`), `softabs_alpha` is the SoftAbs regularization parameter
(default 1.0), `fixedpoint_maxiter` is the maximum number of fixed-point
iterations in the implicit integrator (default 8), and `fixedpoint_tol`
is the convergence tolerance for fixed-point iterations (default
\\10^{-6}\\).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    trajectory length \\L\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute metric \\G(\theta\_{t-1})\\ and its Cholesky factor
        \\L_G\\
    2.  Draw momentum \\p \sim \mathcal{N}(0, G(\theta\_{t-1}))\\
    3.  **For** \\\ell = 1, \ldots, L\\ (generalized leapfrog):
        - Half-step momentum: solve \\\hat{p} = p - \frac{\epsilon}{2}
          \nabla\_\theta H(\theta, p)\\ via fixed-point iteration
        - Full-step position: solve \\\hat{\theta} = \theta +
          \frac{\epsilon}{2}\[G(\theta)^{-1}\hat{p} +
          G(\hat{\theta})^{-1}\hat{p}\]\\ via fixed-point iteration
        - Half-step momentum: \\p' = \hat{p} - \frac{\epsilon}{2}
          \nabla\_\theta H(\hat{\theta}, \hat{p})\\
    4.  Compute \\\log \alpha = H(\theta\_{t-1}, p_0) - H(\theta^\*,
        p^\*)\\
    5.  Accept with probability \\\min(1, \exp(\log \alpha))\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step size. |
| L | 5 | Five generalized leapfrog steps per trajectory. |
| metric | “softabs” | SoftAbs regularization of the Hessian; numerically robust default for arbitrary log-posteriors. |
| softabs_alpha | 1.0 | Sharpness parameter of the SoftAbs map. |
| fixedpoint_maxiter | 8 | At most eight fixed-point iterations per implicit leapfrog step. |
| fixedpoint_tol | 1e-6 | Convergence tolerance for the fixed-point iteration. |

Default `Specs` for RMHMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "RMHMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [Lagrangian Monte Carlo](#lmc)
- [Magnetic HMC](#mhmc)

### Lagrangian Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is typically 60-80%, similar to RMHMC. The explicit integrator avoids fixed-point iterations, potentially allowing larger step sizes. |
| Applications | This algorithm is applicable to the same class of problems as RMHMC (curved parameter spaces, funnel geometries, hierarchical models) but may be more efficient when the metric tensor is expensive to invert. |
| Difficulty | This algorithm is difficult for a user. It has four algorithm specifications and requires a computable metric tensor. The Lagrangian formulation avoids implicit integration but introduces Christoffel symbols. |
| Final Algorithm? | User Discretion. LMC may be used as a final algorithm once the metric mode and step size are tuned. |
| Proposal | Multivariate (Lagrangian dynamics). |

Lagrangian Monte Carlo (LMC) was introduced by [\[83\]](#ref83) as an
alternative to RMHMC that uses the Lagrangian (rather than Hamiltonian)
formulation of mechanics on the statistical manifold. While RMHMC works
with the Hamiltonian \\H(\theta, p)\\ and requires an implicit
generalized leapfrog integrator involving costly fixed-point iterations,
LMC works directly with the Lagrangian \\\mathcal{L}(\theta,
\dot{\theta}) = \frac{1}{2}\dot{\theta}^\top G(\theta)\dot{\theta} +
\log \pi(\theta) - \frac{1}{2}\log\|G(\theta)\|\\ and the Euler-Lagrange
equations. The key advantage is that the resulting integrator can be
made explicit (no fixed-point iterations) by using a carefully designed
splitting scheme, which can significantly reduce the per-step cost when
the metric tensor is expensive. The dynamics follow the same geodesics
as RMHMC, so the two methods sample from the same extended target
distribution.

LMC has four algorithm specifications: `epsilon` is the step size, `L`
is the number of integration steps, `metric` selects the metric tensor
type (`"softabs"`, `"hessian"`, or `"identity"`, default `"softabs"`),
and `softabs_alpha` is the SoftAbs regularization parameter (default
1.0).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, velocity
    \\\dot{\theta}\_0\\, step size \\\epsilon\\, trajectory length \\L\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute metric \\G(\theta\_{t-1})\\
    2.  Draw velocity \\\dot{\theta} \sim \mathcal{N}(0,
        G(\theta\_{t-1})^{-1})\\
    3.  **For** \\\ell = 1, \ldots, L\\ (explicit Lagrangian
        integrator):
        - Compute Christoffel symbols \\\Gamma^k\_{ij}(\theta)\\ from
          \\G(\theta)\\
        - Update velocity: \\\dot{\theta} \leftarrow \dot{\theta} +
          \frac{\epsilon}{2}\bigl\[G(\theta)^{-1}\nabla \log
          \pi(\theta) - \Gamma(\theta, \dot{\theta})\bigr\]\\
        - Update position: \\\theta \leftarrow \theta + \epsilon
          \dot{\theta}\\
        - Recompute and update velocity with new Christoffel symbols
    4.  Compute acceptance probability from the Lagrangian energy
        difference
    5.  Accept with probability \\\min(1, \exp(\Delta \mathcal{L}))\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default                  | Justification                   |
|:--------------|:-------------------------|:--------------------------------|
| epsilon       | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step size.           |
| L             | 5                        | Five Lagrangian leapfrog steps. |
| metric        | “softabs”                | SoftAbs metric.                 |
| softabs_alpha | 1.0                      | SoftAbs sharpness.              |

Default `Specs` for LMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "LMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Riemannian Manifold HMC](#rmhmc)
- [Magnetic HMC](#mhmc)
- [HMC](#hmc)

### Magnetic HMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is typically 60-80%. The magnetic field does not change the energy, so acceptance depends only on the integrator accuracy, as in standard HMC. |
| Applications | This algorithm is applicable to targets with multimodal or twisted geometry. The magnetic field induces non-reversible dynamics that can break the symmetric exploration patterns of standard HMC, improving mixing across modes. |
| Difficulty | This algorithm is difficult for a user. It has five algorithm specifications including the field strength, which requires tuning. The antisymmetric magnetic field matrix must be constructed and integrated. |
| Final Algorithm? | User Discretion. Magnetic HMC may be used as a final algorithm once the field strength and step size are tuned. |
| Proposal | Multivariate (magnetic Hamiltonian dynamics). |

Magnetic Hamiltonian Monte Carlo was introduced by [\[84\]](#ref84) as a
modification of HMC that adds an antisymmetric “magnetic field” matrix
\\F\\ to the equations of motion. The modified Hamiltonian dynamics
become \\\dot{\theta} = M^{-1}p\\ and \\\dot{p} = \nabla \log
\pi(\theta) + F \cdot M^{-1}p\\, where \\F = -F^\top\\ is an
antisymmetric matrix representing the magnetic field. Because \\F\\ is
antisymmetric, \\p^\top F M^{-1} p = 0\\, so the magnetic force does no
work and the Hamiltonian is preserved exactly (up to integrator error).
The magnetic field induces non-reversible, curved trajectories that can
explore the target distribution more efficiently than the straight-line
reversible trajectories of standard HMC. In particular, the curved paths
can spiral around modes and navigate between modes more effectively when
the field strength is properly tuned. The method uses a modified
leapfrog integrator that accounts for the magnetic force term.

Magnetic HMC has five algorithm specifications: `epsilon` is the
leapfrog step size, `L` is the number of leapfrog steps, `metric`
selects the mass matrix type (`"softabs"`, `"hessian"`, or `"identity"`,
default `"softabs"`), `softabs_alpha` is the SoftAbs regularization
parameter (default 1.0), and `field_strength` controls the magnitude of
the antisymmetric magnetic field (default 1.0). The field matrix is
constructed as \\F = s \cdot (J - J^\top)\\ where \\J\\ is a random
matrix and \\s\\ is the field strength.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, magnetic field \\F\\, step
    size \\\epsilon\\, trajectory length \\L\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum \\p \sim \mathcal{N}(0, M)\\
    2.  **For** \\\ell = 1, \ldots, L\\ (magnetic leapfrog):
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\[\nabla \log \pi(\theta) + F M^{-1} p\]\\
        - Full-step position: \\\theta \leftarrow \theta + \epsilon
          M^{-1} p\\
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\[\nabla \log \pi(\theta) + F M^{-1} p\]\\
    3.  Compute \\\log \alpha = H(\theta\_{t-1}, p_0) - H(\theta^\*,
        p^\*)\\
    4.  Accept with probability \\\min(1, \exp(\log \alpha))\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step size. |
| L | 5 | Five magnetic leapfrog steps. |
| metric | “softabs” | SoftAbs metric. |
| softabs_alpha | 1.0 | SoftAbs sharpness. |
| field_strength | 1.0 | Magnetic field strength of unit magnitude. |

Default `Specs` for MHMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "MHMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Riemannian Manifold HMC](#rmhmc)
- [HMC](#hmc)
- [Relativistic Monte Carlo](#relativistic)

### Relativistic Monte Carlo

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The target acceptance rate is typically 65-80%. The relativistic kinetic energy produces heavier tails in the momentum distribution, improving robustness to step size choice. |
| Applications | This algorithm is applicable to smooth, differentiable targets. It is particularly robust for targets with heavy tails or regions of low curvature where standard HMC momentum can overshoot. |
| Difficulty | This algorithm is moderately easy for a beginner. It has three algorithm specifications. The speed-of-light parameter c controls the transition between Newtonian (large c) and relativistic (small c) regimes. |
| Final Algorithm? | Yes. Relativistic MC may be used as a final algorithm. |
| Proposal | Multivariate (relativistic Hamiltonian dynamics). |

Relativistic Monte Carlo was introduced by [\[85\]](#ref85) as a
modification of HMC that replaces the standard quadratic kinetic energy
\\K(p) = \frac{1}{2}p^\top M^{-1}p\\ with a relativistic kinetic energy
\\K(p) = c\sqrt{\\p\\^2 + m^2c^2} - mc^2\\, inspired by special
relativity. The crucial property is that the velocity \\\dot{\theta} =
\nabla_p K = cp / \sqrt{\\p\\^2 + m^2c^2}\\ is bounded by the speed of
light \\c\\, regardless of how large the momentum becomes. This imposes
a natural speed limit on the sampler, preventing the catastrophically
large jumps that can occur in standard HMC when the gradient is large or
the step size is too big. The momentum marginal distribution becomes a
multivariate hyperbolic distribution (heavier tails than Gaussian),
which makes the sampler more robust to targets with heavy tails or
varying curvature. When \\c \to \infty\\, the method reduces to standard
HMC.

Relativistic MC has three algorithm specifications: `epsilon` is the
leapfrog step size, `L` is the number of leapfrog steps, and `c` is the
speed of light parameter (default 10.0). Smaller values of \\c\\ impose
a tighter speed limit and produce more conservative exploration.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    trajectory length \\L\\, speed of light \\c\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw momentum from the relativistic distribution: \\p \sim
        \pi(p) \propto \exp(-c\sqrt{\\p\\^2 + c^2})\\
    2.  **For** \\\ell = 1, \ldots, L\\ (relativistic leapfrog):
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\nabla \log \pi(\theta)\\
        - Full-step position: \\\theta \leftarrow \theta + \epsilon
          \cdot cp / \sqrt{\\p\\^2 + c^2}\\
        - Half-step momentum: \\p \leftarrow p +
          \frac{\epsilon}{2}\nabla \log \pi(\theta)\\
    3.  Compute \\\log \alpha = H(\theta\_{t-1}, p_0) - H(\theta^\*,
        p^\*)\\ with relativistic \\H\\
    4.  Accept with probability \\\min(1, \exp(\log \alpha))\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Stan-style step size. |
| L | 10 | Ten relativistic leapfrog steps. |
| c | 10.0 | Speed-of-light cap; large values recover ordinary HMC. |

Default `Specs` for Relativistic. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "Relativistic", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [HMC](#hmc)
- [Magnetic HMC](#mhmc)
- [Riemannian Manifold HMC](#rmhmc)

## Constraint-handling samplers exploreon restricted domains

Constraint-handling samplers are designed for targets defined on
restricted domains, such as box constraints \\a \leq \theta \leq b\\,
simplex constraints \\\sum_j \theta_j = 1\\, or targets with non-smooth
penalty terms (L1 regularization, total variation). Standard
gradient-based samplers ignore constraints and either require
reparameterization (which distorts the geometry) or reject proposals
outside the feasible set (which wastes computation). The methods in this
section handle constraints natively, either by projecting the dynamics
onto the feasible set or by using proximal operators that account for
non-smooth penalties.

### Projected Langevin

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | There is no acceptance step; all projected proposals are used. The asymptotic bias from discretization vanishes as the step size decreases. |
| Applications | This algorithm is applicable to targets with box or simplex constraints. It is a natural choice when the parameter space is a convex set and gradient information is available. |
| Difficulty | This algorithm is easy for a beginner. It has four algorithm specifications. The constraint type determines the projection operator. |
| Final Algorithm? | User Discretion. The algorithm introduces discretization bias, so it should be used with small step sizes or as a fast constrained-domain explorer. |
| Proposal | Multivariate. |

Projected Langevin dynamics was studied by [\[86\]](#ref86) as a method
for sampling from log-concave distributions on convex sets. The idea is
simple: run an overdamped Langevin step \\\tilde{\theta} = \theta +
\frac{\epsilon}{2}\nabla \log \pi(\theta) + \sqrt{\epsilon}\\\xi\\
(where \\\xi \sim \mathcal{N}(0, I_d)\\) and then project back onto the
constraint set: \\\theta^\* =
\text{Proj}\_{\mathcal{C}}(\tilde{\theta})\\. For box constraints \\\[a,
b\]^d\\, the projection is componentwise clamping \\\theta_j^\* =
\max(a_j, \min(b_j, \tilde{\theta}\_j))\\. For the simplex \\\\\theta :
\theta_j \geq 0, \sum_j \theta_j = 1\\\\, the projection uses the
algorithm of Duchi, Shalev-Shwartz, Singer, and Chandra (2008), which
sorts the components and finds the threshold for the KKT conditions. No
Metropolis correction is applied, so the output has asymptotic bias of
order \\O(\epsilon)\\, but this is often acceptable for constrained
problems where the alternative is to use a transformation that distorts
the geometry.

Projected Langevin has four algorithm specifications: `epsilon` is the
step size, `constraint` is the constraint type (`"box"` or `"simplex"`,
default `"box"`), `lower` is the vector of lower bounds for box
constraints (default 0), and `upper` is the vector of upper bounds for
box constraints (default \\\infty\\).

#### Algorithm

1.  **Initialize** parameters \\\theta_0 \in \mathcal{C}\\, step size
    \\\epsilon\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute gradient \\g = \nabla \log \pi(\theta\_{t-1})\\
    2.  Langevin step: \\\tilde{\theta} = \theta\_{t-1} +
        \frac{\epsilon}{2} g + \sqrt{\epsilon}\\\xi\\, \\\xi \sim
        \mathcal{N}(0, I_d)\\
    3.  Project: \\\theta_t =
        \text{Proj}\_{\mathcal{C}}(\tilde{\theta})\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Step size of the projected Langevin scheme. |
| constraint | “box” | Box constraint by default; switch to “simplex” for compositional parameters. |
| lower | rep(-Inf, LIV) | No lower bound by default. |
| upper | rep(Inf, LIV) | No upper bound by default. |

Default `Specs` for ProjLang. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "ProjLang", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [MALA](#mala)
- [ProxMCMC](#proxmcmc)
- [Non-Reversible Overdamped Langevin](#nrolangevin)

### ProxMCMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The acceptance rate depends on the step size and the strength of the penalty. For moderate regularization, acceptance rates of 40-70% are typical. |
| Applications | This algorithm is applicable to Bayesian models with non-smooth priors such as L1 (Lasso/Laplace) or total variation penalties. It handles the non-differentiable penalty via proximal operators rather than subgradients. |
| Difficulty | This algorithm is moderately difficult for a beginner. It has three algorithm specifications. Understanding proximal operators is required. |
| Final Algorithm? | User Discretion. ProxMCMC may be used as a final algorithm once the step size and regularization strength are tuned. |
| Proposal | Multivariate. |

ProxMCMC was developed by [\[87\]](#ref87) as a framework for sampling
from posterior distributions with non-smooth penalty terms of the form
\\\pi(\theta) \propto \exp(-f(\theta) - \lambda \cdot g(\theta))\\,
where \\f\\ is the smooth log-likelihood and \\g\\ is a non-smooth
penalty (L1 norm, total variation, etc.). Standard gradient-based MCMC
cannot handle \\g\\ because it is not differentiable everywhere.
ProxMCMC replaces the gradient step with a proximal step: instead of
\\\theta + \epsilon \nabla f(\theta)\\, it uses \\\text{prox}\_{\epsilon
\lambda g}(\theta + \epsilon \nabla f(\theta))\\, where the proximal
operator is \\\text{prox}\_{\tau g}(x) = \arg\min_z \\g(z) +
\frac{1}{2\tau}\\z - x\\^2\\\\. For the L1 penalty, the proximal
operator is the soft-thresholding function \\\text{sign}(x) \max(\|x\| -
\tau\lambda, 0)\\; for total variation, it is the TV denoising operator.
The proximal step handles the non-smoothness exactly, and a
Metropolis-Hastings correction ensures the correct stationary
distribution using the Moreau-Yosida envelope of \\g\\ in the acceptance
ratio.

ProxMCMC has three algorithm specifications: `epsilon` is the step size,
`prox` is the penalty type (`"l1"` or `"tv"`, default `"l1"`), and
`lambda_reg` is the regularization strength (default 0.1).

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, step size \\\epsilon\\,
    penalty strength \\\lambda\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Compute smooth gradient \\g = \nabla f(\theta\_{t-1})\\
    2.  Gradient step: \\\tilde{\theta} = \theta\_{t-1} + \epsilon g +
        \sqrt{2\epsilon}\\\xi\\, \\\xi \sim \mathcal{N}(0, I_d)\\
    3.  Proximal step: \\\theta^\* = \text{prox}\_{\epsilon \lambda
        \\\cdot\\}(\tilde{\theta})\\
    4.  Compute acceptance ratio using the Moreau-Yosida envelope
    5.  Accept with probability \\\min(1, \alpha)\\; otherwise retain
        \\\theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| epsilon | 0.1 / LIV^0.25 (~ 0.076) | Step size. |
| prox | “l1” | Soft-thresholding proximity operator (sparsity prior). |
| lambda_reg | 0.1 | Regularization strength. |

Default `Specs` for ProxMCMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "ProxMCMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [Projected Langevin](#projlang)
- [MALA](#mala)
- [Barker Proposal](#barker)

## Samplers for intractable likelihoods and infinite-dimensional targets

Some models have likelihoods that cannot be evaluated pointwise but can
be estimated unbiasedly through simulation, as in state-space models
with latent processes. Pseudo-marginal MCMC (PMCMC) substitutes the true
likelihood with an unbiased particle-filter estimate and still targets
the exact posterior, at the cost of increased variance in the acceptance
probability. The preconditioned Crank-Nicolson (PCN) algorithm addresses
a different challenge: sampling in function spaces (infinite-dimensional
targets) where standard random-walk proposals degenerate as
discretization refines, by using a proposal that is well-defined in the
continuum limit.

### Pseudo-Marginal MCMC

> *New MCMC algorithm in lucifer*

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This algorithm is designed for models where the likelihood is intractable but can be estimated unbiasedly via simulation, such as state-space models, latent variable models, and models requiring Monte Carlo integration. |
| Difficulty | This algorithm is easy for a beginner. It has one algorithm specification and uses the same proposal covariance machinery as Random-Walk Metropolis. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Pseudo-Marginal MCMC (PMCMC) was introduced by [\[89\]](#ref89) and
formalized by [\[88\]](#ref88), who proved that replacing an intractable
likelihood with an unbiased estimate thereof still targets the exact
posterior distribution. The key insight is that in a standard
Metropolis-Hastings algorithm, if the log-posterior at the current
position is never re-evaluated, and the Model function returns a noisy
but unbiased estimate of the log-likelihood, then the resulting Markov
chain has the correct stationary distribution. This makes PMCMC
applicable to problems where the likelihood cannot be computed
analytically, for instance models with latent variables that are
integrated out via importance sampling or particle filters.

PMCMC has one algorithm specification: `B` accepts a list of blocks and
defaults to a single block containing all parameters. The proposal
mechanism is identical to Random-Walk Metropolis (RWM), using a
multivariate normal proposal centered at the current state with the
proposal covariance matrix. The essential difference from RWM is that
PMCMC never re-evaluates the Model at the current parameter values once
a proposal is accepted; it stores the noisy log-posterior estimate from
the acceptance step and reuses it for all subsequent comparisons. This
property is critical: re-evaluating the noisy likelihood at the same
point would introduce bias, because the acceptance ratio would no longer
be an unbiased estimate of the true Metropolis-Hastings ratio.

In practice, the user must ensure that the Model function returns an
unbiased estimate of the log-likelihood on the log scale. Common
strategies include importance sampling with auxiliary variables,
sequential Monte Carlo (particle filtering) for state-space models, and
Monte Carlo integration for marginalizing latent variables. The variance
of the log-likelihood estimate affects mixing: [\[90\]](#ref90)
recommend targeting a variance of approximately 1-3 on the log scale for
the likelihood estimate, as higher variance degrades acceptance rates
and lower variance wastes computational effort on unnecessary precision.

The advantage of PMCMC is that it enables exact Bayesian inference (in
the sense of targeting the correct posterior) for models with
intractable likelihoods, requiring only that the user can construct an
unbiased likelihood estimator. The disadvantage is that noisy likelihood
estimates reduce the effective acceptance rate compared to RWM with an
exact likelihood, and the variance of the estimator must be carefully
controlled. Since PMCMC is not adaptive, it is suitable as a final
algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\; evaluate and store noisy
    log-posterior \\\hat{\ell}\_0 = \log \hat{p}(\theta_0 \mid y)\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Propose \\\theta^\* \sim \mathcal{N}(\theta\_{t-1}, \Sigma)\\
        (as in RWM)
    2.  Evaluate model **once** at \\\theta^\*\\: obtain noisy estimate
        \\\hat{\ell}^\* = \log \hat{p}(\theta^\* \mid y)\\
    3.  Compute \\\alpha = \min\\\bigl(1,\\ \exp(\hat{\ell}^\* -
        \hat{\ell}\_{t-1})\bigr)\\
    4.  With probability \\\alpha\\, accept: \\\theta_t = \theta^\*\\,
        \\\hat{\ell}\_t = \hat{\ell}^\*\\
    5.  Otherwise, retain: \\\theta_t = \theta\_{t-1}\\, \\\hat{\ell}\_t
        = \hat{\ell}\_{t-1}\\ (**do not re-evaluate** at
        \\\theta\_{t-1}\\)
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification                             |
|:--------------|:--------|:------------------------------------------|
| B             | list()  | Single block (no parameter partitioning). |

Default `Specs` for PMCMC. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "PMCMC", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [RWM](#rwm)
- [SMC (via `SMC()`)](#rwm)

### Preconditioned Crank-Nicolson

| Aspect | Description |
|:---|:---|
| Acceptance Rate | The optimal acceptance rate is based on the multivariate normality of the marginal posterior distributions, and ranges from 44% for one parameter to 23.4% for five or more parameters. The observed acceptance rate may be suitable in the interval \[15%,50%\]. |
| Applications | This is a widely applicable, general-purpose algorithm that is most useful with large-dimensional models. |
| Difficulty | This algorithm is easy to use, having only one tuning parameter. |
| Final Algorithm? | Yes. |
| Proposal | Multivariate. |

Preconditioned Crank-Nicolson (pCN) was introduced originally as the PIA
algorithm in [\[91\]](#ref91), and differs only slightly from
Random-Walk Metropolis (RWM). The proposal is a first-order
autoregressive process, rather than a centered random-walk. The pCN
algorithm has an acceptance probability that is invariant to dimension,
making pCN useful for large-dimensional models.

pCN has one algorithm specification, `beta`, which is the autoregressive
parameter in (0,1).

As with RWM, the target acceptance rate is 23.4%. pCN seems to perform
well only with an identity matrix as the proposal covariance matrix.

Since pCN is not adaptive, it is suitable as a final algorithm.

#### Algorithm

1.  **Initialize** parameters \\\theta_0\\, autoregressive parameter
    \\\beta \in (0,1)\\
2.  **For** \\t = 1, \ldots, T\\:
    1.  Draw \\\xi \sim \mathcal{N}(0, \Sigma)\\
    2.  Propose \\\theta^\* = \sqrt{1 - \beta^2}\\ \theta\_{t-1} +
        \beta\\ \xi\\ (first-order autoregressive)
    3.  Compute \\\alpha = \min\\\bigl(1,\\ p(\theta^\* \mid y) /
        p(\theta\_{t-1} \mid y)\bigr)\\
    4.  With probability \\\alpha\\, accept \\\theta_t = \theta^\*\\;
        otherwise \\\theta_t = \theta\_{t-1}\\
3.  **Return** \\\\\theta_1, \ldots, \theta_T\\\\

#### Default specifications

| Specification | Default | Justification |
|:---|:---|:---|
| beta | 0.01 | Autoregressive parameter; small values keep proposals close to the current state, which is essential for high- or infinite-dimensional targets. |

Default `Specs` for pCN. {.table}

#### Example with default specifications

``` r

fit <- lucifer(Model, Data, IV,
               Iterations = 2000, Status = 1000, Thinning = 2,
               Algorithm = "pCN", Specs = NULL,
               Chains = 3)
plot(fit)
ppc <- predict(fit, Model = Model, Data = Data)
plot(ppc, Style = "Fitted", Data = Data)
ppc_dens_overlay(ppc)
plot(fit, ground_truth = ground_truth)
```

#### See Also

- [RAM](#ram)

## References

**\[1\]** Gelman, A., Roberts, G., and Gilks, W. (1996). *Efficient
Metropolis Jumping Rules.* Bayesian Statistics, 5, 599–608.

**\[2\]** Haario, H., Saksman, E., and Tamminen, J. (2001). *An Adaptive
Metropolis Algorithm.* Bernoulli, 7(2), 223–242.
[doi:10.3150/bj/1080222083](https://doi.org/10.3150/bj/1080222083)

**\[3\]** Roberts, G. and Rosenthal, J. (2009). *Examples of Adaptive
MCMC.* Journal of Computational and Graphical Statistics, 18(2),
349–367.
[doi:10.1198/jcgs.2009.06134](https://doi.org/10.1198/jcgs.2009.06134)

**\[4\]** Vihola, M. (2012). *Robust Adaptive Metropolis Algorithm with
Coerced Acceptance Rate.* Statistics and Computing, 22(5), 997–1008.
[doi:10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5)

**\[5\]** Haario, H., Laine, M., Mira, A., and Saksman, E. (2006).
*DRAM: Efficient Adaptive MCMC.* Statistics and Computing, 16, 339–354.
[doi:10.1007/s11222-006-9438-0](https://doi.org/10.1007/s11222-006-9438-0)

**\[6\]** Mira, A. (2001). *On Metropolis-Hastings Algorithms with
Delayed Rejection.* Metron, LIX(3–4), 231–241.

**\[7\]** Christen, J.A. and Fox, C. (2005). *Markov Chain Monte Carlo
Using an Approximation.* Journal of Computational and Graphical
Statistics, 14(4), 795–810.
[doi:10.1198/106186005X76983](https://doi.org/10.1198/106186005X76983)

**\[8\]** Sherlock, C., Golightly, A., and Sherlock, C. (2017).
*Adaptive, Delayed-Acceptance MCMC for Targets with Expensive
Likelihoods.* Journal of Computational and Graphical Statistics, 26(2),
434–444.
[doi:10.1080/10618600.2016.1231064](https://doi.org/10.1080/10618600.2016.1231064)

**\[9\]** Hastings, W.K. (1970). *Monte Carlo Sampling Methods Using
Markov Chains and Their Applications.* Biometrika, 57(1), 97–109.
[doi:10.1093/biomet/57.1.97](https://doi.org/10.1093/biomet/57.1.97)

**\[10\]** Tierney, L. (1994). *Markov Chains for Exploring Posterior
Distributions.* The Annals of Statistics, 22(4), 1701–1762.
[doi:10.1214/aos/1176325750](https://doi.org/10.1214/aos/1176325750)

**\[11\]** Liu, J., Liang, F., and Wong, W. (2000). *The Multiple-Try
Method and Local Optimization in Metropolis Sampling.* Journal of the
American Statistical Association, 95, 121–134.
[doi:10.1080/01621459.2000.10473908](https://doi.org/10.1080/01621459.2000.10473908)

**\[12\]** Craiu, R., Rosenthal, J., and Yang, C. (2009). *Learn From
Thy Neighbor: Parallel-Chain and Regional Adaptive MCMC.* Journal of the
American Statistical Association, 104(488), 1454–1466.
[doi:10.1198/jasa.2009.tm08393](https://doi.org/10.1198/jasa.2009.tm08393)

**\[13\]** Solonen, A., Ollinaho, P., Laine, M., Haario, H., Tamminen,
J., and Jarvinen, H. (2012). *Efficient MCMC for Climate Model Parameter
Estimation: Parallel Adaptive Chains and Early Rejection.* Bayesian
Analysis, 7(2), 1–22.
[doi:10.1214/12-BA714](https://doi.org/10.1214/12-BA714)

**\[14\]** Geman, S. and Geman, D. (1984). *Stochastic Relaxation, Gibbs
Distributions, and the Bayesian Restoration of Images.* IEEE
Transactions on Pattern Analysis and Machine Intelligence, 6(6),
721–741.
[doi:10.1109/TPAMI.1984.4767596](https://doi.org/10.1109/TPAMI.1984.4767596)

**\[15\]** Neal, R.M. (2003). *Slice Sampling.* Annals of Statistics,
31(3), 705–767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461)

**\[16\]** Turchin, V.F. (1971). *On the Computation of Multidimensional
Integrals by the Monte Carlo Method.* Theory of Probability and its
Applications, 16(4), 720–724.
[doi:10.1137/1116083](https://doi.org/10.1137/1116083)

**\[17\]** Bai, Y. (2009). *An Adaptive Directional
Metropolis-within-Gibbs Algorithm.* Technical Report, Department of
Statistics, University of Toronto.

**\[18\]** Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., Teller,
A.H., and Teller, E. (1953). *Equation of State Calculations by Fast
Computing Machines.* Journal of Chemical Physics, 21, 1087–1092.
[doi:10.1063/1.1699114](https://doi.org/10.1063/1.1699114)

**\[19\]** Roberts, G. and Rosenthal, J. (2007). *Coupling and
Ergodicity of Adaptive Markov Chain Monte Carlo Algorithms.* Journal of
Applied Probability, 44, 458–475.
[doi:10.1239/jap/1183667414](https://doi.org/10.1239/jap/1183667414)

**\[20\]** Ritter, C. and Tanner, M. (1992). *Facilitating the Gibbs
Sampler: the Gibbs Stopper and the Griddy-Gibbs Sampler.* Journal of the
American Statistical Association, 87, 861–868.
[doi:10.1080/01621459.1992.10475289](https://doi.org/10.1080/01621459.1992.10475289)

**\[21\]** Geweke, J. and Tanizaki, H. (2001). *Bayesian Estimation of
State-Space Models Using the Metropolis-Hastings Algorithm within Gibbs
Sampling.* Computational Statistics and Data Analysis, 37, 151–170.
[doi:10.1016/S0167-9473(01)00045-9](https://doi.org/10.1016/S0167-9473(01)00045-9)

**\[22\]** Fearnhead, P. (2011). *MCMC for State-Space Models.* In S
Brooks, A Gelman, G Jones, M Xiao-Li (eds.), Handbook of Markov Chain
Monte Carlo, p. 513-530. Chapman & Hall, Boca Raton, FL. ISBN:
978-1-4200-7941-8.

**\[23\]** Dutta, S. (2012). *Multiplicative Random Walk
Metropolis-Hastings on the Real Line.* Sankhya B, 74(2), 315–342.
[doi:10.1007/s13571-012-0040-5](https://doi.org/10.1007/s13571-012-0040-5)

**\[24\]** Zanella, G. (2020). *Informed Proposals for Local MCMC in
Discrete Spaces.* Journal of the American Statistical Association,
115(530), 852–865.
[doi:10.1080/01621459.2019.1585255](https://doi.org/10.1080/01621459.2019.1585255)

**\[25\]** Power, S. and Goldman, J.V. (2019). *Accelerated Sampling on
Discrete Spaces with Non-Reversible Markov Chains.* arXiv preprint,
arXiv:1912.04681.

**\[26\]** Polson, N.G., Scott, J.G., and Windle, J. (2013). *Bayesian
Inference for Logistic Models Using Polya-Gamma Latent Variables.*
Journal of the American Statistical Association, 108(504), 1339–1349.
[doi:10.1080/01621459.2013.829001](https://doi.org/10.1080/01621459.2013.829001)

**\[27\]** Windle, J., Polson, N.G., and Scott, J.G. (2014). *Sampling
Polya-Gamma Random Variates: Alternate and Approximate Techniques.*
arXiv preprint, arXiv:1405.0506.

**\[28\]** Duane, S., Kennedy, A.D., Pendleton, B.J., and Roweth, D.
(1987). *Hybrid Monte Carlo.* Physics Letters B, 195(2), 216–222.
[doi:10.1016/0370-2693(87)91197-X](https://doi.org/10.1016/0370-2693(87)91197-X)

**\[29\]** Neal, R. (2011). *MCMC for Using Hamiltonian Dynamics.* In S
Brooks, A Gelman, G Jones, M Xiao-Li (eds.), Handbook of Markov Chain
Monte Carlo, p. 113-162. Chapman & Hall, Boca Raton, FL. ISBN:
978-1-4200-7941-8.

**\[30\]** Rossky, P., Doll, J., and Friedman, H. (1978). *Brownian
Dynamics as Smart Monte Carlo Discrete Approximations.* Journal of
Chemical Physics, 69, 4628–4633.
[doi:10.1063/1.436415](https://doi.org/10.1063/1.436415)

**\[31\]** Hoffman, M. and Gelman, A. (2014). *The No-U-Turn Sampler:
Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.* Journal of
Machine Learning Research, 15, 1593–1623.

**\[32\]** Robnik, J., De Luca, G.B., Silverstein, E., and Seljak, U.
(2023). *Microcanonical Hamiltonian Monte Carlo.* Journal of Machine
Learning Research, 24(311), 1–34.

**\[33\]** Hoffman, M.D. and Sountsov, P. (2022). *Tuning-Free
Generalized Hamiltonian Monte Carlo.* Proceedings of the 25th
International Conference on Artificial Intelligence and Statistics
(AISTATS), PMLR 151, 7799–7813.

**\[34\]** Nesterov, Y. (2009). *Primal-Dual Subgradient Methods for
Convex Problems.* Mathematical Programming, 120, 221–259.
[doi:10.1007/s10107-007-0149-x](https://doi.org/10.1007/s10107-007-0149-x)

**\[35\]** Biron-Lattes, M., Surjanovic, N., Syed, S., Campbell, T., and
Bouchard-Cote, A. (2024). *autoMALA: Locally Adaptive
Metropolis-Adjusted Langevin Algorithm.* Proceedings of the 27th
International Conference on Artificial Intelligence and Statistics
(AISTATS), PMLR 238, 4600–4608.

**\[36\]** Riou-Durand, L. and Vogrinc, J. (2022). *Metropolis Adjusted
Langevin Trajectories: A Robust Alternative to Hamiltonian Monte Carlo.*
arXiv preprint arXiv:2202.13230.

**\[37\]** Sherlock, C., Sherlock, C., and Sherlock, C. (2023). *Apogee
to Apogee Path Sampler.* arXiv preprint arXiv:2302.10629.

**\[38\]** Bou-Rabee, N., Marsden, A., and Saibaba, K. (2024). *GIST:
Gibbs Self-Tuning for Locally Adaptive Hamiltonian Monte Carlo.* arXiv
preprint arXiv:2404.15253.

**\[39\]** Atchade, Y. (2006). *An Adaptive Version for the Metropolis
Adjusted Langevin Algorithm with a Truncated Drift.* Methodology and
Computing in Applied Probability, 8, 235–254.
[doi:10.1007/s11009-006-9550-0](https://doi.org/10.1007/s11009-006-9550-0)

**\[40\]** Roberts, G. and Tweedie, R. (1996). *Exponential Convergence
of Langevin Distributions and Their Discrete Approximations.* Bernoulli,
2(4), 341–363. [doi:10.2307/3318418](https://doi.org/10.2307/3318418)

**\[41\]** Shaby, B. and Wells, M. (2010). *Exploring an Adaptive
Metropolis Algorithm.* Working Paper, Department of Statistical Science,
Duke University.

**\[42\]** Welling, M. and Teh, Y. (2011). *Bayesian Learning via
Stochastic Gradient Langevin Dynamics.* Proceedings of the 28th
International Conference on Machine Learning (ICML), 681–688.

**\[43\]** Chen, T., Fox, E., and Guestrin, C. (2014). *Stochastic
Gradient Hamiltonian Monte Carlo.* Proceedings of the 31st International
Conference on Machine Learning (ICML), 1683–1691.

**\[44\]** Livingstone, S. and Zanella, G. (2022). *The Barker Proposal:
Combining Robustness and Efficiency in Gradient-Based MCMC.* Journal of
the Royal Statistical Society: Series B, 84(2), 496–523.
[doi:10.1111/rssb.12482](https://doi.org/10.1111/rssb.12482)

**\[45\]** Duncan, A.B., Lelievre, T., and Pavliotis, G.A. (2016).
*Variance Reduction Using Nonreversible Langevin Samplers.* Journal of
Statistical Physics, 163(3), 457–491.
[doi:10.1007/s10955-016-1491-2](https://doi.org/10.1007/s10955-016-1491-2)

**\[46\]** Besag, J. and Green, P.J. (1993). *Spatial Statistics and
Bayesian Computation.* Journal of the Royal Statistical Society, Series
B, 55, 25–37.
[doi:10.1111/j.2517-6161.1993.tb01467.x](https://doi.org/10.1111/j.2517-6161.1993.tb01467.x)

**\[47\]** Neal, R. (1997). *Markov Chain Monte Carlo Methods Based on
Slicing the Density Function.* Technical Report, University of Toronto.

**\[48\]** Tibbits, M., Groendyke, C., Haran, M., and Liechty, J.
(2014). *Automated Factor Slice Sampling.* Journal of Computational and
Graphical Statistics, 23(2), 543–563.
[doi:10.1080/10618600.2013.791193](https://doi.org/10.1080/10618600.2013.791193)

**\[49\]** Murray, I., Adams, R., and MacKay, D. (2010). *Elliptical
Slice Sampling.* Proceedings of the 13th International Conference on
Artificial Intelligence and Statistics (AISTATS), JMLR W&CP 9, 541–548.

**\[50\]** Thompson, M.D. (2011). *Slice Sampling with Multivariate
Steps.* Ph.D. thesis, University of Toronto.
<http://hdl.handle.net/1807/31955>.

**\[51\]** Heiner, M.J., Johnson, S.B., Christensen, J.R., and Dorie, D.
(2024). *Quantile Slice Sampling.* Statistics and Computing, 34, 36.
[doi:10.1007/s11222-023-10338-x](https://doi.org/10.1007/s11222-023-10338-x)

**\[52\]** Li, Y. and Walker, S.G. (2023). *A Latent Slice Sampling
Algorithm.* Computational Statistics and Data Analysis, 179, 107652.
[doi:10.1016/j.csda.2022.107652](https://doi.org/10.1016/j.csda.2022.107652)

**\[53\]** Schar, P., Falck, F., Thomas, R.P., and Ge, H. (2023).
*Gibbsian Polar Slice Sampling.* Proceedings of the 40th International
Conference on Machine Learning (ICML), PMLR 202, 30204–30223.

**\[54\]** Foreman-Mackey, D., Hogg, D., Lang, D., and Goodman, J.
(2013). *emcee: The MCMC Hammer.* Publications of the Astronomical
Society of the Pacific, 125(925), 306–312.
[doi:10.1086/670067](https://doi.org/10.1086/670067)

**\[55\]** Goodman, J. and Weare, J. (2010). *Ensemble Samplers with
Affine Invariance.* Communications in Applied Mathematics and
Computational Science, 5(1), 65–80.
[doi:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)

**\[56\]** Ter Braak, C. (2006). *A Markov Chain Monte Carlo Version of
the Genetic Algorithm Differential Evolution: Easy Bayesian Computing
for Real Parameter Spaces.* Statistics and Computing, 16, 239–249.
[doi:10.1007/s11222-006-8769-1](https://doi.org/10.1007/s11222-006-8769-1)

**\[57\]** Ter Braak, C. and Vrugt, J. (2008). *Differential Evolution
Markov Chain with Snooker Updater and Fewer Chains.* Statistics and
Computing, 18(4), 435–446.
[doi:10.1007/s11222-008-9104-9](https://doi.org/10.1007/s11222-008-9104-9)

**\[58\]** Christen, J.A. and Fox, C. (2010). *A General Purpose
Sampling Algorithm for Continuous Distributions (the t-walk).* Bayesian
Analysis, 5(2), 263–282.
[doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

**\[59\]** Militzer, B. (2023). *Quadratic Monte Carlo.* Physical Review
E, 108(5), 055309.
[doi:10.1103/PhysRevE.108.055309](https://doi.org/10.1103/PhysRevE.108.055309)

**\[60\]** Militzer, B. (2025). *Five Novel Move Types for Ensemble
Monte Carlo.* arXiv preprint, arXiv:2501.02382.

**\[61\]** Karamanis, M. and Beutler, F. (2021). *Ensemble Slice
Sampling: Beyond Linear Inefficiency.* Statistics and Computing, 31, 61.
[doi:10.1007/s11222-021-10038-2](https://doi.org/10.1007/s11222-021-10038-2)

**\[62\]** Vrugt, J.A., ter Braak, C.J.F., Diks, C.G.H., Robinson, B.A.,
Hyman, J.M., and Higdon, D. (2009). *Accelerating Markov Chain Monte
Carlo Simulation by Differential Evolution with Self-Adaptive Randomized
Subspace Sampling.* International Journal of Nonlinear Sciences and
Numerical Simulation, 10(3), 273–290.
[doi:10.1515/IJNSNS.2009.10.3.273](https://doi.org/10.1515/IJNSNS.2009.10.3.273)

**\[63\]** Chen, M. and Schmeiser, B. (1993). *Performance of the Gibbs,
Hit-And-Run and Metropolis Samplers.* Journal of Computational and
Graphical Statistics, 2(3), 251–272.
[doi:10.1080/10618600.1993.10474611](https://doi.org/10.1080/10618600.1993.10474611)

**\[64\]** Garthwaite, P., Fan, Y., and Sisson, S. (2010). *Adaptive
Optimal Scaling of Metropolis-Hastings Algorithms Using the
Robbins-Monro Process.* arXiv preprint, arXiv:1006.3690.

**\[65\]** Gilks, W. and Roberts, G. (1996). *Strategies for Improving
MCMC.* In W Gilks, S Richardson, D Spiegelhalter (eds.), Markov Chain
Monte Carlo in Practice, p. 89-114. Chapman & Hall, Boca Raton, FL.
ISBN: 978-0-412-05551-5.

**\[66\]** Smith, R. (1984). *Efficient Monte Carlo Procedures for
Generating Points Uniformly Distributed Over Bounded Region.* Operations
Research, 32, 1296–1308.
[doi:10.1287/opre.32.6.1296](https://doi.org/10.1287/opre.32.6.1296)

**\[67\]** Boyles, L.B. and Welling, M. (2012). *Refractive Sampling.*
Technical Report, UC Irvine.

**\[68\]** Atchade, Y.F., Roberts, G.O., and Rosenthal, J.S. (2011).
*Towards Optimal Scaling of Metropolis-Coupled Markov Chain Monte
Carlo.* Statistics and Computing, 21(4), 555–568.
[doi:10.1007/s11222-010-9192-1](https://doi.org/10.1007/s11222-010-9192-1)

**\[69\]** Earl, D.J. and Deem, M.W. (2005). *Parallel Tempering:
Theory, Applications, and New Perspectives.* Physical Chemistry Chemical
Physics, 7, 3910–3916.
[doi:10.1039/b509983h](https://doi.org/10.1039/b509983h)

**\[70\]** Frantz, D.D., Freeman, J.D., and Doll, J.L. (1990). *Reducing
Quasi-Ergodic Behavior in Monte Carlo Simulations by J-Walking:
Applications to Atomic Clusters.* Journal of Chemical Physics, 93,
2769–2784. [doi:10.1063/1.458513](https://doi.org/10.1063/1.458513)

**\[71\]** Geyer, C.J. (1991). *Markov Chain Monte Carlo Maximum
Likelihood.* In Keramidas EM (ed.), Computing Science and Statistics:
Proceedings of the 23rd Symposium on the Interface, p. 156–163.
Interface Foundation, Fairfax Station, VA.

**\[72\]** Green, P. (1995). *Reversible Jump Markov Chain Monte Carlo
Computation and Bayesian Model Determination.* Biometrika, 82(4),
711–732.
[doi:10.1093/biomet/82.4.711](https://doi.org/10.1093/biomet/82.4.711)

**\[73\]** Syed, S., Bouchard-Cote, A., Deligiannidis, G., and Doucet,
A. (2022). *Non-Reversible Parallel Tempering: A Scalable Highly
Parallel MCMC Scheme.* Journal of the Royal Statistical Society: Series
B, 84(2), 321–350.
[doi:10.1111/rssb.12464](https://doi.org/10.1111/rssb.12464)

**\[74\]** Miasojedow, B., Moulines, E., and Vihola, M. (2013). *An
Adaptive Parallel Tempering Algorithm.* Journal of Computational and
Graphical Statistics, 22(3), 649–664.
[doi:10.1080/10618600.2013.778779](https://doi.org/10.1080/10618600.2013.778779)

**\[75\]** Marinari, E. and Parisi, G. (1992). *Simulated Tempering: A
New Monte Carlo Scheme.* Europhysics Letters, 19(6), 451–458.
[doi:10.1209/0295-5075/19/6/002](https://doi.org/10.1209/0295-5075/19/6/002)

**\[76\]** Wang, F. and Landau, D.P. (2001). *Efficient, Multiple-Range
Random Walk Algorithm to Calculate the Density of States.* Physical
Review Letters, 86(10), 2050–2053.
[doi:10.1103/PhysRevLett.86.2050](https://doi.org/10.1103/PhysRevLett.86.2050)

**\[77\]** Bierkens, J., Fearnhead, P., and Roberts, G. (2019). *The
Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of
Big Data.* Annals of Statistics, 47(3), 1288–1320.
[doi:10.1214/18-AOS1715](https://doi.org/10.1214/18-AOS1715)

**\[78\]** Bouchard-Cote, A., Vollmer, S.J., and Doucet, A. (2018). *The
Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain
Monte Carlo Method.* Journal of the American Statistical Association,
113(522), 855–867.
[doi:10.1080/01621459.2017.1294075](https://doi.org/10.1080/01621459.2017.1294075)

**\[79\]** Bierkens, J., Kamatani, K., and Roberts, G.O. (2020). *The
Boomerang Sampler.* Proceedings of the 37th International Conference on
Machine Learning (ICML), PMLR 119, 908–918.

**\[80\]** Bou-Rabee, N. and Sanz-Serna, J.M. (2017). *Randomized
Hamiltonian Monte Carlo.* Annals of Applied Probability, 27(4),
2159–2194. [doi:10.1214/16-AAP1255](https://doi.org/10.1214/16-AAP1255)

**\[81\]** Girolami, M. and Calderhead, B. (2011). *Riemann Manifold
Langevin and Hamiltonian Monte Carlo Methods.* Journal of the Royal
Statistical Society: Series B, 73(2), 123–214.
[doi:10.1111/j.1467-9868.2010.00765.x](https://doi.org/10.1111/j.1467-9868.2010.00765.x)

**\[82\]** Betancourt, M. (2013). *A General Metric for Riemannian
Manifold Hamiltonian Monte Carlo.* Geometric Science of Information,
LNCS 8085, 327–334.
[doi:10.1007/978-3-642-40020-9_35](https://doi.org/10.1007/978-3-642-40020-9_35)

**\[83\]** Lan, S., Stathopoulos, V., Mark, B., and Girolami, M. (2015).
*Markov Chain Monte Carlo from Lagrangian Dynamics.* Journal of
Computational and Graphical Statistics, 24(2), 357–378.
[doi:10.1080/10618600.2014.902764](https://doi.org/10.1080/10618600.2014.902764)

**\[84\]** Tripuraneni, N., Rowland, M., Ghahramani, Z., and Turner, R.
(2017). *Magnetic Hamiltonian Monte Carlo.* Proceedings of the 34th
International Conference on Machine Learning (ICML), PMLR 70, 3453–3461.

**\[85\]** Lu, X., Perrone, V., Hasenclever, L., Teh, Y.W., and Vollmer,
S.J. (2017). *Relativistic Monte Carlo.* Proceedings of the 20th
International Conference on Artificial Intelligence and Statistics
(AISTATS), PMLR 54, 1236–1245.

**\[86\]** Bubeck, S., Eldan, R., and Lehec, J. (2018). *Sampling from a
Log-Concave Distribution with Projected Langevin Monte Carlo.* Discrete
and Computational Geometry, 59(4), 757–783.
[doi:10.1007/s00454-018-9992-1](https://doi.org/10.1007/s00454-018-9992-1)

**\[87\]** Zhou, Y., Chen, T., Paisley, J., and Lu, J. (2024). *Proximal
MCMC for Bayesian Inference of Constrained and Regularized Estimation.*
arXiv preprint arXiv:2205.07378.

**\[88\]** Andrieu, C. and Roberts, G.O. (2009). *The Pseudo-Marginal
Approach for Efficient Monte Carlo Computations.* Annals of Statistics,
37(2), 697–725.
[doi:10.1214/07-AOS574](https://doi.org/10.1214/07-AOS574)

**\[89\]** Beaumont, M.A. (2003). *Estimation of Population Growth or
Decline in Genetically Monitored Populations.* Genetics, 164(3),
1139–1160.
[doi:10.1534/genetics.103.012914](https://doi.org/10.1534/genetics.103.012914)

**\[90\]** Doucet, A., Pitt, M.K., Deligiannidis, G., and Kohn, R.
(2015). *Efficient Implementation of Markov Chain Monte Carlo When Using
an Unbiased Likelihood Estimator.* Biometrika, 102(2), 295–313.
[doi:10.1093/biomet/asu075](https://doi.org/10.1093/biomet/asu075)

**\[91\]** Beskos, A., Roberts, G.O., Stuart, A.M., and Voss, J. (2008).
*MCMC Methods for Diffusion Bridges.* Stochastics and Dynamics, 8(3),
319–350.
[doi:10.1142/S0219493708002378](https://doi.org/10.1142/S0219493708002378)
