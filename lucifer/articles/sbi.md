# Simulation-based inference

## Introduction

Many scientific models are defined not by closed-form probability
distributions but by stochastic simulators that encode complex
generative processes. In ecology, an individual-based model might
specify how organisms are born, disperse, compete for resources, and die
through a sequence of random events, yet computing the probability of an
observed population time series under this model is analytically
intractable. The same fundamental problem arises across disciplines: in
cosmology, where the relationship between cosmological parameters and
the cosmic microwave background involves solving Boltzmann equations
with stochastic initial conditions; in epidemiology, where disease
transmission models include contact networks and superspreading events
that preclude analytic likelihood evaluation; and in neuroscience, where
mechanistic models of neural circuits produce spike trains through
stochastic differential equations whose transition densities have no
closed form.

To make this concrete, consider a simple example from population
ecology. The Ricker model describes population dynamics through
\\N\_{t+1} = r N_t \exp(-N_t/K + \epsilon_t)\\, where \\r\\ is the
growth rate, \\K\\ is the carrying capacity, and \\\epsilon_t\\ is
stochastic process noise. Given an observed population time series, one
wants to infer \\r\\ and \\K\\. The likelihood \\p(N_1, \ldots, N_T \mid
r, K)\\ requires marginalizing over all possible realizations of the
noise sequence \\\epsilon_1, \ldots, \epsilon_T\\, which involves a
\\T\\-dimensional integral with no closed form. Yet generating a
simulated time series from the model for any \\(r, K)\\ is trivial:
simply draw the noise terms and iterate the recursion. This asymmetry,
easy to simulate but impossible to evaluate, is the hallmark of
simulator-defined models and the fundamental motivation for
simulation-based inference.

The scale and diversity of these intractable-likelihood problems has
grown enormously in recent years. In genetics, coalescent models with
recombination involve complex genealogical histories that preclude
likelihood computation. In climate science, general circulation models
are multi-scale simulators with millions of interacting components. In
gravitational wave astronomy, waveform templates for binary mergers
depend nonlinearly on 15+ parameters and require expensive numerical
relativity calculations. In all these cases, the simulator embodies
decades of domain expertise and cannot be replaced by a simpler,
analytically tractable model without sacrificing scientific fidelity.
The need for inference methods that work directly with simulators has
driven the development of two complementary approaches: Approximate
Bayesian Computation (ABC) and simulation-based inference (SBI).

Classical Bayesian methods require evaluating the log-likelihood \\\log
p(x \mid \theta)\\ at every iteration, whether one uses Markov chain
Monte Carlo (MCMC), variational inference, or sequential Monte Carlo.
When this quantity is unavailable, these methods simply cannot be
applied directly. Approximate Bayesian Computation (ABC) addressed this
gap by proposing parameters from the prior, simulating synthetic data,
and accepting proposals whose simulations fall within a tolerance ball
\\\rho(x\_{\text{sim}}, x\_{\text{obs}}) \< \epsilon\\ around the
observed data [\[16\]](#ref16). While ABC has been enormously
influential, it suffers from the curse of dimensionality: as the
dimension of the data or summary statistics grows, the fraction of
accepted simulations vanishes, and the tolerance \\\epsilon\\ must be
kept small enough to ensure posterior accuracy yet large enough to
retain a reasonable acceptance rate [\[17\]](#ref17).

Simulation-based inference (SBI), sometimes called likelihood-free
inference or implicit likelihood inference, takes a fundamentally
different approach [\[1\]](#ref1) [\[4\]](#ref4). Rather than discarding
simulations that fall outside a tolerance ball, SBI uses all simulated
pairs \\(\theta_i, x_i)\\ as supervised training data for a neural
density estimator that learns to approximate the posterior distribution
\\p(\theta \mid x)\\, the likelihood function \\p(x \mid \theta)\\, or
the likelihood-to-evidence ratio \\p(x \mid \theta) / p(x)\\ directly.
The result is a full probability distribution, not a kernel-smoothed
approximation, and for amortized methods, posterior inference for new
observations requires only a single forward pass through the trained
network.

Three families of SBI methods have emerged as the dominant paradigms.
Neural Posterior Estimation (NPE) trains a conditional density network
\\q\_\phi(\theta \mid x)\\ to approximate the posterior directly,
yielding amortized inference where the computational cost of training is
incurred once and posterior evaluation for any observation is
essentially free [\[2\]](#ref2). Neural Likelihood Estimation (NLE)
trains a surrogate likelihood \\q\_\phi(x \mid \theta)\\ that can be
combined with any prior via MCMC [\[3\]](#ref3) [\[11\]](#ref11). Neural
Ratio Estimation (NRE) formulates the problem as binary classification
and recovers the likelihood-to-evidence ratio, which also plugs into
MCMC [\[5\]](#ref5).

lucifer’s SBI implementation provides all three families, their
sequential variants (SNPE, SNLE, SNRE), a zero-dependency C++ Mixture
Density Network (MDN) backend with SELU activations and
Cholesky-parametrized covariances, an optional torch backend for Neural
Spline Flows [\[8\]](#ref8), and a comprehensive diagnostic suite (SBC,
expected coverage, C2ST, TARP, posterior predictive checks). The
distinguishing architectural feature, however, is the NLE/NRE-to-MCMC
bridge: the trained neural likelihood or ratio estimator is wrapped as a
standard `Model(parm, Data)` closure that plugs directly into
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
giving the user immediate access to all 62 MCMC algorithms. No
retraining, no adapter code, no Python dependencies. This
interoperability is unique among SBI implementations and means that the
entire diagnostic and post-processing infrastructure of lucifer
(convergence diagnostics, effective sample size, MCSE, LOO-CV) is
available for simulation-based inference out of the box.

By comparison, the sbi Python package [\[23\]](#ref23) provides a
comprehensive PyTorch-based implementation with excellent normalizing
flow support, extensive documentation, and a large user community, but
it requires a Python environment with PyTorch installed and offers no
built-in bridge to classical MCMC samplers beyond a lightweight MCMC
wrapper. BayesFlow [\[22\]](#ref22) focuses on amortized NPE with
invertible neural networks and specializes in repeated inference across
datasets, with an emphasis on summary network learning, but does not
provide NLE or NRE methods. Other notable implementations include the
Julia package ELFI (Engine for Likelihood-Free Inference), which
provides ABC and BOLFI (Bayesian optimization for likelihood-free
inference), and the Python package lampe (Likelihood-free AMortized
Posterior Estimation), which focuses on normalizing flows.

lucifer occupies a complementary position in this landscape. Its core
advantages are threefold: first, a self-contained C++ MDN backend that
requires no external neural network libraries, making installation
trivial on any system with a C++ compiler; second, the NLE/NRE-to-MCMC
bridge that gives access to 62 MCMC algorithms through a standard
`Model()` closure, enabling the full power of classical MCMC theory
(convergence diagnostics, effective sample sizes, parallel chains) for
posterior sampling under a neural surrogate; and third, seamless
connection to the broader lucifer ecosystem, including cross-validation
([`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md),
[`Kfold()`](https://robustecologies.github.io/lucifer/reference/Kfold.md)),
model comparison
([`loo_compare()`](https://robustecologies.github.io/lucifer/reference/loo_compare.md)),
robust Bayesian analysis
([`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)),
and the
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
algorithm benchmarking framework. For R users working on scientific
problems with simulator-defined models, this integrated approach
eliminates the friction of interfacing between separate packages or
languages.

## Mathematical foundations

### The SBI problem

The central goal of Bayesian inference is to compute the posterior
distribution \\p(\theta \mid x_o)\\ given observed data \\x_o\\ and a
model with parameters \\\theta\\. By Bayes’ theorem,

\\p(\theta \mid x_o) = \frac{p(x_o \mid \theta) \\ p(\theta)}{p(x_o)},\\

where \\p(x_o \mid \theta)\\ is the likelihood, \\p(\theta)\\ is the
prior, and \\p(x_o) = \int p(x_o \mid \theta) p(\theta) \\ d\theta\\ is
the marginal likelihood or evidence. In simulation-based inference, the
likelihood \\p(x \mid \theta)\\ cannot be evaluated pointwise; it is
available only implicitly through a simulator \\\mathcal{S}\\ that,
given parameters \\\theta\\, produces synthetic data \\x \sim p(x \mid
\theta)\\ by executing a stochastic program. The three SBI families
differ in what they learn from simulated pairs \\\\(\theta_i,
x_i)\\\_{i=1}^N\\ drawn from the joint distribution \\p(\theta, x) = p(x
\mid \theta) p(\theta)\\.

Throughout what follows, \\\phi\\ denotes the parameters of the neural
network being trained, \\q\_\phi\\ denotes the density parametrized by
the network, \\d\\ denotes the dimensionality of the parameter space
\\\theta \in \mathbb{R}^d\\, and \\m\\ denotes the dimensionality of the
data or summary statistics \\x \in \mathbb{R}^m\\. The notation
\\\sigma(\cdot)\\ is used for the sigmoid function \\\sigma(z) = 1/(1 +
e^{-z})\\, and \\D\_{\text{KL}}(\cdot \\ \cdot)\\ for the
Kullback-Leibler divergence.

### Neural posterior estimation (NPE)

NPE trains a conditional density estimator \\q\_\phi(\theta \mid x)\\ to
approximate the true posterior \\p(\theta \mid x)\\ directly from the
simulated pairs [\[2\]](#ref2). The training objective is the expected
negative log-probability of the training parameters under the estimated
posterior:

\\\mathcal{L}\_{\text{NPE}}(\phi) = -\mathbb{E}\_{p(\theta,
x)}\\\left\[\log q\_\phi(\theta \mid x)\right\] = -\frac{1}{N}
\sum\_{i=1}^N \log q\_\phi(\theta_i \mid x_i).\\

This loss has a clean information-theoretic interpretation. Minimizing
\\\mathcal{L}\_{\text{NPE}}\\ over \\\phi\\ is equivalent to minimizing
the KL divergence \\D\_{\text{KL}}\\\left(p(\theta \mid x) \\\\\\
q\_\phi(\theta \mid x)\right)\\ averaged over the marginal distribution
of \\x\\, up to a constant that does not depend on \\\phi\\
[\[24\]](#ref24). To see this, expand the KL divergence:

\\D\_{\text{KL}}\\\left(p(\theta \mid x) \\\\\\ q\_\phi(\theta \mid
x)\right) = \int p(\theta \mid x) \log \frac{p(\theta \mid
x)}{q\_\phi(\theta \mid x)} \\ d\theta = -H\[p(\theta \mid x)\] - \int
p(\theta \mid x) \log q\_\phi(\theta \mid x) \\ d\theta,\\

where \\H\[p(\theta \mid x)\]\\ is the entropy of the true posterior,
which does not depend on \\\phi\\. Taking the expectation over \\p(x)\\
yields:

\\\mathbb{E}\_{p(x)}\\\left\[D\_{\text{KL}}\\\left(p(\theta \mid x)
\\\\\\ q\_\phi(\theta \mid x)\right)\right\] = \text{const} -
\mathbb{E}\_{p(\theta, x)}\\\left\[\log q\_\phi(\theta \mid x)\right\] =
\text{const} + \mathcal{L}\_{\text{NPE}}(\phi).\\

So minimizing \\\mathcal{L}\_{\text{NPE}}\\ drives \\q\_\phi\\ toward
the true posterior for each \\x\\ simultaneously, without ever needing
to evaluate the posterior density itself.

The default density estimator in lucifer is a Mixture Density Network
(MDN) [\[7\]](#ref7), which models the posterior as a mixture of
Gaussians conditioned on the observation:

\\q\_\phi(\theta \mid x) = \sum\_{k=1}^K \pi_k(x) \\
\mathcal{N}\\\left(\theta;\\ \mu_k(x),\\ \Sigma_k(x)\right),\\

where \\K\\ is the number of mixture components. The mixing weights are
produced by a softmax transformation of raw network outputs:

\\\pi_k(x) = \frac{\exp(\alpha_k(x))}{\sum\_{j=1}^K
\exp(\alpha_j(x))},\\

ensuring \\\pi_k \> 0\\ and \\\sum_k \pi_k = 1\\. The component means
\\\mu_k(x) \in \mathbb{R}^d\\ are direct network outputs. The covariance
matrices are parametrized through their Cholesky factors \\L_k(x)\\,
where \\\Sigma_k(x) = L_k(x) L_k(x)^\top\\, which ensures positive
definiteness by construction. The diagonal elements of \\L_k\\ are
exponentiated to enforce strict positivity, while the off-diagonal
elements are unconstrained. The raw network outputs are thus the
pre-softmax logits \\\alpha_k\\, the mean vectors \\\mu_k\\, and the
lower-triangular Cholesky entries \\L_k\\ for each component, totalling
\\K\\\left(1 + d + \frac{d(d+1)}{2}\right)\\ output units.

This full-covariance parametrization is critical for capturing
correlations between parameters, which diagonal-covariance models would
miss entirely. For a \\d\\-dimensional parameter space with \\K\\
components, the number of Cholesky parameters grows as \\O(K d^2)\\,
which becomes expensive for large \\d\\; this is the primary motivation
for switching to normalizing flows in high dimensions.

The log-density of the mixture is computed via the log-sum-exp trick to
avoid numerical underflow:

\\\log q\_\phi(\theta \mid x) = \log \sum\_{k=1}^K \pi_k(x) \\
\mathcal{N}(\theta; \mu_k, \Sigma_k) = \text{LSE}\_{k=1}^K\\\left\[\log
\pi_k + \log \mathcal{N}(\theta; \mu_k, \Sigma_k)\right\],\\

where \\\text{LSE}\\ denotes the log-sum-exp operation, and the
multivariate normal log-density for each component is:

\\\log \mathcal{N}(\theta; \mu_k, \Sigma_k) = -\frac{d}{2}\log(2\pi) -
\log\|L_k\| - \frac{1}{2}\\L_k^{-1}(\theta - \mu_k)\\^2.\\

The key property of NPE is amortization: once the network is trained,
posterior samples for any new observation \\x_o\\ are obtained by (i)
passing \\x_o\\ through the network to obtain the mixture parameters
\\\\\pi_k, \mu_k, L_k\\\\, (ii) sampling a component index \\k \sim
\text{Categorical}(\pi_1, \ldots, \pi_K)\\, and (iii) drawing \\\theta =
\mu_k + L_k \epsilon\\ where \\\epsilon \sim \mathcal{N}(0, I_d)\\. No
MCMC is required; the entire posterior sampling process is a single
forward pass plus a random draw, which takes milliseconds regardless of
the number of samples requested. This makes NPE ideal for scenarios
where inference must be repeated for many different datasets under the
same model, such as in population-level studies where each individual
produces a separate dataset but the model structure is shared.

### Neural likelihood estimation (NLE)

NLE trains the conditional density in the opposite direction:
\\q\_\phi(x \mid \theta) \approx p(x \mid \theta)\\, learning a
surrogate for the intractable likelihood [\[3\]](#ref3)
[\[11\]](#ref11). The training objective mirrors NPE but with roles
reversed:

\\\mathcal{L}\_{\text{NLE}}(\phi) = -\mathbb{E}\_{p(\theta,
x)}\\\left\[\log q\_\phi(x \mid \theta)\right\] = -\frac{1}{N}
\sum\_{i=1}^N \log q\_\phi(x_i \mid \theta_i).\\

The MDN architecture is the same as in NPE, but the conditioning
variable is now \\\theta\\ and the target variable is \\x\\. The network
takes \\\theta \in \mathbb{R}^d\\ as input and outputs the mixture
parameters for a distribution over \\x \in \mathbb{R}^m\\.

At inference time, the trained neural likelihood is combined with the
prior through Bayes’ theorem. The unnormalized log-posterior, which
serves as the MCMC target density, is:

\\\log p(\theta \mid x_o) \propto \log q\_\phi(x_o \mid \theta) + \log
p(\theta).\\

Here \\\log q\_\phi(x_o \mid \theta)\\ is evaluated by passing the
candidate \\\theta\\ through the trained network and computing the
log-density of the fixed observation \\x_o\\ under the output mixture
distribution. The log-prior \\\log p(\theta)\\ is either supplied by the
user via the `log_prior_fn` argument or approximated from the prior
draws.

This is where lucifer’s architecture creates a unique advantage. The
trained neural likelihood is wrapped as a standard `Model(parm, Data)`
closure that returns the list with `LP` (the log-posterior computed as
above), `Dev` (deviance, \\-2 \cdot \text{LP}\\), `Monitor` (monitored
quantities), `yhat` (predicted values), and `parm` (the parameter
vector, potentially after constraint handling). This closure plugs
directly into
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
with any of the 62 available MCMC algorithms, from simple random-walk
Metropolis to NUTS, Hamiltonian Monte Carlo, slice sampling,
differential evolution, and ensemble methods, without any modification
to the sampler code. No other R or Python package offers this level of
MCMC integration for neural likelihood estimation.

The model closure construction in lucifer proceeds through the following
steps. Given the trained network \\q\_\phi\\, the observed data \\x_o\\,
the prior function, and the standardization parameters (means
\\\bar{\theta}\\, \\\bar{x}\\ and standard deviations \\s\_\theta\\,
\\s_x\\ from the training data), the closure at each MCMC iteration:

1.  Receives the candidate parameter vector \\\theta\\ from the MCMC
    sampler.
2.  Standardizes \\\theta\\ to \\\theta^{\text{std}} = (\theta -
    \bar{\theta}) / s\_\theta\\.
3.  Evaluates \\\log q\_\phi(x_o^{\text{std}} \mid
    \theta^{\text{std}})\\ by computing the MDN mixture log-density of
    the standardized observation.
4.  Adds the log-prior \\\log p(\theta)\\ from either the user-supplied
    `log_prior_fn` or the approximate Gaussian prior.
5.  Returns the standard `LP`, `Dev`, `Monitor`, `yhat`, `parm` list
    expected by
    [`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

NLE is particularly effective for i.i.d. data, because one can train the
network on individual observations and aggregate the log-likelihood at
inference time via \\\log p(x_1, \ldots, x_n \mid \theta) =
\sum\_{i=1}^n \log q\_\phi(x_i \mid \theta)\\. This factorization avoids
the curse of dimensionality that would arise from training on the full
\\n\\-dimensional data vector, and it allows the network to be trained
on a lower-dimensional problem (\\m = 1\\ for scalar observations)
regardless of the sample size \\n\\.

### Neural ratio estimation (NRE)

NRE takes a classification approach to density ratio estimation
[\[5\]](#ref5). The key insight, rooted in the Neyman-Pearson lemma, is
that the likelihood-to-evidence ratio \\r(\theta, x) = p(\theta, x) /
\[p(\theta) p(x)\] = p(x \mid \theta) / p(x)\\ can be recovered from a
binary classifier trained to distinguish joint samples from marginal
samples.

Specifically, one constructs a balanced training set with two classes.
The “joint” class contains pairs \\(\theta_i, x_i) \sim p(\theta, x)\\
where \\x_i\\ was generated by the simulator from \\\theta_i\\, labelled
\\y = 1\\. The “marginal” class contains pairs \\(\theta_j, x_k) \sim
p(\theta) p(x)\\ where the parameter and data are drawn independently,
labelled \\y = 0\\. In practice, marginal pairs are constructed by
simply permuting the \\\theta\\ indices relative to the \\x\\ indices
within a batch, which preserves the correct marginal distributions while
breaking the dependence structure.

A binary classifier \\d\_\phi(\theta, x) = \sigma(f\_\phi(\theta, x))\\
is trained with the standard binary cross-entropy loss:

\\\mathcal{L}\_{\text{NRE}}(\phi) = -\frac{1}{2N} \sum\_{i=1}^N
\Big\[\log d\_\phi(\theta_i, x_i) + \log\big(1 - d\_\phi(\theta_i',
x_i')\big)\Big\],\\

where \\(\theta_i, x_i)\\ are joint pairs and \\(\theta_i', x_i')\\ are
marginal pairs. At the Bayes-optimal classifier, the sigmoid output
equals the class-posterior probability:

\\d\_\phi^\*(\theta, x) = \frac{p(\theta, x)}{p(\theta, x) + p(\theta)
p(x)} = \sigma\\\left(\log \frac{p(\theta, x)}{p(\theta) p(x)}\right),\\

and therefore the logit output converges to the log density ratio:

\\f\_\phi^\*(\theta, x) = \log \frac{p(\theta, x)}{p(\theta)p(x)} = \log
\frac{p(x \mid \theta)}{p(x)} = \log r(\theta, x).\\

Since \\p(x_o)\\ is a constant for fixed observed data, the MCMC target
for NRE becomes:

\\\log p(\theta \mid x_o) \propto f\_\phi(\theta, x_o) + \log
p(\theta),\\

and the same model closure mechanism as NLE is used to bridge into
lucifer’s MCMC algorithms. The \\\log r\\ term replaces the
log-likelihood; the prior term is added explicitly. The connection to
the Neyman-Pearson lemma is direct: the optimal test statistic for
distinguishing joint from marginal samples is precisely the likelihood
ratio, and the trained classifier approximates this optimal test
statistic through empirical risk minimization.

NRE has the simplest training objective among the three families (binary
cross-entropy with a standard classifier architecture) and requires no
density estimation machinery. It is also the most flexible in terms of
classifier architecture, since any binary classifier that accepts
concatenated \\(\theta, x)\\ inputs can serve as the density ratio
estimator. However, the resulting log-ratio surface may be less smooth
than NLE’s log-likelihood surface, sometimes requiring more careful MCMC
tuning or longer chains to achieve comparable effective sample sizes.

### Sequential methods

When the prior is broad relative to the posterior, most simulations from
\\p(\theta)\\ fall in regions of negligible posterior mass and
contribute little information to training. Consider a model with a
10-dimensional parameter space where the prior volume is \\10^6\\ times
the posterior volume: only about 1 in \\10^6\\ prior draws would fall in
the posterior’s high-density region, meaning that with \\N = 10^4\\
training simulations, essentially none of them inform the posterior.
Sequential SBI addresses this inefficiency by iterating the
simulation-training cycle across multiple rounds, progressively
concentrating the proposal distribution near the observed data
[\[2\]](#ref2) [\[9\]](#ref9).

In sequential NPE (SNPE), the procedure operates as follows. In round
\\r = 1\\, simulations are drawn from the prior \\p(\theta)\\ and a
preliminary posterior estimator \\q\_\phi^{(1)}(\theta \mid x)\\ is
trained. In subsequent rounds \\r \geq 2\\, the proposal distribution
\\\tilde{p}\_r(\theta)\\ is derived from the current posterior estimate,
and new simulations are drawn from this focused proposal. The truncated
sequential NPE (TSNPE) variant [\[10\]](#ref10), which is the default in
lucifer, defines the proposal as a truncation of the prior to a region
\\\mathcal{R}\_{r-1}\\ determined by the previous round’s posterior:

\\\tilde{p}\_r(\theta) \propto p(\theta) \cdot \mathbb{1}\[\theta \in
\mathcal{R}\_{r-1}\],\\

where \\\mathcal{R}\_{r-1}\\ is the Cartesian product of intervals
defined using robust bounds based on the interquartile range (IQR) of
the posterior samples from round \\r - 1\\:

\\\mathcal{R}\_{r-1}^{(j)} = \left\[\text{median}(\theta^{(j)}) - 5
\cdot \text{IQR}(\theta^{(j)}),\\ \text{median}(\theta^{(j)}) + 5 \cdot
\text{IQR}(\theta^{(j)})\right\],\\

where \\j\\ indexes the parameter dimensions. The factor of 5 ensures
that the truncation region is generous enough to avoid cutting off
posterior mass (for a Gaussian distribution, \\\pm 5 \cdot \text{IQR}\\
covers roughly \\\pm 6.7\sigma\\, encompassing \\\> 99.99\\\\ of the
mass) while still excluding the vast majority of the prior volume that
contributes nothing to inference.

The TSNPE approach has two important advantages over earlier sequential
methods that used importance weighting. First, it avoids
importance-weighting corrections entirely, which can suffer from high
variance when the proposal and posterior differ substantially. Second,
the truncation is purely geometric and does not modify the loss
function, so the same NPE loss can be used in every round. The only
modification is the restriction of the prior draws to
\\\mathcal{R}\_{r-1}\\, which is implemented as simple rejection
sampling.

Sequential NLE (SNLE) and sequential NRE (SNRE) follow the same
multi-round logic [\[11\]](#ref11) [\[12\]](#ref12). In each round, the
current posterior estimate (obtained via MCMC using the current neural
likelihood or ratio) serves as the proposal for generating new training
data in the next round. For SNLE and SNRE, the surrogate likelihood or
ratio is retrained from scratch on the accumulated data from all rounds,
since the MCMC-based posterior from the current round provides a natural
proposal for the next. The total simulation budget \\N\\ is divided
across \\R\\ rounds, with \\N/R\\ simulations per round by default,
though lucifer allows the per-round budget to be specified explicitly
via the `n_sims_per_round` argument.

### Comparison with ABC

Both ABC and SBI handle intractable likelihoods, but their mechanisms
are fundamentally different, and understanding these differences is
essential for choosing between them.

ABC generates posterior samples by an accept/reject mechanism: draw
\\\theta \sim p(\theta)\\, simulate \\x\_{\text{sim}} \sim p(x \mid
\theta)\\, and accept \\\theta\\ if \\\rho(S(x\_{\text{sim}}), S(x_o))
\< \epsilon\\, where \\\rho\\ is a distance metric, \\S(\cdot)\\
extracts summary statistics, and \\\epsilon\\ is the tolerance. The
accepted \\\theta\\ values approximate samples from \\p(\theta \mid
\rho(S(x), S(x_o)) \< \epsilon)\\, which converges to the true posterior
\\p(\theta \mid x_o)\\ as \\\epsilon \to 0\\ (assuming sufficient
statistics are used). However, the acceptance rate decreases
exponentially with the dimension of the summary statistics, scaling
roughly as \\O(\epsilon^m)\\ where \\m\\ is the summary statistic
dimension [\[17\]](#ref17). For a 10-dimensional summary statistic
vector with \\\epsilon = 0.1\\ (relative to the data range), the
expected acceptance rate is approximately \\10^{-10}\\, requiring
billions of simulations.

SBI uses every simulation as a training example. The neural density
estimator learns to map from observations (or summaries) to posterior
distributions across the entire parameter space simultaneously,
amortizing the information from all simulations into a single set of
network weights. This makes SBI substantially more sample-efficient than
ABC, typically requiring \\10^3\\ to \\10^4\\ simulations where ABC
might need \\10^5\\ to \\10^7\\ for comparable accuracy. SBI also avoids
the need to choose a distance metric \\\rho\\ and tolerance
\\\epsilon\\, though it introduces its own hyperparameters (network
architecture, number of mixture components, learning rate, training
schedule).

The tradeoff is that SBI requires more upfront computational investment
in training the neural network and more expertise in diagnosing the
quality of the learned density estimator, whereas ABC produces samples
immediately (if slowly) and has a simpler theoretical foundation. For
single-use inference on a fixed dataset, the total cost may be
comparable; for repeated inference across many datasets with the same
model structure, NPE’s amortization property gives it a decisive
advantage. lucifer provides both ABC (via the
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
function) and SBI (via
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)),
making it straightforward to compare results from both approaches on the
same problem.

### Convergence properties

Under mild regularity conditions (compact parameter space, continuous
and bounded simulator, universal approximation by the density estimator
class), the NPE posterior \\q\_\phi(\theta \mid x)\\ converges to the
true posterior \\p(\theta \mid x)\\ as the number of simulations \\N \to
\infty\\ and the network capacity grows [\[11\]](#ref11)
[\[24\]](#ref24). Specifically, for mixture density networks with
sufficiently many components \\K\\, the class of MDN densities is dense
in the space of continuous conditional densities on compact domains, and
the empirical risk \\\hat{\mathcal{L}}\_{\text{NPE}}\\ converges to the
population risk \\\mathcal{L}\_{\text{NPE}}\\ at the standard
\\O(1/\sqrt{N})\\ rate, ensuring that the KL divergence between the
learned and true posteriors vanishes.

In practice, convergence depends critically on three factors: the
simulation budget \\N\\ (which determines the statistical error), the
network architecture (which determines the approximation error), and the
optimization procedure (which determines the optimization error). The
total error decomposes as:

\\D\_{\text{KL}}(p \\ q\_\phi) \leq \underbrace{D\_{\text{KL}}(p \\
q\_{\phi^\*})}\_{\text{approximation}} + \underbrace{\|D\_{\text{KL}}(p
\\ q\_\phi) - D\_{\text{KL}}(p \\
q\_{\phi^\*})\|}\_{\text{optimization}} +
\underbrace{\|\hat{\mathcal{L}} - \mathcal{L}\|}\_{\text{estimation}},\\

where \\\phi^\*\\ is the best achievable parameter within the model
class. Increasing the number of mixture components \\K\\ and hidden
units reduces the approximation error; more simulations reduce the
estimation error; better optimization (lower learning rate, more epochs,
gradient clipping) reduces the optimization error.

For NLE and NRE, the additional MCMC step introduces a fourth source of
error: the MCMC mixing error, which arises from finite chain length and
imperfect convergence. However, this error is well-understood and can be
controlled through standard convergence diagnostics (effective sample
size, split-\\\hat{R}\\, trace plots), which is precisely why lucifer’s
MCMC bridge is valuable: the full suite of convergence assessment tools
from classical MCMC applies directly to the neural-likelihood-based
posterior samples.

### Unified view of the three methods

The three SBI families can be understood as different factorizations of
the same joint distribution \\p(\theta, x) = p(x \mid \theta) p(\theta)
= p(\theta \mid x) p(x)\\. NPE directly learns the posterior factor
\\p(\theta \mid x)\\, giving immediate access to posterior samples but
no access to the likelihood. NLE learns the likelihood factor \\p(x \mid
\theta)\\, which must be combined with the prior via MCMC to obtain
posterior samples, but provides a reusable likelihood surrogate that can
be combined with different priors without retraining. NRE learns the
ratio \\p(\theta, x) / \[p(\theta) p(x)\]\\, which equals \\p(x \mid
\theta) / p(x)\\ and differs from the NLE output by a normalizing
constant that cancels in MCMC acceptance ratios.

This factorization perspective reveals a key practical distinction. NPE
is a “one-shot” method: it produces posterior samples directly, but
those samples are tied to the prior used during training. If the user
wants to change the prior (e.g., to perform a sensitivity analysis), NPE
must be retrained from scratch. NLE, by contrast, provides a
prior-independent likelihood surrogate: the same trained neural
likelihood can be combined with any prior via MCMC, simply by changing
the `log_prior_fn` in the model closure. This makes NLE more flexible
for prior sensitivity analysis and for situations where the same
simulator might be used with different prior specifications across
analyses or collaborators. NRE falls between the two: the ratio \\p(x
\mid \theta) / p(x)\\ depends on the marginal \\p(x)\\ which implicitly
depends on the prior, but in practice this dependence is mild and NRE
can be approximately reused across similar priors.

### MDN training algorithm

The MDN training procedure in lucifer follows the standard mini-batch
stochastic gradient descent paradigm with the Adam optimizer
[\[19\]](#ref19). At each epoch, the training data is shuffled and
divided into mini-batches of size `batch_size` (default 256). For each
mini-batch \\\mathcal{B}\\, the forward pass computes the mixture
parameters \\\\\pi_k, \mu_k, L_k\\\\ for each training observation,
evaluates the negative log-mixture-density for each training pair, and
averages to obtain the batch loss:

\\\hat{\mathcal{L}}\_{\mathcal{B}} = -\frac{1}{\|\mathcal{B}\|}
\sum\_{(\theta_i, x_i) \in \mathcal{B}} \log \sum\_{k=1}^K \pi_k(x_i) \\
\mathcal{N}(\theta_i; \mu_k(x_i), L_k(x_i) L_k(x_i)^\top).\\

The backward pass computes gradients \\\nabla\_\phi
\hat{\mathcal{L}}\_{\mathcal{B}}\\ via backpropagation through the SELU
activations, the softmax (for mixing weights), the Cholesky
parametrization (for covariances), and the log-mixture-density formula.
The gradients are clipped to have norm at most `max_grad_norm`, then
passed to Adam, which maintains per-parameter first and second moment
estimates and computes bias-corrected updates.

At the end of each epoch, the validation loss is computed on the
held-out validation set (default 10% of data). If the validation loss is
the lowest observed so far, the current network weights are saved as the
“best” weights. If the validation loss has not improved for `patience`
consecutive epochs (default 20), training terminates and the best
weights are restored. This early stopping mechanism is critical for
preventing overfitting, which is the dominant source of SBI failure for
small simulation budgets.

The entire training procedure, including the forward pass, backward
pass, gradient clipping, Adam updates, and validation evaluation, is
implemented in C++ for the MDN backend, avoiding the overhead of R-level
loops over epochs and mini-batches. The only R-level call during
training is the progress reporting (when `verbose = TRUE`), which is
invoked at most once per epoch. For a typical problem with 10000
simulations, 128-128 hidden layers, 10 mixture components, and batch
size 256, the C++ training loop completes in 1-5 seconds for 200 epochs,
which is negligible compared to the simulation phase for any non-trivial
simulator.

For the NRE classifier, training follows the same epoch-based structure
but with the binary cross-entropy loss instead of the negative
log-mixture-density loss. The training data are organized as pairs of
“joint” and “marginal” batches of equal size, and the gradients from
both batches are combined before the Adam update step. The NRE
classifier is lighter-weight than the MDN (a single logit output rather
than \\K(1 + d + d(d+1)/2)\\ outputs), so training is typically faster
for NRE than for NPE or NLE at the same network depth and width.

## Implementation architecture

lucifer’s SBI module is built on two C++ backends that require no
external dependencies beyond the standard R toolchain. The primary
backend is a Mixture Density Network implemented entirely in C++ with
Armadillo linear algebra, providing the density estimation machinery for
both NPE and NLE. The secondary backend is a binary classifier, also in
C++, that powers NRE. Both backends are compiled as part of the package
installation and require no runtime compilation or external libraries.

The MDN backend uses SELU (Scaled Exponential Linear Unit) activations
[\[20\]](#ref20), defined as \\\text{SELU}(z) = \lambda \cdot
\begin{cases} z & \text{if } z \> 0 \\ \alpha(e^z - 1) & \text{if } z
\leq 0 \end{cases}\\ with \\\lambda \approx 1.0507\\ and \\\alpha
\approx 1.6733\\. SELU activations provide self-normalizing properties:
under certain conditions on the weight initialization, the activations
in each layer converge toward zero mean and unit variance, which
stabilizes training without the need for batch normalization layers.
This is particularly important for the MDN because the output layer must
produce Cholesky factors and means with specific scale properties, and
unstable internal activations would propagate errors to these critical
outputs.

The Adam optimizer [\[19\]](#ref19) handles parameter updates with
configurable learning rate and gradient clipping (controlled via
`max_grad_norm`). Gradient clipping prevents the exploding-gradient
pathology that can occur when the Cholesky factors of the covariance
matrices are poorly conditioned, which happens when the network predicts
near-singular covariances during early training. The clipping rescales
the gradient vector to have norm at most `max_grad_norm` whenever the
raw gradient norm exceeds this threshold:

\\g \leftarrow g \cdot \frac{\text{max\\grad\\norm}}{\max(\\g\\,
\text{max\\grad\\norm})}.\\

Early stopping with patience monitoring on a held-out validation set
prevents overfitting. At each epoch, the validation loss is computed on
a fraction (default 10%) of the data withheld from training. If the
validation loss fails to improve for `patience` consecutive epochs
(default 20), training terminates and the network weights from the best
epoch are restored. This is particularly important because the MDN has
\\O(K d^2)\\ output parameters per mixture component and can easily
memorize small training sets, producing posterior estimates that are
overconfident on the training data but poorly calibrated on new
observations.

The NRE classifier uses the same SELU-activated hidden layers but
terminates in a single logit output rather than the complex mixture
parametrization of the MDN. Training uses binary cross-entropy with the
numerically stable log-sigmoid formulation \\\log \sigma(z) =
-\text{softplus}(-z) = -\log(1 + e^{-z})\\ to avoid overflow when the
classifier becomes highly confident. Both backends support OpenMP
parallelization for the matrix operations in the forward and backward
passes, though the dominant computational cost for most problems is the
simulation step rather than the network training.

A standardization pipeline normalizes both parameter vectors and data
vectors to zero mean and unit variance before training. This is
essential for stable MDN training because the mixture component means
and Cholesky factors are initialized near zero, and unstandardized data
with disparate scales would require the network to learn large-magnitude
outputs before any useful density estimation could begin. The
standardization statistics (means \\\bar{\theta}\\, \\\bar{x}\\ and
standard deviations \\s\_\theta\\, \\s_x\\) are computed from the
training set, stored in the fitted `sbi` object, and applied
automatically during posterior sampling, MCMC model closure evaluation,
and diagnostic computation.

For NPE, posterior sampling applies rejection sampling within the
estimated prior support to suppress artifacts from the MDN’s Gaussian
tails, which can occasionally place mass outside the prior’s domain. The
support bounds are estimated as the range of the training parameters
plus a 10% margin, and samples falling outside these bounds are
discarded and redrawn for up to 10 attempts. For NLE and NRE, the model
closure bridge constructs a `Model(parm, Data)` function that lucifer’s
MCMC dispatcher expects, with the neural log-likelihood or log-ratio as
the LP component and the user’s log-prior added explicitly. If no
`log_prior_fn` is supplied, an approximate log-prior is constructed from
a multivariate Gaussian fitted to the prior draws, with hard bounds
enforcing the prior support.

The simulation pipeline deserves mention because it can be the
computational bottleneck for expensive simulators. In sequential mode
(`n_cores = 1`), the simulation loop is implemented in C++ via
`sbi_simulate_batch_cpp()`, which calls the R prior and simulator
functions from C++ in a tight loop, handles non-finite outputs
gracefully (incrementing a failure counter and continuing), and reports
progress at configurable intervals. For parallel simulation
(`n_cores > 1`), the pipeline switches to R-level
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
which forks the R process across cores. Each forked process
independently draws parameters and runs the simulator, with failed
simulations producing `NULL` results that are filtered before
aggregation. The failure rate is reported at the end, and if more than
99% of simulations fail, an error is thrown with a diagnostic message
suggesting that the user check the simulator and prior functions.

For problems exceeding roughly 20 parameters or requiring more flexible
density estimators, the optional torch backend provides Neural Spline
Flows (NSF) [\[8\]](#ref8) via the torch R package, which provides
native R bindings to libtorch with no Python dependency. The backend is
selected automatically when `backend = "auto"` and the parameter
dimension exceeds 20, or can be forced with `backend = "torch"`. The
[`sbi_network()`](https://robustecologies.github.io/lucifer/reference/sbi_network.md)
constructor allows explicit specification of the architecture type
(`"mdn"` or `"nsf"`) and layer configuration.

The architecture decision between `"mdn"` and `"nsf"` can be summarized
as follows. The MDN backend is preferred when the parameter dimension is
low to moderate (\\d \leq 20\\), when installation simplicity matters
(no torch dependency), and when the posterior is expected to be
well-approximated by a mixture of Gaussians. The NSF backend is
preferred when the parameter dimension is high (\\d \> 20\\), when the
posterior has complex geometry that Gaussian mixtures cannot represent,
or when the user has GPU resources available for faster training. Both
backends produce objects with the same interface, so switching between
them requires changing only the `backend` argument or the
[`sbi_network()`](https://robustecologies.github.io/lucifer/reference/sbi_network.md)
specification, with no changes to downstream analysis code.

A key implementation detail is that the
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
function performs a pilot simulation before training to determine the
dimensionality of both the parameter space and the data/summary space.
It draws a single parameter from the prior, runs the simulator (and
summary function, if provided), and verifies that the output dimensions
match `x_obs`. This pilot check catches dimension mismatches early,
before any expensive simulation or training begins. The function also
validates that the prior and simulator produce finite outputs and that
`log_prior_fn` (if provided) is callable on a prior draw.

### Computational cost breakdown

The total wall-clock time for an SBI run decomposes into three phases:
simulation, training, and posterior sampling (or MCMC). For the C++ MDN
backend with typical settings (10000 simulations, 128-128 hidden layers,
10 components, 200 max epochs), the approximate timings are as follows.
The simulation phase takes \\N \times t\_{\text{sim}}\\, where
\\t\_{\text{sim}}\\ is the per-call time of the simulator; for a
simulator that takes 1 millisecond per call, 10000 simulations take
about 10 seconds; for a simulator taking 1 second per call, the
simulation phase takes nearly 3 hours. The training phase takes
approximately 1-10 seconds for the C++ MDN backend with 10000
simulations (the exact time depends on the output layer size, which
scales as \\O(Kd^2)\\). The posterior sampling phase takes negligible
time for NPE (a single matrix multiplication) or 5-60 seconds for
NLE/NRE MCMC with 5000 iterations.

This breakdown highlights a critical practical point: for any
non-trivial simulator, the simulation phase completely dominates the
computational cost. Network training and posterior sampling are
negligible by comparison. This means that the primary lever for reducing
total computation time is reducing the number of simulations required,
which motivates the use of sequential methods (which achieve better
posterior quality per simulation through focused proposals) and careful
summary statistic selection (which keeps the density estimation problem
low-dimensional, allowing the network to train with fewer examples).

## Worked example 1: normal-normal conjugate model

``` r

library(lucifer)
```

The normal-normal conjugate model provides an analytically tractable
posterior against which SBI accuracy can be verified exactly. The prior
is \\\theta \sim \mathcal{N}(0, 5^2)\\ and the likelihood is \\x_i \mid
\theta \sim \mathcal{N}(\theta, 1)\\ for \\i = 1, \ldots, 10\\, giving a
conjugate posterior \\\theta \mid x \sim \mathcal{N}(\mu\_{\text{post}},
\sigma^2\_{\text{post}})\\ where

\\\sigma^2\_{\text{post}} = \left(\frac{1}{\sigma_0^2} +
\frac{n}{\sigma^2}\right)^{-1} = \left(\frac{1}{25} + 10\right)^{-1}
\approx 0.0996\\

and

\\\mu\_{\text{post}} = \sigma^2\_{\text{post}}
\left(\frac{\mu_0}{\sigma_0^2} + \frac{\sum\_{i=1}^n
x_i}{\sigma^2}\right) = \sigma^2\_{\text{post}} \cdot \sum\_{i=1}^n
x_i.\\

This model is ideal for benchmarking because the simulator is trivial to
implement (a single call to `rnorm`), the posterior has a known closed
form, and the posterior is unimodal and Gaussian, so even a
single-component MDN should be able to approximate it accurately. Any
systematic deviation from the analytical posterior indicates a problem
with the training procedure, the simulation budget, or the
standardization pipeline.

``` r

set.seed(42)

# Define simulator and prior
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)

# Observed data (10 observations)
x_obs_raw <- c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0, 2.4, 1.7, 2.2, 2.6)

# Summary statistics reduce the 10D data to informative low-dimensional
# features. For the normal-normal model, the sample mean and standard
# deviation are sufficient statistics.
summary_fn <- function(x) c(mean(x), sd(x))
x_obs <- summary_fn(x_obs_raw)

# Analytical posterior for comparison
true_sigma2 <- 1 / (1/25 + 10)
true_mu <- true_sigma2 * sum(x_obs_raw)
true_sd <- sqrt(true_sigma2)

cat(sprintf("Analytical posterior: N(%.4f, %.4f)\n", true_mu, true_sd))

# NPE: direct posterior estimation with summary statistics
fit_npe <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NPE",
    summary_fn = summary_fn,
    n_simulations = 10000,
    n_samples = 5000,
    n_components = 3,
    learning_rate = 1e-3,
    max_epochs = 200,
    seed = 42,
    verbose = TRUE
)
```

The
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
function handles the full pipeline: generating simulations from the
prior, applying the summary function to both simulated and observed
data, splitting into training and validation sets, standardizing both
\\\theta\\ and \\x\\ to zero mean and unit variance, training the MDN
with Adam and early stopping, and drawing posterior samples via the
rejection-sampled forward pass. The `seed` argument ensures
reproducibility by setting R’s random number generator state before
simulation begins. Using `summary_fn` is generally recommended when the
raw data dimension exceeds 3-5, since neural density estimators require
exponentially more training data in higher dimensions.

``` r

# Inspect the fit
print(fit_npe)
```

The print method provides a concise one-screen summary: the number of
parameters, simulations, and posterior samples; the network architecture
(hidden layer sizes and number of mixture components); the final
training and validation loss; and a table of posterior summaries
including mean, standard deviation, and 95% credible intervals.

``` r

summary(fit_npe)
```

The summary method extends the print output with details about the
network backend (MDN or torch), the full hidden layer configuration,
training progression (number of epochs trained, best epoch selected by
early stopping, final training and validation losses), and posterior
correlations for multivariate problems. For NLE/NRE fits, it also
includes MCMC diagnostics from the underlying
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md)
call.

``` r

# Compare NPE posterior with analytical truth
cat(sprintf("Analytical posterior mean: %.4f\n", true_mu))
cat(sprintf("NPE posterior mean:       %.4f\n", mean(fit_npe$Posterior[, 1])))
cat(sprintf("Analytical posterior SD:   %.4f\n", true_sd))
cat(sprintf("NPE posterior SD:          %.4f\n", sd(fit_npe$Posterior[, 1])))
cat(sprintf("Absolute error in mean:   %.4f\n",
            abs(mean(fit_npe$Posterior[, 1]) - true_mu)))
```

With 5000 simulations and a 10-component MDN, the NPE posterior should
closely match the analytical result. Deviations of a few hundredths in
the mean and standard deviation are typical and reflect the finite
simulation budget rather than any systematic bias. The errors decrease
approximately as \\O(1/\sqrt{N})\\ with the simulation budget.

``` r

# Posterior density with ground truth
plot(fit_npe, type = "posterior",
     ground_truth = c("theta[1]" = true_mu))
```

The posterior density plot shows the marginal posterior for each
parameter (here just \\\theta_1\\) as a kernel density estimate from the
posterior samples, with a dashed vertical line at the analytical
posterior mean. The density should be approximately Gaussian and
centered near the truth. Substantial asymmetry or multimodality in this
plot for the conjugate normal-normal model would indicate a training
problem.

``` r

# Training loss curve
plot(fit_npe, type = "training")
```

The training curve shows the negative log-likelihood loss on the
training set and validation set across epochs, with a dashed vertical
line at the best epoch selected by early stopping. A well-behaved
training run shows both curves decreasing together initially, with the
validation curve leveling off or increasing slightly as overfitting
begins, at which point training stops and the best-epoch weights are
restored. If the two curves never separate (no overfitting), the network
may benefit from increased capacity (`n_components`, `hidden_layers`) or
the `max_epochs` limit may have been reached prematurely.

For this conjugate model, the posterior is so simple (a single Gaussian)
that a single-component MDN would suffice. The 10-component MDN is
deliberately over-specified to demonstrate that the training procedure
handles excess capacity gracefully: the extra components receive
near-zero mixing weights \\\pi_k\\ and do not distort the posterior.
This robustness to over-specification is an important practical
property, because in real applications the user does not know how many
components are needed and must err on the side of providing too many
rather than too few.

The `Posterior` field of the returned object is a matrix with one row
per sample and one column per parameter. For downstream analysis, these
samples can be used directly to compute any posterior functional: means,
medians, quantiles, highest posterior density intervals, or derived
quantities. For example:

``` r

# 95% highest posterior density interval (equal-tailed)
ci <- quantile(fit_npe$Posterior[, 1], probs = c(0.025, 0.975))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", ci[1], ci[2]))
cat(sprintf("Analytical 95%% CI: [%.4f, %.4f]\n",
            qnorm(0.025, true_mu, true_sd),
            qnorm(0.975, true_mu, true_sd)))

# Posterior probability that theta > 2
prob_gt2 <- mean(fit_npe$Posterior[, 1] > 2)
cat(sprintf("P(theta > 2 | data) = %.4f\n", prob_gt2))
cat(sprintf("Analytical: %.4f\n", 1 - pnorm(2, true_mu, true_sd)))
```

## The NLE-MCMC bridge

The NLE method bridges neural density estimation into lucifer’s full
MCMC infrastructure, and this is the feature that most sharply
distinguishes this implementation from other SBI packages. After
training a neural likelihood \\q\_\phi(x \mid \theta)\\, the surrogate
likelihood is wrapped as a standard `Model(parm, Data)` closure and
passed to
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
giving the user access to every MCMC algorithm in the package. This
means that the full power of adaptive Metropolis, NUTS, Hamiltonian
Monte Carlo, slice sampling, differential evolution, ensemble methods,
and gradient-based samplers is available for posterior exploration under
a neural surrogate likelihood, with all the convergence diagnostics and
post-processing tools that come with standard MCMC output.

The bridge works transparently: the user specifies `method = "NLE"` and
`mcmc_algorithm = "NUTS"` (or any other algorithm), and
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
internally constructs the model closure, sets up the `Data` list with
appropriate `parm.names` and `mon.names`, generates initial values from
a prior draw, creates default `Specs` for the chosen algorithm, and
calls
[`lucifer()`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
The returned `sbi` object contains both the neural network training
results and the full `demonoid` object from the MCMC run, accessible via
`fit$mcmc_fit`.

To understand what the bridge does concretely, consider the model
closure that
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
constructs for NLE. The closure is a function `Model(parm, Data)` that,
at each MCMC iteration, receives a candidate parameter vector `parm`
from the sampler and returns a list with five elements. The `LP` element
contains the log-posterior density, computed as the sum of the neural
log-likelihood \\\log q\_\phi(x_o^{\text{std}} \mid
\theta^{\text{std}})\\ and the log-prior \\\log p(\theta)\\. The `Dev`
element is the deviance \\-2 \cdot \text{LP}\\. The `Monitor` element
contains the parameter values (for recording in the trace). The `yhat`
element contains the predicted observation (here, the expected value of
the neural likelihood). The `parm` element returns the parameter vector,
potentially after applying constraints. This is exactly the same
interface that any hand-coded lucifer model uses, which is why the
entire MCMC infrastructure works without modification.

The default `Specs` for MCMC algorithms are generated automatically
based on the algorithm choice and parameter dimension. For NUTS, the
default `Specs` includes an adaptation length proportional to the number
of iterations, an initial step size, and a target acceptance rate. For
Slice sampling, the default `Specs` includes the slice width. For
adaptive Metropolis methods, the default `Specs` includes initial
proposal scales. The user can override these defaults via the
`mcmc_specs` argument if the automatic settings are not appropriate.

``` r

library(lucifer)

set.seed(42)
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)
x_obs_raw <- c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0, 2.4, 1.7, 2.2, 2.6)
summary_fn <- function(x) c(mean(x), sd(x))
x_obs <- summary_fn(x_obs_raw)

# Define an exact log-prior for better MCMC performance
log_prior <- function(theta) {
    dnorm(theta[1], mean = 0, sd = 5, log = TRUE)
}

# Train the neural likelihood and sample via Slice sampler.
# Slice is more robust than NUTS for neural likelihoods because NUTS
# relies on gradients of the log-likelihood surface, which can be
# unreliable when the MDN has flat regions from wide mixture components.
fit_nle <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NLE",
    summary_fn = summary_fn,
    n_simulations = 10000,
    log_prior_fn = log_prior,
    mcmc_algorithm = "Slice",
    mcmc_iterations = 5000,
    n_components = 5,
    learning_rate = 1e-3,
    verbose = TRUE
)

print(fit_nle)
```

The `log_prior_fn` argument deserves particular attention. When it is
`NULL` (the default), NLE constructs an approximate log-prior by fitting
a multivariate Gaussian to the prior draws, with hard boundary
enforcement at the empirical support of the training data. This
approximation is sufficient for priors that are roughly Gaussian, but it
can introduce bias for heavily skewed, truncated, or multimodal priors.
When the user provides an exact log-prior function, the MCMC target
density uses it directly:

\\\text{LP}(\theta) = \underbrace{\log q\_\phi(x_o \mid
\theta)}\_{\text{neural likelihood}} + \underbrace{\log
p(\theta)}\_{\text{exact prior}},\\

and the only approximation in the entire inference pipeline is the
neural likelihood itself. Providing `log_prior_fn` is strongly
recommended whenever the prior has an analytic form.

``` r

# The mcmc_fit field contains the full demonoid object
cat("MCMC algorithm used:", fit_nle$mcmc_fit$Algorithm, "\n")
cat("Iterations:", fit_nle$mcmc_fit$Iterations, "\n")
```

The bridge works with any of lucifer’s 62 algorithms. To demonstrate the
versatility, the same neural likelihood can be sampled with different
MCMC methods by simply changing the `mcmc_algorithm` argument:

``` r

# Slice sampling on the same neural likelihood
fit_nle_slice <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NLE",
    summary_fn = summary_fn,
    n_simulations = 10000,
    log_prior_fn = log_prior,
    mcmc_algorithm = "Slice",
    mcmc_iterations = 5000,
    n_components = 5,
    learning_rate = 1e-3,
    verbose = TRUE
)

# Adaptive Metropolis-within-Gibbs
fit_nle_amwg <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NLE",
    summary_fn = summary_fn,
    n_simulations = 10000,
    log_prior_fn = log_prior,
    mcmc_algorithm = "AMWG",
    mcmc_iterations = 10000,
    n_components = 5,
    learning_rate = 1e-3,
    verbose = TRUE
)

# NRE with NUTS (often gives sharper posteriors)
fit_nre_nuts <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs_raw,
    method = "NRE",
    n_simulations = 10000,
    log_prior_fn = log_prior,
    mcmc_algorithm = "NUTS",
    mcmc_iterations = 5000,
    verbose = TRUE
)

# Compare posterior summaries across algorithms
cat("NLE+NUTS  mean:", round(mean(fit_nle$Posterior[, 1]), 4), "\n")
cat("NLE+Slice mean:", round(mean(fit_nle_slice$Posterior[, 1]), 4), "\n")
cat("NLE+AMWG  mean:", round(mean(fit_nle_amwg$Posterior[, 1]), 4), "\n")
cat("NRE+NUTS  mean:", round(mean(fit_nre_nuts$Posterior[, 1]), 4), "\n")
```

If all four algorithms have converged, the posterior means should agree
to within Monte Carlo error. Substantial disagreement between algorithms
is a red flag indicating either insufficient network training (the
neural likelihood is a poor surrogate) or insufficient MCMC iterations
(the chains have not converged). The ability to cross-check results
across multiple MCMC algorithms is one of the key practical benefits of
the bridge architecture: algorithmic discrepancies serve as an implicit
diagnostic for neural likelihood quality.

The same bridge architecture applies to NRE. The neural ratio
estimator’s logit output serves as the log-likelihood component in the
MCMC target:

``` r

# NRE with slice sampler. NRE works well on raw data (no summary_fn needed)
# because it learns a ratio rather than a density.
fit_nre <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs_raw,
    method = "NRE",
    n_simulations = 10000,
    log_prior_fn = log_prior,
    mcmc_algorithm = "Slice",
    mcmc_iterations = 5000,
    verbose = TRUE
)

print(fit_nre)
summary(fit_nre)
```

## Worked example 2: Ricker ecological model

The Ricker model is a classic discrete-time population dynamics model
widely used in fisheries science and ecology [\[15\]](#ref15). The
population size \\N_t\\ at time \\t\\ evolves according to

\\N\_{t+1} = r \\ N_t \exp\\\left(-\frac{N_t}{K} + \epsilon_t\right),
\quad \epsilon_t \sim \mathcal{N}(0, \sigma^2),\\

where \\r\\ is the intrinsic growth rate, \\K\\ is the carrying
capacity, and \\\sigma\\ controls process noise. This model generates
complex nonlinear dynamics including stable equilibria, cycles, and
chaos depending on the parameter values, and its likelihood is
analytically intractable due to the stochastic process noise and the
nonlinear state transitions. It is a canonical test problem for
likelihood-free inference methods.

We parameterize on the log scale (\\\log r\\, \\\log K\\) to ensure
positivity and use three summary statistics: the mean, standard
deviation, and lag-1 autocorrelation of the population time series.
These statistics capture the central tendency, variability, and temporal
structure of the dynamics, providing a low-dimensional but informative
representation of the 100-step time series. The mean is informative
about the equilibrium population size (which depends primarily on
\\K\\), the standard deviation captures the amplitude of fluctuations
(which depends on both \\r\\ and \\\sigma\\), and the lag-1
autocorrelation captures the temporal persistence (which depends
primarily on \\r\\ through the stability properties of the deterministic
skeleton). Together, these three statistics provide a reasonable summary
of the dynamics for the purpose of inferring \\r\\ and \\K\\.

The choice of summary statistics is crucial in SBI just as it is in ABC:
insufficient statistics lose information about the parameters, leading
to wider posteriors than necessary, while high-dimensional summaries
make the density estimation problem harder and require more training
simulations. For ecological time series models, a useful strategy is to
start with low-order moments (mean, variance) and temporal correlations
(lag-1, lag-2 autocorrelation), then add distributional statistics
(quantiles, skewness) if the posterior is too diffuse. The `summary_fn`
mechanism in
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
makes it straightforward to iterate on the summary statistic choice
without modifying the simulator itself.

``` r

library(lucifer)

# Ricker population model simulator
ricker_sim <- function(theta) {
    r <- exp(theta[1])     # growth rate on log-scale
    K <- exp(theta[2])     # carrying capacity on log-scale
    sigma <- 0.1           # process noise (known)

    N <- numeric(100)
    N[1] <- K / 2
    for (t in 2:100) {
        N[t] <- N[t - 1] * exp(r * (1 - N[t - 1] / K) +
                                 rnorm(1, 0, sigma))
        N[t] <- max(N[t], 0.01)  # floor to prevent extinction
    }

    # Summary statistics: mean, SD, lag-1 autocorrelation
    c(mean(N), sd(N), cor(N[-1], N[-100]))
}

# Prior on log-parameters
ricker_prior <- function() {
    c(rnorm(1, 0.5, 0.5),   # log(r) ~ N(0.5, 0.5)
      rnorm(1, 5, 1))        # log(K) ~ N(5, 1)
}

# Generate "observed" summary statistics from known true parameters
set.seed(42)
true_theta <- c(log(2), log(150))
x_obs_eco <- ricker_sim(true_theta)

cat(sprintf("True parameters: log(r) = %.3f, log(K) = %.3f\n",
            true_theta[1], true_theta[2]))
cat(sprintf("Observed summaries: mean = %.2f, sd = %.2f, acf1 = %.3f\n",
            x_obs_eco[1], x_obs_eco[2], x_obs_eco[3]))
```

``` r

# NPE inference with 20,000 simulations
fit_eco <- SBI(
    simulator = ricker_sim,
    prior = ricker_prior,
    x_obs = x_obs_eco,
    method = "NPE",
    n_simulations = 20000,
    n_samples = 5000,
    hidden_layers = c(128, 128),
    n_components = 10,
    max_epochs = 200,
    seed = 42,
    verbose = TRUE
)

print(fit_eco)

# Posterior density with ground truth
plot(fit_eco, type = "posterior",
     ground_truth = c("theta[1]" = true_theta[1],
                      "theta[2]" = true_theta[2]))

# Pairwise posterior showing correlation structure
plot(fit_eco, type = "pairs",
     ground_truth = c("theta[1]" = true_theta[1],
                      "theta[2]" = true_theta[2]))

# Training diagnostics
plot(fit_eco, type = "training")
```

The pairwise plot is particularly informative for ecological models
because population dynamics parameters are often strongly correlated. A
higher growth rate can produce similar mean dynamics to a lower carrying
capacity within certain parameter regimes, creating banana-shaped or
otherwise non-Gaussian posterior structures that test the MDN’s ability
to capture complex geometries. The correlation panel in the upper
triangle of the pairs plot reports the Pearson correlation between the
marginal posteriors.

For comparison with ABC, one would use lucifer’s
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
function on the same model. The key difference in computational cost is
instructive: ABC typically requires \\10^5\\ or more simulations with
careful tuning of the tolerance \\\epsilon\\ and distance metric, while
NPE achieves comparable posterior accuracy with \\2 \times 10^4\\
simulations and no tolerance parameter. When the simulator is expensive
(seconds per call, as in many realistic ecological models), this
five-fold reduction in simulation budget translates directly to
wall-clock savings.

### Summary statistics and the summary_fn argument

In the Ricker example above, the simulator itself returns summary
statistics rather than raw data. However, in many applications it is
cleaner to separate the simulator (which produces raw data) from the
summary statistic extraction (which reduces the data to a
fixed-dimensional vector). The `summary_fn` argument in
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
provides this separation. When `summary_fn` is provided, it is applied
to both the simulator output and the observed data `x_obs` before any
density estimation, ensuring consistency between training and inference.

``` r

# Raw simulator: returns the full time series
ricker_raw <- function(theta) {
    r <- exp(theta[1])
    K <- exp(theta[2])
    sigma <- 0.1

    N <- numeric(100)
    N[1] <- K / 2
    for (t in 2:100) {
        N[t] <- N[t - 1] * exp(r * (1 - N[t - 1] / K) +
                                 rnorm(1, 0, sigma))
        N[t] <- max(N[t], 0.01)
    }
    N
}

# Summary function: extracts informative statistics
ricker_summary <- function(x) {
    c(mean(x), sd(x), cor(x[-1], x[-length(x)]),
      quantile(x, 0.25), quantile(x, 0.75))
}

# Observed raw data
set.seed(42)
x_raw_obs <- ricker_raw(c(log(2), log(150)))

# SBI with summary_fn
fit_summary <- SBI(
    simulator = ricker_raw,
    prior = ricker_prior,
    x_obs = x_raw_obs,
    summary_fn = ricker_summary,
    method = "NPE",
    n_simulations = 20000,
    n_samples = 5000,
    verbose = TRUE
)
```

The choice of summary statistics is one of the most consequential design
decisions in SBI (and ABC). Sufficient statistics preserve all
information about the parameters, but for most simulators no sufficient
statistics are known. In practice, one selects statistics that capture
the key features of the data relevant to the parameters of interest:
location and scale statistics for mean/variance parameters, temporal
correlations for dynamical parameters, and distributional features
(quantiles, skewness, kurtosis) for shape parameters. Including too many
statistics increases the dimensionality of the density estimation
problem and can degrade performance; including too few discards
information and widens the posterior. A useful heuristic is to start
with 3-5 well-chosen statistics and increase the number if the posterior
is too diffuse relative to prior expectations.

When the observed data are low-dimensional (e.g., fewer than 20
dimensions), it is often possible to skip summary statistics entirely
and train the density estimator on the raw data. This avoids the
information loss inherent in summary statistic selection and can yield
tighter posteriors, at the cost of requiring more simulations to train
the higher-dimensional density estimator.

### NLE approach for the ecological model

The Ricker model also illustrates the NLE-MCMC bridge in a realistic
setting. By training a neural likelihood on the summary statistics and
then sampling via MCMC, the user obtains posterior samples with full
MCMC convergence diagnostics, which is particularly valuable for
ecological models where posterior correlations and multimodality are
common:

``` r

# NLE with exact log-prior
ricker_log_prior <- function(theta) {
    dnorm(theta[1], 0.5, 0.5, log = TRUE) +
    dnorm(theta[2], 5, 1, log = TRUE)
}

fit_eco_nle <- SBI(
    simulator = ricker_sim,
    prior = ricker_prior,
    x_obs = x_obs_eco,
    method = "NLE",
    n_simulations = 20000,
    log_prior_fn = ricker_log_prior,
    mcmc_algorithm = "NUTS",
    mcmc_iterations = 10000,
    hidden_layers = c(128, 128),
    n_components = 10,
    verbose = TRUE
)

print(fit_eco_nle)

# MCMC diagnostics are available via the mcmc_fit field
summary(fit_eco_nle)

# Compare NPE and NLE posteriors
cat(sprintf("NPE log(r) mean: %.3f, NLE log(r) mean: %.3f\n",
            mean(fit_eco$Posterior[, 1]),
            mean(fit_eco_nle$Posterior[, 1])))
cat(sprintf("NPE log(K) mean: %.3f, NLE log(K) mean: %.3f\n",
            mean(fit_eco$Posterior[, 2]),
            mean(fit_eco_nle$Posterior[, 2])))

# Pairwise comparison
plot(fit_eco_nle, type = "pairs",
     ground_truth = c("theta[1]" = true_theta[1],
                      "theta[2]" = true_theta[2]))
```

Agreement between the NPE and NLE posteriors on the same problem
provides strong evidence that both approximations are reliable.
Disagreement, on the other hand, warrants investigation: run diagnostics
on both fits, increase the simulation budget, and if the discrepancy
persists, trust the NLE result (which has MCMC convergence guarantees
conditional on the neural likelihood being accurate) over the NPE result
(which has no such guarantees).

This two-method cross-check (NPE + NLE on the same model) is a powerful
practical strategy that lucifer makes uniquely easy. In other SBI
implementations, switching between NPE and NLE requires reconfiguring
the inference pipeline and often changing code structure. In lucifer, it
requires changing only the `method` argument from `"NPE"` to `"NLE"` and
specifying `mcmc_algorithm`. The resulting fits share the same
architecture (MDN backend, same hidden layers and components) and differ
only in the direction of conditioning (posterior vs. likelihood) and the
sampling method (forward pass vs. MCMC), providing an informative
orthogonal check on posterior quality.

### Parallel simulation

When the simulator is expensive and many simulations are needed, the
`n_cores` argument enables parallel simulation via
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html). Each
core independently draws parameters from the prior and runs the
simulator, with failed simulations (those producing non-finite values)
silently dropped from the training set:

``` r

fit_parallel <- SBI(
    simulator = ricker_sim,
    prior = ricker_prior,
    x_obs = x_obs_eco,
    method = "NPE",
    n_simulations = 50000,
    n_samples = 5000,
    n_cores = 4,
    verbose = TRUE
)
```

For sequential simulation (the default `n_cores = 1`), the simulation
loop runs in C++ via `sbi_simulate_batch_cpp()`, which calls the R
simulator function from C++ for each pair and handles failure counting
and progress reporting. The parallel path uses R-level `mclapply()`
instead, which is slightly slower per simulation but scales well across
cores. On a typical 8-core machine, using 4-6 cores provides near-linear
speedup for simulators that take more than a few milliseconds per call;
for very fast simulators, the overhead of process forking may offset the
parallelism benefit.

## Worked example 3: two-moons benchmark

The two-moons problem is a standard SBI benchmark specifically designed
to test density estimators on bimodal posteriors [\[25\]](#ref25). The
simulator generates data from a two-dimensional crescent-shaped
distribution conditioned on a two-dimensional parameter, producing a
posterior with two well-separated modes for certain observed data
values. The name “two moons” comes from the shape of the posterior,
which resembles two crescent moons facing each other in the \\(\theta_1,
\theta_2)\\ plane.

This problem is important for SBI development and evaluation for two
reasons. First, multimodal posteriors arise frequently in scientific
applications: in mixture models, in systems with symmetries (where
permuting component labels produces equivalent solutions), and in
nonlinear models with multiple parameter combinations that produce
similar data. A density estimator that cannot represent multimodality
will produce misleading inferential results, either reporting a single
mode that does not correspond to any true mode (the average of the two)
or capturing only one mode and missing the other. Second, the two-moons
problem has known ground truth properties that enable quantitative
assessment of posterior quality through C2ST and calibration
diagnostics.

``` r

library(lucifer)

# Two-moons simulator
two_moons_sim <- function(theta) {
    a <- -abs(theta[1] + theta[2]) / sqrt(2)
    b <- (-theta[1] + theta[2]) / sqrt(2)
    p <- rnorm(1, mean = a, sd = 0.1)
    ang <- runif(1, -pi/2, pi/2)
    c(p + 0.25 * cos(ang), b + 0.25 * sin(ang))
}

two_moons_prior <- function() runif(2, -1, 1)

# Observed data at origin (produces bimodal posterior)
x_obs_moons <- c(0, 0)
```

The number of mixture components \\K\\ in the MDN is the critical
hyperparameter for representing multimodality. With too few components,
the MDN can only approximate a unimodal distribution and will place a
single broad Gaussian between the two modes; with sufficiently many
components, it can allocate separate Gaussian clusters to each mode. The
following comparison illustrates this effect across three settings.

``` r

# K = 3: likely insufficient for bimodality
fit_k3 <- SBI(
    simulator = two_moons_sim,
    prior = two_moons_prior,
    x_obs = x_obs_moons,
    method = "NPE",
    n_simulations = 20000,
    n_samples = 5000,
    n_components = 3,
    hidden_layers = c(128, 128),
    max_epochs = 300,
    seed = 42,
    verbose = TRUE
)

# K = 10: adequate for most bimodal problems
fit_k10 <- SBI(
    simulator = two_moons_sim,
    prior = two_moons_prior,
    x_obs = x_obs_moons,
    method = "NPE",
    n_simulations = 20000,
    n_samples = 5000,
    n_components = 10,
    hidden_layers = c(128, 128),
    max_epochs = 300,
    seed = 42,
    verbose = TRUE
)

# K = 20: high capacity, may overfit with limited data
fit_k20 <- SBI(
    simulator = two_moons_sim,
    prior = two_moons_prior,
    x_obs = x_obs_moons,
    method = "NPE",
    n_simulations = 20000,
    n_samples = 5000,
    n_components = 20,
    hidden_layers = c(256, 256),
    max_epochs = 300,
    seed = 42,
    verbose = TRUE
)

# Compare pairwise posteriors
plot(fit_k3, type = "pairs")
plot(fit_k10, type = "pairs")
plot(fit_k20, type = "pairs")
```

With \\K = 3\\, the posterior will likely appear as a single broad blob
centered between the two moons, failing to capture the bimodal
structure. At \\K = 10\\, the two crescent-shaped modes should be
clearly resolved, with each mode represented by several Gaussian
components that together approximate the crescent shape. Increasing to
\\K = 20\\ with a larger network can provide finer resolution of the
mode shapes but risks overfitting if the simulation budget is
insufficient, manifesting as noisy, fragmented, or spuriously spiky
posterior samples.

For the Neural Spline Flow backend (available via `backend = "torch"`
and `sbi_network(type = "nsf")`), multimodality is handled more
naturally through the flow’s invertible transformations, which can
deform a unimodal base distribution into complex multi-modal shapes
without needing to choose a number of components. This makes flows
particularly attractive for problems where the number and shape of
posterior modes is unknown a priori.

The two-moons problem also provides an instructive test case for NLE.
Because the posterior is bimodal, MCMC samplers used with the neural
likelihood must be capable of crossing between modes. NUTS and other
gradient-based samplers may get trapped in a single mode, while ensemble
methods (like DREAM or Zeus) or parallel tempering can explore both
modes:

``` r

# NLE with slice sampler on two-moons
fit_moons_nle <- SBI(
    simulator = two_moons_sim,
    prior = two_moons_prior,
    x_obs = x_obs_moons,
    method = "NLE",
    n_simulations = 20000,
    log_prior_fn = function(theta) {
        if (all(theta >= -1 & theta <= 1)) 0 else -Inf
    },
    mcmc_algorithm = "Slice",
    mcmc_iterations = 10000,
    hidden_layers = c(128, 128),
    n_components = 15,
    verbose = TRUE
)

plot(fit_moons_nle, type = "pairs")
```

This example highlights why cross-checking NLE results with NPE is
valuable: if NPE clearly shows two modes but NLE with a particular MCMC
algorithm shows only one, the problem is likely the sampler’s inability
to cross the inter-mode valley rather than the neural likelihood itself.
Trying a different MCMC algorithm or running multiple chains with
dispersed initial values can help resolve this.

## Sequential SBI

When the prior is much broader than the posterior, most simulations fall
in regions of negligible posterior density and contribute little useful
information to the neural network’s training. Sequential methods address
this inefficiency by iterating the simulation-training cycle, using each
round’s posterior estimate to focus the next round’s simulations on the
informative region of parameter space. The improvement can be dramatic:
for a problem where the posterior occupies \\10^{-4}\\ of the prior
volume, sequential SBI can achieve comparable accuracy to single-round
SBI with an order of magnitude fewer simulations.

``` r

library(lucifer)

set.seed(42)
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)
x_obs <- c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0, 2.4, 1.7, 2.2, 2.6)

# SNPE with 3 rounds (9000 total simulations = 3000 per round)
fit_snpe <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "SNPE",
    n_simulations = 9000,
    n_rounds = 3,
    n_samples = 5000,
    hidden_layers = c(128, 128),
    n_components = 10,
    max_epochs = 200,
    seed = 42,
    verbose = TRUE
)

print(fit_snpe)
summary(fit_snpe)

# Analytical posterior for comparison
true_sigma2 <- 1 / (1/25 + 10)
true_mu <- true_sigma2 * sum(x_obs)

plot(fit_snpe, type = "posterior",
     ground_truth = c("theta[1]" = true_mu))
```

The SNPE method automatically handles the multi-round logic. In round 1,
simulations are drawn from the prior. After training, the posterior
samples are used to compute IQR-based truncation bounds
\\\mathcal{R}\_1\\. In round 2, simulations are drawn from the prior
restricted to \\\mathcal{R}\_1\\, the network is retrained on the
accumulated data from both rounds, and new truncation bounds
\\\mathcal{R}\_2\\ are computed. Round 3 proceeds similarly. The total
simulation budget (here 9000) is divided equally across rounds (3000 per
round).

``` r

# Training loss across all rounds
plot(fit_snpe, type = "training")
```

The training loss plot for sequential methods shows the loss trajectory
across all rounds concatenated. A characteristic pattern is a jump in
loss at the beginning of each new round (because new training data from
the focused proposal changes the data distribution) followed by rapid
decrease as the network adapts to the more concentrated data. Each
successive round should achieve a lower minimum loss than the previous
one, reflecting the increasing quality of the training data as the
proposal concentrates near the posterior.

For this conjugate normal-normal example, the improvement from
sequential refinement is modest because the prior standard deviation
(\\\sigma_0 = 5\\) is only about 16 times the posterior standard
deviation (\\\sigma\_{\text{post}} \approx 0.316\\). Sequential methods
show their greatest advantage when the prior-to-posterior contraction
ratio is much larger (e.g., 100:1 or more), which is common in
scientific applications where prior knowledge is vague and the data are
highly informative.

Sequential NLE and SNRE follow the same multi-round logic and are
invoked simply by specifying the appropriate `method`:

``` r

# Sequential NLE with 3 rounds and NUTS sampling
fit_snle <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "SNLE",
    n_simulations = 9000,
    n_rounds = 3,
    log_prior_fn = function(theta) dnorm(theta[1], 0, 5, log = TRUE),
    mcmc_algorithm = "NUTS",
    mcmc_iterations = 5000,
    verbose = TRUE
)

print(fit_snle)
```

For SNLE and SNRE, each round involves both a training phase (fitting
the neural likelihood or ratio) and an MCMC phase (sampling the
posterior under the current neural surrogate), which makes these methods
computationally more expensive per round than SNPE. The posterior
samples from MCMC in round \\r\\ are used to define the proposal for
round \\r + 1\\. The wall-clock time for sequential NLE is roughly \\R\\
times the time for single-round NLE (where \\R\\ is the number of
rounds), because each round requires a full network training and MCMC
sampling phase. For expensive simulators where the simulation time
dominates, the overhead of the extra training and MCMC phases is
negligible, and the efficiency gains from focused proposals can be
dramatic.

The `n_sims_per_round` argument allows explicit control over the
per-round simulation budget. This is useful when the user wants to
front-load simulations in the first round (where the prior is broadest
and more simulations are needed to cover the space) and use fewer
simulations in subsequent rounds (where the proposal is focused and the
network needs less data to refine the density estimate):

``` r

# Custom budget: 5000 in round 1, 3000 in round 2, 2000 in round 3
fit_snpe_custom <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "SNPE",
    n_simulations = 10000,
    n_rounds = 3,
    n_sims_per_round = c(5000L, 3000L, 2000L),
    n_samples = 5000,
    verbose = TRUE
)
```

A useful diagnostic for sequential methods is to compare the posterior
from each round against the final posterior or the ground truth. If the
posterior changes dramatically between rounds 2 and 3, this suggests
that additional rounds might be beneficial; if the posterior has
stabilized, further rounds are unlikely to improve results and the
simulation budget is better spent on increasing per-round simulations.
In lucifer, the `training` field of the returned `sbi` object includes
the training loss across all rounds, and the pattern of loss values
across rounds provides a quick qualitative assessment of convergence: if
the final loss in round \\r + 1\\ is substantially lower than in round
\\r\\, the additional round provided useful information; if the losses
are similar, the proposal has converged and additional rounds are
unnecessary.

The interaction between the number of rounds and the per-round
simulation budget deserves careful thought. For a fixed total budget
\\N\\, increasing the number of rounds \\R\\ reduces the per-round
budget \\N/R\\, which can degrade the quality of each round’s network
training. The optimal balance depends on the prior-to-posterior
contraction ratio: problems with extreme contraction (broad prior, sharp
posterior) benefit from more rounds with fewer simulations per round,
because each round’s focused proposal dramatically improves the
signal-to-noise ratio of the training data. Problems with moderate
contraction may perform better with fewer rounds and more simulations
per round, because each round’s network training benefits from a larger
dataset. As a rough guide: 3 rounds is appropriate for contraction
ratios up to 100:1; 5 rounds may be needed for contraction ratios above
1000:1; and single-round NPE often suffices for contraction ratios below
10:1.

## Diagnostics

SBI methods produce approximate posteriors, and unlike MCMC, where
convergence diagnostics such as \\\hat{R}\\ and effective sample size
have well-established theoretical backing, assessing the quality of a
neural density estimator’s output requires specialized calibration
diagnostics. The fundamental question is whether the learned posterior
\\q\_\phi(\theta \mid x)\\ is a faithful representation of the true
posterior \\p(\theta \mid x)\\, and answering this requires generating
additional synthetic data to test the posterior’s calibration
properties. lucifer provides five complementary diagnostic tools, each
addressing a different aspect of posterior quality.

### Simulation-based calibration (SBC)

Simulation-based calibration [\[6\]](#ref6) is the most principled
general-purpose diagnostic for any Bayesian computational method. The
underlying theorem states that if the posterior is correctly calibrated
with respect to the prior predictive distribution, then for parameters
drawn from the prior and data simulated from the model, the rank of the
true parameter among posterior samples should be uniformly distributed.

Formally, for each calibration replicate \\i = 1, \ldots, M\\:

1.  Draw \\\theta^\*\_i \sim p(\theta)\\ from the prior.
2.  Simulate \\x^\*\_i \sim p(x \mid \theta^\*\_i)\\ from the model.
3.  Draw \\L\\ posterior samples \\\theta_1, \ldots, \theta_L \sim
    q\_\phi(\theta \mid x^\*\_i)\\ from the learned posterior.
4.  Compute the rank statistic \\r_i^{(j)} = \sum\_{l=1}^L
    \mathbf{1}\[\theta_l^{(j)} \< \theta_i^{\*(j)}\]\\ for each
    parameter dimension \\j\\.

Under exact calibration, \\r_i^{(j)} \sim \text{Uniform}\\0, 1, \ldots,
L\\\\ for each parameter dimension independently. This is a consequence
of the probability integral transform: if \\\theta^\*\\ and \\\theta_1,
\ldots, \theta_L\\ are exchangeable draws from the same distribution
(which holds when the posterior is correctly calibrated), the rank of
\\\theta^\*\\ among the \\L + 1\\ values is uniformly distributed.

Deviations from uniformity in the rank histogram reveal specific types
of miscalibration that are diagnostically informative. A U-shaped
histogram (excess counts at both tails) indicates the posterior is too
narrow (overconfident), because the true parameter frequently falls
outside the posterior’s bulk, receiving very low or very high ranks. An
inverted-U shape (excess counts in the middle) indicates the posterior
is too broad (underconfident), because the true parameter is more likely
to fall near the center of the posterior than it would under the true
distribution. Left-skewed histograms (excess counts at low ranks)
indicate that the posterior is systematically shifted to the right of
the true value; right-skewed histograms indicate the opposite. A
Kolmogorov-Smirnov test of the rank distribution against Uniform\\(0,
1)\\ provides a formal quantitative assessment, with \\p \< 0.05\\
providing evidence of miscalibration at the 5% significance level.

``` r

library(lucifer)

# Assume fit_npe and simulator/prior from Example 1 are available
sbc_result <- SBC(
    fit_npe,
    simulator = simulator,
    prior = prior,
    n_sims = 300,
    n_posterior_samples = 200,
    verbose = TRUE
)

# KS test p-value: values > 0.05 are consistent with calibration
cat(sprintf("KS p-value: %.4f\n", sbc_result$ks_pvalues[1]))
cat(sprintf("Interpretation: %s\n",
            ifelse(sbc_result$ks_pvalues[1] > 0.05,
                   "consistent with calibration",
                   "evidence of miscalibration")))

# Rank histogram
plot(fit_npe, type = "sbc", sbc = sbc_result)
```

A well-calibrated SBI posterior produces a flat rank histogram with all
bars near the expected count (shown as a dashed red line). The KS
p-value provides a formal test: \\p \> 0.05\\ is consistent with correct
calibration, though with only 300 replicates the test has limited power
against subtle miscalibration. Increasing `n_sims` to 1000 or more
provides a more sensitive test at the cost of additional computation,
since each replicate requires a fresh simulation and posterior sampling
step.

For NLE and NRE methods, SBC is computationally more expensive because
each replicate requires a short MCMC run to obtain posterior samples.
The internal helper uses a lightweight slice sampler with a short chain
(500 iterations or \\2 \times\\ the requested sample size, whichever is
larger) for these calibration runs, trading MCMC quality for speed. If
SBC reveals calibration problems for NLE/NRE, the cause may be either
the neural network (poor likelihood/ratio approximation) or the MCMC
(insufficient convergence in the short calibration chains), and
disentangling these requires examining the MCMC diagnostics of the full
posterior run.

### Expected coverage

Expected coverage checks whether credible intervals from the posterior
have the correct frequentist coverage probability [\[1\]](#ref1). For
each credibility level \\\alpha \in \\0.1, 0.2, \ldots, 0.9\\\\ and each
calibration replicate \\i\\, the diagnostic checks whether the true
parameter \\\theta_i^\*\\ falls within the \\\alpha\\-level equal-tailed
credible interval computed from the posterior samples. A well-calibrated
posterior produces a coverage curve that lies on the diagonal of the P-P
(probability-probability) plot: the 50% credible interval should contain
the true value 50% of the time, the 90% interval 90% of the time, and so
on.

``` r

cov_result <- expected_coverage(
    fit_npe,
    simulator = simulator,
    prior = prior,
    n_sims = 300,
    n_posterior_samples = 500,
    levels = seq(0.1, 0.9, by = 0.1),
    verbose = TRUE
)

# Coverage at key levels
print(cov_result$coverage)

# P-P plot
plot(fit_npe, type = "coverage", coverage = cov_result)
```

Deviations from the diagonal indicate specific problems. If the curve
lies consistently above the diagonal (empirical coverage exceeds nominal
coverage at every level), the posterior is too broad: the intervals are
wider than they need to be, which is conservative but wastes inferential
power. If the curve lies below the diagonal, the posterior is too narrow
and credible intervals are overconfident, which is the more dangerous
failure mode because it can lead to substantively incorrect conclusions.

A common pattern for undertrained MDNs is a “hockey stick” shape where
coverage is approximately correct at low levels (e.g., the 10% and 20%
intervals have correct coverage) but drops below the diagonal at high
levels (the 80% and 90% intervals have insufficient coverage). This
indicates that the bulk of the posterior is well-calibrated but the
tails are too thin, which is a characteristic artifact of MDN models
that have not seen enough training examples in the posterior tails. The
remedy is to increase the simulation budget, which populates the tails
with more training data, or to increase the number of mixture
components, which provides more flexibility for modeling tail behavior.

Expected coverage and SBC are closely related diagnostics that probe
calibration from different angles: SBC checks the full rank distribution
parameter by parameter, while expected coverage specifically targets the
calibration of interval estimates. Both should be examined together for
a comprehensive assessment. When one passes and the other fails, the
failure mode is informative about which aspect of the posterior is
miscalibrated.

### Classifier two-sample test (C2ST)

The classifier two-sample test [\[14\]](#ref14) evaluates posterior
quality by training a binary classifier to distinguish samples from the
approximate posterior \\q\_\phi(\theta \mid x_o)\\ from samples from a
reference distribution. If the classifier achieves cross-validated
accuracy near 0.5, the two sample sets are statistically
indistinguishable, indicating that the approximate posterior is close to
the reference. Accuracy near 1.0 means the distributions are easily
separable, indicating poor approximation.

The test works by combining \\n\\ samples from the approximate posterior
(labelled 1) with \\n\\ samples from the reference distribution
(labelled 0), then training a binary classifier with \\k\\-fold
cross-validation and reporting the mean accuracy across folds. Under the
null hypothesis that both distributions are identical, the expected
accuracy is 0.5 (random guessing). The test statistic is the accuracy
itself, with values significantly above 0.5 indicating distributional
differences.

In lucifer’s implementation, the classifier is the same C++ NRE binary
classifier used for ratio estimation, repurposed as a two-sample test
statistic. The classifier is trained independently on each fold, and the
per-fold accuracies are reported alongside the mean, providing an
uncertainty estimate for the accuracy. A mean accuracy of 0.50-0.55
indicates no detectable difference; 0.55-0.75 indicates moderate
differences that may be practically acceptable; and above 0.75 indicates
large distributional discrepancies requiring investigation.

``` r

# Compare posterior against analytical posterior (conjugate model)
analytical_samples <- matrix(
    rnorm(5000, mean = true_mu, sd = true_sd),
    ncol = 1
)

c2st_result <- C2ST(
    fit_npe,
    reference_samples = analytical_samples,
    n_folds = 5,
    hidden_layers = c(64, 64),
    max_epochs = 50,
    verbose = TRUE
)

cat(sprintf("C2ST accuracy: %.3f\n", c2st_result$accuracy))

# Visualize per-fold accuracies
plot(fit_npe, type = "c2st", c2st = c2st_result)
```

When comparing against the true analytical posterior, a C2ST accuracy
below 0.55 indicates excellent approximation quality, meaning that a
neural classifier cannot reliably distinguish between samples from the
two distributions. Values between 0.55 and 0.75 suggest moderate
distributional differences that may or may not be practically
significant depending on the application. Values above 0.75 indicate
substantial discrepancies that warrant investigation of the network
architecture, simulation budget, or training procedure.

A simpler but less informative use of C2ST is comparing the posterior
against prior draws:

``` r

# Compare against prior (basic sanity check)
c2st_prior <- C2ST(
    fit_npe,
    prior = prior,
    n_folds = 5,
    verbose = TRUE
)

cat(sprintf("C2ST accuracy vs prior: %.3f\n", c2st_prior$accuracy))
```

When comparing against the prior, accuracy should be high (near 1.0),
because a well-trained SBI posterior that has learned anything from the
data should be easily distinguishable from the prior. An accuracy near
0.5 against the prior would indicate that the SBI procedure has failed
completely, producing samples indistinguishable from prior draws.

### Tests of accuracy with random points (TARP)

TARP [\[13\]](#ref13) is a more powerful calibration diagnostic than SBC
that can detect both global and local miscalibration, including
misspecification of posterior correlations that marginal rank histograms
would miss entirely. The key insight is to test calibration using
distances to randomly chosen reference points in parameter space, rather
than using the coordinate-aligned intervals that SBC relies on.

For each of \\M\\ synthetic datasets, TARP: (i) draws a true parameter
\\\theta^\*\\ from the prior and simulates data \\x^\* \sim p(x \mid
\theta^\*)\\; (ii) draws \\L\\ posterior samples \\\theta_1, \ldots,
\theta_L \sim q\_\phi(\theta \mid x^\*)\\; (iii) picks a random
reference point \\\xi\\ from the prior; (iv) computes the distance
\\d(\theta^\*, \xi)\\ from the true parameter to the reference, and the
distances \\d(\theta_l, \xi)\\ from each posterior sample to the
reference; (v) computes the fraction of posterior samples closer to the
reference than the true parameter: \\f = \frac{1}{L}\sum\_{l=1}^L
\mathbf{1}\[d(\theta_l, \xi) \leq d(\theta^\*, \xi)\]\\.

Under correct calibration, \\f\\ should be uniformly distributed on
\\\[0, 1\]\\, because the true parameter is statistically exchangeable
with the posterior samples. The expected coverage probability (ECP) at
nominal level \\\alpha\\ is \\\text{ECP}(\alpha) =
\frac{1}{M}\sum\_{i=1}^M \mathbf{1}\[f_i \leq \alpha\]\\, which should
equal \\\alpha\\ for all \\\alpha\\.

``` r

tarp_result <- TARP(
    fit_npe,
    simulator = simulator,
    prior = prior,
    n_sims = 300,
    n_posterior_samples = 500,
    n_references = 100,
    verbose = TRUE
)

# TARP diagnostic plot
plot(fit_npe, type = "tarp", tarp = tarp_result)
```

The TARP plot shows the ECP curve (blue line) against the diagonal
(dashed line), with a gray confidence band representing \\\pm 2\\
standard errors under the null hypothesis of correct calibration. The
curve should lie within the band for a well-calibrated posterior.
Systematic deviations above or below the diagonal have the same
interpretation as in expected coverage: conservative (too broad) or
overconfident (too narrow) posteriors, respectively.

TARP’s advantage over SBC is that it tests calibration along random
directions in parameter space rather than just along coordinate axes,
making it sensitive to miscalibration in posterior correlations and
other joint properties. For example, consider a two-parameter posterior
where the marginal distributions for \\\theta_1\\ and \\\theta_2\\ are
each correctly calibrated (SBC passes) but the correlation between them
is wrong (e.g., the true posterior has \\\text{Cor}(\theta_1, \theta_2)
= 0.8\\ but the learned posterior has \\\text{Cor}(\theta_1, \theta_2) =
0\\). SBC, which examines each parameter marginally, would detect no
problem. TARP, which tests distances to random reference points (and
therefore probes the joint distribution along arbitrary directions),
would detect the miscalibration because the true parameter would fall
outside the posterior’s credible regions defined by diagonal directions
that probe the correlation structure.

The `n_references` parameter (default 100) controls how many random
reference points are used. Each calibration replicate uses one randomly
selected reference point from this pool. More reference points provide
better coverage of the parameter space directions but do not increase
the computational cost per replicate. The main computational cost is
determined by `n_sims` (the number of synthetic datasets) and
`n_posterior_samples` (the number of posterior samples per dataset). For
a thorough TARP assessment, 500 or more synthetic datasets are
recommended; 300 provides a reasonable quick check.

Using SBC, expected coverage, and TARP together provides the most
comprehensive calibration assessment, with each diagnostic covering a
different aspect of posterior quality: marginal calibration (SBC),
interval calibration (expected coverage), and joint calibration (TARP).

For practical purposes, the computational cost of diagnostics should be
weighed against their informational value. For a quick initial check,
SBC with 200 replicates and 200 posterior samples per replicate takes a
few minutes for NPE and provides a basic sanity check. For a thorough
assessment before publication, 500-1000 replicates with 500 posterior
samples each provides high-powered tests that can detect subtle
miscalibration. For NLE/NRE, each SBC replicate requires a short MCMC
run (500 iterations by default), so the computational cost scales
linearly with `n_sims`; running 300 SBC replicates may take 15-60
minutes depending on the network evaluation cost and the parameter
dimension. TARP is roughly comparable to SBC in computational cost but
provides complementary information about joint calibration, making it
worth the additional investment for multivariate problems where
correlation structure matters.

The C2ST diagnostic is most useful when a reference posterior is
available (either from an analytical solution or from a different,
trusted inference method). When no reference is available, C2ST can be
run against the prior as a basic sanity check, but this provides limited
information beyond confirming that the posterior differs from the prior.
A creative use of C2ST is to compare the outputs of two different SBI
methods (e.g., NPE versus NLE-NUTS) on the same problem; if C2ST
accuracy is near 0.5, both methods agree, providing mutual validation.

PPC is computationally cheap (it requires only forward simulations from
the fitted posterior, with no additional network evaluations) and should
be run routinely as a model consistency check. Unlike the other four
diagnostics, PPC does not test calibration but rather model adequacy: it
answers the question “can the model with these inferred parameters
reproduce the observed data?” rather than “are the inferred parameters
correctly distributed?” Both questions are important for scientific
inference, and answering both provides a more complete picture of the
analysis quality.

In summary, a minimal diagnostic suite for publication-quality SBI
results should include SBC (marginal calibration), expected coverage or
TARP (interval/joint calibration), and PPC (model consistency). Adding
C2ST when a reference is available provides a quantitative measure of
distributional agreement. The full five-diagnostic suite takes 30-60
minutes for typical NPE problems and provides a comprehensive assessment
that addresses all major potential failure modes of neural density
estimation.

### Posterior predictive checks (PPC)

Posterior predictive checks assess model fit by simulating data from the
fitted posterior and comparing the distribution of simulated data (or
summary statistics) with the observed values [\[18\]](#ref18). Unlike
the calibration diagnostics above, PPC does not require generating new
synthetic datasets from the prior; instead, it uses the posterior that
was already computed for the observed data. If the model is
well-specified and the posterior is accurate, the observed data should
look plausible under the posterior predictive distribution
\\p(x\_{\text{rep}} \mid x_o) = \int p(x\_{\text{rep}} \mid \theta)
q\_\phi(\theta \mid x_o) \\ d\theta\\.

``` r

# Generate posterior predictive samples
ppc_data <- sbi_ppc(
    fit_npe,
    simulator = simulator,
    n_sims = 200
)

# Posterior predictive check plot
plot(fit_npe, type = "ppc", ppc_data = ppc_data)
```

The PPC plot shows histograms of each summary statistic (or data
dimension) computed from the 200 posterior predictive replications. For
a well-specified model with an accurate posterior, the observed value of
each statistic should fall in the bulk of the posterior predictive
distribution, not in the tails. Summary statistics whose observed values
fall below the 2.5th or above the 97.5th percentile of the posterior
predictive distribution suggest potential model misspecification (the
simulator cannot reproduce the observed data features) or poor posterior
approximation (the SBI posterior is placing mass on parameter values
that do not generate data resembling the observation).

PPC is particularly valuable for SBI because it directly tests the
inferential pipeline’s self-consistency: parameters drawn from the
learned posterior should produce data that resembles the observation.
This check does not require knowledge of the true posterior and can be
applied to any problem, making it an essential complement to the
calibration diagnostics that require a simulation-based ground truth.

### Diagnostic interpretation guide

The five diagnostics probe different aspects of posterior quality, and
combining them provides a comprehensive picture. The following decision
logic can guide interpretation.

If SBC rank histograms are flat and KS p-values exceed 0.05, and
expected coverage curves lie on the diagonal, the posterior is
well-calibrated in terms of marginal and interval properties. If TARP
also lies on the diagonal, the joint distribution (including
correlations) is well-calibrated. If C2ST accuracy against a known
reference posterior is below 0.55, the distributional agreement is
strong. If PPC shows that the observed data is consistent with the
posterior predictive distribution, the model is well-specified. When all
five diagnostics pass, the SBI posterior can be trusted for scientific
inference.

When diagnostics fail, the pattern of failure is informative. If SBC
shows a U-shaped rank histogram (overconfidence) while expected coverage
shows the curve below the diagonal, the posterior is too narrow. This
can result from insufficient simulation budget (the network has not seen
enough of the parameter space), an overfitted network (too many
components for the data, or insufficient early stopping patience), or a
posterior that is fundamentally more complex than the MDN can represent.
The remedies are, in order: increase `n_simulations`, reduce
`n_components`, and if neither helps, switch to the torch NSF backend.

If SBC shows an inverted-U histogram (underconfidence) and coverage
curves lie above the diagonal, the posterior is too broad. This
typically indicates that the network has not converged during training
(training loss still decreasing at termination) or that the simulation
budget is so small that the network defaults to a diffuse approximation.
The remedy is to increase `max_epochs` or `n_simulations`.

If SBC passes but TARP fails (curve deviates from diagonal), the
marginal distributions are correct but the joint distribution is wrong,
meaning the posterior correlations are misspecified. This can happen
when `n_components` is too low to capture the correlation structure or
when the parameter space has strong nonlinear dependencies. Increasing
`n_components` or switching to a normalizing flow backend typically
resolves this.

If C2ST accuracy is high against a known posterior but SBC appears flat,
the problem may be subtle distributional differences that SBC’s rank
test lacks power to detect with the given `n_sims`. Increasing the
number of SBC replicates to 1000 or more can help, or the C2ST result
may reflect differences that are statistically detectable but
practically irrelevant.

``` r

# Comprehensive diagnostic workflow
# Step 1: SBC (marginal calibration)
sbc <- SBC(fit_npe, simulator, prior, n_sims = 500,
           n_posterior_samples = 200, verbose = TRUE)
plot(fit_npe, type = "sbc", sbc = sbc)

# Step 2: Expected coverage (interval calibration)
cov <- expected_coverage(fit_npe, simulator, prior, n_sims = 500,
                         n_posterior_samples = 500, verbose = TRUE)
plot(fit_npe, type = "coverage", coverage = cov)

# Step 3: TARP (joint calibration including correlations)
tarp <- TARP(fit_npe, simulator, prior, n_sims = 500,
             n_posterior_samples = 500, verbose = TRUE)
plot(fit_npe, type = "tarp", tarp = tarp)

# Step 4: PPC (model consistency)
ppc <- sbi_ppc(fit_npe, simulator, n_sims = 200)
plot(fit_npe, type = "ppc", ppc_data = ppc)

# Step 5: C2ST against reference (if available)
ref_samples <- matrix(rnorm(5000, mean = true_mu, sd = true_sd), ncol = 1)
c2st <- C2ST(fit_npe, reference_samples = ref_samples, verbose = TRUE)
plot(fit_npe, type = "c2st", c2st = c2st)
```

## SBI versus other lucifer methods

Choosing among lucifer’s inference engines depends on the problem’s
characteristics, and making the right choice can dramatically affect
both the quality of the inference and the computational cost. The key
factors are whether the likelihood is available (which determines
whether SBI or ABC is needed), whether the posterior is expected to be
multimodal (which affects the choice of MDN components or flow
architecture), the dimensionality of the parameter space (which
determines whether the MDN or torch backend is appropriate), the
computational cost of the simulator (which determines whether sequential
methods are worthwhile), and whether amortized inference is needed
(which favors NPE). The following table summarizes the key tradeoffs
across the five main engines.

The decision about which engine to use should be guided by a few key
questions, considered in order of priority.

The first and most decisive question is whether the likelihood is
evaluable. If \\\log p(x \mid \theta)\\ can be computed at any parameter
value, MCMC is the gold standard for moderate-dimensional problems (up
to roughly 100 parameters with efficient algorithms like NUTS), offering
exact posterior samples with well-understood convergence properties.
lucifer provides 62 MCMC algorithms spanning random-walk Metropolis,
adaptive methods, gradient-based samplers, slice samplers, ensemble
methods, and more. Variational Bayes provides a faster approximate
alternative when the parameter dimension is very high or exact inference
is unnecessary, and SMC is preferred when the posterior is multimodal or
sequential data assimilation is involved. There is no reason to use SBI
when the likelihood is available, because SBI adds neural network
approximation error that is absent from direct MCMC.

It is worth emphasizing when SBI should not be used, because the
availability of neural density estimators can create a temptation to
apply them even when simpler methods suffice. If the likelihood is
tractable and the parameter dimension is moderate (below 50), MCMC with
NUTS or an adaptive Metropolis method will produce exact posterior
samples with well-characterized convergence properties, and no amount of
neural network sophistication can improve upon exact inference. If the
likelihood is intractable but the parameter dimension is very low (1-2
parameters) and the simulator is fast, ABC-SMC may be simpler to
implement and diagnose than SBI, with fewer hyperparameters to tune. SBI
is most valuable in the regime where the likelihood is intractable, the
parameter dimension is moderate to high (3-50), the simulator is
expensive enough that simulation efficiency matters, or amortized
inference is desired for repeated application to multiple datasets.
Another important use case for SBI is model development: when a
researcher is iteratively refining a simulator model, SBI provides fast
approximate posterior feedback that can guide model development
decisions (e.g., which model parameters are identifiable from the
available data) before committing to the expense of a full MCMC analysis
with a computationally cheaper approximate likelihood.

If the likelihood is not evaluable, the choice narrows to ABC or SBI.
ABC is appropriate for very low-dimensional problems (fewer than roughly
5 summary statistics) where the simulator is fast, a good distance
metric is available, and the user is willing to tune the tolerance
parameter \\\epsilon\\. lucifer’s
[`ABC()`](https://robustecologies.github.io/lucifer/reference/ABC.md)
function provides rejection, MCMC, and SMC-ABC variants with automatic
tolerance selection. For higher-dimensional problems, when the simulator
is expensive and simulation efficiency matters, or when amortized
inference is desired, SBI is preferred.

Within SBI, NPE is the default recommendation: it is the simplest to use
(no MCMC tuning required), provides amortized inference (once trained,
posterior evaluation for new \\x_o\\ is essentially free), and handles
most problems well with the MDN backend. NLE becomes the preferred
method when the user wants to leverage lucifer’s MCMC infrastructure for
posterior exploration, including convergence diagnostics, effective
sample size estimation, and the ability to cross-check results across
multiple MCMC algorithms. NLE is also preferred when the likelihood
factorizes over i.i.d. observations. NRE is useful when the
classification formulation is natural or when one wants the simplest
possible training objective with no density estimation machinery.

Sequential methods (SNPE, SNLE, SNRE) should be used when the prior is
much broader than the posterior and the simulation budget is limited. A
useful heuristic: if the prior standard deviation for each parameter is
more than 10 times the expected posterior standard deviation, sequential
methods will provide substantial efficiency gains.

The following NPE-vs-NLE comparison on the same problem illustrates the
key tradeoff between amortized inference (NPE) and MCMC-backed inference
(NLE):

``` r

library(lucifer)

set.seed(42)
simulator <- function(theta) rnorm(10, mean = theta[1], sd = 1)
prior <- function() rnorm(1, mean = 0, sd = 5)
x_obs <- c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0, 2.4, 1.7, 2.2, 2.6)

# NPE: fast, amortized, no MCMC diagnostics
fit_npe <- SBI(simulator, prior, x_obs, method = "NPE",
               n_simulations = 10000, n_samples = 5000, seed = 42)

# NLE + NUTS: slower, but provides full MCMC diagnostics
fit_nle <- SBI(simulator, prior, x_obs, method = "NLE",
               n_simulations = 10000, mcmc_algorithm = "NUTS",
               mcmc_iterations = 5000,
               log_prior_fn = function(theta) dnorm(theta[1], 0, 5, log = TRUE),
               seed = 42)

# Compare
cat(sprintf("NPE mean: %.4f (%.1f min)\n",
            mean(fit_npe$Posterior[, 1]), fit_npe$Minutes))
cat(sprintf("NLE mean: %.4f (%.1f min)\n",
            mean(fit_nle$Posterior[, 1]), fit_nle$Minutes))

# NLE provides MCMC diagnostics that NPE cannot
if (!is.null(fit_nle$mcmc_fit)) {
    cat("NLE MCMC algorithm:", fit_nle$mcmc_fit$Algorithm, "\n")
}
```

NPE is typically faster because it avoids the MCMC phase entirely, but
NLE provides richer post-hoc analysis through the full `demonoid`
object. For problems where MCMC convergence is a concern or where the
user needs effective sample sizes and convergence diagnostics, NLE’s
additional cost is justified.

To provide concrete timing expectations: for a 1-dimensional conjugate
normal problem with 10000 simulations, NPE completes in approximately
5-15 seconds (simulation + training + sampling), while NLE with 5000
NUTS iterations takes approximately 20-60 seconds (simulation +
training + MCMC). The simulation phase accounts for most of the time
when the simulator is expensive; the network training adds 1-5 seconds,
and the MCMC phase adds 5-30 seconds depending on the algorithm and
number of iterations. For problems with expensive simulators (seconds to
minutes per call), the simulation phase dominates entirely, and the
choice between NPE and NLE has negligible impact on total wall-clock
time.

## Tuning guidance

### Simulation budget

The number of simulations (`n_simulations`) is the single most important
hyperparameter and should be the first thing to increase when
diagnostics reveal problems. As a starting point, 5000 simulations are
sufficient for low-dimensional problems (1-3 parameters) with smooth,
unimodal posterior geometry. For moderate-dimensional problems (4-10
parameters), 10000 to 50000 simulations are typically needed. Problems
with more than 10 parameters, complex posterior structures, or highly
nonlinear simulators may require 100000 or more simulations. The
diagnostics described above (particularly SBC and expected coverage)
should be used to assess whether the simulation budget is sufficient: if
coverage curves deviate systematically from the diagonal, increasing the
simulation budget is the first intervention to try before adjusting
network architecture.

### Network architecture

The `n_components` parameter controls the MDN’s representational
capacity for multimodality and non-Gaussianity. For unimodal,
approximately Gaussian posteriors, \\K = 5\\ is often sufficient and
trains quickly. For posteriors with mild skewness or heavy tails, \\K =
10\\ provides adequate flexibility. For multimodal posteriors or
posteriors with complex geometries, \\K = 15\\ to \\K = 20\\ is
recommended. Setting \\K\\ too high wastes parameters (the output layer
has \\K(1 + d + d(d+1)/2)\\ units) and can lead to overfitting,
particularly with small simulation budgets, manifesting as noisy,
fragmented, or spuriously spiky posterior samples.

The `hidden_layers` parameter controls the network depth and width. The
default `c(128, 128)` (two hidden layers of 128 units each) handles most
problems up to moderate dimension and complexity. For high-dimensional
problems (10+ parameters) or complex nonlinear relationships between
parameters and data, `c(256, 256)` or `c(256, 256, 128)` may be needed.
Deeper networks (3+ hidden layers) rarely help for typical SBI problems
and increase the risk of training instability.

### Learning rate and training schedule

The default learning rate of \\5 \times 10^{-4}\\ for the Adam optimizer
works well for most problems. Increasing to \\10^{-3}\\ can speed
convergence for simple problems but risks training instability for
complex ones, particularly when the Cholesky parametrization produces
near-singular covariance matrices early in training. Decreasing to
\\10^{-4}\\ provides more stable training at the cost of more epochs to
converge.

The `max_epochs` parameter (default 200) sets an upper bound on training
duration, but early stopping with `patience = 20` typically terminates
training well before this limit. If the training loss is still
decreasing noticeably when early stopping triggers (the validation loss
has plateaued but the training loss continues to improve), this suggests
the network is beginning to overfit and the current architecture is
appropriate for the data size. If both training and validation loss are
still decreasing when `max_epochs` is reached, increase the epoch limit
or the patience parameter.

The `max_grad_norm` parameter (default 10.0) controls gradient clipping
and rarely needs adjustment. Reducing it to 5.0 or 1.0 can stabilize
training for pathological problems where the Cholesky factors produce
extreme gradient spikes, but this comes at the cost of slower
convergence. The `batch_size` parameter (default 256) can be reduced to
64 or 128 for small training sets (fewer than 5000 simulations) to
increase the number of gradient updates per epoch. Conversely, for very
large training sets (100000+ simulations), increasing the batch size to
512 or 1024 can speed up training by reducing the number of iterations
per epoch, at the cost of slightly noisier gradient estimates.

The `validation_fraction` parameter (default 0.1) controls what fraction
of the simulations is held out for early stopping monitoring. For small
simulation budgets (fewer than 5000), the default 10% holdout may be too
small to provide a reliable validation signal; consider increasing to
0.15 or 0.2. For large budgets (50000+), the default is more than
adequate, and reducing to 0.05 frees more data for training without
appreciably affecting the early stopping criterion.

The `activation` parameter currently supports only `"selu"`, which is
the recommended activation for the MDN architecture due to its
self-normalizing properties [\[20\]](#ref20). The SELU activation
ensures that the mean and variance of the activations remain stable
across layers without batch normalization, which is particularly
important for the MDN because the output layer must produce well-scaled
Cholesky factors and mixing weights. Alternative activations (ReLU, ELU,
tanh) can work in principle but require more careful initialization and
may need batch normalization layers to prevent activation drift, adding
complexity without clear benefit. SELU has proven robust across a wide
range of SBI problems in practice, from 1-dimensional conjugate models
to 50-dimensional ecological simulators, and there is rarely a reason to
change it.

For users who wish to explore the effect of different training
hyperparameters without modifying source code, the most impactful
parameters to tune (in order of typical importance) are: `n_simulations`
(simulation budget), `n_components` (MDN capacity), `hidden_layers`
(network depth/width), `learning_rate` (optimization step size),
`patience` (early stopping sensitivity), and `max_grad_norm` (gradient
clipping threshold). The default values have been chosen to work well
across a range of problems, but every problem has unique characteristics
that may benefit from tuning.

### When to use sequential methods

Sequential methods are most beneficial when the ratio of prior volume to
posterior volume is large. A useful heuristic: if the prior standard
deviation for each parameter is more than 10 times the expected
posterior standard deviation, sequential methods will provide
substantial efficiency gains over single-round methods. For the typical
scientific application where prior knowledge is vague and the data are
informative, this condition is often met.

The default of 3 rounds with equal simulation budgets per round is a
reasonable starting point for most problems. Increasing to 5 rounds
provides diminishing returns for most problems, as the proposal
converges to the posterior within 2-3 rounds for well-behaved problems.
However, for highly peaked posteriors or problems with strong parameter
correlations, 5 rounds can provide meaningful improvement.

### When to switch to the torch backend

The C++ MDN backend is preferable for most problems up to roughly 20
parameters because it requires no additional dependencies, compiles with
the package, and trains quickly. Consider switching to the torch backend
with Neural Spline Flows when any of the following conditions hold: the
MDN produces poor calibration diagnostics despite adequate simulation
budgets and generous \\K\\, suggesting that the mixture-of-Gaussians
assumption is too restrictive for the posterior geometry; the parameter
dimension exceeds 20, where the MDN’s output layer becomes prohibitively
large (the number of Cholesky parameters grows as \\O(Kd^2)\\); or the
posterior has complex topology (rings, manifold structure, sharp ridges)
that Gaussian mixtures cannot represent efficiently regardless of \\K\\.

### Recommended starting configurations

The following configurations provide reasonable defaults for common
problem categories.

These are starting points, not guarantees. Always run at least SBC and
expected coverage diagnostics after fitting and adjust hyperparameters
if calibration is poor. When in doubt, increase the simulation budget
before adjusting network architecture: insufficient training data is the
most common cause of poor SBI performance.

### Common failure modes

Training loss plateaus at a high value without decreasing indicate that
the network is too small (increase `hidden_layers`) or the learning rate
is too high (decrease `learning_rate`). If only the validation loss
plateaus while the training loss continues to decrease, the network is
overfitting: reduce `n_components` or `hidden_layers`, or increase
`n_simulations`. Posterior samples that cluster at the prior boundaries
suggest that the MDN is placing substantial mass outside the prior’s
domain; the rejection sampling mechanism handles mild cases
automatically, but severe boundary effects require either more
simulations, a reparametrization that moves the posterior away from the
prior boundaries, or the use of `log_prior_fn` to enforce exact prior
bounds. Very high C2ST accuracy when comparing against a known reference
posterior indicates systematic posterior misspecification that requires
investigation of the network architecture, simulation budget, or summary
statistics.

A particularly insidious failure mode is the “oversmooth posterior,”
where the MDN converges to a smooth but overly broad approximation that
passes SBC (because it is conservative) but has poor practical utility.
This manifests as expected coverage curves above the diagonal at all
levels. The root cause is usually insufficient network capacity to
represent the posterior’s structure: the MDN defaults to a wide Gaussian
because it cannot learn the sharp features of the true posterior. The
remedy is to increase both `n_components` and `hidden_layers`, or switch
to the torch NSF backend.

Another common problem is “mode collapse” in multimodal posteriors,
where the MDN captures only one mode and ignores the other. This is
fundamentally a local-optimum problem in the network training: the Adam
optimizer converges to a loss minimum that represents one mode well but
ignores the other. Increasing `n_components` helps by providing more
capacity, but the most reliable solution is to increase `n_simulations`
so that both modes are well-represented in the training data, and to
consider using a normalizing flow backend that avoids the mixture-model
formulation entirely.

For NLE and NRE methods, additional failure modes relate to the MCMC
phase. If the neural likelihood surface has sharp ridges or
discontinuities (which can happen when the network is poorly trained),
gradient-based samplers like NUTS may struggle with divergent
transitions, and non-gradient samplers like Slice or random-walk
Metropolis may be more appropriate. If the MCMC chains show poor mixing
(low effective sample size, high autocorrelation), the neural likelihood
surface may be poorly conditioned for the chosen sampler. Trying
multiple algorithms (which the bridge architecture makes trivial) helps
diagnose whether the issue is the neural likelihood itself or the
sampler-likelihood interaction.

## Practical workflow

Bringing together the concepts from the preceding sections, the
following workflow provides a systematic approach to applying SBI to a
new problem.

The first step is to define and test the simulator. The simulator
function must have the signature `function(theta)` where `theta` is a
numeric vector, and must return a numeric vector of fixed length for
every valid parameter draw from the prior. Before running SBI, verify
that the simulator produces finite, non-degenerate output for a range of
parameter values drawn from the prior. Simulators that crash, produce
`NaN` or `Inf` values, or return vectors of varying length will cause
problems that are difficult to diagnose during training. Common pitfalls
include: division by zero for certain parameter values (add guards or
reparametrize); numerical overflow in exponentials (use log-scale
computation); stochastic processes that can produce zero or negative
populations (add floor guards); and simulators that take excessively
long for certain “pathological” parameter values (add timeout guards). A
useful sanity check is to run 100 prior-predictive simulations and
examine the distribution of outputs:

``` r

# Test simulator robustness
set.seed(42)
n_test <- 100
test_outputs <- vector("list", n_test)
for (i in seq_len(n_test)) {
    theta_test <- prior()
    test_outputs[[i]] <- simulator(theta_test)
}

# Check for failures
n_finite <- sum(vapply(test_outputs, function(x) all(is.finite(x)), logical(1)))
cat(sprintf("%d/%d simulations produced finite output\n", n_finite, n_test))

# Check output dimensionality consistency
dims <- vapply(test_outputs, length, integer(1))
cat(sprintf("Output dimensions: %s\n",
            paste(unique(dims), collapse = ", ")))
```

The second step is to choose summary statistics (if applicable). If the
simulator output is low-dimensional (fewer than roughly 20 dimensions),
try training on the raw output first, as this avoids the information
loss inherent in dimensionality reduction. If the output is
high-dimensional (time series, spatial data, images), choose a set of
3-10 informative summary statistics that capture the features of the
data most relevant to the parameters. Good summary statistics should
vary with the parameters of interest (sensitivity) and vary little with
other sources of stochasticity (stability). Common choices include
moments (mean, variance, skewness), quantiles, autocorrelation functions
for temporal data, spatial correlation functions for spatial data, and
regression coefficients for relational data. The `summary_fn` argument
separates this choice from the simulator definition, making it easy to
experiment with different summary statistic sets and assess their impact
on posterior width and calibration.

The third step is to run a pilot fit with NPE using a moderate
simulation budget (5000-10000 simulations) and examine the training
curves and posterior samples. The training curve should show decreasing
loss on both training and validation sets; if the validation loss never
decreases, the network may be too large or the learning rate too high.
If the training loss plateaus at a high value, the network may be too
small.

``` r

# Pilot fit
fit_pilot <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NPE",
    n_simulations = 5000,
    n_samples = 2000,
    n_components = 10,
    hidden_layers = c(128, 128),
    seed = 42,
    verbose = TRUE
)

# Examine training curve
plot(fit_pilot, type = "training")

# Examine posterior
plot(fit_pilot, type = "posterior")
```

The fourth step is to run calibration diagnostics. SBC and expected
coverage are the minimum; TARP provides additional joint calibration
information. If diagnostics reveal miscalibration, increase
`n_simulations` (first priority), adjust `n_components` or
`hidden_layers` (second priority), or try sequential methods (third
priority).

``` r

# Run SBC on pilot fit
sbc <- SBC(fit_pilot, simulator, prior, n_sims = 300, verbose = TRUE)
plot(fit_pilot, type = "sbc", sbc = sbc)

# If KS p-values are good, increase simulation budget for final fit
if (all(sbc$ks_pvalues > 0.05)) {
    cat("Pilot calibration looks good. Running final fit.\n")
    fit_final <- SBI(
        simulator = simulator,
        prior = prior,
        x_obs = x_obs,
        method = "NPE",
        n_simulations = 20000,
        n_samples = 10000,
        n_components = 10,
        seed = 42,
        verbose = TRUE
    )
} else {
    cat("Calibration issues detected. Consider adjusting hyperparameters.\n")
}
```

The fifth step, applicable to NLE/NRE methods, is to verify MCMC
convergence by examining the `mcmc_fit` component of the returned
object. Trace plots, effective sample sizes, and split-\\\hat{R}\\
values provide standard MCMC diagnostics. Running the same neural
likelihood with two or three different MCMC algorithms and comparing
posterior summaries provides an additional cross-check. This
multi-algorithm comparison is a powerful diagnostic that is unique to
the NLE/NRE-MCMC bridge: if two different MCMC algorithms produce
different posteriors from the same neural likelihood, the problem is
MCMC convergence (try longer chains or different algorithms); if they
agree but disagree with NPE, the problem is the neural density estimator
(try more simulations or different architecture).

The sixth and final step is to report results with appropriate
uncertainty quantification. The posterior samples in `fit$Posterior` can
be used to compute posterior means, medians, credible intervals, and
derived quantities. The
[`summary()`](https://rdrr.io/r/base/summary.html) method provides a
formatted summary table, and the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method with
`type = "posterior"` or `type = "pairs"` provides publication-quality
visualizations. When publishing results, it is good practice to report
the SBI method used, the simulation budget, the network architecture,
and the diagnostic results (at minimum SBC or expected coverage), so
that readers can assess the reliability of the inference.

### Troubleshooting checklist

When SBI produces unsatisfactory results, the following checklist helps
identify and resolve the problem.

If the training loss does not decrease at all, check that the simulator
produces valid output for prior draws, that the observed data `x_obs`
has the same dimensionality as the simulator output, and that the data
are not constant (zero variance in the simulations would make
standardization degenerate). If all simulations produce identical
output, the prior may be too narrow or the simulator may be
deterministic in a way that is inconsistent with the parameter values.

If the training loss decreases but the posterior is clearly wrong (e.g.,
centered far from the truth in a conjugate model), check the
standardization by examining `fit$standardization`: the means and
standard deviations should be finite and non-zero. Also verify that
`x_obs` matches the simulator output convention; a common mistake is to
pass raw data when the simulator returns summary statistics, or vice
versa.

If the posterior is approximately correct but calibration diagnostics
show systematic bias, the simulation budget is likely insufficient.
Double the `n_simulations` parameter and re-run. If calibration
improves, continue doubling until the diagnostics pass. If calibration
does not improve with additional simulations, the network architecture
may be insufficient: try increasing `n_components` or `hidden_layers`.

If NLE/NRE produces poor MCMC mixing (very low effective sample sizes or
chains that get stuck), the neural likelihood/ratio surface may have
pathological features. Try a different MCMC algorithm (e.g., switch from
NUTS to Slice), increase the `mcmc_iterations`, or retrain the neural
network with a smoother architecture (fewer components, wider hidden
layers). The `log_prior_fn` argument is critical for NLE/NRE: always
provide it when an analytical prior is available, as the default
approximate prior can introduce artifacts near the boundaries of the
prior support.

## Advanced: the torch backend

For problems that exceed the MDN’s representational capacity, lucifer
provides an optional Neural Spline Flow (NSF) backend [\[8\]](#ref8)
[\[21\]](#ref21) through the torch R package. NSFs use invertible neural
network layers based on rational-quadratic spline transformations that
can represent complex, multimodal distributions without the parametric
assumptions of Gaussian mixtures. The key idea is to transform a simple
base distribution (typically a standard Gaussian) through a sequence of
learned invertible transformations, each of which can deform the density
in flexible ways while maintaining exact likelihood computation through
the change-of-variables formula.

The torch R package provides native bindings to libtorch with no Python
dependency, making it a self-contained solution for users who need
flow-based density estimation within R.

``` r

# Install torch if needed
# install.packages("torch")
# torch::install_torch()

library(lucifer)

# NPE with Neural Spline Flow on the two-moons problem
fit_nsf <- SBI(
    simulator = two_moons_sim,
    prior = two_moons_prior,
    x_obs = x_obs_moons,
    method = "NPE",
    backend = "torch",
    network = sbi_network(type = "nsf", hidden_layers = c(128, 128)),
    n_simulations = 20000,
    n_samples = 5000,
    max_epochs = 300,
    seed = 42,
    verbose = TRUE
)

print(fit_nsf)
plot(fit_nsf, type = "pairs")
```

The NSF backend automatically handles the flow architecture
construction, including the alternating coupling layers, the spline bin
parameters, and the base distribution. Training uses the same Adam
optimizer and early stopping mechanism as the MDN backend, and the same
standardization pipeline applies. The key advantage of NSF over MDN is
its ability to represent complex, non-Gaussian posterior geometries
without specifying a number of mixture components. The flow transforms a
simple base distribution (standard Gaussian) through a series of learned
invertible maps, and the resulting density can take essentially any
shape, including sharp modes, heavy tails, and complex multimodal
structures. The cost is longer training time (flows have more parameters
per layer than MDNs) and the requirement for the torch package, which
adds a substantial dependency (~500 MB for libtorch).

A pretrained network from a previous
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
call can be passed to a new call via the `pretrained_network` argument,
enabling warm-start training when the user wants to refine a fit with
additional simulations or adjusted hyperparameters:

``` r

# First fit with limited simulations
fit1 <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NPE",
    n_simulations = 5000,
    n_samples = 2000,
    verbose = TRUE
)

# Refine with more simulations, starting from the trained network
fit2 <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NPE",
    n_simulations = 10000,
    n_samples = 5000,
    pretrained_network = fit1$network,
    verbose = TRUE
)
```

This warm-start approach is useful when the initial simulation budget
was insufficient and the user wants to generate additional simulations
without discarding the original training. The network begins from a
reasonable initialization (the previous fit’s weights) rather than
random weights, which typically accelerates convergence on the augmented
dataset. Note that the new training data is generated from the prior
independently of the previous fit, so the total training set consists of
the original simulations plus the new ones, all from the prior. The
standardization statistics are recomputed from the new training data, so
the warm-started network must adapt to potentially different
standardization parameters.

The `pretrained_network` argument also enables a useful workflow for
hyperparameter exploration: train an initial model with conservative
settings, then experiment with different `n_components` or
`hidden_layers` while keeping the simulation budget fixed, by passing
the same simulation data (generated internally by the initial call) to
multiple training runs. While
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
does not currently expose a “reuse simulations” interface directly, the
pretrained network mechanism provides an indirect way to iterate on the
training procedure without regenerating simulations.

For both the MDN and torch backends, the trained network object is a
self-contained data structure that can be saved to disk via
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and reloaded for
later use. This is particularly valuable for expensive simulations:
train the network once, save it, and use it for posterior evaluation,
diagnostics, and sensitivity analysis across multiple R sessions without
retraining.

``` r

# Save the entire SBI fit for later use
saveRDS(fit_npe, file = "sbi_fit_npe.rds")

# Reload in a new session
fit_npe <- readRDS("sbi_fit_npe.rds")

# The network is immediately usable for posterior sampling,
# diagnostics, and as pretrained_network for warm-start training
print(fit_npe)
plot(fit_npe, type = "posterior")

# Pass to a new SBI call as pretrained_network
fit_refined <- SBI(
    simulator = simulator,
    prior = prior,
    x_obs = x_obs,
    method = "NPE",
    n_simulations = 10000,
    pretrained_network = fit_npe$network,
    verbose = TRUE
)
```

## Connection to lucifer’s broader ecosystem

One of the key advantages of implementing SBI within lucifer rather than
as a standalone package is the seamless integration with the rest of the
Bayesian inference ecosystem. Several connections deserve explicit
mention.

For NLE and NRE fits, the `mcmc_fit` component of the returned `sbi`
object is a standard `demonoid` object, which means all of lucifer’s
post-processing functions apply directly. The
[`Consort()`](https://robustecologies.github.io/lucifer/reference/Consort.md)
function can provide algorithm recommendations based on the posterior
geometry. The
[`Rhat()`](https://robustecologies.github.io/lucifer/reference/Rhat.md)
and
[`ESS()`](https://robustecologies.github.io/lucifer/reference/ESS.md)
functions compute convergence diagnostics. The
[`MCSE()`](https://robustecologies.github.io/lucifer/reference/MCSE.md)
function estimates Monte Carlo standard errors for posterior summaries.
If the user runs multiple MCMC chains on the same neural likelihood (by
calling
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
multiple times with the same `pretrained_network` and different seeds),
the chains can be combined and analyzed using split-\\\hat{R}\\ for a
rigorous convergence assessment.

The
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md)
function, which benchmarks multiple MCMC algorithms on a given model,
can also be applied to neural likelihood models. By extracting the model
closure from an NLE fit and passing it to
[`Crucible()`](https://robustecologies.github.io/lucifer/reference/Crucible.md),
the user can systematically compare how different MCMC algorithms
perform on the learned neural likelihood surface. This is particularly
useful for identifying the optimal sampler for a given problem
structure.

For model comparison, the
[`LOO()`](https://robustecologies.github.io/lucifer/reference/LOO.md)
function can be applied to NLE/NRE fits if the neural likelihood model
is constructed to return pointwise log-likelihood values. While this
requires some care (the neural likelihood is an approximation, so LOO-CV
is comparing approximate likelihoods), it provides a principled way to
compare different model structures when the likelihood is intractable.

The SBI module also interacts naturally with lucifer’s
[`RobustBayes()`](https://robustecologies.github.io/lucifer/reference/RobustBayes.md)
framework for sensitivity analysis. By running SBI with different prior
specifications and comparing the resulting posteriors, the user can
assess prior sensitivity without requiring a tractable likelihood. This
is valuable in scientific applications where the prior is uncertain and
robustness to prior choice is important for the credibility of the
conclusions.

``` r

# Example: using lucifer's diagnostics on an NLE fit
fit_nle <- SBI(simulator, prior, x_obs, method = "NLE",
               n_simulations = 10000, mcmc_algorithm = "NUTS",
               mcmc_iterations = 5000, verbose = TRUE)

# Standard MCMC diagnostics apply to the mcmc_fit component
if (!is.null(fit_nle$mcmc_fit)) {
    # Effective sample size
    post <- fit_nle$mcmc_fit$Posterior1
    for (j in seq_len(ncol(post))) {
        cat(sprintf("ESS for theta[%d]: %.0f\n", j,
                    coda::effectiveSize(post[, j])))
    }
}
```

### Reproducibility

The `seed` argument in
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
sets R’s random number generator state before any simulations or network
training begins, ensuring that the results are reproducible across runs.
However, there are two caveats. First, the C++ MDN training involves
floating-point operations whose order may vary slightly across platforms
(due to differences in BLAS implementations or compiler optimizations),
which can cause small numerical differences in the trained network
weights and therefore in the posterior samples, even with the same seed.
These differences are within the statistical noise of the posterior
estimate and do not affect scientific conclusions. Second, when using
`n_cores > 1` for parallel simulation, the order of simulations is not
deterministic across runs (it depends on the OS scheduler), so the
training data will differ even with the same seed. For fully
reproducible results, use `n_cores = 1`.

When reporting SBI results in publications, include the seed value, the
simulation budget, the network architecture (hidden layers and
components), the SBI method used (NPE/NLE/NRE and whether sequential),
the MCMC algorithm (for NLE/NRE), and the diagnostic results (at minimum
SBC p-values or expected coverage curves), so that readers and reviewers
can reproduce the analysis and assess its reliability. For ecological
and biological applications, it is also important to describe the
summary statistics used and justify their choice in terms of the model
parameters of interest.

For NLE/NRE fits, additionally report the MCMC iterations, burn-in,
thinning, and effective sample sizes, following the same standards as
for conventional MCMC analyses. The availability of these standard
diagnostics through lucifer’s MCMC bridge is one of the key advantages
of the NLE/NRE approach: it provides reviewers and readers with the same
convergence assurances they expect from traditional Bayesian analyses,
even when the underlying likelihood is neural-network-based.

## API quick reference

This section provides a concise overview of the SBI module’s public
functions and their key arguments for quick lookup.

### Core inference

The
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
function is the main entry point. It accepts a `simulator` function
(taking a parameter vector, returning a data vector), a `prior` function
(taking no arguments, returning a parameter draw), and the observed data
`x_obs`. The `method` argument selects among `"NPE"`, `"NLE"`, `"NRE"`,
`"SNPE"`, `"SNLE"`, and `"SNRE"`. For NPE, the `n_samples` argument
controls how many posterior samples are drawn via the forward pass; for
NLE/NRE, `mcmc_algorithm` and `mcmc_iterations` control the MCMC phase.
The `log_prior_fn` argument (strongly recommended for NLE/NRE) provides
an exact log-prior density evaluation. The `seed` argument sets R’s RNG
state for reproducibility. The `n_cores` argument enables parallel
simulation via `mclapply()`. The `pretrained_network` argument accepts a
network from a previous
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
call for warm-start training.

### Network specification

The
[`sbi_network()`](https://robustecologies.github.io/lucifer/reference/sbi_network.md)
function creates a network specification with `type` (`"mdn"` or
`"nsf"`), `hidden_layers` (integer vector of layer widths),
`n_components` (number of MDN mixture components), and `activation`
(currently `"selu"` only). This specification is passed to
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
via the `network` argument, or the layer/component settings can be
passed directly to
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
via the `hidden_layers` and `n_components` arguments for convenience.

### Diagnostics

Five diagnostic functions assess posterior quality.
[`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md)
performs simulation-based calibration, producing rank histograms and KS
p-values.
[`expected_coverage()`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md)
computes empirical coverage at multiple credibility levels, producing
P-P plots.
[`C2ST()`](https://robustecologies.github.io/lucifer/reference/C2ST.md)
trains a classifier to distinguish posterior samples from reference
samples, reporting cross-validated accuracy.
[`TARP()`](https://robustecologies.github.io/lucifer/reference/TARP.md)
tests calibration along random directions in parameter space, detecting
joint miscalibration that SBC misses.
[`sbi_ppc()`](https://robustecologies.github.io/lucifer/reference/sbi_ppc.md)
generates posterior predictive replications for model checking.

### S3 methods

The `sbi` class supports three S3 methods.
[`print.sbi()`](https://robustecologies.github.io/lucifer/reference/print.sbi.md)
provides a one-screen summary of the fit.
[`summary.sbi()`](https://robustecologies.github.io/lucifer/reference/summary.sbi.md)
provides extended diagnostics including network architecture, training
progression, and MCMC diagnostics for NLE/NRE.
[`plot.sbi()`](https://robustecologies.github.io/lucifer/reference/plot.sbi.md)
supports eight plot types via the `type` argument: `"posterior"`
(default, marginal densities), `"pairs"` (pairwise scatter/density
matrix), `"training"` (loss curves), `"sbc"` (rank histograms from
[`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md)),
`"coverage"` (P-P plot from
[`expected_coverage()`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md)),
`"ppc"` (posterior predictive histograms from
[`sbi_ppc()`](https://robustecologies.github.io/lucifer/reference/sbi_ppc.md)),
`"tarp"` (ECP curve from
[`TARP()`](https://robustecologies.github.io/lucifer/reference/TARP.md)),
and `"c2st"` (per-fold accuracy bars from
[`C2ST()`](https://robustecologies.github.io/lucifer/reference/C2ST.md)).
The `ground_truth` argument (or its alias `true_values`) adds vertical
dashed lines at known parameter values in posterior and pairs plots.

### Return value structure

The `sbi` object returned by
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
contains the following fields:

The `network` field can be passed to a subsequent
[`SBI()`](https://robustecologies.github.io/lucifer/reference/SBI.md)
call via `pretrained_network` for warm-start training. The
`standardization` field is used internally by the diagnostic functions
([`SBC()`](https://robustecologies.github.io/lucifer/reference/SBC.md),
[`expected_coverage()`](https://robustecologies.github.io/lucifer/reference/expected_coverage.md),
[`TARP()`](https://robustecologies.github.io/lucifer/reference/TARP.md))
and should not be modified. The `simulations` field provides access to
the raw training data, which can be useful for examining the
distribution of simulations across the parameter space or for computing
additional summary statistics post hoc. The `training` field’s `loss`
and `val_loss` vectors are plotted by `plot(fit, type = "training")` and
can be used to manually diagnose training quality or to compare training
curves across different hyperparameter settings.

## References

**\[1\]** Deistler, M., Gloeckler, M., Boelts, J., Karaletsos, T.,
Macke, J.H. et al. (2025). Simulation-based inference: a practical
guide. *arXiv preprint*.
[doi:10.48550/arXiv.2411.17337](https://doi.org/10.48550/arXiv.2411.17337)

**\[2\]** Papamakarios, G. & Murray, I. (2016). Fast epsilon-free
inference of simulation models with Bayesian conditional density
estimation. *Advances in Neural Information Processing Systems* 29.

**\[3\]** Lueckmann, J.M., Goncalves, P.J., Bassetto, G., Oecal, K.,
Nonnenmacher, M. & Macke, J.H. (2017). Flexible statistical inference
for mechanistic models of neural dynamics. *Advances in Neural
Information Processing Systems* 30.

**\[4\]** Cranmer, K., Brehmer, J. & Louppe, G. (2020). The frontier of
simulation-based inference. *Proceedings of the National Academy of
Sciences* 117(48), 30055-30062.
[doi:10.1073/pnas.1912789117](https://doi.org/10.1073/pnas.1912789117)

**\[5\]** Hermans, J., Begy, V. & Louppe, G. (2020). Likelihood-free
MCMC with amortized approximate ratio estimators. *Proceedings of the
37th International Conference on Machine Learning* (ICML).

**\[6\]** Talts, S., Betancourt, M., Simpson, D., Vehtari, A. & Gelman,
A. (2018). Validating Bayesian inference algorithms with
simulation-based calibration. *arXiv preprint*.
[doi:10.48550/arXiv.1804.06788](https://doi.org/10.48550/arXiv.1804.06788)

**\[7\]** Bishop, C.M. (1994). Mixture density networks. *Aston
University Technical Report* NCRG/94/004.

**\[8\]** Durkan, C., Bekasov, A., Murray, I. & Papamakarios, G. (2019).
Neural spline flows. *Advances in Neural Information Processing Systems*
32.

**\[9\]** Greenberg, D., Nonnenmacher, M. & Macke, J.H. (2019).
Automatic posterior transformation for likelihood-free inference.
*Proceedings of the 36th International Conference on Machine Learning*
(ICML).
[doi:10.48550/arXiv.1905.07488](https://doi.org/10.48550/arXiv.1905.07488)

**\[10\]** Deistler, M., Goncalves, P.J. & Macke, J.H. (2022). Truncated
proposals for scalable and hassle-free simulation-based inference.
*Advances in Neural Information Processing Systems* 35.
[doi:10.48550/arXiv.2210.04815](https://doi.org/10.48550/arXiv.2210.04815)

**\[11\]** Papamakarios, G., Sterratt, D. & Murray, I. (2019).
Sequential neural likelihood: fast likelihood-free inference with
autoregressive flows. *Proceedings of the 22nd International Conference
on Artificial Intelligence and Statistics* (AISTATS).
[doi:10.48550/arXiv.1805.07226](https://doi.org/10.48550/arXiv.1805.07226)

**\[12\]** Delaunoy, A., Hermans, J., Rozet, F., Wehenkel, A. & Louppe,
G. (2022). Towards reliable simulation-based inference with balanced
neural ratio estimation. *Advances in Neural Information Processing
Systems* 35.
[doi:10.48550/arXiv.2208.13624](https://doi.org/10.48550/arXiv.2208.13624)

**\[13\]** Lemos, P., Coogan, A., Hezaveh, Y. & Perreault-Levasseur, L.
(2023). Sampling-based accuracy testing of posterior estimators for
general inference. *Proceedings of the 40th International Conference on
Machine Learning* (ICML).
[doi:10.48550/arXiv.2302.03026](https://doi.org/10.48550/arXiv.2302.03026)

**\[14\]** Lopez-Paz, D. & Oquab, M. (2017). Revisiting classifier
two-sample tests. *International Conference on Learning Representations*
(ICLR).

**\[15\]** Wood, S.N. (2010). Statistical inference for noisy nonlinear
ecological dynamic systems. *Nature* 466(7310), 1102-1104.
[doi:10.1038/nature09319](https://doi.org/10.1038/nature09319)

**\[16\]** Beaumont, M.A., Zhang, W. & Balding, D.J. (2002). Approximate
Bayesian computation in population genetics. *Genetics* 162(4),
2025-2035.
[doi:10.1093/genetics/162.4.2025](https://doi.org/10.1093/genetics/162.4.2025)

**\[17\]** Blum, M.G.B., Nunes, M.A., Prangle, D. & Sisson, S.A. (2013).
A comparative review of dimension reduction methods in approximate
Bayesian computation. *Statistical Science* 28(2), 189-208.
[doi:10.1214/12-STS406](https://doi.org/10.1214/12-STS406)

**\[18\]** Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari,
A. & Rubin, D.B. (2013). *Bayesian Data Analysis*. 3rd edition. Chapman
& Hall/CRC. ISBN 978-1439840955.

**\[19\]** Kingma, D.P. & Ba, J. (2015). Adam: a method for stochastic
optimization. *Proceedings of the 3rd International Conference on
Learning Representations* (ICLR).
[doi:10.48550/arXiv.1412.6980](https://doi.org/10.48550/arXiv.1412.6980)

**\[20\]** Klambauer, G., Unterthiner, T., Mayr, A. & Hochreiter, S.
(2017). Self-normalizing neural networks. *Advances in Neural
Information Processing Systems* 30.
[doi:10.48550/arXiv.1706.02515](https://doi.org/10.48550/arXiv.1706.02515)

**\[21\]** Papamakarios, G., Nalisnick, E., Rezende, D.J., Mohamed, S. &
Balaram, L. (2021). Normalizing flows for probabilistic modeling and
inference. *Journal of Machine Learning Research* 22(57), 1-64.

**\[22\]** Radev, S.T., Mertens, U.K., Voss, A., Ardizzone, L. & Kothe,
U. (2020). BayesFlow: learning complex stochastic models with invertible
neural networks. *IEEE Transactions on Neural Networks and Learning
Systems* 33(4), 1452-1466.
[doi:10.1109/TNNLS.2020.3042395](https://doi.org/10.1109/TNNLS.2020.3042395)

**\[23\]** Tejero-Cantero, A., Boelts, J., Deistler, M., Lueckmann,
J.M., Durkan, C., Goncalves, P.J., Greenberg, D.S. & Macke, J.H. (2020).
sbi: a toolkit for simulation-based inference. *Journal of Open Source
Software* 5(52), 2505.
[doi:10.21105/joss.02505](https://doi.org/10.21105/joss.02505)

**\[24\]** Goncalves, P.J., Lueckmann, J.M., Deistler, M., Nonnenmacher,
M., Oecal, K., Bassetto, G., Chintaluri, C., Podlaski, W.F., Haddad,
S.A., Vogels, T.P., Greenberg, D.S. & Macke, J.H. (2020). Training deep
neural density estimators to identify mechanistic models of neural
dynamics. *eLife* 9, e56261.
[doi:10.7554/eLife.56261](https://doi.org/10.7554/eLife.56261)

**\[25\]** Lueckmann, J.M., Boelts, J., Greenberg, D.S., Goncalves, P.J.
& Macke, J.H. (2021). Benchmarking simulation-based inference.
*Proceedings of the 24th International Conference on Artificial
Intelligence and Statistics* (AISTATS).
[doi:10.48550/arXiv.2101.04653](https://doi.org/10.48550/arXiv.2101.04653)
