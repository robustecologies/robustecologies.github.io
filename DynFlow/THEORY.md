# Numerical methods of DynFlow

This document states the mathematics behind every figure that DynFlow draws. It covers how a system written in the formula language is turned into a vector field or a map, which scheme advances it, how perturbations and shared noise enter, how equilibria and their stability are computed, and which tests check each of these against a closed form or a published value. The notation follows the source files, so that a statement here can be traced to the function that implements it; the tests cited in the last section can be run with `node tests/<name>.test.mjs`, and every number quoted there was produced by them.

<br>

## Systems and notation

A system is a list of statements of the form `x' = f(...)` for a flow or `x[n+1] = F(...)` for a map, together with optional `noise`, `aux`, `param`, `init` and `range` lines. Let $x \in \mathbb{R}^d$ be the state, $p \in \mathbb{R}^m$ the parameters and $t$ the time (the iteration index for a map). A flow is then

$$\dot x = f(t, x, p), \qquad t \ge t_0,$$

where $f : \mathbb{R} \times \mathbb{R}^d \times \mathbb{R}^m \to \mathbb{R}^d$ is assembled from the right-hand sides in the order in which the variables are declared, and a map is $x_{n+1} = F(n, x_n, p)$. A line `noise x = g` adds an Itô diffusion coefficient $g_i(t, x, p)$ for the variable $x_i$, so that the system becomes the stochastic differential equation

$$\mathrm{d}x_i = f_i(t, x, p)\,\mathrm{d}t + g_i(t, x, p)\,\mathrm{d}W_i,$$

where the $W_i$ are independent standard Wiener processes; the term `lag(x, tau)` inside an expression denotes $x(t - \tau)$ and turns the system into a delay differential equation. The compiler accepts only numbers, declared names, a closed list of functions and the arithmetic, comparison and logical operators, and it emits a JavaScript function that writes into preallocated arrays; a compiled formula therefore cannot call anything outside that list, so a scene received by link or file can be run without inspecting it.

<br>

## Integration schemes

Flows are advanced by the classical fourth-order Runge-Kutta scheme at a fixed step $h$. A fixed step is chosen deliberately, because an animation advances a whole ensemble in lockstep and the frames must correspond to equal intervals of model time. The price is that stiff systems, such as the Oregonator or the Robertson kinetics listed in janos, are left out of the catalogue: an explicit scheme would likely need a step several orders of magnitude smaller than the time scale of interest, which may make an animation of them impractical.

Stochastic systems are advanced by the Euler-Maruyama scheme,

$$x_{k+1} = x_k + f(t_k, x_k, p)\,h + g(t_k, x_k, p) \odot \Delta W_k, \qquad \Delta W_k \sim \mathcal{N}(0, h I),$$

which converges with strong order 1/2 and weak order 1 for globally Lipschitz coefficients (Kloeden and Platen 1992); the diffusion is interpreted in the sense of Itô throughout, so a Stratonovich model must be written with its drift correction added by hand.

Delay equations are advanced by the same Runge-Kutta scheme, with each delayed value read from a cubic Hermite interpolant of the stored solution and its derivative at the step points (Bellen and Zennaro 2003). Let $s \in [t_a, t_b]$ with $u = (s - t_a)/(t_b - t_a)$; then

$$x(s) \approx (2u^3 - 3u^2 + 1)\,x_a + (u^3 - 2u^2 + u)\,\delta\,\dot x_a + (-2u^3 + 3u^2)\,x_b + (u^3 - u^2)\,\delta\,\dot x_b,$$

where $\delta = t_b - t_a$ and $\dot x_a$, $\dot x_b$ are the right-hand sides evaluated at the stored points. The initial history is constant, $x(s) = x_0$ for $s \le t_0$, and the stored history begins with the point $(t_0, x_0)$ and the right derivative $f(t_0^+, x_0)$; without that first point the interpolant on $(t_0, t_0 + h)$ would be wrong, and the tests below detected this error in an earlier draft of the code. Note that the interpolant is exact for cubic solutions, so on problems whose solution is piecewise polynomial of degree three or less the scheme reproduces the solution to rounding error.

Maps are iterated directly. The doubling map $x \mapsto 2x \bmod 1$ raises a numerical caveat for any figure built on it: in binary floating point each iteration shifts one bit out of the mantissa, so every orbit reaches $x = 0$ within 53 steps. The Kaplan-Yorke entry of the catalogue therefore adds a Gaussian term of standard deviation $10^{-9}$, which re-injects the low-order bits; a perturbation of that size is unlikely to change the attractor at the resolution of any figure.

<br>

## Ensembles and common noise

Every figure is computed from an ensemble of $n$ copies of the system that share the parameters, the time and any forcing, so a particle flow, a density and a stroboscopic cloud are all drawn from one simulator. Noise enters in two ways. A perturbation marked as independent draws a separate increment for every member, which yields an ensemble of independent realisations and, in the limit of many members, the solution of the Fokker-Planck equation. A perturbation marked as common draws one increment per step and applies it to every member, as an environmental forcing does; the ensemble then samples a single realisation of a random dynamical system, and after a long enough transient its cloud is likely to approximate the snapshot (or pullback) attractor of that realisation. The same distinction is available inside formulas: `urand()` and `nrand()` are fresh for every member, whereas `ucommon(k)` and `ncommon(k)` are shared, and the random maps of the catalogue (the Zaslavsky map with a random phase and the random baker's map) are written with the shared draws. The generator is xoshiro128** (Blackman and Vigna 2021) seeded through splitmix32, normal variates come from the Box-Muller transform, and symmetric alpha-stable variates from the method of Chambers, Mallows and Stuck (1976); every scene records its seed, so an exported figure can be reproduced exactly.

<br>

## Perturbations

Perturbations act either on a parameter or on a state. Let $b$ be the base value of a parameter; the modulated value is

$$p(t) = b + A \sin(2\pi t/T + \varphi)$$

for periodic forcing, the sum of two such terms with incommensurate periods for quasiperiodic forcing, $b + r\,(t - t_0)$ clipped to a total change for a ramp, $b + A\,\Theta(t - t_0)$ for a step, and $b + \eta(t)$ for coloured parameter noise, where $\eta$ is an Ornstein-Uhlenbeck process with $\mathrm{d}\eta = -\eta/\tau\,\mathrm{d}t + \sigma\sqrt{2/\tau}\,\mathrm{d}W$ and hence stationary variance $\sigma^2$. State perturbations add $\sigma\,\mathrm{d}W$ (additive), $\sigma x\,\mathrm{d}W$ (multiplicative), $\eta\,\mathrm{d}t$ with $\eta$ an Ornstein-Uhlenbeck process (coloured), $\sigma h^{1/\alpha} S_\alpha$ with $S_\alpha$ a symmetric alpha-stable variate (Lévy), Poisson jumps of fixed or random size or of a fixed fraction of the state, and periodic pulses that remove a fraction of the state or add a fixed amount, as a harvest or a stocking does. The ramp is the device for rate-induced tipping. In the normal form $\dot x = (x + \lambda)^2 - 1$ with $\lambda$ ramped at rate $r$ (Ashwin et al. 2012), the variable $y = x + \lambda$ obeys $\dot y = y^2 - 1 + r$, which has an equilibrium only for $r \le 1$. A ramp faster than $r = 1$ that lasts longer than the escape time therefore tips the system, even though every frozen value of $\lambda$ admits a stable state.

<br>

## Local analysis

The Jacobian $J = \partial f/\partial x$ is computed by central differences with step $10^{-6}\max(1, |x_j|)$ in each coordinate. Its eigenvalues are obtained by reduction to upper Hessenberg form with stabilised elementary transformations followed by the shifted QR iteration of the routine hqr (Wilkinson and Reinsch 1971), which returns real and complex-conjugate pairs without forming the characteristic polynomial. Equilibria of a flow, or fixed points of a map, are the zeros of $G(x) = f(x)$ or $G(x) = F(x) - x$. They are sought by damped Newton iteration from the points of a Halton sequence over the axis box together with the solutions found at the previous parameter value, and solutions closer than $10^{-6}$ of the box size are merged. An equilibrium of a flow is classed as stable when every eigenvalue has negative real part, and a fixed point of a map when every eigenvalue lies inside the unit circle, with a margin of $10^{-7}$; saddles, foci and nodes are distinguished by the signs of the real parts and the presence of imaginary parts. The sweep view continues these branches across a parameter interval, and the phase view applies the same search to the plane spanned by two variables with the others held at their initial values.

The largest Lyapunov exponent quoted in the catalogue is estimated from two nearby orbits renormalised at fixed intervals (Benettin et al. 1980), after a transient; for a flow it is a rate per unit time and for a map a rate per iteration. Finite horizons tend to bias such estimates, which may explain why the Rössler value below sits slightly above the published one.

<br>

## Validation record

Four test files check the engine and the catalogue against results that are known independently of the code, and their output of 21 September 2026 is summarised here. The core tests confirm the moments of the random variates: the normal fourth moment is 2.978 against 3, the alpha-stable variance 1.990 against 2 at $\alpha = 2$, and the Cauchy median 0.997 against 1. They also confirm the order of the Runge-Kutta scheme, with an error ratio of 16.7 on halving the step against 16, and the exact stationary variance of the Euler-Maruyama discretisation of the Ornstein-Uhlenbeck process, $\sigma^2/(\theta(2 - \theta h))$, to within 1.4 percent over 20000 members. Finally, they check the delay scheme on $\dot x = -x(t - 1)$ with unit history, whose method-of-steps solution gives $x(2) = -1/2$ and $x(3) = -1/6$, reproduced to rounding error. The analysis tests check 280 Gaussian matrices of order 2 to 10, for which the eigenvalues reproduce the trace to $1.6 \times 10^{-14}$ and satisfy $|\det(A - \lambda I)|$ below $10^{-14}$ in relative terms. They also recover the Lorenz equilibria, the Hopf point $B = 1 + A^2$ of the Brusselator and the Hopf point $K = 3.75$ of the Rosenzweig-MacArthur model.

The model tests run all 82 catalogue entries and check the claims stated in their descriptions. The largest Lyapunov exponents are 0.9025 for the Lorenz system (reference 0.906), 0.4191 for the Hénon map (Hénon 1976; reference 0.419), 0.4701 for the Lozi map, 0.0739 for the Rössler system, $\ln 2$ to four decimals for the logistic map at $r = 4$, and 0.0205 for the four-species Lotka-Volterra community of Vano et al. (2006), who report 0.0203. Conserved quantities hold to within $10^{-8}$ in relative terms (the Lotka-Volterra first integral, the product $p_1 p_2 p_3$ of the rock-paper-scissors replicator, the pendulum energy), the residence times near the saddles of the May-Leonard system grow as expected for an attracting heteroclinic cycle (May and Leonard 1975), and the SIR model satisfies the final-size relation to $4 \times 10^{-12}$. For the five-species resource competition of Huisman and Weissing (2001), the total content of the first resource converges to the supply, as it must when mortality equals dilution. The browser tests load every model in the studio, run every applicable view, and decode the exported GIF with Pillow and the exported WebM with ffprobe. The checks are those of a figure-making tool and not of a solver library; in particular, no test certifies accuracy for stiff, chaotic or long integrations beyond the horizons stated above, where the fixed-step schemes may lose accuracy before any figure shows it.

<br>

## References

<a id="ref1"></a>**[1]** Ashwin, P., Wieczorek, S., Vitolo, R. and Cox, P. (2012). *Tipping points in open systems: bifurcation, noise-induced and rate-dependent examples in the climate system*. Philosophical Transactions of the Royal Society A 370: 1166-1184. [doi:10.1098/rsta.2011.0306](https://doi.org/10.1098/rsta.2011.0306)

<a id="ref2"></a>**[2]** Bellen, A. and Zennaro, M. (2003). *Numerical Methods for Delay Differential Equations*. Oxford University Press. [doi:10.1093/acprof:oso/9780198506546.001.0001](https://doi.org/10.1093/acprof:oso/9780198506546.001.0001)

<a id="ref3"></a>**[3]** Benettin, G., Galgani, L., Giorgilli, A. and Strelcyn, J.-M. (1980). *Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them*. Meccanica 15: 9-20. [doi:10.1007/BF02128236](https://doi.org/10.1007/BF02128236)

<a id="ref4"></a>**[4]** Blackman, D. and Vigna, S. (2021). *Scrambled linear pseudorandom number generators*. ACM Transactions on Mathematical Software 47: 36. [doi:10.1145/3460772](https://doi.org/10.1145/3460772)

<a id="ref5"></a>**[5]** Chambers, J. M., Mallows, C. L. and Stuck, B. W. (1976). *A method for simulating stable random variables*. Journal of the American Statistical Association 71: 340-344. [doi:10.1080/01621459.1976.10480344](https://doi.org/10.1080/01621459.1976.10480344)

<a id="ref6"></a>**[6]** Hénon, M. (1976). *A two-dimensional mapping with a strange attractor*. Communications in Mathematical Physics 50: 69-77. [doi:10.1007/BF01608556](https://doi.org/10.1007/BF01608556)

<a id="ref7"></a>**[7]** Huisman, J. and Weissing, F. J. (2001). *Fundamental unpredictability in multispecies competition*. The American Naturalist 157: 488-494. [doi:10.1086/319929](https://doi.org/10.1086/319929)

<a id="ref8"></a>**[8]** Kloeden, P. E. and Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations*. Springer. [doi:10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)

<a id="ref9"></a>**[9]** May, R. M. and Leonard, W. J. (1975). *Nonlinear aspects of competition between three species*. SIAM Journal on Applied Mathematics 29: 243-253. [doi:10.1137/0129022](https://doi.org/10.1137/0129022)

<a id="ref10"></a>**[10]** Vano, J. A., Wildenberg, J. C., Anderson, M. B., Noel, J. K. and Sprott, J. C. (2006). *Chaos in low-dimensional Lotka-Volterra models of competition*. Nonlinearity 19: 2391-2404. [doi:10.1088/0951-7715/19/10/006](https://doi.org/10.1088/0951-7715/19/10/006)

<a id="ref11"></a>**[11]** Wilkinson, J. H. and Reinsch, C. (1971). *Handbook for Automatic Computation, Volume II: Linear Algebra*. Springer. [doi:10.1007/978-3-642-86940-2](https://doi.org/10.1007/978-3-642-86940-2)
