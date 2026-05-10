# Levy alpha-stable noise specification

Creates a noise specification for driving an SDE with alpha-stable Levy
increments instead of Gaussian Wiener increments. The
Chambers-Mallows-Stuck (CMS) algorithm generates symmetric or asymmetric
stable variates, with scaling dt^(1/alpha) replacing the Gaussian
sqrt(dt).

## Usage

``` r
levy_noise(alpha = 1.5, beta = 0)
```

## Arguments

- alpha:

  Stability index in (0, 2\]. alpha=2 recovers Gaussian noise. Smaller
  values produce heavier tails.

- beta:

  Skewness parameter in \[-1, 1\]. beta=0 gives symmetric stable noise.

## Value

An S3 object of class `noise_spec` with `type = "levy"`.

## Details

The alpha-stable distribution S(alpha, beta, 1, 0) generalizes the
Gaussian (alpha=2) and Cauchy (alpha=1) distributions. For alpha \< 2,
the distribution has heavy tails (infinite variance for alpha \< 2,
infinite mean for alpha \<= 1). The stability index alpha controls tail
heaviness while the skewness parameter beta controls asymmetry.

The SDE dX = f(X)dt + g(X)dL with alpha-stable Levy noise dL is
discretized as X(n+1) = X(n) + f(X(n))\*dt + g(X(n))*dL*dt^(1/alpha),
where dL is drawn from S(alpha, beta, 1, 0) via the CMS algorithm.

Milstein is not available for Levy-driven SDEs because the Ito isometry
(which underlies the Milstein correction) does not hold for non-Gaussian
noise. Only Euler-Maruyama is valid.

Levy noise can be composed with correlated noise (Cholesky-rotated
alpha-stable increments). It cannot be composed with fBm or colored
noise.

## References

Chambers, J. M., Mallows, C. L., & Stuck, B. W. (1976). A method for
simulating stable random variables. *Journal of the American Statistical
Association*, 71(354), 340-344.
[doi:10.1080/01621459.1976.10480344](https://doi.org/10.1080/01621459.1976.10480344)

Weron, R. (1996). On the Chambers-Mallows-Stuck method for simulating
skewed stable random variables. *Statistics & Probability Letters*,
28(2), 165-171.
[doi:10.1016/0167-7152(95)00113-1](https://doi.org/10.1016/0167-7152%2895%2900113-1)

## See also

[`correlated_noise`](https://robustecologies.github.io/janos/reference/correlated_noise.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- model_spec(
    rhs = list(x ~ -theta * x),
    diffusion = list(x ~ sigma),
    noise = levy_noise(alpha = 1.5, beta = 0),
    state_names = "x",
    parms = list(theta = 1, sigma = 0.3),
    init = c(x = 0)
)
result <- dyn_sim(model, t_max = 100,
                  solver = solver_euler_maruyama(dt = 0.01))
print(result)
plot(result)
} # }
```
