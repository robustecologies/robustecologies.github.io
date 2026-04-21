# Fractional Brownian motion noise specification

Creates a noise specification for driving an SDE with fractional
Brownian motion (fBm) increments. The Hurst parameter H in (0,1)
controls the correlation structure: H \> 0.5 produces persistent
(positively correlated) increments, H \< 0.5 produces anti-persistent
increments, and H = 0.5 recovers standard Brownian motion.

## Usage

``` r
fbm_noise(H = 0.75, method = c("auto", "circulant", "hosking"))
```

## Arguments

- H:

  Hurst parameter in (0, 1). H=0.5 recovers standard Brownian motion.

- method:

  Increment generation method: "auto" (default), "circulant", or
  "hosking". Auto selects circulant when possible, falling back to
  Hosking.

## Value

An S3 object of class `noise_spec` with `type = "fbm"`.

## Details

Fractional Brownian motion B_H(t) is a centered Gaussian process with
covariance E\[B_H(s) B_H(t)\] = 0.5 \* (\|s\|^(2H) + \|t\|^(2H) -
\|t-s\|^(2H)). Unlike standard Brownian motion, fBm has stationary but
(for H != 0.5) correlated increments. The variance scales as Var(B_H(t))
~ t^(2H).

The fBm increments are pre-generated in R using the circulant embedding
method (Wood-Chan, O(n log n)) with the Hosking method as a fallback
when the circulant eigenvalues are negative. The pre-generated increment
vector is passed to the C++ integrator.

Because fBm is not a semimartingale for H != 0.5, the Euler-Maruyama
discretization is a pathwise approximation valid only for additive noise
(g(X) = constant). Milstein is not available. For multiplicative noise
with fBm, the result should be interpreted as a formal approximation.

fBm noise is standalone: it cannot be composed with Levy, correlated, or
colored noise.

## References

Mandelbrot, B. B. & Van Ness, J. W. (1968). Fractional Brownian motions,
fractional noises and applications. *SIAM Review*, 10(4), 422-437.
[doi:10.1137/1010093](https://doi.org/10.1137/1010093)

Wood, A. T. A. & Chan, G. (1994). Simulation of stationary Gaussian
processes in \[0,1\]^d. *Journal of Computational and Graphical
Statistics*, 3(4), 409-432.
[doi:10.2307/1390903](https://doi.org/10.2307/1390903)

Nualart, D. (2006). Stochastic calculus with respect to fractional
Brownian motion. *Annales de la Faculte des Sciences de Toulouse*,
15(1), 63-78.

## See also

[`colored_noise`](https://robustecologies.github.io/janos/reference/colored_noise.md),
[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- model_spec(
    rhs = list(x ~ -theta * x),
    diffusion = list(x ~ sigma),
    noise = fbm_noise(H = 0.75),
    state_names = "x",
    parms = list(theta = 1, sigma = 0.3),
    init = c(x = 0)
)
result <- dyn_sim(model, t_max = 50,
                  solver = solver_euler_maruyama(dt = 0.01))
print(result)
plot(result)
} # }
```
