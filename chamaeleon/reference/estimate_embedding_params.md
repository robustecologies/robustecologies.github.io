# Estimate embedding parameters

Convenience function to estimate both the embedding dimension \\m\\ and
the time delay \\\tau\\ for Takens phase-space reconstruction from a
single scalar time series. The two parameters are estimated
independently from the autocorrelation function (for \\\tau\\) and the
false-nearest- neighbour criterion (for \\m\\); the diagnostic curves
used to make the choice are returned alongside the parameter estimates.

## Usage

``` r
estimate_embedding_params(
  x,
  method_tau = "acf_half",
  max_dim = 10,
  verbose = TRUE
)
```

## Arguments

- x:

  Numeric vector. Univariate time series.

- method_tau:

  Character. Method for \\\tau\\ estimation: `"acf_half"` uses the first
  lag at which the autocorrelation drops below 0.5; `"acf_1e"` uses the
  first lag at which the autocorrelation drops below \\1/e\\.

- max_dim:

  Integer. Maximum embedding dimension tested by the
  false-nearest-neighbour criterion.

- verbose:

  Logical. If `TRUE` (default), print progress messages.

## Value

A list with components `tau` (integer; estimated delay), `m` (integer;
estimated embedding dimension), `acf` (numeric vector; autocorrelation
values used for the \\\tau\\ decision), and `fnn` (numeric vector; FNN
proportion at each tested dimension).

## Details

**Time delay \\\tau\\.** The autocorrelation function (ACF) of the
series is computed for lags up to `min(n/4, 500)`. The first lag at
which the ACF crosses the threshold (`0.5` or \\1/e\\, depending on
`method_tau`) is taken as \\\tau\\. If no crossing is detected, \\\tau\\
falls back to one tenth of the maximum lag.

**Embedding dimension \\m\\.** For each candidate dimension \\m \in \\1,
\ldots, \mathtt{max\\dim}\\\\, the proportion of false nearest
neighbours \\F(m)\\ is computed via
`fnn_proportion_cpp(x, m, tau, rtol = 15, atol = 2)`. The smallest \\m\\
satisfying \\F(m) \< 0.01\\ is selected; if the threshold is never
crossed, `max_dim` is returned.

**Caveats.** For very short series the FNN estimator is noisy; for
multivariate or non-stationary series the user should compute embedding
parameters per channel or per window respectively. The function does not
attempt mutual-information based estimation; for that, call
`takens_embed(x, tau = NULL, method_tau = "mutual_info")` directly.

## References

Fraser AM, Swinney HL (1986). Independent coordinates for strange
attractors from mutual information. Physical Review A 33:1134-1140.
[doi:10.1103/PhysRevA.33.1134](https://doi.org/10.1103/PhysRevA.33.1134)

Kennel MB, Brown R, Abarbanel HDI (1992). Determining embedding
dimension for phase-space reconstruction using a geometrical
construction. Physical Review A 45:3403-3411.
[doi:10.1103/PhysRevA.45.3403](https://doi.org/10.1103/PhysRevA.45.3403)

Takens F (1981). Detecting strange attractors in turbulence. In: Lecture
Notes in Mathematics, 898. Springer.
[doi:10.1007/BFb0091924](https://doi.org/10.1007/BFb0091924)

## See also

[`takens_embed`](https://robustecologies.github.io/chamaeleon/reference/takens_embed.md)
for the embedding itself (which can also auto-estimate \\\tau\\ and
\\m\\ internally);
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for the full chameleon pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- sin(seq(0, 20 * pi, length.out = 1000)) + 0.1 * rnorm(1000)
params <- estimate_embedding_params(x, verbose = FALSE)
print(params)

# Use the estimated parameters for the embedding
embedded <- takens_embed(x, m = params$m, tau = params$tau)
} # }
```
