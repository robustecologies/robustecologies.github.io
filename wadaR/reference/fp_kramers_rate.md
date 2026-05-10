# Kramers escape rate for noise-induced transitions

For a system with metastable states separated by potential barriers,
computes the Kramers escape rate, the mean first-passage time (MFPT),
and the spectral gap of the Fokker-Planck operator. Supports both 1D and
2D models. For 2D systems, the multi-dimensional Kramers formula using
Hessian eigenvalues at critical points is applied.

## Usage

``` r
fp_kramers_rate(
  model,
  xlim,
  ylim = NULL,
  n_grid = 200L,
  epsilon = 0.1,
  verbose = TRUE
)
```

## Arguments

- model:

  A `model_spec` with one or two state variables.

- xlim:

  Numeric vector of length 2 giving the domain for the first state
  variable.

- ylim:

  Numeric vector of length 2 giving the domain for the second state
  variable (required for 2D models, ignored for 1D).

- n_grid:

  Integer number of grid points (default 200). For 2D models, a scalar
  is used for both dimensions.

- epsilon:

  Positive numeric noise intensity (default 0.1).

- verbose:

  Logical; print progress messages (default TRUE).

## Value

An S3 object of class `fp_kramers` containing:

- transitions:

  Data frame of all well-to-well transitions with columns: from_x, to_x,
  barrier_x, barrier_height, rate, mfpt. For 2D, from/to are formatted
  as coordinate strings and saddle_x, saddle_y columns are included.

- spectral_gap:

  Spectral gap from the FP operator (1D only; NA for 2D).

- potential:

  The reconstructed `fp_potential` object (1D) or `NULL` (2D).

- fp:

  The `fp_stationary` or `fp_stationary_2d` object.

- epsilon:

  Noise intensity.

- model:

  The input `model_spec`.

## Details

Kramers (1940) showed that for a particle in a double-well potential
\\V(x)\\ driven by white noise of intensity \\\varepsilon\\, the escape
rate from a well at \\x_s\\ over a barrier at \\x_b\\ is \$\$r_K =
\frac{\sqrt{\|V''(x_s)\| \cdot \|V''(x_b)\|}}{2\pi}
\exp\\\left(-\frac{2\Delta V}{\varepsilon}\right)\$\$ where \\\Delta V =
V(x_b) - V(x_s)\\ is the barrier height. Since \\V''(x) = -f'(x)\\, the
formula involves the absolute values of the drift derivative at well and
barrier locations.

For 2D systems, the generalized Kramers formula (Langer, 1969) is
\$\$r_K =
\frac{\omega\_{\min}\\\|\omega\_{\mathrm{saddle}}\|}{2\pi\gamma}
\exp\\\left(-\frac{\Delta V}{\varepsilon}\right)\$\$ where
\\\omega\_{\min}\\ and \\\omega\_{\mathrm{saddle}}\\ are computed from
the Hessian eigenvalues of the quasi-potential at the minimum and saddle
point, respectively, and \\\gamma\\ is the damping factor (taken as the
geometric mean of the negative Hessian eigenvalues at the saddle). The
quasi-potential is extracted from the 2D stationary density via \\V(x,y)
= -\varepsilon \log p(x,y)\\.

The mean first-passage time is \\\tau = 1/r_K\\. The spectral gap
\\\|\lambda_2\|\\ of the FP operator converges to the sum of forward and
backward Kramers rates as \\\varepsilon \to 0\\, providing an
independent check: for a symmetric double well, \\\|\lambda_2\| \approx
2 r_K\\.

## References

Kramers, H. A. (1940). Brownian motion in a field of force and the
diffusion model of chemical reactions. *Physica*, 7(4), 284-304. DOI:
[doi:10.1016/S0031-8914(40)90098-2](https://doi.org/10.1016/S0031-8914%2840%2990098-2)
.

Langer, J. S. (1969). Statistical theory of the decay of metastable
states. *Annals of Physics*, 54(2), 258-275. DOI:
[doi:10.1016/0003-4916(69)90153-5](https://doi.org/10.1016/0003-4916%2869%2990153-5)
.

Hanggi, P., Talkner, P. & Borkovec, M. (1990). Reaction-rate theory:
fifty years after Kramers. *Reviews of Modern Physics*, 62(2), 251-341.
DOI:
[doi:10.1103/RevModPhys.62.251](https://doi.org/10.1103/RevModPhys.62.251)
.

## See also

[`fp_stationary`](https://robustecologies.github.io/janos/reference/fp_stationary.md),
[`fp_potential`](https://robustecologies.github.io/janos/reference/fp_potential.md),
[`fp_stationary_2d`](https://robustecologies.github.io/janos/reference/fp_stationary_2d.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1D double-well: cubic drift f(x) = x - x^3
dw <- model_spec(
    rhs = list(x ~ x - x^3),
    state_names = "x",
    parms = list(),
    init = c(x = 0.5)
)

kr <- fp_kramers_rate(dw, xlim = c(-2, 2), epsilon = 0.3)

print(kr)
summary(kr)
plot(kr)
plot(kr, type = "rate_vs_epsilon")

# 2D system with two stable equilibria
dw2d <- model_spec(
    rhs = list(x ~ x - x^3, y ~ -y),
    diffusion = list(x ~ 1, y ~ 1),
    state_names = c("x", "y"),
    parms = list(),
    init = c(x = 0.5, y = 0)
)

kr2 <- fp_kramers_rate(dw2d, xlim = c(-2, 2), ylim = c(-2, 2),
                         epsilon = 0.3)
print(kr2)
} # }
```
