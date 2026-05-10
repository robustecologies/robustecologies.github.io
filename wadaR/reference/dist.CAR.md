# Intrinsic conditional autoregressive prior

Log-density and random generation for the intrinsic conditional
autoregressive (CAR) prior used in spatial statistics.

See Details.

See Details.

## Usage

``` r
dcar(x, W, tau, zero.mean = rep(0, length(x)), log = FALSE)

rcar(W, tau)
```

## Arguments

- x:

  vector of spatial random effects of length \\n\\.

- W:

  an \\n \times n\\ binary adjacency matrix (symmetric, with zeros on
  the diagonal).

- tau:

  precision (inverse variance) parameter, which must be positive.

- zero.mean:

  vector of length \\n\\ specifying the centering constraint (typically
  a vector of zeros).

- log:

  logical; if `TRUE`, the logarithm of the density is returned.

## Value

`dcar` gives the (log) density as a scalar. `rcar` generates a single
draw from the CAR prior.

See Details.

See Details.

## Details

The intrinsic CAR prior (Besag, York, and Mollie, 1991) is a pairwise
difference prior for spatially structured random effects on an irregular
lattice. The precision matrix is \\Q = \tau(D - W)\\ where \\D =
\mathrm{diag}(W \mathbf{1})\\ is the degree matrix and \\W\\ is the
adjacency matrix. This prior is improper (rank-deficient by 1 for a
connected graph) and the density is computed using the
pseudo-determinant (product of nonzero eigenvalues).

For identifiability, a sum-to-zero constraint on \\x\\ is conventionally
imposed in the model.

Implementation of `dcar`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

Implementation of `rcar`. Refer to the package vignettes and the cited
references for a complete algorithmic and mathematical description.

## References

Besag, J., York, J., and Mollie, A. (1991). "Bayesian Image Restoration,
with Two Applications in Spatial Statistics (with Discussion)". *Annals
of the Institute of Statistical Mathematics*, 43(1), p. 1–59.

## See also

[`dmvn`](https://robustecologies.github.io/lucifer/reference/dist.Multivariate.Normal.md)

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# 4-node lattice
W <- matrix(c(0,1,1,0, 1,0,0,1, 1,0,0,1, 0,1,1,0), 4, 4)
x <- c(0.5, -0.3, 0.1, -0.3)
d <- dcar(x, W, tau = 1, log = TRUE)
r <- rcar(W, tau = 1)
} # }

if (FALSE) { # \dontrun{
## see package vignettes for full examples involving dcar
} # }
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving rcar
} # }
```
