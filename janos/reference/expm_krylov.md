# Arnoldi-Krylov approximation of exp(t A) v

Approximates the action of the matrix exponential on a vector via
projection onto a Krylov subspace built by Arnoldi iteration. Returns
`exp(t A) v` without ever forming `exp(t A)`.

## Usage

``` r
expm_krylov(A, v, t = 1, m = 30, tol = 1e-10)
```

## Arguments

- A:

  A square numeric matrix.

- v:

  A numeric vector with the same length as the columns of A.

- t:

  Scalar time (default 1).

- m:

  Krylov subspace dimension (default 30).

- tol:

  Tolerance for early termination of the Arnoldi process (default
  1e-10).

## Value

Numeric vector `exp(t A) v`.

## References

Saad, Y. (1992). Analysis of some Krylov subspace approximations to the
matrix exponential operator. *SIAM J. Numer. Anal.* 29(1), 209-228.
[doi:10.1137/0729014](https://doi.org/10.1137/0729014) .

## See also

[`expm_pade`](https://robustecologies.github.io/janos/reference/expm_pade.md),
[`slow_fast_partition`](https://robustecologies.github.io/janos/reference/slow_fast_partition.md).

## Examples

``` r
if (FALSE) { # \dontrun{
A <- matrix(c(-2, 1, 0, -2), 2, 2)
expm_krylov(A, c(1, 0), t = 1)
} # }
```
