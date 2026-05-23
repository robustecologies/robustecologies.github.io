# Matrix exponential via Pade approximation with scaling-and-squaring

Higham 2005 algorithm; degree-13 rational Pade approximant for the
matrix exponential `exp(t * A)`.

## Usage

``` r
expm_pade(A, t = 1)
```

## Arguments

- A:

  A square numeric matrix.

- t:

  Scalar time (default 1).

## Value

The matrix `exp(t * A)`.

## References

Higham, N. J. (2005). The scaling and squaring method for the matrix
exponential revisited. *SIAM J. Matrix Anal. Appl.* 26(4), 1179-1193.
[doi:10.1137/04061101X](https://doi.org/10.1137/04061101X) .

## See also

[`expm_krylov`](https://robustecologies.github.io/janos/reference/expm_krylov.md).

## Examples

``` r
if (FALSE) { # \dontrun{
A <- matrix(c(-2, 1, 0, -2), 2, 2)
expm_pade(A, t = 1)
} # }
```
