# Compute a single term in the tridiagonal series

Computes the term corresponding to multiexponent \\b = (b_1, \ldots,
b\_{n-1})\\ in equation (23).

## Usage

``` r
compute_tridiagonal_term(b, beta)
```

## Arguments

- b:

  Vector of exponents (length n-1)

- beta:

  Vector of consecutive dot products (length n-1)

## Value

The term value
