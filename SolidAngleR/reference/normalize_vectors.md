# Normalize vectors to unit length

Converts a set of vectors to unit vectors by dividing each by its norm.

## Usage

``` r
normalize_vectors(V)
```

## Arguments

- V:

  A matrix where each column is a vector to normalize

## Value

A matrix of the same dimensions with unit-length columns

## Examples

``` r
# Normalize some vectors
V <- cbind(c(3, 4), c(1, 1))
V_unit <- normalize_vectors(V)

# Check norms
apply(V_unit, 2, function(v) sqrt(sum(v^2)))  # All = 1
#> [1] 1 1
```
