# Plot method for cone diagnosis

Visualizes the eigenvalue structure of the associated matrix for a
polyhedral cone. The bar plot shows eigenvalues colored by sign,
providing insight into the cone's geometric properties.

## Usage

``` r
# S3 method for class 'cone_diagnosis'
plot(x, type = "eigenvalues", ...)
```

## Arguments

- x:

  A cone_diagnosis object from
  [`diagnose_cone`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md)

- type:

  Type of plot. Currently only "eigenvalues" is supported.

- ...:

  Additional arguments passed to `barplot`

## Details

The eigenvalues of the associated matrix determine which computational
methods are applicable. Positive eigenvalues (shown in blue) indicate
that the corresponding directions contribute to a positive definite
structure. Negative eigenvalues (shown in red) indicate that
decomposition methods are required.

## Examples

``` r
if (FALSE) { # \dontrun{
# Positive definite cone (orthant)
V <- diag(4)
d <- diagnose_cone(V)
plot(d)

# Non-positive-definite cone
V <- matrix(c(1, 0.9, 0.1,
              0, 0.5, 0.8,
              0, 0, 0.3), nrow = 3)
d <- diagnose_cone(V)
plot(d)
} # }
```
