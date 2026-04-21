# Check if solid angle computation is feasible

Diagnoses whether a cone's solid angle can be computed and suggests the
best method.

## Usage

``` r
diagnose_cone(V)
```

## Arguments

- V:

  Matrix of cone generators

## Value

A list containing diagnostic information

## Examples

``` r
if (FALSE) { # \dontrun{
# Diagnose and analyze a 4D orthogonal cone
V <- diag(4)

# Basic diagnosis
d <- diagnose_cone(V)
print(d)

# Extended summary with spectral analysis
s <- summary(d)
print(s)

# Visualize eigenvalue structure
plot(d)
} # }
```
