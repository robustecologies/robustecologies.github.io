# Summary method for a lyapunov_function object

Produces a detailed numerical diagnostic of a `lyapunov_function`
object: the value of \\V\\ at the equilibrium (which should vanish up to
round-off), a spot check over `n_check` random points in a small ball
around \\x^\*\\ reporting the fraction where \\V \> 0\\ and \\\dot V \<
0\\, and method-specific diagnostics (the condition number of \\P\\ for
quadratic; the Goh weights for the logarithmic Lyapunov function).

## Usage

``` r
# S3 method for class 'lyapunov_function'
summary(object, n_check = 20L, ...)
```

## Arguments

- object:

  A `lyapunov_function` object.

- n_check:

  Integer, number of random points for the spot check. Defaults to 20.

- ...:

  Unused, kept for S3 compatibility.

## Value

A list of class `summary.lyapunov_function` with the spot-check results,
invisibly.

## Details

The spot check evaluates \\V(x)\\ and \\\dot V(x)\\ at `n_check` random
points drawn from a small ball of radius \\0.1 \max(\|x^\*\|, 1)\\
around the equilibrium and reports the fraction satisfying \\V \> 0\\
and \\\dot V \< 0\\. For the Goh method, points are reflected into the
positive orthant so the logarithmic term stays defined.

## References

Goh, B. S. (1977). Global stability in many-species systems. *The
American Naturalist*, 111(977), 135-143.
[doi:10.1086/283144](https://doi.org/10.1086/283144)

## See also

[`lyapunov_function()`](https://robustecologies.github.io/janos/reference/lyapunov_function.md)
— constructor;
[`print.lyapunov_function()`](https://robustecologies.github.io/janos/reference/print.lyapunov_function.md)
— compact header;
[`plot.lyapunov_function()`](https://robustecologies.github.io/janos/reference/plot.lyapunov_function.md)
— level sets and surfaces;
[`analyse_lyapunov()`](https://robustecologies.github.io/janos/reference/analyse_lyapunov.md)
— higher-level family-aware dispatcher.

## Examples

``` r
if (FALSE) { # \dontrun{
A <- diag(c(-1, -2))
summary(lyapunov_function(A))
} # }
```
