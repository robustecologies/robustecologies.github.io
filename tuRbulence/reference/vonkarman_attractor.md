# Create embedded attractor from peaks

Constructs the embedded attractor by creating delay vectors from the
extracted peak sequence.

## Usage

``` r
vonkarman_attractor(peaks, embed_dim = 3L)
```

## Arguments

- peaks:

  Object of class "vonkarman_peaks".

- embed_dim:

  Embedding dimension (default 3).

## Value

Object of class "vonkarman_attractor" with embedded coordinates.

## See also

[`vonkarman_sim()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_sim.md),
[`print.vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/print.vonkarman_attractor.md),
[`summary.vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/summary.vonkarman_attractor.md),
[`plot.vonkarman_attractor()`](https://robustecologies.github.io/tuRbulence/reference/plot.vonkarman_attractor.md),
[`vonkarman_peaks()`](https://robustecologies.github.io/tuRbulence/reference/vonkarman_peaks.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate and extract peaks
sim <- vonkarman_sim(mu = 0.3, n_steps = 100000, seed = 42)
peaks <- vonkarman_peaks(sim, min_separation = 0.1)

# Create embedded attractor
attr <- vonkarman_attractor(peaks, embed_dim = 3)
print(attr)
} # }
```
