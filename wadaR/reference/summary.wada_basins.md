# Summary method for wada_basins objects

Provides extended statistics including basin areas, boundary analysis,
and area ratios for the computed basins of attraction.

## Usage

``` r
# S3 method for class 'wada_basins'
summary(object, ...)
```

## Arguments

- object:

  A wada_basins object.

- ...:

  Additional arguments (ignored).

## Value

A summary.wada_basins object (invisibly) containing:

- resolution:

  Grid resolution

- n_attractors:

  Number of attractors

- basin_counts:

  Number of points in each basin

- basin_fractions:

  Fraction of phase space occupied by each basin

- n_boundary:

  Number of boundary points

- boundary_fraction:

  Fraction of points on boundary

- unclassified:

  Number of unclassified points
