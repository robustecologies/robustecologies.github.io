# Merge basins for Wada analysis

Creates a two-color basin map by merging all basins except one into a
single "merged" basin. This is essential for the merging method of Wada
detection.

## Usage

``` r
merge_basins(basins, keep_basin)
```

## Arguments

- basins:

  Matrix of basin assignments or `wada_basins` object.

- keep_basin:

  Integer. Basin index to keep separate (will be color 1).

## Value

An object of class `merged_basins` containing:

- merged:

  Integer matrix with values 1 (kept basin) and 2 (merged).

- keep_basin:

  Integer. The basin that was kept separate.

- merged_basins:

  Integer vector. Indices of basins that were merged.

- n_attractors:

  Integer. Original number of attractors.

- x_grid:

  Numeric vector. X coordinates (if available from input).

- y_grid:

  Numeric vector. Y coordinates (if available from input).

- x_range:

  Numeric vector. X range (if available from input).

- y_range:

  Numeric vector. Y range (if available from input).

## Details

The merging method for detecting Wada basins works by comparing the
boundaries between a single basin and the union of all other basins. If
the basins are truly Wada, all these boundaries should be identical.

**The merging principle:**

A key insight from Daza et al. (2018) is that a Wada boundary is the
only one that remains unaltered under the action of merging basins. When
we merge \\N_A - 1\\ basins into one, the resulting *slim boundary*
between the kept basin and the merged basin should be identical
regardless of which basin is kept separate.

The merge operation assigns:

- Value 1 to points belonging to the kept basin \\B_k\\

- Value 2 to points belonging to the merged basin \\M_k = \bigcup\_{j
  \neq k} B_j\\

**Slim boundary:**

The boundary between the kept basin and the merged basin is called the
"slim boundary". For Wada basins, all slim boundaries (for different
values of `keep_basin`) are identical. This is verified by computing the
Hausdorff distance between slim boundaries in
[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md).

## References

Daza, A., Wagemakers, A., & Sanjuan, M. A. F. (2018). Ascertaining when
a basin is Wada: The merging method. *Scientific Reports*, 8, 9954.
[doi:10.1038/s41598-018-28119-0](https://doi.org/10.1038/s41598-018-28119-0)

Kennedy, J., & Yorke, J. A. (1991). Basins of Wada. *Physica D:
Nonlinear Phenomena*, 51(1-3), 213-225.
[doi:10.1016/0167-2789(91)90234-Z](https://doi.org/10.1016/0167-2789%2891%2990234-Z)

## See also

[`wada_merging_method`](https://robustecologies.github.io/wadaR/reference/wada_merging_method.md)
for the full Wada detection algorithm,
[`compute_basins`](https://robustecologies.github.io/wadaR/reference/compute_basins.md)
for computing basins,
[`get_boundary`](https://robustecologies.github.io/wadaR/reference/get_boundary.md)
for boundary extraction,
[`hausdorff_distance`](https://robustecologies.github.io/wadaR/reference/hausdorff_distance.md)
for the metric that compares boundaries,
[`print.merged_basins`](https://robustecologies.github.io/wadaR/reference/print.merged_basins.md),
[`summary.merged_basins`](https://robustecologies.github.io/wadaR/reference/summary.merged_basins.md),
[`plot.merged_basins`](https://robustecologies.github.io/wadaR/reference/plot.merged_basins.md)
for the S3 method trio.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example: Visualize merged basins for pendulum                          #
# ===================================================================== #
pendulum <- forced_damped_pendulum(forcing = 1.66)
basins <- compute_basins(pendulum, c(-pi, pi), c(-3, 3), resolution = 200)

# Merge keeping basin 1 separate: B1 vs (B2 + B3)
merged1 <- merge_basins(basins, keep_basin = 1)
plot(merged1)

# Merge keeping basin 2 separate: B2 vs (B1 + B3)
merged2 <- merge_basins(basins, keep_basin = 2)
plot(merged2, colors = c("#2A9D8F", "#E9C46A"))

# Compare slim boundaries (should be identical for Wada basins)
plot(merged1, show_boundary = TRUE)
} # }
```
