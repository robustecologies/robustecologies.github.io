# Algorithm registry information

Returns a data frame summarizing all registered inference algorithms in
lucifer, with their metadata and classification. Useful for discovering
available methods and their properties.

## Usage

``` r
algo_info(category = NULL, subcategory = NULL, requires_gradient = NULL)
```

## Arguments

- category:

  Optional character. Filter by inference category (e.g., `"MCMC"`,
  `"VB"`, `"Laplace"`).

- subcategory:

  Optional character. Filter by subcategory (e.g., `"gradient"`,
  `"ensemble"`, `"multimodal"`).

- requires_gradient:

  Optional logical. If `TRUE`, return only gradient-requiring methods;
  if `FALSE`, only gradient-free methods.

## Value

A data frame with columns: Abbreviation, Full.Name, Category,
Subcategory, Gradient, Discrete, Componentwise, Quality.Tier,
Multimodal.Affinity, Constraint.Affinity, Torch, Eval.Cost.

## Details

Implementation of `algo_info`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# All algorithms
algo_info()

# Only gradient-based MCMC
algo_info(category = "MCMC", subcategory = "gradient")

# Gradient-free methods
algo_info(requires_gradient = FALSE)
} # }
```
