# Extract basin at specific parameter value

Extracts a wada_basins object at a specific parameter value or index
from a bifurcation analysis result.

## Usage

``` r
extract_basin(x, param_value = NULL, index = NULL)
```

## Arguments

- x:

  A bifurcation_result object.

- param_value:

  Numeric. The parameter value (selects closest).

- index:

  Integer. Alternative: specify index directly.

## Value

A basin_result object that can be used with other wadaR functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract basin at gamma = 0.35
basin_035 <- extract_basin(bif_result, param_value = 0.35)

# Use with other functions
plot(basin_035)
entropy <- basin_entropy(basin_035)
} # }
```
