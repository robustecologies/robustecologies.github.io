# Create a network specification for SBI

Defines the neural network architecture used by
[`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md) for
density estimation or classification. The default Mixture Density
Network (MDN) backend requires no additional packages and handles most
problems with up to ~20 parameters. For larger problems or normalizing
flow architectures, install the torch package.

## Usage

``` r
sbi_network(
  type = "mdn",
  hidden_layers = c(128L, 128L),
  n_components = 10L,
  activation = "selu"
)
```

## Arguments

- type:

  Character: `"mdn"` (default) for Mixture Density Network, `"nsf"` for
  Neural Spline Flow (requires torch).

- hidden_layers:

  Integer vector of hidden layer sizes. Default `c(128, 128)`.

- n_components:

  Integer: number of mixture components for MDN. Default 10.

- activation:

  Character: activation function. Currently only `"selu"` is supported.

## Value

A list of class `sbi_network_spec` with the architecture specification.

## Details

Implementation of `sbi_network`. Refer to the package vignettes and the
cited references for a complete algorithmic and mathematical
description.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Default MDN
net <- sbi_network()

# Larger MDN for complex posteriors
net <- sbi_network(hidden_layers = c(256, 256, 128), n_components = 20)

# Neural Spline Flow (requires torch)
net <- sbi_network(type = "nsf", hidden_layers = c(128, 128))
} # }
```
