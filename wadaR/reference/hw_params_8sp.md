# Create parameters for 8-species competition scenario

Create parameters for 8-species competition scenario

## Usage

``` r
hw_params_8sp(
  dilution = 0.25,
  supply = c(10, 10, 10),
  t_invade = 1000,
  N_invade = 0.1
)
```

## Arguments

- dilution:

  Dilution rate D (default 0.25).

- supply:

  Supply concentrations (default c(10, 10, 10)).

- t_invade:

  Invasion time for species 6-8.

- N_invade:

  Initial abundance of invading species.

## Value

List of model parameters.
