# Create a discrete outcome specification

Defines a discrete outcome for systems with categorical final states,
such as multi-species competition models where outcomes are defined by
which species survive.

## Usage

``` r
attractor_outcome(id, label, survivors = NULL, color = NULL)
```

## Arguments

- id:

  Integer or character. Unique identifier for the outcome.

- label:

  Character. Human-readable description of the outcome.

- survivors:

  Character vector. Optional list of surviving species (for competition
  models).

- color:

  Character. Optional color for plotting this outcome.

## Value

A list with class `"attractor_spec"` containing:

- id:

  Outcome identifier

- label:

  Outcome label

- survivors:

  List of survivors (or NULL)

- color:

  Plot color (or NULL)

- type:

  "outcome"

## See also

[`attractor_point`](https://robustecologies.github.io/wadaR/reference/attractor_point.md),
[`compiled_system`](https://robustecologies.github.io/wadaR/reference/compiled_system.md)

## Examples

``` r
# Competition outcomes
outcome1 <- attractor_outcome(id = 1, label = "Species A wins",
                              survivors = "A", color = "#E41A1C")
outcome2 <- attractor_outcome(id = 2, label = "Species B wins",
                              survivors = "B", color = "#377EB8")
outcome3 <- attractor_outcome(id = 3, label = "Coexistence",
                              survivors = c("A", "B"), color = "#4DAF4A")
```
