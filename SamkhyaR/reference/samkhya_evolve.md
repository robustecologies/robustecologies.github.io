# Run the complete evolution of Prakrti

Carries the evolution from the unmanifest equilibrium through to the
gross elements, following the branching derivation of Samkhya-Karika 22
and 25 rather than a single linear chain.

## Usage

``` r
samkhya_evolve(
  guna_dominant = c("sattva", "rajas", "tamas"),
  delta = 0.1,
  verbose = TRUE
)
```

## Arguments

- guna_dominant:

  Which quality predominates when the equilibrium breaks. Passed to
  [`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md).

- delta:

  Magnitude of the disturbance. Passed to
  [`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md).

- verbose:

  Emit progress messages at each step. Default `TRUE`.

## Value

A `samkhya_state` in which all twenty-three evolutes of Prakrti have
manifested.

## Details

The order in which the two branches from ahamkara are taken is
immaterial, since they issue in parallel; this function takes the
vaikrta branch first only because the karika enumerates the elevenfold
set first. Taking the bhutadi branch first yields an identical final
state, which the package test suite asserts.

## References

Isvarakrsna. *Samkhyakarika*, verses 22, 25.

## See also

[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md)
and the individual evolution functions for the step-by-step form,
[`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md)
for what follows completion,
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md),
[`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md).

## Examples

``` r
state <- samkhya_evolve(verbose = FALSE)
state
#> <samkhya_state>
#> Stage: Mahabhutas (the specific)
#> Gunas: sattva = 0.433, rajas = 0.283, tamas = 0.283 (sum = 1.000)
#> Evolutes manifest: 23 of 23
#> Cognitive situation: Visesa (the specific): the gross world stands manifest
#> Viveka in buddhi: Not attained
#> Prakrti acting for this Purusha: Yes
#> Purusha: Witness, neutral, non-agent; neither bound nor released (SK 19, 62)

# All twenty-three evolutes of Prakrti are present.
summary(state)$n_manifest
#> [1] 23
```
