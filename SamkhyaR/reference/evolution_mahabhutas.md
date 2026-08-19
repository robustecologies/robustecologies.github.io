# Evolve the five gross elements

Manifests the five *mahabhutas*, space, air, fire, water and earth,
which arise from the tanmatras and complete the material evolution
(Samkhya-Karika 22 and 38).

## Usage

``` r
evolution_mahabhutas(state, verbose = TRUE)
```

## Arguments

- state:

  A `samkhya_state` in which the tanmatras have manifested.

- verbose:

  Emit a progress message. Default `TRUE`.

## Value

The `samkhya_state` with the five gross elements manifested.

## Details

Samkhya-Karika 38 calls the gross elements the specific, *visesa*, and
says of them that they are tranquil, turbulent and deluding, *santa
ghoras ca mudhas ca*.

Two claims commonly attached to this step do not come from the karika.
The first is that each gross element accumulates the qualities of those
before it, so that space carries sound alone and earth carries all five.
SK 22 states only the order of derivation. Within Samkhya the
accumulation scheme is Vijnanabhiksu's position; Vacaspati Misra holds
instead that each gross element arises from its own tanmatra alone,
which is the derivation recorded in the `parent` column of
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md).
Outside Samkhya the scheme is Vedantic and Puranic, and belongs to the
*pancikarana* framework.

The second is that the gross elements constitute the subtle body. They
do not. SK 40 makes the subtle body the eighteen consisting of mahat,
ahamkara, the eleven organs and the five tanmatras; the gross elements
constitute the gross body alone, which perishes, while the subtle body
transmigrates. See
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md).

## References

Isvarakrsna. *Samkhyakarika*, verses 22, 38-40.

Vacaspati Misra. *Tattvakaumudi*, commentary on verse 22.

## See also

[`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md)
for their material cause,
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md)
for what the subtle body actually comprises,
[`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md),
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md).

## Examples

``` r
state <- samkhya_evolve(verbose = FALSE)
summary(state)$n_manifest
#> [1] 23

# The gross elements are not constituents of the subtle body.
intersect(subset(tattvas, group == "mahabhuta")$iast, linga_sharira()$constituents)
#> character(0)
```
