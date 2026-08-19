# Evolve manas, the eleventh organ

Manifests *manas*, which Samkhya-Karika 27 calls an organ partaking of
the nature of both classes, characterised by *samkalpa*, construction or
determination. Manas is the first member of the elevenfold set issuing
from the sattvic aspect of ahamkara.

## Usage

``` r
evolution_manas(state, verbose = TRUE)
```

## Arguments

- state:

  A `samkhya_state` in which ahamkara has already manifested.

- verbose:

  Emit a progress message. Default `TRUE`.

## Value

The `samkhya_state` with manas manifested.

## Details

Manas belongs to two classifications at once and it is worth keeping
them apart. Causally it is one of the eleven of SK 25, issuing from the
vaikrta aspect of ahamkara together with the ten organs. Functionally it
is one of the three constituents of the internal organ, the
*antahkarana* of SK 33, together with buddhi and ahamkara; the internal
organ is threefold, never twofold. The columns `ekadasaka` and
`antahkarana` of
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md)
record the two memberships separately for this reason.

SK 27 attributes to manas *samkalpaka*, the constructive or
determinative function. It does not attribute *vikalpa* or *savikalpa*,
terms that belong to the Yoga and Buddhist vocabularies respectively and
that are frequently imported into accounts of this verse.

## References

Isvarakrsna. *Samkhyakarika*, verses 25, 27, 33.

## See also

[`evolution_ahamkara()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_ahamkara.md)
for its material cause,
[`evolution_indriyas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_indriyas.md)
for the remaining ten of the elevenfold set,
[`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md)
for the parallel branch,
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md).

## Examples

``` r
state <- init_prakriti() |>
    evolution_buddhi(verbose = FALSE) |>
    evolution_ahamkara(verbose = FALSE) |>
    evolution_manas(verbose = FALSE)

# Threefold internal organ (SK 33), now complete.
summary(state)$antahkarana
#> [1] "mahat"    "ahaṃkāra" "manas"   
```
