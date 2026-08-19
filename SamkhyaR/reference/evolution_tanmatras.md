# Evolve the five tanmatras

Manifests the five *tanmatras*, the bare essences of sound, touch, form,
taste and smell. Samkhya-Karika 25 derives them from the tamasic aspect
of ahamkara, called *bhutadi*, the beginning of the elements.

## Usage

``` r
evolution_tanmatras(state, verbose = TRUE)
```

## Arguments

- state:

  A `samkhya_state` in which ahamkara has already manifested. The organs
  are not required: the two branches issue in parallel.

- verbose:

  Emit a progress message. Default `TRUE`.

## Value

The `samkhya_state` with the five tanmatras manifested.

## Details

The tanmatras descend from ahamkara, not from the organs. Diagrams that
route the organs of knowledge or of action into the tanmatras invert the
derivation of SK 25, in which the elevenfold set and the tanmatras are
two branches from a single cause. This function therefore requires only
ahamkara, and a state carrying ahamkara but no organ at all is a
legitimate argument to it.

Samkhya-Karika 38 reads: *tanmatrany avisesas tebhyo bhutani panca
pancabhyah / ete smrta visesah santa ghoras ca mudhas ca*. The tanmatras
are the non-specific, *avisesa*; from them arise the five gross
elements; and it is *these*, the specific or *visesa*, that are said to
be tranquil, turbulent and deluding. The frequent statement that SK 38
describes the tanmatras as neither soothing nor tormenting nor
stupefying reverses the predication and borrows the triple *priti*,
*apriti*, *visada* from SK 12, where it characterises the three
qualities. The negative characterisation of the tanmatras is a
commentarial inference, not the text of the verse.

The tanmatras belong to the seven of SK 3 that are at once effects and
causes, and to the eighteen constituents of the subtle body of SK 40.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 25, 38, 40.

## See also

[`evolution_ahamkara()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_ahamkara.md)
for their material cause,
[`evolution_mahabhutas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_mahabhutas.md)
for what descends from them,
[`evolution_indriyas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_indriyas.md)
for the parallel branch,
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md).

## Examples

``` r
# The tamasic branch may be taken with no organ manifested at all.
state <- init_prakriti() |>
    evolution_buddhi(verbose = FALSE) |>
    evolution_ahamkara(verbose = FALSE) |>
    evolution_tanmatras(verbose = FALSE)

summary(state)$branches
#> vaikrta bhutadi 
#>       0       5 
```
