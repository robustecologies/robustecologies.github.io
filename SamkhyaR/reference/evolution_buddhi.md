# Evolve buddhi, the cosmic intellect

Manifests *mahat*, the great one, which is buddhi considered in its
cosmic aspect. It is the first evolute of Prakrti and the second
principle of the enumeration. Samkhya-Karika 23 defines it by
*adhyavasaya*, ascertainment or determination, and assigns to it eight
dispositions, four sattvic and four their reverse.

## Usage

``` r
evolution_buddhi(
  state,
  guna_dominant = c("sattva", "rajas", "tamas"),
  delta = 0.1,
  verbose = TRUE
)
```

## Arguments

- state:

  A `samkhya_state`, ordinarily from
  [`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md).

- guna_dominant:

  Which quality comes to predominate when the equilibrium breaks. One of
  `"sattva"` (the default), `"rajas"` or `"tamas"`. Sattva predominates
  in buddhi, whose nature is illumination.

- delta:

  Magnitude of the disturbance, in \\(0, 1)\\. The dominant quality
  gains `delta`, and the other two lose it in proportion to their
  current shares, so the triple remains on the two-simplex. Default 0.1.

- verbose:

  Emit a progress message. Default `TRUE`.

## Value

The `samkhya_state` with mahat manifested and the qualities
redistributed.

## Details

The disturbance is a redistribution, not an increase. Prakrti is
constituted by the three qualities and by nothing besides, so a quality
that comes to predominate does so at the expense of the other two; the
three continue to sum to one. Implementations that add to one quality
without subtracting from the others leave the triple off the simplex and
no longer describe a Samkhya state at all.

Samkhya-Karika 23 gives buddhi's sattvic form as *dharma* (virtue),
*jnana* (knowledge), *viraga* (dispassion) and *aisvarya* (mastery), and
states that the tamasic form is the reverse of this. The reverse of
*viraga* is *raga*, attachment, not the coinage *avairagya* that
circulates in secondary accounts. See
[`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md)
for the eight dispositions and the fiftyfold intellectual creation that
rests on them.

Evolution begins because Purusha is present, not because Purusha acts.
SK 21 compares the conjunction to that of a lame man and a blind man,
who cooperate without either becoming the other; SK 62 denies that
Purusha is ever bound, released or a transmigrant.

## References

Isvarakrsna. *Samkhyakarika*, verses 21-23, 62.

Vacaspati Misra. *Tattvakaumudi*, commentary on verse 23.

## See also

[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md),
[`evolution_ahamkara()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_ahamkara.md)
for the next evolute,
[`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md)
for the dispositions of SK 23 and 43,
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md),
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md).

## Examples

``` r
state <- evolution_buddhi(init_prakriti(), guna_dominant = "sattva")
#> Mahat evolves. The equilibrium of the qualities breaks, sattva predominating (SK 22-23).
state$gunas
#>    sattva     rajas     tamas 
#> 0.4333333 0.2833333 0.2833333 

# The disturbance redistributes; it does not create quality.
sum(state$gunas)
#> [1] 1

# A rajasic intellect is a different distribution of the same one unit.
evolution_buddhi(init_prakriti(), guna_dominant = "rajas", verbose = FALSE)$gunas
#>    sattva     rajas     tamas 
#> 0.2833333 0.4333333 0.2833333 
```
