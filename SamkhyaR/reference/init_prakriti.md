# Initialise Prakrti in its primordial equilibrium

Creates the unmanifest condition of Prakrti, in which the three
qualities stand in perfect equilibrium and nothing has yet evolved. This
is the *avyakta*, the unmanifest, of Samkhya-Karika 10 and 11, also
called *pradhana*, the primary. Purusha is present as witness
throughout, and the state records it as a constant that no function of
this package alters.

## Usage

``` r
init_prakriti(sattva = 1/3, rajas = 1/3, tamas = 1/3)
```

## Arguments

- sattva, rajas, tamas:

  Initial proportions of the three qualities. They must be non-negative
  and must sum to one; the defaults place them in exact equilibrium at
  one third each.

## Value

An object of class `samkhya_state` carrying the guna triple, the single
manifested principle (mulaprakrti itself), the stage of evolution, the
cognitive situation, and the constant description of Purusha.

## Details

Samkhya-Karika 10 contrasts the manifest, which is caused, non-eternal,
non-pervasive, active, plural, dependent, mergent, composite and
subordinate, with the unmanifest, which is the reverse in each respect.
SK 11 adds that both the manifest and the unmanifest are constituted by
the three qualities, are non-discriminating, objective, common,
unconscious and productive, and that Purusha is the reverse of this and
yet like the unmanifest in some of the properties of SK 10.

The equilibrium of the three qualities in the unmanifest condition is
called *samyavastha*. The term is not the karika's own; it belongs to
the later Samkhya of the *Samkhyasutra* and of Vacaspati Misra, and it
is used here because it names the condition compactly and without
ambiguity.

The three qualities are proportions of one substance and therefore sum
to one. Every function in this package preserves that invariant, which
is verified whenever a state is validated.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 10-11, 15-16.

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

## See also

[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
for the complete branching sequence,
[`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md)
for the first evolute,
[`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md),
[`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md),
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md),
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md)
for the enumeration this state indexes into.

## Examples

``` r
state <- init_prakriti()
state
#> <samkhya_state>
#> Stage: Avyakta (unmanifest)
#> Gunas: sattva = 0.333, rajas = 0.333, tamas = 0.333 (sum = 1.000)
#> Evolutes manifest: 0 of 23
#> Cognitive situation: Samyavastha (equilibrium of the three qualities)
#> Viveka in buddhi: Not attained
#> Prakrti acting for this Purusha: Yes
#> Purusha: Witness, neutral, non-agent; neither bound nor released (SK 19, 62)

# The three qualities are in equilibrium and sum to one.
state$gunas
#>    sattva     rajas     tamas 
#> 0.3333333 0.3333333 0.3333333 
sum(state$gunas)
#> [1] 1

# Nothing has evolved yet: only the root itself is present.
summary(state)$n_manifest
#> [1] 0
```
