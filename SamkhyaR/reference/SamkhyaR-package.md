# SamkhyaR: computational exploration of the Samkhya darsana

Tools for working with the enumeration and the causal derivation of
classical Samkhya, one of the six orthodox schools of Indian philosophy.
The package carries the twenty-five principles as a documented data
object, simulates the branching evolution of Prakrti as the
Samkhyakarika states it, renders the derivation as a diagram, and
provides reference objects for the doctrines that surround the
evolution: the means of valid knowledge, the dispositions of buddhi, the
subtle body, and the pre-existence of the effect in its cause.

## Details

Samkhya distinguishes two irreducible realities. Purusha is
consciousness: witness, neutral, spectator, and non-agent (SK 19).
Prakrti is the unconscious material principle, constituted by the three
qualities sattva, rajas and tamas, which in equilibrium constitute the
unmanifest and in disequilibrium the whole manifest world (SK 10-13).
From their conjunction the twenty-three evolutes of Prakrti proceed.

Three points of doctrine govern the design of this package and are
enforced by its tests.

The partition of SK 3 is one uncreated root, seven principles that are
at once effects and causes, sixteen that are effects only, and Purusha
who is neither: \\1 + 7 + 16 + 1 = 25\\. Prakrti therefore has
twenty-three evolutes, since mulaprakrti is not an evolute of itself.

The derivation of SK 22 and 25 branches. From the sattvic aspect of
ahamkara proceed the eleven organs, manas among them; from its tamasic
aspect proceed the five tanmatras; the rajasic aspect supplies the
activity of both and emits nothing. The two branches are parallel, so
the tanmatras do not descend from the organs, and the organs of action
are sattvic in origin.

Nothing happens to Purusha. SK 62 denies that anyone is bound, released
or a transmigrant, and assigns bondage and release to Prakrti in her
manifold supports. Discriminative knowledge arises in buddhi; what ends
at kaivalya is Prakrti's activity on behalf of that Purusha (SK 59,
66-68). Accordingly no function in this package alters the description
of Purusha carried by a `samkhya_state`.

Classical Samkhya is non-theistic. The karika argues for Prakrti and
Purusha without positing a creator, and the theistic reading belongs to
Vijnanabhiksu's sixteenth-century synthesis rather than to the karika
itself.

## Note

This package is the original creation of the author in all conceptual,
theoretical and design aspects. Implementation was assisted by
Anthropic's Claude Code to streamline package development. All original
ideas, creativity and scientific contributions belong to the author, who
maintains full responsibility for the package's correctness and
reliability. Users are encouraged to report any issue through the
package's issue tracker.

## Enumeration

- Data:

  [tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md),
  the twenty-five principles with their structural class, material
  cause, and branch membership.

## Evolution

- Sequence:

  [`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md),
  [`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md),
  [`evolution_ahamkara()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_ahamkara.md),
  [`evolution_manas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_manas.md),
  [`evolution_indriyas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_indriyas.md),
  [`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md),
  [`evolution_mahabhutas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_mahabhutas.md),
  and
  [`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
  for the whole branching sequence in one call.

- Outcome:

  [`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md).

- Methods:

  [`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md),
  [`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md),
  [`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md).

## Doctrinal reference

- Objects:

  [`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md),
  [`samkhya_pramanas()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_pramanas.md),
  [`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md),
  [`satkaryavada()`](https://robustecologies.github.io/SamkhyaR/reference/satkaryavada.md).

## Visualisation

- Diagrams:

  [`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md),
  [`create_samkhya_flowchart_detailed()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart_detailed.md),
  [`print_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/print_samkhya_flowchart.md).

## References

Isvarakrsna. *Samkhyakarika*. Composed c. 350-450 CE; the received text
transmits seventy-two verses, although verse 72 counts the work at
seventy.

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

Larson, G. J., & Bhattacharya, R. S. (Eds.). (1987). *Samkhya: A dualist
tradition in Indian philosophy*. Encyclopedia of Indian Philosophies,
Vol. 4. Princeton, NJ: Princeton University Press.
[doi:10.1515/9781400853533](https://doi.org/10.1515/9781400853533)

## See also

Useful links:

- <https://github.com/robustecologies/SamkhyaR>

- <https://robustecologies.github.io/SamkhyaR>

- Report bugs at <https://github.com/robustecologies/SamkhyaR/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))

Authors:

- Pablo Almaraz <pablo.almaraz@csic.es>
  ([ORCID](https://orcid.org/0000-0003-1416-2695))

## Examples

``` r
# The enumeration, and the partition of Samkhya-Karika 3.
table(tattvas$category)
#> 
#>               avikriti na-prakriti-na-vikriti       prakriti-vikriti 
#>                      1                      1                      7 
#>                 vikara 
#>                     16 

# The complete branching evolution.
state <- samkhya_evolve(verbose = FALSE)
summary(state)
#> <summary.samkhya_state>
#> 
#> Stage: Mahabhutas (the specific)
#> Cognitive situation: Visesa (the specific): the gross world stands manifest
#> 
#> Qualities of Prakrti (SK 12-13)
#>   Sattva = 0.4333   Rajas = 0.2833   Tamas = 0.2833   Sum = 1.0000
#> 
#> Structural partition of what has manifested (SK 3)
#>   Avikriti, uncreated root                : 1 of 1
#>   Prakriti-vikriti, effect and cause both : 7 of 7
#>   Vikara, effect only                     : 16 of 16
#>   Evolutes of Prakrti manifest            : 23 of 23
#> 
#> Branches issuing from ahamkara (SK 25)
#>   Vaikrta, sattvic: the eleven organs     : 11 of 11
#>   Bhutadi, tamasic: the five tanmatras    : 5 of 5
#>   Taijasa, rajasic: supplies activity to both, issues no tattva
#> 
#> Internal organ, threefold (SK 33): mahat, ahaṃkāra, manas
#> Subtle body, eighteenfold (SK 40): 18 of 18 constituents
#> 
#> Viveka in buddhi: Not attained
#> Prakrti acting for this Purusha: Yes
#> Purusha: Unchanged throughout; neither bound nor released (SK 62)

# Discriminative knowledge arises in buddhi; Purusha is unaffected.
final <- get_kaivalya(state, verbose = FALSE)
identical(final$purusha, state$purusha)
#> [1] TRUE
```
