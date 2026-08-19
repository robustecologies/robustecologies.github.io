# The twenty-five tattvas of the Samkhya enumeration

The canonical enumeration of the twenty-five principles (tattva) of
classical Samkhya, carrying for each principle its structural class
under Samkhya-Karika 3, its immediate material cause, the aspect of
ahamkara from which it issues under Samkhya-Karika 25, and its
membership in the internal organ, the elevenfold set and the subtle
body. This table is the single source of truth for the enumeration: the
flowchart builders, the print and summary methods, and the package
vignettes all derive from it, so the numbering cannot drift between
them.

## Usage

``` r
tattvas
```

## Format

A data frame with 25 rows and 13 columns:

- index:

  Integer 1 to 25. The twenty-four material principles are numbered from
  mulaprakrti (1) to prthivi (24) and Purusha is the twenty-fifth
  (pancavimsa), which is the enumeration presupposed by the classical
  commentators and by the Moksadharma.

- iast:

  Sanskrit name in IAST transliteration.

- devanagari:

  Sanskrit name in Devanagari, given in stem form.

- ascii:

  Transliteration without diacritics, for console output.

- english:

  Short English gloss.

- group:

  One of `"mulaprakriti"`, `"buddhi"`, `"ahamkara"`, `"manas"`,
  `"jnanendriya"`, `"karmendriya"`, `"tanmatra"`, `"mahabhuta"`,
  `"purusha"`.

- category:

  Structural class under SK 3: `"avikriti"` (uncreated, mulaprakrti
  alone), `"prakriti-vikriti"` (both creating and created; the seven,
  namely mahat, ahamkara and the five tanmatras), `"vikara"` (created
  only; the sixteen, namely the eleven organs and the five gross
  elements), and `"na-prakriti-na-vikriti"` (neither; Purusha alone).

- parent:

  IAST name of the immediate material cause, or `NA` for mulaprakrti and
  Purusha, which have none.

- guna_aspect:

  `"vaikrta"` for the eleven issuing from the sattvic aspect of
  ahamkara, `"bhutadi"` for the five tanmatras issuing from its tamasic
  aspect, `NA` otherwise (SK 25).

- antahkarana:

  Logical; `TRUE` for buddhi, ahamkara and manas, the threefold internal
  organ of SK 33.

- ekadasaka:

  Logical; `TRUE` for the eleven organs that issue together from the
  vaikrta ahamkara (SK 25).

- linga:

  Logical; `TRUE` for the eighteen constituents of the subtle body (SK
  40).

- sk_verse:

  Karika verse numbers where the principle is treated.

## Source

Isvarakrsna. *Samkhyakarika*, verses 3, 22, 24-28, 33, 38 and 40. Text
and classification as established in Larson (1979) and in the critical
summaries of Larson and Bhattacharya (1987).

## Details

Three structural facts are enforced when the table is built and are
asserted again by the package test suite.

The partition of SK 3 is \\1 + 7 + 16 + 1 = 25\\. Mulaprakrti is
uncreated; mahat, ahamkara and the five tanmatras are simultaneously
effects of what precedes them and causes of what follows; the eleven
organs and the five gross elements are effects only; Purusha is neither
cause nor effect. It follows that Prakrti has exactly **twenty-three**
evolutes, not twenty-four, since mulaprakrti is not an evolute of
itself.

The derivation of SK 22 and 25 is branching, not linear. From ahamkara
the sattvic or *vaikrta* aspect issues the elevenfold set, comprising
manas and the ten organs of knowledge and action; the tamasic or
*bhutadi* aspect issues the five tanmatras. These proceed in parallel:
the tanmatras are not derived from the organs. The rajasic or *taijasa*
aspect supplies the activity by which both proceed and issues no tattva
of its own, which is why no row carries `guna_aspect == "taijasa"`. The
five gross elements then arise from the tanmatras.

The `parent` of each gross element is the corresponding tanmatra alone,
following Vacaspati Misra's one-to-one derivation in the
*Tattvakaumudi*. Vijnanabhiksu instead holds that each gross element
accumulates the qualities of the preceding ones, so that space carries
sound alone and earth carries all five. The accumulation scheme is
frequently attributed to SK 22, which states only the order of
derivation; the reader should treat it as a commentarial position rather
than as the text of the karika.

## References

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

Larson, G. J., & Bhattacharya, R. S. (Eds.). (1987). *Samkhya: A dualist
tradition in Indian philosophy*. Encyclopedia of Indian Philosophies,
Vol. 4. Princeton, NJ: Princeton University Press.
[doi:10.1515/9781400853533](https://doi.org/10.1515/9781400853533)

## See also

[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md)
to begin an evolution from this enumeration,
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
for the complete branching sequence,
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md)
for the eighteen constituents flagged by `linga`,
[`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md)
for the diagram generated from this table.

## Examples

``` r
# The partition of Samkhya-Karika 3.
table(tattvas$category)
#> 
#>               avikriti na-prakriti-na-vikriti       prakriti-vikriti 
#>                      1                      1                      7 
#>                 vikara 
#>                     16 

# Prakrti has twenty-three evolutes, not twenty-four.
sum(tattvas$category %in% c("prakriti-vikriti", "vikara"))
#> [1] 23

# The eleven of SK 25 all issue from the vaikrta aspect of ahamkara.
subset(tattvas, ekadasaka, select = c(index, iast, guna_aspect))
#>    index    iast guna_aspect
#> 4      4   manas     vaikrta
#> 5      5  śrotra     vaikrta
#> 6      6    tvac     vaikrta
#> 7      7  cakṣus     vaikrta
#> 8      8  rasanā     vaikrta
#> 9      9  ghrāṇa     vaikrta
#> 10    10     vāc     vaikrta
#> 11    11    pāṇi     vaikrta
#> 12    12    pāda     vaikrta
#> 13    13    pāyu     vaikrta
#> 14    14 upastha     vaikrta

# The tanmatras issue from the tamasic aspect, in parallel with the eleven.
subset(tattvas, group == "tanmatra", select = c(index, iast, parent, guna_aspect))
#>    index   iast   parent guna_aspect
#> 15    15  śabda ahaṃkāra     bhutadi
#> 16    16 sparśa ahaṃkāra     bhutadi
#> 17    17   rūpa ahaṃkāra     bhutadi
#> 18    18   rasa ahaṃkāra     bhutadi
#> 19    19 gandha ahaṃkāra     bhutadi
```
