# Summarise a Samkhya evolution state

Assembles the structural diagnostics of the state: the partition of the
manifested principles under Samkhya-Karika 3, the composition of the two
branches issuing from ahamkara under SK 25, the membership of the
internal organ under SK 33, the constituents of the subtle body under SK
40, and the proportions of the three qualities.

## Usage

``` r
# S3 method for class 'samkhya_state'
summary(object, ...)
```

## Arguments

- object:

  An object of class `samkhya_state`.

- ...:

  Further arguments, ignored.

## Value

An object of class `summary.samkhya_state`: a list with components
`gunas`, `partition` (a table over the SK 3 classes), `branches` (a
table over the vaikrta and bhutadi aspects), `antahkarana`, `linga`,
`n_manifest`, `n_evolutes`, `viveka`, `prakriti_active`, `jivanmukti`
and `stage`.

## Details

The partition reported is that of SK 3: one uncreated root, seven
principles that are at once effects and causes, sixteen that are effects
only, and Purusha who is neither. The branch table reports how many of
the eleven organs of the vaikrta aspect and how many of the five
tanmatras of the bhutadi aspect have so far manifested; because these
are parallel branches, either may be complete while the other is empty.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 25, 33, 40.

## See also

[`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md),
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md),
[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md),
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md),
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md).

## Examples

``` r
s <- samkhya_evolve(verbose = FALSE)
summary(s)
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

summary(s)$partition
#>         avikriti prakriti-vikriti           vikara 
#>                1                7               16 
```
