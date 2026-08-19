# The three means of valid knowledge

Returns the three *pramanas* that Samkhya admits, in the karika's own
vocabulary, with the definitions of Samkhya-Karika 4 to 6.

## Usage

``` r
samkhya_pramanas()
```

## Value

An object of class `samkhya_pramana`: a list with a `table` component
naming each means, its karika term, the term used by the commentators,
its definition and its verse.

## Details

Samkhya-Karika 4 reads *drstam anumanam aptavacanam ca
sarvapramanasiddhatvat / trividham pramanam istam prameyasiddhih
pramanad dhi*. The three are *drsta*, what is seen, *anumana*,
inference, and *aptavacana*, the word of a reliable authority.

The karika's terms are not the terms usually quoted. *Pratyaksa* and
*sabda*, which appear in most secondary accounts of Samkhya
epistemology, are the vocabulary of Nyaya, adopted by the Samkhya
commentators; the karika itself writes *drsta* and *aptavacana*. Both
sets are reported here so that the difference is visible.

SK 5 defines perception as *prativisayadhyavasaya*, the determination of
each object by its own sense; inference is said to be threefold and to
rest on a mark and on the bearer of the mark; and reliable testimony is
*aptasruti*. SK 6 adds that what lies beyond the senses is established
by inference from the general, and what even that cannot reach is
established by reliable testimony, which is how Samkhya claims to
establish Prakrti and Purusha at all.

Samkhya admits three and not more. Nyaya adds comparison, and Mimamsa
adds presumption and non-apprehension; SK 4's *trividham* is a
deliberate restriction, since the remaining means are held to reduce to
these three.

## References

Isvarakrsna. *Samkhyakarika*, verses 4-6.

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

## See also

[`satkaryavada()`](https://robustecologies.github.io/SamkhyaR/reference/satkaryavada.md)
for the causal doctrine these means establish,
[`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md),
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md).

## Examples

``` r
samkhya_pramanas()
#> <samkhya_pramana>
#> Three means of valid knowledge, and three only (SK 4)
#> 
#> 1. Perception (karika: drsta; commentators: pratyaksa)
#>    Determination of each object by the sense proper to it (SK 5)
#> 2. Inference (karika: anumana; commentators: anumana)
#>    Threefold, resting on a mark and on the bearer of the mark (SK 5)
#> 3. Reliable testimony (karika: aptavacana; commentators: sabda)
#>    The word of one who is reliable; reaches what inference cannot (SK 5-6)
#> 
#> Nyaya adds comparison and Mimamsa adds two more; SK 4 restricts to three.

samkhya_pramanas()$table[, c("karika_term", "commentarial_term")]
#>   karika_term commentarial_term
#> 1       drsta         pratyaksa
#> 2     anumana           anumana
#> 3  aptavacana             sabda
```
