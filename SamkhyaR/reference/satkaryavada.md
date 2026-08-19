# The doctrine that the effect pre-exists in its cause

Returns the five arguments of Samkhya-Karika 9 for *satkaryavada*, the
doctrine that the effect exists in its material cause before it is
produced.

## Usage

``` r
satkaryavada()
```

## Value

An object of class `samkhya_satkarya`: a list with a `table` component
of five rows, each carrying the karika's term, a literal rendering and a
gloss.

## Details

Samkhya-Karika 9 reads *asadakaranad upadanagrahanat sarvasambhava-
bhavat / saktasya sakyakaranat karanabhavac ca sat karyam*. The effect
is existent, because what does not exist cannot be made to be; because a
material cause is resorted to; because not everything arises from
everything; because what is capable produces only what it is capable of;
and because the effect is of the nature of the cause.

The doctrine is what makes the whole evolutionary sequence intelligible.
If the effect did not pre-exist in the cause, the twenty-three evolutes
would be additions to Prakrti rather than manifestations of what she
already contains, and SK 15 and 16's inference from the finitude and
homogeneity of the manifest back to an unmanifest cause would not go
through.

Samkhya's transformation is real, *parinama*: the effect is a genuine
modification of the cause. This distinguishes it from the *vivarta* of
Advaita Vedanta, in which the effect is an appearance that leaves the
cause unaltered, and from the *asatkaryavada* of Nyaya-Vaisesika, in
which the effect is a genuinely new entity.

## References

Isvarakrsna. *Samkhyakarika*, verses 9, 15-16.

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

## See also

[`samkhya_pramanas()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_pramanas.md)
for the means by which the doctrine is established,
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md),
[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md).

## Examples

``` r
satkaryavada()
#> <samkhya_satkarya>
#> Satkaryavada: the effect exists in its cause before production (SK 9)
#> 
#> 1. asadakaranat
#>    Because there is no making of what is not
#>    Production out of nothing is not production
#> 2. upadanagrahanat
#>    Because a material cause is resorted to
#>    The potter takes clay, not milk, to make a pot
#> 3. sarvasambhavabhavat
#>    Because there is no arising of everything from everything
#>    Oil comes from sesame and not from sand
#> 4. saktasya sakyakaranat
#>    Because what is capable produces what it is capable of
#>    A capacity is a capacity for something determinate
#> 5. karanabhavat
#>    Because the effect is of the nature of the cause
#>    Cloth is not other than the threads that constitute it
#> 
#> Transformation is real (parinama), not apparent (vivarta).

satkaryavada()$table$term
#> [1] "asadakaranat"          "upadanagrahanat"       "sarvasambhavabhavat"  
#> [4] "saktasya sakyakaranat" "karanabhavat"         
```
