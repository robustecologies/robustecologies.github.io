# The eight dispositions of buddhi and the fiftyfold intellectual creation

Returns the eight *bhavas* of Samkhya-Karika 23 and 43, and optionally
the fiftyfold *pratyayasarga* of SK 46 to 51 that rests on them.

## Usage

``` r
buddhi_bhavas(pratyayasarga = TRUE)
```

## Arguments

- pratyayasarga:

  Logical; include the fiftyfold intellectual creation. Default `TRUE`.

## Value

An object of class `samkhya_bhava`: a list with components `bhavas` (a
data frame of eight rows) and, when requested, `pratyayasarga` (a data
frame of four rows whose counts sum to fifty).

## Details

Samkhya-Karika 23 gives buddhi's sattvic form as *dharma*, *jnana*,
*viraga* and *aisvarya*, and states that its tamasic form is the
reverse, *tamasam asmad viparyastam*. The reverse of *viraga*,
dispassion, is *raga*, attachment. The form *avairagya*, which
circulates widely, is a coinage: the karika negates the disposition
rather than prefixing its abstract noun.

SK 44 and 45 give what each disposition effects. By virtue there is
ascent, by vice descent, by knowledge release and by its reverse
bondage: *jnanena capavargo viparyayad isyate bandhah*. By dispassion
there is dissolution into Prakrti, by rajasic passion transmigration, by
mastery non-obstruction and by its reverse the contrary.

The fiftyfold intellectual creation of SK 46 to 51 comprises
misconception in five kinds, incapacity in twenty-eight, contentment in
nine and attainment in eight. SK 47 and 48 name the five misconceptions
as *tamas*, *moha*, *mahamoha*, *tamisra* and *andhatamisra*, in eight,
eight, ten, eighteen and eighteen subdivisions respectively. These are
not the five afflictions of Yoga-sutra 2.3, *avidya*, *asmita*, *raga*,
*dvesa* and *abhinivesa*; the identification of the two lists is made by
the commentators, Vacaspati Misra among them, and is not stated by the
karika.

## References

Isvarakrsna. *Samkhyakarika*, verses 23, 43-51.

Vacaspati Misra. *Tattvakaumudi*, commentary on verses 23 and 47.

## See also

[`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md)
for the principle these dispositions inhere in,
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md)
for the body that carries them,
[`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md),
[`samkhya_pramanas()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_pramanas.md).

## Examples

``` r
buddhi_bhavas()
#> <samkhya_bhava>
#> Eight dispositions of buddhi (SK 23, 43-45)
#> 
#>   Sattvic form
#>     dharma      Virtue       Ascent to the higher worlds (SK 44)
#>     jnana       Knowledge    Release (SK 44)
#>     viraga      Dispassion   Dissolution into Prakrti (SK 45)
#>     aisvarya    Mastery      Non-obstruction (SK 45)
#>   Tamasic form
#>     adharma     Vice         Descent to the lower worlds (SK 44)
#>     ajnana      Ignorance    Bondage (SK 44)
#>     raga        Attachment   Transmigration (SK 45)
#>     anaisvarya  Impotence    Obstruction (SK 45)
#> 
#> Intellectual creation, fiftyfold (SK 46-51)
#>   viparyaya   Misconception   5  tamas, moha, mahamoha, tamisra, andhatamisra (SK 47-48)
#>   asakti      Incapacity     28  Injuries of the eleven organs and seventeen failures of buddhi (SK 49)
#>   tusti       Contentment     9  Four internal and five external (SK 50)
#>   siddhi      Attainment      8  Reasoning, hearing, study, three suppressions of pain, and two more (SK 51)
#>               Total          50
#> 
#> The five misconceptions are not the five afflictions of Yoga-sutra 2.3;
#> that identification belongs to the commentators, not to the karika.

# The four sattvic dispositions and their reverses.
buddhi_bhavas()$bhavas[, c("bhava", "form", "effect")]
#>        bhava     form                              effect
#> 1     dharma sattvika Ascent to the higher worlds (SK 44)
#> 2      jnana sattvika                     Release (SK 44)
#> 3     viraga sattvika    Dissolution into Prakrti (SK 45)
#> 4   aisvarya sattvika             Non-obstruction (SK 45)
#> 5    adharma   tamasa Descent to the lower worlds (SK 44)
#> 6     ajnana   tamasa                     Bondage (SK 44)
#> 7       raga   tamasa              Transmigration (SK 45)
#> 8 anaisvarya   tamasa                 Obstruction (SK 45)

# The intellectual creation is fiftyfold.
sum(buddhi_bhavas()$pratyayasarga$n)
#> [1] 50
```
