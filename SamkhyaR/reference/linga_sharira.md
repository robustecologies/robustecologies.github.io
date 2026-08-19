# The subtle body and its eighteen constituents

Returns the constituents of the subtle body, which Samkhya-Karika 40
defines as the eighteen consisting of mahat, ahamkara, the eleven organs
and the five tanmatras. This is the aggregate that transmigrates.

## Usage

``` r
linga_sharira(state = NULL)
```

## Arguments

- state:

  Optionally a `samkhya_state`, in which case the result also reports
  how many of the eighteen have so far manifested.

## Value

An object of class `samkhya_linga`: a list with components
`constituents` (the eighteen IAST names), `table` (the corresponding
rows of
[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md)),
`n` (eighteen) and, when `state` is supplied, `manifest`.

## Details

*Linga-sarira* and *suksma-sarira* are two names for one thing.
Samkhya-Karika 40 reads *purvotpannam asaktam niyatam mahadi-suksma-
paryantam*, arisen before, unattached, constant, extending from mahat to
the subtle ones; SK 39 distinguishes the subtle from the bodies born of
mother and father and from the gross elements. Accounts that make the
two names designate two different aggregates, one running from buddhi to
the tanmatras and another consisting of buddhi, ahamkara and manas,
split a single term of the karika into two and contradict SK 40.

The gross elements are not constituents of the subtle body. The subtle
body transmigrates; the gross body, constituted by the gross elements,
perishes.

SK 41 gives the reason the subtle body must exist: as a picture cannot
stand without a support, nor a shadow without a post, so the subtle
marks cannot subsist without a specific substratum. SK 42 gives its
purpose: it performs its function for the sake of Purusha, like a player
taking on a role, by means of the conjunction of efficient causes and
the dispositions.

## References

Isvarakrsna. *Samkhyakarika*, verses 39-42.

Larson, G. J., & Bhattacharya, R. S. (Eds.). (1987). *Samkhya: A dualist
tradition in Indian philosophy*. Encyclopedia of Indian Philosophies,
Vol. 4. Princeton, NJ: Princeton University Press.
[doi:10.1515/9781400853533](https://doi.org/10.1515/9781400853533)

## See also

[tattvas](https://robustecologies.github.io/SamkhyaR/reference/tattvas.md)
for the `linga` column that flags these principles,
[`buddhi_bhavas()`](https://robustecologies.github.io/SamkhyaR/reference/buddhi_bhavas.md)
for the dispositions the subtle body carries,
[`evolution_mahabhutas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_mahabhutas.md),
[`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md).

## Examples

``` r
linga_sharira()
#> <samkhya_linga>
#> Subtle body, linga-sarira and suksma-sarira being one thing (SK 39-40)
#> 
#>   buddhi       Mahat
#>   ahamkara     Ahamkara
#>   manas        Manas
#>   jnanendriya  Shrotra, Tvac, Chakshus, Rasana, Ghrana
#>   karmendriya  Vac, Pani, Pada, Payu, Upastha
#>   tanmatra     Shabda, Sparsha, Rupa, Rasa, Gandha
#> 
#> Constituents: 18
#> The gross elements are excluded: they constitute the gross body, which perishes.

# Eighteen constituents, none of them a gross element.
length(linga_sharira()$constituents)
#> [1] 18

# Progress within a partly evolved state.
state <- init_prakriti() |>
    evolution_buddhi(verbose = FALSE) |>
    evolution_ahamkara(verbose = FALSE)
linga_sharira(state)$manifest
#> [1] 2
```
