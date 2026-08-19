# Attain kaivalya through discriminative knowledge

Records the arising of discriminative knowledge in buddhi and the
consequent cessation of Prakrti's activity on behalf of this Purusha.
Samkhya-Karika 62 to 68 describe the outcome: not an attainment by
Purusha, who was never bound, but the withdrawal of a Prakrti whose
purpose has been served.

## Usage

``` r
get_kaivalya(state, jivanmukti = TRUE, verbose = TRUE)
```

## Arguments

- state:

  A `samkhya_state` in which all twenty-three evolutes of Prakrti have
  manifested.

- jivanmukti:

  Logical. If `TRUE` (the default), the state records that the body
  persists after knowledge by the momentum of what has already been set
  going, as a potter's wheel continues to turn (SK 67). If `FALSE`, the
  state records the final separation of SK 68.

- verbose:

  Emit a message. Default `TRUE`.

## Value

The `samkhya_state` with `viveka` set to `TRUE`, `prakriti_active` set
to `FALSE`, and `jivanmukti` set as requested. The `purusha` component
is returned unchanged, because nothing happens to Purusha.

## Details

The central claim of this function is negative, and it is the point on
which computational treatments of Samkhya most often go wrong.
Samkhya-Karika 62 reads: *tasman na badhyate 'ddha na mucyate napi
samsarati kascit / samsarati badhyate mucyate ca nanasraya prakrtih*.
Therefore no one is bound, no one is released, and no one transmigrates;
it is Prakrti, in her manifold supports, who transmigrates, is bound and
is released. To set a field on Purusha to the value "liberated" is to
assert exactly what the verse denies.

What does change is Prakrti and buddhi. Discriminative knowledge arises
in buddhi, not in Purusha, and SK 64 gives its content as *nasmi na me
naham*, I am not, nothing is mine, I am not an I. SK 59 compares Prakrti
to a dancer who withdraws from the stage once she has been seen, and SK
61 calls her the most delicate of beings, who having once been seen
never comes again into the sight of that Purusha. SK 66 states that the
purpose of their conjunction is at an end while neither has ceased to
exist.

Embodiment does not stop at the moment of knowledge. SK 67 states that
the body persists through the momentum of impressions already acquired,
as a potter's wheel continues to turn after the potter has stopped
pushing it; this is the condition the later tradition calls
*jivanmukti*. Only at the separation from the body, Prakrti's purpose
being accomplished, is isolation attained both certainly and finally,
*aikantikam atyantikam* (SK 68).

The function requires all twenty-three evolutes because discrimination
is discrimination of Purusha from the whole of Prakrti; a partial
manifestation leaves something undiscriminated. Prakrti has twenty-three
evolutes, not twenty-four: mulaprakrti is not an evolute of itself (SK
3).

## References

Isvarakrsna. *Samkhyakarika*, verses 2, 3, 59, 61-64, 66-68.

Larson, G. J. (1979). *Classical Samkhya: An interpretation of its
history and meaning* (2nd rev. ed.). Delhi: Motilal Banarsidass. ISBN
978-81-208-0503-3.

## See also

[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
for the manifestation that precedes it,
[`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md),
[`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md),
[`linga_sharira()`](https://robustecologies.github.io/SamkhyaR/reference/linga_sharira.md)
for what transmigrates until this point.

## Examples

``` r
state <- samkhya_evolve(verbose = FALSE)
final <- get_kaivalya(state, verbose = FALSE)

# Discriminative knowledge has arisen and Prakrti has ceased to act.
final$viveka
#> [1] TRUE
final$prakriti_active
#> [1] FALSE

# Purusha is returned exactly as it was: nothing happens to Purusha (SK 62).
identical(final$purusha, state$purusha)
#> [1] TRUE

# An incomplete manifestation is refused.
try(get_kaivalya(init_prakriti()))
#> Error : Discrimination is incomplete: 23 of the 23 evolutes of Prakrti have not manifested (mahat, ahaṃkāra, manas, śrotra, tvac). Viveka discriminates Purusha from the whole of Prakrti.
```
