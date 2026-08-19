# Print a Samkhya evolution state

Displays, on a single screen, the condition of Prakrti at the current
point of the evolutionary sequence: the proportions of the three
qualities, how many of the twenty-three evolutes have manifested, the
cognitive situation, and whether discriminative knowledge has arisen in
buddhi.

## Usage

``` r
# S3 method for class 'samkhya_state'
print(x, ...)
```

## Arguments

- x:

  An object of class `samkhya_state`.

- ...:

  Further arguments, ignored.

## Value

`x`, invisibly.

## Details

The display never reports Purusha as bound or as liberated.
Samkhya-Karika 62 is explicit that no one is bound, no one is released
and no one transmigrates; it is Prakrti, in her manifold supports, who
is bound and released. What the state records is therefore the condition
of Prakrti and of buddhi, and the constant description of Purusha given
at SK 19.

## References

Isvarakrsna. *Samkhyakarika*, verses 19, 23, 62.

## See also

[`init_prakriti()`](https://robustecologies.github.io/SamkhyaR/reference/init_prakriti.md)
and
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md)
to construct such an object,
[`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md)
for the fuller diagnostic,
[`plot.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/plot.samkhya_state.md)
for the graphical view,
[`get_kaivalya()`](https://robustecologies.github.io/SamkhyaR/reference/get_kaivalya.md)
for the cessation of Prakrti's activity.

## Examples

``` r
init_prakriti()

print(samkhya_evolve(verbose = FALSE))
#> <samkhya_state>
#> Stage: Mahabhutas (the specific)
#> Gunas: sattva = 0.433, rajas = 0.283, tamas = 0.283 (sum = 1.000)
#> Evolutes manifest: 23 of 23
#> Cognitive situation: Visesa (the specific): the gross world stands manifest
#> Viveka in buddhi: Not attained
#> Prakrti acting for this Purusha: Yes
#> Purusha: Witness, neutral, non-agent; neither bound nor released (SK 19, 62)
```
