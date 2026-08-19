# Plot a Samkhya evolution state

Renders the current state in one of three views: the proportions of the
three qualities, the progress of manifestation across the structural
classes of Samkhya-Karika 3, or the branching derivation of
Samkhya-Karika 22 and 25 as a layered graph.

## Usage

``` r
# S3 method for class 'samkhya_state'
plot(x, type = c("gunas", "tattvas", "hierarchy"), ...)
```

## Arguments

- x:

  An object of class `samkhya_state`.

- type:

  One of `"gunas"` (default), `"tattvas"` or `"hierarchy"`.

- ...:

  Further arguments, ignored.

## Value

A `ggplot` object, returned invisibly.

## Details

The default view is the guna simplex because the three qualities are the
quantity the evolution functions actually change; the other two views
are structural and change only in which principles are shaded as
manifest.

The `"hierarchy"` view draws the derivation as the karika states it:
manas and the ten organs on one branch from ahamkara and the five
tanmatras on another, issuing in parallel rather than in sequence, with
the gross elements descending from the tanmatras alone.

## References

Isvarakrsna. *Samkhyakarika*, verses 3, 12-13, 22, 25.

## See also

[`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md),
[`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md),
[`create_samkhya_flowchart()`](https://robustecologies.github.io/SamkhyaR/reference/create_samkhya_flowchart.md)
for the Mermaid rendering of the same graph.

## Examples

``` r
s <- samkhya_evolve(verbose = FALSE)
plot(s)

plot(s, type = "tattvas")

plot(s, type = "hierarchy")

```
