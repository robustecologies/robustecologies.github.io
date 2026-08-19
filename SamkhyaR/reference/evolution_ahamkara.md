# Evolve ahamkara, the principle of self-appropriation

Manifests *ahamkara*, the I-maker, the third principle. Samkhya-Karika
24 defines it by *abhimana*, appropriation or self-reference. It is the
point at which the evolution divides, since SK 25 assigns three aspects
to it, each with a different issue.

## Usage

``` r
evolution_ahamkara(state, verbose = TRUE)
```

## Arguments

- state:

  A `samkhya_state` in which mahat has already manifested.

- verbose:

  Emit a progress message. Default `TRUE`.

## Value

The `samkhya_state` with ahamkara manifested.

## Details

Samkhya-Karika 25 reads: *sattvika ekadasakah pravartate vaikrtad
ahamkarat / bhutades tanmatrah sa tamasas taijasad ubhayam*. The
elevenfold set proceeds from the sattvic, called *vaikrta*, aspect; the
tanmatras proceed from the tamasic, called *bhutadi*; and both proceed
from the rajasic, called *taijasa*, in the sense that the rajasic aspect
supplies the activity by which either can proceed at all. The taijasa
aspect emits no principle of its own.

Two consequences are frequently lost in secondary accounts and are
enforced by this package. The first is that the elevenfold set comprises
manas *and* the ten organs of knowledge and action together, so the
organs of action are sattvic in origin and not tamasic or rajasic. The
second is that the two branches issue *in parallel*: the tanmatras
descend from ahamkara directly and not from the organs, so
[`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md)
requires only that ahamkara be present.

## References

Isvarakrsna. *Samkhyakarika*, verses 22, 24-25.

Gaudapada. *Samkhyakarikabhasya*, commentary on verse 25.

## See also

[`evolution_buddhi()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_buddhi.md)
for its material cause,
[`evolution_manas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_manas.md)
and
[`evolution_indriyas()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_indriyas.md)
for the vaikrta branch,
[`evolution_tanmatras()`](https://robustecologies.github.io/SamkhyaR/reference/evolution_tanmatras.md)
for the bhutadi branch,
[`samkhya_evolve()`](https://robustecologies.github.io/SamkhyaR/reference/samkhya_evolve.md).

## Examples

``` r
state <- evolution_ahamkara(evolution_buddhi(init_prakriti(), verbose = FALSE),
                            verbose = FALSE)
summary(state)$antahkarana
#> [1] "mahat"    "ahaṃkāra"

# Both branches are now open, and either may be taken first.
summary(state)$branches
#> vaikrta bhutadi 
#>       0       0 
```
