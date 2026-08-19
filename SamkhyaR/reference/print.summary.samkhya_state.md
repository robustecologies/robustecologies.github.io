# Print the summary of a Samkhya evolution state

Print the summary of a Samkhya evolution state

## Usage

``` r
# S3 method for class 'summary.samkhya_state'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.samkhya_state`.

- ...:

  Further arguments, ignored.

## Value

`x`, invisibly.

## See also

[`summary.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/summary.samkhya_state.md),
[`print.samkhya_state()`](https://robustecologies.github.io/SamkhyaR/reference/print.samkhya_state.md).

## Examples

``` r
print(summary(samkhya_evolve(verbose = FALSE)))
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
```
