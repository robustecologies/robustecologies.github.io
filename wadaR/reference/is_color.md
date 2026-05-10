# Check if a string is a valid color specification

Validates whether a string represents a valid color in R. Accepts hex
colors (e.g., "#FF0000"), named colors (e.g., "red"), and other valid
color specifications.

## Usage

``` r
is_color(x)
```

## Arguments

- x:

  Character vector. Strings to check as color specifications.

## Value

Logical vector. TRUE for each valid color, FALSE otherwise.

## Examples

``` r
is_color("#FF0000")
#> [1] TRUE
is_color("red")
#> [1] TRUE
is_color("not_a_color")
#> [1] FALSE
is_color(c("#00FF00", "blue", "invalid"))
#> [1]  TRUE  TRUE FALSE
```
