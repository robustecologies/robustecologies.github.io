# Logical check of a constant

This function provides a logical test of whether or not a vector is a
constant.

## Usage

``` r
is.constant(x)
```

## Arguments

- x:

  A vector.

## Value

A logical result, reporting `TRUE` when a vector is a constant, or
`FALSE` otherwise.

## Details

As opposed to a variable, a constant is a vector in which the elements
contain less than or equal to one unique value.

## See also

[`unique`](https://rdrr.io/r/base/unique.html)

## Examples

``` r
if (FALSE) { # \dontrun{
is.constant(rep(1, 10)) # TRUE
is.constant(1:10) # FALSE
} # }
```
