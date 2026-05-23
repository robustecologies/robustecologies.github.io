# Generate a coloring scheme for patterns

Creates a sequence of colors for tiling patterns using various
distribution methods. This function is designed for geometric art and
Islamic-inspired patterns where color assignment follows specific rules.

## Usage

``` r
pattern_colors(
  n_tiles,
  palette = "alhambra_nazari",
  method = c("cyclic", "random", "gradient", "alternating"),
  seed = NULL
)
```

## Arguments

- n_tiles:

  Integer specifying the number of tiles (elements) to color.

- palette:

  Either a character string specifying a palette name, or a character
  vector of hex color codes. Default is "alhambra_nazari".

- method:

  Character string specifying the coloring method:

  - `"cyclic"`: Colors repeat in sequence (default)

  - `"random"`: Colors assigned randomly from palette

  - `"gradient"`: Smooth interpolation across all tiles

  - `"alternating"`: Alternates between first two colors only

- seed:

  Integer for random seed (only used with method = "random") to ensure
  reproducibility.

## Value

A character vector of `n_tiles` hex color codes.

## Details

Colors are selected from the specified palette using one of four
sampling strategies, each serving a distinct aesthetic purpose in
geometric pattern design. The cyclic method repeats the palette sequence
modularly across all tiles, producing the regular, symmetry-preserving
repetitions characteristic of Islamic geometric art and traditional
Alhambra zellige work. The random method samples from the palette with
replacement using an optional reproducible seed, yielding organic
variation that breaks strict periodicity while maintaining the overall
chromatic character of the source palette. The gradient method
interpolates smoothly between palette colors using
[`colorRampPalette`](https://rdrr.io/r/grDevices/colorRamp.html),
distributing the resulting continuous color ramp evenly across the tile
sequence; this is particularly effective for large patterns where a
progressive color transition conveys spatial structure or
directionality. The alternating method restricts selection to the first
two palette colors in strict alternation, creating checkerboard-like
contrast patterns that emphasize the underlying geometric tessellation.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
if (FALSE) { # \dontrun{
# Cyclic coloring for a 20-tile pattern
cyclic <- pattern_colors(20, "alhambra_nazari", method = "cyclic")

# Random coloring with reproducible seed
random <- pattern_colors(20, "alhambra_nazari", method = "random", seed = 42)

# Gradient across tiles
gradient <- pattern_colors(50, "viridis", method = "gradient")
plot_palette(gradient)

# Alternating two-color pattern
checker <- pattern_colors(20, "alhambra_minimal_blue_white", method = "alternating")
} # }
```
