# Create a shadow color

Generates a shadow variant of a color by reducing its value (brightness)
and slightly desaturating it. This creates realistic shadow effects for
geometric patterns and 3D visualizations.

## Usage

``` r
shadow_color(color, intensity = 0.3)
```

## Arguments

- color:

  Hex color code.

- intensity:

  Numeric value between 0 and 1 specifying the shadow intensity. Higher
  values produce darker, more desaturated shadows. Default is 0.3.

## Value

A hex color code representing the shadow variant.

## Details

Shadows are created by converting the input color from RGB to HSV color
space as defined by Smith (1978), where the three components, hue,
saturation, and value, encode chromatic identity, color purity, and
luminance independently. The function reduces the value (V) component
proportionally to the specified intensity, simulating the physical
attenuation of light in shadow regions. Simultaneously, saturation is
reduced by 30 percent of the intensity parameter, reflecting the
observation that shadows in natural scenes appear less chromatically
saturated than directly illuminated surfaces; this desaturation
approximates the contribution of ambient (diffuse) light that dilutes
spectral purity in occluded areas. The resulting HSV coordinates are
converted back to hexadecimal RGB representation. This approach produces
perceptually convincing shadows because the hue channel is preserved
intact, ensuring that the shadow retains the chromatic identity of the
original surface color.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- "#1E4D8C"  # Alhambra blue

# Create shadows of varying intensity
light_shadow <- shadow_color(base, 0.2)
medium_shadow <- shadow_color(base, 0.4)
dark_shadow <- shadow_color(base, 0.6)

# Compare base and shadow
plot_palette(c(base, shadow_color(base, 0.3)))
} # }
```
