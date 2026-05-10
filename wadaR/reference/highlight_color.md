# Create a highlight color

Generates a highlight variant of a color by increasing its value
(brightness) and reducing saturation. This creates realistic highlight
effects for geometric patterns and 3D visualizations.

## Usage

``` r
highlight_color(color, intensity = 0.3)
```

## Arguments

- color:

  Hex color code.

- intensity:

  Numeric value between 0 and 1 specifying the highlight intensity.
  Higher values produce brighter, more washed-out highlights. Default is
  0.3.

## Value

A hex color code representing the highlight variant.

## Details

The function converts the input color from RGB to HSV color space as
defined by Smith (1978), then increases the value (V) component toward
its maximum of 1.0 proportionally to the intensity parameter,
approximating the effect of specular reflection on a surface. At the
same time, saturation is reduced by 50 percent of the intensity, which
simulates how highlights on physical surfaces tend to desaturate toward
white as the angle of incidence approaches the specular reflection
angle. This dual adjustment, increasing brightness while reducing color
purity, produces highlights that convincingly approximate the Phong
reflection model's specular component without requiring explicit surface
geometry. Because hue is preserved throughout the transformation, the
highlighted color remains recognizably related to the original surface
color, ensuring visual coherence in geometric pattern rendering and 3D
visualizations.

## References

Smith, A. R. (1978). Color gamut transform pairs. *ACM SIGGRAPH Computer
Graphics*, 12(3), 12-19.
[doi:10.1145/965139.807361](https://doi.org/10.1145/965139.807361)

## Examples

``` r
if (FALSE) { # \dontrun{
base <- "#1E4D8C"  # Alhambra blue

# Create highlights of varying intensity
subtle_highlight <- highlight_color(base, 0.2)
medium_highlight <- highlight_color(base, 0.4)
bright_highlight <- highlight_color(base, 0.6)

# Create a 3D effect with shadow and highlight
plot_palette(c(shadow_color(base), base, highlight_color(base)))
} # }
```
