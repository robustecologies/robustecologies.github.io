# With four parameters I can fit an elephant, and with five I can make him wiggle his trunk

This function generates and plots the shape of an elephant using a
Fourier series with four complex parameters, as described by Mayer,
Khairy, and Howard. It can create either a static plot or an animated
GIF with trunk wiggle animation. The function is named after John von
Neumann, the father of computer science.

## Usage

``` r
VonNeumanns_elephant(
  elephant_color = "darkblue",
  line_width = 2,
  trunk_wiggle = 0,
  wiggle_phase = 0,
  animate = FALSE,
  filename = "wiggling_janus.gif",
  n_frames = 30,
  width = 800,
  height = 600,
  fps = 10
)
```

## Arguments

- elephant_color:

  A character string specifying the color for the elephant's body and
  eye. Default is "darkblue".

- line_width:

  A numeric value specifying the line width (thickness) of the
  elephant's outline. Default is 2.

- trunk_wiggle:

  A numeric value controlling the amplitude of trunk wiggle animation.
  Default is 0 (no wiggle). Positive values create upward wiggle,
  negative values downward. Only used when animate=TRUE or for static
  plots with specific wiggle phase.

- wiggle_phase:

  A numeric value representing the phase of the wiggle animation (0 to
  2\*pi). Default is 0. Only used for static plots with trunk wiggle.

- animate:

  A logical value indicating whether to create an animated GIF instead
  of a static plot. Default is FALSE (creates static plot).

- filename:

  A character string for the output GIF filename when animate=TRUE.
  Default is "elephant_wiggle.gif". Ignored when animate=FALSE.

- n_frames:

  An integer specifying the number of frames in the animation when
  animate=TRUE. Default is 30. Ignored when animate=FALSE.

- width:

  An integer specifying GIF width in pixels when animate=TRUE. Default
  is 800. Ignored when animate=FALSE.

- height:

  An integer specifying GIF height in pixels when animate=TRUE. Default
  is 600. Ignored when animate=FALSE.

- fps:

  An integer specifying frames per second for GIF animation when
  animate=TRUE. Default is 10. Must be a factor of 100 (valid values: 1,
  2, 4, 5, 10, 20, 25, 50, 100). Ignored when animate=FALSE.

## Value

When animate=FALSE, creates a static plot and returns NULL invisibly.
When animate=TRUE, creates a GIF file and returns the filename.

## Details

The four-parameter elephant of Mayer, Khairy and Howard (2010)
parametrises the body outline as a closed Fourier curve in the complex
plane, \$\$ z(t) \\=\\ \sum\_{k=-\infty}^{\infty} C_k \\ e^{\mathrm{i} k
t}, \qquad t \in \[0, 2\pi\], \$\$ where only four non-zero complex
coefficients \\C\_{\pm 1}, C\_{\pm 2}, C\_{\pm 3}, C\_{\pm 5}\\ are
needed for a recognisable elephant silhouette: each \\C_k = a_k +
\mathrm{i} b_k\\ contributes a harmonic of frequency \\k\\ with phase
set by the ratio \\b_k / a_k\\. Mayer, Khairy and Howard chose \\(a_1,
b_1) = (50, -30)\\, \\(a_2, b_2) = (18, 8)\\, \\(a_3, b_3) = (12,
-10)\\, \\(a_5, b_5) = (-14, 0)\\, and used \\b_5\\ as the single extra
parameter that "makes the trunk wiggle". The implementation splits the
real and imaginary parts into \\x(t)\\ and \\y(t)\\ separately and
samples the curve on a uniform grid of \\t\\ values. The wiggle is
obtained by letting \\b_5 \mapsto b_5 + A \sin(\phi)\\ with amplitude
\\A\\ and phase \\\phi\\ sweeping \\\[0, 2\pi)\\ across animation
frames. The function is included as an Easter-egg illustration of the
apocryphal John von Neumann quote that gives the title and of Enrico
Fermi's reply to Freeman Dyson during a 1953 meeting at Columbia
University, where the quote was attributed.

## References

Mayer, J., Khairy, K., & Howard, J. (2010). Drawing an elephant with
four complex parameters. *American Journal of Physics*, 78(6), 648-649.
[doi:10.1119/1.3254017](https://doi.org/10.1119/1.3254017)

Dyson, F. (2004). A meeting with Enrico Fermi. *Nature*, 427(6972), 297.
[doi:10.1038/427297a](https://doi.org/10.1038/427297a)

## See also

[`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md) -
the canonical model constructor of janos;
[`dyn_sim()`](https://robustecologies.github.io/janos/reference/dyn_sim.md) -
numerical integrator for janos `model_spec` objects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic static elephant with default settings
VonNeumanns_elephant()

# Static elephant with custom color and line width
VonNeumanns_elephant(elephant_color = "red", line_width = 3)

# Static elephant with trunk wiggle at specific phase
VonNeumanns_elephant(trunk_wiggle = 5, wiggle_phase = pi / 4)

# Static elephant with strong trunk wiggle
VonNeumanns_elephant(
  trunk_wiggle = 10, wiggle_phase = pi / 2,
  elephant_color = "purple", line_width = 2.5
)

# Create animated GIF with default settings
VonNeumanns_elephant(animate = TRUE)

# Create custom animated GIF
VonNeumanns_elephant(
  animate = TRUE,
  filename = "my_elephant.gif",
  trunk_wiggle = 8,
  n_frames = 40,
  fps = 20,
  elephant_color = "darkgreen"
)
} # }
```
