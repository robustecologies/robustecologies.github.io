# koloRo: comprehensive color palettes for R

A comprehensive collection of 282 carefully curated color palettes
organized across 16 categories for scientific visualization,
accessibility, geometric art, and cinematic color grading.

## Note

This package is the original creation of the author in all conceptual,
theoretical, and design aspects. Implementation was assisted by
Anthropic's Claude Code v.2 (Opus 4.5) to streamline package
development. All original ideas, creativity, and scientific
contributions belong to the author, who maintains full responsibility
for the package's correctness and reliability. All the code has been
thoroughly tested, and users are encouraged to report any issues through
the package's issue tracker.

## Palette categories

- RElab:

  Custom palettes from the RElab research group (6)

- scientific:

  Perceptually uniform palettes, viridis, Crameri (17)

- colorblind:

  Colorblind-safe palettes, Okabe-Ito, Paul Tol, Wong (34)

- alhambra:

  Authentic Nasrid dynasty historical pigments (30)

- chameleons:

  Chameleon species coloration patterns (22)

- natural:

  Natural phenomena, ocean, forest, sunset, aurora (41)

- cultural:

  Cultural traditions, Japanese, Persian, Chinese (23)

- artistic:

  Artistic movements, impressionist, art deco, bauhaus (13)

- seasonal:

  Seasonal color schemes (8)

- modern:

  Contemporary design trends, neon, vaporwave, synthwave (12)

- classic:

  Classic color theory, primary, complementary, analogous (16)

- monochrome:

  Grayscale and single-hue scales (11)

- food:

  Food-inspired palettes (11)

- diverging:

  Diverging color scales for bipolar data (10)

- qualitative:

  Qualitative palettes for categorical data (10)

- cinema:

  Film-inspired palettes from Dune cinematography (18)

## Main functions

- [`palettes()`](https://robustecologies.github.io/koloRo/reference/palettes.md):

  Retrieve color palettes by category or name

- [`list_palettes()`](https://robustecologies.github.io/koloRo/reference/list_palettes.md):

  List all available palettes

- [`display_palette()`](https://robustecologies.github.io/koloRo/reference/display_palette.md):

  Visualize a single palette

- [`plot_palette()`](https://robustecologies.github.io/koloRo/reference/plot_palette.md):

  Flexible palette visualization

- [`compare_palettes()`](https://robustecologies.github.io/koloRo/reference/compare_palettes.md):

  Side-by-side palette comparison

- [`palette_ramp()`](https://robustecologies.github.io/koloRo/reference/palette_ramp.md):

  Create continuous color ramps

## Color manipulation

- [`adjust_brightness()`](https://robustecologies.github.io/koloRo/reference/adjust_brightness.md):

  Lighten or darken colors

- [`add_alpha()`](https://robustecologies.github.io/koloRo/reference/add_alpha.md):

  Add transparency to colors

- [`shadow_color()`](https://robustecologies.github.io/koloRo/reference/shadow_color.md):

  Create shadow effects

- [`highlight_color()`](https://robustecologies.github.io/koloRo/reference/highlight_color.md):

  Create highlight effects

- [`pattern_colors()`](https://robustecologies.github.io/koloRo/reference/pattern_colors.md):

  Generate pattern coloring schemes

## ggplot2 integration

- [`scale_color_koloro()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro.md):

  Discrete color scale for ggplot2

- [`scale_fill_koloro()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro.md):

  Discrete fill scale for ggplot2

- [`scale_color_koloro_c()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_c.md):

  Continuous color scale for ggplot2

- [`scale_fill_koloro_c()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro_c.md):

  Continuous fill scale for ggplot2

## See also

Useful links:

- <https://github.com/robustecologies/koloRo>

- <https://robustecologies.github.io/koloRo>

- Report bugs at <https://github.com/robustecologies/koloRo/issues>

## Author

**Maintainer**: Pablo Almaraz <pablo.almaraz@csic.es>
([ORCID](https://orcid.org/0000-0003-1416-2695))

Authors:

- Pablo Almaraz <pablo.almaraz@csic.es>
  ([ORCID](https://orcid.org/0000-0003-1416-2695))
