# Changelog

## koloRo 0.1.2

### Changes

- `ggplot2` promoted from `Suggests` to `Imports`. The 12 exported
  `scale_*_koloro*` functions are no longer guarded by
  [`requireNamespace("ggplot2")`](https://ggplot2.tidyverse.org); the
  dependency is now loaded with the package. No behaviour change for
  users who already had `ggplot2` installed.

- Internal helpers `get_colored_symbol()`, `require_package()` and
  `require_packages()` removed from the export surface (they remain
  available internally via `koloRo:::`). They are package-internal
  utilities and were never intended as part of the public API.

  

## koloRo 0.1.1

### New features

- **Cinema palettes (18)**: Film-inspired color palettes based on the
  cinematography of Denis Villeneuve’s *Dune* (2021) and *Dune: Part
  Two* (2024), with color grading by Greig Fraser and David Cole.
  Includes Arrakis desert landscapes at different times of day
  (`dune_arrakis`, `dune_arrakis_dawn`, `dune_arrakis_night`,
  `dune_deep_desert`), political faction palettes (`dune_atreides`,
  `dune_harkonnen`, `dune_corrino`, `dune_fremen`), key narrative
  elements (`dune_spice`, `dune_sandworm`, `dune_caladan`), Part
  Two-specific scenes (`dune_giedi_prime`, `dune_arena`,
  `dune_southern_desert`, `dune_ritual`, `dune_war`), and overall
  cinematographic color grading (`dune_villeneuve`, `dune_part_two`).

### Bug fixes

- Fixed missing `"chameleons"` category in test suite category vectors.

------------------------------------------------------------------------

## koloRo 0.1.0

### New features

#### Palette collection (264 palettes across 15 categories)

- **RElab palettes**: Custom color schemes from the Robust Ecologies
  Lab, including `RElab_main`, `RElab_primary`, `RElab_extended`,
  `RElab_diverging`, `RElab_sequential`, and `RElab_qualitative`.

- **Scientific visualization palettes**: Complete viridis family
  (`viridis`, `magma`, `plasma`, `inferno`, `cividis`, `turbo`), Fabio
  Crameri’s scientific colormaps (`batlow`, `vik`, `cork`, `broc`), and
  perceptually uniform palettes (`twilight`, `parula`).

- **Colorblind-safe palettes (40+)**: Gold standard palettes including
  Okabe-Ito, Wong (Nature Methods), complete Paul Tol collection
  (`tol_bright`, `tol_vibrant`, `tol_muted`, `tol_light`, `tol_medium`,
  `tol_dark`), IBM colorblind palette, Tableau colorblind, and
  CVD-specific optimized palettes.

- **Alhambra historical palettes (30)**: Authentic color schemes derived
  from spectroscopic analysis of Nasrid dynasty (1232-1492) pigments,
  organized by location (Patio de los Leones, Sala de los Abencerrajes,
  Sala de las Dos Hermanas, Generalife), pigment type (blues, greens,
  reds, golds, whites, darks), and historical period (early, peak, late
  Nasrid).

- **Chameleon species palettes (22)**: Color schemes based on documented
  coloration patterns of chameleon species from scientific literature
  and natural history documentation. Includes African *Chamaeleo*
  species (`chameleon_chamaeleon`, `chameleon_calyptratus`,
  `chameleon_dilepis`, `chameleon_namaquensis`), Malagasy *Furcifer*
  species with panther chameleon locales (`chameleon_pardalis_ambilobe`,
  `chameleon_pardalis_nosybe`, `chameleon_pardalis_tamatave`), carpet
  and jeweled chameleons (`chameleon_lateralis`, `chameleon_campani`),
  horned *Trioceros* species (`chameleon_jacksonii`,
  `chameleon_melleri`, `chameleon_quadricornis`), giant *Calumma*
  species (`chameleon_parsonii`), and pygmy leaf chameleons
  (`chameleon_brookesia`).

- **Natural phenomena palettes (41)**: Oceanic (`ocean`, `deep_sea`,
  `coral_reef`, `bioluminescent`), atmospheric (`sunset`, `aurora`,
  `northern_lights`), terrestrial (`forest`, `desert`, `canyon`), and
  celestial (`galaxy`, `nebula`, `cosmic`) color schemes.

- **Cultural palettes (20+)**: Japanese (`edo`, `zen`, `sakura`),
  Persian (`miniature`, `carpet`, `garden`), Chinese (`imperial`, `ink`,
  `jade`, `porcelain`), Indian (`holi`, `spice`, `monsoon`), African,
  Nordic, and Mediterranean traditions.

- **Artistic movement palettes (15+)**: Historical movements
  (`impressionist`, `renaissance`, `baroque`, `rococo`) and modern
  styles (`art_deco`, `bauhaus`, `pop_art`, `street_art`, `minimalist`).

- **Additional categories**: Seasonal, modern design (`neon`,
  `vaporwave`, `synthwave`, `cyberpunk`), classic color theory,
  monochrome scales, food-inspired, diverging, and qualitative palettes.

#### Core functions

- [`palettes()`](https://robustecologies.github.io/koloRo/reference/palettes.md):
  Main function to retrieve palettes by category or specific name.
- [`list_palettes()`](https://robustecologies.github.io/koloRo/reference/list_palettes.md):
  List all available palettes with category and color count information.
- [`display_palette()`](https://robustecologies.github.io/koloRo/reference/display_palette.md):
  Quick visualization of a single palette.
- [`plot_palette()`](https://robustecologies.github.io/koloRo/reference/plot_palette.md):
  Flexible palette visualization supporting single palettes,
  comparisons, and hex vectors.
- [`compare_palettes()`](https://robustecologies.github.io/koloRo/reference/compare_palettes.md):
  Side-by-side comparison of multiple palettes.
- [`plot_palette_grid()`](https://robustecologies.github.io/koloRo/reference/plot_palette_grid.md):
  Grid layout visualization for multiple palettes.
- [`plot_category()`](https://robustecologies.github.io/koloRo/reference/plot_category.md):
  Display all palettes within a category.
- [`palette_ramp()`](https://robustecologies.github.io/koloRo/reference/palette_ramp.md):
  Create smooth color gradients from any palette.
- [`pattern_colors()`](https://robustecologies.github.io/koloRo/reference/pattern_colors.md):
  Generate color schemes for geometric patterns with cyclic, random,
  gradient, or alternating methods.
- [`palette_info()`](https://robustecologies.github.io/koloRo/reference/palette_info.md):
  Get detailed information about a palette including RGB values.
- [`export_palette()`](https://robustecologies.github.io/koloRo/reference/export_palette.md):
  Export palettes to hex, RGB, or normalized RGB formats.
- [`interpolate_palette()`](https://robustecologies.github.io/koloRo/reference/interpolate_palette.md):
  Smooth interpolation of palettes to any number of colors.

#### Color manipulation functions

- [`adjust_brightness()`](https://robustecologies.github.io/koloRo/reference/adjust_brightness.md):
  Lighten or darken colors.
- [`add_alpha()`](https://robustecologies.github.io/koloRo/reference/add_alpha.md):
  Add transparency to colors.
- [`shadow_color()`](https://robustecologies.github.io/koloRo/reference/shadow_color.md):
  Create shadow effects for 3D illusions.
- [`highlight_color()`](https://robustecologies.github.io/koloRo/reference/highlight_color.md):
  Create highlight effects.
- [`mix_colors()`](https://robustecologies.github.io/koloRo/reference/mix_colors.md):
  Blend two colors at specified ratio.
- [`complementary_color()`](https://robustecologies.github.io/koloRo/reference/complementary_color.md):
  Get the complementary color (opposite on color wheel).
- [`analogous_colors()`](https://robustecologies.github.io/koloRo/reference/analogous_colors.md):
  Generate analogous color schemes.
- [`triadic_colors()`](https://robustecologies.github.io/koloRo/reference/triadic_colors.md):
  Generate triadic color schemes.
- [`desaturate()`](https://robustecologies.github.io/koloRo/reference/desaturate.md):
  Reduce color saturation.
- [`invert_color()`](https://robustecologies.github.io/koloRo/reference/invert_color.md):
  Invert colors.

#### ggplot2 integration

- [`scale_color_koloro()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro.md)
  /
  [`scale_colour_koloro()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro.md):
  Discrete color scales.
- [`scale_fill_koloro()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro.md):
  Discrete fill scales.
- [`scale_color_koloro_c()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_c.md)
  /
  [`scale_colour_koloro_c()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_c.md):
  Continuous color gradients.
- [`scale_fill_koloro_c()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro_c.md):
  Continuous fill gradients.
- [`scale_color_koloro_div()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_div.md)
  /
  [`scale_colour_koloro_div()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_div.md):
  Diverging color scales with customizable midpoint.
- [`scale_fill_koloro_div()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro_div.md):
  Diverging fill scales.
- [`scale_color_koloro_binned()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_binned.md)
  /
  [`scale_colour_koloro_binned()`](https://robustecologies.github.io/koloRo/reference/scale_color_koloro_binned.md):
  Binned (stepped) color scales.
- [`scale_fill_koloro_binned()`](https://robustecologies.github.io/koloRo/reference/scale_fill_koloro_binned.md):
  Binned fill scales.

#### Utility functions

- `get_colored_symbol()`: ANSI-colored Unicode symbols for terminal
  output.
- `require_package()` / `require_packages()`: Interactive package
  dependency management.
- [`is_color()`](https://robustecologies.github.io/koloRo/reference/is_color.md):
  Validate color specifications.
- [`hex_to_rgb()`](https://robustecologies.github.io/koloRo/reference/hex_to_rgb.md)
  /
  [`rgb_to_hex()`](https://robustecologies.github.io/koloRo/reference/rgb_to_hex.md):
  Color format conversion.

#### Documentation

- Six comprehensive vignettes:
  - “Introduction to koloRo” - Quick start guide and overview
  - “Colorblind-safe palettes” - Accessibility best practices
  - “Alhambra palettes” - Historical context and Nasrid dynasty pigments
  - “Chameleon palettes” - Species-based palettes with scientific
    references
  - “Scientific palettes” - Perceptually uniform colormaps for
    oceanography
  - “ggplot2 integration” - Complete guide to ggplot2 scales
