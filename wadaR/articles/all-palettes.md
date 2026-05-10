# Complete palette reference

``` r

library(koloRo)
```

This vignette provides a complete visual reference of all **282
palettes** available in koloRo, organized by category.

## Overview

``` r

lp <- list_palettes()
cat("Total palettes:", nrow(lp), "\n\n")
#> Total palettes: 283
cat("Palettes by category:\n")
#> Palettes by category:
sort(table(lp$category), decreasing = TRUE)
#> 
#>     natural  colorblind    alhambra    cultural  chameleons      cinema 
#>          41          34          30          23          22          18 
#>  scientific     classic    artistic      modern        food  monochrome 
#>          18          16          13          12          11          11 
#>   diverging qualitative    seasonal       RElab 
#>          10          10           8           6
```

## RElab palettes (6)

Custom palettes from the Robust Ecologies Laboratory, designed for
scientific publications.

``` r

relab_pals <- list_palettes()$palette_name[list_palettes()$category == "RElab"]
plot_palette(relab_pals)
```

![](all-palettes_files/figure-html/relab-1.png)

## Scientific palettes (17)

Perceptually uniform colormaps for scientific visualization, including
viridis family ([van der Walt & Smith (2015)](#ref-vanderwalt2015)) and
Fabio Crameri’s scientific colormaps ([Crameri
(2018)](#ref-crameri2018)).

``` r

sci_pals <- list_palettes()$palette_name[list_palettes()$category == "scientific"]
plot_palette(sci_pals)
```

![](all-palettes_files/figure-html/scientific-1.png)

## Colorblind-safe palettes (34)

Palettes designed for accessibility, validated for deuteranopia,
protanopia, and tritanopia. Includes Okabe-Ito ([Okabe & Ito
(2008)](#ref-okabe2008)), Paul Tol ([Tol (2021)](#ref-tol2021)), Bang
Wong ([Wong (2011)](#ref-wong2011)), and other accessibility-focused
palettes.

``` r

cb_pals <- list_palettes()$palette_name[list_palettes()$category == "colorblind"]
n_cb <- length(cb_pals)
plot_palette(cb_pals[1:17])
```

![](all-palettes_files/figure-html/colorblind-1-1.png)

``` r

plot_palette(cb_pals[18:n_cb])
```

![](all-palettes_files/figure-html/colorblind-2-1.png)

## Alhambra palettes (30)

Authentic historical pigments from the 13th-14th century Nasrid dynasty
palace ([Fernández-Puertas (1997)](#ref-fernandez1997)), derived from
spectroscopic analysis of original decorations ([García-Bueno &
Medina-Flórez (2004)](#ref-garcia2004); [Cardell & Navarrete
(2006)](#ref-cardell2006)).

``` r

alh_pals <- list_palettes()$palette_name[list_palettes()$category == "alhambra"]
n_alh <- length(alh_pals)
plot_palette(alh_pals[1:15])
```

![](all-palettes_files/figure-html/alhambra-1-1.png)

``` r

plot_palette(alh_pals[16:n_alh])
```

![](all-palettes_files/figure-html/alhambra-2-1.png)

## Chameleons palettes (22)

Color palettes derived from scientifically documented coloration
patterns of chameleon species worldwide. The collection spans African
*Chamaeleo* species, Malagasy *Furcifer* and *Calumma* genera including
multiple panther chameleon (*Furcifer pardalis*) locales, East African
horned *Trioceros* species, and pygmy *Brookesia* leaf chameleons. All
colors are based on herpetological literature and natural history
documentation ([Stuart-Fox & Moussalli (2008)](#ref-stuart2008);
[Teyssier et al. (2015)](#ref-teyssier2015)).

``` r

cham_pals <- list_palettes()$palette_name[list_palettes()$category == "chameleons"]
n_cham <- length(cham_pals)
plot_palette(cham_pals[1:11])
```

![](all-palettes_files/figure-html/chameleons-1-1.png)

``` r

plot_palette(cham_pals[12:n_cham])
```

![](all-palettes_files/figure-html/chameleons-2-1.png)

## Natural palettes (41)

Colors inspired by natural phenomena: oceans, forests, sunsets, aurora
borealis, deserts, and more.

``` r

nat_pals <- list_palettes()$palette_name[list_palettes()$category == "natural"]
n_nat <- length(nat_pals)
plot_palette(nat_pals[1:14])
```

![](all-palettes_files/figure-html/natural-1-1.png)

``` r

plot_palette(nat_pals[15:28])
```

![](all-palettes_files/figure-html/natural-2-1.png)

``` r

plot_palette(nat_pals[29:n_nat])
```

![](all-palettes_files/figure-html/natural-3-1.png)

## Cultural palettes (23)

Traditional color schemes from world cultures: Japanese, Persian,
Chinese, Indian, Mediterranean, and more.

``` r

cult_pals <- list_palettes()$palette_name[list_palettes()$category == "cultural"]
n_cult <- length(cult_pals)
plot_palette(cult_pals[1:12])
```

![](all-palettes_files/figure-html/cultural-1-1.png)

``` r

plot_palette(cult_pals[13:n_cult])
```

![](all-palettes_files/figure-html/cultural-2-1.png)

## Artistic palettes (13)

Palettes inspired by art movements and famous artists: Impressionism,
Bauhaus, Art Deco, Pop Art, and more.

``` r

art_pals <- list_palettes()$palette_name[list_palettes()$category == "artistic"]
plot_palette(art_pals)
```

![](all-palettes_files/figure-html/artistic-1.png)

## Classic palettes (16)

Timeless color combinations that work across many contexts.

``` r

class_pals <- list_palettes()$palette_name[list_palettes()$category == "classic"]
plot_palette(class_pals)
```

![](all-palettes_files/figure-html/classic-1.png)

## Modern palettes (12)

Contemporary color schemes reflecting current design trends.

``` r

mod_pals <- list_palettes()$palette_name[list_palettes()$category == "modern"]
plot_palette(mod_pals)
```

![](all-palettes_files/figure-html/modern-1.png)

## Food palettes (11)

Appetizing colors inspired by culinary themes.

``` r

food_pals <- list_palettes()$palette_name[list_palettes()$category == "food"]
plot_palette(food_pals)
```

![](all-palettes_files/figure-html/food-1.png)

## Monochrome palettes (11)

Single-hue gradients for elegant, focused visualizations.

``` r

mono_pals <- list_palettes()$palette_name[list_palettes()$category == "monochrome"]
plot_palette(mono_pals)
```

![](all-palettes_files/figure-html/monochrome-1.png)

## Diverging palettes (10)

Bipolar color schemes with neutral midpoints, ideal for data with
meaningful center values.

``` r

div_pals <- list_palettes()$palette_name[list_palettes()$category == "diverging"]
plot_palette(div_pals)
```

![](all-palettes_files/figure-html/diverging-1.png)

## Qualitative palettes (10)

Distinct colors for categorical data where no inherent ordering exists.

``` r

qual_pals <- list_palettes()$palette_name[list_palettes()$category == "qualitative"]
plot_palette(qual_pals)
```

![](all-palettes_files/figure-html/qualitative-1.png)

## Seasonal palettes (8)

Colors evoking the four seasons and holidays.

``` r

seas_pals <- list_palettes()$palette_name[list_palettes()$category == "seasonal"]
plot_palette(seas_pals)
```

![](all-palettes_files/figure-html/seasonal-1.png)

## Cinema palettes (18)

Color palettes inspired by the cinematography of Denis Villeneuve’s
*Dune* (2021) and *Dune: Part Two* (2024), reflecting the color grading
of Greig Fraser and David Cole. The collection spans Arrakis desert
landscapes at different times of day, political faction color
identities, key narrative elements such as spice and sandworms, Part
Two-specific scenes including the Giedi Prime infrared sequence, and
overall cinematographic color grading signatures.

``` r

cin_pals <- list_palettes()$palette_name[list_palettes()$category == "cinema"]
n_cin <- length(cin_pals)
plot_palette(cin_pals[1:9])
```

![](all-palettes_files/figure-html/cinema-1-1.png)

``` r

plot_palette(cin_pals[10:n_cin])
```

![](all-palettes_files/figure-html/cinema-2-1.png)

## Quick reference table

``` r

lp <- list_palettes()
knitr::kable(
  lp[order(lp$category, lp$palette_name), ],
  row.names = FALSE,
  col.names = c("Palette", "Category", "Colors"),
  caption = "All 282 palettes in koloRo"
)
```

| Palette                      | Category    | Colors |
|:-----------------------------|:------------|-------:|
| alhambra_alicatado           | alhambra    |      6 |
| alhambra_artesonado          | alhambra    |      6 |
| alhambra_blues               | alhambra    |      4 |
| alhambra_cool                | alhambra    |      6 |
| alhambra_darks               | alhambra    |      4 |
| alhambra_dos_hermanas        | alhambra    |      5 |
| alhambra_early_nasrid        | alhambra    |      5 |
| alhambra_earth               | alhambra    |      5 |
| alhambra_embajadores         | alhambra    |      6 |
| alhambra_generalife          | alhambra    |      6 |
| alhambra_geometric           | alhambra    |      5 |
| alhambra_golds               | alhambra    |      3 |
| alhambra_greens              | alhambra    |      4 |
| alhambra_late_nasrid         | alhambra    |      5 |
| alhambra_minimal_blue_white  | alhambra    |      2 |
| alhambra_minimal_gold_blue   | alhambra    |      2 |
| alhambra_minimal_green_white | alhambra    |      2 |
| alhambra_mudejar             | alhambra    |      6 |
| alhambra_nazari              | alhambra    |      5 |
| alhambra_patio_leones        | alhambra    |      5 |
| alhambra_peak_nasrid         | alhambra    |      6 |
| alhambra_reds                | alhambra    |      4 |
| alhambra_restoration         | alhambra    |      6 |
| alhambra_sala_abencerrajes   | alhambra    |      5 |
| alhambra_spectrum            | alhambra    |      8 |
| alhambra_sunset              | alhambra    |      6 |
| alhambra_warm                | alhambra    |      6 |
| alhambra_whites              | alhambra    |      3 |
| alhambra_yeseria             | alhambra    |      5 |
| alhambra_zellige             | alhambra    |      6 |
| art_deco                     | artistic    |     10 |
| art_nouveau                  | artistic    |     10 |
| baroque                      | artistic    |     10 |
| bauhaus                      | artistic    |     10 |
| cubist                       | artistic    |     10 |
| impressionist                | artistic    |     10 |
| minimalist                   | artistic    |     10 |
| pop_art                      | artistic    |     10 |
| post_impressionist           | artistic    |     10 |
| renaissance                  | artistic    |     10 |
| rococo                       | artistic    |     10 |
| street_art                   | artistic    |     10 |
| surrealist                   | artistic    |     10 |
| chameleon_brookesia          | chameleons  |      8 |
| chameleon_calyptratus        | chameleons  |      8 |
| chameleon_campani            | chameleons  |      8 |
| chameleon_chamaeleon         | chameleons  |      7 |
| chameleon_dilepis            | chameleons  |      7 |
| chameleon_fischeri           | chameleons  |      8 |
| chameleon_jacksonii          | chameleons  |      8 |
| chameleon_labordi            | chameleons  |      8 |
| chameleon_lateralis          | chameleons  |      8 |
| chameleon_melleri            | chameleons  |      8 |
| chameleon_minor              | chameleons  |      8 |
| chameleon_namaquensis        | chameleons  |      8 |
| chameleon_nasutum            | chameleons  |      8 |
| chameleon_oustaleti          | chameleons  |      8 |
| chameleon_pardalis_ambilobe  | chameleons  |      8 |
| chameleon_pardalis_nosybe    | chameleons  |      8 |
| chameleon_pardalis_tamatave  | chameleons  |      8 |
| chameleon_parsonii           | chameleons  |      8 |
| chameleon_quadricornis       | chameleons  |      8 |
| chameleon_senegalensis       | chameleons  |      7 |
| chameleon_verrucosus         | chameleons  |      8 |
| chameleon_werneri            | chameleons  |      8 |
| dune_arena                   | cinema      |      8 |
| dune_arrakis                 | cinema      |      8 |
| dune_arrakis_dawn            | cinema      |      8 |
| dune_arrakis_night           | cinema      |      8 |
| dune_atreides                | cinema      |      8 |
| dune_caladan                 | cinema      |      8 |
| dune_corrino                 | cinema      |      8 |
| dune_deep_desert             | cinema      |      8 |
| dune_fremen                  | cinema      |      8 |
| dune_giedi_prime             | cinema      |      8 |
| dune_harkonnen               | cinema      |      8 |
| dune_part_two                | cinema      |      8 |
| dune_ritual                  | cinema      |      8 |
| dune_sandworm                | cinema      |      8 |
| dune_southern_desert         | cinema      |      8 |
| dune_spice                   | cinema      |      8 |
| dune_villeneuve              | cinema      |      8 |
| dune_war                     | cinema      |      8 |
| analogous_cool               | classic     |      3 |
| analogous_warm               | classic     |      3 |
| complementary_rb             | classic     |      2 |
| complementary_yp             | classic     |      2 |
| munsell_blue                 | classic     |     10 |
| munsell_green                | classic     |     10 |
| munsell_purple               | classic     |     10 |
| munsell_red                  | classic     |     10 |
| munsell_yellow               | classic     |     10 |
| primary                      | classic     |      3 |
| secondary                    | classic     |      3 |
| split_complementary          | classic     |      3 |
| square                       | classic     |      4 |
| tertiary                     | classic     |      9 |
| tetradic                     | classic     |      4 |
| triadic                      | classic     |      3 |
| cb_dark2                     | colorblind  |      8 |
| cb_paired                    | colorblind  |     10 |
| cb_set2                      | colorblind  |      8 |
| cud_blue_red                 | colorblind  |      2 |
| cud_bluish_green_red         | colorblind  |      2 |
| cud_sky_orange               | colorblind  |      2 |
| cud_vermillion_blue          | colorblind  |      2 |
| deuteranopia_safe            | colorblind  |      8 |
| ibm_colorblind               | colorblind  |      5 |
| ibm_colorblind_extended      | colorblind  |      7 |
| okabe_ito                    | colorblind  |      8 |
| okabe_ito_extended           | colorblind  |     10 |
| protanopia_safe              | colorblind  |      8 |
| r_colorblind_safe            | colorblind  |      8 |
| safe_diverging_brown_teal    | colorblind  |      8 |
| safe_diverging_pink_green    | colorblind  |      8 |
| safe_diverging_purple_green  | colorblind  |      8 |
| safe_qualitative_4           | colorblind  |      4 |
| safe_qualitative_5           | colorblind  |      5 |
| safe_qualitative_6           | colorblind  |      6 |
| safe_qualitative_8           | colorblind  |      8 |
| seaborn_colorblind           | colorblind  |     10 |
| tableau_colorblind10         | colorblind  |     10 |
| tol_bright                   | colorblind  |      7 |
| tol_dark                     | colorblind  |      6 |
| tol_high_contrast            | colorblind  |      3 |
| tol_light                    | colorblind  |      9 |
| tol_medium                   | colorblind  |      6 |
| tol_muted                    | colorblind  |      9 |
| tol_pale                     | colorblind  |      6 |
| tol_vibrant                  | colorblind  |      7 |
| traffic_safe                 | colorblind  |      3 |
| tritanopia_safe              | colorblind  |      8 |
| wong                         | colorblind  |      8 |
| african_sunset               | cultural    |     10 |
| african_textile              | cultural    |     10 |
| african_wildlife             | cultural    |     10 |
| arabian_nights               | cultural    |     10 |
| chinese_imperial             | cultural    |     10 |
| chinese_ink                  | cultural    |     10 |
| chinese_jade                 | cultural    |     10 |
| chinese_porcelain            | cultural    |     10 |
| greek_islands                | cultural    |     10 |
| indian_holi                  | cultural    |     10 |
| indian_monsoon               | cultural    |     10 |
| indian_spice                 | cultural    |     10 |
| japanese_autumn              | cultural    |     10 |
| japanese_edo                 | cultural    |     10 |
| japanese_sakura              | cultural    |     10 |
| japanese_zen                 | cultural    |     10 |
| mediterranean                | cultural    |     10 |
| moroccan                     | cultural    |     10 |
| nordic_aurora                | cultural    |     10 |
| nordic_winter                | cultural    |     10 |
| persian_carpet               | cultural    |     10 |
| persian_garden               | cultural    |     10 |
| persian_miniature            | cultural    |     10 |
| bad_good                     | diverging   |     10 |
| brown_blue_green             | diverging   |     11 |
| cool_warm_diverging          | diverging   |      8 |
| negative_positive            | diverging   |     10 |
| pink_yellow_green            | diverging   |     11 |
| purple_green                 | diverging   |     11 |
| purple_orange                | diverging   |     10 |
| red_yellow_blue              | diverging   |     11 |
| spectral_diverging           | diverging   |     11 |
| vik_diverging                | diverging   |     12 |
| berry                        | food        |     10 |
| candy                        | food        |     10 |
| chocolate                    | food        |     10 |
| citrus                       | food        |     10 |
| coffee                       | food        |     10 |
| herbs                        | food        |     10 |
| ice_cream                    | food        |     10 |
| spices                       | food        |     10 |
| tea                          | food        |     10 |
| tropical                     | food        |     10 |
| wine                         | food        |     10 |
| corporate_blue               | modern      |     10 |
| corporate_gray               | modern      |     10 |
| cyberpunk                    | modern      |     10 |
| instagram_gradient           | modern      |     10 |
| japanese_minimal             | modern      |     10 |
| neon                         | modern      |     10 |
| nordic_minimal               | modern      |     10 |
| pastel                       | modern      |     10 |
| synthwave                    | modern      |     10 |
| tech_gradient                | modern      |     10 |
| tiktok                       | modern      |     10 |
| vaporwave                    | modern      |     10 |
| blue_scale                   | monochrome  |     10 |
| cool_gray                    | monochrome  |     10 |
| grayscale                    | monochrome  |     11 |
| green_scale                  | monochrome  |     10 |
| orange_scale                 | monochrome  |     10 |
| purple_scale                 | monochrome  |     10 |
| red_scale                    | monochrome  |     10 |
| sepia                        | monochrome  |     10 |
| vintage                      | monochrome  |     10 |
| warm_gray                    | monochrome  |     10 |
| yellow_scale                 | monochrome  |     10 |
| abyssal                      | natural     |     10 |
| aurora                       | natural     |     10 |
| bioluminescent               | natural     |     10 |
| boreal                       | natural     |     10 |
| canyon                       | natural     |     10 |
| coral_reef                   | natural     |     10 |
| cosmic                       | natural     |     10 |
| dawn                         | natural     |     10 |
| deep_sea                     | natural     |     10 |
| desert                       | natural     |     10 |
| dusk                         | natural     |     10 |
| earth                        | natural     |     10 |
| fire                         | natural     |     11 |
| forest                       | natural     |     10 |
| galaxy                       | natural     |     10 |
| granite                      | natural     |     10 |
| inferno_natural              | natural     |     10 |
| interstellar                 | natural     |     10 |
| kelp_forest                  | natural     |     10 |
| lava                         | natural     |     11 |
| lightning                    | natural     |     10 |
| midday                       | natural     |     10 |
| midnight                     | natural     |     10 |
| nebula                       | natural     |     10 |
| northern_lights              | natural     |     10 |
| ocean                        | natural     |     10 |
| rainforest                   | natural     |     10 |
| sandstone                    | natural     |     10 |
| savanna                      | natural     |     10 |
| stardust                     | natural     |     10 |
| stone                        | natural     |     10 |
| storm                        | natural     |     10 |
| sunrise                      | natural     |     10 |
| sunset                       | natural     |     11 |
| tidal                        | natural     |     10 |
| tornado                      | natural     |     10 |
| tropical                     | natural     |     10 |
| tundra                       | natural     |     10 |
| twilight_natural             | natural     |     10 |
| volcano                      | natural     |     10 |
| wetland                      | natural     |     10 |
| category10                   | qualitative |     10 |
| category20                   | qualitative |     20 |
| kelly_colors                 | qualitative |     22 |
| paired_12                    | qualitative |     12 |
| safe_qual_8                  | qualitative |      8 |
| set1_9                       | qualitative |      9 |
| set2_8                       | qualitative |      8 |
| set3_12                      | qualitative |     12 |
| tableau10                    | qualitative |     10 |
| tableau20                    | qualitative |     20 |
| RElab_diverging              | RElab       |      9 |
| RElab_extended               | RElab       |      9 |
| RElab_main                   | RElab       |      5 |
| RElab_primary                | RElab       |      5 |
| RElab_qualitative            | RElab       |      8 |
| RElab_sequential             | RElab       |      9 |
|                              | scientific  |     10 |
| batlow                       | scientific  |     10 |
| batlowK                      | scientific  |     10 |
| batlowW                      | scientific  |      9 |
| broc                         | scientific  |     12 |
| cividis                      | scientific  |     10 |
| coolwarm                     | scientific  |      8 |
| cork                         | scientific  |      9 |
| inferno                      | scientific  |     10 |
| magma                        | scientific  |      6 |
| parula                       | scientific  |      9 |
| plasma                       | scientific  |      6 |
| rdylbu                       | scientific  |     11 |
| spectral                     | scientific  |     11 |
| turbo                        | scientific  |      9 |
| twilight                     | scientific  |     12 |
| vik                          | scientific  |     12 |
| viridis                      | scientific  |      6 |
| autumn_foliage               | seasonal    |     10 |
| autumn_harvest               | seasonal    |     10 |
| spring_blossom               | seasonal    |     10 |
| spring_garden                | seasonal    |     10 |
| summer_beach                 | seasonal    |     10 |
| summer_tropical              | seasonal    |     10 |
| winter_frost                 | seasonal    |     10 |
| winter_holiday               | seasonal    |     10 |

All 282 palettes in koloRo {.table}

## Accessing palettes

``` r

# Get a specific palette
palettes(palette = "viridis")

# Get all palettes from a category
palettes(category = "colorblind")

# List all palettes
list_palettes()

# Visualize multiple palettes
plot_palette(c("okabe_ito", "wong", "tol_bright"))

# Get palette info
palette_info("alhambra_nazari")
```

## References

van der Walt, S., & Smith, N. (2015). A better default colormap for
Matplotlib. *Proceedings of SciPy 2015*.
<https://bids.github.io/colormap/>

Crameri, F. (2018). Scientific colour maps. Zenodo.
[doi:10.5281/zenodo.1243862](https://doi.org/10.5281/zenodo.1243862)

Okabe, M., & Ito, K. (2008). Color Universal Design (CUD). *J*Fly\*.
<https://jfly.uni-koeln.de/color/>

Tol, P. (2021). Colour schemes. SRON Technical Note.
<https://personal.sron.nl/~pault/>

Wong, B. (2011). Points of view: Color blindness. *Nature Methods*,
8(6), 441. [doi:10.1038/nmeth.1618](https://doi.org/10.1038/nmeth.1618)

Fernández-Puertas, A. (1997). *The Alhambra: From the ninth century to
Yusuf I*. Saqi Books.

García-Bueno, A., & Medina-Flórez, V.J. (2004). The Nasrid plasterwork
at “qubba Dar al-Manjara l-kubra” in Granada. *Journal of Cultural
Heritage*, 5(1), 75-89.

Cardell, C., & Navarrete, L. (2006). Pigment and plasterwork analyses of
Nasrid polychrome stucco in the Alhambra. *Studies in Conservation*,
51(3), 161-176.

Stuart-Fox, D., & Moussalli, A. (2008). Selection for social signalling
drives the evolution of chameleon colour change. *PLoS Biology*, 6(1),
e25.
[doi:10.1371/journal.pbio.0060025](https://doi.org/10.1371/journal.pbio.0060025)

Teyssier, J., Saenko, S.V., van der Marel, D., & Milinkovitch, M.C.
(2015). Photonic crystals cause active colour change in chameleons.
*Nature Communications*, 6, 6368.
[doi:10.1038/ncomms7368](https://doi.org/10.1038/ncomms7368)

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8    
#>  [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8   
#>  [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Madrid
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] koloRo_0.1.2
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.0    rstudioapi_0.18.0 tools_4.6.0       ragg_1.5.2       
#> [21] bslib_0.10.0      evaluate_1.0.5    yaml_2.3.12       otel_0.2.0       
#> [25] jsonlite_2.0.0    htmlwidgets_1.6.4 rlang_1.2.0       fs_2.1.0
```
