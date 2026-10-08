# Robust Ecologies Lab website

This repository is the website of the Robust Ecologies Lab. GitHub Pages publishes it from the root of the branch `main` at https://robustecologies.github.io. The pages are plain HTML, CSS and JavaScript.

The studio RElabSite writes the pages from the content in this folder. Edit the content with RElabSite and not the `index.html` files, because RElabSite writes them again at each save.

<br>

## Content

| File | Content |
|---|---|
| `site.json` | Title, menu, header icons, footer, home page portraits, comments, statistics |
| `<page>/page.json` | The fields of one page |
| `<page>/body.md` or `body.html` | The long text of a page |
| `404.json` | The page for a missing address |
| `<page>/social.png` | The image that social media show with a link to the page |
| `assets/` | Styles, scripts and fonts, written by RElabSite |
| `projects/` | The SaniVult and GuadalShiftR project sites |
| `LICENSE`, `LICENSE-CC-BY-4.0`, `LICENSES/`, `REUSE.toml` | The licences and the licence of each file (see Licences) |
| `PRIVATE/` | Private material, such as the legacy Quarto site. Git ignores it (`.gitignore`), so it is never published. RElabSite does not read it. |

The repository holds only the content of this folder. To keep another site under the same address, such as the site of an R package, add its folder in RElabSite, Site settings, Publishing.

<br>

## Publishing

RElabSite commits each change to the branch `main` and pushes it, and GitHub Pages then updates the site. No pull request or merge is necessary. The README of RElabSite gives all the steps, the folders to keep and the solutions to common problems.

<br>

## Licences

The repository holds two kinds of work under two licences. Text, figures and images are © Pablo Almaraz under the [Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0, text in `LICENSE-CC-BY-4.0`). Code is under the [GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.html) (GPL-3.0-or-later, text in `LICENSE`): the scripts and style sheets in `assets/`, the live-figure scenes in the `page.json` files and the code shown in posts. The scripts in `assets/js/` are readable source code, and they include the engine of the live figures.

| Material | Licence |
|---|---|
| Text, figures, images and data | CC BY 4.0 |
| Scripts, style sheets, live-figure scenes, code in posts | GPL-3.0-or-later |
| Figures reproduced from published articles | The licence of the article: CC BY 4.0, CC BY-NC-ND 4.0, or the terms of the publisher |
| Fonts in `assets/fonts/` | SIL Open Font License 1.1 (Jost, Libre Franklin) and GUST Font License (TeX Gyre Pagella) |
| `projects/sanivult/` | MIT |
| `projects/guadalshiftr/` | CC BY-NC-SA 4.0 |
| Libraries in `projects/*/site_libs/` | The licence stated in each library |

The repository follows the [REUSE specification](https://reuse.software) of the Free Software Foundation Europe. `REUSE.toml` gives the copyright holder and licence of every file, and `LICENSES/` holds the text of each licence. Logos and trademarks are not licensed. To check the declaration after you add or replace a file, run this command in the root of the repository:

```bash
uvx reuse lint
```
