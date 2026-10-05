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

The repository holds only the content of this folder. To keep another site under the same address, such as the site of an R package, add its folder in RElabSite, Site settings, Publishing.

<br>

## Publishing

RElabSite commits each change to the branch `preview` and pushes it. Merge `preview` into `main` on GitHub to update the site. The README of RElabSite gives all the steps, the folders to keep and the solutions to common problems.
