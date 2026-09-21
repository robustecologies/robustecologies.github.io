# DynFlow

DynFlow builds live figures of dynamical systems for web pages, slides and packages. It has two parts. In the studio (`index.html`), you choose or write a model, perturb it, choose a view and a style, and export the result. The player (`dist/dynflow.js`) is one script file that plays a saved scene anywhere, as a `<dyn-flow>` element or as a standalone page.

DynFlow is plain JavaScript. It has no build step for use, no server requirement and no dependency on R or on any package. It includes a catalogue of 82 models from ecology, chaotic flows, forced oscillators, maps, delay equations, stochastic processes, tipping and epidemics.

<br>

## Start the studio

1. Open `index.html` in a recent browser (Chrome, Edge, Firefox or Safari).
2. Choose a model in the library on the left, or choose a scene in the showcase.
3. Change the parameters, the perturbations, the view and the style in the inspector on the right.
4. Export the figure with the Export button.

The studio also works when you open the file from disk. To use share links and the clipboard, serve the folder:

```bash
cd /path/to/DynFlowApp
python3 -m http.server 8000
# then open http://127.0.0.1:8000
```

<br>

## Keys

| Key | Action |
|---|---|
| Space | Play or pause |
| R | Restart with the same seed |
| N | Restart with a new seed |
| S | Advance one frame |
| F | Full screen |
| P | Save a PNG |

In a 2D view, a click on the figure starts new orbits. In a 3D view, a drag rotates the view.

<br>

## Views

| View | Shows |
|---|---|
| Flow | An ensemble of particles with fading trails, in 2D, rotating 3D or on the simplex of three shares |
| Trajectory | A few long orbits with a colour gradient along the tail |
| Time series | Chosen variables against time |
| Phase plane | Vector field, nullclines, equilibria with their stability, and orbits |
| Sweep | A slow sweep of one parameter over the branches of equilibria (hysteresis) |
| Bifurcation | Iterates or local maxima against one parameter |
| Density | Ensemble density as a heat map, or one variable over time |
| Strobe | The state every period T, or its crossings of a Poincare plane |
| Cobweb | Graphical iteration of a one-dimensional map |

<br>

## Perturbations

A parameter can receive periodic forcing, quasiperiodic forcing, a ramp, a step or coloured noise. A state can receive additive, multiplicative, coloured or alpha-stable noise, random jumps or periodic pulses. Mark a state perturbation as common to give the whole ensemble one realisation of the noise.

<br>

## Write a model

Select Write a model, then type one equation per line:

```
x' = a*x - b*x*y              # differential equation
y[n+1] = r*y*(1 - y)          # map (do not mix with x')
noise x = sigma*x             # Ito diffusion of x
z' = -z + lag(x, tau)         # delay: x(t - tau)
aux f = x/(1 + x)             # helper
param a = 1 [0, 3]            # value and slider range
init x = 0.5                  # initial state
range x = [0, 4]              # axis range
```

Press Apply, or Ctrl + Enter. The Help dialog lists the functions. `THEORY.md` states the numerical methods.

<br>

## Use a scene elsewhere

Save the scene as JSON (Export, Scene file). Then copy `dist/dynflow.js` next to your page and add:

```html
<script src="dynflow.js"></script>
<dyn-flow src="my-scene.json" controls style="height:420px"></dyn-flow>
```

You can place a catalogue model without a scene file:

```html
<dyn-flow model="lorenz" theme="blackboard" controls></dyn-flow>
```

The attributes of the element are `scene`, `src`, `model`, `theme`, `view`, `title`, `controls`, `paused` and `static`. The element starts when it becomes visible. When the reader prefers reduced motion, the element shows a still frame.

In R Markdown, Quarto or pkgdown, emit the same two tags:

```r
htmltools::tagList(
  htmltools::tags$script(src = "dynflow.js"),
  htmltools::HTML('<dyn-flow src="my-scene.json" controls style="display:block;height:420px"></dyn-flow>')
)
```

For slides, use Export, Standalone HTML page. The page contains the player and the scene, so it runs offline and in an iframe. The file `examples/embed.html` shows every way to embed a scene.

<br>

## Exports

| Export | Result |
|---|---|
| PNG | The current frame at the current resolution. Set Style, Resolution to 2x or 3x for print. |
| SVG | Orbits, nullclines and branches as vector paths; particle clouds as an embedded image |
| WebM | A video of the running figure. Convert to MP4 with `ffmpeg -i scene.webm scene.mp4`. |
| GIF | An animation for slides and messages |
| Standalone HTML | One file with the player and the scene |
| Embed code | The script and element tags for a web page |
| Scene file | JSON with the equations, parameters, perturbations, view and style |
| Share link | A URL that holds the compressed scene |

<br>

## Build and test

Edit the files in `src/`, then rebuild the player:

```bash
node tools/build.mjs
```

Run the tests:

```bash
node tests/core.test.mjs       # schemes, formulas, random numbers, perturbations
node tests/analysis.test.mjs   # eigenvalues, equilibria, bifurcation points
node tests/models.test.mjs     # every catalogue model and the claims in its description
node tests/browser.test.mjs    # studio, exports and element in headless Chrome
```

The browser test needs Google Chrome, Python 3, Pillow and ffprobe.

<br>

## Files

| Path | Content |
|---|---|
| `index.html`, `src/studio/` | The studio |
| `src/core/` | Formula compiler, random numbers, simulator, local analysis |
| `src/render/` | Themes and palettes, projections, views, player, exports |
| `src/models/catalogue.js` | The 82 models, each with its source |
| `src/component/dyn-flow.js` | The `<dyn-flow>` element |
| `dist/dynflow.js` | The player in one file, built by `tools/build.mjs` |
| `examples/` | Embedding examples and scene files |
| `tests/` | Numerical and browser tests |
| `vendor/` | KaTeX and fonts |

<br>

## Licence

DynFlow is free software under the GNU General Public License, version 3 or later (`LICENSE`). Copyright (C) 2026 Pablo Almaraz, Robust Ecologies Lab (ORCID [0000-0003-1416-2695](https://orcid.org/0000-0003-1416-2695)).

Third-party files in `vendor/` keep their own licences: KaTeX (MIT, `vendor/katex/LICENSE`), Jost and Libre Franklin (SIL Open Font License 1.1, `vendor/fonts/OFL-1.1.txt`), and TeX Gyre Pagella (GUST Font License, `vendor/fonts/GUST-FONT-LICENSE.txt`).
