# Plot a directed acyclic graph

Produces a DAG visualization of Bayesian model structure, with 7 node
types, automatic classification heuristics, plate-notation grouping
boxes, and stochastic/deterministic edge styling. Two backends are
available: DiagrammeR (default, static Graphviz rendering) and
visNetwork (interactive HTML widget with draggable nodes and tooltips).

## Usage

``` r
plot_dag(
  formula,
  groups = NULL,
  node_types = NULL,
  stochastic = NULL,
  backend = c("diagrammer", "visnetwork"),
  plates = TRUE,
  layout = c("hierarchical", "physics"),
  direction = c("LR", "UD", "DU", "RL"),
  title = NULL,
  greek = TRUE,
  legend = TRUE,
  col = NULL,
  width = NULL,
  height = NULL,
  file = NULL,
  data_nodes = NULL,
  hyper_nodes = NULL,
  interactive = NULL,
  ...
)
```

## Arguments

- formula:

  Model specification in one of four formats:

  Named list

  :   A named list of character vectors encoding parent relationships,
      e.g., `list(y = c("mu", "sigma"), mu = c("mu_0", "sigma_0"))`.

  Formula list

  :   A list of R formulas, e.g.,
      `list(y ~ mu + sigma, mu ~ mu_0 + sigma_0)`.

  Text/LaTeX notation

  :   A single character string with one sampling statement per line in
      the form `LHS ~ Dist(arg1, arg2)`. LaTeX markup is automatically
      stripped, so both plain text and copy-pasted LaTeX from vignettes
      or papers work directly.

  Model function or demonoid

  :   A lucifer `Model` function (signature `function(parm, Data)`) or a
      fitted `demonoid` object returned by
      [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
      The DAG is extracted via static analysis of the function body,
      identifying parameter extractions, distribution calls, and
      deterministic links.

- groups:

  Optional named list mapping group labels to node name vectors, e.g.,
  `list("Observation" = c("y", "sigma"), "Group" = "theta")`. When
  `NULL` (default), groups are inferred from node types.

- node_types:

  Optional named character vector mapping node names to types (one of
  `"data"`, `"parameter"`, `"latent"`, `"hyperparameter"`,
  `"process_error"`, `"obs_error"`, `"missing"`). Overrides
  auto-classification.

- stochastic:

  Controls which edges render as dashed (stochastic) arrows. `NULL`
  (default) auto-detects from node types: outgoing edges from
  `hyperparameter`, `obs_error`, and `process_error` nodes are dashed. A
  character vector overrides with explicit node names. `FALSE` forces
  all edges solid.

- backend:

  Character. Rendering backend: `"diagrammer"` (default, static Graphviz
  via DOT) or `"visnetwork"` (interactive HTML widget).

- plates:

  Logical. Whether to draw plate-notation boxes around groups (default
  `TRUE`). For DiagrammeR, uses Graphviz subgraph clusters; for
  visNetwork, uses JavaScript canvas drawing.

- layout:

  Character. Layout algorithm: `"hierarchical"` (default) or `"physics"`
  (force-directed, visNetwork only).

- direction:

  Character. Layout direction: `"LR"` (left-to-right, default), `"UD"`
  (top-down), `"DU"`, or `"RL"`.

- title:

  Optional character string for the plot title.

- greek:

  Logical. Convert parameter names to Greek letters with subscripts
  (default `TRUE`). For example, `sigma_obs` renders as the Greek sigma
  with subscript "obs". DiagrammeR uses Graphviz HTML labels; visNetwork
  uses Unicode characters. Set `FALSE` for raw names.

- legend:

  Logical. Draw a pictorial legend explaining node shapes, colors, and
  edge styles (default `TRUE`). Currently only applies to the visNetwork
  backend, which uses its built-in legend panel. Silently ignored for
  DiagrammeR. Set `FALSE` to suppress.

- col:

  Optional character vector of 7 colors overriding the default palette
  (order: data, parameter, latent, hyperparameter, process_error,
  obs_error, missing).

- width, height:

  Widget dimensions. Defaults to `NULL` (auto) for DiagrammeR, `"100%"`
  and `"600px"` for visNetwork.

- file:

  Optional file path to save the plot. Supported extensions: `.png`,
  `.pdf`, `.svg` (requires the DiagrammeRsvg and rsvg packages for
  raster/vector export), and `.html` (self-contained widget via
  htmlwidgets).

- data_nodes:

  Deprecated. Use `node_types` instead.

- hyper_nodes:

  Deprecated. Use `node_types` instead.

- interactive:

  Deprecated. Use `backend` instead.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the htmlwidget (class `grViz` for DiagrammeR, class
`visNetwork` for visNetwork).

## Details

**Input formats.** The `formula` argument accepts four complementary
input formats. Named lists and formula lists provide explicit control
over the graph structure. Text/LaTeX notation lets you copy-paste
mathematical model specifications directly from papers or vignettes,
with all LaTeX markup stripped automatically. Passing a lucifer `Model`
function or a fitted `demonoid` object triggers static analysis of the
function body: the parser identifies parameter extractions from
`parm[]`, distribution calls (`dnorm`, `dgamma`, etc.), and
deterministic assignments to reconstruct the DAG. The Model parser works
best with standard lucifer model structure (hyperparameters, priors,
likelihood sections) and may miss edges in models with complex control
flow or unusual coding patterns; for such cases, provide explicit edges
via named lists.

**Node types.** The function classifies nodes into 7 types, each with a
distinct visual encoding designed around three shape families. The *data
family* (rectangles) includes `data` nodes as filled rounded rectangles
in deep navy and `missing` nodes as dashed-border rectangles in red. The
*variable family* (circles) includes `parameter` nodes as filled blue
circles, `latent` nodes as dashed-border teal circles, and
`hyperparameter` nodes as smaller purple circles. The *distribution
family* (ovals/ellipses) includes `process_error` as green ellipses and
`obs_error` as orange ellipses.

**Edge styling** follows the greta convention: solid arrows indicate
deterministic relationships, dashed arrows indicate stochastic
relationships (sampling from a distribution). By default, outgoing edges
from distribution-like nodes (hyperparameter, obs_error, process_error)
are rendered as dashed.

**Auto-classification** proceeds in three layers: (1) topological, where
terminal children become `data` and terminal parents become
`hyperparameter`; (2) naming heuristics, where regex patterns on node
names override topology (e.g., `sigma_obs` maps to `obs_error`); (3)
explicit `node_types`, which take highest priority.

**Plate notation.** When `plates = TRUE`, nodes are enclosed in labeled
grouping boxes. The DiagrammeR backend uses native Graphviz subgraph
clusters (rendered as dashed rectangles with tinted backgrounds). The
visNetwork backend draws plates via JavaScript canvas events.

To save a static image, pass an `.html` file path via `file`. For
PNG/PDF conversion, open the HTML in a browser and use a screenshot
tool, or use the webshot2 package.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Format 1: named list of parent nodes ---
plot_dag(list(y = c("beta", "sigma"),
              beta = c("mu_beta", "tau_beta")))

# --- Format 2: list of R formulas ---
plot_dag(list(y ~ beta + sigma, beta ~ mu_beta + tau_beta))

# --- Format 3: text/LaTeX notation ---
# Plain text notation (one sampling statement per line)
plot_dag("
  y ~ Normal(mu, sigma)
  mu ~ Normal(mu_0, sigma_0)
  sigma ~ HalfCauchy(25)
  mu_0 ~ Normal(0, 100)
  sigma_0 ~ HalfCauchy(5)
")

# LaTeX notation (copy-pasted from a paper or vignette)
plot_dag("
  $$y_i \\sim \\mathcal{N}(\\mu, \\sigma^2)$$
  $$\\mu \\sim \\mathcal{N}(\\mu_0, \\sigma_0^2)$$
  $$\\sigma \\sim \\mathcal{HC}(25)$$
  $$\\mu_0 \\sim \\mathcal{N}(0, 100)$$
  $$\\sigma_0 \\sim \\mathcal{HC}(5)$$
")

# Hierarchical model from LaTeX
plot_dag("
  y_i ~ Normal(theta_i, sigma)
  theta_i ~ Normal(mu, tau)
  mu ~ Normal(0, 100)
  tau ~ HalfCauchy(25)
  sigma ~ HalfCauchy(25)
", title = "Hierarchical normal")

# --- Format 4: lucifer Model function ---
Model <- function(parm, Data) {
    mu <- parm[1]
    sigma <- interval(parm[2], 1e-100, Inf)
    parm[2] <- sigma
    mu.prior <- dnorm(mu, 0, 1000, log = TRUE)
    sigma.prior <- dhalfcauchy(sigma, 25, log = TRUE)
    LL <- sum(dnorm(Data$y, mu, sigma, log = TRUE))
    LP <- LL + mu.prior + sigma.prior
    Monitor <- LP
    yhat <- rnorm(length(Data$y), mu, sigma)
    return(list(LP = LP, Dev = -2 * LL, Monitor = Monitor,
                yhat = yhat, parm = parm))
}
plot_dag(Model, title = "Normal model (auto-extracted)")

# From a fitted demonoid object
# fit <- lucifer(Model, Data, ...)
# plot_dag(fit)

# --- Customization options ---
# State-space model with explicit node types
plot_dag(
  list(y = c("z", "sigma_obs"),
       z = c("z_prev", "sigma_proc", "phi"),
       phi = c("mu_phi", "tau_phi")),
  node_types = c(y = "data", z = "latent", z_prev = "latent",
                 sigma_obs = "obs_error", sigma_proc = "process_error"),
  title = "State-space model"
)

# Hierarchical model with custom plate groups
plot_dag(
  list(y_i = c("theta_i", "sigma"),
       theta_i = c("mu", "tau"),
       mu = c("mu_0", "sigma_0"),
       tau = c("a_0", "b_0")),
  groups = list("Observation" = c("y_i", "sigma"),
                "Group" = "theta_i",
                "Population" = c("mu", "tau"),
                "Hyperpriors" = c("mu_0", "sigma_0", "a_0", "b_0")),
  title = "Hierarchical normal"
)

# Interactive visNetwork backend
plot_dag(list(y ~ beta + sigma), backend = "visnetwork")

# Disable Greek letter rendering
plot_dag(list(y = c("beta", "sigma")), greek = FALSE)

# Hide legend
plot_dag(list(y ~ beta + sigma), legend = FALSE)

# Save to file (multiple formats)
plot_dag(list(y ~ beta + sigma), file = "dag.png")
plot_dag(list(y ~ beta + sigma), file = "dag.pdf",
         width = 800, height = 600)
plot_dag(list(y ~ beta + sigma), file = "dag.svg")
plot_dag(list(y ~ beta + sigma), file = "dag.html")
} # }
```
