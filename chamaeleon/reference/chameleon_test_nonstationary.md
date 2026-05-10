# Test for non-stationary chameleon behavior

Detect time-varying chameleon attractor properties using windowed
analysis. This test is designed for non-autonomous systems where the
attractor structure and its scale-dependent properties may evolve over
time.

## Usage

``` r
chameleon_test_nonstationary(
  x,
  window_size = NULL,
  step_size = NULL,
  min_windows = 5,
  chameleon_index = c("D_trend", "theta_trend", "combined"),
  n_surrogates = 0,
  alpha = 0.05,
  test_stationarity = TRUE,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  Numeric vector. Time series to analyze.

- window_size:

  Integer. Number of observations per window. Default is automatically
  determined as `max(1000, length(x)/10)`.

- step_size:

  Integer. Step between consecutive windows (default `window_size/4` for
  75 percent overlap).

- min_windows:

  Integer. Minimum number of windows required (default 5).

- chameleon_index:

  Character. Which index to track over time:

  "D_trend"

  :   (Default) Spearman correlation of D with log-frequency.

  "theta_trend"

  :   Spearman correlation of theta with log-frequency.

  "combined"

  :   Combined index from both D and theta trends.

- n_surrogates:

  Integer. Surrogates per window for local testing (default 0). Set to 0
  to skip local testing and only compute indices, which significantly
  improves performance. The global test provides statistical inference;
  per-window surrogates are rarely needed.

- alpha:

  Numeric. Significance level (default 0.05).

- test_stationarity:

  Logical. Test whether chameleon character is stationary over time
  (default TRUE).

- parallel:

  Logical. Use parallel computation (default TRUE).

- n_cores:

  Integer. Number of cores (default: detectCores() - 1).

- verbose:

  Logical. Print progress (default FALSE).

- ...:

  Additional arguments passed to window-level analysis.

## Value

Object of class `"chameleon_nonstationary"` containing:

- window_results:

  Data frame with per-window chameleon indices and test results (time,
  index value, p-value, classification).

- global_test:

  Overall test for presence of chameleon behavior anywhere in the time
  series.

- stationarity_test:

  Test for whether chameleon character is constant over time (if
  `test_stationarity = TRUE`).

- change_points:

  Detected change points in chameleon behavior (if any).

- summary_stats:

  Summary statistics of the chameleon index over time.

- params:

  List of analysis parameters.

## Details

**Motivation:**

In non-autonomous dynamical systems, the governing equations depend
explicitly on time: \\\dot{x} = f(x, t)\\. This means the attractor
structure can evolve, and chameleon behavior may appear, disappear, or
change character as the system traverses different dynamical regimes.

Examples include:

- Seasonally forced climate oscillators

- Neural systems with varying input drives

- Ecosystems under environmental change

- Financial markets with regime switching

**Analysis procedure:**

1.  Divide time series into overlapping windows

2.  For each window, compute MEMD decomposition, scale-dependent
    metrics, chameleon index (trend statistic), and optionally perform
    local surrogate test

3.  Analyze temporal evolution of chameleon index

4.  Test for stationarity of chameleon character

5.  Detect change points (if present)

**Stationarity test (H0: constant chameleon behavior):**

Tests whether the chameleon index is stationary using:

- Augmented Dickey-Fuller test for unit root

- Runs test for non-random patterns

- Mann-Kendall test for monotonic trend

**Change point detection:**

Uses CUSUM-based detection to identify abrupt transitions between
chameleon and non-chameleon regimes.

## References

Lucarini V, Faranda D, Wouters J (2012). Universal behaviour of extreme
value statistics for selected observables of dynamical systems. Journal
of Statistical Physics 147:63-73.
[doi:10.1007/s10955-012-0468-z](https://doi.org/10.1007/s10955-012-0468-z)

Freitas ACM, Freitas JM, Todd M (2010). Hitting time statistics and
extreme value theory. Probability Theory and Related Fields 147:675-710.
[doi:10.1007/s00440-009-0221-y](https://doi.org/10.1007/s00440-009-0221-y)

## See also

[`chameleon_test`](https://robustecologies.github.io/chamaeleon/reference/chameleon_test.md)
for stationary chameleon testing;
[`chameleon_analysis`](https://robustecologies.github.io/chamaeleon/reference/chameleon_analysis.md)
for single-window analysis;
[`print.chameleon_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/print.chameleon_nonstationary.md),
[`summary.chameleon_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/summary.chameleon_nonstationary.md),
[`plot.chameleon_nonstationary`](https://robustecologies.github.io/chamaeleon/reference/plot.chameleon_nonstationary.md)
for inspection of the returned object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Synthetic non-stationary series: sinusoid with a regime shift halfway
set.seed(42)
n <- 4000
t <- seq_len(n)
regime <- ifelse(t > n / 2, 1.5, 1.0)
x <- regime * sin(2 * pi * t / 80) + 0.1 * rnorm(n)

# Detect time-evolving chameleon behaviour with sliding windows
result <- chameleon_test_nonstationary(
  x,
  window_size = 1000,
  step_size   = 200,
  chameleon_index = "D_trend",
  n_surrogates    = 0,
  parallel        = FALSE,
  verbose         = FALSE
)

print(result)
summary(result)
plot(result)
} # }
```
