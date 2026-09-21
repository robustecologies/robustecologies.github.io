# ============================ #
# Fibrewise contraction in the sine-forced circle map ####
# ============================ #
#
# The map is F(y, x) = (3 y, x + a sin(2 pi x) + y) mod 1 on the two-torus, with a = 1/8, the
# numerical experiment of Homburg (2012), Figure 1. The script tests five claims and prints one
# JSON line for tools/kb.py verify.
#
# 1. Diffeomorphism and partial hyperbolicity. The fibre derivative is f'(x) = 1 + 2 pi a cos(2 pi x),
#    so at a = 1/8 its extremes are 1 -/+ pi/4 exactly. The script compares a numerical minimum and
#    maximum over a fine grid, and a finite-difference derivative of the fibre map itself, with
#    those closed forms; the finite-difference route does not use the formula for f', so the
#    analytic derivative used elsewhere in the script is validated and not assumed. It then reports
#    whether condition (8) of the source, max f' < min g' = 3, holds, and whether the stronger
#    condition of Lemma 3.1, max{f', 1/f'} near 1 with min g' > 2, holds.
# 2. Fibrewise contraction. Twenty states of one fibre are iterated over random base histories and
#    the circular diameter of their images is recorded. The base orbit of y -> 3y mod 1 is computed
#    exactly from a ternary digit string, because iterating 3y in double precision multiplies the
#    rounding error by three at every step and exhausts the mantissa after about 33 steps, which
#    would make a run of 120 steps a computation on rounding noise.
# 3. Independent route to the contraction rate. Two orbits at initial separation 1e-7 are iterated,
#    and the step ratio log(sep_{n+1} / sep_n) is compared with log f'(x_n). The two are computed
#    from different quantities, a distance between two orbits and a derivative along one. By the
#    mean value theorem the step ratio equals log f'(xi) for some xi between the two orbits, so the
#    discrepancy is bounded by |(log f')'| times the separation, and |(log f')'| = |f''| / f' is at
#    most (2 pi)^2 a / (1 - pi/4) = 23.0 at a = 1/8. Only the steps with separation between 1e-9 and
#    1e-7 are compared, so the linearisation error is at most 23 * 1e-7 = 2.3e-6 and the rounding
#    error, which is about 2.2e-16 divided by the separation, is at most 5e-7.
# 4. Degenerate case. At a = 0 every fibre map is the rotation x -> x + y, an isometry of the
#    circle, so the diameter of any finite set is preserved exactly. Expected change: zero.
# 5. Density of one orbit. Ten thousand points of a single orbit are counted on grids of the torus
#    with 400 and with 1024 cells, which is the left panel of Figure 1 of the source. Only the
#    coarser grid carries an assertion, since the finer one is a report of what the sample reaches
#    at that resolution and not a claim about the orbit.
#
# Tolerances. The closed forms of test 1 are compared at 1e-12, which is above the rounding of a
# grid minimum and far below any error that would matter. The finite-difference derivative is
# compared at 1e-6, the size of the truncation error of a central difference at step 1e-5. Test 3
# is compared at 1e-5, above the bound 2.8e-6 derived above and below the values of log f', which
# range over an interval of width 2.1. Test 4 is compared at 1e-14, a few rounding units of the
# modulo operation. Test 2
# has no closed form and asserts only that the median diameter falls below 1e-6, a threshold chosen
# far above the observed median and far below the initial diameter of 0.95.

set.seed(20260918L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) {
        if (is.character(v)) return(sprintf("\"%s\"", v))
        if (is.logical(v)) return(if (v) "true" else "false")
        if (!is.finite(v)) return("null")
        formatC(v, digits = 6, format = "g")
    }
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

amp <- 1 / 8
fibre <- function(x, y, a = amp) (x + a * sin(2 * pi * x) + y) %% 1
dfibre <- function(x, a = amp) 1 + 2 * pi * a * cos(2 * pi * x)

# Exact orbit of y -> 3y mod 1 from a ternary digit string: shifting the digits is the map.
base_orbit <- function(digits) {
    K <- length(digits); p <- 3^-(1:K)
    vapply(0:(K - 1), function(j) sum(digits[(j + 1):K] * p[1:(K - j)]), numeric(1))
}
# Length of the shortest arc of the circle that contains every point of x.
circ_diameter <- function(x) { xs <- sort(x); 1 - max(c(diff(xs), 1 - (xs[length(xs)] - xs[1]))) }

## 1. Derivative extremes, condition (8) and the refinement of Lemma 3.1 -----------------------
grid <- seq(0, 1, length.out = 2000001)
d_min_num <- min(dfibre(grid)); d_max_num <- max(dfibre(grid))
d_min_exact <- 1 - pi / 4; d_max_exact <- 1 + pi / 4
err_extremes <- max(abs(d_min_num - d_min_exact), abs(d_max_num - d_max_exact))
h <- 1e-5
xs_fd <- seq(0, 1, length.out = 1001)
fd <- (fibre(xs_fd + h, 0) - fibre(xs_fd - h, 0))
fd <- ((fd + 0.5) %% 1 - 0.5) / (2 * h)          # unwrap the modulo before dividing
err_fd <- max(abs(fd - dfibre(xs_fd)))
cond8 <- d_max_exact < 3
cond_lemma31 <- max(d_max_exact, 1 / d_min_exact) < 3
ratio_lemma31 <- max(d_max_exact, 1 / d_min_exact)

## 2. Fibrewise contraction over random base histories ------------------------------------------
n_step <- 120L; n_hist <- 300L; n_state <- 20L
diam <- matrix(NA_real_, nrow = n_hist, ncol = n_step + 1L)
for (r in seq_len(n_hist)) {
    ys <- base_orbit(sample(0:2, n_step + 40L, replace = TRUE))
    x <- (seq_len(n_state) - 0.5) / n_state
    diam[r, 1] <- circ_diameter(x)
    for (s in seq_len(n_step)) { x <- fibre(x, ys[s]); diam[r, s + 1L] <- circ_diameter(x) }
}
final <- diam[, n_step + 1L]
med_final <- median(final)
frac_small <- mean(final < 1e-6)
worst <- max(final)

## 3. Separation ratio against the analytic derivative ------------------------------------------
n_sep <- 200L
ys <- base_orbit(sample(0:2, n_sep + 40L, replace = TRUE))
x1 <- 0.2; x2 <- x1 + 1e-7
sep <- numeric(n_sep); lder <- numeric(n_sep)
for (s in seq_len(n_sep)) {
    sep[s] <- min(abs(x2 - x1), 1 - abs(x2 - x1)); lder[s] <- log(dfibre(x1))
    x1 <- fibre(x1, ys[s]); x2 <- fibre(x2, ys[s])
}
keep <- which(sep[-n_sep] <= 1e-7 & sep[-n_sep] >= 1e-9 & sep[-1] >= 1e-9)
ratio <- log(sep[keep + 1L] / sep[keep])
err_rate <- if (length(keep) >= 10L) max(abs(ratio - lder[keep])) else Inf
n_rate_steps <- length(keep)

## 4. Degenerate case: zero amplitude -----------------------------------------------------------
x <- (seq_len(n_state) - 0.5) / n_state; d0 <- circ_diameter(x)
ys <- base_orbit(sample(0:2, n_step + 40L, replace = TRUE))
for (s in seq_len(n_step)) x <- fibre(x, ys[s], a = 0)
err_isometry <- abs(circ_diameter(x) - d0)

## 5. Coverage of the torus by one orbit --------------------------------------------------------
n_orb <- 10000L
ys <- base_orbit(sample(0:2, n_orb + 40L, replace = TRUE))
xo <- numeric(n_orb); xx <- 0.3
for (s in seq_len(n_orb)) { xo[s] <- xx; xx <- fibre(xx, ys[s]) }
cells20 <- length(unique(floor(ys[1:n_orb] * 20) * 20 + floor(xo * 20)))
cells32 <- length(unique(floor(ys[1:n_orb] * 32) * 32 + floor(xo * 32)))

ok <- c(err_extremes < 1e-12, err_fd < 1e-6, cond8, med_final < 1e-6,
        err_rate < 1e-5, err_isometry < 1e-14, cells20 == 400L)
status <- if (all(ok)) "pass" else "fail"

## Figure ---------------------------------------------------------------------------------------
# Left: the fibre diameter against the number of steps for a sample of histories, with the median
# over all of them, on a logarithmic scale. Right: ten thousand points of one orbit on the torus.
if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    floor_v <- 1e-16
    show <- sample(seq_len(n_hist), 60)
    # t(diam[show, ]) has the steps down its columns, so as.vector runs the steps fastest and the
    # history index slowest; the two index vectors below must follow that order.
    long <- data.frame(step = rep(0:n_step, times = length(show)),
                       value = pmax(as.vector(t(diam[show, ])), floor_v),
                       hist = rep(show, each = n_step + 1L))
    med <- data.frame(step = 0:n_step, value = pmax(apply(diam, 2, median), floor_v))
    p1 <- ggplot(long, aes(step, value, group = hist)) +
        geom_line(colour = "#056796", alpha = 0.18, linewidth = 0.25) +
        geom_line(data = med, aes(step, value), inherit.aes = FALSE,
                  colour = "#be1117", linewidth = 0.8) +
        scale_y_log10() +
        labs(title = "Fibre diameter under iteration",
             subtitle = "Sixty base histories in blue, median of three hundred in red",
             x = "Number of steps", y = "Fibre diameter") +
        kb_theme()
    orb <- data.frame(cy = floor(ys[1:n_orb] * 20), cx = floor(xo * 20))
    tab20 <- as.data.frame(table(orb$cy, orb$cx), stringsAsFactors = FALSE)
    names(tab20) <- c("cy", "cx", "count")
    tab20$y <- (as.numeric(tab20$cy) - 0.5) / 20
    tab20$x <- (as.numeric(tab20$cx) - 0.5) / 20
    p2 <- ggplot(tab20, aes(y, x, fill = count)) +
        geom_tile() +
        scale_fill_viridis_c(name = "Points per cell", option = "mako", direction = -1) +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        labs(title = "One orbit on the torus",
             subtitle = "Ten thousand orbit points on four hundred cells, none empty",
             x = "Base coordinate y", y = "Fibre coordinate x") +
        kb_theme()
    p <- patchwork::wrap_plots(p1, p2, widths = c(1.25, 1)) +
        patchwork::plot_annotation(
            caption = kb_caption(paste(
                "Left: the shortest arc containing twenty states of one fibre, against the number of steps, on a logarithmic scale; it falls to rounding level for most histories and remains of order $10^{-2}$ for a few, which is what a limit without a rate looks like over a fixed horizon.",
                "Right: the same map fills the torus, so the fibrewise contraction on the left and the apparent density on the right belong to one system.",
                "A diameter below $10^{-16}$, which includes the histories whose twenty states became equal in double precision, is drawn at that floor.",
                "Drawn by checks/sine-forced-circle-synchronization.R from the arrays the tests use.")),
            theme = ggplot2::theme(plot.caption = ggplot2::element_text(colour = "grey40", size = 7.5, hjust = 0),
                                   plot.caption.position = "plot"))
    kb_save(p, "sine-forced-circle-synchronization", width = 9.2, height = 4.2)
}

emit("sine-forced-circle-synchronization", status,
     "The fibre maps of the sine-forced circle map at a = 1/8 are diffeomorphisms satisfying the partial hyperbolicity condition, and twenty states of one fibre contract to a point for most base histories while one orbit fills the torus",
     list(min_fibre_derivative = d_min_num, max_fibre_derivative = d_max_num,
          error_against_closed_form = err_extremes, error_finite_difference = err_fd,
          condition_8_holds = cond8, lemma_3_1_ratio = ratio_lemma31,
          lemma_3_1_holds = cond_lemma31,
          median_final_diameter = med_final, fraction_below_1e_6 = frac_small,
          worst_final_diameter = worst, initial_diameter = diam[1, 1],
          error_separation_rate = err_rate, steps_compared_for_rate = n_rate_steps, error_zero_amplitude_isometry = err_isometry,
          cells_visited_of_400 = cells20, cells_visited_of_1024 = cells32))
