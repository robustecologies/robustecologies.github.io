# ============================ #
# Pullback statistical attraction without pullback measure attraction ####
# ============================ #
#
# Claims, with the locators of oljaca2024 (Example 3.5):
#
# The scalar linear equation x' = a(t) x has a(t) = -1 on I_k = [-2^k + 1, -2^(k-1)) for
# k >= 1 and a(t) = 2^(k-1) on J_k = [-2^k, -2^k + 1) for k >= 0, so its process is
# Phi(t, t0)(x0) = x0 exp(int_{t0}^{t} a(s) ds).
#
# 1. The intervals I_k and J_k tile the half line (-infinity, 0).
# 2. Phi(0, -2^k)(x0) = x0 exp(k + 1/2). Reference: the exact integral of the piecewise
#    constant a over the tiling, and, independently, adaptive numerical quadrature.
# 3. For t in the interval Ibar_k = [-2^k + 1, -2^(k-1) - 2k] the process has
#    Phi(0, t)(x0) = x0 exp(k - 1/2 + 2^(k-1) + t), whose supremum over Ibar_k is
#    |x0| exp(-k - 1/2) and tends to 0. The set T of these intervals for k >= 5 therefore
#    pullback attracts every bounded set along times in T.
# 4. T has full density at -infinity: m(T intersect [-s, 0]) / s tends to 1. The local
#    minima of that ratio sit at the right end of each gap, s = 2^k + 2k + 2, where the
#    closed form (2^k - 16 - (k - 4)(k + 6)) / (2^k + 2k + 2) applies; the script compares
#    that closed form with a direct measure computation and takes the infimum over s.
# 5. No nonautonomous set N with liminf m(N(t)) > 0 is pullback attracted to 0 along all
#    times: at t = -2^j the image of N(-2^j) is stretched by exp(j + 1/2), so staying
#    within eps of 0 forces m(N(-2^j)) <= 2 eps exp(-j - 1/2), which tends to 0.
# 6. Degenerate cases: x0 = 0 is fixed, and t = t0 gives the identity.

set.seed(20260915L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) {
        if (is.character(v)) return(sprintf("\"%s\"", v))
        if (is.logical(v)) return(if (v) "true" else "false")
        if (!is.finite(v)) return("null")
        formatC(v, digits = 6, format = "g")
    }
    fields <- paste(sprintf("\"%s\": %s", names(metrics),
                            vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n",
                id, status, summary, fields))
}

# Coefficient a(t) for t < 0, defined by the interval that contains t.
a_of <- function(t) {
    vapply(t, function(u) {
        if (u >= 0) return(0)
        k <- max(0, ceiling(log2(-u)))      # smallest k >= 0 with -2^k <= u
        if (u < -2^k + 1) 2^(k - 1) else -1
    }, numeric(1))
}

# ============================ #
# Test 1: the intervals tile the negative half line ####
# ============================ #
#
# Each point of (-2^K, 0) must belong to exactly one I_k or J_k. The test samples points
# and counts memberships; the count must be one everywhere.

K <- 12L
pts <- -sort(runif(4000L, 1e-6, 2^K - 1e-6))
memberships <- vapply(pts, function(u) {
    inI <- sum(vapply(1:(K + 1), function(k) u >= -2^k + 1 && u < -2^(k - 1), logical(1)))
    inJ <- sum(vapply(0:(K + 1), function(k) u >= -2^k && u < -2^k + 1, logical(1)))
    inI + inJ
}, numeric(1))
tiling_errors <- sum(memberships != 1)

# ============================ #
# Test 2: the exponent at the times -2^k ####
# ============================ #
#
# Exact sum over the tiling: the intervals I_i (i = 1..k) have length 2^(i-1) - 1 and
# coefficient -1, and the intervals J_i (i = 0..k) have length 1 and coefficient 2^(i-1).

# Exact integral of a over [t0, 0], by summing the overlap of each interval of the tiling
# with [t0, 0]. It is independent of the closed forms that it is used to test.
integral_a <- function(t0, kmax = 60L) {
    total <- 0
    for (k in 0:kmax) {
        seg <- function(lo, hi, coef) {
            lo <- max(lo, t0)
            if (hi > lo) total <<- total + coef * (hi - lo)
        }
        seg(-2^k, -2^k + 1, 2^(k - 1))
        if (k >= 1) seg(-2^k + 1, -2^(k - 1), -1)
    }
    total
}
ks <- 1:12
err_exponent <- max(abs(vapply(ks, function(k) integral_a(-2^k), numeric(1)) - (ks + 0.5)))

# Independent route: adaptive quadrature of the piecewise constant coefficient, with the
# jump points supplied as subdivisions. Tolerance 1e-8 for the quadrature.
err_quadrature <- 0
for (k in 1:8) {
    breaks <- sort(unique(c(-2^k, as.vector(rbind(-2^(0:k), -2^(0:k) + 1)), -2^((0:k) - 1), 0)))
    breaks <- breaks[breaks >= -2^k & breaks <= 0]
    value <- sum(vapply(seq_len(length(breaks) - 1), function(i) {
        integrate(function(u) a_of(u), breaks[i], breaks[i + 1],
                  subdivisions = 1000L, rel.tol = 1e-10)$value
    }, numeric(1)))
    err_quadrature <- max(err_quadrature, abs(value - (k + 0.5)))
}

growth_k10 <- exp(integral_a(-2^10))

# ============================ #
# Test 3: contraction along the intervals Ibar_k ####
# ============================ #
#
# For t in Ibar_k the exponent is k - 1/2 + 2^(k-1) + t, which is linear with slope 1, so
# its maximum over Ibar_k is at the right end t = -2^(k-1) - 2k and equals -k - 1/2.

Ibar <- function(k) c(-2^k + 1, -2^(k - 1) - 2 * k)

# Closed form of the exponent for t0 in I_k, from the source.
exponent_from <- function(t0) {
    k <- max(1, ceiling(log2(-t0)))
    stopifnot(t0 >= -2^k + 1, t0 < -2^(k - 1))
    (k - 1 + 0.5) + 2^(k - 1) + t0
}
# The closed form must agree with the exact integral over the tiling.
err_ibar_exact <- max(vapply(5:20, function(k) {
    ts <- seq(Ibar(k)[1], Ibar(k)[2], length.out = 20L)
    max(abs(vapply(ts, exponent_from, numeric(1)) - vapply(ts, integral_a, numeric(1))))
}, numeric(1)))

err_ibar <- 0
sup_ibar <- numeric(0)
for (k in 5:20) {
    ends <- Ibar(k)
    ts <- seq(ends[1], ends[2], length.out = 50L)
    ex <- vapply(ts, exponent_from, numeric(1))
    err_ibar <- max(err_ibar, abs(max(ex) - (-k - 0.5)))
    sup_ibar <- c(sup_ibar, exp(max(ex)))
}
sup_ibar_k5 <- sup_ibar[1]
sup_ibar_k20 <- sup_ibar[length(sup_ibar)]

# Cross-check the exponent formula on Ibar_k against quadrature for one moderate k.
k_test <- 7L
t_test <- mean(Ibar(k_test))
breaks <- sort(unique(c(t_test, as.vector(rbind(-2^(0:k_test), -2^(0:k_test) + 1)),
                        -2^((0:k_test) - 1), 0)))
breaks <- breaks[breaks >= t_test & breaks <= 0]
quad_ibar <- sum(vapply(seq_len(length(breaks) - 1), function(i) {
    integrate(function(u) a_of(u), breaks[i], breaks[i + 1],
              subdivisions = 1000L, rel.tol = 1e-10)$value
}, numeric(1)))
err_ibar_quadrature <- abs(quad_ibar - exponent_from(t_test))

# ============================ #
# Test 4: T has full density at -infinity ####
# ============================ #

measure_T <- function(s, kmax = 60L) {
    total <- 0
    for (k in 5:kmax) {
        ends <- Ibar(k)
        lo <- max(ends[1], -s)
        hi <- ends[2]
        if (hi > lo) total <- total + (hi - lo)
    }
    total
}

# Closed form at the right end of each gap, s = 2^k + 2k + 2, where the covered part is
# exactly the union of Ibar_5, ..., Ibar_k.
density_closed <- function(k) (2^k - 16 - (k - 4) * (k + 6)) / (2^k + 2 * k + 2)
ks_gap <- 5:40
err_density_closed <- max(abs(vapply(ks_gap, function(k) measure_T(2^k + 2 * k + 2) / (2^k + 2 * k + 2),
                                     numeric(1)) - vapply(ks_gap, density_closed, numeric(1))))

# Infimum of the density over s >= S, evaluated on the candidate minima (gap ends) and on a
# dense grid of other times, for two windows.
density_at <- function(s) measure_T(s) / s
grid_s <- function(lo, hi) sort(unique(c(exp(seq(log(lo), log(hi), length.out = 4000L)),
                                        2^(5:40) + 2 * (5:40) + 2, 2^(5:40) - 1)))
window <- function(lo, hi) {
    ss <- grid_s(lo, hi)
    ss <- ss[ss >= lo & ss <= hi]
    min(vapply(ss, density_at, numeric(1)))
}
inf_density_from_2e2 <- window(2e2, 1e12)
inf_density_from_1e6 <- window(1e6, 1e12)
inf_density_from_1e9 <- window(1e9, 1e12)

# ============================ #
# Test 5: no pullback measure attraction ####
# ============================ #
#
# A nonautonomous set N whose image at time 0 lies within eps of the origin must satisfy
# m(N(-2^j)) <= 2 eps exp(-j - 1/2). The bound is reported for eps = 1 and j = 1..20.

eps <- 1
needed_measure <- 2 * eps * exp(-(1:20) - 0.5)
needed_measure_j20 <- needed_measure[20]

# ============================ #
# Test 6: degenerate cases ####
# ============================ #

Phi <- function(t0, x0) x0 * exp(exponent_from(t0))
err_zero <- abs(Phi(mean(Ibar(6L)), 0))
err_same_time <- abs(integral_a(0))              # empty interval, so the process is the identity

ok <- tiling_errors == 0 && err_exponent < 1e-12 && err_quadrature < 1e-8 &&
    err_ibar < 1e-9 && err_ibar_quadrature < 1e-8 && err_ibar_exact < 1e-9 &&
    sup_ibar_k5 < 0.01 && sup_ibar_k20 < 1e-8 &&
    err_density_closed < 1e-9 &&
    inf_density_from_2e2 > 0.5 && inf_density_from_1e6 > 0.99 &&
    inf_density_from_1e9 > 0.9999 &&
    needed_measure_j20 < 1e-8 && err_zero == 0 && err_same_time < 1e-12

emit("alternating-growth-pullback", if (ok) "pass" else "fail",
     "The linear process with alternating growth contracts every bounded set along a set of times of full density at -infinity and stretches it by exp(k + 1/2) at the times -2^k",
     list(tiling_errors = tiling_errors,
          exponent_max_abs_error = err_exponent,
          quadrature_max_abs_error = err_quadrature,
          growth_factor_k10 = growth_k10,
          ibar_exponent_max_abs_error = err_ibar,
          ibar_quadrature_abs_error = err_ibar_quadrature,
          ibar_exact_integral_max_abs_error = err_ibar_exact,
          sup_contraction_k5 = sup_ibar_k5,
          sup_contraction_k20 = sup_ibar_k20,
          density_closed_form_max_abs_error = err_density_closed,
          inf_density_s_above_200 = inf_density_from_2e2,
          inf_density_s_above_1e6 = inf_density_from_1e6,
          inf_density_s_above_1e9 = inf_density_from_1e9,
          required_measure_at_j20 = needed_measure_j20,
          zero_solution_error = err_zero,
          identity_interval_error = err_same_time))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the exponent of the process from time t up to 0, against the logarithm of the elapsed
# time. It is a sawtooth: at the times -2^k it reaches k + 1/2, so the process stretches by
# exp(k + 1/2), while on the intervals Ibar_k it falls to -k - 1/2 and below, so the process
# contracts by exp(-k - 1/2) or more. Right: the share of the window [-s, 0] covered by those
# intervals, which tends to one, so they form a set of times of full density at minus infinity.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    kmax_fig <- 12L
    ts_fig <- -sort(10^seq(log10(2), log10(2^kmax_fig - 1), length.out = 1200))
    saw <- data.frame(t = ts_fig, e = vapply(ts_fig, integral_a, numeric(1)))
    peaks <- data.frame(t = -2^(1:kmax_fig), e = (1:kmax_fig) + 0.5)
    bands <- do.call(rbind, lapply(2:kmax_fig, function(k) {
        ends <- Ibar(k)
        data.frame(xmin = -ends[1], xmax = -ends[2], k = k)
    }))
    bands <- bands[bands$xmax > 0 & bands$xmin > bands$xmax, ]
    p1 <- ggplot(saw, aes(-t, e)) +
        geom_rect(data = bands, aes(xmin = xmax, xmax = xmin, ymin = -Inf, ymax = Inf),
                  inherit.aes = FALSE, fill = "#21918c", alpha = 0.12) +
        geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
        geom_line(colour = "#440154", linewidth = 0.45) +
        geom_point(data = peaks, aes(-t, e), inherit.aes = FALSE, colour = "#B8390E", size = 1.5) +
        annotate("text", x = 2.6, y = 25, hjust = 0, size = 2.6, colour = "#B8390E",
                 label = kb_unicode("Red: the times $-2^k$, where the exponent is $k + 1/2$")) +
        annotate("text", x = 2.6, y = -25, hjust = 0, size = 2.6, colour = "#12726e",
                 label = kb_unicode("Shaded: the intervals $\\bar{I}_k$, where the exponent is at most $-k - 1/2$")) +
        scale_x_log10() +
        coord_cartesian(ylim = c(-32, 32)) +
        labs(title = "Stretching at a thin set of times, contraction on the rest",
             subtitle = "Exponent from time t up to 0; inside the shaded intervals it falls far below the frame",
             x = kb_tex("Elapsed time $t$, backwards from zero"), y = kb_tex("Exponent $\\lambda$")) +
        kb_theme()
    ss <- 2^(5:40) + 2 * (5:40) + 2
    dens <- data.frame(s = ss, share = vapply(ss, density_at, numeric(1)))
    p2 <- ggplot(dens, aes(s, share)) +
        geom_hline(yintercept = 1, colour = "#B8390E", linewidth = 0.4, linetype = "22") +
        geom_line(colour = "#21918c", linewidth = 0.5) +
        geom_point(colour = "#21918c", size = 1.2) +
        scale_x_log10() + scale_y_continuous(limits = c(0, 1.05)) +
        labs(title = "Those intervals have full density at minus infinity",
             subtitle = "Share of the window from -s to 0 that they cover, at the worst instants",
             x = kb_tex("Length $s$ of the window"), y = "Share covered") +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.25, 1)) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "The process contracts every bounded set along a set of times of full density and stretches it without bound along the thin set that remains, which is why pullback attraction along a set of full density is weaker than pullback attraction. In the run recorded by checks/alternating-growth-pullback.R the exponent at the times $-2^k$ agreed with $k + 1/2$ to %.1e, the closed form of the covered share agreed with the exact measure to %.1e, and the share stayed above %.4f for every window longer than $10^9$.",
                err_exponent, err_density_closed, inf_density_from_1e9)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "alternating-growth-pullback", width = 9.4, height = 4.2))
}
