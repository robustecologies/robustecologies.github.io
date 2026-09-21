# ============================ #
# The Lyapunov exponent of the logistic map at r = 4 ####
# ============================ #
#
# Claim. The exponent belongs to the pair formed by a map and an invariant measure, and not to
# the map alone. For T(x) = 4x(1 - x) on [0, 1] the absolutely continuous invariant measure and
# every periodic orbit except the fixed point at the origin give log 2, while the Dirac measure
# at that fixed point gives log 4.
#
# The finite-time exponent has a closed form, derived here and not taken from a source. Write
#
#     H(x) = log(pi sqrt(x (1 - x))).
#
# Then log|T'(x)| = log 2 + H(T(x)) - H(x) for every x in (0, 1), as the algebra shows: with
# T(x) = 4x(1 - x) and 1 - T(x) = (1 - 2x)^2, the right side is log 2 + log 2 + log|1 - 2x|,
# which is log|4 - 8x|. Summing along an orbit makes the middle terms cancel, so
#
#     (1/n) sum_{i<n} log|T'(x_i)| = log 2 + (H(x_n) - H(x_0)) / n.
#
# The deviation of a finite-time exponent from log 2 is therefore known exactly for the orbit in
# hand, it falls like 1/n and not like the 1/sqrt(n) of a generic observable, and it is large for
# an orbit that comes near 0 or 1, where H is singular. The same singularity is why the fixed
# point at the origin escapes the value log 2: the conjugacy to the tent map, x = h(y) with
# h(y) = sin^2(pi y / 2), is a diffeomorphism except at the two endpoints, where h' vanishes.
#
# Method and references, none of which shares the code under test.
#   1. The pointwise identity, by evaluating its two sides, which are different expressions.
#   2. The telescoped form, against the finite-time exponent of computed orbits.
#   3. The periodic orbits, from the exact periodic points k/(2^p - 1) of the tent map carried
#      over by h, at periods 2 to 8.
#   4. The conjugacy, by the identity T(h(y)) = h(tent(y)) at many points.
#   5. The growth of the discrepancy between the two sides of the conjugacy over k steps, whose
#      rate is the exponent measured a second way, by a fit that uses no derivative at all.
#   6. A double-precision orbit against a reference carried at 800 bits, which departs from it
#      while the two exponents stay together.
#
# Degenerate cases: the two fixed points, and the tent map itself, whose derivative has modulus
# two everywhere so that every invariant measure gives log 2.

set.seed(20260920L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

LOG2 <- log(2)
Tmap <- function(x) 4 * x * (1 - x)
Tder <- function(x) 4 - 8 * x
tent <- function(y) 1 - abs(2 * y - 1)
h <- function(y) sin(pi * y / 2)^2
H <- function(x) log(pi * sqrt(x * (1 - x)))

# ============================ #
# Test 1: the pointwise identity ####
# ============================ #
#
# The two sides are built from different expressions, log|4 - 8x| on one and a difference of
# logarithms of the conjugacy on the other, so evaluating both tests the derivation and not the
# typing.

test_identity <- function(n_grid = 200001L) {
    x <- seq(0, 1, length.out = n_grid)
    x <- x[x > 1e-9 & x < 1 - 1e-9]
    lhs <- log(abs(Tder(x)))
    rhs <- LOG2 + H(Tmap(x)) - H(x)
    # Beside the half, T(x) is within rounding of 1 and the subtraction 1 - T(x) inside H loses
    # every digit, so the naive evaluation of the right side is inaccurate there while the
    # identity is not. The stable form uses the factorisation 1 - T(x) = (1 - 2x)^2, which is the
    # one algebraic step of the derivation, and is exact in floating point.
    rhs_stable <- LOG2 + (log(pi) + 0.5 * log(Tmap(x)) + log(abs(1 - 2 * x))) - H(x)
    keep <- is.finite(lhs) & is.finite(rhs) & is.finite(rhs_stable)
    # The size of the cancellation is known in advance: at a distance d from the half the term
    # 1 - T(x) is 4d^2 and is computed with an absolute error of order the machine epsilon, so the
    # half logarithm of it carries an error of order eps / (8 d^2). The residual of the naive form
    # multiplied by d^2 is therefore of the order of the machine epsilon, whatever the window.
    resid <- abs(lhs - rhs)
    d <- abs(x - 0.5)
    away <- keep & d > 1e-2
    win1 <- keep & d > 1e-3 & d < 3e-3
    win2 <- keep & d > 1e-2 & d < 3e-2
    list(n = sum(keep), n_away = sum(away),
         max_abs = max(resid[keep]),
         max_abs_away = max(resid[away]),
         max_abs_stable = max(abs(lhs[keep] - rhs_stable[keep])),
         worst_at = x[keep][which.max(resid[keep])],
         scaled_1 = max(resid[win1] * d[win1]^2),
         scaled_2 = max(resid[win2] * d[win2]^2))
}

# ============================ #
# Test 2: the telescoped deviation ####
# ============================ #

finite_time <- function(x0, n) {
    x <- x0; s <- 0
    for (i in seq_len(n)) { s <- s + log(abs(Tder(x))); x <- Tmap(x) }
    list(lambda = s / n, x_end = x)
}

test_telescope <- function(starts, ns) {
    rows <- list()
    for (x0 in starts) {
        for (n in ns) {
            r <- finite_time(x0, n)
            predicted <- LOG2 + (H(r$x_end) - H(x0)) / n
            rows[[length(rows) + 1L]] <- data.frame(
                x0 = x0, n = n, lambda = r$lambda, predicted = predicted,
                residual = abs(r$lambda - predicted),
                bound = abs(H(r$x_end) - H(x0)) / n, stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 3: the exponent of each invariant measure ####
# ============================ #
#
# The periodic points of the tent map of period dividing p are the rationals k/(2^p - 1), and h
# carries them to the periodic points of the logistic map. The exponent of the measure spread
# evenly over one orbit is the average of log|T'| over its points.

test_measures <- function(p_max = 8L) {
    rows <- list(data.frame(measure = "Dirac at the fixed point 0", points = 1L,
                            lambda = log(abs(Tder(0))), stringsAsFactors = FALSE),
                 data.frame(measure = "Dirac at the fixed point 3/4", points = 1L,
                            lambda = log(abs(Tder(0.75))), stringsAsFactors = FALSE))
    for (p in 2:p_max) {
        ys <- (1:(2^p - 2)) / (2^p - 1)
        xs <- h(ys)
        lam <- vapply(xs, function(x0) finite_time(x0, p)$lambda, numeric(1))
        rows[[length(rows) + 1L]] <- data.frame(
            measure = sprintf("Periodic orbits of period %d", p), points = length(xs),
            lambda = mean(lam), stringsAsFactors = FALSE)
        attr(rows[[length(rows)]], "spread") <- max(abs(lam - LOG2))
    }
    out <- do.call(rbind, rows)
    out$spread <- c(NA, NA, vapply(rows[-(1:2)], function(d) attr(d, "spread"), numeric(1)))
    out
}

# ============================ #
# Test 4 and 5: the conjugacy, and the exponent measured through it ####
# ============================ #

test_conjugacy <- function(n_grid = 40001L, k_max = 44L) {
    y <- seq(0, 1, length.out = n_grid)
    one_step <- max(abs(Tmap(h(y)) - h(tent(y))))
    y0 <- 0.3141592653589793
    a <- h(y0); b <- y0; d <- numeric(k_max)
    for (k in seq_len(k_max)) { a <- Tmap(a); b <- tent(b); d[k] <- abs(a - h(b)) }
    # the discrepancy grows at the rate of the exponent, before it saturates at the size of the
    # interval; the window stops where the discrepancy passes a hundredth
    kk <- which(d > 1e-13 & d < 1e-2)
    slope <- if (length(kk) > 4L) unname(coef(lm(log(d[kk]) ~ kk))[2]) else NA_real_
    list(one_step = one_step, slope = slope, n_fitted = length(kk),
         tent_derivative = max(abs(abs(diff(tent(seq(0, 0.5, length.out = 5001)))) /
                                   diff(seq(0, 0.5, length.out = 5001)) - 2)))
}

# ============================ #
# Test 6: double precision against a reference at 800 bits ####
# ============================ #

test_precision <- function(n = 2000L, prec = 800L) {
    if (!requireNamespace("Rmpfr", quietly = TRUE)) return(NULL)
    x0 <- 0.2345678
    xe <- Rmpfr::mpfr(x0, precBits = prec)
    four <- Rmpfr::mpfr(4L, precBits = prec)
    one <- Rmpfr::mpfr(1L, precBits = prec)
    xd <- x0; se <- Rmpfr::mpfr(0L, precBits = prec); sd_ <- 0
    for (i in seq_len(n)) {
        se <- se + log(abs(four - 8 * xe)); xe <- four * xe * (one - xe)
        sd_ <- sd_ + log(abs(Tder(xd))); xd <- Tmap(xd)
    }
    list(gap = abs(as.numeric(xe) - xd),
         lambda_exact = as.numeric(se) / n, lambda_double = sd_ / n,
         bound_exact = abs(H(as.numeric(xe)) - H(x0)) / n,
         bound_double = abs(H(xd) - H(x0)) / n)
}

# ============================ #
# Run ####
# ============================ #

starts <- c(0.2345678, 0.1111111, 0.618034, 0.8090170, 0.4272198)
ns <- c(10L, 100L, 1000L, 10000L)
t_id <- test_identity()
t_tel <- test_telescope(starts, ns)
t_mu <- test_measures()
t_cj <- test_conjugacy()
t_pr <- test_precision()

# The identity holds; away from the half the naive evaluation confirms it to rounding, and the
# stable evaluation confirms it everywhere. The residual of the naive form beside the half is
# reported and is a cancellation and not a failure.
ok_id <- t_id$max_abs_away < 1e-12 && t_id$max_abs_stable < 1e-13 &&
         t_id$scaled_1 < 1e-15 && t_id$scaled_2 < 1e-15
# The telescoped form is an identity, so the residual is rounding alone; the tolerance allows one
# part in 10^9 of the sum, which a run of 10^4 steps in double precision can accumulate.
ok_tel <- max(t_tel$residual) < 1e-9
# Every exponent within its own certified bound of log 2, the bound being the one the identity
# supplies for that orbit.
# The identity makes this an equality, so the margin allowed is the rounding of the sum, the same
# one part in 10^9 that the test above allows.
ok_bound <- all(abs(t_tel$lambda - LOG2) <= t_tel$bound + 1e-9)
ok_mu <- abs(t_mu$lambda[t_mu$measure == "Dirac at the fixed point 0"] - log(4)) < 1e-12 &&
         abs(t_mu$lambda[t_mu$measure == "Dirac at the fixed point 3/4"] - LOG2) < 1e-12 &&
         max(t_mu$spread, na.rm = TRUE) < 1e-8
ok_cj <- t_cj$one_step < 1e-14 && abs(t_cj$slope - LOG2) < 0.05 * LOG2 && t_cj$tent_derivative < 1e-9
ok_pr <- is.null(t_pr) || (t_pr$gap > 1e-3 &&
                           abs(t_pr$lambda_double - LOG2) <= t_pr$bound_double + 1e-8 &&
                           abs(t_pr$lambda_exact - LOG2) <= t_pr$bound_exact + 1e-8)

status <- if (ok_id && ok_tel && ok_bound && ok_mu && ok_cj && ok_pr) "pass" else "fail"

# ============================ #
# Figure ####
# ============================ #
#
# Left: the finite-time exponent of five orbits against the number of steps, inside the band that
# the closed form certifies for each of them. Right: the exponent of each invariant measure
# tested, which is log 2 for all of them but the Dirac measure at the origin.

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    run <- do.call(rbind, lapply(starts, function(x0) {
        x <- x0; s <- 0; out <- list()
        for (i in seq_len(10000L)) {
            s <- s + log(abs(Tder(x))); x <- Tmap(x)
            if (i %in% unique(round(10^seq(0.7, 4, length.out = 60)))) {
                out[[length(out) + 1L]] <- data.frame(
                    x0 = kb_unicode(sprintf("$x_0 = %.4f$", x0)), n = i, lambda = s / i,
                    bound = abs(H(x) - H(x0)) / i, stringsAsFactors = FALSE)
            }
        }
        do.call(rbind, out)
    }))
    run$lo <- LOG2 - run$bound; run$hi <- LOG2 + run$bound
    pa <- ggplot(run, aes(n, lambda, colour = x0)) +
        geom_ribbon(aes(ymin = lo, ymax = hi, group = x0), fill = "grey80", colour = NA,
                    alpha = 0.35) +
        geom_hline(yintercept = LOG2, colour = "grey35", linewidth = 0.4) +
        geom_line(linewidth = 0.5) +
        scale_x_log10() +
        coord_cartesian(ylim = c(LOG2 - 0.22, LOG2 + 0.22)) +
        scale_colour_viridis_d(option = "viridis", end = 0.85, name = NULL) +
        labs(subtitle = kb_unicode("Finite-time exponent of five orbits, inside the band $(H(x_n) - H(x_0))/n$ that the closed form certifies"),
             x = kb_tex("Steps $n$"), y = kb_tex("$\\lambda_n$")) +
        kb_theme() + guides(colour = guide_legend(nrow = 2))
    bar <- t_mu
    bar$label <- kb_unicode(c("$\\delta$ at 0", "$\\delta$ at 3/4", sprintf("Period %d", 2:8)))
    bar$label <- factor(bar$label, levels = rev(bar$label))
    pb <- ggplot(bar, aes(lambda, label)) +
        geom_vline(xintercept = LOG2, colour = "#056796", linewidth = 0.5, linetype = "22") +
        geom_vline(xintercept = log(4), colour = "#be1117", linewidth = 0.5, linetype = "22") +
        geom_point(size = 2.6, colour = "#21918c") +
        scale_x_continuous(limits = c(0.6, 1.45),
                           breaks = c(LOG2, log(4)), labels = kb_ticks(c("$\\log 2$", "$\\log 4$"))) +
        labs(subtitle = kb_unicode("Exponent of each invariant measure tested"),
             x = kb_tex("$\\lambda$"), y = NULL) +
        kb_theme()
    fig <- (pa | pb) + plot_layout(widths = c(1.25, 1)) +
        plot_annotation(
            title = kb_unicode("The exponent belongs to the measure and not to the map alone"),
            caption = kb_caption(paste(
                "Left: the finite-time exponent $\\lambda_n$ of five orbits of $T(x) = 4x(1-x)$, with the grey",
                "band giving the deviation that the identity $\\lambda_n = \\log 2 + (H(x_n) - H(x_0))/n$",
                "certifies for that orbit, where $H(x) = \\log(\\pi\\sqrt{x(1-x)})$. The band is wide where the",
                "orbit has come near 0 or 1, which is where $H$ is singular. Right: every periodic orbit of",
                "period 2 to 8 and the fixed point at $3/4$ give $\\log 2$, and the Dirac measure at the fixed",
                "point 0 gives $\\log 4$; the conjugacy to the tent map is a diffeomorphism except at the two",
                "endpoints, which is why that one measure escapes. Drawn by",
                "checks/lyapunov-exponent-of-the-logistic-map.R.")),
            theme = kb_theme())
    kb_save(fig, "lyapunov-exponent-of-the-logistic-map", width = 9.2, height = 4.4)
}

emit("lyapunov-exponent-of-the-logistic-map", status,
     "The finite-time exponent of the logistic map at r = 4 is log 2 plus an exact coboundary over n, every periodic orbit and the absolutely continuous measure give log 2, and the Dirac measure at the origin gives log 4",
     list(identity_points = t_id$n,
          identity_max_residual_away_from_one_half = t_id$max_abs_away,
          identity_max_residual_stable_form = t_id$max_abs_stable,
          identity_worst_naive_residual = t_id$max_abs,
          identity_worst_naive_point = t_id$worst_at,
          identity_residual_times_distance_squared_near = t_id$scaled_1,
          identity_residual_times_distance_squared_far = t_id$scaled_2,
          telescope_cases = nrow(t_tel),
          telescope_max_residual = max(t_tel$residual),
          exponent_at_10000_steps = t_tel$lambda[t_tel$n == 10000][1],
          certified_bound_at_10000_steps = t_tel$bound[t_tel$n == 10000][1],
          log_two = LOG2,
          lambda_dirac_at_zero = t_mu$lambda[1],
          log_four = log(4),
          lambda_dirac_at_three_quarters = t_mu$lambda[2],
          periodic_orbits_tested = sum(t_mu$points[-(1:2)]),
          greatest_departure_of_a_periodic_orbit = max(t_mu$spread, na.rm = TRUE),
          conjugacy_one_step_residual = t_cj$one_step,
          exponent_from_the_growth_of_the_discrepancy = t_cj$slope,
          points_fitted = t_cj$n_fitted,
          gap_between_double_and_reference = if (is.null(t_pr)) NA_real_ else t_pr$gap,
          lambda_double_precision = if (is.null(t_pr)) NA_real_ else t_pr$lambda_double,
          lambda_reference = if (is.null(t_pr)) NA_real_ else t_pr$lambda_exact))
