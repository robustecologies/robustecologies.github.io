# ============================ #
# Process, pullback limit and blow-up threshold of x' = -exp(-t^2) x^2 ####
# ============================ #
#
# Claims, with the locators of oljaca2024 (Example 2.4):
#
# 1. The equation x' = -exp(-t^2) x^2 generates the process
#    Phi(t, t0)(x0) = 1 / (1/x0 + int_{t0}^{t} exp(-s^2) ds), and the integral equals
#    (sqrt(pi)/2) (erf(t) - erf(t0)). Reference: a fourth-order Runge-Kutta solution of
#    the equation with step 1e-5, which shares no code with the closed form.
# 2. The closed form satisfies the initial value property Phi(s, s) = id and the cocycle
#    property Phi(t, u) o Phi(u, s) = Phi(t, s).
# 3. The pullback limit of Phi(t, t0)(x0) as t0 -> -infinity is
#    xi(t, x0) = x0 / (1 + (sqrt(pi)/2) x0 (erf(t) + 1)), every xi(., c) solves the
#    equation on the whole line where its denominator does not vanish, and the forward
#    limit of xi(t, x0) as t -> infinity is x0 / (1 + x0 sqrt(pi)).
# 4. Correction to the source. The source fixes an arbitrary x0 > 0 and asserts that
#    A(t) = [xi(t, -x0), xi(t, x0)] is a compact invariant nonautonomous set, hence a
#    forward measure attractor. The lower solution xi(., -x0) escapes to -infinity at the
#    finite time t* with erf(t*) = 2 / (sqrt(pi) x0) - 1 whenever x0 > 1/sqrt(pi), so the
#    construction needs x0 <= 1/sqrt(pi) = 0.5641896 for every fibre to be compact, and
#    x0 < 1/sqrt(pi) for the fibres to be bounded uniformly in t. The script locates the
#    escape time by bisection on the denominator and confirms both regimes.
# 5. Distinct entire solutions stay a positive distance apart uniformly in t, so no
#    solution attracts another in the classical pullback or forward sense.
# 6. Degenerate cases: x0 = 0 gives the zero solution, and t = t0 gives the identity.

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

erf <- function(t) 2 * pnorm(t * sqrt(2)) - 1
g <- function(t) exp(-t^2)

# Closed form of the process, written with the error function.
Phi <- function(t, t0, x0) {
    if (isTRUE(all.equal(x0, 0))) return(rep(0, length(t)))
    x0 / (1 + x0 * (sqrt(pi) / 2) * (erf(t) - erf(t0)))
}

# Independent route: classical fourth-order Runge-Kutta on x' = -exp(-t^2) x^2.
rk4 <- function(t0, x0, t1, h = 2e-5) {
    n <- max(1L, as.integer(ceiling(abs(t1 - t0) / h)))
    h <- (t1 - t0) / n
    f <- function(t, x) -g(t) * x^2
    x <- x0
    t <- t0
    for (i in seq_len(n)) {
        k1 <- f(t, x)
        k2 <- f(t + h / 2, x + h * k1 / 2)
        k3 <- f(t + h / 2, x + h * k2 / 2)
        k4 <- f(t + h, x + h * k3)
        x <- x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t <- t + h
    }
    x
}

# ============================ #
# Test 1: the closed form solves the equation ####
# ============================ #
#
# Tolerance 1e-9: with step 2e-5 the observed error of RK4 over an interval of length at
# most 5 is of order 1e-10, two orders below the tolerance.

cases <- expand.grid(x0 = c(-0.4, 0.1, 2), t0 = c(-3, 0), t = c(-1, 2))
cases <- cases[cases$t > cases$t0, ]
err_rk4 <- 0
for (i in seq_len(nrow(cases))) {
    exact <- Phi(cases$t[i], cases$t0[i], cases$x0[i])
    approx <- rk4(cases$t0[i], cases$x0[i], cases$t[i])
    err_rk4 <- max(err_rk4, abs(exact - approx))
}

# ============================ #
# Test 2: initial value and cocycle properties ####
# ============================ #

err_id <- 0
err_cocycle <- 0
for (i in seq_len(200L)) {
    x0 <- runif(1, -0.4, 2)
    s <- runif(1, -3, 0)
    u <- runif(1, s, 1)
    t <- runif(1, u, 3)
    err_id <- max(err_id, abs(Phi(s, s, x0) - x0))
    err_cocycle <- max(err_cocycle, abs(Phi(t, u, Phi(u, s, x0)) - Phi(t, s, x0)))
}

# ============================ #
# Test 3: pullback limit, entire solutions and forward limit ####
# ============================ #

xi <- function(t, c) c / (1 + (sqrt(pi) / 2) * c * (erf(t) + 1))

# erf(-8) differs from -1 by less than 1e-29, so Phi(t, -8, .) equals the pullback limit
# to machine precision; tolerance 1e-12.
err_pullback <- 0
for (c in c(-0.5, -0.2, 0.05, 0.3, 1.5)) {
    for (t in c(-2, -0.5, 0, 1, 4)) {
        err_pullback <- max(err_pullback, abs(Phi(t, -8, c) - xi(t, c)))
    }
}

# Residual of the equation along xi, by central differences with step 1e-5; the
# truncation error of the difference is of order 1e-10, so the tolerance is 1e-8.
resid_xi <- 0
for (c in c(-0.5, -0.2, 0.05, 0.3, 1.5)) {
    for (t in seq(-4, 4, by = 0.25)) {
        h <- 1e-5
        dxi <- (xi(t + h, c) - xi(t - h, c)) / (2 * h)
        resid_xi <- max(resid_xi, abs(dxi + g(t) * xi(t, c)^2))
    }
}

err_forward <- 0
for (c in c(-0.5, 0.05, 0.3, 1.5)) {
    err_forward <- max(err_forward, abs(xi(20, c) - c / (1 + c * sqrt(pi))))
}

# ============================ #
# Test 4: the threshold 1 / sqrt(pi) for the lower solution ####
# ============================ #
#
# The denominator of xi(t, -x0) is 1 - (sqrt(pi)/2) x0 (erf(t) + 1), which decreases from
# 1 to 1 - sqrt(pi) x0. It vanishes at a finite t exactly when x0 > 1 / sqrt(pi).

threshold <- 1 / sqrt(pi)
denom <- function(t, x0) 1 - (sqrt(pi) / 2) * x0 * (erf(t) + 1)

escape_time <- function(x0) {
    if (denom(40, x0) > 0) return(Inf)
    uniroot(function(t) denom(t, x0), lower = -40, upper = 40, tol = 1e-12)$root
}

t_star_above <- escape_time(1)                    # x0 = 1 > 1 / sqrt(pi)
erf_star <- erf(t_star_above)
erf_star_exact <- 2 / (sqrt(pi) * 1) - 1
sup_below <- max(abs(xi(seq(-40, 40, by = 0.01), -0.5)))          # x0 = 0.5 < threshold
sup_at <- max(abs(xi(seq(-40, 40, by = 0.01), -threshold)))       # x0 = threshold
# For x0 < 1 / sqrt(pi) the lower solution is bounded, and its supremum is the forward
# limit x0 / (1 - sqrt(pi) x0); tolerance 1e-9 for the residual erf tail at t = 40.
sup_below_exact <- 0.5 / (1 - sqrt(pi) * 0.5)
err_sup_below <- abs(sup_below - sup_below_exact)
escape_below <- escape_time(0.5)
escape_boundary <- escape_time(threshold)
# One-sided blow-up: the solution grows without bound as t rises to t_star.
blowup_size <- abs(xi(t_star_above - 1e-8, -1))

# ============================ #
# Test 5: solutions stay apart ####
# ============================ #

cs <- c(-0.4, -0.2, 0, 0.2, 0.4)
tgrid <- seq(-6, 6, by = 0.001)                 # a(t) is constant to machine precision outside
a_of <- function(t) (sqrt(pi) / 2) * (erf(t) + 1)
xi_a <- function(a, c) c / (1 + a * c)

# Write the gap between two entire solutions as a function of a = (sqrt(pi)/2)(erf(t) + 1),
# which increases from 0 to sqrt(pi). Then d(gap)/da = -(xi(c2)^2 - xi(c1)^2), so the gap
# is stationary only where xi(c2) = -xi(c1), that is at a* = -(c1 + c2) / (2 c1 c2). The
# infimum over t is therefore the smallest of the values at a = 0, at a = sqrt(pi) and, when
# it lies in that range, at a*. Tolerance 1e-5: the grid spacing in a is at most 1e-3 and the
# gap is quadratic near an interior minimum.
min_gap <- Inf
err_gap_closed_form <- 0
for (i in seq_along(cs)) {
    for (j in seq_along(cs)) {
        if (i >= j) next
        c1 <- cs[i]
        c2 <- cs[j]
        observed <- min(abs(xi_a(a_of(tgrid), c2) - xi_a(a_of(tgrid), c1)))
        candidates <- c(abs(c2 - c1), abs(xi_a(sqrt(pi), c2) - xi_a(sqrt(pi), c1)))
        if (c1 * c2 != 0) {
            a_star <- -(c1 + c2) / (2 * c1 * c2)
            if (a_star > 0 && a_star < sqrt(pi)) {
                candidates <- c(candidates, abs(xi_a(a_star, c2) - xi_a(a_star, c1)))
            }
        }
        predicted <- min(candidates)
        err_gap_closed_form <- max(err_gap_closed_form, abs(observed - predicted))
        min_gap <- min(min_gap, observed)
    }
}

tgrid <- seq(-40, 40, by = 0.05)
# Fibre length of the candidate attractor A(t) = [xi(t, -x0), xi(t, x0)] with x0 = 0.4.
fibre <- xi(tgrid, 0.4) - xi(tgrid, -0.4)
fibre_min <- min(fibre)
fibre_max <- max(fibre)

# ============================ #
# Test 6: degenerate cases ####
# ============================ #

err_zero <- max(abs(Phi(c(-2, 0, 3), -5, 0)))
err_same_time <- abs(Phi(1.25, 1.25, 0.7) - 0.7)

ok <- err_rk4 < 1e-9 && err_id < 1e-14 && err_cocycle < 1e-12 &&
    err_pullback < 1e-12 && resid_xi < 1e-8 && err_forward < 1e-12 &&
    is.finite(t_star_above) && abs(erf_star - erf_star_exact) < 1e-10 &&
    !is.finite(escape_below) && !is.finite(escape_boundary) &&
    is.finite(sup_below) && err_sup_below < 1e-9 && sup_at > 1e6 && blowup_size > 1e6 &&
    min_gap > 0 && err_gap_closed_form < 1e-5 && fibre_min > 0 &&
    err_zero == 0 && err_same_time == 0

emit("nonautonomous-quadratic-decay", if (ok) "pass" else "fail",
     "The process and pullback limit of x' = -exp(-t^2) x^2 match a Runge-Kutta solution, and the lower fibre escapes in finite time when x0 exceeds 1/sqrt(pi)",
     list(rk4_max_abs_error = err_rk4,
          identity_error = err_id,
          cocycle_max_abs_error = err_cocycle,
          pullback_limit_max_abs_error = err_pullback,
          entire_solution_max_residual = resid_xi,
          forward_limit_max_abs_error = err_forward,
          threshold = threshold,
          escape_time_x0_1 = t_star_above,
          escape_erf_abs_error = abs(erf_star - erf_star_exact),
          blowup_size_near_escape = blowup_size,
          sup_lower_fibre_x0_0.5 = sup_below,
          sup_lower_fibre_closed_form_error = err_sup_below,
          sup_lower_fibre_at_threshold = sup_at,
          min_separation_of_solutions = min_gap,
          separation_closed_form_error = err_gap_closed_form,
          fibre_length_min = fibre_min,
          fibre_length_max = fibre_max,
          zero_solution_error = err_zero,
          same_time_error = err_same_time))

# ============================ #
# Figure ####
# ============================ #
#
# The entire solutions of the equation, which are the fibres of the pullback attractor, over a
# window of time. Bounded entire solutions exist only for c in an interval; the two extreme ones
# bound the attractor, and a solution started above the threshold leaves in finite time.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    tt <- seq(-4, 4, by = 0.01)
    cs <- c(-0.5, -0.3, -0.15, 0.15, 0.3, 0.6)
    fam <- do.call(rbind, lapply(cs, function(c)
        data.frame(t = tt, x = xi(tt, c), c = sprintf("%+.2f", c))))
    upper <- data.frame(t = tt, x = xi(tt, -threshold))
    escaping <- data.frame(t = seq(-3, escape_time(1) - 1e-3, by = 0.005))
    escaping$x <- vapply(escaping$t, function(s) Phi(s, -3, xi(-3, -1)), numeric(1))
    p1 <- ggplot(fam, aes(t, x, colour = c)) +
        geom_line(linewidth = 0.5) +
        geom_line(data = upper, aes(t, x), inherit.aes = FALSE,
                  colour = "#B8390E", linewidth = 0.8, linetype = "22") +
        annotate("text", x = -3.9, y = -4.2, hjust = 0, size = 2.6, colour = "#B8390E",
                 label = kb_unicode("Dashed, at $c = -1/\\sqrt{\\pi}$: the boundary case,\nwhich is unbounded in forward time")) +
        scale_colour_viridis_d(option = "viridis", end = 0.9, name = "Parameter c") +
        coord_cartesian(ylim = c(-5, 1.2)) +
        labs(title = "Entire solutions, the fibres of the pullback attractor",
             subtitle = kb_unicode("Solutions of $x' = -e^{-t^2}x^2$, bounded for $c$ above $-1/\\sqrt{\\pi}$"),
             x = kb_tex("Time $t$"), y = kb_tex("State $x$")) +
        kb_theme()
    p2 <- ggplot(escaping, aes(t, x)) +
        geom_line(colour = "#440154", linewidth = 0.6) +
        geom_vline(xintercept = t_star_above, colour = "#B8390E", linewidth = 0.4, linetype = "22") +
        annotate("text", x = t_star_above - 0.15, y = 6, hjust = 1, size = 2.6, colour = "#B8390E",
                 label = kb_unicode(sprintf("Escape at $t = %.4f$", t_star_above))) +
        coord_cartesian(ylim = c(-8, 8)) +
        labs(title = "Below them, a solution leaves in finite time",
             subtitle = kb_unicode("Started at $t = -3$ from a value below the attractor"),
             x = kb_tex("Time $t$"), y = kb_tex("State $x$")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.2, 1)) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "The damping $e^{-t^2}$ has a finite integral, so the equation contracts only by a finite factor and the pullback attractor is a band of entire solutions and not a single one. In the run recorded by checks/nonautonomous-quadratic-decay.R the process agreed with a Runge-Kutta solution to %.1e, the entire solutions left a residual of %.1e in the equation, and a solution starting below the threshold $-1/\\sqrt{\\pi} = %.4f$ escaped at $t = %.4f$.",
                err_rk4, resid_xi, -threshold, t_star_above)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "nonautonomous-quadratic-decay", width = 9.4, height = 4.2))
}
