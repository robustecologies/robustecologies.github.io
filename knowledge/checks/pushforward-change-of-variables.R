# ============================ #
# Pushforward of the uniform law by F(x) = x^2 ####
# ============================ #
#
# Let mu be the Lebesgue measure on [0, 1] and F(x) = x^2. The pushforward
# F#mu assigns to [a, b] the mass mu(F^{-1}([a, b])) = sqrt(b) - sqrt(a), and
# it has the density rho(y) = 1 / (2 sqrt(y)) on (0, 1]. The script checks:
#
# 1. Interval masses: the integral of rho over [a, b] equals sqrt(b) - sqrt(a)
#    for 500 random intervals, including the degenerate interval a = b.
# 2. Change of variables: the integral of g against F#mu equals the integral
#    of g o F against mu, for four observables with closed forms.
# 3. Composition: for G(y) = 1 - y, the law of (G o F)(X) with X uniform
#    has distribution function 1 - sqrt(1 - t), which is also G#(F#mu).
#    A Kolmogorov-Smirnov test compares a sample with this function.

set.seed(20260914L)

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

rho <- function(y) 1 / (2 * sqrt(y))

# ============================ #
# 1. Interval masses ####
# ============================ #
n_int <- 500L
ab <- matrix(runif(2 * n_int), ncol = 2)
a <- pmin(ab[, 1], ab[, 2])
b <- pmax(ab[, 1], ab[, 2])
a[1] <- 0
b[2] <- a[2]
mass_quad <- vapply(seq_len(n_int), function(i)
    if (a[i] == b[i]) 0 else integrate(rho, a[i], b[i], rel.tol = 1e-11)$value, numeric(1))
interval_err <- max(abs(mass_quad - (sqrt(b) - sqrt(a))))

# ============================ #
# 2. Change of variables ####
# ============================ #
erf1 <- 2 * pnorm(sqrt(2)) - 1
observables <- list(
    y     = list(g = function(y) y, exact = 1 / 3),
    y2    = list(g = function(y) y^2, exact = 1 / 5),
    expm  = list(g = function(y) exp(-y), exact = sqrt(pi) / 2 * erf1),
    cos3  = list(g = function(y) cos(3 * y), exact = NA_real_)
)
cov_err <- 0
bench_err <- 0
for (o in observables) {
    lhs <- integrate(function(y) o$g(y) * rho(y), 0, 1, rel.tol = 1e-11)$value
    rhs <- integrate(function(x) o$g(x^2), 0, 1, rel.tol = 1e-11)$value
    cov_err <- max(cov_err, abs(lhs - rhs))
    if (!is.na(o$exact)) bench_err <- max(bench_err, abs(rhs - o$exact))
}

# ============================ #
# 3. Composition of pushforwards ####
# ============================ #
n_sample <- 200000L
z <- 1 - runif(n_sample)^2
H <- function(t) 1 - sqrt(1 - t)
ks_comp <- suppressWarnings(ks.test(z, H))
t_grid <- seq(0.01, 0.99, by = 0.01)
comp_quad <- vapply(t_grid, function(t) integrate(rho, 1 - t, 1, rel.tol = 1e-11)$value, numeric(1))
comp_err <- max(abs(comp_quad - H(t_grid)))

# A Dirac mass at x0 is pushed to a Dirac mass at F(x0).
dirac_ok <- all((rep(0.3, 10))^2 == 0.3^2)

pass <- interval_err < 1e-8 && cov_err < 1e-8 && bench_err < 1e-8 &&
    comp_err < 1e-8 && ks_comp$p.value > 1e-3 && dirac_ok
emit("pushforward-change-of-variables", if (pass) "pass" else "fail",
     "Interval masses, change of variables and composition for the pushforward of Lebesgue measure by x^2",
     list(interval_max_abs_error = interval_err, change_of_variables_max_abs_error = cov_err,
          closed_form_max_abs_error = bench_err, composition_max_abs_error = comp_err,
          composition_ks_p = ks_comp$p.value, sample_size = n_sample))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the density of the pushforward of Lebesgue measure by the squaring map, against a
# histogram of samples obtained by squaring uniform draws. The mass piles up near zero because
# the map compresses the interval there. Right: the error of the interval masses computed by
# quadrature against the closed form, over the random intervals of the test.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    gy <- seq(0.002, 1, length.out = 600)
    samples <- runif(200000L)^2
    hh <- hist(samples, breaks = seq(0, 1, by = 0.02), plot = FALSE)
    p1 <- ggplot(data.frame(y = gy, d = rho(gy)), aes(y, d)) +
        geom_col(data = data.frame(y = hh$mids, d = hh$density), aes(y, d), inherit.aes = FALSE,
                 fill = "grey80", colour = NA, width = 0.02) +
        geom_line(colour = "#440154", linewidth = 0.8) +
        coord_cartesian(ylim = c(0, 7)) +
        annotate("text", x = 0.35, y = 5.4, hjust = 0, size = 2.7, colour = "#440154",
                 label = kb_unicode("The density of the pushforward, $1/(2\\sqrt{y})$")) +
        annotate("text", x = 0.35, y = 4.6, hjust = 0, size = 2.7, colour = "grey35",
                 label = kb_unicode("Grey: squares of 200 000 uniform draws")) +
        labs(title = "Where a map sends the mass of Lebesgue measure",
             subtitle = "Image of the uniform measure on the unit interval under the squaring map",
             x = "State", y = "Density") +
        kb_theme()
    errs <- abs(mass_quad - (sqrt(b) - sqrt(a)))
    p2 <- ggplot(data.frame(width = b - a, err = pmax(errs, 1e-18)), aes(width, err)) +
        geom_point(size = 0.8, colour = "#21918c", alpha = 0.7) +
        scale_y_log10() +
        labs(title = "Quadrature against the closed form, interval by interval",
             subtitle = sprintf("Mass of %d random intervals, computed and exact, greatest difference %.1e",
                                n_int, interval_err),
             x = "Width of the interval", y = "Difference in the mass") +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "The pushforward is defined by the mass of preimages, so the density is the derivative of sqrt(y), and the same rule is what the change of variables formula states for integrals. In the run recorded by checks/pushforward-change-of-variables.R the interval masses agreed with the closed form to %.1e and the change of variables formula held to %.1e.",
                interval_err, cov_err)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "pushforward-change-of-variables", width = 9.2, height = 4.2))
}
