# ============================ #
# Duality of the Koopman and Perron-Frobenius operators ####
# ============================ #
#
# For a nonsingular map S on a measure space (X, m), the Perron-Frobenius
# operator P on L1(m) and the Koopman operator K g = g o S on L-infinity(m)
# satisfy the duality <P f, g> = <f, K g>, where <f, g> is the integral of
# f g with respect to m. The script checks the identity for the logistic map
# S(x) = 4 x (1 - x) on [0, 1] with m the Lebesgue measure, using the explicit
# form (P f)(x) = [f(xi_minus) + f(xi_plus)] / (4 sqrt(1 - x)).
#
# Sixteen pairs (f, g) are integrated by adaptive quadrature. Four pairs have
# closed forms that serve as an independent benchmark:
#   f = 1,  g = x   : 2/3
#   f = 1,  g = x^2 : 8/15
#   f = f*, g = x   : 1/2   (f* the invariant density, by symmetry)
#   f = f*, g = x^2 : 3/8   (second moment of the arcsine law on [0, 1])
# The script also checks that P preserves the integral of a density.

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

S        <- function(x) 4 * x * (1 - x)
xi_minus <- function(y) -expm1(0.5 * log1p(-y)) / 2
xi_plus  <- function(y) 1 - xi_minus(y)
P_op     <- function(f) function(x) (f(xi_minus(x)) + f(xi_plus(x))) / (4 * sqrt(1 - x))
K_op     <- function(g) function(x) g(S(x))

# Quadrature on [0, 1] split at the discontinuities of the integrand and at
# 1/2, where the two monotone branches of S meet.
integrate01 <- function(h, breaks = numeric(0)) {
    b <- sort(unique(c(0, 0.5, breaks[breaks > 0 & breaks < 1], 1)))
    pieces <- vapply(seq_len(length(b) - 1), function(k)
        integrate(h, b[k], b[k + 1], subdivisions = 2000L, rel.tol = 1e-11)$value,
        numeric(1))
    sum(pieces)
}

densities <- list(
    one    = function(x) rep(1, length(x)),
    linear = function(x) 2 * x,
    square = function(x) 3 * x^2,
    fstar  = function(x) 1 / (pi * sqrt(x * (1 - x)))
)
observables <- list(
    x         = list(g = function(x) x, jumps = numeric(0)),
    x2        = list(g = function(x) x^2, jumps = numeric(0)),
    cos5      = list(g = function(x) cos(5 * x), jumps = numeric(0)),
    indicator = list(g = function(x) as.numeric(x >= 0.1 & x <= 0.4), jumps = c(0.1, 0.4))
)

pairs <- expand.grid(f = names(densities), g = names(observables), stringsAsFactors = FALSE)
pairs$lhs <- NA_real_
pairs$rhs <- NA_real_
for (k in seq_len(nrow(pairs))) {
    f <- densities[[pairs$f[k]]]
    obs <- observables[[pairs$g[k]]]
    lhs_breaks <- obs$jumps
    rhs_breaks <- c(xi_minus(obs$jumps), xi_plus(obs$jumps))
    pairs$lhs[k] <- integrate01(function(x) P_op(f)(x) * obs$g(x), lhs_breaks)
    pairs$rhs[k] <- integrate01(function(x) f(x) * K_op(obs$g)(x), rhs_breaks)
}
duality_err <- max(abs(pairs$lhs - pairs$rhs))

value_of <- function(fn, gn) pairs$lhs[pairs$f == fn & pairs$g == gn]
benchmark_err <- max(abs(c(
    value_of("one", "x") - 2 / 3,
    value_of("one", "x2") - 8 / 15,
    value_of("fstar", "x") - 1 / 2,
    value_of("fstar", "x2") - 3 / 8
)))

mass_err <- max(vapply(densities, function(f)
    abs(integrate01(P_op(f)) - integrate01(f)), numeric(1)))

# The operator is not injective. The two preimages of x satisfy xi_plus = 1 - xi_minus, so
# reflecting a density about 1/2 exchanges them and leaves P f unchanged: P f = P g with
# g(u) = f(1 - u). The test evaluates both sides on a grid for three asymmetric densities
# and confirms that f and g differ, so the identity is not trivial. The comparison is
# pointwise and exact up to rounding, hence the tolerance 1e-12.
reflect <- function(f) function(u) f(1 - u)
grid <- seq(0.001, 0.999, length.out = 999)
asymmetric <- list(linear = densities$linear, square = densities$square,
                   shifted = function(u) 1 + 0.5 * cos(3 * u))
kernel_err <- 0
separation <- Inf
for (f in asymmetric) {
    g <- reflect(f)
    kernel_err <- max(kernel_err, max(abs(P_op(f)(grid) - P_op(g)(grid))))
    separation <- min(separation, max(abs(f(grid) - g(grid))))
}

pass <- duality_err < 1e-7 && benchmark_err < 1e-7 && mass_err < 1e-7 &&
    kernel_err < 1e-12 && separation > 0.1
emit("koopman-pf-duality", if (pass) "pass" else "fail",
     "<P f, g> = <f, g o S> for 16 pairs on the logistic map, with closed-form benchmarks, and P is not injective",
     list(pairs = nrow(pairs), duality_max_abs_error = duality_err,
          benchmark_max_abs_error = benchmark_err, mass_max_abs_error = mass_err,
          reflected_density_max_abs_error = kernel_err,
          smallest_separation_of_the_pair = separation))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the two sides of the duality for the sixteen pairs of the test, against the line where
# they are equal. Right: why the transfer operator cannot be inverted. A density and its
# reflection about one half have the same image, because the two preimages of a point under the
# logistic map add to one, so the operator loses the difference between them.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    pairs$label <- paste(pairs$f, pairs$g, sep = " with ")
    p1 <- ggplot(pairs, aes(lhs, rhs, colour = f, shape = g)) +
        geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
        geom_point(size = 2) +
        scale_colour_viridis_d(option = "viridis", end = 0.85, name = "Density",
                               labels = kb_tex(c("$f^{*}$", "$2x$", "$1$", "$3x^2$"))) +
        scale_shape_manual(values = c(16, 17, 15, 18), name = "Observable",
                           labels = kb_tex(c("$\\cos 5x$", "$1[0.1, 0.4]$", "$x$", "$x^2$"))) +
        labs(title = "The two sides of the duality agree on every pair",
             subtitle = sprintf("Four densities against four observables, greatest difference %.1e", duality_err),
             x = kb_tex("$\\int (\\mathcal{P}f)\\, g \\, dm$"),
             y = kb_tex("$\\int f\\, (\\mathcal{K}g) \\, dm$")) +
        kb_theme() +
        theme(legend.box = "vertical", legend.spacing.y = unit(1, "pt")) +
        guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2))
    gx <- seq(0.001, 0.999, length.out = 600)
    f_lin <- densities$linear
    f_ref <- reflect(f_lin)
    dd <- rbind(
        data.frame(x = gx, y = f_lin(gx), what = "A density, 2x"),
        data.frame(x = gx, y = f_ref(gx), what = "Its reflection, 2(1 - x)"),
        data.frame(x = gx, y = P_op(f_lin)(gx), what = "The image of both under the operator"))
    p2 <- ggplot(dd, aes(x, y, colour = what, linetype = what)) +
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = c("A density, 2x" = "#21918c",
                                       "Its reflection, 2(1 - x)" = "#B8390E",
                                       "The image of both under the operator" = "#440154"),
                            name = NULL) +
        scale_linetype_manual(values = c("solid", "solid", "22"), name = NULL) +
        coord_cartesian(ylim = c(0, 2.6)) +
        labs(title = "The transfer operator is not injective",
             subtitle = sprintf("Both densities have the same image, here to %.1e", kernel_err),
             x = kb_tex("State $x$"), y = kb_tex("Density $f(x)$")) +
        kb_theme() + theme(legend.box = "vertical") +
        guides(colour = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3))
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "The duality is exact and the failure of injectivity is its algebraic shadow: the composition operator is not onto, and the pair is the operator form of irreversibility for a map that folds. In the run recorded by checks/koopman-pf-duality.R the two sides agreed to %.1e over the sixteen pairs, the two images agreed to %.1e, and the two densities themselves differed by %.2f at their farthest point.",
                duality_err, kernel_err, separation)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "koopman-pf-duality", width = 9.6, height = 4.6))
}
