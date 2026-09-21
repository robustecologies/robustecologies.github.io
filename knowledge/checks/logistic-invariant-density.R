# ============================ #
# Invariant density of the logistic map at r = 4 ####
# ============================ #
#
# The map is S(x) = 4 x (1 - x) on [0, 1]. The candidate invariant density is
# f*(x) = 1 / (pi sqrt(x (1 - x))), with distribution function
# F*(x) = (2 / pi) asin(sqrt(x)). The script tests the claim in three
# independent ways and prints one JSON line for tools/kb.py verify.
#
# 1. Fixed point: (P f*)(x) = f*(x) on a grid, where the Perron-Frobenius
#    operator of S is (P f)(x) = [f(xi_minus) + f(xi_plus)] / (4 sqrt(1 - x))
#    and xi_minus, xi_plus = (1 -/+ sqrt(1 - x)) / 2 are the two preimages of x.
# 2. Pushforward of a sample: if X has distribution F*, then S(X) has
#    distribution F*. A one-sample Kolmogorov-Smirnov test compares S(X) with
#    F*; the same test against the uniform law shows that the test has power.
# 3. Ulam's method: the stationary vector of the exact Ulam matrix on n equal
#    bins is compared in l1 with the exact bin masses of F*.

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

S      <- function(x) 4 * x * (1 - x)
f_star <- function(x) 1 / (pi * sqrt(x * (1 - x)))
F_star <- function(x) (2 / pi) * asin(sqrt(x))

# Preimages of y under S. The form -expm1(log1p(-y) / 2) = 1 - sqrt(1 - y)
# avoids the cancellation of the naive expression for small y.
xi_minus <- function(y) -expm1(0.5 * log1p(-y)) / 2
xi_plus  <- function(y) 1 - xi_minus(y)

perron_frobenius <- function(f, x) (f(xi_minus(x)) + f(xi_plus(x))) / (4 * sqrt(1 - x))

# ============================ #
# 1. Fixed point of the Perron-Frobenius operator ####
# ============================ #
grid <- seq(1e-4, 1 - 1e-4, length.out = 20001)
fixed_point_rel_err <- max(abs(perron_frobenius(f_star, grid) / f_star(grid) - 1))
normalisation <- integrate(f_star, 0, 0.5, rel.tol = 1e-10)$value +
    integrate(f_star, 0.5, 1, rel.tol = 1e-10)$value
normalisation_err <- abs(normalisation - 1)

# ============================ #
# 2. Pushforward of a sample drawn from F* ####
# ============================ #
# X = sin(pi U / 2)^2 with U uniform has distribution F*. The generator has
# 32-bit resolution, so a sample of 1e6 draws contains a few ties; the
# asymptotic Kolmogorov-Smirnov p-value is unaffected at this size.
n_sample <- 1000000L
u  <- runif(n_sample)
x0 <- sin(pi * u / 2)^2
y  <- S(x0)
ks_sampler <- suppressWarnings(ks.test(x0, F_star))
ks_push    <- suppressWarnings(ks.test(y, F_star))
ks_uniform <- suppressWarnings(ks.test(y, "punif"))

# ============================ #
# 3. Ulam's method with the exact transition matrix ####
# ============================ #
# Row i of the Ulam matrix is p_ij = m(I_i intersect S^{-1}(I_j)) / m(I_i),
# with m the Lebesgue measure. The preimage of I_j = [c, d] is the union of
# [xi_minus(c), xi_minus(d)] and [xi_plus(d), xi_plus(c)], so the entries are
# lengths of interval intersections and need no sampling.
ulam <- function(n) {
    e  <- seq(0, 1, length.out = n + 1)
    lo <- e[-(n + 1)]
    hi <- e[-1]
    # pmax() drops the dim attribute, so negative overlaps are clamped in place.
    overlap <- function(a2, b2) {
        m <- outer(hi, b2, pmin) - outer(lo, a2, pmax)
        m[m < 0] <- 0
        m
    }
    P <- (overlap(xi_minus(lo), xi_minus(hi)) + overlap(xi_plus(hi), xi_plus(lo))) / (hi - lo)
    row_err <- max(abs(rowSums(P) - 1))
    # Stationary row vector: p (P - I) = 0 with sum(p) = 1.
    A <- t(P) - diag(n)
    A[n, ] <- 1
    b <- c(rep(0, n - 1), 1)
    p <- qr.solve(A, b)
    q <- diff(F_star(e))
    # The stationary vector and the partition are returned as well, so that the figure below
    # draws the same approximation whose error the metrics record.
    list(l1 = sum(abs(p - q)), row_err = row_err, min_p = min(p), p = p, edges = e)
}
bins <- c(50L, 100L, 200L, 400L, 800L)
ulam_runs <- lapply(bins, ulam)
l1 <- vapply(ulam_runs, `[[`, numeric(1), "l1")
row_err <- max(vapply(ulam_runs, `[[`, numeric(1), "row_err"))

# ============================ #
# 4. Affine conjugacy with the quadratic family ####
# ============================ #
# For r != 2, x = 1/2 + (r - 2) y / 4 conjugates S_r(x) = r x (1 - x) to
# Q_a(y) = 1 - a y^2 with a = r (r - 2) / 4: S_r(alpha y + 1/2) = alpha Q_a(y) + 1/2.
r_vals <- runif(2000, 2.01, 4)
y_vals <- runif(2000, -1, 1)
alpha <- (r_vals - 2) / 4
a_vals <- r_vals * (r_vals - 2) / 4
conj_lhs <- r_vals * (alpha * y_vals + 0.5) * (1 - (alpha * y_vals + 0.5))
conj_rhs <- alpha * (1 - a_vals * y_vals^2) + 0.5
conjugacy_err <- max(abs(conj_lhs - conj_rhs))

# ============================ #
# Verdict ####
# ============================ #
pass <- conjugacy_err < 1e-12 &&
    fixed_point_rel_err < 1e-9 &&
    normalisation_err < 1e-6 &&
    ks_sampler$p.value > 1e-3 &&
    ks_push$p.value > 1e-3 &&
    ks_uniform$p.value < 1e-10 &&
    row_err < 1e-10 &&
    all(diff(l1) < 0) &&
    l1[length(l1)] < 0.05

metrics <- list(
    fixed_point_max_rel_error = fixed_point_rel_err,
    normalisation_error = normalisation_err,
    ks_sampler_p = ks_sampler$p.value,
    ks_pushforward_D = unname(ks_push$statistic),
    ks_pushforward_p = ks_push$p.value,
    ks_against_uniform_D = unname(ks_uniform$statistic),
    ulam_row_sum_error = row_err,
    ulam_l1_n50 = l1[1], ulam_l1_n100 = l1[2], ulam_l1_n200 = l1[3],
    ulam_l1_n400 = l1[4], ulam_l1_n800 = l1[5],
    sample_size = n_sample,
    conjugacy_max_abs_error = conjugacy_err
)
emit("logistic-invariant-density", if (pass) "pass" else "fail",
     "f* is a fixed point of P, S pushes F* to F*, Ulam approximations converge to F* in l1, and S_r is conjugate to Q_a",
     metrics)

# ============================ #
# Figure ####
# ============================ #
#
# Left: the stationary vector of the Ulam matrix on 200 cells, drawn as a density, against the
# invariant density f*(x) = 1 / (pi sqrt(x (1 - x))) that the check has shown to be a fixed
# point of the operator. Right: the error of that approximation in the l1 norm, over the five
# partitions the check uses, on logarithmic axes.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    k <- which(bins == 200L)
    u <- ulam_runs[[k]]
    mid <- (u$edges[-1] + u$edges[-length(u$edges)]) / 2
    width <- diff(u$edges)
    dens <- data.frame(x = mid, ulam = u$p / width)
    curve_df <- data.frame(x = seq(1e-4, 1 - 1e-4, length.out = 1200))
    curve_df$f <- f_star(curve_df$x)
    top <- 6
    p1 <- ggplot(dens, aes(x, ulam)) +
        geom_col(width = width, fill = "grey78", colour = NA) +
        geom_line(data = curve_df, aes(x, f), inherit.aes = FALSE,
                  colour = "#440154", linewidth = 0.7) +
        coord_cartesian(ylim = c(0, top)) +
        annotate("text", x = 0.5, y = top * 0.82, hjust = 0.5, size = 2.6, colour = "#440154",
                 label = kb_unicode("The invariant density, $1/(\\pi\\sqrt{x(1 - x)})$")) +
        annotate("text", x = 0.5, y = top * 0.68, hjust = 0.5, size = 2.6, colour = "grey35",
                 label = sprintf("Grey: the stationary vector of the Ulam matrix on %d cells", bins[k])) +
        labs(title = "A matrix approximation of the transfer operator",
             subtitle = "The density is unbounded at both ends, so the cells at the ends hold the excess mass",
             x = "State", y = "Density") +
        kb_theme()
    conv <- data.frame(bins = bins, l1 = l1)
    p2 <- ggplot(conv, aes(bins, l1)) +
        geom_line(colour = "#21918c", linewidth = 0.5) +
        geom_point(colour = "#21918c", size = 1.6) +
        geom_text(aes(label = sprintf("%.4f", l1)), vjust = -0.9, size = 2.4, colour = "grey30") +
        scale_x_log10(breaks = bins) + scale_y_log10() +
        labs(title = "The approximation improves as the partition is refined",
             subtitle = kb_unicode("Distance in the $\\ell^1$ norm between the stationary vector and the exact cell masses"),
             x = kb_tex("Number of cells $N$"), y = kb_tex("Error in the $\\ell^1$ norm")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "Both panels come from the run recorded by checks/logistic-invariant-density.R, in which the arcsine density was a fixed point of the operator to %.1e in relative terms, the rows of every Ulam matrix summed to one to %.1e, and the $\\ell^1$ error fell from %.4f on %d cells to %.4f on %d cells.",
                fixed_point_rel_err, row_err, l1[1], bins[1], l1[length(l1)], bins[length(bins)])),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "logistic-invariant-density", width = 9.2, height = 4.2))
}
