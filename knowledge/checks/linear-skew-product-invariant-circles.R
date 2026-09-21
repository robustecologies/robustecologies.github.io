# ============================ #
# Invariant circles of the linear skew products on the two-torus ####
# ============================ #
#
# Claim. The linear skew product F_{i,j}(y, x) = (i y, x + j y) mod 1 on T^2 preserves the
# function Theta(y, x) = ((i - 1) x - j y) mod 1, so every level set of Theta is invariant and
# F_{i,j} is not topologically transitive. Homburg (2012), section 4, states the invariant family
# in the transposed form "circles parallel to j x = (i - 1) y", which agrees with the above only
# when j = i - 1.
#
# Tests.
#  1. Exact identity. For six pairs (i, j), the circular distance between Theta(F(y, x)) and
#     Theta(y, x) is measured on 2000 random points. The reference is the algebraic identity
#     (i - 1)(x + j y) - j (i y) = (i - 1) x - j y, which the script does not use: it applies the
#     map numerically and evaluates Theta from its definition, so the two sides are computed by
#     different routes and agreement is evidence and not a restatement.
#  2. The transposed form. The same measurement for Phi(y, x) = (j x - (i - 1) y) mod 1, which
#     must drift for the pairs with j != i - 1 and be preserved for the pairs with j = i - 1.
#  3. Orbit test. A single orbit of 200 steps is iterated and the drift of Theta along it is
#     recorded, which tests invariance under composition and not only under one step.
#  4. Degenerate case j = 0. Then F_{i,0}(y, x) = (i y, x) fixes every horizontal circle, so the
#     fibre coordinate is preserved exactly; this is the unperturbed map of Theorem 4.1.
#
# Tolerances. The maps and the invariants are sums and products of integers with numbers in
# [0, 1), so the identity holds in binary floating point up to the rounding of the modulo
# operation: the tolerance for a preserved quantity is 1e-12, far above the observed values and
# far below the drift of a quantity that is not preserved, which is of order 1e-1.

set.seed(20260918L)



emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

F_ij <- function(y, x, i, j) list(y = (i * y) %% 1, x = (x + j * y) %% 1)
theta     <- function(y, x, i, j) ((i - 1) * x - j * y) %% 1
theta_alt <- function(y, x, i, j) (j * x - (i - 1) * y) %% 1
circ_dist <- function(a, b) { d <- abs(a - b) %% 1; pmin(d, 1 - d) }

pairs <- list(c(2, 1), c(3, 1), c(3, 2), c(4, 3), c(3, 0), c(5, 2))
n_pt <- 2000L
rows <- lapply(pairs, function(p) {
    i <- p[1]; j <- p[2]
    y <- runif(n_pt); x <- runif(n_pt)
    z <- F_ij(y, x, i, j)
    data.frame(i = i, j = j,
               drift_theta = max(circ_dist(theta(z$y, z$x, i, j), theta(y, x, i, j))),
               drift_alt   = max(circ_dist(theta_alt(z$y, z$x, i, j), theta_alt(y, x, i, j))),
               agree = j == i - 1)
})
tab <- do.call(rbind, rows)

# Test 3: one orbit of 200 steps for (i, j) = (3, 1).
y <- 0.2718281828; x <- 0.1414213562; th0 <- theta(y, x, 3, 1); drift_orbit <- 0
for (s in seq_len(200L)) {
    z <- F_ij(y, x, 3, 1); y <- z$y; x <- z$x
    drift_orbit <- max(drift_orbit, circ_dist(theta(y, x, 3, 1), th0))
}

# Test 4: degenerate case j = 0.
y0 <- runif(n_pt); x0 <- runif(n_pt)
z0 <- F_ij(y0, x0, 3L, 0L)
drift_fibre <- max(circ_dist(z0$x, x0))

tol <- 1e-12
ok_theta  <- max(tab$drift_theta) < tol
ok_alt    <- all(tab$drift_alt[tab$agree] < tol) && all(tab$drift_alt[!tab$agree] > 0.1)
ok_orbit  <- drift_orbit < tol
ok_fibre  <- drift_fibre < tol
status <- if (ok_theta && ok_alt && ok_orbit && ok_fibre) "pass" else "fail"

# Figure: a single level set of each candidate invariant for (i, j) = (3, 1), together with its
# image under one step of the map. The level set and its image coincide when the function is
# invariant. Both are computed with the same functions the tests above use.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    i <- 3L; j <- 1L; c_lev <- 0.5
    yy <- seq(0, 1, length.out = 200)
    # {(i-1)x - jy = c} has i-1 = 2 branches; {jx - (i-1)y = c} has j = 1 branch.
    set_theta <- data.frame(y = rep(yy, 2),
                            x = c((c_lev + j * yy) / (i - 1), (c_lev + j * yy) / (i - 1) + 0.5) %% 1,
                            panel = "(i - 1)x - jy")
    set_alt <- data.frame(y = yy, x = ((c_lev + (i - 1) * yy) / j) %% 1,
                          panel = "jx - (i - 1)y")
    lev <- rbind(set_theta, set_alt)
    im <- F_ij(lev$y, lev$x, i, j)
    img <- data.frame(y = im$y, x = im$x, panel = lev$panel)
    lev$role <- "Level set"; img$role <- "Its image under one step"
    both <- rbind(lev, img)

    both$role <- factor(both$role, levels = c("Level set", "Its image under one step"))
    p <- ggplot(both, aes(y, x, colour = role, size = role)) +
        geom_point(shape = 16, alpha = 0.9) +
        facet_wrap(~ panel) +
        scale_colour_manual(name = NULL, values = c("Level set" = "#056796",
                                                    "Its image under one step" = "#be1117")) +
        scale_size_manual(name = NULL, values = c("Level set" = 2.1, "Its image under one step" = 0.7)) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        guides(colour = guide_legend(override.aes = list(size = 2.2))) +
        labs(title = "Invariant circles of a linear skew product",
             subtitle = kb_unicode("Level set at $c = 0.5$ of two candidate invariants of $F(y, x) = (3y, x + y) \\bmod 1$, with its image after one step"),
             x = kb_tex("Base coordinate $y$"), y = kb_tex("Fibre coordinate $x$"),
             caption = kb_caption(paste(
                 "Left: the image lies on the level set of $(i - 1)x - jy$, the small red points sitting inside the wider blue ones, so every level set is invariant and the map is not topologically transitive.",
                 "Right: the level set of the transposed form $jx - (i - 1)y$ and its image are different sets, so that function is not invariant for $(i, j) = (3, 1)$.",
                 "Drawn by checks/linear-skew-product-invariant-circles.R from the functions the numerical tests use."))) +
        kb_theme() + theme(panel.spacing.x = unit(1.4, "lines"))
    kb_save(p, "linear-skew-product-invariant-circles")
}

emit("linear-skew-product-invariant-circles", status,
     "(i-1)x - jy is invariant under F_{i,j}, and the transposed form jx - (i-1)y is not unless j = i-1",
     list(max_drift_theta = max(tab$drift_theta),
          max_drift_alt_when_j_eq_i_minus_1 = max(tab$drift_alt[tab$agree]),
          min_drift_alt_otherwise = min(tab$drift_alt[!tab$agree]),
          drift_along_orbit_200_steps = drift_orbit,
          drift_fibre_j_zero = drift_fibre,
          n_pairs = nrow(tab)))
