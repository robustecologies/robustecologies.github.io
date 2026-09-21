# ============================ #
# The asymptotically stable set and the measure attractor of Milnor's square map ####
# ============================ #
#
# Claims, with the locators of milnor1985a (section 1, the example after the correction):
#
# For f(x, y) = (1 - x)(x, y) on the unit square, the edge x = 0 is the only compact
# invariant set that is asymptotically stable, while the omega-limit set of almost every
# point is the origin, so the origin is the unique attractor in the sense of milnor1985.
#
# 1. Every point of the edge is fixed, and the edge is invariant.
# 2. For x0 in (0, 1) the first coordinate obeys x_{n+1} = x_n (1 - x_n), so n x_n tends
#    to 1. Reference: the closed form of the asymptotics of that recurrence, which is
#    independent of the iteration performed here.
# 3. The map multiplies both coordinates by the same factor 1 - x, so every ray through the
#    origin is invariant and y_n / x_n = y_0 / x_0 for every n. The product that gives y_n
#    telescopes, y_n = y_0 prod_{k < n} (1 - x_k) = y_0 x_n / x_0, so n y_n tends to
#    y_0 / x_0. Both the exact ratio and that limit are tested.
# 4. The orbit therefore converges to the origin from every point with x0 in (0, 1], so
#    the omega-limit set is the origin for every point of the square outside the edge,
#    which is a set of full Lebesgue measure.
# 5. The edge attracts a neighbourhood uniformly, since the distance to it is x_n, which
#    decreases monotonically; the origin does not, since points of the edge never move.
# 6. Degenerate cases: a point of the edge is fixed, and x0 = 1 lands on the origin in one
#    step.

set.seed(20260916L)

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

step <- function(v) c((1 - v[1]) * v[1], (1 - v[1]) * v[2])

orbit <- function(v0, n) {
    v <- v0
    xs <- numeric(n); ys <- numeric(n)
    for (i in seq_len(n)) { v <- step(v); xs[i] <- v[1]; ys[i] <- v[2] }
    list(x = xs, y = ys)
}

# ============================ #
# Test 1: the edge is fixed ####
# ============================ #

edge <- seq(0, 1, by = 0.05)
err_edge <- max(abs(vapply(edge, function(y) max(abs(step(c(0, y)) - c(0, y))), numeric(1))))

# ============================ #
# Test 2: the first coordinate decays like 1/n ####
# ============================ #
#
# The recurrence x_{n+1} = x_n (1 - x_n) has n x_n -> 1. The convergence is logarithmically
# slow, since n x_n = 1 + log(n) / n + O(1/n), so the tolerance at n = 2e5 is 1e-3 on the
# product and the residual of the asymptotic form is reported as well.

n <- 200000L
o <- orbit(c(0.7, 0.9), n)
prod_x <- n * o$x[n]
# The expansion of the recurrence is n x_n = 1 - log(n)/n + O(1/n), so at n = 2e5 the
# residual is of order 1/n = 5e-6 and the tolerance is 5e-5.
resid_x <- abs(prod_x - (1 - log(n) / n))

# ============================ #
# Test 3: the second coordinate decays like 1/n ####
# ============================ #

prod_y <- n * o$y[n]
err_limit_y <- abs(prod_y - (0.9 / 0.7) * prod_x)        # exact, since y_n = y_0 x_n / x_0
# The ratio y_n / x_n is invariant, to rounding, along every orbit.
err_ratio <- 0
for (v0 in list(c(0.7, 0.9), c(0.05, 0.2), c(0.9, 0.05))) {
    oo <- orbit(v0, 5000L)
    err_ratio <- max(err_ratio, max(abs(oo$y / oo$x - v0[2] / v0[1])))
}

# ============================ #
# Test 4: convergence to the origin from a grid of the square ####
# ============================ #

# The distance after n steps is x_n sqrt(1 + (y_0/x_0)^2), and n x_n is close to 1, so the
# predicted distance is sqrt(x_0^2 + y_0^2) / (x_0 n) and the test compares with it.
grid <- expand.grid(x = seq(0.02, 1, length.out = 12), y = seq(0, 1, length.out = 12))
steps <- 20000L
far <- 0
err_predicted <- 0
for (i in seq_len(nrow(grid))) {
    v0 <- c(grid$x[i], grid$y[i])
    v <- v0
    for (k in seq_len(steps)) v <- step(v)
    d <- sqrt(sum(v^2))
    far <- max(far, d)
    # x0 = 1 lands on the origin in one step and is the degenerate case of test 6, so the
    # asymptotic law is compared only where the orbit is infinite.
    if (v0[1] < 1) {
        predicted <- sqrt(sum(v0^2)) / (v0[1] * steps)
        err_predicted <- max(err_predicted, abs(d / predicted - 1))
    }
}

# ============================ #
# Test 5: the edge attracts, the origin does not attract the edge ####
# ============================ #
#
# The distance to the edge is the first coordinate, which decreases monotonically, so a
# neighbourhood of the edge is drawn in. A point of the edge stays where it is, so its
# distance to the origin never decreases.

o2 <- orbit(c(0.9, 0.4), 500L)
edge_monotone <- all(diff(o2$x) < 0)
edge_limit <- o2$x[500]
origin_dist_edge <- max(abs(vapply(c(0.3, 0.6, 1), function(y) {
    v <- c(0, y); for (k in seq_len(1000L)) v <- step(v); sqrt(sum(v^2)) - y
}, numeric(1))))

# ============================ #
# Test 6: degenerate cases ####
# ============================ #

err_fixed <- max(abs(step(c(0, 0.37)) - c(0, 0.37)))
err_one <- max(abs(step(c(1, 0.8)) - c(0, 0)))

ok <- err_edge == 0 && abs(prod_x - 1) < 1e-3 && resid_x < 5e-5 &&
    err_limit_y < 1e-12 && err_ratio < 1e-12 &&
    far < 0.01 && err_predicted < 0.01 &&
    edge_monotone && edge_limit < 0.01 && origin_dist_edge == 0 &&
    err_fixed == 0 && err_one == 0

emit("milnor-square-example", if (ok) "pass" else "fail",
     "Almost every orbit of Milnor's square map converges to the origin while the asymptotically stable set is the edge, so the two notions of attractor select different sets",
     list(edge_fixed_error = err_edge,
          n_times_x_at_2e5 = prod_x,
          x_asymptotic_residual = resid_x,
          n_times_y_at_2e5 = prod_y,
          y_limit_closed_form_error = err_limit_y,
          ray_ratio_max_error = err_ratio,
          max_distance_to_origin_after_20000_steps = far,
          distance_closed_form_max_rel_error = err_predicted,
          edge_distance_monotone = edge_monotone,
          edge_distance_after_500_steps = edge_limit,
          origin_distance_change_on_edge = origin_dist_edge,
          fixed_point_error = err_fixed,
          boundary_x_equals_one_error = err_one))

# ============================ #
# Figure ####
# ============================ #
#
# Left: orbits of the map in the square. Every orbit moves along the ray through its starting
# point, since the map scales both coordinates by the same factor, and reaches the origin,
# while the edge x = 0 stands still. Right: the decay of the first coordinate against the
# closed form 1/n, on logarithmic axes.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    starts <- expand.grid(x = c(0.25, 0.5, 0.75, 1), y = c(0.1, 0.35, 0.6, 0.85, 1))
    n_draw <- 60L
    paths <- do.call(rbind, lapply(seq_len(nrow(starts)), function(i) {
        o <- orbit(c(starts$x[i], starts$y[i]), n_draw)
        data.frame(x = c(starts$x[i], o$x), y = c(starts$y[i], o$y), id = i)
    }))
    p1 <- ggplot(paths, aes(x, y, group = id)) +
        geom_path(colour = "#21918c", linewidth = 0.3, alpha = 0.85) +
        geom_point(data = paths[!duplicated(paths$id), ], aes(x, y), inherit.aes = FALSE,
                   colour = "#21918c", size = 0.7) +
        geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), colour = "#B8390E", linewidth = 1.1) +
        annotate("point", x = 0, y = 0, colour = "black", size = 1.6) +
        annotate("text", x = 0.06, y = 0.66, hjust = 0, size = 2.5, colour = "#B8390E",
                 label = kb_unicode("The edge $x = 0$ is fixed pointwise\nand is asymptotically stable")) +
        annotate("text", x = 0.12, y = 0.04, hjust = 0, size = 2.5, colour = "black",
                 label = "Almost every orbit converges to the origin") +
        coord_equal(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02)) +
        labs(title = "Two notions of attractor, selecting different sets",
             subtitle = kb_unicode(sprintf("Orbits of $(x, y)$ mapped to $((1 - x)x, (1 - x)y)$ in the unit square, %d steps each", n_draw)),
             x = "First coordinate", y = "Second coordinate") +
        kb_theme()
    ns <- unique(round(10^seq(0, log10(n), length.out = 160)))
    decay <- data.frame(n = ns, x = o$x[ns])
    p2 <- ggplot(decay, aes(n, x)) +
        geom_line(colour = "#440154", linewidth = 0.6) +
        geom_line(aes(n, 1 / n), colour = "#B8390E", linewidth = 0.4, linetype = "22") +
        annotate("text", x = 10, y = 4e-4, hjust = 0, size = 2.6, colour = "#B8390E",
                 label = kb_unicode("Closed form $1/n$")) +
        scale_x_log10() + scale_y_log10() +
        labs(title = "The first coordinate decays like the reciprocal of the step",
             subtitle = kb_unicode(sprintf("Orbit of $(%.1f, %.1f)$, with $n x_n = %.6f$ after %s steps",
                                           0.7, 0.9, prod_x, format(n, big.mark = " "))),
             x = kb_tex("Step $n$"), y = kb_tex("First coordinate $x_n$")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1, 1)) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "The likely limit set is the origin and the asymptotically stable set is the edge, so a measure attractor and an asymptotically stable attractor need not agree. In the run recorded by checks/milnor-square-example.R every point of a grid of %d initial conditions came within %.1e of the origin after %s steps, the edge did not move at all, and n x_n equalled %.6f.",
                nrow(grid), far, format(steps, big.mark = " "), prod_x)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "milnor-square-example", width = 9.2, height = 4.4))
}
