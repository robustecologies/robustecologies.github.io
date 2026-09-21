# ============================ #
# The recurrent sets of three circle maps and the inclusions between them ####
# ============================ #
#
# Claim. For a continuous map of a compact metric space the recurrent sets are nested,
# L(T) inside Omega(T) inside R(T), where L is the closure of the omega limit sets, Omega the
# nonwandering set and R the chain recurrent set (robinson1999, sec. 2.3). The inclusion of the
# nonwandering set in the chain recurrent set is strict, and the map T(theta) = theta + alpha
# cos^2(theta) separates them by the whole circle: its nonwandering set is the two fixed points
# and its chain recurrent set is all of it.
#
# Method and references. Each quantity is computed and then compared with a closed form derived
# here from the displacement, and the derivations are stated beside the tests.
#
#   Nonwandering set at a grid of width h. Cell [a, a+h] meets its own forward image when
#   T(a) <= a + h and T(a+h) >= a, that is when disp(a) <= h and disp(a+h) >= -h, since T is an
#   increasing homeomorphism of the lift. A later image cannot meet the cell once the first does
#   not, because the motion of these maps is monotone and no orbit winds: the displacement of
#   cos(k theta) changes sign and blocks a turn, and the displacement of cos^2 vanishes at the
#   fixed points, which an orbit approaches and never passes. The computed set is therefore
#   {theta : |disp(theta)| <~ h} and its measure is 4 asin(c) with c = h/alpha for the maps with
#   displacement alpha cos(k theta) and c = sqrt(h/alpha) for alpha cos^2(theta).
#
#   Chain recurrent set at accuracy eps. The closed form is 4 asin(eps/alpha) where the
#   displacement changes sign and 2 pi where it does not; checks/chain-recurrence-of-circle-maps.R
#   verifies it against the strong components of the eps-transition graph.
#
#   Limit set L. Every orbit converges to a fixed point, so L is the fixed point set. The test
#   iterates a grid of starting points and measures how far the last iterate sits from it.
#
#   Alpha limit sets. The maps are homeomorphisms for k*alpha < 1, so the inverse exists and is
#   computed by bisection; the backward orbit converges to the fixed point behind the point.
#
# Degenerate cases: the identity map, every point of which is fixed and therefore nonwandering,
# and a grid so coarse that the whole circle is nonwandering.

set.seed(20260920L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

TWO_PI <- 2 * pi

systems <- list(
    a = list(tex = "$T(\\theta) = \\theta + \\alpha\\cos^2\\theta$", alpha = 0.4, k = NA_integer_,
             disp = function(th, alpha) alpha * cos(th)^2,
             fixed = c(pi / 2, 3 * pi / 2)),
    b = list(tex = "$T(\\theta) = \\theta + \\alpha\\cos\\theta$", alpha = 0.4, k = 1L,
             disp = function(th, alpha) alpha * cos(th),
             fixed = c(pi / 2, 3 * pi / 2)),
    c = list(tex = "$T(\\theta) = \\theta + \\alpha\\cos 3\\theta$", alpha = 0.25, k = 3L,
             disp = function(th, alpha) alpha * cos(3 * th),
             fixed = pi / 6 + (0:5) * pi / 3)
)

lift <- function(sys) function(th) th + sys$disp(th, sys$alpha)

# ============================ #
# The sets ####
# ============================ #

# The fixed points, found from the displacement itself and not taken from the list above. A sign
# change is not enough: the displacement of the squared cosine touches zero without crossing it,
# so the search is for the local minima of its modulus, each refined by a ternary search.
fixed_points_from_the_displacement <- function(sys, n_grid = 200000L, tol = 1e-10) {
    th <- seq(0, TWO_PI, length.out = n_grid + 1L)[-(n_grid + 1L)]
    d <- abs(sys$disp(th, sys$alpha))
    left <- c(d[n_grid], d[-n_grid])
    right <- c(d[-1], d[1])
    cand <- which(d <= left & d <= right & d < 0.25 * max(d))
    step <- TWO_PI / n_grid
    roots <- vapply(cand, function(i) {
        lo <- th[i] - step; hi <- th[i] + step
        for (s in seq_len(200L)) {
            m1 <- lo + (hi - lo) / 3; m2 <- hi - (hi - lo) / 3
            if (abs(sys$disp(m1, sys$alpha)) < abs(sys$disp(m2, sys$alpha))) hi <- m2 else lo <- m1
        }
        (lo + hi) / 2
    }, numeric(1))
    roots <- roots[abs(sys$disp(roots, sys$alpha)) < tol]
    roots <- sort(roots %% TWO_PI)
    keep <- c(TRUE, diff(roots) > 1e-6)
    roots[keep]
}

# The cells that meet one of their own forward images, over the first n_steps iterates.
nonwandering_cells <- function(sys, N, n_steps = 60L) {
    h <- TWO_PI / N
    a <- (0:(N - 1)) * h
    b <- a + h
    Tl <- lift(sys)
    A <- a; B <- b
    hit <- logical(N)
    for (n in seq_len(n_steps)) {
        A <- Tl(A); B <- Tl(B)
        # the image meets the cell modulo one turn when an integer k lies in the window below
        lo <- (A - b) / TWO_PI
        hi <- (B - a) / TWO_PI
        hit <- hit | (floor(hi) >= ceiling(lo))
    }
    list(hit = hit, h = h, measure = sum(hit) * h)
}

closed_form_nonwandering <- function(sys, h) {
    c_val <- if (is.na(sys$k)) sqrt(h / sys$alpha) else h / sys$alpha
    if (c_val >= 1) return(TWO_PI)
    4 * asin(c_val)
}

closed_form_chain_recurrent <- function(sys, eps) {
    if (is.na(sys$k)) return(TWO_PI)
    if (eps >= sys$alpha) return(TWO_PI)
    4 * asin(eps / sys$alpha)
}

# ============================ #
# Test 1: the nonwandering set against its closed form ####
# ============================ #
#
# The grid gives an outer approximation, since a cell is kept whole once any part of it meets an
# image. The bracket allows one cell at each end of each retained arc, which is 2 * (number of
# fixed points) * h.

test_nonwandering <- function(grids) {
    rows <- list()
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        for (N in grids) {
            r <- nonwandering_cells(sys, N)
            exact <- closed_form_nonwandering(sys, r$h)
            rows[[length(rows) + 1L]] <- data.frame(
                system = nm, cells = N, h = r$h, measured = r$measure, exact = exact,
                lower = exact - 2 * length(sys$fixed) * r$h,
                upper = exact + 2 * length(sys$fixed) * r$h,
                stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    out$inside <- out$measured >= out$lower - 1e-12 & out$measured <= out$upper + 1e-12
    out
}

# ============================ #
# Test 2: the rate at which the grid answer falls ####
# ============================ #
#
# The two closed forms give two rates: 4 asin(h/alpha) falls like h, and 4 asin(sqrt(h/alpha))
# like the square root of h, because the fixed points of the squared cosine are degenerate and
# the displacement beside them is quadratic and not linear. The ratio of the measure to the
# matching power of h therefore settles, and the ratio to the other power does not.

test_rate <- function(grids) {
    rows <- list()
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        m <- vapply(grids, function(N) nonwandering_cells(sys, N)$measure, numeric(1))
        h <- TWO_PI / grids
        n <- length(grids)
        rows[[length(rows) + 1L]] <- data.frame(
            system = nm, n_fixed = length(sys$fixed),
            ratio_linear_last = m[n] / h[n],
            ratio_linear_previous = m[n - 1L] / h[n - 1L],
            ratio_sqrt_last = m[n] / sqrt(h[n]),
            predicted_linear = if (is.na(sys$k)) NA_real_ else 4 / sys$alpha,
            predicted_sqrt = if (is.na(sys$k)) 4 / sqrt(sys$alpha) else NA_real_,
            stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 3: the inclusions ####
# ============================ #
#
# Every fixed point sits in a nonwandering cell, every nonwandering cell sits inside the chain
# recurrent set at accuracy h, and the last inclusion is strict for the first map by the whole
# circle. The chain recurrent set is taken from its closed form, which the companion check
# verifies against a graph computation.

test_inclusions <- function(N) {
    rows <- list()
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        r <- nonwandering_cells(sys, N)
        centres <- ((which(r$hit) - 1) + 0.5) * r$h
        fixed <- fixed_points_from_the_displacement(sys)
        in_omega <- vapply(fixed, function(p) {
            i <- floor((p %% TWO_PI) / r$h) + 1L
            r$hit[min(i, N)]
        }, logical(1))
        # A nonwandering cell satisfies disp(a) <= h at its left end and disp(a + h) >= -h at
        # its right end, so by the intermediate value theorem it holds a point with |disp| <= h,
        # which is a point of the chain recurrent set at accuracy h. The centre of such a cell is
        # within h/2 of that point, so its displacement exceeds h by at most L h / 2, where L is
        # the largest slope of the displacement.
        L <- if (is.na(sys$k)) sys$alpha else sys$k * sys$alpha
        in_r <- if (is.na(sys$k)) rep(TRUE, length(centres)) else
            abs(sys$disp(centres, sys$alpha)) <= r$h * (1 + L / 2) + 1e-12
        rows[[length(rows) + 1L]] <- data.frame(
            system = nm, n_fixed = length(fixed),
            fixed_in_omega = all(in_omega),
            omega_in_chain = all(in_r),
            measure_omega = r$measure,
            measure_chain = closed_form_chain_recurrent(sys, r$h),
            stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 4: the limit set and the alpha limit sets ####
# ============================ #
#
# Forward, every orbit reaches a fixed point, so the limit set L is the fixed point set.
# Backward, the inverse is computed by bisection on the increasing lift and the orbit reaches the
# fixed point behind the starting point, which for the second map is the single repelling one.

# The inverse of the lift, by bisection, vectorised over the points.
invert <- function(sys, y, halfwidth = 4) {
    Tl <- lift(sys)
    lo <- y - halfwidth; hi <- y + halfwidth
    for (i in seq_len(90L)) {
        mid <- (lo + hi) / 2
        left <- Tl(mid) < y
        lo <- ifelse(left, mid, lo)
        hi <- ifelse(left, hi, mid)
    }
    (lo + hi) / 2
}

dist_to <- function(v, fixed) {
    vapply(v, function(u) min(abs(((u - fixed + pi) %% TWO_PI) - pi)), numeric(1))
}

test_limit_sets <- function(n_points = 401L, n_iter = 40000L, n_back = 800L) {
    rows <- list()
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        Tl <- lift(sys)
        fixed <- fixed_points_from_the_displacement(sys)
        x <- seq(0, TWO_PI, length.out = n_points + 1L)[-(n_points + 1L)]
        # the starting points are kept away from every fixed point, since a point beside the one
        # ahead of it must first escape that one backwards, and beside a degenerate fixed point
        # the escape is as slow as the approach and would dominate the horizon
        x <- x[dist_to(x, fixed) > 0.2]
        fwd <- x
        for (i in seq_len(n_iter)) fwd <- Tl(fwd)
        bwd <- x; d_half <- NA_real_
        for (i in seq_len(n_back)) {
            bwd <- invert(sys, bwd)
            if (i == n_back %/% 2L) d_half <- max(dist_to(bwd, fixed))
        }
        d_back <- max(dist_to(bwd, fixed))
        # for the second map every backward orbit reaches the one repelling fixed point
        to_repeller <- if (nm == "b") max(abs(((bwd + pi / 2 + pi) %% TWO_PI) - pi)) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
            system = nm, n_tested = length(x), n_fixed = length(fixed),
            max_distance_forward = max(dist_to(fwd, fixed)),
            max_distance_backward = d_back,
            halving_ratio = d_back / d_half,
            n_times_distance_backward = n_back * d_back,
            algebraic_target = 1 / sys$alpha,
            backward_to_repeller = to_repeller, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 5: the degenerate cases ####
# ============================ #

test_degenerate <- function() {
    idm <- systems$b; idm$alpha <- 0; idm$fixed <- numeric(0)
    r_id <- nonwandering_cells(idm, 720L)
    coarse <- nonwandering_cells(systems$b, 8L)     # a cell wider than the largest displacement
    data.frame(identity_measure = r_id$measure, coarse_measure = coarse$measure, target = TWO_PI)
}

# ============================ #
# Run ####
# ============================ #

grids <- c(500L, 1000L, 2000L, 4000L, 8000L)
t_nw <- test_nonwandering(grids)
t_rate <- test_rate(grids)
t_inc <- test_inclusions(4000L)
t_lim <- test_limit_sets()
t_deg <- test_degenerate()

ok_nw <- all(t_nw$inside)
rel_err <- max(abs(t_nw$measured - t_nw$exact) / t_nw$exact)
# The two laws. Where the displacement changes sign the arcs of the nonwandering set are a fixed
# number of cells wide at every grid, 2/(k alpha) of them, so the measure falls like h and the
# ratio to h settles, at a value between the closed form 4/alpha and that plus two cells for each
# arc, which is what rounding a partial cell up at each end costs. The ratio is therefore not
# expected to reach the closed form, and what the test asks is that it sit in that window and not
# move between the two finest grids. At a degenerate fixed point the arc is of order sqrt(h) and
# holds many cells, so the rounding vanishes and the ratio to sqrt(h) does reach its closed form.
rate_window <- function(nm) {
    r <- t_rate[t_rate$system == nm, ]
    r$ratio_linear_last >= r$predicted_linear - 1e-9 &&
        r$ratio_linear_last <= r$predicted_linear + 2 * r$n_fixed &&
        abs(r$ratio_linear_last - r$ratio_linear_previous) <= 2 * r$n_fixed
}
ok_rate <- rate_window("b") && rate_window("c") &&
           abs(t_rate$ratio_sqrt_last[t_rate$system == "a"] - t_rate$predicted_sqrt[t_rate$system == "a"]) <
               0.01 * t_rate$predicted_sqrt[t_rate$system == "a"]
# The inclusion is strict for the first map, and by a factor the grid fixes: the chain recurrent
# set is the whole circle at every accuracy and the nonwandering set at cell width h measures
# 4 asin(sqrt(h/alpha)), so their ratio is 2 pi over that, which is above twenty at this grid.
ok_inc <- all(t_inc$fixed_in_omega) && all(t_inc$omega_in_chain) &&
          t_inc$measure_chain[t_inc$system == "a"] == TWO_PI &&
          t_inc$measure_chain[t_inc$system == "a"] / t_inc$measure_omega[t_inc$system == "a"] > 20
# Forward every orbit reaches a fixed point. Backward the two rates part again: the hyperbolic
# fixed points of the second and third maps collapse the distance to rounding level, while the
# degenerate ones of the first give the algebraic law, the distance falling like 1/(alpha n), so
# the product with the step count approaches 1/alpha from above and the distance halves when the
# step count doubles.
ok_lim <- all(t_lim$max_distance_forward < 2e-3) &&
          all(t_lim$max_distance_backward[t_lim$system %in% c("b", "c")] < 1e-6) &&
          t_lim$backward_to_repeller[t_lim$system == "b"] < 1e-6 &&
          with(t_lim[t_lim$system == "a", ],
               n_times_distance_backward > algebraic_target &&
               n_times_distance_backward < 2 * algebraic_target &&
               halving_ratio > 0.42 && halving_ratio < 0.58)
ok_deg <- abs(t_deg$identity_measure - TWO_PI) < 1e-12 && abs(t_deg$coarse_measure - TWO_PI) < 1e-12

status <- if (ok_nw && ok_rate && ok_inc && ok_lim && ok_deg) "pass" else "fail"

# ============================ #
# Figure ####
# ============================ #
#
# The nonwandering set and the chain recurrent set of the same map at the same scale, against
# the resolution at which each is asked for. The two curves separate for the first map, where the
# nonwandering set falls to two points and the chain recurrent set stays the whole circle.

if (requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages(library(ggplot2))
    labs_tex <- vapply(systems, `[[`, character(1), "tex")
    fine <- exp(seq(log(1e-4), log(0.35), length.out = 240))
    curves <- do.call(rbind, lapply(names(systems), function(nm) {
        sys <- systems[[nm]]
        rbind(
            data.frame(label = labs_tex[[nm]], scale = fine, set = "Nonwandering set",
                       measure = vapply(fine, function(u) closed_form_nonwandering(sys, u), numeric(1)),
                       stringsAsFactors = FALSE),
            data.frame(label = labs_tex[[nm]], scale = fine, set = "Chain recurrent set",
                       measure = vapply(fine, function(u) closed_form_chain_recurrent(sys, u), numeric(1)),
                       stringsAsFactors = FALSE))
    }))
    pts <- data.frame(label = labs_tex[t_nw$system], scale = t_nw$h, measure = t_nw$measured,
                      set = "Nonwandering set", stringsAsFactors = FALSE)
    curves$label <- factor(curves$label, levels = labs_tex)
    pts$label <- factor(pts$label, levels = labs_tex)
    p <- ggplot(curves, aes(scale, measure, colour = set, linewidth = set, alpha = set)) +
        geom_line() +
        geom_point(data = pts, aes(scale, measure), inherit.aes = FALSE,
                   size = 1.9, shape = 21, fill = "white", stroke = 0.7, colour = "#056796") +
        facet_wrap(~ label, nrow = 1, labeller = kb_labeller()) +
        scale_x_log10() +
        scale_y_continuous(limits = c(0, TWO_PI + 0.2),
                           breaks = c(0, pi, TWO_PI),
                           labels = kb_ticks(c("0", "$\\pi$", "$2\\pi$"))) +
        scale_colour_manual(values = c("Nonwandering set" = "#056796",
                                       "Chain recurrent set" = "#be1117"), name = NULL) +
        scale_linewidth_manual(values = c("Nonwandering set" = 0.55,
                                          "Chain recurrent set" = 2.2), name = NULL) +
        scale_alpha_manual(values = c("Nonwandering set" = 1, "Chain recurrent set" = 0.4),
                           name = NULL) +
        labs(title = kb_unicode("Two recurrent sets of the same map, against the scale at which they are asked for"),
             subtitle = kb_unicode("Lebesgue measure of the nonwandering set at cell width $h$ and of the chain recurrent set at jump size $\\varepsilon$"),
             x = kb_tex("Cell width $h$, or jump size $\\varepsilon$"),
             y = kb_tex("Measure"),
             caption = kb_caption(paste(
                 "Thin line: the measure of the set of cells that meet one of their own forward images, which",
                 "is the nonwandering set seen at cell width $h$, with the circles giving the computed value at",
                 "five grids. Wide band: the measure of the chain recurrent set at jump size $\\varepsilon$,",
                 "which the companion check verifies against a graph computation; where the thin line runs",
                 "inside the band the two sets have the same measure. They agree where the",
                 "displacement changes sign and part company on the first map, whose nonwandering set falls to",
                 "its two fixed points like $\\sqrt{h}$ while its chain recurrent set is the whole circle at",
                 "every $\\varepsilon$: a chain may jump past a fixed point and a true orbit may not. Drawn by",
                 "checks/recurrence-hierarchy-of-circle-maps.R with $\\alpha = 0.4$ for the first two maps and",
                 "$\\alpha = 0.25$ for the third."))) +
        kb_theme()
    kb_save(p, "recurrence-hierarchy-of-circle-maps", width = 8.6, height = 4.4)
}

emit("recurrence-hierarchy-of-circle-maps", status,
     "The nonwandering set of a circle map at cell width h has measure 4 asin(h/alpha), or 4 asin(sqrt(h/alpha)) at a degenerate fixed point, and it lies strictly inside the chain recurrent set",
     list(max_relative_error_nonwandering = rel_err,
          n_grid_cases = nrow(t_nw),
          all_within_one_cell_bracket = as.numeric(ok_nw),
          linear_rate_cos_map = t_rate$ratio_linear_last[t_rate$system == "b"],
          linear_rate_cos_map_closed_form = t_rate$predicted_linear[t_rate$system == "b"],
          linear_rate_cos3_map = t_rate$ratio_linear_last[t_rate$system == "c"],
          linear_rate_cos3_map_closed_form = t_rate$predicted_linear[t_rate$system == "c"],
          sqrt_rate_cos_squared_map = t_rate$ratio_sqrt_last[t_rate$system == "a"],
          sqrt_rate_cos_squared_map_closed_form = t_rate$predicted_sqrt[t_rate$system == "a"],
          fixed_points_all_nonwandering = as.numeric(all(t_inc$fixed_in_omega)),
          nonwandering_inside_chain_recurrent = as.numeric(all(t_inc$omega_in_chain)),
          measure_nonwandering_cos_squared_map = t_inc$measure_omega[t_inc$system == "a"],
          measure_chain_recurrent_cos_squared_map = t_inc$measure_chain[t_inc$system == "a"],
          max_distance_to_a_fixed_point_forward = max(t_lim$max_distance_forward),
          max_distance_to_a_fixed_point_backward_hyperbolic = max(t_lim$max_distance_backward[t_lim$system %in% c("b", "c")]),
          backward_product_cos_squared_map = t_lim$n_times_distance_backward[t_lim$system == "a"],
          backward_product_closed_form = t_lim$algebraic_target[t_lim$system == "a"],
          backward_halving_ratio_cos_squared_map = t_lim$halving_ratio[t_lim$system == "a"],
          backward_orbits_to_the_repeller_cos_map = t_lim$backward_to_repeller[t_lim$system == "b"],
          identity_map_measure = t_deg$identity_measure,
          coarse_grid_measure = t_deg$coarse_measure))
