# ============================ #
# Return map and time fractions of the circle flow with a degenerate equilibrium ####
# ============================ #
#
# Claims, with the locators of oljaca2024 (section 3.3, eq. (13) to eq. (20)):
#
# The flow on S^1 x [0, 1], with S^1 identified with [0, 1], is given by
# theta' = y^3 where theta <= y, theta' = 1 where theta > y, and y' = -y^2.
#
# 1. Eq. (16) and eq. (17): the orbits satisfy theta = (x^2 - y^2) / 2 in the region
#    theta <= y, where x is the y-coordinate of the last crossing of the section
#    {theta = 0}, and theta = 1/y + u - 1/u in the region theta > y, where u is the
#    y-coordinate at the crossing of the diagonal theta = y.
# 2. The diagonal is met at u = r(x) = sqrt(1 + x^2) - 1, the section is met again at
#    y = s(u) = 1 / (1 + 1/u - u), and the return map is F = s o r. Reference: a
#    fourth-order Runge-Kutta solution of the piecewise system, which shares no code with
#    the closed forms.
# 3. Eq. (15): F(x) = x^2/2 - 3 x^4/8 + O(x^6) as x tends to 0.
# 4. Consistency of the two descriptions of the fast leg: the time 1 - u taken by theta to
#    run from u to 1 equals the time 1/s(u) - 1/u taken by y to fall from u to s(u).
# 5. The orbit through (0, x) leaves every neighbourhood of the origin, since theta passes
#    through 1/2 on every loop, so the omega-limit set is the whole circle y = 0 and the
#    origin is not an attractor in the topological sense. The fraction of time spent within
#    distance eps of the origin nevertheless converges to 1, which is statistical
#    attraction to the single point (0, 0).
# 6. Degenerate cases: x = 1 on the boundary of the section, and the fixed point at the
#    origin.
#
# Times are handled leg by leg, because the slow leg lasts about 2 / x^2 and adding a
# duration of order one to it would be lost to rounding after a few loops.

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

# Numerically stable forms of the two legs of the return map.
r_map <- function(x) x^2 / (sqrt(1 + x^2) + 1)          # sqrt(1 + x^2) - 1
s_map <- function(u) u / (u + 1 - u^2)                  # 1 / (1 + 1/u - u)
F_map <- function(x) s_map(r_map(x))

t_slow <- function(x) 1 / r_map(x) - 1 / x              # y falls from x to r(x)
t_fast <- function(u) 1 - u                             # theta runs from u to 1

slow_position <- function(tau, x) {
    y <- 1 / (1 / x + tau)
    c((x^2 - y^2) / 2, y)
}
fast_position <- function(tau, u) {
    c(u + tau, 1 / (1 / u + tau))
}
dist_origin <- function(p) sqrt(min(p[1], 1 - p[1])^2 + p[2]^2)

# ============================ #
# Test 1: closed forms against a Runge-Kutta solution ####
# ============================ #
#
# Tolerance 1e-9: the field is discontinuous across the diagonal theta = y, so each leg is
# integrated on its own with the field that governs it, and each leg is smooth. With step
# 1e-5 the observed Runge-Kutta error over a leg is of order 1e-12.

rk4_field <- function(v0, tend, f, h = 1e-5) {
    n <- as.integer(ceiling(tend / h))
    h <- tend / n
    v <- v0
    for (i in seq_len(n)) {
        k1 <- f(v); k2 <- f(v + h * k1 / 2); k3 <- f(v + h * k2 / 2); k4 <- f(v + h * k3)
        v <- v + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    }
    v
}
f_slow <- function(v) c(v[2]^3, -v[2]^2)        # region theta <= y
f_fast <- function(v) c(1, -v[2]^2)             # region theta > y

err_slow_leg <- 0
err_fast_leg <- 0
err_diagonal <- 0
for (x in c(0.4, 0.6, 0.9)) {
    u <- r_map(x)
    end_slow <- rk4_field(c(0, x), t_slow(x), f_slow)
    err_slow_leg <- max(err_slow_leg, abs(end_slow[1] - u), abs(end_slow[2] - u))
    err_diagonal <- max(err_diagonal, abs(end_slow[1] - end_slow[2]))
    end_fast <- rk4_field(c(u, u), t_fast(u), f_fast)
    err_fast_leg <- max(err_fast_leg, abs(end_fast[1] - 1), abs(end_fast[2] - F_map(x)))
}
err_return <- err_fast_leg

# Positions inside each leg, compared with the closed forms.
err_position <- 0
x_ref <- 0.6
u_ref <- r_map(x_ref)
for (frac in c(0.25, 0.75)) {
    tend <- frac * t_slow(x_ref)
    err_position <- max(err_position,
                        max(abs(rk4_field(c(0, x_ref), tend, f_slow) - slow_position(tend, x_ref))))
    tend <- frac * t_fast(u_ref)
    err_position <- max(err_position,
                        max(abs(rk4_field(c(u_ref, u_ref), tend, f_fast) - fast_position(tend, u_ref))))
}

# ============================ #
# Test 2: the series expansion of the return map ####
# ============================ #
#
# The ratio |F(x) - (x^2/2 - 3 x^4/8)| / x^6 must stay bounded and bounded away from zero
# as x tends to 0, which identifies the order of the remainder as exactly six.

xs <- 2^-(2:9)
ratio <- abs(F_map(xs) - (xs^2 / 2 - 3 * xs^4 / 8)) / xs^6
ratio_max <- max(ratio)
ratio_min <- min(ratio)

# ============================ #
# Test 3: consistency of the fast leg ####
# ============================ #

us <- r_map(c(0.05, 0.2, 0.5, 0.9, 1))
err_fast <- max(abs((1 / s_map(us) - 1 / us) - (1 - us)))

# ============================ #
# Test 4: time fractions and the omega-limit set ####
# ============================ #
#
# In-ball time within one loop that starts at (0, x), by bisection on the closed-form
# distance in each leg. Along the slow leg the distance falls monotonically from x to
# u sqrt(2), so the whole leg is inside as soon as x < eps. Along the fast leg the distance
# rises to about 1/2 at the far side of the circle and falls again, so the leg contributes
# an exit window and an entry window.

in_ball_time <- function(x, eps) {
    u <- r_map(x)
    ts <- t_slow(x)
    tf <- t_fast(u)
    d_slow <- function(tau) dist_origin(slow_position(tau, x)) - eps
    d_fast <- function(tau) dist_origin(fast_position(tau, u)) - eps
    slow_part <- if (d_slow(0) < 0) ts else if (d_slow(ts) >= 0) 0 else
        ts - uniroot(d_slow, c(0, ts), tol = 1e-12)$root
    half <- tf / 2
    exit <- if (d_fast(0) >= 0) 0 else if (d_fast(half) < 0) half else
        uniroot(d_fast, c(0, half), tol = 1e-14)$root
    entry <- if (d_fast(tf) >= 0) 0 else if (d_fast(half) < 0) half else
        tf - uniroot(d_fast, c(half, tf), tol = 1e-14)$root
    slow_part + exit + entry
}

eps_ball <- 0.1
x0 <- 0.5
xs_seq <- numeric(0)
t_in <- 0
t_tot <- 0
x <- x0
fractions <- numeric(0)
for (k in 1:7) {
    xs_seq <- c(xs_seq, x)
    t_in <- t_in + in_ball_time(x, eps_ball)
    t_tot <- t_tot + t_slow(x) + t_fast(r_map(x))
    fractions <- c(fractions, t_in / t_tot)
    x <- F_map(x)
}
fraction_final <- fractions[length(fractions)]
fraction_first <- fractions[1]

# Cross-check of one loop by summing an indicator on a fine grid; the grid step is 1e-4 and
# there are three boundary crossings, so the tolerance is 1e-3.
x_chk <- 0.3
u_chk <- r_map(x_chk)
ts_chk <- t_slow(x_chk)
tf_chk <- t_fast(u_chk)
n_slow <- as.integer(ts_chk / 1e-4)
n_fast <- as.integer(tf_chk / 1e-4)
grid_slow <- mean(vapply(seq(0, ts_chk, length.out = n_slow),
                         function(tau) dist_origin(slow_position(tau, x_chk)) < eps_ball,
                         logical(1))) * ts_chk
grid_fast <- mean(vapply(seq(0, tf_chk, length.out = n_fast),
                         function(tau) dist_origin(fast_position(tau, u_chk)) < eps_ball,
                         logical(1))) * tf_chk
err_in_ball <- abs((grid_slow + grid_fast) - in_ball_time(x_chk, eps_ball))

# The orbit passes through theta = 1/2 on every loop, so it leaves the ball of radius 1/4.
max_dist_per_loop <- vapply(xs_seq, function(x) {
    u <- r_map(x)
    dist_origin(fast_position(0.5 - u, u))
}, numeric(1))
min_max_dist <- min(max_dist_per_loop)
sup_y_on_last_loop <- xs_seq[length(xs_seq)]

# ============================ #
# Test 5: degenerate cases ####
# ============================ #

err_boundary <- abs(r_map(1) - (sqrt(2) - 1))
err_fixed <- max(abs(c(r_map(0), s_map(0))))    # the origin is a fixed point of both legs

ok <- err_slow_leg < 1e-9 && err_fast_leg < 1e-9 && err_diagonal < 1e-9 &&
    err_position < 1e-9 &&
    ratio_max < 1 && ratio_min > 0.01 &&
    err_fast < 1e-12 && err_in_ball < 1e-3 &&
    fraction_final > 1 - 1e-12 && fraction_first < 0.9 &&
    min_max_dist > 0.49 && sup_y_on_last_loop < 1e-30 &&
    err_boundary < 1e-15 && err_fixed == 0

emit("circle-flow-degenerate-equilibrium", if (ok) "pass" else "fail",
     "The return map of the circle flow is s(r(x)) with the expansion x^2/2 - 3x^4/8, and the orbit spends an asymptotically full fraction of its time within any ball about the degenerate equilibrium while still crossing the far side of the circle on every loop",
     list(rk4_slow_leg_max_abs_error = err_slow_leg,
          rk4_fast_leg_max_abs_error = err_fast_leg,
          rk4_diagonal_crossing_abs_error = err_diagonal,
          rk4_position_max_abs_error = err_position,
          series_ratio_max = ratio_max,
          series_ratio_min = ratio_min,
          fast_leg_consistency_error = err_fast,
          in_ball_time_grid_abs_error = err_in_ball,
          time_fraction_loop1 = fraction_first,
          time_fraction_loop7 = fraction_final,
          section_value_x2 = xs_seq[2],
          section_value_x4 = xs_seq[4],
          section_value_x7 = xs_seq[7],
          min_far_side_distance = min_max_dist,
          boundary_r_abs_error = err_boundary,
          fixed_point_error = err_fixed))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the running share of time spent within distance 0.1 of the degenerate equilibrium,
# against the greatest distance the orbit reaches on the same loop. The share rises to one while
# the orbit still travels to the far side of the circle on every loop, which is the separation
# between statistical attraction and convergence. Right: the section values, which fall to zero
# like the iterates of the return map, on a logarithmic scale.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    loops <- seq_along(fractions)
    left <- rbind(
        data.frame(loop = loops, value = fractions,
                   what = "Share of time within 0.1 of the equilibrium"),
        data.frame(loop = loops, value = max_dist_per_loop,
                   what = "Greatest distance reached on the loop"))
    p1 <- ggplot(left, aes(loop, value, colour = what, shape = what)) +
        geom_line(linewidth = 0.5) + geom_point(size = 1.7) +
        scale_colour_manual(values = c("Share of time within 0.1 of the equilibrium" = "#21918c",
                                       "Greatest distance reached on the loop" = "#B8390E"),
                            name = NULL) +
        scale_shape_manual(values = c(16, 17), name = NULL) +
        scale_y_continuous(limits = c(0, 1.02)) +
        scale_x_continuous(breaks = loops) +
        labs(title = "Almost all of the time near a point the orbit keeps leaving",
             subtitle = kb_unicode(sprintf("Flow on the cylinder, from the section at $y = %.1f$, over %d loops",
                                           x0, length(loops))),
             x = "Loop", y = kb_tex("Share of time, or distance to $(0,0)$")) +
        kb_theme()
    p2 <- ggplot(data.frame(loop = loops, y = xs_seq), aes(loop, y)) +
        geom_line(colour = "#440154", linewidth = 0.5) +
        geom_point(colour = "#440154", size = 1.7) +
        scale_y_log10() + scale_x_continuous(breaks = loops) +
        labs(title = "The section values under the return map",
             subtitle = "Each loop starts closer to the equilibrium and lasts longer",
             x = "Loop", y = kb_tex("Value $x_k$ on the section")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.35, 1)) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "In the run recorded by checks/circle-flow-degenerate-equilibrium.R the share of time was %.3f on the first loop and %.12f on the seventh, while the orbit reached distance %.2f from the equilibrium on every one of them. The point is therefore a statistical attractor and attracts no orbit in the ordinary sense.",
                fraction_first, fraction_final, min_max_dist)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "circle-flow-degenerate-equilibrium", width = 9.4, height = 4.2))
}
