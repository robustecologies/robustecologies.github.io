# ============================ #
# Return times and pullback overtaking near a hyperbolic homoclinic loop ####
# ============================ #
#
# Claims, with the locators of newman2025 (section 3.5 and section 5.6):
#
# A planar vector field is linear, with f(x, y) = (lam2 x, lam1 y), on the square
# V1 = [0, 1]^2, where lam1 > 0 and lam2 < -lam1; a trajectory that leaves V1 at (z, 1)
# re-enters it at (1, c z) after an excursion of duration tau(z), bounded by M. The
# equilibrium p = (0, 0) then has a homoclinic loop. Write r = -lam2 / lam1 > 1,
# sigma(y) = -log(y) / lam1 for the time that the trajectory entering at (1, y) stays in
# V1, and G(y) = c y^r for the y-coordinate of the next entry. This script uses
# lam1 = 1, lam2 = -1.5, c = 0.5 and tau(z) = 0.3 (1 + z), so M = 0.6.
#
# 1. The linear flow in V1 carries (1, y) to (y^r, 1) in time sigma(y). Reference: a
#    fourth-order Runge-Kutta solution of the linear system, which shares no code with the
#    closed form.
# 2. Eq. (12): sigma_n(y) = sigma(G^n(y)) = -[ (r^n - 1) / (r - 1) log c + r^n log y ] / lam1,
#    where G^n is the n-th iterate of G. Reference: direct iteration of G in logarithms.
# 3. Lemma 5.4: sigma_n(y') - sigma_n(y) = r^n (sigma(y') - sigma(y)).
# 4. Corollary 5.3 and BRT attraction: for 0 < eps <= 1 the trajectory entering at (1, y)
#    lies in V_eps = [0, eps]^2 exactly on the sub-interval [u_n + m2, v_n - m1] of its
#    n-th visit to V1, with m1 = -log(eps)/lam1 and m2 = log(eps)/lam2; the times between
#    consecutive visits to V_eps are therefore bounded by M + m1 + m2, and the durations
#    inside V_eps grow without bound. The point p is BRT attracted from a set of positive
#    measure, hence statistically attracted, and the fraction of time spent in V_eps
#    converges to 1.
# 5. Statement 1 of proposition 3.5: no trajectory that enters V1 at (1, y) converges to p,
#    since it returns to x = 1 at every loop; {p} is therefore not a measure attractor.
# 6. Statement 3 of proposition 3.5: with y1 chosen so that sigma(G(y1)) - sigma(y1) > 2M,
#    with sigma(y2) the midpoint of sigma(G(y1)) and sigma(y1), and with alpha > 0 small,
#    the two sets S1 = {phi_s(1, y) : s in [0, beta], y in [y1, y1 + alpha]} and
#    S2 = {phi_s(1, y) : s in [0, beta], y in [y2 - alpha, y2]} take turns: for every
#    sufficiently large t at least one of phi_t(S1) and phi_t(S2) is contained in V_eps.
#    The script verifies this alternation on a grid of times and initial conditions.
# 7. Boundary case r = 1 (lam2 = -lam1): the differences sigma_n(y') - sigma_n(y) are then
#    constant in n, so the exponential growth of the excursions is lost, while the return
#    times stay bounded.

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

lam1 <- 1
lam2 <- -1.5
r <- -lam2 / lam1
cc <- 0.5
tau <- function(z) 0.3 * (1 + z)
M <- tau(1)

sigma <- function(y) -log(y) / lam1
Gmap <- function(y) cc * y^r
log_Gn <- function(y, n) ((r^n - 1) / (r - 1)) * log(cc) + r^n * log(y)
sigma_n <- function(y, n) -log_Gn(y, n) / lam1

# ============================ #
# Test 1: the linear flow inside the square ####
# ============================ #
#
# Tolerance 1e-9: with step 1e-5 the Runge-Kutta error over an interval of length at most
# 4 is of order 1e-12.

rk4_box <- function(y, h = 1e-5) {
    tend <- sigma(y)
    n <- as.integer(ceiling(tend / h))
    h <- tend / n
    z <- c(1, y)
    f <- function(v) c(lam2 * v[1], lam1 * v[2])
    for (i in seq_len(n)) {
        k1 <- f(z); k2 <- f(z + h * k1 / 2); k3 <- f(z + h * k2 / 2); k4 <- f(z + h * k3)
        z <- z + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    }
    z
}
err_exit <- 0
for (y in c(0.02, 0.1, 0.3, 0.8)) {
    z <- rk4_box(y)
    err_exit <- max(err_exit, abs(z[1] - y^r), abs(z[2] - 1))
}

# ============================ #
# Test 2: the closed form of the iterated return map ####
# ============================ #
#
# Iteration is done in logarithms, so the comparison is exact to rounding; tolerance 1e-10.

err_Gn <- 0
err_sigma_n <- 0
for (y in c(0.05, 0.2, 0.5, 0.9)) {
    lg <- log(y)
    for (n in 0:12) {
        err_Gn <- max(err_Gn, abs(lg - log_Gn(y, n)))
        err_sigma_n <- max(err_sigma_n, abs(-lg / lam1 - sigma_n(y, n)))
        lg <- log(cc) + r * lg
    }
}

# ============================ #
# Test 3: lemma 5.4 ####
# ============================ #

err_lemma54 <- 0
for (n in 0:10) {
    for (pair in list(c(0.1, 0.4), c(0.05, 0.9), c(0.3, 0.31))) {
        lhs <- sigma_n(pair[2], n) - sigma_n(pair[1], n)
        rhs <- r^n * (sigma(pair[2]) - sigma(pair[1]))
        err_lemma54 <- max(err_lemma54, abs(lhs - rhs) / abs(rhs))
    }
}

# ============================ #
# Test 4: itineraries, bounded return times and time fraction ####
# ============================ #
#
# The n-th visit to V1 starts at u_n and ends at v_n = u_n + sigma_n(y); the excursion that
# follows lasts tau(G^n(y)^r / cc ... ) evaluated at the exit abscissa z_n = G^n(y)^r.

itinerary <- function(y, nmax) {
    u <- numeric(nmax + 1L); v <- numeric(nmax + 1L)
    t <- 0
    for (n in 0:nmax) {
        u[n + 1L] <- t
        s <- sigma_n(y, n)
        v[n + 1L] <- t + s
        z <- exp(r * log_Gn(y, n))          # exit abscissa (G^n(y))^r
        t <- v[n + 1L] + tau(z)
    }
    list(u = u, v = v)
}

eps_list <- c(0.5, 0.1, 0.01)
y_start <- 0.3
nmax <- 14L
it <- itinerary(y_start, nmax)
gap_max <- numeric(length(eps_list))
fraction <- numeric(length(eps_list))
window_growth <- numeric(length(eps_list))
for (i in seq_along(eps_list)) {
    eps <- eps_list[i]
    m1 <- -log(eps) / lam1
    m2 <- log(eps) / lam2
    starts <- it$u + m2
    ends <- it$v - m1
    keep <- ends > starts
    starts <- starts[keep]; ends <- ends[keep]
    gap_max[i] <- max(starts[-1] - ends[-length(ends)])
    horizon <- ends[length(ends)]
    inside <- sum(pmax(0, pmin(ends, horizon) - starts))
    fraction[i] <- inside / horizon
    window_growth[i] <- (ends - starts)[length(ends)] / (ends - starts)[1]
}
bound <- M + (-log(eps_list) / lam1) + (log(eps_list) / lam2)
gap_within_bound <- all(gap_max <= bound + 1e-9)

# ============================ #
# Test 5: no convergence to p ####
# ============================ #
#
# At every re-entry the trajectory is at (1, G^n(y)), so its distance to the origin is at
# least 1 for arbitrarily large times.

dist_at_reentry <- min(sqrt(1 + exp(2 * log_Gn(y_start, 0:nmax))))

# ============================ #
# Test 6: the alternation of two sets under pullback ####
# ============================ #

y1 <- 0.2
stopifnot(sigma(Gmap(y1)) - sigma(y1) > 2 * M)
margin <- sigma(Gmap(y1)) - sigma(y1)
y2 <- exp(-(sigma(Gmap(y1)) + sigma(y1)) / 2 * lam1)     # sigma(y2) is the midpoint
stopifnot(y2 > Gmap(y1), y2 < y1)
# Largest alpha for which sigma(G(y1 + alpha)) - sigma(y2 - alpha) > M, by bisection.
f_alpha <- function(a) sigma(Gmap(y1 + a)) - sigma(y2 - a) - M
alpha_max <- uniroot(f_alpha, lower = 1e-9, upper = min(y2 - Gmap(y1), 1 - y1) * 0.999,
                     tol = 1e-12)$root
alpha <- alpha_max / 2
beta <- 0.05

eps_alt <- 0.1
m1 <- -log(eps_alt) / lam1
m2 <- log(eps_alt) / lam2

in_V_eps_interval <- function(y, t, beta, nmax) {
    # TRUE when the whole segment [t, t + beta] of the trajectory entering at (1, y) lies
    # in V_eps, that is when it sits inside one window [u_n + m2, v_n - m1].
    itn <- itinerary(y, nmax)
    s <- itn$u + m2
    e <- itn$v - m1
    any(s <= t & t + beta <= e)
}

ys1 <- seq(y1, y1 + alpha, length.out = 25L)
ys2 <- seq(y2 - alpha, y2, length.out = 25L)
set_inside <- function(ys, t) all(vapply(ys, function(y) in_V_eps_interval(y, t, beta, 20L), logical(1)))

t_grid <- seq(0, 600, by = 0.05)
ok_alt <- vapply(t_grid, function(t) set_inside(ys1, t) || set_inside(ys2, t), logical(1))
last_failure <- if (any(!ok_alt)) max(t_grid[!ok_alt]) else -Inf
frac_covered_after <- mean(ok_alt[t_grid > last_failure])
n_times_after <- sum(t_grid > last_failure)

# ============================ #
# Test 7: the boundary case r = 1 ####
# ============================ #

r_b <- 1
sigma_n_b <- function(y, n) (-log(y) - n * log(cc)) / lam1      # G(y) = c y when r = 1
err_r1 <- 0
for (n in 0:10) {
    lhs <- sigma_n_b(0.4, n) - sigma_n_b(0.1, n)
    rhs <- sigma(0.4) - sigma(0.1)
    err_r1 <- max(err_r1, abs(lhs - rhs))
}

ok <- err_exit < 1e-9 && err_Gn < 1e-10 && err_sigma_n < 1e-10 && err_lemma54 < 1e-12 &&
    gap_within_bound && all(diff(fraction) < 0) && fraction[1] > 0.99 &&
    all(window_growth > 100) && dist_at_reentry >= 1 &&
    is.finite(last_failure) && last_failure < 400 && frac_covered_after == 1 &&
    n_times_after > 3000 && err_r1 < 1e-12

emit("homoclinic-loop-return-times", if (ok) "pass" else "fail",
     "Excursions from a hyperbolic homoclinic loop return in bounded time while the visits grow like r^n, and two sets of initial conditions alternate in the neighbourhood of the equilibrium under backward time",
     list(r = r,
          M = M,
          rk4_exit_max_abs_error = err_exit,
          iterated_map_max_abs_error = err_Gn,
          sigma_n_max_abs_error = err_sigma_n,
          lemma54_max_rel_error = err_lemma54,
          max_return_gap_eps_0.5 = gap_max[1],
          max_return_gap_eps_0.1 = gap_max[2],
          max_return_gap_eps_0.01 = gap_max[3],
          bound_eps_0.01 = bound[3],
          time_fraction_eps_0.5 = fraction[1],
          time_fraction_eps_0.1 = fraction[2],
          time_fraction_eps_0.01 = fraction[3],
          window_growth_eps_0.1 = window_growth[2],
          min_distance_at_reentry = dist_at_reentry,
          alternation_alpha = alpha,
          alternation_last_failure_time = last_failure,
          alternation_times_checked_after = n_times_after,
          alternation_fraction_covered_after = frac_covered_after,
          sigma_gap_at_y1 = margin,
          r_equal_one_max_abs_error = err_r1))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the length of each visit to the neighbourhood of the equilibrium and the gap that
# follows it, against the index of the loop. The visits grow geometrically while the gaps stay
# under the bound M + log(1/eps)/lambda_1 + log(eps)/lambda_2, which is bounded return time.
# Right: the same information as a timeline of the last loops, where the visits are the bars.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    eps_fig <- 0.1
    m1 <- -log(eps_fig) / lam1
    m2 <- log(eps_fig) / lam2
    starts <- it$u + m2
    ends <- it$v - m1
    keep <- ends > starts
    starts <- starts[keep]; ends <- ends[keep]
    bound_fig <- M + m1 + m2
    visits <- data.frame(n = seq_along(starts), length = ends - starts, what = "Visit")
    gaps <- data.frame(n = seq_len(length(starts) - 1L),
                       length = starts[-1] - ends[-length(ends)], what = "Gap between visits")
    p1 <- ggplot(rbind(visits, gaps), aes(n, length, colour = what, shape = what)) +
        geom_hline(yintercept = bound_fig, colour = "#B8390E", linewidth = 0.4, linetype = "22") +
        geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
        annotate("text", x = 1, y = bound_fig * 1.35, hjust = 0, size = 2.6, colour = "#B8390E",
                 label = sprintf("Bound on the gap: %.3f", bound_fig)) +
        scale_colour_manual(values = c("Visit" = "#21918c", "Gap between visits" = "#440154"),
                            name = NULL) +
        scale_shape_manual(values = c(17, 16), name = NULL) +
        scale_y_log10() +
        labs(title = "Visits that grow and excursions that do not",
             subtitle = kb_unicode(sprintf("Neighbourhood of radius %.2f of the equilibrium, eigenvalues $\\lambda_1 = %.1f$ and $\\lambda_2 = %.1f$",
                                           eps_fig, lam1, lam2)),
             x = "Loop", y = kb_tex("Length of the interval, $\\log$ scale")) +
        kb_theme()
    # A zoom on one transition, since at the scale of the whole itinerary the excursion is
    # shorter than the width of a line: the visits last hundreds of time units and the gap
    # between them stays near four.
    j <- length(starts) - 1L
    gap_j <- starts[j + 1L] - ends[j]
    win <- c(ends[j] - 5 * gap_j, starts[j + 1L] + 5 * gap_j)
    bars <- data.frame(xmin = c(win[1], starts[j + 1L]), xmax = c(ends[j], win[2]),
                       label = c(sprintf("Visit %d, of length %.0f", j, ends[j] - starts[j]),
                                 sprintf("Visit %d, of length %.0f", j + 1L, ends[j + 1L] - starts[j + 1L])))
    p2 <- ggplot(bars) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.7, ymax = 1.3),
                  fill = "#21918c", colour = NA) +
        annotate("segment", x = ends[j], xend = starts[j + 1L], y = 1, yend = 1,
                 colour = "#440154", linewidth = 0.7,
                 arrow = grid::arrow(length = unit(0.06, "in"), ends = "both")) +
        annotate("text", x = (ends[j] + starts[j + 1L]) / 2, y = 1.45, size = 2.7,
                 colour = "#440154", label = sprintf("The excursion lasts %.2f", gap_j)) +
        annotate("text", x = win[1] + gap_j * 0.4, y = 1, size = 2.6, colour = "white", hjust = 0,
                 label = sprintf("Visit %d", j)) +
        annotate("text", x = win[2] - gap_j * 0.4, y = 1, size = 2.6, colour = "white", hjust = 1,
                 label = sprintf("Visit %d", j + 1L)) +
        coord_cartesian(xlim = win, ylim = c(0.55, 1.6)) +
        scale_y_continuous(breaks = NULL) +
        labs(title = "One transition, at the scale of the excursion",
             subtitle = sprintf("Between two visits that last %.0f and %.0f time units",
                                ends[j] - starts[j], ends[j + 1L] - starts[j + 1L]),
             x = kb_tex("Time $t$"), y = NULL) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.15, 1)) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "In the run recorded by checks/homoclinic-loop-return-times.R the longest gap at this radius was %.3f against the bound %.3f, the visits grew by a factor of %.1f over %d loops, and the share of time inside the neighbourhood was %.3f. The loop is therefore a bounded-return-time attractor and not a measure attractor, since only the loop itself converges to the equilibrium.",
                gap_max[2], bound[2], window_growth[2], length(starts), fraction[2])),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "homoclinic-loop-return-times", width = 9.4, height = 4.3))
}
