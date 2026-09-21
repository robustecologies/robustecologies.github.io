# ============================ #
# The random attractor of the Ornstein-Uhlenbeck process ####
# ============================ #
#
# Claim, against the definitions of crauel1994 (Def. 2.1, Def. 3.9 and Thm. 3.11):
#
# The scalar equation dX = -a X dt + sigma dW with a > 0 generates a random dynamical system
# whose global random attractor is a single random point,
#
#     A(omega) = { Z(omega) },      Z(omega) = sigma * int_{-infinity}^{0} e^{a s} dW_s ,
#
# invariant in the sense that the cocycle carries the fibre over omega onto the fibre over the
# shifted noise, and pullback attracting: carrying a bounded set forward from earlier and
# earlier noise states collapses it onto that point. The statement is worked out here for this
# equation and is not quoted from the source, which treats the general theory and uses the
# Ornstein-Uhlenbeck process only as an auxiliary device for its infinite dimensional examples.
#
# The solution of the linear equation is explicit,
#
#     phi(t, omega) x = x e^{-a t} + sigma * int_0^t e^{-a (t - s)} dW_s ,
#
# so the pullback image of x from the noise state shifted back by t is
#
#     phi(t, theta_{-t} omega) x = x e^{-a t} + sigma * int_{-t}^{0} e^{a s} dW_s ,
#
# which gives every claim below by inspection. The script tests them against a simulation that
# shares no formula with the solution: an Euler-Maruyama integration driven by the same
# increments of the same noise path.
#
# 1. Contraction: two initial conditions under the same noise separate by exactly the factor
#    e^{-a t}, so the diameter of the image of any bounded set decays at that rate and the
#    attractor cannot contain two points. The simulation reproduces the factor of its own
#    scheme exactly and the continuous-time factor to the order of the step.
# 2. Pullback attraction: the images of a wide interval of initial conditions, carried forward
#    from time -t, converge to one point, which agrees with the stochastic integral.
# 3. Forward evolution does not converge to a point of the state space: the same interval
#    collapses, but the point it collapses onto keeps moving, so the attractor is a random
#    point and not a fixed one. The script records the spread of that point over time.
# 4. Invariance: evolving Z(omega) forward for a time t with the same noise gives the value
#    that the stochastic integral assigns to the shifted noise state, phi(t, omega) Z(omega) =
#    Z(theta_t omega).
# 5. The law of Z is the stationary law of the equation, Gaussian with mean zero and variance
#    sigma^2 / (2 a), which is the invariant measure that the Fokker-Planck equation of
#    checks/ornstein-uhlenbeck-fokker-planck.R has as its stationary density, so the fibre
#    measures delta_{Z(omega)} average to it.
# 6. Zero-noise limit: the variance of Z is sigma^2 / (2 a), so the random attractor collapses
#    onto the deterministic equilibrium as the noise vanishes.

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

a <- 1.3
sigma <- 0.7
dt <- 0.001
Tpast <- 40           # the noise is drawn on [-Tpast, Tfut]
Tfut <- 10

# One realisation of the noise, shared by every test below. The grid carries the increments of
# the Wiener process; the shift theta_t acts on this array by moving the origin.
n_past <- as.integer(Tpast / dt); n_fut <- as.integer(Tfut / dt)
dW <- rnorm(n_past + n_fut, 0, sqrt(dt))
t_grid <- seq(-Tpast, Tfut - dt, by = dt)        # left endpoint of each increment
idx0 <- n_past + 1L                              # first index with t >= 0

# The stochastic integral of the exact solution, as a left endpoint sum, which is the Ito sum.
# Z(theta_r omega) = sigma * int_{-infinity}^{r} e^{a (s - r)} dW_s.
Z_at <- function(r) {
    k <- which(t_grid < r)
    sigma * sum(exp(a * (t_grid[k] - r)) * dW[k])
}
Z0 <- Z_at(0)
# The part of the integral that the finite past leaves out has standard deviation
# sigma e^{-a Tpast} / sqrt(2 a), which is the accuracy of Z0 as a value of the true integral.
truncation_sd <- sigma * exp(-a * Tpast) / sqrt(2 * a)

# Euler-Maruyama from time s to time r along the same increments, started at x.
em <- function(x, s, r) {
    k <- which(t_grid >= s & t_grid < r)
    for (j in k) x <- x - a * x * dt + sigma * dW[j]
    x
}

# ============================ #
# Test 1: contraction under a shared noise ####
# ============================ #
#
# The difference of two solutions solves the noiseless equation, so it is multiplied by
# e^{-a t}; the scheme multiplies it by (1 - a dt)^n, which the simulation must reproduce
# exactly. Two limits decide the tolerances, and both are computed rather than chosen. The
# comparison with the scheme is limited by cancellation: the two trajectories are of order one
# and their difference decays to 1e-5, so the relative error of that difference is of order
# eps * max|X| / |x1 - x2| (1 - a dt)^n, and the assertion allows twenty times that bound. The
# comparison with the continuous factor carries the bias of the scheme itself,
# n |log(1 - a dt) + a dt| = a^2 t dt / 2 in relative terms.

t_c <- 10
n_c <- as.integer(t_c / dt)
X1 <- em(3, 0, t_c); X2 <- em(-2, 0, t_c)
d_sim <- abs(X1 - X2) / 5
d_scheme <- (1 - a * dt)^n_c
d_exact <- exp(-a * t_c)
err_scheme <- abs(d_sim - d_scheme) / d_scheme
err_exact <- abs(d_sim - d_exact) / d_exact
cancellation_bound <- 20 * .Machine$double.eps * max(abs(c(X1, X2))) / (5 * d_scheme)
pred_bias <- a^2 * t_c * dt / 2

# ============================ #
# Test 2: pullback attraction to one point ####
# ============================ #
#
# Carrying a wide interval forward from time -t, the images must collapse onto Z(omega). The
# spread of the images is the diameter of the interval times e^{-a t}, and the distance of the
# common limit from Z0 is of the order of the step of the scheme.

x_set <- c(-50, -10, 0, 7, 50)
pull <- function(t) vapply(x_set, function(x) em(x, -t, 0), numeric(1))
t_set <- c(2, 5, 10, 20)
pulls <- lapply(t_set, pull)
spreads <- vapply(pulls, function(v) diff(range(v)), numeric(1))
# Two floating-point effects bound this comparison and both are computed. Rounding accumulates
# over the n steps of each run, giving a relative error of order n * eps in the images, and the
# spread of the images cancels their common part, giving a further eps * max|x| / spread, which
# at t = 20 is the larger of the two by four orders. The assertion allows a factor of twenty on
# the first effect and five on the second.
spread_bound <- vapply(seq_along(t_set), function(i)
    20 * .Machine$double.eps * max(abs(pulls[[i]])) / spreads[i] +
    5 * (t_set[i] / dt) * .Machine$double.eps, numeric(1))
# The scheme contracts by (1 - a dt)^n, which the images must follow to the cancellation limit,
# and that factor differs from e^{-a t} by the bias a^2 t dt / 2 computed above.
pred_scheme <- diff(range(x_set)) * (1 - a * dt)^as.integer(t_set / dt)
pred_cont <- diff(range(x_set)) * exp(-a * t_set)
spread_rel_each <- abs(spreads - pred_scheme) / pred_scheme
spread_rel <- max(spread_rel_each / spread_bound)      # 1 means at the cancellation limit
spread_rel_cont <- max(abs(spreads - pred_cont) / pred_cont)
pred_bias_20 <- a^2 * max(t_set) * dt / 2
limit_gap <- max(abs(pull(20) - Z0))

# ============================ #
# Test 3: the point of the attractor moves ####
# ============================ #
#
# The same interval carried forward from time 0 also collapses, but onto a point that keeps
# moving with the noise. The spread of Z(theta_r omega) over r is of the order of its own
# standard deviation, sigma / sqrt(2 a), and is not a decaying quantity.

r_grid <- seq(0, Tfut - 0.01, by = 0.05)
Z_path <- vapply(r_grid, Z_at, numeric(1))
stationary_sd <- sigma / sqrt(2 * a)
path_sd <- sd(Z_path)
path_range <- diff(range(Z_path))

# ============================ #
# Test 4: invariance of the attractor ####
# ============================ #
#
# phi(t, omega) Z(omega) = Z(theta_t omega). The left side is the simulation started at Z0 and
# run to t; the right side is the stochastic integral ending at t. They differ by the bias of
# the scheme, of order dt, so the tolerance is 1e-2 of the standard deviation of Z.

inv_err <- vapply(c(1, 3, 7), function(t) abs(em(Z0, 0, t) - Z_at(t)), numeric(1))
inv_rel <- max(inv_err) / stationary_sd

# ============================ #
# Test 5: the law of the random point ####
# ============================ #
#
# Over independent realisations Z is Gaussian with mean zero and variance sigma^2 / (2 a). With
# m = 20000 realisations the standard error of the mean is stationary_sd / sqrt(m) = 3.1e-3 and
# that of the variance about 1 percent, so the tolerances are 4 standard errors on the mean and
# 5 percent on the variance. The integral is truncated at Tpast, which adds nothing at this
# accuracy.

m <- 20000L
n_z <- as.integer(20 / dt)                      # 20 time units of past is e^{-26} of the tail
wts <- exp(a * (seq(-20, -dt, by = dt)))
Zs <- vapply(seq_len(m), function(i) sigma * sum(wts * rnorm(n_z, 0, sqrt(dt))), numeric(1))
z_mean <- mean(Zs); z_var <- var(Zs)
se_mean <- stationary_sd / sqrt(m)
ks_p <- ks.test(Zs[1:2000], "pnorm", 0, stationary_sd)$p.value

# ============================ #
# Test 6: the zero-noise limit ####
# ============================ #
#
# The attractor is a point whose standard deviation is sigma / sqrt(2 a), so it collapses onto
# the deterministic equilibrium at the origin as sigma tends to zero. The check is the closed
# form itself, evaluated at three amplitudes against a simulation of the same integral.

zn <- vapply(c(0.7, 0.07, 0.007), function(s) {
    v <- vapply(seq_len(2000L), function(i) s * sum(wts * rnorm(n_z, 0, sqrt(dt))), numeric(1))
    sd(v)
}, numeric(1))
zn_pred <- c(0.7, 0.07, 0.007) / sqrt(2 * a)
zn_rel <- max(abs(zn - zn_pred) / zn_pred)

ok <- err_scheme < cancellation_bound && err_exact < 1.2 * pred_bias &&
    spread_rel < 1 && spread_rel_cont < 1.2 * pred_bias_20 && limit_gap < 5e-3 &&
    path_sd > 0.2 * stationary_sd && path_range > stationary_sd &&
    inv_rel < 1e-2 &&
    abs(z_mean) < 4 * se_mean && abs(z_var / stationary_sd^2 - 1) < 0.05 && ks_p > 0.001 &&
    zn_rel < 0.05 && truncation_sd < 1e-12

emit("ornstein-uhlenbeck-random-attractor", if (ok) "pass" else "fail",
     "The random attractor of the Ornstein-Uhlenbeck equation is one moving point: a bounded set carried forward from earlier noise collapses onto it, the cocycle maps it to its own shift, and its law is the stationary law of the equation",
     list(contraction_scheme_rel_error = err_scheme,
          contraction_continuous_rel_error = err_exact,
          predicted_scheme_bias = pred_bias,
          cancellation_bound = cancellation_bound,
          pullback_spread_error_in_units_of_the_cancellation_bound = spread_rel,
          pullback_spread_rel_error_at_t20 = spread_rel_each[4],
          cancellation_bound_at_t20 = spread_bound[4],
          pullback_spread_rel_error_continuous = spread_rel_cont,
          predicted_bias_at_t20 = pred_bias_20,
          pullback_spread_at_t20 = spreads[4],
          distance_of_the_limit_from_the_integral = limit_gap,
          moving_point_sd = path_sd,
          moving_point_range = path_range,
          stationary_sd = stationary_sd,
          invariance_error_in_sd = inv_rel,
          law_mean = z_mean,
          law_mean_standard_error = se_mean,
          law_variance_rel_error = abs(z_var / stationary_sd^2 - 1),
          kolmogorov_smirnov_p = ks_p,
          zero_noise_sd_rel_error = zn_rel,
          truncation_sd = truncation_sd))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the pullback picture. Five initial conditions spread over an interval of length 100 are
# carried forward from time -20 under one realisation of the noise, and they collapse onto the
# single point Z(omega) that travels with that realisation. Right: the law of that point over
# independent realisations, against the stationary law of the equation.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    # The whole path of the scheme, rather than its endpoint, for the five initial conditions.
    em_path <- function(x, s, r) {
        k <- which(t_grid >= s & t_grid < r)
        out <- numeric(length(k)); i <- 0L
        for (j in k) { x <- x - a * x * dt + sigma * dW[j]; i <- i + 1L; out[i] <- x }
        data.frame(t = t_grid[k] + dt, x = out)
    }
    keep <- function(d) d[seq(1, nrow(d), by = 40L), ]      # 25 points per time unit, which keeps the file small
    paths <- do.call(rbind, lapply(seq_along(x_set), function(i) {
        d <- keep(em_path(x_set[i], -20, 0)); d$start <- factor(x_set[i]); d
    }))
    r_fine <- seq(-20, 0, by = 0.1)
    attractor <- data.frame(t = r_fine, x = vapply(r_fine, Z_at, numeric(1)))
    p1 <- ggplot(paths, aes(t, x, colour = start)) +
        geom_line(linewidth = 0.4, alpha = 0.9) +
        geom_line(data = attractor, aes(t, x), inherit.aes = FALSE,
                  colour = "black", linewidth = 0.8) +
        annotate("text", x = -19.3, y = 40, hjust = 0, size = 2.6, colour = "grey25",
                 label = "Five initial conditions spread over an interval of length 100") +
        annotate("text", x = -19.3, y = -46, hjust = 0, size = 2.6, colour = "black",
                 label = "Black: the attractor, one point moving with the noise") +
        scale_colour_viridis_d(option = "viridis", end = 0.9, name = "Initial condition") +
        labs(title = "Pullback attraction to a single random point",
             subtitle = bquote("Solutions of" ~ dX == -a*X*dt + sigma*dW ~ "under one realisation, with" ~
                               a == .(a) ~ "and" ~ sigma == .(sigma)),
             x = kb_tex("Time $t$"), y = kb_tex("State $X_t$")) +
        kb_theme()
    p2 <- ggplot(data.frame(z = Zs), aes(z)) +
        geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey80", colour = "grey55",
                       linewidth = 0.2) +
        stat_function(fun = dnorm, args = list(mean = 0, sd = stationary_sd),
                      colour = "#440154", linewidth = 0.8) +
        labs(title = "The law of the point of the attractor",
             subtitle = bquote("Histogram of" ~ .(m) ~ "realisations against the stationary law" ~
                               N(0, sigma^2/(2*a))),
             x = kb_tex("Position $Z$ of the attractor"), y = kb_tex("Density")) +
        kb_theme()
    # At the scale of the initial conditions the attractor is a flat line, so an inset repeats
    # the last five time units at the scale of the attractor itself, where the five images are
    # no longer distinguishable from the point they have collapsed onto.
    zoom <- paths[paths$t > -5, ]
    p_inset <- ggplot(zoom, aes(t, x, colour = start)) +
        geom_line(linewidth = 0.9, alpha = 0.9) +
        geom_line(data = attractor[attractor$t > -5, ], aes(t, x), inherit.aes = FALSE,
                  colour = "black", linewidth = 0.35) +
        scale_colour_viridis_d(option = "viridis", end = 0.9, guide = "none") +
        labs(x = NULL, y = NULL, title = "The last five time units, at the scale of the attractor") +
        kb_theme(base_size = 7) +
        theme(plot.title = element_text(size = 6.5, colour = "grey25"),
              plot.background = element_rect(fill = "white", colour = "grey70"))
    p1 <- p1 + patchwork::inset_element(p_inset, left = 0.30, bottom = 0.06, right = 0.99, top = 0.47)
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1.3, 1)) +
            patchwork::plot_annotation(
                caption = kb_caption(sprintf(
                    "Left: under one fixed realisation of the noise the five images arrive within %.1e of one another at $t = 0$, the contraction $e^{-at}$ predicted for the scheme, and within %.1e of the value of the stochastic integral. Right: over %d independent realisations the point has mean %.4f against a standard error of %.4f and variance within %.1f percent of the closed form. Both panels come from the run recorded by checks/ornstein-uhlenbeck-random-attractor.R.",
                    spreads[4], limit_gap, m, z_mean, se_mean, 100 * abs(z_var / stationary_sd^2 - 1))),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "ornstein-uhlenbeck-random-attractor", width = 9.6, height = 4.3))
}
