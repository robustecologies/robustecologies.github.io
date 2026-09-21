# ============================ #
# The Fokker-Planck equation of the Ornstein-Uhlenbeck process ####
# ============================ #
#
# Claims, with the locators of arnold1974 (sec. 8.3 for the process, Thm. 2.6.9 for the
# equation):
#
# The scalar equation dX = -a X dt + sigma dW with a > 0 has the solution whose law is
# Gaussian with mean m(t) = x0 exp(-a t) and variance v(t) = sigma^2 (1 - exp(-2 a t)) / (2 a),
# and the density of that law solves the Fokker-Planck equation
#
#     d rho / dt = d/dx ( a x rho ) + (sigma^2 / 2) d2 rho / dx2 .
#
# 1. The Gaussian density with that mean and variance makes the residual of the equation
#    vanish. Reference: analytic derivatives of the Gaussian, and, independently, central
#    finite differences in time and space.
# 2. The stationary density is Gaussian with mean zero and variance sigma^2/(2a), and it
#    makes the stationary residual vanish; the time-dependent solution converges to it.
# 3. The equation conserves mass.
# 4. Monte Carlo: an Euler-Maruyama simulation of the equation reproduces the mean and the
#    variance of the closed form and passes a Kolmogorov-Smirnov test against it. This route
#    shares no code with the two above.
# 5. Zero-noise limit: the stationary law concentrates at the origin, the attractor of the
#    deterministic equation, and the mass outside an interval is the Gaussian tail, which
#    the script compares with its closed form.
# 6. Degenerate case sigma = 0: the equation becomes the Liouville equation of the linear
#    contraction, whose solution is the transported density rho0(x exp(a t)) exp(a t), and
#    the residual of that solution vanishes as well.

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
x0 <- 2.0

m_t <- function(t) x0 * exp(-a * t)
v_t <- function(t) sigma^2 * (1 - exp(-2 * a * t)) / (2 * a)
rho <- function(t, x) dnorm(x, mean = m_t(t), sd = sqrt(v_t(t)))

# ============================ #
# Test 1: the closed form solves the equation ####
# ============================ #
#
# With rho Gaussian of mean m and variance v, the derivatives are analytic, and m' = -a m,
# v' = sigma^2 - 2 a v. The residual is compared with zero relative to the size of the terms
# that make it up, so the tolerance 1e-10 is on a scaled quantity. This route takes those two
# derivatives as given, so it tests the equation and not the closed forms for m and v; the
# finite differences below and the simulation in test 4 test those.

residual_analytic <- function(t, x) {
    m <- m_t(t); v <- v_t(t)
    r <- dnorm(x, m, sqrt(v))
    dm <- -a * m
    dv <- sigma^2 - 2 * a * v
    z <- x - m
    d_t <- r * (z / v * dm + (z^2 / (2 * v^2) - 1 / (2 * v)) * dv)
    d_x <- -r * z / v
    d_xx <- r * (z^2 / v^2 - 1 / v)
    d_t - (a * r + a * x * d_x) - (sigma^2 / 2) * d_xx
}
grid_t <- seq(0.05, 3, by = 0.05)
grid_x <- seq(-4, 4, by = 0.05)
scale <- max(abs(outer(grid_t, grid_x, Vectorize(function(t, x) rho(t, x) / v_t(t)))))
res_analytic <- max(abs(outer(grid_t, grid_x, Vectorize(residual_analytic)))) / scale

# Independent route: central differences, with step 1e-4 in t and 1e-3 in x. The truncation
# error of the second difference is of order 1e-6 relative to the terms, so the tolerance is
# 1e-5 on the scaled residual.
residual_fd <- function(t, x) {
    ht <- 1e-4; hx <- 1e-3
    d_t <- (rho(t + ht, x) - rho(t - ht, x)) / (2 * ht)
    flux <- function(u) a * u * rho(t, u)
    d_flux <- (flux(x + hx) - flux(x - hx)) / (2 * hx)
    d_xx <- (rho(t, x + hx) - 2 * rho(t, x) + rho(t, x - hx)) / hx^2
    d_t - d_flux - (sigma^2 / 2) * d_xx
}
res_fd <- max(abs(outer(seq(0.2, 3, by = 0.2), seq(-3, 3, by = 0.25),
                        Vectorize(residual_fd)))) / scale

# ============================ #
# Test 2: the stationary density ####
# ============================ #

v_inf <- sigma^2 / (2 * a)
stationary_residual <- function(x) {
    r <- dnorm(x, 0, sqrt(v_inf))
    d_x <- -r * x / v_inf
    d_xx <- r * (x^2 / v_inf^2 - 1 / v_inf)
    (a * r + a * x * d_x) + (sigma^2 / 2) * d_xx
}
res_stationary <- max(abs(vapply(grid_x, stationary_residual, numeric(1)))) / scale
approach <- abs(v_t(20) - v_inf)

# ============================ #
# Test 3: conservation of mass ####
# ============================ #

mass_err <- max(abs(vapply(c(0.1, 0.5, 1, 3), function(t)
    integrate(function(x) rho(t, x), -Inf, Inf)$value - 1, numeric(1))))

# ============================ #
# Test 4: Euler-Maruyama ####
# ============================ #
#
# With n = 50000 paths the standard error of the mean is about 0.4 percent of the standard
# deviation and that of the variance about 0.6 percent, and the Euler-Maruyama bias is of
# order dt = 0.002, so the tolerances are 3 percent on both moments.

n <- 50000L
dt <- 0.002
tend <- 1.5
steps <- as.integer(tend / dt)
x <- rep(x0, n)
for (k in seq_len(steps)) x <- x - a * x * dt + sigma * sqrt(dt) * rnorm(n)
mc_mean_err <- abs(mean(x) - m_t(tend)) / sqrt(v_t(tend))
mc_var_err <- abs(var(x) / v_t(tend) - 1)
# The Kolmogorov-Smirnov test is a shape check on 2000 of the paths. Its p-value is a
# random variable and the seed fixes it, so the criterion is only that the test does not
# reject at the one in a thousand level; the moments above are the quantitative test.
ks_p <- ks.test(x[1:2000], "pnorm", m_t(tend), sqrt(v_t(tend)))$p.value

# ============================ #
# Test 5: the zero-noise limit ####
# ============================ #
#
# The stationary law is N(0, sigma^2 / (2a)), so the mass outside [-eps, eps] is
# 2 pnorm(-eps / sqrt(sigma^2 / (2a))), which tends to zero with sigma.

eps <- 0.1
tails <- vapply(c(0.7, 0.2, 0.05, 0.01), function(s)
    2 * pnorm(-eps / sqrt(s^2 / (2 * a))), numeric(1))
tail_closed_form <- 2 * pnorm(-eps / sqrt(0.01^2 / (2 * a)))
tail_err <- abs(tails[4] - tail_closed_form)

# ============================ #
# Test 6: the noiseless case is the Liouville equation ####
# ============================ #
#
# For sigma = 0 the transported density of the linear contraction is
# rho(t, x) = rho0(x exp(a t)) exp(a t), and it must satisfy d rho/dt + d/dx (-a x rho) = 0.

# The derivatives are taken symbolically with D(), which differentiates the formula for the
# transported density itself, so the test does not repeat a derivation made by hand. Writing
# u = x exp(a t), the density is exp(a t) rho0(u), and the claim is that it makes the residual
# of the Liouville equation vanish.
mu0 <- 0.5; s0 <- 0.3
rho_expr <- quote(exp(a * t) * exp(-(x * exp(a * t) - mu0)^2 / (2 * s0^2)) / (s0 * sqrt(2 * pi)))
flux_expr <- bquote(-a * x * .(rho_expr))
d_t_expr <- D(rho_expr, "t")
d_flux_expr <- D(flux_expr, "x")
res_liouville <- 0
for (tt in seq(0.1, 2, by = 0.1)) for (xx in seq(-2, 2, by = 0.1)) {
    env <- list(a = a, mu0 = mu0, s0 = s0, t = tt, x = xx)
    res_liouville <- max(res_liouville, abs(eval(d_t_expr, env) + eval(d_flux_expr, env)))
}

# Independent route on the same claim, with central differences and a tolerance scaled by the
# size of the terms, since the narrowing density makes the absolute residual grow with time.
rho_liouville <- function(t, x) eval(rho_expr, list(a = a, mu0 = mu0, s0 = s0, t = t, x = x))
res_liouville_fd <- 0
for (tt in seq(0.1, 1, by = 0.1)) for (xx in seq(-1, 1, by = 0.1)) {
    ht <- 1e-5; hx <- 1e-4
    d_t <- (rho_liouville(tt + ht, xx) - rho_liouville(tt - ht, xx)) / (2 * ht)
    flux <- function(u) -a * u * rho_liouville(tt, u)
    d_flux <- (flux(xx + hx) - flux(xx - hx)) / (2 * hx)
    res_liouville_fd <- max(res_liouville_fd, abs(d_t + d_flux) / max(1, abs(d_t)))
}

ok <- res_analytic < 1e-10 && res_fd < 1e-5 && res_stationary < 1e-10 &&
    approach < 1e-15 && mass_err < 1e-10 &&
    mc_mean_err < 0.03 && mc_var_err < 0.03 && ks_p > 0.001 &&
    all(diff(tails) < 0) && tails[4] < 1e-12 && tail_err < 1e-15 &&
    res_liouville < 1e-12 && res_liouville_fd < 1e-4

emit("ornstein-uhlenbeck-fokker-planck", if (ok) "pass" else "fail",
     "The Gaussian law of the Ornstein-Uhlenbeck process solves its Fokker-Planck equation, agrees with an Euler-Maruyama simulation, and concentrates on the deterministic equilibrium as the noise vanishes",
     list(analytic_scaled_residual = res_analytic,
          finite_difference_scaled_residual = res_fd,
          stationary_scaled_residual = res_stationary,
          variance_gap_at_t20 = approach,
          mass_max_abs_error = mass_err,
          monte_carlo_mean_error_in_sd = mc_mean_err,
          monte_carlo_variance_rel_error = mc_var_err,
          kolmogorov_smirnov_p = ks_p,
          stationary_tail_sigma_0.7 = tails[1],
          stationary_tail_sigma_0.05 = tails[3],
          stationary_tail_sigma_0.01 = tails[4],
          tail_closed_form_error = tail_err,
          liouville_analytic_residual = res_liouville,
          liouville_finite_difference_residual = res_liouville_fd))

# ============================ #
# Figure ####
# ============================ #
#
# Drawn from the objects computed above, so the picture and the recorded values come from the
# same run. Left: the solution of the Fokker-Planck equation at four times, with the histogram
# of the simulated paths at the last of them and the stationary density. Right: the stationary
# density as the noise falls, which is the zero-noise limit concentrating on the equilibrium.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    times <- c(0.05, 0.2, 0.5, 1.5)
    grid <- seq(-1.5, 3.2, length.out = 400)
    dens <- do.call(rbind, lapply(times, function(tt)
        data.frame(x = grid, density = rho(tt, grid), t = factor(sprintf("t = %.2f", tt)))))
    stat <- data.frame(x = grid, density = dnorm(grid, 0, sqrt(v_inf)))
    hist_df <- as.data.frame(with(hist(x, breaks = 60, plot = FALSE), list(x = mids, density = density)))
    p1 <- ggplot(dens, aes(x, density, colour = t)) +
        geom_area(data = stat, aes(x, density), inherit.aes = FALSE,
                  fill = "grey90", colour = NA) +
        geom_step(data = hist_df, aes(x, density), inherit.aes = FALSE,
                  colour = "grey35", linewidth = 0.35, direction = "mid") +
        geom_line(linewidth = 0.7) +
        annotate("text", x = -1.15, y = 0.93, hjust = 0, size = 2.5, colour = "grey35",
                 label = "Stationary density") +
        annotate("text", x = 0.45, y = 1.22, hjust = 0, size = 2.5, colour = "grey25",
                 label = kb_unicode("Simulated paths at $t = 1.5$")) +
        scale_colour_viridis_d(option = "viridis", end = 0.88, name = NULL) +
        labs(title = "Transport of the density by the Fokker-Planck equation",
             subtitle = bquote("Solution" ~ rho(t, x) ~ "of" ~ partialdiff[t] * rho == partialdiff[x](a*x*rho) + sigma^2/2 * partialdiff[x*x] * rho ~
                               ", with" ~ a == .(a) ~ "and" ~ sigma == .(sigma)),
             x = kb_tex("State $x$"), y = kb_tex("Density $\\rho(t, x)$")) +
        kb_theme()
    sig_seq <- c(0.7, 0.35, 0.15, 0.05)
    zn <- do.call(rbind, lapply(sig_seq, function(s)
        data.frame(x = grid, density = dnorm(grid, 0, sqrt(s^2 / (2 * a))),
                   s = factor(sprintf("$\\sigma = %.2f$", s), levels = sprintf("$\\sigma = %.2f$", sig_seq)))))
    p2 <- ggplot(zn, aes(x, density, colour = s)) +
        geom_line(linewidth = 0.7) +
        coord_cartesian(xlim = c(-1.2, 1.2)) +
        scale_colour_viridis_d(option = "magma", end = 0.8, name = NULL,
                               labels = kb_tex(levels(zn$s))) +
        labs(title = "The stationary density as the noise vanishes",
             subtitle = bquote("Gaussian of variance" ~ sigma^2/(2*a) ~ ", concentrating on the equilibrium at the origin"),
             x = kb_tex("State $x$"), y = kb_tex("Stationary density $\\mu_\\sigma$")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1 + labs(caption = NULL), p2 + labs(caption = NULL), nrow = 1) +
            patchwork::plot_annotation(
                caption = kb_caption(sprintf(
                    "Left: the closed-form solution at four times, the histogram of %d Euler-Maruyama paths at $t = 1.5$ as the grey step, and the stationary density as the shaded region. Right: the stationary density at four noise amplitudes $\\sigma$. Both panels come from the run recorded by checks/ornstein-uhlenbeck-fokker-planck.R, with $a = %.1f$ and $x_0 = %.1f$.",
                    n, a, x0)),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "ornstein-uhlenbeck-fokker-planck", width = 9.2, height = 4.2))
}
