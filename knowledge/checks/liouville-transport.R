# ============================ #
# Transport of densities by a flow: the Liouville equation ####
# ============================ #
#
# Let v be a C1 vector field on an interval with flow Phi_t, and let rho0 be
# a C1 initial density. The density of the pushforward (Phi_t)#(rho0 m), with
# m the Lebesgue measure, is rho(t, y) = rho0(Phi_{-t}(y)) |d Phi_{-t}(y) / dy|,
# and it satisfies d rho / dt + d (rho v) / dy = 0. The script checks the
# equation for two flows with closed-form solutions.
#
# 1. Linear contraction v(x) = -a x, Gaussian rho0: the transported density is
#    Gaussian with mean m0 exp(-a t) and standard deviation s0 exp(-a t).
# 2. Logistic growth v(x) = x (1 - x) on (0, 1) with rho0(x) = 6 x (1 - x):
#    Phi_{-t}(y) = y exp(-t) / (1 - y + y exp(-t)), with derivative
#    exp(-t) / (1 - y + y exp(-t))^2.
#
# The residual is computed with symbolic derivatives from D() and, as an
# independent route, with central finite differences. Mass conservation and a
# Kolmogorov-Smirnov test on transported samples complete the check.

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

# Relative residual of d rho / dt + d (rho v) / dy on a grid, both routes.
residuals <- function(rho_expr, v_expr, env, t_values, y_values, h = 1e-5) {
    flux_expr <- substitute(R * V, list(R = rho_expr, V = v_expr))
    d_t <- D(rho_expr, "t")
    d_y <- D(flux_expr, "y")
    rho_fun  <- function(t, y) eval(rho_expr, c(list(t = t, y = y), env))
    flux_fun <- function(t, y) eval(flux_expr, c(list(t = t, y = y), env))
    sym <- fd <- scale <- 0
    for (t in t_values) {
        dt_sym <- eval(d_t, c(list(t = t, y = y_values), env))
        dy_sym <- eval(d_y, c(list(t = t, y = y_values), env))
        dt_fd  <- (rho_fun(t + h, y_values) - rho_fun(t - h, y_values)) / (2 * h)
        dy_fd  <- (flux_fun(t, y_values + h) - flux_fun(t, y_values - h)) / (2 * h)
        sym   <- max(sym, abs(dt_sym + dy_sym))
        fd    <- max(fd, abs(dt_fd + dy_fd))
        scale <- max(scale, abs(dt_sym))
    }
    # A stationary density (v = 0) has no time derivative to scale by, so the
    # absolute residual is returned in that case.
    if (scale == 0) scale <- 1
    c(symbolic = sym / scale, finite_difference = fd / scale)
}

# ============================ #
# 1. Linear contraction ####
# ============================ #
lin_env <- list(a = 0.7, m0 = 1.5, s0 = 0.8)
lin_rho <- quote(exp(-(y - m0 * exp(-a * t))^2 / (2 * (s0 * exp(-a * t))^2)) /
                     (sqrt(2 * pi) * s0 * exp(-a * t)))
lin_v <- quote(-a * y)
lin_t <- c(0.3, 1, 2)
lin_res <- residuals(lin_rho, lin_v, lin_env, lin_t, seq(-1, 3, length.out = 801))
lin_mass <- max(vapply(lin_t, function(t) abs(integrate(function(y)
    eval(lin_rho, c(list(t = t, y = y), lin_env)), -Inf, Inf)$value - 1), numeric(1)))
x0 <- rnorm(100000L, lin_env$m0, lin_env$s0)
xt <- x0 * exp(-lin_env$a * 1)
lin_ks <- suppressWarnings(ks.test(xt, "pnorm", lin_env$m0 * exp(-lin_env$a), lin_env$s0 * exp(-lin_env$a)))

# The degenerate field v = 0 leaves the density unchanged in time.
zero_env <- list(a = 0, m0 = 1.5, s0 = 0.8)
zero_res <- residuals(lin_rho, lin_v, zero_env, lin_t, seq(-1, 3, length.out = 101))

# ============================ #
# 2. Logistic growth ####
# ============================ #
log_env <- list()
log_rho <- quote(6 * (y * exp(-t) / (1 - y + y * exp(-t))) *
                     (1 - y * exp(-t) / (1 - y + y * exp(-t))) *
                     exp(-t) / (1 - y + y * exp(-t))^2)
log_v <- quote(y * (1 - y))
log_t <- c(0.25, 1, 2.5)
log_res <- residuals(log_rho, log_v, log_env, log_t, seq(0.001, 0.999, length.out = 999))
log_mass <- max(vapply(log_t, function(t) abs(integrate(function(y)
    eval(log_rho, list(t = t, y = y)), 0, 1, rel.tol = 1e-11)$value - 1), numeric(1)))
phi <- function(t, x) x * exp(t) / (1 - x + x * exp(t))
u0 <- rbeta(100000L, 2, 2)
ut <- phi(1, u0)
F_t <- function(y) { z <- phi(-1, y); 3 * z^2 - 2 * z^3 }
log_ks <- suppressWarnings(ks.test(ut, F_t))

pass <- lin_res[["symbolic"]] < 1e-10 && lin_res[["finite_difference"]] < 1e-5 &&
    log_res[["symbolic"]] < 1e-10 && log_res[["finite_difference"]] < 1e-5 &&
    zero_res[["symbolic"]] == 0 &&
    lin_mass < 1e-8 && log_mass < 1e-8 &&
    lin_ks$p.value > 1e-3 && log_ks$p.value > 1e-3

emit("liouville-transport", if (pass) "pass" else "fail",
     "Transported densities solve the continuity equation for a linear and a logistic flow",
     list(linear_symbolic_rel_residual = lin_res[["symbolic"]],
          linear_fd_rel_residual = lin_res[["finite_difference"]],
          logistic_symbolic_rel_residual = log_res[["symbolic"]],
          logistic_fd_rel_residual = log_res[["finite_difference"]],
          zero_field_symbolic_residual = zero_res[["symbolic"]],
          linear_mass_error = lin_mass, logistic_mass_error = log_mass,
          linear_ks_p = lin_ks$p.value, logistic_ks_p = log_ks$p.value))

# ============================ #
# Figure ####
# ============================ #
#
# Left: the density transported by the linear flow, which keeps its shape and narrows as the
# flow contracts. Right: the density transported by the logistic flow, whose shape changes
# because the flow is not linear. Both are the exact transported densities whose residuals in
# the continuity equation the check has just measured.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    eval_rho <- function(expr, env, t, ys)
        vapply(ys, function(y) eval(expr, c(env, list(t = t, y = y))), numeric(1))
    ys_lin <- seq(-1, 3, length.out = 500)
    times <- c(0, 0.3, 1, 2)
    lin_df <- do.call(rbind, lapply(times, function(t)
        data.frame(y = ys_lin, rho = eval_rho(lin_rho, lin_env, t, ys_lin),
                   t = factor(sprintf("t = %.1f", t), levels = sprintf("t = %.1f", times)))))
    p1 <- ggplot(lin_df, aes(y, rho, colour = t)) +
        geom_line(linewidth = 0.6) +
        scale_colour_viridis_d(option = "viridis", end = 0.88, name = NULL) +
        labs(title = "A density carried by a linear contraction",
             subtitle = kb_unicode(sprintf("Flow of $y' = -ay$ with $a = %.1f$; the Gaussian keeps its shape and narrows", lin_env$a)),
             x = kb_tex("State $y$"), y = kb_tex("Density $\\rho(t, y)$")) +
        kb_theme()
    ys_log <- seq(0.001, 0.999, length.out = 500)
    log_df <- do.call(rbind, lapply(times, function(t)
        data.frame(y = ys_log, rho = eval_rho(log_rho, log_env, t, ys_log),
                   t = factor(sprintf("t = %.1f", t), levels = sprintf("t = %.1f", times)))))
    p2 <- ggplot(log_df, aes(y, rho, colour = t)) +
        geom_line(linewidth = 0.6) +
        scale_colour_viridis_d(option = "magma", end = 0.82, name = NULL) +
        labs(title = "A density carried by the logistic flow",
             subtitle = kb_tex("Flow of $y' = y(1 - y)$; the shape changes, since the flow is not linear"),
             x = kb_tex("State $y$"), y = kb_tex("Density $\\rho(t, y)$")) +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, nrow = 1) +
            patchwork::plot_annotation(caption = kb_caption(sprintf(
                "Each curve is the density obtained by transporting the initial one along the flow, and each satisfies the continuity equation, which the run recorded by checks/liouville-transport.R verified by symbolic differentiation: the residual was at most %.1e for the linear flow and %.1e for the logistic flow, and the mass stayed one to %.1e.",
                max(lin_res), max(log_res), max(lin_mass))),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "liouville-transport", width = 9.2, height = 4.2))
}
