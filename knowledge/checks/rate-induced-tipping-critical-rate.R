# ============================ #
# Critical rates for a fold with a moving parameter ####
# ============================ #
#
# The note on rate-induced tipping states that a system can fail to keep up with an attractor that
# moves although the frozen system keeps that attractor at every parameter value. This check builds
# the smallest system in which that happens, computes its critical rate and compares it with a
# closed form derived here, so the claim is tested and not restated.
#
# The system is dx/dt = (x - L(t))^2 - 1. For a parameter frozen at L the equilibria are x = L - 1,
# which is stable since the derivative of the right-hand side is 2(x - L) = -2 there, and x = L + 1,
# which is unstable. No value of L destroys either, so the frozen system has no bifurcation and any
# tipping is rate-induced.
#
# 1. Analytic critical rate for a linear ramp. With L(t) = -r t the substitution u = x - L(t) gives
#    du/dt = u^2 - 1 + r, whose equilibria u = -sqrt(1 - r) and u = +sqrt(1 - r) exist if and only
#    if r <= 1 and collide at r = 1. The state therefore tracks for r < 1 and escapes for r > 1, so
#    the critical rate is exactly 1. The script finds it by bisection on the integrated solution and
#    compares with that value. The moving-frame equation is derived by hand, so the script does not
#    use it: it integrates the original equation in the original variable.
# 2. Quasi-static lag. For small r the tracked state sits at u = -sqrt(1 - r) = -1 + r/2 + O(r^2),
#    so the lag behind the frozen equilibrium is r/2 to first order. The script measures the lag at
#    three small rates and compares the ratio lag / r with 1/2.
# 3. Integrator. The solver is validated against a closed form on a problem with the same shape and
#    a known solution: dx/dt = -(x - L(t)) with L(t) = -r t has x(t) = -r t + r + (x_0 - r) e^{-t},
#    which the script compares pointwise. This tests the integration and not the model.
# 4. Asymptotically autonomous shift, which is solvable exactly. With L(t) = -(D/2)(1 + tanh(rt/2))
#    the parameter moves from 0 to -D and the system is autonomous at both ends, which is the
#    setting of the note. In the moving frame the equation is the Riccati equation
#    du/dt = u^2 - 1 + (D r / 4) sech^2(r t / 2), and the substitution u = -psi'/psi turns it into
#    psi_ss + (D / r) sech^2(s) psi = (4 / r^2) psi with s = r t / 2, the Poschl-Teller eigenvalue
#    problem. Tipping is the appearance of a zero of psi, and the threshold is where the decay rate
#    2 / r meets the index lambda fixed by lambda (lambda + 1) = D / r. Solving
#    (2/r)(2/r + 1) = D/r gives the critical rate in closed form,
#        r_c = 4 / (D - 2),
#    so no rate tips a shift of magnitude D <= 2, and the critical rate falls as 4/D for large
#    shifts. The derivation is made here and is not taken from a source; the script tests it by
#    bisecting for r_c at D = 4, 8 and 16, and by checking that a shift of magnitude 1.5 does not
#    tip at rates up to 50.
# 5. End-point tracking. Below the critical rate the solution must arrive at the future equilibrium
#    -D - 1, and above it must escape. The script checks both at r_c times 0.9 and 1.1.
#
# Tolerances. The bisection is run to a width of 1e-6 in r, so the critical rates carry that
# uncertainty and are compared with the closed form at 1e-5. The integrator is compared with the
# closed form at 1e-8, above its own rtol and atol of 1e-10 and far below any difference that would
# matter. The quasi-static ratio is compared at 5e-3, the size of the first neglected term r/8 at
# the largest rate used. Escape is declared when the distance to the moving equilibrium exceeds
# 1e3, which the fold reaches in finite time, and it is measured in that moving frame because the
# state follows the parameter to large values without tipping. Tracking is declared when the final
# state is within 1e-6 of the future equilibrium. The horizon bias of test 1 is compared with pi^2
# at 0.2, two per cent of it, which is above the term of order sqrt(e) that the expansion drops.

set.seed(20260918L)
suppressPackageStartupMessages(library(deSolve))

emit <- function(id, status, summary, metrics) {
    val <- function(v) {
        if (is.character(v)) return(sprintf("\"%s\"", v))
        if (is.logical(v)) return(if (v) "true" else "false")
        if (!is.finite(v)) return("null")
        formatC(v, digits = 6, format = "g")
    }
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

# ---- the shifts and the model ----
ramp <- function(t, r) -r * t
shift <- function(t, r, D) -(D / 2) * (1 + tanh(r * t / 2))
fold_rhs <- function(t, state, parms) {
    L <- parms$path(t, parms)
    list(c(x = (state[["x"]] - L)^2 - 1))
}
integrate_fold <- function(path, times, x0, ..., bound = 1e3) {
    parms <- c(list(path = function(t, p) path(t, ...)), list(...))
    # The root is on the distance to the moving equilibrium and not on the state itself, since the
    # state follows the parameter to large values without tipping.
    out <- ode(y = c(x = x0), times = times, func = fold_rhs, parms = parms,
               method = "lsoda", rtol = 1e-10, atol = 1e-10,
               rootfun = function(t, state, parms) abs(state[["x"]] - parms$path(t, parms)) - bound)
    sol <- as.data.frame(out)
    sol$u <- sol$x - vapply(sol$time, function(t) parms$path(t, parms), numeric(1))
    sol
}
escaped <- function(sol, bound = 1e3) any(!is.finite(sol$u)) || max(sol$u) >= bound * 0.99

critical_rate <- function(path, times, x0_of_r, lo, hi, tol = 1e-6, ...) {
    while (hi - lo > tol) {
        mid <- (lo + hi) / 2
        sol <- integrate_fold(path, times, x0_of_r(mid), r = mid, ...)
        if (escaped(sol)) hi <- mid else lo <- mid
    }
    (lo + hi) / 2
}

# ---- test 3 first: the solver against a closed form ----
linear_rhs <- function(t, state, parms) list(c(x = -(state[["x"]] - ramp(t, parms$r))))
t_lin <- seq(0, 12, by = 0.05); r_lin <- 0.3; x0_lin <- 0.7
num_lin <- as.data.frame(ode(y = c(x = x0_lin), times = t_lin, func = linear_rhs,
                             parms = list(r = r_lin), method = "lsoda", rtol = 1e-10, atol = 1e-10))
exact_lin <- -r_lin * t_lin + r_lin + (x0_lin - r_lin) * exp(-t_lin)
err_solver <- max(abs(num_lin$x - exact_lin))

# ---- test 1: the critical rate of the linear ramp, and the bias a finite horizon imposes ----
# From u(0) = -1 the escape time at r = 1 + e is the integral of du / (u^2 + e) from -1 to infinity,
# which is pi / sqrt(e) - 1 + O(sqrt(e)). A run of length T cannot see a rate below
# 1 + pi^2 / (T + 1)^2, so the measured critical rate carries that bias and must follow it.
horizons <- c(100, 200, 400)
r_by_horizon <- vapply(horizons, function(T)
    critical_rate(ramp, c(0, seq(0.5, T, length.out = 400)), function(r) -1, lo = 0.9, hi = 1.5),
    numeric(1))
bias_constant <- (r_by_horizon - 1) * (horizons + 1)^2        # expected to approach pi^2 = 9.8696
err_bias <- max(abs(bias_constant - pi^2))
r_ramp <- critical_rate(ramp, c(0, seq(0.5, 4000, length.out = 800)), function(r) -1,
                        lo = 0.9, hi = 1.2, tol = 1e-8)
err_ramp <- abs(r_ramp - 1)

# ---- test 2: the quasi-static lag ----
lag_ratio <- vapply(c(0.02, 0.01, 0.005), function(r) {
    sol <- integrate_fold(ramp, c(0, seq(1, 400, length.out = 400)), -1, r = r)
    (sol$u[nrow(sol)] + 1) / r
}, numeric(1))
err_lag <- max(abs(lag_ratio - 0.5))

# ---- tests 4 and 5: the asymptotically autonomous shift ----
t_shift <- c(-60, seq(-40, 60, length.out = 600))
magnitudes <- c(4, 8, 16)
r_closed_form <- 4 / (magnitudes - 2)
r_crit <- vapply(magnitudes, function(D)
    critical_rate(shift, t_shift, function(r) -1, lo = 4 / D * 0.4, hi = 4 / D * 8, D = D), numeric(1))
err_closed_form <- max(abs(r_crit - r_closed_form))
speed_ratio <- r_crit * magnitudes / 4                       # the weaker bound, at least 1
small_shift <- !escaped(integrate_fold(shift, t_shift, -1, r = 50, D = 1.5))   # D <= 2 never tips
# Below the critical rate the state must sit on the frozen equilibrium L(t) - 1 at the end of the
# run. The comparison is with L at the final time and not with its limit -D, because the shift
# itself is still a few parts in a million away from that limit at the end of the horizon.
tracked <- vapply(seq_along(magnitudes), function(k) {
    sol <- integrate_fold(shift, t_shift, -1, r = r_crit[k] * 0.9, D = magnitudes[k])
    last <- nrow(sol)
    abs(sol$x[last] - (shift(sol$time[last], r_crit[k] * 0.9, magnitudes[k]) - 1))
}, numeric(1))
end_point_gap <- abs(shift(max(t_shift), r_crit[1] * 0.9, magnitudes[1]) + magnitudes[1])
tipped <- vapply(seq_along(magnitudes), function(k) {
    sol <- integrate_fold(shift, t_shift, -1, r = r_crit[k] * 1.1, D = magnitudes[k])
    escaped(sol)
}, logical(1))

ok <- err_solver < 1e-8 && err_ramp < 1e-5 && err_bias < 0.2 && err_lag < 5e-3 &&
    err_closed_form < 1e-5 && all(speed_ratio >= 1) && small_shift &&
    max(tracked) < 1e-6 && all(tipped)

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    library(ggplot2); library(patchwork)
    D <- 8; rc <- r_crit[2]
    traj <- do.call(rbind, lapply(c(0.9, 1.1), function(f) {
        sol <- integrate_fold(shift, seq(-30, 30, by = 0.02), -1, r = rc * f, D = D)
        data.frame(time = sol$time, x = pmax(sol$x, -12),
                   case = sprintf("rate %.2f of critical", f))
    }))
    frozen <- data.frame(time = seq(-30, 30, by = 0.05))
    frozen$x <- shift(frozen$time, rc, D) - 1
    p1 <- ggplot(traj, aes(time, x, colour = case)) +
        geom_line(data = frozen, aes(time, x), inherit.aes = FALSE,
                  colour = "grey45", linetype = "dashed", linewidth = 0.5) +
        geom_line(linewidth = 0.7) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07"), name = NULL) +
        coord_cartesian(ylim = c(-12, 4)) +
        labs(title = "Tracking and tipping at the same fold",
             subtitle = kb_unicode(sprintf("Shift of magnitude $D = %d$; critical rate $r_c = %.4f$", D, rc)),
             x = kb_tex("Time $t$"), y = kb_tex("State $x$")) + kb_theme() + theme(legend.position = "top")
    curve <- data.frame(D = seq(2.4, 20, by = 0.1))
    curve$exact <- 4 / (curve$D - 2)
    curve$bound <- 4 / curve$D
    pts <- data.frame(D = magnitudes, r = r_crit)
    p2 <- ggplot(curve, aes(D, exact)) +
        geom_line(colour = "#4A6FA5", linewidth = 0.7) +
        geom_line(aes(y = bound), colour = "grey55", linetype = "dashed", linewidth = 0.5) +
        geom_point(data = pts, aes(D, r), colour = "#FB9E07", size = 2.4) +
        coord_cartesian(ylim = c(0, 4)) +
        labs(title = "Critical rate against the size of the shift",
             subtitle = kb_unicode("Blue: $4/(D - 2)$, exact. Orange: measured. Dashed: $4/D$"),
             x = kb_tex("Magnitude $D$ of the shift"), y = kb_tex("Critical rate $r_c$"),
             caption = kb_caption(paste(
                 "Left: the fold $dx/dt = (x - L(t))^2 - 1$ under $L(t) = -(D/2)(1 + \\tanh(rt/2))$, at nine tenths",
                 "of the critical rate and at eleven tenths; dashed, the moving equilibrium $L(t) - 1$. Right: the",
                 "critical rate of that shift, which the Poschl-Teller solution of the moving-frame Riccati",
                 "equation gives as $4/(D - 2)$, so a shift of magnitude 2 or less never tips.",
                 "Source: checks/rate-induced-tipping-critical-rate.R"))) + kb_theme()
    kb_save(p1 + p2, "rate-induced-tipping-critical-rate", width = 9.2, height = 4.2)
}

emit("rate-induced-tipping-critical-rate", if (ok) "pass" else "fail",
     sprintf("The linear ramp tips at rate %.6f against the closed form 1, and the tanh shift of magnitude D tips at 4/(D - 2), reproduced at D = 4, 8 and 16 with an error of %.1e",
             r_ramp, err_closed_form),
     list(critical_rate_ramp = r_ramp, error_against_closed_form = err_ramp,
          horizon_bias_constant = bias_constant[3], horizon_bias_error = err_bias,
          solver_error = err_solver, quasi_static_lag_ratio = lag_ratio[3], lag_ratio_error = err_lag,
          critical_rate_D4 = r_crit[1], critical_rate_D8 = r_crit[2], critical_rate_D16 = r_crit[3],
          error_against_closed_form_shift = err_closed_form, speed_ratio_D16 = speed_ratio[3],
          no_tipping_below_magnitude_two = small_shift,
          tracking_error_below_critical = max(tracked), end_point_gap_at_horizon = end_point_gap))
