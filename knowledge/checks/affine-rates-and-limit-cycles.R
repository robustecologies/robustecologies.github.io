# ============================ #
# What affine per capita rates exclude in the plane ####
# ============================ #
#
# A planar Kolmogorov system is dx/dt = x S(x, y), dy/dt = y W(x, y) on the open positive quadrant.
# The source of kb/concepts/kolmogorov-system.md reports that Moiseev showed in 1939 that such a
# system with affine S and W admits no limit cycle, which separates the Lotka-Volterra model from
# the saturating models inside the same family. The script does not test that statement in its full
# generality, which is a non-existence claim over all parameters; it tests the two cases in which a
# computation can decide, and it tests that the saturating case behaves otherwise.
#
# 1. Affine rates with self-limitation admit no periodic orbit at all. With
#        S(x, y) = b1 + a11 x + a12 y,   W(x, y) = b2 + a21 x + a22 y,
#    the Dulac function B(x, y) = 1/(xy) gives
#        d/dx (B x S) + d/dy (B y W) = a11 / y + a22 / x,
#    which is strictly negative on the open quadrant when a11 < 0 and a22 < 0, so Dulac's criterion
#    excludes a periodic orbit there. The script computes the left side by numerical differentiation
#    of the vector field itself, which shares no algebra with the hand-derived right side, and
#    compares the two on a grid; it then integrates from twenty initial conditions and measures the
#    amplitude of the last quarter of each trajectory, which must fall to zero.
# 2. The conservative case is the boundary of that argument. With a11 = a22 = 0 the divergence above
#    vanishes identically, Dulac's criterion says nothing, and the classical Volterra system has a
#    continuum of closed orbits and no isolated one. The script verifies the first integral
#        H(x, y) = d x - c log x + b y - a log y
#    of dx/dt = x(a - b y), dy/dt = y(-c + d x) along the trajectories, and verifies that the
#    amplitude varies continuously with the initial condition, which is what distinguishes a
#    continuum of closed orbits from a limit cycle. The first integral is differentiated
#    symbolically with D(), so the test performs the differentiation that could be wrong.
# 3. A saturating rate in the same family produces a limit cycle. The Rosenzweig-MacArthur system
#        S(x, y) = r (1 - x / K) - a y / (1 + a h x),   W(x, y) = e a x / (1 + a h x) - m
#    has the interior equilibrium x* = m / (a (e - m h)), y* = (r / a)(1 - x*/K)(1 + a h x*), which
#    loses stability when x* < (K - 1/(a h)) / 2. With r = 1, a = 1, h = 0.5, e = 0.5 and m = 0.2
#    the equilibrium sits at x* = 0.5 and the threshold is K = 3. The script runs the system at
#    K = 6 and at K = 2.5, one on each side of that analytic threshold, and requires an oscillation
#    of fixed amplitude in the first case and convergence in the second. The amplitude at K = 6 is
#    measured from two initial conditions that are far apart, which is the independent route: a
#    limit cycle is an isolated periodic orbit, so the two must agree, whereas the closed orbits of
#    test 2 must not. The stability of the equilibrium is decided from the eigenvalues of a
#    finite-difference Jacobian and compared with the analytic threshold.
#
# Tolerances. The Dulac divergence is compared at 1e-6, which is the accuracy of a central
# difference of step 1e-5 on a field whose second derivatives are of order one. The first integral
# is required to vary by less than 1e-6 over the trajectory, against a solver run at rtol and atol
# of 1e-11, and its symbolic time derivative is compared with zero at 1e-10. Amplitudes are called
# zero below 1e-6 and nonzero above 1e-2, a gap of four orders that the observed values do not come
# near, and the two amplitudes on the limit cycle are compared at 1e-4, which is the accuracy with
# which a finite integration reaches the cycle.

set.seed(20260919L)
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

amplitude <- function(sol, frac = 0.25) {
    keep <- seq(floor((1 - frac) * nrow(sol)) + 1, nrow(sol))
    max(apply(as.matrix(sol[keep, -1, drop = FALSE]), 2, function(v) max(v) - min(v)))
}

# ---- 1: affine rates with self-limitation ----
lv_par <- list(b = c(1.0, 0.6), A = matrix(c(-0.8, -0.5, 0.4, -0.7), 2, 2, byrow = TRUE))
lv_field <- function(state, p) state * (p$b + as.vector(p$A %*% state))
lv_rhs <- function(t, state, p) list(lv_field(state, p))

test_dulac_affine <- function() {
    p <- lv_par
    h <- 1e-5
    grid <- expand.grid(x = seq(0.05, 3, length.out = 25), y = seq(0.05, 3, length.out = 25))
    # numerical divergence of B * field, with B = 1/(xy)
    div_num <- mapply(function(x, y) {
        g <- function(u) lv_field(u, p) / (u[1] * u[2])
        (g(c(x + h, y))[1] - g(c(x - h, y))[1]) / (2 * h) +
            (g(c(x, y + h))[2] - g(c(x, y - h))[2]) / (2 * h)
    }, grid$x, grid$y)
    div_closed <- p$A[1, 1] / grid$y + p$A[2, 2] / grid$x
    err_div <- max(abs(div_num - div_closed))
    starts <- matrix(runif(40, 0.05, 3), ncol = 2)
    amps <- apply(starts, 1, function(s) {
        sol <- as.data.frame(ode(s, seq(0, 400, by = 0.2), lv_rhs, p,
                                 method = "lsoda", rtol = 1e-11, atol = 1e-11))
        amplitude(sol)
    })
    list(err_div = err_div, max_div = max(div_closed), max_amplitude = max(amps),
         ok = err_div < 1e-6 && max(div_closed) < 0 && max(amps) < 1e-6)
}

# ---- 2: the conservative case, a continuum of closed orbits ----
volterra <- list(a = 1.1, b = 0.8, c = 0.9, d = 0.7)
vo_rhs <- function(t, state, p) list(c(state[1] * (p$a - p$b * state[2]),
                                       state[2] * (-p$c + p$d * state[1])))

test_conservative_centre <- function() {
    p <- volterra
    H <- function(x, y) p$d * x - p$c * log(x) + p$b * y - p$a * log(y)
    # symbolic derivative of H along the field, which must vanish identically
    expr <- quote(d * x - c * log(x) + b * y - a * log(y))
    dHdx <- D(expr, "x"); dHdy <- D(expr, "y")
    pts <- data.frame(x = runif(200, 0.2, 4), y = runif(200, 0.2, 4))
    env <- list2env(p)
    dH_dt <- mapply(function(x, y) {
        assign("x", x, envir = env); assign("y", y, envir = env)
        gx <- eval(dHdx, env); gy <- eval(dHdy, env)
        gx * x * (p$a - p$b * y) + gy * y * (-p$c + p$d * x)
    }, pts$x, pts$y)
    err_symbolic <- max(abs(dH_dt))
    starts <- cbind(c(1.4, 2.0, 2.6), rep(p$a / p$b, 3))
    orbits <- lapply(seq_len(nrow(starts)), function(i)
        as.data.frame(ode(starts[i, ], seq(0, 60, by = 0.02), vo_rhs, p,
                          method = "lsoda", rtol = 1e-11, atol = 1e-11)))
    err_integral <- max(vapply(orbits, function(o) diff(range(H(o[, 2], o[, 3]))), numeric(1)))
    amps <- vapply(orbits, amplitude, numeric(1))
    list(err_symbolic = err_symbolic, err_integral = err_integral,
         amplitude_spread = diff(range(amps)),
         ok = err_symbolic < 1e-10 && err_integral < 1e-6 && diff(range(amps)) > 1e-2)
}

# ---- 3: a saturating rate produces a limit cycle ----
rm_par <- function(K) list(r = 1, K = K, a = 1, h = 0.5, e = 0.5, m = 0.2)
rm_field <- function(state, p) {
    x <- state[1]; y <- state[2]
    fr <- p$a * x / (1 + p$a * p$h * x)
    c(x * (p$r * (1 - x / p$K) - p$a * y / (1 + p$a * p$h * x)),
      y * (p$e * fr - p$m))
}
rm_rhs <- function(t, state, p) list(rm_field(state, p))

test_saturating_cycle <- function() {
    out <- list()
    for (K in c(6, 2.5)) {
        p <- rm_par(K)
        x_star <- p$m / (p$a * (p$e - p$m * p$h))
        y_star <- (p$r / p$a) * (1 - x_star / p$K) * (1 + p$a * p$h * x_star)
        threshold <- (p$K - 1 / (p$a * p$h)) / 2
        h <- 1e-6
        J <- vapply(1:2, function(j) {
            e <- rep(0, 2); e[j] <- h
            (rm_field(c(x_star, y_star) + e, p) - rm_field(c(x_star, y_star) - e, p)) / (2 * h)
        }, numeric(2))
        re_max <- max(Re(eigen(J)$values))
        sols <- lapply(list(c(0.9 * x_star, 0.9 * y_star), c(0.2, 0.2)), function(s)
            as.data.frame(ode(s, seq(0, 3000, by = 0.05), rm_rhs, p,
                              method = "lsoda", rtol = 1e-11, atol = 1e-11)))
        amps <- vapply(sols, amplitude, numeric(1))
        out[[as.character(K)]] <- list(x_star = x_star, threshold = threshold, re_max = re_max,
                                       amps = amps)
    }
    unstable <- out[["6"]]
    stable <- out[["2.5"]]
    list(x_star = unstable$x_star, threshold_unstable = unstable$threshold,
         threshold_stable = stable$threshold,
         re_max_unstable = unstable$re_max, re_max_stable = stable$re_max,
         amplitude_cycle = mean(unstable$amps), amplitude_gap = abs(diff(unstable$amps)),
         amplitude_stable = max(stable$amps),
         ok = unstable$x_star < unstable$threshold && unstable$re_max > 0 &&
             stable$x_star > stable$threshold && stable$re_max < 0 &&
             min(unstable$amps) > 1e-2 && abs(diff(unstable$amps)) < 1e-4 &&
             max(stable$amps) < 1e-6)
}

t1 <- test_dulac_affine()
t2 <- test_conservative_centre()
t3 <- test_saturating_cycle()
ok <- t1$ok && t2$ok && t3$ok

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    lv_orbits <- do.call(rbind, lapply(1:6, function(i) {
        s <- c(runif(1, 0.1, 2.8), runif(1, 0.1, 2.8))
        o <- as.data.frame(ode(s, seq(0, 120, by = 0.05), lv_rhs, lv_par,
                               method = "lsoda", rtol = 1e-11, atol = 1e-11))
        data.frame(x = o[, 2], y = o[, 3], orbit = as.character(i))
    }))
    lv_eq <- solve(lv_par$A, -lv_par$b)
    p1 <- ggplot(lv_orbits, aes(x, y, group = orbit)) +
        geom_path(colour = "#4A6FA5", linewidth = 0.4, alpha = 0.9) +
        geom_point(data = data.frame(x = lv_eq[1], y = lv_eq[2], orbit = "eq"),
                   colour = "grey20", fill = "white", shape = 21, size = 2.2, stroke = 0.8) +
        labs(title = "Affine rates with self-limitation",
             subtitle = "Every orbit reaches the equilibrium",
             x = kb_tex("Prey $x$"), y = kb_tex("Predator $y$")) + kb_theme()
    vo_orbits <- do.call(rbind, lapply(c(1.4, 2.0, 2.6), function(x0) {
        o <- as.data.frame(ode(c(x0, volterra$a / volterra$b), seq(0, 60, by = 0.02), vo_rhs,
                               volterra, method = "lsoda", rtol = 1e-11, atol = 1e-11))
        data.frame(x = o[, 2], y = o[, 3], orbit = sprintf("x(0) = %.1f", x0))
    }))
    p2 <- ggplot(vo_orbits, aes(x, y, colour = orbit)) +
        geom_path(linewidth = 0.4) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07", "#7A5195"), name = NULL) +
        labs(title = "Affine rates without self-limitation",
             subtitle = "Closed orbits, none of them isolated",
             x = kb_tex("Prey $x$"), y = kb_tex("Predator $y$")) + kb_theme()
    rm_orbits <- do.call(rbind, lapply(list(c(0.45, 0.9), c(0.2, 0.2)), function(s) {
        o <- as.data.frame(ode(s, seq(0, 400, by = 0.02), rm_rhs, rm_par(6),
                               method = "lsoda", rtol = 1e-11, atol = 1e-11))
        data.frame(x = o[, 2], y = o[, 3], orbit = sprintf("(%.2f, %.2f)", s[1], s[2]))
    }))
    p3 <- ggplot(rm_orbits, aes(x, y, colour = orbit)) +
        geom_path(linewidth = 0.4, alpha = 0.9) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07"), name = "Initial state") +
        labs(title = "A saturating rate in the same family",
             subtitle = "Two orbits reach one limit cycle",
             x = kb_tex("Prey $x$"), y = kb_tex("Predator $y$")) + kb_theme()
    fig <- (p1 | p2 | p3) +
        plot_annotation(caption = kb_caption(sprintf(
            "All three systems have the form $dx_i/dt = x_i f_i(x)$. Left: $f$ affine with $a_{11} = %.1f$ and $a_{22} = %.1f$, where the Dulac function $1/(xy)$ gives a divergence of at most %.3f on the quadrant drawn, and the largest amplitude left after $t = 400$ is %.1e. Middle: the conservative case $a_{11} = a_{22} = 0$, where the first integral is constant to %.1e along each orbit and the amplitude changes with the initial condition, so no orbit is isolated. Right: a saturating $f$ with $K = 6$, past the analytic instability threshold $K = 3$, where two orbits reach a cycle of amplitude %.3f and the two amplitudes differ by %.1e. Drawn by checks/affine-rates-and-limit-cycles.R.",
            lv_par$A[1, 1], lv_par$A[2, 2], t1$max_div, t1$max_amplitude, t2$err_integral, t3$amplitude_cycle, t3$amplitude_gap)),
            theme = kb_theme())
    kb_save(fig, "affine-rates-and-limit-cycles", width = 9.2, height = 4.0)
}

emit("affine-rates-and-limit-cycles", if (ok) "pass" else "fail",
     "Affine per capita rates in the plane give a negative Dulac divergence or a continuum of closed orbits, while a saturating rate in the same family gives a limit cycle",
     list(err_dulac_divergence = t1$err_div, max_dulac_divergence = t1$max_div,
          amplitude_affine = t1$max_amplitude,
          err_first_integral_symbolic = t2$err_symbolic, drift_first_integral = t2$err_integral,
          amplitude_spread_conservative = t2$amplitude_spread,
          equilibrium_prey = t3$x_star, threshold_unstable = t3$threshold_unstable,
          leading_eigenvalue_unstable = t3$re_max_unstable,
          leading_eigenvalue_stable = t3$re_max_stable,
          amplitude_limit_cycle = t3$amplitude_cycle, amplitude_difference = t3$amplitude_gap,
          amplitude_below_threshold = t3$amplitude_stable))
