# ============================ #
# Time and measure in the modified Bowen example ####
# ============================ #
#
# Claim, with the locators of kleptsyn2006 (Thm. 1 and Thm. 2) and ilyashenko2005 (Thm. 3):
#
# A planar field with a saddle node A and a hyperbolic saddle B joined in a heteroclinic
# cycle has statistical attractor {A, B} and minimal attractor {A}, so a neighbourhood of B
# keeps a share of the time of almost every orbit and loses its share of the transported
# Lebesgue measure.
#
# What this script tests is the mechanism in the normal-form model of that field, obtained
# from Lemma 1 and Lemma 2 of kleptsyn2006 by setting the remainder functions, which tend to
# 1 at the singular points, equal to 1, and by taking the transition maps between the two
# neighbourhoods to be the identity. With the saddle node in the normal form of Lemma 2 with
# parameter a and non-zero eigenvalue -b, and the saddle of Lemma 1 with eigenvalues lam and
# -mu, the section coordinate x near A and its reciprocal u = 1/x obey
#
#   passage time near A   t_A = G(u)/b,        G(u) = u - a log u,
#   Poincare map at A     xi  = x^{-a} exp(-1/x),
#   passage time near B   t_B = log(1/xi)/lam = G(u)/lam,
#   Poincare map at B     x'  = xi^{mu/lam},   so  log u' = (mu/lam) G(u).
#
# The script does not verify the theorems for a general field: the remainders and the
# transition times are dropped, and both are bounded by a constant per loop while G(u) grows
# beyond every bound, which is the sense in which the model is the asymptotic one.
#
# 1. The two Poincare maps and the two passage times are consistent: -log xi computed from
#    the formula for the Poincare map at A agrees with G(u), and the next section coordinate
#    agrees with the recursion for log u'. This tests the transcription of the two lemmas.
# 2. Time: at the end of every complete loop the fraction of time spent near B is exactly
#    (1/lam)/(1/b + 1/lam), because both passage times of a loop are proportional to the same
#    G; at the end of every passage near A that fraction collapses, since the current loop
#    is longer than all the previous ones together. So the fraction has no limit, its upper
#    limit is positive, and B belongs to the statistical attractor.
# 3. Measure: the mass that the transported Lebesgue measure gives the neighbourhood of B is
#    computed exactly, by solving for the two initial conditions whose orbits enter and leave
#    that neighbourhood at the given time, and it falls as the time grows.
# 4. The two readings are tied by Fubini: the time average over [0, T] of the transported
#    mass equals the average over initial conditions of the individual time fractions. The
#    script computes both sides by independent quadratures.

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

a <- 0.3      # parameter of the orbital normal form of the saddle node
b <- 1.0      # magnitude of the non-zero eigenvalue at the saddle node
lam <- 2.0    # unstable eigenvalue at the saddle
mu <- 3.0     # magnitude of the stable eigenvalue at the saddle, mu > lam so the loop contracts

G <- function(u) u - a * log(u)
share_B <- (1 / lam) / (1 / b + 1 / lam)   # 1/3 with these parameters
share_A <- (1 / b) / (1 / b + 1 / lam)

# ============================ #
# Test 1: the two lemmas, transcribed ####
# ============================ #
#
# Both identities are exact in the model, so the tolerance is the rounding of the elementary
# functions involved, taken here as 1e-12 relative.

xs <- c(0.5, 1/3, 0.2, 0.1, 0.05)
err_poincare <- 0; err_recursion <- 0
for (x in xs) {
    u <- 1 / x
    xi <- x^(-a) * exp(-1 / x)                       # Poincare map at the saddle node
    err_poincare <- max(err_poincare, abs(-log(xi) - G(u)) / G(u))
    x_next <- xi^(mu / lam)                          # Poincare map at the saddle
    err_recursion <- max(err_recursion, abs(log(1 / x_next) - (mu / lam) * G(u)) / ((mu / lam) * G(u)))
}

# ============================ #
# The schedule of one orbit ####
# ============================ #
#
# Returns, for the orbit through x on the section at A, the entry and exit times of the
# neighbourhood of B on each loop, and the cumulative time at the end of each loop. The
# recursion is stopped when log u leaves the double range, which happens after four loops.

schedule <- function(x, tmax, maxloop = 60L) {
    u <- 1 / x; tot <- 0
    enter <- numeric(0); leave <- numeric(0); Gs <- numeric(0)
    for (n in seq_len(maxloop)) {
        g <- G(u)
        enter <- c(enter, tot + g / b)
        leave <- c(leave, tot + g / b + g / lam)
        tot <- tot + g / b + g / lam
        Gs <- c(Gs, g)
        if (tot > tmax || !is.finite(tot)) break
        ln_next <- (mu / lam) * g
        if (ln_next > 700) { tot <- Inf; break }
        u <- exp(ln_next)
    }
    list(enter = enter, leave = leave, G = Gs)
}

# time spent in the neighbourhood of B before time t, by the orbit through x
time_in_B <- function(x, t) {
    s <- schedule(x, t)
    sum(pmax(0, pmin(s$leave, t) - pmin(s$enter, t)))
}

# ============================ #
# Test 2: the fraction of time, loop by loop ####
# ============================ #
#
# The predicted values are exact: share_B at the end of each complete loop, and at the end of
# the passage near A of loop n the value (S/lam) / (S (1/b + 1/lam) + G_n/b) with S the sum of
# the previous G. Both are compared with the schedule at 1e-12 relative.

x0 <- 1 / 3
s0 <- schedule(x0, Inf, maxloop = 4L)
Gs <- s0$G
n_loops <- length(Gs)
frac_loop_end <- numeric(n_loops); frac_A_end <- numeric(n_loops)
pred_A_end <- numeric(n_loops)
for (n in seq_len(n_loops)) {
    frac_loop_end[n] <- time_in_B(x0, s0$leave[n]) / s0$leave[n]
    frac_A_end[n] <- time_in_B(x0, s0$enter[n]) / s0$enter[n]
    S <- if (n == 1) 0 else sum(Gs[1:(n - 1)])
    pred_A_end[n] <- (S / lam) / (S * (1 / b + 1 / lam) + Gs[n] / b)
}
err_loop_end <- max(abs(frac_loop_end - share_B))
err_A_end <- max(abs(frac_A_end - pred_A_end))
oscillation <- max(frac_loop_end) - min(frac_A_end)

# ============================ #
# Test 3: the transported mass near the saddle ####
# ============================ #
#
# The entry and the exit time of the neighbourhood of B on loop n are decreasing functions of
# the section coordinate x, so for each loop the set of initial conditions that are in that
# neighbourhood at time t is an interval, whose endpoints are found by bisection to a relative
# accuracy of 1e-12. The mass is their total length divided by the length of the window.

xlo <- 1 / 6; xhi <- 1 / 2
nth_time <- function(x, n, which) {
    s <- schedule(x, Inf, maxloop = n)
    if (length(s$enter) < n) return(Inf)
    if (which == "enter") s$enter[n] else s$leave[n]
}
# For a given loop the set of initial conditions in the neighbourhood of B at time t is
# { x : enter_n(x) <= t <= leave_n(x) }. Both times decrease in x, so the first condition
# holds on an interval ending at xhi and the second on an interval starting at xlo, and the
# set is their intersection. Each endpoint is either a boundary of the window or a root found
# by bisection.
endpoint <- function(n, which, t) {
    f <- function(x) {
        v <- nth_time(x, n, which)
        if (!is.finite(v)) return(1)          # a time beyond the double range exceeds t
        sign(v - t) * min(1, abs(log(v) - log(t)))
    }
    list(flo = f(xlo), fhi = f(xhi), f = f)
}
mass_B <- function(t, maxloop = 6L) {
    total <- 0
    for (n in seq_len(maxloop)) {
        e <- endpoint(n, "enter", t)
        l <- endpoint(n, "leave", t)
        # { enter_n <= t } is [lo, xhi]
        if (e$flo <= 0) lo <- xlo
        else if (e$fhi > 0) next                     # the orbit has not entered anywhere
        else lo <- uniroot(e$f, c(xlo, xhi), tol = .Machine$double.eps^0.75)$root
        # { leave_n >= t } is [xlo, hi]
        if (l$fhi >= 0) hi <- xhi
        else if (l$flo < 0) next                     # the orbit has left everywhere
        else hi <- uniroot(l$f, c(xlo, xhi), tol = .Machine$double.eps^0.75)$root
        if (hi > lo) total <- total + (hi - lo)
    }
    total / (xhi - xlo)
}
ts <- 10^c(2, 4, 6, 8, 12, 16, 20, 25, 30)
mass <- vapply(ts, mass_B, numeric(1))
# the supremum over a block of decades, which is the quantity that must fall
sup_early <- max(mass[ts <= 1e8]); sup_late <- max(mass[ts >= 1e20])

# ============================ #
# Test 4: the same quantity read in the two ways ####
# ============================ #
#
# Fubini: the time average over [0, T] of the transported mass equals the mean over initial
# conditions of the fraction of time each orbit spends near B. The left side is a Riemann sum
# over 4000 times and the right side a Riemann sum over 4000 initial conditions, so the two
# agree only to the resolution of the coarser feature; the tolerance is 2 percent of the value.

Tf <- 200
t_grid <- seq(Tf / 4000, Tf, length.out = 4000L)
lhs <- mean(vapply(t_grid, mass_B, numeric(1)))
x_grid <- seq(xlo, xhi, length.out = 4000L)
rhs <- mean(vapply(x_grid, function(x) time_in_B(x, Tf), numeric(1))) / Tf
fubini_rel <- abs(lhs - rhs) / rhs

ok <- err_poincare < 1e-12 && err_recursion < 1e-12 &&
    err_loop_end < 1e-12 && err_A_end < 1e-12 &&
    abs(frac_loop_end[n_loops] - share_B) < 1e-12 && min(frac_A_end) < 1e-30 &&
    oscillation > 0.3 &&
    sup_late < 0.2 * sup_early && all(mass >= 0) &&
    fubini_rel < 0.02

emit("modified-bowen-example", if (ok) "pass" else "fail",
     "In the normal-form model of the modified Bowen example the share of time near the hyperbolic saddle returns to one third on every loop while the share of the transported measure it holds falls towards zero",
     list(poincare_map_rel_error = err_poincare,
          recursion_rel_error = err_recursion,
          loops_resolved = n_loops,
          fraction_in_B_at_loop_end = frac_loop_end[n_loops],
          closed_form_share_B = share_B,
          fraction_error_at_loop_end = err_loop_end,
          fraction_error_at_passage_end = err_A_end,
          smallest_fraction_at_passage_end = min(frac_A_end),
          oscillation_of_the_fraction = oscillation,
          mass_in_B_at_t_1e2 = mass[1],
          mass_in_B_at_t_1e12 = mass[ts == 1e12],
          mass_in_B_at_t_1e30 = mass[length(mass)],
          sup_mass_up_to_1e8 = sup_early,
          sup_mass_beyond_1e20 = sup_late,
          fubini_relative_gap = fubini_rel))

# ============================ #
# Figure ####
# ============================ #
#
# The separation in one picture, on a logarithmic time axis because the loops grow doubly
# exponentially. Above: the share of time that one orbit has spent near the hyperbolic saddle,
# which returns to the closed form on every loop and collapses in between, so it has no limit.
# Below: the share of the transported Lebesgue measure that the same neighbourhood holds, which
# falls away. The first says that the saddle belongs to the statistical attractor, the second
# that it does not belong to the minimal one.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    # The grid runs to just past the end of the third loop, which is the last one that double
    # precision resolves, and the flat stretches between the spikes are the passages near the
    # saddle node, each longer than everything that came before it.
    lg <- seq(0, 35.4, by = 0.06)
    tt <- 10^lg
    orbit <- data.frame(lg = lg, share = vapply(tt, function(t) time_in_B(x0, t) / t, numeric(1)))
    massd <- data.frame(lg = lg, share = vapply(tt, mass_B, numeric(1)))
    loop_ends <- log10(s0$leave[is.finite(s0$leave)])
    p1 <- ggplot(orbit, aes(lg, share)) +
        geom_vline(xintercept = loop_ends, colour = "grey80", linewidth = 0.3) +
        geom_hline(yintercept = share_B, colour = "#B8390E", linewidth = 0.4, linetype = "22") +
        geom_line(linewidth = 0.5, colour = "#21918c") +
        annotate("text", x = 0.3, y = share_B + 0.045, hjust = 0, size = 2.5, colour = "#B8390E",
                 label = sprintf("Closed form at the end of every loop: %.4f", share_B)) +
        annotate("text", x = 12, y = 0.10, hjust = 0, size = 2.5, colour = "grey35",
                 label = "The flat stretch is one passage near the saddle node,\nlonger than every loop before it together") +
        annotate("text", x = loop_ends[3] - 0.3, y = 0.40, hjust = 1, size = 2.4, colour = "grey45",
                 label = "Grey lines: the end of each loop") +
        scale_y_continuous(limits = c(0, 0.44)) +
        labs(title = "One orbit: the share of time spent near the hyperbolic saddle",
             subtitle = "It returns to the closed form on every loop and collapses in between, so it has no limit",
             x = NULL, y = "Share of the elapsed time") +
        kb_theme()
    floor_mass <- 1e-4
    p2 <- ggplot(massd, aes(lg, pmax(share, floor_mass))) +
        geom_vline(xintercept = loop_ends, colour = "grey80", linewidth = 0.3) +
        geom_line(linewidth = 0.5, colour = "#440154") +
        scale_y_log10(limits = c(floor_mass, 1),
                      breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
                      labels = kb_ticks(c("below $10^{-4}$", "0.001", "0.01", "0.1", "1"))) +
        labs(title = "The ensemble: the share of the transported measure held by the same neighbourhood",
             subtitle = "Each loop returns a smaller share than the loop before it, on a logarithmic scale",
             x = kb_tex("$\\log_{10}$ of the elapsed time"), y = "Share of the measure") +
        kb_theme()
    fig <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, ncol = 1) +
            patchwork::plot_annotation(
                caption = kb_caption(sprintf(
                    "Normal-form model of the field with a saddle node and a hyperbolic saddle, with $a = %.1f$, $b = %.1f$, $\\lambda = %.1f$ and $\\mu = %.1f$. The orbit starts at $x = %.3f$ on the section; the measure is Lebesgue measure on the section interval [%.3f, %.3f], and the share it gives the neighbourhood of the saddle is computed by solving for the two initial conditions that enter and leave it at each time. Both panels come from the run recorded by checks/modified-bowen-example.R, in which the share of time equalled %.6f at every loop end and the share of the measure fell from %.3f at $t = 10^2$ to %.4f at $t = 10^{30}$.",
                    a, b, lam, mu, x0, xlo, xhi, share_B, mass[1], mass[length(mass)])),
                theme = kb_theme())
    } else p1
    invisible(kb_save(fig, "modified-bowen-example", width = 8.6, height = 5.6))
}
