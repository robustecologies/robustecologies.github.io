# ============================ #
# Non-uniform fibre contraction in the sine-forced circle map ####
# ============================ #
#
# The map is F(y, x) = (3 y, x + a sin(2 pi x) + y) mod 1 on the two-torus, with a = 1/8, the
# numerical experiment of Homburg (2012), Figure 1. The notes on invariant graphs claim that this
# class sits in the non-uniform regime of the theory of forced systems: the fibre Lyapunov exponent
# of a typical orbit is negative, and yet no uniform bound |(f^n_y)'(x)| <= c delta^n with
# delta < 1 can hold, because the map carries an invariant measure of positive normal exponent.
# The script tests the four claims that this rests on and prints one JSON line for kb.py verify.
#
# 1. Fixed points with exponents of both signs. At y = 0 the base is fixed, 3 * 0 = 0, and the
#    fibre map is f_0(x) = x + a sin(2 pi x), which fixes x = 0 and x = 1/2. The script finds both
#    from the definition of F, by a bisection on f_0(x) - x that never uses the closed form, and
#    evaluates the fibre derivative there by a central finite difference of the map and by symbolic
#    differentiation with D(). The closed forms are f'(0) = 1 + 2 pi a = 1 + pi/4 and
#    f'(1/2) = 1 - 2 pi a = 1 - pi/4 at a = 1/8, so the normal exponents are log(1 + pi/4) > 0 and
#    log(1 - pi/4) < 0. The positive one is the obstruction: the Dirac measure at (0, 0) is
#    F-invariant and its normal exponent is positive, so the hypothesis of Theorem 1.13 of Sturman
#    and Stark (2000), that every invariant measure on the set has strictly negative exponent,
#    fails on the torus and on every compact invariant set that contains the point.
# 2. The supremum of the n-step exponent. Since every factor of (f^n_y)'(x) = prod f'(x_k) is at
#    most 1 + pi/4, the supremum over the torus of (1/n) log |(f^n_y)'(x)| is at most log(1 + pi/4)
#    for every n, and the fixed point of test 1 attains it. The script checks the bound over a grid
#    of the torus for n = 1 to 50 and checks that the fixed point attains it to machine precision.
#    A uniform bound c delta^n with delta < 1 would force this supremum below log delta < 0, so the
#    equality refutes uniform contraction for every c and every delta.
# 3. The typical exponent is negative, by two routes that share no code. Route A averages
#    log |f'(x_k)| along orbits over random base histories, using the analytic derivative. Route B
#    uses the map alone: two orbits are started 1e-9 apart in the same fibre and the exponent is
#    read from (1/n) log(sep_n / sep_0), which is a difference quotient and not a derivative. By
#    the mean value theorem the two differ by at most |(log f')'| = |f''| / f' <= (2 pi)^2 a /
#    (1 - pi/4) = 23.0 times the separation, which starts at 1e-7, so the linearisation error is at
#    most 2.3e-6 per step. Each history stops before its separation falls below 1e-11, since two
#    orbits closer than that are the same double and the difference loses all its digits; the
#    rounding error of the final separation is then at most 2.2e-16 / 1e-11 = 2.2e-5 in relative
#    terms, and dividing by the number of steps leaves it below 1e-5 for the shortest runs.
#    The base orbit of y -> 3y mod 1 is computed exactly from a ternary digit string, because
#    iterating 3y in double precision multiplies the rounding error by three at every step.
# 4. Degenerate case. At a = 0 every fibre map is a rotation, so every factor is 1 and route A and
#    the supremum are exactly 0. Route B is not exact there: it divides a difference of two numbers
#    of size one by 1e-7, so its floor is the rounding of that difference, 2.2e-16 / 1e-7 = 2.2e-9
#    in relative terms, which the run reproduces at the ninth decimal.
#
# Tolerances. The closed forms of test 1 are compared at 1e-12, above the rounding of a bisection
# on a smooth function and far below any difference that would matter. The finite-difference
# derivative is compared at 1e-6, the truncation error of a central difference at step 1e-5. The
# supremum of test 2 is compared at 1e-12 at the fixed point and the grid is required not to exceed
# it by more than 1e-12. The two routes of test 3 are compared at 1e-4, above the linearisation and
# rounding bounds derived above and far below the value of the exponent itself, which is -0.2. Test
# 4 is compared at 1e-15 for the derivative route and at 1e-6 for the difference route, whose floor
# is the rounding described above.

set.seed(20260918L)

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

# ---- the model ----
amp <- 1 / 8
fibre <- function(x, y, a = amp) (x + a * sin(2 * pi * x) + y) %% 1
dfibre <- function(x, a = amp) 1 + 2 * pi * a * cos(2 * pi * x)
base_orbit <- function(digits) {                       # exact orbit of y -> 3y mod 1
    K <- length(digits); p <- 3^-(1:K)
    vapply(0:(K - 1), function(j) sum(digits[(j + 1):K] * p[1:(K - j)]), numeric(1))
}

# ---- test 1: the two fixed points of the fibre map over the fixed base point ----
fixed_points <- function(a = amp) {
    displacement <- function(x) fibre(x, 0, a) - x     # zero exactly at a fixed point of f_0
    wrapped <- function(x) ((displacement(x) + 0.5) %% 1) - 0.5
    roots <- c(uniroot(wrapped, c(-0.2, 0.2), tol = .Machine$double.eps)$root,
               uniroot(wrapped, c(0.3, 0.7), tol = .Machine$double.eps)$root)
    roots
}
derivative_routes <- function(x, a = amp) {
    h <- 1e-5
    fd <- ((fibre(x + h, 0, a) - fibre(x - h, 0, a) + 0.5) %% 1 - 0.5) / (2 * h)
    sym <- eval(D(quote(x + a * sin(2 * pi * x)), "x"), list(x = x, a = a))
    c(finite_difference = fd, symbolic = sym, closed_form = dfibre(x, a))
}

fp <- fixed_points()
fp_residual <- max(abs(((fibre(fp, 0) - fp + 0.5) %% 1) - 0.5))
d_rep <- derivative_routes(fp[1]); d_att <- derivative_routes(fp[2])
err_fd <- max(abs(d_rep["finite_difference"] - d_rep["closed_form"]),
              abs(d_att["finite_difference"] - d_att["closed_form"]))
err_sym <- max(abs(d_rep["symbolic"] - (1 + pi / 4)), abs(d_att["symbolic"] - (1 - pi / 4)))
lambda_repelling <- log(unname(d_rep["symbolic"]))
lambda_attracting <- log(unname(d_att["symbolic"]))

# ---- test 2: the supremum of the n-step exponent over the torus ----
sup_exact <- log(1 + pi / 4)
n_grid <- 50L
grid_x <- seq(0, 1, length.out = 401)[-401]
grid_y <- seq(0, 1, length.out = 401)[-401]
sup_grid <- vapply(seq_len(n_grid), function(n) {
    best <- -Inf
    for (y0 in grid_y) {
        x <- grid_x; y <- y0; total <- rep(0, length(x))
        for (k in seq_len(n)) { total <- total + log(abs(dfibre(x))); x <- fibre(x, y); y <- (3 * y) %% 1 }
        best <- max(best, max(total) / n)
    }
    best
}, numeric(1))
sup_excess <- max(sup_grid - sup_exact)                       # must not exceed zero
at_fixed_point <- vapply(seq_len(n_grid), function(n) mean(rep(log(dfibre(fp[1])), n)), numeric(1))
sup_attained_err <- max(abs(at_fixed_point - sup_exact))

# ---- test 3: the typical exponent, by two routes ----
typical_exponent <- function(n_step = 400L, n_hist = 200L, sep0 = 1e-7, floor_sep = 1e-11, a = amp) {
    route_a <- numeric(n_hist); route_b <- numeric(n_hist); steps <- integer(n_hist)
    for (r in seq_len(n_hist)) {
        ys <- base_orbit(sample(0:2, n_step + 40L, replace = TRUE))
        x <- runif(1); x2 <- (x + sep0) %% 1
        acc <- 0; sep <- sep0; s <- 0L
        while (s < n_step && sep > floor_sep) {          # stop before the pair collapses to one double
            s <- s + 1L
            acc <- acc + log(abs(dfibre(x, a)))
            x <- fibre(x, ys[s], a); x2 <- fibre(x2, ys[s], a)
            sep <- abs(((x2 - x + 0.5) %% 1) - 0.5)
        }
        steps[r] <- s
        route_a[r] <- acc / s
        route_b[r] <- log(sep / sep0) / s
    }
    list(a = route_a, b = route_b, steps = steps)
}
short <- typical_exponent(n_step = 60L, n_hist = 200L)        # the pair stays linearised and resolved
lambda_a <- mean(short$a); lambda_b <- mean(short$b)
route_gap <- max(abs(short$a - short$b))
median_steps <- median(short$steps)
long <- typical_exponent(n_step = 400L, n_hist = 200L, sep0 = 1e-7, floor_sep = 0)
lambda_long <- mean(long$a)
lambda_se <- sd(long$a) / sqrt(length(long$a))

# ---- test 4: the degenerate case ----
flat <- typical_exponent(n_step = 40L, n_hist = 20L, a = 0)
flat_err_a <- max(abs(flat$a))                                # exact: every factor is 1
flat_err_b <- max(abs(flat$b))                                # limited by the rounding of a difference
flat_sup <- max(abs(log(abs(dfibre(grid_x, 0)))))

ok <- fp_residual < 1e-12 && err_fd < 1e-6 && err_sym < 1e-12 &&
    lambda_repelling > 0 && lambda_attracting < 0 &&
    sup_excess < 1e-12 && sup_attained_err < 1e-12 &&
    lambda_long < -0.1 && route_gap < 1e-4 && flat_err_a < 1e-15 && flat_err_b < 1e-6 && flat_sup < 1e-15

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    library(ggplot2)
    panel <- data.frame(exponent = long$a)
    p1 <- ggplot(panel, aes(x = exponent)) +
        geom_histogram(bins = 30, fill = "#4A6FA5", colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = lambda_long, colour = "#1A1A1A", linewidth = 0.5) +
        geom_vline(xintercept = sup_exact, colour = "#FB9E07", linewidth = 0.7) +
        geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
        labs(title = "Fibre exponents of the sine-forced circle map",
             subtitle = sprintf("400 steps over 200 random base histories; mean %.3f, supremum over the torus %.3f",
                                lambda_long, sup_exact),
             x = kb_tex("Exponent $\\lambda$ of one history"), y = "Histories",
             caption = kb_caption(paste(
                 "Blue: the distribution of $(1/n)\\sum \\log|f'(x_k)|$ over random base histories of the map",
                 "$F(y, x) = (3y, x + \\sin(2\\pi x)/8 + y)$. Black: their mean, which is negative, so a typical",
                 "orbit contracts. Orange: $\\log(1 + \\pi/4)$, the supremum of the same quantity over the torus,",
                 "attained at the repelling fixed point $(0, 0)$ for every $n$. A uniform bound $c\\delta^n$ with",
                 "$\\delta < 1$ would put",
                 "the orange line below zero, so contraction here is non-uniform.",
                 "Source: checks/non-uniform-fibre-contraction.R"))) +
        kb_theme()
    kb_save(p1, "non-uniform-fibre-contraction")
}

emit("non-uniform-fibre-contraction", if (ok) "pass" else "fail",
     sprintf("Typical fibre exponent %.3f is negative while the supremum over the torus is log(1 + pi/4) = %.3f, attained at a repelling fixed point",
             lambda_long, sup_exact),
     list(lambda_typical = lambda_long, lambda_standard_error = lambda_se,
          lambda_repelling_fixed_point = lambda_repelling, lambda_attracting_fixed_point = lambda_attracting,
          sup_exact = sup_exact, sup_grid_excess = sup_excess, sup_attained_error = sup_attained_err,
          fixed_point_residual = fp_residual, finite_difference_error = err_fd, symbolic_error = err_sym,
          two_route_gap = route_gap, two_route_median_steps = median_steps,
          degenerate_error_derivative = max(flat_err_a, flat_sup), degenerate_error_difference = flat_err_b))
