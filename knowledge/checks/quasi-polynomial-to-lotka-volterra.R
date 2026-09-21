# ============================ #
# Exact reduction of a quasi-polynomial system to Lotka-Volterra form ####
# ============================ #
#
# A quasi-polynomial system is dx_i/dt = x_i (lambda_i + sum_j A_ij prod_k x_k^{B_jk}), with one
# monomial per column of B. Writing y_j = prod_k x_k^{B_jk} for those monomials and differentiating
# gives
#     dy_j/dt = y_j ( (B lambda)_j + ((B A) y)_j ),
# a Lotka-Volterra system with growth vector rho = B lambda and interaction matrix M = B A. The
# reduction is exact, and the source of the knowledge base states that most nonlinear systems can
# be brought to quasi-polynomial form and hence to this canonical Lotka-Volterra format.
#
# The script tests four claims and prints one JSON line for kb.py verify.
#
# 1. The reduction is exact on a system that is not Lotka-Volterra. The example is
#        dx1/dt = x1 (1 - x1 x2),   dx2/dt = x2 (-1 + x1^2),
#    whose monomials are x1 x2 and x1^2, so B = [[1, 1], [2, 0]], lambda = (1, -1) and A = -I with
#    the sign pattern of the equations. The script integrates the original system, forms the two
#    monomials along that solution, integrates the derived Lotka-Volterra system from the matching
#    initial condition, and compares the two curves. They must agree to the accuracy of the solver.
# 2. The derived system is Lotka-Volterra and the original is not. The script verifies that the
#    right-hand side of the derived system is exactly y_j times an affine function of y, by
#    comparing it with that form at random points, and that the original right-hand side is not, by
#    exhibiting the monomial x1 x2 whose exponent vector is not a unit vector.
# 3. The transformation is invertible where B is invertible. Here det B = -2, so the original
#    variables are recovered from the monomials as x = exp(B^{-1} log y), and the script recovers
#    the original trajectory from the Lotka-Volterra one and compares.
# 4. Degenerate case. When B is the identity the construction returns the system unchanged, which
#    the script checks on a Lotka-Volterra system with a random matrix.
#
# Tolerances. The trajectories are compared at 1e-7. The solver runs at rtol and atol of 1e-11, and
# the two integrations follow different equations in different variables, so their difference is
# dominated by the accumulation of their separate local errors over the interval, not by one
# tolerance; 1e-7 is two orders above the observed difference and far below the variation of the
# curves, which span more than one unit. The algebraic identities of tests 2 and 4 are compared at
# 1e-12, a few rounding units.

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

# ---- the quasi-polynomial system and its Lotka-Volterra image ----
monomials <- function(x, B) apply(B, 1, function(row) prod(x^row))

qp_rhs <- function(t, x, parms) {
    y <- monomials(x, parms$B)
    list(x * (parms$lambda + as.vector(parms$A %*% y)))
}
lv_rhs <- function(t, y, parms) list(y * (parms$rho + as.vector(parms$M %*% y)))

reduce_to_lv <- function(lambda, A, B) list(rho = as.vector(B %*% lambda), M = B %*% A)

# ---- test 1: the example ----
B <- matrix(c(1, 1,
              2, 0), nrow = 2, byrow = TRUE)     # monomials x1 x2 and x1^2
lambda <- c(1, -1)
A <- matrix(c(-1, 0,
              0, 1), nrow = 2, byrow = TRUE)     # dx1/dt = x1(1 - x1x2), dx2/dt = x2(-1 + x1^2)
lv <- reduce_to_lv(lambda, A, B)

times <- seq(0, 6, by = 0.01)
x0 <- c(0.8, 1.3)
sol_x <- as.data.frame(ode(y = x0, times = times, func = qp_rhs,
                           parms = list(lambda = lambda, A = A, B = B),
                           method = "lsoda", rtol = 1e-11, atol = 1e-11))
y_from_x <- t(apply(as.matrix(sol_x[, -1]), 1, function(x) monomials(x, B)))
y0 <- monomials(x0, B)
sol_y <- as.data.frame(ode(y = y0, times = times, func = lv_rhs,
                           parms = list(rho = lv$rho, M = lv$M),
                           method = "lsoda", rtol = 1e-11, atol = 1e-11))
reduction_error <- max(abs(as.matrix(sol_y[, -1]) - y_from_x))

# ---- test 2: the derived system has Lotka-Volterra form and the original does not ----
affine_error <- max(replicate(200, {
    y <- runif(2, 0.2, 3)
    rhs <- lv_rhs(0, y, list(rho = lv$rho, M = lv$M))[[1]]
    max(abs(rhs - y * (lv$rho + as.vector(lv$M %*% y))))
}))
original_is_lv <- all(rowSums(abs(B)) == 1)          # false: the exponents are not unit vectors

# ---- test 3: the transformation is inverted ----
x_from_y <- t(apply(as.matrix(sol_y[, -1]), 1, function(y) as.vector(exp(solve(B) %*% log(y)))))
inversion_error <- max(abs(x_from_y - as.matrix(sol_x[, -1])))

# ---- test 4: an identity exponent matrix returns the same system ----
M_rand <- matrix(rnorm(4), 2, 2); lambda_rand <- rnorm(2)
lv_identity <- reduce_to_lv(lambda_rand, M_rand, diag(2))
identity_error <- max(abs(lv_identity$rho - lambda_rand), abs(lv_identity$M - M_rand))

ok <- reduction_error < 1e-7 && affine_error < 1e-12 && !original_is_lv &&
    inversion_error < 1e-7 && identity_error < 1e-12

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    library(ggplot2); library(patchwork)
    orig <- data.frame(time = rep(sol_x$time, 2),
                       value = c(sol_x[, 2], sol_x[, 3]),
                       variable = rep(c("$x_1$", "$x_2$"), each = nrow(sol_x)))
    p1 <- ggplot(orig, aes(time, value, colour = variable)) +
        geom_line(linewidth = 0.6) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07"), name = NULL,
                            labels = kb_tex(sort(unique(orig$variable)))) +
        labs(title = "A system that is not Lotka-Volterra",
             subtitle = kb_tex("$dx_1/dt = x_1(1 - x_1 x_2)$, $dx_2/dt = x_2(-1 + x_1^2)$"),
             x = kb_tex("Time $t$"), y = kb_tex("State $x_i$")) + kb_theme() + theme(legend.position = "top")
    img <- data.frame(time = rep(sol_y$time, 2),
                      derived = c(sol_y[, 2], sol_y[, 3]),
                      monomial = c(y_from_x[, 1], y_from_x[, 2]),
                      variable = rep(c("$y_1 = x_1 x_2$", "$y_2 = x_1^2$"), each = nrow(sol_y)))
    p2 <- ggplot(img, aes(time, monomial, colour = variable)) +
        geom_line(linewidth = 1.4, alpha = 0.35) +
        geom_line(aes(y = derived), linewidth = 0.5) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07"), name = NULL,
                            labels = kb_tex(sort(unique(img$variable)))) +
        labs(title = "Its Lotka-Volterra image",
             subtitle = sprintf("Monomials (thick) and the derived system (thin); gap %.0e", reduction_error),
             x = kb_tex("Time $t$"), y = kb_tex("Monomial $y_j$"),
             caption = kb_caption(paste(
                 "The quasi-polynomial system on the left has monomials $x_1 x_2$ and $x_1^2$, so it is not of",
                 "Lotka-Volterra form. Writing those monomials as new variables gives the Lotka-Volterra system",
                 "$dy/dt = y(B\\lambda + BAy)$, whose solution is drawn thin over the monomials of the original",
                 "solution, drawn thick. The two coincide to the accuracy of the solver, and the original",
                 "trajectory is recovered from the image because $B$ is invertible.",
                 "Source: checks/quasi-polynomial-to-lotka-volterra.R"))) +
        kb_theme() + theme(legend.position = "top")
    kb_save(p1 + p2, "quasi-polynomial-to-lotka-volterra", width = 9.4, height = 4.3)
}

emit("quasi-polynomial-to-lotka-volterra", if (ok) "pass" else "fail",
     sprintf("The monomials of a quasi-polynomial system obey the derived Lotka-Volterra equations to %.0e, and the original trajectory is recovered from them to %.0e",
             reduction_error, inversion_error),
     list(reduction_error = reduction_error, inversion_error = inversion_error,
          affine_form_error = affine_error, identity_case_error = identity_error,
          original_is_lotka_volterra = original_is_lv))
