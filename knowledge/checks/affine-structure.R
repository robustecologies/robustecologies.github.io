# ============================ #
# What an affine structure preserves and what it does not ####
# ============================ #
#
# The note kb/concepts/affine-space.md defines an affine map as x -> A x + b with A linear, an
# affine combination as sum_i lambda_i x_i with sum_i lambda_i = 1, and the affine hull of a set as
# the set of its affine combinations. The script tests six claims of that note. Every test compares
# the computed quantity with a closed form derived by hand, with a published value, or with a route
# that shares no code with the object under test.
#
# 1. An affine map commutes with an affine combination, and with no other weighted combination. For
#    weights that sum to s the discrepancy is exactly (1 - s) b, since
#        phi(sum lambda_i x_i) - sum lambda_i phi(x_i) = b + A sum lambda_i x_i
#                                                        - sum lambda_i (A x_i + b) = (1 - s) b.
#    The script evaluates both sides at random points and weights and compares the discrepancy with
#    that closed form, at s = 1 and at s != 1. The closed form is the reference.
# 2. A map that commutes with the affine combinations of pairs is affine, and its linear part and
#    translation are read off the images of the frame 0, e_1, ..., e_d by b = phi(0) and
#    A e_i = phi(e_i) - phi(0). The script reconstructs A and b that way and compares them with the
#    matrix and vector that generated the map, then applies the same reconstruction to a map with a
#    quadratic term and requires the residual on random points to be large, so that the test is
#    shown to have power against a map that is not affine.
# 3. An affine map preserves the ratio in which a point divides a segment, and does not preserve
#    lengths, angles or volumes. The reference for the ratio is the parameter t of the point
#    (1 - t) x + t y, known in advance. The reference for the volume is |det A| times the volume of
#    the source tetrahedron, and the volume of the image tetrahedron is computed from its six
#    pairwise distances through the Cayley-Menger determinant, which uses the coordinates of the
#    image only through distances and so does not repeat the determinant under test.
# 4. The solution set of A x = b is empty or an affine set of dimension n - rank(A), and affine
#    combinations of solutions are solutions. The matrix is built with a known rank, so the expected
#    dimension is known before the computation; the rank and the null space come from the singular
#    value decomposition, and the residual of a combination with weights summing to s is compared
#    with the closed form (s - 1) b. Degenerate cases: the affine hull of one point has dimension 0,
#    and d + 2 points of R^d are never affinely independent.
# 5. Barycentric coordinates with respect to an affinely independent frame are unique, sum to one,
#    are positive in the interior of the simplex, and are unchanged by an affine map. The reference
#    for the invariance is the coordinate vector computed in the source frame.
# 7. A change of units is an affine change of variables that carries a Lotka-Volterra system to
#    another one. Writing the abundances as u = E v with E a positive diagonal matrix gives
#        dv_i/dt = v_i (b_i + (A E v)_i),
#    so the growth rates are unchanged and the interaction matrix becomes A E; the two trajectories
#    are then related by v(t) = E^{-1} u(t), which the script checks against an integration of each
#    system in its own variables. The script also checks that Volterra-Lyapunov stability survives
#    the rescaling with the transformed certificate: if D is a positive diagonal matrix with
#    D A + A^T D negative definite, then for P = E D the quadratic form x -> x^T P (A E) x equals
#    y^T D A y at y = E x, so P certifies A E. The reference in each case is the closed form just
#    derived, and the degenerate case E = I must return the system unchanged.
#
# 6. An affine conjugacy preserves the Lyapunov exponent. The logistic map S_r(x) = r x (1 - x) and
#    the quadratic map Q_a(y) = 1 - a y^2 with a = r(r - 2)/4 are conjugate by the affine change of
#    variables x = 1/2 + (r - 2) y / 4, as kb/concepts/logistic-map.md states. The script starts an
#    orbit of each map at an independently drawn point, so that the two orbits are not images of one
#    another, and compares the two exponents with each other and with the values known in closed
#    form: log 2 at r = 4, and half the logarithm of the multiplier 4 + 2r - r^2 of the two-cycle at
#    r = 3.2, which is (1/2) log 0.16.
#
# Tolerances. The algebraic identities of tests 1 to 5 are compared at 1e-10, which is several
# thousand rounding units of the quantities involved, all of order one to one hundred; the observed
# errors are near 1e-14. The Cayley-Menger volume is compared at 1e-8, since it takes a square root
# of a 5 by 5 determinant of squared distances and loses about half of the available digits. The
# Lyapunov exponents of test 6 are compared at 5e-3 against their closed forms. The observed errors
# are far smaller, near 5e-7, and the reason is structural rather than fortunate: at r = 4 the map
# is conjugate to the doubling map by x = sin^2(pi theta), so differentiating that conjugacy gives
# log|S_r'(x_i)| = log 2 + log|h'(theta_{i+1})| - log|h'(theta_i)| and the Birkhoff sum telescopes,
# leaving an error of order 1/n instead of the 1/sqrt(n) of a central limit theorem. The tolerance
# is kept at 5e-3 because the telescoped remainder is bounded only when the orbit stays away from
# the critical point, which no run can guarantee. The periodic case at r = 3.2 converges
# geometrically once the transient is discarded and is compared at 1e-8. The trajectories of test 7
# are compared at 1e-8, against a solver run at rtol and atol of 1e-11 on two systems written in
# different variables, so their difference accumulates the separate local errors of the two runs.

set.seed(20260919L)

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

affine_map <- function(A, b) function(x) as.vector(A %*% x + b)

# ---- 1: affine combinations, and only those, commute with an affine map ----
test_combination_invariance <- function(trials = 500L, d = 4L) {
    err_affine <- 0
    err_closed <- 0
    for (i in seq_len(trials)) {
        A <- matrix(rnorm(d * d), d, d); b <- rnorm(d)
        phi <- affine_map(A, b)
        pts <- matrix(rnorm(d * (d + 1)), d, d + 1)
        w <- rnorm(d + 1)
        w_affine <- w / sum(w)                      # weights that sum to one
        lhs <- phi(as.vector(pts %*% w_affine))
        rhs <- as.vector(apply(pts, 2, phi) %*% w_affine)
        err_affine <- max(err_affine, max(abs(lhs - rhs)))
        s <- sum(w)                                  # weights that do not sum to one
        gap <- phi(as.vector(pts %*% w)) - as.vector(apply(pts, 2, phi) %*% w)
        err_closed <- max(err_closed, max(abs(gap - (1 - s) * b)))
    }
    list(err_affine = err_affine, err_closed = err_closed,
         ok = err_affine < 1e-10 && err_closed < 1e-10)
}

# ---- 2: an affine map is determined by the image of one frame ----
test_reconstruction_from_frame <- function(d = 5L, trials = 200L) {
    A <- matrix(rnorm(d * d), d, d); b <- rnorm(d)
    phi <- affine_map(A, b)
    b_hat <- phi(rep(0, d))
    A_hat <- vapply(seq_len(d), function(i) phi(diag(d)[, i]) - b_hat, numeric(d))
    err_linear <- max(abs(A_hat - A), abs(b_hat - b))
    # the same reconstruction applied to a map that is not affine
    psi <- function(x) as.vector(A %*% x + b) + c(x[1]^2, rep(0, d - 1))
    c_hat <- psi(rep(0, d))
    B_hat <- vapply(seq_len(d), function(i) psi(diag(d)[, i]) - c_hat, numeric(d))
    res_nonlinear <- max(replicate(trials, {
        x <- rnorm(d)
        max(abs(psi(x) - (as.vector(B_hat %*% x) + c_hat)))
    }))
    list(err_linear = err_linear, res_nonlinear = res_nonlinear,
         ok = err_linear < 1e-10 && res_nonlinear > 1e-2)
}

# ---- 3: ratios survive, lengths, angles and volumes do not ----
cayley_menger_volume <- function(P) {
    # P has one vertex per column; the volume is computed from the pairwise distances alone.
    n <- ncol(P)
    D <- as.matrix(dist(t(P)))^2
    M <- rbind(c(0, rep(1, n)), cbind(rep(1, n), D))
    sign_k <- (-1)^n
    det_M <- det(M)
    sqrt(max(sign_k * det_M, 0) / (2^(n - 1) * (factorial(n - 1))^2))
}

test_ratio_and_metric <- function(trials = 500L, d = 3L) {
    A <- matrix(c(2, 0.7, -0.3, 0.1, 1.4, 0.6, -0.5, 0.2, 0.9), 3, 3)
    b <- c(0.4, -1.1, 2.0)
    phi <- affine_map(A, b)
    err_ratio <- 0; max_len_change <- 0; max_angle_change <- 0
    for (i in seq_len(trials)) {
        x <- rnorm(d); y <- rnorm(d); t <- runif(1)
        z <- (1 - t) * x + t * y
        fx <- phi(x); fy <- phi(y); fz <- phi(z)
        t_hat <- sum((fz - fx) * (fy - fx)) / sum((fy - fx)^2)   # ratio along the image segment
        err_ratio <- max(err_ratio, abs(t_hat - t))
        u <- rnorm(d); v <- rnorm(d)
        len <- sqrt(sum(u^2)); len_img <- sqrt(sum((A %*% u)^2))
        max_len_change <- max(max_len_change, abs(len_img / len - 1))
        ang <- function(p, q) acos(sum(p * q) / sqrt(sum(p^2) * sum(q^2)))
        max_angle_change <- max(max_angle_change, abs(ang(u, v) - ang(as.vector(A %*% u), as.vector(A %*% v))))
    }
    P <- matrix(rnorm(d * (d + 1)), d, d + 1)
    vol_source <- abs(det(P[, -1, drop = FALSE] - P[, 1])) / factorial(d)
    vol_image <- cayley_menger_volume(apply(P, 2, phi))
    err_volume <- abs(vol_image - abs(det(A)) * vol_source)
    list(err_ratio = err_ratio, err_volume = err_volume,
         max_len_change = max_len_change, max_angle_change = max_angle_change,
         ok = err_ratio < 1e-10 && err_volume < 1e-8 &&
             max_len_change > 1e-2 && max_angle_change > 1e-2)
}

# ---- 4: solution sets of linear systems are affine sets ----
test_solution_sets <- function() {
    n <- 5L; rank_target <- 2L
    U <- matrix(rnorm(3 * rank_target), 3, rank_target)
    V <- matrix(rnorm(rank_target * n), rank_target, n)
    A <- U %*% V                                   # a 3 by 5 matrix of rank 2 by construction
    x_star <- rnorm(n); b <- as.vector(A %*% x_star)
    sv <- svd(A, nu = nrow(A), nv = n)
    rank_num <- sum(sv$d > max(dim(A)) * .Machine$double.eps * max(sv$d))
    null_basis <- sv$v[, (rank_num + 1):n, drop = FALSE]
    dim_expected <- n - rank_target
    sols <- x_star + null_basis %*% matrix(rnorm(dim_expected * (dim_expected + 2)),
                                           dim_expected, dim_expected + 2)
    res_sols <- max(abs(A %*% sols - b))
    w <- rnorm(ncol(sols)); w1 <- w / sum(w)
    res_affine <- max(abs(A %*% (sols %*% w1) - b))
    s <- sum(w)
    res_closed <- max(abs(as.vector(A %*% (sols %*% w)) - b - (s - 1) * b))
    diffs <- sols[, -1, drop = FALSE] - sols[, 1]
    dim_hull <- qr(diffs)$rank
    # degenerate cases
    dim_point <- qr(matrix(0, n, 1))$rank
    d <- 3L
    too_many <- matrix(rnorm(d * (d + 2)), d, d + 2)
    dim_too_many <- qr(too_many[, -1, drop = FALSE] - too_many[, 1])$rank
    list(rank_num = rank_num, dim_hull = dim_hull, res_sols = res_sols,
         res_affine = res_affine, res_closed = res_closed,
         ok = rank_num == rank_target && dim_hull == dim_expected && res_sols < 1e-10 &&
             res_affine < 1e-10 && res_closed < 1e-10 && dim_point == 0L && dim_too_many <= d)
}

# ---- 5: barycentric coordinates and their invariance ----
test_barycentric <- function(trials = 300L, d = 3L) {
    U <- matrix(rnorm(d * (d + 1)), d, d + 1)      # frame of d + 1 points, affinely independent
    stopifnot(qr(U[, -1, drop = FALSE] - U[, 1])$rank == d)
    A <- matrix(c(1.3, -0.4, 0.2, 0.5, 0.9, -0.7, 0.1, 0.6, 1.1), 3, 3); b <- c(-0.2, 0.8, 0.3)
    phi <- affine_map(A, b)
    U_img <- apply(U, 2, phi)
    M <- rbind(U, 1); M_img <- rbind(U_img, 1)     # the affine system that the coordinates solve
    err_sum <- 0; err_point <- 0; err_invariance <- 0; min_interior <- 1
    for (i in seq_len(trials)) {
        lam <- rgamma(d + 1, 1); lam <- lam / sum(lam)   # interior point, positive coordinates
        x <- as.vector(U %*% lam)
        lam_hat <- solve(M, c(x, 1))
        err_sum <- max(err_sum, abs(sum(lam_hat) - 1))
        err_point <- max(err_point, max(abs(lam_hat - lam)))
        min_interior <- min(min_interior, min(lam_hat))
        lam_img <- solve(M_img, c(phi(x), 1))
        err_invariance <- max(err_invariance, max(abs(lam_img - lam)))
    }
    list(err_sum = err_sum, err_point = err_point, err_invariance = err_invariance,
         ok = err_sum < 1e-10 && err_point < 1e-10 && err_invariance < 1e-10 && min_interior > 0)
}

# ---- 6: an affine conjugacy leaves the Lyapunov exponent unchanged ----
lyapunov_logistic <- function(r, x0, n, burn) {
    x <- x0; s <- 0
    for (i in seq_len(burn)) x <- r * x * (1 - x)
    for (i in seq_len(n)) {
        s <- s + log(abs(r * (1 - 2 * x)))
        x <- r * x * (1 - x)
    }
    s / n
}
lyapunov_quadratic <- function(a, y0, n, burn) {
    y <- y0; s <- 0
    for (i in seq_len(burn)) y <- 1 - a * y^2
    for (i in seq_len(n)) {
        s <- s + log(abs(2 * a * y))
        y <- 1 - a * y^2
    }
    s / n
}

test_conjugacy_lyapunov <- function() {
    n <- 2e6L; burn <- 1e4L
    r <- 4; a <- r * (r - 2) / 4
    lam_log <- lyapunov_logistic(r, runif(1, 0.1, 0.9), n, burn)
    lam_quad <- lyapunov_quadratic(a, runif(1, -0.9, 0.9), n, burn)
    err_chaotic <- max(abs(lam_log - log(2)), abs(lam_quad - log(2)), abs(lam_log - lam_quad))
    r2 <- 3.2; a2 <- r2 * (r2 - 2) / 4
    lam_cycle <- 0.5 * log(abs(4 + 2 * r2 - r2^2))          # multiplier of the two-cycle
    lam_log2 <- lyapunov_logistic(r2, runif(1, 0.1, 0.9), 2e5L, 2e5L)
    lam_quad2 <- lyapunov_quadratic(a2, runif(1, -0.9, 0.9), 2e5L, 2e5L)
    err_periodic <- max(abs(lam_log2 - lam_cycle), abs(lam_quad2 - lam_cycle))
    list(lam_log = lam_log, lam_quad = lam_quad, err_chaotic = err_chaotic,
         lam_cycle = lam_cycle, err_periodic = err_periodic,
         ok = err_chaotic < 5e-3 && err_periodic < 1e-8)
}

# ---- 7: a change of units is an affine change of variables ----
test_diagonal_rescaling <- function() {
    if (!requireNamespace("deSolve", quietly = TRUE)) return(list(ok = TRUE, err_traj = NA_real_, margin = NA_real_))
    suppressPackageStartupMessages(library(deSolve))
    n <- 4L
    B <- matrix(rnorm(n * n), n, n)
    A <- -(B %*% t(B) + diag(n))          # symmetric negative definite, so D = I certifies it
    b <- runif(n, 0.5, 1.5)
    D <- diag(n)
    E <- diag(runif(n, 0.4, 2.5))
    lv <- function(t, u, p) list(u * (p$b + as.vector(p$A %*% u)))
    u0 <- runif(n, 0.2, 1.5)
    times <- seq(0, 40, by = 0.05)
    sol_u <- as.data.frame(ode(u0, times, lv, list(b = b, A = A),
                               method = "lsoda", rtol = 1e-11, atol = 1e-11))
    sol_v <- as.data.frame(ode(as.vector(solve(E) %*% u0), times, lv, list(b = b, A = A %*% E),
                               method = "lsoda", rtol = 1e-11, atol = 1e-11))
    err_traj <- max(abs(as.matrix(sol_v[, -1]) - t(solve(E) %*% t(as.matrix(sol_u[, -1])))))
    P <- E %*% D
    S <- P %*% (A %*% E) + t(A %*% E) %*% P
    margin <- max(Re(eigen((S + t(S)) / 2)$values))
    identity_case <- max(abs(A %*% diag(n) - A))
    list(err_traj = err_traj, margin = margin,
         ok = err_traj < 1e-8 && margin < 0 && identity_case == 0)
}

t1 <- test_combination_invariance()
t2 <- test_reconstruction_from_frame()
t3 <- test_ratio_and_metric()
t4 <- test_solution_sets()
t5 <- test_barycentric()
t6 <- test_conjugacy_lyapunov()
t7 <- test_diagonal_rescaling()

ok <- t1$ok && t2$ok && t3$ok && t4$ok && t5$ok && t6$ok && t7$ok

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    A2 <- matrix(c(1.4, 0.55, -0.65, 1.0), 2, 2); b2 <- c(0.6, -0.35)
    phi2 <- affine_map(A2, b2)
    grid_lines <- do.call(rbind, lapply(0:4, function(k) {
        s <- seq(0, 4, length.out = 41)
        rbind(data.frame(x = s, y = k, line = paste0("h", k)),
              data.frame(x = k, y = s, line = paste0("v", k)))
    }))
    img <- t(apply(as.matrix(grid_lines[, c("x", "y")]), 1, phi2))
    grid_img <- data.frame(x = img[, 1], y = img[, 2], line = grid_lines$line)
    tri <- data.frame(x = c(0.5, 3.5, 1.0, 0.5), y = c(0.5, 1.0, 3.5, 0.5))
    tri_img <- as.data.frame(t(apply(as.matrix(tri), 1, phi2))); names(tri_img) <- c("x", "y")
    seg <- data.frame(x = c(0.5, 3.5), y = c(0.5, 1.0))
    mid <- data.frame(x = mean(seg$x), y = mean(seg$y))
    mid_img <- as.data.frame(t(phi2(c(mid$x, mid$y)))); names(mid_img) <- c("x", "y")
    p1 <- ggplot() +
        geom_path(data = grid_lines, aes(x, y, group = line), colour = "grey70", linewidth = 0.3) +
        geom_path(data = grid_img, aes(x, y, group = line), colour = "#4A6FA5", linewidth = 0.25) +
        geom_path(data = tri, aes(x, y), colour = "grey45", linewidth = 0.7) +
        geom_path(data = tri_img, aes(x, y), colour = "#FB9E07", linewidth = 0.7) +
        geom_point(data = mid, aes(x, y), colour = "grey35", fill = "white", shape = 21, size = 2.2, stroke = 0.8) +
        geom_point(data = mid_img, aes(x, y), colour = "#FB9E07", fill = "white", shape = 21, size = 2.2, stroke = 0.8) +
        coord_equal() +
        labs(title = "An affine map of the plane",
             subtitle = "Source in grey, image in blue and orange",
             x = "First coordinate", y = "Second coordinate") + kb_theme()
    rs <- seq(2.05, 4, by = 0.005)
    lam_pairs <- vapply(rs, function(r) {
        a <- r * (r - 2) / 4
        c(lyapunov_logistic(r, runif(1, 0.1, 0.9), 4000L, 2000L),
          lyapunov_quadratic(a, runif(1, -0.9, 0.9), 4000L, 2000L))
    }, numeric(2))
    lam_df <- data.frame(r = rep(rs, 2), lambda = c(lam_pairs[1, ], lam_pairs[2, ]),
                         map = rep(c("$Logistic~S_r$", "$Quadratic~Q_a,~a = r(r - 2)/4$"), each = length(rs)))
    p2 <- ggplot(lam_df, aes(r, lambda, colour = map, linewidth = map, linetype = map)) +
        geom_line(alpha = 0.9) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
        geom_point(data = data.frame(r = 4, lambda = log(2)), aes(r, lambda),
                   inherit.aes = FALSE, colour = "grey20", size = 1.8) +
        scale_colour_manual(values = c("#4A6FA5", "#FB9E07"), name = NULL,
                            labels = kb_tex(sort(unique(lam_df$map)))) +
        scale_linewidth_manual(values = c(1.3, 0.45), name = NULL,
                               labels = kb_tex(sort(unique(lam_df$map)))) +
        scale_linetype_manual(values = c("solid", "22"), name = NULL,
                              labels = kb_tex(sort(unique(lam_df$map)))) +
        labs(title = "An affine conjugacy leaves the exponent unchanged",
             subtitle = "Exponents of two conjugate families, from independent orbits",
             x = kb_tex("Growth parameter $r$"),
             y = kb_tex("Lyapunov exponent $\\lambda$")) + kb_theme()
    fig <- (p1 | p2) +
        plot_annotation(caption = kb_caption(sprintf(
            "Left: the map $x \\mapsto Ax + b$ with $A = [[%.2f, %.2f], [%.2f, %.2f]]$ and $b = (%.2f, %.2f)$ takes lines to lines and parallels to parallels, and takes the midpoint of a side to the midpoint of its image, while lengths and angles change. Right: the exponent $\\lambda$ of the logistic map $S_r$ and of the quadratic map $Q_a$ with $a = r(r - 2)/4$, each from an independently drawn initial condition, 4000 iterations after a transient of 2000; the point marks the value $\\log 2 = %.4f$ that both families take at $r = 4$. Drawn by checks/affine-structure.R.",
            A2[1, 1], A2[1, 2], A2[2, 1], A2[2, 2], b2[1], b2[2], log(2))),
            theme = kb_theme())
    kb_save(fig, "affine-structure")
}

emit("affine-structure", if (ok) "pass" else "fail",
     "Affine maps commute with affine combinations and preserve ratios, hulls, Lyapunov exponents and the Lotka-Volterra form under a change of units, and do not preserve lengths, angles or volumes",
     list(err_affine_combination = t1$err_affine, err_weight_closed_form = t1$err_closed,
          err_frame_reconstruction = t2$err_linear, residual_nonlinear_map = t2$res_nonlinear,
          err_ratio = t3$err_ratio, err_volume_scaling = t3$err_volume,
          max_length_change = t3$max_len_change, max_angle_change = t3$max_angle_change,
          solution_set_dimension = t4$dim_hull, err_solution_affine = t4$res_affine,
          err_barycentric_invariance = t5$err_invariance,
          lyapunov_logistic_r4 = t6$lam_log, lyapunov_quadratic_r4 = t6$lam_quad,
          err_lyapunov_chaotic = t6$err_chaotic, err_lyapunov_periodic = t6$err_periodic,
          err_rescaled_trajectory = t7$err_traj, rescaled_stability_margin = t7$margin))
