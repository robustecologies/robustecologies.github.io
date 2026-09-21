# ============================ #
# Koopman modes, generalized Laplace analysis and dynamic mode decomposition on linear systems ####
# ============================ #
#
# Claims, with the locators of budisic2012 (arXiv preprint numbering):
#
# 1. Example 20. For the harmonic oscillator dp1/dt = p2, dp2/dt = -w^2 p1 and the
#    observable F(p) = p, F(p(t)) = exp(i w t) phi1(p0) C1 + exp(-i w t) phi2(p0) C2 with
#    phi1 = (p1 - i p2 / w) / 2, phi2 = (p1 + i p2 / w) / 2, C1 = (1, i w), C2 = (1, -i w).
#    Reference: the closed-form flow matrix, on 200 random initial conditions and times.
#    The in-text formula of the preprint, (p1 -/+ p2 / w) / 2 without the factor i, is
#    evaluated as well and must fail, which documents the misprint in the source note.
# 2. Eq. (44). For the time-h map T of the oscillator, the average
#    (1/K) sum_{k<K} exp(-i w h k) F(T^k p0) equals phi1(p0) C1 plus the remainder
#    phi2(p0) C2 (1 - q^K) / (K (1 - q)) with q = exp(-2 i w h). Reference: this closed form
#    of the finite geometric sum, and its bound 2 |phi2| |C2| / (K |1 - q|).
# 3. Thm. 17. For A = P diag(1.25, 0.8) P^{-1} and F(x) = x, the modes are
#    u_j = phi_j(x) v_j with (phi_1(x), phi_2(x)) = P^{-1} x and v_j the columns of P.
#    The average with weights 1.25^{-k} equals u1 plus u2 (1 - 0.64^K) / (0.36 K), a closed
#    form; after subtracting the exact mode 1.25^k u1, the average with weights 0.8^{-k}
#    equals u2 for every K; without the subtraction it grows without bound.
# 4. Def. 21. Dynamic mode decomposition with r snapshots of a linear map whose data lie in
#    the span of r eigenfunction components: the least-squares residual vanishes, the
#    empirical Ritz values are the eigenvalues, and when the left eigenvectors of the
#    companion matrix are scaled to have first entry 1, the empirical Ritz vectors are
#    phi_j(p) C_j(F). References: the eigenvalues used to build the map and the exact
#    projections P^{-1} p. Tested on a real 4 x 4 map and on the oscillator time-h map.
# 5. Degenerate case, eq. (58) and Remark 22. For the map p -> 0.5 p (p <= 0), p -> 2 p
#    (p > 0), phi1 = min(p, 0) and phi2 = max(p, 0) are eigenfunctions with eigenvalues 0.5
#    and 2, their pointwise product is identically zero, and dynamic mode decomposition with
#    r = 1 from p0 = -1 returns only 0.5 and from p0 = 1 only 2.

set.seed(20260914L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

cnorm <- function(v) sqrt(sum(Mod(v)^2))

# ============================ #
# Test 1: Koopman modes of the harmonic oscillator ####
# ============================ #

w <- 1.7
flow <- function(t) matrix(c(cos(w * t), -w * sin(w * t), sin(w * t) / w, cos(w * t)), 2, 2)
C1 <- c(1, 1i * w)
C2 <- c(1, -1i * w)
phi1 <- function(p) (p[1] - 1i * p[2] / w) / 2
phi2 <- function(p) (p[1] + 1i * p[2] / w) / 2

err1 <- 0
imag1 <- 0
err1_typo <- 0
for (s in seq_len(200L)) {
    p0 <- runif(2, -1, 1)
    t <- runif(1, 0, 10)
    exact <- as.vector(flow(t) %*% p0)
    expansion <- exp(1i * w * t) * phi1(p0) * C1 + exp(-1i * w * t) * phi2(p0) * C2
    err1 <- max(err1, Mod(expansion - exact))
    imag1 <- max(imag1, abs(Im(expansion)))
    typo <- exp(1i * w * t) * (p0[1] - p0[2] / w) / 2 * C1 + exp(-1i * w * t) * (p0[1] + p0[2] / w) / 2 * C2
    err1_typo <- max(err1_typo, Mod(typo - exact))
}
# Entries are of order 1 + w, so rounding stays near 1e-15; 1e-12 leaves a wide margin.
# The misprinted formula differs by sqrt(2) |p2 sin(w t)| / w in the first component,
# which exceeds 0.1 for most of the 200 draws.
pass1 <- err1 < 1e-12 && imag1 < 1e-12 && err1_typo > 0.1

# ============================ #
# Test 2: Fourier average for an eigenvalue on the unit circle ####
# ============================ #

h <- 0.3
Th <- flow(h)
p0 <- c(0.7, -0.4)
Ks <- c(10L, 100L, 1000L, 10000L)
q <- exp(-2i * w * h)
target2 <- phi1(p0) * C1
x <- p0
acc <- c(0 + 0i, 0 + 0i)
err2_closed <- 0
ratio2 <- 0
err2_last <- NA_real_
for (k in 0:(max(Ks) - 1L)) {
    acc <- acc + exp(-1i * w * h * k) * x
    x <- as.vector(Th %*% x)
    K <- k + 1L
    if (K %in% Ks) {
        avg <- acc / K
        remainder <- phi2(p0) * C2 * (1 - q^K) / (K * (1 - q))
        err2_closed <- max(err2_closed, cnorm(avg - target2 - remainder))
        bound <- Mod(phi2(p0)) * cnorm(C2) * 2 / (K * Mod(1 - q))
        ratio2 <- max(ratio2, cnorm(avg - target2) / bound)
        if (K == max(Ks)) err2_last <- cnorm(avg - target2)
    }
}
# 10^4 matrix products and summands of order 1, each with relative rounding near 2e-16,
# accumulate to well below 1e-10; 1e-9 leaves a margin. The bound ratio may not exceed 1
# beyond rounding.
pass2 <- err2_closed < 1e-9 && ratio2 <= 1 + 1e-6 && err2_last < 1e-3

# ============================ #
# Test 3: generalized Laplace analysis with an unstable mode ####
# ============================ #

P <- matrix(c(1, 0.5, 1, -1), 2, 2)
lam <- c(1.25, 0.8)
A <- P %*% diag(lam) %*% solve(P)
x0 <- c(0.3, 0.9)
proj <- solve(P, x0)
u1 <- proj[1] * P[, 1]
u2 <- proj[2] * P[, 2]

K3a <- 200L
x <- x0
acc <- c(0, 0)
err3a <- 0
for (k in 0:(K3a - 1L)) {
    acc <- acc + lam[1]^(-k) * x
    x <- as.vector(A %*% x)
    K <- k + 1L
    closed <- u1 + u2 * (1 - (lam[2] / lam[1])^K) / (K * (1 - lam[2] / lam[1]))
    err3a <- max(err3a, cnorm(acc / K - closed))
}

K3b <- 30L
x <- x0
acc_sub <- c(0, 0)
acc_raw <- c(0, 0)
err3b <- 0
for (k in 0:(K3b - 1L)) {
    acc_sub <- acc_sub + lam[2]^(-k) * (x - lam[1]^k * u1)
    acc_raw <- acc_raw + lam[2]^(-k) * x
    x <- as.vector(A %*% x)
    K <- k + 1L
    err3b <- max(err3b, cnorm(acc_sub / K - u2))
}
raw3b <- cnorm(acc_raw / K3b)
# Weights 1.25^{-k} damp the data, so rounding stays near 1e-14 (tolerance 1e-9). With
# weights 0.8^{-k} the rounding error of A^k x, of order k eps 1.25^k, is amplified by
# 0.8^{-k}; at k = 29 this gives about 1e-8, hence the tolerance 1e-6. The unsubtracted
# average is about |u1| 1.5625^30 / (0.5625 * 30), of order 3e4.
pass3 <- err3a < 1e-9 && err3b < 1e-6 && raw3b > 1e3

# ============================ #
# Test 4: dynamic mode decomposition on exact data ####
# ============================ #

dmd <- function(B) {
    r <- ncol(B) - 1L
    Br <- B[, seq_len(r), drop = FALSE]
    br <- B[, r + 1L]
    cvec <- qr.solve(Br, br)
    eta <- br - as.vector(Br %*% cvec)
    Ar <- matrix(0, r, r)
    if (r > 1L) Ar[cbind(2:r, 1:(r - 1L))] <- 1
    Ar[, r] <- cvec
    eg <- eigen(Ar)
    V <- solve(eg$vectors)
    Vs <- V / V[, 1]
    Z <- Br %*% solve(Vs)
    list(values = eg$values, vectors = Z, residual = cnorm(eta) / cnorm(br), condition = kappa(Br, exact = TRUE))
}

lam4 <- c(0.95, 0.6, -0.4, 0.2)
P4 <- diag(4) + 0.3 * matrix(rnorm(16), 4, 4)
A4 <- P4 %*% diag(lam4) %*% solve(P4)
p4 <- rnorm(4)
B4 <- matrix(0, 4, 5)
B4[, 1] <- p4
for (k in 1:4) B4[, k + 1L] <- as.vector(A4 %*% B4[, k])
res4 <- dmd(B4)
proj4 <- solve(P4, p4)
U4 <- sweep(P4, 2, proj4, `*`)
err4_val <- 0
err4_vec <- 0
for (j in seq_along(lam4)) {
    idx <- which.min(Mod(res4$values - lam4[j]))
    err4_val <- max(err4_val, Mod(res4$values[idx] - lam4[j]))
    err4_vec <- max(err4_vec, cnorm(res4$vectors[, idx] - U4[, j]) / max(apply(U4, 2, cnorm)))
}

Bo <- matrix(0, 2, 3)
Bo[, 1] <- p0
for (k in 1:2) Bo[, k + 1L] <- as.vector(Th %*% Bo[, k])
reso <- dmd(Bo)
mu_o <- c(exp(1i * w * h), exp(-1i * w * h))
modes_o <- cbind(phi1(p0) * C1, phi2(p0) * C2)
err4o_val <- 0
err4o_vec <- 0
for (j in 1:2) {
    idx <- which.min(Mod(reso$values - mu_o[j]))
    err4o_val <- max(err4o_val, Mod(reso$values[idx] - mu_o[j]))
    err4o_vec <- max(err4o_vec, cnorm(reso$vectors[, idx] - modes_o[, j]))
}
# The data are exact up to rounding, so the errors scale with machine epsilon times the
# condition number of the snapshot matrix (reported as a metric, of order 10^2 to 10^3);
# 1e-8 leaves a margin of several orders of magnitude.
pass4 <- res4$residual < 1e-10 && err4_val < 1e-8 && err4_vec < 1e-8 &&
    reso$residual < 1e-10 && err4o_val < 1e-8 && err4o_vec < 1e-8

# ============================ #
# Test 5: degenerate piecewise linear map of eq. (58) ####
# ============================ #

l1 <- 0.5
l2 <- 2
Tpl <- function(p) ifelse(p <= 0, l1 * p, l2 * p)
grid <- seq(-2, 2, length.out = 401)
err5_eig <- max(abs(pmin(Tpl(grid), 0) - l1 * pmin(grid, 0)), abs(pmax(Tpl(grid), 0) - l2 * pmax(grid, 0)))
prod5 <- max(abs(pmin(grid, 0) * pmax(grid, 0)))
neg <- dmd(matrix(c(-1, Tpl(-1)), 1, 2))
pos <- dmd(matrix(c(1, Tpl(1)), 1, 2))
err5_dmd <- max(Mod(neg$values - l1), Mod(neg$vectors[1, 1] - min(-1, 0)), Mod(pos$values - l2), Mod(pos$vectors[1, 1] - max(1, 0)))
# Every operation here is exact in binary floating point (multiplication by 0.5 and 2).
pass5 <- err5_eig < 1e-15 && prod5 == 0 && err5_dmd < 1e-15

status <- if (pass1 && pass2 && pass3 && pass4 && pass5) "pass" else "fail"
emit("koopman-mode-linear-systems", status,
     "Koopman modes of the harmonic oscillator, generalized Laplace analysis and dynamic mode decomposition reproduce the exact modes of linear systems",
     list(oscillator_mode_error = err1, oscillator_imaginary_part = imag1, misprinted_formula_error = err1_typo,
          fourier_average_closed_form_error = err2_closed, fourier_average_bound_ratio = ratio2,
          fourier_average_error_K10000 = err2_last,
          laplace_unstable_mode_error = err3a, laplace_stable_mode_error = err3b, laplace_unsubtracted_norm = raw3b,
          dmd_residual = res4$residual, dmd_value_error = err4_val, dmd_vector_error = err4_vec,
          dmd_snapshot_condition = res4$condition,
          dmd_oscillator_value_error = err4o_val, dmd_oscillator_vector_error = err4o_vec,
          piecewise_eigen_error = err5_eig, piecewise_product_max = prod5, piecewise_dmd_error = err5_dmd))

# ============================ #
# Figure ####
# ============================ #
#
# The generalized Laplace average converging to a Koopman mode, in the two systems of the test.
# The averages are recomputed here with the same maps and the same starting points that the
# check used, and each is drawn against the closed-form bound on its error, which falls like
# one over the number of terms.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    Kmax <- 2000L
    # The oscillator: averaging against the conjugate eigenvalue leaves the other mode as a
    # remainder of size 2 |phi2(p0)| ||C2|| / (K |1 - q|).
    x <- p0; acc <- c(0 + 0i, 0 + 0i); errs_osc <- numeric(Kmax)
    for (k in 0:(Kmax - 1L)) {
        acc <- acc + exp(-1i * w * h * k) * x
        x <- as.vector(Th %*% x)
        errs_osc[k + 1L] <- cnorm(acc / (k + 1L) - target2)
    }
    bound_osc <- Mod(phi2(p0)) * cnorm(C2) * 2 / (seq_len(Kmax) * Mod(1 - q))
    # The unstable linear map: averaging against the larger eigenvalue leaves the smaller mode.
    x <- x0; acc2 <- c(0, 0); errs_lin <- numeric(Kmax)
    for (k in 0:(Kmax - 1L)) {
        acc2 <- acc2 + lam[1]^(-k) * x
        x <- as.vector(A %*% x)
        errs_lin[k + 1L] <- cnorm(acc2 / (k + 1L) - u1)
    }
    df <- rbind(
        data.frame(K = seq_len(Kmax), err = errs_osc, what = "Harmonic oscillator"),
        data.frame(K = seq_len(Kmax), err = errs_lin, what = "Linear map with an unstable mode"))
    bounds <- data.frame(K = seq_len(Kmax), err = bound_osc, what = "Harmonic oscillator")
    fig <- ggplot(df, aes(K, err, colour = what)) +
        geom_line(data = bounds, aes(K, err), inherit.aes = FALSE,
                  colour = "grey55", linewidth = 0.4, linetype = "22") +
        geom_line(linewidth = 0.6) +
        annotate("text", x = 3, y = bound_osc[3] * 1.6, hjust = 0, size = 2.6, colour = "grey40",
                 label = kb_unicode("Grey: the closed-form bound on the remainder, which falls like $1/K$")) +
        scale_x_log10() + scale_y_log10() +
        scale_colour_viridis_d(option = "viridis", end = 0.7, name = NULL) +
        labs(title = "Recovering a Koopman mode from one trajectory",
             subtitle = "Error of the generalized Laplace average against the number of terms, for two linear systems",
             x = kb_tex("Number of terms $K$ in the average"), y = "Distance from the exact mode") +
        kb_theme() +
        labs(caption = kb_caption(sprintf(
            "The average of an observable along one orbit, weighted by the conjugate of an eigenvalue, converges to the Koopman mode of that eigenvalue, and what remains at $K$ terms is the contribution of the other mode divided by $K$. In the run recorded by checks/koopman-mode-linear-systems.R the oscillator average matched its closed form to %.1e and stayed within the bound, and the mode of the unstable eigenvalue was recovered to %.1e after %d terms.",
            err2_closed, err3a, K3a)))
    invisible(kb_save(fig, "koopman-mode-linear-systems", width = 8.6, height = 4.6))
}
