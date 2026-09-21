# ============================ #
# Equilibria of a Lotka-Volterra system as a linear complementarity problem ####
# ============================ #
#
# The generalised Lotka-Volterra system is du_i/dt = u_i (b_i + sum_j a_ij u_j) on the nonnegative
# orthant. A saturated equilibrium is a nonnegative u with b_i + (A u)_i = 0 for every species
# present and b_i + (A u)_i <= 0 for every species absent, so that no absent species can invade.
# Those two conditions are the linear complementarity problem LCP(-A, -b): find u >= 0 with
# -A u - b >= 0 and u' (-A u - b) = 0. When A is Volterra-Lyapunov stable, meaning that there is a
# positive diagonal H with H A + A' H negative definite, the problem has a unique solution and that
# equilibrium attracts every trajectory starting with the same set of species present; the source
# of the knowledge base states this as Theorem 5, citing Takeuchi (1996).
#
# The script tests five claims and prints one JSON line for kb.py verify.
#
# 1. Uniqueness. For a Volterra-Lyapunov stable A built as A = -(M M' + I), which is symmetric and
#    negative definite so that H = I works, the script enumerates all 2^n subsets of species,
#    solves the linear system on each support and keeps the solutions that are positive on the
#    support and non-invasible off it. Exactly one such vector must exist.
# 2. Complementarity. The residuals max |u_i (A u + b)_i| and the violations of u >= 0 and
#    A u + b <= 0 are reported for that vector.
# 3. Dynamics. The system is integrated with deSolve from twenty random positive initial
#    conditions, and every trajectory must reach the vector of test 1. This is the independent
#    route: one side solves a linear complementarity problem, the other integrates the equations,
#    and the two share no code.
# 4. Extinction is predicted, not assumed. A second system is built whose unconstrained solution
#    of A u = -b has negative entries. The LCP solution must then have zeros exactly where the
#    dynamics drives species to extinction, which the script checks against the integration.
# 5. Degenerate cases. With one species the equilibrium is -b/a when positive and zero otherwise,
#    and the script checks both. With b = 0 the only saturated equilibrium is the origin.
#
# Tolerances. The linear algebra is compared at 1e-10, far above the conditioning of these small
# systems and far below the separation between the equilibria. The integration is compared at 1e-6
# after a long run, which is above the solver tolerance of 1e-10 and below the smallest positive
# component of the equilibria used. A species is called extinct in the integration when its
# abundance falls below 1e-8, which no surviving component of these systems approaches.

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

# ---- the model ----
glv_rhs <- function(t, u, parms) list(u * (parms$b + as.vector(parms$A %*% u)))

integrate_glv <- function(A, b, u0, times = c(0, 10^seq(0, 4, length.out = 60))) {
    out <- ode(y = u0, times = times, func = glv_rhs, parms = list(A = A, b = b),
               method = "lsoda", rtol = 1e-10, atol = 1e-12)
    as.vector(out[nrow(out), -1])
}

# Saturated equilibria by enumeration of the supports: the definition, solved directly.
saturated_equilibria <- function(A, b, tol = 1e-10) {
    n <- length(b)
    found <- list()
    for (mask in 0:(2^n - 1)) {
        support <- which(bitwAnd(mask, 2^(seq_len(n) - 1)) > 0)
        u <- rep(0, n)
        if (length(support) > 0) {
            sub <- A[support, support, drop = FALSE]
            if (abs(det(sub)) < 1e-12) next
            u[support] <- solve(sub, -b[support])
            if (any(u[support] <= tol)) next
        }
        rates <- b + as.vector(A %*% u)
        if (any(abs(rates[support]) > tol)) next
        outside <- setdiff(seq_len(n), support)
        if (length(outside) > 0 && any(rates[outside] > tol)) next
        found[[length(found) + 1L]] <- u
    }
    found
}

complementarity_residual <- function(A, b, u) {
    rates <- b + as.vector(A %*% u)
    c(product = max(abs(u * rates)), negative_abundance = max(0, -min(u)),
      invasion = max(0, max(rates)))
}

volterra_lyapunov <- function(A, H = diag(nrow(A))) max(eigen(H %*% A + t(A) %*% H, symmetric = TRUE)$values)

# ---- tests 1 to 3: a system whose interior equilibrium is feasible ----
n <- 5
M <- matrix(rnorm(n * n), n, n)
A <- -(M %*% t(M) + diag(n))                       # symmetric negative definite, so VL-stable
b <- as.vector(-A %*% runif(n, 0.5, 1.5))          # chosen so that the interior solution is positive
vl_margin <- volterra_lyapunov(A)                  # must be negative
eq <- saturated_equilibria(A, b)
n_equilibria <- length(eq)
u_star <- if (n_equilibria == 1) eq[[1]] else rep(NA_real_, n)
residuals <- complementarity_residual(A, b, u_star)
unconstrained <- as.vector(solve(A, -b))
interior_gap <- max(abs(u_star - unconstrained))   # equal when the unconstrained solution is positive

reached <- vapply(seq_len(20), function(k) {
    u0 <- runif(n, 0.05, 3)
    max(abs(integrate_glv(A, b, u0) - u_star))
}, numeric(1))
dynamics_error <- max(reached)

# ---- test 4: a system in which the complementarity problem predicts extinctions ----
b_drop <- b
b_drop[c(2, 4)] <- -abs(b_drop[c(2, 4)]) - 1       # make two species unable to persist
unconstrained_drop <- as.vector(solve(A, -b_drop))
eq_drop <- saturated_equilibria(A, b_drop)
u_drop <- eq_drop[[1]]
predicted_extinct <- which(u_drop <= 1e-12)
final_drop <- integrate_glv(A, b_drop, runif(n, 0.05, 3))
observed_extinct <- which(final_drop < 1e-8)
extinction_agrees <- length(eq_drop) == 1 && identical(predicted_extinct, observed_extinct)
drop_error <- max(abs(final_drop - u_drop))
residuals_drop <- complementarity_residual(A, b_drop, u_drop)

# ---- test 5: degenerate cases ----
one_positive <- saturated_equilibria(matrix(-2, 1, 1), 3)[[1]]      # -b/a = 1.5
one_negative <- saturated_equilibria(matrix(-2, 1, 1), -3)[[1]]     # no positive equilibrium
zero_growth <- saturated_equilibria(A, rep(0, n))
degenerate_error <- max(abs(one_positive - 1.5), abs(one_negative - 0),
                        abs(zero_growth[[1]]))

ok <- vl_margin < 0 && n_equilibria == 1 && max(residuals) < 1e-10 && interior_gap < 1e-10 &&
    dynamics_error < 1e-6 && extinction_agrees && drop_error < 1e-6 &&
    max(residuals_drop) < 1e-10 && length(zero_growth) == 1 && degenerate_error < 1e-12

# ---- figure ----
if (requireNamespace("ggplot2", quietly = TRUE) && file.exists("checks/lib/figure-style.R")) {
    source("checks/lib/figure-style.R")
    library(ggplot2); library(patchwork)
    times <- c(0, 10^seq(-1, 2, length.out = 300))
    traj <- do.call(rbind, lapply(1:3, function(k) {
        u0 <- runif(n, 0.05, 3)
        out <- as.data.frame(ode(y = u0, times = times, func = glv_rhs,
                                 parms = list(A = A, b = b_drop), method = "lsoda",
                                 rtol = 1e-10, atol = 1e-12))
        names(out) <- c("time", paste0("s", seq_len(n)))
        long <- reshape(out, direction = "long", varying = paste0("s", seq_len(n)),
                        v.names = "abundance", timevar = "species", times = seq_len(n))
        long$run <- k
        long
    }))
    traj$species <- factor(traj$species)
    ends <- data.frame(species = factor(seq_len(n)), value = u_drop)
    p1 <- ggplot(traj, aes(time, abundance, colour = species, group = interaction(species, run))) +
        geom_line(linewidth = 0.4, alpha = 0.85) +
        geom_hline(data = ends, aes(yintercept = value, colour = species),
                   linetype = "dashed", linewidth = 0.4) +
        scale_x_log10() +
        labs(title = "Three runs, one saturated equilibrium",
             subtitle = sprintf("Five species; the complementarity problem sets %d of them to zero",
                                length(predicted_extinct)),
             x = kb_tex("Time $t$"), y = kb_tex("Abundance $u_i$")) + kb_theme() + theme(legend.position = "none")
    bars <- data.frame(species = rep(factor(seq_len(n)), 2),
                       value = c(unconstrained_drop, u_drop),
                       route = rep(c("Linear solve", "Complementarity"), each = n))
    p2 <- ggplot(bars, aes(species, value, fill = route)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
        scale_fill_manual(values = c("#4A6FA5", "#FB9E07"), name = NULL) +
        labs(title = "Solving the linear system is not enough",
             subtitle = "The linear solve returns negative abundances",
             x = kb_tex("Species $i$"), y = kb_tex("Abundance $u_i^{*}$ at equilibrium"),
             caption = kb_caption(paste(
                 "Left: three trajectories of $du_i/dt = u_i(b_i + (Au)_i)$ with a Volterra-Lyapunov stable $A$,",
                 "from random positive starts, against the equilibrium that the linear complementarity problem",
                 "returns (dashed). Right: the same problem solved as a linear system, which puts two species at",
                 "negative abundance, beside the complementarity solution, which sets them to zero and leaves the",
                 "rest non-invasible. Source: checks/lotka-volterra-lcp-equilibria.R"))) +
        kb_theme() + theme(legend.position = "top")
    kb_save(p1 + p2, "lotka-volterra-lcp-equilibria", width = 9.4, height = 4.3)
}

emit("lotka-volterra-lcp-equilibria", if (ok) "pass" else "fail",
     sprintf("A Volterra-Lyapunov stable system has exactly one saturated equilibrium, which the complementarity conditions give to %.0e and which twenty random trajectories reach to %.0e",
             max(residuals), dynamics_error),
     list(volterra_lyapunov_margin = vl_margin, saturated_equilibria = n_equilibria,
          complementarity_residual = unname(residuals["product"]),
          invasion_violation = unname(residuals["invasion"]),
          interior_gap_against_linear_solve = interior_gap,
          dynamics_error = dynamics_error, extinctions_predicted = length(predicted_extinct),
          extinctions_agree = extinction_agrees, extinction_dynamics_error = drop_error,
          degenerate_error = degenerate_error))
