# ============================ #
# Chain recurrent set of the circle maps of Norton's example ####
# ============================ #
#
# Claim. For the circle maps T(theta) = theta + alpha * cos(k * theta) with 0 < k*alpha < 1,
# the chain recurrent set at jump size eps < alpha is the open set where the displacement is
# smaller than the jump, {theta : alpha * |cos(k*theta)| < eps}, whose measure is
# 4 * asin(eps/alpha) for every k. For T(theta) = theta + alpha * cos^2(theta) the displacement
# never changes sign, every point returns to itself by going once round the circle, and the
# chain recurrent set is the whole circle at every eps > 0, although only two points are
# periodic. Norton states the two qualitative conclusions for cos^2 and cos in his Example 1
# (norton1995, sec. 4 and sec. 5); the closed form and the case k = 3 are derived below.
#
# Method. The circle is cut into N cells and an eps-transition graph is built: cell j points to
# cell l when some x in cell j and some y in cell l satisfy d(T(x), y) <= eps. The map is an
# increasing homeomorphism of the lift, so the image of a cell is an interval and the out-
# neighbours of a cell form a contiguous arc, which is stored as a pair of lifted indices. The
# chain recurrent cells are the cells that lie on a cycle of this graph, found with an iterative
# Tarjan scan written out below.
#
# Independent references, none of which uses the graph.
#   1. The closed form 4 * asin(eps/alpha), derived from the displacement, bracketed for the
#      discretisation as stated at test_measure_against_closed_form.
#   2. A constructive witness: for a cell the graph calls recurrent, an explicit eps-pseudo-orbit
#      is built from the map itself and every step d(T(x_i), x_{i+1}) < eps is verified.
#   3. A barrier witness: for a cell the graph calls transient, an arc is exhibited on which the
#      displacement exceeds eps with one sign, so no eps-pseudo-orbit crosses it backwards and no
#      return is possible.
#   4. The identity T([-pi/2, pi/2]) = [-pi/2, pi/2] for the maps with k = 1, which gives the
#      omega limit set of that half-circle exactly and separates it from the union of the omega
#      limit sets of its points (norton1995, sec. 4).
#
# Degenerate cases: eps >= alpha, where the two arcs merge and the whole circle is recurrent, and
# alpha = 0, the identity map, whose every point is fixed.

set.seed(20260919L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

TWO_PI <- 2 * pi

# ============================ #
# The systems ####
# ============================ #
#
# Each system carries the displacement T(theta) - theta, its Lipschitz constant, and the name of
# the closed form that the chain recurrent set is compared against.

systems <- list(
    a = list(label = "$T(\\theta) = \\theta + \\alpha\\cos^2\\theta$", alpha = 0.4,
             disp = function(th, alpha) alpha * cos(th)^2,
             lip  = function(alpha) 1 + alpha,          # |1 - alpha sin(2 theta)|
             form = "circle"),
    b = list(label = "$T(\\theta) = \\theta + \\alpha\\cos\\theta$", alpha = 0.4,
             disp = function(th, alpha) alpha * cos(th),
             lip  = function(alpha) 1 + alpha,
             form = "arcsin", k = 1),
    c = list(label = "$T(\\theta) = \\theta + \\alpha\\cos 3\\theta$", alpha = 0.25,
             disp = function(th, alpha) alpha * cos(3 * th),
             lip  = function(alpha) 1 + 3 * alpha,
             form = "arcsin", k = 3)
)

lift_map <- function(sys) function(th) th + sys$disp(th, sys$alpha)

# ============================ #
# The eps-transition graph ####
# ============================ #
#
# Cell j (j = 0, ..., N-1) is the closed interval [j*h, (j+1)*h] with h = 2*pi/N. The map lifts
# to an increasing function, so T(cell j) = [Tlift(j*h), Tlift((j+1)*h)] and the cells within eps
# of that interval are the cells meeting [Tlift(j*h) - eps, Tlift((j+1)*h) + eps]. Their lifted
# indices run from floor((Tlift(a_j) - eps)/h) to ceiling((Tlift(b_j) + eps)/h) - 1, which is an
# outer approximation by at most one cell on each side; the bracket of the measure test allows
# for it.

transition_arcs <- function(sys, N, eps) {
    h <- TWO_PI / N
    a <- (0:(N - 1)) * h
    Tl <- lift_map(sys)
    lo <- floor((Tl(a) - eps) / h)
    hi <- ceiling((Tl(a + h) + eps) / h) - 1
    list(lo = as.integer(lo), hi = as.integer(hi), h = h, N = N, eps = eps)
}

# Iterative Tarjan over the arc adjacency. The neighbours of a cell are generated on demand, so
# no edge list is ever materialised. Components are numbered in reverse topological order, that
# is component 1 is a sink of the condensation, which the Lyapunov construction relies on.
strong_components <- function(g) {
    N <- g$N; lo <- g$lo; hi <- g$hi
    index <- integer(N); low <- integer(N); onstack <- logical(N); comp <- integer(N)
    idx <- 0L; ncomp <- 0L
    stack <- integer(N); sp <- 0L
    callv <- integer(N + 1L); callk <- integer(N + 1L)
    for (root in seq_len(N)) {
        if (index[root] != 0L) next
        cp <- 1L; callv[1L] <- root; callk[1L] <- 0L
        idx <- idx + 1L; index[root] <- idx; low[root] <- idx
        sp <- sp + 1L; stack[sp] <- root; onstack[root] <- TRUE
        while (cp >= 1L) {
            v <- callv[cp]
            k <- callk[cp] + 1L
            if (k <= hi[v] - lo[v] + 1L) {
                callk[cp] <- k
                w <- ((lo[v] + k - 1L) %% N) + 1L   # cell of the lifted index lo[v] + k - 1
                if (index[w] == 0L) {
                    idx <- idx + 1L; index[w] <- idx; low[w] <- idx
                    sp <- sp + 1L; stack[sp] <- w; onstack[w] <- TRUE
                    cp <- cp + 1L; callv[cp] <- w; callk[cp] <- 0L
                } else if (onstack[w] && index[w] < low[v]) {
                    low[v] <- index[w]
                }
            } else {
                if (low[v] == index[v]) {
                    ncomp <- ncomp + 1L
                    repeat {
                        w <- stack[sp]; sp <- sp - 1L; onstack[w] <- FALSE
                        comp[w] <- ncomp
                        if (w == v) break
                    }
                }
                cp <- cp - 1L
                if (cp >= 1L && low[v] < low[callv[cp]]) low[callv[cp]] <- low[v]
            }
        }
    }
    comp
}

# A cell lies on a cycle when its component has more than one cell, or when it has a self-loop.
recurrent_cells <- function(g, comp) {
    N <- g$N
    sizes <- tabulate(comp, nbins = max(comp))
    idx <- seq_len(N)
    self <- (idx - 1L) >= g$lo & (idx - 1L) <= g$hi          # the lifted index of cell i is i-1
    wrap <- (idx - 1L + N) >= g$lo & (idx - 1L + N) <= g$hi  # the same cell one turn ahead
    back <- (idx - 1L - N) >= g$lo & (idx - 1L - N) <= g$hi
    sizes[comp] > 1L | self | wrap | back
}

# The cell width is tied to the jump size, so that the discretisation costs the same fraction of
# eps at every eps and the bracket of the measure test stays equally tight.
cells_for <- function(eps, per_eps = 40) as.integer(round(TWO_PI * per_eps / eps))

chain_recurrent <- function(sys, N, eps) {
    g <- transition_arcs(sys, N, eps)
    comp <- strong_components(g)
    rec <- recurrent_cells(g, comp)
    list(g = g, comp = comp, rec = rec, measure = sum(rec) * g$h)
}

closed_form_measure <- function(sys, eps) {
    if (sys$form == "circle") return(TWO_PI)
    if (eps >= sys$alpha) return(TWO_PI)
    4 * asin(eps / sys$alpha)
}

# ============================ #
# Test 1: the measure against the closed form ####
# ============================ #
#
# Lower bound. If theta is eps-chain recurrent then the cells of one of its eps-pseudo-orbits
# form a closed walk of the graph, so every cell meeting the true set is counted and
# m(R_graph) >= m(R_eps) = 4 asin(eps/alpha).
#
# Upper bound. An edge j -> l gives x in cell j and y in cell l with d(T(x), y) <= eps, so the
# centres satisfy d(T(c_j), c_l) <= eps + (L+1)*h/2 with L the Lipschitz constant of T. A graph
# cycle therefore gives an eps_plus-pseudo-orbit through the centres, every counted cell has its
# centre in R_{eps_plus}, and the cells reach at most h/2 beyond it on each side of each of the
# 2k arcs, so m(R_graph) <= 4 asin(eps_plus/alpha) + 2*k*h.

test_measure_against_closed_form <- function(eps_grid) {
    rows <- list()
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        for (eps in eps_grid) {
            N <- cells_for(eps)
            r <- chain_recurrent(sys, N, eps)
            h <- r$g$h
            eps_plus <- eps + (sys$lip(sys$alpha) + 1) * h / 2
            if (sys$form == "circle") {
                lower <- TWO_PI; upper <- TWO_PI
            } else {
                lower <- closed_form_measure(sys, eps)
                upper <- closed_form_measure(sys, eps_plus) + 2 * sys$k * h
            }
            rows[[length(rows) + 1L]] <- data.frame(
                system = nm, label = sys$label, alpha = sys$alpha, eps = eps, cells = N,
                measured = r$measure, exact = closed_form_measure(sys, eps),
                lower = lower, upper = upper,
                inside = r$measure >= lower - 1e-12 & r$measure <= upper + 1e-12,
                n_classes = sum(tabulate(r$comp[r$rec], nbins = max(r$comp)) > 0L),
                n_large_classes = sum(tabulate(r$comp[r$rec], nbins = max(r$comp)) > 1L),
                stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 2: a constructive eps-pseudo-orbit for a recurrent cell ####
# ============================ #
#
# The graph is consulted only for the starting point. The witness is built from the map itself,
# in the lift, where T is increasing. One step carries a set S to T(S) widened by eps, so the set
# reachable from x0 in n steps is the interval [l_n, r_n] with l_0 = r_0 = x0, l_{n+1} = T(l_n) -
# eps and r_{n+1} = T(r_n) + eps, and x0 is eps-chain recurrent as soon as one of its lifts falls
# inside that interval for some n >= 1. The path is then read off backwards: at each stage the
# admissible predecessors of y_k are T^{-1}([y_k - eps, y_k + eps]) intersected with the interval
# reached at that stage, which is nonempty by construction, and the midpoint is taken. Every step
# of the finished sequence is measured against eps afterwards, so the assertion tests the path
# and not the reasoning that produced it.

invert_lift <- function(sys, y, lo, hi) {
    Tl <- lift_map(sys)
    for (i in seq_len(80L)) {                 # bisection on an increasing lift
        mid <- (lo + hi) / 2
        if (Tl(mid) < y) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
}

pseudo_orbit_back_to_start <- function(sys, x0, eps, max_steps = 4000L) {
    Tl <- lift_map(sys)
    jump <- eps * (1 - 1e-9)                  # strictly admissible deviation
    l <- x0; r <- x0
    ls <- x0; rs <- x0
    hit <- NA_integer_; target <- NA_real_
    for (n in seq_len(max_steps)) {
        l <- Tl(l) - jump; r <- Tl(r) + jump
        ls <- c(ls, l); rs <- c(rs, r)
        m <- round((l - x0) / TWO_PI):round((r - x0) / TWO_PI)
        cand <- x0 + TWO_PI * m
        inside <- cand[cand >= l & cand <= r]
        if (length(inside) > 0L) { hit <- n; target <- inside[which.max(abs(inside - x0))]; break }
        if (r - l >= TWO_PI) { hit <- n; target <- x0 + TWO_PI * round((l - x0) / TWO_PI + 0.5); break }
    }
    if (is.na(hit)) return(NULL)
    y <- numeric(hit + 1L); y[hit + 1L] <- target
    for (k in seq(hit, 1L)) {
        lo_adm <- invert_lift(sys, y[k + 1L] - jump, ls[k] - TWO_PI, rs[k] + TWO_PI)
        hi_adm <- invert_lift(sys, y[k + 1L] + jump, ls[k] - TWO_PI, rs[k] + TWO_PI)
        a <- max(lo_adm, ls[k]); b <- min(hi_adm, rs[k])
        y[k] <- (a + b) / 2
    }
    y
}

pseudo_orbit_max_step <- function(sys, path) {
    if (is.null(path) || length(path) < 2L) return(Inf)
    Tl <- lift_map(sys)
    n <- length(path)
    max(abs(path[-1] - Tl(path[-n])))
}

test_pseudo_orbit_witnesses <- function(eps) {
    rows <- list()
    N <- cells_for(eps)
    for (nm in names(systems)) {
        sys <- systems[[nm]]
        r <- chain_recurrent(sys, N, eps)
        h <- r$g$h
        # The graph is an outer approximation, so its recurrent arcs reach a little beyond the
        # true set at their two ends. The witnesses are taken at the quarter, half and three
        # quarter points of each maximal arc of recurrent cells, which needs no tuned margin.
        runs <- if (all(r$rec)) list(c(1L, N)) else {
            idx <- which(r$rec)
            brk <- which(diff(idx) > 1L)
            Map(c, c(idx[1], idx[brk + 1L]), c(idx[brk], idx[length(idx)]))
        }
        fracs <- if (all(r$rec)) seq(0, 5) / 6 else c(0.25, 0.5, 0.75)
        for (rn in runs) {
            for (f in fracs) {
                x0 <- ((rn[1] + f * (rn[2] - rn[1]) - 1) + 0.5) * h
                path <- pseudo_orbit_back_to_start(sys, x0, eps)
                ok <- !is.null(path) &&
                      abs(((path[length(path)] - x0 + pi) %% TWO_PI) - pi) < 1e-9 &&
                      abs(path[1] - x0) < 1e-9
                rows[[length(rows) + 1L]] <- data.frame(
                    system = nm, x0 = x0, steps = if (is.null(path)) NA_integer_ else length(path) - 1L,
                    max_step = pseudo_orbit_max_step(sys, path), eps = eps, closed = ok,
                    stringsAsFactors = FALSE)
            }
        }
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 3: a barrier for a transient cell ####
# ============================ #
#
# A point can return to itself only if it can travel against the displacement somewhere. On an
# arc where |T(theta) - theta| > eps the displacement keeps one sign and beats every admissible
# jump, so each step advances in the direction of that sign and no eps-pseudo-orbit crosses the
# arc the other way. A transient point therefore carries a certificate with three parts: it lies
# in such an arc, the sign of the displacement is constant there, and a second arc of the
# opposite sign lies elsewhere on the circle, which closes the route round the far side. The
# certificate is read off the displacement and never consults the graph.

barrier_components <- function(sys, eps, n_grid = 60001L) {
    th <- seq(0, TWO_PI, length.out = n_grid)[-n_grid]
    d <- sys$disp(th, sys$alpha)
    blocked <- abs(d) > eps
    # label the cyclic runs of blocked grid points
    lab <- integer(length(th))
    if (!any(blocked)) return(list(th = th, d = d, lab = lab, n = 0L))
    start <- which(blocked & !blocked[c(length(th), seq_len(length(th) - 1L))])
    if (length(start) == 0L) { lab[blocked] <- 1L; return(list(th = th, d = d, lab = lab, n = 1L)) }
    for (i in seq_along(start)) {
        j <- start[i]
        repeat {
            lab[j] <- i
            j <- j %% length(th) + 1L
            if (!blocked[j] || lab[j] != 0L) break
        }
    }
    list(th = th, d = d, lab = lab, n = length(start))
}

test_barrier_witnesses <- function(eps) {
    rows <- list()
    for (nm in c("b", "c")) {
        sys <- systems[[nm]]
        N <- cells_for(eps)
        r <- chain_recurrent(sys, N, eps)
        bc <- barrier_components(sys, eps)
        signs <- vapply(seq_len(bc$n), function(i) sign(mean(bc$d[bc$lab == i])), numeric(1))
        pure <- vapply(seq_len(bc$n), function(i) all(sign(bc$d[bc$lab == i]) == signs[i]), logical(1))
        h <- r$g$h
        centres <- ((which(!r$rec) - 1) + 0.5) * h
        pick <- centres[round(seq(1, length(centres), length.out = min(12L, length(centres))))]
        for (x0 in pick) {
            i <- bc$lab[which.min(abs(bc$th - (x0 %% TWO_PI)))]
            rows[[length(rows) + 1L]] <- data.frame(
                system = nm, x0 = x0, in_barrier = i > 0L,
                sign_constant = if (i > 0L) pure[i] else FALSE,
                margin = abs(sys$disp(x0, sys$alpha)) - eps,
                opposite_exists = if (i > 0L) any(signs == -signs[i]) else FALSE,
                stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 4: the omega limit set of a set is not the union of the omega limit sets ####
# ============================ #
#
# For k = 1 the points -pi/2 and pi/2 are fixed and T is an increasing homeomorphism of the lift,
# so T maps the closed half-circle Y = [-pi/2, pi/2] onto itself, every forward image of Y is Y
# and omega(Y) = Y, of measure pi. Every interior point converges to pi/2, so the union of the
# omega limit sets of the points of Y is the two-point set and has measure zero. The convergence
# carries its own reference. Writing v = pi/2 - theta turns the two maps into v -> v - alpha
# sin(v) and v -> v - alpha sin^2(v), which the script verifies against the maps themselves and
# then iterates in v, where no cancellation occurs. The first has multiplier 1 - alpha at v = 0
# and its ratio v_{n+1}/v_n tends to that value; the second has multiplier 1, the fixed point is
# degenerate, and v_n ~ 1/(alpha n), so the product n * v_n falls to 1/alpha.

test_omega_of_a_set <- function(n_points = 2001L) {
    x0 <- seq(-pi / 2, pi / 2, length.out = n_points)
    rows <- list()
    for (nm in c("a", "b")) {
        sys <- systems[[nm]]
        Tl <- lift_map(sys)
        fixed_err <- max(abs(Tl(c(-pi / 2, pi / 2)) - c(-pi / 2, pi / 2)))
        ends <- c(-pi / 2, pi / 2)
        for (i in seq_len(500L)) ends <- Tl(ends)      # omega_n(Y) has the same endpoints for all n
        x <- x0
        horizon <- if (nm == "a") 16000L else 100L
        dist_at <- numeric(0)
        for (i in seq_len(horizon)) {
            x <- Tl(x)
            if (i %in% c(horizon %/% 4L, horizon %/% 2L, horizon)) dist_at <- c(dist_at, max(abs(x[-c(1, n_points)] - pi / 2)))
        }
        # the recursion conjugated to the distance v = pi/2 - theta, checked against the map at a
        # distance where the subtraction is still exact, then iterated in v alone
        rec <- if (nm == "b") function(v) v - sys$alpha * sin(v) else function(v) v - sys$alpha * sin(v)^2
        vv <- c(1, 0.5, 0.1, 0.01)
        conj_err <- max(abs(rec(vv) - (pi / 2 - Tl(pi / 2 - vv))))
        v <- 0.1; ratio_last <- NA_real_; prod_last <- NA_real_
        steps <- if (nm == "b") 40L else 200000L
        for (i in seq_len(steps)) { vn <- rec(v); if (i == steps) ratio_last <- vn / v; v <- vn }
        prod_last <- steps * v
        rows[[length(rows) + 1L]] <- data.frame(
            system = nm, fixed_err = fixed_err,
            endpoint_err = max(abs(ends - c(-pi / 2, pi / 2))),
            image_span = max(x0) - min(x0), omega_of_Y_measure = pi,
            horizon = horizon, final_distance = dist_at[3], conj_err = conj_err,
            rate_last = ratio_last, rate_target = if (nm == "b") 1 - sys$alpha else 1,
            n_times_distance = if (nm == "a") prod_last else NA_real_,
            algebraic_target = if (nm == "a") 1 / sys$alpha else NA_real_,
            stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 5: the degenerate cases ####
# ============================ #

test_degenerate <- function() {
    sys_b <- systems$b
    eps_big <- sys_b$alpha * 1.0001
    big <- chain_recurrent(sys_b, cells_for(eps_big), eps_big)  # eps above the largest displacement
    idm <- systems$b; idm$alpha <- 0
    ident <- chain_recurrent(idm, cells_for(1e-3), 1e-3)
    data.frame(eps_above_alpha_measure = big$measure,
               identity_measure = ident$measure,
               target = TWO_PI)
}

# ============================ #
# Run ####
# ============================ #

eps_grid <- c(0.2, 0.1, 0.05, 0.025)

t_measure <- test_measure_against_closed_form(eps_grid)
t_orbit   <- test_pseudo_orbit_witnesses(0.05)
t_barrier <- test_barrier_witnesses(0.05)
t_omega   <- test_omega_of_a_set()
t_degen   <- test_degenerate()

ok_measure <- all(t_measure$inside)
rel_err <- with(t_measure[t_measure$system != "a", ], max(abs(measured - exact) / exact))
ok_orbit <- all(t_orbit$closed) && all(t_orbit$max_step < 0.05) && nrow(t_orbit) >= 12L
ok_barrier <- all(t_barrier$in_barrier) && all(t_barrier$sign_constant) &&
              all(t_barrier$opposite_exists) && all(t_barrier$margin > 0)
# omega(Y) = Y: the endpoints are fixed to machine precision and the map is an increasing
# homeomorphism of the lift, so every forward image of Y is Y. The union of the pointwise omega
# limit sets is the two endpoints: the interior collapses to pi/2 at the two derived rates.
ok_omega <- all(t_omega$fixed_err < 1e-15) && all(t_omega$endpoint_err < 1e-12) &&
            all(t_omega$final_distance < 1e-3) && all(t_omega$conj_err < 1e-15) &&
            abs(t_omega$rate_last[t_omega$system == "b"] - 0.6) < 1e-6 &&
            # the remaining correction to n * v_n is of order log(n)/n, so one percent of the
            # limit is the band at n = 2 * 10^5
            abs(t_omega$n_times_distance[t_omega$system == "a"] -
                t_omega$algebraic_target[t_omega$system == "a"]) <
                0.01 * t_omega$algebraic_target[t_omega$system == "a"]
ok_degen <- abs(t_degen$eps_above_alpha_measure - TWO_PI) < 1e-12 &&
            abs(t_degen$identity_measure - TWO_PI) < 1e-12

status <- if (ok_measure && ok_orbit && ok_barrier && ok_omega && ok_degen) "pass" else "fail"

# ============================ #
# Figure ####
# ============================ #
#
# Two panels, both drawn from objects the tests have already produced. The upper panel shows the
# mechanism: the displacement against the angle, the band of width 2 eps that a jump can cover,
# and the cells the graph puts on a cycle. The lower panel shows the measure of that set against
# eps, beside the closed form.

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    labels <- vapply(systems, `[[`, character(1), "label")
    eps_show <- 0.05

    disp_rows <- do.call(rbind, lapply(names(systems), function(nm) {
        sys <- systems[[nm]]
        th <- seq(0, TWO_PI, length.out = 1200)
        data.frame(label = labels[[nm]], theta = th, disp = sys$disp(th, sys$alpha),
                   stringsAsFactors = FALSE)
    }))
    band_rows <- do.call(rbind, lapply(names(systems), function(nm) {
        sys <- systems[[nm]]
        N <- cells_for(eps_show)
        r <- chain_recurrent(sys, N, eps_show)
        h <- r$g$h
        idx <- which(r$rec)
        brk <- which(diff(idx) > 1L)                      # merge the consecutive cells into arcs
        starts <- c(idx[1], idx[brk + 1L]); ends <- c(idx[brk], idx[length(idx)])
        data.frame(label = labels[[nm]], xmin = (starts - 1) * h, xmax = ends * h,
                   stringsAsFactors = FALSE)
    }))
    disp_rows$label <- factor(disp_rows$label, levels = labels)
    band_rows$label <- factor(band_rows$label, levels = labels)

    p_top <- ggplot() +
        geom_rect(data = band_rows,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = "#056796", alpha = 0.22) +
        geom_ribbon(data = disp_rows, aes(x = theta, ymin = -eps_show, ymax = eps_show),
                    fill = "grey55", alpha = 0.3) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
        geom_line(data = disp_rows, aes(theta, disp), colour = "#7a2a2a", linewidth = 0.55) +
        facet_wrap(~ label, nrow = 1, labeller = kb_labeller()) +
        scale_x_continuous(breaks = c(0, pi, TWO_PI), labels = kb_ticks(c("0", "$\\pi$", "$2\\pi$"))) +
        labs(subtitle = kb_unicode("Displacement, the band of width $2\\varepsilon$ at $\\varepsilon = 0.05$, and the cells on a cycle of the graph"),
             x = NULL, y = kb_tex("$T(\\theta) - \\theta$")) +
        kb_theme()

    fine <- seq(0.005, 0.38, length.out = 300)
    curve_rows <- do.call(rbind, lapply(names(systems), function(nm) {
        sys <- systems[[nm]]
        data.frame(label = labels[[nm]], eps = fine,
                   measure = vapply(fine, function(e) closed_form_measure(sys, e), numeric(1)),
                   stringsAsFactors = FALSE)
    }))
    curve_rows$label <- factor(curve_rows$label, levels = labels)
    pts <- t_measure
    pts$label <- factor(pts$label, levels = labels)

    p_bot <- ggplot() +
        geom_line(data = curve_rows, aes(eps, measure, colour = label), linewidth = 0.6) +
        geom_point(data = pts, aes(eps, measured, colour = label),
                   size = 2, shape = 21, fill = "white", stroke = 0.7) +
        scale_colour_manual(values = c("#056796", "#be1117", "#379d08"), name = NULL,
                            labels = kb_tex(as.character(labels))) +
        scale_y_continuous(limits = c(0, TWO_PI + 0.25),
                           breaks = c(0, pi / 2, pi, 3 * pi / 2, TWO_PI),
                           labels = kb_ticks(c("0", "$\\pi/2$", "$\\pi$", "$3\\pi/2$", "$2\\pi$"))) +
        labs(subtitle = kb_unicode("Measure of the chain recurrent set against the jump size $\\varepsilon$"),
             x = kb_tex("Jump size $\\varepsilon$"), y = kb_tex("$m(\\mathcal{R}_\\varepsilon)$")) +
        kb_theme() +
        guides(colour = guide_legend(nrow = 1))

    p <- (p_top / p_bot) + plot_layout(heights = c(1, 1.15)) +
        plot_annotation(
            title = kb_unicode("Chain recurrent set of three circle maps"),
            caption = kb_caption(paste(
                "Upper: the displacement $T(\\theta) - \\theta$ in red, the grey band where its modulus is",
                "below $\\varepsilon = 0.05$, and in blue the cells that lie on a cycle of the",
                "$\\varepsilon$-transition graph. The recurrent cells sit where the displacement enters the",
                "band, except for the first map, whose displacement never changes sign and whose whole circle",
                "returns by going once round. Lower: lines give the closed form,",
                "$4\\arcsin(\\varepsilon/\\alpha)$ where the displacement changes sign and $2\\pi$ where it",
                "does not; circles give the measure counted from the graph, with the cell width held at",
                "$\\varepsilon/40$. Drawn by checks/chain-recurrence-of-circle-maps.R with $\\alpha = 0.4$",
                "for the first two maps and $\\alpha = 0.25$ for the third.")),
            theme = kb_theme())
    kb_save(p, "chain-recurrence-of-circle-maps", width = 7.8, height = 6.2)
}

emit("chain-recurrence-of-circle-maps", status,
     "The chain recurrent set of theta + alpha cos(k theta) has measure 4 asin(eps/alpha), and that of theta + alpha cos^2(theta) is the whole circle",
     list(max_relative_error_measure = rel_err,
          n_measure_cases = nrow(t_measure),
          large_classes_cos_map = t_measure$n_large_classes[t_measure$system == "b"][1],
          large_classes_cos3_map = t_measure$n_large_classes[t_measure$system == "c"][1],
          single_cell_classes_cos_map = with(t_measure[t_measure$system == "b", ], max(n_classes - n_large_classes)),
          all_within_discretisation_bracket = as.numeric(ok_measure),
          n_pseudo_orbit_witnesses = nrow(t_orbit),
          eps_minus_largest_step = 0.05 - max(t_orbit$max_step),
          eps_used_for_witnesses = 0.05,
          min_barrier_margin = min(t_barrier$margin),
          n_barrier_witnesses = nrow(t_barrier),
          omega_of_half_circle_measure = t_omega$omega_of_Y_measure[1],
          omega_endpoint_error = max(t_omega$endpoint_err),
          omega_interior_distance_to_pi_over_2 = max(t_omega$final_distance),
          geometric_rate_cos_map = t_omega$rate_last[t_omega$system == "b"],
          conjugated_recursion_error = max(t_omega$conj_err),
          algebraic_product_cos_squared_map = t_omega$n_times_distance[t_omega$system == "a"],
          identity_map_measure = t_degen$identity_measure,
          eps_above_alpha_measure = t_degen$eps_above_alpha_measure))
