# ============================ #
# Complete Lyapunov function built from attractor-repeller pairs ####
# ============================ #
#
# Claim. For a continuous map of a compact metric space there is a continuous V with three
# properties: V strictly decreases along orbits outside the chain recurrent set, the image of the
# chain recurrent set under V is compact and nowhere dense, and two chain recurrent points share
# a value of V exactly when they belong to the same chain class (norton1995, sec. 6, Thm. 4).
# Norton assembles V from the attractor-repeller pairs as the sum of 2 g_n / 3^n, where g_n is
# zero on the n-th attractor, one on its complementary repeller and strictly decreasing between
# them (norton1995, sec. 6, Lemma). The weights are a base-three expansion with digits 0 and 2,
# so the values taken on the chain recurrent set lie in the middle-thirds Cantor set, which is
# where the second property comes from. Conley makes the same observation for his own weights
# 3^-n, whose digits are 0 and 1 (conley1978a, sec. II.6.4B); with those the values lie in a
# Cantor set that is not the middle-thirds one, and the test below is written for Norton's
# weights, so it asks for no digit 1.
#
# Method. The same three circle maps as in chain-recurrence-of-circle-maps are cut into cells and
# the eps-transition graph is built, but nothing else is shared: the graph is held as an explicit
# edge list and its strong components come from Kosaraju's two-pass scan rather than from the
# Tarjan scan of that script, so the two agree on the chain recurrent set through independent
# code. The attractors used are the forward closures of the recurrent components, which are the
# attractors of the graph in the sense of the definition, and V is assembled from them by the
# formula above.
#
# References that the construction does not share.
#   1. The measure of the chain recurrent set against the closed form 4 asin(eps/alpha), derived
#      in chain-recurrence-of-circle-maps and bracketed here for the discretisation in the same
#      way, which tests the Kosaraju components against analysis.
#   2. The three properties of a complete Lyapunov function, each tested directly on the finished
#      V and on every edge of the graph, and not on the recipe that produced it.
#   3. The membership of every recorded value of V on the chain recurrent set in the Cantor set,
#      tested by expanding the value in base three and requiring no digit one.
#   4. The intersection of the attractor-repeller pairs against the chain recurrent set, which is
#      the statement of Theorem 3 of the source.
#
# Degenerate case: the first map, whose chain recurrent set is the whole circle, so that the
# gradient-like part is empty and V is constant.

set.seed(20260919L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

TWO_PI <- 2 * pi

systems <- list(
    a = list(label = "$T(\\theta) = \\theta + \\alpha\\cos^2\\theta$", alpha = 0.4,
             disp = function(th, alpha) alpha * cos(th)^2, lip = function(alpha) 1 + alpha,
             k = NA_integer_),
    b = list(label = "$T(\\theta) = \\theta + \\alpha\\cos\\theta$", alpha = 0.4,
             disp = function(th, alpha) alpha * cos(th), lip = function(alpha) 1 + alpha, k = 1L),
    c = list(label = "$T(\\theta) = \\theta + \\alpha\\cos 3\\theta$", alpha = 0.25,
             disp = function(th, alpha) alpha * cos(3 * th), lip = function(alpha) 1 + 3 * alpha,
             k = 3L)
)

lift_map <- function(sys) function(th) th + sys$disp(th, sys$alpha)

# ============================ #
# The graph as an explicit edge list ####
# ============================ #
#
# Cell j is [ (j-1) h, j h ] for j = 1, ..., N. An edge j -> l is recorded when some x in cell j
# and some y in cell l satisfy d(T(x), y) <= eps. Both directions are stored in compressed form,
# the out-edges sorted by tail and the in-edges by head, which is what the two passes of Kosaraju
# need.

build_graph <- function(sys, N, eps) {
    h <- TWO_PI / N
    a <- (0:(N - 1)) * h
    Tl <- lift_map(sys)
    lo <- floor((Tl(a) - eps) / h)
    hi <- ceiling((Tl(a + h) + eps) / h) - 1
    len <- hi - lo + 1L
    from <- rep(seq_len(N), len)
    to <- as.integer((rep(lo, len) + sequence(len) - 1L) %% N) + 1L
    ord_out <- order(from); ord_in <- order(to)
    list(N = N, h = h, eps = eps, lo = lo, hi = hi, len = len,
         from = from, to = to,
         out_head = c(0L, cumsum(tabulate(from, nbins = N))), out_adj = to[ord_out],
         in_head  = c(0L, cumsum(tabulate(to,   nbins = N))), in_adj  = from[ord_in])
}

# Kosaraju. The first pass records the order in which the vertices finish on the forward graph;
# the second walks the reverse graph in the opposite order, and each tree it grows is a strong
# component. Components come out in topological order, so a later component never reaches an
# earlier one, and reversing that order puts the sinks first.
kosaraju <- function(g) {
    N <- g$N
    seen <- logical(N); finish <- integer(N); nf <- 0L
    stackv <- integer(N + 1L); stackk <- integer(N + 1L)
    for (s in seq_len(N)) {
        if (seen[s]) next
        sp <- 1L; stackv[1L] <- s; stackk[1L] <- 0L; seen[s] <- TRUE
        while (sp >= 1L) {
            v <- stackv[sp]; k <- stackk[sp] + 1L
            deg <- g$out_head[v + 1L] - g$out_head[v]
            if (k <= deg) {
                stackk[sp] <- k
                w <- g$out_adj[g$out_head[v] + k]
                if (!seen[w]) { seen[w] <- TRUE; sp <- sp + 1L; stackv[sp] <- w; stackk[sp] <- 0L }
            } else { nf <- nf + 1L; finish[nf] <- v; sp <- sp - 1L }
        }
    }
    comp <- integer(N); nc <- 0L
    for (i in rev(seq_len(N))) {
        s <- finish[i]
        if (comp[s] != 0L) next
        nc <- nc + 1L
        sp <- 1L; stackv[1L] <- s; comp[s] <- nc
        while (sp >= 1L) {
            v <- stackv[sp]; sp <- sp - 1L
            for (k in seq_len(g$in_head[v + 1L] - g$in_head[v])) {
                w <- g$in_adj[g$in_head[v] + k]
                if (comp[w] == 0L) { comp[w] <- nc; sp <- sp + 1L; stackv[sp] <- w }
            }
        }
    }
    list(comp = comp, n = nc)                 # component 1 is a source, component n a sink
}

recurrent_cells <- function(g, comp) {
    sizes <- tabulate(comp, nbins = max(comp))
    idx <- seq_len(g$N)
    self <- (idx - 1L) >= g$lo & (idx - 1L) <= g$hi
    wrap <- (idx - 1L + g$N) >= g$lo & (idx - 1L + g$N) <= g$hi
    back <- (idx - 1L - g$N) >= g$lo & (idx - 1L - g$N) <= g$hi
    sizes[comp] > 1L | self | wrap | back
}

# Forward reachability from a set of cells, grown arc by arc.
forward_reach <- function(g, seed) {
    N <- g$N
    seen <- logical(N); seen[seed] <- TRUE
    frontier <- seed
    while (length(frontier) > 0L) {
        len <- g$len[frontier]
        cand <- as.integer((rep(g$lo[frontier], len) + sequence(len) - 1L) %% N) + 1L
        cand <- unique(cand[!seen[cand]])
        seen[cand] <- TRUE
        frontier <- cand
    }
    seen
}

# ============================ #
# The attractor-repeller pairs and the function V ####
# ============================ #
#
# A_n is the forward closure of the n-th recurrent component, which is forward invariant and a
# union of components, so it is an attractor of the graph. Its basin D(A_n) holds the cells whose
# reachable recurrent components all lie inside A_n, and the complementary repeller is everything
# else. The mask of a cell records which recurrent components it reaches and is propagated in one
# pass over the cells taken sinks first, so that a cell is treated after everything it reaches.
# On D(A_n) minus A_n the function g_n is built from the longest remaining path, which decreases
# strictly along every edge of that region and lands at zero on entering A_n.

lyapunov <- function(g, comp, rec) {
    N <- g$N
    rec_comps <- sort(unique(comp[rec]))
    stopifnot(length(rec_comps) <= 30L)           # the mask below holds one bit per component
    nc <- max(comp)
    cells_by_comp <- order(comp, decreasing = TRUE)   # sinks first: Kosaraju numbers sources low
    bit <- integer(nc); bit[rec_comps] <- bitwShiftL(1L, seq_along(rec_comps) - 1L)

    # Which recurrent components each component reaches. A component is settled only after every
    # component it reaches, so the mask is accumulated per component and never per cell, and only
    # the neighbours in strictly later components are read.
    mask_comp <- bit
    for (i in cells_by_comp) {
        ci <- comp[i]
        nb <- as.integer((g$lo[i] + seq_len(g$len[i]) - 1L) %% N) + 1L
        nb <- nb[comp[nb] > ci]
        if (length(nb) > 0L) mask_comp[ci] <- Reduce(bitwOr, mask_comp[comp[nb]], mask_comp[ci])
    }
    mask <- mask_comp[comp]

    # One pair (A_n, A_n*) per recurrent component, with g_n zero on the attractor, one on the
    # complementary repeller and strictly decreasing between them.
    gl <- vector("list", length(rec_comps)); span <- integer(length(rec_comps))
    for (n in seq_along(rec_comps)) {
        seed <- which(comp == rec_comps[n])
        inA <- forward_reach(g, seed)
        mA <- mask[seed[1]]                        # the recurrent components inside A_n
        inD <- bitwAnd(mask, bitwNot(mA)) == 0L
        mid <- inD & !inA                          # the basin outside the attractor
        stopifnot(!any(mid & rec))                 # no chain class survives between A_n and A_n*
        hgt <- integer(N)
        if (any(mid)) {
            for (i in cells_by_comp) {
                if (!mid[i]) next
                nb <- as.integer((g$lo[i] + seq_len(g$len[i]) - 1L) %% N) + 1L
                nb <- nb[mid[nb] & comp[nb] > comp[i]]
                hgt[i] <- if (length(nb) > 0L) 1L + max(hgt[nb]) else 0L
            }
        }
        H <- if (any(mid)) max(hgt[mid]) else 0L
        gn <- numeric(N)
        gn[!inD] <- 1
        gn[mid] <- (hgt[mid] + 1) / (H + 2)
        gl[[n]] <- gn; span[n] <- sum(mid)
    }

    # The order of the pairs in the sum is free, and it decides which pair takes the leading
    # digit; Conley states that the ordering of the components depends on the numbering of the
    # pairs (conley1978a, sec. II.6.4B). Taking them by decreasing size of the region they
    # separate puts the descent across the long transient arcs in the first digits, where a plot
    # can show it, and leaves the pairs that separate a single cell in the tail. Every order
    # gives a complete Lyapunov function; this one gives a legible one.
    ord <- order(span, decreasing = TRUE)
    gmat <- do.call(cbind, gl[ord])
    V <- as.numeric(gmat %*% (2 / 3^seq_along(ord)))
    list(V = V, g = gmat, rec_comps = rec_comps[ord], span = span[ord])
}

# ============================ #
# The three properties, tested on the finished V ####
# ============================ #

clause_strict_decrease <- function(g, V, rec) {
    keep <- !rec[g$from]
    if (!any(keep)) return(list(ok = TRUE, worst = NA_real_, n_edges = 0L))
    slack <- V[g$from[keep]] - V[g$to[keep]]
    list(ok = all(slack > 0), worst = min(slack), n_edges = sum(keep))
}

# V must be constant on each chain class and must take a different value on each. The sum that
# builds V is reassociated by the matrix product, so cells of one class can differ in the last
# bit; the test asks that the spread inside a class be a rounding error and that the gap between
# classes exceed it by three orders of magnitude, which is the statement that survives in
# floating point.
clause_classes <- function(comp, V, rec) {
    cl <- split(V[rec], comp[rec])
    spread <- max(vapply(cl, function(v) diff(range(v)), numeric(1)))
    vals <- vapply(cl, function(v) v[1], numeric(1))
    gap <- if (length(vals) > 1L) min(diff(sort(vals))) else Inf
    list(ok = spread < 1e-12 && length(unique(vals)) == length(vals) &&
              gap > 1e3 * max(spread, .Machine$double.eps),
         spread = spread, n_classes = length(vals), min_gap = gap)
}

# A number of [0, 1] lies in the middle-thirds Cantor set when its base-three expansion uses the
# digits 0 and 2 only. The values here are finite sums of 2/3^n, so 3^depth times a value is an
# integer, and the digits are read from that integer exactly rather than by repeated
# multiplication, which would let rounding turn a 2 into a 1. The residual of the rounding is
# reported beside the digits.
clause_cantor <- function(V, rec, depth) {
    x <- unique(V[rec])
    scaled <- x * 3^depth
    m <- round(scaled)
    worst_digit <- 0
    for (mi in m) {
        r <- mi
        for (i in seq_len(depth)) { d <- r %% 3; r <- (r - d) / 3; worst_digit <- max(worst_digit, min(abs(d), abs(d - 2))) }
        worst_digit <- max(worst_digit, abs(r))      # nothing may be left above the last digit
    }
    list(ok = worst_digit == 0 && max(abs(scaled - m)) < 1e-3,
         worst_digit_distance = worst_digit, rounding_residual = max(abs(scaled - m)))
}

# Theorem 3: the chain recurrent set is the intersection of the unions A + A* over the attractors.
clause_theorem_3 <- function(g, comp, rec, lyap) {
    N <- g$N
    inter <- rep(TRUE, N)
    rec_comps <- lyap$rec_comps
    for (n in seq_along(rec_comps)) {
        gn <- lyap$g[, n]
        inter <- inter & (gn == 0 | gn == 1)       # g_n is 0 on A_n and 1 on A_n*
    }
    list(ok = all(inter == rec), n_wrong = sum(inter != rec))
}

# ============================ #
# Run ####
# ============================ #

eps <- 0.05
per_eps <- 40
results <- list()
for (nm in names(systems)) {
    sys <- systems[[nm]]
    N <- as.integer(round(TWO_PI * per_eps / eps))
    g <- build_graph(sys, N, eps)
    kk <- kosaraju(g)
    rec <- recurrent_cells(g, kk$comp)
    lyap <- lyapunov(g, kk$comp, rec)
    exact <- if (is.na(sys$k)) TWO_PI else 4 * asin(eps / sys$alpha)
    eps_plus <- eps + (sys$lip(sys$alpha) + 1) * g$h / 2
    upper <- if (is.na(sys$k)) TWO_PI else 4 * asin(eps_plus / sys$alpha) + 2 * sys$k * g$h
    results[[nm]] <- list(
        sys = sys, g = g, comp = kk$comp, rec = rec, lyap = lyap,
        measure = sum(rec) * g$h, exact = exact, lower = exact, upper = upper,
        dec = clause_strict_decrease(g, lyap$V, rec),
        cls = clause_classes(kk$comp, lyap$V, rec),
        can = clause_cantor(lyap$V, rec, depth = length(lyap$rec_comps) + 2L),
        th3 = clause_theorem_3(g, kk$comp, rec, lyap))
}

ok_measure <- all(vapply(results, function(r) r$measure >= r$lower - 1e-12 && r$measure <= r$upper + 1e-12, logical(1)))
ok_dec     <- all(vapply(results, function(r) r$dec$ok, logical(1)))
ok_cls     <- all(vapply(results, function(r) r$cls$ok, logical(1)))
ok_can     <- all(vapply(results, function(r) r$can$ok, logical(1)))
ok_th3     <- all(vapply(results, function(r) r$th3$ok, logical(1)))
# The first map has no gradient-like part, so V is constant and no edge leaves the recurrent set.
ok_degen   <- results$a$dec$n_edges == 0L && diff(range(results$a$lyap$V)) == 0 &&
              results$a$cls$n_classes == 1L

status <- if (ok_measure && ok_dec && ok_cls && ok_can && ok_th3 && ok_degen) "pass" else "fail"

# ============================ #
# Figure ####
# ============================ #
#
# The theorem drawn: on the shaded arcs V is flat and the dynamics returns, and between them V
# descends and the dynamics runs one way. Both series come from the objects tested above.

if (requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages(library(ggplot2))
    labels <- vapply(systems, `[[`, character(1), "label")
    vrows <- do.call(rbind, lapply(names(systems), function(nm) {
        r <- results[[nm]]
        data.frame(label = labels[[nm]], theta = ((seq_len(r$g$N) - 1) + 0.5) * r$g$h,
                   V = r$lyap$V, stringsAsFactors = FALSE)
    }))
    brows <- do.call(rbind, lapply(names(systems), function(nm) {
        r <- results[[nm]]
        idx <- which(r$rec); brk <- which(diff(idx) > 1L)
        starts <- c(idx[1], idx[brk + 1L]); ends <- c(idx[brk], idx[length(idx)])
        data.frame(label = labels[[nm]], xmin = (starts - 1) * r$g$h, xmax = ends * r$g$h,
                   stringsAsFactors = FALSE)
    }))
    nrows <- do.call(rbind, lapply(names(systems), function(nm) {
        r <- results[[nm]]
        data.frame(label = labels[[nm]], theta = 0.12,
                   V = 1.04,
                   text = sprintf("%d chain %s", r$cls$n_classes,
                                  if (r$cls$n_classes == 1L) "class" else "classes"),
                   stringsAsFactors = FALSE)
    }))
    for (d in c("vrows", "brows", "nrows")) {
        x <- get(d); x$label <- factor(x$label, levels = labels); assign(d, x)
    }
    p <- ggplot() +
        geom_rect(data = brows, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = "#056796", alpha = 0.22) +
        geom_line(data = vrows, aes(theta, V), colour = "#7a2a2a", linewidth = 0.6) +
        geom_text(data = nrows, aes(theta, V, label = text), hjust = 0, size = 2.7,
                  colour = "grey25") +
        facet_wrap(~ label, nrow = 1, labeller = kb_labeller()) +
        scale_y_continuous(limits = c(-0.02, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        scale_x_continuous(breaks = c(0, pi, TWO_PI), labels = kb_ticks(c("0", "$\\pi$", "$2\\pi$"))) +
        labs(title = kb_unicode("The decomposition of a circle map drawn by its Lyapunov function"),
             subtitle = kb_unicode("Complete Lyapunov function $V$ assembled from the attractor-repeller pairs, with the chain recurrent arcs shaded"),
             x = kb_tex("Angle $\\theta$"), y = kb_tex("$V(\\theta)$"),
             caption = kb_caption(paste(
                 "On a shaded arc $V$ is constant and the points return to themselves along",
                 "$\\varepsilon$-pseudo-orbits; between two arcs $V$ falls at every step and the points run",
                 "one way. The first map has one chain class, the whole circle, so $V$ is constant and",
                 "nothing descends. Values of $V$ on the shaded arcs are sums of $2/3^n$ and so lie in the",
                 "middle-thirds Cantor set, which is where the image of the chain recurrent set gets its",
                 "nowhere density. The pairs are taken in decreasing order of the region each separates,",
                 "which is free to choose and puts the long descents in the leading digits. Each arc has a",
                 "few single-cell classes at its shoulders, which are the outer approximation of the grid;",
                 "they count as chain classes of the graph, and the step at the left edge of the first",
                 "shaded arc of the second panel is one of them. Their measure falls with the cell width.",
                 "Drawn by checks/complete-lyapunov-function-of-circle-maps.R at $\\varepsilon = 0.05$ with",
                 "the cell width $\\varepsilon/40$."))) +
        kb_theme()
    kb_save(p, "complete-lyapunov-function-of-circle-maps", width = 7.8, height = 4.2)
}

emit("complete-lyapunov-function-of-circle-maps", status,
     "A complete Lyapunov function built from the attractor-repeller pairs is flat on each chain class, strictly decreasing elsewhere, and takes values in the Cantor set on the chain recurrent set",
     list(chain_classes_cos_map = results$b$cls$n_classes,
          chain_classes_cos3_map = results$c$cls$n_classes,
          chain_classes_cos_squared_map = results$a$cls$n_classes,
          measure_cos_map = results$b$measure,
          measure_cos_map_closed_form = results$b$exact,
          measure_cos3_map = results$c$measure,
          measure_cos3_map_closed_form = results$c$exact,
          edges_out_of_the_recurrent_set = results$b$dec$n_edges + results$c$dec$n_edges,
          smallest_decrease_per_edge = min(results$b$dec$worst, results$c$dec$worst),
          spread_of_V_within_a_class = max(vapply(results, function(r) r$cls$spread, numeric(1))),
          ratio_of_the_class_gap_to_the_spread = min(results$b$cls$min_gap, results$c$cls$min_gap) /
              max(max(vapply(results, function(r) r$cls$spread, numeric(1))), .Machine$double.eps),
          smallest_gap_between_classes = min(results$b$cls$min_gap, results$c$cls$min_gap),
          worst_base_three_digit_distance = max(vapply(results, function(r) r$can$worst_digit_distance, numeric(1))),
          base_three_rounding_residual = max(vapply(results, function(r) r$can$rounding_residual, numeric(1))),
          cells_separated_by_the_leading_pair_cos_map = results$b$lyap$span[1],
          cells_separated_by_the_leading_pair_cos3_map = results$c$lyap$span[1],
          cells_where_theorem_3_disagrees = sum(vapply(results, function(r) r$th3$n_wrong, numeric(1))),
          constant_V_on_the_first_map = as.numeric(diff(range(results$a$lyap$V)))))
