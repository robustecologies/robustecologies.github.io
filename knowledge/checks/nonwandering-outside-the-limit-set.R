# ============================ #
# A nonwandering point that lies in no omega limit set ####
# ============================ #
#
# Claim. The inclusion of the limit set in the nonwandering set is strict. On the space
#
#     X = {z_m : m >= 1} + {0^j 1 0^inf : j >= 0} + {0^inf},      z_m = 1 0^m 1 0^inf,
#
# a closed shift-invariant subset of the one-sided sequences on two symbols, with the shift T and
# the metric d(x, y) = 2^-k where k is the first index at which x and y differ, the limit set is
# the single fixed point 0^inf while the nonwandering set is {0^inf, 1 0^inf}. The point q =
# 1 0^inf is therefore nonwandering and lies in no omega limit set, and it sits at distance 1
# from the limit set, the largest distance the metric allows. The construction is made here; the
# definitions are those of robinson1999, sec. 2.3, and conley1978a, sec. II.4.2.
#
# Why it works. Every element of X holds at most two symbols 1, so no orbit visits the cylinder
# of words beginning with 1 infinitely often, and therefore no omega limit set meets that
# cylinder; since the cylinder is open and closed, the limit set misses it as well, and q lies in
# it. On the other side, the cylinder [1 0^N] holds the generator z_N together with q, and
# T^(N+1) carries z_N to q, so every neighbourhood of q meets one of its own forward images and
# q is nonwandering.
#
# Method and references. Each element is a finite word with an implicit tail of zeros, so the
# tests are exact string computations and not approximations. Nothing is shared with the
# statement being tested: the sets are recomputed from the definitions, by searching the
# enumerated space for the witnesses that the definitions ask for.
#
#   1. Forward invariance and the count of ones, by enumeration.
#   2. The omega limit set of every element, by iterating to the fixed word.
#   3. The nonwandering property of q, by an exhaustive search for a returning point in every
#      cylinder, which recovers the witness the argument above names.
#   4. The wandering property of every other element, by exhibiting a cylinder that holds that
#      element alone and no return.
#   5. Chain recurrence of q, by an explicit epsilon-pseudo-orbit whose every step is measured
#      against epsilon in the metric.
#
# Degenerate case: the same computation with finitely many generators, where q turns out to be
# wandering, which shows that the infinite family is doing the work.

set.seed(20260920L)

emit <- function(id, status, summary, metrics) {
    val <- function(v) if (is.character(v)) sprintf("\"%s\"", v) else if (!is.finite(v)) "null" else formatC(v, digits = 6, format = "g")
    fields <- paste(sprintf("\"%s\": %s", names(metrics), vapply(metrics, val, character(1))), collapse = ", ")
    cat(sprintf("{\"id\": \"%s\", \"status\": \"%s\", \"summary\": \"%s\", \"metrics\": {%s}}\n", id, status, summary, fields))
}

# ============================ #
# The space ####
# ============================ #
#
# A point is a word of length W with an implicit tail of zeros. W is taken past the last symbol 1
# of every element enumerated, so the word determines the point.

zero_word <- function(W) integer(W)

generator <- function(m, W) {           # z_m = 1 0^m 1 0^inf
    w <- zero_word(W); w[1] <- 1L; w[m + 2L] <- 1L; w
}

single <- function(j, W) {              # 0^j 1 0^inf
    w <- zero_word(W); w[j + 1L] <- 1L; w
}

shift <- function(w) c(w[-1], 0L)

as_key <- function(w) paste(w, collapse = "")

# distance d(x, y) = 2^-k with k the first index at which the two differ, and zero when equal
distance <- function(x, y) {
    i <- which(x != y)
    if (length(i) == 0L) 0 else 2^(-(min(i) - 1L))
}

prefix_equal <- function(x, y, n) all(x[seq_len(n)] == y[seq_len(n)])

build_space <- function(M) {
    W <- M + 4L
    pts <- list()
    for (m in seq_len(M)) pts[[length(pts) + 1L]] <- generator(m, W)
    for (j in 0:(M + 1L)) pts[[length(pts) + 1L]] <- single(j, W)
    pts[[length(pts) + 1L]] <- zero_word(W)
    keys <- vapply(pts, as_key, character(1))
    keep <- !duplicated(keys)
    list(pts = pts[keep], keys = keys[keep], W = W, M = M,
         q = single(0L, W), zero = zero_word(W))
}

# ============================ #
# Test 1: forward invariance and the count of ones ####
# ============================ #
#
# The enumeration is truncated at M generators, and the shift of a generator leaves the family of
# generators for the family of single-symbol points, which the enumeration also holds. A word
# outside the enumeration would show here as a shift with no match.

test_invariance <- function(sp) {
    imgs <- lapply(sp$pts, shift)
    matched <- vapply(imgs, function(w) as_key(w) %in% sp$keys, logical(1))
    ones <- vapply(sp$pts, sum, numeric(1))
    list(n_points = length(sp$pts), all_images_inside = all(matched),
         max_ones = max(ones), n_with_two_ones = sum(ones == 2))
}

# ============================ #
# Test 2: the omega limit set of every element ####
# ============================ #
#
# Every element holds finitely many symbols 1, so the shift reaches the zero word in at most W
# steps and stays. The omega limit set of every element is therefore the single fixed point, and
# the limit set of the system, being the closure of their union, is that point as well.

test_omega <- function(sp) {
    steps <- vapply(sp$pts, function(w) {
        n <- 0L
        while (any(w != 0L) && n <= sp$W + 2L) { w <- shift(w); n <- n + 1L }
        if (any(w != 0L)) NA_integer_ else n
    }, integer(1))
    fixed_zero <- identical(shift(sp$zero), sp$zero)
    list(all_reach_zero = all(!is.na(steps)), max_steps = max(steps, na.rm = TRUE),
         zero_is_fixed = fixed_zero)
}

# ============================ #
# Test 3: q is nonwandering ####
# ============================ #
#
# For each cylinder [1 0^N] the search runs over every enumerated point of that cylinder and
# every positive time, and asks whether the image stays in the cylinder. No witness is supplied
# to the search; the one it finds is compared afterwards with the point the argument names.

test_nonwandering_q <- function(sp, N_max) {
    rows <- list()
    for (N in seq_len(N_max)) {
        len <- N + 1L                        # the word 1 0^N
        inside <- vapply(sp$pts, function(w) prefix_equal(w, sp$q, len), logical(1))
        found_u <- NA_integer_; found_n <- NA_integer_
        for (i in which(inside)) {
            w <- sp$pts[[i]]
            for (n in seq_len(sp$W)) {
                if (prefix_equal(Reduce(function(a, b) shift(a), seq_len(n), w), sp$q, len)) {
                    found_u <- i; found_n <- n; break
                }
            }
            if (!is.na(found_u)) break
        }
        predicted_u <- as_key(generator(N, sp$W))
        rows[[length(rows) + 1L]] <- data.frame(
            N = N, n_in_cylinder = sum(inside), found = !is.na(found_u),
            time = found_n,
            witness_is_the_predicted_generator =
                !is.na(found_u) && identical(sp$keys[found_u], predicted_u),
            predicted_time = N + 1L, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Test 4: every other element is wandering ####
# ============================ #
#
# A point is wandering when some cylinder about it holds no point that returns to that cylinder.
# The search reports the shortest such cylinder, and the two points of the nonwandering set are
# expected to have none.

shortest_wandering_cylinder <- function(sp, idx, max_len) {
    y <- sp$pts[[idx]]
    for (len in seq_len(max_len)) {
        inside <- which(vapply(sp$pts, function(w) prefix_equal(w, y, len), logical(1)))
        returns <- FALSE
        for (i in inside) {
            w <- sp$pts[[i]]
            for (n in seq_len(sp$W)) {
                w <- shift(w)
                if (prefix_equal(w, y, len)) { returns <- TRUE; break }
            }
            if (returns) break
        }
        if (!returns) return(len)
    }
    NA_integer_
}

# The truncation decides how far the search can see. A cylinder of length L about q holds the
# generator z_{L-1}, which the enumeration carries only while L - 1 <= M, so no cylinder longer
# than M + 1 can certify anything about q and the search is capped there. A generator z_m is
# isolated by a cylinder of length m + 3 and a point 0^j 1 0^inf by one of length j + 2, so the
# same cap classifies every generator with m <= M - 2 and every other point with j <= M - 1. The
# argument for the whole family, which no truncation reaches, is in the note.
test_wandering <- function(sp) {
    L <- sp$M + 1L
    lens <- vapply(seq_along(sp$pts), function(i) shortest_wandering_cylinder(sp, i, L), integer(1))
    is_q <- vapply(sp$pts, function(w) identical(w, sp$q), logical(1))
    is_zero <- vapply(sp$pts, function(w) identical(w, sp$zero), logical(1))
    resolved <- vapply(seq_along(sp$pts), function(i) {
        w <- sp$pts[[i]]
        ones <- which(w == 1L)
        if (length(ones) == 2L) (ones[2] - 1L) + 2L <= L else if (length(ones) == 1L) ones[1] + 1L <= L else TRUE
    }, logical(1))
    others <- !is_q & !is_zero & resolved
    list(cap = L,
         q_nonwandering = is.na(lens[is_q]),
         zero_nonwandering = is.na(lens[is_zero]),
         n_others_resolved = sum(others),
         others_wandering = all(!is.na(lens[others])),
         n_nonwandering_among_resolved = sum(is.na(lens[resolved | is_q | is_zero])),
         max_cylinder = max(lens[others], na.rm = TRUE))
}

# ============================ #
# Test 5: q is chain recurrent, by an explicit pseudo-orbit ####
# ============================ #
#
# At accuracy eps = 2^-N the sequence starts at q, jumps to 0^(N+1) 1 0^inf, whose distance to
# the image of q is 2^-(N+1), and then follows the true orbit back to q. Every step is measured
# against eps afterwards.

test_chain_recurrence_of_q <- function(sp, N_max) {
    rows <- list()
    for (N in seq_len(N_max)) {
        eps <- 2^(-N)
        path <- list(sp$q, single(N + 1L, sp$W))
        w <- path[[2]]
        for (i in seq_len(N + 1L)) { w <- shift(w); path[[length(path) + 1L]] <- w }
        steps <- vapply(seq_len(length(path) - 1L),
                        function(i) distance(shift(path[[i]]), path[[i + 1L]]), numeric(1))
        closed <- identical(path[[length(path)]], sp$q)
        rows[[length(rows) + 1L]] <- data.frame(
            N = N, eps = eps, length = length(path) - 1L, max_step = max(steps),
            closed = closed, admissible = max(steps) < eps, stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# ============================ #
# Run ####
# ============================ #

sp <- build_space(14L)
t_inv <- test_invariance(sp)
t_om <- test_omega(sp)
t_nw <- test_nonwandering_q(sp, N_max = 10L)
t_wa <- test_wandering(sp)
t_cr <- test_chain_recurrence_of_q(sp, N_max = 10L)
dist_q_to_limit <- distance(sp$q, sp$zero)

# The degenerate case: with finitely many generators the cylinders below the last one still hold
# a returning point and the ones above it do not, so q is wandering in the truncated space. The
# separation needs the whole family.
sp_small <- build_space(3L)
t_small <- test_nonwandering_q(sp_small, N_max = 6L)

ok_inv <- t_inv$all_images_inside && t_inv$max_ones == 2
ok_om <- t_om$all_reach_zero && t_om$zero_is_fixed
ok_nw <- all(t_nw$found) && all(t_nw$witness_is_the_predicted_generator) &&
         all(t_nw$time == t_nw$predicted_time)
ok_wa <- t_wa$q_nonwandering && t_wa$zero_nonwandering && t_wa$others_wandering &&
         t_wa$n_nonwandering_among_resolved == 2L
ok_cr <- all(t_cr$closed) && all(t_cr$admissible)
ok_sep <- dist_q_to_limit == 1
ok_small <- any(!t_small$found)

status <- if (ok_inv && ok_om && ok_nw && ok_wa && ok_cr && ok_sep && ok_small) "pass" else "fail"

# ============================ #
# Figure ####
# ============================ #
#
# The space drawn as words. Left: the generators and the points they pass through, one row each,
# a filled cell for the symbol 1. Right: the forward orbit of one generator, which walks the
# symbol 1 to the left until it leaves, passing through q on the way to the fixed point.

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    W_show <- 14L
    rows <- list(); labs_row <- character(0); kinds <- character(0)
    for (m in 1:6) {
        rows[[length(rows) + 1L]] <- generator(m, W_show)
        labs_row <- c(labs_row, sprintf("$z_{%d}$", m)); kinds <- c(kinds, "Wandering")
    }
    for (j in 0:4) {
        rows[[length(rows) + 1L]] <- single(j, W_show)
        labs_row <- c(labs_row, if (j == 0) "$q = 1 0^{\\infty}$" else sprintf("$0^{%d} 1 0^{\\infty}$", j))
        kinds <- c(kinds, if (j == 0) "Nonwandering" else "Wandering")
    }
    rows[[length(rows) + 1L]] <- zero_word(W_show)
    labs_row <- c(labs_row, "$0^{\\infty}$")
    kinds <- c(kinds, "Nonwandering and the limit set")
    grid_a <- do.call(rbind, lapply(seq_along(rows), function(i)
        data.frame(row = i, pos = seq_len(W_show) - 1L, v = rows[[i]],
                   kind = kinds[i], stringsAsFactors = FALSE)))
    shown_a <- kb_unicode(labs_row)
    grid_a$label <- factor(shown_a[grid_a$row], levels = rev(shown_a))

    orb <- list(); orb_lab <- character(0)
    w <- generator(6L, W_show)
    for (n in 0:8) {
        orb[[length(orb) + 1L]] <- w
        orb_lab <- c(orb_lab, sprintf("$T^{%d} z_6$", n))
        w <- shift(w)
    }
    grid_b <- do.call(rbind, lapply(seq_along(orb), function(i)
        data.frame(row = i, pos = seq_len(W_show) - 1L, v = orb[[i]], stringsAsFactors = FALSE)))
    shown_b <- kb_unicode(orb_lab)
    grid_b$label <- factor(shown_b[grid_b$row], levels = rev(shown_b))
    marks <- data.frame(label = factor(shown_b[c(8, 9)], levels = rev(shown_b)),
                        pos = W_show - 1.2, text = kb_unicode(c("$= q$", "$= 0^{\\infty}$")))

    pal <- c("Wandering" = "grey60", "Nonwandering" = "#056796",
             "Nonwandering and the limit set" = "#be1117")
    marks_a <- unique(grid_a[, c("label", "kind")])
    marks_a$kind <- factor(marks_a$kind, levels = names(pal))
    marks_a$pos <- -1.6
    pa <- ggplot(grid_a, aes(pos, label)) +
        geom_tile(aes(fill = v == 1L), colour = "white", linewidth = 0.4) +
        geom_point(data = marks_a, aes(pos, label, colour = kind), size = 2.4,
                   inherit.aes = FALSE) +
        scale_fill_manual(values = c("TRUE" = "#21918c", "FALSE" = "grey96"), guide = "none") +
        scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
        scale_x_continuous(breaks = seq(0, W_show - 1, by = 2), limits = c(-2.2, W_show - 0.5)) +
        labs(subtitle = kb_unicode("The space $X$: a filled cell is the symbol 1, and the dot on the left gives the classification"),
             x = kb_tex("Index $i$"), y = NULL) +
        kb_theme() + theme(panel.grid = element_blank()) +
        guides(colour = guide_legend(nrow = 3))
    pb <- ggplot(grid_b, aes(pos, label)) +
        geom_tile(aes(fill = v == 1L), colour = "white", linewidth = 0.4) +
        geom_text(data = marks, aes(pos, label, label = text), inherit.aes = FALSE,
                  size = 2.7, colour = "grey25", hjust = 0) +
        scale_fill_manual(values = c("TRUE" = "#21918c", "FALSE" = "grey96"), guide = "none") +
        scale_x_continuous(breaks = seq(0, W_show - 1, by = 2)) +
        labs(subtitle = kb_unicode("One forward orbit, which passes through $q$ and reaches the fixed point"),
             x = kb_tex("Index $i$"), y = NULL) +
        kb_theme() + theme(panel.grid = element_blank())
    fig <- (pa | pb) + plot_layout(widths = c(1.15, 1)) +
        plot_annotation(
            title = kb_unicode("A nonwandering point outside the limit set"),
            caption = kb_caption(paste(
                "Every element of $X$ holds at most two symbols 1, so no orbit visits the cylinder of words",
                "beginning with 1 more than twice and no omega limit set meets it; every orbit reaches the",
                "fixed point $0^{\\infty}$, so the limit set is that point alone. The point $q = 1 0^{\\infty}$",
                "lies in that cylinder, at distance 1 from the limit set, and is nonwandering all the same:",
                "the cylinder $[1 0^N]$ holds the generator $z_N$, and $T^{N+1}$ takes $z_N$ to $q$. The",
                "nonwandering set is $q$ together with $0^{\\infty}$, and the limit set is $0^{\\infty}$ alone.",
                "Drawn by checks/nonwandering-outside-the-limit-set.R.")),
            theme = kb_theme())
    kb_save(fig, "nonwandering-outside-the-limit-set", width = 9.0, height = 4.0)
}

emit("nonwandering-outside-the-limit-set", status,
     "On a countable compact shift space the limit set is one fixed point and the nonwandering set has two, so the inclusion of the limit set in the nonwandering set is strict",
     list(points_enumerated = t_inv$n_points,
          images_inside_the_space = as.numeric(t_inv$all_images_inside),
          greatest_number_of_ones = t_inv$max_ones,
          every_orbit_reaches_the_fixed_point = as.numeric(t_om$all_reach_zero),
          longest_run_to_the_fixed_point = t_om$max_steps,
          cylinders_tested_for_q = nrow(t_nw),
          witness_always_the_predicted_generator = as.numeric(all(t_nw$witness_is_the_predicted_generator)),
          size_of_the_nonwandering_set = t_wa$n_nonwandering_among_resolved,
          points_classified = t_wa$n_others_resolved + 2,
          cylinder_cap_set_by_the_truncation = t_wa$cap,
          longest_isolating_cylinder = t_wa$max_cylinder,
          size_of_the_limit_set = 1,
          distance_from_q_to_the_limit_set = dist_q_to_limit,
          pseudo_orbits_built = nrow(t_cr),
          largest_step_over_epsilon = max(t_cr$max_step / t_cr$eps),
          cylinders_without_a_witness_in_the_truncated_space = sum(!t_small$found)))
