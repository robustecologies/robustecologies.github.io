# ============================ #
# Double precision along orbits of the logistic map at r = 4 ####
# ============================ #
#
# Claim, in the section "Ontology and epistemology" of perron-frobenius-operator: an orbit of
# S(x) = 4 x (1 - x) computed in double precision departs from the exact orbit by more than 0.1
# after about fifty iterations.
#
# Reference: the exact orbit, carried at 1700 bits, which is more than 500 significant decimal
# digits. The map is conjugate to the tent map, whose slopes are 2 and -2, so a rounding error
# of order 2^-53 grows to order one after about 53 steps. The high precision orbit loses about
# log10(2) = 0.30 digits per step and therefore stays exact to more than 380 digits over the
# 400 steps allowed here.
#
# This check was written in Python and was translated to R, with Rmpfr in place of the decimal
# module of the standard library. The starting values are the decimal strings, parsed once at
# 1700 bits, so the reference orbit differs from the exact decimal orbit by 2^-1700 at the
# start and by less than 2^-1300 at the end, which is far below the threshold tested.

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

if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    emit("logistic-orbit-precision", "error", "The package Rmpfr is not installed", list())
    quit(status = 0)
}

PREC <- 1700L
STARTS <- c("0.1", "0.3", "0.123456789", "0.7", "0.9")
THRESHOLD <- 0.1
MAX_STEPS <- 400L

# The gap between the two orbits at every step, for one starting value.
gaps_from <- function(start) {
    x_double <- as.numeric(start)
    x_exact <- Rmpfr::mpfr(start, precBits = PREC)
    four <- Rmpfr::mpfr(4L, precBits = PREC)
    one <- Rmpfr::mpfr(1L, precBits = PREC)
    gaps <- numeric(MAX_STEPS)
    for (n in seq_len(MAX_STEPS)) {
        x_double <- 4.0 * x_double * (1.0 - x_double)
        x_exact <- four * x_exact * (one - x_exact)
        gaps[n] <- abs(x_double - as.numeric(x_exact))
    }
    gaps
}

gaps <- lapply(STARTS, gaps_from)
names(gaps) <- STARTS
first_over <- vapply(gaps, function(g) {
    k <- which(g > THRESHOLD)
    if (length(k)) k[1] else NA_real_
}, numeric(1))

# 53 bits of precision and a doubling of the error at every step give about 53 steps; the band
# from 40 to 70 allows for the size of the initial rounding error and for the local stretching
# factor, which varies along the orbit.
passed <- all(is.finite(first_over)) && all(first_over >= 40 & first_over <= 70)
low <- min(first_over); high <- max(first_over)

metrics <- as.list(first_over)
names(metrics) <- paste0("steps_from_", gsub(".", "_", STARTS, fixed = TRUE))
emit("logistic-orbit-precision", if (passed) "pass" else "fail",
     sprintf("Double-precision orbits of the logistic map at r = 4 depart from the exact orbits by more than %s after %d to %d iterations",
             format(THRESHOLD), as.integer(low), as.integer(high)),
     metrics)

# ============================ #
# Figure ####
# ============================ #
#
# The gap between the double-precision orbit and the reference orbit, step by step, for the five
# starting values. The vertical axis is logarithmic, so the straight rise is the doubling of the
# error at every step, and the curves flatten at order one because both orbits stay in the unit
# interval and the gap can grow no further.

if (file.exists("checks/lib/figure-style.R") && requireNamespace("ggplot2", quietly = TRUE)) {
    source("checks/lib/figure-style.R")
    floor_gap <- 1e-18
    df <- do.call(rbind, lapply(STARTS, function(s)
        data.frame(n = seq_len(MAX_STEPS), gap = pmax(gaps[[s]], floor_gap), start = s)))
    df <- df[df$n <= 90, ]
    marks <- data.frame(start = STARTS, n = first_over,
                        gap = vapply(STARTS, function(s) gaps[[s]][first_over[[s]]], numeric(1)))
    growth <- data.frame(n = 1:90, gap = 2^(1:90) * .Machine$double.eps / 4)
    growth <- growth[growth$gap <= 2, ]        # clipped to the panel, so no points are dropped silently
    fig <- ggplot(df, aes(n, gap, colour = start)) +
        geom_line(data = growth, aes(n, gap), inherit.aes = FALSE,
                  colour = "grey55", linewidth = 0.4, linetype = "22") +
        geom_hline(yintercept = THRESHOLD, colour = "#B8390E", linewidth = 0.4) +
        geom_line(linewidth = 0.5) +
        geom_point(data = marks, aes(n, gap, colour = start), size = 1.5, show.legend = FALSE) +
        annotate("text", x = 2, y = 2e-1, hjust = 0, size = 2.6, colour = "#B8390E",
                 label = kb_unicode("Threshold $0.1$")) +
        annotate("text", x = 20, y = 1e-9, hjust = 0, size = 2.6, colour = "grey40",
                 label = "Grey: one rounding error of the double, doubled at every step") +
        scale_y_log10(limits = c(floor_gap, 2), breaks = 10^seq(-18, 0, by = 3)) +
        scale_colour_viridis_d(option = "viridis", end = 0.9, name = "Starting value") +
        labs(title = "Where a computed orbit stops being the orbit it started as",
             subtitle = sprintf("Gap between the double-precision orbit and a reference carried at %d bits, for five starting values", PREC),
             x = kb_tex("Step $n$"), y = "Gap between the two orbits") +
        kb_theme() +
        labs(caption = kb_caption(sprintf(
            "The points mark the first step at which the gap passes the threshold, which the run recorded by checks/logistic-orbit-precision.R placed between %d and %d. The map is conjugate to the tent map, whose slope has modulus 2 everywhere, so the initial rounding error of about $2^{-53}$ reaches order one after about 53 steps, which is the grey line. The invariant density of the same map is exact and is checked in checks/logistic-invariant-density.R.",
            as.integer(low), as.integer(high))))
    invisible(kb_save(fig, "logistic-orbit-precision", width = 8.4, height = 4.6))
}
