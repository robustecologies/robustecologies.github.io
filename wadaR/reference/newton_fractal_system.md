# Newton fractal system

Creates a Newton-Raphson iteration map for finding roots of \\z^n - 1 =
0\\ in the complex plane. The resulting basins of attraction form the
classic Newton fractal, which exhibits Wada boundaries for \\n \geq 3\\.

## Usage

``` r
newton_fractal_system(n_roots = 3, relaxation = 1)
```

## Arguments

- n_roots:

  Integer. Degree of the polynomial (number of roots). Must be at least
  3 for Wada basins.

- relaxation:

  Numeric. Relaxation parameter \\\omega\\ for the Newton iteration.
  Default is 1 (standard Newton method). Values \\\omega \neq 1\\
  produce relaxed Newton variants.

## Value

A list containing:

- map:

  Function(state, iter) for the Newton iteration

- roots:

  Complex vector of the n-th roots of unity

- attractors:

  List of attractor specifications at each root

- params:

  List with n_roots and relaxation

- description:

  Human-readable description

## Details

The Newton-Raphson method for solving \\f(z) = z^n - 1 = 0\\ iterates:
\$\$z\_{k+1} = z_k - \omega \frac{z_k^n - 1}{n z_k^{n-1}}\$\$

The \\n\\ roots of unity are located at: \$\$\zeta_k = e^{2\pi i k / n},
\quad k = 0, 1, \ldots, n-1\$\$

These roots are equally spaced on the unit circle. Each root has its own
basin of attraction (the set of initial points that converge to that
root). The basins are separated by Julia set boundaries, which are
fractal.

For \\n \geq 3\\, the Newton fractal exhibits the Wada property: every
boundary point of one basin is simultaneously on the boundary of all
other basins. This was proven by Przytycki (1989).

## References

Przytycki, F. (1989). Remarks on the simple connectedness of basins of
sinks for iterations of rational maps. *Banach Center Publications*,
23(1), 229-235.
[doi:10.4064/-23-1-229-235](https://doi.org/10.4064/-23-1-229-235)

Curry, J. H., Garnett, L., & Sullivan, D. (1983). On the iteration of a
rational function: Computer experiments with Newton's method.
*Communications in Mathematical Physics*, 91(2), 267-277.
[doi:10.1007/BF01211162](https://doi.org/10.1007/BF01211162)

Blanchard, P. (1984). Complex analytic dynamics on the Riemann sphere.
*Bulletin of the American Mathematical Society*, 11(1), 85-141.
[doi:10.1090/S0273-0979-1984-15240-6](https://doi.org/10.1090/S0273-0979-1984-15240-6)

## See also

[`compute_newton_basins`](https://robustecologies.github.io/wadaR/reference/compute_newton_basins.md)
for efficient basin computation,
[`forced_damped_pendulum`](https://robustecologies.github.io/wadaR/reference/forced_damped_pendulum.md)
for a physical Wada system.

## Examples

``` r
if (FALSE) { # \dontrun{
# ===================================================================== #
# Example 1: Classic cubic Newton fractal (z^3 - 1 = 0)                 #
# ===================================================================== #
newton3 <- newton_fractal_system(n_roots = 3)
print(newton3$roots)  # Three roots of unity

# ===================================================================== #
# Example 2: Quintic Newton fractal (z^5 - 1 = 0)                       #
# ===================================================================== #
newton5 <- newton_fractal_system(n_roots = 5)
basins5 <- compute_newton_basins(n_roots = 5, resolution = 500)
plot(basins5, colors = c("#FB9E07", "#170C3A", "#2a1766"), 
     title = "Newton fractal: z^5 - 1 = 0")

# ===================================================================== #
# Example 3: Relaxed Newton method (omega = 0.8)                        #
# ===================================================================== #
newton_relaxed <- newton_fractal_system(n_roots = 3, relaxation = 0.8)
} # }
```
