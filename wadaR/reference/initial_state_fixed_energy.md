# Create initial state function for fixed energy Hamiltonian systems

Generates an initial state specification that constrains initial
conditions to a surface of constant energy for Hamiltonian systems with
separable kinetic energy.

## Usage

``` r
initial_state_fixed_energy(
  energy,
  potential_cpp,
  momentum_direction = c("x", "y", "radial")
)
```

## Arguments

- energy:

  Numeric. Target energy value.

- potential_cpp:

  Character. C++ code computing potential V(x,y). Should use variables
  `x`, `y` and return the potential value.

- momentum_direction:

  Character. Direction of initial momentum: "x" (p_y = 0), "y" (p_x =
  0), or "radial" (outward from origin).

## Value

A list with fields `energy`, `potential_cpp`, and `momentum_direction`
for use in compiled system setup.

## References

Henon, M., & Heiles, C. (1964). The applicability of the third integral
of motion: some numerical experiments. *The Astronomical Journal*, 69,
73-79. [doi:10.1086/109234](https://doi.org/10.1086/109234)

## Examples

``` r
if (FALSE) { # \dontrun{
# Henon-Heiles potential
init_spec <- initial_state_fixed_energy(
    energy = 0.2,
    potential_cpp = "0.5 * (x*x + y*y) + x*x*y - y*y*y/3.0",
    momentum_direction = "x"
)
} # }
```
