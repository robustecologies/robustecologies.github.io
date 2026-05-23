# Simulate 2D incompressible Navier-Stokes equations with the finite difference method

Solves the 2D incompressible Navier-Stokes equations using a finite
difference method with projection method for pressure-velocity coupling.
The equations are: \$\$\frac{\partial \mathbf{u}}{\partial t} +
(\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu
\nabla^2 \mathbf{u} + \mathbf{f}\$\$ \$\$\nabla \cdot \mathbf{u} = 0\$\$

## Usage

``` r
simulate_navier_stokes(
  nx,
  ny,
  lx,
  ly,
  dt,
  nt,
  nu,
  rho = 1,
  initial_condition = NULL,
  boundary_condition = NULL,
  forcing_function = NULL,
  pressure_solver = "jacobi",
  max_iter = NULL,
  tolerance = 0.001,
  save_interval = 1,
  save_final = TRUE,
  verbose = TRUE
)
```

## Arguments

- nx:

  Integer. Number of grid points in x-direction

- ny:

  Integer. Number of grid points in y-direction

- lx:

  Numeric. Domain length in x-direction

- ly:

  Numeric. Domain length in y-direction

- dt:

  Numeric. Time step size

- nt:

  Integer. Number of time steps

- nu:

  Numeric. Kinematic viscosity

- rho:

  Numeric. Fluid density (default = 1.0)

- initial_condition:

  Function. Function(x, y) returning list(u, v) for initial velocity
  field

- boundary_condition:

  List. Boundary conditions with elements:

  - type: "dirichlet", "neumann", "periodic", or "channel"

  - values: list with u_left, u_right, u_top, u_bottom, v_left, etc. The
    "channel" type uses periodic BCs in x (streamwise) and Dirichlet
    no-slip walls in y, which is the correct setup for Poiseuille flow.

- forcing_function:

  Function. External forcing f(x, y, t) returning list(fx, fy)

- pressure_solver:

  Character. Either "jacobi" (default, recommended) or "sor"

- max_iter:

  Integer. Maximum iterations for pressure solver. If NULL (default),
  automatically scales with grid size: max(1000, 20 \* max(nx, ny)).
  This ensures adequate iterations for Jacobi convergence on larger
  grids.

- tolerance:

  Numeric. Convergence tolerance for pressure (default = 1e-3)

- save_interval:

  Integer. Save solution every N steps (default = 1)

- save_final:

  Logical. If TRUE (default), always save the final state even if nt is
  not a multiple of save_interval

- verbose:

  Logical. Print progress messages to console (default = TRUE). Set to
  FALSE to suppress all informational output, which is required for CRAN
  compliance when used in non-interactive contexts.

## Value

A list containing:

- u:

  3D array (nx, ny, n_saves) of x-velocity

- v:

  3D array (nx, ny, n_saves) of y-velocity

- p:

  3D array (nx, ny, n_saves) of pressure

- x:

  Vector of x-coordinates

- y:

  Vector of y-coordinates

- t:

  Vector of time points corresponding to saved frames

- parameters:

  List of simulation parameters

- diagnostics:

  List with convergence info and mass conservation errors

## Details

The solver uses:

- Chorin's projection method for pressure-velocity coupling

- Second-order central differences for diffusion, divergence, and
  gradient

- First-order upwind scheme for convection (stable but introduces
  numerical diffusion of order \\O(\|u\| \Delta x / 2)\\)

- First-order explicit Euler for time integration

- Iterative solver (Jacobi or SOR) for Poisson pressure equation

- Rcpp-optimized computational kernels for performance

Overall accuracy is **first-order** in both space and time: \\O(\Delta
x) + O(\Delta t)\\. The upwind scheme trades accuracy for stability at
high Reynolds numbers.

The reported Reynolds number assumes a characteristic velocity of \\U =
1\\. For flows with different velocity scales, the actual Reynolds
number is \\Re = U \cdot l_x / \nu\\.

## References

Chorin, A. J. (1968). Numerical solution of the Navier-Stokes equations.
*Mathematics of Computation*, 22(104), 745-762. DOI:
[doi:10.1090/S0025-5718-1968-0242392-2](https://doi.org/10.1090/S0025-5718-1968-0242392-2)

Griebel, M., Dornseifer, T., & Neunhoeffer, T. (1998). *Numerical
simulation in fluid dynamics: A practical introduction*. SIAM. ISBN:
978-0-89871-398-5. DOI:
[doi:10.1137/1.9780898719703](https://doi.org/10.1137/1.9780898719703)

Ferziger, J. H., Peric, M., & Street, R. L. (2020). *Computational
methods for fluid dynamics* (4th ed.). Springer. ISBN:
978-3-319-99691-2. DOI:
[doi:10.1007/978-3-319-99693-6](https://doi.org/10.1007/978-3-319-99693-6)

## See also

[`print.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/print.navieRstokes.md),
[`summary.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/summary.navieRstokes.md),
[`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md),
[`flow.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/flow.md),
[`compute_vorticity`](https://robustecologies.github.io/navieRstokes/reference/compute_vorticity.md),
[`vortex_ic`](https://robustecologies.github.io/navieRstokes/reference/vortex_ic.md),
[`pressure_gradient_force`](https://robustecologies.github.io/navieRstokes/reference/pressure_gradient_force.md),
[`shiny_navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/shiny_navieRstokes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Lid-driven cavity flow simulation
result <- simulate_navier_stokes(
  nx = 64, ny = 64, lx = 1.0, ly = 1.0,
  dt = 0.001, nt = 1000, nu = 0.01,
  initial_condition = function(x, y) list(u = 0, v = 0),
  boundary_condition = list(
    type = "dirichlet",
    values = list(
      u_left = 0, u_right = 0, u_top = 1, u_bottom = 0,
      v_left = 0, v_right = 0, v_top = 0, v_bottom = 0
    )
  )
)

# Inspect results with S3 methods
print(result)
summary(result)

# Visualize flow field
plot(result)
plot(result, plot_type = "velocity")
plot(result, plot_type = "vorticity")
plot(result, plot_type = "pressure")

# View specific time step
plot(result, time_index = 50, plot_type = "velocity")
} # }
```
