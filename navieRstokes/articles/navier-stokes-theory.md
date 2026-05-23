# Navier-Stokes equations: mathematical foundations

## Mathematical theory and numerical methods for the Navier-Stokes equations

### Table of contents

1.  [Governing equations](#governing-equations)
2.  [Chorin’s projection method](#chorins-projection-method)
3.  [Spatial discretization](#spatial-discretization)
4.  [Temporal discretization](#temporal-discretization)
5.  [Poisson solver](#poisson-solver)
6.  [Stability analysis](#stability-analysis)
7.  [Boundary conditions](#boundary-conditions)
8.  [Verification and validation](#verification-and-validation)

------------------------------------------------------------------------

### Governing equations

#### The incompressible Navier-Stokes equations

The motion of an incompressible Newtonian fluid is governed by the
Navier-Stokes equations, which express conservation of momentum and mass
[\[5\]](#ref5),[\[8\]](#ref8):

**Momentum equation:** \\ \frac{\partial \mathbf{u}}{\partial t} +
(\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu
\nabla^2 \mathbf{u} + \mathbf{f} \\

**Continuity equation (incompressibility):** \\ \nabla \cdot \mathbf{u}
= 0 \\

where: - \\\mathbf{u} = (u, v)\\ is the velocity field in 2D - \\p\\ is
the pressure field - \\\rho\\ is the constant fluid density - \\\nu =
\mu/\rho\\ is the kinematic viscosity - \\\mathbf{f}\\ is the external
body force per unit mass

#### Physical interpretation

Each term in the momentum equation has a physical meaning:

- \\\frac{\partial \mathbf{u}}{\partial t}\\: Local acceleration
  (unsteady term)
- \\(\mathbf{u} \cdot \nabla)\mathbf{u}\\: Convective acceleration
  (nonlinear term)
- \\-\frac{1}{\rho}\nabla p\\: Pressure gradient force
- \\\nu \nabla^2 \mathbf{u}\\: Viscous diffusion
- \\\mathbf{f}\\: External forcing (gravity, electromagnetic forces,
  etc.)

#### Reynolds number

The dimensionless **Reynolds number** characterizes the ratio of
inertial to viscous forces [\[5\]](#ref5):

\\ Re = \frac{U L}{\nu} \\

where \\U\\ is a characteristic velocity and \\L\\ is a characteristic
length scale. The flow regime depends critically on \\Re\\:

- \\Re \ll 1\\: Stokes flow (viscous-dominated, linear)
- \\Re \sim 1\\: Transitional regime
- \\Re \gg 1\\: Inertial flow (potentially turbulent)

**Package convention:** The
[`simulate_navier_stokes()`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
function reports \\Re = l_x / \nu\\, which assumes a characteristic
velocity \\U = 1\\. For lid-driven cavity with `u_top = 1`, this is
exact. For other flows (e.g., Taylor-Green with amplitude \\A\\), the
actual Reynolds number is \\Re = A \cdot l_x / \nu\\.

------------------------------------------------------------------------

### Chorin’s projection method

#### The pressure-velocity coupling problem

The primary difficulty in solving the incompressible Navier-Stokes
equations is the **pressure-velocity coupling**: the pressure appears as
a constraint to enforce incompressibility \\(\nabla \cdot \mathbf{u} =
0)\\, but there is no explicit evolution equation for \\p\\.

#### Helmholtz-Hodge decomposition

**Theorem 1 (Helmholtz-Hodge decomposition):** Any sufficiently smooth
vector field \\\mathbf{w}\\ defined on a simply-connected domain
\\\Omega\\ with appropriate boundary conditions can be uniquely
decomposed as:

\\ \mathbf{w} = \mathbf{u} + \nabla \phi \\

where \\\nabla \cdot \mathbf{u} = 0\\ (divergence-free part) and
\\\phi\\ is a scalar potential.

This theorem is the foundation of projection methods.

#### Chorin’s fractional-step scheme

Chorin’s method [\[1\]](#ref1) decouples the momentum equation into two
sub-steps:

##### Step 1: Advection-diffusion (predictor step)

Compute an intermediate velocity \\\mathbf{u}^\*\\ by solving the
momentum equation **without** the pressure term:

\\ \frac{\mathbf{u}^\* - \mathbf{u}^n}{\Delta t} = -(\mathbf{u}^n \cdot
\nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n + \mathbf{f}^n \\

This gives: \\ \mathbf{u}^\* = \mathbf{u}^n + \Delta t
\left\[-(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2
\mathbf{u}^n + \mathbf{f}^n\right\] \\

Note that \\\mathbf{u}^\*\\ generally does **not** satisfy the
divergence-free constraint.

##### Step 2: Projection (corrector step)

Project \\\mathbf{u}^\*\\ onto the space of divergence-free vector
fields. The correction is:

\\ \mathbf{u}^{n+1} = \mathbf{u}^\* - \frac{\Delta t}{\rho} \nabla
p^{n+1} \\

Applying the divergence operator to both sides and enforcing \\\nabla
\cdot \mathbf{u}^{n+1} = 0\\:

\\ 0 = \nabla \cdot \mathbf{u}^\* - \frac{\Delta t}{\rho} \nabla^2
p^{n+1} \\

This yields the **pressure Poisson equation**:

\\ \nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^\*
\\

##### Summary of the algorithm

1.  **Predictor**: Solve for \\\mathbf{u}^\*\\ from advection-diffusion
2.  **Poisson**: Solve \\\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla
    \cdot \mathbf{u}^\*\\
3.  **Corrector**: Update \\\mathbf{u}^{n+1} = \mathbf{u}^\* -
    \frac{\Delta t}{\rho} \nabla p^{n+1}\\

**Theorem 2 (Convergence of Chorin’s method):** For smooth solutions,
Chorin’s projection method is **first-order accurate in time** for the
velocity and pressure.

------------------------------------------------------------------------

### Spatial discretization

#### Grid structure

We use a **uniform Cartesian collocated grid**
[\[3\]](#ref3),[\[4\]](#ref4):

\\ x_i = (i-1)\Delta x, \quad i = 1, \ldots, n_x, \quad \Delta x =
\frac{L_x}{n_x - 1} \\ \\ y_j = (j-1)\Delta y, \quad j = 1, \ldots, n_y,
\quad \Delta y = \frac{L_y}{n_y - 1} \\

All variables \\(u, v, p)\\ are stored at the same grid locations
(collocated grid).

**Note:** More sophisticated implementations use **staggered grids**
(MAC grid) where velocities are stored at cell faces and pressure at
cell centers. This eliminates spurious pressure oscillations
(checkerboard modes) but complicates the discretization.

#### Finite difference operators

The discretization follows standard finite difference methods
[\[6\]](#ref6),[\[7\]](#ref7).

##### Second-order central differences

For the diffusion term (Laplacian):

\\ \left.\frac{\partial^2 u}{\partial x^2}\right\|\_{i,j} \approx
\frac{u\_{i+1,j} - 2u\_{i,j} + u\_{i-1,j}}{\Delta x^2} + O(\Delta x^2)
\\

\\ \left.\nabla^2 u\right\|\_{i,j} \approx \frac{u\_{i+1,j} -
2u\_{i,j} + u\_{i-1,j}}{\Delta x^2} + \frac{u\_{i,j+1} - 2u\_{i,j} +
u\_{i,j-1}}{\Delta y^2} + O(\Delta x^2, \Delta y^2) \\

For the gradient (pressure):

\\ \left.\frac{\partial p}{\partial x}\right\|\_{i,j} \approx
\frac{p\_{i+1,j} - p\_{i-1,j}}{2\Delta x} + O(\Delta x^2) \\

For the divergence:

\\ \left.\nabla \cdot \mathbf{u}\right\|\_{i,j} \approx
\frac{u\_{i+1,j} - u\_{i-1,j}}{2\Delta x} + \frac{v\_{i,j+1} -
v\_{i,j-1}}{2\Delta y} + O(\Delta x^2, \Delta y^2) \\

##### First-order upwind differences

For the convective term \\(\mathbf{u} \cdot \nabla)u = u \frac{\partial
u}{\partial x} + v \frac{\partial u}{\partial y}\\, we use **upwind
differencing** for stability:

\\ u\_{i,j} \frac{\partial u}{\partial x}\bigg\|\_{i,j} \approx u\_{i,j}
\times \begin{cases} \frac{u\_{i,j} - u\_{i-1,j}}{\Delta x} & \text{if }
u\_{i,j} \> 0 \text{ (backward)} \\ \frac{u\_{i+1,j} - u\_{i,j}}{\Delta
x} & \text{if } u\_{i,j} \< 0 \text{ (forward)} \end{cases} \\

This choice ensures that information propagates in the direction of the
flow, maintaining stability at high Reynolds numbers.

**Accuracy:** Upwind differencing is **first-order accurate** \\O(\Delta
x)\\ but provides numerical dissipation that stabilizes the convective
term.

**Important:** Because the convection term uses first-order upwind and
time integration uses first-order Euler, the **overall scheme accuracy
is \\O(\Delta x) + O(\Delta t)\\**, despite diffusion/gradient operators
being second-order.

------------------------------------------------------------------------

### Temporal discretization

#### Explicit time stepping

The advection-diffusion step uses **explicit forward Euler** time
integration:

\\ \mathbf{u}^\* = \mathbf{u}^n + \Delta t \cdot
\text{RHS}(\mathbf{u}^n) \\

where RHS includes convection, diffusion, and forcing terms.

**Advantages:** - Simple to implement - Low memory requirements - No
need to solve linear systems for velocity

**Disadvantages:** - Subject to strict stability constraints (CFL
conditions) - Requires small time steps for high Reynolds number or fine
grids

#### Implicit alternatives

For better stability, implicit or semi-implicit schemes can be used: -
**Crank-Nicolson** (second-order): \\\frac{\mathbf{u}^{n+1} -
\mathbf{u}^n}{\Delta t} = \frac{1}{2}\[\text{RHS}(\mathbf{u}^{n+1}) +
\text{RHS}(\mathbf{u}^n)\]\\ - **Implicit Euler** (first-order,
unconditionally stable)

These require solving linear systems at each time step but allow larger
time steps.

------------------------------------------------------------------------

### Poisson solver

#### The discrete Poisson equation

The pressure Poisson equation on the discrete grid is:

\\ \frac{p\_{i+1,j} - 2p\_{i,j} + p\_{i-1,j}}{\Delta x^2} +
\frac{p\_{i,j+1} - 2p\_{i,j} + p\_{i,j-1}}{\Delta y^2} =
\text{RHS}\_{i,j} \\

where \\\text{RHS}\_{i,j} = \frac{\rho}{\Delta t}(\nabla \cdot
\mathbf{u}^\*)\_{i,j}\\.

Rearranging for the unknown \\p\_{i,j}\\:

\\ p\_{i,j} = \frac{\Delta y^2(p\_{i+1,j} + p\_{i-1,j}) + \Delta
x^2(p\_{i,j+1} + p\_{i,j-1}) - \Delta x^2 \Delta y^2 \\
\text{RHS}\_{i,j}}{2(\Delta x^2 + \Delta y^2)} \\

#### Iterative solvers

##### Jacobi method

The Jacobi method [\[3\]](#ref3) updates all points simultaneously using
values from the previous iteration:

\\ p\_{i,j}^{(k+1)} = \frac{\Delta y^2(p\_{i+1,j}^{(k)} +
p\_{i-1,j}^{(k)}) + \Delta x^2(p\_{i,j+1}^{(k)} + p\_{i,j-1}^{(k)}) -
\Delta x^2 \Delta y^2 \\ \text{RHS}\_{i,j}}{2(\Delta x^2 + \Delta y^2)}
\\

**Convergence:** Linear convergence, spectral radius \\\rho \approx 1 -
O(\Delta x^2)\\.

##### Successive over-relaxation (SOR)

SOR [\[3\]](#ref3),[\[4\]](#ref4) accelerates convergence by using
updated values immediately and over-relaxing:

\\ p\_{i,j}^{(k+1)} = (1-\omega) p\_{i,j}^{(k)} + \omega \cdot
p\_{i,j}^{\text{GS}} \\

where \\p\_{i,j}^{\text{GS}}\\ is the Gauss-Seidel update (using
already-updated neighbors) and \\\omega\\ is the relaxation parameter.

**Optimal relaxation parameter:** For a rectangular grid,

\\ \omega\_{\text{opt}} = \frac{2}{1 + \sin\left(\frac{\pi}{\max(n_x,
n_y)}\right)} \\

This choice minimizes the spectral radius and significantly accelerates
convergence.

**Convergence:** With optimal \\\omega\\, SOR converges much faster than
Jacobi, especially for large grids.

#### Practical convergence considerations

In practice, the choice of convergence tolerance depends on several
factors:

1.  **Truncation error balance**: There’s no benefit in solving the
    Poisson equation to accuracy much better than the spatial
    discretization error \\O(\Delta x^2)\\. For a \\64 \times 64\\ grid
    with \\\Delta x \approx 0.016\\, this suggests tolerance \\\sim
    10^{-4}\\ to \\10^{-5}\\.

2.  **Mass conservation**: The pressure correction step ensures \\\nabla
    \cdot \mathbf{u}^{n+1} \approx 0\\. The divergence error is
    proportional to the Poisson solver error. For most applications,
    \\\|\nabla \cdot \mathbf{u}\| \< 10^{-5}\\ is excellent mass
    conservation.

3.  **Computational cost**: Tighter tolerances (e.g., \\10^{-8}\\)
    require many more iterations with diminishing returns for overall
    solution accuracy.

**Recommendation**: Use `tolerance = 1e-5` for most simulations. This
provides good mass conservation while maintaining computational
efficiency.

#### Boundary conditions for pressure

We apply **Neumann boundary conditions** \\\frac{\partial p}{\partial n}
= 0\\ on all boundaries:

\\ p\_{1,j} = p\_{2,j}, \quad p\_{n_x,j} = p\_{n_x-1,j}, \quad p\_{i,1}
= p\_{i,2}, \quad p\_{i,n_y} = p\_{i,n_y-1} \\

This is derived from the normal component of the momentum equation at
the boundary.

**Note:** The Poisson equation with pure Neumann boundary conditions has
a solution unique up to an additive constant (the pressure is determined
only up to a reference value).

------------------------------------------------------------------------

### Stability analysis

#### CFL conditions

Explicit time stepping requires satisfying **Courant-Friedrichs-Lewy
(CFL)** stability conditions [\[3\]](#ref3),[\[6\]](#ref6).

##### Diffusion stability

For the diffusion term \\\nu \nabla^2 \mathbf{u}\\ with explicit Euler,
the stability condition in 2D is:

\\ \Delta t \leq \frac{\Delta x^2}{2 d \nu} \\

where \\d = 2\\ is the spatial dimension. This gives:

\\ \Delta t \leq \frac{\min(\Delta x, \Delta y)^2}{4\nu} \\

**Physical interpretation:** Diffusion propagates information at a rate
proportional to \\\nu / \Delta x^2\\. The time step must resolve this
diffusion timescale.

##### Convection stability

For the convective term \\(\mathbf{u} \cdot \nabla)\mathbf{u}\\, the CFL
condition is:

\\ \text{CFL} = \frac{\|\mathbf{u}\|\_{\max} \Delta t}{\min(\Delta x,
\Delta y)} \< 1 \\

or equivalently:

\\ \Delta t \< \frac{\min(\Delta x, \Delta y)}{\|\mathbf{u}\|\_{\max}}
\\

**Physical interpretation:** A fluid particle should not traverse more
than one grid cell in a single time step.

##### Combined condition

The overall stability requirement is:

\\ \Delta t \< \min\left\\\frac{\Delta x^2}{4\nu}, \\ \frac{\min(\Delta
x, \Delta y)}{\|\mathbf{u}\|\_{\max}}\right\\ \\

At **high Reynolds numbers** (\\Re \gg 1\\), convection dominates and
the convective CFL is more restrictive. At **low Reynolds numbers**
(\\Re \ll 1\\), diffusion dominates.

#### Accuracy vs. stability

The upwind scheme for convection introduces **numerical diffusion** of
order \\O(\|\mathbf{u}\| \Delta x / 2)\\, which acts as an artificial
viscosity. This stabilizes the scheme but reduces accuracy for
under-resolved flows.

**Numerical viscosity estimate:**

\\ \nu\_{\text{num}} \approx \frac{\|\mathbf{u}\| \Delta x}{2} \\

**Effective Reynolds number:**

\\ Re\_{\text{eff}} = \frac{U L}{\nu + \nu\_{\text{num}}} =
\frac{Re}{1 + \frac{\|\mathbf{u}\| \Delta x}{2\nu}} \\

For accurate simulations, we require \\\|\mathbf{u}\| \Delta x \ll
2\nu\\ (well-resolved grid). When this condition is violated, the
numerical diffusion dominates physical diffusion and the effective
viscosity increases significantly.

**Practical guideline:** If \\Re \cdot \Delta x / L \> 1\\, the
numerical diffusion exceeds the physical diffusion. For example, on a
\\64 \times 64\\ grid (\\\Delta x \approx 0.016\\), simulations at \\Re
\> 60\\ will exhibit significant numerical diffusion effects.

------------------------------------------------------------------------

### Boundary conditions

#### Dirichlet boundary conditions (specified values)

For velocity, we directly impose values at the boundary:

\\ u(0, y) = u\_{\text{left}}, \quad u(L_x, y) = u\_{\text{right}},
\quad \text{etc.} \\

**Example:** No-slip walls \\(\mathbf{u} = 0)\\, moving lid \\(u =
U\_{\text{lid}}, v = 0)\\.

#### Neumann boundary conditions (zero gradient)

For outflow or symmetry boundaries:

\\ \frac{\partial u}{\partial n} = 0 \quad \Rightarrow \quad u\_{1,j} =
u\_{2,j}, \quad u\_{n_x,j} = u\_{n_x-1,j} \\

#### Periodic boundary conditions

For periodic domains (e.g., channel flow):

\\ u\_{1,j} = u\_{n_x-1,j}, \quad u\_{n_x,j} = u\_{2,j} \\

#### Pressure boundary conditions

For all boundary types, we impose **Neumann conditions** on pressure:

\\ \frac{\partial p}{\partial n} = 0 \\

This is consistent with the momentum equation at the boundary and
ensures well-posedness of the Poisson problem.

------------------------------------------------------------------------

### Verification and validation

#### Analytical solutions

##### Stokes flow between parallel plates

For low Reynolds number (\\Re \ll 1\\), the momentum equation reduces
to:

\\ 0 = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu
\frac{\partial^2 u}{\partial y^2} \\

With boundary conditions \\u(0) = u(h) = 0\\ and constant pressure
gradient \\\frac{\partial p}{\partial x} = G\\, the solution is:

\\ u(y) = \frac{G}{2\nu\rho} y(h - y) \\

This is the classic **Poiseuille flow** and can be used to verify the
diffusion solver.

##### Taylor-Green vortex

An exact solution to the 2D Navier-Stokes equations (decaying vortices).
The **canonical form** on domain \\\[0, 2\pi\] \times \[0, 2\pi\]\\ is:

\\ u(x, y, t) = -\cos(x) \sin(y) \\ e^{-2\nu t} \\ \\ v(x, y, t) =
\sin(x) \cos(y) \\ e^{-2\nu t} \\

**This package uses a different (but equivalent) sign convention** on
the unit domain \\\[0, 1\] \times \[0, 1\]\\:

\\ u(x, y, t) = A \sin(2\pi x) \cos(2\pi y) \\ e^{-2\nu k^2 t} \\ \\
v(x, y, t) = -A \cos(2\pi x) \sin(2\pi y) \\ e^{-2\nu k^2 t} \\

where \\k = 2\pi\\. Both forms are divergence-free and represent valid
Taylor-Green vortices with opposite rotation sense.

The **kinetic energy** decays exponentially: \\ E(t) = E_0 \\ e^{-4\nu
k^2 t} \\

**Verification:** Simulate this flow and compare with the analytical
solution to verify: - Temporal accuracy (convergence rate in \\\Delta
t\\) - Spatial accuracy (convergence rate in \\\Delta x\\) - Energy
decay rate

#### Benchmark problems

##### Lid-driven cavity

A square cavity with no-slip walls and a moving lid. Well-documented
benchmark results exist for various Reynolds numbers (Ghia et al.,
1982).

**Key diagnostics:** - Velocity profiles along centerlines - Primary
vortex center location - Secondary vortex formation at \\Re \> 100\\

See [\[2\]](#ref2) for benchmark data.

##### Flow past a cylinder

Vortex shedding behind a circular cylinder exhibits periodic behavior
characterized by the **Strouhal number**:

\\ St = \frac{f D}{U} \\

where \\f\\ is the shedding frequency, \\D\\ is the cylinder diameter,
and \\U\\ is the free-stream velocity.

------------------------------------------------------------------------

### Computational complexity

#### Time per iteration

- **Convection/diffusion**: \\O(n_x \cdot n_y)\\ (explicit updates)
- **Poisson solver**: \\O(n\_{\text{iter}} \cdot n_x \cdot n_y)\\ where
  \\n\_{\text{iter}}\\ depends on grid size and solver

For SOR with optimal \\\omega\\, \\n\_{\text{iter}} = O(\sqrt{n_x \cdot
n_y})\\.

**Total per time step:** \\O(n_x \cdot n_y \cdot \sqrt{n_x \cdot n_y}) =
O((n_x n_y)^{3/2})\\

#### Memory requirements

Storage for three 2D fields at multiple time steps:

\\ \text{Memory} \approx 3 \times n_x \times n_y \times
n\_{\text{saves}} \times 8 \text{ bytes} \\

For a \\128 \times 128\\ grid with 100 saved time steps: \\\approx 40\\
MB.

------------------------------------------------------------------------

### Extensions and improvements

#### Higher-order schemes

- **Second-order time stepping**: Use Runge-Kutta or Crank-Nicolson
- **Higher-order spatial schemes**: QUICK, ENO, WENO for convection
- **Multigrid solvers**: Reduce Poisson solver cost to \\O(n_x \cdot
  n_y)\\

#### Advanced physics

- **Variable viscosity**: Temperature-dependent \\\nu(T)\\
- **Non-Newtonian fluids**: Shear-thinning/thickening
- **Multiphase flows**: Level-set or volume-of-fluid methods
- **Free surfaces**: Coupled with interface tracking

#### Turbulence modeling

For high Reynolds numbers where direct simulation is infeasible: -
**Large eddy simulation (LES)**: Resolve large scales, model small
scales - **Reynolds-averaged Navier-Stokes (RANS)**: Ensemble-averaged
equations with turbulence closure models

------------------------------------------------------------------------

### Appendix: Derivation of the Navier-Stokes equations

#### Conservation of momentum

From Newton’s second law applied to a fluid element:

\\ \rho \frac{D\mathbf{u}}{Dt} = \nabla \cdot \mathbf{\sigma} + \rho
\mathbf{f} \\

where \\\frac{D}{Dt} = \frac{\partial}{\partial t} + \mathbf{u} \cdot
\nabla\\ is the material derivative and \\\mathbf{\sigma}\\ is the
stress tensor.

For a **Newtonian fluid**:

\\ \mathbf{\sigma} = -p \mathbf{I} + \mu (\nabla \mathbf{u} + (\nabla
\mathbf{u})^T) \\

For incompressible flow (\\\nabla \cdot \mathbf{u} = 0\\):

\\ \nabla \cdot \mathbf{\sigma} = -\nabla p + \mu \nabla^2 \mathbf{u} \\

Substituting and expanding the material derivative:

\\ \rho \left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot
\nabla)\mathbf{u}\right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho
\mathbf{f} \\

Dividing by \\\rho\\ and using \\\nu = \mu/\rho\\:

\\ \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot
\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u} +
\mathbf{f} \\

#### Conservation of mass

For incompressible flow (\\\rho = \text{const}\\), the continuity
equation reduces to:

\\ \nabla \cdot \mathbf{u} = 0 \\

This is the **divergence-free constraint** that couples the velocity and
pressure.

------------------------------------------------------------------------

### Appendix: Vorticity-streamfunction formulation

An alternative formulation eliminates pressure by taking the curl of the
momentum equation.

#### Vorticity transport equation

Define the **vorticity** \\\omega = \nabla \times \mathbf{u}\\. In 2D:

\\ \omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial
y} \\

Taking the curl of the momentum equation (assuming \\\mathbf{f} = 0\\):

\\ \frac{\partial \omega}{\partial t} + (\mathbf{u} \cdot \nabla)\omega
= \nu \nabla^2 \omega \\

The pressure term vanishes because \\\nabla \times \nabla p = 0\\.

#### Streamfunction

Define the **streamfunction** \\\psi\\ such that:

\\ u = \frac{\partial \psi}{\partial y}, \quad v = -\frac{\partial
\psi}{\partial x} \\

This automatically satisfies \\\nabla \cdot \mathbf{u} = 0\\.

The vorticity is related to the streamfunction by:

\\ \omega = -\nabla^2 \psi \\

#### Coupled system

The equations become:

1.  **Vorticity transport**: \\\frac{\partial \omega}{\partial t} +
    (\mathbf{u} \cdot \nabla)\omega = \nu \nabla^2 \omega\\
2.  **Poisson equation for streamfunction**: \\\nabla^2 \psi = -\omega\\
3.  **Velocity recovery**: \\u = \frac{\partial \psi}{\partial y}, \\ v
    = -\frac{\partial \psi}{\partial x}\\

**Advantages:** No pressure variable, automatic incompressibility.

**Disadvantages:** Difficult to impose boundary conditions, does not
extend easily to 3D.

------------------------------------------------------------------------

### References

**\[1\]** Chorin, A. J. (1968). Numerical solution of the Navier-Stokes
equations. *Mathematics of Computation*, 22(104), 745-762. DOI:
[10.1090/S0025-5718-1968-0242392-2](https://doi.org/10.1090/S0025-5718-1968-0242392-2)

**\[2\]** Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions
for incompressible flow using the Navier-Stokes equations and a
multigrid method. *Journal of Computational Physics*, 48(3), 387-411.
DOI:
[10.1016/0021-9991(82)90058-4](https://doi.org/10.1016/0021-9991(82)90058-4)

**\[3\]** Griebel, M., Dornseifer, T., & Neunhoeffer, T. (1998).
*Numerical simulation in fluid dynamics: A practical introduction*.
SIAM. ISBN: 978-0-89871-398-5. DOI:
[10.1137/1.9780898719703](https://doi.org/10.1137/1.9780898719703)

**\[4\]** Ferziger, J. H., & Peric, M. (2002). *Computational methods
for fluid dynamics* (3rd ed.). Springer. ISBN: 978-3-540-42074-3. DOI:
[10.1007/978-3-642-56026-2](https://doi.org/10.1007/978-3-642-56026-2)

**\[5\]** Kundu, P. K., Cohen, I. M., & Dowling, D. R. (2015). *Fluid
mechanics* (6th ed.). Academic Press. ISBN: 978-0-12-405935-1.

**\[6\]** Strikwerda, J. C. (2004). *Finite difference schemes and
partial differential equations* (2nd ed.). SIAM. ISBN:
978-0-89871-567-5. DOI:
[10.1137/1.9780898717938](https://doi.org/10.1137/1.9780898717938)

**\[7\]** Quarteroni, A., & Valli, A. (1994). *Numerical approximation
of partial differential equations*. Springer. ISBN: 978-3-540-57111-2.
DOI:
[10.1007/978-3-540-85268-1](https://doi.org/10.1007/978-3-540-85268-1)

**\[8\]** Temam, R. (2001). *Navier-Stokes equations: Theory and
numerical analysis*. AMS Chelsea Publishing. ISBN: 978-0-8218-2737-6.
