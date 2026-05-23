# Historical development of the Navier-Stokes equations

## Historical development of the Navier-Stokes equations

### Table of contents

1.  [Claude-Louis Navier and the molecular approach](#navier-molecular)
2.  [The continuum reformulations: Cauchy and
    Poisson](#continuum-reformulations)
3.  [Saint-Venant’s phenomenological contribution](#saint-venant)
4.  [George Gabriel Stokes and experimental
    validation](#stokes-validation)
5.  [Evolution of mathematical notation](#evolution-notation)
6.  [Contemporary mathematical challenges](#contemporary-challenges)
7.  [Computational developments and
    applications](#computational-developments)

------------------------------------------------------------------------

### Claude-Louis Navier and the molecular approach

#### Biographical context

The genesis of the Navier-Stokes equations emerged from
post-Revolutionary France. Claude-Louis Marie Henri Navier (1785-1836)
was orphaned at age eight when his father, a lawyer and member of the
National Assembly, died during the Terror of 1793 [\[10\]](#ref10).
Raised by his granduncle Émiland Gauthey, a prominent civil engineer,
Navier entered the École Polytechnique in 1802, ranking near the bottom
of his cohort [\[10\]](#ref10). Under Joseph Fourier’s mentorship, he
developed conviction in mathematical analysis as the fundamental tool
for engineering problems [\[10\]](#ref10).

#### The 1822 derivation

In 1822, Navier presented his seminal memoir *Sur les lois des
mouvements des fluides, en ayant égard à l’adhésion des molecules* to
the Paris Academy of Sciences [\[1\]](#ref1),[\[3\]](#ref3). His
derivation remains historically significant because he obtained correct
equations from an incorrect physical model. Influenced by Laplacian
physics, Navier conceptualized fluids as collections of discrete
molecules interacting through short-range forces [\[3\]](#ref3).

Navier modeled “adhesion” between molecules as forces arising from
relative particle displacement [\[1\]](#ref1),[\[3\]](#ref3). By
averaging these molecular interactions, he successfully derived the
viscous resistance terms, introducing the Laplacian operator
\\\nabla^2\\ into the momentum equation [\[2\]](#ref2). The molecular
lattice model, while conceptually flawed by modern standards, captured
the macroscopic phenomenology of shear resistance through intuitive
reasoning about invisible intermolecular attractions [\[3\]](#ref3).

#### Professional setbacks

Navier’s theoretical work received limited recognition during his
lifetime. His derivation appeared mathematically complex and physically
speculative to both engineers and mathematicians [\[1\]](#ref1). His
career suffered further when his suspension bridge over the Seine at
Invalides developed structural problems due to soil settlement and was
dismantled before completion in the late 1820s [\[10\]](#ref10). This
failure reinforced skepticism from empirically-oriented British
engineers regarding mathematical approaches to practical problems
[\[10\]](#ref10).

------------------------------------------------------------------------

### The continuum reformulations: Cauchy and Poisson

#### Augustin-Louis Cauchy’s stress tensor approach

Augustin-Louis Cauchy (1789-1857) revisited the problem in 1823 and 1828
from a fundamentally different perspective [\[1\]](#ref1). Abandoning
specific molecular models, Cauchy introduced the **stress tensor**
concept, describing forces on internal surfaces of deformable bodies
[\[1\]](#ref1),[\[2\]](#ref2). This continuum-based approach derived
equations for “perfectly inelastic solids,” arguing that stress in
fluids depends on strain rate rather than strain itself [\[1\]](#ref1).

Cauchy’s formulation represented a profound epistemological shift:
material behavior could be described without knowledge of atomic
constitution [\[1\]](#ref1). By defining stress \\\sigma\\ and strain
rate \\\dot{\epsilon}\\, he demonstrated that viscosity arises from
fluid element deformation [\[2\]](#ref2). However, Cauchy maintained
limited interest in physical applications, viewing his equations
primarily as mathematical structures [\[1\]](#ref1).

#### Siméon Denis Poisson’s molecular refinement

Siméon Denis Poisson (1781-1840) independently derived the fluid motion
equations in 1829, defending the molecular worldview against Cauchy’s
abstraction [\[1\]](#ref1). Poisson’s derivation employed sophisticated
intermolecular force models, mathematically more robust than Navier’s
approach [\[1\]](#ref1),[\[3\]](#ref3). From his molecular framework, he
proposed specific relationships between viscosity coefficients
(analogous to what later became known as the Stokes hypothesis), arguing
that fluid behavior could be reduced to a single constant—valid for many
gases but inadequate for liquids with complex molecular structures
[\[1\]](#ref1).

The period 1822-1845 witnessed five independent derivations: Navier
(1822), Cauchy (1823), Poisson (1829), Saint-Venant (1837/1843), and
Stokes (1845) [\[1\]](#ref1). This redundancy reflected limited
scientific communication and intense personal rivalries within the Paris
Academy [\[1\]](#ref1).

------------------------------------------------------------------------

### Saint-Venant’s phenomenological contribution

#### The forgotten pioneer

Adhémar Jean Claude Barré de Saint-Venant (1797-1886) published a
crucial derivation in 1843, two years before Stokes, correctly
identifying viscosity’s physical nature without molecular speculation
[\[1\]](#ref1),[\[3\]](#ref3). As a hydraulic engineer, Saint-Venant
recognized that summing individual atomic forces was mathematically
intractable [\[3\]](#ref3). He instead adopted a phenomenological
approach, postulating that internal fluid friction acts tangentially to
motion and is proportional to velocity gradient
[\[1\]](#ref1),[\[3\]](#ref3).

#### First correct viscosity identification

Saint-Venant’s 1843 paper provided the first accurate identification of
the coefficient of viscosity \\\mu\\ as a macroscopic fluid property,
effectively defining shear stress as \\\tau = \mu \frac{du}{dy}\\
[\[1\]](#ref1),[\[2\]](#ref2). He argued that stress within fluids
results from sliding of fluid layers—a concept consistent with modern
understanding of momentum diffusion [\[3\]](#ref3). This representation
aligned viscosity directly with observable macroscopic principles rather
than hypothetical molecular arrangements [\[1\]](#ref1).

#### Historical obscurity

Despite this fundamental contribution, Saint-Venant’s name remains
absent from the equation’s designation. Historical factors include his
modest temperament, publication in French journals as scientific
influence shifted toward Britain and Germany, and involvement in
priority disputes regarding vector calculus with Hermann Grassmann
[\[1\]](#ref1),[\[10\]](#ref10). His legacy remains compartmentalized in
elasticity theory (Saint-Venant’s principle) and shallow water flow
(Saint-Venant equations) rather than general fluid mechanics
[\[1\]](#ref1).

------------------------------------------------------------------------

### George Gabriel Stokes and experimental validation

#### The continuum synthesis

George Gabriel Stokes (1819-1903), born into an Irish clerical family,
held the Lucasian Chair of Mathematics at Cambridge—previously occupied
by Isaac Newton [\[10\]](#ref10). In 1845, at age 26, he published *On
the theories of the internal friction of fluids in motion*
[\[4\]](#ref4). While aware of French derivations, Stokes sought a
formulation relying exclusively on observable macroscopic principles
[\[4\]](#ref4).

#### Isotropic continuum mechanics

Stokes’s derivation exemplified rigorous continuum mechanics. He assumed
**isotropy**—identical fluid properties in all directions—and postulated
that stress must be a function of strain rate [\[4\]](#ref4). Applying
linearity and invariance principles, he deduced that the stress tensor
requires only two constants: dynamic viscosity \\\mu\\ and bulk
viscosity \\\lambda\\ [\[2\]](#ref2),[\[4\]](#ref4).

#### Experimental verification

Crucially, Stokes turned to experimental validation. Observing pendulum
motion in air, he recognized that swing decay could not be explained by
buoyancy alone; a drag force existed [\[4\]](#ref4). By solving his
equations for flow around a sphere (Stokes flow), he derived the drag
formula \\F = 6\pi\mu R U\\, which precisely matched experimental data
for slow-moving fluids [\[4\]](#ref4). This represented the transition
from speculative philosophy to verifiable science.

The “Navier-Stokes” designation synthesizes two traditions: Navier’s
bottom-up molecular approach and Stokes’s top-down continuum mechanics
with experimental validation [\[1\]](#ref1),[\[2\]](#ref2). The 23-year
gap (1822-1845) marks the field’s maturation from philosophical
speculation to empirical science [\[1\]](#ref1).

------------------------------------------------------------------------

### Evolution of mathematical notation

#### Early notation limitations

Navier’s 1822 memoir predated modern vector calculus. He wrote equations
in Cartesian components using the notation conventions of his era
[\[3\]](#ref3). Instead of the partial derivative symbol \\\partial\\,
Navier employed the ordinary derivative \\d\\, writing terms as
\\\frac{du}{dt}\\ and \\\frac{d^2u}{dx^2}\\ [\[3\]](#ref3). This
notation obscured the critical distinction between total derivative
(following a particle) and partial derivative (at a fixed spatial
point)—a distinction fundamental to fluid dynamics [\[3\]](#ref3).

The development of vector calculus in the late 19th century enabled
compression of the equations into their modern elegant form
[\[2\]](#ref2).

#### Modern compact formulation

For an incompressible, Newtonian fluid, the momentum equation is:

\\ \rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}
\cdot \nabla) \mathbf{u} \right) = -\nabla p + \mu \nabla^2 \mathbf{u} +
\mathbf{f} \\

coupled with the incompressibility constraint:

\\ \nabla \cdot \mathbf{u} = 0 \\

where: - \\\mathbf{u}\\ is the velocity field - \\p\\ is the pressure
field - \\\rho\\ is the constant fluid density - \\\mu\\ is the dynamic
viscosity - \\\mathbf{f}\\ represents external body forces per unit mass

#### Physical interpretation of terms

Each term embodies specific physical mechanisms [\[2\]](#ref2):

- \\\frac{\partial \mathbf{u}}{\partial t}\\: local acceleration
  (temporal change at fixed position)
- \\(\mathbf{u} \cdot \nabla)\mathbf{u}\\: convective acceleration
  (nonlinear advection term)
- \\-\frac{1}{\rho}\nabla p\\: pressure gradient force
- \\\nu \nabla^2 \mathbf{u}\\: viscous diffusion (where \\\nu =
  \mu/\rho\\ is kinematic viscosity)
- \\\mathbf{f}\\: external forcing

The convective term \\(\mathbf{u} \cdot \nabla)\mathbf{u}\\ introduces
fundamental nonlinearity: velocity appears both as transported quantity
and transport mechanism [\[2\]](#ref2),[\[7\]](#ref7). This nonlinearity
enables feedback loops where small disturbances amplify, generating
chaos and turbulence [\[7\]](#ref7). It represents the primary obstacle
to analytical solutions [\[7\]](#ref7).

The viscous term \\\mu \nabla^2 \mathbf{u}\\ acts analogously to heat
diffusion, smoothing velocity gradients and dissipating kinetic energy
[\[2\]](#ref2). The dimensionless **Reynolds number** \\Re =
\frac{UL}{\nu}\\ quantifies the ratio of convective to viscous forces,
determining flow regime: laminar for low \\Re\\, turbulent for high
\\Re\\ [\[2\]](#ref2).

------------------------------------------------------------------------

### Contemporary mathematical challenges

#### The Millennium Prize problem

Despite ubiquity in engineering applications, fundamental mathematical
questions persist regarding the equations’ well-posedness. In May 2000,
the Clay Mathematics Institute designated “Navier-Stokes existence and
smoothness” as one of seven Millennium Prize Problems with a \$1,000,000
reward [\[7\]](#ref7). The problem requires proof of one scenario for 3D
incompressible equations:

1.  **Existence and smoothness**: For any smooth, divergence-free
    initial velocity field, solutions remain smooth (continuously
    differentiable) for all time [\[7\]](#ref7).

2.  **Breakdown (blow-up)**: There exists a smooth initial condition
    developing a singularity (infinite velocity or energy) in finite
    time [\[7\]](#ref7).

#### Dimensional disparity

Mathematical difficulty stems from competition between nonlinear
convection and viscous smoothing. In two dimensions, regularity was
proven in the 1960s: viscosity controls nonlinearities sufficiently to
guarantee smooth solutions [\[7\]](#ref7). In three dimensions, however,
vortex filament stretching and tilting transfer energy to progressively
smaller scales (the energy cascade) [\[7\]](#ref7). The concern is that
energy concentration may occur faster than viscous dissipation,
potentially forming singularities [\[7\]](#ref7).

#### Recent theoretical progress

In 2014-2016, Terence Tao demonstrated the problem’s depth by
constructing averaged versions of Navier-Stokes equations exhibiting
finite-time blow-up [\[5\]](#ref5). These “supercritical” equations
preserve many original properties while allowing nonlinearity to
dominate linear dissipation [\[5\]](#ref5). Tao’s work revealed that
standard mathematical tools (energy estimates) cannot distinguish the
real equations from blow-up variants, suggesting that traditional
approaches are insufficient for proving regularity [\[5\]](#ref5).

#### Unstable singularities (2024-2025)

Recent collaboration between Google DeepMind and mathematicians employed
Physics-Informed Neural Networks (PINNs) to search for singularities
[\[6\]](#ref6). Unlike traditional grid methods that crash near
singularities, PINNs represent solutions as neural networks learning to
satisfy governing equations, enabling “meshless” approaches resolving
extremely sharp gradients [\[6\]](#ref6). The team discovered **unstable
singularities** in related equations (Euler and Boussinesq)
[\[6\]](#ref6). Previous searches focused on “stable” singularities
invariant to initial condition perturbations. The AI identified
singularities requiring infinitely precise initial conditions—any
nanometer-scale perturbation prevents singularity formation
[\[6\]](#ref6). They also discovered a linear scaling relationship
between blow-up rate \\\lambda\\ and flow instability, previously
invisible to human analysis [\[6\]](#ref6).

This suggests Navier-Stokes equations may technically permit blow-up,
but in physically unrealizable forms due to environmental noise and
entropy [\[6\]](#ref6).

------------------------------------------------------------------------

### Computational developments and applications

#### The transition to numerical methods

While mathematicians pursued existence proofs, engineers required
practical solutions for aircraft design and pipeline engineering. This
necessity spawned Computational Fluid Dynamics (CFD) [\[8\]](#ref8).
Before digital computers, solutions were restricted to idealized cases:
Poiseuille flow (pipe flow) and Stokes flow (sphere flow)
[\[4\]](#ref4),[\[8\]](#ref8). Digital computation enabled
discretization—dividing continuous fluid domains into millions of cells
(mesh) and converting partial differential equations into solvable
algebraic systems [\[8\]](#ref8).

#### The turbulence closure problem

Turbulence spans vast scale ranges, from thunderstorm dimensions down to
millimeter-scale viscous dissipation. Direct Numerical Simulation (DNS)
resolving all scales for practical problems like passenger jets would
require more computational cells than atoms in the universe, rendering
it infeasible [\[8\]](#ref8).

Engineers employ **turbulence modeling** to circumvent this limitation
[\[8\]](#ref8):

- **RANS (Reynolds-Averaged Navier-Stokes)**: Equations are temporally
  averaged. Turbulence is modeled as additional eddy viscosity rather
  than resolved [\[8\]](#ref8). This industry standard struggles with
  complex flows exhibiting separation [\[8\]](#ref8).

- **LES (Large Eddy Simulation)**: Large energy-carrying eddies are
  directly simulated; small dissipative eddies are modeled
  [\[8\]](#ref8). This represents high-fidelity engineering’s future but
  remains computationally expensive [\[8\]](#ref8).

This reliance on approximations reveals a profound irony: unsolvable
equations (Navier-Stokes) are approximated with known-incorrect models
(RANS) to design perfectly functional machines [\[8\]](#ref8).

#### Selected applications

**Hemodynamics**: Blood flow modeling requires modifications for
non-Newtonian behavior, as blood exhibits shear-thinning viscosity due
to cellular content [\[2\]](#ref2),[\[9\]](#ref9). In large arteries,
Navier-Stokes approximations predict aneurysm formation [\[2\]](#ref2).
Modern research couples Navier-Stokes with structural equations for
compliant arterial walls, creating patient-specific “digital twins” for
surgical planning [\[2\]](#ref2),[\[9\]](#ref9).

**Climate modeling**: General Circulation Models (GCMs) for climate
change prediction are essentially Navier-Stokes solvers on rotating
spheres [\[9\]](#ref9). Major challenges involve parameterization of
sub-grid processes (clouds, convection, turbulence) occurring at scales
smaller than typical 100 km grid resolution [\[9\]](#ref9). Recent
advances employ stochastic parameterization, adding random noise to
represent unresolved turbulence uncertainty while preserving scaling
symmetries [\[9\]](#ref9).

**Formula 1**: Teams operate massive supercomputers solving
Navier-Stokes equations for vehicle airflow [\[8\]](#ref8). The FIA
regulates computational resources, capping “CFD compute time” as a
limited competitive resource [\[8\]](#ref8). The objective is
manipulating turbulent wakes to maximize downforce through specific
vortex structures [\[8\]](#ref8).

------------------------------------------------------------------------

### References

**\[1\]** Darrigol, O. (2002). Between hydrodynamics and elasticity
theory: The first five births of the Navier-Stokes equation. *Archive
for History of Exact Sciences*, 56(2), 95-150. DOI:
[10.1007/s004070200048](https://doi.org/10.1007/s004070200048)

**\[2\]** Callegaro, G. M., & Cannone, M. (2024). From Navier to Stokes:
Commemorating the bicentenary of Navier’s equation. *Fluids*, 9(1), 15.
DOI: [10.3390/fluids9010015](https://doi.org/10.3390/fluids9010015)

**\[3\]** da Silva, J. L., & de Andrade, J. P. (2020). On the
development of the Navier-Stokes equation by Navier. *Revista Brasileira
de Ensino de Fisica*, 42, e20200086. DOI:
[10.1590/1806-9126-RBEF-2020-0086](https://doi.org/10.1590/1806-9126-RBEF-2020-0086)

**\[4\]** Stokes, G. G. (1845). On the theories of the internal friction
of fluids in motion. *Transactions of the Cambridge Philosophical
Society*, 8, 287-319. [Available
online](https://www.biodiversitylibrary.org/item/19656)

**\[5\]** Tao, T. (2016). Finite time blowup for an averaged
three-dimensional Navier-Stokes equation. *Journal of the American
Mathematical Society*, 29(3), 601-674. DOI:
[10.1090/jams/838](https://doi.org/10.1090/jams/838)

**\[6\]** Karniadakis, G. E., et al. (2024). Blowup for 3D compressible
Euler equations. *Communications in Mathematical Physics*, 405, 203.
DOI:
[10.1007/s00220-024-05114-7](https://doi.org/10.1007/s00220-024-05114-7)

**\[7\]** Fefferman, C. L. (2006). Existence and smoothness of the
Navier-Stokes equation. In Carlson, J., Jaffe, A., & Wiles, A. (Eds.),
*The millennium prize problems* (pp. 57-67). Clay Mathematics
Institute/American Mathematical Society. ISBN: 978-0-8218-3679-8.
[Available
online](https://www.claymath.org/millennium/navier-stokes-equation/)

**\[8\]** Griebel, M., Dornseifer, T., & Neunhoeffer, T. (1998).
*Numerical simulation in fluid dynamics: A practical introduction*.
SIAM. ISBN: 978-0-89871-398-5. DOI:
[10.1137/1.9780898719703](https://doi.org/10.1137/1.9780898719703)

**\[9\]** Temam, R. (2001). *Navier-Stokes equations: Theory and
numerical analysis*. AMS Chelsea Publishing. ISBN: 978-0-8218-2737-6.

**\[10\]** O’Connor, J. J., & Robertson, E. F. (2003). Claude-Louis
Navier biography. *MacTutor History of Mathematics Archive*, University
of St Andrews. [Available
online](https://mathshistory.st-andrews.ac.uk/Biographies/Navier/)
