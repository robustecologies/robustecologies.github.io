# Changelog

## SolidAngleR 0.1.0

### Major features

#### Implemented methods

- **Generic n-dimensional hypergeometric series**
  - Full implementation for arbitrary dimensions (n \> 3)
  - Recursive generation of weak compositions for series exponents
  - Enables analytical computation for high-dimensional
    positive-definite cones
  - Integrated into
    [`hypergeometric_series()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series.md)
- **C++ Optimization (Rcpp)**
  - Accelerated
    [`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md)
    using RcppArmadillo
  - O(n) Givens rotations implemented in C++
  - Significant performance boost for Monte Carlo simulations
  - Automatic dispatch for inverse transform method
- **General 3D Cone Visualization**
  - New function
    [`plot_cone_3d()`](https://robustecologies.github.io/SolidAngleR/reference/plot_cone_3d.md)
    using plotly
  - Visualizes arbitrary simplicial cones (rays, faces, spherical
    triangle)
  - Interactive exploration of cone geometry
- **Hypergeometric series** (Ribando 2006; Fitisone & Zhou 2023)
  - Full implementation for n=2, n=3 dimensions with closed-form
    solutions
  - Automatic fallback to decomposition for n≥4 (now supports direct
    series for PD cones)
  - Convergence checking with positive-definite criterion
- **Tridiagonal series optimization** (Fitisone & Zhou 2023)
  - Specialized series for cones with tridiagonal Gram matrices
  - Exponential speedup: O(N·n) instead of O(N^(n(n-1)/2))
  - Automatic detection and application
- **Decomposition methods** (Fitisone & Zhou 2023)
  - Brion-Vergne decomposition for arbitrary cones
  - Ensures universal applicability
  - Handles non-positive-definite cases
- **Geometric methods** (Mazonka 2012)
  - Right circular cone formula
  - Polyhedral cone solid angles (Van Oosterom-Strackee formula)
  - Cone segment calculations
  - Intersecting cones computation
- **Multivariate normal integration** (Genz & Bretz)
  - Feasibility domain analysis for ecological interaction matrices
  - Monte Carlo estimation fallback
  - Both normalized and scaled solid angle measures

#### Uniform sampling on n-dimensional spheres and cones

- **O(n) cone sampling algorithm** (Arun & Venkatapathi, 2025)
  - Generates uniform random vectors directly within n-dimensional cones
  - Eliminates exponential cost of rejection sampling in high dimensions
  - Two methods: inverse transform (default) and one-dimensional
    rejection (for numerical stability)
  - Functions:
    [`generate_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_sample.md),
    [`generate_cone_samples()`](https://robustecologies.github.io/SolidAngleR/reference/generate_cone_samples.md),
    [`generate_hollow_cone_sample()`](https://robustecologies.github.io/SolidAngleR/reference/generate_hollow_cone_sample.md)
  - Supporting functions:
    [`theta_to_omega()`](https://robustecologies.github.io/SolidAngleR/reference/theta_to_omega.md),
    [`omega_to_theta()`](https://robustecologies.github.io/SolidAngleR/reference/omega_to_theta.md),
    [`rotate_from_canonical()`](https://robustecologies.github.io/SolidAngleR/reference/rotate_from_canonical.md)
  - Validation function:
    [`verify_cone_uniformity()`](https://robustecologies.github.io/SolidAngleR/reference/verify_cone_uniformity.md)
    with KS and chi-squared tests
- **Comprehensive sampling toolkit**
  - Full sphere sampling:
    [`generate_point_on_sphere()`](https://robustecologies.github.io/SolidAngleR/reference/generate_point_on_sphere.md)
    (Box-Muller/Marsaglia method)
  - Hollow cone sampling for annular regions
  - Arbitrary orientations via efficient O(n) Givens rotations
  - Planar angle generation:
    [`generate_planar_angle_inverse()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_inverse.md),
    [`generate_planar_angle_rejection()`](https://robustecologies.github.io/SolidAngleR/reference/generate_planar_angle_rejection.md)
  - Cost analysis:
    [`rejection_cost()`](https://robustecologies.github.io/SolidAngleR/reference/rejection_cost.md)
    for comparing with naive rejection
- **New vignette: uniform sphere sampling**
  (`vignettes/5.uniform-sphere-sampling.Rmd`)
  - Historical development of sphere sampling techniques (1958-2025)
  - Mathematical foundations with incomplete beta function theory
  - Rigorous statistical validation (KS tests, chi-squared tests)
  - Performance analysis and comparison with rejection sampling
  - Practical applications: Monte Carlo integration, directional
    statistics, ray tracing
  - Comprehensive examples and visualizations
- **Relationship to existing Monte Carlo methods**
  - NEW cone sampling: Generates samples uniformly WITHIN a cone (no
    rejection)
  - Existing
    [`solid_angle_monte_carlo()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_monte_carlo.md):
    Estimates solid angles via rejection sampling
  - Methods are COMPLEMENTARY

#### Main interface functions

- [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md) -
  Main function with automatic method selection
- [`compute_solid_angles()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angles.md) -
  Batch computation for multiple cones
- [`diagnose_cone()`](https://robustecologies.github.io/SolidAngleR/reference/diagnose_cone.md) -
  Diagnostic tool to suggest optimal method
- [`omega()`](https://robustecologies.github.io/SolidAngleR/reference/omega.md) -
  High-level interface for feasibility domains

#### Utility functions

- Matrix operations:
  [`compute_associated_matrix()`](https://robustecologies.github.io/SolidAngleR/reference/compute_associated_matrix.md),
  [`is_positive_definite()`](https://robustecologies.github.io/SolidAngleR/reference/is_positive_definite.md),
  [`is_tridiagonal()`](https://robustecologies.github.io/SolidAngleR/reference/is_tridiagonal.md)
- Vector operations:
  [`normalize_vectors()`](https://robustecologies.github.io/SolidAngleR/reference/normalize_vectors.md),
  [`angle_between()`](https://robustecologies.github.io/SolidAngleR/reference/angle_between.md),
  [`cross_product_3d()`](https://robustecologies.github.io/SolidAngleR/reference/cross_product_3d.md)
- Geometric helpers:
  [`create_tridiagonal_cone()`](https://robustecologies.github.io/SolidAngleR/reference/create_tridiagonal_cone.md),
  [`spherical_triangle_area()`](https://robustecologies.github.io/SolidAngleR/reference/spherical_triangle_area.md),
  [`lhuilier_angle()`](https://robustecologies.github.io/SolidAngleR/reference/lhuilier_angle.md)

#### Documentation

- Comprehensive package documentation with mathematical theory
- Five detailed vignettes:
  - `1.theoretical-analysis`: Theoretical foundations and validation
  - `2.fitisone-zhou-methods`: Decomposition and series methods
  - `3.mazonka-geometric-methods`: Geometric approaches with interactive
    visualizations
  - `4.monte-carlo-methods`: Monte Carlo integration and error analysis
  - `5.uniform-sphere-sampling`: Comprehensive treatment of uniform
    sampling on spheres and cones
- Complete function documentation with equations, examples, and
  references

#### Bug fixes (from pre-release versions)

- **Critical gamma function bug**: Fixed operator precedence error in
  all hypergeometric series implementations
  - Changed `gamma(1 + x/2)` to `gamma((1+x)/2)` throughout
  - Affected:
    [`hypergeometric_series_2d()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_2d.md),
    [`hypergeometric_series_3d()`](https://robustecologies.github.io/SolidAngleR/reference/hypergeometric_series_3d.md),
    [`tridiagonal_series()`](https://robustecologies.github.io/SolidAngleR/reference/tridiagonal_series.md)
  - This bug caused ~10x error in solid angle computations
- **Geometric function normalization**: All geometric methods now return
  normalized values \[0,1\]
  - [`solid_angle_cone()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone.md),
    [`solid_angle_polyhedral()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_polyhedral.md),
    [`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md),
    [`solid_angle_intersecting_cones()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_intersecting_cones.md)
  - Previously returned steradians, now consistent with series methods
- **Names attribute cleanup**: Removed unwanted names from return values
  - Fixed in
    [`compute_solid_angle()`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
    and
    [`solid_angle_decomposition()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_decomposition.md)
- **Vignette compilation**: Fixed all rendering errors
  - Stack overflow prevention (reduced dimension in examples)
  - Polyhedral vertex format (rbind vs cbind)
  - Function renaming (Omega → omega)

### Known limitations

#### For version 0.2.0

- **Decomposition recursion**: Some random cone configurations trigger
  stack overflow
  - Occurs with certain non-positive-definite cone combinations
  - Workaround: Use `method = "series"` for known positive-definite
    cones
  - Will be fixed with improved termination conditions
- **Cone segment relationship**: Mathematical formula for segment/full
  cone relationship needs verification
  - [`solid_angle_cone_segment()`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_cone_segment.md)
    computes correctly
  - Test relationship `omega_segment = omega_full * gamma/(2π)` may be
    oversimplified
  - Requires deeper geometric analysis
- **Tridiagonal convergence**: Some non-orthogonal tridiagonal cones
  converge slowly
  - May require `max_terms > 2000` for certain configurations
  - Results are valid but convergence criterion may not be met
  - Performance optimization planned

#### Future enhancements

- Optimize decomposition algorithm to prevent stack overflow
- Add parallelization for batch computations
- Extend to non-Euclidean geometries (hyperbolic, spherical)
- Add visualization tools for higher dimensions

### Test suite

- **Comprehensive test coverage** with all tests passing
- NEW: 40+ tests for cone sampling functions
  - Basic functionality (unit vectors, dimensions, input validation)
  - Theta-omega mapping and inverse relationships
  - Planar angle generation (inverse and rejection methods)
  - Givens rotation (norm preservation, alignment, angle preservation)
  - Cone sample generation (unit vectors, cone constraints, arbitrary
    orientations)
  - Hollow cone sampling
  - Statistical uniformity (KS tests, chi-squared tests)
  - Rejection cost analysis
  - Edge cases (narrow cones, wide cones, high dimensions, 2D case)
  - Consistency with analytical formulas
- Existing tests for all other major functionality:
  - All analytical formulas (2D, 3D, circular cones)
  - Hypergeometric series convergence
  - Tridiagonal optimization
  - Decomposition methods
  - Matrix operations and validations
  - Batch processing
  - Edge cases and degeneracies

### Dependencies

- Base R (≥ 4.1.0)
- Matrix operations: standard `stats` package
- Multivariate normal integration: `mvtnorm` package
- Suggested: `knitr`, `rmarkdown`, `testthat`, `plotly`, `viridisLite`

### References

#### Core solid angle methods

- Ribando, J. M. (2006). “Measuring Solid Angles Beyond Dimension
  Three.” *Discrete & Computational Geometry*, 36(3):479-487.
- Fitisone, A., & Zhou, Y. (2023). “Solid angle measure of polyhedral
  cones.” arXiv:2304.11102 \[math.CO\].
- Mazonka, O. (2012). “Solid Angle of Conical Surfaces, Polyhedral
  Cones, and Intersecting Spherical Caps.” arXiv:1205.1396 \[math.MG\].
- Van Oosterom, A. & Strackee, J. (1983). “The Solid Angle of a Plane
  Triangle.” *IEEE Trans. Biomed. Eng.*, BME-30(2):125-126.

#### Uniform sphere sampling

- Arun, I., & Venkatapathi, M. (2025). “An O(n) algorithm for generating
  uniform random vectors in n-dimensional cones.” *Sankhyā: The Indian
  Journal of Statistics*, 87-A(2), 327-348. DOI:
  10.1007/s13171-025-00387-9
- Box, G.E.P. & Muller, M.E. (1958). “A note on the generation of random
  normal deviates.” *Annals of Mathematical Statistics*, 29(2), 610-611.
- Marsaglia, G. (1972). “Choosing a point from the surface of a sphere.”
  *Annals of Mathematical Statistics*, 43(2), 645-646.
- Arvo, J. (2001). “Stratified sampling of spherical triangles.”
  *Proceedings of ACM SIGGRAPH 2001*, 437-438.
