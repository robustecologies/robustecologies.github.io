# Decompose simplicial cone

Implements the decomposition method from Theorem 3.3 (second
decomposition) of Fitisone, A., & Zhou (2023), which decomposes a
simplicial cone into cones that either have dimension less than n or
have positive definite associated matrices.

## Usage

``` r
decompose_cone(V, method = "line")
```

## Arguments

- V:

  An n x n matrix where columns are linearly independent unit vectors

- method:

  Method to use: "line" for Theorem 3.3 (default)

## Value

A list containing:

- `cones`: List of matrices, each representing a cone

- `signs`: Vector of signs (+1 or -1) for each cone

- `n_cones`: Number of cones in decomposition

- `all_pd`: Logical indicating if all cones are PD

## Details

By Theorem 3.3, any simplicial cone can be decomposed into finitely many
cones (at most (n-1)! by Corollary 3.4) such that each cone either: (I)
has affine dimension \< n (contributing 0 to solid angle), or (II) has
positive definite associated matrix (Theorem 1.5 applies)

The solid angle of the original cone is: \\\Omega_n(C) = \Sigma s_i
\cdot \Omega_n(C_i)\\

## References

Fitisone, A., & Zhou, Y. (2023). Solid angle measure of polyhedral
cones. arXiv:2304.11102 (math.CO). <https://arxiv.org/abs/2304.11102>

## Examples

``` r
# Decompose a 3D cone
v1 <- c(1, 0, 0)
v2 <- c(-0.5, 0.8, 0.3)
v3 <- c(0.2, 0.3, 0.9)
V <- cbind(v1, v2, v3)
V <- normalize_vectors(V)

decomp <- decompose_cone(V)
print(decomp$n_cones)  # Number of cones
#> [1] 1
print(decomp$all_pd)   # All positive definite?
#> [1] TRUE
```
