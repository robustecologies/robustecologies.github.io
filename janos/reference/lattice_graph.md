# Create a 2D rectangular lattice adjacency matrix

Constructs the adjacency matrix for a rectangular lattice graph with
optional periodic boundary conditions. The resulting matrix encodes
nearest-neighbor connectivity on an `nx` by `ny` grid, where each
interior node connects to its four cardinal neighbors. Edge weights
default to 1.

## Usage

``` r
lattice_graph(nx, ny = NULL, bc = c("none", "periodic"), weights = 1)
```

## Arguments

- nx:

  Integer, number of nodes in the x-direction.

- ny:

  Integer or NULL. Number of nodes in the y-direction. If NULL, a 1D
  chain/ring of `nx` nodes is created.

- bc:

  Character, boundary condition type: `"none"` (open boundaries) or
  `"periodic"` (wraparound). Default: `"none"`.

- weights:

  Numeric scalar, edge weight for all connections. Default: 1.

## Value

A numeric matrix of class `"adjacency_matrix"` with a `layout`
attribute.

## Details

Nodes are numbered in row-major order: node `(i,j)` maps to index
`(j-1)*nx + i` for `i` in `1:nx` and `j` in `1:ny`. With periodic
boundary conditions, the grid wraps around in both dimensions forming a
torus. The returned matrix carries a `layout` attribute with
`type = "lattice"`, `nx`, and `ny` fields to enable automatic 2D heatmap
plotting.

When `ny` is `NULL`, the function returns a 1D chain or ring of `nx`
nodes.

## See also

[`ring_graph`](https://robustecologies.github.io/janos/reference/ring_graph.md),
[`complete_graph`](https://robustecologies.github.io/janos/reference/complete_graph.md),
[`star_graph`](https://robustecologies.github.io/janos/reference/star_graph.md),
[`random_graph`](https://robustecologies.github.io/janos/reference/random_graph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 5x5 periodic lattice
adj <- lattice_graph(5, 5, bc = "periodic")
print(dim(adj))

# 1D chain of 10 nodes
adj1d <- lattice_graph(10)
} # }
```
