# Create a ring (cycle) graph adjacency matrix

Convenience wrapper around
[`lattice_graph`](https://robustecologies.github.io/janos/reference/lattice_graph.md)
for 1D ring or chain topologies.

## Usage

``` r
ring_graph(n, bc = c("periodic", "none"), weights = 1)
```

## Arguments

- n:

  Integer, number of nodes.

- bc:

  Character, `"periodic"` for a ring or `"none"` for an open chain.
  Default: `"periodic"`.

- weights:

  Numeric scalar, edge weight. Default: 1.

## Value

A numeric matrix of class `"adjacency_matrix"`.

## See also

[`lattice_graph`](https://robustecologies.github.io/janos/reference/lattice_graph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
adj <- ring_graph(20, bc = "periodic")
} # }
```
