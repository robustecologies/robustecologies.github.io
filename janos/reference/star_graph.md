# Create a star graph adjacency matrix

Constructs the adjacency matrix for a star graph where a central node
(node 1) is connected to all other nodes, and peripheral nodes have no
mutual connections.

## Usage

``` r
star_graph(n, weights = 1)
```

## Arguments

- n:

  Integer, total number of nodes (including the central node). Must be
  at least 2.

- weights:

  Numeric scalar, edge weight. Default: 1.

## Value

A numeric matrix of class `"adjacency_matrix"`.

## See also

[`lattice_graph`](https://robustecologies.github.io/janos/reference/lattice_graph.md),
[`complete_graph`](https://robustecologies.github.io/janos/reference/complete_graph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
adj <- star_graph(10)
} # }
```
