# Erdos-Renyi random graph adjacency matrix

Constructs the adjacency matrix for a G(n,p) Erdos-Renyi random graph
where each possible edge is included independently with probability `p`.

## Usage

``` r
random_graph(n, p = 0.5, seed = NULL, directed = FALSE, weights = 1)
```

## Arguments

- n:

  Integer, number of nodes. Must be at least 2.

- p:

  Numeric, edge probability in `[0,1]`. Default: 0.5.

- seed:

  Integer, random seed for reproducibility. Default: NULL.

- directed:

  Logical, whether the graph is directed. Default: FALSE.

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
adj <- random_graph(20, p = 0.3, seed = 42)
} # }
```
