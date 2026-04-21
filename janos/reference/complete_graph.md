# Create a complete graph adjacency matrix

Constructs the adjacency matrix for a fully connected graph where every
node is connected to every other node with equal weight.

## Usage

``` r
complete_graph(n, weights = 1)
```

## Arguments

- n:

  Integer, number of nodes. Must be at least 2.

- weights:

  Numeric scalar, edge weight. Default: 1.

## Value

A numeric matrix of class `"adjacency_matrix"`.

## See also

[`lattice_graph`](https://robustecologies.github.io/janos/reference/lattice_graph.md),
[`star_graph`](https://robustecologies.github.io/janos/reference/star_graph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
adj <- complete_graph(5)
} # }
```
