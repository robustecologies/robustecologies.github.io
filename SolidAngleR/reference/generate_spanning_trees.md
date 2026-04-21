# Generate all spanning trees of a complete digraph

Implements the algorithm by Malekmohammadi & Mostafaee (2017) to
systematically generate all directed spanning trees of a complete
directed graph K'\_m with `m` nodes. This is equivalent to finding all
potential vertices of a specific polyhedral cone defined by ratio
constraints.

## Usage

``` r
generate_spanning_trees(m)
```

## Arguments

- m:

  An integer (`>= 2`) representing the number of nodes in the graph.

## Value

A list where each element is a `2 x (m-1)` integer matrix. Each matrix
represents a unique directed spanning tree. In each matrix, a column
`c(u, v)` represents a directed edge from node `u` to node `v`.

## Details

The algorithm works in two main stages. First, it generates all
`m^(m-2)` unique undirected tree structures using their Prüfer sequence
representations. Second, for each structure, it applies all `2^(m-1)`
possible directed orientations to the `m-1` edges. The total number of
generated trees is `2^(m-1) * m^(m-2)`. This method is purely
combinatorial and avoids solving any linear programming problems.

## References

Malekmohammadi, N., & Mostafaee, A. (2017). Obtaining all extreme rays
of a special cone using spanning trees in a complete digraph:
application in DEA. *Journal of the Operational Research Society*,
69(3), 465-472.
[doi:10.1057/s41274-017-0265-9](https://doi.org/10.1057/s41274-017-0265-9)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Basic case for m=3
# Expected number of trees = 2^(3-1) * 3^(3-2) = 4 * 3 = 12
trees_m3 <- generate_spanning_trees(3)
cat("Total trees generated for m=3:", length(trees_m3), "\n")

# Print the first generated tree
print("First tree for m=3:")
print(trees_m3[[1]])

# Example 2: More complex case for m=4
# Expected number of trees = 2^(4-1) * 4^(4-2) = 8 * 16 = 128
trees_m4 <- generate_spanning_trees(4)
cat("\nTotal trees generated for m=4:", length(trees_m4), "\n")

# A spanning tree for a 4-node graph has 3 edges
# Inspect the structure of the 100th tree
print("100th tree for m=4:")
print(trees_m4[[100]])
} # }
```
