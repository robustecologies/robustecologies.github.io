# Enumerate all directed spanning trees of a complete digraph

Generates every directed spanning tree of the complete directed graph
\\K'\_m\\ on \\m\\ nodes, following the combinatorial enumeration of
Malekmohammadi & Mostafaee (2017). The total count is \\2^{m-1} \\
m^{m-2}\\, since by Cayley's formula \\K_m\\ has \\m^{m-2}\\ undirected
spanning trees and each of the \\m - 1\\ edges admits two orientations.

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

Cayley, A. (1889). A theorem on trees. *Quarterly Journal of Pure and
Applied Mathematics*, 23, 376-378.

Malekmohammadi, N., & Mostafaee, A. (2017). Obtaining all extreme rays
of a special cone using spanning trees in a complete digraph:
application in DEA. *Journal of the Operational Research Society*,
69(3), 465-472.
[doi:10.1057/s41274-017-0265-9](https://doi.org/10.1057/s41274-017-0265-9)

## See also

[`solid_angle_3d_from_rays`](https://robustecologies.github.io/SolidAngleR/reference/solid_angle_3d_from_rays.md)
for the consumer that converts enumerated extreme rays into a 3D solid
angle;
[`compute_solid_angle`](https://robustecologies.github.io/SolidAngleR/reference/compute_solid_angle.md)
for the dispatcher that delegates to the relevant backend.

## Examples

``` r
if (FALSE) { # \dontrun{
# m = 3: 2^(3-1) * 3^(3-2) = 12 trees
trees_m3 <- generate_spanning_trees(3)
length(trees_m3)
trees_m3[[1]]

# m = 4: 2^(4-1) * 4^(4-2) = 128 trees
trees_m4 <- generate_spanning_trees(4)
length(trees_m4)
} # }
```
