# Hexbin density of compilation sites on a Robinson world map

Hexbin density of compilation sites on a Robinson world map

## Usage

``` r
vt_plot_world_hex(x)
```

## Arguments

- x:

  A `vt_compilation`.

## Value

A `ggplot`. Falls back to a Cartesian rectangle when
`sf`/`rnaturalearth` are missing.

## References

Pebesma (2018).

## See also

[`vt_plot_world_map()`](https://robustecologies.github.io/VerteTIME/reference/vt_plot_world_map.md)

## Examples

``` r
if (FALSE)  vt_plot_world_hex(vt_ingest_all())  # \dontrun{}
```
