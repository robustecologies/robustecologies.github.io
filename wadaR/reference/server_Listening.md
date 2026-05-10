# Server listening

This function is not intended to be called directly by the user. It is
an internal-only function to prevent cluster problems while using the
`INCA` algorithm in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Usage

``` r
server_Listening(n = 2, port = 19009)
```

## Arguments

- n:

  integer; the number of CPUs. For more information, see
  [`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

- port:

  integer; the port for server listening, defaults to `19009`.

## Value

See Details.

## Details

For the `INCA` algorithm, a server has been built into
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).
The server exchanges information between processes and has been designed
to be portable. The `server_Listening` function is run as a separate
process via the `system` function when `INCA` is selected in
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

Socket connections and the `serialize` function are used as per the Snow
package to update a single proposal covariance matrix given all parallel
chains. The sockets are opened/closed in each process with a small
random sleep time to avoid collisions during connections to the internal
server. Blocking sockets are used to synchronize processes.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md).

## Author

Silvere Vialet-Chabrand

## Examples

``` r
if (FALSE) { # \dontrun{
## see package vignettes for full examples involving server_Listening
} # }
```
