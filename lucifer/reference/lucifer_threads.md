# Configure parallelization thread budget

Sets or queries the thread budget controlling how lucifer distributes
computational resources between inter-chain parallelism (multiple MCMC
chains running simultaneously) and intra-chain parallelism (OpenMP
threads for vectorized distribution math, diagnostics, and PSIS within
each chain).

## Usage

``` r
lucifer_threads(total = NULL, omp = NULL, chains = NULL, strategy = NULL)
```

## Arguments

- total:

  Integer total CPUs to use. If `NULL` (default), detects automatically
  via `parallel::detectCores(logical = FALSE)`.

- omp:

  Integer OpenMP threads per process. Overrides the automatic allocation
  from `strategy`. Set to 1 to disable OpenMP threading within chains.

- chains:

  Integer number of parallel chain workers. Overrides the automatic
  allocation from `strategy`.

- strategy:

  Character string controlling the allocation strategy:

  `"auto"`

  :   Distributes threads evenly between chains and OpenMP. On 8 CPUs
      with 4 chains, each chain gets 2 OpenMP threads. This is the
      recommended default.

  `"chains_only"`

  :   All CPUs to chain parallelism, each chain runs single-threaded
      (legacy behavior). Safest option for systems with threaded BLAS.

  `"omp_only"`

  :   All CPUs to OpenMP within a single chain process. Chains run
      sequentially. Best when running a single chain with large data.

  `"balanced"`

  :   Same as `"auto"`.

## Value

A named list with components `total_cpus`, `omp_threads`,
`chain_workers`, and `strategy` (invisibly when setting, visibly when
querying).

## Details

When called with no arguments, returns the current settings as a named
list without modifying anything. When arguments are provided, updates
the global thread budget and applies the OpenMP thread count immediately
via `omp_set_num_threads()`.

The budget is respected by
[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`Kfold`](https://robustecologies.github.io/lucifer/reference/Kfold.md),
[`SBI`](https://robustecologies.github.io/lucifer/reference/SBI.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md)
(Pathfinder), and
[`ABC`](https://robustecologies.github.io/lucifer/reference/ABC.md) when
using parallel dispatch.

The setting respects the `_R_CHECK_LIMIT_CORES_` environment variable
(used by CRAN) and `OMP_THREAD_LIMIT` when set.

## See also

[`lucifer`](https://robustecologies.github.io/lucifer/reference/lucifer.md),
[`LaplaceApproximation`](https://robustecologies.github.io/lucifer/reference/LaplaceApproximation.md),
[`VariationalBayes`](https://robustecologies.github.io/lucifer/reference/VariationalBayes.md).

## Examples

``` r
# Query current settings
lucifer_threads()
#> $total_cpus
#> [1] 20
#> 
#> $omp_threads
#> [1] 20
#> 
#> $chain_workers
#> [1] 1
#> 
#> $strategy
#> [1] "auto"
#> 

if (FALSE) { # \dontrun{
# Use 8 CPUs with auto strategy (default)
lucifer_threads(total = 8)

# Force all threads to chains, no OpenMP (legacy)
lucifer_threads(total = 8, strategy = "chains_only")

# Manual override: 4 chain workers, 2 OpenMP threads each
lucifer_threads(omp = 2, chains = 4)
} # }
```
