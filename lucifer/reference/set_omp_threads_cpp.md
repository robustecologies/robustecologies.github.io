# Set the number of OpenMP threads at runtime

Sets the number of threads used by OpenMP parallel regions in the C++
backend. Called internally by
[`lucifer_threads`](https://robustecologies.github.io/lucifer/reference/lucifer_threads.md).
Has no effect if compiled without OpenMP support.

## Usage

``` r
set_omp_threads_cpp(n)
```

## Arguments

- n:

  Integer number of threads to use.

## Value

Invisibly returns `NULL`.
