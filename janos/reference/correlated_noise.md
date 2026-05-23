# Correlated Gaussian noise specification

Creates a noise specification for driving an SDE with correlated Wiener
increments. Given an n x n positive-definite covariance matrix Sigma,
the Cholesky factor L = chol(Sigma) is computed so that dW = L \* dZ
where dZ is a vector of independent standard normals.

## Usage

``` r
correlated_noise(Sigma)
```

## Arguments

- Sigma:

  Positive-definite covariance matrix. Its dimension must match the
  number of state variables in the model.

## Value

An S3 object of class `noise_spec` with `type = "correlated"`.

## Details

The standard SDE formulation dX_i = f_i(X)dt + g_i(X)dW_i assumes
independent Wiener processes for each state variable. When the noise
sources are correlated (e.g., shared environmental fluctuations acting
on multiple species), the Wiener increments satisfy E\[dW_i dW_j\] =
Sigma(i,j) dt. The Cholesky factorization Sigma = L L^T allows sampling
correlated increments efficiently as dW = L \* dZ where dZ_i are
independent N(0,1) variates.

Correlated noise can be composed with Levy noise: the Cholesky rotation
is applied to the alpha-stable increments instead of Gaussian
increments. It cannot be composed with fBm or colored noise.

## See also

[`levy_noise`](https://robustecologies.github.io/janos/reference/levy_noise.md),
[`fbm_noise`](https://robustecologies.github.io/janos/reference/fbm_noise.md),
[`colored_noise`](https://robustecologies.github.io/janos/reference/colored_noise.md),
[`system_spec`](https://robustecologies.github.io/janos/reference/system_spec.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
model <- system_spec(
    rhs = list(x ~ -a * x, y ~ -b * y),
    diffusion = list(x ~ sigma, y ~ sigma),
    noise = correlated_noise(Sigma),
    state_names = c("x", "y"),
    parms = list(a = 1, b = 2, sigma = 0.5),
    init = c(x = 1, y = 1)
)
result <- dyn_sim(model, t_max = 10, solver = solver_euler_maruyama())
print(result)
plot(result)
} # }
```
