# Extract per-state C++ RHS expressions from a `model_spec`

Translate the symbolic right-hand sides of a formula-based `model_spec`
into character strings of valid C++ code, one per state variable, ready
to be inlined by external compilers. Designed as a stable
interoperability hook for downstream packages (e.g. `wadaR`) that
operate on the same dynamics under a different code-generation template.

## Usage

``` r
model_spec_rhs_cpp(model)
```

## Arguments

- model:

  A `model_spec` object built from a formula list, either as a
  deterministic ODE (argument `rhs` of
  [`model_spec()`](https://robustecologies.github.io/janos/reference/model_spec.md))
  or as a discrete map (argument `map`). Function-mode `rhs` is not
  supported because the symbolic structure is unavailable.

## Value

A list of class `model_spec_rhs_cpp` with the following entries:

- rhs:

  Character vector of length `n_states`. Each element is a C++
  expression for the time derivative (ODE) or successor (map) of the
  corresponding state variable, using `y[i]` for state references
  (0-based) and `parms[i]` for parameter references (0-based).

- target:

  Character scalar, `"deriv"` for ODE models and `"next_state"` for
  discrete-map models. Identifies the convention an external generator
  should use when assembling assignment statements.

- state_names:

  Character vector of state names, copied from `model$state_names`.

- parm_names:

  Character vector of parameter names, ordered to match the `parms[i]`
  indexing used in the `rhs` strings.

- parms:

  Named list of default parameter values, copied from `model$parms`.

- n_states:

  Integer, copied from `model$n_states`.

- is_map:

  Logical scalar, `TRUE` for discrete-map systems.

- meta:

  Copy of `model$meta`.

## Details

The function consumes the parsed formulas stored on the `model_spec`
object: `rhs_formulas` for ODEs, `map_formulas` for discrete maps. It
constructs 0-based index maps over `state_names` and `names(parms)`,
iterates `parse_formula_rhs()` and `expr_to_cpp()` over each formula,
and returns the resulting character vector without performing any
compilation or template assembly. The output is identical, in identifier
convention, to the strings that janos's internal compilers inline into
RK4/RK45 templates: state references appear as `y[i]`, parameter
references as `parms[i]`, time as `t`, and the constant pi as `M_PI`.

Only deterministic ODE and discrete-map families are accepted.
Stochastic differential equations, delay differential equations,
jump-diffusion processes, partial differential equations, piecewise
deterministic Markov processes, reaction-diffusion master equations and
continuous-time Markov chains are rejected because their right-hand
sides do not reduce to a single character vector of derivative
expressions and the corresponding inlining contracts are not part of
this stable interface. Models constructed with `positive_states` set to
a non-trivial subset are also rejected because positivity clamps are
implemented as guarded statements inside the integrator loop and are not
transferable through a per-state expression vector.

## References

Almaraz, P. (2024). janos: a general-purpose R framework for the
specification, simulation and analysis of dynamical systems.
<https://github.com/robustecologies/janos>.

## See also

[`model_spec`](https://robustecologies.github.io/janos/reference/model_spec.md),
[`dyn_sim`](https://robustecologies.github.io/janos/reference/dyn_sim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
duffing <- model_spec(
    rhs = list(
        x ~ y,
        y ~ -delta * y - x * (x * x - 1) + F * sin(omega * t)
    ),
    state_names = c("x", "y"),
    parms = list(delta = 0.15, F = 0.20, omega = 1.0)
)
rhs <- model_spec_rhs_cpp(duffing)
rhs$rhs
#> [1] "y[1]"
#> [2] "((((-parms[0]) * y[1]) - (y[0] * ((y[0] * y[0]) - 1.0))) +
#>      (parms[1] * std::sin((parms[2] * t))))"
rhs$target
#> [1] "deriv"
} # }
```
