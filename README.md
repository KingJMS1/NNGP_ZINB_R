# ZINB.GP

`ZINB.GP` fits the Bayesian zero-inflated negative-binomial Gaussian-process
model described in [He and Huang (2024)](https://doi.org/10.1016/j.jspi.2023.106098).
The `ZINB_GP()` entry point selects the appropriate spatial, temporal, or
spatial-temporal model from the supplied design matrices and the two
`use_*_gp` switches.

## Installation

After its CRAN release, install the package with:

```r
install.packages("ZINB.GP")
```

Install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("KingJMS1/GP_ZINB_R")
```

## Complete example

This example fits the full model, with spatial and temporal Gaussian-process
(GP) effects in both the zero-inflation and count components. Rows of
`counts` are spatial locations and columns are time points. `make_y_Vs_Vt()`
turns this grid into the response and random-effect design matrices expected
by `ZINB_GP()`.

```r
library(ZINB.GP)

set.seed(1)

# Counts at four locations over three time points.
counts <- matrix(
    c(0, 1, 0,
      2, 0, 1,
      0, 3, 1,
      1, 0, 2),
    nrow = 4,
    byrow = TRUE
)

setup <- make_y_Vs_Vt(counts)
y <- setup$y
Vs <- setup$Vs
Vt <- setup$Vt

# One row per element of y: an intercept and an observation-level covariate.
X <- cbind(
    intercept = 1,
    covariate = scale(seq_along(y))[, 1]
)

# One coordinate pair per row of counts and one time coordinate per column.
coords <- rbind(
    c(0, 0), c(1, 0), c(0, 1), c(1, 1)
)
Ds <- as.matrix(dist(coords))
Dt <- as.matrix(dist(seq_len(ncol(counts))))

fit <- ZINB_GP(
    X = X,
    y = y,
    coords = coords,
    nsim = 1000,
    burn = 200,
    use_count_gp = TRUE,
    use_inflation_gp = TRUE,
    thin = 1,
    save_ypred = TRUE,
    print_progress = TRUE,
    print_iter = 100,
    Vs = Vs,
    Vt = Vt,
    Ds = Ds,
    Dt = Dt
)

# Posterior summaries for the fixed effects.
apply(fit$Alpha, 2, quantile, probs = c(0.025, 0.5, 0.975))
apply(fit$Beta, 2, quantile, probs = c(0.025, 0.5, 0.975))

# With save_ypred = TRUE, rows are saved MCMC iterations and columns are
# observations. These are posterior predictive count draws for this full model.
posterior_mean_count <- colMeans(fit$Y_pred)
posterior_at_risk_probability <- colMeans(fit$at_risk)
```

Use substantially longer chains than this small demonstration for analysis;
assess convergence and effective sample sizes before interpreting estimates.

## Inputs to `ZINB_GP()`

| Argument | Description |
| --- | --- |
| `X` | Numeric fixed-effect design matrix with one row per observation. Include an intercept column if required. |
| `y` | Non-negative integer response vector, with one entry per row of `X`. |
| `coords` | Spatial coordinates for the legacy spatial-temporal implementation. Supply one row per spatial location when using both spatial and temporal GPs. |
| `nsim` | Total MCMC iterations. It must be greater than `burn`; defaults to 5000. |
| `burn` | Number of initial MCMC iterations to discard; defaults to 1000. |
| `use_count_gp` | Logical: include GP random effects in the negative-binomial count component. Defaults to `TRUE`. |
| `use_inflation_gp` | Logical: include GP random effects in the zero-inflation component. Defaults to `FALSE`. At least one of the two GP switches must be `TRUE`. |
| `thin` | Save every `thin`-th post-burn-in iteration; defaults to 1. |
| `kern` | Optional kernel function with arguments `(distance_matrix, length_scale)`. The default is the package squared-exponential kernel. |
| `save_ypred` | Logical: retain posterior predictive draws and at-risk indicators. Defaults to `FALSE`. |
| `print_iter` | Report progress after this many iterations when `print_progress = TRUE`; defaults to 100. |
| `print_progress` | Logical: report MCMC progress with `message()` (suppress with `suppressMessages()`); defaults to `FALSE`. |
| `Vs` | Spatial random-effect design matrix with one row per observation. `make_y_Vs_Vt()` creates it from a location-by-time count matrix. Set to `NULL` when no spatial GP is wanted. |
| `Vt` | Temporal random-effect design matrix with one row per observation. `make_y_Vs_Vt()` creates it from a location-by-time count matrix. Set to `NULL` when no temporal GP is wanted. |
| `Ds` | Square spatial distance matrix, including the baseline location. Its diagonal must be zero and it must correspond to `Vs`. |
| `Dt` | Square temporal distance matrix, including the baseline time point. Its diagonal must be zero and it must correspond to `Vt`. |
| `ltPrior` | Optional list of temporal-length-scale controls: `list(max, mh_sd, a, b)`, where `a` and `b` parameterize the gamma prior and `mh_sd` is the proposal standard deviation. |
| `lsPrior` | Optional list of spatial-length-scale controls: `list(max, mh_sd, a, b)`. |
| `sigmaPrior` | Optional list of inverse-gamma GP-scale prior parameters: `list(a, b)`. |
| `noisePrior` | Optional list of GP noise-ratio controls: `list(a, b, mh_sd)`, where `a` and `b` parameterize the beta prior. |
| `mh_sd_r` | Optional Metropolis-Hastings proposal standard deviation for the negative-binomial dispersion parameter. |

`make_y_Vs_Vt()` uses the first spatial location and first time point as
baselines. Consequently, its `Vs` and `Vt` matrices have one fewer column than
the respective dimensions of `Ds` and `Dt`; do not remove the first row and
column from either distance matrix.

### Selecting the model

Provide `Vs`/`Ds` for spatial effects, `Vt`/`Dt` for temporal effects, or both
for a spatial-temporal model. The GP switches choose the components to which
those effects apply:

| `use_inflation_gp` | `use_count_gp` | Fitted GP effects |
| --- | --- | --- |
| `TRUE` | `TRUE` | Both zero-inflation and count components. |
| `TRUE` | `FALSE` | Zero-inflation component only. |
| `FALSE` | `TRUE` | Count component only. |
| `FALSE` | `FALSE` | Not supported; use a non-GP ZINB model instead. |

## Output

`ZINB_GP()` returns a named list of saved posterior draws. The exact elements
depend on the supplied spatial/temporal matrices and on the selected GP
components. With both spatial and temporal GPs enabled in both components,
the list contains the following elements.

| Element | Description |
| --- | --- |
| `Alpha` | Matrix of fixed-effect draws for the zero-inflation component. |
| `Beta` | Matrix of fixed-effect draws for the negative-binomial count component. |
| `A`, `B` | Matrices of spatial (`A`) and temporal (`B`) random-effect draws for the zero-inflation component. |
| `C`, `D` | Matrices of spatial (`C`) and temporal (`D`) random-effect draws for the count component. |
| `L1s`, `L2s` | Spatial GP length-scale draws for the zero-inflation and count components, respectively. |
| `Sigma1s`, `Sigma2s` | Spatial GP scale draws for the zero-inflation and count components. |
| `Noise1s`, `Noise2s` | Spatial GP kernel-to-noise mixing-ratio draws for the zero-inflation and count components. |
| `L1t`, `L2t` | Temporal GP length-scale draws for the zero-inflation and count components. |
| `Sigma1t`, `Sigma2t` | Temporal GP scale draws for the zero-inflation and count components. |
| `Noise1t`, `Noise2t` | Temporal GP kernel-to-noise mixing-ratio draws for the zero-inflation and count components. |
| `R` | Vector of negative-binomial dispersion-parameter draws. |
| `Y_pred` | Matrix of posterior predictive count draws, included only when `save_ypred = TRUE`. |
| `at_risk` | Matrix of sampled at-risk indicators, included only when `save_ypred = TRUE`. |

When a GP component or spatial/temporal dimension is omitted, the corresponding
random-effect and GP-parameter elements are omitted from the result. Check
`names(fit)` to see the elements returned for a particular configuration.

## Further information

The package reference is available at
<https://kingjms1.github.io/GP_ZINB_R/reference/index.html>.
