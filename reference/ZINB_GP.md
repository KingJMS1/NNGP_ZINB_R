# ZINB_GP

Fits the zero-inflated negative-binomial Gaussian-process model
described in
[doi:10.1016/j.jspi.2023.106098](https://doi.org/10.1016/j.jspi.2023.106098)
. The supplied spatial and temporal design and distance matrices
determine which Gaussian processes are included.

## Usage

``` r
ZINB_GP(
  X,
  y,
  coords,
  nsim = 5000,
  burn = 1000,
  use_count_gp = TRUE,
  use_inflation_gp = FALSE,
  thin = 1,
  kern = NULL,
  save_ypred = FALSE,
  print_iter = 100,
  print_progress = FALSE,
  Vs = NULL,
  Vt = NULL,
  Ds = NULL,
  Dt = NULL,
  ltPrior = NULL,
  lsPrior = NULL,
  sigmaPrior = NULL,
  noisePrior = NULL,
  mh_sd_r = NULL
)
```

## Arguments

- X:

  Fixed-effect design matrix with one row per observation.

- y:

  Non-negative integer count response.

- coords:

  Spatial coordinate matrix with one row per full spatial level,
  including the baseline level omitted from `Vs`.

- nsim:

  Total number of MCMC iterations; must exceed `burn`.

- burn:

  Number of burn-in iterations.

- use_count_gp:

  Whether to include GP random effects in the count component.

- use_inflation_gp:

  Whether to include GP random effects in the zero-inflation component.

- thin:

  Store every `thin`-th iteration after burn-in.

- kern:

  Kernel function accepting a distance matrix and length scale.

- save_ypred:

  Whether to save posterior predictive draws.

- print_iter:

  Report progress every `print_iter` iterations when `print_progress` is
  `TRUE`.

- print_progress:

  Whether to report MCMC progress via
  [`message()`](https://rdrr.io/r/base/message.html); these reports can
  be silenced with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html).

- Vs:

  Spatial random-effect design matrix with one row per observation and
  the baseline spatial column omitted.

- Vt:

  Temporal random-effect design matrix with one row per observation and
  the baseline temporal column omitted.

- Ds:

  Spatial distance matrix for all spatial levels, including the
  baseline; its diagonal must be zero.

- Dt:

  Temporal distance matrix for all temporal levels, including the
  baseline; its diagonal must be zero.

- ltPrior:

  List with `max`, `mh_sd`, `a`, and `b` for temporal length-scale prior
  and proposal controls.

- lsPrior:

  List with `max`, `mh_sd`, `a`, and `b` for spatial length-scale prior
  and proposal controls.

- sigmaPrior:

  List with `a` and `b` inverse-gamma prior parameters for GP scales.

- noisePrior:

  List with `a`, `b`, and `mh_sd` for the GP noise-ratio prior and
  proposal.

- mh_sd_r:

  Proposal standard deviation for the negative-binomial dispersion
  parameter.

## Value

A list containing posterior MCMC draws:

- Alpha:

  Fixed-effect coefficients for the zero-inflation component.

- Beta:

  Fixed-effect coefficients for the count component.

- A, B:

  Spatial and temporal random effects for the zero-inflation component.

- C, D:

  Spatial and temporal random effects for the count component.

- L1t, L2t:

  Temporal GP length scales for the zero-inflation and count components.

- Sigma1t, Sigma2t:

  Temporal GP scale parameters.

- L1s, L2s:

  Spatial GP length scales for the zero-inflation and count components.

- Sigma1s, Sigma2s:

  Spatial GP scale parameters.

- R:

  Negative-binomial dispersion parameter.

- at_risk:

  Latent at-risk indicator draws, included when `save_ypred` is `TRUE`.

- Y_pred:

  Posterior predictive count draws, included when `save_ypred` is
  `TRUE`.

## Details

At least one spatial or temporal GP must be active in the count or
zero-inflation component. Models with no active GP are outside this
entry point and signal an error that points to standard GLM software
instead.

## Examples

``` r
# A small synthetic spatial example.
cells <- expand.grid(spatial = seq_len(4), replicate = seq_len(12))
Vs <- diag(4)[cells$spatial, -1, drop = FALSE]
y <- c(
  0, 0, 0, 0, 15, 1, 0, 0, 6, 0, 0, 0, 1, 6, 0, 0,
  0, 0, 0, 1, 11, 0, 0, 0, 2, 0, 0, 0, 8, 0, 5, 0,
  0, 4, 0, 0, 34, 0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 3
)
X <- cbind("(Intercept)" = 1, x = as.numeric(scale(seq_along(y))))
coords <- rbind(c(0, 0), c(1000, 0), c(0, 1000), c(1000, 1000))

set.seed(1)
fit <- ZINB_GP(
  X = X,
  y = y,
  coords = coords,
  nsim = 5,
  burn = 1,
  thin = 2,
  use_count_gp = TRUE,
  use_inflation_gp = FALSE,
  Vs = Vs,
  Ds = as.matrix(stats::dist(coords))
)
names(fit)
#> [1] "Alpha"   "Beta"    "C"       "L2s"     "Sigma2s" "Noise2s" "R"      
```
