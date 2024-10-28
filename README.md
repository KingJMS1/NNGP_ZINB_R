# ZINB NNGP Bayesian Model

This package implements the model described in "A Framework of Zero-Inflated Bayesian Negative
Binomial Regression Models For Spatiotemporal Data" by Qing He and Hsin-Hsiung Huang (2023). https://doi.org/10.1016/j.jspi.2023.106098. 

This package is a work in progress, feel free to create an issue if you have suggestions or notice any problems.

## Installation Instructions


## Example Use
Detailed examples with full code can be found in the examples folder. For basic usage, the model needs the following:

```
X       Other Predictor variables
y       Zero inflated count response
coords  Spatial coordinates for NNGP
Vs      Spatially varying predictor variables 
        (e.g. one-hot indication of which location this is for varying intercept), 
        wrapped in sparseMatrix from Matrix R package. 
        Will be multiplied by the spatial random effects for prediction.
Vt      Temporal varying predictor variables, wrapped in sparseMatrix from Matrix R package. 
        Will be multiplied by the temporal random effects for prediction.
Ds      Spatial distance matrix, diagonal should be 0, 
        off diagonal is distance between elements i and j in space, inputs to the spatial NNGP kernel
Dt      Temporal distance matirx, diagonal should be 0, 
        off diagonal is distance between elements i and j in time, inputs to the temporal GP kernel
M       How many neighbors to allow in the spatial NNGP algorithm, defaults to 10.
nsim    How long to run MCMC in total, must be greater than burn.
burn    How long to run MCMC before saving samples.
```

Given all of the above, predictions can then be found via:
```r
library(ZINB.GP)
output <- ZINB_NNGP(X, y, coords, Vs, Vt, Ds, Dt, M = M, nsim, burn, save_ypred = TRUE)
predictions <- output$Y_pred
```