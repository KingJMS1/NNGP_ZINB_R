#' ZINB_GP
#' @description Fits the zero-inflated negative-binomial Gaussian-process model
#' described in \doi{10.1016/j.jspi.2023.106098}. The supplied
#' spatial and temporal design and distance matrices determine which Gaussian
#' processes are included.
#'
#' @details At least one spatial or temporal GP must be active in the count or
#'   zero-inflation component. Models with no active GP are outside this entry
#'   point and signal an error that points to standard GLM software instead.
#'
#' @param X Fixed-effect design matrix with one row per observation.
#' @param y Non-negative integer count response.
#' @param coords Spatial coordinate matrix with one row per full spatial level,
#'   including the baseline level omitted from `Vs`.
#' @param nsim Total number of MCMC iterations; must exceed `burn`.
#' @param burn Number of burn-in iterations.
#' @param use_count_gp Whether to include GP random effects in the count component.
#' @param use_inflation_gp Whether to include GP random effects in the zero-inflation component.
#' @param thin Store every `thin`-th iteration after burn-in.
#' @param kern Kernel function accepting a distance matrix and length scale.
#' @param save_ypred Whether to save posterior predictive draws.
#' @param print_iter Report progress every `print_iter` iterations when
#'   `print_progress` is `TRUE`.
#' @param print_progress Whether to report MCMC progress via `message()`; these
#'   reports can be silenced with `suppressMessages()`.
#' @param Vs Spatial random-effect design matrix with one row per observation
#'   and the baseline spatial column omitted.
#' @param Vt Temporal random-effect design matrix with one row per observation
#'   and the baseline temporal column omitted.
#' @param Ds Spatial distance matrix for all spatial levels, including the
#'   baseline; its diagonal must be zero.
#' @param Dt Temporal distance matrix for all temporal levels, including the
#'   baseline; its diagonal must be zero.
#' @param ltPrior List with `max`, `mh_sd`, `a`, and `b` for temporal length-scale prior and proposal controls.
#' @param lsPrior List with `max`, `mh_sd`, `a`, and `b` for spatial length-scale prior and proposal controls.
#' @param sigmaPrior List with `a` and `b` inverse-gamma prior parameters for GP scales.
#' @param noisePrior List with `a`, `b`, and `mh_sd` for the GP noise-ratio prior and proposal.
#' @param mh_sd_r Proposal standard deviation for the negative-binomial dispersion parameter.
#' @return A list containing posterior MCMC draws:
#' \describe{
#'   \item{Alpha}{Fixed-effect coefficients for the zero-inflation component.}
#'   \item{Beta}{Fixed-effect coefficients for the count component.}
#'   \item{A, B}{Spatial and temporal random effects for the zero-inflation component.}
#'   \item{C, D}{Spatial and temporal random effects for the count component.}
#'   \item{L1t, L2t}{Temporal GP length scales for the zero-inflation and count components.}
#'   \item{Sigma1t, Sigma2t}{Temporal GP scale parameters.}
#'   \item{L1s, L2s}{Spatial GP length scales for the zero-inflation and count components.}
#'   \item{Sigma1s, Sigma2s}{Spatial GP scale parameters.}
#'   \item{R}{Negative-binomial dispersion parameter.}
#'   \item{at_risk}{Latent at-risk indicator draws, included when `save_ypred` is `TRUE`.}
#'   \item{Y_pred}{Posterior predictive count draws, included when `save_ypred` is `TRUE`.}
#' }
#'
#' @examples
#' # A small synthetic spatial example.
#' cells <- expand.grid(spatial = seq_len(4), replicate = seq_len(12))
#' Vs <- diag(4)[cells$spatial, -1, drop = FALSE]
#' y <- c(
#'   0, 0, 0, 0, 15, 1, 0, 0, 6, 0, 0, 0, 1, 6, 0, 0,
#'   0, 0, 0, 1, 11, 0, 0, 0, 2, 0, 0, 0, 8, 0, 5, 0,
#'   0, 4, 0, 0, 34, 0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 3
#' )
#' X <- cbind("(Intercept)" = 1, x = as.numeric(scale(seq_along(y))))
#' coords <- rbind(c(0, 0), c(1000, 0), c(0, 1000), c(1000, 1000))
#'
#' set.seed(1)
#' fit <- ZINB_GP(
#'   X = X,
#'   y = y,
#'   coords = coords,
#'   nsim = 5,
#'   burn = 1,
#'   thin = 2,
#'   use_count_gp = TRUE,
#'   use_inflation_gp = FALSE,
#'   Vs = Vs,
#'   Ds = as.matrix(stats::dist(coords))
#' )
#' names(fit)
#' @export
ZINB_GP <- function(X, y, coords, nsim = 5000, burn = 1000, use_count_gp = TRUE, use_inflation_gp = FALSE, thin = 1, kern = NULL, save_ypred = FALSE, print_iter = 100, print_progress = FALSE, Vs = NULL, Vt = NULL, Ds = NULL, Dt = NULL, ltPrior = NULL, lsPrior = NULL, sigmaPrior = NULL, noisePrior = NULL, mh_sd_r = NULL) 
{
    errMsg <- "You must specify at least 1 GP to use. Use optimization GLM software like INLA, MASS, glmmTMB, pscl, etc. to fit the model instead."
    no_gp_design <- is.null(Vs) && is.null(Vt)
    no_gp_component <- !use_count_gp && !use_inflation_gp
    if (no_gp_design || no_gp_component) {
        stop(errMsg)
    }

    validate_zinb_inputs(X, y, nsim, burn, thin)
    if (!is.null(Vs)) {
        validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")
        validate_spatial_coordinates(coords, Vs, Ds)
    }
    if (!is.null(Vt)) {
        validate_gp_design(Vt, Dt, nrow(X), "Vt", "Dt")
    }

    if (is.null(Vs))
    {
        if (is.null(Vt))
        {
            # No GPS
            stop(errMsg)
        }
        else
        {
            # Only temporal GPs
            if (use_count_gp)
            {
                if (use_inflation_gp)
                {
                    # Both types of GP, temporal only
                    results <- ZINB_GP_spatial(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, B = results$A, D = results$C, L1t = results$L1s, Sigma1t = results$Sigma1s, Noise1t = results$Noise1s, L2t = results$L2s, Sigma2t = results$Sigma2s, Noise2t = results$Noise2s, R = results$R)
                    if (save_ypred) {
                        toReturn$Y_pred <- results$Y_pred
                        toReturn$at_risk <- results$at_risk
                    }
                    return(new_zinb_gp_fit(toReturn))
                }
                else
                {
                    # Only count GP, temporal only
                    results <- ZINB_GP_spatial_count(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, D = results$C, L2t = results$L2s, Sigma2t = results$Sigma2s, Noise2t = results$Noise2s, R = results$R)
                    if (save_ypred) {
                        toReturn$Y_pred <- results$Y_pred
                        toReturn$at_risk <- results$at_risk
                    }
                    return(new_zinb_gp_fit(toReturn))
                }
            }
            else
            {
                # No count GP
                if (use_inflation_gp)
                {
                    # Inflation GP only
                    results <- ZINB_GP_spatial_inflation(X, y, Vt, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, sigmaPrior, noisePrior, mh_sd_r, kern)
                    toReturn <- list(Alpha = results$Alpha, Beta = results$Beta, B = results$A, L1t = results$L1s, Sigma1t = results$Sigma1s, Noise1t = results$Noise1s, R = results$R)
                    if (save_ypred) {
                        toReturn$Y_pred <- results$Y_pred
                        toReturn$at_risk <- results$at_risk
                    }
                    return(new_zinb_gp_fit(toReturn))
                }
                else
                {
                    # No GPS
                    stop(errMsg)
                }
            }
        }
    }
    else
    {
        if (is.null(Vt))
        {
            # Only spatial GPs
            if (use_count_gp)
            {
                # Do use count gps
                if (use_inflation_gp)
                {
                    # Both GP types, spatial only
                    return(ZINB_GP_spatial(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # Only count GP, spatial only
                    return(ZINB_GP_spatial_count(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
            }
            else
            {
                # No count GPS
                if (use_inflation_gp)
                {
                    # Inflation gp only, spatial only
                    return(ZINB_GP_spatial_inflation(X, y, Vs, Ds, nsim, burn, thin, save_ypred, print_iter, print_progress, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # No GPS
                    stop(errMsg)
                }
            }
        }
        else
        {
            # All GPs
            if (use_count_gp)
            {
                # Do use count gps
                if (use_inflation_gp)
                {
                    # Full GPs
                    return(ZINB_GP_orig(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # Count GP only
                    return(ZINB_GP_count(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
            }
            else
            {
                # No Count GPS
                if (use_inflation_gp)
                {
                    return(ZINB_GP_inflation(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin, save_ypred, print_iter, print_progress, ltPrior, lsPrior, sigmaPrior, noisePrior, mh_sd_r, kern))
                }
                else
                {
                    # No GPS
                    stop(errMsg)
                }
            }
        }
    }
}
