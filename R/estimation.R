#' Draw a Zero-Inflated Negative-Binomial Response
#'
#' Generates one posterior predictive count for each pair of at-risk and count
#' linear predictors using the parameterization fitted by `ZINB_GP()`.
#'
#' @param eta_at_risk Numeric vector of linear predictors for the Bernoulli
#'   at-risk component.
#' @param eta_count Numeric vector of linear predictors for the
#'   negative-binomial count component.
#' @param r Positive negative-binomial dispersion parameter.
#'
#' @return A non-negative integer vector with one predictive count per linear
#'   predictor pair.
#' @keywords internal
draw_zinb_response <- function(eta_at_risk, eta_count, r) {
    eta_at_risk <- c(eta_at_risk)
    eta_count <- c(eta_count)
    N <- length(eta_at_risk)

    if (length(eta_count) != N) {
        stop("`eta_at_risk` and `eta_count` must have the same length.")
    }

    at_risk <- stats::rbinom(N, 1, sigmoid(eta_at_risk))
    y <- integer(N)

    if (any(at_risk == 1)) {
        y[at_risk == 1] <- stats::rnbinom(
            sum(at_risk == 1),
            size = r,
            prob = sigmoid(-eta_count[at_risk == 1])
        )
    }

    y
}


#' @title Numerically stable sigmoid function
#' @description Computes the inverse-logit transformation while clipping its
#' output away from zero and one.
#' @param eta Numeric vector or matrix of linear predictors.
#' @return `eta` transformed to probabilities in `[1e-6, 1 - 1e-6]`.
sigmoid <- function(eta) {
    eta <- pmin(700, eta)
    return(pmax(1e-6, pmin(1 - 1e-6, exp(eta) / (1 + exp(eta))))) 
}

#' Draw a GP at New Levels Conditional on a Stored Draw
#'
#' @param effect Numeric GP draw at the observed nonbaseline levels.
#' @param distance Augmented distance matrix whose observed levels precede its
#'   new levels.
#' @param length_scale GP length scale.
#' @param sigma GP marginal standard deviation.
#' @param noise_ratio Fraction of variance assigned to the structured kernel.
#' @param kern Kernel function accepting squared distances and a length scale.
#' @return A numeric draw at the new levels.
#' @keywords internal
draw_conditional_gp <- function(
    effect,
    distance,
    length_scale,
    sigma,
    noise_ratio,
    kern
) {
    # The augmented matrix places fitted levels before prediction levels.
    effect <- c(effect)
    distance <- as.matrix(distance)
    n_observed <- length(effect)
    n_new <- nrow(distance) - n_observed

    observed_index <- seq_len(n_observed)
    new_index <- n_observed + seq_len(n_new)
    squared_distance <- distance^2

    # Construct the observed, cross, and new covariance blocks. Independent
    # noise contributes only within the observed and new diagonal blocks.
    K_observed <- sigma^2 * noise_mix(
        kern(
            squared_distance[observed_index, observed_index, drop = FALSE],
            length_scale
        ),
        noise_ratio
    )
    K_cross <- sigma^2 * noise_ratio * kern(
        squared_distance[new_index, observed_index, drop = FALSE],
        length_scale
    )
    K_new <- sigma^2 * noise_mix(
        kern(
            squared_distance[new_index, new_index, drop = FALSE],
            length_scale
        ),
        noise_ratio
    )

    # Apply the Gaussian conditioning rule without explicitly inverting K.
    solved_effect <- solve(K_observed, effect)
    solved_cross <- solve(K_observed, t(K_cross))

    conditional_mean <- c(K_cross %*% solved_effect)
    conditional_covariance <- K_new - K_cross %*% solved_cross

    # Remove small numerical asymmetry before taking the eigendecomposition.
    conditional_covariance <- (
        conditional_covariance + t(conditional_covariance)
    ) / 2

    # Clamp round-off-level negative eigenvalues to obtain a stable MVN draw.
    decomposition <- eigen(conditional_covariance, symmetric = TRUE)
    scales <- sqrt(pmax(decomposition$values, 0))

    conditional_mean + c(
        decomposition$vectors %*% (scales * stats::rnorm(n_new))
    )
}

#' Mark Model Output as a ZINB-GP Fit
#'
#' @param output Model output list.
#' @return The output with its model class attached.
#' @keywords internal
new_zinb_gp_fit <- function(output) {
    class(output) <- c("zinb_gp_fit", "list")
    output
}

#' Predict at New Locations and Times
#'
#' Generates posterior predictive draws at new spatial locations, times, or
#' both. For every retained MCMC iteration, the method draws each active random
#' effect from its GP predictive distribution conditional on the stored fitted
#' random effect, adds it to the corresponding fixed-effect predictor, and
#' draws a zero-inflated negative-binomial response.
#'
#' @param object A fitted zinb_gp_fit object returned by ZINB_GP().
#' @param X Fixed-effect design matrix for the prediction observations.
#' @param Ds_new Augmented spatial distance matrix produced by
#'   make_prediction_inputs(), or NULL for a temporal-only fit.
#' @param Dt_new Augmented temporal distance matrix produced by
#'   make_prediction_inputs(), or NULL for a spatial-only fit.
#' @param Vs_new Spatial design matrix produced by make_prediction_inputs().
#' @param Vt_new Temporal design matrix produced by make_prediction_inputs().
#' @param kern Kernel function used to fit the model. The default is kernel();
#'   supply the fitting kernel again if a custom kernel was used.
#' @param ... Unused.
#' @return An object of class zinb_gp_prediction containing posterior draws
#'   Y_pred, eta_at_risk, and eta_count. Active predicted GP effects are
#'   returned under A, B, C, and D, matching the fitted object.
#' @examples
#' # Fit a small spatial model, then predict at two new locations.
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
#'
#' coords_new <- rbind(c(500, 500), c(250, 750))
#' inputs <- make_prediction_inputs(coords = coords, coords_new = coords_new)
#' X_new <- cbind("(Intercept)" = 1, x = c(-0.5, 0.5))
#' predictions <- predict(
#'   fit,
#'   X = X_new,
#'   Ds_new = inputs$Ds_new,
#'   Vs_new = inputs$Vs_new
#' )
#' predictions$Y_pred
#' @export
#' @importFrom stats predict
predict.zinb_gp_fit <- function(
    object,
    X,
    Ds_new = NULL,
    Dt_new = NULL,
    Vs_new = NULL,
    Vt_new = NULL,
    kern = NULL,
    ...
) {
    # Normalize the fixed-effect design and identify the posterior draw count.
    X <- as.matrix(X)
    n_prediction <- nrow(X)
    n_draws <- nrow(object$Alpha)
    kern <- nullcheck(kern, kernel)

    # Validate the shared fixed-effect and dispersion draws first.
    if (n_prediction < 1L || ncol(X) != ncol(object$Alpha) ||
            any(!is.finite(X))) {
        stop("X must be a finite matrix with one column per fitted fixed effect.")
    }

    if (!identical(dim(object$Alpha), dim(object$Beta)) ||
            length(object$R) != n_draws) {
        stop("The fitted fixed-effect and dispersion draws are not aligned.")
    }

    validate_prediction_dimension <- function(
        distance,
        design,
        effects,
        dimension_name
    ) {
        # A dimension is active when either model component stores its effects.
        active_effects <- effects[vapply(effects, function(name) {
            !is.null(object[[name]])
        }, logical(1))]

        if (length(active_effects) == 0L) {
            if (!is.null(distance) || !is.null(design)) {
                stop("The fitted model has no ", dimension_name, " GP.")
            }
            return(NULL)
        }

        if (is.null(distance) || is.null(design)) {
            stop(
                "Both the new ", dimension_name,
                " distance and design matrices are required."
            )
        }
        design <- as.matrix(design)
        distance <- as.matrix(distance)

        # The design maps observations to the unique new GP levels.
        if (nrow(design) != n_prediction || ncol(design) < 1L ||
                any(!is.finite(design))) {
            stop(
                "The new ", dimension_name,
                " design must have one row per row of X."
            )
        }

        # Components sharing a dimension must use the same fitted GP levels.
        observed_levels <- unique(vapply(active_effects, function(name) {
            ncol(object[[name]])
        }, integer(1)))

        if (length(observed_levels) != 1L) {
            stop("Stored ", dimension_name, " random-effect draws are not aligned.")
        }
        expected_dimension <- observed_levels + ncol(design)

        if (nrow(distance) != expected_dimension ||
                ncol(distance) != expected_dimension ||
                any(!is.finite(distance)) || any(distance < 0) ||
                !isTRUE(all.equal(distance, t(distance)))) {
            stop(
                "The augmented ", dimension_name,
                " distance matrix must be symmetric with one row for every ",
                "observed nonbaseline and new GP level."
            )
        }

        list(distance = distance, design = design)
    }

    # Validate the coordinate inputs for every active model dimension.
    spatial <- validate_prediction_dimension(
        Ds_new, Vs_new, c("A", "C"), "spatial"
    )
    temporal <- validate_prediction_dimension(
        Dt_new, Vt_new, c("B", "D"), "temporal"
    )

    # Match each stored effect to its coordinate dimension and hyperparameters.
    effect_specification <- list(
        A = list(
            dimension = spatial,
            length = "L1s",
            sigma = "Sigma1s",
            noise = "Noise1s"
        ),
        B = list(
            dimension = temporal,
            length = "L1t",
            sigma = "Sigma1t",
            noise = "Noise1t"
        ),
        C = list(
            dimension = spatial,
            length = "L2s",
            sigma = "Sigma2s",
            noise = "Noise2s"
        ),
        D = list(
            dimension = temporal,
            length = "L2t",
            sigma = "Sigma2t",
            noise = "Noise2t"
        )
    )

    # Draw each active random effect at the unique new GP levels.
    predicted_effects <- list()

    for (effect_name in names(effect_specification)) {
        if (is.null(object[[effect_name]])) {
            next
        }

        specification <- effect_specification[[effect_name]]

        predicted_effects[[effect_name]] <- matrix(
            NA_real_,
            nrow = n_draws,
            ncol = ncol(specification$dimension$design)
        )

        for (draw in seq_len(n_draws)) {
            predicted_effects[[effect_name]][draw, ] <- draw_conditional_gp(
                effect = object[[effect_name]][draw, ],
                distance = specification$dimension$distance,
                length_scale = object[[specification$length]][draw],
                sigma = object[[specification$sigma]][draw],
                noise_ratio = object[[specification$noise]][draw],
                kern = kern
            )
        }
    }

    # Allocate one prediction row for every retained posterior iteration.
    eta_at_risk <- eta_count <- matrix(
        NA_real_, nrow = n_draws, ncol = n_prediction
    )
    Y_pred <- matrix(NA_integer_, nrow = n_draws, ncol = n_prediction)

    for (draw in seq_len(n_draws)) {
        # Start with fixed effects, then add each active predicted GP effect.
        eta_at_risk[draw, ] <- c(X %*% object$Alpha[draw, ])
        eta_count[draw, ] <- c(X %*% object$Beta[draw, ])

        if (!is.null(predicted_effects$A)) {
            eta_at_risk[draw, ] <- eta_at_risk[draw, ] +
                c(Vs_new %*% predicted_effects$A[draw, ])
        }

        if (!is.null(predicted_effects$B)) {
            eta_at_risk[draw, ] <- eta_at_risk[draw, ] +
                c(Vt_new %*% predicted_effects$B[draw, ])
        }

        if (!is.null(predicted_effects$C)) {
            eta_count[draw, ] <- eta_count[draw, ] +
                c(Vs_new %*% predicted_effects$C[draw, ])
        }

        if (!is.null(predicted_effects$D)) {
            eta_count[draw, ] <- eta_count[draw, ] +
                c(Vt_new %*% predicted_effects$D[draw, ])
        }

        # Convert the two predictors into a posterior response draw.
        Y_pred[draw, ] <- draw_zinb_response(
            eta_at_risk = eta_at_risk[draw, ],
            eta_count = eta_count[draw, ],
            r = object$R[draw]
        )
    }

    # Return both response predictions and the latent draws that produced them.
    prediction <- c(
        list(
            Y_pred = Y_pred,
            eta_at_risk = eta_at_risk,
            eta_count = eta_count
        ),
        predicted_effects
    )

    class(prediction) <- c("zinb_gp_prediction", "list")
    prediction
}
