#' RMSE across a range of phi via N-fold cross-validation
#'
#' For each fold and each phi value, fits a penalized LASSO model on the
#' training split and computes RMSE on the held-out test split. Results are
#' returned as z-scored RMSE values (per fold) to allow comparison across folds.
#'
#' @param X A numeric matrix of predictor variables, as in \code{glmnet::glmnet()}.
#' @param y A numeric response vector, as in \code{glmnet::glmnet()}.
#' @param penalties A numeric vector of penalty scores in the range [0, 1],
#'   one per column of \code{X}. Typically derived from \code{\link{get_scores}}
#'   and normalized to [0, 1].
#' @param lambda_min A single numeric lambda value, typically from \code{\link{find_lambda}}.
#' @param phi_range A numeric vector of phi values to test. Suggested range [0, 1].
#' @param n_folds An integer number of cross-validation folds. Default is \code{10}.
#'
#' @return A numeric matrix of z-scored RMSE values with dimensions
#'   \code{length(phi_range)} x \code{n_folds}.
#' @export
#'
#' @examples
#' \dontrun{
#'   phi_range <- seq(0, 1, length = 30)
#'   get_rmse_from_xvalid(X, y, penalties,
#'     phi_range = phi_range,
#'     lambda_min = 0.5,
#'     n_folds = 10)
#' }
#'
get_rmse_from_xvalid <- function(
    X, y, penalties, lambda_min, phi_range, n_folds = 10){
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  rmse <- do.call(cbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    do.call(rbind, lapply(phi_range, function(phi){
      lasso_tr <- glmnet::glmnet(
        X[train, ],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
      pred <- stats::predict(lasso_tr, X[test, ])
      rmse <- sqrt(apply((y[test] - pred)^2, 2, mean))
      return(rmse)
    }))
  }))
  rmse <- apply(rmse, 2, function(x)
    (x - mean(x)) / stats::sd(x))
  return(rmse)
}
