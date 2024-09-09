#' RSME for the range of lambda
#' RMSE, via N-fold cross-validation, along a range of phi, given lambda.  Returned as a Z-score
#'
#' @param X matrix as in `glmnet::glmnet()`
#' @param y response array as in `glmnet::glmnet()`
#' @param penalties User-supplied penalty scores.  Range [0-1]
#' @param lambda_min User-supplied lambda value from `find_lambda()`
#' @param phi_range Array of phi ranges to test.  Suggested limits [0,1]
#' @param n_folds Number of cross-validations
#'
#' @return a matrix of n-fold by phi_range z-scored RMSE, from cross-validation glmnet
#' @export
#'
#' @examples
#' #' phirange <- seq(0,1,length = 30)
#'   get_rmse_from_xvalid(X, y, penalties,
#'     phi_range = phirange,
#'     lambda_min = 0.5,
#'     n_folds = 10)
#'
get_rmse_from_xvalid <- function(
    X, y, penalties, lambda_min, phi_range, n_folds = 10){
  asplits <- suppressWarnings(split(sample(1:nrow(X)), 1:n_folds))
  rmse <- do.call(cbind, lapply(names(asplits), function(x){
    train <- unlist(asplits[setdiff(names(asplits), x)])
    test <- unlist(asplits[x])
    do.call(rbind, lapply(phi_range, function(phi){
      lasso_tr <- glmnet(
        X[train,],
        y[train],
        lambda = lambda_min,
        penalty.factor = 1 - penalties * phi)
      pred <- predict(lasso_tr, X[test,])
      rmse <- sqrt(apply((y[test] - pred)^2, 2, mean))
      return(rmse)
    }))
  }))
  rmse <- apply(rmse, 2, function(x)
    (x - mean(x))/sd(x))
  return(rmse)
}
