#' Calculate best lambda: lambda min
#'
#' @param X matrix as in `glmnet`
#' @param y response array as in `glmnet::glmnet()`
#' @param plot plot the resulting cross-validation for glmnet
#'
#' @return A numeric value
#' @export
#'
find_lambda <- function(X, y, plot = F){
  fitcv <- cv.glmnet(
    X, y,
    alpha = 1,
    lambda = NULL)
  if(plot) plot(fitcv, xvar = "lambda", label = T)

  fitcv$lambda.min
}
