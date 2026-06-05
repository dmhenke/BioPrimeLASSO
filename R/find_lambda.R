#' Calculate best lambda (lambda.min)
#'
#' Runs \code{glmnet::cv.glmnet()} with alpha = 1 (LASSO) and returns the
#' lambda value that minimizes cross-validated error.
#'
#' @param X A numeric matrix of predictor variables, as in \code{glmnet::glmnet()}.
#' @param y A numeric response vector, as in \code{glmnet::glmnet()}.
#' @param plot Logical. If \code{TRUE}, plots the cross-validation curve.
#'   Default is \code{FALSE}.
#'
#' @return A single numeric value: \code{lambda.min} from \code{glmnet::cv.glmnet()}.
#' @export
#'
#' @examples
#' \dontrun{
#'   lambda_min <- find_lambda(scale(X), y, plot = FALSE)
#' }
find_lambda <- function(X, y, plot = FALSE){
  fitcv <- glmnet::cv.glmnet(
    X, y,
    alpha = 1,
    lambda = NULL)
  if (plot) plot(fitcv, xvar = "lambda", label = TRUE)

  fitcv$lambda.min
}
