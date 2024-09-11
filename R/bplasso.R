#' Bio-primed LASSO
#'
#' @param X matrix as in `glmnet::glmnet()`
#' @param y response array as in `glmnet::glmnet()`
#' @param scores Gene specific association scores, from `get_scores()`
#' @param n_folds Number of cross-validations
#' @param phi_range  Array of phi ranges to test.  Suggested limits [0,1]
#'
#' @return list of 5 elements
#' 1: "phi" = best phi value as chosen by `find_best_phi_rmse()`
#' 2: "lambda" = lambda value as identified as `lambda.min` by `glmnet::cv.glmnet()`
#' 3: "betas" = data.frame containing baseline lasso betas "betas" and bio-primed lasso betas "betas_pen" for all columns of matrix `x`
#'
#' @export
#'
#' @examples
#'  bplasso(scale(X),
#'   y,
#'   scores,
#'   n_folds = 10,
#'   phi_range = seq(0, 1, length = 30))
#'
bplasso <- function(X, y, scores,
                     n_folds = 10,
                     phi_range = seq(0, 1, length = 30)){
  # Choose lambda
  lambda_min <- find_lambda(X, y, plot = F)
  print(paste("Lambda min:", round(lambda_min,4)))
  # Fit baseline LASSO
  afit <- glmnet(
    X, y,
    alpha = 1,
    lambda = lambda_min)
  betas <- afit$beta[,1]
  if(sum(betas) == 0){
    print("All betas are zero.")
    return(NA)
  }
  penalties <- scores[match(colnames(X), names(scores))]
  names(penalties) <- colnames(X)
  penalties[is.na(penalties)] <- 0
  penalties <- penalties/max(scores)
  # Choose best phi
  rmse <- get_rmse_from_xvalid(
    X, y, penalties, phi_range = phi_range, lambda_min = lambda_min, n_folds = n_folds)
  if(length(unique(dim(rmse) == dim(na.omit(rmse)))) == 2){
    print("Missing values in correlation.")
    return(NA)
  }
  best_phi <- find_best_phi_rmse(rmse, phi_range)
  print(paste("Best phi based on RMSE:", round(best_phi, 4)))
  # Run LASSO with updated lambda & phi
  afit <- glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi)
  betas_pen <- afit$beta[,1]
  return(list(phi = best_phi,
              lambda = lambda_min,
              betas = data.frame(betas, betas_pen)))
}
