#' Identify phi value
#' Cross-validation across range of phi values, find best phi by RMSE
#'
#' @param rmse Correlation matrix along phi value range from `get_rmse_from_xvalid`
#' @param phi_range Array of phi ranges to test.  Suggested limits [0,1]
#'
#' @return A numeric value
#' @export
#'
#' @examples
#' phirange <- seq(0,1,length = 30)
#' rmse <- get_rmse_from_xvalid(X, y, penalties,
#'   phi_range = phirange,
#'   lambda_min = 0.5,
#'   n_folds = 10)
#'   best_phi <- find_best_phi_rmse(rmse, phirange, plot = F)

find_best_phi_rmse <- function(rmse, phi_range){
  median_rmse <- apply(rmse, 1, median)

  aframe <- data.frame(
    phi = phi_range,
    rmse = median_rmse)

  afit <- lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- predict(afit, aframe)
  diff_rmse <- aframe$rmse - preds

  best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]

  return(best_rmse_phi)
}
