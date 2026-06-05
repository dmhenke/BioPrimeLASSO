#' Identify phi value
#'
#' Cross-validation across a range of phi values; selects the best phi by
#' minimum RMSE relative to a linear baseline across the phi range.
#'
#' @param rmse A numeric matrix of z-scored RMSE values (phi x fold) returned
#'   by \code{\link{get_rmse_from_xvalid}}.
#' @param phi_range A numeric vector of phi values tested. Suggested range [0, 1].
#'
#' @return A single numeric value: the phi corresponding to the greatest
#'   negative deviation from the linear trend in median RMSE.
#' @export
#'
#' @examples
#' \dontrun{
#'   phi_range <- seq(0, 1, length = 30)
#'   rmse <- get_rmse_from_xvalid(X, y, penalties,
#'     phi_range = phi_range,
#'     lambda_min = 0.5,
#'     n_folds = 10)
#'   best_phi <- find_best_phi_rmse(rmse, phi_range)
#' }
find_best_phi_rmse <- function(rmse, phi_range){
  median_rmse <- apply(rmse, 1, median)

  aframe <- data.frame(
    phi = phi_range,
    rmse = median_rmse)

  afit <- stats::lm(rmse ~ phi, data = aframe[c(1, nrow(aframe)), ])
  preds <- stats::predict(afit, aframe)
  diff_rmse <- aframe$rmse - preds

  best_rmse_phi <- phi_range[which(diff_rmse == min(diff_rmse))]

  return(best_rmse_phi)
}
