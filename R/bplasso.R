#' Bio-primed LASSO
#'
#' Fits a biologically informed LASSO model that incorporates protein-protein
#' interaction (PPI) network scores as penalty weights during cross-validated
#' regularization. The method prioritizes features that are both statistically
#' relevant and biologically plausible, as measured by their proximity to a
#' gene of interest in a curated PPI network (e.g., STRING DB).
#'
#' @param X A numeric matrix of predictor variables (n x p), as in `glmnet::glmnet()`, typically
#'   scaled omic features (e.g., copy number variation). Rows are observations
#'   (cell lines), columns are features (genes).
#' @param y A numeric response vector, as in `glmnet::glmnet()`, of length n containing the response variable
#'   (e.g., CRISPR dependency scores from DEMETER2 or Chronos).
#' @param scores A numeric vector of biological priority scores, gene-specific association from `get_scores()`, for each
#'   predictor in \code{X}, typically derived from \code{\link{get_scores}}.
#'   Higher scores imply greater biological relevance.
#' @param n_folds An integer specifying the number of cross-validation folds.
#'   Default is \code{10}.
#' @param phi_range A numeric vector defining the grid of phi values to search
#'   over. Phi controls the degree of bio-priming: \code{phi = 0} reduces the
#'   model to standard LASSO; \code{phi = 1} applies maximal biological
#'   weighting. Default is \code{seq(0, 1, length = 30)}.
#'
#' @return A named list of 3 elements:
#' \describe{
#'   \item{phi}{Best phi value as chosen by \code{\link{find_best_phi_rmse}}.}
#'   \item{lambda}{Lambda value identified as \code{lambda.min} by \code{glmnet::cv.glmnet()}.}
#'   \item{betas}{A data.frame containing baseline LASSO coefficients (\code{betas}) and
#'     bio-primed LASSO coefficients (\code{betas_pen}) for all columns of \code{X}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   bplasso(scale(X),
#'     y,
#'     scores,
#'     n_folds = 10,
#'     phi_range = seq(0, 1, length = 30))
#' }
#'
bplasso <- function(X, y, scores,
                     n_folds = 10,
                     phi_range = seq(0, 1, length = 30)){
  # Choose lambda
  lambda_min <- find_lambda(X, y, plot = FALSE)
  message(paste("Lambda min:", round(lambda_min, 4)))
  # Fit baseline LASSO
  afit <- glmnet::glmnet(
    X, y,
    alpha = 1,
    lambda = lambda_min)
  betas <- afit$beta[, 1]
  if (sum(betas) == 0) {
    message("All betas are zero.")
    return(NA)
  }
  penalties <- scores[match(colnames(X), names(scores))]
  names(penalties) <- colnames(X)
  penalties[is.na(penalties)] <- 0
  penalties <- penalties / max(scores)
  # Choose best phi
  rmse <- get_rmse_from_xvalid(
    X, y, penalties, phi_range = phi_range, lambda_min = lambda_min, n_folds = n_folds)
  if (length(unique(dim(rmse) == dim(stats::na.omit(rmse)))) == 2) {
    message("Missing values in RMSE matrix.")
    return(NA)
  }
  best_phi <- find_best_phi_rmse(rmse, phi_range)
  message(paste("Best phi based on RMSE:", round(best_phi, 4)))
  # Run LASSO with updated lambda & phi
  afit <- glmnet::glmnet(
    X,
    y,
    alpha = 1,
    lambda = lambda_min,
    penalty.factor = 1 - penalties * best_phi)
  betas_pen <- afit$beta[, 1]
  return(list(phi = best_phi,
              lambda = lambda_min,
              betas = data.frame(betas, betas_pen)))
}
