#' Extract DLNM Predictions from Model Objects
#'
#' @description
#' Extracts cross-basis coefficients and variance-covariance matrices corresponding to a specific
#' cross-basis object from a fitted regression model (e.g., mixed-effects models fitted via \code{lme4} or \code{mgcv}).
#' Computes exposure-response and lag-response predictions centered at a specified value using the \code{dlnm::crosspred()} function.
#'
#' @param model A fitted regression model object (e.g., \code{lmerMod}, \code{gam}, or \code{glm}) containing
#' cross-basis terms among its fixed coefficients.
#' @param cb_object A cross-basis object created by \code{dlnm::crossbasis()} used in fitting \code{model}.
#' @param cb_name Character string specifying the prefix/name of the cross-basis terms in the model summary.
#' Default is \code{"final_cb"}.
#' @param cen_val Numeric. Centering point (reference value) for exposure-response functions. Default is \code{0}.
#'
#' @return A \code{crosspred} object containing predicted exposure-lag-response associations, overall cumulative
#' associations, relative risks/ratios, and corresponding standard errors/confidence intervals.
#'
#' @importFrom dlnm crosspred
#'
#' @examples
#' \dontrun{
#' library(dlnm)
#'
#' # Mock example using standard linear model and dlnm cross-basis
#' set.seed(123)
#' x <- rnorm(100, mean = 20, sd = 5)
#' y <- rnorm(100, mean = 100, sd = 10)
#'
#' # Create cross-basis matrix
#' cb <- crossbasis(x, lag = 3, argvar = list(fun = "lin"), arglag = list(fun = "integer"))
#' fit <- lm(y ~ cb)
#'
#' # Extract predictions centered at x = 20
#' pred <- get_dlnm_pred(
#'   model = fit,
#'   cb_object = cb,
#'   cb_name = "cb",
#'   cen_val = 20
#' )
#'
#' # Inspect overall effect
#' head(pred$allRRfit)
#' }
#'
#' @export
get_dlnm_pred <- function(model, cb_object, cb_name = "final_cb", cen_val = 0) {
  indices <- grep(paste0("^", cb_name), names(coef(model)))
  if(length(indices) == 0) {
    stop("Errore: Non ho trovato coefficienti che corrispondono a: ", cb_name)}
  beta <- coef(model)[indices]
  V <- vcov(model)[indices, indices]
  names(beta) <- gsub(paste0("^", cb_name), "", names(beta))
  rownames(V) <- colnames(V) <- names(beta)
  pred <- crosspred(cb_object, coef = beta, vcov = V,
                    model.link = "log", cen = cen_val)
  return(pred)}
