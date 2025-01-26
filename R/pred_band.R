#' Prediction Band
#'
#' Calculate a prediction band for a fitted growth curve.
#'
#' @param model a fitted growth model.
#' @param age a vector of ages to calculate the prediction band.
#' @param model a fitted growth model.
#' @param level significance level.
#'
#' @return
#' A data frame containing five columns:
#' \item{age}{age}
#' \item{Lhat}{predicted length}
#' \item{sigma}{growth variability}
#' \item{lower}{lower limit of prediction band}
#' \item{upper}{upper limit of prediction band}
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}/\code{\link{gompertzo}},
#' \code{\link{richards}}/\code{\link{richardso}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}}/\code{\link{vonberto}} are alternative growth models.
#'
#' \code{\link{fishgrowth-package}} gives an overview of the package.
#'
#' @examples
#' # Fit a model
#' init <- list(log_L1=log(20), log_L2=log(70), log_k=log(0.1),
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_had$age, Loto=otoliths_had$len, t1=1, t2=10)
#' model <- vonbert(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#'
#' # Calculate prediction band
#' x <- seq(1, 18, 0.5)
#' pred_band(x, model)
#'
#' @importFrom stats qnorm
#'
#' @export

pred_band <- function(age, model, level)
{
  report <- model$report()
  report$t <- age
  Lhat <- with(report, eval(body(report$curve)))
  sigma <- with(report, sigma_intercept + sigma_slope * Lhat)
  sigma[sigma < 0] <- NA
  lower <- Lhat - qnorm(0.025) * sigma
  upper <- Lhat + qnorm(0.025) * sigma
  data.frame(age, Lhat, sigma, lower, upper)
}
