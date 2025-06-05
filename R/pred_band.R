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
#' @note
#' The variability of length at age (\code{sigma}) increases linearly with
#' length:
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{
#'       sigma_L = alpha + beta * Lhat}
#'
#' This calculation of \code{sigma} is demonstrated in the example below.
#'
#' The \code{lower} and \code{upper} limits of the prediction band are
#' calculated as \eqn{\hat L \pm 1.96\sigma_L}{Lhat +/- 1.96 * sigma_L} at the
#' 95\% significance level.
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}, \code{\link{gompertzo}},
#' \code{\link{richards}}, \code{\link{richardso}}, \code{\link{schnute3}},
#' \code{\link{vonbert}}, and \code{\link{vonberto}} are alternative growth
#' models.
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
#' # Calculate 95% prediction band
#' x <- seq(1, 18, 0.5)
#' band <- pred_band(x, model)
#'
#' # Plot 95% prediction band
#' areaplot::confplot(cbind(lower,upper)~age, band, xlim=c(0,18), ylim=c(0,100),
#'          ylab="len", col="mistyrose")
#' points(len~age, otoliths_had, pch=16, col="#0080a010")
#' lines(Lhat~age, band, lwd=2, col=2)
#' lines(lower~age, band, lty=1, lwd=0.5, col=2)
#' lines(upper~age, band, lty=1, lwd=0.5, col=2)
#'
#' # Calculate sigma by hand
#' report <- model$report()
#' alpha <- report$sigma_intercept
#' beta <- report$sigma_slope
#' Lhat <- band$Lhat
#' alpha + beta * Lhat  # same values as band$sigma calculated by pred_band()
#'
#' @importFrom stats qnorm
#'
#' @export

pred_band <- function(age, model, level=0.95)
{
  report <- model$report()
  report$t <- age
  Lhat <- with(report, eval(body(report$curve)))
  sigma <- with(report, sigma_intercept + sigma_slope * Lhat)
  sigma[sigma < 0] <- NA
  p <- 1 - (1-level)/2
  lower <- Lhat - qnorm(p) * sigma
  upper <- Lhat + qnorm(p) * sigma
  data.frame(age, Lhat, sigma, lower, upper)
}
