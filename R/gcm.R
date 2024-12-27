#' Growth Cessation Model
#'
#' Fit a growth cessation model (GCM) to tags and otoliths.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#'
#' @details
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{L0}, predicted length at age 0
#'   \item \code{log_rmax}, shape parameter that determines the initial slope
#'   \item \code{log_k}, shape parameter that determines how quickly the growth
#'         curve reaches the asymptotic maximum
#'   \item \code{t50}, shape parameter that determines the logistic function
#'         midpoint
#'   \item \code{log_sigma_1}, growth variability at length \code{L_short}
#'   \item \code{log_sigma_2}, growth variability at length \code{L_long}
#'   \item \code{log_age} age at release of tagged individuals (vector)
#' }
#'
#' The \code{data} list contains the following elements:
#' \itemize{
#'   \item \code{Lrel}, length at release of tagged individuals (vector)
#'   \item \code{Lrec}, length at recapture of tagged individuals (vector)
#'   \item \code{liberty}, time at liberty of tagged individuals in years
#'         (vector)
#'   \item \code{Aoto}, age from otoliths (vector)
#'   \item \code{Loto}, length from otoliths (vector)
#'   \item \code{L_short}, length where sd(length) is \code{sigma_1}
#'   \item \code{L_long}, length where sd(length) is \code{sigma_2}
#' }
#'
#' @return
#' TMB model object, produced by \code{\link[RTMB]{MakeADFun}}.
#'
#' @note
#'
#' The growth cessation model (Maunder et al. 2018) predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ r_{\max}\!\left[\,\frac{\log\!\left(1+e^{-k\:\!t_{50}}
#'       \right) \;-\;\log\!\left(1+e^{k(t-t_{50})}\right)}{k}\;+\;t\:\right]}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L_t}
#'
#' where the slope is \eqn{\beta=(\sigma_2-\sigma_1) /
#' (L_\mathrm{long}-L_\mathrm{short})} and the intercept is \eqn{\alpha=\sigma_1
#' - \beta L_\mathrm{short}}.
#'
#' The negative log-likelihood is calculated by comparing the observed and
#' predicted lengths:
#' \preformatted{
#'   nll_Lrel <- -dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE)
#'   nll_Lrec <- -dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE)
#'   nll_Loto <- -dnorm(Loto, Loto_hat, sigma_Loto, TRUE)
#'   nll <- sum(nll_Lrel) + sum(nll_Lrec) + sum(nll_Loto)
#' }
#'
#' @references
#' Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W., Aires‐da‐Silva,
#' A.M., Minte‐Vera, C.V., and Campana, S.E. (2018).
#' The growth cessation model: a growth model for species showing a near
#' cessation in growth with application to bigeye tuna (\emph{Thunnus obesus}).
#' \emph{Marine Biology}, \bold{165}, 76.
#' \doi{10.1007/s00227-018-3336-9}.
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @seealso
#' \code{\link{richards}} and \code{\link{vonbert}} are alternative growth
#' models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @export

gcm <- function(par, data)
{
  wrap <- function(objfun, ...) function(par) objfun(par, ...)
  MakeADFun(wrap(gcm_objfun, data=data), par, silent=TRUE)
}

gcm_curve <- function(t, L0, rmax, k, t50)
{
  L0 + rmax * ((log(1 + exp(-k*t50)) - log(1 + exp(k*(t-t50)))) / k + t)
}

gcm_objfun <- function(par, data)
{
  # Extract parameters
  L0 <- par$L0
  rmax <- exp(par$log_rmax)
  k <- exp(par$log_k)
  t50 <- par$t50
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- exp(par$log_sigma_2)
  age <- exp(par$log_age)

  # Extract data
  Lrel <- data$Lrel
  Lrec <- data$Lrec
  liberty <- data$liberty
  Aoto <- data$Aoto
  Loto <- data$Loto
  L_short <- data$L_short
  L_long <- data$L_long

  # Calculate sigma coefficients (sigma = a + b*age)
  sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)
  sigma_intercept <- sigma_1 - L_short * sigma_slope

  # Calculate Lhat and sigma
  Lrel_hat <- gcm_curve(age, L0, rmax, k, t50)
  Lrec_hat <- gcm_curve(age+liberty, L0, rmax, k, t50)
  Loto_hat <- gcm_curve(Aoto, L0, rmax, k, t50)
  sigma_Lrel <- sigma_intercept + sigma_slope * Lrel_hat
  sigma_Lrec <- sigma_intercept + sigma_slope * Lrec_hat
  sigma_Loto <- sigma_intercept + sigma_slope * Loto_hat

  # Calculate likelihoods
  nll_Lrel <- -dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE)
  nll_Lrec <- -dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE)
  nll_Loto <- -dnorm(Loto, Loto_hat, sigma_Loto, TRUE)
  nll <- sum(nll_Lrel) + sum(nll_Lrec) + sum(nll_Loto)

  # Calculate curve
  age_seq = seq(0, 10, 1/365)  # age 0-10 years, day by day
  curve <- gcm_curve(age_seq, L0, rmax, k, t50)

  # Report quantities of interest
  REPORT(L0)
  REPORT(rmax)
  REPORT(k)
  REPORT(t50)
  REPORT(age)
  REPORT(liberty)
  REPORT(Lrel)
  REPORT(Lrec)
  REPORT(Aoto)
  REPORT(Loto)
  REPORT(Lrel_hat)
  REPORT(Lrec_hat)
  REPORT(Loto_hat)
  REPORT(L_short)
  REPORT(L_long)
  REPORT(sigma_1)
  REPORT(sigma_2)
  REPORT(sigma_Lrel)
  REPORT(sigma_Lrec)
  REPORT(sigma_Loto)
  REPORT(nll_Lrel)
  REPORT(nll_Lrec)
  REPORT(nll_Loto)
  REPORT(age_seq)
  REPORT(curve)
  ADREPORT(curve)

  nll
}
