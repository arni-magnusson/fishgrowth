#' Growth Cessation Model
#'
#' Fit a growth cessation model (GCM) to tags and otoliths.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#' @param t age (vector).
#' @param L0 predicted length at age 0.
#' @param rmax shape parameter that determines the initial slope.
#' @param k shape parameter that determines how quickly the growth curve reaches
#'        the asymptotic maximum.
#' @param t50 shape parameter that determines the logistic function midpoint.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{gcm} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{gcm_curve} and \code{gcm_objfun}
#' are called by the main function to calculate the regression curve and
#' objective function value. The user can also call the auxiliary functions
#' directly for plotting and model exploration.
#'
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
#'   \item \code{log_age}, age at release of tagged individuals (vector)
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
#' The \code{gcm} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{gcm_curve} function returns a numeric vector of predicted length at
#' age.
#'
#' The \code{gcm_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The growth cessation model (Maunder et al. 2018) predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ r_{\max}\!\left[\,\frac{\log\!\left(1+e^{-k\:\!t_{50}}
#'       \right) \;-\;\log\!\left(1+e^{k(t-t_{50})}\right)}{k}\;+\;t\:\right]}{
#'       L0 + rmax * ((log(1 + exp(-k*t50)) - log(1 + exp(k*(t-t50)))) / k + t)}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{
#'       sigma_L = alpha + beta * Lhat}
#'
#' where the slope is \eqn{\beta=(\sigma_2-\sigma_1) /
#' (L_\mathrm{long}-L_\mathrm{short})}{beta = (sigma_2-sigma_1) /
#' (L_long-L_short)} and the intercept is \eqn{\alpha=\sigma_1
#' - \beta L_\mathrm{short}}{alpha = sigma_1 - beta * L_short}.
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
#' Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W., Aires-da-Silva,
#' A.M., Minte-Vera, C.V., and Campana, S.E. (2018).
#' The growth cessation model: a growth model for species showing a near
#' cessation in growth with application to bigeye tuna (\emph{Thunnus obesus}).
#' \emph{Marine Biology}, \bold{165}, 76.
#' \doi{10.1007/s00227-018-3336-9}.
#'
#' @seealso
#' \code{gcm}, \code{\link{richards}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}} are alternative growth models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @examples
#' # Explore initial parameter values
#' plot(len~age, otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' x <- seq(0, 4, 0.1)
#' points(lenRel~I(lenRel/60), tags_ex, col=4)
#' points(lenRec~I(lenRel/60+liberty), tags_ex, col=3)
#' lines(x, gcm_curve(x, L0=20, rmax=120, k=2, t50=0), lty=2)
#'
#' # Prepare parameters and data
#' init <- list(L0=20, log_rmax=log(120), log_k=log(2), t50=0,
#'              log_sigma_1=log(1), log_sigma_2=log(1),
#'              log_age=log(tags_ex$lenRel/60))
#' dat <- list(Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'             liberty=tags_ex$liberty, Aoto=otoliths_ex$age,
#'             Loto=otoliths_ex$len, L_short=30, L_long=60)
#' gcm_objfun(init, dat)
#'
#' # Fit model
#' model <- gcm(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4,iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model, getReportCovariance=FALSE)
#'
#' # Plot results
#' plot(len~age, otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, tags_ex$lenRel, col=4)
#' points(report$age+tags_ex$liberty, tags_ex$lenRec, col=3)
#' lines(report$age_seq, report$curve, lwd=2)
#'
#' # Model summary
#' est <- report[c("L0", "rmax", "k", "t50", "sigma_1", "sigma_2")]
#' est
#' fit[-1]
#' head(summary(sdreport), 6)
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @export

gcm <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  MakeADFun(wrap(gcm_objfun, data=data), par, silent=silent, ...)
}

#' @rdname gcm
#'
#' @export

gcm_curve <- function(t, L0, rmax, k, t50)
{
  L0 + rmax * ((log(1 + exp(-k*t50)) - log(1 + exp(k*(t-t50)))) / k + t)
}

#' @rdname gcm
#'
#' @export

gcm_objfun <- function(par, data)
{
  # Extract parameters
  L0 <- par$L0
  rmax <- exp(par$log_rmax)
  k <- exp(par$log_k)
  t50 <- par$t50
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- exp(par$log_sigma_2)
  age <- tryCatch(exp(par$log_age), errror=as.null)

  # Extract data
  Lrel <- data$Lrel
  Lrec <- data$Lrec
  liberty <- data$liberty
  Aoto <- data$Aoto
  Loto <- data$Loto
  L_short <- data$L_short
  L_long <- data$L_long

  # Calculate sigma coefficients (sigma = a + b*L)
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
  nll_Lrel <- tryCatch(-dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE), error=as.null)
  nll_Lrec <- tryCatch(-dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE), error=as.null)
  nll_Loto <- tryCatch(-dnorm(Loto, Loto_hat, sigma_Loto, TRUE), error=as.null)
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
