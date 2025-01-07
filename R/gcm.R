#' Growth Cessation Model
#'
#' Fit a growth cessation model (GCM) to otoliths and/or tags.
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
#'   \item \code{log_sigma_2} (*), growth variability at length \code{L_long}
#'   \item \code{log_age} (*), age at release of tagged individuals (vector)
#' }
#'
#' *: The parameter \code{log_sigma_2} can be omitted to estimate growth
#' variability that does not vary with length. The parameter vector
#' \code{log_age} can be omitted to fit to otoliths only.
#'
#' The \code{data} list contains the following elements:
#' \itemize{
#'   \item \code{Aoto} (*), age from otoliths (vector)
#'   \item \code{Loto} (*), length from otoliths (vector)
#'   \item \code{Lrel} (*), length at release of tagged individuals (vector)
#'   \item \code{Lrec} (*), length at recapture of tagged individuals (vector)
#'   \item \code{liberty} (*), time at liberty of tagged individuals in years
#'         (vector)
#'   \item \code{L_short}, length where sd(length) is \code{sigma_1}
#'   \item \code{L_long}, length where sd(length) is \code{sigma_2}
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
##'
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
#' \deqn{\hat L_t ~=~ r_{\max}\!\left[\,\frac{\log\left(1+e^{-kt_{50}}
#'       \right) \;-\;\log\left(1+e^{k(t-t_{50})}\right)}{k}\;+\;t\:\right]}{
#'       L0 + rmax * ((log(1 + exp(-k*t50)) - log(1 + exp(k*(t-t50)))) / k + t)}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{
#'       sigma_L = alpha + beta * Lhat}
#'
#' where the slope is \eqn{\beta=(\sigma_2-\sigma_1) /
#' (L_\mathrm{long}-L_\mathrm{short})}{beta = (sigma_2-sigma_1) /
#' (L_long-L_short)} and the intercept is \eqn{\alpha=\sigma_1 - \beta
#' L_\mathrm{short}}{alpha = sigma_1 - beta * L_short}. Alternatively, growth
#' variability can be modelled as a constant
#' (\eqn{\sigma_L=\sigma_1}{sigma_L=sigma_1}) that does not vary with length,
#' see \code{log_sigma_2} above.
#'
#' The negative log-likelihood is calculated by comparing the observed and
#' predicted lengths:
#' \preformatted{
#'   nll_Loto <- -dnorm(Loto, Loto_hat, sigma_Loto, TRUE)
#'   nll_Lrel <- -dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE)
#'   nll_Lrec <- -dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE)
#'   nll <- sum(nll_Loto) + sum(nll_Lrel) + sum(nll_Lrec)
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
#' \code{gcm}, \code{\link{gompertz}},
#' \code{\link{richards}}/\code{\link{richardso}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}}/\code{\link{vonberto}} are alternative growth models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{fishgrowth-package}} gives an overview of the package.
#'
#' @examples
#' # Demonstrate GCM curve shapes
#' plot(0:10, gcm_curve(0:10, L0=0, rmax=20, k=2, t50=2))
#' plot(0:10, gcm_curve(0:10, L0=0, rmax=100, k=2, t50=2))
#'
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
#' dat <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len,
#'             Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'             liberty=tags_ex$liberty, L_short=30, L_long=60)
#' gcm_objfun(init, dat)
#'
#' # Fit model
#' model <- gcm(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
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
#' #############################################################################
#'
#' # Stepwise estimation procedure, described by Maunder et al. (2018)
#' # - estimate L0 and rmax using linear regression on younger fish
#' # - estimate k and t0 using GCM and all data, keeping L0 and rmax fixed
#'
#' # Estimate L0 and rmax
#' plot(otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' fm <- lm(len~age, otoliths_ex)
#' abline(fm)
#' L0 <- coef(fm)[[1]]
#' rmax <- coef(fm)[[2]]
#'
#' # Explore initial parameter values (k, t50, age)
#' t <- seq(0, 4, by=0.1)
#' points(t, gcm_curve(t, L0, rmax, k=3, t50=2), col="gray")
#' points(lenRel~I(lenRel/50), tags_ex, col=4)
#' points(lenRec~I(lenRel/50+liberty), tags_ex, col=3)
#'
#' # Prepare parameters
#' init_step <- list(L0=L0, log_rmax=log(rmax), log_k=log(3), t50=2,
#'                   log_sigma_1=log(1), log_sigma_2=log(1),
#'                   log_age=log(tags_ex$lenRel/50))
#'
#' # Fit model
#' map <- list(L0=factor(NA), log_rmax=factor(NA))  # fix L0 and rmax
#' model_step <- gcm(init_step, dat, map=map)
#' fit_step <- nlminb(model_step$par, model_step$fn, model_step$gr,
#'                    control=list(eval.max=1e4,iter.max=1e4))
#' report_step <- model_step$report()
#' sdreport_step <- sdreport(model_step, getReportCovariance=FALSE)
#'
#' # Plot results
#' plot(len~age, otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' points(report_step$age, tags_ex$lenRel, col=4)
#' points(report_step$age+tags_ex$liberty, tags_ex$lenRec, col=3)
#' lines(report_step$age_seq, report_step$curve, lwd=2)
#'
#' # Model summary
#' est_step <- report_step[c("L0", "rmax", "k", "t50", "sigma_1", "sigma_2")]
#' est_step
#' fit_step[-1]
#' head(summary(sdreport_step), 6)
#'
#' #############################################################################
#'
#' # Fit to otoliths only
#' init_oto <- list(L0=20, log_rmax=log(120), log_k=log(2), t50=0,
#'                  log_sigma_1=log(1), log_sigma_2=log(1))
#' dat_oto <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len,
#'                 L_short=30, L_long=60)
#' model_oto <- gcm(init_oto, dat_oto)
#' fit_oto <- nlminb(model_oto$par, model_oto$fn, model_oto$gr,
#'                   control=list(eval.max=1e4, iter.max=1e4))
#' model_oto$report()[c("L0", "rmax", "k", "t50")]
#'
#' #############################################################################
#'
#' # Fit to tags only
#' init_tags <- list(L0=log(20), log_rmax=log(120), log_k=log(2), t50=0,
#'                   log_sigma_1=log(1), log_sigma_2=log(1),
#'                   log_age=log(tags_ex$lenRel/60))
#' dat_tags <- list(Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'                  liberty=tags_ex$liberty, L_short=30, L_long=60)
#' model_tags <- gcm(init_tags, dat_tags)
#' fit_tags <- nlminb(model_tags$par, model_tags$fn, model_tags$gr,
#'                    control=list(eval.max=1e4, iter.max=1e4))
#' model_tags$report()[c("L0", "rmax", "k", "t50")]
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
  sigma_2 <- if(is.null(par$log_sigma_2)) NULL else exp(par$log_sigma_2)

  # Extract data
  L_short <- data$L_short
  L_long <- data$L_long

  # Calculate sigma coefficients (sigma = a + b*L)
  if(is.null(sigma_2))
  {
    sigma_slope <- 0  # if user did not pass log_sigma_2 then use constant sigma
    sigma_intercept <- sigma_1
  }
  else
  {
    sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)
    sigma_intercept <- sigma_1 - L_short * sigma_slope
  }

  # Initialize likelihood
  nll <- 0

  # Calculate curve
  age_seq = seq(0, 10, 1/365)  # age 0-10 years, day by day
  curve <- gcm_curve(age_seq, L0, rmax, k, t50)

  # Report quantities of interest
  REPORT(L0)
  REPORT(rmax)
  REPORT(k)
  REPORT(t50)
  REPORT(L_short)
  REPORT(L_long)
  REPORT(sigma_1)
  REPORT(sigma_2)
  REPORT(age_seq)
  REPORT(curve)
  ADREPORT(curve)

  # Model includes otolith data
  if(!is.null(data$Aoto) && !is.null(data$Loto))
  {
    # data
    Aoto <- data$Aoto
    Loto <- data$Loto
    # Lhat
    Loto_hat <- gcm_curve(Aoto, L0, rmax, k, t50)
    # sigma
    sigma_Loto <- sigma_intercept + sigma_slope * Loto_hat
    # nll
    nll_Loto <- -dnorm(Loto, Loto_hat, sigma_Loto, TRUE)
    nll <- nll + sum(nll_Loto)
    # report
    REPORT(Aoto)
    REPORT(Loto)
    REPORT(Loto_hat)
    REPORT(sigma_Loto)
    REPORT(nll_Loto)
  }

  # Model includes tagging data
  if(!is.null(par$log_age) && !is.null(data$Lrel) &&
     !is.null(data$Lrec) && !is.null(data$liberty))
  {
    # par
    age <- exp(par$log_age)
    # data
    Lrel <- data$Lrel
    Lrec <- data$Lrec
    liberty <- data$liberty
    # Lhat
    Lrel_hat <- gcm_curve(age, L0, rmax, k, t50)
    Lrec_hat <- gcm_curve(age+liberty, L0, rmax, k, t50)
    # sigma
    sigma_Lrel <- sigma_intercept + sigma_slope * Lrel_hat
    sigma_Lrec <- sigma_intercept + sigma_slope * Lrec_hat
    # nll
    nll_Lrel <- -dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE)
    nll_Lrec <- -dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE)
    nll <- nll + sum(nll_Lrel) + sum(nll_Lrec)
    # report
    REPORT(age)
    REPORT(Lrel)
    REPORT(Lrec)
    REPORT(liberty)
    REPORT(Lrel_hat)
    REPORT(Lrec_hat)
    REPORT(sigma_Lrel)
    REPORT(sigma_Lrec)
    REPORT(nll_Lrel)
    REPORT(nll_Lrec)
  }

  nll
}
