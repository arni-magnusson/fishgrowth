#' Growth Cessation Model
#'
#' Fit a growth cessation model (GCM) to otoliths and/or tags.
#'
#' @param par a parameter list.
#' @param data a data list.
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
#'   \item \code{log_sigma_min}, growth variability at the shortest observed
#'         length in the data
#'   \item \code{log_sigma_max} (*), growth variability at the longest observed
#'         length in the data
#'   \item \code{log_age} (*), age at release of tagged individuals (vector)
#' }
#'
#' *: The parameter \code{log_sigma_max} can be omitted to estimate growth
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
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
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
#' \deqn{\hat L_t ~=~ L_0 ~+~ r_{\max}\!\left[\,\frac{\log\left(1+e^{-kt_{50}}
#'       \right) \;-\;\log\left(1+e^{k(t-t_{50})}\right)}{k}\;+\;t\:\right]}{
#'       Lhat = L0 + rmax * ((log(1 + exp(-k*t50)) - log(1 + exp(k*(t-t50))))
#'       / k + t)}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{sigma = alpha + beta * Lhat}
#'
#' where the intercept is \eqn{\alpha=\sigma_{\min} - \beta L_{\min}}{alpha =
#' sigma_min - beta * L_min}, the slope is
#' \eqn{\beta=(\sigma_{\max}-\sigma_{\min}) / (L_{\max}-L_{\min})}{beta =
#' (sigma_max-sigma_min) / (L_max-L_min)}, and \eqn{L_{\min}}{L_min} and
#' \eqn{L_{\max}}{L_max} are the shortest and longest observed lengths in the
#' data. Alternatively, growth variability can be modelled as a constant
#' \eqn{\sigma_L=\sigma_{\min}}{sigma_L=sigma_min} that does not vary with
#' length, by omitting \code{log_sigma_max} from the parameter list (see above).
#'
#' The negative log-likelihood objective function integrates (sums) the
#' likelihood components from the otoliths and tags,
#'
#' \deqn{f ~=~ \sum_{i=1}^{N_\mathrm{oto}}\, 0.5\log(2\pi) \,+\,
#'       \log\sigma_i \,+\, \frac{(L_i-\hat L_i)^2}{2\sigma_i^2}}{f
#'       = sum(0.5*log(2*pi) + log(sigma_i) + ((L_i-Lhat_i)^2)/(2*sigma_i^2))}
#' \deqn{~~~~ ~+\,~ \sum_{j=1}^{N_\mathrm{tag}}\, 0.5\log(2\pi) \,+\,
#'       \log\sigma_j \,+\, \frac{(L_j-\hat L_j)^2}{2\sigma_j^2}}{
#'       + sum(0.5*log(2*pi) + log(sigma_j) + ((L_j-Lhat_j)^2)/(2*sigma_j^2))}
#' \deqn{~~~~ ~+\,~ \sum_{k=1}^{N_\mathrm{tag}}\, 0.5\log(2\pi) \,+\,
#'       \log\sigma_k \,+\, \frac{(L_k-\hat L_k)^2}{2\sigma_k^2}}{
#'       + sum(0.5*log(2*pi) + log(sigma_k) + ((L_k-Lhat_k)^2)/(2*sigma_k^2))}
#'
#' where \eqn{L_i} are length measurements from the otolith data, \eqn{L_j} are
#' length measurements at tag release, and \eqn{L_k} are length measurements at
#' tag recapture. \if{html,latex}{\eqn{N_\mathrm{oto}} is the number of fish in
#' the otolith data and \eqn{N_\mathrm{tag}} is the number of fish in the
#' tagging data.}
#'
#' @references
#' Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W., Aires-da-Silva,
#' A.M., Minte-Vera, C.V., and Campana, S.E. (2018).
#' The growth cessation model: a growth model for species showing a near
#' cessation in growth with application to bigeye tuna (\emph{Thunnus obesus}).
#' \emph{Marine Biology}, \bold{165}, 76.
#' \doi{10.1007/s00227-018-3336-9}.
#'
#' The \code{\link{fishgrowth-package}} help page includes references describing
#' the parameter estimation method.
#'
#' @seealso
#' \code{gcm}, \code{\link{gompertz}}, \code{\link{gompertzo}},
#' \code{\link{richards}}, \code{\link{richardso}}, \code{\link{schnute3}},
#' \code{\link{vonbert}}, and \code{\link{vonberto}} are alternative growth
#' models.
#'
#' \code{\link{pred_band}} calculates a prediction band for a fitted growth
#' model.
#'
#' \code{\link{otoliths_had}}, \code{\link{otoliths_skj}}, and
#' \code{\link{tags_skj}} are example datasets.
#'
#' \code{\link{fishgrowth-package}} gives an overview of the package.
#'
#' @examples
#' # Model 1: Fit to haddock otoliths
#'
#' # Explore initial parameter values
#' plot(len~age, otoliths_had, xlim=c(0,18), ylim=c(0,105), pch=16,
#'      col="#0080a010")
#' x <- seq(1, 18, 0.1)
#' lines(x, gcm_curve(x, L0=5, rmax=20, k=0.15, t50=0), lty=3)
#'
#' # Prepare parameters and data
#' init <- list(L0=5, log_rmax=log(20), log_k=log(0.15), t50=-1,
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_had$age, Loto=otoliths_had$len)
#' gcm_objfun(init, dat)
#'
#' # Fit model
#' model <- gcm(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' Lhat <- with(report, gcm_curve(x, L0, rmax, k, t50))
#' lines(x, Lhat, lwd=2, col=2)
#' legend("bottomright", c("initial curve","model fit"), col=c(1,2), lty=c(3,1),
#'        lwd=c(1,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L0", "rmax", "k", "t50", "sigma_min", "sigma_max")]
#' fit[-1]
#' summary(sdreport)
#'
#' # Plot 95% prediction band
#' band <- pred_band(x, model)
#' areaplot::confplot(cbind(lower,upper)~age, band, xlim=c(0,18), ylim=c(0,100),
#'          ylab="len", col="mistyrose")
#' points(len~age, otoliths_had, pch=16, col="#0080a010")
#' lines(Lhat~age, band, lwd=2, col=2)
#' lines(lower~age, band, lty=1, lwd=0.5, col=2)
#' lines(upper~age, band, lty=1, lwd=0.5, col=2)
#'
#' #############################################################################
#'
#' # Model 2: Fit to skipjack otoliths and tags
#'
#' # Explore initial parameter values
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' x <- seq(0, 4, 0.1)
#' points(lenRel~I(lenRel/60), tags_skj, col=4)
#' points(lenRec~I(lenRel/60+liberty), tags_skj, col=3)
#' lines(x, gcm_curve(x, L0=20, rmax=120, k=2, t50=0), lty=2)
#' legend("bottomright", c("otoliths","tag releases","tag recaptures",
#'        "initial curve"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,2),
#'        lwd=c(1.2,1.2,1.2,1), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Prepare parameters and data
#' init <- list(L0=20, log_rmax=log(120), log_k=log(4), t50=0,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len,
#'             Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty)
#' gcm_objfun(init, dat)
#'
#' # Fit model
#' model <- gcm(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, gcm_curve(x, L0, rmax, k, t50))
#' lines(x, Lhat, lwd=2)
#' legend("bottomright", c("otoliths","tag releases","tag recaptures",
#'        "model fit"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,1),
#'        lwd=c(1.2,1.2,1.2,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L0", "rmax", "k", "t50", "sigma_min", "sigma_max")]
#' fit[-1]
#' head(summary(sdreport), 6)
#'
#' #############################################################################
#'
#' # Model 3: Stepwise estimation procedure, described by Maunder et al. (2018)
#' # - estimate L0 and rmax using linear regression on younger fish
#' # - estimate k and t50 using GCM and all data, keeping L0 and rmax fixed
#'
#' # Estimate L0 and rmax
#' plot(otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' fm <- lm(len~age, otoliths_skj)
#' abline(fm)
#' L0 <- coef(fm)[[1]]
#' rmax <- coef(fm)[[2]]
#'
#' # Explore initial parameter values (k, t50, age)
#' t <- seq(0, 4, by=0.1)
#' points(t, gcm_curve(t, L0, rmax, k=3, t50=2), col="gray")
#' points(lenRel~I(lenRel/50), tags_skj, col=4)
#' points(lenRec~I(lenRel/50+liberty), tags_skj, col=3)
#' legend("bottomright", c("otoliths","tag releases","tag recaptures",
#'        "linear regression (otoliths)"), col=c(1,4,3,1), pch=c(1,1,1,NA),
#'        lty=c(0,0,0,1), lwd=c(1.2,1.2,1.2,2), bty="n", inset=0.02,
#'        y.intersp=1.25)
#'
#' # Prepare parameters
#' init <- list(L0=L0, log_rmax=log(rmax), log_k=log(3), t50=2,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/50))
#'
#' # Fit model
#' map <- list(L0=factor(NA), log_rmax=factor(NA))  # fix L0 and rmax
#' model <- gcm(init, dat, map=map)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4,iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, gcm_curve(x, L0, rmax, k, t50))
#' lines(x, Lhat, lwd=2)
#' legend("bottomright", c("otoliths","tag releases","tag recaptures",
#'        "model fit"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,1),
#'        lwd=c(1.2,1.2,1.2,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L0", "rmax", "k", "t50", "sigma_min", "sigma_max")]
#' fit[-1]
#' head(summary(sdreport), 6)
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#'
#' @export

gcm <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  if(is.null(par$log_sigma_min))
    stop("'par' list must include 'log_sigma_min'")
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
  sigma_min <- exp(par$log_sigma_min)
  sigma_max <- if(is.null(par$log_sigma_max)) NULL else exp(par$log_sigma_max)

  # Set L_min and L_max to minimum and maximum lengths in data
  L_min <- min(c(data$Loto, data$Lrel, data$Lrec))
  L_max <- max(c(data$Loto, data$Lrel, data$Lrec))

  # Calculate sigma coefficients (sigma = a + b*L)
  if(is.null(sigma_max))
  {
    sigma_slope <- 0  # if user did not pass log_sigma_max then constant sigma
    sigma_intercept <- sigma_min
  }
  else
  {
    sigma_slope <- (sigma_max - sigma_min) / (L_max - L_min)
    sigma_intercept <- sigma_min - L_min * sigma_slope
  }

  # Initialize likelihood
  nll <- 0

  # Report quantities of interest
  type <- "gcm"
  curve <- gcm_curve
  REPORT(type)
  REPORT(curve)
  REPORT(L0)
  REPORT(rmax)
  REPORT(k)
  REPORT(t50)
  REPORT(L_min)
  REPORT(L_max)
  REPORT(sigma_min)
  REPORT(sigma_max)
  REPORT(sigma_intercept)
  REPORT(sigma_slope)

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
