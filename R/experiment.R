#' Experimental Growth Model
#'
#' Fit an experimental growth model from another package.
#'
#' @param par a parameter list.
#' @param data a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @seealso
#' \code{\link[linear]{linear}} describes the experimental growth model.
#'
#' \code{\link{gcm}}, \code{\link{gompertz}}, \code{\link{gompertzo}},
#' \code{\link{richards}}, \code{\link{richardso}}, \code{experiment},
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
#' lines(x, experiment_curve(x, L1=15, L2=70, t1=1, t2=10), lty=3)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(15), log_L2=log(70),
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_had$age, Loto=otoliths_had$len, t1=1, t2=10)
#' experiment_objfun(init, dat)
#'
#' # Fit model
#' model <- experiment(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' Lhat <- with(report, experiment_curve(x, L1, L2, t1, t2))
#' lines(x, Lhat, lwd=2, col=2)
#' legend("bottomright", c("initial curve","model fit"), lty=c(3,1), lwd=c(1,2),
#'        col=c(1,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L1", "L2", "sigma_min", "sigma_max")]
#' fit[-1]
#' summary(sdreport)
#'
#' # Plot 95% prediction band
#' band <- pred_band(x, model)
#' areaplot::confplot(cbind(lower,upper)~age, band, xlim=c(0,18), ylim=c(0,100),
#'          ylab="len", col="mistyrose")
#' points(len~age, otoliths_had, xlim=c(0,18), ylim=c(0,100),
#'        pch=16, col="#0080a010")
#' lines(x, Lhat, lwd=2, col=2)
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
#' lines(x, experiment_curve(x, L1=30, L2=110, t1=0, t2=4), lty=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "initial curve"), lty=c(0,0,0,2), pch=c(1,1,1,NA), lwd=c(1.2,1.2,1.2,1),
#'        col=c(1,4,3,1), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(30), log_L2=log(110),
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len,
#'             Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4)
#' experiment_objfun(init, dat)
#'
#' # Fit model
#' model <- experiment(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, experiment_curve(x, L1, L2, t1, t2))
#' lines(x, Lhat, lwd=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "model fit"), lty=c(0,0,0,1), pch=c(1,1,1,NA), lwd=c(1.2,1.2,1.2,2),
#'        col=c(1,4,3,1), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L1", "L2", "sigma_min", "sigma_max")]
#' fit[-1]
#' head(summary(sdreport), 5)
#'
#' #############################################################################
#'
#' # Model 3: Fit to skipjack otoliths only
#'
#' init <- list(log_L1=log(30), log_L2=log(110),
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4)
#' model <- experiment(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "sigma_min", "sigma_max")]
#'
#' #############################################################################
#'
#' # Model 4: Fit to skipjack otoliths only,
#' # but now estimating constant sigma instead of sigma varying by length
#'
#' # We do this by omitting log_sigma_max
#' init <- list(log_L1=log(30), log_L2=log(110),
#'              log_sigma_min=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4)
#' model <- experiment(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "sigma_min")]
#'
#' #############################################################################
#'
#' # Model 5: Fit to skipjack tags only
#'
#' init <- list(log_L1=log(30), log_L2=log(110),
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4)
#' model <- experiment(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "sigma_min", "sigma_max")]
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#' @importFrom linear linear_curve linear_objfun
#'
#' @export

experiment <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  if(is.null(par$log_sigma_min))
    stop("'par' list must include 'log_sigma_min'")
  if(!is.null(par$log_sigma_max))
  if(is.null(data$t1))
    stop("'data' list must include 't1'")
  if(is.null(data$t2))
    stop("'data' list must include 't2'")
  MakeADFun(wrap(experiment_objfun, data=data), par, silent=silent, ...)
}

#' @rdname experiment
#'
#' @export

experiment_curve <- linear_curve

#' @rdname experiment
#'
#' @export

experiment_objfun <- linear_objfun
