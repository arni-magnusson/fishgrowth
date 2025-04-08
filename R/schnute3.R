#' Schnute Case 3 Model
#'
#' Fit a Schnute Case 3 model to otoliths and/or tags.
#'
#' @param par a parameter list.
#' @param data a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param b shape parameter.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{schnute3} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{schnute3_curve} and
#' \code{schnute3_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}, predicted length at age \code{t1}
#'   \item \code{log_L2}, predicted length at age \code{t2}
#'   \item \code{b}, shape parameter
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
#'   \item \code{t1}, age where predicted length is \code{L1}
#'   \item \code{t2}, age where predicted length is \code{L2}
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
#'
#' @return
#' The \code{schnute3} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{schnute3_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{schnute3_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The Schnute Case 3 model is a special case of the Richards (1959) model,
#' where \eqn{k=0}. If the best model fit of a \code{\link{richards}} model to a
#' particular dataset involves a very small estimated value of \eqn{k}, then the
#' \code{schnute3} model offers a preferable parametrization, as it produces the
#' same curve using fewer parameters.
#'
#' The Schnute Case 3 model (Schnute 1981, Eq. 17) predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ \left[\;L_1^b\;+\;(L_2^b-L_1^b)\,
#'       \frac{t-t_1}{t_2-t_1}\,\right]^{1/b}}{
#'       Lt = (L1^b + (L2^b-L1^b) * (t-t1) / (t2-t1))^(1/b)}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{
#'       sigma_L = alpha + beta * Lhat}
#'
#' where the slope is \eqn{\beta=(\sigma_{\max}-\sigma_{\min}) /
#' (L_{\max}-L_{\min})}{beta = (sigma_max-sigma_min) / (L_max-L_min)}, the
#' intercept is \eqn{\alpha=\sigma_{\min} - \beta L_{\min}}{alpha = sigma_min -
#' beta * L_min}, and \eqn{L_{\min}}{L_min} and \eqn{L_{\max}}{L_max} are the
#' shortest and longest observed lengths in the data. Alternatively, growth
#' variability can be modelled as a constant
#' \eqn{\sigma_L=\sigma_{\min}}{sigma_L=sigma_min} that does not vary with
#' length, by omitting \code{log_sigma_max} from the parameter list (see above).
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
#' Richards, F.J. (1959).
#' A flexible growth function for empirical use.
#' \emph{Journal of Experimental Botany}, \bold{10}, 290-300.
#' \url{https://www.jstor.org/stable/23686557}.
#'
#' Schnute, J. (1981).
#' A versatile growth model with statistically stable parameters.
#' \emph{Canadian Journal of Fisheries and Aquatic Science}, \bold{38},
#' 1128-1140.
#' \doi{10.1139/f81-153}.
#'
#' The \code{\link{fishgrowth-package}} help page includes references describing
#' the parameter estimation method.
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}, \code{\link{gompertzo}},
#' \code{\link{richards}}, \code{\link{richardso}}, \code{schnute3},
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
#' lines(x, schnute3_curve(x, L1=15, L2=70, b=2, t1=1, t2=10), lty=3)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(15), log_L2=log(70), b=2,
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_had$age, Loto=otoliths_had$len, t1=1, t2=10)
#' schnute3_objfun(init, dat)
#'
#' # Fit model
#' model <- schnute3(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' Lhat <- with(report, schnute3_curve(x, L1, L2, b, t1, t2))
#' lines(x, Lhat, lwd=2, col=2)
#' legend("bottomright", c("initial curve","model fit"), col=c(1,2), lty=c(3,1),
#'        lwd=c(1,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L1", "L2", "b", "sigma_min", "sigma_max")]
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
#' lines(x, schnute3_curve(x, L1=25, L2=75, b=3, t1=0, t2=4), lty=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "initial curve"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,2),
#'        lwd=c(1.2,1.2,1.2,1), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(25), log_L2=log(75), b=3,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len,
#'             Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4)
#' schnute3_objfun(init, dat)
#'
#' # Fit model
#' model <- schnute3(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, schnute3_curve(x, L1, L2, b, t1, t2))
#' lines(x, Lhat, lwd=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "model fit"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,1),
#'        lwd=c(1.2,1.2,1.2,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("L1", "L2", "b", "sigma_min", "sigma_max")]
#' fit[-1]
#' head(summary(sdreport), 5)
#'
#' #############################################################################
#'
#' # Model 3: Fit to skipjack otoliths only
#'
#' init <- list(log_L1=log(25), log_L2=log(75), b=3,
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4)
#' model <- schnute3(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "b", "sigma_min", "sigma_max")]
#'
#' #############################################################################
#'
#' # Model 4: Fit to skipjack otoliths only,
#' # but now estimating constant sigma instead of sigma varying by length
#'
#' # We do this by omitting log_sigma_max
#' init <- list(log_L1=log(25), log_L2=log(75), b=3,
#'              log_sigma_min=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4)
#' model <- schnute3(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "b", "sigma_min")]
#'
#' #############################################################################
#'
#' # Model 5: Fit to skipjack tags only
#'
#' init <- list(log_L1=log(25), log_L2=log(75), b=3,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4)
#' model <- schnute3(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "b", "sigma_min", "sigma_max")]
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#'
#' @export

schnute3 <- function(par, data, silent=TRUE, ...)
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
  MakeADFun(wrap(schnute3_objfun, data=data), par, silent=silent, ...)
}

#' @rdname schnute3
#'
#' @export

schnute3_curve <- function(t, L1, L2, b, t1, t2)
{
  (L1^b + (L2^b-L1^b) * (t-t1) / (t2-t1))^(1/b)
}

#' @rdname schnute3
#'
#' @export

schnute3_objfun <- function(par, data)
{
  # Extract parameters
  L1 <- exp(par$log_L1)
  L2 <- exp(par$log_L2)
  b <- par$b
  sigma_min <- exp(par$log_sigma_min)
  sigma_max <- if(is.null(par$log_sigma_max)) NULL else exp(par$log_sigma_max)

  # Extract data
  t1 <- data$t1
  t2 <- data$t2

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
  type <- "schnute3"
  curve <- schnute3_curve
  REPORT(type)
  REPORT(curve)
  REPORT(L1)
  REPORT(L2)
  REPORT(b)
  REPORT(t1)
  REPORT(t2)
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
    Loto_hat <- schnute3_curve(Aoto, L1, L2, b, t1, t2)
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
    Lrel_hat <- schnute3_curve(age, L1, L2, b, t1, t2)
    Lrec_hat <- schnute3_curve(age+liberty, L1, L2, b, t1, t2)
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
