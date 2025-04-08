#' Richards Growth Model (Old Style)
#'
#' Fit a Richards growth model to otoliths and/or tags, using a traditional
#' parametrization.
#'
#' @param par a parameter list.
#' @param data a data list.
#' @param t age (vector).
#' @param Linf asymptotic maximum length.
#' @param k growth coefficient.
#' @param tau location parameter.
#' @param b shape parameter.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{richardso} creates a model object, ready for
#' parameter estimation. The auxiliary functions \code{richardso_curve} and
#' \code{richardso_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_Linf}, asymptotic maximum length
#'   \item \code{log_k}, growth coefficient
#'   \item \code{tau}, location parameter
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
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
#'
#' @return
#' The \code{richardso} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{richardso_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{richardso_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The Schnute parametrization used in \code{\link{richards}} reduces parameter
#' correlation and improves convergence reliability compared to the traditional
#' parametrization used in \code{richardso}. Therefore, the
#' \code{richards} parametrization can be recommended for general usage, as both
#' parametrizations produce the same growth curve. However, there can be some
#' use cases where the traditional parametrization (\code{Linf}, \code{k},
#' \code{tau}, \code{b}) is preferred over the Schnute parametrization
#' (\code{L1}, \code{L2}, \code{k}, \code{b}).
#'
#' The Richards (1959) growth model, as parametrized by Tjørve and Tjørve (2010,
#' Eq. 4), predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ L_\infty\left(1\,-\,\frac{1}{b}\,
#'       e^{-k(t-\tau)}\right)^{\!b}}{
#'       Lt = Linf * (1 - (1/b) * exp(-k*(t-tau)))^b}
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
#' Tjørve, E. and Tjørve, K.M.C. (2010).
#' A unified approach to the Richards-model family for use in growth analyses:
#' Why we need only two model forms.
#' \emph{Journal of Theoretical Biology}, \bold{267}, 417-425.
#'
#' The \code{\link{fishgrowth-package}} help page includes references describing
#' the parameter estimation method.
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}, \code{\link{gompertzo}},
#' \code{\link{richards}}, \code{richardso}, \code{\link{schnute3}},
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
#' lines(x, richardso_curve(x, Linf=100, k=0.1, tau=-1, b=1), lty=3)
#'
#' # Prepare parameters and data
#' init <- list(log_Linf=log(100), log_k=log(0.1), tau=-1, b=1,
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_had$age, Loto=otoliths_had$len)
#' richardso_objfun(init, dat)
#'
#' # Fit model
#' model <- richardso(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' Lhat <- with(report, richardso_curve(x, Linf, k, tau, b))
#' lines(x, Lhat, lwd=2, col=2)
#' legend("bottomright", c("initial curve","model fit"), col=c(1,2), lty=c(3,1),
#'        lwd=c(1,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("Linf", "k", "tau", "b", "sigma_min", "sigma_max")]
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
#' lines(x, richardso_curve(x, Linf=80, k=0.8, tau=-0.5, b=1), lty=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "initial curve"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,2),
#'        lwd=c(1.2,1.2,1.2,1), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Prepare parameters and data
#' init <- list(log_Linf=log(80), log_k=log(0.8), tau=-0.5, b=1,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len,
#'             Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty)
#' richardso_objfun(init, dat)
#'
#' # Fit model
#' model <- richardso(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, richardso_curve(x, Linf, k, tau, b))
#' lines(x, Lhat, lwd=2)
#' legend("bottomright", c("otoliths","tag releases","tac recaptures",
#'        "model fit"), col=c(1,4,3,1), pch=c(1,1,1,NA), lty=c(0,0,0,1),
#'        lwd=c(1.2,1.2,1.2,2), bty="n", inset=0.02, y.intersp=1.25)
#'
#' # Model summary
#' report[c("Linf", "k", "tau", "b", "sigma_min", "sigma_max")]
#' fit[-1]
#' head(summary(sdreport), 6)
#'
#' #############################################################################
#'
#' # Model 3: Fit to skipjack otoliths only
#'
#' init <- list(log_Linf=log(80), log_k=log(0.8), tau=-0.5, b=1,
#'              log_sigma_min=log(3), log_sigma_max=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len)
#' model <- richardso(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("Linf", "k", "tau", "b", "sigma_min", "sigma_max")]
#'
#' #############################################################################
#'
#' # Model 4: Fit to skipjack otoliths only,
#' # but now estimating constant sigma instead of sigma varying by length
#'
#' # We do this by omitting log_sigma_max
#' init <- list(log_Linf=log(80), log_k=log(0.8), tau=-0.5, b=1,
#'              log_sigma_min=log(3))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len)
#' model <- richardso(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("Linf", "k", "tau", "b", "sigma_min")]
#'
#' #############################################################################
#'
#' # Model 5: Fit to skipjack tags only
#'
#' init <- list(log_Linf=log(80), log_k=log(0.8), tau=-0.5, b=1,
#'              log_sigma_min=log(3), log_sigma_max=log(3),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty)
#' model <- richardso(init, dat)           # using 1e3 to keep CRAN checks fast,
#' fit <- nlminb(model$par, model$fn, model$gr,             # but try 1e4 to get
#'               control=list(eval.max=1e3, iter.max=1e3))  # better convergence
#' model$report()[c("Linf", "k", "tau", "b", "sigma_min", "sigma_max")]
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#'
#' @export

richardso <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  if(is.null(par$log_sigma_min))
    stop("'par' list must include 'log_sigma_min'")
  MakeADFun(wrap(richardso_objfun, data=data), par, silent=silent, ...)
}

#' @rdname richardso
#'
#' @export

richardso_curve <- function(t, Linf, k, tau, b)
{
  Linf * (1 - (1/b) * exp(-k*(t-tau)))^b
}

#' @rdname richardso
#'
#' @export

richardso_objfun <- function(par, data)
{
  # Extract parameters
  Linf <- exp(par$log_Linf)
  k <- exp(par$log_k)
  tau <- par$tau
  b <- par$b
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
  type <- "richardso"
  curve <- richardso_curve
  REPORT(type)
  REPORT(curve)
  REPORT(Linf)
  REPORT(k)
  REPORT(tau)
  REPORT(b)
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
    Loto_hat <- richardso_curve(Aoto, Linf, k, tau, b)
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
    Lrel_hat <- richardso_curve(age, Linf, k, tau, b)
    Lrec_hat <- richardso_curve(age+liberty, Linf, k, tau, b)
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
