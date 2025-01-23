#' Gompertz Growth Model
#'
#' Fit a Gompertz growth model to otoliths and/or tags, using the Schnute
#' parametrization.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param k growth coefficient.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{gompertz} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{gompertz_curve} and
#' \code{gompertz_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}, predicted length at age \code{t1}
#'   \item \code{log_L2}, predicted length at age \code{t2}
#'   \item \code{k}, growth coefficient
#'   \item \code{log_sigma_1}, growth variability at length \code{Lshort}
#'   \item \code{log_sigma_2} (*), growth variability at length \code{Llong}
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
#'   \item \code{t1}, age where predicted length is \code{L1}
#'   \item \code{t2}, age where predicted length is \code{L2}
#'   \item \code{Lshort} (*), length where sd(length) is \code{sigma_1}
#'   \item \code{Llong} (*), length where sd(length) is \code{sigma_2}
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only. The reference lengths
#' \code{Lshort} and \code{Llong} are only required when estimating growth
#' variability that varies with length, i.e., when the \code{log_sigma_2}
#' parameter is specified.
#'
#' @return
#' The \code{gompertz} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{gompertz_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{gompertz_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The Schnute parametrization used in \code{gompertz} reduces parameter
#' correlation and improves convergence reliability compared to the traditional
#' parametrization used in \code{\link{gompertzo}}. Therefore, the
#' \code{gompertz} parametrization can be recommended for general usage, as both
#' parametrizations produce the same growth curve. However, there can be some
#' use cases where the traditional parametrization (\code{Linf}, \code{k},
#' \code{tau}) is preferred over the Schnute parametrization (\code{L1},
#' \code{L2}, \code{k}).
#'
#' Gompertz is a special case of the Richards (1959) model, where \eqn{b=0}. If
#' the best model fit of a \code{\link{richards}} model to a particular dataset
#' involves a very small estimated value of \eqn{b}, then the \code{gompertz}
#' model offers a preferable parametrization, as it produces the same curve
#' using fewer parameters.
#'
#' The Gompertz (1825) growth model, as parametrized by Schnute (1981, Eq. 16)
#' predicts length at age as:
#'
#' \deqn{L ~=~ L_1\exp\!\left[\,\log(L_2/L_1)\,
#'       \frac{1-e^{-k(t-t_1)}}{1-e^{-k(t_2-t_1)}}\,\right]}{
#'       L1 * exp(log(L2/L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1))))}
#'
#' The variability of length at age increases linearly with length,
#'
#' \deqn{\sigma_L ~=~ \alpha \,+\, \beta \hat L}{
#'       sigma_L = alpha + beta * Lhat}
#'
#' where the slope is \eqn{\beta=(\sigma_2-\sigma_1) /
#' (L_\mathrm{long}-L_\mathrm{short})}{beta = (sigma_2-sigma_1) /
#' (Llong-Lshort)} and the intercept is \eqn{\alpha=\sigma_1 - \beta
#' L_\mathrm{short}}{alpha = sigma_1 - beta * Lshort}. Alternatively, growth
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
#' Gompertz, B. (1825).
#' On the nature of the function expressive of the law of human mortality, and
#' on a new mode of determining the value of Life Contingencies.
#' \emph{Philosophical Transactions of the Royal Society}, \bold{115}, 513-583.
# \doi{10.1098/rstl.1825.0026}.
#'
#' Schnute, J. (1981).
#' A versatile growth model with statistically stable parameters.
#' \emph{Canadian Journal of Fisheries and Aquatic Science}, \bold{38},
#' 1128-1140.
#' \doi{10.1139/f81-153}.
#'
#' @seealso
#' \code{\link{gcm}}, \code{gompertz}/\code{\link{gompertzo}},
#' \code{\link{richards}}/\code{\link{richards}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}}/\code{\link{vonberto}} are alternative growth models.
#'
#' \code{\link{otoliths_had}}, \code{\link{otoliths_skj}}, and
#' \code{\link{tags_skj}} are example datasets.
#'
#' \code{\link{fishgrowth-package}} gives an overview of the package.
#'
#' @examples
#' # Model 1: Fit to skipjack otoliths and tags
#'
#' # Explore initial parameter values
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' x <- seq(0, 4, 0.1)
#' points(lenRel~I(lenRel/60), tags_skj, col=4)
#' points(lenRec~I(lenRel/60+liberty), tags_skj, col=3)
#' lines(x, gompertz_curve(x, L1=25, L2=75, k=1.2, t1=0, t2=4), lty=2)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(25), log_L2=log(75), k=1.2,
#'              log_sigma_1=log(1), log_sigma_2=log(1),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len,
#'             Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4, Lshort=30, Llong=60)
#' gompertz_objfun(init, dat)
#'
#' # Fit model
#' model <- gompertz(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model)
#'
#' # Plot results
#' plot(len~age, otoliths_skj, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, gompertz_curve(x, L1, L2, k, t1, t2))
#' lines(x, Lhat, lwd=2)
#'
#' # Model summary
#' est <- report[c("L1", "L2", "k", "sigma_1", "sigma_2")]
#' est
#' fit[-1]
#' head(summary(sdreport), 5)
#'
#' #############################################################################
#'
#' # Model 2: Fit to skipjack otoliths only
#'
#' init <- list(log_L1=log(25), log_L2=log(75), k=1.2,
#'              log_sigma_1=log(1), log_sigma_2=log(1))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4,
#'             Lshort=30, Llong=60)
#' model <- gompertz(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "k", "sigma_1", "sigma_2")]
#'
#' #############################################################################
#'
#' # Model 3: Fit to skipjack otoliths only,
#' # but now estimating constant sigma instead of sigma varying by length
#'
#' # We do this by omitting log_sigma_2, Lshort, Llong
#' init <- list(log_L1=log(25), log_L2=log(75), k=1.2,
#'              log_sigma_1=log(1))
#' dat <- list(Aoto=otoliths_skj$age, Loto=otoliths_skj$len, t1=0, t2=4)
#' model <- gompertz(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'                     control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "k", "sigma_1")]
#'
#' #############################################################################
#'
#' # Model 4: Fit to skipjack tags only
#'
#' init <- list(log_L1=log(25), log_L2=log(75), k=1.2,
#'              log_sigma_1=log(1), log_sigma_2=log(1),
#'              log_age=log(tags_skj$lenRel/60))
#' dat <- list(Lrel=tags_skj$lenRel, Lrec=tags_skj$lenRec,
#'             liberty=tags_skj$liberty, t1=0, t2=4, Lshort=30, Llong=60)
#' model <- gompertz(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'                    control=list(eval.max=1e4, iter.max=1e4))
#' model$report()[c("L1", "L2", "k")]
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#'
#' @export

gompertz <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  if(is.null(par$log_sigma_1))
    stop("'par' list must include 'log_sigma_1'")
  if(!is.null(par$log_sigma_2))
  {
    if(is.null(data$Lshort))
      stop("'data' list must include 'Lshort' when 'log_sigma_2' is specified")
    if(is.null(data$Llong))
      stop("'data' list must include 'Llong' when 'log_sigma_2' is specified")
  }
  if(is.null(data$t1))
    stop("'data' list must include 't1'")
  if(is.null(data$t2))
    stop("'data' list must include 't2'")
  MakeADFun(wrap(gompertz_objfun, data=data), par, silent=silent, ...)
}

#' @rdname gompertz
#'
#' @export

gompertz_curve <- function(t, L1, L2, k, t1, t2)
{
  L1 * exp(log(L2/L1) * ((1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))))
}

#' @rdname gompertz
#'
#' @export

gompertz_objfun <- function(par, data)
{
  # Extract parameters
  L1 <- exp(par$log_L1)
  L2 <- exp(par$log_L2)
  k <- par$k
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- if(is.null(par$log_sigma_2)) NULL else exp(par$log_sigma_2)

  # Extract data
  t1 <- data$t1
  t2 <- data$t2
  Lshort <- data$Lshort
  Llong <- data$Llong

  # Calculate sigma coefficients (sigma = a + b*L)
  if(is.null(sigma_2))
  {
    sigma_slope <- 0  # if user did not pass log_sigma_2 then use constant sigma
    sigma_intercept <- sigma_1
  }
  else
  {
    sigma_slope <- (sigma_2 - sigma_1) / (Llong - Lshort)
    sigma_intercept <- sigma_1 - Lshort * sigma_slope
  }

  # Initialize likelihood
  nll <- 0

  # Report quantities of interest
  REPORT(L1)
  REPORT(L2)
  REPORT(k)
  REPORT(t1)
  REPORT(t2)
  REPORT(Lshort)
  REPORT(Llong)
  REPORT(sigma_1)
  REPORT(sigma_2)

  # Model includes otolith data
  if(!is.null(data$Aoto) && !is.null(data$Loto))
  {
    # data
    Aoto <- data$Aoto
    Loto <- data$Loto
    # Lhat
    Loto_hat <- gompertz_curve(Aoto, L1, L2, k, t1, t2)
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
    Lrel_hat <- gompertz_curve(age, L1, L2, k, t1, t2)
    Lrec_hat <- gompertz_curve(age+liberty, L1, L2, k, t1, t2)
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
