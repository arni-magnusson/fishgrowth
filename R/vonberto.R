#' Von Bertalanffy Growth Model (Old Style)
#'
#' Fit a von Bertalanffy growth model to otoliths and/or tags, using a
#' traditional parametrization.
#'
#' @param par a parameter list.
#' @param data a data list.
#' @param t age (vector).
#' @param Linf asymptotic maximum length.
#' @param k growth coefficient.
#' @param t0 age where the predicted length is zero, the x-intercept.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{vonberto} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{vonberto_curve} and
#' \code{vonberto_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_Linf}, asymptotic maximum length
#'   \item \code{log_k}, growth coefficient
#'   \item \code{to}, age where the predicted length is zero, the x-intercept
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
#'   \item \code{Lshort}, length where sd(length) is \code{sigma_1}
#'   \item \code{Llong}, length where sd(length) is \code{sigma_2}
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
#'
#' @return
#' The \code{vonberto} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{vonberto_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{vonberto_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The Schnute-Fournier parametrization used in \code{\link{vonbert}} reduces
#' parameter correlation and improves convergence reliability compared to the
#' traditional parametrization used in \code{vonberto}. Therefore, the
#' \code{vonbert} parametrization can be recommended for general usage, as both
#' parametrizations produce the same growth curve. However, there can be some
#' use cases where the traditional parametrization (\code{Linf}, \code{k},
#' \code{t0}) is preferred over the Schnute-Fournier parametrization (\code{L1},
#' \code{L2}, \code{k}).
#'
#' The von Bertalanffy (1938) growth model, as parametrized by Beverton and Holt
#' (1957), predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ L_\infty\left(1\,-\,e^{-k(t-t_0)}\right)}{
#'       Linf * (1 - exp(-k*(t-t0)))}
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
#' von Bertalanffy, L. (1938).
#' A quantitative theory of organic growth.
#' \emph{Human Biology}, \bold{10}, 181-213.
#' \url{https://www.jstor.org/stable/41447359}.
#'
#' Beverton, R.J.H. and Holt, S.J. (1957).
#' \emph{On the dynamics of exploited fish populations}.
#' London: Her Majesty’s Stationery Office.
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}/\code{\link{gompertzo}},
#' \code{\link{richards}}/\code{\link{richards}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}}/\code{vonberto} are alternative growth models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{fishgrowth-package}} gives an overview of the package.
#'
#' @examples
#' # Explore initial parameter values
#' plot(len~age, otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' x <- seq(0, 4, 0.1)
#' points(lenRel~I(lenRel/60), tags_ex, col=4)
#' points(lenRec~I(lenRel/60+liberty), tags_ex, col=3)
#' lines(x, vonberto_curve(x, Linf=80, k=0.8, t0=-0.5), lty=2)
#'
#' # Prepare parameters and data
#' init <- list(log_Linf=log(80), log_k=log(0.8), t0=-0.5,
#'              log_sigma_1=log(1), log_sigma_2=log(1),
#'              log_age=log(tags_ex$lenRel/60))
#' dat <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len,
#'             Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'             liberty=tags_ex$liberty, Lshort=30, Llong=60)
#' vonberto_objfun(init, dat)
#'
#' # Fit model
#' model <- vonberto(init, dat)
#' fit <- nlminb(model$par, model$fn, model$gr,
#'               control=list(eval.max=1e4, iter.max=1e4))
#' report <- model$report()
#' sdreport <- sdreport(model, getReportCovariance=FALSE)
#'
#' # Plot results
#' plot(len~age, otoliths_ex, xlim=c(0,4), ylim=c(0,100))
#' points(report$age, report$Lrel, col=4)
#' points(report$age+report$liberty, report$Lrec, col=3)
#' Lhat <- with(report, vonberto_curve(x, Linf, k, t0))
#' lines(x, Lhat, lwd=2)
#'
#' # Model summary
#' est <- report[c("Linf", "k", "t0", "sigma_1", "sigma_2")]
#' est
#' fit[-1]
#' head(summary(sdreport), 5)
#'
#' #############################################################################
#'
#' # Fit to otoliths only
#' init_oto <- list(log_Linf=log(80), log_k=log(0.8), t0=-0.5,
#'                  log_sigma_1=log(1), log_sigma_2=log(1))
#' dat_oto <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len,
#'                 Lshort=30, Llong=60)
#' model_oto <- vonberto(init_oto, dat_oto)
#' fit_oto <- nlminb(model_oto$par, model_oto$fn, model_oto$gr,
#'                   control=list(eval.max=1e4, iter.max=1e4))
#' model_oto$report()[c("Linf", "k", "t0")]
#'
#' #############################################################################
#'
#' # Fit to tags only
#' init_tags <- list(log_Linf=log(80), log_k=log(0.8), t0=-0.5,
#'                   log_sigma_1=log(1), log_sigma_2=log(1),
#'                   log_age=log(tags_ex$lenRel/60))
#' dat_tags <- list(Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'                  liberty=tags_ex$liberty, Lshort=30, Llong=60)
#' model_tags <- vonberto(init_tags, dat_tags)
#' fit_tags <- nlminb(model_tags$par, model_tags$fn, model_tags$gr,
#'                    control=list(eval.max=1e4, iter.max=1e4))
#' model_tags$report()[c("Linf", "k", "t0")]
#'
#' @importFrom RTMB dnorm MakeADFun REPORT
#'
#' @export

vonberto <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  MakeADFun(wrap(vonberto_objfun, data=data), par, silent=silent, ...)
}

#' @rdname vonberto
#'
#' @export

vonberto_curve <- function(t, Linf, k, t0)
{
  Linf * (1 - exp(-k*(t-t0)))
}

#' @rdname vonberto
#'
#' @export

vonberto_objfun <- function(par, data)
{
  # Extract parameters
  Linf <- exp(par$log_Linf)
  k <- exp(par$log_k)
  t0 <- par$t0
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- if(is.null(par$log_sigma_2)) NULL else exp(par$log_sigma_2)

  # Extract data
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
  REPORT(Linf)
  REPORT(k)
  REPORT(t0)
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
    Loto_hat <- vonberto_curve(Aoto, Linf, k, t0)
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
    Lrel_hat <- vonberto_curve(age, Linf, k, t0)
    Lrec_hat <- vonberto_curve(age+liberty, Linf, k, t0)
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
