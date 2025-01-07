#' Richards Growth Model
#'
#' Fit a Richards growth model to otoliths and/or tags.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param k growth coefficient.
#' @param b shape parameter.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param silent passed to \code{\link[RTMB]{MakeADFun}}.
#' @param \dots passed to \code{\link[RTMB]{MakeADFun}}.
#'
#' @details
#' The main function \code{richards} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{richards_curve} and
#' \code{richards_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}, predicted length at age \code{t1}
#'   \item \code{log_L2}, predicted length at age \code{t2}
#'   \item \code{log_k}, growth coefficient
#'   \item \code{b}, shape parameter
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
#'   \item \code{t1}, age where predicted length is \code{L1}
#'   \item \code{t2}, age where predicted length is \code{L2}
#'   \item \code{L_short}, length where sd(length) is \code{sigma_1}
#'   \item \code{L_long}, length where sd(length) is \code{sigma_2}
#' }
#'
#' *: The data vectors \code{Aoto} and \code{Loto} can be omitted to fit to
#' tagging data only. The data vectors \code{Lrel}, \code{Lrec}, and
#' \code{liberty} can be omitted to fit to otoliths only.
#'
#' @return
#' The \code{richards} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{richards_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{richards_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The Richards (1959) growth model, as parametrized by Schnute (1981, Eq. 15),
#' predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ \left[\:L_1^b\;+\;(L_2^b-L_1^b)\,
#'       \frac{1-e^{-k(t-t_1)}}{1-e^{-k(t_2-t_1)}}\,\right]^{1/b}}{
#'       (L1^b + (L2^b-L1^b) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1))))^(1/b)}
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
#' @seealso
#' \code{\link{gcm}}, \code{\link{gompertz}}, \code{richards},
#' \code{\link{schnute3}}, and \code{\link{vonbert}}/\code{\link{vonberto}} are
#' alternative growth models.
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
#' lines(x, richards_curve(x, L1=25, L2=75, k=0.8, b=1, t1=0, t2=4), lty=2)
#'
#' # Prepare parameters and data
#' init <- list(log_L1=log(25), log_L2=log(75), log_k=log(0.8), b=1,
#'              log_sigma_1=log(1), log_sigma_2=log(1),
#'              log_age=log(tags_ex$lenRel/60))
#' dat <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len,
#'             Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'             liberty=tags_ex$liberty, t1=0, t2=4, L_short=30, L_long=60)
#' richards_objfun(init, dat)
#'
#' # Fit model
#' model <- richards(init, dat)
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
#' est <- report[c("L1", "L2", "k", "b", "sigma_1", "sigma_2")]
#' est
#' fit[-1]
#' head(summary(sdreport), 6)
#'
#' #############################################################################
#'
#' # Fit to otoliths only
#' init_oto <- list(log_L1=log(25), log_L2=log(75), log_k=log(0.8), b=1,
#'                  log_sigma_1=log(1), log_sigma_2=log(1))
#' dat_oto <- list(Aoto=otoliths_ex$age, Loto=otoliths_ex$len, t1=0, t2=4,
#'                 L_short=30, L_long=60)
#' model_oto <- richards(init_oto, dat_oto)
#' fit_oto <- nlminb(model_oto$par, model_oto$fn, model_oto$gr,
#'                   control=list(eval.max=1e4, iter.max=1e4))
#' model_oto$report()[c("L1", "L2", "k", "b")]
#'
#' #############################################################################
#'
#' # Fit to tags only
#' init_tags <- list(log_L1=log(25), log_L2=log(75), log_k=log(0.8), b=1,
#'                   log_sigma_1=log(1), log_sigma_2=log(1),
#'                   log_age=log(tags_ex$lenRel/60))
#' dat_tags <- list(Lrel=tags_ex$lenRel, Lrec=tags_ex$lenRec,
#'                  liberty=tags_ex$liberty, t1=0, t2=4, L_short=30, L_long=60)
#' model_tags <- richards(init_tags, dat_tags)
#' fit_tags <- nlminb(model_tags$par, model_tags$fn, model_tags$gr,
#'                    control=list(eval.max=1e4, iter.max=1e4))
#' model_tags$report()[c("L1", "L2", "k", "b")]
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @export

richards <- function(par, data, silent=TRUE, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  MakeADFun(wrap(richards_objfun, data=data), par, silent=silent, ...)
}

#' @rdname richards
#'
#' @export

richards_curve <- function(t, L1, L2, k, b, t1, t2)
{
  (L1^b + (L2^b-L1^b) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1))))^(1/b)
}

#' @rdname richards
#'
#' @export

richards_objfun <- function(par, data)
{
  # Extract parameters
  L1 <- exp(par$log_L1)
  L2 <- exp(par$log_L2)
  k <- exp(par$log_k)
  b <- par$b
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- if(is.null(par$log_sigma_2)) NULL else exp(par$log_sigma_2)

  # Extract data
  t1 <- data$t1
  t2 <- data$t2
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
  curve <- richards_curve(age_seq, L1, L2, k, b, t1, t2)

  # Report quantities of interest
  REPORT(L1)
  REPORT(L2)
  REPORT(k)
  REPORT(b)
  REPORT(t1)
  REPORT(t2)
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
    Loto_hat <- richards_curve(Aoto, L1, L2, k, b, t1, t2)
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
    Lrel_hat <- richards_curve(age, L1, L2, k, b, t1, t2)
    Lrec_hat <- richards_curve(age+liberty, L1, L2, k, b, t1, t2)
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
