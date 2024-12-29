#' Von Bertalanffy Growth Model
#'
#' Fit a von Bertalanffy growth model to tags and otoliths.
#'
#' @param par a parameter list.
#' @param data a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param k growth coefficient.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param \dots passed to \code{MakeADFun}.
#'
#' @details
#' The main function \code{vonbert} creates a model object, ready for parameter
#' estimation. The auxiliary functions \code{vonbert_curve} and
#' \code{vonbert_objfun} are called by the main function to calculate the
#' regression curve and objective function value. The user can also call the
#' auxiliary functions directly for plotting and model exploration.
#'
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}, predicted length at age \code{t1}
#'   \item \code{log_L2}, predicted length at age \code{t2}
#'   \item \code{log_k}, growth coefficient
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
#'   \item \code{t1}, age where predicted length is \code{L1}
#'   \item \code{t2}, age where predicted length is \code{L2}
#'   \item \code{L_short}, length where sd(length) is \code{sigma_1}
#'   \item \code{L_long}, length where sd(length) is \code{sigma_2}
#' }
#'
#' @return
#' The \code{vonbert} function returns a TMB model object, produced by
#' \code{\link[RTMB]{MakeADFun}}.
#'
#' The \code{vonbert_curve} function returns a numeric vector of predicted
#' length at age.
#'
#' The \code{vonbert_objfun} function returns the negative log-likelihood as a
#' single number, describing the goodness of fit of \code{par} to the
#' \code{data}.
#'
#' @note
#' The von Bertalanffy (1938) growth model, as parametrized by Schnute and
#' Fournier (1980), predicts length at age as:
#'
#' \deqn{\hat L_t ~=~ L_1 \;+\; (L_2-L_1)\,
#'       \frac{1\,-\,e^{-k(t-t_1)}}{\,1\,-\,e^{-k(t_2-t_1)}\,}}{
#'       L1 + (L2-L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))}
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
#' von Bertalanffy, L. (1938).
#' A quantitative theory of organic growth.
#' \emph{Human Biology}, \bold{10}, 181-213.
#' \url{https://www.jstor.org/stable/41447359}.
#'
#' Schnute, J. and Fournier, D. (1980).
#' A new approach to length-frequency analysis: Growth structure.
#' \emph{Canadian Journal of Fisheries and Aquatic Science}, \bold{37},
#' 1337-1351.
#' \doi{10.1139/f80-172}.
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{richards}}, \code{\link{schnute3}}, and
#' \code{vonbert} are alternative growth models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @export

vonbert <- function(par, data, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  MakeADFun(wrap(vonbert_objfun, data=data), par, ...)
}

#' @rdname vonbert
#'
#' @export

vonbert_curve <- function(t, L1, L2, k, t1, t2)
{
  L1 + (L2-L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))
}

#' @rdname vonbert
#'
#' @export

vonbert_objfun <- function(par, data)
{
  # Extract parameters
  L1 <- exp(par$log_L1)
  L2 <- exp(par$log_L2)
  k <- exp(par$log_k)
  sigma_1 <- exp(par$log_sigma_1)
  sigma_2 <- exp(par$log_sigma_2)
  age <- exp(par$log_age)

  # Extract data
  Lrel <- data$Lrel
  Lrec <- data$Lrec
  liberty <- data$liberty
  Aoto <- data$Aoto
  Loto <- data$Loto
  t1 <- data$t1
  t2 <- data$t2
  L_short <- data$L_short
  L_long <- data$L_long

  # Calculate sigma coefficients (sigma = a + b*age)
  sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)
  sigma_intercept <- sigma_1 - L_short * sigma_slope

  # Calculate Lhat and sigma
  Lrel_hat <- vonbert_curve(age, L1, L2, k, t1, t2)
  Lrec_hat <- vonbert_curve(age+liberty, L1, L2, k, t1, t2)
  Loto_hat <- vonbert_curve(Aoto, L1, L2, k, t1, t2)
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
  curve <- vonbert_curve(age_seq, L1, L2, k, t1, t2)

  # Report quantities of interest
  REPORT(L1)
  REPORT(L2)
  REPORT(k)
  REPORT(age)
  REPORT(liberty)
  REPORT(Lrel)
  REPORT(Lrec)
  REPORT(Aoto)
  REPORT(Loto)
  REPORT(Lrel_hat)
  REPORT(Lrec_hat)
  REPORT(Loto_hat)
  REPORT(t1)
  REPORT(t2)
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
