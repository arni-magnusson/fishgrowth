#' Schnute Case 3 Model
#'
#' Fit a Schnute Case 3 model to tags and otoliths.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#' @param t age (vector).
#' @param L1 predicted length at age \code{t1}.
#' @param L2 predicted length at age \code{t2}.
#' @param b shape parameter.
#' @param t1 age where predicted length is \code{L1}.
#' @param t2 age where predicted length is \code{L2}.
#' @param \dots passed to \code{MakeADFun}.
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
#' \deqn{L ~=~ \left[\;L_1^b\;+\;(L_2^b-L_1^b)\,
#'       \frac{t-t_1}{t_2-t_1}\,\right]^{1/b}}{
#'       (L1^b + (L2^b-L1^b) * (t-t1) / (t2-t1))^(1/b)}
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
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @seealso
#' \code{\link{gcm}}, \code{\link{richards}}, \code{schnute3}, and
#' \code{\link{vonbert}} are alternative growth models.
#'
#' \code{\link{otoliths_ex}} and \code{\link{tags_ex}} are example datasets.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @export

schnute3 <- function(par, data, ...)
{
  wrap <- function(objfun, data)
  {
    function(par) objfun(par, data)
  }
  MakeADFun(wrap(schnute3_objfun, data=data), par, ...)
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
  Lrel_hat <- schnute3_curve(age, L1, L2, b, t1, t2)
  Lrec_hat <- schnute3_curve(age+liberty, L1, L2, b, t1, t2)
  Loto_hat <- schnute3_curve(Aoto, L1, L2, b, t1, t2)
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
  curve <- schnute3_curve(age_seq, L1, L2, b, t1, t2)

  # Report quantities of interest
  REPORT(L1)
  REPORT(L2)
  REPORT(b)
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
