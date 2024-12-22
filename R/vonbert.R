#' von Bertalanffy
#'
#' Fit a von Bertalanffy growth model to tags and otoliths.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#'
#' @details
#' The \code{par} list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}, predicted length at age \code{t1}
#'   \item \code{log_L2}, predicted length at age \code{t2}
#'   \item \code{log_k}, growth coefficient
#'   \item \code{log_sigma_1}, growth variability at length \code{L_short}
#'   \item \code{log_sigma_2}, growth variability at length \code{L_long}
#'   \item \code{log_age} predicted age at release (vector)
#' }
#'
#' The \code{data} list contains the following elements:
#' \itemize{
#'   \item \code{Lrel}, length at release (vector)
#'   \item \code{Lrec}, length at recapture (vector)
#'   \item \code{liberty}, time at liberty in years (vector)
#'   \item \code{Aoto}, age from otoliths (vector)
#'   \item \code{Loto}, length from otoliths (vector)
#'   \item \code{t1}, age where predicted length is \code{L1}
#'   \item \code{t2}, age where predicted length is \code{L2}
#'   \item \code{L_short}, length where sd(length) is \code{sigma_1}
#'   \item \code{L_long}, length where sd(length) is \code{sigma_2}
#' }
#'
#' @return
#' TMB model object, produced by \code{\link[RTMB]{MakeADFun}}.
#'
#' @references
#' von Bertalanffy, L. (1938).
#' A quantitative theory of organic growth.
#' \emph{Human Biology}, \bold{10}, 181-213.
#' \url{https://www.jstor.org/stable/41447359}.
#'
#' Schnute, J. and Fournier, D.. 1980.
#' A new approach to length-frequence analysis: Growth structure.
#' \emph{Canadian Journal of Fisheries and Aquatic Science}, \bold{37},
#' 1337-1351.
#' \doi{10.1139/f80-172}.
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
#'
#' @seealso
#' \code{\link{otoliths_skj}} and \code{\link{tags_skj}} are example datasets.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @export

vonbert <- function(par, data)
{
  wrap <- function(objfun, ...) function(par) objfun(par, ...)
  MakeADFun(wrap(vb_objfun, data=data), par, silent=TRUE)
}

vb_objfun <- function(par, data)
{
  # Extract parameters
  log_L1 <- par$log_L1
  log_L2 <- par$log_L2
  log_k <- par$log_k
  log_sigma_1 <- par$log_sigma_1
  log_sigma_2 <- par$log_sigma_2
  log_age <- par$log_age

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

  # Calculate parameters
  L1 <- exp(log_L1)
  L2 <- exp(log_L2)
  k <- exp(log_k)
  sigma_1 <- exp(log_sigma_1)
  sigma_2 <- exp(log_sigma_2)
  sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)  # s <- a + b*age
  sigma_intercept <- sigma_1 - L_short * sigma_slope
  age <- exp(log_age)

  # Calculate Lhat and sigma
  Lrel_hat <- L1 + (L2-L1) * (1-exp(-k*(age-t1))) / (1-exp(-k*(t2-t1)))
  Lrec_hat <- L1 + (L2-L1) * (1-exp(-k*(age+liberty-t1))) / (1-exp(-k*(t2-t1)))
  Loto_hat <- L1 + (L2-L1) * (1-exp(-k*(Aoto-t1))) / (1-exp(-k*(t2-t1)))
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
  curve <- L1 + (L2-L1) * (1-exp(-k*(age_seq-t1))) / (1-exp(-k*(t2-t1)))

  # Report quantities of interest
  ADREPORT(curve)
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
  REPORT(curve)

  nll
}
