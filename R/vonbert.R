#' von Bertalanffy
#'
#' Fit a von Bertalanffy growth model to tags and otoliths.
#'
#' @param par is a parameter list.
#' @param data is a data list.
#'
#' @details
#' The parameter list contains the following elements:
#' \itemize{
#'   \item \code{log_L1}
#'   \item \code{log_L2}
#'   \item \code{log_k}
#'   \item \code{log_sigma_1}
#'   \item \code{log_sigma_2}
#'   \item \code{log_age} (vector)
#' }
#'
#' The data list contains the following elements:
#' \itemize{
#'   \item \code{Lrel} (vector)
#'   \item \code{Lrec} (vector)
#'   \item \code{liberty} (vector)
#'   \item \code{Aoto} (vector)
#'   \item \code{Loto} (vector)
#'   \item \code{t1}
#'   \item \code{t2}
#'   \item \code{L_short}
#'   \item \code{L_long}
#' }
#'
#' @return
#' List produced by \code{\link[RTMB]{MakeADFun}}.
#'
#' @importFrom RTMB ADREPORT dnorm MakeADFun REPORT
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
