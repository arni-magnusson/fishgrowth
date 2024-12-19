vonbert <- function(par, data)
{
  wrap <- function(objfun, ...) function(par) objfun(par, ...)
  MakeADFun(wrap(vb_objfun, data=data), par, silent=TRUE)
}

vb_objfun <- function(par, data)
{
  # Extract data and parameters
  getAll(data, warn=FALSE)
  getAll(par)

  # Calculate parameters
  Linf <- exp(log_Linf)
  k <- exp(log_k)
  sigma_1 <- exp(log_sigma_1)
  sigma_2 <- exp(log_sigma_2)
  sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)  # s <- a + b*age
  sigma_intercept <- sigma_1 - L_short * sigma_slope
  age <- exp(log_age)

  # Calculate Lhat and sigma
  Lrel_hat <- Linf * (1 - exp(-k * (age - t0)))
  Lrec_hat <- Linf * (1 - exp(-k * (age + liberty - t0)))
  Loto_hat <- Linf * (1 - exp(-k * (Aoto - t0)))
  sigma_Lrel <- sigma_intercept + sigma_slope * Lrel_hat
  sigma_Lrec <- sigma_intercept + sigma_slope * Lrec_hat
  sigma_Loto <- sigma_intercept + sigma_slope * Loto_hat

  # Calculate likelihoods
  nll_Lrel <- -dnorm(Lrel, Lrel_hat, sigma_Lrel, TRUE)
  nll_Lrec <- -dnorm(Lrec, Lrec_hat, sigma_Lrec, TRUE)
  nll_Loto <- -dnorm(Loto, Loto_hat, sigma_Loto, TRUE)
  nll <- sum(nll_Lrel) + sum(nll_Lrec) + sum(nll_Loto)

  # Calculate curve
  age_seq = seq(0, 10, 1/365)  # age 1-10 years, day by day
  curve <- Linf * (1 - exp(-k * (age_seq - t0)))

  # Report quantities of interest
  ADREPORT(curve)
  REPORT(Linf)
  REPORT(k)
  REPORT(t0)
  REPORT(age)
  REPORT(liberty)
  REPORT(Lrel)
  REPORT(Lrec)
  REPORT(Aoto)
  REPORT(Loto)
  REPORT(Lrel_hat)
  REPORT(Lrec_hat)
  REPORT(Loto_hat)
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
