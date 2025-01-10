#' @name fishgrowth-package
#'
#' @aliases fishgrowth
#'
#' @title Fit Growth Curves to Fish Data
#'
#' @description
#' Fit growth models to otolith data and/or tagging data data, using RTMB and
#' maximum likelihood. The otoliths provide direct observed coordinates of age
#' and length. The tagging data provide information about the observed length at
#' release and length at recapture at a later time, where the age at release is
#' unknown and estimated as a vector of parameters.
#'
#' @details
#' \emph{Growth models:}
#' \tabular{ll}{
#'   \code{\link{gcm}}       \tab growth cessation\cr
#'   \code{\link{gompertz}}  \tab Gompertz\cr
#'   \code{\link{gompertzo}} \tab Gompertz (old style)\cr
#'   \code{\link{richards}}  \tab Richards\cr
#'   \code{\link{richardso}} \tab Richards (old style)\cr
#'   \code{\link{schnute3}}  \tab Schnute Case 3\cr
#'   \code{\link{vonbert}}   \tab von Bertalanffy\cr
#'   \code{\link{vonberto}}  \tab von Bertalanffy (old style)
#' }
#' \emph{Example data:}
#' \tabular{ll}{
#'   \code{\link{otoliths_ex}} \tab otoliths\cr
#'   \code{\link{tags_ex}}     \tab tags
#' }
#'
#' @author Arni Magnusson and Mark Maunder.

"_PACKAGE"
