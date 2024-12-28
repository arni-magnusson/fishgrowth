#' @name tao-package
#'
#' @aliases tao
#'
#' @title Growth from Tags and Otoliths
#'
#' @description
#' Fit growth models to tagging data and otolith data, using RTMB and maximum
#' likelihood. The tagging data provide information about the observed length at
#' release and length at recapture at a later time, where the age at release is
#' unknown and estimated as a vector of parameters. The otoliths provide direct
#' observed coordinates of age and length.
#'
#' @details
#' \emph{Growth models:}
#' \tabular{ll}{
#'   \code{\link{gcm}}      \tab growth cessation\cr
#'   \code{\link{richards}} \tab Richards\cr
#'   \code{\link{schnute3}} \tab Schnute Case 3\cr
#'   \code{\link{vonbert}}  \tab von Bertalanffy
#' }
#' \emph{Example data:}
#' \tabular{ll}{
#'   \code{\link{otoliths_ex}} \tab otolith data\cr
#'   \code{\link{tags_ex}}     \tab tagging data
#' }
#'
#' @author Arni Magnusson and Mark Maunder.

"_PACKAGE"
