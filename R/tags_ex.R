#' @docType data
#'
#' @name tags_ex
#'
#' @title Tagging Data
#'
#' @description
#' Simulated tagging data, loosely based on a skipjack tuna dataset analyzed by
#' Macdonald et al. (2022).
#'
#' @usage
#' tags_ex
#'
#' @format
#' Data frame containing three columns:
#' \tabular{ll}{
#'   \code{lenRel}  \tab length at release (cm)\cr
#'   \code{lenRec}  \tab length at recapture (cm)\cr
#'   \code{liberty} \tab time at liberty (years)
#' }
#'
#' @details
#' The simulation code that was used to produce this dataset is included in the
#' package:
#' \preformatted{file.show(system.file(package="tao", "sim/simulate.R"))}
#'
#' @source
#' Macdonald, J., Day, J., Magnusson, A., Maunder, M., Aoki, Y., Matsubara, N.,
#' Tsuda, Y., McKechnie, S., Teears, T., Leroy, B., Castillo-Jordán, C.,
#' Hampton, J., and Hamer, P. (2022).
#' \emph{Review and new analyses of skipjack growth in the Western and Central
#' Pacific Ocean}.
#' Western and Central Pacific Fisheries Commission Report
#' WCPFC-SC18-2022/SA-IP-06.
#' \url{https://meetings.wcpfc.int/node/16254}.
#'
#' @seealso
#' \code{\link{otoliths_ex}} and \code{tags_ex} are example datasets.
#'
#' \code{\link{gcm}}, \code{\link{richards}}, \code{\link{schnute3}}, and
#' \code{\link{vonbert}} are alternative growth models.
#'
#' \code{\link{tao-package}} gives an overview of the package.
#'
#' @examples
#' tags_ex

NA
