#' @name fishgrowth-package
#'
#' @aliases fishgrowth
#'
#' @title Fit Growth Curves to Fish Data
#'
#' @description
#' Fit growth models to otoliths and/or tagging data, using the
#' \code{\link[RTMB]{RTMB}} package and maximum likelihood.
#'
#' The otoliths (or similar measurements of age) provide direct observed
#' coordinates of age and length. The tagging data provide information about the
#' observed length at release and length at recapture at a later time, where the
#' age at release is unknown and estimated as a vector of parameters.
#'
#' The growth models provided by this package can be fitted to otoliths only,
#' tagging data only, or a combination of the two. Growth variability can be
#' modelled as constant or increasing with length.
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
#' \emph{Utilities:}
#' \tabular{ll}{
#'   \code{\link{pred_band}}  \tab prediction band
#' }
#' \emph{Example data:}
#' \tabular{ll}{
#'   \code{\link{otoliths_had}} \tab otoliths (haddock)\cr
#'   \code{\link{otoliths_skj}} \tab otoliths (skipjack)\cr
#'   \code{\link{tags_skj}}     \tab tags (skipjack)
#' }
#'
#' @note
#' The parameter estimation method follows the statistical approach of Maunder
#' et al. (2018), which stems from Aires-da-Silva et al. (2015), Eveson et al.
#' (2004), and Laslett et al. (2002), collectively called the
#' Laslett-Eveson-Polacheck approach.
#'
#' @author Arni Magnusson and Mark Maunder.
#'
#' @references
#' Aires-da-Silva, A.M., Maunder, M.N., Schaefer, K.M., and Fuller, D.W. (2015).
#' Improved growth estimates from integrated analysis of direct aging and
#' tag-recapture data: An illustration with bigeye tuna (\emph{Thunnus obesus})
#' of the eastern Pacific Ocean with implications for management.
#' \emph{Fisheries Research}, \bold{163}, 119--126.
#' \doi{10.1016/j.fishres.2014.04.001}.
#'
#' Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W., Aires-da-Silva,
#' A.M., Minte-Vera, C.V., and Campana, S.E. (2018).
#' The growth cessation model: a growth model for species showing a near
#' cessation in growth with application to bigeye tuna (\emph{Thunnus obesus}).
#' \emph{Marine Biology}, \bold{165}, 76.
#' \doi{10.1007/s00227-018-3336-9}.
#'
#' Eveson, J.P., Laslett, G.M., and Polacheck, T. (2004).
#' An integrated model for growth incorporating tag-recapture, length-frequency,
#' and direct aging data.
#' \emph{Canadian Journal of Fisheries and Aquatic Sciences}, \bold{61},
#' 292--306.
#' \doi{10.1139/f03-163}.
#'
#' Laslett, G.M., Eveson, J.P., and Polacheck, T. (2002).
#' A flexible maximum likelihood approach for fitting growth curves to
#' tag-recapture data.
#' \emph{Canadian Journal of Fisheries and Aquatic Sciences}, \bold{59},
#' 976--986.
#' \doi{10.1139/f02-069}.

"_PACKAGE"
