#' Calculate the MiniMax \emph{p}-Values
#'
#' @description Given a vector of pathway-level MiniMax statistics and a vector
#'    of Beta Distribution parameters, calculate the MiniMax statistics'
#'    \emph{p}-values.
#'
#' @param MiniMax_num A numeric vector of MiniMax statistics
#' @param betaParams_ls A list of the parameters for the Beta Distribution.
#'    These values should be returned by the \code{\link{MiniMax_estBetaParams}}
#'    function.
#'
#' @return A vector of the MiniMax significance levels (\emph{p}-values)
#'    corresponding to the supplied MiniMax statistics.
#'
#' @importFrom stats pbeta
#'
#' @export
#'
#' @examples
#'   # Find the best-fitting paramters based on the MiniMax statistic values
#'   #   under the null distribution
#'   mmBetaParams_ls <- MiniMax_estBetaParams(
#'     MiniMaxNull_num = nullMiniMaxResults_df$MiniMax,
#'     nPlatforms = 3L,
#'     method = "MoM"
#'   )
#'
#'   # Calculate the MiniMax Statistics for each gene set
#'   mmVals_num <- MiniMax_calculateStatistic(
#'     res_df = multiOmicsHighSignalResults_df[, -(1:2)]
#'   )
#'
#'   # Find the p-values corresponding to these statistics
#'   MiniMax_calculatePVal(
#'     MiniMax_num = mmVals_num,
#'     betaParams_ls = mmBetaParams_ls
#'  )
#'
MiniMax_calculatePVal <- function(MiniMax_num, betaParams_ls){

  pbeta(
    q = MiniMax_num,
    shape1 = betaParams_ls[["alpha"]],
    shape2 = betaParams_ls[["beta"]]
  )

}

