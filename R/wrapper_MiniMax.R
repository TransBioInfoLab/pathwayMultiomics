#' Find Gene Set / Pathway Significance across Multi-Omics Data
#'
#' @description Given a data frame of pathway-level \emph{p}-values across
#'    multiple -omics platforms, use the MiniMax technique to assign statistical
#'    significance to concordant or cascading pathway-level biological effects.
#'
#' @param pValues_df A data frame of pathway / gene set \emph{p}-values under
#'    true responses (this data set should contain true biological signal). The
#'    rows correspond to gene sets / pathways, and the columns correspond to the
#'    data platforms for the disease of interest.
#' @param pValuesNull_df A data frame of pathway / gene set \emph{p}-values under
#'    the null hypothesis, most likely constructed from randomly permuting the
#'    response and re-estimating all significance levels (this data set should
#'    NOT contain any true biological signal). As with \code{pValues_df}, the
#'    rows correspond to gene sets / pathways, and the columns correspond to the
#'    data platforms for the disease of interest. NOTE: if this data set is not
#'    provided, only \code{method = "parametric"} will be available.
#' @param orderStat How many platforms should show a biological signal for a
#'    pathway / gene set to have multi-omic "enrichment"? Defaults to 2. See
#'    "Details" for more information.
#' @param method If \code{pValuesNull_df} is provided, which estimation method
#'    will be used to find the parameters of the Beta Distribution? Options are
#'    \code{"parametric"} (no estimation from the data; this should be used only
#'    in cases where no MiniMax statistics under the null hypothesis are
#'    available, such as in the case of pure meta-analysis approaches),
#'    \code{"MLE"} (Maximum Likelihood Estimates), or \code{"MoM"} (Method of
#'    Moments estimates). Using \code{"MLE"} or \code{"MoM"} requires the user
#'    to provide \code{pValuesNull_df}. See "Details" for more information.
#'
#' @details
#'   \strong{Concerning Parameter Estimation Methods:} We currently support 3
#'   options to estimate the parameters of the Beta Distribution. The
#'   "parametric" option does not use the data, and it is therefore the only
#'   option available if \code{pValuesNull_df} is not provided. Instead, it
#'   assumes that the MiniMax statistics will have a Beta \eqn{(k, n + 1 - k)}
#'   distribution, where \eqn{k} is the value of \code{orderStat} and \eqn{n}
#'   has the value \code{nPlatforms}.
#'   See \url{https://en.wikipedia.org/wiki/Order_statistic}.
#'
#'   The next two estimation options make use of the \code{pValuesNull_df} data
#'   frame, which should be calculated by finding the same significance levels
#'   of the statistical tests used on the real data (for each pathway and data
#'   platform), but by using a random permutation of the outcome of interest
#'   instead of the real values; more permutations are better. The "MLE" option
#'   uses the \code{\link[Rfast]{beta.mle}} function to find the Maximum
#'   Likelihood Estimates of \eqn{\alpha} and \eqn{\beta}. The "MoM" option uses
#'   the closed-form Method of Moments estimators of \eqn{\alpha} and
#'   \eqn{\beta} as shown in
#'   \url{https://en.wikipedia.org/wiki/Beta_distribution#Method_of_moments}.
#'
#'   \strong{Concerning Appropriate Order Statistics:} The MiniMax operation is
#'    equivalent to sorting the \emph{p}-values and taking the second smallest.
#'    In our experience, setting this "order statistic" cutoff to 2 is
#'    appropriate for =< 5 data platforms. Biologically, this is equivalent to
#'    saying "if this pathway is dysregulated in at least two data types for
#'    this disease / condition, it is worthy of additional consideration". In
#'    situations where more than 5 data platforms are available for the disease
#'    of interest, we recommend increasing the \code{orderStat} value to 3.
#'
#' @return A copy of the \code{pValues_df} data frame with two additional
#'    columns: \code{MiniMax} (the statistic values for each gene set) and
#'    \code{MiniMaxP} (the \emph{p}-values of these statistics). This data frame
#'    is sorted by ascending MiniMax \emph{p}-value.
#'
#' @export
#'
#'
#' @examples
#'  data("multiOmicsMedSignalResults_df")
#'  data("nullMiniMaxResults_df")
#'
#'  MiniMax(
#'    pValues_df = multiOmicsMedSignalResults_df,
#'    pValuesNull_df = nullMiniMaxResults_df[, -5],
#'    method = "MLE"
#'  )
#'

MiniMax <- function(pValues_df, pValuesNull_df = NULL,
                    orderStat = 2L,
                    method = c("parametric", "MLE", "MoM")){
  # browser()

  ###  Check Inputs  ###
  method = match.arg(method)

  # Extract p-value columns
  resCols_lgl <- vapply(pValues_df, is.numeric, logical(1))
  nPlatforms_int <- sum(resCols_lgl)
  if(orderStat >= nPlatforms_int) {
    stop("The number of platforms must be greater than the order statisic.")
  }


  ###  Estimate Beta Parameters  ###
  if(is.null(pValuesNull_df)){

    if(method != "parametric") {
      warning("This method requires pValuesNull_df. Using parametric estimates.")
    }
    bestParams_ls <- MiniMax_estBetaParams(
      nPlatforms = nPlatforms_int,
      orderStat = orderStat,
      method = "parametric"
    )

  } else {

    resNullCols_lgl <- vapply(pValuesNull_df, is.numeric, logical(1))
    if(nPlatforms_int != sum(resNullCols_lgl)) {
      stop("The number of platforms in pValues_df and pValuesNull_df must be equal.")
    }

    mmNull_num <- MiniMax_calculateStatistic(
      res_df = pValuesNull_df[, resNullCols_lgl, drop = FALSE],
      orderStat = orderStat
    )

    bestParams_ls <- MiniMax_estBetaParams(
      MiniMaxNull_num = mmNull_num,
      nPlatforms = nPlatforms_int,
      orderStat = orderStat,
      method = method
    )

  }


  ###  Calculate the MiniMax Statistic and p-Value  ###
  # To Do: create an error if these columns already exist
  pValues_df[["MiniMax"]] <- MiniMax_calculateStatistic(
    res_df = pValues_df[, resCols_lgl, drop = FALSE],
    orderStat = orderStat
  )
  pValues_df[["MiniMaxP"]] <- MiniMax_calculatePVal(
    MiniMax_num = pValues_df[["MiniMax"]],
    betaParams_ls = bestParams_ls
  )


  ###  Return  ###
  pValues_df[order(pValues_df[["MiniMaxP"]]), ]

}
