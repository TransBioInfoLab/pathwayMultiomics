#' Estimate the Parameters of the MiniMax Statistic's Beta Distribution
#'
#' @description Given a vector of MiniMax statisic values under the null
#'    hypothesis, estimate the parameters of the Beta Distribution which best
#'    fits these values.
#'
#' @param MiniMaxNull_num A numeric vector of MiniMax statistics under the null
#' @param nPlatforms An integer stating how many data platforms are in the
#'    original data.
#' @param orderStat How many platforms should show a biological signal for a
#'    pathway / gene set to have multi-omic "enrichment"? Defaults to 2. See
#'    "Details" for more information.
#' @param method Which estimation method will be used to find the parameters of
#'    the Beta Distribution? Options are \code{"parametric"} (no estimation from
#'    the data), \code{"MLE"} (Maximum Likelihood Estimates), or \code{"MoM"}
#'    (Method of Moments estimates). See "Details" for more information.
#'
#' @details
#'   \strong{Concerning Parameter Estimation Methods:} We currently support 3
#'   options to estimate the parameters of the Beta Distribution. The
#'   "parametric" option does not use the data. Instead, it assumes that the
#'   MiniMax statistics will have a Beta \eqn{(k, n + 1 - k)} distribution,
#'   where \eqn{k} is the value of \code{orderStat} and \eqn{n} has the value
#'   \code{nPlatforms}. See \url{https://en.wikipedia.org/wiki/Order_statistic}.
#'
#'   The next two estimation options make use of the \code{MiniMaxNull_num}
#'   vector, which should be calculated by finding the same significance levels
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
#' @return A list of 3 components: "alpha" and "beta" hold the parameter
#'    estimates of the Beta Distribution, and "method" returns a character
#'    string denoting which estimation method was used.
#'
#' @importFrom stats var
#' @importFrom Rfast beta.mle
#'
#' @export
#'
#'
#' @examples
#'  miniMax_num <- nullMiniMaxResults_df$MiniMax
#'
#'  MiniMax_estBetaParams(miniMax_num, nPlatforms = 3L)
#'  MiniMax_estBetaParams(miniMax_num, nPlatforms = 3L, method = "MoM")
#'  MiniMax_estBetaParams(miniMax_num, nPlatforms = 3L, method = "MLE")
#'

MiniMax_estBetaParams <- function(MiniMaxNull_num,
                                  nPlatforms, orderStat = 2L,
                                  method = c("parametric", "MLE", "MoM")){
  # browser()
  method = match.arg(method)

  switch(method,
    "parametric" = {

      out_ls <- list(
        # See wikipedia.org/wiki/Order_statistic
        alpha = orderStat, beta = nPlatforms + 1 - orderStat,
        method = "Parametric"
      )

    },

    "MoM" = {

      # See wikipedia.org/wiki/Beta_distribution#Method_of_moments
      xBar <- mean(MiniMaxNull_num)
      vBar <- var(MiniMaxNull_num)
      guts <- (xBar * (1 - xBar)) / vBar - 1

      if(guts > 0) {

        out_ls <- list(
          alpha = xBar * guts, beta = (1 - xBar) * guts,
          method = "Method of Moments"
        )

      } else {

        out_ls <- list(
          alpha = NaN, beta = NaN,
          method = "Method of Moments"
        )

      }

    },

    "MLE" = {
      # browser()

      # The Rfast::beta.mle() function can't handle p-values within 1e-09 of
      #   {0, 1}, so detect these and "fudge" them away. Since small p-values
      #   are more valuable than large p-values, we will preserve the ordering
      #   for small values, but truncate the large values (after shifting).
      if(any(MiniMaxNull_num <= 1e-09)){
        MiniMaxNull_num <- MiniMaxNull_num + 1e-09
      }
      bigP_lgl <- MiniMaxNull_num >= (1 - 1e-09)
      if(any(bigP_lgl)){
        MiniMaxNull_num[bigP_lgl] <- 1 - 1e-09
      }

      out_ls <- tryCatch(
        expr = {

          beta_fit <- beta.mle(MiniMaxNull_num)
          list(
            alpha = beta_fit$param[["alpha"]],
            beta  = beta_fit$param[["beta"]],
            method = "Maximum Likelihood"
          )

        },
        error = function(cond){

          message("MLE failed to converge. Please select a different method.")
          list(
            alpha = NA_real_, beta = NA_real_,
            method = "Maximum Likelihood"
          )

        }
      )

    }

  )

  class(out_ls) <- c("MiniMaxParams", "list")
  out_ls

}
