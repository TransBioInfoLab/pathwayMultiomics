#' Calculate the MiniMax Statistics
#'
#' @description Given a data frame of pathway-level \emph{p}-values, calculate
#'    the MiniMax statistic (but not the significance level of this statistic).
#'
#' @param res_df A data frame of \emph{p}-values. The rows correspond to gene
#'    sets / pathways and the columns correspond to a data platform for the
#'    disease of interest.
#' @param orderStat How many platforms should show a biological signal for a
#'    pathway / gene set to have multi-omic "enrichment"? Defaults to 2. See
#'    "Details" for more information.
#'
#' @details The MiniMax statistic is defined as the minimum of all pairwise
#'    maxima of pathway \emph{p}-values. This operation is arithmetically
#'    equivalent to sorting the \emph{p}-values and taking the second smallest.
#'    In our experience, setting this "order statistic" cutoff to 2 is
#'    appropriate for =< 5 data platforms. Biologically, this is equivalent to
#'    saying "if this pathway is dysregulated in at least two data types for
#'    this disease / condition, it is worthy of additional consideration". In
#'    situations where more than 5 data platforms are available for the disease
#'    of interest, we recommend increasing the \code{orderStat} value to 3.
#'
#' @return A vector of the MiniMax statistic values.
#'
#' @export
#'
#' @examples
#'   MiniMax_calculateStatistic(multiOmicsHighSignalResults_df[, -(1:2)])
#'
MiniMax_calculateStatistic <- function(res_df, orderStat = 2L){

  res_mat <- as.matrix(res_df)
  res_num <- numeric(length = nrow(res_mat))

  # Calculate the MiniMax
  for(r in seq_along(res_num)){
    res_num[r] <- sort(res_mat[r, ])[orderStat]
  }

  # Return
  res_num

}

