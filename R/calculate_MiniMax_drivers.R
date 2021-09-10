#' Mark which Platforms are Driving the MiniMax Statistics
#'
#' @description Given a data frame of pathway-level \emph{p}-values, mark which
#'    of the multi-omics platforms are driving the MiniMax statistic.
#'
#' @param res_df A data frame of \emph{p}-values. The rows correspond to gene
#'    sets / pathways and the columns correspond to a data platform for the
#'    disease of interest.
#' @param orderStat How many platforms should show a biological signal for a
#'    pathway / gene set to have multi-omic "enrichment"? Defaults to 2. See
#'    "Details" for more information.
#' @param drivers_char What labels should be given to the driving platforms?
#'    Defaults to the column names of \code{res_df}. If you supply custom labels,
#'    make sure to match them to the column order of \code{res_df}.
#' @param sortLabels Should the driver labels be sorted alphabetically before
#'    concatenation? Defaults to \code{TRUE}; that is, a multi-omics result
#'    driven first by protein expression then by DNA methylation will have the
#'    same label as a result driven first by DNA methylation then by protein
#'    expression. If you would like the magnitude of the \emph{p}-value to set
#'    the label order, then use \code{sortLabels = FALSE}.
#' @param separator What character string should be used to separate the names
#'    of the driving platforms? Defaults to \code{" and "}; for example, if the
#'    platform driver labels are \code{"protein"} and \code{"cnv"}, and if
#'    \code{sortLabels = TRUE}, then the label of drivers would be
#'    \code{"cnv and protein"}.
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
#'    NOTE: this result does not depend on the pathway significance level at all.
#'    This result will simply show you which platforms had the smallest
#'    \emph{p}-values for a particular pathway, even if the MiniMax statistic is
#'    not statistically significant for that pathway. Therefore, we recommend
#'    that this function be used only for interpretation of results post-hoc.
#'
#' @return A vector of the names of the platforms driving the MiniMax statistic
#'    values.
#'
#' @export
#'
#' @examples
#'   MiniMax_calculateDrivers(
#'     multiOmicsHighSignalResults_df[, -(1:2)],
#'     drivers_char = c("cnv", "rnaSeq", "protein")
#'  )
#'
MiniMax_calculateDrivers <- function(res_df,
                                     orderStat = 2L,
                                     drivers_char = colnames(res_df),
                                     sortLabels = TRUE,
                                     separator = " and "){


  res_mat  <- as.matrix(res_df)
  if(orderStat >= ncol(res_mat)) {
    stop("The number of platforms must be greater than the order statisic.")
  }

  nDrivers_int <- seq_len(orderStat)
  res_char <- character(length = nrow(res_mat))

  # Calculate the MiniMax
  for(r in seq_along(res_char)){

    whichPlatforms_idx <- order(res_mat[r, ])[nDrivers_int]
    labels_char <- drivers_char[whichPlatforms_idx]
    if (sortLabels) {
      labels_char <- sort(labels_char)
    }

    res_char[r] <- paste(labels_char, collapse = separator)

  }

  # Return
  res_char

}
