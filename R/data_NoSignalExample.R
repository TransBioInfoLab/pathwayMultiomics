#' Example No-Signal Multi-Omics Pathway Results
#'
#' @description An example of the input data to the MiniMax statistic parameter
#'    estimation functions with no biological (synthetic) signal (the assumed
#'    results of finding pathway significance levels for each data platform
#'    under the null hypothesis)
#'
#' @details This data represents some results from a simulation study we
#'    conducted to show the utility of the MiniMax technique. The results in
#'    this file are from a no-signal scenario (no treatment group, no treated
#'    gene sets / pathways, random binary response).
#'
#'    For the MiniMax statistic to work, the platform-specific p-values must
#'    come from a statistical test with well-controlled Type-I error rates.
#'
#' @format A tibble (a modern data frame created by the \code{tibble} package)
#'    with two columns:
#' \itemize{
#'    \item{\code{terms}}{ the name of a pathway / gene set; this column is
#'       called 'terms' for compatibility with \code{pathwayCollection} objects
#'       from package \code{pathwayPCA}.}
#'    \item{\code{pVal_CNV}}{ for copy number variation data, these are the
#'       p-values for each pathway}
#'    \item{\code{pVal_RNAseq}}{ for RNAseq expression data, these are the
#'       p-values for each pathway}
#'    \item{\code{pVal_Prot}}{ for protein expression data, these are the
#'       p-values for each pathway}
#'    \item{\code{MiniMax}}{ the MiniMax statistic value from 3 platforms}
#' }
#'
"nullMiniMaxResults_df"
