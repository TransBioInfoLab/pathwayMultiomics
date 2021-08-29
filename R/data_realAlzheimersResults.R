#' Real Multi-Omics Data Results for AD
#'
#' @description Real multi-omics data results for Alzheimer's disease using SNP,
#'    DNA methylation, and gene expression data from the ROSMAP study
#'
#' @details This data gives pathway-level statistical significance results from
#'    GSEA under Bacon correction (for SNP analysis), \code{missMethyl} (for DNA
#'    methylation analysis), and \code{fgsea} (for gene expression analysis) for
#'    640 samples from the ROSMAP study using the Broad Institute's C2 CP
#'    collection (2833 pathways with the number of genes less than 200 or more
#'    than 4).
#'
#' @format A tibble (a modern data frame created by the \code{tibble} package)
#'    with eight columns and 2833 rows:
#' \itemize{
#'    \item{\code{pathway}}{ the name of the pathway / gene set in C2 CP}
#'    \item{\code{size}}{ the number of genes in the pathway / gene set}
#'    \item{\code{snpPval}}{ the *p*-values for each pathway after SNP GWAS
#'      analysis}
#'    \item{\code{snpFDR}}{ the false discovery rate for each SNP *p*-value}
#'    \item{\code{dnamPval}}{ the *p*-values for each pathway after DNA
#'      methylation analysis with \code{missMethyl::}}
#'    \item{\code{dnamFDR}}{ the false discovery rate for each DNAm *p*-value}
#'    \item{\code{rnaseqPval}}{ the *p*-values for each pathway after RNAseq
#'      analysis with \code{fgsea::}}
#'    \item{\code{rnaseqFDR}}{ the false discovery rate for each RNAseq *p*-value}
#' }
#'
"alzheimersMultiOmics_df"
