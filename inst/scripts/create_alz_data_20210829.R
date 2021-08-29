# Save ROSMAP Single-Platform Pathway Analysis Results
# Gabriel Odom
# 2021-08-29

library(tidyverse)



######  Single-Platform Results  ##############################################
snp_df <- read_csv(
	file = "inst/extData/alz_data/snp_all_results_bacon_correction_20210809.csv"
) %>%
	select(pathway = TERMS, snpPval = pVal, snpFDR = FDR)

dnaM_df <- read_csv(
	"inst/extData/alz_data/dnaM_c2cp_trimmed_missMethyl_pValues_20210803.csv"
) %>%
	select(pathway, dnamPval = P.DE, dnamFDR = FDR)

rnaSeq_df <- read_csv(
	"inst/extData/alz_data/rnaSeq_c2cp_trimmed_fgsea_pValues_20210810.csv"
) %>%
	select(pathway, size, rnaseqPval = pval, rnaseqFDR = padj)



######  Join and Save  ########################################################
alzheimersMultiOmics_df <-
	snp_df %>%
	full_join(dnaM_df, by = "pathway") %>%
	full_join(rnaSeq_df, by = "pathway") %>%
	select(pathway, size, everything())

usethis::use_data(alzheimersMultiOmics_df)
