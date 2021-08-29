# Create Example Pathway Analysis Results Table
# Gabriel Odom
# 2021-06-11

# This script will take three sets of pathway analysis results on simulated
#   data under low (partition1_delta0.2), medium (partition2_delta0.3), and high
#   (partition4_delta0.5) synthetic signal. The design parameters are:
#   - partition: what proportion of the genes in each pathway received
#     treatment? We used 5 total partitions, so each partition represents an
#     additional 20% of the genes in the gene set (partition1 = 20%, partition2
#     = 40%, ..., partition4 = 80% of the genes in the gene set received
#     treatment).
#   - delta: what was the strength of the treatment in effect size? We centred
#     and scaled all assay / expression values to N(0,1). Therefore, the added
#     signal represents a 0.1, 0.2, ..., 0.5 standard deviation change to the
#     treated samples.
# The data sets were originally stored in miniMax_p_simulation/simAnalysis/
#   under the simulation pathwayPCA_sim_20201202/. These particular files are
#   from simulation 001.
# We will also save a "no signal" data set for examples of the MLE and MoM
#   estimation techniques. This will be a subset of one of the other example
#   data sets, but only with pathway results for untreated pathways.

library(tidyverse)
#



######  Create Low-Signal Example Data  #######################################
CNV_ls <- readRDS(
  "inst/extData/CNV_partition1_delta0.2.RDS"
)
RNAseq_ls <- readRDS(
  "inst/extData/RNAseq_partition1_delta0.2.RDS"
)
Prot_ls <- readRDS(
  "inst/extData/Prot_partition1_delta0.2.RDS"
)

multiOmicsLowSignalResults_df <-
  CNV_ls$pVals_df %>%
  arrange(terms) %>%
  select(terms, treated, rawp) %>%
  rename(pVal_CNV = rawp) %>%
  left_join(
    RNAseq_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_RNAseq = rawp)
  ) %>%
  left_join(
    Prot_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_Prot = rawp)
  )

usethis::use_data(multiOmicsLowSignalResults_df)



######  Create Medium-Signal Example Data  ####################################
CNV_ls <- readRDS(
  "inst/extData/CNV_partition2_delta0.3.RDS"
)
RNAseq_ls <- readRDS(
  "inst/extData/RNAseq_partition2_delta0.3.RDS"
)
Prot_ls <- readRDS(
  "inst/extData/Prot_partition2_delta0.3.RDS"
)

multiOmicsMedSignalResults_df <-
  CNV_ls$pVals_df %>%
  arrange(terms) %>%
  select(terms, treated, rawp) %>%
  rename(pVal_CNV = rawp) %>%
  left_join(
    RNAseq_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_RNAseq = rawp)
  ) %>%
  left_join(
    Prot_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_Prot = rawp)
  )

usethis::use_data(multiOmicsMedSignalResults_df)



######  Create High-Signal Example Data  ######################################
CNV_ls <- readRDS(
  "inst/extData/CNV_partition4_delta0.5.RDS"
)
RNAseq_ls <- readRDS(
  "inst/extData/RNAseq_partition4_delta0.5.RDS"
)
Prot_ls <- readRDS(
  "inst/extData/Prot_partition4_delta0.5.RDS"
)

multiOmicsHighSignalResults_df <-
  CNV_ls$pVals_df %>%
  arrange(terms) %>%
  select(terms, treated, rawp) %>%
  rename(pVal_CNV = rawp) %>%
  left_join(
    RNAseq_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_RNAseq = rawp)
  ) %>%
  left_join(
    Prot_ls$pVals_df %>%
      select(terms, rawp) %>%
      rename(pVal_Prot = rawp)
  )

usethis::use_data(multiOmicsHighSignalResults_df)


######  Create No-Signal Example Data  ########################################
# This data set is to show how the estimation techniques work.
data(multiOmicsLowSignalResults_df)

nullMiniMaxResults_df <-
  multiOmicsLowSignalResults_df %>%
  filter(!treated) %>%
  select(-treated) %>%
  mutate(
    MiniMax = pathwayMiniMax::MiniMax_calculateStatistic(., orderStat = 2L)
  )

usethis::use_data(nullMiniMaxResults_df, overwrite = TRUE)
