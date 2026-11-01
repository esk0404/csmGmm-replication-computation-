# ============================================================
# Example: Breast Cancer Replication Analysis Using csmGmm
# ============================================================
#
# This example demonstrates how to perform a replication
# analysis using two breast cancer GWAS datasets:
#
# Dataset 1: African females (GCST90551892)
# Dataset 2: European females (iCOGS, GCST010098)
#
# Both datasets are based on genome build 37.
#
# The analysis includes:
#   1. Loading GWAS summary statistics
#   2. Preparing and merging datasets
#   3. Generating initial parameter values
#   4. Fitting the csmGmm model using EM
#   5. Identifying replicated SNPs using local FDR
#   6. Generating Manhattan and Z-score plots
#
# ============================================================


# 1. Load required packages
# ------------------------------------------------------------

library(csmGmm)
library(dplyr)
library(readr)


# 2. Download GWAS summary statistics
# ------------------------------------------------------------
#
# Download the following datasets:
#
# iCOGS (GCST010098):
# https://www.ebi.ac.uk/gwas/studies/GCST010098
#
# African females (GCST90551892):
# https://www.ebi.ac.uk/gwas/studies/GCST90551892
#
# Place both files in the data/ directory.
#
# Expected file names:
#   data/GCST010098.tsv
#   data/GCST90551892.tsv.gz


# 3. Read GWAS summary statistics
# ------------------------------------------------------------

iCOGS <- read_tsv(
  "data/GCST010098.tsv",
  col_select = c(
    variant_id,
    chromosome,
    base_pair_location,
    beta,
    standard_error,
    effect_allele,
    other_allele,
    p_value
  ),
  progress = FALSE
)

Africa <- read_tsv(
  "data/GCST90551892.tsv.gz",
  col_select = c(
    variant_id,
    chromosome,
    base_pair_location,
    beta,
    standard_error,
    effect_allele,
    other_allele,
    p_value
  ),
  progress = FALSE
)


# 4. Organize datasets
# ------------------------------------------------------------
#
# Dataset ordering determines the dataset labels in the
# downstream analysis.
#
# Dataset 1: Africa
# Dataset 2: iCOGS

mydatasets <- list(Africa, iCOGS)


# 5. Prepare and merge GWAS datasets
# ------------------------------------------------------------

prepped <- prepare_csmgmm_data_new(
  datasets = mydatasets,
  dataset_names = c("Africa", "iCOGS"),
  z_cap = 8.1,
  strand_aligned = FALSE
)


# 6. Generate initial parameter values
# ------------------------------------------------------------

initParams <- generate_init_lists(2)


# 7. Fit the csmGmm model using EM
# ------------------------------------------------------------

res <- symm_fit_ind_EM(
  testStats = prepped$clean_dat,
  initMuList = initParams$initMuList,
  initPiList = initParams$initPiList,
  sameDirAlt = TRUE,
  eps = 1e-5,
  checkpoint = TRUE
)


# 8. Process local FDR results
# ------------------------------------------------------------
#
# Identify replicated SNPs using a local FDR threshold
# of 0.10.

processed <- process_lfdr_results(
  orig_data = prepped$extra_dat,
  lfdrResults = res$lfdrResults,
  fdr_threshold = 0.10
)


# 9. Save significant SNP results
# ------------------------------------------------------------

if (!dir.exists("results")) {
  dir.create("results")
}

write.csv(
  processed$sig_snps_info,
  "results/csmgmm_rep_results.csv",
  row.names = FALSE
)


# 10. Generate visualizations
# ------------------------------------------------------------

figs <- create_plots(
  plotData = processed$all_snps_info
)

# Manhattan plot
print(figs$manPlot)

# Z-score plot
print(figs$zPlot)
