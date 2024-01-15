
## Author: Junmin Wang
## Date: January 15th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr R package installed.
## This script can calculate the power and Type I error rate of the statistical methods applied to the simulation data.
## Increasing the number of simulations in "sim_corr_val.R" and "sim_uncorr_val.R" will provide a more accurate estimation of power and Type I error rates.

## load libraries
library(dplyr)

## calculate the percentage of simulations with p-value less than the threshold for each protein
## change "methods" to the methods you want the power or Type I error rate to be calculated for
calc_sig_perc <- function(dat, pval.threshold = 0.05, 
                          methods = "sumLimma") {
  dat %>% 
    dplyr::group_by(nSample, nTestPep, protein) %>%
    dplyr::summarise_at(methods,
                 list(perc = function(x) {
                   sum(x < pval.threshold) / length(x)
                 }))
}

## example
## change file path to where the rds file is located.
pval.all.df <- readRDS("path/to/sim_active_mixed_corr_pval_data.rds")

## separate the proteome into active and inactive proteins
pval.all.df.zero <- pval.all.df %>% 
  dplyr::filter(mu == 0)

pval.all.df.posneg <- pval.all.df %>%
  dplyr::filter(mu != 0)

## Type I error rate
type1.err.rate.wide <- calc_sig_perc(pval.all.df.zero, methods = c("sumLimma", 
                                                                   "robRegLimma",
                                                                   "pepSetTestEqCorrMAD",
                                                                   "pepSetTestEqCorrSD",
                                                                   "pepSetTestUneqCorr",
                                                                   "scPepSetTest"))

## power
power.wide <- calc_sig_perc(pval.all.df.posneg, methods = c("sumLimma", 
                                                            "robRegLimma",
                                                            "pepSetTestEqCorrMAD",
                                                            "pepSetTestEqCorrSD",
                                                            "pepSetTestUneqCorr",
                                                            "scPepSetTest"))
