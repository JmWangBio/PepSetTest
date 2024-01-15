
## Author: Junmin Wang
## Date: January 15th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, MASS, and PepSetTest R packages installed.
## Also, please make sure to provide the correct path to "helper.R".
## This script can produce four simulated datasets with rho set to 0: a mixture of proteins (5% active), proteins with equal peptide distributions (5% active), a mixture of proteins (0% active), and proteins with equal peptide distributions (0% active).
## To run 1000 simulations, change "1:1" to "1:1000" in the lapply() functions.
## Beware that running 1000 simulations requires lots of computing power. Using a high performance computing environment is strongly recommended.

## load libraries
library(dplyr)
library(MASS)
library(PepSetTest)

## load helper function
source("path/to/helper.R")

## pval.all.df is a data frame containing the simulation results. Here is what each column stands for:
# protein: Protein ID
# nTestPep: number of peptides within this protein. It could be 3, 10, or 30.
# mu: group mean difference
# nSample: number of samples within each group (i.e., N/2). It could be 3, 15, or 50.
# interPepCor: inter-peptide correlation coefficient used to generate the simulation data.
# Iteration: index of the simulation
# sumLimma: p-value obtained from summation aggregation coupled with LIMMA
# robRegLimma: p-value obtained from robust regression aggregation coupled with LIMMA
# pepSetTestMAD: p-value obtained from peptide set test (no correlation, MAD)
# pepSetTestSD: p-value obtained from peptide set test (no correlation, SD)
# scPepSetTest: p-value obtained from self-contained peptide set test

####################################################
#### simulate a mixture of proteins (5% active) ####
####################################################
nSamples_lst <- c(3, 15, 50) 
interPepCor <- 0
pval.all.df <- data.frame()

for (n in 1:length(nSamples_lst)) {
  print(sprintf("nSamples: %s; interPepCor: %s", 
                nSamples_lst[n],
                interPepCor))
  ## change "1:1" to "1:1000" to run 1000 simulations
  ## Running 1000 simulations takes a long time. Using a high performance computing environment is strongly recommended.
  pval.comb.df <- do.call('rbind', lapply(1:1, function(k) {
    pval.df <- main_sim_active_uncorr(GroupDiff = 0.5,
                                      nTestPeps = c(3, 10, 30),
                                      nTotalPeps = c(4200, 3600, 1200),
                                      inter.pep.cor = interPepCor,
                                      nSamples = nSamples_lst[n],
                                      percDEG = 0.05,
                                      npep.trend = TRUE)
    pval.df$Iteration <- k
    return(pval.df)
  }))
  pval.comb.df$nSample <- nSamples_lst[n]
  pval.comb.df$interPepCor <- interPepCor
  pval.all.df <- rbind(pval.all.df, pval.comb.df)
}

## save result
saveRDS(pval.all.df,
        file = "sim_active_mixed_uncorr_pval_data.rds")

########################################################################
#### simulate proteins with equal peptide distributions (5% active) ####
########################################################################
nSamples_lst <- c(3, 15, 50)
nTestPeps_lst <- c(3, 10, 30)
interPepCor <- 0
pval.all.df <- data.frame()

for (n in 1:length(nSamples_lst)) {
  for (m in 1:length(nTestPeps_lst)) {
    print(sprintf("nSamples: %s; nTestPeps: %s; interPepCor: %s", 
                  nSamples_lst[n],
                  nTestPeps_lst[m],
                  interPepCor))
    ## change "1:1" to "1:1000" to run 1000 simulations
    ## Running 1000 simulations takes a long time. Using a high performance computing environment is strongly recommended.
    pval.comb.df <- do.call('rbind', lapply(1:1, function(k) {
      pval.df <- main_sim_active_uncorr(GroupDiff = 0.5,
                                        nTestPeps = nTestPeps_lst[m],
                                        nTotalPeps = 9000,
                                        inter.pep.cor = interPepCor,
                                        nSamples = nSamples_lst[n],
                                        percDEG = 0.05,
                                        npep.trend = FALSE)
      pval.df$Iteration <- k
      return(pval.df)
    }))
    pval.comb.df$nSample <- nSamples_lst[n]
    pval.comb.df$interPepCor <- interPepCor
    pval.all.df <- rbind(pval.all.df, pval.comb.df)
  }
}

## save result
saveRDS(pval.all.df,
        file = "sim_active_constnp_uncorr_pval_data.rds")

####################################################
#### simulate a mixture of proteins (NO active) ####
####################################################
nSamples_lst <- c(3, 15, 50) 
interPepCor <- 0
pval.all.df <- data.frame()

for (n in 1:length(nSamples_lst)) {
  print(sprintf("nSamples: %s; interPepCor: %s", 
                nSamples_lst[n],
                interPepCor))
  ## change "1:1" to "1:1000" to run 1000 simulations
  ## Running 1000 simulations takes a long time. Using a high performance computing environment is strongly recommended.
  pval.comb.df <- do.call('rbind', lapply(1:1, function(k) {
    pval.df <- main_sim_inactive_uncorr(nTestPeps = c(3, 10, 30),
                                        nTotalPeps = c(4200, 3600, 1200),
                                        inter.pep.cor = interPepCor,
                                        nSamples = nSamples_lst[n],
                                        npep.trend = TRUE)
    pval.df$Iteration <- k
    return(pval.df)
  }))
  pval.comb.df$nSample <- nSamples_lst[n]
  pval.comb.df$interPepCor <- interPepCor
  pval.all.df <- rbind(pval.all.df, pval.comb.df)
}

## save result
saveRDS(pval.all.df,
        file = "sim_inactive_mixed_uncorr_pval_data.rds")

########################################################################
#### simulate proteins with equal peptide distributions (NO active) ####
########################################################################
nSamples_lst <- c(3, 15, 50)
nTestPeps_lst <- c(3, 10, 30)
interPepCor <- 0
pval.all.df <- data.frame()

for (n in 1:length(nSamples_lst)) {
  for (m in 1:length(nTestPeps_lst)) {
    print(sprintf("nSamples: %s; nTestPeps: %s; interPepCor: %s", 
                  nSamples_lst[n],
                  nTestPeps_lst[m],
                  interPepCor))
    ## change "1:1" to "1:1000" to run 1000 simulations
    ## Running 1000 simulations takes a long time. Using a high performance computing environment is strongly recommended.
    pval.comb.df <- do.call('rbind', lapply(1:1, function(k) {
      pval.df <- main_sim_inactive_uncorr(nTestPeps = nTestPeps_lst[m],
                                          nTotalPeps = 9000,
                                          inter.pep.cor = interPepCor,
                                          nSamples = nSamples_lst[n],
                                          npep.trend = FALSE)
      pval.df$Iteration <- k
      return(pval.df)
    }))
    pval.comb.df$nSample <- nSamples_lst[n]
    pval.comb.df$interPepCor <- interPepCor
    pval.all.df <- rbind(pval.all.df, pval.comb.df)
  }
}

## save result
saveRDS(pval.all.df,
        file = "sim_inactive_constnp_uncorr_pval_data.rds")

