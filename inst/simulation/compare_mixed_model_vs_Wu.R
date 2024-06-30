
## Author: Junmin Wang
## Date: March 17th, 2024
## This script compares the mixed model approach and Wu et al.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tidyr, MASS, PepSetTest, and ggplot2 R packages installed.
## Please make sure to provide the correct path to "helper.R" (line 17).
## This script will produce a simulated datasets with rho set to 0.05: proteins with equal peptide distributions (0% active).
## This dataset has 1000 simulations. Beware that running 1000 simulations requires lots of computing power. Using a high performance computing environment is strongly recommended.
## To run fewer simulations, change "1:1000" in the lapply() function.

## load libraries
library(dplyr)
library(tidyr)
library(MASS)
library(PepSetTest)
library(ggplot2)

## load helper function
source("path/to/helper.R")

## simulation
nSamples_lst <- c(3, 15, 50)
nTestPeps_lst <- c(3, 10, 30)
interPepCor <- 0.05
corr.all.df <- data.frame()

for (n in seq_along(nSamples_lst)) {
  for (m in seq_along(nTestPeps_lst)) {
    print(sprintf("nSamples: %s; nTestPeps: %s; interPepCor: %s", 
                  nSamples_lst[n],
                  nTestPeps_lst[m],
                  interPepCor))
    ## Running 1000 simulations takes a long time. Using a high performance computing environment is strongly recommended.
    corr.comb.df <- do.call('rbind', lapply(1:1000, function(k) {
      corr.df <- main_sim_for_corr_estim(nTestPeps = nTestPeps_lst[m],
                                        nTotalPeps = 9000,
                                        inter.pep.cor = interPepCor,
                                        nSamples = nSamples_lst[n],
                                        npep.trend = FALSE)
      corr.df$Iteration <- k
      return(corr.df)
    }))
    corr.comb.df$nSample <- nSamples_lst[n]
    corr.comb.df$interPepCor <- interPepCor
    corr.all.df <- rbind(corr.all.df, corr.comb.df)
  }
}

## make histogram
corr.all.long <- corr.all.df %>%
  dplyr::select(Iteration, nTestPep, nSample, pepSetTestEqCorrSD, pepSetTestUneqCorr) %>%
  pivot_longer(!c(Iteration, nTestPep, nSample), names_to = "method", values_to = "val") %>%
  mutate(Condition = factor(paste0("N = ", nSample * 2, "\n", nTestPep, "-Peptide Proteins")))

p <- ggplot(corr.all.long, aes(x = val, color = factor(ifelse(grepl("Uneq", method), 1, 0)))) +
  geom_density() +
  facet_wrap(~Condition) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.spacing = unit(2, "lines")) +
  scale_color_manual(breaks = c(0, 1), values = c("mediumpurple", "red"),
                     labels = c("Peptide Set Test (Mixed Model)", "Peptide Set Test (Wu et al)")) +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             color = "gray", alpha = 0.8)
