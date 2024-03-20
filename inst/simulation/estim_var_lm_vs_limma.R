
## Author: Junmin Wang
## Date: March 17th, 2024
## This script compares the sample variance and posterior variance estimated by Limma-trend.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tidyr, MASS, PepSetTest, and ggplot2 R packages installed.

## load libraries
library(dplyr)
library(tidyr)
library(MASS)
library(PepSetTest)
library(ggplot2)

## set random state
set.seed(1)

## assign values to parameters
nTestPeps <- c(3, 10, 30)
nTotalPeps <- c(4200, 3600, 1200)
inter.pep.cor <- 0.05
nSamples <- 15

## specify metadata
pepNames <- paste0("Peptide", 1:sum(nTotalPeps))
sampleNames <- paste0("Sample", 1:(2*nSamples))
group <- c(rep("A", nSamples), rep("B", nSamples))
contrasts.par <- "B-A"
pep_mapping_tbl <- data.frame(
  peptide = pepNames,
  protein = c(rep(paste0("Protein", 
                         1:sum(floor(nTotalPeps / nTestPeps))), 
                  times = rep(nTestPeps, times = nTotalPeps / nTestPeps)))
)

## specify distribution parameters
norm.mu <- 0
norm.sigma <- 1

## simulate peptide-wise mean and variance
pepwise.means <- rnorm(sum(nTotalPeps), mean = norm.mu,
                       sd = norm.sigma)
pepwise.vars <- rep(1, sum(nTotalPeps))

## simulate data
AllDat <- do.call('rbind', 
                  lapply(1:length(nTotalPeps), function(i) {
                    do.call('rbind', lapply(1:floor(nTotalPeps[i] / nTestPeps[i]), function(j) {
                      t(MASS::mvrnorm(n = nSamples * 2, 
                                      mu = pepwise.means[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)], 
                                      Sigma = (diag(nTestPeps[i]) * 
                                                 (1-inter.pep.cor) + inter.pep.cor) * 
                                        (sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)]) %*% 
                                           t(sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)])))))
                    }))
                  }))
rownames(AllDat) <- pepNames
colnames(AllDat) <- sampleNames

## aggregate using robust regression
prot_lst <- AggPeps(dat = AllDat,
                    pep_mapping_tbl = pep_mapping_tbl,
                    method = "robreg",
                    logged = TRUE)
prot_dat <- prot_lst$int
NPeptide <- prot_lst$NPeptide

## limma-trend
eBayes.fit <- FitContrasts(prot_dat,
                           contrasts.par,
                           group,
                           logged = TRUE,
                           NPeptide = NPeptide)

## plot histogram of s^2 and s^2.post
s2.df <- data.frame(protein = names(eBayes.fit$s2.post),
                    nTestPep = rep(paste0(c(3, 10, 30), "-Peptide Proteins"),
                                   times = c(1400, 360, 40)),
                    s2_ord = eBayes.fit$sigma^2,
                    s2_post = eBayes.fit$s2.post) %>%
  pivot_longer(!c(protein, nTestPep),
               names_to = "type",
               values_to = "val") %>%
  mutate(nTestPep = factor(nTestPep, levels = unique(nTestPep)))

p <- s2.df %>%
  ggplot(., aes(x = val, color = type)) +
  geom_density() +
  facet_wrap(~nTestPep) +
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
  scale_x_continuous(limits = c(0, 1.25))
