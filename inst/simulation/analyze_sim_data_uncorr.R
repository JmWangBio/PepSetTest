
## Author: Junmin Wang
## Date: March 17th, 2024
## This script takes simulated p-values as input and calculates the power, Type I error rate, FDR, and TPR of the statistical methods applied to the simulation data.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tidyr, and ggplot2 R packages installed.
## Make sure to change the paths of input to where you have saved them (lines 140 - 143).

## load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

## define functions
# calculate the percentage of simulations with p-value less than the threshold for each protein
# change "methods" to the methods you want the power or Type I error rate to be calculated for
calc_sig_perc <- function(dat, pval.threshold = 0.05, 
                          methods = "sumLimma") {
  dat %>% 
    dplyr::group_by(nSample, nTestPep, protein) %>%
    dplyr::summarise_at(methods,
                 list(perc = function(x) {
                   sum(x < pval.threshold) / length(x)
                 }))
}

# adjust p-values to correct for multiple testing
adj_pval <- function(dat, methods) {
  dat %>%
    dplyr::group_by(nSample, Iteration) %>%
    dplyr::mutate_at(methods, 
                     list(bh = function(x) {
                       p.adjust(x, method = "BH")
                     },
                     by = function(x) {
                       p.adjust(x, method = "BY")
                     })) %>%
    dplyr::ungroup()
}

# calculate FDR (average of FDP)
calc_FDR <- function(dat, adjp.threshold, methods) {
  FDP <- dat %>%
    dplyr::group_by(nSample, Iteration, nTestPep) %>%
    dplyr::summarise_at(
      methods,
      list(~ifelse(sum(. < 0.05, na.rm = TRUE) == 0, 0,
                   sum(. < 0.05 & mu == 0, na.rm = TRUE) / sum(. < 0.05, na.rm = TRUE)))
    ) %>%
    ungroup()
  
  FDR <- FDP %>%
    dplyr::group_by(nSample, nTestPep) %>%
    dplyr::summarise_at(methods,
                        list(fdr = function(x) {
                          mean(x)
                        }))
  
  return(FDR)
}

# calculate TPR (average of TPP)
calc_TPR <- function(dat, adjp.threshold, methods) {
  TPP <- dat %>%
    dplyr::group_by(nSample, Iteration, nTestPep) %>%
    dplyr::summarise_at(
      methods,
      list(~ifelse(sum(mu != 0, na.rm = TRUE) == 0, 0,
                   sum(. < 0.05 & mu != 0, na.rm = TRUE) / sum(mu != 0, na.rm = TRUE)))
    ) %>%
    dplyr::ungroup()
  
  TPR <- TPP %>%
    dplyr::group_by(nSample, nTestPep) %>%
    dplyr::summarise_at(methods,
                        list(tpr = function(x) {
                          mean(x)
                        }))
  
  return(TPR)
}

# make box plot
make_boxplot <- function(dat, ylabel, methods) {
  dat.long <- dat %>%
    pivot_longer(!c(nSample, nTestPep, protein),
                 names_to = "method",
                 values_to = "val") %>%
    filter(grepl(paste0(methods, collapse = "|"), method))
  
  suffix <- "_perc"
  ggplot(dat.long %>%
           filter(grepl(suffix, method)) %>%
           mutate(method = factor(method, levels = paste0(methods, suffix)),
                  nTestPep = factor(nTestPep)),
         aes(x = nTestPep, y = val, fill = method)) +
    geom_boxplot(outlier.shape = 0.1, lwd = 0.1) +
    facet_wrap(~nSample, labeller = as_labeller(appender)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 7, margin = margin(r = 5, unit = "pt")),
          legend.title = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.spacing = unit(2, "lines")) +
    labs(x = "Number of Peptides",
         y = ylabel)
}

appender <- function(string) {
  paste0("N = ", 2 * as.numeric(string))
}

# make heatmap
make_heatmap <- function(dat) {
  dat %>%
    ggplot(., aes(x = method, y = nTestPep, fill = val)) +
    geom_tile(show.legend = FALSE) +
    geom_text(aes(label = round(val, 2))) +
    facet_wrap(~nSample, ncol = 1) +
    scale_fill_gradient(low = "white", 
                        high = "lightblue") + 
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(color = "black",
                                     angle = 90),
          axis.text.y = element_text(color = "black")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "# of Peptides")
}

## change file paths to where the simulated p-value files are located.
inactive.constnp.pval <- readRDS("sim_inactive_constnp_uncorr_pval_data.rds")
inactive.mixed.pval <- readRDS("sim_inactive_mixed_uncorr_pval_data.rds")
active.constnp.pval <- readRDS("sim_active_constnp_uncorr_pval_data.rds")
active.mixed.pval <- readRDS("sim_active_mixed_uncorr_pval_data.rds")

## separate the 5% active proteome into active and inactive proteins
active.constnp.pval.zero <- active.constnp.pval %>% 
  dplyr::filter(mu == 0)
active.constnp.pval.posneg <- active.constnp.pval %>% 
  dplyr::filter(mu != 0)
active.mixed.pval.zero <- active.mixed.pval %>% 
  dplyr::filter(mu == 0)
active.mixed.pval.posneg <- active.mixed.pval %>% 
  dplyr::filter(mu != 0)

## Type I error rate
inactive.constnp.pval.t1err <- calc_sig_perc(inactive.constnp.pval, methods = c("sumLimma", 
                                                                                "robRegLimma",
                                                                                "pepSetTestSD",
                                                                                "scPepSetTest"))

inactive.mixed.pval.t1err <- calc_sig_perc(inactive.mixed.pval, methods = c("sumLimma", 
                                                                            "robRegLimma",
                                                                            "pepSetTestSD",
                                                                            "scPepSetTest"))

active.constnp.pval.t1err <- calc_sig_perc(active.constnp.pval.zero, methods = c("sumLimma", 
                                                                                 "robRegLimma",
                                                                                 "pepSetTestMAD",
                                                                                 "pepSetTestSD",
                                                                                 "scPepSetTest"))

active.mixed.pval.t1err <- calc_sig_perc(active.mixed.pval.zero, methods = c("sumLimma", 
                                                                             "robRegLimma",
                                                                             "pepSetTestMAD",
                                                                             "pepSetTestSD",
                                                                             "scPepSetTest"))

## make box plots
p1 <- make_boxplot(inactive.constnp.pval.t1err, ylabel = "Type I Error Rate",
                   methods = c("sumLimma",
                               "robRegLimma",
                               "pepSetTestSD",
                               "scPepSetTest"))

p2 <- make_boxplot(inactive.mixed.pval.t1err, ylabel = "Type I Error Rate",
                   methods = c("sumLimma", 
                               "robRegLimma",
                               "pepSetTestSD",
                               "scPepSetTest"))

p3 <- make_boxplot(active.constnp.pval.t1err, ylabel = "Type I Error Rate",
                   methods = c("sumLimma",
                               "robRegLimma",
                               "pepSetTestMAD",
                               "pepSetTestSD",
                               "scPepSetTest"))

p4 <- make_boxplot(active.mixed.pval.t1err, ylabel = "Type I Error Rate",
                   methods = c("sumLimma", 
                               "robRegLimma",
                               "pepSetTestMAD",
                               "pepSetTestSD",
                               "scPepSetTest"))

## power
active.constnp.pval.power <- calc_sig_perc(active.constnp.pval.posneg, methods = c("sumLimma",
                                                                                   "robRegLimma",
                                                                                   "pepSetTestMAD",
                                                                                   "pepSetTestSD",
                                                                                   "scPepSetTest"))

active.mixed.pval.power <- calc_sig_perc(active.mixed.pval.posneg, methods = c("sumLimma",
                                                                               "robRegLimma",
                                                                               "pepSetTestMAD",
                                                                               "pepSetTestSD",
                                                                               "scPepSetTest"))

## make box plots
p5 <- make_boxplot(active.constnp.pval.power, ylabel = "Power",
                   methods = c("sumLimma",
                               "robRegLimma",
                               "pepSetTestMAD",
                               "pepSetTestSD",
                               "scPepSetTest"))

p6 <- make_boxplot(active.mixed.pval.power, ylabel = "Power",
                   methods = c("sumLimma",
                               "robRegLimma",
                               "pepSetTestMAD",
                               "pepSetTestSD",
                               "scPepSetTest"))

## adjust p-values
adjp.all.df <- adj_pval(active.mixed.pval, methods = c("sumLimma",
                                                 "robRegLimma",
                                                 "pepSetTestMAD",
                                                 "pepSetTestSD",
                                                 "scPepSetTest"))

## FDR
fdr.wide <- calc_FDR(adjp.all.df, adjp.threshold = 0.05,
                     methods = c("sumLimma_bh",
                                 "robRegLimma_bh",
                                 "pepSetTestMAD_by",
                                 "pepSetTestMAD_bh",
                                 "pepSetTestSD_by",
                                 "pepSetTestSD_bh",
                                 "scPepSetTest_bh"))

fdr.long <- fdr.wide %>%
  tidyr::pivot_longer(!c(nSample, nTestPep), 
                      names_to = "method",
                      values_to = "val") %>%
  mutate(nTestPep = factor(nTestPep),
         nSample = factor(paste0("N = ", 2 * nSample)),
         method = factor(method, 
                         levels = c("sumLimma_bh_fdr",
                                    "robRegLimma_bh_fdr",
                                    "pepSetTestMAD_by_fdr",
                                    "pepSetTestMAD_bh_fdr",
                                    "pepSetTestSD_by_fdr",
                                    "pepSetTestSD_bh_fdr",
                                    "scPepSetTest_bh_fdr")))

## TPR
tpr.wide <- calc_TPR(adjp.all.df, adjp.threshold = 0.05,
                     methods = c("sumLimma_bh",
                                 "robRegLimma_bh",
                                 "pepSetTestMAD_by",
                                 "pepSetTestMAD_bh",
                                 "pepSetTestSD_by",
                                 "pepSetTestSD_bh",
                                 "scPepSetTest_bh"))

tpr.long <- tpr.wide %>%
  tidyr::pivot_longer(!c(nSample, nTestPep), 
                      names_to = "method",
                      values_to = "val") %>%
  mutate(nTestPep = factor(nTestPep),
         nSample = factor(paste0("N = ", 2 * nSample)),
         method = factor(method, 
                         levels = c("sumLimma_bh_tpr",
                                    "robRegLimma_bh_tpr",
                                    "pepSetTestMAD_by_tpr",
                                    "pepSetTestMAD_bh_tpr",
                                    "pepSetTestSD_by_tpr",
                                    "pepSetTestSD_bh_tpr",
                                    "scPepSetTest_bh_tpr")))

## make heatmaps
p7 <- make_heatmap(fdr.long)
p8 <- make_heatmap(tpr.long)
