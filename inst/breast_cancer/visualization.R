
## Author: Junmin Wang
## Date: March 17th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the ggplot2, ggrepel, dplyr, tibble, tidyr, and gridExtra R packages installed.
## Before executing this script, run "breast_cancer_dea_train.R" and "breast_cancer_dea_val.R" to generate the required input.
## Make sure to change the paths of input to where you have saved them (lines 91, 92, 95, 96, 99, 100, 103, 104).

## load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(gridExtra)

## define functions
make_volcano <- function(dat) {
  ggplot(data = dat, aes(x = logFC,
                         y = -log10(adj.P.Val),
                         color = ifelse(adj.P.Val < 0.05 | gene == "COL12A1",
                                        "red",
                                        "grey"))) +
    geom_point(size = 0.5) +
    scale_x_continuous(limits = c(-2.3, 2.3)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(-log10(0.95), 10^2),
                       trans = "log10",
                       labels = function(x) {
                         format(x, scientific = FALSE,
                                drop0trailing = TRUE)
                       }) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    ylab(expression("-Log"[10]*"(Adj. P-Value)")) +
    xlab(expression("Log"[2]*"(Poor / Good)")) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(0.05),
               col = "black",
               linetype = "dashed",
               size = 0.3,
               alpha = 0.8) +
    geom_text_repel(data = dat %>%
                      filter(gene == "COL12A1"),
                    aes(logFC, -log10(adj.P.Val), label = gene),
                    size = 3.5, color = "steelblue",
                    segment.color = NA,
                    min.segment.length = unit(0, "line"),
                    nudge_y = 0.05,
                    max.overlaps = Inf)
}

make_violin <- function(dat) {
  ggplot(data = dat %>%
           mutate(group = ifelse(gene == "COL12A1", "COL12A1", "not in\nCOL12A1")),
         aes(x = group, y = logFC, fill = group)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(width = 0.1, fill = "white",
                 color = "black", alpha = 1,
                 outlier.shape = NA) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(x= "",
         y = expression("Log"[2]*"(Poor / Good)")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("mediumpurple",
                                 "lightblue",
                                 "black")) +
    geom_jitter(width = 0.12,
                alpha = 1.0,
                color = "black",
                aes(size = ifelse(group == "COL12A1", 0.2, NA))) +
    scale_size_area(max_size = 0.3) +
    ylim(c(-2.3, 2.3)) +
    annotate("text", x = c(1, 2),
             y = 1, label = c("", ""))
}


## load summation Limma results
bc_train_sumLimma_output <- readRDS("path/to/bc_train_sumLimma_res.rds")  ## change path to where this file is located
bc_val_sumLimma_output <- readRDS("path/to/bc_val_sumLimma_res.rds")  ## change path to where this file is located

## load robust reg. Limma results
bc_train_robRegLimma_output <- readRDS("path/to/bc_train_robRegLimma_res.rds")  ## change path to where this file is located
bc_val_robRegLimma_output <- readRDS("path/to/bc_val_robRegLimma_res.rds")  ## change path to where this file is located

## load peptide set test results
bc_train_pepSetTest_eq_pep_corr_mad_output <- readRDS("path/to/bc_train_pepSetTest_res.rds")  ## change path to where this file is located
bc_val_pepSetTest_eq_pep_corr_mad_output <- readRDS("path/to/bc_val_pepSetTest_res.rds")  ## change path to where this file is located

# load peptide t-test results
bc_train_pep_dea_res <- readRDS("path/to/bc_train_pep_dea_res.rds")  ## change path to where this file is located
bc_val_pep_dea_res <- readRDS("path/to/bc_val_pep_dea_res.rds")  ## change path to where this file is located


## make volcano plots
p1 <- make_volcano(bc_train_sumLimma_output) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  ggtitle("EMC")

p2 <- make_volcano(bc_val_sumLimma_output) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  ggtitle("NKI-AVL + RUMC")

p3 <- make_volcano(bc_train_robRegLimma_output) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  ggtitle("EMC")

p4 <- make_volcano(bc_val_robRegLimma_output) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  ggtitle("NKI-AVL + RUMC")

bc_train_pepSetTest_eq_pep_corr_mad_output_volcano <- bc_train_pepSetTest_eq_pep_corr_mad_output %>%
  mutate(adj.P.Val = ifelse(adj.P.Val > 0.9, 0.9, adj.P.Val))
p5 <- make_volcano(bc_train_pepSetTest_eq_pep_corr_mad_output_volcano) +
  ggtitle("EMC")

bc_val_pepSetTest_eq_pep_corr_mad_output_volcano <- bc_val_pepSetTest_eq_pep_corr_mad_output %>%
  mutate(adj.P.Val = ifelse(adj.P.Val > 0.9, 0.9, adj.P.Val))
p6 <- make_volcano(bc_val_pepSetTest_eq_pep_corr_mad_output_volcano) +
  ggtitle("NKI-AVL + RUMC")


## make bar plot
# assign categories to genes
bc_train_pepSetTest_eq_pep_corr_mad_output$Category <- 0
bc_train_pepSetTest_eq_pep_corr_mad_output[bc_train_pepSetTest_eq_pep_corr_mad_output$gene %in% c("NDUFA6", "NDUFA13", "NDUFA8"), "Category"] <- 1

bc_val_pepSetTest_eq_pep_corr_mad_output$Category <- 0
bc_val_pepSetTest_eq_pep_corr_mad_output[bc_val_pepSetTest_eq_pep_corr_mad_output$gene %in% c("FN1", "COL12A1", "FBN1", "COL6A3", "COL6A1"), "Category"] <- 2

barplot_data <- rbind(bc_train_pepSetTest_eq_pep_corr_mad_output, 
                      bc_val_pepSetTest_eq_pep_corr_mad_output)
barplot_data <- barplot_data[barplot_data$Category > 0, ]

# reorder barplot data by fold change
barplot_data <- barplot_data %>%
  group_by(Category) %>%
  arrange(Category,
          desc(logFC)) %>%
  mutate(gene = factor(gene, levels = rev(unique(gene))))

# plot
p7 <- ggplot(barplot_data, aes(x = logFC, y = gene,
                               fill = ifelse(logFC > 0, "pos", "neg"))) +
  geom_col(width = 0.7) +
  labs(x = expression("Log"[2]*"(Poor / Good)"),
       fill = "",
       color = "") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        legend.position = "none") +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c("lightblue", "khaki"))


## make heat map
heatmap_data <- barplot_data[, c("gene", "Up", "Down", "NPeps")]
heatmap_data <- heatmap_data %>%
  pivot_longer(!c(gene, NPeps), 
               names_to = "Dir",
               values_to = "Num") %>%
  mutate(Perc = ifelse(Dir == "Up", Num / NPeps,
                       - Num / NPeps))

p8 <- ggplot(heatmap_data, aes(x = Dir, y = gene, fill = Perc)) +
  geom_tile(col = "white", size = 0.5) +
  geom_text(aes(label = Num)) +
  scale_fill_gradient2(limits = c(-1, 1),
                       low = "lightblue",
                       mid = "white",
                       midpoint = 0,
                       high = "yellow") +
  labs(x = " ", y = " ") +
  theme_classic() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,
                                    color = alpha("black", 1),
                                    size = 0.5),
        plot.title = element_blank(),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), labels = c("G", "P"))

p9 <- grid.arrange(p7, p8, ncol = 2, nrow = 1,
                  widths = c(1.5, 0.4))


## make violin plots
p10 <- make_violin(bc_train_pep_dea_res) +
  ggtitle("EMC")

p11 <- make_violin(bc_val_pep_dea_res) +
  ggtitle("NKI-AVL + RUMC")
