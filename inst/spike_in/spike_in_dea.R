
## Author: Junmin Wang
## Date: March 17th, 2024
## Acknowledgement: a significant portion of this script is adapted from https://rformassspectrometry.github.io/docs/sec-quant.html
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, QFeatures, msdata, ggplot2, and PepSetTest R packages installed.

## load libraries
library(dplyr)
library(tibble)
library(QFeatures)
library(ggplot2)
library(PepSetTest)

## define functions
plot_volcano <- function(dat, title, yscale, ylimits) {
   ggplot(data = dat, aes(x = logFC, y = -log10(adj.P.Val),
                          color = TP)) +
      geom_point(size = 0.8) +
      labs(x = expression("Log"[2]*"(Condition B / Condition A)"),
           y = expression("-Log"[10]*"(Adj. P-Value)"),
           color = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 10),
            legend.direction = "horizontal",
            legend.text = element_text(size = 10, margin = margin(r = 30, unit = "pt"))) +
      scale_color_manual(breaks = c("Human", "Yeast"),
                         values = c("red", "gray")) +
      scale_y_continuous(trans = yscale,
                         limits = ylimits,
                         labels = function(x) {
                            format(x, scientific = FALSE,
                                   drop0trailing = TRUE)
                         }) +
      scale_x_continuous(limits = c(-8.5, 8.5)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      ggtitle(title)
}

extract_legend <- function(my_ggp) {
   step1 <- ggplot_gtable(ggplot_build(my_ggp))
   step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
   step3 <- step1$grobs[[step2]]
   return(step3)
}

## read data
f <- msdata::quant(pattern = "cptac", full.names = TRUE)
i <- grep("Intensity\\.", names(read.delim(f)))
cptac_se <- readSummarizedExperiment(f, ecol = i, 
                                     fnames = "Sequence", 
                                     sep = "\t")

## clean up sample names and annotate the experiment
colnames(cptac_se) <- sub("I.+\\.", "", colnames(cptac_se))
cptac_se$condition <- sub("_[7-9]", "", colnames(cptac_se))
cptac_se$id <- sub("^.+_", "", colnames(cptac_se))

## annotate features (peptides)
keep_var <- c("Sequence", "Proteins", "Leading.razor.protein", "PEP",
              "Score", "Reverse", "Potential.contaminant")
rowData(cptac_se) <- rowData(cptac_se)[, keep_var]

## convert zeros to NAs
cptac_se <- zeroIsNA(cptac_se)

## remove rows that have 2 or more NAs out of 6
cptac_se <- filterNA(cptac_se, pNA = 1/6)

## create the QFeatures data
cptac <- QFeatures(list(peptides = cptac_se))
colData(cptac) <- colData(cptac_se)

## filter reverse hits, contaminants, and low-confidence identifications
cptac <-
   cptac %>%
   filterFeatures(~ Reverse != "+") %>%
   filterFeatures(~ Potential.contaminant != "+") %>%
   filterFeatures(~ PEP < 0.05)

## log transform peptide abundance
cptac <- logTransform(cptac, i = "peptides",
                      name = "log_peptides")

## normalize peptide abundance
cptac <- normalize(cptac, i = "log_peptides",
                   name = "lognorm_peptides",
                   method = "center.median")

## retrieve log-norm peptides
AllDat <- assay(cptac[["lognorm_peptides"]])
group <- c(rep("A", 3), rep("B", 3))
contrasts.par <- "B-A"
pep_mapping_tbl <- rowData(cptac[["lognorm_peptides"]]) %>%
   as.data.frame(.) %>%
   tibble::rownames_to_column(var = "peptide") %>%
   dplyr::select(peptide, Leading.razor.protein) %>%
   `colnames<-`(c("peptide", "protein"))

## aggregation by summation
sumLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                    contrasts.par = contrasts.par,
                                    group = group,
                                    pep_mapping_tbl = pep_mapping_tbl,
                                    method = "sum", logged = TRUE)

## aggregation by robust regression
robRegLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                       contrasts.par = contrasts.par,
                                       group = group,
                                       pep_mapping_tbl = pep_mapping_tbl,
                                       method = "robreg", logged = TRUE)

## pep set test (mixed model, mad)
pepSetTest_EqCorr_mad_output <- CompPepSetTestWorkflow(dat = AllDat, 
                                                       contrasts.par = contrasts.par,
                                                       group = group,
                                                       pep_mapping_tbl = pep_mapping_tbl,
                                                       stat = "t", 
                                                       correlated = TRUE,
                                                       equal.correlation = TRUE,
                                                       pepC.estim = "mad",
                                                       logged = TRUE)

## add organism labels
sumLimma_output <- sumLimma_output %>% 
   dplyr::mutate(TP = ifelse(grepl("ups", protein), "Human", "Yeast"),
                 TP = factor(TP, levels = c("Human", "Yeast")))

robRegLimma_output <- robRegLimma_output %>% 
   dplyr::mutate(TP = ifelse(grepl("ups", protein), "Human", "Yeast"),
                 TP = factor(TP, levels = c("Human", "Yeast")))

pepSetTest_EqCorr_mad_output <- pepSetTest_EqCorr_mad_output %>% 
   dplyr::mutate(TP = ifelse(grepl("ups", protein), "Human", "Yeast"),
                 TP = factor(TP, levels = c("Human", "Yeast")))

## make plots
p1 <- plot_volcano(sumLimma_output,
                   title = "Protein Limma (Summed\nPeptide Intensity)",
                   yscale = "identity",
                   ylimits = c(0, 3))
p4 <- extract_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- plot_volcano(robRegLimma_output,
                   title = "Protein Limma (Mean Peptide\nIntensity via Robust Regression)",
                   yscale = "identity",
                   ylimits = c(0, 3)) +
   theme(legend.position = "none")

pepSetTest_EqCorr_mad_output <- pepSetTest_EqCorr_mad_output %>%
   mutate(adj.P.Val = ifelse(adj.P.Val > 0.9, 0.9, adj.P.Val))
p3 <- plot_volcano(pepSetTest_EqCorr_mad_output,
                   title = "Peptide Set Test\n(Mixed Model, MAD)",
                   yscale = "log10",
                   ylimits = c(-log10(0.9), 10^2)) +
   theme(legend.position = "none")

p <- gridExtra::grid.arrange(p1, p2, p3, p4, heights = c(4, 1),
                             layout_matrix = rbind(c(1, 2, 3),
                                                   c(4, 4, 4)))
