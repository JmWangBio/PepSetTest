
## Author: Junmin Wang
## Date: June 20th, 2023
## Acknowledgement: a significant portion of this script is adapted from https://rformassspectrometry.github.io/docs/sec-quant.html
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, QFeatures, and PepSetTest R packages installed.

## load libraries
library(dplyr)
library(tibble)
library(QFeatures)
library(PepSetTest)

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
pepSetTest_EqCorr_mad_output <- PepSetTestWorkflow(dat = AllDat, 
                                                      contrasts.par = contrasts.par,
                                                      group = group,
                                                      pep_mapping_tbl = pep_mapping_tbl,
                                                      stats = "t", 
                                                      correlated = TRUE,
                                                      equal.correlation = TRUE,
                                                      pepC.estim = "mad",
                                                      logged = TRUE)

## add organism labels
sumLimma_output <- sumLimma_output %>% 
   dplyr::mutate(TP = grepl("ups", protein)) %>%
   dplyr::mutate(TP = ifelse(TP, "Human", "Yeast")) %>%
   dplyr::mutate(TP = factor(TP, levels = c("Human", "Yeast")))

robRegLimma_output <- robRegLimma_output %>% 
   dplyr::mutate(TP = grepl("ups", protein)) %>%
   dplyr::mutate(TP = ifelse(TP, "Human", "Yeast")) %>%
   dplyr::mutate(TP = factor(TP, levels = c("Human", "Yeast")))

pepSetTest_EqCorr_mad_output <- pepSetTest_EqCorr_mad_output %>% 
   dplyr::mutate(TP = grepl("ups", protein)) %>%
   dplyr::mutate(TP = ifelse(TP, "Human", "Yeast")) %>%
   dplyr::mutate(TP = factor(TP, levels = c("Human", "Yeast")))
