
## Author: Junmin Wang
## Date: June 20th, 2023
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, stringr, and PepSetTest R packages installed.
## In addition, you need the "peptides.txt" files (MaxQuant search output) from De Marchi et al (2016), Mol. Oncol.
## These files can be downloaded from ProteomeXchange with dataset identifiers PXD000484 and PXD000485.
## The former refers to the EMC dataset, and the latter refers to the NKI-AVL + RUMC dataset.
## After downloading these files, change "path/to/peptides.txt" to where "peptides.txt"is located on your computer.

## load libraries
library(dplyr)
library(tibble)
library(stringr)
library(PepSetTest)

## define functions
Normalize <- function (dat.df) {
  dat.df.nor <-	t(t(dat.df) * ( 
    mean(apply(dat.df, 2,
               function(x) median(x, na.rm=T))
    ) / apply(
      dat.df, 2, 
      function(x) median(x, na.rm=T)))
  )
  return(dat.df.nor)
}

## load data
pep_df <- read.delim("path/to/peptides.txt",  ## change path to where "peptides.txt" is located
                     sep = "\t")

## discard contaminants, reverse sequences, and peptides with PEP score < 0.05
pep_df_filtered <- pep_df %>%
  dplyr::filter(PEP  < 0.05,
         !grepl("CON_|REV_", .$Leading.razor.protein))

## extract peptide sequence + protein ID mapping table
PepNames <- pep_df_filtered %>%
  dplyr::select("Sequence", "Leading.razor.protein") %>%
  dplyr::distinct() %>%
  `colnames<-`(c("peptide", "protein"))
  
## extract protein ID + gene name mapping table
ProNames <- pep_df_filtered %>%
  dplyr::select("Leading.razor.protein", "Gene.names") %>%
  dplyr::mutate(Gene.names = sapply(Gene.names, function(x) {
    str_split(x, ';')[[1]][1]
  })) %>%
  dplyr::filter(!duplicated(Leading.razor.protein)) %>%
  `colnames<-`(c("protein", "gene"))

## extract peptide abundance 
dat.m <- pep_df_filtered %>%
  dplyr::select("Sequence", grep("Intensity.", colnames(.), 
                                 value = TRUE)) %>%
  tibble::column_to_rownames(var = "Sequence")

## extract group labels
Group <- as.character(sapply(colnames(dat.m), function(x) {
  stringr::str_split(x, "\\.")[[1]][2]
}))
contrasts.par <- "PD-OR"

## convert 0s to NAs
dat.m[dat.m == 0] <- NA

## remove peptides >= 30% missing values in either group
yvar.px <- rownames(dat.m)
over_miss_peps <- union(rownames(dat.m[rowSums(is.na(dat.m[, grep("OR", colnames(dat.m), value=TRUE)])) / ncol(dat.m[, grep("OR", colnames(dat.m), value=TRUE)]) > 0.3, ]),
                        rownames(dat.m[rowSums(is.na(dat.m[, grep("PD", colnames(dat.m), value=TRUE)])) / ncol(dat.m[, grep("PD", colnames(dat.m), value=TRUE)]) > 0.3, ]))
yvar.px.2 <- setdiff(yvar.px, over_miss_peps)
dat.m.filtered <- dat.m[yvar.px.2, ]
PepNames.filtered <- PepNames[PepNames$peptide %in% yvar.px.2, ]

## normalization
df.nor.filtered <- Normalize(as.matrix(dat.m.filtered))

## aggregation by summation
sumLimma_output <- AggLimmaWorkflow(dat = df.nor.filtered,
                                    contrasts.par = contrasts.par,
                                    group = Group,
                                    pep_mapping_tbl = PepNames.filtered,
                                    method = "sum", logged = FALSE)

sumLimma_output <- sumLimma_output %>%
  dplyr::inner_join(ProNames, by = "protein")

## aggregation by robust regression
robRegLimma_output <- AggLimmaWorkflow(dat = df.nor.filtered,
                                       contrasts.par = contrasts.par,
                                       group = Group,
                                       pep_mapping_tbl = PepNames.filtered,
                                       method = "robreg", logged = FALSE)

robRegLimma_output <- robRegLimma_output %>%
  dplyr::inner_join(ProNames, by = "protein")

## peptide set test
pepSetTest_eq_pep_corr_mad_output <- PepSetTestWorkflow(dat = df.nor.filtered,
                                                      contrasts.par = contrasts.par,
                                                      group = Group,
                                                      pep_mapping_tbl = PepNames.filtered,
                                                      stats = "t", 
                                                      correlated = TRUE,
                                                      equal.correlation = TRUE,
                                                      pepC.estim = "mad",
                                                      logged = FALSE)

pepSetTest_eq_pep_corr_mad_output <- pepSetTest_eq_pep_corr_mad_output %>%
  dplyr::inner_join(ProNames, by = "protein")
