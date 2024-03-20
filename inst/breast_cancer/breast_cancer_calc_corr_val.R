
## Author: Junmin Wang
## Date: March 17th, 2024
## This script calculates the proteome-wide inter-peptide correlations.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, stringr, ggplot2, and latex2exp R packages installed.
## In addition, you need the "peptides.txt" file (MaxQuant search output) from De Marchi et al (2016), Mol. Oncol.
## This file can be downloaded from ProteomeXchange with dataset identifier PXD000485, which refers to the NKI-AVL + RUMC dataset.
## After downloading this file, change "path/to/peptides.txt" to where "peptides.txt"is located on your computer (line 51).

## load libraries
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(latex2exp)

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

make_histogram <- function(dat, col_name) {
  ggplot(dat, aes_string(x = col_name)) +
    geom_histogram(binwidth = 0.1,
                   fill = "lightblue",
                   color = "black") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = TeX("$\\hat{\\bar{\\rho}}_y$"),
         y = "# of Proteins") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.spacing = unit(2, "lines")) +
    geom_vline(xintercept = mean(dat[, col_name, drop = TRUE]), color = "red",
               linetype = "dashed") +
    annotate("text", x = 0.1, y = 200, 
             label = sprintf("Mean: %s", round(mean(dat[, col_name, drop = TRUE]), 2)))
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

## add protein annotation
df.nor.filtered.annot <- as.data.frame(df.nor.filtered) %>%
  rownames_to_column(var = "peptide") %>%
  inner_join(PepNames.filtered, by = "peptide")

## calculate correlation
all_corr_df <- data.frame()
all_prots <- unique(df.nor.filtered.annot$protein)

for (prot in all_prots) {
  tmpdat <- df.nor.filtered.annot %>%
    filter(protein == prot)
  if (nrow(tmpdat) == 1) {
    next
  }
  for (i in (1:(nrow(tmpdat) - 1))) {
    tmpdat_OR_1 <- tmpdat[i, 2:25]
    tmpdat_PD_1 <- tmpdat[i, 26:57]
    tmppep_1 <- tmpdat[i, 1]
    for (j in ((i+1):nrow(tmpdat))) {
      tmpdat_OR_2 <- tmpdat[j, 2:25]
      tmpdat_PD_2 <- tmpdat[j, 26:57]
      tmppep_2 <- tmpdat[j, 1]
      tmp_OR_corr <- cor(as.numeric(tmpdat_OR_1),
                         as.numeric(tmpdat_OR_2),
                         method = "pearson",
                         use = "pairwise.complete.obs")
      tmp_PD_corr <- cor(as.numeric(tmpdat_PD_1),
                         as.numeric(tmpdat_PD_2),
                         method = "pearson",
                         use = "pairwise.complete.obs")
      all_corr_df <- rbind(all_corr_df, 
                           data.frame(Pep1 = tmppep_1,
                                      Pep2 = tmppep_2,
                                      ORCorr = tmp_OR_corr,
                                      PDCorr = tmp_PD_corr,
                                      Prot = prot))
    }
  }
}

## calculate mean correlation for each protein
mean_corr_by_prot_df <- all_corr_df %>%
  group_by(Prot) %>%
  summarise(meanPDCorr = mean(PDCorr, na.rm = TRUE),
            meanORCorr = mean(ORCorr, na.rm = TRUE))

## make plots
p1 <- make_histogram(mean_corr_by_prot_df, col_name = "meanORCorr")
p2 <- make_histogram(mean_corr_by_prot_df, col_name = "meanPDCorr")
