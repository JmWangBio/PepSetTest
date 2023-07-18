
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PepSetTest

<!-- badges: start -->
<!-- badges: end -->

Peptide Set Test (PepSetTest) is a peptide-centric strategy to infer
differentially expressed proteins in LC-MS/MS proteomics data. This test
detects coordinated changes in the expression of peptides originating
from the same protein and compares these changes against the rest of the
peptidome. Compared to traditional aggregation-based approaches, the
peptide set test demonstrates improved statistical power, yet
controlling the Type I error rate correctly in most cases. This test can
be valuable for discovering novel biomarkers and prioritizing drug
targets, especially when the direct application of statistical analysis
to protein data fails to provide substantial insights.

## Installation

You can install PepSetTest like so:

``` r
## If limma is not installed, run the following code:
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("limma")

library(devtools)
install_github("JmWangBio/PepSetTest")
```

Note: please make sure that devtools is installed.

In the rest of this tutorial, I will provide a few examples to
demonstrate how to use this package.

## Example 1

In this example, I will show how a workflow based on the peptide set
test can be directly applied to the peptide-level abundance data. First,
let’s simulate some peptide abundance data. Suppose that in a
hypothetical proteomics study of cancer, 500 peptides in total are
detected. The diseased group (i.e., D) and healthy group (i.e. H) has
three replicates each. Thirty out of 500 peptides are more abundant in
the diseased group than in the healthy group by an order of 2 on the
log2 scale.

``` r
# Load library
library(PepSetTest)

# Generate random peptide data
dat <- matrix(rnorm(3000), ncol = 6)
dat[1:30, 4:6] <- dat[1:30, 4:6] + 2
dat <- 2^dat
colnames(dat) <- paste0("Sample", 1:6)
rownames(dat) <- paste0("Peptide", 1:500)
```

Next, let’s generate the group labels and contrasts. Note that the
contrast should always be expressed as “X-Y”, where X and Y are the
labels of the two groups to be compared.

``` r
# Generate group labels and contrasts
group <- c(rep("H", 3), rep("D", 3))
contrasts.par <- "D-H"
```

Afterwards, let’s simulate a peptide-protein mapping table. Suppose the
500 detected peptides are mapped to 100 proteins. Each protein has 5
peptides.

``` r
# Generate a mapping table
pep_mapping_tbl <- data.frame(peptide = paste0("Peptide", 1:500),
                              protein = paste0("Protein", rep(1:100, each = 5)))
```

Now that all the data are ready, we can run the peptide set test
workflow as follows:

``` r
# Run the workflow
result <- PepSetTestWorkflow(dat, contrasts.par = contrasts.par,
                             group = group,
                             pep_mapping_tbl = pep_mapping_tbl,
                             stats = "t",
                             correlated = FALSE,
                             equal.correlation = FALSE,
                             pepC.estim = "mad",
                             logged = FALSE)
```

The above workflow assumes that the peptides are not correlated and
adopts the sample median absolute deviation (MAD) as the estimator for
the standard deviation of peptides not belonging to the protein of
interest. The table returned gives you the fold change, p-value, and
adjusted p-value of each protein.

The parameters (e.g., “correlated”, “equal.correlation”, “pepC.estim”)
can be adjusted to your preference. For example, if you assume that each
protein has equal inter-peptide correlation coefficients, you may run
the workflow as follows:

``` r
# Run the workflow
result <- PepSetTestWorkflow(dat, contrasts.par = contrasts.par,
                             group = group,
                             pep_mapping_tbl = pep_mapping_tbl,
                             stats = "t",
                             correlated = TRUE,
                             equal.correlation = TRUE,
                             pepC.estim = "mad",
                             logged = FALSE)
```

Check the documentation to learn more about each parameter.

## Example 2

In observational designs, a linear model accounting for co-variates
(e.g., sex, age) may be used to eliminate confounding. Next, I will
demonstrate how the peptide set test can be combined with a
peptide-level linear model to infer differentially expressed proteins.
First, let’s simulate some peptide abundance data. Suppose that in a
hypothetical proteomics study of cancer, 500 peptides in total are
detected. The diseased group (i.e., D) and healthy group (i.e. H) have
six samples each: three samples collected from donors of either sex.
Thirty peptides are more abundant in the diseased group than in the
healthy group by an order of 2 on the log2 scale. The same thirty
peptides are also more abundant in the male group than in the female
group by an order of 1 on the log2 scale.

``` r
# Load library
library(PepSetTest)
library(tidyverse)

# Generate random peptide data
dat <- matrix(rnorm(6000), ncol = 12)
dat[1:30, 7:12] <- dat[1:30, 7:12] + 2
dat[1:30, c(1:3, 7:9)] <- dat[1:30, c(1:3, 7:9)] + 1
dat <- 2^dat
colnames(dat) <- paste0("Sample", 1:12)
rownames(dat) <- paste0("Peptide", 1:500)

# Convert data from wide format to long format
dat.long <- dat %>%
  as.data.frame(.) %>%
  rownames_to_column(var = "peptide") %>%
  pivot_longer(!peptide, names_to = "sample",
               values_to = "abundance")
```

Next, let’s generate the group labels and contrasts.

``` r
# Generate group labels and contrasts
group <- c(rep("H", 6), rep("D", 6))
sex <- rep(rep(c("M", "F"), each = 3), 2)
contrasts.par <- "D-H"
sample_metadata_mapping <- data.frame(sample = paste0("Sample", 1:12),
                                      group = group, 
                                      sex = sex)
dat.long <- dat.long %>%
  inner_join(sample_metadata_mapping, by = "sample")
```

Afterwards, let’s simulate a peptide-protein mapping table.

``` r
# Generate a mapping table
pep_mapping_tbl <- data.frame(peptide = paste0("Peptide", 1:500),
                              protein = paste0("Protein", rep(1:100, each = 5)))
```

Then let’s fit the peptide-level linear model to the data.

``` r
# Fit linear model for every peptide
contrasts_res <- as.data.frame(
  t(sapply(unique(pep_mapping_tbl$peptide), 
           function(x) {
             fit <- lm(abundance ~ group + sex,
                       data = dat.long %>% filter(peptide == x))
             logfc <- -as.numeric(coef(summary(fit))[, "Estimate"]["groupH"])
             t_stat <- -as.numeric(coef(summary(fit))[, "t value"]["groupH"])
             p_val <- as.numeric(coef(summary(fit))[, "Pr(>|t|)"]["groupH"])
             return(c(x, logfc, t_stat, p_val))
           }))
) %>% 
  remove_rownames() %>%
  `colnames<-`(c("feature", "logFC", "t", "P.Value")) %>%
  mutate_at(c("logFC", "t", "P.Value"), as.numeric) %>%
  mutate(adj.P.Val = p.adjust(P.Value, method = "BY"))
```

Finally, let’s run the peptide set test using the t-statistics obtained
from the peptide-level linear model.

``` r
# Run the peptide set test
result <- CompPepSetTest(contrasts_res,
                         pep_mapping_tbl = pep_mapping_tbl,
                         stats = "t",
                         cor_coef = 0,
                         pepC.estim = "mad")
```

The peptide-level data were simulated assuming no inter-peptide
correlation. If peptides are instead assumed to be inter-correlated,
then the correlation coefficients need to be estimated prior to running
the peptide set test using the code as follows:

``` r
# Estimate inter-peptide correlation coefficients
group <- factor(group)
sex <- factor(sex)
design.m <- stats::model.matrix(~ 0 + group + sex)
pep_cors <- EstimInterPepCor(dat, contrasts.par, 
                             design.m,
                             pep_mapping_tbl, 
                             equal.correlation = TRUE,
                             logged = FALSE)

# Run the peptide set test
result <- CompPepSetTest(contrasts_res,
                         pep_mapping_tbl = pep_mapping_tbl,
                         stats = "t",
                         cor_coef = pep_cors,
                         pepC.estim = "mad")
```

## Data Analysis

Scripts used for generating the simulated data and analysing the data
from the spike-in experiment and breast cancer study as shown in the
manuscript are stored in the “inst” folder.
