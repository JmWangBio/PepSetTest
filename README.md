
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

## If statmod is not installed, run the following code:
# install.packages("statmod")

library(devtools)
install_github("JmWangBio/PepSetTest")
```

Note: please make sure that devtools is installed.

Next, I will provide one example to demonstrate how to get started with this package. 
Examples of advanced usage can be found in the vignettes.

## Example 

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

Now that the data are ready, we can run the peptide set test
workflow as follows:

``` r
# Run the workflow
result <- CompPepSetTestWorkflow(dat, contrasts.par = contrasts.par,
                                 group = group,
                                 pep_mapping_tbl = pep_mapping_tbl,
                                 stat = "t",
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
result <- CompPepSetTestWorkflow(dat, contrasts.par = contrasts.par,
                                 group = group,
                                 pep_mapping_tbl = pep_mapping_tbl,
                                 stat = "t",
                                 correlated = TRUE,
                                 equal.correlation = TRUE,
                                 pepC.estim = "mad",
                                 logged = FALSE)
```

Note that without good reason, we should almost always let PepSetTest estimate 
inter-peptide correlations on its own by setting correlated = TRUE. 
Check the documentation to learn more about each parameter.

## Data Analysis

Scripts used for generating the simulated data and analysing the data
from the spike-in experiment and breast cancer study as shown in the
manuscript are stored in the “inst” folder.
