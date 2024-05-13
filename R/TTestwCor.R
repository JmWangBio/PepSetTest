#' Two-sample t-test accounting for inter-peptide correlation
#'
#' Test whether peptides belonging to the same protein are differentially expressed relative to
#' the rest of the peptidome, accounting for inter-peptide correlation.
#' This function is adapted from \code{cameraPR()}
#' in LIMMA R package (Wu and Smyth (2012), \emph{Nucleic Acids Research}).
#'
#' @param statistic a numeric vector of peptide-wise t-statistics.
#' @param index an index vector or a list of index vectors. \code{statistic[index]} selects corresponding rows
#' in the protein of interest, i.e., test set(s).
#' @param inter.pep.cor a numeric vector of inter-peptide correlation coefficients within
#' the protein of interest, i.e., test set(s).
#' @param pepC.estim estimator of the variance of peptide-wise t-statistics not
#' belonging to the protein of interest, i.e., test set. Options include "sd" and "mad".
#' "sd" represents sample standard deviation. "mad" represents sample median
#' absolute deviation.
#'
#' @return \code{TTestwCor} returns a dataframe in which each row represents a protein of
#' interest, i.e., test set. Columns include
#' \item{NPeps}{number of peptides}
#' \item{Correlation}{inter-peptide correlation coefficient}
#' \item{Direction}{direction of change}
#' \item{PValue}{raw p-value}
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Hochberg method}
#' @export
#'
#' @author Junmin Wang
#' @references Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation. \emph{Nucleic Acids Research} 40, e133
#'
#' @examples
#'
#' y <- matrix(rnorm(1000 * 6), 1000, 6)
#' design <- cbind(Intercept = 1, Group = c(0, 0, 0, 1, 1, 1))
#'
#' # First set of 20 genes are genuinely differentially expressed
#' index1 <- 1:20
#' y[index1, 4:6] <- y[index1, 4:6]+1
#'
#' fit <- limma::eBayes(limma::lmFit(y, design))
#' TTestwCor(fit$t[, 2], index = index1,
#' inter.pep.cor = 0,
#' pepC.estim = "sd")
#'
TTestwCor <- function(statistic,
                      index,
                      inter.pep.cor,
                      pepC.estim = c("sd", "mad")) {
  storage.mode(statistic) <- "numeric"
  G <- length(statistic)
  ID <- names(statistic)
  if(G<3) stop("Too few genes in dataset: need at least 3")

  #	Check index
  if(!is.list(index)) index <- list(set1=index)
  nsets <- length(index)

  # Check inter-peptide correlation
  if(length(inter.pep.cor) > 1L) {
    if(length(inter.pep.cor) != nsets)
      stop("Length of inter.pep.cor doesn't match number of proteins")
  }

  # Compute degrees of freedom
  df.pepSetTest <- G-2L

  StatInSetLst <- list()
  StatOutSetLst <- list()
  NPeps <- rep_len(0, nsets)

  for (i in 1:nsets) {
    iset <- index[[i]]
    if(is.character(iset)) iset <- which(ID %in% iset)
    StatInSet <- statistic[iset]
    StatOutSet <- statistic[-iset]
    StatInSetLst[[i]] <- StatInSet
    StatOutSetLst[[i]] <- StatOutSet
    NPeps[i] <- length(StatInSet)
  }

  Down <- Up <- VIFs <- rep_len(0, nsets)

  # Run t-tests sequentially
  for (i in 1:nsets) {
    StatInSet <- StatInSetLst[[i]]
    StatOutSet <- StatOutSetLst[[i]]
    m <- length(StatInSet)
    vif <- 1+(m-1)*inter.pep.cor[i]
    m2 <- G-m
    meanStatInSet <- mean(StatInSet)
    varStatInSet <- stats::var(StatInSet)
    if (is.na(varStatInSet)) {
      varStatInSet <- 0
    }
    if (pepC.estim == "sd") {
      meanStatOutSet <- mean(StatOutSet)
      varStatOutSet <- stats::var(StatOutSet)
    } else if (pepC.estim == "mad") {
      meanStatOutSet <- stats::median(StatOutSet)
      varStatOutSet <- stats::mad(StatOutSet)^2
    }
    delta <- meanStatInSet - meanStatOutSet
    varStatPooled <- ((m - 1) * varStatInSet +
                        (m2 - 1) * varStatOutSet) /
      (m + m2 - 2)
    two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
    Down[i] <- stats::pt(two.sample.t, df=df.pepSetTest)
    Up[i] <- stats::pt(two.sample.t, df=df.pepSetTest, lower.tail=FALSE)
  }
  TwoSided <- 2*pmin(Down,Up)

  #	Assemble into data.frame
  D <- (Down < Up)
  Direction <- rep_len("Up", nsets)
  Direction[D] <- "Down"
  tab <- data.frame(NPeps=NPeps,
                    Correlation=inter.pep.cor,
                    Direction=Direction,
                    PValue=TwoSided,
                    stringsAsFactors=FALSE)
  rownames(tab) <- names(index)

  #	Add adjusted p-value
  if(nsets>1L) tab$adj.P.Val <- stats::p.adjust(tab$PValue,
                                                method="BH")

  tab
}
