#' Estimation of inter-peptide correlation
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' estimate inter-peptide correlation coefficient for each protein via the
#' mixed model approach or approach described in Wu and Smyth (2012), \emph{Nucleic Acids Research}.
#'
#' @inheritParams FitContrasts
#' @inheritParams AggPeps
#' @param design.m design matrix
#' @param equal.correlation Boolean variable indicating whether all pairwise
#' inter-peptide correlation coefficients are assumed to be equal within a protein.
#' If true, the mixed model approach will be applied; otherwise, the approach
#' described in Wu and Smyth (2012), \emph{Nucleic Acids Research} will be applied.
#'
#' @return \code{EstimInterPepCor} returns a numeric vector of inter-peptide
#' correlation coefficients (one value for each protein).
#' @export
#'
#' @references Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation. \emph{Nucleic Acids Research} 40, e133
#'
#' @author Junmin Wang and Steven Novick
#'
#' @examples
#' # Generate random peptide data
#' dat <- 2^matrix(rnorm(540), ncol = 6)
#' colnames(dat) <- paste0("Sample", 1:6)
#' rownames(dat) <- paste0("Peptide", 1:90)
#'
#' # Generate peptide mapping table
#' pep_mapping_tbl <- data.frame(peptide = paste0("Peptide", 1:90),
#' protein = paste0("Protein", rep(1:30, each = 3)))
#'
#' # Generate groups and contrasts
#' group <- c(rep("A", 3), rep("B", 3))
#' contrasts.par <- "B-A"
#'
#' # Generate design matrix
#' group <- factor(group)
#' design.m <- stats::model.matrix(~ 0 + group)
#' 
#' EstimInterPepCor(dat, contrasts.par, design.m, pep_mapping_tbl,
#' equal.correlation = TRUE, logged = FALSE)
#'
EstimInterPepCor <- function(dat, contrasts.par, design.m,
                             pep_mapping_tbl,
                             equal.correlation = FALSE,
                             logged = FALSE) {
  dat <- dat[stats::complete.cases(dat), ]
  if (logged) {
    dat.m <- as.matrix(dat)
  } else {
    dat.m <- as.matrix(log2(dat))
  }
  ## keep peptides with abundance data only
  pep_mapping_tbl <- pep_mapping_tbl[pep_mapping_tbl$peptide %in% rownames(dat.m), ]
  ## extract all proteins
  target_prot_ids <- unique(as.character(pep_mapping_tbl$protein))
  ## obtain peptide indices for each protein
  pep_lst <- lapply(target_prot_ids, function(x) {
    pep_mapping_tbl[pep_mapping_tbl$protein == x, "peptide"]
  })
  names(pep_lst) <- target_prot_ids
  ## calculate correlation
  inter.pep.cor.lst <-
    sapply(target_prot_ids, function(x) {
      if (length(pep_lst[[x]]) >= 2) {
        if (equal.correlation) {
          inter.pep.cor <- FitLmerBySample(dat.m[pep_lst[[x]], ],
                                           design.m)
        } else {
          inter.pep.cor <- limma::interGeneCorrelation(dat.m[pep_lst[[x]], ],
                                                       design.m)$correlation
        }
      } else {
        inter.pep.cor <- NA
      }
      return(inter.pep.cor)
    })
  names(inter.pep.cor.lst) <- target_prot_ids
  return(inter.pep.cor.lst)
}
