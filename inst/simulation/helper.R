
## simulate proteins with correlation (5% active)
main_sim_active_corr <- function(GroupDiff = 0.5, 
                                 nTestPeps = c(3, 10, 30), 
                                 nTotalPeps = c(4200, 3600, 1200), 
                                 inter.pep.cor = 0.05, 
                                 nSamples = 3, 
                                 percDEG = 0.05) {
  ####################
  ## run simulation ##
  ####################
  ## specify metadata
  pepNames <- paste0("Peptide", 1:sum(nTotalPeps))
  sampleNames <- paste0("Sample", 1:(2*nSamples))
  group <- c(rep("A", nSamples), rep("B", nSamples))
  contrasts.par <- "B-A"
  pep_mapping_tbl <- data.frame(
    peptide = pepNames,
    protein = c(rep(paste0("Protein", 
                           1:sum(floor(nTotalPeps / nTestPeps))), 
                    times = rep(nTestPeps, times = nTotalPeps / nTestPeps)))
  )
  
  ## specify distribution parameters
  norm.mu <- 0
  norm.sigma <- 1
  
  ## simulate peptide-wise mean and variance
  pepwise.means <- rnorm(sum(nTotalPeps), mean = norm.mu,
                         sd = norm.sigma)
  pepwise.vars <- rep(1, sum(nTotalPeps))
  
  ## simulate data
  AllDat <- do.call('rbind', 
                    lapply(1:length(nTotalPeps), function(i) {
                      do.call('rbind', lapply(1:floor(nTotalPeps[i] / nTestPeps[i]), function(j) {
                        t(MASS::mvrnorm(n = nSamples * 2, 
                                        mu = pepwise.means[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)], 
                                        Sigma = (diag(nTestPeps[i]) * 
                                                   (1-inter.pep.cor) + inter.pep.cor) * 
                                          (sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)]) %*% 
                                             t(sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)])))))
                      }))
                    }))
  
  ## add group mean difference to x % of peptides
  nDEGs_2 <- floor(percDEG * floor(nTotalPeps / nTestPeps) / 2)
  all_DE_pep_index <- c()
  up_DE_prot_index <- c()
  down_DE_prot_index <- c()
  for (i in 1:length(nDEGs_2)) {
    up_pep_index <- (1:(nDEGs_2[i] * nTestPeps[i])) + cumsum(c(0, nTotalPeps))[i]
    down_pep_index <- ((nDEGs_2[i] * nTestPeps[i] + 1):(2 * nDEGs_2[i] * nTestPeps[i])) + cumsum(c(0, nTotalPeps))[i]
    up_prot_index <- (1:nDEGs_2[i]) + cumsum(c(0, nTotalPeps / nTestPeps))[i]
    down_prot_index <- ((nDEGs_2[i] + 1):(2 * nDEGs_2[i])) + cumsum(c(0, nTotalPeps / nTestPeps))[i]
    AllDat[up_pep_index, (nSamples+1):(2*nSamples)] <- 
      AllDat[up_pep_index, (nSamples+1):(2*nSamples)] + GroupDiff
    AllDat[down_pep_index, (nSamples+1):(2*nSamples)] <- 
      AllDat[down_pep_index, (nSamples+1):(2*nSamples)] - GroupDiff
    all_DE_pep_index <- c(all_DE_pep_index, up_pep_index, down_pep_index)
    up_DE_prot_index <- c(up_DE_prot_index, up_prot_index)
    down_DE_prot_index <- c(down_DE_prot_index, down_prot_index)
  }
  
  rownames(AllDat) <- pepNames
  colnames(AllDat) <- sampleNames
  
  # ## run protein LIMMA with summation
  sumLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                      contrasts.par = contrasts.par,
                                      group = group,
                                      pep_mapping_tbl = pep_mapping_tbl,
                                      method = "sum", 
                                      logged = TRUE)
  
  ## run protein LIMMA with robust regression
  robRegLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                         contrasts.par = contrasts.par,
                                         group = group,
                                         pep_mapping_tbl = pep_mapping_tbl,
                                         method = "robreg", 
                                         logged = TRUE)
  
  ## run peptide set test
  pepSetTest_eq_corr_mad_output <- PepSetTestWorkflow(dat = AllDat, 
                                                      contrasts.par = contrasts.par,
                                                      group = group,
                                                      pep_mapping_tbl = pep_mapping_tbl,
                                                      stats = "t", 
                                                      correlated = TRUE,
                                                      equal.correlation = TRUE,
                                                      pepC.estim = "mad",
                                                      logged = TRUE)
  
  pepSetTest_eq_corr_sd_output <- PepSetTestWorkflow(dat = AllDat, 
                                                     contrasts.par = contrasts.par,
                                                     group = group,
                                                     pep_mapping_tbl = pep_mapping_tbl,
                                                     stats = "t", 
                                                     correlated = TRUE,
                                                     equal.correlation = TRUE,
                                                     pepC.estim = "sd",
                                                     logged = TRUE)
  
  pepSetTest_uneq_corr_output <- PepSetTestWorkflow(dat = AllDat, 
                                                    contrasts.par = contrasts.par,
                                                    group = group,
                                                    pep_mapping_tbl = pep_mapping_tbl,
                                                    stats = "t", 
                                                    correlated = TRUE,
                                                    equal.correlation = FALSE,
                                                    pepC.estim = "sd",
                                                    logged = TRUE)
  
  ## retrieve P values
  proteins <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(protein)
  num.peps <- pepSetTest_eq_corr_mad_output %>% dplyr::arrange(protein) %>% dplyr::pull(NPeps)
  sumLimma.pval <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  robRegLimma.pval <- robRegLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  pepSetTest.eq.corr.mad.pval <- pepSetTest_eq_corr_mad_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  pepSetTest.eq.corr.sd.pval <- pepSetTest_eq_corr_sd_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  pepSetTest.uneq.corr.pval <- pepSetTest_uneq_corr_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  
  pval.df <- data.frame( protein = proteins,
                         nTestPep = num.peps,
                         sumLimma = sumLimma.pval,
                         robRegLimma = robRegLimma.pval,
                         pepSetTestEqCorrMAD = pepSetTest.eq.corr.mad.pval,
                         pepSetTestEqCorrSD = pepSetTest.eq.corr.sd.pval,
                         pepSetTestUneqCorr = pepSetTest.uneq.corr.pval )
  pval.df$mu <- 0
  pval.df[pval.df$protein %in% paste0("Protein", up_DE_prot_index), "mu"] <- GroupDiff
  pval.df[pval.df$protein %in% paste0("Protein", down_DE_prot_index), "mu"] <- -GroupDiff
  return(pval.df)
}


## simulate proteins without correlation (5% active)
main_sim_active_uncorr <- function(GroupDiff = 0.5, 
                                   nTestPeps = c(3, 10, 30), 
                                   nTotalPeps = c(4200, 3600, 1200), 
                                   inter.pep.cor = 0, 
                                   nSamples = 3, 
                                   percDEG = 0.05) {
  ####################
  ## run simulation ##
  ####################
  ## specify metadata
  pepNames <- paste0("Peptide", 1:sum(nTotalPeps))
  sampleNames <- paste0("Sample", 1:(2*nSamples))
  group <- c(rep("A", nSamples), rep("B", nSamples))
  contrasts.par <- "B-A"
  pep_mapping_tbl <- data.frame(
    peptide = pepNames,
    protein = c(rep(paste0("Protein", 
                           1:sum(floor(nTotalPeps / nTestPeps))), 
                    times = rep(nTestPeps, times = nTotalPeps / nTestPeps)))
  )
  
  ## specify distribution parameters
  norm.mu <- 0
  norm.sigma <- 1
  
  ## simulate peptide-wise mean and variance
  pepwise.means <- rnorm(sum(nTotalPeps), mean = norm.mu,
                         sd = norm.sigma)
  pepwise.vars <- rep(1, sum(nTotalPeps))
  
  ## simulate data
  AllDat <- do.call('rbind', 
                    lapply(1:length(nTotalPeps), function(i) {
                      do.call('rbind', lapply(1:floor(nTotalPeps[i] / nTestPeps[i]), function(j) {
                        t(MASS::mvrnorm(n = nSamples * 2, 
                                        mu = pepwise.means[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)], 
                                        Sigma = (diag(nTestPeps[i]) * 
                                                   (1-inter.pep.cor) + inter.pep.cor) * 
                                          (sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)]) %*% 
                                             t(sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)])))))
                      }))
                    }))
  
  ## add group mean difference to x % of peptides
  nDEGs_2 <- floor(percDEG * floor(nTotalPeps / nTestPeps) / 2)
  all_DE_pep_index <- c()
  up_DE_prot_index <- c()
  down_DE_prot_index <- c()
  for (i in 1:length(nDEGs_2)) {
    up_pep_index <- (1:(nDEGs_2[i] * nTestPeps[i])) + cumsum(c(0, nTotalPeps))[i]
    down_pep_index <- ((nDEGs_2[i] * nTestPeps[i] + 1):(2 * nDEGs_2[i] * nTestPeps[i])) + cumsum(c(0, nTotalPeps))[i]
    up_prot_index <- (1:nDEGs_2[i]) + cumsum(c(0, nTotalPeps / nTestPeps))[i]
    down_prot_index <- ((nDEGs_2[i] + 1):(2 * nDEGs_2[i])) + cumsum(c(0, nTotalPeps / nTestPeps))[i]
    AllDat[up_pep_index, (nSamples+1):(2*nSamples)] <- 
      AllDat[up_pep_index, (nSamples+1):(2*nSamples)] + GroupDiff
    AllDat[down_pep_index, (nSamples+1):(2*nSamples)] <- 
      AllDat[down_pep_index, (nSamples+1):(2*nSamples)] - GroupDiff
    all_DE_pep_index <- c(all_DE_pep_index, up_pep_index, down_pep_index)
    up_DE_prot_index <- c(up_DE_prot_index, up_prot_index)
    down_DE_prot_index <- c(down_DE_prot_index, down_prot_index)
  }
  
  rownames(AllDat) <- pepNames
  colnames(AllDat) <- sampleNames
  
  # ## run protein LIMMA with summation
  sumLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                      contrasts.par = contrasts.par,
                                      group = group,
                                      pep_mapping_tbl = pep_mapping_tbl,
                                      method = "sum", 
                                      logged = TRUE)
  
  ## run protein LIMMA with robust regression
  robRegLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                         contrasts.par = contrasts.par,
                                         group = group,
                                         pep_mapping_tbl = pep_mapping_tbl,
                                         method = "robreg", 
                                         logged = TRUE)
  
  ## run peptide set test
  pepSetTest_mad_output <- PepSetTestWorkflow(dat = AllDat, 
                                              contrasts.par = contrasts.par,
                                              group = group,
                                              pep_mapping_tbl = pep_mapping_tbl,
                                              stats = "t", 
                                              correlated = FALSE,
                                              equal.correlation = FALSE,
                                              pepC.estim = "mad",
                                              logged = TRUE)
  
  pepSetTest_sd_output <- PepSetTestWorkflow(dat = AllDat, 
                                             contrasts.par = contrasts.par,
                                             group = group,
                                             pep_mapping_tbl = pep_mapping_tbl,
                                             stats = "t", 
                                             correlated = FALSE,
                                             equal.correlation = FALSE,
                                             pepC.estim = "sd",
                                             logged = TRUE)
  
  ## retrieve P values
  proteins <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(protein)
  num.peps <- pepSetTest_mad_output %>% dplyr::arrange(protein) %>% dplyr::pull(NPeps)
  sumLimma.pval <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  robRegLimma.pval <- robRegLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  pepSetTest.mad.pval <- pepSetTest_mad_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  pepSetTest.sd.pval <- pepSetTest_sd_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  
  pval.df <- data.frame( protein = proteins,
                         nTestPep = num.peps,
                         sumLimma = sumLimma.pval,
                         robRegLimma = robRegLimma.pval,
                         pepSetTestMAD = pepSetTest.mad.pval,
                         pepSetTestSD = pepSetTest.sd.pval )
  pval.df$mu <- 0
  pval.df[pval.df$protein %in% paste0("Protein", up_DE_prot_index), "mu"] <- GroupDiff
  pval.df[pval.df$protein %in% paste0("Protein", down_DE_prot_index), "mu"] <- -GroupDiff
  return(pval.df)
}


## simulate proteins with correlation (NO active)
main_sim_inactive_corr <- function(nTestPeps = c(3, 10, 30), 
                                   nTotalPeps = c(4200, 3600, 1200), 
                                   inter.pep.cor = 0.05, 
                                   nSamples = 3) {
  ####################
  ## run simulation ##
  ####################
  ## specify metadata
  pepNames <- paste0("Peptide", 1:sum(nTotalPeps))
  sampleNames <- paste0("Sample", 1:(2*nSamples))
  group <- c(rep("A", nSamples), rep("B", nSamples))
  contrasts.par <- "B-A"
  pep_mapping_tbl <- data.frame(
    peptide = pepNames,
    protein = c(rep(paste0("Protein", 
                           1:sum(floor(nTotalPeps / nTestPeps))), 
                    times = rep(nTestPeps, times = nTotalPeps / nTestPeps)))
  )
  
  ## specify distribution parameters
  norm.mu <- 0
  norm.sigma <- 1
  
  ## simulate peptide-wise mean and variance
  pepwise.means <- rnorm(sum(nTotalPeps), mean = norm.mu,
                         sd = norm.sigma)
  pepwise.vars <- rep(1, sum(nTotalPeps))
  
  ## simulate data
  AllDat <- do.call('rbind', 
                    lapply(1:length(nTotalPeps), function(i) {
                      do.call('rbind', lapply(1:floor(nTotalPeps[i] / nTestPeps[i]), function(j) {
                        t(MASS::mvrnorm(n = nSamples * 2, 
                                        mu = pepwise.means[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)], 
                                        Sigma = (diag(nTestPeps[i]) * 
                                                   (1-inter.pep.cor) + inter.pep.cor) * 
                                          (sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)]) %*% 
                                             t(sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)])))))
                      }))
                    }))
  rownames(AllDat) <- pepNames
  colnames(AllDat) <- sampleNames
  
  # ## run protein LIMMA with summation
  sumLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                      contrasts.par = contrasts.par,
                                      group = group,
                                      pep_mapping_tbl = pep_mapping_tbl,
                                      method = "sum", 
                                      logged = TRUE)
  
  ## run protein LIMMA with robust regression
  robRegLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                         contrasts.par = contrasts.par,
                                         group = group,
                                         pep_mapping_tbl = pep_mapping_tbl,
                                         method = "robreg", 
                                         logged = TRUE)
  
  ## run peptide set test
  pepSetTest_eq_corr_sd_output <- PepSetTestWorkflow(dat = AllDat, 
                                                     contrasts.par = contrasts.par,
                                                     group = group,
                                                     pep_mapping_tbl = pep_mapping_tbl,
                                                     stats = "t", 
                                                     correlated = TRUE,
                                                     equal.correlation = TRUE,
                                                     pepC.estim = "sd",
                                                     logged = TRUE)
  
  pepSetTest_uneq_corr_output <- PepSetTestWorkflow(dat = AllDat, 
                                                    contrasts.par = contrasts.par,
                                                    group = group,
                                                    pep_mapping_tbl = pep_mapping_tbl,
                                                    stats = "t", 
                                                    correlated = TRUE,
                                                    equal.correlation = FALSE,
                                                    pepC.estim = "sd",
                                                    logged = TRUE)
  
  ## retrieve P values
  proteins <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(protein)
  num.peps <- pepSetTest_eq_corr_sd_output %>% dplyr::arrange(protein) %>% dplyr::pull(NPeps)
  sumLimma.pval <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  robRegLimma.pval <- robRegLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  pepSetTest.eq.corr.sd.pval <- pepSetTest_eq_corr_sd_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  pepSetTest.uneq.corr.pval <- pepSetTest_uneq_corr_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  
  pval.df <- data.frame( protein = proteins,
                         nTestPep = num.peps,
                         sumLimma = sumLimma.pval,
                         robRegLimma = robRegLimma.pval,
                         pepSetTestEqCorrSD = pepSetTest.eq.corr.sd.pval,
                         pepSetTestUneqCorr = pepSetTest.uneq.corr.pval )
  pval.df$mu <- 0
  return(pval.df)
}


## simulate proteins without correlation (NO active)
main_sim_inactive_uncorr <- function(nTestPeps = c(3, 10, 30), 
                                     nTotalPeps = c(4200, 3600, 1200), 
                                     inter.pep.cor = 0, 
                                     nSamples = 3) {
  ####################
  ## run simulation ##
  ####################
  ## specify metadata
  pepNames <- paste0("Peptide", 1:sum(nTotalPeps))
  sampleNames <- paste0("Sample", 1:(2*nSamples))
  group <- c(rep("A", nSamples), rep("B", nSamples))
  contrasts.par <- "B-A"
  pep_mapping_tbl <- data.frame(
    peptide = pepNames,
    protein = c(rep(paste0("Protein", 
                           1:sum(floor(nTotalPeps / nTestPeps))), 
                    times = rep(nTestPeps, times = nTotalPeps / nTestPeps)))
  )
  
  ## specify distribution parameters
  norm.mu <- 0
  norm.sigma <- 1
  
  ## simulate peptide-wise mean and variance
  pepwise.means <- rnorm(sum(nTotalPeps), mean = norm.mu,
                         sd = norm.sigma)
  pepwise.vars <- rep(1, sum(nTotalPeps))
  
  ## simulate data
  AllDat <- do.call('rbind', 
                    lapply(1:length(nTotalPeps), function(i) {
                      do.call('rbind', lapply(1:floor(nTotalPeps[i] / nTestPeps[i]), function(j) {
                        t(MASS::mvrnorm(n = nSamples * 2, 
                                        mu = pepwise.means[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)], 
                                        Sigma = (diag(nTestPeps[i]) * 
                                                   (1-inter.pep.cor) + inter.pep.cor) * 
                                          (sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)]) %*% 
                                             t(sqrt(pepwise.vars[(nTestPeps[i] * (j-1) + 1):(nTestPeps[i] * j)])))))
                      }))
                    }))
  rownames(AllDat) <- pepNames
  colnames(AllDat) <- sampleNames
  
  # ## run protein LIMMA with summation
  sumLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                      contrasts.par = contrasts.par,
                                      group = group,
                                      pep_mapping_tbl = pep_mapping_tbl,
                                      method = "sum", 
                                      logged = TRUE)
  
  ## run protein LIMMA with robust regression
  robRegLimma_output <- AggLimmaWorkflow(dat = AllDat,
                                         contrasts.par = contrasts.par,
                                         group = group,
                                         pep_mapping_tbl = pep_mapping_tbl,
                                         method = "robreg", 
                                         logged = TRUE)
  
  pepSetTest_sd_output <- PepSetTestWorkflow(dat = AllDat, 
                                             contrasts.par = contrasts.par,
                                             group = group,
                                             pep_mapping_tbl = pep_mapping_tbl,
                                             stats = "t", 
                                             correlated = FALSE,
                                             equal.correlation = FALSE,
                                             pepC.estim = "sd",
                                             logged = TRUE)
  
  ## retrieve P values
  proteins <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(protein)
  num.peps <- pepSetTest_sd_output %>% dplyr::arrange(protein) %>% dplyr::pull(NPeps)
  sumLimma.pval <- sumLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  robRegLimma.pval <- robRegLimma_output %>% dplyr::arrange(protein) %>% dplyr::pull(P.Value)
  pepSetTest.sd.pval <- pepSetTest_sd_output %>%
    dplyr::arrange(protein) %>%
    dplyr::pull(PValue)
  
  pval.df <- data.frame( protein = proteins,
                         nTestPep = num.peps,
                         sumLimma = sumLimma.pval,
                         robRegLimma = robRegLimma.pval,
                         pepSetTestSD = pepSetTest.sd.pval )
  pval.df$mu <- 0
  return(pval.df)
}

