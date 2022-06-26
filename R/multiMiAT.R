#' @export
multiMiAT <- function (y, otu.tab, tree, covs = NULL, model = c("bcl", "clm"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W, n.perm = 999, seed = 1) {
  
  Ks <- Whole.Kernels(otu.tab, tree, distance = distance, kernel = kernel) 
  
  if (!("bcl" %in% model) & !("clm" %in% model) & !("crm" %in% model)) {
    stop("Error: Model must have at least one of bcl, clm and crm! And model can be two or three! Please enter the correct omnibus!")
  }

  multiMiRKAT.pvs <- list()
  multiMiRKAT <- MiRKAT.MC <- NULL
  if ("bcl" %in% model){
    re.bcl <- multiMiRKAT_bcl(y = y, Ks = Ks, covs = covs, n.perm = n.perm, seed = seed)
    multiMiRKAT.pvs <- c(multiMiRKAT.pvs, list(multiMiRKAT.bcl.pvs = re.bcl$multiMiRKAT.bcl.pvs))
    multiMiRKAT <- c(multiMiRKAT, re.bcl$multiMiRKAT.bcl.Omnibus.p["multiMiRKAT-N"])
    names(multiMiRKAT)[length(multiMiRKAT)] <- "multiMiRKAT-N"
    
    out.bcl <- MiRKAT_MCN(y, Ks, covs = covs)
    MiRKAT.MC <- c(MiRKAT.MC, round(out.bcl$MiRKAT_MCN, digits = 4))
    names(MiRKAT.MC)[length(MiRKAT.MC)] <- "MiRKAT-MCN"
  }
  
  if ("clm" %in% model){
    re.clm <- multiMiRKAT_clm(y = y, Ks = Ks, covs = covs, n.perm = n.perm, seed = seed)
    multiMiRKAT.pvs <- c(multiMiRKAT.pvs, list(multiMiRKAT.clm.pvs = re.clm$multiMiRKAT.clm.pvs))
    multiMiRKAT <- c(multiMiRKAT, re.clm$multiMiRKAT.clm.Omnibus.p["multiMiRKAT-O"])
    names(multiMiRKAT)[length(multiMiRKAT)] <- "multiMiRKAT-O"
    
    out.clm <- MiRKAT_MCO(y, Ks, covs = covs)
    MiRKAT.MC <- c(MiRKAT.MC, round(out.clm$MiRKAT_MCO, digits = 4))
    names(MiRKAT.MC)[length(MiRKAT.MC)] <- "MiRKAT_MCO"
  }
    
  if ("crm" %in% model){
    re.crm <- multiMiRKAT_crm(y = y, Ks = Ks, covs = covs, n.perm = n.perm, seed = seed)
    multiMiRKAT.pvs <- c(multiMiRKAT.pvs, list(multiMiRKAT.crm.pvs = re.crm$multiMiRKAT.crm.pvs))
    multiMiRKAT <- c(multiMiRKAT, re.crm$multiMiRKAT.crm.Omnibus.p["multiMiRKAT.crm"])
    names(multiMiRKAT)[length(multiMiRKAT)] <- "multiMiRKAT-crm"
  }

  score.test <- as.numeric(round(score.test_bcl(y, otu.tab, covs, W), digits = 4))
  names(score.test) <- "score.test"
  
  multiMiAT.pvs <- c(multiMiRKAT, score.test, MiRKAT.MC)
  multiMiAT.pvs[which(multiMiAT.pvs == 0)] <- 1e-6
  multiMiAT <- round(harmonicmeanp::p.hmp(p = multiMiAT.pvs, w = NULL, L = length(multiMiAT.pvs), w.sum.tolerance = 1e-6, multilevel = F)[1], digits = 4)
  names(multiMiAT) <- "multiMiAT"
  
  multiMiAT <- c(list(multiMiRKAT.pvs$multiMiRKAT.bcl.pvs), list(multiMiRKAT.pvs$multiMiRKAT.clm.pvs), list(score.test), list(MiRKAT.MC), list(multiMiRKAT), list(multiMiAT))
  names(multiMiAT) <- c("multiMiRKAT.bcl.pvs", "multiMiRKAT.clm.pvs", "Score.test", "MiRKAT.MC", "multiMiRKAT", "multiMiAT" )
  
  return(multiMiAT = multiMiAT)
  
}

