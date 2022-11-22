## ----multiMiAT, include = FALSE-----------------------------------------------
multiMiAT <- function (y, otu.tab, tree, covs = NULL, method = c("multiMiRKAT", "score", "MiRKAT-MC"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W, CLR = FALSE, n.perm = 999, seed = 1) {
  
  if (CLR) {
    otu.tab <- t(apply(otu.tab, 1, clr))
  }

  if(class(covs) != "data.frame") {
    stop("Error: data format of covariates must be data frame!")
  }

  Ks <- Whole.Kernels(otu.tab, tree, distance = distance, kernel = kernel) 
  
  if (!("multiMiRKAT" %in% method) & !("score" %in% method) & !("MiRKAT-MC" %in% method)) {
    stop("Error: Method must have at least one of multiMiRKAT, score and MiRKAT-MC! Method can also be any two or all of them! Please enter the correct method!")
  }

  multiMiAT.all <- list()
  multiMiAT.pvs <- NULL
  if ("multiMiRKAT" %in% method){
    re.bcl <- multiMiRKAT_bcl(y = y, Ks = Ks, covs = covs, n.perm = n.perm, seed = seed)
    re.clm <- multiMiRKAT_clm(y = y, Ks = Ks, covs = covs, n.perm = n.perm, seed = seed)
    multiMiRKAT.pvs <- c(list(multiMiRKAT.bcl.pvs = re.bcl$multiMiRKAT.bcl.pvs), list(multiMiRKAT.clm.pvs = re.clm$multiMiRKAT.clm.pvs))
    multiMiRKAT <- c(re.bcl$multiMiRKAT.bcl.Omnibus.p["multiMiRKAT-N"], re.clm$multiMiRKAT.clm.Omnibus.p["multiMiRKAT-O"])
    names(multiMiRKAT.pvs) <- c("multiMiRKAT.bcl.pvs", "multiMiRKAT.clm.pvs")
    names(multiMiRKAT) <- c("multiMiRKAT-N", "multiMiRKAT-O")
    multiMiAT.all <- c(multiMiAT.all, list(multiMiRKAT.pvs))
    names(multiMiAT.all)[length(multiMiAT.all)] <- "multiMiRKAT.pvs"
    multiMiAT.pvs <- c(multiMiAT.pvs, multiMiRKAT)
  }

  if ("score" %in% method){
    score.test <- as.numeric(round(score.test_bcl(y, otu.tab, covs, W), digits = 4))
    names(score.test) <- "score test"
    multiMiAT.pvs <- c(multiMiAT.pvs, score.test)
  }
  
  if ("MiRKAT-MC" %in% method){
    MiRKATMC.bcl <- MiRKAT_MCN(y, Ks, covs = covs)
    MiRKATMC.clm <- MiRKAT_MCO(y, Ks, covs = covs)
    MiRKAT.MC <- c(round(MiRKATMC.bcl$MiRKAT_MCN, digits = 4), round(MiRKATMC.clm$MiRKAT_MCO, digits = 4))
    names(MiRKAT.MC) <- c("MiRKAT-MCN", "MiRKAT-MCO")
    multiMiAT.pvs <- c(multiMiAT.pvs, MiRKAT.MC)
  }
  
  multiMiAT.pvs[which(multiMiAT.pvs == 0)] <- 1e-6
  if (length(method) == 1){
    if (method != "score") {
      multiMiAT.pvs <- c(multiMiAT.pvs, round(harmonicmeanp::p.hmp(p = multiMiAT.pvs, w = NULL, L = length(multiMiAT.pvs), w.sum.tolerance = 1e-6, multilevel = F)[1], digits = 4))
      names(multiMiAT.pvs)[length(multiMiAT.pvs)] <- method
    }
    multiMiAT.all <- c(multiMiAT.all, list(multiMiAT.pvs))
    names(multiMiAT.all)[length(multiMiAT.all)] <- "Combined.methods.pvs"
  } else {
    multiMiAT.omni <- round(harmonicmeanp::p.hmp(p = multiMiAT.pvs, w = NULL, L = length(multiMiAT.pvs), w.sum.tolerance = 1e-6, multilevel = F)[1], digits = 4)
    if (length(method) == 2) {
      names(multiMiAT.omni) <- paste(method[1], "+", method[2], sep = "")
    } else {
      names(multiMiAT.omni) <- "multiMiAT"
    }
    multiMiAT.all <- c(multiMiAT.all, list(multiMiAT.pvs), list(multiMiAT.omni))    
    names(multiMiAT.all)[(length(multiMiAT.all)-1):length(multiMiAT.all)] <- c("Combined.methods.pvs", "Optimal.test.p") 
  }

  return(multiMiAT = multiMiAT.all)
}

