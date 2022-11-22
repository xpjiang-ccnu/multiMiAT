#' @export
Whole.Kernels <- function(otu.tab, tree, distance = c("Jaccard", "BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "EK", "LaK"), sigma.GK = 1, sigma.EK = 1, sigma.LaK = 1) {
  
  jac.k <- bc.k <- u.unif.k <- g.unif.00.k <- g.unif.25.k <- g.unif.50.k <- g.unif.75.k <- w.unif.k <- Ks.names <- NULL
  Ks <- list()
  
  if ("Jaccard" %in% distance) {
    pa.otu.tab <- otu.tab
    for (i in 1:nrow(otu.tab)) {
      ind <- which(otu.tab[i,] > 0)
      pa.otu.tab[i,ind] <- 1
    }
    jac <- as.matrix(bcdist(pa.otu.tab))
    sigma.jac.list <- sigma_value(jac)

    if (is.null(sigma.GK)) {
      sigma.jac.GK = sigma.jac.list$sigma.GK
    } else {
      sigma.jac.GK = sigma.GK
    }

    if (is.null(sigma.EK)) {
      sigma.jac.EK = sigma.jac.list$sigma.EK
    } else {
      sigma.jac.EK = sigma.EK
    }

    if (is.null(sigma.LaK)) {
      sigma.jac.LaK = sigma.jac.list$sigma.LaK
    } else {
      sigma.jac.LaK = sigma.LaK
    }
    
    jac.k <- Kernel.Matrix(jac, kernel = kernel, sigma.GK = sigma.jac.GK, sigma.EK = sigma.jac.EK, sigma.LaK = sigma.jac.LaK)
    Ks.names <- c(Ks.names, "Jaccard")
    Ks <- c(Ks, list(jac.k))
  }

  if ("BC" %in% distance) {
    bc <- as.matrix(bcdist(otu.tab))
    sigma.bc.list <- sigma_value(bc)

    if (is.null(sigma.GK)) {
      sigma.bc.GK = sigma.bc.list$sigma.GK
    } else {
      sigma.bc.GK = sigma.GK
    }

    if (is.null(sigma.EK)) {
      sigma.bc.EK = sigma.bc.list$sigma.EK
    } else {
      sigma.bc.EK = sigma.EK
    }

    if (is.null(sigma.LaK)) {
      sigma.bc.LaK = sigma.bc.list$sigma.LaK
    } else {
      sigma.bc.LaK = sigma.LaK
    } 

    bc.k <- Kernel.Matrix(bc, kernel = kernel, sigma.GK = sigma.bc.GK, sigma.EK = sigma.bc.EK, sigma.LaK = sigma.bc.LaK)
    Ks.names <- c(Ks.names, "BC")
    Ks <- c(Ks, list(bc.k))
  }
  
  phydist <- c("U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac")
  if (TRUE %in% (distance %in% phydist)) {
    unifs <- GUniFrac::GUniFrac(otu.tab, tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
    u.unif <- unifs[, , "d_UW"]
    g.unif.00 <- unifs[, , paste("d_", 0, sep = "")]
    g.unif.25 <- unifs[, , paste("d_", 0.25, sep = "")]
    g.unif.50 <- unifs[, , paste("d_", 0.5, sep = "")]
    g.unif.75 <- unifs[, , paste("d_", 0.75, sep = "")]
    w.unif <- unifs[, , "d_1"]
    for (j in 1:nrow(otu.tab)) {
      ind <- is.na(u.unif[j, ])
      
      if (sum(ind) != 0) u.unif[, ind] <- 0
      ind <- is.na(g.unif.00[j, ])
      
      if (sum(ind) != 0) g.unif.00[, ind] <- 0
      ind <- is.na(g.unif.25[j, ])
      
      if (sum(ind) != 0) g.unif.25[, ind] <- 0
      ind <- is.na(g.unif.50[j, ])
      
      if (sum(ind) != 0) g.unif.50[, ind] <- 0
      ind <- is.na(g.unif.75[j, ])
      
      if (sum(ind) != 0) g.unif.75[, ind] <- 0
      ind <- is.na(w.unif[j, ])
      
      if (sum(ind) != 0) w.unif[, ind] <- 0
    }
    
    if ("U.UniFrac" %in% distance){
      sigma.u.unif.list <- sigma_value(u.unif)

      if (is.null(sigma.GK)) {
        sigma.u.unif.GK = sigma.u.unif.list$sigma.GK
      } else {
        sigma.u.unif.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.u.unif.EK = sigma.u.unif.list$sigma.EK
      } else {
        sigma.u.unif.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.u.unif.LaK = sigma.u.unif.list$sigma.LaK
      } else {
        sigma.u.unif.LaK = sigma.LaK
      }
      
      u.unif.k <- Kernel.Matrix(u.unif, kernel = kernel, sigma.GK = sigma.u.unif.GK, sigma.EK = sigma.u.unif.EK, sigma.LaK = sigma.u.unif.LaK)
      Ks.names <- c(Ks.names, "U.UniFrac")
      Ks <- c(Ks, list(u.unif.k))
    }
    
    if ("G.UniFrac.00" %in% distance){
      sigma.g.unif.00.list <- sigma_value(g.unif.00)

      if (is.null(sigma.GK)) {
        sigma.g.unif.00.GK = sigma.g.unif.00.list$sigma.GK
      } else {
        sigma.g.unif.00.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.g.unif.00.EK = sigma.g.unif.00.list$sigma.EK
      } else {
        sigma.g.unif.00.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.g.unif.00.LaK = sigma.g.unif.00.list$sigma.LaK
      } else {
        sigma.g.unif.00.LaK = sigma.LaK
      }
      
      g.unif.00.k <- Kernel.Matrix(g.unif.00, kernel = kernel, sigma.GK = sigma.g.unif.00.GK, sigma.EK = sigma.g.unif.00.EK, sigma.LaK = sigma.g.unif.00.LaK)
      Ks.names <- c(Ks.names, "G.UniFrac.00")
      Ks <- c(Ks, list(g.unif.00.k))
    }
    
    if ("G.UniFrac.25" %in% distance){
      sigma.g.unif.25.list <- sigma_value(g.unif.25)

      if (is.null(sigma.GK)) {
        sigma.g.unif.25.GK = sigma.g.unif.25.list$sigma.GK
      } else {
        sigma.g.unif.25.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.g.unif.25.EK = sigma.g.unif.25.list$sigma.EK
      } else {
        sigma.g.unif.25.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.g.unif.25.LaK = sigma.g.unif.25.list$sigma.LaK
      } else {
        sigma.g.unif.25.LaK = sigma.LaK
      }
      
      g.unif.25.k <- Kernel.Matrix(g.unif.25, kernel = kernel, sigma.GK = sigma.g.unif.25.GK, sigma.EK = sigma.g.unif.25.EK, sigma.LaK = sigma.g.unif.25.LaK)
      Ks.names <- c(Ks.names, "G.UniFrac.25")
      Ks <- c(Ks, list(g.unif.25.k))
    }
    
    if ("G.UniFrac.50" %in% distance){
      sigma.g.unif.50.list <- sigma_value(g.unif.50)

      if (is.null(sigma.GK)) {
        sigma.g.unif.50.GK = sigma.g.unif.50.list$sigma.GK
      } else {
        sigma.g.unif.50.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.g.unif.50.EK = sigma.g.unif.50.list$sigma.EK
      } else {
        sigma.g.unif.50.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.g.unif.50.LaK = sigma.g.unif.50.list$sigma.LaK
      } else {
        sigma.g.unif.50.LaK = sigma.LaK
      }

      g.unif.50.k <- Kernel.Matrix(g.unif.50, kernel = kernel, sigma.GK = sigma.g.unif.50.GK, sigma.EK = sigma.g.unif.50.EK, sigma.LaK = sigma.g.unif.50.LaK)
      Ks.names <- c(Ks.names, "G.UniFrac.50")
      Ks <- c(Ks, list(g.unif.50.k))
    }
    
    if ("G.UniFrac.75" %in% distance){
      sigma.g.unif.75.list <- sigma_value(g.unif.75)

      if (is.null(sigma.GK)) {
        sigma.g.unif.75.GK = sigma.g.unif.75.list$sigma.GK
      } else {
        sigma.g.unif.75.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.g.unif.75.EK = sigma.g.unif.75.list$sigma.EK
      } else {
        sigma.g.unif.75.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.g.unif.75.LaK = sigma.g.unif.75.list$sigma.LaK
      } else {
        sigma.g.unif.75.LaK = sigma.LaK
      }

      g.unif.75.k <- Kernel.Matrix(g.unif.75, kernel = kernel, sigma.GK = sigma.g.unif.75.GK, sigma.EK = sigma.g.unif.75.EK, sigma.LaK = sigma.g.unif.75.LaK)
      Ks.names <- c(Ks.names, "G.UniFrac.75")
      Ks <- c(Ks, list(g.unif.75.k))
    }
    
    if ("W.UniFrac" %in% distance){
      sigma.w.unif.list <- sigma_value(w.unif)

      if (is.null(sigma.GK)) {
        sigma.w.unif.GK = sigma.w.unif.list$sigma.GK
      } else {
        sigma.w.unif.GK = sigma.GK
      }

      if (is.null(sigma.EK)) {
        sigma.w.unif.EK = sigma.w.unif.list$sigma.EK
      } else {
        sigma.w.unif.EK = sigma.EK
      }

      if (is.null(sigma.LaK)) {
        sigma.w.unif.LaK = sigma.w.unif.list$sigma.LaK
      } else {
        sigma.w.unif.LaK = sigma.LaK
      }

      w.unif.k <- Kernel.Matrix(w.unif, kernel = kernel, sigma.GK = sigma.w.unif.GK, sigma.EK = sigma.w.unif.EK, sigma.LaK = sigma.w.unif.LaK)
      Ks.names <- c(Ks.names, "W.UniFrac")
      Ks <- c(Ks, list(w.unif.k))
    }
  }
  
  names(Ks) <- Ks.names
  return(Ks)
}
