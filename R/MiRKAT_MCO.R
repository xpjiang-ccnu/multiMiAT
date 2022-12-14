#' @export
MiRKAT_MCO <- function(y, Ks, covs = NULL){
  
  if (!is.matrix(y) & !is.vector(y) & !is.integer(y) & !is.factor(y)) {
    stop("Error: y must be a matrix, vector, integer or factor!")
  } 
  
  if (is.matrix(y)){
    if(any(rowSums(y)!=1)) {
      stop("Error: summation of each row must equal to 1!")
    }  
    y.inter <- which(y == 1, arr.ind = T)[,2]
  }
  
  if (is.vector(y)){
    if(!is.null(covs) && length(y) != nrow(covs)) {
      stop("Error: The dimension of the response variable is not equal to that of the covariable!")
    } else {
      y.inter <- as.integer(y)
    }
  }
  
  if (is.integer(y)) {
    y.inter <- y
  }
  
  y.ordered <- as.ordered(y)
  
  if (is.null(covs)) {
    formula = y ~ 1
  } else {
    xnam <- colnames(covs)
    formula <- as.formula(paste("y.ordered ~ ", paste(xnam, collapse= "+")))
  }
  
  inner.POM.null <- function(formula, data = covs){
    fit <- MASS::polr(formula, data, method = 'logistic', na.action = 'na.fail')
    if (is.null(data)) data <- fit$model
    X <- model.matrix(formula, data)
    nSam <- dim(X)[1]
    y <- as.formula(paste0('~ 0 +', as.character(formula)[[2]]))
    y.matrix <- model.matrix(y, data)
    J <- ncol(y.matrix)
    
    pi.hat <- fit$fitted.values
    A <- matrix(1, nrow = J-1, ncol = J-1)
    A[upper.tri(A)] <- 0
    A <- kronecker(A, Diagonal(nSam, 1))
    
    mu.hat <- matrix(A %*% c(pi.hat[, 1: (J-1)]), ncol = J-1)
    y.tilde <- matrix(A %*% c(y.matrix[, 1: (J-1)]), ncol = J-1)
    KY <- tcrossprod(y.tilde - mu.hat)
    
    return(KY)
  }
  
  KY <- inner.POM.null(formula, covs)

  if (is.matrix(Ks)){
    Ks <- list(Ks)
  } else if (!is.list(Ks)){
    stop("Please enter either a single kernel matrix or a list of kernels for Ks.")
  }

  DKAT=function(K,L){
    tr=function(x){return(sum(diag(x)))}
    n=nrow(K)
    I.n=diag(1,n)
    I.1=rep(1,n)
    H=I.n-I.1%*%t(I.1)/n
    K=H%*%K%*%H
    L=H%*%L%*%H
    A=K/tr(K%*%K)  ## standard-version of K
    W=L/tr(L%*%L)
    
    Fstar=tr(A%*%W)
    mean.krv=tr(A)*tr(W)/(n-1)	## mean of DKAT
    
    T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
    Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
    temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
    temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
    temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
    temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
    temp2=temp21*temp22/temp23
    variance.krv=temp1+temp2		## variance of DKAT
    
    T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
    T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
    t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
    t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
    t3=24*(n^2-n-4)*(U*Bs+B*Us)
    t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
    t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
    t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
    t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
    t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
    t8=24*(t81+t82)
    t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
    t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
    t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
    t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
    t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
    t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
    t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
    t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
    t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
    t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
    t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
    t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
    t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
    t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
    t20=-(n-2)*(t201+t202+t203)
    temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
    temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
    mom3=temp31/temp32
    skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of DKAT
    
    m1=mean.krv
    m2=variance.krv
    m3=skewness.krv
    shape=4/m3^2
    scale=sqrt(m2)*m3/2
    location=m1-2*sqrt(m2)/m3
    PIIIpars=list(shape,location,scale)
    pv=1- PearsonDS::ppearsonIII(Fstar, params=PIIIpars)
    # return(list(pv, Fstar, PIIIpars))
    if (is.na(pv)) {pv=1}  ## usually happens if var=0 in the denominator
    return(pv)
  }
  
  pvs <- NULL
  for (i in 1:length(Ks[[1]])) {
    list.kernels <- list()
    for (ii in 1:length(Ks)) { 
      list.kernels <- c(list.kernels, list(Ks[[ii]][[i]]))
    }
    pvs.single.kernel <- mapply(function(x){DKAT(x, KY)}, list.kernels)
    names(pvs.single.kernel) <- names(Ks)
    pvs <- rbind(pvs, pvs.single.kernel)
  }
  
  rownames(pvs) <- names(Ks[[1]])
  pvs.all <- c(pvs)
  pvs.all[which(pvs.all == 0)] <- 1e-6
  p.MiRKAT_MCO <- harmonicmeanp::p.hmp(p = pvs.all, w = NULL, L = length(pvs.all), w.sum.tolerance = 1e-6, multilevel = F)[1]
  
  MiRKAT_MCO <- list(MiRKAT_MCO.pvs = pvs, MiRKAT_MCO = p.MiRKAT_MCO)
  
  return(MiRKAT_MCO)
}

