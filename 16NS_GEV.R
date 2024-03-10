
# Libraries
library(data.table)
library(reshape2)
library(dplyr)
library(SpatialExtremes)
library(ismev)
library(stringr)

# Important functions
# for real data
rm_tie <- function(x){
  if(is.data.frame(x)){
    x0 <- x[,1]
  }else{
    x0 <- x
  }
  x <- as.numeric(x0)
  tie_index <- which(duplicated(x))
  
  if(length(tie_index)>0){
    x[tie_index] <- x[tie_index] + 0.0001
    
    recheck_tie <- which(duplicated(x))
    if(length(recheck_tie)>0){
      data_rm_tie <- rm_tie(x)
    }else{
      data_rm_tie <- x
    }
    return(data_rm_tie)
  }else return(x)
}

# transform function for K-S test
trans.GEV <- function(z_t,mu,sig,xi){
  term <- (z_t-mu)/sig
  z_t_tilda <- (1/xi) * log(1 + xi * term)
  
  return(z_t_tilda)
}

trans.GUM <- function(z_t,mu,sig){
  z_t_tilda <- (z_t - mu)/sig
  return(z_t_tilda)
}

# for parameter outputs
fillNA <- function(mu,sig,xi){
  mu_v <- c(mu,rep(NA,3-length(mu)))
  sig_v <- c(sig,rep(NA,2-length(sig)))
  xi_v <- c(xi,rep(NA,2-length(xi)))
  return(as.numeric(c(mu_v,sig_v,xi_v)))
}

gev.llh <- function(par, x, t, model) {
  if(model=="GEV300"){
    b0 <- par[1] ; b1 <- par[2] ; b2 <- par[3]
    mu <- function(t) {b0+b1*exp(-b2*t)}
    sig <- par[4]  ;  xi<- par[5]
    m <- length(x)
    y <- (x-mu(t))/sig
    if(sig <= 0 || any((1+xi*y)<=0)) return(10^6)
    fn_value <- if(xi!=0){ -m*log(sig)-((1/xi)+1)*sum(log(1+xi*y))-sum((1+xi*y)^(-1/xi)) }else{ return(10^6) } 
    
  }else if(model=="GEV301"){
    b0 <- par[1] ; b1 <- par[2] ; b2 <- par[3]
    c0 <- par[5] ; c1 <- par[6]
    mu <- function(t) {b0+b1*exp(-b2*t)}
    sig <- par[4]  
    xi <- function(t) {c0+c1*t}
    m <- length(x)
    y <- (x-mu(t))/sig
    if(sig <= 0 || any((1+xi(t)*y)<=0)) return(10^6)
    fn_value <- if(any(xi(t)!=0)){ -m*log(sig)-((1/xi(t))+1)*sum(log(1+xi(t)*y))-sum((1+xi(t)*y)^(-1/xi(t))) }else{ return(10^6) } 
    
  }else if(model=="GEV310"){ 
    b0 <- par[1] ;  b1 <- par[2] ;  b2 <- par[3]
    mu <- function(t) {b0+b1*exp(-b2*t)}
    sig0 <- par[4] ;  sig1 <- par[5]
    sig <- function(t) {exp(sig0+sig1*t)}
    xi<- par[6]
    m <- length(x)
    y <- (x-mu(t))/sig(t)
    if(any(sig(t) <= 0) || any((1+xi*y)<=0)) return(10^6)
    fn_value <- if(xi!=0){-sum(log(sig(t))+(1+(1/xi))*log(1+(xi*y))+(1+(xi*y))^(-1/xi)) }else{ return(10^6) }
    
  }else if(model=="GEV311"){
    a0 <- par[1] ; a1 <- par[2] ; a2 <- par[3]
    b0 <- par[4] ; b1 <- par[5]
    c0 <- par[6] ; c1 <- par[7]
    
    mu <- function(t) {a0+a1*exp(-a2*t)}
    sig<- function(t) {exp(b0+(b1*t))}
    xi<- function(t) {c0+(c1*t)}
    m <- length(x)
    y <- (x-mu(t))/sig(t)
    
    if(any(sig(t) <= 0) || any((1+xi(t)*y)<=0)) return(10^6)
    
    fn_value <-
      if(any(xi(t)!=0)){-sum(log(sig(t))+(1+(1/xi(t)))*log(1+(xi(t)*y))+(1+(xi(t)*y))^(-1/xi(t)))
      }else{ return(10^6)
      }
    
  }else if(model=="GD300"){
    b0 <- par[1]  ; b1 <- par[2] ; b2 <- par[3]
    mu <- function(t) {b0+b1*exp(-b2*t)}
    sig <- par[4]
    m <- length(x)
    if(sig<=0) return(10^6)
    y <- (x-mu(t))/sig
    fn_value <- -m*log(sig)-sum(y)-sum(exp(-y))
    
  }else if(model=="GD310"){
    b0 <- par[1] ; b1 <- par[2] ; b2 <- par[3]
    mu <- function(t) {b0+b1*exp(-b2*t)}
    sig0 <- par[4] ;  sig1 <- par[5]
    sig <- function(t) {exp(sig0+sig1*t)}
    m <- length(x)
    if(any(sig(t)<=0)) return(10^6)
    y <- (x-mu(t))/sig(t)
    fn_value <- -sum(log(sig(t))+y+(exp(-y)))
  }
  return(-fn_value)
} 

initial_random_func <- function(a, t = NULL, model){
  t <- t
  if(is.data.frame(a)){
    prec <- a[,1]
  }else{
    prec <- a
  }
  n_int_random <- 20
  
  if(model=="GEV000"){
    est_int <- gev.fit(xdat=prec, show = F)$mle      
    mu_int <- est_int[1]+rnorm(n_int_random,0,1)
    sig_int<- est_int[2]+rnorm(n_int_random,0,0.1)
    xi_int <- est_int[3]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int, sig_int, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, MARGIN=1, FUN=function(x, xdat){
      m <- gev.fit(xdat, muinit=x[1], siginit=x[2], shinit=x[3], show=FALSE)
      result <- c(nllh=m$nllh, x[1], x[2], x[3])
      return(result) 
    })
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[4,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV100"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, show = FALSE)$mle      
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)     
    sig_int<- est_int[3]+ rnorm(n_int_random,0,0.1)
    xi_int <- est_int[4]+rnorm(n_int_random,0,0.1)      
    int_df <- cbind(mu_int0, mu_int1, sig_int, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=1, muinit=x[1:2], siginit=x[3], shinit=x[4], show=FALSE), error=function(e){}) 
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV200"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=c(1,2), show = FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1) 
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)  
    mu_int2 <- est_int[3]+rnorm(n_int_random,0,0.1)
    sig_int<- est_int[4]+ rnorm(n_int_random,0,0.1)
    xi_int <- est_int[5]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=c(1,2), muinit=x[1:3], siginit=x[4], shinit=x[5], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV300"){
    est_int <- gev.fit(xdat=prec, show = FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- rnorm(n_int_random,0,0.5)
    mu_int2 <- rnorm(n_int_random,0,0.1)
    sig_int <- est_int[2]+rnorm(n_int_random,0,0.1)
    xi_int  <- est_int[3]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){}) 
      result <- c(nllh=m$value, x[1], x[2], x[3], x[4], x[5])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV010"){
    est_int <- gev.fit(xdat=prec, ydat=t, sigl=1, siglink=exp, show = FALSE)$mle      
    mu_int <- est_int[1]+rnorm(n_int_random,0,1)    
    sig_int0<- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[3]+rnorm(n_int_random,0,0.1) 
    xi_int <- est_int[4]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int, sig_int0, sig_int1, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, sigl=1, siglink=exp, muinit=x[1], siginit=x[2:3], shinit=x[4], show=FALSE), error=function(e){})
      result <- c(m$nllh, x[1], x[2], x[3], x[4])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV110"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)  
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int <- est_int[5]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, sig_int0, sig_int1, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=1, sigl=1, siglink=exp, muinit=x[1:2], siginit=x[3:4], shinit=x[5], show=FALSE), error=function(e){})
      result <- c(m$nllh, x[1], x[2], x[3], x[4], x[5])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])] 
    
  }else if(model=="GEV210"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1) 
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.01) 
    mu_int2 <- rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.01)
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.01)
    xi_int <- est_int[5]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, t, mul=c(1,2), sigl=1, siglink=exp, muinit=x[1:3], siginit=x[4:5], shinit=x[6], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5], x[6])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[7,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])] 
    
  }else if(model=="GEV310"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1) 
    mu_int2 <- rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,1) 
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.1) 
    xi_int  <- est_int[5]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1, xi_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){})
      result <- c(m$value, x[1], x[2], x[3], x[4], x[5], x[6])
      return(result)
    }) 
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[7,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV001"){
    est_int <- gev.fit(xdat=prec,ydat=t,shl=1,show = FALSE)$mle    
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    sig_int0 <- est_int[2]+ rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[3]+rnorm(n_int_random,0,0.01)
    xi_int1 <- est_int[4]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, sig_int0, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, shl=1, muinit=x[1], siginit=x[2], shinit=x[3:4], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV101"){
    est_int <- gev.fit(xdat=prec,ydat=t,mul=1,shl=1,show = FALSE)$mle    
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int0 <- est_int[3]+rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[5]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, mu_int1, sig_int0, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=1, shl=1, muinit=x[1:2], siginit=x[3], shinit=x[4:5], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV201"){
    est_int <- gev.fit(xdat=prec,ydat=t,mul=c(1,2),shl=1,show = FALSE)$mle    
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    mu_int2 <- est_int[3]+rnorm(n_int_random,0,0.01)
    sig_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[5]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[6]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=c(1,2), shl=1, muinit=x[1:3], siginit=x[4], shinit=x[5:6], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5], x[6])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[7,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV301"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, shl=1, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- rnorm(n_int_random,0,0.01)
    mu_int2 <- est_int[2]+rnorm(n_int_random,0,0.01)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[5]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){}) 
      result <- c(nllh=m$value, x[1], x[2], x[3], x[4], x[5], x[6])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[7,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV011"){
    est_int <- gev.fit(xdat=prec, ydat=t, sigl=1, shl=1, siglink=exp, show = FALSE)$mle      
    mu_int <- est_int[1]+rnorm(n_int_random,0,1)    
    sig_int0<- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[3]+rnorm(n_int_random,0,0.1) 
    xi_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[5]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, sig_int0, sig_int1, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, sigl=1, shl=1, siglink=exp, muinit=x[1], siginit=x[2:3], shinit=x[4:5], show=FALSE), error=function(e){})
      result <- c(m$nllh, x[1], x[2], x[3], x[4], x[5])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GEV111"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=1, sigl=1, shl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)  
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[5]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[6]+rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0,mu_int1, sig_int0, sig_int1, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, ydat=t, mul=1, sigl=1, shl=1, siglink=exp, muinit=x[1:2], siginit=x[3:4], shinit=x[5:6], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5], x[6])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[7,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])] 
    
  }else if(model=="GEV211"){
    est_int <- gev.fit(xdat=prec, ydat=t, mul=c(1,2), sigl=1, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    mu_int2 <- est_int[3]+rnorm(n_int_random,0,0.01)
    sig_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    sig_int1 <- est_int[5]+rnorm(n_int_random,0,0.01)
    xi_int0 <- est_int[6]+rnorm(n_int_random,0,0.01)
    xi_int1 <- rnorm(n_int_random,0,0.01)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gev.fit(xdat, t, mul=c(1,2), sigl=1, shl=1, siglink=exp, muinit=x[1:3], siginit=x[4:5], shinit=x[6:7], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5], x[6], x[7])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[8,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])] 
    
  }else if(model=="GEV311"){
    est_int <- gev.fit(xdat=prec, ydat=t, sigl=1, shl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- rnorm(n_int_random,0,1)
    mu_int2 <- rnorm(n_int_random,0,0.01)
    sig_int0<- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[3]+rnorm(n_int_random,0,0.1)
    xi_int0 <- est_int[4]+rnorm(n_int_random,0,0.1)
    xi_int1 <- est_int[5]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1, xi_int0, xi_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){})
      result <- c(m$value, x[1], x[2], x[3], x[4], x[5], x[6], x[7])
      return(result)
    }) 
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[8,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD000"){
    est_int <- gum.fit(xdat=prec, show = F)$mle       
    mu_int <- est_int[1]+rnorm(n_int_random,0,1)
    sig_int <- est_int[2]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int, sig_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, MARGIN=1, FUN=function(x, xdat){
      m <- gum.fit(xdat, muinit=x[1], siginit=x[2], show=FALSE)
      result <- c(nllh=m$nllh, x[1], x[2])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[3,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD100"){
    est_int <- gum.fit(xdat=prec, ydat=t, mul=1,show = F)$mle   
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)     
    sig_int<- est_int[3]+ rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, sig_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gum.fit(xdat, t, mul=1, muinit=x[1:2], siginit=x[3], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[4,]))]
    min.nllh.int  <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD200"){
    est_int <- gum.fit(xdat=prec,ydat=t,mul=c(1,2),show = F)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)    
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)  
    mu_int2 <- est_int[3]+rnorm(n_int_random,0,0.1)  
    sig_int<- est_int[4]+ rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gum.fit(xdat, t, mul=c(1,2), muinit=x[1:3], siginit=x[4], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD300"){
    est_int <- gum.fit(xdat=prec, show = FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1) 
    mu_int1 <- rnorm(n_int_random,0,0.1)
    mu_int2 <- rnorm(n_int_random,0,0.1)         
    sig_int <- est_int[2]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){}) 
      result <- c(nllh=m$value, x[1], x[2], x[3], x[4])
      return(result) 
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD010"){
    est_int <- gum.fit(xdat=prec, ydat=t, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int <- est_int[1]+rnorm(n_int_random,0,1)     
    sig_int0<- est_int[2]+ rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[3]+ rnorm(n_int_random,0,0.1) 
    int_df <- cbind(mu_int, sig_int0, sig_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gum.fit(xdat, t, sigl=1, siglink=exp, muinit=x[1], siginit=x[2:3], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[4,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD110"){
    est_int <- gum.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1)  
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.1)
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, sig_int0, sig_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gum.fit(xdat, t, mul=1, sigl=1, siglink=exp, muinit=x[1:2], siginit=x[3:4], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[5,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD210"){
    est_int <- gum.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,0.1) 
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1) 
    mu_int2 <- rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,0.1) 
    sig_int1<- est_int[4]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, MARGIN=1, FUN=function(x, xdat, t){
      m <- tryCatch(gum.fit(xdat, t, mul=c(1,2), sigl=1, siglink=exp, muinit=x[1:3], siginit=x[4:5], show=FALSE), error=function(e){})
      result <- c(nllh=m$nllh, x[1], x[2], x[3], x[4], x[5])
      return(result)
    })
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
    
  }else if(model=="GD310"){
    est_int <- gum.fit(xdat=prec, ydat=t, mul=1, sigl=1, siglink=exp, show=FALSE)$mle
    mu_int0 <- est_int[1]+rnorm(n_int_random,0,1) 
    mu_int1 <- est_int[2]+rnorm(n_int_random,0,0.1) 
    mu_int2 <- rnorm(n_int_random,0,0.1)
    sig_int0<- est_int[3]+rnorm(n_int_random,0,1) 
    sig_int1 <- est_int[4]+rnorm(n_int_random,0,0.1)
    int_df <- cbind(mu_int0, mu_int1, mu_int2, sig_int0, sig_int1)
    
    cal_int_all_n <- apply(X=int_df, xdat=prec, t, model, MARGIN=1, FUN=function(x, xdat, t, model){
      int_value <- x
      m <- tryCatch(optim(par=int_value, fn=gev.llh, x=xdat, t=t, model=model, method="BFGS", hessian=T), error=function(e){})
      result <- c(nllh=m$value, x[1], x[2], x[3], x[4], x[5])
      return(result)
    }) 
    if(is.list(cal_int_all_n)) cal_int_all_n <- do.call(cbind, cal_int_all_n)
    cal_int_all_n <- cal_int_all_n[,which(!is.na(cal_int_all_n[6,]))]
    min.nllh.int <- cal_int_all_n[,which.min(cal_int_all_n[1,])]
  }
  return(min.nllh.int)
}

# Model
GEV000 <- function(a, official_int=NULL){
  # set momdel name
  model <- 'GEV000'
  # set initial value
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it){
      tryCatch({maxitCnt=maxitCnt+1;official_int <- initial_random_func(a=a, model=model)[-1]
      Fal <- T }, error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  try({gev.test <- gev.fit(xdat=a, muinit=official_int[1], siginit=official_int[2], shinit=official_int[3], show=FALSE)}, silent=T)
  nllh <- gev.test$nllh     ;n <- length(a)         ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1]     ;sig <- gev.test$mle[2]      ;xi <- gev.test$mle[3]
  se_mu <- gev.test$se[1]   ;se_sig <- gev.test$se[2]    ;se_xi <- gev.test$se[3]
  mu_int <- official_int[1] ;sig_int <- official_int[2]  ;xi_int <- official_int[3]
  BIC <- 2*nllh+log(n)*k    ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  KS <- ks.test(z_t,"pgev",loc=mu,scale=sig,shape=xi)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV100 <- function(a, official_int=NULL){
  # set model name
  model='GEV100'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1)   
  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gev.fit(xdat=a, ydat=ti, mul=1, muinit=official_int[1:2], siginit=official_int[3], shinit=official_int[4], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gev.fit(xdat=a, ydat=ti, mul=1, show=FALSE)
    official_int[1:4] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)                ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:2]     ;sig <- gev.test$mle[3]      ;xi <- gev.test$mle[4]
  se_mu <- gev.test$se[1:2]   ;se_sig <- gev.test$se[3]    ;se_xi <- gev.test$se[4]
  mu_int <- official_int[1:2] ;sig_int <- official_int[3]  ;xi_int <- official_int[4]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*ti[,1]),sig,xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV200 <- function(a, official_int=NULL){
  # set model name
  model = 'GEV200'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=2)   ;  ti[,1] <- seq(1:length(a))     ; ti[,2] <- ti[,1]^2
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <-  initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gev.fit(xdat=a, ydat=ti, mul=c(1,2), muinit=official_int[1:3], siginit=official_int[4], shinit=official_int[5], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gev.fit(xdat=a, ydat=ti, mul=c(1,2), show=FALSE)
    official_int[1:5] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)         ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:3]     ;sig <- gev.test$mle[4]      ;xi <- gev.test$mle[5]
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4]    ;se_xi <- gev.test$se[5]
  mu_int <- official_int[1:3] ;sig_int <- official_int[4]  ;xi_int <- official_int[5]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*ti[,1]+mu[3]*ti[,2]),sig,xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv) 
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV300 <- function(a, official_int=NULL){
  # set model name
  model='GEV300'
  # set initial value
  
  ti <- matrix(nrow=length(a),ncol=1);   ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 100; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  # estimate parameter
  try({gev.test <- optim(par=official_int[1:5], fn=gev.llh, x=a, t=ti, model=model, method="BFGS", hessian=T)}, silent=T)
  try(gev.test$se <- sqrt(diag(solve(gev.test$hessian))), silent = T)
  nllh <- gev.test$value      ;n <- length(a)         ;k <- length(gev.test$par) 
  mu <- gev.test$par[1:3]     ;sig <- gev.test$par[4]      ;xi <- gev.test$par[5]
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4]    ;se_xi <- gev.test$se[5]
  mu_int <- official_int[1:3] ;sig_int <- official_int[4]  ;xi_int <- official_int[5]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1)) 
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*exp(-mu[3]*ti[,1])),sig,xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV010 <- function(a, official_int=NULL){
  # set model name
  model='GEV010'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1)   ;  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({ gev.test <- gev.fit(xdat=a, ydat=ti, sigl=1, siglink=exp, muinit=official_int[1], siginit=official_int[2:3], shinit=official_int[4], show=FALSE)
  }, silent=T)
  if(is.null(gev.test)){
    gev.test <- gev.fit(xdat=a, ydat=ti, sigl=1, siglink=exp, show=FALSE)
    official_int[1:4] <- NA
  }
  
  nllh <- gev.test$nllh     ;n <- length(a)           ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1]     ;sig <- gev.test$mle[2:3]      ;xi <- gev.test$mle[4]
  se_mu <- gev.test$se[1]   ;se_sig <- gev.test$se[2:3]    ;se_xi <- gev.test$se[4]
  mu_int <- official_int[1] ;sig_int <- official_int[2:3]  ;xi_int <- official_int[4]
  BIC <- 2*nllh+log(n)*k    ;AIC <- 2*nllh+2*k             ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu,sig=exp(sig[1]+sig[2]*ti[,1]),xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV110 <- function(a, official_int=NULL){
  #set model name
  model='GEV110'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1)   ;  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gev.fit(xdat=a, ydat=ti, mul=1, sigl=1, siglink=exp, muinit=official_int[1:2], siginit=official_int[3:4], shinit=official_int[5], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gev.test <- gev.fit(xdat=a, ydat=ti, mul=1, sigl=1, siglink=exp, show=FALSE)
    official_int[1:5] <- NA
  }
  
  nllh <- gev.test$nllh       ;n <- length(a)                 ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:2]     ;sig <- gev.test$mle[3:4]     ;xi <- gev.test$mle[5]
  se_mu <- gev.test$se[1:2]   ;se_sig <- gev.test$se[3:4]   ;se_xi <- gev.test$se[5]
  mu_int <- official_int[1:2] ;sig_int <- official_int[3:4] ;xi_int <- official_int[5]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k            ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*ti[,1]),sig=exp(sig[1]+sig[2]*ti[,1]),xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV210 <- function(a, official_int=NULL){
  # set model name 
  model='GEV210'
  # set initial name
  ti <- matrix(nrow=length(a),ncol=2)   ;  ti[,1] <- seq(1:length(a))     ; ti[,2] <- ti[,1]^2
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){})}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gev.fit(xdat=a, ydat=ti, mul=c(1,2), sigl=1, siglink=exp, muinit=official_int[1:3], siginit=official_int[4:5], shinit=official_int[6], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gev.fit(xdat=a, ydat=ti, mul=c(1,2), sigl=1, siglink=exp, show=FALSE)
    official_int[1:6] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)                  ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:3]     ;sig <- gev.test$mle[4:5]      ;xi <- gev.test$mle[6]
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4:5]    ;se_xi <- gev.test$se[6]
  mu_int <- official_int[1:3] ;sig_int <- official_int[4:5]  ;xi_int <- official_int[6]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k             ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*ti[,1]+mu[3]*ti[,2]),sig=exp(sig[1]+sig[2]*ti[,1]),xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GEV310 <- function(a, official_int=NULL){
  # set model name
  model='GEV310'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1);   ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 100; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  # estimate parameter
  try({gev.test <- optim(par = official_int[1:6], fn=gev.llh, x=a, t=ti, model=model, method = "BFGS", hessian = T)}, silent = T)
  try(gev.test$se <- sqrt(diag(solve(gev.test$hessian))),silent = T)
  nllh <- gev.test$value      ;n <- length(a)                ;k <- length(gev.test$par) 
  mu <- gev.test$par[1:3]     ;sig <- gev.test$par[4:5]    ;xi <- gev.test$par[6]
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4:5]  ;se_xi <- gev.test$se[6]
  mu_int <- official_int[1:3] ;sig_int <- official_int[4:5];xi_int <- official_int[6]
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GEV(z_t,mu=(mu[1]+mu[2]*exp(-mu[3]*ti[,1])),sig=exp(sig[1]+sig[2]*ti[,1]),xi)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD000 <- function(a, official_int=NULL){
  # set model name
  model='GD000'
  # set initial value
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1;official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T }, error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  try({gev.test <- gum.fit(xdat=a, muinit=official_int[1], siginit=official_int[2], show=FALSE)}, silent=T)
  nllh <- gev.test$nllh     ;n <- length(a)         ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1]     ;sig <- gev.test$mle[2]      ;xi <- 0
  se_mu <- gev.test$se[1]   ;se_sig <- gev.test$se[2]    ;se_xi <- NA
  mu_int <- official_int[1] ;sig_int <- official_int[2]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k    ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  KS <- ks.test(z_t,"pgev",loc=mu,scale=sig,shape=xi)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv) 
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD100 <- function(a, official_int=NULL){
  # set model name
  model='GD100'
  #set initial value
  ti <- matrix(nrow=length(a),ncol=1)   ;  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){}) }
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gum.fit(xdat=a, ydat=ti, mul=1, muinit=official_int[1:2], siginit=official_int[3], show=FALSE)}, silent=T)
  
  if(is.null(gev.test)){
    gev.test <- gum.fit(xdat=a, ydat=ti, mul=1, show=FALSE)
    official_int[1:3] <- NA
  }
  
  nllh <- gev.test$nllh       ;n <- length(a)         ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:2]     ;sig <- gev.test$mle[3]      ;xi <- 0
  se_mu <- gev.test$se[1:2]   ;se_sig <- gev.test$se[3]    ;se_xi <- NA
  mu_int <- official_int[1:2] ;sig_int <- official_int[3]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*ti[,1]),sig)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD200 <- function(a, official_int=NULL){
  # set model name
  model='GD200'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=2)   ;  ti[,1] <- seq(1:length(a))     ; ti[,2] <- ti[,1]^2
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gum.fit(xdat=a, ydat=ti, mul=c(1,2), muinit=official_int[1:3], siginit=official_int[4], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gum.fit(xdat=a, ydat=ti, mul=c(1,2), show=FALSE)
    official_int[1:4] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)         ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:3]     ;sig <- gev.test$mle[4]      ;xi <- 0
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4]    ;se_xi <- NA
  mu_int <- official_int[1:3] ;sig_int <- official_int[4]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*ti[,1]+mu[3]*ti[,2]),sig)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD300 <- function(a, official_int=NULL){
  # set model name
  model='GD300'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1);   ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 100; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  # estimate parameter
  try({gev.test <- optim(par = official_int[1:4], fn=gev.llh, x=a, t=ti, model=model, method="BFGS", hessian=T)}, silent=T)
  try(gev.test$se <- sqrt(diag(solve(gev.test$hessian))),silent = T)
  nllh <- gev.test$value      ;n <- length(a)         ;k <- length(gev.test$par) 
  mu <- gev.test$par[1:3]     ;sig <- gev.test$par[4]      ;xi <- 0
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4]    ;se_xi <- NA
  mu_int <- official_int[1:3] ;sig_int <- official_int[4]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*exp(-mu[3]*ti[,1])),sig)
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD010 <- function(a, official_int=NULL){
  # set model name
  model='GD010'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1)   ;  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gum.fit(xdat=a, ydat=ti, sigl=1, siglink=exp, muinit=official_int[1], siginit=official_int[2:3], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gum.fit(xdat=a, ydat=ti, sigl=1, siglink=exp, show=FALSE)
    official_int[1:3] <- NA
  }
  nllh <- gev.test$nllh     ;n <- length(a)           ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1]     ;sig <- gev.test$mle[2:3]      ;xi <- 0
  se_mu <- gev.test$se[1]   ;se_sig <- gev.test$se[2:3]    ;se_xi <- NA
  mu_int <- official_int[1] ;sig_int <- official_int[2:3]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k    ;AIC <- 2*nllh+2*k             ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu,sig=exp(sig[1]+sig[2]*ti[,1]))
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv) 
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD110 <- function(a, official_int=NULL){
  # set model name
  model='GD110'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1)   ;  ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gum.fit(xdat=a, ydat=ti, mul=1, sigl=1, siglink=exp, muinit=official_int[1:2], siginit=official_int[3:4], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gum.fit(xdat=a, ydat=ti, mul=1, sigl=1, show=FALSE)
    official_int[1:4] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)          ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:2]     ;sig <- gev.test$mle[3:4]     ;xi <- 0
  se_mu <- gev.test$se[1:2]   ;se_sig <- gev.test$se[3:4]   ;se_xi <- NA
  mu_int <- official_int[1:2] ;sig_int <- official_int[3:4] ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k            ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*ti[,1]),sig=exp(sig[1]+sig[2]*ti[,1]))
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD210 <- function(a, official_int=NULL){
  # set model name
  model='GD210'
  ti <- matrix(nrow=length(a),ncol=2)   ;  ti[,1] <- seq(1:length(a))     ; ti[,2] <- ti[,1]^2
  
  if(is.null(official_int)){
    maxi.it <- 50; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  gev.test <- NULL
  try({gev.test <- gum.fit(xdat=a, ydat=ti, mul=c(1,2), sigl=1, siglink=exp, muinit=official_int[1:3], siginit=official_int[4:5], show=FALSE)}, silent=T)
  if(is.null(gev.test)){
    gev.test <- gum.fit(xdat=a, ydat=ti, mul=c(1,2), sigl=1, siglink=exp, show=FALSE)
    official_int[1:5] <- NA
  }
  nllh <- gev.test$nllh       ;n <- length(a)                  ;k <- length(gev.test$mle) 
  mu <- gev.test$mle[1:3]     ;sig <- gev.test$mle[4:5]      ;xi <- 0
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4:5]    ;se_xi <- NA
  mu_int <- official_int[1:3] ;sig_int <- official_int[4:5]  ;xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k             ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*ti[,1]+mu[3]*ti[,2]),sig=exp(sig[1]+sig[2]*ti[,1]))
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

GD310 <- function(a, official_int=NULL){
  # set model name
  model='GD310'
  # set initial value
  ti <- matrix(nrow=length(a),ncol=1);   ti[,1] <- seq(1:length(a))
  
  if(is.null(official_int)){
    maxi.it <- 100; Fal <- F; maxitCnt <- 0
    while (Fal==F && maxitCnt < maxi.it) {
      tryCatch({maxitCnt=maxitCnt+1; official_int <- initial_random_func(a=a, t=ti, model=model)[-1]
      Fal <- T },error=function(e){} )}
    if(maxitCnt>=maxi.it){ return(NULL) }
  }
  
  # estimate parameter
  try({gev.test <- optim(par = official_int[1:5], fn=gev.llh, x=a, t=ti, model=model, method="BFGS", hessian=T)}, silent=T)
  try(gev.test$se <- sqrt(diag(solve(gev.test$hessian))),silent = T)
  nllh <- gev.test$value      ;n <- length(a)                ;k <- length(gev.test$par) 
  mu <- gev.test$par[1:3]     ;sig <- gev.test$par[4:5]    ;xi <- 0
  se_mu <- gev.test$se[1:3]   ;se_sig <- gev.test$se[4:5]  ;se_xi <- NA
  mu_int <- official_int[1:3] ;sig_int <- official_int[4:5];xi_int <- 0
  BIC <- 2*nllh+log(n)*k      ;AIC <- 2*nllh+2*k           ;AICc <- 2*nllh+2*k+((2*k*(k+1))/(n-k-1))  
  
  z_t <- rm_tie(a)
  z_t_tran <- trans.GUM(z_t,mu=(mu[1]+mu[2]*exp(-mu[3]*ti[,1])),sig=exp(sig[1]+sig[2]*ti[,1]))
  KS <- ks.test(rm_tie(z_t_tran),"pgev",loc=0,scale=1,shape=0)
  KS_ST <- KS$statistic
  KS_pv <- KS$p.value
  
  result <- c(model,fillNA(mu,sig,xi),fillNA(se_mu,se_sig,se_xi),fillNA(mu_int,sig_int,xi_int),
              nllh,BIC,AIC,AICc,KS_ST,KS_pv)  
  
  names(result) <- c("Model",paste0("mu",c(0,1,2)),paste0("sig",c(0,1)),paste0("xi",c(0,1)),
                     paste0("se_mu",c(0,1,2)),paste0("se_sig",c(0,1)),paste0("se_xi",c(0,1)),
                     paste0("mu_int",c(0,1,2)),paste0("sig_int",c(0,1)),paste0("xi_int",c(0,1)),
                     "Nllh","BIC","AIC","AICc","KS_ST","KS_pv")
  
  return(result)
}

fit_all <- function(a){
  table <- NULL
  try(table <- rbind(table,GEV000(a)),silent = T)
  try(table <- rbind(table,GEV100(a)),silent = T)
  try(table <- rbind(table,GEV200(a)),silent = T)
  try(table <- rbind(table,GEV300(a)),silent = T)
  try(table <- rbind(table,GEV010(a)),silent = T)
  try(table <- rbind(table,GEV110(a)),silent = T)
  try(table <- rbind(table,GEV210(a)),silent = T)
  try(table <- rbind(table,GEV310(a)),silent = T)
  try(table <- rbind(table,GD000(a)),silent = T)
  try(table <- rbind(table,GD100(a)),silent = T)
  try(table <- rbind(table,GD200(a)),silent = T)
  try(table <- rbind(table,GD300(a)),silent = T)
  try(table <- rbind(table,GD010(a)),silent = T)
  try(table <- rbind(table,GD110(a)),silent = T)
  try(table <- rbind(table,GD210(a)),silent = T)
  try(table <- rbind(table,GD310(a)),silent = T)
  table <- data.frame(table,stringsAsFactors=F)
  for(i in 2:ncol(table)){
    table[,i] <- as.numeric(table[,i])
  }
  return(table)
}

best_ens_model <- function(a.fit, n.mod, cri){
  #rank_model <- a.fit[order(a.fit[,cri], decreasing = FALSE), ]
  rank_model <- a.fit
  bm_ens <- rank_model[1:n.mod,] 
  bm_ens_split <- split(bm_ens, as.factor(bm_ens$Model))
  bm_ens_RL <- lapply(bm_ens_split, rp=100, nt=100, get_rtlv) %>% do.call(cbind, .)
  
  rm_mod <- NULL
  upper_rl <- 10^4
  lower_rl <- 0
  if(any(bm_ens_RL>=upper_rl) || any(bm_ens_RL<=lower_rl)){
    rm_mod <- colnames(bm_ens_RL)[c(which(bm_ens_RL>=upper_rl), which(bm_ens_RL<=lower_rl))]
  }else{
    return(bm_ens)
  }
  if(length(rm_mod)>0){
    re_rank_model <- rank_model
    save_idx <- NULL
    for (l.rm_mod in 1:length(rm_mod)) { 
      save_idx <- c(save_idx,which(re_rank_model$Model==rm_mod[l.rm_mod]) )
    }
    re_rank_model <- re_rank_model[-save_idx,]
    
    if(n.mod>nrow(re_rank_model)){
      re_choose_ens <- best_ens_model(a.fit=re_rank_model,n.mod=nrow(re_rank_model), cri=cri)  
    }else{
      re_choose_ens <- best_ens_model(a.fit=re_rank_model,n.mod=n.mod,cri=cri)  
    }
  }
  return(re_choose_ens)
} 
