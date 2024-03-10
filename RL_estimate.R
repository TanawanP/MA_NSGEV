
qgev <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE){
  if (!lower.tail) 
    p <- 1 - p
  if (shape == 0) 
    return(loc - scale * log(-log(p)))
  else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}

# Rtlv witout SE
#GEV000
rtlv.GEV000 <- function(x,rp){
  rl <- qgev(1-1/rp,loc=x$mu0,scale=x$sig0,shape=x$xi0)
  return(rl)
}

#GEV100
rtlv.GEV100 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  rl <- qgev(1-1/rp,mu,x$sig0,x$xi0)
  return(rl)
}

#GEV200
rtlv.GEV200 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*(nt^2)
  rl <- qgev(1-1/rp,mu,x$sig0,x$xi0)
  return(rl)
}

#GEV300
rtlv.GEV300 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  rl <- qgev(1-1/rp,mu,x$sig0,x$xi0)
  return(rl)
}

#GEV010
rtlv.GEV010 <- function(x,rp,nt){
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,x$mu0,sig,x$xi0)
  return(rl)
}

#GEV110
rtlv.GEV110 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,x$xi0)
  return(rl)
}

#GEV210
rtlv.GEV210 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*(nt^2)
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,x$xi0)
  return(rl)
}

#GEV310
rtlv.GEV310 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,x$xi0)
  return(rl)
}

#GEV001
rtlv.GEV001 <- function(x,rp,nt){
  xi <- x$xi0+x$xi1*nt
  rl <- x$mu0-(x$sig0/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV101
rtlv.GEV101 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  xi <- x$xi0+x$xi1*nt
  rl <- mu-(x$sig0/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV201
rtlv.GEV201 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*(nt^2)
  rl <- mu-(x$sig0/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV301
rtlv.GEV301 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  xi <- x$xi0+x$xi1*nt
  rl <- mu-(x$sig0/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV011
rtlv.GEV011 <- function(x,rp,nt){
  sig <- exp(x$sig0+x$sig1*nt)
  xi <- x$xi0+x$xi1*nt
  rl <- x$mu0-(sig/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV111
rtlv.GEV111 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  sig <- exp(x$sig0+x$sig1*nt)
  xi <- x$xi0+x$xi1*nt
  rl <- mu-(sig/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV211
rtlv.GEV211 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*(nt^2)
  sig <- exp(x$sig0+x$sig1*nt)
  xi <- x$xi0+x$xi1*nt
  rl <- mu-(sig/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GEV311
rtlv.GEV311 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  sig <- exp(x$sig0+x$sig1*nt)
  xi <- x$xi0+x$xi1*nt
  rl <- mu-(sig/xi)*(1-(-log(1-1/rp))^(-xi))
  return(rl)
}

#GD000
rtlv.GD000 <- function(x,rp){
  rl <- qgev(1-1/rp,x$mu0,x$sig0,0)
  return(rl)
}

#GD100
rtlv.GD100 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  rl <- qgev(1-1/rp,mu,x$sig0,0)
  return(rl)
}

#GD200
rtlv.GD200 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*nt^2
  rl <- qgev(1-1/rp,mu,x$sig0,0)
  return(rl)
}

#GD300
rtlv.GD300 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  rl <- qgev(1-1/rp,mu,x$sig0,0)
  return(rl)
}

#GD010
rtlv.GD010 <- function(x,rp,nt){
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,x$mu0,sig,0)
  return(rl)
}

#GD110
rtlv.GD110 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,0)
  return(rl)
}

#GD210
rtlv.GD210 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*nt+x$mu2*nt^2
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,0)
  return(rl)
}

#GD310
rtlv.GD310 <- function(x,rp,nt){
  mu <- x$mu0+x$mu1*exp(-x$mu2*nt)
  sig <- exp(x$sig0+x$sig1*nt)
  rl <- qgev(1-1/rp,mu,sig,0)
  return(rl)
}

get_rtlv <- function(x,rp,nt){
  if(x$Model=="GEV000"){
    return(rtlv.GEV000(x,rp))
  }else  if(x$Model=="GEV100"){
    return(rtlv.GEV100(x,rp,nt))
  }else if(x$Model=="GEV200"){
    return(rtlv.GEV200(x,rp,nt))
  }else if(x$Model=="GEV300"){
    return(rtlv.GEV300(x,rp,nt))
  }else if(x$Model=="GEV010"){
    return(rtlv.GEV010(x,rp,nt))
  }else if(x$Model=="GEV110"){
    return(rtlv.GEV110(x,rp,nt))
  }else if(x$Model=="GEV210"){
    return(rtlv.GEV210(x,rp,nt))
  }else if(x$Model=="GEV310"){
    return(rtlv.GEV310(x,rp,nt))
  }else if(x$Model=="GD100"){
    return(rtlv.GD100(x,rp,nt))
  }else if(x$Model=="GD000"){
    return(rtlv.GD000(x,rp))
  }else if(x$Model=="GD200"){
    return(rtlv.GD200(x,rp,nt))
  }else if(x$Model=="GD300"){
    return(rtlv.GD300(x,rp,nt))
  }else if(x$Model=="GD010"){
    return(rtlv.GD010(x,rp,nt))
  }else if(x$Model=="GD110"){
    return(rtlv.GD110(x,rp,nt))
  }else if(x$Model=="GD210"){
    return(rtlv.GD210(x,rp,nt))
  }else if(x$Model=="GD310"){
    return(rtlv.GD310(x,rp,nt))
  }
}