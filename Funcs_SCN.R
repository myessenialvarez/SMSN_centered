#################################################################################################
## Functions related to Skew Contaminated Normal (SCN) Distribution
## 
## Last version date: 19/12/2023
## Author: Maria Yessenia Alvarez Gil 
##
#################################################################################################

# Funcion de densidad
DenSCN <-  function(y, mu, sigma2, shape, nu){
  dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
  return(dens)
}

# Funcion de distribucion Acumulada
AcumSCN <- function(y,mu,sigma2,shape,nu){
  n      <- length(y)
  resp   <- matrix(0, n, 1)
  delta  <- shape/sqrt(1 + shape^2)
  SIGMA1 <- matrix(c(sigma2/(nu[2]), -delta * sqrt(sigma2/nu[2]), 
                     -delta * sqrt(sigma2/nu[2]), 1), byrow = TRUE, ncol = 2, nrow = 2)
  SIGMA2 <- matrix(c(sigma2, -delta * sqrt(sigma2), 
                     -delta * sqrt(sigma2), 1), byrow = TRUE, ncol = 2, nrow = 2)
  if (length(mu) == 1) {
    MU <- cbind(rep(mu, n), 0)
  }
  if (length(mu) == n) {
    MU <- cbind(mu, 0)
  }
  Y <- cbind(y, 0)
  
  for (i in 1:n) {
    resp[i] <- 2 * ( nu[1]*mnormt::pmnorm(x = Y[i, ], mean = MU[i, ], varcov = SIGMA1) 
                     + (1-nu[1]) * mnormt::pmnorm(x = Y[i, ], mean = MU[i, ], varcov = SIGMA2)) 
  }
  return(resp)
}

# Diferencia entre acumulada LS y acumulada LI
pSNI_NC<- function(LI, LS, mu, sigma2, shape, nu, type = "SCN"){
  n      <- length(LI)
  resp   <- matrix(0, n, 1)
  delta  <- shape/sqrt(1 + shape^2)
  SIGMA1 <- matrix(c(sigma2/(nu[2]), -delta * sqrt(sigma2/nu[2]), 
                     -delta * sqrt(sigma2/nu[2]), 1), byrow = TRUE, ncol = 2, nrow = 2)
  SIGMA2 <- matrix(c(sigma2, -delta * sqrt(sigma2), 
                     -delta * sqrt(sigma2), 1), byrow = TRUE, ncol = 2, nrow = 2)
  if (length(mu) == 1) {
    MU <- cbind(rep(mu, n), 0)
  }
  if (length(mu) == n) {
    MU <- cbind(mu, 0)
  }
  YI <- cbind(LI, 0)
  YI[LI==-Inf,]<-c(-Inf,-Inf)
  YS <- cbind(LS, 0)
  #YS[LS==Inf,]<-c(Inf,Inf)
  
  for (i in 1:n) {
    resp[i] <- 2 * (( nu[1]*mnormt::pmnorm(x = YS[i, ], mean = MU[i, ], varcov = SIGMA1) 
                      + (1-nu[1]) * mnormt::pmnorm(x = YS[i, ], mean = MU[i, ], varcov = SIGMA2)) -
                      ( nu[1]*mnormt::pmnorm(x = YI[i, ], mean = MU[i, ], varcov = SIGMA1) 
                        + (1-nu[1]) * mnormt::pmnorm(x = YI[i, ], mean = MU[i, ], varcov = SIGMA2)))
  }
  if(length(which(resp <= 0)) > 0) resp[which(resp <= 0)] <- .Machine$double.xmin
  return(resp)
}


# Verosimilanca 
dSCNMod <- function(cc, y, mu, sigma2 = 1, shape, nu, cens="Left", LS=NULL){
  densN<- vector(mode = "numeric", length = length(y))
  densN[cc==0] <- DenSCN(y[cc==0], mu[cc==0], sigma2, shape, nu ) 
  if (cens=="Left"){
    densN[cc==1]<- AcumSCN(y[cc==1], mu[cc==1], sigma2, shape, nu)     
  }
  if (cens=="Right"){
    densN[cc==1]<- 1-AcumSCN(y[cc==1], mu[cc==1], sigma2, shape, nu)     
  }
  if (cens=="Interv"){
    densN[cc==1]<-pSNI_NC(y[cc==1], LS[cc==1], mu[cc==1], sigma2, shape, nu, type = "SCN")
  }
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}

