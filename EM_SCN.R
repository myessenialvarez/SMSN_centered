#################################################################################################
## EM algorithm for estimating the parameters of SCN-CR Model 
## with Centered Error (Interval Censored)
##
## Version based on Victor Hugo Lachos' algorithm
## Last version date: 19/12/2023
## Author: Maria Yessenia Alvarez Gil  
#################################################################################################


EM_SNormC <- function(cc, x, y, beta = NULL, sigma2 = NULL, shape = NULL, nu = NULL, cens="Left", LS=NULL, get.init = TRUE, show.envelope = "FALSE", error = 0.000001, iter.max = 500, family = "SCN"){
 
  p          <- ncol(x)
  n          <- nrow(x)
  
  ## Valores iniciales
  y0        <- y[cc==0]
  x0        <- x[cc==0,]
  reg        <- lm(y0 ~ x0[ , 2:p])
  beta       <- as.vector(coefficients(reg), mode = "numeric")
  sigma2     <- sum((y0 - x0%*%beta)^2)/(n - p)
  shape      <- as.numeric(sign(skewness(y0 - x0%*%beta))*3)
  delta      <- shape / (sqrt(1 + shape^2))
  Delta      <- sqrt(sigma2)*delta
  Gama       <- sigma2 - Delta^2
  muaux      <- x%*%beta
  nu         <- c(0.4,0.7) 
  
  k1         <- (nu[1]/(nu[2]^(0.5))) + 1-nu[1]
  b          <- -sqrt(2/pi)*k1
  
  mu         <- muaux + b*Delta
  
  Lim1  <- Lim2 <- vector(mode = "numeric", length = length(y))
  
  if (cens=="Left")
  {
    Lim1 <- rep(-Inf,n)
    Lim2 <- y
  }
  if (cens=="Right")
  {
    Lim1 <- y
    Lim2 <- rep(Inf,n)
  }
  if (cens=="Interv")
  {
    Lim1 <- y
    Lim2 <- LS
  }
  
  cont       <- 0
  criterio   <- 1
  lkante     <- 1
  ychap      <- y
  n.par      <- p + 2
  
  y1         <- y[cc==1]
  Lim1aux    <- Lim1[cc==1]
  Lim2aux    <- Lim2[cc==1]
  
  while(criterio > error){
    
    cont     <- (cont + 1)
    print(cont)
    print(c(beta,sigma2,shape,nu))
    
    Mtij2    <- 1/(1 + (Delta^2)*(Gama^(-1)))
    Mtij     <- sqrt(Mtij2)
    mutij    <- Mtij2*Delta*(Gama^(-1))*(y - mu)
    A        <- mutij / Mtij
    cnu1     <- (nu[1]*(nu[2])^(0.5))/sqrt(2*pi*(1+shape^2))
    cnu2     <- (1-nu[1])/sqrt(2*pi*(1+shape^2))
    
    u        <- 2/(DenSCN(y, mu, sigma2, shape, nu))*(nu[1]*nu[2]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-0.5)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-0.5)*(y-mu))) 
    E        <- 2/(DenSCN(y, mu, sigma2, shape, nu))*(nu[1]*(nu[2])^(0.5)*dnorm(y, mu, sqrt(sigma2/nu[2]))*dnorm(sqrt(nu[2])*shape*sigma2^(-0.5)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*dnorm(shape*sigma2^(-0.5)*(y-mu)))
    
    S1       <- u
    S2       <- ((mutij + b)*u + Mtij*E)
    S3       <- ((mutij + b)^2*u + Mtij2 + Mtij*(mutij+2*b)*E)
    
    ## Sin censura  
    E00      <- S1               #E[U|Y]
    E01      <- y*S1             #E[UY|Y]
    E02      <- y^2*S1           #E[UY^2|Y]
    E10      <- S2               #E[UT|Y]
    E20      <- S3               #E[UT^2|Y]
    E11      <- y*S2             #E[UTY|Y]
    
    AuxDen    <- 1/(1+shape^2)
    sigma2s   <- sigma2/nu[2]
    sigma2ss  <- (sigma2/nu[2])*AuxDen
    sigma2sss <- sigma2*AuxDen
    
    Aux1     <- pSNI(Lim1,    Lim2, mu,   sigma2s,   shape, NULL, type = "SN")
    Aux11    <- pSNI(Lim1,    Lim2, mu,    sigma2,   shape, NULL, type = "SN")
    Aux2     <- pSNI_NC(Lim1, Lim2, mu,    sigma2,   shape,   nu, type = "SCN")
    Aux3     <- pSNI(Lim1,    Lim2, mu,  sigma2ss,       0, NULL, type = "SN")
    Aux31    <- pSNI(Lim1,    Lim2, mu, sigma2sss,       0, NULL, type = "SN")
    
  
    mu1       <- mu[cc==1]
    np        <- length(mu1)
    aux1MomW1 <- aux2MomW2 <- aux3MomW3 <- aux4MomW4 <- matrix(0, np, 2)
    
    for(j in 1:np){
      A1a           <- MomenSNI2023(mu1[j],    sigma2, shape, NULL, delta=NULL, Lim1=Lim1aux[j], Lim2=Lim2aux[j], type = "SN") 
      A2a           <- MomenSNI2023(mu1[j],   sigma2s, shape, NULL, delta=NULL, Lim1=Lim1aux[j], Lim2=Lim2aux[j], type = "SN")
      A3a           <- MomenSNI2023(mu1[j],  sigma2ss,     0, NULL, delta=NULL, Lim1=Lim1aux[j], Lim2=Lim2aux[j], type = "SN") 
      A4a           <- MomenSNI2023(mu1[j], sigma2sss,     0, NULL, delta=NULL, Lim1=Lim1aux[j], Lim2=Lim2aux[j], type = "SN")
      
      aux1MomW1[j,] <- c(A1a$EUY1, A1a$EUY2)
      aux2MomW2[j,] <- c(A2a$EUY1, A2a$EUY2)
      aux3MomW3[j,] <- c(A3a$EUY1, A3a$EUY2)
      aux4MomW4[j,] <- c(A4a$EUY1, A4a$EUY2)
    }
    
    P1aux <- P2aux  <- P3aux <- P4aux <- P5aux <- P6aux <- WW <- u
    P1aux[cc==1]    <- aux1MomW1[,1]
    P2aux[cc==1]    <- aux1MomW1[,2]
    P3aux[cc==1]    <- aux2MomW2[,1]
    P4aux[cc==1]    <- aux2MomW2[,2]
    P5aux[cc==1]    <- aux3MomW3[,1]
    P6aux[cc==1]    <- aux4MomW4[,1]
   
    AuxT0          <- (2/Aux2)*(cnu1*Aux3 + cnu2*Aux31)
    AuxT1          <- (2/Aux2)*(cnu1*Aux3*P5aux + cnu2*Aux31*P6aux)
    
    # Com censura 
    E00Aux         <- (1/Aux2)*(nu[1]*nu[2]*Aux1+(1-nu[1])*Aux11) 
    E01Aux         <- (1/Aux2)*(nu[1]*nu[2]*Aux1*P3aux+(1-nu[1])*Aux11*P1aux)
    E02Aux         <- (1/Aux2)*(nu[1]*nu[2]*Aux1*P4aux+(1-nu[1])*Aux11*P2aux)
    E10Aux         <- Delta/(Gama + Delta^2)*(E01Aux - E00Aux*mu) + b*E00Aux + Mtij*AuxT0
    E20Aux         <- (Delta/(Gama + Delta^2))^2*(E02Aux - 2*E01Aux*mu + mu^2*E00Aux) + 2*b*(Delta/(Gama + Delta^2))*(E01Aux - E00Aux*mu) + b^2*E00Aux + Mtij2 + Mtij*((Delta/(Gama + Delta^2))*(AuxT1-mu*AuxT0) + 2*b*AuxT0 )
    E11Aux         <- Delta/(Gama + Delta^2)*(E02Aux - E01Aux*mu) + b*E01Aux + Mtij*AuxT1
    
    E00[cc==1]     <- E00Aux[cc==1]        # E[U]
    E01[cc==1]     <- E01Aux[cc==1]        # E[UY]  
    E02[cc==1]     <- E02Aux[cc==1]        # E[UY2]   
    E10[cc==1]     <- E10Aux[cc==1]        # E[UT]
    E20[cc==1]     <- E20Aux[cc==1]        # E[UT2]
    E11[cc==1]     <- E11Aux[cc==1]        # E[UTY]
    
    beta           <- solve(t(x)%*%diag(c(E00))%*%x)%*%(t(x)%*%matrix(E01 - E10*Delta, n, 1))
    muaux          <- x%*%beta
    Delta          <- sum(E11 - E10*muaux)/sum(E20)
    Gama           <- sum(E02 - 2*E01*muaux + E00*muaux^2 + Delta^2*E20 - 2*Delta*E11 + 2*Delta*E10*muaux)/n		
    sigma2         <- Gama + Delta^2
    shape          <- ((sigma2^(-1/2))*Delta)/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
    
    k1             <- (nu[1]/(nu[2]^(0.5))) + 1-nu[1]
    b              <- -sqrt(2/pi)*k1
    
    mu             <- muaux + b*Delta
    
    f<-function(nu){
      # a1 <- exp(nu[1])/(1+exp(nu[1]))
      # a2 <- exp(nu[2])/(1+exp(nu[2]))
      # if (a2 == 0) {
      #   a2 <- .Machine$double.xmin
      # }
      # nuaux = c(a1,a2)
      # ver <- sum(log(dSCNMod(cc, y, mu, sigma2 , shape, nuaux, cens=cens, LS=Lim2)))
 
      ver <- sum(log(dSCNMod(cc, y, mu, sigma2 , shape, nu, cens=cens, LS=Lim2)))
      return(-ver)
    }
    
    
    
    if (nu[2] == 0) {
      nu[2] <- .Machine$double.xmin
    }
    
    
      # nu <- optim(c(0.1, 0.1), f, method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(0.99, 0.99))$par
       
      # Art <- optim(c(0.1,0.1),f, method=c("BFGS"),control=list(maxit=20000))$par


     #   nu3 <- min(round(exp(Art[1])/(1+exp(Art[1]))+0.05,2),0.9)
     #   nu4 <- min(round(exp(Art[2])/(1+exp(Art[2]))+0.05,2),0.9)
     # 
     #  
      # nu <- c(nu3,nu4)
      
       nu <- nu
    
    lk             <- -f(nu)
    logver         <- lk
    criterio       <- abs((lk/lkante - 1))
    lkante         <- lk
    
    if (cont==iter.max){
      criterio     <- error/10
    }
  }
  
  teta_novo        <- matrix(c(beta, sigma2, shape, nu), ncol = 1) # to compute the number of parameters
  
    
  #########################
  ##  Information Matrix
  #########################
  
  sbeta            <- c()
  ssigma2          <- c()
  sshape           <- c()
  MIE              <- matrix(0, p+2, p+2)
  S                <- matrix(0,   1, p+2)
  sigma            <- sqrt(sigma2)
  shape            <- shape
  
  for(i in 1:n){
    sbeta          <- ((1 + shape^2)/sigma2)*(E01[i]*t(x[i,]) - E00[i]*t(x[i,])*mu[i] - Delta*E10[i]*t(x[i,]))
    ssigma2        <- -1/(2*sigma2) + ((1 + shape^2)/(2*sigma2^2))*(E02[i] - 2*E01[i]*mu[i] + (t(mu[i])%*%mu[i])*E00[i]) - ((shape*sqrt(1+shape^2))/(2*sigma^3))*(E11[i] - E10[i]*mu[i])
    sshape         <- shape/(1 + shape^2) - (shape/sigma2)*(E02[i] - 2*E01[i]*mu[i] + E00[i]*(t(mu[i])%*%mu[i])) + ((1 + 2*shape^2)/(sigma*sqrt(1+shape^2)))*(E11[i] - E10[i]*mu[i]) - shape*E20[i]      
    S              <- c(sbeta, ssigma2, sshape)
    MIE1           <- S%*%t(S)
    ind            <- lower.tri(MIE1)
    MIE1[ind]      <- t(MIE1)[ind]
    MIE            <- MIE1 + MIE
  }
  
  se               <- sqrt(diag(solve(MIE)))    # standard errors
  
  if(show.envelope=="TRUE"){
    envelop        <- EnvelopeRMT(teta_novo, y, x, cc, family = "ST")
  }
  
  aic              <- -2*lk +             2*length(teta_novo)
  bic              <- -2*lk +        log(n)*length(teta_novo)
  caic             <- -2*lk +  (log(n) + 1)*length(teta_novo)
  hq               <- -2*lk + 2*log(log(n))*length(teta_novo)
  
  #####################################################
  ### Summary of the model results.
  #####################################################
  
  namesx <- ('(Intercept)     ')
  
  if(ncol(as.matrix(x)) > 1){
    for(i in 2:ncol(as.matrix(x))){
      namesx <- cbind(namesx, paste("x", i-1,"     ", sep = ""))
    }
  }
  
  t.est            <- teta_novo[1:n.par]/se[1:n.par]
  p.est            <- 2*(1 - pt(abs(t.est), n-1))
  
  param            <- round(cbind(teta_novo[1:n.par], se[1:n.par], t.est, p.est), digits = 3)
  namespar         <- colnames(x)
  colx             <- ncol(as.matrix(x))
  
  if(length(namespar) == 0){
    namespar       <- namesx[1:colx]
  }
  
  if(family=="N" || family=="T" ){
    dimnames(param)  <- list(c(namespar, expression(sigma^2)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  if(family=="SN" || family=="ST" || family=="SCN"){
    dimnames(param)  <- list(c(namespar, expression(sigma^2), expression(shape)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  }
  
  cat('\n')
  cat('-------------------------------------------------------\n')
  cat('                  EM estimates and SE                  \n')
  cat('                                                       \n')
  cat('Coefficients:                                          \n')
  print(param)
  cat('-------------------------------------------------------\n')
  
  if (family == "CN" || family == "SCN") {
    cat("nu1 ", round(nu[1], 3), " and nu2", round(nu[2], 3), "\n")
  }
  
  cat('-------------------------------------------------------\n')
  cat('\r \n')
  critFin           <- c(logver, aic, bic, caic, hq)
  critFin           <- round(t(as.matrix(critFin)), digits = 3)
  dimnames(critFin) <- list(c("Value"), c("Loglik", "AIC", "BIC","CAIC", "HQ"))
  cat('\n')
  cat('Model selection criteria                               \n')
  cat('-------------------------------------------------------\n')
  print(critFin)
  cat('-------------------------------------------------------\n')
  cat('\r \n')
  
  return(list(theta=teta_novo, Se=se, nu=nu, iter=cont, logver=logver, AIC=aic,
              BIC=bic, CAIC=caic, HQ=hq, yhat=ychap, E00=E00, E01=E01, E02=E02,
              E10=E10, E20=E20, E11=E11))	
  
  } 
  
  
  
  
  
  
  
  

