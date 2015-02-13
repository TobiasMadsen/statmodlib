#Diverse funktioner til brug for Statistiske modeller 2

#USS,SSD,S,s^2, fishers Z  -------------------------
#Udregningen af SSD
ssd <- function(x){
  var(x)*(length(x)-1)
}

#Udregning af USS
uss <- function(x){
  sum(x*x)
}

#' Fishers Z
#' 
#' Wrapper for atanh
fisherZ <- function(r){
  atanh(r)
}

#' Invers fishers Z
#' 
#' Wrapper for tanh
invFisherZ <- function(z){
  tanh(z)
}

#Udregning sum
#brug R's indbyggede sum

#Estimater i 2D normal----------------

#' Estimer 2D normal fordeling
#' 
#' (x,y) er en observationsrække fra den 2dimensionale 
#' normalfordeling. 
#' 
#' @param S_x summen af x'er
#' @param S_y summen af y'er
#' @param USS_x summen af x^2
#' @param USS_y summen af y^2
#' @param SP_xy summen af x*y
#' @param n stikprøvestørrelsen
estimates2d <- function(S_x,S_y,USS_x,USS_y,SP_xy,n){
  #Se side 25
  mu_x <- S_x/n
  mu_y <- S_y/n
  s2_x <- (USS_x - n*mu_x**2)/(n-1)
  s2_y <- (USS_y - n*mu_y**2)/(n-1)
  rho <- (SP_xy - S_x*S_y/n)/sqrt(s2_x*(n-1))/sqrt(s2_y*(n-1))

  return(c("mean x" = mu_x,
           "mean y" = mu_y,
           "var x"  = s2_x,
           "var y"  = s2_y,
           "rho"    = rho))
}

#' Estimer 2D normal fordeling
#' 
#' (x,y) er en observationsrække fra den 2dimensionale 
#' normalfordeling. 
#' 
#' @param S_x summen af x'er
#' @param S_y summen af y'er
#' @param SSD_x summen af (x-mean(x))^2
#' @param SSD_y summen af (y-mean(y))^2
#' @param SPD_xy summen af (x-mean(x))*(y-mean(y))
#' @param n stikprøvestørrelsen
estimates2dSSD <- function(S_x,S_y,SSD_x,SSD_y,SPD_xy,n){
  mu_x <- S_x/n
  mu_y <- S_y/n
  s2_x <- SSD_x/(n-1)
  s2_y <- SSD_y/(n-1)
  rho <- (SPD_xy)/sqrt(SSD_x)/sqrt(SSD_y)

  return(c("mean x" = mu_x,
           "mean y" = mu_y,
           "var x"  = s2_x,
           "var y"  = s2_y,
           "rho"    = rho))
}

#Konfidensintervaller---------------

#' Konfidensinterval for fishers Z
#' 
#' @param z estimat for fisherZ
#' @param n stikprøvestørrelse
confZ <- function(z,n){
  #Side 34
  z+c(-1,0,1)*qnorm(0.975)/sqrt(n-3)
}

#' Konfidensinterval for rho
#' 
#' @param r estimat for rho
#' @param n stikprøvestørrelse
confRho <- function(r,n){
  #Side 34
  res <- invFisherZ(confZ(fisherZ(r),n))
  names(res) = c("nedre", "estimat", "øvre")
  res
}

#' Konfidensinterval for fælles Z
#' 
#' @param z estimat for fishers Z
#' @param n stikprøvestørrelse
confZfaelles <- function(z,n){
  #Se side 36
  estZ <- sum(z*(n-3))/sum(n-3)
  estZ+c(-1,0,1)*qnorm(0.975)/sqrt(sum(n-3))
}

#' Konfidensinterval for fælles rho
#' 
#' @param r estimat for rho
#' @param n stikprøvestørrelse
confRhofaelles <- function(r,n){
  #Se side 36
  invFisherZ(confZfaelles(fisherZ(r),n))
}

#Test ----------------

#' Test lighed af korrelation k>=2 stikprøver
#'
#' @param z vektor med estimater for fishersZ for hver observationsrække
#' @param n vektor med stikprøvestørrelser
testLighedZ <- function(z,n){
  estZ <- sum(z*(n-3))/sum(n-3)
  X2 <- sum((z-estZ)**2*(n-3))
  p_obs <- 1-pchisq(X2,length(z)-1)
  c(estZ,X2,p_obs)
}

#' Test lighed af korrelation k>=2 stikprøver
#'
#' @param r vektor med estimator for rho for hver observationsrække
#' @param n vektor med stikprøvestørrelser
testLighedRho <- function(rho,n){
  res <- testLighedZ( fisherZ(rho),n )
  names(res) = c("estZ","X2","p_obs")
  res
}

#' Test korrelation = 0 i en stikprøve
#' 
#' @param rho estimat for rho
#' @param n stikprøvestørrelse
testRhoNul <- function(rho,n){
  Test <- sqrt(n-2)*(rho)/sqrt(1-rho**2)
  p_obs <- 2*(1-pt(abs(Test),n-2))
  c(Test,p_obs)
}

#' Hotellings T test
#' 
#' Hotellings T test for H0: (mu_x, mu_y) = (mux0, muy0)
#' 
#' @param mux0 mux0
#' @param muy0 muy0
#' @param muxhat estimat for mu_x
#' @param muyhat estimat for mu_y
#' @param S 2x2 matrix med estimat for fælles kovarians matrix (se formel (2.49))
hotellingsT <- function(mux0,muy0,muxhat,muyhat,S,n){
  S_inv <- solve(S)
  v <- c(muxhat - mux0, muyhat - muy0)
  Test <- n*t(v) %*% S_inv %*% v
  p_obs <- 1-pf(Test*(n-2)/(2*(n-1)),2,n-2)
  c(Test,p_obs)
}