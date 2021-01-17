library(invgamma)
library(ggplot2)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)


setwd("C:/Users/franc/OneDrive - Politecnico di Milano/POLIMI/Bayesian statistics/TA_material/TA11")
Y <- read.csv("data.csv")
Y <- Y[complete.cases(Y),]
head(Y)
dim(Y)

# Compute M, number of Genus, and the vector N, 
# with the numerosities in each Genus:

m <-  length(unique(Y[,3])) 
num_genus <- as.vector(table(Y[,3]))
t <- length(Y)
Y = log((Y[,7]/Y[,5]) / (1 - Y[,7]/Y[,5]))

a0 = 2
b0 = 1
m0 = -5
s20 = 10
alpha0 = 2
beta0 = 1


# MANCA LA PROPOSAL, USO LA PRIOR

ABC_mean <- function(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold){
  
  thetas = rep(NA, m)
  thetas2 = matrix( nrow = max(num_genus), ncol = m)
  Ys = rep(NA, length(Y))
  sigma2 = rep(NA, niter)
  mu = rep(NA, niter)
  tau2 = rep(NA, niter)
  theta = matrix( nrow = niter, ncol = m)
  R = matrix( nrow = niter, ncol = length(Y))
  
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  for (i in 1:niter){
    
    sigma2s <- rinvgamma(1, a0, b0)
    mus <- rnorm(1, m0, s20)
    tau2s <- rinvgamma(1, alpha0, beta0)
    l=1
    
    #thetas2 <- data.frame(matrix(rnorm(m*max(num_genus), mus, tau2s), max(num_genus))) #max(num_genus) righe, m colonne
    
    #apply(thetas2,2,max,na.rm=1)
    
    
    for (j in 1:m){
      thetas[j] <- rnorm(1, mus, tau2s)
      for (k in 1:num_genus[j]){
        Ys[l] <- rnorm(1, thetas[j], sigma2s)
        l= l+1
      }
    }
    
    if (norm(mean(Ys)-mean(Y), type = "2") <= threshold){
      sigma2[i] = sigma2s
      mu[i] = mus
      tau2[i] = tau2s
      theta[i,] = thetas
      for (l in 1:t){
        R[l,] = Ys
      }
      
      nacp = nacp + 1
      
    }
    setTxtProgressBar(pb, i)
    
    
    
  }
  
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(list(sigma2 = sigma2, mu = mu, tau2 = tau2, theta = theta, R = R))
  
}

ABC_mean_kernel_gauss <- function(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold){
  
  thetas = rep(NA, m)
  thetas2 = matrix( nrow = max(num_genus), ncol = m)
  Ys = rep(NA, length(Y))
  sigma2 = rep(NA, niter)
  mu = rep(NA, niter)
  tau2 = rep(NA, niter)
  theta = matrix( nrow = niter, ncol = m)
  R = matrix( nrow = niter, ncol = length(Y))
  
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  for (i in 1:niter){
    
    sigma2s <- rinvgamma(1, a0, b0)
    mus <- rnorm(1, m0, s20)
    tau2s <- rinvgamma(1, alpha0, beta0)
    l=1
    
    #thetas2 <- data.frame(matrix(rnorm(m*max(num_genus), mus, tau2s), max(num_genus))) #max(num_genus) righe, m colonne
    
    #apply(thetas2,2,max,na.rm=1)
    
    
    for (j in 1:m){
      thetas[j] <- rnorm(1, mus, tau2s)
      for (k in 1:num_genus[j]){
        Ys[l] <- rnorm(1, thetas[j], sigma2s)
        l= l+1
      }
    }
    
    u <- runif(1)
    if (u < dnorm(norm(mean(Ys)-mean(Y), type="2")/threshold)/(threshold)){
      sigma2[i] = sigma2s
      mu[i] = mus
      tau2[i] = tau2s
      theta[i,] = thetas
      for (l in 1:t){
        R[l,] = Ys
      }
      nacp = nacp + 1
      
    }
    setTxtProgressBar(pb, i)
    
    
    
  }
  
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(list(sigma2 = sigma2, mu = mu, tau2 = tau2, theta = theta, R = R))
  
}

ABC_mean_kernel_unif <- function(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold){
  
  thetas = rep(NA, m)
  thetas2 = matrix( nrow = max(num_genus), ncol = m)
  Ys = rep(NA, length(Y))
  sigma2 = rep(NA, niter)
  mu = rep(NA, niter)
  tau2 = rep(NA, niter)
  theta = matrix( nrow = niter, ncol = m)
  R = matrix( nrow = niter, ncol = length(Y))
  
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  for (i in 1:niter){
    
    sigma2s <- rinvgamma(1, a0, b0)
    mus <- rnorm(1, m0, s20)
    tau2s <- rinvgamma(1, alpha0, beta0)
    l=1
    
    #thetas2 <- data.frame(matrix(rnorm(m*max(num_genus), mus, tau2s), max(num_genus))) #max(num_genus) righe, m colonne
    
    #apply(thetas2,2,max,na.rm=1)
    
    
    for (j in 1:m){
      thetas[j] <- rnorm(1, mus, tau2s)
      for (k in 1:num_genus[j]){
        Ys[l] <- rnorm(1, thetas[j], sigma2s)
        l= l+1
      }
    }
    #ifelse(abs(norm(mean(Ys)-mean(Y), type="2")/threshold)<=1,0.5,0)/threshold
    u <- runif(1)
    if (u < dunif(norm(mean(Ys)-mean(Y), type="2")/threshold)/(threshold)){
      sigma2[i] = sigma2s
      mu[i] = mus
      tau2[i] = tau2s
      theta[i,] = thetas
      for (l in 1:t){
        R[l,] = Ys
      }
      nacp = nacp + 1
      
    }
    setTxtProgressBar(pb, i)
    
    
    
  }
  
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(list(sigma2 = sigma2, mu = mu, tau2 = tau2, theta = theta, R = R))
  
}

ABC_median <- function(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold){
  
  thetas = rep(NA, m)
  thetas2 = matrix( nrow = max(num_genus), ncol = m)
  Ys = rep(NA, length(Y))
  sigma2 = rep(NA, niter)
  mu = rep(NA, niter)
  tau2 = rep(NA, niter)
  theta = matrix( nrow = niter, ncol = m)
  R = matrix( nrow = niter, ncol = length(Y))
  
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  for (i in 1:niter){
    
    sigma2s <- rinvgamma(1, a0, b0)
    mus <- rnorm(1, m0, s20)
    tau2s <- rinvgamma(1, alpha0, beta0)
    l=1
    
    #thetas2 <- data.frame(matrix(rnorm(m*max(num_genus), mus, tau2s), max(num_genus))) #max(num_genus) righe, m colonne
    
    #apply(thetas2,2,max,na.rm=1)
    
    
    for (j in 1:m){
      thetas[j] <- rnorm(1, mus, tau2s)
      for (k in 1:num_genus[j]){
        Ys[l] <- rnorm(1, thetas[j], sigma2s)
        l= l+1
      }
    }
    
    if (norm(median(Ys)-median(Y), type = "2") <= threshold){
      sigma2[i] = sigma2s
      mu[i] = mus
      tau2[i] = tau2s
      theta[i,] = thetas
      for (l in 1:t){
        R[l,] = Ys
      }
      
      nacp = nacp + 1
      
    }
    setTxtProgressBar(pb, i)
    
    
    
  }
  
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(list(sigma2 = sigma2, mu = mu, tau2 = tau2, theta = theta, R = R))
  
}


set.seed(1234)
niter = 200000

#threshold = 0.15 #migliore 0.14
#Res <- ABC_median(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold)

#threshold = 0.1
#Res <- ABC_mean(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold)

#threshold = 0.075
#Res <- ABC_mean_kernel_gauss(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold)

threshold = 0.05
Res <- ABC_mean_kernel_unif(niter, a0, b0, m0, s20, alpha0, beta0, data, threshold)



# Some plots

# Traceplots for mu, theta, tau2
x11()
par (mfrow=c(1,3))
plot(na.omit(Res$mu), type='l', xlab = "", ylab = "", main = "mu")
plot(na.omit(Res$theta[,1]), type='l', xlab = "", ylab = "", main = "theta[1]")
plot(na.omit(Res$tau2), type='l', xlab = "", ylab = "", main = "tau2")

Rmean <- apply(na.omit(Res$R), 2, mean)

# Values of estimated and true in the first iteration
x11()
par (mfrow=c(2,1))
plot(Rmean, type='l', ylab = "R_obs") #first estimate
plot(Y, type='l', ylab = "R") #ture values

# Values of Y true and estimated
x11()
par (mfrow=c(2,2))
hist(Y, breaks = 40, main = "", xlab = "", ylab = "", col="green" )
plot(density(Y), breaks = 40, main = "", xlab = "", ylab = "" , col="green" )
hist(Rmean, breaks = 40, main = "", xlab = "", ylab = "" , col="red" )
plot(density(Rmean), main = "", xlab = "", ylab = "" , col="red" )

