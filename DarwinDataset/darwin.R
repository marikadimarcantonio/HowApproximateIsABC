library(LearnBayes)
library(mvtnorm)
library(coda)  
library(ggplot2)
library(ggpubr)

############################################
### METROPOLIS-HASTINGS ALGORITHM      #####
############################################
### CAUCHY DISTRIBUTION                #####
############################################

data(darwin)
data = darwin$difference
data

#-------------------------------------------
# prior distribution: product of two independent gaussian distribution
# specify the hyperparameters

mu1 <- 0
sig1 <- 100
mu2 <- 0
sig2 <- sqrt(10)

#Prior as proposal

###############   ABC_median   #####################
#

ABC_median <- function(niter = niter, mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold){
  
  theta_output <- matrix(nrow = niter, ncol = 2)
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  i=0
  N = length(data)
  while(i < niter)
  {
    mu_sim <- rnorm(1, mu1, sig1)
    lambda_sim <- rnorm(1, mu2, sig2)
    Y <- rcauchy(N, mu_sim, abs(lambda_sim)) #devo mettere l'abs perchè la cauchy non vuole valori negativi li!!!
    
    if (norm(median(Y) - median(data),type="2") < threshold){
      nacp = nacp + 1
      theta_output[i,] = c(mu_sim, lambda_sim)
    }
    
    setTxtProgressBar(pb, i)
    i = i + 1
  }
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(theta_output)
}

ABC_median_kernel <- function(niter = niter, mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold){
  
  theta_output <- matrix(nrow = niter, ncol = 2)
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  
  i=0
  N = length(data)
  while(i < niter)
  {
    mu_sim <- rnorm(1, mu1, sig1)
    lambda_sim <- rnorm(1, mu2, sig2)
    Y <- rcauchy(N, mu_sim, abs(lambda_sim)) #devo mettere l'abs perchè la cauchy non vuole valori negativi li!!!
    
    u <- runif(1)
    if (u < dnorm(norm(median(Y)-median(data), type="2")/threshold)/(threshold)){
      nacp = nacp + 1
      theta_output[i,] = c(mu_sim, lambda_sim)
    }
    
    setTxtProgressBar(pb, i)
    i = i + 1
  }
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(theta_output)
}

set.seed(1234)
niter = 200000

threshold = 1.5
theta_post <- ABC_median(niter = niter, mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold)

#threshold = 100
#theta_post <- ABC_median_kernel(niter = niter, mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold)

theta_post <- na.omit(theta_post)

dati = data.frame(mu = theta_post[,1], lambda = theta_post[,2])
x11()  
ggplot(dati, aes(x = mu, y = lambda)) + geom_point(na.rm = TRUE)


x11()
par (mfrow=c(2,2))
plot(dati$mu, type='l', main = "mu")
plot(density(dati$mu), main = "mu", xlab = "" )
plot(dati$lambda, type='l', main = "lambda")
plot(density(dati$lambda), main = "lambda", xlab = "" )



