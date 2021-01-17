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

x <- seq(-1000, 1000, length = 100)
x11()
plot(x, dnorm(x, mean = mu1, sd = sig1), type="l", ylim =c(0,0.006))
points(x, dnorm(x, mean = mu1, sd = 70), type="l", col=2) 


x <- seq(-20, 20, length = 100)
x11()
plot(x, dnorm(x, mean = mu2, sd = sig2), type="l", ylim =c(0,0.2)) 
points(x, dnorm(x, mean = mu2, sd = 2), type="l", col=2) 

#-------------------------------------------
# PROPOSAL

#Uso i valori di sd trovati graficamente

#MAP <- c(24.513453, 2.734407)
MAP <- c(0, 0)
#Sigma <- rbind(c(33.0381982, 0.3279282),c(0.3279282, 0.1366159))
Sigma <- rbind(c(70, 0),c(0, 2)) #like in plots before
# PROPOSAL: Normal(X_current, Sigma)


#-------------------------------------------


ABC <- function(niter, burnin, 
                theta0, Sigma, mu1, mu2, sig1, sig2, data, threshold){
  # niter, burnin: iterations
  # theta0 initial state of the chain
  # Sig: proposal's var-cov matrix
  # mu1, mu2, sig2, sig2: hyperparameters
  
  theta_output <- matrix(nrow = (niter-burnin), ncol = 2)
  # number of accepted moves
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  # Start from theta0
  i=0
  while(i < niter)
  {
    theta_temp <- as.vector(rmvnorm(1, mean = MAP, sigma = Sigma))
    
    # Generate new data using model M with params theta_temp
    mu_sim <- dnorm(theta_temp[1], mu1, sig1)
    lambda_sim <- dnorm(theta_temp[2], mu2, sig2)
    Y <- rcauchy(length(data), mu_sim, lambda_sim) 
    
    if (norm(median(Y) - median(data), type="2") < threshold){
      ## Compute the probability of accepting the transition from th0 to theta_temp
      # First consider the accept/reject ratio. 
      # Numerator:
      acpn <- dnorm(theta_temp[1], mean = mu1, sd = sig1) * dnorm(theta_temp[2], mean = mu2, sd = sig2)
      # Denominator:
      acpd <- dmvnorm(theta_temp, mean = MAP, sigma = Sigma)
      
      acp <- acpn / acpd
      racp <- min(1, acp) 
      
      # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
      u <- runif(1)
      
      ## if u < acp accept the move
      if(u < racp)
      {
        theta0 <- theta_temp
        nacp = nacp + 1
      }
      # otherwise remain in the same state
      
      if(i > burnin )
      {
        theta_output[(i-burnin),] = theta0
      }
      setTxtProgressBar(pb, i)
      i = i + 1
    }
    
  }
  close(pb)
  cat("Acceptance rate =", nacp/niter, "\n")
  return(theta_output)
}

#------------------------------------
# RUN
# fix the seed
set.seed(42)

# fix the parameters of the Metropolis
niter = 5000
burnin = 0
threshold = 160

# At the end the chain will contain
# (niter-burnin) = 10000 observations
theta0 = MAP
theta_post <- ABC(niter = niter, burnin = burnin, theta0 = theta0, Sigma = Sigma, 
                  mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold)

#mcmc_theta_post <- mcmc(theta_post, start = burnin + 1, end = niter)
#summary(mcmc_theta_post) NON FUNZIONA

dati = data.frame(mu = theta_post[,1], lambda = theta_post[,2])
x11()  
ggplot(dati, aes(x = mu, y = lambda)) + geom_point(na.rm = TRUE)


x11()
par (mfrow=c(2,2))
plot(dati$mu, type='l', main = "mu")
plot(density(dati$mu[1:(dim(dati)[1]-1)]), main = "mu", xlab = "" )
plot(dati$lambda, type='l', main = "lambda")
plot(density(dati$lambda[1:(dim(dati)[1]-1)]), main = "lambda", xlab = "" )

