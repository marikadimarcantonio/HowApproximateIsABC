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


#-------------------------------------------
# PROPOSAL

# The target density p(??) is a standard normal density, truncated to the region ?? >1
# i.e. p(??) ??? exp(?????^2/2)I(?? >2)
# The majorizing function is a translated exponential density, m(??) = ?? exp[?????(?? ??? k)],
# with ?? ??? 1.62 and k = 1 (the truncation point).
## p.162 Wiley textbook ##



############ MU PROPOSAL ###########

alpha_mu<- 0.005
k_mu<-0.1
m_mu<-function(x){ifelse(x>=0, alpha_mu*exp(-alpha_mu*(x-k_mu)), alpha_mu*exp(alpha_mu*(x+k_mu)))}

#Plot proposal

x <- seq(-500, 500, length = 100)
x11()
plot(x,dnorm(x, mean = mu1, sd = sig1), type="l",ylim=c(0,0.006),ylab="",main="mu proposal") 
par(new=T)
plot(x, m_mu(x), type="l",col=2,ylim=c(0,0.006),ylab="") 
par(new=F)
legend("topright", legend=c("dnorm", "majorizing function"),
       col=c("black", "red"), lty=1:1, cex=0.8,
       box.lty=0)

#Plot proposal log scale (brutto)

x <- seq(0, 500, length = 100)
x11()
plot(x,dnorm(x, mean = mu1, sd = sig1), type="l", log='x', ylim=c(0,0.006),ylab="",main="mu proposal in log scale") 
lines(x,m_mu(x), type="l", log='x', col=2, ylim=c(0,0.06))
legend("topright", legend=c("dnorm", "majorizing function"),
       col=c("black", "red"), lty=1:1, cex=0.8,
       box.lty=0)


#creo un generatore di numeri random dalla nostra proposal m_theta

library(distr)
dist_mu <-AbscontDistribution(d=m_mu) # signature for a dist with pdf ~ p
rdist_mu <- r(dist_mu)

set.seed(1) # for reproduceable example
X <- rdist_mu(100) # sample from X ~ p
x <- seq(-1000,1000, .01)
x11()
hist(X, freq=F, breaks=50, xlim=c(-1000,1000))
lines(x,m_mu(x),lty=2, col="red")




############ LAMBDA PROPOSAL ###########

alpha_lambda <- 0.15
k_lambda <-0.1
m_lambda <-function(x){ifelse(x>=0, alpha_lambda *exp(-alpha_lambda*(x-k_lambda)), alpha_lambda*exp(alpha_lambda*(x+k_lambda)))}

#Plot proposal

x <- seq(-20, 20, length = 100)
x11()
plot(x,dnorm(x, mean = mu2, sd = sig2), type="l",ylim=c(0,0.15),ylab="",main="lambda proposal") 
par(new=T)
plot(x, m_lambda(x), type="l",col=2,ylim=c(0,0.15),ylab="") 
par(new=F)
legend("topright", legend=c("dnorm", "majorizing function"),
       col=c("black", "red"), lty=1:1, cex=0.8,
       box.lty=0)

#Plot proposal log scale (brutto)

x <- seq(0, 20, length = 100)
x11()
plot(x,dnorm(x, mean = mu2, sd = sig2), type="l", log='x', ylim=c(0,0.15),ylab="",main="lambda proposal in log scale") 
lines(x,m_lambda(x), type="l", log='x', col=2, ylim=c(0,0.15))
legend("topright", legend=c("dnorm", "majorizing function"),
       col=c("black", "red"), lty=1:1, cex=0.8,
       box.lty=0)


#creo un generatore di numeri random dalla nostra proposal m_theta

library(distr)
dist_lambda <-AbscontDistribution(d=m_lambda) # signature for a dist with pdf ~ p
rdist_lambda <- r(dist_lambda)

set.seed(1) # for reproduceable example
X <- rdist_lambda(100) # sample from X ~ p
x <- seq(-50, 50, .01)
x11()
hist(X, freq=F, breaks=50, xlim=c(-50, 50))
lines(x,m_lambda(x),lty=2, col="red")


########################################################


ABC <- function(niter, burnin, 
                      theta0, mu1, mu2, sig1, sig2, data, threshold){
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
    # propose a move to new theta -> theta_temp
    theta_temp <- c(rdist_mu(1), rdist_lambda(1)) #mu and lambda
    #Dovrebbero essere 2 distribuzioni diverse, bisogna calibrare i parametri delle due proposal
    
    # Generate new data using model M with params theta_temp
    mu_sim <- dnorm(theta_temp[1], mu1, sig1)
    lambda_sim <- dnorm(theta_temp[2], mu2, sig2)
    Y <- rcauchy(length(data), mu_sim, lambda_sim) # nel codice ta9 c'era dcauchy
    
    if (norm(Y - data,type="2") < threshold){
      ## Compute the probability of accepting the transition from th0 to theta_temp
      # First consider the accept/reject ratio. 
      # Numerator:
      acpn <- dnorm(theta_temp[1], mean = mu1, sd = sig1) * dnorm(theta_temp[2], mean = mu2, sd = sig2)
      # Denominator:
      acpd <- m_mu(theta_temp[1]) * m_lambda(theta_temp[2])
      
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
threshold = 160 # norm(Y - data,type="2") 167.0283

# At the end the chain will contain
# (niter-burnin) = 10000 observations
#theta0 = MAP
theta0 = c(0,0)
theta_post <- ABC(niter = niter, burnin = burnin, theta0 = theta0, 
                        mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2, data, threshold)


dati = data.frame(mu = theta_post[,1], lambda = theta_post[,2])
x11()  
ggplot(dati, aes(x = mu, y = lambda)) + geom_point(na.rm = TRUE)


x11()
par (mfrow=c(2,2))
plot(dati$mu, type='l', main = "mu")
plot(density(dati$mu[1:(dim(dati)[1]-1)]), main = "mu", xlab = "" )
plot(dati$lambda, type='l', main = "lambda")
plot(density(dati$lambda[1:(dim(dati)[1]-1)]), main = "lambda", xlab = "" )






# ONLY FOR PLOT

#-------------------------------------------
# log-posterior function -> pi(mu, lambda | y1,...yn)

cauchy_log_post <- function(theta, y, mu1, sig1, mu2, sig2){
  mu <- theta[1]
  lambda <- theta[2]
  n <- length(y)
  out <- 0
  
  # log-likelihood for each obs
  for(i in 1:n){
    out <- out - lambda - log(1 + exp(-2*lambda) * (y[i]-mu)^2)
  }
  
  # Add the prior in log scale (because log = T)
  out <- out + dnorm(mu, mean = mu1, sd = sig1, log = T)
  out <- out + dnorm(lambda, mean = mu2, sd = sig2, log = T)
  
  # return the log-posterior
  return(out)
}

#-------------------------------------------
# plot the log-posterior

grid_mu  <- seq(-40, 90, length = 100)
grid_lam <- seq(0.5, 6, length = 100)
expanded_grid <- as.matrix(expand.grid(grid_mu, grid_lam))

log_post <- c()
for(i in 1:nrow(expanded_grid)){
  log_post[i] <- cauchy_log_post(theta = expanded_grid[i,], y = data, 
                                 mu1 = mu1, sig1 = sig1, mu2 = mu2, sig2 = sig2)
}

plot_df <- data.frame(mu = expanded_grid[,1], lam = expanded_grid[,2], lpost = log_post)

#-------------------------------------------
# Visualize the proposal (fixing the position parameter equal to the MAP)

log_proposal <- c()
for(i in 1:nrow(expanded_grid)){
  log_proposal[i] <- dmvnorm(expanded_grid[i,], mean = MAP, sigma = Sigma, log = T)
}


plot_df_2 <- data.frame(mu = expanded_grid[,1], lam = expanded_grid[,2], lpost = log_proposal)
x11()
ggplot(plot_df, mapping = aes(x = mu, y = lam, z = lpost)) + 
  geom_contour() + 
  geom_contour(data = plot_df_2, col = 2, lty = 2) + 
  theme_bw() + 
  xlab("mu") + 
  ylab("lambda") + 
  ggtitle("log-posterior distribution") 










# COMPARISON MH


cauchy_MH <- function(niter, burnin, 
                      theta0, Sigma, mu1, mu2, sig1, sig2){
  # niter, burnin: iterations
  # theta0 initial state of the chain
  # Sig: proposal's var-cov matrix
  # mu1, mu2, sig2, sig2: hyperparameters
  
  # define the matrix that contains the output MCMC sample
  theta_output <- matrix(nrow = (niter-burnin), ncol = 2)
  # number of accepted moves
  nacp = 0  
  
  pb <- txtProgressBar(min = 1, max = niter, initial = 1, style = 3)
  # Start from theta0
  for(i in 1:(niter))
  {
    theta_temp <- as.vector(rmvnorm(1, mean = theta0, sigma = Sigma)) # propose theta_temp: the mean is given by the previous value
    
    ## Compute the logarithm of the probability of accepting the transition from th0 to theta_temp
    # First consider the log of accept/reject ratio. Numerator:
    lacp <- cauchy_log_post(th = theta_temp, y = data, mu1 = mu1, sig1 = sig1, mu2 = mu2, sig2 = sig2)
    # Denominator:
    lacp <- lacp - cauchy_log_post(th = theta0, y = data, mu1 = mu1, sig1 = sig1, mu2 = mu2, sig2 = sig2)
    lacp <- min(0, lacp)  
    # Note: The proposal is symmetrical and therefore does not appear in the value of acceptance / rejection!
    lgu <- log(runif(1))
    
    ## if u < acp accept the move
    if(lgu < lacp)
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
niter = 30000
burnin = 20000

# At the end the chain will contain
# (niter-burnin) = 5000 observations
MAP <- c(24.513453, 2.734407)
Sigma <- rbind(c(33.0381982, 0.3279282),c(0.3279282, 0.1366159))

theta0 = MAP
theta_post <- cauchy_MH(niter = niter, burnin = burnin, theta0 = theta0, Sigma = Sigma, 
                        mu1 = mu1, mu2 = mu2, sig1 = sig1, sig2 = sig2)

# CODA
mcmc_theta_post <- mcmc(theta_post, start = burnin + 1, end = niter)
summary(mcmc_theta_post)

#  The "naive" standard error is the standard error of the mean,
#  which captures simulation error of the mean rather than posterior
#  uncertainty.
#
#    naiveSE = posteriorSD/sqrt(n)
# 
#  The time-series standard error adjust the "naive" standard error for 
#  autocorrelation. ---> more precise estimate

x11()
plot(mcmc_theta_post)

