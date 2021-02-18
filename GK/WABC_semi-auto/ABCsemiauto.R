library(gk)
library(abctools)
library(ggplot2)
library(abc)

################################################################################################

MySemiAutoABC <- function(x, rprior, nobs, nsam, nacc, trasf, plot){
  
  #Sample from the prior
  sim_par <- rprior(nobs)
  #Create an artificial dataset based on the priors simulated
  artificial_dataset <- apply(sim_par, 1, function(x) rgk(nsam,x[1],x[2],x[3],x[4]))
  
  #Compute trasnsformation of the data
  x.tr <- eval(trasf[[1]](x))
  
  trasformazione <- list(function(x){rbind(x, x^2)})
  artificial_dataset.tr <- eval(trasformazione[[1]](artificial_dataset))
  # artificial_dataset.tr <- rbind(artificial_dataset,artificial_dataset^2)
  
  sobs = sort(x)
  simStats = function(theta) sort(rgk(nsam, A=theta[1], B=theta[2], g=theta[3], k=theta[4]))
  
  #Compute summary statistics for parameters
  sumstats = apply(sim_par, 1, simStats)
  sumstats.tr <- eval(trasformazione[[1]](sumstats))
  sumstats.tr <- t(sumstats.tr)
  
  d <- ncol(sim_par)
  p <- ncol(sumstats.tr)
  n <- nrow(sumstats.tr)
  
  sumstats.tr <- cbind(1,sumstats.tr)
  
  parnames <- paste("Parameter", 1:d)  
  varnames <- paste("X", 1:p, sep="")
  B0 <- c()
  B <- matrix(nrow=d, ncol=p) #nsim x ntrasf
  BICs <- c()
  
  #for each parameter (A, B, g, k) create a lm
  for (i in 1:d) { 
    reg <- lm.fit(sumstats.tr,sim_par[,i])
    class(reg) <- "lm" ##Needed to allow BIC calculation below 
    B0[i] <- reg$coefficients[1]
    B[i,] <- reg$coefficients[-1]
    BICs[i] <- AIC(reg, k=log(n))
  }
  
  if (any(is.na(B))) warning("Linear regression problems: ill-conditioning?")
  names(B0) <- parnames
  names(BICs) <- parnames
  rownames(B) <- parnames
  colnames(B) <- varnames
  sa <- list(B0=B0,B=B,BICs=BICs)  
  
  B <- sa$B
  B[is.na(B)] <- 0 ##NAs may exist due to collinearity of sumstats.tr[tobuild,] but can safely be set to zero
  ss.sa <- x.tr %*% t(B) #estimated summary statistics
  obs.sa <-  t(artificial_dataset.tr) %*% t(B) #estimated obs data  
  
  #Do ABC
  abcout.sa = myabc(x=ss.sa, nobs=nobs, M=nacc, sumstat=obs.sa, silent=FALSE)
  return(abcout.sa)
  
}


myabc = function(x, nobs, M, sumstat, silent=FALSE) {
  sumstat_ord = sumstat
  batch_size = 10^4
  nbatches = ceiling(nobs / batch_size)
  last_batch_size = nobs %% batch_size
  if (last_batch_size == 0) { last_batch_size = batch_size }
  if (!silent) { prog_bar = progress::progress_bar$new(total = nbatches, format = "[:bar] :percent eta: :eta") }
  batch_out = myabc_batch(x, param, sumstat_ord, M, nobs)
  samp = batch_out$samp
  if (!silent) { prog_bar$tick() }
  next_batch_size = batch_size
  for (b in 2:nbatches) {
    if (b==nbatches) { next_batch_size = last_batch_size }
    next_samp = myabc_batch(x, param, sumstat_ord, M, nobs)$samp
    samp = rbind(samp, next_samp)
    toacc = order(samp[,5])[1:M]
    samp = samp[toacc,]
    if (!silent) { prog_bar$tick() }
  }
  colnames(samp) = c("A","B", "g", "k", "distance")
  return(samp)
}

myabc_batch = function(x, param, sumstat_ord, nacc, nobs) {
  # Calculate distances
  d = apply(sumstat_ord, 1, function(s) { sum(abs(s-rep(x, nobs)))/nobs })
  # Construct and return output
  toacc = order(d)[1:nacc]
  samp = cbind(sumstat_ord[toacc,], d[toacc])
  list(samp=samp)
}

###################################################################

a = 3
b = 1
g = 2
k = 0.5

nobs <- 10^5 # num of simulations
nsam <- 250
nacc <- 2048 # num of accepted
x <- rgk(nsam, a, b, g, k)

rprior = function(i) {cbind(runif(i,0,5), runif(i,0,5), runif(i,0,5), runif(i,0,5))}

#transformations
trasf <- list(function(x){c(x, x^2)})

ris <- MySemiAutoABC(x, rprior, nobs, nsam, nacc, trasf, plot=TRUE)


improper_uniform_log_density = function(theta) {
  if (theta[2]<0 || theta[4]<0) return(-Inf)
  return(0)
}
mymcmc = function(x, N, model=c("gk", "gh"), logB=FALSE, get_log_prior=improper_uniform_log_density, theta0, Sigma0, t0=100, epsilon=1E-6, silent=FALSE, plotEvery=100) {
  if (!is.numeric(x)) stop("x must be numeric (a vector of observations)")
  if (!silent) { oldask = par(ask=FALSE) } ##Don't ask before progress plots
  output = matrix(nrow=N+1, ncol=4)
  colnames(output) = c("A", ifelse(logB, "log B", "B"), "g", ifelse(model[1]=="gk", "k", "h"))
  if (model[1] == "gk") {
    if (logB) {
      get_log_likelihood = function(theta) sum(dgk(x, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE))
    } else {
      get_log_likelihood = function(theta) sum(dgk(x, theta[1], theta[2], theta[3], theta[4], log=TRUE))
    }
  } else {
    if (logB) {
      get_log_likelihood = function(theta) sum(dgh(x, theta[1], exp(theta[2]), theta[3], theta[4], log=TRUE))
    } else {
      get_log_likelihood = function(theta) sum(dgh(x, theta[1], theta[2], theta[3], theta[4], log=TRUE))
    }
  }
  output[1,] = theta0
  theta = theta0
  Sigma = Sigma0
  C = chol(Sigma0)
  theta_bar = 0*theta0 ##Mean value of theta
  theta_mom2 = 0*Sigma ##2nd moment of theta
  if (!silent) { prog_bar = progress::progress_bar$new(total = N+1, format = "[:bar] :percent eta: :eta") }
  log_prior = get_log_prior(theta)
  log_likelihood = get_log_likelihood(theta)
  if (!silent) { prog_bar$tick() }
  for (i in 1:N) {
    theta_bar = theta_bar*(i-1)/i + theta/i
    theta_mom2 = theta_mom2*(i-1)/i + theta%*%t(theta)/i
    if (i > t0) {
      C = chol(theta_mom2 - theta_bar %*% t(theta_bar) + epsilon*diag(4))
    }
    theta_prop = theta + C %*% stats::rnorm(4)
    log_prior_prop = get_log_prior(theta_prop)
    if (log_prior_prop > -Inf) { ##Skip following if proposal has zero prior
      log_likelihood_prop = get_log_likelihood(theta_prop)
      r = log_prior_prop + log_likelihood_prop - log_prior - log_likelihood
      if (stats::runif(1) < exp(r)) {
        theta = theta_prop
        log_prior = log_prior_prop
        log_likelihood = log_likelihood_prop
      }
    }
    output[i+1,] = theta
    if (!silent && ((i+1) %% plotEvery == 0)) {
      graphics::par(mfrow=c(2,2))
      for (j in 1:4) {
        ylim = range(output[ceiling(i/10):(i+1),j])
        graphics::plot(output[,j], type='s', xlim=c(1,N), ylim=ylim, xlab="MCMC iteration", ylab=colnames(output)[j])
      }
    }
    if (!silent) { prog_bar$tick() }
  }
  if (!silent) { par(ask=oldask) }
  output
}
set.seed(1)
x = rgk(10, A=3, B=1, g=2, k=0.5) ##mall dataset for fast execution
x11()
out = mymcmc(x, N=75000, theta0=c(mean(x),sd(x),0,0), Sigma0=0.1*diag(4))
out = out[50000:75000,1:4]

x11()
par (mfrow=c(2,2))
plot(density(ris[,1]), main='A')
plot(density(ris[,2]), main='B')
plot(density(ris[,3]), main='g')
plot(density(ris[,4]), main='k')

x11()
par (mfrow=c(2,2))
plot(density(out[,1]), main='A')
#legend(3.7, 1.2, legend=c("Target", "Result"),col=c("black", "red"), lty=1:1, cex=0.8)
points(density(ris[,1]),col='red',type='l', main='A')
plot(density(out[,2]), main='B')
points(density(ris[,2]),col='red', type='l', main='A')
plot(density(out[,3]), main='g')
points(density(ris[,3]),col='red', type='l', main='A')
plot(density(out[,4]), main='k')
points(density(ris[,4]),col='red', type='l', main='A')
