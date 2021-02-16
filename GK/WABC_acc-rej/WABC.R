library(gk)
library(Ecdat)
library(ggplot2)

# Parameters
n = 250    
a = 3
b = 1
g = 2
k = 0.5

# Observed data
ys <- rgk(n,a,b,g,k,c=0.8)

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
plot(density(out[,1]), main='A')
plot(density(out[,2]), main='B')
plot(density(out[,3]), main='g')
plot(density(out[,4]), main='k')

# Parameters' check
a_new <- 2.9
b_new <- 0.9
g_new <- 2.1
k_new <- 0.5

y_new <- rgk(n,a_new,b_new,g_new,k_new,c=0.8)
x11()
par (mfrow=c(1,2))
plot(density(ys),main='g-and-k distribution')
plot(density(y_new), main='Observed g-and-k distribution')

N_particles = 2048
N = 2.4*10^6
rprior = function(i) {
  cbind(runif(i,0,10), 
        runif(i,0,10), 
        runif(i,0,10), 
        runif(i,0,10))}

myabc = function(x, N, rprior, M, silent=FALSE) {
  nobs = length(x)
  # Define simStats: a function to simulate one set of summary statistics
  # and sobs: the observed summary statistics
  sobs = sort(x)
  simStats = function(theta) sort(rgk(nobs, A=theta[1], B=theta[2], g=theta[3], k=theta[4]))
  batch_size = 10^4
  nbatches = ceiling(N / batch_size)
  last_batch_size = N %% batch_size
  if (last_batch_size == 0) { last_batch_size = batch_size }
  if (!silent) { prog_bar = progress::progress_bar$new(total = nbatches, format = "[:bar] :percent eta: :eta") }
  batch_out = myabc_batch(sobs, rprior(batch_size), simStats, M)
  samp = batch_out$samp
  if (!silent) { prog_bar$tick() }
  next_batch_size = batch_size
  for (b in 2:nbatches) {
    if (b==nbatches) { next_batch_size = last_batch_size }
    next_samp = myabc_batch(sobs, rprior(next_batch_size), simStats, M)$samp
    samp = rbind(samp, next_samp)
    toacc = order(samp[,5])[1:M]
    samp = samp[toacc,]
    if (!silent) { prog_bar$tick() }
    
  }
  colnames(samp) = c("A","B", "g", "k", "distance")
  return(samp)
}
myabc_batch = function(sobs, priorSims, simStats, M) {
  summaries = apply(priorSims, 1, simStats)
  # Calculate distances
  d = apply(summaries, 2, function(s) { sum(abs(s-sobs))/N })
  # Construct and return output
  toacc = order(d)[1:M]
  samp = cbind(priorSims[toacc,], d[toacc])
  list(samp=samp)
}

## Wasserstein distance ##
risultato<-myabc(ys, N=2.4*10^6, rprior=rprior, M=2048)
## weighted euclidean distance ##
#risultato <- abc(ys,N=2.4*10^6, "gk",rprior=rprior,M=2048,sumstats="all order statistics") 

x11()
par (mfrow=c(2,2))
plot((risultato[,1]), main='A', type='l', xlab='Iterations', ylab= '')
plot((risultato[,2]), main='B', type='l', xlab='Iterations', ylab= '')
plot((risultato[,3]), main='g', type='l', xlab='Iterations', ylab= '')
plot((risultato[,4]), main='k', type='l', xlab='Iterations', ylab= '')

x11()
par (mfrow=c(2,2))
plot(density(out[,1]), main='A')
#legend(3.7, 1.2, legend=c("Target", "Result"),col=c("black", "red"), lty=1:1, cex=0.8)
points(density(risultato[,1]),col='red',type='l', main='A')
plot(density(out[,2]), main='B')
points(density(risultato[,2]),col='red', type='l', main='A')
plot(density(out[,3]), main='g')
points(density(risultato[,3]),col='red', type='l', main='A')
plot(density(out[,4]), main='k')
points(density(risultato[,4]),col='red', type='l', main='A')


#graphics.off()
