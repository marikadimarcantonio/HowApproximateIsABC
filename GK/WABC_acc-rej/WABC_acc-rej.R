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
x11()
out = mcmc(ys, 8000, model="gk", theta0=c(mean(ys),sd(ys),0,0),Sigma0=0.1*diag(4))

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

myabc = function(x, N, rprior, M, sumstats=c("all order statistics", "octiles", "moment estimates"), silent=FALSE) {
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
risultato<-myabc(ys, N=2.4*10^6, rprior=rprior, M=2048,sumstats="all order statistics")
## weighted euclidean distance ##
risultato <- abc(ys,N=2.4*10^6, "gk",rprior=rprior,M=2048,sumstats="all order statistics") 

x11()
par (mfrow=c(2,2))
plot((risultato[,1]), main='A', type='l', xlab='Iterations', ylab= '')
plot((risultato[,2]), main='B', type='l', xlab='Iterations', ylab= '')
plot((risultato[,3]), main='g', type='l', xlab='Iterations', ylab= '')
plot((risultato[,4]), main='k', type='l', xlab='Iterations', ylab= '')

x11()
par (mfrow=c(2,2))
plot(density(risultato[,1]), main='A')
plot(density(risultato[,2]), main='B')
plot(density(risultato[,3]), main='g')
plot(density(risultato[,4]), main='k')

graphics.off()
