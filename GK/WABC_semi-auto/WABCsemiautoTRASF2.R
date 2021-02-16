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
  sumstat_ord = cbind(sort(sumstat[,1]), sort(sumstat[,2]), sort(sumstat[,3]), sort(sumstat[,4]))
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

x11()
par (mfrow=c(2,2))
plot(density(ris[,1]), main='A')
plot(density(ris[,2]), main='B')
plot(density(ris[,3]), main='g')
plot(density(ris[,4]), main='k')

