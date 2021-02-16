library(winference)
library(gk)
library(Rmixmod)

###############################################################################

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



my_wsmc <- function(compute_d, rprior, nthetas, q, proposal, alpha, maxstep){
  
  ncomputed <- c(nthetas)

  thetas <- rprior(nthetas)

  dy <- list()
  for(itheta in 1:nthetas) {
    y_fake <- target$simulate(thetas[itheta,]) 
    distance <- compute_d(y_fake)
    dy[[itheta]] <- list(d = distance, y = y_fake)
  }
  distances <- sapply(X = dy, FUN = function(x) x$d)
  latest_y <- lapply(X = dy, FUN = function(x) x$y)
  threshold <- as.numeric(quantile(distances,probs = q))
  
  istep <- 1
  while (istep < maxstep){
    istep <- istep + 1
    one_step_results <- wsmc_one_step(thetas, distances, latest_y, compute_d, threshold, nthetas, rprior, proposal)
    thetas <- one_step_results$thetas
    distances <- one_step_results$distances
    latest_y <- one_step_results$latest_y
    threshold <- one_step_results$threshold
    
    cat( "  Step: ", istep, "\n", sep = "")
  } 
  return(theta = thetas)
  
}

systematic_resampling_ <- function(weights) {
  .Call('_winference_systematic_resampling_', PACKAGE = 'winference', weights)
}

systematic_resampling <- function(normalized_weights){
  return(systematic_resampling_(normalized_weights))
}

wsmc_one_step <- function(thetas, distances, latest_y, compute_d, threshold, nthetas, rprior, proposal){
  
  weights <- as.numeric(distances <= threshold)
  weights <- weights / sum(weights)
  # fit proposal
  proposal <- myupdate_proposal(thetas, weights, proposal) 
  # systematic resampling
  ancestors <- systematic_resampling(weights)
  thetas <- thetas[ancestors,,drop=F]
  distances <- distances[ancestors]
  old_latest_y <- latest_y
  for (i in 1:(dim(thetas)[1])){
    latest_y[[i]] <- old_latest_y[[ancestors[i]]]
  }
  weights <- rep(1/nthetas, nthetas)
  
  # rejuvenation moves
  acceptrates <- 0
  
  mh_res <- move_step(thetas, distances, latest_y, compute_d, target, threshold, rprior, proposal )
  
  thetas <- mh_res$thetas
  distances <- mh_res$distances
  latest_y <- mh_res$latest_y
  
  # re-fit proposal
  proposal <- myupdate_proposal(thetas, weights, proposal)
  acceptrates <- mh_res$acceptrate
  cat("  acceptance rates:", 100 * acceptrates, "%, threshold =", threshold,
      ", min. dist. =", min(distances), "\n")
  
  nunique <- length(unique(thetas[,1])) # number of unique thetas
  
  if (nunique/nthetas > alpha){
    one_uniform <- runif(1)
    g <- function(epsilon){
      logw <- log(as.numeric(distances <= epsilon))
      w <- exp(logw - max(logw))
      w <- w / sum(w)
      a <- systematic_resampling_given_u(w, one_uniform)
      th <- thetas[a,1]
      return((length(unique(th))-1)/nthetas)
    }
    lower_threshold <- 0
    upper_threshold <- threshold
    if (is.infinite(upper_threshold)){
      upper_threshold <- max(distances)
    }
    opt <- try(optimize(f = function(e) (g(e) - q)^2, interval = c(lower_threshold, upper_threshold)))
    if (!inherits(opt, "try-error"))
      threshold <- opt$minimum
  }
  return(list(threshold = threshold, thetas = thetas, distances = distances, latest_y = latest_y, acceptrates = acceptrates))
}

myupdate_proposal <- function(thetas, weights, proposal){
  proposal$param_prop <- proposal$param_update(thetas, weights)
  return(proposal)
}

move_step <- function(thetas, distances, latest_y, compute_d, target, threshold, rprior, proposal ){
  res_foreach <- list()
  for (i in 1:nrow(thetas)){
    theta <- thetas[i,]
    res_foreach[[i]] <- std_move_step_onetheta(theta, compute_d, target, threshold, rprior, proposal )
  }
  accepts <- sapply(res_foreach, function(x) x$accepted)
  for (i in 1:nrow(thetas)){
    if (accepts[i]){
      thetas[i,] <- res_foreach[[i]]$theta
      distances[i] <- res_foreach[[i]]$distance
      latest_y[[i]] <- res_foreach[[i]]$y
    }
  }
  return(list(acceptrate = mean(accepts), thetas = thetas, distances = distances, latest_y = latest_y))
}

std_move_step_onetheta <- function(theta, compute_d, target, threshold, rprior, proposal){
  theta <- matrix(theta, nrow = 1)
  prior_current <- target$dprior(theta, target$parameters)
  prop_current <- proposal$d(theta, proposal$param_prop)
  theta_prop <- proposal$r(theta, proposal$param_prop)
  prior_proposed <- target$dprior(theta_prop, target$parameters)
  prop_proposed <- proposal$d(theta_prop, proposal$param_prop)
  dproposed <- Inf
  y_prop <- 0
  if (!is.infinite(prior_proposed)){
    y_prop <- target$simulate(theta_prop)
    dproposed <- compute_d(y_prop)
  }
  logratio <- (prior_proposed - prior_current) + (prop_current - prop_proposed)
  logratio <- logratio + log(dproposed < threshold)
  accepted <- (log(runif(1)) < logratio)
  if (accepted){
    return(list(accepted = TRUE, theta = theta_prop, distance = dproposed, y = y_prop,
                nproposals = 1, ncurrent = 0))
  } else {
    return(list(accepted = FALSE, nproposals = 1, ncurrent = 0))
  }
}

###############################################################################

n = 250    
a = 3
b = 1
g = 2
k = 0.5

set.seed(1)

obs <- rgk(n,a,b,g,k,c=0.8)

sort_obs <- sort(obs)

compute_d = function(y){
  sort_y = sort(y)
  mean(abs(sort_y-sort_obs))
}

proposal <-  mixture_rmixmod()
target <- get_gandk()
target$simulate <- function(theta){
  return(matrix(target$robservation(n, theta, target$parameters, target$generate_randomness(n)), nrow = 1))
}

nthetas <- 2048
q <- 0.5
alpha = 0.5
maxstep <- 30
rprior = function(i) {cbind(runif(i,0,5), runif(i,0,5), runif(i,0,5), runif(i,0,5))}

ris <- my_wsmc(compute_d, rprior, nthetas, q, proposal, alpha, maxstep)

risA <- ris[,1]
risB <- ris[,2]
risg <- ris[,3]
risk <- ris[,4]

x11()
par (mfrow=c(2,2))
plot(density(risA), main='A')
plot(density(risB), main='B')
plot(density(risg), main='g')
plot(density(risk), main='k')

set.seed(1)
x <- rgk(10,a,b,g,k,c=0.8)

x11()
out = mymcmc(x, N=75000, theta0=c(mean(x),sd(x),0,0), Sigma0=0.1*diag(4))
out.new <- out[50000:75000,]


x11()
par (mfrow=c(2,2))
plot(density(out.new[,1]), main='A', ylim=c(0,4.2))
points(density(ris[,1]),col='red',type='l', main='A')
plot(density(out.new[,2]), main='B', ylim=c(0,2.1))
points(density(ris[,2]),col='red', type='l', main='B')
plot(density(out.new[,3]), main='g', xlim=c(-2,10.2))
points(density(ris[,3]),col='red', type='l', main='g')
plot(density(out.new[,4]), main='k', ylim=c(0,2.6))
points(density(ris[,4]),col='red', type='l', main='k')

