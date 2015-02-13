####################################################################################
# FUNCTIONS
####################################################################################
rt_Discrete_CDF <- function(mu,y){
  lower <- rep(0,length(mu)); upper <- rep(1,length(mu))
  lower[y == 1] <- pnorm(0,mu[y==1])
  upper[y == 0] <- pnorm(0,mu[y==0])
  trunc.norm <- qnorm(runif(length(mu),lower,upper),mu)
  return(trunc.norm)
}

mixture.like <- function(X,c=.2,p=.15,r=3){
  if (X==0) return(c)
  return((1-c)*dnbinom(X,r,p))
}


do.Gibbs <- function(indat){
  library(msm) # rtnorm 
  library(mnormt) #rmnorm
  rt_Discrete_CDF <- function(mu,y){
    lower <- rep(0,length(mu)); upper <- rep(1,length(mu))
    lower[y == 1] <- pnorm(0,mu[y==1])
    upper[y == 0] <- pnorm(0,mu[y==0])
    trunc.norm <- qnorm(runif(length(mu),lower,upper),mu)
    return(trunc.norm)
  }
  mixture.like <- function(X,c=.2,p=.15,r=3){
    if (X==0) return(c)
    return((1-c)*dnbinom(X,r,p))
  }
  keep.trials <- indat[[1]]
  num.MCMC <- indat[[2]]
  phi.a <- indat[[3]]
  phi.b <- indat[[4]]
  y <- indat[[5]]
  lambda <- 1
  # initialize Gibbs sampler
  foxes <- unique(keep.trials$Freq)
  num.foxes <- length(foxes)
  indiv.trials <- nrow(keep.trials)
  beta.samples <- matrix(0,nrow=num.MCMC,ncol=4)
  up <- rep(Inf,indiv.trials)
  up[keep.trials$Dead ==1] <- 0
  bottom <- rep(-Inf,indiv.trials)
  bottom[keep.trials$Dead ==0] <- 0
  theta <- matrix(0,nrow=num.MCMC,ncol=num.foxes)
  phi <- rep(0,num.MCMC)
  Q <- matrix(0,nrow=nrow(keep.trials),ncol=num.foxes)
  for (i in 1:num.foxes){
    Q[,i] <- as.numeric(keep.trials$Freq == foxes[i])
  }
  
  X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs,keep.trials$Experience,keep.trials$Num_dogs*keep.trials$Experience)
  X.propose <- X
  missing.indicator <- (1:nrow(X))[X[,2] == 0]
  num.missing <- length(missing.indicator)
  Xstar <- matrix(0,num.MCMC,num.missing)
  c0 <- .2
  po <- .15
  r <- 3
  indicator <- rbinom(num.missing,1,1-c0)
  nb.part <- rnbinom(num.missing-sum(1-indicator),r,po)
  indicator[indicator==1] <- nb.part
  
  Xstar[1,] <- indicator
  X[,2][missing.indicator] <- Xstar[1,]
  X[,4] <- X[,2] * X[,3]
  
  accept.rate <- matrix(0,num.MCMC,num.missing)
  
  # Fairly strong prior, the tendency is for the random effects to march towards \infty and -\infty for some
  # fox that die on first day or survive the entire time.
  
  phi <- rep(1,num.MCMC)
  phi[1] <- rgamma(1,phi.a,phi.b)
  # ignoring repeated trials on same fox
  for (i in 2:num.MCMC){
    if (i %% 500 == 0){
      print(i)
      print(Sys.time())
    }
    # Sample Z_it
    mu <- X %*% beta.samples[(i-1),] + Q %*% theta[(i-1),]
    Z <- rt_Discrete_CDF(mu,y)      
    #Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
    # Sample beta & alpha jointly
    Cov.beta <- solve(t(X) %*% X)
    Exp.beta <- Cov.beta %*% t(X) %*% (Z - Q %*% theta[(i-1),])
    beta.samples[i,] <- rmnorm(1,Exp.beta,Cov.beta)
    # Sample theta
    Cov.theta <- solve(t(Q) %*% Q + phi[i-1]*diag(num.foxes))
    Exp.theta <- Cov.theta %*% t(Q) %*% (Z - X %*% beta.samples[i,])
    theta[i,] <- rmnorm(1,Exp.theta,Cov.theta)
    theta[i,] <- theta[i,] - mean(theta[i,])
    # Sample Phi
    phiA <- phi.a + num.foxes/2
    phiB <- phi.b + t(theta[i,]) %*% theta[i,] / 2
    phi[i] <- rgamma(1,phiA,phiB)
    # Sample X*
    # random scan MH
    # randomize order
    rand.indices <- sample(1:num.missing,num.missing)
    num.partitions <- 350
    num.in.partition <- floor(num.missing/num.partitions)
    for (i2 in 1:num.partitions){
      start.at <- (i2-1) * num.in.partition + 1
      end.at <- i2*num.in.partition
      if (i2 == num.partitions) end.at = num.missing
      set.indices <- rand.indices[start.at:end.at]
      num.inset <- length(set.indices)
      change_mag <- rpois(num.inset,lambda)
      change_sign <- rbinom(num.inset,1,.5)
      change_sign[change_sign==0] <- -1
      Xnew <- Xstar[(i-1),][set.indices] + change_mag * change_sign
      Xnew[Xnew < 0] <- Xstar[(i-1),][set.indices][Xnew < 0]
      X.propose[missing.indicator,2][set.indices] <- Xnew
      X.propose[missing.indicator,4][set.indices] <- X.propose[missing.indicator,2][set.indices] * X.propose[missing.indicator,3][set.indices]
      log.alpha.propose <- sum(log(dnorm(Z[missing.indicator][set.indices],X.propose[missing.indicator,][set.indices,] %*% beta.samples[i,] + Q[missing.indicator,][set.indices,] %*% theta[i,],1) * sapply(Xnew,mixture.like)))
      log.alpha.old <- sum(log(dnorm(Z[missing.indicator][set.indices],X[missing.indicator,][set.indices,] %*% beta.samples[i,] + Q[missing.indicator,][set.indices,] %*% theta[i,],1) * sapply(Xstar[(i-1),][set.indices],mixture.like)))
      log.alpha <- log.alpha.propose - log.alpha.old
      if (log.alpha > log(runif(1))){
        accept.and.change <- as.numeric(Xstar[(i-1),set.indices] != Xnew)
        Xstar[i,set.indices] <- Xnew
        X[missing.indicator,2] <- Xstar[i,]
        X[missing.indicator,4] <- X[missing.indicator,2] * X[missing.indicator,3]    
        accept.rate[i,set.indices] <- accept.and.change
      } else {
        Xstar[i,set.indices] <- Xstar[(i-1),set.indices]
      }
    }    
#     change_mag <- rpois(num.missing,lambda)
#     change_sign <- rbinom(num.missing,1,.5)
#     change_sign[change_sign==0] <- -1
#     Xnew <- Xstar[(i-1),] + change_mag * change_sign
#     Xnew[Xnew < 0] <- Xstar[(i-1),Xnew < 0]
#     X.propose[missing.indicator,2] <- Xnew
#     X.propose[missing.indicator,4] <- X.propose[missing.indicator,2] * X.propose[missing.indicator,3]
#     log.alpha.propose <- sum(log(dnorm(Z[missing.indicator],X.propose[missing.indicator,] %*% beta.samples[i,] + Q[missing.indicator,] %*% theta[i,],1) * sapply(Xnew,mixture.like)))
#     log.alpha.old <- sum(log(dnorm(Z[missing.indicator],X[missing.indicator,] %*% beta.samples[i,] + Q[missing.indicator,] %*% theta[i,],1) * sapply(Xstar[(i-1),],mixture.like)))
#     log.alpha <- log.alpha.propose - log.alpha.old
#     if (log.alpha > log(runif(1))){
#       accept.rate[i] <- 1
#       Xstar[i,] <- Xnew
#       X[missing.indicator,2] <- Xstar[i,]
#       X[missing.indicator,4] <- X[missing.indicator,2] * X[missing.indicator,3]    
#     } else {
#       Xstar[i,] <- Xstar[(i-1),]
#     }
    # Old X* from trunc normal prior
#     inv.cov <- rep(beta.samples[i,2] ^2,num.missing) + X[missing.indicator,3] ^2 * beta.samples[i,4]^2 + rep(1/100,num.missing) 
#     exp <- 1/inv.cov * (rep(beta.samples[i,2],num.missing) + X[missing.indicator,3] * beta.samples[i,4]) * 
#       (Z[missing.indicator] - rep(beta.samples[i,1],num.missing) - X[missing.indicator,3] * beta.samples[i,3] - 
#          Q[missing.indicator,] %*% theta[i,])
#     Xstar[i,] <- rtnorm(num.missing,exp,sqrt(1/inv.cov), lower=rep(0,num.missing),upper <- rep(50,num.missing))
#     X[missing.indicator,2] <- Xstar[i,]
#     X[missing.indicator,4] <- X[missing.indicator,2] * X[missing.indicator,3]
  }
  return(list(beta.samples,theta,Xstar,phi,accept.rate))
}

####################################################################################
# Input & Process Data
####################################################################################
####################################################################################
# Assess Convergence
####################################################################################
load('~/Dropbox/FoxData/keep.trials.Rdata')
num.MCMC <- 50000
phi.a <- 5
phi.b <- 5
Survive <- (keep.trials$dead - 1)^2
#head(keep.trials)
keep.trials$Experience <- log(keep.trials$Experience)
indat <- list(keep.trials,num.MCMC,phi.a,phi.b,Survive)
 tmp <- do.Gibbs(indat)
# 
 beta.samples <- tmp[[1]]
 theta <- tmp[[2]]
 Xstar <- tmp[[3]]
 phi <- tmp[[4]]
accept.rate <- tmp[[5]]
# 
# save(beta.samples,theta,Xstar,phi,file='~/Dropbox/FoxData/fullmodel.Rdata')

# keep.trials.log <-keep.trials
# keep.trials.log$Experience <- log(keep.trials.log$Experience)
# indat <- list(keep.trials.log,num.MCMC,phi.a,phi.b,Survive)
# tmp2 <- do.Gibbs(indat)
# 
# beta.samples <- tmp2[[1]]
# theta <- tmp2[[2]]
# Xstar <- tmp2[[3]]
# phi <- tmp2[[4]]
# 
# save(beta.samples,theta,phi,file='~/Dropbox/Fox/Hierarchical-Models/MCMCout.Rdata')
# 
# save(beta.samples,theta,Xstar,phi,file='~/Desktop/fullmodel_log.Rdata')


####################################################################################
# VISUAL ASSESSMENT OF CONVERGENCE
####################################################################################

par(mfcol=c(1,1))
#trace plots
plot(beta.samples[,1],type='l',main='Intercept')

plot(beta.samples[,2],type='l',main='Number of Dogs')

plot(beta.samples[,3],type='l',main = 'Fox Acclimation Time')

plot(beta.samples[,4],type='l',main = 'Number of Dogs * Acclimation Time Interaction')

plot(phi,type='l', main = expression(phi))


#Posterior Means and Credible Intervals
apply(beta.samples,2,mean)
apply(beta.samples,2,quantile,probs=c(.025,.975))

apply(theta,2,mean)
par(mfcol=c(2,2))
foxes <- unique(keep.trials$Freq)
num.foxes <- length(foxes)
for (i in 1:num.foxes){
  plot(theta[,i],type='l',main=paste('Fox id = ',foxes[i])) 
}

par(mfcol=c(1,2))
hist((Xstar[,100]),breaks='FD')

plot(Xstar[,200],type='l')

colMeans(accept.rate)
