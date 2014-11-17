library(msm) # rtnorm
library(mnormt) #rmnorm
in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
#including trial information from the fox that died during a trial
keep.trials <- in.dat[in.dat$Num_dogs >0,c(2,3,4,6)]
# still need to account for missing data on non-trial days
####################################################################
# Model #1
# Simple Probit with dogs only
####################################################################

# initialize Gibbs sampler
num.MCMC <- 5000
indiv.trials <- nrow(keep.trials)
alpha <- rep(0,num.MCMC)
beta <- rep(0,num.MCMC)
up <- rep(Inf,indiv.trials)
up[keep.trials$Dead ==1] <- 0
bottom <- rep(-Inf,indiv.trials)
bottom[keep.trials$Dead ==0] <- 0
X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs)

# Only consider dog effect ignoring experience and repeated trials on same fox
for (i in 2:num.MCMC){
  # Sample Z_it
  mu <- alpha[i] + beta[i] * keep.trials$Num_dogs
  Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
  # Sample beta & alpha jointly
  Cov.beta <- solve(t(X) %*% X )
  Exp.beta <- Cov.beta %*% t(X) %*% Z
  new.vals <- rmnorm(1,Exp.beta,Cov.beta)
  alpha[i+1] <- new.vals[1]
  beta[i+1] <- new.vals[2]
}

#trace plots
par(mfcol=c(1,2))
plot(alpha,type='l')
plot(beta,type='l')

#Posterior Means and Credible Intervals
mean(alpha); quantile(alpha,probs=c(.025,.975))
mean(beta); quantile(beta,probs=c(.025,.975))

#Verify moments near MLE
summary(glm(Survive~X -1,data=keep.trials,family=binomial(link="probit")))


# Compute "Mean" Hazard at each dog value
discrete_hazard <- cbind(sort(unique(keep.trials$Num_dogs)),1-pnorm(mean(alpha) + sort(unique(keep.trials$Num_dogs)) * mean(beta)))
colnames(discrete_hazard) <- c("Number of Dogs",'Hazard')
discrete_hazard

####################################################################
# Model 2
# Probit with dogs and experience
####################################################################

# initialize Gibbs sampler
num.MCMC <- 5000
indiv.trials <- nrow(keep.trials)
beta.samples <- matrix(0,nrow=num.MCMC,ncol=4)
up <- rep(Inf,indiv.trials)
up[keep.trials$Dead ==1] <- 0
bottom <- rep(-Inf,indiv.trials)
bottom[keep.trials$Dead ==0] <- 0
X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs,keep.trials$Experience,keep.trials$Num_dogs*keep.trials$Experience)

# ignoring repeated trials on same fox
for (i in 2:num.MCMC){
  # Sample Z_it
  mu <- X %*% beta.samples[i-1,]
  Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
  # Sample beta & alpha jointly
  Cov.beta <- solve(t(X) %*% X )
  Exp.beta <- Cov.beta %*% t(X) %*% Z
  beta.samples[i,] <- rmnorm(1,Exp.beta,Cov.beta)
}

#trace plots
plot(beta.samples[,1],type='l',main='Intercept')
plot(beta.samples[,2],type='l',main='Number of Dogs')
plot(beta.samples[,3],type='l',main = 'Fox Acclimation Time')
plot(beta.samples[,4],type='l',main = 'Number of Dogs * Acclimation Time Interaction')

#Posterior Means and Credible Intervals
apply(beta.samples,2,mean)
apply(beta.samples,2,quantile,probs=c(.025,.975))

#Verify moments near MLE
summary(glm(Survive~X -1,data=keep.trials,family=binomial(link="probit")))


# Try to simulate different policy regimes
# Compute "Mean" Hazard at each dog value
discrete_hazard <- cbind(X,round(1-pnorm(X %*% apply(beta.samples,2,mean)),5))
discrete_hazard <- discrete_hazard[,-c(1,4)]
colnames(discrete_hazard) <- c("Number of Dogs",'Fox Experience','Discrete Hazard')
discrete_hazard

####################################################################
# Model 3
# Probit with dogs and experience and random effects
####################################################################

# initialize Gibbs sampler
foxes <- unique(keep.trials$Freq)
num.foxes <- length(foxes)
num.MCMC <- 50000
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
# Fairly strong prior, the tendency is for the random effects to march towards \infty and -\infty for some
# fox that die on first day or survive the entire time.
phi.a <- phi.b <- 5
phi <- rep(1,num.MCMC)
# ignoring repeated trials on same fox
for (i in 2:num.MCMC){
  if (i %% 1000 == 0){
    print(i)
    print(Sys.time())
  }
  # Sample Z_it
  mu <- X %*% beta.samples[(i-1),] + Q %*% theta[(i-1),]
  Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
  # Sample beta & alpha jointly
  Cov.beta <- solve(t(X) %*% X )
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
}

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
for (i in 1:num.foxes){
  plot(theta[,i],type='l',main=paste('Fox id = ',foxes[i]))  
}

save(beta.samples,theta,phi,file='~/Dropbox/FoxData/MCMCout.Rdata')



####################################################################
# Policy Analysis
####################################################################
# Regime A.
# no limits on dog density or fox acclimation time
#
# Regime B
# No limits on dog density require at least 1 week of acclimation time before trials
#
# Regime C
# Max dog density and at least 1 week of acclimation time required before trials


## Regime a 
