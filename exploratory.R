library(msm) # rtnorm
library(mnormt) #rmnorm

rt_Discrete_CDF <- function(mu,y){
  lower <- rep(0,length(mu)); upper <- rep(1,length(mu))
  lower[y == 1] <- pnorm(0,mu[y==1])
  upper[y == 0] <- pnorm(0,mu[y==0])
  trunc.norm <- qnorm(runif(length(mu),lower,upper),mu)
  return(trunc.norm)
}


## simulate mixture prior figure
c0 <- .2
n <- 10000
po <- .15
r <- 3
indicator <- rbinom(n,1,c0)


nb.part <- rnbinom(n-sum(indicator),r,po)

hist(c(nb.part,rep(0,sum(indicator))),breaks='FD',prob=T,main='',xlab='Dogs')
# Ignore these next two lines, they're for my working directory.
#setwd("~/Desktop/Homework & Collaboration/STAT 5364/Fox Project")
#in.dat <- read.csv('goodfox_jmp.csv')

in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
#including trial information from the fox that died during a trial
keep.trials <- in.dat[,c(1,2,3,4,6)]

dead <- in.dat[in.dat$Dead == 1,]
max.experience <- aggregate(keep.trials$Experience,list(keep.trials$Freq),max)
alive.times <- max.experience[!max.experience[,1] %in% dead$Freq,][,2]
boxplot(dead$Experience,outline=F,ylim=c(0,max(c(dead$Experience,alive.times))),main='',ylab='Days')
x.jitter <- runif(21,.85,1.15)
points(x.jitter,dead$Experience,pch=16,cex=.8)
x.jitter2 <- runif(6,.85,1.15)
points(x.jitter2,alive.times,pch='X',cex=.7)


####################################################################
####################################################################
# Model #1
# Simple Probit with dogs only
####################################################################
# load(file='~/Dropbox/FoxData/keep.trials.Rdata')
# 
# 
# # initialize Gibbs sampler
# num.MCMC <- 5000
# indiv.trials <- nrow(keep.trials)
# alpha <- rep(0,num.MCMC)
# beta <- rep(0,num.MCMC)
# up <- rep(Inf,indiv.trials)
# up[keep.trials$dead ==1] <- 0
# bottom <- rep(-Inf,indiv.trials)
# bottom[keep.trials$dead ==0] <- 0
# X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs)
# Survive <- (keep.trials$dead - 1)^2
# # Only consider dog effect ignoring experience and repeated trials on same fox
# for (i in 2:num.MCMC){
#   # Sample Z_it
#   if (i %% 100 == 0){
#     print(i)
#     print(Sys.time())
#   }
#   
#   mu <- alpha[i] + beta[i] * keep.trials$Num_dogs
#   Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
#   # Sample beta & alpha jointly
#   Cov.beta <- solve(t(X) %*% X )
#   Exp.beta <- Cov.beta %*% t(X) %*% Z
#   new.vals <- rmnorm(1,Exp.beta,Cov.beta)
#   alpha[i+1] <- new.vals[1]
#   beta[i+1] <- new.vals[2]
# }
# 
# #trace plots
# par(mfcol=c(1,2))
# plot(alpha,type='l')
# plot(beta,type='l')
# 
# #Posterior Means and Credible Intervals
# mean(alpha); quantile(alpha,probs=c(.025,.975))
# mean(beta); quantile(beta,probs=c(.025,.975))
# 
# #Verify moments near MLE
# summary(glm(Survive~X -1,data=keep.trials,family=binomial(link="probit")))


# ####################################################################
# # Model 2
# # Probit with dogs and experience
# ####################################################################
# 
# # initialize Gibbs sampler
# num.MCMC <- 50000
# indiv.trials <- nrow(keep.trials)
# beta.samples <- matrix(0,nrow=num.MCMC,ncol=4)
# up <- rep(Inf,indiv.trials)
# up[keep.trials$dead ==1] <- 0
# bottom <- rep(-Inf,indiv.trials)
# bottom[keep.trials$dead ==0] <- 0
# X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs,keep.trials$Experience,keep.trials$Num_dogs*keep.trials$Experience)
# Survive <- (keep.trials$dead - 1)^2
# y <- Survive
# 
# 
# # ignoring repeated trials on same fox
# for (i in 2:num.MCMC){
#     # Sample Z_it
#   if (i %% 1000 == 0){
#     print(i)
#     print(Sys.time())
#   }
#   # Sample Z_it
#   mu <- X %*% beta.samples[i-1,]
#   #Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
#   Z <- rt_Discrete_CDF(mu,y)
#   # Sample beta & alpha jointly
#   Cov.beta <- solve(t(X) %*% X )
#   Exp.beta <- Cov.beta %*% t(X) %*% Z
#   beta.samples[i,] <- rmnorm(1,Exp.beta,Cov.beta)
# }
# 
# #trace plots
# par(mfcol=c(1,1))
# plot(beta.samples[,1],type='l',main='Intercept')
# plot(beta.samples[,2],type='l',main='Number of Dogs')
# plot(beta.samples[,3],type='l',main = 'Fox Acclimation Time')
# plot(beta.samples[,4],type='l',main = 'Number of Dogs * Acclimation Time Interaction')
# 
# #Posterior Means and Credible Intervals
# apply(beta.samples,2,mean)
# apply(beta.samples,2,quantile,probs=c(.025,.975))
# 
# #Verify moments near MLE
# summary(glm(Survive~X -1,data=keep.trials,family=binomial(link="probit")))
# 
# 

# ####################################################################
# # Model 3
# # Probit with dogs and experience and random effects
# ####################################################################
# 
# # initialize Gibbs sampler
# foxes <- unique(keep.trials$Freq)
# num.foxes <- length(foxes)
# num.MCMC <- 50000
# indiv.trials <- nrow(keep.trials)
# beta.samples <- matrix(0,nrow=num.MCMC,ncol=4)
# up <- rep(Inf,indiv.trials)
# up[keep.trials$Dead ==1] <- 0
# bottom <- rep(-Inf,indiv.trials)
# bottom[keep.trials$Dead ==0] <- 0
# theta <- matrix(0,nrow=num.MCMC,ncol=num.foxes)
# phi <- rep(0,num.MCMC)
# Q <- matrix(0,nrow=nrow(keep.trials),ncol=num.foxes)
# for (i in 1:num.foxes){
#   Q[,i] <- as.numeric(keep.trials$Freq == foxes[i])
# }
# 
# X <- cbind(rep(1,indiv.trials),keep.trials$Num_dogs,keep.trials$Experience,keep.trials$Num_dogs*keep.trials$Experience)
# # Fairly strong prior, the tendency is for the random effects to march towards \infty and -\infty for some
# # fox that die on first day or survive the entire time.
# phi.a <- phi.b <- 5
# phi <- rep(1,num.MCMC)
# # ignoring repeated trials on same fox
# for (i in 2:num.MCMC){
#   if (i %% 1000 == 0){
#     print(i)
#     print(Sys.time())
#   }
#   # Sample Z_it
#   mu <- X %*% beta.samples[(i-1),] + Q %*% theta[(i-1),]
#   Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
#   # Sample beta & alpha jointly
#   Cov.beta <- solve(t(X) %*% X )
#   Exp.beta <- Cov.beta %*% t(X) %*% (Z - Q %*% theta[(i-1),])
#   beta.samples[i,] <- rmnorm(1,Exp.beta,Cov.beta)
#   # Sample theta
#   Cov.theta <- solve(t(Q) %*% Q + phi[i-1]*diag(num.foxes))
#   Exp.theta <- Cov.theta %*% t(Q) %*% (Z - X %*% beta.samples[i,])
#   theta[i,] <- rmnorm(1,Exp.theta,Cov.theta)
#   theta[i,] <- theta[i,] - mean(theta[i,])
#   # Sample Phi
#   phiA <- phi.a + num.foxes/2
#   phiB <- phi.b + t(theta[i,]) %*% theta[i,] / 2
#   phi[i] <- rgamma(1,phiA,phiB)
# }
# 
# par(mfcol=c(1,1))
# #trace plots
# plot(beta.samples[,1],type='l',main='Intercept')
# plot(beta.samples[,2],type='l',main='Number of Dogs')
# plot(beta.samples[,3],type='l',main = 'Fox Acclimation Time')
# plot(beta.samples[,4],type='l',main = 'Number of Dogs * Acclimation Time Interaction')
# plot(phi,type='l', main = expression(phi))
# 
# #Posterior Means and Credible Intervals
# apply(beta.samples,2,mean)
# apply(beta.samples,2,quantile,probs=c(.025,.975))
# 
# apply(theta,2,mean)
# par(mfcol=c(2,2))
# for (i in 1:num.foxes){
#   plot(theta[,i],type='l',main=paste('Fox id = ',foxes[i]))  
# }
# 

####################################################################
# Model 4
# Probit with dogs and experience and random effects
# Include non trial days as zero dogs
####################################################################

# ##IAN PUT THE NEW DATA SET HERE
# keep.trials = transform(keep.trials, Date = as.Date(Date, format = "%m/%d/%Y"))
# foxes = unique(keep.trials$Freq)
# fox.by.day = numeric(0)
# for(i in 1:27){
#   tmp = keep.trials[keep.trials$Freq == foxes[i],, drop = F]
#   dates  = seq(min(tmp$Date) - min(tmp$Experience), max(tmp$Date), by = 1)
#   fox = rep(foxes[i], length(dates))
#   exp = 0:(length(dates) - 1)
#   dog = merge(data.frame(Date = dates), tmp, by = "Date", all.x = T)$Num_dogs
#   dog[is.na(dog)] = 0
#   dead = merge(data.frame(Date = dates), tmp, by = "Date", all.x = T)$Dead
#   dead[is.na(dead)] = 0
#   fox.by.day = rbind(fox.by.day, data.frame(Date = dates, Freq = fox, Experience = exp, Num_dogs = dog, dead = dead))
# }
# 
# 
# ##I'm overwriting the name of the old data set so avoid bugs.
# keep.trials = fox.by.day
# 


# initialize Gibbs sampler
foxes <- unique(keep.trials$Freq)
num.foxes <- length(foxes)
num.MCMC <- 50000
#num.MCMC <- 100
indiv.trials <- nrow(keep.trials)
beta.samples <- matrix(0,nrow=num.MCMC,ncol=4)
log.lik.samples = numeric(num.MCMC)
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
Cov.beta <- solve(t(X) %*% X )
# ignoring repeated trials on same fox
for (i in 2:num.MCMC){
  if (i %% 100 == 0){
    print(i)
    print(Sys.time())
  }
  # Sample Z_it
  mu <- X %*% beta.samples[(i-1),] + Q %*% theta[(i-1),]
  #Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
  Z <- rt_Discrete_CDF(mu,y)
  # Sample beta & alpha jointly
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
  # Compute -2 * log likelihood
  #prob.death = pnorm(X %*%beta.samples[i,] + Q%*%theta[i,])
 # log.lik.samples[i] = -2 * sum(keep.trials$Dead * log(prob.death) + (1 - keep.trials$Dead) * log(1 - prob.death))
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

# This code computes the MAP estimates for the parameters and includes a function to compute P(y_it = 1).

beta.hat = apply(beta.samples,2,mean)
theta.hat = apply(theta,2,mean)
surv.prob = matrix(0, nrow = 26, ncol = 412)
for(i in 1:26){
  for(j in 1:412){
    surv.prob[i,j] = pnorm(theta.hat[i] + X[j,] %*% beta.hat)
  }
}


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
