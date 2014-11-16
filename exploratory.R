library(msm) # rtnorm
library(mnormt) #rmnorm
in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
#only retain data for fox that died during a trial
keep.trials <- in.dat[in.dat$Num_dogs >0,c(2,3,4,6)]

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
  Cov.beta <- solve(t(X) %*% X + diag(2))
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
