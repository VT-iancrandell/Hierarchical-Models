####################################################################################
# FUNCTIONS
####################################################################################
do.Gibbs <- function(indat){
  library(msm) # rtnorm 
  library(mnormt) #rmnorm
  keep.trials <- indat[[1]]
  num.MCMC <- indat[[2]]
  phi.a <- indat[[3]]
  phi.b <- indat[[4]]
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
  missing.indicator <- (1:nrow(X))[X[,2] == 0]
  num.missing <- length(missing.indicator)
  Xstar <- matrix(0,num.MCMC,num.missing)
  Xstar[1,] <- rtnorm(num.missing,0,10,lower=0,upper=50)
  X[,2][missing.indicator] <- Xstar[1,]
  X[,4] <- X[,2] * X[,3]
  # Fairly strong prior, the tendency is for the random effects to march towards \infty and -\infty for some
  # fox that die on first day or survive the entire time.
  
  phi <- rep(1,num.MCMC)
  phi[1] <- rgamma(1,phi.a,phi.b)
  # ignoring repeated trials on same fox
  for (i in 2:num.MCMC){
    if (i %% 100 == 0){
      print(i)
      print(Sys.time())
    }
    # Sample Z_it
    mu <- X %*% beta.samples[(i-1),] + Q %*% theta[(i-1),]
    Z <- mapply(rtnorm,n=rep(1,indiv.trials),mean=mu,sd=rep(1,indiv.trials),lower=bottom,upper=up)
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
    inv.cov <- rep(beta.samples[i,2] ^2,num.missing) + X[missing.indicator,3] ^2 * beta.samples[i,4]^2 + rep(1/100,num.missing) 
    exp <- 1/inv.cov * (rep(beta.samples[i,2],num.missing) + X[missing.indicator,3] * beta.samples[i,4]) * 
      (Z[missing.indicator] - rep(beta.samples[i,1],num.missing) - X[missing.indicator,3] * beta.samples[i,3] - 
         Q[missing.indicator,] %*% theta[i,])
    Xstar[i,] <- rtnorm(num.missing,exp,sqrt(1/inv.cov), lower=rep(0,num.missing),upper <- rep(50,num.missing))
    X[missing.indicator,2] <- Xstar[i,]
    X[missing.indicator,4] <- X[missing.indicator,2] * X[missing.indicator,3]
  }
  return(list(beta.samples,theta,Xstar,phi))
}
####################################################################################
# Input & Process Data
####################################################################################
# Ignore these next two lines, they're for my working directory.
setwd("~/Desktop/Homework & Collaboration/STAT 5364/Fox Project")
in.dat <- read.csv('goodfox_jmp.csv')

in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
#including trial information from the fox that died during a trial
keep.trials <- in.dat[,c(1,2,3,4,6)]

keep.trials = transform(keep.trials, Date = as.Date(Date, format = "%m/%d/%Y"))
foxes = unique(keep.trials$Freq)
fox.by.day = numeric(0)
for(i in 1:27){
  tmp = keep.trials[keep.trials$Freq == foxes[i],, drop = F]
  dates  = seq(min(tmp$Date) - min(tmp$Experience), max(tmp$Date), by = 1)
  fox = rep(foxes[i], length(dates))
  exp = 0:(length(dates) - 1)
  dog = merge(data.frame(Date = dates), tmp, by = "Date", all.x = T)$Num_dogs
  dog[is.na(dog)] = 0
  dead = merge(data.frame(Date = dates), tmp, by = "Date", all.x = T)$Dead
  dead[is.na(dead)] = 0
  fox.by.day = rbind(fox.by.day, data.frame(Date = dates, Freq = fox, Experience = exp, Num_dogs = dog, dead = dead))
}

keep.trials = fox.by.day
# I believe experience =1 implies that is the foxes first day...
keep.trials <- keep.trials[keep.trials$Experience > 0,]

####################################################################################
# Package Data for Parallel MCMC Runs - looking for weird behavior...
####################################################################################
keep.trials = fox.by.day
num.MCMC <- 50000
phi.a <- 5
phi.b <- 5
indat <- list(keep.trials,num.MCMC,phi.a,phi.b)

# time.check <- rep(0,8)
# time.check[1] <- system.time(do.Gibbs(indat))[3]
# for (i in 2:8){
#   num.parallel <- i
#   in.par <- replicate(num.parallel,indat,simplify=F)
#   cl <- makeCluster(num.parallel)
#   time.check[i] <- system.time(output <- clusterApplyLB(cl, in.par, do.Gibbs))[3]
#   stopCluster(cl)
# }


require(parallel)
num.parallel <- detectCores()
in.par <- replicate(num.parallel,indat,simplify=F)
cl <- makeCluster(num.parallel)
output <- clusterApplyLB(cl, in.par, do.Gibbs)
stopCluster(cl)

####################################################################################
# VISUAL ASSESSMENT OF CONVERGENCE
####################################################################################
beta.samples <- theta <- Xstar <- phi<- NULL
for (i in 1:num.parallel){
  beta.samples <- rbind(beta.samples,output[[i]][[1]])
  theta <- rbind(theta,output[[i]][[2]])
  Xstar <- rbind(Xstar,output[[i]][[3]])
  phi <- c(phi,output[[i]][[4]])
}

par(mfcol=c(1,1))
#trace plots
plot(beta.samples[,1],type='l',main='Intercept')
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}

plot(beta.samples[,2],type='l',main='Number of Dogs')
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}

plot(beta.samples[,3],type='l',main = 'Fox Acclimation Time')
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}

plot(beta.samples[,4],type='l',main = 'Number of Dogs * Acclimation Time Interaction')
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}

plot(phi,type='l', main = expression(phi))
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}


#Posterior Means and Credible Intervals
apply(beta.samples,2,mean)
apply(beta.samples,2,quantile,probs=c(.025,.975))

apply(theta,2,mean)
par(mfcol=c(2,2))
for (i in 1:num.foxes){
  plot(theta[,i],type='l',main=paste('Fox id = ',foxes[i])) 
  for (i in 1:num.parallel){
    abline(v=num.MCMC*i)  
  }
}

par(mfcol=c(1,2))
hist(Xstar[,1])

plot(Xstar[,100],type='l')
for (i in 1:num.parallel){
  abline(v=num.MCMC*i)  
}

