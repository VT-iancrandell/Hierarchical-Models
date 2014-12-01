load(file='~/Desktop/fullmodel_log.Rdata')
load(file='~/Dropbox/FoxData/keep.trials.Rdata')

load(file='~/Dropbox/Fox/Hierarchical-Models/MCMCout.Rdata')

#Xstar <- Xstar[1:200000,]
# 
# keep.trials$DateVal <- strptime(keep.trials$Date,'%m/%d/%Y')
# keep.trials$RegA_Exp <- floor(as.numeric(keep.trials$DateVal - strptime('10/1/2002','%m/%d/%Y')))

####################################################################
# Combine MCMC iterations
####################################################################
# #the first 4 come from b0 + b1*exp + b2*dogs - denote model 1, 
# beta.1 <- theta.1 <- Xstar.1 <- phi.1 <- NULL
# for (i in 1:4){
#   beta.1 <- rbind(beta.1,output[[i]][[1]])  
#   theta.1 <- rbind(theta.1,output[[i]][[2]])
# #  Xstar.1 <- rbind(Xstar.1, output[[i]][[3]])
#   phi.1 <- c(phi.1,output[[i]][[3]])
# }
# 
# #the next 4 come from b0 + b1*log(exp) + b2*dogs - denote model 2, 
# beta.2 <- theta.2 <- Xstar.2 <- phi.2 <- NULL
# for (i in 5:8){
#   beta.2 <- rbind(beta.2,output[[i]][[1]])  
#   theta.2 <- rbind(theta.2,output[[i]][[2]])
#   #Xstar.2 <- rbind(Xstar.2, output[[i]][[3]])
#   phi.2 <- c(phi.2,output[[i]][[3]])
# }
####################################################################
# Visualize Convergence
####################################################################
num.MCMC <- length(phi)
par(mfcol=c(1,1))
##########################
plot(beta.samples[,1],type='l',ylab=expression(beta[0]), main  = '')
# for (i in 1:4){
#   abline(v=num.MCMC*i,col='grey')
# }
##########################
plot(beta.samples[,2],type='l',ylab=expression(beta['dogs']), main  = '')
# for (i in 1:6){
#   abline(v=num.MCMC*i,col='grey')
# }
##########################
plot(beta.samples[,3],type='l',ylab=expression(beta['exp']), main  = '')
# for (i in 1:6){
#   abline(v=num.MCMC*i,col='grey')
# }
##########################
plot(beta.samples[,4],type='l',ylab=expression(beta['dogs*exp']), main  = '')
# for (i in 1:6){
#   abline(v=num.MCMC*i,col='grey')
# }
##########################
plot(phi,type='l',ylab=expression(phi), main  = 'Model 1 w missing data')
# for (i in 1:6){
#   abline(v=num.MCMC*i,col='grey')
# }
##########################
num.theta <- 4
theta.vals <- sample(1:27,num.theta)
for (j in 1:num.theta){
  plot(theta[,theta.vals[j]],type='l',ylab=paste('theta = ',theta.vals[j]), main='Model 1 w missing data')
#   for (i in 1:6){
#     abline(v=num.MCMC*i,col='grey')
#   }
}
##########################
# num.xstar <- 4
# xstar.vals <- sample(1:ncol(Xstar.1),num.xstar)
# par(mfcol=c(2,2))
# for (j in 1:num.xstar){
#   plot(Xstar.1[,xstar.vals[j]],type='l',ylab=paste('Xstar = ',xstar.vals[j]), main='Model 1 w missing data')
#   for (i in 1:6){
#     abline(v=num.MCMC*i,col='grey')
#   }
#   hist(Xstar.1[,xstar.vals[j]],main='Model 1 w missing data',xlab='Xstar')
#   
#   plot(Xstar.2[,xstar.vals[j]],type='l',ylab=paste('Xstar = ',xstar.vals[j]), main='Model 2 w missing data')
#   for (i in 1:6){
#     abline(v=num.MCMC*i,col='grey')
#   }  
#   hist(Xstar.2[,xstar.vals[j]],main='Model 1 w missing data',xlab='Xstar')
#   
# }
####################################################################
# Response Surface
####################################################################

x <- seq(0, max(keep.trials$Num_dogs), by = 5)
y <- seq(0, max(log(keep.trials$Experience)), by = .1)
library(plot3D)
grid <- mesh(x, y)
beta.mean <- apply(beta.samples,2,mean)
z    <- with(grid, pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y + beta.mean[4] * x * y))
z2    <- with(grid,beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y + beta.mean[4] * x * y)
z3    <- with(grid, 1-pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y + beta.mean[4] * x * y))

par(mfrow = c(1,1))

persp3D(z = z2, x = x, y = y, ylab= 'Exp',xlab='Dogs',zlab='Latent Survival')
persp3D(z = z, x = x, y = y)
persp3D(z = z, x = x, y = y, facets = FALSE, curtain = TRUE)
persp3D(z = z, x = x, y = y, facets = FALSE)
persp3D(z = z, x = x, y = y, curtain = TRUE,ylab= 'log(Exp)',xlab='Dogs',zlab='Survival Probability')
persp3D(z = z3, x = x, y = y, curtain = TRUE,ylab= 'log(Exp)',xlab='Dogs',zlab='Discrete Hazard')


####################################################################
# Bayes Factor Calculations
####################################################################
# # going to need X & Q
# num.foxes <- length(unique(keep.trials$Freq))
# foxes <- unique(keep.trials$Freq)
# Q <- matrix(0,nrow=nrow(keep.trials),ncol=num.foxes)
# for (i in 1:num.foxes){
#   Q[,i] <- as.numeric(keep.trials$Freq == foxes[i])
# }
# X <- cbind(rep(1,nrow(keep.trials)),keep.trials$Num_dogs,keep.trials$Experience)
# 
# ### CONTINUE HERE
# 
# beta <- beta.1; Xstar <- Xstar.1; theta <- theta.1
# Y <- keep.trials$dead
# compute_likelihood(X,Q,beta,theta,Xstar,Y){
#   num.MCMC <- nrow(beta)
#   num.obs <- nrow(X)
#   X.impute <- X
#   p.d <- numeric(num.MCMC)
#   for (i in 1:num.MCMC){
#     #impute missing values
#     X.impute[X[,2] == 0,2] <- Xstar[i,]
#     probs <- pnorm(X.impute %*% beta[i,] + Q %*% theta[i,])
#     p.d[i] <- prod(probs^Y * (1-probs)^ (1-Y))
#   }
#   
# }
# prob.death = pnorm(X %*%beta.samples[i,] + Q%*%theta[i,])
# 






####################################################################
# Policy Analysis
####################################################################
# Regime A.
# no limits on dog density or fox acclimation time, introduced on Oct. 1  
#
regA <- keep.trials[duplicated(keep.trials$Date) == FALSE,]
tmpX <- keep.trials[keep.trials$Num_dogs == 0,]
pos <- (1:nrow(tmpX))[!duplicated(tmpX$Date)]
X <- cbind(rep(1,nrow(regA)),regA$Num_dogs,log(regA$Experience), log(regA$Experience) * regA$Num_dogs)

surv.prob <- matrix(1,num.MCMC,nrow(X))

X.val <- X
for (j in 1:num.MCMC){
  X.val[X.val[,2] == 0,2] <- Xstar[j,pos]
  X.val[,4] <- X.val[,2] * X.val[,3]
  surv.prob[j,] <- pnorm(X.val %*% beta.samples[j,])
}
cum.surv <- matrix(1,num.MCMC,nrow(X))

cum.surv[,1] <- surv.prob[,1]
for (i in 2:nrow(X)){
  cum.surv[,i] <- cum.surv[,(i-1)] * surv.prob[,i]
}

plot(colMeans(cum.surv),type='l',lwd=2,main="Survival Probability: Regime A",ylim=c(0,1),
     ylab='Survival Prob',xlab='Trial')
lower <- apply(cum.surv,2,quantile,probs=.025)
lines(lower,lwd=2,lty=2)
upper <- apply(cum.surv,2,quantile,probs=.975)
lines(upper,lwd=2,lty=2)


# Regime B
# No limits on dog density require at 2 weeks of acclimation time before trials
#
regB <- regA
tmpX <- keep.trials[keep.trials$Num_dogs == 0,]
pos <- (1:nrow(tmpX))[!duplicated(tmpX$Date)]
#regB$Num_dogs[regB$Experience < 14] <-0
X <- cbind(rep(1,nrow(regB)),regB$Num_dogs,log(regB$Experience), log(regB$Experience) * regB$Num_dogs)

surv.prob <- matrix(1,num.MCMC,nrow(X))

for (j in 1:num.MCMC){
  X.val <- X
  X.val[X.val[,2] == 0,2] <- Xstar[j,pos]
  X.val[1:14,2] <- 0
  X.val[,4] <- X.val[,2] * X.val[,3]
  surv.prob[j,] <- pnorm(X.val %*% beta.samples[j,])  
}

cum.surv <- matrix(1,num.MCMC,nrow(X))

cum.surv[,1] <- surv.prob[,1]
for (i in 2:nrow(X)){
  cum.surv[,i] <- cum.surv[,(i-1)] * surv.prob[,i]
}

plot(colMeans(cum.surv),type='l',lwd=2,main="Survival Probability: Regime B",ylim=c(0,1),
     ylab='Survival Prob',xlab='Trial')
lower <- apply(cum.surv,2,quantile,probs=.025)
lines(lower,lwd=2,lty=2)
upper <- apply(cum.surv,2,quantile,probs=.975)
lines(upper,lwd=2,lty=2)
# Regime C
# Max dog density 400 dogs and 2 week acclimation period
regC <- regA
regC$Num_dogs[regC$Num_dogs > 400] <- 400

X <- cbind(rep(1,nrow(regC)),regC$Num_dogs,log(regC$Experience), log(regC$Experience) * regC$Num_dogs)

surv.prob <- matrix(1,num.MCMC,nrow(X))

X.val <- X
for (j in 1:num.MCMC){
  X.val[X.val[,2] == 0,2] <- Xstar[j,pos]
  X.val[,4] <- X.val[,2] * X.val[,3]
  surv.prob[j,] <- pnorm(X.val %*% beta.samples[j,])
}
cum.surv <- matrix(1,num.MCMC,nrow(X))

cum.surv[,1] <- surv.prob[,1]
for (i in 2:nrow(X)){
  cum.surv[,i] <- cum.surv[,(i-1)] * surv.prob[,i]
}

plot(colMeans(cum.surv),type='l',lwd=2,main="Survival Probability: Regime C",ylim=c(0,1),
     ylab='Survival Prob',xlab='Trial')
lower <- apply(cum.surv,2,quantile,probs=.025)
lines(lower,lwd=2,lty=2)
upper <- apply(cum.surv,2,quantile,probs=.975)
lines(upper,lwd=2,lty=2)
# Regime D
# Max dog density 400 dogs and 2 week acclimation period
regD <- regA
regD$Num_dogs[regD$Num_dogs > 400] <- 400
#regD$Num_dogs[regD$Experience < 14] <-0

X <- cbind(rep(1,nrow(regD)),regD$Num_dogs,log(regD$Experience), log(regD$Experience) * regD$Num_dogs)

surv.prob <- matrix(1,num.MCMC,nrow(X))

for (j in 1:num.MCMC){
  X.val <- X
  X.val[X.val[,2] == 0,2] <- Xstar[j,pos]
  X.val[1:14,2] <- 0
  X.val[,4] <- X.val[,2] * X.val[,3]
  surv.prob[j,] <- pnorm(X.val %*% beta.samples[j,])  
}
cum.surv <- matrix(1,num.MCMC,nrow(X))

cum.surv[,1] <- surv.prob[,1]
for (i in 2:nrow(X)){
  cum.surv[,i] <- cum.surv[,(i-1)] * surv.prob[,i]
}

plot(colMeans(cum.surv),type='l',lwd=2,main="Survival Probability: Regime D",ylim=c(0,1),
     ylab='Survival Prob',xlab='Trial')
lower <- apply(cum.surv,2,quantile,probs=.025)
lines(lower,lwd=2,lty=2)
upper <- apply(cum.surv,2,quantile,probs=.975)
lines(upper,lwd=2,lty=2)