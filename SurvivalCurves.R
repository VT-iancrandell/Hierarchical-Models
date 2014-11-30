load(file='~/Dropbox/Fox/Hierarchical-Models/parallelMCMC_noMissing.Rdata')
load(file='~/Dropbox/FoxData/keep.trials.Rdata')
# 
# keep.trials$DateVal <- strptime(keep.trials$Date,'%m/%d/%Y')
# keep.trials$RegA_Exp <- floor(as.numeric(keep.trials$DateVal - strptime('10/1/2002','%m/%d/%Y')))

####################################################################
# Combine MCMC iterations
####################################################################
#the first 4 come from b0 + b1*exp + b2*dogs - denote model 1, 
beta.1 <- theta.1 <- Xstar.1 <- phi.1 <- NULL
for (i in 1:4){
  beta.1 <- rbind(beta.1,output[[i]][[1]])  
  theta.1 <- rbind(theta.1,output[[i]][[2]])
#  Xstar.1 <- rbind(Xstar.1, output[[i]][[3]])
  phi.1 <- c(phi.1,output[[i]][[3]])
}

#the next 4 come from b0 + b1*log(exp) + b2*dogs - denote model 2, 
beta.2 <- theta.2 <- Xstar.2 <- phi.2 <- NULL
for (i in 5:8){
  beta.2 <- rbind(beta.2,output[[i]][[1]])  
  theta.2 <- rbind(theta.2,output[[i]][[2]])
  #Xstar.2 <- rbind(Xstar.2, output[[i]][[3]])
  phi.2 <- c(phi.2,output[[i]][[3]])
}
####################################################################
# Visualize Convergence
####################################################################
num.MCMC <- length(output[[1]][[3]])
par(mfcol=c(1,2))
##########################
plot(beta.1[,1],type='l',ylab=expression(beta[0]), main  = 'Model 1 w missing data')
for (i in 1:4){
  abline(v=num.MCMC*i,col='grey')
}
plot(beta.2[,1],type='l',ylab=expression(beta[0]), main  = 'Model 2 w missing data')
for (i in 1:4){
  abline(v=num.MCMC*i,col='grey')
}
##########################
plot(beta.1[,2],type='l',ylab=expression(beta['dogs']), main  = 'Model 1 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
plot(beta.2[,2],type='l',ylab=expression(beta['dogs']), main  = 'Model 2 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
##########################
plot(beta.1[,3],type='l',ylab=expression(beta['exp']), main  = 'Model 1 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
plot(beta.2[,3],type='l',ylab=expression(beta['log(exp)']), main  = 'Model 2 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
##########################
plot(phi.1,type='l',ylab=expression(phi), main  = 'Model 1 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
plot(phi.2,type='l',ylab=expression(phi), main  = 'Model 2 w missing data')
for (i in 1:6){
  abline(v=num.MCMC*i,col='grey')
}
##########################
num.theta <- 4
theta.vals <- sample(1:27,num.theta)
for (j in 1:num.theta){
  plot(theta.1[,theta.vals[j]],type='l',ylab=paste('theta = ',theta.vals[j]), main='Model 1 w missing data')
  for (i in 1:6){
    abline(v=num.MCMC*i,col='grey')
  }
  
  plot(theta.2[,theta.vals[j]],type='l',ylab=paste('theta = ',theta.vals[j]), main='Model 2 w missing data')
  for (i in 1:6){
    abline(v=num.MCMC*i,col='grey')
  }  
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
# Bayes Factor Calculations
####################################################################
# going to need X & Q
num.foxes <- length(unique(keep.trials$Freq))
foxes <- unique(keep.trials$Freq)
Q <- matrix(0,nrow=nrow(keep.trials),ncol=num.foxes)
for (i in 1:num.foxes){
  Q[,i] <- as.numeric(keep.trials$Freq == foxes[i])
}
X <- cbind(rep(1,nrow(keep.trials)),keep.trials$Num_dogs,keep.trials$Experience)

### CONTINUE HERE

beta <- beta.1; Xstar <- Xstar.1; theta <- theta.1
Y <- keep.trials$dead
compute_likelihood(X,Q,beta,theta,Xstar,Y){
  num.MCMC <- nrow(beta)
  num.obs <- nrow(X)
  X.impute <- X
  p.d <- numeric(num.MCMC)
  for (i in 1:num.MCMC){
    #impute missing values
    X.impute[X[,2] == 0,2] <- Xstar[i,]
    probs <- pnorm(X.impute %*% beta[i,] + Q %*% theta[i,])
    p.d[i] <- prod(probs^Y * (1-probs)^ (1-Y))
  }
  
}
prob.death = pnorm(X %*%beta.samples[i,] + Q%*%theta[i,])







####################################################################
# Policy Analysis
####################################################################
# Regime A.
# no limits on dog density or fox acclimation time, introduced on Oct. 1  
#
regA <- keep.trials[duplicated(keep.trials$Date) == FALSE,]
X <- cbind(rep(1,nrow(regA)),regA$Num_dogs,regA$RegA_Exp)
discrete_hazard <- 1-pnorm(X %*% apply(beta.samples,2,mean))
survival_prob <- 1-discrete_hazard
for (i in 2:nrow(regA)){
  survival_prob[i] <- survival_prob[i-1] * (1-discrete_hazard[i])
}

plot(survival_prob,type='l',lwd=2,main="Survival Probability: Regime A",ylim=c(0,1),
     ylab='Survival Prob',xlab='Trial')


x <- seq(0, max(regA$Num_dogs), by = 5)
y <- seq(0, max(regA$RegA_Exp), by = 5)
grid <- mesh(x, y)
beta.mean <- apply(beta.samples,2,mean)
z    <- with(grid, pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y))
z2    <- with(grid,1- pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y + beta.mean[4] * x * y))

par(mfrow = c(1,1))

persp3D(z = z2, x = x, y = y, ylab= 'Exp',xlab='Dogs',zlab='hazard')
persp3D(z = z2, x = x, y = y)
persp3D(z = z, x = x, y = y, facets = FALSE, curtain = TRUE)

# ribbons in two directions and larger spaces
ribbon3D(z = z2, x = x, y = y, along = "xy", space = 0.3)
ribbon3D(z = z2, x = x, y = y, along = "y", space = 0.3)

hist3D(z = z, x = x, y = y, border = "black")  




# Regime B
# No limits on dog density require at least 1 week of acclimation time before trials
#
regB <- regA
regB$Num_dogs[regB$RegA_Exp < 14] <-0
X <- cbind(rep(1,nrow(regB)),regB$Num_dogs,regB$RegA_Exp,regB$Num_dogs*regB$RegA_Exp)
discrete_hazard <- 1-pnorm(X %*% apply(beta.samples,2,mean))
survival_prob <- 1-discrete_hazard
for (i in 2:nrow(regA)){
  survival_prob[i] <- survival_prob[i-1] * (1-discrete_hazard[i])
}


lines(survival_prob, lty=2,col='red',lwd=2)

# Regime C
# Max dog density and at least 1 week of acclimation time required before trials


## Regime a 



# load(file='~/Dropbox/FoxData/MCMCout.Rdata')
# in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
# #including trial information from the fox that died during a trial
# keep.trials <- in.dat[in.dat$Num_dogs >0,c(1,3,4)]
# 
# keep.trials$DateVal <- strptime(keep.trials$Date,'%m/%d/%Y')
# keep.trials$RegA_Exp <- floor(as.numeric(keep.trials$DateVal - strptime('10/1/2002','%m/%d/%Y')))
# 
# ####################################################################
# # Policy Analysis
# ####################################################################
# # Regime A.
# # no limits on dog density or fox acclimation time, introduced on Oct. 1  
# #
# regA <- keep.trials[duplicated(keep.trials$Date) == FALSE,]
# X <- cbind(rep(1,nrow(regA)),regA$Num_dogs,regA$RegA_Exp,regA$Num_dogs*regA$RegA_Exp)
# discrete_hazard <- 1-pnorm(X %*% apply(beta.samples,2,mean))
# survival_prob <- 1-discrete_hazard
# for (i in 2:nrow(regA)){
#   survival_prob[i] <- survival_prob[i-1] * (1-discrete_hazard[i])
# }
# 
# plot(survival_prob,type='l',lwd=2,main="Survival Probability: Regime A",ylim=c(0,1),
#      ylab='Survival Prob',xlab='Trial')
# 
# 
# x <- seq(0, max(regA$Num_dogs), by = 5)
# y <- seq(0, max(regA$RegA_Exp), by = 5)
# grid <- mesh(x, y)
# beta.mean <- apply(beta.samples,2,mean)
# z    <- with(grid, pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y))
# z2    <- with(grid,1- pnorm(beta.mean[1] + beta.mean[2] * x + beta.mean[3] * y + beta.mean[4] * x * y))
# 
# par(mfrow = c(1,1))
# 
# persp3D(z = z2, x = x, y = y, ylab= 'Exp',xlab='Dogs',zlab='hazard')
# persp3D(z = z2, x = x, y = y)
# persp3D(z = z, x = x, y = y, facets = FALSE, curtain = TRUE)
# 
# # ribbons in two directions and larger spaces
# ribbon3D(z = z2, x = x, y = y, along = "xy", space = 0.3)
# ribbon3D(z = z2, x = x, y = y, along = "y", space = 0.3)
# 
# hist3D(z = z, x = x, y = y, border = "black")  
# 
# 
# 
# 
# # Regime B
# # No limits on dog density require at least 1 week of acclimation time before trials
# #
# regB <- regA
# regB$Num_dogs[regB$RegA_Exp < 14] <-0
# X <- cbind(rep(1,nrow(regB)),regB$Num_dogs,regB$RegA_Exp,regB$Num_dogs*regB$RegA_Exp)
# discrete_hazard <- 1-pnorm(X %*% apply(beta.samples,2,mean))
# survival_prob <- 1-discrete_hazard
# for (i in 2:nrow(regA)){
#   survival_prob[i] <- survival_prob[i-1] * (1-discrete_hazard[i])
# }
# 
# 
# lines(survival_prob, lty=2,col='red',lwd=2)
# 
# # Regime C
# # Max dog density and at least 1 week of acclimation time required before trials
# 
# 
# ## Regime a 
