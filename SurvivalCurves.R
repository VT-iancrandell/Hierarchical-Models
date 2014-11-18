load(file='~/Dropbox/FoxData/MCMCout.Rdata')
in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
#including trial information from the fox that died during a trial
keep.trials <- in.dat[in.dat$Num_dogs >0,c(1,3,4)]

keep.trials$DateVal <- strptime(keep.trials$Date,'%m/%d/%Y')
keep.trials$RegA_Exp <- floor(as.numeric(keep.trials$DateVal - strptime('10/1/2002','%m/%d/%Y')))

####################################################################
# Policy Analysis
####################################################################
# Regime A.
# no limits on dog density or fox acclimation time, introduced on Oct. 1  
#
regA <- keep.trials[duplicated(keep.trials$Date) == FALSE,]
X <- cbind(rep(1,nrow(regA)),regA$Num_dogs,regA$RegA_Exp,regA$Num_dogs*regA$RegA_Exp)
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
