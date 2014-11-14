in.dat <- read.csv('~/Dropbox/FoxData/goodfox_jmp.csv')
head(in.dat)

num.fox <- length(unique(in.dat$Freq))
head(in.dat,30)
num.trials <- length(unique(in.dat$Date))
