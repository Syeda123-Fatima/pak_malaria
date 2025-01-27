# LOAD REQUIRED LIBRARIES
library(dlnm)
library(mvmeta)
library(splines)
library(MASS)
library(tsModel)
library(lubridate)
library(tidyverse)
library(readr)
library(haven)
library(performance)
library(ggplot2)
library(dplyr)
library(MuMIn)
library(gridExtra)

# LOAD AND SET UP DATA
setwd("your/directory/path")
Lakki <- read.csv("MalariaPak.csv", stringsAsFactors = FALSE)
Lakki$date <- dmy(Lakki$Date)
Lakki$date <- as.Date(Lakki$date, format = "%d/%m/%y")

# SET MODEL PARAMETERS
lag <- 3
arglag <- list(fun="integer")
cb <- crossbasis(Lakki$Tmean, lag=3, argvar=list(fun="thr",thr=22.4), arglag=list(fun="integer"))
plot(cb)

# DEFINE AND RUN GLM MODEL
model <- glm(total_cases ~ cb + ns(Lakki$date, df=2*9) + offset(log(Population)) + ns(Lakki$Relative_Humidity, df=2) + ns(Lakki$Precipitation, df=2), Lakki, family=quasipoisson (link="log"), maxit = 50, na.action="na.exclude")

# DEFINE QAIC CALCULATION FUNCTION
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
} 
fqaic(model)

# PERFORM CROSS PREDICTION
cp <- crosspred(cb, model,at = 22.4:35.3,by=0.1)
pred <- crosspred(cb, model, at = 22.4:35.3, by=0.1) 
allRRfit <- (pred$allRRfit-1)*100
allRRfitlow <- (pred$allRRlow-1)*100
allRRfithigh <- (pred$allRRhigh-1)*100

# SET UP DATA FRAME FOR RISK PREDICTIONS
df <- data.frame(matrix(nrow=length(pred$predvar), ncol=4))
colnames(df) <- c("temp", "RRfit", "RRlow", "RRhigh")
df$temp <- pred$predvar
df$RRfit <- (pred$allRRfit-1)*100
df$RRlow <- (pred$allRRlow-1)*100
df$RRhigh <- (pred$allRRhigh-1)*100

# PLOT EXPOSURE-RESPONSE ASSOCIATION
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
predvar <- quantile(Lakki$EHF,1:99/100,na.rm=T)
xlab <- expression(paste("Temperature (",degree,"C)"))

pdf("Figure 1.pdf",height=4.5,width=10)
layout(matrix(c(1,2),ncol=2,byrow=T))
par(mar=c(2,3,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"3d",ltheta=150,xlab="Temperature (C)",ylab="Lag",zlab="RR", 
     col=gray(0.9), main="Exposure-lag-response")
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.9,4.0),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",main="Overall")

# CROSS-REDUCTION AND SETUP FOR MONTE CARLO SIMULATION
red <- crossreduce(cb,model,at = 22.4:35.3,by=0.1)
coef <- coef(red)
vcov <- vcov(red)
oi <- rep(Lakki$total_cases, length=nrow(Lakki))
oimm <- tapply(Lakki$total_cases,as.numeric(format(Lakki$date,"%j")),
               mean,na.rm=T)[seq(12)]
oiperiod <- sum(oimm)*9

# DEFINE BASELINE PERIOD FOR ATTRIBUTABLE FRACTION ESTIMATION
baselineperiod <- "2014-2022"
baselineseqperiod <- factor(rep(baselineperiod,length.out=12*9))
seqperiod <- factor(c(as.numeric(baselineseqperiod)))
levels(seqperiod) <- c(baselineseqperiod)
temprange <- c("tot")
absrel <- c("abs")

# SET UP MONTE CARLO SIMULATION
nsim <- 1000
ansim_bs <- array(NA,dim=c(length(levels(seqperiod)),length(temprange),
                           length(absrel), nsim+1), 
                  dimnames=list(levels(seqperiod),temprange,absrel,
                                c("est",paste0("sim",seq(nsim)))))

# DEFINE AND INITIALIZE SIMULATION PARAMETERS
argvar=list(fun="thr",thr=22.4)
cenvec <- do.call(onebasis,c(list(x=22.4),argvar))
bvar <- do.call(onebasis,c(list(x=Lakki$Tmean),argvar))
bvarcen <- scale(bvar,center=cenvec,scale=F)

# INDICATOR FOR OPTIMAL TEMPERATURE RANGE FOR MALARIA TRANSMISSION
indheat <- which(Lakki$Tmean >= 22.4 & Lakki$Tmean <= 35.3)

# COMPUTE MONTHLY ATTRIBUTABLE NUMBERS
anbaseline <- (1-exp(-bvarcen%*%coef(red)))*oi
ansim_bs[,"tot","abs",1] <- tapply(anbaseline[indheat],factor(seqperiod[indheat]),sum)

# RUN MONTE CARLO SIMULATIONS
set.seed(13041975)
coefsim <- mvrnorm(nsim,coef,vcov)
for(s in seq(nsim)) {
  anbaseline <- (1-exp(-bvarcen%*%coefsim[s,]))*oi
  ansim_bs[,"tot","abs",s+1] <- tapply(anbaseline[indheat],factor(seqperiod[indheat]),sum)
}

# CALCULATE ESTIMATION INTERVALS FOR ATTRIBUTABLE NUMBERS
estci <- c("est","ci.l","ci.u")
anabs_bs <- afabs_bs <- array(NA,dim=c(length(levels(seqperiod)),
                                       length(estci),length(temprange)), 
                              dimnames=list(levels(seqperiod),estci,temprange))
anabs_bs[,"est",] <- apply(ansim_bs,1:2,mean)
anabs_bs[,"ci.l",] <- apply(ansim_bs,1:2,quantile,0.025)
anabs_bs[,"ci.u",] <- apply(ansim_bs,1:2,quantile,0.975)

# CALCULATE ATTRIBUTABLE FRACTION
afabs_bs[,,] <- anabs_bs[,,]/oiperiod*100
afabs_bs

# LOAD FUTURE PROJECTION DATA FOR DIFFERENT RCP SCENARIOS
rcp4p5 <- read.csv("future_projection_rcp45.csv")
rcp4p5$date <- as.Date(rcp4p5$date, format = "%d/%m/%Y") 
rcp8p5 <- read.csv("future_projection_rcp85.csv")
rcp8p5$date <- as.Date(rcp8p5$date, format = "%d/%m/%Y")
oimoy <- tapply(Lakki$total_cases,as.numeric(format(Lakki$date,"%m")),
                mean,na.rm=T)[seq(12)]
while(any(isna <- is.na(oimoy)))
  oimoy[isna] <- rowMeans(Lag(oimoy,c(-1,1)),na.rm=T)[isna]
oiproj <- rep(oimoy,length=nrow(rcp4p5))

# RUN CROSS-REDUCTION ON FUTURE DATA
red <- crossreduce(cb,model,at = 22.4:35.3,by=0.1)
coef <- coef(red)
vcov <- vcov(red)
oiperiod <- sum(oimm)*9

# SETUP FOR PROJECTION PERIOD
projperiod <- "2044-2052"
projseqperiod <- factor(rep(projperiod,length.out=12*length(seq(2044,2052))))
seqperiod <- factor(c(as.numeric(projseqperiod)))
levels(seqperiod) <- c(projperiod)
temprange <- c("tot")
absrel <- c("abs")

# SETUP RCP SCENARIOS FOR PROJECTION ANALYSIS
rcp <- c(RCP4.5="rcp4p5",RCP8.5="rcp8p5")
nsim <- 1000
ansim <- array(NA,dim=c(length(levels(seqperiod)),length(temprange),
                        length(absrel), length(rcp),nsim+1), 
               dimnames=list(levels(seqperiod),temprange,absrel,
                             names(rcp), c("est",paste0("sim",seq(nsim)))))

# RUN MONTE CARLO SIMULATIONS PER RCP SCENARIO
for (i in seq(rcp)) {
  
  # LOAD PROJECTED TEMPERATURE SERIES FOR CURRENT RCP SCENARIO
  cat("\n\n", names(rcp)[i], "\n")
  Tmeanproj <- get(rcp[[i]])
  argvarproj=list(fun="thr",thr=22.4)
  
  bvar <- do.call(onebasis,c(list(x=Tmeanproj$Tmean),argvarproj))
  cenvec <- do.call(onebasis,c(list(x=22.4),argvarproj))
  bvarcen <- scale(bvar, center = cenvec, scale = FALSE)
  indheat <- which(Lakki$Tmean >= 22.4 & Lakki$Tmean <= 35.3)
  an <- (1-exp(-bvarcen%*%coef(red)))*oiproj
  ansim[,"tot","abs",i,1] <- tapply(an[indheat],factor(seqperiod[indheat]),sum)
  
  # SAMPLE COEFFICIENTS AND RUN SIMULATIONS
  set.seed(13041975)
  coefsim <- mvrnorm(nsim,coef,vcov)
  for(s in seq(nsim)) {
    an <- (1-exp(-bvarcen%*%coefsim[s,]))*oi
    ansim[,"tot","abs",i,s+1] <- tapply(an[indheat],factor(seqperiod[indheat]),sum)
  }
}

# CALCULATE ATTRIBUTABLE NUMBERS AND FRACTIONS FOR RCP SCENARIOS
estci <- c("est","ci.l","ci.u")
anabs <- afabs <- array(NA,dim=c(length(levels(seqperiod)),
                                 length(estci),length(temprange),length(rcp)), 
                        dimnames=list(levels(seqperiod),estci,temprange,names(rcp)))
anabs[,"est",,"RCP4.5"] <- apply(ansim[,,"abs","RCP4.5",],1:1,mean)
anabs[,"ci.l",,"RCP4.5"] <- apply(ansim[,,"abs","RCP4.5",],1:1,quantile,0.025)
anabs[,"ci.u",,"RCP4.5"] <- apply(ansim[,,"abs","RCP4.5",],1:1,quantile,0.975)
anabs[,"est",,"RCP8.5"] <- apply(ansim[,,"abs","RCP8.5",],1:1,mean)
anabs[,"ci.l",,"RCP8.5"] <- apply(ansim[,,"abs","RCP8.5",],1:1,quantile,0.025)
anabs[,"ci.u",,"RCP8.5"] <- apply(ansim[,,"abs","RCP8.5",],1:1,quantile,0.975)
dim(ansim[, , "abs","RCP4.5", ])
afabs[,,,] <- anabs[,,,]/oiperiod*100
afabs
