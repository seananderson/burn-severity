#################################RdNBR
# rm(list=ls())
library(sampling)
# dataDir="D:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv("~/Downloads/ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171029.csv")

library(glmmTMB)
allData$Firemort.BA.p[allData$Firemort.BA.p == 1] <- 0.99
allData$Firemort.BA.p[allData$Firemort.BA.p == 0] <- 0.01
allData$RdNBR <- arm::rescale(allData$RdNBR)
m <- glmmTMB(Firemort.BA.p ~ RdNBR + (1 | FIRE),
  family=list(family="beta",link="logit"), data = allData)
summary(m)

library(rstanarm)
library(betareg)
options(mc.cores = parallel::detectCores())
m <- stan_betareg(Firemort.BA.p ~ RdNBR + as.factor(FIRE),
  data = allData, link = "logit", cores = 1, iter = 800, link.phi = "log")

library(mgcv)
m <- gamm(Firemort.BA.p ~ s(RdNBR) + (1 | FIRE),
  family = mgcv::betar(link = "logit"), data = allData)


names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-RdNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
# png(paste(dataDir,"fitPlots_RdNBR_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR,allData$Firemort.BA.p,pch=16,col="gray",xlab="RdNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR),max(allData$RdNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
# dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-RdNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR,allData$Firemort.trees.p,pch=16,col="gray",xlab="RdNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR),max(allData$RdNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-RdNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR,allData$CHARHT_percMax,pch=16,col="gray",xlab="RdNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR),max(allData$RdNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-RdNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR,allData$BOLESCORCH,pch=16,col="gray",xlab="RdNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR),max(allData$RdNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-RdNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR,allData$CHARCOV,pch=16,col="gray",xlab="RdNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR),max(allData$RdNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_CHARCOV.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
#################################RdNBR_BL

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR_BL
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-RdNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BL_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR_BL,allData$Firemort.BA.p,pch=16,col="gray",xlab="RdNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR_BL),max(allData$RdNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BL_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BL_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR_BL
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-RdNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BL_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR_BL,allData$Firemort.trees.p,pch=16,col="gray",xlab="RdNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR_BL),max(allData$RdNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BL_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BL_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR_BL
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-RdNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BL_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR_BL,allData$CHARHT_percMax,pch=16,col="gray",xlab="RdNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR_BL),max(allData$RdNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BL_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BL_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR_BL
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-RdNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BL_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR_BL,allData$BOLESCORCH,pch=16,col="gray",xlab="RdNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR_BL),max(allData$RdNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BL_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BL_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RdNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RdNBR_BL
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-RdNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RdNBR_BL_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RdNBR_BL,allData$CHARCOV,pch=16,col="gray",xlab="RdNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RdNBR_BL),max(allData$RdNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RdNBR_BL_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RdNBR_BL_CHARCOV.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
#######################dNBR

#################################dNBR
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-dNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR,allData$Firemort.BA.p,pch=16,col="gray",xlab="dNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR),max(allData$dNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-dNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR,allData$Firemort.trees.p,pch=16,col="gray",xlab="dNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR),max(allData$dNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-dNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR,allData$CHARHT_percMax,pch=16,col="gray",xlab="dNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR),max(allData$dNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-dNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR,allData$BOLESCORCH,pch=16,col="gray",xlab="dNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR),max(allData$dNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-dNBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR,allData$CHARCOV,pch=16,col="gray",xlab="dNBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR),max(allData$dNBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_CHARCOV.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
#################################dNBR_BL

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_BL
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-dNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BL_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_BL,allData$Firemort.BA.p,pch=16,col="gray",xlab="dNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_BL),max(allData$dNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BL_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BL_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_BL
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-dNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BL_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_BL,allData$Firemort.trees.p,pch=16,col="gray",xlab="dNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_BL),max(allData$dNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BL_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BL_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_BL
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-dNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BL_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_BL,allData$CHARHT_percMax,pch=16,col="gray",xlab="dNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_BL),max(allData$dNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BL_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BL_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_BL
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-dNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BL_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_BL,allData$BOLESCORCH,pch=16,col="gray",xlab="dNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_BL),max(allData$dNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BL_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BL_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_BL
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-dNBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_BL_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_BL,allData$CHARCOV,pch=16,col="gray",xlab="dNBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_BL),max(allData$dNBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_BL_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_BL_CHARCOV.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################
#################################dNBR_Offset
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-dNBR_Offset)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset,allData$Firemort.BA.p,pch=16,col="gray",xlab="dNBR_Offset",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset),max(allData$dNBR_Offset),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-dNBR_Offset)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset,allData$Firemort.trees.p,pch=16,col="gray",xlab="dNBR_Offset",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset),max(allData$dNBR_Offset),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-dNBR_Offset)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset,allData$CHARHT_percMax,pch=16,col="gray",xlab="dNBR_Offset",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset),max(allData$dNBR_Offset),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-dNBR_Offset)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset,allData$BOLESCORCH,pch=16,col="gray",xlab="dNBR_Offset",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset),max(allData$dNBR_Offset),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-dNBR_Offset)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset,allData$CHARCOV,pch=16,col="gray",xlab="dNBR_Offset",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset),max(allData$dNBR_Offset),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_CHARCOV.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
#################################dNBR_Offset_BL

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset_BL
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-dNBR_Offset_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BL_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset_BL,allData$Firemort.BA.p,pch=16,col="gray",xlab="dNBR_Offset_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset_BL),max(allData$dNBR_Offset_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BL_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset_BL
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-dNBR_Offset_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BL_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset_BL,allData$Firemort.trees.p,pch=16,col="gray",xlab="dNBR_Offset_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset_BL),max(allData$dNBR_Offset_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BL_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset_BL
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-dNBR_Offset_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BL_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset_BL,allData$CHARHT_percMax,pch=16,col="gray",xlab="dNBR_Offset_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset_BL),max(allData$dNBR_Offset_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BL_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset_BL
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-dNBR_Offset_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BL_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset_BL,allData$BOLESCORCH,pch=16,col="gray",xlab="dNBR_Offset_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset_BL),max(allData$dNBR_Offset_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BL_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$dNBR_Offset_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = dNBR_Offset_BL
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-dNBR_Offset_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_dNBR_Offset_BL_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$dNBR_Offset_BL,allData$CHARCOV,pch=16,col="gray",xlab="dNBR_Offset_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$dNBR_Offset_BL),max(allData$dNBR_Offset_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_dNBR_Offset_BL_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_CHARCOV.csv",sep=""),row.names=FALSE)

#################################################################################################
###################################################################################################
checkmofo
##########################################################################################
##########################################################################################
#######################RBR

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-RBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR,allData$Firemort.BA.p,pch=16,col="gray",xlab="RBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR),max(allData$RBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-RBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR,allData$Firemort.trees.p,pch=16,col="gray",xlab="RBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR),max(allData$RBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-RBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR,allData$CHARHT_percMax,pch=16,col="gray",xlab="RBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR),max(allData$RBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-RBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR,allData$BOLESCORCH,pch=16,col="gray",xlab="RBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR),max(allData$RBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-RBR)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR,allData$CHARCOV,pch=16,col="gray",xlab="RBR",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR),max(allData$RBR),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_CHARCOV.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
#################################RBR_BL

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR_BL
#Y = BA killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.BA.p),]
# SET UP MODEL
iForm=as.formula("Firemort.BA.p~1/(1+exp(a*(b-RBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BL_Firemort.BA.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR_BL,allData$Firemort.BA.p,pch=16,col="gray",xlab="RBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR_BL),max(allData$RBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.BA.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.BA.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.BA.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.BA.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BL_Firemort.BA.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BL_Firemort.BA.p.csv",sep=""),row.names=FALSE)

##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR_BL
#Y = trees killed by fire
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$Firemort.trees.p),]
# SET UP MODEL
iForm=as.formula("Firemort.trees.p~1/(1+exp(a*(b-RBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BL_Firemort.trees.p.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR_BL,allData$Firemort.trees.p,pch=16,col="gray",xlab="RBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR_BL),max(allData$RBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$Firemort.trees.p)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$Firemort.trees.p)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$Firemort.trees.p)**2))
  valRSQR=summary(lm(valFits~valData$Firemort.trees.p))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BL_Firemort.trees.p.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BL_Firemort.trees.p.csv",sep=""),row.names=FALSE)

##########################################################################################
##########################################################################################

rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR_BL
#Y = CHARHT_percMax
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARHT_percMax),]
# SET UP MODEL
iForm=as.formula("CHARHT_percMax~1/(1+exp(a*(b-RBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BL_CHARHT_percMax.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR_BL,allData$CHARHT_percMax,pch=16,col="gray",xlab="RBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR_BL),max(allData$RBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=850))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARHT_percMax)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARHT_percMax)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARHT_percMax)**2))
  valRSQR=summary(lm(valFits~valData$CHARHT_percMax))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BL_CHARHT_percMax.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BL_CHARHT_percMax.csv",sep=""),row.names=FALSE)



##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR_BL
#Y = BOLESCORCH
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$BOLESCORCH),]
# SET UP MODEL
iForm=as.formula("BOLESCORCH~1/(1+exp(a*(b-RBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BL_BOLESCORCH.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR_BL,allData$BOLESCORCH,pch=16,col="gray",xlab="RBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR_BL),max(allData$RBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.005,b=375))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$BOLESCORCH)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$BOLESCORCH)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$BOLESCORCH)**2))
  valRSQR=summary(lm(valFits~valData$BOLESCORCH))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BL_BOLESCORCH.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BL_BOLESCORCH.csv",sep=""),row.names=FALSE)


##########################################################################################
##########################################################################################
rm(list=ls())
library(sampling)
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

allData=read.csv(paste(dataDir,"ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171021.csv",sep=""))
names(allData)
dim(allData)
# SET DEFAULTS
nSims=1000  # Number of simulations
nCut=0.6	# Cal/Val stratification probability => 0.6 = 60/40
##### removing the plots that we do not have predictor value
allData = allData[complete.cases(allData$RBR_BL),]
# ...ORDER BY STRATIFYING VARIABLE
allData=allData[order(allData$FIRE),]
# ...CALCULATE STRATUM SIZES
stSizes=as.data.frame(table(allData$FIRE))
stSizes$SSIZE=round(stSizes$Freq*nCut)

#X = RBR_BL
#Y = CHARCOV
#pixel = focal
##### removing the plots that we do not have a response value
allData = allData[complete.cases(allData$CHARCOV),]
# SET UP MODEL
iForm=as.formula("CHARCOV~1/(1+exp(a*(b-RBR_BL)))")
# SET UP OUTPUT DATAFRAME
outCols=c("SIM","nlsA","nlsB","calRMSEnls","calRSQRnls","valRMSEnls","valRSQRnls")
outDF=matrix(data=NA,nrow=nSims,ncol=length(outCols))
outDF=as.data.frame(outDF)
names(outDF)=outCols
# START BASE PLOT FOR NLS FITS
png(paste(dataDir,"fitPlots_RBR_BL_CHARCOV.png",sep=""),height=5,width=5,units="in",res=200)
plot(allData$RBR_BL,allData$CHARCOV,pch=16,col="gray",xlab="RBR_BL",ylab="% Basal area killed by fire",xlim=c(-400,1500),ylim=c(0,1))
xArr=seq(min(allData$RBR_BL),max(allData$RBR_BL),length.out=100)
# START SIMULATIONS
for (sim in seq(nSims))
{
  print(c(sim,"of",nSims))
  flush.console()
  # ...copy dataframe
  intData=allData
  # ...get 60% stratum using a stratified random sample without replacement "srswor"
  s=strata(intData,stratanames="FIRE",size=stSizes$SSIZE,method="srswor")
  calData=getdata(intData,s)
  # ...get rest of data not used for calibration, use this for validation
  valData=intData[which(!(intData$KEY%in%calData$KEY)),]
  # NLS MODEL
  # Fit NLS model on calData
  resNLS=nls(iForm,data=calData,start=list(a=0.003,b=995))
  calRMSE=sqrt(mean((fitted(resNLS)-calData$CHARCOV)**2))
  calRSQR=summary(lm((fitted(resNLS)~calData$CHARCOV)))$adj.r.squared
  # Predict model on valData
  valFits=predict(resNLS,newdata=valData)
  valRMSE=sqrt(mean((valFits-valData$CHARCOV)**2))
  valRSQR=summary(lm(valFits~valData$CHARCOV))$adj.r.squared
  # Get coefficients, store all variables
  calCOEF=coef(resNLS)
  # STORE DATA
  outDF[sim,"SIM"]=sim
  # ...nls model
  outDF[sim,"nlsA"]=calCOEF[1]
  outDF[sim,"nlsB"]=calCOEF[2]
  outDF[sim,"calRMSEnls"]=calRMSE
  outDF[sim,"calRSQRnls"]=calRSQR
  outDF[sim,"valRMSEnls"]=valRMSE
  outDF[sim,"valRSQRnls"]=valRSQR
  # Add lines to plot
  yArr=1/(1+exp(calCOEF[1]*(calCOEF[2]-xArr)))
  lines(xArr,yArr,col=sim)
}
# Close NLS plot
dev.off()

# Make plots of R2 and RMSE histograms comparing calibration and validation
png(paste(dataDir,"fitHists_RBR_BL_CHARCOV.png",sep=""),height=5,width=7,units="in",res=200)
par(mfcol=c(2,2))
hist(outDF$calRSQRnls,main="Calibration",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$valRSQRnls,main="Validation",xlab="R2",xlim=c(0.4,0.8),breaks=20)
hist(outDF$calRMSEnls,main="Calibration",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
hist(outDF$valRMSEnls,main="Validation",xlab="RMSE",xlim=c(0.15,0.35),breaks=20)
dev.off()
# WRITE OUT MODEL FIT DATA
write.csv(outDF,paste(dataDir,"output_CalVal_Stats_RBR_BL_CHARCOV.csv",sep=""),row.names=FALSE)

############################################
############################################
############################################

#Bring in R2 and RMSE tables to combine
rm(list=ls())
dataDir="E:\\Dropbox\\BJH_ComputerFiles\\Work\\BriansRAnalysis\\data\\"

# READ IN DATA, CHECK

## Firemort.BA.p
#RdNBR
Rd_BA=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_Firemort.BA.p.csv",sep=""))
names(Rd_BA)
dim(Rd_BA)
Rd_BA$index=ifelse(Rd_BA$SIM>0,"RdNBR","9999")
Rd_BA$variable=ifelse(Rd_BA$SIM>0,"Firemort.BA.p","9999")
Rd_BA$pixel=ifelse(Rd_BA$SIM>0,"focal","9999")

#RdNBR_BL
Rd_BL_BA=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BL_Firemort.BA.p.csv",sep=""))
names(Rd_BL_BA)
dim(Rd_BL_BA)
Rd_BL_BA$index=ifelse(Rd_BL_BA$SIM>0,"RdNBR_BL","9999")
Rd_BL_BA$variable=ifelse(Rd_BL_BA$SIM>0,"Firemort.BA.p","9999")
Rd_BL_BA$pixel=ifelse(Rd_BL_BA$SIM>0,"BL","9999")

#dNBR
d_BA=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Firemort.BA.p.csv",sep=""))
names(d_BA)
dim(d_BA)
d_BA$index=ifelse(d_BA$SIM>0,"dNBR","9999")
d_BA$variable=ifelse(d_BA$SIM>0,"Firemort.BA.p","9999")
d_BA$pixel=ifelse(d_BA$SIM>0,"focal","9999")

#dNBR_BL
d_BL_BA=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BL_Firemort.BA.p.csv",sep=""))
names(d_BL_BA)
dim(d_BL_BA)
d_BL_BA$index=ifelse(d_BL_BA$SIM>0,"dNBR_BL","9999")
d_BL_BA$variable=ifelse(d_BL_BA$SIM>0,"Firemort.BA.p","9999")
d_BL_BA$pixel=ifelse(d_BL_BA$SIM>0,"BL","9999")


#dNBR_Offset
d_Offset_BA=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_Firemort.BA.p.csv",sep=""))
names(d_Offset_BA)
dim(d_Offset_BA)
d_Offset_BA$index=ifelse(d_Offset_BA$SIM>0,"dNBR_Offset","9999")
d_Offset_BA$variable=ifelse(d_Offset_BA$SIM>0,"Firemort.BA.p","9999")
d_Offset_BA$pixel=ifelse(d_Offset_BA$SIM>0,"focal","9999")

#dNBR_Offset_BL
d_Offset_BL_BA=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_Firemort.BA.p.csv",sep=""))
names(d_Offset_BL_BA)
dim(d_Offset_BL_BA)
d_Offset_BL_BA$index=ifelse(d_Offset_BL_BA$SIM>0,"dNBR_Offset_BL","9999")
d_Offset_BL_BA$variable=ifelse(d_Offset_BL_BA$SIM>0,"Firemort.BA.p","9999")
d_Offset_BL_BA$pixel=ifelse(d_Offset_BL_BA$SIM>0,"BL","9999")

#RBR
RBR_BA=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_Firemort.BA.p.csv",sep=""))
names(RBR_BA)
dim(RBR_BA)
RBR_BA$index=ifelse(RBR_BA$SIM>0,"RBR","9999")
RBR_BA$variable=ifelse(RBR_BA$SIM>0,"Firemort.BA.p","9999")
RBR_BA$pixel=ifelse(RBR_BA$SIM>0,"focal","9999")

#RBR_BL
RBR_BL_BA=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BL_Firemort.BA.p.csv",sep=""))
names(RBR_BL_BA)
dim(RBR_BL_BA)
RBR_BL_BA$index=ifelse(RBR_BL_BA$SIM>0,"RBR_BL","9999")
RBR_BL_BA$variable=ifelse(RBR_BL_BA$SIM>0,"Firemort.BA.p","9999")
RBR_BL_BA$pixel=ifelse(RBR_BL_BA$SIM>0,"BL","9999")

## Firemort.trees.p
#RdNBR
Rd_trees=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_Firemort.trees.p.csv",sep=""))
names(Rd_trees)
dim(Rd_trees)
Rd_trees$index=ifelse(Rd_trees$SIM>0,"RdNBR","9999")
Rd_trees$variable=ifelse(Rd_trees$SIM>0,"Firemort.trees.p","9999")
Rd_trees$pixel=ifelse(Rd_trees$SIM>0,"focal","9999")

#RdNBR_BL
Rd_BL_trees=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BL_Firemort.trees.p.csv",sep=""))
names(Rd_BL_trees)
dim(Rd_BL_trees)
Rd_BL_trees$index=ifelse(Rd_trees$SIM>0,"RdNBR_BL","9999")
Rd_BL_trees$variable=ifelse(Rd_trees$SIM>0,"Firemort.trees.p","9999")
Rd_BL_trees$pixel=ifelse(Rd_BL_trees$SIM>0,"BL","9999")

#dNBR
d_trees=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Firemort.trees.p.csv",sep=""))
names(d_trees)
dim(d_trees)
d_trees$index=ifelse(d_trees$SIM>0,"dNBR","9999")
d_trees$variable=ifelse(d_trees$SIM>0,"Firemort.trees.p","9999")
d_trees$pixel=ifelse(d_trees$SIM>0,"focal","9999")

#dNBR_BL
d_BL_trees=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BL_Firemort.trees.p.csv",sep=""))
names(d_BL_trees)
dim(d_BL_trees)
d_BL_trees$index=ifelse(d_BL_trees$SIM>0,"dNBR_BL","9999")
d_BL_trees$variable=ifelse(d_BL_trees$SIM>0,"Firemort.trees.p","9999")
d_BL_trees$pixel=ifelse(d_BL_trees$SIM>0,"BL","9999")

#dNBR_Offset
d_Offset_trees=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_Firemort.trees.p.csv",sep=""))
names(d_Offset_trees)
dim(d_Offset_trees)
d_Offset_trees$index=ifelse(d_Offset_trees$SIM>0,"dNBR_Offset","9999")
d_Offset_trees$variable=ifelse(d_Offset_trees$SIM>0,"Firemort.trees.p","9999")
d_Offset_trees$pixel=ifelse(d_Offset_trees$SIM>0,"focal","9999")

#dNBR_Offset_BL
d_Offset_BL_trees=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_Firemort.trees.p.csv",sep=""))
names(d_Offset_BL_trees)
dim(d_Offset_BL_trees)
d_Offset_BL_trees$index=ifelse(d_Offset_BL_trees$SIM>0,"dNBR_Offset_BL","9999")
d_Offset_BL_trees$variable=ifelse(d_Offset_BL_trees$SIM>0,"Firemort.trees.p","9999")
d_Offset_BL_trees$pixel=ifelse(d_Offset_BL_trees$SIM>0,"BL","9999")

#RBR
RBR_trees=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_Firemort.trees.p.csv",sep=""))
names(RBR_trees)
dim(RBR_trees)
RBR_trees$index=ifelse(RBR_trees$SIM>0,"RBR","9999")
RBR_trees$variable=ifelse(RBR_trees$SIM>0,"Firemort.trees.p","9999")
RBR_trees$pixel=ifelse(RBR_trees$SIM>0,"focal","9999")

#RBR_BL
RBR_BL_trees=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BL_Firemort.trees.p.csv",sep=""))
names(RBR_BL_trees)
dim(RBR_BL_trees)
RBR_BL_trees$index=ifelse(RBR_BL_trees$SIM>0,"RBR_BL","9999")
RBR_BL_trees$variable=ifelse(RBR_BL_trees$SIM>0,"Firemort.trees.p","9999")
RBR_BL_trees$pixel=ifelse(RBR_BL_trees$SIM>0,"BL","9999")

## CHARHT_percmax
#RdNBR
Rd_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_CHARHT_percMax.csv",sep=""))
names(Rd_CHARHT)
dim(Rd_CHARHT)
Rd_CHARHT$index=ifelse(Rd_trees$SIM>0,"RdNBR","9999")
Rd_CHARHT$variable=ifelse(Rd_trees$SIM>0,"CHARHT","9999")
Rd_CHARHT$pixel=ifelse(Rd_trees$SIM>0,"focal","9999")

#RdNBR_BL
Rd_BL_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BL_CHARHT_percMax.csv",sep=""))
names(Rd_BL_CHARHT)
dim(Rd_BL_CHARHT)
Rd_BL_CHARHT$index=ifelse(Rd_BL_CHARHT$SIM>0,"RdNBR_BL","9999")
Rd_BL_CHARHT$variable=ifelse(Rd_BL_CHARHT$SIM>0,"CHARHT","9999")
Rd_BL_CHARHT$pixel=ifelse(Rd_BL_CHARHT$SIM>0,"BL","9999")

#dNBR
d_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_CHARHT_percMax.csv",sep=""))
names(d_CHARHT)
dim(d_CHARHT)
d_CHARHT$index=ifelse(d_CHARHT$SIM>0,"dNBR","9999")
d_CHARHT$variable=ifelse(d_CHARHT$SIM>0,"CHARHT","9999")
d_CHARHT$pixel=ifelse(d_CHARHT$SIM>0,"focal","9999")

#dNBR_BL
d_BL_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BL_CHARHT_percMax.csv",sep=""))
names(d_BL_CHARHT)
dim(d_BL_CHARHT)
d_BL_CHARHT$index=ifelse(d_BL_CHARHT$SIM>0,"dNBR_BL","9999")
d_BL_CHARHT$variable=ifelse(d_BL_CHARHT$SIM>0,"CHARHT","9999")
d_BL_CHARHT$pixel=ifelse(d_BL_CHARHT$SIM>0,"BL","9999")

#dNBR_Offset
d_Offset_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_CHARHT_percMax.csv",sep=""))
names(d_Offset_CHARHT)
dim(d_Offset_CHARHT)
d_Offset_CHARHT$index=ifelse(d_Offset_CHARHT$SIM>0,"dNBR_Offset","9999")
d_Offset_CHARHT$variable=ifelse(d_Offset_CHARHT$SIM>0,"CHARHT","9999")
d_Offset_CHARHT$pixel=ifelse(d_Offset_CHARHT$SIM>0,"focal","9999")

#dNBR_Offset_BL
d_Offset_BL_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_CHARHT_percMax.csv",sep=""))
names(d_Offset_BL_CHARHT)
dim(d_Offset_BL_CHARHT)
d_Offset_BL_CHARHT$index=ifelse(d_Offset_BL_CHARHT$SIM>0,"dNBR_Offset_BL","9999")
d_Offset_BL_CHARHT$variable=ifelse(d_Offset_BL_CHARHT$SIM>0,"CHARHT","9999")
d_Offset_BL_CHARHT$pixel=ifelse(d_Offset_BL_CHARHT$SIM>0,"BL","9999")

#RBR
RBR_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_CHARHT_percMax.csv",sep=""))
names(RBR_CHARHT)
dim(RBR_CHARHT)
RBR_CHARHT$index=ifelse(RBR_CHARHT$SIM>0,"RBR","9999")
RBR_CHARHT$variable=ifelse(RBR_CHARHT$SIM>0,"CHARHT","9999")
RBR_CHARHT$pixel=ifelse(RBR_CHARHT$SIM>0,"focal","9999")

#RBR_BL
RBR_BL_CHARHT=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BL_CHARHT_percMax.csv",sep=""))
names(RBR_BL_CHARHT)
dim(RBR_BL_CHARHT)
RBR_BL_CHARHT$index=ifelse(RBR_BL_CHARHT$SIM>0,"RBR_BL","9999")
RBR_BL_CHARHT$variable=ifelse(RBR_BL_CHARHT$SIM>0,"CHARHT","9999")
RBR_BL_CHARHT$pixel=ifelse(RBR_BL_CHARHT$SIM>0,"BL","9999")

## BOLESCORCH
#RdNBR
Rd_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BOLESCORCH.csv",sep=""))
names(Rd_BOLESCORCH)
dim(Rd_BOLESCORCH)
Rd_BOLESCORCH$index=ifelse(Rd_trees$SIM>0,"RdNBR","9999")
Rd_BOLESCORCH$variable=ifelse(Rd_trees$SIM>0,"BOLESCORCH","9999")
Rd_BOLESCORCH$pixel=ifelse(Rd_trees$SIM>0,"focal","9999")

#RdNBR_BL
Rd_BL_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BL_BOLESCORCH.csv",sep=""))
names(Rd_BL_BOLESCORCH)
dim(Rd_BL_BOLESCORCH)
Rd_BL_BOLESCORCH$index=ifelse(Rd_BL_BOLESCORCH$SIM>0,"RdNBR_BL","9999")
Rd_BL_BOLESCORCH$variable=ifelse(Rd_BL_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
Rd_BL_BOLESCORCH$pixel=ifelse(Rd_BL_BOLESCORCH$SIM>0,"BL","9999")

#dNBR
d_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BOLESCORCH.csv",sep=""))
names(d_BOLESCORCH)
dim(d_BOLESCORCH)
d_BOLESCORCH$index=ifelse(d_BOLESCORCH$SIM>0,"dNBR","9999")
d_BOLESCORCH$variable=ifelse(d_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
d_BOLESCORCH$pixel=ifelse(d_BOLESCORCH$SIM>0,"focal","9999")

#dNBR_BL
d_BL_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BL_BOLESCORCH.csv",sep=""))
names(d_BL_BOLESCORCH)
dim(d_BL_BOLESCORCH)
d_BL_BOLESCORCH$index=ifelse(d_BL_BOLESCORCH$SIM>0,"dNBR_BL","9999")
d_BL_BOLESCORCH$variable=ifelse(d_BL_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
d_BL_BOLESCORCH$pixel=ifelse(d_BL_BOLESCORCH$SIM>0,"BL","9999")

#dNBR_Offset
d_Offset_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BOLESCORCH.csv",sep=""))
names(d_Offset_BOLESCORCH)
dim(d_Offset_BOLESCORCH)
d_Offset_BOLESCORCH$index=ifelse(d_Offset_BOLESCORCH$SIM>0,"dNBR_Offset","9999")
d_Offset_BOLESCORCH$variable=ifelse(d_Offset_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
d_Offset_BOLESCORCH$pixel=ifelse(d_Offset_BOLESCORCH$SIM>0,"focal","9999")

#dNBR_Offset_BL
d_Offset_BL_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_BOLESCORCH.csv",sep=""))
names(d_Offset_BL_BOLESCORCH)
dim(d_Offset_BL_BOLESCORCH)
d_Offset_BL_BOLESCORCH$index=ifelse(d_Offset_BL_BOLESCORCH$SIM>0,"dNBR_Offset_BL","9999")
d_Offset_BL_BOLESCORCH$variable=ifelse(d_Offset_BL_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
d_Offset_BL_BOLESCORCH$pixel=ifelse(d_Offset_BL_BOLESCORCH$SIM>0,"BL","9999")

#RBR
RBR_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BOLESCORCH.csv",sep=""))
names(RBR_BOLESCORCH)
dim(RBR_BOLESCORCH)
RBR_BOLESCORCH$index=ifelse(RBR_BOLESCORCH$SIM>0,"RBR","9999")
RBR_BOLESCORCH$variable=ifelse(RBR_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
RBR_BOLESCORCH$pixel=ifelse(RBR_BOLESCORCH$SIM>0,"focal","9999")

#RBR_BL
RBR_BL_BOLESCORCH=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BL_BOLESCORCH.csv",sep=""))
names(RBR_BL_BOLESCORCH)
dim(RBR_BL_BOLESCORCH)
RBR_BL_BOLESCORCH$index=ifelse(RBR_BL_BOLESCORCH$SIM>0,"RBR_BL","9999")
RBR_BL_BOLESCORCH$variable=ifelse(RBR_BL_BOLESCORCH$SIM>0,"BOLESCORCH","9999")
RBR_BL_BOLESCORCH$pixel=ifelse(RBR_BL_BOLESCORCH$SIM>0,"BL","9999")

## CHARCOV
#RdNBR
Rd_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_CHARCOV.csv",sep=""))
names(Rd_CHARCOV)
dim(Rd_CHARCOV)
Rd_CHARCOV$index=ifelse(Rd_trees$SIM>0,"RdNBR","9999")
Rd_CHARCOV$variable=ifelse(Rd_trees$SIM>0,"CHARCOV","9999")
Rd_CHARCOV$pixel=ifelse(Rd_trees$SIM>0,"focal","9999")

#RdNBR_BL
Rd_BL_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_RdNBR_BL_CHARCOV.csv",sep=""))
names(Rd_BL_CHARCOV)
dim(Rd_BL_CHARCOV)
Rd_BL_CHARCOV$index=ifelse(Rd_BL_CHARCOV$SIM>0,"RdNBR_BL","9999")
Rd_BL_CHARCOV$variable=ifelse(Rd_BL_CHARCOV$SIM>0,"CHARCOV","9999")
Rd_BL_CHARCOV$pixel=ifelse(Rd_BL_CHARCOV$SIM>0,"BL","9999")

#dNBR
d_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_CHARCOV.csv",sep=""))
names(d_CHARCOV)
dim(d_CHARCOV)
d_CHARCOV$index=ifelse(d_CHARCOV$SIM>0,"dNBR","9999")
d_CHARCOV$variable=ifelse(d_CHARCOV$SIM>0,"CHARCOV","9999")
d_CHARCOV$pixel=ifelse(d_CHARCOV$SIM>0,"focal","9999")

#dNBR_BL
d_BL_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_BL_CHARCOV.csv",sep=""))
names(d_BL_CHARCOV)
dim(d_BL_CHARCOV)
d_BL_CHARCOV$index=ifelse(d_BL_CHARCOV$SIM>0,"dNBR_BL","9999")
d_BL_CHARCOV$variable=ifelse(d_BL_CHARCOV$SIM>0,"CHARCOV","9999")
d_BL_CHARCOV$pixel=ifelse(d_BL_CHARCOV$SIM>0,"BL","9999")

#dNBR_Offset
d_Offset_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_CHARCOV.csv",sep=""))
names(d_Offset_CHARCOV)
dim(d_Offset_CHARCOV)
d_Offset_CHARCOV$index=ifelse(d_Offset_CHARCOV$SIM>0,"dNBR_Offset","9999")
d_Offset_CHARCOV$variable=ifelse(d_Offset_CHARCOV$SIM>0,"CHARCOV","9999")
d_Offset_CHARCOV$pixel=ifelse(d_Offset_CHARCOV$SIM>0,"focal","9999")

#dNBR_Offset_BL
d_Offset_BL_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_dNBR_Offset_BL_CHARCOV.csv",sep=""))
names(d_Offset_BL_CHARCOV)
dim(d_Offset_BL_CHARCOV)
d_Offset_BL_CHARCOV$index=ifelse(d_Offset_BL_CHARCOV$SIM>0,"dNBR_Offset_BL","9999")
d_Offset_BL_CHARCOV$variable=ifelse(d_Offset_BL_CHARCOV$SIM>0,"CHARCOV","9999")
d_Offset_BL_CHARCOV$pixel=ifelse(d_Offset_BL_CHARCOV$SIM>0,"BL","9999")

#RBR
RBR_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_CHARCOV.csv",sep=""))
names(RBR_CHARCOV)
dim(RBR_CHARCOV)
RBR_CHARCOV$index=ifelse(RBR_CHARCOV$SIM>0,"RBR","9999")
RBR_CHARCOV$variable=ifelse(RBR_CHARCOV$SIM>0,"CHARCOV","9999")
RBR_CHARCOV$pixel=ifelse(RBR_CHARCOV$SIM>0,"focal","9999")

#RBR_BL
RBR_BL_CHARCOV=read.csv(paste(dataDir,"output_CalVal_Stats_RBR_BL_CHARCOV.csv",sep=""))
names(RBR_BL_CHARCOV)
dim(RBR_BL_CHARCOV)
RBR_BL_CHARCOV$index=ifelse(RBR_BL_CHARCOV$SIM>0,"RBR_BL","9999")
RBR_BL_CHARCOV$variable=ifelse(RBR_BL_CHARCOV$SIM>0,"CHARCOV","9999")
RBR_BL_CHARCOV$pixel=ifelse(RBR_BL_CHARCOV$SIM>0,"BL","9999")


AllModels = rbind(d_BA,d_trees,d_CHARHT,d_BOLESCORCH,d_CHARCOV,d_BL_BA,d_BL_trees,
                  d_BL_CHARHT,d_BL_BOLESCORCH,d_BL_CHARCOV,Rd_BA,Rd_trees,Rd_CHARHT,
                  Rd_BOLESCORCH,Rd_CHARCOV,Rd_BL_BA,Rd_BL_trees,Rd_BL_CHARHT,Rd_BL_BOLESCORCH,
                  Rd_BL_CHARCOV,d_Offset_BA,d_Offset_trees,d_Offset_CHARHT,d_Offset_BOLESCORCH,
                  d_Offset_CHARCOV,d_Offset_BL_BA,d_Offset_BL_trees,d_Offset_BL_CHARHT,
                  d_Offset_BL_BOLESCORCH,d_Offset_BL_CHARCOV,RBR_BA,RBR_trees,RBR_CHARHT,
                  RBR_BOLESCORCH,RBR_CHARCOV,RBR_BL_BA,RBR_BL_trees,RBR_BL_CHARHT,
                  RBR_BL_BOLESCORCH,RBR_BL_CHARCOV)


#### Making box plots of R2 and RMSE
library(ggplot2)
BA = ggplot(AllModels$RBR_BA)
BA + geom_boxplot()

ggplot((subset(AllModels,variable=="Firemort.BA.p")),aes(x=index, y=valRSQRnls))
geom_boxplot()

boxplot(valRMSEnls~pixel,data=AllModels,ylab="RMSE")
boxplot(valRSQRnls~pixel,data=AllModels,ylab="R^2")

AllModels = subset(AllModels, pixel =="focal")

boxplot(valRSQRnls~index,data=subset(AllModels,variable=="Firemort.BA.p"),ylab="R^2")
boxplot(valRSQRnls~index,data=subset(AllModels,variable=="Firemort.trees.p"),ylab="R^2")
boxplot(valRSQRnls~index,data=subset(AllModels,variable=="CHARHT"),ylab="R^2")
boxplot(valRSQRnls~index,data=subset(AllModels,variable=="BOLESCORCH"),ylab="R^2")
boxplot(valRSQRnls~index,data=subset(AllModels,variable=="CHARCOV"),ylab="R^2")

boxplot(valRMSEnls~index,data=subset(AllModels,variable=="Firemort.BA.p"),ylab="RMSE")
boxplot(valRMSEnls~index,data=subset(AllModels,variable=="Firemort.trees.p"),ylab="RMSE")
boxplot(valRMSEnls~index,data=subset(AllModels,variable=="CHARHT"),ylab="RMSE")
boxplot(valRMSEnls~index,data=subset(AllModels,variable=="BOLESCORCH"),ylab="RMSE")
boxplot(valRMSEnls~index,data=subset(AllModels,variable=="CHARCOV"),ylab="RMSE")
