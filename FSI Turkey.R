library(zoo)
library(dlm)
library(fGarch)
library(xlsx)
library(rmgarch)
library(ggplot2)
library(xts)
library(PerformanceAnalytics)

#Banking Sector

#Banking Dynamic Beta

turkeybdcc<-read.csv(file="TurkeyBankDCC.csv",header=TRUE)
turkeybdcc.zoo<-zoo(turkeybdcc[,-1], order.by=as.Date(strptime(as.character(turkeybdcc[,1]), "%m/%d/%Y")))
turkeybanka.zoo <- turkeybdcc.zoo[, "turkeybanka"]


# 1 Year yield for the risk free rate

turkeybdcc_1y<-read.csv(file="TurkeyBankDCC-3m1y2y.csv",header=TRUE)
turkeybdcc_1y.zoo<-zoo(turkeybdcc_1y[,-1], order.by=as.Date(strptime(as.character(turkeybdcc_1y[,1]), "%m/%d/%Y")))
tr1yyc<-read.csv(file="tr1yyield.csv",header=TRUE)
tr1yyc.zoo<-zoo(tr1yyc[,-1], order.by=as.Date(strptime(as.character(tr1yyc[,1]), "%m/%d/%Y")))
turkeybdcc1yyield<-merge(turkeybdcc_1y.zoo,tr1yyc.zoo)
turkeybdcc1yyield[,3]<-na.approx(na.trim(turkeybdcc1yyield[,3],side="both"))

# 3 month yield for the risk free rate

turkeybdcc_3m<-read.csv(file="TurkeyBankDCC-3m1y2y.csv",header=TRUE)
turkeybdcc_3m.zoo<-zoo(turkeybdcc_3m[,-1], order.by=as.Date(strptime(as.character(turkeybdcc_3m[,1]), "%m/%d/%Y")))
tr3myc<-read.csv(file="tr3monthyield.csv",header=TRUE)
tr3myc.zoo<-zoo(tr3myc[,-1], order.by=as.Date(strptime(as.character(tr3myc[,1]), "%m/%d/%Y")))
turkeybdcc3myield<-merge(turkeybdcc_3m.zoo,tr3myc.zoo)
turkeybdcc3myield[,3]<-na.approx(na.trim(turkeybdcc3myield[,3],side="both"))


# 2 year yield for the risk free rate

turkeybdcc_2y<-read.csv(file="TurkeyBankDCC-3m1y2y.csv",header=TRUE)
turkeybdcc_2y.zoo<-zoo(turkeybdcc_2y[,-1], order.by=as.Date(strptime(as.character(turkeybdcc_2y[,1]), "%m/%d/%Y")))
tr2yyc<-read.csv(file="tr2yyield.csv",header=TRUE)
tr2yyc.zoo<-zoo(tr2yyc[,-1], order.by=as.Date(strptime(as.character(tr2yyc[,1]), "%m/%d/%Y")))
turkeybdcc2yyield<-merge(turkeybdcc_2y.zoo,tr2yyc.zoo)
turkeybdcc2yyield[,3]<-na.approx(na.trim(turkeybdcc2yyield[,3],side="both"))

##Function
## x = banking index
## y = market index
## z = risk free rate
dynamic_beta<-function(xx){
  y<-diff(log(xx[,1])) ##Logged returns of Banking Sector Index 
  x<-diff(log(xx[,2])) #Logged returns of Stock Market Index 
  z<- exp(xx[,3][2:length(xx[,3])]*diff(index(xx[,3][2:length(xx[,3])]))/365)-1
  x<-x-z 
  y<-y-z
  garch.spec<-ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "norm")
  dcc.spec<-dccspec(uspec = multispec( replicate(2, garch.spec) ), dccOrder = c(1,1), distribution = "mvnorm")
  dcc.fit<-dccfit(dcc.spec, data = merge(x,y), fit.control=list(scale=TRUE))
  ## Instead use the info given notes on DCC ... Read the article
  ## The intertemporal capital asset pricing model with dynamic conditional correlation BaliEngleJME
  ## anyway the old one should be like..
  rcv=rcov(dcc.fit)
  beta <- zoo(rcv[1,2,] / rcv[2,2,], order.by=time(x))
  return(beta)
}

dccturkeybeta<-dynamic_beta(turkeybdcc1yyield)
dccturkeybeta_3m<-dynamic_beta(turkeybdcc3myield)
dccturkeybeta_2y<-dynamic_beta(turkeybdcc2yyield)


# Plot Dynamic Betas

plot(dccturkeybeta,xlab="Date", ylab="Dynamic Betas", main="Bank Index Dynamic Betas-1 y. yield")
# png("Dynamic beta 3 month.png", width=8, height=6, units="in", res=600)
plot(dccturkeybeta_3m,xlab="Date", ylab="Dynamic Betas", main="Bank Index Dynamic Betas-3 m. yield")
# dev.off()

# png("Dynamic beta 2 year.png", width=8, height=6, units="in", res=600)
plot(dccturkeybeta_2y,xlab="Date", ylab="Dynamic Betas", main="Bank Index Dynamic Betas-2 y. yield")
# dev.off()

# Dynamic Correlations 

dynamic_cor<-function(xx){
  y<-diff(log(xx[,1])) ##Logged returns of Banking Sector Index 
  x<-diff(log(xx[,2])) #Logged returns of Stock Market Index 
  z <-exp(xx[,3][2:length(xx[,3])]*diff(index(xx[,3][2:length(xx[,3])]))/365)-1
  x<-x-z 
  y<-y-z
  garch.spec<-ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "norm")
  dcc.spec<-dccspec(uspec = multispec( replicate(2, garch.spec) ), dccOrder = c(1,1), distribution = "mvnorm")
  dcc.fit<-dccfit(dcc.spec, data = merge(x,y), fit.control=list(scale=TRUE))
  rcr=zoo(rcor(dcc.fit)[1,2,], order.by=time(x))
  return(rcr)
}


dyncortur<-dynamic_cor(turkeybdcc1yyield)

# Plot Dynamic Correlations

plot(dyncortur)

# Event Analysis


dccturkeybeta_fsi<-window(xts(dccturkeybeta), start="2005-01-11", end="2016-11-17")
dccturkeybeta_fsi_scale<-scale(dccturkeybeta_fsi)
colnames(dccturkeybeta_fsi_scale)<-"dynamic betas"

#Cycles

cycles.dates<-c("2006-05/2006-07",
                "2007-08/2009-12",
                "2013-06/2014-03",
                "2014-12/2015-04")

# Event lists 

risk.dates = c(
  "2007-08-09",
  "2008-09-15",
  "2008-11-12",
  "2011-11-09",
  "2013-05-22",
  "2013-12-18",
  "2014-10-29",
  "2015-12-16",
  "2016-07-15",
  "2016-11-08"
)


risk.labels = c(
  "BNP Paribas",
  "Lehman Brother",
  "TARP Release",
  "European Debt Crisis",
  "1.st Fed Tapering",
  "2.nd Fed Tapering",
  "End of Fed Asset Purchasing Program",
  "Fed Interest Rate Hike and Abusi.Invest.",
  "Coup Attempt",
  "U.S. Presidential Election"
)

chart.TimeSeries(dccturkeybeta_fsi_scale, colorset = "black", legend.loc = "bottomleft", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Dynamic betas",
                 xlab = "Date", main="Dynamic Betas of Banking Index")




# Cmax of Bank Index

cmax_bank<-rollapply(xts(turkeybanka.zoo), align="right", width = 504, FUN = function(x) 1-(xts::last(x)/max(x)))
head(cmax_bank)
tail(cmax_bank)

# Common period

cmax_bank_fsi<-window(cmax_bank, start="2005-01-11", end="2016-11-17")
cmax_bank_fsi_scale<-scale(cmax_bank_fsi)
colnames(cmax_bank_fsi_scale)<-"cmax of bank"

# Event Analysis to save it delete

# png("Banking CMAX.png", width=8, height=6, units="in", res=600)

chart.TimeSeries(cmax_bank_fsi_scale, colorset = "black", legend.loc = "topleft", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "CMAX", 
                 xlab = "Date", main="Banking Index CMAX")

# dev.off()




# Idiosyncratic risk of bank stock prices
# GARCH(1,1) model

dturkeybanka.zoo<-diff(log(turkeybanka.zoo))
acf(as.vector(dturkeybanka.zoo))
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), distribution.model = "norm")

garch.fit_dturkeybanka = ugarchfit(garch11.spec, data = dturkeybanka.zoo, fit.control=list(scale=TRUE))
print(garch.fit_dturkeybanka)
head(sigma(garch.fit_dturkeybanka))
tail(sigma(garch.fit_dturkeybanka))

dturkeybanka_idifsi<-window(sigma(garch.fit_dturkeybanka), start="2005-01-11", end="2016-11-17")
dturkeybanka_idifsi_scale<-scale(dturkeybanka_idifsi)
colnames(dturkeybanka_idifsi_scale)<-"Id. Vol. of Banking Index"

#Event Analysis

# png("Id. Vol. of Banking Index.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(dturkeybanka_idifsi_scale, colorset = "black", legend.loc = "topleft", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Index", 
                 xlab = "Date", main="Banking Index Id. Vol.")
# dev.off()



bank_fsi<-merge(cmax_bank_fsi_scale,dturkeybanka_idifsi_scale,dccturkeybeta_fsi_scale )

head(bank_fsi)
tail(bank_fsi)


# Bank Index

bankfsi<-zoo(bank_fsi[,1])

for(i in 1:length(bankfsi)){
  bankfsi[i]=1/3*bank_fsi[i,1]+1/3*bank_fsi[i,2]+
    1/3*bank_fsi[i,3]
}

bankfsi<-xts(bankfsi)
colnames(bankfsi)<-"bank index"

# Event Analysis


# png("Banking Index_11.29.16.png", width=8, height=6, units="in", res=600)

chart.TimeSeries(bankfsi, colorset = "black", legend.loc = "topleft", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Index", 
                 xlab = "Date", main="Banking Sector Index")

# dev.off()

# Plot

png("Banking Index_plot_11.23.16.png", width=8, height=6, units="in", res=600)

plot(bankfsi, ylab = "Index",  xlab = "Date", main="Banking Sector Index for Turkey")

dev.off()

## Bond Market

## Slope-Curvature

library(readxl)
library(zoo)
library(chron)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xts)
library(quantmod)
library(lubridate)
library(xlsx)

yc<-read.csv(file = "yc.csv", header = TRUE)

head(yc)
tail(yc)
str(yc)
yc.z <- zoo(yc[,-1], order.by=as.Date(strptime(as.character(yc[,1]), "%m/%d/%Y")))

plot(yc.z, main="Yield Curve Data", xlab="Date")

remove_bad_data <- function(dta){
  ## Tricky one.. There must be an iteration.. if there exist no bad data it must stop
  ## Assume that the data is in zoo format.. The assumption is
  ## there should not be jumps more than +- 3 standard deviation
  ## But clearing some of them may make others bad data.. Hence
  ## lets do it iteratively.
  
  ## Do not forget that the best way is to use outlier detection
  ## procedures. If you have time please do this step
  ## by using some of the procedures we used in class
  
  ## Divide the procedure for univarita and multivariate case
  if (NCOL(dta) == 1){
    # add some values if negative values exist.. 
    dum <- min(dta, na.rm=T)
    ifelse(dum < 0, dum <- abs(dum)+1, dum <- 0)
    k_ind <- 1  ## For number of iterations
    a <- time(diff(log(dta+dum)))[diff(log(dta+dum)) < sd(diff(log(dta+dum)), na.rm=T)*-3]
    b <- time(diff(log(dta+dum)))[diff(log(dta+dum)) > 3*sd(diff(log(dta+dum)), na.rm=T)]
    while(length(c(a,b)) > 0 & k_ind <= 4){  ## 4 iterations
      k_ind <- k_ind +1
      dta[c(a,b)] <- NA
      dta <- na.approx(na.trim(dta))
      dum <- min(dta, na.rm=T)
      ifelse(dum < 0, dum <- abs(dum)+1, dum <- 0)
      a <- time(diff(log(dta+dum)))[diff(log(dta+dum)) < sd(diff(log(dta+dum)), na.rm=T)*-3]
      b <- time(diff(log(dta+dum)))[diff(log(dta+dum)) > 3*sd(diff(log(dta+dum)), na.rm=T)]
    }
  } else {  
    for (i in 1:dim(dta)[2]){
      dum <- min(dta[,i], na.rm=T)
      k_ind <- 1
      ifelse(dum < 0, dum <- abs(dum)+1, dum <- 0)
      a <- time(diff(log(dta[,i]+dum)))[diff(log(dta[,i]+dum)) < sd(diff(log(dta[,i]+dum)), na.rm=T)*-3]
      b <- time(diff(log(dta[,i]+dum)))[diff(log(dta[,i]+dum)) > 3*sd(diff(log(dta[,i]+dum)), na.rm=T)]
      while(length(c(a,b)) > 0 & k_ind < 4){  
        k_ind <- k_ind +1
        dta[c(a,b),i] <- NA
        dta[,i] <- na.approx(na.trim(dta[,i]))
        dum <- min(dta[,i], na.rm=T)
        ifelse(dum < 0, dum <- abs(dum)+1, dum <- 0)
        a <- time(diff(log(dta[,i]+dum)))[diff(log(dta[,i]+dum)) <  sd(diff(log(dta[,i]+dum)), na.rm=T)*-3]
        b <- time(diff(log(dta[,i]+dum)))[diff(log(dta[,i]+dum)) > 3*sd(diff(log(dta[,i]+dum)), na.rm=T)]
      }
      
    }
  }
  return(dta)
}

yc2.z <- remove_bad_data(yc.z)

plot(yc2.z, main="Yield Curve Data (Cleaned)", xlab="Date")


## Lets continue with Nelson-Siegel

library(YieldCurve)
maturity.tcmb <- c(1/365, 1/52, 2/52, 1/12, 3/12, 0.5, 9/12, 1,2,3,4,5,10,15,20,30)
NSParameters.d <- Nelson.Siegel(rate=yc2.z,  maturity=maturity.tcmb)
head(NSParameters.d)
tail(NSParameters.d)
plot.zoo(NSParameters.d[,1:3],plot.type ="single",main="Level, Slope and Curvature of the Yield Curve",xlab="Date", ylab="NS Parameters",col = c("red","blue","green"))
legend("top", legend=c("Level","Slope","Curvature"), col=c("red","blue","green"), lty=1, bg="transparent")
grid()



# Level

level<-NSParameters.d[,1]
head(level)
tail(level)
level_fsi<-window(level,start="2005-01-11", end="2016-11-17")
plot(level_fsi)

#Realized volatility of Level by GARCH(1,1) model

garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), distribution.model = "norm")

# Fit the model

garch.fit_level = ugarchfit(garch11.spec, data = level_fsi, fit.control=list(scale=TRUE))
print(garch.fit_level)
plot(sigma(garch.fit_level))


# Realized volatility of Level by GARCH(1,1) model with ARMA(1,1)

garch11_mean11.spec = ugarchspec(mean.model = list(armaOrder = c(1,1)), 
                                 variance.model = list(garchOrder = c(1,1), 
                                                       model = "sGARCH"), distribution.model = "norm")
# Fit the model

garch11mean11.fit_level = ugarchfit(garch11_mean11.spec, data = level_fsi, fit.control=list(scale=TRUE))
print(garch11mean11.fit_level)
plot(sigma(garch11mean11.fit_level))

levelfsi<-xts(sigma(garch11mean11.fit_level))
levelfsi_scale<-scale(levelfsi)
colnames(levelfsi_scale)<-"level"

# Event Analysis

chart.TimeSeries(levelfsi_scale, colorset = "black", legend.loc = "topleft", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Level", xlab = "Date", main="Level of the Yield Curve")

# Slope

slope<-NSParameters.d[,2]
head(slope)
tail(slope)
slope_fsi<-window(slope, start="2005-01-11", end="2016-11-17")
plot(slope_fsi)

# Acf of slope_fsi

acf(as.vector(slope_fsi))

# Pacf of slope_fsi

pacf(as.vector(slope_fsi))



# Fit the model by GARCH(1,1)

garch.fit_slope = ugarchfit(garch11.spec, data = slope_fsi, fit.control=list(scale=TRUE))
print(garch.fit_slope)
plot(sigma(garch.fit_slope))



# Fit the model by GARCH(1,1) with ARMA(1,1)

garch11mean11.fit_slope = ugarchfit(garch11_mean11.spec, data = slope_fsi, fit.control=list(scale=TRUE))
print(garch11mean11.fit_slope)
plot(sigma(garch11mean11.fit_slope))

# GARCH(1,1) with mean arma order(2,0)

garch11_mean20.spec = ugarchspec(mean.model = list(armaOrder = c(2,0)), 
                                 variance.model = list(garchOrder = c(1,1), 
                                                       model = "sGARCH"), distribution.model = "norm")

# Fit the model

garch11mean20.fit_slope = ugarchfit(garch11_mean20.spec, data = slope_fsi, fit.control=list(scale=TRUE))
print(garch11mean20.fit_slope)
plot(sigma(garch11mean20.fit_slope))

slope_fsi<-xts((sigma(garch11mean11.fit_slope)))
slope_fsi_scale<-scale(slope_fsi)
colnames(slope_fsi_scale)<-"slope"


# Event Analysis

# png("Slope of Yield Curve.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(slope_fsi_scale, colorset = "black", legend.loc = "topright", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Slope", 
                 xlab = "Date", main="Slope of the Yield Curve")
# dev.off()

# Curvature

curvature<-NSParameters.d[,3]
curvature_fsi<-window(curvature, start="2005-01-11", end="2016-11-17")
plot(curvature_fsi)


# Acf of curvature_fsi

acf(as.vector(curvature_fsi))

# Pacf of curvature_fsi

pacf(as.vector(curvature_fsi))

#Realized volatility of Curvature by GARCH(1,1) model

# Fit the model

garch.fit_curvature = ugarchfit(garch11.spec, data = curvature_fsi, fit.control=list(scale=TRUE))
print(garch.fit_curvature)
plot(sigma(garch.fit_curvature))

# GARCH(1,1) with mean arma order(1,1)


# Fit the model

garch11mean11.fit_curvature = ugarchfit(garch11_mean11.spec, data = curvature_fsi, fit.control=list(scale=TRUE))
print(garch11mean11.fit_curvature)
plot(sigma(garch11mean11.fit_curvature))

# GARCH(1,1) with mean arma order(2,0)

# Fit the model

garch11mean20.fit_curvature = ugarchfit(garch11_mean20.spec, data = curvature_fsi, fit.control=list(scale=TRUE))
print(garch11mean20.fit_curvature)
plot(sigma(garch11mean20.fit_curvature))

curvature_fsi<-xts((sigma(garch11mean11.fit_curvature)))
curvature_fsi_scale<-scale(curvature_fsi)
colnames(curvature_fsi_scale)<-"curvature"

# Event Analysis

chart.TimeSeries(curvature_fsi_scale, colorset = "black", legend.loc = "topleft", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Curvature", xlab = "Date", main="Curvature of the Yield Curve")

# CDS
turkeycds<-read.csv(file="turkeycds.csv", header=TRUE)
turkeycds.zoo<-zoo(turkeycds[,-1], order.by=as.Date(strptime(as.character(turkeycds[,1]), "%m/%d/%Y")))
head(turkeycds.zoo)
tail(turkeycds.zoo)

# png("Turkey CDS.png", width=8, height=6, units="in", res=600)
plot(turkeycds.zoo,xlab="Date",ylab="CDS", main="Turkey CDS Data")

#dev.off()

plot(diff(turkeycds.zoo))


# Realized Volatility

acf(as.vector(turkeycds.zoo))
pacf(as.vector(turkeycds.zoo))

# Volatility of CDS by GARCH(1,1) model with ARMA(1,1)

acf(as.vector(turkeycds.zoo))
pacf(as.vector(turkeycds.zoo))
garch11mean11.spec = ugarchspec(mean.model = list(armaOrder = c(1,1)), 
                                variance.model = list(garchOrder = c(1,1), 
                                                      model = "sGARCH"), distribution.model = "norm")

# Fit the model
garch.fit_cds11 = ugarchfit(garch11mean11.spec, data = turkeycds.zoo, fit.control=list(scale=TRUE))
print(garch.fit_cds11)

# Volatility of CDS by GARCH(1,1) model with ARMA(1,1)

garch.fit_cds11 = ugarchfit(garch11_mean11.spec, data = turkeycds.zoo, fit.control=list(scale=TRUE))
print(garch.fit_cds11)


# Volatility of CDS by GARCH(1,1) model with ARMA(2,0)

# Fit the model

garch.fit11mean20_cds = ugarchfit(garch11_mean20.spec, data = turkeycds.zoo, fit.control=list(scale=TRUE))
print(garch.fit11mean20_cds)

# Common period

cds_fsi<-window(xts(sigma(garch.fit11mean20_cds)),start="2005-01-11", end="2016-11-17")

cds_fsi_scale<-scale(cds_fsi)
colnames(cds_fsi_scale)<-"cds"

# Event Analysis

# png("CDS.png", width=8, height=6, units="in", res=600)

chart.TimeSeries(cds_fsi_scale, colorset = "black", legend.loc = "topright", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "CDS", 
                 xlab = "Date", main="CDS of Turkey")
# dev.off()

# 10 yr benchmark Government Bond Index

tr10yrbenchmark<-read.csv(file="tr10yrbenchmark.csv", header=TRUE)
tr10yrbenchmark.zoo<-zoo(tr10yrbenchmark[,-1], order.by=as.Date(strptime(as.character(tr10yrbenchmark[,1]), "%m/%d/%Y")))
head(tr10yrbenchmark.zoo)
tail(tr10yrbenchmark.zoo)

# Realized Volatility of 10 yr benchmark Government Bond Index  by GARCH(1,1) model with ARMA(0,0)

pacf(as.vector(tr10yrbenchmark.zoo))

# Fit the model

garch.fit_tr10yr = ugarchfit(garch11.spec, data = tr10yrbenchmark.zoo, fit.control=list(scale=TRUE))

print(garch.fit_tr10yr)


# Volatility of 10 yr benchmark Government Bond Index  by GARCH(1,1) model with ARMA(1,1)

garch.fit_tr10yr11 = ugarchfit(garch11_mean11.spec, data = tr10yrbenchmark.zoo, fit.control=list(scale=TRUE))

print(garch.fit_tr10yr11)

# Volatility of 10 yr benchmark Government Bond Index  by GARCH(1,1) model with ARMA(1,0)


garch11_mean10.spec = ugarchspec(mean.model = list(armaOrder = c(1,0)), 
                                 variance.model = list(garchOrder = c(1,1), 
                                                       model = "sGARCH"), distribution.model = "norm")


garch.fit_tr10yr10 = ugarchfit(garch11_mean10.spec, data = tr10yrbenchmark.zoo, fit.control=list(scale=TRUE))

print(garch.fit_tr10yr10)


# Volatility of 10 yr benchmark Government Bond Index  by GARCH(1,1) model with ARMA(2,0)

# Fit the model

garch.fit20_tr10yr = ugarchfit(garch11_mean20.spec, data = tr10yrbenchmark.zoo, fit.control=list(scale=TRUE))

print(garch.fit20_tr10yr)

plot(garch.fit20_tr10yr, which="all")

plot(sigma(garch.fit20_tr10yr))

# Common Period

tr10yr_fsi<-window(xts(sigma(garch.fit20_tr10yr)), start="2005-01-11", end="2016-11-17")
tr10yr_fsi_scale<-scale(tr10yr_fsi)
colnames(tr10yr_fsi_scale)<-"10 yr. Gov. Bond Index"

#Event Analysis

# png("10 yr Gov. Yield.png", width=8, height=6, units="in", res=600)

chart.TimeSeries(tr10yr_fsi_scale, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Index", xlab = "Date", main="10 yr. Gov. Bond Index")

# dev.off()

# Principal Component Analysis

bond_fsi<-merge(tr10yr_fsi_scale,cds_fsi_scale,slope_fsi_scale)
bond_fsi[,3]<-na.approx(na.trim(bond_fsi[,3], side="both"))
head(bond_fsi)
tail(bond_fsi)

# Bond Index

bondfsi<-zoo(bond_fsi[,1])

for(i in 1:length(bondfsi)){
  bondfsi[i]=1/3*bond_fsi[i,1]+1/3*bond_fsi[i,2]+
    1/3*bond_fsi[i,3]
}

bondfsi<-xts(bondfsi)
colnames(bondfsi)<-"bond index"

# Event Analysis

# png("Bond Index_11.23.16.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(bondfsi, colorset = "black", legend.loc = "topleft", period.areas = cycles.dates, 
                 period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,
                 ylab = "Index", xlab = "Date", main="Bond Market Index ")

# dev.off()

# Equity Market

xu100<-read.csv(file="bist.csv", header=TRUE)
xu100.zoo<-zoo(xu100[,-1], order.by=as.Date(strptime(as.character(xu100[,1]), "%m/%d/%Y")))
plot(xu100.zoo)
plot(diff(log(xu100.zoo)))

garch.fit_xu100 = ugarchfit(garch11_mean11.spec, data = diff(log(xu100.zoo)), fit.control=list(scale=TRUE))

print(garch.fit_xu100)

plot(garch.fit_xu100, which="all")
plot(sigma(garch.fit_xu100))

xu100idfsi<-window(xts(sigma(garch.fit_xu100)), start="2005-01-11", end="2016-11-17")
xu100idfsi_scale<-scale(xu100idfsi)


# Event Analysis


colnames(xu100idfsi_scale)<-"xu100vol"

chart.TimeSeries(xu100idfsi_scale, colorset = "black", legend.loc = "topright", 
                 period.areas = cycles.dates, period.color = "gray90", 
                 event.lines = risk.dates, event.labels = risk.labels, 
                 event.color = "red", lwd = 1.5,ylab = "Id. Vol. of XU100", xlab = "Date", 
                 main="Id. Vol. of XU100 Index")


# CMAX of Bist

cmax_xu100<-rollapply(xts(xu100.zoo), align="right", width = 504, FUN = function(x) 1-(xts::last(x)/max(x)))

str(cmax_xu100)
plot(cmax_xu100)

# Common Period

cmax_xu100_fsi<-window(cmax_xu100, start="2005-01-11", end="2016-11-17")
cmax_xu100_fsi_scale<-scale(cmax_xu100_fsi)
colnames(cmax_xu100_fsi_scale)<-"CMAX of XU100"

# Event Analysis

# png("CMAX of XU100.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(cmax_xu100_fsi_scale, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "CMAX of XU100", xlab = "Date", main="CMAX of Stock Index")
# dev.off()

# Difference between CMAX of Turkey's main stock market index and CMAX of US's main stock market index (SP500).

sp500<-read.csv(file="sp500.csv", header = T)
sp500.zoo<-zoo(sp500[,-1], order.by=as.Date(strptime(as.character(sp500[,1]), "%m/%d/%Y")))
sp500.zoo<-window(sp500.zoo,start=start(xu100.zoo), end=end(xu100.zoo))

# CMAX of sp500

cmax_sp500<-rollapply(xts(sp500.zoo), align="right", width = 504, FUN = function(x) 1-(xts::last(x)/max(x)))


cmax_sp500_fsi<-window(cmax_sp500, start="2005-01-11", end="2016-11-17")

# Difference between xu100 and sp500

cmax_xu100sp500_fsi<-cmax_xu100_fsi-cmax_sp500_fsi


cmax_xu100sp500_fsi_scale<-scale(cmax_xu100sp500_fsi)
colnames(cmax_xu100sp500_fsi_scale)<-" CMAX Diff. "

# Event Analysis

# png("CMAX Difference.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(cmax_xu100sp500_fsi_scale, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "CMAX Diff.", xlab = "Date", main="Difference of CMAX of Stock Indexes")
# dev.off()

# Principal Component Analysis

equity_fsi<-merge(xu100idfsi_scale, cmax_xu100_fsi_scale,cmax_xu100sp500_fsi_scale)
head(equity_fsi)
tail(equity_fsi)



equityfsi<-zoo(equity_fsi[,1])

for(i in 1:length(equityfsi)){
  equityfsi[i]=1/3*equity_fsi[i,1]+1/3*equity_fsi[i,2]+
    1/3*equity_fsi[i,3]
}

# Event Analysis

equityfsi<-xts(equityfsi)
colnames(equityfsi)<-"equity index"

# png("Equity Index_11.29.2016.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(equityfsi, colorset = "black", legend.loc = "topleft", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Index.", xlab = "Date", main="Equity Market Index")
# dev.off()


# Money Market


# 3 month tr libor
tr3libor<-read.csv(file="tr3libor.csv", header=TRUE)
tr3libor.zoo<-zoo(tr3libor[,-1], order.by=as.Date(strptime(as.character(tr3libor[,1]), "%m/%d/%Y")))
head(tr3libor.zoo)
tail(tr3libor.zoo)
plot(tr3libor.zoo)

# Realized volatility

pacf(as.vector(tr3libor.zoo))

# GARCH(1,1) with arma(1,1)

garch.fit_tr3libor = ugarchfit(garch11_mean11.spec, data = tr3libor.zoo, fit.control=list(scale=TRUE))
print(garch.fit_tr3libor)

# Common Period

tr3libor_fsi<-window(xts(sigma(garch.fit_tr3libor)),start="2005-01-11", end="2016-11-17")
tr3libor_fsi_scale<-scale(tr3libor_fsi)

# Event Analysis

colnames(tr3libor_fsi_scale)<-"tr3libor"

# png("TR Libor.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(tr3libor_fsi_scale, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "3 month libor rate", xlab = "Date", main="3 Month Interbank Rate ")
# dev.off()

# Ted Spread (Spread between 3 month libor and 3 month Treasury Yield) 


tr3libor_1.zoo<-window(tr3libor.zoo,start="2005-01-11", end="2016-11-17")
ted<-merge(tr3libor_1.zoo,tr3myc.zoo[2:length(tr3myc.zoo)])
ted[,2]<-na.approx(na.trim(ted[,2],side="both"))
ted<-ted[,1]-ted[,2]
plot(ted)

# Realized volatility with GARCH(1,1) 

garch.fit_ted = ugarchfit(garch11_mean20.spec, data = ted, fit.control=list(scale=TRUE))
print(garch.fit_ted)
plot(sigma(garch.fit_ted))

ted_fsi<-window(xts(sigma(garch.fit_ted)),start="2005-01-11", end="2016-11-17")
ted_fsi_scale<-scale(ted_fsi)

# Event Analysis

colnames(ted_fsi_scale)<-"ted spread"

# png("TED.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(ted_fsi_scale, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "TED Spread", xlab = "Date", main="Spread between 3 month Interbank and 3 month Treasury Yield ")
# dev.off()




money_fsi<-merge(tr3libor_fsi_scale,ted_fsi_scale)


# Money Index

moneyfsi<-zoo(money_fsi[,1])

for (i in 1:length(moneyfsi)){
  moneyfsi[i]=1/2*money_fsi[i,1]+1/2*money_fsi[i,2]
}

# Event Analysis


colnames(moneyfsi)<-"money index"

# png("Money Index_11.29.2016.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(moneyfsi, colorset = "black", legend.loc = "topleft", 
                 period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, 
                 event.labels = risk.labels, event.color = "red", lwd = 1.5,
                 ylab = "Index", xlab = "Date", main="Money Market Index")
# dev.off()


# Foreign Exchange Market

# USD/TRY

tlusd<-read.csv(file="tlusd.csv", header=T, sep=";")
tlusd.zoo<-zoo(tlusd[,-1], order.by=as.Date(strptime(as.character(tlusd[,1]), "%d-%m-%Y")))
tlusd.zoo<-tlusd.zoo[!is.weekend(time(tlusd.zoo))]
tlusd.zoo[tlusd.zoo==0]<-NA
#tlusd.zoo<- na.approx(na.trim(CalculateReturns(tlusd.zoo), side="both"))
tlusd.zoo<-na.approx(na.trim(tlusd.zoo,side="both"))
head(tlusd.zoo)
tail(tlusd.zoo)
plot(tlusd.zoo)


# Realized Volatility of USD/TRY with GARCH(1,1) and ARMA(0,0)

pacf(as.vector(tlusd.zoo))

garch.fit_tlusd = ugarchfit(garch11_mean20.spec, data = tlusd.zoo, fit.control=list(scale=TRUE))
print(garch.fit_tlusd)

plot(sigma(garch.fit_tlusd))

tlusd_fsi<-window(xts(sigma(garch.fit_tlusd)),start="2005-01-11", end="2016-11-17")
tlusd_fsi_scale<-scale(tlusd_fsi)


# EURO/TRY

tleuro<-read.csv(file="tleuro.csv", header=T, sep=";")
tleuro.zoo<-zoo(tleuro[,-1], order.by=as.Date(strptime(as.character(tleuro[,1]), "%d-%m-%Y")))
tleuro.zoo <- tleuro.zoo[!is.weekend(time(tleuro.zoo))]
tleuro.zoo[tleuro.zoo==0]<-NA
#tleuro.zoo<- na.approx(na.trim(CalculateReturns(tleuro.zoo), side="both"))
tleuro.zoo<- na.approx(na.trim(tleuro.zoo,side="both"))
head(tleuro.zoo)
tail(tleuro.zoo)
plot(tleuro.zoo)


# Realized Volatility of EURO/TRY with GARCH(1,1) and ARMA(0,0)

pacf(as.vector(tleuro.zoo))

garch.fit_tleuro = ugarchfit(garch11.spec, data = tleuro.zoo, fit.control=list(scale=TRUE))
print(garch.fit_tleuro)

plot(sigma(garch.fit_tleuro))

tleuro_fsi<-window(xts(sigma(garch.fit_tleuro)),start="2005-01-11", end="2016-11-17")
tleuro_fsi_scale<-scale(tleuro_fsi)



# Index

fx_fsi<-merge(tleuro_fsi_scale,tlusd_fsi_scale)


fxfsi<-zoo(fx_fsi[,1])

for(i in 1:length(fxfsi)){
  fxfsi[i]=1/2*fx_fsi[i,1]+1/2*fx_fsi[i,2]
}

# Event Analysis

fxfsi<-xts(fxfsi)
colnames(fxfsi)<-"fx index"

# png("FX Index_11.29.2016.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(fxfsi, colorset = "black", legend.loc = "topleft", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, event.color = "red", lwd = 1.5,ylab = "Index", xlab = "Date", main="Foreign Exchange Market Index")
# dev.off()

## Constructing Index


fsi<-merge(bankfsi,bondfsi,moneyfsi,equityfsi,fxfsi)


# CISS analysis

# DCC-GARCH

fsi_dcc_1<-xts(fsi)
garch.spec<-ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "norm")

dcc.spec<-dccspec(uspec = multispec( replicate(5, garch.spec) ), dccOrder = c(1,1), distribution = "mvnorm")

dcc.fit<-dccfit(dcc.spec, data =fsi_dcc_1 , fit.control=list(scale=TRUE))

cormat<-rcor(dcc.fit, type="R")
cormatrix<-cormat[1:5,1:5,]
weights<-c(0.2,0.2,0.2,0.2,0.2) # Obtained by Principal Component Analysis


fsi_dcc<-zoo(fsi_dcc_1)

j=1
for (i in 1:length(fsi_dcc[,1])){
  for (j in 1:length(weights))
    fsi_dcc[,j][i]=weights[j]*fsi_dcc_1[,j][i]
  j=j+1
}

S<-matrix(fsi_dcc,nrow=length(fsi_dcc[,1]),ncol=5)
S_tr<-t(S)

ciss<-zoo(fsi_dcc[,1])

for (i in length(fsi_dcc[,1])){
  ciss[i]=(S[i,]%*%cormat[1:5,1:5,i]%*%S_tr[,i])
}

##Function
##x is merge of sub-market indexes
CISS<-function(x){
  garch.spec<-ugarchspec(mean.model = list(armaOrder = c(1,1)), 
                         variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
                         distribution.model = "norm")
  dcc.spec<-dccspec(uspec = multispec(replicate(5, garch.spec)), dccOrder = c(1,1), distribution = "mvnorm")
  dcc.fit<-dccfit(dcc.spec, data =x , fit.control=list(scale=TRUE))
  cormat<-rcor(dcc.fit, type="R")
  cormatrix<-cormat[1:5,1:5,]
  weights<-c(0.2,0.2,0.2,0.2,0.2) # Equal weigths
  y<-xts(x)
  j=1
  for (i in 1:length(y[,1])){
    for (j in 1:length(weights))
      y[,j][i]=weights[j]*x[,j][i]
    j=j+1
  }
  S<-matrix(y,nrow=length(y[,1]),ncol=5)
  S_tr<-t(S)
  ciss<-xts(y[,1])
  for (i in length(y[,1])){
    ciss[i]=(S[i,]%*%cormat[1:5,1:5,i]%*%S_tr[,i])
  }
  return(ciss)
}

ciss_1<-CISS(fsi)


# Event Analysis

ciss<-xts(ciss)
colnames(ciss)<-"ciss" 

# png("CISS_11.29.2016.png", width=8, height=6, units="in", res=600)
chart.TimeSeries(ciss, colorset = "black", legend.loc = "topright", period.areas = cycles.dates, period.color = "gray90", event.lines = risk.dates, event.labels = risk.labels, 
                 event.color = "red", lwd = 1.5,ylab = "CISS", xlab = "Date", main="Financial Stress Index for Turkey")
dev.off()


ciss<-xts(ciss)
colnames(ciss)<-"ciss" 

# png("CISS_plot11.23.2016.png", width=8, height=6, units="in", res=600)
plot(ciss, xlab="Date", ylab="CISS", main = "Financial Stress Index for Turkey")
# dev.off()


