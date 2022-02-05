
library(survival)
library(Formula)
library(Hmisc)
library(SparseM)
library(rms)
library(prodlim)
library(pec)
library(Rcpp)
library(riskRegression)
library(crayon )
library(backports)
library(data.table)
library(plotrix)
library(coxphw)
library(cmprsk)
library(timereg)

data0 = read.csv("C:/Users/sjeon/Desktop/work/dementiaR01/DementiaR01_01212022/01222022 incident_mi.csv")

data = data0[data0$X_mi_m==0,]
cox1 = coxph(Surv(time_to_death, death) ~ as.factor(age_cat) + as.factor(ragender) + bmi_1 + bmi_3 + bmi_4 + adl_dep_count + iadl_diff_count +
               rwalks + rcancre + rdiabe + rhearte + rlunge + as.factor(smoke_cat) + vigact, data=data, x=TRUE )
#, weights=data$baseline_svywgt*/)

predict.cox = 1 - predictSurvProb(cox1, newdata=data, times=365*(1:10)) ##predictSurvProb doesn't work with weights
head(predict.cox)

data$cox.1yr = predict.cox[,1]
data$cox.2yr = predict.cox[,2]
data$cox.3yr = predict.cox[,3]
data$cox.4yr = predict.cox[,4]
data$cox.5yr = predict.cox[,5]
data$cox.10yr = predict.cox[,10]

data$cox.1yr.cll = log(-log(1-data$cox.1yr))
calibrate.cox = coxph(Surv(time_to_death, death) ~ rcs(cox.1yr.cll,3), x=T, data=data)

predict.grid.cox = seq( quantile(data$cox.1yr, probs=0.01, na.rm=TRUE), quantile(data$cox.1yr, probs=0.99, na.rm=TRUE), length=10)
predict.grid.cox.cll = log(-log(1-predict.grid.cox)) 

predict.grid.cox.df <- data.frame(predict.grid.cox)
predict.grid.cox.cll.df <- data.frame(predict.grid.cox.cll)

names(predict.grid.cox.df) <- "cox.1yr"
names(predict.grid.cox.cll.df) <- "cox.1yr.cll"

predict.calibrate.cox <- 1 - predictSurvProb(calibrate.cox,
                                             newdata=predict.grid.cox.cll.df,times=1*365)


plot(predict.grid.cox,predict.calibrate.cox,type="l",lty=1,col="red",
     xlim=c(0,1),ylim=c(0,1), lwd=2,
     xlab = "Predicted probability of 1-year mortality",
     ylab = "Observed probability of 1-year mortality", main="1-year (Austin & Harrell)")
abline(0,1)

par(new=T)
plot(density(data$cox.1yr, na.rm=TRUE),axes=F,xlab=NA,ylab=NA,main="")
axis(side=4)



#######################################################################################
######## OR Can we replicate Grisell's code?
#######################################################################################

xs1 = xs2 = xs5 = xs10 = list(0)

for (i in 1: 20){
  data = data0[which(data0$X_mi_m == i),]
  data = subset(data, select=c(age_cat, ragender, bmi_1,  bmi_3 ,  bmi_4 ,  adl_dep_count , iadl_diff_count,
                               rwalks , rcancre , rdiabe , rhearte , rlunge , smoke_cat, vigact, time_to_death, death, X_mi_m))
  
  data=na.omit(data)
  fitdth <- coxph(Surv(time_to_death, death) ~ as.factor(age_cat) + as.factor(ragender) + bmi_1 + bmi_3 + bmi_4 + adl_dep_count + iadl_diff_count +
                    rwalks + rcancre + rdiabe + rhearte + rlunge + as.factor(smoke_cat) + vigact, data=data, x=TRUE,y=TRUE)
  fitdth #Checked the estimates and se, and they are similar to the ones obtained in SAS PHREG
  
  ####### Calibration plot
  memory.limit(10000000000000)
  set.seed(1234)
  xs1[[i]]=Score(list(" "=fitdth),
                 formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*1, plot="Calibration", na.rm=TRUE)
  xs2[[i]]=Score(list(" "=fitdth),
                 formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*2, plot="Calibration", na.rm=TRUE)
  # xs4[[i]]=Score(list(" "=fitdth),
  #           formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*4, plot="Calibration", na.rm=TRUE)
  xs5[[i]]=Score(list(" "=fitdth),
                 formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*5, plot="Calibration", na.rm=TRUE)
  # xs6[[i]]=Score(list(" "=fitdth),
  #           formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*6, plot="Calibration", na.rm=TRUE)
  # xs8[[i]]=Score(list(" "=fitdth),
  #           formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*8, plot="Calibration", na.rm=TRUE)
  xs10[[i]]=Score(list(" "=fitdth),
                  formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*10, plot="Calibration", na.rm=TRUE)
  # xs12[[i]]=Score(list(" "=fitdth),
  #           formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*12, plot="Calibration", na.rm=TRUE)
  # xs14[[i]]=Score(list(" "=fitdth),
  #           formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=365.25*14, plot="Calibration", na.rm=TRUE)
}
# xsmax=Score(list(" "=fitdth),
#            formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=7779, plot="Calibration", na.rm=TRUE)
# xsseq=Score(list(" "=fitdth),
#             formula=Surv(time_to_death,death)~1,data=data,conf.int=TRUE,times=c(365.25*2,365.25*4,365.25*6,365.25*8,365.25*10,
#                                                                                 65.25*12,365.25*14), plot="Calibration", na.rm=TRUE)
# 
# calfitdth10<-plotCalibration(xs14, cens.method="local") 

#plotCalibration(xs1, cens.method="local", title= "Times=1year")
# title(main="Times=1year")
par(mfrow=c(2,2))
par(cex.axis=1.8, cex.lab=1.8)
par(mar=c(5,7,3,3))
par(mgp=c(5,1,0))
plotCalibration(xs1[[1]], cens.method="local", auc.in.legend =FALSE,
                brier.in.legend = FALSE, col="Blue", cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
                xlab="") 
title(xlab="Predicted Risk", line=3.5)

par(new=T)
plot(density(data$cox.1yr, na.rm=TRUE),axes=F,xlab=NA,ylab=NA,main="", lwd=1.5)
title(main="Year 1", cex.main=2)
mtext(cex=1.7, "AUC 0.710 (0.685, 0.736)", line=-0.9)


##YEAR 2
plotCalibration(xs2[[1]], cens.method="local", auc.in.legend =F,
                brier.in.legend = FALSE, col="Blue", cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
                xlab="") 

title(xlab="Predicted Risk", line=3.5)
par(new=T)
plot(density(data$cox.2yr, na.rm=TRUE),axes=F,xlab=NA,ylab=NA,main="", lwd=1.5)
title(main="Year 2", cex.main=2)
mtext(cex=1.7, "AUC 0.718 (0.699, 0.737)", line=-0.9)


## YEAR 5
plotCalibration(xs5[[1]], cens.method="local", auc.in.legend =F,
                brier.in.legend = FALSE, col="Blue", cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
                xlab="" )

title(xlab="Predicted Risk", line=3.5)
par(new=T)
plot(density(data$cox.5yr, na.rm=TRUE),axes=F,xlab=NA,ylab=NA,main="", lwd=1.5)
title(main="Year 5", cex.main=2)
mtext(cex=1.7, "AUC 0.746 (0.729, 0.762)", line=-0.9)


## YEAR10
plotCalibration(xs10[[1]], cens.method="local", auc.in.legend =F,
                brier.in.legend = FALSE, col="Blue", cex.lab=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, 
                xlab="" )

title(xlab="Predicted Risk", line=3.5)
par(new=T)
plot(density(data$cox.10yr, na.rm=TRUE),axes=F,xlab=NA,ylab=NA,main="", lwd=1.5)
title(main="Year 10", cex.main=2)
mtext(cex=1.7, "AUC 0.830 (0.809, 0.851)", line=-0.9)
