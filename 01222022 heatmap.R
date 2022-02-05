install.packages("dplyr")
install.packages("Rtools")
install.packages("ggtree")
install.packages("ggstance")
install.packages("cowplot")
install.packages("patchwork")

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
library(tidyverse)
library(dplyr)
library(ggtree)
library(ggstance)
library(cowplot)
library(patchwork)

data0 = read.csv("C:/Users/sjeon/Desktop/dementia R01 11012021/mortality_mi_10052021_count_20_long_2.csv")
data0 = subset(data0, select=c("time_to_death", "death", "age_cat", "ragender", "bmi_1", "bmi_3", "bmi_4", "adl_dep_count", "iadl_diff_count",
                 "rwalks", "rcancre", "rdiabe", "rhearte", "rlunge", "smoke_cat", "vigact", "X_mi_m"))
data_original = data0[1:4267,]

pred_cox1=NULL

for(i in 1:20){
  
  data = data0[which(data0$X_mi_m == i),]
  cox1 = coxph(Surv(time_to_death, death) ~ as.factor(age_cat) + as.factor(ragender) + bmi_1 + bmi_3 + bmi_4 + adl_dep_count + iadl_diff_count +
                 rwalks + rcancre + rdiabe + rhearte + rlunge + as.factor(smoke_cat) + vigact, data=data, x=TRUE )
  #, weights=data$baseline_svywgt*/)
  
  x = survfit(cox1, newdata=data)
  pred_cox1= cbind(pred_cox1, predict(cox1, data, type=c("risk")))
}


predict_cox = as.data.frame(cbind(pred_cox1, pred_risk=rowMeans(pred_cox1)))

predict_cox = predict_cox %>% 
              mutate(quantile=ntile(predict_cox[,21], 10))

predict_cox_cov = cbind(predict_cox,data_original )

# set.seed(12345)
# selected_patients = predict_cox_cov %>% group_by(quantile) %>%sample_n(1)
# (selected_patients[1:10, 25:39])
# (selected_patients[1:10, 20:39])
# 
# set.seed(12346)
# selected_patients = predict_cox_cov %>% group_by(quantile) %>%sample_n(1)
# (selected_patients[1:10, 25:38])
# 
# set.seed(123)
# selected_patients = predict_cox_cov %>% group_by(quantile) %>%sample_n(1)
# (selected_patients[1:10, 25:38])

set.seed(28) ###THIS IS IT
set.seed(79) ###THIS IS IT
selected_patients0 = predict_cox_cov %>% group_by(quantile) %>%sample_n(1)
(selected_patients0[1:10, 25:38])
selected_patients0[is.na(selected_patients0)] =0

# Heatmap Plotting

###Practice
# Dummy data
# x <- LETTERS[1:20]
# y <- paste0("var", seq(1,20))
# data <- expand.grid(X=x, Y=y)
# data$Z <- runif(400, 0, 5)
# 
# ggplot(data, aes(X, Y, fill= Z)) + 
#   geom_tile()

selected_patients=selected_patients0[1:10, 25:38]

selected_patients$age_cat = as.factor(selected_patients$age_cat)
selected_patients$smoke_cat = as.factor(selected_patients$smoke_cat)
selected_patients$ragender = selected_patients$ragender-1

patient_mx = model.matrix(~., data=selected_patients)[,-1]
  patient_mx
  patient_mx_ref = cbind(c(1,0,0,0,0,0,0,0,1,0), patient_mx[, 1:7], 
                         c(1,0,0,1,0,0,1,0,0,1), patient_mx[,8:16], 
                         c(0,0,0,1,0,0,1,0,1,0), patient_mx[,17:19] )
  patient_mx_ref[patient_mx_ref ==0 ] = "o"
  patient_mx_ref[patient_mx_ref ==1 ] = "x"

  #12142021 collapse
  coef.mx = diag(cox1$coef) ##fix depending on values in patient_mx###
  linear_pred = patient_mx %*% coef.mx
  linear_pred
  
  age_c = rowSums(linear_pred[, 1:5])
  bmi_c = rowSums(linear_pred[, 7:9])
  smoke_c = rowSums(linear_pred[, 17:18])
#  collapsed_mx = cbind(age_c, linear_pred[,6], bmi_c, linear_pred[, 10:16], smoke_c, linear_pred[,19])
  collapsed_mx = cbind(linear_pred[,16], linear_pred[,13], linear_pred[,15], linear_pred[,14],linear_pred[,19], linear_pred[,12], 
                       linear_pred[,11], linear_pred[,10],  smoke_c, bmi_c, linear_pred[,6], age_c)
  hr_collapsed = exp(collapsed_mx)

  hr_vec = c(hr_collapsed)
  
  Patient = LETTERS[1:10]
  # Predictor = c("Age group",
  #               "Female", "BMI", 
  #               "#of ADL dependencies", "# of IADL difficulties", 
  #               "Walk, difficulty", "Cancer", "Diabetes", "Heart disease", "Lung disease", 
  #               "Smoke", 
  #               "Vigorous activity")
  
  Predictor = c( "Lung disease", "Cancer", "Heart disease", "Diabetes","Vigorous activity",
                 "Walk, difficulty",  "# of IADL difficulties", "#of ADL dependencies",  "Smoke", 
                 "BMI","Female",  "Age group"
                 )

                
  data = expand.grid(X=Patient, Y=Predictor)
  data$HazardRatio = hr_vec
  
  patient_mx = model.matrix(~., data=selected_patients)[,-1]
  patient_mx
  patient_mx_ref = cbind(patient_mx[,16], patient_mx[,13], patient_mx[,15], patient_mx[,14], patient_mx[,19], patient_mx[,12], 
                         patient_mx[,11], patient_mx[,10], 
                         c("Former","Former","Former","Never", "Former","Former","Never","Former", "Never","Current"),
                         c("18.5-25",">=30","25-30", "18.5-25", "25-30", ">=30", "18.5-25", "25-30", "<18.5", "18.5-25"),  
                         patient_mx[, 6], 
                         c("65-69","70-74","75-79","85-89", "75-79", "85-89", "80-84", "80-84", "90+", "70-74")
                         )
  patient_mx_ref[patient_mx_ref =="0" ] = " "
  patient_mx_ref[patient_mx_ref =="1" ] = "Y"
  patient_mx_ref[, 7] = c("0", "0", "1", "0", "2", "1", "4", "1", "1", "5")
  design_vect = c(patient_mx_ref)
#### BEFORE COLLASPING CATEGORICAL VARIABLES IN ONE ROW
  
  #linear_pred_ref = exp(cbind(rep(0,10), linear_pred[, 1:7], rep(0,10), linear_pred[,8:16], rep(0,10), linear_pred[,17:19] ))
  #linear_pred_vec = c(linear_pred_ref)

# Patient = LETTERS[1:10]
# Predictor = c("(ref) Age 65-69", "Age 70-74", "Age 75-79", "Age 80-84", "Age 85-89", "Age 90+",
#               "Female", "BMI<18.5", "(ref) BMI 18.5-25", "BMI 25-30", "BMI>=30", 
#               "ADL Dependence", "IADL dependence", 
#               "Walk, difficulty", "Cancer", "Diabetes", "Heart disease", "Lung disease", "(ref) Never Smoker", "Former Smoker", "Current Smoker", 
#               "Vigorous activity")
# data = expand.grid(X=Patient, Y=Predictor)
# data$HazardRatio = linear_pred_vec
# design_vect= c(patient_mx_ref)
  
data = expand.grid(X=Patient, Y=Predictor)
data$HazardRatio = hr_vec

data_mx = as.matrix(data)
data_mx[, 3] = as.numeric(data_mx[,3])
head(data_mx)

plot = ggplot(data, aes(X, Y, fill= HazardRatio, label=sprintf("%.2f", round(HazardRatio, digits = 1)))) + 
  geom_tile() + 
  scale_fill_gradient2 (low="darkgreen", mid="white", high="#c51b7d", midpoint=1, limits=c(0.694, 4.8)) + 
  ylab("Predictor")+
  xlab("Patient")+ 
  geom_text(size=3.5, aes(label=design_vect))+
  guides(fill = guide_colorbar(order=3))
  
plot

######################################
##### TIME TO EVENT predicted
######################################

pred_time=NULL
for(i in 1:20){
  
  data = data0[which(data0$X_mi_m == i),]
  cox1 = coxph(Surv(time_to_death, death) ~ as.factor(age_cat) + as.factor(ragender) + bmi_1 + bmi_3 + bmi_4 + adl_dep_count + iadl_diff_count +
                 rwalks + rcancre + rdiabe + rhearte + rlunge + as.factor(smoke_cat) + vigact, data=data, x=TRUE )
  #, weights=data$baseline_svywgt*/)
  
  x = survfit(cox1, newdata=selected_patients0)
  pred_time= cbind(pred_time, summary(x)$table[,"median"]/365.25) ## replace this quantile(x) (Quantile.survfit{survival)})
}

qt_25 = quantile(x)$quantile[,1]/365.25
qt_75 = quantile(x)$quantile[,3]/365.25

predict_time = as.data.frame(cbind(pred_time, pred_time=rowMeans(pred_time), qt_25, qt_75))[1:10, ]
predict_time$qt_75[1] = predict_time$pred_time[1]
  
rownames(predict_time) =  LETTERS[1:10]
plot2 = ggplot(predict_time) +
       geom_point(aes(x= LETTERS[1:10], y=predict_time$pred_time), size=3)+
       xlab("Patient") + 
       ylab("Median predicted time to death (years)")+
  #    ylim(0,14)+
      scale_y_continuous(breaks=c(2,4,6,8,10, 12, 14), limits=c(0,15))+ 
      geom_errorbar(aes(ymin=predict_time$qt_25, ymax=predict_time$qt_75, x=LETTERS[1:10]), width=0.15, size=0.5, position=position_dodge(.9))
  
  plot2
   
   plot / plot2 + plot_layout(heights=c(5,2))
   ggsave("C:/Users/sjeon/Desktop/work/dementiaR01/DementiaR01_01212022/01252021heatmap_1.png", height=10)
   