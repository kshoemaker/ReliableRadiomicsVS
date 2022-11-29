### Simulation Design ###
## Week of Jan 18, 2019
  
# setup 
rm(list = ls())
library(MASS)
library(BAMMtools)
library(tidyverse)
library(gtools)


#set.seed(12345678) ## seed good for stages?s
set.seed(1234)
#setwd("~/Google Drive/Classes/CP_Project/LungProjectData")
#setwd("C:/Users/KAShoemaker/Box Sync/Stats + radiomics project")
setwd("~/Documents/GitHub/StatinMedicineProject/CaseStudyData/")
# 78 subjects
hn_data <- read.csv("CTPatientDataForKate.csv", stringsAsFactors = T)
# hn_data <- hn_data[c(-81,-82),]

## load in reliability information 
reliability_data <- read.csv("LungProjectData/TotalSDForKate.csv", header = T, stringsAsFactors = F)
#reliability_data$Control.Total.SD <- log(reliability_data$Control.Total.SD)






hn_aerts <- filter(hn_data,DataType =="A") ## ,  HPV.status == "Negative")
# hn_aerts <- hn_aerts[c(-81,-82),]
# hn_aerts <- hn_data
hn_clinical <- hn_aerts[,1:14]
 hn_radiomics <- hn_aerts[,15:186]
# hn_radiomics <- hn_aerts[,58:100]
hn_radiomics <- hn_radiomics[,which(colnames(hn_radiomics) %in% reliability_data$PreFeat)]


# par(mfrow = c(1,2))
# for(i in 1:64){
#   hist(hn_radiomics[,i], main = paste("Histogram of Variable ",i))
#   if(min(hn_radiomics[,i] < 0)){
#     hist(log(hn_radiomics[,i] - min(hn_radiomics[,i]) + 0.1),main = paste("Histogram of Variable ",i))}
#   else {hist(log(hn_radiomics[,i]),main = paste("Histogram of Variable ",i)) }
# }

#hn_radiomics[,160] <- hn_radiomics[,160]*100000
library(car)
trans_hn_radiomics <- hn_radiomics
for (i in 1:(ncol(hn_radiomics))){
  vec <- hn_radiomics[,i]
  if (min(vec) <= 0){
    fam = "bcnPower"
  } else {
      gam <- 0
      fam <- "bcPower"
    }
  if (mean(vec)<0){ print(paste("vector", i,"has negative mean"))}
  else {
  pT <- powerTransform(vec, family = fam)
    if (pT$roundlam != 1){
    lambda <- pT$roundlam
    gamma <- max(pT$gamma,0)
    trans_hn_radiomics[,i] <- bcnPower(vec, lambda = lambda, gamma = gamma)
    print(paste("transformed vector" , i,"with lambda = ", lambda))
    } else { print(paste("vector",i,"doesn't need transform"))}
  }
  hist(trans_hn_radiomics[,i])
}
# hn_radiomics[,149] <- log(hn_radiomics[,149])
# hn_radiomics[,136] <- log(hn_radiomics[,136])
# hn_radiomics[,132] <- log(hn_radiomics[,132])
# hn_radiomics[,131] <- log(hn_radiomics[,131])
# hn_radiomics[,129] <- log(hn_radiomics[,129])
# hn_radiomics[,128] <- log(hn_radiomics[,128])
# hn_radiomics[,127] <- log(hn_radiomics[,127])
# hn_radiomics[,124] <- log(hn_radiomics[,124])
# hn_radiomics[,123] <- log(hn_radiomics[,123])
# hn_radiomics[,121] <- log(hn_radiomics[,121])
# hn_radiomics[,114] <- log(hn_radiomics[,114])
# hn_radiomics[,113] <- log(hn_radiomics[,113])
# hn_radiomics[,112] <- log(hn_radiomics[,112])
# hn_radiomics[,109] <- log(hn_radiomics[,109])
# hn_radiomics[,92] <- log(hn_radiomics[,92])
# hn_radiomics[,91] <- log(hn_radiomics[,91])
# hn_radiomics[,89] <- log(hn_radiomics[,89])
# hn_radiomics[,88] <- log(hn_radiomics[,88])
# hn_radiomics[,87] <- log(hn_radiomics[,87])
# hn_radiomics[,86] <- log(hn_radiomics[,86] - min(hn_radiomics[,86]) + 0.1)
# hn_radiomics[,84] <- log(hn_radiomics[,84])
# hn_radiomics[,83] <- log(hn_radiomics[,83])
# hn_radiomics[,81] <- log(hn_radiomics[,81])
# hn_radiomics[,69] <- log(hn_radiomics[,69])
# hn_radiomics[,56] <- log(hn_radiomics[,56])
# hn_radiomics[,55] <- log(hn_radiomics[,55])
# hn_radiomics[,52] <- log(hn_radiomics[,52])
# hn_radiomics[,51] <- log(hn_radiomics[,51])
# hn_radiomics[,49] <- log(hn_radiomics[,49])
# hn_radiomics[,48] <- log(hn_radiomics[,48]+0.1)
# hn_radiomics[,47] <- log(hn_radiomics[,47])
# hn_radiomics[,41] <- log(hn_radiomics[,41])
# hn_radiomics[,43] <- log(hn_radiomics[,43]+0.1)
# hn_radiomics[,44] <- log(hn_radiomics[,44]+0.01)
# hn_radiomics[,29] <- log(hn_radiomics[,29])
# hn_radiomics[,15] <- log(hn_radiomics[,15])
# hn_radiomics[,12] <- log(hn_radiomics[,12])
# hn_radiomics[,11] <- log(hn_radiomics[,11])
# hn_radiomics[,9] <- log(hn_radiomics[,9])
# hn_radiomics[,8] <- log(hn_radiomics[,8] + 0.1)
# hn_radiomics[,7] <- log(hn_radiomics[,7])
# hn_radiomics[,6] <- log(hn_radiomics[,6]-min(hn_radiomics[,6])+0.1)
# hn_radiomics[,4] <- log(hn_radiomics[,4]+0.1)
# hn_radiomics[,3] <- log(hn_radiomics[,3]+0.1)
# hn_radiomics[,1] <- log(hn_radiomics[,1])
# hn_radiomics[,142] <- exp(hn_radiomics[,142])
# hn_radiomics[,140] <- exp(hn_radiomics[,140])
# hn_radiomics[,102] <- exp(hn_radiomics[,102])
var_aerts <- apply(hn_radiomics,2,var)

hn_radiomics <- trans_hn_radiomics
hn_radiomics$ShapeVolume <- hn_aerts$ShapeVolume  
hn_radiomics$Age <- hn_aerts$Age
hn_radiomics <- scale(hn_radiomics,T,T)

# N <- N[which(reliability_data$PreFeat %in% colnames(hn_radiomics))]
N <- reliability_data$Control.Total.SD[match(colnames(hn_radiomics), reliability_data$PreFeat)]
# N[41:42] <- mean(N, na.rm = T)
 N[161:162] <- mean(N, na.rm = T)
N <- log(N)
N <- abs(N- max(N))
N <- N/max(N)


## NO GENES
n <- nrow(hn_radiomics)
gene_data <- matrix(rep(0,n*4),ncol = 4)


train_size <- floor(n*3/4)
test_size <- n - train_size
train <- sample(1:n,train_size,replace = F) #seq(1,101,by = 2)  #

            
            
#make groups
groups <- kmeans(hn_clinical$OStime,2)
type <- groups$cluster[train] - 1 
test_type <- groups$cluster[-train] - 1
# trainingOStime <- hn_clinical$OStime[train]
# a <-  max(trainingOStime[groups$cluster == 1])
# b <- min(trainingOStime[groups$cluster == 2])
# border <- mean(a,b)
# type <- groups$cluster - 1 
# test_type <- as.numeric(hn_clinical$OStime[-train] > border)

# what if we do stage?
type2 <- as.integer(hn_clinical$AJCC.stage[train]) -1
test_type2 <- as.integer(hn_clinical$AJCC.stage[-train]) - 1

# what if we do HPV status? 
type3 <- as.integer(hn_clinical$HPV.status[train]) -1
test_type3 <- as.integer(hn_clinical$HPV.status[-train]) - 1


X <- (hn_radiomics[train,])
Xf <- (hn_radiomics[-train,])
Y <- type3
Yf <- test_type3
Z <- as.matrix(gene_data[train,])
Zf <- as.matrix(gene_data[-train,])
N <- t(N)


# library(R.matlab)
# #setwd("C:/Users/KAShoemaker/Google Drive/Classes/CP_Project/ProbitModelCode/")
# setwd("~/Documents/GitHub/StatinMedicineProject/MatLabCode/")
# writeMat("HNModelData_Dec20HPV.mat", X = X, Xf = Xf, Z = Z, Zf = Zf, Y=Y, Yf=Yf, N = N)
# 
# ### logistic fitting?
# library(tidyverse)
# library(caret)
# library(glmnet)
# data <- data.frame(X,Y=as.factor(Y))
# log_fit <- glmnet(X,as.factor(Y), family = "binomial", lambda = 0.08)
# summary(log_fit)
# Yf_predict <- predict(log_fit,Xf, type = "class")
# length(which(Yf_predict == Yf))
# 
# cvlogfit <- cv.glmnet(X,Y, family = "binomial", type.measure = "class")
# cvYPredict <- predict(cvlogfit, Xf, s = "lambda.min", type = "class")
# length(which(cvYPredict == Yf))/length(Yf)
# 
# y <- Y
# x <- X[,161]
# df <- data.frame(x,y)
# just_volume <- glm(y ~ x, family = "binomial", data = df)
# yfvol <- predict(just_volume, newdata = data.frame(x = Xf[,161]), type = "response")
# length((which(Yf == (yfvol > 0.5))))
# 
# # 
# library(survminer)
# library(survival)
# ggsurvplot(survfit(Surv(hn_data$OStime,hn_data$OS)~1),data = hn_data)
# ggsurvplot(survfit(Surv(hn_aerts$OStime,hn_aerts$OS)~1),data = hn_aerts,  title = "Survival Curve",legend = "none",theme = theme_classic(), palette = "skyblue")
# plot(density(hn_aerts$OStime))
# 
# 
# 
# 
# big_set <- c(7,
#              12,
#              41,
#              44,
#              47,
#              52,
#              81,
#              87,
#              121,
#              124,
#              127)
# little_set <- c(52,121,127)
# train_data <- X[,big_set]
# train_data <- as.data.frame(scale(train_data))
# 
# train_data$Y <- as.factor(Y)
# 
# test_data <- Xf[,big_set]
# test_data <- as.data.frame(scale(test_data))
# test_data$Y <- as.factor(Yf)
# 
# 
# ### 
# trial_glm <- glm(Y ~ . , data = train_data,family = binomial(link = "logit"))
# new <- predict(trial_glm,newdata  = test_data, type  = "response")
# plotROC(Yf,new)
# #plot(trial_glm)
# 
# # my_predict <- predict(trial_glm,test_data, type = "response")
# # plotROC(Yf,my_predict)
# # #opt_cutoff <- optimalCutoff(Yf, my_predict)
# # misClassError(Yf, my_predict,threshold = 0.5)
# 
# library(randomForest)
# rf_predict <- randomForest(train_data[,-12],y = train_data$Y)
# rf_out <- predict(rf_predict,test_data[,-12])
# (sum(Yf == rf_out)/length(Yf))
# 
# train_data <- X[,little_set]
# train_data <- as.data.frame(scale(train_data))
# train_data$Y <- as.factor(Y)
# test_data <- Xf[,little_set]
# test_data <- as.data.frame(scale(test_data))
# test_data$Y <- as.factor(Yf)
# 
# rf_predict_small <- randomForest(train_data[,-4],y = train_data$Y)
# rf_out_small <- predict(rf_predict_small,test_data[,-4])
# (sum(Yf == rf_out_small)/length(Yf))
# 
# trial_glm <- glm(Y ~ . , data = train_data,family = binomial(link = "logit"))
# new <- predict(trial_glm,newdata  = test_data, type  = "response")
# plotROC(Yf,new)
