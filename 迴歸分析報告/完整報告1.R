forestdata= read.table("forestfires.csv",sep = ",",header = T)


#將月份和日期字符串變量轉換為數值
forestdata$month <- as.numeric(as.factor(forestdata$month))
forestdata$day <- as.numeric(as.factor(forestdata$day))

# 檢查箱形圖中每個變量的分佈 #可不打
par(mfrow=c(1,1))
boxplot(forestdata$X,main="X")
boxplot(forestdata$Y,main ='Y')
boxplot(forestdata$FFMC, main='FFMC') #outliers
boxplot(forestdata$DMC, main ='DMC') # outliers
boxplot(forestdata$DC, main='DC') # some outliers
boxplot(forestdata$ISI,main='ISI') # outliers
boxplot(forestdata$temp, main='temp') 
boxplot(forestdata$RH,main="RH") # outliers
boxplot(forestdata$wind, main='wind') #
boxplot(forestdata$rain, main='rain')  # heavy outliers...high variability in data
boxplot(forestdata$area, main='area') # heavy outliers..high variability in data
# X，Y，DC，FFMC等變量中也觀察到不對稱性

#可視化數據
par(mfrow=c(2,6))
for (var in 1:(dim(forestdata)[2]-1)){
  d <- density(forestdata[,var])
  plot(d, main = names(forestdata[var]),xlab="")}
#Rain的右偏斜很大，FFMC的左偏斜很大，area的右偏斜率很大，需數值轉換。
#rain轉換數值，大多數值的都是零，也只有八筆資料有rain值，但rain影響沒有很大。
forestdata <- forestdata[,-which(colnames(forestdata)== "rain")] #從模型中刪除

#對area取ln(y+1)減少右偏
forestdata$area <- log(forestdata$area+1)  ###

#由於FFMC是左偏斜的，因此我們將其立方化以對其進行歸一化
forestdata$FFMC<- (forestdata$FFMC^3)   ###

#同質性在這裡我們可以看到數據滿足同質性的假設，其中模型殘差的便異數是隨機分佈的。 
#使用Breush-Pagan檢驗，我們可以看到一個不顯著的結果，表明該數據可以被認為是同質性的。
assumptionsmodel <- lm(area ~ ., data=forestdata)
bptest(assumptionsmodel)
par(mfrow=c(2,2))
plot(assumptionsmodel)

#由於樣本量較大，即使與正常值之間的偏差也很小，因此也會被認為是有意義的，因此可以通過繪圖來評估正常性。
#似乎還沒有正常分佈的殘差有很大的偏差。 這似乎是由於數據庫中存在大量0。 刪除這些0後，我們可以看到殘差變得更正態分佈。 
#儘管這意味著我們損失了大量數據案例，但這是正確使用lm模型所必需的。
assumptionsmodel_all <- lm(area ~ ., data=forestdata)
assumptionsmodel_0 <- lm(area ~ .,data=forestdata[which(forestdata$area>0),]) #刪除所有面積為0
par(mfrow=c(1,2))
hist(assumptionsmodel_all$residuals, main = "areau有0的", xlab = 'Residuals')
abline(v=mean(assumptionsmodel_all$residuals), col='red', lwd=2)
hist(assumptionsmodel_0$residuals,main = "areau沒有0的", xlab = 'Residuals')
abline(v=mean(assumptionsmodel_0$residuals), col='red', lwd=2)
forestdata = forestdata[which(forestdata$area>0),] ###


#基本模型
m <- lm(area ~ ., data=forestdata)
summary(m) 


#-------------------------------------------------------------------------------------------------#
#選擇變數
none = lm(area~1,forestdata)
full = lm(area~.,forestdata)
step(none,scope = list(upper=full,lower=none),direction = 'forward')  # month
step(full,scope = list(upper=full,lower=none),direction = 'backward') # month DMC DC ISI wind RH
step(none,scope = list(upper=full,lower=none),direction = 'both')     # month

library(leaps)
r = leaps(forestdata[,1:11],forestdata$area,method = "adjr2") # month DMC DC ISI RH wind 
r$which[which.max(r$adjr2),]

# model1
m1 = lm(area~month+DMC+DC+ISI+RH+wind,data = forestdata) ; summary(m1)
m2=lm(area~month,data=forestdata);summary(m2)
anova(m1,m2) 

#相關性(交互作用)
par(mfrow=c(1,1))
corrplot(cor(forestdata), method="number", outline = TRUE,type="upper")
# FFMC,ISI=0.74 ， FFMC.temp=0.6 ， DMC.DC=0.5 ，  temp.RH=-0.5 

forestdata$DMC.DC = (forestdata$DMC)*(forestdata$DC)
forestdata$FFMC.ISI = (forestdata$FFMC)*(forestdata$ISI)
forestdata$temp.RH = (forestdata$temp)*(forestdata$RH)
forestdata$FFMC.temp = (forestdata$FFMC)*(forestdata$temp)


###
none = lm(area~1,forestdata)
full = lm(area~.,forestdata) #month DMC DC ISI wind RH
step(none,scope = list(upper=full,lower=none),direction = 'forward')  # month
step(full,scope = list(upper=full,lower=none),direction = 'backward') # month DMC ISI wind DMC*DC temp*RH
step(none,scope = list(upper=full,lower=none),direction = 'both')     # month

r = leaps(forestdata[,-12],forestdata$area,method = "adjr2") # month DMC ISI wind DMC*DC FFMC*ISI
r$which[which.max(r$adjr2),]

# 交互model
m2 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI , data = forestdata) ; summary(m2)
# r2=0.05366，adjr2=0.02465
m3 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI+DMC*DC , data = forestdata) ; summary(m3)
# r2=0.05797，adjr2=0.02536
m4 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI+DMC*DC+temp*RH , data = forestdata) ; summary(m4)
# r2=0.06117，adjr2=0.02114

extractAIC(m2)
extractAIC(m3)
extractAIC(m4)
AIC(m2)
AIC(m3)
AIC(m4)



#-------------------------------------------------------------------------------------------------#
#離群值
Ei = m$residuals
hii = hatvalues(m)
ri = rstandard(m)
ti = rstudent(m)
#影響值
dffitsi = dffits(m)
Di = cooks.distance(m)
dfbetasi = dfbetas(m)

n=270 ;p=11 ;apha=0.05
case1=which(hii>2*p/n);names(case1)=NULL;case1
case2=which(abs(ri)>2);names(case2)=NULL;case2
case3=which(abs(ti)>qt(1-(0.05/2),df=n-p-1,lower.tail = FALSE));names(case3)=NULL;case3
intersect(case1,case2)
intersect(case1,case3) #離群值: 269

case4=which(abs(dffitsi)>2*sqrt(p/n));names(case1)=NULL;case4
case5=which(Di>4/n);names(case2)=NULL;case5
case6=which(abs(dfbetasi)>2/sqrt(n),arr.ind=T);names(case3)=NULL;case6
intersect(case4,case5) #影響值: 269

influencePlot(m,id.n=5)
forestdata2 = forestdata[-269,] #從數據中刪除269


m5 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI , data = forestdata2) ; summary(m5)
# r2=0.05789，adjr2=0.0289
m6 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI+DMC*DC , data = forestdata2) ; summary(m6)
# r2=0.06138，adjr2=0.02876
m7 = lm(area~month+DMC+DC+ISI+RH+wind+FFMC*ISI+DMC*DC+temp*RH , data = forestdata2) ; summary(m7)
# r2=0.06824，adjr2=0.02836

par(mfrow=c(1,1))
cv_all<-cv.lm(data=forestdata2, m7, m=3) 
RMSD_origin<-sqrt(mean((cv_all$cvpred-cv_all$area)^2))
print(paste("All variable model RMSD: ",RMSD_origin))

#下圖顯示了具有不同變量子集的模型如何根據調整後的R2進行比較，其中我們可以在頂部看到最佳子集的變量。
subset.out<-regsubsets(area~ ., data=forestdata2,nbest=1,nvmax=NULL,method="exhaustive")
summary.out<-summary(subset.out)
# 返回最佳模型中使用的變量
bestmodelvariables<- summary.out$which[which.max(summary.out$adjr2),]  
par(mfrow=c(1,1))
# 根據調整後的r2值放所有子集
plot(subset.out,scale="adjr2",main="All Subset Regression")

variables <- names(bestmodelvariables[which(bestmodelvariables==TRUE)])[-1] ;variables
#我們可以看到最佳模型使用以下變量：month, DMC, ISI, wind, DMC.DC, FFMC.ISI, temp.RH。


formula <- noquote(paste("area~",paste(variables, collapse = '+'),sep=""))
modelfinal<-lm(formula,data=forestdata2)  ;summary(modelfinal)
# r2=0.0631，adjr2=0.038




























