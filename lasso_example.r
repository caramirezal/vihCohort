library(lars)
library(glmnet)
library(caret)

data <- read.csv("C08_04paraIL10etc180123.csv",header=TRUE)
dim(data)
head(data)

Y1 <- data$DeltaCD4_W52
X1 <- cbind(data$CD4_S0,data$LPS_pg_ml_W0,data$IL18_pg_ml,data$LeucocitosS0,data$HematocritoS0,
            data$PlaquetasS0,data$GlucosaS0,data$Acido_uricoS0,data$CD4porcentajeS0,data$CD8_S0,data$CD8porcentajeS0,data$CocCD4_CD8_S0,data$LogCV_S0,data$CD4_S12,data$DeltaCV_S8)

# Simple function lasso regression
reg1 <- cv.glmnet(X1,Y1)
plot(reg1)

best_lambda1 <- reg1$lambda.min
fit1 <- predict(reg1,newx=X1, s="lambda.min")

plot(fit1,Y1)
abline(0,1,col="red")
r2 <- cov(fit1,Y1)/(sd(fit1)*sd(Y1)) # rsquared

## Now use: Training and testing sets
# Split dataframe into training & testing sets
inTrain <- createDataPartition(Y1, p = .75, list = FALSE)
Train   <- X1[ inTrain, ] # Training dataset for all model development
Test    <- X1[-inTrain, ] # Final sample for model validation

CV1   <- cv.glmnet(x=Train, y=Y1[inTrain], nfolds=10, type.measure="mse")

# Prediction and mse of Training set
Train_pred1 <- predict(CV1,type="response", newx = Train, s = CV1$lambda)
dim(Train_pred1)

Train_coef1 <- predict(CV1,type="coefficients", newx = Train, s = CV1$lambda)
dim(Train_coef1)

tmse <- apply((Train_pred1-Y1[inTrain])^2,2,mean) # same as mse in Metrics package

# Prediction and mse of Test set
Test_pred1 <- predict(CV1,type="response", newx = Test, s = CV1$lambda)
mse        <- apply((Test_pred1-Y1[-inTrain])^2,2,mean)

# Prediction of Training and Testing set at lambda min
pdf("scatter1.pdf")
Train_pred1lmin <- predict(CV1,type="response", newx = Train, s = "lambda.min")
plot(Train_pred1lmin,Y1[inTrain])
DF <- as.data.frame(cbind(Train_pred1lmin,Y1[inTrain]))
colnames(DF) <- c("Pr","Ob")
abline(fit <- lm(Pr ~ Ob, data=DF), col='red')
legend("topright",legend=paste("R2","=",round(cov(DF$Pr,DF$Ob)/(sd(DF$Pr)*sd(DF$Ob)),digits=2),sep=" "),lwd=1,col="red")
dev.off()

pdf("Model_1.pdf")
plot(log(CV1$lambda),tmse,col="blue",ylim=c(min(c(mse,tmse)),max(c(mse,tmse))),ylab="Mean-Squared Error",xlab="log(lambda)",main="Model 1")
par(new=TRUE)
plot(log(CV1$lambda),mse,col="red",ylim=c(min(c(mse,tmse)),max(c(mse,tmse))),axes=FALSE,ylab="",xlab="")
dev.off()

