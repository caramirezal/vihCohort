library(dplyr)
library(glmnet)
library(caret)

########################################################################
## initial exploratory lasso analysis of raw data

data <- read.csv("data/TablaLPS_IL18_20180217.csv")

data <- data[,sapply(data, is.numeric)]
data <- select(data,-Expediente)
data <- select(data,-IRIS)
input <- as.matrix(select(data,-DeltaCD4_W52))
response <- data$DeltaCD4_W52

## lasso analysis
lasso <- cv.glmnet(input,response)
plot(lasso)

## model fit using lasso
pred.vals <- predict(lasso,
                     newx = input,
                     s = lasso$lambda.min)

plot(response,pred.vals)
abline(0,1,col="red")
coef(lasso)


######################################################################
## drop any variable related to CD4 values

drop <- grepl("CD4",colnames(input))
drop <- colnames(input)[drop]
#drop <- c("CD4_S0","CD4_S52","DeltaCD4W52atleast150")
input <- input[,!(colnames(input)%in%drop)]

lasso <- cv.glmnet(input,response)
cs <- as.vector(coef(lasso))
indexes <- 2:length(cs)
cs <- cs[indexes]
names(cs) <- rownames(coef(lasso))[indexes]
sort(cs,decreasing = TRUE)

coefs <- names(cs[cs!=0])
coefs

############################################################################

## Now use: Training and testing sets
# Split dataframe into training & testing sets
inTrain <- createDataPartition(response, p = .75, list = FALSE)
Train   <- input[ inTrain, ] # Training dataset for all model development
Test    <- input[-inTrain, ] # Final sample for model validation

lasso   <- cv.glmnet(x=Train, y=response[inTrain], nfolds=10, type.measure="mse")
plot(lasso)
coef(lasso)



predicted.Tr <- predict(lasso,
                       newx = Train,
                       s=lasso$lambda.min)
predicted.Te <- predict(lasso,
                        newx=Test,
                        s=lasso$lambda.min)

par(mfrow=c(1,2))
plot(predicted.Tr,response[inTrain],main = "Train")
abline(0,1,col="red")
plot(predicted.Te,response[-inTrain],main = "Test")
abline(0,1,col="red")
par(mfrow=c(1,1))

mean((response[inTrain]-predicted.Tr)^2)
mean((response[-inTrain]-predicted.Te)^2)


lasso.vis <- data.frame(predicted.Tr,response[inTrain],iris[inTrain])
names(lasso.vis) <- c("x","y","outcome")
library(ggplot2)

fit <- predict(lasso,newx = input,s=lasso$lambda.min)
plot(data$LogCV_S08,data$DeltaCD4_W52)
points(data$LogCV_S08,fit,
      xlim=c(min(data$LogCV_S08),max(data$LogCV_S08)),
      ylim=c(min(data$DeltaCD4_W52),max(data$LogCV_S08)),
      col="red")

plot(data$LPS_pg_ml_W0,data$DeltaCD4_W52+)

lm(DeltaCD4_W52~LPS_pg_ml_W0+LogCV_S08+IL18_pg_ml+TNFalfa_pg_ml_W0+CD8porcentajeS12,data)
