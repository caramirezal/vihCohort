library(glmnet)
library(caret)
library(ggplot2)

set.seed(777)

########################################################################
## Data preprocessing

data <- read.csv("data/TablaLPS_IL18_20180217.csv")

data <- data[,sapply(data, is.numeric)]
## droping redundant and non-sense variables
drop <- c("Expediente","IRIS","DeltaCD4_W52","Numero_consecutivo")
drop <- c(drop,grep("CD4",colnames(data),value = TRUE))
drop <- c(drop,grep("CD8porcen",colnames(data),value = TRUE))
## input definition
input <- data[,!(colnames(data)%in%drop)]
input <- as.matrix(input)
## output definition
output <- data$DeltaCD4_W52

#########################################################################
## LASSO implementation

## Define data partition
inTrain <- createDataPartition(output, p = .75, list = FALSE)
Train   <- input[ inTrain, ] # Training dataset for all model development
Test    <- input[-inTrain, ] # Final sample for model validation

## perform LASSO
lasso   <- cv.glmnet(x=Train, y=output[inTrain], nfolds=10,       
                     type.measure="mse")

## plotting MSE vs Lambda
jpeg("figures/MSEvsLambda.jpg")
plot(lasso)
dev.off()

## perform LASSO on train and test sets
predicted.Tr <- predict(lasso,
                        newx = Train,
                        s=lasso$lambda.min)
predicted.Te <- predict(lasso,
                        newx=Test,
                        s=lasso$lambda.min)

## plot cross validation
jpeg("figures/crossValidation.jpg")
par(mfrow=c(1,2))
plot(predicted.Tr,output[inTrain],
     pch=20,
     main = "Train",
     xlab = "input",
     ylab = "output")
abline(0,1,col="red")
plot(predicted.Te,output[-inTrain],
     pch=20,
     main = "Test",
     xlab = "input",
     ylab = "output")
abline(0,1,col="red")
par(mfrow=c(1,1))
dev.off()

## LASSO best ranked variable is correlated to output
jpeg("figures/bestRankedLASSOvariables.jpg")
g <- ggplot(data,aes(LogCV_S08,DeltaCD4_W52,colour=IRIS))
g <- g + geom_point()
g <- g + geom_smooth(method = "lm")
plot(g)
dev.off()
