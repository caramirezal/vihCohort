## Evaluation of variation in the coefficients for lasso

library(glmnet)
library(caret)

#########################################################################
## data preprocessing

## reading data
data <- read.csv("~/scripts/vihCohort/data/TablaLPS_IL18_20180217.csv")
head(data)


## preprocessing raw data
## drop categorical data
data <- data[,sapply(data, is.numeric)]

## list of removed variables
drop <- c("Expediente",
          "IRIS",
          "DeltaCD4_W52",
          "Numero_consecutivo",
          "CD8porcentajeS24",
          "CD4_S52")

## removing list of variables
input <- data[,!(colnames(data)%in%drop)]

## input definition
input <- as.matrix(input)

## output definition
output <- data$DeltaCD4_W52
colnames(input)

#########################################################################
## lasso iteration

## model parameters
n.sim <- 100
ratio <- 0.8
nfolds <- 5


## vector to store coeficients
res <- numeric(ncol(input)+1)

## lasso iterations
for (i in 1:n.sim) {
        ## definition of the train and test data
        inTrain <- createDataPartition(output, p = ratio, list = FALSE)
        Train   <- input[ inTrain, ] 
        
        ## performing lasso
        lasso   <- cv.glmnet(x=Train, y=output[inTrain], nfolds=nfolds,       
                             type.measure="mse")
        
        ## stroing lasso results
        lasso.coef <- coef(lasso)
        lasso.names <- rownames(lasso.coef)
        lasso.coef <- as.vector(lasso.coef)
        names(lasso.coef) <- lasso.names
        res <- res + lasso.coef
}
## calculate average
res <- res/n.sim
## drop intercept coeficiente
res <- res[res >0 & ! names(res) %in% "(Intercept)"]
sort(res,decreasing = TRUE)