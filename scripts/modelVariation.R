## Evaluation of variation in the coefficients for lasso
library(glmnet)
library(dplyr)

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
ratio <- 0.9
nfolds <- 5


## vector to store coeficients
res <- matrix(0,n.sim,ncol(input)+1)

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
        lasso.coef <- as.vector(lasso.coef)
        res[i,] <- lasso.coef
}

colnames(res) <- rownames(coef(lasso))

## saving data
write.csv(res,file = "../data/lassoVariation.csv",row.names=FALSE)
sims <- read.csv("../data/lassoVariation.csv")

## getting top 10 coefficients
sims.summary <- data.frame("mean_lasso"=apply(sims,2,mean),
                           "sd_lasso"=apply(sims,2,sd))
sims.summary <- sims.summary[rownames(sims.summary)!="X.Intercept.",]
sims.summary[order(-abs(sims.summary$means)),][1:10,]
