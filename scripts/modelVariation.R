## Evaluation of variation in the coefficients for lasso
library(glmnet)
library(randomForest)
library(dplyr)
library(caret)

#########################################################################
## data preprocessing

## reading data
data <- read.csv("../data/TablaLPS_IL18_20180217.csv")
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


## matrix to store coeficients
res <- matrix(0,n.sim,ncol(input)+1)

## vector to store lambda
lambda_min <- rep(0,n.sim)
        
        
## lasso iterations
for (i in 1:n.sim) {
        ## definition of the train and test data
        inTrain <- createDataPartition(output, p = ratio, list = FALSE)
        Train   <- input[ inTrain, ] 
        
        ## performing lasso
        lasso   <- cv.glmnet(x=Train, y=output[inTrain],       
                             type.measure="mse")
        
        ## stroing lasso results
        lasso.coef <- coef(lasso)
        lasso.coef <- as.vector(lasso.coef)
        res[i,] <- lasso.coef
        lambda_min[i] <- lasso$lambda.min
}

colnames(res) <- rownames(coef(lasso))

## saving data
#write.csv(res,file = "../data/lassoVariation.csv",row.names=FALSE)
sims <- read.csv("../data/lassoVariation.csv")

## getting top 10 coefficients
sims.summary <- data.frame("mean_lasso"=apply(sims,2,mean),
                           "sd_lasso"=apply(sims,2,sd))
sims.summary <- sims.summary[rownames(sims.summary)!="X.Intercept.",]
sims.summary <- mutate(sims.summary,variable_lasso=row.names(sims.summary))
head(sims.summary)

## getting biggest coefficients in absolute value 
lasso_top <-sims.summary[order(-abs(sims.summary$mean_lasso)),]
lasso_top <- mutate(lasso_top,ranking=1:nrow(lasso_top))
head(lasso_top)

##############################################################################
## Random forest iteration

results_forest <- matrix(0,n.sim,ncol(input))
colnames(results_forest) <- rownames(importance(forest))
results_forest

for (i in 1:n.sim) {
        ## create train partition
        inTrain <- createDataPartition(output, p = ratio, list = FALSE)
        Train   <- input[ inTrain, ]
        
        ## perform random forest
        forest <- randomForest(x = input,y = output,importance = TRUE)
        
        ## store data
        results_forest[i,] <- importance(forest)[,"%IncMSE"]
}

# saving data
#write.csv(results_forest,"../data/importanceForestSims.csv",row.names = FALSE)
results_forest <- read.csv("../data/importanceForestSims.csv")

forest_mean <- apply(results_forest,2,mean)
forest_mean

forest_sd <- apply(results_forest,2,sd)
forest_sd

forest.summary <- data.frame("variable_forest"=names(forest_mean),
                             "mean_forest"=forest_mean,
                             "sd_forest"=forest_sd)
forest.summary <- arrange(forest.summary,desc(mean_forest))
forest.summary <- mutate(forest.summary,ranking=1:nrow(forest.summary))
head(forest.summary)

final_summary <- merge(lasso_top,forest.summary,by = "ranking")
final_summary <- final_summary[1:15,c(1,4,2,3,5,6,7)]
colnames(final_summary) <- c("ranking",
                             "lasso_variable",
                             "coefficient",
                             "sd",
                             "variable_forest",
                             "MSE_increase",
                             "sd")
write.csv(final_summary,"../data/ranking.csv",row.names = FALSE)
