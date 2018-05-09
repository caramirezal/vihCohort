## Machine learning on CD4 reconstitution of VIH patients


#######################################################################
## Data preprocessing

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

##########################################################################

## lasso and ridge regression library
library(glmnet)

## performs and leave one out lasso regression
lasso <- cv.glmnet(x=input,
                   y=output)

## lasso predictions
y_predicted <- predict(lasso,
                       newx = input,
                       s = lasso$lambda.1se,
                       type="response")

## plot lasso predicted versus observed values 
plot(output,
     y_predicted,
     pch=20,
     col="steelblue",
     cex=1.5)
abline(0,1,col="red",lwd=2)

mse_lasso <- mean((y_predicted-output)^2)
mse_lasso

###########################################################################
## performing ridge regression

## perform leave one out cross validation out ridge regression
ridge <- cv.glmnet(x = input, y = output,
                   alpha = 0)

## ridge regression predictions
y_predicted_r <- predict(ridge,
                         newx = input,
                         type="response",
                         ridge$lambda.min)


## observed vs predicted ridge regression values
plot(output, y_predicted_r,
     pch=20, col="steelblue", cex=1.5)
abline(0,1,
       col="red",lwd=2)


## calculating ridge regression mse
mse_ridge <- mean((y_predicted_r-output)^2)
mse_ridge

###########################################################################
## validation set approach using 80% of the data as training set
## for lasso regression

## definition of the train set
n <- nrow(input)
train <- sample(1:n,0.8*n)

## performs 100 simulations 
lasso_sims <- rep(0,100)
for (i in 1:100) {
        ## calculates optimal lambda
        lasso <- cv.glmnet(x=input[train,],y=output[train])
        
        ## assessing fit on test data
        y_predicted_l <- predict(lasso,
                                 newx = input[-train,],
                                 s=lasso$lambda.min,
                                 type="response")
        mse_lasso <- mean((y_predicted_l-output[-train])^2)
        
        ## storing data
        lasso_sims[i] <- mse_lasso
}

hist(lasso_sims,breaks = 500)
