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
                   alpha = 0)   ## alpha must be set to 0 in rr

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
train <- sample(1:n,n*0.8)

## performs 100 simulations 
lasso_sims <- rep(0,100)
for (i in 1:5) {
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

#######################################################################
## Testing multiple train set sizes

## calculate r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

## number of observations
n <- nrow(input)
## number of simulations
n_sims <- 5

## train testing sizes
sizes <- c(0.8,0.9,0.95,1)*n

## matrix to store results of simulations
results <- matrix(0,length(sizes)*n_sims,5)
colnames(results) <- c("mse_lasso","r_sq_lasso",
                       "mse_ridge","r_sq_ridge",
                       "size")
head(results)


## iteration of lasso and ridge regression simulations
for ( i in 1:length(sizes) ) {
        ## train test definition
        train <- sample(1:n,sizes[i])
        for ( j in 1:n_sims ) {
                ##########################################################
                ## lasso 
                
                ## lasso Lambda optimal calculation
                lasso <- cv.glmnet(x=input[train,],y = output[train])
                
                ## calculate lasso predictions
                if ( sizes[i] < 1 ) {   ## sizes < 1 uses test data
                        y_predicted <- predict(lasso, newx=input[-train,],
                                               c=lasso$lambda.min,
                                               type="response")
                } else {   ## sizes = 1 uses all data
                        y_predicted <- predict(lasso, newx=input[train,],
                                               c=lasso$lambda.min,
                                               type="response")
                }
                
                ## lasso mse calculation
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        mse_lasso <- mean((y_predicted-output[-train])^2) 
                } else {   ## sizes = 1 uses all data to predict
                        mse_lasso <- mean((y_predicted-output[train])^2)
                }
                
                ## lasso r^2 calculation
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        r_sq_lasso <- rsq(observed=output[-train],
                                          predicted=y_predicted) 
                } else {   ## sizes = 1 uses all data to predict
                        r_sq_lasso <- rsq(observed=output[train],
                                          predicted=y_predicted)
                }
                
                #########################################################
                ## Ridge regression
                ridge <- cv.glmnet(x=input[train,],y=output[train],
                                   alpha=0)
                
                ## calculate ridge regression
                if ( sizes[i] < 1 ) {   ## sizes < 1 uses test data
                        y_predicted_r <- predict(ridge, newx=input[-train,],
                                               c=lasso$lambda.min,
                                               type="response")
                } else {   ## sizes = 1 uses all data
                        y_predicted_r <- predict(ridge, newx=input[train,],
                                               c=lasso$lambda.min,
                                               type="response")
                }
                
                
                ## ridge regression mse 
                ## lasso mse calculation
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        mse_ridge <- mean((y_predicted_r-output[-train])^2) 
                } else {   ## sizes = 1 uses all data to predict
                        mse_ridge <- mean((y_predicted_r-output[train])^2)
                }
                
                ## Rigde regression r^2 calculation
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        r_sq_ridge <- rsq(observed=output[-train],
                                          predicted=y_predicted_r) 
                } else {   ## sizes = 1 uses all data to predict
                        r_sq_ridge <- rsq(observed=output[train],
                                          predicted=y_predicted_r)
                }
                
                #########################################################
                ## storing results in results matrix
                
                index <- (i-1)*n_sims + j
                results[index, "mse_lasso"] <- mse_lasso
                results[index,"r_sq_lasso"] <- r_sq_lasso
                results[index,"mse_ridge"] <- mse_ridge
                results[index,"r_sq_ridge"] <- r_sq_ridge
                results[index, "size"] <- sizes[i]
        }
}


