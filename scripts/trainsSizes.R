## In this script lasso, ridge and random forest regressions
## are performed in different training sizes in order to select
## the best method and training size.

library(dplyr)
library(reshape2)
library(ggplot2)
library(randomForest)
library(glmnet)        ## lasso and ridge regression


## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

##########################################################################
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


##########################################################################
## initial parameter definitions

## number of observations
n <- nrow(input)
## number of simulations
n_sims <- 100

## train testing sizes
sizes <- c(0.8,0.85,0.9,0.95,1)*n


##########################################################################
## performing simulations

## matrix to store results of simulations
results <- matrix(0,length(sizes)*n_sims,7)
colnames(results) <- c("mse_lasso","r_sq_lasso",
                       "mse_ridge","r_sq_ridge",
                       "mse_forest","r_sq_forest",
                       "size")
head(results)


## iteration of lasso, ridge, and random forest regression simulations
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
                        r_sq_lasso <- r_sq(observed=output[-train],
                                           predicted=y_predicted) 
                } else {   ## sizes = 1 uses all data to predict
                        r_sq_lasso <- r_sq(observed=output[train],
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
                        r_sq_ridge <- r_sq(observed=output[-train],
                                           predicted=y_predicted_r) 
                } else {   ## sizes = 1 uses all data to predict
                        r_sq_ridge <- r_sq(observed=output[train],
                                           predicted=y_predicted_r)
                }
                
                ############################################################
                ## random forest
                ## performs random forest
                forest <- randomForest(output~.,data = input,subset = train)
                
                ## calculate prediction
                if ( sizes[i] < 1 ) {
                        y_predicted_f <- predict(forest,newdata = input[-train,])
                } else {
                        y_predicted_f <- predict(forest,newdata = input[train,])
                }
                
                ## calculate random forest mse 
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        mse_forest <- mean((y_predicted_f-output[-train])^2) 
                } else {   ## sizes = 1 uses all data to predict
                        mse_forest <- mean((y_predicted_f-output[train])^2)
                }
                
                ## random forest r^2 calculation
                if ( sizes[i] < 1 ) { ## sizes < 1 uses test data for predict   
                        r_sq_forest <- r_sq(observed=output[-train],
                                            predicted=y_predicted_f) 
                } else {   ## sizes = 1 uses all data to predict
                        r_sq_forest <- r_sq(observed=output[train],
                                            predicted=y_predicted_f)
                }
                
                
                
                #########################################################
                ## storing results in results matrix
                
                index <- (i-1)*n_sims + j
                results[index, "mse_lasso"] <- mse_lasso
                results[index,"r_sq_lasso"] <- r_sq_lasso
                results[index,"mse_ridge"] <- mse_ridge
                results[index,"r_sq_ridge"] <- r_sq_ridge
                results[index,"mse_forest"] <- mse_forest
                results[index,"r_sq_forest"] <- r_sq_forest
                results[index, "size"] <- sizes[i]
        }
}

write.csv(results,"../data/trainSizesResults.csv",row.names = FALSE)

#######################################################################
## processing results

results.processed <- as.data.frame(results)
results.processed <- mutate(results.processed,size=as.factor(size))
results.processed <- melt(results.processed)
## calculating mean of mse and R^2 values
means <- with(results.processed,
              aggregate(value,list(size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results.processed,
           aggregate(value,list(size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")
results.summary <- mutate(results.summary,size=as.character(size))
results.summary <- mutate(results.summary,size=100*(as.numeric(size)/n))
results.summary

############################################################################
## plotting results

## plotting mse values
mse_summary <- filter(results.summary,grepl("mse",technique))
write.csv(mse_summary,file = "../data/mse_summary.csv",row.names = FALSE) 

theme_set(theme_bw())

jpeg("../figures/trainSizeMSE.jpg")
ggplot(mse_summary, aes(x=size,y=log(mean),colour=technique)) +
        geom_line(size=1.8) +                                    ## not working
        geom_errorbar(aes(ymin=log(mean-sd),ymax=log(mean+sd)), width=.1) +
        geom_point(size=3.5) +
        labs(x="Size", y="Log(MSE)") +
        theme(text = element_text(size=16,face="bold"),
              axis.line = element_line(colour = 'black', size = 0.7))
dev.off()

## plotting results r^2
r_sq_summary <- filter(results.summary,grepl("r_sq",technique))
r_sq_summary

jpeg("../figures/trainSizeRsq.jpg")
ggplot(r_sq_summary, aes(x=size,y=mean,colour=technique)) + 
        geom_line(size=1.8) + 
        geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=1) + 
        geom_point(size=3.5) +
        labs(x="Size", y="R^2") +
        theme(text = element_text(size=16,face = "bold"),       ## change font text  
              axis.line = element_line(colour = 'black', size = 0.7))   ## increase axis width
dev.off()        


