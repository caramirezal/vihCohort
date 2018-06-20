## This script contains the R functions to implement Machine Learning, parameter selection,
## and plotting results in technicalReportML.Rmd 

## Dependencies
library(randomForest)
library(caret)
library(party)
library(gridExtra)
library(reshape2)
library(glmnet)
library(dplyr)
library(e1071)

## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

#########################################################################################
## data preprocessing

## reading data
data <- read.csv("../data/TablaLPS_IL18_20180217.csv")

## drop categorical data
data <- data[,sapply(data, is.numeric)]

## list of removed variables
drop <- c("Expediente",
          "DeltaCD4W52atleast150",
          "DeltaCD4_W52",
          "CD4TCabove200_W052",
          "CD4TCabove350_W052",
          "Numero_consecutivo",
          "CD4_S52")

## removing list of variables
input <- data[,!(colnames(data)%in%drop)]

## input definition
input <- as.matrix(input)

## input to dataframe
input.df <- as.data.frame(input)

## output definition
output <- data$DeltaCD4_W52



##################################################################################
## randomForest evaluation                                                      ##
##################################################################################

## ntree parameter variation in randomForest method

## initial parameters
nsims <- 100
mtryInterval <- seq(2,10,1)

## results is defined data frame to store results
results <- matrix(0,nsims*length(mtryInterval),3)
results <- as.data.frame(results)
colnames(results) <- c("Accuracy","mtry","performance")

for (i in 1:length(mtryInterval)){
        for (j in 1:nsims) {
                
                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                
                ## perform randomForest
                initialTime <- Sys.time()
                rforest <- randomForest(output[intrain]~.,
                                        data=input[intrain,],
                                        mtry=mtryInterval[i]) ## parameter variation
                finalTime <- Sys.time() - initialTime
                
                ## perform predictions
                rforest.pred <- predict(rforest,
                                        newdata = input[-intrain,])
                
                ## calcualte mse values
                mss <- mean((rforest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                mse
                
                ## index to store simulation data
                index <- (nsims*(i-1)) + j
                
                ## storing results
                results[index,1] <- mse
                results[index,2] <- mtryInterval[i]
                results[index,3] <- finalTime
        } 
}

write.csv(results,"../data/rforestMtryVariation.csv",row.names = FALSE)
#results <- read.csv("../figures/rforestMtryVariation.csv")

## plot accuracy test
theme_set(theme_bw())
g <- ggplot(results, aes(x=mtry,y=Accuracy)) + geom_point() 
g <- g + geom_boxplot(aes(group=mtry,fill="steelblue"))
g <- g + theme(legend.position = "none",      ## remove legends
               text = element_text(size=26,face="bold")) 
g <- g + labs(y="MSE")
g <- g + theme(title = element_text(face="bold"),
               axis.title.x = element_text(face="bold"),
               axis.title.y = element_text(face="bold"))

## plot performance time
theme_set(theme_bw())
g1 <- ggplot(results, aes(x=mtry,y=performance)) + geom_point() 
g1 <- g1 + geom_boxplot(aes(group=mtry,       ## x axis cathegorical values  
                            fill=132))        ## 132 -> steelblue color
g1 <- g1 + theme(legend.position = "none",
                 text = element_text(size=26,face="bold"))    
g1 <- g1 + labs(y="Time (secs)")
g1 <- g1 + theme(title = element_text(face="bold"),  
                 axis.title.x = element_text(face="bold"),
                 axis.title.y = element_text(face="bold"))


grid.arrange(g,g1,nrow=1)


#######################################################################################
## ntree parameter variation in randomForest method

## initial parameters
nsims <- 100
ntreeInterval <- seq(50,500,50)

## data to store results
results <- matrix(0,nsims*length(ntreeInterval),3)
results <- as.data.frame(results)
colnames(results) <- c("Accuracy","ntree","performance")

for (i in 1:length(ntreeInterval)){
        for (j in 1:nsims) {
                
                ## perform randomForest
                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                initialTime <- Sys.time()
                forest <- randomForest(output[intrain]~.,
                                       data=input[intrain,],
                                       ntree = ntreeInterval[i]) ## parameter variation
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input[-intrain,])
                mss <- mean((forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## index to store simulation data
                index <- (nsims*(i-1)) + j
                
                ## storing results
                results[index,1] <- mse
                results[index,2] <- ntreeInterval[i]
                results[index,3] <- finalTime
        } 
}

write.csv(results,"../data/rforestNtreeVariation.csv",row.names = FALSE)
#results <- read.csv("../data/rforestNtreeVariation.csv")

## plot Accuracy MSE values
theme_set(theme_bw())
g <- ggplot(results, aes(x=ntree,y=Accuracy)) + geom_point() 
g <- g + geom_boxplot(aes(group=ntree,fill="steelblue"))
g <- g + theme(legend.position = "none",
               text = element_text(size=26,face="bold")) 
g <- g + labs(y="MSE")
g <- g + theme(title = element_text(face="bold"),
               axis.title.x = element_text(face="bold"),
               axis.title.y = element_text(face="bold"))

## plot waiting times
g1 <- ggplot(results, aes(x=ntree,y=performance)) + geom_point() 
g1 <- g1 + geom_boxplot(aes(group=ntree,       ## x axis cathegorical values  
                            fill=132))        ## 132 -> steelblue color
g1 <- g1 + theme(legend.position = "none",
                 text = element_text(size=26,face="bold"))    ## removing legend
g1 <- g1 + labs(y="Time (secs)")
g1 <- g1 + theme(title = element_text(face="bold"),  
                 axis.title.x = element_text(face="bold"),
                 axis.title.y = element_text(face="bold"))


grid.arrange(g,g1,nrow=1)



############################################################################################
## cforest method                                                                         ##
############################################################################################

## mtry parameter variation in cforest method

## initial parameters
nsims <- 100
mtryInterval <- seq(2,10,1)

## results is defined data frame to store results
results <- matrix(0,nsims*length(mtryInterval),3)
results <- as.data.frame(results)
colnames(results) <- c("Accuracy","mtry","performance")

for (i in 1:length(mtryInterval)){
        for (j in 1:nsims) {
                
                ## simulation parameters
                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                cforestControl <- cforest_control(teststat = "quad",
                                                  testtype = "Univ", 
                                                  mincriterion = 0, 
                                                  ntree = 50, 
                                                  mtry = mtryInterval[i], ## variation
                                                  replace = FALSE)
                
                
                ## perform simulation
                initialTime <- Sys.time()
                forest <- cforest(output[intrain]~.,
                                  data=input.df[intrain,],
                                  control = cforestControl)
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input.df[-intrain,])
                mss <- mean((forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## index to store simulation data
                index <- (nsims*(i-1)) + j
                
                ## storing results
                results[index,1] <- mse
                results[index,2] <- mtryInterval[i]
                results[index,3] <- finalTime
        } 
}

#write.csv(results,"../data/cforestMtryVariation.csv",row.names = FALSE)
results <- read.csv("../data/cforestMtryVariation.csv")

## plot accuracy test
theme_set(theme_bw())
g <- ggplot(results, aes(x=mtry,y=Accuracy)) + geom_point() 
g <- g + geom_boxplot(aes(group=mtry,fill="steelblue"))
g <- g + theme(legend.position = "none",      ## remove legends
               text = element_text(size=26,face="bold")) 
g <- g + labs(y="MSE")
g <- g + theme(title = element_text(face="bold"),
               axis.title.x = element_text(face="bold"),
               axis.title.y = element_text(face="bold"))

## plot performance time
theme_set(theme_bw())
g1 <- ggplot(results, aes(x=mtry,y=performance)) + geom_point() 
g1 <- g1 + geom_boxplot(aes(group=mtry,       ## x axis cathegorical values  
                            fill=132))        ## 132 -> steelblue color
g1 <- g1 + theme(legend.position = "none",
                 text = element_text(size=26,face="bold"))    
g1 <- g1 + labs(y="Time (secs)")
g1 <- g1 + theme(title = element_text(face="bold"),  
                 axis.title.x = element_text(face="bold"),
                 axis.title.y = element_text(face="bold"))

grid.arrange(g,g1,nrow=1)


############################################################################################
## ntree parameter variation in cforest

## initial parameters
nsims <- 100
ntreeInterval <- seq(50,500,50)

results <- matrix(0,nsims*length(ntreeInterval),3)
results <- as.data.frame(results)
colnames(results) <- c("Accuracy","ntree","performance")

for (i in 1:length(ntreeInterval)){
        for (j in 1:nsims) {
                
                ## initial settings
                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                cforestControl <- cforest_control(teststat = "quad",
                                                  testtype = "Univ", 
                                                  mincriterion = 0, 
                                                  ntree = ntreeInterval[i], 
                                                  mtry = 3,
                                                  replace = FALSE)
                
                ## perform cforest
                initialTime <- Sys.time()
                forest <- cforest(output[intrain]~.,data=input.df[intrain,],
                                  control = cforestControl)
                finalTime <- Sys.time() - initialTime
                
                ## calculate predictions
                forest.pred <- predict(forest,newdata=input.df[-intrain,])
                
                ## calculate mse
                mss <- mean((forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## index to store simulation data
                index <- (nsims*(i-1)) + j
                
                ## storing results
                results[index,1] <- mse
                results[index,2] <- ntreeInterval[i]
                results[index,3] <- finalTime
        } 
}

write.csv(results,"../data/cforestNtreeVariation.csv",row.names = FALSE)
#results <- read.csv("../data/cforestNtreeVariation.csv")

## plot Accuracy MSE values
g <- ggplot(results, aes(x=ntree,y=Accuracy)) + geom_point() 
g <- g + geom_boxplot(aes(group=ntree,fill="steelblue"))
g <- g + theme(legend.position = "none",
               text = element_text(size = 26,face = "bold")) 
g <- g + labs(y="MSE")
g <- g + theme(title = element_text(face="bold"),
               axis.title.x = element_text(face="bold"),
               axis.title.y = element_text(face="bold"))

## plot waiting times
g1 <- ggplot(results, aes(x=ntree,y=performance)) + geom_point() 
g1 <- g1 + geom_boxplot(aes(group=ntree,       ## x axis cathegorical values  
                            fill=132))        ## 132 -> steelblue color
g1 <- g1 + theme(legend.position = "none",
                 text = element_text(size = 26,face = "bold")) 
g1 <- g1 + labs(y="Time (secs)")
g1 <- g1 + theme(title = element_text(face="bold"),  
                 axis.title.x = element_text(face="bold"),
                 axis.title.y = element_text(face="bold"))

grid.arrange(g,g1,nrow=1)


######################################################################################
## Compare cforest versus randomForest functions for different train
## sizes

nsims <- 100
p = seq(0.5,0.9,0.05)
input.df <- as.data.frame(input)
results <- matrix(0,nsims*length(p),7)
results <- as.data.frame(results)
colnames(results) <- c("Size",
                       "cforest_accuracy",
                       "cforest_rsq",
                       "cforest_time",
                       "rforest_accuracy",
                       "rforest_rsq",
                       "rforest_time")

for ( i in 1:length(p) ) {
        for ( j in 1:nsims ) {
                
                ###########################################################
                ## cforest method
                
                ## simulation parameters
                intrain <- createDataPartition(output,p = p[i],list = FALSE)
                cforestControl <- cforest_control(teststat = "quad",
                                                  testtype = "Univ", 
                                                  mincriterion = 0, 
                                                  ntree = 300, 
                                                  mtry = 6,
                                                  replace = FALSE)
                
                
                ## perform simulation
                initialTime <- Sys.time()
                forest <- cforest(output[intrain]~.,
                                  data=input.df[intrain,],
                                  control = cforestControl)
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input.df[-intrain,])
                mss <- mean((forest.pred-output[-intrain])^2)
                msecforest <- sqrt(mss)
                
                ## calculation of R^2
                rsqcforest <- r_sq(output[-intrain],forest.pred)
                
                index <- (i - 1)*nsims + j
                results$Size[index] <- p[i]
                results$cforest_accuracy[index] <- msecforest
                results$cforest_rsq[index] <- rsqcforest
                results$cforest_time[index] <- finalTime
                
                ########################################################
                ## randomForest method
                initialTime <- Sys.time()
                #intrain <- createDataPartition(output,p = p[i],list = FALSE)
                forest <- randomForest(output[intrain]~.,
                                       data=input[intrain,],
                                       ntree = 400,
                                       mtry=8)
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input.df[-intrain,])
                
                ## calculate accuracy
                mss <- mean((forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## calculate rsq
                rsq <- r_sq(output[-intrain],forest.pred)
                
                ## store results
                results$rforest_accuracy[index] <- mse
                results$rforest_rsq[index] <- rsq
                results$rforest_time[index] <- finalTime
                
        }
}

## stroing results
write.csv(results,"../data/forestTrainSizes.csv",row.names = FALSE)
#results <- read.csv("../data/forestTrainSizes.csv")

## preprocessing previuos results
results <- read.csv("../data/forestTrainSizes.csv")
results <- melt(results,id="Size")

## calculating mean of mse and R^2 values
means <- with(results,aggregate(value,list(Size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results, aggregate(value,list(Size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")


accuracy <- results.summary[grepl("accuracy",results.summary$technique),]


## plotting accuracy
g <- ggplot(accuracy, aes(x=size,y=mean,colour=technique)) 
g <- g +  geom_line()
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g <- g + geom_point(size=3.5)
g <- g + labs(x="Size", y="MSE")
g <- g + theme(text = element_text(size=26,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + theme(legend.position = "none")


time <- results.summary[grepl("time",results.summary$technique),]

## plotting performance time
g1 <- ggplot(time, aes(x=size,y=mean,colour=technique)) 
g1 <- g1 +  geom_line()
g1 <- g1 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g1 <- g1 + geom_point(size=3.5)
g1 <- g1 + labs(x="Size", y="Time (secs)")
g1 <- g1 + theme(text = element_text(size=26,face="bold"),
                 axis.line = element_line(colour = 'black', size = 0.7))
g1 <- g1 + theme(legend.position = "none")

## saving plot
grid.arrange(g,g1,nrow=1)


##################################################################################
##############################################################################
## Random forest iteration
## evaluation of variabibility of importance of variables
## in randomForest method

n.sim <- 100

results_forest <- matrix(0,n.sim,ncol(input))
colnames(results_forest) <- colnames(input)
results_forest

for (i in 1:n.sim) {
        ## create train partition
        intrain <- createDataPartition(output, p = 0.9, list = FALSE)
        
        ## perform random forest
        forest <- randomForest(x = input.df[intrain,],
                               y = output[intrain],importance = TRUE,
                               ntree = 150,
                               mtry = 7)
        
        ## store data
        results_forest[i,] <- importance(forest,scale=FALSE)[,"%IncMSE"]
}

# saving data
write.csv(results_forest,"../data/importancerForest.csv",row.names = FALSE)
#results_forest <- read.csv("../data/importanceForestSims.csv")

## calculating mean and sd 
forest_mean <- apply(results_forest,2,mean)
forest_sd <- apply(results_forest,2,sd)
forest.summary <- data.frame("variable_forest"=names(forest_mean),
                             "mean_forest"=forest_mean,
                             "sd_forest"=forest_sd)
forest.summary <- arrange(forest.summary,desc(mean_forest))
forest.summary <- mutate(forest.summary,ranking=1:nrow(forest.summary))

write.csv(forest.summary,"../data/rankingRForest.csv")


##############################################################################
## Random forest iteration
## evaluation of variabiblity of importance of variables
## in randomForest method

n.sim <- 100

results_cforest <- matrix(0,n.sim,ncol(input))
colnames(results_cforest) <- colnames(input)

for (i in 1:n.sim) {
        ## initial parameters
        intrain <- createDataPartition(output, p = 0.9, list = FALSE)
        cforestControl <- cforest_control(teststat = "quad",
                                          testtype = "Univ", 
                                          mincriterion = 0, 
                                          ntree = 100, 
                                          mtry = 6,
                                          replace = FALSE)
        
        ## perform random forest
        forest <- cforest(output[intrain]~.,
                          data = input.df[intrain,])
        
        ## store data
        results_cforest[i,] <- varimp(forest)
}

# saving data
write.csv(results_forest,"../data/importancecForest.csv",row.names = FALSE)
#results_forest <- read.csv("../data/importanceForestSims.csv")

## calculating mean and sd 
cforest_mean <- apply(results_cforest,2,mean)
cforest_sd <- apply(results_cforest,2,sd)
cforest.summary <- data.frame("variable_cforest"= names(cforest_mean),
                             "mean_cforest"= cforest_mean,
                             "sd_cforest"= cforest_sd)
cforest.summary <- arrange(cforest.summary,desc(mean_cforest))
cforest.summary <- mutate(cforest.summary,ranking=1:nrow(cforest.summary))

write.csv(cforest.summary,"../data/rankingcForest.csv")

rankingVarImp <- merge(forest.summary,cforest.summary)
write.csv(rankingVarImp,"../data/rankingVarImp.csv",row.names = FALSE)

###################################################################################
###################################################################################
## Lasso simulations

## number of simulations
nsims <- 100

## train testing sizes
sizes <- seq(0.5,0.90,0.05)

## matrix to store results of simulations
results <- matrix(0,length(sizes)*nsims,4)
results <- as.data.frame(results)
colnames(results) <- c("mse_lasso","r_sq_lasso",
                       "size","time_lasso")


## iteration of lasso, ridge, and random forest regression simulations
for ( i in 1:length(sizes) ) {
        ## train test definition
        intrain <- createDataPartition(output,p=sizes[i],list = FALSE)
        
        for ( j in 1:nsims ) {
                ##########################################################
                ## lasso 
                
                ## lasso Lambda optimal calculation
                initialTime <- Sys.time()
                lasso <- cv.glmnet(x=input[intrain,],y = output[intrain])
                totalTime <- Sys.time() - initialTime 
                
                y_predicted <- predict(lasso, newx=input[-intrain,],
                                       c=lasso$lambda.min,
                                       type="response")
                
                mss_lasso <- mean((y_predicted-output[-intrain])^2)
                mse_lasso <- sqrt(mss_lasso)
                
                r_sq_lasso <- r_sq(observed=output[-intrain],
                                   predicted=y_predicted)
                
                index <- (i-1)*nsims + j
                results[index, "mse_lasso"] <- mse_lasso
                results[index,"r_sq_lasso"] <- r_sq_lasso
                results[index,"size"] <- sizes[i]
                results[index,"time_lasso"] <- totalTime
        }
}


write.csv(results,"../data/lassoTrainSizes.csv",
          row.names = FALSE)

results.m <- melt(results,id="size")


## calculating mean of mse and R^2 values
means <- with(results.m,aggregate(value,list(size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results.m, aggregate(value,list(size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")
results.summary

## select accuracy values for plot
mse.summary <- results.summary[grepl("mse",results.summary$technique),] 

color <- "steelblue"

## plot mse values
theme_set(theme_light())
g <- ggplot(mse.summary, aes(x=size,y=mean)) 
g <- g +  geom_line(colour=color)
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g <- g + geom_point(size=3.5,colour=color)
g <- g + labs(x="Size", y="MSE")
g <- g + theme(text = element_text(size=22,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + theme(legend.position = "none")


## select accuracy values for plot
rsq.summary <- results.summary[grepl("r_sq",results.summary$technique),] 

## plot rsq values
g1 <- ggplot(rsq.summary, aes(x=size,y=mean)) 
g1 <- g1 +  geom_line(colour=color)
g1 <- g1 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g1 <- g1 + geom_point(size=3.5,colour=color)
g1 <- g1 + labs(x="Size", y="R^2")
g1 <- g1 + theme(text = element_text(size=22,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g1 <- g1 + theme(legend.position = "none")


time.summary <- results.summary[grepl("time",results.summary$technique),]

## plot rsq values
g2 <- ggplot(time.summary, aes(x=size,y=mean)) 
g2 <- g2 +  geom_line(colour=color)
g2 <- g2 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g2 <- g2 + geom_point(size=3.5,colour=color)
g2 <- g2 + labs(x="Size", y="Time (secs)")
g2 <- g2 + theme(text = element_text(size=22,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g2 <- g2 + theme(legend.position = "none")


grid.arrange(g,g1,g2,nrow=1)

##########################################################################
## performing simulations

n_sims <- 100
sizes <- seq(0.5,0.9,0.05)

## matrix to store results of simulations
results <- matrix(0,length(sizes)*n_sims,4)
colnames(results) <- c("mse","r_sq","time","size")
results <- as.data.frame(results)


## iteration of ridge regression simulations
for ( i in 1:length(sizes) ) {
        for ( j in 1:n_sims ) {
                
                train <- createDataPartition(output,p=sizes[i],list=FALSE)
                
                ## Ridge regression
                initialTime <- Sys.time()
                ridge <- cv.glmnet(x=input[train,],y=output[train],
                                   alpha=0)
                totalTime <- Sys.time() - initialTime
                
                ## performing predictions
                y_predicted_r <- predict(ridge, newx=input[-train,],
                                         c=lasso$lambda.min,
                                         type="response")
                
                ## ridge regression mse 
                mss_ridge <- mean((y_predicted_r-output[-train])^2)
                mse_ridge <- sqrt(mss_ridge)
                
                ## Rigde regression r^2 calculation
                r_sq_ridge <- r_sq(observed=output[-train],
                                   predicted=y_predicted_r) 
                
                index <- (i-1)*n_sims + j
                
                results$mse[index] <- mse_ridge
                results$r_sq[index] <- r_sq_ridge
                results$size[index] <- sizes[i]
                results$time[index] <- totalTime
                
        }
}

write.csv(results,"../data/ridgeRegTrainSizes.csv",
          row.names = FALSE)

results.m <- melt(results,id="size")


## calculating mean of mse and R^2 values
means <- with(results.m,aggregate(value,list(size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results.m, aggregate(value,list(size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")
results.summary

## select accuracy values for plot
mse.summary <- results.summary[grepl("mse",results.summary$technique),] 

color <- "steelblue"

## plot mse values
theme_set(theme_light())
g <- ggplot(mse.summary, aes(x=size,y=mean)) 
g <- g +  geom_line(colour=color)
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g <- g + geom_point(size=3.5,colour=color)
g <- g + labs(x="Size", y="MSE")
g <- g + theme(text = element_text(size=22,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + theme(legend.position = "none")


## select accuracy values for plot
rsq.summary <- results.summary[grepl("r_sq",results.summary$technique),] 

## plot rsq values
g1 <- ggplot(rsq.summary, aes(x=size,y=mean)) 
g1 <- g1 +  geom_line(colour=color)
g1 <- g1 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g1 <- g1 + geom_point(size=3.5,colour=color)
g1 <- g1 + labs(x="Size", y="R^2")
g1 <- g1 + theme(text = element_text(size=22,face="bold"),
                 axis.line = element_line(colour = 'black', size = 0.7))
g1 <- g1 + theme(legend.position = "none")


time.summary <- results.summary[grepl("time",results.summary$technique),]

## plot rsq values
g2 <- ggplot(time.summary, aes(x=size,y=mean)) 
g2 <- g2 +  geom_line(colour=color)
g2 <- g2 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01,colour=color) 
g2 <- g2 + geom_point(size=3.5,colour=color)
g2 <- g2 + labs(x="Size", y="Time (secs)")
g2 <- g2 + theme(text = element_text(size=22,face="bold"),
                 axis.line = element_line(colour = 'black', size = 0.7))
g2 <- g2 + theme(legend.position = "none")


grid.arrange(g,g1,g2,nrow=1)

#########################################################################################

nsims <- 300

results <- matrix(0,nsims,3)
results <- as.data.frame(results)
colnames(results) <- c("rforest","lasso","ridge")

for (i in 1:nsims) {
        
        ################################################################################
        ## random forest
        
        ## getting train subset
        intrain <- createDataPartition(output,p = 0.9,list = FALSE)
        
        ## perform randomForest
        rforest <- randomForest(output[intrain]~.,
                                data=input[intrain,],
                                mtry=7,
                                ntree=150)
        ## perform predictions
        rforest.pred <- predict(rforest,
                                newdata = input[-intrain,])
        ## calcualte mse values
        mssForest <- mean((rforest.pred-output[-intrain])^2)
        mseForest <- sqrt(mssForest)
        
        ###############################################################################
        ## lasso
        
        intrain <- createDataPartition(output,p = 0.85,list = FALSE)
        
        ## perform lasso
        lasso <- cv.glmnet(x=input[intrain,],y = output[intrain])
        
        ## predictions
        y_predicted <- predict(lasso, newx=input[-intrain,],
                               c=lasso$lambda.min,
                               type="response")
        
        ## gettind accuracy mse values
        mssLasso <- mean((y_predicted-output[-intrain])^2)
        mseLasso <- sqrt(mssLasso)
        
        ###############################################################################
        ## Ridge regression
        
        train <- createDataPartition(output,p=0.9,list=FALSE)
        
        ## Ridge regression
        ridge <- cv.glmnet(x=input[train,],y=output[train],
                           alpha=0)
        
        ## performing predictions
        y_predicted_r <- predict(ridge, newx=input[-train,],
                                 c=lasso$lambda.min,
                                 type="response")
        
        ## ridge regression mse 
        mssRidge <- mean((y_predicted_r-output[-train])^2)
        mseRidge <- sqrt(mssRidge)
        
        
        ## storing results
        results[i,"rforest"] <- mseForest
        results[i,"lasso"] <- mseLasso
        results[i,"ridge"] <- mseRidge
        
}

write.csv(results,"../data/benchmark.csv")

results.m <- melt(results)

## plot mse values
theme_set(theme_light())
g <- ggplot(results.m, aes(x=variable,y=log(value))) + geom_point() 
g <- g + geom_boxplot(aes(group=variable,fill=variable))
g <- g + theme( # legend.position = "none",
        text = element_text(size = 26,face = "bold")) 
g <- g + labs(x="",y="log(MSE)")
g <- g + theme(title = element_text(face="bold"),
               axis.title.x = element_text(face="bold"),
               axis.title.y = element_text(face="bold"))
g

